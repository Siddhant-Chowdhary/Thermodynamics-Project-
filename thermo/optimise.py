"""Choosing the mixing ratio of two cryoprotectants.

Key thermodynamic point: once ice is present at temperature T, the water
activity (and osmotic pressure) of the remaining liquid is fixed by T alone.
A cryoprotectant cannot change that. What it changes is how MUCH liquid is left
and WHAT that liquid contains. So the objectives here are about composition:

  "salt"      minimise NaCl wt% in the residual liquid (the classic
              salt-dilution mechanism of cryoprotection). Requires a saline base.
  "toxicity"  minimise the worst normalised cryoprotectant exposure,
              max_i (c_i / cap_i), with a user-supplied concentration cap per
              cryoprotectant. Caps must come from the literature and be cited.
  "liquid"    maximise the fraction of the sample still liquid.

The total cryoprotectant dose is fixed (the toxicity budget); only the ratio
between the two agents varies. A grid search is used so the whole trade-off
curve is available for plotting, not just the optimum.
"""

import numpy as np

from .sle import residual_liquid


def design_recipe(solute_a, solute_b, T_store, model, salt_target, caps,
                  base=None, w_max=0.40, n_grid=61, include_dcp=False, min_a=0.0):
    """Minimum-toxicity recipe that meets a thermodynamic protection target.

    This is the well-posed design problem. (The fixed-total-dose version in
    optimise_ratio with objective='toxicity' is NOT: its optimum ratio equals
    caps_a / (caps_a + caps_b) whatever the thermodynamic model, because both
    solutes are conserved in the residual liquid.)

    Here both doses are free. Protection target: the NaCl mass fraction in the
    residual liquid at T_store must not exceed `salt_target` (salt-dilution
    mechanism). Toxicity cost of the added agents:
        cost = (w_a / cap_a)^2 + (w_b / cap_b)^2
    How much NaCl each gram of cryoprotectant removes depends on molar mass and
    on the activity model, so the optimum depends on the thermodynamics.

    salt_target, caps : must be taken from the literature and cited.
    Returns the best (w_a, w_b), its cost, and the full grid for plotting.
    """
    base = dict(base or {})
    if "NaCl" not in base:
        raise ValueError("design_recipe needs NaCl in `base` (e.g. {'NaCl': 0.009})")
    wa = np.linspace(0, w_max, n_grid)
    wb = np.linspace(0, w_max, n_grid)
    salt = np.full((n_grid, n_grid), np.nan)
    cost = (wa[:, None] / caps[solute_a]) ** 2 + (wb[None, :] / caps[solute_b]) ** 2
    for i, a in enumerate(wa):
        for j, b in enumerate(wb):
            comp = dict(base)
            if a > 0:
                comp[solute_a] = a
            if b > 0:
                comp[solute_b] = b
            if sum(comp.values()) >= 0.95:
                continue
            try:
                salt[i, j] = residual_liquid(comp, T_store, model, include_dcp)["composition"]["NaCl"]
            except ValueError:
                continue
    feasible = np.where((salt <= salt_target) & (wa[:, None] >= min_a), cost, np.inf)
    if not np.isfinite(feasible).any():
        raise ValueError("No recipe on the grid meets the salt target; raise w_max")
    i, j = np.unravel_index(np.argmin(feasible), feasible.shape)
    return {"w_a": float(wa[i]), "w_b": float(wb[j]), "cost": float(cost[i, j]),
            "frac_a": float(wa[i] / (wa[i] + wb[j])) if wa[i] + wb[j] > 0 else float("nan"),
            "salt_residual": float(salt[i, j]), "wa": wa, "wb": wb,
            "salt_grid": salt, "cost_grid": cost}


def evaluate_ratio(frac_a, solute_a, solute_b, total_cpa, T_store, model,
                   base=None, include_dcp=False):
    base = dict(base or {})
    comp = dict(base)
    comp[solute_a] = comp.get(solute_a, 0) + total_cpa * frac_a
    comp[solute_b] = comp.get(solute_b, 0) + total_cpa * (1 - frac_a)
    comp = {k: v for k, v in comp.items() if v > 0}
    return residual_liquid(comp, T_store, model, include_dcp)


def optimise_ratio(solute_a, solute_b, total_cpa, T_store, model, objective="salt",
                   base=None, caps=None, n_grid=101, include_dcp=False):
    """Scan frac_a from 0 to 1. Returns best ratio and the full curve.

    total_cpa : total cryoprotectant mass fraction (e.g. 0.20 for 20 wt%)
    base      : other solutes in the medium, e.g. {"NaCl": 0.009} for saline
    caps      : {"glycerol": 0.40, ...} residual mass-fraction caps for "toxicity"
    """
    fracs = np.linspace(0, 1, n_grid)
    scores, results = [], []
    for f in fracs:
        r = evaluate_ratio(f, solute_a, solute_b, total_cpa, T_store, model, base, include_dcp)
        comp = r["composition"]
        if objective == "salt":
            if "NaCl" not in comp:
                raise ValueError("objective='salt' needs NaCl in `base`")
            s = comp["NaCl"]
        elif objective == "toxicity":
            if not caps:
                raise ValueError("objective='toxicity' needs `caps`")
            s = max(comp.get(k, 0) / caps[k] for k in (solute_a, solute_b))
        elif objective == "liquid":
            s = -r["liquid_fraction"]
        elif callable(objective):
            s = objective(r)
        else:
            raise ValueError(f"Unknown objective {objective!r}")
        scores.append(s)
        results.append(r)
    scores = np.array(scores)
    i = int(np.argmin(scores))
    return {"best_frac_a": float(fracs[i]), "best_score": float(scores[i]),
            "fracs": fracs, "scores": scores, "results": results}
