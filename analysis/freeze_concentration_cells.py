"""Freeze concentration outside the cell, what it does inside the cell, and what that means for recipe design.
(Application example; equilibrium thermodynamics only; uses the committed model, nothing refitted.)

OUTSIDE: ice is pure water, so solutes stay in the shrinking liquid until its water activity equals that of ice at T
(thermo.residual_liquid; per-OH polyols + fitted NaCl phi, additive in ln a_w).

INSIDE: water crosses the membrane until its chemical potential is equal on both sides, so intracellular a_w also equals
a_w(ice, T). Assumptions (stated in the report):
  * the cell's own impermeant solutes are represented as NaCl at the isotonic level (0.9 wt% saline), so the cell is
    isotonic before freezing and no further cell parameter is needed (we report cell WATER, not total volume);
  * glycerol is fully penetrating and equilibrated: intracellular molality = extracellular molality at T;
  * sorbitol and NaCl do not cross the membrane;
  * slow cooling: the cell stays at osmotic equilibrium (no intracellular ice).

DESIGN: minimum glycerol vs sorbitol at -15 C under three literature injury criteria:
  (1) NaCl in the unfrozen liquid <= 4.41 wt%   (Lovelock 1953, via Shepard 1976)     -> recipe_v2 criterion
  (2) the same limit for the salt INSIDE the cell (equilibrium cell model above)
  (3) unfrozen fraction of the medium's water >= 25 %  (Mazur, Rall & Rigopoulos 1981: most red cells survived when
      less than 75 % of the water froze, even with > 2 molal salt)

    python analysis/freeze_concentration_cells.py -> predictions/freeze_concentration_cells.csv,
        predictions/recipe_three_criteria.csv, figures/freeze_concentration_cells.png
"""
import csv, sys, warnings
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys as _sys
_sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from analysis import style as S

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from thermo import GroupContributionModel, residual_liquid, ln_aw_ice
from thermo.salt import AdditiveSaltModel
from thermo.constants import solute

warnings.filterwarnings("ignore", message="NaCl molality")
a_oh = float(open(ROOT / "predictions" / "predictions_v1.csv").readline().split("a_OH = ")[1].split()[0])
MODEL = AdditiveSaltModel(GroupContributionModel(a_oh))
M = {k: solute(k)["molar_mass"] / 1000 for k in ("glycerol", "sorbitol", "NaCl")}     # kg/mol
BASE = 0.009
M_S0 = BASE / M["NaCl"] / (1 - BASE)            # isotonic NaCl molality (mol/kg water)
SALT_LIMIT, UNFROZEN_MIN, T_DESIGN = 4.41, 0.25, -15.0

CASES = {
    "saline only": {},
    "sorbitol 9 wt%": {"sorbitol": 0.09},
    "glycerol 6.5 wt%": {"glycerol": 0.065},
    "glycerol 2.75 + sorbitol 5 wt%": {"glycerol": 0.0275, "sorbitol": 0.05},
}


def state(comp, TC):
    """Return unfrozen water fraction, NaCl wt% outside, NaCl-equiv wt% inside, cell water W/W0 at TC (deg C)."""
    T = TC + 273.15
    mf = dict({k: v for k, v in comp.items() if v > 0}, NaCl=BASE)
    r = residual_liquid(mf, T, MODEL, include_dcp=True)
    c = r["composition"]
    w_liq_water = r["liquid_fraction"] * (1 - sum(c.values()))
    unfrozen_water = w_liq_water / (1 - sum(mf.values()))
    m_g = c.get("glycerol", 0.0) / M["glycerol"] / (1 - sum(c.values()))
    target = MODEL.ln_aw(c, T)          # = ln a_w(ice, T) once ice is present; = medium's a_w if no ice yet

    def f(W):
        ms = M_S0 / W
        tot = 1 + ms * M["NaCl"] + m_g * M["glycerol"]
        x = {"NaCl": ms * M["NaCl"] / tot}
        if m_g > 0:
            x["glycerol"] = m_g * M["glycerol"] / tot
        return MODEL.ln_aw(x, T) - target
    W = brentq(f, 1e-3, 1e4, xtol=1e-10)
    ms = M_S0 / W
    salt_in = 100 * ms * M["NaCl"] / (1 + ms * M["NaCl"] + m_g * M["glycerol"])
    return unfrozen_water, 100 * c["NaCl"], salt_in, W


if __name__ == "__main__":
    Ts = np.arange(-1.0, -25.01, -0.5)
    rows, res = [], {k: [] for k in CASES}
    for name, comp in CASES.items():
        for TC in Ts:
            try:
                u, so, si, W = state(comp, TC)
            except ValueError:
                continue
            res[name].append((TC, u, so, si, W))
            rows.append([name, f"{TC:.1f}", f"{u:.4f}", f"{so:.3f}", f"{si:.3f}", f"{W:.4f}"])
    with open(ROOT / "predictions" / "freeze_concentration_cells.csv", "w", newline="") as fh:
        fh.write("# per-OH (committed a_OH) + fitted NaCl, additive ln a_w; equilibrium cell model (assumptions in script)\n")
        w = csv.writer(fh)
        w.writerow(["case", "T_C", "unfrozen_water_fraction", "NaCl_outside_wt", "salt_inside_wt_NaCl_equiv", "cell_water_W_over_W0"])
        w.writerows(rows)
    print(f"At {T_DESIGN} C: unfrozen water | NaCl outside | salt inside | cell water")
    for name in CASES:
        r = [x for x in res[name] if abs(x[0] - T_DESIGN) < 1e-6][0]
        print(f"  {name:32s} {r[1]:6.3f} | {r[2]:6.2f} | {r[3]:6.2f} | {r[4]:6.3f}")

    # ---- design under three criteria
    def need(ws, crit):
        def g(wg):
            u, so, si, _ = state({"glycerol": wg, "sorbitol": ws}, T_DESIGN)
            return {"outside": so - SALT_LIMIT, "inside": si - SALT_LIMIT, "unfrozen": UNFROZEN_MIN - u}[crit]
        if g(0.0) <= 0:
            return 0.0
        if g(0.35) > 0:
            return float("nan")           # needs more than 35 wt% glycerol
        return brentq(g, 0.0, 0.35, xtol=1e-6)

    WS = [i / 100 for i in range(0, 13)]
    design = {c: [need(ws, c) for ws in WS] for c in ("outside", "inside", "unfrozen")}
    with open(ROOT / "predictions" / "recipe_three_criteria.csv", "w", newline="") as fh:
        fh.write(f"# min glycerol (wt%) at {T_DESIGN} C in 0.9 wt% saline, per-OH model; criteria: NaCl outside <= {SALT_LIMIT} wt%, "
                 f"salt inside <= {SALT_LIMIT} wt% (equilibrium cell), unfrozen water >= {UNFROZEN_MIN}\n")
        w = csv.writer(fh)
        w.writerow(["sorbitol_wt", "glycerol_outside_criterion", "glycerol_inside_criterion", "glycerol_unfrozen_criterion"])
        for i, ws in enumerate(WS):
            w.writerow([f"{100*ws:.0f}"] + [f"{100*design[c][i]:.2f}" for c in ("outside", "inside", "unfrozen")])
    print("\nMin glycerol (wt%) at -15 C   sorbitol: " + " ".join(f"{100*x:5.0f}" for x in WS))
    for c in design:
        print(f"  {c:9s}                     " + " ".join(f"{100*x:5.2f}" for x in design[c]))
    for c in design:
        d = design[c]
        print(f"  slope 0->2 wt% sorbitol, {c}: {(d[0]-d[2])/0.02:+.2f} g glycerol saved per g sorbitol (negative = extra glycerol needed)")

    # all criteria at once: inside salt >= outside salt whenever sorbitol is present, so the binding pair is
    # (inside, unfrozen); the optimum sits where the two curves cross
    ws_opt = brentq(lambda ws: need(ws, "inside") - need(ws, "unfrozen"), 0.0, 0.06, xtol=1e-6)
    wg_opt = need(ws_opt, "inside")
    print(f"\nAll three criteria: glycerol alone needs {100*design['unfrozen'][0]:.2f} wt%; optimum sorbitol {100*ws_opt:.2f} wt% "
          f"+ glycerol {100*wg_opt:.2f} wt% ({100*(1-wg_opt/design['unfrozen'][0]):.0f}% less glycerol)")
    with open(ROOT / "predictions" / "recipe_three_criteria.csv", "a", newline="") as fh:
        fh.write(f"# all criteria: optimum sorbitol {100*ws_opt:.2f} wt% + glycerol {100*wg_opt:.2f} wt% "
                 f"(glycerol alone {100*design['unfrozen'][0]:.2f} wt%)\n")

    # ---- figure (2 x 2)
    plt.rcParams.update({"font.size": 6.5, "axes.linewidth": 0.6})
    fig, ax = plt.subplots(2, 2, figsize=(5.1, 3.35))
    ax = ax.ravel()
    styles = {"saline only": (":", S.DATA), "sorbitol 9 wt%": ("-.", S.ALT), "glycerol 6.5 wt%": ("-", S.OURS),
              "glycerol 2.75 + sorbitol 5 wt%": ("--", S.OV)}
    for name, (ls, cc) in styles.items():
        d = np.array(res[name])
        ax[0].plot(d[:, 0], d[:, 1], ls, color=cc, lw=1.4, label=name)
        ax[1].plot(d[:, 0], d[:, 2], ls, color=cc, lw=1.4)
        ax[2].plot(d[:, 0], d[:, 3], ls, color=cc, lw=1.4)
    ax[0].axhspan(0, UNFROZEN_MIN, color=S.ALT, alpha=0.07, lw=0); ax[0].axhline(UNFROZEN_MIN, color=S.ALT, lw=0.7)
    ax[0].text(-16, UNFROZEN_MIN + 0.02, "25 % unfrozen (Mazur 1981)", fontsize=5.0, color="0.35")
    for a in ax[1:3]:
        a.axhspan(SALT_LIMIT, 24, color=S.ALT, alpha=0.07, lw=0); a.axhline(SALT_LIMIT, color=S.ALT, lw=0.7)
        a.text(-16, SALT_LIMIT - 1.8, "haemolysis onset 4.41 wt%", fontsize=5.0, color="0.35")
        a.set_ylim(0, 24)
    for a in ax[:3]:
        a.axvline(T_DESIGN, color="0.75", lw=0.6, ls=":")
        a.invert_xaxis(); a.set_xlabel("Temperature (°C)")
    ax[0].set(ylabel="Unfrozen fraction of water", title="(a) Liquid water left", ylim=(0, 1.02))
    ax[1].set(ylabel="NaCl outside (wt%)", title="(b) Salt outside the cell")
    ax[2].set(ylabel="Salt inside (wt%, NaCl-equiv.)", title="(c) Salt inside the cell (equilibrium)")
    ax[2].text(-0.5, 22, "saline only = sorbitol only", fontsize=5.0, color="0.25")
    ax[0].legend(fontsize=5.0, frameon=False, loc="upper right")
    xs = [100 * x for x in WS]
    for c, ls, lab, cc in (("outside", "-", "salt outside ≤ 4.41 wt%", S.UNIFAC), ("inside", "--", "salt inside ≤ 4.41 wt%", S.ALT),
                           ("unfrozen", ":", "unfrozen water ≥ 25 %", S.OURS)):
        ax[3].plot(xs, [100 * v for v in design[c]], ls, color=cc, lw=1.5, label=lab)
    comb = [100 * max(a, b) for a, b in zip(design["inside"], design["unfrozen"])]
    ax[3].plot(xs, comb, color=S.KM, alpha=0.35, lw=4.0, zorder=0, label="all three (binding)")
    ax[3].plot([100 * ws_opt], [100 * wg_opt], "o", color=S.KM, mec="k", mew=0.6, ms=5.5)
    ax[3].annotate(f"optimum: {100*ws_opt:.1f} + {100*wg_opt:.1f}", (100 * ws_opt, 100 * wg_opt), xytext=(4.5, 3.0),
                   fontsize=6, arrowprops=dict(arrowstyle="-", lw=0.5))
    ax[3].set_ylim(0, 30)
    ax[3].set(xlabel="Sorbitol added (wt%)", ylabel="Minimum glycerol (wt%)", title="(d) Recipe at −15 °C: depends on the criterion")
    ax[3].legend(fontsize=5.2, frameon=False)
    for a in ax:
        a.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(pad=0.3, h_pad=0.8, w_pad=1.0)
    fig.savefig(ROOT / "figures" / "freeze_concentration_cells.png", dpi=300)
    fig.savefig(ROOT / "figures" / "freeze_concentration_cells.pdf")
