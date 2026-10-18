"""Fitting model parameters to freezing point data, with uncertainty.

Data format: a list of (composition_dict, T_freeze_K) pairs, e.g.
    [({"glycerol": 0.10}, 270.6), ({"glycerol": 0.20}, 267.9), ...]
"""

import numpy as np
from scipy.optimize import least_squares

from .activity import GroupContributionModel, MargulesModel
from .sle import freezing_point


def _residuals(model_factory, param, data, include_dcp):
    model = model_factory(param)
    return np.array([freezing_point(c, model, include_dcp) - T for c, T in data])


def fit_parameter(data, model_factory, x0=-1000.0, include_dcp=False):
    """Least-squares fit of a single parameter (J/mol). Returns dict."""
    res = least_squares(lambda p: _residuals(model_factory, p[0], data, include_dcp),
                        x0=[x0], method="lm")
    param = float(res.x[0])
    r = res.fun
    return {"param": param, "residuals_K": r, "rmse_K": float(np.sqrt(np.mean(r ** 2))),
            "model": model_factory(param)}


def fit_group_contribution(data, **kw):
    """Fit a_OH (J/mol per -OH group)."""
    return fit_parameter(data, GroupContributionModel, **kw)


def fit_margules_single(data, solute_name, **kw):
    """Fit one Margules constant for a single solute."""
    return fit_parameter(data, lambda A: MargulesModel({solute_name: A}), **kw)


def bootstrap(data, fitter, n_boot=500, seed=0, **kw):
    """Resample data points with replacement and refit. Returns param distribution."""
    rng = np.random.default_rng(seed)
    n = len(data)
    params = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        sample = [data[i] for i in idx]
        try:
            params.append(fitter(sample, **kw)["param"])
        except Exception:
            continue
    params = np.array(params)
    lo, hi = np.percentile(params, [2.5, 97.5])
    return {"samples": params, "mean": float(params.mean()), "std": float(params.std(ddof=1)),
            "ci95": (float(lo), float(hi))}


def gibbs_duhem_residual(model, T, solute_name, n=200):
    """Max |x_w dln(gamma_w) + x_s dln(gamma_s)| over a binary composition grid.

    Should be ~0 (numerical noise only) for any model derived from a G_ex.
    """
    xs_grid = np.linspace(0.001, 0.5, n)
    lgw, lgs = [], []
    for xs in xs_grid:
        xw = 1 - xs
        lgw.append(model.ln_gamma_water(T, xw, {solute_name: xs}))
        lgs.append(model.ln_gamma_solutes(T, xw, {solute_name: xs})[solute_name])
    d_lgw = np.gradient(np.array(lgw), xs_grid)
    d_lgs = np.gradient(np.array(lgs), xs_grid)
    resid = (1 - xs_grid) * d_lgw + xs_grid * d_lgs
    return float(np.max(np.abs(resid[2:-2])))   # drop edge-difference artefacts
