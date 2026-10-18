"""Run with:  python -m pytest -q   (from the repo root)"""

import math
import numpy as np

from thermo import (IdealModel, MargulesModel, GroupContributionModel, freezing_point,
                    residual_liquid, osmotic_pressure_at_ice, from_wt_percent, mole_fractions)
from thermo.fit import fit_group_contribution, bootstrap, gibbs_duhem_residual
from thermo.optimise import optimise_ratio, design_recipe
from analysis import simulate, find_freezing_point


def test_mole_fraction_10wt_glycerol():
    xw, xs = mole_fractions(from_wt_percent(glycerol=10))
    assert abs(xw - 0.97873) < 1e-4


def test_ideal_freezing_point_matches_hand_calculation():
    Tf = freezing_point(from_wt_percent(glycerol=10), IdealModel())
    assert abs((Tf - 273.15) - (-2.20)) < 0.01


def test_pure_water_freezes_at_Tm():
    assert freezing_point({}, IdealModel()) == 273.15


def test_group_contribution_reduces_to_margules_for_one_solute():
    a = -1500.0
    gc = GroupContributionModel(a)
    mg = MargulesModel({"glycerol": 3 * a})
    c = from_wt_percent(glycerol=20)
    assert abs(freezing_point(c, gc) - freezing_point(c, mg)) < 1e-9


def test_gibbs_duhem_holds():
    assert gibbs_duhem_residual(GroupContributionModel(-2000), 265.0, "sorbitol") < 1e-8


def test_freeze_concentration_mass_balance():
    r = residual_liquid(from_wt_percent(glycerol=10), 263.15, IdealModel())
    # solute conserved: liquid_fraction * residual wt% = starting wt%
    assert abs(r["liquid_fraction"] * r["composition"]["glycerol"] - 0.10) < 1e-9
    assert abs(r["liquid_fraction"] - 0.285) < 0.002


def test_osmotic_pressure_at_minus10():
    assert abs(osmotic_pressure_at_ice(263.15) / 1e6 - 12.2) < 0.1


def test_fit_recovers_known_parameter():
    true = GroupContributionModel(-1500)
    data = [(from_wt_percent(glycerol=w), freezing_point(from_wt_percent(glycerol=w), true))
            for w in (5, 10, 15, 20, 25, 30)]
    assert abs(fit_group_contribution(data)["param"] + 1500) < 1.0


def test_bootstrap_gives_finite_interval():
    rng = np.random.default_rng(1)
    true = GroupContributionModel(-1500)
    data = [(from_wt_percent(glycerol=w),
             freezing_point(from_wt_percent(glycerol=w), true) + rng.normal(0, 0.05))
            for w in (5, 10, 15, 20, 25, 30)]
    b = bootstrap(data, fit_group_contribution, n_boot=100)
    assert b["ci95"][0] < b["mean"] < b["ci95"][1]


def test_analyser_recovers_simulated_freezing_point():
    for Tf in (-1.0, -3.0, -6.0):
        t, T = simulate(Tf, seed=3)
        r = find_freezing_point(t, T, calib_unc_C=0.0)
        assert abs(r["T_freeze_C"] - Tf) < 0.1


def test_fixed_dose_toxicity_optimum_is_set_by_caps_alone():
    # documents WHY design_recipe exists: this formulation ignores thermodynamics
    for model in (IdealModel(), GroupContributionModel(-2500)):
        o = optimise_ratio("glycerol", "sorbitol", 0.20, 263.15, model, objective="toxicity",
                           base={"NaCl": 0.009}, caps={"glycerol": 0.3, "sorbitol": 0.5},
                           n_grid=201)
        assert abs(o["best_frac_a"] - 0.375) < 0.006


def test_design_recipe_depends_on_thermodynamics():
    caps = {"glycerol": 0.30, "sorbitol": 0.40}
    kw = dict(salt_target=0.03, caps=caps, base={"NaCl": 0.009}, w_max=0.35, n_grid=71)
    r_id = design_recipe("glycerol", "sorbitol", 263.15, IdealModel(), **kw)
    r_gc = design_recipe("glycerol", "sorbitol", 263.15, GroupContributionModel(-3000), **kw)
    assert r_id["salt_residual"] <= 0.03 and r_gc["salt_residual"] <= 0.03
    assert abs(r_id["frac_a"] - r_gc["frac_a"]) > 0.05


def test_nacl_correction_reproduces_liquidus():
    from thermo.salt import AdditiveSaltModel
    m = AdditiveSaltModel(IdealModel())
    # literature: 10 wt% -6.56 C, 20 wt% -16.46 C
    assert abs(freezing_point({"NaCl": 0.10}, m, True) - 273.15 + 6.56) < 0.1
    assert abs(freezing_point({"NaCl": 0.20}, m, True) - 273.15 + 16.46) < 0.1


def test_additive_model_without_salt_equals_polyol_model():
    from thermo.salt import AdditiveSaltModel
    g = GroupContributionModel(-2400.0)
    a = freezing_point({"glycerol": 0.2}, g, True)
    b = freezing_point({"glycerol": 0.2}, AdditiveSaltModel(g), True)
    assert abs(a - b) < 1e-9
