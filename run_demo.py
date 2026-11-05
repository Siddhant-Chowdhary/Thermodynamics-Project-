"""End-to-end demo of the pipeline on SYNTHETIC data.

Every number here is generated from an assumed parameter plus noise, only to
show that the pipeline works. None of it is a result. Replace with
data/literature.csv and your measured runs before drawing any conclusion.

    python run_demo.py        -> figures/demo_*.png
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from thermo import IdealModel, GroupContributionModel, freezing_point, residual_liquid, from_wt_percent
from thermo.fit import fit_group_contribution, bootstrap, gibbs_duhem_residual
from thermo.optimise import design_recipe
from analysis import simulate, find_freezing_point

TAG = "SYNTHETIC DATA — pipeline demo, not results"
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
rng = np.random.default_rng(42)
TRUE_A_OH = -1500.0                         # assumed, for generating fake data only
truth = GroupContributionModel(TRUE_A_OH)


def fake(solute, wts, noise=0.05):
    return [(from_wt_percent(**{solute: w}),
             freezing_point(from_wt_percent(**{solute: w}), truth) + rng.normal(0, noise))
            for w in wts]


# 1. fit on "literature" glycerol ------------------------------------------------
gly = fake("glycerol", [5, 10, 15, 20, 25, 30])
fit = fit_group_contribution(gly)
boot = bootstrap(gly, fit_group_contribution, n_boot=300)
model = fit["model"]
print(f"a_OH = {fit['param']:.0f} J/mol per -OH  (95% CI {boot['ci95'][0]:.0f} to {boot['ci95'][1]:.0f}),"
      f"  RMSE {fit['rmse_K']:.3f} K,  Gibbs-Duhem residual {gibbs_duhem_residual(model, 265, 'glycerol'):.1e}")

# 2. blind prediction vs "measured" PG and sorbitol ------------------------------
fig, ax = plt.subplots(figsize=(5.6, 3.6))
w = np.linspace(0.5, 30, 60)
for sol, ls, mk in [("propylene_glycol", "-", "o"), ("glycerol", "--", "s"), ("sorbitol", ":", "^")]:
    ax.plot(w, [freezing_point(from_wt_percent(**{sol: x}), model) - 273.15 for x in w],
            ls, color="black", lw=1.2, label=f"{sol.replace('_', ' ')} (model)")
    pts = gly if sol == "glycerol" else fake(sol, [10, 20, 30])
    ax.plot([100 * list(c.values())[0] for c, _ in pts], [T - 273.15 for _, T in pts],
            mk, mfc="white", mec="black", ms=5, ls="none")
ax.set(xlabel="Cryoprotectant / wt %", ylabel="Freezing point / °C", title=TAG)
ax.title.set_fontsize(7)
ax.legend(frameon=False, fontsize=7)
fig.tight_layout(); fig.savefig("figures/demo_1_prediction.png", dpi=160)

# 3. cooling curve analysis --------------------------------------------------------
t, T = simulate(-2.3, seed=7, supercool_C=3.2)
cc = find_freezing_point(t, T)
fig, ax = plt.subplots(figsize=(5.6, 3.2))
ax.plot(t / 60, T, color="black", lw=0.8)
ax.axhline(cc["T_freeze_C"], color="0.5", ls="--", lw=0.8)
ax.plot(cc["t_freeze_s"] / 60, cc["T_freeze_C"], "o", mfc="white", mec="black")
ax.plot(cc["t_nucleation_s"] / 60, cc["T_supercool_C"], "v", mfc="white", mec="black")
ax.set(xlabel="Time / min", ylabel="Temperature / °C",
       title=f"{TAG}\nfreezing point {cc['T_freeze_C']:.2f} ± {cc['uncertainty_C']:.2f} °C "
             f"(supercooling {cc['supercooling_C']:.1f} °C)")
ax.title.set_fontsize(7)
ax.set_ylim(-8, 5)
fig.tight_layout(); fig.savefig("figures/demo_2_cooling_curve.png", dpi=160)

# 4. freeze concentration -----------------------------------------------------------
fig, ax = plt.subplots(figsize=(5.6, 3.4))
Ts = np.linspace(-1, -20, 60)
for sol, ls in [("propylene_glycol", "-"), ("glycerol", "--"), ("sorbitol", ":")]:
    c0 = from_wt_percent(**{sol: 10})
    fr = []
    for Tc in Ts:
        try:
            fr.append(100 * residual_liquid(c0, Tc + 273.15, model)["liquid_fraction"])
        except ValueError:
            fr.append(np.nan)
    ax.plot(Ts, fr, ls, color="black", lw=1.2, label=f"10 wt% {sol.replace('_', ' ')}")
ax.set(xlabel="Temperature / °C", ylabel="Sample still liquid / %", title=TAG)
ax.title.set_fontsize(7)
ax.invert_xaxis(); ax.legend(frameon=False, fontsize=7)
fig.tight_layout(); fig.savefig("figures/demo_3_freeze_concentration.png", dpi=160)

# 5. recipe design: ideal vs fitted model ------------------------------------------
caps = {"glycerol": 0.30, "sorbitol": 0.40}    # PLACEHOLDER toxicity caps — replace from literature
kw = dict(salt_target=0.03, caps=caps, base={"NaCl": 0.009}, w_max=0.30, n_grid=61)
fig, ax = plt.subplots(figsize=(4.8, 4.2))
for name, m, mk in [("ideal", IdealModel(), "s"), ("fitted", model, "o")]:
    r = design_recipe("glycerol", "sorbitol", 263.15, m, **kw)
    cs = ax.contour(100 * r["wa"], 100 * r["wb"], r["salt_grid"].T, levels=[kw["salt_target"]],
                    colors="black", linestyles="--" if name == "ideal" else "-", linewidths=1)
    ax.plot(100 * r["w_a"], 100 * r["w_b"], mk, mfc="white", mec="black", ms=7,
            label=f"{name}: {100*r['w_a']:.1f}% gly + {100*r['w_b']:.1f}% sorb")
    print(f"recipe ({name}): glycerol {100*r['w_a']:.1f} wt%, sorbitol {100*r['w_b']:.1f} wt%, cost {r['cost']:.3f}")
ax.set(xlabel="Glycerol added / wt %", ylabel="Sorbitol added / wt %",
       title=f"{TAG}\nLines: residual NaCl = 3 wt% at −10 °C; markers: minimum-toxicity recipe")
ax.title.set_fontsize(7)
ax.set_xlim(0, 15); ax.set_ylim(0, 15)
ax.legend(frameon=False, fontsize=7, loc="upper right")
fig.tight_layout(); fig.savefig("figures/demo_4_recipe.png", dpi=160)
print("figures written to figures/")
