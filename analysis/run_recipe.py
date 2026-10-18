"""Design the minimum-toxicity glycerol + sorbitol recipe.

Inputs: data/design_inputs.csv. If any input is still PLACEHOLDER, every output
is stamped as a method demonstration, not a result.

    python analysis/run_recipe.py -> predictions/recipe_v1.csv, figures/recipe_map.png
"""
import csv, sys, time
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from thermo import GroupContributionModel, IdealModel, residual_liquid
from thermo.salt import AdditiveSaltModel
from thermo.optimise import design_recipe

inp = {r["key"]: r for r in csv.DictReader(open(ROOT / "data" / "design_inputs.csv"))}
val = lambda k: float(inp[k]["value"])
placeholder = [k for k, r in inp.items() if r["status"] == "PLACEHOLDER"]
TAG = "METHOD DEMO - placeholder inputs: " + ", ".join(placeholder) if placeholder else "Inputs sourced"
a_oh = float(open(ROOT / "predictions" / "predictions_v1.csv").readline().split("a_OH = ")[1].split()[0])

T = val("T_design_C") + 273.15
base = {"NaCl": val("base_NaCl")}
caps = {"glycerol": val("cap_glycerol"), "sorbitol": val("cap_sorbitol")}
out = []
for label, poly in [("per-OH model", GroupContributionModel(a_oh)), ("ideal polyols", IdealModel())]:
    t0 = time.time()
    model = AdditiveSaltModel(poly)
    r = design_recipe("glycerol", "sorbitol", T, model, val("salt_target"), caps, base=base,
                      w_max=0.40, n_grid=41, include_dcp=True, min_a=val("min_glycerol"))
    out.append([label, f"{100*r['w_a']:.1f}", f"{100*r['w_b']:.1f}", f"{100*r['salt_residual']:.2f}", f"{r['cost']:.3f}"])
    print(f"{label:14} glycerol {100*r['w_a']:.1f} wt%, sorbitol {100*r['w_b']:.1f} wt%, "
          f"residual NaCl {100*r['salt_residual']:.2f} wt%, cost {r['cost']:.3f}  ({time.time()-t0:.0f}s)")
    if label == "per-OH model":
        best = r
none = residual_liquid(base, T, AdditiveSaltModel(IdealModel()), include_dcp=True)
print(f"no cryoprotectant: residual NaCl {100*none['composition']['NaCl']:.1f} wt% at {val('T_design_C')} C")
print(TAG)

with open(ROOT / "predictions" / "recipe_v1.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow([f"# {TAG}"])
    w.writerow(["model", "glycerol_wt", "sorbitol_wt", "residual_NaCl_wt", "cost"]); w.writerows(out)

fig, ax = plt.subplots(figsize=(4.8, 3.8))
S = 100 * best["salt_grid"].T
cs = ax.contour(100 * best["wa"], 100 * best["wb"], S, levels=[2, 5, 8, 12], colors="black", linewidths=0.8)
ax.clabel(cs, fmt="%g wt%% NaCl", fontsize=7)
ax.contour(100 * best["wa"], 100 * best["wb"], best["cost_grid"].T, levels=8, colors="0.7", linewidths=0.6, linestyles="--")
ax.plot(100 * best["w_a"], 100 * best["w_b"], "k*", ms=10)
ax.set(xlabel="Glycerol added / wt %", ylabel="Sorbitol added / wt %")
ax.set_title(TAG, fontsize=6)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout(); fig.savefig(ROOT / "figures" / "recipe_map.png", dpi=170)
