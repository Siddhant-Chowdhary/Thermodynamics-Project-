"""Sensitivity of the glycerol-sparing recipe to the salt target (+-25%) and design temperature, three polyol models.
Uses the committed recipe code (analysis/glycerol_sparing.py) unchanged; only the salt target is varied.
Run from the repo root:  python analysis/recipe_sensitivity.py -> predictions/recipe_sensitivity.csv
"""
import csv, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import glycerol_sparing as gs

ROOT = Path(__file__).resolve().parents[1]
base = gs.SALT_TARGET
rows = []
print(f"base salt target {100*base:.2f} wt%")
print("T (C) | target (wt%) | model | min glycerol, no sorbitol (wt%) | glycerol saved per g sorbitol (0->2 wt%)")
for T in (-10.0, -15.0, -20.0):
    for f in (0.75, 1.0, 1.25):
        gs.SALT_TARGET = base * f
        for mdl in ("per-OH", "ideal", "UNIFAC"):
            g0 = gs.min_glycerol(0.0, T + 273.15, mdl); g2 = gs.min_glycerol(0.02, T + 273.15, mdl)
            ratio = (g0 - g2) / 0.02
            rows.append([T, round(100 * base * f, 2), mdl, round(100 * g0, 2), round(ratio, 2)])
            print(f"{T:5.1f} | {100*base*f:5.2f} | {mdl:6s} | {100*g0:5.2f} | {ratio:.2f}")
gs.SALT_TARGET = base
with open(ROOT / "predictions" / "recipe_sensitivity.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["# salt target varied +-25% around the sourced 4.41 wt% (Lovelock 1953); 0.9 wt% NaCl base; models as in recipe_v2"])
    w.writerow(["T_design_C", "salt_target_wt", "polyol_model", "min_glycerol_no_sorbitol_wt", "glycerol_saved_per_g_sorbitol"])
    w.writerows(rows)
