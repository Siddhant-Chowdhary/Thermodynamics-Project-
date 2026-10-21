"""Blind predictions for ethylene glycol (2 OH, no CH3) and D-mannitol (6 OH, sorbitol isomer),
written 21 Oct BEFORE looking at their CRC rows (only the CRC solute list had been seen).
Same committed a_OH; four models; CRC concentration grid.
    python analysis/predict_crc_extra.py -> predictions/crc_extra_v1.csv
"""
import csv, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from unifac_compare import T_freeze, ln_aw_ice, mole_fracs, R
from scipy.optimize import brentq
import math
ROOT = Path(__file__).resolve().parents[1]
A_OH = float(open(ROOT / "predictions" / "predictions_v1.csv").readline().split("a_OH = ")[1].split()[0])
NOH = {"ethylene_glycol": 2, "mannitol": 6}

def T_model(wt, a):
    names, xs = mole_fracs(wt)
    f = lambda T: math.log(xs[0]) + (1 - xs[0]) * sum(NOH[k] * a / (R * T) * x for k, x in zip(names[1:], xs[1:])) - ln_aw_ice(T)
    return brentq(f, 180.0, 273.15 - 1e-9) - 273.15

grid = [0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 24, 28, 32, 36, 40]
rows = []
for sol, wmax in [("ethylene_glycol", 40), ("mannitol", 18)]:
    for w in [g for g in grid if g <= wmax]:
        wt = {sol: w}
        r = [T_model(wt, A_OH), T_model(wt, 0.0), T_freeze(wt, "UNIFAC"), T_freeze(wt, "UNIFAC-Dortmund")]
        rows.append([sol, w] + [f"{v:.3f}" for v in r])
with open(ROOT / "predictions" / "crc_extra_v1.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow([f"# blind; a_OH = {A_OH} J/mol per OH (from predictions_v1); written before CRC EG/mannitol rows were viewed"])
    w.writerow(["solute", "wt_percent", "T_ours_C", "T_ideal_C", "T_UNIFAC_C", "T_Dortmund_C"])
    w.writerows(rows)
for r in rows:
    if r[1] in (5, 10, 18, 20, 40): print(r)
