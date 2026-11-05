"""Extension E1 - blind predictions (21 Oct) for CRC alcohols and sugars, written BEFORE their CRC
rows were viewed. One-parameter per-OH model with the committed a_OH (glycerol only), ideal solution,
UNIFAC and Dortmund (alcohols only: sugars' ring/ether/carbonyl groups are outside our simple group sets).
    python analysis/predict_e1.py -> predictions/e1_alcohols_sugars_v1.csv
"""
import csv, math, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from unifac_compare import T_freeze, ln_aw_ice, R
from scipy.optimize import brentq
ROOT = Path(__file__).resolve().parents[1]
A_OH = float(open(ROOT / "predictions" / "predictions_v1.csv").readline().split("a_OH = ")[1].split()[0])
SOL = {"methanol": (32.042, 1), "ethanol": (46.069, 1), "1-propanol": (60.096, 1), "2-propanol": (60.096, 1),
       "sucrose": (342.297, 8), "glucose": (180.156, 5), "fructose": (180.156, 5)}
ALCOHOLS = {"methanol", "ethanol", "1-propanol", "2-propanol"}

def T_model(sol, w, a):
    M, n = SOL[sol]
    ns, nw = w / M, (100 - w) / 18.015
    x = ns / (ns + nw); xw = 1 - x
    f = lambda T: math.log(xw) + x * x * n * a / (R * T) - ln_aw_ice(T)
    return brentq(f, 150.0, 273.15 - 1e-9) - 273.15

grid = [0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 24, 28, 32, 36, 40]
rows = []
for sol in SOL:
    for w in grid:
        try:
            r = [T_model(sol, w, A_OH), T_model(sol, w, 0.0)]
        except ValueError:
            continue
        if sol in ALCOHOLS:
            try:
                r += [T_freeze({sol: w}, "UNIFAC"), T_freeze({sol: w}, "UNIFAC-Dortmund")]
            except ValueError:
                r += [float("nan"), float("nan")]
        else:
            r += [float("nan"), float("nan")]
        rows.append([sol, w] + [f"{v:.3f}" for v in r])
with open(ROOT / "predictions" / "e1_alcohols_sugars_v1.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow([f"# Extension E1 blind predictions; a_OH = {A_OH} (glycerol only); written before CRC alcohol/sugar rows were viewed"])
    w.writerow(["solute", "wt_percent", "T_ours_C", "T_ideal_C", "T_UNIFAC_C", "T_Dortmund_C"])
    w.writerows(rows)
print(len(rows), "predictions written")
