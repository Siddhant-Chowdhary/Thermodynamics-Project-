"""Blind test #3 (30 Oct 2026): maltose and lactose (C12H22O11, 8 OH each), CRC 2005 'Concentrative Properties'.
Predictions are written and committed BEFORE the CRC rows for these solutes are entered into data/.
Committed a_OH, no refit. Only the 0.5 wt% rows (0.03 K, uninformative) had been seen when this was written.

    python analysis/predict_disaccharides.py -> predictions/disaccharides_v1.csv
"""
import csv, math, sys
from pathlib import Path
from scipy.optimize import brentq
ROOT = Path(__file__).resolve().parents[1]
R, MW, DH, DCP, TM = 8.314462618, 0.018015, 6009.5, 37.9, 273.15
A_OH = float(open(ROOT / "predictions" / "predictions_v1.csv").readline().split("a_OH = ")[1].split()[0])
MOL = {"maltose": (342.297, 8), "lactose": (342.297, 8)}
GRID = [0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 24, 28, 32, 36, 40, 44, 48]


def ice(T):
    return -(DH / R) * (1 / T - 1 / TM) - (DCP / R) * (TM / T - 1 + math.log(T / TM))


def Tf(w, M, n, a):
    m = w / (100 - w) * 1000 / M
    x = m / (1 / MW + m)
    return brentq(lambda T: math.log(1 - x) + n * a / (R * T) * x * x - ice(T), 150, TM - 1e-9) - TM


if __name__ == "__main__":
    with open(ROOT / "predictions" / "disaccharides_v1.csv", "w", newline="") as fh:
        fh.write(f"# BLIND, committed 30 Oct 2026 before CRC maltose/lactose rows entered; a_OH = {A_OH} J/mol, n_OH = 8, no refit\n")
        w = csv.writer(fh); w.writerow(["solute", "wt_percent", "T_perOH_C", "T_ideal_C"])
        for s, (M, n) in MOL.items():
            for wt in GRID:
                w.writerow([s, wt, f"{Tf(wt, M, n, A_OH):.3f}", f"{Tf(wt, M, n, 0.0):.3f}"])
    print(open(ROOT / "predictions" / "disaccharides_v1.csv").read()[:600])
