"""Blind predictions for Uraji et al., Food Sci. Technol. Int. 2, 38-42 (1996): xylitol, sorbitol, maltitol, lactitol at
30-60 wt%. Written 30 Oct 2026 BEFORE the figure is digitised (no data from that paper seen; only its abstract-level
description). Committed a_OH, no refit. UNIFAC only where group assignment is unambiguous (acyclic polyols).

    python analysis/predict_uraji.py -> predictions/uraji_polyols_v1.csv
"""
import csv, math, sys
from pathlib import Path
from scipy.optimize import brentq
sys.path.insert(0, str(Path(__file__).resolve().parent))
from unifac_compare import UFIP, UFSG, ln_aw_ice, R
from thermo.unifac import UNIFAC

ROOT = Path(__file__).resolve().parents[1]
A_OH = float(open(ROOT / "predictions" / "predictions_v1.csv").readline().split("a_OH = ")[1].split()[0])
MW_W = 18.015
SOL = {"xylitol": (152.146, 5, {2: 2, 3: 3, 14: 5}), "sorbitol": (182.172, 6, {2: 2, 3: 4, 14: 6}),
       "maltitol": (344.31, 9, None), "lactitol": (344.31, 9, None)}


def T_model(w, Mw, n, kind, groups=None):
    nw, ns = (100 - w) / MW_W, w / Mw
    x = ns / (nw + ns)
    def f(T):
        if kind == "perOH":
            lg = n * A_OH / (R * T) * x * x
        elif kind == "ideal":
            lg = 0.0
        else:
            lg = math.log(UNIFAC.from_subgroups(T=T, xs=[1 - x, x], chemgroups=[{16: 1}, groups], version=0,
                                                interaction_data=UFIP, subgroups=UFSG).gammas()[0])
        return math.log(1 - x) + lg - ln_aw_ice(T)
    return brentq(f, 150, 273.15 - 1e-9) - 273.15


rows = []
for s, (Mw, n, g) in SOL.items():
    for w in (30, 40, 50, 60):
        r = [s, w, n, f"{T_model(w, Mw, n, 'perOH'):.2f}", f"{T_model(w, Mw, n, 'ideal'):.2f}",
             f"{T_model(w, Mw, n, 'unifac', g):.2f}" if g else ""]
        rows.append(r); print(r)
with open(ROOT / "predictions" / "uraji_polyols_v1.csv", "w", newline="") as fh:
    fh.write(f"# blind predictions committed 30 Oct 2026 before digitising Uraji 1996; a_OH = {A_OH}\n")
    w = csv.writer(fh); w.writerow(["solute", "wt_percent", "n_OH", "T_perOH_C", "T_ideal_C", "T_UNIFAC_C"]); w.writerows(rows)
