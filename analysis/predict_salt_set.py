"""Pre-registered predictions for the salt + polyol tests (written 21 Oct, BEFORE any data).

Primary model (decided 21 Oct): per-OH polyols + fitted NaCl phi, ADDITIVE IN ln a_w.
Comparison (reported, never switched to): Kleinhans-Mazur summation of binary freezing-point
depressions at equal molality; and ideal polyols + fitted NaCl.
Pass criterion for Shepard 1976 (pre-registered): RMS <= 0.5 K for total solute <= 30 wt%.

    python analysis/predict_salt_set.py -> predictions/salt_set_v1.csv
"""
import csv, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from thermo import GroupContributionModel, IdealModel, freezing_point
from thermo.salt import AdditiveSaltModel
from thermo.constants import solute

a_oh = float(open(ROOT / "predictions" / "predictions_v1.csv").readline().split("a_OH = ")[1].split()[0])
poly = GroupContributionModel(a_oh)
primary = AdditiveSaltModel(poly)
salt_only = AdditiveSaltModel(IdealModel())
ideal_poly = AdditiveSaltModel(IdealModel())


def Tf(comp, m):
    return freezing_point(comp, m, True) - 273.15


def summation(wg, ws):
    ww = 1 - wg - ws
    Mg, Ms = solute("glycerol")["molar_mass"], solute("NaCl")["molar_mass"]
    mg, ms = wg / Mg / ww * 1000, ws / Ms / ww * 1000
    return Tf({"glycerol": mg * Mg / (1000 + mg * Mg)}, poly) + Tf({"NaCl": ms * Ms / (1000 + ms * Ms)}, salt_only)


rows = []
for label, wg, ws in [("LAB session 2: glycerol 20 + NaCl 5", 0.20, 0.05),
                      ("glycerol 10 + NaCl 0.9", 0.10, 0.009),
                      ("glycerol 10 + NaCl 3", 0.10, 0.03)]:
    c = {"glycerol": wg, "NaCl": ws}
    r = [Tf(c, primary), summation(wg, ws), Tf(c, ideal_poly)]
    rows.append([label, f"{100*wg:g}", f"{100*ws:g}"] + [f"{v:.2f}" for v in r])
    print(label, [round(v, 2) for v in r])
with open(ROOT / "predictions" / "salt_set_v1.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["# primary = additive ln a_w; pass criterion for Shepard 1976: RMS <= 0.5 K for total solute <= 30 wt% (pre-registered 21 Oct)"])
    w.writerow(["sample", "glycerol_wt", "NaCl_wt", "T_primary_C", "T_summation_C", "T_ideal_polyol_C"])
    w.writerows(rows)
