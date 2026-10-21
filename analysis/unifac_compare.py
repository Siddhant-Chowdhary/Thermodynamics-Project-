"""Freezing points from original UNIFAC and Dortmund-modified UNIFAC (open-source
`thermo` package), solved with the same SLE equation used in our model.

UNIFAC parameters come from VLE/LLE databanks, never from freezing data, so this
is an independent, external prediction that exists even with no lab work.

Run from the repo root:  python analysis/unifac_compare.py
(The pip package `thermo` shares a name with our local package; this script
deliberately does not import the local one.)
"""
import csv, math
from pathlib import Path
from scipy.optimize import brentq
from thermo.unifac import UNIFAC, DOUFSG, DOUFIP2016, UFSG, UFIP

ROOT = Path(__file__).resolve().parents[1]
R, TM, DH, DCP = 8.314462618, 273.15, 6009.5, 37.9
MW = {"ethylene_glycol": 62.068, "mannitol": 182.172, "water": 18.015, "glycerol": 92.094, "propylene_glycol": 76.094, "sorbitol": 182.172}
# subgroup counts: original UNIFAC (CH3=1, CH2=2, CH=3, OH=14, H2O=16)
ORIG = {"ethylene_glycol": {2: 2, 14: 2}, "mannitol": {2: 2, 3: 4, 14: 6}, "water": {16: 1}, "glycerol": {2: 2, 3: 1, 14: 3},
        "propylene_glycol": {1: 1, 2: 1, 3: 1, 14: 2}, "sorbitol": {2: 2, 3: 4, 14: 6}}
# Dortmund: primary OH(P)=14, secondary OH(S)=81
DORT = {"ethylene_glycol": {2: 2, 14: 2}, "mannitol": {2: 2, 3: 4, 14: 2, 81: 4}, "water": {16: 1}, "glycerol": {2: 2, 3: 1, 14: 2, 81: 1},
        "propylene_glycol": {1: 1, 2: 1, 3: 1, 14: 1, 81: 1}, "sorbitol": {2: 2, 3: 4, 14: 2, 81: 4}}
VERSIONS = {"UNIFAC": (ORIG, 0, UFIP, UFSG), "UNIFAC-Dortmund": (DORT, 1, DOUFIP2016, DOUFSG)}


def ln_aw_ice(T):
    return -(DH / R) * (1 / T - 1 / TM) - (DCP / R) * (TM / T - 1 + math.log(T / TM))


def mole_fracs(wt):  # wt: {solute: wt%}
    n = {"water": (100 - sum(wt.values())) / MW["water"]}
    n.update({k: v / MW[k] for k, v in wt.items()})
    tot = sum(n.values())
    return list(n), [v / tot for v in n.values()]


def T_freeze(wt, version):
    groups, ver, ip, sg = VERSIONS[version]
    names, xs = mole_fracs(wt)
    cg = [groups[k] for k in names]

    def f(T):
        g = UNIFAC.from_subgroups(T=T, xs=xs, chemgroups=cg, version=ver, interaction_data=ip, subgroups=sg)
        return math.log(xs[0] * g.gammas()[0]) - ln_aw_ice(T)
    return brentq(f, 200.0, TM - 1e-9) - 273.15


if __name__ == "__main__":
    lit = [r for r in csv.DictReader(open(ROOT / "data" / "literature.csv"))
           if r["solute"] == "glycerol" and r["source_id"] == "CRC" and float(r["wt_percent"]) >= 2]
    print("Validation on CRC glycerol (where data exist):")
    for v in VERSIONS:
        res = [T_freeze({"glycerol": float(r["wt_percent"])}, v) - float(r["T_freeze_C"]) for r in lit]
        print(f"  {v:16} RMSE {math.sqrt(sum(e*e for e in res)/len(res)):.2f} K, at 40 wt% error {res[-1]:+.2f} K")
    cases = [("propylene_glycol", {"propylene_glycol": w}) for w in (10, 20, 30)] + \
            [("sorbitol", {"sorbitol": w}) for w in (10, 20, 30, 40)] + \
            [("glycerol", {"glycerol": w}) for w in (10, 20, 30)] + \
            [("gly:sorb 1:1", {"glycerol": 10, "sorbitol": 10}), ("gly:sorb 3:1", {"glycerol": 15, "sorbitol": 5})]
    rows = []
    for label, wt in cases:
        a, b = T_freeze(wt, "UNIFAC"), T_freeze(wt, "UNIFAC-Dortmund")
        rows.append([label, ";".join(f"{k}={v}" for k, v in wt.items()), f"{a:.2f}", f"{b:.2f}"])
        print(f"{label:16} {str(wt):40} UNIFAC {a:7.2f}  Dortmund {b:7.2f}")
    with open(ROOT / "predictions" / "unifac_v1.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["case", "composition_wt_percent", "T_UNIFAC_C", "T_UNIFAC_Dortmund_C"])
        w.writerows(rows)
