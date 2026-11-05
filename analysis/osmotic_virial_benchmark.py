"""Benchmark against the cryobiology field's own model: the Elliott et al. molality-based osmotic virial equation,
with the coefficients of Zielinski et al., Cryobiology 69, 305 (2014), Table 3.

    pi = sum_i k_i m_i + sum_i sum_j (B_ii + B_jj)/2 (k_i m_i)(k_j m_j) + (C terms),   ln a_w = -M_w pi

Per-OH and ideal use our ice-equilibrium equation (dH, dCp). The osmotic virial is converted to a freezing point with
Zielinski's OWN relation, dT = [R T_m pi M1/ds] / [1 + R pi M1/ds], ds = 22.00 J/mol/K, M1 = 0.01802 kg/mol, because its
coefficients were fitted with it (using our SLE instead makes it worse, e.g. glycerol 0.46 vs 0.27 K; the fairer choice is kept).

Fairness note (stated in the report): the osmotic virial coefficients were FITTED per solute, largely to the same CRC
freezing-point compilation [their ref. 72] that we use as test data for EG, sucrose etc. On binaries the osmotic virial is
therefore an in-sample reference (one fitted constant per molecule), while per-OH uses ONE constant fitted on glycerol only.
The fair head-to-head is the ternary glycerol + NaCl (Shepard 1976), where both models predict from binaries alone.

    python analysis/osmotic_virial_benchmark.py -> predictions/osmotic_virial_benchmark.csv, figures/osmotic_virial_benchmark.png
"""
import csv, math
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
R, MW_W, DH, DCP, TM = 8.314462618, 0.018015, 6009.5, 37.9, 273.15
A_OH = float(open(ROOT / "predictions" / "predictions_v1.csv").readline().split("a_OH = ")[1].split()[0])
MW = {"glycerol": 92.094, "ethylene_glycol": 62.068, "sucrose": 342.297, "glucose": 180.156, "methanol": 32.042,
      "ethanol": 46.069, "mannitol": 182.172, "NaCl": 58.443}
NOH = {"glycerol": 3, "ethylene_glycol": 2, "sucrose": 8, "glucose": 5, "methanol": 1, "ethanol": 1, "mannitol": 6}
# Zielinski 2014 Table 3: (k, B [kg/mol], C [kg2/mol2]); 'dextrose' = glucose
OV = {"NaCl": (1.678, 0.044, 0.0), "glycerol": (1, 0.023, 0.0), "ethylene_glycol": (1, 0.020, 0.0),
      "sucrose": (1, 0.116, 0.0), "glucose": (1, 0.044, 0.0), "methanol": (1, 0.0036, 0.0),
      "ethanol": (1, 0.012, 0.0004), "mannitol": (1, 0.0, 0.0), "PG": (1, 0.039, 0.0)}


def ln_aw_ice(T):
    return -(DH / R) * (1 / T - 1 / TM) - (DCP / R) * (TM / T - 1 + math.log(T / TM))


def ln_aw_ov(mol):          # mol: {solute: molality}
    km = {s: OV[s][0] * m for s, m in mol.items()}
    pi = sum(km.values())
    for i in km:
        for j in km:
            pi += (OV[i][1] + OV[j][1]) / 2 * km[i] * km[j]
    # third-order term: only for single-solute use here (ethanol); multi-solute C mixing not needed
    if len(km) == 1:
        (s, v), = km.items(); pi += OV[s][2] * v ** 3
    return -MW_W * pi


def ln_aw_perOH(mol, T):
    nw = 1 / MW_W; n = dict(mol); tot = nw + sum(n.values())
    xw = nw / tot
    return math.log(xw) + (1 - xw) * sum(NOH[s] * A_OH / (R * T) * (v / tot) for s, v in n.items())


def ln_aw_ideal(mol):
    nw = 1 / MW_W; return math.log(nw / (nw + sum(mol.values())))


def Tf_ov(mol):
    a = R * (-ln_aw_ov(mol) / MW_W) * 0.01802 / 22.00
    return -TM * a / (1 + a)


def Tf(fn):
    return brentq(lambda T: fn(T) - ln_aw_ice(T), 180.0, TM - 1e-9) - 273.15


lit = list(csv.DictReader(open(ROOT / "data" / "literature.csv")))
lines = ["Part 1, binaries (CRC, 2-40 wt%, ice branch). Osmotic virial is IN-SAMPLE here (fitted to CRC by Zielinski);",
         "per-OH uses one constant from glycerol. RMSE in freezing point, K:",
         "  solute | n | per-OH | osmotic virial | ideal"]
out = []
for s in ("glycerol", "ethylene_glycol", "sucrose", "glucose", "mannitol", "methanol", "ethanol"):
    rows = [r for r in lit if r["solute"] == s and r["source_id"] in ("CRC", "CRC2005") and r["phase_branch"] == "ice"
            and 2 <= float(r["wt_percent"]) <= 40]
    seen = {}
    for r in rows:      # one row per wt% (prefer CRC2005)
        if r["wt_percent"] not in seen or r["source_id"] == "CRC2005":
            seen[r["wt_percent"]] = r
    e = {"perOH": [], "OV": [], "ideal": []}
    for r in seen.values():
        w = float(r["wt_percent"]); m = w / (100 - w) * 1000 / MW[s]; T = float(r["T_freeze_C"])
        p = {"perOH": Tf(lambda TT: ln_aw_perOH({s: m}, TT)), "OV": Tf_ov({s: m}),
             "ideal": Tf(lambda TT: ln_aw_ideal({s: m}))}
        for k in e: e[k].append(p[k] - T)
        out.append(["binary", s, w, "", T] + [round(p[k], 3) for k in ("perOH", "OV", "ideal")])
    rm = {k: math.sqrt(np.mean(np.square(v))) for k, v in e.items()}
    lines.append(f"  {s:15s} | {len(seen):2d} | {rm['perOH']:.2f} | {rm['OV']:.2f} | {rm['ideal']:.2f}")

# Part 2: ternary glycerol + NaCl vs Shepard 1976 (digitised), both models predicting from binaries only
sh = list(csv.DictReader(l for l in open(ROOT / "predictions" / "shepard_test_results.csv") if not l.startswith("#")))
tern = [r for r in sh if r["curve"] not in ("R0", "Rinf")]
res = {"perOH (additive ln a_w)": [], "osmotic virial": [], "OV with our SLE": [], "K-M summation": []}
for r in tern:
    wg, wn = float(r["glycerol_wt"]), float(r["NaCl_wt"]); ww = 100 - wg - wn
    mol = {"glycerol": wg / MW["glycerol"] / ww * 1000, "NaCl": wn / MW["NaCl"] / ww * 1000}
    tov = Tf_ov(mol)
    T = float(r["T_shepard"])
    res["perOH (additive ln a_w)"].append((ww, float(r["T_primary"]) - T))
    res["osmotic virial"].append((ww, tov - T))
    res["K-M summation"].append((ww, float(r["T_summation"]) - T))
    res["OV with our SLE"].append((ww, Tf(lambda TT: ln_aw_ov(mol)) - T))
    out.append(["ternary_shepard", r["curve"], wg, wn, T, r["T_primary"], round(tov, 3), r["T_summation"]])
lines.append("Part 2, glycerol + NaCl ternary vs Shepard 1976 (digitised), predictions from binaries only. RMS (mean signed), K:")
for k, v in res.items():
    lo = [e for w, e in v if w >= 70]; hi = [e for w, e in v if w < 70]
    f = lambda a: f"{math.sqrt(np.mean(np.square(a))):.2f} ({np.mean(a):+.2f})"
    lines.append(f"  {k:24s} solute <=30 wt%: {f(lo)}   solute >30 wt%: {f(hi)}   all: {f([e for _, e in v])}")
txt = "\n".join(lines); print(txt)
with open(ROOT / "predictions" / "osmotic_virial_benchmark.csv", "w", newline="") as fh:
    fh.write("# " + txt.replace("\n", "\n# ") + "\n")
    w = csv.writer(fh)
    w.writerow(["set", "solute_or_curve", "wt_or_glycerol_wt", "NaCl_wt", "T_data_C", "T_perOH_C", "T_osmotic_virial_C", "T_ideal_or_KM_C"])
    w.writerows(out)

plt.rcParams.update({"font.size": 9})
fig, ax = plt.subplots(figsize=(3.8, 3.2))
for (k, v), mk in zip([(k, v) for k, v in res.items() if k != "OV with our SLE"], ("o", "s", "^")):
    v = sorted(v)
    ax.plot([100 - a for a, _ in v], [b for _, b in v], "k" + mk, ms=3.5, mfc="white", ls="none", label=k)
ax.axhline(0, color="k", lw=0.5); ax.axvline(30, color="k", lw=0.5, ls=":")
ax.set(xlabel="Total solute (wt%)", ylabel="Model − Shepard 1976 (K)", title="Glycerol + NaCl: three mixture models")
ax.legend(fontsize=7, frameon=False)
fig.tight_layout(); fig.savefig(ROOT / "figures" / "osmotic_virial_benchmark.png", dpi=300)
