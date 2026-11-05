"""Ternary test: water-glycerol-NaCl liquidus vs Shepard, Goldston & Cocks (1976), Fig. 3 (digitised).

EVALUATION ONLY. a_OH comes from the committed predictions_v1.csv header and the NaCl phi from the
CRC 2005 binary fit; nothing is refitted to Shepard data.

Models
  primary      per-OH glycerol + fitted NaCl phi, additive in ln a_w   (pre-registered 21 Oct)
  summation    Kleinhans-Mazur: sum of binary depressions at equal molality (comparison only)
  ideal_poly   ideal glycerol + fitted NaCl phi
  ideal_all    ideal glycerol + ideal NaCl (nu = 2, phi = 1)

Pre-registered criterion (21 Oct): RMS <= 0.5 K for total solute <= 30 wt% (water >= 70 wt%).

Added after seeing the digitised data (labelled post hoc in the report): a "figure noise floor".
The two edge curves are binaries with good independent data (R = 0: CRC NaCl; R = inf: CRC glycerol,
<= 40 wt%). The scatter of the digitised edge points against CRC measures how far Shepard's smoothed
curves + digitising sit from the truth. A ternary error no larger than that floor cannot be resolved
by this figure.

    python analysis/shepard_test.py -> predictions/shepard_test_results.csv, figures/shepard_test.png
"""
import csv, math, sys, warnings
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from thermo import GroupContributionModel, IdealModel, freezing_point
from thermo.salt import AdditiveSaltModel
from thermo.constants import solute

warnings.filterwarnings("ignore", message="NaCl molality")

a_oh = float(open(ROOT / "predictions" / "predictions_v1.csv").readline().split("a_OH = ")[1].split()[0])
poly = GroupContributionModel(a_oh)
salt = AdditiveSaltModel(IdealModel())                     # fitted NaCl only
MODELS = {
    "primary": AdditiveSaltModel(poly),
    "ideal_poly": AdditiveSaltModel(IdealModel()),
    "ideal_all": AdditiveSaltModel(IdealModel(), nacl_params=(0.0, 0.0, 0.0)),
}
Mg, Ms = solute("glycerol")["molar_mass"], solute("NaCl")["molar_mass"]


def Tf(comp, m):
    comp = {k: v for k, v in comp.items() if v > 0}
    return freezing_point(comp, m, True) - 273.15


def summation(wg, ws):
    ww = 1 - wg - ws
    mg, ms = wg / Mg / ww * 1000, ws / Ms / ww * 1000
    t = 0.0
    if mg > 0:
        t += Tf({"glycerol": mg * Mg / (1000 + mg * Mg)}, poly)
    if ms > 0:
        t += Tf({"NaCl": ms * Ms / (1000 + ms * Ms)}, salt)
    return t


def split(R, w_water):
    tot = (100 - w_water) / 100
    if R == "inf":
        return tot, 0.0
    num, den = (float(x) for x in R.split("/")) if "/" in R else (0.0, 1.0)
    fg = num / (num + den)
    return tot * fg, tot * (1 - fg)


# ---- load digitised points; drop both members of any same-curve pair < 0.5 wt% apart (digitising ambiguity)
pts = [r for r in csv.DictReader(l for l in open(ROOT / "data/digitised/shepard1976_fig3.csv") if not l.lstrip("\"").startswith("#"))
       if r["point_type"] == "curve"]
keep, dropped = [], []
for r in pts:
    twins = [q for q in pts if q is not r and q["curve"] == r["curve"]
             and abs(float(q["wt_water"]) - float(r["wt_water"])) < 0.5]
    (dropped if twins else keep).append(r)

# ---- CRC binaries for the noise floor
lit = list(csv.DictReader(open(ROOT / "data/literature.csv")))
def crc(sol, wmax):
    d = sorted((float(r["wt_percent"]), float(r["T_freeze_C"])) for r in lit
               if r["solute"] == sol and r["source_id"] == "CRC2005" and r["phase_branch"] == "ice")
    return np.array([p for p in d if p[0] <= wmax])
crc_nacl, crc_gly = crc("NaCl", 23.3), crc("glycerol", 40)

rows = []
for r in keep:
    ww = float(r["wt_water"]); T = float(r["T_C"])
    wg, ws = split(r["R_glycerol_over_NaCl"], ww)
    comp = {"glycerol": wg, "NaCl": ws}
    pred = {k: Tf(comp, m) for k, m in MODELS.items()}
    pred["summation"] = summation(wg, ws)
    ref = ""
    if r["curve"] == "R0" and 100 * ws <= crc_nacl[-1, 0]:
        ref = float(np.interp(100 * ws, crc_nacl[:, 0], crc_nacl[:, 1]))
    if r["curve"] == "Rinf" and 100 * wg <= crc_gly[-1, 0]:
        ref = float(np.interp(100 * wg, crc_gly[:, 0], crc_gly[:, 1]))
    rows.append(dict(curve=r["curve"], wt_water=ww, glycerol_wt=100 * wg, NaCl_wt=100 * ws, T_shepard=T,
                     crc=ref, **{f"T_{k}": v for k, v in pred.items()}))

mk = ["primary", "summation", "ideal_poly", "ideal_all"]


def stats(sel, key):
    e = np.array([r[f"T_{key}"] - r["T_shepard"] for r in sel])
    return math.sqrt((e ** 2).mean()), e.mean(), len(e)


ternary = [r for r in rows if r["curve"] not in ("R0", "Rinf")]
lo = [r for r in ternary if r["wt_water"] >= 70]
hi = [r for r in ternary if r["wt_water"] < 70]
edge = [r for r in rows if r["crc"] != ""]
floor_e = np.array([r["T_shepard"] - r["crc"] for r in edge])
floor = math.sqrt((floor_e ** 2).mean())

lines = [f"a_OH = {a_oh} J/mol (committed). Points used {len(rows)}, dropped as digitising duplicates {len(dropped)}.",
         f"Figure noise floor (digitised edge curves vs CRC, n={len(edge)}): RMS {floor:.2f} K, mean {floor_e.mean():+.2f} K",
         "RMS / mean signed error (model - Shepard), K:"]
for name, sel in [("ternary, total solute <= 30 wt% [PRE-REGISTERED]", lo),
                  ("ternary, total solute > 30 wt%", hi), ("binary edges R=0, R=inf", [r for r in rows if r['curve'] in ('R0', 'Rinf')])]:
    lines.append("  " + name)
    for k in mk:
        s = stats(sel, k)
        lines.append(f"    {k:11s} RMS {s[0]:5.2f}  mean {s[1]:+5.2f}  n={s[2]}")
for k in mk:
    e = np.array([r[f"T_{k}"] - r["T_shepard"] for r in ternary])
    lines.append(f"  within +-1.5 K (Shepard's stated precision), all ternary: {k:11s} {np.mean(abs(e) <= 1.5)*100:.0f}% of {len(e)}")
lines.append("Per curve, primary model (RMS, mean signed):")
for c in dict.fromkeys(r["curve"] for r in rows):
    s = stats([r for r in rows if r["curve"] == c], "primary")
    lines.append(f"    {c:7s} RMS {s[0]:5.2f}  mean {s[1]:+5.2f}  n={s[2]}")
txt = "\n".join(lines)
print(txt)

with open(ROOT / "predictions" / "shepard_test_results.csv", "w", newline="") as f:
    f.write("# " + txt.replace("\n", "\n# ") + "\n")
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    for r in rows:
        w.writerow({k: (f"{v:.3f}" if isinstance(v, float) else v) for k, v in r.items()})

# ---- figure: parity, b/w, marker shapes
plt.rcParams.update({"font.size": 9})
fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.9))
marks = {"primary": ("o", "per-OH + NaCl (additive)"), "summation": ("^", "K-M summation"),
         "ideal_poly": ("s", "ideal glycerol + NaCl"), "ideal_all": ("x", "ideal all")}
for k, (m, lab) in marks.items():
    kw = dict(c="k") if m == "x" else dict(facecolors="none", edgecolors="k")
    axs[0].scatter([r["T_shepard"] for r in ternary], [r[f"T_{k}"] for r in ternary], marker=m,
                   s=18, linewidths=0.7, label=lab, **kw)
lim = [-52, 0]
axs[0].plot(lim, lim, "k-", lw=0.6)
axs[0].fill_between(lim, [l - 1.5 for l in lim], [l + 1.5 for l in lim], color="0.85", lw=0, label="±1.5 K (Shepard)")
axs[0].set(xlim=lim, ylim=lim, xlabel="Shepard 1976, digitised T (°C)", ylabel="Predicted T (°C)",
           title="Ternary isopleths (R = 20/80 … 90/10)")
axs[0].legend(fontsize=7, frameon=False, loc="upper left")
ls = ["-", "--", ":", "-."]
for i, c in enumerate(["R20_80", "R50_50", "R70_30", "R90_10"]):
    sel = sorted([r for r in rows if r["curve"] == c], key=lambda r: -r["wt_water"])
    axs[1].plot([r["wt_water"] for r in sel], [r["T_shepard"] - r["T_primary"] for r in sel], "k" + ls[i],
                marker="o", ms=3, mfc="white", lw=0.8, label=c.replace("R", "R = ").replace("_", "/"))
axs[1].axhspan(-floor, floor, color="0.85", lw=0, label=f"figure noise floor ±{floor:.2f} K")
axs[1].axhline(0, color="k", lw=0.5)
axs[1].axvline(70, color="k", lw=0.5, ls=":")
axs[1].text(70.5, 1.85, "← solute ≤ 30 wt%", fontsize=7, ha="right")
axs[1].invert_xaxis()
axs[1].set(xlabel="Water (wt%)", ylabel="Shepard − per-OH model (K)", title="Error along isopleths")
axs[1].legend(fontsize=7, frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.2))
fig.tight_layout()
fig.savefig(ROOT / "figures" / "shepard_test.png", dpi=300)
