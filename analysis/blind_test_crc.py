"""Blind test: ethylene glycol and D-mannitol (CRC 2005) vs predictions committed in
predictions/crc_extra_v1.csv BEFORE the data were viewed. No refitting.

At CRC compositions that are not on the committed grid (mannitol 11, 13, 15 wt%; EG 44-60 wt%)
the models are EVALUATED with the same committed a_OH; these are flagged 'off-grid'.

    python analysis/blind_test_crc.py -> predictions/blind_test_crc_results.csv, figures/blind_test_crc.png
"""
import csv, math, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from predict_crc_extra import T_model, A_OH
from unifac_compare import T_freeze
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]

committed = {(r["solute"], float(r["wt_percent"])): r for r in list(csv.reader(open(ROOT / "predictions" / "crc_extra_v1.csv")))[1:] and
             csv.DictReader(open(ROOT / "predictions" / "crc_extra_v1.csv").readlines()[1:])}
data = [(r["solute"], float(r["wt_percent"]), float(r["T_freeze_C"])) for r in csv.DictReader(open(ROOT / "data" / "literature.csv"))
        if r["solute"] in ("ethylene_glycol", "mannitol")]
MODELS = ["ours", "ideal", "UNIFAC", "Dortmund"]
rows = []
for sol, w, T in data:
    key = (sol, w)
    if key in committed:
        c = committed[key]; p = [float(c["T_ours_C"]), float(c["T_ideal_C"]), float(c["T_UNIFAC_C"]), float(c["T_Dortmund_C"])]; src = "committed"
    else:
        wt = {sol: w}; p = [T_model(wt, A_OH), T_model(wt, 0.0), T_freeze(wt, "UNIFAC"), T_freeze(wt, "UNIFAC-Dortmund")]; src = "off-grid"
    rows.append((sol, w, T, p, src))

def stats(sel):
    out = {}
    for i, m in enumerate(MODELS):
        e = [p[i] - T for _, _, T, p, _ in sel]
        out[m] = (math.sqrt(sum(x * x for x in e) / len(e)), sum(e) / len(e), max(e, key=abs))
    return out

summary = []
print(f"a_OH = {A_OH} J/mol per OH (fitted on glycerol only)")
for sol, cap, lab in [("ethylene_glycol", 40, "EG <= 40 wt%"), ("ethylene_glycol", 60, "EG all (to 60 wt%)"), ("mannitol", 15, "mannitol all (to 15 wt%)")]:
    sel = [r for r in rows if r[0] == sol and r[1] <= cap and r[1] >= 2]
    st = stats(sel)
    print(f"\n{lab}, n={len(sel)} (points >= 2 wt%)")
    for m in MODELS:
        print(f"  {m:9} RMSE {st[m][0]:.2f} K   mean signed error {st[m][1]:+.2f} K   worst {st[m][2]:+.2f} K")
        summary.append([lab, len(sel), m, f"{st[m][0]:.3f}", f"{st[m][1]:+.3f}", f"{st[m][2]:+.3f}"])
with open(ROOT / "predictions" / "blind_test_crc_results.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["subset", "n", "model", "RMSE_K", "mean_signed_error_K", "worst_error_K"]); w.writerows(summary)
    w.writerow([]); w.writerow(["solute", "wt_percent", "T_CRC_C", "T_ours", "T_ideal", "T_UNIFAC", "T_Dortmund", "prediction_source"])
    for sol, wt, T, p, src in rows: w.writerow([sol, wt, T] + [f"{x:.3f}" for x in p] + [src])

plt.rcParams.update({"font.size": 8})
fig, axs = plt.subplots(1, 2, figsize=(8, 3.3))
styles = {"ours": ("k-", None), "ideal": ("k:", None), "UNIFAC": ("k--", None), "Dortmund": ("-.", "0.45")}
for ax, sol, title in zip(axs, ["ethylene_glycol", "mannitol"], ["Ethylene glycol (2 OH) - blind", "D-Mannitol (6 OH) - blind"]):
    sel = sorted([r for r in rows if r[0] == sol], key=lambda r: r[1])
    ws = [r[1] for r in sel]
    for i, m in enumerate(MODELS):
        ls, col = styles[m]
        ax.plot(ws, [r[3][i] for r in sel], ls if col is None else "-.", color=col or "black", lw=1.2 if m == "ours" else 1, label=m if m != "ours" else "per-OH model (this work)")
    ax.plot(ws, [r[2] for r in sel], "o", mfc="white", mec="black", ms=3.5, label="CRC 2005 data")
    ax.set(title=title, xlabel="wt %"); ax.spines[["top", "right"]].set_visible(False)
axs[0].set_ylabel("Freezing point / °C"); axs[0].legend(frameon=False, fontsize=7, loc="lower left")
fig.tight_layout(); fig.savefig(ROOT / "figures" / "blind_test_crc.png", dpi=170)
