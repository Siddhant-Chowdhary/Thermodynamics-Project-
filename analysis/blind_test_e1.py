"""Extension E1 blind test: CRC 2005 alcohols and sugars vs predictions/e1_alcohols_sugars_v1.csv
(committed before the data were viewed). No refitting. Off-grid compositions are evaluated with the
same committed a_OH and flagged.
    python analysis/blind_test_e1.py -> predictions/e1_blind_test_results.csv, figures/e1_blind_test.png
"""
import csv, math, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from predict_e1 import T_model, A_OH, ALCOHOLS
from unifac_compare import T_freeze
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]
lines = open(ROOT / "predictions" / "e1_alcohols_sugars_v1.csv").read().splitlines()[1:]
committed = {(r["solute"], float(r["wt_percent"])): r for r in csv.DictReader(lines)}
SOLS = ["methanol", "ethanol", "1-propanol", "2-propanol", "glucose", "fructose", "sucrose"]
data = [(r["solute"], float(r["wt_percent"]), float(r["T_freeze_C"])) for r in csv.DictReader(open(ROOT / "data" / "literature.csv"))
        if r["solute"] in SOLS and r["source_id"] in ("CRC", "CRC2005")]   # CRC only (Baer 1984 rows added later are not part of E1)
M = ["ours", "ideal", "UNIFAC", "Dortmund"]
rows = []
for sol, w, T in data:
    c = committed.get((sol, w))
    if c:
        p = [float(c["T_ours_C"]), float(c["T_ideal_C"]), float(c["T_UNIFAC_C"]), float(c["T_Dortmund_C"])]; src = "committed"
    else:
        p = [T_model(sol, w, A_OH), T_model(sol, w, 0.0)]
        p += [T_freeze({sol: w}, "UNIFAC"), T_freeze({sol: w}, "UNIFAC-Dortmund")] if sol in ALCOHOLS else [float("nan")] * 2
        src = "off-grid"
    rows.append((sol, w, T, p, src))
out = []
print(f"a_OH = {A_OH} (glycerol only). Points 2-40 wt%. RMSE K (mean signed error K)")
for sol in SOLS:
    sel = [r for r in rows if r[0] == sol and 2 <= r[1] <= 40]
    line = f"  {sol:11} n={len(sel):2}"
    for i, m in enumerate(M):
        e = [r[3][i] - r[2] for r in sel if not math.isnan(r[3][i])]
        if not e: continue
        rm, ms = math.sqrt(sum(x * x for x in e) / len(e)), sum(e) / len(e)
        line += f" | {m} {rm:5.2f} ({ms:+.2f})"; out.append([sol, len(e), m, f"{rm:.3f}", f"{ms:+.3f}"])
    print(line)
with open(ROOT / "predictions" / "e1_blind_test_results.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["solute", "n", "model", "RMSE_K", "mean_signed_error_K"]); w.writerows(out)
    w.writerow([]); w.writerow(["solute", "wt_percent", "T_CRC_C", "T_ours", "T_ideal", "T_UNIFAC", "T_Dortmund", "prediction_source"])
    for sol, wt, T, p, src in rows: w.writerow([sol, wt, T] + [f"{x:.3f}" for x in p] + [src])

plt.rcParams.update({"font.size": 7.5})
fig, axs = plt.subplots(2, 4, figsize=(10, 4.8)); axs = axs.ravel()
for ax, sol in zip(axs, SOLS):
    sel = sorted([r for r in rows if r[0] == sol and r[1] <= 40], key=lambda r: r[1]); ws = [r[1] for r in sel]
    ax.plot(ws, [r[3][1] for r in sel], "k:", lw=1, label="ideal")
    ax.plot(ws, [r[3][0] for r in sel], "k-", lw=1.2, label="per-OH (this work)")
    if sol in ALCOHOLS:
        ax.plot(ws, [r[3][2] for r in sel], "k--", lw=1, label="UNIFAC")
        ax.plot(ws, [r[3][3] for r in sel], "-.", color="0.45", lw=1, label="Dortmund")
    ax.plot(ws, [r[2] for r in sel], "o", mfc="white", mec="black", ms=3, label="CRC 2005")
    ax.set_title(sol + " (blind)"); ax.set_xlabel("wt %"); ax.spines[["top", "right"]].set_visible(False)
axs[0].set_ylabel("Freezing point / °C"); axs[4].set_ylabel("Freezing point / °C")
axs[-1].axis("off"); axs[-1].legend(*axs[0].get_legend_handles_labels(), frameon=False, loc="center")
fig.tight_layout(); fig.savefig(ROOT / "figures" / "e1_blind_test.png", dpi=170)
