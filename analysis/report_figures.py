"""Figures for the 5-page report, drawn only from files already in predictions/ (no new modelling).

    python analysis/report_figures.py -> figures/report_pipeline.png, figures/report_results.png
"""
import csv
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from analysis import style as S

ROOT = Path(__file__).resolve().parents[1]
P = ROOT / "predictions"
S.setup()


def rows(name, skip_hash=True):
    lines = [l for l in open(P / name) if not (skip_hash and l.startswith("#"))]
    return lines


# ------------------------------------------------------------------ pipeline (block diagram)
fig, ax = plt.subplots(figsize=(7.2, 1.55))
ax.set_xlim(0, 100); ax.set_ylim(-1.5, 20); ax.axis("off")
boxes = [
    (1, "Literature data\nCRC glycerol\nfreezing points", "#E8EEF4"),
    (21, "Fit ONE constant\n$a_{OH}$ (per –OH)\nSLE + Margules", "#CFE3F2"),
    (41, "Predictions frozen\nin git before\ntest data seen", "#FCE9C4"),
    (61, "Tests on new data:\n13 solutes, 25 °C,\nglycerol + NaCl", "#D3EFE5"),
    (81, "Design: salt in/out\nof cell, glycerol-\nsparing recipe", "#F3DCEA"),
]
for x, t, c in boxes:
    ax.add_patch(FancyBboxPatch((x, 5), 17, 13, boxstyle="round,pad=0.3,rounding_size=0.8", fc=c, ec="#33475B", lw=0.9))
    ax.text(x + 8.5, 11.5, t, ha="center", va="center", fontsize=7)
for x in (18.6, 38.6, 58.6, 78.6):
    ax.annotate("", xy=(x + 2.2, 11.5), xytext=(x, 11.5), arrowprops=dict(arrowstyle="-|>", lw=1.0, color=S.OURS))
ax.text(69.5, 1.2, "vs ideal, UNIFAC, osmotic\nvirial, Kleinhans–Mazur", ha="center", fontsize=6, style="italic", color="#33475B")
ax.text(89.5, 2.2, "lab check (planned)", ha="center", fontsize=6.3, style="italic", color="#33475B")
fig.tight_layout(pad=0.1); fig.savefig(ROOT / "figures" / "report_pipeline.png", dpi=300); fig.savefig(ROOT / "figures" / "report_pipeline.pdf"); plt.close(fig)

# ------------------------------------------------------------------ results: 3 panels
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.45), gridspec_kw={"width_ratios": [1.35, 1, 1]})

# (a) blind binary tests, RMSE per molecule
rm = {}
for f, key in (("e1_blind_test_results.csv", "solute"), ("blind_test_crc_results.csv", "subset"),
               ("blind_test_sorbitol_pg.csv", "test")):
    lines = rows(f)
    head = lines.index(next(l for l in lines if l.startswith(key)))
    end = next((i for i, l in enumerate(lines[head + 1:], head + 1) if not l.strip()), len(lines))
    for r in csv.DictReader(lines[head:end]):
        val = r.get("RMSE_K") or r.get("RMSE_or_diff_K")
        rm[(r[key], r["model"])] = float(val)
mol = [("EG ≤40 wt%", "EG <= 40 wt%"), ("sucrose", "sucrose"), ("glucose", "glucose"), ("fructose", "fructose"),
       ("mannitol", "mannitol all (to 15 wt%)"), ("2-PrOH", "2-propanol"), ("1-PrOH", "1-propanol"),
       ("MeOH", "methanol"), ("EtOH", "ethanol"), ("sorbitol (conc.)", "sorbitol, Siniti 1999 without suspect x=0.05 point")]
x = np.arange(len(mol))
for model, mk, off in (("ours", "o", -0.18), ("ideal", "s", 0.0), ("UNIFAC", "^", 0.18)):
    ys = [rm.get((k, model), np.nan) for _, k in mol]
    ax[0].plot(x + off, ys, "k" + mk, ms=3.6, mfc="k" if model == "ours" else "white", ls="none", mew=0.7,
               label={"ours": "per-OH (this work)", "ideal": "ideal", "UNIFAC": "UNIFAC"}[model])
ax[0].set_yscale("log"); ax[0].set_ylim(0.01, 30)
ax[0].set_xticks(x); ax[0].set_xticklabels([m for m, _ in mol], rotation=55, ha="right", fontsize=6.5)
ax[0].axvline(8.5, color="k", lw=0.5, ls=":")
ax[0].set_ylabel("RMSE in freezing point (K)")
ax[0].set_title("(a) Blind tests on single solutes", fontsize=7.5)
ax[0].legend(fontsize=6, frameon=False, loc="upper left")

# (b) glycerol + NaCl ternary vs Shepard 1976: error vs total solute
ov = [r for r in csv.DictReader(rows("osmotic_virial_benchmark.csv")) if r["set"] == "ternary_shepard"]
tot = np.array([float(r["wt_or_glycerol_wt"]) + float(r["NaCl_wt"]) for r in ov])
for col, mk, lab in (("T_perOH_C", "o", "per-OH, additive ln $a_w$"), ("T_osmotic_virial_C", "s", "osmotic virial"),
                     ("T_ideal_or_KM_C", "^", "Kleinhans–Mazur")):
    e = np.array([float(r[col]) for r in ov]) - np.array([float(r["T_data_C"]) for r in ov])
    ax[1].plot(tot, e, "k" + mk, ms=2.8, mfc="k" if mk == "o" else "white", mew=0.6, ls="none", label=lab)
ax[1].axhline(0, color="k", lw=0.5); ax[1].axvline(30, color="k", lw=0.5, ls=":")
ax[1].set(xlabel="Total solute (wt%)", ylabel="Model − Shepard 1976 (K)")
ax[1].set_title("(b) Glycerol + NaCl mixtures", fontsize=7.5)
ax[1].legend(fontsize=6, frameon=False, loc="lower left")

# (c) glycerol-sparing curve at -15 C
rec = list(csv.DictReader(rows("recipe_v2_glycerol_sparing.csv")))
for model, ls, lab in (("per-OH", "k-", "per-OH"), ("ideal", "k:", "ideal"), ("UNIFAC", "k--", "UNIFAC")):
    s = [r for r in rec if r["polyol_model"] == model and float(r["T_design_C"]) == -15.0]
    ax[2].plot([float(r["sorbitol_wt"]) for r in s], [float(r["min_glycerol_wt"]) for r in s], ls, lw=1.1, label=lab)
ax[2].axhspan(-0.3, 0.25, color="0.85", lw=0)
ax[2].text(0.2, 0.5, "grey band: zero glycerol\n(unacceptable)", fontsize=5.5, color="0.25")
ax[2].set(xlabel="Sorbitol added (wt%)", ylabel="Minimum glycerol (wt%)", ylim=(-0.3, 9))
ax[2].set_title("(c) Recipe at −15 °C", fontsize=7.5)
ax[2].legend(fontsize=6, frameon=False, loc="upper right")
for a in ax:
    a.spines[["top", "right"]].set_visible(False)
fig.tight_layout(pad=0.3, w_pad=0.8)
fig.savefig(ROOT / "figures" / "report_results.png", dpi=300)
fig.savefig(ROOT / "figures" / "report_results.pdf")

# ------------------------------------------------------------------ results v2 (final report, colour): (c) = post-hoc constants
# maltose (blind test #3) added to panel (a): per-OH and ideal only (no UNIFAC groups for it here, as for sucrose)
for r in csv.DictReader(rows("disaccharides_blind_results.csv")):
    if r["solute"] == "maltose":
        rm[("maltose", {"per-OH": "ours", "ideal": "ideal"}[r["model"]])] = float(r["RMSE_K"])
mol2 = mol[:2] + [("maltose", "maltose")] + mol[2:]
x2 = np.arange(len(mol2))
COL = {"ours": S.OURS, "ideal": S.IDEAL, "UNIFAC": S.UNIFAC}
fig2, bx = plt.subplots(1, 3, figsize=(7.2, 2.55), gridspec_kw={"width_ratios": [1.4, 1, 1]})
bx[0].axvspan(-0.5, 2.5, color=S.OURS, alpha=0.07, lw=0)
bx[0].axvspan(7.5, 10.5, color=S.ALT, alpha=0.07, lw=0)
bx[0].text(1.0, 120, "transfers", ha="center", fontsize=6, color=S.OURS, weight="bold")
bx[0].text(9.0, 120, "fails", ha="center", fontsize=6, color=S.ALT, weight="bold")
for model, mk, off in (("ours", "o", -0.2), ("ideal", "s", 0.0), ("UNIFAC", "^", 0.2)):
    ys = [rm.get((k, model), np.nan) for _, k in mol2]
    bx[0].plot(x2 + off, ys, mk, color=COL[model], ms=4.2 if model == "ours" else 3.6, mec="white" if model == "ours" else COL[model],
               mfc=COL[model], mew=0.4, ls="none", label={"ours": "per-OH (this work)", "ideal": "ideal", "UNIFAC": "UNIFAC"}[model])
bx[0].set_yscale("log"); bx[0].set_ylim(0.01, 300)
bx[0].set_xticks(x2); bx[0].set_xticklabels([m for m, _ in mol2], rotation=55, ha="right", fontsize=6.5)
bx[0].set_ylabel("RMSE in freezing point (K)")
bx[0].set_title("(a) Blind tests, single solutes", fontsize=7.5)
bx[0].legend(fontsize=6, loc="upper center", bbox_to_anchor=(0.5, 1.0), ncol=1)
for col, mk, lab, c in (("T_perOH_C", "o", "per-OH, additive ln $a_w$", S.OURS), ("T_osmotic_virial_C", "s", "osmotic virial", S.OV),
                        ("T_ideal_or_KM_C", "^", "Kleinhans–Mazur", S.KM)):
    e = np.array([float(r[col]) for r in ov]) - np.array([float(r["T_data_C"]) for r in ov])
    bx[1].plot(tot, e, mk, color=c, ms=3.0, mfc=c if mk == "o" else "white", mew=0.7, ls="none", label=lab,
               alpha=0.9, zorder=3 if mk == "o" else 2)
bx[1].axhline(0, color="k", lw=0.5); bx[1].axvline(30, color="0.4", lw=0.6, ls=":")
bx[1].axvspan(30, max(tot) + 3, color=S.OURS, alpha=0.06, lw=0)
bx[1].text(31, bx[1].get_ylim()[1] * 0.85 if bx[1].get_ylim()[1] > 0 else 1, "cryo range", fontsize=5.8, color=S.OURS)
bx[1].set(xlabel="Total solute (wt%)", ylabel="Model − Shepard 1976 (K)")
bx[1].set_title("(b) Blind test, glycerol + NaCl", fontsize=7.5)
bx[1].legend(fontsize=6, loc="lower left")
ph = [r for r in csv.DictReader(rows("posthoc_per_molecule.csv")) if r["analysis"] == "per_molecule"]
mk = {"acyclic polyol": "o", "sugar": "s", "mono-alcohol": "^"}
seen = set()
A_OH = float(open(P / "predictions_v1.csv").readline().split("a_OH = ")[1].split()[0])
for r in ph:
    a = float(r["a_le20wt_or_a2"]) if r["solute_or_window"] != "sorbitol" else float(r["a_all_or_a1"])
    weak = float(r["max_depression_K_or_s"]) < 3.0
    c = S.CLASS[r["class"]]
    lab = r["class"] if r["class"] not in seen else None; seen.add(r["class"])
    bx[2].plot(int(r["n_OH"]), a / 1000, mk[r["class"]], color=c, ms=4.5, mfc="white" if weak else c, mew=0.9, ls="none", label=lab)
short = {"ethylene_glycol": "EG", "glycerol": "gly", "mannitol": "man", "sorbitol": "sorb", "glucose": "glc", "fructose": "fru",
         "sucrose": "suc", "methanol": "MeOH", "ethanol": "EtOH", "1-propanol": "1-PrOH", "2-propanol": "2-PrOH"}
for r in ph:
    a = float(r["a_le20wt_or_a2"]) if r["solute_or_window"] != "sorbitol" else float(r["a_all_or_a1"])
    off = {"ethylene_glycol": (-4, -9), "glycerol": (3, 3), "1-propanol": (3, 3)}.get(r["solute_or_window"], (3, -2))
    bx[2].annotate(short[r["solute_or_window"]], (int(r["n_OH"]), a / 1000), xytext=off, textcoords="offset points", fontsize=5.4)
bx[2].axhline(A_OH / 1000, color=S.OURS, lw=0.8, ls="--")
bx[2].text(5.6, A_OH / 1000 - 0.45, "glycerol constant", fontsize=5.4, color=S.OURS)
bx[2].set(xlabel="Number of OH groups", ylabel="Per-molecule a (kJ mol$^{-1}$ per OH)", ylim=(-7.2, 0), xlim=(0.5, 9.2))
bx[2].set_title("(c) Post hoc: one constant each", fontsize=7.5)
bx[2].legend(fontsize=5.6, loc="lower right")
for a in bx:
    a.spines[["top", "right"]].set_visible(False)
fig2.tight_layout(pad=0.3, w_pad=0.8)
fig2.savefig(ROOT / "figures" / "report_results_v2.png", dpi=300)
fig2.savefig(ROOT / "figures" / "report_results_v2.pdf")
