"""Figure: ideal vs per-OH model vs UNIFAC vs Dortmund, against CRC glycerol, and
the blind predictions for PG and sorbitol. Black-and-white.

    python analysis/plot_model_comparison.py -> figures/model_comparison.png
"""
import csv, math, sys
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parent))
from unifac_compare import T_freeze, ln_aw_ice, mole_fracs, R  # pip `thermo` inside

ROOT = Path(__file__).resolve().parents[1]
pred_hdr = open(ROOT / "predictions" / "predictions_v1.csv").readline()
A_OH = float(pred_hdr.split("a_OH = ")[1].split()[0])
NOH = {"glycerol": 3, "propylene_glycol": 2, "sorbitol": 6}


def T_model(wt, a_oh):
    names, xs = mole_fracs(wt)
    def f(T):
        alph = {k: NOH[k] * a_oh / (R * T) for k in names[1:]}
        lgw = (1 - xs[0]) * sum(alph[k] * x for k, x in zip(names[1:], xs[1:]))
        return math.log(xs[0]) + lgw - ln_aw_ice(T)
    return brentq(f, 180.0, 273.15 - 1e-9) - 273.15

lit = [(float(r["wt_percent"]), float(r["T_freeze_C"])) for r in csv.DictReader(open(ROOT / "data" / "literature.csv"))
       if r["solute"] == "glycerol" and r["source_id"] == "CRC" and float(r["wt_percent"]) >= 2]
w = np.linspace(1, 40, 40)
plt.rcParams.update({"font.size": 8})
fig, axs = plt.subplots(1, 3, figsize=(10, 3.3), sharey=True)
for ax, sol, title in zip(axs, ["glycerol", "propylene_glycol", "sorbitol"],
                          ["Glycerol (fitted, CRC data)", "Propylene glycol (blind)", "Sorbitol (blind)"]):
    ax.plot(w, [T_model({sol: x}, 0.0) for x in w], "k:", lw=1, label="ideal")
    ax.plot(w, [T_model({sol: x}, A_OH) for x in w], "k-", lw=1.4, label="per-OH model (this work)")
    ax.plot(w, [T_freeze({sol: x}, "UNIFAC") for x in w], "k--", lw=1, label="UNIFAC")
    ax.plot(w, [T_freeze({sol: x}, "UNIFAC-Dortmund") for x in w], "-.", color="0.45", lw=1, label="UNIFAC-Dortmund")
    if sol == "glycerol":
        ax.plot(*zip(*lit), "o", mfc="white", mec="black", ms=3.5, label="CRC data")
    ax.set(title=title, xlabel="wt %")
    ax.spines[["top", "right"]].set_visible(False)
axs[0].set_ylabel("Freezing point / °C")
axs[0].legend(frameon=False, fontsize=7, loc="lower left")
fig.tight_layout(); fig.savefig(ROOT / "figures" / "model_comparison.png", dpi=170)
print("saved")
