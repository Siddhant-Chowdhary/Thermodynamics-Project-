"""Shepard, Goldston & Cocks (1976) Fig. 3, reproduced from the original scan, with this work's blind predictions
(per-OH constant fitted on glycerol only, additive ln a_w, fitted NaCl phi; nothing refitted) drawn on top.
Axis calibration from the tick marks of the 300 dpi scan (data/digitised/shepard1976_fig3_scan.png):
x: 90 wt% H2O at 458 px, 10 wt% at 1369.5 px; y: -10 C at 397 px, -50 C at 850.5 px.

    python analysis/shepard_overlay.py -> figures/shepard_overlay.pdf/.png"""
import sys, warnings
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from thermo import GroupContributionModel, freezing_point
from thermo.salt import AdditiveSaltModel
from analysis import style as S
warnings.filterwarnings("ignore")

A_OH = float(open(ROOT / "predictions" / "predictions_v1.csv").readline().split("a_OH = ")[1].split()[0])
MODEL = AdditiveSaltModel(GroupContributionModel(A_OH))
img = plt.imread(ROOT / "data" / "digitised" / "shepard1976_fig3_scan.png")
H, W = img.shape[:2]
x90, x10 = 458.0, 1369.5
y10, y50 = 397.0, 850.5
w_of = lambda px: 90 - (px - x90) * 80 / (x10 - x90)
T_of = lambda py: -10 - (py - y10) * 40 / (y50 - y10)
extent = [w_of(0), w_of(W), T_of(H), T_of(0)]

# end of each isopleth = last composition on Shepard's curve (from the hand-digitised file)
ENDS = {"0": 77.6, "20/80": 73.8, "30/70": 70.8, "40/60": 68.1, "50/50": 58.2, "70/30": 48.1, "80/20": 37.7,
        "90/10": 39.7, "inf": 33.0}


def isopleth(R, w_end):
    fg = 1.0 if R == "inf" else (0.0 if R == "0" else float(R.split("/")[0]) / 100)
    ws = np.linspace(99.0, w_end, 40)
    T = []
    for w in ws:
        tot = (100 - w) / 100
        comp = {k: v for k, v in (("glycerol", tot * fg), ("NaCl", tot * (1 - fg))) if v > 0}
        T.append(freezing_point(comp, MODEL, True) - 273.15)
    return ws, np.array(T)


S.setup(6.5)
fig, ax = plt.subplots(figsize=(3.35, 2.45))
ax.imshow(img, extent=extent, cmap="gray", aspect="auto", interpolation="lanczos")
for i, (R, w_end) in enumerate(ENDS.items()):
    w, T = isopleth(R, w_end)
    ax.plot(w, T, color=S.OURS, lw=1.0, ls=(0, (4, 1.5)), alpha=0.9, label="this work, blind (no refit)" if i == 0 else None)
ax.plot([], [], color="k", lw=1.0, label="Shepard et al. 1976 (original)")
ax.set_xlim(w_of(195), w_of(1640)); ax.set_ylim(T_of(1050), T_of(150))
ax.axis("off")
ax.legend(loc="lower left", fontsize=5.6, frameon=True, framealpha=0.95, edgecolor="0.8", bbox_to_anchor=(0.53, 0.10))
fig.tight_layout(pad=0.05)
fig.savefig(ROOT / "figures" / "shepard_overlay.pdf", dpi=300)
fig.savefig(ROOT / "figures" / "shepard_overlay.png", dpi=300)
