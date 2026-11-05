"""Fig. 1 schematic (colour): freeze concentration around a cell. Illustration only, no data.
    python analysis/schematic.py -> figures/fig_freeze_concentration.pdf/.png"""
import sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle, Ellipse
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from analysis import style as S

S.setup(8)
LIQ, ICE_F, ICE_E, CELL_F, CELL_E, SALT = "#EAF4FB", "#FFFFFF", S.SKY, "#FBE6EF", S.KM, S.ALT
rng = np.random.default_rng(3)


def hexagon(ax, x, y, r, **kw):
    ang = np.linspace(0, 2 * np.pi, 7)[:-1] + np.pi / 6
    ax.add_patch(Polygon(np.c_[x + r * np.cos(ang), y + r * np.sin(ang)], closed=True, **kw))


fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.75))
stages = [("Above freezing", 0.0, 1.00), ("Ice starts to form", 0.45, 0.80), ("Colder: less liquid, saltier", 0.75, 0.62)]
for ax, (title, ice, cell) in zip(axes, stages):
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.set_aspect("equal"); ax.axis("off")
    ax.add_patch(Rectangle((0.2, 0.2), 9.6, 9.6, fc=LIQ, ec="#33475B", lw=1))
    pts = [(x, y) for x in np.linspace(1.2, 8.8, 5) for y in np.linspace(1.2, 8.8, 5)]
    order = rng.permutation(len(pts))
    ice_pts = [pts[i] for i in order[:int(len(pts) * ice)] if np.hypot(pts[i][0] - 5, pts[i][1] - 5) >= 2.4]
    for (x, y) in ice_pts:
        hexagon(ax, x, y, 0.95, fc=ICE_F, ec=ICE_E, lw=1.0)
    keep = []
    for d in rng.uniform(0.6, 9.4, (500, 2)):
        if any(np.hypot(*(d - np.array(p))) < 0.95 for p in ice_pts) or np.hypot(*(d - 5)) < 2.2 * cell + 0.15:
            continue
        keep.append(d)
        if len(keep) == 60:
            break
    keep = np.array(keep)
    ax.plot(keep[:, 0], keep[:, 1], "o", ms=2.4, color=SALT)
    ax.add_patch(Ellipse((5, 5), 4.4 * cell, 3.8 * cell, fc=CELL_F, ec=CELL_E, lw=1.5))
    ax.text(5, 5, "cell", ha="center", va="center", fontsize=8, color=CELL_E, weight="bold")
    ax.set_title(title, fontsize=8.5)
for i in (0, 1):
    fig.patches.append(matplotlib.patches.FancyArrowPatch((0.325 + 0.335 * i, 0.5), (0.345 + 0.335 * i, 0.5),
                       transform=fig.transFigure, arrowstyle="-|>", mutation_scale=10, color=S.OURS))
fig.text(0.5, 0.03, "White hexagons: ice (pure water).   Orange dots: dissolved salt.   Blue: unfrozen liquid.   "
         "The same salt is squeezed into less liquid, and the cell loses water.", ha="center", fontsize=6.8, color="#33475B")
fig.tight_layout(rect=(0, 0.07, 1, 0.97))
fig.savefig(ROOT / "figures" / "fig_freeze_concentration.pdf")
fig.savefig(ROOT / "figures" / "fig_freeze_concentration.png", dpi=300)
