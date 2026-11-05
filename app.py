"""Interactive version of the cell freeze-concentration calculation (report Section 3.3, eqs. 1-8).

Uses the same code as the report (analysis/freeze_concentration_cells.state), so the numbers match it exactly.
    streamlit run app.py
"""
import sys, warnings
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import streamlit as st

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from thermo import freezing_point
from analysis.freeze_concentration_cells import state, MODEL, SALT_LIMIT, UNFROZEN_MIN
from analysis import style as S

warnings.filterwarnings("ignore")
st.set_page_config(page_title="The Ice Decides", layout="wide")

st.title("The Ice Decides")
st.write("Pick a cryoprotectant recipe in 0.9 wt% saline and a temperature. The app solves solid–liquid "
         "equilibrium (one constant per OH group, fitted on glycerol only) and shows how much water is still liquid, "
         "how salty it is outside the cell, and how salty it gets inside the cell. Glycerol enters cells, sorbitol does not.")

c1, c2, c3 = st.columns(3)
g = c1.slider("Glycerol (wt%)", 0.0, 20.0, 6.75, 0.25)
s = c2.slider("Sorbitol (wt%)", 0.0, 15.0, 0.0, 0.25)
T = c3.slider("Temperature (°C)", -30.0, -1.0, -15.0, 0.5)

comp = {k: v / 100 for k, v in (("glycerol", g), ("sorbitol", s)) if v > 0}


@st.cache_data
def tf(g, s):
    mf = dict({k: v / 100 for k, v in (("glycerol", g), ("sorbitol", s)) if v > 0}, NaCl=0.009)
    return freezing_point(mf, MODEL, True) - 273.15


@st.cache_data
def curve(g, s, n=50):
    c = {k: v / 100 for k, v in (("glycerol", g), ("sorbitol", s)) if v > 0}
    Ts = np.linspace(-1.0, -30.0, n)
    return Ts, np.array([state(c, t)[:3] for t in Ts])


uf, out, inside, _ = state(comp, T)
m1, m2, m3, m4 = st.columns(4)
m1.metric("Medium freezes at", f"{tf(g, s):.2f} °C")
m2.metric("Water still liquid", f"{100 * uf:.1f} %", "ok" if uf >= UNFROZEN_MIN else "below 25 %",
          delta_color="normal" if uf >= UNFROZEN_MIN else "inverse", delta_arrow="off")
m3.metric("NaCl outside the cell", f"{out:.2f} wt%", "ok" if out <= SALT_LIMIT else "above 4.41 wt%",
          delta_color="normal" if out <= SALT_LIMIT else "inverse", delta_arrow="off")
m4.metric("Salt inside the cell", f"{inside:.2f} wt%", "ok" if inside <= SALT_LIMIT else "above 4.41 wt%",
          delta_color="normal" if inside <= SALT_LIMIT else "inverse", delta_arrow="off")

S.setup(8)
Ts, mine = curve(g, s)
_, ref = curve(0.0, 0.0)
fig, ax = plt.subplots(1, 3, figsize=(11, 3.1))
titles = ["Water still liquid (fraction)", "NaCl outside the cell (wt%)", "Salt inside the cell (wt%)"]
for k in range(3):
    ax[k].plot(Ts, ref[:, k], ":", color=S.DATA, lw=1.3, label="saline only")
    ax[k].plot(Ts, mine[:, k], "-", color=S.OURS, lw=1.8, label="your recipe")
    ax[k].axvline(T, color="0.6", lw=0.7, ls="--")
    ax[k].set(title=titles[k], xlabel="Temperature (°C)")
    ax[k].invert_xaxis()
    ax[k].spines[["top", "right"]].set_visible(False)
ax[0].axhline(UNFROZEN_MIN, color=S.ALT, lw=0.8)
for k in (1, 2):
    ax[k].axhline(SALT_LIMIT, color=S.ALT, lw=0.8)
    ax[k].set_ylim(0, 25)
ax[0].legend(frameon=False)
fig.tight_layout()
st.pyplot(fig)

st.caption("Orange lines: the three injury limits used in the report (4.41 wt% NaCl outside or inside, Lovelock 1953; "
           "25 % water unfrozen, Mazur et al. 1981). Equilibrium only (slow cooling); the cell's own solutes are "
           "represented as isotonic NaCl. Model, data and tests: "
           "https://github.com/Siddhant-Chowdhary/Thermodynamics-Project-")
