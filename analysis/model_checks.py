"""Model checks on the glycerol fit (21 Oct). Post-hoc diagnostics, NOT blind tests.

(1) Fit-range sensitivity: refit a_OH on CRC glycerol 2-20, 2-30, 2-40 wt%.
(2) Model form: two-parameter per-OH Redlich-Kister/Margules,
      G_ex/RT = x_w x_s [A0 + A1 (x_w - x_s)],  A0 = n_OH a0/RT, A1 = n_OH a1/RT
      => ln gamma_w = x_s^2 [A0 + A1 (3 x_w - x_s)]   (binary)
    fitted on the same glycerol data; residuals compared with the one-parameter model.
EG numbers here are evaluated AFTER the EG data was seen, so they are diagnostics only.

    python analysis/model_checks.py -> predictions/model_checks.csv, figures/model_checks.png
"""
import csv, math, sys
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares, brentq
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from thermo.sle import ln_aw_ice
from thermo.constants import R, M_WATER, solute
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

lit = list(csv.DictReader(open(ROOT / "data" / "literature.csv")))
def series(sol):
    return [(float(r["wt_percent"]), float(r["T_freeze_C"])) for r in lit
            if r["solute"] == sol and r["source_id"] == "CRC2005" and float(r["wt_percent"]) >= 2]
gly, eg = series("glycerol"), series("ethylene_glycol")

def xs_of(sol, w):
    ns = w / solute(sol)["molar_mass"]; nw = (100 - w) / M_WATER
    return ns / (ns + nw)

def Tf(sol, w, a0, a1=0.0):
    x = xs_of(sol, w); n = solute(sol)["n_OH"]; xw = 1 - x
    f = lambda T: math.log(xw) + x * x * (n * a0 / (R * T) + n * a1 / (R * T) * (3 * xw - x)) - ln_aw_ice(T, True)
    return brentq(f, 180, 273.15 - 1e-9) - 273.15

def fit(data, two=False):
    res = lambda p: [Tf("glycerol", w, p[0], p[1] if two else 0.0) - T for w, T in data]
    s = least_squares(res, [-2400.0, 0.0] if two else [-2400.0])
    r = np.array(s.fun); return s.x, math.sqrt(np.mean(r ** 2)), r

def rmse_on(sol, data, a0, a1=0.0, wmax=40):
    e = [Tf(sol, w, a0, a1) - T for w, T in data if w <= wmax]; return math.sqrt(np.mean(np.square(e)))

rows = []
print("(1) fit-range sensitivity (one parameter)")
for wmax in (20, 30, 40):
    p, rm, _ = fit([d for d in gly if d[0] <= wmax])
    eg_r = rmse_on("ethylene_glycol", eg, p[0])
    print(f"  glycerol 2-{wmax} wt%: a_OH = {p[0]:7.1f} J/mol, fit RMSE {rm*1000:4.0f} mK, EG(2-40) RMSE {eg_r:.2f} K")
    rows.append([f"1-param, glycerol 2-{wmax}", f"{p[0]:.1f}", "", f"{rm:.4f}", f"{eg_r:.3f}"])
print("(2) model form on glycerol 2-40 wt%")
p1, rm1, r1 = fit(gly)
p2, rm2, r2 = fit(gly, two=True)
eg2 = rmse_on("ethylene_glycol", eg, p2[0], p2[1])
print(f"  1-param: a_OH = {p1[0]:.1f}; RMSE {rm1*1000:.0f} mK; residuals at 24-40 wt%: {np.round(r1[-5:],3)}")
print(f"  2-param: a0 = {p2[0]:.1f}, a1 = {p2[1]:.1f}; RMSE {rm2*1000:.0f} mK; residuals at 24-40 wt%: {np.round(r2[-5:],3)}")
print(f"  EG(2-40) RMSE with 2-param (post-hoc): {eg2:.2f} K")
rows.append(["2-param, glycerol 2-40", f"{p2[0]:.1f}", f"{p2[1]:.1f}", f"{rm2:.4f}", f"{eg2:.3f}"])
with open(ROOT / "predictions" / "model_checks.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["# post-hoc diagnostics (EG data already seen); not blind tests"])
    w.writerow(["model", "a0_J_per_mol_per_OH", "a1_J_per_mol_per_OH", "glycerol_fit_RMSE_K", "EG_2_40_RMSE_K"]); w.writerows(rows)

fig, ax = plt.subplots(figsize=(4.8, 3.2))
ws = [w for w, _ in gly]
ax.axhline(0, color="0.7", lw=0.8)
ax.plot(ws, r1 * 1000, "ko-", mfc="white", ms=3.5, lw=0.8, label="one parameter")
ax.plot(ws, r2 * 1000, "ks--", ms=3, lw=0.8, label="two parameters")
ax.set(xlabel="Glycerol / wt %", ylabel="Model − CRC / mK")
ax.spines[["top", "right"]].set_visible(False); ax.legend(frameon=False, fontsize=7)
fig.tight_layout(); fig.savefig(ROOT / "figures" / "model_checks.png", dpi=170)
