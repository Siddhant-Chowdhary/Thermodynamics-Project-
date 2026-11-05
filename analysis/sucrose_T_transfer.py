"""Separate the temperature effect from the molecule effect using ONE molecule (sucrose, 8 OH) at two temperatures.

Freezing line (~ -0.1 to -4.9 C): CRC 2005 sucrose freezing points (data/literature.csv).
25 C: Robinson & Stokes, Trans. Faraday Soc. 45, 612 (1949), Table IV (data/robinson_stokes1949_sucrose_25C.csv).
EVALUATION ONLY: committed a_OH, A = n_OH a_OH constant in T (as fixed in PROMPTS.md, 17 Oct). No refit.
Also reports the per-point implied a_OH (a diagnostic, not a fit).
Why it matters: sorbitol fails at 25 C (analysis/sorbitol_25C.py). If sucrose works on the freezing line but fails
at 25 C, the cause is temperature; if it works at both, the sorbitol failure is a molecule effect.
    python analysis/sucrose_T_transfer.py -> predictions/sucrose_T_transfer.csv, figures/sucrose_T_transfer.png
"""
import csv, math
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
R, MW_W, MW_S, NOH = 8.314462618, 0.018015, 0.342297, 8
DH, DCP, TM = 6009.5, 37.9, 273.15
A_OH = float(open(ROOT / "predictions" / "predictions_v1.csv").readline().split("a_OH = ")[1].split()[0])

def ln_aw_ice(T):
    return -(DH / R) * (1 / T - 1 / TM) - (DCP / R) * (TM / T - 1 + math.log(T / TM))
def x_of(m): return m / (1 / MW_W + m)
def phi_model(m, T):
    x = x_of(m); return -(math.log(1 - x) + NOH * A_OH / (R * T) * x * x) / (m * MW_W)
def phi_ideal(m): return -math.log(1 - x_of(m)) / (m * MW_W)
def a_impl(m, phi, T):
    x = x_of(m); return (-phi * m * MW_W - math.log(1 - x)) / (x * x) * R * T / NOH

frz = []
for r in csv.DictReader(open(ROOT / "data" / "literature.csv")):
    if r["solute"] == "sucrose" and r["source_id"] in ("CRC", "CRC2005") and float(r["wt_percent"]) >= 4:
        w = float(r["wt_percent"]); T = float(r["T_freeze_C"]) + 273.15; m = w / (100 - w) / MW_S
        frz.append((m, -ln_aw_ice(T) / (m * MW_W), phi_model(m, T), a_impl(m, -ln_aw_ice(T) / (m * MW_W), T), T))
r25 = []
for r in csv.DictReader(l for l in open(ROOT / "data" / "robinson_stokes1949_sucrose_25C.csv") if not l.startswith("#")):
    m, phi = float(r["m"]), float(r["phi"])
    if m >= 0.3:
        r25.append((m, phi, phi_model(m, 298.15), a_impl(m, phi, 298.15), 298.15))

def summ(rows, lo, hi):
    s = [r for r in rows if lo <= r[0] <= hi]
    e = [r[2] - r[1] for r in s]
    return math.sqrt(np.mean(np.square(e))), np.median([r[3] for r in s]), np.mean([abs(r[1] - 1) for r in s]), len(s)

lines = [f"Sucrose, committed a_OH = {A_OH} J/mol/OH. Points with m >= ~0.12 (freezing) / 0.3 (25 C) to limit rounding noise.",
         "range (mol/kg) | set | RMS err in phi | median implied a_OH | mean |phi-1| | n"]
for lo, hi in ((0.1, 1.0), (1.0, 2.2)):
    for name, rows in (("freezing line", frz), ("25 C", r25)):
        s = summ(rows, lo, hi)
        lines.append(f"  {lo}-{hi} | {name:13s} | {s[0]:.3f} | {s[1]:7.0f} | {s[2]:.3f} | {s[3]}")
s = summ(r25, 2.2, 6.0)
lines.append(f"  2.2-6.0 | 25 C          | {s[0]:.3f} | {s[1]:7.0f} | {s[2]:.3f} | {s[3]}")
txt = "\n".join(lines); print(txt)
with open(ROOT / "predictions" / "sucrose_T_transfer.csv", "w", newline="") as f:
    f.write("# " + txt.replace("\n", "\n# ") + "\n")
    w = csv.writer(f); w.writerow(["set", "m", "phi_data", "phi_perOH", "phi_ideal", "a_OH_implied", "T_K"])
    for name, rows in (("freezing", frz), ("25C", r25)):
        for r in rows:
            w.writerow([name, f"{r[0]:.4f}", f"{r[1]:.4f}", f"{r[2]:.4f}", f"{phi_ideal(r[0]):.4f}", f"{r[3]:.0f}", f"{r[4]:.2f}"])

plt.rcParams.update({"font.size": 9})
fig, ax = plt.subplots(figsize=(3.6, 3.2))
mm = np.linspace(0.05, 6, 100)
ax.plot(mm, [phi_model(m, 298.15) for m in mm], "k-", lw=1, label="per-OH at 25 °C")
ax.plot(mm, [phi_model(m, 271.0) for m in mm], "k--", lw=1, label="per-OH at −2 °C")
ax.plot(mm, [phi_ideal(m) for m in mm], "k:", lw=1, label="ideal")
ax.plot([r[0] for r in r25], [r[1] for r in r25], "ks", ms=3.5, mfc="white", label="Robinson & Stokes 1949, 25 °C")
ax.plot([r[0] for r in frz], [r[1] for r in frz], "k^", ms=3.5, label="CRC freezing points (~0 °C)")
ax.set(xlabel="Sucrose molality (mol/kg)", ylabel="Osmotic coefficient φ", title="Sucrose at two temperatures")
ax.legend(fontsize=6.5, frameon=False)
fig.tight_layout(); fig.savefig(ROOT / "figures" / "sucrose_T_transfer.png", dpi=300)
