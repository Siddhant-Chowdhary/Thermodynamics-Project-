"""POST HOC (done after seeing the blind-test results; not part of the blind test).

1. Per-molecule constant. For each solute, fit its own a (J/mol per OH) to its own freezing points, then ask what
   structural feature explains the spread. The blind-test constant is glycerol's a_OH = -2414.9 J/mol.
2. Gibbs-Helmholtz split. Sucrose is the one molecule with good data at two temperatures (CRC freezing points near
   -1 to -5 C; Robinson & Stokes 1949 at 25 C). Writing the per-OH interaction free energy as a(T) = h - T s gives the
   enthalpy h and entropy s of the water-OH interaction, and explains why the interaction weakens on warming.
   Fitted on matched molality windows so concentration effects do not masquerade as temperature effects.

    python analysis/posthoc_per_molecule.py -> predictions/posthoc_per_molecule.csv, figures/posthoc_per_molecule.png
"""
import csv, math
from pathlib import Path
import numpy as np
from scipy.optimize import brentq, minimize_scalar
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys as _sys
_sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from analysis import style as S

ROOT = Path(__file__).resolve().parents[1]
R, MW_W, DH, DCP, TM = 8.314462618, 0.018015, 6009.5, 37.9, 273.15
A_OH = float(open(ROOT / "predictions" / "predictions_v1.csv").readline().split("a_OH = ")[1].split()[0])
MOL = {  # molar mass g/mol, n_OH, n_C, class
    "ethylene_glycol": (62.068, 2, 2, "acyclic polyol"), "glycerol": (92.094, 3, 3, "acyclic polyol"),
    "mannitol": (182.172, 6, 6, "acyclic polyol"), "sorbitol": (182.172, 6, 6, "acyclic polyol"),
    "glucose": (180.156, 5, 6, "sugar"), "fructose": (180.156, 5, 6, "sugar"), "sucrose": (342.297, 8, 12, "sugar"),
    "methanol": (32.042, 1, 1, "mono-alcohol"), "ethanol": (46.069, 1, 2, "mono-alcohol"),
    "1-propanol": (60.096, 1, 3, "mono-alcohol"), "2-propanol": (60.096, 1, 3, "mono-alcohol"),
}


def ln_aw_ice(T):
    return -(DH / R) * (1 / T - 1 / TM) - (DCP / R) * (TM / T - 1 + math.log(T / TM))


def x_of(m):
    return m / (1 / MW_W + m)


def Tf(m, n, a):
    x = x_of(m)
    return brentq(lambda T: math.log(1 - x) + n * a / (R * T) * x * x - ln_aw_ice(T), 150, TM - 1e-9)


def a_implied(m, ln_aw, T, n):
    x = x_of(m)
    return (ln_aw - math.log(1 - x)) / (x * x) * R * T / n


lit = list(csv.DictReader(open(ROOT / "data" / "literature.csv")))
out, fits = [], {}
print("Per-molecule a (J/mol per OH), least squares on freezing points:")
for s, (Mw, n, nc, cls) in MOL.items():
    if s == "sorbitol":
        rows = [r for r in lit if r["solute"] == s and r["source_id"] == "SINITI1999" and float(r["wt_percent"]) > 40]
        lab = "Siniti 1999, 53-64 wt%"
    else:
        rows = [r for r in lit if r["solute"] == s and r["source_id"] == "CRC2005" and r["phase_branch"] == "ice"
                and 2 <= float(r["wt_percent"]) <= 40]
        lab = "CRC 2005, 2-40 wt%"
    pts = [(float(r["wt_percent"]) / (100 - float(r["wt_percent"])) * 1000 / Mw, float(r["T_freeze_C"]) + 273.15) for r in rows]

    def sse(a, P=pts):
        return sum((Tf(m, n, a) - T) ** 2 for m, T in P)
    a_all = minimize_scalar(sse, bounds=(-15000, 5000), method="bounded").x
    lo = [p for p in pts if p[0] * Mw / 1000 / (1 + p[0] * Mw / 1000) <= 0.20] if s != "sorbitol" else pts
    a_lo = minimize_scalar(lambda a: sse(a, lo), bounds=(-15000, 5000), method="bounded").x if len(lo) >= 3 else float("nan")
    rmse = math.sqrt(sse(a_all) / len(pts))
    depth = max(TM - T for _, T in pts)
    fits[s] = (a_all, a_lo, n, nc, cls, rmse, depth, len(pts))
    out.append(["per_molecule", s, cls, n, nc, f"{a_all:.0f}", f"{a_lo:.0f}", f"{rmse:.3f}", f"{depth:.2f}", len(pts), lab])
    print(f"  {s:16s} {cls:15s} n_OH={n} a={a_all:7.0f} (<=20 wt%: {a_lo:7.0f})  RMSE {rmse:.3f} K  max depression {depth:5.2f} K  n={len(pts)}")

# --- Gibbs-Helmholtz on sucrose: implied a per point at two temperatures, matched molality windows
Mw, n = 342.297, 8
frz = []
for r in lit:
    if r["solute"] == "sucrose" and r["source_id"] == "CRC2005" and float(r["wt_percent"]) >= 4:
        w = float(r["wt_percent"]); T = float(r["T_freeze_C"]) + 273.15; m = w / (100 - w) * 1000 / Mw
        frz.append((m, T, a_implied(m, ln_aw_ice(T), T, n)))
r25 = []
for r in csv.DictReader(l for l in open(ROOT / "data" / "robinson_stokes1949_sucrose_25C.csv") if not l.startswith("#")):
    m, phi = float(r["m"]), float(r["phi"])
    if m >= 0.3:
        r25.append((m, 298.15, a_implied(m, -phi * m * MW_W, 298.15, n)))
print("\nGibbs-Helmholtz, sucrose, a(T) = h - T s (two temperature clusters, matched molality):")
gh = []
for lo, hi in ((0.3, 1.0), (1.0, 2.2)):
    f = [p for p in frz if lo <= p[0] <= hi]; t = [p for p in r25 if lo <= p[0] <= hi]
    T1, a1 = np.mean([p[1] for p in f]), np.median([p[2] for p in f])
    T2, a2 = 298.15, np.median([p[2] for p in t])
    s_ = -(a2 - a1) / (T2 - T1); h_ = a1 + T1 * s_
    gh.append((lo, hi, T1, a1, a2, h_, s_))
    out.append(["gibbs_helmholtz_sucrose", f"{lo}-{hi} mol/kg", "", n, 12, f"{a1:.0f}@{T1:.1f}K", f"{a2:.0f}@298.15K",
                f"h={h_:.0f}", f"s={s_:.2f}", len(f) + len(t), "CRC 2005 + Robinson & Stokes 1949"])
    print(f"  m {lo}-{hi}: a={a1:.0f} at {T1:.1f} K, {a2:.0f} at 298.15 K  ->  h = {h_:.0f} J/mol, s = {s_:.1f} J/mol/K per OH, "
          f"T s at 270 K = {270*s_:.0f}")
with open(ROOT / "predictions" / "posthoc_per_molecule.csv", "w", newline="") as fh:
    fh.write("# POST HOC: per-molecule constants and a Gibbs-Helmholtz split; NOT part of the blind test\n")
    w = csv.writer(fh)
    w.writerow(["analysis", "solute_or_window", "class", "n_OH", "n_C", "a_all_or_a1", "a_le20wt_or_a2", "rmse_K_or_h", "max_depression_K_or_s", "n_points", "data"])
    w.writerows(out)

# --- figure
plt.rcParams.update({"font.size": 7.5, "axes.linewidth": 0.6})
fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.6), gridspec_kw={"width_ratios": [1.25, 1]})
mk = {"acyclic polyol": "o", "sugar": "s", "mono-alcohol": "^"}
done = set()
for s, (a_full, a, n, nc, cls, rmse, depth, npts) in fits.items():   # plot the <= 20 wt% constant (one-constant
    weak = depth < 3.0                                                   # Margules is adequate there for every molecule)          # small depressions -> constant weakly determined
    lab = cls if cls not in done else None; done.add(cls)
    ax[0].plot(n, a, mk[cls], color=S.CLASS[cls], ms=5, mfc="white" if weak else S.CLASS[cls], mew=0.9, ls="none", label=lab)
    ax[0].annotate(s.replace("ethylene_glycol", "EG").replace("-propanol", "-PrOH").replace("methanol", "MeOH").replace("ethanol", "EtOH")
                   + (" (conc.)" if s == "sorbitol" else ""), (n, a), xytext=(4, -2), textcoords="offset points", fontsize=5.6)
ax[0].axhline(A_OH, color=S.OURS, lw=0.8, ls="--")
ax[0].text(3.6, A_OH - 420, "glycerol constant (blind tests)", fontsize=5.8, color="0.35")
ax[0].set(xlabel="Number of OH groups", ylabel="Per-molecule a (J mol$^{-1}$ per OH)", title="(a) One constant per molecule (post hoc)")
ax[0].legend(fontsize=6, frameon=False, loc="lower right")
ax[0].set_ylim(-7200, 0)
ax[0].text(0.02, 0.97, "fits to <= 20 wt% (sorbitol: 53-64 wt%)\nopen symbol: depression < 3 K, weakly determined", transform=ax[0].transAxes, fontsize=5.6, color="0.35", va="top")
frz_w = [p for p in frz if 0.3 <= p[0] <= 2.2]
ax[1].plot([p[1] for p in frz_w], [p[2] for p in frz_w], "^", color=S.OURS, ms=3.5, label="sucrose, freezing points (0.3–2.2 mol/kg)")
r25_w = [p for p in r25 if p[0] <= 2.2]
ax[1].plot([p[1] for p in r25_w], [p[2] for p in r25_w], "s", color=S.ALT, ms=3.5, mfc="white", label="sucrose, 25 °C isopiestic (0.3–2.2 mol/kg)")
Tline = np.linspace(268, 300, 10)
for (lo, hi, T1, a1, a2, h_, s_), ls, cc in zip(gh, (":", "-"), (S.OV, S.KM)):
    ax[1].plot(Tline, h_ - Tline * s_, ls, color=cc, lw=1.1, label=f"a = h − Ts, {lo}–{hi} mol/kg")
ax[1].set(xlabel="Temperature (K)", ylabel="Implied a (J mol$^{-1}$ per OH)", title="(b) Sucrose: a = h − Ts (post hoc)")
ax[1].legend(fontsize=5.6, frameon=False, loc="upper left")
for a_ in ax:
    a_.spines[["top", "right"]].set_visible(False)
fig.tight_layout(pad=0.3, w_pad=1.0)
fig.savefig(ROOT / "figures" / "posthoc_per_molecule.png", dpi=300)
fig.savefig(ROOT / "figures" / "posthoc_per_molecule.pdf")

# single-panel version of (b) for a one-column figure in the report
plt.rcParams.update({"font.size": 6.3})
fig3, cx = plt.subplots(figsize=(2.75, 1.72))
cx.plot([p[1] for p in frz_w], [p[2] for p in frz_w], "^", color=S.OURS, ms=3.4, label="sucrose, freezing points")
cx.plot([p[1] for p in r25_w], [p[2] for p in r25_w], "s", color=S.ALT, ms=3.4, mfc="white", label="sucrose, 25 °C isopiestic")
for (lo, hi, T1, a1, a2, h_, s_), ls, cc in zip(gh, (":", "-"), (S.OV, S.KM)):
    cx.plot(Tline, h_ - Tline * s_, ls, color=cc, lw=1.2, label=f"a = h − Ts ({lo}–{hi} mol/kg)")
cx.set(xlabel="Temperature (K)", ylabel="Implied a (J mol$^{-1}$ per OH)")
cx.legend(fontsize=4.8, frameon=False, loc="upper left", handlelength=1.6, labelspacing=0.25)
cx.spines[["top", "right"]].set_visible(False)
fig3.tight_layout(pad=0.3)
fig3.savefig(ROOT / "figures" / "posthoc_gibbs_helmholtz.png", dpi=300)
fig3.savefig(ROOT / "figures" / "posthoc_gibbs_helmholtz.pdf")
