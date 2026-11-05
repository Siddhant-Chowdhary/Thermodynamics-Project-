"""Temperature-transfer and mixing-rule test with Bower & Robinson, J. Phys. Chem. 67, 1540 (1963), 25 C isopiestic data.

EVALUATION ONLY, nothing refitted. a_OH is the committed value (predictions_v1.csv header). The temperature
assumption was fixed before this data was seen (PROMPTS.md, committed 17 Oct): the interaction energy A = n_OH*a_OH
is temperature-independent, so alpha = A/RT. No prediction file was committed first because, with a_OH and that
assumption fixed, the model has no remaining free choice.

Part A (Table I, binary sorbitol-water at 25 C): osmotic coefficient phi = -ln a_w / (m M_w), data vs
  per-OH (A constant, evaluated at 298.15 K)  <- the model as specified
  per-OH (alpha frozen at 273.15 K)           <- alternative: purely entropic interaction, shown for contrast
  ideal (Raoult, mole fraction), UNIFAC, UNIFAC-Dortmund (pip thermo, 298.15 K)
  plus phi derived from Baer 1984 sorbitol freezing points (~0 C) against the per-OH model on its liquidus,
  to separate a temperature effect from a molecule effect.
Part B (Table II, sorbitol + NaCl + water): the authors' Delta = deviation of the mixture's osmolality from the
  sum of the binary contributions at the mixture's molalities. Our pre-registered mixing rule (additive ln a_w at
  equal molality) assumes Delta = 0. Convert Delta to an equivalent freezing-point error 1.86 K kg/mol * |Delta|
  (assumes the 25 C cross term carries to ~0 C; stated as an assumption) and evaluate it, with the authors'
  eq. 1, at the recipe's residual-liquid compositions at -15 C.

Run from the repo root (pip 'thermo' clash):  python analysis/sorbitol_25C.py
  -> predictions/sorbitol_25C_results.csv, figures/sorbitol_25C.png
"""
import csv, math, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from unifac_compare import VERSIONS, ln_aw_ice, R
from thermo.unifac import UNIFAC
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys as _sys
_sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from analysis import style as S

ROOT = Path(__file__).resolve().parents[1]
MW_W, MW_S = 0.018015, 0.182172          # kg/mol
NOH = 6
A_OH = float(open(ROOT / "predictions" / "predictions_v1.csv").readline().split("a_OH = ")[1].split()[0])
T25 = 298.15
KF = 8.314 * 273.15 ** 2 / 6009.5 * MW_W   # K kg/mol, cryoscopic constant of water (1.86)

rows = list(csv.DictReader(l for l in open(ROOT / "data" / "bower1963_sorbitol_25C.csv") if not l.startswith("#")))
t1 = [r for r in rows if r["table"] == "I"]
t2 = [r for r in rows if r["table"] == "II"]


def xs_of(m):
    return m / (1 / MW_W + m)


def phi_perOH(m, T, T_alpha=None):
    x = xs_of(m); alpha = NOH * A_OH / (R * (T_alpha or T))
    ln_aw = math.log(1 - x) + alpha * x * x          # binary: ln g_w = (1 - x_w) alpha x_s = alpha x_s^2
    return -ln_aw / (m * MW_W)


def phi_ideal(m):
    return -math.log(1 - xs_of(m)) / (m * MW_W)


def phi_unifac(m, T, version):
    groups, ver, ip, sg = VERSIONS[version]
    x = xs_of(m)
    g = UNIFAC.from_subgroups(T=T, xs=[1 - x, x], chemgroups=[groups["water"], groups["sorbitol"]],
                              version=ver, interaction_data=ip, subgroups=sg).gammas()[0]
    return -(math.log(1 - x) + math.log(g)) / (m * MW_W)


def a_per_oh_implied(m, phi, T):
    x = xs_of(m); ln_aw = -phi * m * MW_W
    return (ln_aw - math.log(1 - x)) / (x * x) * R * T / NOH


out = []
lines = [f"a_OH committed = {A_OH} J/mol per OH. Bower & Robinson 1963, 25 C."]
lines.append("Part A, Table I: m (mol/kg) | phi data | per-OH(A const) | per-OH(alpha@273K) | ideal | UNIFAC | Dortmund | implied a_OH (J/mol)")
errs = {k: [] for k in ("perOH", "perOH_alpha273", "ideal", "UNIFAC", "Dortmund")}
for r in t1:
    m, phi = float(r["m_B"]), float(r["phi_B"])
    p = dict(perOH=phi_perOH(m, T25), perOH_alpha273=phi_perOH(m, T25, 273.15), ideal=phi_ideal(m),
             UNIFAC=phi_unifac(m, T25, "UNIFAC"), Dortmund=phi_unifac(m, T25, "UNIFAC-Dortmund"))
    for k in errs:
        errs[k].append(p[k] - phi)
    aimp = a_per_oh_implied(m, phi, T25)
    lines.append(f"  {m:6.4f} | {phi:.4f} | {p['perOH']:.4f} | {p['perOH_alpha273']:.4f} | {p['ideal']:.4f} | "
                 f"{p['UNIFAC']:.4f} | {p['Dortmund']:.4f} | {aimp:7.0f}")
    out.append(dict(part="A_25C", m_sorbitol=m, m_NaCl="", phi_data=phi, **{f"phi_{k}": v for k, v in p.items()},
                    a_OH_implied=aimp, delta_osm="", dT_equiv_K=""))
lines.append("  RMS error in (phi - 1), i.e. in the non-ideal part: " +
             ", ".join(f"{k} {math.sqrt(np.mean(np.square(v))):.3f}" for k, v in errs.items()))
mean_dev = np.mean([abs(float(r['phi_B']) - 1) for r in t1])
lines.append(f"  (mean |phi - 1| in the data = {mean_dev:.3f}, for scale)")

# Baer 1984 at the freezing line: phi from ice equilibrium, and the per-OH model there
lit = [r for r in csv.DictReader(open(ROOT / "data" / "literature.csv")) if r["source_id"] == "BAER1984" and r["solute"] == "sorbitol"]
baer = []
lines.append("Part A2, Baer 1984 (freezing line): m | T (C) | phi data | per-OH at that T | implied a_OH")
for r in lit:
    w = float(r["wt_percent"]); T = float(r["T_freeze_C"]) + 273.15
    m = w / (100 - w) / MW_S
    phi = -ln_aw_ice(T) / (m * MW_W)
    baer.append((m, phi, phi_perOH(m, T)))
    lines.append(f"  {m:6.4f} | {T-273.15:6.3f} | {phi:.4f} | {phi_perOH(m, T):.4f} | {a_per_oh_implied(m, phi, T):7.0f}")
    out.append(dict(part="A2_freezing", m_sorbitol=m, m_NaCl="", phi_data=phi, phi_perOH=phi_perOH(m, T),
                    phi_perOH_alpha273="", phi_ideal=phi_ideal(m), phi_UNIFAC="", phi_Dortmund="",
                    a_OH_implied=a_per_oh_implied(m, phi, T), delta_osm="", dT_equiv_K=""))

# Part B: cross term
lines.append("Part B, Table II: m_sorbitol | m_NaCl | Delta (mol/kg) | Delta / ideal osmolality | equiv. freezing-point error (K)")
for r in t2:
    mB, mC = float(r["m_B"]), float(r["m_C"])
    d = -float(r["neg_delta_over_mBmC"]) * mB * mC
    lines.append(f"  {mB:.4f} | {mC:.4f} | {d:+.4f} | {abs(d)/(2*mC+mB)*100:4.1f}% | {KF*abs(d):.3f}")
    out.append(dict(part="B_ternary", m_sorbitol=mB, m_NaCl=mC, phi_data="", phi_perOH="", phi_perOH_alpha273="",
                    phi_ideal="", phi_UNIFAC="", phi_Dortmund="", a_OH_implied="", delta_osm=d, dT_equiv_K=KF * abs(d)))


def eq1(mB, mC):   # authors' least-squares fit, Delta/(mB mC)
    return -0.05326 + 0.00261 * mC - 0.00235 * mC ** 2 + 0.0134 * mB


# recipe residual liquid at -15 C on the per-OH glycerol-sparing curve: NaCl exactly at target, all solutes
# concentrated by the same factor (target / base NaCl)
rec = [r for r in csv.DictReader(l for l in open(ROOT / "predictions" / "recipe_v2_glycerol_sparing.csv") if not l.startswith("#"))
       if r["polyol_model"] == "per-OH" and float(r["T_design_C"]) == -15.0]
inp = {r["key"]: r for r in csv.DictReader(open(ROOT / "data" / "design_inputs.csv"))}
target, base = float(inp["salt_target"]["value"]), float(inp["base_NaCl"]["value"])
f = target / base
lines.append("Part B2, recipe residual liquid at -15 C (per-OH curve), sorbitol-NaCl cross term via eq. 1 "
             "(sorbitol molality beyond the data range 0.22-1.13 mol/kg is extrapolation):")
lines.append("  sorbitol wt% | glycerol wt% | m_sorb | m_NaCl (liquid) | Delta | equiv. error (K)")
for r in rec:
    ws, wg = float(r["sorbitol_wt"]) / 100, float(r["min_glycerol_wt"]) / 100
    if ws == 0:
        continue
    gl = {"s": ws * f, "g": wg * f, "n": base * f}
    ww = 1 - sum(gl.values())
    mB, mC = gl["s"] / MW_S / ww, gl["n"] / 0.058443 / ww
    d = eq1(mB, mC) * mB * mC
    lines.append(f"  {100*ws:5.1f} | {100*wg:5.2f} | {mB:5.2f} | {mC:5.2f} | {d:+.3f} | {KF*abs(d):.2f}")

txt = "\n".join(lines)
print(txt)
with open(ROOT / "predictions" / "sorbitol_25C_results.csv", "w", newline="") as fh:
    fh.write("# " + txt.replace("\n", "\n# ") + "\n")
    w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
    w.writeheader()
    for r in out:
        w.writerow({k: (f"{v:.5g}" if isinstance(v, float) else v) for k, v in r.items()})

# figure: phi vs m, b/w
plt.rcParams.update({"font.size": 6.5, "axes.linewidth": 0.6})
fig, ax = plt.subplots(1, 2, figsize=(4.3, 1.97))
mm = np.linspace(0.05, 3.9, 120)
ax[0].plot(mm, [phi_perOH(m, T25) for m in mm], "-", color=S.OURS, lw=1.4, label="per-OH, A const (model)")
ax[0].plot(mm, [phi_perOH(m, T25, 273.15) for m in mm], "--", color=S.OURS, lw=1.1, label="per-OH, α fixed at 273 K")
ax[0].plot(mm, [phi_unifac(m, T25, "UNIFAC") for m in mm], "-.", color=S.UNIFAC, lw=1.3, label="UNIFAC")
ax[0].plot(mm, [phi_ideal(m) for m in mm], ":", color=S.IDEAL, lw=1.4, label="ideal")
ax[0].plot([float(r["m_B"]) for r in t1], [float(r["phi_B"]) for r in t1], "s", color=S.DATA, ms=3.5, mfc="white", mew=1.0, label="Bower & Robinson 1963, 25 °C")
ax[0].plot([b[0] for b in baer], [b[1] for b in baer], "^", color=S.ALT, ms=3, label="Baer 1984, at freezing (~0 °C)")
ax[0].set(xlabel="Sorbitol molality (mol/kg)", ylabel="Osmotic coefficient φ", title="Binary sorbitol–water")
ax[0].set_ylim(0.86, 1.62)
ax[0].legend(fontsize=4.9, frameon=False, loc="upper left")
mBs = [float(r["m_B"]) for r in t2]; mCs = [float(r["m_C"]) for r in t2]
dT = [KF * float(r["neg_delta_over_mBmC"]) * float(r["m_B"]) * float(r["m_C"]) for r in t2]
for (mk, lo, hi), cc in zip((("o", 0, 1), ("s", 1, 3), ("^", 3, 4.5), ("D", 4.5, 7)), (S.SKY, S.OURS, S.OV, S.ALT)):
    sel = [i for i, c in enumerate(mCs) if lo <= c < hi]
    ax[1].plot([mBs[i] for i in sel], [dT[i] for i in sel], mk, color=cc, ms=3.5, mfc=cc, mec="white", mew=0.4, ls="none",
               label=f"NaCl {lo}–{hi} mol/kg")
ax[1].set(xlabel="Sorbitol molality (mol/kg)", ylabel="Neglected cross term as ΔT (K)",
          title="Sorbitol + NaCl: cost of additivity")
ax[1].legend(fontsize=5, frameon=False)
fig.tight_layout()
fig.savefig(ROOT / "figures" / "sorbitol_25C.png", dpi=300)
fig.savefig(ROOT / "figures" / "sorbitol_25C.pdf")

# single-panel (osmotic coefficient only) for a one-column figure in the report; the cross term is reported in the text
plt.rcParams.update({"font.size": 6.3})
fig1, bx = plt.subplots(figsize=(2.75, 1.72))
bx.plot(mm, [phi_perOH(m, T25) for m in mm], "-", color=S.OURS, lw=1.4, label="per-OH, glycerol constant")
bx.plot(mm, [phi_perOH(m, T25, 273.15) for m in mm], "--", color=S.OURS, lw=1.1, label="per-OH, α fixed at 273 K")
bx.plot(mm, [phi_unifac(m, T25, "UNIFAC") for m in mm], "-.", color=S.UNIFAC, lw=1.3, label="UNIFAC")
bx.plot(mm, [phi_ideal(m) for m in mm], ":", color=S.IDEAL, lw=1.4, label="ideal")
bx.plot([float(r["m_B"]) for r in t1], [float(r["phi_B"]) for r in t1], "s", color=S.DATA, ms=3.2, mfc="white", mew=0.9,
        label="Bower & Robinson 1963, 25 °C")
bx.plot([b[0] for b in baer], [b[1] for b in baer], "^", color=S.ALT, ms=3.2, label="Baer 1984, near 0 °C (flagged)")
bx.set(xlabel="Sorbitol molality (mol/kg)", ylabel="Osmotic coefficient φ", ylim=(0.86, 1.62))
bx.legend(fontsize=4.8, frameon=False, loc="upper left", handlelength=1.6, labelspacing=0.25)
bx.spines[["top", "right"]].set_visible(False)
fig1.tight_layout(pad=0.2)
fig1.savefig(ROOT / "figures" / "sorbitol_25C_phi.pdf")
fig1.savefig(ROOT / "figures" / "sorbitol_25C_phi.png", dpi=300)
