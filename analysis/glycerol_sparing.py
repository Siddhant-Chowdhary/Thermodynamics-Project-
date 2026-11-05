"""Recipe as a glycerol-sparing curve (decided 21 Oct with Sid).

For each sorbitol dose added to 0.9 wt% saline, find the MINIMUM glycerol dose that keeps the NaCl
mass fraction in the residual (unfrozen) liquid at or below the salt target at the design temperature.
No toxicity weights (cannot be sourced). Zero-glycerol recipes are reported but marked biologically
unacceptable: sorbitol does not penetrate cells (Pegg 2007).

Inputs (data/design_inputs.csv): salt target 4.41 wt% NaCl (Lovelock 1953 via Shepard 1976);
design T -15 C (start of the lethal zone, Gao & Critser 2000) with -10 and -20 C as a sweep.
Robustness: polyol activity from (i) per-OH model (this work), (ii) ideal, (iii) UNIFAC; NaCl always from
the fitted phi(m) (CRC 2005), combined additively in ln a_w (pre-registered 21 Oct).

Self-contained (does not import the local 'thermo' package, because the pip 'thermo' is used for UNIFAC).
    python analysis/glycerol_sparing.py -> predictions/recipe_v2_glycerol_sparing.csv, figures/glycerol_sparing.png
"""
import csv, math, sys
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0, str(Path(__file__).resolve().parent))
from unifac_compare import ln_aw_ice, R, ORIG, UFIP, UFSG
from thermo.unifac import UNIFAC
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
MW = {"water": 18.015, "glycerol": 92.094, "sorbitol": 182.172, "NaCl": 58.443}
NOH = {"glycerol": 3, "sorbitol": 6}
A_OH = float(open(ROOT / "predictions" / "predictions_v1.csv").readline().split("a_OH = ")[1].split()[0])
inp = {r["key"]: r for r in csv.DictReader(open(ROOT / "data" / "design_inputs.csv"))}
SALT_TARGET = float(inp["salt_target"]["value"]); BASE_NACL = float(inp["base_NaCl"]["value"])

# NaCl phi(m), refitted here exactly as in thermo/salt.py (CRC 2005 rows, linear least squares)
rows = [r for r in csv.DictReader(open(ROOT / "data" / "literature.csv")) if r["solute"] == "NaCl" and r["source_id"] == "CRC2005"]
X, y = [], []
for r in rows:
    w = float(r["wt_percent"]); m = w / (100 - w) * 1000 / MW["NaCl"]
    phi = -ln_aw_ice(float(r["T_freeze_C"]) + 273.15) / (2 * m * MW["water"] / 1000); s = math.sqrt(m)
    X.append([-s / (1 + s), m, m * m]); y.append(phi - 1)
PHI = np.linalg.lstsq(np.array(X), np.array(y), rcond=None)[0]

def ln_aw(mass, T, model):
    """mass: dict of solute masses (g) incl. 'water'. Additive polyols + NaCl."""
    nw = mass["water"] / MW["water"]
    polys = {k: mass[k] / MW[k] for k in ("glycerol", "sorbitol") if mass.get(k, 0) > 0}
    val = 0.0
    if polys:
        tot = nw + sum(polys.values()); xw = nw / tot; xs = {k: v / tot for k, v in polys.items()}
        if model == "per-OH":
            lg = (1 - xw) * sum(NOH[k] * A_OH / (R * T) * x for k, x in xs.items())
        elif model == "ideal":
            lg = 0.0
        else:
            names = list(xs); cg = [ORIG["water"]] + [ORIG[k] for k in names]
            lg = math.log(UNIFAC.from_subgroups(T=T, xs=[xw] + [xs[k] for k in names], chemgroups=cg, version=0,
                                                interaction_data=UFIP, subgroups=UFSG).gammas()[0])
        val += math.log(xw) + lg
    if mass.get("NaCl", 0) > 0:
        m = mass["NaCl"] / MW["NaCl"] / (mass["water"] / 1000); s = math.sqrt(m)
        val += -2 * m * MW["water"] / 1000 * (1 - PHI[0] * s / (1 + s) + PHI[1] * m + PHI[2] * m * m)
    return val

def residual_salt(wg, ws, T, model):
    """NaCl mass fraction in the residual liquid at T (per 100 g of starting solution)."""
    base = {"glycerol": 100 * wg, "sorbitol": 100 * ws, "NaCl": 100 * BASE_NACL}
    w0 = 100 - sum(base.values())
    target = ln_aw_ice(T)
    f = lambda frac: ln_aw(dict(base, water=w0 * (1 - frac)), T, model) - target
    if f(0.0) <= 0:      # no ice at this T
        frac = 0.0
    else:
        frac = brentq(f, 0.0, 0.999, xtol=1e-10)
    liquid = w0 * (1 - frac) + sum(base.values())
    return base["NaCl"] / liquid

def min_glycerol(ws, T, model):
    g = lambda wg: residual_salt(wg, ws, T, model) - SALT_TARGET
    if g(0.0) <= 0: return 0.0
    return brentq(g, 0.0, 0.45, xtol=1e-5)

if __name__ == "__main__":
    MODELS = ["per-OH", "ideal", "UNIFAC"]; TS = [-10.0, -15.0, -20.0]; WS = [i / 100 for i in range(0, 13)]
    out = []
    for T in TS:
        for mdl in MODELS:
            for ws in WS:
                wg = min_glycerol(ws, T + 273.15, mdl); out.append([T, mdl, round(100 * ws, 1), round(100 * wg, 2)])
    with open(ROOT / "predictions" / "recipe_v2_glycerol_sparing.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([f"# min glycerol (wt%) for residual NaCl <= {100*SALT_TARGET:.2f} wt% in {100*BASE_NACL:.1f} wt% saline; a_OH={A_OH}; committed before lab session 2"])
        w.writerow(["T_design_C", "polyol_model", "sorbitol_wt", "min_glycerol_wt"]); w.writerows(out)
    def get(T, mdl): return [r[3] for r in out if r[0] == T and r[1] == mdl]
    for T in TS:
        print(f"T = {T:5.1f} C | " + " | ".join(f"{m}: " + ", ".join(f"{s:.0f}:{g:.1f}" for s, g in zip([100*x for x in WS][:7], get(T, m)[:7])) for m in MODELS))
    # replacement ratio near zero sorbitol at -15 C (glycerol saved per g sorbitol)
    for mdl in MODELS:
        g = get(-15.0, mdl); print(f"  replacement ratio at -15 C ({mdl}): {(g[0]-g[2])/2:.2f} g glycerol saved per g sorbitol (0->2 wt%)")
    # lab recipe point: sorbitol 5 wt% on the -15 C per-OH curve
    wg5 = min_glycerol(0.05, 258.15, "per-OH")
    print(f"Proposed lab recipe point: glycerol {100*wg5:.2f} wt% + sorbitol 5 wt% + NaCl 0.9 wt%")
    # its freezing point by each model (no ice yet)
    for mdl in MODELS:
        base = {"glycerol": 100 * wg5, "sorbitol": 5.0, "NaCl": 0.9}; base["water"] = 100 - sum(base.values())
        Tf = brentq(lambda T: ln_aw(base, T, mdl) - ln_aw_ice(T), 200, 273.15 - 1e-9) - 273.15
        print(f"   predicted freezing point ({mdl}): {Tf:.2f} C")
    # figure
    plt.rcParams.update({"font.size": 8})
    fig, axs = plt.subplots(1, 2, figsize=(8, 3.3), sharey=True)
    xs = [100 * x for x in WS]
    for mdl, ls in zip(MODELS, ["k-", "k:", "k--"]):
        axs[0].plot(xs, get(-15.0, mdl), ls, lw=1.3 if mdl == "per-OH" else 1, label=mdl if mdl != "per-OH" else "per-OH model (this work)")
    for T, ls in zip(TS, ["k:", "k-", "k--"]):
        axs[1].plot(xs, get(T, "per-OH"), ls, lw=1.2, label=f"{T:.0f} °C")
    for ax, t in zip(axs, ["Design T = −15 °C: three polyol models", "Per-OH model: design temperature sweep"]):
        ax.axhspan(-0.5, 0.3, color="0.85", lw=0)
        ax.set(title=t, xlabel="Sorbitol added / wt %"); ax.spines[["top", "right"]].set_visible(False); ax.legend(frameon=False, fontsize=7)
    axs[0].set_ylabel("Minimum glycerol / wt %")
    axs[0].text(6.5, 0.6, "grey band: no penetrating agent\n(biologically unacceptable)", fontsize=6.5, color="0.3")
    fig.suptitle(f"Glycerol needed to keep residual NaCl ≤ {100*SALT_TARGET:.2f} wt% (0.9 wt% saline base)", fontsize=8.5)
    fig.tight_layout(); fig.savefig(ROOT / "figures" / "glycerol_sparing.png", dpi=170)
