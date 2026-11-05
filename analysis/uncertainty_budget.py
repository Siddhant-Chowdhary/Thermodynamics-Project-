"""Uncertainty budget for a measured freezing point (GUM-style, standard uncertainties added in quadrature).

Inputs: data/lab/uncertainty_inputs.csv. Rows marked ASSUMED are typical values for a PRE-LAB PREVIEW and must be
replaced with measured ones after probe calibration and each run. Outputs carry a METHOD DEMO stamp while any
input is ASSUMED.

Components (u = standard uncertainty, K):
  resolution       DS18B20 step / sqrt(12)
  calibration      sd/sqrt(n) of ice-point readings, plus how well the ice bath realises 0 C
  nonlinearity     offset change between 0 C and the sample's T (one-point calibration) / sqrt(3)
  composition      weighing: balance readability on solute and water masses -> u(w) -> |dT/dw| u(w)
  purity           water in the supplied solute (rectangular) -> |dT/dw| w p / sqrt(3)
  supercooling     ice formed at nucleation, f = cp dT_sc / dH_fus, concentrates the liquid; we CORRECT for it
                   (shift computed with the per-OH model) and carry a fraction of the correction as uncertainty
  peak picking     scatter of raw points at the post-recalescence maximum
  repeatability    standard uncertainty of the mean of two runs = |difference| / 2
Combined u_c = quadrature sum; expanded U = 2 u_c (~95 %).

    python analysis/uncertainty_budget.py -> predictions/uncertainty_budget_preview.csv
"""
import csv, math, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from thermo import GroupContributionModel, freezing_point
from thermo.salt import AdditiveSaltModel

a_oh = float(open(ROOT / "predictions" / "predictions_v1.csv").readline().split("a_OH = ")[1].split()[0])
MODEL = AdditiveSaltModel(GroupContributionModel(a_oh))
CP_W, DH_W = 4.18, 333.6        # J/g/K, J/g

inp = {r["key"]: r for r in csv.DictReader(l for l in open(ROOT / "data" / "lab" / "uncertainty_inputs.csv") if not l.startswith("#"))}
v = {k: float(r["value"]) for k, r in inp.items()}
assumed = sorted(k for k, r in inp.items() if r["status"] == "ASSUMED")

SAMPLES = {"glycerol 30": {"glycerol": 0.30}, "PG 30": {"propylene_glycol": 0.30}, "sorbitol 40": {"sorbitol": 0.40},
           "glycerol 15 + sorbitol 15": {"glycerol": 0.15, "sorbitol": 0.15},
           "glycerol 20 + NaCl 5": {"glycerol": 0.20, "NaCl": 0.05}}


def Tf(comp):
    return freezing_point(comp, MODEL, True) - 273.15


def budget(comp):
    w = sum(comp.values()); T0 = Tf(comp)
    # sensitivity to total solute mass fraction (all solutes scaled together)
    h = 1e-4
    dTdw = (Tf({k: x * (w + h) / w for k, x in comp.items()}) - Tf({k: x * (w - h) / w for k, x in comp.items()})) / (2 * h)
    M = v["sample_mass"]; ms = w * M; mw = M - ms; uw_mass = v["balance_resolution"] / math.sqrt(3)
    u_w = math.hypot(mw / M ** 2 * uw_mass, ms / M ** 2 * uw_mass)
    # supercooling: fraction of WATER frozen at nucleation; liquid concentrates; correction = shift in T
    f = CP_W * v["supercooling"] / DH_W
    w_liq = ms / (ms + mw * (1 - f))
    shift = Tf({k: x * w_liq / w for k, x in comp.items()}) - T0
    u = {
        "resolution": v["probe_resolution"] / math.sqrt(12),
        "calibration": math.hypot(v["cal_scatter_sd"] / math.sqrt(v["cal_n_readings"]), v["ice_point_realisation"]),
        "nonlinearity": v["probe_nonlinearity"] * min(1.0, abs(T0) / 10.0) / math.sqrt(3),
        "composition": abs(dTdw) * u_w,
        "purity": abs(dTdw) * w * v["purity_water_fraction"] / math.sqrt(3),
        "supercooling": abs(shift) * v["supercool_corr_unc_frac"],
        "peak picking": v["peak_pick_sd"],
        "repeatability": v["duplicate_difference"] / 2,
    }
    uc = math.sqrt(sum(x * x for x in u.values()))
    return T0, dTdw, shift, u, uc


rows = []
stamp = "METHOD DEMO (assumed inputs: " + ", ".join(assumed) + ")" if assumed else "all inputs measured/spec"
print(stamp)
for name, comp in SAMPLES.items():
    T0, dTdw, shift, u, uc = budget(comp)
    big = max(u, key=u.get)
    print(f"\n{name}: predicted Tf {T0:.2f} C, dT/dw {dTdw:.1f} K per unit mass fraction, supercooling shift {shift:+.3f} K "
          f"(add {-shift:+.3f} K to the measured peak)")
    for k, x in u.items():
        print(f"   {k:14s} {x:.3f} K")
    print(f"   combined u_c = {uc:.3f} K, expanded U (k=2) = {2*uc:.2f} K; largest term: {big}")
    rows.append([name, f"{T0:.2f}", f"{dTdw:.1f}", f"{-shift:+.3f}"] + [f"{u[k]:.3f}" for k in u] + [f"{uc:.3f}", f"{2*uc:.2f}", big])

with open(ROOT / "predictions" / "uncertainty_budget_preview.csv", "w", newline="") as fh:
    fh.write(f"# {stamp}\n")
    w = csv.writer(fh)
    w.writerow(["sample", "T_pred_C", "dT_dw_K", "supercool_correction_K"] + [f"u_{k.replace(' ', '_')}_K" for k in u] +
               ["u_combined_K", "U_expanded_k2_K", "largest_term"])
    w.writerows(rows)
