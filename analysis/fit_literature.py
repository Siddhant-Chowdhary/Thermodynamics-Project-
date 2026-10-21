"""Fit a_OH on literature glycerol data (CRC, 2-40 wt%) and write blind predictions.

    python analysis/fit_literature.py  -> predictions/predictions_v1.csv, figures/fit_glycerol.png

Predictions for PG, sorbitol and mixtures are made with NO data on those solutes.
Commit the CSV before any measurement: the commit timestamp is the blind-test proof.
"""
import csv, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from thermo import IdealModel, GroupContributionModel, freezing_point, from_wt_percent
from thermo.fit import fit_group_contribution, bootstrap, gibbs_duhem_residual

DCP = True
rows = list(csv.DictReader(open(ROOT / "data" / "literature.csv")))
data = [({"glycerol": float(r["wt_percent"]) / 100}, float(r["T_freeze_C"]) + 273.15)
        for r in rows if r["solute"] == "glycerol" and r["source_id"] in ("CRC", "CRC2005")
        and float(r["wt_percent"]) >= 2.0]

fit = fit_group_contribution(data, include_dcp=DCP)
boot = bootstrap(data, fit_group_contribution, n_boot=400, include_dcp=DCP)
a = fit["param"]
print(f"n = {len(data)} CRC glycerol points (2-40 wt%)")
print(f"a_OH = {a:.1f} J/mol per OH   (A_glycerol = {3*a:.0f} J/mol)")
print(f"bootstrap 95% CI a_OH: {boot['ci95'][0]:.1f} to {boot['ci95'][1]:.1f}")
print(f"RMSE = {fit['rmse_K']*1000:.0f} mK, max |res| = {np.max(np.abs(fit['residuals_K']))*1000:.0f} mK")
print(f"Gibbs-Duhem residual = {gibbs_duhem_residual(fit['model'], 265, 'glycerol'):.1e}")
print("residuals by wt% (K):", [(round(c['glycerol']*100), round(r, 3)) for (c, _), r in zip(data, fit['residuals_K'])])

cases = [("propylene_glycol", {"propylene_glycol": w}) for w in (10, 20, 30)] + \
        [("sorbitol", {"sorbitol": w}) for w in (10, 20, 30)] + \
        [("glycerol", {"glycerol": w}) for w in (10, 20, 30)] + \
        [("gly:sorb 1:1", {"glycerol": 10, "sorbitol": 10}),
         ("gly:sorb 3:1", {"glycerol": 15, "sorbitol": 5})]
ideal = IdealModel()
samples = boot["samples"]
out = []
for label, wt in cases:
    comp = from_wt_percent(**wt)
    Tp = freezing_point(comp, fit["model"], DCP) - 273.15
    Ts = np.array([freezing_point(comp, GroupContributionModel(s), DCP) - 273.15 for s in samples])
    Ti = freezing_point(comp, ideal, DCP) - 273.15
    lo, hi = np.percentile(Ts, [2.5, 97.5])
    out.append([label, ";".join(f"{k}={v}" for k, v in wt.items()), f"{Tp:.2f}", f"{lo:.2f}", f"{hi:.2f}", f"{Ti:.2f}"])
    print(f"{label:16} {str(wt):40} pred {Tp:7.2f} C  [{lo:.2f}, {hi:.2f}]  ideal {Ti:7.2f}")

with open(ROOT / "predictions" / "predictions_v1.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow([f"# a_OH = {a:.1f} J/mol per OH, 95% CI {boot['ci95'][0]:.1f}..{boot['ci95'][1]:.1f}; fitted on CRC glycerol 2-40 wt%; dCp term on"])
    w.writerow(["case", "composition_wt_percent", "T_pred_C", "ci95_low_C", "ci95_high_C", "T_ideal_C"])
    w.writerows(out)

# figure: glycerol fit
ws = np.linspace(0.5, 45, 90)
fig, ax = plt.subplots(figsize=(5.4, 3.5))
ax.plot([c["glycerol"]*100 for c, _ in data], [T-273.15 for _, T in data], "o", mfc="white", mec="black", ms=4, label="CRC (fitted)")
ax.plot(ws, [freezing_point(from_wt_percent(glycerol=x), fit["model"], DCP)-273.15 for x in ws], "k-", lw=1.2, label=f"per-OH model, a_OH = {a:.0f} J/mol")
ax.plot(ws, [freezing_point(from_wt_percent(glycerol=x), ideal, DCP)-273.15 for x in ws], "k:", lw=1, label="ideal solution")
ax.set(xlabel="Glycerol / wt %", ylabel="Freezing point / °C")
ax.spines[["top", "right"]].set_visible(False); ax.legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(ROOT / "figures" / "fit_glycerol.png", dpi=170)
