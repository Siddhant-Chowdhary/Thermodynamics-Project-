"""Sanity checks on data/literature.csv.

For every ice-branch point: ideal-solution freezing point, deviation, and the
apparent osmotic coefficient phi = -ln a_w / (nu * m * M_w), where a_w is fixed
by the measured T via ice equilibrium. phi ~ 1 means near-ideal; a smooth phi
vs m within one source is a quick transcription check.
"""
import csv, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from thermo.sle import freezing_point, ln_aw_ice
from thermo.activity import IdealModel
from thermo.constants import solute, M_WATER

path = Path(__file__).resolve().parents[1] / "data" / "literature.csv"
rows = list(csv.DictReader(open(path)))
model = IdealModel()
print(f"{'solute':12}{'wt%':>6}{'src':>10}{'T_lit':>8}{'T_ideal':>9}{'dev':>7}{'m':>7}{'phi':>7}")
for r in rows:
    if r["phase_branch"] != "ice":
        continue
    s, w, T = r["solute"], float(r["wt_percent"]), float(r["T_freeze_C"])
    Ti = freezing_point({s: w / 100}, model, include_dcp=True) - 273.15
    m = w / (100 - w) * 1000 / solute(s)["molar_mass"]
    phi = -ln_aw_ice(T + 273.15, include_dcp=True) / (solute(s)["nu"] * m * M_WATER / 1000)
    print(f"{s[:12]:12}{w:6.1f}{r['source_id']:>10}{T:8.2f}{Ti:9.2f}{T - Ti:7.2f}{m:7.3f}{phi:7.3f}")
