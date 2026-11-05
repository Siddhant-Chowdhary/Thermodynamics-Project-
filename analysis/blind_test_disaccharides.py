"""Score blind test #3: CRC 2005 maltose and lactose vs predictions/disaccharides_v1.csv (committed before the data rows
were entered). No refit.  python analysis/blind_test_disaccharides.py -> predictions/disaccharides_blind_results.csv"""
import csv, math
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
P = {(r["solute"], float(r["wt_percent"])): r for r in csv.DictReader(
    l for l in open(ROOT / "predictions" / "disaccharides_v1.csv") if not l.startswith("#"))}
lit = [r for r in csv.DictReader(open(ROOT / "data" / "literature.csv")) if r["solute"] in ("maltose", "lactose")]
out = []
for s in ("maltose", "lactose"):
    e = {"per-OH": [], "ideal": []}; wmax = 0
    for r in lit:
        w = float(r["wt_percent"])
        if r["solute"] == s and (s, w) in P:
            T = float(r["T_freeze_C"]); wmax = max(wmax, w)
            e["per-OH"].append(float(P[(s, w)]["T_perOH_C"]) - T); e["ideal"].append(float(P[(s, w)]["T_ideal_C"]) - T)
    for m, v in e.items():
        out.append([s, m, len(v), f"{wmax:g}", f"{math.sqrt(sum(x*x for x in v)/len(v)):.3f}", f"{max(map(abs, v)):.3f}"])
with open(ROOT / "predictions" / "disaccharides_blind_results.csv", "w", newline="") as fh:
    fh.write("# blind test #3, CRC 2005 pp. 8-67/8-68; scored only at grid points present in both files\n")
    w = csv.writer(fh); w.writerow(["solute", "model", "n", "max_wt", "RMSE_K", "max_abs_K"]); w.writerows(out)
for r in out: print(r)
