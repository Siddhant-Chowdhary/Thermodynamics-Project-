"""Blind predictions (all four models) for the short lab set, written BEFORE measuring.
    python analysis/predict_lab_set.py -> predictions/lab_set_v2.csv
"""
import csv, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from unifac_compare import T_freeze
from plot_model_comparison import T_model, A_OH
ROOT = Path(__file__).resolve().parents[1]
SET = [("glycerol 30 (method check vs CRC)", {"glycerol": 30}),
       ("PG 30", {"propylene_glycol": 30}),
       ("sorbitol 40", {"sorbitol": 40}),
       ("glycerol 15 + sorbitol 15", {"glycerol": 15, "sorbitol": 15}),
       ("PG 20 (if time)", {"propylene_glycol": 20})]
rows = []
for label, wt in SET:
    r = [T_model(wt, A_OH), T_model(wt, 0.0), T_freeze(wt, "UNIFAC"), T_freeze(wt, "UNIFAC-Dortmund")]
    rows.append([label] + [f"{v:.2f}" for v in r] + [f"{max(r)-min(r):.1f}"])
    print(f"{label:34} ours {r[0]:7.2f} ideal {r[1]:7.2f} UNIFAC {r[2]:7.2f} Dortmund {r[3]:7.2f}  spread {max(r)-min(r):.1f}")
with open(ROOT / "predictions" / "lab_set_v2.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["sample", "T_ours_C", "T_ideal_C", "T_UNIFAC_C", "T_Dortmund_C", "spread_K"]); w.writerows(rows)
