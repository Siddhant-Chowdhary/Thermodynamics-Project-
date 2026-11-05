"""Blind tests with the committed a_OH (predictions_v1, 16 Oct; model frozen, no refit):
 (a) sorbitol vs Baer & Baldwin 1984 (cryoscope, 1-11 wt%) and Siniti et al. 1999 (DSC, 35-64 wt%)
 (b) propylene glycol vs the osmotic-virial correlation of Zielinski et al. 2014 (Table 3: PG B = 0.039 molal^-1),
     which was FITTED to PG data (their refs 5, 47, 72, 75). This is a fitted correlation, not raw data.
 (c) cross-check of Baer's cryoscope: Baer sucrose vs CRC 2005 sucrose.
At compositions not on the committed grid the frozen model is evaluated (flagged).
    python analysis/blind_test_sorbitol_pg.py -> predictions/blind_test_sorbitol_pg.csv
"""
import csv, math, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_model_comparison import T_model, A_OH
from unifac_compare import T_freeze
ROOT = Path(__file__).resolve().parents[1]
lit = list(csv.DictReader(open(ROOT / "data" / "literature.csv")))

def zielinski_dT(m, B, k=1.0):
    """Zielinski et al. 2014 eqs (2), (9): osmolality -> freezing point depression."""
    R, M1, ds, Tm = 8.314, 0.01802, 22.00, 273.15
    pi = k * m + B * (k * m) ** 2
    return R * Tm * pi * M1 / ds / (1 + R * pi * M1 / ds)

rows = []
def stats(label, pts):
    out = []
    for i, name in enumerate(["ours", "ideal", "UNIFAC", "Dortmund"]):
        e = [p[i] - T for _, T, p in pts]
        out.append((name, math.sqrt(sum(x * x for x in e) / len(e)), sum(e) / len(e)))
    print(f"{label} (n={len(pts)}): " + " | ".join(f"{n} RMSE {r:.2f} K ({s:+.2f})" for n, r, s in out))
    for n, r, s in out: rows.append([label, len(pts), n, f"{r:.3f}", f"{s:+.3f}"])

def models(sol, w):
    return [T_model({sol: w}, A_OH), T_model({sol: w}, 0.0), T_freeze({sol: w}, "UNIFAC"), T_freeze({sol: w}, "UNIFAC-Dortmund")]

for sid, lab in [("BAER1984", "sorbitol, Baer 1984 cryoscope, 1-11 wt%"), ("SINITI1999", "sorbitol, Siniti 1999 DSC, 35-64 wt%")]:
    pts = [(float(r["wt_percent"]), float(r["T_freeze_C"]), models("sorbitol", float(r["wt_percent"]))) for r in lit
           if r["solute"] == "sorbitol" and r["source_id"] == sid]
    stats(lab, pts)
    if sid == "SINITI1999":
        stats("sorbitol, Siniti 1999 without suspect x=0.05 point", [p for p in pts if p[0] > 40])
    for w, T, p in pts: print(f"   {w:6.2f} wt%  data {T:7.2f}  ours {p[0]:7.2f}  ideal {p[1]:7.2f}  UNIFAC {p[2]:7.2f}  Dortmund {p[3]:7.2f}")

print("\nPG: committed predictions vs Zielinski 2014 PG correlation (fitted to PG data)")
for w in (10, 20, 30, 40):
    m = w / (100 - w) * 1000 / 76.094
    Tz = -zielinski_dT(m, 0.039)
    p = models("propylene_glycol", w)
    print(f"   PG {w} wt%: Zielinski {Tz:7.2f} | ours {p[0]:7.2f} ideal {p[1]:7.2f} UNIFAC {p[2]:7.2f} Dortmund {p[3]:7.2f}")
    rows.append([f"PG {w} wt% vs Zielinski correlation", 1, "ours-minus-Zielinski", f"{p[0]-Tz:+.3f}", ""])
print("\nSanity: Zielinski glycerol correlation vs CRC 40 wt%:", round(-zielinski_dT(7.239, 0.023), 2), "(CRC -15.50)")

crc = {float(r["wt_percent"]): float(r["T_freeze_C"]) for r in lit if r["solute"] == "sucrose" and r["source_id"] == "CRC2005"}
baer = {float(r["wt_percent"]): float(r["T_freeze_C"]) for r in lit if r["solute"] == "sucrose" and r["source_id"] == "BAER1984"}
common = [w for w in baer if w in crc]
print("\nBaer cryoscope check, sucrose vs CRC:", [(w, baer[w], crc[w]) for w in sorted(common)])
with open(ROOT / "predictions" / "blind_test_sorbitol_pg.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["test", "n", "model", "RMSE_or_diff_K", "mean_signed_error_K"]); w.writerows(rows)
