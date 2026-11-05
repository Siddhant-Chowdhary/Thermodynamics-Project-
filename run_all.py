"""Rebuild every figure and results table in one command, and check that the committed predictions reproduce.

    python run_all.py            # ~ a few minutes

How it protects the blind tests: everything runs in a TEMPORARY COPY of the repo, so no committed file is ever
overwritten while running. Afterwards:
  * figures/ and the evaluation result tables are copied back into the repo;
  * committed PREDICTION files are never copied back; the regenerated copy is only compared byte-for-byte with the
    committed one and reported as "reproduces" or "DIFFERS".
"""
import filecmp, shutil, subprocess, sys, tempfile, time
from pathlib import Path

ROOT = Path(__file__).resolve().parent

# committed, pre-registered prediction files: never overwritten, only checked
# Known, explained differences when a prediction file is regenerated with today's code/data:
EXPECTED_DIFF = {
    "salt_set_v1.csv": "NaCl phi was refit on CRC 2005 rows (commit 35eae3f) AFTER this file was committed; one value moves 0.01 K",
    "recipe_v1.csv": "superseded METHOD DEMO from 17 Oct, made before the NaCl refit and sourced inputs",
}
PREDICTIONS = ["disaccharides_v1.csv", "predictions_v1.csv", "unifac_v1.csv", "lab_set_v2.csv", "salt_set_v1.csv", "crc_extra_v1.csv",
               "e1_alcohols_sugars_v1.csv", "recipe_v2_glycerol_sparing.csv", "recipe_v1.csv", "lab_set_v3_glycerol_nacl.csv"]

# (script, what it does). Order matters only for the first one.
STEPS = [
    ("analysis/fit_literature.py", "fit a_OH on CRC glycerol (prediction file checked, not replaced)"),
    ("analysis/unifac_compare.py", "UNIFAC / Dortmund on glycerol"),
    ("analysis/plot_model_comparison.py", "figure: model comparison"),
    ("analysis/model_checks.py", "a_OH stability, second-parameter check"),
    ("analysis/predict_crc_extra.py", "re-create blind predictions (check only)"),
    ("analysis/predict_e1.py", "re-create E1 blind predictions (check only)"),
    ("analysis/predict_lab_set.py", "re-create lab predictions (check only)"),
    ("analysis/predict_salt_set.py", "re-create salt predictions (check only)"),
    ("analysis/blind_test_crc.py", "blind test: ethylene glycol, mannitol"),
    ("analysis/blind_test_e1.py", "blind test: alcohols and sugars"),
    ("analysis/blind_test_sorbitol_pg.py", "sorbitol (Baer, Siniti) and PG"),
    ("analysis/shepard_test.py", "glycerol + NaCl ternary vs Shepard 1976"),
    ("analysis/sorbitol_25C.py", "sorbitol at 25 C; sorbitol + NaCl cross term"),
    ("analysis/sucrose_T_transfer.py", "sucrose at two temperatures"),
    ("analysis/osmotic_virial_benchmark.py", "osmotic virial benchmark (Zielinski 2014)"),
    ("analysis/glycerol_sparing.py", "recipe: glycerol-sparing curve (prediction file checked, not replaced)"),
    ("analysis/recipe_sensitivity.py", "recipe sensitivity to salt target"),
    ("analysis/run_recipe.py", "old weighted-cost recipe demo (METHOD DEMO)"),
    ("analysis/freeze_concentration_cells.py", "freeze concentration, equilibrium cell model, three-criteria recipe"),
    ("analysis/posthoc_per_molecule.py", "POST HOC: per-molecule constants, Gibbs-Helmholtz split"),
    ("analysis/predict_disaccharides.py", "re-create maltose/lactose blind predictions (check only)"),
    ("analysis/blind_test_disaccharides.py", "blind test #3: maltose, lactose"),
    ("analysis/schematic.py", "Fig. 1 schematic (illustration)"),
    ("analysis/shepard_overlay.py", "Shepard 1976 Fig. 3 with blind predictions overlaid"),
    ("analysis/report_figures.py", "report figures (pipeline, results)"),
    ("analysis/worked_calculations.py", "calculations workbook values (hand = code asserted)"),
    ("analysis/uncertainty_budget.py", "lab uncertainty budget (preview until inputs are measured)"),
]


def main():
    tmp = Path(tempfile.mkdtemp(prefix="cryo_sle_run_"))
    work = tmp / "repo"
    shutil.copytree(ROOT, work, ignore=shutil.ignore_patterns(".git", "__pycache__", "*.zip"))
    failed = []
    print(f"Running in a temporary copy: {work}\n")
    for script, what in STEPS:
        t0 = time.time()
        r = subprocess.run([sys.executable, script], cwd=work, capture_output=True, text=True)
        status = "ok" if r.returncode == 0 else "FAILED"
        print(f"[{status:6s}] {time.time()-t0:5.1f}s  {script:42s} {what}")
        if r.returncode != 0:
            failed.append(script)
            print("         " + r.stderr.strip().splitlines()[-1] if r.stderr.strip() else "")

    print("\nCommitted predictions (never overwritten; regenerated copy compared byte-for-byte):")
    for name in PREDICTIONS:
        a, b = ROOT / "predictions" / name, work / "predictions" / name
        same = a.exists() and b.exists() and filecmp.cmp(a, b, shallow=False)
        note = "" if same else ("  (expected: " + EXPECTED_DIFF[name] + ")" if name in EXPECTED_DIFF else "  <-- UNEXPECTED, investigate")
        print(f"  {'reproduces' if same else 'DIFFERS   '}  {name}{note}")

    print("\nCopying back figures and evaluation results:")
    n = 0
    for f in sorted((work / "figures").glob("*.png")):
        shutil.copy2(f, ROOT / "figures" / f.name); n += 1
    for f in sorted((work / "predictions").glob("*.csv")):
        if f.name not in PREDICTIONS:
            changed = not (ROOT / "predictions" / f.name).exists() or not filecmp.cmp(f, ROOT / "predictions" / f.name, shallow=False)
            shutil.copy2(f, ROOT / "predictions" / f.name)
            print(f"  {'updated  ' if changed else 'unchanged'}  predictions/{f.name}")
    print(f"  {n} figures written to figures/")
    shutil.rmtree(tmp, ignore_errors=True)
    if failed:
        print(f"\n{len(failed)} step(s) failed: {', '.join(failed)}"); sys.exit(1)
    print("\nAll steps ran.")


if __name__ == "__main__":
    main()
