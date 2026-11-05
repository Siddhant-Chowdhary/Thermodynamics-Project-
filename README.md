# The Ice Decides

CL207 (Chemical Engineering Thermodynamics) Project
Siddhant Chowdhary, 25B0347.

The report is [`report/report_final.pdf`](report/report_final.pdf). This repo has all the code and data behind it.

## What this project is about

The theme involves using thermodynamic concepts in biological systems via cryoprotectants.
When cells are frozen slowly, ice forms outside them. Ice is pure water, so all the salt gets squeezed into
the little liquid that is left, and that salty liquid damages the cells. Glycerol protects them but it is toxic,
so you want to use as little as possible.

I used solid–liquid equilibrium from the course to predict how salty that liquid gets. The only unknown is the
activity coefficient of water. I modelled it with one number per OH group, fitted that number on glycerol only,
and then tested it on other molecules and on glycerol + salt mixtures without changing it.

## **Try yourself: https://the-ice-decides.streamlit.app**

I wanted people to see the thermodynamics, not just read it, so I built a small app where you play cryobiologist.
Slide glycerol, sorbitol and temperature and watch the ice take over: how much water stays liquid and how salty it
gets outside and inside the cell. Try 9 wt% sorbitol with no glycerol and see the outside stay safe while the
inside turns to brine. It runs the same code as the report, so every number is real. Also, the graphs are interactive and keep changing as we change the conditions, so sure do have a go!! 

**Easiest: in Google Colab (nothing to install)**

Open this link and click Runtime > Run all:
https://colab.research.google.com/github/Siddhant-Chowdhary/Thermodynamics-Project-/blob/main/run_in_colab.ipynb

The notebook is the file `run_in_colab.ipynb` in the main folder of this repo. It downloads the repo, installs
everything, rebuilds all the figures and tables, runs the tests and shows the figures. It takes a few minutes.
You don't need to upload any files as all the files are extracted from this github repo.


**Interactive version:** `streamlit run app.py` opens a small page where you pick glycerol, sorbitol and
temperature and see how much water stays liquid and how salty it gets inside and outside the cell.

## Where each figure in the report comes from

All scripts are in `analysis/`. Figures go to `figures/`, result tables to `predictions/`.

| Report | Script |
|---|---|
| Fig. 1 (freezing sketch) | `schematic.py` |
| Fig. 2 (workflow) | `report_figures.py` |
| Fig. 3 and Table 2 (blind tests) | `fit_literature.py`, `blind_test_crc.py`, `blind_test_e1.py`, `blind_test_disaccharides.py`, `blind_test_sorbitol_pg.py`, `unifac_compare.py`, `shepard_test.py`, `osmotic_virial_benchmark.py`, `posthoc_per_molecule.py`, then `report_figures.py` draws it |
| Fig. 4 (Shepard's phase diagram with my predictions on top) | `shepard_overlay.py` |
| Fig. 5 (salt outside and inside the cell, recipes) | `freeze_concentration_cells.py` |
| Fig. 6 (enthalpy and entropy from sucrose) | `posthoc_per_molecule.py`, `sucrose_T_transfer.py` |
| Fig. 7 (sorbitol at 25 °C) | `sorbitol_25C.py` |
| Table 1 (planned lab samples) | `predict_lab_set.py`, `predict_salt_set.py`, `uncertainty_budget.py` |
| Every equation worked out with numbers | `worked_calculations.py` (it checks my hand calculation against the code) |

## How I kept the tests honest

1. I fitted the one constant on glycerol only (`fit_literature.py`). It is saved at the top of
   `predictions/predictions_v1.csv`.
2. Before looking at any test data, I saved the predictions in `predictions/` and committed them. I never edited
   those files afterwards.
3. My first upload of this repo (21 Oct 2026) already has `predictions_v1.csv`, `unifac_v1.csv`,
   `crc_extra_v1.csv`, `lab_set_v2.csv` and `salt_set_v1.csv`, and they are the same as now. Later rounds are
   `e1_alcohols_sugars_v1.csv`, `lab_set_v3_glycerol_nacl.csv`, `disaccharides_v1.csv` and
   `uraji_polyols_v1.csv`.
4. Anything I fitted after seeing results is labelled "post hoc".

## Data

`data/literature.csv` has every literature number I used, with its source and page. Most are from the CRC
Handbook. Shepard et al. (1976) only give their glycerol + NaCl data as a graph, so I read the points off the
graph myself with WebPlotDigitizer (`data/digitised/`). The 25 °C data are from Bower & Robinson (1963) and
Robinson & Stokes (1949).

## Folders

```
thermo/        the model: SLE solver, activity coefficients, NaCl, fitting
analysis/      one script per analysis (table above)
data/          literature data with sources, digitised graph
predictions/   saved blind predictions and test results
figures/       all figures
report/        the report and a worked-calculations PDF
arduino/       temperature logger for the planned freezing-point experiment
tests/         unit tests (python -m pytest -q)
```

## Note

There is a pip package called `thermo` (used only for the UNIFAC comparison) with the same name as my `thermo/`
folder. Run scripts from this main folder, the way `run_all.py` does, and it works.

The results are for slow cooling (equilibrium). They tell you what conditions a cell sees, not whether it
survives.
