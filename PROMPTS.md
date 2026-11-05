# Prompts for Claude Code — paste one at a time, in order

How to work with it:
- Start every session with prompt 0, so it reads CLAUDE.md and the current state.
- For big tasks (3, 4, 6, 9), say "plan first, show me the plan, then build". Approve the plan before it writes code.
- After each task: "run the tests, show me the diff summary, commit with a clear message".
- If it proposes a modelling assumption, it should ask you first (CLAUDE.md says so). If it doesn't, stop it.
- Anything in [brackets] is for you to fill in.

---

## 0. Orientation (every new session)
```
Read CLAUDE.md, README.md and git log. Run the tests. Then tell me in under 10 lines: what's done, what's next,
and anything in the repo that contradicts CLAUDE.md. Don't change any files.
```

## 1. Push to GitHub
```
Create a private GitHub repo called cryo-sle with gh and push all commits. Then show me the commit list with
dates as GitHub records them, and the URL of the predictions_v1.csv commit (this is the blind-test evidence).
```

## 2. Enter library data
```
I'm going to paste freezing-point tables from [source, edition, page]. Add them to data/literature.csv with
source, table_or_page, quality and phase_branch filled in. For the CRC glycerol rows that already exist, compare
against what I paste and set verified_vs_original=yes where they match; list any mismatches. Then run
analysis/check_literature.py and tell me if any new source looks internally inconsistent (non-smooth phi, values
warmer than ideal). Don't refit anything.

[paste table here]
```

## 3. Blind test against literature data (the core result)
```
Plan first. Write analysis/blind_test.py that compares the COMMITTED predictions (predictions_v1.csv for our model,
unifac_v1.csv for UNIFAC and Dortmund, plus the ideal solution) against PG and sorbitol freezing points in
data/literature.csv that were added after the predictions commit. Don't refit a_OH. Where literature
concentrations don't match predicted ones, evaluate each model at the literature concentration using the
committed a_OH from the predictions_v1.csv header; say clearly in the code that this is evaluation, not fitting.
Output: (a) a parity plot, measured vs predicted, one marker shape per model, y=x line, black-and-white;
(b) a table of RMSE and mean signed error per model per solute; (c) a short plain-text verdict on whether
per-OH transfer works for PG and for sorbitol, and in which direction it fails if it does.
Save to figures/blind_test_parity.png and predictions/blind_test_results.csv. Add a test that fails if
blind_test.py changes any file in predictions/ that existed before it ran.
```

## 4. Salt test against Shepard 1976
```
Plan first. Using the water–NaCl–glycerol freezing points from Shepard, Goldston & Cocks, Cryobiology 13, 9 (1976),
which I'll paste [paste data: composition + T], test the AdditiveSaltModel (per-OH glycerol + fitted NaCl phi) and
two baselines: ideal glycerol + fitted NaCl, and ideal everything. No refitting. Parity plot + RMSE table. Then
tell me at what glycerol:salt ratio or concentration the additive assumption starts to fail.
```

## 5. Temperature-transfer test (sorbitol at 25 °C)
```
Plan first. From the isopiestic sorbitol–water data at 25 °C [paste: molality, osmotic coefficient or a_w],
compute ln a_w predicted by our per-OH model with the committed a_OH. Our model has A independent of T, so
alpha = A/RT changes with T. Compare against the data and against the ideal solution. Explain in plain words
what agreement or disagreement means for the assumption that A is temperature-independent.
Figure + table, no refitting.
```

## 6. Recipe with real inputs
```
Update data/design_inputs.csv with these values and citations: [cap_glycerol = …, source; cap_sorbitol = …, source;
salt_target = …, source; T_design_C = …, justification; min_glycerol = …]. Set status to SOURCED for each one I
give a source for. Then run analysis/run_recipe.py. Add a sensitivity analysis: vary each cap and the salt target
by ±50% and T_design by ±5 °C, one at a time, and show how the optimum doses move (tornado-style table + one
figure). Compare the recipe from our model vs ideal polyols vs UNIFAC-based polyols to show whether the
thermodynamic model changes the design. Commit the recipe prediction as predictions/recipe_v2.csv BEFORE any lab
measurement of it.
```

## 7. Lab day
Before:
```
Walk me through calibrating the 3 DS18B20 probes in ice water using arduino/serial_to_csv.py. Write a small
script that reads a calibration CSV and stores each probe's offset in data/raw/probe_offsets.csv.
```
After each run:
```
New run files are in data/raw/ [list names]. Apply the probe offsets, run find_freezing_point on each sample
trace, plot each cooling curve with the detected freezing point marked, and add the results to
data/lab_results.csv (sample, composition by actual weighed masses, T_f, uncertainty, run id, notes).
Flag any duplicate pair that differs by more than 0.3 °C.
```
After the session:
```
Plan first. Write analysis/lab_parity.py comparing lab_results.csv against the committed lab_set_v2.csv (all four
models). Glycerol 30 wt% is the method check against CRC: report its offset and say whether to correct for it or
not, and why. Parity plot + error table, no refitting.
```

## 8. Aspen cross-check
```
I ran Aspen Properties with [UNIFAC / UNIF-DMD] for water + [solute] and exported [gamma_water vs x at T = …]
[paste or file path]. Compare against what analysis/unifac_compare.py computes with the thermo package at the same
T and x. Report the max difference in gamma_w and in the resulting freezing point. If they disagree by more than
0.1 K in T_f, find out why (group assignment, parameter table version) before we use either.
```

## 9. Final figures and reproducibility
```
Plan first. Make run_all.py that regenerates every figure and results table from data/ and predictions/ in one
command, without overwriting any committed prediction file. The final five figures: (1) freeze-concentration
schematic (exists in cryo-report), (2) glycerol fit with ideal and UNIFAC lines, (3) blind-test parity plot,
(4) residual NaCl vs storage T for saline-only, glycerol-only and the designed recipe, (5) recipe map with the
salt-target boundary and optimum. Same style everywhere: black-and-white, 8–9 pt fonts, units on axes, sized for
a two-column paper. Update README with one-command reproduction and the results table.
```

## 10. Report writing
```
Open report_paper.tex [path]. Replace placeholder results section by section, starting with Results, then
Discussion, Methods, Introduction, Abstract last. Every number you write must come from a file in the repo: put
a LaTeX comment next to it naming the file. If a number doesn't exist yet, write \todo{…} instead of inventing one.
Keep claims narrow: cite Zargarzadeh & Elliott 2019 as closest prior work. Plain, student-researcher tone,
no hype words.
```

## 11. Hostile review (do this before Checkpoint 3 and before submission)
```
Act as a strict thermodynamics professor grading this project for the best-student prize. Read the report draft,
CLAUDE.md and the results files. List the 10 weakest points: claims not backed by data, circular tests,
uncertainty not propagated, assumptions not stated, figures that mislead, and anything I can't defend in a viva.
For each one: severity, the exact location, and the smallest fix. Don't fix anything yet.
```

## 12. Viva prep
```
Generate the 15 hardest questions a professor would ask about this project, with short answers grounded in our
actual numbers and files. Include: why T alone fixes a_w with ice present; per-OH vs per-molecule Margules; why
UNIFAC fails here; how we know the predictions were blind; measurement uncertainty vs effect size; why the
fixed-dose optimisation is degenerate; what the model can't say about cell survival.
```
