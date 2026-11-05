# CLAUDE.md — cryo-sle (CL207 optional project, IIT Bombay)

Owner: Siddhant (Sid) Chowdhary, 25B0347. Course: CL207 Chemical Engineering Thermodynamics.
Submission deadline: **6 Nov 2026**. Checkpoint 3 with the professor: ~27 Oct. Marks: 10 bonus + best-student prize consideration.
Prize criterion: "clear understanding of chemical engineering thermodynamics, applied through the project".

## The project in one paragraph
When cells freeze, ice is pure water, so salt concentrates in the shrinking liquid and damages cells. Cryoprotectants
(glycerol, sorbitol, propylene glycol) dilute the salt but are toxic. We model freezing with solid–liquid equilibrium:
ln(x_w γ_w) = −(ΔH_fus/R)(1/T − 1/T_m) − (ΔCp/R)(T_m/T − 1 + ln(T/T_m)). γ_w comes from a one-constant excess-Gibbs model
whose water–solute interaction is written **per –OH group**: A_i = n_OH,i × a_OH. Fit a_OH once on glycerol; predict PG,
sorbitol and mixtures blind; then design the lowest-toxicity glycerol + sorbitol recipe that keeps residual NaCl below a target.

## Rules (non-negotiable)
1. **Never invent data or citations.** Every number in data/ has a source and page. Unknown = clearly marked PLACEHOLDER.
2. **Blind tests stay blind.** Never refit a_OH on PG/sorbitol/mixture data. Never edit files in predictions/ after they are
   committed; new predictions go in a new versioned file. Commit + push predictions BEFORE comparing to any new data.
3. **Ask Sid before** any modelling assumption that changes conclusions, any change to research objectives, or choosing
   between viable options. Drafting, refactoring, plotting and tests don't need approval.
4. **Verify literature claims before stating them.** Say "not found in this form", never "never done".
5. Bootstrap CIs on a_OH reflect data scatter only, not model error. Don't present them as prediction uncertainty.
6. Figures: COLOUR since 5 Nov (Sid): Okabe-Ito palette in analysis/style.py, one colour per model everywhere; keep marker
   shapes/line styles so greyscale prints still read. Units on every axis. Placeholder-based outputs
   must carry a visible "METHOD DEMO" stamp.
7. Commit after every meaningful step with a clear message. Run `python -m pytest -q` before committing.
8. Sid prefers brief, direct answers: plain explanation first, then technical. Push back on weak reasoning.

## Decisions already made (don't re-litigate)
- Solutes: propylene glycol (2 OH), glycerol (3 OH), sorbitol (6 OH). DMSO/trehalose dropped.
- Fit set: CRC "Concentrative Properties" glycerol, 2–40 wt% (points < 2 wt% dropped: rounding). Lange/Lane-1925 values at
  10–20 wt% are suspect (warmer than ideal) and not used for fitting.
- NaCl + polyols: **additive** in ln a_w at equal molality, no cross term (Kleinhans & Mazur 2007 support). NaCl via empirical
  φ(m) = 1 − a√m/(1+√m) + b m + c m² fitted to the NaCl ice liquidus (thermo/salt.py). Shepard 1976 ternary data is a TEST, not a fit.
- Recipe (revised 21 Oct): glycerol-sparing curve = min glycerol for each sorbitol dose s.t. residual NaCl <= 4.41 wt% at -15 C
  (sweep -10/-20 C), 0.9 wt% saline base; no toxicity weights; zero-glycerol end unacceptable (sorbitol non-penetrating). The old
  weighted-cost optimiser (thermo/optimise.py) and the degenerate fixed-dose result stay as a documented side finding.
- Lab is optional: one ~3 h session requested (glycerol 30, PG 30, sorbitol 40, glycerol 15 + sorbitol 15 wt%, duplicates).
  The project must stand on lab-free tests if the lab doesn't happen.
- Novelty is narrow. Closest prior work: Zargarzadeh & Elliott, J Phys Chem B 123, 1099 (2019) — Margules on water/glycerol SLE. Must cite.

- **Decided 21 Oct (pre-registered, before any salt data):** primary mixing rule = additive in ln a_w (follows directly from the
  SLE/activity framework of the course); Kleinhans–Mazur freezing-point summation reported as a comparison only, never switched to.
  Shepard 1976 pass criterion: RMS <= 0.5 K for total solute <= 30 wt%; above that, report where it fails. Lab session 2 salt sample:
  glycerol 20 wt% + NaCl 5 wt% (predictions/salt_set_v1.csv). NaCl data used provisionally (unverified) until the printed CRC check.

## Results so far (real, reproducible)
| Item | Value | File |
|---|---|---|
| a_OH | −2414.9 J/mol per OH (bootstrap 95% CI −2461 to −2395) | predictions/predictions_v1.csv header |
| Fit quality, CRC glycerol 2–40 wt% | RMSE 0.017 K, max 0.041 K | analysis/fit_literature.py |
| UNIFAC on same glycerol data | RMSE 2.1 K (+5.8 K at 40 wt%); Dortmund 2.4 K (+6.9 K). Both give γ_w > 1 (wrong sign). Aspen confirmation pending | analysis/unifac_compare.py |
| NaCl φ fit | a,b,c = 0.2568, 0.0500, 0.0066; RMSE 0.023 K to eutectic | thermo/salt.py |
| No-CPA residual NaCl at −10 °C | 14.0 wt% | analysis/run_recipe.py |
| Recipe | Glycerol-sparing curve (see below); salt target and design T sourced | analysis/glycerol_sparing.py |

**Blind test #1 (21 Oct, CRC 2005, no refitting):** ethylene glycol (2 OH, no CH3), 2-40 wt%: RMSE ours 0.83 K vs ideal 2.49,
UNIFAC 3.99, Dortmund 3.43 (ours 3-5x better; under-predicts depression above ~30 wt%, worst +2.7 K at 40 wt%). D-mannitol (6 OH),
2-15 wt%: all models within 0.03-0.26 K; inconclusive (non-ideality too small at its solubility limit). Glycerol CRC rows verified exactly;
NaCl now CRC 2005 (17 rows to 22 wt%), 23.3 wt% Wikipedia row excluded; NaCl refit a,b,c = 0.2602, 0.0523, 0.0061.
Source: CRC Handbook of Chemistry and Physics, Internet Version 2005, D. R. Lide (ed.), section 8 pp. 8-58 to 8-84.

**Extension E1 blind test (21 Oct, CRC 2005, predictions committed first, no refit):** RMSE ours / ideal, 2-40 wt%:
methanol 2.20/4.84, ethanol 2.57/4.12, 1-propanol 0.89/0.50, 2-propanol 0.36/0.56, glucose 0.12/0.23, fructose 0.14/0.15,
sucrose 0.02/0.35 (UNIFAC/Dortmund worse for all alcohols). Up to 20 wt% per-OH beats ideal for every molecule except 1-propanol (tie).
Post-hoc per-molecule A/n_OH (2-20 wt%): EG -2473, glycerol -2426, sucrose -2127 (transfer holds); mannitol -915, glucose -1427,
fructose -1058 (weakly determined, small depressions); monoalcohols -4113 (MeOH), -4980 (EtOH), -2413 (1-PrOH), -6495 (2-PrOH):
scattered and non-monotonic in alkyl size => a linear CHx group term is NOT supported; two-group model dropped. Scope statement:
per-OH transfer works for polyols with ~one OH per carbon; it does not describe monohydric alcohols. PG (2 OH, 3 C, CH3) is in between.

**Sorbitol / PG tests (21 Oct, committed a_OH, no refit):** sorbitol 1-11 wt% (Baer & Baldwin, J Dairy Sci 67, 2860, 1984,
cryoscope; their sucrose matches CRC within 0.02 K): RMSE ours 0.07 K, ideal 0.12, UNIFAC 0.17, Dortmund 0.24. Sorbitol 53-64 wt%
(Siniti et al., Thermochim Acta 325, 171, 1999, DSC 'tentative'): ours over-predicts depression by 3.5-5.5 K (UNIFAC best, 1.6 K);
x=0.05 point suspect. Consistent with post-hoc A/n_OH ~ -1000 to -1500 for C6 polyols: per-OH likely overstates for large polyols at
high concentration. PG vs Zielinski 2014 osmotic-virial correlation (fitted to PG data, not raw data): ours under-predicts depression by
0.07/0.37/1.16/2.95 K at 10/20/30/40 wt%; ideal 0.2/1.0/2.9/6.5 K; UNIFAC worse. Sorbitol eutectic ~ -40 C at x 0.16-0.20 (66-71 wt%):
40 wt% lab sample is safely in the ice region.

**Shepard 1976 ternary test (22 Oct, Fig. 3 digitised by Sid, data/digitised/shepard1976_fig3.csv, analysis/shepard_test.py, no refit):**
pre-registered (solute <= 30 wt%): primary RMS 0.69 K -> MISSES 0.5 K target (K-M summation 0.59, ideal poly 0.89, ideal all 1.76).
Post hoc, labelled: figure noise floor (digitised edge curves vs CRC binaries) 0.34 K; per-curve biases non-monotonic in R (digitising).
Solute > 30 wt%: primary 1.23 K vs summation 3.39, ideal poly 4.64, ideal all 5.99. 89% of ternary points within Shepard's +-1.5 K
(summation 76%). Report: narrow miss of the strict target; additive-ln a_w is best by ~3x where cryo matters (> 30 wt%).
Don't call it inconclusive (Sid, 22 Oct); never re-state the criterion as if pre-registered. Salt target sourced: 4.41 wt% NaCl (Lovelock 1953 via Shepard 1976 p. 15).

**Sorbitol 25 C test (22 Oct, Bower & Robinson, J Phys Chem 67, 1540, 1963, NBS isopiestic; data/bower1963_sorbitol_25C.csv;
analysis/sorbitol_25C.py; committed a_OH, A const in T as fixed in PROMPTS.md 17 Oct; no refit):** data phi-1 = 0.008-0.074 for
0.9-3.7 mol/kg. Per-OH predicts 0.085-0.31 (about 4-5x too non-ideal); implied a_OH at 25 C = -413 to -744 J/mol (vs -2415).
RMS error in phi: per-OH 0.163, alpha-fixed 0.183, ideal 0.060, UNIFAC 0.113, Dortmund 0.716. => per-OH with glycerol's a_OH does
NOT transfer to sorbitol at 25 C. Can't yet separate T effect from molecule effect (needs glycerol 25 C osmotic data).
**Baer 1984 sorbitol flagged:** its implied phi is 1.14 even at 0.055 mol/kg (must -> 1 at infinite dilution), i.e. a ~13%
systematic excess depression (likely sample/purity or concentration basis). Our earlier 'sorbitol 1-11 wt% works' win over ideal is
NOT trustworthy; report Baer as data-quality-flagged. Consistent picture now: per-OH overstates C6 polyols (sorbitol 25 C, Siniti,
mannitol A/n_OH -915). Works for EG/glycerol (C2-C3), sucrose. Recipe replacement ratio still stands because it was shown
across per-OH / ideal / UNIFAC (ideal is closest to 25 C sorbitol data).
Ternary (Table II): neglected sorbitol-NaCl cross term (salting-in, Delta<0) = 0.6-5% of osmolality; as freezing-point error
0.01-0.15 K for NaCl <= 2 mol/kg, up to 1.2 K at 5.8 mol/kg. At the recipe's -15 C residual liquid (eq. 1, partly extrapolated):
0.05-0.13 K for recipes with glycerol. Supports the additive mixing rule where the recipe operates.

**Sucrose at two temperatures (22 Oct, analysis/sucrose_T_transfer.py; CRC freezing points + Robinson & Stokes, Trans Faraday Soc
45, 612, 1949, Table IV, 25 C; no refit):** per-OH fits sucrose on the freezing line (RMS err in phi 0.010-0.014; implied a_OH
-2160 to -2270) and is moderately too non-ideal at 25 C (0.029-0.090; implied a_OH -1640 to -2010). => a real but modest
temperature effect (~20-25% weaker interaction at 25 C, beyond the 1/T already in alpha). Sorbitol at 25 C implies -413 to -744,
i.e. 3-4x weaker than sucrose at the SAME T => the sorbitol failure is mainly a MOLECULE effect, not temperature.
(The 1949 paper omits glycerol data; it points to Scatchard, Hamer & Wood, JACS 60, 3061, 1938, not obtained.)

**Recipe sensitivity (22 Oct, analysis/recipe_sensitivity.py -> predictions/recipe_sensitivity.csv):** salt target +-25%
(3.31/4.41/5.51 wt%) x T (-10/-15/-20 C) x 3 models. Glycerol dose moves a lot (-15 C per-OH: 9.2 / 6.5 / 4.9 wt%), but the
replacement ratio barely moves: 0.68-0.91 g glycerol saved per g sorbitol over ALL 27 cases. Headline design result = the ratio.

**Osmotic virial benchmark (22 Oct, analysis/osmotic_virial_benchmark.py; Zielinski 2014 Table 3 coefficients; OV converted to
T with Zielinski's own relation, ds = 22.00):** binaries (OV is in-sample, fitted per molecule to CRC): RMSE per-OH / OV / ideal:
glycerol 0.02/0.27/1.31, EG 0.83/0.52/2.49, sucrose 0.02/0.02/0.35, glucose 0.12/0.01/0.23, mannitol 0.06/0.02/0.03,
MeOH 2.20/1.29/4.84, EtOH 2.57/0.91/4.12. Fair head-to-head = Shepard ternary (both predict from binaries): <=30 wt% solute
per-OH 0.69, OV 0.69, K-M 0.59; >30 wt%: per-OH 1.23, OV 5.75 (2.66 if put through our SLE), K-M 3.39.
=> one transferable constant matches per-molecule fits at low concentration and beats both cryobiology models when concentrated.

**Recipe redefined 21 Oct (Sid approved): glycerol-sparing curve.** Min glycerol to keep residual NaCl <= 4.41 wt% in 0.9 wt% saline at
-15 C (Gao & Critser lethal zone; sweep -10/-20 C). No toxicity weights. -15 C, no sorbitol: glycerol 6.5 (per-OH) / 7.6 (ideal) /
8.4 (UNIFAC) wt%. Replacement ratio robust across models: 0.74-0.86 g glycerol saved per g sorbitol. Zero-glycerol end marked unacceptable.
Proposed lab point glycerol 2.75 + sorbitol 5 + NaCl 0.9 wt%: predicted Tf -1.79/-1.73/-1.65 C (models differ by only 0.14 K, so a
freezing-point measurement barely tests the recipe; a residual-salt measurement would). predictions/recipe_v2_glycerol_sparing.csv.

Committed blind predictions: predictions_v1.csv (16 Oct), unifac_v1.csv, lab_set_v2.csv (18 Oct).

## Code map
- thermo/ — constants, composition, activity (IdealModel, MargulesModel, GroupContributionModel), sle (freezing_point,
  residual_liquid; uses model.ln_aw(...) if the model defines it), fit (least squares, bootstrap, Gibbs–Duhem), optimise
  (design_recipe with min_a option; optimise_ratio kept to document the degenerate formulation), salt (NaCl φ, AdditiveSaltModel).
- analysis/ — check_literature, fit_literature, unifac_compare, plot_model_comparison, predict_lab_set, run_recipe,
  cooling_curve (find_freezing_point), simulate_curve.
- arduino/ — ds18b20_logger.ino (3 probes, one data pin, 4.7 kΩ pull-up), serial_to_csv.py.
- data/ — literature.csv (quality-flagged; verified_vs_original = no until checked against print), LITERATURE_NOTES.md,
  design_inputs.csv, raw/ (lab CSVs: YYYYMMDD_runN_solute_wt.csv).
- **Gotcha:** the pip package `thermo` (Caleb Bell, used for UNIFAC) has the same name as our local package. Scripts that use
  UNIFAC (analysis/unifac_compare.py, plot_model_comparison.py, predict_lab_set.py) must be run as scripts from the repo root
  and must not import the local thermo. Use include_dcp=True everywhere for results.

## Commands
```
pip install -r requirements.txt
python run_all.py                        # rebuild ALL figures + result tables in a temp copy; checks committed predictions reproduce
python -m pytest -q                      # 14 tests
python analysis/fit_literature.py        # fit + predictions_v1 (do NOT rerun into the same file after new data exists)
python analysis/unifac_compare.py
python analysis/plot_model_comparison.py
python analysis/run_recipe.py
```

## Report must-includes (from Sid)
- Introduction: quote, in quotation marks, cited to Pegg, Methods Mol Biol 368, 39 (2007) (abstract): "Cryoprotectants, simply by
  increasing the total concentration of all solutes in the system, reduce the amount of ice formed at any given temperature; but to be
  biologically acceptable they must be able to penetrate into the cells and have low toxicity." Link it to freeze concentration and to the
  minimum-glycerol constraint (sorbitol is non-penetrating).
- Ours vs UNIFAC paragraph (Sid, 22 Oct): UNIFAC = general group contribution, many parameters fitted mostly on VLE near/above room T,
  nothing fitted to this system; gives gamma_w > 1 (wrong sign) for water-glycerol, under-predicts depression by up to ~6 K.
  Ours = one group (-OH), one constant fitted on glycerol freezing data; right sign and size for EG/glycerol/sucrose; fails for
  sorbitol and monoalcohols. Trade-off: breadth vs accuracy in a narrow class. Include the UNIFAC numbers from the results tables.

## Still to do (in order)
1. Push to GitHub (private now, public before submission).
2-4. DONE (library data entered; sorbitol/PG tests; Shepard ternary test).
5. DONE 22 Oct: sorbitol 25 C test (fails) + sucrose two-T test (sorbitol failure is mainly molecule, modest T effect).
6. DONE: recipe with sourced inputs (21 Oct) + sensitivity (22 Oct). Osmotic virial benchmark DONE 22 Oct.
7. Lab analysis (if lab happens): calibrate probes on water, find_freezing_point, parity plot vs lab_set_v2.csv.
8. Aspen cross-check of UNIFAC numbers.
9. run_all.py DONE 22 Oct (all 18 steps run; every result table reproduces byte-for-byte; 6/8 prediction files reproduce,
   2 differ for documented reasons). Still: final figure set; report (LaTeX report_paper.tex in the separate cryo-report folder); README.
10. (Sid, 22 Oct) Post-hoc section, clearly LABELLED post hoc and kept separate from blind tests: per-molecule A/n_OH refit to show
    what was learned after the sorbitol failure (scope: works for EG, glycerol, sucrose; fails for sorbitol, monoalcohols).
11. (Sid, 22 Oct, after final report) Colab notebook: clone repo, pip install, run run_all.py, show figures; Sid runs it himself
    and keeps screenshots / shows prof live, as evidence the code runs end to end.

## Done 23 Oct
- report/report_paper.tex/.pdf: Checkpoint 3 draft, 5 pages in prof's structure (abstract 199 words). report/checkpoint3_summary.pdf: one-page handout
  with 5 feedback questions for the prof. data/NOVELTY_SEARCH.md: closest prior work = Marcolli & Peter 2005 (polyol UNIFAC at 298 K),
  Peres & Macedo 1997 (OHring proximity group, sugars), Powell-Palm & Consiglio 2026 review (suggests GC screening, doesn't apply it).
  analysis/uncertainty_budget.py + data/lab/uncertainty_inputs.csv: preview U(k=2) ~0.25 K per sample (ASSUMED inputs, METHOD DEMO);
  supercooling correction ~ +0.2-0.3 K for 2 K supercooling must be applied to the measured peak.
- Left for Sid: Arduino build + ice-water calibration + dry run; DWSIM UNIFAC check; friend cold-read; GitHub upload; Checkpoint 3.

## Done 30 Oct (no lab; Sid approved: freeze-concentration, cell example, post hoc, Gibbs-Helmholtz, web tool, journal restyle)
- analysis/freeze_concentration_cells.py: unfrozen fraction + salt outside/inside an equilibrium cell (cell solutes = isotonic NaCl,
  glycerol penetrates, sorbitol/NaCl excluded). -15 C: saline 3.9% water liquid, 18.8 wt% NaCl; sorbitol 9 wt% -> outside 4.28 but
  inside 18.8 (= saline); glycerol 6.5 -> 4.43 both. Recipe under 3 criteria (predictions/recipe_three_criteria.csv): outside salt
  +0.76 g gly saved / g sorb; inside salt -2.10 (sorbitol COSTS glycerol); unfrozen >= 25% (Mazur 1981) +0.66; all three: optimum
  sorbitol 2.06 + glycerol 10.86 wt% vs glycerol alone 12.22 (11% less). HEADLINE CHANGE: sparing ratio depends on injury criterion.
- analysis/posthoc_per_molecule.py (POST HOC): per-molecule a <=20 wt%: EG -2473, gly -2426, suc -2127, glc -1392, fru -1086,
  man -915, sorb(conc) -1269, MeOH -4113, EtOH -4980, 1-PrOH -2413, 2-PrOH -6495. Gibbs-Helmholtz on sucrose: h = -6.3 to -7.6 kJ/mol,
  s = -15 to -20 J/mol/K per OH.
- predictions/lab_set_v3_glycerol_nacl.csv: glycerol 25 + NaCl 5 committed 30 Oct before any measurement (-11.71 C).
- web/index.html -> artifact "Cryo Recipe Explorer" (JS model matches Python to all printed digits at -15 C).
- report/report_final.tex/.pdf: journal two-column, 5 pages (abstract 194 words). run_all: all steps run, all results reproduce.

## Checkpoint 3 must-mention (Sid, 22 Oct)
- How the Shepard data was obtained: no tables in the paper; Fig. 3 hand-digitised by Sid (WebPlotDigitizer, 600 dpi, 9 isopleths,
  76 points); calibration vs NaCl-water eutectic (77.2 wt%, -21.6 C read vs 76.7, -21.1); 4 duplicate clicks dropped by fixed rule;
  figure noise floor 0.34 K vs CRC binaries; nothing refitted. Shows the effort behind the ternary test.

## Done 5 Nov (container was reset; repo restored from cryo-sle-full-30oct.zip, git history before that is lost)
- Re-applied 30 Oct work from the session log: report cost/time/feasibility paragraph, calculations workbook
  (analysis/worked_calculations.py + report/calculations_workbook*), blind test #3 maltose/lactose
  (predictions/disaccharides_v1.csv regenerated by the committed script BEFORE the CRC rows were re-entered).
  Maltose to 40 wt%: RMSE per-OH 0.029 K vs ideal 0.229 (blind). Lactose (to 8 wt%) too dilute to discriminate.
- All report figures recoloured (analysis/style.py), Fig. 1 schematic now analysis/schematic.py (vector PDF).
- report_final.tex: Methods 2.1 rewritten as 'Thermodynamic framework' (eqs 2-6: G^ex, ln gamma_w and ln gamma_i with
  Gibbs-Duhem, NaCl phi, additive mixture eq, K_f check, osmotic coefficient, inversion for a, Gibbs-Helmholtz);
  explicit heat-melt-cool expressions in intro; maltose row in Table 2; abstract 199 words; 5 pages (guideline: strict).
