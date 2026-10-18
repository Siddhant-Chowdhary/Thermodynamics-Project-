# Literature freezing-point data: status and notes (Phase 1)

Last updated: 2 Oct 2026. Data file: `data/literature.csv`. Check script: `python analysis/check_literature.py`.

Every value in the CSV was read from the source named in its row. Nothing is typed in from memory.
The `verified_vs_original` column is `no` for every row until it has been checked against the printed source.

## Coverage

| Solute | Rows | Source | Quality | Usable for |
|---|---|---|---|---|
| Glycerol | 21 (0.5–40 wt%) | CRC "Concentrative Properties of Aqueous Solutions" (unofficial online copy) | handbook, 0.01 °C | **fitting a_OH** after the print copy is checked |
| Glycerol | 6 ice-branch rows (10–60 wt%) + 4 glycerol-branch rows | Lange's 10th ed., via Wikipedia | secondary, 0.1 °C | cross-check only (see below) |
| NaCl | 18 (0.5–23.3 wt%) | Wikipedia NaCl data page (no citation on page; probably CRC) | secondary | check against CRC, then use |
| Propylene glycol | 6 (10–60 wt%) | Engineering ToolBox | **low**: rounded to 1 °C, no source cited | sanity check only, **not for the blind test** |
| Sorbitol | **0** | none found open-access | | **gap: needs the library** |

## Findings from the checks

1. **The CRC glycerol rows look internally sound.** The molality column matches wt% → mol/kg for the rows checked. The apparent osmotic coefficient φ rises smoothly from 1.00 at 2 wt% to 1.28 at 40 wt%. Below 2 wt%, φ looks wrong (0.69–0.88), but that comes from rounding: a 0.01 °C step is a large share of a 0.07 °C depression. Drop these points from the fit or give them little weight.
2. **Lange's (Lane 1925) low-concentration glycerol values are suspect.** At 10 wt% Lange gives −1.6 °C, which is warmer than the ideal value (−2.19 °C) and implies φ = 0.72. That is implausible for glycerol. At 20 wt%, Lange and CRC also disagree (−4.8 vs −5.46 °C). The two sources agree at 40 wt%. Above 40 wt%, Lange's φ is non-monotonic (1.28 → 1.27 → 1.41), so it needs checking against the original Lane 1925 table.
   → CRC covers only up to 40 wt%. For data above that, read Lane 1925 itself rather than a copy.
3. **NaCl is clearly non-ideal, and the current code treats it as ideal ions (ν = 2).** φ is about 0.92 at low concentration and about 1.26 at the eutectic. The ideal model gets T_f wrong by up to 0.18 °C below 10 wt% and by 5 °C at 23 wt%. The salt-dilution target depends directly on this, so NaCl needs its own activity correction, such as an osmotic-coefficient fit or interpolating the data directly. It cannot be left ideal.
4. **Early warning for the per-OH hypothesis.** This is based on low-quality data, so it isn't conclusive. Compared at similar molality, the Engineering ToolBox PG values imply φ ≈ 1.47 at 5.6 mol/kg, against glycerol's 1.22 at 5.1 mol/kg. The per-OH model makes non-ideality scale with n_OH, so it predicts PG (2 OH) should be *less* non-ideal than glycerol (3 OH). If good PG data shows the same thing, the CH3 group matters and a one-parameter OH model will under-predict PG. That would itself be a reportable result. Get proper PG data before Phase 2.

## Library list (IITB): what to pull, in priority order

1. **CRC Handbook**, "Concentrative Properties of Aqueous Solutions": check the glycerol rows; transcribe the NaCl, 1,2-propanediol and sorbitol/D-glucitol entries (if listed); record the edition and page numbers.
2. **Sorbitol** (the critical gap):
   - "Freezing Points of Bulking Agents Used in Manufacture of Low-Calorie Frozen Desserts", *J. Dairy Sci.* 1984 (ScienceDirect pii S0022030284816471).
   - "Study of the respective binary phase diagrams of sorbitol with mannitol, maltitol and water", *Thermochim. Acta* 1998 (pii S0040603198005760).
   - Uraji et al., *Food Sci. Technol. Int.* 2, 38–42 (1996). Open access, but the polyol data appear only in figures (30–60 wt%), so this is a last resort via digitising.
3. **Glycerol, full range:** L.B. Lane, *Ind. Eng. Chem.* 17, 924 (1925), doi:10.1021/ie50189a017.
4. **Propylene glycol:** ASHRAE Handbook *Fundamentals* (secondary-coolant properties chapter), or Dow "Dowfrost" engineering tables. For the cryo literature, look at Boutron & Kaufmann (1979), water–1,2-propanediol, *Cryobiology*.
5. **Must-read for novelty and benchmarking:** L. Zargarzadeh & J.A.W. Elliott, "Comparison of the Osmotic Virial Equation with the Margules Activity Model for Solid–Liquid Equilibrium", *J. Phys. Chem. B* 123, 1099–1107 (2019), doi:10.1021/acs.jpcb.8b09285. They fit one- and two-parameter Margules models to water/glycerol SLE. **Our glycerol Margules fit is therefore not new and must cite this paper.** The per-OH *transfer* to other polyols is still the distinct part.
6. **Ternary data with salt (could become a test set):**
   - Shepard, Goldston & Cocks, "The H2O–NaCl–glycerol phase diagram and its application in cryobiology", *Cryobiology* 13, 9–23 (1976), doi:10.1016/0011-2240(76)90154-1.
   - Kleinhans & Mazur, *Cryobiology* 54, 212–222 (2007), doi:10.1016/j.cryobiol.2007.01.007 (PMC2730495). They show that adding the binary glycerol and NaCl depressions reproduces the ternary diagram up to about 55 wt%. Their methods section lists the binary data sources they used.

## Pages that could not be read from here

- Melinder, KTH doctoral thesis (diva-portal): blocked. It tabulates PG, glycerol and NaCl freezing points.
- J. Dairy Sci. 1984 bulking-agents paper: HTTP 403.
- Kleinhans & Mazur full text (PMC / Europe PMC): captcha, then rate-limited.
- The CRC online copy is cut off alphabetically at 1-propanol, so its NaCl, PG and sorbitol entries were not reachable.
