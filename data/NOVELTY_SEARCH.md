# Novelty search record

Searched 23 Oct 2026 (web search plus the reference lists of papers already in hand). This is not exhaustive. The claim we can make is
"not found in this form in the literature we searched", never "never done". Before submission, Sid should repeat the
Google Scholar part ("Cited by" on Zargarzadeh & Elliott 2019 and Marcolli & Peter 2005) and add the date to this file.

## Search terms used
1. group contribution model water activity polyols hydroxyl group freezing point depression
2. Water Activity in Polyol Systems, J. Chem. Eng. Data (Ninni / Meirelles)
3. cryoprotectant mixture design glycerol sorbitol freeze concentration salt thermodynamic model
4. UNIFAC sugars polyols aqueous solutions modified hydroxyl group parameters (Peres & Macedo, Spiliotis & Tassios)
5. activity coefficient model cryobiology predict freezing point new cryoprotectant without measurement
6. hydroxyl group interaction parameter Margules transfer polyols freezing point cryoprotectant

## Closest hits and how this work differs

| Work | What it does | Difference from this project |
|---|---|---|
| Zargarzadeh & Elliott, J. Phys. Chem. B 123, 1099 (2019), doi:10.1021/acs.jpcb.8b09285 | 1- and 2-parameter Margules and osmotic virial models fitted to water/glycerol (and other) SLE | Fitted per molecule. No transfer between molecules, no blind tests, no design. **Our glycerol fit is not new and cites this paper.** |
| Marcolli & Peter, Atmos. Chem. Phys. 5, 1545 (2005), doi:10.5194/acp-5-1545-2005 | New UNIFAC parameterisation for polyols (C2–C7) from water-activity data at 298 K. Notes that standard UNIFAC overestimates a_w for polyols with 2 or more OH groups. | Closest group-contribution precedent. Theirs is a multi-parameter UNIFAC fitted at 298 K for atmospheric aerosols. Ours is one constant fitted on freezing data of one molecule, tested blind on freezing points and the cryobiology ternary. Their UNIFAC finding agrees with ours (UNIFAC gets water–glycerol wrong). |
| Peres & Macedo, Fluid Phase Equilib. 139, 47 (1997), doi:10.1016/S0378-3812(97)00196-9 | Modified UNIFAC for sugar solutions, adding an "OHring" group for proximity effects between hydroxyls | Supports our explanation of the sorbitol failure (neighbouring OH groups interact). Sugars, not cryoprotection. |
| Fysun et al., Chem. Ing. Tech. 87, 1327 (2015), doi:10.1002/cite.20140034 | Prediction of water activity in aqueous polyol solutions | Only the abstract page was reachable, so the method is unknown. Read it before claiming anything about it. |
| Elliott et al., J. Phys. Chem. B 111, 1775 (2007); Zielinski et al., Cryobiology 69, 305 (2014) | Osmotic virial equation; coefficients fitted for each solute; mixtures predicted from binaries | Needs data for every molecule. We benchmark against it (Shepard ternary). |
| Kleinhans & Mazur, Cryobiology 54, 212 (2007) | Synthesises ternary diagrams by summing binary freezing-point depressions | We benchmark against it. Our rule adds ln a_w instead. |
| Powell-Palm & Consiglio, arXiv:2601.12573 (2026), book chapter | Review of cryopreservation thermodynamics. Says group-contribution models (UNIFAC, AIOMFAC) could support predictive screening of CPA mixtures, but does not apply them to freezing points | Shows the idea is seen as useful but not yet done in this form. Good for motivation. |

## Conclusion (wording for the report)
Group-contribution descriptions of polyol water activity exist (Marcolli & Peter 2005; Peres & Macedo 1997), as do Margules
fits of water–glycerol SLE (Zargarzadeh & Elliott 2019) and per-solute cryobiology models (osmotic virial, Kleinhans–Mazur).
We did not find, in the literature searched, the following combination:
- one hydroxyl constant fitted on a single cryoprotectant's freezing data;
- that constant tested blind against freezing points of other solutes and the glycerol–NaCl ternary, using the cryobiology
  field's own models as benchmarks;
- the constant then used for a glycerol-sparing cryoprotectant design.
