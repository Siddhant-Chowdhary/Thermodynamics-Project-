"""Physical constants and solute data.

Water properties are standard values. Solute molar masses follow from molecular
formulae. Anything used for a quantitative claim in the report should still be
checked against, and cited from, a primary source.
"""

R = 8.314462618          # J mol^-1 K^-1

# ---- water ----
T_M_WATER = 273.15       # K, normal melting point
DH_FUS_WATER = 6009.5    # J mol^-1, enthalpy of fusion at T_m
DCP_FUS_WATER = 37.9     # J mol^-1 K^-1, Cp(liquid) - Cp(ice) near 0 C
M_WATER = 18.015         # g mol^-1
V_WATER = 18.07e-6       # m^3 mol^-1, molar volume of liquid water

# ---- solutes ----
# molar_mass: g/mol, n_OH: hydroxyl groups, nu: particles per formula unit in
# solution (2 for fully dissociated NaCl, treated as ideal ions).
SOLUTES = {
    "propylene_glycol": {"formula": "C3H8O2",   "molar_mass": 76.094,  "n_OH": 2, "nu": 1},
    "glycerol":         {"formula": "C3H8O3",   "molar_mass": 92.094,  "n_OH": 3, "nu": 1},
    "sorbitol":         {"formula": "C6H14O6",  "molar_mass": 182.172, "n_OH": 6, "nu": 1},
    "ethylene_glycol":  {"formula": "C2H6O2",   "molar_mass": 62.068,  "n_OH": 2, "nu": 1},
    "mannitol":         {"formula": "C6H14O6",  "molar_mass": 182.172, "n_OH": 6, "nu": 1},
    "sucrose":          {"formula": "C12H22O11", "molar_mass": 342.297, "n_OH": 8, "nu": 1},
    "NaCl":             {"formula": "NaCl",     "molar_mass": 58.443,  "n_OH": 0, "nu": 2},
}


def solute(name):
    try:
        return SOLUTES[name]
    except KeyError as exc:
        raise KeyError(f"Unknown solute '{name}'. Known: {sorted(SOLUTES)}") from exc
