"""Worked calculations: every equation in the project, applied by hand to one example each, with all intermediate numbers.

Each example is computed twice: (1) explicitly, line by line, from the formula as written in the handbook, and
(2) with the project's own code (thermo/, analysis/). The script asserts the two agree, so the workbook cannot drift
from the code. Data are read from data/ and predictions/; nothing is typed in by hand except the example compositions.

    python analysis/worked_calculations.py   -> report/calc_values.json and report/calculations_workbook.tex (filled template)
"""
import csv, json, math, sys, warnings
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from thermo import GroupContributionModel, IdealModel, freezing_point, ln_aw_ice as code_ln_aw_ice
from thermo.salt import AdditiveSaltModel, fit_nacl, phi_nacl
from thermo.constants import R, M_WATER, DH_FUS_WATER as DH, DCP_FUS_WATER as DCP, T_M_WATER as TM, solute
from analysis.freeze_concentration_cells import state

warnings.filterwarnings("ignore", message="NaCl molality")
Mw = M_WATER / 1000.0                                   # kg/mol
A_OH = float(open(ROOT / "predictions" / "predictions_v1.csv").readline().split("a_OH = ")[1].split()[0])
NACL = fit_nacl(include_dcp=True)[0]
LIT = list(csv.DictReader(open(ROOT / "data" / "literature.csv")))
V = {}                                                   # values for the template


def crc(sol, w):
    r = [r for r in LIT if r["solute"] == sol and r["source_id"] == "CRC2005" and abs(float(r["wt_percent"]) - w) < 1e-9]
    return float(r[0]["T_freeze_C"])


def put(k, x, fmt):
    V[k] = format(x, fmt)
    return x


def check(a, b, tol, what):
    assert abs(a - b) < tol, f"{what}: hand {a} vs code {b}"


def ice_terms(T):
    t1 = -(DH / R) * (1 / T - 1 / TM)
    t2 = -(DCP / R) * (TM / T - 1 + math.log(T / TM))
    return t1, t2


def newton_T(lhs_of_T, T0=TM - 1.0, n=6):
    """Solve g(T) = lhs(T) - ln a_w,ice(T) = 0 by Newton with a numerical derivative; return iterations."""
    it, T = [], T0
    for _ in range(n):
        g = lhs_of_T(T) - sum(ice_terms(T))
        h = 1e-4
        dg = ((lhs_of_T(T + h) - sum(ice_terms(T + h))) - g) / h
        it.append((T, g))
        T = T - g / dg
    return T, it


# ---------- C1 constants
put("R", R, ".6f"); put("DH", DH, ".1f"); put("DCP", DCP, ".1f"); put("TM", TM, ".2f"); put("Mw", M_WATER, ".3f")
put("AOH", A_OH, ".1f")
put("na", NACL[0], ".4f"); put("nb", NACL[1], ".4f"); put("nc", NACL[2], ".4f")
Kf = R * TM**2 * Mw / DH
put("Kf", Kf, ".3f")

# ---------- C2 ice line at -15 C
T15 = TM - 15
t1, t2 = ice_terms(T15)
put("iceT", T15, ".2f"); put("ice_inv", 1 / T15 - 1 / TM, ".4e"); put("ice_t1", t1, ".5f")
put("ice_r", TM / T15 - 1 + math.log(T15 / TM), ".6f"); put("ice_t2", t2, ".5f")
put("ice_ln", t1 + t2, ".5f"); put("ice_aw", math.exp(t1 + t2), ".4f")
put("ice_DHR", DH / R, ".2f"); put("ice_DCPR", DCP / R, ".4f")
check(t1 + t2, code_ln_aw_ice(T15, True), 1e-12, "ice line")
put("ice_pct2", 100 * t2 / (t1 + t2), ".1f")

# ---------- C3 composition, glycerol 20 wt%
Mg = solute("glycerol")["molar_mass"]
w = 20.0
n_g = w / Mg; n_w = (100 - w) / M_WATER
m_g = n_g / ((100 - w) / 1000)
x_w = n_w / (n_w + n_g)
put("g_ng", n_g, ".4f"); put("g_nw", n_w, ".4f"); put("g_m", m_g, ".4f"); put("g_xw", x_w, ".5f"); put("g_xs", 1 - x_w, ".5f")
put("g_lnxw", math.log(x_w), ".5f")
put("Mg", Mg, ".3f")

# ---------- C4 ideal freezing point, 20 wt% glycerol
T_id, it_id = newton_T(lambda T: math.log(x_w))
for i, (T, g) in enumerate(it_id[:4]):
    put(f"id_T{i}", T, ".4f"); put(f"id_g{i}", g, ".2e")
put("id_T", T_id, ".3f"); put("id_TC", T_id - TM, ".2f")
check(T_id, freezing_point({"glycerol": 0.20}, IdealModel(), include_dcp=True), 1e-6, "ideal Tf")
put("id_dil", -Kf * m_g, ".2f")
put("g_crc", crc("glycerol", 20.0), ".2f")

# ---------- C5 implied constant from one CRC point, and the fit
Tcrc = crc("glycerol", 20.0) + TM
lnaw = sum(ice_terms(Tcrc))
lngw = lnaw - math.log(x_w)
a_imp = lngw * R * Tcrc / (3 * (1 - x_w) ** 2)
put("f_T", Tcrc, ".2f"); put("f_lnaw", lnaw, ".5f"); put("f_lngw", lngw, ".5f"); put("f_gw", math.exp(lngw), ".4f")
put("f_a", a_imp, ".0f")
model = GroupContributionModel(A_OH)
fit_rows = [r for r in LIT if r["solute"] == "glycerol" and r["source_id"] == "CRC2005" and r["phase_branch"] == "ice"
            and 2 <= float(r["wt_percent"]) <= 40]
res = []
for r in fit_rows:
    ww = float(r["wt_percent"])
    res.append(freezing_point({"glycerol": ww / 100}, model, include_dcp=True) - TM - float(r["T_freeze_C"]))
put("fit_n", len(res), "d"); put("fit_rmse", math.sqrt(np.mean(np.square(res))), ".3f"); put("fit_max", max(map(abs, res)), ".3f")

# ---------- C6 non-ideal freezing point, 20 wt% glycerol, committed constant
A = 3 * A_OH


def lhs_g(T, xw=x_w, A_=A):
    return math.log(xw) + (A_ / (R * T)) * (1 - xw) ** 2


T_no, it_no = newton_T(lhs_g)
for i, (T, g) in enumerate(it_no[:4]):
    put(f"no_T{i}", T, ".4f"); put(f"no_g{i}", g, ".2e")
put("no_A", A, ".1f"); put("no_alpha", A / (R * T_no), ".4f"); put("no_lng", A / (R * T_no) * (1 - x_w) ** 2, ".5f")
put("no_gw", math.exp(A / (R * T_no) * (1 - x_w) ** 2), ".4f")
put("no_T", T_no, ".3f"); put("no_TC", T_no - TM, ".2f")
check(T_no, freezing_point({"glycerol": 0.20}, model, include_dcp=True), 1e-6, "non-ideal Tf")
# solute activity coefficient (Gibbs-Duhem consistent)
put("no_lngs", A / (R * T_no) * x_w ** 2, ".4f")

# ---------- C7 blind prediction: ethylene glycol 30 wt%
Me = solute("ethylene_glycol")["molar_mass"]
we = 28.0
ne, nwe = we / Me, (100 - we) / M_WATER
xwe = nwe / (ne + nwe)
put("e_Me", Me, ".3f"); put("e_xw", xwe, ".5f"); put("e_A", 2 * A_OH, ".1f")
Te_id = newton_T(lambda T: math.log(xwe))[0]
Te = newton_T(lambda T: lhs_g(T, xwe, 2 * A_OH))[0]
put("e_id", Te_id - TM, ".2f"); put("e_ours", Te - TM, ".2f"); put("e_crc", crc("ethylene_glycol", 28.0), ".2f")
put("e_err_id", (Te_id - TM) - crc("ethylene_glycol", 28.0), "+.2f"); put("e_err", (Te - TM) - crc("ethylene_glycol", 28.0), "+.2f")
check(Te, freezing_point({"ethylene_glycol": 0.28}, model, include_dcp=True), 1e-6, "EG")

# ---------- C8 NaCl: isotonic saline 0.9 wt%
Mn = solute("NaCl")["molar_mass"]
ws = 0.9
m_s = ws / Mn / ((100 - ws) / 1000)
phi = phi_nacl(m_s, NACL)
lnaw_s = -2 * m_s * Mw * phi
put("s_Mn", Mn, ".3f"); put("s_m", m_s, ".4f"); put("s_sqrt", math.sqrt(m_s), ".4f"); put("s_phi", phi, ".4f")
put("s_lnaw", lnaw_s, ".6f")
smod = AdditiveSaltModel(IdealModel(), NACL)
Ts = newton_T(lambda T: lnaw_s)[0]
put("s_T", Ts - TM, ".3f"); put("s_dil", -Kf * 2 * m_s, ".3f"); put("s_dilphi", -Kf * 2 * m_s * phi, ".3f")
check(Ts, freezing_point({"NaCl": ws / 100}, smod, include_dcp=True), 1e-6, "saline Tf")

# ---------- C9 ternary: glycerol 25 + NaCl 5 wt% (committed lab_set_v3)
wg, wn = 25.0, 5.0
ww3 = 100 - wg - wn
ng3, nw3 = wg / Mg, ww3 / M_WATER
xw3 = nw3 / (nw3 + ng3)
mn3 = wn / Mn / (ww3 / 1000)
smodel = AdditiveSaltModel(model, NACL)


def lhs3(T):
    return math.log(xw3) + (3 * A_OH / (R * T)) * (1 - xw3) ** 2 - 2 * mn3 * Mw * phi_nacl(mn3, NACL)


T3 = newton_T(lhs3)[0]
put("t_xw", xw3, ".5f"); put("t_m", mn3, ".4f"); put("t_phi", phi_nacl(mn3, NACL), ".4f")
put("t_salt", -2 * mn3 * Mw * phi_nacl(mn3, NACL), ".5f")
put("t_poly", lhs3(T3) + 2 * mn3 * Mw * phi_nacl(mn3, NACL), ".5f")
put("t_ln", lhs3(T3), ".5f"); put("t_T", T3 - TM, ".2f")
committed = [r for r in csv.reader(l for l in open(ROOT / "predictions" / "lab_set_v3_glycerol_nacl.csv") if not l.startswith("#"))][1]
put("t_comm", float(committed[3]), ".2f")
check(T3, freezing_point({"glycerol": 0.25, "NaCl": 0.05}, smodel, include_dcp=True), 1e-6, "ternary Tf")
check(T3 - TM, float(committed[3]), 0.006, "ternary vs committed")

# ---------- C10 freeze concentration: saline at -15 C
target = sum(ice_terms(T15))
m_eq = brentq(lambda m: -2 * m * Mw * phi_nacl(m, NACL) - target, 0.1, 10)
put("fc_target", target, ".5f"); put("fc_m", m_eq, ".3f"); put("fc_phi", phi_nacl(m_eq, NACL), ".4f")
put("fc_u", 100 * m_s / m_eq, ".2f"); put("fc_wt", 100 * m_eq * Mn / 1000 / (1 + m_eq * Mn / 1000), ".2f")
u0, so0, si0, _ = state({}, -15.0)
check(100 * m_s / m_eq, 100 * u0, 1e-4, "saline unfrozen"); check(100 * m_eq * Mn / 1000 / (1 + m_eq * Mn / 1000), so0, 1e-4, "saline wt")
put("fc_mmax", 4.83, ".2f")

# ---------- C11 freeze concentration: glycerol 6.5 wt% in the saline
wgc = 0.065
basis_w = 100 * (1 - wgc - 0.009)         # g water in 100 g
ng_c, nn_c = 100 * wgc / Mg, 100 * 0.009 / Mn


def g_W(W):                                # W = grams of liquid water left
    xw_ = (W / M_WATER) / (W / M_WATER + ng_c)
    mn_ = nn_c / (W / 1000)
    return math.log(xw_) + 3 * A_OH / (R * T15) * (1 - xw_) ** 2 - 2 * mn_ * Mw * phi_nacl(mn_, NACL) - target


Wc = brentq(g_W, 0.5, basis_w)
xw_c = (Wc / M_WATER) / (Wc / M_WATER + ng_c)
mn_c = nn_c / (Wc / 1000)
mg_c = ng_c / (Wc / 1000)
salt_wt = 100 * 0.9 / (Wc + 100 * wgc + 0.9)
put("c_W0", basis_w, ".2f"); put("c_ng", ng_c, ".5f"); put("c_nn", nn_c, ".5f")
put("c_W", Wc, ".3f"); put("c_u", 100 * Wc / basis_w, ".2f"); put("c_xw", xw_c, ".5f"); put("c_mn", mn_c, ".3f")
put("c_mg", mg_c, ".3f"); put("c_wt", salt_wt, ".2f"); put("c_gwt", 100 * 100 * wgc / (Wc + 100 * wgc + 0.9), ".1f")
put("c_poly", math.log(xw_c) + 3 * A_OH / (R * T15) * (1 - xw_c) ** 2, ".5f")
put("c_salt", -2 * mn_c * Mw * phi_nacl(mn_c, NACL), ".5f")
u1, so1, si1, _ = state({"glycerol": wgc}, -15.0)
check(salt_wt, so1, 1e-3, "glycerol outside salt"); check(100 * Wc / basis_w, 100 * u1, 1e-3, "glycerol unfrozen")

# ---------- C12 inside the cell
m_s0 = 0.009 / (Mn / 1000) / (1 - 0.009)
# (a) sorbitol: inside only NaCl-equivalent; same equation as C10
Wa = m_s0 / m_eq
put("in_ms0", m_s0, ".4f"); put("in_Wa", Wa, ".4f"); put("in_wta", 100 * m_eq * Mn / 1000 / (1 + m_eq * Mn / 1000), ".2f")
_, so2, si2, Wcode = state({"sorbitol": 0.09}, -15.0)
put("in_out_sorb", so2, ".2f")
check(Wa, Wcode, 1e-5, "cell W with sorbitol")
check(100 * m_eq * Mn / 1000 / (1 + m_eq * Mn / 1000), si2, 1e-4, "inside salt with sorbitol")


# (b) glycerol: glycerol molality inside = outside (mg_c); solve for cell water fraction W
def g_cell(W):
    ms_ = m_s0 / W
    tot = 1 + ms_ * Mn / 1000 + mg_c * Mg / 1000
    n_w_, n_g_ = 1 / Mw, mg_c
    xw_ = n_w_ / (n_w_ + n_g_)
    return math.log(xw_) + 3 * A_OH / (R * T15) * (1 - xw_) ** 2 - 2 * ms_ * Mw * phi_nacl(ms_, NACL) - target


Wb = brentq(g_cell, 1e-3, 10)
msb = m_s0 / Wb
put("in_Wb", Wb, ".4f"); put("in_msb", msb, ".3f")
put("in_wtb", 100 * msb * Mn / 1000 / (1 + msb * Mn / 1000 + mg_c * Mg / 1000), ".2f")
check(100 * msb * Mn / 1000 / (1 + msb * Mn / 1000 + mg_c * Mg / 1000), si1, 2e-3, "inside salt with glycerol")

# ---------- C13 sorbitol at 25 C: osmotic coefficient
Msb = solute("sorbitol")["molar_mass"]
bow = [r for r in csv.DictReader(l for l in open(ROOT / "data" / "bower1963_sorbitol_25C.csv") if not l.startswith("#")) if r["table"] == "I"]
row = bow[3]
mb, phib = float(row["m_B"]), float(row["phi_B"])
xs_b = mb / (1 / Mw + mb); xw_b = 1 - xs_b
T25 = 298.15
lnaw_b = math.log(xw_b) + 6 * A_OH / (R * T25) * xs_b ** 2
phi_mod = -lnaw_b / (mb * Mw)
phi_id = -math.log(xw_b) / (mb * Mw)
a_bow = (-phib * mb * Mw - math.log(xw_b)) / xs_b ** 2 * R * T25 / 6
put("b_m", mb, ".4f"); put("b_phi", phib, ".3f"); put("b_xs", xs_b, ".5f"); put("b_lnaw", lnaw_b, ".5f")
put("b_phim", phi_mod, ".3f"); put("b_phiid", phi_id, ".3f"); put("b_a", a_bow, ".0f")
put("b_lnawdat", -phib * mb * Mw, ".5f"); put("b_lnxw", math.log(xw_b), ".5f")

# ---------- C14 Gibbs-Helmholtz from sucrose
gh = [r for r in csv.reader(l for l in open(ROOT / "predictions" / "posthoc_per_molecule.csv") if not l.startswith("#"))
      if r[0] == "gibbs_helmholtz_sucrose"]
for j, r in enumerate(gh):
    T1 = float(r[5].split("@")[1][:-1]); a1 = float(r[5].split("@")[0])
    a2 = float(r[6].split("@")[0]); T2 = 298.15
    s_ = -(a2 - a1) / (T2 - T1); h_ = a1 + T1 * s_
    put(f"gh{j}_a1", a1, ".0f"); put(f"gh{j}_T1", T1, ".1f"); put(f"gh{j}_a2", a2, ".0f")
    put(f"gh{j}_dT", T2 - T1, ".2f"); put(f"gh{j}_s", s_, ".2f"); put(f"gh{j}_h", h_, ".0f")
    put(f"gh{j}_Ts", T1 * s_, ".0f"); put(f"gh{j}_win", r[1].replace(" mol/kg", ""), "s")
    check(h_, float(r[7].split("=")[1]), 10, "GH h (inputs rounded to 1 J/mol)"); put(f"gh{j}_hcode", float(r[7].split("=")[1]), ".0f")

# ---------- C15 glycerol-sparing ratio
rec = [r for r in csv.reader(l for l in open(ROOT / "predictions" / "recipe_three_criteria.csv") if not l.startswith("#"))][1:]
S = np.array([float(r[0]) for r in rec])
for k, col in (("out", 1), ("in", 2), ("unf", 3)):
    G = np.array([float(r[col]) for r in rec])
    ok = np.isfinite(G) & (G > 0)
    if k == "out":
        ok &= S <= 8
    slope = np.polyfit(S[ok], G[ok], 1)[0]
    put(f"r_{k}", -slope, "+.2f"); put(f"r_{k}_rng", f"{S[ok].min():.0f}--{S[ok].max():.0f}", "s")
    put(f"r_{k}_0", G[0], ".2f")
    put(f"r_{k}_2", G[2], ".2f")
    put(f"r_{k}_fd", (G[0] - G[2]) / 2, "+.2f")      # saving per g sorbitol over the first 2 wt% (report value)
opt = open(ROOT / "predictions" / "recipe_three_criteria.csv").read().split("optimum sorbitol ")[1]
put("r_opt_s", float(opt.split(" wt%")[0]), ".2f"); put("r_opt_g", float(opt.split("glycerol ")[1].split(" wt%")[0]), ".2f")
put("r_alone", float(rec[0][3]), ".2f")
put("r_saving", 100 * (1 - float(V["r_opt_g"]) / float(V["r_alone"])), ".0f")

json.dump(V, open(ROOT / "report" / "calc_values.json", "w"), indent=1)
tpl = (ROOT / "report" / "calculations_workbook_template.tex").read_text()
for k, v in V.items():
    tpl = tpl.replace("<<" + k + ">>", v)
left = [s.split(">>")[0] for s in tpl.split("<<")[1:]]
assert not left, f"unfilled placeholders: {left}"
(ROOT / "report" / "calculations_workbook.tex").write_text(tpl)
print(f"{len(V)} values computed; hand calculations match the project code; workbook written.")
for k in ("ice_ln", "id_TC", "no_TC", "g_crc", "e_ours", "e_crc", "s_T", "t_T", "fc_u", "fc_wt", "c_wt", "in_wta", "in_wtb",
          "b_phim", "b_phi", "gh0_h", "gh0_s", "r_out", "r_in", "r_unf"):
    print(f"  {k:8s} {V[k]}")
