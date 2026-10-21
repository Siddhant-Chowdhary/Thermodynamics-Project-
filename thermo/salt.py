"""Non-ideal NaCl and the additive salt + polyol water-activity model.

NaCl is described by an empirical osmotic coefficient fitted to the NaCl-water
ice liquidus (data/literature.csv):

    phi(m) = 1 - a*sqrt(m)/(1 + sqrt(m)) + b*m + c*m^2      (m in mol/kg water)
    ln a_w(NaCl part) = -nu * m * M_w * phi(m),   nu = 2

The sqrt term has the Debye-Hueckel shape (phi < 1 when dilute); b, c capture
the rise at high m. Each data point sits at its own liquidus T, so phi carries
the temperature along the liquidus (assumption: phi(m) is used at nearby T).

Mixtures (assumption agreed 18 Oct): contributions are ADDITIVE at equal
molality, i.e. no salt-polyol cross term:

    ln a_w(mix) = ln a_w(polyols alone in the same water) + ln a_w(NaCl alone in the same water)

Kleinhans & Mazur (Cryobiology 54, 212, 2007) found this additivity reproduces
glycerol + NaCl ternary freezing data well up to ~55 wt%. We TEST it against
Shepard et al. (1976) rather than fitting to it.
"""
import csv, math
from pathlib import Path
import numpy as np

from .constants import M_WATER, solute
from .sle import ln_aw_ice

M_MAX_FIT = 4.83  # mol/kg, top of the CRC 2005 NaCl data used in the fit (22 wt%)
_DATA = Path(__file__).resolve().parents[1] / "data" / "literature.csv"


def phi_nacl(m, params):
    a, b, c = params
    s = math.sqrt(m)
    return 1.0 - a * s / (1.0 + s) + b * m + c * m * m


def fit_nacl(path=_DATA, include_dcp=True, max_wt=23.3):
    """Fit (a, b, c) by linear least squares to the NaCl ice-liquidus rows."""
    rows = [r for r in csv.DictReader(open(path))
            if r["solute"] == "NaCl" and r["phase_branch"] == "ice" and r["source_id"] == "CRC2005"
            and float(r["wt_percent"]) <= max_wt]
    M = solute("NaCl")["molar_mass"]
    X, y, pts = [], [], []
    for r in rows:
        w = float(r["wt_percent"]); T = float(r["T_freeze_C"]) + 273.15
        m = w / (100 - w) * 1000 / M
        phi = -ln_aw_ice(T, include_dcp) / (2 * m * M_WATER / 1000)
        s = math.sqrt(m)
        X.append([-s / (1 + s), m, m * m]); y.append(phi - 1.0); pts.append((m, phi, T))
    params, *_ = np.linalg.lstsq(np.array(X), np.array(y), rcond=None)
    return tuple(float(p) for p in params), pts


class AdditiveSaltModel:
    """Wraps a polyol activity model; adds non-ideal NaCl additively."""

    def __init__(self, polyol_model, nacl_params=None, include_dcp=True):
        self.polyol = polyol_model
        self.params = nacl_params or fit_nacl(include_dcp=include_dcp)[0]
        self.name = f"{polyol_model.name} + non-ideal NaCl (additive)"

    # used by thermo.sle.ln_aw_solution when present
    def ln_aw(self, massfrac, T):
        w_water = 1.0 - sum(massfrac.values())
        n_w = w_water / M_WATER
        polys = {k: v / solute(k)["molar_mass"] for k, v in massfrac.items() if k != "NaCl" and v > 0}
        ln_aw = 0.0
        if polys:
            tot = n_w + sum(polys.values())
            x_w = n_w / tot
            xs = {k: v / tot for k, v in polys.items()}
            ln_aw += math.log(x_w) + self.polyol.ln_gamma_water(T, x_w, xs)
        w_salt = massfrac.get("NaCl", 0.0)
        if w_salt > 0:
            m = (w_salt / solute("NaCl")["molar_mass"]) / (w_water / 1000.0)   # mol/kg water
            if m > M_MAX_FIT and not getattr(self, "_warned", False):
                import warnings
                warnings.warn(f"NaCl molality {m:.2f} mol/kg is above the fitted range ({M_MAX_FIT} mol/kg)")
                self._warned = True
            ln_aw += -2 * m * (M_WATER / 1000.0) * phi_nacl(m, self.params)
        return ln_aw
