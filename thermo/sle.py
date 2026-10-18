"""Solid-liquid equilibrium between pure ice and an aqueous solution.

Equilibrium condition (fugacity of water equal in ice and solution), with the
fusion free energy obtained from the heat-melt-cool construction:

    ln(x_w gamma_w) = -(dH/R)(1/T - 1/T_m)                        [constant dH]
                      - (dCp/R)(T_m/T - 1 + ln(T/T_m))            [optional dCp term]

The right-hand side depends only on T. So once ice is present, the water
activity of the liquid in contact with it is fixed by temperature alone,
whatever the solutes are.
"""

import math
import numpy as np
from scipy.optimize import brentq

from .constants import R, T_M_WATER, DH_FUS_WATER, DCP_FUS_WATER, V_WATER
from .composition import mole_fractions, validate


def ln_aw_ice(T, include_dcp=False):
    """ln(water activity) of a solution in equilibrium with ice at T."""
    val = -(DH_FUS_WATER / R) * (1.0 / T - 1.0 / T_M_WATER)
    if include_dcp:
        val -= (DCP_FUS_WATER / R) * (T_M_WATER / T - 1.0 + math.log(T / T_M_WATER))
    return val


def ln_aw_solution(massfrac, T, model):
    if hasattr(model, "ln_aw"):                 # composite models (e.g. salt.AdditiveSaltModel)
        return model.ln_aw(massfrac, T)
    x_w, xs = mole_fractions(massfrac)
    return math.log(x_w) + model.ln_gamma_water(T, x_w, xs)


def freezing_point(massfrac, model, include_dcp=False, T_low=150.0):
    """Liquidus temperature (K) of a solution of given composition."""
    validate(massfrac)
    if sum(massfrac.values()) == 0:
        return T_M_WATER

    def f(T):
        return ln_aw_solution(massfrac, T, model) - ln_aw_ice(T, include_dcp)

    T_high = T_M_WATER - 1e-9
    if f(T_high) >= 0:
        raise ValueError("No freezing point found below T_m for this composition/model")
    if f(T_low) <= 0:
        raise ValueError(f"Liquidus lies below {T_low} K or the model has left its "
                         "valid range for this composition")
    return brentq(f, T_low, T_high, xtol=1e-9)


def residual_liquid(massfrac0, T, model, include_dcp=False):
    """Freeze-concentration at temperature T for a starting composition.

    Ice is pure water, so solute stays in the shrinking liquid. Returns a dict:
      liquid_fraction  : mass of liquid / starting mass
      ice_fraction     : mass of ice / starting mass
      composition      : mass fractions of each solute in the residual liquid
      total_solute_wt  : total solute wt% in the residual liquid
    """
    Tf0 = freezing_point(massfrac0, model, include_dcp)
    if T >= Tf0:
        comp = dict(massfrac0)
        return {"liquid_fraction": 1.0, "ice_fraction": 0.0, "composition": comp,
                "total_solute_wt": 100 * sum(comp.values()), "T_freeze0": Tf0}

    m_solutes = dict(massfrac0)
    m_water0 = 1.0 - sum(m_solutes.values())

    def liquid_comp(phi):                       # phi = fraction of water frozen
        m_w = m_water0 * (1.0 - phi)
        total = m_w + sum(m_solutes.values())
        return {k: v / total for k, v in m_solutes.items()}, total

    def g(phi):
        comp, _ = liquid_comp(phi)
        return freezing_point(comp, model, include_dcp) - T

    # bracket by stepping up the frozen fraction until the liquid's freezing
    # point drops below T; stop if the model leaves its valid range first
    phi_lo, phi_hi = 0.0, None
    for phi_try in np.linspace(0.02, 0.99, 98):
        try:
            if g(phi_try) < 0:
                phi_hi = phi_try
                break
            phi_lo = phi_try
        except ValueError:
            break
    if phi_hi is None:
        raise ValueError("T is outside the range where the ice liquidus applies "
                         "(below the eutectic or glass region for this composition)")
    phi = brentq(g, phi_lo, phi_hi, xtol=1e-10)
    comp, liq = liquid_comp(phi)
    return {"liquid_fraction": liq, "ice_fraction": 1.0 - liq, "composition": comp,
            "total_solute_wt": 100 * sum(comp.values()), "T_freeze0": Tf0}


def osmotic_pressure_at_ice(T, include_dcp=False):
    """Osmotic pressure (Pa) of any liquid in equilibrium with ice at T.

    Independent of which solutes are present: equilibrium with ice fixes a_w.
    """
    return -R * T * ln_aw_ice(T, include_dcp) / V_WATER


def osmotic_pressure(massfrac, T, model):
    """Osmotic pressure (Pa) of an unfrozen solution at T."""
    return -R * T * ln_aw_solution(massfrac, T, model) / V_WATER
