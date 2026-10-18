"""Conversions between mass fractions and mole fractions.

A composition is a dict {solute_name: mass_fraction}; water makes up the rest.
Example: {"glycerol": 0.10} is 10 wt% glycerol in water.
Dissociating solutes (NaCl) contribute nu moles of particles per formula unit.
"""

from .constants import M_WATER, solute


def validate(massfrac):
    total = sum(massfrac.values())
    if any(w < 0 for w in massfrac.values()):
        raise ValueError("Mass fractions must be non-negative")
    if total >= 1:
        raise ValueError(f"Solute mass fractions sum to {total:.3f}; must be below 1")
    return 1.0 - total


def moles_per_gram(massfrac):
    """Moles of water and of each solute's particles per 1 g of solution."""
    w_water = validate(massfrac)
    n_w = w_water / M_WATER
    n = {name: w * solute(name)["nu"] / solute(name)["molar_mass"]
         for name, w in massfrac.items() if w > 0}
    return n_w, n


def mole_fractions(massfrac):
    """Return (x_water, {solute: x_solute}) on a particle basis."""
    n_w, n = moles_per_gram(massfrac)
    total = n_w + sum(n.values())
    return n_w / total, {k: v / total for k, v in n.items()}


def from_wt_percent(**wt_percent):
    """Convenience: from_wt_percent(glycerol=10) -> {'glycerol': 0.10}."""
    return {k: v / 100.0 for k, v in wt_percent.items()}
