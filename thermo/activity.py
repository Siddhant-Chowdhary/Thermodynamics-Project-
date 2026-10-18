"""Activity coefficient models.

All models are written as an excess Gibbs energy containing only water-solute
pair interactions (solute-solute interactions neglected):

    G_ex / RT = sum_i  alpha_i * x_i * x_w,       alpha_i = A_i / (R T)

Differentiating n*G_ex/RT gives

    ln gamma_w = (1 - x_w) * sum_i alpha_i x_i
    ln gamma_k = alpha_k * x_w - x_w * sum_i alpha_i x_i

For a single solute this reduces exactly to the one-constant Margules model
(RT ln gamma_w = A x_s^2). The models differ only in how A_i is set:

    IdealModel               A_i = 0
    MargulesModel            A_i given per solute (fitted molecule by molecule)
    GroupContributionModel   A_i = n_OH,i * a_OH  (one parameter shared by all polyols)

The group-contribution model is what allows parameters fitted on one polyol to
predict another, and to predict mixtures, with no new parameters.
"""

from .constants import R, solute


class _PairModel:
    name = "pair model"

    def A(self, name):
        raise NotImplementedError

    def alphas(self, T, xs):
        return {k: self.A(k) / (R * T) for k in xs}

    def ln_gamma_water(self, T, x_w, xs):
        a = self.alphas(T, xs)
        return (1.0 - x_w) * sum(a[k] * xs[k] for k in xs)

    def ln_gamma_solutes(self, T, x_w, xs):
        a = self.alphas(T, xs)
        s = sum(a[k] * xs[k] for k in xs)
        return {k: a[k] * x_w - x_w * s for k in xs}


class IdealModel(_PairModel):
    name = "ideal"

    def A(self, name):
        return 0.0


class MargulesModel(_PairModel):
    """One Margules constant (J/mol) per solute, e.g. {'glycerol': -8000}."""
    name = "Margules (per molecule)"

    def __init__(self, A_by_solute):
        self.A_by_solute = dict(A_by_solute)

    def A(self, name):
        return self.A_by_solute.get(name, 0.0)


class GroupContributionModel(_PairModel):
    """A_i = n_OH,i * a_OH. Ions (n_OH = 0) behave ideally."""
    name = "group contribution (per -OH)"

    def __init__(self, a_OH):
        self.a_OH = float(a_OH)

    def A(self, name):
        return solute(name)["n_OH"] * self.a_OH
