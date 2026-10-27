"""Eutaxy, extremality and the additive Hermite--Humbert constant gamma_K.

For a unary form alpha (totally positive, minimum m(alpha)) put
gamma(alpha) = m(alpha) / N(alpha)^(1/m), the invariant of Leibak's additive
theory; gamma_K = sup gamma.  The gradient of log N at alpha is the functional
alpha' -> Tr(alpha' alpha^{-1}), and the gradients of the minimal vectors are
alpha' -> Tr(alpha' x^2).  So alpha is

* eutactic if alpha^{-1} = sum_x rho_x x^2 with all rho_x > 0, i.e. alpha^{-1}
  lies in the interior of the Voronoi cone cone{x^2 : x in S(alpha)};
* weakly eutactic if it lies in the closed cone.

For a perfect form the cone is full dimensional, with facet normals beta_F
(Tr(beta_F x^2) >= 0 on S(alpha), = 0 on F), so eutaxy is the exact test
Tr(beta_F alpha^{-1}) > 0 for every facet F.  Voronoi's theorem in this setting
(extreme <=> perfect and eutactic, the argument of Bavard for a minimum of linear
forms divided by a strictly log-concave homogeneous function) is recorded in
docs/THEORY.md, section 4, which also records what Leibak's thesis (Thm 2.1, checked
against the text) does and does not cover: its eutaxy statement is the weak one, vacuous
for unary forms.

gamma is maximal at a perfect form, so gamma_K = max over the perfect classes of
N(alpha)^(-1/m) (minimum 1): gamma_K^m = 1 / min N(alpha) is rational and computed
exactly.  The maximiser must be eutactic (a consistency check).
"""
from fractions import Fraction

# gamma_m^m for the lattice Hermite constants known exactly (m <= 8)
HERMITE_POW = {1: Fraction(1), 2: Fraction(4, 3), 3: Fraction(2), 4: Fraction(4), 5: Fraction(8),
               6: Fraction(64, 3), 7: Fraction(64), 8: Fraction(256)}


def eutaxy(O, P):
    """'eutactic', 'weakly' (on the boundary of the cone) or 'no', for a PerfectForm P."""
    inv = O.inverse(P.alpha)
    vals = [O.trace_pair(beta, inv) for beta, _ in P.facets]
    if all(v > 0 for v in vals):
        return 'eutactic'
    if all(v >= 0 for v in vals):
        return 'weakly'
    return 'no'


def gamma_pow(P):
    """gamma(alpha)^m = 1 / N(alpha) for a perfect form of minimum 1 (exact)."""
    return 1 / Fraction(P.norm)


def summary(O, V):
    """Exact gamma_K^m, the maximising classes and the eutaxy of every class."""
    eu = [eutaxy(O, P) for P in V.classes]
    g = [gamma_pow(P) for P in V.classes]
    gmax = max(g)
    argmax = [i for i, v in enumerate(g) if v == gmax]
    out = {
        'eutaxy': eu,
        'n_extreme': sum(e == 'eutactic' for e in eu),
        'n_weakly_eutactic': sum(e == 'weakly' for e in eu),
        'gamma_pow': gmax,                  # gamma_K^m
        'gamma': float(gmax) ** (1.0 / O.m),
        'argmax': argmax,
        'argmax_eutactic': all(eu[i] == 'eutactic' for i in argmax),
    }
    hp = HERMITE_POW.get(O.m)
    if hp is not None:
        # Leibak's bound gamma_K <= gamma_m D^(1/m), i.e. gamma_K^m <= gamma_m^m D
        out['leibak_ratio'] = gmax / (hp * O.disc)     # <= 1; = 1 iff the bound is attained
    return out
