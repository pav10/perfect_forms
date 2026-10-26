"""The exact neighbour walk of Voronoi's algorithm, for any linear space of forms.

A space of forms is given by linear maps  coordinates -> rational Gram matrix  on a
fixed Z-lattice.  Along alpha + t beta (alpha of minimum 1, beta vanishing on the
minimal vectors of a face F and positive on the others) the neighbour is at
rho = min over v with q_beta(v) < 0 of (q_alpha(v) - 1)/(-q_beta(v)).  It is found by
the classical iteration: for any t with alpha + t beta positive definite and minimum
< 1, the vectors below 1 give an upper bound for rho; iterate until the minimum is
exactly 1.  Positive definiteness is decided by exact Gauss reduction; the minimum
by the integer Fincke--Pohst of qfmin.py.  Floating point only estimates the
boundary of positivity, to choose the starting t.
"""
from fractions import Fraction

import numpy as np

from . import qfmin


class WalkError(ArithmeticError):
    pass


def rational_between(lo, hi):
    """A rational strictly between lo < hi, of small height."""
    lo, hi = Fraction(lo), Fraction(hi)
    mid = (lo + hi) / 2
    B = 16
    while True:
        r = mid.limit_denominator(B)
        if lo < r < hi:
            return r
        B *= 16


def positivity_bound(Ga, Gb):
    """Float estimate of sup{t > 0 : Ga + t Gb positive definite} (Ga PD), or None if
    the ray stays positive definite.  Ga is LLL-reduced first (an exact unimodular
    change of basis applied to both matrices) to keep the float problem well
    conditioned."""
    den, Gi = qfmin.integral_gram(Ga)
    U, Gr = qfmin.lll(Gi)
    n = len(Ga)
    Bu = [[sum(U[k][i] * Gb[k][l] * U[l][j] for k in range(n) for l in range(n)) for j in range(n)]
          for i in range(n)]
    A = np.array([[float(x) for x in r] for r in Gr])
    B = np.array([[float(x * den) for x in r] for r in Bu])
    d = np.sqrt(np.diag(A))
    A = A / np.outer(d, d)
    B = B / np.outer(d, d)
    L = np.linalg.cholesky(A)
    Li = np.linalg.inv(L)
    M = Li @ B @ Li.T
    lam = np.linalg.eigvalsh((M + M.T) / 2).min()
    if lam >= 0:
        return None
    return -1.0 / lam


def walk(Ga, Gb, qa, qb, stats=None, tmax=None):
    """Exact rho for the ray Ga + t Gb (Gram matrices, Fractions); qa(x), qb(x) give the
    exact values of the two forms at a lattice vector x.  tmax: an estimate of the
    positivity boundary (computed from the Gram matrices if not given).  Returns rho."""
    n = len(Ga)
    if tmax is None:
        tmax = positivity_bound(Ga, Gb)
    if tmax is None:
        raise WalkError('walk: the ray stays positive definite')
    if not tmax > 0:
        raise WalkError('walk: bad estimate of the positivity boundary')

    def G_at(t):
        return [[Ga[i][j] + t * Gb[i][j] for j in range(n)] for i in range(n)]

    def below(G):
        if stats is not None:
            stats['enumerations'] = stats.get('enumerations', 0) + 1
        return qfmin.short_vectors(G, 1, strict=True)

    lo = Fraction(0)
    t = rational_between(0, tmax)
    steps = 0
    while True:
        G = G_at(t)
        if not qfmin.is_positive_definite(G):
            t = rational_between(lo, t)
            continue
        S = below(G)
        if S:
            break
        lo = t
        steps += 1
        # tmax is only an estimate: approach it, and after many steps also allow t
        # beyond it (positivity is re-checked exactly)
        t = rational_between(t, tmax) if t < tmax and steps < 60 else 2 * t
    while S:
        tn = t
        for _, x in S:
            a, b = qa(x), qb(x)
            if b >= 0:
                raise WalkError('walk: vector below 1 with q_beta(x) >= 0')
            tn = min(tn, (a - 1) / (-b))
        t = tn
        S = below(G_at(t))
    return t
