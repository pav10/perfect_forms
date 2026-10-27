"""Exact short-vector enumeration for positive definite rational quadratic forms.

The Gram matrix is scaled to integers and LLL-reduced (PARI's qflllgram; the
transformation is checked to be unimodular, so LLL only affects speed).  The
enumeration is Fincke--Pohst on the exact rational Gauss decomposition

    q(y) = sum_i d_i (y_i + sum_{j>i} mu_ij y_j)^2,

rewritten with integers only: with delta_i the common denominator of row i,
t_i = delta_i y_i + sum_{j>i} delta_i mu_ij y_j is an integer and the term is
W_i t_i^2 / L for integers W_i, L.  Each loop range is then an exact integer
interval (integer square root), so no vector is ever missed and no floating point
is used at all.  This generalises the 3-variable enumerator ``vu_enum`` of
voronoi_unary.gp to any dimension.
"""
from fractions import Fraction
from math import gcd, isqrt

from .pari import pari


def _lcm(a, b):
    return a * b // gcd(a, b)


def gauss_reduction(G):
    """Exact rational Gauss reduction (PARI's qfgaussred convention): returns (d, mu)
    with q(y) = sum_i d[i] (y_i + sum_{j>i} mu[i][j] y_j)^2, or None if G is not
    positive definite."""
    n = len(G)
    A = [[Fraction(x) for x in row] for row in G]
    d = [Fraction(0)] * n
    mu = [[Fraction(0)] * n for _ in range(n)]
    for i in range(n):
        di = A[i][i]
        if di <= 0:
            return None
        d[i] = di
        for j in range(i + 1, n):
            mu[i][j] = A[i][j] / di
        for j in range(i + 1, n):
            if A[i][j]:
                f = A[i][j]
                for k in range(j, n):
                    A[j][k] -= f * mu[i][k]
                    if k != j:
                        A[k][j] = A[j][k]
    return d, mu


def is_positive_definite(G):
    return gauss_reduction(G) is not None


def integral_gram(G):
    """(den, Gi) with Gi = den * G an integer matrix, den > 0 minimal."""
    den = 1
    for row in G:
        for x in row:
            den = _lcm(den, Fraction(x).denominator)
    return den, [[int(Fraction(x) * den) for x in row] for row in G]


def lll(Gi):
    """(U, U^T Gi U) with U unimodular (columns = new basis), Gi integral and PD."""
    n = len(Gi)
    Up = pari.qflllgram(pari.matrix(n, n, [x for row in Gi for x in row]))
    if int(pari.matsize(Up)[1]) != n:
        raise ArithmeticError('qflllgram: matrix not positive definite')
    U = [[int(Up[i, j]) for j in range(n)] for i in range(n)]
    GU = [[sum(Gi[i][k] * U[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    Gr = [[sum(U[k][i] * GU[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    return U, Gr


class _Enumerator:
    """Integer data of the Gauss reduction of an integral PD matrix."""

    def __init__(self, Gi):
        red = gauss_reduction(Gi)
        if red is None:
            raise ArithmeticError('form is not positive definite')
        d, mu = red
        n = self.n = len(Gi)
        self.delta = []
        self.a = []
        w = []
        for i in range(n):
            de = 1
            for j in range(i + 1, n):
                de = _lcm(de, mu[i][j].denominator)
            self.delta.append(de)
            self.a.append([0] * (i + 1) + [int(mu[i][j] * de) for j in range(i + 1, n)])
            w.append(d[i] / (de * de))
        L = 1
        for x in w:
            L = _lcm(L, x.denominator)
        self.L = L
        self.W = [int(x * L) for x in w]

    def run(self, bound):
        """All y != 0, one per +- pair (last nonzero coordinate > 0), with
        sum_i W_i t_i^2 <= bound (an integer).  Returns list of (value, y)."""
        n, W, a, delta = self.n, self.W, self.a, self.delta
        out = []
        y = [0] * n

        def rec(i, budget, sym):
            ai = a[i]
            s = 0
            for j in range(i + 1, n):
                if y[j]:
                    s += ai[j] * y[j]
            Wi, de = W[i], delta[i]
            T = isqrt(budget // Wi)
            lo = -((T + s) // de)          # ceil((-T - s) / de)
            hi = (T - s) // de             # floor((T - s) / de)
            if sym and lo < 0:
                lo = 0
            for v in range(lo, hi + 1):
                t = de * v + s
                rem = budget - Wi * t * t
                if rem < 0:
                    continue
                y[i] = v
                if i == 0:
                    if not (sym and v == 0):
                        out.append((bound - rem, tuple(y)))
                else:
                    rec(i - 1, rem, sym and v == 0)
            y[i] = 0

        rec(n - 1, bound, True)
        return out


def _normalize_sign(x):
    for t in x:
        if t:
            return x if t > 0 else tuple(-u for u in x)
    return x


def short_vectors(G, B, strict=False):
    """All x != 0 (one per +- pair, first nonzero coordinate positive) with
    x^T G x <= B (< B if strict), for a rational PD matrix G.  Returns a list of
    (value, x) with exact Fraction values, sorted by value then x."""
    den, Gi = integral_gram(G)
    U, Gr = lll(Gi)
    E = _Enumerator(Gr)
    Bs = Fraction(B) * den * E.L
    bound = Bs.numerator // Bs.denominator
    n = len(G)
    out = []
    for val, yv in E.run(bound):
        q = Fraction(val, den * E.L)
        if q > B or (strict and q == B):
            continue
        x = tuple(sum(U[i][j] * yv[j] for j in range(n)) for i in range(n))
        out.append((q, _normalize_sign(x)))
    out.sort()
    return out


def minimum(G):
    """(min, minimal vectors) of a rational PD matrix, exactly; one vector per +- pair."""
    den, Gi = integral_gram(G)
    U, Gr = lll(Gi)
    E = _Enumerator(Gr)
    n = len(G)
    B = min(Gr[i][i] for i in range(n)) * E.L
    vecs = E.run(B)
    m = min(v for v, _ in vecs)
    out = []
    for val, yv in vecs:
        if val == m:
            x = tuple(sum(U[i][j] * yv[j] for j in range(n)) for i in range(n))
            out.append(_normalize_sign(x))
    out.sort()
    return Fraction(m, den * E.L), out


def qeval(G, x):
    n = len(x)
    return sum(G[i][j] * x[i] * x[j] for i in range(n) for j in range(n))
