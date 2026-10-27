"""Perfect unary forms over number fields of signature (1, 1) (complex cubic fields), plan E2.

Here d = r1 + r2 = 2: a form is a = (a_1, a_2) in R^2_{>0} with value
a_1 sigma_1(x)^2 + 2 a_2 |sigma_2(x)|^2 on x in O_K (manuscript eq. (2.1), e_sigma = 2 at the
complex place), and dK_1(m_O)/O^x is a circle (n_K = mu_K = number of edges).

Everything lives in K through the real embedding sigma_1: |sigma_2(x)|^2 = N(x)/sigma_1(x), so
q(x) = (sigma_1(x^2), 2 sigma_1(N(x)/x)), perfect forms have a_i = sigma_1(A_i) with A_i in K, and
the value of a form at x is sigma_1(A_1 x^2 + 2 A_2 N(x)/x).  Signs are exact and cheap:
sigma_1(y) has the sign of N(y), since N(y) = sigma_1(y) |sigma_2(y)|^2.

Minimal vectors: candidates from a floating point Fincke--Pohst with a widened bound (floating
point only widens the search), then every value is compared exactly in K.  Units act by
(A_1, A_2) -> (A_1 eps^2, A_2 N(eps)/eps); equivalence is decided exactly (square root of
A_1'/A_1 in K, then the second coordinate).  Certificates: the circle has as many vertices as
edges, and the side-pairing units generate O^x/+-1 (index 1 in bnfinit's unit group).

Benchmark: Komatsu--Watanabe, IJNT 10 (2014), Thms 3.1-3.2: for v^3 + l v - 1 = 0, l >= 2,
4 l^3 + 27 squarefree, n_K = 1 with S(a_0) = {+-1, +-v}, a_0 = ((v^-1 - 1)/l, (1 - v^2)/(2l)).
"""
from fractions import Fraction

import numpy as np

from . import linalg
from .pari import pari
from .field import member, to_frac


class ComplexCubicError(ArithmeticError):
    pass


class ComplexCubic:
    def __init__(self, pol):
        self.pol = pari(pol) if isinstance(pol, str) else pol
        nf = self.nf = pari.nfinit(self.pol)
        if [int(t) for t in member(nf, 'sign')] != [1, 1]:
            raise ValueError('need signature (1, 1)')
        m = self.m = 3
        zk = member(nf, 'zk')
        self.basis_alg = [pari.nfbasistoalg(nf, pari.Col([int(i == j) for i in range(m)])) for j in range(m)]
        self.mult = [[tuple(int(t) for t in pari.nfalgtobasis(nf, pari.nfeltmul(nf, pari.Col([int(i == j) for i in range(m)]),
                                                                                  pari.Col([int(i == k) for i in range(m)]))))
                      for k in range(m)] for j in range(m)]
        self.one = (1, 0, 0)
        roots = pari.polroots(self.pol, precision=256)
        real = [r for r in roots if abs(float(pari.imag(r))) < 1e-30]
        cplx = [r for r in roots if float(pari.imag(r)) > 1e-30]
        self.r1 = real[0]
        self.r2 = cplx[0]
        self.e1 = np.array([float(pari.real(pari.subst(pari.lift(b), 'x', self.r1))) for b in self.basis_alg])
        self.e2 = np.array([complex(float(pari.real(z)), float(pari.imag(z)))
                            for z in (pari.subst(pari.lift(b), 'x', self.r2) for b in self.basis_alg)])
        self.disc = int(member(nf, 'disc'))
        self._bnf = None

    # --- exact arithmetic in K (coordinates in nf.zk) -------------------------

    def mul(self, x, y):
        out = [0, 0, 0]
        for j, xj in enumerate(x):
            if xj:
                for k, yk in enumerate(y):
                    if yk:
                        p = xj * yk
                        for l, c in enumerate(self.mult[j][k]):
                            if c:
                                out[l] += p * c
        return tuple(out)

    def mult_matrix(self, x):
        cols = [self.mul(x, tuple(int(i == k) for i in range(3))) for k in range(3)]
        return [[cols[k][i] for k in range(3)] for i in range(3)]

    def norm(self, x):
        return linalg.det(self.mult_matrix(x))

    def inverse(self, x):
        return tuple(linalg.solve(self.mult_matrix(x), self.one))

    def sign(self, y):
        """Exact sign of sigma_1(y)."""
        if not any(y):
            return 0
        n = self.norm(y)
        return 1 if n > 0 else -1

    def s(self, x):
        return self.mul(x, x)

    def t(self, x):
        """N(x)/x (so sigma_1(t(x)) = |sigma_2(x)|^2)."""
        n = self.norm(x)
        return tuple(n * c for c in self.inverse(x))

    def value(self, A, x):
        """A_1 x^2 + 2 A_2 N(x)/x in K; its real embedding is the value of the form at x."""
        A1, A2 = A
        u, v = self.mul(A1, self.s(x)), self.mul(A2, self.t(x))
        return tuple(Fraction(a) + 2 * Fraction(b) for a, b in zip(u, v))

    def to_alg(self, x):
        return pari.nfbasistoalg(self.nf, pari.Col([pari(Fraction(t).numerator) / Fraction(t).denominator for t in x]))

    def from_alg(self, a):
        return tuple(to_frac(t) for t in pari.nfalgtobasis(self.nf, a))

    def real(self, y):
        return float(sum(float(t) * e for t, e in zip(y, self.e1)))

    # --- minimal vectors ------------------------------------------------------

    def gram_float(self, a1, a2):
        """Gram matrix (floats) of x -> a1 sigma_1(x)^2 + 2 a2 |sigma_2(x)|^2 on Z^3."""
        G = a1 * np.outer(self.e1, self.e1) + 2 * a2 * np.real(np.outer(self.e2, np.conj(self.e2)))
        return (G + G.T) / 2

    def candidates(self, A, bound):
        """All x != 0 (one per +- pair) whose value is <= bound, plus possibly more (floating point
        Fincke--Pohst with a widened bound)."""
        a1, a2 = self.real(A[0]), self.real(A[1])
        G = self.gram_float(a1, a2)
        Gp = pari.matrix(3, 3, [pari(float(x)) for x in G.flatten()])
        U = pari.qflllgram(Gp)
        U = [[int(U[i, j]) for j in range(3)] for i in range(3)]
        if abs(linalg.det(U)) != 1:
            raise ComplexCubicError('LLL not unimodular')
        Um = np.array(U, dtype=float)
        Gr = Um.T @ G @ Um
        B = float(bound) * (1 + 1e-7) + 1e-12
        V = pari.qfminim(pari.matrix(3, 3, [pari(float(x)) for x in Gr.flatten()]), B, None, 2)[2]
        out = []
        for j in range(int(pari.matsize(V)[1])):
            y = [int(V[i, j]) for i in range(3)]
            x = tuple(sum(U[i][k] * y[k] for k in range(3)) for i in range(3))
            out.append(x)
        return out

    def minimum(self, A, guess=None):
        """(min as an element of K, minimal vectors), exactly among the candidates."""
        if guess is None:
            guess = min(self.real(self.value(A, x)) for x in [(1, 0, 0), (0, 1, 0), (0, 0, 1)])
        C = self.candidates(A, guess)
        best, vecs = None, []
        for x in C:
            v = self.value(A, x)
            if best is None:
                best, vecs = v, [x]
                continue
            d = self.sign(tuple(p - q for p, q in zip(v, best)))
            if d < 0:
                best, vecs = v, [x]
            elif d == 0:
                vecs.append(x)
        return best, [_norm_sign(x) for x in vecs]


def _norm_sign(x):
    for t in x:
        if t:
            return tuple(x) if t > 0 else tuple(-u for u in x)
    return tuple(x)


class ComplexCubicVoronoi:
    """Voronoi's algorithm on the circle dK_1(m_O)/O^x for a field of signature (1, 1)."""

    def __init__(self, pol, max_classes=100000):
        self.K = ComplexCubic(pol)
        self.max_classes = max_classes
        self.classes = []        # (A, minimal vectors, extreme pair)
        self.edges = []          # (class i, endpoint index, class j, eps)

    def _below(self, A, one=Fraction(1)):
        """Minimal vectors of A with value < 1 (exactly)."""
        K = self.K
        out = []
        for x in K.candidates(A, 1):
            v = K.value(A, x)
            if K.sign(tuple(Fraction(p) - (one if i == 0 else 0) for i, p in enumerate(v))) < 0:
                out.append(x)
        return out

    def _scale_to_min_one(self, A):
        K = self.K
        mn, S = K.minimum(A)
        inv = K.inverse(mn)
        return (K.mul(A[0], inv), K.mul(A[1], inv)), S

    def _q(self, x):
        K = self.K
        return K.s(x), tuple(2 * c for c in K.t(x))

    def _pair(self, S):
        """The two extreme minimal vectors: extremes of the ratio q_1/q_2 (exactly)."""
        K = self.K

        def less(x, y):     # q1(x)/q2(x) < q1(y)/q2(y)  <=>  q1(x) q2(y) - q1(y) q2(x) < 0
            qx, qy = self._q(x), self._q(y)
            d = tuple(a - b for a, b in zip(K.mul(qx[0], qy[1]), K.mul(qy[0], qx[1])))
            return K.sign(d) < 0
        lo = hi = S[0]
        for x in S[1:]:
            if less(x, lo):
                lo = x
            if less(hi, x):
                hi = x
        if lo == hi:
            raise ComplexCubicError('form is not perfect')
        return lo, hi

    def direction(self, x_keep, x_drop=None):
        """beta with <beta, q(x_keep)> = 0, oriented so that <beta, q(x_drop)> > 0."""
        K = self.K
        q_keep = self._q(x_keep)
        B = (q_keep[1], tuple(-c for c in q_keep[0]))
        if x_drop is not None:
            q_drop = self._q(x_drop)
            val = tuple(a + b for a, b in zip(K.mul(B[0], q_drop[0]), K.mul(B[1], q_drop[1])))
            if K.sign(val) < 0:
                B = (tuple(-c for c in B[0]), tuple(-c for c in B[1]))
        return B

    def walk(self, A, B):
        """First form A + t B (t > 0, possibly t in K) at which a vector not minimal for A attains
        the minimum 1; <B, q(x)> = 0 on the minimal vectors kept along the edge."""
        K = self.K

        def at(t):
            return (tuple(Fraction(a) + c for a, c in zip(A[0], K.mul(t, B[0]))),
                    tuple(Fraction(a) + c for a, c in zip(A[1], K.mul(t, B[1]))))

        def pos(F):
            return K.sign(F[0]) > 0 and K.sign(F[1]) > 0

        a1, a2, b1, b2 = K.real(A[0]), K.real(A[1]), K.real(B[0]), K.real(B[1])
        neg = [-a / b for a, b in ((a1, b1), (a2, b2)) if b < 0]
        if not neg:
            raise ComplexCubicError('unbounded direction')
        tmax = Fraction(min(neg))
        lo, t = Fraction(0), tmax / 2
        while True:
            F = at((t, 0, 0))
            if not pos(F):
                t = (lo + t) / 2
                continue
            S = self._below(F)
            if S:
                break
            lo, t = t, (t + tmax) / 2
        best = (t, Fraction(0), Fraction(0))
        while S:
            for x in S:
                qx = self._q(x)
                aval = tuple(p + q for p, q in zip(K.mul(A[0], qx[0]), K.mul(A[1], qx[1])))
                bval = tuple(p + q for p, q in zip(K.mul(B[0], qx[0]), K.mul(B[1], qx[1])))
                if K.sign(bval) >= 0:
                    raise ComplexCubicError('vector below 1 with <beta, q(x)> >= 0')
                num = tuple(Fraction(c) - (1 if i == 0 else 0) for i, c in enumerate(aval))
                tx = K.mul(num, K.inverse(tuple(-c for c in bval)))
                if K.sign(tuple(p - q for p, q in zip(tx, best))) < 0:
                    best = tx
            F = at(best)
            S = self._below(F)
        return F

    def _equiv(self, A, P):
        """eps in O^x with A = P eps-action, or None."""
        K = self.K
        u = K.mul(A[0], K.inverse(P[0]))
        if any(Fraction(c).denominator != 1 for c in u) or abs(K.norm(u)) != 1:
            return None
        tvar = pari.varhigher('t')
        R = pari.nfroots(K.nf, tvar ** 2 - K.to_alg(u))
        if len(R) == 0:
            return None
        eps = K.from_alg(R[0])
        if any(c.denominator != 1 for c in eps):
            return None
        for e in (eps, tuple(-c for c in eps)):
            if K.mul(P[1], K.t(e)) == tuple(Fraction(c) for c in A[1]):
                return tuple(int(c) for c in e)
        return None

    def run(self):
        K = self.K
        A = ((Fraction(1), Fraction(0), Fraction(0)), (Fraction(1), Fraction(0), Fraction(0)))
        A, S = self._scale_to_min_one(A)
        # ascend: walk until two minimal vectors with independent q
        while True:
            try:
                lo, hi = self._pair(S)
                break
            except ComplexCubicError:
                A = self.walk(A, self.direction(S[0]))
                mn, S = K.minimum(A)
        self.classes.append((A, S, (lo, hi)))
        i = 0
        while i < len(self.classes):
            A, S, (lo, hi) = self.classes[i]
            for k, (keep, drop) in enumerate([(lo, hi), (hi, lo)]):
                nb = self.walk(A, self.direction(keep, drop))
                j, eps = None, None
                for jj, (P, _, _) in enumerate(self.classes):
                    e = self._equiv(nb, P)
                    if e is not None:
                        j, eps = jj, e
                        break
                if j is None:
                    mn, Sn = K.minimum(nb)
                    if K.sign(tuple(Fraction(c) - (1 if t == 0 else 0) for t, c in enumerate(mn))) != 0:
                        raise ComplexCubicError('neighbour has minimum != 1')
                    j, eps = len(self.classes), K.one
                    self.classes.append((nb, Sn, self._pair(Sn)))
                    if j > self.max_classes:
                        raise ComplexCubicError('too many classes')
                self.edges.append((i, k, j, eps))
            i += 1
        return self

    def unit_index(self):
        K = self.K
        if K._bnf is None:
            pari.setrand(1)
            K._bnf = pari.bnfinit(K.pol, 1)
        units = [e for (_, _, _, e) in self.edges]
        ex = []
        for u in units:
            e = pari.bnfisunit(K._bnf, K.to_alg(u))
            ex.append(int(e[0]))
        from math import gcd
        g = 0
        for e in ex:
            g = gcd(g, e)
        return g
