"""Voronoi's algorithm for binary Humbert forms in the additive normalisation (plan H1).

Over a totally real field K of degree m with maximal order O, a binary form is
A = [[a, b], [b, c]] in Sym_2(K (x) R), totally positive definite, with value
q_A(v) = Tr_{K/Q}(v^T A v) = Tr(a x^2 + 2 b x y + c y^2) on v = (x, y) in O^2
(Koecher, Leibak's additive theory, Gunnells--Yasaki).  On the Z-basis
(w_1 e_1, ..., w_m e_1, w_1 e_2, ..., w_m e_2) of O^2 this is a positive definite
2m-ary Z-form; a perfect form of minimum 1 has a, b, c in K (it solves a rational
linear system), so everything is exact.  The space of forms has dimension 3m.

Equivalence is modulo GL_2(O): A' ~ A iff A' = g^T A g, g in GL_2(O) (homothety is
fixed by minimum 1).  It is decided exactly with PARI's qfisom applied to the m
Z-forms Tr(theta^k v^T A v), k = 0..m-1, for a totally positive generator theta of
K: a Z-isometry preserving all of them preserves the K-valued form v^T A w, hence is
O-linear.  Automorphism groups (stabilisers) come from qfauto on the same forms.

Labels: complete lists of perfect forms are CERTIFIED computations (exact walk,
exact equivalence); there is no Euler-characteristic certificate here yet (the
quotient is an orbifold with cusps; see README).  Additive and multiplicative
(Icaza/BCIO) normalisations differ and are never mixed (plan H1).
"""
import time
from fractions import Fraction

from . import linalg, polyface, qfmin, walk
from .pari import pari


class BinaryError(ArithmeticError):
    pass


def _normalize(v):
    for t in v:
        if t:
            return tuple(v) if t > 0 else tuple(-u for u in v)
    return tuple(v)


class BinaryForms:
    """The space Sym_2(K) of binary forms over the maximal order O (field.Order)."""

    def __init__(self, O):
        self.O = O
        m = self.m = O.m
        self.n = 2 * m          # rank of the Z-lattice O^2
        self.N = 3 * m          # dimension of the space of forms
        # a totally positive generator theta of K (for the isometry test)
        k = 0
        while True:
            theta = tuple([k] + [1] + [0] * (m - 2)) if m > 1 else (1,)
            if m == 1 or (O.totally_positive(theta) and self._generates(theta)):
                break
            k += 1
        self.theta = theta
        self.theta_pows = [O.power(theta, j) for j in range(m)]

    def _generates(self, t):
        O = self.O
        P = pari.minpoly(O.to_alg(t))
        return int(pari.poldegree(P)) == O.m

    # --- coordinates: A = (a_1..a_m, b_1..b_m, c_1..c_m) -------------------

    def split(self, A):
        m = self.m
        return A[:m], A[m:2 * m], A[2 * m:]

    def gram(self, A):
        O, m = self.O, self.m
        a, b, c = self.split(A)
        Ga, Gb, Gc = O.gram(a), O.gram(b), O.gram(c)
        G = [[Fraction(0)] * (2 * m) for _ in range(2 * m)]
        for j in range(m):
            for k in range(m):
                G[j][k] = Ga[j][k]
                G[j][m + k] = G[m + k][j] = Gb[j][k]
                G[m + j][m + k] = Gc[j][k]
        return G

    def vec(self, v):
        m = self.m
        return tuple(v[:m]), tuple(v[m:])

    def row(self, v):
        O = self.O
        x, y = self.vec(v)
        T = O.T
        m = self.m
        sx, sy, xy = O.sq(x), O.sq(y), O.mul(x, y)
        r1 = [sum(sx[p] * T[l][p] for p in range(m)) for l in range(m)]
        r2 = [2 * sum(xy[p] * T[l][p] for p in range(m)) for l in range(m)]
        r3 = [sum(sy[p] * T[l][p] for p in range(m)) for l in range(m)]
        return tuple(r1 + r2 + r3)

    def value(self, A, v):
        return linalg.dot(self.row(v), A)

    def det(self, A):
        O = self.O
        a, b, c = self.split(A)
        return tuple(s - t for s, t in zip(O.mul(a, c), O.mul(b, b)))

    def inverse_row(self, A):
        """The row vector of the functional H -> Tr(tr(A^{-1} H)) (gradient of log N det)."""
        O, m, T = self.O, self.m, self.O.T
        a, b, c = self.split(A)
        dinv = O.inverse(self.det(A))
        p, q, s = O.mul(c, dinv), tuple(-t for t in O.mul(b, dinv)), O.mul(a, dinv)
        r1 = [sum(p[k] * T[l][k] for k in range(m)) for l in range(m)]
        r2 = [2 * sum(q[k] * T[l][k] for k in range(m)) for l in range(m)]
        r3 = [sum(s[k] * T[l][k] for k in range(m)) for l in range(m)]
        return tuple(r1 + r2 + r3)

    def scaled(self, lam, A):
        """lam * A for lam in K (coordinates)."""
        O = self.O
        a, b, c = self.split(A)
        return tuple(O.mul(lam, a)) + tuple(O.mul(lam, b)) + tuple(O.mul(lam, c))

    def iso_forms(self, A):
        """The integral Z-forms Tr(theta^k v^T A v), k = 0..m-1, with a common scale."""
        Gs = [self.gram(self.scaled(t, A)) for t in self.theta_pows]
        den = linalg.common_denominator([x for G in Gs for r in G for x in r])
        n = self.n
        return [pari.matrix(n, n, [int(x * den) for r in G for x in r]) for G in Gs], den


class BinaryPerfectForm:
    def __init__(self, B, A, minvecs=None):
        self.A = tuple(Fraction(t) for t in A)
        if minvecs is None:
            mn, minvecs = qfmin.minimum(B.gram(self.A))
            if mn != 1:
                raise BinaryError('minimum %s != 1' % mn)
        self.minvecs = [_normalize(v) for v in minvecs]
        self.index = {v: i for i, v in enumerate(self.minvecs)}
        self.rows = [B.row(v) for v in self.minvecs]
        if linalg.rank(self.rows) != B.N:
            raise BinaryError('form is not perfect')
        self.facets = polyface.facets(self.rows)
        self.normdet = B.O.norm(B.det(self.A))
        self.forms, self.den = B.iso_forms(self.A)
        aut = pari.qfauto(pari(self.forms))
        self.aut_order = int(aut[0])
        self.neighbours = [None] * len(self.facets)

    @property
    def key(self):
        return (self.normdet, len(self.minvecs))


class BinaryVoronoi:
    def __init__(self, O, max_classes=10 ** 5, log=None):
        self.O = O
        self.B = BinaryForms(O)
        self.max_classes = max_classes
        self.log = log
        self.classes = []
        self.by_key = {}
        self.stats = {'walks': 0, 'enumerations': 0, 'isometry_tests': 0}

    def walk(self, A, beta):
        B = self.B
        self.stats['walks'] += 1
        rho = walk.walk(B.gram(A), B.gram(beta), lambda v: B.value(A, v), lambda v: B.value(beta, v), self.stats)
        return tuple(Fraction(a) + rho * b for a, b in zip(A, beta))

    def initial_perfect_form(self):
        B = self.B
        m = B.m
        one = [Fraction(1)] + [Fraction(0)] * (m - 1)
        A = tuple(one + [Fraction(0)] * m + one)
        mn, S = qfmin.minimum(B.gram(A))
        A = tuple(t / mn for t in A)
        while True:
            rows = [B.row(v) for v in S]
            if linalg.rank(rows) == B.N:
                return A
            beta = linalg.kernel(rows)[0]
            A = self.walk(A, beta)
            mn, S = qfmin.minimum(B.gram(A))
            if mn != 1:
                raise BinaryError('ascent: minimum %s != 1' % mn)

    def _identify(self, A):
        B = self.B
        mn, S = qfmin.minimum(B.gram(A))
        if mn != 1:
            raise BinaryError('neighbour: minimum %s != 1' % mn)
        key = (B.O.norm(B.det(A)), len(S))
        forms = None
        for j in self.by_key.get(key, ()):
            P = self.classes[j]
            if forms is None:
                forms, den = B.iso_forms(A)
            if den != P.den:
                continue
            self.stats['isometry_tests'] += 1
            # g with A = g^T P.A g  (as Z-forms: forms = g~ P.forms g)
            g = pari.qfisom(pari(forms), pari(P.forms))
            if g != 0:
                n = B.n
                g = [[int(g[i, k]) for k in range(n)] for i in range(n)]
                for F, G in zip(forms, P.forms):
                    if pari.matrix(n, n, [x for r in g for x in r]).mattranspose() * G * \
                            pari.matrix(n, n, [x for r in g for x in r]) != F:
                        raise BinaryError('qfisom returned a non-isometry')
                return j, g
        P = BinaryPerfectForm(B, A, S)
        j = len(self.classes)
        self.classes.append(P)
        self.by_key.setdefault(key, []).append(j)
        if j >= self.max_classes:
            raise BinaryError('too many classes')
        n = B.n
        return j, [[int(i == k) for k in range(n)] for i in range(n)]

    def run(self):
        t0 = time.time()
        A = self.initial_perfect_form()
        self._identify(A)
        i = 0
        while i < len(self.classes):
            P = self.classes[i]
            for f, (beta, F) in enumerate(P.facets):
                nb = self.walk(P.A, beta)
                P.neighbours[f] = self._identify(nb)
            i += 1
            if self.log:
                self.log('%d classes processed, %d known, %.1fs' % (i, len(self.classes), time.time() - t0))
        self.time = time.time() - t0
        return self

    # --- invariants -------------------------------------------------------

    def eutaxy(self, P):
        r = self.B.inverse_row(P.A)
        vals = [linalg.dot(beta, r) for beta, _ in P.facets]
        if all(v > 0 for v in vals):
            return 'eutactic'
        if all(v >= 0 for v in vals):
            return 'weakly'
        return 'no'

    def summary(self):
        out = []
        m = self.B.m
        for i, P in enumerate(self.classes):
            g2m = 1 / Fraction(P.normdet)          # gamma(A)^(2m) = m(A)^(2m) / N(det A), m(A) = 1
            out.append({'index': i, 'A': [str(t) for t in P.A], 'n_minvecs': len(P.minvecs),
                        'n_facets': len(P.facets), 'aut_order': P.aut_order,
                        'normdet': str(P.normdet), 'gamma_pow': str(g2m), 'gamma': float(g2m) ** (1 / (2 * m)),
                        'eutaxy': self.eutaxy(P), 'neighbours': sorted(set(j for j, _ in P.neighbours))})
        return out


def edge_balance(V):
    """Certificate: for classes i, j let N_ij be the number of facets of the Voronoi cone
    of A_i whose neighbour lies in class j.  Counting the GL_2(O)-orbits of edges
    {A_i, A_j} with their stabilisers gives N_ij / |Aut A_i| = N_ji / |Aut A_j|
    (sum over the edge orbits of 1/|Stab(edge)|).  Returns the list of violations."""
    from collections import Counter
    N = [Counter(j for j, _ in P.neighbours) for P in V.classes]
    bad = []
    for i, P in enumerate(V.classes):
        for j in range(len(V.classes)):
            lhs = Fraction(N[i][j], P.aut_order)
            rhs = Fraction(N[j][i], V.classes[j].aut_order)
            if lhs != rhs:
                bad.append((i, j, lhs, rhs))
    return bad
