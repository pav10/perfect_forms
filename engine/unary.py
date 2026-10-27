"""Voronoi's algorithm for perfect unary forms over a totally real field, exactly.

A unary form is a totally positive alpha in K; its value on x in O is
Tr(alpha x^2), a positive definite m-ary Z-form.  A perfect form of minimum 1 is
an element of K (Leibak; manuscript Lemma 5.2), so every vertex of the Ryshkov
polyhedron K_1(m_O) is represented exactly.  Two forms are equivalent iff
alpha' = alpha eta^2 with eta in O^x.

The algorithm (a generalisation of voronoi_unary.gp from cubic maximal orders to
any degree m and any order):

1. ascend from alpha = 1/m along directions orthogonal to the squared minimal
   vectors until they have rank m (a perfect form);
2. at a perfect form, compute the face lattice of its Voronoi cone (polyface.py);
   for every facet walk along the edge of K_1(m_O) to the neighbour (``walk``);
3. identify the neighbour with a known class exactly (``Voronoi._identify``) or
   record a new class; repeat until all neighbours are identified.

Exactness: minimal vectors by the integer Fincke--Pohst of qfmin.py; positive
definiteness by exact Gauss reduction; equivalence by exact arithmetic in K
(the unit eta is found by pairing minimal vectors and checked by
alpha' == alpha eta^2).  The units of bnfinit enter no decision: they are only used
to keep representatives small and to compute the unit-index certificate.

Certificates (certify.py): the alternating sum of the global face numbers of the
torus complex is 0, and the side-pairing units generate O^x / +-1 (index 1 in
bnf's unit group), which proves that bnf's units are fundamental.
"""
import time
from fractions import Fraction

import numpy as np

from . import linalg, polyface, qfmin, walk


class VoronoiError(ArithmeticError):
    pass


def _neg(x):
    return tuple(-t for t in x)


def normalize(x):
    """Representative of {x, -x}: first nonzero coordinate positive."""
    for t in x:
        if t:
            return tuple(x) if t > 0 else _neg(x)
    return tuple(x)


class PerfectForm:
    """A perfect form alpha (minimum 1) with its minimal vectors and face lattice."""

    def __init__(self, O, alpha, minvecs=None):
        self.alpha = tuple(Fraction(t) for t in alpha)
        if minvecs is None:
            mn, minvecs = qfmin.minimum(O.gram(self.alpha))
            if mn != 1:
                raise VoronoiError('minimum %s != 1' % mn)
        self.minvecs = [normalize(x) for x in minvecs]
        self.index = {x: i for i, x in enumerate(self.minvecs)}
        self.rows = [O.form_row(x) for x in self.minvecs]
        if linalg.rank(self.rows) != O.m:
            raise VoronoiError('form is not perfect')
        self.facets = polyface.facets(self.rows)
        self.faces = polyface.face_lattice(self.rows, self.facets)
        if not polyface.check_polytope_euler(self.faces, O.m):
            raise VoronoiError('local face lattice fails Euler--Poincare')
        self.vertices = sorted(i for F in self.faces[0] for i in F)
        self.key = O.ideal_key(self.alpha)
        self.norm = O.norm(self.alpha)
        self.neighbours = [None] * len(self.facets)      # (class j, eta)
        self.face_id = {}                                 # frozenset of vectors -> (dim, idx)
        for d, L in self.faces.items():
            for t, F in enumerate(L):
                self.face_id[frozenset(self.minvecs[i] for i in F)] = (d, t)

    @property
    def degree(self):
        return len(self.facets)


class Voronoi:
    """Voronoi's algorithm modulo units for the order O (field.Order)."""

    def __init__(self, O, max_classes=10 ** 6, reduce=True, log=None):
        self.O = O
        self.max_classes = max_classes
        self.reduce = reduce
        self.log = log
        self.classes = []
        self.by_key = {}
        self.stats = {'walks': 0, 'enumerations': 0, 'equivalence_tests': 0}
        self._units = None
        if reduce and O.m > 1:
            self._setup_reduction()

    # ------------------------------------------------------------------
    # size reduction modulo unit squares (heuristic: affects heights only)

    def _setup_reduction(self):
        O = self.O
        self._units = O.reduced_units()
        L = np.array([O.log_embed(u) for u in self._units])          # (m-1) x m
        self._Lpinv = np.linalg.pinv(2 * L.T)                          # m x (m-1) -> least squares
        self._unit_inv = [O.as_int(O.inverse(u)) for u in self._units]

    def reduce_form(self, alpha):
        """(alpha', w) with alpha = alpha' w^2 and alpha' of balanced embeddings."""
        O = self.O
        if self._units is None:
            return alpha, O.one
        la = O.log_embed(alpha)
        la = la - la.mean()
        k = np.rint(self._Lpinv @ la).astype(int)
        w = O.one
        for u, ui, e in zip(self._units, self._unit_inv, k):
            e = int(e)
            if e > 0:
                w = O.mul(w, O.power(u, e))
            elif e < 0:
                w = O.mul(w, O.power(ui, -e))
        if w == O.one:
            return alpha, w
        winv2 = O.inverse(O.mul(w, w))
        return O.mul(alpha, winv2), w

    # ------------------------------------------------------------------

    def walk(self, alpha, beta, F):
        """From the form alpha (minimum 1) along beta (Tr(beta x^2) = 0 for x in F, > 0
        for the other minimal vectors; beta not totally nonnegative) to the first form
        alpha + rho beta at which a vector outside F attains the minimum 1 (walk.py)."""
        O = self.O
        self.stats['walks'] += 1
        # estimate of the positivity boundary min_{sigma_i(beta) < 0} sigma_i(alpha)/|sigma_i(beta)|,
        # in log space (coordinates can exceed the float range when the regulator is large)
        ea, eb = O.embed_hp(alpha), O.embed_hp(beta)
        neg = [float(-a / b) for a, b in zip(ea, eb) if b < 0]
        if not neg:
            raise VoronoiError('walk: unbounded direction')
        try:
            rho = walk.walk(O.gram(alpha), O.gram(beta),
                            lambda x: O.trace_pair(alpha, O.sq(x)), lambda x: O.trace_pair(beta, O.sq(x)),
                            self.stats, tmax=min(neg))
        except walk.WalkError as e:
            raise VoronoiError(str(e))
        return tuple(Fraction(a) + rho * b for a, b in zip(alpha, beta))

    def initial_perfect_form(self):
        O = self.O
        alpha = tuple([Fraction(1, O.m)] + [Fraction(0)] * (O.m - 1))
        mn, S = qfmin.minimum(O.gram(alpha))
        alpha = tuple(a / mn for a in alpha)
        S = [normalize(x) for x in S]
        while True:
            rows = [O.form_row(x) for x in S]
            if linalg.rank(rows) == O.m:
                return alpha, S
            beta = linalg.kernel(rows)[0]
            alpha = self.walk(alpha, beta, S)
            mn, S = qfmin.minimum(O.gram(alpha))
            if mn != 1:
                raise VoronoiError('ascent: minimum %s != 1' % mn)
            S = [normalize(x) for x in S]

    # ------------------------------------------------------------------
    # equivalence

    def _find_eta(self, alpha, Fvecs, P):
        """eta in O^x with alpha = P.alpha eta^2, or None.  If such eta exists, eta x lies
        in S(P.alpha) (up to sign) for every x in S(alpha), in particular for x in F."""
        O = self.O
        self.stats['equivalence_tests'] += 1
        x1 = min(Fvecs, key=lambda v: (sum(abs(t) for t in v), v))
        n1 = abs(O.norm(x1))
        for y in P.minvecs:
            if abs(O.norm(y)) != n1:
                continue
            eta = O.div(y, x1)
            if not O.integral(eta):
                continue
            eta = O.as_int(eta)
            if abs(O.norm(eta)) != 1:
                continue
            if O.mul(P.alpha, O.sq(eta)) == alpha:
                return eta
        return None

    def _identify(self, alpha, Fvecs):
        """(class index, eta) with alpha = classes[j].alpha * eta^2; creates a new class
        if alpha is equivalent to none of the known ones."""
        O = self.O
        key = O.ideal_key(alpha)
        for j in self.by_key.get(key, ()):
            eta = self._find_eta(alpha, Fvecs, self.classes[j])
            if eta is not None:
                return j, eta
        if self.reduce:
            ared, w = self.reduce_form(alpha)
        else:
            ared, w = alpha, O.one
        P = PerfectForm(O, ared)
        if P.key != key:
            raise VoronoiError('reduction changed the ideal key')
        j = len(self.classes)
        self.classes.append(P)
        self.by_key.setdefault(key, []).append(j)
        if len(self.classes) > self.max_classes:
            raise VoronoiError('more than %d classes' % self.max_classes)
        return j, w

    # ------------------------------------------------------------------

    def run(self):
        O = self.O
        t0 = time.time()
        alpha, S = self.initial_perfect_form()
        self._identify(alpha, S)
        i = 0
        while i < len(self.classes):
            P = self.classes[i]
            for f, (beta, F) in enumerate(P.facets):
                Fvecs = [P.minvecs[k] for k in F]
                nb = self.walk(P.alpha, beta, Fvecs)
                j, eta = self._identify(nb, Fvecs)
                P.neighbours[f] = (j, eta)
            i += 1
            if self.log and i % 100 == 0:
                self.log('%d classes processed, %d known, %.1fs' % (i, len(self.classes), time.time() - t0))
        self.time = time.time() - t0
        self._glue()
        return self

    # ------------------------------------------------------------------
    # global face complex of the torus  dK_1(m_O) / O^x

    def _glue(self):
        """Union-find of the local faces of all classes along the side pairings.

        A local face is (class i, frozenset of minimal vectors).  Across the facet F of
        class i the neighbour equals classes[j] * eta^2, and x -> eta x maps the faces
        G subset F of P(alpha_i) onto faces of P(alpha_j).  Faces of a polyhedral
        manifold are glued exactly by such facet adjacencies, so the classes of the
        union-find are the cells of the torus complex (dual numbering: dimension d of
        the Voronoi-cone face is the codimension of the face of K_1(m_O))."""
        O = self.O
        parent = {}

        def find(a):
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a

        def union(a, b):
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[max(ra, rb)] = min(ra, rb)

        for i, P in enumerate(self.classes):
            for d, L in P.faces.items():
                for t in range(len(L)):
                    parent[(i, d, t)] = (i, d, t)
        for i, P in enumerate(self.classes):
            for f, (beta, F) in enumerate(P.facets):
                j, eta = P.neighbours[f]
                Q = self.classes[j]
                M = O.mult_matrix(eta)
                img = {}
                for k in F:
                    x = P.minvecs[k]
                    y = normalize(tuple(linalg.matvec(M, x)))
                    if y not in Q.index:
                        raise VoronoiError('side pairing does not map minimal vectors')
                    img[k] = y
                for d, L in P.faces.items():
                    for t, G in enumerate(L):
                        if G <= F:
                            key = frozenset(img[k] for k in G)
                            dt = Q.face_id.get(key)
                            if dt is None or dt[0] != d:
                                raise VoronoiError('side pairing does not map faces to faces')
                            union((i, d, t), (j, d, dt[1]))
        self.cells = {}
        for a in parent:
            r = find(a)
            self.cells.setdefault(r[1], {}).setdefault(r, []).append(a)

    # ------------------------------------------------------------------
    # invariants

    def f_vector(self):
        """(f_0, ..., f_{m-1}) of the dual complex: f_d = number of d-faces of the
        Voronoi-cone complex modulo units; f_0 = mu (extreme squares), f_{m-1} = n."""
        m = self.O.m
        f = [len(self.cells.get(d, {})) for d in range(m - 1)]
        return f + [len(self.classes)]

    def euler(self):
        return sum((-1) ** d * v for d, v in enumerate(self.f_vector()))

    def side_pairing_units(self):
        return [eta for P in self.classes for (j, eta) in P.neighbours]

    def degrees(self):
        return sorted(P.degree for P in self.classes)

    def facet_sizes(self):
        """Number of vertices of each facet of K_1(m_O) mod units (for m = 3: the
        polygon sizes): occurrences of each global vertex of the dual complex."""
        return sorted(len(v) for v in self.cells.get(0, {}).values())

    def extreme_square_classes(self):
        """Representatives of the extreme squares mod units, keyed by the ideal xO
        (exact: xO = yO iff y/x in O^x).  An independent count of mu."""
        O = self.O
        reps = {}
        for P in self.classes:
            for k in P.vertices:
                x = P.minvecs[k]
                reps.setdefault(O.ideal_key(x), x)
        return reps
