"""Relative minima of the lattice psi(O) modulo units (Voronoi 1896, Buchmann).

x in O \\ 0 is a relative minimum if no y in O \\ {0, +-x} has |sigma_i(y)| <=
|sigma_i(x)| for all i (manuscript Def. 6.2).  Two minima x, y are adjacent if the
box B(x v y) = {z : |z_i| <= max(|x_i|, |y_i|)} contains no lattice point other
than 0, +-x, +-y.  The adjacency graph is connected (follow a generic path of box
shapes a in R^m_{>0}: the minimiser of max_i |z_i|/a_i changes only at ties of two
minima, which are adjacent), and units act freely on it, so a walk from x = 1
modulo units finds all ell_O classes, and the units relating neighbours to
class representatives generate O^x / +-1.

Neighbours of x lie in the box B* = prod [-b_k, b_k], b_k = |sigma_k(N_k(x))|,
where N_k(x) minimises |sigma_k| subject to |sigma_j| < |sigma_j(x)| (j != k):
if y is adjacent to x and |y_k| > |x_k|, then |y_k| <= b_k, since otherwise
N_k(x) lies in B(x v y).

Lattice points in a box are enumerated exactly as {y : Tr(beta y^2) <= R} with
beta totally positive (checked exactly) and R a rational upper bound of
sum_j sigma_j(beta) b_j^2; all comparisons |sigma_i(y)| vs |sigma_i(z)| are exact
(field.Order.compare_abs: certified floating point, else Sturm sequences).  This
code shares only the field arithmetic and the enumerator with unary.py.
"""
from fractions import Fraction

import numpy as np

from . import linalg, qfmin
from .pari import pari


class RelMinError(ArithmeticError):
    pass


def _neg(x):
    return tuple(-t for t in x)


class _Pts:
    """Lattice points with float embeddings and certified error bounds, for fast
    exact comparisons of |sigma_i|."""

    def __init__(self, O, pts):
        self.O = O
        self.pts = list(pts)
        X = np.array([[float(t) for t in p] for p in self.pts]) if self.pts else np.zeros((0, O.m))
        self.ab = np.abs(X @ O.E.T)
        self.err = 1e-12 * (np.abs(X) @ np.abs(O.E).T) + 1e-300
        self.big = [not all(abs(t) < 2 ** 50 for t in p) for p in self.pts]

    def less(self, a, b, i):
        """|sigma_i(p_a)| < |sigma_i(p_b)| (exactly)."""
        if not (self.big[a] or self.big[b]):
            if self.ab[a, i] + self.err[a, i] < self.ab[b, i] - self.err[b, i]:
                return True
            if self.ab[a, i] - self.err[a, i] > self.ab[b, i] + self.err[b, i]:
                return False
        return self.O.compare_abs(self.pts[a], self.pts[b], i) < 0


class RelativeMinima:
    def __init__(self, O, max_classes=10 ** 6):
        self.O = O
        self.max_classes = max_classes
        self.reps = []           # class representatives
        self.by_key = {}
        self.neighbours = []     # per class: list of (class j, unit eps) with nb = eps * rep_j
        self.stats = {'enumerations': 0, 'points': 0}

    # --- boxes -----------------------------------------------------------

    def _weight(self, targets):
        """Totally positive beta with sigma_j(beta) ~ targets[j] > 0 (heuristic quality;
        total positivity is checked exactly).  The linear system E c = targets is solved
        in PARI at a precision adapted to the spread of the targets."""
        O = self.O
        lt = np.log2(np.array(targets, dtype=float))
        spread = float(lt.max() - lt.min())
        p = max(8, int(np.ceil(-lt.min() + np.log2(O.m * np.abs(O.E).max()))) + 8)
        bits = int(spread + p + 128)
        while True:
            c = O.solve_embeddings(targets, bits)
            emax = max([int(pari.exponent(t)) for t in c if t != 0] or [0])
            if p + emax + 32 > bits:
                bits = p + emax + 96
                c = O.solve_embeddings(targets, bits)
            beta = tuple(Fraction(int(pari.round(t * pari(2) ** p)), 2 ** p) for t in c)
            if any(beta) and O.totally_positive(beta):
                return beta
            p += 16
            bits += 64
            if p > 20000:
                raise RelMinError('could not build a positive weight')

    def _enum(self, beta, R):
        self.stats['enumerations'] += 1
        pts = [x for _, x in qfmin.short_vectors(self.O.gram(beta), R)]
        self.stats['points'] += len(pts)
        return pts

    def _upper(self, beta, w, j):
        """Rational upper bound of sigma_j(beta w^2) (w in O)."""
        O = self.O
        y = O.mul(beta, O.sq(w))
        den = linalg.common_denominator(y)
        yi = tuple(int(t * den) for t in y)
        fb = O._float_bound(yi, j)
        if fb is not None:
            return Fraction(fb[0] + 2 * fb[1] + 1e-300) / den
        return O.trace_pair(beta, O.sq(w))            # Tr >= each (positive) term

    def single_neighbour(self, x, k):
        """N_k(x): the y minimising |sigma_k(y)| subject to |sigma_j(y)| < |sigma_j(x)|
        for all j != k."""
        O = self.O
        ax = np.abs(O.embed_accurate(x))
        D = abs(O.disc) ** 0.5
        T = 2 * D / np.prod([ax[j] for j in range(O.m) if j != k])   # Minkowski: a point exists
        targets = [1 / ax[j] ** 2 for j in range(O.m)]
        targets[k] = 1 / T ** 2
        beta = self._weight(targets)
        base = O.trace_pair(beta, O.sq(x))
        s = Fraction(1)
        while True:
            pts = self._enum(beta, base + s)
            P = _Pts(O, [x] + [p for p in pts if p != x and p != _neg(x)])
            best = None
            for a in range(1, len(P.pts)):
                if all(P.less(a, 0, j) for j in range(O.m) if j != k):
                    if best is None or P.less(a, best, k):
                        best = a
            if best is not None:
                y = P.pts[best]
                # guaranteed complete if sigma_k(beta y^2) <= s
                if self._le_at(O.mul(beta, O.sq(y)), s, k):
                    return y
            s *= 4

    def _le_at(self, w, s, k):
        """sigma_k(w) <= s for w in K (Fraction coordinates), exactly."""
        O = self.O
        diff = [Fraction(-t) for t in w]
        diff[0] += s                         # first basis vector is 1
        den = linalg.common_denominator(diff)
        di = tuple(int(t * den) for t in diff)
        return O.sign_at(di, k) >= 0

    def neighbours_of(self, x):
        O = self.O
        Ns = [self.single_neighbour(x, k) for k in range(O.m)]
        ab = [abs(O.embed_accurate(Ns[k])[k]) for k in range(O.m)]
        beta = self._weight([1 / b ** 2 for b in ab])
        R = sum(self._upper(beta, Ns[k], k) for k in range(O.m))
        pts = self._enum(beta, R)
        cand = [p for p in pts if p != x and p != _neg(x)]
        P = _Pts(O, [x] + cand)
        n = len(P.pts)
        out = []
        for a in range(1, n):
            # x not in B(y): some coordinate of y is smaller than that of x
            if not any(P.less(a, 0, i) for i in range(O.m)):
                continue
            ok = True
            for z in range(1, n):
                if z == a:
                    continue
                # z in B(x v y) ?
                inside = True
                for i in range(O.m):
                    big = 0 if P.less(a, 0, i) else a            # index of the larger of x, y
                    if P.less(big, z, i):
                        inside = False
                        break
                if inside:
                    ok = False
                    break
            if ok:
                out.append(P.pts[a])
        for k in range(O.m):
            if Ns[k] not in out and _neg(Ns[k]) not in out:
                raise RelMinError('single-direction neighbour not adjacent')
        return out

    # --- walk modulo units --------------------------------------------------

    def _reduce(self, x):
        O = self.O
        units = O.reduced_units()
        if not units:
            return x
        L = np.array([O.log_embed(u) for u in units])
        lx = O.log_embed(x)
        lx = lx - lx.mean()
        k = np.rint(np.linalg.pinv(L.T) @ lx).astype(int)
        w = O.one
        for u, e in zip(units, k):
            if e:
                w = O.mul(w, O.power(u, -int(e)))
        return O.as_int(O.mul(x, w))

    def _identify(self, y):
        O = self.O
        key = O.ideal_key(y)
        j = self.by_key.get(key)
        if j is not None:
            eps = O.div(y, self.reps[j])
            if not O.is_unit(eps):
                raise RelMinError('same ideal but quotient is not a unit')
            return j, O.as_int(eps)
        r = self._reduce(y)
        j = len(self.reps)
        self.reps.append(r)
        self.by_key[key] = j
        self.neighbours.append(None)
        if j >= self.max_classes:
            raise RelMinError('too many classes')
        return j, O.as_int(O.div(y, r))

    def run(self):
        O = self.O
        self._identify(O.one)
        i = 0
        while i < len(self.reps):
            nbs = self.neighbours_of(self.reps[i])
            self.neighbours[i] = [self._identify(y) for y in nbs]
            i += 1
        return self

    def count(self):
        return len(self.reps)

    def side_pairing_units(self):
        return [eps for L in self.neighbours for (_, eps) in L]

    def is_relative_minimum(self, x):
        """Direct exact test (for cross-checks): enumerate the box B(x)."""
        O = self.O
        ax = np.abs(O.embed_accurate(x))
        beta = self._weight([1 / a ** 2 for a in ax])
        R = sum(self._upper(beta, x, k) for k in range(O.m))
        pts = [p for p in self._enum(beta, R) if p != x and p != _neg(x)]
        P = _Pts(O, [x] + pts)
        for a in range(1, len(P.pts)):
            if not any(P.less(0, a, i) for i in range(O.m)):
                return False
        return True
