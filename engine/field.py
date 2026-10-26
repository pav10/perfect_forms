"""Orders in totally real number fields, with exact arithmetic.

Conventions
-----------
* An order O is given by a Z-basis w_1 = 1, w_2, ..., w_m.  For the maximal
  order this is PARI's integral basis nf.zk; for any other order it is a basis
  given in nf.zk coordinates (as in pav10/sails).
* Elements of K are coordinate vectors in this basis: tuples of int (elements of
  O) or of Fraction (elements of K).
* All decisions (membership, norms, units, squares, positivity, equality) are
  exact.  Floating point embeddings are used only to choose small
  representatives modulo units and to pick starting points of searches.
"""
import math
from fractions import Fraction

import numpy as np

from . import linalg
from .pari import pari

_MEMBERS = {}


def member(obj, name):
    f = _MEMBERS.get(name)
    if f is None:
        f = _MEMBERS[name] = pari('(o)->o.%s' % name)
    return f(obj)


def to_frac(g):
    return Fraction(int(pari.numerator(g)), int(pari.denominator(g)))


def parse_poly(pol):
    """Coefficient list (ascending, leading 1 included) or a polynomial string -> PARI t_POL."""
    if isinstance(pol, str):
        return pari(pol)
    return pari.Pol(list(reversed([int(c) for c in pol])))


class Order:
    """An order O of a totally real field K = Q[x]/(f)."""

    def __init__(self, pol, basis=None, label=None):
        self.pol = parse_poly(pol)
        self.label = label
        nf = self.nf = pari.nfinit(self.pol)
        r1, r2 = [int(t) for t in member(nf, 'sign')]
        self.signature = (r1, r2)
        if r2:
            raise NotImplementedError('only totally real fields are implemented so far')
        m = self.m = r1
        self.field_disc = int(member(nf, 'disc'))
        zk = member(nf, 'zk')
        if basis is None:
            self.basis_zk = [[int(i == j) for i in range(m)] for j in range(m)]
        else:
            self.basis_zk = [[int(t) for t in v] for v in basis]
            if self.basis_zk[0] != [1] + [0] * (m - 1):
                raise ValueError('the first basis vector of the order must be 1')
        self.maximal = basis is None
        # Mo: O-coordinates -> zk coordinates (columns = basis vectors); Minv: back
        self._Mo = pari.matrix(m, m, [self.basis_zk[j][i] for i in range(m) for j in range(m)])
        self._Minv = self._Mo ** -1
        self.index = abs(int(pari.matdet(self._Mo)))
        if self.index == 0:
            raise ValueError('order basis is not of full rank')
        self.disc = self.index ** 2 * self.field_disc
        self.basis_alg = [pari.nfbasistoalg(nf, pari.Col(self.basis_zk[j])) for j in range(m)]
        # multiplication table w_j w_k = sum_l mult[j][k][l] w_l
        self.mult = [[None] * m for _ in range(m)]
        for j in range(m):
            for k in range(j, m):
                v = self._from_zk(pari.nfeltmul(nf, pari.Col(self.basis_zk[j]), pari.Col(self.basis_zk[k])))
                if any(t.denominator != 1 for t in v):
                    raise ValueError('the given basis does not span a ring')
                self.mult[j][k] = self.mult[k][j] = tuple(int(t) for t in v)
        tr = [int(pari.nfelttrace(nf, pari.Col(self.basis_zk[l]))) for l in range(m)]
        self.tr = tuple(tr)
        # trace matrix T[j][k] = Tr(w_j w_k) and T3[l][j][k] = Tr(w_l w_j w_k)
        self.T = [[sum(self.mult[j][k][p] * tr[p] for p in range(m)) for k in range(m)] for j in range(m)]
        self.T3 = [[[sum(self.mult[j][k][p] * self.T[l][p] for p in range(m)) for k in range(m)]
                    for j in range(m)] for l in range(m)]
        self.one = tuple([1] + [0] * (m - 1))
        # float embeddings sigma_i(w_j), real places ordered by the root (increasing)
        roots = sorted(pari.polrootsreal(self.pol, precision=256), key=lambda r: float(r))
        self.roots = roots
        self.E = np.array([[float(pari.subst(pari.lift(self.basis_alg[j]), 'x', r)) for j in range(m)]
                           for r in roots])
        self._bnf = None
        self._units = None

    # ------------------------------------------------------------------
    # coordinates

    def _from_zk(self, v):
        """zk coordinates (PARI column, possibly rational) -> O coordinates (Fractions)."""
        v = pari.nfalgtobasis(self.nf, v)
        w = self._Minv * v
        return tuple(to_frac(t) for t in w)

    def to_zk(self, x):
        return self._Mo * pari.Col([pari(Fraction(t).numerator) / Fraction(t).denominator for t in x])

    def to_alg(self, x):
        return pari.nfbasistoalg(self.nf, self.to_zk(x))

    def from_alg(self, a):
        return self._from_zk(a)

    def integral(self, x):
        """Is x (Fraction coordinates) in O?"""
        return all(Fraction(t).denominator == 1 for t in x)

    def as_int(self, x):
        if not self.integral(x):
            raise ValueError('element not in the order')
        return tuple(int(t) for t in x)

    def poly_str(self, x):
        return str(pari.lift(self.to_alg(x)))

    # ------------------------------------------------------------------
    # arithmetic (exact, in Python)

    def mul(self, x, y):
        m = self.m
        out = [0] * m
        for j, xj in enumerate(x):
            if xj:
                for k, yk in enumerate(y):
                    if yk:
                        p = xj * yk
                        mjk = self.mult[j][k]
                        for l in range(m):
                            if mjk[l]:
                                out[l] += p * mjk[l]
        return tuple(out)

    def sq(self, x):
        return self.mul(x, x)

    def mult_matrix(self, x):
        """Matrix M with coords(x y) = M coords(y)."""
        m = self.m
        cols = [self.mul(x, tuple(int(i == k) for i in range(m))) for k in range(m)]
        return [[cols[k][i] for k in range(m)] for i in range(m)]

    def norm(self, x):
        return linalg.det(self.mult_matrix(x))

    def trace(self, x):
        return sum(t * c for t, c in zip(x, self.tr))

    def trace_pair(self, x, y):
        """Tr(x y)."""
        return sum(x[j] * self.T[j][k] * y[k] for j in range(self.m) if x[j] for k in range(self.m) if y[k])

    def inverse(self, x):
        return tuple(linalg.solve(self.mult_matrix(x), self.one))

    def div(self, x, y):
        return tuple(linalg.solve(self.mult_matrix(y), x))

    def power(self, x, e):
        if e < 0:
            x, e = self.inverse(x), -e
        r = self.one
        while e:
            if e & 1:
                r = self.mul(r, x)
            x = self.mul(x, x)
            e >>= 1
        return r

    def is_unit(self, x):
        """x in O^x (exact)."""
        return self.integral(x) and abs(self.norm(x)) == 1

    def totally_positive(self, x):
        """Exact: the characteristic polynomial of a totally real element has only
        positive roots iff its coefficients strictly alternate in sign."""
        cp = pari.charpoly(self.to_alg(x))
        c = [to_frac(t) for t in pari.Vec(cp)]           # leading coefficient first
        return all(ci != 0 and (ci > 0) == (i % 2 == 0) for i, ci in enumerate(c))

    def sqrt(self, x):
        """A square root of x in K (coordinates), or None if x is not a square in K."""
        t = pari.varhigher('t')
        R = pari.nfroots(self.nf, t ** 2 - self.to_alg(x))
        if len(R) == 0:
            return None
        return self.from_alg(R[0])

    def gram(self, alpha):
        """Gram matrix of x -> Tr(alpha x^2) on the basis of O (Fractions)."""
        m = self.m
        G = [[Fraction(0)] * m for _ in range(m)]
        for l, al in enumerate(alpha):
            if al:
                T3l = self.T3[l]
                for j in range(m):
                    for k in range(m):
                        G[j][k] += al * T3l[j][k]
        return G

    def form_row(self, x):
        """Integer vector r with Tr(alpha x^2) = r . alpha for every alpha."""
        s = self.sq(x)
        return tuple(sum(s[p] * self.T[l][p] for p in range(self.m)) for l in range(self.m))

    def ideal_key(self, alpha):
        """Canonical key of the Z-module alpha O (exact): (den, HNF of den M_alpha).
        alpha O = beta O  iff  beta / alpha in O^x."""
        M = self.mult_matrix(alpha)
        den = linalg.common_denominator([t for row in M for t in row])
        m = self.m
        H = pari.mathnf(pari.matrix(m, m, [int(t * den) for row in M for t in row]))
        return (den, tuple(int(H[i, j]) for i in range(m) for j in range(m)))

    # ------------------------------------------------------------------
    # embeddings (floating point: heuristics only)

    def embed(self, x):
        return self.E @ np.array([float(t) for t in x])

    def log_embed(self, x):
        """log|sigma_i(x)| in floating point (heuristic).  Falls back to PARI at a
        precision adapted to the size of the coordinates (huge units, tiny embeddings)."""
        bits = max(max(Fraction(t).numerator.bit_length(), Fraction(t).denominator.bit_length()) for t in x)
        if bits < 900:
            e = self.embed(x)
            scale = np.abs(self.E) @ np.abs(np.array([float(t) for t in x]))
            if np.all(np.abs(e) > 1e-8 * scale):
                return np.log(np.abs(e))
        a = pari.lift(self.to_alg(x))
        prec = max(2048, 2 * self.m * bits + 256)
        return np.array([float(pari.log(pari.abs(pari.subst(a, 'x', r)))) for r in self._hp_roots(prec)])

    def _hp_roots(self, bits=2048):
        if getattr(self, '_hproots_bits', 0) < bits:
            self._hproots = sorted(pari.polrootsreal(self.pol, precision=bits), key=lambda r: float(r))
            self._hproots_bits = bits
        return self._hproots

    # ------------------------------------------------------------------
    # class group and units (PARI's bnfinit: conditional on GRH unless certified)

    @property
    def bnf(self):
        if self._bnf is None:
            pari.setrand(1)
            self._bnf = pari.bnfinit(self.pol, 1)
        return self._bnf

    def class_number(self):
        return int(member(self.bnf, 'no'))

    def fundamental_units(self):
        """Fundamental units of O (coordinates), from bnfinit's units of O_K.

        For a non-maximal order the unit group is the kernel of
        O_K^x -> (O_K / f O_K)^x / (O / f O_K)^x, f the exponent of O_K/O
        (method of pav10/sails, field.py)."""
        if self._units is None:
            fu = list(member(self.bnf, 'fu'))
            if self.maximal or self.index == 1:
                self.unit_index = 1
                units = fu
                self.unit_exponents = [[int(i == j) for i in range(len(fu))] for j in range(len(fu))]
            else:
                units = self._order_units(fu)
            self._units = [self.as_int(self.from_alg(u)) for u in units]
        return self._units

    def reduced_units(self):
        """A basis of O^x / +-1 that is LLL-reduced in log space (heuristic use only:
        size reduction of representatives)."""
        if getattr(self, '_red_units', None) is None:
            units = self.fundamental_units()
            r = len(units)
            if r <= 1:
                self._red_units = list(units)
            else:
                L = np.array([self.log_embed(u) for u in units])
                G = L @ L.T
                U = pari.qflllgram(pari.matrix(r, r, [pari(float(x)) for x in G.flatten()]))
                U = [[int(U[i, j]) for j in range(r)] for i in range(r)]
                if abs(linalg.det(U)) != 1:
                    raise ArithmeticError('LLL transformation not unimodular')
                self._red_units = []
                for j in range(r):
                    u = self.one
                    for i in range(r):
                        if U[i][j]:
                            u = self.mul(u, self.power(units[i], U[i][j]))
                    self._red_units.append(self.as_int(u))
        return self._red_units

    def embed_hp(self, x):
        """sigma_i(x) as PARI reals at a precision adapted to the coordinates (heuristic
        use: estimates only)."""
        bits = max(max(Fraction(t).numerator.bit_length(), Fraction(t).denominator.bit_length()) for t in x)
        a = pari.lift(self.to_alg(x))
        prec = max(256, 2 * self.m * bits + 256)
        return [pari.subst(a, 'x', r) for r in self._hp_roots(prec)]

    def solve_embeddings(self, targets, bits):
        """Coordinates c (PARI reals, `bits` of precision) with sigma_i(sum c_j w_j) =
        targets[i] (targets: floats; heuristic use)."""
        roots = self._hp_roots(max(256, bits))
        m = self.m
        E = pari.matrix(m, m, [pari.subst(pari.lift(self.basis_alg[j]), 'x', roots[i])
                                for i in range(m) for j in range(m)])
        t = pari.Col([pari(Fraction(float(v)).numerator) / Fraction(float(v)).denominator for v in targets])
        return list(pari.matsolve(E, t))

    def embed_accurate(self, x):
        """sigma_i(x) as floats with ~1e-8 relative accuracy (high precision fallback)."""
        e = self.embed(x)
        scale = np.abs(self.E) @ np.abs(np.array([float(t) for t in x]))
        if np.all(np.abs(e) > 1e-6 * scale):
            return e
        a = pari.lift(self.to_alg(x))
        return np.array([float(pari.subst(a, 'x', r)) for r in self._hp_roots()])

    # ------------------------------------------------------------------
    # exact signs at real places

    def _isolating_intervals(self):
        """Disjoint rational intervals [a_i, b_i], each containing exactly the i-th real
        root of the defining polynomial (certified with Sturm sequences)."""
        if getattr(self, '_iso', None) is None:
            r = [Fraction(str(pari.bestappr(t, 2 ** 200))) for t in self._hp_roots(1024)]
            gaps = [r[i + 1] - r[i] for i in range(len(r) - 1)]
            d = (min(gaps) / 4) if gaps else Fraction(1)
            iso = []
            for t in r:
                a, b = t - d, t + d
                if int(pari.polsturm(self.pol, [pari(a), pari(b)])) != 1:
                    raise ArithmeticError('root isolation failed')
                iso.append([a, b])
            self._iso = iso
        return self._iso

    def _float_bound(self, w, i):
        """(value, error bound) of sigma_i(w) in double precision, or None if the
        coordinates are not exactly representable.  The bound 1e-12 * sum |w_k E_ik| is
        ~700 times the worst-case rounding error of an m-term dot product (m <= 10)
        with entries correct to one ulp."""
        if not all(isinstance(t, int) and abs(t) < 2 ** 50 for t in w):
            return None
        row = self.E[i]
        v = 0.0
        s = 0.0
        for t, e in zip(w, row):
            if t:
                v += t * e
                s += abs(t * e)
        return v, 1e-12 * s

    def sign_at(self, w, i):
        """Exact sign (-1, 0, 1) of sigma_i(w), places ordered by the real root."""
        if not any(w):
            return 0
        fb = self._float_bound(w, i)
        if fb is not None and abs(fb[0]) > fb[1]:
            return 1 if fb[0] > 0 else -1
        g = pari.lift(self.to_alg(w))
        a, b = self._iso_refine(i, g)
        ga = to_frac(pari.subst(g, 'x', pari(a)))
        return 1 if ga > 0 else -1

    def _iso_refine(self, i, g):
        """Shrink the isolating interval of root i until g has no root in it."""
        iso = self._isolating_intervals()
        a, b = iso[i]
        fa = to_frac(pari.subst(self.pol, 'x', pari(a)))
        while int(pari.polsturm(g, [pari(a), pari(b)])) != 0:
            c = (a + b) / 2
            fc = to_frac(pari.subst(self.pol, 'x', pari(c)))
            if (fc > 0) == (fa > 0):
                a, fa = c, fc
            else:
                b = c
        iso[i] = [a, b]
        return a, b

    def compare_abs(self, x, y, i):
        """Exact sign of |sigma_i(x)| - |sigma_i(y)| (0 iff x = +-y for x, y in O)."""
        return self.sign_at(tuple(a - b for a, b in zip(self.sq(x), self.sq(y))), i)

    def bnf_regulator(self):
        return float(member(self.bnf, 'reg'))

    def regulator(self):
        """Regulator of O^x (float): bnf regulator times [O_K^x : O^x]."""
        self.fundamental_units()
        return self.bnf_regulator() * self.unit_index

    def certify_bnf(self):
        return int(pari.bnfcertify(self.bnf)) == 1

    def unit_exponents_of(self, units):
        """Exponent vectors of units (coordinates in O) w.r.t. bnf's fundamental units of
        O_K (torsion dropped).  Raises if bnfisunit fails."""
        out = []
        for u in units:
            e = pari.bnfisunit(self.bnf, self.to_alg(u))
            if len(e) == 0:
                raise ArithmeticError('bnfisunit: not a unit of the bnf unit group')
            out.append([int(e[i]) for i in range(len(e) - 1)])
        return out

    def subgroup_index(self, units):
        """[O^x / +-1 : <units>] (0 if the units do not have full rank), computed with
        the units of bnfinit.  Exact once bnf's units are known to be fundamental, and
        index 1 proves that they are (Voronoi side-pairing argument)."""
        r = self.m - 1
        if r == 0:
            return 1
        if not units:
            return 0
        self.fundamental_units()
        ex = self.unit_exponents_of(units)
        M = pari.matrix(r, len(ex), [ex[j][i] for i in range(r) for j in range(len(ex))])
        H = pari.mathnf(M)
        if int(pari.matsize(H)[1]) != r:
            return 0
        sub = abs(int(pari.matdet(H)))
        # index of O^x in O_K^x
        return sub // self.unit_index if sub % self.unit_index == 0 else Fraction(sub, self.unit_index)

    def _order_units(self, fu_max):
        """Fundamental units of a non-maximal order from those of O_K (adapted from
        pav10/sails, sails/field.py: NumberFieldData._order_units)."""
        from itertools import product
        n = self.m
        nf = self.nf
        Mo, Minv = self._Mo, self._Minv
        f = int(pari.denominator(Minv))                      # exponent of O_K / O
        r = len(fu_max)
        bid = pari.idealstar(nf, f, 1)
        cyc = [int(t) for t in member(bid, 'cyc')]
        mm = len(cyc)
        if mm == 0:
            self.unit_index = 1
            self.unit_exponents = [[int(i == j) for i in range(r)] for j in range(r)]
            return list(fu_max)
        A = [[int(t) for t in pari.ideallog(nf, u, bid)] for u in fu_max]
        Lf = f * Minv
        U, V, D = pari.matsnf(Lf, 1)
        diag = [int(D[i, i]) for i in range(n)]
        Uinv = U ** -1
        logs = []
        for k in product(*[range(d) for d in diag]):
            xo = Uinv * pari.Col(list(k))
            x = Mo * xo
            if math.gcd(int(pari.nfeltnorm(nf, x)), f) != 1:
                continue
            logs.append([int(t) for t in pari.ideallog(nf, x, bid)])
        if logs:
            H = pari.mathnfmodid(pari.matconcat(pari.Vec([pari.Col(l) for l in logs])), cyc)
        else:
            H = pari.matdiagonal(cyc)
        cols = [pari.Col(a) for a in A] + [pari.Col([-H[i, j] for i in range(mm)]) for j in range(mm)]
        Ker = pari.matkerint(pari.matconcat(pari.Vec(cols)))
        k = int(pari.matsize(Ker)[1])
        E = pari.mathnf(pari.matrix(r, k, [Ker[i, j] for i in range(r) for j in range(k)]))
        if int(pari.matsize(E)[1]) != r:
            raise ArithmeticError('unit lattice of the order has wrong rank')
        self.unit_index = abs(int(pari.matdet(E)))
        self.unit_exponents = [[int(E[i, j]) for i in range(r)] for j in range(r)]
        return [pari.nfbasistoalg(nf, pari.nffactorback(nf, list(fu_max), [int(E[i, j]) for i in range(r)]))
                for j in range(r)]

    # ------------------------------------------------------------------

    def polredabs(self):
        return str(pari.polredabs(self.pol))

    def describe(self):
        return {'pol': str(self.pol), 'degree': self.m, 'signature': list(self.signature),
                'field_disc': self.field_disc, 'disc': self.disc, 'index': self.index,
                'basis_zk': self.basis_zk}
