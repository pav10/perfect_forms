"""Regression benchmarks (plan, section 7).  Run:  python tests/test_engine.py [--slow]

Every check must pass before any run is reported.  Sources of the expected values:
* cubic rows: perfect-unary-cubic.tex (Obs. 5.4-5.6, the table after Obs. 5.6) and
  the certified GP census data/legacy/census_cubic_D20000_gp.jsonl (voronoi_unary.gp);
* real quadratic: Yasaki, JTNB 25 (2013): h_d = 1 for squarefree d = (2k+1)^2 + 1
  (proof of Cor. 5.1: l_d = 1, then Thm. 5.1), the examples Q(sqrt 19) and Q(sqrt 23) of
  Sec. 5.4 and Cor. 3.1 (two or three minimal vectors up to sign);
* regulator of 3.3.49.1 from the LMFDB.
"""
import itertools
import json
import os
import random
import sys
import time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from engine import extreme, qfmin                      # noqa: E402
from engine.certify import certificates                # noqa: E402
from engine.field import Order                         # noqa: E402
from engine.pari import pari                           # noqa: E402
from engine.relmin import RelativeMinima               # noqa: E402
from engine.unary import Voronoi                       # noqa: E402


def run(pol):
    O = Order(pol)
    V = Voronoi(O).run()
    c = certificates(O, V)
    assert c['certified'], (pol, c)
    f = V.f_vector()
    return O, V, f


def test_enumerator_bruteforce():
    random.seed(1)
    for _ in range(200):
        n = random.choice([2, 3, 4])
        A = [[random.randint(-3, 3) for _ in range(n)] for _ in range(n)]
        G = [[Fraction(sum(A[k][i] * A[k][j] for k in range(n)) + (n if i == j else 0), random.choice([1, 2, 3]))
              for j in range(n)] for i in range(n)]
        G = [[G[min(i, j)][max(i, j)] for j in range(n)] for i in range(n)]
        if not qfmin.is_positive_definite(G):
            continue
        B = Fraction(random.randint(5, 30), random.randint(1, 3))
        R = 8
        bf = []
        for x in itertools.product(range(-R, R + 1), repeat=n):
            if any(x) and next(t for t in x if t) > 0:
                q = qfmin.qeval(G, x)
                if q <= B:
                    bf.append((q, x))
        sv = qfmin.short_vectors(G, B)
        if all(max(abs(t) for t in x) < R for _, x in sv):
            assert sorted(sv) == sorted(bf), (G, B)
            if bf:
                m, mv = qfmin.minimum(G)
                assert m == min(q for q, _ in bf)
                assert sorted(mv) == sorted(x for q, x in bf if q == m)


CUBIC = {   # D: (polynomial, mu, n, degree multiset or None)
    49: ('x^3-x^2-2*x+1', 1, 2, {3: 2}),
    257: ('x^3-x^2-4*x+3', 2, 4, {3: 4}),
    404: ('x^3-x^2-5*x-1', 4, 7, {3: 6, 4: 1}),
    473: ('x^3-5*x-1', 2, 4, None),
    1708: ('x^3-x^2-8*x-2', 12, 21, None),
    6856: ('x^3-x^2-19*x-1', 23, 46, None),
    6997: ('x^3-x^2-13*x+14', 13, 26, None),
    16689: ('x^3-x^2-26*x+24', 77, 152, None),
}


def test_cubic_benchmarks():
    for D, (pol, mu, n, degs) in CUBIC.items():
        O, V, f = run(pol)
        assert O.disc == D
        assert (f[0], f[2]) == (mu, n), (D, f)
        if degs:
            from collections import Counter
            assert dict(Counter(V.degrees())) == degs, (D, V.degrees())


def test_ennola_shanks():
    discs = [49, 257, 697, 1489, 2777, 4729, 7537, 11417, 16609, 23377, 32009]
    for a, D in zip(range(2, 13), discs):
        O, V, f = run('x^3+(%d)*x^2-(%d)*x-1' % (a - 1, a))
        assert O.disc == D
        assert (f[0], f[2]) == ((2, 4) if a == 3 else (1, 2)), (a, f)
    for a in range(-1, 12):
        O, V, f = run('x^3-(%d)*x^2-(%d)*x-1' % (a, a + 3))
        assert (f[0], f[2]) == (1, 2), (a, f)


def test_large_cubic():
    O, V, f = run('x^3-x^2-221*x+29')           # manuscript table: D = 676829, (mu, n) = (136, 272)
    assert (O.disc, f[0], f[2]) == (676829, 136, 272)
    assert set(V.degrees()) == {3}


def test_yasaki_quadratic_family():
    """h_d = 1 for squarefree d = (2k+1)^2 + 1 (Yasaki); for d = 2 also n = mu = 1."""
    for k in range(0, 12):
        d = (2 * k + 1) ** 2 + 1
        if not int(pari.issquarefree(d)):
            continue
        O, V, f = run('x^2-%d' % d)
        assert f == [1, 1], (d, f)


def test_yasaki_examples():
    """Yasaki, JTNB 25 (2013), Sec. 5.4: Q(sqrt 19) has four classes and Q(sqrt 23) two, with the
    minimal vectors listed there (a + b omega, omega = sqrt d); each set is the set of minimal
    vectors (up to sign) of exactly one class, and the classes found are pairwise distinct."""
    from engine import linalg
    from engine.unary import PerfectForm, normalize
    examples = {19: [[(1, 0), (4, 1)], [(4, 1), (9, 2), (13, 3)], [(13, 3), (48, 11), (61, 14)],
                     [(61, 14), (170, 39)]],
                23: [[(1, 0), (5, 1)], [(5, 1), (24, 5)]]}
    for d, sets in examples.items():
        O, V, f = run('x^2-%d' % d)
        assert f[1] == len(sets), (d, f)
        hit = []
        for S in sets:
            S = sorted(normalize(x) for x in S)
            P = PerfectForm(O, linalg.solve([O.form_row(x) for x in S[:2]], [1, 1]))
            assert sorted(P.minvecs) == S, (d, S, P.minvecs)
            js = [j for j in V.by_key.get(O.ideal_key(P.alpha), ())
                  if V._find_eta(P.alpha, P.minvecs, V.classes[j]) is not None]
            assert len(js) == 1, (d, S, js)
            hit += js
        assert sorted(hit) == list(range(len(sets))), (d, hit)


def test_regulator_lmfdb():
    O, V, f = run('x^3-x^2-2*x+1')
    assert abs(O.regulator() - 0.525454682123) < 1e-10


def test_higher_degree_euler():
    for pol in ['x^4-x^3-3*x^2+x+1', 'x^4-5*x^2+5', 'x^4-x^3-4*x^2+4*x+1', 'x^5-x^4-4*x^3+3*x^2+3*x-1']:
        run(pol)


def test_relmin_crosscheck():
    for pol in ['x^2-46', 'x^3-x^2-5*x-1', 'x^3-x^2-8*x-2', 'x^3-x^2-17*x+18', 'x^4-x^3-3*x^2+x+1']:
        O, V, f = run(pol)
        R = RelativeMinima(O).run()
        keys = set(R.by_key)
        assert all(k in keys for k in V.extreme_square_classes()), pol
        assert O.subgroup_index(R.side_pairing_units()) == 1
        assert f[0] <= R.count()
        for x in V.extreme_square_classes().values():
            assert R.is_relative_minimum(x)


def test_quadratic_relmin_continued_fraction():
    from engine.relmin import quadratic_relmin_count
    for D in range(5, 200):
        if not int(pari.isfundamental(D)):
            continue
        O = Order('x^2-x-%d' % ((D - 1) // 4) if D % 4 == 1 else 'x^2-%d' % (D // 4))
        assert RelativeMinima(O).run().count() == quadratic_relmin_count(D), D


def test_gamma():
    O, V, f = run('x^2-3')                 # Q(sqrt 3) attains Leibak's bound: gamma^2 = (4/3) 12 = 16
    s = extreme.summary(O, V)
    assert s['gamma_pow'] == 16 and s['leibak_ratio'] == 1 and s['argmax_eutactic']
    O, V, f = run('x^3-x^2-8*x-2')         # D = 1708: the A_3 case gamma^3 = 2 D
    s = extreme.summary(O, V)
    assert s['gamma_pow'] == 2 * 1708 and s['argmax_eutactic']


def test_leibak_equality():
    """docs/THEORY.md 4: Leibak's bound is attained iff (m = 2) O has a basis u, v with
    N(u) + N(v) + N(u - v) = 0, (m = 3) some beta >> 0 with (beta) = p^2 d^-1, 2O = p q^2, has
    Tr(beta x^2) even; in both cases 2 ramifies.  Real quadratic D < 1200 and five cubic fields
    (attained: 316, 940; 2O = p q^2 but not attained: 568; 2O = q^3: 148; 2 unramified: 49)."""
    sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'analysis'))
    import leibak_equality as le
    att = 0
    for D in range(5, 1200):
        if not int(pari.isfundamental(D)) or D < 0:
            continue
        O, V, f = run(_quadratic_pol(D))
        a = extreme.summary(O, V)['leibak_ratio'] == 1
        assert a == le.quadratic(D), D
        att += a
    assert att == 19, att
    for pol, a in (('x^3-x^2-4*x+2', True), ('x^3-7*x-4', True), ('x^3-x^2-6*x-2', False),
                   ('x^3-x^2-3*x+1', False), ('x^3-x^2-2*x+1', False)):
        O, V, f = run(pol)
        assert (extreme.summary(O, V)['leibak_ratio'] == 1) == a == le.cubic(pol), (pol, O.disc)


def test_one_square_class():
    """docs/THEORY.md 7 (manuscript Sec. 9).  Lemma 7.2: an integral perfect form with unit
    minimal vectors is I_3 and generates the codifferent (N(alpha) D = 1); one such class in
    the simplest cubic fields with O = Z[theta_t] (D = 49, 81, 169), none in 3.3.148.1.
    Lemma 7.5: minimal triples are bases or degenerate (checked inside counts); degenerate
    classes: D = 49 and 169 (ell = 1) one, D = 81 none, D = 148 two.  Lemma 7.6: for D = 169,
    eps = N_k(1) and 1 - eps are units and eps is a root of P_1(x - 1)."""
    sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'analysis'))
    import one_square_class as osc
    import minkowski_voronoi as mv
    for pol, k in (('x^3-x^2-2*x+1', 1), ('x^3-3*x-1', 1), ('x^3-x^2-4*x-1', 1), ('x^3-x^2-3*x+1', 0)):
        O, V, f = run(pol)
        res = osc.integral_unit_forms(O, [(P.alpha, P.minvecs) for P in V.classes])
        assert len(res) == k, (pol, res)
        assert all(r['identity'] and r['basis'] and r['N_alpha_D'] == '1' for r in res), (pol, res)
    for pol, ell, dg in (('x^3-x^2-2*x+1', 1, 1), ('x^3-x^2-4*x-1', 1, 1), ('x^3-3*x-1', 2, 0),
                         ('x^3-x^2-3*x+1', 3, 2)):
        c = mv.counts(pol)
        assert (c['ell'], c['V'], c['V_degenerate'], c['euler']) == (ell, 2 * ell, dg, 0), (pol, c)
    O = Order('x^3-x^2-4*x-1')
    M = RelativeMinima(O)
    for k in range(3):
        e = M.single_neighbour(O.one, k)
        if O.embed_accurate(e)[k] < 0:
            e = tuple(-t for t in e)
        one_minus = tuple(-e[i] + (i == 0) for i in range(3))
        assert abs(O.norm(e)) == 1 and abs(O.norm(one_minus)) == 1
        assert str(pari.charpoly(O.to_alg(e))) == 'x^3 - 4*x^2 + x + 1'


def test_binary_additive():
    """Binary forms (additive normalisation): equivalence under random GL_2(O) transforms,
    Galois pairing, and the edge-balance certificate.  The class counts agree with
    Gunnells--Yasaki Table 1 (test_gunnells_yasaki_binary)."""
    from engine.binary import BinaryVoronoi, edge_balance, mass_formula
    expected = {'x^2-x-1': 2, 'x^2-2': 2, 'x^2-3': 3, 'x^2-x-3': 9}
    for pol, k in expected.items():
        O = Order(pol)
        V = BinaryVoronoi(O).run()
        assert len(V.classes) == k, (pol, len(V.classes))
        assert edge_balance(V) == [], pol
        assert mass_formula(V)[0] == 0, pol
        assert all(V.eutaxy(P) != 'no' for P in V.classes if P.normdet == min(Q.normdet for Q in V.classes))
    random.seed(5)
    O = Order('x^2-x-3')
    V = BinaryVoronoi(O).run()
    B = V.B
    for i, P in enumerate(V.classes):
        while True:
            M = [[tuple(random.randint(-2, 2) for _ in range(2)) for _ in range(2)] for _ in range(2)]
            d = tuple(s - t for s, t in zip(O.mul(M[0][0], M[1][1]), O.mul(M[0][1], M[1][0])))
            if abs(O.norm(d)) == 1:
                break
        a, b, c = B.split(P.A)
        Aa = [[a, b], [b, c]]

        def entry(i_, j_):
            out = tuple([Fraction(0)] * 2)
            for k in range(2):
                for l in range(2):
                    out = tuple(s + t for s, t in zip(out, O.mul(O.mul(M[k][i_], Aa[k][l]), M[l][j_])))
            return out
        A2 = entry(0, 0) + entry(0, 1) + entry(1, 1)
        n0 = len(V.classes)
        j, _ = V._identify(A2)
        assert j == i and len(V.classes) == n0


# Gunnells--Yasaki, Perfect forms over totally real number fields (arXiv:0908.3130), Table 1:
# number of GL_2(O)-classes of perfect binary forms over Q(sqrt d), by discriminant D, as printed.
# The entry for D = 184 (4306) disagrees with this engine's certified 4308 (data/ERRATA.md); the
# tests below use only D <= 44.
GY_TABLE = {5: 2, 8: 2, 12: 3, 13: 9, 17: 34, 21: 16, 24: 22, 28: 30, 29: 39, 33: 149, 37: 97, 40: 52,
            41: 144, 44: 66, 53: 93, 56: 144, 57: 515, 60: 89, 61: 206, 65: 301, 76: 404, 88: 645, 92: 343,
            104: 326, 120: 524, 124: 1473, 136: 1135, 140: 487, 152: 995, 156: 941, 168: 918, 172: 2683,
            184: 4306, 188: 1358, 204: 1992, 220: 2397, 232: 3042, 236: 2977, 248: 2360, 264: 3139}


def _quadratic_pol(D):
    return 'x^2-x-%d' % ((D - 1) // 4) if D % 4 == 1 else 'x^2-%d' % (D // 4)


def test_gunnells_yasaki_binary(slow=False):
    """Binary perfect forms (additive) vs Gunnells--Yasaki Table 1; the slow run covers
    every D <= 44, the fast one the smallest fields."""
    from engine.binary import BinaryVoronoi, mass_formula
    Ds = [D for D in GY_TABLE if D <= (44 if slow else 13)]
    for D in Ds:
        V = BinaryVoronoi(Order(_quadratic_pol(D))).run()
        assert len(V.classes) == GY_TABLE[D], (D, len(V.classes))
        assert mass_formula(V)[0] == 0, D


def test_leibak_porter_ling_and_yasaki():
    """Leibak--Porter--Ling, JTNB 35 (2023), Thm 1.2: for real quadratic K = Q(sqrt d),
    n_K = 1 iff d = m^2+1 (m odd), m^2-1 (m even), m^2+4 (m odd) or m^2-4 (m > 3 odd).
    Yasaki, JTNB 25 (2013), Thm 5.1: l_d/2 <= n_K <= l_d with l_d the period length of the
    continued fraction of c_d (sqrt d, or (-1 + sqrt d)/2 for d = 1 mod 4), and Cor. 3.1."""
    import math
    from engine.relmin import quadratic_relmin_count

    def lpl(d):
        for m in range(1, math.isqrt(d) + 3):
            if (d == m * m + 1 and m % 2) or (d == m * m - 1 and m % 2 == 0) or \
                    (d == m * m + 4 and m % 2) or (d == m * m - 4 and m % 2 and m > 3):
                return True
        return False
    for d in range(2, 400):
        if not int(pari.issquarefree(d)):
            continue
        D = d if d % 4 == 1 else 4 * d
        O = Order(_quadratic_pol(D))
        V = Voronoi(O).run()
        n = len(V.classes)
        assert (n == 1) == lpl(d), (d, n)
        assert all(len(P.minvecs) in (2, 3) for P in V.classes), d       # Yasaki, Cor. 3.1
        ell = quadratic_relmin_count(D)
        assert ell <= 2 * n and n <= ell, (d, n, ell)


def test_simplest_cubics_lpl():
    """Leibak--Porter--Ling, Thm 1.5: n_K = 2 for x^3 - t x^2 - (t+3) x - 1 (t >= 0) when O_K is
    monogenic (Z[theta] = O_K)."""
    for t in range(0, 12):
        pol = 'x^3-(%d)*x^2-(%d)*x-1' % (t, t + 3)
        if int(pari.nfdisc(pari(pol))) != int(pari.poldisc(pari(pol))):
            continue
        O, V, f = run(pol)
        assert f[2] == 2, (t, f)


def test_komatsu_watanabe():
    """Komatsu--Watanabe, IJNT 10 (2014), Thms 3.1-3.2: for K = Q(v), v^3 + l v - 1 = 0, l >= 2 and
    4 l^3 + 27 squarefree, n_K = 1, the perfect form is a_0 = ((v^-1 - 1)/l, (1 - v^2)/(2l)) up to
    units, and S(a_0) = {+-1, +-v}."""
    from engine.complexcubic import ComplexCubicVoronoi
    for l in (2, 4, 5, 7, 10):
        V = ComplexCubicVoronoi('x^3+%d*x-1' % l).run()
        assert len(V.classes) == 1 and len(V.edges) == 2 and V.unit_index() == 1, l
        K = V.K
        A, S, _ = V.classes[0]
        v = K.from_alg(pari.Mod(pari('x'), K.pol))
        assert sorted(S) == sorted([K.one, tuple(int(c) for c in v)]), (l, S)
        a1 = pari('(v,l)->lift((1/v - 1)/l)')(pari.Mod(pari('x'), K.pol), l)
        a2 = pari('(v,l)->lift((1 - v^2)/(2*l))')(pari.Mod(pari('x'), K.pol), l)
        assert pari.lift(K.to_alg(A[0])) == a1 and pari.lift(K.to_alg(A[1])) == a2, l


def test_census_legacy(slow=False):
    """All 832 totally real cubic fields with D <= 20000 against the certified GP census."""
    if not slow:
        return
    leg = {}
    with open(os.path.join(os.path.dirname(HERE), 'data', 'legacy', 'census_cubic_D20000_gp.jsonl')) as fh:
        for line in fh:
            r = json.loads(line)
            leg[str(pari.polredabs(pari(r['pol'])))] = r
    L = list(pari('nflist("C3",[1,20000],0)')) + list(pari('nflist("S3",[1,20000],0)'))
    assert len(L) == 832 == len(leg)
    mu1 = 0
    for p in L:
        pol = str(pari.polredabs(p))
        O, V, f = run(pol)
        r = leg[pol]
        assert (f[2], f[1], f[0], V.degrees()) == (r['n'], r['E'], r['mu'], sorted(r['degs'])), pol
        mu1 += f[0] == 1
    assert mu1 == 27


def main():
    slow = '--slow' in sys.argv
    tests = [(k, v) for k, v in sorted(globals().items()) if k.startswith('test_') and callable(v)]
    failed = 0
    for name, fn in tests:
        t = time.time()
        try:
            if name in ('test_census_legacy', 'test_gunnells_yasaki_binary'):
                fn(slow=slow)
            else:
                fn()
            print('PASS %-32s %6.1fs' % (name, time.time() - t))
        except Exception as e:      # noqa: BLE001
            failed += 1
            print('FAIL %-32s %r' % (name, e))
    print('%d tests, %d failed' % (len(tests), failed))
    sys.exit(1 if failed else 0)


if __name__ == '__main__':
    main()
