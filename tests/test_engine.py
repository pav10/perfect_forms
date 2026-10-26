"""Regression benchmarks (plan, section 7).  Run:  python tests/test_engine.py [--slow]

Every check must pass before any run is reported.  Sources of the expected values:
* cubic rows: perfect-unary-cubic.tex (Obs. 5.4-5.6, the table after Obs. 5.6) and
  the certified GP census data/legacy/census_cubic_D20000_gp.jsonl (voronoi_unary.gp);
* real quadratic: Yasaki, Thm. (manuscript Sec. 1.2): h_d = 1 for squarefree
  d = (2k+1)^2 + 1;
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


def test_binary_additive():
    """Binary forms (additive normalisation): equivalence under random GL_2(O) transforms,
    Galois pairing, and the edge-balance certificate.  Class counts are regression values
    of this code (not yet compared with Ong / Leibak / Gunnells--Yasaki [REF?])."""
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
            if name == 'test_census_legacy':
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
