"""Independent check of one field of a binary census (used for the Gunnells--Yasaki
discrepancy at D = 184, data/ERRATA.md).

Recomputes the GL_2(O_K)-classes of perfect binary forms over Q(sqrt d), groups them by
invariants that do not use qfisom (N(det A), number of minimal vectors, |Aut A|, the common
denominator of the trace forms, and the number of vectors of each trace form up to 4 times
its minimum), and runs qfisom in both directions on every pair that no invariant separates.

usage: python analysis/binary_verify.py d
"""
import collections
import itertools
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from engine.binary import BinaryVoronoi, edge_balance, mass_formula   # noqa: E402
from engine.field import Order                                         # noqa: E402
from engine.pari import pari                                           # noqa: E402


def counts(P, bound=4):
    out = []
    for F in P.forms:
        mn = int(pari.qfminim(F, 0, 0, 2)[1])
        out.append((mn, int(pari.qfminim(F, bound * mn, 0, 2)[0])))
    return tuple(out)


def main(d):
    d = int(d)
    pol = 'x^2-x-%d' % ((d - 1) // 4) if d % 4 == 1 else 'x^2-%d' % d
    t = time.time()
    V = BinaryVoronoi(Order(pol)).run()
    print('d = %d: %d classes (%.0fs); edge balance %s; mass %s' % (
        d, len(V.classes), time.time() - t, not edge_balance(V), mass_formula(V)[0]), flush=True)
    groups = collections.defaultdict(list)
    for j, P in enumerate(V.classes):
        groups[(P.normdet, len(P.minvecs), P.aut_order, P.den, counts(P))].append(j)
    print('invariant groups: %d; sizes %s' % (
        len(groups), dict(sorted(collections.Counter(len(g) for g in groups.values()).items()))), flush=True)
    pairs = iso = 0
    for g in groups.values():
        for i, j in itertools.combinations(g, 2):
            pairs += 1
            a, b = V.classes[i], V.classes[j]
            if pari.qfisom(pari(a.forms), pari(b.forms)) != 0 or pari.qfisom(pari(b.forms), pari(a.forms)) != 0:
                iso += 1
                print('isometric pair', i, j, flush=True)
    print('pairs not separated by invariants: %d; isometric: %d' % (pairs, iso))


if __name__ == '__main__':
    main(*sys.argv[1:])
