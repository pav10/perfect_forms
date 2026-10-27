"""The Minkowski--Voronoi complex of O modulo units (Karpenkov--Ustinov, JNT 176 (2017), Sec. 2).

Its 2-faces are the relative minima, its edges the adjacent pairs, and its vertices the
minimal triples {x, y, z}: the box Pi(x, y, z) = {w : |w_i| <= max(|x_i|, |y_i|, |z_i|)}
contains no lattice point besides 0, +-x, +-y, +-z, and each of x, y, z attains the
maximum in some coordinate (KU, Def. 2.4).  It is the "tropical" counterpart of Voronoi's
complex of perfect forms (2-faces = extreme squares, vertices = perfect forms): a perfect
form is a linear functional sum a_i X_i supporting {(sigma_i(y)^2)} at >= 3 points, a
minimal triple a box max_i X_i/c_i supporting {(|sigma_i(y)|)} at 3 points.

For each class of relative minima x (relmin.py), the triples through x are counted among
the lattice points of the box B* that relmin.py enumerates to find the neighbours of x
(Pi(x, y, z) lies in B* for neighbours y, z).  Counts modulo units: F = ell, E = (sum of
the neighbour counts)/2, V = (sum of the triple counts)/3.

usage: python analysis/minkowski_voronoi.py CENSUS.jsonl[.xz] [LIMIT (0 = all)] [OUT.jsonl] [JOBS]
"""
import json
import lzma
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from engine.field import Order                        # noqa: E402
from engine.relmin import RelativeMinima, _Pts, _neg  # noqa: E402


def triples_at(M, x):
    """Number of minimal triples {x, y, z} (y, z neighbours of x up to sign)."""
    O = M.O
    Ns = [M.single_neighbour(x, k) for k in range(O.m)]
    ab = [abs(O.embed_accurate(Ns[k])[k]) for k in range(O.m)]
    beta = M._weight([1 / b ** 2 for b in ab])
    R = sum(M._upper(beta, Ns[k], k) for k in range(O.m))
    pts = M._enum(beta, R)
    nbs = M.neighbours_of(x)
    S = {x, _neg(x)} | set(nbs) | {_neg(y) for y in nbs}
    cand = [p for p in pts if p not in S]
    P = _Pts(O, [x] + nbs + cand)          # 0 = x, 1..k = neighbours, then the rest
    k = len(nbs)
    m = O.m

    def arg_max(idx, i):
        b = idx[0]
        for a in idx[1:]:
            if P.less(b, a, i):
                b = a
        return b

    count = 0
    for a in range(1, k + 1):
        for b in range(a + 1, k + 1):
            idx = (0, a, b)
            top = [arg_max(idx, i) for i in range(m)]
            if set(top) != set(idx):            # each point attains a maximum
                continue
            empty = True
            for z in list(range(1, k + 1)) + list(range(k + 1, len(P.pts))):
                if z in idx:
                    continue
                if all(not P.less(top[i], z, i) for i in range(m)):
                    empty = False
                    break
            count += empty
    return count, k


def counts(pol):
    O = Order(pol)
    M = RelativeMinima(O).run()
    tri, deg = 0, 0
    for x in M.reps:
        t, k = triples_at(M, x)
        tri += t
        deg += k
    ell = M.count()
    if deg % 2 or tri % 3:
        raise ArithmeticError('incidence counts not divisible: %d, %d' % (deg, tri))
    V, E, F = tri // 3, deg // 2, ell
    return {'ell': ell, 'V': V, 'E': E, 'F': F, 'euler': V - E + F}


def _one(r):
    t = time.time()
    c = counts(r['pol'])
    c.update(label=r['label'], D=r['disc'], n=r['n'], mu=r['mu'], seconds=round(time.time() - t, 2))
    c['ok'] = c['ell'] == r.get('ell', c['ell']) and c['euler'] == 0 and c['V'] == 2 * c['ell']
    return c


def main(path, limit=0, out=None, jobs=1):
    import multiprocessing as mp
    op = lzma.open if path.endswith('.xz') else open
    with op(path, 'rt') as fh:
        recs = [json.loads(l) for l in fh if l.strip()]
    recs = [r for r in recs if r['status'] == 'CERTIFIED']
    if int(limit):
        recs = recs[:int(limit)]
    fh = open(out, 'w') if out else None
    bad = 0
    with mp.get_context('fork').Pool(int(jobs)) as pool:
        for c in pool.imap_unordered(_one, recs, chunksize=2):
            if not c['ok']:
                bad += 1
                print('!!', json.dumps(c), flush=True)
            if fh:
                fh.write(json.dumps(c) + '\n')
                fh.flush()
    print('%d fields; V = 2 ell, E = 3 ell, chi = 0 fails in %d' % (len(recs), bad))


if __name__ == '__main__':
    main(*sys.argv[1:])
