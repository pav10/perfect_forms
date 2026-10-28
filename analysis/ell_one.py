"""ell_O = 1 test (analysis/FINDINGS.md section 9): all relative minima are units iff all
Minkowski--Voronoi neighbours of 1 are units (unit translations map neighbours to
neighbours, and the adjacency graph is connected).  Expensive for large D (the neighbours of
1 lie in long boxes); FINDINGS uses it only on the fields that pass the necessary condition
"at most 2 sail facets modulo units" (manuscript Prop. 7.1 with V = 1).

usage: python analysis/ell_one.py FIELDS.jsonl.xz DMAX OUT.jsonl JOBS
       python analysis/ell_one.py --structure CANDIDATES.json
The second form adds to each candidate record the Minkowski--Voronoi data (ell, number of
minimal triple classes V and degenerate ones, analysis/minkowski_voronoi.py) and, for each
direction k, the neighbour eps = N_k(1) normalised by sigma_k(eps) > 0: the signs of its
other two conjugates, whether eps and 1 - eps are units, whether 1 - eps is adjacent to 1,
and the characteristic polynomial of eps (docs/THEORY.md Lemma 7.6: P_t(x - 1)).
"""
import json
import lzma
import multiprocessing as mp
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from engine.field import Order              # noqa: E402
from engine.relmin import RelativeMinima    # noqa: E402


def one(row):
    t = time.time()
    O = Order(row['pol'])
    M = RelativeMinima(O)
    nb = M.neighbours_of(O.one)
    units = all(abs(O.norm(y)) == 1 for y in nb)
    return {'label': row['label'], 'D': row['D'], 'pol': row['pol'], 'deg1': len(nb), 'ell1': units,
            'seconds': round(time.time() - t, 3)}


def structure(row):
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from minkowski_voronoi import counts
    from engine.pari import pari
    O = Order(row['pol'])
    M = RelativeMinima(O)
    nb = M.neighbours_of(O.one)
    nbs = set(nb) | {tuple(-t for t in y) for y in nb}
    dirs = []
    for k in range(O.m):
        e = M.single_neighbour(O.one, k)
        if O.embed_accurate(e)[k] < 0:
            e = tuple(-t for t in e)
        om = tuple((i == 0) - e[i] for i in range(O.m))
        dirs.append({'signs': ''.join('+' if O.sign_at(e, j) > 0 else '-' for j in range(O.m) if j != k),
                     'eps_unit': abs(O.norm(e)) == 1, 'one_minus_eps_unit': abs(O.norm(om)) == 1,
                     'one_minus_eps_adjacent': om in nbs, 'charpoly': str(pari.charpoly(O.to_alg(e)))})
    c = counts(row['pol'])
    return dict(row, mv={'ell': c['ell'], 'V': c['V'], 'V_degenerate': c['V_degenerate']}, directions=dirs)


if __name__ == '__main__' and sys.argv[1] == '--structure':
    rows = json.load(open(sys.argv[2]))
    rows = [structure(r) for r in rows]
    with open(sys.argv[2], 'w') as fh:
        json.dump(rows, fh, indent=1)
    print('%d candidates; ell = 1: %d' % (len(rows), sum(r['ell1'] for r in rows)))
    for r in rows:
        print(r['D'], r['ell1'], r['mv'], [(d['signs'], d['eps_unit'], d['one_minus_eps_unit'],
                                          d['one_minus_eps_adjacent'], d['charpoly']) for d in r['directions']][:1])
elif __name__ == '__main__':
    rows = [json.loads(l) for l in lzma.open(sys.argv[1], 'rt')]
    rows = [r for r in rows if r['D'] <= int(sys.argv[2])]
    with open(sys.argv[3], 'a') as fh, mp.get_context('fork').Pool(int(sys.argv[4]), maxtasksperchild=200) as pool:
        for rec in pool.imap_unordered(one, rows, chunksize=8):
            fh.write(json.dumps(rec) + '\n')
            fh.flush()
