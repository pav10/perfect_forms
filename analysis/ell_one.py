"""ell_O = 1 test (analysis/FINDINGS.md section 9): all relative minima are units iff all
Minkowski--Voronoi neighbours of 1 are units (unit translations map neighbours to
neighbours, and the adjacency graph is connected).  Expensive for large D (the neighbours of
1 lie in long boxes); FINDINGS uses it only on the fields that pass the necessary condition
"at most 2 sail facets modulo units" (manuscript Prop. 7.1 with V = 1).

usage: python analysis/ell_one.py FIELDS.jsonl.xz DMAX OUT.jsonl JOBS
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


if __name__ == '__main__':
    rows = [json.loads(l) for l in lzma.open(sys.argv[1], 'rt')]
    rows = [r for r in rows if r['D'] <= int(sys.argv[2])]
    with open(sys.argv[3], 'a') as fh, mp.get_context('fork').Pool(int(sys.argv[4]), maxtasksperchild=200) as pool:
        for rec in pool.imap_unordered(one, rows, chunksize=8):
            fh.write(json.dumps(rec) + '\n')
            fh.flush()
