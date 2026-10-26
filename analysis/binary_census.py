"""Census of perfect binary forms (additive normalisation) over real quadratic fields.

For each squarefree d in the range: the GL_2(O_K)-classes of perfect binary forms
Tr(v^T A v) (engine/binary.py), with both certificates (edge balance, mass formula), the
number of classes up to Gal(K/Q) as well, the eutactic (= extreme) classes and the exact
maximum of gamma(A)^(2m) = 1/N(det A) (minimum 1).  Records are written as JSONL.

usage: python analysis/binary_census.py D_MIN D_MAX OUT.jsonl [--jobs J]
(d ranges over squarefree integers, K = Q(sqrt d); the counts for D <= 44 are checked
against Gunnells--Yasaki, Table 1, in tests/test_engine.py)
"""
import json
import os
import sys
import time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from engine.binary import BinaryVoronoi, edge_balance, mass_formula   # noqa: E402
from engine.field import Order                                         # noqa: E402
from engine.pari import pari                                           # noqa: E402
from engine.record import provenance                                   # noqa: E402


def galois_orbits(V):
    O, B = V.O, V.B
    gens = [s for s in pari.nfgaloisconj(O.nf) if s != pari('x')]
    if not gens:
        return len(V.classes)
    sigma = gens[0]

    def conj(z):
        return O.from_alg(pari.Mod(pari.subst(pari.lift(O.to_alg(z)), 'x', sigma), O.pol))

    seen, orbits = set(), 0
    n0 = len(V.classes)
    for i, P in enumerate(V.classes[:n0]):
        if i in seen:
            continue
        a, b, c = B.split(P.A)
        j, _ = V._identify(tuple(conj(a)) + tuple(conj(b)) + tuple(conj(c)))
        if len(V.classes) != n0:
            raise ArithmeticError('Galois conjugate of a perfect form not found among the classes')
        seen |= {i, j}
        orbits += 1
    return orbits


def one(d):
    t0 = time.time()
    pol = 'x^2-x-%d' % ((d - 1) // 4) if d % 4 == 1 else 'x^2-%d' % d
    O = Order(pol)
    V = BinaryVoronoi(O).run()
    eb = edge_balance(V)
    mass, table = mass_formula(V)
    S = V.summary()
    gmax = max(Fraction(s['gamma_pow']) for s in S)
    rec = {
        'd': d, 'D': O.disc, 'pol': O.polredabs(), 'class_number': O.class_number(),
        'n_perfect': len(V.classes),
        'n_perfect_mod_galois': galois_orbits(V),
        'n_eutactic': sum(s['eutaxy'] == 'eutactic' for s in S),
        'n_weakly_eutactic': sum(s['eutaxy'] == 'weakly' for s in S),
        'gamma_pow_max': str(gmax), 'gamma_max': float(gmax) ** 0.25,
        'cells_by_dim': {d_: len(v) for d_, v in sorted(table.items())},
        'certificates': {'edge_balance_ok': not eb, 'mass': str(mass), 'mass_ok': mass == 0},
        'classes': S,
        'stats': V.stats,
    }
    rec['status'] = 'CERTIFIED' if (not eb and mass == 0) else 'QUARANTINE'
    rec['provenance'] = provenance(time.time() - t0)
    return rec


def main():
    import argparse
    import multiprocessing as mp
    p = argparse.ArgumentParser()
    p.add_argument('dmin', type=int)
    p.add_argument('dmax', type=int)
    p.add_argument('out')
    p.add_argument('--jobs', type=int, default=1)
    a = p.parse_args()
    done = set()
    if os.path.exists(a.out):
        with open(a.out) as fh:
            done = {json.loads(l)['d'] for l in fh if l.strip()}
    ds = [d for d in range(max(2, a.dmin), a.dmax + 1) if int(pari.issquarefree(d)) and d not in done]
    with open(a.out, 'a') as fh, mp.get_context('fork').Pool(a.jobs) as pool:
        for rec in pool.imap_unordered(one, ds):
            fh.write(json.dumps(rec) + '\n')
            fh.flush()
            print(rec['d'], rec['n_perfect'], rec['n_perfect_mod_galois'], rec['status'],
                  rec['provenance']['seconds'], flush=True)


if __name__ == '__main__':
    main()
