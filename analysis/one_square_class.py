"""Integral perfect forms with unit minimal vectors (cubic census check of
docs/THEORY.md Sec. 5, Lemma 5.2).

Lemma 5.2 (PROVED): if alpha is totally positive, lies in the codifferent O^v, has
min Tr(alpha x^2) = 1 on O \\ 0 and three minimal vectors u_1, u_2, u_3 that are units
with distinct squares, then the Gram matrix (Tr(alpha u_j u_k)) is the identity,
{u_1, u_2, u_3} is a Z-basis of O and N(alpha) D = 1, so alpha generates O^v.  Two such
forms differ by a totally positive unit, hence there are at most [U+ : U^2] classes of
them modulo units (the census counts perfect forms modulo O^x acting by alpha -> alpha u^2).

For every CERTIFIED record this script rescales each perfect form to minimum 1, tests
alpha in O^v (Tr(alpha w_j) in Z), and for the integral ones whose minimal vectors are
units checks the three conclusions and the bound on the number of classes.

usage: python analysis/one_square_class.py CENSUS.jsonl[.xz]
"""
import json
import lzma
import os
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from engine.field import Order        # noqa: E402
from engine import linalg             # noqa: E402


def records(path):
    op = lzma.open if path.endswith('.xz') else open
    with op(path, 'rt') as fh:
        for line in fh:
            if line.strip():
                yield json.loads(line)


def tp_mod_squares(O):
    """[U+ : U^2] = 8 / #(sign vectors of units), for m = 3."""
    units = O.fundamental_units()
    signs = set()
    for e in range(8):
        w = O.one
        for i, u in enumerate(units):
            if e >> i & 1:
                w = O.mul(w, u)
        s = tuple(O.sign_at(w, i) for i in range(O.m))
        signs.add(s)
        signs.add(tuple(-t for t in s))
    return 8 // len(signs)


def check(r):
    O = Order(r['pol'])
    D = O.disc
    out = []
    for f in r['forms']:
        alpha = [Fraction(t) for t in f['alpha_coords']]
        x0 = tuple(f['minvecs'][0])
        mn = O.trace_pair(alpha, O.sq(x0))
        alpha = [t / mn for t in alpha]
        integral = all(O.trace_pair(alpha, w).denominator == 1 for w in
                       [tuple(int(i == j) for i in range(O.m)) for j in range(O.m)])
        units = all(abs(n) == 1 for n in f['minvec_norms'])
        if not (integral and units):
            continue
        G = O.gram(alpha)
        V = [tuple(v) for v in f['minvecs']]
        GV = [[sum(V[a][j] * G[j][k] * V[b][k] for j in range(O.m) for k in range(O.m))
               for b in range(len(V))] for a in range(len(V))]
        identity = len(V) == 3 and all(GV[a][b] == (a == b) for a in range(3) for b in range(3))
        basis = len(V) == 3 and abs(linalg.det([list(v) for v in V])) == 1
        nd = abs(O.norm(alpha)) * D
        out.append({'identity': identity, 'basis': basis, 'N_alpha_D': str(nd)})
    return out, (tp_mod_squares(O) if out else None)


def main(path):
    n_fields = n_forms = 0
    bad = []
    fields = []
    for r in records(path):
        if r['status'] != 'CERTIFIED':
            continue
        n_fields += 1
        res, idx = check(r)
        if not res:
            continue
        n_forms += len(res)
        ok = all(t['identity'] and t['basis'] and t['N_alpha_D'] == '1' for t in res) and len(res) <= idx
        fields.append((r['disc'], len(res), idx))
        if not ok:
            bad.append((r['label'], res, idx))
    print('%d certified fields; %d integral perfect forms with unit minimal vectors, in %d fields' % (
        n_fields, n_forms, len(fields)))
    print('Gram = identity, minimal vectors a Z-basis, N(alpha) D = 1 and #classes <= [U+ : U^2]: %s' % (
        'all' if not bad else 'FAILS for %s' % bad))
    print('(D, #forms, [U+ : U^2]) for the first 30: %s' % fields[:30])
    by = {}
    for _, k, idx in fields:
        by[(k, idx)] = by.get((k, idx), 0) + 1
    print('fields by (#forms, [U+ : U^2]): %s' % dict(sorted(by.items())))


if __name__ == '__main__':
    main(*sys.argv[1:])
