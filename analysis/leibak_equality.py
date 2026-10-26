"""When is Leibak's bound gamma_K^m <= gamma_m^m D attained?  (docs/THEORY.md, section 4)

PROVED: for 2 <= m <= 8, gamma_K^m = gamma_m^m D iff some beta >> 0 in the codifferent d^-1
has N(beta d) = det(Lambda_m) (3 for m = 2, 4 for m = 3) and Tr(beta x^2) in 2Z for all
x in O; then 2 ramifies in K.  For m = 2 this is: O has a Z-basis (u, v) with
N(u) + N(v) + N(u - v) = 0; for m = 3: 2O = p q^2 and (beta) = p^2 d^-1.

This script compares the criterion with the exact `leibak_ratio` of census records:
* degree 2: exact and unconditional (qfbsolve on the indefinite forms [A, 2(A+C), C] with
  A^2 + AC + C^2 = D/4);
* degree 3: the principal ideal test uses bnfinit (GRH-conditional; a cross-check only).

usage: python analysis/leibak_equality.py CENSUS.jsonl[.xz] [...]
"""
import collections
import itertools
import json
import lzma
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from engine.pari import pari        # noqa: E402


def quadratic(D):
    """Is there a Z-basis (u, v) of O_K, disc K = D, with N(u) + N(v) + N(u - v) = 0?  Such a
    basis has norm form [A, 2(A + C), C], so D = 4(A^2 + AC + C^2), and conversely such a
    form is a norm form of O_K iff it is GL_2(Z)-equivalent to +-(x^2 - (D/4) y^2), i.e.
    properly represents 1 or -1."""
    if D % 4:
        return False
    d = D // 4
    s = math.isqrt(4 * d // 3) + 1
    for A in range(-s, s + 1):
        for C in range(-s, s + 1):
            if A * A + A * C + C * C == d:
                f = pari('Qfb(%d,%d,%d)' % (A, 2 * (A + C), C))
                if any(len(pari.qfbsolve(f, n)) for n in (1, -1)):
                    return True
    return False


def cubic(pol):
    """Is there beta >> 0 with (beta) = p^2 d^-1 (2O = p q^2) and Tr(beta w^2) even for the
    integral basis w?  Units are taken modulo squares (8 classes), which preserve both
    conditions."""
    bnf = pari.bnfinit(pari(pol), 1)
    nf = pari('(b)->b.nf')(bnf)
    dec = pari.idealprimedec(nf, 2)
    if sorted((int(P[2]), int(P[3])) for P in dec) != [(1, 1), (2, 1)]:
        return False
    P = [P for P in dec if int(P[2]) == 1][0]
    I = pari.idealmul(nf, pari.idealpow(nf, P, 2), pari.idealinv(nf, pari('(nf)->nf.diff')(nf)))
    e, g = pari.bnfisprincipal(bnf, I)
    if any(int(c) for c in e):
        return False
    u1, u2 = list(pari('(b)->b.fu')(bnf))
    w = pari('(nf)->nf.zk')(nf)
    for s, a, b in itertools.product((1, -1), (0, 1), (0, 1)):
        beta = s * pari.nfeltmul(nf, g, pari.nfeltmul(nf, pari.nfeltpow(nf, u1, a), pari.nfeltpow(nf, u2, b)))
        if any(int(x) != 1 for x in pari.nfeltsign(nf, beta)):
            continue
        ba = pari.nfbasistoalg(nf, beta)
        if all(int(pari.nfelttrace(nf, pari.nfalgtobasis(nf, ba * wi * wi))) % 2 == 0 for wi in w):
            return True
    return False


def main(paths):
    for p in paths:
        op = lzma.open if p.endswith('.xz') else open
        c, bad, ram = collections.Counter(), [], collections.Counter()
        with op(p, 'rt') as fh:
            for line in fh:
                r = json.loads(line)
                if r['status'] != 'CERTIFIED':
                    continue
                att = r['leibak_ratio'] == '1'
                crit = quadratic(r['disc']) if r['degree'] == 2 else cubic(r['pol']) if r['degree'] == 3 else None
                c[(att, crit)] += 1
                if att != crit:
                    bad.append(r['label'])
                if att:
                    ram[r['disc'] % 16] += 1
        print('%s: (attained, criterion) -> fields %s; disagreements %s; D mod 16 of the attained: %s' % (
            os.path.basename(p), dict(c), bad[:20], dict(ram)))


if __name__ == '__main__':
    main(sys.argv[1:])
