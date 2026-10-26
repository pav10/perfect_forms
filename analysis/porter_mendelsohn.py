"""Compare n_K with the bound of Porter--Mendelsohn (IJNT 2024; arXiv:2208.03304), Thm 1:

    n_K <= e^{2 n rho} gamma_K^n / |Delta_K|^{1/2},     rho = rho_inf(Lambda_K),

where Lambda_K = Log(O_K^x) = {(log|sigma_i(u)|)_i} in H = {sum x_i = 0} and
rho_inf(Lambda) = max_{x in H} min_{l in Lambda} ||x - l||_inf (their Section 1), and with their
conjectured improvement e^{2 n rho} gamma_K^n / |Delta_K| (plan B4).  gamma_K^n is exact (records);
rho is numerical: Reg/2 in degree 2, and for degree 3 the maximum of the distance function over
a fine grid of the fundamental domain refined by local optimisation (label: numerical, not
certified).

usage: python analysis/porter_mendelsohn.py CENSUS.jsonl[.xz] [OUT.jsonl]
"""
import json
import lzma
import math
import os
import sys
from fractions import Fraction

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from engine.field import Order        # noqa: E402


def covering_radius_inf(B):
    """rho_inf of the lattice spanned by the rows of B (r x n, rows in H), r = 1 or 2."""
    B = np.asarray(B, dtype=float)
    r = B.shape[0]
    if r == 1:
        return 0.5 * np.abs(B[0]).max()
    # Lagrange (Gauss) reduction of the 2-dimensional basis (Euclidean)
    b1, b2 = B[0], B[1]
    while True:
        if b1 @ b1 > b2 @ b2:
            b1, b2 = b2, b1
        mu = round((b1 @ b2) / (b1 @ b1))
        if mu == 0:
            break
        b2 = b2 - mu * b1
    pts = np.array([i * b1 + j * b2 for i in range(-3, 4) for j in range(-3, 4)])

    def f(st):
        x = st[0] * b1 + st[1] * b2
        return np.abs(pts - x).max(axis=1).min()

    N = 120
    g = np.linspace(0, 1, N, endpoint=False)
    S, T = np.meshgrid(g, g, indexing='ij')
    X = S.reshape(-1, 1) * b1 + T.reshape(-1, 1) * b2                   # (N^2, n)
    F = np.abs(X[:, None, :] - pts[None, :, :]).max(axis=2).min(axis=1)
    k = int(F.argmax())
    best, arg = float(F[k]), (S.reshape(-1)[k], T.reshape(-1)[k])
    # local refinement (pattern search)
    step = 1.0 / N
    s, t = arg
    while step > 1e-12:
        moved = False
        for ds, dt in ((step, 0), (-step, 0), (0, step), (0, -step), (step, step), (-step, -step),
                       (step, -step), (-step, step)):
            v = f((s + ds, t + dt))
            if v > best:
                best, s, t, moved = v, s + ds, t + dt, True
        if not moved:
            step /= 2
    return best


def main(path, out=None):
    op = lzma.open if path.endswith('.xz') else open
    with op(path, 'rt') as fh:
        recs = [json.loads(l) for l in fh if l.strip()]
    rows = []
    for r in recs:
        O = Order(r['pol'])
        n = O.m
        B = [O.log_embed(u) for u in O.reduced_units()]
        rho = covering_radius_inf(B)
        g = Fraction(r['gamma_pow'])                    # gamma_K^n, exact
        D = r['disc']
        log_rhs = 2 * n * rho + math.log(g) - 0.5 * math.log(D)
        log_conj = 2 * n * rho + math.log(g) - math.log(D)
        rows.append({'label': r['label'], 'D': D, 'n': r['n'], 'Reg': r['regulator'], 'rho_inf': float(rho),
                     'log_bound': log_rhs, 'log_conjectured_bound': log_conj,
                     'slack_log': log_rhs - math.log(r['n']),
                     'conjecture_ok': bool(math.log(r['n']) <= log_conj)})
    if out:
        with open(out, 'w') as fh:
            for row in rows:
                fh.write(json.dumps(row) + '\n')
    ok = sum(math.log(x['n']) <= x['log_bound'] for x in rows)
    cj = sum(x['conjecture_ok'] for x in rows)
    print('%d fields: PM bound holds in %d; conjectured bound (/|Delta|) holds in %d' % (len(rows), ok, cj))
    sl = sorted(x['slack_log'] for x in rows)
    print('log(bound/n_K): min %.2f, median %.2f, max %.2f' % (sl[0], sl[len(sl) // 2], sl[-1]))
    bad = sorted((x for x in rows if not x['conjecture_ok']), key=lambda x: x['D'])
    print('conjecture fails at:', [(x['label'], x['n'], round(math.exp(x['log_conjectured_bound']), 3)) for x in bad[:20]])
    return rows


if __name__ == '__main__':
    main(*sys.argv[1:])
