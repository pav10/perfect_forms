"""Statistics of a cubic census (plan B2-B5, C1-C4), printed as Markdown.

usage: python analysis/cubic_census.py CENSUS.jsonl[.xz] [...]
"""
import json
import lzma
import math
import statistics
import sys
from collections import Counter
from fractions import Fraction


def read(paths):
    out = {}
    for p in paths:
        op = lzma.open if p.endswith('.xz') else open
        with op(p, 'rt') as fh:
            for line in fh:
                if line.strip():
                    r = json.loads(line)
                    if r['status'] == 'CERTIFIED':
                        out[r['label']] = r
    return sorted(out.values(), key=lambda r: r['disc'])


def q(xs, f='%.3f'):
    xs = sorted(xs)
    if not xs:
        return '-'
    return (f + ' [' + f + ', ' + f + ']') % (statistics.median(xs), xs[0], xs[-1])


def main(paths):
    R = read(paths)
    print('## Cubic census: %d certified fields, D = %d..%d\n' % (len(R), R[0]['disc'], R[-1]['disc']))
    ns = Counter(r['n'] for r in R)
    print('* n_K values: smallest %s; n_K = 1: %d fields, n_K = 3: %d fields' % (
        sorted(ns)[:12], ns.get(1, 0), ns.get(3, 0)))
    print('* mu_K = 1: %d fields; max n_K = %d at %s' % (
        sum(r['mu'] == 1 for r in R), max(ns), max(R, key=lambda r: r['n'])['label']))
    print('* every certificate passed (records are CERTIFIED only); relative minima checked: %d' % (
        sum('relmin' in r['certificates'] for r in R)))
    sand = [r for r in R if 'ell' in r]
    lo = sum(4 * r['regulator'] / math.log(2 * r['disc']) ** 2 <= r['mu'] for r in sand)
    hi = sum(r['ell'] <= 4 * r['regulator'] / math.log(2) ** 2 for r in sand)
    print('* regulator sandwich [M, Thm 6.5]: lower bound holds in %d/%d, mu <= ell in %d/%d, upper bound in %d/%d\n' % (
        lo, len(sand), sum(r['mu'] <= r['ell'] for r in sand), len(sand), hi, len(sand)))

    print('### mu/Reg and ell/Reg by discriminant window (Reg > 10)\n')
    print('| D window | fields | mu/Reg, h=1 | mu/Reg, h>1 | ell/Reg | n/ell | simplicial |')
    print('|---|---|---|---|---|---|---|')
    edges = [0, 2 * 10 ** 4, 4 * 10 ** 4, 6 * 10 ** 4, 8 * 10 ** 4, 10 ** 5, 10 ** 9]
    for a, b in zip(edges, edges[1:]):
        W = [r for r in R if a < r['disc'] <= b]
        if not W:
            continue
        Wr = [r for r in W if r['regulator'] > 10]
        h1 = [r['mu'] / r['regulator'] for r in Wr if r['class_number'] == 1]
        hg = [r['mu'] / r['regulator'] for r in Wr if r['class_number'] > 1]
        el = [r['ell'] / r['regulator'] for r in Wr if 'ell' in r]
        ne = [r['n'] / r['ell'] for r in Wr if 'ell' in r]
        simp = sum(set(map(int, r['degrees'])) == {3} for r in W) / len(W)
        print('| (%d, %d] | %d | %s | %s | %s | %s | %.1f%% |' % (a, min(b, R[-1]['disc']), len(W), q(h1), q(hg), q(el), q(ne), 100 * simp))
    cyc = [r for r in R if math.isqrt(r['disc']) ** 2 == r['disc']]
    ncy = [r for r in R if math.isqrt(r['disc']) ** 2 != r['disc']]
    print('\n* cyclic (D a square): %d fields, mu/Reg %s; non-Galois: mu/Reg %s (Reg > 10)' % (
        len(cyc), q([r['mu'] / r['regulator'] for r in cyc if r['regulator'] > 10]),
        q([r['mu'] / r['regulator'] for r in ncy if r['regulator'] > 10])))
    degs = Counter(int(d) for r in R for d, c in r['degrees'].items() for _ in range(c))
    print('* vertex degrees over all classes: %s; maximal degree 6 in %d fields' % (
        dict(sorted(degs.items())), sum(max(map(int, r['degrees'])) == 6 for r in R)))
    att = [r for r in R if r.get('leibak_ratio') == '1']
    print('* gamma_K = gamma_3 D^(1/3) (Leibak bound attained): %d fields; D mod 4 of these: %s' % (
        len(att), dict(Counter(r['disc'] % 4 for r in att))))
    lr = [float(Fraction(r['leibak_ratio'])) for r in R]
    print('* gamma_K^3 / (2 D): %s' % q(lr))
    ext = [r['n_extreme'] / r['n'] for r in R]
    print('* fraction of perfect classes that are extreme (eutactic): %s; weakly eutactic classes in %d fields' % (
        q(ext), sum(r['n_weakly_eutactic'] > 0 for r in R)))


if __name__ == '__main__':
    main(sys.argv[1:])
