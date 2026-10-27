"""Compare perfect-form invariants (n, mu, ell) with sail invariants of pav10/sails.

Sail data: the summary file of the sails releases (read-only), one line per field with,
per sign class sigma (representatives of GF(2)^3 / sgn(O^x)), the number of
indecomposables and of sail facets modulo totally positive units, and S.  Summed over
the classes these count orbits of the full unit group O^x on all 2^m orthants, the
normalisation of n_K, mu_K and ell_O.

usage: python analysis/sails_link.py CENSUS.jsonl[.xz] SAILS_SUMMARY.xz
"""
import ast
import json
import lzma
import math
import statistics
import sys


def read_sails(path):
    out = {}
    with lzma.open(path, 'rt') as fh:
        header = fh.readline().rstrip('\n').split('|')
        for line in fh:
            row = dict(zip(header, line.rstrip('\n').split('|')))
            ints = lambda k: [int(t) for t in row[k].strip('[]').split(',')]
            sig = ast.literal_eval(row['sail signatures'])
            out[row['label']] = {
                'classes': len(sig), 'tp_index': 0,
                'ind': ints('# indecomposables'), 'fac': ints('# facets'), 'S': ints('S'),
                'maxnorm': ints('max norm of indecomposables'),
            }
    return out


def read_census(path):
    op = lzma.open if path.endswith('.xz') else open
    with op(path, 'rt') as fh:
        return [json.loads(l) for l in fh if l.strip()]


def corr(x, y):
    mx, my = statistics.mean(x), statistics.mean(y)
    sx = math.sqrt(sum((a - mx) ** 2 for a in x))
    sy = math.sqrt(sum((b - my) ** 2 for b in y))
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / (sx * sy)


def main(census, sails):
    S = read_sails(sails)
    rows = []
    for r in read_census(census):
        s = S.get(r['label'])
        if s is None:
            continue
        rows.append({
            'label': r['label'], 'D': r['disc'], 'h': r['class_number'], 'R': r['regulator'],
            'n': r['n'], 'mu': r['mu'], 'ell': r.get('ell'),
            'classes': s['classes'],
            'F_tp': s['fac'][0], 'F_all': sum(s['fac']),
            'I_tp': s['ind'][0], 'I_all': sum(s['ind']),
            'S_tp': s['S'][0], 'S_all': sum(s['S']),
        })
    print('%d fields joined' % len(rows))
    keys = ['n', 'mu', 'ell', 'F_tp', 'F_all', 'I_tp', 'I_all', 'S_tp', 'S_all', 'R']
    print('correlations (Pearson) with n, mu, ell:')
    for k in keys:
        xs = [r[k] for r in rows if r['ell'] is not None]
        print('  %-6s  n: %.3f  mu: %.3f  ell: %.3f' % (
            k, corr([r['n'] for r in rows if r['ell'] is not None], xs),
            corr([r['mu'] for r in rows if r['ell'] is not None], xs),
            corr([r['ell'] for r in rows if r['ell'] is not None], xs)))
    for a, b in [('mu', 'F_all'), ('n', 'F_all'), ('ell', 'F_all'), ('mu', 'F_tp'), ('n', 'I_all'), ('ell', 'I_all'),
                 ('mu', 'S_all'), ('n', 'S_all')]:
        q = [r[a] / r[b] for r in rows if r[a] is not None and r[b]]
        print('%s/%s: min %.3f median %.3f max %.3f' % (a, b, min(q), statistics.median(q), max(q)))
    return rows


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
