"""Compare a binary census (analysis/binary_census.py) with the literature and with
pav10/hermite-humbert (read-only).

* Gunnells--Yasaki, Table 1 (arXiv:0908.3130): number of GL_2(O_K)-classes of perfect binary
  forms by discriminant D (the table is the GY_TABLE of tests/test_engine.py).
* Leibak's bound gamma_{2,K}^4 <= gamma_4^4 D^2 = 4 D^2 (thesis Thm 2.2) and its equality
  cases.
* The multiplicative Hermite--Humbert constant gamma*_{K,2} of pav10/hermite-humbert
  (Icaza's normalisation: S[v] = S1[v1] S2[v2], gamma*(S) = m(S)/sqrt(det S1 det S2),
  certified there with closed forms).  By AM-GM, S1[v1] + S2[v2] >= 2 sqrt(S1[v1] S2[v2]),
  so the additive constant satisfies gamma_{2,K} >= 2 sqrt(gamma*_{K,2}) (Leibak's thesis,
  Prop. 2.1, with n = 2), i.e. gamma_{2,K}^4 >= 16 gamma*^2.  Checked numerically
  (gamma* is given there to 25 digits; gamma_{2,K}^4 is exact here).

usage: python analysis/binary_compare.py CENSUS.jsonl [HERMITE_HUMBERT_TSV]
(default TSV: $HH_PATH or ../hermite-humbert/results/gamma_K2.tsv)
"""
import json
import os
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'tests'))

from test_engine import GY_TABLE       # noqa: E402


def read_hh(path):
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        head = fh.readline().rstrip('\n').split('\t')
        for line in fh:
            row = dict(zip(head, line.rstrip('\n').split('\t')))
            if row.get('status', '').startswith('certified'):
                out[int(row['D'])] = (float(row['gamma']), row['closed_form'])
    return out


def main(census, hh_path=None):
    hh_path = hh_path or os.environ.get('HH_PATH', os.path.join(os.path.dirname(ROOT), 'hermite-humbert',
                                                                 'results', 'gamma_K2.tsv'))
    with open(census) as fh:
        R = sorted((json.loads(l) for l in fh if l.strip()), key=lambda r: r['d'])
    bad = [r['d'] for r in R if r['status'] != 'CERTIFIED']
    R = [r for r in R if r['status'] == 'CERTIFIED']
    print('%d certified fields (d = %d..%d); not certified: %s' % (len(R), R[0]['d'], R[-1]['d'], bad))

    gy = [(r['D'], r['n_perfect'], GY_TABLE[r['D']]) for r in R if r['D'] in GY_TABLE]
    print('Gunnells--Yasaki Table 1: %d discriminants compared, %d agree; disagreements: %s' % (
        len(gy), sum(a == b for _, a, b in gy), [t for t in gy if t[1] != t[2]]))
    missing = sorted(set(GY_TABLE) - {r['D'] for r in R})
    if missing:
        print('  GY entries not in this census: %s' % missing)

    eq = [r for r in R if Fraction(r['gamma_pow_max']) == 4 * r['D'] ** 2]
    over = [r['d'] for r in R if Fraction(r['gamma_pow_max']) > 4 * r['D'] ** 2]
    ram2 = [r for r in R if r['D'] % 4 == 0]
    print("Leibak's bound gamma^4 <= 4 D^2: violated for %s; equality for %d fields, all with 2 | D: %s; "
          "fields with 2 | D: %d, of which equality in %d" % (
              over or 'none', len(eq), all(r['D'] % 4 == 0 for r in eq), len(ram2),
              sum(Fraction(r['gamma_pow_max']) == 4 * r['D'] ** 2 for r in ram2)))
    print('  equality: d = %s' % [r['d'] for r in eq])
    print('  2 | D without equality: d = %s' % [r['d'] for r in ram2 if r not in eq])

    hh = read_hh(hh_path)
    rows = []
    for r in R:
        if r['d'] in hh:
            gstar, cf = hh[r['d']]
            g4 = Fraction(r['gamma_pow_max'])
            rows.append((r['d'], float(g4), 16 * gstar ** 2, float(g4) / (16 * gstar ** 2), cf))
    if rows:
        ok = sum(a >= b * (1 - 1e-12) for _, a, b, _, _ in rows)
        print('hermite-humbert (%s): %d fields with both constants; gamma_add^4 >= 16 gamma_mult^2 in %d' % (
            os.path.basename(hh_path), len(rows), ok))
        print('  ratio gamma_add^4 / (16 gamma_mult^2): min %.4f (d = %d), max %.4f (d = %d)' % (
            min(t[3] for t in rows), min(rows, key=lambda t: t[3])[0],
            max(t[3] for t in rows), max(rows, key=lambda t: t[3])[0]))
        tight = [t for t in rows if abs(t[3] - 1) < 1e-9]
        print('  equality (ratio 1): d = %s' % [t[0] for t in tight])
    else:
        print('hermite-humbert table not found at %s' % hh_path)
    return R


if __name__ == '__main__':
    main(*sys.argv[1:])
