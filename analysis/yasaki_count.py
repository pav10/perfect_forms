"""Yasaki's global count (JTNB 25 (2013), Sec. 1 and 4): for the 121,580 real quadratic fields
Q(sqrt d), 2 <= d <= 200000 squarefree, there are 10,732,735 perfect unary forms in total
(h_d = GL_1(O)-orbits of homothety classes = n_K).

Reads certified census records (any number of files; the union must contain every such
field exactly once), and checks along the way Yasaki's Thm 5.1 (l_d/2 <= h_d <= l_d) and
Leibak--Porter--Ling's Thm 1.2 (h_d = 1 iff d is of type T1-T4) on the whole range.

usage: python analysis/yasaki_count.py CENSUS.jsonl[.xz] [...]
"""
import json
import lzma
import math
import os
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from engine.pari import pari        # noqa: E402

D_MAX = 200000
YASAKI_FIELDS, YASAKI_FORMS = 121580, 10732735


def lpl(d):
    """Leibak--Porter--Ling, Thm 1.2: h_d = 1 iff d = m^2+1 (m odd), m^2-1 (m even),
    m^2+4 (m odd) or m^2-4 (m > 3 odd)."""
    m = math.isqrt(d)
    for k in (m - 1, m, m + 1, m + 2):
        if k < 1:
            continue
        if (d == k * k + 1 and k % 2) or (d == k * k - 1 and k % 2 == 0) or \
                (d == k * k + 4 and k % 2) or (d == k * k - 4 and k % 2 and k > 3):
            return True
    return False


def main(paths):
    want = {}
    for d in range(2, D_MAX + 1):
        if int(pari.issquarefree(d)):
            want[d if d % 4 == 1 else 4 * d] = d
    got = {}
    for p in paths:
        op = lzma.open if p.endswith('.xz') else open
        with op(p, 'rt') as fh:
            for line in fh:
                r = json.loads(line)
                if r['degree'] != 2 or r['disc'] not in want:
                    continue
                if r['status'] != 'CERTIFIED':
                    raise ValueError('%s is not CERTIFIED' % r['label'])
                if r['disc'] in got and got[r['disc']]['n'] != r['n']:
                    raise ValueError('%s: two records disagree' % r['label'])
                got[r['disc']] = r
    missing = sorted(set(want) - set(got))
    total = sum(r['n'] for r in got.values())
    bounds = sum(r['ell'] <= 2 * r['n'] and r['n'] <= r['ell'] for r in got.values())
    lpl_ok = sum((r['n'] == 1) == lpl(want[D]) for D, r in got.items())
    by_res = Counter()
    for D, r in got.items():
        by_res[want[D] % 4] += r['n']
    print('fields: %d of %d (Yasaki: %d); missing: %s' % (len(got), len(want), YASAKI_FIELDS, missing[:10]))
    print('total number of perfect forms: %d (Yasaki: %d) -> %s' % (
        total, YASAKI_FORMS, 'AGREES' if total == YASAKI_FORMS and not missing else 'DIFFERS'))
    print('Yasaki Thm 5.1 holds in %d fields; LPL Thm 1.2 in %d; h_d = 1 in %d fields' % (
        bounds, lpl_ok, sum(r['n'] == 1 for r in got.values())))
    print('forms by d mod 4: %s; max h_d = %d at d = %d' % (
        dict(sorted(by_res.items())), max(r['n'] for r in got.values()),
        want[max(got, key=lambda D: got[D]['n'])]))
    return total, missing


if __name__ == '__main__':
    main(sys.argv[1:])
