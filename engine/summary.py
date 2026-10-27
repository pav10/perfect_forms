"""Summary tables of result files (JSONL records, plain or .xz)."""
import json
import lzma
import statistics
from collections import Counter


def read_records(paths):
    for p in paths:
        op = lzma.open if p.endswith('.xz') else open
        with op(p, 'rt') as fh:
            for line in fh:
                line = line.strip()
                if line:
                    yield json.loads(line)


def line(r):
    return '%s|%d|%s|%d|%.6f|%d|%d|%d|%s|%s' % (
        r['label'], r['disc'], r['pol'], r['class_number'], r['regulator'], r['n'], r['e'], r['mu'],
        r['f_vector'], r['gamma_pow'])


def main(paths, what='overview'):
    recs = sorted(read_records(paths), key=lambda r: (r['disc'], r['label'] or ''))
    if what == 'lines':
        print('label|disc|pol|h|regulator|n|e|mu|f_vector|gamma^m')
        for r in recs:
            print(line(r))
        return
    if not recs:
        print('no records')
        return
    m = recs[0]['degree']
    cert = sum(r['status'] == 'CERTIFIED' for r in recs)
    print('records: %d (degree %d, disc %d..%d), certified: %d' % (len(recs), m, recs[0]['disc'], recs[-1]['disc'], cert))
    ns = [r['n'] for r in recs]
    mus = [r['mu'] for r in recs]
    print('n_K:  min %d, median %s, max %d (at %s)' % (min(ns), statistics.median(ns), max(ns),
                                                    max(recs, key=lambda r: r['n'])['label']))
    print('mu_K: min %d, max %d; mu = 1 for %d fields' % (min(mus), max(mus), sum(x == 1 for x in mus)))
    print('n_K distribution (smallest values):', sorted(Counter(ns).items())[:12])
    simp = sum(set(map(int, r['degrees'])) == {m} for r in recs)
    print('simplicial (every vertex of degree %d): %d (%.1f%%)' % (m, simp, 100.0 * simp / len(recs)))
    rat = [r['mu'] / r['regulator'] for r in recs if r['class_number'] == 1 and r['regulator'] > 10]
    if rat:
        print('mu/Reg (h = 1, Reg > 10): %d fields, median %.3f, range [%.3f, %.3f]' % (
            len(rat), statistics.median(rat), min(rat), max(rat)))
    att = [r['label'] for r in recs if r.get('leibak_ratio') == '1']
    print('Leibak bound gamma_K = gamma_m D^(1/m) attained: %d fields %s' % (len(att), att[:20]))
