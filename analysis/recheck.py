"""Recompute census records with the current code and compare them field by field.

Every field of the record except `provenance` and `stats` (timings and counters) must
agree.  Used to (a) replace records whose provenance says the code was uncommitted
(`+dirty`) and (b) spot-check records computed by earlier commits.

usage: python analysis/recheck.py CENSUS.jsonl OUT_REPORT.json [--sample K] [--seed S]
           [--replace-dirty FIXED.jsonl] [--jobs J]
"""
import json
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from engine.record import compute       # noqa: E402

SKIP = {'provenance', 'stats'}


def _recompute(r):
    new = compute(r['pol'], label=r['label'], forms='forms' in r, relmin='ell' in r and r['degree'] != 2)
    new = json.loads(json.dumps(new))          # same normalisation as the stored record
    diffs = [k for k in set(r) | set(new) if k not in SKIP and r.get(k) != new.get(k)]
    return r['label'], new, diffs


def main():
    import argparse
    import multiprocessing as mp
    p = argparse.ArgumentParser()
    p.add_argument('census')
    p.add_argument('report')
    p.add_argument('--sample', type=int, default=0)
    p.add_argument('--seed', type=int, default=1)
    p.add_argument('--replace-dirty')
    p.add_argument('--jobs', type=int, default=1)
    a = p.parse_args()
    with open(a.census) as fh:
        recs = [json.loads(l) for l in fh if l.strip()]
    dirty = [r for r in recs if r['provenance']['git'].endswith('+dirty')]
    clean = [r for r in recs if not r['provenance']['git'].endswith('+dirty')]
    random.seed(a.seed)
    sample = random.sample(clean, min(a.sample, len(clean))) if a.sample else []
    todo = dirty + sample if a.replace_dirty else sample
    report = {'records': len(recs), 'dirty': len(dirty), 'sample': len(sample), 'mismatches': []}
    fixed = {}
    with mp.get_context('fork').Pool(a.jobs) as pool:
        for label, new, diffs in pool.imap_unordered(_recompute, todo):
            if diffs:
                report['mismatches'].append({'label': label, 'keys': sorted(diffs)})
            fixed[label] = new
    if a.replace_dirty:
        dl = {r['label'] for r in dirty}
        with open(a.replace_dirty, 'w') as fh:
            for r in recs:
                fh.write(json.dumps(fixed[r['label']] if r['label'] in dl else r) + '\n')
    with open(a.report, 'w') as fh:
        json.dump(report, fh, indent=1)
    print(json.dumps({k: (v if k != 'mismatches' else len(v)) for k, v in report.items()}))


if __name__ == '__main__':
    main()
