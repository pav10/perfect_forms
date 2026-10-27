"""Command line interface:  python -m engine COMMAND ...

  field POL       one field (polynomial in x), prints the record
  fields          build a cross-checked field list (JSONL)
  run             a resumable, parallel batch over a field list
  summary         tables from result files
"""
import argparse
import json
import multiprocessing as mp
import os
import signal
import sys
import time
import traceback


def _compute(args):
    label, pol, opts = args
    from .record import compute

    def handler(signum, frame):
        raise TimeoutError('timeout after %ss' % opts['timeout'])

    if opts['timeout']:
        signal.signal(signal.SIGALRM, handler)
        signal.alarm(opts['timeout'])
    try:
        rec = compute(pol, label=label, forms=opts['forms'], relmin=opts['relmin'])
        return label, pol, rec, None
    except Exception as e:      # noqa: BLE001 - reported in the errors file
        return label, pol, None, '%s: %s' % (type(e).__name__, e) + (
            '\n' + traceback.format_exc() if opts['traceback'] else '')
    finally:
        if opts['timeout']:
            signal.alarm(0)


def _done_labels(*paths):
    done = set()
    for p in paths:
        if os.path.exists(p):
            with open(p) as fh:
                for line in fh:
                    line = line.strip()
                    if line:
                        try:
                            done.add(json.loads(line)['label'])
                        except (ValueError, KeyError):
                            pass
    return done


def cmd_field(a):
    from .record import compute
    rec = compute(a.pol, label=a.label, forms=not a.no_forms, relmin=a.relmin)
    if a.json:
        print(json.dumps(rec))
    else:
        keys = ['label', 'pol', 'disc', 'class_number', 'regulator', 'n', 'e', 'mu', 'f_vector',
                'degrees', 'facet_sizes', 'extreme_square_norms', 'gamma_pow', 'gamma', 'n_extreme',
                'leibak_ratio', 'ell', 'status']
        for k in keys:
            if k in rec:
                print('%-22s %s' % (k, rec[k]))
        print('%-22s %s' % ('certificates', rec['certificates']))
        print('%-22s %s' % ('seconds', rec['provenance']['seconds']))


def cmd_fields(a):
    from . import fieldlists
    lst, info = fieldlists.build(a.degree, a.disc_max, a.disc_min, mirror=a.mirror)
    fieldlists.write(lst, a.out)
    print('%d fields written to %s; %s' % (len(lst), a.out, info), file=sys.stderr)


def cmd_run(a):
    from . import fieldlists
    fields = fieldlists.read(a.fields, a.disc_min, a.disc_max)
    if a.limit:
        fields = fields[:a.limit]
    qpath, epath = a.out + '.quarantine', a.out + '.errors'
    done = _done_labels(a.out, qpath) | (set() if a.retry_errors else _done_labels(epath))
    todo = [(lab, pol, {'timeout': a.timeout, 'forms': not a.no_forms, 'relmin': a.relmin,
                        'traceback': a.traceback}) for lab, pol, _ in fields if lab not in done]
    print('%d fields, %d done, %d to do' % (len(fields), len(fields) - len(todo), len(todo)), file=sys.stderr)
    t0 = time.time()
    k = nerr = nq = 0
    ctx = mp.get_context('fork')
    with open(a.out, 'a') as out, open(qpath, 'a') as qf, open(epath, 'a') as ef, \
            ctx.Pool(a.jobs, maxtasksperchild=200) as pool:
        for label, pol, rec, err in pool.imap_unordered(_compute, todo, chunksize=1):
            k += 1
            if err is not None:
                nerr += 1
                ef.write(json.dumps({'label': label, 'pol': pol, 'error': err}) + '\n')
                ef.flush()
            elif rec['status'] == 'CERTIFIED':
                out.write(json.dumps(rec) + '\n')
                out.flush()
            else:
                nq += 1
                qf.write(json.dumps(rec) + '\n')
                qf.flush()
            if a.progress and k % a.progress == 0:
                print('%d/%d done (%d errors, %d quarantined) in %.0fs' % (k, len(todo), nerr, nq, time.time() - t0),
                      file=sys.stderr)
    print('finished: %d computed, %d errors, %d quarantined, %.0fs' % (k, nerr, nq, time.time() - t0), file=sys.stderr)


def cmd_summary(a):
    from . import summary
    summary.main(a.files, a.what)


def main(argv=None):
    p = argparse.ArgumentParser(prog='python -m engine', description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='cmd', required=True)

    q = sub.add_parser('field', help='one field')
    q.add_argument('pol', help='defining polynomial in x, e.g. "x^3-x^2-2*x+1"')
    q.add_argument('--label')
    q.add_argument('--json', action='store_true', help='print the full JSON record')
    q.add_argument('--no-forms', action='store_true', help='omit the list of perfect forms')
    q.add_argument('--relmin', action='store_true', help='also run the relative-minima cross-check')
    q.set_defaults(func=cmd_field)

    q = sub.add_parser('fields', help='build a cross-checked field list')
    q.add_argument('--degree', type=int, required=True)
    q.add_argument('--disc-max', type=int, required=True)
    q.add_argument('--disc-min', type=int, default=1)
    q.add_argument('--mirror', help='LMFDB mirror file degN.txt (degrees >= 3)')
    q.add_argument('--out', required=True)
    q.set_defaults(func=cmd_fields)

    q = sub.add_parser('run', help='batch over a field list (resumable)')
    q.add_argument('--fields', required=True)
    q.add_argument('--out', required=True, help='JSONL of certified records; also OUT.quarantine, OUT.errors')
    q.add_argument('--disc-min', type=int)
    q.add_argument('--disc-max', type=int)
    q.add_argument('--limit', type=int)
    q.add_argument('--jobs', type=int, default=os.cpu_count())
    q.add_argument('--timeout', type=int, default=0, help='seconds per field (0 = none)')
    q.add_argument('--no-forms', action='store_true')
    q.add_argument('--relmin', action='store_true')
    q.add_argument('--retry-errors', action='store_true')
    q.add_argument('--traceback', action='store_true')
    q.add_argument('--progress', type=int, default=0)
    q.set_defaults(func=cmd_run)

    q = sub.add_parser('summary', help='summary tables of result files')
    q.add_argument('files', nargs='+')
    q.add_argument('--what', default='overview', choices=['overview', 'lines'])
    q.set_defaults(func=cmd_summary)

    a = p.parse_args(argv)
    a.func(a)


if __name__ == '__main__':
    main()
