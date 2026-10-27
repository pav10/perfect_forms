"""Field lists (plan, B2): totally real fields of bounded discriminant.

Two independent sources are compared before a list is used:

* PARI's nflist (cubic: groups C3 and S3, signature 0, i.e. totally real), and
* the LMFDB's totally real fields as mirrored by rvisser7/indecomposables
  (``totally_real_fields/degN.txt``: ``lmfdb_index|coeffs|monogenic`` with the
  coefficients ascending and the leading 1 omitted), which also gives the LMFDB
  labels n.n.D.i.

Real quadratic fields are generated directly (fundamental discriminants).
Each list line is a JSON object {"label", "pol", "D"} with pol in polredabs form.
"""
import json
import os
import urllib.request

from .pari import pari

MIRROR = 'https://raw.githubusercontent.com/rvisser7/indecomposables/main/totally_real_fields/deg%d.txt'


def fetch_mirror(degree, dest):
    urllib.request.urlretrieve(MIRROR % degree, dest)
    return dest


def read_mirror(path, disc_max=None, disc_min=1):
    """[(label, polredabs string, D)] from a mirror file, sorted by (D, index)."""
    out = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('lmfdb') or '|' not in line:
                continue
            idx, co, *_ = line.split('|')
            coeffs = [int(t) for t in co.split(',')] + [1]
            n = len(coeffs) - 1
            pol = pari.Pol(list(reversed(coeffs)))
            D = int(pari.nfdisc(pol))
            if D < disc_min:
                continue
            if disc_max is not None and D > disc_max:
                break
            out.append(('%d.%d.%d.%s' % (n, n, D, idx), str(pari.polredabs(pol)), D))
    return out


def quadratic(disc_max, disc_min=1):
    """Real quadratic fields Q(sqrt d) with fundamental discriminant D <= disc_max."""
    out = []
    for D in range(max(5, disc_min), disc_max + 1):
        if int(pari.isfundamental(D)):
            pol = str(pari.polredabs(pari('x^2 - %d' % D)))
            out.append(('2.2.%d.1' % D, pol, D))
    return out


def nflist_cubic(disc_max, disc_min=1):
    L = list(pari('nflist("C3",[%d,%d],0)' % (disc_min, disc_max))) + \
        list(pari('nflist("S3",[%d,%d],0)' % (disc_min, disc_max)))
    return sorted((int(pari.nfdisc(p)), str(pari.polredabs(p))) for p in L)


def build(degree, disc_max, disc_min=1, mirror=None):
    """The list, cross-checked: for degree 3 nflist and the mirror must agree field by
    field (same polredabs polynomials with the same discriminants)."""
    if degree == 2:
        return quadratic(disc_max, disc_min), {'source': 'isfundamental', 'checked': 'n/a'}
    if mirror is None:
        raise ValueError('need the LMFDB mirror file for degree %d' % degree)
    lst = read_mirror(mirror, disc_max, disc_min)
    info = {'source': 'LMFDB (rvisser7 mirror)', 'count': len(lst)}
    if degree == 3:
        nfl = nflist_cubic(disc_max, disc_min)
        mine = sorted((D, pol) for _, pol, D in lst)
        if nfl != mine:
            raise ArithmeticError('nflist and LMFDB lists differ (%d vs %d fields)' % (len(nfl), len(mine)))
        info['checked'] = 'nflist C3+S3: identical (%d fields)' % len(nfl)
    else:
        info['checked'] = 'none (LMFDB only)'
    return lst, info


def write(lst, path):
    with open(path, 'w') as fh:
        for label, pol, D in lst:
            fh.write(json.dumps({'label': label, 'pol': pol, 'D': D}) + '\n')


def read(path, disc_min=None, disc_max=None):
    import lzma
    out = []
    with (lzma.open(path, 'rt') if path.endswith('.xz') else open(path)) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            r = json.loads(line)
            if disc_min is not None and r['D'] < disc_min:
                continue
            if disc_max is not None and r['D'] > disc_max:
                continue
            out.append((r['label'], r['pol'], r['D']))
    return out
