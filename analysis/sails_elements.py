"""Element-level comparison of sail vertices, relative minima and extreme squares.

For each field, computes with this engine the extreme squares (Voronoi) and the
relative minima (relmin.py), and with pav10/sails (imported read-only from SAILS_PATH,
default ../sails; no bytecode is written there) the sails of all sign classes.  All
elements are compared modulo units through the exact key HNF(xO).

Reports per field: V_all (sail vertices mod O^x over all sign classes), F_all, ell, mu, n,
and the inclusions between the three sets.

usage: python analysis/sails_elements.py FIELDLIST.jsonl[.xz] DISC_MAX [OUT.jsonl] [--jobs J] [--disc-min D]
       python analysis/sails_elements.py --census CENSUS.jsonl[.xz] OUT.jsonl [JOBS]   (extreme squares from the records)
       python analysis/sails_elements.py --facets CENSUS.jsonl[.xz] CENSUS_OUT.jsonl OUT.jsonl
           (facets and barycentric coordinates of the extreme squares that are not sail vertices)
"""
import json
import os
import sys
import time

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.append(os.environ.get('SAILS_PATH', os.path.join(os.path.dirname(os.path.dirname(HERE)), 'sails')))

from engine import fieldlists                     # noqa: E402
from engine.field import Order                    # noqa: E402
from engine.pari import pari                      # noqa: E402
from engine.relmin import RelativeMinima          # noqa: E402
from engine.unary import Voronoi                  # noqa: E402
from sails import compute_field                   # noqa: E402


def sail_vertex_keys(O, pol):
    """Sail data of all sign classes, keyed by HNF(xO): vertices, lattice points on edges
    (ridges) and all lattice points of facets; the number of facets and of edges (ridges,
    each shared by two facets) modulo units; and the lattice barycentres of the triangular
    facets at integer distance 1 only (`barycentres` in the output; German points in
    general, any k and distance, come from --facets)."""
    coeffs = [int(c) for c in reversed(list(pari.Vec(pari(pol))))]
    r = compute_field(coeffs).record()
    basis = [pari(b.replace('a', 'x')) for b in r['integral_basis']]
    cache = {}

    def key(v):
        v = tuple(v)
        if v not in cache:
            a = sum(int(c) * b for c, b in zip(v, basis))
            x = O.as_int(O.from_alg(pari.Mod(a, O.pol)))
            cache[v] = (O.ideal_key(x), x)
        return cache[v]

    verts, edge_pts, lattice, bary = {}, {}, {}, {}
    facets = ridges = 0
    for orth in r['orthants']:
        facets += orth['num_facets']
        for F in orth['facets']:
            ridges += len(F['ridges'])
            for v in F['vertices']:
                k, x = key(v)
                verts.setdefault(k, x)
            for rd in F['ridges']:
                for v in rd['points']:
                    k, x = key(v)
                    edge_pts.setdefault(k, x)
            for v in F['lattice_points']:
                k, x = key(v)
                lattice.setdefault(k, x)
            if len(F['vertices']) == 3 and F['d'] == 1:
                vs = [tuple(v) for v in F['vertices']]
                for v in F['lattice_points']:
                    if tuple(v) not in vs and all(3 * v[i] == sum(w[i] for w in vs) for i in range(len(v))):
                        k, x = key(v)
                        bary.setdefault(k, x)
    return verts, facets, lattice, edge_pts, ridges // 2, bary


def one(args):
    label, pol, D = args
    t = time.time()
    O = Order(pol)
    V = Voronoi(O).run()
    f = V.f_vector()
    R = RelativeMinima(O).run()
    rk = set(R.by_key)
    ek = set(V.extreme_square_classes())
    vk, F_all, lat, edg, E_all, bary = sail_vertex_keys(O, pol)
    vks = set(vk)
    pos = {'vertex': 0, 'edge': 0, 'facet': 0, 'off': 0}
    for k in ek:
        pos['vertex' if k in vks else 'edge' if k in edg else 'facet' if k in lat else 'off'] += 1
    return {'label': label, 'D': D, 'h': O.class_number(), 'R': O.regulator(),
            'n': f[-1], 'e': f[1], 'mu': f[0], 'ell': len(rk), 'V_all': len(vks),
            'F_all': F_all, 'E_all': E_all, 'L_all': len(lat), 'ext_position': pos,
            'barycentres': len(bary), 'barycentres_extreme': len(set(bary) & ek),
            'ext_vertex_or_barycentre': ek <= (vks | set(bary)),
            'relmin_eq_vertices_or_barycentres': rk == (vks | set(bary)),
            'relmin_in_vertices': rk <= vks, 'vertices_in_relmin': vks <= rk,
            'ext_in_vertices': ek <= vks, 'relmin_on_sail': rk <= set(lat),
            'relmin_not_vertex': len(rk - vks), 'vertex_not_relmin': len(vks - rk),
            'seconds': round(time.time() - t, 2)}


_CENSUS = {}


def one_from_census(rec):
    """Like one(), but the extreme squares are read from a census record (forms with the
    flag 'extreme') and ell from the record; relative minima are not recomputed."""
    t = time.time()
    O = Order(rec['pol'])
    ek = {}
    for F in rec['forms']:
        for x, e in zip(F['minvecs'], F['extreme']):
            if e:
                ek.setdefault(O.ideal_key(tuple(x)), tuple(x))
    if len(ek) != rec['mu']:
        raise ArithmeticError('census record: extreme squares do not match mu')
    vk, F_all, lat, edg, E_all, bary = sail_vertex_keys(O, rec['pol'])
    vks = set(vk)
    pos = {'vertex': 0, 'edge': 0, 'facet': 0, 'off': 0}
    for k in ek:
        pos['vertex' if k in vks else 'edge' if k in edg else 'facet' if k in lat else 'off'] += 1
    ekk = set(ek)
    return {'label': rec['label'], 'D': rec['disc'], 'h': rec['class_number'], 'R': rec['regulator'],
            'n': rec['n'], 'e': rec['e'], 'mu': rec['mu'], 'ell': rec.get('ell'), 'V_all': len(vks),
            'F_all': F_all, 'E_all': E_all, 'L_all': len(lat), 'ext_position': pos,
            'barycentres': len(bary), 'barycentres_extreme': len(set(bary) & ekk),
            'ext_vertex_or_barycentre': ekk <= (vks | set(bary)),
            'ext_in_vertices': ekk <= vks, 'seconds': round(time.time() - t, 2)}


def facet_data(rec):
    """For each extreme square of a census record that is not a sail vertex: the facets of
    the sails containing it (number of vertices, integer distance d, number of lattice
    points) and its exact barycentric coordinates when the facet is a triangle."""
    from fractions import Fraction
    from engine import linalg
    O = Order(rec['pol'])
    ek = {}
    for F in rec['forms']:
        for x, e in zip(F['minvecs'], F['extreme']):
            if e:
                ek.setdefault(O.ideal_key(tuple(x)), tuple(x))
    coeffs = [int(c) for c in reversed(list(pari.Vec(pari(rec['pol']))))]
    r = compute_field(coeffs).record()
    basis = [pari(b.replace('a', 'x')) for b in r['integral_basis']]

    def key(v):
        a = sum(int(c) * b for c, b in zip(v, basis))
        return O.ideal_key(O.as_int(O.from_alg(pari.Mod(a, O.pol))))
    verts = {key(v) for orth in r['orthants'] for F in orth['facets'] for v in F['vertices']}
    out = []
    for orth in r['orthants']:
        for F in orth['facets']:
            for v in F['lattice_points']:
                k = key(v)
                if k in ek and k not in verts:
                    row = {'x': list(ek[k]), 'norm': str(O.norm(ek[k])), 'facet_vertices': len(F['vertices']),
                           'd': F['d'], 'lattice_points': len(F['lattice_points'])}
                    if len(F['vertices']) == 3:
                        M = [[Fraction(F['vertices'][j][i]) for j in range(3)] for i in range(3)]
                        row['barycentric'] = [str(t) for t in linalg.solve(M, [Fraction(t) for t in v])]
                    out.append(row)
    return {'label': rec['label'], 'D': rec['disc'], 'non_vertex_extreme': out}


def main_facets(census, se_out, out):
    """facet_data() for the fields of a --census output with extreme squares off the vertices."""
    import lzma
    with open(se_out) as fh:
        flagged = {json.loads(l)['label'] for l in fh if json.loads(l)['ext_position']['facet'] or
                   json.loads(l)['ext_position']['edge'] or json.loads(l)['ext_position']['off']}
    op = lzma.open if census.endswith('.xz') else open
    with op(census, 'rt') as fh, open(out, 'w') as fo:
        for line in fh:
            rec = json.loads(line)
            if rec['label'] in flagged:
                fo.write(json.dumps(facet_data(rec)) + '\n')


def main_census(census, out, jobs=1):
    import lzma
    import multiprocessing as mp
    op = lzma.open if census.endswith('.xz') else open
    with op(census, 'rt') as fh:
        recs = [json.loads(l) for l in fh if l.strip()]
    done = set()
    if os.path.exists(out):
        with open(out) as fh:
            done = {json.loads(l)['label'] for l in fh if l.strip()}
    todo = [r for r in recs if r['label'] not in done and r['status'] == 'CERTIFIED']
    with open(out, 'a') as fh, mp.get_context('fork').Pool(int(jobs), maxtasksperchild=100) as pool:
        for rec in pool.imap_unordered(one_from_census, todo, chunksize=4):
            fh.write(json.dumps(rec) + '\n')
            fh.flush()


def main(fields, disc_max, out=None, jobs=1, disc_min=None):
    import multiprocessing as mp
    rows = fieldlists.read(fields, disc_min=int(disc_min) if disc_min else None, disc_max=int(disc_max))
    done = set()
    if out and os.path.exists(out):
        with open(out) as fh:
            done = {json.loads(l)['label'] for l in fh if l.strip()}
    todo = [r for r in rows if r[0] not in done]
    fh = open(out, 'a') if out else None
    with mp.get_context('fork').Pool(int(jobs), maxtasksperchild=100) as pool:
        for rec in pool.imap_unordered(one, todo, chunksize=1):
            line = json.dumps(rec)
            if fh:
                fh.write(line + '\n')
                fh.flush()
            else:
                print(line, flush=True)


if __name__ == '__main__':
    import argparse
    if len(sys.argv) > 1 and sys.argv[1] == '--facets':
        # python analysis/sails_elements.py --facets CENSUS.jsonl[.xz] CENSUS_OUT.jsonl OUT.jsonl
        main_facets(sys.argv[2], sys.argv[3], sys.argv[4])
        sys.exit(0)
    if len(sys.argv) > 1 and sys.argv[1] == '--census':
        # python analysis/sails_elements.py --census CENSUS.jsonl[.xz] OUT.jsonl [JOBS]
        main_census(sys.argv[2], sys.argv[3], sys.argv[4] if len(sys.argv) > 4 else 1)
        sys.exit(0)
    p = argparse.ArgumentParser()
    p.add_argument('fields')
    p.add_argument('disc_max')
    p.add_argument('out', nargs='?')
    p.add_argument('--jobs', default=1)
    p.add_argument('--disc-min')
    a = p.parse_args()
    main(a.fields, a.disc_max, a.out, a.jobs, a.disc_min)
