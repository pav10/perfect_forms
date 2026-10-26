"""Element-level comparison of sail vertices, relative minima and extreme squares.

For each field, computes with this engine the extreme squares (Voronoi) and the
relative minima (relmin.py), and with pav10/sails (imported read-only from SAILS_PATH,
default ../sails; no bytecode is written there) the sails of all sign classes.  All
elements are compared modulo units through the exact key HNF(xO).

Reports per field: V_all (sail vertices mod O^x over all sign classes), F_all, ell, mu, n,
and the inclusions between the three sets.

usage: python analysis/sails_elements.py FIELDLIST.jsonl[.xz] DISC_MAX [OUT.jsonl]
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
    (ridges) and all lattice points of facets; and the number of facets and of edges
    (ridges, each shared by two facets) modulo units."""
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

    verts, edge_pts, lattice = {}, {}, {}
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
    return verts, facets, lattice, edge_pts, ridges // 2


def main(fields, disc_max, out=None):
    rows = fieldlists.read(fields, disc_max=int(disc_max))
    fh = open(out, 'w') if out else None
    tot = {'relmin_in_vertices': 0, 'vertices_in_relmin': 0, 'ext_in_vertices': 0, 'relmin_on_sail': 0, 'n': 0}
    for label, pol, D in rows:
        t = time.time()
        O = Order(pol)
        V = Voronoi(O).run()
        f = V.f_vector()
        R = RelativeMinima(O).run()
        rk = set(R.by_key)
        ek = set(V.extreme_square_classes())
        vk, F_all, lat, edg, E_all = sail_vertex_keys(O, pol)
        vks = set(vk)
        pos = {'vertex': 0, 'edge': 0, 'facet': 0, 'off': 0}
        for k in ek:
            pos['vertex' if k in vks else 'edge' if k in edg else 'facet' if k in lat else 'off'] += 1
        rec = {'label': label, 'D': D, 'n': f[-1], 'e': f[1], 'mu': f[0], 'ell': len(rk), 'V_all': len(vks),
               'F_all': F_all, 'E_all': E_all, 'L_all': len(lat), 'ext_position': pos,
               'relmin_in_vertices': rk <= vks, 'vertices_in_relmin': vks <= rk,
               'ext_in_vertices': ek <= vks, 'relmin_on_sail': rk <= set(lat),
               'relmin_not_vertex': len(rk - vks), 'vertex_not_relmin': len(vks - rk),
               'seconds': round(time.time() - t, 2)}
        for k in ('relmin_in_vertices', 'vertices_in_relmin', 'ext_in_vertices', 'relmin_on_sail'):
            tot[k] += rec[k]
        tot['n'] += 1
        print(json.dumps(rec), flush=True)
        if fh:
            fh.write(json.dumps(rec) + '\n')
    print('TOTAL', json.dumps(tot))


if __name__ == '__main__':
    main(*sys.argv[1:])
