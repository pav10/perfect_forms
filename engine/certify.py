"""Certificates of a Voronoi computation (plan, section 1.1).

A record is CERTIFIED only if all of the following hold; otherwise it is
quarantined (written to a separate file) and never repaired by hand.

* euler: the alternating sum of the face numbers of the torus complex
  dK_1(m_O)/O^x is 0 (it is a (d-1)-torus, manuscript Prop. 7.1);
* unit_index: the side-pairing units generate a subgroup of index 1 in the unit
  group of bnfinit.  Since they generate O^x/+-1 unconditionally (manuscript
  Prop. 5.3), this also proves that bnfinit's units are fundamental (no GRH);
* mu_keys: the number of classes of extreme squares counted by the gluing of the
  face complex equals the number of distinct ideals xO over the extreme squares x
  (an independent exact count of mu);
* edges: sum of vertex degrees = 2 * (number of edges);
* local Euler--Poincare of every Voronoi cone and all gluing maps (checked while
  running; a failure raises VoronoiError).
"""


def certificates(O, V):
    f = V.f_vector()
    m = O.m
    euler = V.euler()
    idx = O.subgroup_index(V.side_pairing_units())
    mu_keys = len(V.extreme_square_classes())
    edges = f[m - 2] if m >= 2 else 0
    out = {
        'euler': euler,
        'euler_ok': euler == 0,
        'unit_index': int(idx) if idx == int(idx) else str(idx),
        'unit_index_ok': idx == 1,
        'mu_keys': mu_keys,
        'mu_keys_ok': mu_keys == f[0],
        'edges_ok': sum(V.degrees()) == 2 * edges,
    }
    out['certified'] = all(out[k] for k in ('euler_ok', 'unit_index_ok', 'mu_keys_ok', 'edges_ok'))
    return out


def relmin_certificates(O, V, R):
    """Cross-check with the independent relative-minima enumerator (relmin.py):
    every extreme square is a relative minimum (manuscript Lemma 6.3), so
    mu <= ell, and the relative-minima side pairings also generate the units."""
    keys = set(R.by_key)
    ext = V.extreme_square_classes()
    idx = O.subgroup_index(R.side_pairing_units())
    out = {
        'ell': R.count(),
        'extreme_squares_are_relmin': all(k in keys for k in ext),
        'relmin_unit_index': int(idx) if idx == int(idx) else str(idx),
    }
    out['ok'] = out['extreme_squares_are_relmin'] and idx == 1 and len(ext) <= R.count()
    return out
