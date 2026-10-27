"""Exact face lattice of the Voronoi cone of a perfect form.

At a perfect form alpha (minimum 1) the minimal vectors x give integer vectors
r(x) with Tr(alpha x^2) = r(x) . alpha = 1.  The Voronoi cone is cone{r(x)}; it
is full dimensional (perfection) and its cross-section r(x) . alpha = 1 is an
(m-1)-polytope P(alpha), the dual polytope conv{q(x)} of the plan.  Faces are
represented by the set of indices of the points lying on them (non-vertex points
included), facets also by a primitive integer normal beta with
r(x) . beta >= 0 for all x and = 0 exactly on the facet.  Everything is exact
integer linear algebra; this replaces the polygon code (vu_edges) of
voronoi_unary.gp, which handled m = 3 only.
"""
from itertools import combinations

from . import linalg


def facets(points):
    """Facets of cone(points) (points: integer vectors spanning R^m and lying on an
    affine hyperplane missing 0).  Returns a list of (beta, frozenset of indices),
    sorted by the index sets."""
    k = len(points)
    m = len(points[0])
    if linalg.rank(points) != m:
        raise ValueError('points do not span: form is not perfect')
    found = {}
    for S in combinations(range(k), m - 1):
        if any(set(S) <= F for F in found):
            continue
        ker = linalg.kernel([points[i] for i in S])
        if len(ker) != 1:
            continue
        b = ker[0]
        vals = [linalg.dot(b, p) for p in points]
        if all(v >= 0 for v in vals):
            pass
        elif all(v <= 0 for v in vals):
            b = tuple(-t for t in b)
            vals = [-v for v in vals]
        else:
            continue
        F = frozenset(i for i, v in enumerate(vals) if v == 0)
        if F not in found:
            found[F] = b
    return sorted(((b, F) for F, b in found.items()), key=lambda t: sorted(t[1]))


def face_lattice(points, facet_list=None):
    """All proper nonempty faces of P = conv(points) on the hyperplane, by dimension:
    {dim: sorted list of frozensets of point indices}, dim = 0 .. m-2.  The faces are
    the nonempty intersections of facets (closure under intersection)."""
    m = len(points[0])
    if facet_list is None:
        facet_list = facets(points)
    faces = {F for _, F in facet_list}
    frontier = set(faces)
    while frontier:
        new = set()
        for A in frontier:
            for B in faces:
                C = A & B
                if C and C not in faces and C not in new:
                    new.add(C)
        faces |= new
        frontier = new
    out = {d: [] for d in range(m - 1)}
    for F in faces:
        d = linalg.rank([points[i] for i in F]) - 1
        out[d].append(F)
    for d in out:
        out[d].sort(key=sorted)
    return out


def check_polytope_euler(faces, m):
    """Euler--Poincare for the boundary of an (m-1)-polytope:
    sum_{j=0}^{m-2} (-1)^j f_j = 1 - (-1)^(m-1)."""
    s = sum((-1) ** d * len(faces[d]) for d in faces)
    return s == 1 - (-1) ** (m - 1)
