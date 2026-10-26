# Findings

Labels as in CLAUDE.md: **PROVED**, **CERTIFIED** (certified computation; the range is
stated), **CONJECTURE**.  [M] = paper/perfect-unary-cubic.tex.

## 1. Perfect forms and sails of totally real cubic fields (plan: Yasaki-type link)

Notation (all counts modulo the full unit group O^x, over all 2^3 orthants; for the
sails this is the sum over the sign classes of pav10/sails of the counts modulo totally
positive units): n = n_K, mu = mu_K (extreme squares), ell = ell_O (relative minima),
V, E, F = vertices, edges, facets of the sails (Klein polyhedra) of all sign classes.

Data: `data/derived/sails_vs_perfect_deg3_D2e4.jsonl.xz` (all 832 fields with D <= 20000),
produced by `analysis/sails_elements.py` (this engine for n, mu, ell; the sails code of
pav10/sails, run read-only, for the sails; elements compared modulo units by the exact
key HNF(xO)).

* **PROVED.**  Every sail vertex is a relative minimum, so V <= ell.
  *Proof.*  Let v be a vertex of the Klein polyhedron of its orthant and y != +-v a
  lattice point with |sigma_i(y)| <= |sigma_i(v)| for all i (strictly, since equality at
  one place forces y = +-v).  Then v + y and v - y lie in the open orthant of v, and v is
  their midpoint, so v is not a vertex.  []
* **PROVED** ([M, Lemma 6.3]).  Every extreme square is a relative minimum: mu <= ell.
* **CERTIFIED** (D <= 20000).  mu <= V <= ell in every field; hence (Euler bound
  [M, Thm 4.3]) n <= 2 mu <= 2 V.  V = ell in 658 of 832 fields; V/ell in [0.5, 1].
* **CERTIFIED** (D <= 20000).  Of the 15,090 classes of extreme squares, 15,083 are sail
  vertices.  The other 7 (one in each of 3.3.4001.1, 13877.1, 14388.1, 14408.1, 15124.1,
  17144.1, 17581.1) are each the barycentre (v1 + v2 + v3)/3 of a triangular facet of the
  sail of their sign class at integer distance 1; none lies on a sail edge.  Example:
  3.3.4001.1 (x^3 - 11x - 7), x = -t^2 + t + 12 of norm 13: a totally positive unit
  multiple of x is (1 + v2 + v3)/3 for the facet {1, v2, v3} of the totally positive sail,
  yet q(x) is a vertex of the squared sail (the centroid of the three perfect forms at x
  has x as its only minimal vector).  So "extreme squares are sail vertices" is false.
* **CERTIFIED** (D <= 20000).  Every extreme square is a sail vertex or the barycentre of a
  triangular sail facet at integer distance 1 (all 832 fields).  Such barycentres are
  rarely extreme: of the 190 classes of triangular distance-1 facets containing their
  barycentre, 7 give extreme squares.
* **CONJECTURE.**  For every totally real cubic field, every extreme square is a vertex of
  the sail of its sign class or the barycentre of a triangular sail facet at integer
  distance 1.
* **PROVED** (relative minima on a sail).  A lattice point x in the open orthant C is a
  relative minimum iff no w != 0 has x + w and x - w both in C, iff
  K ∩ (2x - K) ∩ L = {x} for the Klein polyhedron K of C.  Consequently (i) a lattice point
  in the relative interior of a sail edge is never a relative minimum (its two lattice
  neighbours on the edge are x +- w); (ii) a lattice point x in the relative interior of a
  facet F is a relative minimum iff F ∩ (2x - F) contains no lattice point besides x
  (if x +- w are both in C, then c.w = 0 for the facet functional c, since no lattice point
  of C lies below the facet plane).
* **CERTIFIED** (D <= 20000).  The relative minima that are not sail vertices lie in the
  interior of triangular facets at distance 1; besides barycentres they occur e.g. at
  barycentric coordinates (1/5, 2/5, 2/5) (3.3.404.1, a 5-point triangle) and
  (1/7, 3/7, 3/7) (3.3.1257.1, a 6-point triangle).  "Relative minima = vertices and
  barycentres" holds in 710 of 832 fields only.
* **CERTIFIED** (D <= 20000).  Every relative minimum lies on the sail of its sign class
  (it is a lattice point of some facet).
* **CERTIFIED** (D <= 20000), ratios: mu/V in [0.22, 1] (median 0.61); n/V in [0.44, 2]
  (median 1.2); n/F in [0.33, 1.53] (median 0.91); n/E in [0.19, 0.87].
  n <= F fails in 109 fields, n <= E holds in all.
* The 27 fields with mu = 1 fall into three sail types: (V, F) = (1, 2) for the 10
  cyclic ones (D = 49, 81, 169, 361, 1369, 3969, 6241, 9409, 13689, 19321), (3, 4) for
  the 7 Ennola fields a = 4..10 (D = 697, ..., 16609), and (4, 4) for the other 10
  (D = 985, 2597, 5521, 6809, 7249, 10309, 13801, 14945, 16357, 17609).

## 2. Perfect unary forms, cubic census D <= 20000 (CERTIFIED)

Reproduces [M, Obs. 5.6] exactly (27 fields with mu = 1; max n = 152 at D = 16689; 68.0%
simplicial; mu/Reg median 0.616 in [0.238, 0.948] over h = 1, Reg > 10).  New:

* n_K takes the values 2, 4, 5, 6, 7, 8, ...; n_K = 1 and n_K = 3 do not occur.
* Leibak's bound gamma_K^3 <= 2 D is attained (some trace lattice is similar to A_3) in
  56 fields, all with D = 0 mod 4 (3.3.316.1, 940.1, 1708.1, ...).
* The three large fields of [M, Sec. 5] are reproduced: D = 676829, 58271449, 223164633
  give (mu, n) = (136, 272), (546, 1092), (1767, 3534), all simplicial.

## 3. Binary forms over real quadratic fields, additive normalisation (CERTIFIED)

**Agreement with the literature.**  Every count computed so far equals Gunnells--Yasaki,
*Perfect forms over totally real number fields* (arXiv:0908.3130), Table 1 (GL_2(O)-classes
by discriminant D): D = 5, 8, 12, 13, 17, 24, 28, 40, 44 give 2, 2, 3, 9, 34, 22, 30, 52, 66;
and Leibak's thesis Thm 3.1 (exactly 22 classes over Q(sqrt 6)).  Leibak's Thm 2.2 remark
(the bound gamma_{2,K} <= gamma_4 D^(1/2) is attained for d = 2, 3, 6) agrees with
gamma^4 = (2D)^2 below.

Perfect binary forms Tr(v^T A v) modulo GL_2(O_K) (engine/binary.py; certificates: edge
balance and the mass formula sum (-1)^dim/|Stab| = 0):

| d | 5 | 2 | 3 | 13 | 17 | 6 | 7 | 10 (h = 2) | 11 |
|---|---|---|---|---|---|---|---|---|---|
| perfect classes mod GL_2(O_K) | 2 | 2 | 3 | 9 | 34 | 22 | 30 | 52 | 66 |
| ... up to Galois as well | 2 | 2 | 3 | 6 | | 13 | 17 | 29 | 36 |

Classes with equal invariants are Galois conjugates (Gal(K/Q) is not in GL_2(O_K)).
For d = 2 (mod 4) and 3 (mod 4) in this range the maximum is gamma^4 = (2D)^2, i.e.
gamma_{2,K} = gamma_4 D^(1/2) (a trace lattice similar to D_4); not for d = 5, 13.

## 4. The floor: is n_K >= 2 for totally real cubic fields? (plan C1)

* **CERTIFIED.**  n_K >= 2 for every totally real cubic field with D <= 68297 (the census
  so far); the values n_K = 1 and n_K = 3 never occur.
* **PROVED** (reduction).  If n_K = 1, then by the Euler identity n = 2 mu - sum(deg - 3)
  and n >= mu/2 [M, Thm 4.3] either
  (a) mu = 1 and the unique vertex has degree 4: the torus is one quadrilateral (V, E, F) =
      (1, 2, 1); or
  (b) mu = 2 and the unique vertex has degree 6: two triangles (V, E, F) = (1, 3, 2).
  In case (b) the trace lattice (O, Tr(alpha x^2)) has 6 pairs of minimal vectors, so it is
  similar to A_3 (the unique 3-dimensional lattice with kissing number 12), and gamma_K
  attains Leibak's bound gamma_K^3 = 2 D.
  In case (a), after multiplying alpha by x^2 for an extreme square x, the extreme minimal
  vectors are +-1, +-e1, +-e2, +-e1 e2 for a basis e1, e2 of O^x/+-1 (the face of K_1 at
  the extreme square 1 has the four vertices alpha, alpha e1^2, alpha e2^2, alpha e1^2 e2^2,
  and the side pairings of the one-square torus are the two generators), and the Voronoi
  polygon is the quadrilateral 1, e1^2, e1^2 e2^2, e2^2 (in this cyclic order).  Hence its
  diagonals meet: mu + (1 - mu) a b = nu a + (1 - nu) b in K (x) R for a = e1^2, b = e2^2
  and some mu, nu in (0, 1) -- a rational linear relation among 1, a, b, ab with
  coefficient sum 0.  Moreover |Tr(alpha y)| <= 1/2 for y = e1, e2, e1 e2, e1^2 e2, e1 e2^2
  (from Tr(alpha (z +- w)^2) >= 1 for the pairs z, w of minimal vectors).
* Open: exclude (a) and (b).  Positivity of the embeddings alone does not exclude the
  relation in (a) (it defines a curve through (1, 1) in the (log a, log b) plane), so a
  proof must use the arithmetic of the units or the minimality of alpha.

## 5. Real quadratic fields (plan B1)

Census of all real quadratic fields with D <= 10^5 (running; numbers below for the first
24,023 fields, all certified), with ell = period length of the continued fraction of
omega_d (= the relative minima count; checked against relmin.py for D < 200).

* **CERTIFIED.**  n_K = 1 exactly for the d of Leibak--Porter--Ling's types T1-T4
  (JTNB 35 (2023), Thm 1.2) -- 322 fields; this independently confirms their theorem in
  this range and validates the engine.
* **CERTIFIED.**  Yasaki's bound l_d/2 <= n_K <= l_d (JTNB 25 (2013), as quoted by
  Porter--Mendelsohn) holds in every field; both ends are attained (n = l in 2142 fields,
  n = l/2 in 842).  n/l has median 0.79.
* ell/Reg has median 0.842 for Reg > 20: Levy's constant 12 log 2/pi^2 = 0.8428 of the
  continued fraction algorithm.

## 6. The cubic analogue of Yasaki's bound (plan Goal 1)

In degree 2, ell is the continued fraction period and l/2 <= n <= l.  In degree 3, ell is
Voronoi's period (relative minima mod units).  Cubic census D <= 10^5 (4804 fields):

* **PROVED.**  n <= 2 mu <= 2 ell (Euler bound [M, Thm 4.3] and [M, Lemma 6.3]).
* **CERTIFIED** (D <= 10^5).  (4/9) ell <= n <= 2 ell; n/ell has median 1.16 in every
  discriminant window (0.444 .. 2.0); n < ell/2 only in 4 fields (n = 4, ell = 9:
  D = 11197, 45684, 49897, 61825).  mu/ell in [2/9, 1], median 0.59.
* ell/Reg drifts slowly: medians 1.02, 0.97, 0.95, 0.94, 0.93 over D windows of width 2*10^4
  (Reg > 10) -- a "cubic Levy constant" for Voronoi's algorithm, if it converges.
* **CONJECTURE.**  c ell <= n_K for an absolute c > 0 (the data suggest c = 4/9).
