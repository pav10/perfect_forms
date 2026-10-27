# Findings

Labels as in CLAUDE.md: **PROVED**, **CERTIFIED** (certified computation; the range is
stated), **CONJECTURE**.  [M] = paper/perfect-unary-cubic.tex.

## 1. Perfect forms and sails of totally real cubic fields (plan: Yasaki-type link)

Notation (all counts modulo the full unit group O^x, over all 2^3 orthants; for the
sails this is the sum over the sign classes of pav10/sails of the counts modulo totally
positive units): n = n_K, mu = mu_K (extreme squares), ell = ell_O (relative minima),
V, E, F = vertices, edges, facets of the sails (Klein polyhedra) of all sign classes.

Data: `data/derived/sails_vs_perfect_deg3_D2e4.jsonl.xz` (all 832 fields with D <= 20000,
with relative minima as elements) and `sails_vs_perfect_deg3_D1e5.jsonl.xz` (all 4804
fields with D <= 10^5; extreme squares from the census records, `--census`), produced by
`analysis/sails_elements.py` (this engine for n, mu, ell; the sails code of pav10/sails,
run read-only, for the sails; elements compared modulo units by the exact key HNF(xO));
`sails_nonvertex_extreme_deg3_D1e5.jsonl.xz` lists the facets of the extreme squares that
are not sail vertices (`--facets`).

**Summary, D <= 10^5 (4804 fields, CERTIFIED).**  mu <= V <= ell in every field (V = ell in
3054); n <= 2V and n <= E in every field, n <= F fails in 333.  Of the 202,071 classes of
extreme squares, 201,948 are sail vertices and 123 lie in the relative interior of a
triangular facet at integer distance 1 (none on an edge): 121 at the barycentre of a
triangle with 4 lattice points, and 2 (3.3.61108.1, 3.3.75425.1) at barycentric
coordinates (1/5, 2/5, 2/5) of a triangle with 5 lattice points.  Of the 2512 classes of
4-point distance-1 triangles, 121 have an extreme barycentre.  Ratios (median [min, max]):
mu/V 0.59 [0.22, 1], n/V 1.17 [0.44, 2], n/F 0.89 [0.33, 1.53], n/E 0.51 [0.19, 0.87],
stable across discriminant windows of width 2*10^4 (mu/V 0.61, 0.59, 0.59, 0.59, 0.59).
The statements below marked D <= 20000 use the relative minima as elements, which the
D <= 10^5 run does not recompute.

*Literature (checked against the texts).*  [By] V. A. Bykovskii, Local minima of lattices
and vertices of Klein polyhedra, Funct. Anal. Appl. 40 (2006) 56-57; [Ge] O. N. German,
Klein polyhedra and relative minima of lattices, Math. Notes 79 (2006) 505-510.  Both work
with a lattice in R^n and the Klein polyhedron K of an orthant; "local" = "relative"
minimum (for psi(O), K totally real, the definitions agree: |sigma_i(y)| = |sigma_i(x)|
forces y = +-x).  [Ge]'s exceptional lattices (the orbit L^3_{1/2}) have points on the
coordinate axes, so psi(O) is never exceptional, and psi(O) has no nonzero point on a
coordinate plane (N(y) != 0), so [Ge, Cor. 3] applies to every sign class.

* **PROVED** ([By, Theorem]; also [Ge, Cor. 1]).  Every sail vertex is a relative minimum,
  so V <= ell.  (Our proof: if y != +-v has |sigma_i(y)| <= |sigma_i(v)| for all i, then
  v +- y lie in the open orthant of v and v is their midpoint.)  In dimension 2 the two
  sets coincide ([By], "well known"; [Ge, Prop. 1]): for real quadratic fields the sail
  vertices of the two sign classes mod units are the ell convergents.
* **PROVED** ([Ge, Lemma 1]).  A lattice point x in the open orthant C is a relative minimum
  iff no w != 0 has x + w and x - w both in C (equivalently 2x is not a sum of two lattice
  points of C other than x, x).  Consequently a lattice point in the relative interior of
  a sail edge is never a relative minimum (its two lattice neighbours on the edge are
  x +- w).
* **PROVED** ([Ge, Thm 1]).  Every relative minimum lies on the sail of its sign class.
  (CERTIFIED independently for D <= 20000: each is a lattice point of some facet.)
* **PROVED** ([Ge, Thm 2 and Cor. 3]).  The relative minima in an orthant are exactly the
  sail vertices and German's points Xi: v = (k(u1 + u2) + u3)/(2k + 1), k >= 1, in a
  triangular face F = conv(u1, u2, u3) whose lattice points are u1, u2, u3 and the k points
  u3 + (m/k)(v - u3), m = 1..k (on the median from u3; v is the last).  Barycentric
  coordinates (k, k, 1)/(2k + 1): the barycentre for k = 1 (4 lattice points on F),
  (2/5, 2/5, 1/5) for k = 2 (5 points), (3/7, 3/7, 1/7) for k = 3 (6 points).
* **CERTIFIED** (D <= 20000; 174 fields have relative minima that are not sail vertices).
  These are German points with k = 1 (148 classes), 2 (49), 3 (14), 4 (4), 5 (2), 6 (2)
  and 7 (1; 3.3.19089.1, a 10-point triangle); e.g. k = 2 in 3.3.404.1, k = 3 in
  3.3.1257.1.  Their faces have integer distance 1, except for 17 barycentres (k = 1) on
  faces at distance 3 (3.3.3356.1, 6556.1, 7528.1, 9281.1, 10172.1, 10641.1, 10712.1,
  10949.1, 12092.1, 12284.1, 14165.1, 14316.1, 14360.1, 14769.1, 16689.1, 17581.1,
  19605.1).  Relative minima = sail vertices and barycentres (k = 1) in 763 of 832
  fields.
* **PROVED** ([Ge, Cor. 3] and [M, Lemma 6.3]).  Every extreme square is a sail vertex or a
  German point.
* **CERTIFIED** (D <= 20000).  Of the 15,090 classes of extreme squares, 15,083 are sail
  vertices.  The other 7 (one in each of 3.3.4001.1, 13877.1, 14388.1, 14408.1, 15124.1,
  17144.1, 17581.1) are German points with k = 1, the barycentre (v1 + v2 + v3)/3 of a
  triangular facet of the sail of their sign class at integer distance 1.  Example:
  3.3.4001.1 (x^3 - 11x - 7), x = -t^2 + t + 12 of norm 13: a totally positive unit
  multiple of x is (1 + v2 + v3)/3 for the facet {1, v2, v3} of the totally positive sail,
  yet q(x) is a vertex of the squared sail (the centroid of the three perfect forms at x
  has x as its only minimal vector).  So "extreme squares are sail vertices" is false.
  Such barycentres are rarely extreme: of the 190 classes of 4-point triangular
  distance-1 facets, 7 give extreme squares.
* **CERTIFIED** (D <= 10^5, summary above).  Extreme German points with k = 2 occur
  (3.3.61108.1, 3.3.75425.1), so the earlier conjecture "every extreme square is a sail
  vertex or a barycentre" is false (data/ERRATA.md); no extreme square with k >= 3 occurs.
* Open: can an extreme square be a German point with k >= 3, or lie on a face at integer
  distance > 1?  Relative minima of both kinds occur (above); for D <= 10^5 the 123
  non-vertex extreme squares all have k <= 2 and distance 1.
* **CERTIFIED** (D <= 20000), ratios: mu/V in [0.22, 1] (median 0.61); n/V in [0.44, 2]
  (median 1.2); n/F in [0.33, 1.53] (median 0.91); n/E in [0.19, 0.87].
  n <= F fails in 109 fields, n <= E holds in all.
* The 27 fields with mu = 1 fall into three sail types: (V, F) = (1, 2) for the 10
  cyclic ones (D = 49, 81, 169, 361, 1369, 3969, 6241, 9409, 13689, 19321), (3, 4) for
  the 7 Ennola fields a = 4..10 (D = 697, ..., 16609), and (4, 4) for the other 10
  (D = 985, 2597, 5521, 6809, 7249, 10309, 13801, 14945, 16357, 17609).

* **PROVED** (Bykovskii, German, Euler; manuscript Thm 7.4).  ell = V + G, where G is the
  number of German points mod units; each sail facet carries at most one German point
  (k is fixed by the number of lattice points of the facet, and for k >= 2 the interior
  points lie on one median), so G <= F <= 2V (Euler on the sail tori: V - E + F = 0, every
  facet has >= 3 sides).  Hence ell/3 <= V <= ell and, with the regulator sandwich,
  4 Reg/(3 (log 2D)^2) <= V <= 4 Reg/(log 2)^2, 2 Reg/(3 (log 2D)^2) <= F <= 8 Reg/(log 2)^2,
  and n <= 2 mu <= 2 ell <= 6V.  (For sails the volume argument alone gives facets
  O(sqrt(D) Reg); this is linear in Reg with an absolute constant, plus a lower bound.)
  CERTIFIED (D <= 10^5): E = V + F, F <= 2V, V >= ell/2 in every field.
* **PROVED** (manuscript Lemma 7.6).  A totally positive relative minimum is indecomposable
  (if x = y + z with y, z >> 0, y is dominated by x).  So totally positive extreme squares,
  sail vertices and German points are indecomposable.
* **CERTIFIED** (D <= 20000, all 832 fields).  Every relative minimum is an extreme square or
  adjacent to one (the extreme squares dominate the adjacency graph of relative minima);
  the stronger "every minimal triple contains an extreme square" fails for 135 of the
  50,176 minimal triples.  This is the local statement behind a possible lower bound
  n >= c ell (manuscript, open problem 4); the graph degrees are not bounded (26 occurs).

### Three hulls of one point set

Let P = {(|sigma_i(y)|)_i : y in O \ 0} in R^3_{>0} (all orthants folded together);
O^x acts on it by the diagonal maps diag(|sigma_i(eps)|).

1. *Klein sails*: for each orthant C, the boundary of conv(O cap C); vertices V, facets F.
2. *The squared sail* Q = conv(P^2) + R^3_{>=0} (P^2 = coordinatewise squares).  The Ryshkov
   polyhedron {a : sum a_i p_i >= 1 for all p in P^2} is its polar, so the extreme squares
   are the vertices of Q and the perfect forms its facets (a perfect form alpha is the facet
   {p : sum sigma_i(alpha) p_i = 1}): dQ/O^x is a polyhedral torus with (V, E, F) =
   (mu, e, n).
3. *The Minkowski staircase* dS, S = P + R^3_{>=0} (Karpenkov--Ustinov, JNT 176 (2017),
   Def. 2.8): its outer corners are the relative minima, its inner corners the minimal
   triples {x, y, z} (the box max(|x_i|, |y_i|, |z_i|) contains no other lattice point
   and each of x, y, z attains one of the maxima).  These are the vertices of the
   Minkowski--Voronoi complex MV(O) of KU, Def. 2.5 (2-faces = relative minima, edges =
   adjacent pairs).  A minimal triple is a box -- a tropical linear form max_i X_i/c_i --
   supporting P at three points, as a perfect form is a linear form supporting P^2.

* **PROVED.**  Every vertex of the folded hull conv(P) + R^3_{>=0} is an extreme square.
  *Proof.*  If x^2 is not extreme, q(x) = (sigma_i(x)^2) lies in conv(P^2 \ {q(x)}) +
  R^3_{>=0}: sigma_i(x)^2 >= sum_k lambda_k p_{k,i}^2 >= (sum_k lambda_k p_{k,i})^2 for each
  i (Jensen), so |sigma(x)| lies in conv(P \ {|sigma(x)|}) + R^3_{>=0} and is not a
  vertex.  []  So: vertices of the folded hull are in {extreme squares} and in {sail
  vertices}, and both are in {relative minima}.
* **PROVED.**  MV(O)/O^x has (V, E, F) = (2 ell, 3 ell, ell).
  *Proof.*  For totally real K, |sigma_i(y)| = |sigma_i(z)| forces y = +-z (sigma_i(y/z) =
  +-1), so P is in general position in KU's sense and every minimal triple has exactly
  three neighbouring triples (KU, Sec. 2.3.1): 2E = 3V.  Radial projection maps dS
  homeomorphically onto the trace-zero plane H (S is closed under scaling up, and every
  ray t*e^h enters S since the box of e^h contains lattice points for t large and none for
  N(t e^h) < 1), and units act on H by the translations Log(eps); so dS/O^x is the torus
  H/Log(O^x) and V - E + F = 0 with F = ell.  []
* **CERTIFIED** (D <= 20000, `analysis/minkowski_voronoi.py`,
  `data/derived/minkowski_voronoi_deg3_D2e4.jsonl.xz`).  Counting the minimal triples
  directly gives (V, E, F) = (2 ell, 3 ell, ell) in all 832 fields.
* **PROVED.**  n <= 2 mu <= 2 ell = V(MV): the number of perfect forms is at most the number
  of vertices of the Minkowski--Voronoi complex.  In degree 2, MV(O)/O^x is a circle with
  ell vertices and Yasaki's Thm 5.1 reads V(MV)/2 <= n <= V(MV); in degree 3 the data give
  (2/9) V(MV) <= n <= V(MV) (section 6).

## 2. Perfect unary forms, cubic census D <= 20000 (CERTIFIED)

Reproduces [M, Obs. 5.6] exactly (27 fields with mu = 1; max n = 152 at D = 16689; 68.0%
simplicial; mu/Reg median 0.616 in [0.238, 0.948] over h = 1, Reg > 10).  New:

* n_K takes the values 2, 4, 5, 6, 7, 8, ...; n_K = 1 and n_K = 3 do not occur.
* Leibak's bound gamma_K^3 <= 2 D is attained (some trace lattice is similar to A_3) in
  56 fields, all with D = 0 mod 4 (3.3.316.1, 940.1, 1708.1, ...); explained in section 7
  (2O = p q^2 is necessary).
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

Census of all 30,394 real quadratic fields with D <= 10^5 (`data/census/deg2_D1e5.jsonl.xz`,
all CERTIFIED; 1,673,343 classes of perfect forms in total, max n = 522 at 2.2.95209.1).
ell = period length of the continued fraction of omega = (D mod 2 + sqrt D)/2; omega differs
from Yasaki's c_d by an integer, so ell = Yasaki's l_d (and = the relative minima count,
checked against relmin.py for D < 200).

* **CERTIFIED.**  n_K = 1 exactly for the d of Leibak--Porter--Ling's types T1-T4
  (JTNB 35 (2023), Thm 1.2) -- 364 fields; this independently confirms their theorem in
  this range and validates the engine.
* **CERTIFIED.**  Yasaki's bound l_d/2 <= h_d <= l_d (JTNB 25 (2013), Thm 5.1; checked against
  the text) holds in every field; both ends are attained (n = l in 2497 fields, n = l/2 in
  974).  n/l has median 0.79, the same (0.79) for d = 1, 2, 3 mod 4.
* **CERTIFIED.**  Yasaki's examples (Sec. 5.4): Q(sqrt 19) has 4 classes, Q(sqrt 23) has 2, with
  exactly the minimal vectors listed there (tests/test_engine.py); every class for d < 400
  has 2 or 3 minimal vectors (his Cor. 3.1).
* ell/Reg has median 0.841 over the 23,673 fields with Reg > 20: Levy's constant
  12 log 2/pi^2 = 0.8428 of the continued fraction algorithm.  n/Reg has median 0.664.
* n = mu in every field (Euler characteristic 0 of the circle).  The fraction of eutactic
  classes per field has median 0.97; no weakly eutactic class occurs.
* Leibak's bound gamma_K^2 = (4/3) D is attained in 983 fields, all with d = 3 mod 4; the
  exact criterion is section 7.
* **CERTIFIED** (Yasaki's global count; `analysis/yasaki_count.py` on `deg2_D1e5.jsonl.xz`
  and `deg2_d2e5_rest.jsonl.xz`).  The 121,580 fields Q(sqrt d), 2 <= d <= 200000
  squarefree, carry exactly 10,732,735 classes of perfect unary forms in total, the number
  stated by Yasaki (JTNB 25 (2013), Secs. 1 and 4; the only number his paper tabulates,
  the rest being plots), from independent code: every one of the 121,580 records is
  CERTIFIED.  By d mod 4: 3,648,331 (d = 1), 3,247,134 (d = 2), 3,837,270 (d = 3);
  max h_d = 888 at d = 196771.  On this whole range his Thm 5.1 (l_d/2 <= h_d <= l_d) and
  Leibak--Porter--Ling's Thm 1.2 hold in every field (h_d = 1 in 685 fields, exactly the
  types T1-T4).

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

## 7. Equality in Leibak's bound (plan E4)

* **PROVED** (docs/THEORY.md section 4).  For 2 <= m <= 8, gamma_K^m = gamma_m^m D iff some
  beta >> 0 in the codifferent d^-1 has N(beta d) = det Lambda_m (3, 4, 4, 4, 3, 2, 1) and
  Tr(beta x^2) even for all x in O; then 2 ramifies in K.  For m = 2: iff O has a Z-basis
  (u, v) with N(u) + N(v) + N(u - v) = 0 (so d = 3 mod 4); for m = 3: then 2O = p q^2 and
  (beta) = p^2 d^-1.
* **CERTIFIED** (D <= 10^5, `analysis/leibak_equality.py`).  The criterion and the exact
  gamma_K of the census agree in all 30,394 quadratic fields (983 attain the bound; the
  quadratic test is exact and unconditional) and all 4804 cubic fields (306 attain it, every
  one with 2O = p q^2; 959 further fields with 2O = p q^2 do not; the cubic principal ideal
  test uses bnfinit, so this half of the cross-check is GRH-conditional).
* **CERTIFIED** (D <= 10^5).  Every field attaining the bound has D = 12 mod 16, in degree 2
  (forced by the theorem) and in degree 3 (where 2O = p q^2 alone allows v_2(D) = 3 as well).
* Consequence for section 4: in case (b) of the n_K = 1 reduction the trace lattice is
  similar to A_3, so a cubic field with n_K = 1 of type (b) has 2O = p q^2 and p^2 d^-1
  narrowly principal.  Together with the observed D mod 16 this puts type (b) in a thin
  2-adic family; case (a) is not affected.

## 8. The Porter--Mendelsohn bound (plan B4)

Porter--Mendelsohn (IJNT 2024; arXiv:2208.03304), Thm 1: n_K <= e^{2 n rho} gamma_K^n /
|Delta_K|^{1/2}, rho = rho_inf(Log O^x); their conjectured improvement replaces
|Delta_K|^{1/2} by |Delta_K|.  `analysis/porter_mendelsohn.py`; gamma_K^n exact, rho
numerical (Reg/2 in degree 2), so these are numerical checks, not certificates.

* Real quadratic, D <= 10^5 (30,394 fields, `data/derived/porter_mendelsohn_deg2_D1e5.jsonl.xz`):
  the bound and the conjectured bound hold in every field.  log(bound/n_K) has median 102
  (min 1.77 at D = 5); the conjectured bound is tightest for D = 5, 8, 13, 12, 17
  (log slack 0.96, 1.76, 2.39, 2.92, 3.21).
* Totally real cubic, D <= 10^5 (4804 fields, `data/derived/porter_mendelsohn_deg3_D1e5.jsonl.xz`):
  the bound and the conjectured bound hold in every field.  log(bound/n_K) has median 38
  (min 4.18 at D = 49); the conjectured bound is tightest for D = 49, 81, 169, 148, 257 (log
  slack 2.24, 2.68, 3.85, 4.23, 4.78), i.e. for the smallest fields, as in degree 2.
