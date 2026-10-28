# Theory behind the engine

Labels: **PROVED** (a proof is given here or in `paper/perfect-unary-cubic.tex`, cited
as [M]), **CERTIFIED** (a certified computation), **CONJECTURE** (heuristic or
conjecture).  References to [M] use its numbering (theorem counter per section).

## 1. Coordinates

K is totally real of degree m, O = O_K with Z-basis w_1 = 1, ..., w_m (PARI's
`nf.zk`).  A unary form is a totally positive alpha in K ⊗ R, with value
q_alpha(x) = Tr(alpha x^2) on x in O; its Gram matrix Tr(alpha w_j w_k) is rational
when alpha is in K.  For x in O let r(x) be the integer vector with
Tr(alpha x^2) = r(x) · alpha (alpha in w-coordinates); r(x) = T sq(x) with T the
trace matrix and sq(x) the coordinates of x^2.

* **PROVED** ([M, Lemma 5.2]).  A perfect form of minimum 1 lies in K: it solves the
  rational system r(x_i) · alpha = 1 for m minimal vectors with independent r(x_i).
* **PROVED** ([M, Lemma 5.2, Lemma 3.3(b)]).  alpha and alpha' are equivalent (same orbit of
  O^x, homothety fixed by minimum 1) iff alpha' = alpha eta^2 with eta in O^x.

## 2. Voronoi's algorithm (engine/unary.py)

1. Ascent: from alpha = 1/m (minimum 1, minimal vectors +-1) move along beta with
   r(x) · beta = 0 on the current minimal vectors until a new vector reaches value 1;
   repeat until the r(x) span R^m.
2. At a perfect form alpha the Voronoi cone is cone{r(x) : x in S(alpha)}; its facets
   (polyface.py, exact) are the edges of K_1(m_O) at alpha.  For the facet F with
   inward normal beta, the neighbour is alpha + rho beta with
   rho = min over x with r(x)·beta < 0 of (r(x)·alpha - 1)/(-r(x)·beta), found by
   the classical iteration (walk): any t with alpha + t beta positive definite and
   minimum < 1 gives the upper bound min over its vectors below 1; iterate until the
   minimum is exactly 1.  Positivity is decided by exact Gauss reduction, the
   minimum by the integer Fincke--Pohst of qfmin.py.  A float estimate of the
   positivity boundary only chooses the starting t.
3. Identification: the key of alpha' is the HNF of the Z-module alpha' O (exact;
   equal keys iff alpha'/alpha in O^x).  For each known class with the same key, take
   x in F (a minimal vector of alpha'); if alpha' = alpha_j eta^2 then eta x is (up to
   sign) a minimal vector y of alpha_j, so eta is among the quotients y/x.  Each
   candidate is tested exactly (eta in O, N(eta) = +-1, alpha_j eta^2 = alpha').
   **PROVED**: the test is complete and correct.
4. New classes are stored after size reduction modulo unit squares (rounding in a
   log-space LLL-reduced unit basis: heuristic, affects heights only).

## 3. Certificates (engine/certify.py)

* **PROVED** ([M, Prop. 5.3]).  The side-pairing units (eta of step 3, over all
  classes and facets) generate O^x/{+-1}.  Hence index 1 of the subgroup they
  generate in bnfinit's unit group proves bnfinit's units fundamental without GRH.
* **PROVED** ([M, Prop. 7.1]).  dK_1(m_O)/O^x is an (m-1)-torus with a finite regular
  CW structure, so the alternating sum of its face numbers is 0.
* Face numbers are computed by gluing the local face lattices: the local face G of
  the cone of class i, G inside the facet F, is identified with the face eta G of class
  j when the neighbour across F is alpha_j eta^2.  **PROVED**: these identifications
  generate all of them, since the dual complex is a polyhedral manifold, in which the
  cells around a face of codimension >= 2 are connected through facets containing it
  (strong connectivity of the link, a sphere of dimension >= 1), and a facet lies in
  exactly two top cells.  Every gluing map is checked to send faces to faces of the
  same dimension; failure aborts the run.
* Independent count of mu: the number of distinct ideals xO over the extreme squares
  x must equal f_0 of the glued complex.

## 4. Eutaxy, extremality, gamma_K (engine/extreme.py)

gamma(alpha) = m(alpha)/N(alpha)^(1/m) (Leibak's additive normalisation, scale invariant);
gamma_K = sup gamma.  Write <beta, y> = Tr(beta y).  alpha (minimum 1) is *eutactic*
if alpha^{-1} = sum_{x in S} rho_x x^2 with all rho_x > 0.

**PROVED** (Voronoi's theorem, unary additive case).  alpha is a strict local maximum
of gamma modulo scaling iff alpha is perfect and eutactic; a form that is not both is
not even a local maximum.

*Proof.*  Note sum rho_x = sum rho_x <alpha, x^2> = <alpha, alpha^{-1}> = m.
(<=) For small beta, m(alpha + beta) = 1 + e(beta) with e(beta) = min_S <beta, x^2>,
and log N(alpha + beta) = log N(alpha) + <beta, alpha^{-1}> + O(|beta|^2).  So
log gamma(alpha + beta) - log gamma(alpha) <= phi(beta) + O(|beta|^2) with
phi(beta) = e(beta) - (1/m) <beta, alpha^{-1}>.  By eutaxy
<beta, alpha^{-1}> = sum rho_x <beta, x^2> >= m e(beta), so phi <= 0, with equality
iff <beta, x^2> = e(beta) for all x in S, i.e. (perfection) beta in R alpha.  phi is
piecewise linear, so on a complement of R alpha it is <= -c |beta| with c > 0, which
beats the O(|beta|^2) term.
(=>) If alpha^{-1} is not in cone{x^2}, separate: <beta, x^2> >= 0 on S and
<beta, alpha^{-1}> < 0; along alpha + t beta the minimum stays >= 1 and N decreases
to first order.  If alpha^{-1} is in the cone but only on its relative boundary, or
alpha is not perfect, there is beta not in R alpha with <beta, x^2> >= 0 on S and
<beta, alpha^{-1}> = 0 (a supporting functional of the face containing alpha^{-1},
or a functional vanishing on all x^2); then the minimum stays >= 1 and
log N(alpha + t beta) = log N(alpha) - (t^2/2) sum (sigma_i(beta)/sigma_i(alpha))^2
+ O(t^3) decreases.  In all cases gamma strictly increases.  []

For a perfect form the cone is full dimensional with facet normals beta_F, so eutaxy
is the exact test <beta_F, alpha^{-1}> > 0 for all facets F.

*Literature (checked against the text).*  Leibak's thesis (2006), Thm 2.1 (= [Lei05a,
Thm 5]): an extreme Humbert tuple is perfect, and if its rank is > 1 it is *weakly*
eutactic.  Weak eutaxy (thesis, Def. 2.2) allows the coefficients rho to depend on the
place, so it is vacuous for unary forms (the abstract's footnote: the claim for weak eutaxy
holds iff the number of variables is at least 2).  The unary equivalence above, with
eutaxy in the strong sense (one coefficient rho_x per minimal vector, as for Coulangeon),
is therefore not in the thesis; the proof above is self-contained.  Thesis Thm 2.3: for
totally real K a critical unary form exists.

**PROVED.**  gamma_K = max over the perfect classes of N(alpha)^(-1/m), so
gamma_K^m = 1/min N(alpha) is rational.  Indeed gamma = N^(-1/m) on dK_1(m_O), a
compact torus modulo units; on each face -log N is convex, so gamma is convex there
and its maximum over a face is attained at a vertex.  The maximiser is extreme, so it
must be eutactic (checked in every record).

**PROVED** (Leibak's bound).  det(Tr(alpha w_j w_k)) = N(alpha) disc(O), so the Hermite
invariant of the Z-lattice (O, q_alpha) is gamma(alpha)/D^(1/m) <= gamma_m, i.e.
gamma_K^m <= gamma_m^m D, with equality iff some trace lattice is similar to a
critical lattice of dimension m (A_2, A_3, D_4, ...).  The records store the exact
ratio `leibak_ratio` = gamma_K^m / (gamma_m^m D).  This is Leibak's bound (thesis Thm 2.2
= [Lei05a, Thm 6]: gamma_{m,K} <= gamma_{m[K:Q]} |disc K|^{1/[K:Q]}), which he notes is
attained for gamma_{1,Q(sqrt 3)}, gamma_{2,Q(sqrt 2)}, gamma_{2,Q(sqrt 3)} and gamma_{2,Q(sqrt 6)};
our computations agree (unary Q(sqrt 3): ratio 1; binary: gamma^4 = (2D)^2 for d = 2, 3, 6).

**PROVED** (equality in Leibak's bound, unary forms).  Let 2 <= m <= 8 and let Lambda_m be
the critical lattice of dimension m (A_2, A_3, D_4, D_5, E_6, E_7, E_8, unique up to
similarity), scaled to be even with minimum 2; det Lambda_m = d_m = 3, 4, 4, 4, 3, 2, 1.
Let d be the different of O = O_K, so d^-1 = {beta : Tr(beta O) in Z}.

(i) gamma_K^m = gamma_m^m D iff there is beta >> 0 in d^-1 with N(beta d) = d_m and
    Tr(beta x^2) in 2Z for all x in O.
(ii) In that case 2 ramifies in K, and O/beta d is isomorphic to Lambda_m^#/Lambda_m.
(iii) m = 3: in that case 2O = p q^2 and beta d = p^2 (p, q of degree 1).
(iv) m = 2: gamma_K^2 = (4/3) D iff O has a Z-basis (u, v) with N(u) + N(v) + N(u - v) = 0.
     Then D = 4(A^2 + AC + C^2) with A = N(u), C = N(v), so D = 12 mod 16, i.e.
     K = Q(sqrt d) with d = 3 mod 4, and 3 is not inert.

*Proof.*  (i) Equality holds iff some q_alpha attains the Hermite inequality, i.e. the
lattice (O, Tr(alpha x y)) is similar to Lambda_m.  Its Gram matrix in an integral basis
is rational and that of Lambda_m is primitive with minimum 2, so the similarity factor
lambda = min(q_alpha)/2 is rational, and beta = alpha/lambda gives an even integral lattice
(O, Tr(beta x y)) of determinant N(beta) D = N(beta d) = d_m; integrality means beta in
d^-1.  Conversely, an even positive lattice L of dimension m and determinant d_m has
min L >= 2 and min L <= gamma_m d_m^(1/m) = 2 (Lambda_m has minimum 2), so L is critical.
(ii) If 2 is unramified, O/2O is reduced, so the squares span O modulo 2O, and
z -> Tr(beta z) mod 2 vanishes on O: beta/2 in d^-1.  Then N((beta/2) d) = d_m / 2^m is an
integer, which fails for 2 <= m <= 8.  The dual of L = (O, Tr(beta x y)) is beta^-1 d^-1, so
L^#/L = beta^-1 d^-1 / O, isomorphic to O/beta d.
(iii) A_3^#/A_3 is cyclic of order 4.  An ideal I of norm 4 with O/I cyclic is p^2 with p
of degree 1 and 2 not in p^2, i.e. p unramified; with (ii), 2O = p q^2.
(iv) For a basis (u, v) the forms Q_a = a_1 L_1^2 + a_2 L_2^2, L_i(x, y) = sigma_i(x u + y v),
a_i > 0, are exactly the forms Tr(alpha (x u + y v)^2), alpha >> 0.  They are the positive
forms apolar to F = L_1 L_2 = N(x u + y v) = A x^2 + B x y + C y^2 (the apolarity condition
a C - b B/2 + c A = 0 for Q = a x^2 + b x y + c y^2 is linear and holds for L_1^2 and L_2^2), so
their roots tau in the upper half plane fill the geodesic A |tau|^2 + B Re tau + C = 0.  The
hexagonal form x^2 + x y + y^2 has root rho = e^(2 pi i/3), and rho lies on that geodesic iff
A - B/2 + C = 0, i.e. N(u - v) = A - B + C = -(A + C).  Every basis change of O moves rho
through SL_2(Z) rho, so some q_alpha is hexagonal iff some basis satisfies the condition.
The discriminant of the norm form is D = B^2 - 4AC = 4(A^2 + AC + C^2); A^2 + AC + C^2 is
never 2 mod 4, so d = D/4 = 3 mod 4, and by (ii) the ideal beta d has norm 3.  []

Examples: Q(sqrt 3) with (u, v) = (1, 2 + sqrt 3) (norms 1, 1, -2); Q(sqrt 7) with
(1, 3 + sqrt 7) (norms 1, 2, -3).  Leibak's thesis (Thm 2.2) records only that the bound
is attained for Q(sqrt 3) (unary) and a few binary cases; the characterisation above is
not in the papers in the project files.  `analysis/leibak_equality.py` checks (i)-(iv)
against the exact `leibak_ratio` of the censuses (FINDINGS section 7).

## 5. Relative minima (engine/relmin.py)

x in O \ 0 is a relative minimum if no y in O \ {0, +-x} has |sigma_i(y)| <= |sigma_i(x)|
for all i [M, Def. 6.2]; ties |sigma_i(y)| = |sigma_i(x)| force y = +-x.  Minima x, y
are *adjacent* if B(x v y) = {z : |z_i| <= max(|x_i|, |y_i|)} contains no lattice point
besides 0, +-x, +-y.

**PROVED.**  The adjacency graph is connected.  For a in R^m_{>0} let x(a) minimise
max_i |z_i|/a_i over z in O \ 0.  Every minimum x is x(a) on an open set of a (near
a = |sigma(x)|), and along a generic path of a's the minimiser changes only where two
minima tie with no third lattice point in the common box, i.e. between adjacent minima.

**PROVED.**  Neighbours of x lie in the box prod [-b_k, b_k], b_k = |sigma_k(N_k(x))|,
where N_k(x) minimises |sigma_k| subject to |sigma_j| < |sigma_j(x)| (j != k), and each
N_k(x) is adjacent to x.  (If y is adjacent to x and |y_k| > b_k, N_k(x) lies in
B(x v y).)

All lattice points of a box are enumerated exactly as {y : Tr(beta y^2) <= R} with beta
totally positive (checked exactly) and R a rational upper bound of the box's value;
comparisons of |sigma_i| are exact (certified double precision with a 700-fold margin
over the worst-case rounding error, otherwise Sturm sequences on an isolating
interval of the root).  The walk modulo units uses the exact key HNF(xO).

**PROVED** ([M, Lemma 6.3]).  Every extreme square is a relative minimum, so mu <= ell;
each record with `--relmin` checks this class by class.

## 6. Exactness statement

All decisions (minimum, minimal vectors, positivity, facets, equivalence, units,
eutaxy, signs) use exact integer/rational arithmetic, PARI's exact routines (HNF,
`nfroots`, `polsturm`, `bnfisunit` on given units) or floating point with a proved
error bound.  Floating point is otherwise used only to choose starting points and
representatives.  PARI's `bnfinit` units are used for size reduction and for the unit
index certificate only.

## 7. Indecomposables of totally real cubic fields

K totally real cubic, O = O_K, D = disc K, O+ the totally positive elements, U = O^x, U+
the totally positive units, U^2 the squares of units, O^v = {beta : Tr(beta O) in Z} the
codifferent (N(O^v) = 1/D).  x in O+ is *indecomposable* if it is not y + z with
y, z in O+; Ind is the set of indecomposables (1 in Ind).  P = conv(O+) is the Klein
polyhedron of the totally positive orthant; U+ acts freely on its boundary (the sail)
with quotient a 2-torus with finitely many cells (manuscript Prop. 7.1).  Every facet F
is a compact polygon with an outer functional L_F(x) = Tr(gamma_F x), gamma_F >> 0 in K,
L_F >= 1 on P and L_F = 1 on F.  The cone determinant of a triangle conv(v1, v2, v3)
with v_i in O is [O : Z v1 + Z v2 + Z v3].

**PROVED** (Lemma 7.1).  (a) Totally positive units are indecomposable.  (b) Sail vertices
are indecomposable.  (c) If x in O+ and L_F(x) < 2 for some facet F, then x is
indecomposable.  *Proof.*  (a) y < u coordinatewise forces N(y) < 1.  (b) If v = y + z,
then y and y + 2z lie in P (z is in the recession cone) and v is their midpoint.
(c) L_F(y + z) = L_F(y) + L_F(z) >= 2.  []

**PROVED** (Lemma 7.2, integral perfect forms with unit minimal vectors).  Let alpha >> 0
lie in O^v, with Tr(alpha x^2) >= 1 for all x in O \ 0 and Tr(alpha u_i^2) = 1 for three
units u_1, u_2, u_3 with distinct squares.  Then Tr(alpha u_j u_k) = delta_jk, {u_1, u_2, u_3}
is a Z-basis of O, and N(alpha) D = 1, i.e. alpha O = O^v.  Consequently any two such
alpha differ by a totally positive unit.
*Proof.*  b(x, y) = Tr(alpha x y) is integral on O (xy in O, alpha in O^v), and its
determinant on O is N(alpha) D, a positive integer.  For j != k, u_j +- u_k != 0, so
1 <= b(u_j +- u_k, u_j +- u_k) = 2 +- 2 b(u_j, u_k), hence |b(u_j, u_k)| <= 1/2, hence
b(u_j, u_k) = 0.  The u_j are orthonormal, so independent, and L = Z u_1 + Z u_2 + Z u_3 has
determinant 1 = [O : L]^2 N(alpha) D.  Both factors are positive integers.  If alpha' is
another one, alpha' O = alpha O, so alpha'/alpha is a unit, totally positive.  []

In the census language: an integral perfect form whose minimal vectors are units is
the cubic lattice I_3 with an orthonormal basis of units; there are at most [U+ : U^2]
classes of them.  `analysis/one_square_class.py` checks this on the cubic census
(FINDINGS section 9).

**PROVED** (Lemma 7.3, unit-square facets).  Let F = conv(u_1^2, u_2^2, u_3^2) be a
triangular facet of P, u_i units, with cone determinant 1.  Then gamma_F satisfies the
hypotheses of Lemma 7.2: L_F = 1 on the Z-basis u_1^2, u_2^2, u_3^2 of O, so L_F(O) is in
Z, i.e. gamma_F in O^v; and Tr(gamma_F x^2) = L_F(x^2) >= 1 because x^2 in O+ lies in P.
Hence gamma_F O = O^v, and all such facets lie in one U+-orbit: if gamma_F' = eta gamma_F
with eta in U+, then F' = eta^-1 F, since F_{eta gamma} = {y in P : Tr(gamma eta y) = 1}
= eta^-1 F_gamma.

**PROVED** (Theorem 7.4).  In every totally real cubic field some indecomposable is not
the square of a unit; i.e. Ind / U^2 has at least two elements.
*Proof.*  Suppose Ind = U^2 (one class, the class of 1).
1. U+ = U^2: totally positive units are indecomposable (7.1a), hence squares.
2. The vertex set of P is U+ and V = 1: vertices are indecomposable (7.1b), hence units;
   the vertex set is a non-empty U+-invariant subset of the single orbit U+.
3. Facets contain no lattice points other than their vertices: such a point p is in O+
   with L_F(p) = 1, so indecomposable (7.1c), and not a unit, since units are vertices.
4. Euler on the sail torus: V - E + F = 0 with V = 1 and 2E = sum_k k F_k (F_k facet
   classes with k vertices) gives sum_k (k - 2) F_k = 2: two triangle classes, or one
   quadrilateral class.
5. No quadrilateral.  An empty lattice quadrilateral (step 3; the lattice points of the
   facet plane form a 2-dimensional affine lattice) is a parallelogram: by Pick it has
   area 1, its two triangles across a diagonal are unimodular, and in coordinates
   p_2 = (0,0), p_1 = (1,0), p_3 = (0,1) the fourth vertex is (1 + k, 1 - k), convex only
   for k = 0.  On the torus with V = 1, E = 2, F = 1 the gluing is a b a^-1 b^-1
   (orientable, and units preserve orientation), so the vertices are 1, a, ab, b in
   cyclic order with a, b in U+.  Then 1 + ab = a + b, i.e. (1 - a)(1 - b) = 0: impossible.
6. Both triangles have cone determinant 1.  Otherwise the half-open parallelepiped
   {sum lambda_i v_i : 0 <= lambda_i < 1} of a triangle F = conv(v1, v2, v3) contains a
   lattice point x != 0.  It is totally positive, and so is x' = v1 + v2 + v3 - x, so
   1 <= L_F(x) = 3 - L_F(x') <= 2.  If L_F(x) < 2, x is indecomposable (7.1c) and not a
   unit: if L_F(x) > 1 it is not on the sail, and if L_F(x) = 1 it is a lattice point of F
   other than a vertex.  If L_F(x) = 2, x' is a lattice point of F other than a vertex,
   against step 3.
7. By steps 1, 2 and 6 both triangle classes satisfy Lemma 7.3 (their vertices are in
   U+ = U^2), so they lie in one U+-orbit, against step 4.  []

Remarks.  (i) The perfect forms enter through Lemma 7.2: in the situation of Theorem 7.4
the squared sail coincides with P (both are conv(U+) + R^3_{>=0}), each facet is a
perfect form, and Lemma 7.2 says it would be the integral form I_3; two classes of I_3
facets cannot coexist.  (ii) The real quadratic analogue fails, as it must (Q(sqrt 5) has
Ind = U^2): the sail modulo U+ is a circle, which needs only one edge.  (iii) "Every
indecomposable is a unit" (Ind = U+): steps 2-6 go through with U+ in place of U^2
(FINDINGS section 9: S+ = 2, two unimodular triangles conv(1, a, ab), conv(1, b, ab),
b = -a + m + n ab with m + n >= 3 by convexity), and Theorem 7.4 adds that U+ != U^2, so
the units do not have all signatures.  This version is open.  (iv) Degree >= 4: Lemmas
7.2 and 7.3 hold verbatim (Gram matrix I_m, alpha O = O^v), but steps 4-6 do not: empty
lattice simplices of dimension >= 3 need not be unimodular, the parallelepiped argument
only gives 1 <= L_F(x) <= m - 1, and a torus of dimension >= 3 can have a single cell.

**PROVED** (Lemma 7.5, minimal triples; any order O of a totally real cubic field).  A
minimal triple {x, y, z} of the Minkowski-Voronoi complex (manuscript Prop. 7.7) is either
a Z-basis of O or, after changing signs, z = x + y.
*Proof.*  Let m_i = max of |x_i|, |y_i|, |z_i|; each of x, y, z attains a different m_i,
say |x_1| = m_1, |y_2| = m_2, |z_3| = m_3, and all other entries are smaller in absolute
value (equality |a_i| = |b_i| forces a = +-b).  The closed box B = prod [-m_i, m_i]
contains no lattice point other than 0, +-x, +-y, +-z.
Independent case.  Let L = Z x + Z y + Z z.  Minkowski: 8 m_1 m_2 m_3 <= 8 covol(O), and
|det(x, y, z)| <= 4 m_1 m_2 m_3 (the determinant is affine in each entry, so on the box
of matrices with |entries| <= m it is maximal at sign matrices, where it is at most 4),
so [O : L] <= 4.  O/L has no element of order 2: such an element is represented by
(a x + b y + c z)/2 with a, b, c in {0, 1} not all 0, with any signs.  One or two nonzero
terms already lie strictly inside B.  With three terms, coordinate i of x + s y + t z
reaches 2 m_i only if the three terms in it have one sign, which fixes (s, t) uniquely
for each i.  So of the four sign choices at least one has all |coordinates| < 2 m_i, and
half of it is a lattice point strictly inside B.  O/L has no element of order 3: it is
represented by (a x + b y + c z)/3 with a, b, c in {-1, 0, 1}, and
|a x_i + b y_i + c z_i| < 3 m_i.  So [O : L] is prime to 6 and at most 4, hence 1.
Dependent case.  x, y are adjacent, so (x +- y)/2 is not in O (it lies strictly inside
B(x v y)), and {x, y} is a basis of O cap (Q x + Q y) by the same argument in dimension 2.
Write z = p x + q y with p, q in Z; changing signs of x, y, p >= 0 and q >= 0.  q = 0 or
p = 0 makes z a multiple of x or y, not a relative minimum unless z = +-x, +-y.  If
p, q >= 1 and (p, q) != (1, 1), then x + y = (q - 1)/(p + q - 1) x + (p - 1)/(p + q - 1) y +
1/(p + q - 1) z is a convex combination with at least two positive weights, hence strictly
inside B, so x + y in {+-x, +-y, +-z}: impossible.  So z = x + y.  []

`analysis/minkowski_voronoi.py` checks |det(x, y, z)| in {0, 1} for every minimal triple
and counts the degenerate classes: for the 832 fields with D <= 20000, 12,903 of the 50,176
classes are degenerate (at most half in each field; none only for D = 81).  In each of the 26 fields with ell = 1 (D <= 10^6) one of the two triple classes is
a basis and the other is degenerate, z = +-x +- y with x, y, z units (through 1 it is
{1, eps, 1 - eps}, eps = N_k(1)); D = 81 (ell = 2) has four basis classes and no degenerate
one.

**PROVED** (Lemma 7.6, sign pattern).  Let eps be a unit with |eps_1| > 1 > |eps_2|, |eps_3|,
signs (+, +, -), such that 1 - eps is also a unit.  Then eps is a root of
x^3 - (t + 3) x^2 + t x + 1 = P_t(x - 1), P_t(x) = x^3 - t x^2 - (t + 3) x - 1, so K is the
simplest cubic field K_t and Z[theta_t] is contained in O with theta_t = eps - 1.
*Proof.*  Write the minimal polynomial as x^3 - A x^2 + B x - C.  C = N(eps) < 0, so
C = -1; N(1 - eps) = (1 - eps_1)(1 - eps_2)(1 - eps_3) < 0 equals the value at 1, so
1 - A + B + 1 = -1, B = A - 3.  With t = A - 3 this is P_t(x - 1).  []
