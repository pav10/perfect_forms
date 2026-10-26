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
* **CONJECTURE.**  For every totally real cubic field, every extreme square is a vertex of
  the sail of its sign class or the barycentre of a triangular sail facet at integer
  distance 1.  (Then mu <= V + #{such facets} and n <= 2 mu.)
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

Perfect binary forms Tr(v^T A v) modulo GL_2(O_K) (engine/binary.py; certificates: edge
balance and the mass formula sum (-1)^dim/|Stab| = 0):

| d | 5 | 2 | 3 | 13 | 17 | 6 | 7 | 10 (h = 2) |
|---|---|---|---|---|---|---|---|---|
| perfect classes mod GL_2(O_K) | 2 | 2 | 3 | 9 | 34 | 22 | 30 | 52 |

Classes with equal invariants are Galois conjugates (Gal(K/Q) is not in GL_2(O_K)).  The
comparison with Ong (d = 2, 3, 5), Leibak (d = 6) and Gunnells--Yasaki is pending: their
tables are not in the project files [REF?].
