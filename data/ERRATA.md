# Errata and discrepancies

Rule (CLAUDE.md, rule 4): when computations disagree, the exact certified one wins, and
the discrepancy is recorded here.

## 2026-10: Python engine vs voronoi_unary.gp (cubic census D <= 20000)

No discrepancy.  `engine/` (exact, any degree) reproduces every one of the 832 records of
`data/legacy/census_cubic_D20000_gp.jsonl` (produced by `engine/voronoi_unary.gp`): n, E,
mu, the vertex-degree multiset and the norms of the extreme squares, with Euler
characteristic 0 and unit index 1 in both codes.

## Earlier draft of the manuscript (recorded in its Remark 5.7)

| field | earlier (box enumeration + LP, superseded) | certified |
|---|---|---|
| D = 473, x^3-x^2-10x+1 (used Z[theta], index 3) | (mu, n) = (6, 12) | (2, 4) |
| D = 6856 | mu = 3 | (mu, n) = (23, 46) |

The box-enumeration code (`perfect_unary_sail.*`) is not in this repository; see
`attic/README.md`.

## Komatsu--Watanabe reference

An earlier session wrongly declared Komatsu--Watanabe, *Perfect unary forms over certain
cubic fields*, Int. J. Number Theory 10 (2014) 1337--1342, non-existent and removed it from
the manuscript.  The paper exists (it is in the project files; verified online, DOI
10.1142/S1793042114500304).  The citation and the introduction paragraph are reinstated in
`paper/perfect-unary-cubic.tex`.

## 2026-10: a conjecture of analysis/FINDINGS.md refuted by the D <= 10^5 data

FINDINGS section 1 conjectured (from all 832 fields with D <= 20000) that every extreme
square is a sail vertex or the barycentre of a triangular sail facet at integer distance 1.
The comparison over D <= 10^5 (`analysis/sails_elements.py --census`, then `--facets`)
finds 3.3.61108.1, x^3 - x^2 - 63x - 167: the extreme square with coordinates (8, 0, -1) in
the integral basis (norm 53) lies in the triangle with vertices (1, 0, 0), (5, 1, 0),
(28, -2, -5) (distance 1, five lattice points) at barycentric coordinates (2/5, 2/5, 1/5).
The second counterexample in this range is 3.3.75425.1, of the same type; the other 121
extreme squares off the sail vertices are barycentres.  The certified data for D <= 20000
are unchanged; the conjecture is replaced by the weaker statement in FINDINGS section 1.

## 2026-10: relative minima off the sail vertices (FINDINGS section 1)

An earlier version of FINDINGS section 1 stated as CERTIFIED (D <= 20000) that the relative
minima that are not sail vertices lie in triangular facets at integer distance 1, with
barycentric coordinates (1/3, 1/3, 1/3), (1/5, 2/5, 2/5) or (1/7, 3/7, 3/7).  A complete
classification of these 220 classes (174 fields) shows German points (O. N. German, Math.
Notes 79 (2006), Thm 2) with k = 1..7, and 17 barycentres on faces at integer distance 3
(listed in FINDINGS).  The distance-1 claim had not been checked as stated:
`analysis/sails_elements.py` collects barycentres of distance-1 triangles only.  The number
710 of 832 ("relative minima = vertices and barycentres") referred to that collection and
also required every collected barycentre to be a relative minimum; the corrected statement
(non-vertex relative minima all have k = 1) holds in 763 fields.

## 2026-10: Gunnells--Yasaki, Table 1, D = 184

Gunnells--Yasaki, *Perfect forms over totally real number fields* (arXiv:0908.3130), Table 1,
gives N_D = 4306 GL_2(O)-classes of perfect binary forms over Q(sqrt 46) (D = 184, h = 1).
`engine/binary.py` finds 4308 classes, with both certificates (edge balance and the mass
formula).  An independent recheck (`analysis/binary_verify.py 46`, run twice, including once after a
machine restart) reproduced 4308. That recheck grouped the classes by invariants that do
not use qfisom: N(det A), the number of minimal vectors, |Aut A|, the denominator, and
vector counts of the trace forms. It then ran qfisom in both directions on all 5400 pairs
the invariants do not separate: no two classes are isometric.  Every class is an exactly
perfect form of minimum 1, so there are at least 4308 inequivalent classes, and the
certificates indicate the list is complete.  All other 39 entries of Table 1 agree with
this engine.  The certified count 4308 is used (rule 4).

