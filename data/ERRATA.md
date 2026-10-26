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
The certified data for D <= 20000 are unchanged; the conjecture is replaced by the weaker
statement in FINDINGS section 1.
