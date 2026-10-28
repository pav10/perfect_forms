# perfect_forms — perfect unary forms over number fields, exactly

Exact Voronoi theory of perfect unary (Humbert) forms over totally real number fields:
the number n_K of perfect forms modulo homothety and units, the number mu_K of extreme
squares (facets of the Ryshkov polyhedron), the whole face complex of the torus
dK_1(m_O)/O^x, the additive Hermite--Humbert constant gamma_K (exactly), and the relative
minima of O as an independent cross-check.  The research plan is [`PLAN.md`](PLAN.md); the
manuscript is [`paper/perfect-unary-cubic.tex`](paper/perfect-unary-cubic.tex); the
mathematics behind every step is in [`docs/THEORY.md`](docs/THEORY.md).

* Python on top of PARI/GP (cypari2).  Any totally real degree (tested 2--5).
* Every decision is exact; floating point only chooses starting points (`docs/THEORY.md` §6).
* Every record carries certificates (Euler characteristic, unit index, an independent
  count of mu) and provenance, and is labelled `CERTIFIED` or quarantined.
* The earlier cubic GP code [`engine/voronoi_unary.gp`](engine/voronoi_unary.gp) is kept as
  the reference implementation; the new engine reproduces all 832 records of its census.

## Installation

```bash
pip install -r requirements.txt        # cypari2 (brings PARI), numpy
python tests/test_engine.py            # 18 checks, ~2 min;  --slow adds the 832-field census
```

## Quick start

```bash
python -m engine field "x^3-x^2-5*x-1" --label 3.3.404.1 --relmin
#   n 7, e 11, mu 4, f_vector [4, 11, 7], degrees {3: 6, 4: 1}, gamma_pow 1616/3, ell 5, CERTIFIED

python -m engine fields --degree 3 --disc-max 1000000 --mirror deg3.txt --out fields/deg3_D1e6.jsonl
python -m engine run --fields fields/deg3_D1e6.jsonl.xz --disc-max 20000 \
    --out data/census/deg3_D2e4.jsonl --jobs 4 --relmin --progress 100
python -m engine summary data/census/deg3_D2e4.jsonl.xz
```

## Commands

| command | purpose |
|---|---|
| `field POL [--label L] [--json] [--no-forms] [--relmin]` | one field (polynomial in x); prints the main invariants, or the full JSON record |
| `fields --degree N --disc-max X [--disc-min Y] [--mirror FILE] --out FILE` | field list as JSONL `{"label","pol","D"}`; degree 2 generated from fundamental discriminants; degree >= 3 from the LMFDB mirror of [rvisser7/indecomposables](https://github.com/rvisser7/indecomposables) (`totally_real_fields/degN.txt`), and for degree 3 checked field by field against PARI's `nflist` (C3 and S3) |
| `run --fields FILE --out OUT [--disc-min/--disc-max] [--limit K] [--jobs J] [--timeout S] [--no-forms] [--relmin] [--retry-errors] [--progress N]` | parallel batch, resumable (labels already in OUT, OUT.quarantine or OUT.errors are skipped); certified records go to OUT, the others to OUT.quarantine, failures to OUT.errors |
| `summary FILE... [--what overview\|lines]` | overview statistics, or one line per field |

## The record (one JSON object per field)

| key | meaning |
|---|---|
| `label`, `pol`, `pol_input`, `degree`, `signature`, `disc`, `order_index` | field data (`pol` in polredabs form; `order_index` 1 = maximal order) |
| `class_number`, `regulator`, `fundamental_units`, `integral_basis` | from `bnfinit` (units certified by the unit-index certificate) |
| `n`, `e`, `mu`, `f_vector` | n_K, edges, mu_K and the face vector (f_0, ..., f_{m-1}) of the dual complex: f_0 = mu (extreme squares), f_{m-1} = n (perfect forms); for m = 3, (mu, e, n) |
| `degrees` | {vertex degree: count}: number of facets of each Voronoi cone (= edges at the vertex) |
| `facet_sizes` | {number of vertices: count} for the facets of K_1(m_O) (for m = 3 the polygon sizes) |
| `extreme_square_norms` | \|N(x)\| for one extreme square x per class |
| `gamma_pow`, `gamma`, `gamma_argmax` | gamma_K^m (exact rational), gamma_K, the maximising classes |
| `n_extreme`, `n_weakly_eutactic` | number of perfect classes that are eutactic (= extreme), and on the boundary of the cone |
| `leibak_ratio` | gamma_K^m / (gamma_m^m D), exact; 1 iff Leibak's bound is attained |
| `ell` | (with `--relmin`) number of classes of relative minima |
| `forms` | (unless `--no-forms`) per class: `alpha` (and `alpha_coords` in the integral basis), `norm` N(alpha), `degree`, `eutaxy`, `minvecs` (coordinates), `extreme` (vertex of the Voronoi polytope), `minvec_norms`, `neighbours` (class indices across the facets) |
| `certificates` | `euler` (must be 0), `unit_index` (must be 1), `mu_keys` (independent count of mu), `edges_ok`, and with `--relmin`: `relmin` {`ell`, `extreme_squares_are_relmin`, `relmin_unit_index`} |
| `status` | `CERTIFIED` or `QUARANTINE` |
| `stats`, `provenance` | counters; git hash, PARI/cypari2/Python versions, precision, seed, seconds |

## Validation (`tests/test_engine.py`, plan §7)

| check | result |
|---|---|
| exact enumerator vs brute force (random rational forms, dim 2--4) | identical |
| cubic D = 49, 257, 404 (degrees 3^6 4^1), 473, 1708, 6856, 6997, 16689 | (mu, n) = (1,2), (2,4), (4,7), (2,4), (12,21), (23,46), (13,26), (77,152) |
| Ennola a = 2..12, Shanks a = -1..11 | mu = 1 except Ennola a = 3 (mu = 2); discriminants as in [M, Obs. 5.4] |
| D = 676829 (x^3-x^2-221x+29) | (mu, n) = (136, 272), all simplicial |
| census D <= 20000 (`--slow`) | 832 fields, 27 with mu = 1, identical to the GP census |
| real quadratic, d = (2k+1)^2+1 squarefree (Yasaki, proof of Cor. 5.1) | n = mu = 1 |
| Q(sqrt 19), Q(sqrt 23) (Yasaki, Sec. 5.4) | 4 and 2 classes; each listed set of minimal vectors is that of exactly one class |
| regulator of 3.3.49.1 vs LMFDB | 0.525454682123, agrees to 1e-10 |
| degree 4 and 5 examples | Euler 0, unit index 1 |
| relative minima | extreme squares are relative minima, mu <= ell, unit index 1 |
| gamma_K | Q(sqrt 3): gamma^2 = 16 (Leibak's bound attained); D = 1708: gamma^3 = 2D |
| equality in Leibak's bound (docs/THEORY.md 4) | criterion = exact gamma_K for all real quadratic D < 1200 (19 equalities) and five cubic fields |

| binary forms over Q(sqrt d) vs Gunnells--Yasaki Table 1 (`--slow`: all D <= 44) | identical class counts (2, 2, 3, 9, 34, 16, 22, 30, 39, 149, 97, 52, 144, 66); Leibak's 22 classes over Q(sqrt 6); edge balance and mass formula |
| real quadratic, d < 400: Leibak--Porter--Ling Thm 1.2, Yasaki Thm 5.1 and Cor. 3.1 | n_K = 1 exactly for the types T1--T4; l_d/2 <= n_K <= l_d; 2 or 3 minimal vectors per form |
| simplest cubics with monogenic O_K (Leibak--Porter--Ling Thm 1.5) | n_K = 2 |
| complex cubic v^3 + l v - 1 (Komatsu--Watanabe Thms 3.1--3.2), l = 2, 4, 5, 7, 10 | n_K = 1, the form is a_0, S(a_0) = {+-1, +-v}, unit index 1 |

Checked outside the test suite (needs the census data): Yasaki's global count, 10,732,735
perfect forms over the 121,580 fields Q(sqrt d), 2 <= d <= 200000 squarefree (JTNB 25
(2013), Sec. 4), is reproduced exactly (`analysis/yasaki_count.py`, analysis/FINDINGS.md
section 5).  Not covered: the multiplicative (Pohst--Wagner) binary benchmark, which
belongs to pav10/hermite-humbert.

## Binary forms (additive normalisation, plan H1)

`engine/binary.py`: perfect binary Humbert forms Tr(v^T A v), A in Sym_2(K), over the
maximal order of a totally real field, modulo GL_2(O_K) and homothety (Koecher/Leibak
additive normalisation; the multiplicative Icaza/BCIO normalisation of pav10/hermite-humbert
is a different theory and is never mixed with this one).

```python
from engine.field import Order
from engine.binary import BinaryVoronoi, edge_balance, mass_formula
V = BinaryVoronoi(Order('x^2-x-3')).run()        # Q(sqrt 13)
len(V.classes)                                   # 9 perfect classes mod GL_2(O_K)
edge_balance(V), mass_formula(V)[0]              # [], 0   (certificates)
V.summary()                                      # per class: minimal vectors, |Aut|, N(det A), gamma, eutaxy
```

Exact equivalence and stabilisers: PARI's `qfisom`/`qfauto` on the Z-forms
Tr(theta^k v^T A v), k < m (theta a totally positive generator), which forces O_K-linearity.
Certificates: edge balance N_ij/|Aut A_i| = N_ji/|Aut A_j|, and the mass formula
sum over interior cells mod GL_2(O_K) of (-1)^dim/|Stab| = chi_orb(GL_2(O_K)) = 0.
Results so far are in [`analysis/FINDINGS.md`](analysis/FINDINGS.md) §3.

## Analysis

* [`analysis/FINDINGS.md`](analysis/FINDINGS.md): labelled results (PROVED / CERTIFIED /
  CONJECTURE), in particular the comparison of perfect forms with the sails of pav10/sails.
* `analysis/sails_link.py CENSUS SAILS_SUMMARY.xz`: counts vs the sails summary data.
* `analysis/sails_elements.py FIELDS DISC_MAX [OUT] [--jobs J]`: element-level comparison
  (extreme squares, relative minima, sail vertices, modulo units); imports pav10/sails
  read-only from `SAILS_PATH` (default `../sails`) without writing bytecode there.
  `--census CENSUS OUT [JOBS]` reads the extreme squares from census records instead;
  `--facets CENSUS CENSUS_OUT OUT` lists the facets (shape, distance, exact barycentric
  coordinates) of the extreme squares that are not sail vertices.
* `analysis/minkowski_voronoi.py CENSUS [LIMIT] [OUT] [JOBS]`: cells (V, E, F) of the
  Minkowski--Voronoi complex modulo units (relative minima, adjacent pairs, minimal triples),
  and `V_degenerate`, the classes of linearly dependent minimal triples (every other minimal
  triple is checked to be a Z-basis of O, docs/THEORY.md Lemma 7.5).
* `analysis/cubic_census.py CENSUS...`: census statistics (Markdown).
* `analysis/binary_census.py D_MIN D_MAX OUT [--jobs J]`: certified census of perfect binary
  forms over Q(sqrt d).
* `analysis/binary_verify.py d`: independent recheck of one binary census field (classes
  separated by qfisom-free invariants, qfisom both ways on the rest).
* `analysis/binary_compare.py CENSUS [HH_TSV]`: a binary census against Gunnells--Yasaki
  Table 1, Leibak's bound gamma^4 <= 4 D^2, and the multiplicative constant of
  pav10/hermite-humbert (read-only; gamma_add >= 2 sqrt(gamma_mult) by AM-GM).
* `analysis/recheck.py CENSUS REPORT [--sample K] [--replace-dirty FIXED] [--jobs J]`:
  recompute records and compare (provenance repair, spot checks).
* `analysis/porter_mendelsohn.py CENSUS [OUT]`: n_K against the Porter--Mendelsohn bound
  (covering radius numerical).
* `analysis/leibak_equality.py CENSUS...`: the exact criterion for equality in Leibak's
  bound (docs/THEORY.md section 4) against the census values.
* `analysis/ell_one.py FIELDS DMAX OUT JOBS`: whether all relative minima are units
  (ell = 1); see analysis/FINDINGS.md section 9 (indecomposables).
* `analysis/one_square_class.py CENSUS`: the integral perfect forms (in the codifferent)
  whose minimal vectors are units, with the checks of docs/THEORY.md Lemma 7.2 (Gram matrix
  I_3, Z-basis, N(alpha) D = 1, at most [U+ : U^2] classes); FINDINGS section 9.
* `analysis/yasaki_count.py CENSUS...`: Yasaki's global count (10,732,735 perfect forms over
  the 121,580 fields Q(sqrt d), d <= 200000 squarefree) from certified records, with his
  Thm 5.1 and Leibak--Porter--Ling Thm 1.2 on the whole range.

## Data

| file | content |
|---|---|
| `fields/deg3_D1e6.jsonl.xz` | the 54,600 totally real cubic fields with D <= 10^6 (LMFDB labels; identical to nflist) |
| `fields/deg2_D4e5.jsonl.xz`, `deg4_D1e6`, `deg5_D1e7` | real quadratic (D <= 4*10^5), quartic (D <= 10^6, 13,073), quintic (D <= 10^7, 9,494) |
| `data/census/deg3_D2e4.jsonl.xz` | certified records (with forms and relative minima) of all 832 cubic fields with D <= 20000 |
| `data/census/deg3_D1e5.jsonl.xz` | certified records (with forms and relative minima) of all 4804 cubic fields with D <= 10^5; provenance re-checked, see `deg3_D1e5.recheck.json` |
| `data/census/deg2_D1e5.jsonl.xz` | certified records (no forms; ell from the continued fraction) of all 30,394 real quadratic fields with D <= 10^5 |
| `data/census/binary_deg2.jsonl.xz` | certified binary census (perfect binary forms mod GL_2(O_K), additive normalisation) for all squarefree d <= 100 and d = 105, 113, 129, 141; per class: A, |Aut|, minimal vectors, gamma, eutaxy |
| `data/census/deg2_d2e5_rest.jsonl.xz` | the other 91,186 fields Q(sqrt d), d <= 200000 squarefree (D up to 799,996), same format; with the previous file, Yasaki's range |
| `data/legacy/census_cubic_D20000_gp.jsonl` | output of `voronoi_unary.gp` on the same fields |
| `data/derived/sails_vs_perfect_deg3_D2e4.jsonl.xz` | per field: n, e, mu, ell, sail vertices/edges/facets (all sign classes, mod units), positions of the extreme squares |
| `data/derived/sails_vs_perfect_deg3_D1e5.jsonl.xz` | the same for all 4804 cubic fields with D <= 10^5 (`--census`; no relative minima as elements) |
| `data/derived/sails_nonvertex_extreme_deg3_D1e5.jsonl.xz` | the 123 extreme squares (D <= 10^5) that are not sail vertices: facet shape, distance, exact barycentric coordinates |
| `data/derived/minkowski_voronoi_deg3_D2e4.jsonl.xz` | (V, E, F) of the Minkowski--Voronoi complex mod units and the number `V_degenerate` of degenerate minimal triple classes, all 832 cubic fields with D <= 20000 |
| `data/derived/porter_mendelsohn_deg{2,3}_D1e5.jsonl.xz` | n_K against the Porter--Mendelsohn bound (covering radius numerical, 6 digits) |
| `data/ERRATA.md` | discrepancies between computations, and corrections |

## Status against the plan

| item | status |
|---|---|
| P0: KW citation, skeleton, port of the GP code, cubic tests | done (`paper/` fixed; `attic/` holds nothing yet) |
| E1: degree m (maximal orders) | done, tested m = 2..5 |
| E3: relative minima | done (`--relmin`), cross-checks mu <= ell |
| E4: eutaxy, gamma_K | done, exact |
| E2: non-totally-real signatures | signature (1,1) done (`engine/complexcubic.py`, Komatsu--Watanabe reproduced); (2,1), (1,2), (0,3) not started |
| B2: cubic census | D <= 10^5 done (4804 fields, analysis/FINDINGS.md); D <= 10^6 needs a go-ahead |
| B1: real quadratic census | D <= 10^5 done (30,394 fields): LPL Thm 1.2 and Yasaki Thm 5.1 confirmed in every field; d <= 200000 squarefree (121,580 fields) done: Yasaki's total of 10,732,735 perfect forms reproduced exactly |
| H1: binary forms | additive normalisation: census of all squarefree d <= 100 (126,321 classes, certified); Gunnells--Yasaki Table 1 reproduced except D = 184 (4308 vs 4306, ERRATA); rank-2 conjectures in analysis/FINDINGS.md section 3 |
| sails link (cubic) | D <= 10^5 done (analysis/FINDINGS.md §1): three hulls, Minkowski--Voronoi complex (V, E, F) = (2 ell, 3 ell, ell) |
| indecomposables (cubic) | PROVED (docs/THEORY.md §7, manuscript Thm 9.4): no totally real cubic field has only squares of units as indecomposables; "only units" and the classification of ell = 1 open (FINDINGS §9) |

## Layout

```
engine/      field.py (orders, exact arithmetic, exact signs), qfmin.py (exact LLL + Fincke-Pohst),
             polyface.py (face lattices), unary.py (Voronoi), certify.py, extreme.py (eutaxy, gamma_K),
             relmin.py (relative minima), walk.py (exact neighbour walk), binary.py (binary forms),
             complexcubic.py (signature (1,1)),
             record.py, fieldlists.py, summary.py, cli.py,
             voronoi_unary.gp (reference cubic GP code)
fields/      field lists        data/      results        analysis/  scripts reading data/
docs/        THEORY.md          paper/     manuscript     attic/     superseded code
tests/       test_engine.py     PLAN.md    research plan  CLAUDE.md  rules for Claude sessions
```
