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
python tests/test_engine.py            # 10 checks, ~2 min;  --slow adds the 832-field census
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
| real quadratic, d = (2k+1)^2+1 squarefree (Yasaki) | n = mu = 1 |
| regulator of 3.3.49.1 vs LMFDB | 0.525454682123, agrees to 1e-10 |
| degree 4 and 5 examples | Euler 0, unit index 1 |
| relative minima | extreme squares are relative minima, mu <= ell, unit index 1 |
| gamma_K | Q(sqrt 3): gamma^2 = 16 (Leibak's bound attained); D = 1708: gamma^3 = 2D |

Not yet covered: complex cubic fields (Komatsu--Watanabe, plan E2: needs non-totally-real
signatures), Yasaki's tables and the Leibak--Porter--Ling list (reference data not in the
project files), binary forms (plan H1).

## Data

| file | content |
|---|---|
| `fields/deg3_D1e6.jsonl.xz` | the 54,600 totally real cubic fields with D <= 10^6 (LMFDB labels; identical to nflist) |
| `fields/deg2_D4e5.jsonl.xz`, `deg4_D1e6`, `deg5_D1e7` | real quadratic (D <= 4*10^5), quartic (D <= 10^6, 13,073), quintic (D <= 10^7, 9,494) |
| `data/census/deg3_D2e4.jsonl.xz` | certified records (with forms and relative minima) of all 832 cubic fields with D <= 20000 |
| `data/legacy/census_cubic_D20000_gp.jsonl` | output of `voronoi_unary.gp` on the same fields |
| `data/ERRATA.md` | discrepancies between computations, and corrections |

## Status against the plan

| item | status |
|---|---|
| P0: KW citation, skeleton, port of the GP code, cubic tests | done (`paper/` fixed; `attic/` holds nothing yet) |
| E1: degree m (maximal orders) | done, tested m = 2..5 |
| E3: relative minima | done (`--relmin`), cross-checks mu <= ell |
| E4: eutaxy, gamma_K | done, exact |
| E2: non-totally-real signatures | not started |
| B2: cubic census | D <= 20000 done; larger ranges next |
| H1: binary forms | planned (additive normalisation first) |

## Layout

```
engine/      field.py (orders, exact arithmetic, exact signs), qfmin.py (exact LLL + Fincke-Pohst),
             polyface.py (face lattices), unary.py (Voronoi), certify.py, extreme.py (eutaxy, gamma_K),
             relmin.py (relative minima), record.py, fieldlists.py, summary.py, cli.py,
             voronoi_unary.gp (reference cubic GP code)
fields/      field lists        data/      results        analysis/  scripts reading data/
docs/        THEORY.md          paper/     manuscript     attic/     superseded code
tests/       test_engine.py     PLAN.md    research plan  CLAUDE.md  rules for Claude sessions
```
