# Perfect unary forms over number fields — computational plan

Scope: a Claude Code project that (1) establishes Yasaki-type bounds for perfect unary forms over totally real cubic fields and beyond, (2) classifies small perfect and extreme forms, and (3) brings in the project's other threads (volume sandwich, higher-rank lower bounds, Hermite–Humbert constants). It is designed so that every module extends to higher degree and to rank n ≥ 2.

Notation: K a number field with signature (r1, r2), d = r1 + r2, O an order. n_K is the number of perfect unary forms modulo homothety and O^×. μ_K is the number of unit-classes of extreme squares, i.e. facets of the Ryshkov polyhedron. ℓ_O is the number of unit-classes of relative minima of O. Reg, D are the regulator and discriminant of O.

---

## 0. Starting point

**Proved, in `perfect-unary-cubic.tex`:**
- Squared-sail duality.
- Euler bound for d = 3: ½μ_K ≤ n_K ≤ 2μ_K, with the defect identity n_K = 2μ_K − Σ_v(deg v − 3).
- Regulator sandwich for totally real cubic orders: 4Reg/(log 2D)² ≤ μ_K ≤ ℓ_O ≤ 4Reg/(log 2)².
- The degree-m analogue of the sandwich.
- n_K ≥ 2mμ_K/τ_m.

**Computed and certified:** all 832 totally real cubic fields with D ≤ 20000, plus three fields with n_K up to 3534. The code is `voronoi_unary.gp`, with an exact Voronoi algorithm, an Euler check and a unit-index certificate.

**Errors to fix before anything else:**
1. **Komatsu–Watanabe is real.** "Perfect unary forms over certain cubic fields", Int. J. Number Theory 10 (2014) 1337–1342, is in the project files. It proves n_K = 1 for the complex cubic fields Q(v), v³ + ℓv − 1 = 0, with 4ℓ³ + 27 squarefree. Last turn I wrongly declared it non-existent and deleted it. Reinstate the citation and the intro paragraph; theorem numbers in the framework section can cite KW Thm 2.2 (quoting WYH Thm 2.4) and OY Thm 5.1, which are verified in KW's text.
2. **The old box-enumeration code** (`perfect_unary_sail.*`) is superseded. Move it to `attic/` and never use it as evidence.

**Open mathematical items flagged in the manuscript:**
- Corollary 6.8(c), on orders Z + pO_K, still needs an audit.
- Buchmann's period-length papers need checking for prior art on ℓ_O ≪ Reg.
- The Gunnells–Yasaki Thm 5.3 sign-genericity gap is unresolved.

---

## 1. Repository and invariants

```
engine/      exact Voronoi for unary forms (PARI/GP core + Python driver)
  unary.gp           current cubic code, generalised (§2)
  polyface.py        face lattice of dual polytopes (Sage/pplpy, exact rationals)
  certify.py         Euler characteristic, unit-index, LP extremality checks
  relmin.gp          independent relative-minima enumerator (Voronoi 1896 / Buchmann)
  binary/            rank-2 Humbert Voronoi (§6)
fields/      field lists: nflist + LMFDB cross-check, one JSONL per (degree, signature)
data/        results, one record per (field, order): schema in §1.2
analysis/    notebooks/scripts producing tables and figures; never edit data/
tests/       regression benchmarks (§7); must pass before any run is reported
paper/       LaTeX (manuscript + successors)
attic/       superseded code, kept for provenance only
CLAUDE.md    the rules in §8
```

### 1.1 Non-negotiable invariants
- **Exactness.** Every decision is made in exact arithmetic: minimality, equivalence, adjacency, extremality and eutaxy. Floating point may only widen search ranges, or serve as a pre-filter whose result is re-verified exactly.
- **Certificates.** Every record carries:
  - the alternating face count Σ(−1)^i f_i = 0;
  - the index of the side-pairing units in the unit group, which must be 1. This also certifies the unit group without GRH.

  A record failing either check is quarantined, not repaired.
- **Provenance.** Each record stores the code git hash, the PARI/Sage versions, real precision, wall time and the polynomial in `polredabs` form.
- **Three labels.** Results in tables and writing are labelled PROVED, CERTIFIED-COMPUTATION or CONJECTURE/HEURISTIC, and never mixed.

### 1.2 Record schema (JSONL)
- **Field data:** field id, polynomial, signature, D, order index, h, Reg, fundamental units.
- **Counts:** n_K, e_K, μ_K, ℓ_O.
- **Combinatorics:** the vertex-degree multiset and face-size multiset.
- **Forms:** for each perfect form, the element α ∈ K, its extreme minimal vectors, and the norms of those vectors.
- **Constants:** γ_K, the covering radius ρ_∞(Λ_K) and the certificates.
- **Provenance:** as in §1.1.

---

## 2. Engine (prerequisite for everything)

**E1. Generalise the cubic code to totally real degree m and to orders.**
- *Perfect forms.* By Leibak (Cor. 1.1 of the thesis), a perfect form is Tr(αx²) with α ∈ K totally positive.
- *Neighbours.* Use Voronoi neighbours via the facets of the dual polytope conv{q(x) : x ∈ S(a)}, an (m−1)-polytope. Its face lattice is computed exactly in `polyface.py`; the cubic code had polygons only.
- *Short vectors.* Use exact enumeration from the rational Gauss reduction after LLL, generalising the 3-variable enumerator to m variables.
- *Acceptance:* reproduces all §7 benchmarks for m = 2 and 3.

**E2. Non-totally-real signatures.** The coordinates of perfect forms now lie in the Galois closure. Represent them exactly in a number field and decide signs by certified interval arithmetic, raising precision until the interval excludes 0.
- *Acceptance:* reproduces Komatsu–Watanabe, n_K = 1 for v³ + ℓv − 1 with ℓ = 2, …, 20 and 4ℓ³ + 27 squarefree. Also handles the d = 3 signatures (2,1), (1,2), (0,3) with Euler checks.

**E3. Independent relative-minima enumerator** (`relmin.gp`). This is Voronoi's 1896 algorithm on the genuine lattice ψ(O). It gives ℓ_O, and Lemma 6.3 then predicts that every extreme square found by E1 appears among the relative minima. This cross-checks two independent codes.

**E4. Extremality and eutaxy.** A unary form is extreme if it is a local maximum of γ(a) = m(a)/N(a)^{1/m}. Implement the Voronoi-type test of Leibak: extreme ⇔ perfect and eutactic. Eutaxy means that (1/σ_i(α))_i lies in the relative interior of cone{q(x) : x ∈ S(a)}, decided by an exact LP. Then γ_K is the maximum of γ over the finitely many perfect classes, so **γ_K is computed exactly**.
- *Caveat:* first state precisely, from Leibak's thesis (Thm 2.1, Thm 2.3), which unary statement is proved there.

---

## 3. Goal 1 — Yasaki-type bounds, cubic and beyond

**B1. Baseline in degree 2.** Reproduce Yasaki's bounds in his normalisation, and the h_d = 1 classification of Leibak–Porter–Ling, for squarefree d ≤ 10⁵ using E1 with m = 2. This validates the engine against a known theorem. It also measures, for d = 2, the ratio n/ℓ and the dependence of the period on Reg.

**B2. Cubic census to large discriminant.**
- *Stages:* D ≤ 10⁵, then D ≤ 10⁶ (about 7·10⁴ fields), then D ≤ 10⁷ only if cost permits.
- *Field lists:* `nflist("S3"/"C3", [1, X], 0)`, cross-checked against LMFDB counts.
- *Outputs:* n, μ, ℓ, the degree profile, γ_K and ρ_∞ for every field.
- *Cost:* the D ≤ 20000 census ran in under a minute, and cost grows roughly linearly in n_K.

**B3. Sharp constants.**
- *Lower bound:* determine empirically inf μ/(Reg/(log D)²) and sup μ/Reg over the census, and compare with the proved 4 and 4/(log 2)². Test the refined lower bound 4Reg/(3 log γ_K)², which E4 makes computable.
- *Upper bound:* the target is a provable improvement of the constant in μ ≤ ℓ ≤ c·Reg. Use sharper packing of the separated sign classes in Lemma 6.4, and compare ℓ_O/Reg with Buchmann.

**B4. Test Porter–Mendelsohn.** They prove n_K ≤ e^{2nρ_∞}γ_K^n/|Δ|^{1/2}, and conjecture the improvement to /|Δ|. Compute the right-hand sides exactly from B2 and report the slack. Search for counterexamples to the conjectured form, which is plausible since n_K ≍ Reg while their bound is exponential in ρ_∞.

**B5. Limit law of μ/Reg**, a "cubic Lévy constant". Compute its distribution by discriminant window, separately for h = 1 and h > 1 and for cyclic versus non-Galois fields. Test whether it is concentrated, its mean and whether it drifts; the census median is 0.62 and three large fields give about 0.47. Also measure the fraction of simplicial fields as D → ∞.

**B6. Unboundedness, the open problem.**
- *Orders:* compute μ for O_p = Z + pO_K for p ≤ 10³ in a few fields, and check numerically the claim μ ≫ p^{1/3}/(log p)² of Cor. 6.8(c).
- *Maximal orders:* search parametric families (cubic analogues of Yasaki's long-period families, and families built from products of small-norm elements) for Reg/(log D)² → ∞ with visible parametric structure in the Voronoi output. A family with closed-form perfect forms would give an unconditional sup n_K = ∞.

**B7. Degree 4 and 5, totally real.** Run small censuses (LMFDB ranges) and test the m-dimensional sandwich and n ≥ 2mμ/τ_m. Record the full face vector (f_0, …, f_{m−2}) to test the Dehn–Sommerville relations of the torus complex. The aim is a d ≥ 4 upper bound for n_K in terms of μ_K, which the manuscript leaves open.

---

## 4. Goal 2 — classification of small perfect and extreme forms

**C1. The floor.**
- *Observed:* in the census, μ_K = 1 always comes with n_K = 2, the one-hexagon torus (27 fields with D ≤ 20000). μ_K = 2 always gives n_K = 4.
- *Question:* is n_K = 1 impossible over totally real cubic fields? It would require one perfect form with exactly four extreme squares tiling the torus by one quadrilateral. This contrasts with Komatsu–Watanabe, where n_K = 1 for complex cubics.
- *Task:* prove or refute n_K ≥ 2. The computation supplies the list of possible local configurations.

**C2. Fields with μ_K = 1**, the cubic analogue of Leibak–Porter–Ling. Produce the complete list up to D ≤ 10⁶ and identify the families: Ennola, Shanks, cyclic cubics, and others. Many of the 27 have square discriminant (49, 81, 169, 361, …), i.e. are cyclic. For those, test whether the C3 Galois action fixes the hexagon. Then conjecture a classification, and prove the finiteness or infinitude of the non-family members. The proved constraint is Reg ≤ ¼(log 2D)².

**C3. Atlas of combinatorial types.** For each field, record the torus tiling: vertex degrees 3–6, face sizes, and the isomorphism type of the CW structure on T². Then:
- tabulate which tilings occur for small n_K, for instance all fields with n_K ≤ 8;
- find the smallest fields realising each vertex degree 4, 5, 6 and each face size;
- decide whether the maximal degree 6, six extreme minimal vectors up to sign, is attained by the hexagonal (A₂-type) trace form only.

**C4. Extreme forms and γ_K.**
- For every census field, list the perfect classes that are eutactic, hence extreme, and the global maximiser γ_K.
- Tabulate γ_K against Leibak's bound γ_K ≤ γ_m D^{1/m} and Porter–Mendelsohn's bound γ_K ≤ m e^{2ρ_∞}.
- Determine the fields where the bound γ_m D^{1/m} is attained.

**C5. Small minimal vectors.** Record the norms of extreme squares, giving the distribution of |N(x)| of extreme squares relative to √D. Test whether every non-unit extreme square has norm ≤ c·D^{1/3}, which connects to the A-reducibility of Leibak–Porter–Ling. Do not develop unit-reducibility as a topic: record only the A-reducibility constant per field as a by-product.

---

## 5. Thread 2: volume of the fundamental domain

By Koecher, F_K is the union of the Voronoi cones of the perfect classes. The engine therefore gives Vol(F_K¹) **exactly** as a finite sum of simplex volumes after triangulating each cone.

**V1.** Compute Vol(F_K¹) and the minimum and maximum Voronoi-cone volumes for the census. Compare with:
- Porter–Mendelsohn's upper bound 1/n!;
- the lower bound |Δ|^{1/2}/(n! e^{2nρ}γ^n) per cone;
- the heuristic Siegel-type normalisation.

This settles numerically how far apart the two sides of the volume sandwich are, and whether the missing lower bound on Vol(F_K¹) has a chance.

**V2.** Test the proposed relation Vol(X_B⁺/Λ_K) = R_K with ratio to |Δ|^{1/2} ≈ h_K, and the identification of F_K¹ with Yatsyna's polytope K(f), both heuristic in memory. Report as HEURISTIC with data.

---

## 6. Higher rank: Thread 1 and Hermite–Humbert constants

**H1. Rank-2 engine over real quadratic fields first.**
- *Algorithm:* binary Humbert-form Voronoi in Leibak's additive normalisation.
- *Validation:* Leibak's complete list of binary perfect forms over Q(√6), and Ong's lists for Q(√2), Q(√3), Q(√5).
- *Rigor:* Baeza–Coulangeon–Icaza–O'Ryan and Pohst–Wagner work in the **multiplicative** (Icaza) normalisation γ_{K,2}. Their extreme-form lists for Q(√2), Q(√3), Q(√5), Q(√6), Q(√13), Q(√21) are a separate benchmark. The additive and multiplicative theories differ (BCIO exhibit perfect tuples with coefficients outside K), so keep both labelled and never conflate them.

**H2. Initial perfect forms from the unary census.** By Leibak's Thm 1.1, if ax² is perfect over a totally real K and g is a perfect Z-form, with one of the two lattices of E-type, then ag is perfect over K. The unary data of §3 thus seeds rank-n Voronoi for every census field. Use the e-type-lattices skill for which g are E-type.

**H3. Binary and ternary counts over cubic fields.** Run the rank-2 Voronoi for small totally real cubic fields to get P_{K,2}, and P_{K,3} if feasible. Compare with:
- the Thread 1 lower bound log P_{K,n} ≥ (d/2)n² log n − O_K(n²);
- Moshchevitin's lower bound over Z, adapted via Kitaoka;
- the number of perfect forms obtained by tensor constructions.

This measures how much of the count tensor constructions account for.

**H4. Gunnells–Yasaki Thm 5.3.** For every perfect form produced in H1–H3, test sign-genericity of the minimal vectors. Search systematically for a counterexample to the genericity assumption for m ≥ 3, or for evidence that it always holds.

---

## 7. Regression benchmarks (all must pass; `tests/`)

| case | expected |
|---|---|
| cubic, D = 49, 257, 404 | (μ, n) = (1, 2), (2, 4), (4, 7); 404 has degrees 3⁶4¹ |
| cubic, D = 473, 1708, 6856, 6997 | (2, 4), (12, 21), (23, 46), (13, 26) |
| cubic, D = 16689 | (77, 152) |
| Ennola a = 2…12; Shanks a = −1…11 | μ = 1 except Ennola a = 3 (μ = 2) |
| census D ≤ 20000 | 832 fields, 27 with μ = 1, all certificates pass |
| complex cubic v³ + ℓv − 1 (Komatsu–Watanabe) | n = 1, S(a₀) = {±1, ±v} |
| real quadratic | Yasaki's tabulated values; Leibak–Porter–Ling h_d = 1 list |
| binary over Q(√6) | Leibak's list (additive) |
| binary over Q(√13) | Pohst–Wagner γ_{K,2} (multiplicative) |
| unit groups | regulator agrees with LMFDB to 10⁻¹⁰ (e.g. 3.3.49.1: 0.525454682123) |

---

## 8. Rules for Claude Code (`CLAUDE.md`)

1. **Certify before reporting.** Never report a count without passing certificates; never repair a failing record by hand.
2. **Labels.** Every claim in `paper/` or a summary carries PROVED, CERTIFIED or CONJECTURE.
3. **References.**
   - Cite only papers present in the project files or verified online.
   - Theorem numbers must be checked against the source text in the project files.
   - Unverified items are marked `[REF?]`.
   - Never declare a reference non-existent without searching the project files first. This is the mistake made with Komatsu–Watanabe.
4. **Disagreements.** When old and new computations disagree, the exact certified one wins. Record the discrepancy in `data/ERRATA.md`.
5. **Determinism.** Seeds and precision are fixed; every run is reproducible from the record's provenance.
6. **Performance.** Profile before optimising. Parallelise over fields only; the per-field algorithm stays sequential and exact.

---

## 9. Phases and deliverables

| phase | content | exit criterion | output |
|---|---|---|---|
| P0 | fix KW citation, repo skeleton, port `voronoi_unary.gp`, tests | §7 cubic rows pass | — |
| P1 | E1, E3, E4 (degree m, orders, relmin, eutaxy) | m = 2 and m = 3 benchmarks; ℓ ≥ μ on the census | — |
| P2 | B1, B2 (D ≤ 10⁶), C1–C5, B3–B5 | census complete, atlas tables | Paper A: Euler bound, regulator sandwich, census, γ_K tables, Porter–Mendelsohn comparison |
| P3 | E2, B6, B7, V1–V2 | Komatsu–Watanabe reproduced; degree-4/5 data | Paper B: classification of small and extreme unary forms; higher degree |
| P4 | H1–H4 | rank-2 benchmarks pass | Paper C: higher rank via unary seeds; Thread 1 data; Gunnells–Yasaki test |

**Most interesting expected outcomes, in order:**
1. An unconditional proof that sup n_K = ∞ via an explicit family.
2. A proof or refutation of n_K ≥ 2 for totally real cubics.
3. The limiting law of μ/Reg.
4. The first exact tables of the unary Hermite–Humbert constant γ_K for cubic fields.
5. A d ≥ 4 Euler-type upper bound for n_K from face-vector data.
6. Rank-n perfect-form counts over cubic fields seeded by unary forms.
