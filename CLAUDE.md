# perfect_forms — notes for Claude sessions

- Push only to `pav10/perfect_forms` (working branch `claude/sharp-clarke-v0s20b`,
  default branch `main`).  The repositories `pav10/sails` and `pav10/hermite-humbert`
  are **read-only** for this project: read their code and data, never commit, push
  or open PRs there.
- `PLAN.md` is the research plan; `README.md` documents commands, record schema and
  status.  Keep both in sync with the code.
- Focus: perfect unary forms over **maximal orders** of totally real fields (the
  number n_K of perfect forms first).  The `basis=` hook of `engine.field.Order` for
  non-maximal orders is not a current focus: do not build on it or report data from it
  unless asked.
- Tests before every push: `python tests/test_engine.py` (fast, ~2 min);
  `python tests/test_engine.py --slow` also re-derives the 832-field cubic census
  against `data/legacy/census_cubic_D20000_gp.jsonl`.
- Runs expected to take much more than an hour need the user's go-ahead.

## Rules (PLAN.md section 8)

1. **Certify before reporting.**  Never report a count from a record whose status is
   not `CERTIFIED` (Euler characteristic 0, unit index 1, independent mu count, edge
   count).  Never repair a failing record by hand: it stays in `*.quarantine`.
2. **Labels.**  Every claim in `paper/`, `docs/` or a summary carries one of PROVED,
   CERTIFIED (certified computation) or CONJECTURE (conjecture/heuristic), never mixed.
3. **References.**  Cite only papers present in the project files or verified online;
   theorem numbers must be checked against the source text; unverified items are
   marked `[REF?]`.  Never declare a reference non-existent without searching the
   project files first (the Komatsu--Watanabe mistake).
4. **Disagreements.**  When old and new computations disagree the exact certified one
   wins; record the discrepancy in `data/ERRATA.md`.
5. **Determinism.**  PARI's seed is fixed (`setrand(1)`); every record stores its
   provenance (git hash, PARI/cypari2/Python versions, precision, wall time).
6. **Performance.**  Profile before optimising.  Parallelise over fields only; the
   per-field algorithm stays sequential and exact.

## Data conventions

- `fields/`: field lists (JSONL `{"label","pol","D"}`, xz-compressed), built with
  `python -m engine fields` and cross-checked (cubic: nflist == LMFDB field by field).
- `data/census/`: certified records (`*.jsonl.xz`); `analysis/` reads them and never
  edits them.  `data/legacy/`: results of the earlier GP code, kept for regression.
- The manuscript `paper/perfect-unary-cubic.tex` has no build here (no LaTeX in the
  sandbox); if LaTeX becomes available build with `pdflatex` twice and commit the PDF.
