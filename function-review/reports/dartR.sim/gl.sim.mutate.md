# Review: gl.sim.mutate (dartR.sim)

- Family mode: analysis (applies random mutations to SNP genotypes)
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: cc29b1f (`origin/dev`); dartR.base f9f1be8 (`origin/dev`)
- Datasets: `testset.gl` (loci with no calls removed: 250 individuals, 252 loci, 12.5% missing); `testset.gs`; `bandicoot.gl` (documented example); `glSim(500, 20000)` for timing
- Baseline: `tests/testthat/test-gl.sim.mutate.R` (snapshot captured pre-review, 16 expectations, all pass)

## Verdict

**Standards: Needs work** — the whole object is converted to a matrix once per mutation, SilicoDArT data are not refused, locus metrics are left stale, and the standard skeleton is missing.
**Spec: Needs work** — when the number of mutations drawn is 0, the loop still runs twice, so every call makes 1–2 mutations whatever `mut.rate` is.

What works: the mutation model itself. A homozygote becomes a heterozygote, and a heterozygote becomes either homozygote. No genotype jumps from 0 to 2, missing genotypes are left alone, and at high rates the number of mutations matches `mut.rate` (10.3 observed against 11.0 expected at 1e-4).

## Findings

**F1 [HIGH, confidence: high] — mutations occur when none are drawn (principle: model correctness)**
`R/gl.sim.mutate.r:32` — `for (ii in 1:nm)` with `nm = 0` iterates over `c(1, 0)`, so two mutation events are applied.
Failure scenario: `mut.rate = 0` changes 1–2 genotypes in every call (10 of 10 calls on `testset.gl`). At the default `1e-6`, `testset.gl` should get 0.13 mutations per call and gets 1.7, 13 times too many. A simulation that calls the function once per generation adds about 2 mutations per generation regardless of the rate, which inflates diversity and new-allele counts in small datasets.
Proposed change: `for (ii in seq_len(nm))`. **Consequence: output changes whenever no mutation is drawn — no genotype changes, as the rate implies.**

**F2 [MEDIUM, confidence: high] — SilicoDArT data get impossible genotypes (FS4, DAT1)**
`R/gl.sim.mutate.r:40–44` — presence/absence scores (0/1) are treated as SNP genotypes, so a mutated 1 becomes 0 or 2.
Failure scenario: `testset.gs` at `mut.rate = 1e-3` returns 24 scores of 2 in a SilicoDArT object.
Proposed change: `utils.check.datatype(x, accept = "SNP")`. **Consequence: SilicoDArT input, which runs today, stops with an error.**

**F3 [MEDIUM, confidence: high] — locus metrics are stale after mutation (DAT4)**
`R/gl.sim.mutate.r:26–52` — genotypes change but `loc.metrics` and `loc.metrics.flags` do not; on `testset.gl` the `maf` flag stays `TRUE`.
Failure scenario: a function that trusts the flag (for example a MAF filter) uses allele frequencies from before the mutations.
Proposed change: reset the flags with `utils.reset.flags()` when at least one genotype changed.

**F4 [MEDIUM, confidence: high] — the whole object is converted to a matrix for every mutation (DAT6 proposed, STY2)**
`R/gl.sim.mutate.r:35` — `as.matrix(x)[ri, rl]` converts every individual and locus to read one genotype.
Failure scenario: 500 individuals × 20,000 loci at `mut.rate = 2e-6` (40 mutations) takes 7.9 s, and the time grows with the number of mutations times the size of the object.
Proposed change: read the genotype from the mutated individual only (`as.matrix(x[ri, ])`, which the code already builds). The random draws happen in the same order, so seeded results do not change.

**F5 [LOW, confidence: high] — `mut.rate` is not checked (FS5)**
`R/gl.sim.mutate.r:31` — `mut.rate = 2` gives "NA/NaN argument" (with an `rbinom` warning); `mut.rate = "a"` gives "invalid arguments".
Proposed change: stop unless `mut.rate` is a single number between 0 and 1.

**F6 [LOW, confidence: high] — standard skeleton missing (FS2, FS3, FS8, FS9)**
`R/gl.sim.mutate.r:26–52` — no `verbose` argument, no start/end messages, no history entry; commented `fbm` code remains.
Proposed change: add `verbose = NULL`, `utils.flag.start()`, completion message, history, a summary of the number of mutations applied at `verbose >= 3`; remove the dead code. Adds one argument (`verbose`) at the end of the signature; the dartR GUI call (`x`, `mut.rate`) is not affected.

**F7 [LOW, confidence: high] — documentation (DOC1, DOC5, DOC7 proposed)**
`R/gl.sim.mutate.r:1–24` — the model is not described (which genotypes a mutation produces; that events at missing genotypes are dropped, so the realised rate is lower by the missing fraction). `@param x` says "Name of the genlight object"; there is no `@family`, no `@details`, and `@author` lacks the Author(s)/Custodian parts.
Proposed change: describe the model in `@details`; fix the tags.

## Proposed changes

1. No mutations when none are drawn (F1). **Consequence: output changes whenever the draw is 0; at low rates most calls now change nothing.**
2. SNP data only (F2). **Consequence: SilicoDArT input stops with an error.**
3. Reset locus-metric flags after mutation (F3).
4. Read only the mutated individual, not the whole matrix (F4). Seeded results unchanged.
5. Validate `mut.rate` (F5).
6. Standard skeleton (F6).
7. Documentation (F7).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run. PLT, DEP: not applicable.
- Spec: `mut.rate = 0`; default rate; rate calibration at 1e-3/1e-4; genotype transitions; missing data; SilicoDArT; invalid `mut.rate`; metadata and flags; documented example; timing on a 500 × 20,000 object — run.
- Callers: no internal callers; none in sibling `dartR.*` packages; dartr2shiny `shiny_fun/Fun_gl.sim.mutate.R` calls `gl.sim.mutate(x, mut.rate)` from a numeric input — read.
- FBM path (DAT6): not tested; the `fbm` argument is commented out. `x@gen[[ri]] <-` writes SNPbin slots and would not work on an FBM-backed object; recorded, not proposed.
- Google Group / GitHub issues: not searched.

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | consequence approved |
| 2 | approved | Luis | consequence approved |
| 3 | approved | Luis |  |
| 4 | approved | Luis |  |
| 5 | approved | Luis |  |
| 6 | approved | Luis |  |
| 7 | approved | Luis |  |

## Addendum notes (not changed)

- The `dev` manifest had `gl.diagnostics.sim` twice and #53/#54 still `pr-open` after the #54 merge conflict was resolved; rewritten on this branch (bookkeeping only).

## Outcome

| Change | Evidence | Snapshot result |
|---|---|---|
| 1 | test "no mutations when none are drawn" | flipped: `mut.rate = 0` changes 0 genotypes in 10 of 10 calls (was 1-2 in every call); default rate averages < 0.5 per call (was 1.7) |
| 2 | test "SilicoDArT input stops" | flipped: error (was scores of 2) |
| 3 | test "flags reset and history added" | flipped: `maf` and frequency flags `FALSE` after mutation; existing `loc.metrics` columns unchanged (`utils.reset.flags()` adds missing metric columns as `NA`) |
| 4 | test "seeded results unchanged"; old and new code on 20 seeds | unchanged: identical genotype matrices for 20 of 20 seeds; seed 11 gives 122 changed genotypes in both. 500 x 20,000 at `mut.rate = 2e-6`: 0.044 s (was 7.9 s) |
| 5 | test "invalid mut.rate stops with a clear error" | flipped: 2, -1 and "a" stop with "between 0 and 1" (were "NA/NaN argument" / "invalid arguments") |
| 6 | same test | new: one history entry; summary at `verbose = 3`; silent at `verbose = 0` |
| 7 | `devtools::document()`; `@family` adds `gl.sim.mutate` to the "Other simulation functions" links | docs only |

Unchanged and passing: transitions (no 0 <-> 2 jumps), missing genotypes untouched, rate calibration at 1e-3, documented example on `bandicoot.gl`. Tests: `test-gl.sim.mutate.R` 22 expectations pass; full suite 221 pass. `R CMD check` on the tracked files: 0 errors; 1 warning (installed packages built under R 4.4.3, local environment); 1 note (timestamps). NEWS.md updated.
PR: pending.

```json
{
  "function": "gl.sim.mutate",
  "package": "dartR.sim",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "cc29b1f",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "principle:model-correctness", "status": "approved", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "FS4", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DAT4", "status": "approved", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "STY2", "status": "approved", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "FS5", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "FS3", "status": "approved", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved", "change": 7}
  ],
  "coverage_skipped": ["FBM path: fbm argument disabled; SNPbin writes incompatible, not proposed", "forum/issues search: not run"],
  "status": "pr-open",
  "pr": null
}
```
