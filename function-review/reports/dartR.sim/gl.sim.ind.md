# Review: gl.sim.ind (dartR.sim)

- Family mode: analysis (simulates individuals from the allele frequencies of a genlight)
- Date: 2026-09-23
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 8e81e6d (`origin/dev`)
- Datasets: `testset.gl` (full, and with all-`NA` loci removed by `gl.filter.allna`); `testset.gs` (all-`NA` loci removed); `testset.gl` loci replicated 40 times (10,080 loci) for timing
- Baseline: `tests/testthat/test-gl.sim.ind.R` (snapshot captured pre-review, 13 expectations, all pass)

## Verdict

**Standards: Needs work** — the function skips the standard skeleton: it has no `verbose`, no start or end flag, no data-type check and no input validation, and it drops the locus metadata.
**Spec: Needs work** — any locus with no calls stops the run, which is why the function's own example fails `R CMD check` on `testset.gl`. SilicoDArT input is silently turned into diploid SNP data.

What works: genotypes are drawn under Hardy-Weinberg equilibrium from the input's alternative-allele frequencies. In 2000 simulated individuals, allele frequencies are within 0.013 of the input and heterozygosity equals 2pq (0.038 and 0.038). Ploidy, locus names, alleles and positions are kept.

## Findings

**F1 [HIGH, confidence: high] — loci without calls stop the function (DOC5)**
`R/gl.sim.ind.r:57–61` — `colMeans(..., na.rm = TRUE)` gives `NaN` for a locus with no calls, and `sample()` stops with "NA in probability vector". `@details` documents the crash and tells users to remove such loci first.
Failure scenario: `gl.sim.ind(testset.gl, n = 10)`, the first line of the function's own example, stops (`testset.gl` has 3 such loci among 255). `R CMD check` on dartR.sim therefore fails at this example. Proposed change: give those loci `NA` genotypes, keeping the loci aligned with `x`, and report the count at `verbose >= 2`.
Correction (Phase C): an earlier version of this finding said `gl.report.nall()` stops on the same data and would need an `na.rm` fix. That is wrong: `gl.report.nall()` removes loci that are all missing in any population before calling `gl.sim.ind()`, and it runs on `testset.gl` both before and after this change.

**F2 [HIGH, confidence: high] — SilicoDArT input is silently turned into diploid SNP data (DAT1, FS4)**
`R/gl.sim.ind.r:56–75` — there is no data-type check. Presence/absence scores (0/1) are halved as if they were allele dosages and written into a diploid object.
Failure scenario: `gl.sim.ind(gl.filter.allna(testset.gs), n = 5)` returns a ploidy-2 genlight with genotypes 0/1/2 that mean nothing biologically.
Proposed change: `utils.check.datatype(x, accept = "SNP")` at the top. **Consequence: SilicoDArT input stops with an error instead of returning meaningless genotypes.**

**F3 [MEDIUM, confidence: high] — locus metadata and chromosome are dropped (DAT2, DAT5)**
`R/gl.sim.ind.r:65–76` — the new object is built without `chromosome`, `@other$loc.metrics` or `loc.metrics.flags`.
Failure scenario: a genlight with chromosomes gives a simulated object with `@chromosome = NULL`, so per-chromosome analyses (LD, `gl.keep.loc` by chromosome) fail. `loc.metrics` has 0 rows (`gl.compliance.check()` rebuilds a minimal one).
Proposed change: carry `chromosome`, copy `loc.metrics` from `x` and reset the metrics flags (`utils.reset.flags()`), so recalculable metrics are not taken as describing the simulated individuals.

**F4 [MEDIUM, confidence: high] — one `sample()` call per genotype makes the function slow (STY2)**
`R/gl.sim.ind.r:58–61` — `apply()` over an n × L matrix calls `sample()` once per cell.
Failure scenario: 1000 individuals × 10,080 loci take 46.3 s; the equivalent draw `rbinom(n * L, 2, p)` takes 0.17 s. Two independent allele draws per locus give the same Hardy-Weinberg genotype distribution. `gl.report.nall()` calls this function `reps × levels` times.
Proposed change: draw genotypes with `rbinom()`. **Consequence: seeded outputs change (same distribution, different random draws).**

**F5 [MEDIUM, confidence: high] — standard skeleton missing (FS2, FS3, FS8, FS9, VRB1)**
`R/gl.sim.ind.r:48–81` — no `verbose` argument, no `gl.check.verbosity()`, no `utils.flag.start()`, no completion message, and no history entry on the returned genlight.
Failure scenario: `gl.set.verbosity()` has no effect, and the object carries no record of how it was made.
Proposed change: add `verbose = NULL` as the last argument, plus the FS2/FS3/FS9 lines and a history entry. Adding an argument at the end does not break existing calls.

**F6 [LOW, confidence: high] — `n` is not validated (FS5)**
`R/gl.sim.ind.r:58` — `n = 0` stops with "subscript out of bounds"; a non-integer `n` is silently truncated by `matrix()`.
Proposed change: stop unless `n` is a single whole number ≥ 1.

**F7 [LOW, confidence: high] — documentation (DOC1, DOC5, DOC6 and DOC7 proposed, style)**
`R/gl.sim.ind.r:1–47`
- `@details` promises a crash that F1 removes, and uses promotional wording ("The beauty of the function is, that it is lightning fast").
- `@param x` says "Name of".
- No `@family`; `@author` lacks the `Author(s):` / `Custodian:` structure (DOC7).
- A commented-out `fbm` parameter is left in the roxygen and the signature.
- The model (Hardy-Weinberg, linkage equilibrium, alternative-allele frequencies) is not stated.
Proposed change: rewrite `@details`, fix the tags, remove the `fbm` leftovers.

## Proposed changes

1. Loci without calls get `NA` genotypes; count reported at `verbose >= 2` (F1). **Consequence: calls that stop today return an object; the example and `R CMD check` pass.**
2. Accept SNP data only (F2). **Consequence: SilicoDArT input errors instead of returning meaningless diploid genotypes.**
3. Carry `chromosome` and `loc.metrics` (flags reset) (F3).
4. Draw genotypes with `rbinom()` (F4). **Consequence: seeded outputs change; the distribution is unchanged.**
5. Add `verbose = NULL`, FS skeleton and history (F5).
6. Validate `n` (F6).
7. Documentation (F7).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run. PLT and DEP: not applicable.
- Spec: structure, allele frequencies and Hardy-Weinberg on 2000 simulated individuals; all-`NA` loci; SilicoDArT; non-genlight input; `n = 0` and `n = 1`; `rbind` of two outputs (no duplicate names); timing at 200 and 10,080 loci — run with `devtools::load_all()`.
- Callers: `gl.report.nall()` (dartR.sim) — run on `testset.gl` after the change (it filters loci with no calls first). No other `dartR.*` package or dartr2shiny calls `gl.sim.ind()`.
- FBM path (DAT6): SKIPPED — the `fbm` argument is commented out and `as.matrix(x)` densifies; not evaluated on an FBM-backed object.
- Google Group / GitHub issues: not searched.

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis |  |
| 2 | approved | Luis |  |
| 3 | approved | Luis |  |
| 4 | approved | Luis |  |
| 5 | approved | Luis |  |
| 6 | approved | Luis |  |
| 7 | approved | Luis |  |

## Outcome

| Change | Evidence | Snapshot result |
|---|---|---|
| 1 | test "loci without calls get NA genotypes"; the example runs on `testset.gl` | flipped: 3 loci with no calls return `NA`, all other genotypes non-missing |
| 2 | test "SilicoDArT input stops" | flipped: error from `utils.check.datatype()` |
| 3 | test "chromosome and loc.metrics are carried" | flipped: chromosome kept; `loc.metrics` has `nLoc` rows; recalculable flags `FALSE` |
| 4 | structure and Hardy-Weinberg test unchanged and passing | 1000 x 10,080 loci: 1.19 s (was 46.3 s) |
| 5 | test "history and verbosity" | new |
| 6 | test "n is validated" | flipped: `n = 0` and `n = 2.5` give "whole number" errors |
| 7 | `devtools::document()`; `man/gl.sim.ind.Rd` regenerated | drift example: 0, 32, 48 fixed loci at generations 1, 25, 50 |

Tests: `test-gl.sim.ind.R` 19 expectations pass; full suite passes. NEWS.md updated.
PR: #50.

```json
{
  "function": "gl.sim.ind",
  "package": "dartR.sim",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "8e81e6d",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DAT1", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DAT2", "status": "approved", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "STY2", "status": "approved", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "FS2", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "FS5", "status": "approved", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved", "change": 7}
  ],
  "coverage_skipped": ["FBM path: fbm argument disabled", "forum/issues search: not run"],
  "status": "pr-open",
  "pr": 50
}
```
