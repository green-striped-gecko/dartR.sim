# Review: gl.report.nall (dartR.sim)

- Family mode: report (allele-count rarefaction curve with per-population points)
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 94a7379 (`origin/dev`); dartR.base f9f1be8 (`origin/dev`); the installed dartR.sim (used by the parallel workers) built from 94a7379
- Datasets: `platypus.gl` (3 populations, structured), `testset.gl` (30 populations, missing data), `possums.gl` (documented example), `testset.gs`
- Baseline: `tests/testthat/test-gl.report.nall.R` (snapshot captured pre-review, 16 expectations, all pass)

## Verdict

**Standards: Needs work** — `verbose = 0` still prints, SilicoDArT input fails inside a worker with an unrelated message, one simulation job crashes, and `reps = 0` runs two replicates.
**Spec: Needs work** — the documentation says the curve comes from subsampling individuals, but it comes from individuals simulated under Hardy–Weinberg from the pooled allele frequencies. The two agree closely on the reference data except at the largest sample sizes, where the simulated curve never reaches 1.

What works: the input object is returned untouched and nothing is added to its history; results are computed and returned whether or not the plot is shown (PLT3); per-population allele counts are correct.

## Findings

**F1 [MEDIUM, confidence: high] — the curve is simulated, not subsampled (DOC5)**
`R/gl.report.nall.r:5–8, 33–40, 156–161` — each replicate calls `gl.sim.ind(x, n)`, which draws n new individuals under Hardy–Weinberg and linkage equilibrium from the pooled allele frequencies. The title, `@description` and `@details` describe "subsampling individuals from the pooled set", and `@param x` says missing data prevent subsampling.
Failure scenario: a reader interprets a population below the ribbon as having fewer alleles than a random subsample of the same size would have. The ribbon is the expectation under a panmictic population with the pooled frequencies, so structure between populations is not part of it. Compared with true subsampling (10 replicates), the curves differ by 0.01 on `platypus.gl` (0.850 vs 0.837 at 5 individuals) and by up to 0.04 on `testset.gl`. At the full sample the simulated curve stays below 1 (0.961 at 250 individuals on `testset.gl`; a subsample of everyone is 1.0 by definition).
Proposed change: documentation only — describe the null model (simulated panmictic individuals from the pooled frequencies), explain that the curve need not reach 1, restate the interpretation guide in those terms, and drop the claim that missing data must be removed (allele frequencies use the called genotypes). Replacing the simulation with true subsampling would change results and is not proposed.

**F2 [MEDIUM, confidence: high] — `verbose = 0` is not silent (VRB3)**
`R/gl.report.nall.r:96, 123` — `gl.filter.allna(x, by.pop = TRUE)` runs without `verbose`, and the default `plot.colors.pop = gl.colors("dis")` prints its own start/end messages.
Failure scenario: `gl.report.nall(x, verbose = 0)` prints 8 lines ("Starting gl.filter.allna", …, "Completed: gl.colors").
Proposed change: pass `verbose = 0` to `gl.filter.allna()` and use `gl.colors("dis", verbose = 0)` as the default. The dartR GUI passes its own colours and is not affected.

**F3 [MEDIUM, confidence: high] — SilicoDArT input fails after the cluster starts (FS4)**
`R/gl.report.nall.r:119` — `utils.check.datatype()` accepts SilicoDArT, then `gl.sim.ind()` refuses it inside a worker.
Failure scenario: `testset.gs` gives "task 1 failed - Fatal Error: inappropriate object passed to function, found SilicoDArT expecting SNP" after the cluster has started.
Proposed change: `utils.check.datatype(x, accept = "SNP")` before any work. **Consequence: SilicoDArT input stops at the start with a clear error (it already fails today, later).**

**F4 [MEDIUM, confidence: high] — one simulation job crashes; `reps = 0` runs twice; arguments unchecked (FS5)**
`R/gl.report.nall.r:173, 178–190`
- `simlevels = 5, reps = 1` (one job): `foreach(.combine = rbind)` returns a vector, and `colnames()` fails with "'names' attribute [2] must be the same length as the vector [1]".
- `reps = 0`: `1:reps` is `c(1, 0)`, so two replicates run.
- `simlevels`, `reps` and `ncores` are not checked (fractional or negative values reach `expand.grid()`/`makeCluster()`).
Proposed change: bind results so a single job keeps two columns; stop unless `simlevels` are whole numbers ≥ 1 and `reps`, `ncores` whole numbers ≥ 1. **Consequence: `reps = 0`, which runs today, stops with an error.**

**F5 [LOW, confidence: high] — a cluster is started even for one core (principle: efficiency)**
`R/gl.report.nall.r:176–186` — `makeCluster(ncores)` runs on every call. On `possums.gl` the call takes 4.3 s elapsed for 0.2 s of computation.
Proposed change: with `ncores = 1`, run the replicates in the current session without a cluster. Results for a given seed then follow the session's random numbers; with `ncores > 1` nothing changes.

**F6 [LOW, confidence: high] — documentation and messages (DOC1, DOC5, DOC7 proposed, FS3)**
`R/gl.report.nall.r:1–116`
- `@param ncores` says "[default 10]"; the default is 2.
- `@return` says `sim` is a `data.frame`; it is a tibble.
- `utils.flag.start()` gets the outdated `build = "v.2023.3"`.
- No `@family`; `@author` lacks the Author(s)/Custodian parts; `@param x` lacks the family wording; the description uses a curly apostrophe (DOC6).
Proposed change: correct the defaults and the return type (return `sim` as a `data.frame`, matching the documentation); drop `build`; fix the tags. **Consequence: `sim` changes class from tibble to `data.frame` (same columns and values).**

## Proposed changes

1. Document the simulated null model and its interpretation (F1). Documentation only.
2. `verbose = 0` silent (F2).
3. SNP data checked at the start (F3). **Consequence: SilicoDArT input stops at the start with a clear error.**
4. Single job works; `simlevels`, `reps`, `ncores` validated (F4). **Consequence: `reps = 0` and invalid values stop with an error.**
5. No cluster for `ncores = 1` (F5).
6. Documentation fixes; `sim` returned as a `data.frame` (F6). **Consequence: `sim` is a `data.frame`, not a tibble.**

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run. Report family: input untouched and no history appended (checked); PLT3 checked (results returned with `plot.display = FALSE`).
- Spec: simulated curve against true subsampling on `platypus.gl` and `testset.gl`; full-sample value; per-population points; documented example; SilicoDArT; missing data; `reps = 0`, one job, `simlevels` above `nInd`; verbosity; timing — run.
- Callers: dartr2shiny `shiny_fun/Fun_gl.report.nall.R` passes `x`, `simlevels`, `reps`, `ncores`, `plot.theme` and an evaluated colour vector; none of the proposed changes alter that call. No callers in sibling `dartR.*` packages.
- The parallel workers load the installed dartR.sim, so tests use the installed `gl.sim.ind()`; installed from 94a7379 for this review.
- FBM path (DAT6): not tested.
- Google Group / GitHub issues: not searched.

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis |  |
| 2 | approved | Luis |  |
| 3 | approved | Luis | consequence approved |
| 4 | approved | Luis | consequence approved |
| 5 | approved | Luis |  |
| 6 | approved | Luis | consequence approved |

## Addendum notes (not changed)

- The `dev` manifest again had a duplicate `gl.sim.offspring` row and #55–#57 marked `pr-open` after their merges; rewritten on this branch (bookkeeping only).

## Outcome

| Change | Evidence | Snapshot result |
|---|---|---|
| 1 | `man/gl.report.nall.Rd` regenerated; test "curve is simulated and stays below 1 at the full sample" | docs only; curve values unchanged |
| 2 | test "verbose = 0 is silent" | flipped: no output (was 8 lines) |
| 3 | test "SilicoDArT stops at the start" | flipped: "SilicoDArT" error before any cluster starts (was "task 1 failed" from a worker) |
| 4 | test "single job works; arguments checked" | flipped: `simlevels = 5, reps = 1` returns one row (was a `names` error); `reps = 0`, fractional `simlevels`, `ncores = 0` stop with named errors (`reps = 0` ran two replicates) |
| 5 | test "parallel and single-core runs agree" | new: `ncores = 1` runs in the session; means agree with `ncores = 2` within 0.03, points identical. `possums.gl`, 3 sample sizes × 10 replicates: 0.55 s at 1 core, 3.7 s at 2 cores (cluster start-up) |
| 6 | test "sim is a data.frame"; `devtools::document()` | flipped: `sim` is a `data.frame` with columns `Npop, mnall, low, high` (was a tibble); `ncores` default documented as 2; `build` dropped |

Unchanged and passing: returned structure, input untouched, results with `plot.display = FALSE`, per-population points; documented example at `verbose = 3`. Tests: `test-gl.report.nall.R` 22 expectations pass; full suite 314 pass. `R CMD check` on the tracked files: 0 errors; 1 warning (installed packages built under R 4.4.3, local environment); 1 note (timestamps). NEWS.md updated.
PR: #59.

```json
{
  "function": "gl.report.nall",
  "package": "dartR.sim",
  "family": "report",
  "skill_version": "2.0.0",
  "commit": "94a7379",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "VRB3", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "FS4", "status": "approved", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "approved", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "principle:efficiency", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved", "change": 6}
  ],
  "coverage_skipped": ["FBM path: not tested", "forum/issues search: not run"],
  "status": "pr-open",
  "pr": 59
}
```
