# Review: gl.sim.apply (dartR.sim)

- Family mode: analysis (utility: applies a function over `gl.sim.WF.run()` output)
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0. The function was written by the same model on 2026-09-23; this review is not independent of its author.
- Package commit: eb39994 (`origin/dev`); dartR.base f9f1be8 (`origin/dev`)
- Datasets: `gl.sim.WF.run()` default example, 2 iterations (generations 1 and 10); `testset.gl`
- Baseline: `tests/testthat/test-gl.sim.apply.R` (24 existing expectations plus 5 added for this review, all pass)

## Verdict

**Standards: Ready** — the skeleton, input checks, verbosity and documentation follow the conventions; the input is not modified.
**Spec: Needs work** — three result shapes that `fun` can reasonably return end in an R error after every call to `fun` has run, and a repeated generation silently replaces an earlier result.

What works: tagging by `sim.vars$generation` rather than list names; binding of vectors and consistent data frames; nested lists for other results; checks before `fun` runs; errors in `fun` reported with iteration and generation.

## Findings

**F1 [MEDIUM, confidence: high] — `fun` returning `NULL` stops the run (DOC5)**
`R/gl.sim.apply.r:140–141` — `attr(res, "iteration") <- it` fails on `NULL`.
Failure scenario: a function called for its side effect (saving a file, printing a plot) that returns `NULL` or `invisible(NULL)` stops with "attempt to set an attribute on NULL" after the first generation.
Proposed change: accept `NULL` results: they add no rows to a bound data frame and stay `NULL` in the nested list.

**F2 [MEDIUM, confidence: high] — results with different columns fail at the end (DOC5)**
`R/gl.sim.apply.r:158–177` — results are bound with base `rbind()`, which needs identical columns. The binding runs after every call to `fun`.
Failure scenario: `fun` returns data frames whose columns differ between generations (a population column that appears only after a split, a statistic that is absent when a population is extinct), a mix of data frames and vectors, or data frames with row names in some generations only. The run stops with "names do not match previous names" or "numbers of columns of arguments do not match", and all results are lost.
Proposed change: bind with `data.table::rbindlist(use.names = TRUE, fill = TRUE)` (already imported), so missing columns are filled with `NA`; return a `data.frame`. Document it.

**F3 [LOW, confidence: high] — a repeated generation overwrites the earlier result; repeated `iteration` values run twice (FS5)**
`R/gl.sim.apply.r:96–102, 142` — results are stored by the name `generation_<n>`.
Failure scenario: an iteration whose list holds generation 10 twice (a list assembled by hand, or two runs combined) returns one result for it, silently. `iteration = c(1, 1)` runs `fun` twice on iteration 1 and keeps one.
Proposed change: stop before `fun` runs when an iteration has repeated generations, naming them; use `unique(iteration)`.

## Proposed changes

1. `NULL` results allowed (F1).
2. Bind with missing columns filled with `NA` (F2). **Consequence: results that stop with an error today are returned, with `NA` where a column is absent.**
3. Repeated generations stop with an error; repeated `iteration` values are used once (F3). **Consequence: input with a repeated generation, which runs today, stops with an error.**

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run. DAT, PLT: not applicable (the input is not modified; no plot).
- Spec: `NULL` results; data frames with differing columns; data frames mixed with vectors; row names in some generations only; repeated generations and iterations; factor vectors; `verbose = 3`; genlight objects without `sim.vars` — run. Whether `gl.sim.WF.run()` stores a generation twice: `every_gen` of 3, 4 and 7 with 10 generations gives no repeats.
- Callers: none in `R/` or sibling `dartR.*` packages; dartr2shiny does not use it.
- Google Group / GitHub issues: not applicable (added 2026-09-23).

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis |  |
| 2 | approved | Luis | consequence approved |
| 3 | approved | Luis | consequence approved |

## Outcome

| Change | Evidence | Snapshot result |
|---|---|---|
| 1 | test "fun returning NULL is allowed" | flipped: nested list with `NULL` elements (was "attempt to set an attribute on NULL"); `NULL` mixed with vectors adds no rows |
| 2 | test "results with different columns are bound with NA" | flipped: columns `iteration, generation, a, b` with `NA` where absent (was "names do not match previous names"); data frames mixed with vectors bound (was "numbers of columns of arguments do not match") |
| 3 | test "repeated generations stop; repeated iterations run once" | flipped: "iteration 1, generation 10 appears more than once" before `fun` runs (was one result kept); `iteration = c(1, 1)` runs once |

Unchanged and passing: the 24 existing expectations (vector and data-frame binding, nested lists, flat lists and single genlights, tagging from `sim.vars`, checks, located errors). Tests: `test-gl.sim.apply.R` 35 expectations pass; full suite 358 pass. `R CMD check` on the tracked files: 0 errors; 1 warning (installed packages built under R 4.4.3, local environment); 1 note (timestamps). NEWS.md updated.
PR: pending.

```json
{
  "function": "gl.sim.apply",
  "package": "dartR.sim",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "eb39994",
  "verdict_standards": "ready",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 2},
    {"id": "F3", "severity": "LOW", "confidence": "high", "rule": "FS5", "status": "approved", "change": 3}
  ],
  "coverage_skipped": ["forum/issues search: not applicable, function added 2026-09-23"],
  "status": "pr-open",
  "pr": null
}
```
