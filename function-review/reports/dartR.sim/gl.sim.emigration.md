# Review: gl.sim.emigration (dartR.sim)

- Family mode: analysis (moves individuals between populations of a genlight object)
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 5effe05 (`origin/dev`); dartR.base f9f1be8 (`origin/dev`, includes the `rbind.dartR` fix of dartR.base PR #422)
- Datasets: `testset.gl` subsets (EmmacMDBCond, EmmacMDBCudg, EmmacMDBForb; loci with no calls removed); `testset.gl` in full; `possums.gl` (documented example)
- Baseline: `tests/testthat/test-gl.sim.emigration.R` (snapshot captured pre-review, 26 expectations, all pass)

## Verdict

**Standards: Needs work** — inputs are not checked, so most mistakes end in errors raised deep inside R; `ind.metrics$pop` is not updated for individuals that move; the standard skeleton is missing.
**Spec: Needs work** — with a probability matrix (`emi.m`), individuals move in the opposite direction to the one the matrix gives. Emigration is applied one population pair at a time, so individuals can move twice, and any population that empties stops the function.

What works: deterministic emigration (`emi.table`) moves the requested numbers in the documented direction, and genotypes, `ind.metrics` rows and `latlon` stay with their individuals now that `rbind()` keeps individual metadata (dartR.base PR #422).

## Findings

**F1 [HIGH, confidence: high] — `emi.m` moves individuals in the reverse direction (principle: model correctness)**
`R/gl.sim.emigration.r:80, 90` — the probabilistic branch stores migrants as `migs[from, to]`, but the moving loop reads `migs[to, from]`, the convention used by `emi.table`.
Failure scenario: two populations; `emi.m` sends the emigrants of EmmacMDBCond to EmmacMDBCudg, and those of EmmacMDBCudg stay (diagonal). With `perc.mig = 0.5`, no EmmacMDBCond individual moves and 4 EmmacMDBCudg individuals move to EmmacMDBCond. The number of emigrants is also taken from the wrong population. With asymmetric migration (source–sink, stepping stone, unequal sizes) the simulated gene flow is reversed. The dartR GUI offers `emi.m` as the "Probabilistic emigration matrix".
Proposed change: store as `migs[to, from]`, so both modes use the documented from = column, to = row. **Consequence: output of `emi.m` runs changes; individuals move in the direction given by the matrix.**

**F2 [MEDIUM, confidence: high] — individuals move twice and emigration order matters (DOC5)**
`R/gl.sim.emigration.r:87–97` — population pairs are processed in turn, and emigrants are drawn from the population as it stands at that point, immigrants included. The documentation mentions this as a limitation.
Failure scenario: three populations, `emi.table` moving 5 from Cond to Cudg and 5 from Cudg to Forb. Three of Cond's emigrants end in Forb, so Cudg sends only 2 of its own individuals and Forb receives Cond genotypes the table never asked for. With symmetric `emi.m`, `perc.mig = 0.2` on two populations of 10, the mean fraction of individuals away from their original population is 0.183 over 300 runs, because some movers return.
Proposed change: draw every population's emigrants from its residents at the start, then move them all at once, relabelling their population. Each individual moves at most once, and the realised numbers equal the table. Relabelling in place replaces the per-pair `rbind()` calls (0.39 s per run for `possums.gl`'s 10 × 10 example). **Consequence: seeded outputs change; individuals move at most once; the order of individuals in the output changes (still grouped by population).**

**F3 [MEDIUM, confidence: high] — a population that empties stops the function (principle: model correctness)**
`R/gl.sim.emigration.r:92–94` — removing the last individual from a population fails with "Subsetting resulted in zero individuals"; asking for more emigrants than a population holds fails with "cannot take a sample larger than the population".
Failure scenario: `emi.table` moving all 10 individuals of a population stops; `perc.mig = 0.99` stops at random, depending on the draw. On `testset.gl`, whose populations include two of size 1, the documented example pattern (one migrant between every pair) stops.
Proposed change: allow a population to empty, with a warning at `verbose >= 1`: in genlight output it disappears from `popNames()`; in list output its element is `NULL`, and `NULL` elements are accepted as input so repeated calls can reuse the same matrices. Stop with a clear error, naming the population, when `emi.table` asks for more emigrants than it holds.
(Wording corrected in Phase C: the first draft said an emptied population keeps its level. Keeping an empty level makes `seppop()` warn and `nPop()` count it, so it is dropped instead.)

**F4 [MEDIUM, confidence: high] — moved individuals keep their old population in `ind.metrics` (DAT2)**
`R/gl.sim.emigration.r:98–99` — `pop()` is updated for movers, `ind.metrics$pop` is not.
Failure scenario: after moving 5 individuals, 5 rows of `ind.metrics` name a population that differs from `pop()`. `gl.reassign.pop(as.pop = "pop")` or any function reading `ind.metrics$pop` puts the migrants back in their source population.
Proposed change: set `ind.metrics$pop` to the new population when the column exists.

**F5 [MEDIUM, confidence: high] — inputs are not checked (FS4, FS5)**
`R/gl.sim.emigration.r:57–84`
- Neither `emi.table` nor `perc.mig` + `emi.m` given: "'x' must be an array of at least two dimensions".
- An `emi.m` column of zeros (a population with no destinations): "NA in probability vector"; the natural reading is that nobody leaves.
- An `emi.table` of the wrong dimension (3 × 3 for 2 populations) runs silently using its top-left corner; a non-square `emi.m` gives "subscript out of bounds".
- `perc.mig` is not checked to be between 0 and 1.
- A list of length 1 is treated as a single genlight (`seppop` fails); an unnamed list loses its population names (`levels(pop())` is `NULL`).
Proposed change: `utils.check.datatype()` on genlight input (each element for lists); stop unless one mode is fully specified; check square matrices of the number of populations, non-negative whole numbers in `emi.table`, `perc.mig` in [0, 1]; treat an all-zero `emi.m` column as no emigration; accept data frames (the GUI reads the matrices with `read.csv()`); decide list or genlight by class, not length; name unnamed list elements by their population.

**F6 [LOW, confidence: high] — standard skeleton missing (FS2, FS3, FS8, FS9, VRB2)**
`R/gl.sim.emigration.r:50–109` — no `verbose` argument, no start/end messages; the history gains one `rbind.dartR(...)` entry per pair moved but none for `gl.sim.emigration`; `if (!is.null(pop.size))` guards nothing; commented `fbm` code remains.
Proposed change: add `verbose = NULL`, `utils.flag.start()`, completion message and one history entry (genlight output); summary of moves at `verbose >= 3`; remove dead code. Adds one argument (`verbose`) at the end of the signature.

**F7 [LOW, confidence: high] — documentation (DOC1, DOC5, DOC7 proposed)**
`R/gl.sim.emigration.r:1–48`
- `perc.mig` is called a "Percentage", but it is a proportion: `perc.mig = 10` sends everyone. The dartR GUI slider already offers it as a fraction from 0 to 1.
- The description of moving twice and population order must follow change 2 if approved.
- No `@family`; `@author` lacks the Author(s)/Custodian parts.
Proposed change: document `perc.mig` as a proportion in [0, 1], keeping the argument name; describe the model; add the tags.

## Proposed changes

1. `emi.m` moves individuals from column to row, as documented (F1). **Consequence: output of `emi.m` runs changes; migration follows the direction in the matrix.**
2. Emigrants drawn from residents at the start and moved all at once, by relabelling (F2). **Consequence: seeded outputs change; each individual moves at most once; output order changes (still grouped by population).**
3. Populations may empty; clear error when `emi.table` asks for more emigrants than a population holds (F3).
4. `ind.metrics$pop` follows `pop()` (F4).
5. Input checks, including data frames from the GUI, all-zero `emi.m` columns and list handling (F5). **Consequence: wrong-dimension tables and out-of-range `perc.mig`, accepted today, stop with an error.**
6. Standard skeleton and one history entry (F6).
7. Documentation (F7).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run. PLT, DEP: not applicable.
- Spec: direction of `emi.table` and `emi.m`; moving twice (3-population chain); realised migration fraction (300 runs); emptying a population; too many emigrants; wrong dimensions; zero columns; data-frame input; list input (named, unnamed, length 1); genotype, `ind.metrics`, `latlon` and `loc.metrics` integrity; history; timing — run.
- Callers: no internal callers; no calls in sibling `dartR.*` packages; dartr2shiny `shiny_fun/Fun_gl.sim.emigration.R` passes `perc.mig` (slider 0–1), and `emi.m`/`emi.table` as data frames from `read.csv(row.names = 1)` — read. The added `verbose` argument does not affect that call.
- Results depend on the `rbind.dartR` fix (dartR.base PR #422, merged); with an older dartR.base, `ind.metrics` is lost for some inputs (see the #422 downstream tests).
- FBM path (DAT6): not applicable (the `fbm` argument is commented out).
- Google Group / GitHub issues: not searched.

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | consequence approved |
| 2 | approved | Luis | consequence approved |
| 3 | approved | Luis | empty-population handling detailed in Phase C (see F3) |
| 4 | approved | Luis |  |
| 5 | approved | Luis | consequence approved |
| 6 | approved | Luis |  |
| 7 | approved | Luis |  |

## Outcome

| Change | Evidence | Snapshot result |
|---|---|---|
| 1 | test "emi.m moves individuals from column to row" | flipped: with emigrants of EmmacMDBCond sent to EmmacMDBCudg, Cond individuals move and none from Cudg (was the reverse); realised fraction moving 0.203 for `perc.mig = 0.2` (was 0.183) |
| 2 | test "individuals move at most once" | flipped: Cond -> Cudg 5 and Cudg -> Forb 5 give 0 Cond individuals in Forb (was 3); output grouped by population. `possums.gl` example: 0.011 s per run (was 0.39 s) |
| 3 | test "populations can empty" | flipped: emptying EmmacMDBCond returns 20 individuals in EmmacMDBCudg with a warning (was an error); `perc.mig = 0.99` runs; list output has a `NULL` element that is accepted back. Too many emigrants: "asks for more emigrants than individuals in population(s): <pop> (<asked> of <size>)" |
| 4 | test "ind.metrics$pop follows pop()" | flipped: 0 mismatches (was 5) |
| 5 | tests "input errors", "list input", "data.frame matrices" | flipped: wrong dimensions, `perc.mig = 10`, fractional `emi.table`, non-genlight list elements stop with named errors; zero `emi.m` column moves nobody (was an error); `list(x)` returns a list; unnamed list elements named after their population (were `NULL`). Data frames still accepted |
| 6 | test "history and messages" | new: one `gl.sim.emigration` history entry (no `rbind.dartR` entries); summary at `verbose = 3`; silent at `verbose = 0` |
| 7 | `devtools::document()`; `@family` adds `gl.sim.emigration` to the "Other simulation functions" links in the other `.Rd` files | docs only |

Unchanged and passing: `emi.table` direction and counts; genotypes, `ind.metrics` ids and `latlon` stay with their individuals. Stepping-stone run on all of `testset.gl` at `verbose = 3` completes, and the result goes through `gl.report.heterozygosity()`. Tests: `test-gl.sim.emigration.R` 43 expectations pass; full suite 201 pass. `R CMD check` on the tracked files: 0 errors; 1 warning (installed packages built under R 4.4.3, local environment); 1 note (timestamps). NEWS.md updated.
PR: pending.

```json
{
  "function": "gl.sim.emigration",
  "package": "dartR.sim",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "5effe05",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "principle:model-correctness", "status": "approved", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "principle:model-correctness", "status": "approved", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DAT2", "status": "approved", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "FS3", "status": "approved", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 7}
  ],
  "coverage_skipped": ["FBM path: fbm argument disabled", "forum/issues search: not run"],
  "status": "pr-open",
  "pr": null
}
```
