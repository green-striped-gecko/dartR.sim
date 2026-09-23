# Review: gl.diagnostics.sim (dartR.sim)

- Family mode: report (plots simulated He and FST against theoretical expectations)
- Date: 2026-09-23
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: e4210df (`origin/dev`, includes the fixes to `gl.sim.WF.table`, `gl.sim.WF.run` and `gl.sim.ind`)
- Datasets: `gl.sim.WF.run()` outputs (the example: 2 populations of 10; a theory run: 2 populations of 50, `number_transfers = 1` every generation, 500 neutral loci, 200 generations × 12 replicates); an exact identity-by-descent recursion computed independently
- Baseline: `tests/testthat/test-gl.diagnostics.sim.R` (snapshot captured pre-review, 7 expectations, all pass)

## Verdict

**Standards: Needs work** — there is no input validation, `verbose = 0` still prints messages, and the function returns only the plot, so the numbers behind it cannot be tested or reused.
**Spec: Rework** — both theoretical curves are wrong for the model `gl.sim.WF.run` simulates. The expected FST is about half the correct value (0.059 against 0.111), and the expected He curve starts from the wrong point whenever the first stored generation is not 0. The function exists to validate the simulator, and it currently flags a correct simulator as wrong.

What works: the observed curves. Mean Nei FST from the simulations (0.114 at equilibrium) matches the exact value for the simulated model (0.111), which also confirms the migration fix in #48.

## Findings

**F1 [HIGH, confidence: high] — the expected FST is about half the correct value (principle: model correctness)**
`R/gl.diagnostics.sim.r:186–188` — the line is `1 / (4 Ne m (n/(n-1))^2 + 1)` with n = 2, i.e. `1 / (16 Ne m + 1)`, and m = `number_transfers / transfer_each_gen / N`. The vertical "equilibrium" line is a heuristic (twice a half-life).
Failure scenario: 2 populations of 50, one individual exchanged each way every generation (m = 0.02, Ne = N = 50, verified for this simulator in #48):

| Source | FST at equilibrium |
|---|---|
| Simulation (Nei FST, generations 101–200, 12 replicates) | 0.114 |
| Exact identity-by-descent recursion for this model | 0.111 |
| Plotted expectation, `1/(16Nm+1)` | 0.059 |
| `1/(1 + 4Nm·n/(n−1))` (m = immigrant fraction from other demes) | 0.111 |

The recursion is: drift within demes (F0 ← 1/2N + (1 − 1/2N) F0), then symmetric exchange ((1−m)² + m² of gene pairs stay together). It also reproduces the approach to equilibrium (0.091 at generation 20 against 0.099 observed). The coded formula uses Takahata's parametrisation, in which m counts migrants drawn from a pool that includes the own deme; the simulator's m counts only immigrants from other demes.
Proposed change: plot the expected FST trajectory from the exact recursion (a curve over generations, starting from the observed FST at the first stored generation), instead of the constant line and the heuristic vertical line. For more than 2 populations with `all_connected`, use the n-deme version (an immigrant comes from each other deme with probability m/(n−1)). For `line` and `circle`, stop with an informative error, since the island recursion does not apply. **Consequence: the plotted expectation changes; the observed curve is unchanged.**

**F2 [HIGH, confidence: high] — the expected He curve starts from the wrong point (principle: model correctness)**
`R/gl.diagnostics.sim.r:110–112` — the expectation is `He_first × (1 − 1/2Ne)^t` with t = the absolute generation. He_first is the observed He at the first *stored* generation, but it is decayed as if it were generation 0.
Failure scenario: phase 1 (10 generations) then phase 2, first stored generation 11: observed He 0.275, expected 0.156. The observed curve then looks far above theory. In the default case (first stored generation 1), the curve is one generation off.
Proposed change: use t − t_first as the exponent. **Consequence: the plotted He expectation changes.**

**F3 [MEDIUM, confidence: high] — simulation variables are parsed with the wrong quote and index (principle: failure path)**
`R/gl.diagnostics.sim.r:172–184`
- `population_size_phase2` is stripped of `'` only. Sizes quoted in the CSV (`"10 10"`, as in `sim_variables.csv`) become `NA`, and the expected FST line disappears with a ggplot warning.
- `number_transfers_phase2` and `transfer_each_gen_phase2` are stored once per dispersal pair but indexed by population numbers (`[pops_fst]`).
Failure scenario: `sim_variables.csv` with `population_size_phase2 = "10 10"`: no expected FST line, and no error.
Proposed change: strip both quote types; read the transfer values for the pair `pops_fst` (or the single value when all pairs share it).

**F4 [MEDIUM, confidence: high] — inputs are not validated (FS5)**
`R/gl.diagnostics.sim.r:73–87` — no check that `x` is `gl.sim.WF.run` output, that `iteration` exists and is non-empty, that there are at least 2 populations, or that `pop_he`/`pops_fst` exist.
Failure scenario: a single-population run returns a plot with empty FST panels (NA warnings only); `iteration = 2` of a 1-iteration run gives "subscript out of bounds"; `pops_fst = c(1, 3)` with 2 populations gives an empty curve.
Proposed change: check these up front and stop with messages that name the problem.

**F5 [MEDIUM, confidence: high] — only the plot is returned (PLT3)**
`R/gl.diagnostics.sim.r:226` — `invisible(p3)` is the only output. Observed and expected values cannot be tabulated, tested or combined across iterations (e.g. with `gl.sim.apply()`).
Proposed change: return `invisible(list(plot = p3, he = <data frame gen/observed/expected by Ne>, fst = <data frame gen/observed/expected>))`. **Consequence: the return type changes from a patchwork to a list; code that uses the return value as a plot must use `$plot`.**

**F6 [LOW, confidence: high] — verbosity leaks and outdated flag (VRB1, FS3)**
`R/gl.diagnostics.sim.r:83, 116` — `dartR.base::gl.colors("dis")` runs at its own default verbosity, so `verbose = 0` still prints "Starting gl.colors…"; `utils.flag.start(build = "Jody")`.
Proposed change: pass `verbose = 0` to `gl.colors`; drop `build =`.

**F7 [LOW, confidence: high] — documentation (DOC1, DOC5, DOC7 proposed)**
`R/gl.diagnostics.sim.r:1–57`
- `@details` gives the FST formula of F1.
- It does not say that five He curves are drawn for Ne to 2Ne, or that the He expectation assumes an isolated deme (migration also slows He loss).
- No `@family`; `@author` lacks the `Author(s):` part; the Crow & Kimura reference repeats its title; `@param Ne` does not say it can be one value per population.
Proposed change: rewrite `@details` for the new expectations and fix the tags.

**F8 [HIGH, confidence: high, out of scope: gl.sim.WF.run] — population labels are scrambled with 10 or more populations**
`R/utils.sims.r` `store()` and `R/gl.sim.WF.run.r` — `pop()` is set from character labels, so the factor levels sort as "1", "10", "11", "12", "2", …; then `popNames<-` renames those levels by position to "1", "2", "3", ….
Failure scenario: 12 populations, no dispersal: every individual born in population 10 is labelled "2", those from 11 are "3", those from 2 are "5", and so on (only population 1 is right). Every per-population analysis of a run with ≥ 10 populations uses the wrong labels.
Proposed change (separate PR, `gl.sim.WF.run`): build the factor with `levels = as.character(pops_vector)` in `store()`.

## Proposed changes

1. Expected FST as a trajectory from the exact identity recursion (2 populations or `all_connected`; error for `line`/`circle`) (F1). **Consequence: the plotted FST expectation changes (about twice the old line for 2 populations).**
2. Expected He decays from the first stored generation (F2). **Consequence: the plotted He expectation changes.**
3. Parse population sizes and transfer values correctly (F3).
4. Validate inputs (F4).
5. Return `list(plot, he, fst)` (F5). **Consequence: the return value is a list, not a patchwork.**
6. Quiet `gl.colors`; drop `build =` (F6).
7. Documentation (F7).
8. Out of scope: fix population labels for ≥ 10 populations in `gl.sim.WF.run` (F8), as a separate PR.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run. DAT: read-only; `x` is not modified.
- Spec: the example; sizes quoted in the CSV; a phase-1 run; one population; iteration and `pops_fst` out of range; missing `Ne`; FST theory against the simulation (12 replicates × 200 generations, 6 min) and against an independent exact recursion — run.
- He expectation with migration: not simulated; the limitation is documented under F7.
- `line`/`circle` dispersal: not simulated; F1 proposes an error for them.
- `plot.file` saving: read, not run.
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
| 8 | approved as separate PR | Luis | next PR, gl.sim.WF.run |
| A1 | approved | Luis | addendum: reshape2 and stringr removed from Imports (unused after change 1/5) |

## Outcome

| Change | Evidence | Snapshot result |
|---|---|---|
| 1 | test "expected FST follows the island recursion"; theory run (6 replicates, 2 populations of 50, m = 0.02) | expected 0.0998 / 0.1162 / 0.1200 at generations 21 / 41 / 101 against observed 0.0958 / 0.1125 / 0.1147; 3 populations `all_connected`: expected 0.082, observed 0.079 at generation 100. `line` stops with an error |
| 2 | test "expected He starts at the first stored generation" | flipped: expected = observed at generation 11 (was 0.156 against 0.275) |
| 3 | test "population sizes quoted in the CSV are read" | flipped: expected FST not `NA` |
| 4 | test "invalid inputs stop with clear errors" | flipped: one population, a bad iteration and a bad `pops_fst` stop with named errors |
| 5 | test "returns plot and tables" | flipped: `list(plot, he, fst)` |
| 6 | test "verbose = 0 prints nothing" | `gl.colors()` call removed (its colours were never used by the plot); `build =` dropped |
| 7 | `devtools::document()`; `man/gl.diagnostics.sim.Rd` regenerated | new `@details` |

A1: `R CMD check` dependencies NOTE resolved.
Tests: `test-gl.diagnostics.sim.R` 14 expectations pass; full suite passes. Callers: none in other `dartR.*` packages or dartr2shiny. NEWS.md updated.
PR: #51.

```json
{
  "function": "gl.diagnostics.sim",
  "package": "dartR.sim",
  "family": "report",
  "skill_version": "2.0.0",
  "commit": "e4210df",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "principle:model-correctness", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "principle:model-correctness", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "principle:failure-path", "status": "approved", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "approved", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "PLT3", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "VRB1", "status": "approved", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved", "change": 7},
    {"id": "F8", "severity": "HIGH", "confidence": "high", "rule": "out-of-scope:gl.sim.WF.run", "status": "approved", "change": 8}
  ],
  "coverage_skipped": ["He with migration: not simulated", "line/circle dispersal: not simulated", "plot.file saving: read only", "forum/issues search: not run"],
  "status": "pr-open",
  "pr": 51
}
```
