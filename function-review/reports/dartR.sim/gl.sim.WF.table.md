# Review: gl.sim.WF.table (dartR.sim)

- Family mode: analysis (builds the input table for `gl.sim.WF.run`; no genlight output)
- Date: 2026-09-23
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: fd727f3 (`dev_luis`, level with `origin/dev` 55c9126 in content)
- Datasets: `inst/extdata/ref_variables.csv` (dartR.sim), `testset.gl` (first 200 loci, with synthetic `@chromosome = "1"` and sorted `@position`, because `testset.gl` has no chromosome slot), synthetic recombination maps written to `tempfile()`
- Baseline: `tests/testthat/test-gl.sim.WF.table.R` (snapshot captured pre-review, 20 expectations, all pass)

## Verdict

**Standards: Needs work** — the FS skeleton is present, but fatal errors are raised with an empty message, there is no parameter validation, and the quote-and-eval block is copied three times.
**Spec: Rework** — three documented input routes (`...` overrides, `real_freq = TRUE` alone, a recombination map whose intervals differ from `chunk_bp`) crash or return a wrong table, and `loci_deleterious` above `chunk_number` returns a different number of loci than requested.

What works: the default CSV route is reproducible under `seed`, and the neutral, advantageous and mutation counts match the request.

## Findings

**F1 [HIGH, confidence: high] — `...` override corrupts `chromosome_name` (DOC5)**
`R/gl.sim.WF.table.r:162` — when any argument is passed through `...`, every string variable is wrapped in single quotes again. The distribution variables are cleaned by the `gsub()` at lines 186–191, but `chromosome_name` is not, so the CSV value `"1"` becomes the literal string `"\"1\""`.
Failure scenario: `gl.sim.WF.table(file_var = fv, x = x, interactive_vars = FALSE, real_loc = TRUE)` stops with "Chromosome name is not in the genlight object" although chromosome 1 exists; `gl.sim.WF.table(..., chunk_number = 20)` returns `chr_name = "\"1\""` for every locus. The same breaks the recombination-map and targets-of-selection lookups. Workaround today: also pass `chromosome_name = "1"`.
Proposed change: strip double quotes from `chromosome_name` alongside the distribution variables, and wrap only the values that were overridden.

**F2 [HIGH, confidence: high] — `real_freq = TRUE` with `real_loc = FALSE` is broken (DOC5)**
`R/gl.sim.WF.table.r:510, 560, 616, 667` — `s`, `h`, `q` and `type` for the real loci are set only when `real_loc == TRUE`. With `real_freq = TRUE` alone, those rows stay `NA`.
Failure scenario: with the default CSV (gamma `s`, equation `q`), the cap step at line 691 indexes with `NA` and the call stops with "missing values are not allowed in subscripted assignments of data frames". If both distributions are `"equal"`, the call succeeds but returns 200 rows with `type = NA`, which `gl.sim.WF.run` (line 187) does not recognise as neutral or real, so the real frequencies are never used. The Shiny tooltip documents this route as supported.
Proposed change: gate the real-locus assignments on `real_loc == TRUE | real_freq == TRUE`.

**F3 [HIGH, confidence: high] — recombination map is assumed to have `chunk_bp`-sized intervals (DOC5)**
`R/gl.sim.WF.table.r:208, 416, 429` — map midpoints are generated as `seq(chunk_bp / 2, chr_length, chunk_bp)` and per-locus `c` divides by `chunk_bp`; the map's own `from`/`to` columns are ignored except for the last `to`. Line 208 also sets the whole row (including `from`/`to`) to 0 when `cM` is `NA`.
Failure scenario: (a) 50 kb map with `chunk_bp = 1e5`: stops with "'vec' must be sorted non-decreasingly". (b) 200 kb map, 0 cM over the first 10 Mb and 4 cM per interval after: the returned table puts 1.96 Morgans inside the first 10 Mb, where the map has none — silently wrong recombination. (c) `NA` in the last `cM`: `chr_length` becomes 0 and the call stops with "wrong sign in 'by' argument". The documented example maps (`fly_recom_map.csv`) are not shipped (F10), so users build their own and hit these cases.
Proposed change: take midpoints and interval lengths from the map's `from`/`to`; replace `NA` in the `cM` column only; place neutral loci over the map's length.

**F4 [HIGH, confidence: high] — `loci_deleterious` is rounded per chunk (DOC5)**
`R/gl.sim.WF.table.r:321` — when `loci_deleterious >= chunk_number`, each chunk gets `round(loci_deleterious / chunk_number)` targets.
Failure scenario: with 100 chunks, requesting 149 gives 100 deleterious loci, 150 gives 200, 250 gives 200. The mismatch is silent; the value written back into the returned `ref_vars` is still the requested number.
Proposed change: give each chunk `floor(n / chunk_number)` targets and add one to a random sample of `n %% chunk_number` chunks, so the total equals the request.

**F5 [MEDIUM, confidence: high] — value caps are silent, undocumented and keyed to the deleterious settings (VRB4 proposed rule, DOC5)**
`R/gl.sim.WF.table.r:690–699` — `q > 0.5` is set to 0.5 when `q_distribution_del != "equal"`, and `s > 1` to 0.99 / `s < -0.5` to -0.5 when `s_distribution_del != "equal"`, but both caps apply to every locus class.
Failure scenario: `q_distribution_adv = "equal", q_adv = 0.8` returns advantageous `q = 0.5`, because the deleterious method is `"equation"`. `s_distribution_del = "equal"` with exponential advantageous `s` (`exp_rate = 1`) returns advantageous `s` down to -6.58 (123 of 200 loci beyond -0.5), while the same advantageous settings with gamma deleterious `s` are capped at -0.5. None of the caps is documented or reported.
Proposed change: cap each class according to its own distribution setting, document the caps in `@details`, and report the number of capped loci at `verbose >= 1`.

**F6 [MEDIUM, confidence: high] — fatal errors carry no message (FS5, VRB2)**
`R/gl.sim.WF.table.r:201, 224, 232, 238, 407` — the pattern `message(error("...")); stop()` prints the text as a message and then raises an error with an empty message.
Failure scenario: inside `tryCatch()`, `suppressMessages()`, a Shiny app or dartr2shiny, the user or caller sees only "Error:" with no reason.
Proposed change: `stop(error("..."))`.

**F7 [MEDIUM, confidence: high] — unknown `...` names are silently ignored; no input validation (API2 proposed rule, FS5)**
`R/gl.sim.WF.table.r:153–159` — names in `...` that are not variables in the CSV match nothing and are dropped. There is no check that `file_var` exists when `interactive_vars = FALSE`, or that `x` is a genlight.
Failure scenario: a typo `chunk_numbr = 20` runs with the default 100 chunks and no warning. Omitting `file_var` with `interactive_vars = FALSE` fails with R's "argument "file_var" is missing" from inside `read.csv`.
Proposed change: stop with the list of unknown names; check `file_var` and `x` at the top of the function.

**F8 [LOW, confidence: medium] — `sample()` on a single candidate position (principle: R `sample()` length-1 trap)**
`R/gl.sim.WF.table.r:336, 382` — when `seq(start, end, by = sample_resolution)` yields one value (an interval in a targets file shorter than `sample_resolution`), `sample(v, size = 1)` draws from `1:v` instead of returning `v`. When `sample_resolution` rounds to 0 (many targets in short intervals), `seq()` stops with "invalid '(to - from)/by'".
Failure scenario: a targets file mixing a 1 Mb interval with a 5 kb interval places the 5 kb interval's target anywhere between 1 and its start position. Not reproduced here because no targets file ships with the package.
Proposed change: index-sample (`v[sample.int(length(v), size)]`) and floor `sample_resolution` at 1.

**F9 [LOW, confidence: high] — `seed` resets the session RNG (principle: no hidden side effects)**
`R/gl.sim.WF.table.r:93–95` — `set.seed(seed)` changes the global random stream for everything the user runs afterwards.
Failure scenario: two calls to `gl.sim.WF.run` after `gl.sim.WF.table(seed = 1)` in a replicate loop start from a stream the user did not set. This is common practice across dartRverse, so the proposal is documentation only.
Proposed change: state in `@param seed` that it sets the session seed.

**F10 [LOW, confidence: high] — documentation does not match the package (DOC1, DOC2, DOC5, DOC6 proposed, DOC7 proposed)**
`R/gl.sim.WF.table.r:1–81`
- `fly_recom_map.csv` and `fly_targets_of_selection.csv` (line 58–59) are not shipped by dartR.sim, dartR.data or dartR.base; the required columns (`Chr`, `from`, `to`, `cM`; `chr_name`, `start`, `end`, `targets`) are not documented anywhere.
- Line 46 points to `ref_variables.csv` in dartR.data; the example (line 71) uses the dartR.sim copy.
- Line 21 URL (`georges.biomatix.org/dartR`) is dead.
- Line 48 uses curly quotes (DOC6); `@param verbose` text differs from the DOC2 standard; `@author` has no `Author(s):` part (DOC7).
- `@param x` says "Name of the genlight object"; it takes the object. Line 52 says `x$position and x$chromosome`; the code requires `@chromosome` to contain `chromosome_name`, which `testset.gl` does not have.
Failure scenario: a user following `@details` cannot find the example files or learn the file format, and builds a map that hits F3.
Proposed change: ship the two example files in `inst/extdata` (or drop the references), document the column formats, and correct the items above.

**F11 [INFO, confidence: high] — standards housekeeping (FS3, STY1)**
`R/gl.sim.WF.table.r:102–104, 112–125, 162–175, 131–136, 145–149, 178–182` — `utils.flag.start(build = "Jody")` uses the outdated `build=` argument; the quote-wrapping and assign-and-eval blocks are repeated three times.
Failure scenario: the F1 defect exists because a fix in one copy was not mirrored in another.
Proposed change: drop `build=`; move quoting and evaluation into one internal helper.

## Proposed changes

1. Strip quotes from `chromosome_name` after evaluation and wrap only overridden values (F1).
2. Assign `s`, `h`, `q` and `type` to real loci when `real_loc` or `real_freq` is `TRUE` (F2).
3. Use the map's own `from`/`to` for midpoints and interval length, replace `NA` only in `cM`, and span neutral loci over the map length (F3). **Consequence: numerical output (`c`, `loc_cM`, `loc_bp`) changes for any recombination map whose interval size differs from `chunk_bp`.**
4. Distribute `loci_deleterious` so the total equals the request (F4). **Consequence: the number and positions of deleterious loci change whenever `loci_deleterious >= chunk_number` and is not a multiple of it; seeded results change.**
5. Cap `q` and `s` per locus class using that class's own distribution setting; document the caps; report the capped count at `verbose >= 1` (F5). **Consequence: returned `q`/`s` values change when deleterious and advantageous settings differ (e.g. advantageous `q_adv = 0.8` is no longer cut to 0.5 under an equation deleterious `q`; exponential advantageous `s` is capped at -0.5 regardless of the deleterious setting).**
6. Replace `message(error(...)); stop()` with `stop(error(...))` (F6).
7. Stop on unknown `...` names; validate `file_var` and `x` up front (F7). **Consequence: calls with misspelt or obsolete `...` arguments that run today will error.**
8. Index-safe sampling and `sample_resolution >= 1` (F8).
9. Document that `seed` sets the session seed (F9).
10. Documentation fixes listed under F10, including shipping or removing the example map/targets files; run `devtools::document()` (F10).
11. Drop `build=` and factor the quote/eval code into one helper (F11).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run. DAT and PLT: not applicable (no genlight or plot returned).
- Spec: behaviour vs roxygen and Shiny tooltips on the CSV route, `...` route, `real_loc`, `real_freq`, recombination maps (50 kb, 100 kb, 200 kb, `NA` cM), deleterious/advantageous/mutation counts, caps — run with `devtools::load_all()`.
- Interactive route (`interactive_vars = TRUE`, Shiny): SKIPPED — needs a browser session; F1's quoting logic in that branch was read, not run.
- Targets-of-selection file route: SKIPPED as a run — no example file ships (F10); F8 is from reading the code.
- Numerical check of `q_equilibrium` against Crow & Kimura: SKIPPED — it lives in `utils.sims.r`, out of scope for this function.
- Downstream use in `gl.sim.WF.run`: read for F2 only.
- FBM path (DAT6): not applicable — only `@chromosome`, `@position` and `nLoc()` are read from `x`.
- Google Group / GitHub issues: not searched in this pass.

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis |  |
| 2 | approved | Luis |  |
| 3 | approved | Luis |  |
| 4 | approved | Luis |  |
| 5 | approved | Luis |  |
| 6 | approved | Luis |  |
| 7 | approved | Luis | as proposed (strict: error, not warning) |
| 8 | approved | Luis |  |
| 9 | approved | Luis |  |
| 10 | approved | Luis |  |
| 11 | approved | Luis |  |

## Addendum findings (found while applying; not applied, need approval)

**F12 [MEDIUM, confidence: high] — map rate assigned by midpoint, gaps inherit neighbouring rates (DOC5)**
`R/gl.sim.WF.table.r` (RECOMBINATION MAP block) — each locus takes the rate of the interval whose midpoint is at or below it, so a locus in the first half of an interval uses the previous interval's rate, and loci in gaps between map intervals use a neighbour's rate. Kept as is in change 3 to stay within the approved scope (it reproduces the old default output exactly).
Failure scenario: `fly_recom_map.csv`, chromosome 2L: the table accumulates 0.579 Morgans against 0.553 Morgans in the map (+4.7%).
Proposed change: assign each locus to the interval whose `from`–`to` contains it; zero rate outside intervals. **Consequence: `c` and `loc_cM` change for every run with a map file; default runs (uniform `chunk_cM`) do not change.**

**F13 [out of scope: gl.sim.WF.run] — real frequencies with all-`NA` loci**
With `real_freq = TRUE` and a genlight containing loci with no calls (e.g. `testset.gl[, 1:200]`), `gl.sim.WF.run` stops with "NA in probability vector". The table from this function is valid; after `gl.filter.callrate(threshold = 1)` both real routes run end to end. Recommend covering this when `gl.sim.WF.run` is reviewed.

## Outcome

| Change | Evidence | Snapshot result |
|---|---|---|
| 1 | test "`...` override keeps chromosome_name clean" | flipped: `chr_name` `"\"1\""` → `"1"`; `real_loc` via `...` returns 200 real loci |
| 2 | test "real_freq = TRUE with real_loc = FALSE"; `gl.sim.WF.run` end to end on complete-call data | flipped: error → 200 real loci, no `NA` |
| 3 | test "map intervals of any size"; before/after S5 (aligned map) identical, S6/S8 (fly map) changed | flipped: 0 Morgans in the first 10 Mb of the 200 kb map (was 1.96); 50 kb map with `NA` runs. Fly 2L now spans 23.1 Mb (was 9.95 Mb) |
| 4 | test "loci_deleterious total equals the request" | flipped: 149/150/250 → 149/150/250 (was 100/200/200); multiples of `chunk_number` unchanged |
| 5 | test "caps follow each class's own setting"; before/after S2, S9 identical | flipped: advantageous `s` ≥ -0.5 with equal deleterious `s`; `q_adv = 0.8` kept. Also: neutral `q_neutral > 0.5` is no longer cut to 0.5 under a non-equal deleterious `q` (same rule: neutral loci have no distribution setting) |
| 6 | test "errors carry their message" | new |
| 7 | test "unknown ... argument stops" | flipped: silent → error naming `chunk_numbr` |
| 8 | before/after S7 (fly targets) identical — `sample.int` draws the same stream | no diff |
| 9, 10 | `devtools::document()`; `man/gl.sim.WF.table.Rd` regenerated; fly files in `inst/extdata`, test "fly example files ship and run" | new |
| 11 | `utils.wf.ref.values()` in `R/utils.sims.r`; `build=` dropped | returned `ref_vars` values lose the extra quote wrapping (`'"gamma"'` → `"gamma"`); `gl.sim.WF.run` reads only `q_neutral`, `real_freq`, `real_loc` from it |

Before/after comparison on nine seeded scenarios (default; mixed selected/mutation loci; `real_loc`; 99 deleterious; aligned 100 kb map; fly map; fly targets; fly map + targets; equation `h`/`q` with equal advantageous `s`): identical `reference` (by `all.equal`) except the two fly-map scenarios, which fall under change 3.
Tests: `tests/testthat/test-gl.sim.WF.table.R` 28 expectations pass. Callers checked: `dartR.captive` (`utils.classes.diagnostics.relatedness.r`, passes `file_var`, `x`, `interactive_vars`) and `dartr2shiny` template (`file_var`, `interactive_vars`) — no `...` use, unaffected. NEWS.md created with the entry.
PR: #47.

```json
{
  "function": "gl.sim.WF.table",
  "package": "dartR.sim",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "fd727f3",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 2},
    {"id": "F3", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 3},
    {"id": "F4", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "VRB4", "status": "approved", "change": 5},
    {"id": "F6", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "approved", "change": 6},
    {"id": "F7", "severity": "MEDIUM", "confidence": "high", "rule": "API2", "status": "approved", "change": 7},
    {"id": "F8", "severity": "LOW", "confidence": "medium", "rule": "principle:sample-length-1", "status": "approved", "change": 8},
    {"id": "F9", "severity": "LOW", "confidence": "high", "rule": "principle:no-hidden-side-effects", "status": "approved", "change": 9},
    {"id": "F10", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved", "change": 10},
    {"id": "F11", "severity": "INFO", "confidence": "high", "rule": "FS3", "status": "approved", "change": 11}
  ],
  "coverage_skipped": ["interactive Shiny route: needs browser", "targets file route: no example file", "q_equilibrium numerics: in utils.sims.r", "forum/issues search: not run"],
  "status": "pr-open",
  "pr": 47
}
```
