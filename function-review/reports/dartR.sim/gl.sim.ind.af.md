# Review: gl.sim.ind.af (dartR.sim)

- Family mode: analysis (simulates genotypes from per-population allele frequencies)
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 13b8a3f (`origin/dev`)
- Datasets: `platypus.gl` (documented example); synthetic frequency tables (two populations listed Z before A with frequencies 0.9 and 0.1; one population at 0.3 for 500 loci)
- Baseline: `tests/testthat/test-gl.sim.ind.af.R` (snapshot captured pre-review, 27 expectations, all pass)

## Verdict

**Standards: Needs work** — the output is a plain `genlight` with a malformed locus-metrics table, inputs are only partly checked, and the standard skeleton (verbosity, messages, history) is missing.
**Spec: Needs work** — frequencies and sizes are matched to populations by position, not by name. When populations are not listed in alphabetical order, or `popn` is a factor, some populations are simulated with another population's frequencies or size, and nothing warns about it.

What works: Hardy–Weinberg proportions (0.490 / 0.420 / 0.090 observed for p = 0.3, n = 2000) and the documented example, whose `gl.allele.freq()` table is already sorted alphabetically.

## Findings

**F1 [HIGH, confidence: high] — frequencies go to the wrong population when populations are not in alphabetical order (principle: model correctness)**
`R/gl.sim.ind.af.r:82, 122, 193` — population names and sizes follow `unique(df$popn)` (order of appearance), but `split(df, df$popn)` returns populations in sorted order, and `df_pops[[i]]` is taken by position.
Failure scenario: a table listing population Z (frequency 0.9) before A (frequency 0.1), `pop.sizes = c(Z = 100, A = 10)`. The individuals labelled Z have a mean frequency of 0.100 and those labelled A 0.902, so the two populations have swapped frequencies. Any user-built table not sorted alphabetically is affected; `gl.allele.freq()` output is sorted, so the documented workflow is not.
Proposed change: convert `popn` to character on entry and take each population's frequencies by name (`df_pops[[pop_name]]`). **Consequence: output changes for tables not in alphabetical order (populations get their own frequencies).**

**F2 [HIGH, confidence: high] — factor `popn` with named `pop.sizes` gives populations the wrong sizes (principle: model correctness)**
`R/gl.sim.ind.af.r:108` — `pop.sizes[pops_in_df]` indexes by a factor, which R treats as the integer level codes, not the names.
Failure scenario: `popn` a factor with levels `c("Z", "A")`, `pop.sizes = c(A = 10, Z = 100)`: Z gets 10 individuals and A gets 100. `gl.allele.freq()` returns `popn` as a factor, so named sizes on its output are exposed whenever the names are not in level order. A related effect: unused factor levels (a table subset from a larger one) stop with "All populations must provide the same set of loci. Populations failing this check: Q", which points to the loci rather than the unused level.
Proposed change: same as F1 — `popn` as character on entry removes both the level-code indexing and the unused levels. **Consequence: output changes for factor `popn` with named sizes (populations get the sizes asked for).**

**F3 [MEDIUM, confidence: high] — "frequency of the first allele" describes the wrong allele (DOC5)**
`R/gl.sim.ind.af.r:11–15, 27, 40` — the genotype counts copies of the allele whose frequency is supplied. In a `genlight`, 0/1/2 counts the alternate (second) allele, and `gl.allele.freq()`'s `frequency` is the alternate-allele frequency. The code is consistent with that; the documentation says "first allele".
Failure scenario: a user supplies reference-allele frequencies as the docs ask. The simulated data have the frequencies mirrored (0.8 becomes 0.2 in `gl.alf()`); heterozygosity and FST are unaffected, but any frequency-based comparison with the source data is inverted.
Proposed change: documentation only — say "frequency of the alternate allele (the allele counted in the genotype, as returned by `gl.allele.freq()`)".

**F4 [MEDIUM, confidence: high] — returned object lacks standard metadata; locus metrics table is malformed (DAT2, DAT5)**
`R/gl.sim.ind.af.r:228–250` — the output is a plain `genlight`, not a `dartR` object. `utils.reset.flags()` builds `loc.metrics` from scratch, so its first column is named `array(NA, nLoc(x))` and every metric is `NA`. `ind.metrics` has `fid`, `iid`, `sex`, `phenotype` but no `id` or `pop`, which dartR functions look up.
Failure scenario: downstream functions that join on `ind.metrics$id` or read `loc.metrics$AlleleID` find nothing; `gl.compliance.check()` must be run first to fill in the missing parts.
Proposed change: return a `dartR` object; build `loc.metrics` with `AlleleID` (locus names) before resetting flags, so no malformed column appears; add `id` and `pop` to `ind.metrics`, keeping `fid`, `iid`, `sex` and `phenotype` for back-compatibility.

**F5 [MEDIUM, confidence: high] — `pop.sizes` and duplicate rows are not checked (FS5)**
`R/gl.sim.ind.af.r:86–90, 124–136, 190`
- `pop.sizes = c(2.7, 3.9)` is truncated to 2 and 3 without a message.
- `pop.sizes = c(NA, 3)` stops with "missing value where TRUE/FALSE needed".
- A population with the same locus twice passes the locus check (`setequal`), and the second row is silently dropped.
Proposed change: stop unless `pop.sizes` are whole numbers ≥ 1 with no `NA`; stop when a population has duplicate locus rows, naming them. **Consequence: fractional sizes and tables with duplicate rows, accepted today, now stop with an error.**

**F6 [LOW, confidence: high] — standard skeleton missing (FS2, FS3, FS8, FS9, VRB2)**
`R/gl.sim.ind.af.r:62–257` — no `verbose` argument, no start/end messages, no history entry; errors use plain `stop()` instead of `stop(error(...))`. `exists("utils.reset.flags")` is always true because the function is imported, so the guard does nothing. Commented-out `fbm` and `check.loci` code remain.
Proposed change: add `verbose = NULL` with `gl.check.verbosity()`, `utils.flag.start()`, completion message and history; use `error()`; drop the dead guard and commented code. Adds one argument (`verbose`) at the end of the signature; no callers in `dartR.*` siblings or dartr2shiny.

**F7 [LOW, confidence: medium] — C++ compiled on every new session (principle: runtime dependency)**
`R/gl.sim.ind.af.r:142–180` — `make_chr` and `make_geno` are compiled by `Rcpp::cppFunction()` at run time. The first call in a session takes 2.9 s here, most of it compiling, and the function fails on machines without a C++ toolchain (Windows without Rtools). Both helpers do what one vectorised line of R does: draw a haplotype matrix as `runif() < q` and add the two haplotypes.
Failure scenario: a Windows user without Rtools cannot run the function at all.
Proposed change: replace both helpers with vectorised R (drawing 40 million alleles takes about 0.8 s, similar to the current run). **Consequence: results for a given seed change (the random stream is consumed in a different order); the distribution is the same.**

**F8 [LOW, confidence: high] — documentation gaps (DOC1, DOC5, DOC7 proposed)**
`R/gl.sim.ind.af.r:1–60`
- `@return` says IDs are `"0_<popIndex>_<i>"`; they use the population name (`"0_SEVERN_ABOVE_1"`).
- `@details` says sex is assigned "in alternating blocks"; it alternates per individual (m, f, m, …).
- No `@family`, no `@param verbose`, and `@author` lacks the Author(s)/Custodian parts.
Proposed change: correct the text; add the tags; regenerate `man/`.

## Proposed changes

1. Match frequencies and sizes to populations by name, with `popn` as character (F1, F2). **Consequence: output changes for tables not in alphabetical order and for factor `popn` with named sizes; populations get their own frequencies and sizes.**
2. Describe `frequency` as the alternate-allele frequency (F3). Documentation only.
3. Return a `dartR` object with a proper `loc.metrics` (`AlleleID`) and `ind.metrics` `id`/`pop` added; existing columns kept (F4).
4. Validate `pop.sizes` (whole numbers ≥ 1, no `NA`) and reject duplicate locus rows within a population (F5). **Consequence: fractional sizes and duplicate rows, accepted today, now stop with an error.**
5. Standard skeleton: `verbose`, messages, history, `error()`, clean-up (F6).
6. Replace the run-time C++ with vectorised R (F7). **Consequence: results for a given seed change; the distribution does not.**
7. Documentation fixes (F8).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run. PLT: not applicable. DEP: `Rcpp` and `data.table` are in Imports, so no DEP1 guard is needed (DEP3's Suggests policy is a proposed rule; not raised).
- Spec: documented example on `platypus.gl`; Hardy–Weinberg proportions; population order (character and factor `popn`, named and unnamed sizes); unused factor levels; fractional/`NA` sizes; duplicate rows; metadata; first-call timing — run.
- Callers: no internal callers in `R/`; no calls in sibling `dartR.*` packages; dartr2shiny holds a source copy (`input_generator/dartR.sim/gl.sim.ind.af.r`) with no argument specification — read.
- FBM path (DAT6): not applicable (the `fbm` argument is commented out).
- Toolchain failure (F7) on Windows: not reproduced; inferred from `Rcpp::cppFunction()` requiring a compiler.
- Google Group / GitHub issues: not searched.

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | consequence approved |
| 2 | approved | Luis |  |
| 3 | approved | Luis |  |
| 4 | approved | Luis | consequence approved |
| 5 | approved | Luis |  |
| 6 | approved | Luis | consequence approved |
| 7 | approved | Luis |  |

## Addendum notes (not changed)

- `utils.reset.flags()` (dartR.base) creates `loc.metrics.flags` from scratch with a first column named `array(NA, 1)` when the object has none. Every function that builds a new object and resets flags gets this column; it is a dartR.base issue, not changed here.

## Outcome

| Change | Evidence | Snapshot result |
|---|---|---|
| 1 | tests "populations get their own frequencies", "factor popn with named pop.sizes matches by name", "unused factor levels are ignored" | flipped: Z listed before A now has mean frequency > 0.8 (was 0.100); factor `popn` with `c(A = 10, Z = 100)` gives Z = 100, A = 10 (was swapped); unused level `Q` is ignored (was a locus error) |
| 2 | `man/gl.sim.ind.af.Rd` regenerated | docs only |
| 3 | test "metadata of the returned object" | flipped: `dartR` object; `loc.metrics$AlleleID` equals `locNames()`, no `array(NA, nLoc(x))` column; `ind.metrics` columns `id, pop, fid, iid, sex, phenotype` |
| 4 | tests "pop.sizes checks", "duplicated locus rows stop" | flipped: `c(2.7, 3.9)`, `c(NA, 3)`, `c(0, 3)` stop with "whole numbers"; a duplicate row stops naming `Z-l1` |
| 5 | test "verbose messages"; history entry | new: `verbose = 3` prints the summary, `verbose = 0` is silent; one history entry |
| 6 | `make_chr`/`make_geno` removed; `importFrom(Rcpp, cppFunction)` dropped from NAMESPACE (other files call `Rcpp::` directly) | seeded genotypes change; Hardy-Weinberg test (0.49/0.42/0.09) and example frequency test pass. First call on the example: 0.03 s (was 2.9 s, mostly compiling). 1000 individuals x 20,000 loci: 2.2 s (was 1.6 s once compiled) |
| 7 | `devtools::document()`; `@family` adds `gl.sim.ind.af` to the "Other simulation functions" links in six other `.Rd` files | docs only |

Unchanged and passing: documented example structure (150 x 383, IDs `0_<pop>_<i>`), Hardy-Weinberg proportions, sex alternation, `loc.all` placeholder. The output goes through `gl.filter.callrate()` and `gl.report.heterozygosity()` without errors. Tests: `test-gl.sim.ind.af.R` 41 expectations pass; full suite 178 pass. `R CMD check` on the tracked files: 0 errors; 1 warning (installed packages built under R 4.4.3, local environment); 1 note (timestamps). NEWS.md updated.
PR: #54.

```json
{
  "function": "gl.sim.ind.af",
  "package": "dartR.sim",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "13b8a3f",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "principle:model-correctness", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "principle:model-correctness", "status": "approved", "change": 1},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 2},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DAT2", "status": "approved", "change": 3},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "approved", "change": 4},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "FS3", "status": "approved", "change": 5},
    {"id": "F7", "severity": "LOW", "confidence": "medium", "rule": "principle:runtime-dependency", "status": "approved", "change": 6},
    {"id": "F8", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved", "change": 7}
  ],
  "coverage_skipped": ["FBM path: fbm argument disabled", "Windows toolchain failure: inferred, not reproduced", "forum/issues search: not run"],
  "status": "pr-open",
  "pr": 54
}
```
