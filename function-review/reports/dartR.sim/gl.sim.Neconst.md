# Review: gl.sim.Neconst (dartR.sim)

- Family mode: analysis (simulates SNP genotypes of a population at mutation–drift equilibrium)
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: cc29b1f (`origin/dev`); dartR.base f9f1be8 (`origin/dev`)
- Datasets: simulated output only (the function takes no genlight input): 20 individuals × 20,000 loci for the spectrum; 20 × 5,000 for `mutation_rate`; 100 × 100,000 for timing
- Baseline: `tests/testthat/test-gl.sim.Neconst.R` (snapshot captured pre-review, 22 expectations, all pass)

## Verdict

**Standards: Needs work** — inputs are not checked, `verbose` is accepted but never used, the history records `gl.compliance.check()` instead of this call, and dead code remains.
**Spec: Needs work** — about a fifth of the requested SNPs come back monomorphic, and the spectrum of the rest is distorted at the tails. The documentation presents `mutation_rate` as a driver of the simulation, but at realistic values it has no effect on the output.

What works: the model choice. Wright's Beta(4Nμ, 4Nμ) stationary distribution, conditioned on polymorphism, gives the neutral spectrum (∝ 1/i + 1/(2N − i)) when 4Nμ is small and the correct heterozygosity when it is large.

## Findings

**F1 [HIGH, confidence: high] — about 19% of loci are monomorphic and the spectrum is distorted (principle: model correctness)**
`R/gl.sim.Neconst.r:46–64` — frequencies are drawn on a grid of 4N + 1 points (`num_bins <- 2 * (ninds * 2) + 1`, commented "for plotting"), but a population of N diploids has 2N gene copies, and genotypes are then resampled binomially. Frequencies near 1/(4N) often give no copy of the rare allele in the 2N copies.
Failure scenario: `gl.sim.Neconst(20, 20000)` returns 3,749 monomorphic loci (18.7%). Among polymorphic loci, singletons are 18% below the neutral expectation (1,609 against 1,959) and the middle of the spectrum is 6% above it. The documented example computes `gl.sfs()` on this output, and any spectrum-based statistic (Tajima's D, θ estimators, MAF filters) inherits the distortion.
Proposed change: draw each locus's allele count k (1 … 2N − 1) with weight ∝ [k/2N · (1 − k/2N)]^(θ − 1), θ = 4 · `ninds` · `mutation_rate`, and place the k alleles at random among the 2N gene copies. A prototype gives 0 monomorphic loci, singletons at 2,444 against 2,411 expected, the middle of the spectrum at 0.996 of expectation, and heterozygosity at θ = 4 of 0.444 (Beta expectation 0.444).
(Corrected in Phase C: the first draft also claimed 2.2 s against 4.7 s for 100 × 100,000. That timing left out `gl.compliance.check()`, which takes 3.6 s of the run in both versions; the full function takes 4.9 s after the change against 4.7 s before.) **Consequence: output changes: every locus is polymorphic and the spectrum matches the neutral expectation.**

**F2 [MEDIUM, confidence: high] — `mutation_rate` and the title describe a different function (DOC5)**
`R/gl.sim.Neconst.r:1–16` — the title says "constant mutation rate" (the function name says constant Ne). `mutation_rate` is described as the per-generation rate, but only θ = 4 · `ninds` · `mutation_rate` enters, and while θ is much smaller than 1 the output does not depend on it.
Failure scenario: `mutation_rate = 1e-8` and `1e-5` give the same heterozygosity (0.195); only values near 0.01 or above change the output. A user varying `mutation_rate` over realistic values sees no effect and cannot tell why. `ninds` is both the sample and the effective population size, which the docs do not say.
Proposed change: documentation only — title "Simulates a population at mutation–drift equilibrium with constant Ne"; describe θ, the conditioning on polymorphism, that `ninds` is Ne and the sample, and that small θ gives the neutral spectrum whatever `mutation_rate` is.

**F3 [MEDIUM, confidence: high] — inputs are not checked (FS5)**
`R/gl.sim.Neconst.r:26–57`
- `mutation_rate = 0` or `ninds = 0`: "too few positive probabilities".
- `mutation_rate = -1`: "NA in probability vector" with a `dbeta` warning.
- `nlocs = 100.7` is truncated to 100 loci without a message.
Proposed change: stop unless `ninds` is a whole number ≥ 2, `nlocs` a whole number ≥ 1, and `mutation_rate` a single number > 0.

**F4 [LOW, confidence: high] — history and locus metadata (FS8, DAT5)**
`R/gl.sim.Neconst.r:68–75` — the only history entry is `gl.compliance.check(x = inds, verbose = 0)`, so the call that made the data is not recorded. `loc.metrics` has a first column named `array(NA, nLoc(x))` (created by `utils.reset.flags()` inside the compliance check); `ind.metrics` has `id` but no `pop`.
Proposed change: record `match.call()` as the history; build `loc.metrics` with `AlleleID` before the compliance check; add `pop` to `ind.metrics`. Individual names ("1", "2", …), locus names, population name ("pop1") and alleles ("A/C") stay as they are.

**F5 [LOW, confidence: high] — standard skeleton missing; dead code (FS2, FS3, FS9, DOC2)**
`R/gl.sim.Neconst.r:26–80` — `verbose` defaults to 0 and is never used (`verbose = 5` prints nothing); no start/end messages. `L` is computed and never used; `pop_size` and the inner `allele_frequency_distribution()` only rename values; commented `fbm` code remains.
Proposed change: `verbose = NULL` with `gl.check.verbosity()`, `utils.flag.start()`, completion message, a summary at `verbose >= 3`; remove the dead code. **Consequence: with the default verbosity (2) the function prints "Starting…/Completed…" messages; it is silent today.**

**F6 [LOW, confidence: high] — documentation tags (DOC1, DOC7 proposed)**
`R/gl.sim.Neconst.r:1–24` — no `@name`, no `@family`; `verbose` does not use the standard text; parameters do not end with `[required]`/`[default]` consistently; `@author` lacks the Author(s)/Custodian parts.
Proposed change: fix the tags.

## Proposed changes

1. Allele counts drawn on the 2N gene copies, conditioned on polymorphism, with alleles placed at random (F1). **Consequence: output changes; every locus is polymorphic and the spectrum matches the neutral expectation.**
2. Correct title and describe θ, Ne and the role of `mutation_rate` (F2). Documentation only.
3. Validate `ninds`, `nlocs`, `mutation_rate` (F3). **Consequence: fractional `nlocs`, accepted today, stops with an error.**
4. History records this call; `loc.metrics` with `AlleleID`; `ind.metrics$pop` (F4).
5. Standard skeleton; dead code removed (F5). **Consequence: start/end messages print at the default verbosity.**
6. Documentation tags (F6).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run. PLT: not applicable. DEP: `dartR.popgen` is imported for `gl.sfs`, which is used only in the example; not raised.
- Spec: site frequency spectrum against the neutral expectation (20 × 20,000); monomorphic fraction; `mutation_rate` from 1e-8 to 0.05; invalid inputs; metadata; verbosity; timing; prototype of change 1 (spectrum, heterozygosity against the Beta expectation, timing) — run.
- Callers: no internal callers; none in sibling `dartR.*` packages; dartr2shiny holds a source copy only (`input_generator/dartR.sim/gl.sim.Neconst.r`), with no GUI module — read.
- Input-data checks (FS4, DAT1–DAT4): not applicable; the function takes no genlight input.
- FBM path (DAT6): not applicable (the `fbm` argument is commented out).
- Google Group / GitHub issues: not searched.

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | consequence approved |
| 2 | approved | Luis |  |
| 3 | approved | Luis | consequence approved |
| 4 | approved | Luis |  |
| 5 | approved | Luis | consequence approved |
| 6 | approved | Luis |  |

## Addendum notes (not changed)

- `gl.compliance.check()` (dartR.base) takes 3.6 s of the 4.9 s for 100 × 100,000; drawing genotypes takes 0.9 s. Not changed here.
- `utils.check.datatype()` takes SNP data whose genotypes are all 0 or 1 and that has no `loc.all` for SilicoDArT (seen with `ninds = 3, nlocs = 1`). Setting `loc.all` before the compliance check avoids it; the heuristic is in dartR.base (also noted in the `gl.sim.offspring` review).

## Outcome

| Change | Evidence | Snapshot result |
|---|---|---|
| 1 | test "every locus is polymorphic and the spectrum is neutral" | flipped: 0 monomorphic loci in 20 × 20,000 (was 3,749); singletons 2,444/2,365 against 2,411 expected (were 1,609/1,623 against 1,959); middle of the spectrum 0.996 of expectation (was 1.06). Heterozygosity at θ = 4: 0.444 (Beta expectation 4/9). Genotypes are placed by drawing each individual's two copies without replacement (`rhyper()`), so Hardy–Weinberg holds (heterozygotes 0.2353 against 0.2353 expected). Run time 4.9 s for 100 × 100,000 (was 4.7 s) |
| 2 | `man/gl.sim.Neconst.Rd` regenerated | docs only |
| 3 | test "inputs are checked" | flipped: `mutation_rate` 0 or -1, `ninds` 0 or 10.5, `nlocs = 100.7` stop with named errors |
| 4 | test "metadata and history" | flipped: history is the `gl.sim.Neconst()` call (was `gl.compliance.check()`); `loc.metrics$AlleleID` present, no `array(NA, nLoc(x))` column; `ind.metrics` has `id`, `pop`. Names, population and alleles unchanged |
| 5 | same test | flipped: summary at `verbose = 3` (was silent at any level); silent at `verbose = 0` |
| 6 | `devtools::document()`; `@family` adds `gl.sim.Neconst` to the "Other simulation functions" links | docs only |

Tests: `test-gl.sim.Neconst.R` 29 expectations pass; full suite 228 pass. `R CMD check` on the tracked files: 0 errors; 1 warning (installed packages built under R 4.4.3, local environment); 1 note (timestamps). NEWS.md updated.
PR: pending.

```json
{
  "function": "gl.sim.Neconst",
  "package": "dartR.sim",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "cc29b1f",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "principle:model-correctness", "status": "approved", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "approved", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "FS8", "status": "approved", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "FS2", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved", "change": 6}
  ],
  "coverage_skipped": ["input-data checks: no genlight input", "FBM path: fbm argument disabled", "forum/issues search: not run"],
  "status": "pr-open",
  "pr": null
}
```
