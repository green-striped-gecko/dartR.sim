# Review: gl.sim.WF.run (dartR.sim)

- Family mode: analysis (forward Wright-Fisher simulator; returns lists of genlight objects). Reviewed together with its engine in `R/utils.sims.r` (`reproduction()`, `recomb()`, `selection_fun()`, `migration()`, `store()`).
- Date: 2026-09-23
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 0f77f7f (`dev_luis`; includes the `gl.sim.WF.table` fixes of PR #47, which this function consumes)
- Datasets: `sim_variables.csv` and `ref_variables.csv` (dartR.sim); `testset.gl` (full, and filtered to call rate 1 for the real-data routes); synthetic parent pairs for recombination; synthetic offspring pools for selection
- Baseline: `tests/testthat/test-gl.sim.WF.run.R` (snapshot captured pre-review, 21 expectations, all pass, ~20 s)

## Verdict

**Standards: Needs work** — errors are raised with an empty message, C++ is recompiled inside the generation loop, the `...` override logic is duplicated and fragile, and several failure paths depend on `verbose`.
**Spec: Rework** — two core parts of the model are wrong under default settings. Recombination does not follow the map: each meiosis has at most one fresh crossover, and crossovers accumulate across siblings. Advantageous selection is switched off whenever `local_adap` is left at `NULL`. Mutation, real-frequency initialisation, generation labels and extinction handling also produce wrong output or crash.

What works: genetic drift. Neutral heterozygosity decays at the rate expected for Ne = N (Ne estimated at 50.1 for N = 50, 20 replicates, 50 generations), and a lower `variance_offspring` lowers Ne as expected (25.0). Deleterious selection acts in the right direction with roughly the right strength. Outputs are diploid genlights with complete locus and individual metrics.

## Population-genetic model, as implemented

For reference, this is the model the code implements (numbers from the runs listed under Coverage):

| Process | Implementation | Assessment |
|---|---|---|
| Initial generation | Each chromosome is drawn allele by allele with P(allele "1") = `q` — linkage and Hardy-Weinberg equilibrium | Correct; matches `@description` |
| Mating | `N/2` monogamous pairs (males = first half of rows, females = second half); family size ~ negative binomial(`mu = number_offspring`, `size = variance_offspring`) | Correct. Ne = N with Poisson families; Ne falls with `variance_offspring` (the NB size parameter: smaller value = more variance). Not documented |
| Next generation | Exactly `N/2` males and `N/2` females sampled without replacement from the offspring pool | Correct; enforces an equal sex ratio |
| Recombination | Per meiosis, events ~ Poisson(λ = ceiling(map length in Morgans)); each event places a chiasma between loci i and i+1 with probability `c[i] / λ`, otherwise none | Design correct (mean crossovers = map length, no interference). **Implementation wrong (F1)** |
| Selection (fitness) | Multiplicative across loci: `1-s` for "11", `1-hs` for heterozygotes, `1` for "00"; advantageous loci have negative `s` | Correct |
| Selection, relative model | Viability: parents sampled without replacement with probability proportional to fitness | Weakens selection when the offspring pool is small (F11) |
| Selection, absolute model | An offspring survives if fitness > U(0, `genetic_load`), so P(survive) = min(1, w / `genetic_load`) | Consistent, but every individual with w ≥ `genetic_load` survives, so there is no selection above that fitness. Not documented (F16) |
| Local / clinal adaptation | Advantageous `s` set to 0 outside `local_adap` populations; multiplied by a decreasing factor along the cline | **Broken (F2, F3)** |
| Migration | Symmetric swaps of `number_transfers` individuals (split between sexes) for every ordered population pair every `transfer_each_gen` generations | Each unordered pair is swapped twice (F10) |
| Mutation | Each offspring gets one new "1" allele at probability `mut_rate`, at a locus drawn from a pool of mutation loci; a locus leaves the pool when mutated and returns when the allele is lost | Pool also takes neutral loci (F6) |

## Findings

**F1 [BLOCKER, confidence: high] — recombination does not follow the map (principle: model correctness)**
`R/utils.sims.r:207–242` (`reproduction()`), with `recomb()` at `:264–282`.
Three defects combine:
(a) A crossover is attempted only when the Poisson draw is `> 1`, so meioses with exactly one event get none.
(b) `for (event in males_recom_events)` loops over the single value of the count, not over `seq_len(count)`, so at most one chiasma is placed per meiosis.
(c) `male_chromosomes` and `female_chromosomes` are overwritten by the recombinant pair and reused for the next sibling, so crossovers accumulate along the family: sibling k carries the breakpoints of siblings 1…k.
Failure scenario: a parent with haplotypes all-0 and all-1, 300 pairs × 10 offspring (the default family size):

| Map length | Expected crossovers per gamete | 1st sibling | 5th | 10th | Mean |
|---|---|---|---|---|---|
| 0.99 M | 0.99 | 0.25 | 1.36 | 2.49 | 1.61 |
| 4.95 M | 4.95 | 0.93 | 4.51 | 8.70 | 5.11 |
| 9.99 M (default table) | 9.99 | 1.00 | 4.98 | 9.85 | 5.96 |

Linkage disequilibrium, haplotype-block structure, Hill–Robertson interference, background selection, and the rate at which neutral loci decouple from selected loci are all wrong. Siblings differ systematically in recombination.
Proposed change: for each gamete, draw k ~ Poisson(λ) and apply k chiasmata (`seq_len(k)`, including k = 1) to the parent's own chromosomes, never to the previous sibling's recombinants. **Consequence: every simulation with `recombination = TRUE` (the default) produces different genotypes; linkage-related statistics change.**

**F2 [BLOCKER, confidence: high] — advantageous selection is switched off when `local_adap` is `NULL` (principle: model correctness)**
`R/gl.sim.WF.run.r:276–277, 717–729` — `gsub('"', "", NULL)` returns `character(0)`, and `as.numeric()` of that is `numeric(0)`, which is not `NULL`. So `!is.null(local_adap)` is always `TRUE`, `pops_local` becomes every population, and advantageous `s` is set to 0 everywhere. The clinal branch can never run.
Failure scenario: 40 advantageous loci, `s = 0.1`, `h = 0.5`, `q = 0.1`, N = 500, 20 generations: expected q = 0.226 (deterministic); simulated 0.098 (unchanged from the start, as under neutrality). With `local_adap = 1` passed explicitly (one population) it reaches 0.179. Deleterious loci are unaffected.
Proposed change: treat a zero-length `local_adap`/`clinal_adap` as `NULL`. **Consequence: every run with `selection = TRUE` and advantageous loci changes; advantageous alleles now respond to selection.**

**F3 [HIGH, confidence: high] — clinal adaptation indexes populations incorrectly (principle: model correctness)**
`R/gl.sim.WF.run.r:743–763` — `clinal_s[y]` uses the population number instead of its position in the cline; `reference_clinal` has one element per cline population but is indexed by every population number.
Failure scenario: 3 populations with `clinal_adap = "2 3"`: population 1 gets population 2's coefficients, population 2 gets `NA` `s` (the fitness is `NA`, so the sampling step crashes), and population 3 is out of bounds. With more than 1/(`clinal_strength`/100) + 1 populations in the cline, the multiplier goes negative and turns advantageous alleles deleterious. This branch is currently unreachable (F2).
Proposed change: index the multiplier by position in the cline, keep one table per population, give populations outside the cline advantageous `s = 0` (as `local_adap` does), and floor the multiplier at 0. **Custodian decision needed:** is `s = 0` outside the cline the intended semantics? **Consequence: clinal runs start to apply selection (today they are neutral for advantageous loci).**

**F4 [HIGH, confidence: high] — returned generation labels are wrong (DOC5)**
`R/gl.sim.WF.run.r:242, 540–543, 911–912, 970–979` — slots are filled by `count_store`, but names are assigned positionally from `gen_store`. With `phase1 = TRUE` and `store_phase1 = FALSE`, the phase-1 generations in `gen_store` are never stored, so names shift. `count_store` is not reset between iterations, so iteration 2 onwards writes past the end of its list.
Failure scenario: phase 1 (10 generations) + phase 2 (10), `every_gen = 5`: names `generation_1, generation_6, generation_11` hold generations 11, 16, 20. `number_iterations = 2`: iteration 2's elements are named `NA`. Any downstream analysis that reads the list names (such as a He-by-generation plot) is mislabelled; only `@other$sim.vars$generation` is correct.
Proposed change: reset `count_store` per iteration and name each element from the generation actually stored.

**F5 [HIGH, confidence: high] — extinction handling crashes or overwrites results (principle: failure path)**
`R/gl.sim.WF.run.r:783–832`
(a) The extinction branch runs only when `verbose >= 2`. At lower verbosity the loop continues into `sample()` and stops with "cannot take a sample larger than the population".
(b) At verbose 2, an extinction in generation 1 stores the founders, which have no parent columns: "column not found: [V5]".
(c) A later extinction writes into slot `count_store` without incrementing it, overwriting the last stored generation: stored generations 1 and 3 were returned as 1 and 4, labelled `generation_1, generation_3`.
(d) The test compares each population's offspring count with the whole `population_size` vector, so with unequal sizes a population can be declared extinct because of another population's size.
Failure scenario: `number_offspring_phase2 = 1`, or `real_pop_size = TRUE` with `testset.gl` (populations of 1–2 individuals): the run crashes at the default `verbose`.
Proposed change: detect extinction regardless of `verbose`; compare each population with its own size; return what was stored so far plus a flag or message at `verbose >= 1`; do not overwrite stored generations.

**F6 [HIGH, confidence: high] — mutation recycles non-mutation loci and crashes when the pool empties (principle: model correctness)**
`R/gl.sim.WF.run.r:686–705, 874–906`
(a) After each generation every locus with allele-"1" count 0 is added to the mutation pool, including neutral, real, deleterious and advantageous loci.
(b) When the pool is empty, the `break` runs only at `verbose >= 2`; otherwise `sample(integer(0), 1)` stops with "invalid first argument".
(c) `sample(pool, 1)` with one locus left draws from `1:pool`, so it mutates an arbitrary locus.
Failure scenario: `q_neutral = 0`, no mutation loci, `mutation = TRUE`, `mut_rate = 1`: after 10 generations 46 of 100 neutral loci are polymorphic, although neutral loci were never meant to mutate. The same run at `verbose = 0` crashes.
Proposed change: recycle only loci whose type is `mutation_neu`, `mutation_del` or `mutation_adv`; skip mutation when the pool is empty regardless of `verbose`; sample by index. **Consequence: mutation runs change; neutral and standing-variation loci no longer receive new mutations.**

**F7 [HIGH, confidence: high] — real allele frequencies are loaded with flipped coding, fail on missing data, and assume sorted loci (DAT5, principle: model correctness)**
`R/gl.sim.WF.run.r:300–311, 455–487`
(a) `gl.alf()$alf1` is the reference-allele frequency (genotype 0), but it is used as the probability of allele "1", which `store()` writes as genotype 2.
(b) A locus with no calls in a population gives an `NA` probability, and the run stops with "NA in probability vector".
(c) With `real_loc = TRUE`, frequencies are taken in the genlight's locus order, but the reference table assigns real loci by sorted position.
Failure scenario: one population of `testset.gl` (call rate 1), 200 simulated individuals: the correlation between real and simulated alternative-allele frequencies is −1.00. `testset.gl[, 1:200]` crashes. A genlight whose loci are not sorted by position gets frequencies on the wrong loci.
Proposed change: use `alf2`; replace `NA` frequencies with the frequency across all populations (and `q_neutral` if still `NA`), reporting the count at `verbose >= 2`; order frequencies by position when `real_loc = TRUE`. **Consequence: real-frequency runs change (the coding flips back).**

**F8 [HIGH, confidence: high] — the interactive route crashes at the first reproduction (principle: failure path)**
`R/gl.sim.WF.run.r:104`; `R/utils.sims.r` `interactive_sim_run()` — `replace_parents` is set to `NULL` and only the CSV defines it; the Shiny app does not return it. `sample(..., replace = NULL)` stops with "invalid 'replace' argument".
Failure scenario: `gl.sim.WF.run(ref_table = rt)` with the default `interactive_vars = TRUE`: the app returns, then the run stops. Not run here (needs a browser); the `sample()` behaviour was checked directly.
Proposed change: default `replace_parents` to `FALSE` when the app does not supply it (and add it to the app).

**F9 [MEDIUM, confidence: high] — the phase 1 → phase 2 founder sampling clones individuals (principle: model correctness)**
`R/gl.sim.WF.run.r:546–574` — phase-2 founders are drawn from phase-1 populations with `replace = TRUE`, even when phase 1 has enough individuals.
Failure scenario: phase 1 N = 100 → phase 2 N = 100: in the first phase-2 generation, 10 of 33 fathers are clones that mated with more than one female (3 of 18 when N2 = 50). Clones inflate relatedness and drift at the founding event beyond what the change in N implies.
Proposed change: sample without replacement when each sex has enough individuals; use replacement only when phase 2 is larger than phase 1. **Consequence: phase-1 runs change.**

**F10 [MEDIUM, confidence: high] — migration swaps each population pair twice (DOC5)**
`R/gl.sim.WF.run.r:591–651` — `all_connected`, `line` and `circle` list both (i, j) and (j, i), and `migration()` already swaps in both directions.
Failure scenario: 3 populations, `all_connected`, `number_transfers = 1`: each population receives 4 immigrants per event (8, 7 and 7 with `number_transfers = 2`). `line`: 2, 4, 2. The documented meaning is "Number of dispersing individuals". Separately, with `real_pop_size = TRUE`, odd real sizes are used unadjusted (e.g. 11 instead of 12), so the female row range starts inside the male half.
Proposed change: process each unordered pair once, document that each migration event swaps `number_transfers` individuals in each direction per connected pair, and use the even-adjusted sizes. **Custodian decision needed** on the intended meaning. **Consequence: the migration rate halves for the built-in dispersal types.**

**F11 [MEDIUM, confidence: high] — relative selection is weakened when the offspring pool is small (principle: model correctness)**
`R/gl.sim.WF.run.r:857–869` — weighted sampling without replacement does not give inclusion probabilities proportional to fitness; as the pool approaches N, everyone is chosen.
Failure scenario: one locus, `s = 0.2`, `h = 0.5`, q = 0.5, N = 2000: expected q after selection 0.476. Pool of 1.2 × N gives 0.493 (60% of the change lost); 5 × N gives 0.479; 20 × N gives 0.472. The default `number_offspring = 10` gives about 5 × N, but `number_offspring` of 2–3 nearly removes selection.
Proposed change: document the dependence and warn at `verbose >= 1` when the pool is under 3 × N. An alternative (a Bernoulli viability step with probability w / max(w), then random sampling) changes results and is left for a separate decision.

**F12 [MEDIUM, confidence: high] — `...` overrides are fragile (API2 proposed rule)**
`R/gl.sim.WF.run.r:153–179` — values are pasted into R code and evaluated. `local_adap` and `clinal_adap` are not quoted, so `local_adap = "1 2"` stops with "unexpected numeric constant". Unknown names are ignored (`gen_number_phase = 3` runs with the default 10), and repeated names misalign the assignment ("replacement has 4 rows, data has 3").
Proposed change: the same approach as `gl.sim.WF.table` in PR #47 — a helper that strips quotes from character variables, evaluates the rest, and stops on unknown names. **Consequence: calls with misspelt `...` names start to error.**

**F13 [MEDIUM, confidence: high] — errors carry no message; the dependency guard returns −1 (FS5, VRB2, DEP1)**
`R/gl.sim.WF.run.r:96–101, 208–227, 389–404` — `message(error(...)); stop()` produces an empty condition message. The `stringi` guard returns `-1` instead of stopping; `stringi` is already in Imports, so the guard is dead code.
Failure scenario: a mismatch between `real_freq` in the table and in the run gives "Error:" inside `tryCatch` or Shiny.
Proposed change: `stop(error(...))`; remove the dead guard.

**F14 [MEDIUM, confidence: high] — C++ is recompiled inside the generation loop (STY2, DEP3 proposed rule)**
`R/gl.sim.WF.run.r:417, 883`; `R/utils.sims.r:307, 393` — `Rcpp::cppFunction()` runs for every iteration (`make_chr`), every generation (`make_freqs`), every generation × population (`make_fit` in `selection_fun()`) and every stored generation (`make_geno`). Even with the compile cache, each call costs ~0.07 s. `reproduction()` grows its data frame with `rbind` in a loop, and the real-frequency initialisation rebuilds the same probability list for every individual.
Failure scenario: 10 populations × 1000 generations with selection spends ~12 minutes in `cppFunction()` lookups alone.
Proposed change: compile each C++ function once per `gl.sim.WF.run()` call and pass it down; build offspring lists and bind once; hoist the loop-invariant probability list.

**F15 [LOW, confidence: high] — `phase1` with `real_pops` / `real_pop_size` is broken (principle: failure path)**
`R/gl.sim.WF.run.r:284–289, 338–345, 389–404` — `real_pops` updates `number_pops` but not `number_pops_phase1`, so the consistency check stops the run. With `real_pop_size`, the phase-1 branch sets `population_size_phase1` but never `population_size`.
Failure scenario: `phase1 = TRUE, real_pops = TRUE, real_pop_size = TRUE` stops with "Number of entries for population sizes do not agree…".
Proposed change: set the phase-1 counts and sizes from `x` as for phase 2.

**F16 [LOW, confidence: high] — documentation gaps and minor output details (DOC1, DOC2, DOC5, DOC6 and DOC7 proposed)**
`R/gl.sim.WF.run.r:1–55` and minor code details:
- No `@details`: the model assumptions in the table above are undocumented (monogamy, equal sex ratio, viability selection on offspring, absolute-model survival = min(1, w/`genetic_load`), `variance_offspring` as the NB size parameter, migration semantics, the mutation model).
- `@return` does not describe the structure (a list per iteration of genlights per stored generation, with `@other$sim.vars`).
- The wiki URL (line 13) points to the retired dartR repository; the curly apostrophe on line 7 (DOC6); the verbose text differs from DOC2; `@author` has no `Author(s):` part (DOC7); `@param x` says "Name of".
- `sample_percent` rounds up to an even number (50% of 50 stores 26).
- `density_mutations_per_cm` (stored as `del_ind_cM`) counts advantageous and mutation loci as deleterious.
- `gen_store` duplicates the last generation when `(generations − 1) %% every_gen == 0`.
Proposed change: add `@details` describing the model, fix the items above, and document the rounding.

**F17 [INFO, confidence: high] — housekeeping (FS3, STY1)**
`R/gl.sim.WF.run.r:87–89, 115–179, 330` — outdated `build = "Jody"`; quote/eval blocks triplicated; `iteration %% 1 == 0` is always true.
Proposed change: drop `build=`; covered by change 12's helper.

## Proposed changes

1. Recombination: k ~ Poisson(λ) chiasmata per gamete, including k = 1, applied to the parent's own chromosomes (F1). **Consequence: all runs with `recombination = TRUE` change; linkage statistics change.**
2. Treat an empty `local_adap`/`clinal_adap` as `NULL` (F2). **Consequence: advantageous loci respond to selection in all runs with selection on.**
3. Fix the clinal indexing; populations outside the cline get advantageous `s = 0`; floor the multiplier at 0 (F3). **Consequence: clinal runs start to apply selection.**
4. Reset `count_store` per iteration and name elements from the stored generation (F4).
5. Extinction: detect at any `verbose`, compare per population, keep stored generations, report at `verbose >= 1` (F5).
6. Mutation: recycle only mutation-type loci; skip when the pool is empty; sample by index (F6). **Consequence: mutation runs change.**
7. Real frequencies: use `alf2`, fill `NA` from the pooled frequency then `q_neutral`, order by position (F7). **Consequence: real-frequency runs change.**
8. Default `replace_parents = FALSE` when not supplied; add it to the app (F8).
9. Founder sampling without replacement when possible (F9). **Consequence: phase-1 runs change.**
10. Migration: one swap per unordered pair, document the migrants per event, use even-adjusted real sizes (F10). **Consequence: the migration rate halves for built-in dispersal types.**
11. Document the pool-size dependence of relative selection and warn below 3 × N (F11).
12. Shared helper for variables and `...` overrides; unknown names stop (F12, F17). **Consequence: misspelt `...` names start to error.**
13. `stop(error(...))`; remove the dead `stringi` guard (F13).
14. Compile C++ once per call; bind offspring once; hoist the invariant probability list (F14).
15. Phase 1 with `real_pops`/`real_pop_size` (F15).
16. Documentation (model `@details`, `@return`, URL, DOC2, DOC6, DOC7) and minor output details (F16, F17).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run. PLT: not applicable.
- Spec, empirical (with `devtools::load_all()`): recombination (direct calls to `reproduction()`, 3 map lengths); selection (a 20-generation run against the deterministic expectation; a single-generation pool experiment for sampling intensity; the absolute-model survival curve); drift (He decay, 20 replicates × 2 offspring variances); migration (the dispersal loop replicated around `migration()`, 3 populations, 2 types × 2 transfer numbers); mutation (pool recycling, empty pool); phase switch (founder clones); generation labels; extinction (generation 1 and later); real frequencies (coding, `NA`); `...` overrides.
- Interactive Shiny route: SKIPPED as a run — needs a browser; F8 is from reading the code plus the `sample(replace = NULL)` behaviour.
- `file_dispersal` route and `gl.sim.create_dispersal()` semantics: read, not run. Whether a user dispersal file is meant to be directional affects change 10.
- `gl.sim.WF.table` inputs: reviewed separately (PR #47).
- Performance at scale: only the per-call `cppFunction()` overhead was measured (0.073 s); no profiling at large N.
- FBM path (DAT6): not applicable — `x` is only read for `pop()`, locus names and frequencies.
- Google Group / GitHub issues: not searched in this pass.

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis |  |
| 2 | approved | Luis |  |
| 3 | approved | Luis | populations outside the cline get advantageous s = 0 (as proposed) |
| 4 | approved | Luis |  |
| 5 | approved | Luis |  |
| 6 | approved | Luis |  |
| 7 | approved | Luis |  |
| 8 | approved | Luis |  |
| 9 | approved | Luis |  |
| 10 | approved | Luis | one swap per unordered pair (as proposed) |
| 11 | approved | Luis | documentation + warning only |
| 12 | approved | Luis | strict: unknown names error |
| 13 | approved | Luis |  |
| 14 | approved | Luis |  |
| 15 | approved | Luis |  |
| 16 | approved | Luis |  |
| F18 | approved | Luis | addendum |
| F19 | approved | Luis | addendum; touches gl.sim.create_dispersal |

## Addendum findings (found while applying; approved and applied)

**F18 [MEDIUM, confidence: high] — the mutation pool is not reset between iterations (principle: model correctness)**
`R/gl.sim.WF.run.r` — `mutation_loci_location` is built once before the iteration loop and consumed inside it, so iteration 2 starts with the pool that iteration 1 left behind. Iterations are not independent replicates for mutation runs.
Proposed change: reset the pool at the start of each iteration. **Consequence: mutation runs with `number_iterations > 1` change from iteration 2 onwards.**

**F19 [MEDIUM, confidence: high] — `file_dispersal` tables from `gl.sim.create_dispersal()` still swap each pair twice (DOC5)**
`R/gl.sim.create.dispersal.r` writes both (i, j) and (j, i). Change 10 was applied to the built-in dispersal types only, because a user file may list directions on purpose.
Proposed change: have `gl.sim.create_dispersal()` write each unordered pair once and document that a row means a two-way swap. **Consequence: migration halves for runs that use a generated dispersal file.** This touches another function, so it may belong in its own review.

**F20 [HIGH, confidence: high] — population labels scrambled with 10 or more populations (found in the gl.diagnostics.sim review as F8; approved there as a separate PR)**
`R/utils.sims.r` `store()` — `pop()` was set from character labels, so the factor levels sorted as "1", "10", "11", "2", …; `gl.sim.WF.run()` then renamed the levels by position.
Failure scenario: 12 populations without dispersal: 220 of 230 individuals carried another population's label (population 10 labelled "2", 2 labelled "5", …).
Change: the factor is built with `levels = as.character(p_vector)`. Test "population labels are right with 10 or more populations" (fails before, passes after); `real_pops = TRUE` with 12 populations of `testset.gl` labels every individual with its birth population's name. **Consequence: labels change for runs with 10 or more populations.**

## Addendum findings (post-review, from the gl.sim.create_dispersal review)

Found on 2026-09-24 while reviewing `gl.sim.create_dispersal` (PR #58); reviewed on `origin/dev` 085a0b2. Baseline: three tests in `tests/testthat/test-gl.sim.WF.run.R` marked (F21), which run the current `migration()` through a copy of the dispersal loop.

**F21 [HIGH, confidence: high] — which sexes migrate depends on the phase setting and on the other rows, not on each row (principle: model correctness)**
`R/gl.sim.WF.run.r:403–410, 515–520, 600–620`; `R/utils.sims.r:163–230` — two switches (move males, move females) are set once per phase from the phase's `number_transfers` and carried from row to row; `migration()` flips both after every row with 1 transfer, and moves one male and one female when both are on.
Failure scenario (10 individuals per population, immigrants counted by sex):
- Phase set to 1, a dispersal-file row of 3: 2 males and 0 females move per direction each generation.
- Phase set to 2, rows of 1, 2, 2: in generation 2 the row of 1 moves a male and a female and the two rows of 2 move nobody; in generation 3 the row of 1 moves nobody and the others move 2 each.
- Every row at 1 with two pairs (3 populations in a line, no file needed): pair 1 swaps a male and pair 2 a female in every generation (4 of 4). With an even number of pairs each pair always swaps the same sex; sexes alternate only with an odd number.
Proposed change: decide sexes per row from the row's own value. 2 or more: ceil(n/2) males and floor(n/2) females, as now. 1: one individual, alternating male/female from one dispersal event to the next for that pair (a per-row state kept across generations and reset at each phase). 0: no transfer. Remove the phase-level switches. Move the dispersal loop into an internal helper so it can be tested directly. **Consequence: seeded outputs change for runs with one transfer per pair and more than one pair, and for dispersal files whose values differ from the phase setting. Runs with one pair, or with 2 or more transfers everywhere, are unchanged.**

| Change | Decision | By | Note |
|---|---|---|---|
| F21 | approved | Luis | consequence approved (2026-09-24) |

F21 outcome: `migration()` takes the row's `n_transfer` and a per-pair `next_male` state; the dispersal loop is the internal helper `dispersal_event()` (`R/utils.sims.r`); the phase-level switches are removed. Tests (F21) flipped: a row of 3 moves 2 males and 1 female (was 2 and 0); rows of 1, 2, 2 move 1, 2, 2 in every generation (was 2, 0, 0 then 0, 2, 2); a row of 0 moves nobody; two pairs of 1 each alternate male, female, male, female (were fixed at male and female); no transfer outside dispersal generations. All earlier `gl.sim.WF.run` tests unchanged (single-pair seeded outputs identical). End to end: 4 populations with a file of 0, 1, 2, 3, 1, 2 runs; the pair set to 0 exchanges no fathers. Tests: `test-gl.sim.WF.run.R` 48 expectations pass; full suite 347 pass. PR: pending.

## Outcome

| Change | Evidence | Snapshot result |
|---|---|---|
| 1 | test "crossovers follow the map"; `gamete()` replaces repeated `recomb()` calls (all chiasmata drawn at once, pairs at the same interval cancel) | flipped. 0.99 M map: mean 1.00 crossovers (1st sibling 0.99, 10th 1.01; was 0.25 → 2.49); P(0) 0.373 vs Poisson 0.372; r between loci 1 and 50 = 0.316 vs Haldane 0.312. 9.99 M: 9.93 (was 5.96) |
| 2 | test "advantageous selection acts" | flipped: q = 0.210 after 20 generations (expected 0.226; was 0.098) |
| 3 | test "clinal adaptation along part of the populations" | `clinal_adap = "2 3"` with 3 populations runs |
| 4 | test "generation labels match stored generation" | flipped: `generation_11/16/20` with phase 1; iteration 2 named `generation_1/10` |
| 5 | test "extinction keeps stored generations" | flipped: no crash at `verbose = 0`; message at `verbose >= 1` |
| 6 | test "only mutation loci mutate" | flipped: 0 of 100 neutral loci polymorphic (was 46); pool-empty runs at `verbose = 0` |
| 7 | test "real_freq uses the alternative allele" | flipped: correlation +1.00 (was −1.00); data with no-call loci run |
| 8 | test "CSV without replace_parents runs" | default applied. Adding a `replace_parents` input to the Shiny app was NOT done: it needs a browser to test; the default covers the crash |
| 9 | test "phase-2 founders are not cloned" | flipped: 0 fathers with more than one mate (was 10 of 33) |
| 10 | test "each connected pair migrates once" | flipped: ~2 immigrant parents per population per generation (was 4) |
| 11 | test "weak relative selection is reported" | new warning |
| 12 | test "... overrides" | `local_adap = "1 2"` works; `gen_number_phase = 3` errors |
| 13 | test "errors carry their message" | new |
| 14 | C++ compiled once per session (`utils.wf.cpp()`); offspring bound once; real-frequency draw vectorised | default example 0.56 s (0.42 s before the review; 1.09 s with correct recombination before `gamete()`) |
| 15 | test "phase 1 with real populations and sizes" | runs |
| 16 | `devtools::document()`; `man/gl.sim.WF.run.Rd` regenerated | model `@details`, `@return`, DOC2/DOC6/DOC7 fixed. Code details listed in F16 (sample rounding, `del_ind_cM`, duplicate last generation) documented or left as is; not changed |

| F18 | test "mutation pool is reset for each iteration" | iterations 2-6 average 4.4 polymorphic mutation loci in generation 1 (2.0 without the reset) |
| F19 | test "gl.sim.create_dispersal writes each pair once"; `@param number_transfers` rewritten | 3 populations, all_connected: 3 rows (was 6) |

Unchanged: default-example structure test and the drift test (Ne = N) still pass. `recomb()` is no longer called but is kept.
Tests: `test-gl.sim.WF.run.R` 36 expectations and `test-gl.sim.WF.table.R` 28 pass. Callers: `dartR.captive` (`utils.classes.diagnostics.relatedness.r`) passes `number_iterations`, `every_gen`, `sample_percent`, `gen_number_phase2`, all valid; dartr2shiny passes CSV-based arguments. NEWS.md updated.
PR: #48.

```json
{
  "function": "gl.sim.WF.run",
  "package": "dartR.sim",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "0f77f7f",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "BLOCKER", "confidence": "high", "rule": "principle:model-correctness", "status": "approved", "change": 1},
    {"id": "F2", "severity": "BLOCKER", "confidence": "high", "rule": "principle:model-correctness", "status": "approved", "change": 2},
    {"id": "F3", "severity": "HIGH", "confidence": "high", "rule": "principle:model-correctness", "status": "approved", "change": 3},
    {"id": "F4", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 4},
    {"id": "F5", "severity": "HIGH", "confidence": "high", "rule": "principle:failure-path", "status": "approved", "change": 5},
    {"id": "F6", "severity": "HIGH", "confidence": "high", "rule": "principle:model-correctness", "status": "approved", "change": 6},
    {"id": "F7", "severity": "HIGH", "confidence": "high", "rule": "DAT5", "status": "approved", "change": 7},
    {"id": "F8", "severity": "HIGH", "confidence": "high", "rule": "principle:failure-path", "status": "approved", "change": 8},
    {"id": "F9", "severity": "MEDIUM", "confidence": "high", "rule": "principle:model-correctness", "status": "approved", "change": 9},
    {"id": "F10", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 10},
    {"id": "F11", "severity": "MEDIUM", "confidence": "high", "rule": "principle:model-correctness", "status": "approved", "change": 11},
    {"id": "F12", "severity": "MEDIUM", "confidence": "high", "rule": "API2", "status": "approved", "change": 12},
    {"id": "F13", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "approved", "change": 13},
    {"id": "F14", "severity": "MEDIUM", "confidence": "high", "rule": "STY2", "status": "approved", "change": 14},
    {"id": "F15", "severity": "LOW", "confidence": "high", "rule": "principle:failure-path", "status": "approved", "change": 15},
    {"id": "F16", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved", "change": 16},
    {"id": "F17", "severity": "INFO", "confidence": "high", "rule": "FS3", "status": "approved", "change": 16},
    {"id": "F18", "severity": "MEDIUM", "confidence": "high", "rule": "principle:model-correctness", "status": "approved", "change": "addendum"},
    {"id": "F19", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": "addendum"}
  ],
  "coverage_skipped": ["interactive Shiny route: needs browser", "file_dispersal route: read only", "large-N profiling: not run", "forum/issues search: not run"],
  "status": "pr-open",
  "pr": 48
}
```
