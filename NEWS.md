# dartR.sim 1.2.2.9000

## gl.report.nall

* The documentation describes the null model: the curve comes from
  individuals simulated under Hardy-Weinberg from the pooled allele
  frequencies (`gl.sim.ind()`), not from subsampling the real individuals,
  so it need not reach 1 at the full sample size.
* `verbose = 0` is silent (it printed the messages of `gl.filter.allna()` and
  `gl.colors()`).
* SilicoDArT input stops at the start with a clear error (it failed inside a
  parallel worker).
* A single simulation job (one sample size, `reps = 1`) no longer crashes;
  `simlevels`, `reps` and `ncores` are validated (`reps = 0` ran two
  replicates and now stops with an error).
* `ncores = 1` runs without starting a cluster (0.55 s instead of 3.7 s for
  the `possums.gl` example).
* `sim` is returned as a `data.frame`, as documented (it was a tibble).

## gl.sim.Neconst

* Allele counts are drawn on the population's 2N gene copies, conditioned on
  polymorphism, and placed at random among individuals. Before, frequencies
  came from a grid of 4N + 1 points and were resampled, so about 19% of the
  requested SNPs were monomorphic and singletons were 18% short of the
  neutral spectrum. **Output changes: every locus is polymorphic and the
  spectrum matches the neutral expectation.**
* Title and documentation describe the model: theta = 4 * ninds *
  mutation_rate, `ninds` is Ne, and `mutation_rate` changes the output only
  when theta approaches 1.
* Inputs are validated (fractional `nlocs` now stops with an error).
* `verbose` follows the dartR default: start/end messages print at the
  default verbosity. The history records this call; `loc.metrics` holds
  `AlleleID` and `ind.metrics` gains `pop`.
## gl.sim.mutate

* No mutations are applied when none are drawn. Before, a `1:0` loop applied
  1-2 mutations in every such call, so `mut.rate = 0` still mutated and the
  default rate gave `testset.gl` about 1.7 mutations per call instead of
  0.13. **Output changes whenever the draw is 0.**
* SilicoDArT input stops with an error (it got presence/absence scores of 2).
* Locus metrics flags are reset after mutation.
* Only the mutated individual is converted, not the whole matrix: 40
  mutations on 500 individuals x 20,000 loci take 0.04 s instead of 7.9 s.
  Seeded results are unchanged.
* `mut.rate` is validated; new `verbose` argument; history is recorded.
## gl.sim.emigration

* `emi.m` moves individuals from column to row, as documented. Before, the
  direction was reversed: a matrix sending the emigrants of A to B moved B
  individuals into A. **Output of every `emi.m` run changes**; `emi.table`
  runs keep their direction.
* Emigrants are drawn from the residents at the start of the call and all
  move at once, so an individual moves at most once (before, a migrant could
  move again in the same call). **Seeded outputs change**, and the order of
  individuals changes (still grouped by population). About 35 times faster
  (0.39 s -> 0.011 s for the `possums.gl` example).
* A population can lose all its individuals (before, the function stopped).
  Asking for more emigrants than a population holds stops with an error
  naming it.
* `ind.metrics$pop` is updated for moved individuals.
* Inputs are validated: matrices of the wrong dimension and `perc.mig`
  outside 0-1 now stop with an error; data frames are accepted; an all-zero
  `emi.m` column means nobody leaves; a list of length 1 is handled as a
  list, and unnamed list elements are named after their population.
* `perc.mig` is documented as a proportion (it always was one).
* New `verbose` argument; one history entry per call.

## gl.sim.ind.af

* Frequencies and sizes are matched to populations by name. Before, a table
  not sorted alphabetically gave some populations another population's
  frequencies (Z listed before A: Z simulated at 0.10 instead of 0.90), and a
  factor `popn` with named `pop.sizes` swapped the sizes. **Output changes
  for those inputs**; sorted tables such as `gl.allele.freq()` output are
  unaffected.
* `frequency` is documented as the alternate-allele frequency (the allele
  counted in the genotypes), which is what the code has always used.
* Returns a `dartR` object; `loc.metrics` holds `AlleleID` (it had a column
  named `array(NA, nLoc(x))`), and `ind.metrics` gains `id` and `pop`.
* Fractional or `NA` `pop.sizes`, and duplicate population-locus rows, stop
  with an error (they were truncated or dropped silently).
* Genotypes are drawn with `rbinom()` in R instead of C++ compiled at run
  time: no compiler is needed and the first call no longer takes ~3 s.
  **Seeded outputs change**; the genotype distribution is unchanged.
* New `verbose` argument; history is recorded.
## gl.sim.offspring

* Loci are inherited independently. Before, `ifelse()` recycled the random
  draws, so some pairs of loci were always passed on together (a mother
  heterozygous at loci 1 and 11 passed the same allele at both in every
  offspring). **Seeded outputs change**, including the sibs simulated by
  `dartR.captive::gl.sim.relatedness()`.
* Each mother has exactly `noffpermother` offspring (before, mothers were
  drawn at random and could get none). **Family sizes change.**
* Parents must be SNP data with the same loci in the same order; otherwise
  the function stops.
* The offspring are a `dartR` object with the mothers' locus metadata, and
  `@other$ind.metrics` records sex, mother and father. `@other$sex` is kept.

## gl.diagnostics.sim

* The expected FST is now a curve over generations from the exact
  identity-by-descent recursion of the simulated island model. The old line,
  1/(16 Ne m + 1), was about half the correct value (0.059 against 0.111 for
  2 populations of 50 exchanging one individual per generation; simulated
  0.114). **The plotted expectation changes.** Dispersal types "line" and
  "circle", and dispersal files, stop with an error.
* The expected He decays from the first stored generation (it was offset
  when phase 1 was used). **The plotted expectation changes.**
* **Returns a list** (`plot`, `he`, `fst`) instead of the plot alone; use
  `$plot` for the plot.
* Population sizes quoted in the CSV are read; inputs are validated;
  `verbose = 0` is silent.

## gl.sim.ind

* Loci with no calls get `NA` genotypes instead of stopping the function
  (its own example failed on `testset.gl`).
* SilicoDArT input now stops with an error; before, it was returned as
  meaningless diploid genotypes.
* The output keeps the chromosome and locus metrics of the input (metrics
  flags reset).
* Genotypes are drawn with `rbinom()`: about 40 times faster (46 s -> 1.2 s
  for 1000 individuals x 10,080 loci). **Seeded outputs change**; the
  Hardy-Weinberg genotype distribution is unchanged.
* New `verbose` argument; history is recorded; `n` is validated.

## gl.sim.apply (new)

* `gl.sim.apply()` runs any function over the output of `gl.sim.WF.run()`
  and tags each result with its iteration and generation (read from
  `@other$sim.vars$generation`). Data frames and vectors are bound into one
  tidy data frame; other results are returned as a nested, tagged list.

## gl.sim.WF.table

* Arguments passed through `...` no longer put quotes around
  `chromosome_name`; before, any `...` override made `real_loc = TRUE` and
  the map/targets chromosome lookups fail unless `chromosome_name` was also
  passed.
* `real_freq = TRUE` with `real_loc = FALSE` now returns loci typed `"real"`.
  Before, it stopped with an error or returned those loci with `NA` values.
* Recombination maps with intervals of any size are now supported. Each
  interval's own `from`/`to` is used, and neutral loci are spread over the
  whole mapped chromosome. **Output changes** for maps whose intervals are not
  aligned to `chunk_bp` from position 1 (including `fly_recom_map.csv`).
  `NA` in `cM` now marks an interval as not recombining instead of breaking
  the call.
* `loci_deleterious >= chunk_number` now returns exactly the number of loci
  requested (before, 149 gave 100 and 150 gave 200). **Output changes** when
  `loci_deleterious` is not a multiple of `chunk_number`.
* Caps on `q` (0.5) and `s` (0.99, -0.5) now apply to each locus class
  according to that class's own distribution setting, are documented, and
  are reported at `verbose >= 1`. **Output changes** when deleterious and
  advantageous settings differ.
* Names passed through `...` that are not simulation variables now stop the
  function with an error. Before, they were ignored.
* Errors now carry their message (`stop(error(...))`).
* `fly_recom_map.csv` and `fly_targets_of_selection.csv` now ship in
  `inst/extdata`.

## gl.sim.WF.run

* With 10 or more populations, individuals are now labelled with the
  population they belong to. Before, the labels were scrambled: for
  example, individuals of population 10 were labelled "2" and those of
  population 2 were labelled "5". **Population labels change for runs with
  10 or more populations.**

* Recombination now follows the map. Each gamete gets a Poisson number of
  crossovers with mean equal to the map length, independently of its
  siblings. Before, a meiosis had at most one new crossover (none when the
  Poisson draw was 1), and crossovers accumulated from one sibling to the
  next. **Output changes for every run with `recombination = TRUE`.**
* Advantageous selection now acts when `local_adap` is `NULL`. Before, it
  was silently switched off in every population. **Output changes for runs
  with selection and advantageous loci.**
* Clinal adaptation indexes populations by their position in the cline;
  populations outside the cline get advantageous `s = 0`; the multiplier is
  floored at 0. **Output changes for clinal runs.**
* Returned elements are named by the generation they hold. Before, names
  were shifted with `phase1 = TRUE` and were `NA` for iterations after the
  first.
* Extinction is detected at any `verbose` and for each population; the
  iteration stops and the generations stored so far are returned. Before,
  it crashed or overwrote the last stored generation.
* Mutation only uses loci of type `mutation_neu`, `mutation_del` and
  `mutation_adv`; an empty pool is skipped. Before, lost neutral and
  selected loci also received mutations. **Output changes for mutation
  runs.**
* Real allele frequencies (`real_freq = TRUE`) use the alternative allele,
  fill missing frequencies from the pooled frequency, and follow locus
  position. Before, the coding was flipped. **Output changes for
  real-frequency runs.**
* Phase-2 founders are sampled without replacement when phase 1 has enough
  individuals. **Output changes for runs with `phase1 = TRUE`.**
* Each connected pair of populations now exchanges individuals once per
  migration event (it was twice). **Migration rate halves for dispersal
  types `line`, `circle` and `all_connected`.**
* A warning is printed when relative selection is weakened by a small
  offspring pool.
* Names passed through `...` that are not simulation variables now stop the
  function; `local_adap = "1 2"` and `clinal_adap` work through `...`.
* `phase1 = TRUE` works with `real_pops` and `real_pop_size`; the
  interactive route no longer fails on `replace_parents`.
* The pool of mutation loci is reset at the start of each iteration, so
  iterations are independent replicates. **Output changes for mutation runs
  with `number_iterations > 1`.**
* C++ helpers are compiled once per session.
* `@details` now describes the model.

## gl.sim.create_dispersal

* Each connected pair of populations is written once, since
  `gl.sim.WF.run()` swaps individuals in both directions for each row.
  **Migration halves for runs that use a newly generated dispersal file.**
