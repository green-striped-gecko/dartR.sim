#' @name gl.sim.WF.run
#' @title Runs Wright-Fisher simulations
#' @family simulation functions
#' @description
#' This function simulates populations made up of diploid organisms that 
#' reproduce in non-overlapping generations. Each individual has a pair of 
#' homologous chromosomes that contains interspersed selected and neutral loci. 
#' For the initial generation, the genotype for each individual's chromosomes is
#' randomly drawn from distributions at linkage equilibrium and in 
#' Hardy-Weinberg equilibrium. 
#' 
#' The simulations will take a little longer the first time you use the
#' function gl.sim.WF.run() in an R session because C++ functions must be
#' compiled.
#' @param file_var Path of the variables file 'sim_variables.csv' (see details) 
#' [required if interactive_vars = FALSE].
#' @param ref_table Reference table created by the function 
#' \code{\link{gl.sim.WF.table}} [required].
#' @param x Genlight object containing the SNP data to extract
#' values for some simulation variables (see details) [default NULL].
#' @param file_dispersal Path of the file with the dispersal table created with
#' the function \code{\link{gl.sim.create_dispersal}} [default NULL]. 
#' @param number_iterations Number of iterations of the simulations [default 1].
#' @param every_gen Generation interval at which simulations should be stored in
#' a genlight object [default 10].
#' @param sample_percent Percentage of individuals, from the total population, 
#' to sample and save in the genlight object every generation. The number is
#' rounded up to an even number, e.g. 50 percent of 50 individuals stores 26
#' [default 50].
#' @param store_phase1 Whether to store simulations of phase 1 in genlight
#' objects [default FALSE].
#' @param store_founders Whether to store the founders of the first phase
#' simulated, before they reproduce, as "generation_0" [default FALSE].
#' Their ind.metrics hold F_founder, the proportion of loci whose two copies
#' are identical by descent (0 unless real_inbreeding = TRUE), and
#' F_founder_map, the same proportion measured along the map. Founders
#' have no parents (pat and mat are NA).
#' @param store_pedigree Whether to return the pedigree of every individual
#' of every generation, sampled or not, as the attribute "pedigree" of each
#' iteration: a data frame with id (as in indNames of the stored genlights),
#' pat, mat (NA for founders), generation (0 for founders), pop, F_founder
#' (founders only, NA otherwise) and in_population (FALSE for offspring
#' sampled by sample_families that did not join the population). Read it with
#' attr(res[["iteration_1"]], "pedigree") [default FALSE].
#' @param interactive_vars Run a shiny app to input interactively the values of
#' simulations variables [default TRUE].
#' @param seed Set the seed for the simulations. This calls set.seed(), so it
#' also sets the random number stream of the R session for any code run
#' afterwards [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' brief progress messages; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @param ... Any simulation variable of 'sim_variables.csv' and its value,
#' e.g. gen_number_phase2 = 20 or local_adap = "1 2". The value replaces the
#' one in the csv file or the Shiny app. Names that are not simulation
#' variables stop the function with an error.
#' @details
#' Values for the simulation variables can be entered in a Shiny app
#' (interactive_vars = TRUE) or in the file 'sim_variables.csv', which can be
#' found by typing in the R console:
#' system.file('extdata', 'sim_variables.csv', package = 'dartR.sim').
#' 
#' The model, generation by generation:
#' \itemize{
#' \item Migration. Every transfer_each_gen generations, each connected pair
#' of populations swaps number_transfers individuals in each direction (half
#' males, half females; a single transfer alternates between sexes). With
#' dispersal_type "line", "circle" or "all_connected", each pair of
#' connected populations is processed once. A file from
#' \code{\link{gl.sim.create_dispersal}} is used row by row.
#' \item Reproduction. N/2 monogamous pairs are formed at random, except that
#' a proportion sib_mating of the pairs are full siblings (see Inbreeding
#' below). The number
#' of offspring per pair follows a negative binomial distribution with mean
#' number_offspring and size variance_offspring: large values of
#' variance_offspring give Poisson family sizes (Ne close to N), small values
#' increase the variance of family size and lower Ne. Sex is assigned at
#' random.
#' \item Recombination. For each gamete, the number of recombination events
#' is Poisson with mean equal to the map length in Morgans rounded up; each
#' event places a crossover between two loci with probability proportional
#' to the recombination rate c of the reference table, so the mean number of
#' crossovers per gamete equals the map length. There is no interference.
#' recombination_males = FALSE switches recombination off in males.
#' \item Mutation. Each offspring receives, with probability mut_rate, one new
#' allele at a locus drawn from the loci of type mutation_neu, mutation_del
#' or mutation_adv that are not segregating. A locus returns to this pool
#' when its new allele is lost.
#' \item Selection. Fitness is multiplicative across loci: 1 - s for the
#' homozygote of the simulated allele, 1 - hs for heterozygotes and 1
#' otherwise; advantageous alleles have negative s. With the "relative"
#' model, the next generation is sampled from the offspring without
#' replacement with probability proportional to fitness; this weakens
#' selection when the offspring pool is not much larger than N (a warning is
#' printed at verbose >= 1 when it is less than 3 times N). With the
#' "absolute" model, an offspring survives with probability
#' min(1, fitness / genetic_load), so there is no selection among offspring
#' whose fitness is at least genetic_load. local_adap lists the populations
#' where advantageous alleles are under selection (s = 0 elsewhere).
#' clinal_adap gives the first and last populations of a cline along which
#' advantageous s is multiplied by 1, 1 - k, 1 - 2k, ... with
#' k = clinal_strength / 100 (floored at 0); outside the cline s = 0.
#' \item Next generation. N/2 males and N/2 females are sampled from the
#' offspring. If there are too few of either sex, the population is extinct:
#' the iteration stops and the generations stored so far are returned.
#' }
#' With phase1 = TRUE, phase 2 starts from individuals sampled from the
#' phase-1 populations, without replacement unless phase 2 is larger.
#' 
#' If a genlight object is used (real_pops, real_pop_size, real_loc,
#' real_freq or real_inbreeding), the simulated allele is the alternative
#' allele of the genlight (genotype 2). Loci without calls in a population
#' take the frequency across all populations. real_pop_size sets the sizes
#' of the first phase simulated, rounded up to even numbers.
#' 
#' Inbreeding. sib_mating_phase1 and sib_mating_phase2 set the proportion of
#' pairs that are full siblings (brother and sister with the same father and
#' mother), one value or one per population (space delimited and in
#' quotes); NULL or 0 is random mating. At equilibrium the inbreeding
#' coefficient is F = b / (4 - 3b) for a proportion b of sib matings. Sib
#' matings need females with a brother among the parents, so when families
#' are small (low number_offspring) fewer sib matings than asked for can take
#' place and a warning is printed at verbose >= 1. With real_inbreeding =
#' TRUE, F = 1 - Ho / He is estimated in each population of x (on the loci
#' of chromosome_name when real_loc = TRUE; negative values are set to 0)
#' and x must have as many populations as the first phase simulated. The
#' founders are then made inbred: along the chromosome they alternate
#' between stretches whose two copies are identical by descent (a fraction F
#' of the map, mean length 25 cM, as in the offspring of full siblings) and
#' stretches drawn independently. In a phase whose sib_mating is NULL, the
#' proportion of sib matings is set to b = 4F / (1 + 3F), which keeps F in
#' the following generations; a sib_mating value that is set is used
#' instead. Founders are unrelated to each other, so the first generation
#' (their offspring) has F close to 0; F returns close to the value in x from
#' the second generation. To store only generations at equilibrium, use
#' phase 1 as a burn-in (phase1 = TRUE, store_phase1 = FALSE).
#' inbreeding_founders (one value or one per population, in the order of
#' the populations) gives the founders' F directly and replaces the value
#' estimated from x, with or without real_inbreeding; negative values are
#' set to 0. It is useful because heterozygote dropout inflates F in the
#' full locus set: estimate F on high call-rate loci and pass it here while
#' x keeps all loci for the allele frequencies.
#' 
#' Differentiation. With real_freq = TRUE, each population's founders are
#' drawn from its sample frequencies in x, which carry sampling noise and
#' so inflate FST among founders. real_freq_shrink = "auto" shrinks them
#' toward their mean, q' = mean + lambda (q - mean), with lambda chosen so
#' the between-population variance of the frequencies equals that of x after
#' removing the variance expected from sampling; an FST estimator that
#' corrects for sample size (e.g. Hudson's) then gives about the same value
#' for the founders as for x, whatever the population size. A number sets
#' lambda directly; NULL (default) does not shrink. The lambda used is
#' stored in sim.vars$freq_shrink_lambda. Drift after the founders raises
#' FST again, faster in small populations.
#' 
#' Migration. With real_migration = TRUE, the number of individuals each
#' pair of populations swaps per generation is set from x's FST (Hudson's,
#' mean over pairs of populations). Each pair swaps number_transfers
#' individuals in each direction, so with n populations and dispersal_type
#' "all_connected" the island model gives an equilibrium
#' FST = 1 / (1 + 4 n T) for T transfers per generation (Takahata, 1983)
#' when Ne = N. FST depends on Ne m rather than N m, so
#' T = (1 / FST - 1) / (4 n) x N / Ne (mean over populations of N / Ne,
#' with Ne the expected Ne, see below) replaces number_transfers. T can be
#' fractional: each pair moves floor(T) or floor(T) + 1 individuals, the latter with probability T - floor(T). It
#' applies to the phases with dispersal, which must be "all_connected" (no
#' dispersal file). T is stored in sim.vars$migrants_real.
#' 
#' Effective population size. Ne = 4N / (Vk + 2), with Vk the variance in
#' the number of offspring an individual leaves in the next generation.
#' With parents sampled without replacement Ne / N = k / (k + 1), and with
#' replace_parents = TRUE (parents mate a Poisson number of times, which
#' produces half siblings) Ne / N = k / (2k + 1), at most 1/2; k is
#' variance_offspring, one value or one per population. ne_phase1 and
#' ne_phase2 set a target Ne for each population, for example estimated
#' from x with gl.LDNe() (per population; with small samples its estimate
#' can be infinite and is then of no use), and variance_offspring is set to
#' reach it. Family-size variance can only lower Ne, so a target above N
#' (N/2 with replace_parents = TRUE) is capped there with a warning; raise
#' the population size instead. With real_pop_size = TRUE the populations
#' are as small as x's samples; to simulate larger populations while storing
#' samples of x's sizes, set population_size_phase2 and use
#' real_sample_size = TRUE, which replaces sample_percent. The expected Ne is
#' stored in sim.vars$ne_expected, which \code{\link{gl.diagnostics.sim}}
#' uses by default.
#' 
#' Family-structured samples. A sample made of a few families (clutches,
#' litters) makes unrelated pairs look related, because allele frequencies
#' come from few families. sample_families gives the sizes of the full-sib
#' families in each stored sample (e.g. "9 15 9"; ";" between
#' populations, one group for all populations). For each size a mating
#' pair with that many offspring is chosen and that many of its offspring
#' are sampled from the offspring pool, before the next generation is drawn
#' (so, like hatchlings, they may not survive to reproduce); the rest of the
#' sample is drawn at random from the population. sample_parents = TRUE
#' also stores each family's two parents, from the previous generation;
#' they count against the sample size. The sample size is x's with
#' real_sample_size = TRUE (exactly, with or without families), otherwise
#' set by sample_percent. Family sampling starts at generation 1, as the
#' founders have no parents; number_offspring must be at least the largest
#' family, otherwise large families are rare and a warning is printed.
#' ind.metrics gains sample_role ("family", "parent" or "random") and
#' family (the parents' ids).
#' 
#' Planted crosses. To recreate a mating design with shared parents (e.g. 4
#' litters from 2 sires x 2 dams), sample_crosses gives one "SxD" token per
#' family of sample_families, where the same label is the same individual
#' (e.g. "S1xD1 S1xD2 S2xD1 S2xD2"; ";" between populations);
#' sample_design = "2x2" is shorthand for every cross of 2 sires and 2
#' dams. Sires and dams are drawn at random from the previous generation and
#' each cross makes a family of exactly the given size, with the usual
#' meiosis. The planted offspring are only sampled: they do not join the
#' population, so drift and gene flow are unchanged. With sample_parents =
#' TRUE each parent is sampled once. ind.metrics also gains sample_cross
#' (the offspring's cross) and parent_label (the parent's label).
#' @return A list with one element per iteration ("iteration_1", ...). Each
#' is a list of genlight objects, one per stored generation, named by the
#' generation they hold ("generation_1", ...). Each genlight has the
#' reference table in @other$loc.metrics, sex and parents in
#' @other$ind.metrics, and the simulation variables (including the
#' generation) in @other$sim.vars.
#' @author Author(s): Luis Mijangos. Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' ref_table <- gl.sim.WF.table(file_var=system.file("extdata", 
#' "ref_variables.csv", package = "dartR.sim"),interactive_vars = FALSE)
#' 
#' res_sim <- gl.sim.WF.run(file_var = system.file("extdata",
#'  "sim_variables.csv", package ="dartR.sim"),ref_table=ref_table,
#'  interactive_vars = FALSE)
#' @seealso \code{\link{gl.sim.WF.table}}
#' @import stats
#' @import shiny
#' @export

gl.sim.WF.run <- function(file_var,
                          ref_table,
                          x = NULL,
                          file_dispersal = NULL,
                          number_iterations = 1,
                          every_gen = 10,
                          sample_percent = 50,
                          store_phase1 = FALSE,
                          store_founders = FALSE,
                          store_pedigree = FALSE,
                          interactive_vars = TRUE,
                          seed = NULL,
                          verbose = NULL,
                          ...) {
  
    
    # -------------------------------
    # SET SEED FOR REPRODUCIBILITY
    # -------------------------------
    if (!is.null(seed)) {
      set.seed(seed)
    }
    
    # -------------------------------
    # SET VERBOSITY LEVEL FOR MESSAGES
    # -------------------------------
    verbose <- gl.check.verbosity(verbose)
    
    # -------------------------------
    # FLAG THE START OF THE FUNCTION EXECUTION
    # -------------------------------
    funname <- match.call()[[1]]
    utils.flag.start(func = funname,
                     verbose = verbose)
    
    # -------------------------------
    # CHECK INPUTS
    # -------------------------------
    if (interactive_vars == FALSE &&
        (missing(file_var) || !file.exists(file_var))) {
      stop(error("  When interactive_vars = FALSE, file_var must be the path",
                 "to an existing 'sim_variables.csv' file\n"))
    }
    if (!is.null(x) && !is(x, "genlight")) {
      stop(error("  x must be a genlight object\n"))
    }
    
    # -------------------------------
    # RETRIEVE SIMULATION VARIABLES
    # -------------------------------
    # Variables come from the Shiny app or from the CSV file; either way
    # sim_vars is a two-column table (variable, value) of character values
    if (interactive_vars == TRUE) {
      sim_vars <- interactive_sim_run()
    } else {
      sim_vars <- suppressWarnings(read.csv(file_var))
      sim_vars <- sim_vars[, c("variable", "value")]
    }
    ## Variables absent from the Shiny app (replace_parents) or from
    ## older sim_variables.csv files take their
    ## defaults: parents sampled without replacement, random mating and no
    ## inbreeding from the real data
    defaults <- c(replace_parents = "FALSE", sib_mating_phase1 = "NULL",
                  sib_mating_phase2 = "NULL", real_inbreeding = "FALSE",
                  real_freq_shrink = "NULL", real_migration = "FALSE",
                  ne_phase1 = "NULL", ne_phase2 = "NULL",
                  real_sample_size = "FALSE",
                  inbreeding_founders = "NULL",
                  sample_families = "NULL", sample_parents = "FALSE",
                  sample_crosses = "NULL", sample_design = "NULL")
    missing_vars <- setdiff(names(defaults), sim_vars$variable)
    if (length(missing_vars) > 0) {
      sim_vars <- rbind(sim_vars,
                        data.frame(variable = missing_vars,
                                   value = unname(defaults[missing_vars])))
    }
    sim_vars <- sim_vars[order(sim_vars$variable), ]
    
    # -------------------------------
    # OVERRIDE VARIABLES WITH ADDITIONAL ARGUMENTS (if any)
    # -------------------------------
    # Names that are not simulation variables are refused, so a misspelt
    # name cannot be silently ignored
    input_list <- list(...)
    if (length(input_list) > 0) {
      unknown <- setdiff(names(input_list), sim_vars$variable)
      if (is.null(names(input_list)) || any(names(input_list) == "") ||
          anyDuplicated(names(input_list)) > 0 || length(unknown) > 0) {
        stop(error("  Arguments passed through ... must be named simulation",
                   "variables, each given once. Unknown:",
                   paste(unknown, collapse = ", "), "\n"))
      }
      for (var in names(input_list)) {
        sim_vars[sim_vars$variable == var, "value"] <-
          paste(as.character(input_list[[var]]), collapse = " ")
      }
    }
    
    # Create one R variable per simulation variable in this environment.
    # These variables are character strings; lists of values (population
    # sizes, local_adap, clinal_adap) are space delimited
    list2env(
      utils.wf.ref.values(
        sim_vars,
        char_vars = c("chromosome_name", "dispersal_type_phase1",
                      "dispersal_type_phase2", "natural_selection_model",
                      "population_size_phase1", "population_size_phase2",
                      "local_adap", "clinal_adap", "sib_mating_phase1",
                      "sib_mating_phase2", "real_freq_shrink",
                      "variance_offspring_phase1", "variance_offspring_phase2",
                      "ne_phase1", "ne_phase2", "inbreeding_founders",
                      "sample_families", "sample_crosses", "sample_design")),
      envir = environment())
    
    # -------------------------------
    # EXTRACT REFERENCE TABLE INFORMATION
    # -------------------------------
    reference <- ref_table$reference
    ref_vars <- ref_table$ref_vars
    
    # Identify positions of neutral and selected loci in the reference table.
    neutral_loci_location <- which(reference$type == "neutral" |
                                     reference$type == "real")
    adv_loci <- which(reference$type=="mutation_adv" | reference$type=="advantageous")
    mutation_loci_adv <- which(reference$type == "mutation_adv")
    mutation_loci_del <- which(reference$type == "mutation_del")
    mutation_loci_neu <- which(reference$type == "mutation_neu")
    
    # Combine mutation loci positions and sort them.
    mutation_loci_location <- c(mutation_loci_adv, mutation_loci_del, mutation_loci_neu)
    mutation_loci_location <- mutation_loci_location[order(mutation_loci_location)]
    mutation_loci_types <- mutation_loci_location
    
    # Identify loci with real data.
    real <- which(reference$type == "real")
    
    # Convert the neutral allele frequency to numeric.
    q_neutral <- as.numeric(ref_vars[ref_vars$variable=="q_neutral", "value"])
    
    # Check consistency between simulation variables and reference table variables.
    real_freq_table <- ref_vars[ref_vars$variable=="real_freq", "value"]
    if (real_freq_table != real_freq) {
      stop(error("  The value for the real_freq parameter was set differently",
                 "in the simulations and in the creation of the reference",
                 "table. They should be the same. Please check it.\n"))
    }
    
    real_loc_table <- ref_vars[ref_vars$variable=="real_loc", "value"]
    if (real_loc_table != real_loc) {
      stop(error("  The value for the real_loc parameter was set differently",
                 "in the simulations and in the creation of the reference",
                 "table. They should be the same. Please check it.\n"))
    }
    
    # Ensure that if real dataset values are required, the 'x' parameter is provided.
    if ((real_pops == TRUE | real_pop_size == TRUE | real_loc == TRUE |
         real_freq == TRUE | real_inbreeding == TRUE | real_migration == TRUE) &&
         is.null(x)) {
      stop(error("  The real dataset to extract information is missing\n"))
    }
    
    # If phase1 is disabled, set its generation count to zero.
    if (phase1 == FALSE) {
      gen_number_phase1 <- 0
    }
    
    # -------------------------------
    # CALCULATE TOTAL NUMBER OF GENERATIONS
    # -------------------------------
    number_generations <- gen_number_phase1 + gen_number_phase2
    
    # Define at which generations to store output genlight objects.
    gen_store <- c(seq(1, number_generations, every_gen), number_generations)
    ## Each iteration is a list of genlight objects named by the generation
    ## they hold
    final_res <- rep(list(list()), number_iterations)
    
    # -------------------------------
    # SET UP LOCI AND RECOMBINATION MAP
    # -------------------------------
    loci_number <- nrow(reference)
    recombination_map <- reference[, c("c", "loc_bp", "loc_cM")]
    
    # Adjust the recombination map so that the overall recombination probability
    # matches the number of recombination events per meiosis.
    recom_event <- ceiling(sum(recombination_map[, "c"], na.rm = TRUE))
    recombination_map[loci_number + 1, 1] <- recom_event - sum(recombination_map[, 1])
    recombination_map[loci_number + 1, 2] <- recombination_map[loci_number, 2]
    recombination_map[loci_number + 1, 3] <- recombination_map[loci_number, 3]
    
    # Prepare a map for the plink format with chromosome, locus, cM and bp positions.
    plink_map <- as.data.frame(matrix(nrow = nrow(reference), ncol = 4))
    plink_map[, 1] <- reference$chr_name
    plink_map[, 2] <- rownames(reference)
    plink_map[, 3] <- reference$loc_cM
    plink_map[, 4] <- reference$loc_bp
    
    # -------------------------------
    # SPLIT LISTS OF VALUES
    # -------------------------------
    # Space-delimited strings become numeric vectors. An empty local_adap or
    # clinal_adap stays NULL: numeric(0) is not NULL and used to switch on
    # local adaptation in no population, silencing advantageous selection
    split_num <- function(v) {
      v <- as.numeric(unlist(strsplit(trimws(v), " +")))
      if (length(v) == 0) NULL else v
    }
    population_size_phase2 <- split_num(population_size_phase2)
    population_size_phase1 <- split_num(population_size_phase1)
    local_adap <- split_num(local_adap)
    clinal_adap <- split_num(clinal_adap)
    sib_mating_phase1 <- split_num(sib_mating_phase1)
    sib_mating_phase2 <- split_num(sib_mating_phase2)
    variance_offspring_phase1 <- split_num(variance_offspring_phase1)
    variance_offspring_phase2 <- split_num(variance_offspring_phase2)
    ne_phase1 <- split_num(ne_phase1)
    ne_phase2 <- split_num(ne_phase2)
    inbreeding_founders <- split_num(inbreeding_founders)

    # -------------------------------
    # DETERMINE NUMBER OF POPULATIONS
    # -------------------------------
    if (real_pops == TRUE) {
      number_pops_phase1 <- number_pops_phase2 <- nPop(x)
    }
    number_pops <- if (phase1 == TRUE) number_pops_phase1 else number_pops_phase2
    
    ## Real census sizes (rounded up to even numbers) replace the sizes of
    ## the first phase simulated
    if (real_pop_size == TRUE) {
      real_sizes <- unname(unlist(table(pop(x))))
      real_sizes <- (real_sizes %% 2 != 0) + real_sizes
      if (phase1 == TRUE) {
        population_size_phase1 <- real_sizes
      } else {
        population_size_phase2 <- real_sizes
      }
    }

    # -------------------------------
    # EFFECTIVE POPULATION SIZE
    # -------------------------------
    # variance_offspring is one value or one per population. With ne_phase*
    # (one value or one per population) it is set so the expected Ne equals
    # the target (see ne_ratio()); Ne above its maximum for the census size
    # (N, or N/2 with replace_parents = TRUE) is capped there with a warning.
    # ne_expected_phase* is the expected Ne of each population
    resolve_ne <- function(k, ne, sizes, phase) {
      sizes <- rep_len(sizes, number_pops)
      if (anyNA(k) || any(k <= 0) || !length(k) %in% c(1, number_pops)) {
        stop(error("  variance_offspring_", phase, " must be positive, one",
                   " value or one per population\n", sep = ""))
      }
      k <- rep_len(k, number_pops)
      if (!is.null(ne)) {
        if (anyNA(ne) || any(ne <= 0) || !length(ne) %in% c(1, number_pops)) {
          stop(error("  ne_", phase, " must be positive, one value or one",
                     " per population\n", sep = ""))
        }
        ne <- rep_len(ne, number_pops)
        k <- k_for_ne(ne / sizes, replace_parents)
        if (any(is.infinite(k)) && verbose >= 1) {
          max_ne <- sizes * ne_ratio(1e6, replace_parents)
          cat(warn("  Warning: ne_", phase, " is above the largest Ne the ",
                   "census size allows in population(s) ",
                   paste(which(is.infinite(k)), collapse = ", "), " (",
                   paste(round(max_ne[is.infinite(k)]), collapse = ", "),
                   if (replace_parents) "; N/2 with replace_parents = TRUE",
                   "); increase population_size_", phase, "\n", sep = ""))
        }
        k[is.infinite(k)] <- 1e6
      }
      list(k = k, ne = sizes * ne_ratio(k, replace_parents))
    }
    res_ne <- resolve_ne(variance_offspring_phase2, ne_phase2,
                         population_size_phase2, "phase2")
    variance_offspring_phase2 <- res_ne$k
    ne_expected_phase2 <- res_ne$ne
    ne_expected_phase1 <- NULL
    if (phase1 == TRUE) {
      res_ne <- resolve_ne(variance_offspring_phase1, ne_phase1,
                           population_size_phase1, "phase1")
      variance_offspring_phase1 <- res_ne$k
      ne_expected_phase1 <- res_ne$ne
    }
    if (verbose >= 2) {
      message(report("  Expected Ne of each population in phase 2:",
                     paste(round(ne_expected_phase2), collapse = " "), "\n"))
    }

    # Stored sample of each population: x's sample sizes with
    # real_sample_size = TRUE, instead of sample_percent
    sample_sizes_real <- NULL
    if (real_sample_size == TRUE) {
      if (nPop(x) != number_pops) {
        stop(error("  real_sample_size needs as many populations in x as",
                   "populations simulated; use real_pops = TRUE\n"))
      }
      sample_sizes_real <- unname(unlist(table(pop(x))))
      smallest <- pmin(rep_len(population_size_phase2, number_pops),
                       if (phase1 == TRUE)
                         rep_len(population_size_phase1, number_pops) else Inf)
      if (any(sample_sizes_real > smallest)) {
        stop(error("  real_sample_size: the sample of x is larger than the",
                   "population simulated in population(s)",
                   paste(which(sample_sizes_real > smallest), collapse = ", "),
                   "\n"))
      }
    }
    
    # Family-structured samples: sizes of the full-sib families sampled in
    # each population, space delimited, with ";" between populations (one
    # group is used for every population; an empty group means no
    # families). Applies from generation 1 (founders have no parents)
    family_sizes <- NULL
    family_crosses <- NULL
    if ((!is.null(sample_crosses) || !is.null(sample_design)) &&
        is.null(sample_families)) {
        stop(error("  sample_crosses and sample_design need sample_families",
                   "(the size of each family)\n"))
    }
    if (!is.null(sample_families)) {
      groups <- trimws(strsplit(sample_families, ";")[[1]])
      family_sizes <- lapply(groups, function(g) {
        suppressWarnings(as.numeric(unlist(strsplit(g, " +"))))
      })
      if (length(family_sizes) == 1) {
        family_sizes <- rep(family_sizes, number_pops)
      }
      if (length(family_sizes) != number_pops ||
          any(vapply(family_sizes, function(v) {
            anyNA(v) || any(v < 2 | v %% 1 != 0)
          }, logical(1)))) {
        stop(error("  sample_families must give whole family sizes of at",
                   "least 2, one group for all populations or one group per",
                   "population separated by \";\"\n"))
      }
      
      # Planted crosses: one "SxD" token per family (the same label is the
      # same individual), ";" between populations; sample_design "axb" is
      # shorthand for every cross of a sires and b dams
      if (!is.null(sample_crosses) && !is.null(sample_design)) {
        stop(error("  Use sample_crosses or sample_design, not both\n"))
      }
      if (!is.null(sample_design)) {
        sample_crosses <- paste(vapply(
          trimws(strsplit(sample_design, ";")[[1]]),
          function(d) paste(expand_design(d), collapse = " "), character(1)),
          collapse = ";")
      }
      if (!is.null(sample_crosses)) {
        family_crosses <- lapply(trimws(strsplit(sample_crosses, ";")[[1]]),
                                 function(g) unlist(strsplit(g, " +")))
        if (length(family_crosses) == 1) {
          family_crosses <- rep(family_crosses, number_pops)
        }
        valid <- length(family_crosses) == number_pops &&
          all(vapply(seq_len(number_pops), function(i) {
            cr <- family_crosses[[i]]
            length(cr) == length(family_sizes[[i]]) &&
              all(grepl("^[^x]+x[^x]+$", cr))
          }, logical(1)))
        if (!valid) {
          stop(error("  sample_crosses (or sample_design) must give one",
                     "\"SxD\" cross per family of sample_families, e.g.",
                     "\"S1xD1 S1xD2 S2xD1 S2xD2\" or sample_design = \"2x2\"",
                     "for four families\n"))
        }
      }
      if (!is.null(sample_sizes_real)) {
        n_parents <- if (!is.null(family_crosses)) {
          vapply(family_crosses, function(cr) {
            length(unique(unlist(strsplit(cr, "x"))))
          }, numeric(1))
        } else 2 * vapply(family_sizes, length, numeric(1))
        need <- vapply(family_sizes, sum, numeric(1)) +
          if (sample_parents == TRUE) n_parents else 0
        if (any(need > sample_sizes_real)) {
          stop(error("  sample_families: families",
                     if (sample_parents == TRUE) "and their parents",
                     "exceed the sample size of x in population(s)",
                     paste(which(need > sample_sizes_real), collapse = ", "),
                     "\n"))
        }
      }
      max_size <- max(unlist(family_sizes), 0)
      n_off <- min(number_offspring_phase2,
                   if (phase1 == TRUE) number_offspring_phase1 else Inf)
      if (is.null(family_crosses) && max_size > n_off && verbose >= 1) {
        cat(warn("  Warning: families of", max_size, "are sampled but the",
                 "mean number of offspring per mating is", n_off, "; large",
                 "families will be rare. Increase number_offspring\n"))
      }
    }

    # -------------------------------
    # EXTRACT FREQUENCY INFORMATION FROM THE REAL DATA (IF APPLICABLE)
    # -------------------------------
    # The simulated allele "1" is stored as genotype 2, so it takes the
    # frequency of the alternative allele (alf2). Loci without calls in a
    # population take the frequency across all populations, then q_neutral.
    # With real_loc = TRUE the reference table orders real loci by position,
    # so frequencies are ordered by position too
    pop_list_freq <- rep(NA, number_pops)
    freq_shrink_lambda <- NULL
    if (real_freq == TRUE) {
      x_freq <- x
      if (real_loc == TRUE) {
        loc_to_keep <- which(as.character(x@chromosome) == chromosome_name)
        loc_to_keep <- loc_to_keep[order(x@position[loc_to_keep])]
        x_freq <- x[, loc_to_keep]
        x_freq@other$loc.metrics <- x@other$loc.metrics[loc_to_keep, , drop = FALSE]
      }
      # Alternative-allele frequency (as gl.alf()$alf2), computed here so
      # it does not depend on the gl.alf() signature of the installed
      # dartR.base (the CRAN release has no verbose argument)
      alt_freq <- function(g) colMeans(as.matrix(g), na.rm = TRUE) / 2
      freq_pooled <- alt_freq(x_freq)
      freq_pooled[is.na(freq_pooled)] <- q_neutral
      pop_list_freq <- lapply(seppop(x_freq), function(p) {
        f <- alt_freq(p)
        f[is.na(f)] <- freq_pooled[is.na(f)]
        return(f)
      })

      # Differentiation among founders: frequencies are shrunk toward
      # their mean so that the sampling noise of x does not inflate FST
      # (see shrink_freq())
      if (!is.null(real_freq_shrink) && length(pop_list_freq) > 1) {
        shrink <- if (real_freq_shrink == "auto") "auto" else
          suppressWarnings(as.numeric(real_freq_shrink))
        if (!identical(shrink, "auto") &&
            (is.na(shrink) || shrink < 0 || shrink > 1)) {
          stop(error("  real_freq_shrink must be NULL, \"auto\" or a number",
                     "between 0 and 1\n"))
        }
        n_genotyped <- lapply(seppop(x_freq), function(p) {
          colSums(!is.na(as.matrix(p)))
        })
        shrunk <- shrink_freq(pop_list_freq, n_genotyped, shrink)
        pop_list_freq <- shrunk$freq
        freq_shrink_lambda <- shrunk$lambda
        if (verbose >= 2) {
          message(report("  Population frequencies shrunk toward their mean,",
                         "lambda =", round(freq_shrink_lambda, 3), "\n"))
        }
      }
    }

    # -------------------------------
    # INBREEDING FROM THE REAL DATA AND PROPORTION OF SIB MATINGS
    # -------------------------------
    # F = 1 - Ho / He of each real population, on the loci of chromosome_name
    # when real_loc = TRUE. It makes the founders inbred and, in a phase
    # whose sib_mating is not set, it sets the proportion of full-sib
    # matings that keeps F at equilibrium, F = b / (4 - 3b), so
    # b = 4F / (1 + 3F)
    F_real <- NULL
    sib_real <- NULL
    # inbreeding_founders (one value or one per population) replaces the F
    # estimated from x, e.g. F estimated on high call-rate loci, as
    # heterozygote dropout inflates F in the full locus set
    if (real_inbreeding == TRUE && is.null(inbreeding_founders)) {
      if (nPop(x) != number_pops) {
        stop(error("  real_inbreeding needs as many populations in x as",
                   "populations simulated in the first phase (", nPop(x),
                   "in x,", number_pops, "simulated); use real_pops = TRUE\n"))
      }
      x_inb <- x
      if (real_loc == TRUE) {
        x_inb <- x[, which(as.character(x@chromosome) == chromosome_name)]
      }
      F_real <- inbreeding_real(x_inb)
      if (anyNA(F_real) && verbose >= 1) {
        cat(warn("  Warning: F cannot be estimated in population(s)",
                 paste(names(F_real)[is.na(F_real)], collapse = ", "),
                 "(no polymorphic loci); it is set to 0\n"))
      }
      F_real[is.na(F_real)] <- 0
      if (any(F_real < 0) && verbose >= 1) {
        cat(warn("  Warning: F is negative in population(s)",
                 paste(names(F_real)[F_real < 0], collapse = ", "),
                 "(excess of heterozygotes); it is set to 0\n"))
      }
      F_real <- pmax(F_real, 0)
    }
    if (!is.null(inbreeding_founders)) {
      if (anyNA(inbreeding_founders) || any(inbreeding_founders > 1) ||
          !length(inbreeding_founders) %in% c(1, number_pops)) {
        stop(error("  inbreeding_founders must be values up to 1, one value",
                   "or one per population\n"))
      }
      if (real_inbreeding == TRUE && verbose >= 1) {
        cat(report("  inbreeding_founders is set, so F is not estimated",
                   "from x\n"))
      }
      if (any(inbreeding_founders < 0) && verbose >= 1) {
        cat(warn("  Warning: negative values of inbreeding_founders are",
                 "set to 0\n"))
      }
      F_real <- pmax(rep_len(inbreeding_founders, number_pops), 0)
      names(F_real) <- if (real_pops == TRUE) popNames(x) else
        as.character(seq_len(number_pops))
    }
    if (!is.null(F_real)) {
      sib_real <- unname(4 * F_real / (1 + 3 * F_real))
      if (verbose >= 2) {
        message(report("  Inbreeding of the founders (F) and proportion of sib",
                       "matings:", paste0(names(F_real), " F = ",
                                          round(F_real, 3), ", sib = ",
                                          round(sib_real, 3),
                                          collapse = "; "), "\n"))
      }
    }

    # A phase's sib_mating is one value or one per population. When it is
    # not set it is the value from the founders' F (real_inbreeding or
    # inbreeding_founders) or 0 (random mating); a value that is set wins over the real data
    resolve_sib <- function(v, n_pops, phase) {
      if (is.null(v)) {
        return(if (is.null(sib_real)) rep(0, n_pops) else sib_real)
      }
      if (anyNA(v) || any(v < 0 | v > 1) || !length(v) %in% c(1, n_pops)) {
        stop(error("  sib_mating_", phase, " must be proportions between 0 ",
                   "and 1, one value or one per population\n", sep = ""))
      }
      if (!is.null(sib_real) && verbose >= 1) {
        cat(report("  sib_mating_", phase, " is set, so the proportion of sib ",
                   "matings set from the founders' F is not used in ", phase, "\n",
                   sep = ""))
      }
      return(rep_len(v, n_pops))
    }
    sib_mating_phase2 <- resolve_sib(sib_mating_phase2, number_pops,
                                     "phase2")
    if (phase1 == TRUE) {
      sib_mating_phase1 <- resolve_sib(sib_mating_phase1, number_pops,
                                       "phase1")
    }

    # -------------------------------
    # MIGRATION FROM THE REAL DATA
    # -------------------------------
    # Each connected pair swaps number_transfers individuals in each
    # direction, so with n populations "all_connected" is an island model
    # with Nm = (n - 1) * number_transfers immigrants per generation and
    # equilibrium FST = 1 / (1 + 4 Nm n / (n - 1)) = 1 / (1 + 4 n T). The
    # number of transfers per generation T that keeps x's FST (Hudson's) is
    # then T = (1 / FST - 1) / (4 n), whatever the population size. It can
    # be fractional: each pair moves floor(T) or floor(T) + 1 individuals,
    # the latter with probability T - floor(T)
    migrants_real <- NULL
    if (real_migration == TRUE) {
      if (nPop(x) < 2) {
        stop(error("  real_migration needs at least two populations in x\n"))
      }
      if (!is.null(file_dispersal)) {
        stop(error("  real_migration cannot be used with a dispersal file\n"))
      }
      disp_phases <- c(dispersal_phase2, if (phase1 == TRUE) dispersal_phase1)
      disp_types <- c(dispersal_type_phase2,
                      if (phase1 == TRUE) dispersal_type_phase1)
      if (any(disp_phases & disp_types != "all_connected")) {
        stop(error("  real_migration follows the island model and needs",
                   "dispersal_type = 'all_connected'\n"))
      }
      if (!any(disp_phases) && verbose >= 1) {
        cat(warn("  Warning: real_migration has no effect because dispersal",
                 "is off in every phase\n"))
      }
      x_mig <- x
      if (real_loc == TRUE) {
        x_mig <- x[, which(as.character(x@chromosome) == chromosome_name)]
      }
      fst_real <- fst_hudson(x_mig)
      if (!is.finite(fst_real) || fst_real < 0.001) {
        if (verbose >= 1) {
          cat(warn("  Warning: FST in x is", round(fst_real, 4), "; it is",
                   "set to 0.001 to compute the number of migrants\n"))
        }
        fst_real <- 0.001
      }
      # FST is set by Ne m, not N m, so T is multiplied by N / Ne (mean over
      # populations) of each phase
      migrants_base <- (1 / fst_real - 1) / (4 * number_pops)
      scale_ne <- function(sizes, ne) mean(rep_len(sizes, number_pops) / ne)
      migrants_real_phase2 <- migrants_base *
        scale_ne(population_size_phase2, ne_expected_phase2)
      if (phase1 == TRUE) {
        migrants_real_phase1 <- migrants_base *
          scale_ne(population_size_phase1, ne_expected_phase1)
      }
      if (verbose >= 2) {
        message(report("  FST in x =", round(fst_real, 4), "; individuals",
                       "transferred per pair of populations and generation",
                       "in phase 2 =", round(migrants_real_phase2, 3), "\n"))
      }
    }
    
    # -------------------------------
    # CALCULATE MUTATION DENSITY
    # -------------------------------
    # This calculation is based on the average proportion of heterozygotes per locus,
    # the number of loci, and the recombination map's length.
    freq_deleterious <- reference[-as.numeric(neutral_loci_location),]
    freq_deleterious_b <- mean(2 * (freq_deleterious$q) * (1 - freq_deleterious$q))
    density_mutations_per_cm <- (freq_deleterious_b * nrow(freq_deleterious)) /
      (recombination_map[loci_number, "loc_cM"] * 100)
    
    # -------------------------------
    # STORE A GENERATION AS A GENLIGHT OBJECT
    # -------------------------------
    # Samples sample_percent of each population (rounded up to an even
    # number, half of each sex) and stores it with the simulation variables.
    # disp_pairs is the dispersal table of the generation (NULL without
    # dispersal). Variables such as population_size are those of the
    # current phase when the function is called
    last_families <- NULL
    store_generation <- function(p_list, generation, iteration, disp_pairs,
                                 pool = NULL, parents = NULL) {
      families <- NULL
      if (!is.null(family_sizes) && generation >= 1 && !is.null(pool)) {
        # Family-structured sample: n = x's sample size, sample_percent or
        # the whole population
        n_sample <- if (!is.null(sample_sizes_real)) sample_sizes_real else
          if (sample_percent < 100) {
            n_s <- round(population_size * (sample_percent / 100))
            (n_s %% 2 != 0) + n_s
          } else population_size
        families <- lapply(pops_vector, function(x) {
          if (!is.null(family_crosses)) {
            sample_planted(p_list[[x]], parents[[x]], n_sample[x],
                           family_sizes[[x]], family_crosses[[x]],
                           sample_parents, generation, x, recom_event,
                           recombination, recombination_males,
                           recombination_map, loci_number)
          } else {
            sample_with_families(p_list[[x]], pool[[x]], parents[[x]],
                                 n_sample[x], family_sizes[[x]],
                                 sample_parents)
          }
        })
        # Kept for the pedigree of sampled offspring (see the generation loop)
        last_families <<- families
        if (any(vapply(families, function(f) f$short, logical(1))) &&
            verbose >= 1) {
          cat(warn("  Warning: generation", generation, "has fewer or smaller",
                   "families than sample_families asks for; they are",
                   "completed with random individuals\n"))
        }
        pop_list_temp <- lapply(families, function(f) f$sample)
        population_size_temp <- vapply(pop_list_temp, nrow, numeric(1))
      } else if (!is.null(sample_sizes_real)) {
        population_size_temp <- sample_sizes_real
        pop_list_temp <- lapply(pops_vector, function(x) {
          sample_sexes(p_list[[x]], population_size_temp[x])
        })
      } else if (sample_percent < 100) {
        population_size_temp <- round(population_size * (sample_percent / 100))
        population_size_temp <- (population_size_temp %% 2 != 0) + population_size_temp
        pop_list_temp <- lapply(pops_vector, function(x) {
          rbind(
            p_list[[x]][sample(which(p_list[[x]]$V1 == "Male"), size = population_size_temp[x] / 2),],
            p_list[[x]][sample(which(p_list[[x]]$V1 == "Female"), size = population_size_temp[x] / 2),]
          )
        })
      } else {
        population_size_temp <- population_size
        pop_list_temp <- p_list
      }
      
      # Combine and format simulation and reference variables for storage.
      s_vars_temp <- rbind(ref_vars, sim_vars)
      s_vars_temp <- setNames(data.frame(t(s_vars_temp[,-1])), s_vars_temp[, 1])
      s_vars_temp$generation <- generation
      s_vars_temp$iteration <- iteration
      s_vars_temp$seed <- seed
      s_vars_temp$del_ind_cM <- density_mutations_per_cm
      s_vars_temp$sample_percent <- sample_percent
      s_vars_temp$file_dispersal <- file_dispersal
      s_vars_temp$freq_shrink_lambda <- freq_shrink_lambda
      
      if (!is.null(disp_pairs)) {
        s_vars_temp$number_transfers_phase2 <- paste(disp_pairs$number_transfers, collapse = " ")  
        s_vars_temp$transfer_each_gen_phase2 <- paste(disp_pairs$transfer_each_gen, collapse = " ") 
      }
      s_vars_temp$ne_expected <- paste(round(ne_expected, 2), collapse = " ")
      s_vars_temp$variance_offspring_used <- paste(variance_offspring,
                                                   collapse = " ")
      s_vars_temp$population_size_used <- paste(population_size,
                                               collapse = " ")
      # With real_migration, the mean number of transfers per event
      # (phase 2, also for generation 0, as gl.diagnostics.sim() reads the
      # first stored generation)
      if (!is.null(migrants_real)) {
        s_vars_temp$migrants_real <- migrants_real
        s_vars_temp$number_transfers_phase2 <-
          migrants_real_phase2 * transfer_each_gen_phase2
      }
      
      res <- store(
        p_vector = pops_vector,
        p_size = population_size_temp,
        p_list = pop_list_temp,
        n_loc_1 = loci_number,
        ref = reference,
        p_map = plink_map,
        s_vars = s_vars_temp,
        g = generation
      )
      
      # Role of each individual in a family-structured sample
      if (!is.null(families)) {
        ids <- unlist(lapply(families, function(f) f$sample$id))
        m <- match(indNames(res), ids)
        res@other$ind.metrics$sample_role <-
          unlist(lapply(families, function(f) f$role))[m]
        res@other$ind.metrics$family <-
          unlist(lapply(families, function(f) f$family))[m]
        if (!is.null(family_crosses)) {
          res@other$ind.metrics$sample_cross <-
            unlist(lapply(families, function(f) f$cross))[m]
          res@other$ind.metrics$parent_label <-
            unlist(lapply(families, function(f) f$label))[m]
        }
      }
      
      # Assign population names to the stored object.
      if (real_pops == TRUE) {
        popNames(res) <- popNames(x)
      } else {
        popNames(res) <- as.character(pops_vector)
      }
      return(res)
    }
    
    # The weak-selection warning is printed once per run
    warn_pool <- TRUE
    # So is the warning that too few females have a brother for sib_mating
    warn_sib <- TRUE
    
    # -------------------------------
    # START SIMULATION ITERATION LOOP
    # -------------------------------
    for (iteration in 1:number_iterations) {
      # Each iteration starts with the full pool of mutation loci, so
      # iterations are independent replicates
      mutation_loci_location <- mutation_loci_types
      
      if (iteration %% 1 == 0 & verbose >= 2) {
        message(report(" Starting iteration =", iteration, "\n"))
      }
      
      # -------------------------------
      # SETUP VARIABLES FOR PHASE 1 (IF APPLICABLE)
      # -------------------------------
      if (phase1 == TRUE) {
        population_size <- population_size_phase1
        
        # Assign phase 1 specific simulation parameters.
        selection <- selection_phase1
        dispersal <- dispersal_phase1
        dispersal_type <- dispersal_type_phase1
        number_transfers <- number_transfers_phase1
        transfer_each_gen <- transfer_each_gen_phase1
        variance_offspring <- variance_offspring_phase1
        ne_expected <- ne_expected_phase1
        if (real_migration == TRUE) migrants_real <- migrants_real_phase1
        number_offspring <- number_offspring_phase1
        sib_mating <- sib_mating_phase1
        
        store_values <- store_phase1
        

        # Sex of the next single-individual transfer of each dispersal pair,
        # reset at the start of each phase (see dispersal_event())
        next_male <- NULL
        
      } else {
        # -------------------------------
        # SETUP VARIABLES FOR PHASE 2
        # -------------------------------
        population_size <- population_size_phase2
        variance_offspring <- variance_offspring_phase2
        ne_expected <- ne_expected_phase2
        if (real_migration == TRUE) migrants_real <- migrants_real_phase2
      }
      
      # -------------------------------
      # ERROR CHECKS ON POPULATION NUMBERS BETWEEN PHASES
      # -------------------------------
      if (phase1 == TRUE & number_pops_phase1 != number_pops_phase2) {
        stop(error("  Number of populations in phase 1 and phase 2 must be the",
                   "same\n"))
      }
      
      if (length(population_size_phase2) != number_pops_phase2) {
        stop(error("  Number of entries for population sizes do not agree with",
                   "the number of populations for phase 2\n"))
      }
      
      if (length(population_size_phase1) != number_pops_phase1 & phase1 == TRUE) {
        stop(error("  Number of entries for population sizes do not agree with",
                   "the number of populations for phase 1\n"))
      }
      
      # -------------------------------
      # INITIALIZE POPULATIONS
      # -------------------------------
      if (verbose >= 2) {
        message(report("  Initialising populations\n"))
      }
      
      # Generate chromosomes for all individuals (each individual has two chromosomes).
      chr_temp <- utils.wf.cpp()$make_chr(j = sum(population_size) * 2, q = reference$q)
      # Split chromosomes among populations.
      chr_pops_temps <- split(chr_temp, rep(1:number_pops, (c(population_size) * 2)))
      chr_pops <- lapply(chr_pops_temps, split, c(1:2))
      
      pops_vector <- 1:number_pops
      pop_list <- as.list(pops_vector)
      # Proportion of loci (founder_F) and of the map (founder_F_map)
      # identical by descent in each founder, by id. For the map, each locus
      # stands for the map between the midpoints to its neighbours
      founder_F <- founder_F_map <- c()
      mid <- (head(reference$loc_cM, -1) + reference$loc_cM[-1]) / 2
      locus_map <- diff(c(min(reference$loc_cM), mid, max(reference$loc_cM)))
      
      # Create the population data frames with sex, population number and chromosomes.
      for (pop_n in pops_vector) {
        pop <- as.data.frame(matrix(ncol = 4, nrow = population_size[pop_n]))
        pop[, 1] <- rep(c("Male", "Female"), each = population_size[pop_n] / 2)
        pop[, 2] <- pop_n   # Population identifier
        pop[, 3] <- chr_pops[[pop_n]][1]  # First chromosome
        pop[, 4] <- chr_pops[[pop_n]][2]  # Second chromosome
        pop$id <- paste0("0_",pop_n,"_",1:nrow(pop)) # ID
        
        # Real loci take the real allele frequencies of this population
        if (real_freq == TRUE) {
          q_real <- pop_list_freq[[pop_n]]
          for (individual_pop in 1:population_size[pop_n]) {
            stringi::stri_sub_all(pop[individual_pop, 3], from = real, length = 1) <-
              as.character(as.integer(runif(length(q_real)) < q_real))
            stringi::stri_sub_all(pop[individual_pop, 4], from = real, length = 1) <-
              as.character(as.integer(runif(length(q_real)) < q_real))
          }
        }

        # Inbred founders: in the stretches of the map that are identical by
        # descent, the second chromosome copies the first
        if (!is.null(F_real) && F_real[pop_n] > 0) {
          for (individual_pop in 1:population_size[pop_n]) {
            ibd <- which(ibd_loci(reference$loc_cM, F_real[pop_n]))
            founder_F[pop[individual_pop, "id"]] <- length(ibd) / loci_number
            founder_F_map[pop[individual_pop, "id"]] <-
              if (sum(locus_map) > 0) sum(locus_map[ibd]) / sum(locus_map) else
                length(ibd) / loci_number
            if (length(ibd) > 0) {
              chr_2 <- strsplit(pop[individual_pop, 4], "")[[1]]
              chr_2[ibd] <- strsplit(pop[individual_pop, 3], "")[[1]][ibd]
              pop[individual_pop, 4] <- paste(chr_2, collapse = "")
            }
          }
        }

        # Save the initialized population in the list.
        pop_list[[pop_n]] <- pop
      }
      
      # If only one population exists, disable dispersal.
      if (length(pop_list) == 1) {
        dispersal <- FALSE
      }
      
      # Pedigree of every individual of every generation (store_pedigree)
      pedigree <- list()
      record_pedigree <- function(p_list, generation) {
        p <- rbindlist(p_list, fill = TRUE)
        pop_labels <- if (real_pops == TRUE) popNames(x) else
          as.character(pops_vector)
        data.frame(
          id = p$id,
          pat = if (is.null(p$V5)) NA_character_ else as.character(p$V5),
          mat = if (is.null(p$V6)) NA_character_ else as.character(p$V6),
          generation = generation,
          pop = pop_labels[as.numeric(p$V2)],
          F_founder = if (generation == 0) {
            v <- unname(founder_F[p$id])
            if (is.null(v)) v <- rep(0, nrow(p))
            v[is.na(v)] <- 0
            v
          } else NA_real_,
          in_population = TRUE,
          stringsAsFactors = FALSE
        )
      }
      if (store_pedigree == TRUE) {
        pedigree[["0"]] <- record_pedigree(pop_list, 0)
      }
      
      # Founders, before they reproduce, stored as generation 0. They have
      # no parents; F_founder is their proportion of loci identical by descent
      if (store_founders == TRUE) {
        founders <- lapply(pop_list, function(p) {
          p$V5 <- NA
          p$V6 <- NA
          p
        })
        gen0 <- store_generation(founders, 0, iteration, NULL)
        by_id <- function(v) {
          out <- rep(0, nInd(gen0))
          if (length(v) > 0) {
            out <- unname(v[indNames(gen0)])
            out[is.na(out)] <- 0
          }
          out
        }
        gen0@other$ind.metrics$F_founder <- by_id(founder_F)
        gen0@other$ind.metrics$F_founder_map <- by_id(founder_F_map)
        final_res[[iteration]][["generation_0"]] <- gen0
      }
      
      # -------------------------------
      # START GENERATION LOOP
      # -------------------------------
      for (generation in 1:number_generations) {
        if (phase1 == TRUE & generation == 1 & verbose >= 2) {
          message(report(" Starting phase 1\n"))
        }
        if (generation %% 10 == 0 & verbose >= 2) {
          message(report("  Starting generation =", generation, "\n"))
        }
        
        # -------------------------------
        # SWITCH FROM PHASE 1 TO PHASE 2 (IF APPLICABLE)
        # -------------------------------
        if (generation == (gen_number_phase1 + 1)) {
          if (verbose >= 2) {
            message(report(" Starting phase 2\n"))
          }
          
          # Update simulation parameters for phase 2.
          selection <- selection_phase2
          dispersal <- dispersal_phase2
          dispersal_type <- dispersal_type_phase2
          number_transfers <- number_transfers_phase2
          transfer_each_gen <- transfer_each_gen_phase2
          variance_offspring <- variance_offspring_phase2
          ne_expected <- ne_expected_phase2
          if (real_migration == TRUE) migrants_real <- migrants_real_phase2
          number_offspring <- number_offspring_phase2
          dispersal <- dispersal_phase2
          sib_mating <- sib_mating_phase2
          population_size <- population_size_phase2
          
          store_values <- TRUE
          
          # Reset the sex of the next single-individual transfer of each pair
          next_male <- NULL
          

          # Phase-2 founders are sampled from phase 1: from one phase-1
          # population for all (same_line = TRUE) or each from its own.
          # Sampling is without replacement unless phase 2 needs more
          # individuals of a sex than phase 1 has, so founders are not
          # cloned
          if (phase1 == TRUE) {
            sample_founders <- function(source, size) {
              males <- which(source$V1 == "Male")
              females <- which(source$V1 == "Female")
              rbind(
                source[males[sample.int(length(males), size = size / 2,
                                        replace = size / 2 > length(males))], ],
                source[females[sample.int(length(females), size = size / 2,
                                          replace = size / 2 > length(females))], ]
              )
            }
            pop_sample <- if (same_line == TRUE) sample(pops_vector, 1) else NULL
            pop_list <- lapply(pops_vector, function(x) {
              source_pop <- if (same_line == TRUE) pop_sample else x
              pop_temp <- sample_founders(pop_list[[source_pop]], population_size[x])
              pop_temp$V2 <- x
              return(pop_temp)
            })
          }
        }
        
        # -------------------------------
        # DISPERSAL PHASE
        # -------------------------------
        if (number_pops == 1) {
          dispersal <- FALSE
        }
        if (dispersal == TRUE) {
          if (is.null(file_dispersal)) {
            # Determine dispersal pairs based on dispersal type.
            if (dispersal_type == "all_connected") {
              dispersal_pairs <- as.data.frame(expand.grid(pops_vector, pops_vector))
              dispersal_pairs$same_pop <- dispersal_pairs$Var1 == dispersal_pairs$Var2
              dispersal_pairs <- dispersal_pairs[which(dispersal_pairs$same_pop == FALSE),]
              colnames(dispersal_pairs) <- c("pop1", "pop2", "same_pop")
            }
            if (dispersal_type == "line") {
              dispersal_pairs <- as.data.frame(rbind(
                cbind(head(pops_vector, -1), pops_vector[-1]),
                cbind(pops_vector[-1], head(pops_vector, -1))
              ))
              colnames(dispersal_pairs) <- c("pop1", "pop2")
            }
            if (dispersal_type == "circle") {
              dispersal_pairs <- as.data.frame(rbind(
                cbind(pops_vector, c(pops_vector[-1], pops_vector[1])),
                cbind(c(pops_vector[-1], pops_vector[1]), pops_vector)
              ))
              colnames(dispersal_pairs) <- c("pop1", "pop2")
            }
            
            # migration() swaps individuals in both directions, so each
            # connected pair is processed once
            pair_key <- paste(pmin(dispersal_pairs$pop1, dispersal_pairs$pop2),
                              pmax(dispersal_pairs$pop1, dispersal_pairs$pop2))
            dispersal_pairs <- dispersal_pairs[!duplicated(pair_key), ]
            
            # Set additional dispersal parameters.
            dispersal_pairs$number_transfers <- number_transfers
            dispersal_pairs$transfer_each_gen <- transfer_each_gen
            
            # real_migration: transfers per event (T per generation times the
            # generations between events), fractional part drawn per pair
            if (!is.null(migrants_real)) {
              t_event <- migrants_real * transfer_each_gen
              dispersal_pairs$number_transfers <- floor(t_event) +
                rbinom(nrow(dispersal_pairs), 1, t_event - floor(t_event))
            }
            
          } else {
            # If a dispersal file is provided, read the file.
            dispersal_pairs <- suppressWarnings(read.csv(file_dispersal))
          }
          
          # Population sizes for each pair (even-adjusted, as simulated)
          dispersal_pairs$size_pop1 <- population_size[dispersal_pairs$pop1]
          dispersal_pairs$size_pop2 <- population_size[dispersal_pairs$pop2]
          
          # Process each dispersal pair; which sexes move is decided per row
          # from its number_transfers
          res <- dispersal_event(pop_list, dispersal_pairs, generation,
                                 next_male)
          pop_list <- res[[1]]
          next_male <- res[[2]]
        }
        
        # Parents of this generation's offspring (for sample_parents)
        parents_list <- pop_list
        
        # -------------------------------
        # REPRODUCTION PHASE
        # -------------------------------
        offspring_list <- lapply(pops_vector, function(x) {
          tmp_rep <- reproduction(
            pop = pop_list[[x]],
            pop_number = x,
            pop_size = population_size[x],
            var_off = variance_offspring[x],
            num_off = number_offspring,
            r_event = recom_event,
            recom = recombination,
            r_males = recombination_males,
            r_map_1 = recombination_map,
            n_loc = loci_number,
            gen = generation,
            rep_parents = replace_parents,
            sib = sib_mating[x]
          )
          if (!is.null(tmp_rep)) {
            tmp_rep$id <- paste0(generation, "_", x, "_", 1:nrow(tmp_rep))
          }
          return(tmp_rep)
        })

        sib_short <- vapply(offspring_list, function(o) {
          isTRUE(attr(o, "sib_short"))
        }, logical(1))
        if (any(sib_short) & warn_sib == TRUE & verbose >= 1) {
          cat(warn("  Warning: in generation", generation, "fewer females",
                   "had a brother than the proportion of sib matings",
                   "(sib_mating), so fewer sib matings took place.",
                   "Increase number_offspring to have larger families.\n"))
          warn_sib <- FALSE
        }

        # Relative selection samples parents without replacement in
        # proportion to fitness; when the offspring pool is not much larger
        # than N, most offspring are kept and selection is weakened
        if (selection == TRUE & natural_selection_model == "relative" &
            warn_pool == TRUE & verbose >= 1) {
          pool_ratio <- sapply(pops_vector, function(x) {
            NROW(offspring_list[[x]]) / population_size[x]
          })
          if (any(pool_ratio < 3)) {
            cat(warn("  Warning: the offspring pool is less than 3 times the",
                     "population size (minimum ratio", round(min(pool_ratio), 2),
                     "), which weakens relative selection. Increase",
                     "number_offspring for stronger selection.\n"))
            warn_pool <- FALSE
          }
        }
        
        # -------------------------------
        # MUTATION PHASE
        # -------------------------------
        if (mutation == TRUE) {
          for (off_pop in 1:length(offspring_list)) {
            offspring_pop <- offspring_list[[off_pop]]
            if (is.null(offspring_pop)) {
              next
            }
            # Add a uniform random number to each offspring for mutation decision.
            offspring_pop$runif <- runif(nrow(offspring_pop))
            
            for (offspring_ind in 1:nrow(offspring_pop)) {
              if (length(mutation_loci_location) == 0) {
                if (verbose >= 2) {
                  message(important("  No more locus to mutate\n"))
                }
                break()
              }
              
              # Determine if a mutation occurs based on the mutation rate.
              if (offspring_pop[offspring_ind, "runif"] < mut_rate) {
                # Sample by index: sample() on a single number n draws from 1:n
                locus_to_mutate <- mutation_loci_location[
                  sample.int(length(mutation_loci_location), 1)]
                # Remove the mutated locus from available mutation loci.
                mutation_loci_location <- mutation_loci_location[-which(mutation_loci_location == locus_to_mutate)]
                chromosomes <- c(offspring_pop[offspring_ind, 3], offspring_pop[offspring_ind, 4])
                chr_to_mutate <- sample(1:2, 1)
                chr_to_mutate_b <- chromosomes[chr_to_mutate]
                # Mutate the selected locus by setting it to "1".
                substr(chr_to_mutate_b, as.numeric(locus_to_mutate),
                       as.numeric(locus_to_mutate)) <- "1"
                chromosomes[chr_to_mutate] <- chr_to_mutate_b
                offspring_pop[offspring_ind, 3] <- chromosomes[1]
                offspring_pop[offspring_ind, 4] <- chromosomes[2]
              } else {
                next()
              }
            }
            offspring_list[[off_pop]] <- offspring_pop
          }
        }
        
        # -------------------------------
        # SELECTION PHASE
        # -------------------------------
        if (selection == TRUE) {
          if (!is.null(local_adap)) {
            pops_local <- setdiff(pops_vector, local_adap)
            # Create a copy of the reference parameters for each population.
            reference_local <- replicate(length(pops_vector), 
                                         reference[, c("s", "h")], 
                                         simplify = FALSE)
            
            # Set selection coefficient s to 0 for populations not under local adaptation.
            reference_local[pops_local] <- lapply(reference_local[pops_local], 
                                                  function(y) {
                                                    y[adv_loci, "s"] <- 0
                                                    return(y)
                                                  })
            
            # Apply selection using the modified reference for local adaptation.
            offspring_list <- lapply(pops_vector, function(x) {
              selection_fun(
                offspring = offspring_list[[x]],
                h = reference_local[[x]][, "h"],
                s = reference_local[[x]][, "s"],
                sel_model = natural_selection_model,
                g_load = genetic_load
              )
            })
            
          } else if (!is.null(clinal_adap)) {
            # Populations clinal_adap[1] to clinal_adap[2] form the cline.
            # Along it, advantageous s is multiplied by 1, 1 - k, 1 - 2k, ...
            # (k = clinal_strength / 100, floored at 0); populations outside
            # the cline get advantageous s = 0, as in local adaptation
            pops_clinal <- seq(clinal_adap[1], clinal_adap[2])
            clinal_s <- pmax(0, 1 - (seq_along(pops_clinal) - 1) *
                               (clinal_strength / 100))
            reference_clinal <- lapply(pops_vector, function(y) {
              ref_pop <- reference[, c("s", "h")]
              position_cline <- match(y, pops_clinal)
              ref_pop[adv_loci, "s"] <- if (is.na(position_cline)) {
                0
              } else {
                ref_pop[adv_loci, "s"] * clinal_s[position_cline]
              }
              return(ref_pop)
            })
            
            offspring_list <- lapply(pops_vector, function(x) {
              selection_fun(
                offspring = offspring_list[[x]],
                h = reference_clinal[[x]][, "h"],
                s = reference_clinal[[x]][, "s"],
                sel_model = natural_selection_model,
                g_load = genetic_load
              )
            })
            
          } else {
            # Apply selection using the original reference parameters.
            offspring_list <- lapply(pops_vector, function(x) {
              selection_fun(
                offspring = offspring_list[[x]],
                h = reference[, "h"],
                s = reference[, "s"],
                sel_model = natural_selection_model,
                g_load = genetic_load
              )
            })
          }
        }
        
        # -------------------------------
        # SAMPLING OF THE NEXT GENERATION
        # -------------------------------
        # A population goes extinct when its offspring (after selection)
        # include fewer males or fewer females than half its size. The
        # iteration stops; generations stored so far are kept
        test_extinction <- sapply(pops_vector, function(x) {
          sexes <- offspring_list[[x]]$V1
          sum(sexes == "Male") < population_size[x] / 2 |
            sum(sexes == "Female") < population_size[x] / 2
        })
        
        if (any(test_extinction)) {
          if (verbose >= 1) {
            cat(important("  Population(s)",
                          paste(pops_vector[test_extinction], collapse = ", "),
                          "became extinct at generation", generation,
                          "of iteration", iteration, ". Stopping this",
                          "iteration; generations stored so far are kept.\n"))
          }
          break()
        }
        
        # -------------------------------
        # SAMPLING PARENTS FOR NEXT GENERATION (WITHOUT SELECTION)
        # -------------------------------
        if (selection == FALSE) {
          pop_list <- lapply(pops_vector, function(x) {
            rbind(
              offspring_list[[x]][sample(which(offspring_list[[x]]$V1 == "Male"), size = population_size[x] / 2),],
              offspring_list[[x]][sample(which(offspring_list[[x]]$V1 == "Female"), size = population_size[x] / 2),]
            )
          })
        }
        
        # -------------------------------
        # SAMPLING PARENTS WITH ABSOLUTE OR RELATIVE SELECTION
        # -------------------------------
        if (selection == TRUE & natural_selection_model == "absolute") {
          pop_list <- lapply(pops_vector, function(x) {
            rbind(
              offspring_list[[x]][sample(which(offspring_list[[x]]$V1 == "Male"), size = population_size[x] / 2),],
              offspring_list[[x]][sample(which(offspring_list[[x]]$V1 == "Female"), size = population_size[x] / 2),]
            )
          })
        }
        if (selection == TRUE & natural_selection_model == "relative") {
          # Selection is performed in proportion to relative fitness for each sex.
          pop_list <- lapply(pops_vector, function(x) {
            males_pop <- offspring_list[[x]][which(offspring_list[[x]]$V1 == "Male"),]
            females_pop <- offspring_list[[x]][which(offspring_list[[x]]$V1 == "Female"),]
            rbind(
              males_pop[sample(row.names(males_pop), size = (population_size[x] / 2),
                               prob = (males_pop$relative_fitness * 2)), ],
              females_pop[sample(row.names(females_pop), size = (population_size[x] / 2),
                                 prob = (females_pop$relative_fitness * 2)), ]
            )
          })
        }
        
        # -------------------------------
        # RECYCLE MUTATIONS FROM ELIMINATED LOCI
        # -------------------------------
        if (mutation == TRUE) {
          # Combine all populations.
          pops_merge <- rbindlist(pop_list)
          pops_seqs <- c(pops_merge$V3, pops_merge$V4)
          

          # Get frequencies for each locus.
          freqs <- utils.wf.cpp()$make_freqs(pops_seqs)
          # Mutation loci whose new allele was lost return to the pool; other
          # loci (neutral, real, selected) never receive mutations
          mutation_eliminated <- intersect(which(freqs == 0),
                                           mutation_loci_types)
          mutation_loci_location <- union(mutation_loci_location, mutation_eliminated)
        }
        
        # -------------------------------
        # STORE GENERATION RESULTS INTO GENLIGHT OBJECTS
        # -------------------------------
        if (store_pedigree == TRUE) {
          pedigree[[as.character(generation)]] <-
            record_pedigree(pop_list, generation)
        }
        
        if (generation %in% gen_store & store_values == TRUE) {
          gen_name <- paste0("generation_", generation)
          disp_pairs <- if (dispersal == TRUE) dispersal_pairs else NULL
          final_res[[iteration]][[gen_name]] <-
            store_generation(pop_list, generation, iteration, disp_pairs,
                             pool = offspring_list, parents = parents_list)
          # Sampled offspring that did not join the population (they were
          # not drawn as parents of the next generation) join the pedigree
          sampled <- final_res[[iteration]][[gen_name]]
          if (store_pedigree == TRUE &&
              !is.null(sampled@other$ind.metrics$sample_role)) {
            im <- sampled@other$ind.metrics
            in_pool <- setdiff(indNames(sampled)[im$sample_role == "family"],
                               pedigree[[as.character(generation)]]$id)
            if (length(in_pool) > 0) {
              fam_all <- do.call(rbind, lapply(last_families,
                                                function(f) f$sample))
              extra <- record_pedigree(
                list(fam_all[fam_all$id %in% in_pool, ]), generation)
              extra$in_population <- FALSE
              pedigree[[as.character(generation)]] <-
                rbind(pedigree[[as.character(generation)]], extra)
            }
          }
        }
      }  # End generation loop
      if (store_pedigree == TRUE) {
        attr(final_res[[iteration]], "pedigree") <-
          do.call(rbind, unname(pedigree))
      }
    }  # End iteration loop
    
    # -------------------------------
    # FINALIZE RESULTS
    # -------------------------------
    # Name the elements of the final results list; each iteration's
    # elements are already named by the generation they hold
    names(final_res) <- paste0("iteration_", 1:number_iterations)
    
    # -------------------------------
    # FLAG THE END OF THE FUNCTION EXECUTION
    # -------------------------------
    if (verbose >= 1) {
      message(report("Completed:", funname, "\n"))
    }
    
    # Return the final simulation results invisibly.
    return(invisible(final_res))
  }
