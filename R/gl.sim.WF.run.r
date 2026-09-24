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
#' \item Reproduction. N/2 monogamous pairs are formed at random. The number
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
#' If a genlight object is used (real_pops, real_pop_size, real_loc or
#' real_freq), the simulated allele is the alternative allele of the
#' genlight (genotype 2). Loci without calls in a population take the
#' frequency across all populations. real_pop_size sets the sizes of the
#' first phase simulated, rounded up to even numbers.
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
    ## The Shiny app does not ask for replace_parents; parents are then
    ## sampled without replacement, as in sim_variables.csv
    if (!"replace_parents" %in% sim_vars$variable) {
      sim_vars <- rbind(sim_vars,
                        data.frame(variable = "replace_parents",
                                   value = "FALSE"))
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
                      "local_adap", "clinal_adap")),
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
         real_freq == TRUE) && is.null(x)) {
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
    # EXTRACT FREQUENCY INFORMATION FROM THE REAL DATA (IF APPLICABLE)
    # -------------------------------
    # The simulated allele "1" is stored as genotype 2, so it takes the
    # frequency of the alternative allele (alf2). Loci without calls in a
    # population take the frequency across all populations, then q_neutral.
    # With real_loc = TRUE the reference table orders real loci by position,
    # so frequencies are ordered by position too
    pop_list_freq <- rep(NA, number_pops)
    if (real_freq == TRUE) {
      x_freq <- x
      if (real_loc == TRUE) {
        loc_to_keep <- which(as.character(x@chromosome) == chromosome_name)
        loc_to_keep <- loc_to_keep[order(x@position[loc_to_keep])]
        x_freq <- x[, loc_to_keep]
        x_freq@other$loc.metrics <- x@other$loc.metrics[loc_to_keep, , drop = FALSE]
      }
      freq_pooled <- suppressMessages(gl.alf(x_freq, verbose = 0))$alf2
      freq_pooled[is.na(freq_pooled)] <- q_neutral
      pop_list_freq <- lapply(seppop(x_freq), function(p) {
        f <- suppressMessages(gl.alf(p, verbose = 0))$alf2
        f[is.na(f)] <- freq_pooled[is.na(f)]
        return(f)
      })
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
    
    # The weak-selection warning is printed once per run
    warn_pool <- TRUE
    
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
        number_offspring <- number_offspring_phase1
        
        store_values <- store_phase1
        

        # Sex of the next single-individual transfer of each dispersal pair,
        # reset at the start of each phase (see dispersal_event())
        next_male <- NULL
        
      } else {
        # -------------------------------
        # SETUP VARIABLES FOR PHASE 2
        # -------------------------------
        population_size <- population_size_phase2
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
        
        # Save the initialized population in the list.
        pop_list[[pop_n]] <- pop
      }
      
      # If only one population exists, disable dispersal.
      if (length(pop_list) == 1) {
        dispersal <- FALSE
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
          number_offspring <- number_offspring_phase2
          dispersal <- dispersal_phase2
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
        
        # -------------------------------
        # REPRODUCTION PHASE
        # -------------------------------
        offspring_list <- lapply(pops_vector, function(x) {
          tmp_rep <- reproduction(
            pop = pop_list[[x]],
            pop_number = x,
            pop_size = population_size[x],
            var_off = variance_offspring,
            num_off = number_offspring,
            r_event = recom_event,
            recom = recombination,
            r_males = recombination_males,
            r_map_1 = recombination_map,
            n_loc = loci_number,
            gen = generation,
            rep_parents = replace_parents
          )
          if (!is.null(tmp_rep)) {
            tmp_rep$id <- paste0(generation, "_", x, "_", 1:nrow(tmp_rep))
          }
          return(tmp_rep)
        })
        
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
        if (generation %in% gen_store & store_values == TRUE) {
          gen_name <- paste0("generation_", generation)
          
          # Subsample individuals if sample_percent is less than 100.
          if (sample_percent < 100) {
            population_size_temp <- round(population_size * (sample_percent / 100))
            population_size_temp <- (population_size_temp %% 2 != 0) + population_size_temp
            pop_list_temp <- lapply(pops_vector, function(x) {
              rbind(
                pop_list[[x]][sample(which(pop_list[[x]]$V1 == "Male"), size = population_size_temp[x] / 2),],
                pop_list[[x]][sample(which(pop_list[[x]]$V1 == "Female"), size = population_size_temp[x] / 2),]
              )
            })
          } else {
            population_size_temp <- population_size
            pop_list_temp <- pop_list
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
          
          if (dispersal == TRUE) {
            s_vars_temp$number_transfers_phase2 <- paste(dispersal_pairs$number_transfers, collapse = " ")  
            s_vars_temp$transfer_each_gen_phase2 <- paste(dispersal_pairs$transfer_each_gen, collapse = " ") 
          }
          
          # Store the generation output.
          final_res[[iteration]][[gen_name]] <- store(
            p_vector = pops_vector,
            p_size = population_size_temp,
            p_list = pop_list_temp,
            n_loc_1 = loci_number,
            ref = reference,
            p_map = plink_map,
            s_vars = s_vars_temp,
            g = generation
          )
          
          # Assign population names to the stored object.
          if (real_pops == TRUE) {
            popNames(final_res[[iteration]][[gen_name]]) <- popNames(x)
          } else {
            popNames(final_res[[iteration]][[gen_name]]) <- as.character(pops_vector)
          }
        }
      }  # End generation loop
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
