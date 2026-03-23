#' @name gl.sim.WF.run
#' @title Runs Wright-Fisher simulations
#' @description
#' This function simulates populations made up of diploid organisms that 
#' reproduce in non-overlapping generations. Each individual has a pair of 
#' homologous chromosomes that contains interspersed selected and neutral loci. 
#' For the initial generation, the genotype for each individual’s chromosomes is
#' randomly drawn from distributions at linkage equilibrium and in 
#' Hardy-Weinberg equilibrium. 
#' 
#' See documentation and tutorial for a complete description of the simulations.
#' These documents can be accessed at 
#' https://github.com/green-striped-gecko/dartR/wiki/Simulations-tutorial
#' 
#' Take into account that the simulations will take a little longer the
#' first time you use the function gl.sim.WF.run() because C++ functions must
#' be compiled.
#' @param file_var Path of the variables file 'sim_variables.csv' (see details) 
#' [required if interactive_vars = FALSE].
#' @param ref_table Reference table created by the function 
#' \code{\link{gl.sim.WF.table}} [required].
#' @param x Name of the genlight object containing the SNP data to extract
#' values for some simulation variables (see details) [default NULL].
#' @param file_dispersal Path of the file with the dispersal table created with
#' the function \code{\link{gl.sim.create_dispersal}} [default NULL]. 
#' @param number_iterations Number of iterations of the simulations [default 1].
#' @param every_gen Generation interval at which simulations should be stored in
#' a genlight object [default 10].
#' @param sample_percent Percentage of individuals, from the total population, 
#' to sample and save in the genlight object every generation [default 50].
#' @param store_phase1 Whether to store simulations of phase 1 in genlight
#' objects [default FALSE].
#' @param interactive_vars Run a shiny app to input interactively the values of
#' simulations variables [default TRUE].
#' @param seed Set the seed for the simulations [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @param ... Any variable and its value can be added separately within the 
#' function, will be changed over the input value supplied by the csv file. See 
#' tutorial. 
#' @return Returns genlight objects with simulated data.
#' @author Custodian: Luis Mijangos
#' @examples
#' ref_table <- gl.sim.WF.table(file_var=system.file("extdata", 
#' "ref_variables.csv", package = "dartR.sim"),interactive_vars = FALSE)
#' 
#' res_sim <- gl.sim.WF.run(file_var = system.file("extdata",
#'  "sim_variables.csv", package ="dartR.sim"),ref_table=ref_table,
#'  interactive_vars = FALSE)
#' @seealso \code{\link{gl.sim.WF.table}}
#' @family simulation functions
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
                     build = "Jody",
                     verbose = verbose)
    
    # -------------------------------
    # CHECK FOR REQUIRED PACKAGES
    # -------------------------------
    pkg <- "stringi"
    if (!(requireNamespace(pkg, quietly = TRUE))) {
      message(error(
        "Package",
        pkg,
        " needed for this function to work. Please install it.\n"
      ))
      return(-1)
    }
    
    replace_parents <- NULL
    
    # -------------------------------
    # RETRIEVE SIMULATION VARIABLES
    # -------------------------------
    # If interactive_vars is TRUE, launch the shiny app to retrieve variables interactively.
    if (interactive_vars == TRUE) {
      
      sim_vars <- interactive_sim_run()
      
      # Wrap specific variable values in quotes to treat them as strings.
      sim_vars[sim_vars$variable=="population_size_phase2" ,"value"] <-
        paste0("'", sim_vars[sim_vars$variable=="population_size_phase2" ,"value"], "'")
      sim_vars[sim_vars$variable=="population_size_phase1" ,"value"] <-
        paste0("'", sim_vars[sim_vars$variable=="population_size_phase1" ,"value"], "'")
      sim_vars[sim_vars$variable=="dispersal_type_phase2" ,"value"] <- 
        paste0("'", sim_vars[sim_vars$variable=="dispersal_type_phase2" ,"value"], "'")
      sim_vars[sim_vars$variable=="dispersal_type_phase1" ,"value"] <- 
        paste0("'", sim_vars[sim_vars$variable=="dispersal_type_phase1" ,"value"], "'")
      sim_vars[sim_vars$variable=="natural_selection_model" ,"value"] <- 
        paste0("'", sim_vars[sim_vars$variable=="natural_selection_model" ,"value"], "'")
      
      # Order the simulation variables alphabetically by variable name.
      sim_vars <- sim_vars[order(sim_vars$variable),]
      
      # Create assignment strings for each variable (e.g., var <- value).
      vars_assign <- unlist(unname(
        mapply(paste, sim_vars$variable, "<-", sim_vars$value, SIMPLIFY = F)
      ))
      
      # Evaluate the assignments to create the variables in the environment.
      eval(parse(text = vars_assign))
      
    } else {
      # If interactive_vars is FALSE, read variables from a CSV file.
      sim_vars <- suppressWarnings(read.csv(file_var))
      sim_vars <- sim_vars[, 2:3]  # Use only the variable names and values columns.
      
      sim_vars <- sim_vars[order(sim_vars$variable),]
      
      vars_assign <- unlist(unname(
        mapply(paste, sim_vars$variable, "<-", sim_vars$value, SIMPLIFY = F)
      ))
      eval(parse(text = vars_assign))
    }
    
    # -------------------------------
    # OVERRIDE VARIABLES WITH ADDITIONAL ARGUMENTS (if any)
    # -------------------------------
    input_list <- list(...)
    
    if (length(input_list) > 0) {
      sim_vars <- sim_vars[order(sim_vars$variable),]
      input_list <- input_list[order(names(input_list))]
      
      # Identify which simulation variables are provided in the additional input.
      val_change <- which(sim_vars$variable %in% names(input_list))
      sim_vars[val_change, "value"] <- unlist(input_list)
      
      # Ensure that some specific variable values are kept as strings.
      sim_vars[sim_vars$variable=="population_size_phase2", "value"] <-
        paste0("'", sim_vars[sim_vars$variable=="population_size_phase2", "value"], "'")
      sim_vars[sim_vars$variable=="population_size_phase1", "value"] <-
        paste0("'", sim_vars[sim_vars$variable=="population_size_phase1", "value"], "'")
      sim_vars[sim_vars$variable=="dispersal_type_phase2", "value"] <- 
        paste0("'", sim_vars[sim_vars$variable=="dispersal_type_phase2", "value"], "'")
      sim_vars[sim_vars$variable=="dispersal_type_phase1", "value"] <- 
        paste0("'", sim_vars[sim_vars$variable=="dispersal_type_phase1", "value"], "'")
      sim_vars[sim_vars$variable=="natural_selection_model", "value"] <- 
        paste0("'", sim_vars[sim_vars$variable=="natural_selection_model", "value"], "'")
      
      vars_assign <- unlist(unname(
        mapply(paste, sim_vars$variable, "<-", sim_vars$value, SIMPLIFY = F)
      ))
      eval(parse(text = vars_assign))
    }
    
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
    
    # Identify loci with real data.
    real <- which(reference$type == "real")
    
    # Convert the neutral allele frequency to numeric.
    q_neutral <- as.numeric(ref_vars[ref_vars$variable=="q_neutral", "value"])
    
    # Check consistency between simulation variables and reference table variables.
    real_freq_table <- ref_vars[ref_vars$variable=="real_freq", "value"]
    if (real_freq_table != real_freq) {
      message(error(
        "  The value for the real_freq parameter was set differently in the simulations
   and in the creation of the reference table. They should be the same. 
   Please check it.\n"))
      stop()
    }
    
    real_loc_table <- ref_vars[ref_vars$variable=="real_loc", "value"]
    if (real_loc_table != real_loc) {
      message(error("  The value for the real_loc parameter was set differently in 
                the simulations and in the creation of the reference table. 
                They should be the same. Please check it.\n"))
      stop()
    }
    
    # Ensure that if real dataset values are required, the 'x' parameter is provided.
    if ((real_pops == TRUE | real_pop_size == TRUE | real_loc == TRUE | 
         real_freq == TRUE) && is.null(x)) {
      message(error(" The real dataset to extract information is missing\n"))
      stop()
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
    final_res <- rep(list(as.list(rep(NA, length(gen_store)))), number_iterations)
    
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
    # CLEAN UP STRING VARIABLES
    # -------------------------------
    # Remove extra quotes from character parameters.
    dispersal_type_phase2 <- gsub('\"', "", dispersal_type_phase2, fixed = TRUE)
    dispersal_type_phase1 <- gsub('\"', "", dispersal_type_phase1, fixed = TRUE)
    natural_selection_model <- gsub('\"', "", natural_selection_model, fixed = TRUE)
    chromosome_name <- gsub('\"', "", chromosome_name, fixed = TRUE)
    population_size_phase2 <- gsub('\"', "", population_size_phase2, fixed = TRUE)
    population_size_phase2 <- as.numeric(unlist(strsplit(population_size_phase2, " ")))
    population_size_phase1 <- gsub('\"', "", population_size_phase1, fixed = TRUE)
    population_size_phase1 <- as.numeric(unlist(strsplit(population_size_phase1, " ")))
    local_adap <- gsub('\"', "", local_adap, fixed = TRUE)
    local_adap <- as.numeric(unlist(strsplit(local_adap, " ")))
    clinal_adap <- gsub('\"', "", clinal_adap, fixed = TRUE)
    clinal_adap <- as.numeric(unlist(strsplit(clinal_adap, " ")))
    
    # -------------------------------
    # DETERMINE NUMBER OF POPULATIONS
    # -------------------------------
    if (phase1 == TRUE & real_pops == FALSE) {
      number_pops <- number_pops_phase1
    }
    if (phase1 == TRUE & real_pops == TRUE & !is.null(x)) {
      number_pops <- nPop(x)
    }
    if (phase1 == FALSE & real_pops == FALSE) {
      number_pops <- number_pops_phase2
    }
    if (phase1 == FALSE & real_pops == TRUE & !is.null(x)) {
      number_pops <- number_pops_phase2 <- nPop(x)
    }
    
    # -------------------------------
    # EXTRACT FREQUENCY INFORMATION FROM THE REAL DATA (IF APPLICABLE)
    # -------------------------------
    if (real_freq == TRUE & !is.null(x) & real_loc == TRUE) {
      pop_list_freq_temp <- seppop(x)
      loc_to_keep <- locNames(pop_list_freq_temp[[1]])[
        which(pop_list_freq_temp[[1]]$chromosome == chromosome_name)
      ]
      pop_list_freq_temp <- lapply(pop_list_freq_temp, gl.keep.loc, loc.list = loc_to_keep, verbose = 0)
      pop_list_freq <- lapply(pop_list_freq_temp, gl.alf)
    } 
    if (real_freq == TRUE & !is.null(x) & real_loc == FALSE) {
      pop_list_freq_temp <- seppop(x)
      pop_list_freq <- lapply(pop_list_freq_temp, gl.alf)
    }
    if (real_freq == FALSE) {
      pop_list_freq <- rep(NA, number_pops)
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
    # START SIMULATION ITERATION LOOP
    # -------------------------------
    for (iteration in 1:number_iterations) {
      if (iteration %% 1 == 0 & verbose >= 2) {
        message(report(" Starting iteration =", iteration, "\n"))
      }
      
      # -------------------------------
      # SETUP VARIABLES FOR PHASE 1 (IF APPLICABLE)
      # -------------------------------
      if (phase1 == TRUE) {
        if (real_pop_size == TRUE & !is.null(x)) {
          # Extract population sizes from the real data object.
          population_size_phase1 <- unname(unlist(table(pop(x))))
          # Ensure even population sizes by adjusting odd numbers.
          population_size_phase1 <- (population_size_phase1 %% 2 != 0) + population_size_phase1
        } else{
          population_size <- population_size_phase1
        }
        
        # Assign phase 1 specific simulation parameters.
        selection <- selection_phase1
        dispersal <- dispersal_phase1
        dispersal_type <- dispersal_type_phase1
        number_transfers <- number_transfers_phase1
        transfer_each_gen <- transfer_each_gen_phase1
        variance_offspring <- variance_offspring_phase1
        number_offspring <- number_offspring_phase1
        
        store_values <- store_phase1
        
        # Initialize counters for storing phase 1 results if required.
        if (store_phase1 == TRUE) {
          gen <- 0
          count_store <- 0
        }
        
        # Decide which sex(es) will be transferred based on number_transfers.
        if (number_transfers >= 2) {
          maletran <- TRUE
          femaletran <- TRUE
        } else if (number_transfers == 1) {
          maletran <- TRUE
          femaletran <- FALSE
        }
        
      } else {
        # -------------------------------
        # SETUP VARIABLES FOR PHASE 2
        # -------------------------------
        if (real_pop_size == TRUE & !is.null(x)) {
          population_size_phase2 <- unname(unlist(table(pop(x))))
          population_size_phase2 <- (population_size_phase2 %% 2 != 0) + population_size_phase2
          population_size <- population_size_phase2
        } else{
          population_size <- population_size_phase2
        }
      }
      
      # -------------------------------
      # ERROR CHECKS ON POPULATION NUMBERS BETWEEN PHASES
      # -------------------------------
      if (phase1 == TRUE & number_pops_phase1 != number_pops_phase2) {
        message(error("  Number of populations in phase 1 and phase 2 must be the same\n"))
        stop()
      }
      
      if (length(population_size_phase2) != number_pops_phase2) {
        message(error("  Number of entries for population sizes do not agree with 
           the number of populations for phase 2\n"))
        stop()
      }
      
      if (length(population_size_phase1) != number_pops_phase1 & phase1 == TRUE) {
        message(error("  Number of entries for population sizes do not agree with 
                  the number of populations for phase 1\n"))
        stop()
      }
      
      # -------------------------------
      # INITIALIZE POPULATIONS
      # -------------------------------
      if (verbose >= 2) {
        message(report("  Initialising populations\n"))
      }
      
      # Create a dummy function to bypass package checking.
      make_chr <- function(){}  
      
      # Define a C++ function to create chromosomes based on allele frequencies.
      Rcpp::cppFunction(plugins="cpp11",
                        
                        'StringVector make_chr(int j, NumericVector q) {
          StringVector out(j);
          int size = 1;
          IntegerVector x = IntegerVector::create(1, 0);
          bool rep = false;
          for (int i = 0; i < j; i++) {
            std::ostringstream temp;
            for (int z = 0; z < q.length(); z++) {
              NumericVector p = NumericVector::create(q[z], 1 - q[z]);
              temp << sample(x, size, rep, p);
            }
            out[i] = temp.str();
          }
          return out;
        }'
      )
      
      # Generate chromosomes for all individuals (each individual has two chromosomes).
      chr_temp <- make_chr(j = sum(population_size) * 2, q = reference$q)
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
        
        # If real frequency data is provided and location information is used:
        if (real_freq == TRUE & real_loc == TRUE) {
          for (individual_pop in 1:population_size[pop_n]) {
            q_prob_t <- pop_list_freq[[pop_n]]$alf1
            q_prob_t2 <- cbind(q_prob_t, 1 - q_prob_t)
            q_prob <- split(q_prob_t2, row(q_prob_t2))
            
            # Update chromosome 1 using allele frequency probabilities.
            stringi::stri_sub_all(pop[individual_pop, 3], from = real, length = 1) <- 
              mapply(function(y) { sample(x = c(1, 0), size = 1, prob = y, replace = FALSE) },
                     q_prob, USE.NAMES = FALSE)
            
            # Update chromosome 2 using allele frequency probabilities.
            stringi::stri_sub_all(pop[individual_pop, 4], from = real, length = 1) <- 
              mapply(function(y) { sample(x = c(1, 0), size = 1, prob = y, replace = FALSE) },
                     q_prob, USE.NAMES = FALSE)
          }
        }
        
        # If only real frequency (without location) is used.
        if (real_freq == TRUE & real_loc == FALSE) {
          for (individual_pop in 1:population_size[pop_n]) {
            q_prob_t <- pop_list_freq[[pop_n]]$alf1
            q_prob_t2 <- cbind(q_prob_t, 1 - q_prob_t)
            q_prob <- split(q_prob_t2, row(q_prob_t2))
            
            stringi::stri_sub_all(pop[individual_pop, 3], from = real, length = 1) <- 
              mapply(function(y) { sample(x = c(1, 0), size = 1, prob = y, replace = FALSE) },
                     q_prob, USE.NAMES = FALSE)
            stringi::stri_sub_all(pop[individual_pop, 4], from = real, length = 1) <- 
              mapply(function(y) { sample(x = c(1, 0), size = 1, prob = y, replace = FALSE) },
                     q_prob, USE.NAMES = FALSE)
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
          
          # Set transfer flags for phase 2.
          if (number_transfers >= 2) {
            maletran <- TRUE
            femaletran <- TRUE
          } else if (number_transfers == 1) {
            maletran <- TRUE
            femaletran <- FALSE
          }
          
          # Initialize counters for storing phase 2 results.
          if (!exists("count_store")) {
            count_store <- 0
          }
          gen <- 0
          
          # Resample populations if the simulation requires it.
          if (phase1 == TRUE) {
            if (same_line == TRUE) {
              # Sample from one population and apply to all.
              pop_sample <- sample(pops_vector, 1)
              pop_list_temp <- lapply(pops_vector, function(x) {
                pop_temp <- rbind(
                  pop_list[[pop_sample]][sample(which(pop_list[[pop_sample]]$V1 == "Male"),
                                                size = population_size[x] / 2, replace = TRUE),],
                  pop_list[[pop_sample]][sample(which(pop_list[[pop_sample]]$V1 == "Female"),
                                                size = population_size[x] / 2, replace = TRUE),]
                )
                pop_temp$V2 <- x
                return(pop_temp)
              })
              pop_list <- pop_list_temp
            }
            if (same_line == FALSE) {
              pop_list <- lapply(pops_vector, function(x) {
                pop_temp <- rbind(
                  pop_list[[x]][sample(which(pop_list[[x]]$V1 == "Male"),
                                       size = population_size[x] / 2, replace = TRUE),],
                  pop_list[[x]][sample(which(pop_list[[x]]$V1 == "Female"),
                                       size = population_size[x] / 2, replace = TRUE),]
                )
                pop_temp$V2 <- x
                return(pop_temp)
              })
            }
          }
        }
        
        # Increment generation counter if storing values.
        if (store_values == TRUE) {
          gen <- gen + 1
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
            
            # Set additional dispersal parameters.
            dispersal_pairs$number_transfers <- number_transfers
            dispersal_pairs$transfer_each_gen <- transfer_each_gen
            
          } else {
            # If a dispersal file is provided, read the file.
            dispersal_pairs <- suppressWarnings(read.csv(file_dispersal))
          }
          
          # Define population sizes for each pair.
          if (real_pop_size == TRUE) {
            dispersal_pairs$size_pop1 <- unname(unlist(table(pop(x))))[dispersal_pairs$pop1]
            dispersal_pairs$size_pop2 <- unname(unlist(table(pop(x))))[dispersal_pairs$pop2]
          } else {
            dispersal_pairs$size_pop1 <- population_size[dispersal_pairs$pop1]
            dispersal_pairs$size_pop2 <- population_size[dispersal_pairs$pop2]
          }
          
          # Process each dispersal pair.
          for (dis_pair in 1:nrow(dispersal_pairs)) {
            res <- migration(
              population1 = pop_list[[dispersal_pairs[dis_pair, "pop1"]]],
              population2 = pop_list[[dispersal_pairs[dis_pair, "pop2"]]],
              gen = generation,
              size_pop1 = dispersal_pairs$size_pop1[dis_pair],
              size_pop2 = dispersal_pairs$size_pop2[dis_pair],
              trans_gen = dispersal_pairs$transfer_each_gen[dis_pair],
              n_transfer = dispersal_pairs$number_transfers[dis_pair],
              male_tran = maletran,
              female_tran = femaletran
            )
            
            # Update populations after migration.
            pop_list[[dispersal_pairs[dis_pair, "pop1"]]] <- res[[1]]
            pop_list[[dispersal_pairs[dis_pair, "pop1"]]]$V2 <- dispersal_pairs[dis_pair, "pop1"]
            pop_list[[dispersal_pairs[dis_pair, "pop2"]]] <- res[[2]]
            pop_list[[dispersal_pairs[dis_pair, "pop2"]]]$V2 <- dispersal_pairs[dis_pair, "pop2"]
            maletran <- res[[3]]
            femaletran <- res[[4]]
          }
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
          tmp_rep$id <- paste0(generation,"_",x,"_",1:nrow(tmp_rep))
          return(tmp_rep)
        })
        
        # -------------------------------
        # MUTATION PHASE
        # -------------------------------
        if (mutation == TRUE) {
          for (off_pop in 1:length(offspring_list)) {
            offspring_pop <- offspring_list[[off_pop]]
            # Add a uniform random number to each offspring for mutation decision.
            offspring_pop$runif <- runif(nrow(offspring_pop))
            
            for (offspring_ind in 1:nrow(offspring_pop)) {
              if (length(mutation_loci_location) == 0 & verbose >= 2) {
                message(important("  No more locus to mutate\n"))
                break()
              }
              
              # Determine if a mutation occurs based on the mutation rate.
              if (offspring_pop[offspring_ind, "runif"] < mut_rate) {
                locus_to_mutate <- sample(mutation_loci_location, 1)
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
            pops_clinal <- seq(clinal_adap[1], clinal_adap[2], 1)
            reference_clinal_temp <- replicate(length(pops_vector), 
                                               reference[, c("s", "h")], 
                                               simplify = FALSE)
            # Calculate clinal selection coefficients.
            clinal_s <- 1 - c(0, (1:(length(pops_clinal)-1) * (clinal_strength / 100)))
            reference_clinal <- lapply(pops_clinal, function(y) {
              reference_clinal_temp[[y]][adv_loci, "s"] <-  
                reference_clinal_temp[[y]][adv_loci, "s"] * clinal_s[y]
              return(reference_clinal_temp[[y]])
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
        # Check if any population went extinct (e.g., insufficient males or females).
        test_extinction <- unlist(lapply(pops_vector, function(x) {
          length(which(offspring_list[[x]]$V1 == "Male")) < population_size / 2 |
            length(which(offspring_list[[x]]$V1 == "Female")) < population_size / 2
        }))
        
        if (any(test_extinction == TRUE) & verbose >= 2) {
          message(important(" One Population became EXTINCT at generation", generation, "\n"))
          message(important("  Breaking this iteration and passing to the next iteration", "\n"))
          
          if (sample_percent != 100) {
            population_size_temp <- round(population_size * (sample_percent / 100))
            # Ensure even numbers after subsampling.
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
          
          # Combine simulation and reference variables for saving.
          s_vars_temp <- rbind(ref_vars, sim_vars)
          s_vars_temp <- setNames(data.frame(t(s_vars_temp[,-1])), s_vars_temp[, 1])
          s_vars_temp$generation <- generation
          s_vars_temp$iteration <- iteration
          s_vars_temp$seed <- seed
          s_vars_temp$del_ind_cM <- density_mutations_per_cm
          s_vars_temp$sample_percent <- sample_percent
          s_vars_temp$file_dispersal <- file_dispersal
          
          # Store the current generation's results.
          final_res[[iteration]][[count_store]] <- store(
            p_vector = pops_vector,
            p_size = population_size_temp,
            p_list = pop_list_temp,
            n_loc_1 = loci_number,
            ref = reference,
            p_map = plink_map,
            s_vars = s_vars_temp
          )
          
          if (real_pops == TRUE) {
            popNames(final_res[[iteration]][[count_store]]) <- popNames(x)
          }
          break()  # End the current iteration if extinction occurs.
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
          
          # Define a dummy function to bypass package checking.
          make_freqs <- function(){}  
          
          # Define a C++ function to compute allele frequencies across loci.
          Rcpp::cppFunction(plugins="cpp11",
                            "NumericVector make_freqs(StringVector seqs) {
              int seqN = seqs.length();
              int locN = strlen(seqs(0));
              NumericMatrix freq_mat = NumericMatrix(seqN, locN);
              NumericVector out(locN);
              for (int i = 0; i < seqN; i++) {
                for (int j = 0; j < locN; j++) {
                  freq_mat(i, j) = seqs(i)[j] - '0';
                }
              }
              for (int y = 0; y < locN; y++){
                out[y] = sum(freq_mat(_, y));
              }
              return out;
            }"
          )
          
          # Get frequencies for each locus.
          freqs <- make_freqs(pops_seqs)
          deleterious_eliminated <- which(freqs == 0)
          # Recycle loci with no deleterious mutations.
          mutation_loci_location <- union(mutation_loci_location, deleterious_eliminated)
        }
        
        # -------------------------------
        # STORE GENERATION RESULTS INTO GENLIGHT OBJECTS
        # -------------------------------
        if (generation %in% gen_store & exists("count_store")) {
          count_store <- count_store + 1
          
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
          final_res[[iteration]][[count_store]] <- store(
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
            popNames(final_res[[iteration]][[count_store]]) <- popNames(x)
          } else {
            popNames(final_res[[iteration]][[count_store]]) <- as.character(pops_vector)
          }
        }
      }  # End generation loop
    }  # End iteration loop
    
    # -------------------------------
    # FINALIZE RESULTS
    # -------------------------------
    # Name the elements of the final results list.
    names(final_res) <- paste0("iteration_", 1:number_iterations)
    final_res <- lapply(final_res, function(x) {
      names(x) <- paste0("generation_", gen_store)
      return(x)
    })
    
    # Remove any NA entries from the results.
    final_res <- lapply(final_res, function(x) {
      x[!is.na(x)]
    })
    
    # -------------------------------
    # FLAG THE END OF THE FUNCTION EXECUTION
    # -------------------------------
    if (verbose >= 1) {
      message(report("Completed:", funname, "\n"))
    }
    
    # Return the final simulation results invisibly.
    return(invisible(final_res))
  }
