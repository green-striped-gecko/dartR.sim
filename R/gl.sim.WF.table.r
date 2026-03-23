#' @name gl.sim.WF.table
#' @title Creates the reference table for running gl.sim.WF.run
#' @description
#' This function creates a reference table to be used as input for the function
#'  \code{\link{gl.sim.WF.run}}. The created table has eight columns with the 
#'  following information for each locus to be simulated:
#' \itemize{ 
#' \item q - initial frequency.
#' \item h - dominance coefficient.
#' \item s - selection coefficient.
#' \item c - recombination rate.
#' \item loc_bp - chromosome location in base pairs.
#' \item loc_cM - chromosome location in centiMorgans.
#' \item chr_name - chromosome name.
#' \item type - SNP type.
#' } 
#' 
#' The reference table can be further modified as required. 
#' 
#' See documentation and tutorial for a complete description of the simulations.
#' These documents can be accessed at http://georges.biomatix.org/dartR 
#' 
#' @param file_var Path of the variables file 'ref_variables.csv' (see details) 
#' [required if interactive_vars = FALSE].
#' @param x Name of the genlight object containing the SNP data to extract
#' values for some simulation variables (see details) [default NULL].
#' @param file_targets_sel Path of the file with the targets for selection (see 
#' details) [default NULL].
#' @param file_r_map Path of the file with the recombination map (see details)
#' [default NULL].
#' @param interactive_vars Run a shiny app to input interactively the values of
#'  simulation variables [default TRUE].
#' @param seed Set the seed for the simulations [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @param ... Any variable and its value can be added separately within the 
#' function, will be changed over the input value supplied by the csv file. See 
#' tutorial. 
#' @details
#' Values for the variables to create the reference table can be submitted into 
#' the function interactively through a Shiny app if interactive_vars = TRUE. 
#' Optionally, if interactive_vars = FALSE, values for variables can be 
#' submitted by using the csv file 'ref_variables.csv' which can be found by 
#' typing in the R console:
#'  system.file('extdata', 'ref_variables.csv', package ='dartR.data').
#'  
#' The values of the variables can be modified using the third column (“value”) 
#' of this file. 
#' 
#' If a genlight object is used as input for some of the simulation variables, 
#' this function access the information stored in the slots x$position and 
#' x$chromosome.
#' 
#' Examples of the format required for the recombination map file and the 
#' targets for selection file can be found by typing in the R console:
#' \itemize{ 
#' \item system.file('extdata', 'fly_recom_map.csv', package ='dartR.data')
#' \item system.file('extdata', 'fly_targets_of_selection.csv', package ='dartR.data')
#' }
#' 
#' To show further information of the variables in interactive mode, it might be
#'  necessary to call first: 'library(shinyBS)' for the information to be 
#'  displayed.
#' @return Returns a list with the reference table used as input for the function
#'  \code{\link{gl.sim.WF.run}} and a table with the values variables used to 
#'  create the reference table.
#' @author Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' ref_table <- gl.sim.WF.table(file_var=system.file("extdata", 
#' "ref_variables.csv", package = "dartR.sim"),interactive_vars = FALSE)
#' 
#' res_sim <- gl.sim.WF.run(file_var = system.file("extdata",
#'  "sim_variables.csv", package ="dartR.sim"),ref_table=ref_table,
#'  interactive_vars = FALSE)
#'  
#' @seealso \code{\link{gl.sim.WF.run}}
#' @family simulation functions
#' @rawNamespace import(fields, except = flame)
#' @export

gl.sim.WF.table <- function(file_var, 
                            x = NULL, 
                            file_targets_sel = NULL, 
                            file_r_map = NULL,
                            interactive_vars = TRUE, 
                            seed = NULL,
                            verbose = NULL,
                            ...) {
  
  ## Set the seed if one is provided, to ensure reproducible simulations
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  ## Set the verbosity level using a helper function
  verbose <- gl.check.verbosity(verbose)
  
  ## Flag the start of the function execution (for logging purposes)
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "Jody",
                   verbose = verbose)
  
  ## LOADING VARIABLES
  ##### If running interactively, open a shiny app to input simulation variables #####
  if (interactive_vars == TRUE) {
    ref_vars <- interactive_reference()
    
    # For certain variables, wrap the values in single quotes (for later evaluation)
    ref_vars[ref_vars$variable=="chromosome_name", "value"] <-
      paste0("'", ref_vars[ref_vars$variable=="chromosome_name", "value"], "'")
    ref_vars[ref_vars$variable=="h_distribution_del", "value"] <-
      paste0("'", ref_vars[ref_vars$variable=="h_distribution_del", "value"], "'")
    ref_vars[ref_vars$variable=="s_distribution_del", "value"] <-
      paste0("'", ref_vars[ref_vars$variable=="s_distribution_del", "value"], "'")
    ref_vars[ref_vars$variable=="q_distribution_del", "value"] <-
      paste0("'", ref_vars[ref_vars$variable=="q_distribution_del", "value"], "'")
    ref_vars[ref_vars$variable=="h_distribution_adv", "value"] <-
      paste0("'", ref_vars[ref_vars$variable=="h_distribution_adv", "value"], "'")
    ref_vars[ref_vars$variable=="s_distribution_adv", "value"] <-
      paste0("'", ref_vars[ref_vars$variable=="s_distribution_adv", "value"], "'")
    ref_vars[ref_vars$variable=="q_distribution_adv", "value"] <-
      paste0("'", ref_vars[ref_vars$variable=="q_distribution_adv", "value"], "'")
    
    ## Order the variables alphabetically
    ref_vars <- ref_vars[order(ref_vars$variable),]
    
    ## Create assignment strings (e.g., variable <- value) for each simulation variable
    vars_assign <- unlist(unname(
      mapply(paste, ref_vars$variable, "<-",
             ref_vars$value, SIMPLIFY = F)
    ))
    ## Evaluate the assignments to set variables in the environment
    eval(parse(text = vars_assign))
    
  } else {
    ## If not interactive, read the variables from a CSV file
    ref_vars <- suppressWarnings(read.csv(file_var))
    ref_vars <- ref_vars[, 2:3]  # select only the variable names and their values
    ref_vars <- ref_vars[order(ref_vars$variable),]
    
    ## Create and evaluate assignment strings from the CSV file values
    vars_assign <- unlist(unname(
      mapply(paste, ref_vars$variable, "<-",
             ref_vars$value, SIMPLIFY = F)
    ))
    eval(parse(text = vars_assign))
  }
  
  ## Process additional input arguments passed through ...
  input_list <- list(...)
  if (length(input_list) > 0) {
    ref_vars <- ref_vars[order(ref_vars$variable),]
    input_list <- input_list[order(names(input_list))]
    ## Identify which variables in ref_vars are being overridden
    val_change <- which(ref_vars$variable %in% names(input_list))
    ref_vars[val_change, "value"] <- unlist(input_list)
    
    ## Ensure certain variables have values wrapped in single quotes
    ref_vars[ref_vars$variable=="chromosome_name", "value"] <-
      paste0("'", ref_vars[ref_vars$variable=="chromosome_name", "value"], "'")
    ref_vars[ref_vars$variable=="h_distribution_del", "value"] <-
      paste0("'", ref_vars[ref_vars$variable=="h_distribution_del", "value"], "'")
    ref_vars[ref_vars$variable=="s_distribution_del", "value"] <-
      paste0("'", ref_vars[ref_vars$variable=="s_distribution_del", "value"], "'")
    ref_vars[ref_vars$variable=="q_distribution_del", "value"] <-
      paste0("'", ref_vars[ref_vars$variable=="q_distribution_del", "value"], "'")
    ref_vars[ref_vars$variable=="h_distribution_adv", "value"] <-
      paste0("'", ref_vars[ref_vars$variable=="h_distribution_adv", "value"], "'")
    ref_vars[ref_vars$variable=="s_distribution_adv", "value"] <-
      paste0("'", ref_vars[ref_vars$variable=="s_distribution_adv", "value"], "'")
    ref_vars[ref_vars$variable=="q_distribution_adv", "value"] <-
      paste0("'", ref_vars[ref_vars$variable=="q_distribution_adv", "value"], "'")
    
    ## Reassign the variables with updated values
    vars_assign <- unlist(unname(
      mapply(paste, ref_vars$variable, "<-",
             ref_vars$value, SIMPLIFY = F)
    ))
    eval(parse(text = vars_assign))
  }
  
  ## Clean up any extra quotation marks from the distribution parameters
  s_distribution_del <- gsub('\"', "", s_distribution_del, fixed = TRUE)
  s_distribution_adv <- gsub('\"', "", s_distribution_adv, fixed = TRUE)
  h_distribution_del <- gsub('\"', "", h_distribution_del, fixed = TRUE)
  h_distribution_adv <- gsub('\"', "", h_distribution_adv, fixed = TRUE)
  q_distribution_del <- gsub('\"', "", q_distribution_del, fixed = TRUE)
  q_distribution_adv <- gsub('\"', "", q_distribution_adv, fixed = TRUE)
  
  ##### LOADING INFORMATION #####
  ## RECOMBINATION MAP: load recombination data if provided
  if (!is.null(file_r_map)) {
    map <- read.csv(file_r_map)
    map$Chr <- as.character(map$Chr)
    
    ## Check that the chromosome name is present in the recombination map
    if (!chromosome_name %in% map$Chr) {
      message(error("  Chromosome name is not in the recombination map file\n"))
      stop()
    }
    ## Subset the map for the given chromosome and adjust units for cM
    map <- map[which(map$Chr == chromosome_name),]
    map$cM <- map$cM / 100
    ## Replace any NA values with 0
    map[is.na(map$cM), ] <- 0
  } else {
    ## If no recombination map file is provided, create a default map using chunk data
    map <- as.data.frame(matrix(nrow = chunk_number))
    map[, 1] <- chunk_cM / 100
    colnames(map) <- "cM"
  }
  
  ## TARGETS OF SELECTION: load targets file if provided
  targets_temp <- NULL
  if (!is.null(file_targets_sel)) {
    targets_temp <- read.csv(file_targets_sel)
    targets_temp$chr_name <- as.character(targets_temp$chr_name)
    
    ## Check that the chromosome name exists in the targets file
    if (!chromosome_name %in% targets_temp$chr_name) {
      message(error("  Chromosome name is not in the targets of selection file\n"))
      stop()
    }
    targets_temp <- targets_temp[which(targets_temp$chr_name == chromosome_name),]
  }
  
  ## REAL DATA: ensure that if real location or frequency info is needed, the dataset is provided
  if ((real_loc == TRUE | real_freq == TRUE) && is.null(x)) {
    message(error(" The real dataset to extract information is missing\n"))
    stop()
  }
  location_real_temp <- NULL
  if (real_loc == TRUE & !is.null(x)) {
    if (!chromosome_name %in% x@chromosome) {
      message(error("  Chromosome name is not in the genlight object\n"))
      stop()
    }
    ## Extract chromosome and position data from the genlight object
    location_real_temp <- as.data.frame(cbind(as.character(x$chromosome), x$position))
    colnames(location_real_temp) <- c("chr", "pos")
    
    location_real_temp <- location_real_temp[location_real_temp$chr == chromosome_name, ]
    location_real_temp <- as.numeric(location_real_temp[, "pos"])
    location_real_temp <- location_real_temp[order(location_real_temp)]
  }
  
  ##### CHROMOSOME LENGTH #####
  ## Determine chromosome length based on recombination map or default chunks
  if (!is.null(file_r_map)) {
    chr_length <- tail(map$to, 1)
  } else {
    chr_length <- chunk_number * chunk_bp
  }
  
  ## Adjust chromosome length and chunk size if real locations are provided
  if (real_loc == TRUE & !is.null(x) & is.null(file_r_map)) {
    chr_length <- tail(location_real_temp, 1)
    chunk_bp <- chr_length / chunk_number
  }
  
  ## Similarly adjust for targets of selection if provided and no recombination map
  if (!is.null(file_targets_sel) & is.null(file_r_map)) {
    chr_length <- tail(targets_temp$end, 1)
    chunk_bp <- chr_length / chunk_number
  }
  
  ##### LOCATIONS ##########
  ## Process real dataset locations
  location_real_bp <- NULL
  if (real_loc == TRUE & !is.null(x)) {
    location_real_bp <- location_real_temp
    location_real_bp <- location_real_bp[order(location_real_bp)]
    ## Round to the nearest ten and add 1 to each location
    location_real_bp <- round(location_real_bp, -1)
    location_real_bp <- location_real_bp + 1
  } 
  
  ## If only real frequencies are provided, generate a sequence of positions
  if (real_loc == FALSE & real_freq == TRUE & !is.null(x)) {
    location_real_temp <- round(seq(chunk_bp / (nLoc(x) + 1),
                                    (chunk_number * chunk_bp),
                                    chunk_bp / nLoc(x)))
    location_real_bp <- sample(location_real_temp, size = nLoc(x))
    location_real_bp <- location_real_bp[order(location_real_bp)]
    location_real_bp <- round(location_real_bp, -1)
    location_real_bp <- location_real_bp + 1
  }
  
  ## Neutral loci simulations: create locations for neutral loci across chunks
  location_neutral_bp <- NULL
  if (chunk_neutral_loci > 0) {
    location_neutral_bp <- round(seq(chunk_bp / (chunk_neutral_loci + 1),
                                     (chunk_number * chunk_bp),
                                     chunk_bp / chunk_neutral_loci))
    location_neutral_bp <- round(location_neutral_bp, -1)
    location_neutral_bp <- location_neutral_bp + 2
  }
  
  ## Deleterious loci: set up locations either from the targets file or default values
  location_deleterious_bp <- NULL
  if (!is.null(file_targets_sel) | loci_deleterious > 0) {
    if (!is.null(file_targets_sel)) {
      del <- targets_temp
      del$targets <- ceiling(del$targets * (deleterious_factor / 100))
      del$distance <-  del$end - del$start
    } 
    ## If no targets file, create default intervals and assign targets
    if (is.null(file_targets_sel)) {
      del <- as.data.frame(matrix(nrow = chunk_number, ncol = 3))
      colnames(del) <- c("start", "end", "targets")
      del$start <- seq(1, chr_length, (chr_length / chunk_number))
      del$end <- seq((chr_length / chunk_number), chr_length, (chr_length / chunk_number))
      if (loci_deleterious < chunk_number) {
        row_targets <- sample(1:chunk_number, size = loci_deleterious)
        del[row_targets, "targets"] <- 1
        del[is.na(del$targets), "targets"] <- 0
      } else {
        del$targets <- round(loci_deleterious / chunk_number)
      }
      del$distance <-  del$end - del$start
    }
    sample_resolution <- round(mean(del$distance) / max(del$targets) / 10)
    ## For each interval, sample positions for deleterious loci
    for (i in 1:nrow(del)) {
      location_deleterious_temp <- mapply(
        FUN = function(a, b) {
          seq(from = a, to = b, by = sample_resolution)
        },
        a = unname(unlist(del[i, "start"])),
        b = unname(unlist(del[i, "end"]))
      )
      location_deleterious_temp <- as.vector(round(location_deleterious_temp))
      location_deleterious_temp <- sample(location_deleterious_temp, size = del[i, "targets"])
      location_deleterious_bp <- c(location_deleterious_bp, location_deleterious_temp)
    }
    location_deleterious_bp <- location_deleterious_bp[order(location_deleterious_bp)]
    ## Round positions and adjust with an offset
    location_deleterious_bp <- round(location_deleterious_bp, -1)
    location_deleterious_bp <- location_deleterious_bp + 3
    loci_deleterious <- length(location_deleterious_bp)
  }
  
  ## Advantageous loci: sample positions randomly and adjust them
  location_advantageous_bp <- NULL
  if (loci_advantageous > 0) {
    location_advantageous_bp <- sample(1:chr_length, loci_advantageous)
    location_advantageous_bp <- location_advantageous_bp[order(location_advantageous_bp)]
    location_advantageous_bp <- round(location_advantageous_bp, -1)
    location_advantageous_bp <- location_advantageous_bp + 4
  }
  
  ## Mutations: determine mutation loci positions for neutral, deleterious, and advantageous mutations
  loci_mutation <- loci_mut_neu + loci_mut_del + loci_mut_adv
  location_mutations_bp <- NULL
  if (!is.null(file_targets_sel) | (loci_mutation > 0)) {
    if (!is.null(file_targets_sel)) {
      mutations <- targets_temp
      mutations$targets <- ceiling(mutations$targets * (mutations_factor / 100))
      mutations$distance <-  mutations$end - mutations$start
    } else {
      mutations <- as.data.frame(matrix(nrow = chunk_number, ncol = 3))
      colnames(mutations) <- c("start", "end", "targets")
      mutations$start <- seq(1, chr_length, (chr_length / chunk_number))
      mutations$end <- seq((chr_length / chunk_number), chr_length, (chr_length / chunk_number))
      mutations$targets <- ceiling(loci_mutation / chunk_number)
      mutations$distance <-  mutations$end - mutations$start
    }
    sample_resolution <- round(mean(mutations$distance) / max(mutations$targets) / 20)
    ## Sample mutation positions for each interval
    for (i in 1:nrow(mutations)) {
      location_mutations_temp <- mapply(
        FUN = function(a, b) {
          seq(from = a, to = b, by = sample_resolution)
        },
        a = unname(unlist(mutations[i, "start"])),
        b = unname(unlist(mutations[i, "end"]))
      )
      location_mutations_temp <- as.vector(round(location_mutations_temp))
      location_mutations_temp <- sample(location_mutations_temp, size = mutations[i, "targets"])
      location_mutations_bp <- c(location_mutations_bp, location_mutations_temp)
    }
    location_mutations_bp <- location_mutations_bp[order(location_mutations_bp)]
    location_mutations_bp <- round(location_mutations_bp, -1)
    location_mutations_bp <- location_mutations_bp + 5
    location_mutations_bp <- location_mutations_bp[sample(1:length(location_mutations_bp), size = loci_mutation)]
  }
  
  ## Ensure uniqueness and update counts for deleterious, advantageous, and mutation loci
  location_deleterious_bp <- unique(location_deleterious_bp)
  loci_deleterious <- length(location_deleterious_bp)
  location_advantageous_bp <- unique(location_advantageous_bp)
  loci_advantageous <- length(location_advantageous_bp)
  location_mutations_bp <- unique(location_mutations_bp)
  loci_mutation <- length(location_mutations_bp)
  
  ## Combine all loci positions and order them
  location_loci_bp <- c(location_real_bp, location_neutral_bp,
                        location_deleterious_bp, location_advantageous_bp,
                        location_mutations_bp)
  location_loci_bp <- location_loci_bp[order(location_loci_bp)]
  
  ## Check that the number of loci exceeds the number of genome chunks
  if (chunk_number > length(location_loci_bp)) {
    message(error("  Number of loci should be more than the number of genome chunks\n"))
    stop()
  }
  
  total_loci <- length(location_loci_bp)
  
  ##### RECOMBINATION MAP #####
  ## Build a recombination map by cross-multiplying loci positions with the map data
  recombination_map_temp <- map
  recombination_map_temp$midpoint <- seq(chunk_bp / 2, chr_length, chunk_bp)[1:nrow(recombination_map_temp)]
  recombination_temp <- unlist(lapply(location_loci_bp, findInterval, vec = as.numeric(paste(unlist(recombination_map_temp$midpoint)))))
  
  ## Correct intervals that fall before the first midpoint
  recombination_temp[recombination_temp == 0] <- 1
  recombination_2 <- recombination_map_temp[recombination_temp, "cM"]
  recombination_map <- as.data.frame(cbind(location_loci_bp, recombination_2))
  recombination_map$c <- NA
  
  ## Calculate recombination rates for each interval (except the last)
  for (target_row in 1:(nrow(recombination_map) - 1)) {
    recombination_map[target_row, "c"] <- ((recombination_map[target_row + 1, "location_loci_bp"] - 
                                              recombination_map[target_row, "location_loci_bp"]) * 
                                             recombination_map[target_row, "recombination_2"]) / chunk_bp
  }
  ## Set recombination rate for the last locus to zero (to avoid function crash)
  recombination_map[nrow(recombination_map), "c"] <- 0
  recombination_map$accum <- cumsum(recombination_map[, "c"])
  
  ## Identify rows corresponding to different types of loci in the recombination map
  location_neutral_row <- NULL
  if (chunk_neutral_loci > 0) {
    location_neutral_row <- lapply(location_neutral_bp, function(x) {
      which(recombination_map$location_loci == x)
    })
    location_neutral_row <- unname(unlist(location_neutral_row))
  }
  
  location_real_row <- NULL
  if (real_loc == TRUE | real_freq == TRUE) {
    location_real_row <- lapply(location_real_bp, function(x) {
      which(recombination_map$location_loci == x)
    })
    location_real_row <- unname(unlist(location_real_row))
  }
  
  location_deleterious_row <- NULL
  if (loci_deleterious > 0) {
    location_deleterious_row <- lapply(location_deleterious_bp, function(x) {
      which(recombination_map$location_loci == x)
    })
    location_deleterious_row <- unname(unlist(location_deleterious_row))
  }
  
  location_advantageous_row <- NULL
  if (loci_advantageous > 0) {
    location_advantageous_row <- lapply(location_advantageous_bp, function(x) {
      which(recombination_map$location_loci == x)
    })
    location_advantageous_row <- unname(unlist(location_advantageous_row))
  }
  
  location_mutations_row <- NULL
  if (loci_mutation > 0) {
    location_mutations_row  <- lapply(location_mutations_bp, function(x) {
      which(recombination_map$location_loci == x)
    })
    location_mutations_row <- unname(unlist(location_mutations_row))
  }
  
  ##### REFERENCE TABLE ########
  ## Split mutation rows into categories: neutral, deleterious, and advantageous
  mut_tmp <- location_mutations_row
  if (loci_mut_neu > 0) {
    neu_tmp <- sample(1:length(mut_tmp), size = loci_mut_neu)
    location_mut_neu_row <- mut_tmp[neu_tmp]
    location_mut_neu_row <- location_mut_neu_row[order(location_mut_neu_row)]
    mut_tmp <- mut_tmp[-neu_tmp]
  }
  
  if (loci_mut_del > 0) {
    del_tmp <- sample(1:length(mut_tmp), size = loci_mut_del)
    location_mut_del_row <- mut_tmp[del_tmp]
    location_mut_del_row <- location_mut_del_row[order(location_mut_del_row)]
    mut_tmp <- mut_tmp[-del_tmp]
  }
  
  if (loci_mut_adv > 0) {
    adv_tmp <- sample(1:length(mut_tmp), size = loci_mut_adv)
    location_mut_adv_row <- mut_tmp[adv_tmp]
    location_mut_adv_row <- location_mut_adv_row[order(location_mut_adv_row)]
    mut_tmp <- mut_tmp[-adv_tmp]
  }
  
  ## Initialize vectors for selection coefficient (s), dominance (h), and frequency (q)
  s <- rep(NA, total_loci)
  h <- rep(NA, total_loci)
  q <- rep(NA, total_loci)
  
  ##### SELECTION COEFFICIENT (s) #####
  ## Assign selection coefficients based on locus type
  if (chunk_neutral_loci > 0) {
    s[location_neutral_row] <- 0
  }
  if (real_loc == TRUE) {
    s[location_real_row] <- 0
  }
  if (loci_deleterious > 0) {
    if (s_distribution_del == "equal") {
      s[location_deleterious_row] <- s_del
    }
    if (s_distribution_del == "gamma") {
      s[location_deleterious_row] <- rgamma(loci_deleterious, shape = gamma_shape, scale = gamma_scale)
    }
    if (s_distribution_del == "log_normal") {
      s[location_deleterious_row] <- rlnorm(loci_deleterious, meanlog = log(log_mean), sdlog = log(log_sd))
    }
  }
  if (loci_advantageous > 0) {
    if (s_distribution_adv == "equal") {
      s[location_advantageous_row] <- s_adv * -1
    }
    if (s_distribution_adv == "exponential") {
      s[location_advantageous_row] <- rexp(loci_advantageous, rate = exp_rate) * -1
    }
  }
  if (loci_mut_neu > 0) {
    s[location_mut_neu_row] <- 0
  }
  if (loci_mut_del > 0) {
    if (s_distribution_del == "equal") {
      s[location_mut_del_row] <- s_del
    }
    if (s_distribution_del == "gamma") {
      s[location_mut_del_row] <- rgamma(loci_mut_del, shape = gamma_shape, scale = gamma_scale)
    }
    if (s_distribution_del == "log_normal") {
      s[location_mut_del_row] <- rlnorm(loci_mut_del, meanlog = log(log_mean), sdlog = log(log_sd))
    }
  }
  if (loci_mut_adv > 0) {
    if (s_distribution_adv == "equal") {
      s[location_mut_adv_row] <- s_adv * -1
    }
    if (s_distribution_adv == "exponential") {
      s[location_mut_adv_row] <- rexp(loci_mut_adv, rate = exp_rate) * -1
    }
  }
  
  ##### DOMINANCE (h) #####
  ## Assign dominance values based on locus type and distribution
  if (chunk_neutral_loci > 0) {
    h[location_neutral_row] <- 0
  }
  if (real_loc == TRUE) {
    h[location_real_row] <- 0
  }
  if (loci_deleterious > 0) {
    if (h_distribution_del == "equal") {
      h[location_deleterious_row] <- h_del
    }
    if (h_distribution_del == "normal") {
      h[location_deleterious_row] <- rnorm(loci_deleterious, mean = h_mean_del, sd = h_sd_del)
    }
    if (h_distribution_del == "equation") {
      h[location_deleterious_row] <- 1 / ((1 / h_intercept_del) - (-1 * h_rate_del * abs(s[location_deleterious_row])))
    }
  }
  if (loci_advantageous > 0) {
    if (h_distribution_adv == "equal") {
      h[location_advantageous_row] <- h_adv
    }
    if (h_distribution_adv == "normal") {
      h[location_advantageous_row] <- rnorm(loci_advantageous, mean = h_mean_adv, sd = h_sd_adv)
    }
    if (h_distribution_adv == "equation") {
      h[location_advantageous_row] <- 1 / ((1 / h_intercept_adv) - (-1 * h_rate_adv * abs(s[location_advantageous_row])))
    }
  }
  if (loci_mut_neu > 0) {
    h[location_mut_neu_row] <- 0
  }
  if (loci_mut_del > 0) {
    if (h_distribution_del == "equal") {
      h[location_mut_del_row] <- h_del
    }
    if (h_distribution_del == "normal") {
      h[location_mut_del_row] <- rnorm(loci_mut_del, mean = h_mean_del, sd = h_sd_del)
    }
    if (h_distribution_del == "equation") {
      h[location_mut_del_row] <- 1 / ((1 / h_intercept_del) - (-1 * h_rate_del * abs(s[location_mut_del_row])))
    }
  }
  if (loci_mut_adv > 0) {
    if (h_distribution_adv == "equal") {
      h[location_mut_adv_row] <- h_adv
    }
    if (h_distribution_adv == "normal") {
      h[location_mut_adv_row] <- rnorm(loci_mut_adv, mean = h_mean_adv, sd = h_sd_adv)
    }
    if (h_distribution_adv == "equation") {
      h[location_mut_adv_row] <- 1 / ((1 / h_intercept_adv) - (-1 * h_rate_adv * abs(s[location_mut_adv_row])))
    }
  }
  
  ###### INITIAL FREQUENCY (q) #####
  ## Assign initial allele frequency values based on locus type
  if (chunk_neutral_loci > 0) {
    q[location_neutral_row] <- q_neutral
  }
  if (real_loc == TRUE) {
    q[location_real_row] <- q_neutral
  }
  if (loci_deleterious > 0) {
    if (q_distribution_del == "equal") {
      q[location_deleterious_row] <- q_del
    }
    if (q_distribution_del == "equation") {
      a <- abs(s[location_deleterious_row]) * (1 - (2 * h[location_deleterious_row]))
      b <- (h[location_deleterious_row] * abs(s[location_deleterious_row])) * (1 + q_equation_del)
      c <- rep.int(-(q_equation_del), times = loci_deleterious)
      df_q <- as.data.frame(cbind(a, b, c))
      ## Solve the quadratic equation for equilibrium frequency based on Crow & Kimura (page 260)
      q[location_deleterious_row] <- mapply(q_equilibrium, a = df_q$a, b = df_q$b, c = df_q$c, USE.NAMES = F)
    }
  }
  if (loci_advantageous > 0) {
    if (q_distribution_adv == "equal") {
      q[location_advantageous_row] <- q_adv
    }
    if (q_distribution_adv == "equation") {
      a <- abs(s[location_advantageous_row]) * (1 - (2 * h[location_advantageous_row]))
      b <- (h[location_advantageous_row] * abs(s[location_advantageous_row])) * (1 + q_equation_adv)
      c <- rep.int(-(q_equation_adv), times = loci_advantageous)
      df_q <- as.data.frame(cbind(a, b, c))
      q[location_advantageous_row] <- mapply(q_equilibrium, a = df_q$a, b = df_q$b, c = df_q$c, USE.NAMES = F)
    }
  }
  if (loci_mut_neu > 0) {
    q[location_mut_neu_row] <- 0
  }
  if (loci_mut_del > 0) {
    q[location_mut_del_row] <- 0
  }
  if (loci_mut_adv > 0) {
    q[location_mut_adv_row] <- 0 
  }
  
  ## Create the reference table data frame using the computed values
  reference <- as.data.frame(matrix(nrow = total_loci))
  reference$q <- q
  reference$h <- h
  reference$s <- s
  reference$c <- recombination_map[1:total_loci, "c"]
  reference$loc_bp <- recombination_map[1:total_loci, "location_loci_bp"]
  reference$loc_cM <- recombination_map[1:total_loci, "accum"]
  reference$chr_name <- chromosome_name
  reference$type <- NA
  reference <- reference[, -1]  # remove the first temporary column
  
  ## Label each row in the reference table based on the type of locus
  if (real_loc == TRUE) {
    reference[location_real_row, "type"] <- "real"
  }
  if (chunk_neutral_loci > 0) {
    reference[location_neutral_row, "type"] <- "neutral"
  }
  if (loci_deleterious > 0) {
    reference[location_deleterious_row, "type"] <- "deleterious"
  }
  if (loci_advantageous > 0) {
    reference[location_advantageous_row, "type"] <- "advantageous"
  }
  if (loci_mut_neu > 0) {
    reference[location_mut_neu_row, "type"] <- "mutation_neu"
  }
  if (loci_mut_del > 0) {
    reference[location_mut_del_row, "type"] <- "mutation_del"
  }
  if (loci_mut_adv > 0) {
    reference[location_mut_adv_row, "type"] <- "mutation_adv"
  }
  
  ## Adjust values: cap q at 0.5 for non-equal distributions and adjust extreme s values
  if (q_distribution_del != "equal") {
    q_more_than_point5 <- as.numeric(row.names(reference[reference$q > 0.5, ]))
    reference[q_more_than_point5, "q"] <- 0.5
  }
  if (s_distribution_del != "equal") {
    s_more_than_one <- as.numeric(row.names(reference[reference$s > 1, ]))
    reference[s_more_than_one, "s"] <- 0.99
    s_less_minus_one <- as.numeric(row.names(reference[reference$s < -0.5, ]))
    reference[s_less_minus_one, "s"] <- -0.5
  }
  
  ## Prepare the result list containing the reference table and the variable values table
  ref_res <- list(reference, ref_vars)
  names(ref_res) <- c("reference", "ref_vars")
  
  ##### END OF FUNCTION #####
  
  ## Flag the end of the function execution
  if (verbose >= 1) {
    message(report("Completed:", funname, "\n"))
  }
  
  ## Return the results invisibly
  return(invisible(ref_res))
}
