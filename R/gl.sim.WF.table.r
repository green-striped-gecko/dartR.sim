#' @name gl.sim.WF.table
#' @title Creates the reference table for running gl.sim.WF.run
#' @family simulation functions
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
#' @param file_var Path of the variables file 'ref_variables.csv' (see details) 
#' [required if interactive_vars = FALSE].
#' @param x Genlight object containing the SNP data to extract
#' values for some simulation variables (see details) [default NULL].
#' @param file_targets_sel Path of the file with the targets for selection (see 
#' details) [default NULL].
#' @param file_r_map Path of the file with the recombination map (see details)
#' [default NULL].
#' @param interactive_vars Run a shiny app to input interactively the values of
#'  simulation variables [default TRUE].
#' @param seed Set the seed for the simulations. This calls set.seed(), so it
#' also sets the random number stream of the R session for any code run
#' afterwards [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' brief progress messages; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @param ... Any simulation variable of 'ref_variables.csv' and its value,
#' e.g. chunk_number = 20. The value replaces the one in the csv file or the
#' Shiny app. Names that are not simulation variables stop the function with
#' an error.
#' @details
#' Values for the variables to create the reference table can be submitted into 
#' the function interactively through a Shiny app if interactive_vars = TRUE. 
#' Optionally, if interactive_vars = FALSE, values for variables can be 
#' submitted by using the csv file 'ref_variables.csv' which can be found by 
#' typing in the R console:
#'  system.file('extdata', 'ref_variables.csv', package = 'dartR.sim').
#'  
#' The values of the variables can be modified using the third column ("value") 
#' of this file. 
#' 
#' If a genlight object is used as input (real_loc = TRUE or real_freq = TRUE),
#' this function uses the slots x@position and x@chromosome. x@chromosome must
#' contain the value of the variable chromosome_name.
#' 
#' The recombination map file needs the columns Chr (chromosome name), from and
#' to (start and end of each interval in bp) and cM (centiMorgans in the
#' interval). Intervals can have any size; intervals with cM = NA are treated
#' as not recombining. The targets of selection file needs the columns
#' chr_name, start and end (in bp) and targets (number of targets of selection
#' in the region). Examples of both files can be found by typing in the R
#' console:
#' \itemize{ 
#' \item system.file('extdata', 'fly_recom_map.csv', package = 'dartR.sim')
#' \item system.file('extdata', 'fly_targets_of_selection.csv', package = 'dartR.sim')
#' }
#' 
#' Values drawn from distributions or equations are capped: q at 0.5,
#' deleterious s at 0.99 and advantageous s at -0.5. Each class is capped only
#' when its own distribution setting is not "equal". The number of capped loci
#' is reported when verbose >= 1.
#' 
#' To show further information of the variables in interactive mode, it might be
#'  necessary to call first: 'library(shinyBS)' for the information to be 
#'  displayed.
#' @return Returns a list with the reference table used as input for the function
#'  \code{\link{gl.sim.WF.run}} and a table with the values variables used to 
#'  create the reference table.
#' @author Author(s): Luis Mijangos. Custodian: Luis Mijangos -- Post to
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
                   verbose = verbose)

  ## CHECK INPUTS
  if (interactive_vars == FALSE &&
      (missing(file_var) || !file.exists(file_var))) {
    stop(error("  When interactive_vars = FALSE, file_var must be the path to",
               "an existing 'ref_variables.csv' file\n"))
  }
  if (!is.null(x) && !is(x, "genlight")) {
    stop(error("  x must be a genlight object\n"))
  }

  ## LOADING VARIABLES
  ## Variables come from the Shiny app or from the CSV file; either way
  ## ref_vars is a two-column table (variable, value) of character values
  if (interactive_vars == TRUE) {
    ref_vars <- interactive_reference()
  } else {
    ref_vars <- suppressWarnings(read.csv(file_var))
    ref_vars <- ref_vars[, c("variable", "value")]
  }
  ref_vars <- ref_vars[order(ref_vars$variable), ]

  ## Values passed through ... replace the values in ref_vars. Names that are
  ## not simulation variables are refused, so a misspelt name cannot be
  ## silently ignored
  input_list <- list(...)
  if (length(input_list) > 0) {
    unknown <- setdiff(names(input_list), ref_vars$variable)
    if (is.null(names(input_list)) || any(names(input_list) == "") ||
        length(unknown) > 0) {
      stop(error("  Arguments passed through ... must be named simulation",
                 "variables. Unknown:", paste(unknown, collapse = ", "),
                 "\n"))
    }
    for (var in names(input_list)) {
      ref_vars[ref_vars$variable == var, "value"] <-
        as.character(input_list[[var]])
    }
  }

  ## Create one R variable per simulation variable in this environment
  list2env(utils.wf.ref.values(ref_vars), envir = environment())

  ##### LOADING INFORMATION #####
  ## RECOMBINATION MAP: load recombination data if provided
  if (!is.null(file_r_map)) {
    map <- read.csv(file_r_map, check.names = FALSE)
    ## Header names can carry stray spaces (e.g. "to " in fly_recom_map.csv)
    colnames(map) <- trimws(colnames(map))
    map$Chr <- as.character(map$Chr)

    ## Check that the chromosome name is present in the recombination map
    if (!chromosome_name %in% map$Chr) {
      stop(error("  Chromosome name is not in the recombination map file\n"))
    }
    ## Subset the map for the given chromosome and convert cM to Morgans
    map <- map[which(map$Chr == chromosome_name), ]
    map <- map[order(map$from), ]
    ## Intervals without an estimate are treated as not recombining
    map$cM[is.na(map$cM)] <- 0
    map$cM <- map$cM / 100
  }

  ## TARGETS OF SELECTION: load targets file if provided
  targets_temp <- NULL
  if (!is.null(file_targets_sel)) {
    targets_temp <- read.csv(file_targets_sel)
    targets_temp$chr_name <- as.character(targets_temp$chr_name)

    ## Check that the chromosome name exists in the targets file
    if (!chromosome_name %in% targets_temp$chr_name) {
      stop(error("  Chromosome name is not in the targets of selection file\n"))
    }
    targets_temp <- targets_temp[which(targets_temp$chr_name == chromosome_name),]
  }

  ## REAL DATA: ensure that if real location or frequency info is needed, the dataset is provided
  if ((real_loc == TRUE | real_freq == TRUE) && is.null(x)) {
    stop(error("  The real dataset to extract information is missing\n"))
  }
  location_real_temp <- NULL
  if (real_loc == TRUE & !is.null(x)) {
    if (!chromosome_name %in% x@chromosome) {
      stop(error("  Chromosome name is not in the genlight object\n"))
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
    ## Neutral loci are spread over the whole mapped chromosome, whatever
    ## the size of the map intervals
    chunk_bp <- chr_length / chunk_number
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
        ## Spread the loci evenly and give the remainder to random chunks,
        ## so the total equals loci_deleterious
        del$targets <- floor(loci_deleterious / chunk_number)
        extra <- loci_deleterious %% chunk_number
        if (extra > 0) {
          row_targets <- sample(1:chunk_number, size = extra)
          del[row_targets, "targets"] <- del[row_targets, "targets"] + 1
        }
      }
      del$distance <-  del$end - del$start
    }
    sample_resolution <- max(1, round(mean(del$distance) / max(del$targets) / 10))
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
      ## Sample by index: sample() on a single number n would draw from 1:n
      location_deleterious_temp <- location_deleterious_temp[
        sample.int(length(location_deleterious_temp), size = del[i, "targets"])]
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
    sample_resolution <- max(1, round(mean(mutations$distance) / max(mutations$targets) / 20))
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
      location_mutations_temp <- location_mutations_temp[
        sample.int(length(location_mutations_temp), size = mutations[i, "targets"])]
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
    stop(error("  Number of loci should be more than the number of genome chunks\n"))
  }
  
  total_loci <- length(location_loci_bp)
  
  ##### RECOMBINATION MAP #####
  ## Without a map file, the chromosome is split into chunk_number intervals
  ## of chunk_bp bp, each with chunk_cM cM. This is built here, after chunk_bp
  ## may have been rescaled to the real data or the targets file
  if (is.null(file_r_map)) {
    map <- data.frame(from = (0:(chunk_number - 1)) * chunk_bp + 1,
                      to = (1:chunk_number) * chunk_bp,
                      cM = chunk_cM / 100)
  }
  ## Each map interval has its own length and midpoint, taken from its from
  ## and to columns, so maps need not have chunk_bp-sized intervals
  map_length <- map$to - map$from + 1
  map_midpoint <- map$to - map_length / 2
  ## Each locus takes the rate (Morgans per bp) of the interval whose midpoint
  ## is the closest one at or below it
  recombination_temp <- findInterval(location_loci_bp, map_midpoint)
  
  ## Correct intervals that fall before the first midpoint
  recombination_temp[recombination_temp == 0] <- 1
  recombination_2 <- map[recombination_temp, "cM"] /
    map_length[recombination_temp]
  recombination_map <- as.data.frame(cbind(location_loci_bp, recombination_2))
  recombination_map$c <- NA
  
  ## Recombination between each locus and the next one (except the last)
  for (target_row in 1:(nrow(recombination_map) - 1)) {
    recombination_map[target_row, "c"] <- (recombination_map[target_row + 1, "location_loci_bp"] - 
                                              recombination_map[target_row, "location_loci_bp"]) * 
                                             recombination_map[target_row, "recombination_2"]
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
  if (real_loc == TRUE | real_freq == TRUE) {
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
  if (real_loc == TRUE | real_freq == TRUE) {
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
  if (real_loc == TRUE | real_freq == TRUE) {
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
  if (real_loc == TRUE | real_freq == TRUE) {
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
  
  ## Cap values drawn from distributions or equations: q at 0.5, deleterious
  ## s at 0.99 and advantageous s at -0.5. Each class is capped only when its
  ## own setting is not "equal", so values set by the user are kept
  del_rows <- which(reference$type %in% c("deleterious", "mutation_del"))
  adv_rows <- which(reference$type %in% c("advantageous", "mutation_adv"))
  cap_q <- c(if (q_distribution_del != "equal") del_rows,
             if (q_distribution_adv != "equal") adv_rows)
  cap_q <- cap_q[reference$q[cap_q] > 0.5]
  reference[cap_q, "q"] <- 0.5
  cap_s_del <- if (s_distribution_del != "equal") del_rows
  cap_s_del <- cap_s_del[reference$s[cap_s_del] > 1]
  reference[cap_s_del, "s"] <- 0.99
  cap_s_adv <- if (s_distribution_adv != "equal") adv_rows
  cap_s_adv <- cap_s_adv[reference$s[cap_s_adv] < -0.5]
  reference[cap_s_adv, "s"] <- -0.5
  n_capped <- length(cap_q) + length(cap_s_del) + length(cap_s_adv)
  if (verbose >= 1 && n_capped > 0) {
    cat(warn("  Values capped: q > 0.5 set to 0.5 in", length(cap_q),
             "loci; deleterious s > 1 set to 0.99 in", length(cap_s_del),
             "loci; advantageous s < -0.5 set to -0.5 in", length(cap_s_adv),
             "loci\n"))
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
