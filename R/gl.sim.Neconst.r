#' @name gl.sim.Neconst
#' @title Simulates a population at mutation-drift equilibrium with constant Ne
#' @family simulation functions
#' @description This function simulates SNP genotypes of a population of
#' constant effective size at mutation-drift equilibrium, using Wright's beta
#' distribution for allele frequencies.
#' @param ninds Number of individuals in the population, which is also its
#' effective population size (Ne); a whole number of at least 2 [required].
#' @param nlocs Number of loci (SNPs) to simulate, a whole number of at least 1
#' [required].
#' @param mutation_rate Mutation rate per locus per generation, a number
#' greater than 0 [default 1e-8].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' brief progress messages; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @return A genlight object representing the simulated population.
#' @details At mutation-drift equilibrium under symmetric mutation between two
#' alleles, the allele frequency p of a locus follows Wright's beta
#' distribution, with density proportional to (p(1 - p))^(theta - 1), where
#' theta = 4 * ninds * mutation_rate.
#'
#' The population has 2 * ninds gene copies. For each locus, the number of
#' copies k of one allele is drawn from 1 to 2 * ninds - 1 with probability
#' proportional to (k/2N (1 - k/2N))^(theta - 1), i.e. conditioned on the locus
#' being polymorphic, and the k copies are placed at random among the gene
#' copies of the individuals. Every locus is therefore a SNP.
#'
#' When theta is much smaller than 1 (e.g. mutation_rate = 1e-8 for any
#' realistic ninds), the spectrum is the neutral site frequency spectrum,
#' proportional to 1/k + 1/(2N - k), and does not depend on mutation_rate.
#' mutation_rate changes the output only when theta approaches 1 or more, when
#' intermediate frequencies become more common.
#'
#' Individuals are named "1", "2", ..., loci "Loc1", "Loc2", ..., all
#' individuals belong to population "pop1", and alleles are set to the
#' placeholder "A/C".
#' @author Author(s): Bernd Gruber. Custodian: Bernd Gruber -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @importFrom stats dbeta rhyper
#' @importFrom dartR.popgen gl.sfs
#' @export
#' @examples
#' # Simulate a population with 50 individuals and 4000 loci
#' gg <- gl.sim.Neconst(ninds = 50, nlocs = 4000, mutation_rate = 1e-8, verbose = 0)
#' dartR.popgen::gl.sfs(gg)

gl.sim.Neconst <- function(ninds,
                           nlocs,
                           mutation_rate = 1e-8,
                           verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   verbose = verbose)

  # FUNCTION SPECIFIC ERROR CHECKING
  is_whole <- function(v, min) {
    is.numeric(v) && length(v) == 1 && !is.na(v) && v >= min && v %% 1 == 0
  }
  if (!is_whole(ninds, 2)) {
    stop(error("  ninds must be a whole number of at least 2\n"))
  }
  if (!is_whole(nlocs, 1)) {
    stop(error("  nlocs must be a whole number of at least 1\n"))
  }
  if (!is.numeric(mutation_rate) || length(mutation_rate) != 1 ||
      is.na(mutation_rate) || mutation_rate <= 0) {
    stop(error("  mutation_rate must be a single number greater than 0\n"))
  }

  # DO THE JOB
  n_copies <- 2 * ninds
  theta <- 4 * ninds * mutation_rate

  # Number of copies of one allele at each locus, from Wright's beta
  # distribution on the 2N gene copies, conditioned on polymorphism
  k <- seq_len(n_copies - 1)
  p <- k / n_copies
  w <- dbeta(p, theta, theta)
  count <- sample(k, nlocs, replace = TRUE, prob = w)

  # Place the copies at random among the gene copies: each individual in turn
  # draws its two copies without replacement from the copies not yet assigned
  genmat <- matrix(0L, nrow = ninds, ncol = nlocs)
  left_allele <- count
  left_total <- n_copies
  for (i in seq_len(ninds)) {
    g <- rhyper(nlocs, left_allele, left_total - left_allele, 2)
    genmat[i, ] <- g
    left_allele <- left_allele - g
    left_total <- left_total - 2
  }
  
  loc_names <- paste0("Loc", seq_len(nlocs))
  inds <- new(
    "dartR",
    gen = genmat,
    ind.names = as.character(seq_len(ninds)),
    loc.names = loc_names,
    ploidy = rep(2, ninds)
  )
  # Alleles set before the compliance check, so that SNP data whose genotypes
  # happen to be all 0 or 1 are not taken for SilicoDArT
  inds@loc.all <- rep("A/C", nlocs)
  inds@other$loc.metrics <- data.frame(AlleleID = loc_names)
  inds <- gl.compliance.check(inds, verbose = 0)
  inds@other$ind.metrics$pop <- as.character(pop(inds))

  if (verbose >= 3) {
    cat(report("  Simulated", ninds, "individuals at", nlocs,
               "polymorphic loci (theta =", signif(theta, 3), ")\n"))
  }

  # ADD TO HISTORY
  inds@other$history <- list(match.call())

  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  return(inds)
}
