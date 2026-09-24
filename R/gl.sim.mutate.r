#' @name gl.sim.mutate
#' @title Simulates mutations within a genlight object
#' @family simulation functions
#' @description
#' This function is intended to be used within the simulation framework of
#' dartR. It adds the ability to add a constant mutation rate across all loci.
#' Only works for biallelic SNP data.
#'
#' @details
#' Each of the nInd * nLoc * 2 allele copies mutates with probability
#' mut.rate, so the number of mutation events is drawn from a binomial
#' distribution, and each event is placed at a random individual and locus.
#' A mutation changes one allele: a homozygote (0 or 2) becomes a heterozygote
#' (1), and a heterozygote becomes either homozygote with equal probability.
#'
#' Events that fall on a missing genotype are dropped, so the realised rate is
#' lower than mut.rate by the proportion of missing data. In principle
#' 'double mutations' at the same locus of the same individual can occur, but
#' should be rare.
#'
#' If any genotype changes, the locus metrics flags are reset so that metrics
#' are recalculated by the functions that use them.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param mut.rate Constant mutation rate per allele copy, a number between 0
#' and 1 [default 1e-6].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' brief progress messages; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @return Returns a genlight object with the applied mutations
#' @export
#' @author Author(s): Bernd Gruber. Custodian: Bernd Gruber -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' b2 <- gl.sim.mutate(bandicoot.gl,mut.rate=1e-4 )
#' #check the mutations that have occurred
#' table(as.matrix(bandicoot.gl), as.matrix(b2))

gl.sim.mutate <- function(x,
                          mut.rate = 1e-06,
                          verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   verbose = verbose)

  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, accept = "SNP", verbose = 0)

  # FUNCTION SPECIFIC ERROR CHECKING
  if (!is.numeric(mut.rate) || length(mut.rate) != 1 || is.na(mut.rate) ||
      mut.rate < 0 || mut.rate > 1) {
    stop(error("  mut.rate must be a single number between 0 and 1\n"))
  }

  # DO THE JOB
  nm <- rbinom(1, nInd(x) * nLoc(x) * 2, mut.rate)
  n_applied <- 0
  for (ii in seq_len(nm)) {
    ri <- sample(1:nInd(x), 1)
    rl <- sample(1:nLoc(x), 1)
    # Genotypes of the mutated individual only (not the whole matrix)
    xx <- as.integer(x@gen[[ri]])
    cs <- xx[rl]
    if (!is.na(cs)) {
      if (!cs %% 2)
        nv <- 1
      else
        nv <- sample(c(0, 2), 1)
      xx[rl] <- nv
      x@gen[[ri]] <- new("SNPbin", snp = xx, ploidy = 2L)
      n_applied <- n_applied + 1
    }
  }

  if (n_applied > 0) {
    x <- utils.reset.flags(x, verbose = 0)
  }

  if (verbose >= 3) {
    cat(report("  Applied", n_applied, "mutations;", nm - n_applied,
               "fell on missing genotypes\n"))
  }

  # ADD TO HISTORY
  nh <- length(x@other$history)
  x@other$history[[nh + 1]] <- match.call()

  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  return(x)
}
