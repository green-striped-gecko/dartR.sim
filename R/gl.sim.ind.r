#' @name gl.sim.ind
#' @title Simulates individuals based on allele frequencies
#' @family simulation functions
#' @description
#' This function simulates individuals based on the allele frequencies of a
#' genlight object. The output is a genlight object with the same loci as the
#' input genlight object.
#'
#' @param x Genlight object containing the SNP data [required].
#' @param n Number of individuals that should be simulated [default 50].
#' @param popname A population name for the simulated individuals
#' [default "pop1"].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' brief progress messages; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @details
#' For each locus, the frequency p of the alternative allele (genotype 2) is
#' estimated from x, ignoring missing data. Each simulated genotype is the
#' number of alternative alleles in two independent draws, so genotypes are
#' in Hardy-Weinberg proportions and loci are in linkage equilibrium. Loci
#' with no calls in x get missing genotypes (NA).
#'
#' The simulated object keeps the locus names, alleles, positions,
#' chromosomes and locus metrics of x; the metrics flags are reset, so
#' recalculable metrics are recalculated when needed.
#'
#' The simulated object has the same flavour as x: if the genotypes of x are
#' held in a file-backed matrix (@fbm, see \code{dartR.base::gl.gen2fbm}),
#' the simulated genotypes are written to a new file-backed matrix, otherwise
#' they are stored in @gen. In the file-backed matrix missing genotypes are
#' stored as the raw code 3, which is read back as NA.
#'
#' The function can be used to simulate populations for sampling designs or
#' for power analysis. The example below explores drift by simulating several
#' generations, each from the allele frequencies of the previous one.
#' @return A genlight object with n individuals, file-backed (@fbm) if x is
#' file-backed and gen-backed (@gen) otherwise.
#' @author Author(s): Bernd Gruber. Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' glsim <- gl.sim.ind(testset.gl, n=10, popname='sims')
#' glsim
#' ###Simulate drift over 10 generation
#' # assuming a bottleneck of only 10 individuals
#' # [ignoring effect of mating and mutation]
#' # Simulate 20 individuals with no structure and 50 SNP loci
#' founder <- glSim(n.ind = 20, n.snp.nonstruc = 50, ploidy=2)
#' #number of fixed loci in the first generation
#' res <- sum(colMeans(as.matrix(founder), na.rm=TRUE) %%2 ==0)
#' simgl <- founder
#' #49 generations of only 10 individuals
#' for (i in 2:50) {
#'    simgl <- gl.sim.ind(simgl, n=10, popname='sims', verbose = 0)
#'    res[i]<- sum(colMeans(as.matrix(simgl), na.rm=TRUE) %%2 ==0)
#' }
#' plot(1:50, res, type='b', xlab='generation', ylab='# fixed loci')
#' @export

gl.sim.ind <- function(x,
                       n = 50,
                       popname = "pop1",
                       verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  funname <- utils.flag.start(func = funname,
                   verbose = verbose)
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)
  
  # FUNCTION SPECIFIC ERROR CHECKING
  if (!is.numeric(n) || length(n) != 1 || is.na(n) || n < 1 || n %% 1 != 0) {
    stop(error("  n must be a single whole number of at least 1\n"))
  }
  
  # DO THE JOB
  # Frequency of the alternative allele (genotype 2) of each locus
  # (for an FBM-backed x without decoding the genotype matrix)
  alf <- utils.sim.alf(x)
  no_calls <- is.nan(alf)
  if (any(no_calls) && verbose >= 2) {
    cat(warn("  Warning:", sum(no_calls), "loci have no calls; their",
             "simulated genotypes are NA\n"))
  }
  
  # Each genotype is the number of alternative alleles in two independent
  # draws, i.e. Hardy-Weinberg proportions (1-p)^2, 2p(1-p), p^2
  simind <- matrix(NA_integer_, nrow = n, ncol = length(alf))
  simind[, !no_calls] <- rbinom(n * sum(!no_calls), size = 2,
                                prob = rep(alf[!no_calls], each = n))
  
  # Same flavour as x: FBM-backed if x is, gen-backed otherwise
  glsim <- utils.sim.dartR(simind,
                           fbm = utils.sim.is.fbm(x),
                           ind.names = as.character(1:n),
                           loc.names = locNames(x),
                           pop = rep(popname, n))
  glsim@loc.all <- x@loc.all
  if (!is.null(x@position)) {
    position(glsim) <- position(x)
  }
  if (!is.null(x@chromosome)) {
    chromosome(glsim) <- chromosome(x)
  }
  
  # Locus metadata come from x; the flags are reset because metrics such as
  # call rate describe x, not the simulated individuals
  if (!is.null(x@other$loc.metrics)) {
    glsim@other$loc.metrics <- x@other$loc.metrics
  }
  glsim <- utils.reset.flags(glsim, verbose = 0)
  
  # ADD TO HISTORY
  glsim@other$history <- list(utils.history.call(match.call()))
  
  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  return(glsim)
}
