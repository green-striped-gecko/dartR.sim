#' @name gl.sim.offspring
#' @title Simulates offspring based on alleles provided by parents
#' @family simulation functions
#' @description
#' This function takes a population (or a single individual) of fathers and
#' of mothers (provided as genlight objects) and simulates offspring by
#' Mendelian inheritance. It can be used to simulate population dynamics and
#' check their effect on allele frequencies and number of alleles, or to
#' simulate the relatedness of siblings and compare it with the relatedness
#' found in a population.
#'
#' @param fathers Genlight object of potential fathers [required].
#' @param mothers Genlight object of mothers [required].
#' @param noffpermother Number of offspring per mother [required].
#' @param sexratio The sex ratio of simulated offspring 
#' (females / (females + males); 1 equals 100 percent females) [default 0.5].
#' @param popname Population name of the returned genlight object 
#' [default "offspring"].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' brief progress messages; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @details
#' Each mother has exactly noffpermother offspring. The father of each
#' offspring is drawn at random from fathers, so offspring of the same mother
#' are full siblings only when a single father is given; otherwise most are
#' half siblings.
#' 
#' Each parent passes on one of its two alleles at every locus with
#' probability 1/2, independently across loci (no linkage). Missing
#' genotypes in a parent give missing genotypes in its offspring.
#' 
#' fathers and mothers must be SNP data with the same loci in the same
#' order. The offspring keep the locus metadata of the mothers (the metrics
#' flags are reset), and @other$ind.metrics records each offspring's sex,
#' mother and father (their individual names).
#' @return A genlight object with nInd(mothers) * noffpermother individuals.
#' @importFrom stats runif
#' @export
#' @author Author(s): Bernd Gruber. Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' #Simulate 10 potential fathers
#' gl.fathers <- glSim(10, 20, ploidy=2)
#' #Simulate 10 potential mothers
#' gl.mothers <- glSim(10, 20, ploidy=2)
#' res <- gl.sim.offspring(gl.fathers, gl.mothers, 2, sexratio=0.5)

gl.sim.offspring <- function(fathers,
                             mothers,
                             noffpermother,
                             sexratio = 0.5, 
                             popname = "offspring",
                             verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   verbose = verbose)
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(fathers, accept = "SNP", verbose = 0)
  datatype <- utils.check.datatype(mothers, accept = "SNP", verbose = 0)
  
  # FUNCTION SPECIFIC ERROR CHECKING
  if (!identical(locNames(fathers), locNames(mothers))) {
    stop(error("  fathers and mothers must have the same loci in the same",
               "order\n"))
  }
  if (!is.numeric(noffpermother) || length(noffpermother) != 1 ||
      is.na(noffpermother) || noffpermother < 1 || noffpermother %% 1 != 0) {
    stop(error("  noffpermother must be a single whole number of at least",
               "1\n"))
  }
  if (!is.numeric(sexratio) || length(sexratio) != 1 || is.na(sexratio) ||
      sexratio < 0 || sexratio > 1) {
    stop(error("  sexratio must be a single number between 0 and 1\n"))
  }
  
  # DO THE JOB
  # Each mother has exactly noffpermother offspring; each offspring has a
  # father drawn at random
  noff <- nInd(mothers) * noffpermother
  mother <- rep(seq_len(nInd(mothers)), each = noffpermother)
  father <- sample(seq_len(nInd(fathers)), noff, replace = TRUE)
  
  if (verbose >= 2 &&
      (anyNA(as.matrix(mothers)) || anyNA(as.matrix(fathers)))) {
    cat(warn("  Warning: You have missing data in your genlight object.",
             "This most likely will cause unwanted structure in your",
             "offspring. Best to remove or impute missing values.\n"))
  }
  
  # Gamete of each parent at each locus: a heterozygote passes on either
  # allele with probability 1/2 (coded 0 or 2 so that the offspring genotype
  # is the mean of the two gametes); each heterozygous genotype gets its own
  # draw, so loci are inherited independently
  gamete <- function(parents, which_parent) {
    g <- as.matrix(parents)[which_parent, , drop = FALSE]
    het <- which(g == 1)
    g[het] <- sample(c(0, 2), length(het), replace = TRUE)
    return(g)
  }
  offmat <- (gamete(mothers, mother) + gamete(fathers, father)) / 2
  
  gl2 <-
    new(
      "dartR",
      gen = offmat,
      ploidy = 2,
      ind.names = paste0("Po_", 1:noff),
      loc.names = locNames(mothers),
      pop = rep(popname, noff)
    )
  if (length(mothers@loc.all) == nLoc(mothers)) {
    gl2@loc.all <- mothers@loc.all
  }
  if (!is.null(mothers@position)) {
    position(gl2) <- position(mothers)
  }
  if (!is.null(mothers@chromosome)) {
    chromosome(gl2) <- chromosome(mothers)
  }
  if (!is.null(mothers@other$loc.metrics)) {
    gl2@other$loc.metrics <- mothers@other$loc.metrics
  }
  gl2 <- utils.reset.flags(gl2, verbose = 0)
  
  # Sex and parents of each offspring
  # (parents are named by indNames, or by position if they have no names)
  parent_names <- function(parents) {
    nm <- indNames(parents)
    if (is.null(nm)) as.character(seq_len(nInd(parents))) else nm
  }
  sr <- factor(ifelse(runif(noff) < sexratio, "female", "male"),
               levels = c("female", "male"))
  gl2@other$sex <- sr
  gl2@other$ind.metrics <- data.frame(id = indNames(gl2),
                                      pop = popname,
                                      sex = sr,
                                      mother = parent_names(mothers)[mother],
                                      father = parent_names(fathers)[father])
  
  # ADD TO HISTORY
  gl2@other$history <- list(match.call())
  
  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  return(gl2)
}
