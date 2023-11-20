#' @name gl.sim.offspring
#' @title Simulates offspring based on alleles provided by parents
#' @description
#' This takes a population (or a single individual) of fathers (provided as a
#' genlight object) and mother(s) and simulates offspring based on 'random'
#'  mating. It can be used to simulate population dynamics and check the effect
#'  of those dynamics and allele frequencies, number of alleles. Another
#'  application is to simulate relatedness of siblings and compare it to actual
#'  relatedness found in the population to determine kinship.
#'
#' @param fathers Genlight object of potential fathers [required].
#' @param mothers Genlight object of potential mothers simulated [required].
#' @param noffpermother Number of offspring per mother [required].
#' @param sexratio The sex ratio of simulated offspring 
#' (females / females +males, 1 equals 100 percent females) [default 0.5.].
#' @param popname population name of the returned genlight object 
#' [default offspring]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @return A genlight object with n individuals.
#' @importFrom stats runif
#' @export
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
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
  
    noff <- nInd(mothers) * noffpermother
    mother <- sample(1:nInd(mothers), noff, replace = T)
    father <- sample(1:nInd(fathers), noff, replace = T)
    
    
    if (sum( c(is.na(as.matrix(mothers)), is.na(as.matrix(fathers))))>0 & verbose >= 2) 
      message(warn("Warning: You have missing data in your genlight object.\nThis most likely will cause unwanted structure in you offspring.\nBest to remove or impute missing values."))
    mmat <- as.matrix(mothers)[mother, ]
    mhet <- sum(mmat == 1, na.rm=T)
    if (!is.na(mhet)) {
        mother.half <-
            ifelse(mmat == 1, sample(c(0, 2), mhet, replace = T), mmat)
    }
    
    fmat <- as.matrix(fathers)[father, ]
    fhet <- sum(fmat == 1, na.rm = T)
    if (!is.na(fhet)) {
        father.half <-
            ifelse(fmat == 1, sample(c(0, 2), fhet, replace = T), fmat)
    } 
    
    offmat <- (mother.half + father.half) / 2
    gl2 <-
        new(
            "genlight",
            gen = offmat,
            ind.names = paste0("Po_", 1:noff),
            loc.names = locNames(mothers),
            ploidy = rep(2, nrow(offmat))
        )
    pop(gl2)<- rep(popname, nInd(gl2))
    
    # set sex ratio
    sr <-
        factor(ifelse(runif(nInd(gl2)) < sexratio, "female", "male"))
    gl2@other$sex <- sr
    return(gl2)
}
