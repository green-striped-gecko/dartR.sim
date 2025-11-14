#' @title Simulate a population with constant mutation rate
#' @description This function simulates a population with a constant mutation 
#' rate using a beta distribution for allele frequencies.
#' @param ninds Number of individuals in the population [required].
#' @param nlocs Number of loci in the population [required].
#' @param mutation_rate Mutation rate per generation (default is 1e-8) 
#' [default 1e-8].
# @param fbm If TRUE, the genlight object will be converted to a filebacked
# large matrix format, which is faster if the dataset is large
# [default FALSE, because still in a testing phase].
#' @param verbose Verbosity level (default is 0).
#' @return A genlight object representing the simulated population.
#' @details The function generates a genlight object with the specified number 
#' of individuals and loci, simulating allele frequencies based on a beta 
#' distribution. The mutation rate is used to calculate theta, which in turn is 
#' the parameter in the beta function.
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#' @importFrom stats dbeta
#' @importFrom dartR.popgen gl.sfs
#' @export
#' @examples
#' # Simulate a population with 50 individuals and 4000 loci
#' gg <- gl.sim.Neconst(ninds = 50, nlocs = 4000, mutation_rate = 1e-8, verbose = 0)
#' dartR.popgen::gl.sfs(gg)

gl.sim.Neconst <- function(ninds,
                           nlocs,
                           mutation_rate = 1e-8,
                           # fbm = FALSE,
                           verbose = 0) {
  # Function to calculate allele frequency distribution
  allele_frequency_distribution <- function(pop_size, mutation_rate, num_bins = 100) {
    # Theta parameter
    theta <- 4 * pop_size * mutation_rate
    
    # Generate the beta distribution
    x <- seq(0, 1, length.out = num_bins)
    y <- dbeta(x, theta, theta)  # Beta distribution for allele frequencies
    
    return(data.frame(frequency = x, density = y))
  }
  
  # Parameters
  pop_size <- ninds      # Effective population size
  mutation_rate <- mutation_rate # Mutation rate per generation
  num_bins <- 2 * (ninds * 2) + 1       # Number of frequency bins for plotting
  
  # Calculate allele frequency distribution
  freq_dist <- allele_frequency_distribution(pop_size, mutation_rate, num_bins)
  
  freq_dist$density[1] <- freq_dist$density[num_bins] <- 0
  
  # create a genlight object based on freq_dist
  allele_freq <- sample(freq_dist$frequency,
                        nlocs,
                        replace = TRUE,
                        prob = freq_dist$density)

  # Generate genotypes
  genmat <- matrix(0, nrow = ninds, ncol = nlocs)
  
  for (i in 1:nlocs) {
    genmat[, i] <- rbinom(ninds , 2, 1 - allele_freq[i])
  }
  
  L <- (nlocs / (ninds * 4 * mutation_rate))
  
  inds <- new(
    "dartR",
    gen = genmat,
    ind.names = 1:ninds,
    ploidy = rep(2, ninds)
  )

  inds <- gl.compliance.check(inds, verbose = 0)
  
  # if (fbm) inds <- gl.gen2fbm(inds)
  
  return(inds)
}
