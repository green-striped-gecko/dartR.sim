#' @name gl.report.nall
#' @title
#' Allele Count Reporting and Rarefaction Curve Simulation
#'
#' @description
#' The function reports allele counts for each population and simulates a
#' rarefaction curve. This curve demonstrates the proportion of the total 
#' allelic diversity captured as increasing numbers of individuals are 
#' subsampled from the combined dataset.
#'
#' @param x Name of the genlight/dartR object containing the SNP data.
#' The object needs to have no missing data as subsampling from missing data
#' is not possible. So we recommend to filter by callrate using a threshold
#' of 1 [required].
#' @param simlevels a vector that defines the different levels the combined
#' population should be subsampled [default seq(1,nInd(x),5)].
#' @param reps number of replicates per subsample level [default 10].
#' @param plot.colors.pop A color palette for population plots or a list with
#' as many colors as there are populations in the dataset 
#' [default gl.colors("dis")].
#' @param ncores number of cores to be used for parallel processing 
#' [default 10].
#' @param plot.display Specify if plot is to be produced [default TRUE].
#' @param plot.theme User specified theme [default theme_dartR()].
#' @param plot.dir Directory to save the plot RDS files [default as specified
#' by the global working directory or tempdir()]
#' @param plot.file Filename (minus extension) for the RDS plot file
#' [Required for plot save]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' 
#' @details
#' To assess the effect of sample size on allelic diversity, the function 
#' simulates subsampling at various levels (as specified by the `simlevels` 
#' argument) and replicates this process (`reps` times). The maximum number 
#' of alleles is determined from a pooled population of all individuals, and 
#' simulation results are normalized to this maximum value. Finally, the 
#' function aggregates the simulation outputs to yield the mean, minimum, 
#' and maximum proportions of allele retention for each sample size, and 
#' produces a plot of the rarefaction curve and the allele counts for each 
#' population.
#'
#' @return 
#' A list with three elements:
#' \itemize{ 
#' \item 1. A data frame with the number of simulated alleles for each simlevel
#' (mean, min and mix)
#' \item 2. A data frame with the number of alleles of each subpopulation
#' \item 3. A ggplot object with the rarefaction curve of alleles for all
#' individuals combined
#' }
#'
#' @importFrom ggrepel geom_text_repel
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel makeCluster stopCluster
#' @importFrom dplyr group_by summarise
#' @importFrom doParallel registerDoParallel
#' @export
#' @author Bernd Gruber (bernd.gruber@@canberra.edu.au)
#' @examples
#' \dontrun{
#' dummy <- gl.report.nall(possums.gl[c(1:5,31:35),], simlevels=seq(1,10,3),
#' reps=5, ncores=2)
#' }

########################## packages needed
#library(parallel)
#library(foreach)
#library(ggrepel)
#library(dartRverse)  #needs dartR.base and dartR.sim

gl.report.nall <- function(x,
                           simlevels = seq(1, nInd(x), 5),
                           reps = 10,
                           plot.colors.pop = gl.colors("dis"),
                           ncores = 2,
                           plot.display = TRUE,
                           plot.theme = theme_dartR(),
                           plot.dir = NULL,
                           plot.file = NULL,
                           verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  if (verbose == 0) {
    plot.display <- FALSE
  }
  
  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir, verbose = 0)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "v.2023.3",
                   verbose = verbose)
  
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)
  
  # FUNCTION SPECIFIC ERROR CHECKING
  
  x <- gl.filter.allna(x,by.pop = TRUE)
  
  #  initial functions needed
  gl.report.nall.pop <- function(x, x2 = NULL) {
    if (!is.null(x2)) {
      pops <- list(pop1 = x, pop2 = x2)
    } else{
      pops <- seppop(x)
    }
    
    res <- data.frame(popname = NA,
                      Npop = NA,
                      N.all = NA)
    
    for (i in 1:length(pops)) {
      p <- colMeans(as.matrix(pops[[i]]), na.rm = TRUE) / 2
      res[i, 1] = names(pops)[i]
      res[i, 2] = nInd(pops[[i]])
      res[i, 3] = sum((p == 1) + (p == 0) + (2 * ((p > 0) & (p < 1))))
    }
    return(res)
  }
  
  ss <- function(sample) {
    foundersim <- gl.sim.ind(x, sample, popname = "foundersim")
    nasim <- gl.report.nall.pop(foundersim)
    res <- c(nasim$Npop, nasim$N.all)
    return(res)
  }
  
  Npop <- Nallsim <- mnall <- low <- high <- N.all <- popname <- ip <- NULL
  ##############################################################
  ### find maxnumber of alleles in the population
  onePop <- x
  pop(onePop) <- rep("A", nInd(onePop))
  maxnall <- gl.report.nall.pop(onePop)$N.all
  sims <- expand.grid(rep = 1:reps, Npop = simlevels)
  
  # run simulations in parallel
  
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  sims2 <- foreach::foreach(ip = 1:nrow(sims),
                            .combine = rbind,
                            .packages = "dartR.sim") %dopar%
    {
      simres <- ss(sims[ip, 2])
      return(simres)
    }
  
  parallel::stopCluster(cl)
  sims2 <- as.data.frame(sims2)
  colnames(sims2) <- c("Npop", "Nallsim")
  
  sims2$Nallsim <- sims2$Nallsim / maxnall
  
  df <- sims2 |>
    dplyr::group_by(Npop) |>
    dplyr::summarise(
      mnall = mean(Nallsim),
      low = min(Nallsim),
      high = max(Nallsim)
    )
  
  ###calculate the points if needed
  # number of alleles for each pop of interest....
  
  nall.pop <- gl.report.nall.pop(x)
  nall.pop$N.all <- nall.pop$N.all / maxnall
  
  # printing plots and reports assigning colors to populations
  if (is(plot.colors.pop, "function")) {
    colors_pops <- plot.colors.pop(length(levels(pop(x))))
  }
  
  if (!is(plot.colors.pop, "function")) {
    colors_pops <- plot.colors.pop
  }
  
  p1 <- ggplot(df, aes(x = Npop, y = mnall)) +
    geom_line() +
    geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2) +
    xlab("Number of sampled individuals") +
    ylab("Percentage of alleles retained \n compared to all individuals") +
    ylim(c(0.5, 1)) +
    geom_text_repel(data = nall.pop,
                    aes(x = Npop, y = N.all, label = popname),
                    hjust = -1) +
    geom_point(data = nall.pop,
               aes(x = Npop, y = N.all, colour = popname),
               size = 6,
               pch = 16) +
    guides(color = guide_legend(title = "Population")) +
    scale_color_manual(values = colors_pops) +
    plot.theme
  
  # PRINTING OUTPUTS
  
  if (plot.display) {
    print(p1)
  }
  
  if (!is.null(plot.file)) {
    tmp <- utils.plot.save(p1,
                           dir = plot.dir,
                           file = plot.file,
                           verbose = verbose)
  }
  
  # FLAG SCRIPT END
  
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  # RETURN OUTPUTS (df, points and plot as list)
  return(list(sim = df, points = nall.pop, p1 = p1))
  
}