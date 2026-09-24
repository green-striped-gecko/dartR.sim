#' @name gl.report.nall
#' @title
#' Report allelic retention and simulate a rarefaction curve
#'
#' @family simulation functions
#' @description
#' This function reports per-population allele counts and simulates a
#' rarefaction-style curve showing the proportion of the dataset's total allelic
#' diversity expected as progressively more individuals are sampled.
#'
#' @param x Name of the genlight/dartR object containing the SNP data. Loci
#' with no calls in any one population are removed before the analysis;
#' otherwise missing data are allowed, as allele frequencies are computed from
#' the called genotypes [required].
#' @param simlevels A vector of whole numbers of at least 1 with the sample
#' sizes (numbers of individuals) to simulate [default seq(1,nInd(x),5)].
#' @param reps Number of replicate simulated samples per sample size, a whole
#' number of at least 1 [default 10].
#' @param plot.colors.pop A color palette for population plots or a list with
#' as many colors as there are populations in the dataset 
#' [default gl.colors("dis")].
#' @param ncores Number of cores to be used for parallel processing, a whole
#' number of at least 1. With 1, the simulations run in the current R session
#' without starting a cluster [default 2].
#' @param plot.display Specify if plot is to be produced [default TRUE].
#' @param plot.theme A `ggplot2` theme object for styling the plot
#'  [default theme_dartR()].
#' @param plot.dir Directory to save the plot RDS files [default as specified
#' by the global working directory or tempdir()].
#' @param plot.file Filename (minus extension) for the RDS plot file
#' [Required for plot save].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' 
#' @details
#' The function estimates how sampling effort affects observed allelic diversity.
#' All individuals are pooled into a single group, and the allele frequencies of
#' the pool are computed. For each sample size in `simlevels`, `reps` samples of
#' that many individuals are simulated from the pooled frequencies with
#' \code{\link{gl.sim.ind}}, i.e. under Hardy-Weinberg and linkage equilibrium.
#' The curve is therefore the expectation for a panmictic population with the
#' allele frequencies of the whole dataset, not a subsampling of the real
#' individuals: structure between populations is not part of it, and alleles
#' that are rare in the pool can be missed even at the full sample size, so the
#' curve need not reach 1.
#'
#' The number of alleles counts 1 for a fixed locus and 2 for a polymorphic
#' one. The maximum is the count in the pooled real data; all simulation
#' outputs and per-population observations are divided by this maximum and
#' expressed as a proportion of alleles retained.
#'
#' For each target sample size, replicated samples are aggregated to yield
#' the mean, minimum, and maximum proportions of alleles retained. A plot is
#' produced showing (i) the mean rarefaction curve with an uncertainty ribbon
#' (min–max across replicates) and (ii) points for each empirical population at
#' its observed sample size and retained proportion.
#'
#' How to use the output
#'
#' - Assess genetic diversity and sampling sufficiency. The curve indicates
#'   how quickly allelic diversity accumulates with additional individuals, and
#'   where diminishing returns begin.
#' - Interpret population points relative to the curve.
#'   \itemize{
#'     \item Above the curve: population retains more allelic diversity than
#'     expected for its sample size (e.g., unusually high diversity or more
#'     private/low-frequency alleles).
#'     \item On/within the ribbon: diversity consistent with a sample of that
#'     size from a panmictic population with the pooled allele frequencies.
#'     \item Below the curve: population retains fewer alleles than expected,
#'     suggesting reduced diversity (e.g., drift, bottleneck), uneven missingness,
#'     or data-quality issues.
#'   }
#'
#' @return 
#' A list with three elements (the input object is not modified):
#' \itemize{
#'   \item `sim`: `data.frame` with columns `Npop` (sample size),
#'   `mnall` (mean proportion retained), `low` (minimum), and `high` (maximum)
#'   across replicates.
#'   \item `points`: `data.frame` with observed per-population values at their
#'   actual sample sizes (columns include `popname`, `Npop`, and scaled `N.all`).
#'   \item `p1`: a `ggplot` object showing the rarefaction curve, uncertainty
#'   ribbon, and per-population points.
#' }
#'
#' @importFrom ggrepel geom_text_repel
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel makeCluster stopCluster
#' @importFrom dplyr group_by summarise
#' @importFrom doParallel registerDoParallel
#' @importFrom methods is
#' @export
#' @author Author(s): Bernd Gruber. Custodian: Bernd Gruber -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' \donttest{
# if (isTRUE(getOption("dartR_fbm"))) possums.gl <- gl.gen2fbm(possums.gl)
#' dummy <- gl.report.nall(possums.gl[c(1:5,31:35),], simlevels=seq(1,10,3),
#' reps=5, ncores=2)
#' }

gl.report.nall <- function(x,
                           simlevels = seq(1, nInd(x), 5),
                           reps = 10,
                           plot.colors.pop = gl.colors("dis", verbose = 0),
                           ncores = 2,
                           plot.display = TRUE,
                           plot.theme = theme_dartR(),
                           plot.dir = NULL,
                           plot.file = NULL,
                           verbose = NULL) {
  # --- SET VERBOSITY (also suppress plot if verbose == 0)
  verbose <- gl.check.verbosity(verbose)
  if (verbose == 0) {
    plot.display <- FALSE
  }
  
  # --- SET/VALIDATE WORKING DIRECTORY FOR OUTPUTS
  plot.dir <- gl.check.wd(plot.dir, verbose = 0)
  
  # --- FLAG SCRIPT START (for logging/build info)
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   verbose = verbose)
  
  # --- CHECK INPUT DATATYPE (the simulation needs SNP data)
  datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)
  
  # --- FUNCTION SPECIFIC ERROR CHECKING
  is_whole <- function(v, min) {
    is.numeric(v) && length(v) >= 1 && !anyNA(v) && all(v >= min) &&
      all(v %% 1 == 0)
  }
  if (!is_whole(simlevels, 1)) {
    stop(error("  simlevels must be whole numbers of at least 1\n"))
  }
  if (!is_whole(reps, 1) || length(reps) != 1) {
    stop(error("  reps must be a single whole number of at least 1\n"))
  }
  if (!is_whole(ncores, 1) || length(ncores) != 1) {
    stop(error("  ncores must be a single whole number of at least 1\n"))
  }
  
  # --- FUNCTION-SPECIFIC PRE-FILTERING
  # Remove loci that are all NA within populations to avoid artificial inflation/deflation
  x <- gl.filter.allna(x, by.pop = TRUE, verbose = 0)
  
  # --- Helper: compute allele counts per population ---------------------------
  # Returns a data.frame with popname, Npop (sample size), and N.all (total alleles)
  gl.report.nall.pop <- function(x, x2 = NULL) {
    # Either use two provided pops (x, x2) or split x by population factor
    if (!is.null(x2)) {
      pops <- list(pop1 = x, pop2 = x2)
    } else{
      pops <- seppop(x)
    }
    
    res <- data.frame(popname = NA,
                      Npop = NA,
                      N.all = NA)
    
    for (i in 1:length(pops)) {
      # Per-locus allele frequency proxy from genotype means (ignoring NAs)
      p <- colMeans(as.matrix(pops[[i]]), na.rm = TRUE) / 2
      
      # Store population name and size
      res[i, 1] = names(pops)[i]
      res[i, 2] = nInd(pops[[i]])
      
      # Count alleles per locus:
      #   1 allele if fixed (p==0 or p==1), 2 alleles if polymorphic (0<p<1)
      res[i, 3] = sum((p == 1) + (p == 0) + (2 * ((p > 0) & (p < 1))))
    }
    return(res)
  }
  
  # --- Helper: single simulation draw for a given sample size -----------------
  # Draw 'sample' individuals from pooled x; return Npop and N.all
  ss <- function(sample) {
    foundersim <- gl.sim.ind(x, sample, popname = "foundersim", verbose = 0)
    nasim <- gl.report.nall.pop(foundersim)
    res <- c(nasim$Npop, nasim$N.all)
    return(res)
  }
  
  # Predeclare symbols used in dplyr/ggplot NSE to avoid NOTE in R CMD check
  Npop <- Nallsim <- mnall <- low <- high <- N.all <- popname <- ip <- NULL
  
  # --- Reference maximum: total alleles in the pooled dataset -----------------
  # Collapse all individuals into a single population and compute its N.all
  onePop <- x
  pop(onePop) <- rep("A", nInd(onePop))
  maxnall <- gl.report.nall.pop(onePop)$N.all
  
  # Grid of simulation jobs: all combinations of replicate and sample size
  sims <- expand.grid(rep = seq_len(reps), Npop = simlevels)
  
  # --- SIMULATION (parallel unless ncores == 1) --------------------------------
  if (ncores == 1) {
    sims2 <- lapply(sims$Npop, ss)
  } else {
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    sims2 <- foreach::foreach(ip = seq_len(nrow(sims)),
                              .packages = "dartR.sim") %dopar%
      {
        # For each row in 'sims', run a single draw at the requested Npop
        simres <- ss(sims[ip, 2])
        return(simres)
      }
    parallel::stopCluster(cl)
  }
  
  # Collect and label simulation results (one row per draw, also when there
  # is a single draw)
  sims2 <- as.data.frame(do.call(rbind, sims2))
  colnames(sims2) <- c("Npop", "Nallsim")
  
  # Scale simulated allele counts by the pooled maximum to get proportions
  sims2$Nallsim <- sims2$Nallsim / maxnall
  
  # Summarise replicates at each sample size: mean and range (min–max)
  df <- sims2 |>
    dplyr::group_by(Npop) |>
    dplyr::summarise(
      mnall = mean(Nallsim),
      low = min(Nallsim),
      high = max(Nallsim)
    ) |>
    as.data.frame()
  
  # --- Observed per-population points -----------------------------------------
  # Compute observed (scaled) allele counts for each empirical population
  nall.pop <- gl.report.nall.pop(x)
  nall.pop$N.all <- nall.pop$N.all / maxnall
  
  # Preserve the order of populations as they appear in the data
  # pop_order <- unique(as.character(pop(x))) 
  pop_order <- levels(pop(x))
  
  # --- Resolve colors for populations -----------------------------------------
  if (is(plot.colors.pop, "function")) {
    # Generate colors from function given the number of populations
    colors_pops <- plot.colors.pop(length(levels(pop(x))))
  }
  if (!is(plot.colors.pop, "function")) {
    # Use provided vector of colors directly
    colors_pops <- plot.colors.pop
  }
  # Name colors to match population order for consistent mapping
  colors_pops <- setNames(colors_pops, pop_order)
  
  # --- Build plot --------------------------------------------------------------
  ymin <- min(c(df$low, nall.pop$N.all)) * 0.95
  p1 <- ggplot(df, aes(x = Npop, y = mnall)) +
    # Mean rarefaction curve
    geom_line() +
    # Uncertainty band across replicates (min–max)
    geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2) +
    # Axis labels
    xlab("Number of sampled individuals") +
    ylab("Proportion of alleles retained \n relative to pooled maximum") +
    # Clamp y-axis to focus on informative range
    ylim(c(ymin, 1)) +
    # Labels for population points
    geom_text_repel(data = nall.pop,
                    aes(x = Npop, y = N.all, label = popname),
                    hjust = -1) +
    # Observed per-population points
    geom_point(data = nall.pop,
               aes(x = Npop, y = N.all, colour = popname),
               size = 6,
               pch = 16,
               show.legend = FALSE) +
    # Legend title
    guides(color = guide_legend(title = "Population")) +
    # Color mapping with preserved order
    scale_color_manual(values = colors_pops,
                       breaks = pop_order,
                       limits = pop_order) +
    # Apply user-provided theme
    plot.theme
  
  # --- PRINT/SAVE OUTPUTS ------------------------------------------------------
  if (plot.display) {
    print(p1)
  }
  
  if (!is.null(plot.file)) {
    tmp <- utils.plot.save(p1,
                           dir = plot.dir,
                           file = plot.file,
                           verbose = verbose)
  }
  
  # --- FLAG SCRIPT END ---------------------------------------------------------
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  # --- RETURN RESULTS ----------------------------------------------------------
  # sim    : summary table of mean/min/max by sample size
  # points : observed per-population scaled allele counts
  # p1     : ggplot object for the rarefaction curve and points
  return(list(sim = df, points = nall.pop, p1 = p1))
}
