#' @name gl.report.nall
#' @title
#' Report allelic retention and simulate a rarefaction curve
#'
#' @description
#' This function reports per-population allele counts and simulates a
#' rarefaction-style curve showing the proportion of the dataset’s total allelic
#' diversity captured as progressively more individuals are sampled.
#'
#' @param x Name of the genlight/dartR object containing the SNP data.
#' The object needs to have no missing data as subsampling from missing data
#' is not possible. So we recommend to filter by callrate using a threshold
#' of 1 [required].
#' @param simlevels A vector that defines the different levels the combined
#' population should be subsampled [default seq(1,nInd(x),5)].
#' @param reps Number of replicate subsamples per sample size [default 10].
#' @param plot.colors.pop A color palette for population plots or a list with
#' as many colors as there are populations in the dataset 
#' [default gl.colors("dis")].
#' @param ncores Number of cores to be used for parallel processing 
#' [default 10].
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
#' The function estimates how sampling effort affects observed allelic diversity
#' by repeatedly subsampling individuals from the pooled set of all
#' individuals at user-defined sample sizes (`simlevels`), with each subsample
#' replicated (`reps` times). The maximum attainable allele count is first
#' determined by pooling all individuals into a single group; all simulation
#' outputs and per-population observations are then normalized to this pooled
#' maximum and expressed as a proportion of alleles retained.
#'
#' For each target sample size, replicated subsamples are aggregated to yield
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
#'     \item On/within the ribbon: diversity consistent with random sampling
#'     from the pooled dataset at that size.
#'     \item Below the curve: population retains fewer alleles than expected,
#'     suggesting reduced diversity (e.g., drift, bottleneck), uneven missingness,
#'     or data-quality issues.
#'   }
#'
#' @return 
#' A list with three elements:
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
#' @author Custodian: Bernd Gruber -- Post to
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
                           plot.colors.pop = gl.colors("dis"),
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
                   build = "v.2023.3",
                   verbose = verbose)
  
  # --- CHECK INPUT DATATYPE (returns a label; stored but not used later)
  datatype <- utils.check.datatype(x, verbose = verbose)
  
  # --- FUNCTION-SPECIFIC PRE-FILTERING
  # Remove loci that are all NA within populations to avoid artificial inflation/deflation
  x <- gl.filter.allna(x, by.pop = TRUE)
  
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
    foundersim <- gl.sim.ind(x, sample, popname = "foundersim")
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
  sims <- expand.grid(rep = 1:reps, Npop = simlevels)
  
  # --- PARALLEL SIMULATION -----------------------------------------------------
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  sims2 <- foreach::foreach(ip = 1:nrow(sims),
                            .combine = rbind,
                            .packages = "dartR.sim") %dopar%
    {
      # For each row in 'sims', run a single draw at the requested Npop
      simres <- ss(sims[ip, 2])
      return(simres)
    }
  parallel::stopCluster(cl)
  
  # Collect and label simulation results
  sims2 <- as.data.frame(sims2)
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
    )
  
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
