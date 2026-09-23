#' @name gl.diagnostics.sim
#' @title Comparing simulations against theoretical expectations
#' @family simulation functions
#' @param x Output from function \code{\link{gl.sim.WF.run}} [required].
#' @param Ne Effective population size used for the theoretical expectations:
#' one value, or one value per population [required].
#' @param iteration Iteration number to analyse [default 1].
#' @param pop_he Population (position in popNames) in which the rate of loss
#' of heterozygosity is compared against theoretical expectations
#' [default 1].
#' @param pops_fst Pair of populations (positions in popNames) in which FST is
#' compared against theoretical expectations [default c(1,2)].
#' @param plot_theme User specified theme [default theme_dartR()].
#' @param plot.file Name for the RDS binary file to save (base name only,
#' exclude extension) [default NULL].
#' @param plot.dir Directory in which to save files [default = working
#' directory].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' brief progress messages; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @details
#' Two plots compare the simulations against theoretical expectations:
#' 
#' \enumerate{
#' \item Heterozygosity. The expected He under drift in an isolated
#' population (Crow & Kimura, 1970, p. 329) is
#' 
#' He(t) = He(t0) (1 - 1/(2Ne))^(t - t0),
#' 
#' where t0 is the first stored generation and He(t0) the observed He in
#' population pop_he at t0. Five curves are drawn, for Ne, 1.25 Ne, ... 2 Ne,
#' to show how sensitive the expectation is to Ne. Migration slows the loss
#' of heterozygosity, so with dispersal the observed He is expected to fall
#' more slowly than these curves.
#' 
#' \item FST. The expected FST between the populations in pops_fst is
#' computed generation by generation with the identity-by-descent recursion
#' of the island model simulated by \code{\link{gl.sim.WF.run}}: in each
#' generation, genes are exchanged between populations (a fraction m of
#' each population are immigrants, m = number_transfers / N in generations
#' with dispersal, coming equally from the other populations), then drift
#' acts within populations (Ne). FST = (F0 - F1) / (1 - F1), where F0 and F1
#' are the probabilities of identity by descent of two genes within and
#' between populations. The recursion starts at the first stored generation
#' of phase 2 from the observed FST. At equilibrium, for n populations and
#' small m, this is close to 1 / (1 + 4 Ne m n / (n - 1)) (Takahata, 1983,
#' with m counting only immigrants from other populations). The expectation
#' is available without dispersal (m = 0) and for dispersal_type
#' "all_connected"; "line", "circle" and dispersal files are not island
#' models and stop the function.
#' }
#' Observed FST is Nei's estimator (hierfstat::pairwise.neifst).
#' @return Invisibly, a list with the plot (plot), a data frame of observed
#' and expected He by generation and Ne (he), and a data frame of observed
#' and expected FST by generation (fst). The plot is printed.
#' @author Author(s): Luis Mijangos. Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' \donttest{
#' ref_table <- gl.sim.WF.table(file_var=system.file('extdata', 
#' 'ref_variables.csv', package = 'dartR.sim'),interactive_vars = FALSE)
#' 
#' res_sim <- gl.sim.WF.run(file_var = system.file('extdata',
#'  'sim_variables.csv', package ='dartR.sim'),ref_table=ref_table,
#'  interactive_vars = FALSE,number_pops_phase2=2,population_size_phase2="10 10")
#'  
#'  res <- gl.diagnostics.sim(x=res_sim, Ne=10)
#'  head(res$fst)
#'  }
#'@references
#'\itemize{
#'\item Crow JF, Kimura M. An introduction to population genetics theory.
#' Harper and Row, New York. 1970.
#'\item Takahata N. Gene identity and genetic differentiation of populations in 
#'the finite island model. Genetics. 1983;104(3):497-512.
#'  }
#' @export

gl.diagnostics.sim <- function(x,
                               Ne,
                               iteration = 1,
                               pop_he = 1,
                               pops_fst = c(1, 2),
                               plot_theme = theme_dartR(),
                               plot.file = NULL,
                               plot.dir = NULL,
                               verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir, verbose = 0)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   verbose = verbose)
  
  # CHECK INPUTS
  if (!is.list(x) || !is.numeric(iteration) || length(iteration) != 1 ||
      !iteration %in% seq_along(x)) {
    stop(error("  iteration must be one iteration of x, which must be the",
               "output of gl.sim.WF.run()\n"))
  }
  x <- x[[iteration]]
  if (length(x) == 0) {
    stop(error("  Iteration", iteration, "has no stored generations (the",
               "populations went extinct)\n"))
  }
  if (!all(vapply(x, function(g) {
    is(g, "genlight") && !is.null(g@other$sim.vars$generation)
  }, logical(1)))) {
    stop(error("  x must be the output of gl.sim.WF.run()\n"))
  }
  n_pops <- nPop(x[[1]])
  if (n_pops < 2) {
    stop(error("  At least two populations are needed to compare FST\n"))
  }
  if (!all(c(pop_he, pops_fst) %in% seq_len(n_pops)) ||
      length(pop_he) != 1 || length(pops_fst) != 2 ||
      pops_fst[1] == pops_fst[2]) {
    stop(error("  pop_he must be one population and pops_fst two different",
               "populations, given as positions between 1 and", n_pops,
               "\n"))
  }
  
  # Simulation variables are stored as text, possibly quoted and space
  # delimited
  sim_num <- function(v) {
    as.numeric(unlist(strsplit(trimws(gsub("[\"']", "", v)), " +")))
  }
  sim_vars <- x[[1]]@other$sim.vars
  Ne <- sim_num(Ne)
  Ne <- if (length(Ne) == 1) rep(Ne, n_pops) else Ne
  if (length(Ne) != n_pops || anyNA(Ne)) {
    stop(error("  Ne must be one value or one value per population\n"))
  }
  
  # Dispersal settings of phase 2 (the island recursion needs all_connected)
  dispersal <- as.logical(gsub("[\"']", "", sim_vars$dispersal_phase2))
  if (dispersal == TRUE) {
    dispersal_type <- gsub("[\"']", "", sim_vars$dispersal_type_phase2)
    if (!is.null(sim_vars$file_dispersal) || dispersal_type != "all_connected") {
      stop(error("  The expected FST follows the island model and needs",
                 "dispersal_type = 'all_connected' (not a dispersal file,",
                 "'line' or 'circle')\n"))
    }
    number_transfers <- unique(sim_num(sim_vars$number_transfers_phase2))
    transfer_each_gen <- unique(sim_num(sim_vars$transfer_each_gen_phase2))
  } else {
    number_transfers <- 0
    transfer_each_gen <- 1
  }
  population_size <- sim_num(sim_vars$population_size_phase2)
  
  # DO THE JOB
  lab <- gen <- He <- value <- variable <- fst_obs <- expected <- NULL
  
  sep_pops <- lapply(x, seppop)
  generations_sim <- vapply(x, function(y) {
    as.numeric(y@other$sim.vars$generation)
  }, numeric(1))
  
  ####################### He #######################
  
  he_obs <- vapply(sep_pops, function(y) mean(gl.He(y[[pop_he]])),
                   numeric(1))
  
  # Expected He decays from the observed He at the first stored generation
  Ne_he <- seq(Ne[pop_he], Ne[pop_he] * 2, Ne[pop_he] / 4)
  first <- which.min(generations_sim)
  expected_het_3 <- do.call(rbind, lapply(Ne_he, function(ne) {
    data.frame(He = he_obs[first] *
                 (1 - 1 / (2 * ne)) ^ (generations_sim - generations_sim[first]),
               Ne = as.character(ne),
               gen = generations_sim)
  }))
  expected_het_3$lab <- factor(paste("Ne ", expected_het_3$Ne),
                               levels = paste("Ne ", Ne_he))
  
  he_pop <- data.frame(gen = generations_sim, value = he_obs,
                       variable = paste0("pop", pop_he))
  
  p1 <- ggplot(data=expected_het_3, aes(x=gen,y=He,color=lab)) +
    geom_line(linewidth=0.75,linetype = "dashed") +
    geom_line(data=he_pop,aes(x=gen,y=value,color=variable),linewidth=1.5) +
    labs(x="Generations",
         y="He", 
         title=paste("Rate of loss of heterozygosity\nacross generations population",
                     paste(pop_he,collapse = " ") ))+ 
    plot_theme +
    theme(legend.title=element_blank())
  
  ####################### FST #######################
  
  fst_gen <- vapply(sep_pops, function(y) {
    merge_pop <- Reduce(rbind, y[pops_fst])
    temp <- hierfstat::genind2hierfstat(gl2gi(merge_pop, verbose = 0))
    hierfstat::pairwise.neifst(temp)[1, 2]
  }, numeric(1))
  
  # Island-model recursion of identity by descent within (F0) and between
  # (F1) populations: migration (a fraction m of immigrants, from each of the
  # other n - 1 populations equally), then drift. It starts at the first
  # stored generation of phase 2, with F1 = 0 and F0 = observed FST
  phase1 <- as.logical(gsub("[\"']", "", sim_vars$phase1))
  start_phase2 <- if (isTRUE(phase1)) {
    sim_num(sim_vars$gen_number_phase1) + 1
  } else {
    1
  }
  stored_phase2 <- which(generations_sim >= start_phase2)
  start <- stored_phase2[which.min(generations_sim[stored_phase2])]
  Ne_fst <- mean(Ne[pops_fst])
  N_fst <- mean(population_size[pops_fst])
  m <- (n_pops - 1) * number_transfers / N_fst
  same <- (1 - m) ^ 2 + m ^ 2 / (n_pops - 1)
  p_same_between <- 2 * (1 - m) * m / (n_pops - 1) +
    (n_pops - 2) * m ^ 2 / (n_pops - 1) ^ 2
  F0 <- fst_gen[start]
  F1 <- 0
  gens_expected <- generations_sim[start]:max(generations_sim)
  fst_expected <- numeric(length(gens_expected))
  fst_expected[1] <- F0
  for (i in seq_along(gens_expected)[-1]) {
    g <- gens_expected[i]
    # gl.sim.WF.run moves individuals when g %% transfer_each_gen == 0
    if (m > 0 && g != 1 && g %% transfer_each_gen == 0) {
      F0_new <- same * F0 + (1 - same) * F1
      F1 <- p_same_between * F0 + (1 - p_same_between) * F1
      F0 <- F0_new
    }
    F0 <- 1 / (2 * Ne_fst) + (1 - 1 / (2 * Ne_fst)) * F0
    fst_expected[i] <- (F0 - F1) / (1 - F1)
  }
  
  generations_fst <- data.frame(gen = generations_sim, fst_obs = fst_gen)
  expected_fst <- data.frame(gen = gens_expected, expected = fst_expected)
  
  p2 <- ggplot(generations_fst) +
    geom_line(aes(x = gen, y = fst_obs, color = "Fst observed"),
              linewidth = 1) +
    geom_line(data = expected_fst,
              aes(x = gen, y = expected, color = "Fst expected"),
              linewidth = 1, linetype = "dashed") +
    labs(x = "Generations", y = "Fst",
         title = paste("Fst between populations:",
                       paste(pops_fst, collapse = " "))) +
    scale_color_manual(values = c("Fst observed" = "deeppink",
                                  "Fst expected" = "chartreuse4")) +
    plot_theme +
    theme(legend.title = element_blank())
  
  # PRINTING OUTPUTS
  # using package patchwork
  p3 <- (p1 / p2)
  print(p3)
  
  # Optionally save the plot ---------------------
  
  if(!is.null(plot.file)){
    tmp <- utils.plot.save(p3,
                           dir=plot.dir,
                           file=plot.file,
                           verbose=verbose)
  }
  
  he_table <- merge(data.frame(gen = generations_sim, observed = he_obs),
                    data.frame(gen = expected_het_3$gen,
                               Ne = as.numeric(expected_het_3$Ne),
                               expected = expected_het_3$He))
  fst_table <- merge(generations_fst, expected_fst, all = TRUE)
  names(fst_table) <- c("gen", "observed", "expected")
  
  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }
  
  # RETURN
  return(invisible(list(plot = p3, he = he_table, fst = fst_table)))
}
