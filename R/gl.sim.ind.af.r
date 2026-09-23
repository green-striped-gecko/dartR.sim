#' @name gl.sim.ind.af
#' @title Simulates diploid genotypes from per-population allele frequencies
#' @family simulation functions
#' @description
#' This function generates a diploid SNP dataset by sampling genotypes for a
#' specified number of individuals per population from user-provided allele
#' frequencies. The result is returned as a genlight object with population
#' and individual metadata.
#'
#' @details
#' The input df must have three columns: population name, locus name, and the
#' frequency of the alternate allele for that population-locus combination.
#' The alternate allele is the allele counted in the genotypes of a genlight
#' object, and its frequency is the column frequency returned by
#' gl.allele.freq(x, by = "popxloc").
#'
#' For each individual, the genotype at each locus is the number of copies of
#' the alternate allele in two gametes drawn independently with the supplied
#' frequency. The procedure assumes Hardy-Weinberg proportions and linkage
#' equilibrium (loci are sampled independently and there is no
#' within-population structure beyond the supplied allele frequencies).
#'
#' Populations are matched to their frequencies and sizes by name, so the rows
#' of df can be in any order. Every population must provide each locus exactly
#' once.
#'
#' Individual IDs have the form "0_<population name>_<i>". Sex alternates
#' between individuals ("m", "f", "m", ...) and a placeholder phenotype is set
#' to "control" for all individuals; both are stored in
#' other$ind.metrics, together with id, pop, fid (population) and iid
#' (individual ID). Locus alleles are set to the placeholder "G/C", and
#' other$loc.metrics holds the locus names (AlleleID) with all metrics flags
#' reset.
#'
#' @param df A data.frame with three columns: (1) population name,
#'   (2) locus name, and (3) frequency of the alternate allele (numeric in
#'   [0, 1]) [required].
#' @param pop.sizes A vector of whole numbers of at least 1, with the number of
#'   individuals to simulate in each population. If unnamed, one element per
#'   population in the order in which the populations first appear in df; if
#'   named, the names must match the population names in df [required].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' brief progress messages; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#'
#' @return A genlight object with sum(pop.sizes) individuals and one locus per
#' locus in df, with diploid SNP genotypes (0, 1, 2 copies of the alternate
#' allele), pop() set to population names, and the metadata described in
#' Details.
#'
#' @author Author(s): Luis Mijangos. Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @importFrom stats rbinom
#' @export
#'
#' @examples
#' t1 <- gl.filter.callrate(platypus.gl,threshold = 1, mono.rm = TRUE)
#' r1 <- gl.allele.freq(t1, by='popxloc' )
#' r2 <- r1[,c("popn",'locus',"frequency")]
#' res <- gl.sim.ind.af(df = r2, pop.sizes= c(50,50,50))

gl.sim.ind.af <- function(df,
                          pop.sizes,
                          verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   verbose = verbose)

  # FUNCTION SPECIFIC ERROR CHECKING
  if (!is.data.frame(df) || ncol(df) != 3L) {
    stop(error("  df must be a data.frame with exactly three columns:",
               "population, locus, frequency.\n"))
  }
  colnames(df) <- c("popn", "locus", "frequency")
  # Work with names as character: a factor would be indexed by its level
  # codes and would carry unused levels
  df$popn <- as.character(df$popn)
  df$locus <- as.character(df$locus)

  if (!is.numeric(df$frequency) || anyNA(df$frequency)) {
    stop(error("  df$frequency must be numeric with no NAs.\n"))
  }
  if (any(df$frequency < 0 | df$frequency > 1)) {
    stop(error("  All allele frequencies must be in [0, 1].\n"))
  }

  # Unique populations in the order they appear
  pops_in_df <- unique(df$popn)
  n_pops <- length(pops_in_df)

  if (!is.numeric(pop.sizes) || anyNA(pop.sizes) || any(pop.sizes < 1) ||
      any(pop.sizes %% 1 != 0)) {
    stop(error("  pop.sizes must be whole numbers of at least 1, with no",
               "NAs.\n"))
  }

  if (is.null(names(pop.sizes))) {
    # Unnamed: order corresponds to unique(df$popn)
    if (length(pop.sizes) != n_pops) {
      stop(error("  Length of pop.sizes (", length(pop.sizes),
                 ") must equal the number of unique populations in df (",
                 n_pops, "). Alternatively, provide a named vector whose",
                 "names match df$popn.\n"))
    }
    names(pop.sizes) <- pops_in_df
  } else {
    missing_names <- setdiff(pops_in_df, names(pop.sizes))
    extra_names <- setdiff(names(pop.sizes), pops_in_df)
    if (length(missing_names)) {
      stop(error("  Missing sizes for populations:",
                 paste(missing_names, collapse = ", "), "\n"))
    }
    if (length(extra_names)) {
      stop(error("  Unknown population names in pop.sizes:",
                 paste(extra_names, collapse = ", "), "\n"))
    }
  }
  pop.sizes <- as.integer(pop.sizes[pops_in_df])

  dup <- duplicated(df[, c("popn", "locus")])
  if (any(dup)) {
    stop(error("  Each population must provide each locus once. Duplicated",
               "population-locus rows:",
               paste(df$popn[dup], df$locus[dup], sep = "-", collapse = ", "),
               "\n"))
  }

  # Canonical locus order from the whole df
  loci_all <- unique(df$locus)
  n_loci <- length(loci_all)

  df_pops <- split(df, df$popn)
  bad <- vapply(df_pops, function(d) {
    !setequal(d$locus, loci_all)
  }, logical(1))
  if (any(bad)) {
    stop(error("  All populations must provide the same set of loci.",
               "Populations failing this check:",
               paste(names(df_pops)[bad], collapse = ", "), "\n"))
  }

  # DO THE JOB
  # Genotypes of each population, taken by name (split() returns the
  # populations in sorted order, not in order of appearance). The genotype
  # is the sum of two independent gametes, each carrying the alternate allele
  # with probability q, i.e. a binomial draw of size 2.
  geno_list <- lapply(seq_len(n_pops), function(i) {
    d <- df_pops[[pops_in_df[i]]]
    q <- d$frequency[match(loci_all, d$locus)]
    n_ind <- pop.sizes[i]
    matrix(rbinom(n_ind * n_loci, size = 2, prob = rep(q, each = n_ind)),
           nrow = n_ind, ncol = n_loci)
  })
  geno <- do.call(rbind, geno_list)

  pop_vec <- rep(pops_in_df, times = pop.sizes)
  id_vec <- unlist(lapply(seq_len(n_pops), function(i) {
    sprintf("0_%s_%d", pops_in_df[i], seq_len(pop.sizes[i]))
  }))
  # Sex alternates within each population, starting with a male
  sex_vec <- unlist(lapply(pop.sizes, function(n) {
    rep(c("m", "f"), length.out = n)
  }))

  gl <- new(
    "dartR",
    gen = geno,
    ploidy = 2,
    ind.names = id_vec,
    loc.names = loci_all,
    pop = pop_vec
  )

  # Placeholder locus allele labels
  gl@loc.all <- rep("G/C", n_loci)

  # Locus metrics: locus names, with all metrics flagged for recalculation
  gl@other$loc.metrics <- data.frame(AlleleID = loci_all)
  gl <- utils.reset.flags(gl, verbose = 0)

  gl@other$ind.metrics <- data.frame(
    id = id_vec,
    pop = pop_vec,
    fid = pop_vec,
    iid = id_vec,
    sex = factor(sex_vec),
    phenotype = factor(rep("control", nrow(geno))),
    stringsAsFactors = FALSE
  )

  if (verbose >= 3) {
    cat(report("  Simulated", nInd(gl), "individuals in", n_pops,
               "populations at", n_loci, "loci\n"))
  }

  # ADD TO HISTORY
  gl@other$history <- list(match.call())

  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  return(gl)
}
