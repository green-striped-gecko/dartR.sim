#' @name gl.sim.emigration
#' @title Simulates emigration between populations
#' @family simulation functions
#' @description
#' A function that allows to exchange individuals of populations within a
#' genlight object (=simulate emigration between populations).
#'
#' @details
#' There are two ways to specify emigration. If an emi.table is provided (a
#' square matrix of dimension of the populations that specifies the emigration
#' from column x to row y), then emigration is deterministic in terms of numbers
#' of individuals as specified in the table. If perc.mig and emi.m are provided,
#' then emigration is probabilistic: each individual emigrates with probability
#' perc.mig, and the population it moves to is drawn from the relative
#' probabilities in the column of emi.m of its population (from = column,
#' to = row). A column of zeros means that no individual leaves that
#' population.
#'
#' If the diagonal of emi.m is non zero, emigrants can be assigned to their own
#' population, where they stay. So most often you want to set the diagonal of
#' the emi.m matrix to zero. The diagonal of emi.table is ignored.
#'
#' Emigrants are drawn at random from the individuals present in each
#' population at the start of the call, and all move at the same time, so an
#' individual moves at most once per call. The returned object is grouped by
#' population; the population of each moved individual is updated in pop() and
#' in other$ind.metrics$pop.
#'
#' A population can lose all its individuals. In a returned genlight object it
#' no longer appears in popNames(); in a returned list its element is NULL.
#' A list with NULL elements is accepted as input, so repeated calls can use
#' the same matrices.
#'
#' @param x A genlight object, or a list of genlight objects with one element
#' per population (NULL for an empty population) [required].
#' @param perc.mig Proportion of individuals that emigrate, between 0 and 1
#' (each individual emigrates with this probability) [default NULL].
#' @param emi.m Probabilistic emigration matrix (emigrate from = column,
#' to = row), a square matrix with one row and one column per population
#' [default NULL].
#' @param emi.table If provided, emi.m and perc.mig are ignored. Deterministic
#' emigration as specified in the matrix (a square matrix of dimension of the
#' number of populations). e.g. an entry in the 'emi.table[2,1]<- 5' means that
#' five individuals emigrate from population 1 to population 2 (from=columns and
#'  to=row) [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' brief progress messages; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @return A list or a single [depends on the input] genlight object, where
#' emigration between population has happened
#' @author Author(s): Bernd Gruber. Custodian: Bernd Gruber -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' x <- possums.gl
#' #one individual moves from every population to
#' #every other population
#' emi.tab <- matrix(1, nrow=nPop(x), ncol=nPop(x))
#' diag(emi.tab)<- 0
#' np <- gl.sim.emigration(x, emi.table=emi.tab)
#' np
#' @importFrom stats rbinom
#' @export

gl.sim.emigration <- function(x,
                              perc.mig = NULL,
                              emi.m = NULL,
                              emi.table = NULL,
                              verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   verbose = verbose)

  # CHECK DATATYPE
  is_list <- is.list(x) && !is(x, "genlight")
  if (is_list) {
    if (length(x) == 0) {
      stop(error("  x must be a genlight object or a non-empty list of",
                 "genlight objects\n"))
    }
    present <- !vapply(x, is.null, logical(1))
    if (!any(present) ||
        !all(vapply(x[present], is, logical(1), "genlight"))) {
      stop(error("  Every element of the list x must be a genlight object",
                 "(or NULL for an empty population)\n"))
    }
    for (g in x[present]) {
      datatype <- utils.check.datatype(g, verbose = 0)
    }
    loc_ref <- locNames(x[present][[1]])
    if (!all(vapply(x[present], function(g) {
      identical(locNames(g), loc_ref)
    }, logical(1)))) {
      stop(error("  All genlight objects in x must have the same loci in the",
                 "same order\n"))
    }
  } else {
    datatype <- utils.check.datatype(x, verbose = 0)
  }

  # Populations and the population of each individual
  if (is_list) {
    pn <- names(x)
    if (is.null(pn)) pn <- rep("", length(x))
    # Unnamed elements are named after their population
    for (i in which(pn == "" | is.na(pn))) {
      lv <- if (present[i]) unique(as.character(pop(x[[i]]))) else NULL
      pn[i] <- if (length(lv) == 1 && !is.na(lv)) lv else paste0("pop", i)
    }
    if (anyDuplicated(pn)) {
      stop(error("  The names of the list x must be unique\n"))
    }
    gl <- do.call("rbind", unname(x[present]))
    cur <- rep(pn[present], vapply(x[present], nInd, integer(1)))
  } else {
    pn <- popNames(x)
    gl <- x
    cur <- as.character(pop(x))
  }
  n.pops <- length(pn)
  sizes <- as.vector(table(factor(cur, levels = pn)))

  # FUNCTION SPECIFIC ERROR CHECKING
  check_matrix <- function(m, name) {
    if (is.data.frame(m)) m <- as.matrix(m)
    if (!is.matrix(m) || !is.numeric(m) || anyNA(m) || any(m < 0)) {
      stop(error("  ", name, " must be a numeric matrix with no NAs or",
                 "negative values\n"))
    }
    if (nrow(m) != n.pops || ncol(m) != n.pops) {
      stop(error("  ", name, " must be a square matrix with one row and one",
                 "column per population (", n.pops, "); it is ", nrow(m),
                 " x ", ncol(m), "\n"))
    }
    dimnames(m) <- NULL
    return(m)
  }

  if (!is.null(emi.table)) {
    migs <- check_matrix(emi.table, "emi.table")
    if (any(migs %% 1 != 0)) {
      stop(error("  emi.table must contain whole numbers of individuals\n"))
    }
    diag(migs) <- 0 # staying is not emigration
    too_many <- colSums(migs) > sizes
    if (any(too_many)) {
      stop(error("  emi.table asks for more emigrants than individuals in",
                 "population(s):",
                 paste0(pn[too_many], " (", colSums(migs)[too_many], " of ",
                        sizes[too_many], ")", collapse = ", "), "\n"))
    }
  } else {
    if (is.null(perc.mig) || is.null(emi.m)) {
      stop(error("  Provide either emi.table, or both perc.mig and emi.m\n"))
    }
    if (!is.numeric(perc.mig) || length(perc.mig) != 1 || is.na(perc.mig) ||
        perc.mig < 0 || perc.mig > 1) {
      stop(error("  perc.mig must be a single proportion between 0 and 1\n"))
    }
    emi.m <- check_matrix(emi.m, "emi.m")
    # Number of emigrants of each population (from = column) and their
    # destinations (to = row), drawn from the column's relative probabilities
    migs <- matrix(0, nrow = n.pops, ncol = n.pops)
    for (from in seq_len(n.pops)) {
      if (sum(emi.m[, from]) == 0) next # nobody leaves this population
      n_emi <- rbinom(1, sizes[from], perc.mig)
      if (n_emi == 0) next
      dest <- sample.int(n.pops, n_emi, replace = TRUE, prob = emi.m[, from])
      migs[, from] <- tabulate(dest, nbins = n.pops)
    }
    diag(migs) <- 0 # emigrants assigned to their own population stay
  }

  # DO THE JOB
  # Emigrants are drawn from the residents at the start and all move at once,
  # so an individual moves at most once
  new_pop <- cur
  for (from in seq_len(n.pops)) {
    n_out <- sum(migs[, from])
    if (n_out == 0) next
    residents <- which(cur == pn[from])
    movers <- residents[sample.int(length(residents), n_out)]
    new_pop[movers] <- rep(pn, times = migs[, from])
  }

  # Group individuals by population (stable, so the input order is kept
  # within each population)
  ord <- order(factor(new_pop, levels = pn))
  gl <- gl[ord]
  pop(gl) <- factor(new_pop[ord], levels = pn)
  if (!is.null(gl@other$ind.metrics$pop)) {
    gl@other$ind.metrics$pop <- new_pop[ord]
  }

  new_sizes <- as.vector(table(factor(new_pop, levels = pn)))
  emptied <- pn[new_sizes == 0]
  if (verbose >= 1 && length(emptied) > 0) {
    cat(warn("  Warning: population(s) with no individuals left:",
             paste(emptied, collapse = ", "), "\n"))
  }
  if (verbose >= 3) {
    cat(report("  Moved", sum(migs), "individuals between", n.pops,
               "populations\n"))
  }

  # return list or single genlight object (depending on the input)
  if (is_list) {
    xout <- setNames(vector("list", n.pops), pn)
    for (i in which(new_sizes > 0)) {
      xout[[i]] <- gl[new_pop[ord] == pn[i]]
    }
  } else {
    # ADD TO HISTORY
    nh <- length(gl@other$history)
    gl@other$history[[nh + 1]] <- match.call()
    xout <- gl
  }

  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  return(xout)
}
