#' @name gl.sim.apply
#' @title Applies a function to every generation of simulation outputs
#' @family simulation functions
#' @description
#' This function runs any function (for example a dartRverse function) on
#' every genlight object returned by \code{\link{gl.sim.WF.run}} and tags each
#' result with its iteration and generation.
#' @param x Output of \code{\link{gl.sim.WF.run}} (a list of iterations, each
#' a list of genlight objects), a list of genlight objects from one
#' iteration, or a single genlight object [required].
#' @param fun Function to apply to each genlight object, given as a function
#' or its name. For several steps, use an anonymous function, e.g.
#' fun = function(g) gl.pcoa(g, nfactors = 3) [required].
#' @param ... Further arguments passed to fun. They are not checked here;
#' fun checks its own arguments.
#' @param iteration Iterations to use, as numbers (positions in x)
#' [default NULL, all iterations].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' brief progress messages; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @details
#' The generation of each genlight object is read from
#' x@other$sim.vars$generation, which gl.sim.WF.run() stores, not from the
#' names of the list. The iteration is the position of the element in x. If
#' any genlight object lacks sim.vars$generation, the function stops before
#' running fun.
#'
#' How results are returned depends on what fun returns:
#' \itemize{
#' \item Data frames are bound into one data frame, with the columns
#' iteration and generation first. If the data frames have row names other
#' than 1, 2, 3, ..., they are kept in a column called name.
#' \item Vectors (numeric, character or logical) become a data frame with the
#' columns iteration, generation, name (the names of the vector, or its
#' positions) and value.
#' \item If any result is of another type (a genlight, a plot, a dist, a
#' list, ...), nothing is bound: the results are returned as a list of
#' iterations, each a list of generations, named "iteration_1", ... and
#' "generation_1", ...; each result also carries the attributes iteration
#' and generation.
#' }
#' A tidy data frame makes it easy to summarise across replicates, e.g. the
#' mean He per generation over iterations.
#'
#' Iterations that went extinct before any generation was stored are empty
#' and are skipped.
#' @return A data frame or a nested list of results (see details).
#' @author Author(s): Luis Mijangos. Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' ref_table <- gl.sim.WF.table(file_var = system.file("extdata",
#'   "ref_variables.csv", package = "dartR.sim"), interactive_vars = FALSE,
#'   verbose = 0)
#' res_sim <- gl.sim.WF.run(file_var = system.file("extdata",
#'   "sim_variables.csv", package = "dartR.sim"), ref_table = ref_table,
#'   interactive_vars = FALSE, number_iterations = 2, verbose = 0)
#' # expected heterozygosity of every locus, tagged by iteration and generation
#' he <- gl.sim.apply(res_sim, gl.He, verbose = 0)
#' head(he)
#' # mean He per generation across iterations
#' aggregate(value ~ generation, data = he, FUN = mean)
#' @seealso \code{\link{gl.sim.WF.run}}
#' @export

gl.sim.apply <- function(x,
                         fun,
                         ...,
                         iteration = NULL,
                         verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   verbose = verbose)

  # CHECK INPUTS
  fun <- match.fun(fun)

  ## Bring x to the shape of gl.sim.WF.run(): a list of iterations, each a
  ## list of genlight objects
  if (is(x, "genlight")) {
    x <- list(list(x))
  } else if (is.list(x) && length(x) > 0 &&
             all(vapply(x, is, logical(1), "genlight"))) {
    x <- list(x)
  }
  if (!is.list(x) || length(x) == 0 ||
      !all(vapply(x, is.list, logical(1)))) {
    stop(error("  x must be the output of gl.sim.WF.run(), a list of",
               "genlight objects or a genlight object\n"))
  }

  iterations <- seq_along(x)
  if (!is.null(iteration)) {
    if (!is.numeric(iteration) || any(!iteration %in% iterations)) {
      stop(error("  iteration must be numbers between 1 and", length(x),
                 "\n"))
    }
    iterations <- iteration
  }

  ## Every element must be a genlight with sim.vars$generation, so fun is
  ## not run on data it cannot be tagged for
  problems <- NULL
  for (it in iterations) {
    for (i in seq_along(x[[it]])) {
      g <- x[[it]][[i]]
      if (!is(g, "genlight") || is.null(g@other$sim.vars$generation)) {
        problems <- c(problems, paste0("iteration ", it, ", element ", i))
      }
    }
  }
  if (!is.null(problems)) {
    stop(error("  These elements are not genlight objects with",
               "@other$sim.vars$generation (as stored by gl.sim.WF.run()):",
               paste(problems, collapse = "; "), "\n"))
  }

  # DO THE JOB
  results <- list()
  for (it in iterations) {
    if (length(x[[it]]) == 0) {
      if (verbose >= 2) {
        cat(report("  Iteration", it, "has no stored generations; skipped\n"))
      }
      next
    }
    res_it <- list()
    for (g in x[[it]]) {
      gen <- as.numeric(g@other$sim.vars$generation)
      res <- tryCatch(
        fun(g, ...),
        error = function(e) {
          stop(error("  fun failed at iteration", it, ", generation", gen,
                     ":", conditionMessage(e), "\n"), call. = FALSE)
        }
      )
      attr(res, "iteration") <- it
      attr(res, "generation") <- gen
      res_it[[paste0("generation_", gen)]] <- res
    }
    results[[paste0("iteration_", it)]] <- res_it
    if (verbose >= 2) {
      cat(report("  Iteration", it, ":", length(res_it),
                 "generations processed\n"))
    }
  }

  ## Bind into one data frame when every result is a data frame or a vector
  flat <- unlist(results, recursive = FALSE, use.names = FALSE)
  tabular <- length(flat) > 0 &&
    all(vapply(flat, function(r) {
      is.data.frame(r) || (is.atomic(r) && is.null(dim(r)))
    }, logical(1)))

  if (tabular) {
    results <- do.call(rbind, lapply(flat, function(r) {
      tags <- data.frame(iteration = attr(r, "iteration"),
                         generation = attr(r, "generation"))
      if (is.data.frame(r)) {
        attr(r, "iteration") <- NULL
        attr(r, "generation") <- NULL
        rn <- rownames(r)
        if (!identical(rn, as.character(seq_len(nrow(r))))) {
          r <- cbind(name = rn, r)
        }
        if (nrow(r) == 0) {
          return(NULL)
        }
        return(cbind(tags[rep(1, nrow(r)), ], r))
      }
      nm <- if (is.null(names(r))) seq_along(r) else names(r)
      data.frame(tags[rep(1, length(r)), ], name = nm,
                 value = as.vector(r))
    }))
    rownames(results) <- NULL
  }

  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  return(results)
}
