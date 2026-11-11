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
#' any genlight object lacks sim.vars$generation, or an iteration holds the
#' same generation twice, the function stops before running fun. Repeated
#' values in iteration are used once.
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
#' When results are bound, a column missing from some results is filled with
#' NA (for example a statistic that is absent in generations where a
#' population is extinct). NULL results (from a fun run for its side effect,
#' such as saving a file) add no rows to a bound data frame and are NULL in
#' the nested list.
#'
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
    iterations <- unique(iteration)
  }

  ## Every element must be a genlight with sim.vars$generation, and each
  ## generation must appear once per iteration, so fun is not run on data it
  ## cannot be tagged for
  problems <- NULL
  for (it in iterations) {
    gens <- NULL
    for (i in seq_along(x[[it]])) {
      g <- x[[it]][[i]]
      if (!is(g, "genlight") || is.null(g@other$sim.vars$generation)) {
        problems <- c(problems, paste0("iteration ", it, ", element ", i))
      } else {
        gens <- c(gens, as.numeric(g@other$sim.vars$generation))
      }
    }
    if (anyDuplicated(gens)) {
      problems <- c(problems,
                    paste0("iteration ", it, ", generation ",
                           paste(unique(gens[duplicated(gens)]),
                                 collapse = ", "), " appears more than once"))
    }
  }
  if (!is.null(problems)) {
    stop(error("  These elements are not genlight objects with",
               "@other$sim.vars$generation (as stored by gl.sim.WF.run()),",
               "or repeat a generation:",
               paste(problems, collapse = "; "), "\n"))
  }

  # DO THE JOB
  results <- list()
  records <- list() # iteration, generation and result of each call
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
      # NULL results (functions run for their side effect) are kept as NULL
      if (!is.null(res)) {
        attr(res, "iteration") <- it
        attr(res, "generation") <- gen
      }
      res_it[paste0("generation_", gen)] <- list(res)
      records[[length(records) + 1]] <- list(it = it, gen = gen, res = res)
    }
    results[[paste0("iteration_", it)]] <- res_it
    if (verbose >= 2) {
      cat(report("  Iteration", it, ":", length(res_it),
                 "generations processed\n"))
    }
  }

  ## Bind into one data frame when every (non-NULL) result is a data frame or
  ## a vector; columns missing from some results are filled with NA
  non_null <- Filter(function(r) !is.null(r$res), records)
  tabular <- length(non_null) > 0 &&
    all(vapply(non_null, function(r) {
      is.data.frame(r$res) || (is.atomic(r$res) && is.null(dim(r$res)))
    }, logical(1)))

  if (tabular) {
    pieces <- lapply(non_null, function(rec) {
      r <- rec$res
      attr(r, "iteration") <- NULL
      attr(r, "generation") <- NULL
      if (is.data.frame(r)) {
        if (nrow(r) == 0) {
          return(NULL)
        }
        rn <- rownames(r)
        if (!identical(rn, as.character(seq_len(nrow(r))))) {
          r <- cbind(name = rn, r)
        }
        return(cbind(data.frame(iteration = rep(rec$it, nrow(r)),
                                generation = rep(rec$gen, nrow(r))), r))
      }
      nm <- if (is.null(names(r))) seq_along(r) else names(r)
      data.frame(iteration = rep(rec$it, length(r)),
                 generation = rep(rec$gen, length(r)),
                 name = nm, value = as.vector(r))
    })
    results <- as.data.frame(data.table::rbindlist(pieces, use.names = TRUE,
                                                   fill = TRUE))
    rownames(results) <- NULL
  }

  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  return(results)
}
