# Tests for gl.sim.apply

fv_ref <- system.file("extdata", "ref_variables.csv", package = "dartR.sim")
fv_sim <- system.file("extdata", "sim_variables.csv", package = "dartR.sim")
rt <- gl.sim.WF.table(file_var = fv_ref, interactive_vars = FALSE,
                      verbose = 0, seed = 1)
sims <- gl.sim.WF.run(file_var = fv_sim, ref_table = rt,
                      interactive_vars = FALSE, verbose = 0, seed = 1,
                      number_iterations = 2)

test_that("vectors are bound with iteration, generation, name, value", {
  he <- gl.sim.apply(sims, gl.He, verbose = 0)
  expect_s3_class(he, "data.frame")
  expect_named(he, c("iteration", "generation", "name", "value"))
  expect_equal(nrow(he), 2 * 2 * 100)
  expect_setequal(unique(he$generation), c(1, 10))
  g <- sims[[2]][["generation_10"]]
  expect_equal(he$value[he$iteration == 2 & he$generation == 10],
               unname(gl.He(g)))
  expect_equal(he$name[1:3], locNames(sims[[1]][[1]])[1:3])
})

test_that("data frames are bound; ... reaches fun", {
  f <- function(g, k) data.frame(n = nInd(g) * k)
  res <- gl.sim.apply(sims, f, k = 2, verbose = 0)
  expect_named(res, c("iteration", "generation", "n"))
  expect_equal(nrow(res), 4)
  expect_true(all(res$n == 2 * 26))
})

test_that("other results come back as a tagged nested list", {
  res <- gl.sim.apply(sims, function(g) g[1:5, ], verbose = 0)
  expect_named(res, c("iteration_1", "iteration_2"))
  expect_named(res[[1]], c("generation_1", "generation_10"))
  expect_s4_class(res[[2]][[2]], "genlight")
  expect_equal(attr(res[[2]][[2]], "generation"), 10)
  expect_equal(attr(res[[2]][[2]], "iteration"), 2)
})

test_that("flat lists, single genlights and iteration subsets", {
  expect_equal(unique(gl.sim.apply(sims[[2]], gl.He, verbose = 0)$iteration),
               1)
  expect_equal(nrow(gl.sim.apply(sims[[1]][[1]], nInd, verbose = 0)), 1)
  res <- gl.sim.apply(sims, "nInd", iteration = 2, verbose = 0)
  expect_equal(unique(res$iteration), 2)
  expect_error(gl.sim.apply(sims, nInd, iteration = 3, verbose = 0),
               "iteration")
})

test_that("generation comes from sim.vars, not from list names", {
  s <- sims
  names(s[[1]]) <- c("a", "b")
  res <- gl.sim.apply(s, nInd, verbose = 0)
  expect_equal(res$generation[res$iteration == 1], c(1, 10))
})

test_that("elements without sim.vars$generation stop before fun runs", {
  s <- sims
  s[[1]][[2]]@other$sim.vars <- NULL
  expect_error(gl.sim.apply(s, nInd, verbose = 0),
               "iteration 1, element 2")
  expect_error(gl.sim.apply(list(list(1:3)), nInd, verbose = 0),
               "iteration 1, element 1")
  expect_error(gl.sim.apply(1:3, nInd, verbose = 0), "x must be")
})

test_that("empty iterations are skipped; errors in fun are located", {
  s <- sims
  s[[1]] <- list()
  res <- gl.sim.apply(s, nInd, verbose = 0)
  expect_equal(unique(res$iteration), 2)
  expect_error(gl.sim.apply(sims, function(g) stop("boom"), verbose = 0),
               "iteration 1 , generation 1 : boom")
})
