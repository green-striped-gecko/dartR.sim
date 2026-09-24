# Characterization tests for gl.report.nall
# Baseline snapshotted before review (review-gl.report.nall), bugs included.
# Assertions marked [approved diff] were flipped in Phase C under the finding
# named in the test title; see function-review/reports/dartR.sim/.
# The simulation runs on parallel workers, which load the installed
# dartR.sim; results depend on the installed gl.sim.ind().

x <- gl.filter.allna(platypus.gl, by.pop = TRUE, verbose = 0)
run_nall <- function(...) {
  suppressWarnings(gl.report.nall(x, ncores = 1, plot.display = FALSE, ...))
}

test_that("returned structure; input untouched; results without plotting", {
  x0 <- x
  set.seed(1)
  r <- run_nall(simlevels = c(5, 20, 61), reps = 3, verbose = 0)
  expect_named(r, c("sim", "points", "p1"))
  expect_identical(r$sim$Npop, c(5, 20, 61))
  expect_true(all(r$sim$low <= r$sim$mnall & r$sim$mnall <= r$sim$high))
  expect_false(is.unsorted(r$sim$mnall))
  expect_identical(r$points$popname, popNames(x))
  expect_identical(r$points$Npop, as.vector(table(pop(x))))
  expect_true(all(r$points$N.all > 0.5 & r$points$N.all <= 1))
  expect_s3_class(r$p1, "ggplot")
  expect_identical(x, x0)
})

test_that("curve is simulated and stays below 1 at the full sample (documented, F1)", {
  set.seed(1)
  r <- run_nall(simlevels = nInd(x), reps = 3, verbose = 0)
  # a subsample of every individual would hold every allele (1.0)
  expect_lt(r$sim$mnall, 0.99)
})

test_that("[approved diff] verbose = 0 is silent (F2)", {
  expect_silent(run_nall(simlevels = 5, reps = 2, verbose = 0))
})

test_that("[approved diff] SilicoDArT stops at the start (F3)", {
  gs <- gl.filter.allna(testset.gs, verbose = 0)
  expect_error(gl.report.nall(gs, simlevels = 5, reps = 2, ncores = 1,
                              plot.display = FALSE, verbose = 0),
               "SilicoDArT")
})

test_that("[approved diff] single job works; arguments checked (F4)", {
  set.seed(1)
  r <- run_nall(simlevels = 5, reps = 1, verbose = 0)
  expect_identical(nrow(r$sim), 1L)
  expect_identical(r$sim$low, r$sim$high)
  expect_error(run_nall(simlevels = 5, reps = 0, verbose = 0), "reps")
  expect_error(run_nall(simlevels = c(5, 2.5), reps = 2, verbose = 0),
               "simlevels")
  expect_error(run_nall(simlevels = 5, reps = 2, ncores = 0, verbose = 0),
               "ncores")
})

test_that("[approved diff] parallel and single-core runs agree (F5)", {
  set.seed(1)
  r1 <- run_nall(simlevels = c(5, 40), reps = 4, verbose = 0)
  r2 <- suppressWarnings(gl.report.nall(x, simlevels = c(5, 40), reps = 4,
                                        ncores = 2, plot.display = FALSE,
                                        verbose = 0))
  expect_identical(r1$sim$Npop, r2$sim$Npop)
  expect_equal(r1$sim$mnall, r2$sim$mnall, tolerance = 0.03)
  expect_identical(r1$points, r2$points)
})

test_that("[approved diff] sim is a data.frame (F6)", {
  set.seed(1)
  r <- run_nall(simlevels = 5, reps = 2, verbose = 0)
  expect_identical(class(r$sim), "data.frame")
  expect_identical(colnames(r$sim), c("Npop", "mnall", "low", "high"))
})
