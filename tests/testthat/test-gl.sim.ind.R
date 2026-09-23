# Characterization tests for gl.sim.ind
# Baseline snapshotted before review (review-gl.sim.ind), bugs included.
# Assertions marked [approved diff] were flipped in Phase C under the finding
# named in the test title; see function-review/reports/dartR.sim/.

x_clean <- gl.filter.allna(testset.gl, verbose = 0)

test_that("output structure and allele frequencies", {
  set.seed(1)
  s <- gl.sim.ind(x_clean, n = 2000, popname = "sims", verbose = 0)
  expect_s4_class(s, "genlight")
  expect_identical(c(nInd(s), nLoc(s)), c(2000L, nLoc(x_clean)))
  expect_true(all(ploidy(s) == 2))
  expect_identical(levels(pop(s)), "sims")
  expect_identical(locNames(s), locNames(x_clean))
  expect_identical(s@loc.all, x_clean@loc.all)
  p_in <- colMeans(as.matrix(x_clean), na.rm = TRUE) / 2
  p_out <- colMeans(as.matrix(s)) / 2
  expect_lt(max(abs(p_out - p_in)), 0.04)
  # Hardy-Weinberg: observed heterozygosity close to 2pq
  expect_equal(mean(as.matrix(s) == 1), mean(2 * p_out * (1 - p_out)),
               tolerance = 0.05)
})

test_that("[approved diff] loci without calls get NA genotypes (F1)", {
  s <- gl.sim.ind(testset.gl, n = 5, verbose = 0)
  no_calls <- is.nan(colMeans(as.matrix(testset.gl), na.rm = TRUE))
  expect_identical(nLoc(s), nLoc(testset.gl))
  expect_true(all(is.na(as.matrix(s)[, no_calls])))
  expect_false(anyNA(as.matrix(s)[, !no_calls]))
})

test_that("[approved diff] SilicoDArT input stops (F2)", {
  gs <- gl.filter.allna(testset.gs, verbose = 0)
  expect_error(gl.sim.ind(gs, n = 5, verbose = 0))
})

test_that("[approved diff] chromosome and loc.metrics are carried (F3)", {
  x <- x_clean
  x@chromosome <- factor(rep(c("1", "2"), length.out = nLoc(x)))
  s <- gl.sim.ind(x, n = 5, verbose = 0)
  expect_identical(as.character(s@chromosome), as.character(x@chromosome))
  expect_equal(nrow(s@other$loc.metrics), nLoc(x))
  expect_false(any(unlist(s@other$loc.metrics.flags[
    c("monomorphs", "OneRatioRef", "OneRatioSnp", "AvgPIC")]), na.rm = TRUE))
})

test_that("[approved diff] history and verbosity (F5)", {
  s <- gl.sim.ind(x_clean, n = 5, verbose = 0)
  expect_length(s@other$history, 1)
  expect_output(gl.sim.ind(x_clean, n = 5, verbose = 1), "Completed")
})

test_that("[approved diff] n is validated (F6)", {
  expect_error(gl.sim.ind(x_clean, n = 0, verbose = 0), "whole number")
  expect_error(gl.sim.ind(x_clean, n = 2.5, verbose = 0), "whole number")
})
