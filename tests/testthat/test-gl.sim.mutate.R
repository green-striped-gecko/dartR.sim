# Characterization tests for gl.sim.mutate
# Baseline snapshotted before review (review-gl.sim.mutate), bugs included.
# Assertions marked [approved diff] were flipped in Phase C under the finding
# named in the test title; see function-review/reports/dartR.sim/.

x <- gl.filter.allna(testset.gl, verbose = 0)
n_diff <- function(a, b) sum(as.matrix(a) != as.matrix(b), na.rm = TRUE)

test_that("mutations change one allele and leave missing data alone", {
  set.seed(5)
  b <- gl.sim.mutate(x, mut.rate = 1e-3, verbose = 0)
  tab <- table(before = as.matrix(x), after = as.matrix(b))
  # homozygotes become heterozygotes; nothing jumps 0 <-> 2
  expect_identical(tab["0", "2"], 0L)
  expect_identical(tab["2", "0"], 0L)
  expect_gt(tab["0", "1"] + tab["2", "1"], 0L)
  expect_identical(is.na(as.matrix(b)), is.na(as.matrix(x)))
  expect_identical(indNames(b), indNames(x))
  expect_identical(locNames(b), locNames(x))
  expect_identical(pop(b), pop(x))
  expect_true(all(ploidy(b) == 2))
})

test_that("number of mutations follows mut.rate at high rates", {
  # expected: nInd * nLoc * 2 * mut.rate, less the share at missing genotypes
  expected <- nInd(x) * nLoc(x) * 2 * 1e-3 * mean(!is.na(as.matrix(x)))
  set.seed(3)
  d <- replicate(10, n_diff(x, gl.sim.mutate(x, mut.rate = 1e-3, verbose = 0)))
  expect_equal(mean(d), expected, tolerance = 0.15)
})

test_that("[approved diff] no mutations when none are drawn (F1)", {
  set.seed(1)
  d <- replicate(10, n_diff(x, gl.sim.mutate(x, mut.rate = 0, verbose = 0)))
  expect_true(all(d == 0))
  # default rate: expected 0.13 mutations per call on this dataset
  set.seed(2)
  d <- replicate(20, n_diff(x, gl.sim.mutate(x, verbose = 0)))
  expect_lt(mean(d), 0.5)
})

test_that("[approved diff] SilicoDArT input stops (F2)", {
  gs <- gl.filter.allna(testset.gs, verbose = 0)
  expect_error(gl.sim.mutate(gs, mut.rate = 1e-3, verbose = 0))
})

test_that("[approved diff] invalid mut.rate stops with a clear error (F5)", {
  expect_error(gl.sim.mutate(x, mut.rate = 2, verbose = 0), "between 0 and 1")
  expect_error(gl.sim.mutate(x, mut.rate = -1, verbose = 0), "between 0 and 1")
  expect_error(gl.sim.mutate(x, mut.rate = "a", verbose = 0), "between 0 and 1")
})

test_that("[approved diff] flags reset and history added (F3, F6)", {
  set.seed(5)
  b <- gl.sim.mutate(x, mut.rate = 1e-3, verbose = 0)
  # existing metrics kept; reset.flags may add missing metric columns as NA
  cols <- colnames(x@other$loc.metrics)
  expect_identical(b@other$loc.metrics[, cols], x@other$loc.metrics[, cols])
  expect_false(any(unlist(b@other$loc.metrics.flags[
    c("maf", "FreqHets", "FreqHomRef", "FreqHomSnp")])))
  expect_identical(length(b@other$history), length(x@other$history) + 1L)
  expect_output(gl.sim.mutate(x, mut.rate = 1e-3, verbose = 3),
                "Applied [0-9]+ mutations")
  expect_silent(gl.sim.mutate(x, mut.rate = 1e-3, verbose = 0))
})

test_that("seeded results unchanged by reading one individual (F4)", {
  # snapshot of the pre-review output for this seed
  set.seed(11)
  b <- gl.sim.mutate(x, mut.rate = 1e-3, verbose = 0)
  expect_identical(n_diff(x, b), 122L)
})
