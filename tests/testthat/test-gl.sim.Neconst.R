# Characterization tests for gl.sim.Neconst
# Baseline snapshotted before review (review-gl.sim.Neconst), bugs included.
# Assertions marked [approved diff] were flipped in Phase C under the finding
# named in the test title; see function-review/reports/dartR.sim/.

# Unfolded site frequency spectrum: loci by count of the counted allele
sfs <- function(g) {
  k <- colSums(as.matrix(g))
  as.vector(table(factor(k, levels = 0:(2 * nInd(g)))))
}
n_ind <- 20L

test_that("output structure", {
  set.seed(1)
  g <- gl.sim.Neconst(ninds = n_ind, nlocs = 500, verbose = 0)
  expect_s4_class(g, "dartR")
  expect_identical(c(nInd(g), nLoc(g)), c(n_ind, 500L))
  expect_true(all(ploidy(g) == 2))
  expect_false(anyNA(as.matrix(g)))
  expect_identical(nPop(g), 1L)
})

test_that("[approved diff] every locus is polymorphic and the spectrum is neutral (F1)", {
  set.seed(1)
  g <- gl.sim.Neconst(ninds = n_ind, nlocs = 20000, verbose = 0)
  s <- sfs(g)
  expect_identical(s[1] + s[length(s)], 0L)
  poly <- s[2:(2 * n_ind)]
  i <- seq_len(2 * n_ind - 1)
  ex <- 1 / i + 1 / (2 * n_ind - i)
  ex <- ex / sum(ex) * sum(poly)
  # singletons and the middle of the spectrum within 5% of expectation
  expect_equal(poly[1] + poly[2 * n_ind - 1], ex[1] + ex[2 * n_ind - 1],
               tolerance = 0.05)
  mid <- 5:35
  expect_equal(sum(poly[mid]), sum(ex[mid]), tolerance = 0.05)
  # Hardy-Weinberg within the population (sampling without replacement)
  p <- colSums(as.matrix(g)) / (2 * n_ind)
  expect_equal(mean(as.matrix(g) == 1),
               mean(2 * p * (1 - p) * 2 * n_ind / (2 * n_ind - 1)),
               tolerance = 0.02)
})

test_that("mutation_rate only matters when 4 * ninds * mutation_rate is large", {
  he <- function(g) {
    p <- colMeans(as.matrix(g)) / 2
    mean(2 * p * (1 - p))
  }
  set.seed(1)
  a <- he(gl.sim.Neconst(n_ind, 5000, mutation_rate = 1e-8, verbose = 0))
  set.seed(1)
  b <- he(gl.sim.Neconst(n_ind, 5000, mutation_rate = 1e-5, verbose = 0))
  set.seed(1)
  c <- he(gl.sim.Neconst(n_ind, 5000, mutation_rate = 0.05, verbose = 0))
  expect_equal(a, b, tolerance = 0.01)
  # theta = 4: expected heterozygosity of Beta(4, 4) is 4/9
  expect_equal(c, 4 / 9, tolerance = 0.02)
})

test_that("[approved diff] inputs are checked (F3)", {
  expect_error(gl.sim.Neconst(10, 100, mutation_rate = 0, verbose = 0),
               "greater than 0")
  expect_error(gl.sim.Neconst(10, 100, mutation_rate = -1, verbose = 0),
               "greater than 0")
  expect_error(gl.sim.Neconst(0, 100, verbose = 0), "at least 2")
  expect_error(gl.sim.Neconst(10.5, 100, verbose = 0), "at least 2")
  expect_error(gl.sim.Neconst(10, 100.7, verbose = 0), "at least 1")
})

test_that("[approved diff] metadata and history (F4, F5)", {
  set.seed(1)
  g <- gl.sim.Neconst(ninds = 10, nlocs = 50, verbose = 0)
  expect_identical(indNames(g)[1:2], c("1", "2"))
  expect_identical(locNames(g)[1:2], c("Loc1", "Loc2"))
  expect_identical(popNames(g), "pop1")
  expect_identical(unique(g@loc.all), "A/C")
  expect_identical(colnames(g@other$ind.metrics), c("id", "pop"))
  expect_identical(g@other$ind.metrics$pop, rep("pop1", 10))
  expect_identical(g@other$loc.metrics$AlleleID, locNames(g))
  expect_false("array(NA, nLoc(x))" %in% colnames(g@other$loc.metrics))
  expect_length(g@other$history, 1)
  expect_match(deparse(g@other$history[[1]])[1], "gl.sim.Neconst")
  expect_output(gl.sim.Neconst(ninds = 10, nlocs = 50, verbose = 3),
                "Simulated 10 individuals at 50 polymorphic loci")
  expect_silent(gl.sim.Neconst(ninds = 10, nlocs = 50, verbose = 0))
  # a single locus whose genotypes are all 0 or 1 is still SNP data
  set.seed(2)
  expect_no_error(gl.sim.Neconst(3, 1, verbose = 0))
})
