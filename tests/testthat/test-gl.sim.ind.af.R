# Characterization tests for gl.sim.ind.af
# Baseline snapshotted before review (review-gl.sim.ind.af), bugs included.
# Assertions marked [approved diff] were flipped in Phase C under the finding
# named in the test title; see function-review/reports/dartR.sim/.

# Two populations with contrasting frequencies, listed Z before A
# (order of appearance differs from alphabetical order).
n_loc <- 1000L
two_pops <- data.frame(
  popn = rep(c("Z", "A"), each = n_loc),
  locus = rep(paste0("l", seq_len(n_loc)), 2),
  frequency = rep(c(0.9, 0.1), each = n_loc)
)
pop_freq <- function(g) {
  tapply(rowMeans(as.matrix(g)) / 2, as.character(pop(g)), mean)
}

test_that("documented example: structure and allele frequencies", {
  t1 <- gl.filter.callrate(platypus.gl, threshold = 1, mono.rm = TRUE,
                           verbose = 0)
  r1 <- gl.allele.freq(t1, by = "popxloc", verbose = 0)
  r2 <- r1[, c("popn", "locus", "frequency")]
  set.seed(1)
  res <- gl.sim.ind.af(df = r2, pop.sizes = c(50, 50, 50), verbose = 0)
  expect_s4_class(res, "genlight")
  expect_s4_class(res, "dartR")
  expect_identical(c(nInd(res), nLoc(res)), c(150L, nLoc(t1)))
  expect_true(all(ploidy(res) == 2))
  expect_identical(locNames(res), locNames(t1))
  expect_identical(as.vector(table(pop(res))), c(50L, 50L, 50L))
  expect_identical(indNames(res)[1], "0_SEVERN_ABOVE_1")
  # simulated frequencies track the input per population
  p_in <- r2$frequency[r2$popn == "SEVERN_ABOVE"]
  p_out <- colMeans(as.matrix(res[pop(res) == "SEVERN_ABOVE"])) / 2
  expect_lt(mean(abs(p_out - p_in)), 0.06)
})

test_that("Hardy-Weinberg genotype proportions", {
  d <- data.frame(popn = "p", locus = paste0("l", 1:500), frequency = 0.3)
  set.seed(2)
  m <- as.matrix(gl.sim.ind.af(d, pop.sizes = 2000, verbose = 0))
  expect_equal(c(mean(m == 0), mean(m == 1), mean(m == 2)),
               c(0.49, 0.42, 0.09), tolerance = 0.02)
})

test_that("[approved diff] populations get their own frequencies (F1)", {
  set.seed(1)
  g <- gl.sim.ind.af(two_pops, pop.sizes = c(Z = 100, A = 10), verbose = 0)
  expect_identical(as.vector(table(pop(g))[c("Z", "A")]), c(100L, 10L))
  f <- pop_freq(g)
  expect_gt(f[["Z"]], 0.8)
  expect_lt(f[["A"]], 0.2)
  # unnamed sizes follow order of appearance
  g <- gl.sim.ind.af(two_pops, pop.sizes = c(100, 10), verbose = 0)
  expect_identical(as.vector(table(pop(g))[c("Z", "A")]), c(100L, 10L))
  expect_gt(pop_freq(g)[["Z"]], 0.8)
})

test_that("[approved diff] factor popn with named pop.sizes matches by name (F2)", {
  d <- two_pops
  d$popn <- factor(d$popn, levels = c("Z", "A"))
  set.seed(1)
  g <- gl.sim.ind.af(d, pop.sizes = c(A = 10, Z = 100), verbose = 0)
  expect_identical(as.vector(table(pop(g))[c("Z", "A")]), c(100L, 10L))
  f <- pop_freq(g)
  expect_gt(f[["Z"]], 0.8)
  expect_lt(f[["A"]], 0.2)
})

test_that("[approved diff] unused factor levels are ignored (F2)", {
  d <- two_pops
  d$popn <- factor(d$popn, levels = c("A", "Z", "Q"))
  g <- gl.sim.ind.af(d, pop.sizes = c(10, 100), verbose = 0)
  expect_identical(as.vector(table(pop(g))[c("Z", "A")]), c(10L, 100L))
})

test_that("[approved diff] pop.sizes checks (F5)", {
  d <- two_pops
  expect_error(gl.sim.ind.af(d, pop.sizes = c(2.7, 3.9), verbose = 0),
               "whole numbers")
  expect_error(gl.sim.ind.af(d, pop.sizes = c(NA, 3), verbose = 0),
               "whole numbers")
  expect_error(gl.sim.ind.af(d, pop.sizes = c(0, 3), verbose = 0),
               "whole numbers")
  expect_error(gl.sim.ind.af(d, pop.sizes = c(3, 3, 3), verbose = 0),
               "must equal")
  expect_error(gl.sim.ind.af(d, pop.sizes = c(Z = 3, B = 3), verbose = 0),
               "Missing sizes")
  expect_error(gl.sim.ind.af(d[, 1:2], pop.sizes = c(3, 3), verbose = 0),
               "three columns")
  d$frequency[1] <- 1.2
  expect_error(gl.sim.ind.af(d, pop.sizes = c(3, 3), verbose = 0),
               "\\[0, 1\\]")
})

test_that("[approved diff] duplicated locus rows stop (F5)", {
  d <- rbind(two_pops, data.frame(popn = "Z", locus = "l1", frequency = 0))
  expect_error(gl.sim.ind.af(d, pop.sizes = c(3, 3), verbose = 0),
               "Z-l1")
})

test_that("[approved diff] metadata of the returned object (F4, F6)", {
  set.seed(1)
  g <- gl.sim.ind.af(two_pops, pop.sizes = c(3, 2), verbose = 0)
  expect_s4_class(g, "dartR")
  expect_identical(colnames(g@other$ind.metrics),
                   c("id", "pop", "fid", "iid", "sex", "phenotype"))
  expect_identical(g@other$ind.metrics$id, indNames(g))
  expect_identical(g@other$ind.metrics$pop, as.character(pop(g)))
  expect_identical(as.character(g@other$ind.metrics$sex),
                   c("m", "f", "m", "m", "f"))
  expect_identical(indNames(g)[4], "0_A_1")
  expect_identical(nrow(g@other$loc.metrics), n_loc)
  expect_identical(g@other$loc.metrics$AlleleID, locNames(g))
  expect_false("array(NA, nLoc(x))" %in% colnames(g@other$loc.metrics))
  expect_true(all(is.na(g@other$loc.metrics$CallRate)))
  expect_identical(unique(g@loc.all), "G/C")
  expect_length(g@other$history, 1)
  expect_null(g@chromosome)
})

test_that("[approved diff] verbose messages (F6)", {
  expect_output(gl.sim.ind.af(two_pops, pop.sizes = c(3, 2), verbose = 3),
                "Simulated 5 individuals in 2 populations at 1000 loci")
  expect_silent(gl.sim.ind.af(two_pops, pop.sizes = c(3, 2), verbose = 0))
})
