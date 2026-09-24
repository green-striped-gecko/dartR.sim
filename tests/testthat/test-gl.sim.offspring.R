# Characterization tests for gl.sim.offspring
# Baseline snapshotted before review (review-gl.sim.offspring), bugs included.
# Assertions marked [approved diff] were flipped in Phase C under the finding
# named in the test title; see function-review/reports/dartR.sim/.

mk_parent <- function(g) {
  new("dartR", gen = matrix(g, nrow = 1), ploidy = 2, ind.names = "p",
      loc.names = paste0("L", seq_along(g)),
      loc.all = rep("A/G", length(g)))
}
x_off <- gl.impute(gl.filter.allna(testset.gl, verbose = 0),
                   method = "random", verbose = 0)

test_that("Mendelian proportions and sex ratio", {
  het <- mk_parent(rep(1L, 20))
  set.seed(4)
  o <- gl.sim.offspring(het, het, 4000, verbose = 0)
  p <- prop.table(table(factor(as.matrix(o), levels = 0:2)))
  expect_equal(as.vector(p), c(0.25, 0.5, 0.25), tolerance = 0.02)
  expect_true(all(ploidy(o) == 2))
  expect_identical(locNames(o), locNames(het))
  s <- gl.sim.offspring(x_off[6:10, ], x_off[1:5, ], 200, sexratio = 0.3,
                        verbose = 0)
  expect_equal(mean(s@other$sex == "female"), 0.3, tolerance = 0.05)
  expect_identical(nInd(s), 1000L)
})

test_that("[approved diff] loci are inherited independently (F1)", {
  g <- rep(0L, 20)
  g[c(1, 11)] <- 1L
  set.seed(1)
  m <- as.matrix(gl.sim.offspring(mk_parent(rep(0L, 20)), mk_parent(g), 2000,
                                  verbose = 0))
  expect_equal(mean(m[, 1] == m[, 11]), 0.5, tolerance = 0.05)
  # real mother: no pair of heterozygous loci perfectly linked
  mo <- x_off[1, ]
  fa <- x_off[2, ]
  gm <- as.matrix(x_off)[1:2, ]
  loci <- which(gm[1, ] == 1 & gm[2, ] != 1)
  set.seed(5)
  o <- as.matrix(gl.sim.offspring(fa, mo, 1000, verbose = 0))[, loci]
  r <- cor(sweep(o, 2, gm[2, loci] / 2))
  expect_equal(sum(abs(r[upper.tri(r)]) > 0.99), 0)
})

test_that("[approved diff] each mother has noffpermother offspring (F2)", {
  o <- gl.sim.offspring(x_off[6:10, ], x_off[1:5, ], 4, verbose = 0)
  expect_equal(as.vector(table(o@other$ind.metrics$mother)), rep(4, 5))
  expect_true(all(o@other$ind.metrics$father %in% indNames(x_off)[6:10]))
})

test_that("[approved diff] metadata and parentage are carried (F4)", {
  o <- gl.sim.offspring(x_off[6:10, ], x_off[1:5, ], 2, verbose = 0)
  expect_s4_class(o, "dartR")
  expect_identical(o@loc.all, x_off@loc.all)
  expect_identical(position(o), position(x_off))
  expect_equal(nrow(o@other$loc.metrics), nLoc(x_off))
  expect_named(o@other$ind.metrics, c("id", "pop", "sex", "mother", "father"))
  expect_identical(o@other$ind.metrics$sex, o@other$sex)
  expect_length(o@other$history, 1)
  expect_equal(nInd(rbind(x_off[1:3, ], o[1:2, ])), 5)
})

test_that("[approved diff] inputs are checked (F3)", {
  fa <- x_off[6:10, sample(nLoc(x_off))]
  expect_error(gl.sim.offspring(fa, x_off[1:5, ], 2, verbose = 0),
               "same loci")
  expect_error(gl.sim.offspring(testset.gs[1:3, ], testset.gs[4:6, ], 2,
                                verbose = 0))
  expect_error(gl.sim.offspring(x_off[6:10, ], x_off[1:5, ], 0, verbose = 0),
               "noffpermother")
  expect_error(gl.sim.offspring(x_off[6:10, ], x_off[1:5, ], 2,
                                sexratio = 2, verbose = 0), "sexratio")
})
