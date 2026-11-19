# Characterization tests for gl.sim.emigration
# Baseline snapshotted before review (review-gl.sim.emigration), bugs included.
# Assertions marked [approved diff] were flipped in Phase C under the finding
# named in the test title; see function-review/reports/dartR.sim/.
# Needs dartR.base with the rbind.dartR fix (dartR.base PR #422).

x2 <- gl.filter.allna(
  gl.keep.pop(testset.gl, c("EmmacMDBCond", "EmmacMDBCudg"), verbose = 0),
  verbose = 0
)
x3 <- gl.filter.allna(
  gl.keep.pop(testset.gl, c("EmmacMDBCond", "EmmacMDBCudg", "EmmacMDBForb"),
              verbose = 0),
  verbose = 0
)
# Cross-table of original population (rows) against population after
# emigration (columns), in the input's population order
moves <- function(before, after) {
  orig <- setNames(as.character(pop(before)), indNames(before))
  table(factor(orig[indNames(after)], levels = popNames(before)),
        factor(as.character(pop(after)), levels = popNames(before)))
}

test_that("emi.table moves individuals from column to row", {
  t <- matrix(0, 2, 2)
  t[2, 1] <- 5 # 5 from EmmacMDBCond to EmmacMDBCudg
  set.seed(1)
  r <- gl.sim.emigration(x2, emi.table = t, verbose = 0)
  expect_identical(as.vector(moves(x2, r)), c(5L, 0L, 5L, 10L))
  expect_identical(nInd(r), nInd(x2))
  expect_identical(nLoc(r), nLoc(x2))
  expect_identical(levels(pop(r)), levels(pop(x2)))
})

test_that("genotypes and individual metadata stay with individuals", {
  set.seed(1)
  r <- gl.sim.emigration(x3, emi.table = matrix(1, 3, 3), verbose = 0)
  expect_identical(as.matrix(r)[indNames(x3), ], as.matrix(x3))
  expect_identical(nrow(r@other$ind.metrics), nInd(r))
  expect_identical(r@other$ind.metrics$id, indNames(r))
  m <- match(indNames(r), indNames(x3))
  expect_equal(unname(as.matrix(r@other$latlon)),
               unname(as.matrix(x3@other$latlon[m, ])))
  expect_identical(nrow(r@other$loc.metrics), nLoc(x3))
})

test_that("[approved diff] emi.m moves individuals from column to row (F1)", {
  # emi.m: emigrants of EmmacMDBCond go to EmmacMDBCudg; emigrants of
  # EmmacMDBCudg stay (diagonal)
  e <- matrix(c(0, 1, 0, 1), 2, 2)
  set.seed(1)
  r <- gl.sim.emigration(x2, perc.mig = 0.5, emi.m = e, verbose = 0)
  mv <- moves(x2, r)
  expect_gt(mv[1, 2], 0L) # EmmacMDBCond individuals moved
  expect_identical(mv[2, 1], 0L) # nobody from EmmacMDBCudg moved
  # expected fraction moving follows perc.mig
  set.seed(3)
  f <- replicate(300, {
    r <- gl.sim.emigration(x2, perc.mig = 0.2,
                           emi.m = matrix(c(0, 1, 1, 0), 2, 2), verbose = 0)
    mv <- moves(x2, r)
    (mv[1, 2] + mv[2, 1]) / nInd(x2)
  })
  expect_equal(mean(f), 0.2, tolerance = 0.05)
})

test_that("[approved diff] ind.metrics$pop follows pop() (F4)", {
  t <- matrix(0, 2, 2)
  t[2, 1] <- 5
  set.seed(1)
  r <- gl.sim.emigration(x2, emi.table = t, verbose = 0)
  expect_identical(r@other$ind.metrics$pop, as.character(pop(r)))
})

test_that("[approved diff] individuals move at most once (F2)", {
  t <- matrix(0, 3, 3)
  t[2, 1] <- 5 # Cond -> Cudg
  t[3, 2] <- 5 # Cudg -> Forb
  set.seed(1)
  r <- gl.sim.emigration(x3, emi.table = t, verbose = 0)
  mv <- moves(x3, r)
  expect_identical(mv[1, 3], 0L) # no Cond individual ended in Forb
  expect_identical(c(mv[1, 2], mv[2, 3]), c(5L, 5L))
  # output grouped by population
  expect_false(is.unsorted(as.integer(pop(r))))
})

test_that("[approved diff] populations can empty (F3)", {
  t <- matrix(0, 2, 2)
  t[2, 1] <- 10
  r <- gl.sim.emigration(x2, emi.table = t, verbose = 0)
  expect_identical(popNames(r), "EmmacMDBCudg")
  expect_identical(nInd(r), 20L)
  expect_output(gl.sim.emigration(x2, emi.table = t, verbose = 1),
                "no individuals left: EmmacMDBCond")
  set.seed(1)
  expect_no_error(gl.sim.emigration(x2, perc.mig = 0.99,
                                    emi.m = matrix(c(0, 1, 1, 0), 2, 2),
                                    verbose = 0))
  # list output keeps an empty element, and accepts it back
  rl <- gl.sim.emigration(seppop(x2), emi.table = t, verbose = 0)
  expect_null(rl$EmmacMDBCond)
  back <- matrix(0, 2, 2)
  back[1, 2] <- 3
  rl2 <- gl.sim.emigration(rl, emi.table = back, verbose = 0)
  expect_identical(unname(sapply(rl2, nInd)), c(3L, 17L))
})

test_that("[approved diff] input errors (F3, F5)", {
  t <- matrix(0, 2, 2)
  t[2, 1] <- 1000
  expect_error(gl.sim.emigration(x2, emi.table = t, verbose = 0),
               "more emigrants than individuals")
  expect_error(gl.sim.emigration(x2, verbose = 0), "Provide either")
  expect_error(gl.sim.emigration(x2, emi.table = matrix(1, 3, 3),
                                 verbose = 0), "square matrix")
  expect_error(gl.sim.emigration(x2, perc.mig = 0.2, emi.m = matrix(1, 3, 2),
                                 verbose = 0), "square matrix")
  expect_error(gl.sim.emigration(x2, perc.mig = 10,
                                 emi.m = matrix(c(0, 1, 1, 0), 2, 2),
                                 verbose = 0), "between 0 and 1")
  expect_error(gl.sim.emigration(x2, emi.table = matrix(0.5, 2, 2),
                                 verbose = 0), "whole numbers")
  expect_error(gl.sim.emigration(list(x2, "a"), emi.table = matrix(0, 2, 2),
                                 verbose = 0), "genlight")
  # a zero column in emi.m: nobody leaves that population
  set.seed(1)
  r <- gl.sim.emigration(x2, perc.mig = 0.5, emi.m = matrix(0, 2, 2),
                         verbose = 0)
  expect_identical(sum(moves(x2, r)) - sum(diag(moves(x2, r))), 0L)
  # a list of length 1 is a list
  r1 <- gl.sim.emigration(list(x2), emi.table = matrix(0, 1, 1), verbose = 0)
  expect_type(r1, "list")
})

test_that("baseline: data.frame matrices (dartR GUI) are accepted", {
  t <- matrix(0, 2, 2)
  t[2, 1] <- 5
  expect_identical(nInd(gl.sim.emigration(x2, emi.table = as.data.frame(t),
                                          verbose = 0)),
                   nInd(x2))
  set.seed(1)
  r <- gl.sim.emigration(x2, perc.mig = 0.2,
                         emi.m = as.data.frame(matrix(c(0, 1, 1, 0), 2, 2)),
                         verbose = 0)
  expect_identical(nInd(r), nInd(x2))
})

test_that("[approved diff] list input (F5)", {
  t <- matrix(0, 2, 2)
  t[2, 1] <- 5
  set.seed(1)
  rl <- gl.sim.emigration(seppop(x2), emi.table = t, verbose = 0)
  expect_type(rl, "list")
  expect_identical(names(rl), popNames(x2))
  expect_identical(unname(sapply(rl, nInd)), c(5L, 15L))
  expect_identical(levels(pop(rl[[2]])), "EmmacMDBCudg")
  # unnamed list: elements are named after their population
  set.seed(1)
  ru <- gl.sim.emigration(unname(seppop(x2)), emi.table = t, verbose = 0)
  expect_identical(names(ru), popNames(x2))
  expect_identical(levels(pop(ru[[1]])), "EmmacMDBCond")
})

test_that("[approved diff] history and messages (F6)", {
  set.seed(1)
  r <- gl.sim.emigration(x2, emi.table = matrix(c(0, 1, 1, 0), 2, 2),
                         verbose = 0)
  calls <- vapply(r@other$history, function(h) deparse(h)[1], character(1))
  expect_identical(sum(grepl("gl.sim.emigration", calls)), 1L)
  expect_identical(length(r@other$history), length(x2@other$history) + 1L)
  expect_output(gl.sim.emigration(x2, emi.table = matrix(c(0, 1, 1, 0), 2, 2),
                                  verbose = 3),
                "Moved 2 individuals between 2 populations")
  expect_silent(gl.sim.emigration(x2, emi.table = matrix(c(0, 1, 1, 0), 2, 2),
                                  verbose = 0))
})
