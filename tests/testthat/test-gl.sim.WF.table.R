# Characterization tests for gl.sim.WF.table
# Baseline snapshotted before review (review-gl.sim.WF.table), bugs included.
# Assertions marked [approved diff] were flipped in Phase C under the finding
# named in the test title; see function-review/reports/dartR.sim/.

fv <- system.file("extdata", "ref_variables.csv", package = "dartR.sim")
wf_tab <- function(...) {
  gl.sim.WF.table(file_var = fv, interactive_vars = FALSE, verbose = 0, ...)
}
chr1_gl <- function() {
  x <- testset.gl[, 1:200]
  set.seed(3)
  x@chromosome <- factor(rep("1", nLoc(x)))
  x@position <- sort(sample(1:1e7, nLoc(x)))
  x
}

test_that("default CSV gives 100 neutral loci with fixed layout", {
  res <- wf_tab(seed = 1)
  ref <- res$reference
  expect_named(res, c("reference", "ref_vars"))
  expect_identical(dim(ref), c(100L, 8L))
  expect_identical(colnames(ref),
                   c("q", "h", "s", "c", "loc_bp", "loc_cM", "chr_name", "type"))
  expect_true(all(ref$type == "neutral"))
  expect_identical(unique(ref$chr_name), "1")
  expect_equal(ref$loc_bp[1:3], c(50002, 150002, 250002))
  expect_equal(max(ref$loc_cM), 9.9)
  expect_equal(ref$c[100], 0)
  expect_identical(res, wf_tab(seed = 1))
})

test_that("deleterious, advantageous and mutation loci counts", {
  ref <- wf_tab(seed = 1, chromosome_name = "1", loci_deleterious = 500,
                loci_advantageous = 50, loci_mut_neu = 50, loci_mut_del = 50,
                loci_mut_adv = 50)$reference
  tab <- table(ref$type)
  expect_equal(as.vector(tab[c("deleterious", "advantageous", "neutral",
                               "mutation_neu", "mutation_del",
                               "mutation_adv")]),
               c(500, 50, 100, 50, 50, 50))
  expect_true(all(ref$q <= 0.5))
  expect_true(all(ref$s[ref$type == "advantageous"] >= -0.5))
})

test_that("[approved diff] ... override keeps chromosome_name clean (F1)", {
  ref <- wf_tab(seed = 1, chunk_number = 20)$reference
  expect_identical(unique(ref$chr_name), "1")
  x <- chr1_gl()
  ref <- wf_tab(seed = 1, x = x, real_loc = TRUE)$reference
  expect_equal(sum(ref$type == "real"), 200)
})

test_that("[approved diff] unknown ... argument stops (F7)", {
  expect_error(wf_tab(seed = 1, chunk_numbr = 20), "chunk_numbr")
  expect_error(gl.sim.WF.table(interactive_vars = FALSE, verbose = 0),
               "file_var")
  expect_error(wf_tab(x = 1:10), "genlight")
})

test_that("real_loc = TRUE uses genlight positions", {
  x <- chr1_gl()
  ref <- wf_tab(seed = 1, x = x, real_loc = TRUE,
                chromosome_name = "1")$reference
  expect_equal(as.vector(table(ref$type)[c("neutral", "real")]), c(100, 200))
})

test_that("[approved diff] real_freq = TRUE with real_loc = FALSE (F2)", {
  x <- chr1_gl()
  ref <- wf_tab(seed = 1, x = x, real_freq = TRUE)$reference
  expect_equal(sum(ref$type == "real"), 200)
  expect_false(anyNA(ref[, c("q", "h", "s", "type")]))
})

test_that("[approved diff] caps follow each class's own setting (F5)", {
  ref <- wf_tab(seed = 1, loci_advantageous = 200,
                s_distribution_del = "equal", exp_rate = 1)$reference
  expect_gte(min(ref$s[ref$type == "advantageous"]), -0.5)
  ref <- wf_tab(seed = 1, loci_advantageous = 50,
                q_distribution_adv = "equal", q_adv = 0.8)$reference
  expect_true(all(ref$q[ref$type == "advantageous"] == 0.8))
})

test_that("[approved diff] map intervals of any size are placed correctly (F3)", {
  map <- data.frame(Chr = "1", from = seq(1, 2e7, 2e5),
                    to = seq(2e5, 2e7, 2e5),
                    cM = c(rep(0, 50), rep(4, 50)))
  mf <- tempfile(fileext = ".csv")
  write.csv(map, mf, row.names = FALSE)
  ref <- wf_tab(seed = 1, file_r_map = mf)$reference
  # first 10 Mb carry 0 cM in the map
  expect_equal(ref$loc_cM[max(which(ref$loc_bp < 1e7))], 0)
  expect_gt(max(ref$loc_bp), 1.9e7)
  map$from <- seq(1, 5e6, 5e4)
  map$to <- seq(5e4, 5e6, 5e4)
  map$cM[100] <- NA
  write.csv(map, mf, row.names = FALSE)
  expect_s3_class(wf_tab(seed = 1, file_r_map = mf)$reference, "data.frame")
})

test_that("[approved diff] loci_deleterious total equals the request (F4)", {
  n_del <- function(n) {
    ref <- wf_tab(seed = 1, loci_deleterious = n)$reference
    sum(ref$type == "deleterious")
  }
  expect_identical(vapply(c(99, 149, 150, 250), n_del, numeric(1)),
                   c(99, 149, 150, 250))
})

test_that("fly example files ship and run", {
  fm <- system.file("extdata", "fly_recom_map.csv", package = "dartR.sim")
  ft <- system.file("extdata", "fly_targets_of_selection.csv",
                    package = "dartR.sim")
  ref <- wf_tab(seed = 1, chromosome_name = "2L", file_r_map = fm,
                file_targets_sel = ft)$reference
  expect_true(all(c("neutral", "deleterious") %in% ref$type))
})

test_that("errors carry their message (F6)", {
  expect_error(wf_tab(real_loc = TRUE), "real dataset")
})
