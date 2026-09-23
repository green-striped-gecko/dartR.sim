# Characterization tests for gl.sim.WF.run
# Baseline snapshotted before review (review-gl.sim.WF.run), bugs included.
# Assertions marked [approved diff] were flipped in Phase C under the finding
# named in the test title; see function-review/reports/dartR.sim/.

fv_ref <- system.file("extdata", "ref_variables.csv", package = "dartR.sim")
fv_sim <- system.file("extdata", "sim_variables.csv", package = "dartR.sim")
wf_ref <- function(...) {
  gl.sim.WF.table(file_var = fv_ref, interactive_vars = FALSE, verbose = 0,
                  seed = 1, ...)
}
wf_run <- function(ref_table, ...) {
  gl.sim.WF.run(file_var = fv_sim, ref_table = ref_table,
                interactive_vars = FALSE, verbose = 0, ...)
}
true_gen <- function(res_it) {
  unname(sapply(res_it, function(g) g@other$sim.vars$generation))
}

test_that("default example: structure and seeded output", {
  res <- wf_run(wf_ref(), seed = 1)
  expect_named(res, "iteration_1")
  expect_named(res[[1]], c("generation_1", "generation_10"))
  g <- res[[1]][["generation_10"]]
  expect_s4_class(g, "genlight")
  expect_identical(c(nInd(g), nLoc(g)), c(26L, 100L))
  expect_true(all(ploidy(g) == 2))
  expect_identical(nrow(g@other$loc.metrics), nLoc(g))
  expect_identical(nrow(g@other$ind.metrics), nInd(g))
  expect_equal(true_gen(res[[1]]), c(1, 10))
  expect_identical(res, wf_run(wf_ref(), seed = 1))
})

test_that("drift: neutral He decays at rate consistent with Ne = N", {
  rt <- wf_ref(chunk_neutral_loci = 5)
  he <- function(g) {
    p <- colMeans(as.matrix(g)) / 2
    mean(2 * p * (1 - p))
  }
  h <- sapply(1:5, function(i) {
    r <- wf_run(rt, seed = 100 + i, every_gen = 30, sample_percent = 100,
                gen_number_phase2 = 31, population_size_phase2 = "50",
                dispersal_phase2 = FALSE)
    c(he(r[[1]][[1]]), he(r[[1]][[2]]))
  })
  ratio <- mean(h[2, ]) / mean(h[1, ])
  # expected (1 - 1/100)^30 = 0.74
  expect_gt(ratio, 0.65)
  expect_lt(ratio, 0.82)
})

test_that("[approved diff] crossovers follow the map, siblings independent (F1)", {
  rt <- wf_ref(chunk_number = 10, chunk_cM = 50, chunk_neutral_loci = 10)
  ref <- rt$reference
  L <- nrow(ref)
  rmap <- ref[, c("c", "loc_bp", "loc_cM")]
  ev <- ceiling(sum(rmap$c))
  rmap[L + 1, 1] <- ev - sum(rmap$c)
  rmap[L + 1, 2:3] <- rmap[L, 2:3]
  n <- 200
  pop <- data.frame(V1 = rep(c("Male", "Female"), each = n / 2), V2 = 1,
                    V3 = strrep("0", L), V4 = strrep("1", L),
                    id = paste0("i", 1:n))
  set.seed(2)
  off <- reproduction(pop, 1, n, var_off = 1e6, num_off = 10, r_event = ev,
                      recom = TRUE, r_males = TRUE, r_map_1 = rmap,
                      n_loc = L, gen = 1, rep_parents = FALSE)
  xo <- sapply(off$V3, function(h) sum(diff(utf8ToInt(h)) != 0))
  sib <- ave(seq_along(off$V5), off$V5, FUN = seq_along)
  # visible crossovers between adjacent loci: sum of Haldane r = 4.74
  expected <- sum(0.5 * (1 - exp(-2 * rmap$c[1:L])))
  expect_equal(mean(xo), expected, tolerance = 0.05)
  expect_equal(mean(xo[sib == 1]), mean(xo[sib == 10]), tolerance = 0.1)
})

test_that("[approved diff] advantageous selection acts when local_adap NULL (F2)", {
  rt <- wf_ref(chunk_number = 10, chunk_neutral_loci = 2,
               loci_advantageous = 40, s_distribution_adv = "equal",
               s_adv = 0.1, q_distribution_adv = "equal", q_adv = 0.1,
               h_distribution_adv = "equal", h_adv = 0.5)
  r <- wf_run(rt, seed = 3, sample_percent = 100, selection_phase2 = TRUE,
              population_size_phase2 = "500", gen_number_phase2 = 20)
  g <- r[[1]][[length(r[[1]])]]
  adv <- rt$reference$type == "advantageous"
  q_adv <- mean(colMeans(as.matrix(g)[, adv]) / 2)
  # deterministic expectation after 20 generations is 0.226
  expect_gt(q_adv, 0.18)
})

test_that("[approved diff] clinal adaptation along part of the populations (F3)", {
  rt <- wf_ref(loci_advantageous = 20)
  expect_s4_class(
    wf_run(rt, seed = 1, number_pops_phase2 = 3,
           population_size_phase2 = "50 50 50", selection_phase2 = TRUE,
           clinal_adap = "2 3")[[1]][[1]],
    "genlight")
})

test_that("[approved diff] generation labels match stored generation (F4)", {
  rt <- wf_ref()
  r <- wf_run(rt, seed = 1, phase1 = TRUE, every_gen = 5)
  expect_named(r[[1]], c("generation_11", "generation_16", "generation_20"))
  expect_equal(true_gen(r[[1]]), c(11, 16, 20))
  r <- wf_run(rt, seed = 1, number_iterations = 2)
  expect_named(r[[2]], c("generation_1", "generation_10"))
  expect_equal(true_gen(r[[2]]), c(1, 10))
})

test_that("[approved diff] extinction keeps stored generations (F5)", {
  r <- wf_run(wf_ref(), seed = 1, number_offspring_phase2 = 1)
  expect_length(r[[1]], 0)
  expect_output(
    gl.sim.WF.run(file_var = fv_sim, ref_table = wf_ref(),
                  interactive_vars = FALSE, verbose = 1, seed = 2,
                  number_offspring_phase2 = 2.2, every_gen = 2,
                  gen_number_phase2 = 20, sample_percent = 100),
    "extinct")
})

test_that("[approved diff] only mutation loci mutate; empty pool is skipped (F6)", {
  rt <- wf_ref(q_neutral = 0)
  r <- wf_run(rt, seed = 1, mutation = TRUE, mut_rate = 1,
              sample_percent = 100)
  g <- r[[1]][[length(r[[1]])]]
  expect_equal(sum(colSums(as.matrix(g)) > 0), 0)
  rt <- wf_ref(q_neutral = 0, loci_mut_neu = 20)
  r <- wf_run(rt, seed = 1, mutation = TRUE, mut_rate = 1,
              sample_percent = 100)
  g <- r[[1]][[length(r[[1]])]]
  poly <- colSums(as.matrix(g)) > 0
  expect_true(all(g@other$loc.metrics$type[poly] == "mutation_neu"))
})

test_that("[approved diff] real_freq uses the alternative allele (F7)", {
  x <- gl.filter.callrate(testset.gl, threshold = 1, verbose = 0)
  x <- gl.keep.pop(x, pop.list = popNames(x)[1], verbose = 0)
  rt <- wf_ref(x = x, real_freq = TRUE)
  r <- wf_run(rt, x = x, seed = 1, real_freq = TRUE, real_pops = TRUE,
              sample_percent = 100, population_size_phase2 = "200",
              gen_number_phase2 = 1)
  g <- r[[1]][[1]]
  real <- g@other$loc.metrics$type == "real"
  sim_alt <- colMeans(as.matrix(g)[, real]) / 2
  real_alt <- colMeans(as.matrix(x), na.rm = TRUE) / 2
  expect_gt(cor(real_alt, sim_alt), 0.9)
  # loci without calls no longer stop the run
  x2 <- gl.keep.pop(testset.gl, pop.list = popNames(testset.gl)[1],
                    verbose = 0)[, 1:200]
  x2@chromosome <- factor(rep("1", nLoc(x2)))
  rt2 <- wf_ref(x = x2, real_freq = TRUE)
  expect_s4_class(
    wf_run(rt2, x = x2, seed = 1, real_freq = TRUE, real_pops = TRUE,
           gen_number_phase2 = 1)[[1]][[1]],
    "genlight")
})

test_that("[approved diff] CSV without replace_parents runs (F8)", {
  v <- read.csv(fv_sim)
  v <- v[v$variable != "replace_parents", ]
  f <- tempfile(fileext = ".csv")
  write.csv(v, f, row.names = FALSE)
  expect_s4_class(
    gl.sim.WF.run(file_var = f, ref_table = wf_ref(),
                  interactive_vars = FALSE, verbose = 0,
                  seed = 1)[[1]][[1]],
    "genlight")
})

test_that("[approved diff] phase-2 founders are not cloned (F9)", {
  r <- wf_run(wf_ref(), seed = 1, phase1 = TRUE, every_gen = 1,
              sample_percent = 100, dispersal_phase2 = FALSE,
              population_size_phase2 = "100")
  im <- r[[1]][[1]]@other$ind.metrics
  expect_equal(sum(tapply(im$mat, im$pat, function(v) length(unique(v))) > 1),
               0)
})

test_that("[approved diff] each connected pair migrates once (F10)", {
  r <- wf_run(wf_ref(), seed = 1, number_pops_phase2 = 3,
              population_size_phase2 = "40 40 40", every_gen = 1,
              sample_percent = 100, gen_number_phase2 = 11)
  immigrants <- sapply(r[[1]][-1], function(g) {
    im <- g@other$ind.metrics
    birth <- sub("^[0-9]+_([0-9]+)_.*$", "\\1", c(im$pat, im$mat))
    here <- rep(as.character(pop(g)), 2)
    parents <- unique(c(im$pat, im$mat)[birth != here])
    length(parents) / nPop(g)
  })
  # 2 per population per generation (was 4 before the fix)
  expect_lt(mean(immigrants), 3)
})

test_that("[approved diff] weak relative selection is reported (F11)", {
  expect_output(
    gl.sim.WF.run(file_var = fv_sim, ref_table = wf_ref(),
                  interactive_vars = FALSE, verbose = 1, seed = 1,
                  selection_phase2 = TRUE, number_offspring_phase2 = 3),
    "offspring pool")
})

test_that("[approved diff] ... overrides: lists work, unknown names stop (F12)", {
  expect_s4_class(
    wf_run(wf_ref(loci_advantageous = 20), seed = 1,
           number_pops_phase2 = 3, population_size_phase2 = "50 50 50",
           selection_phase2 = TRUE, local_adap = "1 2")[[1]][[1]],
    "genlight")
  expect_error(wf_run(wf_ref(), seed = 1, gen_number_phase = 3),
               "gen_number_phase")
})

test_that("errors carry their message (F13)", {
  expect_error(wf_run(wf_ref(), seed = 1, real_freq = TRUE), "real_freq")
})

test_that("[approved diff] phase 1 with real populations and sizes (F15)", {
  x <- gl.keep.pop(testset.gl, pop.list = popNames(testset.gl)[1:3],
                   verbose = 0)
  r <- wf_run(wf_ref(), x = x, seed = 1, phase1 = TRUE, real_pops = TRUE,
              real_pop_size = TRUE, population_size_phase2 = "20 20 20",
              number_offspring_phase1 = 20,
              number_pops_phase2 = 3)
  expect_equal(nPop(r[[1]][[1]]), 3)
})

test_that("[addendum] mutation pool is reset for each iteration (F18)", {
  # nearly every offspring becomes a parent, so each iteration's first
  # generation keeps most of the 5 new mutations if the pool is full
  rt <- wf_ref(loci_mut_neu = 5)
  r <- wf_run(rt, seed = 1, mutation = TRUE, mut_rate = 1,
              sample_percent = 100, number_iterations = 6, every_gen = 1,
              gen_number_phase2 = 1, population_size_phase2 = "500",
              number_offspring_phase2 = 2.2, dispersal_phase2 = FALSE)
  mut <- rt$reference$type == "mutation_neu"
  n_poly <- sapply(r, function(it) {
    sum(colSums(as.matrix(it[[1]])[, mut]) > 0)
  })
  # 4.4 with the reset; 2.0 without it (pool left over from iteration 1)
  expect_gt(mean(n_poly[-1]), 3.5)
})

test_that("[addendum] gl.sim.create_dispersal writes each pair once (F19)", {
  f <- "disp.csv"
  gl.sim.create_dispersal(number_pops = 3, outpath = tempdir(),
                          outfile = f, verbose = 0)
  d <- read.csv(file.path(tempdir(), f))
  expect_equal(nrow(d), 3)
  expect_false(anyDuplicated(paste(pmin(d$pop1, d$pop2),
                                   pmax(d$pop1, d$pop2))) > 0)
})

test_that("population labels are right with 10 or more populations (F20)", {
  r <- wf_run(wf_ref(), seed = 1, number_pops_phase2 = 12,
              population_size_phase2 = paste(c(10, rep(20, 11)), collapse = " "),
              dispersal_phase2 = FALSE, sample_percent = 100, every_gen = 1,
              gen_number_phase2 = 2)
  g <- r[[1]][[2]]
  born <- sapply(strsplit(g@other$ind.metrics$pat, "_"), "[", 2)
  # without dispersal, every parent was born in its offspring's population
  expect_identical(as.character(pop(g)), born)
  expect_identical(popNames(g), as.character(1:12))
  expect_equal(as.vector(table(pop(g))), c(10, rep(20, 11)))
})
