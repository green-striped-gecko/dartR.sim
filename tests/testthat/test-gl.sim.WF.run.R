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

# --- Dispersal: which sexes migrate (addendum F21) ---------------------------
# Populations as simulated: males in the first half, females in the second;
# V3 labels each individual with its population of origin.
disp_pop <- function(p, n = 10) {
  data.frame(V1 = rep(c("Male", "Female"), each = n / 2), V2 = p,
             V3 = paste0("p", p, "_", seq_len(n)), V4 = "x")
}
disp_table <- function(p1, p2, nt) {
  data.frame(pop1 = p1, pop2 = p2, number_transfers = nt,
             transfer_each_gen = 1)
}
# Individuals arriving in pop1 of each row, by sex, for each generation and
# row, using the dispersal step of gl.sim.WF.run() (dispersal_event()).
# Rows must not share populations, so arrivals belong to one row.
disp_moves <- function(pairs, gens) {
  pop_list <- lapply(seq_len(max(unlist(pairs[, 1:2]))), disp_pop)
  pairs$size_pop1 <- 10
  pairs$size_pop2 <- 10
  next_male <- NULL
  out <- NULL
  for (generation in gens) {
    before <- lapply(pop_list, function(p) p$V3)
    res <- dispersal_event(pop_list, pairs, generation, next_male)
    pop_list <- res[[1]]
    next_male <- res[[2]]
    for (i in seq_len(nrow(pairs))) {
      p1 <- pop_list[[pairs$pop1[i]]]
      new <- p1[!p1$V3 %in% before[[pairs$pop1[i]]], ]
      out <- rbind(out, data.frame(gen = generation, row = i,
                                   males = sum(new$V1 == "Male"),
                                   females = sum(new$V1 == "Female")))
    }
  }
  out
}

test_that("[approved diff] a row of 3 moves 2 males and 1 female (F21)", {
  set.seed(1)
  m <- disp_moves(disp_table(1, 2, 3), gens = 2:3)
  expect_identical(m$males, c(2L, 2L))
  expect_identical(m$females, c(1L, 1L))
})

test_that("[approved diff] rows move their own number whatever the others are (F21)", {
  set.seed(1)
  m <- disp_moves(disp_table(c(1, 3, 5), c(2, 4, 6), c(1, 2, 2)), gens = 2:3)
  expect_identical(m$males + m$females, c(1L, 2L, 2L, 1L, 2L, 2L))
  # rows of 2: one male and one female each
  expect_true(all(m$males[m$row > 1] == 1 & m$females[m$row > 1] == 1))
  # row of 0 moves nobody
  set.seed(1)
  m0 <- disp_moves(disp_table(c(1, 3), c(2, 4), c(0, 2)), gens = 2:3)
  expect_identical(m0$males[m0$row == 1] + m0$females[m0$row == 1], c(0L, 0L))
})

test_that("[approved diff] each pair of 1 alternates sexes over time (F21)", {
  set.seed(1)
  m <- disp_moves(disp_table(c(1, 3), c(2, 4), 1), gens = 2:5)
  expect_identical(m$males[m$row == 1], c(1L, 0L, 1L, 0L))
  expect_identical(m$males[m$row == 2], c(1L, 0L, 1L, 0L))
  expect_true(all(m$males + m$females == 1))
  # no transfer outside dispersal generations
  set.seed(1)
  t5 <- disp_table(1, 2, 2)
  t5$transfer_each_gen <- 5
  m5 <- disp_moves(t5, gens = 2:6)
  expect_identical(m5$males + m5$females, c(0L, 0L, 0L, 2L, 0L))
})

# ---- Inbreeding: sib mating and real_inbreeding ----

fis_pop <- function(g) {
  m <- as.matrix(g)
  q <- colMeans(m) / 2
  1 - sum(colMeans(m == 1)) / sum(2 * q * (1 - q))
}

test_that("sib_pairs pairs the asked proportion of full siblings", {
  # 50 families of 4 brothers and 4 sisters
  n <- 400
  fam <- rep(1:50, times = 2, each = 4)
  pop <- data.frame(V1 = rep(c("Male", "Female"), each = n / 2),
                    V5 = paste0("f", fam), V6 = paste0("m", fam))
  rownames(pop) <- paste0("r", 1:n)
  set.seed(1)
  realised <- replicate(50, {
    p <- sib_pairs(pop, n, sib = 0.3, rep_parents = FALSE)
    expect_false(anyDuplicated(p$males) > 0)
    expect_setequal(p$males, rownames(pop)[1:(n / 2)])
    mean(pop[p$males, "V5"] == pop[p$females, "V5"])
  })
  # random pairing adds 1/50 sib pairs by chance
  expect_equal(mean(realised), 0.3 + 0.7 / 50, tolerance = 0.05)
  # no families: no sib pairs and flagged short
  pop_nofam <- pop
  pop_nofam$V5 <- pop_nofam$V6 <- NA
  p <- sib_pairs(pop_nofam, n, sib = 0.3, rep_parents = FALSE)
  expect_true(p$short)
})

test_that("reproduction with sib = 0 keeps the random-mating stream", {
  rt <- wf_ref()
  gens <- function(r) lapply(r[[1]], function(g) g@gen)
  expect_identical(gens(wf_run(rt, seed = 1)),
                   gens(wf_run(rt, seed = 1, sib_mating_phase2 = 0)))
})

test_that("ibd_loci: IBD fraction is F, in stretches along the map", {
  cm <- seq(0, 10, length.out = 2000)
  set.seed(1)
  f <- replicate(500, mean(ibd_loci(cm, 0.2)))
  expect_equal(mean(f), 0.2, tolerance = 0.05)
  set.seed(1)
  runs <- rle(ibd_loci(cm, 0.2))
  # mean IBD stretch close to 25 cM = 100 loci here
  expect_gt(mean(runs$lengths[runs$values]), 40)
  expect_false(any(ibd_loci(cm, 0)))
  expect_true(all(ibd_loci(cm, 1)))
})

test_that("inbreeding_real is 1 - Ho / He per population", {
  x <- gl.keep.pop(testset.gl, pop.list = popNames(testset.gl)[1:2],
                   verbose = 0)
  f <- inbreeding_real(x)
  expect_named(f, popNames(x))
  g <- as.matrix(seppop(x)[[1]])
  n <- colSums(!is.na(g))
  q <- colMeans(g, na.rm = TRUE) / 2
  he <- 2 * q * (1 - q) * 2 * n / (2 * n - 1)
  ho <- colMeans(g == 1, na.rm = TRUE)
  expect_equal(unname(f[1]), 1 - sum(ho, na.rm = TRUE) / sum(he, na.rm = TRUE))
})

test_that("sib_mating raises F towards b / (4 - 3b)", {
  rt <- wf_ref(chunk_neutral_loci = 5)
  f <- sapply(1:2, function(i) {
    r <- wf_run(rt, seed = i, sib_mating_phase2 = 0.3,
                population_size_phase2 = "200", gen_number_phase2 = 12,
                every_gen = 12, sample_percent = 100,
                dispersal_phase2 = FALSE)
    fis_pop(r[[1]][[length(r[[1]])]])
  })
  # expected 0.3 / 3.1 = 0.097
  expect_gt(mean(f), 0.05)
  expect_lt(mean(f), 0.15)
})

test_that("sib_mating is validated", {
  rt <- wf_ref()
  expect_error(wf_run(rt, seed = 1, sib_mating_phase2 = 2), "sib_mating")
  expect_error(wf_run(rt, seed = 1, sib_mating_phase2 = "0.1 0.2"),
               "sib_mating")
  expect_error(wf_run(rt, seed = 1, real_inbreeding = TRUE), "missing")
})

test_that("real_inbreeding: population check, estimate used, set value wins", {
  x <- gl.keep.pop(testset.gl, pop.list = popNames(testset.gl)[1:2],
                   verbose = 0)
  x <- gl.filter.callrate(x, threshold = 1, verbose = 0)
  rt <- wf_ref(x = x, real_freq = TRUE)
  expect_error(wf_run(rt, x = x, seed = 1, real_freq = TRUE,
                      real_inbreeding = TRUE), "real_pops")
  expect_s4_class(
    wf_run(rt, x = x, seed = 1, real_freq = TRUE, real_pops = TRUE,
           real_inbreeding = TRUE, population_size_phase2 = "20 20",
           gen_number_phase2 = 2)[[1]][[1]],
    "genlight")
  expect_output(
    gl.sim.WF.run(file_var = fv_sim, ref_table = rt, x = x,
                  interactive_vars = FALSE, verbose = 1, seed = 1,
                  real_freq = TRUE, real_pops = TRUE, real_inbreeding = TRUE,
                  sib_mating_phase2 = 0, population_size_phase2 = "20 20",
                  gen_number_phase2 = 2),
    "is not used")
})

test_that("CSV without the inbreeding variables runs", {
  v <- read.csv(fv_sim)
  v <- v[!v$variable %in% c("sib_mating_phase1", "sib_mating_phase2",
                            "real_inbreeding"), ]
  f <- tempfile(fileext = ".csv")
  write.csv(v, f, row.names = FALSE)
  rt <- wf_ref()
  expect_identical(
    gl.sim.WF.run(file_var = f, ref_table = rt, interactive_vars = FALSE,
                  verbose = 0, seed = 1)[[1]][[2]]@gen,
    wf_run(rt, seed = 1)[[1]][[2]]@gen)
})

# ---- store_founders ----

test_that("store_founders off leaves the output unchanged", {
  rt <- wf_ref()
  a <- wf_run(rt, seed = 1)
  b <- wf_run(rt, seed = 1, store_founders = FALSE)
  expect_identical(a, b)
  expect_false("generation_0" %in% names(a[[1]]))
})

test_that("store_founders stores generation 0, parents of generation 1", {
  rt <- wf_ref()
  r <- wf_run(rt, seed = 1, store_founders = TRUE, every_gen = 1,
              sample_percent = 100, gen_number_phase2 = 2,
              population_size_phase2 = "20")
  expect_identical(names(r[[1]])[1], "generation_0")
  g0 <- r[[1]][["generation_0"]]
  g1 <- r[[1]][["generation_1"]]
  expect_identical(nInd(g0), 20L)
  expect_identical(g0@other$sim.vars$generation, 0)
  expect_true(all(is.na(g0@other$ind.metrics$pat)))
  expect_true(all(is.na(g0@other$ind.metrics$mat)))
  expect_identical(g0@other$ind.metrics$F_founder, rep(0, 20))
  expect_identical(g0@other$ind.metrics$F_founder_map, rep(0, 20))
  expect_identical(names(g0@other$ind.metrics)[3:4], c("pat", "mat"))
  # every parent of generation 1 is a stored founder
  expect_true(all(c(g1@other$ind.metrics$pat, g1@other$ind.metrics$mat) %in%
                    indNames(g0)))
  # later generations are those of a run without founders stored
  r2 <- wf_run(rt, seed = 1, every_gen = 1, sample_percent = 100,
               gen_number_phase2 = 2, population_size_phase2 = "20")
  expect_identical(r[[1]][["generation_2"]]@gen, r2[[1]][["generation_2"]]@gen)
})

test_that("store_founders with real_inbreeding: F_founder is the IBD share", {
  # platypus.gl has F > 0 in its three populations (testset.gl has F < 0)
  x <- gl.filter.callrate(platypus.gl, threshold = 0.9, verbose = 0)
  x <- gl.filter.monomorphs(x, verbose = 0)
  rt <- wf_ref(x = x, real_freq = TRUE, chunk_number = 20)
  r <- wf_run(rt, x = x, seed = 3, real_freq = TRUE, real_pops = TRUE,
              real_inbreeding = TRUE, sib_mating_phase2 = 0,
              population_size_phase2 = "200 200 200", gen_number_phase2 = 1,
              store_founders = TRUE, sample_percent = 100)
  g0 <- r[[1]][["generation_0"]]
  f <- g0@other$ind.metrics$F_founder
  expect_true(all(f >= 0 & f <= 1))
  F_x <- pmax(inbreeding_real(x), 0)
  expect_true(all(F_x > 0))
  expect_equal(as.numeric(tapply(f, pop(g0), mean)), unname(F_x),
               tolerance = 0.3)
  # homozygosity follows F_founder: IBD loci are homozygous
  expect_equal(g0@other$ind.metrics$F_founder_map, f, tolerance = 0.1)
  het <- rowMeans(as.matrix(g0) == 1)
  expect_lt(cor(f, het), -0.3)
})

# ---- real_freq_shrink ----

hudson_fst <- function(g, loci = TRUE) {
  ps <- seppop(g)
  a <- as.matrix(ps[[1]])[, loci, drop = FALSE]
  b <- as.matrix(ps[[2]])[, loci, drop = FALSE]
  n1 <- 2 * colSums(!is.na(a))
  n2 <- 2 * colSums(!is.na(b))
  p1 <- colMeans(a, na.rm = TRUE) / 2
  p2 <- colMeans(b, na.rm = TRUE) / 2
  num <- (p1 - p2)^2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
  den <- p1 * (1 - p2) + p2 * (1 - p1)
  ok <- is.finite(num) & is.finite(den)
  sum(num[ok]) / sum(den[ok])
}

test_that("shrink_freq shrinks toward the mean; auto removes sampling noise", {
  set.seed(1)
  L <- 2000
  p <- runif(L, 0.1, 0.9)
  # no true differentiation: samples of 10 individuals from the same p
  freq <- lapply(1:2, function(i) rbinom(L, 20, p) / 20)
  n <- lapply(1:2, function(i) rep(10, L))
  s <- shrink_freq(freq, n, "auto")
  expect_lt(s$lambda, 0.3)
  half <- shrink_freq(freq, n, 0.5)
  expect_equal(half$freq[[1]] - half$freq[[2]],
               0.5 * (freq[[1]] - freq[[2]]))
  expect_equal(half$freq[[1]] + half$freq[[2]], freq[[1]] + freq[[2]])
})

test_that("real_freq_shrink = auto matches founder FST to x's FST", {
  x <- gl.filter.callrate(platypus.gl, threshold = 0.9, verbose = 0)
  x <- gl.filter.monomorphs(x, verbose = 0)
  x <- gl.keep.pop(x, pop.list = popNames(x)[c(1, 3)], verbose = 0)
  rt <- wf_ref(x = x, real_freq = TRUE, chunk_neutral_loci = 0)
  fst <- function(shrink, seed) {
    r <- wf_run(rt, x = x, seed = seed, real_freq = TRUE, real_pops = TRUE,
                real_pop_size = TRUE, real_freq_shrink = shrink,
                store_founders = TRUE, sample_percent = 100,
                gen_number_phase2 = 1)
    g0 <- r[[1]][["generation_0"]]
    c(hudson_fst(g0, g0@other$loc.metrics$type == "real"),
      as.numeric(g0@other$sim.vars$freq_shrink_lambda))
  }
  auto <- rowMeans(sapply(1:4, function(s) fst("auto", s)))
  none <- mean(sapply(1:4, function(s) fst("NULL", s)[1]))
  target <- hudson_fst(x)
  expect_gt(none, target)
  expect_equal(auto[1], target, tolerance = 0.2)
  expect_gt(auto[2], 0)
  expect_lt(auto[2], 1)
  expect_error(fst("2", 1), "real_freq_shrink")
})

# ---- real_migration ----

test_that("fst_hudson matches the pairwise Hudson estimator", {
  x <- gl.filter.callrate(platypus.gl, threshold = 0.9, verbose = 0)
  x2 <- gl.keep.pop(x, pop.list = popNames(x)[1:2], verbose = 0)
  expect_equal(fst_hudson(x2), hudson_fst(x2))
})

test_that("real_migration sets T = (1/FST - 1) / (4n) and holds FST", {
  x <- gl.filter.callrate(platypus.gl, threshold = 0.9, verbose = 0)
  x <- gl.filter.monomorphs(x, verbose = 0)
  rt <- wf_ref(x = x, real_freq = TRUE, chunk_neutral_loci = 0)
  run <- function(seed, ...) {
    wf_run(rt, x = x, seed = seed, real_freq = TRUE, real_pops = TRUE,
           real_freq_shrink = "auto", dispersal_phase2 = TRUE,
           population_size_phase2 = "100 100 100", sample_percent = 100,
           gen_number_phase2 = 20, every_gen = 19, ...)
  }
  r <- run(1, real_migration = TRUE)
  g <- r[[1]][[length(r[[1]])]]
  target <- fst_hudson(x)
  expect_equal(as.numeric(g@other$sim.vars$migrants_real),
               (1 / target - 1) / (4 * 3), tolerance = 1e-4)
  f <- mean(sapply(1:3, function(s) {
    g <- run(s, real_migration = TRUE)[[1]][[2]]
    fst_hudson(g[, g@other$loc.metrics$type == "real"])
  }))
  expect_equal(f, target, tolerance = 0.3)
  expect_error(run(1, real_migration = TRUE, dispersal_type_phase2 = "line"),
               "all_connected")
  expect_error(wf_run(wf_ref(), seed = 1, real_migration = TRUE), "missing")
})

# ---- Effective population size ----

test_that("ne_ratio and k_for_ne are inverse; maxima 1 and 1/2", {
  for (rp in c(FALSE, TRUE)) {
    k <- c(0.1, 1, 10)
    expect_equal(k_for_ne(ne_ratio(k, rp), rp), k)
  }
  expect_equal(ne_ratio(1e6, TRUE), 0.5, tolerance = 1e-5)
  expect_true(is.infinite(k_for_ne(0.6, TRUE)))
})

test_that("ne_phase2 sets variance_offspring and the stored Ne", {
  rt <- wf_ref()
  sv <- function(r) r[[1]][[1]]@other$sim.vars
  r <- wf_run(rt, seed = 1, population_size_phase2 = "100 100",
              number_pops_phase2 = 2, ne_phase2 = "30 60")
  expect_identical(sv(r)$ne_expected, "30 60")
  expect_equal(as.numeric(strsplit(sv(r)$variance_offspring_used, " ")[[1]]),
               c(30 / 70, 60 / 40))
  expect_identical(sv(r)$population_size_used, "100 100")
  # with replacement Ne is at most N/2: capped, with a warning
  expect_output(
    r2 <- gl.sim.WF.run(file_var = fv_sim, ref_table = rt,
                        interactive_vars = FALSE, verbose = 1, seed = 1,
                        population_size_phase2 = "100",
                        replace_parents = TRUE, ne_phase2 = 80),
    "above the largest Ne")
  expect_equal(as.numeric(sv(r2)$ne_expected), 50, tolerance = 0.01)
  expect_error(wf_run(rt, seed = 1, ne_phase2 = "-5"), "ne_phase2")
  expect_error(wf_run(rt, seed = 1, variance_offspring_phase2 = "1 2 3"),
               "variance_offspring_phase2")
})

test_that("ne_phase2: heterozygosity follows gl.diagnostics.sim's Ne", {
  rt <- wf_ref(chunk_neutral_loci = 5)
  ratio <- sapply(1:3, function(s) {
    r <- wf_run(rt, seed = s, population_size_phase2 = "100",
                replace_parents = TRUE, ne_phase2 = 30, sample_percent = 100,
                gen_number_phase2 = 21, every_gen = 20,
                dispersal_phase2 = FALSE)
    # a single population: FST cannot be computed, so compare He directly
    he <- function(g) mean(gl.He(g))
    g <- r[[1]]
    c(he(g[[2]]) / he(g[[1]]))
  })
  expected <- (1 - 1 / 60)^20
  expect_equal(mean(ratio), expected, tolerance = 0.05)
})

test_that("real_sample_size stores x's sample sizes; real_migration scales by N/Ne", {
  x <- gl.filter.callrate(platypus.gl, threshold = 0.9, verbose = 0)
  x <- gl.filter.monomorphs(x, verbose = 0)
  rt <- wf_ref(x = x, real_freq = TRUE, chunk_neutral_loci = 0)
  run <- function(sizes = "100 100 100", ...) {
    wf_run(rt, x = x, seed = 1, real_freq = TRUE, real_pops = TRUE,
           population_size_phase2 = sizes, gen_number_phase2 = 2,
           dispersal_phase2 = TRUE, real_migration = TRUE, ...)
  }
  r <- run(real_sample_size = TRUE, ne_phase2 = "25 50 50",
           store_founders = TRUE)
  expect_equal(as.vector(table(pop(r[[1]][["generation_0"]]))),
               c(24L, 18L, 42L))
  expect_equal(as.vector(table(pop(r[[1]][[2]]))), c(24L, 18L, 42L))
  base <- (1 / fst_hudson(x) - 1) / 12
  expect_equal(as.numeric(r[[1]][[2]]@other$sim.vars$migrants_real),
               base * mean(100 / c(25, 50, 50)))
  # generation 0 holds the phase-2 rate too (read by gl.diagnostics.sim)
  expect_equal(
    as.numeric(r[[1]][["generation_0"]]@other$sim.vars$number_transfers_phase2),
    base * mean(100 / c(25, 50, 50)))
  expect_error(run(sizes = "20 20 20", real_sample_size = TRUE), "larger")
})

test_that("gl.diagnostics.sim takes Ne from sim.vars by default", {
  rt <- wf_ref()
  r <- wf_run(rt, seed = 1, number_pops_phase2 = 2,
              population_size_phase2 = "40 40", ne_phase2 = 30,
              gen_number_phase2 = 5, every_gen = 2)
  pdf(NULL)
  on.exit(dev.off())
  d1 <- gl.diagnostics.sim(r, verbose = 0)
  d2 <- gl.diagnostics.sim(r, Ne = 30, verbose = 0)
  expect_equal(d1$he, d2$he)
  expect_equal(d1$fst, d2$fst)
})

# ---- store_pedigree ----

test_that("store_pedigree returns every individual, sampled or not", {
  rt <- wf_ref()
  r <- wf_run(rt, seed = 1, number_pops_phase2 = 2,
              population_size_phase2 = "20 30", gen_number_phase2 = 3,
              every_gen = 1, sample_percent = 50, replace_parents = TRUE,
              store_founders = TRUE, store_pedigree = TRUE)
  ped <- attr(r[[1]], "pedigree")
  expect_named(ped, c("id", "pat", "mat", "generation", "pop", "F_founder"))
  # 4 generations (0-3) of 50 individuals
  expect_identical(nrow(ped), 200L)
  expect_false(anyDuplicated(ped$id) > 0)
  expect_true(all(is.na(ped$pat[ped$generation == 0])))
  expect_true(all(ped$F_founder[ped$generation == 0] == 0))
  expect_true(all(is.na(ped$F_founder[ped$generation > 0])))
  # every parent is in the pedigree, one generation earlier
  kids <- ped[ped$generation > 0, ]
  gen_of <- setNames(ped$generation, ped$id)
  expect_true(all(gen_of[kids$pat] == kids$generation - 1))
  expect_true(all(gen_of[kids$mat] == kids$generation - 1))
  # stored (sampled) individuals are a subset, with matching parents
  for (g in r[[1]]) {
    m <- match(indNames(g), ped$id)
    expect_false(anyNA(m))
    expect_identical(as.character(g@other$ind.metrics$pat), ped$pat[m])
  }
  # off by default: no attribute, output unchanged
  expect_null(attr(wf_run(rt, seed = 1)[[1]], "pedigree"))
})

# ---- inbreeding_founders ----

test_that("inbreeding_founders sets the founders' F, with or without x", {
  rt <- wf_ref(chunk_number = 20)
  f0 <- function(r) {
    g0 <- r[[1]][["generation_0"]]
    tapply(g0@other$ind.metrics$F_founder, pop(g0), mean)
  }
  # no genlight needed
  r <- wf_run(rt, seed = 1, number_pops_phase2 = 2,
              population_size_phase2 = "200 200", gen_number_phase2 = 1,
              inbreeding_founders = "0.3 0", sib_mating_phase2 = 0,
              store_founders = TRUE, sample_percent = 100)
  expect_equal(unname(f0(r)[1]), 0.3, tolerance = 0.2)
  expect_equal(unname(f0(r)[2]), 0)
  # replaces the estimate from x, and sets the sib-mating rate
  x <- gl.filter.callrate(platypus.gl, threshold = 0.9, verbose = 0)
  x <- gl.filter.monomorphs(x, verbose = 0)
  rtx <- wf_ref(x = x, real_freq = TRUE, chunk_number = 20)
  expect_output(
    rx <- gl.sim.WF.run(file_var = fv_sim, ref_table = rtx, x = x,
                        interactive_vars = FALSE, verbose = 1, seed = 1,
                        real_freq = TRUE, real_pops = TRUE,
                        real_inbreeding = TRUE, inbreeding_founders = -0.1,
                        population_size_phase2 = "50 50 50",
                        gen_number_phase2 = 1, store_founders = TRUE,
                        sample_percent = 100),
    "not estimated")
  expect_true(all(rx[[1]][["generation_0"]]@other$ind.metrics$F_founder == 0))
  expect_error(wf_run(rt, seed = 1, inbreeding_founders = "0.1 0.2 0.3"),
               "inbreeding_founders")
})
