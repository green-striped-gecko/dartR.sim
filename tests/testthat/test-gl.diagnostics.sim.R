# Characterization tests for gl.diagnostics.sim
# Baseline snapshotted before review (review-gl.diagnostics.sim). Assertions
# marked [approved diff] were flipped in Phase C under the finding named in
# the test title; see function-review/reports/dartR.sim/.

fv_ref <- system.file("extdata", "ref_variables.csv", package = "dartR.sim")
fv_sim <- system.file("extdata", "sim_variables.csv", package = "dartR.sim")
rt_d <- gl.sim.WF.table(file_var = fv_ref, interactive_vars = FALSE,
                        verbose = 0, seed = 1)
res_d <- gl.sim.WF.run(file_var = fv_sim, ref_table = rt_d,
                       interactive_vars = FALSE, verbose = 0, seed = 1,
                       number_pops_phase2 = 2,
                       population_size_phase2 = "10 10")
diag_quiet <- function(...) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  suppressWarnings(suppressMessages(gl.diagnostics.sim(..., verbose = 0)))
}

test_that("[approved diff] returns plot and tables (F5)", {
  res <- diag_quiet(res_d, Ne = 10)
  expect_named(res, c("plot", "he", "fst"))
  expect_s3_class(res$plot, "patchwork")
  expect_named(res$fst, c("gen", "observed", "expected"))
  expect_true(all(c("gen", "observed", "Ne", "expected") %in% names(res$he)))
})

test_that("[approved diff] expected FST follows the island recursion (F1)", {
  r <- gl.sim.WF.run(file_var = fv_sim, ref_table = rt_d,
                     interactive_vars = FALSE, verbose = 0, seed = 1,
                     number_pops_phase2 = 2, population_size_phase2 = "10 10",
                     gen_number_phase2 = 60, every_gen = 10)
  fst <- diag_quiet(r, Ne = 10)$fst
  # the recursion starts from the observed value
  expect_equal(fst$expected[1], fst$observed[1])
  # equilibrium close to 1 / (1 + 4 Ne m n / (n - 1)), Ne = 10, m = 0.1
  expect_equal(tail(fst$expected, 1), 1 / (1 + 4 * 10 * 0.1 * 2),
               tolerance = 0.1)
  line <- gl.sim.WF.run(file_var = fv_sim, ref_table = rt_d,
                        interactive_vars = FALSE, verbose = 0, seed = 1,
                        number_pops_phase2 = 3,
                        population_size_phase2 = "10 10 10",
                        dispersal_type_phase2 = "line")
  expect_error(diag_quiet(line, Ne = 10), "all_connected")
})

test_that("[approved diff] population sizes quoted in the CSV are read (F3)", {
  v <- read.csv(fv_sim)
  v$value[v$variable == "number_pops_phase2"] <- "2"
  v$value[v$variable == "population_size_phase2"] <- '"10 10"'
  f <- tempfile(fileext = ".csv")
  write.csv(v, f, row.names = FALSE)
  r <- gl.sim.WF.run(file_var = f, ref_table = rt_d,
                     interactive_vars = FALSE, verbose = 0, seed = 1)
  expect_false(anyNA(diag_quiet(r, Ne = 10)$fst$expected))
})

test_that("[approved diff] expected He starts at the first stored generation (F2)", {
  r <- gl.sim.WF.run(file_var = fv_sim, ref_table = rt_d,
                     interactive_vars = FALSE, verbose = 0, seed = 1,
                     number_pops_phase2 = 2, number_pops_phase1 = 2,
                     population_size_phase1 = "10 10",
                     population_size_phase2 = "10 10", phase1 = TRUE,
                     every_gen = 2, sample_percent = 100)
  he <- diag_quiet(r, Ne = 10)$he
  first <- he[he$gen == min(he$gen), ]
  expect_equal(min(he$gen), 11)
  expect_equal(first$expected, first$observed)
})

test_that("[approved diff] invalid inputs stop with clear errors (F4)", {
  one_pop <- gl.sim.WF.run(file_var = fv_sim, ref_table = rt_d,
                           interactive_vars = FALSE, verbose = 0, seed = 1)
  expect_error(diag_quiet(one_pop, Ne = 50), "two populations")
  expect_error(diag_quiet(res_d, Ne = 10, iteration = 2), "iteration")
  expect_error(diag_quiet(res_d, Ne = 10, pops_fst = c(1, 3)), "pops_fst")
})

test_that("[approved diff] verbose = 0 prints nothing (F6)", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  out <- utils::capture.output(
    suppressWarnings(gl.diagnostics.sim(res_d, Ne = 10, verbose = 0)),
    type = "message")
  expect_length(out, 0)
})
