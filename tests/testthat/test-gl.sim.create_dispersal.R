# Characterization tests for gl.sim.create_dispersal
# Baseline snapshotted before review (review-gl.sim.create_dispersal), bugs
# included. Assertions marked [approved diff] were flipped in Phase C under
# the finding named in the test title; see function-review/reports/dartR.sim/.

td <- tempfile("dispersal")
dir.create(td)
# Write the table to td and read it back
disp <- function(...) {
  gl.sim.create_dispersal(..., outpath = td, verbose = 0)
  utils::read.csv(file.path(td, "dispersal_table.csv"))
}
pairs_of <- function(d) paste(d$pop1, d$pop2, sep = "-")

test_that("each connected pair is written once, per dispersal type", {
  expect_identical(pairs_of(disp(number_pops = 4)),
                   c("2-1", "3-1", "3-2", "4-1", "4-2", "4-3"))
  expect_identical(pairs_of(disp(number_pops = 4, dispersal_type = "line")),
                   c("1-2", "2-3", "3-4"))
  expect_identical(pairs_of(disp(number_pops = 4, dispersal_type = "circle")),
                   c("1-2", "2-3", "3-4", "4-1"))
  expect_identical(pairs_of(disp(number_pops = 2, dispersal_type = "circle")),
                   "1-2")
})

test_that("file columns and values as read by gl.sim.WF.run", {
  d <- disp(number_pops = 3, dispersal_type = "line", number_transfers = 2,
            transfer_each_gen = 5)
  expect_identical(colnames(d),
                   c("pop1", "pop2", "number_transfers", "transfer_each_gen"))
  expect_true(all(d$number_transfers == 2))
  expect_true(all(d$transfer_each_gen == 5))
  expect_true(file.exists(file.path(td, "dispersal_table.csv")))
})

test_that("[approved diff] invalid arguments stop (F1)", {
  expect_error(disp(number_pops = 1), "at least 2")
  expect_error(disp(number_pops = 1, dispersal_type = "circle"), "at least 2")
  expect_error(disp(number_pops = 0), "at least 2")
  expect_error(disp(number_pops = 2.5), "at least 2")
  expect_error(disp(number_pops = 3, dispersal_type = "ring"),
               "must be one of")
  expect_error(disp(number_pops = 3, transfer_each_gen = 0), "at least 1")
  expect_error(disp(number_pops = 3, number_transfers = -1), "at least 0")
  expect_error(disp(number_pops = 3, number_transfers = 1.5), "at least 0")
  # 0 transfers is allowed (switches a pair off)
  expect_true(all(disp(number_pops = 3, number_transfers = 0)$
                    number_transfers == 0))
})

test_that("[approved diff] missing outpath falls back to tempdir() (F2)", {
  f <- file.path(tempdir(), "disp_fallback.csv")
  unlink(f)
  expect_output(gl.sim.create_dispersal(3, outpath = file.path(td, "nope"),
                                        outfile = "disp_fallback.csv",
                                        verbose = 1),
                "does not exist")
  expect_true(file.exists(f))
})

test_that("[approved diff] returns the table invisibly; cat messages (F3, F4)", {
  expect_invisible(gl.sim.create_dispersal(3, outpath = td, verbose = 0))
  r <- gl.sim.create_dispersal(3, outpath = td, verbose = 0)
  expect_equal(r, disp(number_pops = 3))
  out <- capture.output(gl.sim.create_dispersal(3, outpath = td, verbose = 2))
  expect_true(any(grepl("dispersal_table.csv\\s*$", out)))
  expect_false(any(grepl("dispersal_table.csv/", out)))
  expect_silent(gl.sim.create_dispersal(3, outpath = td, verbose = 0))
})
