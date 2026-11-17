# File-backed objects: gl.He, gl.Ho, gl.report.nall and the per-population
# helpers read the genotypes in blocks of loci
# (dartR.base::utils.fbm.colwise). Results and messages must be identical to
# those of the gen-backed object for any block size, and no call may decode
# the whole FBM.

fbm_pair <- function() {
  x <- gl.compliance.check(testset.gl[1:40, 1:80], verbose = 0)
  m <- as.matrix(x)
  m[, 5] <- NA                               # an all-NA locus
  m[pop(x) == levels(pop(x))[1], 7] <- NA    # all-NA in one population
  x@gen <- new("genlight", m, ploidy = ploidy(x))@gen
  list(x = x, f = suppressWarnings(gl.gen2fbm(x, verbose = 0)))
}

calls <- list(
  He = list(gl.He),
  Ho = list(gl.Ho),
  alf = list(utils.sim.alf),
  fst = list(fst_hudson),
  inbreeding = list(inbreeding_real),
  nall = list(gl.report.nall, simlevels = c(2, 6), reps = 2, ncores = 1,
              plot.display = FALSE, verbose = 3)
)

run <- function(FUN, x, ...) {
  set.seed(1)
  out <- utils::capture.output(val <- suppressWarnings(FUN(x, ...)))
  if (is.list(val) && !is.data.frame(val)) {
    val[vapply(val, inherits, logical(1), "ggplot")] <- NULL
  }
  list(val = val, out = out)
}

test_that("gen- and file-backed objects give identical results", {
  skip_if_not_installed("bigstatsr")
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  old <- options(bigstatsr.block.sizeGB = getOption("bigstatsr.block.sizeGB"))
  on.exit(options(old), add = TRUE)
  xs <- fbm_pair()
  for (bs in c(1, 1e-6)) {                   # 1e-6: one locus per block
    options(bigstatsr.block.sizeGB = bs)
    for (nm in names(calls)) {
      cl <- calls[[nm]]
      g <- do.call(run, c(list(cl[[1]], xs$x), cl[-1]))
      f <- do.call(run, c(list(cl[[1]], xs$f), cl[-1]))
      expect_identical(f$val, g$val, label = paste(nm, bs))
      expect_identical(f$out, g$out, label = paste(nm, bs, "output"))
    }
  }
})

test_that("file-backed objects are never decoded as a whole", {
  skip_if_not_installed("bigstatsr")
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  xs <- fbm_pair()
  local_mocked_bindings(.fbm_as_matrix = function(x) {
    stop("the whole FBM was decoded")
  }, .package = "dartR.base")
  for (nm in names(calls)) {
    cl <- calls[[nm]]
    set.seed(1)
    expect_no_error(utils::capture.output(suppressWarnings(
      do.call(cl[[1]], c(list(xs$f), cl[-1])))), message = nm)
  }
})
