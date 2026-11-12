# Checks that the FBM-flavoured result b of a simulation function is a valid
# file-backed dartR object that matches the gen-flavoured result a of the same
# call with the same seed, and that its missing genotypes are stored as the
# raw FBM code 3 (decoded as NA)
expect_fbm_twin <- function(a, b) {
  expect_null(a@fbm)
  expect_s4_class(b@fbm, "FBM.code256")
  expect_length(b@gen, 0)
  expect_true(methods::validObject(b))
  expect_identical(c(nInd(b), nLoc(b)), c(nInd(a), nLoc(a)))
  expect_identical(indNames(b), indNames(a))
  expect_identical(locNames(b), locNames(a))
  expect_identical(as.character(pop(b)), as.character(pop(a)))
  expect_identical(b@loc.all, a@loc.all)
  expect_identical(sort(names(b@other)), sort(names(a@other)))
  expect_equal(unname(as.matrix(b)), unname(as.matrix(a)),
               ignore_attr = TRUE)
  raw <- b@fbm$copy(code = 0:255)[]
  expect_identical(unname(is.na(as.matrix(b))), unname(raw == 3))
}
