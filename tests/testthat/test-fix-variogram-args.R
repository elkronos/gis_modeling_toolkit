# tests/testthat/test-fix-variogram-args.R
# ---------------------------------------------------------------------------
# estimate_sac_range() checks cutoff, range_frac, n_max, seed and
# keep_directional_fits, and the arguments added with them.
#
# A bad value used to come back as an unrelated refusal (cutoff = 0, -0.5, NA
# and c(0.2, 0.5) all returned NA with "no variogram model could be fitted
# (singular fits)"; range_frac = -1 returned "fitted range exceeds the largest
# lag fitted") or as a bare R error that named no argument ("missing value
# where TRUE/FALSE needed" for range_frac = NA and n_max = NA, "non-numeric
# argument to binary operator" for cutoff = "half").
# ---------------------------------------------------------------------------

fva_field <- function(n = 120, seed = 1) {
  set.seed(seed)
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000)
  D <- as.matrix(stats::dist(cbind(x, y)))
  z <- as.numeric(t(chol(0.8 * exp(-D / 100) + diag(0.2, n))) %*% rnorm(n))
  sf::st_as_sf(data.frame(x = 5e5 + x, y = 5e6 + y, z = z),
               coords = c("x", "y"), crs = 32632)
}

test_that("a bad cutoff is an error that names the argument", {
  skip_if_not_installed("gstat")
  pts <- fva_field()
  for (bad in list(0, -0.5, NA, NA_real_, "half", c(0.2, 0.5), 3, Inf, NULL, TRUE))
    expect_error(estimate_sac_range(pts, "z", cutoff = bad),
                 "^estimate_sac_range\\(\\): `cutoff` must be a single number above 0 and at most 1",
                 info = paste(deparse(bad), collapse = " "))
  # The message shows what was passed.
  expect_error(estimate_sac_range(pts, "z", cutoff = 3), "got 3\\.$")
  expect_error(estimate_sac_range(pts, "z", cutoff = c(0.2, 0.5)), "got numeric of length 2\\.$")
  expect_error(estimate_sac_range(pts, "z", cutoff = "half"), "got character \\(half\\)\\.$")
})

test_that("a bad range_frac, n_max, seed or keep_directional_fits is an error that names it", {
  skip_if_not_installed("gstat")
  pts <- fva_field()
  for (bad in list(NA, -1, 0, Inf, "1", c(1, 2)))
    expect_error(estimate_sac_range(pts, "z", range_frac = bad),
                 "^estimate_sac_range\\(\\): `range_frac` must be a single finite number above 0",
                 info = paste(deparse(bad), collapse = " "))
  for (bad in list(NA, 10, 29.9, "many", c(100, 200), -Inf))
    expect_error(estimate_sac_range(pts, "z", n_max = bad),
                 "^estimate_sac_range\\(\\): `n_max` must be a single number of at least 30",
                 info = paste(deparse(bad), collapse = " "))
  for (bad in list(NA, "y", 1, c(TRUE, FALSE), NULL))
    expect_error(estimate_sac_range(pts, "z", keep_directional_fits = bad),
                 "^estimate_sac_range\\(\\): `keep_directional_fits` must be TRUE or FALSE",
                 info = paste(deparse(bad), collapse = " "))
  for (bad in list("a", NA, c(1, 2), Inf))
    expect_error(estimate_sac_range(pts, "z", n_max = 60, seed = bad),
                 "^estimate_sac_range\\(\\): `seed` must be NULL or a single whole number",
                 info = paste(deparse(bad), collapse = " "))
})

test_that("the new arguments are checked the same way", {
  skip_if_not_installed("gstat")
  pts <- fva_field()
  expect_error(estimate_sac_range(pts, "z", lags = "log"), "'arg' should be one of")
  for (bad in list(4, 101, 7.5, NA, "15", c(10, 15)))
    expect_error(estimate_sac_range(pts, "z", n_lags = bad),
                 "^estimate_sac_range\\(\\): `n_lags` must be a single whole number from 5 to 100",
                 info = paste(deparse(bad), collapse = " "))
  for (bad in list(0, -0.1, 1.5, NA, "0.05", c(0.05, 0.1)))
    expect_error(estimate_sac_range(pts, "z", structure_alpha = bad),
                 "^estimate_sac_range\\(\\): `structure_alpha` must be a single number above 0 and at most 1",
                 info = paste(deparse(bad), collapse = " "))
})

test_that("the values at the edges of what is allowed still run", {
  skip_if_not_installed("gstat")
  pts <- fva_field()
  base <- estimate_sac_range(pts, "z")
  expect_s3_class(base, "sac_range")
  # No subsample is drawn at any of these, so the answer is the default one.
  expect_identical(estimate_sac_range(pts, "z", n_max = Inf), base)
  expect_identical(estimate_sac_range(pts, "z", n_max = 5000.7), base)
  expect_identical(estimate_sac_range(pts, "z", seed = NULL), base)
  # A fraction is floored: 60.9 points is 60 points.
  r60 <- estimate_sac_range(pts, "z", n_max = 60.9)
  expect_identical(attr(r60, "n_used"), 60L)
  expect_identical(attr(r60, "n_subsampled"), 60L)
  expect_s3_class(estimate_sac_range(pts, "z", cutoff = 1), "sac_range")
  expect_s3_class(estimate_sac_range(pts, "z", n_lags = 5, lags = "equal"), "sac_range")
  expect_s3_class(estimate_sac_range(pts, "z", structure_alpha = 1), "sac_range")
})
