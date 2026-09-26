# tests/testthat/test-review2-fold-separation.R
# ---------------------------------------------------------------------------
# Second review round: fold_separation().  Each test fails on the code before
# the fix.  Helpers are in helper-review2-folds.R.
# ---------------------------------------------------------------------------

# ---- fold_separation() -----------------------------------------------------

test_that("fold_separation labels folds by the fold_id a cv_*() result carries", {
  set.seed(1)
  d <- r2_pts(runif(150, 0, 1000), runif(150, 0, 1000), a = rnorm(150))
  d$z <- d$a + rnorm(150)
  f <- make_folds(d, k = 5, method = "block_kfold", seed = 1)
  d$z[f$folds[[1]]$test] <- NA               # fold 1 has nothing to score
  cv <- suppressWarnings(r2_quiet(cv_spatial(d, "z", "a", fit_fn = r2_fit, folds = f)))
  s <- fold_separation(cv$folds, d)
  expect_identical(s$fold, as.integer(cv$fold_metrics$fold))   # 2 3 4 5, not 1 2 3 4
  m <- merge(as.data.frame(cv$fold_metrics), as.data.frame(s), by = "fold")
  expect_identical(m$n_test.x, m$n_test.y)
})

test_that("fold_separation measures a recorded range in the CRS it was recorded in", {
  set.seed(1)
  d <- r2_pts(7e5 + runif(80, 0, 2000), 3.95e6 + runif(80, 0, 2000), crs = 32617)
  f <- make_folds(d, k = 4, method = "block_kfold", seed = 1)
  f$params$sac_range <- structure(600, class = "sac_range", crs = sf::st_crs(32617))
  a <- fold_separation(f, d)
  b <- fold_separation(f, sf::st_transform(d, 2264))   # the same layer in US feet
  expect_equal(b$within_range, a$within_range, tolerance = 1e-6)
  expect_equal(attr(b, "crs"), "EPSG:32617")
  # The folds' own CRS is used when the range carries none.
  f$params$sac_range <- 600
  b2 <- fold_separation(f, sf::st_transform(d, 2264))
  expect_equal(b2$within_range, a$within_range, tolerance = 1e-6)
  # A bare number stays in data_sf's own units, as documented.
  b3 <- fold_separation(f$folds, sf::st_transform(d, 2264), sac = 600)
  expect_equal(attr(b3, "crs"), "EPSG:2264")
  expect_output(print(a), "in EPSG:32617 units")
  expect_error(fold_separation(f, d, sac = units::set_units(0.6, km)),
               "`sac` must be a plain number")
})
