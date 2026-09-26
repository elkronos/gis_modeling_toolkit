# tests/testthat/test-review-resolution.R
# ---------------------------------------------------------------------------
# Regressions from the review of the resolution step: resolution_profile()
# and determine_optimal_levels().
# ---------------------------------------------------------------------------

# A hand-made sac_range, so the variogram-based columns exist without gstat
# and without the time an estimate takes.
rr_sac <- function(pts, range = 1500, nugget = 1, psill = 1) {
  vm <- data.frame(model = c("Nug", "Exp"), psill = c(nugget, psill),
                   range = c(0, range / 3), stringsAsFactors = FALSE)
  structure(range, class = c("sac_range", "numeric"), variogram_model = vm,
            crs = sf::st_crs(pts))
}


test_that("one missing response or predictor value does not switch the profile to the raw response", {
  # The response is the predictor plus white noise, and the predictor carries
  # a strong east-west trend: the residuals have no structure, the raw
  # response a great deal.  One NA used to make the first lm.fit() on every
  # row fail, and the profile then scored the raw response.
  set.seed(10)
  n <- 300
  d <- data.frame(x = runif(n, 0, 10000), y = runif(n, 0, 10000))
  d$p <- d$x / 1000 + rnorm(n)
  d$z <- 3 * d$p + rnorm(n)
  pts <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
  sac <- rr_sac(pts)
  clean <- resolution_profile(pts, "z", "p", levels = c(10, 20, 40), sac = sac)
  expect_identical(attr(clean, "variable"), "residuals")
  for (col in c("z", "p")) {
    holed <- pts
    holed[[col]][17] <- NA
    lines <- capture_spatialkit_log(
      prof <- resolution_profile(holed, "z", "p", levels = c(10, 20, 40), sac = sac))
    expect_identical(attr(prof, "variable"), "residuals", info = col)
    expect_true(log_has(lines, "1 of 300 row"), info = col)
    expect_false(log_has(lines, "OLS fit on `predictor_vars` failed"), info = col)
    # One row fewer of the same residuals: the RSS barely moves.  The raw
    # response's RSS was several times larger.
    expect_equal(prof$rss, clean$rss, tolerance = 0.05, info = col)
    expect_equal(prof$cp, clean$cp, tolerance = 0.05, info = col)
  }
})
