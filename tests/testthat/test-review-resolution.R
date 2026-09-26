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


test_that("a layer larger than sample_n is bounded, judged and scored on all its points", {
  # 1200 points on a 10 km square, fitted on a 300-point subsample.  A range
  # of 1200 puts the floor at ceiling(1e8 / 1200^2) = 70 cells, which the
  # layer supports (floor(1200 / 9) = 133) and the subsample alone did not
  # (floor(300 / 9) = 33): the profile used to call these data unsupported,
  # run its ladder from 2 to 33, and hand a count of at most 33 to the
  # tessellation of all 1200 points.
  set.seed(3)
  N <- 1200
  d <- data.frame(x = runif(N, 0, 10000), y = runif(N, 0, 10000), z = rnorm(N))
  pts <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
  sac <- rr_sac(pts, range = 1200, nugget = 0.5, psill = 1)
  lines <- capture_spatialkit_log(
    prof <- resolution_profile(pts, "z", sac = sac, n_levels = 4, nstart = 3,
                               sample_n = 300))
  b <- attr(prof, "bounds")
  expect_identical(b$n, 1200L)
  expect_identical(b$n_sample, 300L)
  expect_identical(b$ceiling, 133L)
  expect_identical(b$ceiling_from, "min_cell_n")
  expect_true(b$supported)
  expect_false(log_has(lines, "cannot support"))
  expect_identical(min(prof$levels), b$floor)
  expect_identical(max(prof$levels), 133L)
  expect_output(print(prof), "on 1200 points")
  # The bounds do not move with sample_n.
  b2 <- attr(resolution_profile(pts, "z", sac = sac, n_levels = 4, nstart = 3,
                                sample_n = 600), "bounds")
  expect_identical(b2[c("floor", "ceiling", "supported", "n")],
                   b[c("floor", "ceiling", "supported", "n")])
  # Reliability is for cells holding the layer's points, not the subsample's.
  vm  <- attr(sac, "variogram_model")
  cf  <- spatialkit:::.vgm_correlation_fn(vm)
  bb  <- sf::st_bbox(pts)
  rbV <- spatialkit:::.rbar_rect(cf, bb[["xmax"]] - bb[["xmin"]], bb[["ymax"]] - bb[["ymin"]])
  expect_equal(prof$reliability,
               vapply(prof$levels, spatialkit:::.reliability_at, numeric(1),
                      area = b$area, n_total = 1200, nugget = 0.5, psill = 1,
                      cor_fn = cf, rbar_V = rbV))
  # Cp: the subsample's RSS with its own optimism added back, plus the
  # variance of cell means built from all 1200 points.  Every subsample cell
  # holds a scored row here, so L_m = L.
  expect_equal(prof$cp, prof$rss / 300 + 0.5 * (prof$levels / 300 + prof$levels / 1200))
})
