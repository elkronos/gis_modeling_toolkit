# ---------------------------------------------------------------------------
# Tests for estimate_sac_range() OLS residual alignment guard
# ---------------------------------------------------------------------------

test_that("estimate_sac_range falls back to raw response when residual length mismatches", {
  skip_if_not_installed("gstat")

  set.seed(42)
  n <- 60
  pts <- sf::st_sf(
    resp = rnorm(n),
    pred = rnorm(n),
    geometry = sf::st_sfc(
      lapply(seq_len(n), function(i) sf::st_point(c(i * 100, 0))),
      crs = 32632
    )
  )

  # Normal case: residual vector aligns, no warning expected about mismatch
  expect_no_warning(
    estimate_sac_range(pts, response_var = "resp", predictor_vars = "pred"),
    message = "OLS residual length"
  )

  # Inject NAs in predictor — na.exclude should keep alignment, no mismatch warning
  pts$pred[c(5, 10)] <- NA
  expect_no_warning(
    estimate_sac_range(pts, response_var = "resp", predictor_vars = "pred"),
    message = "OLS residual length"
  )
})
