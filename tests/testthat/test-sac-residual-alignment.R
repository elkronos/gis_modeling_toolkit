# ---------------------------------------------------------------------------
# Tests for estimate_sac_range() OLS residual alignment guard.
#
# When predictor_vars are supplied the variogram is fitted to OLS residuals
# rather than the raw response.  If the residual vector ever came back a
# different length from the data, pairing it with the coordinates would model
# the wrong value at every location -- so the function checks the length and
# falls back to the raw response.
#
# Asserting only that the guard does NOT fire leaves it untested: the whole
# block could be deleted and this file would stay green.  The mocked test
# below drives the guard itself.
# ---------------------------------------------------------------------------

# Captured before any mock is installed, so the mock can delegate to it.
.real_lm <- stats::lm

.sac_pts <- function(n = 60, seed = 42) {
  set.seed(seed)
  sf::st_sf(
    resp = rnorm(n),
    pred = rnorm(n),
    geometry = sf::st_sfc(
      lapply(seq_len(n), function(i) sf::st_point(c(i * 100, 0))),
      crs = 32632
    )
  )
}

test_that("estimate_sac_range uses OLS residuals when they align", {
  skip_if_not_installed("gstat")
  pts <- .sac_pts()

  lines <- capture_spatialkit_log(
    with_pred <- estimate_sac_range(pts, response_var = "resp",
                                    predictor_vars = "pred")
  )
  expect_false(log_has(lines, "OLS residual length"))

  # Residuals of resp ~ pred are not resp, so the two calls must disagree --
  # otherwise "no warning" would be satisfiable by ignoring predictor_vars.
  raw <- estimate_sac_range(pts, response_var = "resp")
  expect_false(isTRUE(all.equal(as.numeric(with_pred), as.numeric(raw))))
})

test_that("na.exclude keeps residuals aligned when a predictor has NAs", {
  skip_if_not_installed("gstat")
  pts <- .sac_pts()
  pts$pred[c(5, 10)] <- NA

  lines <- capture_spatialkit_log(
    est <- estimate_sac_range(pts, response_var = "resp",
                              predictor_vars = "pred")
  )
  # na.exclude pads the residual vector back to nrow(pts); na.omit would not,
  # and the guard would fire.
  expect_false(log_has(lines, "OLS residual length"))
  expect_true(is.na(est) || is.finite(as.numeric(est)))
})

test_that("a misaligned residual vector is detected and the raw response used", {
  # The positive case.  A length mismatch cannot be produced through the public
  # interface -- which is the point of the guard being a backstop -- so lm() is
  # mocked to drop one residual.
  skip_if_not_installed("gstat")
  pts <- .sac_pts()

  raw <- estimate_sac_range(pts, response_var = "resp")

  lines <- capture_spatialkit_log({
    local_mocked_bindings(
      lm = function(...) {
        f <- .real_lm(...)
        f$residuals <- f$residuals[-1L]
        f
      },
      .package = "stats"
    )
    fallback <- estimate_sac_range(pts, response_var = "resp",
                                   predictor_vars = "pred")
  })

  expect_true(log_has(lines, "OLS residual length"))
  expect_true(log_has(lines, "using raw response"))
  # It really did fall back: the answer is the raw-response answer, not the
  # residual one (which the first test showed differs).
  expect_equal(as.numeric(fallback), as.numeric(raw))
})
