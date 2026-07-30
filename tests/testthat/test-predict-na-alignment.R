# ---------------------------------------------------------------------------
# Tests for NA-safe prediction alignment (.clean_row_mask / .expand_predictions)
# ---------------------------------------------------------------------------

test_that(".clean_row_mask identifies NA and non-finite rows", {
  pts <- sf::st_sf(
    y = c(1, NA, 3, 4, Inf),
    x1 = c(10, 20, 30, NA, 50),
    geometry = sf::st_sfc(
      sf::st_point(c(0, 0)), sf::st_point(c(1, 0)),
      sf::st_point(c(0, 1)), sf::st_point(c(1, 1)),
      sf::st_point(c(0.5, 0.5)),
      crs = 32632
    )
  )

  mask <- spatialkit:::.clean_row_mask(pts, "y", "x1", require_response = TRUE)
  # Row 1: ok, Row 2: NA in y, Row 3: ok, Row 4: NA in x1, Row 5: Inf in y
  expect_equal(unname(mask), c(TRUE, FALSE, TRUE, FALSE, FALSE))
})


test_that(".expand_predictions fills dropped rows with NA", {
  preds <- c(10, 30)
  clean_idx <- c(TRUE, FALSE, TRUE, FALSE, FALSE)
  result <- spatialkit:::.expand_predictions(preds, clean_idx, 5L)

  expect_length(result, 5)
  expect_equal(result[1], 10)
  expect_equal(result[3], 30)
  expect_true(is.na(result[2]))
  expect_true(is.na(result[4]))
  expect_true(is.na(result[5]))
})


test_that(".expand_predictions is a no-op when all rows are clean", {
  preds <- c(1, 2, 3)
  clean_idx <- c(TRUE, TRUE, TRUE)
  result <- spatialkit:::.expand_predictions(preds, clean_idx, 3L)
  expect_identical(result, preds)
})


test_that(".clean_row_mask works when response column is absent (out-of-sample)", {
  pts <- sf::st_sf(
    x1 = c(10, NA, 30),
    geometry = sf::st_sfc(
      sf::st_point(c(0, 0)), sf::st_point(c(1, 0)),
      sf::st_point(c(0, 1)),
      crs = 32632
    )
  )

  # Response "y" is not present — with require_response = FALSE,

  # .clean_row_mask should only check predictors.
  mask <- spatialkit:::.clean_row_mask(pts, "y", "x1", require_response = FALSE)
  expect_equal(unname(mask), c(TRUE, FALSE, TRUE))
})


test_that(".clean_row_mask ignores response NAs when require_response = FALSE (regression)", {
  # This is the critical scenario: the user reuses their full dataset (which

  # contains the response column with NAs) for prediction.  Before the fix,
  # rows with NA in the response were silently dropped even though the response
  # is not needed for prediction.
  pts <- sf::st_sf(
    y  = c(NA, NA, 3),
    x1 = c(10, 20, 30),
    geometry = sf::st_sfc(
      sf::st_point(c(0, 0)), sf::st_point(c(1, 0)),
      sf::st_point(c(0, 1)),
      crs = 32632
    )
  )

  # With require_response = FALSE, all three rows should survive because
  # predictor x1 is complete — the NAs in y must be ignored.
  mask <- spatialkit:::.clean_row_mask(pts, "y", "x1", require_response = FALSE)
  expect_equal(unname(mask), c(TRUE, TRUE, TRUE))

  # With require_response = TRUE (training mode), the NAs in y should drop rows.
  mask_train <- spatialkit:::.clean_row_mask(pts, "y", "x1", require_response = TRUE)
  expect_equal(unname(mask_train), c(FALSE, FALSE, TRUE))
})


test_that("prep_model_data succeeds without response when require_response = FALSE", {
  pts <- sf::st_sf(
    x1 = c(10, 20, 30),
    geometry = sf::st_sfc(
      sf::st_point(c(0, 0)), sf::st_point(c(1, 0)),
      sf::st_point(c(0, 1)),
      crs = 32632
    )
  )

  # Should NOT error: response column "y" is absent but not required.
  result <- prep_model_data(pts, "y", "x1", require_response = FALSE)
  expect_equal(nrow(result), 3L)

  # Should still error if a predictor is missing.
  expect_error(
    prep_model_data(pts, "y", "missing_col", require_response = FALSE),
    "missing required column"
  )
})


test_that("prep_model_data ignores response NAs when require_response = FALSE (regression)", {
  # Reused full dataset for prediction — response column present but has NAs.
  pts <- sf::st_sf(
    y  = c(NA, NA, 3),
    x1 = c(10, 20, 30),
    geometry = sf::st_sfc(
      sf::st_point(c(0, 0)), sf::st_point(c(1, 0)),
      sf::st_point(c(0, 1)),
      crs = 32632
    )
  )

  # All 3 rows should survive: NAs in y must be ignored in prediction mode.
  result <- prep_model_data(pts, "y", "x1", require_response = FALSE)
  expect_equal(nrow(result), 3L)

  # In training mode, only row 3 (no NA in y) should survive.
  result_train <- prep_model_data(pts, "y", "x1", require_response = TRUE)
  expect_equal(nrow(result_train), 1L)
})


test_that("prep_model_data still requires response when require_response = TRUE (default)", {
  pts <- sf::st_sf(
    x1 = c(10, 20, 30),
    geometry = sf::st_sfc(
      sf::st_point(c(0, 0)), sf::st_point(c(1, 0)),
      sf::st_point(c(0, 1)),
      crs = 32632
    )
  )

  expect_error(
    prep_model_data(pts, "y", "x1"),
    "missing required column"
  )
})


test_that("predict.bayesian_fit() errors early when predictor columns are missing from newdata", {
  # Build a minimal bayesian_fit-like object with the fields predict() checks.
  fake_fit <- structure(
    list(
      engine         = structure(list(), class = "brmsfit"),
      data_sf        = sf::st_sf(
        y  = c(1, 2, 3),
        x1 = c(4, 5, 6),
        x2 = c(7, 8, 9),
        geometry = sf::st_sfc(
          sf::st_point(c(0, 0)),
          sf::st_point(c(1, 0)),
          sf::st_point(c(0, 1)),
          crs = 32632
        )
      ),
      response_var   = "y",
      predictor_vars = c("x1", "x2")
    ),
    class = c("bayesian_fit", "spatial_fit")
  )

  # newdata missing "x2" — should error immediately with a clear message.
  bad_nd <- sf::st_sf(
    x1 = c(10, 20),
    geometry = sf::st_sfc(
      sf::st_point(c(0, 0)),
      sf::st_point(c(1, 0)),
      crs = 32632
    )
  )

  expect_error(
    predict(fake_fit, newdata = bad_nd),
    "missing required predictor column.*x2"
  )
})
