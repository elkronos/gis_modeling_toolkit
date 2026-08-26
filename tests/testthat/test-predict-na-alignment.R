# ---------------------------------------------------------------------------
# Tests for NA-safe prediction alignment
# (prep_model_data() row cleaning + .expand_predictions)
#
# Every predict method now attaches an `..orig_row_id..` sentinel column,
# hands newdata to prep_model_data(require_response = FALSE), and derives
# `clean_idx` from the sentinel values that survived.  That is the single
# source of truth for "which rows were dropped"; the tests below assert
# against it rather than against a private mask helper.
# ---------------------------------------------------------------------------

test_that("prep_model_data drops NA and non-finite rows", {
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
  pts$row <- seq_len(5L)

  # Row 1: ok, Row 2: NA in y, Row 3: ok, Row 4: NA in x1, Row 5: Inf in y
  clean <- prep_model_data(pts, "y", "x1", require_response = TRUE)
  expect_equal(nrow(clean), 2L)
  expect_equal(clean$row, c(1L, 3L))

  # Prediction mode ignores the response entirely, so only x1 removes a row.
  clean_pred <- prep_model_data(pts, "y", "x1", require_response = FALSE)
  expect_equal(nrow(clean_pred), 4L)
  expect_equal(clean_pred$row, c(1L, 2L, 3L, 5L))
})


test_that("the ..orig_row_id.. sentinel recovers exactly which rows survived", {
  # This is the mechanism every predict method uses: stamp the sentinel,
  # clean, then read the surviving ids back.  Asserting it here keeps the
  # regression intent of the deleted .clean_row_mask tests alive against the
  # single implementation, with no optional backend required.
  pts <- sf::st_sf(
    y  = c(1, NA, 3, 4, Inf),
    x1 = c(10, 20, 30, NA, 50),
    geometry = sf::st_sfc(
      sf::st_point(c(0, 0)), sf::st_point(c(1, 0)),
      sf::st_point(c(0, 1)), sf::st_point(c(1, 1)),
      sf::st_point(c(0.5, 0.5)),
      crs = 32632
    )
  )
  n_orig <- nrow(pts)
  pts$..orig_row_id.. <- seq_len(n_orig)

  cleaned   <- prep_model_data(pts, "y", "x1", require_response = FALSE)
  clean_idx <- seq_len(n_orig) %in% cleaned$..orig_row_id..
  expect_equal(clean_idx, c(TRUE, TRUE, TRUE, FALSE, TRUE))

  # ... and predictions on the clean subset realign to their original slots.
  preds_clean <- c(10, 20, 30, 50)
  out <- spatialkit:::.expand_predictions(preds_clean, clean_idx, n_orig)
  expect_length(out, n_orig)
  expect_equal(out, c(10, 20, 30, NA, 50))

  # The sentinel is inert: it never becomes a modelling column and never
  # changes which rows are kept.
  expect_equal(nrow(cleaned),
               nrow(prep_model_data(pts[, c("y", "x1")], "y", "x1",
                                    require_response = FALSE)))
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


test_that(".expand_predictions errors rather than recycling on a length mismatch", {
  # 4 clean slots and a length-2 `preds` used to fill 1 2 . 1 . 2 and hand
  # that back as a prediction vector.
  expect_error(
    spatialkit:::.expand_predictions(c(1, 2), c(TRUE, TRUE, FALSE, TRUE, TRUE), 5L),
    "2 prediction value\\(s\\) for 4 row\\(s\\)"
  )
  # The equal-length fast path must be guarded too: n_clean == n_orig but
  # preds is short.
  expect_error(
    spatialkit:::.expand_predictions(c(1, 2), rep(TRUE, 4L), 4L),
    "cannot be aligned"
  )
  # Too many predictions is equally an error.
  expect_error(
    spatialkit:::.expand_predictions(c(1, 2, 3), c(TRUE, FALSE, TRUE), 3L),
    "3 prediction value\\(s\\) for 2 row\\(s\\)"
  )
})


test_that("predict() realigns to the original row positions end to end", {
  # The regression the sentinel path exists for: a newdata row carrying an NA
  # predictor must come back as NA *in its own position*, not shift the rest
  # of the vector up.
  skip_if_not_installed("ranger")
  set.seed(11)
  n <- 60
  train <- sf::st_as_sf(
    data.frame(x = runif(n, 0, 100), y = runif(n, 0, 100), x1 = rnorm(n)),
    coords = c("x", "y"), crs = 32632
  )
  train$z <- 3 * train$x1 + rnorm(n, 0, 0.1)
  fit <- fit_rf_model(train, "z", "x1", num_trees = 60, seed = 1)

  nd <- sf::st_as_sf(
    data.frame(x = c(10, 20, 30, 40, 50), y = c(10, 20, 30, 40, 50),
               x1 = c(-1, NA, 0, Inf, 1)),
    coords = c("x", "y"), crs = 32632
  )
  p <- suppressWarnings(predict(fit, newdata = nd))

  expect_length(p, 5L)
  expect_equal(which(is.na(p)), c(2L, 4L))

  # The surviving predictions equal what the same rows get on their own.
  p_ok <- suppressWarnings(predict(fit, newdata = nd[c(1L, 3L, 5L), ]))
  expect_equal(unname(p[c(1L, 3L, 5L)]), unname(p_ok))
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
