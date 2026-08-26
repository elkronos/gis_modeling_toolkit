# ===========================================================================
# .compute_reg_metrics(): the shared metric calculator behind summary(),
# model_metrics() and every cross-validation path. The y_train_mean tests
# matter because R-squared against the WRONG baseline is the classic way to
# report a flattering number: out-of-sample R-squared must be measured
# against the training mean, not the test fold's own mean.
# ===========================================================================

test_that(".compute_reg_metrics matches hand-computed values", {
  y    <- c(1, 2, 3, 4, 5)
  yhat <- c(1.1, 2.2, 2.9, 4.1, 4.8)
  met  <- spatialkit:::.compute_reg_metrics(y, yhat)

  expect_s3_class(met, "data.frame")
  expect_named(met, c("n", "RMSE", "MAE", "MAPE", "SMAPE", "R2", "Adj_R2"))
  expect_identical(met$n, 5L)

  # Errors are -0.1, -0.2, +0.1, -0.1, +0.2, so:
  #   RSS = 0.01 + 0.04 + 0.01 + 0.01 + 0.04 = 0.11
  #   TSS = sum((y - 3)^2) = 4 + 1 + 0 + 1 + 4 = 10
  expect_equal(met$RMSE, sqrt(0.11 / 5))                       # 0.1483240
  expect_equal(met$MAE, (0.1 + 0.2 + 0.1 + 0.1 + 0.2) / 5)     # 0.14
  expect_equal(met$R2, 1 - 0.11 / 10)                          # 0.989

  # MAPE = mean(|e| / |y|) * 100
  expect_equal(met$MAPE,
               mean(c(0.1 / 1, 0.2 / 2, 0.1 / 3, 0.1 / 4, 0.2 / 5)) * 100)
  # SMAPE = mean(2 |e| / (|y| + |yhat|)) * 100
  expect_equal(met$SMAPE,
               mean(c(0.2 / 2.1, 0.4 / 4.2, 0.2 / 5.9,
                      0.2 / 8.1, 0.4 / 9.8)) * 100)

  # Adj R2 is suppressed unless `p` is given.
  expect_true(is.na(met$Adj_R2))
})


test_that(".compute_reg_metrics computes Adj R2 only when p leaves df to spare", {
  y    <- c(1, 2, 3, 4, 5)
  yhat <- c(1.1, 2.2, 2.9, 4.1, 4.8)
  r2   <- 1 - 0.11 / 10

  # n = 5, p = 1  ->  1 - (1 - R2) * (n - 1) / (n - p - 1)
  m1 <- spatialkit:::.compute_reg_metrics(y, yhat, p = 1)
  expect_equal(m1$Adj_R2, 1 - (1 - r2) * 4 / 3)
  expect_lt(m1$Adj_R2, m1$R2)          # the penalty always costs something

  m2 <- spatialkit:::.compute_reg_metrics(y, yhat, p = 2)
  expect_equal(m2$Adj_R2, 1 - (1 - r2) * 4 / 2)
  expect_lt(m2$Adj_R2, m1$Adj_R2)      # ... and more with more predictors

  # n == p + 1 would divide by zero; n < p + 1 would flip the sign.  Both
  # return NA rather than a number.
  expect_true(is.na(spatialkit:::.compute_reg_metrics(y, yhat, p = 4)$Adj_R2))
  expect_true(is.na(spatialkit:::.compute_reg_metrics(y, yhat, p = 9)$Adj_R2))
})


test_that(".compute_reg_metrics returns an all-NA row for an empty comparison", {
  # .cv_overall_metrics() relies on this shape when a fold produced nothing:
  # n = 0 and every metric NA, with the columns still present so rbind() of
  # fold results does not ragged out.
  empty <- spatialkit:::.compute_reg_metrics(numeric(0), numeric(0))
  expect_identical(empty$n, 0L)
  expect_named(empty, c("n", "RMSE", "MAE", "MAPE", "SMAPE", "R2", "Adj_R2"))
  expect_true(all(is.na(unlist(empty[, -1L]))))

  # Same when every row is non-finite rather than absent.
  allna <- spatialkit:::.compute_reg_metrics(c(NA, NaN, Inf), c(1, 2, 3))
  expect_identical(allna$n, 0L)
  expect_true(all(is.na(unlist(allna[, -1L]))))

  # A non-finite PREDICTION drops its row too, and the rest still score.
  partial <- spatialkit:::.compute_reg_metrics(c(1, 2, 3), c(1, NA, 3))
  expect_identical(partial$n, 2L)
  expect_equal(partial$RMSE, 0)
})


test_that(".compute_reg_metrics degrades the ratio metrics rather than dividing by zero", {
  # MAPE is undefined at y = 0 and SMAPE at y = yhat = 0.  Both must come back
  # NA instead of Inf/NaN, which would poison a mean over folds.
  zeros <- spatialkit:::.compute_reg_metrics(c(0, 0, 0), c(0.1, -0.1, 0.2))
  expect_true(is.na(zeros$MAPE))
  expect_false(is.na(zeros$SMAPE))     # denominator is |y| + |yhat| > 0

  both_zero <- spatialkit:::.compute_reg_metrics(c(0, 0), c(0, 0))
  expect_true(is.na(both_zero$MAPE))
  expect_true(is.na(both_zero$SMAPE))
  expect_equal(both_zero$RMSE, 0)

  # A constant response has zero TSS, so R2 is undefined rather than -Inf.
  constant <- spatialkit:::.compute_reg_metrics(rep(5, 4), c(5, 5.1, 4.9, 5))
  expect_true(is.na(constant$R2))
  expect_true(is.finite(constant$RMSE))
})


test_that(".compute_reg_metrics rejects a y_train_mean of the wrong length", {
  # Recycling a length-2 baseline against 5 observations silently produced a
  # wrong R2; it is now an error naming both lengths.
  y    <- c(1, 2, 3, 4, 5)
  yhat <- c(1.1, 2.2, 2.9, 4.1, 4.8)
  expect_error(
    spatialkit:::.compute_reg_metrics(y, yhat, y_train_mean = c(1, 2)),
    "must be length 1 .* or length 5 .*; got 2"
  )
  # Length 1 (scalar baseline) and length n (per observation) are both fine --
  # the check must not have become "length n only".
  expect_s3_class(spatialkit:::.compute_reg_metrics(y, yhat, y_train_mean = 3),
                  "data.frame")
  expect_s3_class(
    spatialkit:::.compute_reg_metrics(y, yhat, y_train_mean = rep(3, 5)),
    "data.frame")
})


test_that(".compute_reg_metrics uses scalar y_train_mean as baseline", {
  y    <- c(1, 2, 3, 4, 5)
  yhat <- c(1.1, 2.2, 2.9, 4.1, 4.8)
  train_mean <- 2.5
  met <- spatialkit:::.compute_reg_metrics(y, yhat, y_train_mean = train_mean)
  # Manually compute expected R² with the scalar baseline
  rss <- sum((y - yhat)^2)
  tss <- sum((y - train_mean)^2)
  expect_equal(met$R2, 1 - rss / tss)
})


test_that(".compute_reg_metrics uses per-observation y_train_mean vector as baseline", {
  y    <- c(1, 2, 3, 4, 5)
  yhat <- c(1.1, 2.2, 2.9, 4.1, 4.8)
  ytm_vec <- c(2.0, 2.5, 3.0, 3.5, 4.0)
  met <- spatialkit:::.compute_reg_metrics(y, yhat, y_train_mean = ytm_vec)
  rss <- sum((y - yhat)^2)
  tss <- sum((y - ytm_vec)^2)
  expect_equal(met$R2, 1 - rss / tss)
})


test_that(".compute_reg_metrics filters NA with per-observation y_train_mean", {
  y    <- c(1, 2, NA, 4, 5)
  yhat <- c(1.1, 2.2, 2.9, 4.1, 4.8)
  ytm_vec <- c(2.0, 2.5, 3.0, 3.5, 4.0)
  met <- spatialkit:::.compute_reg_metrics(y, yhat, y_train_mean = ytm_vec)
  ok <- c(TRUE, TRUE, FALSE, TRUE, TRUE)
  rss <- sum((y[ok] - yhat[ok])^2)
  tss <- sum((y[ok] - ytm_vec[ok])^2)
  expect_equal(met$n, 4L)
  expect_equal(met$R2, 1 - rss / tss)
})
