# ===========================================================================
# .compute_reg_metrics(): the shared metric calculator behind summary(),
# model_metrics() and every cross-validation path. The y_train_mean tests
# matter because R-squared against the WRONG baseline is the classic way to
# report a flattering number: out-of-sample R-squared must be measured
# against the training mean, not the test fold's own mean.
# ===========================================================================

test_that(".compute_reg_metrics returns correct structure", {
  y    <- c(1, 2, 3, 4, 5)
  yhat <- c(1.1, 2.2, 2.9, 4.1, 4.8)
  met  <- spatialkit:::.compute_reg_metrics(y, yhat)
  expect_s3_class(met, "data.frame")
  expect_true(all(c("n", "RMSE", "MAE", "R2") %in% names(met)))
  expect_equal(met$n, 5L)
  expect_true(met$R2 > 0.9)
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
