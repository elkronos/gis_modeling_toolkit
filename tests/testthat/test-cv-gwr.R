# cv_gwr() had one call site in the whole suite -- an expect_error() on
# argument validation -- and its functional path ran only as an unasserted
# internal step of one compare_models_cv() test.  Nothing checked that the
# number it returns describes the predictions it returns.

test_that("cv_gwr() pools the fold predictions it returns", {
  skip_if_not_installed("GWmodel")
  skip_if_not_installed("sp")

  set.seed(21); n <- 90
  d <- sf::st_as_sf(
    data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
               elev = rnorm(n)),
    coords = c("x", "y"), crs = 32632)
  d$price <- 10 + 2 * d$elev + rnorm(n, 0, 0.5)

  folds <- make_folds(d, k = 3, method = "block_kfold", seed = 7)
  cv <- suppressWarnings(suppressMessages(
    cv_gwr(d, "price", "elev", folds = folds, bandwidth = 30)))

  # Every fold ran, and every row was scored exactly once.
  expect_equal(cv$n_folds_attempted, 3L)
  expect_equal(cv$n_folds_succeeded, 3L)
  expect_equal(nrow(cv$fold_metrics), 3L)
  expect_true(all(cv$fold_status$status == "ok"))
  expect_equal(nrow(cv$predictions), n)
  expect_equal(sort(cv$predictions$..row_id), seq_len(n))
  expect_equal(cv$overall$n_pred, n)
  expect_equal(cv$n_dropped, 0L)
  expect_length(cv$orphan_rows, 0L)
  expect_equal(cv$n_unknown_ids, 0L)

  # The pooled metrics must be the metrics OF those predictions, not a mean of
  # fold means: recompute them by hand from $predictions.
  p <- cv$predictions
  expect_equal(cv$overall$RMSE, sqrt(mean((p$y - p$yhat)^2)))
  expect_equal(cv$overall$MAE,  mean(abs(p$y - p$yhat)))
  expect_equal(cv$overall$R2,
               1 - sum((p$y - p$yhat)^2) / sum((p$y - p$y_train_mean)^2))

  # Documented: Adj_R2 is NA for pooled CV, because the pooled predictions come
  # from k separately fitted models with no single parameter count.
  expect_true(is.na(cv$overall$Adj_R2))

  # The folds it reports are the folds it was given.
  key <- function(f) lapply(f, function(z) sort(z$test))
  expect_equal(key(cv$folds), key(folds$folds))

  # A recognisable signal must be recovered: the response is 2 * elev plus
  # small noise, so a cross-validated R2 near zero would mean the backend ran
  # but learned nothing.
  expect_gt(cv$overall$R2, 0.8)
})
