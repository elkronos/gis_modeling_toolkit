# ===========================================================================
# Spatial weight matrices for residual_morans_i().
#
# Two failure modes, both silent: a sparse matrix being dropped in favour of
# the kNN fallback (giving a different statistic than the one asked for),
# and non-row-standardised weights being treated as though they were.
# See test-morans-variance.R for the variance formula itself.
# ===========================================================================

test_that("residual_morans_i accepts sparse Matrix weights (no kNN fallback)", {
  skip_if_not_installed("sf")
  skip_if_not_installed("Matrix")

  set.seed(21)
  n <- 30
  pts <- sf::st_as_sf(
    data.frame(x = runif(n), y = runif(n), resp = rnorm(n)),
    coords = c("x", "y"), crs = 32631
  )
  resid_vec <- rnorm(n)
  fake_fit <- structure(
    list(data_sf = pts, residuals = resid_vec, engine = list()),
    class = "spatial_fit"
  )
  registerS3method("residuals", "spatial_fit",
                   function(object, ...) object$residuals)

  # Row-standardised weight matrix (3 neighbours each, weight 1/3)
  W <- matrix(0, n, n)
  for (i in seq_len(n)) {
    nbrs <- sample(setdiff(seq_len(n), i), 3)
    W[i, nbrs] <- 1 / 3
  }
  W_sparse <- Matrix::Matrix(W, sparse = TRUE)

  res_dense  <- residual_morans_i(fake_fit, weights = W)
  res_sparse <- residual_morans_i(fake_fit, weights = W_sparse)

  # Sparse weights must be USED (identical results), not silently replaced
  # by the kNN fallback (which would give different I / sd / z).
  expect_equal(res_sparse$observed, res_dense$observed, tolerance = 1e-12)
  expect_equal(res_sparse$sd,       res_dense$sd,       tolerance = 1e-12)
  expect_equal(res_sparse$z,        res_dense$z,        tolerance = 1e-12)
  expect_equal(res_sparse$p_value,  res_dense$p_value,  tolerance = 1e-12)
})


test_that("residual_morans_i warns on non-row-standardised user weights", {
  skip_if_not_installed("sf")
  skip_if_not_installed("FNN")

  # Build a minimal spatial_fit stub with residuals and coordinates
  n <- 20
  set.seed(42)
  coords_mat <- cbind(runif(n), runif(n))
  pts <- sf::st_as_sf(
    data.frame(x = coords_mat[, 1], y = coords_mat[, 2], resp = rnorm(n)),
    coords = c("x", "y"), crs = 32631
  )

  fake_fit <- structure(
    list(
      data_sf  = pts,
      residuals = rnorm(n),
      engine   = list()
    ),
    class = c("spatial_fit")
  )
  residuals.spatial_fit <- function(object, ...) object$residuals
  registerS3method("residuals", "spatial_fit", residuals.spatial_fit)

  # A row-standardised matrix should NOT trigger the row-standardisation warning
  W_good <- matrix(0, n, n)
  for (i in seq_len(n)) {
    nbrs <- sample(setdiff(seq_len(n), i), 4)
    W_good[i, nbrs] <- 1 / 4
  }
  result <- residual_morans_i(fake_fit, weights = W_good)
  expect_true(is.list(result))
  expect_true(is.finite(result$observed))

  # A non-row-standardised (binary) matrix SHOULD trigger a logger warning.
  # The spatialkit logger appends to a temp file (configured in .onLoad),
  # so we read new log entries from that file after the call.
  W_bad <- matrix(0, n, n)
  for (i in seq_len(n)) {
    nbrs <- sample(setdiff(seq_len(n), i), 4)
    W_bad[i, nbrs] <- 1
  }

  log_file <- file.path(tempdir(), "spatialkit_model_log.log")
  # Record the file size before so we can isolate new entries
  before_size <- if (file.exists(log_file)) file.info(log_file)$size else 0L

  residual_morans_i(fake_fit, weights = W_bad)

  # Read any new content appended to the log file after the call
  if (file.exists(log_file) && file.info(log_file)$size > before_size) {
    con <- file(log_file, open = "rb")
    on.exit(close(con), add = TRUE)
    if (before_size > 0) seek(con, where = before_size)
    new_log <- readLines(con, warn = FALSE)
  } else {
    new_log <- character(0)
  }

  expect_true(any(grepl("not row-standardised", new_log)))
})
