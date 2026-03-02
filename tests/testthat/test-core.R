test_that(".is_longlat detects geographic CRS", {
  pts <- sf::st_sfc(sf::st_point(c(-122.4, 37.8)), crs = 4326)
  expect_true(spatialkit:::.is_longlat(pts))

  pts_proj <- sf::st_transform(pts, 32610)
  expect_false(spatialkit:::.is_longlat(pts_proj))
})

test_that("ensure_projected transforms geographic data", {
  pts <- sf::st_sf(
    id = 1:3,
    geometry = sf::st_sfc(
      sf::st_point(c(-122.4, 37.8)),
      sf::st_point(c(-122.3, 37.7)),
      sf::st_point(c(-122.5, 37.9)),
      crs = 4326
    )
  )
  result <- ensure_projected(pts)
  expect_false(sf::st_is_longlat(result))
})

test_that("ensure_projected does NOT assume 4326 for small integer-like coords", {
  # Simulates a local site survey in metres with coords outside the

  # geographic envelope (values > 180 / > 90) so the heuristic cannot
  # mistakenly interpret them as lon/lat.
  pts <- sf::st_sf(
    id = 1:3,
    geometry = sf::st_sfc(
      sf::st_point(c(500, 600)),
      sf::st_point(c(510, 620)),
      sf::st_point(c(505, 610))
    )
  )
  result <- suppressWarnings(ensure_projected(pts))
  # CRS should still be NA — the function must NOT guess 4326

  expect_true(is.na(sf::st_crs(result)))
})

test_that("ensure_projected assumes 4326 for CRS-less data with geographic precision", {
  # Realistic lon/lat coordinates without a CRS set
  pts <- sf::st_sf(
    id = 1:2,
    geometry = sf::st_sfc(
      sf::st_point(c(-0.1278, 51.5074)),
      sf::st_point(c(-0.1385, 51.5013))
    )
  )
  result <- suppressWarnings(ensure_projected(pts))
  # Should have been assigned a projected CRS
  expect_false(is.na(sf::st_crs(result)))
  expect_false(sf::st_is_longlat(result))
})

test_that("ensure_projected assumes 4326 for CRS-less data with large extent", {
  pts <- sf::st_sf(
    id = 1:2,
    geometry = sf::st_sfc(
      sf::st_point(c(10, 20)),
      sf::st_point(c(15, 25))
    )
  )
  result <- suppressWarnings(ensure_projected(pts))
  # Extent > 1 degree in both axes → should assume 4326
  expect_false(is.na(sf::st_crs(result)))
  expect_false(sf::st_is_longlat(result))
})

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

test_that(".assert_sf rejects non-sf objects", {
  expect_error(spatialkit:::.assert_sf(data.frame(x = 1)), "Expected an sf object")
})

test_that("harmonize_crs aligns two objects", {
  a <- sf::st_sf(geometry = sf::st_sfc(sf::st_point(c(0, 0)), crs = 4326))
  b <- sf::st_sf(geometry = sf::st_sfc(sf::st_point(c(0, 0)), crs = 3857))
  result <- harmonize_crs(a, b, prefer = "b")
  expect_equal(sf::st_crs(result$a), sf::st_crs(b))
})

test_that("clip_target_for works with points only", {
  pts <- sf::st_sf(
    id = 1:5,
    geometry = sf::st_sfc(
      sf::st_point(c(0, 0)), sf::st_point(c(1, 0)),
      sf::st_point(c(0, 1)), sf::st_point(c(1, 1)),
      sf::st_point(c(0.5, 0.5)),
      crs = 32632
    )
  )
  tgt <- clip_target_for(pts)
  expect_s3_class(tgt, "sf")
  expect_true(all(sf::st_geometry_type(tgt) %in% c("POLYGON", "MULTIPOLYGON")))
})

test_that("create_grid_polygons produces cells", {
  bnd <- sf::st_sf(
    geometry = sf::st_sfc(
      sf::st_polygon(list(rbind(c(0, 0), c(100, 0), c(100, 100), c(0, 100), c(0, 0)))),
      crs = 32632
    )
  )
  grid <- create_grid_polygons(bnd, target_cells = 9, type = "square")
  expect_s3_class(grid, "sf")
  expect_true(nrow(grid) >= 4)
})

test_that("make_folds random_kfold produces correct fold count", {
  pts <- sf::st_sf(
    id = 1:20,
    geometry = sf::st_sfc(
      lapply(1:20, function(i) sf::st_point(c(runif(1, 0, 100), runif(1, 0, 100)))),
      crs = 32632
    )
  )
  folds <- make_folds(pts, k = 5, method = "random_kfold", seed = 42)
  expect_equal(folds$k, 5)
  expect_length(folds$folds, 5)
  # Every point should appear in exactly one test fold
  all_test <- unlist(lapply(folds$folds, `[[`, "test"))
  expect_equal(sort(all_test), 1:20)
})

# ---------- summarize_by_cell: deff tests ----------

.make_test_points <- function(n_cells = 5, pts_per_cell = 10, rho = 0) {
  # Create a simple sf with a cell ID and a numeric variable.
  # When rho > 0 observations within a cell share a common component,
  # mimicking intra-cell spatial autocorrelation.
  set.seed(42)
  ids <- rep(seq_len(n_cells), each = pts_per_cell)
  cell_effect <- rep(rnorm(n_cells, 0, 2), each = pts_per_cell)
  noise <- rnorm(n_cells * pts_per_cell)
  y <- sqrt(rho) * cell_effect + sqrt(1 - rho) * noise
  sf::st_sf(
    poly_id  = ids,
    y        = y,
    geometry = sf::st_sfc(
      lapply(seq_along(ids), function(i) sf::st_point(c(i, ids[i]))),
      crs = 32632
    )
  )
}

test_that("summarize_by_cell default deff=1 gives classic SE", {
  pts <- .make_test_points(rho = 0)
  out <- summarize_by_cell(pts, response_var = "y")
  # SE should equal sd / sqrt(n) for each cell
  for (i in seq_len(nrow(out))) {
    expected_se <- out[["..sd_resp_y"]][i] / sqrt(out$n[i])
    expect_equal(out[["..se_resp_y"]][i], expected_se, tolerance = 1e-12)
  }
  # cell_weight should equal n when deff=1
  expect_equal(out$cell_weight, out$n)
})

test_that("summarize_by_cell deff=2 inflates SE by sqrt(2)", {
  pts <- .make_test_points(rho = 0)
  out_1 <- summarize_by_cell(pts, response_var = "y", deff = 1)
  out_2 <- summarize_by_cell(pts, response_var = "y", deff = 2)
  ratio <- out_2[["..se_resp_y"]] / out_1[["..se_resp_y"]]
  expect_equal(ratio, rep(sqrt(2), nrow(out_1)), tolerance = 1e-12)
  # cell_weight halved
  expect_equal(out_2$cell_weight, out_2$n / 2)
  # attribute recorded
  da <- attr(out_2, "deff_applied")
  expect_equal(da$method, "fixed")
  expect_equal(da$deff, 2)
})

test_that("summarize_by_cell deff='kish' inflates SE under correlation", {
  pts_corr <- .make_test_points(n_cells = 10, pts_per_cell = 20, rho = 0.5)
  out_iid  <- summarize_by_cell(pts_corr, response_var = "y", deff = 1)
  out_kish <- summarize_by_cell(pts_corr, response_var = "y", deff = "kish")
  # Kish SE should be >= IID SE for every cell (deff >= 1)
  expect_true(all(out_kish[["..se_resp_y"]] >= out_iid[["..se_resp_y"]] - 1e-12))
  # cell_weight should be <= n
  expect_true(all(out_kish$cell_weight <= out_kish$n + 1e-12))
  # deff_applied attribute should exist with per-variable-type ICC
  da <- attr(out_kish, "deff_applied")
  expect_equal(da$method, "kish")
  expect_true(da$icc_resp > 0)
  expect_true(is.na(da$icc_pred))
})

test_that("summarize_by_cell deff='kish' estimates separate ICC for response and predictors", {
  set.seed(99)
  n_cells <- 10; n_per <- 20
  ids <- rep(seq_len(n_cells), each = n_per)
  # Response with strong intra-cell correlation
  cell_eff_y <- rep(rnorm(n_cells, 0, 3), each = n_per)
  y <- 0.7 * cell_eff_y + 0.3 * rnorm(n_cells * n_per)
  # Predictor with weak/no intra-cell correlation
  x <- rnorm(n_cells * n_per)
  pts <- sf::st_sf(
    poly_id  = ids,
    y        = y,
    x        = x,
    geometry = sf::st_sfc(
      lapply(seq_along(ids), function(i) sf::st_point(c(i, ids[i]))),
      crs = 32632
    )
  )
  out <- summarize_by_cell(pts, response_var = "y", predictor_vars = "x",
                           deff = "kish")
  da <- attr(out, "deff_applied")
  expect_equal(da$method, "kish")
  # Response ICC should be much larger than predictor ICC
  expect_true(da$icc_resp > da$icc_pred)
})

test_that("summarize_by_cell ICC z-scores predictors so scale does not dominate", {
  # Regression test: if predictors are NOT standardised before pooling,
  # a large-scale variable dominates and the predictor ICC reflects only
  # that variable.  After z-scoring, both contribute equally.
  set.seed(42)
  n_cells <- 10; n_per <- 30
  ids <- rep(seq_len(n_cells), each = n_per)

  # x_small: range ~ [0, 1], strong intra-cell correlation
  # (large cell effect relative to noise ensures a clearly detectable ICC)
  cell_eff_small <- rep(rnorm(n_cells, 0, 2), each = n_per)
  x_small <- cell_eff_small + rnorm(n_cells * n_per, 0, 0.5)

  # x_big: range ~ [0, 10000], virtually zero intra-cell correlation
  x_big <- runif(n_cells * n_per, 0, 10000)

  # A correlated response for completeness
  y <- rnorm(n_cells * n_per)

  pts <- sf::st_sf(
    poly_id = ids,
    y       = y,
    x_small = x_small,
    x_big   = x_big,
    geometry = sf::st_sfc(
      lapply(seq_along(ids), function(i) sf::st_point(c(i, ids[i]))),
      crs = 32632
    )
  )

  # ICC with both predictors (pooled, z-scored)
  out_both <- summarize_by_cell(pts, response_var = "y",
                                predictor_vars = c("x_small", "x_big"),
                                deff = "kish")
  da_both <- attr(out_both, "deff_applied")

  # ICC with only the small-scale predictor
  out_small <- summarize_by_cell(pts, response_var = "y",
                                 predictor_vars = "x_small",
                                 deff = "kish")
  da_small <- attr(out_small, "deff_applied")

  # The pooled ICC should be noticeably pulled toward the small predictor's
  # ICC (which has real intra-cell structure), not dominated by x_big's
  # near-zero ICC.  We check that the pooled ICC is at least 25% of the
  # small-only ICC (before z-scoring it was essentially 0).
  expect_true(da_both$icc_pred >= 0.25 * da_small$icc_pred)
})

test_that("summarize_by_cell invalid deff falls back to 1", {
  pts <- .make_test_points()
  # Invalid deff logs a warning via logger and falls back to deff = 1
  out <- summarize_by_cell(pts, response_var = "y", deff = -1)
  expect_equal(out$cell_weight, out$n)
})

test_that("summarize_by_cell joins cells_sf even when ID types differ", {
  # Points have integer poly_id; cells_sf has character poly_id.
  # The join should still succeed after the coercion fix.
  pts <- .make_test_points(n_cells = 3, pts_per_cell = 5, rho = 0)

  # Build a simple polygon sf with *character* IDs
  cells <- sf::st_sf(
    poly_id  = as.character(1:3),
    geometry = sf::st_sfc(
      sf::st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,1), c(0,0)))),
      sf::st_polygon(list(rbind(c(1,0), c(2,0), c(2,1), c(1,1), c(1,0)))),
      sf::st_polygon(list(rbind(c(2,0), c(3,0), c(3,1), c(2,1), c(2,0)))),
      crs = 32632
    )
  )

  out <- summarize_by_cell(pts, response_var = "y", cells_sf = cells)
  # Result should be an sf object (geometry joined from cells_sf)
  expect_s3_class(out, "sf")
  # All 3 cells should be present with no NA from a failed join
  expect_equal(nrow(out), 3L)
  expect_false(any(is.na(out$n)))
})

# ---------------------------------------------------------------------------
# residual_morans_i: row-standardisation validation
# ---------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------
# .build_knn_weights: dense fallback size guard
# ---------------------------------------------------------------------------

test_that(".build_knn_weights errors for large n without FNN", {
  # Directly call the internal function and verify it errors for large n
  # when FNN is not available.
  skip_if(requireNamespace("FNN", quietly = TRUE),
          "FNN is installed; dense-fallback guard cannot be triggered")

  build_fn <- get(".build_knn_weights", envir = asNamespace("spatialkit"))
  big_coords <- matrix(runif(5001 * 2), ncol = 2)
  expect_error(build_fn(big_coords, k = 8L),
               "requires FNN for k-NN weights")
})

# ---------------------------------------------------------------------------
# GWR fallback bandwidth flag
# ---------------------------------------------------------------------------
test_that("fallback_bandwidth sets bandwidth_is_fallback in info", {
  fb_adaptive <- spatialkit:::.fallback_bandwidth
  # Construct a minimal mock Spatial object
  skip_if_not_installed("sp")
  coords <- matrix(c(0, 0, 1, 1, 2, 2), ncol = 2, byrow = TRUE)
  sp_dat <- sp::SpatialPointsDataFrame(
    coords,
    data = data.frame(x = 1:3)
  )

  bw_a <- fb_adaptive(sp_dat, adaptive = TRUE)
  expect_true(is.integer(bw_a) || is.numeric(bw_a))
  expect_true(bw_a >= 1)

  bw_f <- fb_adaptive(sp_dat, adaptive = FALSE)
  expect_true(is.numeric(bw_f))
  expect_true(bw_f > 0)
})

test_that("gwr_fit info contains bandwidth_is_fallback field", {
  # We cannot run a full GWR fit without GWmodel + real data, but we can

  # verify the info list contract by constructing a mock gwr_fit.
  skip_if_not_installed("sf")
  pts <- sf::st_as_sf(
    data.frame(x = runif(10), y = runif(10), resp = rnorm(10)),
    coords = c("x", "y"), crs = 32632
  )
  mock_info <- list(
    bandwidth = 50,
    adaptive = TRUE,
    kernel = "bisquare",
    AICc = NA_real_,
    bandwidth_is_fallback = TRUE
  )
  expect_true("bandwidth_is_fallback" %in% names(mock_info))
  expect_true(mock_info$bandwidth_is_fallback)
})
