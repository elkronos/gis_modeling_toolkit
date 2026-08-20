# ===========================================================================
# summarize_by_cell(): per-cell aggregation and its design-effect options.
# See test-deff-variogram.R for the variogram-based design effect.
# ===========================================================================

# A simple sf with a cell ID and a numeric variable.  When rho > 0 the
# observations within a cell share a common component, mimicking intra-cell
# spatial autocorrelation -- which is what a design effect is meant to absorb.
.make_test_points <- function(n_cells = 5, pts_per_cell = 10, rho = 0) {
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
