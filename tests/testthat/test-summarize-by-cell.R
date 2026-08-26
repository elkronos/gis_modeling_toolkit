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

  # --- The property z-scoring exists for: unit invariance -------------------
  # Re-express either predictor in units 10^6 times larger.  Standardising
  # before pooling makes the pooled ICC EXACTLY unchanged.  Without it the
  # rescaled column carries essentially all the pooled variance, so the answer
  # would be that column's own ICC.
  .rescaled <- function(small_mult, big_mult) {
    p <- pts
    p$x_small <- p$x_small * small_mult
    p$x_big   <- p$x_big   * big_mult
    attr(summarize_by_cell(p, response_var = "y",
                           predictor_vars = c("x_small", "x_big"),
                           deff = "kish"), "deff_applied")$icc_pred
  }
  expect_equal(.rescaled(1, 1e6), da_both$icc_pred, tolerance = 1e-10)
  expect_equal(.rescaled(1e6, 1), da_both$icc_pred, tolerance = 1e-10)

  # --- Two-sided bracket on the pooled value --------------------------------
  # Derivation.  Both columns are z-scored, so each contributes unit total
  # variance and the cross term vanishes for independent columns:
  #     SSB_pooled = (SSB_small + SSB_big) / 2  ~  SSB_small / 2
  # while the cluster size n0 doubles (each cell now holds 2 * n_per values).
  # The variance component is sigma_b^2 = (MSB - MSW) / n0, so halving the
  # numerator and doubling n0 divides it by ~4; with sigma_w^2 ~ 1 the pooled
  # ICC lands near rho_small / 4.  The bracket spans 0.15-0.40 of it: an
  # un-z-scored pooling gives ~0 (below), and ignoring x_big entirely gives
  # 1.0 (above).
  expect_gte(da_both$icc_pred, 0.15 * da_small$icc_pred)
  expect_lte(da_both$icc_pred, 0.40 * da_small$icc_pred)

  # And it is genuinely pulled by BOTH: strictly under the small-only ICC.
  expect_lt(da_both$icc_pred, da_small$icc_pred)
})


test_that("summarize_by_cell invalid deff falls back to 1", {
  pts <- .make_test_points()

  # The fallback is logged through logger, which raises no R condition, so it
  # has to be captured (see helper-logging.R) -- expect_warning() would pass
  # whether or not anything was emitted.
  lines <- capture_spatialkit_log(
    out <- summarize_by_cell(pts, response_var = "y", deff = -1)
  )
  expect_true(log_has(lines, "deff must be >= 1"))

  # ... and the fallback really is deff = 1: weights untouched, SE the classic
  # sd / sqrt(n), and no design effect recorded.
  expect_equal(out$cell_weight, out$n)
  expect_equal(out[["..se_resp_y"]], out[["..sd_resp_y"]] / sqrt(out$n),
               tolerance = 1e-12)
  expect_null(attr(out, "deff_applied"))

  # A non-numeric, non-keyword value takes the same branch.
  lines2 <- capture_spatialkit_log(
    out2 <- summarize_by_cell(pts, response_var = "y", deff = "nonsense")
  )
  expect_true(log_has(lines2, "deff must be >= 1"))
  expect_null(attr(out2, "deff_applied"))
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
# deff = "variogram", end to end
#
# The internals are covered in isolation by test-deff-variogram.R; what is
# asserted here is the WIRING -- that summarize_by_cell() builds the
# correlation function from the supplied fit, passes deff_max_n through, and
# applies the per-cell design effect to the SEs and the cell weights.
# ---------------------------------------------------------------------------

# Six cells strung out along x, each holding 25 points drawn from a field with
# a 40-unit exponential correlation range, so within-cell correlation is
# strong and the design effect is well above 1.
.vgm_cell_points <- function(n_cells = 6L, per = 25L, seed = 77) {
  set.seed(seed)
  n  <- n_cells * per
  cx <- rep(seq(0, 500, length.out = n_cells), each = per)
  x  <- cx + runif(n, 0, 60)
  y  <- runif(n, 0, 60)
  d  <- as.matrix(stats::dist(cbind(x, y)))
  z  <- as.numeric(t(chol(exp(-d / 40) + diag(1e-6, n))) %*% rnorm(n))
  sf::st_as_sf(
    data.frame(x = x, y = y, z = z,
               poly_id = rep(seq_len(n_cells), each = per)),
    coords = c("x", "y"), crs = 32632)
}

# sum(R) / n with R from the fitted variogram: the definition the per-cell
# design effect is supposed to implement, computed here from the coordinates
# directly.
.hand_deff <- function(pts, sac, id_col = "poly_id") {
  cor_fn <- spatialkit:::.vgm_correlation_fn(attr(sac, "variogram_model"))
  xy  <- sf::st_coordinates(pts)[, 1:2, drop = FALSE]
  ids <- unique(sf::st_drop_geometry(pts)[[id_col]])
  vapply(ids, function(id) {
    idx <- which(sf::st_drop_geometry(pts)[[id_col]] == id)
    d   <- as.matrix(stats::dist(xy[idx, , drop = FALSE]))
    R   <- matrix(as.numeric(cor_fn(as.numeric(d))), nrow(d), ncol(d))
    diag(R) <- 1
    min(max(sum(R) / length(idx), 1), length(idx))
  }, numeric(1))
}

test_that("summarize_by_cell(deff = 'variogram') applies sum(R)/n per cell", {
  skip_if_not_installed("gstat")
  pts <- .vgm_cell_points()
  sac <- suppressWarnings(estimate_sac_range(pts, "z"))
  expect_false(is.null(attr(sac, "variogram_model")))

  out <- summarize_by_cell(pts, response_var = "z", deff = "variogram",
                           sac = sac)
  da  <- attr(out, "deff_applied")

  expect_equal(da$method, "variogram")
  expect_equal(da$max_n, 500L)
  expect_length(da$deff, nrow(out))
  # The per-cell design effects are the hand-computed sum(R)/n, not a constant
  # and not Kish's 1 + (n-1)*rho.
  expect_equal(unname(da$deff), unname(.hand_deff(pts, sac)), tolerance = 1e-8)
  expect_true(all(da$deff > 1))

  # SE scales as sqrt(deff): the ..se_ closures run with deff = 1 and the
  # result is rescaled afterwards, so this equality is the whole mechanism.
  iid <- summarize_by_cell(pts, response_var = "z", deff = 1)
  expect_equal(out[["..se_resp_z"]], iid[["..se_resp_z"]] * sqrt(da$deff),
               tolerance = 1e-10)
  # ... and the effective sample size is n / deff.
  expect_equal(out$cell_weight, out$n / da$deff, tolerance = 1e-10)
  expect_true(all(out$cell_weight < out$n))

  # The means and SDs are untouched by the design effect.
  expect_equal(out[["resp_mean_z"]], iid[["resp_mean_z"]])
  expect_equal(out[["..sd_resp_z"]], iid[["..sd_resp_z"]])
})


test_that("summarize_by_cell passes deff_max_n through to the subsampler", {
  skip_if_not_installed("gstat")
  pts <- .vgm_cell_points()
  sac <- suppressWarnings(estimate_sac_range(pts, "z"))

  full    <- summarize_by_cell(pts, response_var = "z", deff = "variogram",
                               sac = sac)
  capped  <- summarize_by_cell(pts, response_var = "z", deff = "variogram",
                               sac = sac, deff_max_n = 10L)
  da_cap  <- attr(capped, "deff_applied")

  expect_equal(da_cap$max_n, 10L)
  # A design effect cannot exceed the number of points it was computed on, so
  # capping the correlation matrix at 10 caps every cell's deff at 10.
  expect_true(all(da_cap$deff <= 10))
  expect_true(all(attr(full, "deff_applied")$deff > 10))
  # A cap larger than every cell changes nothing.
  wide <- summarize_by_cell(pts, response_var = "z", deff = "variogram",
                            sac = sac, deff_max_n = 5000L)
  expect_equal(attr(wide, "deff_applied")$deff,
               attr(full, "deff_applied")$deff)
})


test_that("summarize_by_cell falls back to deff = 1 with no fitted variogram", {
  # Reachable whenever gstat is unavailable, and also -- deterministically,
  # without uninstalling anything -- when nothing supplies a variogram to fit.
  pts <- .vgm_cell_points()

  lines <- capture_spatialkit_log(
    out <- summarize_by_cell(pts, predictor_vars = "z", deff = "variogram")
  )
  expect_true(log_has(lines, "requires a fitted variogram model"))
  expect_true(log_has(lines, "Falling back to deff = 1"))
  expect_null(attr(out, "deff_applied"))

  iid <- summarize_by_cell(pts, predictor_vars = "z", deff = 1)
  expect_equal(out[["..se_pred_z"]], iid[["..se_pred_z"]])
  expect_equal(out$cell_weight, out$n)

  # Same when the estimation attempt itself fails (a broken or absent gstat).
  lines2 <- capture_spatialkit_log({
    local_mocked_bindings(
      estimate_sac_range = function(...) stop("gstat unavailable"))
    out2 <- summarize_by_cell(pts, response_var = "z", deff = "variogram")
  })
  expect_true(log_has(lines2, "requires a fitted variogram model"))
  expect_null(attr(out2, "deff_applied"))
})


test_that("summarize_by_cell(deff = 'variogram') coerces non-POINT geometry first", {
  # st_coordinates() returns one row per VERTEX, so a MULTIPOINT layer produced
  # a coordinate matrix twice as tall as the id vector and fed the wrong points
  # into every cell's correlation matrix.
  skip_if_not_installed("gstat")
  pts <- .vgm_cell_points()
  xy  <- sf::st_coordinates(pts)

  # Each feature becomes a two-vertex MULTIPOINT centred on its original
  # position, so the coerced centroids are the original points exactly.
  multi <- pts
  sf::st_geometry(multi) <- sf::st_sfc(
    lapply(seq_len(nrow(pts)), function(i) {
      sf::st_multipoint(rbind(c(xy[i, 1] - 3, xy[i, 2]),
                              c(xy[i, 1] + 3, xy[i, 2])))
    }), crs = 32632)

  coerced <- coerce_to_points(multi, "auto")
  expect_equal(nrow(sf::st_coordinates(coerced)), nrow(pts))
  expect_equal(unname(sf::st_coordinates(coerced)[, 1:2]), unname(xy[, 1:2]),
               tolerance = 1e-6)

  sac <- suppressWarnings(estimate_sac_range(pts, "z"))
  point_out <- summarize_by_cell(pts, response_var = "z", deff = "variogram",
                                 sac = sac)
  multi_out <- summarize_by_cell(multi, response_var = "z", deff = "variogram",
                                 sac = sac)

  expect_equal(attr(multi_out, "deff_applied")$deff,
               attr(point_out, "deff_applied")$deff, tolerance = 1e-8)
  expect_equal(multi_out$cell_weight, point_out$cell_weight, tolerance = 1e-10)
  expect_equal(multi_out[["..se_resp_z"]], point_out[["..se_resp_z"]],
               tolerance = 1e-10)
})


test_that("the variogram deff survives the cells_sf join, realigned per cell", {
  skip_if_not_installed("gstat")
  pts <- .vgm_cell_points(n_cells = 4L, per = 20L)
  sac <- suppressWarnings(estimate_sac_range(pts, "z"))

  plain <- summarize_by_cell(pts, response_var = "z", deff = "variogram",
                             sac = sac)

  # Cell polygons, including one with NO observations, in a different order
  # from the summarised rows.
  bb <- sf::st_bbox(pts)
  mk <- function(i) {
    x0 <- as.numeric(bb["xmin"]) + (i - 1) * 200
    sf::st_polygon(list(rbind(c(x0, -100), c(x0 + 200, -100),
                              c(x0 + 200, 400), c(x0, 400), c(x0, -100))))
  }
  cells <- sf::st_sf(poly_id = c(4L, 3L, 2L, 1L, 99L),
                     geometry = sf::st_sfc(lapply(1:5, mk), crs = 32632))

  joined <- summarize_by_cell(pts, response_var = "z", deff = "variogram",
                              sac = sac, cells_sf = cells)
  da <- attr(joined, "deff_applied")

  # dplyr::left_join() rebuilds attributes from its `x` template, which used to
  # drop "deff_applied" entirely.
  expect_false(is.null(da))
  expect_equal(da$method, "variogram")
  expect_s3_class(joined, "sf")
  expect_equal(nrow(joined), 5L)

  # The per-cell vector is realigned to the JOINED row order, with NA for the
  # cell that has no observations.
  expect_length(da$deff, 5L)
  lookup <- stats::setNames(attr(plain, "deff_applied")$deff, plain$poly_id)
  expect_equal(unname(da$deff), unname(lookup[as.character(joined$poly_id)]))
  expect_true(is.na(da$deff[joined$poly_id == 99L]))
  expect_false(anyNA(da$deff[joined$poly_id != 99L]))
})
