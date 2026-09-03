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


test_that("summarize_by_cell deff=2 inflates SE by more than sqrt(2)", {
  # sqrt(deff) is only HALF the correction.  deff inflates Var(xbar) to
  # sigma^2 * deff / n, but under the same clustering the ordinary sample
  # variance is biased low by exactly the amount that inflation assumes:
  #     E[s^2] = sigma^2 (n - deff) / (n - 1).
  # Applying sqrt(deff) alone therefore leaves the second error in place and
  # the two compound.  Measured 95% CI coverage of the sqrt(deff)-only form at
  # n = 20: 0.905 at rho = 0.3, 0.796 at rho = 0.6, 0.632 at rho = 0.8; with
  # the sqrt((n-1)/(n-deff)) rescale, 0.952 / 0.952 / 0.953 against a nominal
  # 0.95.
  pts <- .make_test_points(rho = 0)
  out_1 <- summarize_by_cell(pts, response_var = "y", deff = 1)
  out_2 <- summarize_by_cell(pts, response_var = "y", deff = 2)
  ratio <- out_2[["..se_resp_y"]] / out_1[["..se_resp_y"]]
  expect_equal(ratio, sqrt(2) * sqrt((out_1$n - 1) / (out_1$n - 2)),
               tolerance = 1e-12)
  # Strictly above the sqrt(deff)-only factor, which is the regression guard.
  expect_true(all(ratio > sqrt(2)))
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


test_that("summarize_by_cell deff='kish' applies Kish's exact 1 + (n-1) rho", {
  # "Kish SE >= IID SE" holds for ANY design effect that grows with rho, so it
  # does not pin the formula: 1 + n*rho, 1 + (n+1)*rho and (1 + rho)^n all pass
  # it.  Kish's deff is specifically
  #
  #     deff_i = 1 + (n_i - 1) * rho,     se_i = sd_i / sqrt(n_i / deff_i),
  #
  # floored at 1.  Cells of DIFFERENT sizes, so an off-by-one in the (n - 1)
  # cannot be absorbed into a rescaled rho: the ratio between the true and a
  # mis-specified deff varies from cell to cell.
  set.seed(4242)
  n_per <- c(2L, 3L, 5L, 9L, 14L, 20L)
  ids   <- rep(seq_along(n_per), times = n_per)
  n     <- length(ids)
  y     <- rep(rnorm(length(n_per), 0, 3), times = n_per) + rnorm(n)
  pts   <- sf::st_sf(
    poly_id  = ids,
    y        = y,
    geometry = sf::st_sfc(lapply(seq_len(n),
                                 function(i) sf::st_point(c(i, ids[i]))),
                          crs = 32632)
  )

  out <- summarize_by_cell(pts, response_var = "y", deff = "kish")
  rho <- attr(out, "deff_applied")$icc_resp
  expect_gt(rho, 0)          # otherwise deff is floored at 1 and pins nothing
  expect_equal(out$n, n_per)

  deff  <- pmax(1, 1 + (out$n - 1) * rho)
  # TWO corrections: sqrt(deff/n) for the mean's inflated variance, and
  # sqrt((n-1)/(n-deff)) because s^2 itself is biased low by the same
  # clustering (E[s^2] = sigma^2 (n - deff)/(n - 1)).  See .se_with_deff().
  se_ok <- out[["..sd_resp_y"]] * sqrt(deff / out$n) *
    sqrt((out$n - 1) / (out$n - deff))
  expect_equal(out[["..se_resp_y"]], se_ok, tolerance = 1e-12)

  # cell_weight is the effective sample size n / deff, from the same formula.
  expect_equal(out$cell_weight, out$n / deff, tolerance = 1e-12)

  # And the (n - 1) inside deff is load-bearing: the nearest plausible
  # alternative gives visibly different standard errors on this design, so the
  # assertion above is not satisfiable by an off-by-one.
  d_off <- pmax(1, 1 + out$n * rho)
  se_off_by_one <- out[["..sd_resp_y"]] * sqrt(d_off / out$n) *
    sqrt((out$n - 1) / pmax(out$n - d_off, .Machine$double.eps))
  expect_false(isTRUE(all.equal(se_ok, se_off_by_one, tolerance = 1e-6)))

  # Regression guard: the sqrt(deff)-only form is strictly smaller wherever
  # deff > 1, which is every cell with more than one observation here.
  se_old <- out[["..sd_resp_y"]] / sqrt(out$n / deff)
  bigger <- deff > 1
  expect_true(any(bigger))
  expect_true(all(se_ok[bigger] > se_old[bigger]))
})

test_that("the Kish standard error covers at the nominal rate", {
  # The property the formula exists for.  The design-effect SE targets
  # Var(cell mean) UNCONDITIONALLY -- the cluster effect is part of the
  # sampling variability -- so the estimand is the grand mean, not the realised
  # cluster mean.  Without the s^2 correction this ran at 0.80 for rho = 0.6.
  skip_on_cran()
  n_cells <- 20L; n_per <- 20L; rho <- 0.6
  hit <- logical(0)
  for (r in 1:150) {
    set.seed(9000 + r)
    cid <- rep(seq_len(n_cells), each = n_per)
    u   <- rnorm(n_cells)
    v   <- sqrt(rho) * u[cid] + sqrt(1 - rho) * rnorm(n_cells * n_per)
    pts <- sf::st_sf(
      poly_id  = cid, y = v,
      geometry = sf::st_sfc(lapply(seq_along(cid),
                                   function(i) sf::st_point(c(i, cid[i]))),
                            crs = 32632))
    out <- summarize_by_cell(pts, response_var = "y", deff = "kish")
    hit <- c(hit, abs(out$resp_mean_y) <=
               stats::qt(0.975, n_per - 1) * out[["..se_resp_y"]])
  }
  expect_gt(mean(hit, na.rm = TRUE), 0.90)
  expect_lt(mean(hit, na.rm = TRUE), 0.99)
})

test_that("summarize_by_cell reports NA spread for a one-observation cell", {
  # A cell holding a single point has no within-cell variation, so its sd and
  # se are undefined rather than zero -- the contract the README leans on when
  # it explains why one-seed-per-point Voronoi cells are an interpolation and
  # not an aggregation.  Zero would read as "measured with perfect precision"
  # and would sail straight into any weighting that reads these columns.
  set.seed(808)
  n_per <- c(1L, 1L, 4L, 6L)
  ids   <- rep(seq_along(n_per), times = n_per)
  n     <- length(ids)
  pts   <- sf::st_sf(
    poly_id  = ids,
    y        = rnorm(n),
    e        = rnorm(n),
    geometry = sf::st_sfc(lapply(seq_len(n),
                                 function(i) sf::st_point(c(i, ids[i]))),
                          crs = 32632)
  )

  for (dv in list(1, "kish")) {
    out  <- summarize_by_cell(pts, response_var = "y", predictor_vars = "e",
                              deff = dv)
    lone <- out$n == 1L
    multi <- !lone
    expect_equal(sum(lone), 2L)
    expect_gt(sum(multi), 0L)

    for (col in c("..sd_resp_y", "..se_resp_y", "..sd_pred_e", "..se_pred_e")) {
      expect_true(all(is.na(out[[col]][lone])),
                  info = paste(col, "deff =", dv))
      # Not NA everywhere -- the columns do carry numbers where there is
      # something to measure, so the assertion above is about the lone cells.
      expect_true(all(is.finite(out[[col]][multi])),
                  info = paste(col, "deff =", dv))
      expect_type(out[[col]], "double")
    }
    # The mean is still reported: it is the SPREAD that is undefined, not the
    # summary itself.
    expect_true(all(is.finite(out$resp_mean_y)))
  }
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

  # --- The pooled value is the (weighted) MEAN of the per-column ICCs -------
  # Each z-scored column contributes its own cells as separate groups to one
  # pooled one-way ANOVA, so with equal cell sizes the pooled ICC is the
  # average of the per-column ICCs.  x_big has essentially no intra-cell
  # correlation, so the pooled value sits near rho_small / 2.
  #
  # An earlier version of this test bracketed rho_small / 4 and derived it
  # from pooling the two columns under the SAME cell labels.  That derivation
  # was correct about the code and wrong about the statistic: independent
  # per-column cell effects average away in a shared cell mean, and the
  # between-cell sum of squares shrinks by about 1/m.  Measured with four
  # columns at true rho = 0.5, the old pooled ICC was 0.12 (per-column ANOVA
  # 0.495), so every predictor SE was ~44% too small.
  out_big <- summarize_by_cell(pts, response_var = "y",
                               predictor_vars = "x_big", deff = "kish")
  icc_big <- attr(out_big, "deff_applied")$icc_pred
  if (is.null(icc_big)) icc_big <- 0            # deff not applied when ICC <= 0
  expect_equal(da_both$icc_pred, (da_small$icc_pred + icc_big) / 2,
               tolerance = 0.06)
  # Regression guard against the shared-label pooling: that gave ~rho/4.
  expect_gt(da_both$icc_pred, 0.35 * da_small$icc_pred)

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

  # The ..se_ closures run with deff = 1 and the result is rescaled afterwards,
  # so this equality is the whole mechanism.  The factor is the same one
  # .se_with_deff() applies on the Kish path -- sqrt(deff) for the mean's
  # inflated variance AND sqrt((n-1)/(n-deff)) for the downward bias in s^2 --
  # so the two paths cannot drift apart.
  iid  <- summarize_by_cell(pts, response_var = "z", deff = 1)
  infl <- sqrt(da$deff) * sqrt((out$n - 1) / (out$n - da$deff))
  expect_equal(out[["..se_resp_z"]], iid[["..se_resp_z"]] * infl,
               tolerance = 1e-10)
  # Regression guard: strictly above the sqrt(deff)-only rescale.
  expect_true(all(infl > sqrt(da$deff)))
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
  # `deff_max_n` caps the SUBSAMPLE the correlation matrix is built on, not the
  # answer.  The design effect reported is still the CELL's -- estimated as
  # 1 + (n_i - 1) * Rbar with Rbar taken from the subsample -- so it is bounded
  # by the cell size, not by max_n, and stays close to the un-subsampled value.
  # The old code returned sum(R)/n_used, i.e. the design effect of a cell of
  # max_n points, understating a large cell by roughly n_i/max_n.
  n_cell <- capped$n
  expect_true(all(da_cap$deff <= n_cell))
  expect_true(any(da_cap$deff > 10))           # the old form could not exceed 10
  expect_equal(unname(da_cap$deff), unname(attr(full, "deff_applied")$deff),
               tolerance = 0.35)
  # A cap larger than every cell changes nothing at all.
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


# ---------------------------------------------------------------------------
# Fourth-pass regressions: design effects are built from what a column
# actually observes, in the CRS the variogram was fitted in.
# ---------------------------------------------------------------------------

.vgm_na_fixture <- function() {
  # 4 cells x 10 points on a projected grid; a fitted exponential variogram
  # with a 60-unit range so within-cell pairs are strongly correlated.
  set.seed(31)
  ids <- rep(1:4, each = 10)
  xy  <- cbind(runif(40, 0, 100) + 500 * ((ids - 1) %% 2),
               runif(40, 0, 100) + 500 * ((ids - 1) %/% 2))
  pts <- sf::st_sf(poly_id = ids, val = rnorm(40),
                   geometry = sf::st_sfc(lapply(seq_len(40), function(i)
                     sf::st_point(xy[i, ])), crs = 32632))
  vgm <- gstat::vgm(psill = 1, model = "Exp", range = 20, nugget = 0)
  sac <- structure(60, class = c("sac_range", "numeric"),
                   variogram_model = vgm, crs = sf::st_crs(32632))
  list(pts = pts, sac = sac)
}

test_that("NA-response rows do not move the variogram design effect or SE", {
  # A row with no response contributes nothing to the response mean, so it
  # cannot contribute to that mean's clustering either.  Previously the cell's
  # deff was formed at the ROW count and the SE rescaled by it: adding 8
  # NA-response rows to a 2-observation cell moved its SE from 8.46 to 26.38.
  skip_if_not_installed("gstat")
  fx <- .vgm_na_fixture()
  a  <- suppressMessages(summarize_by_cell(fx$pts, response_var = "val",
                                           deff = "variogram", sac = fx$sac))
  # Blank out 8 of the 10 responses in cell 1 -- the locations stay.
  p2 <- fx$pts; p2$val[p2$poly_id == 1][3:10] <- NA
  b  <- suppressMessages(summarize_by_cell(p2, response_var = "val",
                                           deff = "variogram", sac = fx$sac))
  # Reference: the same 2 observed rows as their own 2-point cell.
  p3 <- fx$pts[!(fx$pts$poly_id == 1 & seq_len(40) %in% 3:10), ]
  cc <- suppressMessages(summarize_by_cell(p3, response_var = "val",
                                           deff = "variogram", sac = fx$sac))
  i_b <- which(b$poly_id == 1); i_c <- which(cc$poly_id == 1)
  expect_equal(b[["..se_resp_val"]][i_b], cc[["..se_resp_val"]][i_c], tolerance = 1e-10)
  expect_equal(b$cell_weight[i_b],        cc$cell_weight[i_c],        tolerance = 1e-10)
  expect_equal(attr(b, "deff_applied")$deff[i_b], attr(cc, "deff_applied")$deff[i_c],
               tolerance = 1e-10)
  # `n` still counts rows; cell_weight counts what the response observed.
  expect_equal(b$n[i_b], 10L)
  expect_lt(b$cell_weight[i_b], 2 + 1e-9)
  # Untouched cells are identical between the two runs.
  expect_equal(b[["..se_resp_val"]][-i_b], a[["..se_resp_val"]][-which(a$poly_id == 1)])
})

test_that("cell_weight is the effective count of the response, not of rows", {
  set.seed(5)
  d <- data.frame(poly_id = rep(1:3, each = 10), val = rnorm(30),
                  x = runif(30, 0, 100), y = runif(30, 0, 100))
  d$val[d$poly_id == 1][1:7] <- NA           # 3 finite responses in cell 1
  pts <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
  out <- summarize_by_cell(pts, response_var = "val", deff = 2)
  i <- which(out$poly_id == 1)
  expect_equal(out$n[i], 10L)
  expect_equal(out$cell_weight[i], 3 / 2)     # was 10 / 2 = 5
  expect_equal(out$cell_weight[-i], rep(10 / 2, 2))
})

test_that("the variogram design effect is evaluated in the variogram's CRS", {
  # A range is a length in the CRS the variogram was fitted in (metres).
  # Evaluating the correlation function at distances taken from lon/lat input
  # (degrees) saturated every pair: measured SE ratio of 354 between UTM and
  # EPSG:4326 input for the SAME points, cells and variogram.
  skip_if_not_installed("gstat")
  fx <- .vgm_na_fixture()
  utm <- fx$pts
  geo <- sf::st_transform(utm, 4326)
  a <- suppressMessages(summarize_by_cell(utm, response_var = "val",
                                          deff = "variogram", sac = fx$sac))
  b <- suppressMessages(summarize_by_cell(geo, response_var = "val",
                                          deff = "variogram", sac = fx$sac))
  expect_equal(b[["..se_resp_val"]], a[["..se_resp_val"]], tolerance = 1e-6)
  expect_equal(attr(b, "deff_applied")$deff, attr(a, "deff_applied")$deff,
               tolerance = 1e-6)
  # And the recorded CRS is the variogram's, whatever the input was in.
  expect_equal(attr(b, "deff_applied")$crs, sf::st_crs(32632))
})

test_that("deff = 'kish' tolerates rows with no cell id", {
  # keep_unassigned = TRUE output carries NA ids; those rows belong to no
  # group and must be left out of the ICC rather than crash it.
  set.seed(8)
  d <- data.frame(poly_id = c(rep(1:3, each = 8), NA, NA), val = rnorm(26),
                  x = runif(26, 0, 100), y = runif(26, 0, 100))
  pts <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
  expect_no_error(out <- suppressMessages(
    summarize_by_cell(pts, response_var = "val", deff = "kish")))
  expect_true(all(is.finite(out[["..se_resp_val"]][!is.na(out$poly_id)])))
})

test_that("the pooled ICC groups by (variable, cell), not by cell alone", {
  # Pooling m z-scored columns under one cell label averages independent cell
  # effects away and shrinks the ICC by ~1/m.  With m = 4 columns at a common
  # true rho, the pooled estimate must sit near rho, not rho / 4.
  set.seed(77)
  n_cells <- 40; n_per <- 12; rho <- 0.5
  ids <- rep(seq_len(n_cells), each = n_per)
  mk  <- function() sqrt(rho) * rnorm(n_cells)[ids] + sqrt(1 - rho) * rnorm(n_cells * n_per)
  d <- data.frame(poly_id = ids, y = rnorm(length(ids)),
                  p1 = mk(), p2 = mk(), p3 = mk(), p4 = mk(),
                  x = runif(length(ids), 0, 100), y0 = runif(length(ids), 0, 100))
  pts <- sf::st_as_sf(d, coords = c("x", "y0"), crs = 32632)
  out <- suppressMessages(summarize_by_cell(pts, response_var = "y",
                                            predictor_vars = c("p1", "p2", "p3", "p4"),
                                            deff = "kish"))
  icc <- attr(out, "deff_applied")$icc_pred
  expect_gt(icc, 0.35)                        # rho/4 = 0.125 would fail this
  expect_lt(icc, 0.65)
})
