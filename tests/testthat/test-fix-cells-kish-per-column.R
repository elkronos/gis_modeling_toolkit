# tests/testthat/test-fix-cells-kish-per-column.R
# ---------------------------------------------------------------------------
# summarize_by_cell(deff = "kish"): one ICC per summarised column.
#
# One ICC used to be pooled over every predictor listed, so a column's
# standard error depended on which other predictors were named beside it.
# On 60 cells of 12 with an ICC-0 and an ICC-0.87 predictor the pooled 0.42
# multiplied both standard errors by 3.13, where their own ICCs give 1.00
# and 8.96; with no response, `cell_weight` was built from the pooled value
# as well.
# ---------------------------------------------------------------------------

# 60 cells of 12: a white-noise predictor, a strongly clustered one, and an
# unclustered response.
.fk_points <- function(seed = 3, k = 60L, per = 12L) {
  set.seed(seed)
  n <- k * per
  cell <- rep(seq_len(k), each = per)
  d <- data.frame(x = 4e5 + runif(n, 0, 1e5), y = 4.6e6 + runif(n, 0, 1e5),
                  poly_id = cell)
  d$resp  <- rnorm(n)
  d$noise <- rnorm(n)
  d$clust <- rnorm(k)[cell] + rnorm(n, 0, 0.33)
  sf::st_as_sf(d, coords = c("x", "y"), crs = 32615)
}

# The same ANOVA estimator on one column, by hand.
.fk_icc <- function(v, g) {
  g <- factor(g); ni <- tabulate(g); k <- nlevels(g); N <- length(v)
  gm <- tapply(v, g, mean)
  msb <- sum(ni * (gm - mean(v))^2) / (k - 1)
  msw <- sum((v - gm[g])^2) / (N - k)
  n0 <- (N - sum(ni^2) / N) / (k - 1)
  max(0, min(1, (msb - msw) / (msb + (n0 - 1) * msw)))
}

test_that("each predictor's standard error, effective n and interval use its own ICC", {
  pts <- .fk_points()
  both   <- summarize_by_cell(pts, "resp", c("noise", "clust"), deff = "kish",
                              conf_level = 0.95, estimate = "mean")
  aloneN <- summarize_by_cell(pts, "resp", "noise", deff = "kish", conf_level = 0.95,
                              estimate = "mean")
  aloneC <- summarize_by_cell(pts, "resp", "clust", deff = "kish", conf_level = 0.95,
                              estimate = "mean")

  # Adding or removing a predictor leaves the other's columns as they were.
  for (stat in c("..se_", "..neff_", "..df_", "..ci_lo_", "..ci_hi_")) {
    expect_equal(both[[paste0(stat, "pred_noise")]], aloneN[[paste0(stat, "pred_noise")]],
                 tolerance = 1e-12, info = stat)
    expect_equal(both[[paste0(stat, "pred_clust")]], aloneC[[paste0(stat, "pred_clust")]],
                 tolerance = 1e-12, info = stat)
  }
  # In either order.
  flipped <- summarize_by_cell(pts, "resp", c("clust", "noise"), deff = "kish",
                               conf_level = 0.95, estimate = "mean")
  expect_equal(flipped[["..se_pred_noise"]], both[["..se_pred_noise"]], tolerance = 1e-12)
  expect_equal(flipped[["..se_pred_clust"]], both[["..se_pred_clust"]], tolerance = 1e-12)
  # The response is untouched by any of it.
  only <- summarize_by_cell(pts, "resp", deff = "kish", conf_level = 0.95, estimate = "mean")
  expect_equal(both[["..se_resp_resp"]], only[["..se_resp_resp"]], tolerance = 1e-12)
  expect_equal(both$cell_weight, only$cell_weight, tolerance = 1e-12)

  # The numbers: the unclustered predictor's SE is the naive one, the
  # clustered one's is about nine times it.  Pooled, both were 3.13 times.
  naive <- summarize_by_cell(pts, "resp", c("noise", "clust"), estimate = "mean")
  expect_equal(both[["..se_pred_noise"]], naive[["..se_pred_noise"]], tolerance = 1e-12)
  ratio_c <- stats::median(both[["..se_pred_clust"]] / naive[["..se_pred_clust"]])
  expect_gt(ratio_c, 7)
  expect_lt(ratio_c, 11)
})

test_that("the ICCs are recorded per predictor column, by name", {
  pts <- .fk_points()
  d <- sf::st_drop_geometry(pts)
  both <- summarize_by_cell(pts, "resp", c("noise", "clust"), deff = "kish", estimate = "mean")
  icc <- attr(both, "icc")
  expect_named(icc$pred, c("noise", "clust"))
  expect_equal(icc$pred[["noise"]], .fk_icc(d$noise, d$poly_id), tolerance = 1e-10)
  expect_equal(icc$pred[["clust"]], .fk_icc(d$clust, d$poly_id), tolerance = 1e-10)
  expect_lt(icc$pred[["noise"]], 0.02)
  expect_gt(icc$pred[["clust"]], 0.8)
  expect_length(icc$resp, 1L)
  # The same vector on the applied record.
  da <- attr(both, "deff_applied")
  expect_identical(da$icc_pred, icc$pred)
  expect_identical(da$method, "kish")
  # One predictor: still named.
  one <- summarize_by_cell(pts, "resp", "clust", deff = "kish", estimate = "mean")
  expect_named(attr(one, "icc")$pred, "clust")
  # No predictor, no response column of that type: NA, as before.
  expect_identical(attr(summarize_by_cell(pts, "resp", deff = "kish", estimate = "mean"),
                        "icc")$pred,
                   NA_real_)
})

test_that("with no response, cell_weight uses the first predictor's own ICC", {
  pts <- .fk_points()
  # Primary variable `noise` (ICC 0): its weight is its count, whatever else
  # is listed.  It fell from 12 to 2.12 when `clust` was listed after it.
  nb <- summarize_by_cell(pts, predictor_vars = c("noise", "clust"), deff = "kish",
                          estimate = "mean")
  n1 <- summarize_by_cell(pts, predictor_vars = "noise", deff = "kish", estimate = "mean")
  expect_equal(nb$cell_weight, n1$cell_weight, tolerance = 1e-12)
  expect_equal(nb$cell_weight, rep(12, 60))
  # The clustered predictor is still corrected, and the rows say a
  # correction was applied.
  expect_true(all(nb$deff_applied))
  expect_identical(unique(nb$deff_method), "kish")
  # Primary variable `clust`: its own, strong, ICC sets the weight.
  cb <- summarize_by_cell(pts, predictor_vars = c("clust", "noise"), deff = "kish",
                          estimate = "mean")
  c1 <- summarize_by_cell(pts, predictor_vars = "clust", deff = "kish", estimate = "mean")
  expect_equal(cb$cell_weight, c1$cell_weight, tolerance = 1e-12)
  expect_true(all(cb$cell_weight < 2))
})

test_that("a column with missing values keeps its own ICC and count", {
  pts <- .fk_points()
  pts$clust[c(1, 2, 14, 30)] <- NA
  both  <- summarize_by_cell(pts, "resp", c("noise", "clust"), deff = "kish",
                             conf_level = 0.9, estimate = "mean")
  alone <- summarize_by_cell(pts, "resp", "clust", deff = "kish", conf_level = 0.9,
                             estimate = "mean")
  expect_equal(both[["..se_pred_clust"]], alone[["..se_pred_clust"]], tolerance = 1e-12)
  expect_equal(both[["..neff_pred_clust"]], alone[["..neff_pred_clust"]], tolerance = 1e-12)
})
