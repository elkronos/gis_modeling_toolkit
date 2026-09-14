# tests/testthat/test-select-on.R
# ---------------------------------------------------------------------------
# select_on = "split": choose on one spatially blocked half, estimate on the
# other, so what is estimated afterwards is not post-selection.
# ---------------------------------------------------------------------------

so_field <- function(n = 400, seed = 2) {
  set.seed(seed)
  xy <- data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000))
  D  <- as.matrix(stats::dist(xy))
  S  <- as.numeric(t(chol(exp(-D / 100) + diag(1e-8, n))) %*% rnorm(n))
  xy$w <- rnorm(n); xy$noise <- rnorm(n)
  xy$z <- 2 * xy$w + S + rnorm(n, sd = 0.5)
  sf::st_as_sf(xy, coords = c("x", "y"), crs = 32632)
}
so_fit <- function(tr, vars) lm_spatial_fit(tr, "z", vars)

test_that(".spatial_half_split makes two disjoint, exhaustive, spatially blocked halves", {
  pts <- so_field(300)
  sp  <- spatialkit:::.spatial_half_split(pts, seed = 5, caller = "test")
  expect_s3_class(sp, "spatialkit_split")
  expect_length(intersect(sp$selection, sp$estimation), 0L)
  expect_setequal(c(sp$selection, sp$estimation), seq_len(300))
  expect_identical(sp$method, "block_kfold")
  # Blocked, not interleaved: the two halves are the two folds of
  # make_folds(k = 2, method = "block_kfold") under the same seed.
  f <- make_folds(pts, k = 2, method = "block_kfold", seed = 5)
  expect_identical(sp$selection, sort(as.integer(f$assignment$row_id[f$assignment$fold == 1L])))
  # Reproducible; a different seed gives a different split.
  expect_identical(spatialkit:::.spatial_half_split(pts, seed = 5, caller = "t")$selection,
                   sp$selection)
  expect_false(identical(spatialkit:::.spatial_half_split(pts, seed = 6, caller = "t")$selection,
                         sp$selection))
  # A pre-existing ..row_id column does not change the positions returned.
  pts2 <- pts; pts2$..row_id <- 1000L + seq_len(300)
  expect_identical(spatialkit:::.spatial_half_split(pts2, seed = 5, caller = "t")$selection,
                   sp$selection)
  expect_error(spatialkit:::.spatial_half_split(pts[1:10, ], caller = "f"),
               "at least 20 points")
  expect_output(print(sp), "Spatial half-split")
})

test_that("determine_optimal_levels(select_on = 'split') selects on one half and returns both", {
  pts <- so_field(300)
  out <- determine_optimal_levels(pts, max_levels = 8, response_var = "z",
                                  predictor_vars = "w", select_on = "split")
  sp <- attr(out, "split")
  expect_s3_class(sp, "spatialkit_split")
  expect_setequal(c(sp$selection, sp$estimation), seq_len(300))
  # The selection is the one made on the selection half alone.
  half <- determine_optimal_levels(pts[sp$selection, ], max_levels = 8,
                                   response_var = "z", predictor_vars = "w")
  expect_identical(as.integer(out), as.integer(half))
  # The default path is unchanged and carries no split.
  all_pts <- determine_optimal_levels(pts, max_levels = 8, response_var = "z",
                                      predictor_vars = "w")
  expect_null(attr(all_pts, "split"))
  # The geometric path with a split: same integer vector, plus the attribute.
  geo <- determine_optimal_levels(pts, max_levels = 6, select_on = "split")
  expect_type(geo, "integer")
  expect_s3_class(attr(geo, "split"), "spatialkit_split")
  expect_identical(as.integer(geo),
                   as.integer(determine_optimal_levels(pts[attr(geo, "split")$selection, ],
                                                       max_levels = 6)))
  expect_error(determine_optimal_levels(pts, select_on = "half"), "'arg' should be one of")
})

test_that("resolution_profile(select_on = 'split') profiles the selection half", {
  skip_if_not_installed("gstat")
  pts <- so_field(400)
  prof <- resolution_profile(pts, response_var = "z", predictor_vars = "w",
                             n_levels = 5, select_on = "split")
  sp <- attr(prof, "split")
  expect_s3_class(sp, "spatialkit_split")
  expect_setequal(c(sp$selection, sp$estimation), seq_len(400))
  expect_identical(attr(prof, "bounds")$n, length(sp$selection))
  expect_output(print(prof), "estimate on the other")
  expect_null(attr(resolution_profile(pts, n_levels = 4), "split"))
})

test_that("select_features_forward(select_on = 'split') scores the selected set on the hold-out half", {
  pts <- so_field(400)
  sel <- select_features_forward(pts, "z", c("w", "noise"), so_fit, k = 3,
                                 seed = 1, quiet = TRUE, select_on = "split")
  expect_identical(sel$params$select_on, "split")
  sp <- sel$split
  expect_s3_class(sp, "spatialkit_split")
  expect_setequal(c(sp$selection, sp$estimation), seq_len(400))
  expect_true("w" %in% sel$selected)
  expect_true(is.finite(sel$score))
  expect_true(is.finite(sel$score_holdout))
  # score_holdout is exactly the selected set fitted on the selection half
  # and scored on the estimation half against the selection half's mean.
  fit <- so_fit(pts[sp$selection, ], sel$selected)
  yh  <- as.numeric(predict(fit, newdata = pts[sp$estimation, ]))
  met <- spatialkit:::.compute_reg_metrics(pts$z[sp$estimation], yh,
                                           y_train_mean = mean(pts$z[sp$selection]))
  expect_equal(sel$score_holdout, met$RMSE)
  # The sweep itself ran on the selection half: same result as calling it on
  # that half directly.
  direct <- select_features_forward(pts[sp$selection, ], "z", c("w", "noise"), so_fit,
                                    k = 3, seed = 1, quiet = TRUE)
  expect_identical(sel$selected, direct$selected)
  expect_equal(sel$score, direct$score)
  # Default: no split, no hold-out score.
  all_sel <- select_features_forward(pts, "z", c("w", "noise"), so_fit, k = 3,
                                     seed = 1, quiet = TRUE)
  expect_null(all_sel$split)
  expect_true(is.na(all_sel$score_holdout))
  expect_identical(all_sel$params$select_on, "all")
})

test_that("a hold-out score that cannot be computed is NA with a warning, not an error", {
  pts <- so_field(300)
  bad_fit <- function(tr, vars) {
    fit <- lm_spatial_fit(tr, "z", vars)
    class(fit) <- c("so_broken", class(fit))
    fit
  }
  # Predicts inside the sweep (test folds of ~65 rows) and fails on the
  # hold-out half (~150 rows), so the selection succeeds and only the
  # hold-out score is lost.
  registerS3method("predict", "so_broken",
                   function(object, newdata = NULL, ...) {
                     if (!is.null(newdata) && nrow(newdata) > 100L) stop("no predictions today")
                     predict.lmsurf_fit(object, newdata = newdata, ...)
                   })
  expect_warning(
    sel <- select_features_forward(pts, "z", c("w", "noise"), bad_fit, k = 3,
                                   seed = 1, quiet = TRUE, select_on = "split"),
    "hold-out score could not be computed")
  expect_true(is.na(sel$score_holdout))
  expect_true(length(sel$selected) >= 1L)
})
