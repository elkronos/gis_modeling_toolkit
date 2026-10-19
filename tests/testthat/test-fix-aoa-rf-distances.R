# tests/testthat/test-fix-aoa-rf-distances.R
# ---------------------------------------------------------------------------
# Two defects in the reference distances of area_of_applicability().
#
# (a) The dense path (no FNN) built squared distances as
#     |a|^2 + |b|^2 - 2ab, which for an exact copy of a reference row is a
#     few ulp away from 0: clamped to 0 for most copies, about 1e-8 for the
#     rest.  Against a threshold at or near 0 that put exact copies of
#     training rows outside the AOA and kept a threshold that should have
#     been 0 just above it.  The distance returned is now recomputed
#     exactly, so the two paths agree.
#
# (b) With folds = NULL a training row that repeats another's predictor
#     values was its "nearest other row" at distance 0.  Sixty sites visited
#     four times gave a threshold of 0 and 30 of 200 new points inside (the
#     30 exact copies), with a log line only.  A repeat is the same point of
#     predictor space: the reference is the nearest row with different
#     values.  Where folds split the repeats the threshold is 0 by
#     definition, and that is a warning with a class.
# ---------------------------------------------------------------------------

# 60 sites visited 4 times, and 200 new rows of which the first 30 are exact
# copies of a training row.  The layout the defect was reported on.
dist_twins <- function() {
  set.seed(1)
  site <- data.frame(a = stats::rnorm(60, 50, 10), b = stats::rnorm(60, 1000, 300),
                     c = stats::rnorm(60, 0.01, 0.002))
  tr <- site[rep(1:60, each = 4), ]
  nw <- rbind(site[1:30, ],
              data.frame(a = stats::rnorm(170, 50, 10),
                         b = stats::rnorm(170, 1000, 300),
                         c = stats::rnorm(170, 0.01, 0.002)))
  rownames(tr) <- NULL; rownames(nw) <- NULL
  list(tr = tr, nw = nw, site = rep(1:60, each = 4))
}
dist_pv <- c("a", "b", "c")


# --- (a) the dense path ------------------------------------------------------

test_that("the dense path gives an exact copy of a reference row distance 0", {
  set.seed(2)
  D <- scale(cbind(stats::rnorm(200, 50, 10), stats::rnorm(200, 1000, 300),
                   stats::rnorm(200, 0.01, 0.002), stats::rnorm(200)))
  attr(D, "scaled:center") <- NULL; attr(D, "scaled:scale") <- NULL
  Q <- D[sample(200, 120), , drop = FALSE]
  d <- .aoa_min_dist(Q, D, use_fnn = FALSE)
  # The fast formula alone left about one copy in ten near 1e-8.
  expect_identical(d, rep(0, 120))
  # Whatever the block size.
  expect_identical(.aoa_min_dist(Q, D, use_fnn = FALSE, chunk_size = 7L),
                   rep(0, 120))
  # And with repeated reference rows.
  expect_identical(.aoa_min_dist(Q, rbind(D, D[1:50, ]), use_fnn = FALSE),
                   rep(0, 120))
})

test_that("the dense path finds a copy among near-coincident reference rows", {
  set.seed(3)
  D <- matrix(stats::rnorm(400), 100, 4)
  # A second row within 1e-9 of row 1: the fast formula cannot tell them
  # apart, so the exact pass has to look at both.
  D2 <- rbind(D, D[1, ] + c(1e-9, 0, 0, 0))
  expect_identical(.aoa_min_dist(D[1, , drop = FALSE], D2, use_fnn = FALSE), 0)
  expect_identical(.aoa_min_dist(D2[101, , drop = FALSE], D2, use_fnn = FALSE), 0)
  # The leave-one-out distance between the two is their true separation.
  loo <- .aoa_min_dist(D2, D2, use_fnn = FALSE, exclude = seq_len(101))
  expect_equal(loo[c(1, 101)], c(1e-9, 1e-9), tolerance = 1e-5)
})

test_that("the dense and FNN paths agree to rounding without duplicates", {
  skip_if_not_installed("FNN")
  set.seed(4)
  D <- matrix(stats::rnorm(300 * 5), 300, 5)
  Q <- matrix(stats::rnorm(80 * 5), 80, 5)
  expect_equal(.aoa_min_dist(Q, D, use_fnn = FALSE),
               .aoa_min_dist(Q, D, use_fnn = TRUE), tolerance = 1e-12)
  expect_equal(.aoa_train_dist(D, use_fnn = FALSE),
               .aoa_train_dist(D, use_fnn = TRUE), tolerance = 1e-12)
  lab <- rep(1:5, length.out = 300)
  sp <- .aoa_fold_splits(lab, 300)
  expect_equal(.aoa_train_dist(D, splits = sp, use_fnn = FALSE),
               .aoa_train_dist(D, splits = sp, use_fnn = TRUE), tolerance = 1e-12)
})

test_that("with twins split across folds the two paths give the same zeros and flags", {
  skip_if_not_installed("FNN")
  L <- dist_twins()
  set.seed(9)
  lab <- sample(rep(1:5, length.out = 240))
  run <- function(fnn) suppressWarnings(
    area_of_applicability(L$nw, train_sf = L$tr, predictor_vars = dist_pv,
                          use_fnn = fnn, folds = lab))
  f <- run(TRUE); d <- run(FALSE)
  # The dense path used to report 216 zeros where FNN had 236, copies' DI up
  # to 2.65e-08, and 26 of the 30 copies inside.
  expect_identical(d$train_DI == 0, f$train_DI == 0)
  expect_identical(sum(d$train_DI == 0), 236L)
  expect_identical(d$threshold, 0)
  expect_identical(f$threshold, 0)
  expect_identical(d$aoa$DI[1:30], rep(0, 30))
  expect_identical(d$aoa$AOA, f$aoa$AOA)
  expect_identical(sum(d$aoa$AOA[1:30]), 30L)
  expect_equal(d$train_DI, f$train_DI, tolerance = 1e-12)
  expect_equal(d$aoa$DI, f$aoa$DI, tolerance = 1e-12)
})

test_that("a threshold of exactly 0 is not lost on the dense path", {
  # 100 sites visited twice; for 78 of them the two visits are in different
  # folds, so 156 of 200 training rows (78%) have their twin on the training
  # side and a DI of 0.  Q3 is then 0 and so is the threshold.  The dense
  # path used to lose about a tenth of those zeros to rounding, which took
  # the share below three quarters: the threshold came out near 1.6e-08 and
  # the zero-threshold caution never fired.
  set.seed(6)
  site <- data.frame(a = stats::rnorm(100, 50, 10), b = stats::rnorm(100, 1000, 300),
                     c = stats::rnorm(100, 0.01, 0.002), d = stats::rnorm(100))
  tr <- site[rep(1:100, each = 2), ]
  rownames(tr) <- NULL
  first  <- rep(1:5, length.out = 100)
  second <- ifelse(seq_len(100) <= 78, first %% 5 + 1, first)
  lab <- as.vector(rbind(first, second))
  nw <- data.frame(a = stats::rnorm(50, 50, 10), b = stats::rnorm(50, 1000, 300),
                   c = stats::rnorm(50, 0.01, 0.002), d = stats::rnorm(50))
  for (fnn in c(FALSE, if (requireNamespace("FNN", quietly = TRUE)) TRUE)) {
    expect_warning(
      res <- area_of_applicability(nw, train_sf = tr,
                                   predictor_vars = c("a", "b", "c", "d"),
                                   use_fnn = fnn, folds = lab),
      class = "spatialkit_aoa_zero_threshold")
    expect_identical(sum(res$train_DI == 0), 156L)
    expect_identical(res$threshold, 0)
  }
})


# --- (b) repeated rows without folds -----------------------------------------

test_that("repeat visits no longer set the threshold to 0 without folds", {
  L <- dist_twins()
  for (fnn in c(FALSE, if (requireNamespace("FNN", quietly = TRUE)) TRUE)) {
    # No warning of any kind: the case is handled, not degenerate.
    res <- expect_no_warning(
      area_of_applicability(L$nw, train_sf = L$tr, predictor_vars = dist_pv,
                            use_fnn = fnn))
    ded <- area_of_applicability(L$nw, train_sf = unique(L$tr),
                                 predictor_vars = dist_pv, use_fnn = fnn)
    # It was: threshold 0, 30 of 200 inside (26 on the dense path).
    expect_gt(res$threshold, 0)
    expect_gt(res$n_inside, 180L)
    expect_identical(res$n_inside, 188L)
    expect_identical(res$n_duplicates, 180L)
    expect_true(all(res$train_DI > 0))
    expect_length(res$train_DI, 240L)
    # Every visit to a site carries its site's reference distance.
    expect_identical(res$train_DI, res$train_DI[rep(seq(1, 240, by = 4), each = 4)])
    # The same decisions as the 60 distinct sites give, and the same
    # threshold as a distance.  The scales differ by known constants only:
    # the standard deviations are still taken over all 240 rows (n - 1 is
    # 239, not 59, so every scaled distance is sqrt(239 / 236) times the
    # distinct-site one), and the normaliser is still the mean pairwise
    # distance over all 240 rows, so that a prediction point's DI does not
    # depend on whether folds were passed.
    expect_identical(res$aoa$AOA, ded$aoa$AOA)
    expect_equal(res$threshold * res$normalizer * sqrt(236 / 239),
                 ded$threshold * ded$normalizer, tolerance = 1e-12)
    expect_identical(ded$n_duplicates, 0L)
  }
})

test_that("the distinct-row reference is logged and printed, not warned", {
  L <- dist_twins()
  lines <- capture_spatialkit_log(
    res <- area_of_applicability(L$nw, train_sf = L$tr, predictor_vars = dist_pv))
  expect_true(log_has(lines, "180 of 240 training rows repeat another row's predictor values exactly \\(60 distinct points\\)"))
  expect_false(log_has(lines, "threshold is 0"))
  out <- utils::capture.output(print(res))
  expect_true(any(grepl("180 of 240 rows repeat another's values", out, fixed = TRUE)))
})

test_that("a layer with no repeated row is computed on its rows as given", {
  skip_if_not_installed("FNN")
  set.seed(8)
  tr <- data.frame(a = stats::rnorm(150), b = stats::rnorm(150, 5, 2))
  nw <- data.frame(a = stats::rnorm(40), b = stats::rnorm(40, 5, 2))
  res <- area_of_applicability(nw, train_sf = tr, predictor_vars = c("a", "b"),
                               weights = NULL)
  expect_identical(res$n_duplicates, 0L)
  # By hand, the way it was always done: nearest other row over the mean
  # pairwise distance.
  X  <- as.matrix(tr)
  sc <- .aoa_scaling(X)
  Z  <- .aoa_apply_scaling(X, sc)
  nn <- as.numeric(FNN::get.knn(Z, k = 1L)$nn.dist[, 1L])
  expect_identical(res$train_DI, nn / mean(as.numeric(stats::dist(Z))))
  expect_identical(res$threshold, .aoa_threshold(res$train_DI))
  g <- .aoa_row_groups(Z)
  expect_identical(g$keep, seq_len(150))
  expect_identical(g$n_dup, 0L)
})

test_that(".aoa_row_groups groups exactly equal rows and nothing else", {
  X <- rbind(c(1, 2), c(3, 4), c(1, 2), c(1, 2 + 1e-12), c(3, 4), c(0, 0))
  g <- .aoa_row_groups(X)
  expect_identical(g$keep, c(1L, 2L, 4L, 6L))
  expect_identical(g$map, c(1L, 2L, 1L, 3L, 2L, 4L))
  expect_identical(g$n_dup, 2L)
  one <- .aoa_row_groups(matrix(5, 1, 3))
  expect_identical(one$n_dup, 0L)
})

test_that("folds that split the repeats warn with a class; folds that keep them do not", {
  L <- dist_twins()
  set.seed(9)
  lab <- sample(rep(1:5, length.out = 240))
  cond <- NULL
  res <- withCallingHandlers(
    area_of_applicability(L$nw, train_sf = L$tr, predictor_vars = dist_pv,
                          folds = lab),
    spatialkit_aoa_zero_threshold = function(w) {
      cond <<- w; invokeRestart("muffleWarning")
    })
  expect_s3_class(cond, "spatialkit_aoa_zero_threshold")
  expect_s3_class(cond, "warning")
  expect_match(conditionMessage(cond), "threshold is 0 because 236 of 240 training rows")
  expect_match(conditionMessage(cond), "leave_location_out")
  expect_identical(res$threshold, 0)
  expect_true(is.na(res$n_duplicates))
  # tryCatch() by class is what a script would write.
  caught <- tryCatch(
    area_of_applicability(L$nw, train_sf = L$tr, predictor_vars = dist_pv,
                          folds = lab),
    spatialkit_aoa_zero_threshold = function(w) "caught")
  expect_identical(caught, "caught")
  # One fold per site keeps every visit together: an ordinary threshold.
  by_site <- expect_no_warning(
    area_of_applicability(L$nw, train_sf = L$tr, predictor_vars = dist_pv,
                          folds = L$site %% 5))
  expect_gt(by_site$threshold, 0)
  # A threshold the caller supplied is theirs: nothing to warn about.
  expect_no_warning(
    area_of_applicability(L$nw, train_sf = L$tr, predictor_vars = dist_pv,
                          folds = lab, threshold = 0.4))
})
