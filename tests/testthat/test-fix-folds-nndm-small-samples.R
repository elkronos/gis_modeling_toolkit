# tests/testthat/test-fix-folds-nndm-small-samples.R
# ---------------------------------------------------------------------------
# The out-of-reach warning of make_folds(method = "nndm") on small samples.
#
# Its tolerance, 3 * sqrt(1/n + 1/m) for n folds and m prediction points, is
# 1 or more whenever n * m / (n + m) <= 9, and a deficit is a difference of
# two shares, so it cannot exceed 1.  For nine folds or fewer, and for up to
# eighteen with as many prediction points, plain leave-one-out therefore came
# back without a word however near the prediction points were.  There the
# warning now rests on complete separation: every prediction point nearer to
# the data than any sample is to its nearest other sample.  Each test that
# expects the warning fails on the code before the fix.
# ---------------------------------------------------------------------------

ns_pts <- function(xy) sf::st_as_sf(data.frame(x = 5e5 + xy[, 1], y = 5e6 + xy[, 2]),
                                    coords = c("x", "y"), crs = 32632)
ns_collect <- function(expr) {
  w <- list()
  val <- withCallingHandlers(expr, warning = function(c) {
    w[[length(w) + 1L]] <<- c; invokeRestart("muffleWarning")
  })
  list(value = val,
       reach = Filter(function(c) inherits(c, "spatialkit_nndm_out_of_reach"), w))
}

test_that("twelve samples with prediction points 5 m away raise the warning", {
  set.seed(3)
  xy <- cbind(runif(12, 0, 1e4), runif(12, 0, 1e4))
  u <- ns_pts(xy)
  near <- ns_pts(xy + matrix(rnorm(24, 0, 5), 12))
  res <- r2_quiet(ns_collect(make_folds(u, method = "nndm", prediction_points = near)))
  p <- res$value$params

  # Plain leave-one-out, two hundred times too far.
  expect_identical(p$n_removed_total, 0L)
  expect_gt(p$realised_median / p$target_median, 50)
  # The share rule cannot fire here: its tolerance is above any deficit.
  expect_equal(p$ecdf_deficit_tol, 3 * sqrt(1 / 12 + 1 / 12))
  expect_gt(p$ecdf_deficit_tol, 1)
  expect_identical(p$loo_ecdf_deficit, 1)
  # Every prediction distance is below every leave-one-out distance.
  expect_lt(max(p$target_distances), min(p$realised_distances))
  expect_true(p$target_separated)

  # It used to say nothing.
  expect_true(p$target_out_of_reach)
  expect_length(res$reach, 1L)
  cnd <- res$reach[[1L]]
  expect_true(cnd$separated)
  expect_equal(cnd$deficit, 1)
  expect_equal(cnd$tolerance, p$ecdf_deficit_tol)
  msg <- conditionMessage(cnd)
  expect_match(msg, "nearer to the data than leave-one-out can reach")
  expect_match(msg, "plain leave-one-out, no neighbour excluded")
  expect_match(msg, "Every one of the 12 prediction points is nearer to the data")
  expect_match(msg, "With 12 folds and 12 prediction points")
  expect_match(msg, "pessimistic for these prediction locations, an upper bound")
  # The numbers in the sentence are the ones in the result.
  expect_match(msg, sprintf("(%s at most)", format(signif(max(p$target_distances), 4))),
               fixed = TRUE)
  expect_match(msg, sprintf("(%s at least)", format(signif(min(p$realised_distances), 4))),
               fixed = TRUE)

  # fold_separation() reads the folds against the target they carry.
  s <- fold_separation(res$value, u)
  expect_true(attr(s, "target")$nndm_out_of_reach)
  expect_identical(attr(s, "target")$reading, "farther")
  out <- gsub("\\s+", " ", paste(utils::capture.output(print(s)), collapse = " "))
  expect_match(out, "NNDM can only lengthen a held-out distance")
})

test_that("a 3 x 3 sample predicted between its nodes is out of reach", {
  # Nine nodes 3333 apart; no point of the square is farther than 2357 from
  # one.  The tolerance for nine folds is above 1 whatever the target.
  g <- as.matrix(expand.grid((1:3 - 0.5) * 1e4 / 3, (1:3 - 0.5) * 1e4 / 3))
  set.seed(5)
  res <- r2_quiet(ns_collect(make_folds(
    ns_pts(g), method = "nndm",
    prediction_points = ns_pts(cbind(runif(1500, 0, 1e4), runif(1500, 0, 1e4))))))
  p <- res$value$params
  expect_gt(p$ecdf_deficit_tol, 1)
  expect_true(p$target_separated)
  expect_true(p$target_out_of_reach)
  expect_length(res$reach, 1L)
  expect_match(conditionMessage(res$reach[[1L]]),
               "Every one of the 1500 prediction points is nearer to the data")
})

test_that("twelve samples with a grid over their extent raise nothing", {
  for (seed in c(201, 202, 203, 204)) {
    set.seed(seed)
    u <- ns_pts(cbind(runif(12, 0, 1e4), runif(12, 0, 1e4)))
    grid <- sf::st_as_sf(sf::st_make_grid(u, n = c(20, 20), what = "centers"))
    res <- r2_quiet(ns_collect(make_folds(u, method = "nndm", prediction_points = grid)))
    p <- res$value$params
    expect_length(res$reach, 0L)
    expect_false(p$target_separated)
    expect_false(p$target_out_of_reach)
    expect_lt(p$loo_ecdf_deficit, 1)
  }
})

test_that("a partial shortfall on a small sample is recorded, not warned about", {
  # Prediction points at 0.3 of the mean spacing from twelve samples: most
  # are nearer than leave-one-out reaches, not all.  With 12 against 12 that
  # is within what two small sets of distances differ by; the deficit and
  # the two medians are in the result for the caller to read.
  set.seed(9048)
  xy <- cbind(runif(12, 0, 1e4), runif(12, 0, 1e4))
  a <- runif(12, 0, 2 * pi); sp <- 1e4 / sqrt(12)
  res <- r2_quiet(ns_collect(make_folds(
    ns_pts(xy), method = "nndm",
    prediction_points = ns_pts(xy + 0.3 * sp * cbind(cos(a), sin(a))))))
  p <- res$value$params
  expect_false(p$target_separated)
  expect_gt(p$loo_ecdf_deficit, 0.5)
  expect_lt(p$loo_ecdf_deficit, 1)
  expect_gt(p$ecdf_deficit_tol, 1)
  expect_false(p$target_out_of_reach)
  expect_length(res$reach, 0L)
  expect_gt(p$realised_median, p$target_median)
})

test_that("a handful of distances is too few to call a separation", {
  # Four samples and two prediction points beside two of them: the ordering
  # is complete, and one of only choose(6, 4) = 15, so it is recorded and
  # nothing is raised.
  xy <- cbind(c(1000, 9000, 1000, 9000), c(1000, 1000, 9000, 9500))
  u <- ns_pts(xy)
  two <- ns_pts(xy[1:2, ] + 5)
  res <- r2_quiet(ns_collect(make_folds(u, method = "nndm", prediction_points = two)))
  p <- res$value$params
  expect_true(p$target_separated)
  expect_false(p$target_out_of_reach)
  expect_length(res$reach, 0L)

  # Seventeen prediction points beside three samples are enough:
  # choose(20, 3) = 1140.
  tri <- cbind(c(1000, 9000, 5000), c(1000, 1500, 9000))
  set.seed(4)
  near <- ns_pts(tri[rep(1:3, length.out = 17), ] + matrix(rnorm(34, 0, 5), 17))
  res <- r2_quiet(ns_collect(make_folds(ns_pts(tri), method = "nndm", prediction_points = near)))
  expect_true(res$value$params$target_separated)
  expect_true(res$value$params$target_out_of_reach)
  expect_length(res$reach, 1L)
  # Sixteen are one short: choose(19, 3) = 969.
  res <- r2_quiet(ns_collect(make_folds(ns_pts(tri), method = "nndm",
                                        prediction_points = near[1:16, ])))
  expect_true(res$value$params$target_separated)
  expect_false(res$value$params$target_out_of_reach)
  expect_length(res$reach, 0L)
})

test_that("on larger samples the share rule decides as before", {
  # 200 against 200, 5 m away: past the tolerance, and separated too.
  set.seed(1)
  xy <- cbind(runif(200, 0, 1e4), runif(200, 0, 1e4))
  res <- r2_quiet(ns_collect(make_folds(
    ns_pts(xy), method = "nndm",
    prediction_points = ns_pts(xy + matrix(rnorm(400, 0, 5), 200)))))
  cnd <- res$reach[[1L]]
  expect_lt(cnd$tolerance, 1)
  expect_gt(cnd$deficit, cnd$tolerance)
  # The message is the share rule's, with its tolerance.
  expect_match(conditionMessage(cnd), "past the 0.300 allowed for 200 folds and 200 prediction points")
  expect_no_match(conditionMessage(cnd), "Every one of the")
  # A uniform sample with a grid target: neither rule fires.
  set.seed(101)
  u <- ns_pts(cbind(runif(300, 0, 1e4), runif(300, 0, 1e4)))
  grid <- sf::st_as_sf(sf::st_make_grid(u, n = c(40, 40), what = "centers"))
  res <- r2_quiet(ns_collect(make_folds(u, method = "nndm", prediction_points = grid)))
  expect_length(res$reach, 0L)
  expect_false(res$value$params$target_separated)
  expect_false(res$value$params$target_out_of_reach)
})
