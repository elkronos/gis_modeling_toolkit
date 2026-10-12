# tests/testthat/test-fix-folds-nndm-reach.R
# ---------------------------------------------------------------------------
# make_folds(method = "nndm") returned plain leave-one-out, with no warning
# and max_ecdf_excess = 0, when the prediction points sat nearer to the data
# than leave-one-out can reach.  Excluding neighbours only lengthens a
# held-out distance, so that side of the match is out of reach by
# construction; it is now measured, recorded and warned about.  Each test
# fails on the code before the fix.
# ---------------------------------------------------------------------------

nr_pts <- function(xy) sf::st_as_sf(data.frame(x = 5e5 + xy[, 1], y = 5e6 + xy[, 2]),
                                    coords = c("x", "y"), crs = 32632)
nr_grid <- function(pts, m = 40) {
  bb <- sf::st_bbox(pts)
  sf::st_as_sf(expand.grid(x = seq(bb[["xmin"]], bb[["xmax"]], length.out = m),
                           y = seq(bb[["ymin"]], bb[["ymax"]], length.out = m)),
               coords = c("x", "y"), crs = sf::st_crs(pts))
}
nr_collect <- function(expr) {
  w <- list()
  val <- withCallingHandlers(expr, warning = function(c) {
    w[[length(w) + 1L]] <<- c; invokeRestart("muffleWarning")
  })
  list(value = val,
       reach = Filter(function(c) inherits(c, "spatialkit_nndm_out_of_reach"), w),
       all = vapply(w, conditionMessage, ""))
}

test_that("prediction points 5 m from the samples raise the classed warning", {
  set.seed(1)
  xy <- cbind(runif(200, 0, 1e4), runif(200, 0, 1e4))
  u <- nr_pts(xy)
  near <- nr_pts(xy + matrix(rnorm(400, 0, 5), 200))
  res <- r2_quiet(nr_collect(make_folds(u, method = "nndm", prediction_points = near)))
  p <- res$value$params

  # What came back is plain leave-one-out, a hundred times too far.
  expect_identical(p$n_removed_total, 0L)
  expect_true(all(lengths(lapply(res$value$folds, `[[`, "train")) == 199L))
  expect_gt(p$realised_median / p$target_median, 20)
  expect_identical(p$max_ecdf_excess, 0)

  # It used to say nothing.  Now: one warning of its own class ...
  expect_length(res$reach, 1L)
  cnd <- res$reach[[1L]]
  expect_match(conditionMessage(cnd), "nearer to the data than leave-one-out can reach")
  expect_match(conditionMessage(cnd), "plain leave-one-out, no neighbour excluded")
  expect_match(conditionMessage(cnd), "pessimistic for these prediction locations, an upper bound")
  expect_equal(cnd$deficit, p$loo_ecdf_deficit)
  expect_equal(cnd$tolerance, p$ecdf_deficit_tol)
  # ... and a record in the result.
  expect_true(p$target_out_of_reach)
  expect_gt(p$loo_ecdf_deficit, 0.9)
  # Nothing was excluded, so the folds' deficit is leave-one-out's.
  expect_equal(p$max_ecdf_deficit, p$loo_ecdf_deficit)
  expect_equal(p$ecdf_deficit_tol, 3 * sqrt(1 / 200 + 1 / 200))
  # Nothing was excluded, so the folds are as far from the target as plain
  # leave-one-out is.
  expect_equal(p$wasserstein, p$wasserstein_loo)
  expect_gt(p$wasserstein, 10 * p$target_median)

  # The deficit is the mirror of max_ecdf_excess, recomputed here by hand.
  at <- sort(unique(c(p$target_distances, p$realised_distances)))
  by_hand <- max(stats::ecdf(p$target_distances)(at) - stats::ecdf(p$realised_distances)(at))
  expect_equal(p$max_ecdf_deficit, by_hand)
})

test_that("a regular sample predicted between its nodes is out of reach too", {
  g <- as.matrix(expand.grid(seq(250, 9750, by = 500), seq(250, 9750, by = 500)))
  set.seed(5)
  res <- r2_quiet(nr_collect(make_folds(
    nr_pts(g), method = "nndm",
    prediction_points = nr_pts(cbind(runif(1500, 0, 1e4), runif(1500, 0, 1e4))))))
  expect_length(res$reach, 1L)
  expect_true(res$value$params$target_out_of_reach)
})

test_that("layouts whose folds match their target raise no such warning", {
  # A uniform sample with a grid over its extent: plain leave-one-out is the
  # right answer, and one fold's worth of deficit (1/n) would fire on every
  # such layout.
  for (seed in c(101, 102, 103)) {
    set.seed(seed)
    u <- nr_pts(cbind(runif(300, 0, 1e4), runif(300, 0, 1e4)))
    res <- r2_quiet(nr_collect(make_folds(u, method = "nndm", prediction_points = nr_grid(u))))
    p <- res$value$params
    expect_length(res$reach, 0L)
    expect_false(p$target_out_of_reach)
    expect_gt(p$loo_ecdf_deficit, 1 / 300)          # the 1/n rule would have fired
    expect_lt(p$loo_ecdf_deficit, p$ecdf_deficit_tol)
  }
  # A clustered sample with the grid, the layout NNDM is for: neighbours are
  # excluded, and the folds end closer to the target than leave-one-out was.
  set.seed(2001)
  c0 <- cbind(runif(8, 1500, 8500), runif(8, 1500, 8500)); i <- sample(8, 200, TRUE)
  cl <- nr_pts(c0[i, ] + matrix(rnorm(400, 0, 400), 200))
  res <- r2_quiet(nr_collect(make_folds(cl, method = "nndm", prediction_points = nr_grid(cl))))
  p <- res$value$params
  expect_length(res$reach, 0L)
  expect_false(p$target_out_of_reach)
  expect_gt(p$n_removed_total, 0L)
  expect_lt(p$wasserstein, p$wasserstein_loo)
})

test_that("an overshoot of the matching is recorded, and is not called out of reach", {
  # Two tight clusters and a grid over the square between them.  Plain
  # leave-one-out is far too close here, so nothing is out of reach; but
  # excluding a point's own cluster moves it out to the other one, past most
  # of the prediction distances.  That is the granularity of the matching,
  # a different thing from prediction points nearer than leave-one-out.
  set.seed(13)
  ctr <- cbind(c(250, 750), c(250, 750)); i <- rep(1:2, each = 40)
  two <- nr_pts(cbind(ctr[i, 1] + rnorm(80, 0, 40), ctr[i, 2] + rnorm(80, 0, 40)))
  grid <- sf::st_as_sf(expand.grid(x = 5e5 + seq(50, 950, by = 75), y = 5e6 + seq(50, 950, by = 75)),
                       coords = c("x", "y"), crs = 32632)
  res <- r2_quiet(nr_collect(make_folds(two, method = "nndm", prediction_points = grid)))
  p <- res$value$params
  expect_length(res$reach, 0L)
  expect_false(p$target_out_of_reach)
  expect_lt(p$loo_ecdf_deficit, p$ecdf_deficit_tol)
  expect_gt(p$max_ecdf_deficit, p$loo_ecdf_deficit)
  expect_gt(p$n_removed_total, 0L)
  # fold_separation() still reads the folds as they are, without blaming
  # the prediction points.
  s <- fold_separation(res$value, two)
  expect_false(attr(s, "target")$nndm_out_of_reach)
  out <- gsub("\\s+", " ", paste(utils::capture.output(print(s)), collapse = " "))
  expect_no_match(out, "NNDM can only lengthen")
})

test_that("fold_separation() reads NNDM folds against the target they carry", {
  set.seed(1)
  xy <- cbind(runif(200, 0, 1e4), runif(200, 0, 1e4))
  u <- nr_pts(xy)
  far <- r2_quiet(suppressWarnings(make_folds(
    u, method = "nndm", prediction_points = nr_pts(xy + matrix(rnorm(400, 0, 5), 200)))))
  s <- fold_separation(far, u, sac = 2000)
  tg <- attr(s, "target")
  expect_identical(tg$reading, "farther")
  expect_identical(tg$basis, "all sites")
  expect_true(tg$nndm_out_of_reach)
  expect_equal(tg$target_all_median, far$params$target_median)
  out <- gsub("\\s+", " ", paste(utils::capture.output(print(s)), collapse = " "))
  # It used to say these folds "reproduce the distances from the prediction
  # points to the data", for every NNDM fold set.
  expect_no_match(out, "reproduce the distances")
  expect_match(out, "farther from their training data than the prediction points do")
  expect_match(out, "NNDM can only lengthen a held-out distance")

  set.seed(101)
  v <- nr_pts(cbind(runif(300, 0, 1e4), runif(300, 0, 1e4)))
  ok <- r2_quiet(make_folds(v, method = "nndm", prediction_points = nr_grid(v)))
  expect_identical(attr(fold_separation(ok, v, sac = 2000), "target")$reading, "matched")
  expect_identical(attr(fold_separation(ok, v), "target")$reading, "matched")
})
