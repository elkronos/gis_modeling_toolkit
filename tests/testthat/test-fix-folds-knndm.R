# tests/testthat/test-fix-folds-knndm.R
# ---------------------------------------------------------------------------
# make_folds(method = "knndm"): k folds matched to the distances from the
# prediction points to the data (Linnenbrink et al. 2024).  Distance-matched
# folds existed only as leave-one-out, refused above 5,000 rows.  The
# candidates are a random partition and clusterings of the samples packed
# into k folds; the one whose held-out distances are closest to the target
# (smallest Wasserstein distance W) is returned.
# ---------------------------------------------------------------------------

kn_pts <- function(xy) sf::st_as_sf(data.frame(x = 5e5 + xy[, 1], y = 5e6 + xy[, 2]),
                                    coords = c("x", "y"), crs = 32632)
kn_box_grid <- function(lo = 0, hi = 1e4, m = 30)
  kn_pts(as.matrix(expand.grid(seq(lo, hi, length.out = m), seq(lo, hi, length.out = m))))
kn_uniform <- function(n, seed) { set.seed(seed); cbind(runif(n, 0, 1e4), runif(n, 0, 1e4)) }
kn_clustered <- function(n, seed, nc = 12) {
  set.seed(seed)
  c0 <- cbind(runif(nc, 1000, 9000), runif(nc, 1000, 9000)); i <- sample(nc, n, TRUE)
  pmin(pmax(c0[i, ] + matrix(rnorm(2 * n, 0, 600), n), 0), 1e4)
}
kn_collect <- function(expr) {
  w <- list()
  val <- withCallingHandlers(expr, warning = function(c) {
    w[[length(w) + 1L]] <<- c; invokeRestart("muffleWarning")
  })
  list(value = val, warnings = w,
       classes = vapply(w, function(c) class(c)[1L], ""))
}
kn_ns <- function(name) get(name, envir = asNamespace("spatialkit"))

test_that("knndm returns k folds that partition the rows", {
  pts <- kn_pts(kn_clustered(300, 1))
  f <- r2_quiet(make_folds(pts, k = 5, method = "knndm",
                           prediction_points = kn_box_grid(), seed = 1))
  expect_identical(f$method, "knndm")
  expect_identical(f$k, 5L)
  expect_length(f$folds, 5L)
  tests <- unlist(lapply(f$folds, `[[`, "test"))
  expect_setequal(tests, seq_len(300))
  expect_false(anyDuplicated(tests) > 0)
  for (s in f$folds) {
    expect_length(intersect(s$train, s$test), 0L)
    expect_setequal(c(s$train, s$test), seq_len(300))
  }
  # `assignment` says the same thing as the splits.
  expect_identical(f$assignment$row_id, seq_len(300))
  for (j in 1:5)
    expect_setequal(f$assignment$row_id[f$assignment$fold == j], f$folds[[j]]$test)
  p <- f$params
  expect_identical(p$k_requested, 5L)
  expect_identical(p$fold_sizes, tabulate(f$assignment$fold, 5L))
  expect_length(p$realised_distances, 300L)
  expect_length(p$target_distances, 900L)
  expect_identical(p$crs, "EPSG:32632")
  expect_identical(p$clustering, "ward")
})

test_that("the realised distances are each point's distance to the nearest point of another fold", {
  xy <- kn_clustered(200, 2)
  f <- r2_quiet(make_folds(kn_pts(xy), k = 4, method = "knndm",
                           prediction_points = kn_box_grid(), seed = 2))
  lab <- f$assignment$fold
  D <- as.matrix(stats::dist(xy))
  by_hand <- vapply(seq_len(200), function(i) min(D[i, lab != lab[i]]), 0)
  expect_equal(f$params$realised_distances, by_hand)
  expect_equal(f$params$realised_median, stats::median(by_hand))
  # The target is the distance from each prediction point to its nearest sample.
  g <- sf::st_coordinates(kn_box_grid()) - cbind(rep(5e5, 900), rep(5e6, 900))
  tg <- vapply(seq_len(900), function(i) sqrt(min((xy[, 1] - g[i, 1])^2 + (xy[, 2] - g[i, 2])^2)), 0)
  expect_equal(f$params$target_distances, tg)
})

test_that("W of the folds returned is never above W of its own random candidate", {
  for (seed in 1:4) for (layout in c("uniform", "clustered")) {
    xy <- if (layout == "uniform") kn_uniform(250, seed) else kn_clustered(250, seed)
    f <- r2_quiet(suppressWarnings(make_folds(
      kn_pts(xy), k = 5, method = "knndm", prediction_points = kn_box_grid(),
      seed = seed)))
    p <- f$params; info <- paste(layout, seed)
    expect_lte(p$wasserstein, p$wasserstein_random, label = info)
    # The random partition is the first candidate, and is always admissible.
    expect_identical(p$candidates$candidate[1L], "random", info = info)
    expect_true(p$candidates$admissible[1L], info = info)
    expect_equal(p$candidates$W[1L], p$wasserstein_random, info = info)
    # What is returned is the smallest W among the admissible candidates.
    expect_equal(p$wasserstein, min(p$candidates$W[p$candidates$admissible]), info = info)
    # W is the Wasserstein distance between the two sets of distances.
    expect_equal(p$wasserstein,
                 kn_ns(".fold_sep_compare")(p$realised_distances, p$target_distances)$w1,
                 info = info)
    # A candidate that was set aside was not scored.
    expect_true(all(is.na(p$candidates$W[!p$candidates$admissible])), info = info)
  }
})

test_that("clustering brings the folds closer to the target than a random partition on clustered samples", {
  # 12 clusters predicted onto a grid over the box: a random partition keeps
  # a neighbour of every held-out point in training.
  ratio <- vapply(1:4, function(seed) {
    p <- r2_quiet(make_folds(kn_pts(kn_clustered(300, seed)), k = 5, method = "knndm",
                             prediction_points = kn_box_grid(), seed = seed))$params
    expect_identical(p$chosen, "clusters")
    expect_true(p$n_clusters >= 5L)
    p$wasserstein / p$wasserstein_random
  }, 0)
  expect_lt(max(ratio), 0.9)
})

test_that("the random candidate is the partition random_kfold builds from the same seed", {
  # Prediction points 5 m from the samples: no clustering comes closer than
  # the random partition, which is what comes back.
  xy <- kn_uniform(150, 3)
  set.seed(30); near <- kn_pts(xy + matrix(rnorm(300, 0, 5), 150))
  f <- r2_quiet(suppressWarnings(make_folds(kn_pts(xy), k = 5, method = "knndm",
                                            prediction_points = near, seed = 11)))
  expect_identical(f$params$chosen, "random")
  expect_true(is.na(f$params$n_clusters))
  r <- make_folds(kn_pts(xy), k = 5, method = "random_kfold", seed = 11)
  expect_identical(f$assignment$fold, r$assignment$fold)
  expect_equal(f$params$wasserstein, f$params$wasserstein_random)
  # And for a fold count that does not divide n.
  set.seed(5); lab <- kn_ns(".random_fold_labels")(23L, 4L)
  r2 <- make_folds(kn_pts(xy[1:23, ]), k = 4, method = "random_kfold", seed = 5)
  expect_identical(lab, r2$assignment$fold)
})

test_that("fold sizes stay within balance_tol, which is a limit for knndm", {
  # Four clusters of very different sizes: few clusters per fold is uneven.
  set.seed(4)
  c0 <- cbind(c(2000, 8000, 2000, 8000), c(2000, 2000, 8000, 8000))
  i <- sample(4, 240, TRUE, prob = c(0.55, 0.25, 0.15, 0.05))
  pts <- kn_pts(c0[i, ] + matrix(rnorm(480, 0, 400), 240))
  far <- kn_box_grid(-1e4, 2e4)
  loose <- r2_quiet(suppressWarnings(make_folds(pts, k = 3, method = "knndm",
                                                prediction_points = far, seed = 1)))
  expect_lte(loose$params$balance_ratio, 3)
  expect_identical(loose$params$balance_tol, 3)
  tight <- r2_quiet(suppressWarnings(make_folds(pts, k = 3, method = "knndm",
                                                prediction_points = far, seed = 1,
                                                balance_tol = 1.1)))
  expect_lte(tight$params$balance_ratio, 1.1)
  expect_equal(tight$params$balance_ratio,
               max(tight$params$fold_sizes) / min(tight$params$fold_sizes))
  # The limit set candidates aside, and those are the ones past it.
  cd <- tight$params$candidates
  cl <- cd[cd$candidate == "clusters" & !is.na(cd$balance_ratio), ]
  expect_true(any(!cl$admissible))
  expect_true(all(cl$balance_ratio[cl$admissible] <= 1.1))
  # With the limit lifted the uneven candidates are scored, and the match is
  # at least as close.
  free <- r2_quiet(suppressWarnings(make_folds(pts, k = 3, method = "knndm",
                                               prediction_points = far, seed = 1,
                                               balance_tol = Inf)))
  expect_lte(free$params$wasserstein, tight$params$wasserstein)
  # On the default tolerance over several clustered draws.
  for (seed in 1:5) {
    p <- r2_quiet(suppressWarnings(make_folds(
      kn_pts(kn_clustered(200, seed, nc = 5)), k = 5, method = "knndm",
      prediction_points = kn_box_grid(), seed = seed)))$params
    expect_lte(p$balance_ratio, 3, label = paste("seed", seed))
  }
})

test_that("min_train caps the largest fold, and a value no k folds can meet is an error", {
  pts <- kn_pts(kn_clustered(200, 6, nc = 3))
  far <- kn_box_grid(-1e4, 2e4)
  # Two folds at the default 0.5: each fold is exactly half.
  two <- r2_quiet(suppressWarnings(make_folds(pts, k = 2, method = "knndm",
                                              prediction_points = far, seed = 1)))
  expect_identical(sort(two$params$fold_sizes), c(100L, 100L))
  # A lower floor admits uneven halves.
  low <- r2_quiet(suppressWarnings(make_folds(pts, k = 2, method = "knndm",
                                              prediction_points = far, seed = 1,
                                              min_train = 0.3)))
  expect_lte(max(low$params$fold_sizes), 140L)
  expect_gt(sum(low$params$candidates$admissible), sum(two$params$candidates$admissible))
  expect_lte(low$params$wasserstein, two$params$wasserstein)
  expect_error(make_folds(pts, k = 5, method = "knndm", prediction_points = far,
                          min_train = 0.9),
               "min_train = 0.9 cannot be met with k = 5 folds.*min_train <= 0.8")
  expect_error(make_folds(pts, k = 5, method = "knndm", prediction_points = far,
                          min_train = 1),
               "`min_train` must be a single number in \\(0, 1\\)")
})

test_that("a seed makes knndm reproducible, and the row order of the result is the input's", {
  pts <- kn_pts(kn_clustered(250, 7))
  a <- r2_quiet(make_folds(pts, k = 5, method = "knndm", prediction_points = kn_box_grid(), seed = 9))
  b <- r2_quiet(make_folds(pts, k = 5, method = "knndm", prediction_points = kn_box_grid(), seed = 9))
  expect_identical(a$assignment, b$assignment)
  expect_identical(a$params$candidates, b$params$candidates)
  # The session RNG is left where it was.
  set.seed(1); r1 <- runif(1)
  set.seed(1); invisible(r2_quiet(make_folds(pts, k = 5, method = "knndm",
                                             prediction_points = kn_box_grid(), seed = 9)))
  expect_identical(runif(1), r1)
  # Row IDs that are not 1..n travel through the splits.
  pts$..row_id <- 1000L + seq_len(250)
  d <- r2_quiet(make_folds(pts, k = 5, method = "knndm", prediction_points = kn_box_grid(), seed = 9))
  expect_setequal(unlist(lapply(d$folds, `[[`, "test")), 1000L + seq_len(250))
  expect_identical(d$assignment$fold, a$assignment$fold)
})

test_that("prediction points nearer than any fold can reach raise the out-of-reach warning", {
  xy <- kn_uniform(200, 8)
  set.seed(80); near <- kn_pts(xy + matrix(rnorm(400, 0, 5), 200))
  res <- r2_quiet(kn_collect(make_folds(kn_pts(xy), k = 5, method = "knndm",
                                        prediction_points = near, seed = 1)))
  expect_identical(res$classes, "spatialkit_nndm_out_of_reach")
  cnd <- res$warnings[[1L]]; p <- res$value$params
  expect_match(conditionMessage(cnd), "nearer to the data than any fold can hold a point out")
  # What comes back is the random partition or a candidate as fine as one:
  # nothing coarser can be nearer.
  expect_match(conditionMessage(cnd), "the 5 folds returned \\((the random partition|19[0-9] clusters)")
  expect_true(p$chosen == "random" || p$n_clusters >= 190L)
  expect_match(conditionMessage(cnd), "pessimistic for these prediction locations, an upper bound")
  expect_true(p$target_out_of_reach)
  expect_false(p$short_of_target)
  expect_equal(cnd$deficit, p$loo_ecdf_deficit)
  expect_equal(cnd$tolerance, p$ecdf_deficit_tol)
  expect_equal(p$ecdf_deficit_tol, 3 * sqrt(1 / 200 + 1 / 200))
  # The same rule, on the same two sets of distances, as method = "nndm".
  nn <- r2_quiet(suppressWarnings(make_folds(kn_pts(xy), method = "nndm",
                                             prediction_points = near)))$params
  expect_equal(p$loo_ecdf_deficit, nn$loo_ecdf_deficit)
  expect_equal(p$wasserstein_loo, nn$wasserstein_loo)
  expect_identical(p$target_separated, nn$target_separated)
  # A k-fold partition trains on fewer points than leave-one-out, so it is
  # no nearer to this target.
  expect_gte(p$wasserstein, p$wasserstein_loo)
})

test_that("a grid over evenly spread samples raises no warning", {
  for (seed in 1:5) {
    pts <- kn_pts(kn_uniform(300, seed))
    bb <- sf::st_bbox(pts)
    grid <- sf::st_as_sf(expand.grid(x = seq(bb[["xmin"]], bb[["xmax"]], length.out = 30),
                                     y = seq(bb[["ymin"]], bb[["ymax"]], length.out = 30)),
                         coords = c("x", "y"), crs = 32632)
    res <- r2_quiet(kn_collect(make_folds(pts, k = 5, method = "knndm",
                                          prediction_points = grid, seed = seed)))
    expect_length(res$warnings, 0L)
    expect_false(res$value$params$target_out_of_reach)
    expect_false(res$value$params$short_of_target)
    # The folds sit within a third of the target mean of the target (0.05
    # to 0.22 over these draws; a random partition on clustered samples is
    # at 0.7).
    expect_lt(res$value$params$wasserstein / res$value$params$target_mean, 1 / 3)
  }
})

test_that("prediction points well beyond the samples raise the short-of-target warning", {
  pts <- kn_pts(kn_clustered(300, 9))
  res <- r2_quiet(kn_collect(make_folds(pts, k = 5, method = "knndm",
                                        prediction_points = kn_box_grid(-1e4, 2e4),
                                        seed = 1)))
  expect_identical(res$classes, "spatialkit_knndm_short_of_target")
  cnd <- res$warnings[[1L]]; p <- res$value$params
  expect_true(p$short_of_target)
  expect_false(p$target_out_of_reach)
  expect_match(conditionMessage(cnd), "still short of the target")
  expect_match(conditionMessage(cnd), "tends to read low for these prediction locations")
  expect_match(conditionMessage(cnd), "upper bound on performance")
  # Five folds of twelve clusters: the candidate kept is the coarsest the
  # two limits allow, and the message says which limits.
  expect_match(conditionMessage(cnd), "the coarsest clustering allowed")
  expect_match(conditionMessage(cnd), "min_train = 0.5")
  expect_match(conditionMessage(cnd), "Fewer folds hold out larger regions")
  expect_equal(cnd$held_out_mean, p$realised_mean)
  expect_equal(cnd$target_mean, p$target_mean)
  expect_equal(p$realised_mean, mean(p$realised_distances))
  expect_equal(p$target_mean, mean(p$target_distances))
  # The rule: the mean held-out distance is more than a quarter below the
  # target's, and that shift is at least half of W.
  expect_lt(p$realised_mean, 0.75 * p$target_mean)
  expect_gte(p$target_mean - p$realised_mean, 0.5 * p$wasserstein)
  # It is still the closest candidate.
  expect_lte(p$wasserstein, p$wasserstein_random)
  # A wider area is matched by coarser clusters than the samples' own extent.
  own <- r2_quiet(make_folds(pts, k = 5, method = "knndm",
                             prediction_points = kn_box_grid(), seed = 1))$params
  expect_false(own$short_of_target)
  expect_lt(p$n_clusters, own$n_clusters)
  expect_gt(p$realised_median, own$realised_median)
})

test_that("k is required, raised from 1 and lowered to n with the classed warning", {
  pts <- kn_pts(kn_uniform(8, 10))
  grid <- kn_box_grid(m = 6)
  expect_error(make_folds(pts, method = "knndm", prediction_points = grid),
               "`k` \\(the number of folds\\) is required for method = \"knndm\"")
  one <- r2_quiet(suppressWarnings(make_folds(pts, k = 1, method = "knndm",
                                              prediction_points = grid, seed = 1)))
  expect_identical(one$k, 2L)
  res <- r2_quiet(kn_collect(make_folds(pts, k = 20, method = "knndm",
                                        prediction_points = grid, seed = 1)))
  red <- Filter(function(c) inherits(c, "spatialkit_k_reduced"), res$warnings)
  expect_length(red, 1L)
  expect_identical(red[[1L]]$k_requested, 20L)
  expect_identical(red[[1L]]$k, 8L)
  expect_match(conditionMessage(red[[1L]]), "returning 8 folds, not the k = 20 asked for")
  expect_identical(res$value$k, 8L)
  expect_identical(res$value$params$k_requested, 20L)
  # One point per fold is leave-one-out: the held-out distance is the
  # distance to the nearest other sample.
  D <- as.matrix(stats::dist(sf::st_coordinates(pts))); diag(D) <- Inf
  expect_equal(res$value$params$realised_distances, unname(apply(D, 1, min)))
  expect_error(make_folds(kn_pts(kn_uniform(2, 1)), k = 2, method = "knndm",
                          prediction_points = grid),
               "need at least 3 points")
})

test_that("prediction_points are required, and cells are reduced to points", {
  pts <- kn_pts(kn_clustered(120, 11))
  expect_error(make_folds(pts, k = 4, method = "knndm"),
               "make_folds\\(knndm\\): `prediction_points` is required")
  cells <- sf::st_make_grid(sf::st_as_sfc(sf::st_bbox(pts)), n = c(12, 12))
  centres <- sf::st_centroid(cells)
  a <- r2_quiet(suppressWarnings(make_folds(pts, k = 4, method = "knndm",
                                            prediction_points = cells, seed = 2)))
  b <- r2_quiet(suppressWarnings(make_folds(pts, k = 4, method = "knndm",
                                            prediction_points = centres, seed = 2)))
  expect_equal(a$params$target_distances, b$params$target_distances)
  expect_identical(a$assignment$fold, b$assignment$fold)
  expect_identical(a$params$n_prediction_points, 144L)
  # The arguments of the other methods are refused or unused, not half-read.
  expect_error(make_folds(pts, k = 4, method = "knndm", prediction_points = centres,
                          max_test = 10),
               "`max_test` is only used by the leave-one-out methods")
  expect_error(make_folds(pts, k = 4, method = "knndm", prediction_points = centres,
                          blocks = cells),
               "`blocks` is only used by method = \"block_kfold\"")
})

test_that("above the Ward limit the samples are grouped by k-means and Ward carries on", {
  core <- kn_ns(".knndm_folds")
  xy <- kn_clustered(600, 12)
  g <- sf::st_coordinates(kn_box_grid()) - cbind(rep(5e5, 900), rep(5e6, 900))
  D <- sqrt(outer(g[, 1], xy[, 1], "-")^2 + outer(g[, 2], xy[, 2], "-")^2)
  tg <- apply(D, 1, min)
  # Without FNN: the search falls back to sf.  No skip, so this path runs on
  # every machine.
  set.seed(3); small <- core(xy, tg, k = 5, ward_max = 150L, use_fnn = FALSE)
  expect_identical(small$clustering, "kmeans+ward")
  expect_setequal(unique(small$fold), 1:5)
  cd <- small$candidates
  expect_identical(cd$candidate[1L], "random")
  expect_equal(cd$W[small$chosen], min(cd$W[cd$admissible]))
  expect_lte(cd$W[small$chosen], small$w_random)
  # Cluster counts above the 150 groups were tried too (k-means runs of
  # their own), and none has more clusters than distinct points.
  expect_true(any(cd$n_clusters > 150L, na.rm = TRUE))
  expect_true(all(cd$n_clusters < 600L, na.rm = TRUE))
  expect_lt(cd$W[small$chosen], small$w_random)
  # The held-out distances returned are those of the folds returned.
  Dx <- as.matrix(stats::dist(xy))
  by_hand <- vapply(seq_len(600), function(i) min(Dx[i, small$fold != small$fold[i]]), 0)
  expect_equal(small$realised, by_hand)
  # Ward on all points for comparison: the same kind of answer.
  set.seed(3); full <- core(xy, tg, k = 5, use_fnn = FALSE)
  expect_identical(full$clustering, "ward")
  expect_lt(abs(cd$W[small$chosen] - full$candidates$W[full$chosen]) / mean(tg), 0.15)
})

test_that("the FNN and sf searches give the same folds", {
  skip_if_not_installed("FNN")
  core <- kn_ns(".knndm_folds")
  xy <- kn_clustered(500, 13)
  g <- sf::st_coordinates(kn_box_grid()) - cbind(rep(5e5, 900), rep(5e6, 900))
  tg <- FNN::get.knnx(xy, g, k = 1)$nn.dist[, 1]
  for (wm in c(4000L, 120L)) {
    set.seed(4); a <- core(xy, tg, k = 5, ward_max = wm, use_fnn = TRUE)
    set.seed(4); b <- core(xy, tg, k = 5, ward_max = wm, use_fnn = FALSE)
    expect_identical(a$fold, b$fold, info = paste("ward_max", wm))
    expect_equal(a$candidates$W, b$candidates$W, info = paste("ward_max", wm))
    expect_equal(a$realised, b$realised, info = paste("ward_max", wm))
    expect_equal(a$loo_distances, b$loo_distances, info = paste("ward_max", wm))
  }
})

test_that("the nearest-other-point distance without FNN is exact, duplicates included", {
  loo <- kn_ns(".loo_dist_xy")
  for (n in c(2L, 3L, 5L, 64L, 97L)) {
    set.seed(n); xy <- cbind(runif(n), runif(n))
    D <- as.matrix(stats::dist(xy)); diag(D) <- Inf
    expect_equal(loo(xy, use_fnn = FALSE), unname(apply(D, 1, min)), info = paste("n", n))
  }
  # A co-located pair is at distance 0 from each other.
  set.seed(1); xy <- cbind(runif(20), runif(20)); xy <- rbind(xy, xy[c(3, 3, 7), ])
  D <- as.matrix(stats::dist(xy)); diag(D) <- Inf
  expect_equal(loo(xy, use_fnn = FALSE), unname(apply(D, 1, min)))
  skip_if_not_installed("FNN")
  expect_equal(loo(xy, use_fnn = TRUE), unname(apply(D, 1, min)))
})

test_that("Ward's method continued from groups reproduces the cuts of the full tree", {
  set.seed(14); xy <- cbind(runif(300), runif(300))
  full <- stats::hclust(stats::dist(xy), method = "ward.D2")
  grp  <- stats::cutree(full, k = 60)
  part <- kn_ns(".ward_from_groups")(xy, grp)
  for (q in c(2L, 3L, 7L, 20L, 45L, 59L)) {
    a <- stats::cutree(full, k = q); b <- stats::cutree(part, k = q)[grp]
    # The same partition, whatever the labels.
    expect_identical(length(unique(paste(a, b))), q, info = paste("q", q))
  }
})

test_that("largest-first packing gives each block to the fold with the fewest points", {
  pack <- kn_ns(".pack_largest_first")
  set.seed(1); p <- pack(c(5L, 4L, 3L, 3L, 1L), 2L)
  expect_identical(p$order, c(1L, 2L, 3L, 4L, 5L))
  # 5 | 4, then 3 to the fold at 4, 3 to the fold at 5, 1 to either fold at 8.
  expect_identical(sort(p$loads), c(8L, 8L))
  expect_true(p$fold_of[1L] != p$fold_of[2L])
  expect_identical(p$fold_of[3L], p$fold_of[2L])
  expect_identical(p$fold_of[4L], p$fold_of[1L])
  expect_identical(as.integer(tapply(c(5L, 4L, 3L, 3L, 1L), factor(p$fold_of, 1:2), sum)),
                   p$loads)
  # One block per fold when there are as many folds as blocks.
  set.seed(2); q <- pack(c(2L, 9L, 4L), 3L)
  expect_setequal(q$fold_of, 1:3)
  expect_identical(sort(q$loads), c(2L, 4L, 9L))
})

test_that("co-located samples and few distinct locations do not break knndm", {
  set.seed(15)
  base <- cbind(runif(30, 0, 1e4), runif(30, 0, 1e4))
  xy <- base[sample(30, 120, TRUE), ]            # repeated measurements
  f <- r2_quiet(suppressWarnings(make_folds(kn_pts(xy), k = 4, method = "knndm",
                                            prediction_points = kn_box_grid(), seed = 1)))
  expect_setequal(unlist(lapply(f$folds, `[[`, "test")), seq_len(120))
  expect_true(all(is.finite(f$params$realised_distances)))
  expect_true(all(f$params$candidates$n_clusters < 30L, na.rm = TRUE))
  # Fewer distinct locations than folds leaves the random partition only.
  two <- kn_pts(rbind(base[rep(1, 6), ], base[rep(2, 6), ]))
  g <- r2_quiet(suppressWarnings(make_folds(two, k = 4, method = "knndm",
                                            prediction_points = kn_box_grid(m = 6), seed = 1)))
  expect_identical(g$params$chosen, "random")
  expect_true(is.na(g$params$clustering))
  expect_identical(nrow(g$params$candidates), 1L)
})

test_that("fold_separation reads knndm folds against the target they carry", {
  pts <- kn_pts(kn_clustered(300, 16))
  f <- r2_quiet(make_folds(pts, k = 5, method = "knndm",
                           prediction_points = kn_box_grid(), seed = 1))
  fs <- r2_quiet(fold_separation(f, pts))
  tg <- attr(fs, "target")
  expect_identical(attr(fs, "method"), "knndm")
  expect_identical(tg$source, "the prediction points the kNNDM folds were built for")
  expect_identical(tg$basis, "all sites")
  expect_identical(tg$n_prediction_points, 900L)
  expect_equal(tg$w1, f$params$wasserstein)
  expect_equal(tg$held_out_mean, f$params$realised_mean)
  expect_false(tg$nndm_out_of_reach)
  expect_true(tg$reading %in% c("closer", "farther", "matched", "mixed"))
  out <- paste(capture.output(print(fs)), collapse = " ")
  expect_match(out, "Fold separation: knndm, 5 fold")
  expect_match(out, "the prediction points the kNNDM folds were built for")
  # Out of reach: the print says why no fold design gets there.
  xy <- kn_uniform(150, 17)
  set.seed(170); near <- kn_pts(xy + matrix(rnorm(300, 0, 5), 150))
  nf <- r2_quiet(suppressWarnings(make_folds(kn_pts(xy), k = 5, method = "knndm",
                                             prediction_points = near, seed = 1)))
  ns <- r2_quiet(fold_separation(nf, kn_pts(xy)))
  expect_true(attr(ns, "target")$nndm_out_of_reach)
  expect_identical(attr(ns, "target")$reading, "farther")
  expect_match(gsub("\\s+", " ", paste(capture.output(print(ns)), collapse = " ")),
               "kNNDM cannot match prediction points this near to the data")
})

test_that("knndm folds run through cv_spatial, and on even samples score like random folds", {
  # The mean over draws, not each draw: a single draw differs by several
  # percent either way.
  knn_fit <- function(train_sf, ...) new_spatial_fit(
    subclass = "kn_knnsurf_fit",
    engine = list(xy = sf::st_coordinates(train_sf), y = train_sf$z),
    formula = z ~ 1, response_var = "z", predictor_vars = character(0),
    data_sf = train_sf)
  predict.kn_knnsurf_fit <- function(object, newdata = NULL, ...) {
    e <- object$engine
    q <- if (is.null(newdata)) e$xy else sf::st_coordinates(newdata)
    vapply(seq_len(nrow(q)), function(i) {
      d <- (e$xy[, 1] - q[i, 1])^2 + (e$xy[, 2] - q[i, 2])^2
      mean(e$y[order(d)[seq_len(min(5L, length(d)))]])
    }, 0)
  }
  registerS3method("predict", "kn_knnsurf_fit", predict.kn_knnsurf_fit)
  ratio <- vapply(1:6, function(seed) {
    set.seed(seed); n <- 200
    xy <- cbind(runif(n, 0, 1000), runif(n, 0, 1000))
    z <- as.numeric(t(chol(exp(-as.matrix(stats::dist(xy)) / 60) + diag(1e-8, n))) %*% rnorm(n)) +
      rnorm(n, sd = 0.3)
    pts <- sf::st_as_sf(data.frame(x = 5e5 + xy[, 1], y = 5e6 + xy[, 2], z = z),
                        coords = c("x", "y"), crs = 32632)
    grid <- kn_pts(as.matrix(expand.grid(seq(17, 983, length.out = 25),
                                         seq(17, 983, length.out = 25))))
    fk <- r2_quiet(make_folds(pts, k = 5, method = "knndm", prediction_points = grid, seed = seed))
    fr <- make_folds(pts, k = 5, method = "random_kfold", seed = seed)
    rmse <- function(f) r2_quiet(suppressWarnings(cv_spatial(
      pts, "z", character(0), fit_fn = knn_fit, folds = f, seed = 1)))$overall$RMSE
    ck <- r2_quiet(suppressWarnings(cv_spatial(pts, "z", character(0), fit_fn = knn_fit,
                                               folds = fk, seed = 1)))
    expect_identical(ck$n_folds_succeeded, 5L)
    ck$overall$RMSE / rmse(fr)
  }, 0)
  expect_lt(abs(mean(ratio) - 1), 0.05)
})

test_that("the messages that used to stop at the nndm row limit now name knndm", {
  set.seed(18)
  big <- kn_pts(cbind(runif(5001, 0, 1e5), runif(5001, 0, 1e5)))
  grid <- kn_box_grid(0, 1e5, m = 20)
  expect_error(make_folds(big, method = "nndm", prediction_points = grid),
               "exceeds the safety threshold of 5000.*method = \"knndm\", which matches k folds")
  other <- kn_ns(".block_other_method")
  expect_match(other(5001L), "method = \"knndm\" does not depend on a block size")
  expect_match(other(5001L), "takes at most 5000 rows, and this layer has 5001")
  expect_match(other(5000L), "method = \"nndm\" does not depend on a block size")
  # And the method that was pointed to runs on that layer.
  f <- r2_quiet(suppressWarnings(make_folds(big, k = 5, method = "knndm",
                                            prediction_points = grid, seed = 1)))
  expect_identical(f$k, 5L)
  expect_identical(f$params$clustering, "kmeans+ward")
  expect_lte(f$params$wasserstein, f$params$wasserstein_random)
  expect_lte(f$params$balance_ratio, 3)
  expect_setequal(unlist(lapply(f$folds, `[[`, "test")), seq_len(5001))
})
