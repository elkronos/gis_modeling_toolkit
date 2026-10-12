# tests/testthat/test-fix-folds-max-test.R
# ---------------------------------------------------------------------------
# make_folds(max_test = ) for the two leave-one-out methods: folds are built
# for a seeded random subset of the held-out points while every row stays
# available for training.  The model is refitted max_test times and the
# folds hold 4 * n * max_test bytes of row indices in place of 4 * n^2.
# ---------------------------------------------------------------------------

mt_pts <- function(n, seed, side = 1e4) {
  set.seed(seed)
  sf::st_as_sf(data.frame(x = 5e5 + runif(n, 0, side), y = 5e6 + runif(n, 0, side),
                          a = rnorm(n), z = rnorm(n)),
               coords = c("x", "y"), crs = 32632)
}
mt_grid <- function(pts, m = 20) {
  bb <- sf::st_bbox(pts)
  sf::st_as_sf(expand.grid(x = seq(bb[["xmin"]], bb[["xmax"]], length.out = m),
                           y = seq(bb[["ymin"]], bb[["ymax"]], length.out = m)),
               coords = c("x", "y"), crs = sf::st_crs(pts))
}

test_that("buffered_loo with max_test builds that many of the full design's folds", {
  pts <- mt_pts(120, 1)
  full <- r2_quiet(make_folds(pts, method = "buffered_loo", buffer = 1500))
  sub  <- r2_quiet(make_folds(pts, method = "buffered_loo", buffer = 1500,
                              max_test = 25, seed = 7))
  expect_identical(sub$k, 25L)
  expect_length(sub$folds, 25L)
  expect_identical(sub$params$max_test, 25L)
  held <- vapply(sub$folds, function(s) s$test, 0L)
  expect_false(anyDuplicated(held) > 0)
  expect_identical(held, sort(held))
  # Each fold is the fold the full design has for that point: every row is
  # still available for training.
  for (j in seq_along(held)) {
    expect_identical(sub$folds[[j]]$train, full$folds[[held[j]]]$train)
    expect_false(held[j] %in% sub$folds[[j]]$train)
  }
  # Rows that are not held out have no fold.
  expect_identical(sub$assignment$row_id, seq_len(120))
  expect_identical(which(!is.na(sub$assignment$fold)), held)
  expect_identical(sub$assignment$fold[held], seq_len(25))
  expect_identical(sum(is.na(sub$assignment$fold)), 95L)
  # The subset is the seed's.
  again <- r2_quiet(make_folds(pts, method = "buffered_loo", buffer = 1500,
                               max_test = 25, seed = 7))
  expect_identical(again$assignment, sub$assignment)
  other <- r2_quiet(make_folds(pts, method = "buffered_loo", buffer = 1500,
                               max_test = 25, seed = 8))
  expect_false(identical(other$assignment$fold, sub$assignment$fold))
})

test_that("nndm with max_test matches over every point and builds a subset of the folds", {
  pts <- mt_pts(150, 2)
  # A grid three times as wide as the samples, so that neighbours are excluded.
  bb <- sf::st_bbox(pts)
  wide <- sf::st_as_sf(expand.grid(x = seq(bb[["xmin"]] - 1e4, bb[["xmax"]] + 1e4, length.out = 20),
                                   y = seq(bb[["ymin"]] - 1e4, bb[["ymax"]] + 1e4, length.out = 20)),
                       coords = c("x", "y"), crs = 32632)
  full <- r2_quiet(suppressWarnings(make_folds(pts, method = "nndm", prediction_points = wide)))
  sub  <- r2_quiet(suppressWarnings(make_folds(pts, method = "nndm", prediction_points = wide,
                                               max_test = 40, seed = 3)))
  expect_gt(full$params$n_removed_total, 0L)
  expect_identical(sub$k, 40L)
  expect_length(sub$folds, 40L)
  held <- vapply(sub$folds, function(s) s$test, 0L)
  expect_identical(length(unique(held)), 40L)
  for (j in seq_along(held))
    expect_identical(sub$folds[[j]]$train, full$folds[[held[j]]]$train)
  expect_identical(which(!is.na(sub$assignment$fold)), held)
  # The matching, and everything `params` says about it, is the full design's.
  same <- setdiff(names(full$params), "seed")
  expect_identical(sub$params[same], full$params[same])
  expect_length(sub$params$realised_distances, 150L)
  expect_identical(sub$params$max_test, 40L)
  expect_null(full$params$max_test)
  # Without a subset nothing about the result has changed.
  expect_identical(full$k, 150L)
  expect_identical(full$assignment$fold, seq_len(150))
})

test_that("a max_test at or above n holds every point out", {
  pts <- mt_pts(40, 3)
  full <- r2_quiet(make_folds(pts, method = "buffered_loo", buffer = 800))
  all_ <- r2_quiet(make_folds(pts, method = "buffered_loo", buffer = 800, max_test = 40))
  more <- r2_quiet(make_folds(pts, method = "buffered_loo", buffer = 800, max_test = 500))
  expect_identical(all_$folds, full$folds)
  expect_identical(more$folds, full$folds)
  expect_identical(more$assignment, full$assignment)
  expect_identical(more$k, 40L)
})

test_that("max_test is validated, and refused for the k-fold methods", {
  pts <- mt_pts(60, 4)
  for (bad in list(1, 2.5, -3, NA_real_, c(5, 6), "10", Inf))
    expect_error(make_folds(pts, method = "buffered_loo", buffer = 500, max_test = bad),
                 "`max_test` must be a single whole number >= 2")
  for (m in c("random_kfold", "block_kfold", "leave_location_out", "knndm"))
    expect_error(make_folds(pts, k = 5, method = m, max_test = 10,
                            prediction_points = mt_grid(pts)),
                 "`max_test` is only used by the leave-one-out methods", info = m)
})

test_that("max_test lets buffered_loo run where its storage limit refuses the full design", {
  pts <- mt_pts(13000, 5, side = 2e5)
  # 4 * 13000^2 bytes is 0.63 GB of row indices: refused, and the message
  # names the way through with a number that works.
  err <- tryCatch(make_folds(pts, method = "buffered_loo", buffer = 2000),
                  error = function(e) conditionMessage(e))
  expect_match(err, "would produce 13000 leave-one-out splits holding about 0.63 GB")
  expect_match(err, "pass `max_test` to build and score a random subset")
  ok <- as.numeric(sub(".*max_test = ([0-9]+) keeps the splits under 0.6 GB.*", "\\1", err))
  expect_true(is.finite(ok) && ok >= 2 && 4 * 13000 * ok / 1024^3 <= 0.6)
  sub <- r2_quiet(make_folds(pts, method = "buffered_loo", buffer = 2000,
                             max_test = 30, seed = 1))
  expect_identical(sub$k, 30L)
  expect_true(all(vapply(sub$folds, function(s) length(s$train), 0L) > 12000L))
  expect_true(all(vapply(sub$folds, function(s) !(s$test %in% s$train), NA)))
  # The largest subset the message offered is accepted by the size check.
  expect_lte(4 * 13000 * ok / 1024^3, 0.6)
})

test_that("folds built with max_test are scored and measured as they are", {
  pts <- mt_pts(90, 6)
  sub <- r2_quiet(make_folds(pts, method = "buffered_loo", buffer = 1200,
                             max_test = 20, seed = 2))
  cv <- r2_quiet(suppressWarnings(cv_spatial(pts, "z", "a", fit_fn = r2_fit,
                                             folds = sub, seed = 1)))
  expect_identical(cv$n_folds_attempted, 20L)
  expect_identical(cv$n_folds_succeeded, 20L)
  held <- vapply(sub$folds, function(s) s$test, 0L)
  fs <- r2_quiet(fold_separation(sub, pts))
  expect_identical(nrow(fs), 20L)
  expect_identical(sum(fs$n_test), 20L)
  expect_true(all(fs$min_dist > 1200))
  expect_identical(attr(fs, "n_unknown_ids"), 0L)
  # nndm folds carry their target, and the subset is read against it.
  nf <- r2_quiet(suppressWarnings(make_folds(pts, method = "nndm",
                                             prediction_points = mt_grid(pts),
                                             max_test = 30, seed = 2)))
  ns <- r2_quiet(fold_separation(nf, pts))
  expect_identical(nrow(ns), 30L)
  expect_identical(attr(ns, "target")$source,
                   "the prediction points the NNDM folds were built for")
})
