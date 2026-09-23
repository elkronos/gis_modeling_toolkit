# fold_separation() measures what blocking is for: how far a held-out point
# ended up from the training data.  The existing check compares the block size
# against the range, which is a statement about the design; these tests pin
# the measurement of the result.

make_pts <- function(n = 200, seed = 1) {
  set.seed(seed)
  sf::st_as_sf(data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000)),
               coords = c("x", "y"), crs = 32632)
}

test_that("blocked folds separate the hold-out further than random ones", {
  pts <- make_pts()
  rnd <- make_folds(pts, k = 4, method = "random_kfold", seed = 1)
  blk <- make_folds(pts, k = 4, method = "block_kfold", block_size = 250, seed = 1)

  s_rnd <- fold_separation(rnd, pts)
  s_blk <- fold_separation(blk, pts)

  expect_s3_class(s_rnd, "fold_separation")
  expect_equal(nrow(s_rnd), 4L)
  expect_equal(sum(s_rnd$n_test), nrow(pts))       # every row held out once
  expect_true(all(s_rnd$n_train + s_rnd$n_test == nrow(pts)))

  # The substantive claim: blocking moves the distances up.
  expect_gt(stats::median(s_blk$median_dist), stats::median(s_rnd$median_dist))
  expect_gt(min(s_blk$min_dist), min(s_rnd$min_dist))

  # Random folds have no blocks; blocked folds report how many each holds.
  expect_true(all(is.na(s_rnd$n_blocks)))
  expect_true(all(s_blk$n_blocks >= 1L))
  expect_equal(sum(s_blk$n_blocks), blk$params$blocks_used)
})

test_that("the distances are the real nearest-training distances", {
  pts <- make_pts(n = 60, seed = 4)
  f <- make_folds(pts, k = 3, method = "random_kfold", seed = 2)
  s <- fold_separation(f, pts)

  # Recompute fold 1 by hand from the coordinates, with no package help.
  xy <- sf::st_coordinates(pts)
  te <- f$folds[[1]]$test; tr <- f$folds[[1]]$train
  d  <- apply(xy[te, , drop = FALSE], 1L, function(p)
    min(sqrt((xy[tr, 1] - p[1])^2 + (xy[tr, 2] - p[2])^2)))
  expect_equal(s$min_dist[1], min(d))
  expect_equal(s$median_dist[1], stats::median(d))
})

test_that("a supplied range gives the share of the hold-out inside it", {
  pts <- make_pts()
  f <- make_folds(pts, k = 4, method = "random_kfold", seed = 1)
  s <- fold_separation(f, pts, sac = 200)
  expect_false(any(is.na(s$within_range)))
  expect_true(all(s$within_range >= 0 & s$within_range <= 1))
  # On random folds over a 1000-unit square, every held-out point has a
  # training neighbour far closer than 200 units.
  expect_equal(unique(s$within_range), 1)
  # Without a range the column is present but empty, and the print method
  # then omits it rather than showing a column of NAs.
  s0 <- fold_separation(f, pts)
  expect_true(all(is.na(s0$within_range)))
  expect_false(grepl("within_range", paste(utils::capture.output(print(s0)),
                                           collapse = " ")))
  expect_true(grepl("within_range", paste(utils::capture.output(print(s)),
                                          collapse = " ")))
})

test_that("row identity follows ..row_id, not row position", {
  pts <- make_pts(n = 80, seed = 5)
  f <- make_folds(pts, k = 4, method = "block_kfold", block_size = 300, seed = 1)
  base <- fold_separation(f, pts)

  # Same layer, rows shuffled, carrying the ids make_folds assigned.  The
  # measurement must not move: matching by position instead would pair the
  # wrong points and change every distance.
  pts_id <- pts; pts_id$..row_id <- seq_len(nrow(pts_id))
  shuffled <- pts_id[sample(nrow(pts_id)), , drop = FALSE]
  expect_equal(fold_separation(f, shuffled)$min_dist, base$min_dist)
  expect_equal(fold_separation(f, shuffled)$median_dist, base$median_dist)
})

test_that("it takes a bare list of splits and refuses anything else", {
  pts <- make_pts(n = 50, seed = 6)
  f <- make_folds(pts, k = 3, method = "random_kfold", seed = 1)
  expect_equal(fold_separation(f$folds, pts)$min_dist,
               fold_separation(f, pts)$min_dist)
  expect_error(fold_separation(f, "not sf"), "must be an sf object")
  expect_error(fold_separation(list(), pts), "make_folds\\(\\) result")
  expect_error(fold_separation(list(list(bad = 1)), pts), "make_folds\\(\\) result")
})

test_that("fold entries naming rows the layer does not have are skipped", {
  pts <- make_pts(n = 60, seed = 7)
  f <- make_folds(pts, k = 3, method = "random_kfold", seed = 1)
  s <- fold_separation(f, pts[1:40, , drop = FALSE])
  expect_gt(attr(s, "n_unknown_ids"), 0L)
  expect_true(all(is.finite(s$n_test)))
})
