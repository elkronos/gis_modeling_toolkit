# tests/testthat/test-fix-folds-separation-target.R
# ---------------------------------------------------------------------------
# fold_separation() printed "this score is optimistic" whenever more than
# half the held-out points had a training point inside the range.  That share
# cannot give a direction: the same 100 percent goes with an unbiased score
# on evenly sited data and with a score 20 percent low on clustered data.
# Without a prediction target it now reports distances and no verdict; with
# one it compares the two sets of distances and words the reading by the
# signed gap.  Each test fails on the code before the fix.
# ---------------------------------------------------------------------------

st_pts <- function(xy) sf::st_as_sf(data.frame(x = 2e5 + xy[, 1], y = 4e6 + xy[, 2]),
                                    coords = c("x", "y"), crs = 32617)
st_even <- function(n = 300, w = 3e5, h = 3e5, seed = 1) {
  set.seed(seed); st_pts(cbind(runif(n, 0, w), runif(n, 0, h)))
}
st_clustered <- function(n = 300, w = 3e5, h = 3e5, seed = 1) {
  set.seed(seed)
  c0 <- cbind(runif(8, 0.1 * w, 0.9 * w), runif(8, 0.1 * h, 0.9 * h)); i <- sample(8, n, TRUE)
  xy <- c0[i, ] + cbind(rnorm(n, 0, 0.04 * w), rnorm(n, 0, 0.04 * w))
  st_pts(cbind(pmin(pmax(xy[, 1], 0), w), pmin(pmax(xy[, 2], 0), h)))
}
st_grid <- function(pts, m = 30) {
  bb <- sf::st_bbox(pts)
  sf::st_as_sf(expand.grid(x = seq(bb[["xmin"]], bb[["xmax"]], length.out = m),
                           y = seq(bb[["ymin"]], bb[["ymax"]], length.out = m)),
               coords = c("x", "y"), crs = sf::st_crs(pts))
}
st_flat <- function(x) gsub("\\s+", " ", paste(utils::capture.output(print(x)), collapse = " "))
st_q <- function(expr) r2_quiet(suppressWarnings(suppressMessages(expr)))
# Nearest-neighbour distance from each row of a to the rows of b, by hand.
st_nn <- function(a, b) apply(a, 1L, function(p) min(sqrt((b[, 1] - p[1])^2 + (b[, 2] - p[2])^2)))

test_that("without a prediction target there are distances and no verdict", {
  pts <- st_even(150, 8e5, 3e5, seed = 101)
  rnd <- st_q(make_folds(pts, k = 5, method = "random_kfold", seed = 1))
  blk <- st_q(make_folds(pts, k = 5, method = "block_kfold", block_size = 1e5, seed = 1))
  buf <- st_q(make_folds(pts, method = "buffered_loo", buffer = 3e4))
  for (f in list(rnd, blk, buf)) {
    s <- fold_separation(f, pts, sac = 1e5)
    expect_null(attr(s, "target"))
    expect_true(all(is.na(s$target_median_dist)) && all(is.na(s$dist_gap)))
    out <- st_flat(s)
    # The old print called all three "optimistic" and named a remedy; on this
    # layout random folds are about unbiased and wider blocks read higher.
    expect_no_match(out, "optimistic")
    expect_no_match(out, "leak")
    expect_no_match(out, "widen the")
    expect_no_match(out, "use blocked or buffered folds")
    expect_match(out, "No prediction target is given, so these are distances and no conclusion about the score")
    expect_match(out, "Pass `prediction_points` to compare")
  }
  # No range either: still the distances and the same closing line.
  expect_match(st_flat(fold_separation(rnd, pts)), "No prediction target is given")
})

test_that("within_range is the correlation-weighted share, not a count inside the range", {
  pts <- st_even(120, seed = 7)
  f <- st_q(make_folds(pts, k = 4, method = "random_kfold", seed = 2))
  rng <- 6e4
  s <- fold_separation(f, pts, sac = rng)
  xy <- sf::st_coordinates(pts)
  for (j in 1:4) {
    d <- st_nn(xy[f$folds[[j]]$test, , drop = FALSE], xy[f$folds[[j]]$train, , drop = FALSE])
    expect_equal(s$within_range[j], mean(exp(-3 * d / rng)))
    # Every one of these points is inside the range, which the old column
    # reported as 1 for each fold whatever the distances were.
    expect_true(all(d < rng))
    expect_lt(s$within_range[j], 1)
  }
  expect_identical(attr(s, "range_model"), "Exp")
  d_all <- unlist(lapply(f$folds, function(sp) st_nn(xy[sp$test, , drop = FALSE], xy[sp$train, , drop = FALSE])))
  expect_equal(attr(s, "pooled")$within_range, mean(exp(-3 * d_all / rng)))
  expect_equal(attr(s, "pooled")$median_dist, stats::median(d_all))
  # A point at a twentieth of the range weighs more than one at nine tenths.
  expect_gt(exp(-3 * 0.05), 5 * exp(-3 * 0.9))

  # A range that carries a spherical fit is weighted by the spherical model.
  sph <- structure(rng, class = c("sac_range", "numeric"), crs = sf::st_crs(pts),
                   variogram_model = data.frame(model = c("Nug", "Sph"), psill = c(0.1, 1),
                                                range = c(0, rng)))
  s2 <- fold_separation(f, pts, sac = sph)
  expect_identical(attr(s2, "range_model"), "Sph")
  d1 <- st_nn(xy[f$folds[[1]]$test, , drop = FALSE], xy[f$folds[[1]]$train, , drop = FALSE])
  expect_equal(s2$within_range[1], mean(1 - 1.5 * (d1 / rng) + 0.5 * (d1 / rng)^3))
})

test_that("the three clear cases read matched, farther and closer", {
  # Evenly sited points, a grid over the same box as the target.
  even <- st_even(seed = 11); grid <- st_grid(even)
  rnd <- st_q(make_folds(even, k = 5, method = "random_kfold", seed = 1))
  blk <- st_q(make_folds(even, k = 5, method = "block_kfold", block_size = 1e5, seed = 1))
  # With a range the reading is on the correlation scale, without one on the
  # mean distance; the clear cases agree on both.
  for (sac in list(1e5, NULL)) {
    a <- attr(fold_separation(rnd, even, sac = sac, prediction_points = grid), "target")
    expect_identical(a$reading, "matched")
    expect_identical(a$scale, if (is.null(sac)) "distance" else "correlation")
    b <- attr(fold_separation(blk, even, sac = sac, prediction_points = grid), "target")
    expect_identical(b$reading, "farther")
    expect_gt(b$mean_diff, 0)
  }
  # Clustered points, a grid over the whole box: random folds hold a point
  # out next to its cluster while the map is mostly predicted between them.
  cl <- st_clustered(seed = 12); gcl <- st_grid(cl)
  rcl <- st_q(make_folds(cl, k = 5, method = "random_kfold", seed = 1))
  for (sac in list(1e5, NULL)) {
    d <- attr(fold_separation(rcl, cl, sac = sac, prediction_points = gcl), "target")
    expect_identical(d$reading, "closer")
    expect_lt(d$mean_diff, 0)
  }

  # The print words each as a statement about distances.
  expect_match(st_flat(fold_separation(rnd, even, sac = 1e5, prediction_points = grid)),
               "about as far from their training data as the prediction points do")
  out_far <- st_flat(fold_separation(blk, even, sac = 1e5, prediction_points = grid))
  expect_match(out_far, "farther from their training data than the prediction points do, so the cross-validated error tends to read high")
  out_near <- st_flat(fold_separation(rcl, cl, sac = 1e5, prediction_points = gcl))
  expect_match(out_near, "closer to their training data than the prediction points do, so the cross-validated error tends to read low")
  for (o in c(out_far, out_near)) {
    expect_no_match(o, "optimistic"); expect_no_match(o, "leak")
    expect_match(o, "Held-out minus prediction: [+-][0-9.e+]+ in mean distance \\(W1 [0-9.e+]+\\); [+-][0-9.]+ on the correlation scale")
  }
  # The old print said "optimistic: widen the blocks" of the blocked design
  # that already reads high on this layout.
  expect_no_match(out_far, "widen")
})

test_that("the target distances are the ones a hand computation gives", {
  pts <- st_even(90, seed = 21); grid <- st_grid(pts, m = 12)
  f <- st_q(make_folds(pts, k = 3, method = "random_kfold", seed = 4))
  s <- fold_separation(f, pts, prediction_points = grid)
  tg <- attr(s, "target")
  xy <- sf::st_coordinates(pts); gxy <- sf::st_coordinates(grid)
  expect_identical(tg$n_prediction_points, 144L)
  expect_identical(tg$basis, "each fold's training data")
  expect_identical(tg$source, "prediction_points")
  # Each fold's prediction distances are to that fold's own training rows.
  per <- lapply(f$folds, function(sp) list(
    cv = st_nn(xy[sp$test, , drop = FALSE], xy[sp$train, , drop = FALSE]),
    tg = st_nn(gxy, xy[sp$train, , drop = FALSE])))
  for (j in 1:3) {
    expect_equal(s$target_median_dist[j], stats::median(per[[j]]$tg))
    expect_equal(s$dist_gap[j], mean(per[[j]]$cv) - mean(per[[j]]$tg))
  }
  # Pooled: folds count by their share of the held-out points.
  wts <- vapply(per, function(p) length(p$cv), 1) / nrow(pts)
  expect_equal(tg$held_out_mean, mean(unlist(lapply(per, `[[`, "cv"))))
  expect_equal(tg$target_mean, sum(wts * vapply(per, function(p) mean(p$tg), 1)))
  expect_equal(tg$mean_diff, tg$held_out_mean - tg$target_mean)
  # The distances to all sites sit beside them, and are shorter: each fold
  # trains on two thirds of the sites.
  all_d <- st_nn(gxy, xy)
  expect_equal(tg$target_all_median, stats::median(all_d))
  expect_equal(tg$target_all_mean, mean(all_d))
  expect_lt(tg$target_all_mean, tg$target_mean)
  # W1 is the unsigned area between the two distributions, never below the
  # signed one.
  expect_gte(tg$w1, abs(tg$mean_diff))
  expect_true(is.na(tg$cor_gap))
})

test_that("prediction_points takes points, cells or an sfc, in any CRS", {
  pts <- st_even(90, seed = 22); grid <- st_grid(pts, m = 10)
  f <- st_q(make_folds(pts, k = 3, method = "random_kfold", seed = 4))
  base <- attr(fold_separation(f, pts, prediction_points = grid), "target")
  expect_equal(attr(fold_separation(f, pts, prediction_points = sf::st_geometry(grid)), "target")$mean_diff,
               base$mean_diff)
  ll <- attr(fold_separation(f, pts, prediction_points = sf::st_transform(grid, 4326)), "target")
  expect_equal(ll$mean_diff, base$mean_diff, tolerance = 1e-6)
  # Grid cells are reduced to one point each.
  cells <- sf::st_as_sf(sf::st_make_grid(pts, n = c(10, 10)))
  expect_identical(attr(st_q(fold_separation(f, pts, prediction_points = cells)), "target")$n_prediction_points,
                   100L)
  expect_error(fold_separation(f, pts, prediction_points = data.frame(x = 1)),
               "`prediction_points` must be an sf or sfc layer")
  expect_error(fold_separation(f, pts, prediction_points = grid[0, ]),
               "`prediction_points` has no rows")
})

test_that("leave-one-out schemes are measured against all sites", {
  pts <- st_even(80, seed = 23); grid <- st_grid(pts, m = 10)
  loo <- st_q(make_folds(pts, method = "buffered_loo", buffer = 2e4))
  s <- fold_separation(loo, pts, prediction_points = grid)
  tg <- attr(s, "target")
  expect_identical(tg$basis, "all sites")
  xy <- sf::st_coordinates(pts)
  all_d <- st_nn(sf::st_coordinates(grid), xy)
  expect_equal(tg$target_mean, mean(all_d))
  expect_equal(unique(s$target_median_dist), stats::median(all_d))
  # Eighty one-point folds are not printed row by row.
  out <- utils::capture.output(print(s))
  expect_true(any(grepl("and 70 more fold\\(s\\)", out)))
  expect_lt(length(out), 40L)
})

test_that("crossing distributions read as mixed, and no reading names a verdict", {
  rd <- spatialkit:::.fold_sep_reading
  expect_identical(rd(-0.30, 0.30, 0.1), "closer")
  expect_identical(rd(0.30, 0.31, 0.1), "farther")
  expect_identical(rd(0.02, 0.05, 0.1), "matched")
  # A net shift that is the remainder of two opposite ones names no side.
  expect_identical(rd(-0.11, 0.50, 0.1), "mixed")
  expect_identical(rd(0.02, 0.40, 0.1), "mixed")
  expect_true(is.na(rd(NA_real_, 0.1, 0.1)))

  # Two samples that cross: equal means, different spreads.
  cmp <- spatialkit:::.fold_sep_compare(c(1, 9), c(4, 6))
  expect_equal(cmp$diff, 0)
  expect_equal(cmp$w1, 3)
  # One wholly to one side of the other: the two areas are equal.
  cmp2 <- spatialkit:::.fold_sep_compare(c(1, 2, 3), c(11, 12, 13))
  expect_equal(cmp2$diff, -10); expect_equal(cmp2$w1, 10)
  # Weights: the second value counted three times is a sample of four.
  expect_equal(spatialkit:::.fold_sep_compare(1, c(0, 4), wb = c(1, 3))$mean_b, 3)

  adv <- spatialkit:::.fold_separation_advice
  tg <- function(reading) list(reading = reading, basis = "each fold's training data",
                               target_median = 10, target_all_median = 9, mean_diff = -3,
                               w1 = 4, cor_gap = -0.2, cor_w1 = 0.2)
  for (rdg in c("closer", "farther", "matched", "mixed")) for (m in c("block_kfold", "nndm")) {
    txt <- adv(tg(rdg), m)
    expect_no_match(txt, "optimistic"); expect_no_match(txt, "leak")
    expect_no_match(txt, "widen"); expect_no_match(txt, "use blocked")
  }
  expect_match(adv(tg("mixed"), "block_kfold"), "give no single direction")
  expect_match(adv(NULL, "block_kfold"), "No prediction target is given")
  expect_match(adv(tg("matched"), "random_kfold"), "do not make a score unbiased by themselves")
})
