# tests/testthat/test-fix-folds-zero-target.R
# ---------------------------------------------------------------------------
# fold_separation() with a prediction target it could not put a reading on.
#
# Prediction points on top of the sites (a model scored where it was sampled)
# have a prediction mean of 0.  On the distance scale the gap was divided by
# that mean, so it was infinite, the reading was NA, and print() closed with
# "No prediction target is given ... Pass `prediction_points`" under a header
# that named the target.  NNDM folds measured in another CRS than the one
# their target distances are in had the target dropped and printed the same
# sentence.  Each test fails on the code before the fix.
# ---------------------------------------------------------------------------

zt_pts <- function(n = 80, seed = 7) {
  set.seed(seed)
  sf::st_as_sf(data.frame(x = 5e5 + runif(n, 0, 1e4), y = 5e6 + runif(n, 0, 1e4)),
               coords = c("x", "y"), crs = 32632)
}
zt_flat <- function(x) gsub("\\s+", " ", paste(utils::capture.output(print(x)), collapse = " "))
zt_q <- function(expr) r2_quiet(suppressWarnings(suppressMessages(expr)))

test_that("prediction points on the samples read 'farther' without a range", {
  u <- zt_pts()
  nn <- zt_q(make_folds(u, method = "nndm", prediction_points = u))
  expect_identical(nn$params$target_median, 0)
  expect_true(nn$params$target_out_of_reach)
  designs <- list(
    # The target the NNDM folds carry, and the same points passed again.
    nndm_carried = list(folds = nn, pp = NULL),
    nndm_passed  = list(folds = nn, pp = u),
    buffered_loo = list(folds = zt_q(make_folds(u, method = "buffered_loo", buffer = 800)),
                        pp = u),
    random_60    = list(folds = zt_q(make_folds(u, k = 60, method = "random_kfold", seed = 1)),
                        pp = u))
  for (nm in names(designs)) {
    s <- fold_separation(designs[[nm]]$folds, u, prediction_points = designs[[nm]]$pp)
    tg <- attr(s, "target")
    expect_identical(tg$basis, "all sites", info = nm)
    expect_identical(tg$scale, "distance", info = nm)
    # Every prediction point is a site, so every prediction distance is 0 ...
    expect_identical(tg$target_mean, 0, info = nm)
    expect_identical(tg$target_median, 0, info = nm)
    # ... and the held-out distances lie wholly above them: the gap is the
    # held-out mean, and it is all of the difference between the two.
    expect_gt(tg$mean_diff, 0)
    expect_equal(tg$mean_diff, tg$held_out_mean, info = nm)
    expect_equal(tg$w1, tg$mean_diff, info = nm)
    # It was NA, from a gap of Inf.
    expect_identical(tg$reading, "farther", info = nm)
    out <- zt_flat(s)
    expect_match(out, "prediction target: 80 point\\(s\\)", info = nm)
    # The closing lines used to deny, under that header, that there was one.
    expect_no_match(out, "No prediction target is given", info = nm)
    expect_no_match(out, "Pass `prediction_points`", info = nm)
    expect_match(out, "The prediction points sit a median of 0 from the nearest of all sites",
                 info = nm)
    expect_match(out, "farther from their training data than the prediction points do, so the cross-validated error tends to read high",
                 info = nm)
  }
  # The NNDM folds say why nothing could be done about it.
  expect_match(zt_flat(fold_separation(nn, u)), "NNDM can only lengthen a held-out distance")
  # With a range the reading was already right; the two scales now agree.
  expect_identical(attr(fold_separation(nn, u, sac = 3000), "target")$reading, "farther")
})

test_that("held-out and prediction distances that are all 0 read 'matched'", {
  # Twenty locations measured twice, the two measurements in opposite folds:
  # every held-out point has its twin in the training data, and every
  # prediction point (the sites again) lies on a training point too.
  set.seed(9)
  xy <- cbind(runif(20, 0, 1e4), runif(20, 0, 1e4))
  d <- sf::st_as_sf(data.frame(x = 5e5 + rep(xy[, 1], 2), y = 5e6 + rep(xy[, 2], 2)),
                    coords = c("x", "y"), crs = 32632)
  splits <- list(list(train = 21:40, test = 1:20), list(train = 1:20, test = 21:40))
  s <- fold_separation(splits, d, prediction_points = d)
  tg <- attr(s, "target")
  expect_identical(tg$held_out_mean, 0)
  expect_identical(tg$target_mean, 0)
  expect_identical(tg$mean_diff, 0)
  expect_identical(tg$w1, 0)
  expect_identical(tg$reading, "matched")       # was NA, from 0 / 0
  out <- zt_flat(s)
  expect_no_match(out, "No prediction target is given")
  expect_match(out, "about as far from their training data as the prediction points do")
})

test_that("the reading rule takes a tolerance of 0", {
  rd <- spatialkit:::.fold_sep_reading
  # Identical distributions are matched whatever the tolerance.
  expect_identical(rd(0, 0, 0), "matched")
  expect_identical(rd(0, 0, 0.1), "matched")
  # With no tolerance any net shift that is most of the difference has a side.
  expect_identical(rd(589.1, 589.1, 0), "farther")
  expect_identical(rd(-3, 3, 0), "closer")
  # Equal means and different spreads still cross.
  expect_identical(rd(0, 3, 0), "mixed")
  # The cases with a tolerance are read as before.
  expect_identical(rd(-0.30, 0.30, 0.1), "closer")
  expect_identical(rd(0.02, 0.05, 0.1), "matched")
  expect_identical(rd(-0.11, 0.50, 0.1), "mixed")
})

test_that("a share of the prediction mean is read the same in distance units", {
  # The distance-scale reading is now taken on the distances themselves,
  # against a quarter of the prediction mean, in place of ratios to that
  # mean.  Where the mean is above 0 the two are the same rule.
  set.seed(31)
  pts <- sf::st_as_sf(data.frame(x = 2e5 + runif(200, 0, 3e5), y = 4e6 + runif(200, 0, 3e5)),
                      coords = c("x", "y"), crs = 32617)
  grid <- sf::st_as_sf(sf::st_make_grid(pts, n = c(20, 20), what = "centers"))
  rd <- spatialkit:::.fold_sep_reading
  for (f in list(zt_q(make_folds(pts, k = 5, method = "random_kfold", seed = 1)),
                 zt_q(make_folds(pts, k = 5, method = "block_kfold", block_size = 1e5, seed = 1)),
                 zt_q(make_folds(pts, method = "buffered_loo", buffer = 2e4)))) {
    tg <- attr(fold_separation(f, pts, prediction_points = grid), "target")
    expect_gt(tg$target_mean, 0)
    expect_identical(tg$tolerance, 0.25)
    expect_identical(tg$reading,
                     rd(tg$mean_diff / tg$target_mean, tg$w1 / tg$target_mean, 0.25))
  }
})

test_that("NNDM folds measured in another CRS say their target is not compared", {
  u <- zt_pts(150, seed = 3)
  g <- sf::st_as_sf(sf::st_make_grid(u, n = c(25, 25), what = "centers"))
  f <- zt_q(make_folds(u, method = "nndm", prediction_points = g))
  expect_identical(f$params$crs, "EPSG:32632")

  # Measured where the folds were built, the carried target is compared.
  s0 <- fold_separation(f, u)
  expect_type(attr(s0, "target"), "list")
  expect_null(attr(s0, "not_compared"))

  # The same rows in another projected CRS: the carried distances are lengths
  # in EPSG:32632 and cannot be set beside distances in EPSG:3857.
  u2 <- sf::st_transform(u, 3857)
  sac_other <- structure(3000, class = c("sac_range", "numeric"), crs = sf::st_crs(3857))
  for (s in list(zt_q(fold_separation(f, u2)), zt_q(fold_separation(f, u, sac = sac_other)))) {
    expect_identical(attr(s, "crs"), "EPSG:3857")
    # NULL, and not the sentence below: attr() matches names partially, so
    # the reason is kept under a name that does not begin with "target".
    expect_null(attr(s, "target"))
    note <- attr(s, "not_compared")
    expect_type(note, "character")
    expect_match(note, "625 prediction points")
    expect_match(note, "EPSG:32632")
    expect_match(note, "EPSG:3857")
    out <- zt_flat(s)
    # These folds were built for a target; the print used to say none was given.
    expect_no_match(out, "No prediction target is given")
    expect_no_match(out, "prediction target: ")
    expect_match(out, "so the two are not compared")
    expect_match(out, "Pass `prediction_points` to have the target measured here")
  }
  # Following that advice gives the reading the folds have where they were built.
  s2 <- fold_separation(f, u2, prediction_points = g)
  expect_null(attr(s2, "not_compared"))
  expect_identical(attr(s2, "target")$reading, attr(s0, "target")$reading)
  expect_no_match(zt_flat(s2), "not compared")

  # Other fold schemes carry no target, and still say none is given.
  r <- zt_q(make_folds(u, k = 5, method = "random_kfold", seed = 1))
  sr <- fold_separation(r, u2)
  expect_null(attr(sr, "not_compared"))
  expect_match(zt_flat(sr), "No prediction target is given")
})

test_that("the closing lines never deny a target that exists", {
  adv <- spatialkit:::.fold_separation_advice
  none <- adv(NULL, "block_kfold")
  expect_match(none, "No prediction target is given")
  # A target that was not compared, with the reason.
  noted <- adv(NULL, "nndm", note = "These NNDM folds carry X in A and the distances here are in B, so the two are not compared.")
  expect_no_match(noted, "No prediction target is given")
  expect_match(noted, "so the two are not compared\\. Pass `prediction_points` to have the target measured here")
  # A target that was measured but has no reading: the numbers, and no
  # sentence saying there is no target.
  unread <- adv(list(reading = NA_character_, basis = "all sites", target_median = 0,
                     target_all_median = 0, mean_diff = 12, w1 = 12, cor_gap = NA_real_,
                     cor_w1 = NA_real_), "buffered_loo")
  expect_no_match(unread, "No prediction target is given")
  expect_no_match(unread, "Pass `prediction_points`")
  expect_match(unread, "The prediction points sit a median of 0 from the nearest of all sites")
  expect_match(unread, "no reading is drawn")
})
