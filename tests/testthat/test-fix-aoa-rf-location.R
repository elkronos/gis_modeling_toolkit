# tests/testthat/test-fix-aoa-rf-location.R
# ---------------------------------------------------------------------------
# area_of_applicability() and models whose predictions depend on location.
#
# The index measured location only for a forest fitted with
# include_coords = TRUE.  A GWR's coefficients and a Bayesian fit's Gaussian
# process are functions of location by construction, and for them the index
# was computed on the covariates alone: a GWR with a slope rising eastward,
# trained on x 0-1000, had 13 times the in-extent RMSE at x 1500-2500 and
# 84% of those points inside the AOA.
#
# Location now has an index of its own for those classes (and under
# coords = "yes"): nearest training point in map space, thresholded by the
# same rule and folds, and a point is inside only if inside on both.  It is
# kept apart from the covariates because as two more equally weighted
# columns it is diluted: with 12 covariates 89% of points 200-500 units
# beyond the extent still read as inside.
# ---------------------------------------------------------------------------

loc_pts <- function(m, x0, x1, p = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  d <- data.frame(x = stats::runif(m, x0, x1), y = stats::runif(m, 0, 1000))
  for (j in seq_len(p)) d[[paste0("a", j)]] <- stats::rnorm(m)
  d$z <- (1 + d$x / 250) * d$a1 + stats::rnorm(m, 0, 0.3)
  sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
}

# A stand-in carrying a class and the fields area_of_applicability() reads,
# for the classes that cannot be fitted in a unit test.
loc_fake_fit <- function(train, vars, cls) {
  structure(list(data_sf = train, predictor_vars = vars, response_var = "z",
                 info = list(), n = nrow(train)),
            class = c(cls, "spatial_fit"))
}


test_that("a GWR fit has location measured: far points are outside", {
  skip_if_not_installed("sp")
  trn <- loc_pts(400, 0, 1000, seed = 3)
  gf  <- suppressWarnings(fit_gwr_model(trn, "z", "a1", bandwidth = 40))
  share <- function(x0, x1) {
    nd <- loc_pts(300, x0, x1)
    a  <- area_of_applicability(nd, model = gf)
    list(a = a, inside = mean(a$aoa$AOA),
         inside_loc = mean(a$aoa$DI_location <= a$location$threshold),
         inside_cov = mean(a$aoa$DI <= a$threshold))
  }
  within <- share(0, 1000)
  far    <- share(1500, 2500)
  farther <- share(4000, 5000)
  near   <- share(1000, 1200)

  # Inside the extent the location index lets nearly every point through,
  # and the AOA is the covariate index's, a few points apart.
  expect_gt(within$inside_loc, 0.9)
  expect_gt(within$inside, within$inside_cov - 0.1)
  expect_gt(within$inside, 0.75)
  # Beyond it the covariates still say "inside" for most points (that was
  # the whole answer before), and location says otherwise.
  expect_gt(far$inside_cov, 0.7)
  expect_lt(far$inside, 0.1)
  expect_lt(farther$inside, 0.1)
  expect_lt(near$inside, 0.5)

  a <- far$a
  expect_identical(a$params$coords, "location index")
  expect_true(all(c("DI", "DI_location", "AOA") %in% names(a$aoa)))
  expect_true(is.list(a$location))
  expect_gt(a$location$threshold, 0)
  expect_length(a$location$train_DI, 400L)
  expect_identical(a$location$n_outside, 300L)
  expect_identical(a$n_new, a$n_inside + a$n_outside + a$n_na)
  expect_identical(a$aoa$AOA,
                   a$aoa$DI <= a$threshold &
                     a$aoa$DI_location <= a$location$threshold)
  # The cut-off as a distance is of the order of the spacing of 400 points
  # on a 1000 x 1000 square, not of the extent.
  cut_m <- a$location$threshold * a$location$normalizer
  expect_gt(cut_m, 20); expect_lt(cut_m, 200)
  # The geometry handed back is the caller's.
  expect_false(any(c("..x", "..y") %in% names(a$aoa)))

  out <- utils::capture.output(print(a))
  expect_true(any(grepl("location    : measured on its own", out, fixed = TRUE)))
  expect_true(any(grepl("on location", out, fixed = TRUE)))

  # coords = "no" is the covariates alone, as before, with a logged caution.
  lines <- capture_spatialkit_log(
    off <- area_of_applicability(loc_pts(300, 1500, 2500, seed = 30), model = gf,
                                 coords = "no"))
  expect_true(log_has(lines, "coords = \"no\", so location is not measured"))
  expect_null(off$location)
  expect_false("DI_location" %in% names(off$aoa))
  expect_gt(mean(off$aoa$AOA), 0.7)
  expect_identical(off$params$coords, "none")
})

test_that("an unsampled gap inside the bounding box is outside for a GWR", {
  skip_if_not_installed("sp")
  set.seed(1)
  f  <- function(x) 3 + 2 * sin(2 * pi * x / 800)
  xx <- c(stats::runif(200, 0, 300), stats::runif(200, 700, 1000))
  aa <- stats::rnorm(400)
  trn <- sf::st_as_sf(
    data.frame(x = xx, y = stats::runif(400, 0, 1000), a = aa,
               z = f(xx) * aa + stats::rnorm(400, 0, 0.3)),
    coords = c("x", "y"), crs = 32632)
  gf <- suppressWarnings(fit_gwr_model(trn, "z", "a", bandwidth = 40))
  gap <- sf::st_as_sf(
    data.frame(x = stats::runif(200, 400, 600), y = stats::runif(200, 0, 1000),
               a = stats::rnorm(200)), coords = c("x", "y"), crs = 32632)
  a <- area_of_applicability(gap, model = gf)
  expect_lt(mean(a$aoa$AOA), 0.1)

  # The documented surface workflow: the grid takes its covariates from the
  # nearest observation, so every cell's covariate DI is 0 and every cell
  # used to be inside, the gap included.
  surf <- suppressWarnings(suppressMessages(
    predict_surface(gf, n_cells = 900, covariates = trn)))
  s <- area_of_applicability(surf, model = gf)
  gx <- sf::st_coordinates(suppressWarnings(sf::st_centroid(sf::st_geometry(s$aoa))))[, 1]
  expect_true(all(s$aoa$DI == 0))
  in_gap <- gx > 400 & gx < 600
  expect_gt(sum(in_gap), 50L)
  expect_lt(mean(s$aoa$AOA[in_gap]), 0.1)
  expect_gt(mean(s$aoa$AOA[gx < 250 | gx > 750]), 0.8)
})

test_that("the location index is not diluted by the number of covariates", {
  for (p in c(1L, 12L)) {
    trn <- loc_pts(400, 0, 1000, p = p, seed = 100 * p + 1)
    pv  <- paste0("a", seq_len(p))
    inside <- function(x0, x1)
      mean(area_of_applicability(loc_pts(300, x0, x1, p = p), train_sf = trn,
                                 predictor_vars = pv, coords = "yes")$aoa$AOA)
    expect_gt(inside(0, 1000), 0.8)
    # As two more equally weighted columns: 82-89% inside at 1200-1500 and
    # 7-15% at 1500-2500 with 8 to 12 covariates.
    expect_lt(inside(1200, 1500), 0.05)
    expect_lt(inside(1500, 2500), 0.05)
    expect_lt(inside(4000, 5000), 0.05)
  }
})

test_that("a bayesian_fit and a gwr_fit are location models under coords = \"auto\"", {
  trn <- loc_pts(200, 0, 1000, seed = 5)
  far <- loc_pts(100, 3000, 4000, seed = 6)
  for (cls in c("gwr_fit", "bayesian_fit")) {
    fit <- loc_fake_fit(trn, "a1", cls)
    a <- area_of_applicability(far, model = fit)
    expect_identical(a$params$coords, "location index")
    expect_identical(a$n_inside, 0L)
    expect_identical(a$location$n_outside, 100L)
  }
  # Any other fit, and no fit at all, are measured on the covariates alone
  # unless asked.
  for (fit in list(loc_fake_fit(trn, "a1", "lmsurf_fit"), NULL)) {
    a <- if (is.null(fit))
      area_of_applicability(far, train_sf = trn, predictor_vars = "a1")
    else area_of_applicability(far, model = fit)
    expect_identical(a$params$coords, "none")
    expect_null(a$location)
    expect_false("DI_location" %in% names(a$aoa))
    expect_gt(a$n_inside, 50L)
  }
  yes <- area_of_applicability(far, model = loc_fake_fit(trn, "a1", "lmsurf_fit"),
                               coords = "yes")
  expect_identical(yes$n_inside, 0L)
})

test_that("the location index uses the folds and handles shared locations", {
  trn <- loc_pts(300, 0, 1000, seed = 7)
  nd  <- loc_pts(200, 0, 1300, seed = 8)
  fit <- loc_fake_fit(trn, "a1", "gwr_fit")
  fo  <- suppressWarnings(suppressMessages(
    make_folds(trn, k = 5, method = "block_kfold", seed = 1)))
  loo <- area_of_applicability(nd, model = fit)
  blk <- suppressWarnings(area_of_applicability(nd, model = fit, folds = fo))
  # A blocked hold-out is further from its training points than a
  # leave-one-out neighbour is, in map space as in predictor space.
  expect_gt(blk$location$threshold, loo$location$threshold)
  expect_identical(blk$aoa$DI_location, loo$aoa$DI_location)
  expect_gte(blk$n_inside, loo$n_inside)

  # Four visits to each of 75 sites: without folds a site's other visits are
  # the same location, not its neighbours.
  site <- loc_pts(75, 0, 1000, seed = 9)
  rep4 <- site[rep(1:75, each = 4), ]
  rep4$a1 <- stats::rnorm(300)
  fit4 <- loc_fake_fit(rep4, "a1", "gwr_fit")
  a4 <- expect_no_warning(area_of_applicability(nd, model = fit4))
  expect_gt(a4$location$threshold, 0)
  expect_identical(a4$location$n_duplicates, 225L)
  expect_identical(a4$n_duplicates, 0L)
  # Folds that split the visits: the location threshold is 0, and that is
  # the classed warning.
  set.seed(10)
  lab <- sample(rep(1:5, length.out = 300))
  expect_warning(
    z <- area_of_applicability(nd, model = fit4, folds = lab),
    class = "spatialkit_aoa_zero_threshold")
  expect_identical(z$location$threshold, 0)
})

test_that("the location index needs geometry, and a projected CRS", {
  trn <- loc_pts(100, 0, 1000, seed = 11)
  nd  <- loc_pts(20, 0, 1000, seed = 12)
  fit <- loc_fake_fit(trn, "a1", "gwr_fit")
  expect_error(area_of_applicability(sf::st_drop_geometry(nd), model = fit),
               "location is measured.*`newdata` is not.*coords = \"no\"")
  expect_s3_class(
    area_of_applicability(sf::st_drop_geometry(nd), model = fit, coords = "no"),
    "aoa")
  ll_tr <- sf::st_transform(trn, 4326); ll_nd <- sf::st_transform(nd, 4326)
  expect_error(area_of_applicability(ll_nd, train_sf = ll_tr,
                                     predictor_vars = "a1", coords = "yes"),
               "geographic CRS")
  # newdata in another CRS is brought into the training one for the
  # distances, and handed back as it came.
  a_ll <- area_of_applicability(ll_nd, model = fit)
  a_pr <- area_of_applicability(nd, model = fit)
  expect_equal(a_ll$aoa$DI_location, a_pr$aoa$DI_location, tolerance = 1e-6)
  expect_identical(sf::st_crs(a_ll$aoa), sf::st_crs(ll_nd))
  # A row with no geometry cannot be placed: inside on the covariates, not
  # classified overall.
  nd2 <- nd
  sf::st_geometry(nd2)[1] <- sf::st_point()
  a2 <- area_of_applicability(nd2, model = fit)
  expect_true(is.na(a2$aoa$DI_location[1]))
  expect_false(is.na(a2$aoa$DI[1]))
  expect_identical(a2$n_new, a2$n_inside + a2$n_outside + a2$n_na)
  expect_error(area_of_applicability(nd, model = fit, coords = "sometimes"),
               "should be one of")
})

test_that("a forest with the coordinates as predictors keeps its own path", {
  skip_if_not_installed("ranger")
  trn <- loc_pts(200, 0, 1000, seed = 13)
  far <- loc_pts(100, 3000, 4000, seed = 14)
  rf  <- suppressWarnings(suppressMessages(
    fit_rf_model(trn, "z", "a1", num_trees = 100, include_coords = TRUE)))
  a <- area_of_applicability(far, model = rf)
  expect_identical(a$params$coords, "predictors")
  expect_identical(a$predictor_vars, c("a1", "..x", "..y"))
  expect_null(a$location)
  expect_false("DI_location" %in% names(a$aoa))
  expect_identical(a$n_inside, 0L)
  expect_identical(area_of_applicability(far, model = rf, coords = "yes")$aoa$DI,
                   a$aoa$DI)
  # coords = "no" turns it off for a forest too, and says so.
  lines <- capture_spatialkit_log(
    off <- area_of_applicability(far, model = rf, coords = "no"))
  expect_true(log_has(lines, "coords = \"no\", so location is not measured"))
  expect_identical(off$predictor_vars, "a1")
})

test_that("with blocked folds a gap narrower than the hold-out distance reads as inside", {
  # Two sampled strips, x 0-300 and 700-1000, and nothing between them.
  # Without folds the location cut-off is the spacing of the points, and the
  # middle of the gap is outside.  With blocked folds it is the distance
  # from a held-out block to its training data, which here is wider than
  # half the gap, so the same points read as inside: the cross-validated
  # error was earned on hold-outs that far from their training rows.  The
  # help says so; it used to say the index marks unsampled gaps, without
  # condition.
  set.seed(1)
  xx <- c(stats::runif(200, 0, 300), stats::runif(200, 700, 1000))
  trn <- sf::st_as_sf(
    data.frame(x = xx, y = stats::runif(400, 0, 1000), a1 = stats::rnorm(400),
               z = stats::rnorm(400)), coords = c("x", "y"), crs = 32632)
  gap <- sf::st_as_sf(
    data.frame(x = stats::runif(200, 400, 600), y = stats::runif(200, 0, 1000),
               a1 = stats::rnorm(200)), coords = c("x", "y"), crs = 32632)
  fit <- loc_fake_fit(trn, "a1", "gwr_fit")
  fo  <- suppressWarnings(suppressMessages(
    make_folds(trn, k = 5, method = "block_kfold", seed = 1)))
  loo <- area_of_applicability(gap, model = fit)
  blk <- suppressWarnings(area_of_applicability(gap, model = fit, folds = fo))
  cut <- function(a) a$location$threshold * a$location$normalizer
  inside_loc <- function(a) mean(a$aoa$DI_location <= a$location$threshold)
  # Every gap point is at least 100 units from the nearest strip.
  expect_lt(cut(loo), 100)
  expect_identical(inside_loc(loo), 0)
  expect_gt(cut(blk), 100)
  expect_gt(inside_loc(blk), 0.5)
  # The distances themselves do not depend on the folds; the cut-off does.
  expect_identical(blk$aoa$DI_location, loo$aoa$DI_location)
})
