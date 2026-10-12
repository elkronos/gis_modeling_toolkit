# tests/testthat/test-fix-cells-kriging.R
# ---------------------------------------------------------------------------
# kriging_adequacy(): a response that is not numeric, replicates that show
# more variation than the model's nugget, and what kriging buys in a cell,
# measured against the plain mean on the same variogram (plain_var, kr_gain,
# kr_shift).
# ---------------------------------------------------------------------------

.fg_field <- function(n = 200, seed = 42) {
  set.seed(seed)
  x <- 5e5 + runif(n, 0, 10000); y <- 5e6 + runif(n, 0, 10000)
  d <- as.matrix(stats::dist(cbind(x, y)))
  lat <- as.numeric(t(chol(exp(-d / 2500) + diag(0.2, n))) %*% rnorm(n))
  sf::st_as_sf(data.frame(x = x, y = y, lat = lat), coords = c("x", "y"), crs = 32632)
}
.fg_setup <- function() {
  pts <- .fg_field()
  bnd <- sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(pts)))
  cells <- create_grid_polygons(bnd, target_cells = 9, type = "square")
  list(pts = pts, cells = cells,
       asg = assign_features_to_polygons(pts, cells))
}
.fg_sac <- function(nugget = 0.2, psill = 1, range = 2500) structure(
  3 * range, class = c("sac_range", "numeric"),
  variogram_model = if (nugget > 0)
    gstat::vgm(psill = psill, model = "Exp", range = range, nugget = nugget)
  else gstat::vgm(psill = psill, model = "Exp", range = range),
  crs = sf::st_crs(32632))

.fg_run <- function(expr) {
  w <- character(0); cls <- list(); conds <- list()
  val <- withCallingHandlers(
    expr,
    warning = function(c) {
      w <<- c(w, conditionMessage(c)); cls[[length(cls) + 1L]] <<- class(c)
      conds[[length(conds) + 1L]] <<- c
      invokeRestart("muffleWarning")
    })
  list(value = val, warnings = w, classes = cls, conditions = conds)
}


# --- the response's type ------------------------------------------------------

test_that("a factor, character or date response is an error naming the column and its class", {
  skip_if_not_installed("gstat")
  s <- .fg_setup()
  asg <- s$asg
  lat <- asg$lat
  asg$cat3 <- cut(lat, stats::quantile(lat, c(0, 1/3, 2/3, 1)), include.lowest = TRUE,
                  labels = c("low", "mid", "high"))
  asg$ord  <- factor(as.character(asg$cat3), levels = c("low", "mid", "high"),
                     ordered = TRUE)
  asg$chr  <- as.character(asg$cat3)
  asg$num_chr <- as.character(round(lat, 3))
  asg$day  <- as.Date("2024-01-01") + round(50 * lat)

  # A factor used to be kriged on its level codes, and the answer changed with
  # the order of the levels, with nothing said.
  expect_error(kriging_adequacy(asg, "cat3", s$cells, k = 4),
               "response 'cat3' is a factor, and kriging needs a numeric or logical variable")
  expect_error(kriging_adequacy(asg, "cat3", s$cells, k = 4),
               "arbitrary relabelling of its levels")
  expect_error(kriging_adequacy(asg, "cat3", s$cells, k = 4),
               "pass an indicator column \\(cat3 == \"<level>\"\\)")
  expect_error(kriging_adequacy(asg, "ord", s$cells, k = 4),
               "response 'ord' is a factor")
  # Labels used to end in "fewer than 10 points with a finite response".
  expect_error(kriging_adequacy(asg, "chr", s$cells, k = 4),
               "response 'chr' is of class character")
  # Numbers held as text, and dates, were converted without a word.
  expect_error(kriging_adequacy(asg, "num_chr", s$cells, k = 4),
               "response 'num_chr' is of class character.*not converted for you")
  expect_error(kriging_adequacy(asg, "day", s$cells, k = 4),
               "response 'day' is of class Date")
  # With `sac` supplied estimate_sac_range() is never called, so the check
  # cannot be left to it.
  expect_error(kriging_adequacy(asg, "cat3", s$cells, k = 4, sac = .fg_sac()),
               "response 'cat3' is a factor")
})

test_that("a logical response runs and its mean is the cell share; numeric is unchanged", {
  skip_if_not_installed("gstat")
  s <- .fg_setup()
  asg <- s$asg
  asg$pos <- asg$lat > 0
  ka <- suppressWarnings(kriging_adequacy(asg, "pos", s$cells, k = 4, sac = .fg_sac()))
  d <- sf::st_drop_geometry(asg)
  share <- as.numeric(tapply(d$pos, factor(d$poly_id, levels = ka$poly_id), mean))
  expect_equal(ka$mean, share, tolerance = 1e-12)
  # (a cell with no point has no mean)
  expect_true(all(ka$mean >= 0 & ka$mean <= 1, na.rm = TRUE))
  expect_gt(sum(is.finite(ka$mean)), 5L)
  # An integer column is numeric.
  asg$cnt <- as.integer(round(10 * asg$lat))
  expect_s3_class(suppressWarnings(kriging_adequacy(asg, "cnt", s$cells, k = 4,
                                                    sac = .fg_sac())), "kriging_adequacy")
  # A numeric response is kriged as before.
  kn <- suppressWarnings(kriging_adequacy(asg, "lat", s$cells, k = 4, sac = .fg_sac()))
  expect_equal(kn$mean, as.numeric(tapply(d$lat, factor(d$poly_id, levels = kn$poly_id),
                                          mean)), tolerance = 1e-12)
})


# --- replicates against the nugget -------------------------------------------

# `n_loc` stations, each visited `m` times: a field with no nugget of its own
# (partial sill 1, range parameter 2500), and a measurement error of variance
# `me_var` on every visit.
.fg_visits <- function(n_loc = 200, m = 2, me_var = 0.5, seed = 31) {
  set.seed(seed)
  x <- 5e5 + runif(n_loc, 0, 10000); y <- 5e6 + runif(n_loc, 0, 10000)
  d <- as.matrix(stats::dist(cbind(x, y)))
  s <- as.numeric(t(chol(exp(-d / 2500) + diag(1e-8, n_loc))) %*% rnorm(n_loc))
  v <- data.frame(x = rep(x, m), y = rep(y, m), station = rep(seq_len(n_loc), m))
  v$lat <- s[v$station] + rnorm(nrow(v), 0, sqrt(me_var))
  sf::st_as_sf(v, coords = c("x", "y"), crs = 32632)
}

test_that("a replicate variance above the nugget is used as the measurement error, and said", {
  skip_if_not_installed("gstat")
  v <- .fg_visits()
  bnd <- sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(v)))
  cells <- create_grid_polygons(bnd, target_cells = 9, type = "square")
  asg <- assign_features_to_polygons(v, cells)
  d <- sf::st_drop_geometry(asg)
  loc_mean <- as.numeric(tapply(d$lat, d$station, mean))
  sw2 <- sum((d$lat - loc_mean[d$station])^2) / 200          # pooled, 200 d.f.
  lab <- rep(1:4, length.out = 200)[d$station]               # a station is held out whole

  # A model with no nugget, and replicates that differ: the replicate variance
  # used to be capped at the nugget (here 0), so the station means were
  # kriged as exact and "(0 of 0)" was all that was said.
  res <- .fg_run(kriging_adequacy(asg, "lat", cells, sac = .fg_sac(nugget = 0), folds = lab))
  hit <- vapply(res$classes, function(k) "spatialkit_replicate_var_above_nugget" %in% k,
                logical(1))
  expect_identical(sum(hit), 1L)
  expect_match(res$warnings[hit], "400 point\\(s\\) share a location with another")
  expect_match(res$warnings[hit], sprintf("variance is %.3g (200 degrees of freedom)", sw2),
               fixed = TRUE)
  expect_match(res$warnings[hit], "above the variogram's nugget of 0")
  expect_match(res$warnings[hit], "is used as the measurement-error variance")
  expect_match(res$warnings[hit],
               sprintf("Refit the variogram with a nugget of at least %.3g", sw2), fixed = TRUE)
  expect_equal(res$conditions[[which(hit)]]$within_var, sw2, tolerance = 1e-10)
  expect_identical(res$conditions[[which(hit)]]$nugget, 0)
  rp <- attr(res$value, "replicates")
  expect_identical(rp$n_points, 400L)
  expect_equal(rp$within_var, sw2, tolerance = 1e-10)
  expect_identical(rp$df, 200L)
  expect_identical(rp$nugget, 0)
  expect_equal(rp$nugget_used, sw2, tolerance = 1e-10)
  expect_true(rp$above_nugget)
  expect_output(print(res$value),
                "replicates: pooled within-location variance .* is above the nugget 0 and is used")

  # It is used, not only reported: kriging the 200 means with error variance
  # sw2 / 2 is kriging all 400 visits with sw2 as their measurement error.
  vm0 <- attr(.fg_sac(nugget = 0), "variogram_model")
  bk <- gstat::krige(lat ~ 1, asg, cells, model = vm0, weights = rep(1 / sw2, nrow(asg)),
                     nmax = 1000, debug.level = 0)
  full <- suppressWarnings(kriging_adequacy(asg, "lat", cells, sac = .fg_sac(nugget = 0),
                                            folds = lab, nmax = 1000))
  expect_equal(full$kr_pred, as.numeric(bk$var1.pred), tolerance = 1e-6)
  expect_equal(full$kr_var, as.numeric(bk$var1.var), tolerance = 1e-6)

  # And the cross-validation statistic is nearer 1 for it.  Kriging the
  # station means as exact under the same model is what the cap did: over
  # ten layers like this one the statistic read 2.96-4.42 that way and
  # 0.84-1.15 with the replicate variance used.
  means <- asg[!duplicated(asg$station), ]
  means$lat <- loc_mean[means$station]
  exact <- suppressWarnings(kriging_adequacy(means, "lat", cells, sac = .fg_sac(nugget = 0),
                                             folds = lab[!duplicated(d$station)]))
  expect_null(attr(exact, "replicates"))
  z_new <- attr(res$value, "cv")$zscore_var
  z_cap <- attr(exact, "cv")$zscore_var
  expect_gt(z_cap, 2)
  expect_lt(abs(z_new - 1), abs(z_cap - 1))
  expect_gt(z_new, 0.75); expect_lt(z_new, 1.3)
  # The kriging variances grow with it: they were too small.
  expect_true(all(res$value$kr_var > exact$kr_var))
})

test_that("the replicate variance is used the same way on both sides of the nugget", {
  skip_if_not_installed("gstat")
  v <- .fg_visits()
  bnd <- sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(v)))
  cells <- create_grid_polygons(bnd, target_cells = 9, type = "square")
  asg <- assign_features_to_polygons(v, cells)
  d <- sf::st_drop_geometry(asg)
  loc_mean <- as.numeric(tapply(d$lat, d$station, mean))
  sw2 <- sum((d$lat - loc_mean[d$station])^2) / 200
  lab <- rep(1:4, length.out = 200)[d$station]

  # A nugget with room for it: the plain warning, no class, no print line,
  # and the rest of the nugget stays with the station as micro-scale variation.
  ok <- .fg_run(kriging_adequacy(asg, "lat", cells, sac = .fg_sac(nugget = 2 * sw2), folds = lab))
  shared <- grepl("share a location", ok$warnings)
  expect_identical(sum(shared), 1L)
  expect_false("spatialkit_replicate_var_above_nugget" %in% ok$classes[[which(shared)]])
  expect_match(ok$warnings[shared], "differs between replicates")
  rp <- attr(ok$value, "replicates")
  expect_false(rp$above_nugget)
  expect_equal(rp$nugget_used, 2 * sw2, tolerance = 1e-10)
  expect_false(any(grepl("replicates:", utils::capture.output(print(ok$value)))))

  # A nugget a thousandth below the replicate variance and one a thousandth
  # above it give the same kriging to that order: nothing jumps where the
  # replicate variance crosses the nugget.
  below <- suppressWarnings(kriging_adequacy(asg, "lat", cells, folds = lab,
                                             sac = .fg_sac(nugget = sw2 * 0.999)))
  above <- suppressWarnings(kriging_adequacy(asg, "lat", cells, folds = lab,
                                             sac = .fg_sac(nugget = sw2 * 1.001)))
  expect_true(attr(below, "replicates")$above_nugget)
  expect_false(attr(above, "replicates")$above_nugget)
  expect_equal(below$kr_var, above$kr_var, tolerance = 5e-3)
  expect_equal(below$plain_var, above$plain_var, tolerance = 5e-3)
  expect_equal(attr(below, "cv")$zscore_var, attr(above, "cv")$zscore_var, tolerance = 5e-3)

  # No replicates: no record.
  one <- asg[!duplicated(asg$station), ]
  none <- suppressWarnings(kriging_adequacy(one, "lat", cells, sac = .fg_sac(),
                                            folds = lab[!duplicated(d$station)]))
  expect_null(attr(none, "replicates"))
})


# --- kriging against the plain mean, on one variogram ------------------------

# The plain mean's error variance against the cell's own average, written out
# with whole matrices: Var(sum a_j z_j - mean of the field over the cell's
# 500 discretisation points).  `e` is each location's error variance.
.fg_plain_var <- function(xy, a, e, cell, range, psill) {
  geo <- sf::as_Spatial(sf::st_set_crs(sf::st_geometry(cell), sf::NA_crs_))
  G <- sp::coordinates(sp::spsample(geo, n = 500, type = "regular", offset = c(0.5, 0.5)))
  cross <- function(A, B) sqrt(outer(A[, 1], B[, 1], "-")^2 + outer(A[, 2], B[, 2], "-")^2)
  Cs <- function(h) psill * exp(-h / range)
  as.numeric(t(a) %*% Cs(cross(xy, xy)) %*% a) + sum(a^2 * e) -
    2 * sum(a * rowMeans(Cs(cross(xy, G)))) + mean(Cs(cross(G, G)))
}

test_that("plain_var is the plain mean's variance under the variogram, and kr_var never exceeds it", {
  skip_if_not_installed("gstat")
  s <- .fg_setup()
  ka <- suppressWarnings(kriging_adequacy(s$asg, "lat", s$cells, sac = .fg_sac(), k = 4))
  df <- sf::st_drop_geometry(ka)
  # The comparison with s^2/n is gone; the two columns that replace it are in.
  expect_false("kr_exceeds_design" %in% names(df))
  expect_true(all(c("plain_var", "kr_gain", "kr_shift") %in% names(df)))

  pop <- df$n >= 1L
  expect_gt(sum(pop), 8L)
  expect_gt(sum(!pop), 0L)                    # this grid has empty cells
  expect_true(any(df$n == 1L))                # and a cell with one point
  expect_true(all(is.finite(df$plain_var[pop]) & df$plain_var[pop] > 0))
  # Ordinary kriging is the best unbiased combination of a neighbourhood
  # that holds the cell's own points, and the plain mean is one of them.
  expect_true(all(df$kr_var[pop] <= df$plain_var[pop] * (1 + 1e-6)))
  expect_equal(df$kr_gain[pop], 1 - df$kr_var[pop] / df$plain_var[pop], tolerance = 1e-12)
  expect_true(all(df$kr_gain[pop] > 0 & df$kr_gain[pop] < 1))
  # Defined for a cell with one point, where s^2/n is not.
  one <- df$n == 1L
  expect_true(all(is.na(df$se[one])))
  expect_true(all(is.finite(df$kr_gain[one]) & is.finite(df$kr_shift[one])))
  # A cell with no point has no plain mean to compare with.
  expect_true(all(is.na(df$plain_var[!pop]) & is.na(df$kr_gain[!pop]) &
                    is.na(df$kr_shift[!pop])))
  expect_true(all(is.finite(df$kr_pred[!pop])))

  # plain_var against the definition written out, for three cells.
  d <- sf::st_drop_geometry(s$asg); xy <- sf::st_coordinates(s$asg)[, 1:2]
  for (i in utils::head(which(df$n >= 2L), 2L)) {
    own <- d$poly_id == df$poly_id[i]
    expect_equal(df$plain_var[i],
                 .fg_plain_var(xy[own, , drop = FALSE], rep(1 / sum(own), sum(own)),
                               rep(0.2, sum(own)), s$cells[i, ], range = 2500, psill = 1),
                 tolerance = 1e-8)
  }
  i <- which(one)[1L]; own <- d$poly_id == df$poly_id[i]
  expect_equal(df$plain_var[i], .fg_plain_var(xy[own, , drop = FALSE], 1, 0.2, s$cells[i, ],
                                              range = 2500, psill = 1), tolerance = 1e-8)

  # Cells whose neighbourhood had to grow to hold all their own points.
  small <- suppressWarnings(kriging_adequacy(s$asg, "lat", s$cells, sac = .fg_sac(), k = 4,
                                             nmax = 10))
  sm <- sf::st_drop_geometry(small)
  grown <- !is.na(sm$kr_n_used) & sm$kr_n_used > 10L
  expect_gt(sum(grown), 3L)
  expect_true(all(sm$kr_var[grown] <= sm$plain_var[grown] * (1 + 1e-6)))
  expect_true(all(sm$kr_gain[sm$n >= 1L] > 0))
  # plain_var is a property of the cell's own points: the neighbourhood does
  # not enter it.
  expect_equal(sm$plain_var, df$plain_var, tolerance = 1e-12)
})

test_that("kr_gain does not move with the scale of the sill, and kr_shift is in the model's units", {
  skip_if_not_installed("gstat")
  s <- .fg_setup()
  base <- suppressWarnings(kriging_adequacy(s$asg, "lat", s$cells, sac = .fg_sac(), k = 4))
  # The same data, the same kriged estimates, sill and nugget multiplied by f.
  # The comparison with s^2/n flagged 0, 1 and 8 of these 9 cells at f = 1/3,
  # 1 and 3.
  for (f in c(1 / 3, 3)) {
    sc <- suppressWarnings(kriging_adequacy(s$asg, "lat", s$cells, k = 4,
                                            sac = .fg_sac(nugget = 0.2 * f, psill = f)))
    expect_equal(sc$kr_pred, base$kr_pred, tolerance = 1e-8)
    expect_equal(sc$kr_gain, base$kr_gain, tolerance = 1e-8)
    expect_equal(sc$kr_ratio, base$kr_ratio, tolerance = 1e-8)
    expect_equal(sc$plain_var, f * base$plain_var, tolerance = 1e-8)
    expect_equal(sc$kr_var, f * base$kr_var, tolerance = 1e-6)
    # A sill three times too large makes every shift look sqrt(3) smaller:
    # the shift is a check of the variogram, like the CV statistic.
    expect_equal(sc$kr_shift * sqrt(f), base$kr_shift, tolerance = 1e-5)
  }
})

test_that("kr_shift is the kriged-minus-plain difference over its own model standard deviation", {
  skip_if_not_installed("gstat")
  s <- .fg_setup()
  # Every point in every neighbourhood, so the kriging can be written out.
  ka <- suppressWarnings(kriging_adequacy(s$asg, "lat", s$cells, sac = .fg_sac(), k = 4,
                                          nmax = 1000))
  df <- sf::st_drop_geometry(ka)
  pop <- df$n >= 1L
  expect_equal(df$kr_shift[pop],
               (df$kr_pred[pop] - df$mean[pop]) / sqrt(df$plain_var[pop] - df$kr_var[pop]),
               tolerance = 1e-12)
  # The variance of (kriged - plain) from the kriging weights themselves is
  # plain_var - kr_var: the kriging error is uncorrelated with the difference.
  d <- sf::st_drop_geometry(s$asg); xy <- sf::st_coordinates(s$asg)[, 1:2]; n <- nrow(xy)
  cross <- function(A, B) sqrt(outer(A[, 1], B[, 1], "-")^2 + outer(A[, 2], B[, 2], "-")^2)
  K <- exp(-cross(xy, xy) / 2500) + diag(0.2, n)
  geo <- sf::as_Spatial(sf::st_set_crs(sf::st_geometry(s$cells), sf::NA_crs_))
  for (i in utils::head(which(pop), 3L)) {
    G <- sp::coordinates(sp::spsample(geo[i], n = 500, type = "regular", offset = c(0.5, 0.5)))
    k <- rowMeans(exp(-cross(xy, G) / 2500))
    lam <- solve(rbind(cbind(K, 1), c(rep(1, n), 0)), c(k, 1))[seq_len(n)]
    a <- as.numeric(d$poly_id == df$poly_id[i]); a <- a / sum(a)
    expect_equal(sum(lam * d$lat), df$kr_pred[i], tolerance = 1e-6)
    expect_equal(as.numeric(t(lam - a) %*% K %*% (lam - a)),
                 df$plain_var[i] - df$kr_var[i], tolerance = 1e-4)
  }
})

test_that("where kriging returns the plain mean the gain is 0 and the shift is not defined", {
  skip_if_not_installed("gstat")
  s <- .fg_setup()
  # One cell holding every point, under a pure nugget: the best unbiased
  # combination of exchangeable observations is their plain mean, so the two
  # estimates and the two variances coincide, and their difference over its
  # standard deviation is rounding error over rounding error.
  one <- sf::st_sf(poly_id = 1L, geometry = sf::st_as_sfc(sf::st_bbox(s$pts)))
  asg <- s$asg; asg$poly_id <- 1L
  nug <- structure(7500, class = c("sac_range", "numeric"),
                   variogram_model = gstat::vgm(psill = 0.7, model = "Nug", range = 0),
                   crs = sf::st_crs(32632))
  ka <- suppressWarnings(kriging_adequacy(asg, "lat", one, sac = nug, k = 4, nmax = 1000))
  expect_equal(ka$kr_pred, ka$mean, tolerance = 1e-8)
  expect_equal(ka$plain_var, 0.7 / 200, tolerance = 1e-12)
  expect_equal(ka$kr_var, ka$plain_var, tolerance = 1e-5)
  expect_gte(ka$kr_gain, 0)
  expect_lt(ka$kr_gain, 1e-5)
  expect_true(is.na(ka$kr_shift))
  expect_output(print(ka), "kr_gain")
  expect_false(any(grepl("kriged minus plain mean", utils::capture.output(print(ka)))))
})

test_that("with repeat visits the plain mean's variance uses the error variances the kriging uses", {
  skip_if_not_installed("gstat")
  v <- .fg_visits(n_loc = 90, m = 3, me_var = 0.1, seed = 8)
  bnd <- sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(v)))
  cells <- create_grid_polygons(bnd, target_cells = 9, type = "square")
  asg <- assign_features_to_polygons(v, cells)
  d <- sf::st_drop_geometry(asg)
  loc_mean <- as.numeric(tapply(d$lat, d$station, mean))
  sw2 <- sum((d$lat - loc_mean[d$station])^2) / 180
  first <- !duplicated(d$station); xy <- sf::st_coordinates(asg)[first, 1:2]
  cell_of <- d$poly_id[first]
  for (nug in c(0, 0.4)) {                     # the replicates above, and below, the nugget
    ka <- suppressWarnings(kriging_adequacy(asg, "lat", cells, sac = .fg_sac(nugget = nug),
                                            k = 4))
    df <- sf::st_drop_geometry(ka)
    pop <- df$n >= 1L
    expect_true(all(df$kr_var[pop] <= df$plain_var[pop] * (1 + 1e-6)), info = nug)
    # A station's three visits share what the nugget holds beyond the
    # replicate variance and average the rest.
    e <- max(nug - sw2, 0) + sw2 / 3
    for (i in utils::head(which(df$n >= 6L), 2L)) {
      own <- cell_of == df$poly_id[i]
      expect_equal(df$plain_var[i],
                   .fg_plain_var(xy[own, , drop = FALSE], rep(1 / sum(own), sum(own)),
                                 rep(e, sum(own)), cells[i, ], range = 2500, psill = 1),
                   tolerance = 1e-8, info = nug)
    }
  }
})


# --- what print() says of the comparison -------------------------------------

test_that("print() reports the gain, the shift on the model's scale and the naive SE against plain_var", {
  skip_if_not_installed("gstat")
  s <- .fg_setup()
  ka <- suppressWarnings(kriging_adequacy(s$asg, "lat", s$cells, sac = .fg_sac(), k = 4))
  df <- sf::st_drop_geometry(ka)
  out <- utils::capture.output(print(ka))
  # The line that read as a count of cells where kriging does not pay is gone.
  expect_false(any(grepl("kriging variance above s^2/n", out, fixed = TRUE)))
  gain <- out[grepl("kr_gain", out, fixed = TRUE)]
  expect_length(gain, 1L)
  expect_match(gain, "share of the plain mean's error variance that kriging removes", fixed = TRUE)
  expect_match(gain, sprintf("median %.2f, range %.2f-%.2f over the %d cell(s) with points and a kriged estimate",
                             stats::median(df$kr_gain, na.rm = TRUE),
                             min(df$kr_gain, na.rm = TRUE), max(df$kr_gain, na.rm = TRUE),
                             sum(is.finite(df$kr_gain))), fixed = TRUE)
  sh <- out[grepl("kriged minus plain mean", out, fixed = TRUE)]
  expect_match(sh, "sqrt(plain_var - kr_var)", fixed = TRUE)
  expect_match(sh, sprintf("|shift| > 1 in %d of %d cell(s), > 2 in %d",
                           sum(abs(df$kr_shift) > 1, na.rm = TRUE), sum(is.finite(df$kr_shift)),
                           sum(abs(df$kr_shift) > 2, na.rm = TRUE)), fixed = TRUE)
  # The naive standard error is set against the same model, not against kr_var.
  two <- is.finite(df$se)
  nv <- out[grepl("naive s^2/n", out, fixed = TRUE)]
  expect_length(nv, 1L)
  expect_match(nv, sprintf("in %d of the %d cell(s) with two or more points",
                           sum(df$se[two]^2 < df$plain_var[two]), sum(two)), fixed = TRUE)
})
