# tests/testthat/test-review3-S4-range-kriging.R
# ---------------------------------------------------------------------------
# Regressions from the third review of estimate_sac_range(), its print() and
# plot() methods, and kriging_adequacy().  Each test names the finding it
# closes.
# ---------------------------------------------------------------------------

# n points on a 1000 m square in EPSG:32632, a white-noise covariate `cov`,
# and z = 0.5 * cov plus an exponential field with range parameter `rp`
# (effective range 3 * rp) and a 0.2 nugget; rp = NULL gives white noise.
.r3s4_field <- function(seed, n = 300, rp = 10) {
  set.seed(seed)
  d <- data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000))
  if (is.null(rp)) {
    d$z <- rnorm(n); d$cov <- rnorm(n)
  } else {
    D <- as.matrix(stats::dist(d)); d$cov <- rnorm(n)
    d$z <- 0.5 * d$cov + as.numeric(t(chol(exp(-D / rp) + diag(0.2, n))) %*% rnorm(n))
  }
  sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
}

# The white-noise REML refusal of test-sac-range.R (seed 7, 0.27 m), fitted
# once for the tests that read it.
.r3s4_wn <- local({
  val <- NULL
  function() {
    if (is.null(val)) {
      set.seed(7)
      n <- 300
      d <- data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
                      z = rnorm(n), w = rnorm(n))
      pts <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
      lines <- capture_spatialkit_log(
        r <- suppressWarnings(estimate_sac_range(pts, "z", "w", detrend = "reml")))
      val <<- list(r = r, lines = lines, xy = as.matrix(d[, c("x", "y")]))
    }
    val
  }
})

.r3s4_warnings <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  })
  list(value = val, warnings = w)
}

.r3s4_first_lag <- function(r) {
  vg <- attr(r, "variogram")
  min(vg$dist[vg$np > 0])
}


# ---------------------------------------------------------------------------
# S4-RANGE-KRIGING-1: the REML floor is set by the point pairs, not the bins
# ---------------------------------------------------------------------------

test_that("a REML range the point pairs support is returned, and the same at any cutoff", {
  # A true effective range of 30 m on 300 points: the REML fit returned 21.0
  # m and it was refused as below the first bin of the diagnostic variogram
  # (30 m at cutoff = 0.5), while at cutoff = 0.1 (first bin 6 m) the same
  # REML number came back.
  skip_if_not_installed("gstat")
  skip_if_not_installed("nlme")
  pts <- .r3s4_field(9, n = 300, rp = 10)
  r5 <- estimate_sac_range(pts, "z", "cov", detrend = "reml")
  r1 <- estimate_sac_range(pts, "z", "cov", detrend = "reml", cutoff = 0.1)
  expect_true(is.finite(r5))
  expect_equal(as.numeric(r1), as.numeric(r5))
  # Shorter than the first bin, which alone used to refuse it ...
  expect_lt(as.numeric(r5), .r3s4_first_lag(r5))
  # ... and longer than the distance within which 30 pairs of its points lie.
  d30 <- sort(as.numeric(stats::dist(sf::st_coordinates(pts))))[30]
  expect_gt(as.numeric(r5), d30)
  expect_null(attr(r5, "range_floor"))
})

test_that("a REML range below 30 pairs is refused, says so, and keeps the REML numbers", {
  skip_if_not_installed("gstat")
  skip_if_not_installed("nlme")
  wn <- .r3s4_wn()
  r <- wn$r
  expect_true(is.na(r))
  expect_s3_class(r, "sac_range")
  expect_identical(attr(r, "rejected_reason"), "fitted range is below the shortest lag fitted")
  # The floor is the 30th-shortest distance between the 300 points the fit
  # used (all of them: no duplicates, n <= reml_max_n), below the first lag.
  d30 <- sort(as.numeric(stats::dist(wn$xy)))[30]
  expect_equal(attr(r, "range_floor"), d30)
  expect_lt(attr(r, "range_floor"), .r3s4_first_lag(r))
  expect_lt(attr(r, "rejected_range"), attr(r, "range_floor"))
  # The warning names that reference, not a white-noise verdict.
  expect_true(log_has(wn$lines, "fewer than 30 pairs of the 300 points the REML fit used"))
  expect_true(log_has(wn$lines, "within which 30 of its point pairs lie"))
  expect_false(log_has(wn$lines, "no spatial structure"))
  # The refusal carries the REML fit's own numbers, as the success does.
  rm <- attr(r, "reml")
  expect_type(rm, "list")
  expect_identical(rm$n_used, 300L)
  expect_false(rm$subsampled)
  expect_true(all(c("nugget_prop", "sigma2") %in% names(rm)))
})

test_that("the gstat-path floor refusal points at a smaller cutoff", {
  # Not reachable on ordinary data (the binned fit did not go below the first
  # lag on white noise or on 60 m fields), so the fit is stubbed to return a
  # 6 m effective range on a field whose variogram rises normally.
  skip_if_not_installed("gstat")
  pts <- .r3s4_field(2, n = 150, rp = 100)
  local_mocked_bindings(
    fit.variogram = function(object, model, ...) {
      m <- gstat::vgm(psill = 1, model = "Exp", range = 2, nugget = 0.2)
      attr(m, "singular") <- FALSE
      attr(m, "SSErr") <- 1
      m
    },
    .package = "gstat")
  lines <- capture_spatialkit_log(r <- estimate_sac_range(pts, "z"))
  expect_true(is.na(r))
  expect_identical(attr(r, "rejected_reason"), "fitted range is below the shortest lag fitted")
  expect_equal(attr(r, "range_floor"), .r3s4_first_lag(r))
  expect_null(attr(r, "reml"))
  expect_true(log_has(lines, "shorter than the shortest lag the variogram resolves"))
  expect_true(log_has(lines, "re-run with a smaller `cutoff`"))
})


# ---------------------------------------------------------------------------
# S4-RANGE-KRIGING-12: no "supply predictor_vars" for a detrended estimate
# ---------------------------------------------------------------------------

test_that("the past-the-lags advice fits whether the variogram was detrended", {
  skip_if_not_installed("gstat")
  set.seed(1); n <- 300
  d <- data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000))
  d$w <- rnorm(n); d$z <- (d$x - 5e5) / 200 + 0.3 * d$w + rnorm(n, sd = 0.3)
  pts <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
  lines <- capture_spatialkit_log(
    r <- suppressWarnings(estimate_sac_range(pts, "z", predictor_vars = "w")))
  skip_if(!identical(attr(r, "rejected_reason"),
                     "fitted range exceeds the largest lag fitted"),
          "the trend did not carry the fit past the lags on this platform")
  expect_true(log_has(lines, "try predictors that carry the trend"))
  expect_true(log_has(lines, "detrend = \"reml\""))
  expect_false(log_has(lines, "supply `predictor_vars`"))
  # Without predictors the advice to supply them stands.
  lines0 <- capture_spatialkit_log(
    r0 <- suppressWarnings(estimate_sac_range(pts, "z")))
  skip_if(!identical(attr(r0, "rejected_reason"),
                     "fitted range exceeds the largest lag fitted"),
          "the trend did not carry the fit past the lags on this platform")
  expect_true(log_has(lines0, "supply `predictor_vars` to detrend"))
})


# ---------------------------------------------------------------------------
# S4-RANGE-KRIGING-10: print() on a CRS-less estimate says what it is in
# ---------------------------------------------------------------------------

test_that("print() names the units of a CRS-less estimate and what was modelled", {
  skip_if_not_installed("gstat")
  pts <- .r3s4_field(3, n = 150, rp = 100)
  r <- suppressWarnings(estimate_sac_range(sf::st_set_crs(pts, NA), "z"))
  out <- utils::capture.output(print(r))
  expect_true(any(grepl("in the coordinate units of a layer with no CRS; variogram of the response itself",
                        out, fixed = TRUE)))
  # A hand-built object that records no CRS at all prints the number only.
  bare <- structure(300, class = c("sac_range", "numeric"))
  expect_identical(utils::capture.output(print(bare)), "300 ")
})


# ---------------------------------------------------------------------------
# S4-RANGE-KRIGING-11: plot() captions the floor refusal with its numbers
# ---------------------------------------------------------------------------

test_that("plot() captions the floor refusal with the refused range and the floor", {
  skip_if_not_installed("gstat")
  skip_if_not_installed("nlme")
  skip_if_not_installed("ggplot2")
  r <- .r3s4_wn()$r
  sub <- gsub("\n", " ", plot(r)$labels$subtitle, fixed = TRUE)
  expect_match(sub, sprintf("the fitted range (%.3g) is below %.3g", attr(r, "rejected_range"),
                            attr(r, "range_floor")), fixed = TRUE)
  expect_match(sub, "30 pairs", fixed = TRUE)
  # On the variogram path (or an object without the floor recorded) it is
  # the first lag that is named.
  gs <- r
  attr(gs, "range_floor") <- NULL
  attr(gs, "reml") <- NULL
  sub2 <- gsub("\n", " ", plot(gs)$labels$subtitle, fixed = TRUE)
  expect_match(sub2, sprintf("is below the shortest lag fitted (%.3g)", .r3s4_first_lag(r)),
               fixed = TRUE)
})


# ---------------------------------------------------------------------------
# kriging_adequacy(): S4-RANGE-KRIGING-8, -9 and S3-ASSIGNMENT-3
# ---------------------------------------------------------------------------

.r3s4_ka_setup <- function() {
  set.seed(1)
  n <- 200
  x <- 5e5 + runif(n, 0, 1000); y <- 5e6 + runif(n, 0, 1000)
  d <- as.matrix(stats::dist(cbind(x, y)))
  z <- as.numeric(t(chol(0.8 * exp(-d / 100) + diag(0.2, n))) %*% rnorm(n))
  pts <- sf::st_as_sf(data.frame(x = x, y = y, z = z), coords = c("x", "y"), crs = 32632)
  bnd <- sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(pts)))
  cells <- create_grid_polygons(bnd, target_cells = 16, type = "square")
  list(pts = pts, cells = cells, asg = assign_features_to_polygons(pts, cells))
}

test_that("kriging_adequacy() says why there is no variogram model, and whose", {
  skip_if_not_installed("gstat")
  s <- .r3s4_ka_setup()
  sing <- structure(NA_real_, class = c("sac_range", "numeric"),
                    rejected_reason = "no variogram model could be fitted (singular fits)")
  expect_error(kriging_adequacy(s$asg, "z", s$cells, sac = sing),
               "`sac` carries no variogram model.*no variogram model could be fitted \\(singular fits\\)")
  # No `sac` passed: the estimate made here is named, not the argument.
  flat <- s$asg; flat$z <- 1
  err <- tryCatch(suppressWarnings(kriging_adequacy(flat, "z", s$cells)),
                  error = conditionMessage)
  expect_match(err, "the variogram estimated here from 'z' carries no variogram model",
               fixed = TRUE)
  # The estimate's own reason, now that its early NA returns carry one.
  expect_match(err, "nothing to krige with: the response is constant", fixed = TRUE)
  expect_false(grepl("`sac`", err, fixed = TRUE))
})

test_that("kriging_adequacy() takes CRS-less cells, or points, to be in the other layer's CRS", {
  skip_if_not_installed("gstat")
  s <- .r3s4_ka_setup()
  # A variogram fitted in kilometres: the points go into that CRS, and the
  # cells must follow from the points' own metres, not be labelled km.
  sac_km <- estimate_sac_range(sf::st_transform(s$pts, "+proj=utm +zone=32 +units=km"), "z")
  base <- kriging_adequacy(s$asg, "z", s$cells, sac = sac_km, k = 4)
  w <- .r3s4_warnings(kriging_adequacy(s$asg, "z", sf::st_set_crs(s$cells, NA),
                                       sac = sac_km, k = 4))
  expect_true(any(grepl("`cells_sf` has no CRS", w$warnings, fixed = TRUE)))
  expect_equal(w$value$kr_var, base$kr_var)
  expect_equal(w$value$n, base$n)
  # The reverse: points without a CRS, cells with one.
  w2 <- .r3s4_warnings(kriging_adequacy(sf::st_set_crs(s$asg, NA), "z", s$cells,
                                        sac = sac_km, k = 4))
  expect_true(any(grepl("`assigned_points_sf` has no CRS", w2$warnings, fixed = TRUE)))
  expect_equal(w2$value$kr_var, base$kr_var)
})

test_that("kriging_adequacy() matches cell IDs as summarize_by_cell() does", {
  skip_if_not_installed("gstat")
  sq <- function(x0, y0, s) sf::st_polygon(list(rbind(c(x0, y0), c(x0 + s, y0),
                                                      c(x0 + s, y0 + s), c(x0, y0 + s), c(x0, y0))))
  grid <- sf::st_sf(poly_id = 1:16, geometry = sf::st_sfc(lapply(0:15, function(k)
    sq(5e5 + (k %% 4) * 100, 5e6 + (k %/% 4) * 100, 100)), crs = 32632))
  set.seed(11); n <- 200; x <- runif(n, 0, 400); y <- runif(n, 0, 400)
  z <- as.numeric(t(chol(exp(-as.matrix(stats::dist(cbind(x, y))) / 80) + diag(1e-6, n))) %*% rnorm(n))
  pts <- sf::st_as_sf(data.frame(x = 5e5 + x, y = 5e6 + y, z = z), coords = c("x", "y"), crs = 32632)
  sac <- structure(240, class = c("sac_range", "numeric"),
                   variogram_model = gstat::vgm(1, "Exp", 80, 0.05), crs = sf::st_crs(32632))

  # Integer cell IDs from 1e5 up against the assigned layer's IDs read back
  # as double: "1e+05" missed "100000", and those cells came back with n = 0.
  g2 <- grid; g2$poly_id <- g2$poly_id * 100000L
  a2 <- assign_features_to_polygons(pts, g2); a2$poly_id <- as.double(a2$poly_id)
  sm <- suppressWarnings(summarize_by_cell(a2, "z", cells_sf = g2))
  w <- .r3s4_warnings(kriging_adequacy(a2, "z", g2, sac = sac))
  expect_false(any(grepl("match no", w$warnings)))
  k2 <- w$value
  expect_equal(sum(k2$n), n)
  expect_equal(k2$n, sm$n[match(as.character(k2$poly_id), sm$poly_id)])
  expect_true(all(k2$n > 0L))

  # Cells keyed by 'id', which summarize_by_cell() joins.
  gid <- grid; names(gid)[names(gid) == "poly_id"] <- "id"
  a <- assign_features_to_polygons(pts, gid)
  kid <- kriging_adequacy(a, "z", gid, sac = sac)
  expect_s3_class(kid, "kriging_adequacy")
  expect_true("id" %in% names(kid))
  expect_equal(sum(kid$n), n)

  # A point whose ID matches no cell is said, not dropped silently.
  a3 <- assign_features_to_polygons(pts, grid)
  w3 <- .r3s4_warnings(kriging_adequacy(a3, "z", grid[-1, ], sac = sac))
  expect_true(any(grepl(sprintf("%d of the %d point\\(s\\) with a cell ID carry one that matches no",
                                sum(a3$poly_id == 1L), n), w3$warnings)))
  expect_equal(sum(w3$value$n), n - sum(a3$poly_id == 1L))

  # Neither layer has an ID column it knows: both lists are named.
  none <- grid; names(none)[names(none) == "poly_id"] <- "zone"
  expect_error(kriging_adequacy(a3, "z", none, sac = sac),
               "`cells_sf`. Looked for: poly_id, polygon_id, id, cell_id, grid_id", fixed = TRUE)
})
