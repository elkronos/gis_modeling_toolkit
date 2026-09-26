# tests/testthat/test-review-resolution.R
# ---------------------------------------------------------------------------
# Regressions from the review of the resolution step: resolution_profile()
# and determine_optimal_levels().
# ---------------------------------------------------------------------------

# A hand-made sac_range, so the variogram-based columns exist without gstat
# and without the time an estimate takes.
rr_sac <- function(pts, range = 1500, nugget = 1, psill = 1) {
  vm <- data.frame(model = c("Nug", "Exp"), psill = c(nugget, psill),
                   range = c(0, range / 3), stringsAsFactors = FALSE)
  structure(range, class = c("sac_range", "numeric"), variogram_model = vm,
            crs = sf::st_crs(pts))
}


test_that("one missing response or predictor value does not switch the profile to the raw response", {
  # The response is the predictor plus white noise, and the predictor carries
  # a strong east-west trend: the residuals have no structure, the raw
  # response a great deal.  One NA used to make the first lm.fit() on every
  # row fail, and the profile then scored the raw response.
  set.seed(10)
  n <- 300
  d <- data.frame(x = runif(n, 0, 10000), y = runif(n, 0, 10000))
  d$p <- d$x / 1000 + rnorm(n)
  d$z <- 3 * d$p + rnorm(n)
  pts <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
  sac <- rr_sac(pts)
  clean <- resolution_profile(pts, "z", "p", levels = c(10, 20, 40), sac = sac)
  expect_identical(attr(clean, "variable"), "residuals")
  for (col in c("z", "p")) {
    holed <- pts
    holed[[col]][17] <- NA
    lines <- capture_spatialkit_log(
      prof <- resolution_profile(holed, "z", "p", levels = c(10, 20, 40), sac = sac))
    expect_identical(attr(prof, "variable"), "residuals", info = col)
    expect_true(log_has(lines, "1 of 300 row"), info = col)
    expect_false(log_has(lines, "OLS fit on `predictor_vars` failed"), info = col)
    # One row fewer of the same residuals: the RSS barely moves.  The raw
    # response's RSS was several times larger.
    expect_equal(prof$rss, clean$rss, tolerance = 0.05, info = col)
    expect_equal(prof$cp, clean$cp, tolerance = 0.05, info = col)
  }
})


test_that("a layer larger than sample_n is bounded, judged and scored on all its points", {
  # 1200 points on a 10 km square, fitted on a 300-point subsample.  A range
  # of 1200 puts the floor at ceiling(1e8 / 1200^2) = 70 cells, which the
  # layer supports (floor(1200 / 9) = 133) and the subsample alone did not
  # (floor(300 / 9) = 33): the profile used to call these data unsupported,
  # run its ladder from 2 to 33, and hand a count of at most 33 to the
  # tessellation of all 1200 points.
  set.seed(3)
  N <- 1200
  d <- data.frame(x = runif(N, 0, 10000), y = runif(N, 0, 10000), z = rnorm(N))
  pts <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
  sac <- rr_sac(pts, range = 1200, nugget = 0.5, psill = 1)
  lines <- capture_spatialkit_log(
    prof <- resolution_profile(pts, "z", sac = sac, n_levels = 4, nstart = 3,
                               sample_n = 300))
  b <- attr(prof, "bounds")
  expect_identical(b$n, 1200L)
  expect_identical(b$n_sample, 300L)
  expect_identical(b$ceiling, 133L)
  expect_identical(b$ceiling_from, "min_cell_n")
  expect_true(b$supported)
  expect_false(log_has(lines, "cannot support"))
  expect_identical(min(prof$levels), b$floor)
  expect_identical(max(prof$levels), 133L)
  expect_output(print(prof), "on 1200 points")
  # The bounds do not move with sample_n.
  b2 <- attr(resolution_profile(pts, "z", sac = sac, n_levels = 4, nstart = 3,
                                sample_n = 600), "bounds")
  expect_identical(b2[c("floor", "ceiling", "supported", "n")],
                   b[c("floor", "ceiling", "supported", "n")])
  # Reliability is for cells holding the layer's points, not the subsample's.
  vm  <- attr(sac, "variogram_model")
  cf  <- spatialkit:::.vgm_correlation_fn(vm)
  bb  <- sf::st_bbox(pts)
  rbV <- spatialkit:::.rbar_rect(cf, bb[["xmax"]] - bb[["xmin"]], bb[["ymax"]] - bb[["ymin"]])
  expect_equal(prof$reliability,
               vapply(prof$levels, spatialkit:::.reliability_at, numeric(1),
                      area = b$area, n_total = 1200, nugget = 0.5, psill = 1,
                      cor_fn = cf, rbar_V = rbV))
  # Cp: the subsample's RSS with its own optimism added back, plus the
  # variance of cell means built from all 1200 points.  Every subsample cell
  # holds a scored row here, so L_m = L.
  expect_equal(prof$cp, prof$rss / 300 + 0.5 * (prof$levels / 300 + prof$levels / 1200))
})


test_that("select_on = 'split' chooses a count for the layer it is applied to", {
  # Four clusters at the corners of a 10 km square.  Every spatial half holds
  # two of them, and determine_optimal_levels(select_on = "split") used to
  # run on that half and answer 2, which the documented workflow then
  # applied to a tessellation of all four clusters.
  set.seed(1)
  cen <- expand.grid(cx = c(1000, 9000), cy = c(1000, 9000))
  g <- rep(1:4, each = 60)
  d <- data.frame(x = rnorm(240, cen$cx[g], 300), y = rnorm(240, cen$cy[g], 300),
                  a = rnorm(240), b = rnorm(240))
  d$z <- d$a + rnorm(240)
  pts <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
  all <- suppressWarnings(determine_optimal_levels(
    pts, max_levels = 12, response_var = "z", predictor_vars = c("a", "b"), set_seed = 1))
  spl <- suppressWarnings(determine_optimal_levels(
    pts, max_levels = 12, response_var = "z", predictor_vars = c("a", "b"),
    select_on = "split", set_seed = 1))
  expect_identical(as.integer(all[1]), 4L)
  expect_identical(as.integer(spl[1]), 4L)
})


test_that("select_on = 'split' reads the selection half's response on the whole layer's cells", {
  # The Moran pass is the only step that reads the response.  Record what it
  # is handed: the selection half's rows, labelled by k-means cells fitted to
  # every point.  On the half's own partition every level would label the
  # half's rows with all k cells; on the whole layer's, a spatial half falls
  # in only some of them.
  set.seed(8)
  n <- 200
  d <- data.frame(x = runif(n, 0, 5000), y = runif(n, 0, 5000), w = rnorm(n))
  d$z <- 2 * d$w + rnorm(n)
  pts <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
  seen <- list()
  local_mocked_bindings(
    .morans_i_for_k = function(xy, response, predictors, cluster_ids) {
      seen[[length(seen) + 1L]] <<- list(resp = response, cl = cluster_ids)
      c(I = NA_real_, z = NA_real_)
    },
    .package = "spatialkit")
  out <- suppressWarnings(determine_optimal_levels(
    pts, max_levels = 12, response_var = "z", predictor_vars = "w",
    criterion = "morans_i", select_on = "split", set_seed = 2))
  sel <- attr(out, "split")$selection
  expect_true(length(seen) > 0L)
  for (s in seen) expect_identical(s$resp, pts$z[sel])
  ks <- vapply(seen, function(s) max(s$cl), integer(1))
  used <- vapply(seen, function(s) length(unique(s$cl)), integer(1))
  expect_true(any(used < ks))
})


test_that("a split profile is bounded on the whole layer and reads one half's response", {
  # Three clusters of 210, 60 and 30 points.  Profiled on its selection half,
  # the split used to take the half's hull (1.6e6 against 4.6e7 m^2), the
  # half's point count and the half's clusters, and so bounded and judged a
  # different layer from the one the count was applied to.
  set.seed(10)
  mk <- function(n, cx, cy) data.frame(x = rnorm(n, cx, 150), y = rnorm(n, cy, 150))
  d <- rbind(mk(210, 2000, 2000), mk(60, 8000, 3000), mk(30, 5000, 8000))
  d$z <- 0.0002 * d$x + rnorm(300)
  pts <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
  sac <- rr_sac(pts, range = 1500, nugget = 1, psill = 0.5)
  prof <- function(p, ...) resolution_profile(p, "z", sac = sac, n_levels = 5, nstart = 5, ...)
  all <- prof(pts)
  spl <- suppressWarnings(prof(pts, select_on = "split"))
  keys <- c("floor", "ceiling", "ceiling_from", "supported", "area", "n", "n_distinct")
  expect_identical(attr(spl, "bounds")[keys], attr(all, "bounds")[keys])
  expect_identical(spl$levels, all$levels)
  expect_identical(spl$wss, all$wss)
  # The estimation half's response is never read; the selection half's is.
  s <- attr(spl, "split")
  pts_e <- pts; pts_e$z[s$estimation] <- rnorm(length(s$estimation), 50, 20)
  spl_e <- suppressWarnings(prof(pts_e, select_on = "split"))
  expect_identical(spl_e[c("rss", "cp", "moran_z")], spl[c("rss", "cp", "moran_z")])
  pts_s <- pts; pts_s$z[s$selection] <- rnorm(length(s$selection), 50, 20)
  spl_s <- suppressWarnings(prof(pts_s, select_on = "split"))
  expect_false(isTRUE(all.equal(spl_s$rss, spl$rss)))
})


test_that("a WSS curve that falls like c / k has no elbow", {
  # The linear-axis chord rule answers sqrt(k_min * k_max) on c / k, here 4,
  # and that used to be returned as the knee.
  eb <- spatialkit:::.elbow_from_wss(1000 / (1:12))
  expect_false(eb$structured)
  expect_true(all(abs(eb$diagnostics$sag) < 1e-12))
  # A bend at the cluster count on the same scale is one.
  bent <- spatialkit:::.elbow_from_wss(c(1000, 400, 150, 60 * 4 / (4:12)))
  expect_true(bent$structured)
  expect_identical(bent$knee_k, 4L)
})


test_that("determine_optimal_levels() warns when uniform points have no elbow", {
  set.seed(3)
  n <- 400
  pts <- sf::st_as_sf(data.frame(x = runif(n, 0, 10000), y = runif(n, 0, 10000)),
                      coords = c("x", "y"), crs = 32632)
  for (ml in c(12, 40)) {
    expect_warning(k <- determine_optimal_levels(pts, max_levels = ml),
                   "no elbow.*set by the ladder")
    expect_type(k, "integer")
  }
})


test_that("determine_optimal_levels() finds 2 to 8 separated clusters without a warning", {
  # At the default max_levels = 12 the linear chord rule answered 3 or 4 for
  # eight clusters: the fall of the between-cluster WSS outweighed the knee.
  for (K in c(2L, 3L, 4L, 8L)) {
    set.seed(K)
    ctr <- expand.grid(x = seq(0, by = 2500, length.out = 4), y = c(0, 2500))[seq_len(K), ]
    g <- rep(seq_len(K), each = 40)
    pts <- sf::st_as_sf(data.frame(x = ctr$x[g] + rnorm(40 * K, sd = 100),
                                   y = ctr$y[g] + rnorm(40 * K, sd = 100)),
                        coords = c("x", "y"), crs = 32632)
    expect_no_warning(k <- determine_optimal_levels(pts, max_levels = 12))
    expect_identical(k[1], K, info = paste("K =", K))
  }
})


test_that("a geometry-only profile of uniform points names no cell count", {
  set.seed(5)
  n <- 600
  pts <- sf::st_as_sf(data.frame(x = runif(n, 0, 10000), y = runif(n, 0, 10000)),
                      coords = c("x", "y"), crs = 32632)
  prof <- resolution_profile(pts, n_levels = 8, nstart = 5)
  expect_true(all(is.na(prof$elbow)))
  expect_output(print(prof), "elbow       : none")
  expect_error(select_resolution(prof, "elbow"), "no elbow")
  # build_tessellation() used to read the chord rule's sqrt(first x last
  # level) off this profile, 16 to 18 cells on any large uniform layer.
  bnd <- sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(pts)))
  expect_error(build_tessellation(pts, boundary = bnd, method = "hex",
                                  approx_n_cells = prof, quiet = TRUE),
               "no cluster structure")
  expect_error(get_voronoi_seeds(bnd, method = "kmeans", n = prof,
                                 sample_points = pts, set_seed = 1),
               "no cluster structure")
})
