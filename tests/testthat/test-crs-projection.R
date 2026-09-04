# ===========================================================================
# CRS detection, projection, and harmonisation.
# Choosing a projection is upstream of everything the package measures:
# variogram ranges, block sizes, GWR bandwidths and GP length-scales are
# all in CRS units. See also test-crs-selection.R for the wide-extent rules.
# ===========================================================================

test_that(".is_longlat detects geographic CRS", {
  pts <- sf::st_sfc(sf::st_point(c(-122.4, 37.8)), crs = 4326)
  expect_true(spatialkit:::.is_longlat(pts))

  pts_proj <- sf::st_transform(pts, 32610)
  expect_false(spatialkit:::.is_longlat(pts_proj))
})


test_that("ensure_projected transforms geographic data", {
  pts <- sf::st_sf(
    id = 1:3,
    geometry = sf::st_sfc(
      sf::st_point(c(-122.4, 37.8)),
      sf::st_point(c(-122.3, 37.7)),
      sf::st_point(c(-122.5, 37.9)),
      crs = 4326
    )
  )
  result <- ensure_projected(pts)
  expect_false(sf::st_is_longlat(result))
})


test_that("ensure_projected does NOT assume 4326 for small integer-like coords", {
  # Simulates a local site survey in metres with coords outside the

  # geographic envelope (values > 180 / > 90) so the heuristic cannot
  # mistakenly interpret them as lon/lat.
  pts <- sf::st_sf(
    id = 1:3,
    geometry = sf::st_sfc(
      sf::st_point(c(500, 600)),
      sf::st_point(c(510, 620)),
      sf::st_point(c(505, 610))
    )
  )
  result <- suppressWarnings(ensure_projected(pts))
  # CRS should still be NA — the function must NOT guess 4326

  expect_true(is.na(sf::st_crs(result)))
})


test_that("ensure_projected assumes 4326 for CRS-less data with geographic precision", {
  # Realistic lon/lat coordinates without a CRS set
  pts <- sf::st_sf(
    id = 1:2,
    geometry = sf::st_sfc(
      sf::st_point(c(-0.1278, 51.5074)),
      sf::st_point(c(-0.1385, 51.5013))
    )
  )
  result <- suppressWarnings(ensure_projected(pts))
  # Should have been assigned a projected CRS
  expect_false(is.na(sf::st_crs(result)))
  expect_false(sf::st_is_longlat(result))
})


test_that("ensure_projected assumes 4326 for CRS-less data with large extent", {
  pts <- sf::st_sf(
    id = 1:2,
    geometry = sf::st_sfc(
      sf::st_point(c(10, 20)),
      sf::st_point(c(15, 25))
    )
  )
  result <- suppressWarnings(ensure_projected(pts))
  # Extent > 1 degree in both axes → should assume 4326
  expect_false(is.na(sf::st_crs(result)))
  expect_false(sf::st_is_longlat(result))
})


# ---------------------------------------------------------------------------
# harmonize_crs()
#
# The test point is (10 E, 50 N), NOT the origin.  (0, 0) is the origin of both
# EPSG:4326 and EPSG:3857, so its coordinates do not move under the transform
# -- an implementation that called st_set_crs() instead of st_transform() (the
# exact hazard R/crs-geometry.R warns about) would relabel the CRS and pass a
# test written on that point identically.
#
# The expected projected coordinates are computed from the Web Mercator
# forward equations rather than from sf, so the assertion does not check
# st_transform() against itself:
#     x = R * lambda,  y = R * ln(tan(pi/4 + phi/2)),  R = 6378137
# ---------------------------------------------------------------------------

.wm_R    <- 6378137
.wm_x    <- .wm_R * (10 * pi / 180)
.wm_y    <- .wm_R * log(tan(pi / 4 + (50 * pi / 180) / 2))
.ll_pt   <- function() sf::st_sf(id = 1L,
                                 geometry = sf::st_sfc(sf::st_point(c(10, 50)),
                                                       crs = 4326))
.wm_pt   <- function() sf::st_sf(id = 2L,
                                 geometry = sf::st_sfc(sf::st_point(c(.wm_x, .wm_y)),
                                                       crs = 3857))

test_that("harmonize_crs(prefer = 'b') reprojects a onto b's CRS", {
  a <- .ll_pt(); b <- .wm_pt()
  result <- harmonize_crs(a, b, prefer = "b")

  expect_equal(sf::st_crs(result$a), sf::st_crs(3857))
  expect_equal(sf::st_crs(result$b), sf::st_crs(3857))

  # `a` MOVED: metres, not degrees, matching the hand-computed projection.
  expect_equal(unname(sf::st_coordinates(result$a)[1, 1:2]),
               c(.wm_x, .wm_y), tolerance = 1e-6)
  # `b` was already in the target CRS and must be untouched.
  expect_equal(unname(sf::st_coordinates(result$b)[1, 1:2]),
               c(.wm_x, .wm_y), tolerance = 1e-9)
})

test_that("harmonize_crs(prefer = 'a') reprojects b onto a's CRS", {
  a <- .ll_pt(); b <- .wm_pt()
  result <- harmonize_crs(a, b, prefer = "a")

  expect_equal(sf::st_crs(result$a), sf::st_crs(4326))
  expect_equal(sf::st_crs(result$b), sf::st_crs(4326))
  expect_equal(unname(sf::st_coordinates(result$b)[1, 1:2]), c(10, 50),
               tolerance = 1e-7)
  expect_equal(unname(sf::st_coordinates(result$a)[1, 1:2]), c(10, 50),
               tolerance = 1e-9)
})

test_that("harmonize_crs(target_crs) moves BOTH inputs to a third CRS", {
  a <- .ll_pt(); b <- .wm_pt()
  result <- harmonize_crs(a, b, target_crs = 32632)   # UTM 32N

  expect_equal(sf::st_crs(result$a), sf::st_crs(32632))
  expect_equal(sf::st_crs(result$b), sf::st_crs(32632))
  # Both describe the same place, so they must land on the same coordinates.
  expect_equal(unname(sf::st_coordinates(result$a)[1, 1:2]),
               unname(sf::st_coordinates(result$b)[1, 1:2]), tolerance = 1e-3)
  # ... and that place is not either input's original coordinates.
  expect_gt(abs(sf::st_coordinates(result$a)[1, 1] - 10), 1000)
})

test_that("harmonize_crs stamps a CRS-less input and says so", {
  # st_set_crs() relabels without moving anything, so the assumption has to be
  # announced.  .log_warn() raises no R condition -- see helper-logging.R.
  b <- .wm_pt()
  a_na <- sf::st_sf(id = 1L, geometry = sf::st_sfc(sf::st_point(c(10, 50))))

  expect_warning(lines <- capture_spatialkit_log(res <- harmonize_crs(a_na, b)),
                 "WITHOUT reprojection")
  expect_true(log_has(lines, "`a` has no CRS"))
  expect_true(log_has(lines, "WITHOUT reprojection"))
  expect_equal(sf::st_crs(res$a), sf::st_crs(3857))
  # Stamped, not transformed: the coordinates are exactly as supplied.
  expect_equal(unname(sf::st_coordinates(res$a)[1, 1:2]), c(10, 50))

  # The mirror branch, with the CRS-less object second.
  expect_warning(lines_b <- capture_spatialkit_log(res_b <- harmonize_crs(b, a_na)),
                 "WITHOUT reprojection")
  expect_true(log_has(lines_b, "`b` has no CRS"))
  expect_equal(sf::st_crs(res_b$b), sf::st_crs(3857))
  expect_equal(unname(sf::st_coordinates(res_b$b)[1, 1:2]), c(10, 50))

  # Neither input has a CRS: nothing to align, nothing to warn about.
  b_na <- sf::st_sf(id = 2L, geometry = sf::st_sfc(sf::st_point(c(1, 2))))
  lines_none <- capture_spatialkit_log(res_none <- harmonize_crs(a_na, b_na))
  expect_false(log_has(lines_none, "has no CRS"))
  expect_true(is.na(sf::st_crs(res_none$a)))
  expect_true(is.na(sf::st_crs(res_none$b)))
})

test_that("harmonize_crs on_transform_error controls a failed transform", {
  a <- .ll_pt(); b <- .wm_pt()
  # A syntactically parseable but geometrically impossible projection: PROJ
  # rejects |lat_0| > 90, so st_transform() fails rather than returning
  # something wrong.
  bad_crs <- paste("+proj=omerc +lat_0=91 +lonc=0 +alpha=0 +k=1",
                   "+x_0=0 +y_0=0 +datum=WGS84")

  # Default: surface the failure.
  expect_error(
    suppressWarnings(harmonize_crs(a, b, target_crs = bad_crs)),
    "st_transform\\(\\) failed"
  )

  # Opt-in override: no error, but the unsafe fallback is announced.
  lines <- capture_spatialkit_log(
    suppressWarnings(harmonize_crs(a, b, target_crs = bad_crs,
                                   on_transform_error = "set_crs"))
  )
  expect_true(log_has(lines, "coordinates are NOT reprojected"))
})

test_that("harmonize_crs rejects non-spatial input", {
  a <- .ll_pt()
  expect_error(harmonize_crs(data.frame(x = 1), a), "`a` must be sf or sfc")
  expect_error(harmonize_crs(a, data.frame(x = 1)), "`b` must be sf or sfc")
})


test_that("ensure_projected refuses to assume lon/lat for a small integer-grid survey", {
  # With no CRS at all, the only evidence available is the numbers themselves.
  # Falling inside the lon/lat envelope is necessary but nowhere near
  # sufficient: a local site survey in metres sits inside it too, and stamping
  # EPSG:4326 on one silently teleports it to the Gulf of Guinea and then
  # "projects" it.  So a second test is required -- the data must look
  # POSITIVELY geographic, either by spanning more than a degree or by
  # carrying fractional-degree precision.
  nocrs <- function(x, y)
    sf::st_as_sf(data.frame(x = x, y = y), coords = c("x", "y"))

  # A 1 m quadrat frame on a whole-metre grid: extent 1 in both axes, no
  # fractional part anywhere.  Neither test is met, so nothing is assumed.
  quad <- nocrs(c(0, 0, 1, 1), c(0, 1, 0, 1))
  lines <- capture_spatialkit_log(out <- ensure_projected(quad))
  expect_true(log_has(lines, "Not assuming EPSG:4326"))
  expect_true(log_has(lines, "lack the decimal precision or extent"))
  expect_true(is.na(sf::st_crs(out)))
  expect_equal(sf::st_coordinates(out), sf::st_coordinates(quad))

  # (a) More than a degree in ONE axis is enough to tip it the other way, so
  #     the extent test is doing real work and is not simply always false.
  wide_x <- nocrs(c(0, 0, 2, 2), c(0, 1, 0, 1))
  expect_warning(lx <- capture_spatialkit_log(ox <- ensure_projected(wide_x)),
                 "assuming EPSG:4326")
  expect_true(log_has(lx, "assuming EPSG:4326"))
  expect_false(is.na(sf::st_crs(ox)))
  expect_false(sf::st_is_longlat(ox))

  wide_y <- nocrs(c(0, 0, 1, 1), c(0, 2, 0, 2))
  expect_false(is.na(sf::st_crs(
    suppressWarnings(ensure_projected(wide_y)))))

  # (b) So is fractional-degree precision within the same 1-degree box.
  fine <- nocrs(c(9.10, 9.20, 9.30, 9.15), c(48.70, 48.75, 48.80, 48.72))
  lf <- capture_spatialkit_log(of <- suppressWarnings(ensure_projected(fine)))
  expect_true(log_has(lf, "assuming EPSG:4326"))
  expect_equal(sf::st_crs(of)$epsg, 32632L)     # the containing UTM zone

  # An explicit CRS short-circuits the whole heuristic: the quadrat is left
  # exactly where it is, with no guessing and no warning about it.
  stamped <- quad
  sf::st_crs(stamped) <- sf::st_crs(32632)
  quiet <- capture_spatialkit_log(os <- ensure_projected(stamped))
  expect_false(log_has(quiet, "assuming EPSG:4326"))
  expect_equal(sf::st_crs(os), sf::st_crs(32632))
  expect_equal(sf::st_coordinates(os), sf::st_coordinates(quad))
})


# ---------------------------------------------------------------------------
# One interpretation of CRS-less coordinates, applied everywhere
# ---------------------------------------------------------------------------

.crsless_lonlat_pts <- function(n = 100, seed = 3) {
  # Decimal-degree-looking coordinates with NO CRS: a layer that the lon/lat
  # heuristic takes as EPSG:4326.
  set.seed(seed)
  d <- data.frame(x = runif(n, 24.1, 25.9), y = runif(n, 40.1, 41.4),
                  elev = rnorm(n))
  d$price <- 10 + 3 * d$elev + 2 * sin(d$x * 4) + rnorm(n, 0, 0.3)
  sf::st_as_sf(d, coords = c("x", "y"), crs = NA)
}

test_that("predict() places CRS-less newdata where the fit placed the training rows", {
  # prep_model_data() assumed lon/lat and PROJECTED the training data; every
  # predict() method then passed the fit's CRS as a target, and the target
  # branch used to STAMP it onto the raw numbers.  The same rows sat in two
  # different places, and predict(fit, newdata = training rows) disagreed with
  # fitted(fit) by up to one response SD (R2 0.98 in-sample, 0.64 via newdata).
  skip_if_not_installed("ranger")
  pts <- .crsless_lonlat_pts()
  expect_warning(
    fit <- fit_rf_model(pts, "price", "elev", include_coords = TRUE,
                        num_trees = 200, seed = 1),
    "assuming EPSG:4326")
  # data_sf carries the assumption the fit was built under ...
  expect_identical(attr(fit$data_sf, "crs_assumed"), "EPSG:4326")
  # ... and newdata drawn from the raw, CRS-less rows lands where they did.
  # The replayed assumption is announced too.
  expect_warning(p_raw <- predict(fit, newdata = pts), "interpreting it as EPSG:4326")
  p_fit <- predict(fit, newdata = fit$data_sf)
  expect_equal(p_raw, p_fit, tolerance = 1e-10)
  # A SUBSET whose own bounding box would not trigger the heuristic still
  # gets the training-time interpretation, not a stamp.
  sub <- pts[1:5, ]
  expect_equal(suppressWarnings(predict(fit, newdata = sub)), p_fit[1:5], tolerance = 1e-10)
})

test_that("ensure_projected() interprets CRS-less lon/lat the same with and without a target", {
  pts  <- .crsless_lonlat_pts()
  expect_warning(auto <- ensure_projected(pts), "assuming EPSG:4326")   # no target
  expect_warning(tgt  <- ensure_projected(pts, target_crs = sf::st_crs(auto)),
                 "taken as EPSG:4326")                                   # target given
  expect_equal(sf::st_coordinates(tgt), sf::st_coordinates(auto), tolerance = 1e-9)
  expect_identical(attr(tgt, "crs_assumed"), "EPSG:4326")
  # Coordinates that do NOT look like lon/lat -- inside the envelope but with
  # neither decimal precision nor an extent above one degree, the site-survey
  # case the heuristic exists to leave alone -- are still stamped, as before.
  loc <- sf::st_as_sf(data.frame(x = c(0, 1, 0, 1), y = c(0, 0, 1, 1), v = 1),
                      coords = c("x", "y"), crs = NA)
  st  <- ensure_projected(loc, target_crs = sf::st_crs(32632))
  expect_equal(sf::st_coordinates(st),
               cbind(X = c(0, 1, 0, 1), Y = c(0, 0, 1, 1)), ignore_attr = TRUE)
  expect_null(attr(st, "crs_assumed"))
})

test_that("CRS-less LINESTRINGs inside the lon/lat envelope get midpoints", {
  # The temporary projection used for 'halfway along the line' was reversed
  # with st_transform(<midpoints>, NA_crs_), which errors: 'crs not found'.
  set.seed(4)
  n <- 12
  ls <- lapply(seq_len(n), function(i) {
    x0 <- runif(1, 24, 26); y0 <- runif(1, 40, 41)
    sf::st_linestring(rbind(c(x0, y0), c(x0 + 0.01, y0 + 0.02)))
  })
  lsf <- sf::st_sf(v = rnorm(n), geometry = sf::st_sfc(ls, crs = sf::NA_crs_))
  expect_no_error(mp <- suppressWarnings(coerce_to_points(lsf, "line_midpoint")))
  expect_true(all(sf::st_geometry_type(mp) == "POINT"))
  expect_true(is.na(sf::st_crs(mp)))                 # output space == input space
  xy <- sf::st_coordinates(mp)
  expect_true(all(xy[, 1] > 24 & xy[, 1] < 26.1))   # back in the input's numbers
  expect_no_error(suppressWarnings(prep_model_data(lsf, "v", character(0))))
})

test_that("hex and square grids accept the CRS-less input voronoi accepts", {
  # The points were left CRS-less while create_grid_polygons() projected the
  # boundary on its own; st_intersects() then refused the pair.
  pts <- .crsless_lonlat_pts(60)
  bnd <- sf::st_as_sfc(sf::st_bbox(pts))
  bnd <- sf::st_set_crs(bnd, NA)
  for (m in c("voronoi", "hex", "square")) {
    expect_no_error(t <- suppressWarnings(
      build_tessellation(pts, boundary = bnd, method = m,
                         approx_n_cells = 12, quiet = TRUE)), message = m)
    expect_gt(nrow(t$cells), 0)
  }
})
