# tests/testthat/test-fix-crs-unlocated-readers.R
# ---------------------------------------------------------------------------
# The CRS readers read the rows that have a location.  A coordinate that is
# not finite or lies beyond 1e150 is not one.
#
# BUG-31: one infinite coordinate made a layer without a CRS read as planar
# (the lon/lat test returned FALSE and no box when the layer's box was not
# finite), and a boundary without a CRS with such a vertex was given the
# points' geographic CRS in silence.
# BUG-136: a finite coordinate beyond 1e150 put the box off the lon/lat
# envelope, with the same result.
# BUG-137: the local CRS chosen for a lon/lat layer read such rows: a
# longitude of 1e200 stopped ensure_projected() with sf's bare "invalid crs",
# a latitude of 1e200 moved the projection of every other row.
# ---------------------------------------------------------------------------

# The value, and the messages of the R warnings, of an expression; the log
# lines are not R warnings.
with_warnings <- function(expr) {
  w <- character(0)
  v <- withCallingHandlers(expr, warning = function(c) {
    w <<- c(w, conditionMessage(c))
    invokeRestart("muffleWarning")
  })
  list(value = v, warnings = w)
}
pts_of <- function(d, crs = NA) sf::st_as_sf(d, coords = c("x", "y"), crs = crs)
xy_of  <- function(x) unname(sf::st_coordinates(x)[, 1:2, drop = FALSE])

deg  <- data.frame(x = c(-93.5, -93.1, -92.8, -93.3, -92.9), y = c(44.9, 45.1, 44.7, 45.3, 44.8))
five <- pts_of(deg)
# Built with st_as_sf(), so the layer's box is the one of its coordinates.
with_row <- function(d, at) pts_of(rbind(d, data.frame(x = at[1L], y = at[2L])))
# The 300 points of the scan, near (-80.4, 45.6), with row 1 replaced.
near_ontario <- function() {
  set.seed(4402)
  data.frame(x = stats::runif(300, -81, -79.7), y = stats::runif(300, 45.15, 46.05),
             z = stats::rnorm(300))
}
row1_at <- function(at, crs = NA) {
  d <- near_ontario()
  p <- pts_of(d, crs)
  g <- sf::st_geometry(p)
  g[1L] <- sf::st_sfc(if (is.null(at)) sf::st_point() else sf::st_point(at), crs = sf::st_crs(g))
  sf::st_set_geometry(p, g)
}
unlocated <- list(inf = c(Inf, 45), ninf = c(-Inf, 45), far_lon = c(1e200, 45),
                  nfar_lon = c(-1e200, 45), far_lat = c(-80, 1e200))
ring <- function(m) sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(m))))

# ---- the shared helper -----------------------------------------------------

test_that(".unlocated_features() reads x and y of every vertex against the bound", {
  g <- sf::st_sfc(sf::st_point(c(1, 2)), sf::st_point(), sf::st_point(c(Inf, 3)),
                  sf::st_point(c(1, 1e200)), sf::st_point(c(1e150, -1e150)))
  expect_identical(.unlocated_features(g), c(FALSE, FALSE, TRUE, TRUE, FALSE))
  # Z is not read.
  gz <- sf::st_sfc(sf::st_point(c(1, 2, Inf)), sf::st_point(c(Inf, 2, 0)))
  expect_identical(.unlocated_features(gz), c(FALSE, TRUE))
  sq <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0))
  far <- sq; far[3L, 1L] <- 1e200
  gp <- sf::st_sfc(sf::st_polygon(list(sq)), sf::st_polygon(list(far)),
                   sf::st_multipolygon(list(list(sq), list(far + 2))),
                   sf::st_linestring(rbind(c(0, 0), c(Inf, 1))),
                   sf::st_multipoint(rbind(c(0, 0), c(1, -Inf))),
                   sf::st_geometrycollection(list(sf::st_point(c(0, 0)),
                                                  sf::st_linestring(rbind(c(0, 0), c(1, 1e151))))),
                   sf::st_polygon())
  expect_identical(.unlocated_features(gp), c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE))
  expect_identical(.unlocated_features(sf::st_sfc()), logical(0))
})

test_that(".empty_unlocated() empties those rows in their own type and dimension", {
  sq <- rbind(c(0, 0, 5), c(1, 0, 5), c(1, 1, 5), c(0, 0, 5))
  far <- sq; far[2L, 2L] <- Inf
  g <- sf::st_sfc(sf::st_polygon(list(sq)), sf::st_polygon(list(far)), crs = 4326)
  e <- .empty_unlocated(g)
  expect_s3_class(e, "sfc_POLYGON")
  expect_identical(sf::st_is_empty(e), c(FALSE, TRUE))
  expect_identical(class(e[[2L]]), class(g[[2L]]))
  expect_identical(e[[1L]], g[[1L]])
  expect_identical(sf::st_crs(e), sf::st_crs(g))
  p <- sf::st_sf(v = 1:3, geometry = sf::st_sfc(sf::st_point(c(1, 2, 3)), sf::st_point(c(1e200, 2, 3)),
                                                sf::st_point(c(4, 5, 6))))
  ep <- .empty_unlocated(p)
  expect_identical(ep$v, 1:3)
  expect_identical(sf::st_is_empty(ep), c(FALSE, TRUE, FALSE))
  expect_identical(class(sf::st_geometry(ep)[[2L]]), c("XYZ", "POINT", "sfg"))
  # A layer whose rows all have a location comes back as it is.
  expect_identical(.empty_unlocated(five), five)
  # Rows can be named.
  expect_identical(sf::st_is_empty(.empty_unlocated(five, 2L)), c(FALSE, TRUE, FALSE, FALSE, FALSE))
  # Mixed types keep each row's type.
  gm <- sf::st_sfc(sf::st_point(c(Inf, 0)), sf::st_linestring(rbind(c(0, 0), c(1, 1e200))),
                   sf::st_point(c(1, 1)))
  em <- .empty_unlocated(gm)
  expect_identical(as.character(sf::st_geometry_type(em)), c("POINT", "LINESTRING", "POINT"))
  expect_identical(sf::st_is_empty(em), c(TRUE, TRUE, FALSE))
})

# ---- BUG-31 and BUG-136: the lon/lat test of a layer without a CRS ----------

test_that("the lon/lat test reads the box of the coordinates that have a location", {
  ref <- .looks_like_lonlat(five)
  expect_true(ref$lonlat)
  for (at in unlocated) {
    x <- with_row(deg, at)
    expect_false(all(is.finite(sf::st_bbox(x))) && all(abs(sf::st_bbox(x)) <= 1e150))
    expect_identical(.looks_like_lonlat(x), ref, label = paste(at, collapse = ", "))
  }
  # Under a degree across, where the precision test decides: no "probable
  # complete loss of accuracy in modulus" from 1e200 %% 1.
  small <- data.frame(x = c(-93.25, -93.21, -93.13), y = c(44.95, 44.98, 44.91))
  for (at in unlocated) {
    expect_no_warning(r <- .looks_like_lonlat(with_row(small, at)))
    expect_identical(r, .looks_like_lonlat(pts_of(small)))
  }
  # Planar metres keep their own located box, not NULL.
  m <- data.frame(x = c(4e5, 4.6e5, 4.3e5), y = c(4.9e6, 5e6, 4.95e6))
  for (at in list(c(Inf, 4.9e6), c(4e5, 1e200))) {
    r <- .looks_like_lonlat(with_row(m, at))
    expect_false(r$lonlat)
    expect_identical(r$bb, sf::st_bbox(pts_of(m)))
  }
  # A coordinate of exactly 1e150 has a location (as .values_at_xy() reads
  # it): the box reaches it.
  r <- .looks_like_lonlat(pts_of(rbind(m, data.frame(x = c(1e150, 4e5, Inf), y = c(4.9e6, -1e150, 5e6)))))
  expect_false(r$lonlat)
  expect_identical(unname(as.numeric(r$bb)), c(4e5, -1e150, 1e150, 5e6))
  r <- .looks_like_lonlat(pts_of(rbind(m, data.frame(x = c(-1e150, Inf), y = c(1e150, 5e6)))))
  expect_identical(unname(as.numeric(r$bb)), c(-1e150, 4.9e6, 4.6e5, 1e150))
  # No coordinate with a location: as before.
  expect_identical(.looks_like_lonlat(pts_of(data.frame(x = c(Inf, -Inf), y = c(45, 45)))),
                   list(lonlat = FALSE, bb = NULL))
  expect_identical(.looks_like_lonlat(pts_of(data.frame(x = c(1e200, 3), y = c(45, -1e300)))),
                   list(lonlat = FALSE, bb = NULL))
  expect_identical(.looks_like_lonlat(five[0, ]), list(lonlat = FALSE, bb = NULL))
})

test_that("ensure_projected() reads a layer without a CRS as the same layer with the row empty", {
  for (at in unlocated) {
    lab <- paste(at, collapse = ", ")
    got  <- with_warnings(ensure_projected(row1_at(at)))
    want <- with_warnings(ensure_projected(row1_at(NULL)))
    expect_identical(sf::st_crs(got$value), sf::st_crs(32617), label = lab)
    expect_identical(attr(got$value, "crs_assumed"), "EPSG:4326", label = lab)
    expect_identical(attr(got$value, "crs_choice"), attr(want$value, "crs_choice"), label = lab)
    expect_identical(xy_of(got$value)[-1L, ], xy_of(want$value)[-1L, ], label = lab)
    expect_identical(got$warnings, want$warnings, label = lab)
    # The row stays without a location.
    expect_false(any(is.finite(xy_of(got$value)[1L, ])), label = lab)
  }
  # The five points and a sixth at Inf: the five's CRS, coordinates and warning.
  six  <- with_warnings(ensure_projected(with_row(deg, c(Inf, 45))))
  ref  <- with_warnings(ensure_projected(five))
  expect_identical(sf::st_crs(six$value), sf::st_crs(32615))
  expect_identical(xy_of(six$value)[1:5, ], xy_of(ref$value))
  expect_identical(six$warnings, ref$warnings)
  expect_length(six$warnings, 1L)
})

test_that("with a target CRS, the readers reproject such a layer, and do not stamp it", {
  b32615 <- sf::st_transform(sf::st_set_crs(five, 4326), 32615)
  for (at in unlocated) {
    lab <- paste(at, collapse = ", ")
    x <- with_row(deg, at)
    calls <- list(
      function(l) ensure_projected(l, target_crs = 32615),
      function(l) .transform_or_stamp(l, sf::st_crs(32615), "points", "probe"),
      function(l) harmonize_crs(l, b32615)$a,
      function(l) harmonize_crs(l, b32615, prefer = "b")$a)
    for (f in calls) {
      got <- with_warnings(f(x))
      ref <- with_warnings(f(five))
      expect_identical(sf::st_crs(got$value), sf::st_crs(32615), label = lab)
      expect_identical(xy_of(got$value)[1:5, ], xy_of(ref$value), label = lab)
      expect_identical(got$warnings, ref$warnings, label = lab)
      expect_false(any(grepl("do not look like lon/lat", got$warnings)), label = lab)
    }
  }
  # No coordinate with a location: stamped, as before.
  allinf <- pts_of(data.frame(x = c(Inf, -Inf), y = c(45, 45)))
  r <- with_warnings(.transform_or_stamp(allinf, sf::st_crs(32615), "points", "probe"))
  expect_identical(sf::st_crs(r$value), sf::st_crs(32615))
  expect_true(any(grepl("do not look like lon/lat", r$warnings)))
  e <- with_warnings(ensure_projected(allinf))
  expect_true(is.na(sf::st_crs(e$value)))
  expect_identical(attr(e$value, "crs_assumed"), "none")
  expect_length(e$warnings, 0L)
})

test_that("a boundary without a CRS is read on the vertices that have a location", {
  bm     <- ring(rbind(c(4e5, 4.9e6), c(4.6e5, 4.9e6), c(4.6e5, 5e6), c(4e5, 5e6), c(4e5, 4.9e6)))
  bm_inf <- ring(rbind(c(4e5, 4.9e6), c(4.6e5, 4.9e6), c(Inf, 5e6), c(4e5, 5e6), c(4e5, 4.9e6)))
  bm_far <- ring(rbind(c(4e5, 4.9e6), c(4.6e5, 4.9e6), c(4.6e5, 1e200), c(4e5, 5e6), c(4e5, 4.9e6)))
  bd     <- ring(rbind(c(-94, 44), c(-92, 44), c(-92, 46), c(-94, 46), c(-94, 44)))
  bd_inf <- ring(rbind(c(-94, 44), c(-92, 44), c(Inf, 46), c(-92, 46), c(-94, 46), c(-94, 44)))
  ll <- sf::st_crs(4326)
  msg_of <- function(b) tryCatch(.crsless_boundary_as_lonlat(b, ll, "probe"),
                                 error = function(e) conditionMessage(e))
  expect_match(msg_of(bm), "are not lon/lat", fixed = TRUE)
  expect_identical(msg_of(bm_inf), msg_of(bm))
  expect_identical(msg_of(bm_far), msg_of(bm))
  ref <- with_warnings(.crsless_boundary_as_lonlat(bd, ll, "probe"))
  got <- with_warnings(.crsless_boundary_as_lonlat(bd_inf, ll, "probe"))
  expect_identical(sf::st_crs(got$value), ll)
  expect_identical(got$warnings, ref$warnings)
  expect_length(got$warnings, 1L)
  # cross-validation's boundary reader
  ref <- with_warnings(.cv_boundary_crs(bd, "probe"))
  got <- with_warnings(.cv_boundary_crs(bd_inf, "probe"))
  expect_identical(sf::st_crs(got$value), ll)
  expect_identical(got$warnings, ref$warnings)
  # Through a builder (it stopped with "`boundary` must be polygonal"): a
  # row with a vertex at Inf has no location as a whole, so the builder
  # leaves it out as given; beside the boundary without that vertex it is
  # refused with that boundary's message, and alone, as a boundary with no
  # polygon that has a location.
  set.seed(3)
  pts <- sf::st_as_sf(data.frame(x = stats::runif(20, -93.9, -92.1), y = stats::runif(20, 44.1, 45.9)),
                      coords = c("x", "y"), crs = 4326)
  err_of <- function(b) tryCatch(build_tessellation(pts, method = "voronoi", boundary = b),
                                 error = function(e) conditionMessage(e))
  bm_and_inf <- sf::st_sf(geometry = c(sf::st_geometry(bm), sf::st_geometry(bm_inf)))
  expect_identical(err_of(bm_and_inf), err_of(bm))
  expect_match(err_of(bm_inf), "`boundary` has no polygon with a location", fixed = TRUE)
})

test_that("a boundary without a CRS and with no located coordinate is refused by name", {
  ll <- sf::st_crs(4326)
  sq <- rbind(c(Inf, 0), c(-Inf, 0), c(Inf, 1), c(Inf, 0))
  # No rows, every coordinate infinite, every coordinate finite and beyond
  # 1e150: the advice fits each.
  far <- rbind(c(1e200, 1e200), c(2e200, 1e200), c(2e200, 2e200), c(1e200, 1e200))
  for (b in list(ring(rbind(c(-94, 44), c(-92, 44), c(-92, 46), c(-94, 44)))[0, ], ring(sq), ring(far))) {
    e <- tryCatch(.crsless_boundary_as_lonlat(b, ll, "probe"), error = function(e) e)
    expect_s3_class(e, "error")
    expect_null(e$call)
    expect_identical(conditionMessage(e), paste0(
      "probe(): `boundary` has no CRS and no coordinate that is finite and within ",
      "1e150, so it cannot be read as lon/lat beside `points_sf`, which is in a ",
      "geographic CRS (EPSG:4326). Pass a boundary with coordinates that are finite ",
      "and within 1e150, or set its CRS with sf::st_crs()."))
    e <- tryCatch(.crsless_boundary_as_lonlat(b, ll, "probe", assumed = TRUE), error = function(e) e)
    expect_match(conditionMessage(e), "which has no CRS either and was taken as lon/lat (EPSG:4326)",
                 fixed = TRUE)
  }
})

test_that("public routes read a layer without a CRS with an infinite row as the layer with a CRS", {
  set.seed(40)
  d <- data.frame(x = stats::runif(40, -93.6, -92.6), y = stats::runif(40, 44.6, 45.4),
                  z = stats::rnorm(40), a = stats::rnorm(40))
  d$x[7L] <- Inf
  p   <- pts_of(d)
  ref <- with_warnings(prep_model_data(sf::st_set_crs(p, 4326), "z", "a"))
  got <- with_warnings(prep_model_data(p, "z", "a"))
  expect_identical(sf::st_crs(got$value), sf::st_crs(ref$value))
  expect_identical(nrow(got$value), 39L)
  expect_identical(xy_of(got$value), xy_of(ref$value))
  expect_identical(attr(got$value, "crs_assumed"), "EPSG:4326")
  expect_true(any(grepl("assuming EPSG:4326", got$warnings, fixed = TRUE)))
  # It stopped with GEOS's "encountered NaN/Inf numbers".
  ref <- with_warnings(build_tessellation(sf::st_set_crs(p, 4326), method = "voronoi"))
  got <- with_warnings(build_tessellation(p, method = "voronoi"))
  expect_identical(sf::st_geometry(got$value$cells), sf::st_geometry(ref$value$cells))
  expect_identical(got$value$index, ref$value$index)
  expect_identical(sf::st_crs(got$value$cells), sf::st_crs(32615))
})

# ---- BUG-137: the chooser of a lon/lat layer ---------------------------------

test_that("the chooser reads a lon/lat layer as the same layer with the unlocated row empty", {
  want <- with_warnings(ensure_projected(row1_at(NULL, 4326)))
  expect_identical(sf::st_crs(want$value), sf::st_crs(32617))
  expect_length(want$warnings, 0L)
  for (at in unlocated) {
    lab <- paste(at, collapse = ", ")
    # A longitude of 1e200 stopped with "invalid crs: +proj=laea ... +lon_0=-2.95e190".
    expect_no_error(got <- with_warnings(ensure_projected(row1_at(at, 4326))))
    expect_identical(sf::st_crs(got$value), sf::st_crs(want$value), label = lab)
    expect_identical(attr(got$value, "crs_choice"), attr(want$value, "crs_choice"), label = lab)
    expect_identical(xy_of(got$value)[-1L, ], xy_of(want$value)[-1L, ], label = lab)
    # No bare sf, GDAL or base R condition.
    expect_length(got$warnings, 0L)
    expect_false(any(is.finite(xy_of(got$value)[1L, ])), label = lab)
  }
  # purpose = "area": a latitude of 1e200 moved the standard parallels.
  want <- ensure_projected(row1_at(NULL, 4326), purpose = "area")
  for (at in unlocated) {
    got <- with_warnings(ensure_projected(row1_at(at, 4326), purpose = "area"))
    expect_identical(sf::st_crs(got$value), sf::st_crs(want), label = paste(at, collapse = ", "))
    expect_identical(xy_of(got$value)[-1L, ], xy_of(want)[-1L, ])
    expect_length(got$warnings, 0L)
  }
})

test_that("the centre and the wrap test read the located rows only", {
  x <- row1_at(c(1e200, 45), 4326)
  expect_identical(.lonlat_centre(x), .lonlat_centre(row1_at(NULL, 4326)))
  # A layer across 180 degrees: the equal-area projection of the layer with
  # the row emptied.
  set.seed(9)
  d <- data.frame(x = c(stats::runif(30, 176, 180), stats::runif(30, -180, -176)),
                  y = stats::runif(60, -20, -15))
  wrap <- function(at) {
    p <- pts_of(d, 4326)
    g <- sf::st_geometry(p)
    g[1L] <- sf::st_sfc(if (is.null(at)) sf::st_point() else sf::st_point(at), crs = 4326)
    sf::st_set_geometry(p, g)
  }
  want <- .pick_local_projected_crs(wrap(NULL))
  expect_match(want$crs$input, "+proj=laea", fixed = TRUE)
  for (at in unlocated)
    expect_identical(.pick_local_projected_crs(wrap(at)), want, label = paste(at, collapse = ", "))
  # A polygon layer: the feature with a vertex beyond the bound takes no part.
  sq <- function(x0, y0) rbind(c(x0, y0), c(x0 + 1, y0), c(x0 + 1, y0 + 1), c(x0, y0 + 1), c(x0, y0))
  bad <- sq(-99, 40); bad[3L, 2L] <- 1e200
  polys <- sf::st_sfc(lapply(c(-104, -102, -100), function(x0) sf::st_polygon(list(sq(x0, 40)))),
                      crs = 4326)
  with_bad  <- c(polys, sf::st_sfc(sf::st_polygon(list(bad)), crs = 4326))
  with_none <- c(polys, sf::st_sfc(sf::st_polygon(), crs = 4326))
  expect_identical(.pick_local_projected_crs(with_bad), .pick_local_projected_crs(with_none))
})

test_that("a lon/lat layer whose rows all have a location is chosen for as before", {
  # The chooser reads the layer as it came: .empty_unlocated() returns it.
  x <- row1_at(c(-80.2, 45.5), 4326)
  expect_identical(.empty_unlocated(x), x)
  expect_identical(.pick_local_projected_crs(x)$crs, sf::st_crs(32617))
})

test_that("resolution_profile() runs on a lon/lat layer with a longitude of 1e200", {
  # It stopped with sf's bare "invalid crs".
  skip_on_cran()
  quiet <- function(expr) suppressWarnings(suppressMessages(expr))
  want <- quiet(resolution_profile(row1_at(NULL, 4326), "z", sac = NULL, n_levels = 3))
  got  <- quiet(resolution_profile(row1_at(c(1e200, 45), 4326), "z", sac = NULL, n_levels = 3))
  expect_identical(nrow(got), nrow(want))
  expect_equal(as.data.frame(got), as.data.frame(want))
})

# ---- BUG-137: the rows stay without a location after the move -----------------

# Four half-degree squares near (-92, 45) in degrees, the fourth with vertex 3
# at `at`, or empty when `at` is NULL.
squares_at <- function(at, crs = 4326) {
  sq <- function(x0, y0) rbind(c(x0, y0), c(x0 + 0.5, y0), c(x0 + 0.5, y0 + 0.5), c(x0, y0 + 0.5), c(x0, y0))
  b <- sq(-93.9, 45.3)
  if (!is.null(at)) b[3L, ] <- at
  sf::st_sf(v = 1:4, geometry = sf::st_sfc(
    sf::st_polygon(list(sq(-93.2, 45.1))), sf::st_polygon(list(sq(-92.2, 44.1))),
    sf::st_polygon(list(sq(-91.2, 45.6))),
    if (is.null(at)) sf::st_polygon() else sf::st_polygon(list(b)), crs = crs))
}
lines_at <- function(at, crs = 4326) {
  sf::st_sfc(sf::st_linestring(rbind(c(-93, 45), c(-92.5, 45.2))),
             if (is.null(at)) sf::st_linestring() else
               sf::st_linestring(rbind(c(-93, 44), at, c(-92, 44.6))),
             sf::st_linestring(rbind(c(-91, 45), c(-90.5, 45.2))), crs = crs)
}
far_vertex <- list(far_lon = c(1e200, 45.8), far_lat = c(-93.4, 1e200), inf = c(Inf, 45.8),
                   nfar_lon = c(-1e200, 45.8))

test_that("moved from lon/lat, a polygon or line with a vertex without a location comes back empty", {
  # sf::st_transform() drops the vertex PROJ cannot move and keeps the rest:
  # the square came back as a finite triangle.
  quiet <- function(expr) suppressWarnings(expr)
  for (at in far_vertex) {
    lab <- paste(at, collapse = ", ")
    for (case in list(
      list(make = squares_at, f = function(l) ensure_projected(l)),
      list(make = squares_at, f = function(l) ensure_projected(l, purpose = "area")),
      list(make = squares_at, f = function(l) ensure_projected(l, target_crs = 32615)),
      list(make = lines_at,   f = function(l) ensure_projected(l)),
      list(make = lines_at,   f = function(l) ensure_projected(l, target_crs = 5070)))) {
      got  <- case$f(case$make(at))
      want <- case$f(case$make(NULL))
      g <- sf::st_geometry(got)
      w <- sf::st_geometry(want)
      bad <- if (inherits(got, "sf")) 4L else 2L
      expect_identical(sf::st_crs(got), sf::st_crs(want), label = lab)
      expect_identical(attr(got, "crs_choice"), attr(want, "crs_choice"), label = lab)
      expect_true(sf::st_is_empty(g[bad]), label = lab)
      expect_identical(class(g[[bad]]), class(w[[bad]]), label = lab)
      expect_identical(unclass(g)[-bad], unclass(w)[-bad], label = lab)
    }
    # Without a CRS: taken as EPSG:4326 (the lon/lat test reads the other
    # vertices), projected as the layer with EPSG:4326 set.
    for (f in list(function(l) ensure_projected(l),
                   function(l) ensure_projected(l, target_crs = 32615))) {
      got  <- quiet(f(squares_at(at, sf::NA_crs_)))
      want <- quiet(f(squares_at(at)))
      expect_identical(attr(got, "crs_assumed"), "EPSG:4326", label = lab)
      expect_identical(sf::st_crs(got), sf::st_crs(want), label = lab)
      expect_true(sf::st_is_empty(got)[4L], label = lab)
      expect_identical(sf::st_geometry(got)[1:3], sf::st_geometry(want)[1:3], label = lab)
    }
  }
  # A point row comes back as an empty point.
  got <- ensure_projected(row1_at(c(1e200, 45), 4326))
  expect_identical(sf::st_is_empty(got), c(TRUE, rep(FALSE, 299L)))
  # Rows that all have a location come back as the transform gives them.
  x <- squares_at(c(-93.4, 45.8))
  got <- ensure_projected(x)
  attr(got, "crs_choice") <- NULL
  expect_identical(got, sf::st_transform(x, sf::st_crs(got)))
  expect_identical(ensure_projected(x, target_crs = 32615), sf::st_transform(x, 32615))
})

test_that("a lon/lat layer with no row that has a location says so and falls back", {
  one <- squares_at(c(Inf, 45.8))[4L, ]
  lines <- capture_spatialkit_log(got <- ensure_projected(one))
  expect_identical(sf::st_crs(got), sf::st_crs(3857))
  expect_true(sf::st_is_empty(got))
  expect_true(log_has(lines, "no row has a location"))
  expect_false(log_has(lines, "NaN"))
  lines <- capture_spatialkit_log(got <- ensure_projected(one, purpose = "area"))
  expect_match(sf::st_crs(got)$input, "+proj=eqearth", fixed = TRUE)
  expect_true(log_has(lines, "no row has a location.*Equal Earth"))
  # Every row empty: the same.
  lines <- capture_spatialkit_log(got <- ensure_projected(squares_at(NULL)[4L, ]))
  expect_identical(sf::st_crs(got), sf::st_crs(3857))
  expect_true(log_has(lines, "no row has a location"))
})
