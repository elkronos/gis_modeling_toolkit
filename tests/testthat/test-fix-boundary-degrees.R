# tests/testthat/test-fix-boundary-degrees.R
# ---------------------------------------------------------------------------
# Boundaries in degrees that the seeders and the builders misread.
#
# In longitude and latitude sf refuses to sample a region whose bounding box
# is over 10^4 times its area, so random and k-means seeding stopped with
# sf's bare "sampling box is N times larger than sampling region" on two
# small squares on two continents, which they sample once projected.  Such a
# boundary is now sampled piece by piece, in proportion to the areas of its
# pieces; a piece sf still refuses stops the call with an error that names
# `boundary` and says to transform it to a projected CRS.
#
# A lon/lat boundary with a vertex at a latitude beyond 90 is no place on the
# Earth.  s2 read it as a shape folded over the pole: the seeders placed
# seeds on the far meridian, the grid in degrees laid cells there, and the
# Voronoi builder, clip_target_for() and make_folds() read it in silence or
# stopped with a false message.  Each reader now refuses it with the class
# "spatialkit_cells_degenerate" (site "boundary").
# ---------------------------------------------------------------------------

quiet <- function(expr) suppressWarnings(suppressMessages(expr))
sq <- function(x0, y0, d)
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + d, y0), c(x0 + d, y0 + d), c(x0, y0 + d), c(x0, y0))))
layer <- function(..., crs = 4326) sf::st_sf(geometry = sf::st_sfc(list(...), crs = crs))
# Is every point of `pts` in `region`, on the sphere?
all_in <- function(pts, region)
  all(lengths(spatialkit:::.with_s2(sf::st_intersects(pts, region))) > 0L)
with_s2_setting <- function(on, code) {
  old <- sf::sf_use_s2()
  on.exit(suppressMessages(sf::sf_use_s2(old)), add = TRUE)
  suppressMessages(sf::sf_use_s2(on))
  code
}

# Two 0.05 degree squares on two continents, and three spread over 100
# degrees of longitude and 20 of latitude: each is valid under s2, and the
# box of each boundary is over 10^4 times its area (2e5 and 1.3e5 times).
far2 <- layer(sq(-94, 39, 0.05), sq(-60, 70, 0.05))
far3 <- layer(sq(-50, 10, 0.1), sq(0, 20, 0.05), sq(50, 30, 0.05))

test_that("a lon/lat boundary of small pieces far apart is sampled piece by piece", {
  for (on in c(TRUE, FALSE)) with_s2_setting(on, {
    for (b in list(far2, far3)) {
      lab <- sprintf("s2 %s, %d pieces", on, nrow(b))
      s <- quiet(voronoi_seeds_random(b, 10, set_seed = 1))
      expect_identical(nrow(s), 10L, label = lab)
      expect_true(all_in(s, b), label = lab)
      expect_identical(quiet(voronoi_seeds_random(b, 10, set_seed = 1)), s, label = lab)
      r <- quiet(get_voronoi_seeds(b, n = 10, method = "random", set_seed = 1))
      expect_identical(nrow(r), 10L, label = lab)
      expect_true(all_in(r, b), label = lab)
      expect_identical(quiet(get_voronoi_seeds(b, n = 10, method = "random", set_seed = 1)),
                       r, label = lab)
      k <- quiet(get_voronoi_seeds(b, n = 6, method = "kmeans", set_seed = 1))
      expect_identical(nrow(k), 6L, label = lab)
      expect_true(all_in(k, b), label = lab)
      expect_identical(quiet(get_voronoi_seeds(b, n = 6, method = "kmeans", set_seed = 1)),
                       k, label = lab)
    }
  })
})

test_that("the pieces get seeds in proportion to their areas", {
  # The first square has about four times the area of each of the others.
  s <- quiet(voronoi_seeds_random(far3, 3000, set_seed = 1))
  piece <- vapply(spatialkit:::.with_s2(sf::st_intersects(s, far3)), `[`, integer(1), 1L)
  n_in <- tabulate(piece, 3L)
  a <- spatialkit:::.with_s2(as.numeric(sf::st_area(far3)))
  expect_identical(sum(n_in), 3000L)
  # A chi-squared test of the counts against the shares of the areas, at a
  # fixed seed; the counts are random, so a p-value is the measure of fit.
  expect_gt(stats::chisq.test(n_in, p = a / sum(a))$p.value, 0.01)
})

test_that("a boundary sf refuses even piece by piece stops in the caller's name", {
  # One thin row whose own box is about 2e5 times its area.
  thin <- layer(sf::st_polygon(list(rbind(c(0, 0), c(10, 10), c(10, 10.0001), c(0, 0)))))
  expect_true(spatialkit:::.with_s2(sf::st_is_valid(thin)))
  for (f in list(
    function() voronoi_seeds_random(thin, 5, set_seed = 1),
    function() get_voronoi_seeds(thin, n = 5, method = "random", set_seed = 1),
    function() get_voronoi_seeds(thin, n = 5, method = "kmeans", set_seed = 1))) {
    e <- tryCatch(quiet(f()), error = function(e) e)
    expect_s3_class(e, "error")
    expect_match(conditionMessage(e),
                 "^(voronoi_seeds_random|get_voronoi_seeds)\\(\\): sf refuses to sample `boundary`")
    expect_match(conditionMessage(e), "over 10,000 times its area", fixed = TRUE)
    expect_match(conditionMessage(e), "Transform `boundary` to a projected CRS", fixed = TRUE)
  }
})

test_that("a boundary sf samples whole never takes the piece-by-piece path", {
  # The boundaries sampled before keep their seeds: the new path runs only
  # where sf refused, so it is mocked to fail here.
  testthat::local_mocked_bindings(
    .sample_by_piece = function(...) stop("the piece-by-piece path ran"),
    .package = "spatialkit")
  hawaii <- layer(sq(-160, 22, 0.05), sq(-157, 21, 0.05), sq(-155.5, 19.5, 0.05))
  near   <- layer(sq(10, 45, 0.05), sq(11, 45, 0.05))
  proj   <- sf::st_transform(far2, 3857)
  for (b in list(hawaii, near, proj)) {
    expect_identical(nrow(quiet(voronoi_seeds_random(b, 10, set_seed = 1))), 10L)
    expect_identical(nrow(quiet(get_voronoi_seeds(b, n = 10, method = "random",
                                                  set_seed = 1))), 10L)
  }
  # The k-means cloud (2000 points) on the lon/lat boundaries only: sf sizes
  # a projected draw by the ratio of the box to the area, which for `proj`
  # asks for billions of points.
  for (b in list(hawaii, near))
    expect_identical(nrow(quiet(get_voronoi_seeds(b, n = 5, method = "kmeans",
                                                  set_seed = 1))), 5L)
})

# A row past latitude 90, alone and beside a valid row, and points inside the
# part of it below the pole.
past_pole <- sf::st_polygon(list(rbind(c(-94, 89.5), c(-93.8, 89.5), c(-93.8, 95),
                                       c(-94, 95), c(-94, 89.5))))
beyond_alone  <- layer(past_pole)
beyond_beside <- layer(sq(-94.5, 39.5, 0.4), past_pole)
pts_pole <- sf::st_sf(geometry = sf::st_sfc(sf::st_point(c(-93.9, 89.6)),
                                            sf::st_point(c(-93.85, 89.8)),
                                            sf::st_point(c(-93.95, 89.7)), crs = 4326))
calls_on <- function(b) list(
  vsr       = function() voronoi_seeds_random(b, 5, set_seed = 1),
  random    = function() get_voronoi_seeds(b, n = 5, method = "random", set_seed = 1),
  kmeans    = function() get_voronoi_seeds(b, n = 5, method = "kmeans", set_seed = 1),
  grid_deg  = function() create_grid_polygons(b, cellsize = 0.05, crs = 4326, quiet = TRUE),
  grid_open = function() create_grid_polygons(b, cellsize = 0.05, crs = 4326, clip = FALSE,
                                              quiet = TRUE),
  grid_proj = function() create_grid_polygons(b, target_cells = 10, quiet = TRUE),
  cached    = function() create_grid_polygons_cached(b, target_cells = 10, quiet = TRUE),
  voronoi   = function() create_voronoi_polygons(pts_pole, boundary = b, quiet = TRUE),
  bt_vor    = function() build_tessellation(pts_pole, boundary = b, method = "voronoi",
                                            quiet = TRUE),
  bt_sq     = function() build_tessellation(pts_pole, boundary = b, method = "square",
                                            approx_n_cells = 5, quiet = TRUE),
  bt_hex    = function() build_tessellation(pts_pole, boundary = b, method = "hex",
                                            approx_n_cells = 5, quiet = TRUE),
  clip      = function() clip_target_for(pts_pole, boundary = b, quiet = TRUE),
  folds     = function() make_folds(pts_pole, method = "block_kfold", k = 2, boundary = b,
                                    seed = 1))

test_that("a lon/lat boundary with a latitude beyond 90 is refused by every reader", {
  for (on in c(TRUE, FALSE)) with_s2_setting(on, {
    for (b in list(beyond_alone, beyond_beside)) {
      row <- nrow(b)
      calls <- calls_on(b)
      for (nm in names(calls)) {
        lab <- sprintf("%s, s2 %s, %d rows", nm, on, row)
        e <- tryCatch(quiet(calls[[nm]]()), error = function(e) e)
        expect_s3_class(e, "spatialkit_cells_degenerate")
        expect_identical(e$site, "boundary", label = lab)
        expect_identical(e$rows, row, label = lab)
        expect_identical(e$n, row, label = lab)
        expect_match(conditionMessage(e), sprintf(
          "row %d of the %d polygons? of `boundary` has a vertex at a latitude beyond 90",
          row, row), label = lab)
        expect_match(conditionMessage(e), "latitude 95", fixed = TRUE, label = lab)
        expect_match(conditionMessage(e), "check the CRS of `boundary`", fixed = TRUE,
                     label = lab)
      }
    }
  })
})

test_that("a lon/lat boundary that reaches latitude 90 and not beyond is read as given", {
  at_pole <- layer(sf::st_polygon(list(rbind(c(-94, 89.5), c(-93.8, 89.5), c(-93.8, 90),
                                             c(-94, 90), c(-94, 89.5)))))
  s <- quiet(voronoi_seeds_random(at_pole, 5, set_seed = 1))
  expect_true(all(sf::st_coordinates(s)[, "X"] > -94 & sf::st_coordinates(s)[, "X"] < -93.8))
  expect_true(all_in(s, at_pole))
  g <- quiet(create_grid_polygons(at_pole, cellsize = 0.05, crs = 4326))
  cx <- sf::st_coordinates(quiet(sf::st_centroid(sf::st_geometry(g))))[, "X"]
  expect_true(all(cx > -94 & cx < -93.8))
  calls <- calls_on(at_pole)
  for (nm in names(calls))
    expect_false(inherits(tryCatch(quiet(calls[[nm]]()), error = function(e) e), "error"),
                 label = nm)
})

test_that("a boundary without a CRS stamped with a geographic `crs` is tested for its latitude", {
  # The coordinates as given are those of `past_pole`; stamped with
  # EPSG:4326 by the call, they were read folded over the pole: 320 cells in
  # degrees from -160 to 20 E, and "must be polygonal" from the square grid.
  bna <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(20, 89), c(20.2, 89), c(20.2, 93), c(20, 93), c(20, 89))))))
  pts <- sf::st_sf(geometry = sf::st_sfc(sf::st_point(c(20.05, 89.3)),
                                         sf::st_point(c(20.1, 89.5)),
                                         sf::st_point(c(20.15, 89.6)), crs = 4326))
  calls <- list(
    grid_deg  = function() create_grid_polygons(bna, cellsize = 0.05, crs = 4326, quiet = TRUE),
    grid_open = function() create_grid_polygons(bna, cellsize = 0.05, crs = 4326, clip = FALSE,
                                                quiet = TRUE),
    grid_n    = function() create_grid_polygons(bna, target_cells = 10, crs = 4326, quiet = TRUE),
    bt_sq     = function() build_tessellation(pts, boundary = bna, method = "square",
                                              approx_n_cells = 5, crs = 4326, quiet = TRUE),
    bt_hex    = function() build_tessellation(pts, boundary = bna, method = "hex",
                                              approx_n_cells = 5, crs = 4326, quiet = TRUE),
    bt_vor    = function() build_tessellation(pts, boundary = bna, method = "voronoi",
                                              crs = 4326, quiet = TRUE))
  for (nm in names(calls)) {
    e <- tryCatch({ quiet(calls[[nm]]()); NULL }, error = function(e) e)
    expect_true(inherits(e, "spatialkit_cells_degenerate"), label = nm)
    expect_identical(e$cause, "latitude", label = nm)
    expect_identical(e$site, "boundary", label = nm)
    expect_match(if (is.null(e)) "<returned>" else conditionMessage(e),
                 "(row 1: latitude 93)", fixed = TRUE, label = nm)
  }
})

test_that("a latitude just beyond 90 is printed beyond 90", {
  # Rounded to six digits, 90 + 1e-9 read "latitude 90" in an error that
  # says the vertex lies beyond 90.
  for (top in c(90 + 1e-9, 90 + 2e-14, -(90 + 1e-9))) {
    ring <- if (top > 0) rbind(c(20, 89), c(20.2, 89), c(20.2, top), c(20, 90), c(20, 89))
            else rbind(c(20, -89), c(20, -90), c(20.2, top), c(20.2, -89), c(20, -89))
    b <- layer(sf::st_polygon(list(ring)))
    e <- tryCatch(quiet(voronoi_seeds_random(b, 5, set_seed = 1)), error = function(e) e)
    expect_s3_class(e, "spatialkit_cells_degenerate")
    shown <- as.numeric(sub("^latitude ", "", e$reasons))
    expect_gt(abs(shown), 90)
    # The fewest digits that show it beyond 90, so the printed value is the
    # latitude to within the last digit shown.
    expect_equal(shown, top, tolerance = 1e-12)
  }
})

test_that("an sf error other than the box refusal keeps its own path", {
  # Only sf's refusal of a box over 10^4 times the area is caught and drawn
  # piece by piece; any other error of sf::st_sample() is raised as it is.
  testthat::local_mocked_bindings(
    st_sample = function(...) stop("another sf error"), .package = "sf")
  b <- layer(sq(10, 45, 0.5))
  for (f in list(
    function() voronoi_seeds_random(b, 5, set_seed = 1),
    function() get_voronoi_seeds(b, n = 5, method = "random", set_seed = 1))) {
    e <- tryCatch(quiet(f()), error = function(e) e)
    expect_identical(conditionMessage(e), "another sf error")
  }
})
