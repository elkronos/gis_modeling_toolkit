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
  # where sf refused, or where sf's draw in the box is above 10^6 points
  # (test-fix-sample-draw-size.R), so it is mocked to fail here.
  testthat::local_mocked_bindings(
    .sample_by_piece = function(...) stop("the piece-by-piece path ran"),
    .package = "spatialkit")
  hawaii <- layer(sq(-160, 22, 0.05), sq(-157, 21, 0.05), sq(-155.5, 19.5, 0.05))
  near   <- layer(sq(10, 45, 0.05), sq(11, 45, 0.05))
  # The box of `proj` is 1.8e5 times its area: 5 seeds take 9.2e5 points in
  # the box, under the cap, and 10 seeds take 1.8e6, above it.
  proj   <- sf::st_transform(far2, 3857)
  for (b in list(hawaii, near)) {
    expect_identical(nrow(quiet(voronoi_seeds_random(b, 10, set_seed = 1))), 10L)
    expect_identical(nrow(quiet(get_voronoi_seeds(b, n = 10, method = "random",
                                                  set_seed = 1))), 10L)
  }
  expect_identical(nrow(quiet(voronoi_seeds_random(proj, 5, set_seed = 1))), 5L)
  expect_identical(nrow(quiet(get_voronoi_seeds(proj, n = 5, method = "random",
                                                set_seed = 1))), 5L)
  expect_error(quiet(voronoi_seeds_random(proj, 10, set_seed = 1)),
               "the piece-by-piece path ran", fixed = TRUE)
  # The k-means cloud (2000 points): on `near` it is drawn whole.  The box
  # of `hawaii` is 1548 times its area, so its cloud takes 3.1e6 points in
  # the box, and the cloud on `proj` 3.7e8: both are above the cap and go
  # piece by piece.
  expect_identical(nrow(quiet(get_voronoi_seeds(near, n = 5, method = "kmeans",
                                                set_seed = 1))), 5L)
  for (b in list(hawaii, proj))
    expect_error(quiet(get_voronoi_seeds(b, n = 5, method = "kmeans", set_seed = 1)),
                 "the piece-by-piece path ran", fixed = TRUE)
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

# ---------------------------------------------------------------------------
# A lon/lat boundary whose bounding box has no height: every vertex on one
# parallel.  s2 reads its edges as great circles, so the ring holds ground,
# and the seeders refuse it, since sf draws in its box.  Projected, the
# parallel is a curve and the ring has an area, so the refusal says to
# transform `boundary` to a projected CRS, a polar one for a ring around a
# pole.  The tests seed each boundary so transformed, so the advice holds.
# ---------------------------------------------------------------------------

ring_layer <- function(xy) layer(sf::st_polygon(list(xy)))
# Three vertices from 10 to 12 E at 60 N, and rings around each pole.
on_60n  <- ring_layer(rbind(c(10, 60), c(11, 60), c(12, 60), c(10, 60)))
north87 <- ring_layer(rbind(c(0, 87), c(90, 87), c(180, 87), c(-90, 87), c(0, 87)))
south87 <- ring_layer(rbind(c(0, -87), c(-90, -87), c(180, -87), c(90, -87), c(0, -87)))
projected_remedy <- paste0(
  "Transform `boundary` to a projected CRS in which a parallel is a curve, such as the ",
  "one ensure_projected() picks for it (for a ring around a pole, a polar one such as ",
  "EPSG:3413 or 3031), to seed it.")
seeders_on <- function(b) list(
  vsr    = function() voronoi_seeds_random(b, 5, set_seed = 1),
  random = function() get_voronoi_seeds(b, n = 5, method = "random", set_seed = 1),
  kmeans = function() get_voronoi_seeds(b, n = 5, method = "kmeans", set_seed = 1))
err_msg <- function(f) {
  e <- tryCatch({ quiet(f()); NULL }, error = function(e) e)
  if (inherits(e, "error")) conditionMessage(e) else "<returned>"
}

test_that("the no-height refusal of a lon/lat ring on a parallel says to project it", {
  # With s2 off as with it on: get_voronoi_seeds() repaired the rows in the
  # plane with s2 off, which turned these rings into lines, and refused them
  # as rows s2 refuses.
  for (on in c(TRUE, FALSE)) with_s2_setting(on, {
    for (nm in c("on_60n", "north87", "south87")) {
      b <- get(nm)
      expect_true(spatialkit:::.with_s2(sf::st_is_valid(b)), label = nm)
      calls <- seeders_on(b)
      for (f in names(calls)) {
        lab <- sprintf("%s, %s, s2 %s", nm, f, on)
        msg <- err_msg(calls[[f]])
        expect_match(msg, paste0(
          "^(voronoi_seeds_random|get_voronoi_seeds)\\(\\): `boundary` cannot be sampled ",
          "\\(its bounding box, read from its coordinates, has no height\\), so no ",
          "(seed|point of the k-means cloud) can be drawn in it\\. "), label = lab)
        expect_true(endsWith(msg, projected_remedy), label = lab)
      }
    }
  })
})

test_that("the projected CRS the refusal names seeds the same boundary", {
  # The CRS clip_target_for() picks for the ring at 60 N (a UTM zone), the
  # polar stereographic CRS of each pole, and the CRS ensure_projected()
  # picks for each ring, which the refusal names.
  pts <- sf::st_sf(geometry = sf::st_sfc(sf::st_point(c(11, 60.001)), crs = 4326))
  utm <- sf::st_crs(quiet(clip_target_for(pts, boundary = on_60n, quiet = TRUE)))
  expect_false(sf::st_is_longlat(utm))
  picked <- function(b) sf::st_crs(quiet(ensure_projected(b)))
  for (case in list(list(on_60n, utm), list(north87, sf::st_crs(3413)),
                    list(south87, sf::st_crs(3031)), list(on_60n, picked(on_60n)),
                    list(north87, picked(north87)), list(south87, picked(south87)))) {
    p <- sf::st_transform(case[[1]], case[[2]])
    bb <- sf::st_bbox(p)
    expect_gt(as.numeric(bb[["ymax"]] - bb[["ymin"]]), 0)
    expect_gt(as.numeric(sf::st_area(p)), 0)
    calls <- seeders_on(p)
    for (f in names(calls)) {
      s <- quiet(calls[[f]]())
      expect_identical(nrow(s), 5L, label = f)
      # Seeds fall in the ring, allowing 1 cm for the snap of a k-means
      # centre onto its edge.
      d <- as.numeric(sf::st_distance(s, p))
      expect_true(all(d < 0.01), label = f)
    }
  }
})

test_that("a box with no height that projection cannot fix offers no remedy", {
  # The equator and a meridian are great circles: a ring on one holds no
  # ground on the sphere, and projected it is still a line.  A projected
  # region with no height has no area.
  on_equator  <- ring_layer(rbind(c(10, 0), c(11, 0), c(12, 0), c(10, 0)))
  on_meridian <- ring_layer(rbind(c(10, 60), c(10, 61), c(10, 62), c(10, 60)))
  for (b in list(on_equator, on_meridian)) {
    msg <- err_msg(function() voronoi_seeds_random(b, 5, set_seed = 1))
    expect_match(msg, "^voronoi_seeds_random\\(\\): `boundary` cannot be sampled .*\\.$")
    expect_false(grepl("projected CRS", msg, fixed = TRUE))
  }
  flat <- sf::st_sfc(sf::st_polygon(list(rbind(c(5e5, 5e6), c(6e5, 5e6), c(7e5, 5e6),
                                               c(5e5, 5e6)))), crs = 32617)
  # A projected box with no height at a northing within 90 of 0, which reads
  # as a latitude other than 0 and 90 (the CRS, not the number, decides),
  # a box that is one point, and a ring with every vertex on a pole.
  near_zero <- sf::st_sfc(sf::st_polygon(list(rbind(c(5e5, 50), c(6e5, 50), c(7e5, 50),
                                                    c(5e5, 50)))), crs = 32631)
  one_point <- sf::st_sfc(sf::st_polygon(list(rbind(c(10, 60), c(10, 60), c(10, 60),
                                                    c(10, 60)))), crs = 4326)
  on_pole <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 90), c(90, 90), c(180, 90),
                                                  c(0, 90)))), crs = 4326)
  for (g in list(flat, near_zero, one_point, on_pole)) {
    e <- tryCatch(spatialkit:::.sample_box_check(g, "probe", "no seed can be drawn in it"),
                  error = function(e) e)
    expect_identical(conditionMessage(e), paste0(
      "probe(): `boundary` cannot be sampled (its bounding box, read from its coordinates, ",
      "has no height), so no seed can be drawn in it."))
  }
  # In Web Mercator, World Mercator and plate carree a parallel stays
  # straight, so the ring on it still has no height: the remedy names a CRS
  # in which a parallel is a curve, as the one ensure_projected() picks.
  for (crs in c(3857, 3395, 4087)) {
    bb <- sf::st_bbox(sf::st_transform(on_60n, crs))
    expect_equal(as.numeric(bb[["ymax"]] - bb[["ymin"]]), 0, label = crs)
  }
})

test_that("a polygon with no height among others is named as one of the polygons", {
  # A small square in Kansas beside the ring at 60 N: the box of the whole
  # is over 10^4 times its area, so it is sampled piece by piece, and the
  # ring is the piece whose box has no height.  The box of `boundary` has a
  # height.
  b <- layer(sq(-94, 39, 0.05), sf::st_polygon(list(rbind(c(10, 60), c(11, 60), c(12, 60),
                                                         c(10, 60)))))
  bb <- sf::st_bbox(b)
  expect_gt(as.numeric(bb[["ymax"]] - bb[["ymin"]]), 0)
  calls <- seeders_on(b)
  for (f in names(calls)) {
    msg <- err_msg(calls[[f]])
    expect_match(msg, paste0(
      "^(voronoi_seeds_random|get_voronoi_seeds)\\(\\): one of the polygons of `boundary` ",
      "cannot be sampled \\(its bounding box, read from its coordinates, has no height\\)"),
      label = f)
    expect_true(endsWith(msg, projected_remedy), label = f)
  }
  # Projected, both pieces have an area and the boundary is seeded.
  p <- sf::st_transform(b, sf::st_crs(quiet(ensure_projected(b))))
  for (f in names(seeders_on(p)))
    expect_identical(nrow(quiet(seeders_on(p)[[f]]())), 5L, label = f)
})

test_that("the help of the seeders states the remedy of the no-height refusal", {
  root <- testthat::test_path("..", "..")
  flat <- function(x)
    paste(unlist(lapply(x, function(el) if (is.character(el)) el else flat(el))), collapse = "")
  page_text <- function(page) {
    file <- file.path(root, "man", paste0(page, ".Rd"))
    rd <- if (file.exists(file)) tools::parse_Rd(file)
          else tools::Rd_db("spatialkit")[[paste0(page, ".Rd")]]
    expect_s3_class(rd, "Rd")
    gsub("[[:space:]]+", " ", flat(rd))
  }
  # The rules stand once, in "Sampling inside the boundary" of
  # ?get_voronoi_seeds; ?voronoi_seeds_random draws its seeds as that section
  # describes and points to it.
  expect_true(grepl(paste0(
    "A boundary whose bounding box has no height or no width (three vertices on one ",
    "parallel) is refused. In longitude and latitude s2 reads the edges as great ",
    "circles, which leave the bounding box read from the coordinates, and sf (without ",
    "lwgeom) draws in that box: a boundary whose box has less than half its ",
    "area on the sphere (three vertices on or near one parallel) is refused too. That ",
    "error, and the error for a box with no height whose vertices lie on one parallel ",
    "other than the equator, say to transform boundary to a projected CRS in which a ",
    "parallel is a curve, such as the one ensure_projected() picks for it (a polar one, ",
    "such as EPSG:3413 or 3031, for a ring around a pole), where the region has an area ",
    "and lies in its box."), page_text("get_voronoi_seeds"), fixed = TRUE),
    label = "get_voronoi_seeds")
  expect_true(grepl(paste0(
    "The seeds are drawn as \"Sampling inside the boundary\" in get_voronoi_seeds() ",
    "describes"), page_text("voronoi_seeds_random"), fixed = TRUE),
    label = "voronoi_seeds_random")
})

# ---------------------------------------------------------------------------
# A lon/lat region that lies mostly outside the box read from its
# coordinates.  s2 reads the edges as great circles, which leave the box,
# and sf (without lwgeom) draws in the box and keeps the points in the
# region.  The s2 union of a ring on the 30th parallel, or on the 45th at
# 100 W, moves its latitudes by a few units in the last place, so its box
# has a height of about 1e-14 degrees and passed the test of a box with no
# height, and sf drew without end; so did a sliver 1e-7 degrees tall.
# Which latitudes passed differed between macOS and Linux.  sf's draw is
# mocked to stop here, so a regression fails and does not hang.
# ---------------------------------------------------------------------------

on_30n   <- ring_layer(rbind(c(10, 30), c(11, 30), c(12, 30), c(10, 30)))
on_45n_w <- ring_layer(rbind(c(-100, 45), c(-99, 45), c(-98, 45), c(-100, 45)))
sliver_7 <- ring_layer(rbind(c(10, 60), c(11, 60.0000001), c(12, 60), c(10, 60)))
# A triangle 120 degrees wide at 80 N: its box has 0.42 times its area, so
# sf drew in less than half of it.
tri_80n  <- ring_layer(rbind(c(0, 80), c(60, 79), c(120, 80), c(0, 80)))
# sf's draw, mocked: a region whose box holds half its area or more, or
# that sf refuses before any draw (its box is over 10^4 times its area), goes
# to sf; any other stops here, where sf drew without end.
mocked_draw <- function(x, size, ...) {
  plan <- spatialkit:::.sample_draw_plan(x, size)
  if (any(plan$sf_refuses) || all(plan$ratio >= 0.5))
    return(suppressWarnings(sf::st_sample(x, size = size, ...)))
  stop("sf::st_sample() was reached on a region outside its box", call. = FALSE)
}
short_box <- paste0(
  "^(voronoi_seeds_random|get_voronoi_seeds)\\(\\): `boundary` cannot be sampled ",
  "\\(its edges are great circles that leave its bounding box, read from its ",
  "coordinates, which sf draws in: the box has [0-9.e+-]+ times its area on the ",
  "sphere\\), so no (seed|point of the k-means cloud) can be drawn in it\\. ")

test_that("a lon/lat region whose box has less than half its area is refused at once", {
  # The test of the plan, on the region the seeders sample, first.
  for (nm in c("on_30n", "on_45n_w", "sliver_7", "tri_80n")) {
    u <- quiet(spatialkit:::.boundary_union(get(nm), caller = "probe", task = "t",
                                            all = TRUE, nothing = "x"))
    plan <- spatialkit:::.sample_draw_plan(u, 5)
    expect_lt(plan$ratio, 0.5, label = nm)
    e <- tryCatch(spatialkit:::.sample_box_short(plan, "probe", "x"), error = function(e) e)
    expect_s3_class(e, "error")
  }
  testthat::local_mocked_bindings(.st_sample_quiet = mocked_draw)
  for (on in c(TRUE, FALSE)) with_s2_setting(on, {
    for (nm in c("on_30n", "on_45n_w", "sliver_7", "tri_80n")) {
      calls <- seeders_on(get(nm))
      for (f in names(calls)) {
        lab <- sprintf("%s, %s, s2 %s", nm, f, on)
        msg <- err_msg(calls[[f]])
        # Either refusal: the union of some of these rings keeps its
        # latitudes exact, and its box has no height.
        expect_true(grepl(short_box, msg) || grepl("has no height\\)", msg), label = lab)
        expect_true(endsWith(msg, projected_remedy), label = lab)
      }
    }
  })
})

test_that("the projected CRS ensure_projected() picks seeds those regions", {
  for (nm in c("on_30n", "on_45n_w", "sliver_7", "tri_80n")) {
    b <- get(nm)
    p <- sf::st_transform(b, sf::st_crs(quiet(ensure_projected(b))))
    calls <- seeders_on(p)
    for (f in names(calls)) {
      s <- quiet(calls[[f]]())
      expect_identical(nrow(s), 5L, label = paste(nm, f))
      expect_true(all(as.numeric(sf::st_distance(s, p)) < 0.01), label = paste(nm, f))
    }
  }
})

test_that("a piece whose box has less than half its area is named as one of the polygons", {
  # The sliver beside a small square in Kansas: the box of the whole is over
  # 10^4 times its area, so it is sampled piece by piece.
  b <- layer(sq(-94, 39, 0.05), sf::st_polygon(list(rbind(
    c(10, 60), c(11, 60.0000001), c(12, 60), c(10, 60)))))
  testthat::local_mocked_bindings(.st_sample_quiet = mocked_draw)
  calls <- seeders_on(b)
  for (f in names(calls)) {
    msg <- err_msg(calls[[f]])
    expect_match(msg, paste0(
      "^(voronoi_seeds_random|get_voronoi_seeds)\\(\\): one of the polygons of `boundary` ",
      "cannot be sampled \\(its edges are great circles that leave its bounding box"), label = f)
    expect_true(endsWith(msg, projected_remedy), label = f)
  }
})

test_that("with lwgeom, which sf draws in the box of great circles with, nothing is refused", {
  u <- quiet(spatialkit:::.boundary_union(tri_80n, caller = "probe", task = "t",
                                          nothing = "x"))
  plan <- spatialkit:::.sample_draw_plan(u, 5)
  testthat::local_mocked_bindings(.sf_box_from_coordinates = function() FALSE)
  expect_null(spatialkit:::.sample_box_short(plan, "probe", "x"))
  # A projected box is never read so.
  expect_null(spatialkit:::.sample_box_short(list(longlat = FALSE, ratio = 0.1), "probe", "x"))
})

test_that("a region whose box holds half of it or more is drawn as before", {
  # A 0.01 degree cell, a square of one degree at 60 N, and 60 degrees from
  # the equator to 60 N (its box has 0.98 times its area).
  for (g in list(sq(10, 60, 0.01), sq(10, 60, 1),
                 sf::st_polygon(list(rbind(c(0, 0), c(60, 0), c(60, 60), c(0, 60), c(0, 0)))))) {
    b <- layer(g)
    plan <- spatialkit:::.sample_draw_plan(sf::st_geometry(b), 5)
    expect_gte(plan$ratio, 0.5)
    expect_null(spatialkit:::.sample_box_short(plan, "probe", "x"))
    s <- quiet(voronoi_seeds_random(b, 5, set_seed = 1))
    expect_identical(nrow(s), 5L)
    expect_true(all_in(s, b))
  }
})

test_that("a lon/lat ring on one meridian, whose union is a line, is refused", {
  # Its union on the sphere is a MULTILINESTRING whose box has a width of
  # rounding size, which sf sampled as lines: voronoi_seeds_random()
  # returned MULTIPOINT seeds, two of them empty.
  on_meridian <- ring_layer(rbind(c(10, -80), c(10, 0), c(10, 80), c(10, -80)))
  testthat::local_mocked_bindings(.st_sample_quiet = mocked_draw)
  for (on in c(TRUE, FALSE)) with_s2_setting(on, {
    expect_identical(
      err_msg(function() voronoi_seeds_random(on_meridian, 5, set_seed = 1)),
      paste0("voronoi_seeds_random(): `boundary` has no polygon with a location (the union ",
             "of its rows leaves no polygon with an area), so no seed can be drawn in it."))
    # get_voronoi_seeds() repairs every row, which s2 does not refuse here.
    expect_identical(
      err_msg(function() get_voronoi_seeds(on_meridian, n = 5, method = "random", set_seed = 1)),
      paste0("get_voronoi_seeds(): `boundary` has no polygon with a location (the repair ",
             "of its rows leaves no polygon with an area), so no seed can be drawn in it."))
  })
})
