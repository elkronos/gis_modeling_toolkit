# tests/testthat/test-fix-tessellation-msgs.R
# ---------------------------------------------------------------------------
# Messages about cells that named the wrong function or left out a size.
#
# BUG-167: create_voronoi_polygons(), called directly, counted the points
# and seeds in no cell, and logged the points it assigned to the nearest
# cell, in the name of build_tessellation().  They now name the function
# called.
# BUG-264: clip_target_for() with no boundary, on points at one location,
# said "degenerate bbox; using small buffer around points." for a disk of
# radius 1 unit of the CRS.  The message now gives that radius in the
# layer's unit and the ways to another size; the target is unchanged.
# BUG-271: points closer together than about 3e-153 were buffered by a
# share of their spread whose square is not a normal double, and GEOS gave
# a MULTIPOINT, an empty LINESTRING or a disk of subnormal area; a box
# 1e-320 on a side came back as a POINT.  Such points now get the disk of
# radius 1 unit of the CRS, with the message that gives their spread.
# ---------------------------------------------------------------------------

.tm_box <- function(w) sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(
  c(xmin = 5e5, ymin = 4e6, xmax = 5e5 + w, ymax = 4e6 + w), crs = sf::st_crs(32617))))

.tm_warnings <- function(expr) {
  w <- list()
  value <- withCallingHandlers(expr, warning = function(c) {
    w[[length(w) + 1L]] <<- c
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = w)
}

.tm_unassigned <- function(r) {
  un <- Filter(function(c) inherits(c, "spatialkit_unassigned"), r$warnings)
  vapply(un, conditionMessage, character(1))
}

test_that("the points in no cell are counted in the name of the function called", {
  withr::local_seed(21)
  p <- sf::st_as_sf(data.frame(x = 5e5 + stats::runif(30, 0, 5000),
                               y = 4e6 + stats::runif(30, 0, 5000)),
                    coords = c("x", "y"), crs = 32617)
  bnd <- .tm_box(2500)
  direct <- .tm_warnings(create_voronoi_polygons(p, boundary = bnd, quiet = TRUE))
  via    <- .tm_warnings(build_tessellation(p, boundary = bnd, method = "voronoi", quiet = TRUE))
  n_out  <- sum(is.na(direct$value$index))
  expect_gt(n_out, 0L)
  head_of <- function(fn) sprintf("%s(): %d of 30 point(s) (%.1f%%) fall in no cell and are NA",
                                  fn, n_out, 100 * n_out / 30)
  expect_length(.tm_unassigned(direct), 1L)
  expect_true(startsWith(.tm_unassigned(direct), head_of("create_voronoi_polygons")))
  expect_length(.tm_unassigned(via), 1L)
  expect_true(startsWith(.tm_unassigned(via), head_of("build_tessellation")))
  # The rest of the text, and the result, are the same on both routes.
  expect_identical(sub("^[a-z_]+\\(\\): ", "", .tm_unassigned(direct)),
                   sub("^[a-z_]+\\(\\): ", "", .tm_unassigned(via)))
  expect_identical(direct$value$index, via$value$index)
})

test_that("the seeds without a cell are counted in the name of the function called", {
  s <- get_voronoi_seeds(method = "provided", seeds = sf::st_as_sf(
    data.frame(x = 5e5 + c(500, 1500, 500, 1500, 6500), y = 4e6 + c(500, 500, 1500, 1500, 1000)),
    coords = c("x", "y"), crs = 32617))
  bnd <- .tm_box(2500)
  direct <- .tm_warnings(create_voronoi_polygons(s, boundary = bnd, quiet = TRUE))
  via    <- .tm_warnings(build_tessellation(s, boundary = bnd, method = "voronoi", quiet = TRUE))
  expect_identical(is.na(direct$value$index), c(FALSE, FALSE, FALSE, FALSE, TRUE))
  expect_identical(direct$value$index, via$value$index)
  expect_length(.tm_unassigned(direct), 1L)
  expect_true(startsWith(.tm_unassigned(direct), paste0(
    "create_voronoi_polygons(): 1 of 5 seed(s) (20.0%) have no cell and are NA in `index`")))
  expect_length(.tm_unassigned(via), 1L)
  expect_true(startsWith(.tm_unassigned(via), paste0(
    "build_tessellation(): 1 of 5 seed(s) (20.0%) have no cell and are NA in `index`")))
})

test_that("the log line on points assigned to the nearest cell names the function called", {
  U <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(0, 0), c(1000, 0), c(1000, 1000), c(550, 1000), c(550, 300),
    c(450, 300), c(450, 1000), c(0, 1000), c(0, 0)))), crs = 32617))
  # The sixth point lies 1e-4 above the notch's floor, outside the boundary
  # and a hair from the fifth point's cell, so it is assigned to that cell.
  p <- sf::st_as_sf(data.frame(x = c(150, 850, 150, 850, 500, 500),
                               y = c(150, 150, 650, 650, 299.99995, 300.0001)),
                    coords = c("x", "y"), crs = 32617)
  line <- "point\\(s\\) sitting just outside every cell"
  direct <- capture_spatialkit_log(
    r1 <- suppressWarnings(create_voronoi_polygons(p, boundary = U, quiet = TRUE)))
  via <- capture_spatialkit_log(
    r2 <- suppressWarnings(build_tessellation(p, boundary = U, method = "voronoi", quiet = TRUE)))
  expect_identical(r1$params$snapped$which, 6L)
  expect_identical(r2$params$snapped$which, 6L)
  expect_true(log_has(direct, paste0("create_voronoi_polygons\\(\\): 1 ", line)))
  expect_false(log_has(direct, paste0("build_tessellation\\(\\): 1 ", line)))
  expect_true(log_has(via, paste0("build_tessellation\\(\\): 1 ", line)))
})

test_that("a cover point assigned to the nearest cell, or in no cell, is named in the function called", {
  # Four seeds in a 1 km square, and a cover of four points: one 1e-4 m
  # outside the square, a hair from a cell, and one 2 km outside it.  The
  # cover is checked against the cells, and its log line and warning named
  # build_tessellation() on a direct call to create_voronoi_polygons(), and
  # create_voronoi_polygons() through build_tessellation(method = "voronoi").
  U <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(0, 0), c(1000, 0), c(1000, 1000), c(0, 1000), c(0, 0)))), crs = 32617))
  s <- get_voronoi_seeds(method = "provided", seeds = sf::st_as_sf(
    data.frame(x = c(250, 750, 250, 750), y = c(250, 250, 750, 750)),
    coords = c("x", "y"), crs = 32617))
  cov <- sf::st_as_sf(data.frame(x = c(100, 500, 1000.0001, 3000), y = c(100, 500, 500, 500)),
                      coords = c("x", "y"), crs = 32617)
  line <- "point\\(s\\) sitting just outside every cell"
  direct <- NULL; via <- NULL
  log_d <- capture_spatialkit_log(
    direct <- .tm_warnings(create_voronoi_polygons(s, boundary = U, cover = cov, quiet = TRUE)))
  log_v <- capture_spatialkit_log(
    via <- .tm_warnings(build_tessellation(s, boundary = U, cover = cov, method = "voronoi",
                                           quiet = TRUE)))
  expect_identical(direct$value$params$cover_outside, 1L)
  expect_identical(via$value$params$cover_outside, 1L)
  expect_true(log_has(log_d, paste0("create_voronoi_polygons\\(\\): 1 ", line)))
  expect_false(log_has(log_d, paste0("build_tessellation\\(\\): 1 ", line)))
  expect_true(log_has(log_v, paste0("build_tessellation\\(\\): 1 ", line)))
  expect_false(log_has(log_v, paste0("create_voronoi_polygons\\(\\): 1 ", line)))
  head_of <- function(fn) paste0(fn, "(): 1 of the 4 `cover` point(s) (25.0%) fall outside every cell")
  expect_length(.tm_unassigned(direct), 1L)
  expect_true(startsWith(.tm_unassigned(direct), head_of("create_voronoi_polygons")))
  expect_length(.tm_unassigned(via), 1L)
  expect_true(startsWith(.tm_unassigned(via), head_of("build_tessellation")))
  expect_identical(sub("^[a-z_]+\\(\\): ", "", .tm_unassigned(direct)),
                   sub("^[a-z_]+\\(\\): ", "", .tm_unassigned(via)))
})

# Two identical points; `crs` as sf reads it.
.tm_one_place <- function(crs, xy) sf::st_as_sf(
  data.frame(x = rep(xy[1], 2), y = rep(xy[2], 2)), coords = c("x", "y"), crs = crs)

.tm_msgs <- function(expr) {
  m <- character()
  value <- withCallingHandlers(expr, message = function(c) {
    m <<- c(m, conditionMessage(c))
    invokeRestart("muffleMessage")
  })
  list(value = value, msgs = m)
}

test_that("clip_target_for() names the radius and its unit for points at one location", {
  km <- "+proj=utm +zone=17 +units=km +datum=WGS84"
  to_km <- "+proj=utm +zone=17 +to_meter=1000 +datum=WGS84"
  cases <- list(
    list(crs = 32617, xy = c(500000, 4e6), unit = "1 metre (1 unit of the CRS)",
         of = "the CRS"),
    list(crs = km, xy = c(500, 4000), unit = "1 kilometre (1 unit of the CRS)",
         of = "the CRS"),
    list(crs = to_km, xy = c(500, 4000), unit = "1 kilometre (1 unit of the CRS)",
         of = "the CRS"),
    # With no CRS the radius is in the units of the coordinates, and the
    # message names no CRS as their unit.
    list(crs = NA, xy = c(500000, 4e6),
         unit = "1 unit of the coordinates (the points have no CRS)", of = "the coordinates"))
  for (cs in cases) {
    pts <- .tm_one_place(cs$crs, cs$xy)
    r <- .tm_msgs(clip_target_for(pts))
    expect_identical(r$msgs, paste0(
      "clip_target_for(): degenerate bbox: the points sit at one location, so the box has no ",
      "extent to take a share of; using a disk of radius ", cs$unit, " around it. Pass ",
      "`boundary` for a target of real size, or an `expand` of 1 or more, which is taken ",
      "as the radius in the units of ", cs$of, ".\n"), info = format(cs$crs))
    # The target is the disk it was: radius 1 in the CRS's unit, the CRS kept.
    g <- sf::st_geometry(r$value)
    expect_identical(sf::st_crs(g), sf::st_crs(pts))
    expect_equal(as.numeric(sf::st_area(g)), pi, tolerance = 1e-3, info = format(cs$crs))
    bb <- sf::st_bbox(g)
    expect_equal(as.numeric(bb[c("xmin", "xmax")]), cs$xy[1] + c(-1, 1), tolerance = 1e-12)
    # A fraction `expand` takes a share of no extent: the same disk and words.
    r2 <- .tm_msgs(clip_target_for(pts, expand = 0.3))
    expect_identical(r2$msgs, r$msgs)
    expect_identical(sf::st_geometry(r2$value), g)
    # An `expand` of 1 or more is the radius, with the message it always had.
    r3 <- .tm_msgs(clip_target_for(pts, expand = 5))
    expect_identical(r3$msgs, "clip_target_for(): degenerate bbox; using small buffer around points.\n")
    expect_equal(as.numeric(sf::st_area(r3$value)), 25 * pi, tolerance = 1e-3)
  }
  # quiet = TRUE says nothing, as before.
  expect_identical(.tm_msgs(clip_target_for(.tm_one_place(32617, c(500000, 4e6)),
                                            quiet = TRUE))$msgs, character())
})

test_that("clip_target_for() says how far apart points are whose spread has no share", {
  # Two points 5e-324 m apart: the box has an extent, but 5 percent of it
  # is 0, so the target is the disk of radius 1 unit of the CRS around them,
  # and the message gives their spread in place of "one location".
  pts <- sf::st_as_sf(data.frame(x = c(0, 5e-324), y = c(0, 0)),
                      coords = c("x", "y"), crs = 32617)
  r <- .tm_msgs(clip_target_for(pts))
  expect_identical(r$msgs, paste0(
    "clip_target_for(): degenerate bbox: the points lie within ", format(5e-324),
    " of each other, so the box has no extent to take a share of; using a disk of ",
    "radius 1 metre (1 unit of the CRS) around them. Pass `boundary` for a target of ",
    "real size, or an `expand` of 1 or more, which is taken as the radius in the units ",
    "of the CRS.\n"))
  expect_equal(as.numeric(sf::st_area(r$value)), pi, tolerance = 1e-3)
})

test_that("clip_target_for() gives the disk of radius 1 to points whose spread squares below doubles", {
  # The share of the spread the buffer would take (5 percent) is below
  # sqrt(.Machine$double.xmin), about 1.5e-154, so its square is not a
  # normal double: the points have no extent to take a share of.  The disk
  # of radius 1 m is drawn at the origin, so its area is pi to the
  # tolerance of GEOS's 30 segments a quadrant, as in the tests above.
  # Each spread is the diagonal of the points' box, formatted as given.
  line <- function(d) sf::st_as_sf(data.frame(x = c(0, d), y = c(0, 0)),
                                   coords = c("x", "y"), crs = 32617)
  box <- function(d) sf::st_as_sf(data.frame(x = c(0, d), y = c(0, d)),
                                  coords = c("x", "y"), crs = 32617)
  cases <- list(list(pts = line(1e-320), spread = 1e-320),
                list(pts = line(1e-310), spread = 1e-310),
                list(pts = line(1e-160), spread = 1e-160),
                list(pts = line(1e-155), spread = 1e-155),
                list(pts = box(1e-320), spread = 1e-320 * sqrt(1 + 1)),
                list(pts = box(1e-160), spread = 1e-160 * sqrt(1 + 1)))
  for (cs in cases) {
    lab <- format(cs$spread)
    r <- .tm_msgs(clip_target_for(cs$pts))
    expect_identical(r$msgs, paste0(
      "clip_target_for(): degenerate bbox: the points lie within ", format(cs$spread),
      " of each other, so the box has no extent to take a share of; using a disk of ",
      "radius 1 metre (1 unit of the CRS) around them. Pass `boundary` for a target of ",
      "real size, or an `expand` of 1 or more, which is taken as the radius in the units ",
      "of the CRS.\n"), info = lab)
    g <- sf::st_geometry(r$value)
    expect_identical(as.character(sf::st_geometry_type(g)), "POLYGON", info = lab)
    expect_equal(as.numeric(sf::st_area(g)), pi, tolerance = 1e-3, info = lab)
  }
  # A spread of 1e-150, whose share squares to a normal double, keeps the
  # buffer and the box it had: the same geometry st_buffer() and the
  # bounding box give, with the same message, or none for the box.
  p <- line(1e-150)
  r <- .tm_msgs(clip_target_for(p))
  expect_identical(r$msgs, "clip_target_for(): degenerate bbox; using small buffer around points.\n")
  want <- spatialkit:::.safe_make_valid(suppressWarnings(
    sf::st_buffer(sf::st_union(sf::st_geometry(p)), dist = 0.05 * 1e-150)))
  expect_identical(sf::st_geometry(r$value), want)
  expect_gt(as.numeric(sf::st_area(r$value)), 0)
  p <- box(1e-150)
  r <- .tm_msgs(clip_target_for(p))
  expect_identical(r$msgs, character())
  expect_identical(sf::st_geometry(r$value), spatialkit:::.safe_make_valid(
    sf::st_as_sfc(sf::st_bbox(p))))
})

test_that("clip_target_for() keeps its message and target for points on one line", {
  pts <- sf::st_as_sf(data.frame(x = 5e5 + c(0, 100, 200), y = rep(4e6, 3)),
                      coords = c("x", "y"), crs = 32617)
  r <- .tm_msgs(clip_target_for(pts))
  expect_identical(r$msgs, "clip_target_for(): degenerate bbox; using small buffer around points.\n")
  # Disks of 5 percent of the 200 m extent round the three points, unioned.
  expect_equal(as.numeric(sf::st_area(r$value)), 942.0472, tolerance = 1e-5)
})

test_that("?clip_target_for says what the target is for points at one location", {
  file <- testthat::test_path("..", "..", "man", "clip_target_for.Rd")
  rd <- if (file.exists(file)) tools::parse_Rd(file)
        else tools::Rd_db("spatialkit")[["clip_target_for.Rd"]]
  txt <- gsub("[[:space:]]+", " ", paste(unlist(rd), collapse = ""))
  expect_match(txt, paste0("When every point sits at one location, or the share of their spread ",
                           "the buffer would take (5 percent, or a fractional expand) is below ",
                           "that length, the box has no extent to take ",
                           "a share of, and the target is a disk of radius 1 unit of the CRS ",
                           "around it (1 m in a metre CRS, 1 km in a kilometre CRS), unless ",
                           "expand is 1 or more, which is then the radius; the message says so. ",
                           "Pass boundary for a target of real size."), fixed = TRUE)
})
