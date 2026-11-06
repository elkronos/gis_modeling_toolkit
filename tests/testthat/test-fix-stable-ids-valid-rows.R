# tests/testthat/test-fix-stable-ids-valid-rows.R
# ---------------------------------------------------------------------------
# ensure_stable_poly_id(make_valid = TRUE) on a layer in degrees with s2 on.
#
# The rows returned were repaired with sf::st_make_valid(), which in
# longitude and latitude with s2 on is s2's rebuild of every row, valid or
# not: valid rows came back with their corners moved by up to 7e-8 degrees
# and a POLYGON as a MULTIPOLYGON, where with s2 off (GEOS's repair) or in a
# projected CRS no corner moved.  Only the rows s2 or GEOS calls invalid
# are rebuilt now (.make_valid_rows()), and the others are returned as they
# arrived.  The sort copy is repaired as before, so the IDs
# do not change.
# ---------------------------------------------------------------------------

.sv_s2 <- function(on, env = parent.frame()) {
  old <- sf::sf_use_s2()
  suppressMessages(sf::sf_use_s2(on))
  withr::defer(suppressMessages(sf::sf_use_s2(old)), envir = env)
}
.sv_q <- function(x) suppressWarnings(suppressMessages(x))

# North Carolina's first 8 counties in EPSG:4326, and two grids over them
# returned in degrees, every row valid under s2 and under GEOS.
.sv_layers <- function() {
  nc <- sf::st_transform(sf::st_read(system.file("shape/nc.shp", package = "sf"),
                                     quiet = TRUE)[1:8, ], 4326)
  old <- sf::sf_use_s2()
  suppressMessages(sf::sf_use_s2(TRUE))
  on.exit(suppressMessages(sf::sf_use_s2(old)))
  list(nc8     = nc,
       grid_cs = .sv_q(create_grid_polygons(nc, cellsize = 0.1, crs = 4326, quiet = TRUE)),
       grid_tc = .sv_q(create_grid_polygons(nc, target_cells = 40, crs = 4326, quiet = TRUE)))
}

# The largest distance between a corner of `a` and the nearest corner of `b`,
# both ways (Hausdorff distance over the corners), for two sfg.
.sv_corner_gap <- function(a, b) {
  x <- sf::st_coordinates(sf::st_sfc(a))[, 1:2, drop = FALSE]
  y <- sf::st_coordinates(sf::st_sfc(b))[, 1:2, drop = FALSE]
  d <- outer(x[, 1], y[, 1], "-")^2 + outer(x[, 2], y[, 2], "-")^2
  max(sqrt(c(apply(d, 1, min), apply(d, 2, min))))
}

# Each input row's ID, read through the column `orig` the input carries.
.sv_ids <- function(out) out$poly_id[order(out$orig)]

# Whether any two rows share interior, in a projected CRS by GEOS.
.sv_overlap <- function(x) {
  old <- sf::sf_use_s2()
  suppressMessages(sf::sf_use_s2(FALSE))
  on.exit(suppressMessages(sf::sf_use_s2(old)))
  p  <- sf::st_transform(sf::st_geometry(x), 32617)
  rl <- sf::st_relate(p, p, pattern = "2********")
  any(vapply(seq_along(rl), function(i) any(rl[[i]] != i), logical(1)))
}


test_that("valid rows in degrees come back as they arrived with s2 on, and with no corner moved with s2 off", {
  skip_if_not(file.exists(system.file("shape/nc.shp", package = "sf")))
  L <- .sv_layers()
  for (nm in names(L)) {
    x <- L[[nm]]
    x$orig <- seq_len(nrow(x))
    expect_true(local({ .sv_s2(TRUE);  all(sf::st_is_valid(x)) }), info = nm)
    expect_true(local({ .sv_s2(FALSE); all(sf::st_is_valid(x)) }), info = nm)
    on  <- local({ .sv_s2(TRUE);  .sv_q(ensure_stable_poly_id(x)) })
    off <- local({ .sv_s2(FALSE); .sv_q(ensure_stable_poly_id(x)) })
    from_on  <- on$orig
    from_off <- off$orig
    # With s2 on every row is the input's row, to the bit, and of its type.
    expect_identical(sf::st_geometry(on), sf::st_geometry(x)[from_on], info = nm)
    expect_identical(as.character(sf::st_geometry_type(on, by_geometry = TRUE)),
                     as.character(sf::st_geometry_type(x, by_geometry = TRUE))[from_on],
                     info = nm)
    # Under both settings no corner moves.
    gap_on  <- max(mapply(.sv_corner_gap, sf::st_geometry(on),
                          sf::st_geometry(x)[from_on]))
    gap_off <- max(mapply(.sv_corner_gap, sf::st_geometry(off),
                          sf::st_geometry(x)[from_off]))
    expect_identical(gap_on, 0, info = nm)
    expect_identical(gap_off, 0, info = nm)
    # The same IDs under both settings, and no row overlaps another.
    expect_identical(.sv_ids(on), .sv_ids(off), info = nm)
    expect_false(.sv_overlap(on), info = nm)
    expect_false(.sv_overlap(off), info = nm)
  }
})

test_that("the IDs are those the code before the change gave", {
  skip_if_not(file.exists(system.file("shape/nc.shp", package = "sf")))
  L <- lapply(.sv_layers(), function(x) { x$orig <- seq_len(nrow(x)); x })
  .sv_s2(TRUE)
  # Pinned: the IDs of the 8 counties before the change, under either setting.
  expect_identical(.sv_ids(.sv_q(ensure_stable_poly_id(L$nc8))),
                   c(1L, 2L, 3L, 8L, 4L, 5L, 7L, 6L))
  # The code before the change rebuilt every row with s2 first, and then the
  # sort copy of the rebuilt rows.  Rows rebuilt beforehand are valid under
  # s2, so they now pass through as they are and their sort copy is the one
  # the code before the change measured: its IDs are those it gave.
  for (nm in names(L)) {
    x <- L[[nm]]
    before <- .sv_q(ensure_stable_poly_id(sf::st_make_valid(x)))
    expect_true(all(sf::st_is_valid(sf::st_make_valid(x))), info = nm)
    expect_identical(.sv_ids(.sv_q(ensure_stable_poly_id(x))), .sv_ids(before), info = nm)
  }
})

test_that("a row s2 refuses is still rebuilt with s2 on, and the rows beside it are not", {
  # The 0.5 degree square with one corner listed twice, a bowtie, a valid square.
  sq  <- sf::st_polygon(list(rbind(c(-94, 39), c(-93.5, 39), c(-93.5, 39),
                                   c(-93.5, 39.5), c(-94, 39.5), c(-94, 39))))
  bow <- sf::st_polygon(list(rbind(c(-92, 39), c(-91.5, 39.5), c(-91.5, 39),
                                   c(-92, 39.5), c(-92, 39))))
  ok  <- sf::st_polygon(list(rbind(c(-90.0001234, 39.0004321), c(-89.4998765, 39.0004321),
                                   c(-89.4998765, 39.4995678), c(-90.0001234, 39.4995678),
                                   c(-90.0001234, 39.0004321))))
  x <- sf::st_sf(v = 1:3, geometry = sf::st_sfc(sq, bow, ok, crs = 4326))
  .sv_s2(TRUE)
  expect_identical(sf::st_is_valid(x), c(FALSE, FALSE, TRUE))
  out <- .sv_q(ensure_stable_poly_id(x))
  g <- sf::st_geometry(out)[order(out$v)]
  rebuilt <- suppressWarnings(sf::st_make_valid(sf::st_geometry(x)))
  # The square with a corner listed twice gets s2's rebuild, as before, and
  # comes back valid under s2.  s2's rebuild leaves the bowtie's crossing
  # edges, so the bowtie gets GEOS's repair of its coordinates, as with s2
  # off (BUG-195; this line pinned s2's rebuild before).
  expect_identical(g[[1]], rebuilt[[1]])
  geos <- suppressWarnings(sf::st_make_valid(sf::st_set_crs(sf::st_geometry(x), NA)))
  expect_identical(g[[2]], geos[[2]])
  expect_true(sf::st_is_valid(g[1]))
  # The valid row comes back as it arrived.
  expect_identical(g[[3]], ok)
  # The repair itself, on an sfc and on an sf, s2 on: only the refused rows change.
  h <- spatialkit:::.make_valid_rows(sf::st_geometry(x))
  expect_identical(h[[3]], ok)
  expect_identical(h[[1]], rebuilt[[1]])
  expect_identical(sf::st_geometry(spatialkit:::.make_valid_rows(x)), h)
  # A layer of valid rows is returned unchanged.
  expect_identical(spatialkit:::.make_valid_rows(x[3, ]), x[3, ])
})

test_that("a row s2 calls valid and GEOS does not is still rebuilt with s2 on", {
  # s2 reads rings without their orientation, so it calls valid a hole outside
  # its shell, a part inside another part and a ring of collinear points,
  # which GEOS calls invalid.  These rows are rebuilt as any invalid row is,
  # and come back valid under both engines; the valid row beside them comes
  # back as it arrived.
  P  <- function(...) sf::st_polygon(list(...))
  sq <- function(x0, y0, w) rbind(c(x0, y0), c(x0 + w, y0), c(x0 + w, y0 + w),
                                  c(x0, y0 + w), c(x0, y0))
  ok <- P(sq(-80.0001234, 30.0004321, 0.3))
  cases <- list(
    hole_outside  = P(sq(-90, 39, 0.5), sq(-89, 39, 0.1)),
    nested_shells = sf::st_multipolygon(list(list(sq(-90, 39, 0.5)),
                                             list(sq(-89.8, 39.2, 0.1)))),
    line_ring     = P(rbind(c(-90, 39), c(-89, 39), c(-88, 39), c(-90, 39))))
  geos_valid <- function(g) sf::st_is_valid(sf::st_set_crs(sf::st_geometry(g), NA))
  for (nm in names(cases)) {
    x <- sf::st_sf(v = 1:2, geometry = sf::st_sfc(cases[[nm]], ok, crs = 4326))
    on <- local({
      .sv_s2(TRUE)
      expect_identical(sf::st_is_valid(x), c(TRUE, TRUE), info = nm)
      .sv_q(ensure_stable_poly_id(x))
    })
    expect_identical(geos_valid(x), c(FALSE, TRUE), info = nm)
    g <- sf::st_geometry(on)[order(on$v)]
    expect_identical(geos_valid(g), c(TRUE, TRUE), info = nm)
    expect_true(local({ .sv_s2(TRUE); all(sf::st_is_valid(g)) }), info = nm)
    expect_identical(g[[2]], ok, info = nm)
    # The rebuilt row is s2's rebuild of the row, as the code before gave.
    expect_identical(g[[1]], local({ .sv_s2(TRUE); suppressWarnings(
      sf::st_make_valid(sf::st_geometry(x)))[[1]] }), info = nm)
  }
  # The hole outside its shell, projected: the same area under both settings,
  # and GEOS can cut it.
  x <- sf::st_sf(v = 1:2, geometry = sf::st_sfc(cases$hole_outside, ok, crs = 4326))
  area <- function(s) local({
    .sv_s2(s)
    o <- .sv_q(ensure_stable_poly_id(x))
    p <- sf::st_transform(sf::st_geometry(o)[order(o$v)][1], 32615)
    cut <- suppressWarnings(sf::st_intersection(p, sf::st_as_sfc(sf::st_bbox(p))))
    c(as.numeric(sf::st_area(p)), as.numeric(sf::st_area(cut)))
  })
  a_on <- area(TRUE)
  a_off <- area(FALSE)
  expect_equal(a_on, a_off, tolerance = 1e-6)
  expect_equal(a_on[1], a_on[2], tolerance = 1e-9)
})

test_that("a projected layer and a layer in degrees with s2 off keep GEOS's repair of every row", {
  sq <- sf::st_polygon(list(rbind(c(-94, 39), c(-93.5, 39), c(-93.5, 39.5),
                                  c(-94, 39.5), c(-94, 39))))
  bow <- sf::st_polygon(list(rbind(c(-92, 39), c(-91.5, 39.5), c(-91.5, 39),
                                   c(-92, 39.5), c(-92, 39))))
  ll <- sf::st_sf(v = 1:2, geometry = sf::st_sfc(sq, bow, crs = 4326))
  pr <- sf::st_transform(ll, 32615)
  local({
    .sv_s2(TRUE)
    expect_identical(spatialkit:::.make_valid_rows(pr),
                     suppressWarnings(sf::st_make_valid(pr)))
  })
  local({
    .sv_s2(FALSE)
    expect_identical(spatialkit:::.make_valid_rows(ll),
                     suppressWarnings(sf::st_make_valid(ll)))
  })
})

test_that("create_voronoi_polygons() on seeds in degrees gives every seed its own cell under both settings", {
  # A square with a notch cut from the top, in degrees, and five seeds in a
  # row across the top: the middle one stands in the notch, outside the
  # boundary, and its own cell is the strip below the notch, not the cell
  # nearest it.
  k <- function(x, y) cbind(x / 7.3e3 - 79.3141592653, y / 7.3e3 + 35.2718281828)
  notch <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(k(
    c(0, 1000, 1000, 550, 550, 450, 450, 0, 0),
    c(0, 0, 1000, 1000, 300, 300, 1000, 1000, 0)))), crs = 4326))
  xy <- k(c(100, 440, 500, 560, 900), c(950, 950, 950, 960, 940))
  seeds <- sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2]), coords = c("x", "y"), crs = 4326)
  seeds$seed_id <- seq_len(nrow(seeds))
  seeds$method  <- "provided"
  for (s in c(TRUE, FALSE)) {
    tv <- local({ .sv_s2(s); .sv_q(create_voronoi_polygons(seeds, boundary = notch, quiet = TRUE)) })
    pts <- sf::st_transform(seeds, sf::st_crs(tv$cells))
    expect_false(lengths(sf::st_intersects(pts[3, ], sf::st_union(sf::st_transform(notch, sf::st_crs(tv$cells))))) > 0L,
                 info = paste("s2", s))
    own <- sf::st_nearest_feature(suppressWarnings(sf::st_point_on_surface(sf::st_geometry(tv$cells))), pts)
    expect_identical(tv$index[own], tv$cells$cell_id, info = paste("s2", s))
    expect_false(anyNA(tv$index), info = paste("s2", s))
  }
})

test_that("create_grid_polygons_cached() in degrees with s2 on returns the builder's cells, corners and types", {
  skip_if_not(file.exists(system.file("shape/nc.shp", package = "sf")))
  nc <- .sv_layers()$nc8
  .sv_s2(TRUE)
  g <- .sv_q(create_grid_polygons_cached(nc, cellsize = 0.1, crs = 4326, quiet = TRUE,
                                         cache_env = new.env()))
  u <- .sv_q(create_grid_polygons(nc, cellsize = 0.1, crs = 4326, quiet = TRUE))
  expect_identical(nrow(g), nrow(u))
  ctr <- function(z) sf::st_centroid(sf::st_transform(sf::st_geometry(z), 32617))
  m <- sf::st_nearest_feature(ctr(g), ctr(u))
  expect_identical(sort(m), seq_len(nrow(u)))
  gap <- max(mapply(.sv_corner_gap, sf::st_geometry(g), sf::st_geometry(u)[m]))
  expect_lt(gap, 1e-9)
  expect_identical(as.character(sf::st_geometry_type(g, by_geometry = TRUE)),
                   as.character(sf::st_geometry_type(u, by_geometry = TRUE))[m])
})

test_that("the help says only the rows s2 or GEOS calls invalid are rebuilt", {
  page <- function(f) {
    file <- test_path("..", "..", "man", f)
    rd <- if (file.exists(file)) tools::parse_Rd(file)
          else tools::Rd_db("spatialkit")[[f]]
    expect_s3_class(rd, "Rd")
    gsub("[[:space:]]+", " ", paste(unlist(rd), collapse = ""))
  }
  e <- page("ensure_stable_poly_id.Rd")
  expect_match(e, paste0("only the rows s2 or GEOS does not call valid are ",
                         "rebuilt, and the other rows are returned as they arrived"),
               fixed = TRUE)
  expect_no_match(e, "only the rows sf::st_is_valid() does not call valid", fixed = TRUE)
  cgp <- page("create_grid_polygons_cached.Rd")
  expect_match(cgp, "rebuilds the cells s2 or GEOS calls invalid with s2's rebuild", fixed = TRUE)
  expect_no_match(cgp, "repairs cells in degrees with s2 or with GEOS", fixed = TRUE)
})
