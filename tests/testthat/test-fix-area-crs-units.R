# tests/testthat/test-fix-area-crs-units.R
# ---------------------------------------------------------------------------
# sf::st_area() gives a projected layer's area in its CRS's unit squared when
# sf can name the unit, and in square metres when it cannot and the CRS
# carries a +to_meter: there it multiplies GEOS's area by to_meter^2.
# summarize_by_cell(area = TRUE) returned that as `cell_area`, which its help
# gives in the CRS's unit squared, and resolution_profile() divided the
# hull's area in square metres by the square of a range in CRS units.  Both
# now read the planar area in the CRS's unit squared (.planar_area()).
# summarize_by_cell(area = TRUE) also stopped with GEOS's bare message on
# cells that carry M throughout; they are now measured in x and y.
#
# The layer: 12 squares of 1000 m in EPSG:32615, 240 points (seed 2381).
# The scaled copies multiply every coordinate by `by` and state the unit as
# +to_meter = 1 / by, so they describe the same ground: by = 1e-6 (a unit of
# 1000 km), by = 0.304801 (a unit of about 3.28 m) and by = 2000 (a unit of
# 0.5 mm), units sf has no name for.  kriging_adequacy()'s sliver guard,
# make_folds()'s scale of supplied blocks and the strips create_grid_polygons()
# drops after its clip read areas beside lengths in the CRS's unit too.
# ---------------------------------------------------------------------------

.au_store <- new.env(parent = emptyenv())
.au_layer <- function() {
  if (is.null(.au_store$layer)) {
    crs <- sf::st_crs(32615)
    sq <- lapply(seq_len(12L), function(i) {
      x0 <- 5e5 + 1000 * ((i - 1L) %% 4L); y0 <- 4e6 + 1000 * ((i - 1L) %/% 4L)
      sf::st_polygon(list(cbind(x0 + c(0, 1000, 1000, 0, 0),
                                y0 + c(0, 0, 1000, 1000, 0))))
    })
    cells <- sf::st_sf(poly_id = sprintf("c%02d", 1:12),
                       geometry = sf::st_sfc(sq, crs = crs))
    d <- local({
      withr::local_seed(2381, .rng_kind = "Mersenne-Twister",
                        .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
      n <- 240
      x <- stats::runif(n, 0, 4000); y <- stats::runif(n, 0, 3000)
      data.frame(x = 5e5 + x, y = 4e6 + y,
                 z = sin(x / 300) + cos(y / 250) + stats::rnorm(n, 0, 0.3))
    })
    pts <- sf::st_as_sf(d, coords = c("x", "y"), crs = crs)
    .au_store$layer <- list(cells = cells, pts = pts)
  }
  .au_store$layer
}

# The layer with every coordinate multiplied by `by`, in a CRS whose unit is
# 1 / by metres.
.au_scale <- function(x, by) {
  p4 <- sf::st_crs(x)$proj4string
  p4 <- trimws(gsub("\\+units=[^ ]+", "", p4))
  crs <- sf::st_crs(paste(p4, sprintf("+to_meter=%.17g", 1 / by)))
  g <- sf::st_geometry(x) * by
  sf::st_geometry(x) <- sf::st_set_crs(g, crs)
  x
}

# An sfc of polygons with an M (or a Z and an M) value at every vertex.
.au_add_m <- function(g, dim = c("XYM", "XYZM")) {
  dim <- match.arg(dim)
  sf::st_sfc(lapply(g, function(p) {
    sf::st_polygon(lapply(p, function(r) {
      if (dim == "XYM") cbind(r, 3)
      else cbind(r, 0, 3)
    }), dim = dim)
  }), crs = sf::st_crs(g))
}

.au_cells_summary <- function(pts, cells) {
  a <- assign_features_to_polygons(pts, cells, "poly_id")
  summarize_by_cell(a, "z", cells_sf = cells, area = TRUE)
}

test_that("the scaled CRSs have a unit sf does not name, so sf::st_area() is in square metres", {
  L <- .au_layer()
  for (by in c(1e-6, 0.304801, 2000)) {
    s <- .au_scale(L$cells, by)
    expect_null(sf::st_crs(s)$ud_unit)
    expect_equal(sf::st_crs(s)$to_meter, 1 / by)
    expect_identical(as.character(units(sf::st_area(s))), "m^2")
  }
})

test_that("summarize_by_cell(area = TRUE) gives cell_area in the CRS's unit squared", {
  L <- .au_layer()
  base <- .au_cells_summary(L$pts, L$cells)
  expect_equal(base$cell_area, rep(1e6, 12))
  for (by in c(1e-6, 0.304801, 2000)) {
    s <- .au_cells_summary(.au_scale(L$pts, by), .au_scale(L$cells, by))
    # The shoelace area of a square of side 1000 * by.
    side <- 1000 * by
    expect_equal(s$cell_area, rep(side^2, 12), tolerance = 1e-9)
    expect_equal(s$cell_area, base$cell_area * by^2, tolerance = 1e-9)
    expect_identical(s$n, base$n)
    expect_equal(s$n_per_area, base$n_per_area / by^2, tolerance = 1e-9)
    expect_equal(s$n_per_area, as.numeric(s$n) / s$cell_area)
  }
})

test_that("resolution_profile() measures the hull and the range floor in the CRS's unit", {
  L <- .au_layer()
  run <- function(pts) {
    sac <- suppressWarnings(estimate_sac_range(pts, "z"))
    suppressWarnings(resolution_profile(pts, "z", sac = sac, seed = 7L))
  }
  base <- run(L$pts)
  hull <- sf::st_convex_hull(sf::st_union(sf::st_geometry(L$pts)))
  expect_identical(attr(base, "bounds")$area, as.numeric(sf::st_area(hull)))
  # A range of about 2 km over a hull of about 11 square km: a floor of 3.
  expect_gt(attr(base, "bounds")$floor, 2L)
  for (by in c(0.304801, 1e-6, 2000)) {
    pts <- .au_scale(L$pts, by)
    p <- run(pts)
    expect_equal(attr(p, "bounds")$area, attr(base, "bounds")$area * by^2, tolerance = 1e-9)
    expect_equal(attr(p, "bounds")$range, attr(base, "bounds")$range * by, tolerance = 1e-6)
    expect_true(is.finite(attr(p, "bounds")$floor))
    expect_identical(attr(p, "bounds")$floor, attr(base, "bounds")$floor)
  }
})

test_that("areas in metres and in a unit sf names are sf::st_area()'s as before", {
  L <- .au_layer()
  base <- .au_cells_summary(L$pts, L$cells)
  expect_identical(base$cell_area, as.numeric(sf::st_area(L$cells)))
  p <- suppressWarnings(resolution_profile(L$pts, "z", seed = 7L))
  hull <- sf::st_convex_hull(sf::st_union(sf::st_geometry(L$pts)))
  expect_identical(attr(p, "bounds")$area, as.numeric(sf::st_area(hull)))
  # US survey feet (EPSG:2264): sf names the unit and gives its square.
  ft_cells <- sf::st_transform(L$cells, 2264)
  ft_pts <- sf::st_transform(L$pts, 2264)
  s_ft <- suppressWarnings(.au_cells_summary(ft_pts, ft_cells))
  expect_identical(s_ft$cell_area, as.numeric(sf::st_area(ft_cells)))
  p_ft <- suppressWarnings(resolution_profile(ft_pts, "z", seed = 7L))
  hull_ft <- sf::st_convex_hull(sf::st_union(sf::st_geometry(ft_pts)))
  expect_identical(attr(p_ft, "bounds")$area, as.numeric(sf::st_area(hull_ft)))
})

test_that(".planar_area() changes sf::st_area() only where sf converted to square metres", {
  L <- .au_layer()
  expect_identical(.planar_area(L$cells), as.numeric(sf::st_area(L$cells)))
  expect_identical(.planar_area(sf::st_geometry(L$cells)),
                   as.numeric(sf::st_area(L$cells)))
  ft <- sf::st_transform(L$cells, 2264)
  expect_identical(.planar_area(ft), as.numeric(sf::st_area(ft)))
  # Lon/lat keeps its geodesic square metres; no CRS keeps GEOS's area.
  ll <- sf::st_transform(L$cells, 4326)
  expect_identical(.planar_area(ll), as.numeric(sf::st_area(ll)))
  na <- sf::st_set_crs(L$cells, NA)
  expect_identical(.planar_area(na), as.numeric(sf::st_area(na)))
  # A unit sf does not name: GEOS's area in that unit squared.
  for (by in c(1e-6, 0.304801, 2000)) {
    sc <- .au_scale(L$cells, by)
    expect_identical(.planar_area(sc), as.numeric(sf::st_area(sf::st_set_crs(sc, NA))))
    expect_equal(.planar_area(sc), as.numeric(sf::st_area(sc)) * by^2, tolerance = 1e-12)
  }
})

test_that("summarize_by_cell(area = TRUE) measures cells that carry M in x and y", {
  L <- .au_layer()
  base <- .au_cells_summary(L$pts, L$cells)
  for (dim in c("XYM", "XYZM")) {
    cm <- L$cells
    sf::st_geometry(cm) <- .au_add_m(sf::st_geometry(L$cells), dim)
    a <- assign_features_to_polygons(L$pts, L$cells, "poly_id")
    s <- summarize_by_cell(a, "z", cells_sf = cm, area = TRUE)
    expect_identical(s$cell_area, base$cell_area)
    expect_identical(s$n_per_area, base$n_per_area)
    expect_identical(attr(s, "area_error"), attr(base, "area_error"))
    # The cells come back with their own geometry, M included.
    expect_identical(class(sf::st_geometry(s)[[1]])[1], dim)
    expect_identical(sf::st_geometry(s), sf::st_geometry(cm))
  }
})

test_that("kriging_adequacy()'s sliver guard reads the box and the area in one unit", {
  skip_if_not_installed("gstat")
  L <- .au_layer()
  # A diagonal strip of 1 m by 1500 m: a box of about 1.1e6 over an area of
  # 1500, a ratio of 750 in every unit.
  strip <- sf::st_sfc(sf::st_polygon(list(cbind(
    5e5 + c(0, 1, 1061.66, 1060.66, 0), 4e6 + c(0, 0, 1060.66, 1060.66, 0)))),
    crs = sf::st_crs(32615))
  strip <- sf::st_sf(poly_id = "s", geometry = strip)
  r0 <- .cell_box_ratio(strip)
  expect_gt(r0, 500)
  expect_equal(.cell_box_ratio(L$cells), rep(1, 12))
  for (by in c(1e-6, 0.304801, 2000)) {
    expect_equal(.cell_box_ratio(.au_scale(strip, by)), r0, tolerance = 1e-9)
    expect_equal(.cell_box_ratio(.au_scale(L$cells, by)), rep(1, 12), tolerance = 1e-9)
  }
  run <- function(pts, cells) {
    a <- assign_features_to_polygons(pts, cells, "poly_id")
    sac <- suppressWarnings(estimate_sac_range(pts, "z"))
    kriging_adequacy(a, "z", cells, sac = sac, k = 3L, seed = 1L)
  }
  base <- run(L$pts, L$cells)
  kc <- grep("^kr_", names(base), value = TRUE)
  expect_false(anyNA(sf::st_drop_geometry(base)[kc]))
  # In a unit of 0.5 mm every square was left out as a sliver (a ratio of
  # 4e6); it is kriged as in metres.
  expect_no_warning(s <- run(.au_scale(L$pts, 2000), .au_scale(L$cells, 2000)))
  expect_equal(sf::st_drop_geometry(s)[kc], sf::st_drop_geometry(base)[kc],
               tolerance = 1e-6)
  expect_identical(attr(s, "cells_left_out"), attr(base, "cells_left_out"))
})

test_that("make_folds() compares supplied blocks with the range in the CRS's unit", {
  L <- .au_layer()
  blocks_of <- function(cells) {
    g <- sf::st_geometry(cells)
    sf::st_sf(geometry = sf::st_sfc(lapply(c(1, 3, 9, 11), function(i)
      sf::st_union(g[c(i, i + 1)])[[1]]), crs = sf::st_crs(cells)))
  }
  run <- function(pts, cells) {
    w <- character(0)
    f <- withCallingHandlers(
      make_folds(pts, k = 2, method = "block_kfold", blocks = blocks_of(cells),
                 auto_range = TRUE, response_var = "z", seed = 1L),
      warning = function(cnd) {
        w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
      })
    list(f = f, w = grep("median scale", w, value = TRUE))
  }
  # Blocks of 2000 m by 1000 m (a scale of 1414.2) under a range of about
  # 1986 m: the warning says so, and gives the size the range asks for.
  base <- run(L$pts, L$cells)
  expect_length(base$w, 1L)
  # The range comes from a variogram fit, whose last bits move with the
  # platform: the figures are read from the fit, not pinned.
  r0 <- as.numeric(base$f$params$sac_range)
  expect_gt(r0, 1500); expect_lt(r0, 2500)
  expect_match(base$w, sprintf(
    "median scale 1414.2) are smaller than the autocorrelation range (%.1f)", r0),
    fixed = TRUE)
  expect_match(base$w, sprintf("supply blocks at least %.0f units across.", ceiling(r0)),
               fixed = TRUE)
  for (by in c(1e-6, 0.304801, 2000)) {
    s <- run(.au_scale(L$pts, by), .au_scale(L$cells, by))
    expect_identical(s$f$assignment, base$f$assignment)
    # In a unit of 1000 km the blocks were read as 1414 units across and no
    # warning came.
    expect_length(s$w, 1L)
    # The message prints its figures with one decimal, so in a unit of
    # 1000 km they read 0.0; only the larger units are read back here.
    if (by < 0.01) next
    scale <- as.numeric(sub(".*median scale ([^)]+)\\).*", "\\1", s$w))
    range <- as.numeric(sub(".*autocorrelation range \\(([^)]+)\\)\\..*", "\\1", s$w))
    across <- as.numeric(sub(".*supply blocks at least ([^ ]+) units across\\.$", "\\1", s$w))
    expect_equal(scale, sqrt(2) * 1000 * by, tolerance = 5e-4)
    expect_equal(range, as.numeric(s$f$params$sac_range), tolerance = 5e-4)
    expect_gte(across, as.numeric(s$f$params$sac_range))
    expect_lte(across, as.numeric(s$f$params$sac_range) * 1.01)
  }
})

test_that("create_grid_polygons()'s clip drops an edge strip in every unit", {
  # One lattice cell of 1000 m and a part of it 0.5 margins wide along its
  # left edge: the strip a clip leaves where the boundary's edge lies on the
  # cell's edge but for rounding.  It was dropped in metres and kept in a
  # unit of 3.28 m, where sf::st_area() gave square metres against a bound in
  # the unit squared.
  for (by in c(1, 0.304801)) {
    crs <- if (by == 1) sf::st_crs(32615) else
      sf::st_crs(paste(trimws(gsub("\\+units=[^ ]+", "", sf::st_crs(32615)$proj4string)),
                       sprintf("+to_meter=%.17g", 1 / by)))
    x0 <- 5e5 * by; y0 <- 4e6 * by; side <- 1000 * by
    sq <- function(x, y, w, h) sf::st_polygon(list(cbind(x + c(0, w, w, 0, 0),
                                                         y + c(0, 0, h, h, 0))))
    lattice <- sf::st_sf(poly_id = 1L, geometry = sf::st_sfc(sq(x0, y0, side, side), crs = crs))
    bb <- sf::st_bbox(lattice)
    tol <- .edge_tol(bb)
    parts <- sf::st_sf(poly_id = c(1L, 1L),
                       geometry = sf::st_sfc(sq(x0, y0, tol / 2, side),
                                             sq(x0 + side / 2, y0, side / 2, side), crs = crs))
    kept <- .drop_edge_slivers(parts, lattice, tol)
    expect_identical(nrow(kept), 1L, info = paste("by", by))
  }
})
