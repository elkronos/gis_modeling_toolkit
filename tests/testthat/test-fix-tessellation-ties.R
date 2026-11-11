# The tessellation under the package's tie rule (BUG-106, BUG-108, BUG-109,
# BUG-112, BUG-201): a point on a cell edge, or snapped as near two cells,
# takes the smaller cell_id; co-circular points are triangulated in a frame
# that does not move with the unit or the origin; near-coincident points are
# merged by that frame; a side within rounding of a whole number of cells
# gets that number.  So the cells and the index are the same in another unit
# or with another origin.

.tt_crs_m    <- function() sf::st_crs(26915)
.tt_crs_usft <- function() sf::st_crs("+proj=utm +zone=15 +datum=NAD83 +units=us-ft +no_defs")
.tt_tm <- function(x0, y0, units = "m") sf::st_crs(sprintf(paste0(
  "+proj=tmerc +lat_0=0 +lon_0=-93 +k=0.9996 +x_0=%s +y_0=%s +datum=NAD83 ",
  "+units=%s +no_defs"), x0, y0, units))
.tt_crs <- function() list(m = .tt_crs_m(), usft = .tt_crs_usft(), m01 = .tt_tm(500000.1, 0.1))
# Metres in each frame's unit.
.tt_unit <- c(m = 1, usft = 3937 / 1200, m01 = 1)
.tt_frames <- function(x) lapply(.tt_crs(), function(cr) sf::st_transform(x, cr))

# A lattice of nx by ny points at 250 m from (450000, 4550000).
.tt_lattice <- function(nx, ny, jitter = 0) {
  g <- expand.grid(i = 0:(nx - 1), j = 0:(ny - 1))
  x <- 4.5e5 + g$i * 250; y <- 4.55e6 + g$j * 250
  if (jitter > 0) {
    x <- x + withr::with_seed(99, stats::runif(length(x), -jitter, jitter))
    y <- y + withr::with_seed(98, stats::runif(length(y), -jitter, jitter))
  }
  sf::st_as_sf(data.frame(x = x, y = y, z = seq_along(x) + 0), coords = c("x", "y"),
               crs = .tt_crs_m())
}

.tt_quiet <- function(expr) suppressWarnings(suppressMessages(
  logger::with_log_threshold(expr, threshold = logger::FATAL, namespace = "spatialkit",
                             index = 2)))

.tt_rect <- function(x0, y0, w, h) sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
  c(x0, y0), c(x0 + w, y0), c(x0 + w, y0 + h), c(x0, y0 + h), c(x0, y0)))), crs = .tt_crs_m()))

# Each cell keyed by its centroid in metres, to the centimetre.
.tt_key <- function(cells) {
  xy <- suppressWarnings(sf::st_coordinates(sf::st_centroid(
    sf::st_transform(sf::st_geometry(cells), .tt_crs_m()))))
  paste(round(xy[, 1L], 2), round(xy[, 2L], 2))
}
# The cell of each row, by that key.
.tt_row_cells <- function(tz) .tt_key(tz$cells)[match(tz$index, tz$cells$cell_id)]

test_that("build_tessellation() indexes a lattice point on a cell edge to the smaller cell_id in every frame", {
  lat <- .tt_lattice(40, 30)
  box <- sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(lat)))
  F <- .tt_frames(lat); B <- .tt_frames(box)
  for (meth in c("square", "hex")) {
    tz <- lapply(names(F), function(f) .tt_quiet(build_tessellation(
      F[[f]], boundary = B[[f]], method = meth, cellsize = 1000 * .tt_unit[[f]], quiet = TRUE)))
    names(tz) <- names(F)
    # The rule, read in metres: of the cells within a micrometre of the
    # point, the smallest cell_id.
    near <- sf::st_is_within_distance(lat, tz$m$cells, dist = 1e-6)
    want <- vapply(near, function(r) min(tz$m$cells$cell_id[r]), integer(1))
    expect_identical(tz$m$index, want, info = meth)
    expect_true(any(lengths(near) > 1L), info = meth)
    # Before, 80 (US feet) and 226 (origin moved) rows of the square grid,
    # and 0 and 156 of the hex grid, were in another cell than in metres.
    for (f in names(F)) {
      expect_setequal(.tt_key(tz[[f]]$cells), .tt_key(tz$m$cells))
      expect_identical(.tt_row_cells(tz[[f]]), .tt_row_cells(tz$m), info = paste(meth, f))
      expect_identical(tz[[f]]$params$snapped$n, 0L, info = paste(meth, f))
    }
  }
})

test_that("build_tessellation() indexes a jittered lattice by the cell that holds each point, as before", {
  lat <- .tt_lattice(40, 30, jitter = 60)
  box <- sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(.tt_lattice(40, 30))))
  F <- .tt_frames(lat); B <- .tt_frames(box)
  for (meth in c("square", "hex")) {
    tz <- lapply(names(F), function(f) .tt_quiet(build_tessellation(
      F[[f]], boundary = B[[f]], method = meth, cellsize = 1000 * .tt_unit[[f]], quiet = TRUE)))
    names(tz) <- names(F)
    hit <- sf::st_intersects(lat, tz$m$cells)
    want <- vapply(hit, function(r) if (length(r)) tz$m$cells$cell_id[r[1L]] else NA_integer_,
                   integer(1))
    # A point outside every cell by at most a thousandth of the cell width
    # (1 m) is snapped to the nearest cell; one further out is NA.
    out <- which(is.na(want))
    near <- sf::st_nearest_feature(lat[out, ], tz$m$cells)
    gap <- as.numeric(sf::st_distance(lat[out, ], tz$m$cells[near, ], by_element = TRUE))
    want[out] <- ifelse(gap <= 1, tz$m$cells$cell_id[near], NA_integer_)
    expect_identical(tz$m$index, want, info = meth)
    expect_identical(tz$m$params$snapped$which, out[gap <= 1], info = meth)
    expect_identical(sum(is.na(tz$m$index)), 67L, info = meth)
    for (f in names(F))
      expect_identical(.tt_row_cells(tz[[f]]), .tt_row_cells(tz$m), info = paste(meth, f))
  }
})

test_that("build_tessellation() snaps a point as near two cells to the smaller cell_id in every frame", {
  # A jittered lattice inside a 9100 x 6100 m box, and 180 points 0.1 to 0.85 m
  # outside it, each beside a corner two boundary cells share.
  box <- .tt_rect(4.5e5, 4.55e6, 9100, 6100)
  xs <- 4.5e5 + 1000 * (1:9); ys <- 4.55e6 + 1000 * (1:6)
  d <- c(0.1, 0.25, 0.4, 0.55, 0.7, 0.85)
  out <- rbind(data.frame(x = rep(xs, each = 6), y = 4.55e6 - d),
               data.frame(x = rep(xs, each = 6), y = 4.55e6 + 6100 + d),
               data.frame(x = 4.5e5 - d, y = rep(ys, each = 6)),
               data.frame(x = 4.5e5 + 9100 + d, y = rep(ys, each = 6)))
  out$z <- 0
  op <- sf::st_as_sf(out, coords = c("x", "y"), crs = .tt_crs_m())
  # One point 0.1 m outside a single cell, and one 5 m outside (beyond the
  # thousandth of the cell width).
  extra <- sf::st_as_sf(data.frame(x = c(4.5e5 + 500, 4.5e5 + 500), y = c(4.55e6 - 0.1, 4.55e6 - 5),
                                   z = 0), coords = c("x", "y"), crs = .tt_crs_m())
  inner <- .tt_lattice(37, 25, jitter = 60)
  inner <- inner[lengths(sf::st_intersects(inner, box)) > 0L, ]
  pts <- rbind(inner, op, extra)
  io <- nrow(inner) + seq_len(nrow(op))
  F <- .tt_frames(pts); B <- .tt_frames(box)
  tz <- lapply(names(F), function(f) .tt_quiet(build_tessellation(
    F[[f]], boundary = B[[f]], method = "square", cellsize = 1000 * .tt_unit[[f]], quiet = TRUE)))
  names(tz) <- names(F)
  m <- tz$m
  # Each snapped point in the smaller cell_id of the two cells at its corner.
  dd <- sf::st_distance(op, m$cells)
  two <- t(apply(matrix(as.numeric(dd), nrow(op)), 1L, function(r) order(r)[1:2]))
  expect_identical(m$index[io], as.integer(pmin(m$cells$cell_id[two[, 1L]],
                                                m$cells$cell_id[two[, 2L]])))
  # The single cell's point snaps to it; the far point stays NA.
  n <- nrow(pts)
  expect_identical(m$index[n - 1L],
                   m$cells$cell_id[lengths(sf::st_intersects(m$cells, sf::st_buffer(extra[1, ], 0.2))) > 0L])
  expect_true(is.na(m$index[n]))
  expect_identical(m$params$snapped$n, 181L)
  expect_equal(m$params$snapped$distance[seq_along(io)], rep(d, 30), tolerance = 1e-6)
  # Before, 4 (US feet) and 8 (origin moved) of the 180 went to another cell
  # than in metres.
  for (f in names(F)) {
    expect_setequal(.tt_key(tz[[f]]$cells), .tt_key(m$cells))
    expect_identical(.tt_row_cells(tz[[f]]), .tt_row_cells(m), info = f)
    expect_identical(tz[[f]]$params$snapped$n, 181L, info = f)
    expect_identical(tz[[f]]$params$snapped$which, m$params$snapped$which, info = f)
    expect_equal(tz[[f]]$params$snapped$distance / .tt_unit[[f]], m$params$snapped$distance,
                 tolerance = 1e-6, info = f)
  }
})

test_that("a snapped point within the relative part of the tie goes to the smaller cell_id", {
  # Two cells of a CRS-less layer near the origin, where the relative part
  # of the tie (one part in ten million of the distance) is as large as the
  # margin of an edge.  A point 9.9e-4 below their shared corner, a hair to
  # the right, is outside both and nearer cell 2 by `delta`.
  sq <- function(x0, y0, x1, y1)
    sf::st_polygon(list(rbind(c(x0, y0), c(x1, y0), c(x1, y1), c(x0, y1), c(x0, y0))))
  cells <- sf::st_sf(cell_id = 1:2, geometry = sf::st_sfc(sq(-1, 0, 0, 0.5), sq(0, 0, 1, 0.5)))
  tol <- .edge_tol(sf::st_bbox(cells))
  d <- 9.9e-4
  for (delta in c(0.5, 1.5, 2.5) * 1e-10) {
    p <- sf::st_sf(z = 1, geometry = sf::st_sfc(sf::st_point(c(sqrt(2 * d * delta), -d))))
    dd <- as.numeric(sf::st_distance(p, cells))
    ix <- .tt_quiet(.build_point_cell_index(p, cells))
    expect_identical(attr(ix, "snapped")$n, 1L)
    expect_equal(attr(ix, "snapped")$distance, d, tolerance = 1e-12)
    if (delta == 1.5e-10) {
      # Beyond the margin alone, within the margin plus the relative part:
      # a tie only through the relative part.
      expect_gt(dd[1] - dd[2], tol)
      expect_lte(dd[1] - dd[2], tol + dd[2] * .TIE_TOL)
    }
    # Within the tie (0.5e-10, 1.5e-10): the smaller cell_id.  Beyond it
    # (2.5e-10): the nearer cell.
    expect_identical(as.integer(ix), if (delta < 2e-10) 1L else 2L, info = format(delta))
  }
})

# Each triangle as the sorted rows at its corners.
.tt_tri_keys <- function(cells, l) {
  xy <- sf::st_coordinates(l)[, 1:2]; k <- paste(xy[, 1], xy[, 2])
  vapply(sf::st_geometry(cells), function(g) {
    r <- sf::st_coordinates(g)[1:3, 1:2]
    paste(sort(match(paste(r[, 1], r[, 2]), k)), collapse = "-")
  }, "")
}

test_that("build_tessellation(method = \"triangles\") splits a lattice the same way in every frame", {
  skip_if_not_installed("geometry")
  for (dims in list(c(10L, 10L), c(37L, 25L))) {
    lat <- .tt_lattice(dims[1], dims[2])
    F <- .tt_frames(lat)
    tz <- lapply(F, function(l) .tt_quiet(build_tessellation(l, method = "triangles", quiet = TRUE)))
    k <- mapply(function(t, l) .tt_tri_keys(t$cells, l), tz, F, SIMPLIFY = FALSE)
    n_tri <- 2L * (dims[1] - 1L) * (dims[2] - 1L)
    # Before, 182 and 1827 triangles with the origin moved, 20 and 99 of them
    # flat, and 66 and 860 triangles in US feet not in the metre set.
    for (f in names(F)) {
      a <- as.numeric(sf::st_area(tz[[f]]$cells)) / .tt_unit[[f]]^2
      expect_length(a, n_tri)
      expect_gt(min(a), 0.99 * 250^2 / 2)
      expect_identical(k[[f]], k$m, info = paste(dims[1], f))
      expect_identical(k[[f]][tz[[f]]$index], k$m[tz$m$index], info = paste(dims[1], f))
    }
  }
})

test_that("build_tessellation(method = \"triangles\") on a jittered lattice gives what it gave", {
  skip_if_not_installed("geometry")
  for (dims in list(c(10L, 10L), c(37L, 25L))) {
    lat <- .tt_lattice(dims[1], dims[2], jitter = 60)
    F <- .tt_frames(lat)
    tz <- lapply(F, function(l) .tt_quiet(build_tessellation(l, method = "triangles", quiet = TRUE)))
    k <- mapply(function(t, l) .tt_tri_keys(t$cells, l), tz, F, SIMPLIFY = FALSE)
    for (f in names(F)) {
      expect_length(k[[f]], if (dims[1] == 10L) 186L else 1828L)
      expect_setequal(k[[f]], k$m)
    }
    # The triangles qhull draws on the coordinates centred as before.
    xy <- sf::st_coordinates(lat)
    ctr <- (apply(xy, 2L, min) + apply(xy, 2L, max)) / 2
    raw <- geometry::delaunayn(sweep(xy, 2L, ctr))
    expect_setequal(k$m, apply(raw, 1L, function(r) paste(sort(r), collapse = "-")))
  }
})

test_that("the GEOS fallback of the triangulation splits a lattice the same way in every frame", {
  lat <- .tt_lattice(10, 10)
  F <- .tt_frames(lat)
  tri <- lapply(F, function(l) .triangulate_in_frame(sf::st_coordinates(l)[, 1:2], sf::st_crs(l)))
  k <- mapply(.tt_tri_keys, tri, F, SIMPLIFY = FALSE)
  for (f in names(F)) {
    a <- as.numeric(sf::st_area(tri[[f]])) / .tt_unit[[f]]^2
    expect_length(a, 162L)
    expect_gt(min(a), 0.99 * 250^2 / 2)
    expect_identical(k[[f]], k$m, info = f)
    # Each triangle is drawn on the points' own coordinates.
    expect_false(anyNA(unlist(strsplit(k[[f]], "-"))), info = f)
  }
})

test_that(".dedup_points() merges copies a hair apart in every unit and origin", {
  x <- withr::with_seed(4, 4.5e5 + stats::runif(60, 0, 5000))
  y <- withr::with_seed(5, 4.55e6 + stats::runif(60, 0, 5000))
  crs_km <- sf::st_crs("+proj=utm +zone=15 +datum=NAD83 +units=km +no_defs")
  for (eps in c(0, 3e-8, 1e-9, 1e-6, 1)) {
    pm <- sf::st_as_sf(data.frame(x = c(x, x[1:20] + eps), y = c(y, y[1:20]), z = 0),
                       coords = c("x", "y"), crs = .tt_crs_m())
    F <- c(.tt_frames(pm), list(km = sf::st_transform(pm, crs_km)))
    kept <- lapply(F, function(l) {
      l$row <- seq_len(nrow(l))
      .dedup_points(l)$row
    })
    vor <- vapply(F, function(l) nrow(.tt_quiet(build_tessellation(l, method = "voronoi", quiet = TRUE))$cells),
                  integer(1))
    # Before, copies 3e-8 m apart were 80 rows in metres and 65 in
    # kilometres; 1e-9 m apart 80, 62, 79 and 60.
    for (f in names(F)) expect_identical(kept[[f]], kept$m, info = paste(eps, f))
    expect_true(all(vor == vor[["m"]]), info = format(eps))
    if (eps %in% c(0, 3e-8, 1e-9)) expect_identical(kept$m, 1:60, info = format(eps))
    if (eps == 1) expect_identical(kept$m, 1:80)
  }
})

test_that(".dedup_points() keys MULTIPOINT features and rows without a location as before", {
  mp <- sf::st_multipoint(rbind(c(0, 0), c(1, 1), c(2, 2)))
  g <- sf::st_sfc(mp, mp, sf::st_point(c(5, 5)), sf::st_point(c(5, 5 + 1e-12)),
                  sf::st_point(), sf::st_point(), sf::st_point(c(Inf, 1)),
                  sf::st_point(c(1e200, 1)), sf::st_point(c(1e200, 1)), crs = 32632)
  out <- .dedup_points(sf::st_sf(id = seq_along(g), geometry = g))
  expect_identical(out$id, c(1L, 3L, 5L, 7L, 8L))
})

test_that("create_grid_polygons(cellsize =) gives a side that is a whole number of cells that number in every frame", {
  box <- .tt_rect(4.5e5, 4.55e6, 9000, 6000)
  B <- .tt_frames(box)
  # Before: square 54 / 54 / 60 cells, the smallest 5e-8 m^2 with the origin
  # moved; hex 76 / 76 / 79; square without the clip 54 / 54 / 60.
  for (type in c("square", "hex")) for (clip in c(TRUE, FALSE)) {
    gr <- lapply(names(B), function(f) .tt_quiet(create_grid_polygons(
      B[[f]], type = type, cellsize = 1000 * .tt_unit[[f]], clip = clip, quiet = TRUE)))
    names(gr) <- names(B)
    want <- if (type == "square") 54L else if (clip) 76L else 84L
    for (f in names(B)) {
      a <- as.numeric(sf::st_area(gr[[f]])) / .tt_unit[[f]]^2
      expect_identical(nrow(gr[[f]]), want, info = paste(type, clip, f))
      expect_gt(min(a), 2e5)
      expect_setequal(.tt_key(gr[[f]]), .tt_key(gr$m))
    }
  }
  # The square side derived from target_cells: 1000 m here.  Without the
  # clip no sliver drop hides an extra row or column: with the origin moved
  # the width is 9.000000000000087 sides, which sf counts as 10.
  for (clip in c(TRUE, FALSE)) {
    gr <- lapply(B, function(b) .tt_quiet(create_grid_polygons(
      b, type = "square", target_cells = 54, clip = clip, quiet = TRUE)))
    for (f in names(B)) {
      expect_identical(nrow(gr[[f]]), 54L, info = paste(clip, f))
      expect_setequal(.tt_key(gr[[f]]), .tt_key(gr$m))
    }
  }
  # build_tessellation() on a 37 x 25 lattice with its own box (54 / 54 / 70
  # cells before).
  lat <- .tt_lattice(37, 25)
  F <- .tt_frames(lat)
  tz <- lapply(names(F), function(f) .tt_quiet(build_tessellation(
    F[[f]], boundary = sf::st_as_sf(sf::st_as_sfc(sf::st_bbox(F[[f]]))), method = "square",
    cellsize = 1000 * .tt_unit[[f]], quiet = TRUE)))
  names(tz) <- names(F)
  for (f in names(F)) {
    expect_identical(nrow(tz[[f]]$cells), 54L, info = f)
    expect_identical(.tt_row_cells(tz[[f]]), .tt_row_cells(tz$m), info = f)
    expect_false(anyNA(tz[[f]]$index))
  }
})

test_that("create_grid_polygons(cellsize =) counts a side as whole only within the margin of an edge", {
  bb <- sf::st_bbox(.tt_rect(4.5e5, 4.55e6, 9000, 6000))
  m <- .edge_tol(bb)
  expect_equal(m, 4.556e-4)
  # A side over a whole number of cells by no more than the margin is that
  # number; one over by more gets another cell, however many cells it holds.
  expect_identical(.grid_side_counts(c(9000 + m / 2, 6000), c(1000, 1000), bb), c(9L, 6L))
  expect_identical(.grid_side_counts(c(9000 + 2 * m, 6000), c(1000, 1000), bb), c(10L, 6L))
  expect_identical(.grid_side_counts(c(1e6 + 0.05, 20), c(20, 20), bb), c(50001L, 1L))
  expect_identical(.grid_side_counts(c(8999.9995, 6000), c(1000, 1000), bb), c(9L, 6L))
  expect_identical(.grid_side_counts(c(0.5, 0.5), c(1000, 1000), bb), c(1L, 1L))
  # Boxes 0.2 mm (within the margin of 0.455 mm) and 0.8 mm (beyond it) over
  # 9000 m.  Points on the east edge, inside the boundary, are in a cell and
  # are not snapped.  Before the margin was absolute, the 0.8 mm box had 54
  # cells and its 59 edge points were snapped.
  for (extra in c(2e-4, 8e-4)) {
    b <- .tt_rect(4.5e5, 4.55e6, 9000 + extra, 6000)
    pts <- sf::st_as_sf(data.frame(x = 4.5e5 + 9000 + extra, y = 4.55e6 + 100 + (0:58) * 100,
                                   z = 1), coords = c("x", "y"), crs = .tt_crs_m())
    for (clip in c(TRUE, FALSE)) {
      tz <- .tt_quiet(build_tessellation(pts, boundary = b, method = "square", cellsize = 1000,
                                         clip = clip, quiet = TRUE))
      expect_identical(nrow(tz$cells), if (extra < m) 54L else 60L, info = paste(extra, clip))
      expect_false(anyNA(tz$index), info = paste(extra, clip))
      expect_identical(tz$params$snapped$n, 0L, info = paste(extra, clip))
      if (extra > m) {
        # The narrow column holds them.
        x0 <- vapply(sf::st_geometry(tz$cells)[match(tz$index, tz$cells$cell_id)],
                     function(g) sf::st_bbox(g)[["xmin"]], 0)
        expect_equal(x0, rep(4.5e5 + 9000, length(x0)), tolerance = 1e-12, info = paste(extra, clip))
      }
    }
  }
  # A transect of (1e6 + 0.05) x 20 m at 20 m: 0.05 m over 50000 cells,
  # within one part in ten million of that number but a hundred times the
  # margin.  Its east edge is in the 50001st cell.  Under a relative count
  # the grid had 50000 cells, and the three points below were NA.
  b <- .tt_rect(1e5, 4.55e6, 1e6 + 0.05, 20)
  pts <- sf::st_as_sf(data.frame(x = 1e5 + 1e6 + 0.05 - c(0, 0.001, 0.003), y = 4.55e6 + 10,
                                 z = 1), coords = c("x", "y"), crs = .tt_crs_m())
  tz <- .tt_quiet(build_tessellation(pts, boundary = b, method = "square", cellsize = 20,
                                     clip = FALSE, quiet = TRUE))
  expect_identical(nrow(tz$cells), 50001L)
  expect_false(anyNA(tz$index))
  expect_identical(tz$params$snapped$n, 0L)
})

test_that("create_grid_polygons() drops the slivers of cells on an inner edge of the boundary, and keeps narrow parts", {
  # An L-shaped boundary whose inner corner lies on grid lines.
  L <- sf::st_sf(geometry = sf::st_union(rbind(.tt_rect(4.5e5, 4.55e6, 9000, 3000),
                                              .tt_rect(4.5e5, 4.55e6, 4000, 6000))))
  B <- .tt_frames(L)
  for (type in c("square", "hex")) {
    gr <- lapply(names(B), function(f) .tt_quiet(create_grid_polygons(
      B[[f]], type = type, cellsize = 1000 * .tt_unit[[f]], quiet = TRUE)))
    names(gr) <- names(B)
    for (f in names(B)) {
      a <- as.numeric(sf::st_area(gr[[f]])) / .tt_unit[[f]]^2
      expect_gt(min(a), 1e3)
      expect_setequal(.tt_key(gr[[f]]), .tt_key(gr$m))
    }
  }
  # A side one metre over a whole number of cells keeps its column of parts
  # one metre wide.
  wide <- .tt_rect(4.5e5, 4.55e6, 9001, 6000)
  gr <- .tt_quiet(create_grid_polygons(wide, type = "square", cellsize = 1000, quiet = TRUE))
  expect_identical(nrow(gr), 60L)
  expect_equal(min(as.numeric(sf::st_area(gr))), 1000, tolerance = 1e-6)
  # A control box that is no whole number of cells: as before.
  for (type in c("square", "hex")) {
    gr <- .tt_quiet(create_grid_polygons(.tt_rect(4.5e5, 4.55e6, 9100, 6100), type = type,
                                         cellsize = 1000, quiet = TRUE))
    expect_identical(nrow(gr), if (type == "square") 70L else 80L)
  }
})

test_that("the help pages state the tie rules of the tessellation", {
  page <- function(name) {
    file <- test_path("..", "..", "man", name)
    rd <- if (file.exists(file)) tools::parse_Rd(file)
          else tools::Rd_db("spatialkit")[[name]]
    expect_s3_class(rd, "Rd")
    gsub("[[:space:]]+", " ", paste(unlist(rd), collapse = ""))
  }
  bt <- page("build_tessellation.Rd")
  expect_match(bt, paste0("a point in two or more cells (on a shared edge, or within that ",
                          "margin of it) takes the smallest cell_id among them"), fixed = TRUE)
  expect_match(bt, paste0("a point the same distance from two cells, to one part in ten ",
                          "million, takes the smaller cell_id"), fixed = TRUE)
  expect_match(bt, "so they are split the same way in any unit or origin", fixed = TRUE)
  expect_match(bt, "points within about a billionth of the extent of each other count as coincident",
               fixed = TRUE)
  expect_match(page("create_voronoi_polygons.Rd"),
               "points within about a billionth of the extent of each other count as coincident",
               fixed = TRUE)
  expect_match(page("create_grid_polygons.Rd"),
               paste0("exceeds a whole number of square cells, of this size or of the size ",
                      "derived from target_cells, by no more than a margin of one part in ",
                      "1e10 of the largest coordinate of the boundary"), fixed = TRUE)
  expect_match(bt, paste0("exceeds a whole number of square cells by no more than a margin of ",
                          "one part in 1e10 of the largest coordinate of the boundary"),
               fixed = TRUE)
})
