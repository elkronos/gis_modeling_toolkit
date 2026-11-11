# Two points that differ only in the last bits of their coordinates (one
# station recorded twice) are the same distance from every grid point and
# every location asked for (BUG-216).  Point covariates of predict_surface()
# and the scored locations of cell_values_at() are ranked by column, then by
# row read the same way on y, with a step that is never less than one part
# in 1e10 of the largest coordinate, so those bits, which move with the unit
# and the origin, do not choose between the two.  predict_surface() takes
# the first of them in row order; cell_values_at() fills from the cell first
# in ID order, as it does for rows at one location.

.cp_crs_m    <- function() sf::st_crs(26915)
.cp_crs_usft <- function() sf::st_crs("+proj=utm +zone=15 +datum=NAD83 +units=us-ft +no_defs")
.cp_crs_m01  <- function() sf::st_crs(paste0("+proj=tmerc +lat_0=0 +lon_0=-93 +k=0.9996 ",
                                            "+x_0=500000.1 +y_0=0.1 +datum=NAD83 +units=m +no_defs"))
.cp_frames <- function(pm) list(m = pm, usft = sf::st_transform(pm, .cp_crs_usft()),
                                m01 = sf::st_transform(pm, .cp_crs_m01()))
.cp_unit <- c(m = 1, usft = 3937 / 1200, m01 = 1)

.cp_quiet <- function(expr) suppressWarnings(suppressMessages(
  logger::with_log_threshold(expr, threshold = logger::FATAL, namespace = "spatialkit",
                             index = 2)))

.cp_x0 <- 4.5e5
.cp_y0 <- 4.55e6

.cp_pts <- function(d) sf::st_as_sf(d, coords = c("x", "y"), crs = .cp_crs_m())

# Two points at position i near (450500, 4550500), the second 5.82e-11 m east
# of the first and dy south of it (one ulp of y is 9.31e-10 m there), w 1 and
# 2.  `far` adds points that widen the layer: "one" a point 6.4 km away, w 3;
# "lattice" a 5 x 4 lattice at 250 m from 5 km away, w 3 and on, so the
# search goes to FNN where it is installed.
.cp_pair <- function(i, dy, far = c("none", "one", "lattice")) {
  far <- match.arg(far)
  xa <- .cp_x0 + 500 + i * 0.37; ya <- .cp_y0 + 500 + i * 0.29
  d <- data.frame(x = c(xa, xa + 5.82e-11), y = c(ya, ya + dy), w = c(1, 2))
  if (far == "one")
    d <- rbind(d, data.frame(x = .cp_x0 + 5000, y = .cp_y0 + 5000, w = 3))
  if (far == "lattice") {
    g <- expand.grid(i = 0:4, j = 0:3)
    d <- rbind(d, data.frame(x = .cp_x0 + 5000 + 250 * g$i, y = .cp_y0 + 5000 + 250 * g$j,
                             w = 2 + seq_len(nrow(g))))
  }
  .cp_pts(d)
}

# The layouts whose row moved between the frames before: two points only,
# 1.86e-9 m apart in y (the points' own box is last-bits small, so the last
# bits of x parted them into two columns), and with points far away,
# 9.31e-10 m apart in y (one column, and the raw y chose).
.cp_cases <- function() list(
  list(i = 1, dy = -2 * 9.31e-10, far = "none"),
  list(i = 3, dy = -2 * 9.31e-10, far = "none"),
  list(i = 1, dy = -9.31e-10, far = "one"),
  list(i = 3, dy = -9.31e-10, far = "one"),
  list(i = 6, dy = -9.31e-10, far = "one"),
  list(i = 1, dy = -9.31e-10, far = "lattice"),
  list(i = 3, dy = -9.31e-10, far = "lattice"),
  list(i = 6, dy = -9.31e-10, far = "lattice"))

# The search of .surface_nearest_rows() for one grid point, in each frame,
# for the layer in the given row order: the row it takes.
.cp_rows <- function(case, perm = NULL) {
  L  <- .cp_frames(.cp_pair(case$i, case$dy, case$far))
  gp <- .cp_frames(.cp_pts(data.frame(x = .cp_x0 + 100, y = .cp_y0 + 100)))
  vapply(names(L), function(f) {
    lay <- if (is.null(perm)) L[[f]] else L[[f]][perm, ]
    .surface_nearest_rows(sf::st_geometry(gp[[f]]), lay)
  }, integer(1))
}

test_that(".surface_nearest_rows() gives a tie between two points within the last bits of each other to the first in row order in every frame", {
  for (case in .cp_cases()) {
    lab <- sprintf("i %g, dy %g, far %s", case$i, case$dy, case$far)
    n <- nrow(.cp_pair(case$i, case$dy, case$far))
    # Before, the rows moved between the frames: 1, 1, 2 for the pair alone,
    # 2, 2, 1 with a point or a lattice far away.
    expect_identical(.cp_rows(case), c(m = 1L, usft = 1L, m01 = 1L), info = lab)
    # The layer reversed: its first of the two is the second point.
    expect_identical(.cp_rows(case, perm = n:1), c(m = n - 1L, usft = n - 1L, m01 = n - 1L),
                     info = lab)
    # Shuffled: the first of the two in the shuffled order.
    perm <- withr::with_seed(case$i, sample(n))
    want <- min(which(perm %in% 1:2))
    expect_identical(.cp_rows(case, perm = perm), c(m = want, usft = want, m01 = want),
                     info = lab)
  }
  # Without FNN the search is exact, and gives the same rows.
  local_mocked_bindings(.has_backend = function(pkg, ...) FALSE)
  for (case in .cp_cases()) {
    lab <- sprintf("no FNN: i %g, dy %g, far %s", case$i, case$dy, case$far)
    expect_identical(.cp_rows(case), c(m = 1L, usft = 1L, m01 = 1L), info = lab)
  }
})

test_that("predict_surface() takes the same covariates beside two points within the last bits of each other in every frame", {
  g <- expand.grid(i = 0:19, j = 0:14)
  x <- .cp_x0 + g$i * 250; y <- .cp_y0 + g$j * 250
  tr <- .cp_frames(.cp_pts(data.frame(
    x = x, y = y,
    z = withr::with_seed(1, sin((x - min(x)) / 1500) + stats::rnorm(length(x), 0, 0.5)),
    w = withr::with_seed(2, (x - min(x)) / 1000 + stats::rnorm(length(x))))))
  gp <- .cp_frames(.cp_pts(data.frame(x = .cp_x0 + c(100, 900, 2000),
                                      y = .cp_y0 + c(100, 300, 2000))))
  surf <- function(f, cov) {
    fit <- lm_spatial_fit(tr[[f]], "z", "w")
    .cp_quiet(predict_surface(fit, grid = gp[[f]], covariates = cov[, "w"]))
  }
  for (case in .cp_cases()[1:5]) {
    lab <- sprintf("i %g, dy %g, far %s", case$i, case$dy, case$far)
    L <- .cp_frames(.cp_pair(case$i, case$dy, case$far))
    ref <- surf("m", L$m)
    for (f in names(L)) {
      p <- surf(f, L[[f]])
      # Every grid point is nearer the pair than the far point.
      want <- c(1, 1, 1)
      # Before, w was 2 and .pred 0.7308 with the origin moved, 1 and 0.7378
      # in metres and US feet, at the pair alone.
      expect_identical(p$w, want, info = paste(lab, f))
      # z ~ w fitted on the same numbers in each frame; to rounding.
      expect_equal(p$.pred, ref$.pred, tolerance = 1e-12, info = paste(lab, f))
      # The layer reversed takes the second point.
      r <- surf(f, L[[f]][rev(seq_len(nrow(L[[f]]))), ])
      expect_identical(r$w, c(2, 2, 2), info = paste(lab, f))
    }
  }
})

# Scored rows of cell_values_at() in two 1 km cells side by side: one
# station on their shared edge, recorded twice, the first record (row 1) in
# cell 2 and the second, one ulp (9.31e-10 m) north of it, in cell 1; ten
# more rows in the south of each cell.  A location 2.5 km north of the
# station lies in no cell and is filled from the nearest scored location.
.cp_cells <- function() {
  sq <- function(a, b) sf::st_polygon(list(rbind(c(a, b), c(a + 1000, b), c(a + 1000, b + 1000),
                                                c(a, b + 1000), c(a, b))))
  sf::st_sf(poly_id = 1:2, geometry = sf::st_sfc(sq(.cp_x0, .cp_y0), sq(.cp_x0 + 1000, .cp_y0),
                                                 crs = .cp_crs_m()))
}

.cp_station <- function(i) {
  sx <- .cp_x0 + 1000; sy <- .cp_y0 + 500 + i * 0.29
  k <- 0:9
  rest <- data.frame(x = c(.cp_x0 + 50 + 90 * k, .cp_x0 + 1050 + 90 * k),
                     y = .cp_y0 + 50 + 30 * c(k, k),
                     poly_id = rep(1:2, each = 10L),
                     z = c(10 + (k - 4.5) / 10, 20 + (k - 4.5) / 10))
  list(pts = .cp_pts(rbind(data.frame(x = c(sx, sx), y = c(sy, sy + 9.31e-10),
                                      poly_id = c(2L, 1L), z = c(20, 10)), rest)),
       q = .cp_pts(data.frame(x = sx, y = sy + 2500)))
}

.cp_fill <- function(lay, perm = NULL) {
  pts <- .cp_frames(lay$pts); q <- .cp_frames(lay$q); cells <- .cp_frames(.cp_cells())
  lapply(stats::setNames(names(pts), names(pts)), function(f) {
    p <- if (is.null(perm)) pts[[f]] else pts[[f]][perm, ]
    v <- .cp_quiet(summarize_by_cell(p, response_var = "z", cells_sf = cells[[f]]))
    .cp_quiet(cell_values_at(q[[f]], v, p, "z"))
  })
}

.cp_check_fill <- function(res, lab) {
  for (f in names(res)) {
    # Before, the fill came from cell 2 in metres and US feet and from cell 1
    # with the origin moved, at positions 1 and 4.
    expect_identical(res[[f]]$filled_from, 1L, info = paste(lab, f))
    # The cell's value is pooled from coordinates read in each frame.
    expect_equal(res[[f]]$value, res$m$value, tolerance = 1e-9, info = paste(lab, f))
    expect_equal(res[[f]]$dist_nearest / .cp_unit[[f]], res$m$dist_nearest,
                 tolerance = 1e-9, info = paste(lab, f))
  }
}

test_that("cell_values_at() fills beside two scored locations within the last bits of each other from the cell first in ID order in every frame", {
  # The two records at one point, rows at one location: the cell first in
  # ID order, as the help says.
  same <- .cp_station(1)
  sf::st_geometry(same$pts)[[2]] <- sf::st_geometry(same$pts)[[1]]
  ref <- .cp_fill(same)$m
  expect_identical(ref$filled_from, 1L)
  for (i in c(1, 4)) {
    lay <- .cp_station(i)
    res <- .cp_fill(lay)
    .cp_check_fill(res, sprintf("position %g", i))
    expect_equal(res$m$value, ref$value, tolerance = 1e-9)
    n <- nrow(lay$pts)
    .cp_check_fill(.cp_fill(lay, perm = n:1), sprintf("position %g, reversed", i))
    .cp_check_fill(.cp_fill(lay, perm = withr::with_seed(i, sample(n))),
                   sprintf("position %g, shuffled", i))
  }
  # Without FNN the search is exact, and gives the same fill.
  local_mocked_bindings(.has_backend = function(pkg, ...) FALSE)
  for (i in c(1, 4)) .cp_check_fill(.cp_fill(.cp_station(i)), sprintf("no FNN: position %g", i))
})

test_that("cell_values_at() fills beside two scored locations within the last bits of each other with a scored row far away", {
  # The two records alone in cells 1 and 2 (row 1 in cell 1), and a third
  # scored row in a cell 5 km away, so the scored locations share a column
  # and the raw y chose between them before.
  sq <- function(a, b) sf::st_polygon(list(rbind(c(a, b), c(a + 1000, b), c(a + 1000, b + 1000),
                                                c(a, b + 1000), c(a, b))))
  cells <- .cp_frames(sf::st_sf(poly_id = 1:3, geometry = sf::st_sfc(
    sq(.cp_x0, .cp_y0), sq(.cp_x0 + 1000, .cp_y0), sq(.cp_x0 + 5000, .cp_y0 + 5000),
    crs = .cp_crs_m())))
  q <- .cp_frames(.cp_pts(data.frame(x = .cp_x0 + 600, y = .cp_y0 + 3000)))
  for (i in c(3, 5, 6)) {
    xa <- .cp_x0 + 1000; ya <- .cp_y0 + 500 + i * 0.29
    pts <- .cp_frames(.cp_pts(data.frame(x = c(xa + 5.82e-11, xa, .cp_x0 + 5500),
                                         y = c(ya - 9.31e-10, ya, .cp_y0 + 5500),
                                         poly_id = c(2L, 1L, 3L), z = c(20, 10, 30))))
    for (perm in list(1:3, 3:1, c(2L, 3L, 1L))) {
      res <- lapply(names(pts), function(f) {
        p <- pts[[f]][perm, ]
        v <- .cp_quiet(summarize_by_cell(p, response_var = "z", cells_sf = cells[[f]]))
        .cp_quiet(cell_values_at(q[[f]], v, p, "z"))$filled_from
      })
      # Before, cell 2 in metres and US feet and cell 1 with the origin moved.
      expect_identical(unlist(res), c(1L, 1L, 1L),
                       info = sprintf("position %g, rows %s", i, paste(perm, collapse = " ")))
    }
  }
})

test_that(".near_points_extent() is the diagonal, or the margin of an edge over .TIE_TOL", {
  # A layer wider than a thousandth of its largest coordinate: the diagonal.
  expect_identical(.near_points_extent(c(4.5e5, 4.53e5), c(4.55e6, 4.554e6)), 5000)
  # Two points within the last bits of each other: the floor, the margin of
  # .edge_tol(), which is four units in the last place of the largest
  # coordinate where the diagonal is last-bits small (BUG-225: no longer one
  # part in 1e10 of it, 455 m of extent and a step of 0.455 mm).  The two
  # points, one unit in the last place apart in x and in y, are still one
  # column.
  ulp <- 2^(floor(log2(4.55e6)) - 52)
  ex <- .near_points_extent(c(4.5e5, 4.5e5 + 5.82e-11), c(4.55e6, 4.55e6 - 9.31e-10))
  expect_equal(ex, 4 * ulp / .TIE_TOL, tolerance = 1e-12)
  expect_identical(.near_columns(c(4.5e5, 4.5e5 + 5.82e-11), c(0, 0), ex), c(1L, 1L))
  expect_identical(.near_columns(c(4.55e6, 4.55e6 - 9.31e-10), c(0, 0), ex), c(1L, 1L))
  # A layer at the origin, and one point.
  expect_identical(.near_points_extent(c(0, 0), c(0, 0)), 0)
  expect_identical(.near_points_extent(5, 7), 4 * 2^(floor(log2(7)) - 52) / .TIE_TOL)
})

# Cells 1 and 2 side by side and cell 3 5 km east of cell 1, and the fill
# of a location 2.5 km north of cell 1 in each frame from scored rows `d`
# (columns x, y, poly_id, z) in the given row order.
.cp_fill3 <- function(d, perm) {
  sq <- function(a, b) sf::st_polygon(list(rbind(c(a, b), c(a + 1000, b), c(a + 1000, b + 1000),
                                                c(a, b + 1000), c(a, b))))
  cells <- .cp_frames(sf::st_sf(poly_id = 1:3, geometry = sf::st_sfc(
    sq(.cp_x0, .cp_y0), sq(.cp_x0 + 1000, .cp_y0), sq(.cp_x0 + 5000, .cp_y0), crs = .cp_crs_m())))
  q <- .cp_frames(.cp_pts(data.frame(x = .cp_x0 + 600, y = .cp_y0 + 3000)))
  pts <- .cp_frames(.cp_pts(d[perm, ]))
  vapply(names(pts), function(f) {
    v <- .cp_quiet(summarize_by_cell(pts[[f]], response_var = "z", cells_sf = cells[[f]]))
    .cp_quiet(cell_values_at(q[[f]], v, pts[[f]], "z"))$filled_from
  }, integer(1))
}

test_that("cell_values_at() fills from the cell first in ID order where only the two records of one station are scored", {
  # The two records alone, row 1 in cell 2 and row 2 in cell 1, 5.82e-11 m
  # east and one ulp south of row 1: the scored locations' own box is
  # last-bits small, so only the floor of the extent puts them in one column
  # and one row.  Before (and with the diagonal alone as the extent), cell 2
  # in metres and US feet and cell 1 with the origin moved.  Two scored
  # locations are searched exactly with or without FNN; both are run.
  for (no_fnn in c(FALSE, TRUE)) {
    if (no_fnn) local_mocked_bindings(.has_backend = function(pkg, ...) FALSE)
    for (i in c(1, 4)) {
      xa <- .cp_x0 + 1000 + i * 0.37; ya <- .cp_y0 + 500 + i * 0.29
      d <- data.frame(x = c(xa, xa + 5.82e-11), y = c(ya, ya - 9.31e-10),
                      poly_id = c(2L, 1L), z = c(20, 10))
      for (perm in list(1:2, 2:1))
        expect_identical(.cp_fill3(d, perm), c(m = 1L, usft = 1L, m01 = 1L),
                         info = sprintf("no FNN %s, position %g, rows %s", no_fnn, i,
                                        paste(perm, collapse = " ")))
    }
  }
})

test_that("cell_values_at() numbers two locations of one column and one row by the first cell among the rows of each", {
  # Location A holds rows in cells 3 and 1; location B, 5.82e-11 m east and
  # one ulp south of A, a row in cell 2; a fourth row lies in cell 3, 4.5 km
  # away.  The first cell in ID order among A's rows is 1, before B's 2, so
  # the fill comes from cell 1, as for rows at one location.  Numbered by
  # the last cell of each location instead (3 for A), B came first and the
  # fill came from cell 2, and from cell 1 with the origin moved at
  # position 3; before the fix the same.
  for (i in c(1, 3)) {
    xa <- .cp_x0 + 1000 + i * 0.37; ya <- .cp_y0 + 500 + i * 0.29
    d <- data.frame(x = c(xa, xa, xa + 5.82e-11, .cp_x0 + 5500),
                    y = c(ya, ya, ya - 9.31e-10, .cp_y0 + 500),
                    poly_id = c(3L, 1L, 2L, 3L), z = c(30, 10, 20, 31))
    for (perm in list(1:4, 4:1, c(3L, 1L, 4L, 2L)))
      expect_identical(.cp_fill3(d, perm), c(m = 1L, usft = 1L, m01 = 1L),
                       info = sprintf("position %g, rows %s", i, paste(perm, collapse = " ")))
  }
})

test_that("?predict_surface states which grid points beside two points within the last bits of each other take the same covariates in every frame", {
  file <- testthat::test_path("..", "..", "man", "predict_surface.Rd")
  rd <- if (file.exists(file)) tools::parse_Rd(file) else tools::Rd_db("spatialkit")[["predict_surface.Rd"]]
  txt <- gsub("[[:space:]]+", " ", paste(unlist(rd), collapse = ""))
  # A grid point that ties the two by the rule takes the same row in every
  # frame; one a few centimetres from them need not tie, and the nearer by
  # its coordinates, which move with the unit and the origin, takes it.
  expect_match(txt, paste0("as does a grid point the same distance from two points that differ ",
                           "only in the last bits of their coordinates"), fixed = TRUE)
  expect_match(txt, paste0("closer to such a pair than about twenty million times the distance ",
                           "between them (2 cm for points one billionth of a metre apart) need ",
                           "not be the same distance from both by that rule, and takes the ",
                           "nearer as its coordinates give it, which can change with the unit ",
                           "or the origin"), fixed = TRUE)
  expect_no_match(txt, "or over two points that differ only in the last bits", fixed = TRUE)
})
