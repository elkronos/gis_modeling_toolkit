# A point a hair from two records of one station that differ in the last
# bits of their coordinates (BUG-268).  The nearest-centre search counted
# two centres as a tie when their squared distances agreed to one part in
# ten million, a margin relative to the least squared distance alone.  Two
# records `delta` apart have squared distances from a point `r` away that
# differ by about 2 r delta, within that margin only while r is more than
# about twenty million times delta (2 cm for records one billionth of a
# metre apart).  Closer, the nearer by the last bits was taken, and those
# bits move with the unit and the origin.  The tie now adds the length
# margin of the centres' box (.edge_tol()): a squared distance d2 ties with
# the least, best, when d2 <= best * (1 + .TIE_TOL) + tol * (2 sqrt(best) + tol).

.np_crs_m    <- function() sf::st_crs(26915)
.np_crs_usft <- function() sf::st_crs("+proj=utm +zone=15 +datum=NAD83 +units=us-ft +no_defs")
.np_crs_m01  <- function() sf::st_crs(paste0("+proj=tmerc +lat_0=0 +lon_0=-93 +k=0.9996 ",
                                            "+x_0=500000.1 +y_0=0.1 +datum=NAD83 +units=m +no_defs"))
.np_frames <- function(pm) list(m = pm, usft = sf::st_transform(pm, .np_crs_usft()),
                                m01 = sf::st_transform(pm, .np_crs_m01()))
.np_unit <- c(m = 1, usft = 3937 / 1200, m01 = 1)

.np_quiet <- function(expr) suppressWarnings(suppressMessages(
  logger::with_log_threshold(expr, threshold = logger::FATAL, namespace = "spatialkit",
                             index = 2)))

.np_pts <- function(d) sf::st_as_sf(d, coords = c("x", "y"), crs = .np_crs_m())

# One unit in the last place of y near 4.55e6 m.
.np_uy <- 9.313225746154785e-10

# Two records of one station at position i near (450500, 4550500), the
# second one ulp north of the first, w 1 and 2.  `far` adds points that
# widen the layer: "one" a point 7 km away, w 3; "lattice" a 5 x 4 lattice
# at 250 m from 5 km away, w 3 and on, so that there are more than 16
# points and the FNN shortlist runs where FNN is installed.
.np_pair <- function(i, far = c("none", "one", "lattice")) {
  far <- match.arg(far)
  xa <- 450500 + i * 0.37; ya <- 4550500 + i * 0.29
  d <- data.frame(x = c(xa, xa), y = c(ya, ya + .np_uy), w = c(1, 2))
  if (far == "one")
    d <- rbind(d, data.frame(x = 455500, y = 4555500, w = 3))
  if (far == "lattice") {
    g <- expand.grid(i = 0:4, j = 0:3)
    d <- rbind(d, data.frame(x = 455500 + 250 * g$i, y = 4555500 + 250 * g$j,
                             w = 2 + seq_len(nrow(g))))
  }
  list(layer = .np_pts(d), x = xa, y = ya)
}

# Eight grid points on a circle of radius r about the first record.
.np_circle <- function(x, y, r) {
  th <- seq(0, 2 * pi, length.out = 9)[-9]
  .np_pts(data.frame(x = x + r * cos(th), y = y + r * sin(th)))
}

# Positions whose grid points took another row between the frames before
# the fix, at 1 mm and at 2 cm, with and without points far away
# (repro/positions.R of the batch: 16 and 32 of 40 positions moved).
.np_positions <- c(5, 10, 14, 18)

.np_rows <- function(i, r, far, perm = NULL) {
  p  <- .np_pair(i, far)
  L  <- .np_frames(p$layer)
  G  <- .np_frames(.np_circle(p$x, p$y, r))
  vapply(names(L), function(f) {
    lay <- if (is.null(perm)) L[[f]] else L[[f]][perm, ]
    .surface_nearest_rows(sf::st_geometry(G[[f]]), lay)
  }, integer(8))
}

test_that(".surface_nearest_rows() gives a grid point a hair from two records of one station the first of them in every frame", {
  for (no_fnn in c(FALSE, TRUE)) {
    if (no_fnn) local_mocked_bindings(.has_backend = function(pkg, ...) FALSE)
    for (far in c("none", "one", "lattice")) for (r in c(1e-3, 2e-2)) for (i in .np_positions) {
      lab <- sprintf("no FNN %s, far %s, r %g, position %g", no_fnn, far, r, i)
      n <- nrow(.np_pair(i, far)$layer)
      # Before, the grid points nearer the second record by the last bits
      # took row 2 in one frame and row 1 in another.
      rows <- .np_rows(i, r, far)
      expect_true(all(rows == 1L), info = lab)
      # The layer reversed: the first of the two is now the second record,
      # at row n - 1.
      expect_true(all(.np_rows(i, r, far, perm = n:1) == n - 1L), info = lab)
    }
  }
})

test_that("predict_surface() takes the same covariates a hair from two records of one station in every frame", {
  g <- expand.grid(i = 0:19, j = 0:14)
  x <- 4.5e5 + g$i * 250; y <- 4.55e6 + g$j * 250
  tr <- .np_frames(.np_pts(data.frame(
    x = x, y = y,
    z = withr::with_seed(1, sin((x - min(x)) / 1500) + stats::rnorm(length(x), 0, 0.5)),
    w = withr::with_seed(2, (x - min(x)) / 1000 + stats::rnorm(length(x))))))
  for (far in c("none", "lattice")) for (i in .np_positions[1:2]) for (r in c(1e-3, 2e-2)) {
    lab <- sprintf("far %s, position %g, r %g", far, i, r)
    p <- .np_pair(i, far)
    L <- .np_frames(p$layer)
    G <- .np_frames(.np_circle(p$x, p$y, r))
    surf <- function(f, cov) {
      fit <- lm_spatial_fit(tr[[f]], "z", "w")
      .np_quiet(predict_surface(fit, grid = G[[f]], covariates = cov[, "w"]))
    }
    ref <- surf("m", L$m)
    for (f in names(L)) {
      s <- surf(f, L[[f]])
      # Before, some grid points took w 2 in one frame and w 1 in another.
      expect_identical(s$w, rep(1, 8), info = paste(lab, f))
      # z ~ w fitted on the same numbers in each frame; to rounding.
      expect_equal(s$.pred, ref$.pred, tolerance = 1e-12, info = paste(lab, f))
    }
  }
})

# Scored rows of cell_values_at() in two 1 km cells side by side, and a
# station 1.5 km north of their shared edge, outside both, recorded twice:
# the first record in cell 2, the second, one ulp north, in cell 1; ten
# more rows in the south of each cell, so 22 scored locations, more than
# 16, and the FNN shortlist runs where FNN is installed.  Locations 1 mm
# from the station lie in no cell and are filled from the nearest scored
# location (the probe of BUG-268).
.np_station <- function(i) {
  x0 <- 450000; y0 <- 4550000
  sq <- function(a, b) sf::st_polygon(list(rbind(c(a, b), c(a + 1000, b), c(a + 1000, b + 1000),
                                                c(a, b + 1000), c(a, b))))
  cells <- sf::st_sf(poly_id = 1:2, geometry = sf::st_sfc(sq(x0, y0), sq(x0 + 1000, y0),
                                                          crs = .np_crs_m()))
  sx <- x0 + 1000 + i * 0.37; sy <- y0 + 2500 + i * 0.29
  k <- 0:9
  d <- rbind(data.frame(x = c(sx, sx), y = c(sy, sy + .np_uy), poly_id = c(2L, 1L), z = c(20, 10)),
             data.frame(x = c(x0 + 50 + 90 * k, x0 + 1050 + 90 * k), y = y0 + 50 + 30 * c(k, k),
                        poly_id = rep(1:2, each = 10L),
                        z = c(10 + (k - 4.5) / 10, 20 + (k - 4.5) / 10)))
  list(pts = .np_pts(d), cells = cells, q = .np_circle(sx, sy, 1e-3))
}

.np_fill <- function(lay, perm = NULL) {
  pts <- .np_frames(lay$pts); q <- .np_frames(lay$q); cells <- .np_frames(lay$cells)
  lapply(stats::setNames(names(pts), names(pts)), function(f) {
    p <- if (is.null(perm)) pts[[f]] else pts[[f]][perm, ]
    v <- .np_quiet(summarize_by_cell(p, response_var = "z", cells_sf = cells[[f]]))
    .np_quiet(cell_values_at(q[[f]], v, p, "z"))
  })
}

.np_check_fill <- function(res, lab) {
  for (f in names(res)) {
    # Before, locations south of the station were filled from cell 2 in
    # metres and US feet and from cell 1 with the origin moved, at
    # positions 5, 9 and 11.  The two records tie; the cell first in ID
    # order fills, as for rows at one location.
    expect_identical(res[[f]]$filled_from, rep(1L, 8), info = paste(lab, f))
    # The cell's value is pooled from coordinates read in each frame.
    expect_equal(res[[f]]$value, res$m$value, tolerance = 1e-9, info = paste(lab, f))
    # A distance of 1 mm from coordinates near 4.55e6 m that PROJ moved to
    # another frame and back: each coordinate is rounded there by about
    # 1e-9 m, a millionth of the distance, so it is compared to 1e-5 of
    # it (1e-8 m), not to the last bits.
    expect_equal(res[[f]]$dist_nearest / .np_unit[[f]], res$m$dist_nearest,
                 tolerance = 1e-5, info = paste(lab, f))
  }
}

test_that("cell_values_at() fills a location a hair from two records of one station from the cell first in ID order in every frame", {
  for (no_fnn in c(FALSE, TRUE)) {
    if (no_fnn) local_mocked_bindings(.has_backend = function(pkg, ...) FALSE)
    for (i in c(5, 9, 11)) {
      lay <- .np_station(i)
      lab <- sprintf("no FNN %s, position %g", no_fnn, i)
      .np_check_fill(.np_fill(lay), lab)
      n <- nrow(lay$pts)
      .np_check_fill(.np_fill(lay, perm = n:1), paste(lab, "reversed"))
    }
  }
})

test_that("the margin closes the gap at the pair and is the relative rule to the bit at 0", {
  xa <- 450505; ya <- 4550505
  cx <- c(xa, xa); cy <- c(ya, ya + .np_uy)
  # A point 1 mm south: nearer the first by the last bits; 1 mm north,
  # nearer the second.  Both tie, and take centre 1.
  px <- c(xa, xa); py <- c(ya - 1e-3, ya + 1e-3)
  expect_identical(.nearest_exact(px, py, cx, cy), c(1L, 1L))
  # With no margin the rule is the relative one, which parted them.
  expect_identical(.nearest_exact(px, py, cx, cy, tol = 0), c(1L, 2L))
  expect_identical(.tie_bound(4, 0), 4 * (1 + .TIE_TOL))
  expect_identical(.tie_bound(4, 0.5), 4 * (1 + .TIE_TOL) + 0.5 * (2 * 2 + 0.5))
  # The margin is .edge_tol() of the centres that have a location.
  expect_identical(.centres_tol(cx, cy), .edge_tol(.xy_box(cx, cy)))
  expect_identical(.centres_tol(c(cx, NA), c(cy, 1)), .edge_tol(.xy_box(cx, cy)))
  expect_identical(.centres_tol(NA_real_, NA_real_), 0)
  # Centres reached through FNN and through the shortlist of
  # .values_at_nearest(), one ulp apart, among more than 16.
  skip_if_not_installed("FNN")
  ax <- c(cx, 455500 + 250 * (0:19 %% 5)); ay <- c(cy, 4555500 + 250 * (0:19 %/% 5))
  expect_identical(.nearest_centre(px, py, ax, ay, use_fnn = TRUE), c(1L, 1L))
  expect_identical(.values_at_nearest(px, py, ax, ay, TRUE), c(1L, 1L))
  expect_identical(.values_at_nearest(px, py, ax, ay, FALSE), c(1L, 1L))
})

test_that("the three searches give one answer on a 200 x 200 lattice, with and without FNN", {
  skip_if_not_installed("FNN")
  q <- expand.grid(x = seq(0, 1000, length.out = 200), y = seq(0, 1000, length.out = 200))
  qx <- q$x + 4.5e5; qy <- q$y + 4.55e6
  # Random centres, lattice centres on the points' own steps, and lattice
  # centres with every other one recorded twice, one ulp north (added at the
  # coordinates' size, where the ulp is 9.31e-10 m).
  rc <- withr::with_seed(268, cbind(stats::runif(40, 0, 1000), stats::runif(40, 0, 1000)))
  rc <- cbind(rc[, 1L] + 4.5e5, rc[, 2L] + 4.55e6)
  lc <- as.matrix(expand.grid(seq(0, 1000, length.out = 8), seq(0, 1000, length.out = 6)))
  lc <- cbind(lc[, 1L] + 4.5e5, lc[, 2L] + 4.55e6)
  half <- lc[c(TRUE, FALSE), , drop = FALSE]
  dc <- rbind(lc, cbind(half[, 1L], half[, 2L] + .np_uy))
  expect_true(all(dc[-seq_len(nrow(lc)), 2L] != half[, 2L]))
  for (nm in c("random", "lattice", "lattice twice")) {
    cen <- switch(nm, random = rc, lattice = lc, `lattice twice` = dc)
    cx <- cen[, 1L]; cy <- cen[, 2L]
    exact <- .nearest_exact(qx, qy, cx, cy)
    expect_identical(.nearest_centre(qx, qy, cx, cy, use_fnn = TRUE), exact, info = nm)
    expect_identical(.nearest_centre(qx, qy, cx, cy, use_fnn = FALSE), exact, info = nm)
    expect_identical(.values_at_nearest(qx, qy, cx, cy, TRUE), exact, info = nm)
    expect_identical(.values_at_nearest(qx, qy, cx, cy, FALSE), exact, info = nm)
  }
  # A centre recorded twice never takes a point as its second record.
  twice <- nrow(lc) + seq_len(nrow(half))
  expect_false(any(.nearest_exact(qx, qy, dc[, 1L], dc[, 2L]) %in% twice))
})

test_that("?predict_surface and ?cell_values_at state the margin of the point tie", {
  read_rd <- function(name) {
    file <- testthat::test_path("..", "..", "man", paste0(name, ".Rd"))
    rd <- if (file.exists(file)) tools::parse_Rd(file) else tools::Rd_db("spatialkit")[[paste0(name, ".Rd")]]
    gsub("[[:space:]]+", " ", paste(unlist(rd), collapse = ""))
  }
  rule <- paste0("no more than the square of the least distance lengthened by a margin, ",
                 "plus one part in ten million of the least squared distance")
  ps <- read_rd("predict_surface")
  expect_match(ps, rule, fixed = TRUE)
  expect_match(ps, paste0("as does a grid point near two points that lie within that margin of ",
                          "each other"), fixed = TRUE)
  expect_no_match(ps, "twenty million", fixed = TRUE)
  cv <- read_rd("cell_values_at")
  expect_match(cv, rule, fixed = TRUE)
  expect_match(cv, "while they lie within the margin of each other", fixed = TRUE)
  for (nm in c("resolution_profile", "make_folds", "summarize_by_cell",
               "assign_features_to_polygons"))
    expect_match(read_rd(nm), "lengthened by a margin", fixed = TRUE, info = nm)
})

test_that("three records of one station take the first number among more than 16 centres, the shortlist included", {
  # A at the query circle's centre, B one ulp north and C one ulp east, in
  # orders that put A after the others, and 20 lattice centres 5 km away, so
  # the FNN shortlist of .values_at_nearest() holds two of the three at
  # first.  Its furthest member then lies past the relative bound but within
  # the margin, and a stop test without the margin ended the search before
  # the lowest-numbered record came in.
  ux <- 5.820766091346741e-11
  lx <- 455500 + 250 * (0:19 %% 5); ly <- 4555500 + 250 * (0:19 %/% 5)
  th <- seq(0, 2 * pi, length.out = 9)[-9]
  for (i in c(3, 17, 29)) {
    xa <- 450500 + i * 0.37; ya <- 4550500 + i * 0.29
    sx <- c(xa, xa, xa + ux); sy <- c(ya, ya + .np_uy, ya)
    for (ord in list(c(2, 3, 1), c(3, 2, 1))) for (r in c(1e-3, 2e-2)) {
      lab <- sprintf("position %g, order %s, r %g", i, paste(ord, collapse = ""), r)
      cx <- c(sx[ord], lx); cy <- c(sy[ord], ly)
      px <- xa + r * cos(th); py <- ya + r * sin(th)
      expect_identical(.nearest_exact(px, py, cx, cy), rep(1L, 8), info = lab)
      expect_identical(.nearest_centre(px, py, cx, cy, use_fnn = FALSE), rep(1L, 8), info = lab)
      expect_identical(.values_at_nearest(px, py, cx, cy, FALSE), rep(1L, 8), info = lab)
      if (requireNamespace("FNN", quietly = TRUE)) {
        expect_identical(.nearest_centre(px, py, cx, cy, use_fnn = TRUE), rep(1L, 8), info = lab)
        expect_identical(.values_at_nearest(px, py, cx, cy, TRUE), rep(1L, 8), info = lab)
      }
    }
  }
})

test_that("with no margin the bound stays finite where a squared distance overflows", {
  # Every centre at (0, 0): the margin is 0.  A point at 1e200 has a squared
  # distance of Inf; 0 * Inf would make the bound NaN and the centre NA.
  expect_identical(.tie_bound(Inf, 0), Inf)
  expect_identical(.nearest_exact(c(1e200, 1), c(0, 0), c(0, 0), c(0, 0)), c(1L, 1L))
})

# Two distinct centres 600 m apart in EPSG:26915 and a location off their
# bisector toward the second by `off` metres, 316 m from both: its distances
# differ by about 1.9 off.  The margin is one part in 1e10 of the largest
# coordinate: 4.55e-4 m at the UTM northing, 2.6e-7 m with the origin next
# to the data.
.np_bisector <- function() {
  x0 <- 450000; y0 <- 4550000
  list(x0 = x0, y0 = y0, cx = c(x0 + 700, x0 + 1300), cy = c(y0 + 2500, y0 + 2500),
       off = c(1e-5, 1e-4, 2e-4))
}

test_that("the margin moves with the origin, as the help states", {
  b <- .np_bisector()
  px <- b$x0 + 1000 + b$off; py <- rep(b$y0 + 2600, length(b$off))
  # Far from the origin the distances tie and the lower-numbered centre takes
  # the location; with the origin next to the data (an exact subtraction)
  # they do not, and the nearer does.  Differences of 1.9e-5 to 3.8e-4 m lie
  # well inside the first margin and far outside the second, so this does not
  # turn on the last bits.
  expect_identical(.nearest_centre(px, py, b$cx, b$cy, use_fnn = FALSE), c(1L, 1L, 1L))
  expect_identical(.nearest_centre(px - b$x0, py - b$y0, b$cx - b$x0, b$cy - b$y0,
                                   use_fnn = FALSE), c(2L, 2L, 2L))
  expect_gt(.centres_tol(b$cx, b$cy), 1.9 * max(b$off))
  expect_lt(.centres_tol(b$cx - b$x0, b$cy - b$y0), 1.9 * min(b$off))
  # cell_values_at() through PROJ: EPSG:26915 and the same transverse
  # Mercator with its false origin next to the data.
  local <- sf::st_crs(paste("+proj=tmerc +lat_0=0 +lon_0=-93 +k=0.9996 +x_0=51000",
                            "+y_0=-4549000 +datum=NAD83 +units=m +no_defs"))
  sq <- function(a, c0) sf::st_polygon(list(rbind(c(a, c0), c(a + 1000, c0), c(a + 1000, c0 + 1000),
                                                 c(a, c0 + 1000), c(a, c0))))
  cells <- sf::st_sf(poly_id = 1:2, geometry = sf::st_sfc(sq(b$x0, b$y0), sq(b$x0 + 1000, b$y0),
                                                          crs = .np_crs_m()))
  k <- 0:9
  d <- rbind(data.frame(x = b$cx, y = b$cy, poly_id = 1:2, z = c(10, 20)),
             data.frame(x = c(b$x0 + 50 + 90 * k, b$x0 + 1050 + 90 * k),
                        y = b$y0 + 50 + 30 * c(k, k), poly_id = rep(1:2, each = 10L),
                        z = c(10 + (k - 4.5) / 10, 20 + (k - 4.5) / 10)))
  P <- .np_pts(d); Q <- .np_pts(data.frame(x = px, y = py))
  fill <- function(crs) {
    p <- if (crs == .np_crs_m()) P else sf::st_transform(P, crs)
    q <- if (crs == .np_crs_m()) Q else sf::st_transform(Q, crs)
    cl <- if (crs == .np_crs_m()) cells else sf::st_transform(cells, crs)
    v <- .np_quiet(summarize_by_cell(p, response_var = "z", cells_sf = cl))
    .np_quiet(cell_values_at(q, v, p, "z"))$filled_from
  }
  expect_identical(fill(.np_crs_m()), c(1L, 1L, 1L))
  expect_identical(fill(local), c(2L, 2L, 2L))
})

test_that("the help pages and the vignette state the margin and its dependence on the origin", {
  read_rd <- function(name) {
    file <- testthat::test_path("..", "..", "man", paste0(name, ".Rd"))
    rd <- if (file.exists(file)) tools::parse_Rd(file) else tools::Rd_db("spatialkit")[[paste0(name, ".Rd")]]
    gsub("[[:space:]]+", " ", paste(unlist(rd), collapse = ""))
  }
  cv <- read_rd("cell_values_at")
  expect_match(cv, paste0("whose distances to two of them differ by less than the margin with one ",
                          "origin and by more with another is a tie with the first"), fixed = TRUE)
  expect_match(read_rd("resolution_profile"),
               paste0("a row whose distances to two centres differ by less than the margin with ",
                      "one origin and by more with another goes to the lower-numbered centre ",
                      "with the first and to the nearer with the second"), fixed = TRUE)
  expect_match(read_rd("predict_surface"),
               paste0("Both margins grow with the largest coordinate (about 0.46 mm at a UTM ",
                      "northing near 4,550 km"), fixed = TRUE)
  vig <- testthat::test_path("..", "..", "vignettes", "diagnostics.Rmd")
  skip_if_not(file.exists(vig), "the vignette sources are not here")
  txt <- gsub("[[:space:]]+", " ", paste(readLines(vig, warn = FALSE), collapse = " "))
  expect_no_match(txt, "whose squares agree to one part in ten million count as the same",
                  fixed = TRUE)
  expect_match(txt, "lengthened by a small margin", fixed = TRUE)
})
