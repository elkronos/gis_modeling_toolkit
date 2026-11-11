# predict_surface() with line and polygon covariates (BUG-200): a grid point
# on the edge two covariate polygons share, or at the shared end of two
# lines, is at distance 0 from both.  It takes the feature whose
# bounding-box centre is in the lowest column, then lowest in y, then the
# lower row, so the surface is the same in another unit, with another
# origin and with the rows of `covariates` in another order.

.se_crs_m    <- function() sf::st_crs(26915)
.se_crs_usft <- function() sf::st_crs("+proj=utm +zone=15 +datum=NAD83 +units=us-ft +no_defs")
.se_crs_m01  <- function() sf::st_crs(paste0("+proj=tmerc +lat_0=0 +lon_0=-93 +k=0.9996 ",
                                            "+x_0=500000.1 +y_0=0.1 +datum=NAD83 +units=m +no_defs"))
.se_frames <- function(pm) list(m = pm, usft = sf::st_transform(pm, .se_crs_usft()),
                                m01 = sf::st_transform(pm, .se_crs_m01()))

.se_quiet <- function(expr) suppressWarnings(suppressMessages(
  logger::with_log_threshold(expr, threshold = logger::FATAL, namespace = "spatialkit",
                             index = 2)))

# Training points: a 20 x 15 lattice at 250 m with z and w.
.se_train <- function() {
  g <- expand.grid(i = 0:19, j = 0:14)
  x <- 4.5e5 + g$i * 250; y <- 4.55e6 + g$j * 250
  z <- withr::with_seed(1, sin((x - min(x)) / 1500) + cos((y - min(y)) / 1250) +
                          stats::rnorm(length(x), 0, 0.5))
  w <- withr::with_seed(2, (x - min(x)) / 1000 + stats::rnorm(length(x)))
  sf::st_as_sf(data.frame(x = x, y = y, z = z, w = w), coords = c("x", "y"),
               crs = .se_crs_m())
}

.se_square <- function(x0, y0, s) sf::st_polygon(list(rbind(
  c(x0, y0), c(x0 + s, y0), c(x0 + s, y0 + s), c(x0, y0 + s), c(x0, y0))))

# 10 x 8 square cells of 500 m, w the row number, numbered a row at a time.
.se_cells <- function() {
  b <- expand.grid(ci = 0:9, cj = 0:7)
  sf::st_sf(w = seq_len(nrow(b)) + 0,
            geometry = sf::st_sfc(lapply(seq_len(nrow(b)), function(k)
              .se_square(4.5e5 + 500 * b$ci[k], 4.55e6 + 500 * b$cj[k], 500)),
              crs = .se_crs_m()))
}

# 72 grid points on the shared vertical edges of the cells, and the row of
# the cell to the left of each.
.se_edge_grid <- function() {
  g <- expand.grid(i = 1:9, j = 0:7)
  pts <- sf::st_as_sf(data.frame(x = 4.5e5 + 500 * g$i, y = 4.55e6 + 250 + 500 * g$j),
                      coords = c("x", "y"), crs = .se_crs_m())
  list(pts = pts, left = (g$i - 1) + 10 * g$j + 1)
}

.se_surface <- function(train, grid, cov) {
  fit <- lm_spatial_fit(train, "z", "w")
  .se_quiet(predict_surface(fit, grid = grid, covariates = cov[, "w"]))
}

test_that("predict_surface() gives a grid point on a shared cell edge the cell to its left in every frame", {
  cells <- .se_frames(.se_cells())
  eg <- .se_edge_grid()
  grid <- .se_frames(eg$pts)
  train <- .se_frames(.se_train())
  s <- lapply(names(cells), function(f) .se_surface(train[[f]], grid[[f]], cells[[f]]))
  names(s) <- names(cells)
  for (f in names(cells)) {
    expect_identical(s[[f]]$w, eg$left, info = f)
    expect_identical(s[[f]]$w, s$m$w, info = f)
  }
  # The fits differ by the last bits of their coordinates only through the
  # covariate, which is now the same; .pred is identical across frames.
  expect_identical(s$usft$.pred, s$m$.pred)
  expect_identical(s$m01$.pred, s$m$.pred)
})

test_that("predict_surface()'s edge rule does not read the row order of line or polygon covariates", {
  cells <- .se_cells()
  eg <- .se_edge_grid()
  train <- .se_train()
  a <- .se_surface(train, eg$pts, cells)
  for (s in 1:3) {
    o <- withr::with_seed(s, sample(nrow(cells)))
    b <- .se_surface(train, eg$pts, cells[o, ])
    expect_identical(b$w, a$w, info = paste("shuffle", s))
    expect_identical(b$.pred, a$.pred, info = paste("shuffle", s))
  }
  # In the frame with the origin moved as well.
  cm <- .se_frames(cells)$m01
  gm <- .se_frames(eg$pts)$m01
  tm <- .se_frames(train)$m01
  o <- withr::with_seed(7, sample(nrow(cells)))
  expect_identical(.se_surface(tm, gm, cm[o, ])$w, eg$left)
})

test_that("predict_surface() gives a grid point at the shared end of two lines the line to its left", {
  g <- expand.grid(ci = 0:9, cj = 0:7)
  segs <- lapply(seq_len(nrow(g)), function(k) sf::st_linestring(rbind(
    c(4.5e5 + 500 * g$ci[k], 4.55e6 + 250 + 500 * g$cj[k]),
    c(4.5e5 + 500 * (g$ci[k] + 1), 4.55e6 + 250 + 500 * g$cj[k]))))
  ln <- .se_frames(sf::st_sf(w = seq_len(nrow(g)) + 0,
                             geometry = sf::st_sfc(segs, crs = .se_crs_m())))
  eg <- .se_edge_grid()
  grid <- .se_frames(eg$pts)
  train <- .se_frames(.se_train())
  for (f in names(ln)) {
    expect_identical(.se_surface(train[[f]], grid[[f]], ln[[f]])$w, eg$left, info = f)
    o <- withr::with_seed(1, sample(nrow(ln[[f]])))
    expect_identical(.se_surface(train[[f]], grid[[f]], ln[[f]][o, ])$w, eg$left,
                     info = paste(f, "shuffled"))
  }
  # Two lines on a vertical, one above the other, meeting at the grid
  # point: the lower one, whatever the row order.
  v <- sf::st_sf(w = c(1, 2), geometry = sf::st_sfc(
    sf::st_linestring(rbind(c(4.5e5, 4.55e6 + 500), c(4.5e5, 4.55e6 + 1000))),
    sf::st_linestring(rbind(c(4.5e5, 4.55e6), c(4.5e5, 4.55e6 + 500))),
    crs = .se_crs_m()))
  pt <- sf::st_as_sf(data.frame(x = 4.5e5, y = 4.55e6 + 500), coords = c("x", "y"),
                     crs = .se_crs_m())
  vf <- .se_frames(v); pf <- .se_frames(pt)
  for (f in names(vf)) {
    expect_identical(.se_surface(train[[f]], pf[[f]], vf[[f]])$w, 2, info = f)
    expect_identical(.se_surface(train[[f]], pf[[f]], vf[[f]][2:1, ])$w, 2, info = f)
  }
})

test_that("predict_surface() settles a tie between polygons a grid point lies outside of", {
  # Points 100 m below the cells, on the lines of their shared edges: the
  # same distance from the two cells above.
  cells <- .se_frames(.se_cells())
  gi <- 1:9
  pts <- .se_frames(sf::st_as_sf(data.frame(x = 4.5e5 + 500 * gi, y = 4.55e6 - 100),
                                 coords = c("x", "y"), crs = .se_crs_m()))
  train <- .se_frames(.se_train())
  for (f in names(cells)) {
    expect_identical(.se_surface(train[[f]], pts[[f]], cells[[f]])$w, gi + 0, info = f)
    o <- withr::with_seed(2, sample(nrow(cells[[f]])))
    expect_identical(.se_surface(train[[f]], pts[[f]], cells[[f]][o, ])$w, gi + 0,
                     info = paste(f, "shuffled"))
  }
})

test_that("predict_surface() takes the nearest line or polygon as before where none tie", {
  cells <- .se_cells()
  train <- .se_train()
  # Grid points inside the cells, at their centres and off them.
  g <- expand.grid(i = 0:9, j = 0:7)
  inside <- sf::st_as_sf(data.frame(x = 4.5e5 + 500 * g$i + rep_len(c(250, 30, 470), 80),
                                    y = 4.55e6 + 500 * g$j + rep_len(c(250, 470, 30), 80)),
                         coords = c("x", "y"), crs = .se_crs_m())
  s <- .se_surface(train, inside, cells)
  expect_identical(s$w, cells$w[sf::st_nearest_feature(inside, cells)])
  expect_identical(s$w, seq_len(nrow(g)) + 0)
  # Random points against random discs that do not touch (a point inside
  # two overlapping discs is at distance 0 from both, a tie), and a point
  # 1 m nearer one of two cells than the other.
  ctr <- withr::with_seed(8, data.frame(x = stats::runif(60, 0, 5000),
                                        y = stats::runif(60, 0, 5000)))
  blobs <- sf::st_sf(w = withr::with_seed(9, stats::rnorm(60)),
                     geometry = sf::st_buffer(sf::st_as_sfc(sf::st_as_sf(
                       ctr, coords = c("x", "y"), crs = .se_crs_m())),
                       withr::with_seed(10, stats::runif(60, 20, 300))))
  blobs <- blobs[lengths(sf::st_intersects(blobs)) == 1L, ]
  gr <- withr::with_seed(11, sf::st_as_sf(
    data.frame(x = stats::runif(500, -500, 5500), y = stats::runif(500, -500, 5500)),
    coords = c("x", "y"), crs = .se_crs_m()))
  pts <- withr::with_seed(12, sf::st_as_sf(
    data.frame(x = stats::runif(200, 0, 5000), y = stats::runif(200, 0, 5000),
               z = stats::rnorm(200), w = stats::rnorm(200)),
    coords = c("x", "y"), crs = .se_crs_m()))
  s2 <- .se_surface(pts, gr, blobs)
  expect_identical(s2$w, blobs$w[sf::st_nearest_feature(gr, blobs)])
  near <- sf::st_as_sf(data.frame(x = 4.5e5 + 501, y = 4.55e6 - 100),
                       coords = c("x", "y"), crs = .se_crs_m())
  expect_identical(.se_surface(train, near, cells)$w, 2)
})

test_that("predict_surface() leaves an empty grid point and a layer in degrees as they were", {
  cells <- .se_cells()
  eg <- .se_edge_grid()
  train <- .se_train()
  # An empty grid point keeps NA; an empty covariate row is never taken.
  ge <- eg$pts
  sf::st_geometry(ge)[[3]] <- sf::st_point()
  ce <- cells
  sf::st_geometry(ce)[[2]] <- sf::st_polygon()
  s <- .se_surface(train, ge, ce)
  expect_true(is.na(s$w[3]))
  want <- eg$left
  want[want == 2] <- 3       # the empty cell's edge point takes its other side
  expect_identical(s$w[-3], want[-3])
  # Cells in degrees keep sf::st_nearest_feature()'s pick.
  cd <- sf::st_transform(cells, 4326)
  gd <- sf::st_transform(eg$pts, 4326)
  td <- sf::st_transform(train, 4326)
  fit <- lm_spatial_fit(td, "z", "w")
  sd <- .se_quiet(predict_surface(fit, grid = gd, covariates = cd[, "w"]))
  fcrs <- sf::st_crs(fit$data_sf)
  expect_identical(sd$w, cd$w[sf::st_nearest_feature(
    sf::st_transform(gd, fcrs), sf::st_transform(cd, fcrs))])
})

test_that(".surface_nearest_shapes() ranks features by box-centre column, then y, then row", {
  # Two coincident squares: the lower row wins.  A square and a line on its
  # left edge: the line's centre is in the lower column.
  sq <- sf::st_sfc(.se_square(0, 0, 10), .se_square(0, 0, 10),
                   sf::st_linestring(rbind(c(0, 0), c(0, 10))), crs = .se_crs_m())
  p <- sf::st_sfc(sf::st_point(c(0, 5)), sf::st_point(c(5, 5)), sf::st_point(c(10, 5)),
                  crs = .se_crs_m())
  expect_identical(.surface_nearest_shapes(p, sq), c(3L, 1L, 1L))
  expect_identical(.surface_nearest_shapes(p, sq[c(2, 1, 3)]), c(3L, 1L, 1L))
  expect_identical(.surface_nearest_shapes(p, sq[c(3, 2, 1)]), c(1L, 2L, 2L))
})

test_that("predict_surface() settles a tie in an oblique direction", {
  # Squares of 100 m spaced 500 m apart in a row, and grid points 250 m
  # below the gaps, each the same distance from a corner of the square on
  # either side, along a slant.  The candidates' square around a point holds
  # both corners; each point takes the square to its left.
  i  <- 0:9
  sq <- .se_frames(sf::st_sf(w = i + 1, geometry = sf::st_sfc(lapply(i, function(k)
    .se_square(4.5e5 + 500 * k, 4.55e6, 100)), crs = .se_crs_m())))
  pts <- .se_frames(sf::st_as_sf(data.frame(x = 4.5e5 + 500 * i[-10] + 300, y = 4.55e6 - 250),
                                 coords = c("x", "y"), crs = .se_crs_m()))
  for (f in names(sq)) {
    g <- sf::st_geometry(sq[[f]])
    p <- sf::st_geometry(pts[[f]])
    expect_identical(.surface_nearest_shapes(p, g), 1:9, info = f)
    for (s in 1:3) {
      o <- withr::with_seed(s, sample(10))
      expect_identical(o[.surface_nearest_shapes(p, g[o])], 1:9, info = paste(f, s))
    }
  }
})

test_that("predict_surface() gives a grid point without a location NA covariates whatever the layer", {
  # A POINT grid kept as given, one point at x = 1e200: lines, polygons and
  # a layer in degrees gave it the covariates of a row GEOS gave it, and a
  # finite .pred, where ?predict_surface states NA.
  cells <- .se_cells()
  lines <- sf::st_sf(w = cells$w,
                     geometry = sf::st_cast(sf::st_geometry(cells), "LINESTRING"))
  cpts  <- suppressWarnings(sf::st_centroid(cells))
  train <- .se_train()
  grid  <- sf::st_as_sf(data.frame(x = c(4.5e5 + 100, 1e200, 4.5e5 + 300),
                                   y = 4.55e6 + 100), coords = c("x", "y"),
                        crs = .se_crs_m())
  for (cv in list(cells, lines, cpts)) {
    s <- .se_surface(train, grid, cv)
    expect_identical(s$w, c(1, NA, 1))
    expect_true(is.na(s$.pred[2]))
    expect_true(all(is.finite(s$.pred[c(1, 3)])))
  }
  td <- sf::st_transform(train, 4326)
  xy <- sf::st_coordinates(sf::st_transform(grid[c(1, 3), ], 4326))
  gd <- sf::st_as_sf(data.frame(x = c(xy[1, 1], 1e200, xy[2, 1]),
                                y = c(xy[1, 2], xy[1, 2], xy[2, 2])),
                     coords = c("x", "y"), crs = 4326)
  for (cv in list(cells, cpts)) {
    s <- .se_surface(td, gd, sf::st_transform(cv, 4326))
    expect_identical(s$w, c(1, NA, 1))
    expect_true(is.na(s$.pred[2]))
  }
  # And at the level of the search.
  pu <- sf::st_sfc(sf::st_point(c(1e200, 4.55e6)), sf::st_point(c(4.5e5 + 1, 4.55e6 + 1)),
                   crs = .se_crs_m())
  expect_identical(.surface_nearest_shapes(pu, sf::st_geometry(cells)), c(NA, 1L))
})

test_that("predict_surface() reads the y of the box centres by the column rule", {
  # A polygon with a hole and the island that fills it have one box centre,
  # which after a move of the origin differs in its last bits in y only.  A
  # grid point on the hole's ring is at distance 0 from both and takes the
  # first of the two in row order in every frame.  A third cell 3 km away
  # widens the layer; the BUG-212 test below takes the donut and island
  # alone.
  ring  <- function(x0, y0, s) rbind(c(x0, y0), c(x0 + s, y0), c(x0 + s, y0 + s),
                                     c(x0, y0 + s), c(x0, y0))
  x0 <- 4.5e5; y0 <- 4.55e6
  donut  <- sf::st_polygon(list(ring(x0, y0, 1000), ring(x0 + 250, y0 + 250, 500)[5:1, ]))
  island <- .se_square(x0 + 250, y0 + 250, 500)
  far    <- .se_square(x0 + 3000, y0 + 2000, 500)
  di <- .se_frames(sf::st_sf(w = c(1, 2, 3), geometry = sf::st_sfc(donut, island, far,
                                                                   crs = .se_crs_m())))
  pd <- .se_frames(sf::st_as_sf(data.frame(x = x0 + c(250, 500, 250), y = y0 + c(500, 750, 250)),
                                coords = c("x", "y"), crs = .se_crs_m()))
  train <- .se_frames(.se_train())
  for (f in names(di)) {
    g <- sf::st_geometry(di[[f]])
    expect_identical(.surface_nearest_shapes(sf::st_geometry(pd[[f]]), g), c(1L, 1L, 1L),
                     info = f)
    expect_identical(.surface_nearest_shapes(sf::st_geometry(pd[[f]]), g[c(2, 1, 3)]),
                     c(1L, 1L, 1L), info = f)
    expect_identical(.surface_nearest_shapes(sf::st_geometry(pd[[f]]), g[c(3, 2, 1)]),
                     c(2L, 2L, 2L), info = f)
    expect_identical(.se_surface(train[[f]], pd[[f]], di[[f]])$w, c(1, 1, 1), info = f)
  }
})

test_that("predict_surface() reads a tie far from the layer by the relative reach", {
  # A grid point 100 km below two lines that meet above it.  The left line
  # is 5 mm farther than the right one: one part in 2e7, within one part in
  # ten million and beyond the absolute margin (about 0.5 mm here), so a
  # tie, which the left line takes.  At 20 mm (two parts in ten million) the
  # right line is the nearest.
  x0 <- 4.5e5; y0 <- 4.55e6
  two <- function(dy) sf::st_sfc(
    sf::st_linestring(rbind(c(x0 - 1000, y0 + dy), c(x0, y0 + dy))),
    sf::st_linestring(rbind(c(x0, y0), c(x0 + 1000, y0))), crs = .se_crs_m())
  p <- sf::st_sfc(sf::st_point(c(x0, y0 - 1e5)), crs = .se_crs_m())
  expect_identical(.surface_nearest_shapes(p, two(0.005)), 1L)
  expect_identical(.surface_nearest_shapes(p, two(0.005)[2:1]), 2L)
  expect_identical(.surface_nearest_shapes(p, two(0.02)), 2L)
  expect_identical(.surface_nearest_shapes(p, two(0.02)[2:1]), 1L)
})

test_that("predict_surface() sends a layer in degrees to sf::st_nearest_feature()", {
  local_mocked_bindings(.surface_nearest_shapes = function(...) stop("shapes rule called"))
  cd <- sf::st_transform(.se_cells(), 4326)
  td <- sf::st_transform(.se_train(), 4326)
  gd <- sf::st_transform(.se_edge_grid()$pts, 4326)
  s  <- .se_surface(td, gd, cd)
  expect_identical(nrow(s), 72L)
  expect_false(anyNA(s$w))
  # A projected layer of polygons goes to the shapes rule.
  expect_error(.se_surface(.se_train(), .se_edge_grid()$pts, .se_cells()),
               "shapes rule called")
})

test_that("predict_surface() measures few pairs for grid points far from the layer", {
  # 900 cells of 100 m and 40 grid points 10 to 20 km away: the square
  # around each point holds every cell, and the box distance keeps those
  # that can tie before any distance to a feature is taken.
  b <- expand.grid(ci = 0:29, cj = 0:29)
  cells <- sf::st_sfc(lapply(seq_len(nrow(b)), function(k)
    .se_square(4.5e5 + 100 * b$ci[k], 4.55e6 + 100 * b$cj[k], 100)), crs = .se_crs_m())
  a <- seq(0, 2 * pi, length.out = 41)[-41]
  rad <- 1e4 + 250 * seq_along(a)
  pts <- sf::st_sfc(lapply(seq_along(a), function(k) sf::st_point(c(
    4.5e5 + 1500 + rad[k] * cos(a[k]), 4.55e6 + 1500 + rad[k] * sin(a[k])))),
    crs = .se_crs_m())
  want <- .surface_nearest_shapes(pts, cells)
  # Each pick is within the reach of the least distance to any cell.
  dall <- matrix(as.numeric(sf::st_distance(pts, cells)), length(pts))
  dmin <- apply(dall, 1L, min)
  dgot <- dall[cbind(seq_along(pts), want)]
  expect_true(all(dgot <= dmin * (1 + .TIE_TOL) + .edge_tol(sf::st_bbox(cells))))
  measured <- 0
  real <- sf::st_distance
  local_mocked_bindings(st_distance = function(x, y, ...) {
    measured <<- measured + max(length(x), length(y))
    real(x, y, ...)
  }, .package = "sf")
  got <- .surface_nearest_shapes(pts, cells)
  expect_identical(got, want)
  # One distance per point for the nearest, and a few per point that ties.
  expect_lt(measured, 10 * length(pts))
})

test_that("predict_surface() settles a tie the same way in every frame where all box centres coincide", {
  # BUG-212: a layer of only a polygon with a hole and the island that fills
  # it.  The two have one box centre, which after a move of the origin
  # differs in its last bits in y, and the centres' own bounding box is
  # last-bits small.  Read against that box, the bits parted the centres and
  # a grid point on the hole's ring took the island with the origin moved
  # and the polygon in metres and US feet.  Read against the layer's
  # diagonal, the centres are one column and one row, and the tie goes to
  # the first in row order in every frame.
  ring  <- function(x0, y0, s) rbind(c(x0, y0), c(x0 + s, y0), c(x0 + s, y0 + s),
                                     c(x0, y0 + s), c(x0, y0))
  x0 <- 4.5e5; y0 <- 4.55e6
  donut  <- sf::st_polygon(list(ring(x0, y0, 1000), ring(x0 + 250, y0 + 250, 500)[5:1, ]))
  island <- .se_square(x0 + 250, y0 + 250, 500)
  di <- .se_frames(sf::st_sf(w = c(1, 2), geometry = sf::st_sfc(donut, island,
                                                                crs = .se_crs_m())))
  # Left edge, top edge, bottom-left corner, right edge, bottom edge.
  pd <- .se_frames(sf::st_as_sf(data.frame(x = x0 + c(250, 500, 250, 750, 600),
                                           y = y0 + c(500, 750, 250, 400, 250)),
                                coords = c("x", "y"), crs = .se_crs_m()))
  train <- .se_frames(.se_train())
  s <- r <- list()
  for (f in names(di)) {
    g <- sf::st_geometry(di[[f]])
    p <- sf::st_geometry(pd[[f]])
    expect_identical(.surface_nearest_shapes(p, g), rep(1L, 5L), info = f)
    expect_identical(.surface_nearest_shapes(p, g[2:1]), rep(1L, 5L), info = f)
    s[[f]] <- .se_surface(train[[f]], pd[[f]], di[[f]])
    r[[f]] <- .se_surface(train[[f]], pd[[f]], di[[f]][2:1, ])
    expect_identical(s[[f]]$w, rep(1, 5L), info = f)
    expect_identical(r[[f]]$w, rep(2, 5L), info = f)
  }
  expect_identical(s$usft$.pred, s$m$.pred)
  expect_identical(s$m01$.pred, s$m$.pred)
  expect_identical(r$usft$.pred, r$m$.pred)
  expect_identical(r$m01$.pred, r$m$.pred)
})

test_that("predict_surface() settles a tie the same way in every frame where coincident box centres part in x", {
  # BUG-212: two lines crossing at their midpoints, and a ring of 32 sides
  # with the disc that fills its hole.  Each pair has one box centre, which
  # in US feet or with the origin moved differs in its last bits in x.  Read
  # against the centres' own box, those bits made two columns and the grid
  # point took the second feature in those frames.
  x0 <- 4.5e5; y0 <- 4.55e6
  circ <- function(r) {
    a <- 2 * pi * (0:32) / 32
    a[33] <- 0
    cbind(x0 + 500 + r * cos(a), y0 + 500 + r * sin(a))
  }
  h <- sf::st_linestring(rbind(c(x0, y0 + 500), c(x0 + 1000, y0 + 500)))
  v <- sf::st_linestring(rbind(c(x0 + 500, y0), c(x0 + 500, y0 + 1000)))
  ann  <- sf::st_polygon(list(circ(500), circ(250)[33:1, ]))
  disc <- sf::st_polygon(list(circ(250)))
  # The crossing, and eight vertices of the hole's ring.
  cases <- list(lines = list(g = sf::st_sfc(h, v, crs = .se_crs_m()),
                             p = cbind(x0 + 500, y0 + 500)),
                ring  = list(g = sf::st_sfc(ann, disc, crs = .se_crs_m()),
                             p = circ(250)[seq(1, 29, by = 4), ]))
  train <- .se_frames(.se_train())
  for (k in names(cases)) {
    n <- nrow(cases[[k]]$p)
    L <- .se_frames(sf::st_sf(w = c(1, 2), geometry = cases[[k]]$g))
    P <- .se_frames(sf::st_as_sf(data.frame(x = cases[[k]]$p[, 1], y = cases[[k]]$p[, 2]),
                                 coords = c("x", "y"), crs = .se_crs_m()))
    s <- r <- list()
    for (f in names(L)) {
      g <- sf::st_geometry(L[[f]])
      p <- sf::st_geometry(P[[f]])
      expect_identical(.surface_nearest_shapes(p, g), rep(1L, n), info = paste(k, f))
      expect_identical(.surface_nearest_shapes(p, g[2:1]), rep(1L, n), info = paste(k, f))
      s[[f]] <- .se_surface(train[[f]], P[[f]], L[[f]])
      r[[f]] <- .se_surface(train[[f]], P[[f]], L[[f]][2:1, ])
      expect_identical(s[[f]]$w, rep(1, n), info = paste(k, f))
      expect_identical(r[[f]]$w, rep(2, n), info = paste(k, f))
    }
    expect_identical(s$usft$.pred, s$m$.pred, info = k)
    expect_identical(s$m01$.pred, s$m$.pred, info = k)
  }
})

test_that(".near_columns() reads its step against the extent it is given", {
  # Two first coordinates 1.2e-10 apart, and the points' own box no larger:
  # without an extent the step is more than one part in ten million of the
  # box's diagonal, so two columns, as before.  Against the diagonal of a
  # 1000 m box they are one column.
  x <- c(1e6, 1e6 + 1.2e-10); y <- c(5e6, 5e6)
  expect_identical(.near_columns(x, y), 1:2)
  expect_identical(.near_columns(x, y, sqrt(2) * 1000), c(1L, 1L))
  expect_identical(.near_columns(c(1e6, 1e6), y), c(1L, 1L))
  # The default is the diagonal of the points' own bounding box.
  xs <- withr::with_seed(3, c(stats::runif(40, 0, 1000), 500 + c(0, 1e-5, 2e-4)))
  ys <- withr::with_seed(4, stats::runif(43, 0, 800))
  ex <- sqrt(diff(range(xs))^2 + diff(range(ys))^2)
  expect_identical(.near_columns(xs, ys), .near_columns(xs, ys, ex))
  expect_identical(.near_columns(xs, ys)[41:43], c(1L, 1L, 2L) +
                     .near_columns(xs, ys)[41] - 1L)
  # A larger extent joins steps the points' own box keeps apart.
  expect_identical(length(unique(.near_columns(xs, ys, 1e4)[41:43])), 1L)
})

test_that("?predict_surface states the tie rule for lines and polygons", {
  file <- test_path("..", "..", "man", "predict_surface.Rd")
  rd <- if (file.exists(file)) tools::parse_Rd(file)
        else tools::Rd_db("spatialkit")[["predict_surface.Rd"]]
  page <- gsub("[[:space:]]+", " ", paste(unlist(rd), collapse = ""))
  expect_match(page, paste0(
    "A grid point the same distance from several features of a layer of lines or ",
    "polygons, or of one that mixes them with points"), fixed = TRUE)
  expect_match(page, paste0(
    "takes the one whose bounding-box centre is lowest in x, then lowest in y, x values ",
    "and y values within one part in ten million of the diagonal of the bounding box ",
    "of covariates counting as one, then the first in row order. For a layer in degrees ",
    "the nearest feature is sf::st_nearest_feature()'s"), fixed = TRUE)
  expect_no_match(page, "For lines and polygons, and for points in degrees", fixed = TRUE)
})
