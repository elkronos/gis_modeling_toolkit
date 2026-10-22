# tests/testthat/test-fix-crs-polygon-order.R
# ---------------------------------------------------------------------------
# The CRS chosen for a polygon layer, and the area distortion measured on
# one, do not depend on the order of its features.  .crs_distance_error()
# thins the vertices of a few detailed outlines, and .crs_area_error()
# measures at most 200 of a layer's polygons; both samples were taken at even
# spacing along the rows.  The vertices are now taken along the distinct
# vertices sorted by coordinate, as the points already were, and the polygons
# along the polygons sorted by their bounding boxes and then by the means of
# their vertices, with the polygons that reach furthest along x and y added.
# ---------------------------------------------------------------------------

quiet <- function(expr) suppressWarnings(suppressMessages(expr))
square <- function(x0, y0, h)
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + h, y0), c(x0 + h, y0 + h), c(x0, y0 + h), c(x0, y0))))
# Seven row orders besides the layer's own.
row_orders <- function(n) { set.seed(99); lapply(1:7, function(s) sample(n)) }

test_that("a polygon layer of 300 features in seven row orders gives one CRS and one area error", {
  # 300 squares of 0.25 degrees over the conterminous United States.  The CRS
  # was already the same in every order (one point per feature, sorted); the
  # area error of the projected layer measured 200 of the 300 by row and
  # gave seven figures in eight orders.
  set.seed(7)
  ll <- sf::st_sf(id = 1:300, geometry = sf::st_sfc(lapply(1:300, function(i)
    square(runif(1, -122, -70), runif(1, 26, 48), 0.25)), crs = 4326))
  d0 <- quiet(ensure_projected(ll))
  a0 <- quiet(ensure_projected(ll, purpose = "area"))
  e0 <- .crs_area_error(a0)
  expect_true(is.finite(e0))
  for (p in row_orders(300)) {
    d <- quiet(ensure_projected(ll[p, ]))
    a <- quiet(ensure_projected(ll[p, ], purpose = "area"))
    expect_true(sf::st_crs(d) == sf::st_crs(d0))
    expect_true(sf::st_crs(a) == sf::st_crs(a0))
    expect_identical(attr(a, "crs_choice")$distance_error, attr(a0, "crs_choice")$distance_error)
    expect_identical(.crs_area_error(a), e0)
    expect_identical(.crs_area_error(d), .crs_area_error(d0))
  }
})

test_that("cells near the area tolerance are refused, or given a density, in every row order alike", {
  # 299 cells of 20 km within 200 km of the central meridian of UTM zone 14
  # and one 700 km east of it: 1.2 percent over all 300, 0.15 percent
  # without the far one.  By row, the far cell was among the 200 measured in
  # some orders only, and summarize_by_cell(area = TRUE) refused the cells in
  # those orders and gave a density in the others.
  set.seed(2)
  xs <- c(500000 + runif(299, -2e5, 2e5), 500000 + 7e5); ys <- c(runif(299, 3.6e6, 4.6e6), 4.0e6)
  cells <- sf::st_sf(cell_id = 1:300, geometry = sf::st_sfc(lapply(1:300, function(i)
    square(xs[i], ys[i], 2e4)), crs = 32614))
  set.seed(3)
  pts <- do.call(rbind, lapply(1:300, function(i) data.frame(cell_id = i,
    x = xs[i] + runif(3, 0, 2e4), y = ys[i] + runif(3, 0, 2e4), z = rnorm(3))))
  pts <- sf::st_as_sf(pts, coords = c("x", "y"), crs = 32614)
  e0 <- .crs_area_error(cells)
  expect_gt(e0, .area_error_tol)
  for (p in row_orders(300)) {
    expect_identical(.crs_area_error(cells[p, ]), e0)
    expect_error(quiet(summarize_by_cell(pts, "z", id_col = "cell_id", cells_sf = cells[p, ],
                                         area = TRUE, estimate = "mean")),
                 "`area = TRUE` refused")
    lines <- capture_spatialkit_log(quiet(ensure_projected(cells[p, ], purpose = "area")))
    expect_true(log_has(lines, "distorts areas across this extent by up to"))
  }
  # The same cells with the far one moved in: a density in every order, and
  # the same area_error attribute.
  near <- cells
  sf::st_geometry(near)[300] <- sf::st_sfc(square(500000 + 1e5, 4.0e6, 2e4), crs = 32614)
  s0 <- quiet(summarize_by_cell(pts, "z", id_col = "cell_id", cells_sf = near, area = TRUE,
                                estimate = "mean"))
  expect_lt(attr(s0, "area_error"), .area_error_tol)
  for (p in row_orders(300)) {
    s <- quiet(summarize_by_cell(pts, "z", id_col = "cell_id", cells_sf = near[p, ], area = TRUE,
                                 estimate = "mean"))
    expect_identical(attr(s, "area_error"), attr(s0, "area_error"))
  }
})

test_that("a few detailed outlines get one CRS, and one score per candidate, in any feature order", {
  # Five outlines of 200 vertices each over the conterminous United States:
  # fewer than 40 features, so their 1,005 vertices are scored too, thinned
  # to 160 first.  Thinned by row, the eight orders below gave an Albers CRS
  # three times and a Lambert azimuthal one five times, and the scores of
  # every candidate moved.
  blob <- function(lon0, lat0, r, nv, wig) {
    a  <- seq(0, 2 * pi, length.out = nv + 1L)[-(nv + 1L)]
    rr <- r * (1 + wig * sin(5 * a + lon0) + 0.5 * wig * cos(11 * a + lat0))
    m  <- cbind(lon0 + rr * cos(a) / cos(lat0 * pi / 180), lat0 + rr * sin(a))
    sf::st_polygon(list(rbind(m, m[1, ])))
  }
  set.seed(49)
  lon <- runif(5, -118, -72); lat <- runif(5, 27, 48); r <- runif(5, 1, 3)
  outl <- sf::st_sf(id = 1:5, geometry = sf::st_sfc(lapply(1:5, function(i)
    blob(lon[i], lat[i], r[i], 200, 0.15)), crs = 4326))
  expect_gt(nrow(sf::st_coordinates(outl)), 160)
  x0 <- quiet(ensure_projected(outl))
  orders <- list(c(1, 4, 2, 3, 5), c(3, 2, 4, 5, 1), 5:1, c(2, 4, 5, 3, 1), c(5, 4, 2, 3, 1),
                 c(1, 3, 2, 5, 4), c(4, 1, 5, 2, 3))
  for (p in orders) {
    x <- quiet(ensure_projected(outl[p, ]))
    expect_true(sf::st_crs(x) == sf::st_crs(x0), info = paste(p, collapse = ""))
    expect_identical(attr(x, "crs_choice")$distance_error, attr(x0, "crs_choice")$distance_error,
                     info = paste(p, collapse = ""))
  }
  # A ring listed from another vertex is the same outline, and is scored the
  # same.  Its closing vertex repeats another vertex, and counted twice it
  # moved the vertices kept: started 50 vertices on, the rings got a Lambert
  # azimuthal CRS where the outlines as built got an Albers one.
  rot <- outl
  sf::st_geometry(rot) <- sf::st_sfc(lapply(sf::st_geometry(outl), function(g) {
    m <- g[[1]][-nrow(g[[1]]), ]; m <- m[c(51:nrow(m), 1:50), ]
    sf::st_polygon(list(rbind(m, m[1, ])))
  }), crs = 4326)
  xr <- quiet(ensure_projected(rot))
  expect_true(sf::st_crs(xr) == sf::st_crs(x0))
  expect_identical(attr(xr, "crs_choice")$distance_error, attr(x0, "crs_choice")$distance_error)
  for (code in c(32613, 5070))
    expect_identical(.crs_distance_error(rot, sf::st_crs(code)),
                     .crs_distance_error(outl, sf::st_crs(code)))
})

test_that("a thinned outline is not densified, so a border two features share scores the same in either order", {
  # Two boxes that share a side, 17 copies of each: 34 features and 170
  # vertices, so the outlines are thinned, to their 6 distinct vertices.  Of
  # a vertex on the shared side, the copy kept is the one of whichever box
  # comes first.  Were the vertices kept densified as if they were
  # neighbours along a ring, the points added would follow that copy: in Web
  # Mercator the error would be 0.417 with the left box first and 0.414 with
  # the right one.
  box <- function(x0, y0, w, h)
    sf::st_polygon(list(rbind(c(x0, y0), c(x0 + w, y0), c(x0 + w, y0 + h), c(x0, y0 + h), c(x0, y0))))
  a <- box(-110, 30, 12, 15); b <- box(-98, 30, 14, 15)
  lay <- function(gs) sf::st_sf(id = seq_along(gs), geometry = sf::st_sfc(gs, crs = 4326))
  ab <- lay(c(rep(list(a), 17), rep(list(b), 17)))
  ba <- lay(c(rep(list(b), 17), rep(list(a), 17)))
  expect_gt(nrow(sf::st_coordinates(ab)), 160)
  for (code in c(3857, 32613, 5070)) {
    e <- .crs_distance_error(ab, sf::st_crs(code))
    expect_true(is.finite(e))
    expect_identical(.crs_distance_error(ba, sf::st_crs(code)), e, info = paste("EPSG", code))
  }
})

test_that("the help and NEWS say how the area check and the outline thinning sample", {
  root <- testthat::test_path("..", "..")
  if (!(file.exists(file.path(root, "DESCRIPTION")) && dir.exists(file.path(root, "man"))))
    skip("the help's source is not here to read")
  read <- function(...) gsub("[[:space:]]+", " ",
                             paste(readLines(file.path(root, ...), warn = FALSE), collapse = " "))
  rd <- read("man", "summarize_by_cell.Rd")
  expect_match(rd, paste0("spread of planar-to-geodesic area ratios over the cells (over 200 of ",
                          "them when there are more, taken at even spacing along the cells sorted ",
                          "by their bounding boxes in that CRS and then by the means of their ",
                          "vertices, together with the cells that reach furthest along x and y ",
                          "either way, so the order of the rows does not matter)"), fixed = TRUE)
  news <- read("NEWS.md")
  expect_match(news, paste0("The projection chosen for a few detailed outlines, and the area check ",
                            "on more than 200 polygons, depended on the order of the features."),
               fixed = TRUE)
  expect_match(news, paste0("The vertices are now thinned along the outlines' distinct vertices ",
                            "sorted by coordinate"), fixed = TRUE)
  expect_match(news, paste0("The polygons are taken along the polygons sorted by their bounding ",
                            "boxes, and polygons that share a box, as the Delaunay triangles of ",
                            "gridded points do, by the means of their vertices; the polygons that ",
                            "reach furthest along x and y, either way, are always measured."),
               fixed = TRUE)
})

test_that("cells that share a bounding box get one area error, and one decision, in any row order", {
  # The Delaunay triangles of a 16 x 11 lattice of points: 300 triangles,
  # each sharing its bounding box with the other half of its square.  Sorted
  # on the box alone they kept their row order among themselves: in the 21
  # orders below the area error took 12 values, and summarize_by_cell(area =
  # TRUE) refused the cells in two orders and gave a density in the others.
  lat <- expand.grid(x = seq(283000, 283000 + 15 * 6e4, by = 6e4),
                     y = seq(3.6e6, 3.6e6 + 10 * 6e4, by = 6e4))
  set.seed(7); lat$z <- rnorm(nrow(lat))
  lat <- sf::st_as_sf(lat, coords = c("x", "y"), crs = 32614)
  cells <- quiet(build_tessellation(lat, method = "triangles", quiet = TRUE))$cells
  expect_equal(nrow(cells), 300L)
  bb <- t(vapply(sf::st_geometry(cells), function(s) as.numeric(sf::st_bbox(s)), numeric(4)))
  expect_true(all(duplicated(bb) | duplicated(bb, fromLast = TRUE)))
  asg <- quiet(assign_features_to_polygons(lat, cells, polygon_id_col = "cell_id", quiet = TRUE))
  outcome <- function(o) tryCatch(
    attr(quiet(summarize_by_cell(asg, "z", id_col = "cell_id", cells_sf = cells[o, ],
                                 area = TRUE, estimate = "mean")), "area_error"),
    error = function(e) if (grepl("`area = TRUE` refused", conditionMessage(e))) "refused" else
      conditionMessage(e))
  set.seed(11); ords <- c(list(seq_len(300)), lapply(1:20, function(s) sample(300)))
  e0 <- .crs_area_error(cells)
  r0 <- outcome(ords[[1L]])
  expect_identical(r0, e0)
  for (o in ords[-1L]) {
    expect_identical(.crs_area_error(cells[o, ]), e0)
    expect_identical(outcome(o), r0)
  }
})

test_that("square, hexagonal and split-square grids of more than 200 cells get one area error in any row order", {
  # The cells of one column of a square grid share their xmin, and the two
  # triangles of a split square, or a square and the diamond inscribed in
  # it, share their whole bounding box.  A sort that stops short of a key
  # telling them apart leaves them in their row order, and with them which
  # of them are measured: on the box alone the split squares and the
  # diamonds below gave several figures in twelve orders, on xmin alone
  # every layer here did.
  grid <- function(x0, y0, nx, ny, w, crs, square = TRUE) {
    g <- sf::st_make_grid(sf::st_as_sfc(sf::st_bbox(c(xmin = x0, ymin = y0, xmax = x0 + nx * w,
                                                       ymax = y0 + ny * w), crs = sf::st_crs(crs))),
                          n = c(nx, ny), square = square)
    sf::st_sf(id = seq_along(g), geometry = g)
  }
  split <- function(x0, y0, nx, ny, w, crs) {
    out <- list()
    for (i in 0:(nx - 1)) for (j in 0:(ny - 1)) {
      a <- c(x0 + i * w, y0 + j * w); b <- a + c(w, 0); cc <- a + c(w, w); d <- a + c(0, w)
      out[[length(out) + 1L]] <- sf::st_polygon(list(rbind(a, b, cc, a)))
      out[[length(out) + 1L]] <- sf::st_polygon(list(rbind(a, cc, d, a)))
    }
    sf::st_sf(id = seq_along(out), geometry = sf::st_sfc(out, crs = crs))
  }
  diamond <- function(x0, y0, nx, ny, w, crs) {
    out <- list()
    for (i in 0:(nx - 1)) for (j in 0:(ny - 1)) {
      a <- c(x0 + i * w * 1.5, y0 + j * w * 1.5)
      out[[length(out) + 1L]] <- square(a[1], a[2], w)
      out[[length(out) + 1L]] <- sf::st_polygon(list(rbind(a + c(w / 2, 0), a + c(w, w / 2),
                                                           a + c(w / 2, w), a + c(0, w / 2),
                                                           a + c(w / 2, 0))))
    }
    sf::st_sf(id = seq_along(out), geometry = sf::st_sfc(out, crs = crs))
  }
  layers <- list(
    "square grid, EPSG 3857"     = grid(-1.1e7, 4.0e6, 20, 15, 5e4, 3857),
    "square grid, UTM 14"        = grid(2e5, 3.6e6, 20, 15, 3e4, 32614),
    "hexagons, EPSG 3857"        = grid(-1.1e7, 4.0e6, 20, 16, 5e4, 3857, square = FALSE),
    "split squares, EPSG 3857"   = split(-1.1e7, 4.0e6, 15, 10, 5e4, 3857),
    "split squares, UTM 14"      = split(3e5, 3.6e6, 15, 10, 4e4, 32614),
    "square and diamond, EPSG 3857" = diamond(-1.1e7, 4.0e6, 15, 10, 5e4, 3857))
  set.seed(99); ords <- lapply(1:12, function(s) sample(1e4))
  for (nm in names(layers)) {
    lay <- layers[[nm]]
    expect_gt(nrow(lay), 200)
    e0 <- .crs_area_error(lay)
    expect_true(is.finite(e0), info = nm)
    for (o in ords)
      expect_identical(.crs_area_error(lay[order(o[seq_len(nrow(lay))]), ]), e0, info = nm)
  }
})

test_that("the cell that reaches furthest north is measured wherever it lies along x", {
  # 299 cells of 20 km between 36N and 42N in a conic CRS, whose area scale
  # changes with latitude, and one more at 53.5N: over all 300 cells the
  # ratios spread by 7.5 percent.  Taken at even spacing along the cells
  # sorted by their boxes, xmin first, the far cell was left out at five of
  # these longitudes in every row order; the 200 measured then spread by 0.3
  # percent, and summarize_by_cell(area = TRUE) gave a density.
  lcc <- sf::st_crs("+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39 +lon_0=-96 +datum=WGS84 +units=m +no_defs")
  set.seed(12)
  ll <- sf::st_as_sf(data.frame(x = runif(299, -104, -88), y = runif(299, 36, 42)),
                     coords = c("x", "y"), crs = 4326)
  xy <- sf::st_coordinates(sf::st_transform(ll, lcc))
  far_ll <- sf::st_sfc(lapply(-103:-89, function(l) sf::st_point(c(l, 53.5))), crs = 4326)
  far <- sf::st_coordinates(sf::st_transform(far_ll, lcc))
  layer <- function(k) sf::st_sf(cell_id = 1:300, geometry = sf::st_sfc(c(
    lapply(1:299, function(i) square(xy[i, 1], xy[i, 2], 2e4)),
    list(square(far[k, 1], far[k, 2], 2e4))), crs = lcc))
  for (k in seq_len(nrow(far))) {
    cells <- layer(k)
    expect_gt(.crs_area_error(cells, max_n = 1e6), 0.07)
    expect_gt(.crs_area_error(cells), 0.07, label = paste("far cell at", -104 + k))
  }
  # The far cell at 100W, in eight row orders: refused in every one.
  cells <- layer(4L)
  set.seed(3)
  pts <- do.call(rbind, lapply(1:300, function(i) {
    b <- sf::st_bbox(sf::st_geometry(cells)[i])
    data.frame(cell_id = i, x = b[["xmin"]] + runif(3, 0, 2e4), y = b[["ymin"]] + runif(3, 0, 2e4),
               z = rnorm(3))
  }))
  pts <- sf::st_as_sf(pts, coords = c("x", "y"), crs = lcc)
  for (p in c(list(seq_len(300)), row_orders(300)))
    expect_error(quiet(summarize_by_cell(pts, "z", id_col = "cell_id", cells_sf = cells[p, ],
                                         area = TRUE, estimate = "mean")),
                 "`area = TRUE` refused")
})

test_that("a thinned outline scores as the same outline does unthinned", {
  # One outline over the conterminous United States, of 2,000 vertices
  # (thinned to 160) and of 150 (scored whole).  The vertices kept have to
  # span the outline: kept at one end of it, the 160 westmost, the score in
  # UTM 14 was 1.9 percent against 2.6 for the outline of 150 vertices, and
  # in Web Mercator 28 percent against 51.
  blob <- function(nv) {
    a  <- seq(0, 2 * pi, length.out = nv + 1L)[-(nv + 1L)]
    rr <- 11 * (1 + 0.1 * sin(5 * a - 97) + 0.05 * cos(11 * a + 38))
    m  <- cbind(-97 + rr * cos(a) / cos(38 * pi / 180), 38 + rr * sin(a))
    sf::st_sf(id = 1, geometry = sf::st_sfc(sf::st_polygon(list(rbind(m, m[1, ]))), crs = 4326))
  }
  dense <- blob(2000); coarse <- blob(150)
  expect_gt(nrow(sf::st_coordinates(dense)), 160)
  expect_lt(nrow(sf::st_coordinates(coarse)), 160)
  for (code in c(32614, 5070, 3857)) {
    e_d <- .crs_distance_error(dense, sf::st_crs(code))
    e_c <- .crs_distance_error(coarse, sf::st_crs(code))
    expect_lt(abs(e_d / e_c - 1), 0.05, label = paste("EPSG", code))
  }
})
