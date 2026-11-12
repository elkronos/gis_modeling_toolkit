# assign_features_to_polygons() under the package's tie rule (BUG-107,
# BUG-120): a point on the edge two polygons share is in both in every unit
# and origin, the smallest-area rule compares bounding-box centres within the
# same margin and then the mean of the vertices, and a feature whose polygons
# are equal on every key is given none, with a warning, whatever the row
# order of `polygons_sf`.

.at_crs_m    <- function() sf::st_crs(26915)
.at_crs_usft <- function() sf::st_crs("+proj=utm +zone=15 +datum=NAD83 +units=us-ft +no_defs")
.at_crs_m01  <- function() sf::st_crs(paste0(
  "+proj=tmerc +lat_0=0 +lon_0=-93 +k=0.9996 +x_0=500000.1 +y_0=0.1 +datum=NAD83 ",
  "+units=m +no_defs"))
.at_frames <- function(pm) list(m = pm, usft = sf::st_transform(pm, .at_crs_usft()),
                                m01 = sf::st_transform(pm, .at_crs_m01()))
.at_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

# A lattice of nx by ny points at 250 m from (450000, 4550000); `r` is the row.
.at_lattice <- function(nx, ny, jitter = 0) {
  g <- expand.grid(i = 0:(nx - 1), j = 0:(ny - 1))
  x <- 4.5e5 + g$i * 250; y <- 4.55e6 + g$j * 250
  if (jitter > 0) {
    x <- x + withr::with_seed(99, stats::runif(length(x), -jitter, jitter))
    y <- y + withr::with_seed(98, stats::runif(length(y), -jitter, jitter))
  }
  sf::st_as_sf(data.frame(x = x, y = y, r = seq_along(x)), coords = c("x", "y"),
               crs = .at_crs_m())
}

# The cell of each row, named by the cell's centroid in metres, in the order
# of `r`.
.at_cells_of <- function(a, cells) {
  c <- suppressWarnings(sf::st_coordinates(sf::st_centroid(
    sf::st_transform(sf::st_geometry(cells), .at_crs_m()))))
  key <- paste(round(c[, 1], 3), round(c[, 2], 3))
  key[match(a$poly_id, cells$poly_id)][order(a$r)]
}

# The 10 x 7 grid of square cells over the lattice's box widened by 250 m:
# cells 950 m wide, every fifth column edge on a column of the lattice.
.at_box <- function() {
  bb <- sf::st_bbox(.at_lattice(37, 25))
  bb[1:2] <- bb[1:2] - 250; bb[3:4] <- bb[3:4] + 250
  bb
}
.at_grid <- function(crs) .at_quiet(create_grid_polygons(
  sf::st_as_sf(sf::st_transform(sf::st_as_sfc(.at_box()), crs)), n = c(10, 7),
  type = "square", quiet = TRUE))

test_that("a lattice point on a grid edge goes to one cell in every unit and origin", {
  skip_on_cran()
  for (jit in c(0, 60)) {
    fr <- .at_frames(.at_lattice(37, 25, jitter = jit))
    grid_m <- .at_grid(.at_crs_m())
    for (tb in c("smallest_area", "first")) {
      # (a) the grid built in each frame over the box moved there.
      built <- lapply(fr, function(l) {
        g <- .at_grid(sf::st_crs(l))
        .at_cells_of(.at_quiet(assign_features_to_polygons(l, g, tie_break = tb)), g)
      })
      # (b) the grid built once in metres and moved with the points.
      moved <- lapply(fr, function(l) {
        g <- if (isTRUE(sf::st_crs(l) == .at_crs_m())) grid_m
             else sf::st_transform(grid_m, sf::st_crs(l))
        .at_cells_of(.at_quiet(assign_features_to_polygons(l, g, tie_break = tb)), g)
      })
      info <- sprintf("jitter %g, %s", jit, tb)
      expect_false(anyNA(built$m), info = info)
      expect_length(built$m, 925L)
      for (f in c("usft", "m01")) {
        expect_identical(built[[f]], built$m, info = paste(info, f, "built"))
        expect_identical(moved[[f]], built$m, info = paste(info, f, "moved"))
      }
    }
  }
})

test_that("the lattice on the edges is assigned by the stated rule, and the control as the polygons hold it", {
  skip_on_cran()
  g <- .at_grid(.at_crs_m())
  lat <- .at_lattice(37, 25)
  a <- .at_quiet(assign_features_to_polygons(lat, g))
  # The cells are axis-aligned rectangles in metres: a point is in every cell
  # whose box holds it (to 1e-6 m), and the rule takes the smallest area,
  # then the lowest, then the leftmost centre.
  xy <- sf::st_coordinates(lat)
  bb <- t(vapply(sf::st_geometry(g), function(x) as.numeric(sf::st_bbox(x)), numeric(4)))
  ar <- signif(as.numeric(sf::st_area(g)), 9)
  cand <- lapply(seq_len(nrow(xy)), function(k)
    which(bb[, 1] - 1e-6 <= xy[k, 1] & xy[k, 1] <= bb[, 3] + 1e-6 &
            bb[, 2] - 1e-6 <= xy[k, 2] & xy[k, 2] <= bb[, 4] + 1e-6))
  want <- vapply(cand, function(cc)
    cc[order(ar[cc], (bb[cc, 2] + bb[cc, 4]) / 2, (bb[cc, 1] + bb[cc, 3]) / 2)[1L]], integer(1))
  expect_true(any(lengths(cand) > 1L))
  expect_identical(as.integer(a$poly_id[order(a$r)]), as.integer(g$poly_id[want]))
  # The rows on an edge are counted in "ties".
  expect_identical(attr(a, "ties")$which, which(lengths(cand) > 1L))
  # The jittered control: no point within the margin of an edge, so every
  # point is in the one cell sf::st_intersects() finds, as before.
  ctl <- .at_lattice(37, 25, jitter = 60)
  hits <- sf::st_intersects(ctl, g)
  expect_true(all(lengths(hits) == 1L))
  b <- .at_quiet(assign_features_to_polygons(ctl, g))
  expect_identical(as.integer(b$poly_id[order(b$r)]), as.integer(g$poly_id[unlist(hits)]))
  expect_identical(attr(b, "ties")$n, 0L)
})

# Two 100 m cells sharing the edge x = 500100 in EPSG:32617.
.at_two_cells <- function(x0 = 5e5, y0 = 4e6, crs = 32617) {
  P <- function(m) sf::st_polygon(list(m))
  sf::st_sf(pid = c("C1", "C2"), geometry = sf::st_sfc(
    P(rbind(c(x0, y0), c(x0 + 100, y0), c(x0 + 100, y0 + 100), c(x0, y0 + 100), c(x0, y0))),
    P(rbind(c(x0 + 100, y0), c(x0 + 200, y0), c(x0 + 200, y0 + 100), c(x0 + 100, y0 + 100),
            c(x0 + 100, y0))), crs = crs))
}

test_that("a point within the margin of a shared edge is in both cells, and one beyond it in one", {
  cells <- .at_two_cells()
  # The margin: one part in 1e10 of the largest coordinate (4e-4 m here), at
  # most a millionth of the diagonal of the cells' box (2.24e-4 m for the
  # 200 m x 100 m pair).  A point 2e-4 m east of the edge is within it; one
  # 3e-4 m east of it, which the margin of the largest coordinate alone
  # took in, lies in C2 alone (BUG-225).
  dx <- c(-1e-5, 0, 1e-5, 2e-4, 1e-3, -1e-3, 3e-4)
  pts <- sf::st_sf(r = seq_along(dx), geometry = sf::st_sfc(
    lapply(dx, function(d) sf::st_point(c(500100 + d, 4000050))), crs = 32617))
  for (ord in list(1:2, 2:1)) {
    a <- .at_quiet(assign_features_to_polygons(pts, cells[ord, ], polygon_id_col = "pid"))
    expect_identical(a$pid, c("C1", "C1", "C1", "C1", "C2", "C1", "C2"),
                     info = paste(ord, collapse = ""))
    expect_identical(attr(a, "ties")$which, 1:4)
    f <- .at_quiet(assign_features_to_polygons(pts, cells[ord, ], polygon_id_col = "pid",
                                               tie_break = "first"))
    expect_identical(f$pid, c(rep(cells$pid[ord][1], 4), "C2", "C1", "C2"))
  }
  # Another predicate is applied as given, with no margin.
  w <- .at_quiet(assign_features_to_polygons(pts, cells, polygon_id_col = "pid",
                                             predicate = sf::st_within,
                                             keep_unassigned = TRUE))
  expect_identical(w$pid, c("C1", NA, "C2", "C2", "C2", "C1", "C2"))
})

test_that("in degrees with s2 on a point a hair either side of a meridian edge is a tie", {
  skip_if_not(sf::sf_use_s2())
  P <- function(m) sf::st_polygon(list(m))
  cells <- sf::st_sf(pid = c("W", "E"), geometry = sf::st_sfc(
    P(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0))),
    P(rbind(c(1, 0), c(2, 0), c(2, 1), c(1, 1), c(1, 0))), crs = 4326))
  dx <- c(-1e-12, 0, 1e-12, 1e-6)
  pts <- sf::st_sf(r = 1:4, geometry = sf::st_sfc(
    lapply(dx, function(d) sf::st_point(c(1 + d, 0.5))), crs = 4326))
  for (ord in list(1:2, 2:1)) {
    a <- .at_quiet(assign_features_to_polygons(pts, cells[ord, ], polygon_id_col = "pid"))
    expect_identical(a$pid, c("W", "W", "W", "E"), info = paste(ord, collapse = ""))
  }
})

# The cases of BUG-120, EPSG:32617.
.at_P <- function(m) sf::st_polygon(list(m))
.at_tri <- function() sf::st_sf(pid = c("T1", "T2"), geometry = sf::st_sfc(
  .at_P(rbind(c(0, 0), c(100, 0), c(100, 100), c(0, 0))),
  .at_P(rbind(c(0, 0), c(100, 100), c(0, 100), c(0, 0))), crs = 32617))
.at_sq <- function() .at_P(rbind(c(0, 0), c(100, 0), c(100, 100), c(0, 100), c(0, 0)))
.at_rh <- function() .at_P(rbind(c(-50, 50), c(50, 0), c(150, 50), c(50, 100), c(-50, 50)))
.at_pts <- function() sf::st_sf(r = 1:3, geometry = sf::st_sfc(
  sf::st_point(c(50, 50)), sf::st_point(c(20, 20)), sf::st_point(c(80, 80)), crs = 32617))
.at_feat <- function() sf::st_sf(r = 1L, geometry = sf::st_sfc(
  .at_P(rbind(c(40, 40), c(60, 40), c(60, 60), c(40, 60), c(40, 40))), crs = 32617))

test_that("the two triangles of a split square give one triangle per feature in both row orders", {
  tri <- .at_tri()
  for (ord in list(1:2, 2:1)) {
    w <- character(0)
    a <- withCallingHandlers(
      assign_features_to_polygons(.at_pts(), tri[ord, ], polygon_id_col = "pid"),
      warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") })
    # The lower-right triangle's vertices average to (66.7, 33.3), the other's
    # to (33.3, 66.7): the lower mean wins.
    expect_identical(a$pid, c("T1", "T1", "T1"), info = paste(ord, collapse = ""))
    expect_identical(attr(a, "ties")$which, 1:3)
    expect_null(attr(a, "ties")$unresolved)
    expect_length(w, 0L)
    b <- assign_features_to_polygons(.at_feat(), tri[ord, ], polygon_id_col = "pid")
    expect_identical(b$pid, "T1", info = paste(ord, collapse = ""))
    expect_identical(attr(b, "ties")$which, 1L)
  }
  # "first" keeps the first polygon in row order.
  expect_identical(.at_quiet(assign_features_to_polygons(.at_pts(), tri[2:1, ],
                     polygon_id_col = "pid", tie_break = "first"))$pid, rep("T2", 3))
  expect_identical(.at_quiet(assign_features_to_polygons(.at_feat(), tri[2:1, ],
                     polygon_id_col = "pid", tie_break = "first"))$pid, "T2")
})

test_that("polygons equal on every key leave the feature without a polygon, with a warning, in both row orders", {
  cases <- list(
    rhombus = sf::st_sf(pid = c("square", "rhombus"),
                        geometry = sf::st_sfc(.at_sq(), .at_rh(), crs = 32617)),
    copy = sf::st_sf(pid = c("A", "B"), geometry = sf::st_sfc(.at_sq(), .at_sq(), crs = 32617)))
  for (nm in names(cases)) for (ord in list(1:2, 2:1)) {
    polys <- cases[[nm]][ord, ]
    info <- paste(nm, paste(ord, collapse = ""))
    w <- NULL
    a <- withCallingHandlers(
      assign_features_to_polygons(.at_pts(), polys, polygon_id_col = "pid",
                                  keep_unassigned = TRUE),
      spatialkit_ties_unresolved = function(cnd) { w <<- cnd; invokeRestart("muffleWarning") })
    expect_identical(a$pid, rep(NA_character_, 3), info = info)
    expect_s3_class(w, "warning")
    expect_identical(w$n, 3L); expect_identical(w$total, 3L)
    expect_match(conditionMessage(w), paste0("3 of 3 feature\\(s\\) match two or more polygons ",
                                             "that the smallest-area rule"))
    expect_match(conditionMessage(w), "kept with an NA ID", fixed = TRUE)
    ties <- attr(a, "ties")
    expect_identical(ties$which, 1:3, info = info)
    expect_identical(ties$unresolved, 1:3, info = info)
    # They matched polygons, so they are not counted as in none.
    expect_identical(attr(a, "unassigned")$n, 0L, info = info)
    # Dropped with keep_unassigned = FALSE, with the ties warning alone.
    ws <- character(0)
    d <- withCallingHandlers(
      assign_features_to_polygons(.at_pts(), polys, polygon_id_col = "pid"),
      warning = function(cnd) { ws <<- c(ws, class(cnd)[1]); invokeRestart("muffleWarning") })
    expect_identical(nrow(d), 0L, info = info)
    expect_identical(ws, "spatialkit_ties_unresolved", info = info)
    expect_identical(attr(d, "ties")$unresolved, 1:3, info = info)
    # The largest-overlap join, likewise.
    w <- NULL
    b <- withCallingHandlers(
      assign_features_to_polygons(.at_feat(), polys, polygon_id_col = "pid",
                                  keep_unassigned = TRUE),
      spatialkit_ties_unresolved = function(cnd) { w <<- cnd; invokeRestart("muffleWarning") })
    expect_identical(b$pid, NA_character_, info = info)
    expect_identical(w$n, 1L)
    expect_identical(attr(b, "ties")$unresolved, 1L)
    expect_identical(attr(b, "unassigned")$n, 0L)
    # "first" keeps the first polygon in row order and warns of nothing.
    expect_no_warning(f <- assign_features_to_polygons(.at_pts(), polys, polygon_id_col = "pid",
                                                       tie_break = "first"))
    expect_identical(f$pid, rep(polys$pid[1], 3), info = info)
  }
  # A point in one of the two only is assigned, and the count is of the rest.
  pts <- rbind(.at_pts(), sf::st_sf(r = 4L, geometry = sf::st_sfc(sf::st_point(c(140, 50)),
                                                                 crs = 32617)))
  expect_warning(a <- assign_features_to_polygons(pts, cases$rhombus, polygon_id_col = "pid",
                                                  keep_unassigned = TRUE),
                 "3 of 4 feature", class = "spatialkit_ties_unresolved")
  expect_identical(a$pid, c(NA, NA, NA, "rhombus"))
})

test_that("grid cells that share an edge are told apart by their centres, as before", {
  cells <- .at_two_cells(0, 0)
  pe <- sf::st_sf(r = 1L, geometry = sf::st_sfc(sf::st_point(c(100, 50)), crs = 32617))
  for (ord in list(1:2, 2:1))
    expect_identical(assign_features_to_polygons(pe, cells[ord, ], polygon_id_col = "pid")$pid,
                     "C1")
  # A 20 m square across the edge: equal overlaps, the left cell.
  sq <- sf::st_sf(r = 1L, geometry = sf::st_sfc(.at_P(rbind(
    c(90, 40), c(110, 40), c(110, 60), c(90, 60), c(90, 40))), crs = 32617))
  for (ord in list(1:2, 2:1))
    expect_identical(assign_features_to_polygons(sq, cells[ord, ], polygon_id_col = "pid")$pid,
                     "C1")
})

test_that(".tie_rule_rows() keeps the least area, then keys within the margin, by group", {
  key <- rbind(c(10, 5, 1, 0, 0),        # group 1: area decides
               c(20, 0, 0, 0, 0),
               c(10, 5, 2, 0, 0),        # group 2: centre y within tol, x decides
               c(10, 5 + 1e-9, 1, 0, 0),
               c(10, 5, 1, 3, 3),        # group 3: equal on every key
               c(10, 5, 1, 3, 3),
               c(Inf, NA, NA, NA, NA))   # group 4: alone
  grp <- c(1, 1, 2, 2, 3, 3, 4)
  expect_identical(.tie_rule_rows(key, grp, tol = 1e-6),
                   c(TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE))
  # Without the margin the lower centre y wins.
  expect_identical(.tie_rule_rows(key, grp, tol = 0)[3:4], c(TRUE, FALSE))
})

test_that("?assign_features_to_polygons states the rule", {
  root <- testthat::test_path("..", "..")
  if (!(file.exists(file.path(root, "DESCRIPTION")) && dir.exists(file.path(root, "man"))))
    skip("the help's source is not here to read")
  read <- function(...) gsub("[[:space:]]+", " ",
                             paste(readLines(file.path(root, ...), warn = FALSE), collapse = " "))
  rd <- read("man", "assign_features_to_polygons.Rd")
  expect_match(rd, "then the one whose vertices' mean is lowest, then leftmost", fixed = TRUE)
  expect_match(rd, "\\code{\"spatialkit_ties_unresolved\"}", fixed = TRUE)
  expect_match(rd, "a POINT or MULTIPOINT feature with a point within a margin of a polygon is in it",
               fixed = TRUE)
  # The centres and means are compared within the margin itself, also where
  # the predicate's margin is none (degrees with s2 off, another predicate).
  expect_match(rd, paste0("centres and means within one part in 1e10 of the largest coordinate ",
                          "of the polygons, at most a millionth of the diagonal of their ",
                          "bounding box and at least four times the spacing of doubles at that ",
                          "coordinate (with any \\code{predicate}, and in longitude and ",
                          "latitude whatever \\code{sf::sf_use_s2()} says)"), fixed = TRUE)
  expect_false(grepl("within the margin of \\code{predicate}", rd, fixed = TRUE))
})

test_that("a MULTIPOINT within the margin of an edge goes where the POINT goes; a line gets no margin", {
  cells <- .at_two_cells()
  e <- 500100 + 1e-5
  geoms <- sf::st_sfc(
    sf::st_point(c(e, 4000050)),
    sf::st_multipoint(matrix(c(e, 4000050), 1)),
    # One point within the margin of the edge, one inside C2: in both cells.
    sf::st_multipoint(rbind(c(e, 4000050), c(500150, 4000050))),
    # A line 1e-5 m east of the edge lies in C2 alone, as sf::st_intersects() says.
    sf::st_linestring(rbind(c(e, 4000010), c(e, 4000090))),
    sf::st_point(c(500150, 4000050)), crs = 32617)
  f <- sf::st_sf(r = 1:5, geometry = geoms)
  for (ord in list(1:2, 2:1)) {
    info <- paste(ord, collapse = "")
    # Three features and two cells (the cells are buffered by the margin),
    a <- .at_quiet(assign_features_to_polygons(f, cells[ord, ], polygon_id_col = "pid"))
    expect_identical(a$pid, c("C1", "C1", "C1", "C2", "C2"), info = info)
    expect_identical(attr(a, "ties")$which, 1:3, info = info)
    # and one feature at a time (the feature is buffered).
    one <- vapply(1:5, function(k) .at_quiet(assign_features_to_polygons(
      f[k, ], cells[ord, ], polygon_id_col = "pid"))$pid, character(1))
    expect_identical(one, c("C1", "C1", "C1", "C2", "C2"), info = info)
  }
})

test_that("among cells of one area the lowest centre wins before the leftmost", {
  # A lower-right and an upper-left 100 m square meeting at the corner
  # (100, 100): the lower-right centre (150, 50) is the lowest, the
  # upper-left one (50, 150) the leftmost.
  B <- function(x, y) .at_P(rbind(c(x, y), c(x + 100, y), c(x + 100, y + 100), c(x, y + 100),
                                  c(x, y)))
  cells <- sf::st_sf(pid = c("Lr", "Ul"), geometry = sf::st_sfc(B(100, 0), B(0, 100), crs = 32617))
  pt <- sf::st_sf(r = 1L, geometry = sf::st_sfc(sf::st_point(c(100, 100)), crs = 32617))
  # A 20 m square about the corner overlaps each by 100 m^2.
  sq <- sf::st_sf(r = 1L, geometry = sf::st_sfc(.at_P(rbind(
    c(90, 90), c(110, 90), c(110, 110), c(90, 110), c(90, 90))), crs = 32617))
  for (ord in list(1:2, 2:1)) {
    info <- paste(ord, collapse = "")
    a <- assign_features_to_polygons(pt, cells[ord, ], polygon_id_col = "pid")
    expect_identical(a$pid, "Lr", info = info)
    expect_identical(attr(a, "ties")$which, 1L, info = info)
    b <- assign_features_to_polygons(sq, cells[ord, ], polygon_id_col = "pid")
    expect_identical(b$pid, "Lr", info = info)
    expect_identical(attr(b, "ties")$which, 1L, info = info)
  }
})

test_that("a point in a polygon whose buffer GEOS draws wrongly keeps its exact polygon", {
  # A bowtie ring (invalid): the point (10, 50) in its left lobe is in it by
  # sf::st_intersects() and in no polygon after st_buffer(, 1e-6).  With more
  # points than polygons the polygons are buffered for the margin.
  bow <- .at_P(rbind(c(0, 0), c(100, 100), c(100, 0), c(0, 100), c(0, 0)))
  sq  <- .at_P(rbind(c(200, 0), c(300, 0), c(300, 100), c(200, 100), c(200, 0)))
  pl  <- sf::st_sf(pid = c("bow", "sq"), geometry = sf::st_sfc(bow, sq, crs = 32617))
  pts <- sf::st_sf(r = 1:3, geometry = sf::st_sfc(
    sf::st_point(c(10, 50)), sf::st_point(c(90, 50)), sf::st_point(c(250, 50)), crs = 32617))
  expect_identical(lengths(sf::st_intersects(pts, pl)), c(1L, 1L, 1L))
  expect_identical(lengths(sf::st_intersects(pts, sf::st_buffer(pl, 1e-6)))[1], 0L)
  for (ord in list(1:2, 2:1)) {
    expect_no_warning(a <- suppressMessages(assign_features_to_polygons(
      pts, pl[ord, ], polygon_id_col = "pid", keep_unassigned = TRUE)))
    expect_identical(a$pid, c("bow", "bow", "sq"), info = paste(ord, collapse = ""))
  }
})

test_that("centres closer than the margin are equal in degrees with s2 on and off", {
  # Two unit squares in degrees, one moved north by 1e-11 degrees: the
  # centres differ by less than one part in 1e10 of the largest coordinate,
  # so the point in both is unresolved whatever sf::sf_use_s2() says.
  sq <- .at_P(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0)))
  ov <- sf::st_sf(pid = c("lo", "hi"), geometry = sf::st_sfc(
    sq, sq + c(0, 1e-11), crs = 4326))
  pt <- sf::st_sf(r = 1L, geometry = sf::st_sfc(sf::st_point(c(0.5, 0.5)), crs = 4326))
  old <- sf::sf_use_s2()
  on.exit(suppressMessages(sf::sf_use_s2(old)), add = TRUE)
  for (s2 in c(TRUE, FALSE)) {
    suppressMessages(sf::sf_use_s2(s2))
    for (ord in list(1:2, 2:1)) {
      info <- paste("s2", s2, paste(ord, collapse = ""))
      expect_warning(a <- suppressMessages(assign_features_to_polygons(
        pt, ov[ord, ], polygon_id_col = "pid", keep_unassigned = TRUE)),
        class = "spatialkit_ties_unresolved")
      expect_identical(a$pid, NA_character_, info = info)
    }
  }
})
