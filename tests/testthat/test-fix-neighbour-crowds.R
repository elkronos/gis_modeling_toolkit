# BUG-30: a point that holds more cells' centroids than a cell has neighbours
# is found with the reach of a tie, not by exact equality.
#
# A cell whose centroid shares its point with eight other cells or more takes
# the first eight of them without a search.  The point was found by exact
# equality.  In another unit, or with the origin moved, rounding parts such
# centroids by a few units in the last place, and the cells were then
# searched for among themselves by those distances: which eight each took,
# and with them every pooled value, followed the unit and the origin.  Now
# centroids within a ten-millionth of the largest coordinate of one another,
# in x and in y, are one point, and its cells take each other by number.

fnc_paths <- function() c(FALSE, if (requireNamespace("FNN", quietly = TRUE)) TRUE)
fnc_tag   <- function(f) if (f) "with FNN" else "exact search"
fnc_quiet <- function(expr)
  logger::with_log_threshold(expr, threshold = logger::FATAL,
                             namespace = "spatialkit", index = 2)

# A 10 x 10 lattice of one-row cells at 250 m from (1e5, 1e5), `cc` cells
# of 1, 2, ..., cc rows with every row on the point `pt`, and a one-row cell
# at each row of `extra`.
fnc_crowd <- function(pt, cc, origin = c(1e5, 1e5), extra = matrix(numeric(0), 0L, 2L)) {
  lat <- unname(as.matrix(expand.grid(0:9 * 250, 0:9 * 250)))
  lat <- lat + rep(origin, each = nrow(lat))
  xy <- rbind(lat, matrix(rep(pt, each = sum(seq_len(cc))), ncol = 2), extra)
  gi <- c(seq_len(nrow(lat)), nrow(lat) + rep(seq_len(cc), seq_len(cc)),
          nrow(lat) + cc + seq_len(nrow(extra)))
  list(xy = xy, gi = gi, K = max(gi), crowd = nrow(lat) + seq_len(cc))
}
fnc_frames <- list(feet = function(xy) xy / 0.3048,
                   `US feet` = function(xy) xy * 3937 / 1200,
                   `origin moved by 0.1` = function(xy) xy + 0.1)


test_that("ten or more cells on one point keep their neighbours and values in another unit or origin", {
  cases <- list(
    list("twelve cells on a lattice point", c(1e5 + 1000, 1e5 + 1000), 12L, 5),
    list("twelve cells between lattice points", c(1e5 + 1125, 1e5 + 1125), 12L, 62),
    list("twenty cells on a lattice point", c(1e5 + 1000, 1e5 + 1000), 20L, 25))
  for (case in cases) {
    L <- fnc_crowd(case[[2L]], case[[3L]])
    y <- withr::with_seed(case[[4L]], stats::rnorm(nrow(L$xy)))
    A <- .pool_fit(y, L$gi, L$K, L$xy, use_fnn = FALSE)
    expect_identical(nrow(unique(A$cen[L$crowd, ])), 1L, info = case[[1L]])
    if (case[[3L]] == 20L) {
      expect_identical(A$status, "ok", info = case[[1L]])
      expect_gt(A$tau2, 0)
    }
    parted <- 0L
    for (f in fnc_paths()) for (nm in names(fnc_frames)) {
      tag <- sprintf("%s, %s, %s", case[[1L]], nm, fnc_tag(f))
      B <- .pool_fit(y, L$gi, L$K, fnc_frames[[nm]](L$xy), use_fnn = f)
      parted <- parted + (nrow(unique(B$cen[L$crowd, ])) > 1L)
      expect_identical(B$nb, A$nb, info = tag)
      for (v in c("value", "target", "w", "tau2", "loo_pooled"))
        expect_equal(B[[v]], A[[v]], tolerance = 1e-12, info = paste(tag, v))
      expect_identical(B$status, A$status, info = tag)
    }
    expect_gt(parted, 0L)
  }
})

test_that("a crowd beside centroids that chain on beyond the reach keeps its neighbours and values in another unit or origin", {
  # The crowd of twelve between lattice points at UTM magnitudes (a reach of
  # about 0.46 m), and cells 0.3 m and 0.6 m from it along x, or 0.3 m on
  # either side: a chain of pairs within the reach that spans more than it.
  # At half the reach the crowd is one point and the two cells are apart.
  P <- c(4.5e5 + 1125, 4.55e6 + 1125)
  sides <- list(`0.3 m and 0.6 m along x` = rbind(P + c(0.3, 0), P + c(0.6, 0)),
                `0.3 m on either side` = rbind(P + c(0.3, 0), P - c(0.3, 0)))
  for (nm in names(sides)) {
    L <- fnc_crowd(P, 12L, origin = c(4.5e5, 4.55e6), extra = sides[[nm]])
    y <- withr::with_seed(5, stats::rnorm(nrow(L$xy)))
    A <- .pool_fit(y, L$gi, L$K, L$xy, use_fnn = FALSE)
    expect_identical(nrow(unique(A$cen[L$crowd, ])), 1L, info = nm)
    expect_lt(.TIE_TOL * max(abs(L$xy)), 0.6)
    expect_gt(.TIE_TOL * max(abs(L$xy)), 0.3)
    for (i in L$crowd)
      expect_identical(A$nb[i, ], setdiff(L$crowd, i)[1:8], info = sprintf("%s, region %d", nm, i))
    parted <- 0L
    for (f in fnc_paths()) for (fr in names(fnc_frames)) {
      tag <- sprintf("%s, %s, %s", nm, fr, fnc_tag(f))
      B <- .pool_fit(y, L$gi, L$K, fnc_frames[[fr]](L$xy), use_fnn = f)
      parted <- parted + (nrow(unique(B$cen[L$crowd, ])) > 1L)
      expect_identical(B$nb, A$nb, info = tag)
      for (v in c("value", "target", "w", "tau2", "loo_pooled"))
        expect_equal(B[[v]], A[[v]], tolerance = 1e-12, info = paste(tag, v))
    }
    expect_gt(parted, 0L)
  }
})

test_that("a crowd parted by rounding takes its first neighbours by number", {
  # Twelve centroids on (5e5, 5e6), regions 1 to 12, beside a 6 x 6 lattice
  # at 100 m, regions 13 to 48, whose point 15, region 27, is that point too.
  # Regions 1 to 4 are pushed 2 units in the last place along x and 7 and 10
  # 1 unit down along y, so the lattice's column at x = 5e5 comes between
  # parts of the crowd in coordinate order.
  lat <- unname(as.matrix(expand.grid(-2:3 * 100, -2:3 * 100))) +
    rep(c(5e5, 5e6), each = 36L)
  on  <- which(lat[, 1] == 5e5 & lat[, 2] == 5e6)
  expect_identical(on, 15L)
  crowd <- matrix(c(5e5, 5e6), 12L, 2L, byrow = TRUE)
  cen0 <- rbind(crowd, lat)
  cen  <- cen0
  ulp  <- function(v, n) v * (1 + n * .Machine$double.eps)
  cen[1:4, 1]      <- ulp(cen[1:4, 1], 2)
  cen[c(7, 10), 2] <- ulp(cen[c(7, 10), 2], -1)
  expect_identical(nrow(unique(cen[1:12, ])), 3L)
  ox <- order(cen[, 1], cen[, 2])
  expect_false(all(diff(sort(match(1:12, ox))) == 1L))
  point <- c(1:12, 12L + on)
  for (f in fnc_paths()) {
    a <- .region_neighbours(cen0, 8L, f)
    b <- .region_neighbours(cen, 8L, f)
    expect_identical(b, a, info = fnc_tag(f))
    for (i in point)
      expect_identical(b[i, ], setdiff(point, i)[1:8], info = sprintf("%s, region %d", fnc_tag(f), i))
  }
})

test_that("summarize_by_cell() keeps the pooled values of a crowd in another unit or origin", {
  pooled <- function(layer) sf::st_drop_geometry(fnc_quiet(suppressWarnings(
    summarize_by_cell(layer, "z", id_col = "cell_id", estimate = "pooled"))))
  usft <- sf::st_crs("+proj=utm +zone=15 +datum=NAD83 +units=us-ft +no_defs")
  # UTM zone 15 with its false easting and northing moved by 0.1 m.
  moved <- sf::st_crs(paste("+proj=tmerc +lat_0=0 +lon_0=-93 +k=0.9996 +x_0=500000.1",
                            "+y_0=0.1 +datum=NAD83 +units=m +no_defs"))
  L  <- fnc_crowd(c(4.5e5 + 1000, 4.55e6 + 1000), 12L, origin = c(4.5e5, 4.55e6))
  z  <- withr::with_seed(5, stats::rnorm(nrow(L$xy)))
  pm <- sf::st_as_sf(data.frame(x = L$xy[, 1], y = L$xy[, 2], cell_id = L$gi, z = z),
                     coords = c("x", "y"), crs = 26915)
  ncen <- function(layer) {
    xy <- unname(sf::st_coordinates(layer))
    cen <- .pool_geometry(layer$cell_id, L$K, xy, tabulate(layer$cell_id, L$K), FALSE)$cen
    nrow(unique(cen[L$crowd, ]))
  }
  expect_identical(ncen(pm), 1L)
  a <- pooled(pm)
  cols <- grep("^\\.\\.pooled", names(a), value = TRUE)
  expect_true("..pooled_resp_z" %in% cols)
  parted <- 0L
  for (crs in list(usft, moved)) {
    b <- sf::st_transform(pm, crs)
    parted <- parted + (ncen(b) > 1L)
    bb <- pooled(b)
    expect_identical(bb$cell_id, a$cell_id, info = crs$input)
    for (cc in cols) expect_equal(bb[[cc]], a[[cc]], tolerance = 1e-12, info = paste(crs$input, cc))
  }
  expect_gt(parted, 0L)
})

test_that("centroids within the reach count as one point, and centroids beyond it are searched", {
  # A 10 x 10 lattice at 250 m at UTM magnitudes, regions 1 to 100, and
  # twelve regions, 101 to 112, from lattice point 45 along x.  The reach is
  # a ten-millionth of the largest coordinate, about 0.46 m here.
  lat <- unname(as.matrix(expand.grid(0:9 * 250, 0:9 * 250))) +
    rep(c(4.5e5, 4.55e6), each = 100L)
  line <- function(step) cbind(lat[45L, 1] + step * 0:11, lat[45L, 2])
  near <- rbind(lat, line(0.04))
  expect_lt(0.44, .TIE_TOL * max(abs(near)))
  point <- c(45L, 101:112)
  # The twelve, 0.04 m apart, and lattice point 45 are one point of
  # thirteen: each takes the first eight others by number, in metres and in
  # feet.
  for (f in fnc_paths()) for (s in c(1, 1 / 0.3048)) {
    nb <- .region_neighbours(near * s, 8L, f)
    for (i in point)
      expect_identical(nb[i, ], setdiff(point, i)[1:8],
                       info = sprintf("%s, scale %g, region %d", fnc_tag(f), s, i))
  }
  # One metre apart they are searched for: region 112 takes the nearest.
  # So are they 0.3 m apart, each within the reach of the next, since the
  # chain spans 3.3 m.
  for (step in c(1, 0.3)) for (f in fnc_paths()) {
    nb <- .region_neighbours(rbind(lat, line(step)), 8L, f)
    expect_identical(nb[112L, ], 111:104, info = sprintf("%s, %g m", fnc_tag(f), step))
    expect_identical(nb[101L, ], c(45L, 102:108), info = sprintf("%s, %g m", fnc_tag(f), step))
  }
})

test_that("with the origin moved far the reach shrinks, and centroids within it without coinciding are searched", {
  # The twelve centroids 0.04 m apart from lattice point 45 at UTM
  # magnitudes are one point of thirteen.  With the origin at the corner of
  # the layer the reach is about 0.000225 m, the twelve lie beyond it, and
  # each takes the nearest.  The help and NEWS state this limit.
  lat <- unname(as.matrix(expand.grid(0:9 * 250, 0:9 * 250))) +
    rep(c(4.5e5, 4.55e6), each = 100L)
  near <- rbind(lat, cbind(lat[45L, 1] + 0.04 * 0:11, lat[45L, 2]))
  local <- near - rep(c(4.5e5, 4.55e6), each = nrow(near))
  expect_lt(.TIE_TOL * max(abs(local)), 0.04)
  for (f in fnc_paths()) {
    expect_identical(.region_neighbours(near, 8L, f)[112L, ], c(45L, 101:107), info = fnc_tag(f))
    expect_identical(.region_neighbours(local, 8L, f)[112L, ], 111:104, info = fnc_tag(f))
  }
  root <- testthat::test_path("..", "..")
  skip_if_not(file.exists(file.path(root, "DESCRIPTION")) &&
                file.exists(file.path(root, "NEWS.md")), "not run from the source tree")
  squash <- function(s) gsub("[[:space:]]+", " ", paste(s, collapse = " "))
  rd   <- squash(readLines(file.path(root, "man", "summarize_by_cell.Rd")))
  news <- squash(readLines(file.path(root, "NEWS.md")))
  for (txt in list(rd, news)) {
    expect_match(txt, "with the origin moved far", fixed = TRUE)
    expect_match(txt, "taken again at half", fixed = TRUE)
  }
  expect_false(grepl("get the same values, in any unit or origin", news, fixed = TRUE))
  expect_false(grepl("no other centroid lies that close to any of them", rd, fixed = TRUE))
})

test_that(".shared_points() joins chains within the reach that span no more than it", {
  # The rule written out: centroids that chains of pairs within h in x and
  # in y link, while the chain spans no more than h along x and along y;
  # the centroids of a wider chain are linked again at h / 2, and so on,
  # until no chain is wider than the reach it was linked at.
  rule <- function(x, y, h) {
    lab <- match(paste(x, y), paste(x, y))
    at <- seq_along(x); top <- length(x)
    while (length(at) > 1L && h > 0) {
      xa <- x[at]; ya <- y[at]; n <- length(at)
      near <- abs(outer(xa, xa, "-")) <= h & abs(outer(ya, ya, "-")) <= h
      ch <- seq_len(n)
      repeat {
        c2 <- vapply(seq_len(n), function(i) min(ch[near[i, ]]), integer(1))
        if (identical(c2, ch)) break
        ch <- c2
      }
      wide <- logical(n)
      for (g in unique(ch)) {
        i <- which(ch == g)
        wide[i] <- diff(range(xa[i])) > h || diff(range(ya[i])) > h
      }
      lab[at[!wide]] <- top + ch[!wide]
      top <- top + n; at <- at[wide]; h <- h / 2
    }
    lab
  }
  same <- function(p, q) identical(outer(p, p, "=="), outer(q, q, "=="))
  # With no reach a point is the centroids that are equal.
  x <- c(3, 1, 3, 2, 1, 3); y <- c(1, 1, 1, 5, 1, 2)
  p <- .shared_points(x, y, 0)
  expect_true(same(p, c(1, 2, 1, 3, 2, 4)))
  expect_identical(.shared_points(1, 1, 0), 1L)
  expect_identical(.shared_points(numeric(0), numeric(0), 0), integer(0))
  expect_identical(.shared_points(c(1, 1), c(1, 1), 0), c(1L, 1L))
  expect_identical(length(unique(.shared_points(c(1, 2), c(1, 1), 0.5))), 2L)
  # A gap just above the reach parts two centroids, along x or along y; at
  # the reach they stay together.
  h <- 0.25
  expect_identical(length(unique(.shared_points(c(0, h * (1 + 1e-9)), c(0, 0), h))), 2L)
  expect_identical(length(unique(.shared_points(c(0, 0), c(0, h * (1 + 1e-9)), h))), 2L)
  expect_identical(length(unique(.shared_points(c(0, h), c(0, h), h))), 1L)
  expect_identical(length(unique(.shared_points(c(0, h), c(0, 0), h))), 1L)
  expect_identical(length(unique(.shared_points(c(0, 0), c(0, h), h))), 1L)
  # A pair within the reach in x but not in y is two points, though other
  # centroids lie between them in order of x.
  expect_identical(length(unique(.shared_points(c(0, 0.1, 0.05, 0.05), c(0, 0, 9, 9), h))), 2L)
  # Four centroids whose gaps along x and along y are all h, but no two of
  # which lie within h of each other in x and in y: four points.
  expect_identical(length(unique(.shared_points(0:3, c(2, 4, 1, 3), 1))), 4L)
  # A chain of pairs within the reach that spans more than it is no point,
  # but the equal centroids on it stay together.
  p <- .shared_points(c(0, 0.2, 0.4, 0.4), c(0, 0, 0, 0), h)
  expect_true(same(p, c(1, 2, 3, 3)))
  expect_identical(length(unique(.shared_points(c(0, 0.1, 0.2), c(0, 0, 0), h))), 1L)
  # Its centroids are linked again at half the reach, and again, so
  # centroids parted by rounding on it stay one point.
  x <- c(5, 5 * (1 + 2 * .Machine$double.eps), 5, 5.2, 5.4, 5.45)
  p <- .shared_points(x, rep(0, 6), h)
  expect_true(same(p, c(1, 1, 1, 2, 3, 3)))
  # With 5.1 added the chain at half the reach still spans 0.2, more than
  # it: at a quarter of the reach the crowd is one point.
  p <- .shared_points(c(x, 5.1), rep(0, 7), h)
  expect_true(same(p, c(1, 1, 1, 2, 3, 3, 4)))
  # A slightly rotated lattice at 20 h: every centroid its own point.
  g <- as.matrix(expand.grid(0:29, 0:29)) * 20 * h
  g <- g %*% matrix(c(cos(0.01), sin(0.01), -sin(0.01), cos(0.01)), 2)
  expect_identical(length(unique(.shared_points(g[, 1], g[, 2], h))), 900L)
  # Random centroids, with near and equal copies, and copies parted by
  # rounding: the rule's points.
  for (seed in 1:8) {
    b <- withr::with_seed(seed, cbind(stats::runif(150, 0, 10), stats::runif(150, 0, 10)))
    d <- withr::with_seed(seed + 100, b[sample(150, 100, replace = TRUE), ] +
                            stats::runif(200, -0.15, 0.15))
    r <- d[1:30, ] * (1 + withr::with_seed(seed + 200, sample(-3:3, 60, replace = TRUE)) *
                        .Machine$double.eps)
    v <- rbind(b, d, d[1:20, ], r)
    expect_true(same(.shared_points(v[, 1], v[, 2], h), rule(v[, 1], v[, 2], h)),
                info = paste("seed", seed))
  }
})
