# BUG-15: on gridded locations the variogram's short-lag pairs and the
# pooling's neighbours follow one tie rule.
#
# estimate_sac_range() pairs each location with its nearest other location,
# and on a lattice two to four are equally near.  FNN's kd-tree and the
# search in base R took different ones (80 pairs against 143 on the frozen
# fixture's 12 x 12 lattice), and Cp's noise variance, held to the
# semivariance of those pairs, followed.  The pooling's searches counted only
# exact ties, and centroids of a lattice that tie exactly in one unit differ
# in their last bits in another, so the same layer in metres and in feet, or
# with its origin moved, gave some cells other neighbours.  Now two squared
# distances within one part in ten million (.TIE_TOL) of each other are a
# tie, which goes to the lower x (x values within that part of the distance
# counting as equal), then the lower y, then the lower number.  Between
# regions that number is the region's own, so regions whose centroids
# coincide in one unit and part by rounding in another keep their
# neighbours.

fnt_paths <- function() c(FALSE, if (requireNamespace("FNN", quietly = TRUE)) TRUE)
fnt_tag   <- function(f) if (f) "with FNN" else "exact search"

# The package's functions as they run where FNN cannot be loaded, for the
# rest of the calling frame.
fnt_without_fnn <- function(env = parent.frame())
  testthat::local_mocked_bindings(!!!core_frozen_without("FNN"), .package = "spatialkit",
                                  .env = env)
fnt_quiet <- function(expr)
  logger::with_log_threshold(expr, threshold = logger::FATAL,
                             namespace = "spatialkit", index = 2)

# The frozen fixture's lattice: 12 x 12 locations 80 m apart, the rows in
# coordinate order, as estimate_sac_range() hands them on.
fnt_lattice <- function() {
  g <- expand.grid(gx = seq(40, 920, by = 80), gy = seq(40, 920, by = 80))
  z <- withr::with_seed(1, sin(g$gx / 200) + cos(g$gy / 150) + stats::rnorm(nrow(g), sd = 0.3))
  xy <- cbind(5e5 + g$gx, 5e6 + g$gy)
  o <- order(xy[, 1], xy[, 2])
  list(xy = unname(xy[o, ]), z = z[o])
}

# The rule written out, for one location at a time: of the candidates `s`,
# those within one part in 1e7 of the least squared distance; of them those
# whose x is within 1e-7 of the distance of the least x; of those the ones
# whose y is; and of those the first in coordinate order.
fnt_pick <- function(s, d2, xy) {
  m <- min(d2[s])
  s <- s[d2[s] <= m * (1 + 1e-7)]
  s <- s[xy[s, 1] <= min(xy[s, 1]) + 1e-7 * sqrt(m)]
  s <- s[xy[s, 2] <= min(xy[s, 2]) + 1e-7 * sqrt(m)]
  s[order(xy[s, 1], xy[s, 2])][1L]
}
# Each location's nearest other location by the rule, and the pairs of
# locations that makes, each once, in the order of the locations.
fnt_rule_pairs <- function(xy) {
  n <- nrow(xy)
  j <- vapply(seq_len(n), function(i) {
    d2 <- (xy[, 1] - xy[i, 1])^2 + (xy[, 2] - xy[i, 2])^2
    fnt_pick(seq_len(n)[-i], d2, xy)
  }, integer(1))
  pr <- cbind(pmin(seq_len(n), j), pmax(seq_len(n), j))
  pr[!duplicated(pr), , drop = FALSE]
}

# Regions of a k-means fit at `K` centres, each row in the region of its
# nearest centre, numbered 1 to K.
fnt_regions <- function(xy, K, seed) {
  withr::local_seed(seed, .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion",
                    .rng_sample_kind = "Rejection")
  km <- suppressWarnings(stats::kmeans(xy, K, nstart = 3, iter.max = 100))
  g  <- .nearest_centre(xy[, 1], xy[, 2], km$centers[, 1], km$centers[, 2])
  match(g, sort(unique(g)))
}


test_that("the tolerance of a tie is one part in ten million", {
  expect_identical(.TIE_TOL, 1e-7)
})

test_that("on a lattice the short-lag pairs follow the rule, with and without FNN", {
  L <- fnt_lattice(); xy <- L$xy; z <- L$z
  a <- .sac_short_lag(xy, z)
  expect_identical(a$pairs, 143L)
  expect_identical(a$dist, 80)
  pr <- fnt_rule_pairs(xy)
  expect_identical(a$pairs, nrow(pr))
  expect_identical(a$gamma, mean(0.5 * (z[pr[, 1]] - z[pr[, 2]])^2))
  # Without FNN: the same list, to the last bit.
  expect_identical(local({ fnt_without_fnn(); .sac_short_lag(xy, z) }), a)
  # The rows in another order: the same pairs, summed in another order.
  p <- withr::with_seed(2, sample(nrow(xy)))
  r <- .sac_short_lag(xy[p, ], z[p])
  expect_identical(r$pairs, a$pairs)
  expect_equal(r$gamma, a$gamma, tolerance = 1e-14)
  expect_identical(r$dist, a$dist)
  expect_identical(local({ fnt_without_fnn(); .sac_short_lag(xy[p, ], z[p]) }), r)
  # In feet, and with the origin moved: the same pairs in the same order.
  for (case in list(list("feet", xy / 0.3048, 80 / 0.3048),
                    list("origin moved by 0.1", xy + 0.1, 80),
                    list("origin at (0, 0)", sweep(xy, 2L, c(5e5, 5e6)), 80))) {
    s <- .sac_short_lag(case[[2L]], z)
    expect_identical(s$pairs, a$pairs, info = case[[1L]])
    expect_identical(s$gamma, a$gamma, info = case[[1L]])
    expect_equal(s$dist, case[[3L]], tolerance = 1e-12, info = case[[1L]])
    expect_identical(local({ fnt_without_fnn(); .sac_short_lag(case[[2L]], z) }), s,
                     info = case[[1L]])
  }
  # The rule prefers the lower x, so with the axes swapped other pairs come
  # out; they are the rule's on the swapped layer.
  sw <- xy[, 2:1]
  o  <- order(sw[, 1], sw[, 2])
  b  <- .sac_short_lag(sw[o, ], z[o])
  ps <- fnt_rule_pairs(sw[o, ])
  expect_identical(b$pairs, nrow(ps))
  expect_identical(b$gamma, mean(0.5 * (z[o][ps[, 1]] - z[o][ps[, 2]])^2))
})

test_that("without FNN the short-lag search is not stopped by the work limit", {
  L <- fnt_lattice()
  a <- .sac_short_lag(L$xy, L$z)
  withr::local_options(spatialkit.exact_pairs_max = 0)
  expect_identical(local({ fnt_without_fnn(); .sac_short_lag(L$xy, L$z) }), a)
})

test_that("estimate_sac_range() gives a lattice the same short-lag pairs with and without FNN", {
  skip_if_not_installed("gstat")
  L <- fnt_lattice()
  layer <- sf::st_as_sf(data.frame(x = L$xy[, 1], y = L$xy[, 2], z = L$z),
                        coords = c("x", "y"), crs = 32632)
  fit <- function() fnt_quiet(suppressWarnings(estimate_sac_range(layer, "z")))
  s1 <- fit()
  s2 <- local({ fnt_without_fnn(); fit() })
  at <- c("short_lag_pairs", "short_lag_semivariance", "short_lag_dist")
  expect_identical(attributes(s1)[at], attributes(s2)[at])
  expect_identical(attr(s1, "short_lag_pairs"), 143L)
  expect_identical(attr(s1, "short_lag_dist"), 80)
})

test_that("the frozen fixture holds one record of the lattice calls, with and without FNN", {
  base <- readRDS(test_path("fixtures", "core-frozen.rds"))
  lattice <- c("profile_lattice", "profile_lattice_density", "profile_integer",
               "select_lattice_cp", "summary_lattice", "select_lattice_density_cp",
               "summary_lattice_density", "select_integer_cp", "summary_integer")
  expect_true(all(lattice %in% names(base$records)))
  expect_false(any(lattice %in% names(base$without$records)))
  for (nm in c("profile_lattice", "profile_lattice_density", "profile_integer")) {
    sac <- attr(base$records[[nm]]$value, "sac")
    expect_identical(attr(sac, "short_lag_pairs"), 143L, info = nm)
  }
})

test_that("on a lattice the pooling's neighbours and values do not move with the unit or the origin", {
  # A 250 m lattice of 60 x 50 locations in whole metres, in 1,500 regions
  # from k-means.
  latm <- unname(as.matrix(expand.grid(x = 0:59 * 250, y = 0:49 * 250)))
  y  <- withr::with_seed(31, sin(latm[, 1] / 3200) + cos(latm[, 2] / 2500) +
                           stats::rnorm(nrow(latm), 0, 0.5))
  gi <- fnt_regions(latm, 1500L, 3L)
  K  <- max(gi)
  A  <- .pool_fit(y, gi, K, latm, use_fnn = FALSE)
  keep <- setdiff(names(A), "cen")
  maps <- list(feet = latm / 0.3048, `US feet` = latm * 3937 / 1200, km = latm / 1000,
               `origin moved by (0.1, 0.1)` = latm + 0.1,
               `origin moved by (5e5, 5.1e6)` = sweep(latm, 2L, c(5e5, 5.1e6), "+"))
  for (f in fnt_paths()) {
    expect_identical(.pool_fit(y, gi, K, latm, use_fnn = f), A, info = fnt_tag(f))
    for (nm in names(maps)) {
      B <- .pool_fit(y, gi, K, maps[[nm]], use_fnn = f)
      expect_identical(B[keep], A[keep], info = sprintf("%s, %s", nm, fnt_tag(f)))
    }
  }
})

test_that("summarize_by_cell() gives one layer the same pooled values in metres and in US feet", {
  g <- expand.grid(i = 0:59, j = 0:59)
  z <- withr::with_seed(8, sin(g$i / 6) + cos(g$j / 5) + stats::rnorm(nrow(g), 0, 0.5))
  pm <- sf::st_as_sf(data.frame(x = 4.5e5 + g$i * 250, y = 4.55e6 + g$j * 250,
                                cell_id = (g$j %/% 2) * 30 + g$i %/% 2 + 1, z = z),
                     coords = c("x", "y"), crs = 26915)
  pf <- sf::st_transform(pm, sf::st_crs("+proj=utm +zone=15 +datum=NAD83 +units=us-ft +no_defs"))
  pooled <- function(layer) sf::st_drop_geometry(fnt_quiet(suppressWarnings(
    summarize_by_cell(layer, "z", id_col = "cell_id", estimate = "pooled"))))
  a <- pooled(pm); b <- pooled(pf)
  cols <- grep("^\\.\\.pooled", names(a), value = TRUE)
  expect_true("..pooled_resp_z" %in% cols)
  expect_identical(a$cell_id, b$cell_id)
  # Column by column: the record of the fit names the CRS, which differs.
  for (cc in cols) expect_identical(b[[cc]], a[[cc]], info = cc)
})

test_that("squared distances within one part in ten million are a tie, which goes to the lower x", {
  near <- function(cen, f) .region_neighbours(cen, 1L, f)[1L, 1L]
  for (f in fnt_paths()) {
    tag <- fnt_tag(f)
    # From (0, 0): (1, 0) and (0, 1 + 1e-9) tie, and the lower x is taken.
    expect_identical(near(rbind(c(0, 0), c(1, 0), c(0, 1 + 1e-9), c(5, 5)), f), 3L, info = tag)
    # One part in 5e7 apart, beyond the margin of one part in 1e9 that the
    # kd-tree's answer is checked by.
    expect_identical(near(rbind(c(0, 0), c(1, 0), c(0, 1 + 1e-8), c(5, 5)), f), 3L, info = tag)
    # Further by more than the tolerance: the nearer is taken.
    expect_identical(near(rbind(c(0, 0), c(1, 0), c(0, 1 + 1e-6), c(5, 5)), f), 2L, info = tag)
    # The nearest other location, on the locations in coordinate order.
    u <- rbind(c(0, 0), c(0, 1 + 1e-9), c(1, 0), c(5, 5))
    expect_identical(.nearest_location_search(u, 1L, f), 2L, info = tag)
    u <- rbind(c(0, 0), c(0, 1 + 1e-6), c(1, 0), c(5, 5))
    expect_identical(.nearest_location_search(u, 1L, f), 3L, info = tag)
    # A row alone at its location reads the nearest other location.
    xy <- rbind(c(0, 0), c(1, 0), c(0, 1 + 1e-9), c(5, 5))
    expect_identical(.nearest_other(xy, 1:4, 1L, f), 3L, info = tag)
  }
})

test_that("x values within the tolerance of the distance count as equal, and the lower y is taken", {
  # From (0, 0), (0, 1) and (1e-12, -1) tie; their x agree to within 1e-7 of
  # the distance, so the lower y decides.  By x alone (0, 1) would be taken.
  cen <- rbind(c(0, 0), c(0, 1), c(1e-12, -1), c(5, 5))
  for (f in fnt_paths()) {
    expect_identical(.region_neighbours(cen, 1L, f)[1L, 1L], 3L, info = fnt_tag(f))
    expect_identical(.nearest_other(cen, 1:4, 1L, f), 3L, info = fnt_tag(f))
    u <- cen[order(cen[, 1], cen[, 2]), ]
    expect_identical(u[.nearest_location_search(u, which(u[, 2] == 0 & u[, 1] == 0), f), ],
                     c(1e-12, -1), info = fnt_tag(f))
  }
  # x values that differ by more than that: the lower x is taken.
  cen2 <- rbind(c(0, 0), c(0, 1), c(1e-6, -1), c(5, 5))
  for (f in fnt_paths())
    expect_identical(.region_neighbours(cen2, 1L, f)[1L, 1L], 2L, info = fnt_tag(f))
})

test_that("the kd-tree's search and the exact search pick by one rule, and the partition keeps exact ties", {
  # A 0.7 lattice from (0.1, 0.3), moved to UTM-sized coordinates and put
  # in feet: no tie in it is exact.
  g  <- unname(as.matrix(expand.grid(seq(0.1, by = 0.7, length.out = 25),
                                     seq(0.3, by = 0.7, length.out = 20))))
  u  <- sweep(g, 2L, c(5e5, 5.1e6), "+") / 0.3048
  u  <- u[order(u[, 1], u[, 2]), ]
  rows <- seq_len(nrow(u))
  exact <- .knn_exact_rows(u[, 1], u[, 2], rows, 8L, tol = .TIE_TOL)
  # The first neighbour of an interior location is the one to its left.
  inner <- which(u[, 1] > min(u[, 1]))
  expect_true(all(u[exact[inner, 1], 1] < u[inner, 1]))
  skip_if_not_installed("FNN")
  expect_identical(.knn_tied_rows(u, rows, 8L, tol = .TIE_TOL), exact)
  expect_identical(.knn_tied_rows(u, rows, 1L, tol = .TIE_TOL),
                   .knn_exact_rows(u[, 1], u[, 2], rows, 1L, tol = .TIE_TOL))
  # With no tolerance, as the partition's searches call them, only exact ties
  # tie: from (0, 0), (1, 0) is nearer than (0, 1 + 1e-9).
  v <- rbind(c(0, 0), c(0, 1 + 1e-9), c(1, 0), c(5, 5))
  expect_identical(.knn_exact_rows(v[, 1], v[, 2], 1L, 1L), matrix(3L))
  expect_identical(.knn_tied_rows(v, 1L, 1L), matrix(3L))
  for (f in fnt_paths())
    expect_identical(.knn_locations(v, 1L, use_fnn = f)$idx[1L, 1L], 3L, info = fnt_tag(f))
})


# ---- Regions whose centroids coincide, and the parts of the rule ------------------

# Each region's neighbours by the rule, one at a time, the last tie to the
# lower region number.
fnt_neighbours <- function(cen, k = 8L) {
  K <- nrow(cen); k <- min(as.integer(k), K - 1L)
  picks <- vapply(seq_len(K), function(j) {
    d2 <- (cen[, 1] - cen[j, 1])^2 + (cen[, 2] - cen[j, 2])^2
    left <- seq_len(K)[-j]; out <- integer(k)
    for (i in seq_len(k)) {
      m <- min(d2[left])
      s <- left[d2[left] <= m * (1 + 1e-7)]
      s <- s[cen[s, 1] <= min(cen[s, 1]) + 1e-7 * sqrt(m)]
      s <- s[cen[s, 2] <= min(cen[s, 2]) + 1e-7 * sqrt(m)]
      out[i] <- min(s)
      left <- left[left != out[i]]
    }
    out
  }, integer(k))
  matrix(as.vector(picks), K, k, byrow = TRUE)
}
fnt_sets_moved <- function(A, B)
  sum(vapply(seq_len(nrow(A)), function(i) !setequal(A[i, ], B[i, ]), TRUE))

test_that("of regions still tied the lower number is taken, not the lower place in coordinate order", {
  # From (0, 0) regions 2 and 3 are at one squared distance, with one x and
  # y within the tolerance of the distance.  Region 3 comes first in
  # coordinate order; region 2 is taken.  Numbered the other way round, the
  # other one is.
  cen <- rbind(c(0, 0), c(1, 1e-13), c(1, 0), c(5, 5))
  for (f in fnt_paths()) {
    expect_identical(.region_neighbours(cen, 1L, f)[1L, 1L], 2L, info = fnt_tag(f))
    expect_identical(.region_neighbours(cen[c(1, 3, 2, 4), ], 1L, f)[1L, 1L], 2L,
                     info = fnt_tag(f))
  }
  # The step itself: the lowest `num`, the position by default.
  d2 <- matrix(c(1, 1), 1L); cand <- matrix(1:2, 1L)
  expect_identical(.knn_pick_rows(d2, cand, 1L, .TIE_TOL, c(1, 1), c(0, 1e-13)), matrix(1L))
  expect_identical(.knn_pick_rows(d2, cand, 1L, .TIE_TOL, c(1, 1), c(0, 1e-13), num = 2:1),
                   matrix(2L))
})

test_that("regions whose centroids coincide to rounding get the rule's neighbours on both paths", {
  # Two pairs and a triple of centroids a few units in the last place apart,
  # among 60 others; in the first pair and the triple the lower number has
  # the higher coordinates.
  base <- withr::with_seed(15, cbind(stats::runif(60, 1e5, 2e5), stats::runif(60, 1e5, 2e5)))
  ulp <- function(v, n) v * (1 + n * .Machine$double.eps)
  for (p in list(c(3, 41), c(17, 8), c(25, 26, 50))) {
    base[p[-1L], ] <- base[rep(p[1L], length(p) - 1L), ]
    base[p[1L], ] <- ulp(base[p[1L], ], 2)            # the first one further up
  }
  base[c(26, 50), 2] <- ulp(base[c(26, 50), 2], c(-1, 3))
  expect_false(anyDuplicated(base) > 0L)
  for (f in fnt_paths())
    for (k in c(8L, 1L))
      expect_identical(.region_neighbours(base, k, f), fnt_neighbours(base, k),
                       info = sprintf("%s, %d neighbours", fnt_tag(f), k))
})

test_that("two cells whose centroids coincide keep their neighbours and values in another unit or origin", {
  # A 12 x 12 lattice of one-row cells at 250 m and a cell of seven rows on
  # one of its points.  In whole metres the two centroids coincide; in feet
  # or with the origin moved by 0.1 m rounding parts them.
  lat <- unname(as.matrix(expand.grid(0:11 * 250, 0:11 * 250))) + 1e5
  xy <- rbind(lat, matrix(rep(lat[40L, ], each = 7L), ncol = 2))
  gi <- c(seq_len(nrow(lat)), rep(nrow(lat) + 1L, 7L))
  K  <- max(gi)
  y  <- withr::with_seed(47, stats::rnorm(nrow(xy)))
  A  <- .pool_fit(y, gi, K, xy, use_fnn = FALSE)
  expect_identical(unname(A$cen[40L, ]), unname(A$cen[K, ]))
  keep <- setdiff(names(A), "cen")
  maps <- list(feet = xy / 0.3048, `US feet` = xy * 3937 / 1200,
               `origin moved by 0.1` = xy + 0.1)
  B <- .pool_fit(y, gi, K, maps[["origin moved by 0.1"]], use_fnn = FALSE)
  expect_false(identical(unname(B$cen[40L, ]), unname(B$cen[K, ])))
  for (f in fnt_paths()) for (nm in names(maps)) {
    B <- .pool_fit(y, gi, K, maps[[nm]], use_fnn = f)
    expect_identical(B[keep], A[keep], info = sprintf("%s, %s", nm, fnt_tag(f)))
  }
})

test_that("summarize_by_cell() keeps the pooled values of cells that share a centroid in another unit or origin", {
  pooled <- function(layer) sf::st_drop_geometry(fnt_quiet(suppressWarnings(
    summarize_by_cell(layer, "z", id_col = "cell_id", estimate = "pooled"))))
  same <- function(a, b, tag) {
    cols <- grep("^\\.\\.pooled", names(a), value = TRUE)
    expect_true("..pooled_resp_z" %in% cols, info = tag)
    expect_identical(b$cell_id, a$cell_id, info = tag)
    for (cc in cols) expect_identical(b[[cc]], a[[cc]], info = paste(tag, cc))
  }
  # Whether the centroids of two cells part in a layer.
  parted <- function(layer, i, j) {
    xy <- unname(sf::st_coordinates(layer)); gi <- layer$cell_id
    cen <- .pool_geometry(gi, max(gi), xy, tabulate(gi), FALSE)$cen
    !identical(unname(cen[i, ]), unname(cen[j, ]))
  }
  usft <- sf::st_crs("+proj=utm +zone=15 +datum=NAD83 +units=us-ft +no_defs")
  # UTM zone 15 with its false easting and northing moved by 0.1 m.
  moved <- sf::st_crs(paste("+proj=tmerc +lat_0=0 +lon_0=-93 +k=0.9996 +x_0=500000.1",
                            "+y_0=0.1 +datum=NAD83 +units=m +no_defs"))
  # The layout of test-core-pool-edge.R at UTM-sized coordinates: a grid of
  # four-row cells and two cross-shaped cells centred on one point.
  g  <- expand.grid(cx = seq(100, 900, by = 200), cy = seq(100, 700, by = 200))
  x  <- c(rep(g$cx, each = 4L) + c(-20, 20, -20, 20), 240, 260, 250, 250)
  y  <- c(rep(g$cy, each = 4L) + c(-20, -20, 20, 20), 250, 250, 230, 270)
  cell <- c(rep(seq_len(nrow(g)), each = 4L), 21L, 21L, 22L, 22L)
  z  <- withr::with_seed(3007, sin(x / 300) + cos(y / 250) + stats::rnorm(length(x), 0, 0.3))
  pm <- sf::st_as_sf(data.frame(x = 4.5e5 + x, y = 4.55e6 + y, cell_id = cell, z = z),
                     coords = c("x", "y"), crs = 26915)
  expect_false(parted(pm, 21L, 22L))
  a <- pooled(pm)
  for (crs in list(usft, moved)) {
    b <- sf::st_transform(pm, crs)
    same(a, pooled(b), crs$input)
  }
  expect_true(parted(sf::st_transform(pm, moved), 21L, 22L))
  # A 30 x 30 lattice of one-row cells at 250 m, and a cross-shaped cell of
  # four rows centred on each of 20 of its points.
  gl <- expand.grid(i = 0:29, j = 0:29)
  pts <- withr::with_seed(5, sample(nrow(gl), 20L))
  x <- 250 * gl$i; y <- 250 * gl$j; cell <- seq_len(nrow(gl))
  for (t in seq_along(pts)) {
    p <- pts[t]
    x <- c(x, x[p] - 50, x[p] + 50, x[p], x[p]); y <- c(y, y[p], y[p], y[p] - 50, y[p] + 50)
    cell <- c(cell, rep(nrow(gl) + t, 4L))
  }
  z  <- withr::with_seed(6, sin(x / 1500) + cos(y / 1200) + stats::rnorm(length(x), 0, 0.5))
  pm <- sf::st_as_sf(data.frame(x = 4.5e5 + x, y = 4.55e6 + y, cell_id = cell, z = z),
                     coords = c("x", "y"), crs = 26915)
  pf <- sf::st_transform(pm, usft)
  same(pooled(pm), pooled(pf), "lattice with cross-shaped cells, US feet")
  expect_true(any(vapply(seq_along(pts), function(t) parted(pf, pts[t], nrow(gl) + t), NA)))
})

test_that("up to nine cells on one point keep their neighbours in another unit or origin", {
  # Cells of 1, 2, ... rows, every row on one point, beside a 10 x 10 lattice
  # of one-row cells: in whole metres their centroids are that point, in feet
  # or with the origin moved some part by rounding.  With nine cells on the
  # point each takes all eight others, in whatever order; the values agree
  # to rounding.
  lat <- unname(as.matrix(expand.grid(0:9 * 250, 0:9 * 250))) + 1e5
  for (case in list(list("nine cells between lattice points", c(1e5 + 1125, 1e5 + 1125), 9L),
                    list("eight cells on a lattice point", c(1e5 + 1000, 1e5 + 1000), 8L))) {
    cc <- case[[3L]]
    xy <- rbind(lat, matrix(rep(case[[2L]], each = sum(seq_len(cc))), ncol = 2))
    gi <- c(seq_len(nrow(lat)), nrow(lat) + rep(seq_len(cc), seq_len(cc)))
    K  <- max(gi)
    y  <- withr::with_seed(50 + cc, stats::rnorm(nrow(xy)))
    A  <- .pool_fit(y, gi, K, xy, use_fnn = FALSE)
    crowd <- nrow(lat) + seq_len(cc)
    expect_identical(nrow(unique(A$cen[crowd, ])), 1L, info = case[[1L]])
    parted <- 0L
    for (f in fnt_paths())
      for (xy2 in list(xy / 0.3048, xy * 3937 / 1200, xy + 0.1)) {
        B <- .pool_fit(y, gi, K, xy2, use_fnn = f)
        parted <- parted + (nrow(unique(B$cen[crowd, ])) > 1L)
        expect_identical(fnt_sets_moved(A$nb, B$nb), 0L, info = case[[1L]])
        expect_equal(B$value, A$value, tolerance = 1e-12, info = case[[1L]])
      }
    expect_gt(parted, 0L)
  }
})

test_that("the tolerance on x and y is a fraction of the distance, at any scale", {
  # From (0, 0), (0, s) and (1e-8 s, -s) tie, and their x agree to within
  # 1e-7 of the distance s, so the lower y is taken, at any s.
  for (s in c(1e-3, 1, 1e3)) for (f in fnt_paths()) {
    cen <- rbind(c(0, 0), c(0, s), c(1e-8 * s, -s), c(5 * s, 5 * s))
    tag <- sprintf("s = %g, %s", s, fnt_tag(f))
    expect_identical(.region_neighbours(cen, 1L, f)[1L, 1L], 3L, info = tag)
    u <- cen[order(cen[, 1], cen[, 2]), ]
    expect_identical(u[.nearest_location_search(u, which(u[, 1] == 0 & u[, 2] == 0), f), ],
                     c(1e-8 * s, -s), info = tag)
  }
})

test_that("a kd-tree shortlist that ends inside the tolerance is made longer", {
  # From (0, 0): A at squared distance 1, B at 1 + 2e-8 with the lowest x of
  # the first four, C at 1 + 3e-8, D at 1 + 4e-8, and E at 1 + 5e-8 with the
  # lowest x of all.  All five tie and E is taken.  The kd-tree's first
  # shortlist holds A to D; its furthest, D, is beyond the pick B by more
  # than one part in 1e9 and by less than the tolerance.
  at <- function(r2, deg) sqrt(r2) * c(cos(deg * pi / 180), sin(deg * pi / 180))
  pts <- rbind(c(0, 0), at(1, 10), at(1 + 2e-8, 60), at(1 + 3e-8, 20), at(1 + 4e-8, 30),
               at(1 + 5e-8, 170), c(9, 9), c(-9, 9), c(9, -9), c(-9, -9))
  o <- order(pts[, 1], pts[, 2]); u <- pts[o, ]
  q <- which(o == 1L); e <- which(o == 6L)
  for (f in fnt_paths()) {
    expect_identical(.nearest_location_search(u, q, f), e, info = fnt_tag(f))
    expect_identical(.region_neighbours(pts, 1L, f)[1L, 1L], 6L, info = fnt_tag(f))
  }
})

test_that("squared distances two parts in ten million apart are not a tie", {
  for (f in fnt_paths())
    expect_identical(.region_neighbours(rbind(c(0, 0), c(1, 0), c(0, 1 + 1e-7), c(5, 5)), 1L, f)[1L, 1L],
                     2L, info = fnt_tag(f))
})

test_that("without FNN the short-lag pairs never call FNN", {
  skip_if_not_installed("FNN")
  L <- fnt_lattice()
  want <- .sac_short_lag(L$xy, L$z)
  got <- local({
    fnt_without_fnn()
    testthat::local_mocked_bindings(get.knn = function(...) stop("FNN called"),
                                    get.knnx = function(...) stop("FNN called"),
                                    .package = "FNN")
    .sac_short_lag(L$xy, L$z)
  })
  expect_identical(got, want)
})
