# tests/testthat/test-fix-partition-ties.R
# ---------------------------------------------------------------------------
# On gridded locations many distances tie exactly, and the last bits of the
# coordinates, which move with the unit and the origin, settled those ties
# in the density partition: its neighbour searches and its nearest-centre
# search counted exact ties alone, the first k-means++ centre was drawn from
# weights equal but for their last bits, the locations were put in order by
# first coordinates equal but for their last bits, and of two starts with
# sums equal to rounding the rounding chose.  So the regions
# resolution_profile() draws, and its rss and loo, changed with the unit or
# the origin.  The partition now takes the pooling's rule (.TIE_TOL), and
# the assignment's edge rule reaches as far as the tie between two centres.
# ---------------------------------------------------------------------------

.pt_seed <- function(seed, env = parent.frame())
  withr::local_seed(seed, .local_envir = env, .rng_kind = "Mersenne-Twister",
                    .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")

.pt_quiet <- function(expr) suppressWarnings(suppressMessages(
  logger::with_log_threshold(expr, threshold = logger::FATAL, namespace = "spatialkit",
                             index = 2)))

# Rows in another cell, in the worse direction: 0 when the two partitions are
# the same whatever their labels.
.pt_moved <- function(a, b) {
  t <- table(a, b)
  length(a) - min(sum(apply(t, 1L, max)), sum(apply(t, 2L, max)))
}

.pt_crs_m    <- function() sf::st_crs(26915)
.pt_crs_usft <- function() sf::st_crs("+proj=utm +zone=15 +datum=NAD83 +units=us-ft +no_defs")
.pt_crs_m01  <- function() sf::st_crs(paste(
  "+proj=tmerc +lat_0=0 +lon_0=-93 +k=0.9996 +x_0=500000.1 +y_0=0.1 +datum=NAD83",
  "+units=m +no_defs"))

# A layer in EPSG:26915 with a smooth response and noise.
.pt_layer <- function(x, y) {
  z <- local({
    .pt_seed(1)
    sin((x - min(x)) / 1500) + cos((y - min(y)) / 1250) + stats::rnorm(length(x), 0, 0.5)
  })
  sf::st_as_sf(data.frame(x = x, y = y, z = z, row = seq_along(x)), coords = c("x", "y"),
               crs = .pt_crs_m())
}
.pt_lattice <- function() {
  g <- expand.grid(i = 0:39, j = 0:29)
  .pt_layer(4.5e5 + g$i * 250, 4.55e6 + g$j * 250)
}
# The layer in metres, in UTM 15N in US feet, and in UTM 15N with the false
# easting and northing moved by 0.1 m.
.pt_frames <- function(layer)
  list(m = layer, usft = sf::st_transform(layer, .pt_crs_usft()),
       m01 = sf::st_transform(layer, .pt_crs_m01()))

# Lattices on which most rows tie, as coordinates in metres.
.pt_lattices <- function() list(
  square = { g <- expand.grid(i = 0:39, j = 0:29); cbind(4.5e5 + g$i * 250, 4.55e6 + g$j * 250) },
  triangular = { g <- expand.grid(i = 0:19, j = 0:19)
                 cbind(4e5 + (g$i + (g$j %% 2) / 2) * 200, 4.6e6 + g$j * 200 * sqrt(3) / 2) },
  repeated = { g <- expand.grid(i = 0:19, j = 0:19)
               b <- cbind(4.2e5 + g$i * 300, 4.7e6 + g$j * 300)
               rbind(b, b[local({ .pt_seed(5); sample.int(400L, 120L) }), ]) })


# ---- The item's test ----------------------------------------------------------------

test_that("a lattice gives the same regions, rss and loo in metres, in feet and with the origin moved", {
  frames <- .pt_frames(.pt_lattice())
  xy <- lapply(frames, function(l) unname(sf::st_coordinates(l)))
  # The frames differ in their last bits: after the change of origin the 40
  # columns of the lattice have hundreds of distinct x.
  expect_identical(length(unique(xy$m[, 1L])), 40L)
  expect_gt(length(unique(xy$m01[, 1L])), 40L)
  prof <- lapply(frames, function(l) .pt_quiet(resolution_profile(l, "z")))
  for (r in c("usft", "m01")) {
    expect_identical(prof[[r]]$levels, prof$m$levels, info = r)
    expect_lt(max(abs(prof[[r]]$loo - prof$m$loo)), 1e-12)
    expect_lt(max(abs(prof[[r]]$rss - prof$m$rss) / prof$m$rss), 1e-12)
  }
  # The cells of the level, row by row: each row's nearest centre in its
  # frame's coordinates.
  cells <- lapply(names(frames), function(r) {
    cen <- attr(prof[[r]], "centres")[[1L]]
    .nearest_centre(xy[[r]][, 1L], xy[[r]][, 2L], cen[, 1L], cen[, 2L])
  })
  expect_identical(cells[[2L]], cells[[1L]])
  expect_identical(cells[[3L]], cells[[1L]])
})

test_that("a profile that draws a subsample draws the same locations in every frame", {
  # Without a response the automatic ladder is fitted to a subsample of
  # `sample_n` rows, drawn from the rows in coordinate order.  After the change
  # of origin the first coordinates of a column differ in their last bits, so
  # an exact order put the column's rows in an order those bits set and the
  # draw took other locations: other cells, and another wss.
  g <- expand.grid(i = 0:19, j = 0:19)
  frames <- .pt_frames(.pt_layer(4.5e5 + g$i * 250, 4.55e6 + g$j * 250))
  xy <- lapply(frames, function(l) unname(sf::st_coordinates(l)))
  prof <- lapply(frames, function(l) .pt_quiet(resolution_profile(l, sample_n = 300L)))
  expect_identical(attr(prof$m, "bounds")$n_sample, 300L)
  for (r in c("usft", "m01")) expect_identical(prof[[r]]$levels, prof$m$levels, info = r)
  # Level by level, among the levels all three frames have: the wss as a
  # fraction of the total sum of squares, and the cells row by row.
  common <- Reduce(intersect, lapply(prof, function(p) p$levels))
  expect_gt(length(common), 5L)
  for (k in common) {
    at <- vapply(prof, function(p) match(k, p$levels), integer(1L))
    tss <- vapply(xy, function(a) sum(sweep(a, 2L, colMeans(a))^2), numeric(1L))
    w <- vapply(names(prof), function(r) prof[[r]]$wss[at[[r]]], numeric(1L)) / tss
    cells <- lapply(names(frames), function(r) {
      cen <- attr(prof[[r]], "centres")[[at[[r]]]]
      .nearest_centre(xy[[r]][, 1L], xy[[r]][, 2L], cen[, 1L], cen[, 2L])
    })
    for (j in 2:3) {
      r <- names(frames)[j]
      expect_lt(abs(w[[r]] - w[["m"]]), 1e-9, label = paste("wss at", k, r))
      expect_identical(.pt_moved(cells[[1L]], cells[[j]]), 0L, info = paste(k, r))
    }
  }
})

test_that("the help bounds the rounding the cells survive by the spacing of the locations", {
  # The tie rule is relative to the squared distances, so the rounding a tie
  # survives is set by the distance between neighbouring locations (about a
  # forty-millionth of it), not by the extent.  The item's lattice in feet:
  # to 5 decimals the rounding is 6e-9 of the spacing and the cells are
  # those of the lattice in metres; to 4 decimals it is 6e-8, which is 1.2e-9
  # of the extent, and the cells move.
  g <- expand.grid(i = 0:39, j = 0:29)
  m <- cbind(4.5e5 + g$i * 250, 4.55e6 + g$j * 250)
  fit <- function(a) local({ .pt_seed(1); .density_kmeans(a, 25L) })$cluster
  base <- fit(m)
  ft5 <- round(m / 0.3048, 5); ft4 <- round(m / 0.3048, 4)
  expect_lt(max(abs(ft5 * 0.3048 - m)) / 250, 1e-8)
  expect_identical(.pt_moved(base, fit(ft5)), 0L)
  ex <- sqrt(sum(apply(m, 2L, function(v) diff(range(v)))^2))
  expect_lt(max(abs(ft4 * 0.3048 - m)) / ex, 1e-8)
  expect_gt(.pt_moved(base, fit(ft4)), 0L)
  # And the help says so: the source page from the source tree, the installed
  # page under R CMD check, which runs the tests where there is no man/.
  file <- test_path("..", "..", "man", "resolution_profile.Rd")
  rd <- if (file.exists(file)) tools::parse_Rd(file)
        else tools::Rd_db("spatialkit")[["resolution_profile.Rd"]]
  expect_s3_class(rd, "Rd")
  txt <- gsub("[[:space:]]+", " ", paste(unlist(rd), collapse = ""))
  expect_match(txt, "forty-millionth of the distance between neighbouring locations",
               fixed = TRUE)
  expect_no_match(txt, "ten-millionth of the extent", fixed = TRUE)
})

test_that("a uniform layer gives the profile it gave, in every frame", {
  # Its distances do not tie, so nothing in it changes: the figures are
  # those of the profile before the tie rule.
  layer <- local({
    .pt_seed(3)
    .pt_layer(4.5e5 + stats::runif(1200, 0, 1e4), 4.55e6 + stats::runif(1200, 0, 7500))
  })
  prof <- lapply(.pt_frames(layer), function(l) .pt_quiet(resolution_profile(l, "z")))
  for (r in c("usft", "m01")) {
    expect_identical(prof[[r]]$loo, prof$m$loo, info = r)
    expect_identical(prof[[r]]$rss, prof$m$rss, info = r)
  }
  expect_equal(prof$m$loo, 0.3100867, tolerance = 1e-6)
  expect_equal(prof$m$rss, 271.3873, tolerance = 1e-6)
})

test_that("the density partition of a lattice is the same in five frames, in any row order and without FNN", {
  for (nm in names(.pt_lattices())) {
    m <- .pt_lattices()[[nm]]
    base <- sf::st_as_sf(data.frame(x = m[, 1L], y = m[, 2L]), coords = c("x", "y"),
                         crs = .pt_crs_m())
    xy <- c(lapply(.pt_frames(base), function(l) unname(sf::st_coordinates(l))),
            list(scaled = m / 0.3048, shifted = sweep(m, 2L, c(0.1, -0.37), "+")))
    perm <- local({ .pt_seed(11); sample.int(nrow(m)) })
    for (k in c(7L, 25L)) for (init in c("kmeans++", "random")) {
      info <- paste(nm, k, init)
      fit <- lapply(xy, function(a) { .pt_seed(2); .density_kmeans(a, k, init = init) })
      for (r in names(xy)[-1L])
        expect_identical(.pt_moved(fit$m$cluster, fit[[r]]$cluster), 0L, info = paste(info, r))
      pf <- local({ .pt_seed(2); .density_kmeans(m[perm, ], k, init = init) })
      expect_identical(.pt_moved(fit$m$cluster[perm], pf$cluster), 0L, info = paste(info, "rows"))
      nf <- local({ .pt_seed(2); .density_kmeans(m, k, init = init, use_fnn = FALSE) })
      expect_identical(.pt_moved(fit$m$cluster, nf$cluster), 0L, info = paste(info, "no FNN"))
    }
  }
})


# ---- The pieces --------------------------------------------------------------------

test_that("a tie between centres goes to the lower number, by a rule of two passes", {
  # From (0, 0): centre 2 at 10 (1 + 2e-10), centre 3 at 10.  A tie, to 2.
  cx <- c(100, 10 * (1 + 2e-10), 0); cy <- c(100, 0, 10)
  expect_identical(.nearest_exact(0, 0, cx, cy), 2L)
  # At 10 (1 + 2e-7) the squared distances differ by 4e-7: the nearer, 3.
  cx[2L] <- 10 * (1 + 2e-7)
  expect_identical(.nearest_exact(0, 0, cx, cy), 3L)
  # Squared distances 1 + 1.5e-7, 1 + 0.8e-7 and 1: centres 2 and 3 tie with
  # the least, centre 1 does not, so centre 2.  A single pass that keeps a
  # centre only when it is nearer by more than the tolerance keeps 3.
  r <- sqrt(c(1 + 1.5e-7, 1 + 0.8e-7, 1))
  expect_identical(.nearest_exact(0, 0, r, c(0, 0, 0)), 2L)
  # The kd-tree path gives the same, among 40 centres.
  far <- local({ .pt_seed(4); stats::runif(37, 50, 90) })
  for (f in c(TRUE, FALSE))
    expect_identical(.nearest_centre(0, 0, c(r, far), c(0, 0, 0, far), use_fnn = f), 2L,
                     info = paste("use_fnn", f))
})

test_that("the neighbours of a lattice's locations are the same with and without FNN, whatever the last bits", {
  skip_if_not_installed("FNN")
  for (nm in names(.pt_lattices())) {
    m <- .pt_lattices()[[nm]]
    noisy <- m + local({ .pt_seed(8); matrix(stats::runif(length(m), -1, 1) * 1e-9, ncol = 2L) })
    for (a in list(m, noisy)) {
      u <- .locations_by_column(.partition_locations(a))$u
      for (k in c(4L, 8L)) {
        f <- .knn_locations(u, k, use_fnn = TRUE)
        e <- .knn_locations(u, k, use_fnn = FALSE)
        expect_identical(t(apply(f$idx, 1L, sort)), t(apply(e$idx, 1L, sort)), info = nm)
        expect_identical(f$d, e$d, info = nm)
      }
    }
    # The neighbour sets, as rows of the layer, do not move with the noise.
    sets <- lapply(list(m, noisy), function(a) {
      loc <- .locations_by_column(.partition_locations(a))
      row_of <- integer(nrow(loc$u)); row_of[loc$gid] <- loc$order
      idx <- .knn_locations(loc$u, 8L)$idx
      out <- t(apply(matrix(row_of[idx], nrow(idx)), 1L, sort))
      out[order(row_of), , drop = FALSE]
    })
    if (nm != "repeated") expect_identical(sets[[2L]], sets[[1L]], info = nm)
  }
})

test_that("weights equal but for their last bits draw the same first k-means++ centre", {
  n <- 300L
  x <- seq_len(n) / n; y <- rep(0, n)
  w1 <- rep(1, n)
  w2 <- w1 * (1 + 1e-12 * local({ .pt_seed(9); stats::runif(n) }))
  for (seed in 1:5) {
    a <- local({ .pt_seed(seed); .wkmeanspp(x, y, w1, 1L) })
    b <- local({ .pt_seed(seed); .wkmeanspp(x, y, w2, 1L) })
    expect_identical(b, a, info = paste("seed", seed))
  }
  # Weights no two of which are within one part in ten million draw what
  # sample.int() draws on them.
  w3 <- local({ .pt_seed(10); stats::rexp(n) })
  expect_gt(min(diff(sort(w3)) / sort(w3)[-n]), 1e-7)
  for (seed in 1:5) {
    a <- local({ .pt_seed(seed); .wkmeanspp(x, y, w3, 1L) })
    b <- local({ .pt_seed(seed); sample.int(n, 1L, prob = w3) })
    expect_identical(a, b, info = paste("seed", seed))
  }
})

test_that("first coordinates equal but for their last bits are one column when the locations are put in order", {
  g <- expand.grid(i = 0:9, j = 0:7)
  m <- cbind(4.5e5 + g$i * 250, 4.55e6 + g$j * 250)
  noisy <- m
  noisy[, 1L] <- noisy[, 1L] + local({ .pt_seed(12); stats::runif(nrow(m), -1, 1) }) * 1e-9
  expect_gt(length(unique(noisy[, 1L])), 10L)
  rows <- lapply(list(m, noisy), function(a) {
    loc <- .locations_by_column(.partition_locations(a))
    row_of <- integer(nrow(loc$u)); row_of[loc$gid] <- loc$order
    row_of
  })
  expect_identical(rows[[2L]], rows[[1L]])
  # So the partition is the same.
  for (init in c("kmeans++", "random")) {
    a <- local({ .pt_seed(13); .density_kmeans(m, 6L, init = init) })
    b <- local({ .pt_seed(13); .density_kmeans(noisy, 6L, init = init) })
    expect_identical(b$cluster, a$cluster, info = init)
  }
  # Where no two first coordinates are that close, the order is the exact one.
  set_xy <- local({ .pt_seed(14); cbind(stats::runif(50), stats::runif(50)) })
  loc <- .partition_locations(set_xy)
  expect_identical(.locations_by_column(loc), loc)
})

test_that("of two starts whose sums agree to rounding the first is kept, in any unit", {
  # Random starts 1 and 3 on this lattice end in mirror images with sums equal
  # to rounding; which was kept rested on the last bits.
  g <- expand.grid(i = 0:39, j = 0:29)
  m <- cbind(4.5e5 + g$i * 250, 4.55e6 + g$j * 250)
  a <- local({ .pt_seed(3); .density_kmeans(m, 7L, init = "random") })
  b <- local({ .pt_seed(3); .density_kmeans(m / 0.3048, 7L, init = "random") })
  expect_identical(.pt_moved(a$cluster, b$cluster), 0L)
  expect_identical(b$iter, a$iter)
})


# ---- The assignment reaches as far as the tie -------------------------------------

test_that("a row two centres tie for is assigned where the profile scored it, beyond the margin of an edge", {
  # Two clusters 50 km apart in a transverse Mercator with coordinates near
  # 0, so the margin of an edge is a few micrometres, and one row P 1e-4 m
  # from the bisector of the two centres, on the side of centre 2.  Its
  # squared distances agree to about 1.5e-8: a tie, which the profile gives
  # to centre 1.  The polygons hold P in the cell of centre 2.
  crs <- sf::st_crs(paste("+proj=tmerc +lat_0=41 +lon_0=-93 +k=1 +x_0=0 +y_0=0",
                          "+datum=NAD83 +units=m +no_defs"))
  cl <- local({
    .pt_seed(1)
    rbind(cbind(stats::rnorm(200, -25000, 2000), stats::rnorm(200, 0, 2000)),
          cbind(stats::rnorm(200, 25000, 2000), stats::rnorm(200, 0, 2000)))
  })
  mk <- function(p) {
    xy <- rbind(cl, p)
    z <- local({ .pt_seed(2); stats::rnorm(nrow(xy)) })
    sf::st_as_sf(data.frame(x = xy[, 1L], y = xy[, 2L], z = z, row = seq_len(nrow(xy))),
                 coords = c("x", "y"), crs = crs)
  }
  run <- function(pts) {
    prof <- .pt_quiet(resolution_profile(pts, "z", partition = "density", levels = 2:4,
                                         seed = 123, min_cell_n = 5L, estimate = "pooled"))
    .pt_quiet(select_resolution(prof, criterion = "size", obs_per_cell = 200))
  }
  # Put P where it is 1e-4 m from the bisector of the centres it leads to.
  p <- c(0, 6000)
  for (it in 1:4) {
    cen <- run(mk(p))$count_rule$scored$centres
    mid <- colMeans(cen)
    dir <- (cen[2L, ] - cen[1L, ]) / sqrt(sum((cen[2L, ] - cen[1L, ])^2))
    p <- p - (sum((p - mid) * dir) - 1e-4) * dir
  }
  pts <- mk(p)
  sel <- run(pts)
  cen <- sel$count_rule$scored$centres
  d2 <- (p[1L] - cen[, 1L])^2 + (p[2L] - cen[, 2L])^2
  expect_lt(d2[2L], d2[1L])
  expect_lt((d2[1L] - d2[2L]) / d2[2L], 1e-7)
  expect_identical(.nearest_centre(p[1L], p[2L], cen[, 1L], cen[, 2L]), 1L)
  # P is farther from the edge than the margin.
  dist <- abs(sum((p - colMeans(cen)) * (cen[2L, ] - cen[1L, ]))) / sqrt(sum((cen[2L, ] - cen[1L, ])^2))
  expect_gt(dist, 10 * .edge_tol(sf::st_bbox(pts)))
  seeds <- .pt_quiet(get_voronoi_seeds(method = "kmeans", n = sel, sample_points = pts))
  tess  <- .pt_quiet(build_tessellation(seeds, method = "voronoi", quiet = TRUE))
  asg   <- .pt_quiet(assign_features_to_polygons(pts, tess$cells, polygon_id_col = "cell_id",
                                                 quiet = TRUE))
  xy <- unname(sf::st_coordinates(pts))
  scored <- tess$index[.nearest_centre(xy[, 1L], xy[, 2L], cen[, 1L], cen[, 2L])]
  placed <- asg$cell_id[match(pts$row, asg$row)]
  expect_identical(placed, scored)
  # Without the record the polygons put P in the other cell.
  bare <- tess$cells; attr(bare, "count_rule") <- NULL
  alone <- .pt_quiet(assign_features_to_polygons(pts, bare, polygon_id_col = "cell_id",
                                                 quiet = TRUE))
  expect_false(identical(alone$cell_id[match(nrow(pts), alone$row)], scored[nrow(pts)]))
  # And the summary's check reads match.
  out <- .pt_quiet(summarize_by_cell(asg, "z", id_col = "cell_id", cells_sf = tess$cells,
                                     estimate = "pooled"))
  expect_identical(attr(out, "count_rule")$check, list(estimate = "match", cells = "match"))
  # Cells drawn in the same projection in US feet: the margin's extension is
  # measured in the record's CRS and scaled to the join's by the ratio of the
  # centres' separations, so P is still placed where it was scored, with the
  # features in metres and in US feet.
  crs_ft <- sf::st_crs(paste("+proj=tmerc +lat_0=41 +lon_0=-93 +k=1 +x_0=0 +y_0=0",
                             "+datum=NAD83 +units=us-ft +no_defs"))
  tess_ft <- .pt_quiet(build_tessellation(seeds, method = "voronoi", quiet = TRUE,
                                          crs = crs_ft))
  expect_true(sf::st_crs(tess_ft$cells) == crs_ft)
  scored_ft <- tess_ft$index[.nearest_centre(xy[, 1L], xy[, 2L], cen[, 1L], cen[, 2L])]
  for (f in list(pts, sf::st_transform(pts, crs_ft))) {
    asg_ft <- .pt_quiet(assign_features_to_polygons(f, tess_ft$cells,
                                                    polygon_id_col = "cell_id", quiet = TRUE))
    expect_identical(asg_ft$cell_id[match(pts$row, asg_ft$row)], scored_ft)
    out_ft <- .pt_quiet(summarize_by_cell(asg_ft, "z", id_col = "cell_id",
                                          cells_sf = tess_ft$cells, estimate = "pooled"))
    expect_identical(attr(out_ft, "count_rule")$check,
                     list(estimate = "match", cells = "match"))
  }
})
