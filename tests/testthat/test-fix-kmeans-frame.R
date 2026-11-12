# tests/testthat/test-fix-kmeans-frame.R
# ---------------------------------------------------------------------------
# On gridded locations many points are equally far from two centres, and
# stats::kmeans() (Hartigan and Wong) settles each such tie by the last bits
# of the coordinates, which move with the unit and the origin.  So the
# regions of the kind "kmeans" (resolution_profile(), the seeding functions,
# the refit of select_resolution(), the sweep of determine_optimal_levels())
# changed with the unit or the origin.  stats::kmeans() is now handed the
# coordinates centred on the middle of their bounding box, divided by its
# diagonal and rounded to 2^-30 of it (.kmeans_frame()), and the centres and
# sums of squares are taken back in the caller's coordinates (.kmeans_back()).
# ---------------------------------------------------------------------------

.kf_seed <- function(seed, env = parent.frame())
  withr::local_seed(seed, .local_envir = env, .rng_kind = "Mersenne-Twister",
                    .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")

.kf_quiet <- function(expr) suppressWarnings(suppressMessages(
  logger::with_log_threshold(expr, threshold = logger::FATAL, namespace = "spatialkit",
                             index = 2)))

# Rows in another cell, in the worse direction: 0 when the two partitions are
# the same whatever their labels.
.kf_moved <- function(a, b) {
  t <- table(a, b)
  length(a) - min(sum(apply(t, 1L, max)), sum(apply(t, 2L, max)))
}

.kf_crs <- function(nm) sf::st_crs(switch(nm,
  m    = "EPSG:26915",
  usft = "+proj=utm +zone=15 +datum=NAD83 +units=us-ft +no_defs",
  ft   = "+proj=utm +zone=15 +datum=NAD83 +units=ft +no_defs",
  km   = "+proj=utm +zone=15 +datum=NAD83 +units=km +no_defs",
  m01  = paste("+proj=tmerc +lat_0=0 +lon_0=-93 +k=0.9996 +x_0=500000.1 +y_0=0.1",
               "+datum=NAD83 +units=m +no_defs"),
  m37  = paste("+proj=tmerc +lat_0=0 +lon_0=-93 +k=0.9996 +x_0=500137.37 +y_0=-2000.01",
               "+datum=NAD83 +units=m +no_defs")))

# A layer in EPSG:26915 with a smooth response and noise.
.kf_layer <- function(m) {
  x <- m[, 1L]; y <- m[, 2L]
  z <- local({
    .kf_seed(1)
    sin((x - min(x)) / 1500) + cos((y - min(y)) / 1250) + stats::rnorm(length(x), 0, 0.5)
  })
  sf::st_as_sf(data.frame(x = x, y = y, z = z), coords = c("x", "y"), crs = .kf_crs("m"))
}
.kf_grid <- function(nx, ny, sp, x0 = 4.5e5, y0 = 4.55e6) {
  g <- expand.grid(i = seq_len(nx) - 1L, j = seq_len(ny) - 1L)
  cbind(x0 + g$i * sp, y0 + g$j * sp)
}
# The layer in metres, in UTM 15N in US feet, and in UTM 15N with the false
# easting and northing moved by 0.1 m.
.kf_frames <- function(layer, frames = c("m", "usft", "m01"))
  stats::setNames(lapply(frames, function(f)
    if (f == "m") layer else sf::st_transform(layer, .kf_crs(f))), frames)
.kf_xy <- function(frames) lapply(frames, function(l) unname(sf::st_coordinates(l)[, 1:2]))
# Points of another frame back in metres.
.kf_in_m <- function(s) unname(sf::st_coordinates(sf::st_transform(s, .kf_crs("m")))[, 1:2])


# ---- The item's test ----------------------------------------------------------------

test_that("a lattice gives the same k-means regions, rss and loo in metres, in feet and with the origin moved", {
  frames <- .kf_frames(.kf_layer(.kf_grid(40, 30, 250)))
  xy <- .kf_xy(frames)
  expect_gt(length(unique(xy$m01[, 1L])), 40L)
  prof <- lapply(frames, function(l) .kf_quiet(resolution_profile(l, "z", partition = "kmeans")))
  for (r in c("usft", "m01")) {
    expect_identical(prof[[r]]$levels, prof$m$levels, info = r)
    expect_lt(max(abs(prof[[r]]$loo - prof$m$loo)), 1e-12)
    expect_lt(max(abs(prof[[r]]$rss - prof$m$rss) / prof$m$rss), 1e-12)
  }
  # The cells, row by row, numbered as the profile numbers them.
  cells <- lapply(names(frames), function(r) {
    cen <- attr(prof[[r]], "centres")[[1L]]
    .nearest_centre(xy[[r]][, 1L], xy[[r]][, 2L], cen[, 1L], cen[, 2L])
  })
  expect_identical(cells[[2L]], cells[[1L]])
  expect_identical(cells[[3L]], cells[[1L]])
})

test_that("the k-means seeds of a lattice are the same in metres, in feet and with the origin moved", {
  frames <- .kf_frames(.kf_layer(.kf_grid(40, 30, 250)))
  for (init in list(NULL, "kmeans++")) {
    info <- if (is.null(init)) "random starts" else init
    v <- lapply(frames, function(l)
      .kf_quiet(voronoi_seeds_kmeans(l, k = 40, partition = "kmeans", init = init)))
    g <- lapply(frames, function(l)
      .kf_quiet(get_voronoi_seeds(method = "kmeans", n = 40, sample_points = l,
                                  partition = "kmeans", init = init, set_seed = 3)))
    for (r in c("usft", "m01")) {
      expect_lt(max(abs(.kf_in_m(v[[r]]) - .kf_in_m(v$m))), 1e-6, label = paste(info, r))
      expect_lt(max(abs(.kf_in_m(g[[r]]) - .kf_in_m(g$m))), 1e-6, label = paste(info, r))
      expect_identical(attr(g[[r]], "kmeans")$cluster, attr(g$m, "kmeans")$cluster,
                       info = paste(info, r))
      expect_identical(attr(g[[r]], "kmeans")$size, attr(g$m, "kmeans")$size,
                       info = paste(info, r))
    }
  }
})

test_that(".kmeans_best() gives a lattice the same clusters in every frame at 5, 25 and 60 centres", {
  xy <- .kf_xy(.kf_frames(.kf_layer(.kf_grid(40, 30, 250))))
  for (k in c(5L, 25L, 60L)) for (s in 1:2) {
    cl <- lapply(xy, function(a) { .kf_seed(s); .kmeans_best(a, k)$km$cluster })
    expect_identical(cl$usft, cl$m, info = paste(k, s))
    expect_identical(cl$m01, cl$m, info = paste(k, s))
  }
})

test_that("lattices that span a kilometre or more get the same k-means cells in five other frames", {
  lattices <- list(.kf_grid(20, 15, 100, 3.2e5, 4.1e6), .kf_grid(25, 25, 250, 4.5e5, 4.1e6),
                   .kf_grid(60, 8, 30.48, 3.2e5, 4.55e6), .kf_grid(12, 12, 1000, 4.5e5, 4.55e6))
  for (i in seq_along(lattices)) {
    m <- lattices[[i]]
    expect_gte(sqrt(sum(apply(m, 2L, function(v) diff(range(v)))^2)), 1000)
    xy <- .kf_xy(.kf_frames(.kf_layer(m), c("m", "usft", "ft", "km", "m01", "m37")))
    for (k in c(5L, 25L)) {
      kp <- lapply(xy, function(a) { .kf_seed(1); .partition_centres(a, k, sort_rows = TRUE)$cluster })
      rn <- lapply(xy, function(a) {
        .kf_seed(1)
        suppressWarnings(.partition_centres(a, k, init = "random", nstart = 10L,
                                            iter.max = 100L, one_call = TRUE))$cluster
      })
      for (r in names(xy)[-1L]) {
        expect_identical(kp[[r]], kp$m, info = paste(i, k, r, "kmeans++"))
        expect_identical(rn[[r]], rn$m, info = paste(i, k, r, "random"))
      }
    }
  }
})

test_that("the help states how far the kind \"kmeans\" is the same in another unit or origin", {
  # The source page from the source tree, the installed page under R CMD check,
  # which runs the tests where there is no man/.
  txt <- function(f) {
    file <- test_path("..", "..", "man", f)
    rd <- if (file.exists(file)) tools::parse_Rd(file) else tools::Rd_db("spatialkit")[[f]]
    expect_s3_class(rd, "Rd")
    gsub("[[:space:]]+", " ", paste(unlist(rd), collapse = ""))
  }
  rp <- txt("resolution_profile.Rd")
  expect_match(rp, "rounded to 2^{-30} of it", fixed = TRUE)
  expect_match(rp, paste0("The cells are then the same in another unit or with another origin ",
                          "unless the change moves a coordinate across the middle of a rounding step"),
               fixed = TRUE)
  expect_match(rp, "does more often the smaller the extent", fixed = TRUE)
  for (f in c("get_voronoi_seeds.Rd", "voronoi_seeds_kmeans.Rd")) {
    s <- txt(f)
    expect_match(s, "rounded to 2^-30 of it, so that the seeds do not move with the unit or the origin",
                 fixed = TRUE, info = f)
    expect_match(s, "span less than a kilometre a change of unit or origin can still move them",
                 fixed = TRUE, info = f)
  }
})


# ---- The other places the kind runs ---------------------------------------------------

test_that("select_resolution() refits a k-means level of a lattice to the same seeds in every frame", {
  frames <- .kf_frames(.kf_layer(.kf_grid(40, 30, 250)))
  sel <- lapply(frames, function(l) .kf_quiet({
    p <- resolution_profile(l, "z", estimate = NULL, obs_per_cell = NULL, sac = NULL,
                            partition = "kmeans", sample_n = 600)
    s <- select_resolution(p, criterion = "size", obs_per_cell = 10)
    get_voronoi_seeds(n = s, method = "kmeans")
  }))
  expect_false(is.null(attr(sel$m, "refit")))
  for (r in c("usft", "m01")) {
    expect_lt(max(abs(.kf_in_m(sel[[r]]) - .kf_in_m(sel$m))), 1e-6, label = r)
    a <- attr(sel[[r]], "refit"); b <- attr(sel$m, "refit")
    expect_identical(a$iter, b$iter, info = r)
    expect_lt(abs(a$wss / a$wss_scored - b$wss / b$wss_scored), 1e-12)
  }
})

test_that("determine_optimal_levels() scores the same counts on a lattice in every frame", {
  # Twenty blocks of a 100 m lattice, 2 km apart: an elbow at 20, whose
  # neighbourhood the model-aware criterion scores.
  g <- expand.grid(i = 0:5, j = 0:5, bx = 0:4, by = 0:3)
  m <- cbind(4.5e5 + g$i * 100 + g$bx * 2000, 4.55e6 + g$j * 100 + g$by * 2000)
  b <- sf::st_as_sf(data.frame(x = m[, 1L], y = m[, 2L],
                               z = local({ .kf_seed(2); stats::rnorm(nrow(m)) })),
                    coords = c("x", "y"), crs = .kf_crs("m"))
  d <- lapply(.kf_frames(b), function(l)
    .kf_quiet(determine_optimal_levels(l, max_levels = 40, criterion = "morans_i",
                                       response_var = "z")))
  for (r in c("usft", "m01")) {
    expect_identical(as.integer(d[[r]]), as.integer(d$m), info = r)
    expect_identical(attr(d[[r]], "diagnostics")$moran_z, attr(d$m, "diagnostics")$moran_z,
                     info = r)
  }
})


# ---- What stays as it was --------------------------------------------------------------

test_that("a layer with no near ties gets the clusters stats::kmeans() gives it, and the same sums", {
  xy <- local({ .kf_seed(4); cbind(x = 4.5e5 + stats::runif(500, 0, 9000),
                                   y = 4.55e6 + stats::runif(500, 0, 6000)) })
  # The run of the seeding functions: one call from random starts.
  fit <- local({ .kf_seed(9); .partition_centres(xy, 12L, init = "random", nstart = 10L,
                                                 iter.max = 100L, one_call = TRUE) })
  km  <- local({ .kf_seed(9); stats::kmeans(xy, centers = 12L, nstart = 10L, iter.max = 100L) })
  expect_identical(fit$cluster, km$cluster)
  expect_identical(dimnames(fit$km$centers), dimnames(km$centers))
  expect_equal(fit$km$centers, km$centers, tolerance = 1e-12)
  expect_equal(fit$withinss, km$withinss, tolerance = 1e-10)
  expect_equal(fit$wss, km$tot.withinss, tolerance = 1e-12)
  expect_equal(fit$km$totss, km$totss, tolerance = 1e-12)
  # k-means++ starts: the starts earlier versions drew, on the raw coordinates.
  kb <- local({ .kf_seed(9); .kmeans_best(xy, 12L, nstart = 5L) })
  old <- local({
    .kf_seed(9)
    fits <- lapply(1:5, function(r) suppressWarnings(stats::kmeans(xy, centers = .kmeanspp_centers(xy, 12L),
                                                                   iter.max = 50L, nstart = 1L)))
    fits[[which.min(vapply(fits, function(f) f$tot.withinss, numeric(1L)))]]
  })
  expect_identical(.kf_moved(kb$km$cluster, old$cluster), 0L)
  expect_equal(kb$wss, old$tot.withinss, tolerance = 1e-12)
})

test_that("stats::kmeans() still refuses more centres than distinct points", {
  xy <- local({ .kf_seed(5); m <- cbind(stats::runif(5, 0, 100), stats::runif(5, 0, 100))
                m[rep(1:5, 4), ] })
  expect_error(local({ .kf_seed(1); .partition_centres(xy, 6L, init = "random", nstart = 10L,
                                                       one_call = TRUE) }),
               "more cluster centers than distinct data points")
  expect_null(local({ .kf_seed(1); .kmeans_best(xy, 6L, init = "random") }))
  expect_null(local({ .kf_seed(1); suppressWarnings(.kmeans_best(xy, 6L)) }))
  fit <- local({ .kf_seed(1); .partition_centres(xy, 5L, init = "random", nstart = 10L,
                                                 one_call = TRUE) })
  expect_identical(sort(fit$size), rep(4L, 5L))
  expect_identical(fit$wss, 0)
})


test_that("locations closer than a rounding step stay distinct, so every count that ran still runs", {
  # Forty locations, each twice, a micrometre apart: 3.5e-10 of the extent,
  # under the rounding step.  Rounded, some twins would be one point, and
  # stats::kmeans() would refuse counts the callers count as supported.
  m <- local({ .kf_seed(3); cbind(4.5e5 + stats::runif(40, 0, 2000),
                                  4.55e6 + stats::runif(40, 0, 2000)) })
  m <- rbind(m, m + 1e-6)
  b <- sf::st_as_sf(data.frame(x = m[, 1L], y = m[, 2L], z = local({ .kf_seed(4); stats::rnorm(80) })),
                    coords = c("x", "y"), crs = .kf_crs("m"))
  z <- .kmeans_frame(m)$z
  expect_identical(nrow(unique(z)), 80L)
  # Only the twins keep unrounded coordinates.
  fr <- .kmeans_frame(m[1:40, , drop = FALSE])
  expect_identical(fr$z, round(fr$z * 2^30) / 2^30)
  for (init in list(NULL, "kmeans++")) {
    s <- .kf_quiet(voronoi_seeds_kmeans(b, k = 70, partition = "kmeans", init = init))
    expect_identical(nrow(s), 70L)
    g <- .kf_quiet(get_voronoi_seeds(method = "kmeans", n = 70, sample_points = b,
                                     partition = "kmeans", init = init))
    expect_identical(nrow(g), 70L)
  }
  p <- .kf_quiet(resolution_profile(b, "z", estimate = NULL, obs_per_cell = NULL, sac = NULL,
                                    partition = "kmeans", min_cell_n = 1))
  expect_identical(max(p$levels), 79L)
  expect_false(anyNA(p$wss))
})


test_that("a lattice with repeat visits gets the same k-means cells in every frame, every row on the grid", {
  # Every location two or three times: the rule that keeps distinct locations
  # unrounded counts locations, not rows, so no row of a repeat visit is
  # left off the grid.
  g <- .kf_grid(40, 30, 250)
  m <- g[c(seq_len(nrow(g)), seq_len(nrow(g)), seq(1L, nrow(g), by = 2L)), ]
  xy <- .kf_xy(.kf_frames(.kf_layer(m)))
  for (r in names(xy)) {
    z <- .kmeans_frame(xy[[r]])$z
    expect_identical(z, round(z / 2^-30) * 2^-30, info = r)
  }
  for (k in c(5L, 25L, 60L)) {
    kp <- lapply(xy, function(a) { .kf_seed(1); .partition_centres(a, k, sort_rows = TRUE)$cluster })
    rn <- lapply(xy, function(a) {
      .kf_seed(1)
      suppressWarnings(.partition_centres(a, k, init = "random", nstart = 10L,
                                          iter.max = 100L, one_call = TRUE))$cluster
    })
    for (r in c("usft", "m01")) {
      expect_identical(kp[[r]], kp$m, info = paste(k, r, "kmeans++"))
      expect_identical(rn[[r]], rn$m, info = paste(k, r, "random"))
    }
  }
})

test_that("a layer at a single location gets one seed there under partition = \"kmeans\"", {
  one <- sf::st_as_sf(data.frame(x = rep(4.5e5, 10), y = rep(4.55e6, 10), z = 1:10),
                      coords = c("x", "y"), crs = .kf_crs("m"))
  z <- .kmeans_frame(sf::st_coordinates(one)[, 1:2])$z
  expect_identical(unname(z), matrix(0, 10L, 2L))
  expect_identical(.n_locations(sf::st_coordinates(one)), 1L)
  for (init in list(NULL, "kmeans++")) {
    info <- if (is.null(init)) "random starts" else init
    v <- .kf_quiet(voronoi_seeds_kmeans(one, k = 3, partition = "kmeans", init = init))
    g <- .kf_quiet(get_voronoi_seeds(method = "kmeans", n = 3, sample_points = one,
                                     partition = "kmeans", init = init))
    for (s in list(v, g)) {
      expect_identical(nrow(s), 1L, info = info)
      expect_equal(unname(sf::st_coordinates(s)[1L, 1:2]), c(4.5e5, 4.55e6), info = info)
    }
  }
})

test_that("the count of distinct locations is the same in metres, kilometres and feet", {
  # Forty locations, each twice, 2e-8 m apart.  The callers counted them
  # after rounding to 8 or 10 decimals of the layer's unit: 80 in metres
  # and fewer in kilometres, so the ladder's ceiling and the number of seeds
  # changed with the unit.
  m <- local({ .kf_seed(3); cbind(4.5e5 + stats::runif(40, 0, 2000),
                                  4.55e6 + stats::runif(40, 0, 2000)) })
  m <- rbind(m, m + 2e-8)
  b <- sf::st_as_sf(data.frame(x = m[, 1L], y = m[, 2L], z = local({ .kf_seed(4); stats::rnorm(80) })),
                    coords = c("x", "y"), crs = .kf_crs("m"))
  frames <- .kf_frames(b, c("m", "km", "usft"))
  for (r in names(frames)) {
    l <- frames[[r]]
    expect_identical(.n_locations(sf::st_coordinates(l)), 80L, info = r)
    expect_identical(nrow(.kf_quiet(voronoi_seeds_kmeans(l, k = 79, partition = "kmeans"))), 79L,
                     info = r)
    expect_identical(nrow(.kf_quiet(get_voronoi_seeds(method = "kmeans", n = 79, sample_points = l,
                                                      partition = "kmeans"))), 79L, info = r)
    for (kind in c("kmeans", "density")) {
      p <- .kf_quiet(resolution_profile(l, "z", estimate = NULL, obs_per_cell = NULL, sac = NULL,
                                        partition = kind, min_cell_n = 1))
      expect_identical(attr(p, "bounds")$n_distinct, 80L, info = paste(r, kind))
      expect_identical(attr(p, "bounds")$ceiling, 79L, info = paste(r, kind))
      expect_identical(max(p$levels), 79L, info = paste(r, kind))
    }
    d <- .kf_quiet(determine_optimal_levels(l, max_levels = 79, criterion = "morans_i",
                                            response_var = "z"))
    expect_length(attr(d, "diagnostics")$wss, 79L)
  }
})


# ---- The pieces --------------------------------------------------------------------

test_that("the frame of a lattice is the same in another unit and with the origin moved", {
  xy <- .kf_xy(.kf_frames(.kf_layer(.kf_grid(40, 30, 250)), c("m", "usft", "km", "m01")))
  z <- lapply(xy, function(a) .kmeans_frame(a)$z)
  for (r in names(xy)[-1L]) expect_identical(z[[r]], z$m, info = r)
  # Centred on the middle of the bounding box, scaled by its diagonal, on a
  # grid of 2^-30.
  expect_equal(apply(z$m, 2L, function(v) (min(v) + max(v)) / 2), c(0, 0))
  expect_equal(sqrt(sum(apply(z$m, 2L, function(v) diff(range(v)))^2)), 1, tolerance = 1e-9)
  expect_identical(z$m, round(z$m * 2^30) / 2^30)
  # The step the item fixes: some coordinate lies on 2^-30 and not on 2^-29.
  expect_identical(.KMEANS_STEP, 2^-30)
  expect_true(any(round(z$m * 2^29) != z$m * 2^29))
  # Points in the caller's coordinates map onto the same grid.
  fr <- .kmeans_frame(xy$m)
  expect_identical(fr$to(xy$m[1:7, , drop = FALSE]), fr$z[1:7, , drop = FALSE])
})

test_that(".kmeans_back() reports a fit in the frame in the caller's coordinates", {
  xy <- local({ .kf_seed(6); cbind(x = stats::runif(80, 0, 5e3), y = stats::runif(80, 0, 3e3)) })
  fr <- .kmeans_frame(xy)
  km <- local({ .kf_seed(2); stats::kmeans(fr$z, centers = 4L, nstart = 3L) })
  back <- .kmeans_back(km, xy)
  expect_identical(back$cluster, km$cluster)
  expect_identical(back$size, km$size)
  mu <- rowsum(xy, km$cluster) / as.vector(table(km$cluster))
  expect_equal(unname(back$centers), unname(mu), tolerance = 1e-12)
  expect_identical(dimnames(back$centers), dimnames(km$centers))
  w <- vapply(1:4, function(j) sum(sweep(xy[km$cluster == j, , drop = FALSE], 2L, mu[j, ])^2),
              numeric(1L))
  expect_equal(back$withinss, w, tolerance = 1e-12)
  expect_equal(back$tot.withinss, sum(w), tolerance = 1e-12)
  expect_equal(back$betweenss, back$totss - sum(w), tolerance = 1e-12)
})
