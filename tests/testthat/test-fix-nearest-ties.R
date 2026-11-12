# Nearest-point picks under the package's tie rule (BUG-110, BUG-111,
# BUG-172): a point the same distance from several others, or midpoints
# whose x differ only in their last bits, are settled by one rule that does
# not read those bits, so the result is the same in another unit or with
# another origin.

.nt_crs_m    <- function() sf::st_crs(26915)
.nt_crs_usft <- function() sf::st_crs("+proj=utm +zone=15 +datum=NAD83 +units=us-ft +no_defs")
.nt_crs_m01  <- function() sf::st_crs(paste0("+proj=tmerc +lat_0=0 +lon_0=-93 +k=0.9996 ",
                                            "+x_0=500000.1 +y_0=0.1 +datum=NAD83 +units=m +no_defs"))
.nt_frames <- function(pm) list(m = pm, usft = sf::st_transform(pm, .nt_crs_usft()),
                                m01 = sf::st_transform(pm, .nt_crs_m01()))

# A lattice of nx by ny points at 250 m, with a response z and a covariate w.
.nt_lattice <- function(nx, ny, jitter = 0) {
  g <- expand.grid(i = 0:(nx - 1), j = 0:(ny - 1))
  x <- 4.5e5 + g$i * 250; y <- 4.55e6 + g$j * 250
  if (jitter > 0) {
    x <- x + withr::with_seed(99, stats::runif(length(x), -jitter, jitter))
    y <- y + withr::with_seed(98, stats::runif(length(y), -jitter, jitter))
  }
  z <- withr::with_seed(1, sin((x - min(x)) / 1500) + cos((y - min(y)) / 1250) +
                          stats::rnorm(length(x), 0, 0.5))
  w <- withr::with_seed(2, (x - min(x)) / 1000 + stats::rnorm(length(x)))
  sf::st_as_sf(data.frame(x = x, y = y, z = z, w = w, i = g$i, j = g$j),
               coords = c("x", "y"), crs = .nt_crs_m())
}

.nt_quiet <- function(expr) suppressWarnings(suppressMessages(
  logger::with_log_threshold(expr, threshold = logger::FATAL, namespace = "spatialkit",
                             index = 2)))

.nt_surface_w <- function(cov, grid) {
  vapply(names(cov), function(f) {
    fit <- lm_spatial_fit(cov[[f]], "z", "w")
    .nt_quiet(predict_surface(fit, grid = grid[[f]], covariates = cov[[f]][, "w"]))$w
  }, numeric(nrow(grid$m)))
}

test_that("predict_surface() fills a grid point between lattice points from one of them in every frame", {
  cov <- .nt_frames(.nt_lattice(20, 15))
  g <- expand.grid(i = 0:18, j = 0:13)
  grid <- .nt_frames(sf::st_as_sf(
    data.frame(x = 4.5e5 + 125 + g$i * 250, y = 4.55e6 + 125 + g$j * 250),
    coords = c("x", "y"), crs = .nt_crs_m()))
  # Each grid point is 176.8 m from four covariate points.  The tie goes to
  # the lowest column, then the lowest y: the point below and to the left.
  want <- cov$m$w[g$i + 20L * g$j + 1L]
  w <- .nt_surface_w(cov, grid)
  # Before, 147 of 266 points took another row in US feet and 189 with the
  # origin moved.
  expect_identical(w[, "m"], want)
  expect_identical(w[, "usft"], want)
  expect_identical(w[, "m01"], want)
  # Without FNN the search is exact, and gives the same rows.
  testthat::local_mocked_bindings(.has_backend = function(pkg) FALSE)
  expect_identical(.nt_surface_w(cov, grid), w)
})

test_that("predict_surface()'s tie rule does not read the row order of `covariates`", {
  cov <- .nt_lattice(20, 15)
  g <- expand.grid(i = 0:18, j = 0:13)
  grid <- sf::st_as_sf(data.frame(x = 4.5e5 + 125 + g$i * 250, y = 4.55e6 + 125 + g$j * 250),
                       coords = c("x", "y"), crs = .nt_crs_m())
  fit <- lm_spatial_fit(cov, "z", "w")
  a <- .nt_quiet(predict_surface(fit, grid = grid, covariates = cov[, "w"]))
  o <- withr::with_seed(4, sample(nrow(cov)))
  b <- .nt_quiet(predict_surface(fit, grid = grid, covariates = cov[o, "w"]))
  expect_identical(b$w, a$w)
  expect_identical(b$.pred, a$.pred)
})

test_that("predict_surface() takes the nearest covariate point as before where no two tie", {
  # The jittered lattice of the scan: no grid point is the same distance from
  # two covariate points, so every row is sf::st_nearest_feature()'s.
  cov <- .nt_frames(.nt_lattice(20, 15, jitter = 60))
  g <- expand.grid(i = 0:18, j = 0:13)
  grid <- .nt_frames(sf::st_as_sf(
    data.frame(x = 4.5e5 + 125 + g$i * 250, y = 4.55e6 + 125 + g$j * 250),
    coords = c("x", "y"), crs = .nt_crs_m()))
  w <- .nt_surface_w(cov, grid)
  want <- cov$m$w[sf::st_nearest_feature(grid$m, cov$m)]
  for (f in names(cov)) expect_identical(w[, f], want, info = f)
  # Random points against random covariates.
  pts <- withr::with_seed(5, sf::st_as_sf(
    data.frame(x = stats::runif(300, 0, 5000), y = stats::runif(300, 0, 5000),
               z = stats::rnorm(300), w = stats::rnorm(300)),
    coords = c("x", "y"), crs = .nt_crs_m()))
  gr <- withr::with_seed(6, sf::st_as_sf(
    data.frame(x = stats::runif(500, 0, 5000), y = stats::runif(500, 0, 5000)),
    coords = c("x", "y"), crs = .nt_crs_m()))
  fit <- lm_spatial_fit(pts, "z", "w")
  s <- .nt_quiet(predict_surface(fit, grid = gr, covariates = pts[, "w"]))
  expect_identical(s$w, pts$w[sf::st_nearest_feature(gr, pts)])
  # An empty covariate row is never the nearest; an empty grid point gets NA.
  pe <- pts
  sf::st_geometry(pe)[[1]] <- sf::st_point()
  ge <- gr
  sf::st_geometry(ge)[[2]] <- sf::st_point()
  s2 <- .nt_quiet(predict_surface(fit, grid = ge, covariates = pe[, "w"]))
  expect_identical(s2$w, pe$w[sf::st_nearest_feature(ge, pe)])
  expect_true(is.na(s2$w[2]))
  expect_false(any(s2$w[-2] == pe$w[1]))
})

test_that("predict_surface() reads polygon covariates with sf::st_nearest_feature() as before", {
  cov <- .nt_lattice(6, 5)
  cells <- sf::st_buffer(cov, 125, endCapStyle = "SQUARE")
  gr <- sf::st_as_sf(data.frame(x = 4.5e5 + c(10, 260, 510), y = 4.55e6 + c(30, 280, 530)),
                     coords = c("x", "y"), crs = .nt_crs_m())
  fit <- lm_spatial_fit(cov, "z", "w")
  s <- .nt_quiet(predict_surface(fit, grid = gr, covariates = cells[, "w"]))
  expect_identical(s$w, cells$w[sf::st_nearest_feature(gr, cells)])
})

test_that("predict_surface() reads point covariates in degrees with sf::st_nearest_feature() as before", {
  # A 0.01 degree lattice with the grid at its midpoints: each grid point is
  # the same distance from four covariate points in degrees, but not on the
  # sphere, where the two to the north are nearer.  The help states that
  # points in degrees keep sf::st_nearest_feature()'s pick.
  g <- expand.grid(i = 0:19, j = 0:14)
  cov <- sf::st_as_sf(data.frame(x = -93 + g$i * 0.01, y = 41 + g$j * 0.01,
                                 z = withr::with_seed(1, stats::rnorm(300)),
                                 w = withr::with_seed(2, stats::rnorm(300))),
                      coords = c("x", "y"), crs = 4326)
  gg <- expand.grid(i = 0:18, j = 0:13)
  grid <- sf::st_as_sf(data.frame(x = -93 + 0.005 + gg$i * 0.01, y = 41 + 0.005 + gg$j * 0.01),
                       coords = c("x", "y"), crs = 4326)
  fit <- lm_spatial_fit(cov, "z", "w")
  s <- .nt_quiet(predict_surface(fit, grid = grid, covariates = cov[, "w"]))
  nf <- sf::st_nearest_feature(grid, cov)
  expect_identical(s$w, cov$w[nf])
  # The planar rule would take the point below and to the left instead.
  expect_true(all(g$j[nf] == gg$j + 1L))
  expect_false(any(s$w == cov$w[gg$i + 20L * gg$j + 1L]))
})

test_that("the kmeans record of scored seeds gives a tied point to the lower seed_id in every frame", {
  # A cloud on a lattice and seeds on every second point of it: a cloud point
  # between two or four seeds is the same distance from each.
  g  <- expand.grid(i = 0:19, j = 0:14)
  cl <- sf::st_as_sf(data.frame(x = 4.5e5 + g$i * 250, y = 4.55e6 + g$j * 250),
                     coords = c("x", "y"), crs = .nt_crs_m())
  s  <- expand.grid(i = seq(0, 19, by = 2), j = seq(0, 14, by = 2))
  sd <- sf::st_as_sf(data.frame(x = 4.5e5 + s$i * 250, y = 4.55e6 + s$j * 250),
                     coords = c("x", "y"), crs = .nt_crs_m())
  C <- .nt_frames(cl); S <- .nt_frames(sd)
  rec <- lapply(names(C), function(f) .scored_seeds_record(C[[f]], S[[f]], 5L))
  names(rec) <- names(C)
  # The nearest seed, a tie to the lower seed_id, worked out on the lattice
  # indices: the seed at the lower even i, then the lower even j.
  want <- match(paste(2L * (g$i %/% 2L), 2L * (g$j %/% 2L)), paste(s$i, s$j))
  # Before, 81, 107 and 112 rows went to another seed.
  for (f in names(C)) {
    expect_identical(rec[[f]]$cluster, want, info = f)
    expect_identical(rec[[f]]$size, rec$m$size, info = f)
    expect_identical(rec[[f]]$rows, seq_len(nrow(cl)), info = f)
  }
  usft <- 3937 / 1200
  expect_equal(rec$usft$tot_withinss / rec$m$tot_withinss, usft^2, tolerance = 1e-9)
  expect_equal(rec$m01$withinss, rec$m$withinss, tolerance = 1e-9)
  # The same without FNN.
  testthat::local_mocked_bindings(.has_backend = function(pkg) FALSE)
  expect_identical(.scored_seeds_record(C$usft, S$usft, 5L)$cluster, want)
})

test_that("the kmeans record of scored seeds never assigns a point to a seed without a location", {
  cl <- sf::st_as_sf(data.frame(x = c(0, 10, 20, 30, 40), y = 0), coords = c("x", "y"),
                     crs = .nt_crs_m())
  sd <- sf::st_sfc(sf::st_point(c(0, 0)), sf::st_point(), sf::st_point(c(40, 0)),
                   crs = .nt_crs_m())
  rec <- .scored_seeds_record(cl, sd, 5L)
  # The point at 20 is the same distance from seeds 1 and 3 and goes to 1.
  expect_identical(rec$cluster, c(1L, 1L, 1L, 3L, 3L))
  expect_identical(rec$size, c(3L, 0L, 2L))
  expect_identical(rec$withinss, c(500, 0, 100))
})

test_that("get_voronoi_seeds()'s record for a profile's centres holds the profile's cells in every frame", {
  F <- .nt_frames(.nt_lattice(40, 30))
  out <- lapply(F, function(l) {
    prof <- .nt_quiet(resolution_profile(l, "z", seed = 1))
    s <- .nt_quiet(get_voronoi_seeds(method = "kmeans", n = prof, sample_points = l))
    cen <- .scored_centres(attr(s, "count_rule"))
    xy <- sf::st_coordinates(l)
    list(rec = attr(s, "kmeans"),
         rule = .nearest_centre(xy[, 1], xy[, 2], cen[, 1], cen[, 2], use_fnn = FALSE))
  })
  # Before, the record left the profile's cells in 4, 1 and 4 rows, and its
  # cluster differed from the one in metres in 5 rows (US feet) and 2
  # (origin moved).
  for (f in names(F)) {
    expect_identical(out[[f]]$rec$cluster, out[[f]]$rule, info = f)
    expect_identical(out[[f]]$rec$cluster, out$m$rec$cluster, info = f)
    expect_identical(out[[f]]$rec$size, out$m$rec$size, info = f)
  }
  expect_equal(out$usft$rec$tot_withinss / out$m$rec$tot_withinss, (3937 / 1200)^2,
               tolerance = 1e-9)
  expect_equal(out$m01$rec$tot_withinss, out$m$rec$tot_withinss, tolerance = 1e-9)
})

test_that("coerce_to_points() keeps the lower of two equal parts whose midpoints share x at every origin", {
  # Part a runs up and to the right; part b, 1000 higher, runs to the left
  # with another slope.  Both are 100.04 long and have midpoint x 30.15, but
  # reached from other coordinates the two x values differ in their last bits
  # once the origin moves, and those bits chose the part: 33 of 200 origins
  # took the upper part.
  a <- rbind(c(0.13, 0.29), c(60.17, 80.31))
  L <- sqrt(sum(diff(a)^2)); dy <- sqrt(L^2 - 30^2)
  b <- rbind(c(45.15, 1000), c(15.15, 1000 + dy))
  withr::local_seed(3, .rng_kind = "Mersenne-Twister",
                    .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
  res <- vapply(1:100, function(k) {
    o <- c(stats::runif(1, 1e5, 9e5) + stats::runif(1),
           stats::runif(1, 1e6, 9e6) + stats::runif(1))
    pa <- sweep(a, 2, o, "+"); pb <- sweep(b, 2, o, "+")
    vapply(list(list(pa, pb), list(pb, pa)), function(parts) {
      x <- sf::st_sf(id = 1, geometry = sf::st_sfc(sf::st_multilinestring(parts), crs = 32617))
      sf::st_coordinates(coerce_to_points(x, "auto"))[1, 2] - o[2]
    }, numeric(1))
  }, numeric(2))
  # The midpoint of the lower part, at y 40.3, in both part orders.
  expect_true(all(abs(res - 40.3) < 1e-6))
})

test_that("coerce_to_points() reads tied midpoints' x by the box's diagonal, not by a fixed rounding", {
  # Two parts of one length whose midpoint x values, 30.1500005001 and
  # 30.1500004999, differ by 2e-10 across the sixth decimal: one column to one
  # part in ten million of the midpoints' box (about 1000), so the lower y,
  # part a, is kept in both part orders.  Ordered by the exact x, or by x
  # rounded to six decimals, part b was taken.
  a <- rbind(c(0, 0), c(60.3000010002, 80))
  L <- sqrt(sum(diff(a)^2)); dy <- sqrt(L^2 - 60^2)
  b <- rbind(c(30.1500004999 - 30, 1000), c(30.1500004999 + 30, 1000 + dy))
  for (parts in list(list(a, b), list(b, a))) {
    x <- sf::st_sf(id = 1, geometry = sf::st_sfc(sf::st_multilinestring(parts), crs = 32617))
    xy <- sf::st_coordinates(coerce_to_points(x, "auto"))[1, 1:2]
    expect_equal(unname(xy), c(30.1500005001, 40), tolerance = 1e-12)
  }
})

test_that("coerce_to_points() keeps the southern of two equal parts in degrees whose midpoints share x", {
  # Parts symmetric about the central meridian of UTM zone 17 and about the
  # equator: one length in the temporary projection and one midpoint
  # longitude, -81, reached from other coordinates.  Before, 8 of 60 such
  # features took the northern part.
  withr::local_seed(11, .rng_kind = "Mersenne-Twister",
                    .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
  todeg <- function(m)
    sf::st_coordinates(sf::st_transform(sf::st_sfc(sf::st_multipoint(m), crs = 32617), 4326))[, 1:2]
  res <- vapply(1:60, function(k) {
    p <- stats::runif(1, 100, 5000); n1 <- stats::runif(1, 1e5, 3e6)
    n2 <- n1 + stats::runif(1, -3000, 3000)
    g1 <- todeg(rbind(c(5e5 - p, n1), c(5e5 + p, n2)))
    g2 <- todeg(rbind(c(5e5 + p, -n1), c(5e5 - p, -n2)))
    vapply(list(list(g1, g2), list(g2, g1)), function(parts) {
      x <- sf::st_sf(id = 1, geometry = sf::st_sfc(sf::st_multilinestring(parts), crs = 4326))
      xy <- sf::st_coordinates(suppressWarnings(coerce_to_points(x, "auto")))[1, 1:2]
      c(xy[[2]] < 0, abs(xy[[1]] + 81) < 1e-6)
    }, logical(2))
  }, logical(4))
  expect_true(all(res))
})

test_that("the help pages state the tie rules", {
  page <- function(name) {
    file <- test_path("..", "..", "man", name)
    rd <- if (file.exists(file)) tools::parse_Rd(file)
          else tools::Rd_db("spatialkit")[[name]]
    expect_s3_class(rd, "Rd")
    gsub("[[:space:]]+", " ", paste(unlist(rd), collapse = ""))
  }
  ps <- page("predict_surface.Rd")
  expect_match(ps, paste0("Distances count as equal when the squared distance is no more than the ",
                          "square of the least distance lengthened by a margin, plus one part in ",
                          "ten million of the least squared distance. Among them the point lowest ",
                          "in x is taken, then the one lowest in y, then the first in row order. ",
                          "x values count as one when they lie within one part in ten million of ",
                          "the diagonal of the points' bounding box, or within the margin where ",
                          "that is more; y values are read the same way."), fixed = TRUE)
  expect_match(ps, paste0("It is one part in 1e10 of the largest coordinate of the features that ",
                          "have a location, at most a millionth of the diagonal of their bounding ",
                          "box and at least four times the spacing of doubles at that ",
                          "coordinate."), fixed = TRUE)
  # sf::st_nearest_feature() searches a layer in degrees with s2, not GEOS,
  # while s2 is on.
  expect_match(ps, paste0("order its search meets the features (GEOS, or s2 while ",
                          "sf::sf_use_s2() is on)"), fixed = TRUE)
  expect_match(page("coerce_to_points.Rd"),
               paste0("x values within one part in ten million of the diagonal of the ",
                      "midpoints' bounding box count as one, so y decides"), fixed = TRUE)
  expect_match(page("get_voronoi_seeds.Rd"),
               paste0("A point the same distance from two seeds (squared distances within ",
                      "one part in ten million) goes to the lower seed_id"), fixed = TRUE)
})
