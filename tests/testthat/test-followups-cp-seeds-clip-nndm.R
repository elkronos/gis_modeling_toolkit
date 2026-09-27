# tests/testthat/test-followups-cp-seeds-clip-nndm.R
# ---------------------------------------------------------------------------
# resolution_profile() gives Cp a noise variance when the variogram has no
# nugget to give; the seeding functions hand on the partition the profile
# scored rather than a k-means of their own; a clipped cell whose overlap
# with the boundary is an area plus a line keeps its area; and the NNDM
# min_train warning states the statistic that fired it.
# ---------------------------------------------------------------------------

.fx_pts <- function(n = 300, seed = 1) {
  set.seed(seed)
  d <- data.frame(x = 5e5 + runif(n, 0, 3000), y = 5e6 + runif(n, 0, 2000))
  d$z <- sin((d$x - 5e5) / 400) + rnorm(n, 0, 0.3)
  sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
}

.fx_warnings <- function(expr) {
  ws <- character(0)
  value <- withCallingHandlers(
    suppressMessages(expr),
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(value = value, warnings = ws)
}


# ---- Cp's noise variance without a nugget -----------------------------------

test_that("with no nugget, cp takes Mallows' residual mean square of the finest level", {
  d <- .fx_pts()
  # A range alone: no variogram model, so no nugget.
  out <- .fx_warnings(resolution_profile(d, "z", sac = 400, n_levels = 6, nstart = 5))
  p <- out$value
  fb <- grep("`cp` takes its noise variance from the residual mean square", out$warnings,
             value = TRUE)
  expect_length(fb, 1L)
  expect_match(fb, "`sac` carries no variogram model and so no nugget, so `cp`", fixed = TRUE)
  expect_true(all(is.finite(p$cp)))
  expect_true(all(is.na(p$reliability)))

  cn <- attr(p, "cp_noise")
  expect_identical(cn$source, "finest-level residual mean square")
  # Every row is scored and every cell holds one (300 points, no subsample),
  # so m = n and L_m = L; the finest level whose cells hold two rows on
  # average is the last one with n - L >= L.
  n <- nrow(d); L <- p$levels
  j <- max(which(n - L >= L))
  expect_identical(cn$level, as.integer(L[j]))
  expect_equal(cn$value, p$rss[j] / (n - L[j]))
  expect_equal(p$cp, p$rss / n + cn$value * 2 * L / n)

  expect_output(print(p), "none usable \\(reliability is NA; cp's noise variance")
  expect_s3_class(select_resolution(p, "cp"), "resolution_selection")
})

test_that("a nugget the variogram does give is still the one cp uses", {
  d <- .fx_pts()
  vm <- data.frame(model = c("Nug", "Exp"), psill = c(0.3, 0.7), range = c(0, 200),
                   stringsAsFactors = FALSE)
  sac <- structure(600, class = c("sac_range", "numeric"), variogram_model = vm,
                   crs = sf::st_crs(32632))
  p <- suppressWarnings(suppressMessages(
    resolution_profile(d, "z", sac = sac, n_levels = 6, nstart = 5)))
  cn <- attr(p, "cp_noise")
  expect_identical(cn$source, "variogram nugget")
  expect_equal(cn$value, 0.3)
  expect_equal(p$cp, p$rss / nrow(d) + 0.3 * 2 * p$levels / nrow(d))
})

test_that("a geometry-only profile has no cp and no noise record", {
  p <- resolution_profile(.fx_pts(), n_levels = 4, nstart = 2)
  expect_true(all(is.na(p$cp)))
  expect_null(attr(p, "cp_noise"))
})


# ---- The scored partition as seeds ------------------------------------------

test_that("select_resolution() carries the centres of the partition it scored", {
  d <- .fx_pts()
  p <- suppressWarnings(suppressMessages(
    resolution_profile(d, "z", sac = 400, n_levels = 6, nstart = 5)))
  sel <- select_resolution(p, "cp")
  s <- sel$seeds
  expect_s3_class(s, "sf")
  expect_identical(nrow(s), as.integer(sel$best))
  expect_identical(s$seed_id, seq_len(sel$best))
  expect_equal(sf::st_crs(s), sf::st_crs(d))
  expect_equal(unname(sf::st_coordinates(s)[, 1:2]),
               attr(p, "centres")[[as.character(sel$best)]])
})

test_that("get_voronoi_seeds() and voronoi_seeds_kmeans() return the scored centres", {
  d <- .fx_pts()
  p <- suppressWarnings(suppressMessages(
    resolution_profile(d, "z", sac = 400, n_levels = 6, nstart = 5)))
  sel <- select_resolution(p, "cp")
  cen <- attr(p, "centres")[[as.character(sel$best)]]

  s1 <- get_voronoi_seeds(method = "kmeans", n = sel, sample_points = d)
  expect_equal(unname(sf::st_coordinates(s1)[, 1:2]), cen)
  expect_identical(s1$method, rep("kmeans", nrow(s1)))
  expect_match(attr(s1, "n_from"), "select_resolution\\(criterion = \"cp\"\\)")
  # Each sample point falls to its nearest seed, and each seed is the mean of
  # the points that fall to it: the seeds' Voronoi cells are the partition
  # the profile scored, not a partition near it.
  km <- attr(s1, "kmeans")
  expect_identical(km$rows, seq_len(nrow(d)))
  expect_true(is.na(km$iter))
  expect_identical(km$nstart, 5L)
  xy <- sf::st_coordinates(d)
  means <- cbind(tapply(xy[, 1], km$cluster, mean), tapply(xy[, 2], km$cluster, mean))
  expect_equal(unname(means), cen, tolerance = 1e-8)

  # The same from the profile itself, which is read at its default criterion.
  s2 <- get_voronoi_seeds(method = "kmeans", n = p, sample_points = d)
  expect_equal(unname(sf::st_coordinates(s2)[, 1:2]),
               unname(sf::st_coordinates(select_resolution(p)$seeds)[, 1:2]))

  # voronoi_seeds_kmeans() takes the selection too, and returns the seeds in
  # the points' CRS.
  s3 <- voronoi_seeds_kmeans(d, k = sel)
  expect_equal(unname(sf::st_coordinates(s3)[, 1:2]), cen)
  s4 <- voronoi_seeds_kmeans(sf::st_transform(d, 4326), k = sel)
  expect_identical(sf::st_crs(s4)$epsg, 4326L)
  expect_equal(unname(sf::st_coordinates(sf::st_transform(s4, 32632))[, 1:2]), cen,
               tolerance = 1e-6)
})

test_that("a count as a number, or a selection without seeds, still runs a k-means", {
  d <- .fx_pts()
  p <- suppressWarnings(suppressMessages(
    resolution_profile(d, "z", sac = 400, n_levels = 6, nstart = 5)))
  sel <- select_resolution(p, "cp")
  fresh <- get_voronoi_seeds(method = "kmeans", n = sel$best, sample_points = d,
                             set_seed = 1)
  expect_false(is.na(attr(fresh, "kmeans")$iter))
  old <- sel
  old$seeds <- NULL                 # a selection made before seeds were kept
  s <- get_voronoi_seeds(method = "kmeans", n = old, sample_points = d, set_seed = 1)
  expect_false(is.na(attr(s, "kmeans")$iter))
  expect_identical(nrow(s), as.integer(sel$best))
  expect_identical(nrow(voronoi_seeds_kmeans(d, k = old)), as.integer(sel$best))
})


# ---- A clipped cell that is an area plus a line -----------------------------

.fx_sq <- function(x0, x1, y0, y1)
  sf::st_polygon(list(rbind(c(x0, y0), c(x1, y0), c(x1, y1), c(x0, y1), c(x0, y0))))

test_that("a grid cell over the inner corner of an L keeps its area and its points", {
  # The L: a bar up the left side and one across the top.  The cell [0,2]^2
  # meets it in the left bar's foot plus the line y = 2, 1 <= x <= 2, where
  # it touches the top bar from below: a GEOMETRYCOLLECTION, which used to be
  # dropped whole with the ten points in it.
  L <- sf::st_union(sf::st_sfc(.fx_sq(0, 1, 0, 2), .fx_sq(0, 4, 2, 4)))
  bnd <- sf::st_sf(geometry = sf::st_sfc(L[[1]], crs = 32632))
  set.seed(1)
  pts <- sf::st_as_sf(data.frame(x = c(runif(10, 0.05, 0.95), runif(20, 0.1, 3.9)),
                                 y = c(runif(10, 0.1, 1.9), runif(20, 2.1, 3.9))),
                      coords = c("x", "y"), crs = 32632)
  tt <- build_tessellation(pts, boundary = bnd, method = "square", cellsize = 2,
                           quiet = TRUE)
  expect_identical(nrow(tt$cells), 3L)
  expect_false(anyNA(tt$index))
  expect_true(all(as.character(sf::st_geometry_type(tt$cells)) %in%
                    c("POLYGON", "MULTIPOLYGON")))
  expect_equal(sum(as.numeric(sf::st_area(tt$cells))), as.numeric(sf::st_area(bnd)))
})

test_that(".keep_polygonal() keeps a collection's area and drops what has none", {
  gc_area <- sf::st_geometrycollection(list(
    .fx_sq(0, 1, 0, 1), sf::st_linestring(rbind(c(1, 1), c(2, 1)))))
  gc_two  <- sf::st_geometrycollection(list(
    .fx_sq(0, 1, 0, 1), .fx_sq(2, 3, 0, 1), sf::st_point(c(5, 5))))
  gc_none <- sf::st_geometrycollection(list(
    sf::st_linestring(rbind(c(0, 0), c(1, 0))), sf::st_point(c(3, 3))))
  x <- sf::st_sf(id = 1:5,
                 geometry = sf::st_sfc(gc_area, .fx_sq(0, 1, 0, 1), gc_none,
                                       sf::st_linestring(rbind(c(0, 0), c(1, 1))),
                                       gc_two, crs = 32632))
  out <- spatialkit:::.keep_polygonal(x)
  expect_identical(out$id, c(1L, 2L, 5L))
  expect_identical(as.character(sf::st_geometry_type(out)),
                   c("POLYGON", "POLYGON", "MULTIPOLYGON"))
  expect_equal(as.numeric(sf::st_area(out)), c(1, 1, 2))
  expect_equal(sf::st_crs(out), sf::st_crs(x))
})


# ---- The NNDM min_train warning ---------------------------------------------

test_that("the nndm min_train warning names the ECDF excess that fired it", {
  set.seed(1)
  mk <- function(x, y) sf::st_as_sf(data.frame(x = x, y = y), coords = c("x", "y"),
                                    crs = 3857)
  p <- mk(rnorm(100, 10000, 800), rnorm(100, 10000, 800))
  g <- mk(rep(seq(0, 20000, by = 1000), 21), rep(seq(0, 20000, by = 1000), each = 21))
  out <- .fx_warnings(make_folds(p, method = "nndm", prediction_points = g))
  w <- grep("stopped the distance matching", out$warnings, value = TRUE)
  expect_length(w, 1L)
  expect_match(w, "At distances up to phi \\([0-9.e+]+\\) the share of folds")
  expect_match(w, "runs above the target's share by as much as 0\\.[0-9]{3} \\(more than one fold's worth, 1/100\\)")
  expect_match(w, "need not differ the same way", fixed = TRUE)
})
