# ===========================================================================
# R/seeding.R: get_voronoi_seeds(), voronoi_seeds_kmeans(),
# voronoi_seeds_random() and .robust_st_sample().
#
# Three exported seed generators with zero coverage.  Two properties matter
# most and neither is visible from a single call:
#
#   * reproducibility -- the same seed must give the same seeds, or a
#     tessellation cannot be reproduced from a script;
#   * RNG restoration -- seeding must not leak into the caller's stream.
#     test-fold-methods.R documents a .with_seed() on.exit collision that
#     silently replaced the caller's RNG state once already in this codebase.
# ===========================================================================

.seed_sq <- function(x0 = 0, y0 = 0, x1 = 100, y1 = 100, crs = 32632) {
  sf::st_sf(geometry = sf::st_sfc(
    sf::st_polygon(list(rbind(c(x0, y0), c(x1, y0), c(x1, y1), c(x0, y1),
                              c(x0, y0)))),
    crs = crs))
}

.seed_points <- function(n = 50, seed = 4, crs = 32632) {
  set.seed(seed)
  sf::st_as_sf(
    data.frame(x = runif(n, 0, 100), y = runif(n, 0, 100)),
    coords = c("x", "y"), crs = crs
  )
}

# before/after around a seeded call, with the caller's stream pinned.
.rng_untouched <- function(expr) {
  set.seed(20250101); expected <- stats::runif(3)
  set.seed(20250101); force(expr)
  isTRUE(all.equal(stats::runif(3), expected))
}


# ---------------------------------------------------------------------------
# Shared output contract
# ---------------------------------------------------------------------------

test_that("all three seed generators emit the same seed_id / method columns", {
  bnd <- .seed_sq()
  pts <- .seed_points()

  cases <- list(
    random    = get_voronoi_seeds(bnd, method = "random", n = 5, set_seed = 1),
    kmeans    = get_voronoi_seeds(bnd, method = "kmeans", n = 4, set_seed = 2),
    provided  = get_voronoi_seeds(seeds = pts[1:3, ], method = "provided"),
    vs_kmeans = voronoi_seeds_kmeans(pts, k = 6, set_seed = 3),
    vs_random = voronoi_seeds_random(bnd, k = 7, set_seed = 3)
  )
  expected_method <- c(random = "random", kmeans = "kmeans",
                       provided = "provided", vs_kmeans = "kmeans",
                       vs_random = "random")

  for (nm in names(cases)) {
    s <- cases[[nm]]
    expect_s3_class(s, "sf")
    expect_true(all(c("seed_id", "method") %in% names(s)), info = nm)
    # seed_id is a dense 1..n key, so the three functions are interchangeable.
    expect_equal(s$seed_id, seq_len(nrow(s)), info = nm)
    expect_equal(unique(s$method), unname(expected_method[nm]), info = nm)
    expect_true(all(as.character(sf::st_geometry_type(s, by_geometry = TRUE)) ==
                      "POINT"), info = nm)
    expect_equal(sf::st_crs(s), sf::st_crs(32632), info = nm)
  }

  expect_equal(nrow(cases$random), 5L)
  expect_equal(nrow(cases$kmeans), 4L)
  expect_equal(nrow(cases$provided), 3L)
  expect_equal(nrow(cases$vs_kmeans), 6L)
  expect_equal(nrow(cases$vs_random), 7L)
})


# ---------------------------------------------------------------------------
# get_voronoi_seeds()
# ---------------------------------------------------------------------------

test_that("get_voronoi_seeds(random) is reproducible and stays inside the boundary", {
  bnd <- .seed_sq()

  a <- get_voronoi_seeds(bnd, method = "random", n = 8, set_seed = 11)
  b <- get_voronoi_seeds(bnd, method = "random", n = 8, set_seed = 11)
  expect_equal(sf::st_coordinates(a), sf::st_coordinates(b))

  c_ <- get_voronoi_seeds(bnd, method = "random", n = 8, set_seed = 12)
  expect_false(isTRUE(all.equal(sf::st_coordinates(a), sf::st_coordinates(c_))))

  # Sampled within the polygon, not merely within its bbox by luck.
  expect_true(all(lengths(sf::st_within(a, bnd)) > 0L))
})


test_that("get_voronoi_seeds(kmeans) clusters the sampling cloud", {
  bnd <- .seed_sq()

  a <- get_voronoi_seeds(bnd, method = "kmeans", n = 6, set_seed = 21)
  b <- get_voronoi_seeds(bnd, method = "kmeans", n = 6, set_seed = 21)
  expect_equal(sf::st_coordinates(a), sf::st_coordinates(b))
  expect_equal(nrow(a), 6L)

  # Centres of a uniform cloud over a square land inside the square.
  expect_true(all(lengths(sf::st_within(a, bnd)) > 0L))

  # An explicit sampling cloud is used in place of the internal one, so two
  # clearly separated clusters give two clearly separated centres.
  cloud <- sf::st_as_sf(
    data.frame(x = c(rep(5, 20), rep(95, 20)) + rep(c(-1, 1), 20),
               y = c(rep(5, 20), rep(95, 20))),
    coords = c("x", "y"), crs = 32632)
  centres <- get_voronoi_seeds(bnd, method = "kmeans", n = 2,
                               sample_points = cloud, set_seed = 5)
  xy <- sf::st_coordinates(centres)
  expect_equal(nrow(centres), 2L)
  expect_equal(sort(round(xy[, 1])), c(5, 95))
})


test_that("get_voronoi_seeds(provided) returns every supplied seed and warns about n", {
  seeds <- .seed_points(4, seed = 6)

  out <- get_voronoi_seeds(seeds = seeds, method = "provided")
  expect_equal(nrow(out), 4L)
  expect_equal(sf::st_coordinates(out)[, 1:2],
               sf::st_coordinates(seeds)[, 1:2])

  # `n` is documented as IGNORED here; a mismatch is announced rather than
  # silently truncating or padding.
  lines <- capture_spatialkit_log(
    out_n <- get_voronoi_seeds(seeds = seeds, method = "provided", n = 9)
  )
  expect_true(log_has(lines, "method = 'provided' ignores 'n'"))
  expect_equal(nrow(out_n), 4L)

  # A matching n says nothing.
  quiet <- capture_spatialkit_log(
    get_voronoi_seeds(seeds = seeds, method = "provided", n = 4))
  expect_false(log_has(quiet, "ignores 'n'"))
})


test_that("get_voronoi_seeds validates its inputs", {
  bnd <- .seed_sq()
  pts <- .seed_points(10)

  expect_error(get_voronoi_seeds(bnd, method = "random"), "'n' is required")
  expect_error(get_voronoi_seeds(bnd, method = "kmeans"), "'n' is required")
  expect_error(get_voronoi_seeds(method = "provided"), "provide 'seeds'")
  expect_error(get_voronoi_seeds(method = "random", n = 3),
               "random seeds require 'boundary'")
  expect_error(get_voronoi_seeds(method = "kmeans", n = 3),
               "k-means seeds require 'boundary'")
  # Geometry types are checked through .assert_sf().
  expect_error(get_voronoi_seeds(pts, method = "random", n = 3),
               "`boundary` geometry must be one of: POLYGON")
  expect_error(get_voronoi_seeds(seeds = bnd, method = "provided"),
               "`seeds` geometry must be one of: POINT")
  expect_error(get_voronoi_seeds(bnd, method = "kmeans", n = 3,
                                 sample_points = bnd),
               "`sample_points` geometry must be one of: POINT")
  expect_error(get_voronoi_seeds(bnd, method = "nonsense", n = 3),
               "'arg' should be one of")
})


# ---------------------------------------------------------------------------
# voronoi_seeds_kmeans() / voronoi_seeds_random()
# ---------------------------------------------------------------------------

test_that("voronoi_seeds_kmeans reproduces its centres and validates input", {
  pts <- .seed_points(60, seed = 7)

  a <- voronoi_seeds_kmeans(pts, k = 5, set_seed = 13)
  b <- voronoi_seeds_kmeans(pts, k = 5, set_seed = 13)
  expect_equal(sf::st_coordinates(a), sf::st_coordinates(b))
  expect_false(isTRUE(all.equal(sf::st_coordinates(a),
                                sf::st_coordinates(
                                  voronoi_seeds_kmeans(pts, k = 5,
                                                       set_seed = 14)))))

  # The siblings now validate spatial input the same way get_voronoi_seeds()
  # does, instead of failing somewhere inside st_coordinates().
  expect_error(voronoi_seeds_kmeans(.seed_sq(), k = 3),
               "`points_sf` geometry must be one of: POINT")
  expect_error(voronoi_seeds_kmeans(sf::st_drop_geometry(pts), k = 3),
               "Expected an sf object")
  # A zero-row layer never reaches the clustering code: .assert_sf() rejects it
  # first, because an empty geometry set has no type to check.  Without the
  # regexp this expectation passed on that message and looked like coverage of
  # the "nothing to cluster" branch, which it is not -- see the next test.
  expect_error(voronoi_seeds_kmeans(pts[0, , drop = FALSE], k = 3),
               "`points_sf` geometry must be one of: POINT")
})

test_that("voronoi_seeds_kmeans drops unusable points and reports having none", {
  # st_coordinates() yields one ALL-NA row per EMPTY POINT rather than zero
  # rows, so a layer of empty points has nrow() > 0, passes .assert_sf() (its
  # geometry type is POINT), and used to reach stats::kmeans() as a matrix of
  # NAs -- which fails with "NA/NaN/Inf in foreign function call (arg 1)".
  ll <- sf::st_crs(32632)

  # (a) An empty point INSIDE a populated layer is dropped, and the rest
  #     clusters exactly as if it had never been there.
  set.seed(3)
  n  <- 20
  xy <- data.frame(x = runif(n, 0, 100), y = runif(n, 0, 100))
  good <- sf::st_as_sf(xy, coords = c("x", "y"), crs = ll)
  mixed <- rbind(
    good,
    sf::st_sf(geometry = sf::st_sfc(sf::st_point(), crs = ll))
  )
  expect_equal(nrow(mixed), n + 1L)
  expect_true(any(sf::st_is_empty(mixed)))

  lines <- capture_spatialkit_log(out <- voronoi_seeds_kmeans(mixed, k = 3,
                                                              set_seed = 11))
  expect_true(log_has(lines, "dropping 1 point\\(s\\) with empty or non-finite"))
  expect_equal(nrow(out), 3L)
  expect_true(all(is.finite(sf::st_coordinates(out))))
  # Identical to clustering the clean layer: the empty row contributed nothing.
  expect_equal(sf::st_coordinates(out),
               sf::st_coordinates(voronoi_seeds_kmeans(good, k = 3,
                                                       set_seed = 11)))

  # (b) A layer of NOTHING BUT empty points is the only way to reach the
  #     "nothing to cluster" branch -- a zero-row layer is stopped earlier.
  allempty <- sf::st_sf(
    id = 1:3,
    geometry = sf::st_sfc(sf::st_point(), sf::st_point(), sf::st_point(),
                          crs = ll))
  expect_equal(nrow(allempty), 3L)
  expect_error(voronoi_seeds_kmeans(allempty, k = 2),
               "has no usable coordinates; nothing to cluster")
})


test_that("voronoi_seeds_kmeans clamps k to the distinct positions and says so", {
  # 12 rows but only 3 distinct positions: k-means cannot produce 10 centres,
  # and stats::kmeans() would error rather than explain.
  dup <- sf::st_as_sf(
    data.frame(x = rep(c(1, 2, 3), 4), y = rep(c(1, 2, 3), 4)),
    coords = c("x", "y"), crs = 32632)

  lines <- capture_spatialkit_log(out <- voronoi_seeds_kmeans(dup, k = 10))
  expect_true(log_has(lines, "only 3 unique positions among 12 point\\(s\\)"))
  expect_true(log_has(lines, "clamping"))
  expect_equal(nrow(out), 3L)
  expect_equal(out$seed_id, 1:3)

  # No clamp, no message.
  quiet <- capture_spatialkit_log(ok <- voronoi_seeds_kmeans(dup, k = 3))
  expect_false(log_has(quiet, "clamping"))
  expect_equal(nrow(ok), 3L)
})


test_that("voronoi_seeds_kmeans projects lon/lat input and returns it in the input CRS", {
  # k-means in degrees is distance-distorted; the centres must nevertheless
  # come back in the CRS they were supplied in.
  set.seed(8)
  ll <- sf::st_as_sf(
    data.frame(lon = runif(60, 9.0, 9.4), lat = runif(60, 48.6, 48.9)),
    coords = c("lon", "lat"), crs = 4326)

  out <- voronoi_seeds_kmeans(ll, k = 4, set_seed = 15)
  expect_equal(sf::st_crs(out), sf::st_crs(4326))
  expect_equal(nrow(out), 4L)
  xy <- sf::st_coordinates(out)
  expect_true(all(xy[, 1] > 8.9 & xy[, 1] < 9.5))
  expect_true(all(xy[, 2] > 48.5 & xy[, 2] < 49.0))
})


test_that("voronoi_seeds_random accepts sf or sfc and stays inside the boundary", {
  bnd <- .seed_sq()

  a <- voronoi_seeds_random(bnd, k = 9, set_seed = 17)
  b <- voronoi_seeds_random(sf::st_geometry(bnd), k = 9, set_seed = 17)

  expect_equal(nrow(a), 9L)
  expect_equal(sf::st_coordinates(a), sf::st_coordinates(b))
  expect_equal(sf::st_crs(a), sf::st_crs(bnd))
  expect_true(all(lengths(sf::st_within(a, bnd)) > 0L))

  expect_error(voronoi_seeds_random(.seed_points(5), k = 3),
               "`boundary` geometry must be one of: POLYGON")
  expect_error(voronoi_seeds_random(list(), k = 3), "Expected an sf object")
})


# ---------------------------------------------------------------------------
# RNG restoration
# ---------------------------------------------------------------------------

test_that("every seeded generator leaves the caller's RNG stream alone", {
  bnd <- .seed_sq()
  pts <- .seed_points(40, seed = 2)

  expect_true(.rng_untouched(
    get_voronoi_seeds(bnd, method = "random", n = 5, set_seed = 1)))
  expect_true(.rng_untouched(
    get_voronoi_seeds(bnd, method = "kmeans", n = 4, set_seed = 1)))
  expect_true(.rng_untouched(voronoi_seeds_kmeans(pts, k = 4, set_seed = 1)))
  expect_true(.rng_untouched(voronoi_seeds_random(bnd, k = 4, set_seed = 1)))

  # The guard on the guard: a generator that consumed RNG without restoring it
  # WOULD be caught -- set_seed = NULL is documented as "no seeding", and it
  # draws from the caller's stream.
  expect_false(.rng_untouched(
    get_voronoi_seeds(bnd, method = "random", n = 5, set_seed = NULL)))
})


test_that("seeding restores an absent .Random.seed rather than inventing one", {
  bnd <- .seed_sq()
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    rm(".Random.seed", envir = globalenv())

  invisible(voronoi_seeds_random(bnd, k = 3, set_seed = 1))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))

  # Restore a stream for the rest of the file.
  set.seed(1)
})


# ---------------------------------------------------------------------------
# .robust_st_sample()
# ---------------------------------------------------------------------------

test_that(".robust_st_sample returns exactly the requested number of points", {
  geom <- sf::st_geometry(.seed_sq())

  for (n in c(1L, 13L, 100L)) {
    set.seed(31)
    pts <- spatialkit:::.robust_st_sample(geom, n)
    expect_s3_class(pts, "sfc")
    expect_length(pts, n)
    expect_true(all(as.character(sf::st_geometry_type(pts)) == "POINT"))
  }

  # st_sample(type = "random") is only approximately exact for polygons, which
  # is the reason for the padding loop; the result must still be exact.
  set.seed(32)
  awkward <- sf::st_geometry(sf::st_buffer(
    sf::st_sfc(sf::st_point(c(0, 0)), crs = 32632), 1))
  expect_length(spatialkit:::.robust_st_sample(awkward, 25L), 25L)
})


# ---------------------------------------------------------------------------
# Every method returns seeds in the BOUNDARY's CRS
# ---------------------------------------------------------------------------

test_that("get_voronoi_seeds returns seeds in the boundary CRS for every method", {
  # NEWS states this for all three methods, and it is what makes the seeds
  # droppable straight into build_tessellation() alongside the boundary.  Two
  # of the branches happen to produce boundary-CRS seeds already (they sample
  # or align against the boundary), so only "provided" -- where the caller
  # hands over seeds carrying whatever CRS they were stored in -- shows the
  # difference.  All three are asserted because all three are promised.
  bnd <- .seed_sq()                                   # EPSG:32632
  expect_equal(sf::st_crs(bnd), sf::st_crs(32632))

  pts32 <- .seed_points(12, seed = 4)                 # same CRS ...
  pts_ll <- sf::st_transform(pts32, 4326)             # ... and a lon/lat copy
  expect_equal(sf::st_crs(pts_ll)$epsg, 4326L)

  cases <- list(
    provided = get_voronoi_seeds(bnd, method = "provided", seeds = pts_ll),
    random   = get_voronoi_seeds(bnd, method = "random", n = 5, set_seed = 1),
    kmeans   = get_voronoi_seeds(bnd, method = "kmeans", n = 4,
                                 sample_points = pts_ll, set_seed = 1),
    kmeans_internal = get_voronoi_seeds(bnd, method = "kmeans", n = 4,
                                        set_seed = 1)
  )
  for (nm in names(cases))
    expect_equal(sf::st_crs(cases[[nm]]), sf::st_crs(bnd), info = nm)

  # Aligning means REPROJECTING, not restamping: the provided seeds come back
  # at the metric coordinates the degrees denote, not at the degrees relabelled.
  got <- sf::st_coordinates(cases$provided)
  expect_equal(unname(got), unname(sf::st_coordinates(pts32)), tolerance = 1e-6)
  # ... and they are NOT the degrees they arrived as, merely relabelled.
  expect_false(isTRUE(all.equal(unname(got),
                                unname(sf::st_coordinates(pts_ll)),
                                tolerance = 1e-3)))

  # And with no boundary at all there is nothing to align to, so the seeds keep
  # their own CRS rather than silently losing it.
  expect_equal(sf::st_crs(get_voronoi_seeds(method = "provided",
                                            seeds = pts_ll)),
               sf::st_crs(4326))
})


# ---------------------------------------------------------------------------
# k-means seed counts are clamped, not fatal
# ---------------------------------------------------------------------------

test_that("k-means seeding clamps n to what k-means can actually produce", {
  # stats::kmeans() refuses BOTH k > distinct rows and k >= nrow(x), the latter
  # with "number of cluster centres must lie between 1 and nrow(x)".  So
  # n = nrow(sample_points) -- the obvious "one seed per point" request -- used
  # to die on a raw kmeans error rather than clamp, despite `@return`
  # promising "at most n".
  bnd   <- .seed_sq()
  cloud <- .seed_points(12, seed = 4)
  N     <- nrow(cloud)
  expect_equal(N, 12L)

  # Below the ceiling: exactly what was asked for, and no clamp message.
  quiet <- capture_spatialkit_log(
    under <- get_voronoi_seeds(bnd, method = "kmeans", n = N - 1L,
                               sample_points = cloud, set_seed = 1))
  expect_equal(nrow(under), N - 1L)
  expect_false(log_has(quiet, "clamping"))
  expect_equal(nrow(voronoi_seeds_kmeans(cloud, k = N - 1L, set_seed = 1)),
               N - 1L)

  # AT the ceiling and ABOVE it: clamped to nrow - 1, with a message naming the
  # count actually used.  No error either way.
  for (k in c(N, N + 3L)) {
    lines <- capture_spatialkit_log(
      at <- get_voronoi_seeds(bnd, method = "kmeans", n = k,
                              sample_points = cloud, set_seed = 1))
    expect_equal(nrow(at), N - 1L, info = paste("get_voronoi_seeds n =", k))
    expect_equal(at$seed_id, seq_len(N - 1L))
    expect_true(log_has(lines, "clamping"), info = paste("n =", k))

    lines2 <- capture_spatialkit_log(
      at2 <- voronoi_seeds_kmeans(cloud, k = k, set_seed = 1))
    expect_equal(nrow(at2), N - 1L, info = paste("voronoi_seeds_kmeans k =", k))
    expect_true(log_has(lines2, "clamping to 11"), info = paste("k =", k))
  }
})
