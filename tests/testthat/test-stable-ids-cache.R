# ===========================================================================
# R/stable-ids-cache.R: ensure_stable_poly_id(), .cache_key(),
# create_grid_polygons_cached() and clear_grid_cache().
#
# "Stable" and "cached" are property claims -- order-invariance of the IDs,
# and identity on a cache hit -- that are cheap to check and were entirely
# unverified.
#
# Every cached call below passes its own `cache_env`, so nothing here touches
# the package-level .gmt_cache that a user session would share.
# ===========================================================================

.sic_sq <- function(x0, y0, x1, y1) {
  sf::st_polygon(list(rbind(c(x0, y0), c(x1, y0), c(x1, y1), c(x0, y1),
                            c(x0, y0))))
}

.sic_boundary <- function(side = 100, crs = 32632) {
  sf::st_sf(geometry = sf::st_sfc(.sic_sq(0, 0, side, side), crs = crs))
}

# Four unit squares on a 2 x 2 layout, deliberately NOT in sorted order.
.sic_polys <- function(crs = 32632) {
  sf::st_sf(v = 1:4,
            geometry = sf::st_sfc(.sic_sq(2, 2, 3, 3),   # v = 1, top-right
                                  .sic_sq(0, 0, 1, 1),   # v = 2, bottom-left
                                  .sic_sq(2, 0, 3, 1),   # v = 3, bottom-right
                                  .sic_sq(0, 2, 1, 3),   # v = 4, top-left
                                  crs = crs))
}


# ---------------------------------------------------------------------------
# ensure_stable_poly_id()
# ---------------------------------------------------------------------------

test_that("ensure_stable_poly_id assigns IDs by position, not by input order", {
  polys <- .sic_polys()
  # Sort in the native CRS so the expected order is the plain (x, y) one; the
  # default lon/lat sort key adds meridian convergence, which is correct but
  # not what this assertion is about.
  out <- ensure_stable_poly_id(polys, transform_for_sort = NULL)

  expect_s3_class(out, "sf")
  expect_equal(out$poly_id, 1:4)
  # Sorted by representative x, then y:
  #   (0,0) -> (0,2) -> (2,0) -> (2,2), i.e. v = 2, 4, 3, 1
  expect_equal(out$v, c(2L, 4L, 3L, 1L))

  # The claim in the name: shuffle the rows and the SAME polygon keeps the
  # SAME id.  This is what makes a cached grid joinable across sessions.
  set.seed(1)
  for (i in 1:5) {
    shuffled <- polys[sample(nrow(polys)), , drop = FALSE]
    out_s <- ensure_stable_poly_id(shuffled, transform_for_sort = NULL)
    expect_equal(out_s$v, out$v)
    expect_equal(out_s$poly_id, out$poly_id)
    expect_equal(sf::st_coordinates(out_s)[, 1:2],
                 sf::st_coordinates(out)[, 1:2])
  }

  # The default lon/lat sort key is stable too, just possibly a different
  # permutation.
  d1 <- ensure_stable_poly_id(polys)
  d2 <- ensure_stable_poly_id(polys[c(3L, 1L, 4L, 2L), , drop = FALSE])
  expect_equal(d1$v, d2$v)
  expect_equal(d1$poly_id, 1:4)
})


test_that("ensure_stable_poly_id honours id_col, method and sfc input", {
  polys <- .sic_polys()

  named <- ensure_stable_poly_id(polys, id_col = "cell_id",
                                 transform_for_sort = NULL)
  expect_true("cell_id" %in% names(named))
  expect_false("poly_id" %in% names(named))
  expect_equal(named$cell_id, 1:4)

  # All three representative-point strategies agree for axis-aligned squares,
  # where centroid, interior point and bbox centre coincide.
  ref <- ensure_stable_poly_id(polys, transform_for_sort = NULL)$v
  for (m in c("centroid", "surface_point", "bbox_center")) {
    got <- suppressWarnings(
      ensure_stable_poly_id(polys, method = m, transform_for_sort = NULL))
    expect_equal(got$v, ref, info = m)
  }
  expect_error(ensure_stable_poly_id(polys, method = "nope"),
               "'arg' should be one of")

  # An sfc is promoted to sf and given the id column.
  from_sfc <- ensure_stable_poly_id(sf::st_geometry(polys),
                                    transform_for_sort = NULL)
  expect_s3_class(from_sfc, "sf")
  expect_equal(from_sfc$poly_id, 1:4)
})


test_that("ensure_stable_poly_id drops non-polygon rows and refuses an empty result", {
  polys <- .sic_polys()
  mixed <- rbind(polys,
                 sf::st_sf(v = 5L, geometry = sf::st_sfc(sf::st_point(c(9, 9)),
                                                         crs = 32632)))

  expect_warning(lines <- capture_spatialkit_log(
    out <- ensure_stable_poly_id(mixed, transform_for_sort = NULL)),
    "dropping 1 non-polygon row")
  expect_true(log_has(lines, "dropping 1 non-polygon row\\(s\\) \\(POINT: 1\\)"))
  expect_equal(nrow(out), 4L)
  expect_equal(out$poly_id, 1:4)
  expect_false(5L %in% out$v)

  points_only <- sf::st_sf(v = 1L,
                           geometry = sf::st_sfc(sf::st_point(c(1, 1)),
                                                 crs = 32632))
  expect_error(suppressWarnings(ensure_stable_poly_id(points_only)),
               "no polygon rows found")
  expect_error(ensure_stable_poly_id(data.frame(x = 1)),
               "must be an sf/sfc object")
})


# ---------------------------------------------------------------------------
# .cache_key()
# ---------------------------------------------------------------------------

test_that(".cache_key separates keys that must not collide", {
  ck   <- spatialkit:::.cache_key
  bnd  <- .sic_boundary(100)
  bnd2 <- .sic_boundary(50)

  # A key is a single string, always -- exists() below needs exactly one.
  key <- ck(bnd, "square", 25)
  expect_type(key, "character")
  expect_length(key, 1L)
  expect_true(nzchar(key))
  expect_identical(key, ck(bnd, "square", 25))          # deterministic

  # The regression: as.integer() TRUNCATED target_cells, so 25.2 and 25.7
  # hashed to one key and the second call silently got the first one's grid.
  expect_false(ck(bnd, "square", 25.2) == ck(bnd, "square", 25.7))
  expect_false(ck(bnd, "square", 25.2) == ck(bnd, "square", 25))
  expect_identical(ck(bnd, "square", 25.2), ck(bnd, "square", 25.2))

  # A NULL target_cells used to collapse paste0() to character(0), which
  # crashes exists() downstream.  It now hashes like any other value.
  null_key <- ck(bnd, "square", NULL)
  expect_length(null_key, 1L)
  expect_false(null_key == ck(bnd, "square", 25))

  # Type, geometry and CRS all participate.
  expect_false(ck(bnd, "square", 25) == ck(bnd, "hex", 25))
  expect_false(ck(bnd, "square", 25) == ck(bnd2, "square", 25))
  expect_false(ck(bnd, "square", 25) ==
                 ck(sf::st_transform(bnd, 3857), "square", 25))

  # Extra grid arguments matter, but their ORDER does not -- the same call
  # written two ways must hit the same cache entry.
  expect_false(ck(bnd, "square", 25, clip = TRUE) ==
                 ck(bnd, "square", 25, clip = FALSE))
  expect_identical(ck(bnd, "square", 25, clip = TRUE, cellsize = 2),
                   ck(bnd, "square", 25, cellsize = 2, clip = TRUE))
})


# ---------------------------------------------------------------------------
# create_grid_polygons_cached() / clear_grid_cache()
# ---------------------------------------------------------------------------

test_that("create_grid_polygons_cached defaults to the same type as create_grid_polygons", {
  # The two defaults disagreed -- c("hex", "square") here against
  # c("square", "hex") there -- so the same call returned a DIFFERENT
  # tessellation depending on which one you reached for.
  bnd <- .sic_boundary(100)
  env <- new.env(parent = emptyenv())

  cached <- create_grid_polygons_cached(bnd, target_cells = 16,
                                        cache_env = env)
  direct <- create_grid_polygons(bnd, target_cells = 16, quiet = TRUE)

  expect_equal(nrow(cached), nrow(direct))
  expect_equal(nrow(cached), 16L)                       # a 4 x 4 square grid
  expect_equal(sum(as.numeric(sf::st_area(cached))),
               sum(as.numeric(sf::st_area(direct))), tolerance = 1e-6)

  # An explicit type = "square" is the same thing; "hex" is visibly not.
  expect_equal(nrow(create_grid_polygons_cached(bnd, target_cells = 16,
                                                type = "square",
                                                cache_env = env)), 16L)
  hex <- create_grid_polygons_cached(bnd, target_cells = 16, type = "hex",
                                     cache_env = env)
  expect_false(nrow(hex) == nrow(cached))

  # The cached result always carries a stable poly_id.
  expect_true("poly_id" %in% names(cached))
  expect_equal(sort(cached$poly_id), seq_len(nrow(cached)))
})


test_that("create_grid_polygons_cached returns the stored object on a hit", {
  bnd <- .sic_boundary(100)
  env <- new.env(parent = emptyenv())

  first <- create_grid_polygons_cached(bnd, target_cells = 9, cache_env = env)
  expect_equal(length(ls(env, all.names = TRUE)), 1L)

  # Tag the stored value.  A genuine cache hit hands back the TAGGED object; a
  # recomputation (which is deterministic here, so identical() alone could not
  # tell the difference) would not carry the tag.
  key <- ls(env, all.names = TRUE)[1L]
  tagged <- get(key, envir = env)
  tagged$..cache_probe <- seq_len(nrow(tagged))
  assign(key, tagged, envir = env)

  second <- create_grid_polygons_cached(bnd, target_cells = 9, cache_env = env)
  expect_true("..cache_probe" %in% names(second))
  expect_equal(length(ls(env, all.names = TRUE)), 1L)   # no second entry

  # A different target_cells is a different entry, not a hit.
  create_grid_polygons_cached(bnd, target_cells = 16, cache_env = env)
  expect_equal(length(ls(env, all.names = TRUE)), 2L)
})


test_that("create_grid_polygons_cached keys nearby fractional targets apart", {
  # The user-visible form of the .cache_key() truncation bug: two calls, two
  # entries, and the second call must not be served the first one's grid.
  bnd <- .sic_boundary(100)
  env <- new.env(parent = emptyenv())

  a <- create_grid_polygons_cached(bnd, target_cells = 25.2, cache_env = env)
  key_a <- ls(env, all.names = TRUE)[1L]
  probe <- get(key_a, envir = env)
  probe$..from <- "25.2"
  assign(key_a, probe, envir = env)

  b <- create_grid_polygons_cached(bnd, target_cells = 25.7, cache_env = env)
  expect_equal(length(ls(env, all.names = TRUE)), 2L)
  expect_false("..from" %in% names(b))                  # not the 25.2 entry

  # A NULL target_cells (grid sized by cellsize instead) must not crash the
  # exists() lookup.
  n1 <- create_grid_polygons_cached(bnd, target_cells = NULL, cellsize = 25,
                                    cache_env = env)
  expect_equal(nrow(n1), 16L)
  expect_equal(length(ls(env, all.names = TRUE)), 3L)
  n2 <- create_grid_polygons_cached(bnd, target_cells = NULL, cellsize = 25,
                                    cache_env = env)
  expect_equal(length(ls(env, all.names = TRUE)), 3L)   # a hit, not a fourth
  expect_equal(nrow(n2), nrow(n1))
})


test_that("clear_grid_cache empties the environment and reports the count", {
  bnd <- .sic_boundary(100)
  env <- new.env(parent = emptyenv())

  for (tc in c(4, 9, 16))
    create_grid_polygons_cached(bnd, target_cells = tc, cache_env = env)
  expect_equal(length(ls(env, all.names = TRUE)), 3L)

  removed <- clear_grid_cache(env)
  expect_equal(removed, 3L)
  expect_equal(length(ls(env, all.names = TRUE)), 0L)

  # Idempotent, and it returns invisibly.
  expect_equal(clear_grid_cache(env), 0L)
  expect_invisible(clear_grid_cache(env))
})


test_that("create_grid_polygons_cached validates the boundary", {
  env <- new.env(parent = emptyenv())
  expect_error(create_grid_polygons_cached(data.frame(x = 1), target_cells = 9,
                                           cache_env = env),
               "must be sf/sfc POLYGON/MULTIPOLYGON")
  expect_error(create_grid_polygons_cached(.sic_boundary(), target_cells = 9,
                                           type = "triangles",
                                           cache_env = env),
               "'arg' should be one of")

  # An sfc boundary is accepted and produces the same grid as the sf form.
  bnd <- .sic_boundary(100)
  from_sfc <- create_grid_polygons_cached(sf::st_geometry(bnd),
                                          target_cells = 9, cache_env = env)
  from_sf  <- create_grid_polygons_cached(bnd, target_cells = 9,
                                          cache_env = env)
  expect_equal(nrow(from_sfc), nrow(from_sf))
})
