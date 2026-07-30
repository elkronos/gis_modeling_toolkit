# ---------------------------------------------------------------------------
# Regression tests for smaller fixes:
#  - .dedup_points() must stay feature-aligned for MULTIPOINT input
#  - .fallback_bandwidth() adaptive clamp must respect [10, 0.9n]
# ---------------------------------------------------------------------------

test_that(".dedup_points removes duplicate POINT features (fast path)", {
  g <- sf::st_sfc(
    sf::st_point(c(0, 0)),
    sf::st_point(c(1, 1)),
    sf::st_point(c(0, 0)),  # duplicate of feature 1
    crs = 32632
  )
  out <- spatialkit:::.dedup_points(g)
  expect_length(out, 2L)
})

test_that(".dedup_points handles MULTIPOINT without misaligning features", {
  mp  <- sf::st_multipoint(rbind(c(0, 0), c(1, 1), c(2, 2)))
  g <- sf::st_sfc(
    mp,                          # feature 1
    mp,                          # feature 2: duplicate of feature 1
    sf::st_point(c(5, 5)),       # feature 3
    sf::st_point(c(5, 5)),       # feature 4: duplicate of feature 3
    sf::st_point(c(9, 9)),       # feature 5: unique
    crs = 32632
  )
  # Before the fix, st_coordinates() returned one row per *vertex* (7 rows
  # for 5 features), so the dedup mask was misaligned with features.
  out <- spatialkit:::.dedup_points(g)
  expect_length(out, 3L)
  gtypes <- as.character(sf::st_geometry_type(out, by_geometry = TRUE))
  expect_equal(sort(gtypes), c("MULTIPOINT", "POINT", "POINT"))

  # sf-object variant must subset rows, not vertices
  sf_obj <- sf::st_sf(id = 1:5, geometry = g)
  out_sf <- spatialkit:::.dedup_points(sf_obj)
  expect_s3_class(out_sf, "sf")
  expect_equal(nrow(out_sf), 3L)
  expect_equal(out_sf$id, c(1L, 3L, 5L))
})

test_that("create_voronoi_polygons accepts MULTIPOINT input", {
  g <- sf::st_sfc(
    sf::st_multipoint(rbind(c(0, 0), c(10, 0))),
    sf::st_multipoint(rbind(c(0, 10), c(10, 10))),
    sf::st_point(c(5, 5)),
    crs = 32632
  )
  pts <- sf::st_sf(geometry = g)
  res <- suppressMessages(create_voronoi_polygons(pts, quiet = TRUE))
  expect_true(is.list(res))
  expect_s3_class(res$cells, "sf")
  expect_gt(nrow(res$cells), 0L)
})

test_that(".fallback_bandwidth adaptive clamp respects [10, 0.9n]", {
  skip_if_not_installed("sp")
  fb <- spatialkit:::.fallback_bandwidth

  make_sp <- function(n) {
    sp::SpatialPointsDataFrame(
      cbind(seq_len(n), seq_len(n)),
      data = data.frame(x = seq_len(n))
    )
  }

  # Tiny n: lower clamp of 10 applies (old code returned floor(0.9n) < 10)
  expect_equal(fb(make_sp(3), adaptive = TRUE), 10L)
  # Mid n: 0.9n binds below the 50-neighbour target
  expect_equal(fb(make_sp(30), adaptive = TRUE), 27L)
  # Large n: the ~50-neighbour target binds
  expect_equal(fb(make_sp(200), adaptive = TRUE), 50L)
})
