# tests/testthat/test-cell-count.R
# ---------------------------------------------------------------------------
# The level-selection step's answer is accepted where a cell count is
# required: build_tessellation(approx_n_cells = ) and get_voronoi_seeds(n = )
# take a determine_optimal_levels() vector, a select_resolution() result or a
# resolution_profile(), and record on the output which count was used and
# what chose it.
# ---------------------------------------------------------------------------

cc_points <- function(n = 300, seed = 5) {
  set.seed(seed)
  cx <- runif(8, 100, 900); cy <- runif(8, 100, 900)
  i  <- sample(8, n, replace = TRUE)
  sf::st_as_sf(data.frame(x = cx[i] + rnorm(n, 0, 40), y = cy[i] + rnorm(n, 0, 40),
                          z = rnorm(n)),
               coords = c("x", "y"), crs = 32632)
}
cc_bnd <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
  c(0, 0), c(1000, 0), c(1000, 1000), c(0, 1000), c(0, 0)))), crs = 32632))
resolve <- spatialkit:::.resolve_cell_count


test_that("a plain number passes through with no provenance", {
  r <- resolve(12, "n", "f")
  expect_equal(r$n, 12); expect_null(r$from)
  expect_null(resolve(NULL, "n", "f"))
})

test_that("a determine_optimal_levels() vector is read at its first candidate", {
  pts <- cc_points()
  lv  <- determine_optimal_levels(pts, max_levels = 12, top_n = 3, set_seed = 1)
  expect_true(is.integer(lv) && length(lv) >= 1L)
  lines <- capture_spatialkit_log(r <- resolve(lv, "approx_n_cells", "build_tessellation"),
                                  level = logger::INFO)
  expect_equal(r$n, as.numeric(lv[1]))
  if (length(lv) > 1L) {
    expect_match(r$from, sprintf("first of %d ranked candidates", length(lv)))
    expect_true(log_has(lines, "using the first"))
  }
})

test_that("a select_resolution() result and a profile are read at their best level", {
  skip_if_not_installed("gstat")
  pts <- cc_points()
  prof <- suppressWarnings(resolution_profile(pts, response_var = "z", n_levels = 8, seed = 1))
  sel  <- select_resolution(prof, criterion = "elbow")
  r1 <- resolve(sel, "n", "get_voronoi_seeds")
  expect_equal(r1$n, as.numeric(sel$best))
  expect_match(r1$from, "select_resolution\\(criterion = \"elbow\"\\)")
  # The profile itself goes through select_resolution()'s default criterion
  # (cp needs a response -- present here -- and a usable variogram).
  r2 <- tryCatch(resolve(prof, "n", "get_voronoi_seeds"), error = function(e) e)
  if (!inherits(r2, "error")) {
    expect_equal(r2$n, as.numeric(select_resolution(prof)$best))
    expect_match(r2$from, "resolution_profile\\(\\) read with select_resolution")
  } else {
    expect_match(conditionMessage(r2), "NA at every level")
  }
})

test_that("build_tessellation() and get_voronoi_seeds() take the vector and record its source", {
  pts <- cc_points()
  lv  <- c(9L, 12L, 6L)
  t_hex <- build_tessellation(pts, boundary = cc_bnd, method = "hex",
                              approx_n_cells = lv, quiet = TRUE)
  expect_equal(t_hex$params$approx_n_cells, 9)
  expect_match(t_hex$params$approx_n_cells_from, "first of 3 ranked candidates \\(9, 12, 6\\)")
  ref <- build_tessellation(pts, boundary = cc_bnd, method = "hex",
                            approx_n_cells = 9, quiet = TRUE)
  expect_equal(nrow(t_hex$cells), nrow(ref$cells))
  expect_null(ref$params$approx_n_cells_from)

  seeds <- get_voronoi_seeds(cc_bnd, method = "kmeans", n = lv, sample_points = pts,
                             set_seed = 1)
  expect_equal(nrow(seeds), 9L)
  expect_match(attr(seeds, "n_from"), "first of 3 ranked candidates")
  ref_s <- get_voronoi_seeds(cc_bnd, method = "kmeans", n = 9, sample_points = pts,
                             set_seed = 1)
  expect_null(attr(ref_s, "n_from"))
  expect_equal(sf::st_coordinates(seeds), sf::st_coordinates(ref_s))

  # The pipeline end to end, no number carried by hand.
  lv2   <- determine_optimal_levels(pts, max_levels = 12, top_n = 3, set_seed = 1)
  seeds2 <- get_voronoi_seeds(cc_bnd, method = "kmeans", n = lv2, sample_points = pts,
                              set_seed = 1)
  expect_equal(nrow(seeds2), as.integer(lv2[1]))
  vor <- build_tessellation(seeds2, boundary = cc_bnd, method = "voronoi", quiet = TRUE)
  expect_equal(nrow(vor$cells), as.integer(lv2[1]))
})

test_that("what is not a count is refused with a message naming the shapes accepted", {
  pts <- cc_points()
  expect_error(build_tessellation(pts, boundary = cc_bnd, method = "hex",
                                  approx_n_cells = "9", quiet = TRUE),
               "must be a number, the integer vector determine_optimal_levels\\(\\) returns")
  expect_error(get_voronoi_seeds(cc_bnd, method = "random", n = 0, set_seed = 1),
               "must resolve to a positive number of cells; got 0")
  expect_error(get_voronoi_seeds(cc_bnd, method = "random", n = list(a = 1), set_seed = 1),
               "got an object of class list")
})
