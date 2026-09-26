# tests/testthat/test-review2-resolution.R
# ---------------------------------------------------------------------------
# Regressions from the second review of the resolution step:
# resolution_profile() and determine_optimal_levels().
# ---------------------------------------------------------------------------

# A hand-made sac_range, so the variogram-based columns exist without gstat.
r2_sac <- function(range, nugget = 1, psill = 1, crs = sf::st_crs(32632), ...) {
  vm <- data.frame(model = c("Nug", "Exp"), psill = c(nugget, psill),
                   range = c(0, range / 3), stringsAsFactors = FALSE)
  structure(range, class = c("sac_range", "numeric"), variogram_model = vm,
            crs = crs, ...)
}

r2_pts <- function(n = 300, seed = 1, ext = 1000, x0 = 5e5, y0 = 5e6, crs = 32632) {
  set.seed(seed)
  d <- data.frame(x = x0 + runif(n, 0, ext), y = y0 + runif(n, 0, ext), w = rnorm(n))
  d$z <- d$w + rnorm(n)
  sf::st_as_sf(d, coords = c("x", "y"), crs = crs)
}

# Every R warning an expression raises, muffled.
r2_warnings <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(expr, warning = function(c) {
    w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning")
  })
  list(value = val, warnings = w)
}


test_that("one empty point does not turn determine_optimal_levels() into 1", {
  set.seed(1)
  pts <- sf::st_as_sf(
    data.frame(x = 5e5 + c(runif(25, 0, 10), runif(25, 90, 100)),
               y = 5e6 + c(runif(25, 0, 10), runif(25, 90, 100))),
    coords = c("x", "y"), crs = 32632)
  clean <- determine_optimal_levels(pts, max_levels = 6)
  holed <- rbind(pts[1:10, ], sf::st_sf(geometry = sf::st_sfc(sf::st_point(), crs = 32632)),
                 pts[11:50, ])
  out <- r2_warnings(determine_optimal_levels(holed, max_levels = 6))
  expect_identical(out$value, clean)
  expect_true(any(grepl("dropping 1 point", out$warnings)))
  # Under a split the positions index the layer as passed: the empty row is
  # in neither half, and every other row is in one.
  set.seed(2)
  d <- data.frame(x = runif(80, 0, 1000), y = runif(80, 0, 1000), w = rnorm(80))
  d$z <- d$w + rnorm(80)
  p80 <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
  p81 <- rbind(p80[1:40, ], sf::st_sf(w = 0, z = 0,
                                      geometry = sf::st_sfc(sf::st_point(), crs = 32632)),
               p80[41:80, ])
  sp <- attr(suppressWarnings(determine_optimal_levels(p81, max_levels = 6,
                                                       select_on = "split")), "split")
  expect_setequal(c(sp$selection, sp$estimation), setdiff(seq_len(81), 41L))
})


test_that("a few missing predictor values do not take Moran's z away from determine_optimal_levels()", {
  # Three NAs in 400 rows used to make every cell mean holding one NA and the
  # whole window fall back to the geometric ranking.
  pts <- r2_pts(400, seed = 2)
  holed <- pts
  holed$w[c(5, 50, 150)] <- NA
  dol <- function(p) suppressWarnings(determine_optimal_levels(
    p, max_levels = 40, response_var = "z", predictor_vars = "w", criterion = "morans_i"))
  a <- dol(pts); b <- dol(holed)
  da <- attr(a, "diagnostics"); db <- attr(b, "diagnostics")
  expect_false(is.null(db))
  ks <- da$eval_ks[is.finite(da$moran_z[da$eval_ks])]
  expect_gt(length(ks), 0L)
  expect_true(all(is.finite(db$moran_z[ks])))
  expect_equal(db$moran_z[ks], da$moran_z[ks], tolerance = 0.2)
  lines <- capture_spatialkit_log(dol(holed))
  expect_true(log_has(lines, "3 of 400 row"))
})


test_that("a misspelt column is an error in determine_optimal_levels(), and a missing one a warning", {
  pts <- r2_pts(200)
  expect_error(determine_optimal_levels(pts, response_var = "Z", predictor_vars = "w"),
               "column 'Z' not found")
  expect_error(determine_optimal_levels(pts, response_var = "z", predictor_vars = c("w", "Elev")),
               "Elev.*not found")
  out <- r2_warnings(determine_optimal_levels(pts, max_levels = 6, response_var = "z",
                                              criterion = "morans_i"))
  expect_true(any(grepl("predictor_vars were not given", out$warnings)))
  expect_null(attr(out$value, "diagnostics"))
})
