# tests/testthat/test-fix-resolution-usability.R
# ---------------------------------------------------------------------------
# Small usability faults of the count functions: a constant response got a
# message about cell sizes, determine_optimal_levels() refused a response
# without predictors, a range floor equal to the ceiling left a ladder of one
# level, and gstat's console print came through resolution_profile().
# ---------------------------------------------------------------------------

fus_catch <- function(expr) {
  ws <- list()
  value <- withCallingHandlers(expr, warning = function(w) {
    ws[[length(ws) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = ws, classes = unlist(lapply(ws, class)),
       text = vapply(ws, conditionMessage, character(1)))
}
fus_pts <- function(n, seed = 1) {
  set.seed(seed)
  sf::st_as_sf(data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), z = rnorm(n)),
               coords = c("x", "y"), crs = 32632)
}


test_that("a constant response is called constant, by the profile and by select_resolution()", {
  pts <- fus_pts(200)
  pts$z <- 5
  out <- fus_catch(resolution_profile(pts, "z", n_levels = 5, nstart = 2))
  p <- out$value
  expect_true("spatialkit_constant_response" %in% out$classes)
  expect_true(any(grepl("the response is constant over the 200 scored row", out$text)))
  expect_true(any(grepl("`cp`, `moran_i` and `moran_z` are NA at every level", out$text,
                        fixed = TRUE)))
  expect_true(attr(p, "one_cell")$constant)
  expect_true(all(is.na(p$cp)))
  expect_true(all(is.na(p$moran_i)) && all(is.na(p$moran_z)))
  # The message used to ask for "a level whose cells hold two scored rows",
  # which every level of this profile has.
  e <- tryCatch(select_resolution(p, "cp"), error = function(e) conditionMessage(e))
  expect_match(e, "`cp` is NA at every level because the response is constant")
  expect_false(grepl("two scored rows", e))
  # A response that varies is not called constant.
  v <- fus_catch(resolution_profile(fus_pts(200), "z", sac = 300, n_levels = 5, nstart = 2))
  expect_false("spatialkit_constant_response" %in% v$classes)
  expect_false(attr(v$value, "one_cell")$constant)
})

test_that("determine_optimal_levels() scores a response without predictors", {
  pts <- fus_pts(400, seed = 5)
  out <- fus_catch(determine_optimal_levels(pts, max_levels = 40, response_var = "z",
                                            criterion = "morans_i"))
  # It fell back to the geometric ranking with "needs both response_var and
  # predictor_vars", though the statistic runs on an intercept alone.
  expect_false(any(grepl("falling back to geometric", out$text)))
  # The ranking is by a diagnostic, and that comes with it as a classed
  # warning (see test-fix-resolution-moran-diagnostic.R).
  expect_true("spatialkit_diagnostic_ranking" %in% out$classes)
  d <- attr(out$value, "diagnostics")
  expect_false(is.null(d))
  expect_true(any(is.finite(d$moran_z)))
  expect_true(all(out$value %in% d$eval_ks[is.finite(d$moran_z[d$eval_ks])]))
  # A response alone does not change the default criterion.
  geo <- fus_catch(determine_optimal_levels(pts, max_levels = 40, response_var = "z"))
  expect_null(attr(geo$value, "diagnostics"))
  expect_identical(as.integer(geo$value),
                   as.integer(fus_catch(determine_optimal_levels(pts, max_levels = 40))$value))
  # Without a response there is still nothing to score, and that is said.
  none <- fus_catch(determine_optimal_levels(pts, max_levels = 40, criterion = "morans_i"))
  expect_true(any(grepl("criterion = 'morans_i' needs a response_var, and none was given",
                        none$text)))
  expect_null(attr(none$value, "diagnostics"))
})

test_that("a range floor equal to the ceiling does not leave a ladder of one level", {
  pts <- fus_pts(90)
  area <- as.numeric(sf::st_area(sf::st_convex_hull(sf::st_union(sf::st_geometry(pts)))))
  # 90 points support 10 cells; a range that puts the floor at 10 as well.
  r <- sqrt(area / 9.5)
  vm <- data.frame(model = c("Nug", "Exp"), psill = c(0.3, 1), range = c(0, r / 3),
                   stringsAsFactors = FALSE)
  sac <- structure(r, class = c("sac_range", "numeric"), variogram_model = vm,
                   crs = sf::st_crs(pts))
  for (rf in c(TRUE, FALSE)) {
    p <- suppressWarnings(resolution_profile(pts, "z", sac = sac, nstart = 2, range_floor = rf))
    b <- attr(p, "bounds")
    expect_identical(b$floor, 10L)
    expect_identical(b$ceiling, 10L)
    expect_true(b$supported)
    expect_identical(min(p$levels), 2L)
    expect_gt(nrow(p), 3L)
    expect_false(identical(select_resolution(p, "cp")$edge, "the only level scored"))
  }
})

test_that("resolution_profile() lets nothing through to the console from the variogram fit", {
  skip_if_not_installed("gstat")
  # A layer on which gstat's fit is singular and the semivariance is far from
  # the scale of the distances: fit.variogram() then print()s "a possible
  # solution MIGHT be to scale semivariances and/or distances" each time.
  set.seed(11)
  n <- 300
  sing <- sf::st_as_sf(data.frame(x = 5e5 + runif(n, 0, 1e-3), y = 5e6 + runif(n, 0, 1e-3),
                                  z = rnorm(n) * 1e6),
                       coords = c("x", "y"), crs = 32632)
  lines <- capture_spatialkit_log(
    printed <- utils::capture.output(
      p <- suppressWarnings(resolution_profile(sing, "z", n_levels = 3, nstart = 1))))
  expect_length(printed, 0L)
  expect_s3_class(p, "resolution_profile")
  # What was held back is in the log, when there was anything.
  if (any(grepl("printed", lines)))
    expect_true(any(grepl("held back here: .*scale semivariances", lines)))
})

test_that("a response that is constant to within rounding scores nothing, as the message says", {
  # 1e6 + 1e-9 x noise: called constant, and then given a finite `cp`
  # computed from its rounding error, and 2 cells by select_resolution().
  near <- fus_pts(200)
  near$z <- 1e6 + 1e-9 * near$z
  out <- fus_catch(resolution_profile(near, "z", n_levels = 5, nstart = 2))
  p <- out$value
  expect_true("spatialkit_constant_response" %in% out$classes)
  expect_true(attr(p, "one_cell")$constant)
  expect_true(all(is.na(p$cp)))
  expect_true(all(is.na(p$moran_i)) && all(is.na(p$moran_z)))
  expect_null(attr(p, "cp_noise"))
  expect_true(is.na(attr(p, "one_cell")$cp))
  # No noise variance is read off rounding error, and none is announced.
  expect_false(any(grepl("takes its noise variance", out$text)))
  expect_error(select_resolution(p, "cp"), "because the response is constant")
  expect_error(select_resolution(p, "moran_z", diagnostic = TRUE),
               "because the response is constant")

  # A response the predictors fit exactly leaves residuals that are rounding
  # error on the response's scale.  Against their own sum of squares they
  # never looked constant, and Cp took a noise variance of 7.9e-31.
  lin <- fus_pts(200)
  lin$w <- rnorm(200)
  lin$z <- 3 * lin$w + 2
  ol <- fus_catch(resolution_profile(lin, "z", "w", n_levels = 5, nstart = 2))
  expect_true("spatialkit_constant_response" %in% ol$classes)
  expect_true(any(grepl("the OLS residual of the response on `predictor_vars` is constant",
                        ol$text, fixed = TRUE)))
  expect_true(all(is.na(ol$value$cp)))
  expect_error(select_resolution(ol$value, "cp"),
               "because the OLS residual of the response on the predictors is constant")

  # A small scale is not a constant: 1e-9 x noise is scored as any response.
  tiny <- fus_pts(200)
  tiny$z <- 1e-9 * tiny$z
  ot <- fus_catch(resolution_profile(tiny, "z", sac = 300, n_levels = 5, nstart = 2))
  expect_false("spatialkit_constant_response" %in% ot$classes)
  expect_false(attr(ot$value, "one_cell")$constant)
  expect_true(all(is.finite(ot$value$cp)))
})
