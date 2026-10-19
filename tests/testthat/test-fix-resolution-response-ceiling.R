# tests/testthat/test-fix-resolution-response-ceiling.R
# ---------------------------------------------------------------------------
# The support ceiling of resolution_profile() counts the rows that carry a
# response.  It counted every point, so a layer with the response missing on
# three rows in four got four times the cells its scored rows support.
# ---------------------------------------------------------------------------

frc_field <- function(n, seed, a = 100, nugget = 0.3, side = 1000) {
  set.seed(seed)
  xy <- data.frame(x = runif(n, 0, side), y = runif(n, 0, side))
  D  <- as.matrix(stats::dist(xy))
  xy$z <- as.numeric(t(chol(exp(-D / a) + diag(1e-8, n))) %*% rnorm(n)) +
    rnorm(n, sd = sqrt(nugget))
  sf::st_as_sf(xy, coords = c("x", "y"), crs = 32632)
}

# Every warning an expression raises, with its classes, none of them left to
# reach the reporter.
frc_catch <- function(expr) {
  ws <- list()
  value <- withCallingHandlers(expr, warning = function(w) {
    ws[[length(ws) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = ws,
       classes = unlist(lapply(ws, class)),
       text = vapply(ws, conditionMessage, character(1)))
}

# A variogram given by hand, so the bounds do not depend on a fit.
frc_sac <- function(pts, range = 240, nugget = 0.3, psill = 1) {
  vm <- data.frame(model = c("Nug", "Exp"), psill = c(nugget, psill),
                   range = c(0, range / 3), stringsAsFactors = FALSE)
  structure(range, class = c("sac_range", "numeric"), variogram_model = vm,
            crs = sf::st_crs(pts))
}


test_that("the ceiling is floor(rows with a response / min_cell_n), and says so", {
  pts <- frc_field(900, seed = 1)
  sac <- frc_sac(pts)
  set.seed(11)
  miss <- sample(900, 675)
  na <- pts
  na$z[miss] <- NA

  out <- frc_catch(resolution_profile(na, "z", sac = sac, nstart = 3, estimate = NULL,
                                      obs_per_cell = NULL, partition = "kmeans"))
  p <- out$value
  b <- attr(p, "bounds")
  # 225 scored rows at min_cell_n = 9 support 25 cells; every point would
  # have given 100, at 2.25 scored rows per cell.
  expect_identical(b$ceiling, 25L)
  expect_identical(b$ceiling_from, "min_cell_n (rows with a response)")
  expect_identical(b$n, 900L)
  expect_identical(b$n_response, 225L)
  expect_identical(max(p$levels), 25L)
  expect_lte(select_resolution(p, "cp")$best, 25L)

  # Said by an R warning with a class, which carries both counts.
  expect_true("spatialkit_response_ceiling" %in% out$classes)
  w <- out$warnings[[which(vapply(out$warnings, inherits, logical(1),
                                  "spatialkit_response_ceiling"))[1L]]]
  expect_identical(w$n_response, 225L)
  expect_identical(w$ceiling, 25L)
  expect_match(conditionMessage(w), "225 of the 900 points carry a response")
  expect_match(conditionMessage(w), "25 cells, not the 100")

  # The same ceiling as the 225 rows passed alone.
  alone <- frc_catch(resolution_profile(na[-miss, ], "z", sac = sac, nstart = 3,
                                        estimate = NULL, obs_per_cell = NULL,
                                        partition = "kmeans"))
  expect_identical(attr(alone$value, "bounds")$ceiling, 25L)
  expect_identical(attr(alone$value, "bounds")$ceiling_from, "min_cell_n")
  expect_false("spatialkit_response_ceiling" %in% alone$classes)

  # The print names the rows the ceiling counts, and the edge label too.
  expect_output(print(p), "on 900 points \\(225 with a response\\)")
  expect_output(print(p), "ceiling 25 from min_cell_n = 9 on the 225 rows with a response")
  edge <- spatialkit:::.ladder_edge(25L, p$levels, p$levels, b)
  expect_identical(edge, "the support ceiling (rows with a response / min_cell_n)")
})

test_that("without a response, or with none on any row, the ceiling counts every point", {
  pts <- frc_field(300, seed = 2)
  geo <- frc_catch(resolution_profile(pts, n_levels = 5, nstart = 2))
  bg <- attr(geo$value, "bounds")
  expect_identical(bg$ceiling, 33L)
  expect_identical(bg$ceiling_from, "min_cell_n")
  expect_true(is.na(bg$n_response))
  expect_false("spatialkit_response_ceiling" %in% geo$classes)

  # No row has a response: nothing can be scored, and the geometry columns
  # keep the ladder a geometry-only profile has.
  none <- pts
  none$z <- NA_real_
  out <- frc_catch(resolution_profile(none, "z", n_levels = 5, nstart = 2))
  bn <- attr(out$value, "bounds")
  expect_identical(bn$n_response, 0L)
  expect_identical(bn$ceiling, 33L)
  expect_identical(out$value$levels, geo$value$levels)
  expect_true(all(is.na(out$value$cp)))
  expect_false("spatialkit_response_ceiling" %in% out$classes)
  expect_false("spatialkit_not_supported" %in% out$classes)
})

test_that("on a subsample the ceiling is held to two scored rows per cell", {
  set.seed(3)
  n <- 4000
  d <- data.frame(x = runif(n, 0, 2000), y = runif(n, 0, 2000), z = rnorm(n))
  d$z[sample(n, 3000)] <- NA
  pts <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
  out <- frc_catch(resolution_profile(pts, "z", sac = 300, sample_n = 600,
                                      n_levels = 5, nstart = 2))
  b <- attr(out$value, "bounds")
  # About 150 of the 600 subsample rows are scored; half the subsample (300)
  # left one scored row in two cells.
  expect_true(is.finite(b$n_scored) && b$n_scored < 600L)
  expect_lte(b$ceiling, floor(b$n_scored / 2))
  expect_identical(b$ceiling_from, "sample_n")
  expect_output(print(out$value),
                sprintf("from the %d scored rows of the 600-point subsample", b$n_scored))
  expect_identical(spatialkit:::.ladder_edge(b$ceiling, out$value$levels, out$value$levels, b),
                   "the subsample's ceiling (two scored subsample rows per cell; raise sample_n)")
  # With a response on every row the cap is half the subsample, as it was.
  set.seed(3)
  full <- sf::st_as_sf(data.frame(x = runif(n, 0, 2000), y = runif(n, 0, 2000), z = rnorm(n)),
                       coords = c("x", "y"), crs = 32632)
  bf <- attr(frc_catch(resolution_profile(full, "z", sac = 300, sample_n = 600,
                                          n_levels = 4, nstart = 1))$value, "bounds")
  expect_identical(bf$n_scored, 600L)
  expect_identical(bf$ceiling, 300L)
  expect_identical(bf$ceiling_from, "sample_n")
})

test_that("scored rows per occupied cell are reported, and a concentrated response is named", {
  # The response on the western quarter of the extent alone: n / L is then
  # not the scored rows per occupied cell, because three cells in four hold
  # none.
  set.seed(5)
  n <- 900
  d <- data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), z = rnorm(n))
  d$z[d$x >= 250] <- NA
  pts <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
  n_resp <- sum(is.finite(d$z))

  wide <- frc_catch(resolution_profile(pts, "z", sac = 300, levels = c(25, 100), nstart = 3))
  p <- wide$value
  expect_true(all(c("cells_resp", "cell_resp_mean") %in% names(p)))
  at100 <- p[p$levels == 100L, ]
  # About a quarter of the 100 cells hold a scored row, with about 9 each.
  expect_lt(at100$cells_resp, 45L)
  expect_gt(at100$cells_resp, 15L)
  expect_gt(at100$cell_resp_mean, 6)
  expect_lt(at100$cell_resp_mean, 13)
  expect_equal(p$cell_resp_mean, n_resp / p$cells_resp)

  # The default ladder stops at floor(n_resp / 9), which is conservative
  # here, and the warning says so and names the way round it.
  def <- frc_catch(resolution_profile(pts, "z", sac = 300, n_levels = 5, nstart = 2))
  expect_identical(attr(def$value, "bounds")$ceiling, as.integer(floor(n_resp / 9)))
  w <- def$text[vapply(def$warnings, inherits, logical(1), "spatialkit_response_ceiling")]
  expect_length(w, 1L)
  expect_match(w, "lower than the occupied cells need")
  expect_match(w, "`cell_resp_mean`", fixed = TRUE)
  expect_match(w, "a profile of the rows with a response alone")
  # A geometry-only profile leaves both columns NA.
  geo <- resolution_profile(pts, levels = c(5, 10), nstart = 1)
  expect_true(all(is.na(geo$cells_resp)) && all(is.na(geo$cell_resp_mean)))
})

test_that("the warning names the bound in force: sample_n on a subsample, not min_cell_n", {
  # 4,000 points, 1,000 scored, cells fitted to 600: both ceilings are
  # sample_n's (half the scored subsample rows, and half the subsample), and
  # the message used to read "not the 300 that every point at min_cell_n = 9
  # would give", where min_cell_n allows 111 and 444.
  set.seed(3)
  n <- 4000
  d <- data.frame(x = runif(n, 0, 2000), y = runif(n, 0, 2000), z = rnorm(n))
  d$z[sample(n, 3000)] <- NA
  pts <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
  out <- frc_catch(resolution_profile(pts, "z", sac = 300, sample_n = 600,
                                      n_levels = 5, nstart = 2))
  b <- attr(out$value, "bounds")
  w <- out$warnings[[which(vapply(out$warnings, inherits, logical(1),
                                  "spatialkit_response_ceiling"))[1L]]]
  msg <- conditionMessage(w)
  expect_match(msg, sprintf("%d of the 600 points of the subsample carry a response", b$n_scored))
  expect_match(msg, "(1000 of the layer's 4000)", fixed = TRUE)
  expect_match(msg, sprintf("the ceiling is %d cells, two scored subsample rows per cell",
                            b$ceiling))
  expect_match(msg, "not the 300 that half the subsample would give")
  expect_match(msg, "`sample_n` sets this ceiling, not `min_cell_n`", fixed = TRUE)
  expect_match(msg, "the layer's 1000 rows with a response support 111 cells at min_cell_n = 9")
  expect_false(grepl("every point at min_cell_n", msg))
  # `cell_resp_mean` counts subsample rows, and the message says by how much
  # the layer's cells hold more.
  expect_match(msg, "`cell_resp_mean` counts the scored rows of the subsample", fixed = TRUE)
  expect_match(msg, sprintf("about %.3g times as many", 1000 / b$n_scored), fixed = TRUE)
  expect_equal(out$value$cell_resp_mean, b$n_scored / out$value$cells_resp)
  expect_identical(w$ceiling_from, "sample_n")
  expect_identical(w$n_sample, 600L)
  expect_identical(w$n_scored, b$n_scored)
  expect_identical(w$n_response, 1000L)

  # A subsample that does not bind: min_cell_n is named, and the rows
  # `cell_resp_mean` counts are still said to be the subsample's.
  set.seed(8)
  n <- 3000
  d <- data.frame(x = runif(n, 0, 2000), y = runif(n, 0, 2000), z = rnorm(n))
  d$z[sample(n, 2700)] <- NA
  sub <- frc_catch(resolution_profile(sf::st_as_sf(d, coords = c("x", "y"), crs = 32632), "z",
                                      sac = 300, sample_n = 1500, n_levels = 4, nstart = 1))
  expect_identical(attr(sub$value, "bounds")$ceiling_from, "min_cell_n (rows with a response)")
  m2 <- sub$text[vapply(sub$warnings, inherits, logical(1), "spatialkit_response_ceiling")]
  expect_match(m2, "300 of the 3000 points carry a response, so the ceiling is 33 cells, not the 333 that every point at min_cell_n = 9 would give")
  expect_match(m2, "rows of the subsample; the cells of the layer hold about")
})

test_that("rows short of a predictor are not said to lack a response", {
  # The response is complete; a predictor is missing on half the rows.  The
  # warning read "150 of the 300 points carry a response".
  set.seed(5)
  n <- 300
  d <- data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), z = rnorm(n), w = rnorm(n))
  d$w[sample(n, 150)] <- NA
  pts <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
  out <- frc_catch(resolution_profile(pts, "z", "w", sac = 300, n_levels = 5, nstart = 2))
  b <- attr(out$value, "bounds")
  expect_identical(b$n_response, 150L)
  expect_true(b$predictors)
  expect_identical(b$ceiling, 16L)
  msg <- out$text[vapply(out$warnings, inherits, logical(1), "spatialkit_response_ceiling")]
  expect_length(msg, 1L)
  expect_match(msg, "150 of the 300 points carry a response and every predictor, so the ceiling is 16 cells")
  expect_match(msg, "built from the rows with a response and every predictor")
  expect_false(grepl("carry a response,", msg, fixed = TRUE))
  expect_output(print(out$value), "on 300 points \\(150 with a response and every predictor\\)")
  expect_output(print(out$value),
                "from min_cell_n = 9 on the 150 rows with a response and every predictor")
  # And the bound an optimum at that ceiling sits on is named the same way.
  expect_identical(spatialkit:::.ladder_edge(b$ceiling, out$value$levels, out$value$levels, b),
                   "the support ceiling (rows with a response and every predictor / min_cell_n)")
  # Without predictors the wording is the response's alone.
  expect_false(attr(resolution_profile(pts, n_levels = 3, nstart = 1), "bounds")$predictors)
})
