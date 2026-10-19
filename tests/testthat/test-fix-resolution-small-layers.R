# tests/testthat/test-fix-resolution-small-layers.R
# ---------------------------------------------------------------------------
# A layer too small for a ladder of cell counts is not reported as supported:
# resolution_profile() used to raise floor(n / min_cell_n) = 0 or 1 to a
# ceiling of 2 "from min_cell_n", supported = TRUE, with nothing said, and
# determine_optimal_levels() returned 1 for two points and for none.
# ---------------------------------------------------------------------------

fsl_catch <- function(expr) {
  ws <- list()
  value <- withCallingHandlers(expr, warning = function(w) {
    ws[[length(ws) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = ws, classes = unlist(lapply(ws, class)),
       text = vapply(ws, conditionMessage, character(1)))
}
fsl_pts <- function(n, seed = 1) {
  set.seed(seed)
  sf::st_as_sf(data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), z = rnorm(n)),
               coords = c("x", "y"), crs = 32632)
}


test_that("fewer than two cells at min_cell_n is supported = FALSE, with a classed warning", {
  for (n in c(3L, 10L, 17L)) {
    out <- fsl_catch(resolution_profile(fsl_pts(n), "z", nstart = 2, estimate = NULL,
                                        obs_per_cell = NULL, sac = NULL, partition = "kmeans"))
    b <- attr(out$value, "bounds")
    expect_false(b$supported)
    expect_identical(b$ceiling_from, "fewer than two supported cells")
    expect_identical(b$ceiling, 2L)
    expect_identical(out$value$levels, 2L)       # the one level is still shown
    w <- out$warnings[vapply(out$warnings, inherits, logical(1), "spatialkit_not_supported")]
    expect_length(w, 1L)
    expect_match(conditionMessage(w[[1L]]),
                 sprintf("%d points at min_cell_n = 9 support %d cell", n, n %/% 9L))
    expect_match(conditionMessage(w[[1L]]), "`supported` is FALSE")
    expect_identical(w[[1L]]$reason, "fewer than two supported cells")
    expect_output(print(out$value), "fewer than two cells meet min_cell_n = 9")
    expect_output(print(out$value), "no ladder: not supported")
  }
  # Eighteen points are two cells of nine: a ladder of one level, supported.
  ok <- fsl_catch(resolution_profile(fsl_pts(18L), "z", nstart = 2, estimate = NULL,
                                     obs_per_cell = NULL, sac = NULL, partition = "kmeans"))
  expect_true(attr(ok$value, "bounds")$supported)
  expect_identical(attr(ok$value, "bounds")$ceiling_from, "min_cell_n")
  expect_false("spatialkit_not_supported" %in% ok$classes)
})

test_that("a min_cell_n above half the points is the same case, with or without a response", {
  pts <- fsl_pts(200L)
  for (resp in list(NULL, "z")) {
    out <- fsl_catch(resolution_profile(pts, resp, min_cell_n = 150, nstart = 2,
                                        estimate = NULL, obs_per_cell = NULL, sac = NULL,
                                        partition = "kmeans"))
    b <- attr(out$value, "bounds")
    expect_false(b$supported)
    expect_identical(b$ceiling_from, "fewer than two supported cells")
    expect_true("spatialkit_not_supported" %in% out$classes)
    expect_true(any(grepl("200 points at min_cell_n = 150 support 1 cell, and it takes two",
                          out$text)))
    expect_true(any(grepl("Lower `min_cell_n` (at most 100 here)", out$text, fixed = TRUE)))
  }
  # The rows with a response are what is counted: 200 points, 12 scored.
  few <- pts
  few$z[13:200] <- NA
  out <- fsl_catch(resolution_profile(few, "z", nstart = 2, estimate = NULL,
                                      obs_per_cell = NULL, sac = NULL, partition = "kmeans"))
  expect_false(attr(out$value, "bounds")$supported)
  expect_true(any(grepl("12 rows with a response at min_cell_n = 9 support 1 cell", out$text)))
})

test_that("points on one location have no ladder either", {
  one <- sf::st_as_sf(data.frame(x = rep(500, 30), y = rep(500, 30), z = rnorm(30)),
                      coords = c("x", "y"), crs = 32632)
  out <- fsl_catch(resolution_profile(one, nstart = 2))
  b <- attr(out$value, "bounds")
  expect_false(b$supported)
  expect_identical(b$ceiling_from, "fewer than two distinct locations")
  expect_true(any(grepl("the 30 points sit on 1 distinct location", out$text)))
  expect_true("spatialkit_not_supported" %in% out$classes)
})

test_that("determine_optimal_levels() does not hand back 1 for a layer with no ladder", {
  pts <- fsl_pts(40L)
  # No rows: an error.  It returned 1.
  expect_error(determine_optimal_levels(pts[0, ]), "has no rows")
  # One or two points: 1, said with the class every count that is not a
  # finding carries, or NA, or an error, as asked.
  for (n in 1:2) {
    out <- fsl_catch(determine_optimal_levels(pts[seq_len(n), ]))
    expect_identical(as.integer(out$value), 1L)
    expect_true("spatialkit_no_elbow" %in% out$classes)
    expect_match(out$text[1L], "too few for a ladder of cell counts")
    expect_true(is.na(determine_optimal_levels(pts[seq_len(n), ], on_no_elbow = "na")))
    expect_error(determine_optimal_levels(pts[seq_len(n), ], on_no_elbow = "error"),
                 class = "spatialkit_no_elbow")
  }
  # Every point on one location.
  one <- sf::st_as_sf(data.frame(x = rep(500, 12), y = rep(500, 12)),
                      coords = c("x", "y"), crs = 32632)
  out <- fsl_catch(determine_optimal_levels(one))
  expect_identical(as.integer(out$value), 1L)
  expect_true("spatialkit_no_elbow" %in% out$classes)
  expect_match(out$text[1L], "sit on 1 distinct location")
})

test_that("the too-few warning gives advice that can be followed", {
  # One row with a response of 200: "1 rows ... support 0 cells ... Lower
  # `min_cell_n` (at most 0 here)", and min_cell_n cannot go below 1.
  one <- fsl_pts(200L)
  one$z[-1] <- NA
  out <- fsl_catch(resolution_profile(one, "z", n_levels = 5, nstart = 2))
  msg <- out$text[vapply(out$warnings, inherits, logical(1), "spatialkit_not_supported")]
  expect_length(msg, 1L)
  expect_match(msg, "1 row with a response at min_cell_n = 9 supports 0 cells, and it takes two")
  expect_match(msg, "Two cells take a response on at least two rows, whatever `min_cell_n` is",
               fixed = TRUE)
  expect_false(grepl("at most 0", msg))
  expect_false(grepl("1 rows", msg))
  # Two rows: a min_cell_n of 1 would give two cells, and that is what is said.
  two <- fsl_pts(200L)
  two$z[-(1:2)] <- NA
  m2 <- fsl_catch(resolution_profile(two, "z", n_levels = 5, nstart = 2))$text
  expect_true(any(grepl("2 rows with a response at min_cell_n = 9 support 0 cells", m2)))
  expect_true(any(grepl("Lower `min_cell_n` (at most 1 here), or add rows with a response.",
                        m2, fixed = TRUE)))
})

test_that("levels the caller asked for are said to be scored as given", {
  # Ten points, levels = 2:4: the warning said "the profile scores 2 cells"
  # over a table of three levels, and the print "no ladder".
  out <- fsl_catch(resolution_profile(fsl_pts(10L), "z", levels = 2:4, nstart = 2))
  p <- out$value
  expect_identical(p$levels, 2:4)
  expect_false(attr(p, "bounds")$supported)
  msg <- out$text[vapply(out$warnings, inherits, logical(1), "spatialkit_not_supported")]
  expect_match(msg, "10 points at min_cell_n = 9 support 1 cell, and it takes two, so no cell count meets `min_cell_n`",
               fixed = TRUE)
  expect_match(msg, "The 3 levels asked for in `levels` (2 to 4 cells) are scored as given",
               fixed = TRUE)
  expect_false(grepl("The profile scores 2 cells", msg))
  expect_false(grepl("no ladder", msg))
  shown <- utils::capture.output(print(p))
  expect_true(any(grepl("levels as given: not supported", shown)))
  expect_false(any(grepl("no ladder", shown)))
  # One level asked for.
  m1 <- fsl_catch(resolution_profile(fsl_pts(10L), "z", levels = 3, nstart = 2))$text
  expect_true(any(grepl("The 1 level asked for in `levels` (3 cells) is scored as given",
                        m1, fixed = TRUE)))
  # A min_cell_n no level can meet, with levels given.
  big <- fsl_catch(resolution_profile(fsl_pts(200L), "z", min_cell_n = 150,
                                      levels = c(2, 4, 8, 12, 16), nstart = 2))
  expect_true(any(grepl("The 5 levels asked for in `levels` (2 to 16 cells) are scored as given",
                        big$text, fixed = TRUE)))
  expect_output(print(big$value), "levels as given: not supported")
})
