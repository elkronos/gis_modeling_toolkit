# tests/testthat/test-fix-resolution-bound-warning.R
# ---------------------------------------------------------------------------
# A count that a bound of the ladder chose reaches the caller of the
# tessellation functions as an R warning with a class, and so does
# `supported = FALSE`.  Both were log lines and a free-text attribute:
# tryCatch(), withCallingHandlers() and options(warn = 2) saw nothing.
# ---------------------------------------------------------------------------

fbw_catch <- function(expr) {
  ws <- list()
  value <- withCallingHandlers(suppressMessages(expr), warning = function(w) {
    ws[[length(ws) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = ws, classes = unlist(lapply(ws, class)),
       text = vapply(ws, conditionMessage, character(1)))
}
# A smooth trend with little noise: Cp keeps falling to the support ceiling.
fbw_pts <- function(n = 300, seed = 1) {
  set.seed(seed)
  d <- data.frame(x = 5e5 + runif(n, 0, 3000), y = 5e6 + runif(n, 0, 2000))
  d$z <- sin((d$x - 5e5) / 400) + rnorm(n, 0, 0.3)
  sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
}
fbw_bnd <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
  c(5e5, 5e6), c(5e5 + 3000, 5e6), c(5e5 + 3000, 5e6 + 2000), c(5e5, 5e6 + 2000),
  c(5e5, 5e6)))), crs = 32632))
fbw_bound <- function(out)
  out$warnings[vapply(out$warnings, inherits, logical(1), "spatialkit_bound_chosen")]


test_that("an optimum at the support ceiling is a classed warning in every function that takes the count", {
  d <- fbw_pts()
  p <- fbw_catch(resolution_profile(d, "z", sac = 400, n_levels = 6, nstart = 5,
                                    estimate = NULL, obs_per_cell = NULL,
                                    partition = "kmeans"))$value
  sel <- select_resolution(p, "cp")
  expect_true(sel$at_ceiling)
  expect_identical(sel$edge, "the support ceiling (n / min_cell_n)")

  calls <- list(
    function(x) get_voronoi_seeds(method = "kmeans", n = x, sample_points = d,
                                  partition = "kmeans"),
    function(x) voronoi_seeds_kmeans(d, k = x, partition = "kmeans"),
    function(x) build_tessellation(d, boundary = fbw_bnd, method = "hex",
                                   approx_n_cells = x, quiet = TRUE))
  for (f in calls) for (x in list(sel, p)) {
    out <- fbw_catch(f(x))
    w <- fbw_bound(out)
    expect_length(w, 1L)
    expect_identical(w[[1L]]$edge, sel$edge)
    expect_identical(w[[1L]]$criterion, "cp")
    expect_identical(w[[1L]]$n, sel$best)
    expect_match(conditionMessage(w[[1L]]),
                 sprintf("is %d cells, the optimum of `cp`, and that optimum is the support ceiling",
                         sel$best))
    expect_match(conditionMessage(w[[1L]]), "the bound chose the count, not the criterion")
    expect_match(conditionMessage(w[[1L]]), "Lower min_cell_n")
  }
  # Caught by class, and an error under warn = 2 semantics.
  expect_identical(tryCatch(voronoi_seeds_kmeans(d, k = sel, partition = "kmeans"),
                            spatialkit_bound_chosen = function(w) "caught"), "caught")
  expect_error(withCallingHandlers(
    get_voronoi_seeds(method = "kmeans", n = sel, sample_points = d, partition = "kmeans"),
    warning = function(w) stop("escalated: ", class(w)[1L])),
    "escalated: spatialkit_bound_chosen")
  # What chose the count rides on the resolved value as fields, too.
  r <- suppressWarnings(spatialkit:::.resolve_cell_count(sel, "n", "f"))
  expect_identical(r$edge, sel$edge)
  expect_identical(r$criterion, "cp")
  expect_true(r$at_ceiling)
})

test_that("an interior optimum, a plain number and a ranked vector raise no such warning", {
  d <- fbw_pts()
  p <- fbw_catch(resolution_profile(d, "z", sac = 400, n_levels = 6, nstart = 5))$value
  # An interior optimum, made by hand: the third level of the ladder.
  mid <- p
  mid$cp <- abs(seq_along(mid$levels) - 3) + 1
  sel <- select_resolution(mid, "cp")
  expect_identical(sel$best, mid$levels[3L])
  expect_true(is.na(sel$edge))
  for (x in list(sel, 12, c(9L, 12L, 6L))) {
    out <- fbw_catch(get_voronoi_seeds(method = "kmeans", n = x, sample_points = d,
                                       set_seed = 1))
    expect_length(fbw_bound(out), 0L)
  }
})

test_that("a reliability pick at the range floor is a bound choosing, and is said", {
  d <- fbw_pts(400)
  vm <- data.frame(model = c("Nug", "Exp"), psill = c(0.3, 0.7), range = c(0, 200),
                   stringsAsFactors = FALSE)
  sac <- structure(600, class = c("sac_range", "numeric"), variogram_model = vm,
                   crs = sf::st_crs(32632))
  p <- fbw_catch(resolution_profile(d, "z", sac = sac, n_levels = 8, nstart = 3))$value
  rel <- select_resolution(p, "reliability")
  expect_true(rel$at_floor)
  out <- fbw_catch(get_voronoi_seeds(method = "kmeans", n = rel, sample_points = d))
  w <- fbw_bound(out)
  expect_length(w, 1L)
  expect_match(conditionMessage(w[[1L]]), "range floor")
  expect_match(conditionMessage(w[[1L]]), "Fewer cells would be wider than the range")
})

test_that("supported = FALSE is a classed warning, not a log line alone", {
  set.seed(2)
  n <- 60
  d <- data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), z = rnorm(n))
  pts <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
  # A range of 120 on a square kilometre needs about 60 cells; 60 points at
  # nine a cell support 6.
  vm <- data.frame(model = c("Nug", "Exp"), psill = c(0.2, 1), range = c(0, 40),
                   stringsAsFactors = FALSE)
  sac <- structure(120, class = c("sac_range", "numeric"), variogram_model = vm,
                   crs = sf::st_crs(pts))
  lines <- capture_spatialkit_log(
    out <- fbw_catch(resolution_profile(pts, "z", sac = sac, n_levels = 5, nstart = 2)))
  expect_false(attr(out$value, "bounds")$supported)
  w <- out$warnings[vapply(out$warnings, inherits, logical(1), "spatialkit_not_supported")]
  expect_length(w, 1L)
  expect_match(conditionMessage(w[[1L]]), "cannot support a tessellation")
  expect_identical(w[[1L]]$reason, "floor above ceiling")
  # Still logged, once.
  expect_identical(sum(grepl("cannot support a tessellation", lines)), 1L)
  expect_identical(tryCatch(resolution_profile(pts, "z", sac = sac, n_levels = 5, nstart = 2),
                            spatialkit_not_supported = function(w) "caught"), "caught")
})
