# tests/testthat/test-fix-resolution-no-structure.R
# ---------------------------------------------------------------------------
# A count that is not a finding about the data reaches the caller as a flag
# and a classed condition, not as warning text alone: a Cp optimum that one
# cell for the whole layer matches (structured = FALSE), and a
# determine_optimal_levels() count read off a WSS curve with no elbow.
# ---------------------------------------------------------------------------

fns_catch <- function(expr) {
  ws <- list()
  value <- withCallingHandlers(expr, warning = function(w) {
    ws[[length(ws) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = ws, classes = unlist(lapply(ws, class)),
       text = vapply(ws, conditionMessage, character(1)))
}
fns_noise <- function(n = 600, seed = 7) {
  set.seed(seed)
  sf::st_as_sf(data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), z = rnorm(n)),
               coords = c("x", "y"), crs = 32632)
}
fns_field <- function(n = 400, seed = 2, a = 100, nugget = 0.3) {
  set.seed(seed)
  xy <- data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000))
  D  <- as.matrix(stats::dist(xy))
  xy$z <- as.numeric(t(chol(exp(-D / a) + diag(1e-8, n))) %*% rnorm(n)) +
    rnorm(n, sd = sqrt(nugget))
  sf::st_as_sf(xy, coords = c("x", "y"), crs = 32632)
}
fns_sac <- function(pts, range = 240, nugget = 0.3, psill = 1, value = range, ...) {
  vm <- data.frame(model = c("Nug", "Exp"), psill = c(nugget, psill),
                   range = c(0, range / 3), stringsAsFactors = FALSE)
  structure(value, class = c("sac_range", "numeric"), variogram_model = vm,
            crs = sf::st_crs(pts), ...)
}


test_that("the profile keeps what one cell scores, the level every other has to beat", {
  pts <- fns_noise()
  p <- fns_catch(resolution_profile(pts, "z", sac = 300, n_levels = 8, nstart = 3))$value
  one <- attr(p, "one_cell")
  z <- pts$z
  expect_identical(one$n, 600L)
  expect_equal(one$rss, sum((z - mean(z))^2))
  # Cp by the formula of the levels, at L = 1, with their noise variance.
  tau2 <- attr(p, "cp_noise")$value
  expect_equal(one$cp, one$rss / 600 + tau2 * 2 * 1 / 600)
  expect_false(one$constant)
  expect_gt(one$rss, max(p$rss))
  expect_output(print(p), "one cell    : rss")
  # No response, no reference.
  expect_null(attr(resolution_profile(pts, n_levels = 3, nstart = 1), "one_cell"))
})

test_that("select_resolution() flags a Cp optimum that one cell matches", {
  # White noise: no level predicts better than one mean for the layer.
  pts <- fns_noise()
  p <- fns_catch(resolution_profile(pts, "z", sac = 300, n_levels = 8, nstart = 3,
                                    estimate = NULL, obs_per_cell = NULL,
                                    partition = "kmeans"))$value
  s <- select_resolution(p, "cp")
  expect_identical(s$structured, FALSE)
  expect_equal(s$one_cell, attr(p, "one_cell")$cp)
  expect_lte(s$one_cell, s$value * 1.02)
  expect_match(s$note, "Cp at one cell")
  expect_match(s$note, "the count is not identified")
  expect_output(print(s), "not identified")
  expect_identical(s$best, p$levels[which.min(p$cp)])       # the pick itself is unchanged
  expect_identical(summary(p)$structured[summary(p)$criterion == "cp"], FALSE)
  expect_output(print(summary(p)), "one cell for the whole layer does as well")

  # A field with structure: the best level beats one cell by far more than
  # the tolerance.
  f <- fns_field()
  pf <- fns_catch(resolution_profile(f, "z", sac = fns_sac(f), n_levels = 8, nstart = 3,
                                     estimate = NULL, obs_per_cell = NULL,
                                     partition = "kmeans"))$value
  sf_ <- select_resolution(pf, "cp")
  expect_identical(sf_$structured, TRUE)
  expect_gt(sf_$one_cell, sf_$value * 1.02)
  expect_null(sf_$note)

  # A `cp` column that is not the one the profile computed has no reference.
  edited <- pf
  edited$cp <- rev(edited$cp)
  expect_true(is.na(select_resolution(edited, "cp")$structured))
  # The criteria that are not Cp: an elbow is structure by definition, and
  # reliability says nothing either way.
  expect_true(is.na(select_resolution(pf, "reliability")$structured))
})

test_that("a count with no structure behind it is a classed warning where it becomes cells", {
  pts <- fns_noise()
  p <- fns_catch(resolution_profile(pts, "z", sac = 300, n_levels = 8, nstart = 3,
                                    estimate = NULL, obs_per_cell = NULL,
                                    partition = "kmeans"))$value
  for (x in list(p, select_resolution(p, "cp"))) {
    out <- fns_catch(get_voronoi_seeds(method = "kmeans", n = x, sample_points = pts,
                                       partition = "kmeans"))
    expect_true("spatialkit_no_structure" %in% out$classes)
    w <- out$warnings[[which(vapply(out$warnings, inherits, logical(1),
                                    "spatialkit_no_structure"))[1L]]]
    expect_identical(w$criterion, "cp")
    expect_match(conditionMessage(w), "get_voronoi_seeds\\(\\): `n` is [0-9]+ cells")
    expect_match(conditionMessage(w), "the count is not identified")
    expect_identical(nrow(out$value), as.integer(w$n))
  }
  expect_identical(tryCatch(voronoi_seeds_kmeans(pts, k = p, partition = "kmeans"),
                            spatialkit_no_structure = function(w) "caught"), "caught")
  # A structured field raises none.
  f <- fns_field()
  pf <- fns_catch(resolution_profile(f, "z", sac = fns_sac(f), n_levels = 8, nstart = 3,
                                     estimate = NULL, obs_per_cell = NULL,
                                     partition = "kmeans"))$value
  out <- fns_catch(get_voronoi_seeds(method = "kmeans", n = pf, sample_points = f,
                                     partition = "kmeans"))
  expect_false("spatialkit_no_structure" %in% out$classes)
})

test_that("a variogram refused for showing no spatial structure gives Cp no nugget, and the flag", {
  # The model such a variogram still carries has the nugget a spurious short
  # range left: here 0.05 on a variance of 1.  Cp must not take it.
  pts <- fns_noise()
  sac <- fns_sac(pts, range = 30, nugget = 0.05, psill = 0.95, value = NA_real_,
                 rejected_reason = "no spatial structure detected at the nearest-neighbour scale")
  out <- fns_catch(resolution_profile(pts, "z", sac = sac, n_levels = 8, nstart = 3))
  p <- out$value
  expect_true(any(grepl("its nugget \\(0.05\\) is not a fitted value either", out$text)))
  expect_identical(attr(p, "cp_noise")$source, "finest-level residual mean square")
  expect_gt(attr(p, "cp_noise")$value, 0.8)
  expect_null(attr(p, "variogram"))
  expect_true(all(is.na(p$reliability)))
  s <- select_resolution(p, "cp")
  expect_identical(s$structured, FALSE)
  expect_true(any(grepl("the profile's variogram found no spatial structure", s$note)))
})

test_that("determine_optimal_levels() says a count is not a finding with a class and the numbers", {
  set.seed(5)
  uni <- sf::st_as_sf(data.frame(x = runif(400, 0, 1000), y = runif(400, 0, 1000)),
                      coords = c("x", "y"), crs = 32632)
  out <- fns_catch(determine_optimal_levels(uni))
  expect_true("spatialkit_no_elbow" %in% out$classes)
  w <- out$warnings[[1L]]
  expect_s3_class(w, "spatialkit_no_elbow")
  expect_identical(w$k, out$value[1L])
  expect_identical(w$structured, FALSE)
  expect_lt(w$sag, log(1.25))
  expect_match(conditionMessage(w), "the WSS curve has no elbow")
  # The count itself, and its shape, are what they were.
  expect_type(out$value, "integer")
  expect_null(attributes(out$value))
  expect_identical(tryCatch(determine_optimal_levels(uni),
                            spatialkit_no_elbow = function(w) "caught"), "caught")

  # "na": no count where there is no finding, with the flag and the sag.
  na <- fns_catch(determine_optimal_levels(uni, on_no_elbow = "na"))
  expect_length(na$warnings, 0L)
  expect_true(is.na(na$value))
  expect_type(na$value, "integer")
  expect_identical(attr(na$value, "structured"), FALSE)
  expect_equal(attr(na$value, "sag"), w$sag)
  # "error": a classed error.
  e <- tryCatch(determine_optimal_levels(uni, on_no_elbow = "error"), error = function(e) e)
  expect_s3_class(e, "spatialkit_no_elbow")
  expect_identical(e$k, out$value[1L])
  expect_error(determine_optimal_levels(uni, on_no_elbow = "nope"), "should be one of")
  # NA is not a cell count, and the tessellation functions say where it came from.
  expect_error(get_voronoi_seeds(method = "kmeans", n = na$value, sample_points = uni),
               "on_no_elbow = \"na\"", fixed = TRUE)

  # Eight separated clusters: an elbow, no condition, the same answer under
  # every setting.
  set.seed(6)
  ctr <- cbind(runif(8, 0, 1000), runif(8, 0, 1000))
  g <- rep(1:8, each = 50)
  cl <- sf::st_as_sf(data.frame(x = ctr[g, 1] + rnorm(400, 0, 8), y = ctr[g, 2] + rnorm(400, 0, 8)),
                     coords = c("x", "y"), crs = 32632)
  ok <- fns_catch(determine_optimal_levels(cl, max_levels = 14))
  expect_length(ok$warnings, 0L)
  expect_identical(ok$value[1L], 8L)
  expect_null(attributes(ok$value))
  expect_identical(determine_optimal_levels(cl, max_levels = 14, on_no_elbow = "error"), ok$value)
  # Under "na" the count carries the flag too, so one test reads both outcomes.
  k8 <- determine_optimal_levels(cl, max_levels = 14, on_no_elbow = "na")
  expect_identical(as.integer(k8), ok$value)
  expect_identical(attr(k8, "structured"), TRUE)
  expect_gte(attr(k8, "sag"), log(1.25))
})

test_that("the model-aware diagnostics carry the flag and the sag", {
  set.seed(5)
  n <- 400
  uni <- sf::st_as_sf(data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000),
                                 z = rnorm(n), w = rnorm(n)),
                      coords = c("x", "y"), crs = 32632)
  out <- fns_catch(determine_optimal_levels(uni, max_levels = 40, response_var = "z",
                                            predictor_vars = "w", criterion = "morans_i"))
  d <- attr(out$value, "diagnostics")
  expect_false(is.null(d))
  expect_identical(d$structured, FALSE)
  expect_lt(d$sag, log(1.25))
  expect_true("spatialkit_no_elbow" %in% out$classes)
  # Under "na" the window the ladder placed is not scored at all.
  na <- determine_optimal_levels(uni, max_levels = 40, response_var = "z",
                                 predictor_vars = "w", criterion = "morans_i",
                                 on_no_elbow = "na")
  expect_true(is.na(na))
  expect_identical(attr(na, "structured"), FALSE)
})
