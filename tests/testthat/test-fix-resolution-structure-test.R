# tests/testthat/test-fix-resolution-structure-test.R
# ---------------------------------------------------------------------------
# `structured`, the flag select_resolution() puts on a Cp optimum, needs the
# cells to beat one mean for the layer by more than chance: an F test of the
# cell means at the optimum, adjusted for the levels searched, beside the
# comparison of Cp at one cell with Cp at the optimum.  That comparison alone
# called white noise structured on every layer whose Cp had a noise variance
# of 0, and on one small layer in seven.
# ---------------------------------------------------------------------------

fst_catch <- function(expr) {
  ws <- list()
  value <- withCallingHandlers(expr, warning = function(w) {
    ws[[length(ws) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = ws, classes = unlist(lapply(ws, class)),
       text = vapply(ws, conditionMessage, character(1)))
}
fst_noise <- function(n, seed) {
  set.seed(seed)
  sf::st_as_sf(data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), z = rnorm(n)),
               coords = c("x", "y"), crs = 32632)
}
fst_field <- function(n = 400, seed = 2, a = 100, nugget = 0.3) {
  set.seed(seed)
  xy <- data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000))
  D  <- as.matrix(stats::dist(xy))
  xy$z <- as.numeric(t(chol(exp(-D / a) + diag(1e-8, n))) %*% rnorm(n)) +
    rnorm(n, sd = sqrt(nugget))
  sf::st_as_sf(xy, coords = c("x", "y"), crs = 32632)
}
fst_sac <- function(pts, range = 240, nugget = 0.3, psill = 1) {
  vm <- data.frame(model = c("Nug", "Exp"), psill = c(nugget, psill),
                   range = c(0, range / 3), stringsAsFactors = FALSE)
  structure(range, class = c("sac_range", "numeric"), variogram_model = vm,
            crs = sf::st_crs(pts))
}


test_that("the profile keeps the F test of each level's cells against one mean", {
  pts <- fst_noise(300, 11)
  p <- fst_catch(resolution_profile(pts, "z", sac = 300, n_levels = 6, nstart = 2))$value
  one <- attr(p, "one_cell")
  expect_named(one$f_p, as.character(p$levels))
  expect_true(all(one$f_p >= 0 & one$f_p <= 1))
  # The cells of a level are the Voronoi cells of its centres, so the test
  # can be redone from outside: a one-way analysis of variance of the
  # response on the cells.
  xy <- sf::st_coordinates(pts)
  for (L in p$levels[c(1L, 3L, length(p$levels))]) {
    cen <- attr(p, "centres")[[as.character(L)]]
    cl <- max.col(-(outer(xy[, 1], cen[, 1], "-")^2 + outer(xy[, 2], cen[, 2], "-")^2))
    aov_p <- stats::anova(stats::lm(pts$z ~ factor(cl)))[["Pr(>F)"]][1L]
    expect_equal(one$f_p[[as.character(L)]], aov_p, tolerance = 1e-8)
  }
  # No response, no test.
  expect_null(attr(resolution_profile(pts, n_levels = 3, nstart = 1), "one_cell")$f_p)
})

test_that("a Cp with no noise variance does not call white noise structured", {
  # A nugget of 0: Cp is RSS / n, which falls with every finer level, so one
  # cell is never within `tol` of the optimum.  `structured` was TRUE on 10
  # of 10 such layers, in the call that warns Cp has no penalty, and the
  # tessellation functions raised no spatialkit_no_structure.
  for (seed in 41:45) {
    pts <- fst_noise(400, seed)
    out <- fst_catch(resolution_profile(pts, "z", sac = fst_sac(pts, 90, nugget = 0),
                                        n_levels = 8, nstart = 2))
    expect_true(any(grepl("the variogram's nugget is 0", out$text)))
    expect_identical(attr(out$value, "cp_noise")$value, 0)
    s <- select_resolution(out$value, "cp")
    # The comparison of Cp at one cell with the optimum passes, as it must.
    expect_gt(s$one_cell, s$value * 1.02)
    expect_identical(s$structured, FALSE)
    expect_gte(s$structure_p, 0.05)
    expect_match(s$note, "by no more than chance allows")
    expect_match(s$note, "Cp's noise variance is 0")
    expect_match(s$note, "the count is not identified")
  }
  # And the count reaches the tessellation functions with the classed warning.
  seeds <- fst_catch(get_voronoi_seeds(method = "kmeans", n = s, sample_points = pts))
  expect_true("spatialkit_no_structure" %in% seeds$classes)
  expect_output(print(summary(out$value)), "one cell for the whole layer does as well")
})

test_that("with no noise variance a field with structure is still called structured, and why", {
  f <- fst_field()
  p <- fst_catch(resolution_profile(f, "z", sac = fst_sac(f, 240, nugget = 0),
                                    n_levels = 8, nstart = 2))$value
  s <- select_resolution(p, "cp")
  expect_identical(s$structured, TRUE)
  expect_lt(s$structure_p, 1e-6)
  # The Cp comparison says nothing there, and the result says what the flag
  # rests on.
  expect_match(s$note, "Cp's noise variance is 0")
  expect_match(s$note, "rests on the F test of the cells against one mean alone")
  seeds <- fst_catch(get_voronoi_seeds(method = "kmeans", n = s, sample_points = f))
  expect_false("spatialkit_no_structure" %in% seeds$classes)
})

test_that("on small white-noise layers the flag is held to the level of its test", {
  # Forty layers of 60 points, noise variance from the finest level.  Cp at
  # one cell is above the optimum by more than `tol` on about one in five,
  # which alone used to make `structured` TRUE (27 of 200 such layers).
  res <- t(vapply(1:40, function(i) {
    p <- fst_catch(resolution_profile(fst_noise(60, 300 + i), "z", sac = 300, nstart = 2))$value
    s <- select_resolution(p, "cp")
    c(structured = isTRUE(s$structured), tol_alone = s$one_cell > s$value * 1.02,
      p = s$structure_p, levels = sum(is.finite(p$cp)))
  }, numeric(4)))
  expect_gte(sum(res[, "tol_alone"]), 6)
  expect_lte(sum(res[, "structured"]), 5)
  expect_lt(sum(res[, "structured"]), sum(res[, "tol_alone"]))
  # TRUE only where the adjusted p-value is under 5 percent, and the
  # adjustment is by the levels searched.
  expect_true(all(res[res[, "structured"] == 1, "p"] < 0.05))
  expect_true(all(res[res[, "tol_alone"] == 1 & res[, "structured"] == 0, "p"] >= 0.05))

  # One of the layers the Cp comparison passes and the test does not.
  i <- which(res[, "tol_alone"] == 1 & res[, "structured"] == 0)[1L]
  pts <- fst_noise(60, 300 + i)
  p <- fst_catch(resolution_profile(pts, "z", sac = 300, nstart = 2))$value
  s <- select_resolution(p, "cp")
  expect_equal(s$structure_p,
               min(1, attr(p, "one_cell")$f_p[[as.character(s$best)]] * sum(is.finite(p$cp))))
  expect_match(s$note, sprintf("the cell means at %d cells differ from one mean", s$best))
  expect_match(s$note, sprintf("adjusting for the %d levels searched", sum(is.finite(p$cp))))
  expect_match(s$note, "a gap that chance produces on a layer of this size")
  expect_output(print(s), "not identified")
  w <- fst_catch(get_voronoi_seeds(method = "kmeans", n = s, sample_points = pts))
  expect_true("spatialkit_no_structure" %in% w$classes)
})

test_that("a field with structure passes both, on a small layer too", {
  f <- fst_field(300, seed = 21)
  p <- fst_catch(resolution_profile(f, "z", sac = 300, n_levels = 8, nstart = 2))$value
  s <- select_resolution(p, "cp")
  expect_identical(s$structured, TRUE)
  expect_lt(s$structure_p, 0.05)
  expect_gt(s$one_cell, s$value * 1.02)
  expect_null(s$note)
  # The criteria that are not Cp carry no p-value.
  expect_true(is.na(select_resolution(p, "moran_z", diagnostic = TRUE)$structure_p))
})

test_that("without the profile's F test the flag is not asserted", {
  f <- fst_field()
  p <- fst_catch(resolution_profile(f, "z", sac = fst_sac(f), n_levels = 8, nstart = 2))$value
  expect_identical(select_resolution(p, "cp")$structured, TRUE)
  bare <- p
  one <- attr(bare, "one_cell")
  one$f_p <- NULL
  attr(bare, "one_cell") <- one
  s <- select_resolution(bare, "cp")
  expect_true(is.na(s$structured))
  expect_true(is.na(s$structure_p))
  expect_match(s$note, "carries no F test")
  # Where Cp at one cell is within `tol`, that alone still says FALSE.
  n <- fst_noise(600, 7)
  pn <- fst_catch(resolution_profile(n, "z", sac = 300, n_levels = 8, nstart = 3))$value
  one <- attr(pn, "one_cell")
  one$f_p <- NULL
  attr(pn, "one_cell") <- one
  expect_identical(select_resolution(pn, "cp")$structured, FALSE)
})
