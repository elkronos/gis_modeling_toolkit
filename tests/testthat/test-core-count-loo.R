# tests/testthat/test-core-count-loo.R
# ---------------------------------------------------------------------------
# resolution_profile(estimate = ): the leave-one-out error of the estimate the
# summary reports, on every scored row of the layer, in the cells the seeds
# will build.
#
# What a profile has without an estimate is held identical() to the same call
# with one, conditions included.  The five columns are held to the reference
# functions of helper-pool-reference.R on the regions .nearest_centre() gives
# from the profile's centres.
#
# Every call names `partition`, `nstart`, `estimate`, `obs_per_cell` and
# `sac`.  Unless a test says otherwise the regions are of the density kind
# from 3 k-means++ starts.  A layer is drawn under the seed its test states:
# the locations (x, then y), then the response.  A backend is switched off by
# replacing .has_backend().
# ---------------------------------------------------------------------------

# ---- Layers and helpers ------------------------------------------------------

.co_seed <- function(seed, env = parent.frame())
  withr::local_seed(seed, .local_envir = env, .rng_kind = "Mersenne-Twister",
                    .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")

.co_points <- function(x, y, z, crs = 32617)
  sf::st_as_sf(data.frame(x = x, y = y, z = z), coords = c("x", "y"), crs = crs)

# `n` locations uniform on a square of side `side`, with the field of the
# engine's tests: sin(x / 2e4) + cos(y / 1.5e4) plus normal noise of standard
# deviation 1.5.
.co_layer <- function(n, seed, side = 1e5) {
  d <- local({
    .co_seed(seed)
    x <- stats::runif(n, 0, side); y <- stats::runif(n, 0, side)
    data.frame(x = x, y = y, z = sin(x / 2e4) + cos(y / 1.5e4) + stats::rnorm(n, 0, 1.5))
  })
  .co_points(d$x, d$y, d$z)
}

# The lattice of the engine's property tests: 60 x 50 locations at a spacing
# of 0.7 from (0.1, 0.3).  None of its coordinates is exact in binary and
# every location has neighbours at exactly equal distances.
.co_lattice <- function() {
  xy <- as.matrix(expand.grid(x = seq(0.1, by = 0.7, length.out = 60),
                              y = seq(0.3, by = 0.7, length.out = 50)))
  z <- local({
    .co_seed(31)
    sin(xy[, 1] / 9) + cos(xy[, 2] / 7) + stats::rnorm(nrow(xy), 0, 0.5)
  })
  .co_points(xy[, 1], xy[, 2], z)
}

.co_profile <- function(pts, ..., response = "z", partition = "density", nstart = 3,
                        sac = "auto", estimate = "pooled", obs_per_cell = 4.5)
  resolution_profile(pts, response, ..., partition = partition, nstart = nstart, sac = sac,
                     estimate = estimate, obs_per_cell = obs_per_cell)

# A call with its warnings kept and muffled.
.co_catch <- function(expr) {
  got <- list()
  value <- withCallingHandlers(expr, warning = function(w) {
    got[[length(got) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = got)
}
.co_classes <- function(r) vapply(r$warnings, function(w) class(w)[1L], character(1))
.co_texts   <- function(r) vapply(r$warnings, conditionMessage, character(1))
.co_of <- function(r, class) Filter(function(w) inherits(w, class), r$warnings)
.co_fields <- function(w) unclass(w)[setdiff(names(unclass(w)), c("message", "call"))]
.co_error <- function(expr) tryCatch({ suppressWarnings(expr); NULL }, error = function(e) e)

.co_without_fnn <- function(env = parent.frame())
  testthat::local_mocked_bindings(
    .has_backend = function(p)
      if (identical(p, "FNN")) FALSE else requireNamespace(p, quietly = TRUE),
    .env = env)

.co_new <- c("scored_per_cell", "loo", "loo_se", "share_alone", "own_weight")
.co_xy  <- function(pts) unname(sf::st_coordinates(pts)[, 1:2, drop = FALSE])

# The regions of a level as the profile scores them: each scored row in the
# cell of its nearest centre, the centres that hold none dropped and the rest
# numbered in the order of their centres.
.co_regions <- function(xy, centres) {
  g <- .nearest_centre(xy[, 1], xy[, 2], centres[, 1], centres[, 2])
  u <- sort(unique(g))
  list(gi = match(g, u), K = length(u), u = u)
}

# What the reference functions give on those regions, with the rows in the
# order the layer has them.
.co_reference <- function(xy, y, centres) {
  rg  <- .co_regions(xy, centres)
  nn1 <- pool_nearest_other(xy)
  fr  <- fit_regions(rg$gi, rg$gi, y, xy, xy, nn1)
  N   <- length(y)
  list(K = fr$K, N = N,
       loo_mean = loo_plain(y, rg$gi, oneway(y, rg$gi, rg$K), nn1),
       loo_mean_se = stats::sd(fr$eP^2) / sqrt(N),
       loo_pooled = fr$loo_L, loo_pooled_se = fr$loo_L_se,
       own_weight = stats::median(fr$w), share_alone = sum(fr$n == 1L) / N,
       one_row = mean(fr$n == 1L), sigma2 = fr$sigma2, tau2 = fr$tau2)
}

# Hold a profile's five columns and its record to the reference, level by
# level, to a relative 1e-10: the profile's rows are in coordinate order, so
# its sums run in another order than the reference's.
.co_expect_reference <- function(p, xy, y, estimate, what = "") {
  est <- attr(p, "estimate")
  cen <- attr(p, "centres")
  for (i in seq_len(nrow(p))) {
    at  <- paste(what, estimate, p$levels[i], "cells")
    ref <- .co_reference(xy, y, cen[[as.character(p$levels[i])]])
    expect_equal(p$loo[i], if (estimate == "mean") ref$loo_mean else ref$loo_pooled,
                 tolerance = 1e-10, info = at)
    expect_equal(p$loo_se[i], if (estimate == "mean") ref$loo_mean_se else ref$loo_pooled_se,
                 tolerance = 1e-10, info = at)
    if (estimate == "pooled")
      expect_equal(p$own_weight[i], ref$own_weight, tolerance = 1e-10, info = at)
    else expect_identical(p$own_weight[i], NA_real_, info = at)
    expect_identical(p$share_alone[i], ref$share_alone, info = at)
    expect_identical(p$scored_per_cell[i], ref$N / ref$K, info = at)
    # The share of regions with one row, which the experiment's summary keeps.
    expect_equal(p$share_alone[i] * ref$N / ref$K, ref$one_row, tolerance = 1e-12, info = at)
    expect_identical(est$cells[[i]], ref$K, info = at)
    expect_equal(est$sigma2[[i]], ref$sigma2, tolerance = 1e-10, info = at)
    if (estimate == "pooled") expect_equal(est$tau2[[i]], ref$tau2, tolerance = 1e-8, info = at)
    else expect_identical(est$tau2[[i]], NA_real_, info = at)
  }
  invisible(NULL)
}

.co_store <- new.env(parent = emptyenv())
.co_keep <- function(key, expr) {
  if (is.null(.co_store[[key]])) .co_store[[key]] <- expr
  .co_store[[key]]
}
.co_600  <- function() .co_keep("l600", .co_layer(600, 401))
.co_1200 <- function() .co_keep("l1200", .co_layer(1200, 402))


# ---- An estimate changes nothing a profile has without one --------------------

# The same call with and without an estimate: every column and attribute of
# the call without is in the call with, identical; every condition of the
# call without is raised again, same class and text, same order; and what the
# estimate adds is five columns, one attribute and conditions of `extra`.
.co_expect_unchanged <- function(call, what, extra = character(0)) {
  base <- .co_catch(call(NULL))
  cols <- names(base$value)
  at   <- attributes(base$value)
  for (est in c("pooled", "mean")) {
    info <- paste(what, est)
    r <- .co_catch(call(est))
    expect_identical(names(r$value), c(cols, .co_new), info = info)
    for (col in cols) expect_identical(r$value[[col]], base$value[[col]], info = paste(info, col))
    expect_identical(setdiff(names(attributes(r$value)), names(at)), "estimate", info = info)
    for (a in setdiff(names(at), "names"))
      expect_identical(attr(r$value, a, exact = TRUE), at[[a]], info = paste(info, a))
    added <- .co_classes(r) %in% extra
    expect_identical(.co_classes(r)[!added], .co_classes(base), info = info)
    expect_identical(.co_texts(r)[!added], .co_texts(base), info = info)
    expect_identical(sort(unique(.co_classes(r)[added])), sort(extra), info = info)
    expect_identical(attr(r$value, "estimate")$value, est, info = info)
  }
  invisible(base)
}

test_that("an estimate adds five columns and one attribute, and changes nothing else", {
  skip_if_not_installed("gstat")
  below <- .co_layer(400, 403)
  above <- .co_1200()
  for (kind in c("kmeans", "density")) {
    # 5 starts for the kind "kmeans", to keep the test short; the number is
    # named either way.
    ns <- if (kind == "kmeans") 5 else 3
    # The ladder the function builds, on a layer below `sample_n` and on one above it.
    .co_expect_unchanged(function(est)
      .co_profile(below, partition = kind, nstart = ns, estimate = est, obs_per_cell = NULL,
                  sac = NULL, n_levels = 6, sample_n = 1500),
      paste(kind, "automatic, below sample_n"))
    .co_expect_unchanged(function(est)
      .co_profile(above, partition = kind, nstart = ns, estimate = est, obs_per_cell = NULL,
                  sac = NULL, n_levels = 6, sample_n = 400),
      paste(kind, "automatic, above sample_n"))
    # A ladder given as `levels`.
    .co_expect_unchanged(function(est)
      .co_profile(below, partition = kind, nstart = ns, estimate = est, obs_per_cell = NULL,
                  sac = NULL, levels = c(8L, 30L, 120L), sample_n = 1500),
      paste(kind, "levels, below sample_n"))
    .co_expect_unchanged(function(est)
      .co_profile(above, partition = kind, nstart = ns, estimate = est, obs_per_cell = NULL,
                  sac = NULL, levels = c(8L, 30L, 120L), sample_n = 400),
      paste(kind, "levels, above sample_n"))
    # The response read on one half.
    .co_expect_unchanged(function(est)
      .co_profile(above, partition = kind, nstart = ns, estimate = est, obs_per_cell = NULL,
                  sac = NULL, n_levels = 5, sample_n = 400, select_on = "split"),
      paste(kind, "split"))
  }
})

test_that("an estimate changes nothing on a constant response, a supplied sac or missing values", {
  skip_if_not_installed("gstat")
  pts <- .co_layer(400, 403)
  flat <- pts
  flat$z <- 0.1
  base <- .co_expect_unchanged(function(est)
    .co_profile(flat, estimate = est, obs_per_cell = NULL, sac = NULL, n_levels = 5),
    "constant")
  expect_true("spatialkit_constant_response" %in% .co_classes(base))

  .co_expect_unchanged(function(est)
    .co_profile(pts, estimate = est, obs_per_cell = NULL, sac = 20000, n_levels = 5,
                range_floor = TRUE),
    "a range as sac")

  holes <- pts
  holes$z[c(3, 40, 41, 200:260)] <- NA
  base <- .co_expect_unchanged(function(est)
    .co_profile(holes, estimate = est, obs_per_cell = NULL, sac = NULL, n_levels = 5,
                min_cell_n = 5),
    "missing values")
  expect_true("spatialkit_response_ceiling" %in% .co_classes(base))

  # Rows that share locations: one further condition, of a new class.
  twice <- pts
  sf::st_geometry(twice) <- sf::st_geometry(pts)[rep(1:200, 2)]
  .co_expect_unchanged(function(est)
    .co_profile(twice, estimate = est, obs_per_cell = NULL, sac = NULL, n_levels = 5),
    "repeated locations", extra = "spatialkit_repeated_locations")
})

test_that("the estimate leaves the caller's random-number stream as the call without it does", {
  pts <- .co_600()
  after <- function(est) {
    .co_seed(77)
    suppressWarnings(.co_profile(pts, estimate = est, obs_per_cell = 9))
    stats::runif(3)
  }
  expect_identical(after("pooled"), after(NULL))
  expect_identical(after("mean"), after(NULL))
  expect_identical(after(NULL), local({ .co_seed(77); stats::runif(3) }))
})


# ---- The five columns, against the reference functions ------------------------

test_that("loo, its standard error and the weights are the reference functions' at each level", {
  skip_if_not_installed("FNN")
  pts <- .co_600()
  xy  <- .co_xy(pts)
  for (est in c("mean", "pooled")) {
    for (kind in c("density", "kmeans")) {
      p <- .co_profile(pts, partition = kind, estimate = est, obs_per_cell = c(18, 9, 4.5, 2))
      expect_identical(p$levels, c(33L, 67L, 133L, 300L))
      expect_true(all(is.finite(p$loo)))
      .co_expect_reference(p, xy, pts$z, est, kind)
      # On a ladder of sizes without a split the layer's count is the fit's.
      expect_equal(p$scored_per_cell, p$cell_resp_mean, tolerance = 1e-12)
      # There are rows alone in their cell at the finest size, so the rule for them is reached.
      expect_gt(p$share_alone[4L], 0.05)
      # The record.
      e <- attr(p, "estimate")
      expect_identical(names(e),
                       c("value", "neighbours", "response", "n", "rows", "features", "seed",
                         "constant", "loo_one_cell", "loo_one_cell_se", "loo_levels", "sigma2",
                         "tau2", "cells", "f_p"))
      expect_identical(e$value, est)
      expect_identical(e$neighbours, if (est == "pooled") 8L else NA_integer_)
      expect_identical(e$response, "z")
      expect_identical(e$n, 600L)
      expect_identical(e$rows, "all")
      expect_identical(e$features, "points")
      expect_identical(e$seed, 123L)
      expect_false(e$constant)
      expect_identical(e$loo_levels, stats::setNames(p$loo, p$levels))
      # One mean for all rows: each row predicted by the mean of the others.
      e1 <- (pts$z - mean(pts$z)) * 600 / 599
      expect_equal(e$loo_one_cell, mean(e1^2), tolerance = 1e-12)
      expect_equal(e$loo_one_cell_se, stats::sd(e1^2) / sqrt(600), tolerance = 1e-12)
      expect_identical(names(e$f_p), as.character(p$levels))
      expect_true(all(is.finite(e$f_p) & e$f_p >= 0 & e$f_p <= 1))
    }
  }
  # The same centres, so the same cells, whichever estimate is scored.
  expect_identical(attr(.co_profile(pts, estimate = "mean", obs_per_cell = 9), "centres"),
                   attr(.co_profile(pts, estimate = "pooled", obs_per_cell = 9), "centres"))
})

test_that("on the automatic ladder loo is the reference's on the ladder's cells", {
  skip_if_not_installed("FNN")
  skip_if_not_installed("gstat")
  pts <- .co_600()
  xy  <- .co_xy(pts)
  for (est in c("mean", "pooled")) {
    p <- suppressWarnings(.co_profile(pts, partition = "kmeans", nstart = 5, estimate = est,
                                      obs_per_cell = NULL, sac = NULL, n_levels = 6))
    expect_identical(p$levels[1L], 2L)
    .co_expect_reference(p, xy, pts$z, est, "automatic")
  }
})

test_that("on a layer above sample_n loo reads every scored row of the layer", {
  skip_if_not_installed("FNN")
  skip_if_not_installed("gstat")
  pts <- .co_1200()
  xy  <- .co_xy(pts)
  for (est in c("mean", "pooled")) {
    p <- suppressWarnings(.co_profile(pts, estimate = est, obs_per_cell = NULL, sac = NULL,
                                      n_levels = 5, sample_n = 300))
    b <- attr(p, "bounds")
    expect_identical(b$n_sample, 300L)
    expect_identical(attr(p, "estimate")$n, 1200L)
    expect_identical(attr(p, "estimate")$rows, "all")
    # The reference on all 1,200 rows, in the cells of the centres fitted to 300.
    .co_expect_reference(p, xy, pts$z, est, "above sample_n")
    # `cell_resp_mean` counts the subsample's rows and `scored_per_cell` the layer's.
    expect_equal(p$scored_per_cell[1L], 600, tolerance = 0.2)
    expect_equal(p$cell_resp_mean[1L], 150, tolerance = 0.2)
    expect_true(all(p$scored_per_cell > 2.5 * p$cell_resp_mean))
    # The record is the layer's as well.  One mean for all 1,200 rows, where
    # `one_cell` holds the subsample's 300.
    e   <- attr(p, "estimate")
    one <- attr(p, "one_cell")
    expect_identical(one$n, 300L)
    e1 <- (pts$z - mean(pts$z)) * 1200 / 1199
    expect_equal(e$loo_one_cell, mean(e1^2), tolerance = 1e-12)
    expect_equal(e$loo_one_cell_se, stats::sd(e1^2) / sqrt(1200), tolerance = 1e-12)
    expect_false(isTRUE(all.equal(e$loo_one_cell, one$rss / 300 * (300 / 299)^2, tolerance = 1e-6)))
    # The F test of each level's cells against one mean: the cells of the
    # nearest centres on every row, where `one_cell$f_p` reads the fit's
    # clusters on the subsample.
    f_p <- vapply(seq_len(nrow(p)), function(i) {
      rg  <- .co_regions(xy, attr(p, "centres")[[i]])
      rss <- sum((pts$z - stats::ave(pts$z, rg$gi))^2)
      tss <- sum((pts$z - mean(pts$z))^2)
      stats::pf(((tss - rss) / (rg$K - 1)) / (rss / (1200 - rg$K)), rg$K - 1, 1200 - rg$K,
                lower.tail = FALSE)
    }, numeric(1))
    expect_equal(unname(e$f_p), f_p, tolerance = 1e-8)
    expect_identical(names(e$f_p), as.character(p$levels))
    expect_false(isTRUE(all.equal(unname(e$f_p), unname(one$f_p), tolerance = 1e-6)))
    expect_false(e$constant)
  }
})

test_that("the estimate's reference is constant when the layer is, not when the subsample is", {
  skip_if_not_installed("gstat")
  # 2,000 rows, a response of 0 on all but six: a subsample of 300 that holds
  # none of the six is constant, and the layer is not.
  pts <- .co_layer(2000, 408)
  pts$z <- 0
  pts$z[local({ .co_seed(409); sample.int(2000L, 6L) })] <- c(1, 2, 3, 1, 2, 3)
  call <- function(seed, est = "pooled")
    .co_catch(.co_profile(pts, partition = "kmeans", nstart = 2, estimate = est,
                          obs_per_cell = NULL, sac = NULL, n_levels = 4, sample_n = 300,
                          seed = seed))
  found <- NA_integer_
  for (seed in 1:30) {
    r <- call(seed)
    if ("spatialkit_constant_response" %in% .co_classes(r)) { found <- seed; break }
  }
  expect_false(is.na(found))
  expect_false(identical(found, 123L))
  p <- r$value
  # The present columns follow the subsample: `cp` is NA at every level.
  expect_true(attr(p, "one_cell")$constant)
  expect_true(all(is.na(p$cp)))
  # The estimate's follow the layer.
  e <- attr(p, "estimate")
  expect_false(e$constant)
  expect_identical(e$n, 2000L)
  expect_true(all(is.finite(p$loo)))
  expect_true(is.finite(e$loo_one_cell))
  e1 <- (pts$z - mean(pts$z)) * 2000 / 1999
  expect_equal(e$loo_one_cell, mean(e1^2), tolerance = 1e-12)
  expect_false(any(grepl("the response is constant", grep("^  estimate", utils::capture.output(print(p)),
                                                          value = TRUE), fixed = TRUE)))
  # The record holds the seed the profile was built under, whatever it is.
  expect_identical(e$seed, found)
  expect_identical(attr(call(77L)$value, "estimate")$seed, 77L)
  none <- local({ .co_seed(410); attr(call(NULL)$value, "estimate") })
  expect_true("seed" %in% names(none))
  expect_null(none$seed)
})

test_that("a level where no cell holds two scored rows has no loo", {
  pts <- .co_layer(40, 404)
  xy  <- .co_xy(pts)
  # Every row its own centre.
  for (est in c("mean", "pooled")) {
    lv <- .estimate_at_level(xy, pts$z, xy, est, use_fnn = FALSE, cache = new.env())
    expect_identical(lv$cells, 40L)
    expect_identical(lv$scored_per_cell, 1)
    expect_identical(lv$share_alone, 1)
    expect_identical(lv$loo, NA_real_)
    expect_identical(lv$loo_se, NA_real_)
    expect_identical(lv$own_weight, NA_real_)
  }
  # One centre: both estimates give the error of one mean.
  one <- matrix(c(5e4, 5e4, 9e9, 9e9), 2, 2, byrow = TRUE)
  e1  <- (pts$z - mean(pts$z)) * 40 / 39
  for (est in c("mean", "pooled")) {
    lv <- .estimate_at_level(xy, pts$z, one, est, use_fnn = FALSE, cache = new.env())
    expect_identical(lv$cells, 1L)
    expect_equal(lv$loo, mean(e1^2), tolerance = 1e-12)
    expect_identical(lv$share_alone, 0)
    expect_identical(lv$own_weight, if (est == "pooled") 1 else NA_real_)
  }
})


# ---- The rows loo reads --------------------------------------------------------

test_that("under select_on = \"split\" only the selection half is scored", {
  skip_if_not_installed("FNN")
  pts <- .co_1200()
  xy  <- .co_xy(pts)
  for (est in c("mean", "pooled")) {
    p <- suppressWarnings(.co_profile(pts, estimate = est, obs_per_cell = 4.5,
                                      select_on = "split"))
    sel <- attr(p, "split")$selection
    e   <- attr(p, "estimate")
    expect_identical(e$rows, "selection half")
    expect_identical(e$n, length(sel))
    expect_lt(e$n, 800L)
    expect_gt(e$n, 400L)
    # The size counts the rows of both halves: the cells are built on every point.
    expect_identical(p$levels, 267L)
    expect_identical(attr(p, "bounds")$n_scored, length(sel))
    expect_identical(attr(p, "bounds")$n_response, 1200L)
    # The reference on the selection rows alone, in the cells of every point.
    .co_expect_reference(p, xy[sel, , drop = FALSE], pts$z[sel], est, "split")
    # The halves are spatial blocks, so the selection rows fill about half
    # of the cells: fewer cells hold a scored row, at about the size each.
    expect_lt(e$cells[[1L]], 0.7 * 267)
    expect_gt(e$cells[[1L]], 0.3 * 267)
    expect_identical(p$scored_per_cell, length(sel) / e$cells[[1L]])
    expect_gt(p$scored_per_cell, 2.25)
    expect_lt(p$scored_per_cell, 6)
    out <- utils::capture.output(print(p))
    expect_true(any(grepl(sprintf("loo on %d rows of the selection half; one cell: loo ",
                                  length(sel)), out, fixed = TRUE)))
  }
})

test_that("rows without a response stay in the geometry and out of loo", {
  skip_if_not_installed("FNN")
  pts <- .co_600()
  pts$z[seq(1, 600, by = 4)] <- NA
  ok  <- is.finite(pts$z)
  xy  <- .co_xy(pts)
  for (est in c("mean", "pooled")) {
    p <- .co_profile(pts, estimate = est, obs_per_cell = c(9, 3))
    expect_identical(attr(p, "estimate")$n, 450L)
    expect_identical(p$levels, c(50L, 150L))
    expect_identical(p$cell_n_median > 0, c(TRUE, TRUE))
    .co_expect_reference(p, xy[ok, , drop = FALSE], pts$z[ok], est, "missing")
  }
})

test_that("the same layer with its rows permuted gives the same profile", {
  skip_if_not_installed("gstat")
  pts  <- .co_1200()
  perm <- local({ .co_seed(405); sample.int(1200L) })
  for (est in c("mean", "pooled")) {
    expect_identical(.co_profile(pts[perm, ], estimate = est, obs_per_cell = c(9, 4.5, 2)),
                     .co_profile(pts, estimate = est, obs_per_cell = c(9, 4.5, 2)), info = est)
    # A subsample is drawn from the rows in coordinate order.
    expect_identical(
      suppressWarnings(.co_profile(pts[perm, ], partition = "kmeans", nstart = 5, estimate = est,
                                   obs_per_cell = NULL, sac = NULL, n_levels = 5, sample_n = 400)),
      suppressWarnings(.co_profile(pts, partition = "kmeans", nstart = 5, estimate = est,
                                   obs_per_cell = NULL, sac = NULL, n_levels = 5, sample_n = 400)),
      info = est)
  }
})

test_that("on a lattice the profile does not depend on the order of the rows", {
  lat  <- .co_lattice()
  perm <- local({ .co_seed(406); sample.int(nrow(lat)) })
  for (est in c("mean", "pooled")) {
    a <- .co_profile(lat, estimate = est, obs_per_cell = c(4.5, 2))
    expect_identical(a$levels, c(667L, 1500L))
    # Rows alone in their cell, whose nearest other row ties on a lattice.
    expect_gt(a$share_alone[2L], 0.05)
    expect_true(all(is.finite(a$loo)))
    expect_identical(.co_profile(lat[perm, ], estimate = est, obs_per_cell = c(4.5, 2)), a,
                     info = est)
  }
})

test_that("a lon/lat layer is scored as the same layer projected first", {
  ll <- sf::st_transform(.co_600(), 4326)
  a  <- suppressWarnings(.co_profile(ll, estimate = "pooled", obs_per_cell = 9))
  b  <- suppressWarnings(.co_profile(ensure_projected(ll), estimate = "pooled", obs_per_cell = 9))
  expect_identical(a, b)
  expect_true(is.finite(a$loo))
})

test_that("features that are not points are scored at the point each is reduced to", {
  pts  <- .co_layer(300, 407)
  sq   <- sf::st_buffer(pts, 50, endCapStyle = "SQUARE")
  expect_false(any(sf::st_geometry_type(sq) == "POINT"))
  p    <- suppressWarnings(.co_profile(sq, estimate = "pooled", obs_per_cell = 9))
  expect_identical(attr(p, "estimate")$features, "other")
  expect_identical(attr(p, "estimate")$n, 300L)
  as_points <- suppressWarnings(.co_profile(coerce_to_points(sq, "auto"), estimate = "pooled",
                                            obs_per_cell = 9))
  expect_identical(attr(as_points, "estimate")$features, "points")
  expect_identical(p$loo, as_points$loo)
  expect_identical(attr(p, "centres"), attr(as_points, "centres"))
})


# ---- A constant response, and repeated locations ------------------------------

test_that("a constant response has no loo, and the present warning says so in its present words", {
  pts <- .co_600()
  pts$z <- 0.1
  for (est in c("mean", "pooled")) {
    r <- .co_catch(.co_profile(pts, estimate = est, obs_per_cell = c(18, 4.5)))
    expect_identical(.co_classes(r), "spatialkit_constant_response", info = est)
    expect_identical(.co_texts(r),
                     sprintf(paste0("resolution_profile(): the response is constant over the 600 ",
                                    "scored row(s) (its total sum of squares, %.3g, is rounding ",
                                    "error at most), so no cell count fits it better than another: ",
                                    "`cp`, `moran_i` and `moran_z` are NA at every level."),
                             attr(r$value, "one_cell")$rss), info = est)
    p <- r$value
    expect_identical(p$loo, c(NA_real_, NA_real_))
    expect_identical(p$loo_se, c(NA_real_, NA_real_))
    expect_identical(p$own_weight, c(NA_real_, NA_real_))
    e <- attr(p, "estimate")
    expect_true(e$constant)
    expect_identical(e$loo_one_cell, NA_real_)
    expect_identical(e$loo_one_cell_se, NA_real_)
    # What does not read the values is computed.
    expect_true(all(is.finite(p$scored_per_cell)))
    expect_true(all(is.finite(p$share_alone)))
    out <- utils::capture.output(print(p))
    expect_true(any(grepl("one cell: loo NA (the response is constant)", out, fixed = TRUE)))
    # The warning is the one the call raises without an estimate.
    none <- .co_catch(.co_profile(pts, estimate = NULL, obs_per_cell = c(18, 4.5)))
    expect_identical(.co_texts(none), .co_texts(r), info = est)
  }
  # A response that varies by rounding error alone is constant by the same rule.
  pts$z <- 1e6 + 1e-9 * sin(seq_len(600))
  p <- suppressWarnings(.co_profile(pts, estimate = "pooled", obs_per_cell = 4.5))
  expect_true(attr(p, "estimate")$constant)
  expect_identical(p$loo, NA_real_)
})

test_that("rows that share locations raise spatialkit_repeated_locations once", {
  # 200 stations, three visits each.
  pts <- .co_600()
  sf::st_geometry(pts) <- sf::st_geometry(pts)[rep(1:200, 3)]
  for (est in c("mean", "pooled")) {
    r <- .co_catch(.co_profile(pts, estimate = est, obs_per_cell = c(18, 9)))
    expect_identical(.co_classes(r), "spatialkit_repeated_locations", info = est)
    w <- r$warnings[[1L]]
    expect_identical(conditionMessage(w),
                     paste0("resolution_profile(): 600 of the 600 rows with a value of 'z' share a ",
                            "location with another row. The leave-one-out error leaves one row out ",
                            "and keeps the others at its location, so it, and an interval built ",
                            "from it, describes a further observation at a place already observed, ",
                            "not an observation at a new place."))
    expect_identical(.co_fields(w), list(n = 600L, total = 600L, column = "z"))
    expect_null(conditionCall(w))
    expect_true(all(is.finite(r$value$loo)), info = est)
  }
  # It is said after the levels are scored, of a `loo` the profile holds.
  events <- character(0)
  real <- .partition_centres
  force(real)
  local({
    testthat::local_mocked_bindings(.partition_centres = function(...) {
      events <<- c(events, "fit")
      real(...)
    })
    withCallingHandlers(.co_profile(pts, estimate = "pooled", obs_per_cell = c(18, 9)),
                        warning = function(w) {
                          events <<- c(events, class(w)[1L])
                          invokeRestart("muffleWarning")
                        })
  })
  expect_identical(events, c("fit", "fit", "spatialkit_repeated_locations"))
  # Some rows share and some do not: the rows that share are counted.
  part <- .co_600()
  sf::st_geometry(part)[1:10] <- sf::st_geometry(part)[11:20]
  r <- .co_catch(.co_profile(part, estimate = "pooled", obs_per_cell = 9))
  expect_identical(.co_of(r, "spatialkit_repeated_locations")[[1L]]$n, 20L)
  expect_identical(.co_of(r, "spatialkit_repeated_locations")[[1L]]$total, 600L)
  # `total` is the rows with a value, not the points of the layer.
  holes <- part
  holes$z[501:600] <- NA
  r <- .co_catch(.co_profile(holes, estimate = "pooled", obs_per_cell = 9))
  w <- .co_of(r, "spatialkit_repeated_locations")[[1L]]
  expect_identical(.co_fields(w), list(n = 20L, total = 500L, column = "z"))
  expect_match(conditionMessage(w), "20 of the 500 rows with a value of 'z' share a location",
               fixed = TRUE)
  # A row that shares its location with a row without a response shares it
  # with no scored row.
  part$z[11:20] <- NA
  r <- .co_catch(.co_profile(part, estimate = "pooled", obs_per_cell = 9))
  expect_length(.co_of(r, "spatialkit_repeated_locations"), 0L)
  # Without an estimate nothing is left out, and nothing is said.
  r <- .co_catch(.co_profile(pts, estimate = NULL, obs_per_cell = c(18, 9)))
  expect_length(r$warnings, 0L)
  # Nor on a constant response, whose loo is NA.
  pts$z <- 3.7
  r <- .co_catch(.co_profile(pts, estimate = "pooled", obs_per_cell = 18))
  expect_identical(.co_classes(r), "spatialkit_constant_response")
})

test_that("with repeated locations each row is left out alone", {
  skip_if_not_installed("FNN")
  # Two visits to each of 300 stations: a row's nearest other row is the
  # other visit, which pool_nearest_other() does not serve, so the pooled
  # value, which needs no nearest row, is the one held to the reference.
  pts <- .co_600()
  sf::st_geometry(pts) <- sf::st_geometry(pts)[rep(1:300, 2)]
  xy  <- .co_xy(pts)
  p   <- suppressWarnings(.co_profile(pts, estimate = "pooled", obs_per_cell = c(9, 3)))
  for (i in 1:2) {
    rg <- .co_regions(xy, attr(p, "centres")[[i]])
    fr <- fit_regions(rg$gi, rg$gi, pts$z, xy, xy, nn1 = seq_len(600))
    expect_equal(p$loo[i], fr$loo_L, tolerance = 1e-10)
    expect_equal(p$own_weight[i], stats::median(fr$w), tolerance = 1e-10)
    # Both visits of a station are in one cell, so no row is alone in its cell.
    expect_identical(p$share_alone[i], 0)
  }
})


# ---- Without FNN ---------------------------------------------------------------

test_that("without FNN every column of the profile is the run's with it", {
  skip_if_not_installed("FNN")
  pts <- .co_1200()
  for (kind in c("kmeans", "density")) for (est in c("mean", "pooled")) {
    what <- paste(kind, est)
    with <- .co_profile(pts, partition = kind, estimate = est, obs_per_cell = c(9, 4.5, 2))
    without <- local({
      .co_without_fnn()
      # Nothing reaches the kd-tree: every step is told the answer.
      testthat::local_mocked_bindings(
        get.knn  = function(...) stop("FNN::get.knn() was reached"),
        get.knnx = function(...) stop("FNN::get.knnx() was reached"), .package = "FNN")
      .co_profile(pts, partition = kind, estimate = est, obs_per_cell = c(9, 4.5, 2))
    })
    expect_lt(max(with$levels), 5000L)
    for (col in setdiff(names(with), c("moran_i", "moran_z")))
      expect_identical(without[[col]], with[[col]], info = paste(what, col))
    # Their weights are then held as a dense matrix.
    expect_equal(without$moran_i, with$moran_i, tolerance = 1e-10, info = what)
    expect_equal(without$moran_z, with$moran_z, tolerance = 1e-10, info = what)
    expect_identical(attr(without, "estimate"), attr(with, "estimate"), info = what)
    expect_identical(attr(without, "centres"), attr(with, "centres"), info = what)
    expect_identical(attr(without, "bounds"), attr(with, "bounds"), info = what)
    expect_gt(with$share_alone[3L], 0)
  }
})

test_that("a profile asks for each backend once and hands the answers down", {
  pts <- .co_600()
  asked <- character(0)
  counting <- function(expr) {
    asked <<- character(0)
    testthat::local_mocked_bindings(.has_backend = function(pkg) {
      asked <<- c(asked, pkg)
      requireNamespace(pkg, quietly = TRUE)
    })
    suppressWarnings(expr)
    sort(asked)
  }
  for (est in list("mean", "pooled", NULL)) {
    # A ladder of sizes, three levels.
    expect_identical(counting(.co_profile(pts, estimate = est, obs_per_cell = c(9, 4.5, 2))),
                     c("FNN", "Matrix"))
    # A ladder given as `levels`, of the kind "kmeans".
    expect_identical(counting(.co_profile(pts, partition = "kmeans", nstart = 2, estimate = est,
                                          obs_per_cell = NULL, sac = 20000,
                                          levels = c(12L, 40L, 90L))),
                     c("FNN", "Matrix"))
  }
  # Each is asked for when a step first needs the answer, so a profile none
  # of whose steps reads one asks for neither: a geometry-only profile of the
  # kind "kmeans" loads no optional package, as it never has.
  expect_identical(counting(.co_profile(pts, response = NULL, partition = "kmeans", nstart = 2,
                                        estimate = NULL, obs_per_cell = NULL, sac = NULL,
                                        n_levels = 3)),
                   character(0))
  expect_identical(counting(.co_profile(pts, response = NULL, partition = "kmeans", nstart = 2,
                                        estimate = NULL, obs_per_cell = c(18, 9), sac = "auto")),
                   character(0))
  # The density kind reads FNN for its fits, and nothing reads Matrix.
  expect_identical(counting(.co_profile(pts, response = NULL, estimate = NULL,
                                        obs_per_cell = c(18, 9))),
                   "FNN")
  # A response on three cells has no Moran step, which is what reads both.
  expect_identical(counting(.co_profile(pts, partition = "kmeans", nstart = 2, estimate = NULL,
                                        obs_per_cell = NULL, sac = 20000, levels = 3L)),
                   character(0))
  # An estimate reads FNN for each row's nearest centre.
  expect_identical(counting(.co_profile(pts, partition = "kmeans", nstart = 2, estimate = "pooled",
                                        obs_per_cell = NULL, sac = 20000, levels = 3L)),
                   "FNN")
  # A call the fit guard stops has asked only what the guard reads: nothing
  # for the kind "kmeans" that passes, and FNN for the condition it stops with.
  withr::local_options(spatialkit.fit_work_max = 1)
  expect_identical(counting(tryCatch(.co_profile(pts, response = NULL, partition = "kmeans",
                                                 nstart = 2, estimate = NULL, obs_per_cell = 9),
                                     error = function(e) NULL)),
                   "FNN")
})

test_that("without FNN the present ladder's steps are told so as well", {
  skip_if_not_installed("FNN")
  skip_if_not_installed("gstat")
  pts <- .co_600()
  call <- function() suppressWarnings(
    .co_profile(pts, estimate = "mean", obs_per_cell = NULL, sac = 20000, levels = c(20L, 90L)))
  with <- call()
  without <- local({
    .co_without_fnn()
    testthat::local_mocked_bindings(
      get.knn  = function(...) stop("FNN::get.knn() was reached"),
      get.knnx = function(...) stop("FNN::get.knnx() was reached"), .package = "FNN")
    call()
  })
  for (col in setdiff(names(with), c("moran_i", "moran_z")))
    expect_identical(without[[col]], with[[col]], info = col)
  expect_equal(without$moran_z, with$moran_z, tolerance = 1e-10)
})

test_that("without FNN the search of each row's nearest centre is bounded", {
  pts <- .co_600()
  call <- function(...) {
    .co_without_fnn()
    testthat::local_mocked_bindings(.partition_centres = function(...) stop("a fit was reached"))
    .co_error(.co_profile(pts, partition = "kmeans", nstart = 2, ...))
  }
  withr::local_options(spatialkit.exact_pairs_max = 600 * 133 - 1)
  e <- call(estimate = "pooled", obs_per_cell = c(18, 4.5))
  expect_identical(class(e), c("spatialkit_work_limit", "error", "condition"))
  # The scored rows times the centres of the finest level, before any fit.
  expect_identical(.co_fields(e), list(work = 600 * 133, limit = 600 * 133 - 1,
                                       what = "nearest centre"))
  expect_identical(conditionMessage(e),
                   paste0("resolution_profile(): finding the nearest of 133 centres for each of 600 ",
                          "rows would compare 8e+04 pairs without the FNN package, above the limit ",
                          "of 8e+04 (options(spatialkit.exact_pairs_max = )). Install FNN: the ",
                          "result is the same."))
  expect_null(conditionCall(e))
  # At the limit the search runs: the fit is reached.
  withr::local_options(spatialkit.exact_pairs_max = 600 * 133)
  expect_identical(conditionMessage(call(estimate = "pooled", obs_per_cell = c(18, 4.5))),
                   "a fit was reached")
  # Without an estimate no row is sent to a centre.
  withr::local_options(spatialkit.exact_pairs_max = 0)
  expect_identical(conditionMessage(call(estimate = NULL, obs_per_cell = c(18, 4.5))),
                   "a fit was reached")
  # On the ladder the function builds, the same search is bounded the same way.
  e <- call(estimate = "pooled", obs_per_cell = NULL, sac = 20000, levels = c(10L, 40L))
  expect_s3_class(e, "spatialkit_work_limit")
  expect_identical(e$what, "nearest centre")
  expect_identical(e$work, 600 * 40)
  # The rows counted are the scored rows, not the points: a response on half
  # of the rows, 300 rows in the 67 cells their size gives.
  half <- pts
  half$z[seq(2, 600, by = 2)] <- NA
  e <- local({
    .co_without_fnn()
    testthat::local_mocked_bindings(.partition_centres = function(...) stop("a fit was reached"))
    .co_error(.co_profile(half, partition = "kmeans", nstart = 2, estimate = "pooled",
                          obs_per_cell = 4.5))
  })
  expect_identical(.co_fields(e), list(work = 300 * 67, limit = 0, what = "nearest centre"))
  expect_match(conditionMessage(e), "finding the nearest of 67 centres for each of 300 rows",
               fixed = TRUE)
})


test_that("the nearest other location of a location is searched once for the whole profile", {
  skip_if_not_installed("FNN")
  pts  <- .co_600()
  real <- .nearest_location_search
  force(real)
  asked <- list()
  spy <- function(u, q, ...) {
    asked[[length(asked) + 1L]] <<- q
    real(u, q, ...)
  }
  run <- function(est, fnn = TRUE) {
    testthat::local_mocked_bindings(.nearest_location_search = spy)
    if (!fnn) .co_without_fnn()
    .co_profile(pts, estimate = est, obs_per_cell = c(9, 4.5, 2))
  }
  # With FNN: one search, of every location, before the sweep.
  with <- run("mean")
  expect_length(asked, 1L)
  expect_identical(sort(asked[[1L]]), 1:600)
  expect_gt(with$share_alone[3L], 0.05)
  # Without FNN: on first need, level by level, and no location twice.
  asked <- list()
  without <- run("mean", fnn = FALSE)
  expect_gt(length(asked), 1L)
  expect_identical(anyDuplicated(unlist(asked)), 0L)
  expect_lt(length(unlist(asked)), 600L)
  expect_identical(without$loo, with$loo)
  # A pooled value is predicted from its target: no nearest row is looked for.
  asked <- list()
  run("pooled")
  expect_length(asked, 0L)
})


# ---- The criterion against the error at new locations ---------------------------

test_that("loo tracks the squared error at new sites drawn like the data", {
  skip_on_cran()
  # Twenty fields.  1,600 sites uniform on a 1000 x 1000 square, x then y; an
  # exponential covariance with range parameter 100 (an effective range of
  # 300), sill 1 and nugget 0.1.  The first 600 sites are the layer and the
  # other 1,000 are new.
  nearest_row <- function(from, to)
    vapply(seq_len(nrow(from)), function(i)
      which.min((to[, 1] - from[i, 1])^2 + (to[, 2] - from[i, 2])^2), integer(1))
  ratio <- array(NA_real_, c(20L, 2L, 2L),
                 dimnames = list(NULL, c("mean", "pooled"), c("4.5", "2")))
  for (r in 1:20) {
    d <- local({
      .co_seed(1000 + r)
      x <- stats::runif(1600, 0, 1000); y <- stats::runif(1600, 0, 1000)
      D <- as.matrix(stats::dist(cbind(x, y)))
      z <- t(chol(exp(-D / 100) + diag(1e-10, 1600))) %*% stats::rnorm(1600) +
        stats::rnorm(1600, 0, sqrt(0.1))
      list(xy = cbind(x, y), z = as.numeric(z))
    })
    tr <- 1:600; te <- 601:1600
    pts <- .co_points(d$xy[tr, 1], d$xy[tr, 2], d$z[tr])
    for (size in c(4.5, 2)) for (est in c("mean", "pooled")) {
      p <- resolution_profile(pts, "z", estimate = est, obs_per_cell = size,
                              partition = "density", nstart = 3, sac = "auto", seed = 123)
      cen <- attr(p, "centres")[[1L]]
      rg  <- .co_regions(d$xy[tr, ], cen)
      fit <- .pool_fit(d$z[tr], rg$gi, rg$K, d$xy[tr, ], want = est)
      val <- if (est == "mean") fit$m else fit$value
      # A new site takes the value of the cell of its nearest centre; one
      # whose cell holds no row, that of the cell of its nearest scored row.
      at <- match(.nearest_centre(d$xy[te, 1], d$xy[te, 2], cen[, 1], cen[, 2]), rg$u)
      if (anyNA(at)) {
        empty <- which(is.na(at))
        at[empty] <- rg$gi[nearest_row(d$xy[te[empty], , drop = FALSE], d$xy[tr, ])]
      }
      ratio[r, est, as.character(size)] <- p$loo / mean((d$z[te] - val[at])^2) - 1
    }
  }
  expect_false(anyNA(ratio))
  mean_ratio <- apply(ratio, c(2, 3), mean)
  # The run this band was set from gave -2.9 and +0.7 percent for plain means
  # and -3.3 and -1.2 for pooled values, at 4.5 and 2 rows a cell.
  expect_true(all(abs(mean_ratio) < 0.10), info = paste(format(mean_ratio, digits = 3),
                                                        collapse = " "))
  # Single fields vary by more than the band (-17 to +15 percent in that
  # run), which is why the mean is asserted and no single field is.
})
