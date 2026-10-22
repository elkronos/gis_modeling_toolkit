# tests/testthat/test-fix-cp-predictor-path.R
# ---------------------------------------------------------------------------
# The path of resolution_profile() that has no estimate (predictor_vars, a
# geometry-only profile, `estimate = NULL`), read by Cp:
#
#   * Cp's noise variance never is a nugget at the fit's lower bound, and it
#     is held to the semivariance of nearest-neighbour pairs; attr(,
#     "cp_noise") says which source was used;
#   * the two variogram warnings carry classes;
#   * a cp selection carries a record of how the count was set, so a summary
#     of pooled values on its cells says the count was chosen for plain
#     means;
#   * the level chosen off a profile without an estimate fitted to a
#     subsample is refitted to every point before its seeds are handed on;
#     a profile with an estimate, scored on every row in the cells of the
#     centres fitted, hands those centres on;
#   * the print header names the kind fitted to a subsample;
#   * on simulated fields the cp count's held-out error is within a stated
#     regret of the best level's;
#   * the help and NEWS say what changed and what the ways back give.
#
# Every profile names `partition`, `nstart`, `estimate` and `obs_per_cell`;
# a variogram a test needs is made by hand, so that no number rests on how
# gstat's optimiser ends on a platform, except in the one test that reads
# the profile's own estimate.
# ---------------------------------------------------------------------------

.bp_seed <- function(seed, env = parent.frame())
  withr::local_seed(seed, .local_envir = env, .rng_kind = "Mersenne-Twister",
                    .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")

# `n` points on a 1000 m square in EPSG:32632, a response with a smooth
# trend plus noise of standard deviation `sd`, and two predictors.  Drawn
# under `seed`: x, y, the noise, then the predictors.
.bp_layer <- function(n, seed, sd = 0.5) {
  d <- local({
    .bp_seed(seed)
    x <- stats::runif(n, 0, 1000); y <- stats::runif(n, 0, 1000)
    e <- stats::rnorm(n, 0, sd)
    x1 <- stats::rnorm(n); x2 <- stats::rnorm(n)
    data.frame(x = 5e5 + x, y = 5e6 + y, x1 = x1, x2 = x2,
               z = 2 + 0.5 * x1 - 0.3 * x2 + sin(x / 150) + cos(y / 120) + e)
  })
  sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
}

# A variogram made by hand, with the attributes resolution_profile() reads as
# estimate_sac_range() sets them.  `range` is the effective range (NA for a
# refused one, with `reason`); `converged` goes on the model; `at_bound`,
# `short_lag` and `detrended` are left off when NULL.
.bp_sac <- function(range = 300, nugget = 0.3, psill = 1, reason = NULL, converged = TRUE,
                    at_bound = NULL, short_lag = NULL, short_dist = 25, detrended = TRUE) {
  vm <- data.frame(model = c("Nug", "Exp"), psill = c(nugget, psill),
                   range = c(0, if (is.na(range)) 3000 else range / 3),
                   stringsAsFactors = FALSE)
  if (!is.null(converged)) attr(vm, "converged") <- converged
  s <- structure(range, class = c("sac_range", "numeric"), variogram_model = vm,
                 nugget = nugget, detrended = detrended, crs = sf::st_crs(32632))
  if (!is.null(reason)) {
    attr(s, "rejected_range") <- 9000
    attr(s, "rejected_reason") <- reason
  }
  if (!is.null(at_bound)) attr(s, "nugget_at_bound") <- at_bound
  if (!is.null(short_lag)) {
    attr(s, "short_lag_semivariance") <- short_lag
    attr(s, "short_lag_pairs") <- 200L
    attr(s, "short_lag_dist") <- short_dist
  }
  s
}

# A profile of the predictor path, with its warnings kept.
.bp_profile <- function(pts, sac, ..., predictors = c("x1", "x2"), n_levels = 6L,
                        partition = "kmeans", nstart = 3L)
  .bp_catch(resolution_profile(pts, "z", predictor_vars = predictors, sac = sac,
                               n_levels = n_levels, partition = partition, nstart = nstart,
                               estimate = NULL, obs_per_cell = NULL, ...))

.bp_catch <- function(expr) {
  got <- list()
  value <- withCallingHandlers(expr, warning = function(w) {
    got[[length(got) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = got)
}
.bp_of <- function(r, class) Filter(function(w) inherits(w, class), r$warnings)
.bp_text <- function(r) vapply(r$warnings, conditionMessage, character(1))

# Cp of the profile at every level with the noise variance `tau2`, by the
# formula of ?resolution_profile on a profile with no subsample: every level's
# cells hold a scored row, so RSS / n + 2 tau2 L / n.
.bp_cp <- function(p, tau2) {
  n <- attr(p, "one_cell")$n
  p$rss / n + tau2 * 2 * p$levels / n
}

# The residual mean square of the finest level whose cells hold two scored
# rows on average, as the profile reads it: RSS / (n - L).
.bp_rms <- function(p) {
  n <- attr(p, "one_cell")$n
  ok <- is.finite(p$rss) & (n - p$levels) >= p$levels
  j <- which(ok)[which.max(p$levels[ok])]
  list(level = as.integer(p$levels[j]), value = p$rss[j] / (n - p$levels[j]))
}


# ---- A nugget at the fit's lower bound ----------------------------------------------

test_that("Cp does not take a nugget at the fit's lower bound as its noise variance", {
  pts <- .bp_layer(400, 2101)
  # A nugget of 0, marked as at its bound as estimate_sac_range() marks it,
  # with a range the variogram identifies: Cp used to take 0, and fall with
  # every finer level to the ceiling.
  for (sac in list(.bp_sac(nugget = 0, at_bound = TRUE),
                   # Unmarked: judged by the same rule, 1e-4 of the sill.
                   .bp_sac(nugget = 0),
                   .bp_sac(nugget = 6e-7, psill = 0.99))) {
    r <- .bp_profile(pts, sac)
    p <- r$value
    cn <- attr(p, "cp_noise")
    expect_identical(cn$source, "finest-level residual mean square")
    rms <- .bp_rms(p)
    expect_identical(cn$level, rms$level)
    expect_equal(cn$value, rms$value)
    expect_gt(cn$value, 0)
    expect_equal(p$cp, .bp_cp(p, cn$value))
    expect_true(all(p$cp > p$rss / attr(p, "one_cell")$n))
    expect_equal(attr(p, "one_cell")$cp,
                 attr(p, "one_cell")$rss / 400 + cn$value * 2 / 400)
    # The range still serves `reliability`, and the variogram is the sac's.
    expect_true(all(is.finite(p$reliability)))
    expect_equal(attr(p, "variogram")$nugget, sac_nugget(sac))
    # One classed warning says why, and no warning without a class is raised.
    w <- .bp_of(r, "spatialkit_cp_noise")
    expect_length(w, 1L)
    expect_identical(w[[1L]]$source, "finest-level residual mean square")
    expect_equal(w[[1L]]$value, cn$value)
    expect_identical(w[[1L]]$level, cn$level)
    expect_identical(w[[1L]]$uncapped_source, "finest-level residual mean square")
    expect_true(is.na(w[[1L]]$cap))
    expect_match(conditionMessage(w[[1L]]), "is at the fit's lower bound, which is not an estimate",
                 fixed = TRUE)
    expect_match(conditionMessage(w[[1L]]),
                 "so `cp` takes its noise variance from the residual mean square of the finest level",
                 fixed = TRUE)
    expect_false(any(grepl("reliability` needs", .bp_text(r), fixed = TRUE)))
    expect_length(r$warnings, 1L)
    expect_false(any(grepl("nugget is 0", .bp_text(r), fixed = TRUE)))
    # The print says the nugget is at its bound and where Cp's noise came from.
    expect_match(paste(utils::capture.output(print(p)), collapse = "\n"),
                 sprintf(paste0("(at the fit's lower bound), partial sill %s, range 300; cp's ",
                                "noise variance %s is the residual mean square at %d cells"),
                         format(signif(attr(p, "variogram")$psill, 3)),
                         format(signif(cn$value, 3)), cn$level),
                 fixed = TRUE)
  }

  # The same nugget behind a range refused as past the lags: Cp takes
  # neither, and the warning on the range says so.
  r <- .bp_profile(pts, .bp_sac(NA_real_, nugget = 0, at_bound = TRUE,
                                reason = "fitted range exceeds the largest lag fitted"))
  p <- r$value
  expect_identical(attr(p, "cp_noise")$source, "finest-level residual mean square")
  expect_true(all(is.na(p$reliability)))
  u <- .bp_of(r, "spatialkit_variogram_unusable")
  expect_length(u, 1L)
  expect_identical(u[[1L]]$nugget_used, FALSE)
  expect_identical(u[[1L]]$converged, TRUE)
  expect_match(conditionMessage(u[[1L]]),
               "its nugget (0) is at the fit's lower bound, which is not an estimate", fixed = TRUE)
  expect_false(any(grepl("still takes its nugget", .bp_text(r), fixed = TRUE)))
  expect_length(.bp_of(r, "spatialkit_cp_noise"), 1L)

  # A nugget above 1e-4 of the sill is taken as it stands.
  r <- .bp_profile(pts, .bp_sac(nugget = 0.2))
  expect_identical(attr(r$value, "cp_noise"),
                   list(source = "variogram nugget", value = 0.2, level = NA_integer_))
  expect_length(r$warnings, 0L)
})

test_that("on simulated fields with no noise Cp never takes the nugget estimate_sac_range() puts at its bound", {
  skip_if_not_installed("gstat")
  # Smooth fields with no noise: the true nugget is 0, and the fit's nugget
  # often runs into its bound.  Every profile on the predictor path reads its
  # own variogram.
  at_bound <- 0L
  for (seed in 1:6) {
    pts <- local({
      .bp_seed(3300 + seed)
      n <- 300
      x <- stats::runif(n, 0, 1000); y <- stats::runif(n, 0, 1000)
      D <- as.matrix(stats::dist(cbind(x, y)))
      f <- as.numeric(t(chol(exp(-D / 120) + diag(1e-8, n))) %*% stats::rnorm(n))
      x1 <- stats::rnorm(n); x2 <- stats::rnorm(n)
      sf::st_as_sf(data.frame(x = 5e5 + x, y = 5e6 + y, x1 = x1, x2 = x2,
                              z = 1 + 0.5 * x1 - 0.3 * x2 + f),
                   coords = c("x", "y"), crs = 32632)
    })
    p <- suppressWarnings(resolution_profile(pts, "z", predictor_vars = c("x1", "x2"),
                                             n_levels = 6L, nstart = 3L, partition = "kmeans",
                                             estimate = NULL, obs_per_cell = NULL))
    sac <- attr(p, "sac")
    cn <- attr(p, "cp_noise")
    if (isTRUE(attr(sac, "nugget_at_bound"))) {
      at_bound <- at_bound + 1L
      expect_false(identical(cn$source, "variogram nugget"))
    }
    # Whatever its source, the noise variance is positive and no more than
    # the semivariance of nearest-neighbour pairs.
    if (!is.null(cn)) {
      expect_gt(cn$value, 0)
      expect_lte(cn$value, attr(sac, "short_lag_semivariance"))
    }
  }
  skip_if_not(at_bound > 0L, "no fit put the nugget at its bound on this platform")
  expect_gt(at_bound, 0L)
})


# ---- The semivariance of nearest-neighbour pairs ---------------------------------------

test_that("Cp's noise variance is held to the short-lag semivariance, and cp_noise says which source", {
  pts <- .bp_layer(400, 2101)
  # A nugget above the semivariance of the closest pairs.
  r <- .bp_profile(pts, .bp_sac(nugget = 0.5, short_lag = 0.3, short_dist = 25))
  p <- r$value
  cn <- attr(p, "cp_noise")
  expect_identical(cn, list(source = "short-lag semivariance", value = 0.3, level = NA_integer_,
                            cap = 0.3,
                            uncapped = list(source = "variogram nugget", value = 0.5,
                                            level = NA_integer_)))
  expect_equal(p$cp, .bp_cp(p, 0.3))
  expect_equal(attr(p, "one_cell")$cp, attr(p, "one_cell")$rss / 400 + 0.3 * 2 / 400)
  w <- .bp_of(r, "spatialkit_cp_noise")
  expect_length(w, 1L)
  expect_length(r$warnings, 1L)
  expect_identical(w[[1L]][c("source", "value", "level", "cap", "uncapped_source",
                             "uncapped_value")],
                   list(source = "short-lag semivariance", value = 0.3, level = NA_integer_,
                        cap = 0.3, uncapped_source = "variogram nugget", uncapped_value = 0.5))
  expect_match(conditionMessage(w[[1L]]),
               paste0("the variogram's nugget (0.5) is above the semivariance of ",
                      "nearest-neighbour pairs (0.3 at a median separation of 25), which ",
                      "bounds the nugget from above, so `cp` takes that semivariance"),
               fixed = TRUE)
  expect_match(paste(utils::capture.output(print(p)), collapse = "\n"),
               paste0("nugget 0.5, partial sill 1, range 300; cp's noise variance 0.3 is the ",
                      "semivariance of nearest-neighbour pairs, below the nugget (0.5)"),
               fixed = TRUE)
  # A semivariance above the nugget holds nothing, and neither does one of 0.
  for (sl in c(0.8, 0)) {
    r <- .bp_profile(pts, .bp_sac(nugget = 0.5, short_lag = sl))
    expect_identical(attr(r$value, "cp_noise"),
                     list(source = "variogram nugget", value = 0.5, level = NA_integer_))
    expect_equal(r$value$cp, .bp_cp(r$value, 0.5))
    expect_length(r$warnings, 0L)
  }

  # The residual mean square of the finest level, held the same way: a model
  # that did not converge gives Cp no nugget.
  r <- .bp_profile(pts, .bp_sac(NA_real_, nugget = 0.5, converged = FALSE, short_lag = 0.05,
                                reason = "fitted range exceeds the largest lag fitted"))
  p <- r$value
  cn <- attr(p, "cp_noise")
  rms <- .bp_rms(p)
  expect_gt(rms$value, 0.05)
  expect_identical(cn[c("source", "value", "level", "cap")],
                   list(source = "short-lag semivariance", value = 0.05, level = NA_integer_,
                        cap = 0.05))
  expect_identical(cn$uncapped$source, "finest-level residual mean square")
  expect_identical(cn$uncapped$level, rms$level)
  expect_equal(cn$uncapped$value, rms$value)
  expect_equal(p$cp, .bp_cp(p, 0.05))
  w <- .bp_of(r, "spatialkit_cp_noise")
  expect_length(w, 1L)
  expect_identical(w[[1L]]$uncapped_source, "finest-level residual mean square")
  expect_match(conditionMessage(w[[1L]]),
               paste0("so `cp` takes its noise variance from the semivariance of ",
                      "nearest-neighbour pairs (0.05 at a median separation of 25)"),
               fixed = TRUE)
  expect_match(conditionMessage(w[[1L]]),
               sprintf("in place of the residual mean square of the finest level (%d cells: ",
                       rms$level), fixed = TRUE)
  expect_match(paste(utils::capture.output(print(p)), collapse = "\n"),
               sprintf(paste0("none usable (reliability is NA; cp's noise variance 0.05 is the ",
                              "semivariance of nearest-neighbour pairs, below the residual mean ",
                              "square at %d cells (%s))"),
                       rms$level, format(signif(rms$value, 3))),
               fixed = TRUE)
  # select_resolution() reads the column computed with it.
  s <- select_resolution(p, criterion = "cp")
  expect_identical(s$value, min(p$cp))
  expect_equal(s$one_cell, attr(p, "one_cell")$cp)
})


# ---- The two variogram warnings ---------------------------------------------------------

test_that("the two variogram warnings of the profile have classes and fields", {
  pts <- .bp_layer(400, 2101)
  past <- "fitted range exceeds the largest lag fitted"
  short <- "fitted range is below the shortest lag fitted"
  cases <- list(
    list(sac = .bp_sac(NA_real_, reason = past, converged = FALSE), used = FALSE, conv = FALSE),
    list(sac = .bp_sac(NA_real_, reason = short, converged = TRUE), used = FALSE, conv = TRUE),
    list(sac = .bp_sac(NA_real_, reason = past, converged = TRUE), used = TRUE, conv = TRUE),
    list(sac = .bp_sac(NA_real_, reason = past, converged = NULL), used = TRUE, conv = NA))
  for (cs in cases) {
    r <- .bp_profile(pts, cs$sac)
    u <- .bp_of(r, "spatialkit_variogram_unusable")
    expect_length(u, 1L)
    expect_identical(u[[1L]]$reason, attr(cs$sac, "rejected_reason"))
    expect_identical(u[[1L]]$converged, cs$conv)
    expect_identical(u[[1L]]$nugget, 0.3)
    expect_identical(u[[1L]]$nugget_used, cs$used)
    expect_match(conditionMessage(u[[1L]]), "reports no usable range", fixed = TRUE)
    # Where the nugget is not used, the noise variance comes from the
    # finest level, and that is the second warning.
    expect_length(.bp_of(r, "spatialkit_cp_noise"), if (cs$used) 0L else 1L)
    expect_identical(attr(r$value, "cp_noise")$source,
                     if (cs$used) "variogram nugget" else "finest-level residual mean square")
    # Every warning of the profile has a class now.
    expect_true(all(vapply(r$warnings, function(w)
      !inherits(w, "simpleWarning") && any(startsWith(class(w), "spatialkit_")), logical(1))))
  }
  # A handler for the class catches each, as the help says.
  got <- tryCatch(resolution_profile(pts, "z", predictor_vars = c("x1", "x2"),
                                     sac = .bp_sac(NA_real_, reason = past, converged = FALSE),
                                     n_levels = 4L, nstart = 2L, partition = "kmeans",
                                     estimate = NULL, obs_per_cell = NULL),
                  spatialkit_variogram_unusable = function(w) "unusable")
  expect_identical(got, "unusable")
  got <- withCallingHandlers(
    tryCatch(resolution_profile(pts, "z", predictor_vars = c("x1", "x2"),
                                sac = .bp_sac(NA_real_, reason = past, converged = FALSE),
                                n_levels = 4L, nstart = 2L, partition = "kmeans",
                                estimate = NULL, obs_per_cell = NULL),
             spatialkit_cp_noise = function(w) "noise"),
    spatialkit_variogram_unusable = function(w) invokeRestart("muffleWarning"))
  expect_identical(got, "noise")
})


# ---- The record of a cp count -------------------------------------------------------------

test_that("a count chosen by cp off the predictor path carries a record, and pooled values on it are warned about", {
  pts <- .bp_layer(500, 2102)
  p <- .bp_profile(pts, .bp_sac(nugget = 0.3), partition = "density")$value
  expect_null(attr(p, "estimate"))
  sel <- select_resolution(p, criterion = "cp")
  rec <- sel$count_rule
  expect_s3_class(rec, "spatialkit_count_rule")
  expect_identical(rec$rule, "cp")
  expect_identical(rec$estimate, "mean")
  expect_identical(rec$cells, sel$best)
  expect_null(rec$scored)
  expect_identical(rec$response, NA_character_)
  expect_identical(rec$n, 500L)
  expect_identical(rec$n_fit, 500L)
  expect_identical(rec$scored_cells, NA)
  expect_identical(rec$partition$kind, "density")
  # The selection gains the record alone, after its seeds.
  expect_identical(names(sel),
                   c("best", "flat", "criterion", "value", "at_ceiling", "at_floor", "edge",
                     "n_levels", "values", "structured", "one_cell", "structure_p",
                     "range_floor", "note", "seeds", "count_rule", "partition"))
  on_seeds <- rec; on_seeds$scored_cells <- TRUE
  expect_identical(attr(sel$seeds, "count_rule"), on_seeds)

  # From the selection, and from the profile passed as the count itself.
  for (count in list(sel, p)) {
    seeds <- suppressWarnings(get_voronoi_seeds(method = "kmeans", n = count, sample_points = pts))
    expect_identical(attr(seeds, "count_rule"), on_seeds)
    tess <- build_tessellation(seeds, method = "voronoi", quiet = TRUE)
    expect_identical(tess$params$count_rule, on_seeds)
    asg <- assign_features_to_polygons(pts, tess$cells, polygon_id_col = "cell_id", quiet = TRUE)
    pooled <- .bp_catch(summarize_by_cell(asg, "z", id_col = "cell_id", cells_sf = tess$cells,
                                          estimate = "pooled"))
    w <- .bp_of(pooled, "spatialkit_estimate_mismatch")
    expect_length(w, 1L)
    expect_identical(w[[1L]]$chosen_for, "mean")
    expect_identical(w[[1L]]$reported, "pooled")
    expect_identical(w[[1L]]$rule, "cp")
    expect_identical(w[[1L]]$cells, sel$best)
    expect_match(conditionMessage(w[[1L]]),
                 sprintf("the number of cells (%d, set by cp) was chosen for the plain mean",
                         sel$best), fixed = TRUE)
    expect_identical(attr(pooled$value, "count_rule")$check,
                     list(estimate = "mismatch",
                          cells = "not checked: the profile scored no estimate"))
    plain <- .bp_catch(summarize_by_cell(asg, "z", id_col = "cell_id", cells_sf = tess$cells,
                                         estimate = "mean"))
    expect_length(.bp_of(plain, "spatialkit_estimate_mismatch"), 0L)
    expect_identical(attr(plain$value, "count_rule")$check$estimate, "match")
  }
  expect_identical(attr(suppressWarnings(voronoi_seeds_kmeans(pts, sel)), "count_rule"),
                   on_seeds)
  # The record's print says what the count was chosen for.
  expect_match(paste(utils::capture.output(print(rec)), collapse = "\n"),
               "chosen for  : the plain mean", fixed = TRUE)

  # A diagnostic read off the same profile still carries none.
  z <- suppressWarnings(select_resolution(p, criterion = "moran_z", diagnostic = TRUE))
  expect_null(z$count_rule)
  expect_null(attr(z$seeds, "count_rule"))
})


# ---- The level chosen off a subsample, refitted to every point ------------------------------

test_that("the level chosen off a profile fitted to a subsample is refitted to every point before its seeds are handed on", {
  pts <- .bp_layer(1200, 2103)
  for (kind in c("kmeans", "density")) {
    info <- kind
    p <- suppressWarnings(
      resolution_profile(pts, "z", predictor_vars = c("x1", "x2"), sac = .bp_sac(nugget = 0.3),
                         sample_n = 300, n_levels = 6L, partition = kind, nstart = 3L,
                         estimate = NULL, obs_per_cell = NULL))
    b <- attr(p, "bounds")
    expect_identical(b$n_sample, 300L, info = info)
    rows <- attr(p, "refit_rows")
    expect_named(rows, c("xy", "iter_max"))
    expect_identical(dim(rows$xy), c(1200L, 2L), info = info)
    expect_identical(rows$iter_max, if (kind == "kmeans") 50L else 100L, info = info)
    # The coordinates of every point, as the profile holds them.
    xy <- unname(sf::st_coordinates(pts)[, 1:2])
    expect_identical(rows$xy, xy[order(xy[, 1], xy[, 2]), ], info = info)

    withr::local_seed(99)
    rng <- get(".Random.seed", envir = globalenv())
    sel <- select_resolution(p, criterion = "cp")
    # No random number is drawn.
    expect_identical(get(".Random.seed", envir = globalenv()), rng, info = info)
    scored <- attr(p, "centres")[[as.character(sel$best)]]
    got <- unname(sf::st_coordinates(sel$seeds)[, 1:2])
    expect_false(isTRUE(all.equal(got, unname(scored))), info = info)
    # The kind's own iteration from the centres scored, on every point.
    want <- if (kind == "kmeans")
      stats::kmeans(rows$xy, centers = scored, iter.max = 50L)$centers
    else
      .density_kmeans(rows$xy, nrow(scored), power = 1, k_nn = 8L, iter.max = 100L,
                      centers = scored)$centers
    expect_equal(got, unname(want[, 1:2]), tolerance = 1e-12, info = info)
    rf <- attr(sel$seeds, "refit")
    expect_named(rf, c("rows", "subsample", "kind", "iter", "iter_max", "converged",
                       "wss_scored", "wss"))
    expect_identical(rf$rows, 1200L, info = info)
    expect_identical(rf$subsample, 300L, info = info)
    expect_identical(rf$kind, kind, info = info)
    expect_true(rf$converged, info = info)
    # Every point in the cells of the centres scored, about each cell's mean,
    # against the refitted cells.
    a0 <- .nearest_exact(rows$xy[, 1], rows$xy[, 2], scored[, 1], scored[, 2])
    m0 <- rowsum(rows$xy, a0) / as.numeric(table(a0))
    expect_equal(rf$wss_scored, sum((rows$xy - m0[as.character(a0), ])^2), info = info)
    if (kind == "kmeans") expect_lte(rf$wss, rf$wss_scored)
    # The record says these are other cells built from the count.
    expect_identical(attr(sel$seeds, "count_rule")$scored_cells, FALSE, info = info)
    expect_identical(sel$count_rule$scored_cells, NA, info = info)
    # The seeding functions hand on the refitted centres, with what was done.
    seeds <- suppressWarnings(get_voronoi_seeds(method = "kmeans", n = sel, sample_points = pts))
    expect_equal(unname(sf::st_coordinates(seeds)[, 1:2]), got, tolerance = 1e-9, info = info)
    expect_identical(attr(seeds, "refit"), rf, info = info)
    expect_identical(attr(seeds, "count_rule")$scored_cells, FALSE, info = info)
    seeds_p <- suppressWarnings(get_voronoi_seeds(method = "kmeans", n = p, sample_points = pts))
    expect_identical(attr(seeds_p, "refit"), rf, info = info)
    seeds_k <- suppressWarnings(voronoi_seeds_kmeans(pts, sel))
    expect_identical(attr(seeds_k, "refit"), rf, info = info)
    expect_identical(attr(seeds_k, "count_rule")$scored_cells, FALSE, info = info)
    # The print of the selection says so.
    expect_match(paste(utils::capture.output(print(sel)), collapse = " "),
                 "seeds       : refitted to all 1200 points from the centres scored on a",
                 fixed = TRUE)
  }

  # summary() and plot() read no seeds and refit nothing; a selection does.
  n_refit <- 0L
  real <- .refit_level
  local_mocked_bindings(.refit_level = function(...) {
    n_refit <<- n_refit + 1L
    real(...)
  })
  invisible(summary(p))
  expect_identical(n_refit, 0L)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    invisible(plot(p))
    expect_identical(n_refit, 0L)
  }
  bnd <- sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(pts)))
  invisible(suppressWarnings(build_tessellation(pts, boundary = bnd, method = "hex",
                                                approx_n_cells = p, quiet = TRUE)))
  expect_identical(n_refit, 0L)
  invisible(select_resolution(p, criterion = "cp"))
  expect_identical(n_refit, 1L)
})

test_that("a refit above the work limit is not run, and a profile on every point has nothing to refit", {
  pts <- .bp_layer(1200, 2103)
  p <- suppressWarnings(
    resolution_profile(pts, "z", predictor_vars = c("x1", "x2"), sac = .bp_sac(nugget = 0.3),
                       sample_n = 300, n_levels = 6L, partition = "kmeans", nstart = 3L,
                       estimate = NULL, obs_per_cell = NULL))
  withr::local_options(spatialkit.fit_work_max = 10)
  sel <- select_resolution(p, criterion = "cp")
  expect_null(attr(sel$seeds, "refit"))
  expect_equal(unname(sf::st_coordinates(sel$seeds)[, 1:2]),
               unname(attr(p, "centres")[[as.character(sel$best)]]))
  expect_identical(attr(sel$seeds, "count_rule")$scored_cells, TRUE)
  withr::local_options(spatialkit.fit_work_max = NULL)

  # Fitted to every point: no coordinates kept, the centres scored handed on.
  q <- suppressWarnings(
    resolution_profile(pts[1:500, ], "z", predictor_vars = c("x1", "x2"),
                       sac = .bp_sac(nugget = 0.3), n_levels = 6L, partition = "kmeans",
                       nstart = 3L, estimate = NULL, obs_per_cell = NULL))
  expect_null(attr(q, "refit_rows"))
  s <- select_resolution(q, criterion = "cp")
  expect_null(attr(s$seeds, "refit"))
  expect_equal(unname(sf::st_coordinates(s$seeds)[, 1:2]),
               unname(attr(q, "centres")[[as.character(s$best)]]))
})

test_that("a density refit reads the profile's density_power and density_k", {
  pts <- .bp_layer(1200, 2103)
  p <- suppressWarnings(
    resolution_profile(pts, "z", predictor_vars = c("x1", "x2"), sac = .bp_sac(nugget = 0.3),
                       sample_n = 300, n_levels = 6L, partition = "density",
                       density_power = 0.5, density_k = 4, nstart = 3L,
                       estimate = NULL, obs_per_cell = NULL))
  expect_identical(attr(p, "partition")$density_power, 0.5)
  expect_identical(attr(p, "partition")$density_k, 4L)
  sel <- select_resolution(p, criterion = "cp")
  rows <- attr(p, "refit_rows")
  scored <- attr(p, "centres")[[as.character(sel$best)]]
  got <- unname(sf::st_coordinates(sel$seeds)[, 1:2])
  refit <- function(power, k_nn)
    unname(.density_kmeans(rows$xy, nrow(scored), power = power, k_nn = k_nn,
                           iter.max = 100L, centers = scored)$centers)
  expect_equal(got, refit(0.5, 4L), tolerance = 1e-12)
  # The weights of the default settings give other centres, so the two
  # settings are read.
  expect_false(isTRUE(all.equal(got, refit(1, 4L), tolerance = 1e-8)))
  expect_false(isTRUE(all.equal(got, refit(0.5, 8L), tolerance = 1e-8)))
  expect_identical(attr(sel$seeds, "refit")$kind, "density")
})

test_that("where Hartigan and Wong refuse the centres scored, Lloyd's iteration with every weight one refits them, and the log says so", {
  pts <- .bp_layer(1200, 2103)
  p <- suppressWarnings(
    resolution_profile(pts, "z", predictor_vars = c("x1", "x2"), sac = .bp_sac(nugget = 0.3),
                       sample_n = 300, n_levels = 6L, partition = "kmeans", nstart = 3L,
                       estimate = NULL, obs_per_cell = NULL))
  L <- as.character(.without_refit(select_resolution(p, criterion = "cp"))$best)
  # Two centres in one place: stats::kmeans() refuses such a start.
  cen <- attr(p, "centres")
  cen[[L]][2L, ] <- cen[[L]][1L, ]
  attr(p, "centres") <- cen
  expect_error(stats::kmeans(attr(p, "refit_rows")$xy, centers = cen[[L]], iter.max = 50L),
               "not distinct")
  lines <- capture_spatialkit_log(sel <- select_resolution(p, criterion = "cp"))
  expect_identical(as.character(sel$best), L)
  want <- .density_kmeans(attr(p, "refit_rows")$xy, nrow(cen[[L]]), power = 0, k_nn = 8L,
                          iter.max = 50L, centers = cen[[L]])$centers
  expect_equal(unname(sf::st_coordinates(sel$seeds)[, 1:2]), unname(want), tolerance = 1e-12)
  # Not the weights of the density kind.
  dens <- .density_kmeans(attr(p, "refit_rows")$xy, nrow(cen[[L]]), power = 1, k_nn = 8L,
                          iter.max = 50L, centers = cen[[L]])$centers
  expect_false(isTRUE(all.equal(unname(sf::st_coordinates(sel$seeds)[, 1:2]), unname(dens),
                                tolerance = 1e-8)))
  expect_identical(attr(sel$seeds, "refit")$kind, "kmeans")
  expect_true(any(grepl(sprintf(paste0("Hartigan and Wong's algorithm refuses the %s centres ",
                                       "scored as a start, so they are refitted to the 1200 ",
                                       "points by Lloyd's iteration in its place"), L),
                        lines, fixed = TRUE)))
})

test_that("without FNN a density refit whose weights would compare more pairs than the limit is not run", {
  pts <- .bp_layer(1200, 2103)
  p <- suppressWarnings(
    resolution_profile(pts, "z", predictor_vars = c("x1", "x2"), sac = .bp_sac(nugget = 0.3),
                       sample_n = 300, n_levels = 6L, partition = "density", nstart = 3L,
                       estimate = NULL, obs_per_cell = NULL))
  local_mocked_bindings(.has_backend = function(pkg) FALSE)
  rows <- attr(p, "refit_rows")
  # 1200 locations: 1.44 million pairs.
  withr::local_options(spatialkit.exact_pairs_max = 1.4e6)
  lines <- capture_spatialkit_log(sel <- select_resolution(p, criterion = "cp"))
  scored <- attr(p, "centres")[[as.character(sel$best)]]
  expect_null(attr(sel$seeds, "refit"))
  expect_equal(unname(sf::st_coordinates(sel$seeds)[, 1:2]), unname(scored))
  expect_identical(attr(sel$seeds, "count_rule")$scored_cells, TRUE)
  expect_true(any(grepl("would pass a limit of the work a fit may do", lines, fixed = TRUE)))
  # Just above the pairs, the exact search runs.
  withr::local_options(spatialkit.exact_pairs_max = 1.5e6)
  sel <- select_resolution(p, criterion = "cp")
  expect_identical(attr(sel$seeds, "refit")$rows, 1200L)
  expect_equal(unname(sf::st_coordinates(sel$seeds)[, 1:2]),
               unname(.density_kmeans(rows$xy, nrow(scored), power = 1, k_nn = 8L,
                                      iter.max = 100L, use_fnn = FALSE,
                                      centers = scored)$centers),
               tolerance = 1e-12)
})

test_that("a profile saved before the refit keeps no coordinates: its centres are handed on as scored, and the log says why", {
  pts <- .bp_layer(1200, 2103)
  p <- suppressWarnings(
    resolution_profile(pts, "z", predictor_vars = c("x1", "x2"), sac = .bp_sac(nugget = 0.3),
                       sample_n = 300, n_levels = 6L, partition = "kmeans", nstart = 3L,
                       estimate = NULL, obs_per_cell = NULL))
  attr(p, "refit_rows") <- NULL
  lines <- capture_spatialkit_log(sel <- select_resolution(p, criterion = "cp"))
  expect_null(attr(sel$seeds, "refit"))
  expect_equal(unname(sf::st_coordinates(sel$seeds)[, 1:2]),
               unname(attr(p, "centres")[[as.character(sel$best)]]))
  expect_identical(attr(sel$seeds, "count_rule")$scored_cells, TRUE)
  expect_true(any(grepl(paste0("this profile was fitted to a subsample of 300 of 1200 points and ",
                               "keeps no coordinates of its points"), lines, fixed = TRUE)))
  # summary() and plot() read no seeds and say nothing of it.
  quiet <- capture_spatialkit_log(invisible(summary(p)))
  expect_false(any(grepl("keeps no coordinates", quiet, fixed = TRUE)))
})


# ---- The bound and the fallback, at their edges ----------------------------------------------

test_that("a nugget at its bound with no level to take a noise variance from leaves cp NA at every level, with a classed warning", {
  pts <- .bp_layer(400, 2101)
  # Every level holds fewer than two scored rows a cell.
  r <- .bp_catch(resolution_profile(pts, "z", predictor_vars = c("x1", "x2"),
                                    sac = .bp_sac(nugget = 0, at_bound = TRUE),
                                    levels = c(250, 300), partition = "kmeans", nstart = 3L,
                                    estimate = NULL, obs_per_cell = NULL))
  p <- r$value
  expect_identical(p$levels, c(250L, 300L))
  expect_true(all(is.na(p$cp)))
  expect_true(is.na(attr(p, "one_cell")$cp))
  expect_null(attr(p, "cp_noise"))
  # The range still serves `reliability`.
  expect_true(all(is.finite(p$reliability)))
  expect_length(r$warnings, 1L)
  w <- .bp_of(r, "spatialkit_cp_noise")
  expect_length(w, 1L)
  expect_true(is.na(w[[1L]]$source))
  expect_true(is.na(w[[1L]]$value))
  expect_true(is.na(w[[1L]]$level))
  expect_identical(w[[1L]]$uncapped_source, "variogram nugget")
  expect_identical(w[[1L]]$uncapped_value, 0)
  expect_match(conditionMessage(w[[1L]]),
               paste0("no level's cells hold at least two scored rows on average to take a ",
                      "noise variance from, so `cp` is NA at every level"), fixed = TRUE)
  # The print marks the nugget and says why `cp` is NA, as the warning does
  # (BUG-07): it printed the nugget with no mark and nothing on `cp`.
  out <- utils::capture.output(print(p))
  expect_true(any(grepl(paste0("  variogram   : nugget 0 (at the fit's lower bound), partial sill 1, ",
                               "range 300; cp is NA at every level: no level's cells hold two scored ",
                               "rows on average"), out, fixed = TRUE)))
  # A geometry-only profile with the same `sac` has no `cp` to speak of: the
  # mark alone.
  g <- suppressWarnings(resolution_profile(pts, NULL, sac = .bp_sac(nugget = 0, at_bound = TRUE),
                                           levels = c(250, 300), partition = "kmeans",
                                           nstart = 3L, estimate = NULL, obs_per_cell = NULL))
  expect_null(attr(g, "cp_noise"))
  out <- utils::capture.output(print(g))
  expect_true(any(grepl("  variogram   : nugget 0 (at the fit's lower bound), partial sill 1, range 300",
                        out, fixed = TRUE)))
  expect_false(any(grepl("cp is NA at every level", out, fixed = TRUE)))
  # A constant response with the same `sac`: the mark alone, since `cp` is NA
  # there for the reason the warning on a constant response gives, not for
  # want of a level to take a noise variance from.
  const <- pts; const$z <- 5
  k <- .bp_catch(resolution_profile(const, "z", sac = .bp_sac(nugget = 0, at_bound = TRUE),
                                    levels = c(250, 300), partition = "kmeans", nstart = 3L,
                                    estimate = NULL, obs_per_cell = NULL))
  expect_length(.bp_of(k, "spatialkit_constant_response"), 1L)
  expect_null(attr(k$value, "cp_noise"))
  expect_true(attr(k$value, "one_cell")$constant)
  out <- utils::capture.output(print(k$value))
  expect_true(any(grepl("  variogram   : nugget 0 (at the fit's lower bound), partial sill 1, range 300",
                        out, fixed = TRUE)))
  expect_false(any(grepl("cp is NA at every level", out, fixed = TRUE)))
})

test_that("the bound is a nugget of at most 1e-4 of the sill: just under it the nugget is refused, just over it taken", {
  pts <- .bp_layer(400, 2101)
  # Unmarked, and marked FALSE: the rule decides either way.
  for (at_bound in list(NULL, FALSE)) {
    under <- .bp_profile(pts, .bp_sac(nugget = 0.99e-4, psill = 1, at_bound = at_bound))
    cn <- attr(under$value, "cp_noise")
    expect_identical(cn$source, "finest-level residual mean square")
    expect_identical(cn$level, .bp_rms(under$value)$level)
    expect_gt(cn$value, 0.1)
    expect_length(.bp_of(under, "spatialkit_cp_noise"), 1L)
    over <- .bp_profile(pts, .bp_sac(nugget = 1.01e-4, psill = 1, at_bound = at_bound))
    expect_identical(attr(over$value, "cp_noise"),
                     list(source = "variogram nugget", value = 1.01e-4, level = NA_integer_))
    expect_length(over$warnings, 0L)
  }
  # A mark says the bound whatever the share.
  marked <- .bp_profile(pts, .bp_sac(nugget = 1.01e-4, psill = 1, at_bound = TRUE))
  expect_identical(attr(marked$value, "cp_noise")$source, "finest-level residual mean square")
})


# ---- The count against held-out error -------------------------------------------------------

test_that("on simulated fields with a known nugget and range the cp count's held-out error is within 15 percent of the best level's", {
  skip_if_not_installed("gstat")
  # The item's test (GAP-11).  Each layer: a Gaussian field of exponential
  # covariance with the given range on a 10 km square, interpolated from a
  # 30 x 30 grid, plus white noise of the given variance (the nugget) and
  # two predictors.  2,500 points: 2,000 to profile, more than `sample_n`,
  # on the default predictor path, and 500 held out.  A level's held-out
  # error is that of the OLS fit plus its cell's mean training residual, the
  # cells being those of the centres the profile scored at that level.  The
  # stated regret: the cp count's error at most 15 percent above the least
  # error of any level.  On these six layers it is 1.0, 0, 0, 2.7, 0 and 9.9
  # percent.
  sim <- function(seed, nugget, range) {
    .bp_seed(seed)
    n <- 2500L; ext <- 1e4
    x <- stats::runif(n, 0, ext); y <- stats::runif(n, 0, ext)
    g <- as.matrix(expand.grid(gx = seq(0, ext, length.out = 30),
                               gy = seq(0, ext, length.out = 30)))
    f0 <- as.numeric(t(chol(exp(-as.matrix(stats::dist(g)) / range) + diag(1e-6, nrow(g)))) %*%
                       stats::rnorm(nrow(g)))
    nn <- FNN::get.knnx(g, cbind(x, y), 4L)
    w <- 1 / (nn$nn.dist + 1)
    f <- rowSums(w * matrix(f0[nn$nn.index], ncol = 4L)) / rowSums(w)
    x1 <- stats::rnorm(n); x2 <- stats::rnorm(n)
    z <- 1 + 0.5 * x1 - 0.3 * x2 + f + stats::rnorm(n, sd = sqrt(nugget))
    test <- sample.int(n, 500L)
    sf::st_as_sf(data.frame(x = 5e5 + x, y = 5e6 + y, x1 = x1, x2 = x2, z = z,
                            held_out = seq_len(n) %in% test),
                 coords = c("x", "y"), crs = 32632)
  }
  layers <- list(c(nugget = 0.5, range = 1500), c(0.1, 3000), c(0, 2000), c(1, 800),
                 c(0.02, 4000), c(0.3, 2500))
  regret <- numeric(0)
  for (i in seq_along(layers)) {
    pts <- sim(4400 + i, layers[[i]][[1L]], layers[[i]][[2L]])
    tr <- pts[!pts$held_out, ]; te <- pts[pts$held_out, ]
    p <- suppressWarnings(resolution_profile(tr, "z", c("x1", "x2")))
    expect_identical(attr(p, "bounds")$n_sample, 1500L)
    sel <- suppressWarnings(select_resolution(p, criterion = "cp"))
    fit <- stats::lm(z ~ x1 + x2, data = sf::st_drop_geometry(tr))
    res <- stats::residuals(fit)
    pred <- stats::predict(fit, newdata = sf::st_drop_geometry(te))
    xtr <- sf::st_coordinates(tr)[, 1:2]; xte <- sf::st_coordinates(te)[, 1:2]
    err <- vapply(as.character(p$levels), function(L) {
      cen <- attr(p, "centres")[[L]]
      a <- FNN::get.knnx(cen, xtr, 1L)$nn.index[, 1L]
      b <- FNN::get.knnx(cen, xte, 1L)$nn.index[, 1L]
      m <- tapply(res, factor(a, levels = seq_len(nrow(cen))), mean)
      m[is.na(m)] <- 0
      mean((te$z - pred - m[b])^2)
    }, numeric(1))
    regret[i] <- err[[as.character(sel$best)]] / min(err) - 1
  }
  expect_true(all(regret <= 0.15), info = paste(format(round(100 * regret, 1)), collapse = ", "))
})


test_that("a profile with an estimate fitted to a subsample hands on the centres it scored, and the summary checks them", {
  # The estimate is scored on every row of the layer in the cells of the
  # centres fitted to the subsample, so those are the cells its record
  # describes: nothing is kept to refit them, and they are handed on as they
  # are (the builder's first version refitted them, and the summary could
  # then no longer check the cells).
  pts <- .bp_layer(1200, 2103)
  for (est in c("pooled", "mean")) {
    info <- est
    p <- suppressWarnings(
      resolution_profile(pts, "z", sample_n = 300, n_levels = 6L, partition = "density",
                         nstart = 3L, sac = "auto", estimate = est, obs_per_cell = NULL))
    expect_identical(attr(p, "estimate")$value, est, info = info)
    expect_identical(attr(p, "bounds")$n_sample, 300L, info = info)
    expect_identical(attr(p, "estimate")$n, 1200L, info = info)
    expect_null(attr(p, "refit_rows"))
    lines <- capture_spatialkit_log(
      sel <- suppressWarnings(select_resolution(p, criterion = "loo")))
    # A profile with an estimate is not taken for one saved before the refit.
    expect_false(any(grepl("keeps no coordinates", lines, fixed = TRUE)), info = info)
    expect_null(attr(sel$seeds, "refit"))
    scored <- unname(attr(p, "centres")[[as.character(sel$best)]])
    expect_equal(unname(sf::st_coordinates(sel$seeds)[, 1:2]), scored, info = info)
    rec <- sel$count_rule
    # The record keeps what the profile scored, at the centres it scored,
    # and the seeds are those centres.
    expect_identical(rec$scored$loo, p$loo[match(sel$best, p$levels)], info = info)
    expect_equal(rec$scored$centres, scored, info = info)
    expect_identical(attr(sel$seeds, "count_rule")$scored_cells, TRUE, info = info)
    seeds <- suppressWarnings(get_voronoi_seeds(method = "kmeans", n = sel, sample_points = pts))
    expect_null(attr(seeds, "refit"))
    expect_identical(attr(seeds, "count_rule")$scored_cells, TRUE, info = info)
    tess <- build_tessellation(seeds, method = "voronoi", quiet = TRUE)
    asg <- assign_features_to_polygons(pts, tess$cells, polygon_id_col = "cell_id",
                                       quiet = TRUE)
    out <- .bp_catch(summarize_by_cell(asg, "z", id_col = "cell_id", cells_sf = tess$cells,
                                       estimate = est))
    expect_identical(attr(out$value, "count_rule")$check,
                     list(estimate = "match", cells = "match"), info = info)
    expect_length(.bp_of(out, "spatialkit_cells_differ"), 0L)
    expect_length(.bp_of(out, "spatialkit_estimate_mismatch"), 0L)
  }
})


# ---- The print header -----------------------------------------------------------------------

test_that("the print of a profile fitted to a subsample names the kind that was fitted", {
  pts <- .bp_layer(600, 2104)
  for (kind in c("density", "kmeans")) {
    p <- suppressWarnings(
      resolution_profile(pts, "z", predictor_vars = c("x1", "x2"), sac = .bp_sac(nugget = 0.3),
                         sample_n = 200, n_levels = 4L, partition = kind, nstart = 2L,
                         estimate = NULL, obs_per_cell = NULL))
    first <- utils::capture.output(print(p))[1L]
    expect_identical(first,
                     sprintf("Resolution profile: %d levels on 600 points (%s fitted to a subsample of 200)",
                             nrow(p),
                             if (kind == "density") "density-weighted k-means" else "k-means"))
  }
  # And with rows that carry no response.
  pts$z[1:50] <- NA
  p <- suppressWarnings(
    resolution_profile(pts, "z", predictor_vars = c("x1", "x2"), sac = .bp_sac(nugget = 0.3),
                       sample_n = 200, n_levels = 4L, partition = "density", nstart = 2L,
                       estimate = NULL, obs_per_cell = NULL))
  expect_match(utils::capture.output(print(p))[1L],
               "(550 with a response and every predictor; density-weighted k-means fitted to a subsample of 200)",
               fixed = TRUE)
})


# ---- The help -----------------------------------------------------------------------------

test_that("the help lists the two classes, the source of Cp's noise variance and the refit", {
  # A help page's text, markup dropped, on one line: the file in man/ from the
  # source tree, the installed page otherwise.
  flat <- function(x)
    paste(unlist(lapply(x, function(el) if (is.character(el)) el else flat(el))), collapse = "")
  rd <- function(page) {
    root <- testthat::test_path("..", "..")
    file <- file.path(root, "man", paste0(page, ".Rd"))
    db <- if (file.exists(file)) tools::parse_Rd(file)
          else tools::Rd_db("spatialkit")[[paste0(page, ".Rd")]]
    gsub("[[:space:]]+", " ", flat(db))
  }
  rp <- rd("resolution_profile")
  expect_match(rp, "spatialkit_variogram_unusable", fixed = TRUE)
  expect_match(rp, "spatialkit_cp_noise", fixed = TRUE)
  expect_match(rp, "The remaining warnings of these functions (points dropped, a variogram of another variable) are R warnings without a class",
               fixed = TRUE)
  expect_false(grepl("a zero nugget, a variogram of another variable", rp, fixed = TRUE))
  expect_match(rp, "is the fit's lower bound and not an estimate", fixed = TRUE)
  expect_match(rp, "held to the semivariance of nearest-neighbour pairs", fixed = TRUE)
  expect_match(rp, "\"short-lag semivariance\"", fixed = TRUE)
  expect_match(rp, "refit_rows", fixed = TRUE)
  sr <- rd("select_resolution")
  expect_match(sr, "refitted to every point", fixed = TRUE)
  expect_match(sr, "attr(, \"refit\")", fixed = TRUE)
  expect_match(sr, "selection off a profile without an estimate has count_rule as well",
               fixed = TRUE)
  # Which profiles are refitted, and which are not.
  expect_match(sr, "On a profile without an estimate whose centres were fitted to a subsample",
               fixed = TRUE)
  expect_match(sr, "which every such profile this version builds carries; on one built before it the centres scored are the seeds, and a log line says so",
               fixed = TRUE)
  expect_match(sr, "A profile with an estimate scores it on every row of the layer in the cells of the centres fitted to the subsample, so its seeds are those centres",
               fixed = TRUE)
  expect_match(sr, "Lloyd's iteration with every weight one runs in its place, and a log line says so",
               fixed = TRUE)
  expect_match(rp, "a profile with an estimate scores it on every row in the cells of the centres fitted, and hands those centres on",
               fixed = TRUE)
  expect_match(rp, "when the centres were fitted to a subsample and the profile has no estimate, refit_rows",
               fixed = TRUE)
  expect_match(rp, "gives cp its nugget (held to the semivariance of nearest-neighbour pairs, and not at all when it is at the fit's lower bound), with a warning",
               fixed = TRUE)
  expect_match(rp, "so on a layer of more than sample_n points the count read off the profile can move with seed and with sample_n",
               fixed = TRUE)
})

test_that("NEWS and ?spatialkit say what the ways back to the earlier results no longer give", {
  root <- testthat::test_path("..", "..")
  # ?spatialkit: the paragraph on the two choices.
  file <- file.path(root, "man", "spatialkit-package.Rd")
  flat <- function(x)
    paste(unlist(lapply(x, function(el) if (is.character(el)) el else flat(el))), collapse = "")
  db <- if (file.exists(file)) tools::parse_Rd(file)
        else tools::Rd_db("spatialkit")[["spatialkit-package.Rd"]]
  pk <- gsub("[[:space:]]+", " ", flat(db))
  # Whether `text` says `x`, failing with `x` alone and not the whole text.
  says <- function(text, x) expect_true(grepl(x, text, fixed = TRUE), info = x)
  never <- function(text, x) expect_false(grepl(x, text, fixed = TRUE), info = x)
  says(pk, "the profile it returned before, with two exceptions.")
  says(pk, "so cp and the count it chooses differ wherever either applies, and nothing gives the earlier noise variance back.")
  says(pk, "is refitted to every point before its seeds are handed on; the centres the profile scored are attr(profile, \"centres\")[[as.character(sel$best)]].")
  never(pk, "the profile it returned before. The seeding functions")

  news_file <- file.path(root, "NEWS.md")
  if (!file.exists(news_file)) news_file <- system.file("NEWS.md", package = "spatialkit")
  if (!nzchar(news_file) || !file.exists(news_file)) skip("NEWS.md is not there to read")
  lines <- readLines(news_file, warn = FALSE, encoding = "UTF-8")
  top <- grep("^# ", lines)
  dev <- lines[seq.int(top[1L], top[2L] - 1L)]
  breaking <- dev[seq.int(grep("^## Breaking changes", dev), grep("^## New features", dev) - 1L)]
  bc <- gsub("[[:space:]]+", " ", paste(breaking, collapse = " "))
  # The promises that no longer hold are gone.
  never(bc, "On a profile without an estimate nothing changes.")
  never(bc, "returns the centres it scored, as before. ")
  never(bc, "For the count Cp gave before, build the profile")
  # Each of the three entries names the exceptions, and two entries say
  # what changed and what, if anything, gives the earlier result.
  says(bc, "gives Cp's count, with the noise variance of the entry on Cp below.")
  says(bc, "obs_per_cell = NULL` gives the previous ladder and columns, except `cp` wherever the entry on Cp's noise variance below applies")
  says(bc, "returns the centres it scored, as before, except off a profile without an estimate fitted to a subsample, whose chosen level is refitted to every point")
  says(bc, "Cp's noise variance is then the variogram's nugget, as in 2.0.0, wherever that nugget is above the fit's lower bound and no more than the semivariance of nearest-neighbour pairs.")
  says(bc, "and no argument gives the earlier noise variance back.")
  says(bc, "are not the centres 2.0.0 returned, `partition = \"kmeans\"` and `estimate = NULL, obs_per_cell = NULL` included. The centres the profile scored, which 2.0.0 returned, are `attr(profile, \"centres\")[[as.character(sel$best)]]`")
})
