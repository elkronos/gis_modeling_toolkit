# tests/testthat/test-core-switch.R
# ---------------------------------------------------------------------------
# The defaults after the switch, read off the functions themselves.
#
# The default path end to end, with nothing named beyond the data and the
# response: a profile of one level at round(n / 4.5) cells of the density
# kind, with `cp` and `reliability` NA and no variogram fitted; a selection
# by size; seeds, cells and the assignment; a summary with pooled columns
# and a record whose check reads "match".
#
# Each way back gives the object the baseline recorded before the switch
# (fixtures/core-frozen.rds, written by dev/baseline-core-frozen.R): the
# profile with estimate = NULL, obs_per_cell = NULL, sac = NULL and
# partition = "kmeans"; criterion = "cp" on it; the seeding functions with
# partition = "kmeans"; the summary with estimate = "mean".  The same calls
# as they were recorded, with nothing named, now give other objects.
#
# What the seeding functions say of the kind when a call names none: nothing
# beside scored centres of either kind and any density setting, and nothing
# beside method = "random" or "provided".  And Moran's z on a present ladder
# above its limit without FNN: NA with a warning, where it stopped before.
# ---------------------------------------------------------------------------

.sw_seed <- function(seed, env = parent.frame())
  withr::local_seed(seed, .local_envir = env, .rng_kind = "Mersenne-Twister",
                    .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")

# A call with its warnings and messages kept and muffled.
.sw_catch <- function(expr) {
  warnings <- list()
  messages <- character(0)
  value <- withCallingHandlers(expr,
    warning = function(w) {
      warnings[[length(warnings) + 1L]] <<- w
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  list(value = value, warnings = warnings, messages = messages)
}
.sw_of <- function(r, class) Filter(function(w) inherits(w, class), r$warnings)
.sw_text <- function(r) vapply(r$warnings, conditionMessage, character(1))

# The layer of the end-to-end tests of the record: seed 7001; `u`, then `v`,
# 3,000 values each, uniform on 0 to 1e5; then the response,
# sin(u / 2e4) + cos(v / 1.5e4) plus normal noise of standard deviation 1.5;
# the points at (5e5 + u, 5e6 + v) in EPSG:32617.
.sw_store <- new.env(parent = emptyenv())
.sw_layer <- function() {
  if (is.null(.sw_store$pts)) {
    d <- local({
      .sw_seed(7001)
      u <- stats::runif(3000, 0, 1e5)
      v <- stats::runif(3000, 0, 1e5)
      z <- sin(u / 2e4) + cos(v / 1.5e4) + stats::rnorm(3000, 0, 1.5)
      data.frame(x = 5e5 + u, y = 5e6 + v, z = z)
    })
    .sw_store$pts <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32617)
  }
  .sw_store$pts
}


# ---- The default path --------------------------------------------------------------

test_that("the default path is pooled values at 4.5 rows a cell, in density cells, with no variogram", {
  pts <- .sw_layer()
  # The profile, with nothing named but the data and the response.
  log <- NULL
  P <- NULL
  log <- capture_spatialkit_log(P <- .sw_catch(resolution_profile(pts, "z")),
                                level = logger::WARN)
  prof <- P$value
  expect_length(P$warnings, 0L)
  expect_length(P$messages, 0L)
  expect_false(any(grepl("variogram", log, ignore.case = TRUE)))
  # One level, at round(n / 4.5) cells of the density kind from 3 k-means++
  # starts, scored for the pooled value.
  expect_identical(nrow(prof), 1L)
  expect_identical(prof$levels, 667L)
  expect_identical(prof$levels, as.integer(round(3000 / 4.5)))
  expect_identical(attr(prof, "estimate")$value, "pooled")
  expect_identical(attr(prof, "partition")[c("kind", "init", "density_power", "density_k")],
                   list(kind = "density", init = "kmeans++", density_power = 1,
                        density_k = 8L))
  expect_identical(attr(prof, "nstart"), 3L)
  expect_true(is.finite(prof$loo))
  # No variogram: `cp` and `reliability` are NA, and the profile says why.
  expect_identical(prof$cp, NA_real_)
  expect_identical(prof$reliability, NA_real_)
  expect_null(attr(prof, "variogram"))
  expect_null(attr(prof, "sac"))
  expect_identical(attr(prof, "bounds")$range, NA_real_)
  expect_identical(attr(prof, "cp_noise")$source, "none")
  # The profile is the one the values the defaults resolve to give when they
  # are named.
  expect_identical(
    resolution_profile(pts, "z", estimate = "pooled", obs_per_cell = 4.5, sac = "auto",
                       partition = "density", nstart = 3),
    prof)

  # The selection: by the size, with no rule named.
  S <- .sw_catch(select_resolution(prof))
  sel <- S$value
  expect_length(S$warnings, 0L)
  expect_identical(sel$criterion, "size")
  expect_identical(sel$best, 667L)
  expect_identical(sel$count_rule$rule, "size")
  expect_identical(sel$count_rule$estimate, "pooled")
  expect_identical(sel$count_rule$partition$kind, "density")
  expect_identical(sel$count_rule$partition$nstart, 3L)
  expect_identical(sel, select_resolution(prof, criterion = "size"))

  # The seeds, the cells and the assignment, by the calls a user writes.
  G <- .sw_catch(get_voronoi_seeds(method = "kmeans", n = sel, sample_points = pts))
  expect_length(G$warnings, 0L)
  seeds <- G$value
  expect_identical(nrow(seeds), 667L)
  expect_identical(attr(seeds, "partition")$kind, "density")
  expect_identical(attr(seeds, "count_rule")$scored_cells, TRUE)
  tess <- build_tessellation(seeds, method = "voronoi", quiet = TRUE)
  asg  <- assign_features_to_polygons(pts, tess$cells, polygon_id_col = "cell_id", quiet = TRUE)
  expect_identical(nrow(asg), 3000L)
  # Every row in the cell of its nearest scored centre.
  rec <- sel$count_rule
  cen <- unname(sf::st_coordinates(sel$seeds)[, 1:2])
  xy  <- unname(sf::st_coordinates(sf::st_transform(pts, rec$crs))[, 1:2])
  nearest <- .nearest_centre(xy[, 1L], xy[, 2L], cen[, 1L], cen[, 2L], .has_backend("FNN"))
  expect_identical(asg$cell_id, tess$index[nearest])

  # The summary, with no estimate named: pooled columns beside the plain
  # ones, and a record whose check reads "match".
  V <- .sw_catch(summarize_by_cell(asg, "z", id_col = "cell_id", cells_sf = tess$cells))
  out <- V$value
  expect_length(V$warnings, 0L)
  expect_true(all(c("resp_mean_z", "..pooled_resp_z", "..pooled_w_resp_z",
                    "..pooled_target_resp_z") %in% names(out)))
  expect_false(any(c("..obs_lo_resp_z", "..obs_hi_resp_z", "..ci_lo_resp_z") %in% names(out)))
  expect_false(is.null(attr(out, "pooled")))
  got <- attr(out, "count_rule")
  expect_identical(got$check, list(estimate = "match", cells = "match"))
  expect_identical(got$rule, "size")
  # The summary's leave-one-out error is the profile's at that level.
  own <- attr(out, "pooled")$columns$resp_z
  expect_lt(abs(own$loo - prof$loo) / prof$loo, 1e-10)
  expect_identical(out, suppressWarnings(
    summarize_by_cell(asg, "z", id_col = "cell_id", cells_sf = tess$cells,
                      estimate = "pooled")))

  # voronoi_seeds_kmeans() given the selection returns the same centres.
  K <- .sw_catch(voronoi_seeds_kmeans(pts, sel))
  expect_length(K$warnings, 0L)
  expect_equal(partition_xy(K$value), partition_xy(seeds), tolerance = 1e-12)
})

test_that("a count given as a number is density cells from 3 k-means++ starts", {
  pts <- .sw_layer()[1:600, ]
  a <- get_voronoi_seeds(method = "kmeans", n = 40, sample_points = pts, set_seed = 5)
  expect_identical(a, get_voronoi_seeds(method = "kmeans", n = 40, sample_points = pts,
                                        set_seed = 5, partition = "density", kmeans_nstart = 3))
  expect_identical(attr(a, "kmeans")$nstart, 3L)
  b <- voronoi_seeds_kmeans(pts, 40, set_seed = 5)
  expect_identical(b, voronoi_seeds_kmeans(pts, 40, set_seed = 5, partition = "density",
                                           nstart = 3))
  expect_identical(attr(b, "partition")$kind, "density")
  # A number of starts that is not whole is refused under the density kind,
  # which a call that names no kind now has.
  expect_error(get_voronoi_seeds(method = "kmeans", n = 40, sample_points = pts,
                                 kmeans_nstart = 2.5),
               "`kmeans_nstart` must be NULL or a single positive whole number of starts; got 2.5.",
               fixed = TRUE)
  expect_s3_class(suppressWarnings(
    get_voronoi_seeds(method = "kmeans", n = 40, sample_points = pts, kmeans_nstart = 2.5,
                      partition = "kmeans")), "sf")
})


# ---- The ways back -------------------------------------------------------------------

.sw_fixture <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) cache <<- readRDS(test_path("fixtures", "core-frozen.rds"))
    cache
  }
})

# A record against the baseline's: by the rule for any platform, and on the
# baseline's own platform identical in every part.
.sw_expect_baseline <- function(got, want, nm) {
  d <- core_frozen_differences(got, want)
  expect(length(d) == 0L,
         paste(c(sprintf("`%s` is not what it was on the baseline.", nm), d), collapse = "\n"))
  if (core_frozen_same_platform(core_frozen_platform(), .sw_fixture()$info))
    for (part in c("value", "conditions", "output", "printed", "log", "rng"))
      expect_identical(got[[part]], want[[part]], label = sprintf("the %s of `%s`", part, nm),
                       expected.label = "the baseline's")
}

# Calls of the fixture, in the order the list runs them, by name; `old`
# names those to run as the baseline recorded them.
.sw_calls <- function(names, old = character(0)) {
  base <- .sw_fixture()
  calls <- core_frozen_calls()
  calls <- calls[intersect(names(calls), names)]
  for (nm in old) calls[[nm]] <- str2lang(base$calls_recorded[[nm]])
  calls
}

test_that("each way back gives the object the baseline recorded", {
  skip_on_cran()
  skip_if_not_installed("gstat")
  base <- .sw_fixture()
  ways <- c("profile_default", "select_default_cp", "select_default_unnamed",
            "seeds_number", "seeds_selection", "kmeans_number", "kmeans_selection",
            "cells_number", "assigned_number", "values_default", "values_cells_sf")
  # Each way back names the old default.
  text <- core_frozen_call_text(.sw_calls(ways))
  expect_match(text$profile_default, paste0("estimate = NULL, obs_per_cell = NULL, ",
                                            "sac = NULL, partition = \"kmeans\""), fixed = TRUE)
  expect_match(text$select_default_unnamed, "criterion = \"cp\"", fixed = TRUE)
  for (nm in c("seeds_number", "seeds_selection", "kmeans_number", "kmeans_selection"))
    expect_match(text[[nm]], "partition = \"kmeans\"", fixed = TRUE, info = nm)
  for (nm in c("values_default", "values_cells_sf"))
    expect_match(text[[nm]], "estimate = \"mean\"", fixed = TRUE, info = nm)

  now <- core_frozen_run(base$data, .sw_calls(ways))
  for (nm in ways) .sw_expect_baseline(now[[nm]], base$records[[nm]], nm)
})

test_that("the calls as the baseline recorded them, with nothing named, give other objects", {
  skip_on_cran()
  skip_if_not_installed("gstat")
  base <- .sw_fixture()
  moved <- c("profile_default", "seeds_number", "kmeans_number", "values_default")
  run <- core_frozen_run(base$data, .sw_calls(c(moved, "cells_number", "assigned_number"),
                                              old = moved))
  for (nm in moved)
    expect_false(isTRUE(all.equal(run[[nm]]$value, base$records[[nm]]$value)), info = nm)
  # The profile is one level of pooled values in density cells.
  p <- run$profile_default$value
  expect_identical(nrow(p), 1L)
  expect_identical(attr(p, "estimate")$value, "pooled")
  expect_identical(attr(p, "partition")$kind, "density")
  # The seeds are density seeds; the summary adds the pooled columns to the
  # plain ones it had, which are unchanged.
  expect_identical(attr(run$seeds_number$value, "partition")$kind, "density")
  expect_identical(attr(run$kmeans_number$value, "partition")$kind, "density")
  v <- run$values_default$value
  old <- base$records$values_default$value
  expect_true("..pooled_resp_z" %in% names(v))
  expect_true(all(names(old) %in% names(v)))
})

test_that("calls that read no default the switch moved give the baseline's object with nothing named", {
  skip_on_cran()
  skip_if_not_installed("gstat")
  base <- .sw_fixture()
  # Scored centres of the density kind passed as the count, and the random
  # and provided seedings: their text names no kind, as when recorded.
  same <- c("seeds_density_selection", "seeds_density_profile", "kmeans_density_selection",
            "kmeans_density_profile")
  expect_identical(core_frozen_call_text(.sw_calls(same)), base$calls_recorded[same])
  random <- c("seeds_random", "seeds_provided", "seeds_random_lonlat")
  calls <- .sw_calls(c("profile_density", "select_density_cp", "profile_density_clustered",
                       same, random), old = random)
  now <- core_frozen_run(base$data, calls)
  for (nm in c(same, random)) .sw_expect_baseline(now[[nm]], base$records[[nm]], nm)
  for (nm in random) expect_length(now[[nm]]$conditions, 0L)
})


# ---- What the seeding functions say of the kind ----------------------------------------

.sw_pts <- function()
  partition_points_sf(partition_clustered_xy(n_bg = 100L, n_cl = 300L, seed = 29L))

# Four profiles of scored centres at two levels from two starts, every
# setting the switch moved named: the k-means kind, and the density kind at
# its defaults, at power 0.5 and with 12 neighbours.
.sw_profiles <- function() {
  if (is.null(.sw_store$profiles)) {
    pts <- .sw_pts()
    prof <- function(...)
      resolution_profile(pts, levels = c(6L, 14L), nstart = 2, estimate = NULL,
                         obs_per_cell = NULL, sac = NULL, ...)
    .sw_store$profiles <- list(
      kmeans  = prof(partition = "kmeans"),
      density = prof(partition = "density"),
      power   = prof(partition = "density", density_power = 0.5),
      k_nn    = prof(partition = "density", density_k = 12))
  }
  .sw_store$profiles
}
# A selection at 14 cells that carries the level's centres and nothing a
# warning about the count could be raised for.
.sw_sel <- function(p, L = 14L) structure(
  list(best = L, seeds = spatialkit:::.profile_seeds(p, L), criterion = "cp",
       edge = NA_character_, structured = TRUE), class = "resolution_selection")

test_that("a selection passed as the count by a call that names no partition argument raises nothing", {
  pts <- .sw_pts()
  profiles <- .sw_profiles()
  via <- list(function(sel, ...) voronoi_seeds_kmeans(pts, sel, ...),
              function(sel, ...) get_voronoi_seeds(method = "kmeans", n = sel,
                                                   sample_points = pts, ...))
  unused <- function(r) length(.sw_of(r, "spatialkit_partition_unused"))
  for (nm in names(profiles)) {
    sel <- .sw_sel(profiles[[nm]])
    for (f in via) {
      r <- .sw_catch(f(sel))
      expect_length(r$warnings, 0L)
      expect_equal(partition_xy(r$value), attr(profiles[[nm]], "centres")[["14"]],
                   tolerance = 1e-12)
      expect_identical(attr(r$value, "partition"),
                       attr(profiles[[nm]], "partition")[c("kind", "init", "density_power",
                                                           "density_k")])
    }
  }
  # Beside the profile built at power 0.5: the power it has, named, raises
  # nothing; another power or neighbour count is said to be unused.
  sel <- .sw_sel(profiles$power)
  for (f in via) {
    expect_length(.sw_catch(f(sel, density_power = 0.5))$warnings, 0L)
    for (other in list(list(density_power = 0.7), list(density_k = 4))) {
      r <- .sw_catch(do.call(f, c(list(sel), other)))
      expect_identical(unused(r), 1L)
      expect_length(r$warnings, 1L)
    }
  }
  # Beside the profile built with 12 neighbours.
  sel <- .sw_sel(profiles$k_nn)
  for (f in via) {
    expect_length(.sw_catch(f(sel, density_k = 12))$warnings, 0L)
    expect_identical(unused(.sw_catch(f(sel, density_k = 8))), 1L)
  }
  # Beside the k-means profile: a density setting supplied is read by
  # nothing, and the k-means kind named is the kind that built it.
  sel <- .sw_sel(profiles$kmeans)
  for (f in via) {
    r <- .sw_catch(f(sel, density_power = 0.5))
    expect_identical(unused(r), 1L)
    expect_length(r$warnings, 1L)
    expect_length(.sw_catch(f(sel, partition = "kmeans"))$warnings, 0L)
  }
})

test_that("a selection made by select_resolution() is handed on in silence whatever its settings", {
  pts <- .sw_pts()
  xy <- partition_xy(pts)
  pts$z <- local({
    .sw_seed(30)
    sin(xy[, 1L] / 150) + stats::rnorm(nrow(xy), sd = 0.5)
  })
  for (args in list(list(density_power = 0.5), list(density_k = 12))) {
    prof <- suppressWarnings(do.call(resolution_profile,
                    c(list(pts, "z", levels = c(6L, 14L), nstart = 2, estimate = NULL,
                           obs_per_cell = NULL, sac = NULL, partition = "density"), args)))
    expect_identical(attr(prof, "partition")[[names(args)]],
                     if (identical(names(args), "density_k")) 12L else 0.5)
    sel <- select_resolution(prof, criterion = "cp")
    for (r in list(.sw_catch(voronoi_seeds_kmeans(pts, sel)),
                   .sw_catch(get_voronoi_seeds(method = "kmeans", n = sel, sample_points = pts)))) {
      # Only what the count itself may be told: an optimum at an end of the
      # levels, or a response without structure.
      classes <- vapply(r$warnings, function(w) class(w)[1L], character(1))
      expect_true(all(classes %in% c("spatialkit_bound_chosen", "spatialkit_no_structure")))
      expect_length(.sw_of(r, "spatialkit_partition_unused"), 0L)
      expect_equal(partition_xy(r$value), attr(prof, "centres")[[as.character(sel$best)]],
                   tolerance = 1e-12)
      expect_identical(attr(r$value, "partition")[[names(args)]],
                       attr(prof, "partition")[[names(args)]])
    }
  }
})

test_that("a random or provided seeding that names no kind raises nothing", {
  pts <- .sw_pts()
  bnd <- sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(pts)))
  random   <- function(...) get_voronoi_seeds(bnd, method = "random", n = 5, set_seed = 2, ...)
  provided <- function(...) get_voronoi_seeds(method = "provided", seeds = pts[1:4, ], ...)
  for (case in list(list(random, "random"), list(provided, "provided"))) {
    f <- case[[1L]]; method <- case[[2L]]
    base <- .sw_catch(f())
    expect_length(base$warnings, 0L)
    # A density setting supplied is read by nothing, and the call is told.
    r <- .sw_catch(f(density_power = 0.5))
    expect_identical(.sw_text(r), sprintf(paste0(
      "get_voronoi_seeds(): `density_power` is read by method = \"kmeans\" only; method = ",
      "\"%s\" ignores it."), method))
    expect_identical(r$value, base$value)
    # The k-means kind named is nothing to say.
    expect_length(.sw_catch(f(partition = "kmeans"))$warnings, 0L)
    expect_identical(f(partition = "kmeans"), base$value)
  }
})


# ---- Moran's z above its limit on a present ladder -----------------------------------

test_that("a present ladder above the limit of Moran's z returns without FNN or Matrix", {
  skip_if_not_installed("gstat")
  pts <- local({
    d <- local({
      .sw_seed(323)
      x <- stats::runif(720, 0, 1e5); y <- stats::runif(720, 0, 1e5)
      data.frame(x = x, y = y, z = sin(x / 2e4) + cos(y / 1.5e4) + stats::rnorm(720, 0, 1.5))
    })
    sf::st_as_sf(d, coords = c("x", "y"), crs = 32617)
  })
  for (pkg in c("FNN", "Matrix")) {
    r <- local({
      testthat::local_mocked_bindings(.KNN_DENSE_MAX = 159L)
      testthat::local_mocked_bindings(
        .has_backend = function(p) if (p %in% pkg) FALSE else requireNamespace(p, quietly = TRUE))
      .sw_catch(resolution_profile(pts, "z", levels = 160, sample_n = 720, estimate = NULL,
                                   obs_per_cell = NULL, sac = NULL, partition = "kmeans",
                                   init = "random", nstart = 1))
    })
    expect_s3_class(r$value, "resolution_profile")
    w <- .sw_of(r, "spatialkit_moran_skipped")
    expect_length(w, 1L)
    expect_identical(w[[1L]]$missing, pkg)
    expect_identical(r$value$levels, 160L)
    expect_identical(r$value$moran_i, NA_real_, info = pkg)
    expect_identical(r$value$moran_z, NA_real_, info = pkg)
    expect_true(is.finite(r$value$rss))
  }
})
