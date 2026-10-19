# tests/testthat/test-core-defaults.R
# ---------------------------------------------------------------------------
# The defaults that were switched together, at their values after the
# switch: eleven formal defaults in five functions, and the value of one
# internal function.  A call that names none of them has pooled values at
# 4.5 rows per cell, in cells of the density kind fitted from 3 k-means++
# starts, with no variogram fitted for a ladder of sizes.  The values they
# had before are the ways back, and test-core-switch.R holds each way back
# to the object it returned.
# ---------------------------------------------------------------------------

test_that("the eleven defaults that switch together hold their switched values", {
  profile <- formals(resolution_profile)
  seeds   <- formals(get_voronoi_seeds)
  kmeans  <- formals(voronoi_seeds_kmeans)
  select  <- formals(select_resolution)
  values  <- formals(summarize_by_cell)

  # 1 to 3: the estimate follows the response, the ladder follows the
  # estimate, and a variogram is fitted on a present ladder only.
  expect_true("estimate" %in% names(profile))
  expect_identical(profile$estimate, "auto")
  expect_true("obs_per_cell" %in% names(profile))
  expect_identical(profile$obs_per_cell, "auto")
  expect_true("sac" %in% names(profile))
  expect_identical(profile$sac, "auto")

  # 4: a profile is read by the rule that goes with its estimate.  "auto"
  # moved to the front and nothing else in the vector changed.
  expect_identical(eval(select$criterion),
                   c("auto", "cp", "reliability", "elbow", "moran_z", "loo", "size"))
  expect_identical(eval(select$criterion)[1L], "auto")

  # 5: the summary reports pooled values beside the plain mean.
  expect_identical(eval(values$estimate), c("pooled", "mean"))
  expect_identical(eval(values$estimate)[1L], "pooled")

  # 6 to 8: the cells are of the density kind in all three functions.
  for (f in list(profile, seeds, kmeans)) {
    expect_identical(eval(f$partition), c("density", "kmeans"))
    expect_identical(eval(f$partition)[1L], "density")
  }

  # 9 to 11: the number of starts is the kind's own.
  expect_true("nstart" %in% names(profile))
  expect_null(profile$nstart)
  expect_true("kmeans_nstart" %in% names(seeds))
  expect_null(seeds$kmeans_nstart)
  expect_true("nstart" %in% names(kmeans))
  expect_null(kmeans$nstart)
})

test_that("Moran's z is left out without its packages on every ladder", {
  expect_identical(.moran_skip_any_ladder(), TRUE)
  expect_identical(formals(.moran_skip_any_ladder), NULL)
})

test_that("the formals around the switched ones are as they were", {
  # The switch edits the eleven defaults and nothing else in a signature.
  expect_identical(
    names(formals(resolution_profile)),
    c("data_sf", "response_var", "predictor_vars", "levels", "n_levels", "min_cell_n",
      "sample_n", "nstart", "seed", "sac", "select_on", "range_floor", "partition",
      "init", "density_power", "density_k", "iter_max", "estimate", "obs_per_cell"))
  expect_identical(
    names(formals(get_voronoi_seeds)),
    c("boundary", "method", "n", "seeds", "sample_points", "kmeans_nstart", "kmeans_iter",
      "set_seed", "partition", "init", "density_power", "density_k"))
  expect_identical(
    names(formals(voronoi_seeds_kmeans)),
    c("points_sf", "k", "set_seed", "nstart", "partition", "init", "density_power",
      "density_k", "iter_max"))
  expect_identical(
    names(formals(summarize_by_cell)),
    c("assigned_points_sf", "response_var", "predictor_vars", "id_col", "agg_funs",
      "cells_sf", "deff", "sac", "deff_max_n", "quiet", "conf_level", "area", "estimate",
      "obs_level"))
  # Not switched (section 6): the arguments an interval, the design effect,
  # the present ladder and the density kind read keep their values.
  p <- formals(resolution_profile)
  expect_identical(p$n_levels, 20L)
  expect_identical(p$min_cell_n, 9L)
  expect_identical(p$sample_n, 1500L)
  expect_identical(p$seed, 123L)
  expect_identical(p$range_floor, FALSE)
  expect_identical(eval(p$init), c("kmeans++", "random"))
  expect_identical(p$density_power, 1)
  expect_identical(p$density_k, 8L)
  expect_null(p$iter_max)
  v <- formals(summarize_by_cell)
  expect_null(v$conf_level)
  expect_null(v$obs_level)
  expect_identical(v$deff, 1)
  expect_null(v$sac)
  expect_identical(formals(get_voronoi_seeds)$kmeans_iter, 100)
  expect_null(formals(get_voronoi_seeds)$init)
  expect_null(formals(voronoi_seeds_kmeans)$init)
  expect_identical(formals(voronoi_seeds_kmeans)$set_seed, 456)
  expect_identical(eval(formals(get_voronoi_seeds)$method),
                   c("kmeans", "random", "provided"))
})

test_that("the values the defaults switched to are values their arguments take", {
  # Each value one of the eleven switched to is a value its argument took
  # before the switch, so the switch edited defaults and nothing else.
  expect_true("auto" %in% eval(formals(select_resolution)$criterion))
  expect_true("pooled" %in% eval(formals(summarize_by_cell)$estimate))
  for (f in list(resolution_profile, get_voronoi_seeds, voronoi_seeds_kmeans))
    expect_true("density" %in% eval(formals(f)$partition))
  expect_identical(.NSTART_DENSITY, 3L)
  expect_identical(.OBS_PER_CELL_DEFAULT, 4.5)
  # The arguments this stage's rules read sit after the present ones.
  expect_identical(names(formals(select_resolution)),
                   c("profile", "criterion", "tol", "range_floor", "diagnostic", "obs_per_cell"))
  expect_null(formals(select_resolution)$obs_per_cell)
  expect_identical(formals(select_resolution)$tol, 0.02)
  expect_identical(formals(select_resolution)$range_floor, TRUE)
  expect_identical(formals(select_resolution)$diagnostic, FALSE)
})
