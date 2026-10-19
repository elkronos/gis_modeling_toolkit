# tests/testthat/test-core-partition-convergence.R
# ---------------------------------------------------------------------------
# The iteration limit.  A start of the density kind runs for at most 100
# steps of Lloyd's iteration in the profile and the seeding functions (the
# `iter_max` and `kmeans_iter` arguments set another limit).  A start that
# was kept and had not converged at the limit is returned with the centres
# of its last step, and that is said: a warning of class
# "spatialkit_not_converged", and `converged` in the record.
# ---------------------------------------------------------------------------

.pv_engine <- function(xy, k, seed, iter.max, init = "kmeans++", nstart = 1L) {
  set.seed(seed)
  spatialkit:::.partition_centres(xy, k, partition = "density", init = init,
                                  nstart = nstart, iter.max = iter.max)
}

# Two single starts on 1,200 uniform points that take more than 40 steps to
# stop moving and fewer than 100 (measured: 59 and 48), one from a random
# start and one from a k-means++ start.
.pv_cases <- list(
  random   = list(k = 20L, seed = 3L,  init = "random"),
  kmeanspp = list(k = 10L, seed = 20L, init = "kmeans++"))


test_that("the engine says whether the start it kept stopped before the limit", {
  xy <- partition_clustered_xy(n_bg = 150L, n_cl = 350L, seed = 7L)
  full <- .pv_engine(xy, 30L, 4L, 500L)
  expect_true(full$converged)
  expect_gt(full$iter, 3L)
  # Stopped at exactly the step it needs: converged.  One step fewer: not.
  at <- .pv_engine(xy, 30L, 4L, full$iter)
  expect_true(at$converged)
  expect_identical(at$centers, full$centers)
  short <- .pv_engine(xy, 30L, 4L, full$iter - 1L)
  expect_false(short$converged)
  expect_identical(short$iter, full$iter - 1L)
  cut <- .pv_engine(xy, 30L, 4L, 3L)
  expect_false(cut$converged)
  expect_identical(cut$iter, 3L)
  # The cut-off centres are those of the last step, with their own cells: a
  # partition, of a higher weighted sum than the converged one.
  expect_identical(cut$cluster, as.integer(partition_nearest(xy, cut$centers)))
  expect_gt(cut$objective, full$objective)
  # The k-means kind: stats::kmeans()'s own verdict, with the iterations held
  # to the limit (it reports one more for a fit it cut off).
  set.seed(4)
  km <- spatialkit:::.partition_centres(partition_sorted(xy), 30L, partition = "kmeans",
                                        init = "kmeans++", nstart = 2L, iter.max = 1L)
  expect_false(km$converged)
  expect_identical(km$iter, 1L)
  expect_identical(km$km$ifault, 2L)
  set.seed(4)
  km <- spatialkit:::.partition_centres(partition_sorted(xy), 30L, partition = "kmeans",
                                        init = "kmeans++", nstart = 2L)
  expect_true(km$converged)
  expect_identical(km$km$ifault, 0L)
})

test_that("the two test cases take more than 40 steps and fewer than 100", {
  xy <- partition_uniform_xy(1200L, seed = 41L)
  for (nm in names(.pv_cases)) {
    cs <- .pv_cases[[nm]]
    full <- .pv_engine(xy, cs$k, cs$seed, 100L, cs$init)
    expect_true(full$converged, label = nm)
    expect_gt(full$iter, 45L, label = nm)
    expect_lt(full$iter, 100L, label = nm)
    # Cut off at 40 or at 10 steps the centres are elsewhere, by tens of
    # units on a square of side 1000, and points are in other cells.
    for (lim in c(40L, 10L)) {
      cut <- .pv_engine(xy, cs$k, cs$seed, lim, cs$init)
      expect_false(cut$converged, label = nm)
      expect_gt(max(abs(cut$centers - full$centers)), 10, label = nm)
      expect_gt(sum(cut$cluster != full$cluster), 20L, label = nm)
    }
  }
})

test_that("voronoi_seeds_kmeans() runs a density start for up to 100 steps, or iter_max", {
  xy  <- partition_uniform_xy(1200L, seed = 41L)
  pts <- partition_points_sf(xy)
  for (nm in names(.pv_cases)) {
    cs <- .pv_cases[[nm]]
    full <- .pv_engine(xy, cs$k, cs$seed, 100L, cs$init)
    r <- partition_warnings(voronoi_seeds_kmeans(pts, cs$k, partition = "density", init = cs$init,
                                                 nstart = 1, set_seed = cs$seed))
    expect_length(r$warnings, 0L)
    expect_identical(partition_xy(r$value), full$centers, label = nm)
    expect_identical(attr(r$value, "partition")[c("iter", "iter_max", "converged")],
                     list(iter = full$iter, iter_max = 100L, converged = TRUE), label = nm)
    # A limit of 40: the engine's centres at 40 steps, and a warning.
    cut <- .pv_engine(xy, cs$k, cs$seed, 40L, cs$init)
    r <- partition_warnings(voronoi_seeds_kmeans(pts, cs$k, partition = "density", init = cs$init,
                                                 nstart = 1, set_seed = cs$seed, iter_max = 40))
    expect_identical(partition_xy(r$value), cut$centers, label = nm)
    expect_length(r$warnings, 1L)
    w <- r$warnings[[1L]]
    expect_s3_class(w, "spatialkit_not_converged")
    expect_identical(w[c("level", "iter", "iter_max")],
                     list(level = cs$k, iter = 40L, iter_max = 40L))
    expect_match(conditionMessage(w), sprintf(
      "voronoi_seeds_kmeans(): at %d cells the start that was kept had not converged after 40 steps of Lloyd's iteration (partition = \"density\")",
      cs$k), fixed = TRUE)
    expect_match(conditionMessage(w), "Raise `iter_max`", fixed = TRUE)
    expect_identical(attr(r$value, "partition")[c("iter", "iter_max", "converged")],
                     list(iter = 40L, iter_max = 40L, converged = FALSE), label = nm)
    # A limit above the default changes nothing for a start that converges.
    expect_identical(partition_xy(voronoi_seeds_kmeans(pts, cs$k, partition = "density",
                                                       init = cs$init, nstart = 1,
                                                       set_seed = cs$seed, iter_max = 400)),
                     full$centers, label = nm)
  }
})

test_that("get_voronoi_seeds() runs a density start for kmeans_iter steps, 100 by default", {
  xy  <- partition_uniform_xy(1200L, seed = 41L)
  pts <- partition_points_sf(xy)
  for (nm in names(.pv_cases)) {
    cs <- .pv_cases[[nm]]
    call <- function(...) get_voronoi_seeds(method = "kmeans", n = cs$k, sample_points = pts,
                                            partition = "density", init = cs$init,
                                            kmeans_nstart = 1, set_seed = cs$seed, ...)
    full <- .pv_engine(xy, cs$k, cs$seed, 100L, cs$init)
    r <- partition_warnings(call())
    expect_length(r$warnings, 0L)
    expect_identical(partition_xy(r$value), full$centers, label = nm)
    expect_true(attr(r$value, "kmeans")$converged, label = nm)
    expect_identical(attr(r$value, "kmeans")$iter, full$iter, label = nm)
    cut <- .pv_engine(xy, cs$k, cs$seed, 40L, cs$init)
    r <- partition_warnings(call(kmeans_iter = 40))
    expect_identical(partition_xy(r$value), cut$centers, label = nm)
    expect_length(r$warnings, 1L)
    w <- r$warnings[[1L]]
    expect_s3_class(w, "spatialkit_not_converged")
    expect_identical(w[c("level", "iter", "iter_max")],
                     list(level = cs$k, iter = 40L, iter_max = 40L))
    expect_match(conditionMessage(w), "get_voronoi_seeds(): at", fixed = TRUE)
    expect_match(conditionMessage(w), "Raise `kmeans_iter`", fixed = TRUE)
    expect_false(attr(r$value, "kmeans")$converged, label = nm)
    expect_identical(attr(r$value, "kmeans")$iter, 40L, label = nm)
    expect_identical(attr(r$value, "partition")[c("iter", "iter_max", "converged")],
                     list(iter = 40L, iter_max = 40L, converged = FALSE), label = nm)
  }
})

test_that("resolution_profile() runs a density start for up to 100 steps, or iter_max", {
  # Geometry only, so nothing draws from the stream before the sweep.
  xy  <- partition_uniform_xy(1200L, seed = 41L)
  pts <- partition_points_sf(xy)
  for (nm in names(.pv_cases)) {
    cs <- .pv_cases[[nm]]
    key <- as.character(cs$k)
    call <- function(...) resolution_profile(pts, levels = cs$k, nstart = 1, seed = cs$seed,
                                             partition = "density", init = cs$init, ...)
    full <- .pv_engine(xy, cs$k, cs$seed, 100L, cs$init)
    r <- partition_warnings(call())
    expect_length(partition_of_class(r, "spatialkit_not_converged"), 0L)
    expect_identical(attr(r$value, "centres")[[key]], full$centers, label = nm)
    rec <- attr(r$value, "partition")
    expect_identical(rec$iter_max, 100L)
    expect_identical(rec$iter, stats::setNames(full$iter, key))
    expect_identical(rec$converged, stats::setNames(TRUE, key))
    cut <- .pv_engine(xy, cs$k, cs$seed, 40L, cs$init)
    r <- partition_warnings(call(iter_max = 40))
    expect_identical(attr(r$value, "centres")[[key]], cut$centers, label = nm)
    w <- partition_of_class(r, "spatialkit_not_converged")
    expect_length(w, 1L)
    expect_identical(w[[1L]][c("level", "iter", "iter_max")],
                     list(level = cs$k, iter = 40L, iter_max = 40L))
    expect_match(conditionMessage(w[[1L]]), sprintf(
      "resolution_profile(): at 1 of the 1 levels (%d cells) the start that was kept had not converged after 40 steps",
      cs$k), fixed = TRUE)
    rec <- attr(r$value, "partition")
    expect_identical(rec$iter_max, 40L)
    expect_identical(rec$iter, stats::setNames(40L, key))
    expect_identical(rec$converged, stats::setNames(FALSE, key))
  }
})

test_that("a profile names the levels that were cut off, in one warning", {
  xy  <- partition_clustered_xy(n_bg = 150L, n_cl = 350L, seed = 7L)
  pts <- partition_points_sf(xy)
  lv  <- c(2L, 9L, 30L, 60L)
  # The sweep, level after level in one stream, as the profile runs it.
  sweep <- function(limit) {
    set.seed(11)
    lapply(lv, function(L) spatialkit:::.partition_centres(xy, L, partition = "density",
                                                           nstart = 2L, iter.max = limit))
  }
  limit <- 9L
  eng <- sweep(limit)
  conv <- vapply(eng, function(e) e$converged, logical(1))
  # Some levels stop within 9 steps and some do not.
  expect_true(any(conv) && any(!conv))
  r <- partition_warnings(resolution_profile(pts, levels = lv, nstart = 2, seed = 11,
                                             partition = "density", iter_max = limit))
  w <- partition_of_class(r, "spatialkit_not_converged")
  expect_length(w, 1L)
  expect_identical(w[[1L]]$level, lv[!conv])
  expect_identical(w[[1L]]$iter, rep(limit, sum(!conv)))
  expect_identical(w[[1L]]$iter_max, limit)
  expect_match(conditionMessage(w[[1L]]),
               sprintf("at %d of the 4 levels (%s cells)", sum(!conv),
                       paste(lv[!conv], collapse = ", ")), fixed = TRUE)
  rec <- attr(r$value, "partition")
  expect_identical(rec$converged, stats::setNames(conv, lv))
  expect_identical(rec$iter, stats::setNames(vapply(eng, function(e) e$iter, integer(1)), lv))
  # The centres are the engine's at that limit, cut off or not.
  expect_identical(unname(attr(r$value, "centres")), lapply(eng, function(e) e$centers))
  # With room to finish there is nothing to say.
  r <- partition_warnings(resolution_profile(pts, levels = lv, nstart = 2, seed = 11,
                                             partition = "density"))
  expect_length(partition_of_class(r, "spatialkit_not_converged"), 0L)
  expect_true(all(attr(r$value, "partition")$converged))
})

test_that("the k-means kind from a starting rule that is not the function's own says it too", {
  xy  <- partition_clustered_xy(n_bg = 150L, n_cl = 350L, seed = 7L)
  pts <- partition_points_sf(xy)
  # A seeding function from k-means++ starts.
  r <- partition_warnings(voronoi_seeds_kmeans(pts, 40, init = "kmeans++", nstart = 2,
                                               iter_max = 1, partition = "kmeans"))
  expect_length(r$warnings, 1L)
  w <- r$warnings[[1L]]
  expect_s3_class(w, "spatialkit_not_converged")
  expect_identical(w[c("level", "iter", "iter_max")], list(level = 40L, iter = 1L, iter_max = 1L))
  expect_match(conditionMessage(w),
               "had not converged after 1 iteration of k-means (partition = \"kmeans\")",
               fixed = TRUE)
  expect_identical(attr(r$value, "partition")[c("iter", "iter_max", "converged")],
                   list(iter = 1L, iter_max = 1L, converged = FALSE))
  r <- partition_warnings(get_voronoi_seeds(method = "kmeans", n = 40, sample_points = pts,
                                            init = "kmeans++", kmeans_nstart = 2,
                                            kmeans_iter = 1, set_seed = 2, partition = "kmeans"))
  expect_length(partition_of_class(r, "spatialkit_not_converged"), 1L)
  expect_false(attr(r$value, "kmeans")$converged)
  # A profile from random starts.
  r <- partition_warnings(resolution_profile(pts, levels = c(10L, 40L), nstart = 2,
                                             init = "random", iter_max = 1, estimate = NULL,
                                             obs_per_cell = NULL, sac = NULL,
                                             partition = "kmeans"))
  w <- partition_of_class(r, "spatialkit_not_converged")
  expect_length(w, 1L)
  expect_identical(w[[1L]]$level, c(10L, 40L))
  expect_identical(attr(r$value, "partition")$converged, c("10" = FALSE, "40" = FALSE))
  expect_identical(attr(r$value, "partition")$iter_max, 1L)
  # Left to their own limits these converge and say nothing.
  expect_length(partition_warnings(voronoi_seeds_kmeans(pts, 40, init = "kmeans++",
                                                        partition = "kmeans"))$warnings, 0L)
  r <- partition_warnings(resolution_profile(pts, levels = c(10L, 40L), nstart = 2,
                                             init = "random", estimate = NULL,
                                             obs_per_cell = NULL, sac = NULL,
                                             partition = "kmeans"))
  expect_length(partition_of_class(r, "spatialkit_not_converged"), 0L)
  expect_identical(attr(r$value, "partition")$iter_max, 100L)
})

test_that("the runs the functions have always made report as they always did", {
  xy  <- partition_clustered_xy(n_bg = 150L, n_cl = 350L, seed = 7L)
  pts <- partition_points_sf(xy)
  # The seeding functions' k-means run is one stats::kmeans() call, and a limit
  # that cuts it off gets that function's warning and no record.
  r <- partition_warnings(voronoi_seeds_kmeans(pts, 40, nstart = 2, iter_max = 1,
                                               partition = "kmeans"))
  expect_gt(length(r$warnings), 0L)
  expect_true(all(vapply(r$warnings, function(w)
    grepl("did not converge in 1 iteration", conditionMessage(w), fixed = TRUE), logical(1))))
  expect_length(partition_of_class(r, "spatialkit_not_converged"), 0L)
  expect_null(attr(r$value, "partition"))
  set.seed(456)
  km <- suppressWarnings(stats::kmeans(partition_xy(pts), centers = 40, iter.max = 1, nstart = 2L))
  expect_identical(partition_xy(r$value), unname(km$centers))
  r <- partition_warnings(get_voronoi_seeds(method = "kmeans", n = 40, sample_points = pts,
                                            kmeans_nstart = 2, kmeans_iter = 1, set_seed = 2,
                                            partition = "kmeans"))
  expect_gt(length(r$warnings), 0L)
  expect_length(partition_of_class(r, "spatialkit_not_converged"), 0L)
  expect_named(attr(r$value, "kmeans"),
               c("cluster", "rows", "size", "withinss", "tot_withinss", "iter", "nstart"))

  # The default profile has never reported a k-means start that ran out of
  # iterations, and does not now: every start is stood in for by one marked
  # as cut off, and nothing is raised or recorded.  Naming a limit is asking,
  # and is answered.
  real <- spatialkit:::.partition_centres
  local_mocked_bindings(.partition_centres = function(...) {
    fit <- real(...)
    fit$converged <- FALSE
    fit
  })
  r <- partition_warnings(resolution_profile(pts, levels = c(5L, 12L), nstart = 2,
                                             estimate = NULL, obs_per_cell = NULL, sac = NULL,
                                             partition = "kmeans"))
  expect_length(partition_of_class(r, "spatialkit_not_converged"), 0L)
  expect_null(attr(r$value, "partition"))
  r <- partition_warnings(resolution_profile(pts, levels = c(5L, 12L), nstart = 2,
                                             iter_max = 50, estimate = NULL,
                                             obs_per_cell = NULL, sac = NULL,
                                             partition = "kmeans"))
  w <- partition_of_class(r, "spatialkit_not_converged")
  expect_length(w, 1L)
  expect_identical(w[[1L]]$level, c(5L, 12L))
  expect_identical(w[[1L]]$iter_max, 50L)
  expect_null(attr(r$value, "partition"))
})

test_that("iter_max is NULL or a positive whole number", {
  pts <- partition_points_sf(partition_clustered_xy(n_bg = 40L, n_cl = 80L, seed = 28L))
  calls <- list(
    resolution_profile = function(...) resolution_profile(pts, n_levels = 3, nstart = 2, ...),
    voronoi_seeds_kmeans = function(...) voronoi_seeds_kmeans(pts, 5, ...))
  for (nm in names(calls)) {
    f <- calls[[nm]]
    for (bad in list(0, -1, 2.5, NA, "a", c(10, 20), Inf))
      expect_error(f(iter_max = bad),
                   sprintf("%s(): `iter_max` must be NULL or a single positive whole number of iterations", nm),
                   fixed = TRUE, info = nm)
    expect_error(f(iter_max = 0), "and at least 1; got 0", fixed = TRUE, info = nm)
    expect_no_error(f(iter_max = 200))
    expect_no_error(f(partition = "density", iter_max = 200))
  }
  # The limit is used: a density start held to one step is not the full run.
  one  <- partition_quietly(voronoi_seeds_kmeans(pts, 8, partition = "density", iter_max = 1))
  full <- voronoi_seeds_kmeans(pts, 8, partition = "density")
  expect_gt(max(abs(partition_xy(one) - partition_xy(full))), 1)
})
