# tests/testthat/test-core-partition-pipeline.R
# ---------------------------------------------------------------------------
# The density kind from end to end: resolution_profile() builds its ladder
# with it, select_resolution() picks a level, the seeding functions hand on
# that level's centres, build_tessellation() draws their Voronoi cells,
# assign_features_to_polygons() puts every observation in one of them and
# summarize_by_cell() counts them.  The regions that come out are the ones
# the profile scored, and every result says what kind built them.
# Also here: the arguments' validation and the classed warnings for scored
# centres of another kind and for a centre that holds no point.
# ---------------------------------------------------------------------------

.pl_pts <- function(n_bg = 150L, n_cl = 450L, seed = 21L, crs = 32632,
                    origin = c(5e5, 5e6), side = 4000) {
  xy <- partition_clustered_xy(n_bg = n_bg, n_cl = n_cl, clusters = 4L, side = side,
                               origin = origin, seed = seed)
  cleanup <- spatialkit:::.with_seed(seed + 1L)
  on.exit(cleanup(), add = TRUE)
  z <- sin(xy[, 1L] / 700) + cos(xy[, 2L] / 900) + stats::rnorm(nrow(xy), sd = 0.8)
  partition_points_sf(xy, crs = crs, z = z)
}

.pl_record <- list(kind = "density", init = "kmeans++", density_power = 1, density_k = 8L)

# The pipeline after the profile, with the checks every case shares.
.pl_run <- function(pts, prof, n, level) {
  seeds <- partition_quietly(get_voronoi_seeds(method = "kmeans", n = n, sample_points = pts))
  expect_identical(nrow(seeds), as.integer(level))
  expect_identical(attr(seeds, "partition"), .pl_record)
  tess <- build_tessellation(seeds, method = "voronoi", quiet = TRUE)
  expect_identical(tess$params$partition, .pl_record)
  expect_identical(tess$params$boundary_source, "seeds' cover")
  expect_identical(nrow(tess$cells), as.integer(level))
  # One cell per seed, and each seed in its own cell.
  expect_identical(sort(tess$index), sort(tess$cells$cell_id))

  assigned <- assign_features_to_polygons(pts, tess$cells, polygon_id_col = "cell_id",
                                          quiet = TRUE)
  # Every observation in exactly one cell.
  expect_identical(nrow(assigned), nrow(pts))
  expect_false(anyNA(assigned$cell_id))
  expect_true(all(assigned$cell_id %in% tess$cells$cell_id))
  inside <- lengths(sf::st_intersects(pts, tess$cells))
  expect_true(all(inside >= 1L))
  expect_lte(sum(inside > 1L), 0L)                # none on a shared edge here

  summ <- summarize_by_cell(assigned, response_var = "z", id_col = "cell_id", estimate = "mean")
  expect_identical(sum(summ$n), nrow(pts))
  expect_identical(nrow(summ), as.integer(level))
  list(seeds = seeds, tess = tess, assigned = assigned, summ = summ)
}


test_that("the density kind runs from profile to summary, with the counts the profile reported", {
  pts  <- .pl_pts()
  prof <- partition_quietly(resolution_profile(pts, "z", n_levels = 8, nstart = 4,
                                               partition = "density", estimate = NULL,
                                               obs_per_cell = NULL, sac = NULL))
  expect_identical(attr(prof, "partition")[1:4], .pl_record)
  expect_named(attr(prof, "partition")$objective, as.character(prof$levels))
  expect_true(all(is.finite(attr(prof, "partition")$objective)))
  expect_true(all(diff(attr(prof, "partition")$objective) < 0))
  # Every column a k-means profile has, and every level scored.
  expect_true(all(is.finite(prof$wss)) && all(is.finite(prof$cp)) && all(is.finite(prof$rss)))
  expect_true(all(diff(prof$wss) < 0))
  expect_output(print(prof), "4 density-weighted k-means++ starts (power 1, 8 neighbours) per level",
                fixed = TRUE)

  sel <- select_resolution(prof, criterion = "cp")
  expect_identical(sel$partition, .pl_record)
  expect_identical(attr(sel$seeds, "partition"), .pl_record)
  level <- sel$best
  out <- .pl_run(pts, prof, sel, level)

  # The seeds are the profile's centres at that level ...
  cen <- attr(prof, "centres")[[as.character(level)]]
  expect_equal(partition_xy(out$seeds), cen, tolerance = 1e-12)
  # ... so the cells hold the points the profile counted: its smallest and
  # median cell, and the whole set of counts of the scored partition.
  counts <- as.integer(table(out$assigned$cell_id))
  row <- prof[prof$levels == level, ]
  expect_identical(min(counts), row$cell_n_min)
  expect_equal(stats::median(counts), row$cell_n_median)
  scored <- tabulate(partition_nearest(partition_xy(pts), cen), nbins = level)
  expect_identical(sort(counts), sort(scored))
  expect_identical(sort(as.integer(out$summ$n)), sort(scored))
  # Point by point: an observation is in the cell of its nearest centre.
  # `index` gives the cell each seed sits in.
  nearest <- partition_nearest(partition_xy(pts), partition_xy(out$seeds))
  expect_identical(out$assigned$cell_id, out$tess$index[nearest])
  # And the residual sum of squares the profile scored is the cells'.
  rss <- sum(tapply(pts$z, out$assigned$cell_id, function(v) sum((v - mean(v))^2)))
  expect_equal(rss, row$rss, tolerance = 1e-10)
  # The record of the seeds describes the cells: sizes and sum of squares.
  km <- attr(out$seeds, "kmeans")
  expect_identical(sort(km$size), sort(scored))
  expect_equal(km$tot_withinss, row$wss, tolerance = 1e-10)
})

test_that("every level of a density ladder gives the cells it was scored on", {
  pts  <- .pl_pts(seed = 22L)
  prof <- resolution_profile(pts, levels = c(3L, 12L, 40L), nstart = 3, partition = "density")
  xy <- partition_xy(pts)
  for (i in seq_len(nrow(prof))) {
    L <- prof$levels[i]
    seeds <- spatialkit:::.profile_seeds(prof, L)
    expect_identical(attr(seeds, "partition"), .pl_record)
    s <- voronoi_seeds_kmeans(pts, structure(list(best = L, seeds = seeds, criterion = "cp",
                                                  edge = NA_character_, structured = TRUE),
                                             class = "resolution_selection"))
    expect_identical(attr(s, "partition"), .pl_record)
    tess <- build_tessellation(s, method = "voronoi", quiet = TRUE)
    expect_false(anyNA(tess$index))
    counts <- tabulate(match(
      assign_features_to_polygons(pts, tess$cells, polygon_id_col = "cell_id",
                                  quiet = TRUE)$cell_id,
      tess$cells$cell_id), nbins = L)
    expect_identical(min(counts), prof$cell_n_min[i])
    expect_equal(stats::median(counts), prof$cell_n_median[i])
  }
})

test_that("a profile passed as the count, and lon/lat points, go through the same way", {
  pts  <- .pl_pts(n_bg = 100L, n_cl = 300L, seed = 23L)
  prof <- partition_quietly(resolution_profile(pts, "z", n_levels = 6, nstart = 3,
                                               partition = "density", estimate = NULL,
                                               obs_per_cell = NULL, sac = NULL))
  level <- select_resolution(prof, criterion = "cp")$best
  out <- .pl_run(pts, prof, prof, level)
  expect_match(attr(out$seeds, "n_from"), "resolution_profile() read with", fixed = TRUE)

  # Lon/lat points: the profile works in a projected CRS, the seeds come back
  # in degrees, and the cells still hold every observation.
  ll <- sf::st_transform(pts, 4326)
  prof_ll <- partition_quietly(resolution_profile(ll, "z", n_levels = 6, nstart = 3,
                                                  partition = "density", estimate = NULL,
                                                  obs_per_cell = NULL, sac = NULL))
  sel_ll <- select_resolution(prof_ll, criterion = "cp")
  seeds <- partition_quietly(get_voronoi_seeds(method = "kmeans", n = sel_ll, sample_points = ll))
  expect_true(sf::st_is_longlat(seeds))
  expect_identical(attr(seeds, "partition"), .pl_record)
  tess <- partition_quietly(build_tessellation(seeds, method = "voronoi", quiet = TRUE))
  expect_identical(tess$params$partition, .pl_record)
  assigned <- partition_quietly(assign_features_to_polygons(
    sf::st_transform(ll, sf::st_crs(tess$cells)), tess$cells, polygon_id_col = "cell_id",
    quiet = TRUE))
  expect_identical(nrow(assigned), nrow(ll))
  expect_false(anyNA(assigned$cell_id))
  counts <- as.integer(table(assigned$cell_id))
  row <- prof_ll[prof_ll$levels == sel_ll$best, ]
  expect_identical(min(counts), row$cell_n_min)
  expect_equal(stats::median(counts), row$cell_n_median)
})

test_that("fresh density seeds go through the tessellation with their record", {
  pts <- .pl_pts(seed = 24L)
  rec_rnd <- list(kind = "density", init = "random", density_power = 0.5, density_k = 5L)
  s <- voronoi_seeds_kmeans(pts, 30, partition = "density", init = "random",
                            density_power = 0.5, density_k = 5, nstart = 10)
  # What built the seeds, and how the start that was kept ended: the steps
  # it took, the limit, and whether it stopped before the limit.
  full <- attr(s, "partition")
  expect_named(full, c("kind", "init", "density_power", "density_k",
                       "iter", "iter_max", "converged"))
  expect_identical(full[1:4], rec_rnd)
  expect_identical(full$iter_max, 100L)
  expect_true(full$converged)
  expect_true(is.integer(full$iter) && full$iter >= 1L && full$iter < 100L)
  expect_identical(unique(s$method), "kmeans")
  expect_identical(s$seed_id, 1:30)
  expect_s3_class(attr(s, "cover"), "sfc")
  tess <- build_tessellation(s, method = "voronoi", quiet = TRUE)
  expect_identical(tess$params$partition, full)
  expect_identical(create_voronoi_polygons(s, quiet = TRUE)$params$partition, full)
  expect_false(anyNA(tess$index))
  # The record survives a transform and a subset, as the cover does.
  expect_identical(attr(sf::st_transform(s, 3857), "partition"), full)
  expect_identical(attr(s[1:10, ], "partition"), full)
  # A CRS given to the tessellation moves the seeds; the record still arrives.
  expect_identical(build_tessellation(s, method = "voronoi", crs = 3857,
                                      quiet = TRUE)$params$partition, full)
  # So it does when build_tessellation() projects the seeds itself: seeds in
  # lon/lat, and seeds with no CRS whose coordinates are degrees.
  ll <- sf::st_transform(pts, 4326)
  s_ll <- voronoi_seeds_kmeans(ll, 30, partition = "density", init = "random",
                               density_power = 0.5, density_k = 5, nstart = 10)
  expect_identical(attr(s_ll, "partition")[1:4], rec_rnd)
  expect_identical(partition_quietly(build_tessellation(s_ll, method = "voronoi",
                                                        quiet = TRUE))$params$partition,
                   attr(s_ll, "partition"))
  s_na <- partition_quietly(voronoi_seeds_kmeans(sf::st_set_crs(ll, NA), 30,
                                                 partition = "density", init = "random",
                                                 density_power = 0.5, density_k = 5, nstart = 10))
  expect_true(is.na(sf::st_crs(s_na)))
  expect_identical(attr(s_na, "partition")[1:4], rec_rnd)
  expect_identical(partition_quietly(build_tessellation(s_na, method = "voronoi",
                                                        quiet = TRUE))$params$partition,
                   attr(s_na, "partition"))
  # Read by one helper wherever it sits.
  expect_identical(spatialkit:::.partition_record(s), rec_rnd)
  expect_identical(spatialkit:::.partition_record(tess), rec_rnd)
  prof <- resolution_profile(pts, levels = c(5L, 9L), nstart = 2, partition = "density",
                             estimate = NULL, obs_per_cell = NULL, sac = NULL)
  expect_identical(spatialkit:::.partition_record(prof), .pl_record)

  # get_voronoi_seeds(): the record of the run describes the seeds' cells.
  g <- get_voronoi_seeds(method = "kmeans", n = 30, sample_points = pts,
                         partition = "density", set_seed = 4, kmeans_nstart = 10)
  km <- attr(g, "kmeans")
  expect_named(km, c("cluster", "rows", "size", "withinss", "tot_withinss", "iter", "nstart",
                     "objective", "converged"))
  expect_true(km$converged)
  expect_identical(attr(g, "partition"),
                   c(.pl_record, list(iter = km$iter, iter_max = 100L, converged = TRUE)))
  xy <- partition_xy(pts)
  expect_identical(km$cluster, as.integer(partition_nearest(xy, partition_xy(g))))
  expect_identical(km$size, tabulate(km$cluster, nbins = 30L))
  mu <- rowsum(xy, km$cluster) / km$size
  expect_equal(km$tot_withinss, sum((xy - mu[km$cluster, ])^2), tolerance = 1e-10)
  expect_gt(km$objective, 0)
  expect_identical(km$nstart, 10L)
})

test_that("a profile fitted to a subsample still tessellates the whole layer", {
  pts  <- .pl_pts(n_bg = 400L, n_cl = 1200L, seed = 25L)
  prof <- partition_quietly(resolution_profile(pts, "z", n_levels = 5, nstart = 2,
                                               sample_n = 500, partition = "density",
                                               estimate = NULL, obs_per_cell = NULL,
                                               sac = NULL))
  expect_identical(attr(prof, "bounds")$n_sample, 500L)
  sel <- select_resolution(prof, criterion = "cp")
  seeds <- partition_quietly(get_voronoi_seeds(method = "kmeans", n = sel, sample_points = pts))
  tess <- build_tessellation(seeds, method = "voronoi", quiet = TRUE)
  assigned <- assign_features_to_polygons(pts, tess$cells, polygon_id_col = "cell_id",
                                          quiet = TRUE)
  expect_identical(nrow(assigned), nrow(pts))
  expect_false(anyNA(assigned$cell_id))
  expect_identical(sum(summarize_by_cell(assigned, "z", id_col = "cell_id",
                                         estimate = "mean")$n), nrow(pts))
  expect_identical(tess$params$partition, .pl_record)
})


# ---- The profile's columns under the density kind ---------------------------------

test_that("reliability does not read the cells, and Cp counts them as it does for k-means", {
  pts <- .pl_pts(seed = 26L)
  # One variogram for both profiles, written by hand so that it has a range
  # and a nugget whatever is installed: only the cells differ.
  vm <- data.frame(model = c("Nug", "Exp"), psill = c(0.6, 0.5), range = c(0, 300))
  sac <- structure(900, class = c("sac_range", "numeric"), variogram_model = vm,
                   crs = sf::st_crs(pts), detrended = FALSE)
  lv  <- c(4L, 10L, 25L, 60L)
  p_km  <- partition_quietly(resolution_profile(pts, "z", levels = lv, nstart = 3, sac = sac,
                                                estimate = NULL, obs_per_cell = NULL,
                                                partition = "kmeans"))
  p_den <- partition_quietly(resolution_profile(pts, "z", levels = lv, nstart = 3, sac = sac,
                                                partition = "density", estimate = NULL,
                                                obs_per_cell = NULL))
  expect_true(all(is.finite(p_km$reliability)))
  expect_identical(p_den$reliability, p_km$reliability)
  expect_identical(attr(p_den, "bounds"), attr(p_km, "bounds"))
  expect_identical(attr(p_den, "one_cell")[c("n", "rss", "cp", "constant")],
                   attr(p_km, "one_cell")[c("n", "rss", "cp", "constant")])
  # Cp is the formula of the k-means kind on the density kind's cells, with
  # the variogram's nugget.
  cn <- attr(p_den, "cp_noise")
  expect_identical(cn$source, "variogram nugget")
  expect_identical(cn$value, 0.6)
  n <- nrow(pts)
  expect_equal(p_den$cp, p_den$rss / n + cn$value * 2 * lv / n, tolerance = 1e-12)
  # The test of the cells against one mean is the F test on these cells.
  f_p <- attr(p_den, "one_cell")$f_p
  tss <- attr(p_den, "one_cell")$rss
  expect_equal(unname(f_p),
               stats::pf(((tss - p_den$rss) / (lv - 1)) / (p_den$rss / (n - lv)),
                         lv - 1, n - lv, lower.tail = FALSE), tolerance = 1e-10)
  expect_false(isTRUE(all.equal(p_den$rss, p_km$rss)))
  # k-means minimises the unweighted sum of squares, so the density kind's is
  # no lower at any level, and its weighted sum is not that sum.
  expect_true(all(p_den$wss > p_km$wss))
  expect_true(all(attr(p_den, "partition")$objective != p_den$wss))
})

test_that("the elbow is read off the unweighted curve, and uniform points have none", {
  cleanup <- spatialkit:::.with_seed(30L)
  xy <- cbind(stats::runif(600, 0, 1000), stats::runif(600, 0, 1000))
  cleanup()
  u <- partition_points_sf(xy)
  p <- resolution_profile(u, n_levels = 8, nstart = 3, partition = "density")
  expect_true(all(is.na(p$elbow)))
  expect_error(select_resolution(p, "elbow"), "no elbow")
  # Clustered points have one, and it is the sag of the profile's own wss.
  cl <- .pl_pts(seed = 27L)
  pc <- resolution_profile(cl, n_levels = 8, nstart = 3, partition = "density")
  expect_true(any(is.finite(pc$elbow)))
  xs  <- partition_xy(cl)
  tss <- sum(sweep(xs, 2L, colMeans(xs))^2)
  rd  <- spatialkit:::.elbow_read(c(1L, pc$levels), c(tss, pc$wss))
  expect_identical(pc$elbow, rd$sag[-1L])
  sel <- select_resolution(pc, "elbow")
  expect_identical(sel$partition, .pl_record)
})


test_that("a level with an empty cell is reported, and Cp counts the occupied cells", {
  # A centre of the density kind that loses its points stays where it is.
  # The engine is stood in for here, so that the level has an empty centre
  # whatever the starts do: its partition into four cells is handed back
  # with a fifth centre that holds no point.  The last two tests of this
  # file reach the same warning with nothing stood in for.
  pts <- .pl_pts(n_bg = 60L, n_cl = 120L, seed = 31L)
  real <- spatialkit:::.partition_centres
  local_mocked_bindings(.partition_centres = function(xy, k, ...) {
    fit <- real(xy, k - 1L, ...)
    fit$centers  <- rbind(fit$centers, c(min(xy[, 1L]) - 1e4, min(xy[, 2L]) - 1e4))
    fit$size     <- c(fit$size, 0L)
    fit$withinss <- c(fit$withinss, 0)
    fit
  })
  vm <- data.frame(model = c("Nug", "Exp"), psill = c(0.6, 0.5), range = c(0, 300))
  sac <- structure(900, class = c("sac_range", "numeric"), variogram_model = vm,
                   crs = sf::st_crs(pts), detrended = FALSE)
  w <- NULL
  prof <- withCallingHandlers(
    resolution_profile(pts, "z", levels = 5L, nstart = 2, sac = sac, partition = "density"),
    spatialkit_empty_cells = function(c) { w <<- c; invokeRestart("muffleWarning") },
    warning = function(c) invokeRestart("muffleWarning"))
  expect_s3_class(w, "spatialkit_empty_cells")
  expect_identical(w$level, 5L)
  expect_identical(w$empty, 1L)
  expect_match(conditionMessage(w), "at 5 cells, 1 of the centres hold no point", fixed = TRUE)
  expect_identical(prof$cell_n_min, 0L)
  expect_identical(prof$cells_resp, 4L)
  # Four cell means, not five, in the variance term.
  n <- nrow(pts)
  expect_equal(prof$cp, prof$rss / n + 0.6 * (4 / n + 4 / n), tolerance = 1e-12)
})


# ---- Arguments -----------------------------------------------------------------------

test_that("the partition arguments are validated like their neighbours", {
  pts <- .pl_pts(n_bg = 40L, n_cl = 80L, seed = 28L)
  calls <- list(
    resolution_profile = function(...) resolution_profile(pts, n_levels = 3, nstart = 2, ...),
    voronoi_seeds_kmeans = function(...) voronoi_seeds_kmeans(pts, 5, ...),
    get_voronoi_seeds = function(...) get_voronoi_seeds(method = "kmeans", n = 5,
                                                        sample_points = pts, set_seed = 1, ...))
  for (nm in names(calls)) {
    f <- calls[[nm]]
    expect_error(f(partition = "equal"), "should be one of", info = nm)
    expect_error(f(init = "furthest"), "should be one of", info = nm)
    expect_error(f(partition = "density", density_power = -0.1),
                 sprintf("%s(): `density_power` must be a single number from 0 to 2 and at least 0", nm),
                 fixed = TRUE, info = nm)
    expect_error(f(partition = "density", density_power = 2.5),
                 "`density_power` must be a single number from 0 to 2 and at most 2",
                 fixed = TRUE, info = nm)
    expect_error(f(partition = "density", density_power = "1"),
                 "`density_power` must be a single number from 0 to 2; got character",
                 fixed = TRUE, info = nm)
    expect_error(f(partition = "density", density_power = c(1, 2)),
                 "`density_power` must be a single number from 0 to 2; got numeric of length 2",
                 fixed = TRUE, info = nm)
    expect_error(f(partition = "density", density_power = NA_real_),
                 "`density_power` must be a single number from 0 to 2", fixed = TRUE, info = nm)
    expect_error(f(partition = "density", density_k = 0),
                 "`density_k` must be a single positive whole number of neighbours and at least 1",
                 fixed = TRUE, info = nm)
    expect_error(f(partition = "density", density_k = 2.5),
                 "`density_k` must be a single positive whole number of neighbours; got 2.5",
                 fixed = TRUE, info = nm)
    expect_error(f(partition = "density", density_k = Inf),
                 "`density_k` must be a single positive whole number of neighbours",
                 fixed = TRUE, info = nm)
    # The ends of the range are allowed.
    expect_no_error(f(partition = "density", density_power = 0, density_k = 1))
    expect_no_error(f(partition = "density", density_power = 2, density_k = 500))
    # Set where they are not read: said, and the result is the k-means
    # kind's without them.
    expect_warning(a <- f(density_power = 0.5, partition = "kmeans"),
                   sprintf("%s(): `density_power` is read by partition = \"density\" only", nm),
                   fixed = TRUE)
    expect_identical(a, f(partition = "kmeans"))
    expect_warning(f(density_power = 0.5, density_k = 4, partition = "kmeans"),
                   "`density_power` and `density_k` are read by partition = \"density\" only",
                   fixed = TRUE)
  }
  expect_no_error(voronoi_seeds_kmeans(pts, 5, init = NULL))
  # A value of `init` that is not a string names the function and the
  # argument, where match.arg() names neither.
  for (nm in names(calls)) {
    f <- calls[[nm]]
    expect_error(f(init = 1), sprintf("%s(): `init` should be one of \"kmeans++\", \"random\"", nm),
                 fixed = TRUE, info = nm)
    expect_error(f(init = 1), "got numeric", fixed = TRUE, info = nm)
    expect_error(f(init = "furthest"), "got \"furthest\"", fixed = TRUE, info = nm)
    expect_error(f(init = c("random", "kmeans++")), "got character of length 2", fixed = TRUE,
                 info = nm)
    expect_error(f(init = NA_character_), "`init` should be one of", fixed = TRUE, info = nm)
  }
  expect_error(voronoi_seeds_kmeans(pts, 5, init = 1), "\"random\", or NULL; got numeric",
               fixed = TRUE)
  # As match.arg() does, a leading part of a name is the name.
  expect_identical(voronoi_seeds_kmeans(pts, 5, init = "rand", partition = "kmeans"),
                   voronoi_seeds_kmeans(pts, 5, partition = "kmeans"))
  expect_identical(attr(voronoi_seeds_kmeans(pts, 5, init = "kmeans"), "partition")$init,
                   "kmeans++")
})

test_that("methods that draw no k-means say that the partition arguments are not theirs", {
  pts <- .pl_pts(n_bg = 40L, n_cl = 80L, seed = 28L)
  bnd <- sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(pts)))
  expect_warning(r <- get_voronoi_seeds(bnd, method = "random", n = 5, set_seed = 2,
                                        partition = "density"),
                 "`partition` and `init` are read by method = \"kmeans\" only", fixed = TRUE)
  expect_identical(r, get_voronoi_seeds(bnd, method = "random", n = 5, set_seed = 2))
  expect_null(attr(r, "partition"))
  expect_warning(get_voronoi_seeds(method = "provided", seeds = pts[1:4, ], init = "kmeans++"),
                 "method = \"provided\" ignores them", fixed = TRUE)
})

test_that("scored centres are returned whatever kind the call names, and the call is told", {
  pts <- .pl_pts(n_bg = 100L, n_cl = 300L, seed = 29L)
  p_km  <- resolution_profile(pts, levels = c(6L, 14L), nstart = 2, estimate = NULL,
                              obs_per_cell = NULL, sac = NULL, partition = "kmeans")
  p_den <- resolution_profile(pts, levels = c(6L, 14L), nstart = 2, partition = "density",
                              estimate = NULL, obs_per_cell = NULL, sac = NULL)
  sel_of <- function(p, L) structure(
    list(best = L, seeds = spatialkit:::.profile_seeds(p, L), criterion = "cp",
         edge = NA_character_, structured = TRUE), class = "resolution_selection")
  km14 <- sel_of(p_km, 14L); den14 <- sel_of(p_den, 14L)

  # A k-means selection with the density kind asked for: k-means centres.
  w <- tryCatch(voronoi_seeds_kmeans(pts, km14, partition = "density"), warning = identity)
  expect_s3_class(w, "spatialkit_partition_unused")
  expect_identical(w$kind, "kmeans")
  expect_identical(w$init, "kmeans++")
  s <- suppressWarnings(voronoi_seeds_kmeans(pts, km14, partition = "density"))
  expect_equal(partition_xy(s), attr(p_km, "centres")[["14"]], tolerance = 1e-12)
  expect_null(attr(s, "partition"))
  # A density selection with k-means named, or another starting rule.
  expect_warning(g <- get_voronoi_seeds(method = "kmeans", n = den14, sample_points = pts,
                                        partition = "kmeans"),
                 class = "spatialkit_partition_unused")
  expect_identical(attr(g, "partition"), .pl_record)
  expect_equal(partition_xy(g), attr(p_den, "centres")[["14"]], tolerance = 1e-12)
  expect_warning(voronoi_seeds_kmeans(pts, den14, init = "random"),
                 "built with partition = \"density\" from k-means++ starts (power 1, 8 neighbours)",
                 fixed = TRUE)
  # The density kind with another power than the profile's.
  expect_warning(voronoi_seeds_kmeans(pts, den14, partition = "density", density_power = 0.5),
                 class = "spatialkit_partition_unused")
  expect_warning(voronoi_seeds_kmeans(pts, den14, partition = "density", density_k = 4),
                 class = "spatialkit_partition_unused")
  # Unnamed, or named as built: nothing to say.
  expect_no_warning(voronoi_seeds_kmeans(pts, den14))
  expect_no_warning(voronoi_seeds_kmeans(pts, den14, partition = "density", init = "kmeans++"))
  expect_no_warning(get_voronoi_seeds(method = "kmeans", n = km14, sample_points = pts))
  # A count as a number is a fresh run of the kind asked for.
  fresh <- voronoi_seeds_kmeans(pts, 14, partition = "density", nstart = 10)
  expect_false(isTRUE(all.equal(partition_xy(fresh), attr(p_km, "centres")[["14"]])))
})

test_that("the partition labels read as the print method prints them", {
  lab <- spatialkit:::.partition_label
  expect_identical(lab(list(kind = "density", init = "kmeans++", density_power = 1,
                            density_k = 8L), 25L),
                   "25 density-weighted k-means++ starts (power 1, 8 neighbours)")
  expect_identical(lab(list(kind = "kmeans", init = "random"), 5L), "5 random starts")
  expect_identical(lab(list(kind = "kmeans", init = "random"), 1L), "1 random start")
  expect_identical(lab(list(kind = "density", init = "random", density_power = 0.5,
                            density_k = 3L)),
                   "density-weighted random starts (power 0.5, 3 neighbours)")
  expect_identical(lab(list(kind = "density", init = "kmeans++", density_power = 2,
                            density_k = 1L), 3L),
                   "3 density-weighted k-means++ starts (power 2, 1 neighbour)")
  pts <- .pl_pts(n_bg = 40L, n_cl = 80L, seed = 28L)
  expect_output(print(resolution_profile(pts, n_levels = 3, nstart = 2, init = "random",
                                         estimate = NULL, obs_per_cell = NULL, sac = NULL,
                                         partition = "kmeans")),
                "geometry only; 2 random starts per level", fixed = TRUE)
  # A selection names the kind of its cells when the profile recorded one,
  # and prints as it always did when it did not.
  den <- partition_quietly(resolution_profile(pts, "z", n_levels = 4, nstart = 2,
                                              partition = "density", estimate = NULL,
                                              obs_per_cell = NULL, sac = NULL))
  expect_output(print(select_resolution(den, criterion = "cp")),
                "  cells       : density-weighted k-means++ starts (power 1, 8 neighbours)",
                fixed = TRUE)
  km <- partition_quietly(resolution_profile(pts, "z", n_levels = 4, nstart = 2,
                                             estimate = NULL, obs_per_cell = NULL, sac = NULL,
                                             partition = "kmeans"))
  expect_false(any(grepl("cells       :",
                         capture.output(print(select_resolution(km, criterion = "cp"))),
                         fixed = TRUE)))
})

test_that("a seeding function says when a density seed holds no point", {
  pts <- .pl_pts(n_bg = 60L, n_cl = 120L, seed = 31L)
  real <- spatialkit:::.partition_centres
  local_mocked_bindings(.partition_centres = function(xy, k, ...) {
    fit <- real(xy, k - 1L, ...)
    fit$centers  <- rbind(fit$centers, c(min(xy[, 1L]) - 1e4, min(xy[, 2L]) - 1e4))
    fit$size     <- c(fit$size, 0L)
    fit$withinss <- c(fit$withinss, 0)
    fit
  })
  w <- tryCatch(voronoi_seeds_kmeans(pts, 6, partition = "density"), warning = identity)
  expect_s3_class(w, "spatialkit_empty_cells")
  expect_identical(w$empty, 1L)
  expect_match(conditionMessage(w), "1 of the 6 seeds hold no point", fixed = TRUE)
  expect_warning(g <- get_voronoi_seeds(method = "kmeans", n = 6, sample_points = pts,
                                        partition = "density", set_seed = 1),
                 class = "spatialkit_empty_cells")
  expect_identical(nrow(g), 6L)
  expect_identical(attr(g, "kmeans")$size[6L], 0L)
})


# ---- A centre with no point, with nothing stood in for -----------------------------

# At counts from about a fifth of the points up a start of the density kind,
# a random one most often, can end with a centre that holds no point (see
# "Conditions" in the help of resolution_profile()).  800 points with 28
# dense places and 400 centres from single random starts: measured under
# seeds 1 to 30, two runs with one such centre (seeds 7 and 23).  The seeds
# are searched for and not written down, so that the tests do not hang on
# one platform's arithmetic.
.pl_empty_layer <- function() {
  xy <- partition_towns_xy(n = 800L, seed = 5L)
  cleanup <- spatialkit:::.with_seed(6L)
  on.exit(cleanup(), add = TRUE)
  list(xy = xy, pts = partition_points_sf(xy, z = stats::rnorm(nrow(xy))))
}

test_that("a real density run that leaves a centre empty is reported and still tessellates", {
  lay <- .pl_empty_layer()
  xy <- lay$xy; pts <- lay$pts
  k <- 400L
  runs <- lapply(1:30, function(s) partition_warnings(
    voronoi_seeds_kmeans(pts, k, partition = "density", init = "random", nstart = 1, set_seed = s)))
  hit <- which(vapply(runs, function(r)
    length(partition_of_class(r, "spatialkit_empty_cells")) > 0L, logical(1)))
  expect_gt(length(hit), 0L)
  skip_if(length(hit) == 0L, "no start of the 30 left a centre empty here")
  s <- hit[1L]
  r <- runs[[s]]
  # Nothing else is raised, and the warning counts the seeds that are no
  # point's nearest.
  expect_length(r$warnings, 1L)
  w <- r$warnings[[1L]]
  seeds <- r$value
  expect_identical(nrow(seeds), k)
  expect_true(all(is.finite(partition_xy(seeds))))
  sizes <- tabulate(partition_nearest(xy, partition_xy(seeds)), nbins = k)
  expect_gte(w$empty, 1L)
  expect_identical(w$empty, sum(sizes == 0L))
  expect_identical(w$level, k)
  expect_match(conditionMessage(w),
               sprintf("voronoi_seeds_kmeans(): %d of the 400 seeds hold no point", w$empty),
               fixed = TRUE)
  expect_match(conditionMessage(w),
               sprintf("so %d of their Voronoi cells will hold a point", k - w$empty),
               fixed = TRUE)
  # It is the engine's partition under that seed, empty centre and all.
  set.seed(s)
  eng <- spatialkit:::.partition_centres(xy, k, partition = "density", init = "random",
                                         nstart = 1L, iter.max = 100L)
  expect_identical(partition_xy(seeds), eng$centers)
  expect_identical(eng$size, sizes)
  expect_true(eng$converged)

  # The tessellation has a cell for every seed.  Every observation is in
  # exactly one, and the cells that hold one are the occupied centres'.
  tess <- build_tessellation(seeds, method = "voronoi", quiet = TRUE)
  expect_identical(nrow(tess$cells), k)
  expect_identical(tess$params$partition, attr(seeds, "partition"))
  assigned <- assign_features_to_polygons(pts, tess$cells, polygon_id_col = "cell_id",
                                          quiet = TRUE)
  expect_identical(nrow(assigned), nrow(pts))
  expect_false(anyNA(assigned$cell_id))
  expect_true(all(lengths(sf::st_intersects(pts, tess$cells)) == 1L))
  expect_identical(length(unique(assigned$cell_id)), k - w$empty)
  expect_identical(sort(as.integer(table(assigned$cell_id))), sort(sizes[sizes > 0L]))
  summ <- summarize_by_cell(assigned, response_var = "z", id_col = "cell_id", estimate = "mean")
  expect_identical(sum(summ$n), nrow(pts))
  expect_identical(nrow(summ), k - w$empty)

  # get_voronoi_seeds() under the same seed: the same seeds, the same
  # warning, and a size of 0 in the record of the run.
  rg <- partition_warnings(get_voronoi_seeds(method = "kmeans", n = k, sample_points = pts,
                                             partition = "density", init = "random",
                                             kmeans_nstart = 1, set_seed = s))
  expect_length(rg$warnings, 1L)
  expect_s3_class(rg$warnings[[1L]], "spatialkit_empty_cells")
  expect_identical(rg$warnings[[1L]]$empty, w$empty)
  expect_identical(partition_xy(rg$value), partition_xy(seeds))
  expect_identical(attr(rg$value, "kmeans")$size, sizes)
})

test_that("a real profile level with an empty centre is reported, and its seeds carry it", {
  lay <- .pl_empty_layer()
  xy <- lay$xy; pts <- lay$pts
  # Measured under seeds 1 to 12: one empty centre at 400 cells under seeds
  # 2 and 10.
  runs <- lapply(1:12, function(s) partition_warnings(
    resolution_profile(pts, levels = c(200L, 400L), nstart = 1, seed = s, min_cell_n = 1,
                       partition = "density", init = "random")))
  hit <- which(vapply(runs, function(r)
    length(partition_of_class(r, "spatialkit_empty_cells")) > 0L, logical(1)))
  expect_gt(length(hit), 0L)
  skip_if(length(hit) == 0L, "no profile of the 12 had a level with an empty centre here")
  r <- runs[[hit[1L]]]
  prof <- r$value
  for (w in partition_of_class(r, "spatialkit_empty_cells")) {
    row <- prof[prof$levels == w$level, ]
    cen <- attr(prof, "centres")[[as.character(w$level)]]
    expect_identical(nrow(cen), w$level)
    sizes <- tabulate(partition_nearest(xy, cen), nbins = w$level)
    expect_gte(w$empty, 1L)
    expect_identical(w$empty, sum(sizes == 0L))
    expect_identical(row$cell_n_min, 0L)
    expect_match(conditionMessage(w),
                 sprintf("resolution_profile(): at %d cells, %d of the centres hold no point",
                         w$level, w$empty), fixed = TRUE)
    # The level's seeds are all its centres; their cells hold every point,
    # in as many cells as the level has occupied ones.
    seeds <- voronoi_seeds_kmeans(pts, structure(
      list(best = w$level, seeds = spatialkit:::.profile_seeds(prof, w$level), criterion = "cp",
           edge = NA_character_, structured = TRUE), class = "resolution_selection"))
    expect_identical(nrow(seeds), w$level)
    tess <- build_tessellation(seeds, method = "voronoi", quiet = TRUE)
    assigned <- assign_features_to_polygons(pts, tess$cells, polygon_id_col = "cell_id",
                                            quiet = TRUE)
    expect_false(anyNA(assigned$cell_id))
    expect_identical(length(unique(assigned$cell_id)), w$level - w$empty)
  }
  # A level with no empty centre has no such warning and a smallest cell of
  # at least one point.
  quiet_levels <- setdiff(prof$levels, vapply(partition_of_class(r, "spatialkit_empty_cells"),
                                              function(w) w$level, integer(1)))
  expect_true(all(prof$cell_n_min[prof$levels %in% quiet_levels] >= 1L))
})
