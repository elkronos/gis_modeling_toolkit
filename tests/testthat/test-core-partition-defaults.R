# tests/testthat/test-core-partition-defaults.R
# ---------------------------------------------------------------------------
# The partition arguments leave every default result as it was.
#
# resolution_profile() has always scored the best of `nstart` k-means++
# starts; get_voronoi_seeds(method = "kmeans") and voronoi_seeds_kmeans()
# have always run stats::kmeans() from random starts.  Naming a function's
# own kind and starting rule returns the object a call without the arguments
# returns, and that object is pinned here to the computation it has always
# been: nothing is added to it, not even a record of the kind.
# ---------------------------------------------------------------------------

.pd_pts <- function(n = 240L, seed = 2L) {
  xy <- partition_clustered_xy(n_bg = n %/% 3L, n_cl = n - n %/% 3L, clusters = 3L,
                               side = 2000, origin = c(5e5, 5e6), seed = seed)
  cleanup <- spatialkit:::.with_seed(seed + 100L)
  on.exit(cleanup(), add = TRUE)
  z <- sin(xy[, 1L] / 300) + cos(xy[, 2L] / 400) + stats::rnorm(nrow(xy), sd = 0.5)
  partition_points_sf(xy, z = z)
}

# .kmeans_best() as it stood before the partition engine: the best of
# `nstart` runs of stats::kmeans(), each from k-means++ seeds.
.pd_old_best <- function(xy, k, nstart = 25L, iter.max = 50L) {
  best <- NULL
  for (r in seq_len(nstart)) {
    km <- try(suppressWarnings(stats::kmeans(xy, centers = spatialkit:::.kmeanspp_centers(xy, k),
                                             iter.max = iter.max, nstart = 1L)),
              silent = TRUE)
    if (inherits(km, "try-error")) next
    if (is.null(best) || km$tot.withinss < best$tot.withinss) best <- km
  }
  best
}


test_that("the engine's k-means kind from k-means++ starts is .kmeans_best() bit for bit", {
  xy <- partition_sorted(partition_clustered_xy(n_bg = 100L, n_cl = 200L, seed = 5L))
  for (k in c(2L, 9L, 40L)) {
    set.seed(31); old <- .pd_old_best(xy, k, nstart = 6L)
    set.seed(31); kb  <- spatialkit:::.kmeans_best(xy, k, nstart = 6L)
    set.seed(31); eng <- spatialkit:::.partition_centres(xy, k, partition = "kmeans",
                                                        init = "kmeans++", nstart = 6L)
    expect_identical(kb$km, old)
    expect_identical(eng$km, old)
    expect_identical(eng$centers, unname(old$centers))
    expect_identical(eng$cluster, old$cluster)
    expect_identical(eng$wss, old$tot.withinss)
    expect_identical(eng$objective, old$tot.withinss)
    expect_identical(eng$withinss, as.numeric(old$withinss))
    expect_identical(eng$spread, kb$spread)
  }
})

test_that("resolution_profile() with its own kind named is the profile without the arguments", {
  pts <- .pd_pts()
  plain <- partition_quietly(resolution_profile(pts, "z", n_levels = 6, nstart = 5))
  named <- partition_quietly(resolution_profile(pts, "z", n_levels = 6, nstart = 5,
                                                partition = "kmeans", init = "kmeans++"))
  expect_identical(named, plain)
  # The density arguments at their defaults change nothing either.
  expect_identical(partition_quietly(resolution_profile(pts, "z", n_levels = 6, nstart = 5,
                                                        density_power = 1, density_k = 8L)),
                   plain)
  # No record of the kind: the profile is the object it always was.
  expect_null(attr(plain, "partition"))
  expect_identical(capture.output(print(named)), capture.output(print(plain)))
  expect_match(capture.output(print(plain)), "5 k-means++ restarts per level",
               fixed = TRUE, all = FALSE)

  geo  <- resolution_profile(pts, n_levels = 5, nstart = 4)
  expect_identical(resolution_profile(pts, n_levels = 5, nstart = 4,
                                      partition = "kmeans", init = "kmeans++"), geo)
  expect_null(attr(geo, "partition"))
  # `iter_max = NULL` is the kind's own limit, 50 for this one; neither it
  # nor the limit named leaves a record on the default kind.
  expect_identical(resolution_profile(pts, n_levels = 5, nstart = 4, iter_max = NULL), geo)
  expect_identical(resolution_profile(pts, n_levels = 5, nstart = 4, iter_max = 50), geo)
  # The default vector of `init`, passed as it is written, is the default.
  expect_identical(resolution_profile(pts, n_levels = 5, nstart = 4,
                                      init = c("kmeans++", "random")), geo)
})

test_that("the default profile's centres are those of the k-means++ run it always made", {
  # Geometry only, so nothing draws from the stream before the sweep, and the
  # first level's centres are the best of `nstart` k-means++ starts on the
  # rows in coordinate order under the profile's seed.
  pts <- .pd_pts(n = 150L, seed = 4L)
  xs  <- partition_sorted(partition_xy(pts))
  prof <- resolution_profile(pts, levels = 7L, nstart = 6, seed = 19)
  set.seed(19); old <- .pd_old_best(xs, 7L, nstart = 6L)
  expect_identical(attr(prof, "centres")[["7"]], unname(old$centers))
  expect_identical(prof$wss, old$tot.withinss)
})

test_that("select_resolution() and the scored seeds of a default profile carry no record", {
  pts  <- .pd_pts()
  prof <- partition_quietly(resolution_profile(pts, "z", n_levels = 6, nstart = 5))
  sel  <- select_resolution(prof)
  expect_false("partition" %in% names(sel))
  expect_null(attr(sel$seeds, "partition"))
  s_plain <- partition_quietly(voronoi_seeds_kmeans(pts, sel))
  s_named <- partition_quietly(voronoi_seeds_kmeans(pts, sel, partition = "kmeans"))
  expect_identical(s_named, s_plain)
  expect_null(attr(s_plain, "partition"))
  g_plain <- partition_quietly(get_voronoi_seeds(method = "kmeans", n = sel, sample_points = pts))
  g_named <- partition_quietly(get_voronoi_seeds(method = "kmeans", n = sel, sample_points = pts,
                                                 partition = "kmeans", init = "kmeans++"))
  expect_identical(g_named, g_plain)
  expect_null(attr(g_plain, "partition"))
  # Naming the kind the profile was built with is not a mismatch.
  expect_no_warning(withCallingHandlers(
    get_voronoi_seeds(method = "kmeans", n = sel, sample_points = pts, partition = "kmeans"),
    spatialkit_bound_chosen = function(w) invokeRestart("muffleWarning"),
    spatialkit_no_structure = function(w) invokeRestart("muffleWarning")))
})

test_that("voronoi_seeds_kmeans() with its own kind named is the call without the arguments", {
  pts <- .pd_pts()
  plain <- voronoi_seeds_kmeans(pts, 11)
  expect_identical(voronoi_seeds_kmeans(pts, 11, partition = "kmeans"), plain)
  expect_identical(voronoi_seeds_kmeans(pts, 11, partition = "kmeans", init = "random"), plain)
  expect_identical(voronoi_seeds_kmeans(pts, 11, init = NULL, density_power = 1, density_k = 8),
                   plain)
  expect_null(attr(plain, "partition"))
  # `iter_max = NULL` is the kind's own limit, 50 for this one.
  expect_identical(voronoi_seeds_kmeans(pts, 11, iter_max = NULL), plain)
  expect_identical(voronoi_seeds_kmeans(pts, 11, iter_max = 50), plain)
  # The run it has always made: stats::kmeans() from 10 random starts under
  # seed 456, on the rows as they come.
  set.seed(456)
  km <- stats::kmeans(sf::st_coordinates(pts)[, 1:2], centers = 11, iter.max = 50, nstart = 10L)
  expect_identical(partition_xy(plain), unname(km$centers))
  # Other seeds and start counts as well.
  set.seed(7)
  km3 <- stats::kmeans(sf::st_coordinates(pts)[, 1:2], centers = 30, iter.max = 50, nstart = 3L)
  expect_identical(partition_xy(voronoi_seeds_kmeans(pts, 30, set_seed = 7, nstart = 3,
                                                     partition = "kmeans")),
                   unname(km3$centers))
})

test_that("for the seeding functions k-means++ starts are a new run, and say so", {
  # Their own starting rule is random, so `init = "kmeans++"` is not the
  # default under another name: it gives other seeds and records the rule.
  pts <- .pd_pts()
  plain <- voronoi_seeds_kmeans(pts, 11)
  pp    <- voronoi_seeds_kmeans(pts, 11, init = "kmeans++")
  expect_false(isTRUE(all.equal(partition_xy(pp), partition_xy(plain))))
  # The best of `nstart` k-means++ starts on the rows in coordinate order.
  xs <- partition_sorted(sf::st_coordinates(pts)[, 1:2])
  set.seed(456); old <- .pd_old_best(xs, 11L, nstart = 10L, iter.max = 50)
  expect_identical(partition_xy(pp), unname(old$centers))
  # The record: the kind and the rule, and how the kept start ended.
  expect_identical(old$ifault, 0L)
  expect_identical(attr(pp, "partition"),
                   list(kind = "kmeans", init = "kmeans++",
                        density_power = NA_real_, density_k = NA_integer_,
                        iter = as.integer(old$iter), iter_max = 50L, converged = TRUE))
})

test_that("get_voronoi_seeds() with its own kind named is the call without the arguments", {
  pts <- .pd_pts()
  plain <- get_voronoi_seeds(method = "kmeans", n = 9, sample_points = pts, set_seed = 3)
  expect_identical(get_voronoi_seeds(method = "kmeans", n = 9, sample_points = pts, set_seed = 3,
                                     partition = "kmeans"), plain)
  expect_identical(get_voronoi_seeds(method = "kmeans", n = 9, sample_points = pts, set_seed = 3,
                                     partition = "kmeans", init = "random"), plain)
  expect_null(attr(plain, "partition"))
  expect_named(attr(plain, "kmeans"),
               c("cluster", "rows", "size", "withinss", "tot_withinss", "iter", "nstart"))
  set.seed(3)
  km <- stats::kmeans(x = sf::st_coordinates(pts)[, 1:2], centers = 9, iter.max = 100,
                      nstart = 10)
  expect_identical(partition_xy(plain), unname(km$centers))
  expect_identical(attr(plain, "kmeans")$cluster, as.integer(km$cluster))
  expect_identical(attr(plain, "kmeans")$withinss, as.numeric(km$withinss))
  expect_identical(attr(plain, "kmeans")$tot_withinss, as.numeric(km$tot.withinss))
  expect_identical(attr(plain, "kmeans")$iter, as.integer(km$iter))

  # A cloud drawn inside a boundary, and the other two methods.
  bnd <- sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(pts)))
  expect_identical(get_voronoi_seeds(bnd, method = "kmeans", n = 5, set_seed = 8,
                                     partition = "kmeans"),
                   get_voronoi_seeds(bnd, method = "kmeans", n = 5, set_seed = 8))
  expect_identical(get_voronoi_seeds(bnd, method = "random", n = 5, set_seed = 8,
                                     partition = "kmeans"),
                   get_voronoi_seeds(bnd, method = "random", n = 5, set_seed = 8))
  expect_identical(get_voronoi_seeds(method = "provided", seeds = pts[1:4, ],
                                     partition = "kmeans"),
                   get_voronoi_seeds(method = "provided", seeds = pts[1:4, ]))
})

test_that("a tessellation of default seeds has the params it always had", {
  pts  <- .pd_pts()
  tess <- build_tessellation(voronoi_seeds_kmeans(pts, 11), method = "voronoi", quiet = TRUE)
  expect_named(tess$params, c("clip", "expand", "keep_duplicates", "snapped", "n_outside",
                              "boundary_source", "cover_outside"))
  cvp <- create_voronoi_polygons(voronoi_seeds_kmeans(pts, 11), quiet = TRUE)
  expect_false("partition" %in% names(cvp$params))
  # The reader of the record says plain k-means of the default seeds.  Their
  # cells carry nothing that tells them from the Voronoi cells of any other
  # points, so of the tessellation it says nothing.
  expect_identical(spatialkit:::.partition_record(voronoi_seeds_kmeans(pts, 11)),
                   list(kind = "kmeans", init = NA_character_,
                        density_power = NA_real_, density_k = NA_integer_))
  expect_identical(spatialkit:::.partition_record(tess),
                   list(kind = NA_character_, init = NA_character_,
                        density_power = NA_real_, density_k = NA_integer_))
})

test_that("the reader of the record gives no kind to what k-means did not build", {
  pts <- .pd_pts()
  bnd <- sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(pts)))
  none <- list(kind = NA_character_, init = NA_character_,
               density_power = NA_real_, density_k = NA_integer_)
  rec <- spatialkit:::.partition_record
  expect_identical(rec(get_voronoi_seeds(bnd, method = "random", n = 5, set_seed = 1)), none)
  expect_identical(rec(voronoi_seeds_random(bnd, 5, set_seed = 1)), none)
  expect_identical(rec(get_voronoi_seeds(method = "provided", seeds = pts[1:4, ])), none)
  expect_identical(rec(pts), none)
  expect_identical(rec(partition_quietly(build_tessellation(pts, boundary = bnd, method = "hex",
                                                            approx_n_cells = 20, quiet = TRUE))),
                   none)
  expect_identical(rec(partition_quietly(build_tessellation(pts, method = "voronoi",
                                                            quiet = TRUE))), none)
  expect_identical(rec(NULL), none)
  expect_identical(rec(list(a = 1)), none)
  # What k-means built with no record is plain k-means: seeds, a profile and
  # a selection.
  prof <- partition_quietly(resolution_profile(pts, "z", n_levels = 4, nstart = 2))
  expect_identical(rec(prof)$kind, "kmeans")
  expect_identical(rec(select_resolution(prof))$kind, "kmeans")
  expect_identical(rec(get_voronoi_seeds(method = "kmeans", n = 5, sample_points = pts,
                                         set_seed = 1))$kind, "kmeans")
})

test_that("determine_optimal_levels() takes no partition argument and keeps plain k-means", {
  expect_false(any(c("partition", "init", "density_power", "density_k") %in%
                     names(formals(determine_optimal_levels))))
  # Its sweep is the best of 25 k-means++ starts at every k, one after the
  # other under its seed, on the rows in coordinate order: the WSS curve it
  # reports is that computation.
  pts <- .pd_pts(n = 300L, seed = 6L)
  out <- partition_quietly(determine_optimal_levels(pts, max_levels = 40L, set_seed = 5,
                                                    response_var = "z",
                                                    criterion = "morans_i"))
  diag <- attr(out, "diagnostics")
  expect_false(is.null(diag))
  xy <- partition_xy(pts)
  xs <- xy[order(xy[, 1L], xy[, 2L], pts$z), , drop = FALSE]
  set.seed(5)
  wss <- c(sum(sweep(xs, 2L, colMeans(xs))^2),
           vapply(2:length(diag$wss), function(k) .pd_old_best(xs, k)$tot.withinss, numeric(1)))
  expect_equal(diag$wss, wss, tolerance = 1e-12)
  expect_identical(diag$wss[-1L], wss[-1L])
})

test_that("the argument order of the three functions is unchanged ahead of the new ones", {
  new_args <- c("partition", "init", "density_power", "density_k")
  # The five partition arguments of resolution_profile(), found by name: they
  # stand together in this order, straight after the arguments the function
  # had before them, and `estimate` and `obs_per_cell` follow `iter_max`.
  prof_args <- names(formals(resolution_profile))
  at <- match(c(new_args, "iter_max"), prof_args)
  expect_identical(at, at[1L] + 0:4)
  expect_identical(prof_args[at[1L] - 1L], "range_floor")
  expect_identical(prof_args[at[5L] + 1:2], c("estimate", "obs_per_cell"))
  expect_identical(utils::tail(names(formals(get_voronoi_seeds)), 4L), new_args)
  expect_identical(utils::tail(names(formals(voronoi_seeds_kmeans)), 5L),
                   c(new_args, "iter_max"))
  expect_null(formals(resolution_profile)$iter_max)
  expect_null(formals(voronoi_seeds_kmeans)$iter_max)
  expect_identical(names(formals(voronoi_seeds_kmeans))[1:4],
                   c("points_sf", "k", "set_seed", "nstart"))
  expect_identical(eval(formals(resolution_profile)$init), c("kmeans++", "random"))
  expect_identical(eval(formals(resolution_profile)$partition), c("kmeans", "density"))
  expect_null(formals(voronoi_seeds_kmeans)$init)
  expect_null(formals(get_voronoi_seeds)$init)
})
