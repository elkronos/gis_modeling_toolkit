# tests/testthat/test-core-partition-arguments.R
# ---------------------------------------------------------------------------
# The partition arguments reach the engine.  For each public function the
# centres at a power, a neighbour count and a starting rule other than the
# default are the engine's at those settings under the same seed, and are
# not the default's: an argument that is recorded and then dropped on the
# way down would pass every check of the record.
# Also here: what the functions say of an argument that nothing reads.
# ---------------------------------------------------------------------------

.pa_xy <- function() partition_clustered_xy(n_bg = 200L, n_cl = 500L, seed = 17L)

# The engine at the settings a public function is given, under its seed.
.pa_engine <- function(xy, k, seed, nstart, settings) {
  set.seed(seed)
  do.call(spatialkit:::.partition_centres,
          c(list(xy, k, partition = "density", nstart = nstart, iter.max = 100L), settings))
}

# Settings other than the defaults (power 1, 8 neighbours, k-means++ starts),
# one at a time and together.
.pa_settings <- list(
  power  = list(density_power = 0.4),
  k_nn   = list(density_k = 3),
  random = list(init = "random"),
  all    = list(density_power = 1.6, density_k = 20, init = "random"))

# What the record of a call at those settings says.
.pa_record <- function(settings) {
  rec <- list(kind = "density", init = "kmeans++", density_power = 1, density_k = 8L)
  if (!is.null(settings$init)) rec$init <- settings$init
  if (!is.null(settings$density_power)) rec$density_power <- settings$density_power
  if (!is.null(settings$density_k)) rec$density_k <- as.integer(settings$density_k)
  rec
}


test_that("the engine's settings change its centres", {
  # The premise of the three tests below: each setting gives another
  # partition than the default's, by far more than rounding.
  xy <- .pa_xy()
  base <- .pa_engine(xy, 25L, 6L, 3L, list())
  for (nm in names(.pa_settings)) {
    other <- .pa_engine(xy, 25L, 6L, 3L, .pa_settings[[nm]])
    expect_gt(max(abs(other$centers - base$centers)), 1, label = nm)
    expect_false(identical(other$cluster, base$cluster), label = nm)
  }
  # The power and the neighbour count act through the weights.
  expect_false(isTRUE(all.equal(.pa_engine(xy, 25L, 6L, 3L, list(density_power = 0.4))$weights,
                                base$weights)))
  expect_false(isTRUE(all.equal(.pa_engine(xy, 25L, 6L, 3L, list(density_k = 3))$weights,
                                base$weights)))
})

test_that("a random start of the density kind is k locations drawn in coordinate order", {
  # One start, one step of Lloyd's iteration: the centres are then the
  # weighted means of the Voronoi cells of the start, which pins the start to
  # sample.int() over the distinct locations in coordinate order.
  xy <- .pa_xy()
  xs <- partition_sorted(xy)
  for (k in c(4L, 30L)) for (seed in c(2L, 9L)) {
    set.seed(seed)
    got <- spatialkit:::.partition_centres(xy, k, partition = "density", init = "random",
                                           nstart = 1L, iter.max = 1L)
    set.seed(seed)
    start <- xs[sample.int(nrow(xs), k), , drop = FALSE]
    cl <- partition_nearest(xs, start)
    w  <- got$weights[order(xy[, 1L], xy[, 2L])]
    ref <- rowsum(xs * w, cl) / as.numeric(rowsum(w, cl))
    expect_identical(sort(unique(cl)), seq_len(k))       # every start location has its own point
    expect_equal(got$centers, ref, tolerance = 1e-10, ignore_attr = TRUE)
    # A k-means++ start under the same seed is another start.
    set.seed(seed)
    pp <- spatialkit:::.partition_centres(xy, k, partition = "density", init = "kmeans++",
                                          nstart = 1L, iter.max = 1L)
    expect_gt(max(abs(pp$centers - got$centers)), 1)
  }
})

test_that("voronoi_seeds_kmeans() hands density_power, density_k and init to the engine", {
  xy  <- .pa_xy()
  pts <- partition_points_sf(xy)
  k <- 25L
  base <- voronoi_seeds_kmeans(pts, k, partition = "density", nstart = 3, set_seed = 6)
  expect_identical(partition_xy(base), .pa_engine(xy, k, 6L, 3L, list())$centers)
  for (nm in names(.pa_settings)) {
    st <- .pa_settings[[nm]]
    s <- do.call(voronoi_seeds_kmeans,
                 c(list(pts, k, partition = "density", nstart = 3, set_seed = 6), st))
    expect_identical(partition_xy(s), .pa_engine(xy, k, 6L, 3L, st)$centers, label = nm)
    expect_gt(max(abs(partition_xy(s) - partition_xy(base))), 1, label = nm)
    expect_identical(attr(s, "partition")[1:4], .pa_record(st), label = nm)
  }
  # The default number of starts is 10.
  expect_identical(partition_xy(voronoi_seeds_kmeans(pts, k, partition = "density", set_seed = 6)),
                   .pa_engine(xy, k, 6L, 10L, list())$centers)
})

test_that("get_voronoi_seeds() hands density_power, density_k and init to the engine", {
  xy  <- .pa_xy()
  pts <- partition_points_sf(xy)
  k <- 25L
  call <- function(...) get_voronoi_seeds(method = "kmeans", n = k, sample_points = pts,
                                          partition = "density", kmeans_nstart = 3,
                                          set_seed = 6, ...)
  base <- call()
  expect_identical(partition_xy(base), .pa_engine(xy, k, 6L, 3L, list())$centers)
  for (nm in names(.pa_settings)) {
    st <- .pa_settings[[nm]]
    s <- do.call(call, st)
    eng <- .pa_engine(xy, k, 6L, 3L, st)
    expect_identical(partition_xy(s), eng$centers, label = nm)
    expect_gt(max(abs(partition_xy(s) - partition_xy(base))), 1, label = nm)
    expect_identical(attr(s, "partition")[1:4], .pa_record(st), label = nm)
    # The record of the run is the engine's at those settings as well.
    expect_identical(attr(s, "kmeans")$cluster, eng$cluster, label = nm)
    expect_identical(attr(s, "kmeans")$objective, eng$objective, label = nm)
    expect_identical(attr(s, "kmeans")$iter, eng$iter, label = nm)
  }
  # The default number of starts is 10.
  expect_identical(partition_xy(get_voronoi_seeds(method = "kmeans", n = k, sample_points = pts,
                                                  partition = "density", set_seed = 6)),
                   .pa_engine(xy, k, 6L, 10L, list())$centers)
})

test_that("resolution_profile() hands density_power, density_k and init to the engine", {
  # Geometry only, so nothing draws from the stream before the sweep.
  xy  <- .pa_xy()
  pts <- partition_points_sf(xy)
  k <- 25L
  call <- function(...) resolution_profile(pts, levels = k, nstart = 3, seed = 6,
                                           partition = "density", ...)
  base <- call()
  expect_identical(attr(base, "centres")[[as.character(k)]],
                   .pa_engine(xy, k, 6L, 3L, list())$centers)
  for (nm in names(.pa_settings)) {
    st <- .pa_settings[[nm]]
    p <- do.call(call, st)
    eng <- .pa_engine(xy, k, 6L, 3L, st)
    cen <- attr(p, "centres")[[as.character(k)]]
    expect_identical(cen, eng$centers, label = nm)
    expect_gt(max(abs(cen - attr(base, "centres")[[as.character(k)]])), 1, label = nm)
    expect_identical(attr(p, "partition")[1:4], .pa_record(st), label = nm)
    expect_identical(p$wss, eng$wss, label = nm)
    expect_identical(unname(attr(p, "partition")$objective), eng$objective, label = nm)
    expect_identical(p$cell_n_min, min(eng$size), label = nm)
    # The selection's seeds are those centres, with the record.
    seeds <- spatialkit:::.profile_seeds(p, k)
    expect_identical(partition_xy(seeds), eng$centers, label = nm)
    expect_identical(attr(seeds, "partition"), .pa_record(st), label = nm)
  }
  # The second level of a ladder is built with the same settings: its starts
  # follow the first level's in the stream.
  st <- .pa_settings$all
  p2 <- do.call(resolution_profile, c(list(pts, levels = c(8L, k), nstart = 2, seed = 6,
                                           partition = "density"), st))
  set.seed(6)
  e1 <- do.call(spatialkit:::.partition_centres,
                c(list(xy, 8L, partition = "density", nstart = 2L, iter.max = 100L), st))
  e2 <- do.call(spatialkit:::.partition_centres,
                c(list(xy, k, partition = "density", nstart = 2L, iter.max = 100L), st))
  expect_identical(attr(p2, "centres")[["8"]], e1$centers)
  expect_identical(attr(p2, "centres")[[as.character(k)]], e2$centers)
})


# ---- Arguments nothing reads -----------------------------------------------------

.pa_scored <- function() {
  xy <- .pa_xy()
  cleanup <- spatialkit:::.with_seed(18L)
  z <- sin(xy[, 1L] / 150) + stats::rnorm(nrow(xy), sd = 0.5)
  cleanup()
  pts <- partition_points_sf(xy, z = z)
  sel_of <- function(p, L) structure(
    list(best = L, seeds = spatialkit:::.profile_seeds(p, L), criterion = "cp",
         edge = NA_character_, structured = TRUE), class = "resolution_selection")
  p_km  <- partition_quietly(resolution_profile(pts, "z", levels = c(6L, 14L), nstart = 2))
  p_den <- partition_quietly(resolution_profile(pts, "z", levels = c(6L, 14L), nstart = 2,
                                                partition = "density", density_power = 0.5))
  list(pts = pts, p_km = p_km, p_den = p_den, km = sel_of(p_km, 14L), den = sel_of(p_den, 14L))
}

test_that("a density setting named beside scored density centres is not called ignored", {
  sc <- .pa_scored()
  rec <- list(kind = "density", init = "kmeans++", density_power = 0.5, density_k = 8L)
  # The power the profile was built with, named with no `partition`: the
  # centres returned are density ones at that power, and nothing is wrong.
  for (f in list(function(...) voronoi_seeds_kmeans(sc$pts, sc$den, ...),
                 function(...) get_voronoi_seeds(method = "kmeans", n = sc$den,
                                                 sample_points = sc$pts, ...))) {
    r <- partition_warnings(f(density_power = 0.5))
    expect_length(r$warnings, 0L)
    expect_identical(attr(r$value, "partition"), rec)
    expect_equal(partition_xy(r$value), attr(sc$p_den, "centres")[["14"]], tolerance = 1e-12)
    # Another power, or another neighbour count, than the profile's: the
    # scored centres are returned and the call is told its setting is unused.
    for (other in list(list(density_power = 0.7), list(density_k = 4))) {
      r <- partition_warnings(do.call(f, other))
      expect_length(r$warnings, 1L)
      expect_s3_class(r$warnings[[1L]], "spatialkit_partition_unused")
      expect_match(conditionMessage(r$warnings[[1L]]),
                   "built with partition = \"density\" from k-means++ starts (power 0.5, 8 neighbours)",
                   fixed = TRUE)
      expect_identical(attr(r$value, "partition"), rec)
    }
    # With `partition = "density"` named, the defaults of the two settings
    # are named with it, and the power is not the profile's.
    r <- partition_warnings(f(partition = "density"))
    expect_length(partition_of_class(r, "spatialkit_partition_unused"), 1L)
    expect_length(partition_warnings(f(partition = "density", density_power = 0.5))$warnings, 0L)
  }
})

test_that("a density setting named beside k-means centres, scored or fresh, is said to be ignored", {
  sc <- .pa_scored()
  text <- "`density_power` is read by partition = \"density\" only, and partition = \"kmeans\" ignores it"
  r <- partition_warnings(voronoi_seeds_kmeans(sc$pts, sc$km, density_power = 0.7))
  expect_length(r$warnings, 1L)
  expect_match(conditionMessage(r$warnings[[1L]]), paste0("voronoi_seeds_kmeans(): ", text),
               fixed = TRUE)
  expect_null(attr(r$value, "partition"))
  r <- partition_warnings(get_voronoi_seeds(method = "kmeans", n = sc$km, sample_points = sc$pts,
                                            density_power = 0.7))
  expect_length(r$warnings, 1L)
  expect_match(conditionMessage(r$warnings[[1L]]), paste0("get_voronoi_seeds(): ", text),
               fixed = TRUE)
  # A fresh run: once, and the seeds are the default's.
  r <- partition_warnings(voronoi_seeds_kmeans(sc$pts, 9, density_power = 0.7))
  expect_length(r$warnings, 1L)
  expect_match(conditionMessage(r$warnings[[1L]]), text, fixed = TRUE)
  expect_identical(r$value, voronoi_seeds_kmeans(sc$pts, 9))
})

test_that("get_voronoi_seeds() checks its two counts before a new kind reads them", {
  pts <- partition_points_sf(.pa_xy())
  for (extra in list(list(partition = "density"), list(init = "kmeans++"))) {
    f <- function(...) do.call(get_voronoi_seeds,
                               c(list(method = "kmeans", n = 5, sample_points = pts, set_seed = 1),
                                 extra, list(...)))
    for (bad in list(0, NA, "a", 2.5, c(1, 2), Inf, -3)) {
      expect_error(f(kmeans_nstart = bad),
                   "get_voronoi_seeds(): `kmeans_nstart` must be NULL or a single positive whole number of starts",
                   fixed = TRUE)
      expect_error(f(kmeans_iter = bad),
                   "get_voronoi_seeds(): `kmeans_iter` must be a single positive whole number of iterations",
                   fixed = TRUE)
    }
    # What is recorded is what ran.
    g <- f(kmeans_nstart = 2, kmeans_iter = 60)
    expect_identical(attr(g, "kmeans")$nstart, 2L)
    expect_identical(attr(g, "partition")$iter_max, 60L)
  }
  # The run the function has always made hands both to stats::kmeans() as
  # they are, as it always did: a start count of 2.5 is that function's two
  # starts, and is not refused here.
  set.seed(1)
  km <- stats::kmeans(x = partition_xy(pts), centers = 5, iter.max = 100, nstart = 2.5)
  g <- get_voronoi_seeds(method = "kmeans", n = 5, sample_points = pts, set_seed = 1,
                         kmeans_nstart = 2.5)
  expect_identical(partition_xy(g), unname(km$centers))
  expect_identical(attr(g, "kmeans")$nstart, 2L)
})

test_that("summary() of a profile names the kind of its cells when the profile records one", {
  sc <- .pa_scored()
  sm <- suppressWarnings(summary(sc$p_den))
  expect_identical(attr(sm, "partition"),
                   list(kind = "density", init = "kmeans++", density_power = 0.5, density_k = 8L))
  expect_output(print(sm), "  cells: density-weighted k-means++ starts (power 0.5, 8 neighbours)",
                fixed = TRUE)
  # A k-means profile's summary is the object it always was, and prints so.
  sm_km <- suppressWarnings(summary(sc$p_km))
  expect_null(attr(sm_km, "partition"))
  expect_false(any(grepl("cells:", capture.output(print(sm_km)), fixed = TRUE)))
})
