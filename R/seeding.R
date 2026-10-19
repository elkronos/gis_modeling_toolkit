#' Generate seed points for Voronoi tessellation
#'
#' Creates an sf POINT layer of "seed" locations. Multiple strategies are
#' supported: user-provided points, uniform random sampling within a boundary,
#' or k-means clustering of a sampling cloud.
#'
#' Seeds lie inside the ground they stand for (k-means seeds inside the
#' hull of the points they were fitted to, random seeds inside the
#' boundary), so Voronoi cells clipped to the seeds' own hull would leave
#' its edge uncovered. The result therefore
#' carries that ground as `attr(, "cover")`: the `boundary` when one was
#' given, otherwise the convex hull of `sample_points`, otherwise (for the
#' scored centres of a [resolution_profile()] or a [select_resolution()]
#' result passed as `n`) the hull of the points the profile was computed on.
#' [build_tessellation()] and
#' [create_voronoi_polygons()] derive their default boundary from it, so
#' `build_tessellation(seeds, method = "voronoi")` holds every observation
#' without a `boundary`. Seeds with no cover are tessellated over their own
#' hull unless a `boundary` or `cover` is passed then: `method = "provided"`
#' seeds given with neither `boundary` nor `sample_points` have none, and no
#' warning is raised for them, since they may be the observation sites
#' themselves. The attribute survives
#' `sf::st_transform()` and row subsetting; `rbind()` drops it.
#'
#' A k-means seed is a mean of the points of its cluster, weighted under the
#' density kind, so it can lie outside a `boundary` that is not convex, in a
#' bay or a bend of a coast, even when every point lies inside. Tessellated
#' within that boundary, its cell is built if any of it lies inside, but the
#' seed itself lies in no cell: [build_tessellation()] warns that it falls
#' in no cell (class `spatialkit_unassigned`) and gives it `NA` in `index`.
#' `seeds` in [cell_values_at()] says how to give such a seed its cell's
#' ID.
#'
#' @param boundary Optional polygonal sf object defining the sampling area.
#' @param method One of "kmeans", "random", "provided".
#' @param n Integer; number of seeds to return. Required for
#'   `method = "kmeans"` and `method = "random"`. **Ignored** for
#'   `method = "provided"`, where every row of `seeds` is returned; a mismatch
#'   between `n` and `nrow(seeds)` is reported as a warning.
#'
#'   For `method = "kmeans"` it is an upper bound only:
#'   k-means cannot produce more centres than there are distinct positions in
#'   the sampling cloud, nor as many centres as there are rows. When `n`
#'   exceeds either ceiling it is clamped, with a warning naming the count
#'   actually used. `n = nrow(sample_points)` is the common case, and yields
#'   `nrow(sample_points) - 1` seeds. Check `nrow()` on the result; do not
#'   assume `n`.
#'
#'   Besides a number, `n` accepts what the level-selection step returned:
#'   the integer vector of ranked candidates from [determine_optimal_levels()]
#'   (its first element is used), a [select_resolution()] result (its
#'   `$best`), or a [resolution_profile()]. A profile built without an
#'   `estimate` is read with `select_resolution()` by `cp`. A profile built
#'   with `estimate` is read by `criterion = "auto"`, the rule that goes
#'   with its estimate: a size for pooled values and the lowest
#'   leave-one-out error for plain means. The size is the profile's own when
#'   it was built at one, and 4.5 rows per cell otherwise, and the count is
#'   the profile's level for that size, or the built level nearest to
#'   `max(2, round(n / obs_per_cell))` when the profile has none for it
#'   (with a warning of class `"spatialkit_size_not_reached"`).
#'   `select_resolution(profile)` with no `criterion` reads such a profile
#'   the same way, since its default criterion is `"auto"`. A profile built
#'   with `obs_per_cell` and without an `estimate` is read by `cp` like any
#'   profile without one; pass `select_resolution(profile, criterion =
#'   "size")` for the count its size gives. The result then
#'   carries `attr(, "n_from")`
#'   saying which. A count read off a selection or a profile comes with what
#'   the selection found: a warning of class `"spatialkit_bound_chosen"` when
#'   the optimum sits on a bound of the ladder, so the bound chose the count,
#'   and one of class `"spatialkit_no_structure"` when one cell for the whole
#'   layer scores as well. A count set by a size has no optimum, so no bound
#'   can choose it; its `"spatialkit_no_structure"` warning says that the
#'   count stands and that the cell values will say little that one mean
#'   does not. A profile without an estimate and with no `cp` to read is an
#'   error (class `"spatialkit_no_default_criterion"` when another criterion
#'   is finite); see [resolution_profile()].
#'
#'   With `method = "kmeans"`, a selection or profile also carries the
#'   centres of the partition the profile scored at that count, and those
#'   are returned instead of a new k-means: a k-means partition, weighted or
#'   not, is the Voronoi partition of its centres, so their cells are the
#'   cells the criteria judged. `sample_points`, when given, is then only
#'   assigned to them (for `attr(, "kmeans")`) and gives the cover, and
#'   `boundary` only
#'   sets the CRS and the cover (see "Details");
#'   `kmeans_nstart`, `kmeans_iter` and `set_seed` are not used, and neither
#'   are `partition`, `init` and the density settings: the centres are of
#'   the kind the profile was built with, and naming another here is a
#'   warning of class `"spatialkit_partition_unused"`. A call names what it
#'   supplies: a `partition` or an `init` passed in the call; with
#'   `partition = "density"` passed, both density settings at their values,
#'   defaults included; without it, a density setting other than its
#'   default. A call that supplies none of the four is told nothing,
#'   whichever kind and settings the profile was built with. Pass the
#'   count as a number (`n = sel$best`) for a fresh k-means instead.
#' @param seeds sf POINT object of user-provided seeds (method = "provided").
#' @param sample_points Optional sf POINT cloud for k-means clustering. Only
#'   the first two coordinate columns are clustered, so a Z or M dimension does
#'   not join the distance calculation and dominate it; rows with empty or
#'   non-finite coordinates are dropped with a warning, so they never reach
#'   `stats::kmeans()`, which fails on them without naming a cause. A lon/lat
#'   cloud, or one with no CRS whose coordinates look like lon/lat (the
#'   heuristic [ensure_projected()] applies, with its warning), is projected
#'   before clustering.
#' @param kmeans_nstart The number of starts, of which the best is kept.
#'   `NULL` (default) is the kind's own number: 3 for `partition =
#'   "density"` and 10 for `partition = "kmeans"`, the default of earlier
#'   versions. `kmeans_nstart = 10` gives the density kind the starts it had
#'   when 10 was the default. Under `partition = "kmeans"` from random starts
#'   it is `nstart` of [stats::kmeans()]. A fresh partition is not the one
#'   [resolution_profile()] scored a count on; pass the selection or profile
#'   itself as `n` for that one (see `n`). Under `partition = "density"`,
#'   and with `init = "kmeans++"`, it has to be `NULL` or a positive whole
#'   number, so a number that is not whole, which `partition = "kmeans"`
#'   takes, stops a call that names no kind. It is not read by `method =
#'   "random"` or `"provided"`.
#' @param kmeans_iter Integer; iter.max for kmeans(). Default 100. Under
#'   `partition = "density"` it is the most steps of Lloyd's iteration a
#'   start takes. Under `partition = "density"`, and with `init =
#'   "kmeans++"`, it has to be a positive whole number, and a kept start that
#'   was cut off at it is a warning of class `"spatialkit_not_converged"`
#'   (`partition = "kmeans"` from random starts keeps the warning
#'   [stats::kmeans()] gives).
#' @param set_seed Optional integer RNG seed.
#' @param partition For `method = "kmeans"`: what the seeds minimise.
#'   `"density"` (default) is the sum of squared distances from the points
#'   of the cloud to their nearest seed, each point weighted by its local
#'   density to the power `density_power`: compact cells, many small ones
#'   where the points are dense and a few large ones where they are sparse.
#'   At the default power they hold about the same number of points each
#'   once `n` is above the number of dense places (clusters) in the cloud;
#'   below that the counts are about as uneven as under `"kmeans"`.
#'   `"kmeans"` is the same sum unweighted: compact cells, with several
#'   times the points in a cell of a dense place that a cell of a sparse
#'   place holds. It is the run of earlier versions, whose default it was:
#'   [stats::kmeans()] from 10 random starts on the rows as they come. The
#'   density kind follows the
#'   density of `sample_points`; a cloud this function draws inside
#'   `boundary` is uniform and has none to follow. Without the 'FNN' package
#'   the density kind finds its weights by comparing every pair of distinct
#'   locations, so its time grows with the square of their number; the seeds
#'   are the same, and nothing stops the call. See "The partition" in
#'   [resolution_profile()] for what each kind does and how.
#' @param init For `method = "kmeans"`: how each start is drawn, `"random"`
#'   (`n` distinct locations at random, as [stats::kmeans()] draws them) or
#'   `"kmeans++"` (the seeds one at a time, each in proportion to its squared
#'   distance from those already drawn, times its weight under `partition =
#'   "density"`). `NULL` (default) is each kind's own rule: `"random"` for
#'   `partition = "kmeans"`, the run this function has always made, and
#'   `"kmeans++"` for `partition = "density"`. Under `partition = "kmeans"`
#'   from random starts the cloud is clustered in the order its rows come
#'   in; with `init = "kmeans++"`, and under `partition = "density"`, the
#'   rows are put in coordinate order first, so the seeds do not depend on
#'   their order.
#' @param density_power For `partition = "density"`: the power the local
#'   density is raised to for a point's weight, from 0 to 2. Default 1, which
#'   evens out the cell counts once the cells outnumber the dense places; 0
#'   weights every point alike. Setting it under `partition = "kmeans"` is a
#'   warning.
#' @param density_k For `partition = "density"`: the number of nearest
#'   neighbours the local density is read from. Default 8.
#' @return An sf POINT object with seed_id and method columns. With
#'   `method = "kmeans"` it also carries `attr(, "kmeans")`, the run behind
#'   the seeds: `cluster` (the `seed_id` each clustered cloud point was
#'   assigned to), `rows` (those points' positions in the cloud, since rows
#'   with unusable coordinates are dropped first), `size` (points per seed),
#'   `withinss` and `tot_withinss` (the within-cluster sums of squares),
#'   `iter` and `nstart`. For the scored centres of a selection or profile
#'   (see `n`) it describes `sample_points` assigned to them, each point to
#'   its nearest seed, with `iter` `NA` and `nstart` the profile's restarts,
#'   and is absent when no `sample_points` were given. Under `partition =
#'   "density"` `withinss` and `tot_withinss` count every point once, about
#'   its cell's mean. A run this call makes with `partition = "density"` or
#'   `init = "kmeans++"` also records `converged`, whether the kept start
#'   stopped before `kmeans_iter`, and under `partition = "density"`
#'   `objective`, the weighted sum the seeds minimise; scored centres have
#'   neither. Seeds placed with `partition = "density"` or `init =
#'   "kmeans++"`, and the scored centres of a profile of the density kind or
#'   from random starts, carry `attr(, "partition")`: a
#'   list with `kind`, `init`, `density_power` and `density_k`, and for
#'   seeds this call placed `iter`, `iter_max` and `converged` (the
#'   iterations the kept start took, the limit, and whether it stopped
#'   before the limit). [build_tessellation()] returns the list as
#'   `params$partition`. Seeds with no such
#'   record are plain k-means: from random starts when this call placed
#'   them, or the scored centres of a k-means profile from k-means++ starts.
#'   Every method also
#'   attaches `attr(, "cover")` when the call says what ground the seeds stand
#'   for (see "Details"): an sfc holding the boundary's polygons or the hull
#'   of the clustered points, in its own CRS.
#'
#'   When `n` is a [select_resolution()] result that carries `count_rule`,
#'   or a profile built with an `estimate`, the seeds carry that record as
#'   `attr(, "count_rule")`: how the count was set, for which cell value, and
#'   what the profile scored at it. Its `scored_cells` is `TRUE` for the
#'   scored centres `method = "kmeans"` hands on, and `FALSE` for seeds this
#'   call places itself at that count: a k-means of its own when the
#'   selection holds no centres at the count, and every `method = "random"`
#'   seeding. `method = "provided"` does not read `n` and writes no record;
#'   seeds that already carry one (`sel$seeds`) keep it, marked `FALSE` when
#'   their number is not the record's count. [build_tessellation()] and
#'   [assign_features_to_polygons()] hand the record on, and
#'   [summarize_by_cell()] checks its result against it. A count given as a
#'   number, or read off a profile without an estimate by `cp` or `elbow`,
#'   has no record.
#' @family tessellation
#' @examples
#' library(sf)
#' bnd <- st_sf(geometry = st_sfc(st_polygon(list(rbind(
#'   c(0, 0), c(100, 0), c(100, 100), c(0, 100), c(0, 0)
#' ))), crs = 32632))
#' get_voronoi_seeds(bnd, method = "random", n = 5, set_seed = 1)
#'
#' # Seeds from a cloud with one dense place: density-weighted k-means by
#' # default, from 3 starts; partition = "kmeans" gives the k-means run of
#' # earlier versions.  The points each seed's cell holds, under each kind,
#' # with more seeds than dense places:
#' set.seed(2)
#' cloud <- st_as_sf(data.frame(x = c(rnorm(150, 30, 5), runif(50, 0, 100)),
#'                              y = c(rnorm(150, 30, 5), runif(50, 0, 100))),
#'                   coords = c("x", "y"), crs = 32632)
#' dens <- get_voronoi_seeds(bnd, method = "kmeans", n = 20,
#'                           sample_points = cloud, set_seed = 1)
#' km   <- get_voronoi_seeds(bnd, method = "kmeans", n = 20,
#'                           sample_points = cloud, set_seed = 1,
#'                           partition = "kmeans")
#' rbind(density = sort(attr(dens, "kmeans")$size),
#'       kmeans  = sort(attr(km, "kmeans")$size))   # points per seed
#' @export
get_voronoi_seeds <- function(boundary = NULL,
                              method = c("kmeans", "random", "provided"),
                              n = NULL,
                              seeds = NULL,
                              sample_points = NULL,
                              kmeans_nstart = NULL,
                              kmeans_iter = 100,
                              set_seed = NULL,
                              partition = c("density", "kmeans"),
                              init = NULL,
                              density_power = 1,
                              density_k = 8) {
  method <- match.arg(method)
  partition_given <- !missing(partition)
  # The density settings the caller supplied, as against those left at
  # their defaults: read where the kind in force is this function's own
  # and not one the call named.
  density_given <- c(if (!missing(density_power)) "density_power",
                     if (!missing(density_k)) "density_k")
  partition <- match.arg(partition)

  # --- validation ---
  # The kind of partition the k-means route builds and how its starts are
  # drawn (R/partition.R).  A density setting given where no density
  # partition is built is reported below, once it is known what the call
  # returns: scored centres passed as `n` may be density ones.
  part <- .seed_partition_args(partition, init, density_power, density_k,
                               "get_voronoi_seeds")
  # Said of a kind the call named, not of the function's own: a random or
  # provided seeding that names no kind has asked for nothing it does not get.
  if (!identical(method, "kmeans") &&
      ((partition_given && identical(partition, "density")) || !is.null(init)))
    .warn_and_log(paste0("get_voronoi_seeds(): `partition` and `init` are read by ",
                         "method = \"kmeans\" only; method = \"%s\" ignores them."), method)
  # The other paths, which are not the run this function has always made,
  # check the two counts before the engine reads them: a start count of 0
  # ran one start and recorded none, and NA stopped inside rep().  A NULL
  # number of starts is the kind's own, and is replaced where the argument is
  # read, here and at the fit below, and nowhere before.  method = "random"
  # and "provided" never read it, and neither do scored centres passed beside
  # this function's own run (k-means from random starts), so a call on those
  # paths does not evaluate the argument.
  if (identical(method, "kmeans") && !part$present) {
    kmeans_nstart <- .check_whole_count(.seed_nstart(kmeans_nstart, part$kind),
                                        "kmeans_nstart", "get_voronoi_seeds",
                                        "starts", null_ok = TRUE)
    kmeans_iter   <- .check_whole_count(kmeans_iter, "kmeans_iter", "get_voronoi_seeds",
                                        "iterations")
  }
  if (!is.null(boundary)) .assert_sf(boundary, c("POLYGON", "MULTIPOLYGON"), "boundary")
  # The level-selection step's own answer is accepted as `n`.
  n_arg  <- n
  n_res  <- .resolve_cell_count(n, "n", "get_voronoi_seeds")
  n      <- if (is.null(n_res)) NULL else n_res$n
  n_from <- if (is.null(n_res)) NULL else n_res$from
  if (method %in% c("random", "kmeans") && is.null(n))
    stop("get_voronoi_seeds(): argument 'n' is required for method = '", method,
         "'.", call. = FALSE)
  if (method == "provided") {
    if (is.null(seeds))
      stop("get_voronoi_seeds(): provide 'seeds' (sf POINT) for method = 'provided'.",
           call. = FALSE)
    .assert_sf(seeds, "POINT", "seeds")
    if (!is.null(n) && !identical(as.integer(n), nrow(seeds))) {
      .warn_and_log(
        "get_voronoi_seeds(): method = 'provided' ignores 'n'; returning all %d row(s) of 'seeds' rather than the %d requested.",
        nrow(seeds), as.integer(n)
      )
    }
  }
  if (!is.null(sample_points)) .assert_sf(sample_points, "POINT", "sample_points")

  # A select_resolution() result or a resolution_profile() as `n` carries the
  # centres of the partition the profile scored at that count.  A k-means
  # partition is the Voronoi partition of its centres, so those centres give
  # back the very cells the criteria judged; a new k-means (a different
  # algorithm, different starts, a different cloud) gave other cells, and
  # the count was then defended on a partition nobody built.
  if (identical(method, "kmeans") && !is.null(n_res$seeds)) {
    # The ground these centres stand for: the boundary, else the points they
    # are for, else the hull of the points the profile was computed on, which
    # a profile keeps and a selection's centres carry as their cover.
    scored_cover <- if (!is.null(boundary)) boundary
                    else if (!is.null(sample_points)) sample_points
                    else if (inherits(n_arg, "resolution_profile"))
                      attr(n_arg, "hull", exact = TRUE)
                    else NULL
    .scored_partition_unused(n_res$seeds, part, part$init_given, partition_given,
                             arg = "n", caller = "get_voronoi_seeds",
                             density_given = density_given)
    out <- .use_scored_seeds(n_res$seeds,
                             like  = if (!is.null(boundary)) boundary else sample_points,
                             cloud = sample_points, caller = "get_voronoi_seeds",
                             cover = scored_cover, count_rule = n_res$count_rule)
    cov <- attr(out, "cover", exact = TRUE)
    rec <- attr(out, "partition", exact = TRUE)
    rule <- .count_rule_get(out)
    out <- .align_crs(out, boundary)
    # Kept across the CRS alignment whatever sf does with custom attributes.
    attr(out, "cover") <- cov
    if (!is.null(rec)) attr(out, "partition") <- rec
    attr(out, "n_from") <- n_from
    return(.count_rule_put(out, rule))
  }

  # No scored centres: `density_power` or `density_k` set under partition =
  # "kmeans" is read by nothing in this call.
  .warn_density_unused(part, "get_voronoi_seeds")
  # The same for a random or provided seeding with the density kind in force
  # as this function's own: the warning above speaks of a kind the call did
  # not get, and cannot fire here.
  if (!identical(method, "kmeans") && identical(part$kind, "density") &&
      !partition_given && length(density_given))
    .warn_and_log(paste0("get_voronoi_seeds(): %s %s read by method = \"kmeans\" only; ",
                         "method = \"%s\" ignores %s."),
                  paste(sprintf("`%s`", density_given), collapse = " and "),
                  if (length(density_given) > 1L) "are" else "is", method,
                  if (length(density_given) > 1L) "them" else "it")

  # --- RNG handling (always restored) ---
  cleanup <- .with_seed(set_seed)
  on.exit(cleanup(), add = TRUE)

  # On the sphere for a lon/lat boundary, whatever sf_use_s2() says: with it
  # off, sf printed "although coordinates are longitude/latitude, st_union
  # assumes that they are planar" on every random or k-means seeding.
  boundary_union <- function(b) {
    b <- .safe_make_valid(b)
    if (.is_longlat(b)) .with_s2(sf::st_union(b)) else sf::st_union(b)
  }

  out <- switch(
    method,

    "provided" = {
      s <- seeds
      s$seed_id <- seq_len(nrow(s))
      s$method  <- "provided"
      # Seeds that carry the record of a count are handed on with it (the
      # scored centres of a selection, `sel$seeds`).  Fewer or more of them
      # than the count are not the centres the profile scored, and the
      # record says so.
      carried <- .count_rule_get(s)
      if (!is.null(carried) && !identical(nrow(s), as.integer(carried$cells)))
        s <- .count_rule_put(s, carried, scored_cells = FALSE)
      s
    },

    "random" = {
      if (is.null(boundary))
        stop("get_voronoi_seeds(): random seeds require 'boundary'.", call. = FALSE)
      b <- boundary_union(boundary)
      pts <- .robust_st_sample(b, n)
      pts_sfc <- sf::st_sfc(pts, crs = sf::st_crs(boundary))
      sf::st_sf(seed_id = seq_along(pts_sfc), method = "random",
                geometry = pts_sfc)
    },

    "kmeans" = {
      if (!is.null(sample_points)) {
        cloud <- if (!is.null(boundary)) .align_crs(sample_points, boundary) else sample_points
      } else {
        if (is.null(boundary))
          stop("get_voronoi_seeds(): k-means seeds require 'boundary' (or provide 'sample_points').",
               call. = FALSE)
        b <- boundary_union(boundary)
        cloud_n <- .kmeans_cloud_size(n)
        cloud_geom <- .robust_st_sample(b, cloud_n)
        cloud_sfc <- sf::st_sfc(cloud_geom, crs = sf::st_crs(boundary))
        cloud <- sf::st_sf(geometry = cloud_sfc)
      }

      # If cloud is in lon/lat, project to a local CRS before k-means so
      # that clustering is distance-faithful (k-means in degrees is
      # distorted except in very small areas).
      # A cloud with no CRS gets the lon/lat heuristic every other entry point
      # applies (ensure_projected() warns when it fires).  .is_longlat() alone
      # is FALSE for a missing CRS, so CRS-less degrees were clustered as if
      # they were metres -- 800 points 89 km wide and 111 km tall were split
      # east-west instead of north-south -- while build_tessellation()
      # projected the very same points.
      cloud_crs   <- sf::st_crs(cloud)
      cloud_for_km <- cloud
      cloud_is_ll <- .is_longlat(cloud) ||
        (is.na(cloud_crs) && isTRUE(.looks_like_lonlat(cloud)$lonlat))
      if (cloud_is_ll) {
        cloud_for_km <- ensure_projected(cloud)
        .log_info("get_voronoi_seeds(kmeans): projecting cloud from lon/lat before k-means clustering.")
      }
      xy <- sf::st_coordinates(cloud_for_km)
      # Only the map coordinates.  st_coordinates() on POINT Z returns X, Y, Z,
      # and kmeans() would then cluster in 3-D with elevation dominating the
      # distance -- the returned "seeds" were 3-D points whose XY was wrong
      # (max |dx| 18.4 against a 2-D reference).  voronoi_seeds_kmeans(), the
      # sibling that does the same job, uses the first two columns.
      if (ncol(xy) >= 2L) xy <- xy[, 1:2, drop = FALSE]
      # And only usable rows.  An EMPTY POINT survives as an all-NA row and
      # kmeans() dies with "NA/NaN/Inf in foreign function call (arg 1)", a
      # message naming nothing; voronoi_seeds_kmeans() drops such rows with a
      # warning, so do the same here.
      bad_xy <- !stats::complete.cases(xy) |
                !is.finite(xy[, 1L]) | !is.finite(xy[, 2L])
      if (any(bad_xy)) {
        .warn_and_log("get_voronoi_seeds(kmeans): dropping %d point(s) with empty or non-finite coordinates.",
                      sum(bad_xy))
        xy <- xy[!bad_xy, , drop = FALSE]
      }
      if (nrow(xy) < 2L)
        stop("get_voronoi_seeds(kmeans): fewer than 2 usable points in the sampling cloud.",
             call. = FALSE)

      # Two separate ceilings, both of which stats::kmeans() enforces with a
      # raw message the caller cannot act on.  `n_uniq` is "more cluster
      # centers than distinct data points"; nrow(xy) - 1L is "number of cluster
      # centres must lie between 1 and nrow(x)", which fires at k == nrow(x)
      # exactly -- so n = nrow(sample_points) used to die on a raw kmeans error.
      n_uniq <- nrow(unique(round(xy, 10)))
      k_max  <- min(n_uniq, nrow(xy) - 1L)
      k_use  <- max(1L, min(as.integer(n), k_max))
      if (k_use < n) {
        .warn_and_log("get_voronoi_seeds(kmeans): requested %d seeds but the sampling cloud supports at most %d (%d unique position(s) among %d point(s)); clamping.",
                  as.integer(n), k_use, n_uniq, nrow(xy))
      }

      # One engine for every kind (.partition_centres()).  The default is
      # the run this function has always made, stats::kmeans() from
      # `kmeans_nstart` random starts on the cloud's rows as they come; any
      # other kind or starting rule works on the rows in coordinate order,
      # so its seeds do not depend on the order of the cloud.
      # `kmeans_nstart` is read by the fit when it first needs it, as it
      # always was; NULL is then the kind's own number.
      fit <- .partition_centres(xy, k_use, partition = part$kind, init = part$init,
                                nstart = .seed_nstart(kmeans_nstart, part$kind),
                                iter.max = kmeans_iter,
                                density_power = density_power, density_k = density_k,
                                sort_rows = !part$present, one_call = part$present)
      if (is.null(fit))
        stop(sprintf(paste0("get_voronoi_seeds(kmeans): no start gave a partition of ",
                            "the %d point(s) into %d cells."), nrow(xy), k_use),
             call. = FALSE)
      if (identical(part$kind, "density"))
        .warn_empty_seeds(fit$size, "get_voronoi_seeds")
      # The default run's own warning from stats::kmeans() stands as it is;
      # the other paths say it with a class.
      if (!part$present && isFALSE(fit$converged))
        .warn_not_converged("get_voronoi_seeds", k_use, fit$iter, kmeans_iter,
                            part$kind, "kmeans_iter")
      run_record <- list(iter = as.integer(fit$iter), iter_max = as.integer(kmeans_iter),
                         converged = isTRUE(fit$converged))
      centers <- .seed_centre_matrix(fit, xy)
      centers_sfc <- sf::st_sfc(
        lapply(seq_len(nrow(centers)), function(i) sf::st_point(centers[i, ])),
        crs = sf::st_crs(cloud_for_km)
      )
      s <- sf::st_sf(seed_id = seq_len(k_use), method = "kmeans", geometry = centers_sfc)
      # Back to the cloud's own coordinates if we projected for k-means: its
      # CRS, or for a CRS-less cloud the degrees it was given in, CRS-less
      # again (st_transform() to a missing CRS is an error).
      if (cloud_is_ll)
        sf::st_geometry(s) <- .back_to_input_crs(sf::st_geometry(s), cloud_for_km,
                                                 cloud_crs)
      # Which cloud points fed which seed, and how tight each cluster is:
      # `cluster` is one seed_id per clustered cloud row (`rows` gives those
      # rows' positions in the cloud, since unusable rows were dropped).
      km_record <- list(
        cluster      = as.integer(fit$cluster),
        rows         = which(!bad_xy),
        size         = as.integer(fit$size),
        withinss     = as.numeric(fit$withinss),
        tot_withinss = as.numeric(fit$wss),
        iter         = as.integer(fit$iter),
        nstart       = as.integer(.seed_nstart(kmeans_nstart, part$kind))
      )
      # The density kind minimises a weighted sum, which `tot_withinss` (every
      # point counted once, about its cell's mean) is not.
      if (identical(part$kind, "density")) km_record$objective <- fit$objective
      # Whether the kept start stopped before `kmeans_iter`.  Not added to
      # the default run's record, which is the list it always was.
      if (!part$present) km_record$converged <- run_record$converged
      s
    }
  )

  # Single alignment point for every branch: the seeds are returned in the
  # boundary's CRS whenever one is available.  (The branches above no longer
  # align individually, which made this block unreachable.)
  out <- .align_crs(out, boundary)
  # Where `n` came from, when it was a level-selection result rather than a
  # number; the count itself is nrow(out) (or fewer, see `n`).
  if (!is.null(n_from)) attr(out, "n_from") <- n_from
  # The k-means run behind the seeds, on the kmeans path only.
  if (identical(method, "kmeans")) attr(out, "kmeans") <- km_record
  # The kind that built the regions, when it is not this function's own
  # default (see `partition`); seeds with no record are plain k-means.
  if (identical(method, "kmeans") && !part$present)
    attr(out, "partition") <- c(part[.partition_fields()], run_record)
  # How the count was set, when a selection or a profile gave it: these
  # seeds were placed here, at that count, and are not the centres the
  # profile scored.  `method = "provided"` does not read `n`, and writes none.
  if (!identical(method, "provided"))
    out <- .count_rule_put(out, n_res$count_rule, scored_cells = FALSE)
  # The ground the seeds stand for, set after the alignment so it is kept
  # whatever a transform does with custom attributes: the boundary when one
  # was given (random seeds, a cloud sampled in it, and a cloud the caller
  # placed in it), otherwise the hull of the clustered cloud.
  .set_cover(out, if (!is.null(boundary)) boundary else sample_points)
}


#' Attach the cover of the ground a layer of seeds stands for
#'
#' Seeds lie inside the data, or inside the boundary they were drawn in, so
#' the hull of the seeds is smaller than the ground their Voronoi cells have
#' to hold.  The cover is that ground: the polygons of a boundary, or the
#' convex hull of the points that were clustered (see \code{.cover_of()}),
#' in its own CRS.  \code{create_voronoi_polygons()} reads it when no
#' boundary is passed.
#'
#' @param seeds The sf POINT layer of seeds.
#' @param src An sf/sfc layer of polygons or points, or \code{NULL} when the
#'   call does not say what ground the seeds stand for.
#' @return \code{seeds}, with \code{attr(, "cover")} set when \code{src}
#'   gives one.
#' @keywords internal
#' @noRd
.set_cover <- function(seeds, src) {
  cov <- tryCatch(.cover_of(src), error = function(e) NULL)
  if (!is.null(cov)) attr(seeds, "cover") <- cov
  seeds
}


#' Sample points inside a geometry, retrying when exact sampling fails
#'
#' Tries exact sampling first, then falls back to iterative padding to
#' reach the desired count.
#'
#' @param geom An sfc geometry to sample within.
#' @param n Integer; desired number of sample points.
#' @return An sfc_POINT of length n.
#' @keywords internal
#' @noRd
.robust_st_sample <- function(geom, n) {
  # st_sample() sizes its draw from st_area(), which on lon/lat input needs
  # lwgeom (not a dependency) when sf_use_s2() is FALSE: random and k-means
  # seeding on a lon/lat boundary failed with "package lwgeom required".
  if (.is_longlat(geom) && !isTRUE(sf::sf_use_s2()))
    return(.with_s2(.robust_st_sample(geom, n)))
  # Without lwgeom (not a dependency) sf warns "coordinate ranges not computed
  # along great circles; install package lwgeom to get rid of this warning"
  # on every lon/lat draw, so every random or k-means seeding on a lon/lat
  # boundary raised it, once or twice, in the ordinary case.  Only that
  # warning is muffled; the draw is the same.
  st_sample_quiet <- function(...) withCallingHandlers(
    sf::st_sample(...),
    warning = function(w)
      if (grepl("coordinate ranges not computed along great circles",
                conditionMessage(w), fixed = TRUE))
        invokeRestart("muffleWarning"))
  pts <- try(st_sample_quiet(geom, size = n, type = "random", exact = TRUE),
             silent = TRUE)
  if (inherits(pts, "try-error")) {
    pts <- st_sample_quiet(geom, size = n, type = "random")
  }
  # Pad or trim to exactly n
  max_attempts <- 10L
  attempt <- 0L
  while (length(pts) < n && attempt < max_attempts) {
    attempt <- attempt + 1L
    extra <- st_sample_quiet(geom, size = n - length(pts), type = "random")
    pts <- c(pts, extra)
  }
  if (length(pts) > n) pts <- pts[seq_len(n)]
  if (length(pts) < n)
    .log_warn(".robust_st_sample(): could only generate %d of %d requested points.",
              length(pts), n)
  pts
}


#' K-means seed generation from point coordinates
#'
#' Places `k` seed points at k-means cluster centres of the observed
#' coordinates, so seeds (and the Voronoi cells built from them) follow the
#' sampling density: clusters of observations attract seeds, empty ground gets
#' none. Reach for this when you want cells that follow the data, so that
#' counts per cell vary far less than on a fixed grid over clustered points
#' (the density kind, the default, evens them out once there are more cells
#' than clusters of points; plain k-means, `partition = "kmeans"`, minimises
#' the spread of points around each centre and does not equalise them),
#' which is what keeps per-cell aggregates in
#' [summarize_by_cell()] from resting on one or two observations. Use
#' [voronoi_seeds_random()] instead when you want coverage of the study area
#' rather than of the data, and [get_voronoi_seeds()] to pick between them by
#' name.
#'
#' Lon/lat input is projected first so the k-means distances are metric and not
#' degrees; so is input with no CRS whose coordinates look like lon/lat (the
#' heuristic [ensure_projected()] applies, with its warning), and the seeds
#' come back in the input's own coordinates. Rows with empty or non-finite
#' coordinates are dropped with a warning, and `k` is clamped to the number of
#' distinct positions.
#'
#' A count given as a number gets a fresh partition, which is not the
#' partition [resolution_profile()] or [determine_optimal_levels()] scored
#' that count on: those fit their own run, the profile to every row or to a
#' subsample, and keep the best of their own starts. Pass the
#' [select_resolution()] result or the profile itself as `k` to get the
#' scored one: the centres of the scored partition are then returned (in
#' the CRS of `points_sf`, with `attr(, "n_from")`), and no k-means is run.
#'
#' `partition` says what the seeds minimise. `"kmeans"` is the sum of
#' squared distances from the points to their nearest seed, from `nstart`
#' random starts by [stats::kmeans()] (Hartigan-Wong) unless `init` says
#' otherwise; it was the default of earlier versions. The seeds
#' follow the square root of the point density, so the cells are compact,
#' smaller where the points are dense, and a cell in a dense place still
#' holds several times the points of a cell in a sparse place. `"density"`
#' weights each point by its local density (the number of points in the disc
#' that reaches its `density_k`-th nearest neighbour, over the disc's area)
#' to the power `density_power`, and is the default. At the default power
#' of 1 the seeds follow the point density itself, and the cells are
#' compact: many small cells where the points are dense, a few large ones
#' where they are sparse. They
#' hold about the same number of points each once `k` is above the number
#' of dense places (clusters of points). With fewer cells than that a cell
#' takes in whole clusters, large and small, and the counts are about as
#' uneven as under `"kmeans"`. The density kind is the best of
#' `nstart` starts, each refined by at most 100 steps of Lloyd's iteration
#' (`iter_max` sets another limit, and a kept start that was cut off at the
#' limit is a warning of class `"spatialkit_not_converged"`);
#' it does not depend on the order of the rows or on the unit of the
#' coordinates, and points at one location share a cell. A seed can end
#' with no point at counts from about a fifth of the points up, mostly from
#' `init = "random"`; that is a warning of class `"spatialkit_empty_cells"`.
#' "The partition" in
#' [resolution_profile()] gives the details and a simulated comparison of
#' the two kinds.
#'
#' The centres lie inside the points, so their own hull is smaller than the
#' data's. The result carries the convex hull of `points_sf` as
#' `attr(, "cover")`, and [build_tessellation()] clips the cells to that hull
#' (grown by 2 percent of its diagonal) when no `boundary` is passed, so every
#' point of `points_sf` falls in a cell.
#'
#' @param points_sf An sf object with POINT geometries.
#' @param k Integer; requested number of clusters, treated as an upper bound
#'   only. It is clamped, with a warning, to whichever is smaller of the number
#'   of distinct point positions and `nrow(points_sf) - 1`, because k-means can
#'   produce neither more centres than there are distinct points nor as many
#'   centres as there are rows. Check `nrow()` on the result. A
#'   [select_resolution()] result or a [resolution_profile()] is also
#'   accepted: its scored centres are returned when it carries them,
#'   whatever `partition`, `init` and density settings the call names (see
#'   `n` in [get_voronoi_seeds()] for which arguments a call is taken to
#'   name), and
#'   otherwise its count is used. A profile is read as [get_voronoi_seeds()]
#'   reads one passed as `n`: by `cp` when it was built without an
#'   `estimate`, and by `select_resolution(criterion = "auto")` when it was
#'   built with one (a size for pooled values, the lowest leave-one-out
#'   error for plain means). Such a count raises the warnings and the
#'   error described under `n` in [get_voronoi_seeds()] (classes
#'   `"spatialkit_bound_chosen"`, `"spatialkit_no_structure"` and
#'   `"spatialkit_no_default_criterion"`).
#' @param set_seed Optional integer RNG seed. Default 456, so a call gives the
#'   same seeds every time whatever the session's random-number state; an
#'   outer [set.seed()] does not change them, and the caller's random-number
#'   stream is left as it was. Pass `NULL` to draw the k-means starts from
#'   the session's stream instead (the default of [get_voronoi_seeds()]).
#' @param nstart The number of starts, of which the best is kept. `NULL`
#'   (default) is the kind's own number: 3 for `partition = "density"` and
#'   10 for `partition = "kmeans"`, the default of earlier versions. A number
#'   is run under either kind and either `init`: `nstart = 10` gives the
#'   density kind the starts it had when 10 was the default.
#' @param partition What the seeds minimise. `"density"` (default): the sum
#'   of squared distances from the points to their nearest seed, each point
#'   weighted by its local density to the power `density_power`; compact
#'   cells that, at the default power, hold about the same number of points
#'   each once `k` is above the number of dense places, and counts about as
#'   uneven as `"kmeans"` gives below that. `"kmeans"`: the same sum
#'   unweighted; compact cells whose counts are several times higher in
#'   dense places than in sparse ones. It is the run of earlier versions,
#'   whose default it was: [stats::kmeans()] from 10 random starts on the
#'   rows as they come. Without the 'FNN' package the density kind
#'   finds its weights by comparing every pair of distinct locations, so its
#'   time grows with the square of their number; the seeds are the same, and
#'   nothing stops the call. See "Details".
#' @param init How each start is drawn: `"random"` (`k` distinct locations
#'   at random, as [stats::kmeans()] draws them) or `"kmeans++"` (the seeds
#'   one at a time, each in proportion to its squared distance from those
#'   already drawn, times its weight under `partition = "density"`). `NULL`
#'   (default) is each kind's own rule: `"random"` for `partition =
#'   "kmeans"`, the run this function has always made, and `"kmeans++"` for
#'   `partition = "density"`. Under `partition = "kmeans"` from random
#'   starts the points are clustered in the order their rows come in; with
#'   `init = "kmeans++"`, and under `partition = "density"`, the rows are put
#'   in coordinate order first, so the seeds do not depend on their order.
#' @param density_power For `partition = "density"`: the power the local
#'   density is raised to for a point's weight, from 0 to 2. Default 1, which
#'   evens out the cell counts once the cells outnumber the dense places; 0
#'   weights every point alike, and values
#'   above 1 give the cells of dense places fewer points than those of
#'   sparse places. Setting it under `partition = "kmeans"` is a warning.
#' @param density_k For `partition = "density"`: the number of nearest
#'   neighbours the local density is read from. Default 8.
#' @param iter_max The most iterations a start runs for: a single positive
#'   whole number, or `NULL` (default) for the kind's own limit, 50 for
#'   `partition = "kmeans"` and 100 steps of Lloyd's iteration for
#'   `partition = "density"`. Under `partition = "density"`, and with `init =
#'   "kmeans++"`, a kept start that was cut off at the limit is a warning of
#'   class `"spatialkit_not_converged"`; `partition = "kmeans"` from random
#'   starts keeps the warning [stats::kmeans()] gives.
#' @return An sf object of **at most** `k` cluster-centre POINTs (fewer when
#'   `k` exceeds the number of distinct positions), with `seed_id` and
#'   `method = "kmeans"` columns matching [get_voronoi_seeds()], and
#'   `attr(, "cover")`, the convex hull of `points_sf` (an sfc in its CRS).
#'   Seeds placed with `partition = "density"` or `init = "kmeans++"`, and
#'   the scored centres of a profile of the density kind or from random
#'   starts, also carry `attr(, "partition")`: a list with
#'   `kind`, `init`, `density_power` and `density_k`, and for seeds this
#'   call placed `iter`, `iter_max` and `converged` (the iterations the kept
#'   start took, the limit, and whether it stopped before the limit).
#'   [build_tessellation()] returns the list as `params$partition`. Seeds
#'   with no such record are plain k-means: from random starts when this call
#'   placed them, or the scored centres of a k-means profile from k-means++
#'   starts. When `k` is a
#'   [select_resolution()] result that carries `count_rule`, or a profile
#'   built with an `estimate`, the seeds carry that record as
#'   `attr(, "count_rule")`, with `scored_cells` `TRUE` for the scored
#'   centres and `FALSE` for centres this call fits itself at that count
#'   (see [get_voronoi_seeds()]).
#' @family tessellation
#' @examples
#' library(sf)
#' set.seed(1)
#' pts <- st_as_sf(
#'   data.frame(x = 5e5 + runif(100, 0, 1000), y = 5e6 + runif(100, 0, 1000)),
#'   coords = c("x", "y"), crs = 32632
#' )
#' seeds <- voronoi_seeds_kmeans(pts, k = 8)   # density-weighted, 3 starts
#' nrow(seeds)   # at most 8: one seed per non-empty cluster
#' seeds         # the cluster centres, as an sf POINT layer in the points' CRS
#' # The k-means seeds of earlier versions: 10 random starts.
#' voronoi_seeds_kmeans(pts, k = 8, partition = "kmeans")
#' @export
voronoi_seeds_kmeans <- function(points_sf, k, set_seed = 456, nstart = NULL,
                                 partition = c("density", "kmeans"), init = NULL,
                                 density_power = 1, density_k = 8, iter_max = NULL) {
  .assert_sf(points_sf, "POINT", "points_sf")
  partition_given <- !missing(partition)
  # The density settings the caller supplied (see get_voronoi_seeds()).
  density_given <- c(if (!missing(density_power)) "density_power",
                     if (!missing(density_k)) "density_k")
  partition <- match.arg(partition)
  # How the count was set, when a selection or a profile gives it.
  count_rule <- NULL
  # The kind of partition and how its starts are drawn (R/partition.R).  A
  # density setting given where no density partition is built is reported
  # below, once it is known what the call returns: scored centres passed as
  # `k` may be density ones.
  part <- .seed_partition_args(partition, init, density_power, density_k,
                               "voronoi_seeds_kmeans")
  if (!is.null(iter_max))
    iter_max <- .check_whole_count(iter_max, "iter_max", "voronoi_seeds_kmeans",
                                   "iterations", null_ok = TRUE)
  # A select_resolution() result or a resolution_profile() as `k`: the
  # centres of the partition it scored at the chosen count when it carries
  # them (as get_voronoi_seeds() does), otherwise the count alone.
  if (inherits(k, c("resolution_selection", "resolution_profile"))) {
    k_res <- .resolve_cell_count(k, "k", "voronoi_seeds_kmeans")
    if (!is.null(k_res$seeds)) {
      .scored_partition_unused(k_res$seeds, part, part$init_given, partition_given,
                               arg = "k", caller = "voronoi_seeds_kmeans",
                               density_given = density_given)
      out <- .use_scored_seeds(k_res$seeds, like = points_sf, cloud = NULL,
                               caller = "voronoi_seeds_kmeans", cover = points_sf,
                               count_rule = k_res$count_rule)
      attr(out, "n_from") <- k_res$from
      return(out)
    }
    k <- k_res$n
    count_rule <- k_res$count_rule
  }
  # No scored centres: `density_power` or `density_k` set under partition =
  # "kmeans" is read by nothing in this call.
  .warn_density_unused(part, "voronoi_seeds_kmeans")
  # `k` was never validated here, unlike get_voronoi_seeds(), which routes
  # `n` through .resolve_cell_count(): k = NA or a length-2 vector aborted on
  # the clamp below, k above .Machine$integer.max became NA through
  # as.integer(), and k = 0 or negative silently returned ONE seed -- which
  # the @return section's "at most k" does not describe.
  .check_scalar(k, "k", "voronoi_seeds_kmeans", min = 1,
                max = .Machine$integer.max,
                what = "a single positive number of seeds")

  # NULL is the kind's own number of starts.  It is replaced here, where the
  # argument is first read: a call that returned scored centres above never
  # read it, and does not evaluate it.
  nstart <- .seed_nstart(nstart, part$kind)
  .check_scalar(nstart, "nstart", "voronoi_seeds_kmeans", min = 1,
                max = .Machine$integer.max,
                what = "NULL or a single positive number of k-means restarts")

  # Project to metric CRS if lon/lat to make k-means distance-faithful.  A
  # CRS-less layer gets the lon/lat heuristic ensure_projected() applies
  # everywhere else (with its warning); .is_longlat() is FALSE for a missing
  # CRS, so CRS-less degrees used to be clustered as if they were metres.
  crs_in     <- sf::st_crs(points_sf)
  pts_for_km <- points_sf
  pts_is_ll  <- .is_longlat(points_sf) ||
    (is.na(crs_in) && isTRUE(.looks_like_lonlat(points_sf)$lonlat))
  if (pts_is_ll) {
    pts_for_km <- ensure_projected(points_sf)
  }
  coords <- sf::st_coordinates(pts_for_km)
  if (ncol(coords) >= 2L) coords <- coords[, 1:2, drop = FALSE]
  # st_coordinates() yields one all-NA row per EMPTY POINT rather than zero
  # rows, so a row-count check alone would let empty geometries through to
  # kmeans(), which fails with "NA/NaN/Inf in foreign function call".
  finite_rows <- stats::complete.cases(coords) &
    apply(coords, 1L, function(r) all(is.finite(r)))
  n_drop <- sum(!finite_rows)
  if (n_drop > 0L) {
    .warn_and_log("voronoi_seeds_kmeans(): dropping %d point(s) with empty or non-finite coordinates.",
              n_drop)
    coords <- coords[finite_rows, , drop = FALSE]
  }
  n <- nrow(coords)
  if (n == 0L)
    stop("voronoi_seeds_kmeans(): `points_sf` has no usable coordinates; ",
         "nothing to cluster.", call. = FALSE)

  # As in get_voronoi_seeds(): stats::kmeans() refuses k > distinct rows AND
  # k >= nrow(x), the latter with "number of cluster centres must lie between 1
  # and nrow(x)".  k = nrow(points_sf) therefore used to error rather than
  # clamp, despite `@return` promising "at most k".
  n_uniq <- nrow(unique(round(coords, 10)))
  k_max  <- min(n_uniq, n - 1L)
  k_use  <- max(1L, min(as.integer(k), k_max))
  if (k_use < k) {
    .warn_and_log("voronoi_seeds_kmeans(): requested %d seeds but only %d unique positions among %d point(s); clamping to %d.",
              as.integer(k), n_uniq, n, k_use)
  }

  cleanup <- .with_seed(set_seed)
  on.exit(cleanup(), add = TRUE)
  # One engine for every kind (.partition_centres()).  The default is the run
  # this function has always made, stats::kmeans() from `nstart` random
  # starts on the rows as they come; any other kind or starting rule works on
  # the rows in coordinate order, so its seeds do not depend on the order of
  # `points_sf`.  Iterations per start, unless `iter_max` sets them: 50 for
  # k-means, as always, and 100 Lloyd steps for the density kind (see
  # resolution_profile()).
  limit <- if (!is.null(iter_max)) iter_max
           else if (identical(part$kind, "density")) 100L else 50
  fit <- .partition_centres(coords, k_use, partition = part$kind, init = part$init,
                            nstart = as.integer(nstart), iter.max = limit,
                            density_power = density_power, density_k = density_k,
                            sort_rows = !part$present, one_call = part$present)
  if (is.null(fit))
    stop(sprintf(paste0("voronoi_seeds_kmeans(): no start gave a partition of the %d ",
                        "point(s) into %d cells."), n, k_use), call. = FALSE)
  if (identical(part$kind, "density"))
    .warn_empty_seeds(fit$size, "voronoi_seeds_kmeans")
  # The default run's own warning from stats::kmeans() stands as it is; the
  # other paths say it with a class.
  if (!part$present && isFALSE(fit$converged))
    .warn_not_converged("voronoi_seeds_kmeans", k_use, fit$iter, limit, part$kind,
                        "iter_max")
  cent <- as.data.frame(.seed_centre_matrix(fit, coords)); names(cent) <- c("x", "y")
  result <- sf::st_as_sf(cent, coords = c("x", "y"), crs = sf::st_crs(pts_for_km))
  # Back to the input's own coordinates if we projected: its CRS, or for a
  # CRS-less layer the degrees it was given in, CRS-less again
  # (st_transform() to a missing CRS is an error).
  if (pts_is_ll) {
    sf::st_geometry(result) <- .back_to_input_crs(sf::st_geometry(result),
                                                  pts_for_km, crs_in)
  }
  # Same output contract as get_voronoi_seeds(), so the three seeding
  # functions are drop-in interchangeable.
  result$seed_id <- seq_len(nrow(result))
  result$method  <- "kmeans"
  # The kind that built the regions, when it is not this function's own
  # default (see `partition`); seeds with no record are plain k-means.
  if (!part$present)
    attr(result, "partition") <- c(part[.partition_fields()],
                                   list(iter = as.integer(fit$iter),
                                        iter_max = as.integer(limit),
                                        converged = isTRUE(fit$converged)))
  # How the count was set, when a selection or a profile gave it and its
  # scored centres were not at hand: these centres were fitted here.
  result <- .count_rule_put(result, count_rule, scored_cells = FALSE)
  # The ground the centres stand for: the hull of the points clustered.
  .set_cover(result, points_sf)
}


#' Random seed generation within a polygonal boundary
#'
#' Draws `k` seed points uniformly at random inside `boundary`, ignoring where
#' the observations are. Reach for this when the cells should cover the study
#' area evenly (so that sparsely sampled ground still gets its own cells and
#' is visibly under-sampled in the results) instead of concentrating
#' resolution where the data already are, which is what
#' [voronoi_seeds_kmeans()] does. It is also the honest choice for a null or
#' sensitivity comparison: re-running an analysis over several random seedings
#' shows how much of a result depends on one particular tessellation.
#'
#' Seeds are drawn uniformly inside the polygon. When a draw falls short of `k`
#' it is topped up with further uniform draws from the same polygon, so the
#' result has exactly `k` seeds; only a geometry that still yields too few
#' after ten top-ups returns fewer, and that shortfall is logged.
#'
#' @param boundary An sf or sfc polygonal object.
#' @param k Integer; number of random seeds.
#' @param set_seed Optional integer RNG seed. Default `NULL`: the seeds are
#'   drawn from the session's random-number stream, so consecutive calls give
#'   different seedings and [set.seed()] before a call makes it reproducible.
#'   Pass a number to get the same seeds whatever that stream holds; the
#'   caller's stream is then left as it was. The default used to be 456,
#'   which made every call return the same "random" seeding, even inside a
#'   loop over [set.seed()].
#' @return An sf object of `k` random POINTs (fewer only in the degenerate
#'   case above), with `seed_id` and `method = "random"` columns matching
#'   [get_voronoi_seeds()], and `attr(, "cover")`, the polygons of `boundary`:
#'   [build_tessellation()] clips the cells to it (grown by 2 percent of its
#'   diagonal) when no `boundary` is passed there, so the cells fill the area
#'   the seeds were drawn in and not only the seeds' own hull.
#' @family tessellation
#' @examples
#' library(sf)
#' bnd <- st_sf(geometry = st_sfc(st_polygon(list(rbind(
#'   c(0, 0), c(100, 0), c(100, 100), c(0, 100), c(0, 0)
#' ))), crs = 32632))
#' set.seed(1)
#' seeds <- voronoi_seeds_random(bnd, k = 10)
#' nrow(seeds)   # 10
#' seeds
#' # Another call is another seeding; set_seed pins one.
#' identical(st_coordinates(voronoi_seeds_random(bnd, k = 10, set_seed = 7)),
#'           st_coordinates(voronoi_seeds_random(bnd, k = 10, set_seed = 7)))
#' @export
voronoi_seeds_random <- function(boundary, k, set_seed = NULL) {
  # `k` was never validated here, unlike get_voronoi_seeds(), which routes
  # `n` through .resolve_cell_count(): k = NA or a length-2 vector aborted on
  # the clamp below, k above .Machine$integer.max became NA through
  # as.integer(), and k = 0 or negative silently returned ONE seed -- which
  # the @return section's "at most k" does not describe.
  .check_scalar(k, "k", "voronoi_seeds_random", min = 1,
                max = .Machine$integer.max,
                what = "a single positive number of seeds")
  # `@param boundary` documents sf *or* sfc, and .assert_sf() only accepts sf.
  if (inherits(boundary, "sfc")) boundary <- sf::st_as_sf(boundary)
  .assert_sf(boundary, c("POLYGON", "MULTIPOLYGON"), "boundary")

  cleanup <- .with_seed(set_seed)
  on.exit(cleanup(), add = TRUE)
  # On the sphere for lon/lat, as in get_voronoi_seeds(): with s2 off sf
  # otherwise printed its planar-union message on every call.
  geom <- if (.is_longlat(boundary)) .with_s2(sf::st_union(boundary)) else
    sf::st_union(boundary)
  pts <- .robust_st_sample(geom, k)
  out <- sf::st_sf(geometry = pts) |> sf::st_set_crs(sf::st_crs(boundary))
  # Same output contract as get_voronoi_seeds(), so the three seeding
  # functions are drop-in interchangeable.
  out$seed_id <- seq_len(nrow(out))
  out$method  <- "random"
  # The ground the seeds stand for: the area they were drawn in.
  .set_cover(out, boundary)
}


#' Size of the point cloud k-means seeding draws inside the boundary
#'
#' Fifty candidate points per requested seed, and never fewer than 2,000.
#' Computed in double precision: \code{50L * as.integer(n)} overflowed to
#' \code{NA} (with a warning) for \code{n} above 42,949,672, and
#' \code{max(2000L, NA)} is \code{NA}, which \code{st_sample()} then refused
#' with an unrelated message.
#'
#' @param n Number of seeds requested.
#' @return A single positive number.
#' @keywords internal
#' @noRd
.kmeans_cloud_size <- function(n) {
  max(2000, 50 * as.numeric(n))
}


#' Seeds from the partition a resolution profile scored
#'
#' \code{resolution_profile()} keeps, for every level, the centres of the
#' best-of-\code{nstart} k-means++ partition its criteria were computed on.
#' A k-means partition is the Voronoi partition of its centres, so these are
#' the seeds of exactly the cells that were scored.
#'
#' @param seeds The centres, from \code{.resolve_cell_count()}'s \code{seeds}
#'   (an sf POINT layer in the profile's projected CRS, with
#'   \code{attr(, "nstart")}).
#' @param like The layer whose CRS the seeds are returned in (the boundary or
#'   the points being seeded), or \code{NULL} for the profile's own CRS.
#' @param cloud Optional sf POINT layer to assign to the seeds for
#'   \code{attr(, "kmeans")}; \code{NULL} for none.
#' @param caller Name for the log line.
#' @param cover The ground the centres stand for, as \code{.set_cover()}
#'   takes it: the boundary, the points being seeded, or the hull the profile
#'   kept of the points it was computed on.  \code{NULL} falls back to a
#'   cover the centres themselves carry (\code{attr(seeds, "cover")}), and
#'   to none when there is none.
#' @param count_rule The record of how the count was set
#'   (\code{.resolve_cell_count()}'s \code{count_rule}), or \code{NULL} for a
#'   count that has none.
#' @return An sf POINT layer with \code{seed_id} and \code{method = "kmeans"}
#'   columns, \code{attr(, "kmeans")} when \code{cloud} is given,
#'   \code{attr(, "cover")} when a cover is known, and
#'   \code{attr(, "count_rule")}, marked \code{scored_cells = TRUE}, when
#'   \code{count_rule} is a record: these are the centres the profile scored.
#' @keywords internal
#' @noRd
.use_scored_seeds <- function(seeds, like, cloud, caller, cover = NULL, count_rule = NULL) {
  nstart <- attr(seeds, "nstart", exact = TRUE)
  # The kind that built the scored partition, when the profile recorded one.
  rec <- attr(seeds, "partition", exact = TRUE)
  if (is.null(cover)) cover <- attr(seeds, "cover", exact = TRUE)
  moved  <- .seeds_in_crs_of(seeds, like)
  out <- sf::st_sf(seed_id = seq_len(nrow(moved)), method = "kmeans",
                   geometry = sf::st_geometry(moved))
  if (is.null(rec))
    .log_info(paste0("%s(): returning the %d centres of the partition resolution_profile() ",
                     "scored at that count (the best of %s k-means++ restarts) rather than ",
                     "running a new k-means, so their Voronoi cells are the cells the ",
                     "criteria judged. Pass the count as a number for a fresh k-means."),
              caller, nrow(out), if (is.null(nstart)) "its" else format(nstart))
  else
    .log_info(paste0("%s(): returning the %d centres of the partition resolution_profile() ",
                     "scored at that count (partition = \"%s\", the best of %s) rather ",
                     "than running a new one, so their Voronoi cells are the cells the ",
                     "criteria judged. Pass the count as a number for a fresh run."),
              caller, nrow(out), rec$kind, .partition_label(rec, nstart))
  if (!is.null(cloud)) attr(out, "kmeans") <- .scored_seeds_record(cloud, seeds, nstart)
  if (!is.null(rec)) attr(out, "partition") <- rec
  # These are the centres the profile scored at the count, so their Voronoi
  # cells are the cells its numbers describe, and the record says so.
  out <- .count_rule_put(out, count_rule, scored_cells = TRUE)
  # The centres lie inside the points they were fitted to; their cells have
  # to reach the points' own hull, not the centres'.
  .set_cover(out, cover)
}


#' The number of starts a seeding function runs
#'
#' \code{NULL} is the kind's own number: the 10 the seeding functions have
#' always run for \code{partition = "kmeans"}, whichever starting rule, and
#' \code{.NSTART_DENSITY} for \code{partition = "density"}.  Anything else
#' is returned as it is, for the caller's check.  It is called where the
#' argument is first read, so that a call on a path that never reads the
#' number of starts does not evaluate it.
#'
#' @param nstart The argument: \code{kmeans_nstart} of
#'   \code{get_voronoi_seeds()}, \code{nstart} of
#'   \code{voronoi_seeds_kmeans()}.
#' @param kind The kind in force, \code{"kmeans"} or \code{"density"}.
#' @return \code{nstart}, or the kind's own number when it is \code{NULL}.
#' @keywords internal
#' @noRd
.seed_nstart <- function(nstart, kind) {
  if (!is.null(nstart)) return(nstart)
  if (identical(kind, "density")) .NSTART_DENSITY else 10
}


#' The partition arguments of a seeding function, resolved
#'
#' The seeding functions have always run \code{stats::kmeans()} from random
#' starts, so \code{init = NULL}, their default, means each kind's own rule:
#' random starts for \code{partition = "kmeans"}, which leaves every existing
#' call as it was, and k-means++ for \code{partition = "density"}.
#'
#' @param partition The kind, already matched.
#' @param init \code{NULL}, \code{"kmeans++"} or \code{"random"}.
#' @param density_power,density_k,caller See \code{.check_partition_args()}.
#' @return The list of \code{.check_partition_args()}, with \code{present}
#'   (whether this is the run the function has always made) and
#'   \code{init_given} (\code{init} as matched, \code{NULL} when it was not
#'   set).  A density setting that \code{partition = "kmeans"} does not read
#'   is not reported here: the caller does that with
#'   \code{.warn_density_unused()} once it knows the call runs a k-means,
#'   since scored centres passed as the count may be density ones.
#' @keywords internal
#' @noRd
.seed_partition_args <- function(partition, init, density_power, density_k, caller) {
  init <- .match_init(init, caller, null_ok = TRUE)
  rule <- if (!is.null(init)) init
          else if (identical(partition, "density")) "kmeans++" else "random"
  part <- .check_partition_args(partition, rule, density_power, density_k, caller,
                                warn_unused = FALSE)
  part$init_given <- init
  part$present <- identical(part$kind, "kmeans") && identical(rule, "random")
  part
}


#' Say that scored centres are returned in place of the kind asked for
#'
#' A selection or a profile passed as the count carries the centres of the
#' partition it scored, and those are what the seeding functions return.  A
#' \code{partition}, an \code{init} or a density setting the same call names
#' is then not used; when it differs from what built the scored partition,
#' the caller would otherwise get regions of one kind having asked for
#' another.
#'
#' The kind in force in a call is the one its \code{partition} resolves to:
#' the kind the caller supplied, or else the function's own.  The density
#' settings a call names follow from it.  With the density kind supplied,
#' both, at the values in force: naming the kind names its settings.  With
#' the density kind in force as the function's own, those the caller supplied
#' (\code{density_given}).  With the k-means kind, those not at their
#' defaults.  So a call that names no partition argument is told nothing,
#' whichever kind and settings the profile was built with.
#'
#' @param seeds The scored centres.
#' @param part The call's partition arguments, as
#'   \code{.seed_partition_args()} resolved them.
#' @param init The call's \code{init}, which may be \code{NULL}.
#' @param partition_given Whether \code{partition} was set in the call.
#' @param arg,caller Names for the message.
#' @param density_given Which of \code{density_power} and \code{density_k}
#'   the caller supplied, by name.
#' @return \code{NULL}, invisibly; a warning of class
#'   \code{spatialkit_partition_unused} when the two differ.
#' @keywords internal
#' @noRd
.scored_partition_unused <- function(seeds, part, init, partition_given, arg, caller,
                                     density_given = character(0)) {
  rec <- attr(seeds, "partition", exact = TRUE)
  # A profile with no record is k-means from k-means++ starts.
  kind <- if (is.null(rec)) "kmeans" else rec$kind
  rule <- if (is.null(rec)) "kmeans++" else rec$init
  # The density kind in force as the function's own, and not as a kind the
  # call named.
  own_density <- identical(part$kind, "density") && !isTRUE(partition_given)
  # The density settings the call names.
  named <- if (own_density) density_given
           else if (identical(part$kind, "density")) c("density_power", "density_k")
           else part$unused
  differs <- (isTRUE(partition_given) && !identical(part$kind, kind)) ||
    (!is.null(init) && !identical(init, rule)) ||
    # Scored density centres, with another power or neighbour count named.
    (identical(kind, "density") &&
       (("density_power" %in% named &&
           !isTRUE(all.equal(part$asked$density_power, rec$density_power))) ||
          ("density_k" %in% named &&
             !identical(as.integer(part$asked$density_k), as.integer(rec$density_k))))) ||
    # Scored k-means centres, with a density setting supplied to a function
    # whose own kind is the density kind: the setting is read by nothing.
    (!identical(kind, "density") && own_density && length(density_given) > 0L)
  if (differs)
    .warn_classed("spatialkit_partition_unused",
                  paste0("%s(): `%s` carries the centres of the partition the profile ",
                         "scored, built with partition = \"%s\" from %s starts%s, and ",
                         "those are returned: the `partition`, `init` and density ",
                         "settings of this call are not used. Build the profile with ",
                         "them, or pass the count as a number for a fresh run."),
                  caller, arg, kind, if (identical(rule, "random")) "random" else "k-means++",
                  if (identical(kind, "density"))
                    sprintf(" (power %s, %d neighbours)", format(rec$density_power),
                            as.integer(rec$density_k))
                  else "",
                  data = list(kind = kind, init = rule))
  # Scored k-means centres with a density setting named and no kind: the
  # setting is read by nothing, as in a fresh k-means run.
  else if (!identical(kind, "density"))
    .warn_density_unused(part, caller)
  invisible(NULL)
}


#' The centres of an engine fit as the matrix the seeds are built from
#'
#' A \code{stats::kmeans()} fit hands over its own \code{centers}, row and
#' column names included, as the seeding functions always read it.  The
#' density kind's centres get the same names, so the layers of seeds have one
#' structure whatever built them.
#'
#' @param fit A \code{.partition_centres()} result.
#' @param xy The coordinate matrix that was clustered.
#' @return A matrix of centres with \code{1..k} as row names.
#' @keywords internal
#' @noRd
.seed_centre_matrix <- function(fit, xy) {
  if (!is.null(fit$km)) return(fit$km$centers)
  cen <- fit$centers
  dimnames(cen) <- list(as.character(seq_len(nrow(cen))), colnames(xy)[1:2])
  cen
}


#' Scored seeds in the CRS of the layer they are for
#'
#' As \code{.align_crs()}, and for a CRS-less layer the coordinates it is
#' in: degrees when they look like lon/lat (\code{.looks_like_lonlat()}, the
#' heuristic the profile's own \code{ensure_projected()} applied), otherwise
#' the profile's own numbers, CRS-less.
#'
#' @keywords internal
#' @noRd
.seeds_in_crs_of <- function(seeds, like) {
  if (is.null(like)) return(seeds)
  to <- sf::st_crs(like); from <- sf::st_crs(seeds)
  if (is.na(from)) return(seeds)
  if (is.na(to)) {
    if (isTRUE(.looks_like_lonlat(like)$lonlat))
      seeds <- sf::st_transform(seeds, 4326)
    return(sf::st_set_crs(seeds, NA))
  }
  if (from == to) seeds else sf::st_transform(seeds, to)
}


#' Which scored seed each cloud point falls to
#'
#' The \code{attr(, "kmeans")} record of \code{get_voronoi_seeds()} for seeds
#' taken from a resolution profile: each usable cloud point is assigned to its
#' nearest seed in the profile's projected CRS (the seed whose Voronoi cell
#' holds it), with the sizes and within-cell sums of squares that assignment
#' gives.  \code{iter} is \code{NA}: no k-means ran here.
#'
#' @keywords internal
#' @noRd
.scored_seeds_record <- function(cloud, seeds, nstart) {
  scrs <- sf::st_crs(seeds); ccrs <- sf::st_crs(cloud)
  if (!is.na(scrs)) {
    if (is.na(ccrs))
      cloud <- if (isTRUE(.looks_like_lonlat(cloud)$lonlat))
        sf::st_transform(sf::st_set_crs(cloud, 4326), scrs)
      else sf::st_set_crs(cloud, scrs)
    else if (ccrs != scrs)
      cloud <- sf::st_transform(cloud, scrs)
  } else if (!is.na(ccrs)) {
    cloud <- sf::st_set_crs(cloud, NA)
  }
  xy <- sf::st_coordinates(cloud)
  xy <- if (ncol(xy) >= 2L) xy[, 1:2, drop = FALSE] else matrix(NA_real_, nrow(cloud), 2L)
  ok <- stats::complete.cases(xy) & is.finite(xy[, 1L]) & is.finite(xy[, 2L])
  sxy <- sf::st_coordinates(seeds)[, 1:2, drop = FALSE]
  k   <- nrow(sxy)
  cl  <- integer(0); d2 <- numeric(0)
  if (any(ok)) {
    cl <- as.integer(sf::st_nearest_feature(cloud[ok, ], seeds))
    d2 <- rowSums((xy[ok, , drop = FALSE] - sxy[cl, , drop = FALSE])^2)
  }
  # Seeds of the density kind are weighted means, which do not sit at the
  # mean of their cell's points; the within-cell sum of squares is then taken
  # about that mean, as the profile and a fresh run report it.  A k-means
  # seed is the mean of the cell it was fitted to.
  if (identical(attr(seeds, "partition", exact = TRUE)$kind, "density") && length(cl)) {
    xo <- xy[ok, , drop = FALSE]
    mu <- rowsum(xo, cl) / as.numeric(table(cl))
    d2 <- rowSums((xo - mu[as.character(cl), , drop = FALSE])^2)
  }
  withinss <- as.numeric(tapply(d2, factor(cl, levels = seq_len(k)), sum, default = 0))
  list(
    cluster      = cl,
    rows         = which(ok),
    size         = tabulate(cl, nbins = k),
    withinss     = withinss,
    tot_withinss = sum(withinss),
    iter         = NA_integer_,
    nstart       = if (is.null(nstart)) NA_integer_ else as.integer(nstart)
  )
}
