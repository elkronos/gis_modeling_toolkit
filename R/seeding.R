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
#' given, otherwise the convex hull of `sample_points` (for lon/lat points,
#' the hull of their coordinates in degrees, with the longitudes of points on
#' both sides of 180 degrees taken from 0 to 360, -170 as 190), otherwise
#' (for the scored centres of a [resolution_profile()] or a
#' [select_resolution()] result passed as `n`) the hull of the points the
#' profile was computed on.
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
#' within that boundary, its cell is built if any of it lies inside: the part
#' of the study area nearer the seed than any other seed.
#' [build_tessellation()] gives the seed that cell in `index`, as it gives
#' every seed the cell it generates, so `seeds$cell_id <- tess$index` holds
#' each seed's cell for [cell_values_at()]. A seed with every part of the
#' study area nearer another seed has no cell, `NA` in `index`, and a
#' warning of class `spatialkit_unassigned` gives the number.
#'
#' @param boundary Optional polygonal sf object defining the sampling area.
#'   One that mixes XY and XYZ (or XYM) geometries is read without Z and M,
#'   since sf reads a layer in the dimension of its first geometry. So is one whose
#'   geometries carry M (XYM or XYZM), which GEOS refuses. Where seeds or
#'   the k-means cloud are sampled inside it, its rows are repaired with
#'   [sf::st_make_valid()] before they are unioned. In
#'   longitude and latitude the union runs with s2, and a row s2 still
#'   refuses (crossing edges, a vertex listed twice) is repaired again, with
#'   crossing edges split as GEOS splits them in a projected CRS; a row s2
#'   refuses even then stops the call with an error of class
#'   `"spatialkit_cells_degenerate"` that names it and says to transform
#'   `boundary` to a projected CRS.
#'   A row that is empty, or has a coordinate that is not finite or lies
#'   beyond 1e150 (read as given, and again after any move into another
#'   CRS), covers no ground and is left out; a boundary with no such row
#'   left is refused with an error that says so. A boundary in a geographic
#'   CRS with a vertex at a latitude beyond 90 degrees, which is no place on
#'   the Earth, is refused with an error of class
#'   `"spatialkit_cells_degenerate"` that names its rows and says to check
#'   its CRS.
#'   A row the repair leaves without a polygon (a ring with no height:
#'   every vertex on one line, or in degrees on one parallel) is left out
#'   as well, and a boundary with no polygon left is refused. In
#'   longitude and latitude, a boundary whose bounding box is over 10,000
#'   times its area (small pieces far apart), which sf does not sample whole,
#'   is sampled piece by piece, the points shared among the pieces in
#'   proportion to their areas; a piece whose own box is over 10,000 times
#'   its area stops the call with an error that says to transform
#'   `boundary` to a projected CRS. A boundary whose bounding box has no
#'   height or no width (three vertices on one parallel) is refused.
#'   sf draws the points in the bounding box and keeps those inside, so it
#'   draws the number asked for times the ratio of the box to the area. In
#'   any CRS, a boundary that would take more than 10^6 points drawn in its
#'   box, or 10 times the points asked for if more (small pieces far apart,
#'   a thin strip at a slant), is sampled piece by piece in the same way
#'   (or whole, in draws of at most 10^6 points in its box, where that takes
#'   less time), and a piece that would still take more is drawn in such
#'   draws. Where the points drawn in all would exceed 2e7 (or 10 times the
#'   points asked for, if more), the call stops before any draw with an
#'   error that gives the ratio and says to pass `sample_points`
#'   (`method = "kmeans"`, whose cloud is `max(2000, 50 * n)` points) or to
#'   ask for fewer seeds.
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
#'   (its first element is used; its `NA` is refused with a note on where it
#'   came from), a [select_resolution()] result (its
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
#'   cells the criteria judged. On a profile without an estimate whose
#'   centres were fitted to a subsample, [select_resolution()] refits them
#'   to every point first, and those centres are returned, with
#'   `attr(, "refit")`: their cells are near the cells judged but not the
#'   same. `sample_points`, when given, is then only
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
#'   Methods `"kmeans"` and `"random"` do not read it, and ignore it with a
#'   warning; so does a call that passes `seed =`, which R matches to
#'   `seeds` (the seed of the random stream is `set_seed`).
#'   Seeds in another CRS than `boundary` are returned in the boundary's CRS;
#'   seeds that mix XY and XYZ (or XYM) points are then returned in XY, with a
#'   warning of class `"spatialkit_zm_dropped"`, since sf reads a layer in
#'   the dimension of its first geometry.
#' @param sample_points Optional sf POINT cloud for k-means clustering.
#'   Method `"random"` places its seeds in `boundary` and does not read it;
#'   it is ignored there with a warning. Only
#'   the first two coordinate columns are clustered, so a Z or M dimension does
#'   not join the distance calculation and dominate it; rows without a
#'   location (an empty point, or a coordinate that is not finite or lies
#'   beyond 1e150, read as given before any move) are dropped with a warning,
#'   so they never reach `stats::kmeans()`, which fails on them without
#'   naming a cause. A lon/lat
#'   cloud, or one with no CRS whose coordinates look like lon/lat (the
#'   heuristic [ensure_projected()] applies, with its warning), is projected
#'   before clustering. A cloud that mixes XY and XYZ (or XYM) points is read
#'   without Z and M before it is projected or clustered, since sf reads a
#'   layer in the dimension of its first geometry.
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
#'   [stats::kmeans()] from 10 random starts on the rows as they come, with
#'   the coordinates centred, divided by the diagonal of their bounding box
#'   and rounded to 2^-30 of it, so that the seeds do not move with the unit
#'   or the origin; on gridded locations they can differ from those of
#'   earlier versions, and where the locations span less than a kilometre a
#'   change of unit or origin can still move them (see "The partition" in
#'   [resolution_profile()]). The
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
#'   and is absent when no `sample_points` were given. A point the same
#'   distance from two seeds (squared distances within one part in ten
#'   million) goes to the lower `seed_id`, the rule the profile scores its
#'   cells by, so the record is the same in another unit or with another
#'   origin. Under `partition = "density"` `withinss` and `tot_withinss`
#'   count every point once, about its cell's mean. A run this call makes with `partition = "density"` or
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
#'   or a profile built with an `estimate` or read by `cp`, the seeds carry
#'   that record as
#'   `attr(, "count_rule")`: how the count was set, for which cell value, and
#'   what the profile scored at it. Its `scored_cells` is `TRUE` for the
#'   scored centres `method = "kmeans"` hands on, and `FALSE` for centres
#'   refitted to every point (see `n`) and for seeds this
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
  # An argument the call passes and the method does not read is reported, as
  # `partition`, `init`, `n` and the density settings are above and below:
  # `seeds` under method = "kmeans" or "random", and `sample_points` under
  # method = "random", were set aside with nothing said.  `seed = 1` matches
  # `seeds` by its prefix (the seed argument is `set_seed`), so such a call
  # drew from the session's stream in silence.  A call the method refuses
  # below is told that, not this: a random seeding without a boundary, and a
  # k-means seeding with neither a boundary nor `sample_points` nor scored
  # centres as `n`.
  refused_below <- is.null(boundary) &&
    (identical(method, "random") ||
       (identical(method, "kmeans") && is.null(sample_points) && is.null(n_res$seeds)))
  if (!identical(method, "provided") && !refused_below) {
    if (!is.null(seeds))
      .warn_and_log(paste0("get_voronoi_seeds(): `seeds` is read by method = \"provided\" ",
                           "only; method = \"%s\" ignores it."), method)
    if (identical(method, "random") && !is.null(sample_points))
      .warn_and_log(paste0("get_voronoi_seeds(): `sample_points` is not read by method = ",
                           "\"random\", which places its seeds in `boundary`; it is ignored."))
  }
  # A layer that mixes XY and XYZ (or XYM) geometries is read in XY
  # (.xy_if_mixed()), before any union, transform or read of its
  # coordinates: sf reads a layer in the dimension of its first geometry, so
  # the cloud was clustered on garbled coordinates or stopped, and GEOS or s2
  # stopped on the boundary or sampled it on other rings.  So is a boundary
  # whose polygons carry M (.xy_if_mixed_or_m()), on which GEOS stopped.  The
  # seeds are new points and the cover is this function's copy, so nothing
  # is said.  Sample points in an sfc_GEOMETRY column are read as the same
  # points in sfc_POINT (.typed_points()): sf::st_coordinates() refuses that
  # class.
  if (!is.null(boundary)) boundary <- .xy_if_mixed_or_m(boundary)
  if (!is.null(sample_points))
    sample_points <- .xy_if_mixed(.typed_points(sample_points))
  # A sample point with a coordinate that is not finite or lies beyond 1e150
  # has no location.  It is read as given, before any move, and emptied
  # (.empty_unlocated()), so that it is dropped below as an empty point is:
  # moved into the boundary's CRS, a northing of 1e200 in EPSG:32616 came
  # back finite in EPSG:32617 and was clustered, and in the working CRS a
  # coordinate of 1e200 left k-means with no partition.
  if (!is.null(sample_points)) sample_points <- .empty_unlocated(sample_points)

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
    refit <- attr(out, "refit", exact = TRUE)
    rule <- .count_rule_get(out)
    out <- .align_crs(out, boundary, caller = "get_voronoi_seeds", arg = "n",
                      label = "the seeds of `n`")
    # Kept across the CRS alignment whatever sf does with custom attributes.
    attr(out, "cover") <- cov
    if (!is.null(rec)) attr(out, "partition") <- rec
    if (!is.null(refit)) attr(out, "refit") <- refit
    attr(out, "n_from") <- n_from
    return(.count_rule_put(out, rule))
  }

  # No scored centres: `density_power` or `density_k` set under partition =
  # "kmeans" is read by nothing in this call.  A k-means seeding is told so
  # in the words of the kind.  A random or provided seeding reads neither
  # setting whatever the kind, so it is told in the words of the method: the
  # kind's words said that the density kind would have read the setting,
  # which a random or provided seeding does not build.  The settings named
  # are the ones each kind reported before: under the k-means kind those not
  # at their defaults, under the density kind as this function's own those
  # the call supplied (with the density kind named, the warning about
  # `partition` and `init` above is the one said).
  if (identical(method, "kmeans")) {
    .warn_density_unused(part, "get_voronoi_seeds")
  } else {
    unread <- if (identical(part$kind, "kmeans")) part$unused
              else if (!partition_given) density_given
              else character(0)
    if (length(unread))
      .warn_and_log(paste0("get_voronoi_seeds(): %s %s read by method = \"kmeans\" only; ",
                           "method = \"%s\" ignores %s."),
                    paste(sprintf("`%s`", unread), collapse = " and "),
                    if (length(unread) > 1L) "are" else "is", method,
                    if (length(unread) > 1L) "them" else "it")
  }

  # --- RNG handling (always restored) ---
  cleanup <- .with_seed(set_seed)
  on.exit(cleanup(), add = TRUE)

  # Every row repaired, then unioned, on the sphere for a lon/lat boundary
  # whatever sf_use_s2() says (.boundary_union()).  In degrees
  # .safe_make_valid() is s2's rebuild, which leaves crossing edges as they
  # are, and the union of a bowtie stopped with s2's bare "Loop 0 is not
  # valid: Edge 0 crosses edge 2"; the rows s2 still refuses are rebuilt with
  # those edges split, as GEOS does in a projected CRS, and a row s2 refuses
  # even then stops the call with an error that names it.
  task <- if (identical(method, "kmeans")) "the sampling of the k-means cloud inside it"
          else "the sampling of random seeds inside it"
  nothing <- if (identical(method, "kmeans")) "no point of the k-means cloud can be drawn in it"
             else "no seed can be drawn in it"
  boundary_union <- function(b)
    .boundary_union(b, caller = "get_voronoi_seeds", task = task, all = TRUE,
                    nothing = nothing)
  # The rows of `boundary` sampled below are those with a location, read as
  # given before the repair (.boundary_located()): a row that is empty or
  # has a vertex at Inf or beyond 1e150 covers no ground, and one with a
  # latitude beyond 90 is refused.  An empty row came out of the repair as
  # an empty line, and the union of nothing stopped the sampling with R's
  # "no applicable method for 'st_geometry'"; a row with a vertex at 1e200
  # took every seed.  The seeds' cover is these rows.
  located <- function(b) .boundary_located(b, "get_voronoi_seeds", nothing)

  out <- switch(
    method,

    "provided" = {
      s <- seeds
      # The seeds are returned in the boundary's CRS (.align_crs() below).
      # Moved there, a layer that mixes XY and XYZ (or XYM) points stopped
      # sf::st_transform() with "OGR error" when its first point had Z, so
      # the transformed copy is returned in XY with a warning of class
      # "spatialkit_zm_dropped" (.xy_if_mixed_warn()).  Seeds that are not
      # moved come back as they were given, Z included.
      if (!is.null(boundary) && !is.na(sf::st_crs(s)) && !is.na(sf::st_crs(boundary)) &&
          sf::st_crs(s) != sf::st_crs(boundary))
        s <- .xy_if_mixed_warn(s, "get_voronoi_seeds", "seeds")
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
      boundary <- located(boundary)
      b <- boundary_union(boundary)
      pts <- .robust_st_sample(b, n, "get_voronoi_seeds", nothing,
                                remedy = .sample_remedy("random"))
      pts_sfc <- sf::st_sfc(pts, crs = sf::st_crs(boundary))
      sf::st_sf(seed_id = seq_along(pts_sfc), method = "random",
                geometry = pts_sfc)
    },

    "kmeans" = {
      if (!is.null(sample_points)) {
        cloud <- if (!is.null(boundary))
          .align_crs(sample_points, boundary, caller = "get_voronoi_seeds",
                     arg = "sample_points")
        else sample_points
      } else {
        if (is.null(boundary))
          stop("get_voronoi_seeds(): k-means seeds require 'boundary' (or provide 'sample_points').",
               call. = FALSE)
        boundary <- located(boundary)
        b <- boundary_union(boundary)
        cloud_n <- .kmeans_cloud_size(n)
        cloud_geom <- .robust_st_sample(b, cloud_n, "get_voronoi_seeds", nothing,
                                        remedy = .sample_remedy("kmeans"))
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
      # warning, so do the same here.  A coordinate beyond 1e150 has no
      # location either (.values_at_xy()), read again after the projection.
      bad_xy <- is.na(.values_at_xy(cloud_for_km)[, 1L])
      if (any(bad_xy)) {
        .warn_and_log("get_voronoi_seeds(kmeans): dropping %d point(s) without a location (%s).",
                      sum(bad_xy), .UNLOCATED_WHAT)
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
      # The distinct positions are counted as stats::kmeans() is handed them
      # (.n_locations()), so the count does not move with the unit.
      n_uniq <- .n_locations(xy)
      k_max  <- min(n_uniq, nrow(xy) - 1L)
      k_use  <- max(1L, min(as.integer(n), k_max))
      if (k_use < n) {
        # The request as given, not as.integer(n) (see voronoi_seeds_kmeans()).
        .warn_and_log("get_voronoi_seeds(kmeans): requested %s seeds but the sampling cloud supports at most %d (%d unique position(s) among %d point(s)); clamping.",
                  format(n, scientific = FALSE), k_use, n_uniq, nrow(xy))
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
  # Only provided seeds can be moved here: the other branches build them in
  # the boundary's CRS.
  out <- .align_crs(out, boundary, caller = "get_voronoi_seeds", arg = "seeds")
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
#' The box the draw is made in is read from the coordinates first.  With no
#' width or no height (a lon/lat ring of three vertices on one parallel,
#' which s2 reads as a sliver with an area) \code{sf::st_sample()} drew in
#' it without end, so the call stops in the caller's name instead.
#'
#' \code{sf::st_sample()} draws \code{n} times the ratio of the bounding box
#' to the area of points in the box and keeps those inside, at about 470
#' bytes of memory per point drawn.  The size of that draw is read first
#' (\code{.sample_draw_plan()}), without touching the random stream.  Where
#' it is at most \code{.sample_draw_cap(n)} points (10^6, or 10 \code{n} if
#' larger), the geometry is sampled whole, as before, to the bit.  Above
#' that, as for small pieces far apart or a thin strip at a slant, in any
#' CRS, the geometry is sampled piece by piece (\code{.sample_by_piece()}),
#' and a piece whose own draw is still above the cap is drawn in successive
#' draws of at most 10^6 points in its box (\code{.sample_bounded()}).
#' Where the whole takes fewer points in its box than its pieces would in
#' theirs, with each call of sf counted as 5000 points, the whole is drawn
#' in those bounded draws instead (many small pieces whose draw is just
#' above the cap).  sf's draw of 2.8e9 points stopped the call with R's
#' "vector memory limit ... reached", which names no cause.
#'
#' In longitude and latitude sf refuses to sample a region whose bounding
#' box is over 10^4 times its area ("sampling box is N times larger than
#' sampling region").  Only that refusal is caught: the geometry is then
#' sampled piece by piece as above.  A piece sf still refuses stops the call
#' with an error in the caller's name that gives the cause and says to
#' transform `boundary` to a projected CRS.
#'
#' @param geom An sfc geometry to sample within.
#' @param n Integer; desired number of sample points.
#' @param caller Name of the exported function, for the errors.
#' @param nothing What cannot be drawn ("no seed can be drawn in it").
#' @param by_piece \code{FALSE} for a piece of the piece-by-piece draw,
#'   which stops where sf refuses it and is drawn in bounded draws where its
#'   own draw is above the cap.
#' @param remedy The sentence that ends the refusal of a draw above
#'   \code{.sample_total_cap()}: what the caller can change.
#' @param plan The draw plan of \code{geom} for \code{n}, where the caller
#'   has it (a piece of \code{.sample_by_piece()}); \code{NULL} to read it
#'   here.
#' @return An sfc_POINT of length n.
#' @keywords internal
#' @noRd
.robust_st_sample <- function(geom, n, caller = "spatialkit",
                              nothing = "no point can be drawn in it", by_piece = TRUE,
                              remedy = .sample_remedy("random"), plan = NULL) {
  # st_sample() sizes its draw from st_area(), which on lon/lat input needs
  # lwgeom (not a dependency) when sf_use_s2() is FALSE: random and k-means
  # seeding on a lon/lat boundary failed with "package lwgeom required".
  if (.is_longlat(geom) && !isTRUE(sf::sf_use_s2()))
    return(.with_s2(.robust_st_sample(geom, n, caller, nothing, by_piece, remedy, plan)))
  .sample_box_check(geom, caller, nothing)
  # The draw sf makes is bounded: above the cap the geometry goes piece by
  # piece, and a piece in bounded draws.  A lon/lat region sf refuses is left
  # to sf's refusal below, which draws nothing.
  if (is.null(plan)) plan <- .sample_draw_plan(geom, n)
  if (isTRUE(plan$draws > .sample_draw_cap(n)) && !any(plan$sf_refuses))
    return(if (isTRUE(by_piece))
             .sample_by_piece(geom, n, NULL, caller, nothing, by_piece, remedy, whole = plan)
           else .sample_bounded(geom, n, plan$ratio, caller))
  # sf's refusal of a lon/lat region whose box is over 10^4 times its area.
  # It stops before any number is drawn, so the stream is where it was.
  box_refused <- function(e)
    grepl("sampling box is .* times larger than sampling region", conditionMessage(e))
  pts <- try(.st_sample_quiet(geom, size = n, type = "random", exact = TRUE),
             silent = TRUE)
  if (inherits(pts, "try-error")) {
    pts <- tryCatch(.st_sample_quiet(geom, size = n, type = "random"),
                    error = function(e) if (box_refused(e)) e else stop(e))
    if (inherits(pts, "error"))
      return(.sample_by_piece(geom, n, pts, caller, nothing, by_piece, remedy))
  }
  # Pad or trim to exactly n
  max_attempts <- 10L
  attempt <- 0L
  while (length(pts) < n && attempt < max_attempts) {
    attempt <- attempt + 1L
    extra <- .st_sample_quiet(geom, size = n - length(pts), type = "random")
    pts <- c(pts, extra)
  }
  if (length(pts) > n) pts <- pts[seq_len(n)]
  if (length(pts) < n)
    .log_warn(".robust_st_sample(): could only generate %d of %d requested points.",
              length(pts), n)
  pts
}

#' sf::st_sample() without sf's lwgeom warning
#'
#' Without lwgeom (not a dependency) sf warns "coordinate ranges not
#' computed along great circles; install package lwgeom to get rid of this
#' warning" on every lon/lat draw, so every random or k-means seeding on a
#' lon/lat boundary raised it, once or twice, in the ordinary case.  Only
#' that warning is muffled; the draw is the same.
#'
#' @param ... Passed to \code{sf::st_sample()}.
#' @return What \code{sf::st_sample()} returns.
#' @keywords internal
#' @noRd
.st_sample_quiet <- function(...) withCallingHandlers(
  sf::st_sample(...),
  warning = function(w)
    if (grepl("coordinate ranges not computed along great circles",
              conditionMessage(w), fixed = TRUE))
      invokeRestart("muffleWarning"))

#' Points sf may draw in a bounding box in one draw, and in all
#'
#' One draw of 10^6 points in the box takes about 0.5 GB of points (1.3 GB
#' at peak) and two to four seconds; 2e7 in all take about a minute.  The
#' 10 \code{n} keeps a large \code{n} on a region that fills its box as it
#' was.
#'
#' @param n Number of points asked for.
#' @return A single number.
#' @keywords internal
#' @noRd
.sample_draw_cap  <- function(n) max(1e6, 10 * as.numeric(n))
.sample_total_cap <- function(n) max(2e7, 10 * as.numeric(n))

#' The cost of a draw, in points drawn in the box
#'
#' One call of \code{sf::st_sample()} costs about as much time as 5000
#' points drawn in the box (16 ms against 3.5 microseconds a point,
#' measured on a projected square).  \code{.sample_by_piece()} compares
#' the whole with its pieces by this.
#'
#' @param points Points drawn in the box.
#' @param calls Calls of \code{sf::st_sample()}.
#' @return A single number.
#' @keywords internal
#' @noRd
.sample_cost <- function(points, calls) as.numeric(points) + 5000 * as.numeric(calls)

#' The last sentence of the refusal of a draw above the cap
#'
#' @param what \code{"kmeans"} for the k-means cloud of
#'   \code{get_voronoi_seeds()}, \code{"random"} for the seeds of
#'   \code{get_voronoi_seeds()}, \code{"k"} for those of
#'   \code{voronoi_seeds_random()}.
#' @return A sentence.
#' @keywords internal
#' @noRd
.sample_remedy <- function(what) {
  fill <- "give a boundary whose polygons fill more of their bounding boxes"
  switch(what,
    kmeans = paste0("Pass `sample_points` (points inside `boundary`, clustered in ",
                    "place of a cloud drawn in it), or ", fill, "."),
    random = paste0("Ask for fewer seeds (`n`), or ", fill, "."),
    k      = paste0("Ask for fewer seeds (`k`), or ", fill, "."))
}

#' The draw sf::st_sample() makes to sample a region
#'
#' sf draws \code{round(n * a0 / a1)} points in the bounding box and keeps
#' those inside.  Here a0 and a1 are as sf takes them: in a projected CRS
#' (or none) the area of the box over the planar area, with the ratio 1
#' where either is not finite and positive (sf then draws \code{n});
#' in longitude and latitude the area of the band between the box's
#' parallels and meridians on the sphere over the s2 area, and sf refuses
#' the region when the ratio is over 10^4 (or not a number).  Nothing is
#' drawn from the random stream.
#'
#' @param geom An sfc.
#' @param n Number of points asked for.
#' @param each \code{TRUE} to read each element of \code{geom} on its own
#'   (the pieces of \code{.sample_by_piece()}); \code{FALSE} for the whole.
#' @return A list: \code{longlat}; \code{a0}, \code{a1} and \code{ratio}
#'   (one per element with \code{each}, else one); \code{draws}, \code{n}
#'   times the ratio (the whole only, else \code{NA}); \code{sf_refuses}, per
#'   ratio, \code{TRUE} where sf refuses the lon/lat region.
#' @keywords internal
#' @noRd
.sample_draw_plan <- function(geom, n, each = FALSE) {
  geom <- sf::st_geometry(geom)
  longlat <- isTRUE(sf::st_is_longlat(geom))
  bbs <- if (isTRUE(each)) lapply(geom, sf::st_bbox) else list(sf::st_bbox(geom))
  if (longlat) {
    # sf's own expression, so that the ratio, and its refusal over 10^4,
    # match sf's to the bit; R is s2::s2_earth_radius_meters().  With lwgeom
    # (not a dependency) sf takes the box of the geometry segmentized along
    # great circles instead, which this does not follow.
    R <- 6371010
    toRad <- pi / 180
    a0 <- vapply(bbs, function(bb) {
      h1 <- sin(bb["ymax"] * toRad)
      h2 <- sin(bb["ymin"] * toRad)
      as.numeric(2 * pi * R^2 * (h1 - h2) * (bb["xmax"] - bb["xmin"]) / 360)
    }, numeric(1))
    a1 <- as.numeric(.with_s2(sf::st_area(geom)))
    if (!isTRUE(each)) a1 <- sum(a1)
    ratio <- a0 / a1
    sf_refuses <- !(ratio <= 1e4)
  } else {
    # The box polygons as sf::st_as_sfc() builds them, all in one sfc, so
    # their areas are sf's and are read in one call.
    boxes <- lapply(bbs, function(bb)
      sf::st_polygon(list(matrix(as.double(bb)[c(1, 2, 3, 2, 3, 4, 1, 4, 1, 2)],
                                 ncol = 2, byrow = TRUE))))
    a0 <- as.numeric(sf::st_area(sf::st_sfc(boxes, crs = sf::st_crs(geom))))
    a1 <- as.numeric(sf::st_area(geom))
    if (!isTRUE(each)) a1 <- sum(a1)
    ratio <- ifelse(is.finite(a0) & is.finite(a1) & a0 > 0 & a1 > 0, a0 / a1, 1)
    sf_refuses <- rep(FALSE, length(ratio))
  }
  list(longlat = longlat, a0 = a0, a1 = a1, ratio = ratio,
       draws = if (isTRUE(each)) NA_real_ else as.numeric(n) * ratio,
       sf_refuses = sf_refuses)
}

#' Refuse a region whose bounding box has no width or no height
#'
#' \code{sf::st_sample()} draws in the bounding box of the region, read from
#' its coordinates, and keeps the points inside.  A lon/lat ring of three
#' vertices on one parallel is valid under s2 (its edges are great circles,
#' so it holds a sliver of ground) and its box has no height: no point drawn
#' in the box falls inside, and the draw ran without end.
#' create_grid_polygons() refuses such a box ("non-positive extent").
#'
#' @param geom An sfc.
#' @param caller,nothing As \code{.robust_st_sample()}.
#' @return \code{NULL}, invisibly, when the box has a width and a height;
#'   otherwise it stops.
#' @keywords internal
#' @noRd
.sample_box_check <- function(geom, caller, nothing) {
  bb <- sf::st_bbox(geom)
  w  <- as.numeric(bb[["xmax"]] - bb[["xmin"]])
  h  <- as.numeric(bb[["ymax"]] - bb[["ymin"]])
  if (!(is.finite(h) && h > 0) || !(is.finite(w) && w > 0))
    stop(sprintf(paste0("%s(): `boundary` cannot be sampled (its bounding box, read from ",
                        "its coordinates, has no %s), so %s."),
                 caller, if (!(is.finite(h) && h > 0)) "height" else "width", nothing),
         call. = FALSE)
  invisible(NULL)
}

#' Sample a region piece by piece
#'
#' See \code{.robust_st_sample()}.  It runs where sf refuses a lon/lat
#' region whole (\code{parent} is that refusal), and where sf's draw for the
#' whole region is above \code{.sample_draw_cap()} (\code{parent} is
#' \code{NULL}), in any CRS.  The geometry is cast to its polygons, and
#' \code{n} is shared among them by their areas (on the sphere in
#' longitude and latitude, planar otherwise) in one multinomial draw from
#' the stream, so the points stay uniform over the whole; a single polygon
#' takes all \code{n} with no draw.  Each piece is then sampled by
#' \code{.robust_st_sample()}, in the order of the pieces.
#'
#' Before any draw, the points to be drawn in the boxes of the pieces are
#' added up from their areas: \code{n} times the sum of the box areas over
#' the area of the whole.  Above \code{.sample_total_cap(n)} the call stops
#' in the caller's name with the box-to-area ratio, the draws needed and
#' \code{remedy}.  A lon/lat piece sf refuses is left out of that sum: if
#' it is drawn, it stops the call with sf's refusal, as before.
#'
#' Where the draw for the whole is above the cap (\code{whole} given), the
#' whole is drawn in bounded draws (\code{.sample_bounded()}) instead of
#' piece by piece when its points in the box are within
#' \code{.sample_total_cap(n)} and cost less (\code{.sample_cost()}) than
#' the pieces' points and one call of sf for each piece that takes a point,
#' or when the pieces' points are above that cap: many small pieces whose
#' whole draw is just above the cap take a few draws, not one call each.
#'
#' @param geom,n,caller,nothing,by_piece,remedy As
#'   \code{.robust_st_sample()}.
#' @param parent sf's refusal of the box of \code{geom}, or \code{NULL}.
#' @param whole \code{.sample_draw_plan(geom, n)}, where \code{parent} is
#'   \code{NULL}.
#' @return An sfc_POINT of length n.
#' @keywords internal
#' @noRd
.sample_by_piece <- function(geom, n, parent, caller, nothing, by_piece,
                             remedy = .sample_remedy("random"), whole = NULL) {
  pieces <- if (isTRUE(by_piece)) suppressWarnings(sf::st_cast(
    sf::st_collection_extract(geom, "POLYGON", warn = FALSE), "POLYGON"))
  if (!is.null(parent) && (!isTRUE(by_piece) || length(pieces) < 2L)) {
    said <- sub(";$", "", trimws(strsplit(conditionMessage(parent), "\n",
                                         fixed = TRUE)[[1L]][1L]))
    stop(sprintf(paste0(
      "%s(): sf refuses to sample `boundary` in longitude and latitude: the bounding ",
      "box of %s is over 10,000 times %s area (sf: \"%s\"), so %s. Transform `boundary` ",
      "to a projected CRS, such as the one ensure_projected() picks for it ",
      "(sf::st_transform()), and call %s() again."),
      caller, if (isTRUE(by_piece)) "`boundary`" else "one of its polygons",
      if (isTRUE(by_piece)) "its" else "that polygon's", said, nothing, caller),
      call. = FALSE)
  }
  plan <- .sample_draw_plan(pieces, n, each = TRUE)
  area <- plan$a1
  if (is.null(parent) && !is.null(whole) &&
      .sample_whole_bounded(whole, plan, n))
    return(.sample_bounded(geom, n, whole$ratio, caller))
  .sample_total_check(plan, n, length(pieces), caller, nothing, remedy)
  take <- if (length(pieces) == 1L) as.numeric(n)
          else as.vector(stats::rmultinom(1L, n, area / sum(area)))
  # Each piece's plan is read from the plan of all pieces.
  got  <- lapply(which(take > 0L), function(i)
    .robust_st_sample(pieces[i], take[i], caller, nothing, by_piece = FALSE,
                      remedy = remedy,
                      plan = list(longlat = plan$longlat, a0 = plan$a0[i], a1 = plan$a1[i],
                                  ratio = plan$ratio[i],
                                  draws = as.numeric(take[i]) * plan$ratio[i],
                                  sf_refuses = plan$sf_refuses[i])))
  do.call(c, got)
}

#' Points the pieces of a region take in their boxes, in expectation
#'
#' Each piece takes \code{n a1 / sum(a1)} points in expectation and draws
#' its ratio times that in its box.  A lon/lat piece sf refuses is left out.
#'
#' @param plan \code{.sample_draw_plan(pieces, n, each = TRUE)}.
#' @param n Number of points asked for.
#' @return A single number; \code{NA} where the area is not finite and
#'   positive.
#' @keywords internal
#' @noRd
.sample_pieces_need <- function(plan, n) {
  total_area <- sum(plan$a1)
  if (!(is.finite(total_area) && total_area > 0)) return(NA_real_)
  drawn <- !plan$sf_refuses
  as.numeric(n) * sum((plan$a1 * plan$ratio)[drawn]) / total_area
}

#' Whether to draw the whole in bounded draws rather than piece by piece
#'
#' See \code{.sample_by_piece()}.
#'
#' @param whole \code{.sample_draw_plan(geom, n)}.
#' @param plan \code{.sample_draw_plan(pieces, n, each = TRUE)}.
#' @param n Number of points asked for.
#' @return \code{TRUE} or \code{FALSE}.
#' @keywords internal
#' @noRd
.sample_whole_bounded <- function(whole, plan, n) {
  if (isTRUE(any(whole$sf_refuses)) || !isTRUE(whole$draws <= .sample_total_cap(n)))
    return(FALSE)
  need <- .sample_pieces_need(plan, n)
  if (!isTRUE(need <= .sample_total_cap(n))) return(TRUE)
  calls <- min(length(plan$ratio), as.numeric(n))
  .sample_cost(whole$draws, ceiling(whole$draws / 1e6)) <= .sample_cost(need, calls)
}

#' Refuse a piece-by-piece draw whose points in all are above the cap
#'
#' See \code{.sample_by_piece()}.
#'
#' @param plan \code{.sample_draw_plan(pieces, n, each = TRUE)}.
#' @param n Number of points asked for.
#' @param n_pieces Number of pieces.
#' @param caller,nothing,remedy As \code{.robust_st_sample()}.
#' @return \code{NULL}, invisibly, at or under the cap; otherwise it stops.
#' @keywords internal
#' @noRd
.sample_total_check <- function(plan, n, n_pieces, caller, nothing, remedy) {
  need <- .sample_pieces_need(plan, n)
  cap  <- .sample_total_cap(n)
  if (!isTRUE(need > cap)) return(invisible(NULL))
  drawn <- !plan$sf_refuses
  big <- which(drawn)[which.max((plan$a1 * plan$ratio)[drawn])]
  num <- function(x) format(signif(x, 3), big.mark = ",", scientific = FALSE, trim = TRUE)
  one <- n_pieces == 1L
  stop(sprintf(paste0(
    "%s(): `boundary` cannot be sampled: the bounding box of %s is %s times %s area, ",
    "so drawing %s points inside %s takes about %s points drawn in %s, more than the ",
    "%s drawn at most, so %s. %s"),
    caller, if (one) "`boundary`" else "one of its polygons", num(plan$ratio[big]),
    if (one) "its" else "that polygon's", num(n), if (one) "it" else "`boundary`",
    num(need), if (one) "the box" else "the boxes of its polygons", num(cap),
    nothing, remedy), call. = FALSE)
}

#' Sample a region in draws of at most 10^6 points in its box
#'
#' sf draws \code{size} times the box-to-area ratio of points in the box, so
#' each draw asks for at most 10^6 / ratio points (a fraction where the
#' ratio is above 10^6, which sf rounds to that many points in the box),
#' with \code{exact = FALSE}; the points are kept in the order drawn until
#' \code{n} are taken.  \code{.sample_by_piece()} has checked the points in
#' all against \code{.sample_total_cap()} first.  A ratio that is not
#' finite and at least 1 is refused: a region with no area would draw
#' without end.
#'
#' @param geom An sfc of one polygon, or the region whole.
#' @param n Number of points.
#' @param ratio Its box-to-area ratio, from \code{.sample_draw_plan()}.
#' @param caller Name for the log line.
#' @return An sfc_POINT of length n (fewer, with a log line, if the draws
#'   fall short).
#' @keywords internal
#' @noRd
.sample_bounded <- function(geom, n, ratio, caller) {
  if (!(length(ratio) == 1L && is.finite(ratio) && ratio >= 1))
    stop(sprintf("%s(): internal error: a bounded draw needs a box-to-area ratio of at least 1, not %s.",
                 caller, format(ratio)), call. = FALSE)
  per_draw <- 1e6 / ratio
  got <- list()
  have <- 0
  calls <- 0
  max_calls <- 4 * ceiling(as.numeric(n) / per_draw) + 10
  while (have < n && calls < max_calls) {
    calls <- calls + 1
    p <- .st_sample_quiet(geom, size = min(as.numeric(n) - have, per_draw),
                          type = "random", exact = FALSE, warn_if_not_integer = FALSE)
    got[[length(got) + 1L]] <- p
    have <- have + length(p)
  }
  pts <- do.call(c, got)
  if (length(pts) > n) pts <- pts[seq_len(n)]
  if (length(pts) < n)
    .log_warn("%s(): could only generate %d of %d requested points.",
              caller, length(pts), n)
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
#' come back in the input's own coordinates. Rows without a location (an
#' empty point, or a coordinate that is not finite or lies beyond 1e150, read
#' as given before the projection) are dropped with a warning, and `k` is
#' clamped to the number of distinct positions or one fewer than the number
#' of points, whichever is smaller, and not below 1.
#'
#' A count given as a number gets a fresh partition, which is not the
#' partition [resolution_profile()] or [determine_optimal_levels()] scored
#' that count on: those fit their own run, the profile to every row or to a
#' subsample, and keep the best of their own starts. Pass the
#' [select_resolution()] result or the profile itself as `k` to get the
#' scored one: the centres of the scored partition are then returned (in
#' the CRS of `points_sf`, with `attr(, "n_from")`), and no k-means is run.
#' On a profile without an estimate whose centres were fitted to a
#' subsample, they are those centres refitted to every point by
#' [select_resolution()], with `attr(, "refit")`.
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
#' @param points_sf An sf object with POINT geometries. Points that mix XY
#'   and XYZ (or XYM) are read without Z and M before they are projected or
#'   clustered, since sf reads a layer in the dimension of its first geometry.
#' @param k Integer; requested number of clusters, treated as an upper bound
#'   only. It is clamped, with a warning, to whichever is smaller of the number
#'   of distinct point positions and one fewer than the number of points with
#'   a location (rows without one are dropped first), and not below 1, because
#'   k-means can produce neither more centres than there are distinct points
#'   nor as many centres as there are rows. Check `nrow()` on the result. A
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
#'   rows as they come, with the coordinates centred, divided by the
#'   diagonal of their bounding box and rounded to 2^-30 of it, so that the
#'   seeds do not move with the unit or the origin; on gridded locations
#'   they can differ from those of earlier versions, and where the locations
#'   span less than a kilometre a change of unit or origin can still move
#'   them (see "The partition" in [resolution_profile()]). Without the 'FNN' package the density kind
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
#'   `k` is clamped, see `k`), with `seed_id` and
#'   `method = "kmeans"` columns matching [get_voronoi_seeds()], and
#'   `attr(, "cover")`, the convex hull of the coordinates of `points_sf`
#'   (an sfc in its CRS): for lon/lat points, the hull in degrees, with the
#'   longitudes of points on both sides of 180 degrees taken from 0 to 360
#'   (-170 as 190).
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
#'   built with an `estimate` or read by `cp`, the seeds carry that record as
#'   `attr(, "count_rule")`, with `scored_cells` `TRUE` for the scored
#'   centres and `FALSE` for refitted centres and for centres this call fits
#'   itself at that count (see [get_voronoi_seeds()]).
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
  # Points that mix XY and XYZ (or XYM) are read in XY (.xy_if_mixed()),
  # before the projection and the read of their coordinates, as in
  # get_voronoi_seeds(); the seeds are new points, so nothing is said.  So
  # are points in an sfc_GEOMETRY column, as the same points in sfc_POINT
  # (.typed_points()): sf::st_coordinates() refuses that class.
  points_sf <- .xy_if_mixed(.typed_points(points_sf))
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
  # The NA determine_optimal_levels() returns where the points name no
  # count is refused with the note get_voronoi_seeds() gives it.
  if (is.numeric(k) && length(k) == 1L && is.na(k) && nzchar(.no_count_note(k)))
    stop(sprintf("voronoi_seeds_kmeans(): `k` must be a single positive number of seeds; got NA.%s",
                 .no_count_note(k)), call. = FALSE)
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
  # A point with a coordinate that is not finite or lies beyond 1e150 has no
  # location.  It is read as given, before the projection, and emptied
  # (.empty_unlocated()), so that it takes no part in the lon/lat test or the
  # choice of the CRS and is dropped below as an empty point is: a northing
  # of 1e200 left k-means with no partition.
  points_sf  <- .empty_unlocated(points_sf)
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
  # kmeans(), which fails with "NA/NaN/Inf in foreign function call".  A
  # coordinate beyond 1e150 has no location either (.values_at_xy()), read
  # again after the projection.
  finite_rows <- !is.na(.values_at_xy(pts_for_km)[, 1L])
  n_drop <- sum(!finite_rows)
  if (n_drop > 0L) {
    .warn_and_log("voronoi_seeds_kmeans(): dropping %d point(s) without a location (%s).",
              n_drop, .UNLOCATED_WHAT)
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
  n_uniq <- .n_locations(coords)
  k_max  <- min(n_uniq, n - 1L)
  k_use  <- max(1L, min(as.integer(k), k_max))
  if (k_use < k) {
    # The request as given: as.integer() read k = 9.5 as "requested 9 seeds",
    # a request equal to the bound it was clamped to.
    .warn_and_log("voronoi_seeds_kmeans(): requested %s seeds but the points support at most %d (%d unique position(s) among %d point(s)); clamping.",
              format(k, scientific = FALSE), k_use, n_uniq, n)
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
#' @param boundary An sf or sfc polygonal object. One that mixes XY and XYZ
#'   (or XYM) geometries is read without Z and M, since sf reads a layer in
#'   the dimension of its first geometry. So is one whose
#'   geometries carry M (XYM or XYZM), which GEOS refuses. Rows that GEOS
#'   (in a projected CRS) or s2 (in longitude and latitude) refuses are
#'   repaired before the rows are unioned, as in [get_voronoi_seeds()];
#'   valid rows are used as they are.
#'   A row that is empty, or has a coordinate that is not finite or lies
#'   beyond 1e150 (read as given, and again after any move into another
#'   CRS), covers no ground and is left out; a boundary with no such row
#'   left is refused with an error that says so. A boundary in a geographic
#'   CRS with a vertex at a latitude beyond 90 degrees, which is no place on
#'   the Earth, is refused with an error of class
#'   `"spatialkit_cells_degenerate"` that names its rows and says to check
#'   its CRS.
#'   A row the repair leaves without a polygon (a ring with no height:
#'   every vertex on one line, or in degrees on one parallel) is left out
#'   as well, and a boundary with no polygon left is refused. In
#'   longitude and latitude, a boundary whose bounding box is over 10,000
#'   times its area (small pieces far apart), which sf does not sample whole,
#'   is sampled piece by piece, the points shared among the pieces in
#'   proportion to their areas; a piece whose own box is over 10,000 times
#'   its area stops the call with an error that says to transform
#'   `boundary` to a projected CRS. A boundary whose bounding box has no
#'   height or no width (three vertices on one parallel) is refused.
#'   sf draws the points in the bounding box and keeps those inside, so it
#'   draws `k` times the ratio of the box to the area. In any CRS, a
#'   boundary that would take more than 10^6 points drawn in its box, or
#'   10 `k` if more (small pieces far apart, a thin strip at a slant), is
#'   sampled piece by piece in the same way (or whole, in draws of at most
#'   10^6 points in its box, where that takes less time), and a piece that
#'   would still take more is drawn in such draws. Where the points drawn
#'   in all would exceed 2e7 (or 10 `k`, if more), the call stops before
#'   any draw with an error that gives the ratio and says to ask for fewer
#'   seeds.
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
  # The NA determine_optimal_levels() returns where the points name no
  # count is refused with the note get_voronoi_seeds() gives it.
  if (is.numeric(k) && length(k) == 1L && is.na(k) && nzchar(.no_count_note(k)))
    stop(sprintf("voronoi_seeds_random(): `k` must be a single positive number of seeds; got NA.%s",
                 .no_count_note(k)), call. = FALSE)
  .check_scalar(k, "k", "voronoi_seeds_random", min = 1,
                max = .Machine$integer.max,
                what = "a single positive number of seeds")
  # `@param boundary` documents sf *or* sfc, and .assert_sf() only accepts sf.
  if (inherits(boundary, "sfc")) boundary <- sf::st_as_sf(boundary)
  .assert_sf(boundary, c("POLYGON", "MULTIPOLYGON"), "boundary")
  # A boundary that mixes XY and XYZ (or XYM) polygons, or whose polygons
  # carry M, is read in XY (.xy_if_mixed_or_m()), as in get_voronoi_seeds():
  # GEOS stopped on it, and s2 sampled other rings of a mixed one in silence.
  boundary <- .xy_if_mixed_or_m(boundary)

  cleanup <- .with_seed(set_seed)
  on.exit(cleanup(), add = TRUE)
  # On the sphere for lon/lat, as in get_voronoi_seeds(): with s2 off sf
  # otherwise printed its planar-union message on every call.  The rows were
  # unioned as they came, so a boundary get_voronoi_seeds() repairs and
  # samples stopped here with GEOS's or s2's bare message (a vertex listed
  # twice in degrees, a bowtie among valid rows).  The rows GEOS or s2
  # refuses are repaired first (.boundary_union()); valid rows are unioned as
  # they arrived, so their seeds are those they were.
  # Only the rows with a location, read as given before the repair
  # (.boundary_located(), as in get_voronoi_seeds()); a latitude beyond 90
  # is refused.
  task    <- "the sampling of random seeds inside it"
  nothing <- "no seed can be drawn in it"
  boundary <- .boundary_located(boundary, "voronoi_seeds_random", nothing)
  geom <- .boundary_union(boundary, caller = "voronoi_seeds_random", task = task,
                          nothing = nothing)
  pts <- .robust_st_sample(geom, k, "voronoi_seeds_random", nothing,
                           remedy = .sample_remedy("k"))
  out <- sf::st_sf(geometry = pts) |> sf::st_set_crs(sf::st_crs(boundary))
  # Same output contract as get_voronoi_seeds(), so the three seeding
  # functions are drop-in interchangeable.
  out$seed_id <- seq_len(nrow(out))
  out$method  <- "random"
  # The ground the seeds stand for: the area they were drawn in.
  .set_cover(out, boundary)
}


#' Repair the rows of a boundary before they are unioned
#'
#' get_voronoi_seeds(), voronoi_seeds_random(), the key of
#' create_grid_polygons_cached(), make_folds(method = "block_kfold") and a
#' grid laid in degrees by create_grid_polygons() union `boundary` or clip to
#' it.  A row GEOS calls invalid (a bowtie among valid rows) stopped the
#' union with GEOS's bare "TopologyException", and in longitude and latitude,
#' where sf hands the union to s2, a row s2 refuses (a vertex listed twice, a
#' figure eight, crossing edges) stopped it with s2's bare "Loop 0 is not
#' valid".
#'
#' In a projected CRS (or none), a row is repaired with
#' \code{.safe_make_valid()}, which splits crossing edges there.  In
#' longitude and latitude a row s2 refuses is repaired with it too, under s2
#' whatever \code{sf::sf_use_s2()} says, and the rows s2 still refuses are
#' repaired one by one (\code{.s2_mend()}): first without the vertices that
#' repeat the one before them, which moves no area, then rebuilt with
#' crossing edges split (\code{.split_crossing_edges()}), which
#' \code{sf::st_make_valid()} does not do in degrees.  A row s2 refuses even
#' then is passed to \code{refuse}, which stops.
#'
#' @param b The boundary, sf or sfc.
#' @param refuse A function of \code{parent} (s2's error), \code{rows},
#'   \code{reasons}, \code{n} and \code{s2_forced} that raises the caller's
#'   error; \code{NULL} for the error \code{.boundary_refused()} builds,
#'   which names \code{caller}.
#' @param caller Name of the exported function, for that error.
#' @param task What cannot run in degrees, for that error.
#' @param all \code{TRUE} to put every row through \code{.safe_make_valid()}
#'   first, under the engine of the session, as get_voronoi_seeds() and
#'   create_grid_polygons() always have; \code{FALSE} to repair only the rows
#'   the engine of the union refuses (GEOS's validity in a projected CRS,
#'   s2's conversion in degrees), so a boundary of valid rows is unioned as it
#'   arrived, to the bit.
#' @param nothing What cannot be done in `boundary` when the repair leaves no
#'   polygon in it ("no seed can be drawn in it"), for the error
#'   \code{.stop_boundary_unlocated()} raises then.
#' @return \code{b}, of the same class, with the rows repaired.  A row the
#'   repair leaves without a polygon (a ring with no height, which the repair
#'   empties or turns into a line) is left out, when \code{b} had a polygon
#'   before the repair.  With no polygon left the call stops: the union
#'   handed the samplers nothing, and they stopped with R's "no applicable
#'   method for 'st_geometry'".
#' @keywords internal
#' @noRd
.boundary_mend <- function(b, refuse = NULL, caller = NULL, task = NULL, all = FALSE,
                           nothing = "nothing can be built or drawn in it") {
  g  <- sf::st_geometry(b)
  ll <- .is_longlat(g)
  had <- any(.has_polygon(g))
  if (isTRUE(all)) {
    g <- .safe_make_valid(g)
  } else if (!ll && length(g)) {
    bad <- which(!(suppressWarnings(sf::st_is_valid(g)) %in% TRUE))
    if (length(bad)) g[bad] <- .safe_make_valid(g[bad])
  }
  if (ll && length(g)) {
    s2_forced <- !isTRUE(sf::sf_use_s2())
    g <- .with_s2({
      # Without `all`, the rows s2 refuses get the repair the others skip,
      # s2's rebuild, as get_voronoi_seeds() gives every row with s2 on.
      if (!isTRUE(all)) {
        bad <- .s2_refused(g)
        if (length(bad)) g[bad] <- .safe_make_valid(g[bad])
      }
      m <- .s2_mend(g, split = TRUE)
      if (length(m$refused)) {
        sub    <- m$g[m$refused]
        parent <- tryCatch({ suppressMessages(sf::st_as_s2(sub)); NULL },
                           error = function(e) e)
        if (is.null(parent))
          parent <- simpleError("s2 refuses the polygon")
        why <- .s2_refused_rows(sub)
        reasons <- if (!is.null(why) && length(why$rows) == length(m$refused))
          why$reasons else rep(conditionMessage(parent), length(m$refused))
        if (is.null(refuse))
          .boundary_refused(parent, m$refused, reasons, length(g), caller, task, s2_forced)
        refuse(parent = parent, rows = m$refused, reasons = reasons, n = length(g),
               s2_forced = s2_forced)
      }
      m$g
    })
  }
  # A row the repair leaves without a polygon (an empty collection, or a
  # line for a ring with no height) covers no ground and is left out, so a
  # square beside it is unioned as the square alone; with no polygon left
  # the call stops.
  keep <- if (had && length(g)) .has_polygon(g) else rep(TRUE, length(g))
  if (!any(keep))
    .stop_boundary_unlocated(if (is.null(caller)) "spatialkit" else caller, nothing,
                             if (ll) "s2" else "geos")
  if (!inherits(b, "sf")) return(if (all(keep)) g else g[keep])
  if (!identical(g, sf::st_geometry(b))) sf::st_geometry(b) <- g
  if (!all(keep)) b <- b[keep, ]
  b
}

#' The union of a boundary's rows, each repaired first
#'
#' \code{.boundary_mend()}, then \code{sf::st_union()}, on the sphere for a
#' boundary in longitude and latitude whatever \code{sf::sf_use_s2()} says:
#' with it off, sf printed "although coordinates are longitude/latitude,
#' st_union assumes that they are planar" on every call.
#'
#' @inheritParams .boundary_mend
#' @return An sfc of one geometry.
#' @keywords internal
#' @noRd
.boundary_union <- function(b, refuse = NULL, caller = NULL, task = NULL, all = FALSE,
                            nothing = "nothing can be built or drawn in it") {
  g <- sf::st_geometry(.boundary_mend(b, refuse, caller, task, all, nothing))
  if (.is_longlat(g)) .with_s2(sf::st_union(g)) else sf::st_union(g)
}

#' The error for a boundary row s2 refuses even repaired, or that lies
#' beyond a pole
#'
#' Of class \code{"spatialkit_cells_degenerate"}, built by
#' \code{.raise_cells_degenerate()}: it names `boundary`, its rows and s2's
#' reason for each, and says to transform `boundary` to a projected CRS,
#' where GEOS repairs it.  It carries \code{site = "boundary"}.
#'
#' With \code{cause = "latitude"} the rows have a vertex at a latitude
#' beyond 90 degrees in a geographic CRS (\code{.boundary_located()}).  Such
#' a vertex is no place on the Earth: s2 read the row as a shape over the
#' pole, on the far meridian, and no transform places it.  The error then
#' names the rows and their latitudes and says to check the CRS of
#' `boundary`; it carries \code{cause = "latitude"} as well, so a caller
#' that raises the s2 refusal again in its own words passes this one on as
#' it is.
#'
#' @param parent,rows,reasons,n,s2_forced As \code{.boundary_mend()} gives
#'   them; for \code{cause = "latitude"}, \code{reasons} gives each row's
#'   latitude.
#' @param caller Name of the exported function.
#' @param task What cannot run in degrees; for \code{cause = "latitude"},
#'   what cannot be done in `boundary` ("no seed can be drawn in it").
#' @param cause \code{"s2"} or \code{"latitude"}.
#' @return Does not return.
#' @keywords internal
#' @noRd
.boundary_refused <- function(parent, rows, reasons, n, caller, task, s2_forced = FALSE,
                              cause = c("s2", "latitude")) {
  cause <- match.arg(cause)
  if (cause == "latitude") {
    rows_text <- if (length(rows) == 1L) sprintf("row %d", rows)
                 else if (length(rows) <= 6L)
                   sprintf("rows %s and %d", paste(rows[-length(rows)], collapse = ", "),
                           rows[length(rows)])
                 else sprintf("rows %s and %d more", paste(rows[1:6], collapse = ", "),
                              length(rows) - 6L)
    shown <- seq_len(min(length(rows), 6L))
    msg <- sprintf(paste0(
      "%s(): %s of the %d %s of `boundary` %s a vertex at a latitude beyond 90 degrees ",
      "(%s), which is no place on the Earth, so %s. Coordinates in metres stamped with ",
      "a geographic CRS, or x and y swapped, are the usual cause: check the CRS of ",
      "`boundary` (sf::st_crs()) and the order of its coordinates, and call %s() again."),
      caller, rows_text, n, if (n == 1L) "polygon" else "polygons",
      if (length(rows) == 1L) "has" else "have",
      paste(sprintf("row %d: %s", rows[shown], reasons[shown]), collapse = "; "),
      task, caller)
    stop(errorCondition(msg, class = "spatialkit_cells_degenerate", call = NULL,
                        rows = rows, n = n, reasons = reasons, parent = parent,
                        repaired = FALSE, s2_forced = FALSE, site = "boundary",
                        cause = "latitude"))
  }
  .raise_cells_degenerate(
    parent, rows, reasons, n, caller, "`boundary`", task, again = NULL,
    of = sprintf("the %d %s of `boundary`", n, if (n == 1L) "polygon" else "polygons"),
    s2_forced = s2_forced, source = "", site = "boundary",
    remedy = sprintf(paste0(
      "Transform `boundary` to a projected CRS, such as the one ensure_projected() picks ",
      "for it (sf::st_transform()), and call %s() again; it is then repaired there with ",
      "sf::st_make_valid(), which splits crossing edges."), caller))
}


# Why `boundary` has no polygon with a location, for the error
# .stop_boundary_unlocated() raises: read as given, after a move, or after
# the repair of the rows the engine of the union refuses.
.BOUNDARY_UNLOCATED_WHY <- c(
  given = paste0("every row is empty, or has a coordinate that is not finite or lies ",
                 "beyond 1e150"),
  moved = "every row is empty after its move to the CRS the call works in, which cannot place them",
  s2    = "the repair of the rows s2 refuses leaves no polygon with an area",
  geos  = "the repair of the rows GEOS refuses leaves no polygon with an area")

#' The error for a boundary with no polygon that has a location
#'
#' @param caller Name of the exported function.
#' @param nothing What cannot be done in `boundary` ("no seed can be drawn
#'   in it").
#' @param why One of the names of \code{.BOUNDARY_UNLOCATED_WHY}.
#' @return Does not return.
#' @keywords internal
#' @noRd
.stop_boundary_unlocated <- function(caller, nothing, why) {
  stop(sprintf("%s(): `boundary` has no polygon with a location (%s), so %s.",
               caller, .BOUNDARY_UNLOCATED_WHY[[why]], nothing), call. = FALSE)
}

#' The rows of a boundary that have a location
#'
#' A row of `boundary` that is empty, or that has a vertex with a coordinate
#' that is not finite or lies beyond 1e150 (\code{.unlocated_features()}),
#' covers no ground.  Such rows are left out, read on the coordinates as
#' given (\code{moved = FALSE}), before any move and any repair, so the
#' repair never turns an empty polygon into an empty line and a builder
#' reads a square beside such a row as it reads the square alone.  A move
#' empties a row it cannot place (\code{.transform_keep_unlocated()}), so
#' the test runs again after it (\code{moved = TRUE}).  With no row left the
#' call stops (\code{.stop_boundary_unlocated()}).
#'
#' Read as given, a boundary in a geographic CRS with a vertex at a latitude
#' beyond 90 degrees, on a row that has a location, is refused first
#' (\code{.boundary_refused(cause = "latitude")}): s2 read such a row as a
#' shape folded over the pole, so seeds and cells were placed on the far
#' meridian.
#'
#' Read as given, a row the repair would leave without a polygon (a ring
#' with no height; \code{.repair_leaves_no_polygon()}) is left out too, and
#' with no polygon left the call stops with the reason that the repair
#' leaves none.  Left in, it set the CRS a grid was laid in, entered the
#' cover of the seeds, and came out of each builder's repair as a line.
#' A boundary without a CRS that a call stamps with a geographic CRS is
#' read again as given after the stamp, so its latitude is tested.
#'
#' @param b The boundary, sf or sfc, or \code{NULL}.
#' @param caller Name of the exported function.
#' @param nothing What cannot be done in `boundary`, for the errors.
#' @param moved \code{TRUE} after a move of `boundary`.
#' @return \code{b} itself when every row has a location; otherwise its rows
#'   that have one.
#' @keywords internal
#' @noRd
.boundary_located <- function(b, caller, nothing, moved = FALSE) {
  if (is.null(b)) return(b)
  g <- sf::st_geometry(b)
  if (!length(g)) return(b)
  unloc <- .unlocated_features(g)
  if (!isTRUE(moved) && isTRUE(sf::st_crs(g)$IsGeographic)) {
    far <- which(.beyond_pole_features(g) & !unloc)
    if (length(far)) {
      lat <- vapply(far, function(i) {
        y <- sf::st_coordinates(g[i])[, "Y"]
        y[which.max(ifelse(is.finite(y), abs(y), -Inf))]
      }, numeric(1))
      # Six digits, or as many more as show the latitude beyond 90: rounded
      # to six, 90 + 1e-9 read "latitude 90".
      shown <- vapply(lat, function(v) {
        for (d in 6:17) {
          s <- format(v, digits = d)
          if (abs(as.numeric(s)) > 90) break
        }
        s
      }, character(1))
      .boundary_refused(simpleError("a vertex at a latitude beyond 90 degrees"), far,
                        sprintf("latitude %s", shown), length(g),
                        caller, nothing, cause = "latitude")
    }
  }
  keep <- !unloc & !sf::st_is_empty(g)
  if (!any(keep))
    .stop_boundary_unlocated(caller, nothing, if (isTRUE(moved)) "moved" else "given")
  # Read as given, a row whose repair leaves no polygon (a ring with no
  # height) covers no ground either, and is left out before any move, so the
  # CRS a call picks, the cover of the seeds and the cells are those of the
  # other rows.  With no polygon left the call stops.
  if (!isTRUE(moved)) {
    flat <- .repair_leaves_no_polygon(g, which(keep))
    if (any(flat)) {
      keep <- keep & !flat
      if (!any(keep))
        .stop_boundary_unlocated(caller, nothing, if (.is_longlat(g)) "s2" else "geos")
    }
  }
  if (all(keep)) return(b)
  if (inherits(b, "sf")) b[keep, ] else g[keep]
}

#' Which rows of a boundary its repair leaves without a polygon
#'
#' A row that holds a polygon and that the engine of the boundary's union
#' refuses (s2 in longitude and latitude, as \code{.boundary_mend()} judges
#' it, GEOS otherwise) is repaired, and a row the repair leaves without a
#' polygon (a ring with no height: an empty collection, or a line) is
#' marked.  Rows the engine takes are not repaired, so a boundary of valid
#' rows costs one validity test.
#'
#' @param g An sfc.
#' @param idx The rows to test (those with a location).
#' @return A logical vector, one element per row of \code{g}.
#' @keywords internal
#' @noRd
.repair_leaves_no_polygon <- function(g, idx = seq_along(g)) {
  out <- logical(length(g))
  idx <- idx[.has_polygon(g[idx])]
  if (!length(idx)) return(out)
  if (.is_longlat(g)) {
    .with_s2({
      bad <- idx[.s2_refused(g[idx])]
      if (length(bad)) out[bad] <- !.has_polygon(.safe_make_valid(g[bad]))
    })
  } else {
    bad <- idx[!(suppressWarnings(sf::st_is_valid(g[idx])) %in% TRUE)]
    if (length(bad)) out[bad] <- !.has_polygon(.safe_make_valid(g[bad]))
  }
  out
}

#' Repair the rows of a boundary and leave out those the repair empties
#'
#' The builders repair `boundary` with \code{.safe_make_valid()}.  A row
#' with no height (every vertex on one line, or in degrees on one parallel)
#' comes out of the repair as a line or an empty collection: beside a square
#' the hex and square grids refused the boundary as "not polygonal", the
#' clip target held the line, and alone the Voronoi cells held none of the
#' seeds.  A row that held a polygon before the repair and holds none after
#' it is left out, as \code{.boundary_mend()} leaves it out for the seeders;
#' with no polygon left the call stops (\code{.stop_boundary_unlocated()}).
#' A row that was not a polygon as given is kept, as before, for the type
#' tests of the callers.
#'
#' @param b The boundary, sf or sfc.
#' @param caller Name of the exported function.
#' @param nothing What cannot be done in `boundary`, for the error.
#' @return \code{.safe_make_valid(b)}, without the rows the repair empties.
#' @keywords internal
#' @noRd
.boundary_repaired <- function(b, caller, nothing) {
  before <- .has_polygon(sf::st_geometry(b))
  b <- .safe_make_valid(b)
  g <- sf::st_geometry(b)
  after <- .has_polygon(g)
  lost <- before & !after
  if (!any(lost)) return(b)
  if (!any(after))
    .stop_boundary_unlocated(caller, nothing,
                             if (.is_longlat(g) && isTRUE(sf::sf_use_s2())) "s2" else "geos")
  if (inherits(b, "sf")) b[!lost, ] else g[!lost]
}

#' Whether each row holds a polygon that is not empty
#'
#' The repair of a row s2 or GEOS refuses can leave an empty collection or a
#' line (a ring with no height).
#'
#' @param g An sfc.
#' @return A logical vector, one element per row.
#' @keywords internal
#' @noRd
.has_polygon <- function(g) {
  poly <- function(s) {
    if (inherits(s, "GEOMETRYCOLLECTION"))
      return(any(vapply(s, poly, logical(1))))
    inherits(s, c("POLYGON", "MULTIPOLYGON")) && length(s) > 0L
  }
  vapply(unclass(g), poly, logical(1))
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
#' the seeds of exactly the cells that were scored.  On a profile without an
#' estimate fitted to a subsample, \code{select_resolution()} refits the
#' chosen level to every point first (\code{attr(seeds, "refit")}), and the
#' seeds are then those centres, whose cells are not the cells scored.
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
#'   \code{attr(, "cover")} when a cover is known, \code{attr(, "refit")}
#'   when the centres were refitted to every point, and
#'   \code{attr(, "count_rule")} when \code{count_rule} is a record, marked
#'   \code{scored_cells = TRUE} when these are the centres the profile scored
#'   and \code{FALSE} when they were refitted.
#' @keywords internal
#' @noRd
.use_scored_seeds <- function(seeds, like, cloud, caller, cover = NULL, count_rule = NULL) {
  nstart <- attr(seeds, "nstart", exact = TRUE)
  # The kind that built the scored partition, when the profile recorded one.
  rec <- attr(seeds, "partition", exact = TRUE)
  # The centres refitted to every point, when the profile was fitted to a
  # subsample (select_resolution()).
  refit <- attr(seeds, "refit", exact = TRUE)
  if (is.null(cover)) cover <- attr(seeds, "cover", exact = TRUE)
  moved  <- .seeds_in_crs_of(seeds, like)
  out <- sf::st_sf(seed_id = seq_len(nrow(moved)), method = "kmeans",
                   geometry = sf::st_geometry(moved))
  if (!is.null(refit))
    .log_info(paste0("%s(): returning the %d centres of the partition resolution_profile() ",
                     "scored at that count on a subsample of %d points, refitted to all %d ",
                     "points by select_resolution(), rather than running a new fit; their ",
                     "Voronoi cells are near the cells the criteria judged but not the same. ",
                     "Pass the count as a number for a fresh fit."),
              caller, nrow(out), as.integer(refit$subsample), as.integer(refit$rows))
  else if (is.null(rec))
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
  if (!is.null(refit)) attr(out, "refit") <- refit
  # These are the centres the profile scored at the count, so their Voronoi
  # cells are the cells its numbers describe, and the record says so; or the
  # centres refitted to every point, which are other cells built from the
  # count.
  out <- .count_rule_put(out, count_rule, scored_cells = is.null(refit))
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
    # A projected seed with a coordinate beyond .LONLAT_MOVE_MAX comes out
    # empty (.lonlat_of_placed()): at x = 1e20 in EPSG:3857 the move never
    # came back from PROJ.
    if (isTRUE(.looks_like_lonlat(like)$lonlat))
      seeds <- .lonlat_of_placed(seeds)
    return(sf::st_set_crs(seeds, NA))
  }
  # The same for a `like` in degrees with its own CRS.
  if (from == to) seeds else .lonlat_of_placed(seeds, to)
}


#' Which scored seed each cloud point falls to
#'
#' The \code{attr(, "kmeans")} record of \code{get_voronoi_seeds()} for seeds
#' taken from a resolution profile: each usable cloud point is assigned to its
#' nearest seed in the profile's projected CRS (the seed whose Voronoi cell
#' holds it), by the profile's rule (\code{.nearest_centre()}, a tie to the
#' lower \code{seed_id}), with the sizes and within-cell sums of squares that
#' assignment gives.  \code{iter} is \code{NA}: no k-means ran here.
#'
#' @keywords internal
#' @noRd
.scored_seeds_record <- function(cloud, seeds, nstart) {
  # A cloud point without a location, read as given, is emptied before the
  # move, which can put a northing of 1e200 at a finite place, and is left
  # out below with the empty points.
  cloud <- .empty_unlocated(cloud)
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
  ok <- !is.na(.values_at_xy(cloud)[, 1L])
  sxy <- sf::st_coordinates(seeds)[, 1:2, drop = FALSE]
  k   <- nrow(sxy)
  cl  <- integer(0); d2 <- numeric(0)
  # The seed of each point by the profile's rule (.nearest_centre()), with
  # the seeds numbered by seed_id: squared distances within .TIE_TOL of the
  # least are a tie, and a tie goes to the lower-numbered seed.
  # sf::st_nearest_feature() settled a point the same distance from two
  # seeds by the last bits of its distances, which move with the unit and
  # the origin, and not always as the profile scored it.  A seed without a
  # location is never the nearest.
  live <- which(is.finite(sxy[, 1L]) & is.finite(sxy[, 2L]))
  if (!length(live)) ok[] <- FALSE
  if (any(ok)) {
    cl <- live[.nearest_centre(xy[ok, 1L], xy[ok, 2L], sxy[live, 1L], sxy[live, 2L],
                               use_fnn = .has_backend("FNN"))]
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
