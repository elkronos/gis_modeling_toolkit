# =============================================================================
# Block-kriging adequacy: what a variogram-based aggregator would deliver per
# cell, and whether its variance can be believed -- computed beside the plain
# cell means, changing none of them.
# =============================================================================

#' Block-kriging adequacy diagnostics for a set of cells
#'
#' \code{\link{summarize_by_cell}()} aggregates by plain means inside cell
#' boundaries.  A block-kriging aggregator would instead weight observations
#' by their spatial correlation with the cell, and would return an estimate
#' and a variance for every cell, thin or empty.  Whether that is worth
#' having on a given layer is a question with a measurable answer, and this
#' function measures it, changing no cell value: for every cell it reports
#' the block-kriging estimate and variance implied by a fitted variogram,
#' that variance as a share of the variance the cell's mean would have with
#' no data at all, and, where the cell has points,
#' whether it exceeds the design-based variance of the plain mean,
#' \eqn{s^2/n}; and it scores the variogram itself by blocked
#' cross-validation.
#'
#' @section Reading the columns:
#' \describe{
#'   \item{\code{kr_ratio}}{The block-kriging variance over the cell's prior
#'     variance, in \eqn{[0, 1]}.  The prior variance is the variance the
#'     cell's mean would have with no data at all, \eqn{\bar C(B,B)}: the
#'     covariance averaged over pairs of points in the cell, on the
#'     discretisation \pkg{gstat} block-kriges with, and without the nugget,
#'     which averages out over a block (\pkg{gstat} leaves it out of the
#'     block variance too).  Each cell has its own: a cell's mean varies less
#'     than a single point does, and far less once the cell is wider than the
#'     range, so the point sill is not the scale.  It is the coverage score,
#'     and it needs no hand-set threshold in metres or point counts: as it
#'     approaches 1 the estimate carries almost no information from the data
#'     about the cell and is reverting to the estimated mean.  Ordinary
#'     kriging adds the variance of that estimated mean, so a cell the data do
#'     not reach comes out at or above its prior variance and reads 1.  A
#'     cell at 0.05 is well determined; a cell at 0.8 is mostly prior.
#'     \code{NA} when the model is a pure nugget, where a cell mean has no
#'     prior variance to be a share of.}
#'   \item{\code{kr_exceeds_design}}{\code{TRUE} where the kriging variance is
#'     larger than \eqn{s^2/n} from the cell's own points: kriging is not
#'     earning its keep there, and that is said per cell instead of
#'     globally.  \code{NA} for cells with fewer than two points, where
#'     \eqn{s^2} does not exist.}
#'   \item{\code{kr_shift}}{The kriged estimate minus the plain mean, in
#'     units of the plain mean's standard error (\code{NA} where that is not
#'     defined).  How much the aggregator would move the value, against the
#'     precision the value has.}
#' }
#'
#' @section The cross-validation statistic:
#' The kriging variance is only as good as the variogram.  With blocked
#' folds, each held-out point is kriged from the training folds and the
#' standardised error \eqn{(z - \hat z) / \sqrt{kv}} is recorded; the
#' variance of those over all held-out points should be about 1.  Its
#' departure from 1 measures how badly the kriging variance is understated
#' (or overstated): 1.4 means the variances are roughly 40 percent too
#' small.  Measured on simulated exponential fields (n = 300 on a 1000-unit
#' extent, sill 1, nugget 0.2, five blocked folds, eight draws per
#' configuration): 0.93--1.07 with the true variogram, 0.85--1.01 with the
#' variogram estimated from the same points.  With the nugget understated
#' tenfold it moved only to 0.95--1.24.  That is a property of the folds and
#' not a weakness of the statistic: under blocked folds every held-out point
#' is far from the training data, where the kriging variance is close to the
#' sill whatever the nugget, so the blocked statistic checks the sill and
#' range.  To check the nugget, pass random folds
#' (\code{make_folds(method = "random_kfold")}) as \code{folds}: the
#' held-out points are then close to their neighbours, where the nugget
#' decides the variance.  The statistic is computed fold by fold with
#' \code{gstat::krige()} on the splits \code{\link{make_folds}()} built, each
#' held-out point kriged from that split's own training set, so the folds
#' carry the same separation the package uses everywhere else:
#' \code{"buffered_loo"} and \code{"nndm"} keep the points they exclude
#' around each held-out one out of its kriging, and \code{print()} names the
#' scheme that ran.  A vector of fold labels is run as k-fold, each fold
#' kriged from all the others.
#'
#' @section What it said about block kriging as an aggregator:
#' On the same simulated fields, with 16, 36 and 64 square cells: under
#' uniform sampling the kriged and plain cell means differed by more than
#' one standard error in 11--24 percent of cells and by more than two in
#' 0--7 percent, the kriging variance exceeded \eqn{s^2/n} in 10--26 percent
#' of populated cells, and at most one cell was empty.  Under clustered
#' sampling (eight clusters of 60-unit spread) the two aggregators parted:
#' shifts above one standard error in 34--63 percent of cells and above two
#' in 9--25 percent, the kriging variance below \eqn{s^2/n} in only
#' 13--52 percent of populated cells, and 3--27 of the cells empty, each with
#' a kriged estimate and variance where the plain mean has nothing.
#' So a block-kriging aggregator earns its place on clustered layers and
#' rarely on uniform ones, and this function says which kind a layer is.
#'
#' @section What this needs:
#' A variogram model.  \code{sac} is a \code{\link{estimate_sac_range}()}
#' result carrying one (\code{attr(, "variogram_model")}); when \code{NULL}
#' it is estimated here from the response.  The model families are the ones
#' the package interprets elsewhere: exponential, spherical and Gaussian
#' components with a nugget.  Anything else is refused by name.  A
#' model whose range was not identified (an \code{NA} estimate with the
#' model attached) is used with a warning that says why it was refused, and
#' the reason is kept as \code{attr(, "rejected_reason")}: a range past the
#' fitted lags means the sill was never reached and the ratios rest on an
#' extrapolation; a fit that did not converge stopped wherever the optimiser
#' halted; a variogram that falls with distance, or a range below the
#' shortest lag, describes the data poorly at some lags.
#'
#' The model has to be of the response itself.  A variogram of residuals
#' (\code{estimate_sac_range(predictor_vars = ...)}, \code{attr(,
#' "detrended")} \code{TRUE}) leaves out the variance the predictors
#' explain, while the response is kriged here without them, so
#' \code{kr_var}, \code{kr_ratio} and the cross-validation statistic come out
#' too small (the statistic at 4.3--5.2 against 0.67--1.53 with a spatially
#' structured covariate); it is used with a warning.  The points are put in
#' the CRS the variogram was fitted in (\code{attr(sac, "crs")}), because
#' its range is a length in that CRS's units.  Requires \pkg{gstat}.
#'
#' @section The kriging neighbourhood:
#' \pkg{gstat} kriges a cell from the \code{nmax} locations nearest its
#' centre.  A cell holding more locations than that would be estimated from
#' its middle alone, which describes the middle rather than the cell (in one
#' simulated case a cell of 1,500 points came out 0.43 off at
#' \code{nmax = 50}, against 0.07 from all of them and 0.02 for its plain
#' mean).  So wherever the \code{nmax} locations nearest a cell's centre
#' leave out any of the cell's own locations, the cell is kriged from all of
#' its own locations plus the \code{nmax} nearest outside it;
#' \code{kr_n_used} says how many locations each cell was kriged from.  The
#' cost of a kriging system grows with the cube of its size (about 1 s at
#' 2,000 locations and 17 s at 5,000), so a cell that would need more than
#' \code{max_neighbours} is left out with a warning, its \code{kr_} columns
#' \code{NA}.
#'
#' \pkg{gstat} also discretises each cell into 500 points on a regular grid
#' laid over the cell's whole bounding box, keeping those inside, so the
#' memory a cell takes grows with the ratio of that box to its area: about
#' 54 MB more for a thin diagonal strip at a ratio of 708, and 592 MB at
#' 7,072.  A cell whose bounding box exceeds its area more than
#' \code{max_box_ratio} times (a sliver, parts far apart, a cell that is
#' mostly hole) is left out the same way.  \code{attr(, "cells_left_out")}
#' counts both kinds, and \code{print()} says how many cells have no
#' estimate and why.
#'
#' @section Repeat measurements at one location:
#' Two observations at the same coordinates (visits to a station, records
#' geocoded to one address) make a kriging system singular, because
#' \pkg{gstat} gives them the full sill, nugget included, as their
#' covariance, as if they were one observation.  The kriging and its
#' cross-validation therefore use one observation per location, the mean of
#' its replicates, with a warning; \code{n} and \code{mean} still count every
#' point.  How much of the nugget \eqn{c_0} a mean of \eqn{m} replicates
#' keeps is read off the replicates: their pooled within-location variance
#' \eqn{s_w^2}, capped at \eqn{c_0}, is the part that differs from visit to
#' visit and averages down, and the rest is micro-scale variation the visits
#' share, so the mean carries error variance \eqn{c_0 - s_w^2 + s_w^2/m}.
#' That goes to \pkg{gstat} as a known measurement error (its
#' \code{weights}) on the model with the nugget set to zero.  For replicates
#' that differ only by measurement error this is exactly the kriging of every
#' observation, and for identical replicates it is the kriging of one; a
#' location seen once is kriged as before.  In the cross-validation a
#' location is held out whole, under the fold of its first row, and its
#' error variance is part of its standardised error.  A cell or held-out
#' location \pkg{gstat} still cannot krige is reported \code{NA} with a
#' warning, and \code{print()} says how many.
#'
#' @param assigned_points_sf Points with a cell identifier column, as
#'   \code{\link{assign_features_to_polygons}()} returns.
#' @param response_var The response column.
#' @param cells_sf The cell polygons, with the matching ID column.
#' @param id_col Preferred name of the ID column.  Default \code{"poly_id"}.
#' @param sac Optional \code{sac_range} carrying a variogram model.
#' @param folds Optional \code{\link{make_folds}()} result on
#'   \code{assigned_points_sf} for the cross-validation statistic; built here
#'   with \code{block_kfold} when \code{NULL}.
#' @param k,seed Folds and seed for that construction.
#' @param nmax The number of neighbours each kriging system uses
#'   (\code{gstat}'s \code{nmax}): the locations nearest the cell's centre,
#'   or, for a cell those leave some of its own locations out of, all of its
#'   own plus this many outside it (see "The kriging neighbourhood").  The
#'   cross-validation kriges each held-out point from its \code{nmax} nearest
#'   training locations.  Default 50.
#' @param max_neighbours The largest kriging system a cell is given when its
#'   neighbourhood has to grow to hold all its own locations; a cell that
#'   would need more is left out (\code{kr_} columns \code{NA}) with a
#'   warning.  Never below \code{nmax}.  Default 2000.
#' @param max_box_ratio A cell whose bounding box is more than this many
#'   times its area is left out (\code{kr_} columns \code{NA}) with a
#'   warning, because \pkg{gstat}'s discretisation of it costs memory in
#'   proportion.  Default 1000.
#' @param quiet Suppress progress messages.  Default \code{TRUE}.
#' @return An \code{sf} object of class \code{"kriging_adequacy"}, one row
#'   per cell with the cell geometry and: the ID column, \code{n} (points in
#'   the cell), \code{mean} (the plain mean), \code{se} (its naive standard
#'   error), \code{kr_pred}, \code{kr_var}, \code{kr_ratio},
#'   \code{kr_exceeds_design}, \code{kr_shift} and \code{kr_n_used} (how many
#'   locations the cell was kriged from; \code{NA} for a cell left out).
#'   Attributes:
#'   \code{variogram} (the model frame), \code{sill}, \code{nugget},
#'   \code{range}, \code{range_identified}, \code{rejected_reason} (why the
#'   range was refused, from \code{sac}; \code{NA} when it was identified),
#'   \code{cv} (a list:
#'   \code{zscore_var}, \code{zscore_mean}, \code{rmse}, \code{n_pred},
#'   \code{k}, \code{method}), \code{nmax}, \code{n_points} (the points
#'   used), \code{n_locations} (the distinct locations among them, which
#'   the kriging used) and \code{cells_left_out} (counts of cells left out,
#'   named \code{shape} for \code{max_box_ratio} and \code{size} for
#'   \code{max_neighbours}).
#' @references
#' Cressie, N. (1993). \emph{Statistics for Spatial Data}, revised edition.
#' Wiley. (Block kriging, chapter 3; cross-validation of the kriging
#' variance, section 2.6.4.)
#' @family aggregation
#' @examples
#' if (requireNamespace("gstat", quietly = TRUE)) {
#'   library(sf)
#'   set.seed(1)
#'   n <- 200
#'   x <- 5e5 + runif(n, 0, 1000); y <- 5e6 + runif(n, 0, 1000)
#'   d <- as.matrix(dist(cbind(x, y)))
#'   z <- as.numeric(t(chol(0.8 * exp(-d / 100) + diag(0.2, n))) %*% rnorm(n))
#'   pts <- st_as_sf(data.frame(x = x, y = y, z = z), coords = c("x", "y"), crs = 32632)
#'   bnd <- st_sf(geometry = st_as_sfc(st_bbox(pts)))
#'   cells <- create_grid_polygons(bnd, target_cells = 16, type = "square")
#'   asg <- assign_features_to_polygons(pts, cells)
#'   ka <- kriging_adequacy(asg, "z", cells, k = 4)
#'   print(ka)                   # the report; print() because only a block's
#'                               # last value shows on its own
#'   attr(ka, "cv")$zscore_var   # about 1 when the variogram is right
#' }
#' @export
kriging_adequacy <- function(assigned_points_sf, response_var, cells_sf,
                             id_col = "poly_id", sac = NULL, folds = NULL,
                             k = 5L, seed = 123L, nmax = 50L,
                             max_neighbours = 2000L, max_box_ratio = 1000,
                             quiet = TRUE) {
  .msg <- function(...) if (!quiet) message(...)
  if (!requireNamespace("gstat", quietly = TRUE))
    stop("kriging_adequacy(): package 'gstat' is required.", call. = FALSE)
  .assert_sf(assigned_points_sf, c("POINT", "MULTIPOINT"), "assigned_points_sf",
             caller = "kriging_adequacy")
  .assert_sf(cells_sf, c("POLYGON", "MULTIPOLYGON"), "cells_sf",
             caller = "kriging_adequacy")
  if (!is.character(response_var) || length(response_var) != 1L ||
      !(response_var %in% names(assigned_points_sf)))
    stop("kriging_adequacy(): `response_var` must name a column of `assigned_points_sf`.",
         call. = FALSE)
  if (!is.numeric(nmax) || length(nmax) != 1L || !is.finite(nmax) || nmax < 1)
    stop("kriging_adequacy(): `nmax` must be a single positive number.", call. = FALSE)

  # --- the ID column, as summarize_by_cell() finds it ---
  id_candidates <- unique(c(id_col, "poly_id", "polygon_id", "cell_id"))
  id_pts   <- id_candidates[id_candidates %in% names(assigned_points_sf)]
  id_cells <- id_candidates[id_candidates %in% names(cells_sf)]
  if (!length(id_pts) || !length(id_cells))
    stop("kriging_adequacy(): could not find a cell ID column in both layers. ",
         "Looked for: ", paste(id_candidates, collapse = ", "), ".", call. = FALSE)
  id_pts <- id_pts[[1L]]; id_cells <- id_cells[[1L]]

  if (!is.numeric(max_neighbours) || length(max_neighbours) != 1L ||
      !is.finite(max_neighbours) || max_neighbours < 1)
    stop("kriging_adequacy(): `max_neighbours` must be a single positive number.", call. = FALSE)
  if (!is.numeric(max_box_ratio) || length(max_box_ratio) != 1L ||
      is.na(max_box_ratio) || max_box_ratio < 1)
    stop("kriging_adequacy(): `max_box_ratio` must be a single number of at least 1.",
         call. = FALSE)

  # --- geometry: projected points, cells in the same CRS ---
  pts <- assigned_points_sf
  if (!all(sf::st_geometry_type(pts, by_geometry = TRUE) == "POINT"))
    pts <- coerce_to_points(pts, "auto")
  # A supplied variogram's range is a length in the CRS it was fitted in, so
  # the points go into that CRS, as summarize_by_cell() puts them.  Projected
  # anew here instead, a sac fitted in metres met points in km (or US feet)
  # with nothing said: the variance CV statistic came out at 3.05 against 0.85.
  sac_crs <- if (!is.null(sac)) attr(sac, "crs") else NULL
  pts <- if (inherits(sac_crs, "crs") && !is.na(sac_crs))
    .transform_or_stamp(pts, sac_crs, what = "assigned_points_sf",
                        caller = "kriging_adequacy")
  else ensure_projected(pts)
  # Row identities survive the completeness filter below, so folds built on
  # the layer as passed (or here, on the kept rows) map back by `..row_id`.
  if (!("..row_id" %in% names(pts))) pts$..row_id <- seq_len(nrow(pts))
  cells <- .align_crs(cells_sf, pts)
  cells <- .safe_make_valid(cells)
  z <- suppressWarnings(as.numeric(sf::st_drop_geometry(pts)[[response_var]]))
  ok <- is.finite(z) & !sf::st_is_empty(pts)
  if (sum(ok) < 10L)
    stop("kriging_adequacy(): fewer than 10 points with a finite response; nothing ",
         "to krige from.", call. = FALSE)
  if (any(!ok))
    .log_warn("kriging_adequacy(): dropping %d point(s) with a missing response or empty geometry.",
              sum(!ok))
  kp <- pts[ok, , drop = FALSE]
  kp$..z <- z[ok]

  # --- the variogram model ---
  if (is.null(sac)) {
    .msg("kriging_adequacy(): estimating the variogram of ", response_var, " ...")
    sac <- estimate_sac_range(kp, "..z", seed = seed)
  }
  vm <- attr(sac, "variogram_model")
  if (is.null(vm) || !is.data.frame(vm))
    stop("kriging_adequacy(): `sac` carries no variogram model (estimate_sac_range() ",
         "returned a bare NA), so there is nothing to krige with.", call. = FALSE)
  fam <- as.character(vm$model)
  bad_fam <- setdiff(fam, c("Nug", "Exp", "Sph", "Gau"))
  if (length(bad_fam))
    stop("kriging_adequacy(): the variogram model has a component of family ",
         paste(sQuote(bad_fam), collapse = ", "), ", which this package does not ",
         "interpret (exponential, spherical and Gaussian, plus a nugget).",
         call. = FALSE)
  nugget <- sum(vm$psill[fam == "Nug"], na.rm = TRUE)
  sill   <- sum(vm$psill, na.rm = TRUE)
  if (!is.finite(sill) || sill <= 0)
    stop("kriging_adequacy(): the variogram model has no positive sill.", call. = FALSE)
  range_identified <- is.finite(suppressWarnings(as.numeric(sac)))
  # Why the range was refused, kept on the result: "sill never reached" was
  # said of every refusal, and a non-converged fit, a variogram that falls
  # with distance and one that runs past the fitted lags are different
  # findings about the model the variances rest on.
  rejected_reason <- if (range_identified) NA_character_ else
    as.character(attr(sac, "rejected_reason") %||% "no effective range")[1L]
  if (!range_identified)
    .warn_and_log(paste0("kriging_adequacy(): the variogram's range was not ",
                         "identified (%s), so the model the kriging variances ",
                         "and ratios rest on is one the data did not pin down."),
                  rejected_reason)
  # A variogram of residuals on predictors describes what the predictors left
  # over, and it is the response itself that is kriged here, with no
  # predictors: the variance they explained is missing from the model.
  # Measured with a spatially structured covariate: the variance of the
  # standardised CV errors went from 0.67-1.53 to 4.3-5.2, and no cell was
  # flagged as kriging worse than its plain mean where 2-7 had been.
  if (isTRUE(attr(sac, "detrended")))
    .warn_and_log(paste0("kriging_adequacy(): `sac` is the variogram of the ",
                         "residuals on predictors (detrend = \"%s\"), but the ",
                         "response itself is kriged here, without them; the ",
                         "variance they explain is missing from the model, so ",
                         "kr_var, kr_ratio and the cross-validation statistic ",
                         "understate the uncertainty. Pass a variogram of the ",
                         "response (estimate_sac_range() without predictor_vars) ",
                         "or sac = NULL."),
                  as.character(attr(sac, "detrend_method") %||% "ols")[1L])

  # --- one observation per location for the kriging systems ---
  # gstat gives two observations at distance zero the full sill, nugget
  # included, as their covariance, as if they were one, so repeat visits to a
  # station make every kriging system holding them singular, and gstat
  # returns NA and says so only at debug.level > 0.  Each location is kriged
  # from the mean of its replicates instead (the plain cell means below still
  # use every point).  The replicates say how much of the nugget that mean
  # keeps: their pooled within-location variance is the part that differs
  # between visits (capped at the nugget) and averages down by the count; the
  # rest is micro-scale variation the visits share.  That goes to gstat as a
  # known measurement error (`weights` = 1 / variance) on the model with its
  # nugget zeroed: a location seen once gets the whole nugget, as before;
  # replicates differing only by measurement error give the kriging of every
  # observation; identical replicates give the kriging of one.
  xy  <- sf::st_coordinates(kp)[, 1:2, drop = FALSE]
  key <- paste(xy[, 1], xy[, 2])
  loc <- match(key, unique(key))
  first <- !duplicated(loc)
  kk <- kp[first, , drop = FALSE]
  vk <- vm
  w <- NULL
  err_var <- rep(0, nrow(kk))    # already in vk's nugget when w is NULL
  if (any(!first)) {
    m_loc <- tabulate(loc)
    kk$..z <- as.numeric(rowsum(kp$..z, loc)) / m_loc
    me <- min(sum((kp$..z - kk$..z[loc])^2) / sum(m_loc - 1L), nugget)
    err_var <- nugget - me + me / m_loc
    vk$psill[fam == "Nug"] <- 0
    w <- 1 / err_var
    .warn_and_log(paste0("kriging_adequacy(): %d point(s) share a location with another; ",
                         "kriging from the mean at each of the %d distinct locations, with ",
                         "the part of the nugget that differs between replicates (%.3g of ",
                         "%.3g) divided by their count."),
                  sum(m_loc[loc] > 1L), nrow(kk), me, nugget)
  }

  ids_pts <- as.character(sf::st_drop_geometry(kp)[[id_pts]])
  ids_cells <- as.character(sf::st_drop_geometry(cells)[[id_cells]])

  # --- cells gstat cannot discretise at a bounded cost ---
  # gstat discretises a polygon block with spsample(n = 500, type =
  # "regular"), which lays its grid over the feature's whole bounding box and
  # keeps the points inside, so the memory it takes grows with the ratio of
  # that box to the cell's area: +54 MB for a diagonal strip at 708, +592 MB
  # at 7,072, and gigabytes beyond.  A sliver, a multi-part cell with parts
  # far apart or a cell that is mostly hole is left out instead, and counted.
  box_ratio <- .cell_box_ratio(cells)
  sliver <- !(is.finite(box_ratio) & box_ratio <= max_box_ratio)
  if (any(sliver))
    .warn_and_log(paste0("kriging_adequacy(): %d cell(s) left out (kr_ columns NA): ",
                         "the bounding box of each is more than max_box_ratio = %s ",
                         "times its area (up to %s), and gstat discretises a block ",
                         "over its whole bounding box, at a memory cost in proportion."),
                  sum(sliver), format(max_box_ratio),
                  format(signif(max(box_ratio[sliver]), 3)))

  # --- each cell's kriging neighbourhood ---
  # gstat takes the nmax locations nearest a block's centre, so a cell holding
  # more than that was kriged from its middle alone: one 60 km cell of 1,500
  # points came out 0.43 off at nmax = 50 against 0.07 from all of them, and
  # the plain mean it is compared with was 0.02 off.  Where the nmax nearest
  # the centre leave out any of the cell's own locations, it is kriged from
  # all of them plus the nmax nearest outside it; elsewhere gstat's own
  # neighbourhood is used unchanged.
  cap <- max(as.integer(max_neighbours), as.integer(nmax))
  loc_cell <- match(ids_pts[first], ids_cells)
  nb <- .cell_neighbourhoods(cells, kk, loc_cell, as.integer(nmax), cap, skip = sliver)
  if (any(nb$enlarged))
    .log_info(paste0("kriging_adequacy(): %d cell(s) hold locations the nmax = %d ",
                     "nearest their centre leave out; each is kriged from all of its ",
                     "own locations plus the nmax nearest outside it (largest ",
                     "system: %d)."),
              sum(nb$enlarged), as.integer(nmax), max(nb$n_used[nb$enlarged]))
  if (any(nb$too_big))
    .warn_and_log(paste0("kriging_adequacy(): %d cell(s) left out (kr_ columns NA): ",
                         "kriging one from all of its own locations plus the nmax ",
                         "nearest outside it takes a system of up to %d, above ",
                         "max_neighbours = %d (the cost grows with the cube of the ",
                         "size). Raise max_neighbours to krige them."),
                  sum(nb$too_big), max(nb$size[nb$too_big]), cap)

  # --- block kriging onto the cells ---
  .msg("kriging_adequacy(): block kriging onto ", nrow(cells), " cells ...")
  .krige_cells <- function(loc, nd, wt, nmx)
    tryCatch(
      gstat::krige(..z ~ 1, locations = loc, newdata = nd, model = vk,
                   weights = wt, nmax = nmx, debug.level = 0),
      error = function(e)
        stop("kriging_adequacy(): gstat::krige() failed: ", conditionMessage(e),
             call. = FALSE))
  kr_pred <- kr_var <- rep(NA_real_, nrow(cells))
  batch <- !sliver & !nb$enlarged & !nb$too_big
  if (any(batch)) {
    bk <- .krige_cells(kk, cells[batch, , drop = FALSE], w, as.integer(nmax))
    kr_pred[batch] <- suppressWarnings(as.numeric(bk$var1.pred))
    kr_var[batch]  <- suppressWarnings(as.numeric(bk$var1.var))
  }
  for (ci in which(nb$enlarged)) {
    s <- nb$sel[[ci]]
    bk <- .krige_cells(kk[s, , drop = FALSE], cells[ci, , drop = FALSE], w[s], Inf)
    kr_pred[ci] <- suppressWarnings(as.numeric(bk$var1.pred))
    kr_var[ci]  <- suppressWarnings(as.numeric(bk$var1.var))
  }
  kr_var[is.finite(kr_var) & kr_var < 0] <- 0
  # gstat answers a system it cannot solve with NA and, at debug.level 0,
  # nothing else.
  kriged <- batch | nb$enlarged
  n_na <- sum(kriged & (!is.finite(kr_pred) | !is.finite(kr_var)))
  if (n_na)
    .warn_and_log(paste0("kriging_adequacy(): gstat::krige() returned no estimate for %d ",
                         "of %d cell(s) (a kriging system it could not solve); their ",
                         "kr_ columns are NA."), n_na, sum(kriged))
  # The scale for kr_ratio.  kr_var is the variance of a cell MEAN, and with
  # no data it levels off at the cell's own C(B,B), not at the point sill:
  # divided by the sill, an empty cell far from every datum read about 0.1,
  # and a large empty cell ranked below small populated ones.  Not computed
  # for a sliver, whose discretisation is what it was left out to avoid.
  prior_var <- rep(NA_real_, nrow(cells))
  if (any(!sliver))
    prior_var[!sliver] <- .block_prior_var(cells[!sliver, , drop = FALSE], vm)

  # --- the plain means and their design-based variance, per cell ---
  n_c  <- as.integer(table(factor(ids_pts, levels = ids_cells)))
  m_c  <- tapply(kp$..z, factor(ids_pts, levels = ids_cells), mean)
  v_c  <- tapply(kp$..z, factor(ids_pts, levels = ids_cells), stats::var)
  s2_n <- as.numeric(v_c) / n_c
  s2_n[n_c < 2L] <- NA_real_
  se_c <- sqrt(s2_n)

  out <- cells[, id_cells, drop = FALSE]
  out$n <- n_c
  out$mean <- as.numeric(m_c)
  out$se <- se_c
  out$kr_pred <- kr_pred
  out$kr_var <- kr_var
  out$kr_ratio <- ifelse(prior_var > 0, pmin(pmax(kr_var / prior_var, 0), 1), NA_real_)
  out$kr_exceeds_design <- ifelse(is.finite(s2_n), kr_var > s2_n, NA)
  out$kr_shift <- ifelse(is.finite(se_c) & se_c > 0, (kr_pred - out$mean) / se_c, NA_real_)

  out$kr_n_used <- nb$n_used

  # --- the cross-validation statistic ---
  if (is.null(folds)) {
    folds <- make_folds(kp, k = k, method = "block_kfold", seed = seed)
    method <- "block_kfold"
  } else {
    method <- if (is.list(folds) && !is.null(folds$method)) folds$method else "supplied labels"
  }
  # One split per fold, in locations: the ones held out, and the ones they are
  # kriged from.  A make_folds() result is run on its own training sets, which
  # for buffered_loo and nndm leave out the points around each held-out one;
  # reduced to fold labels, as they were, both ran as plain leave-one-out
  # (a 250 m buffer: variance 1.02 and RMSE 0.757, identical to 1:n, against
  # 0.873 and 0.944 with the buffer kept).
  splits <- .location_splits(folds, kp, first)
  held <- sort(unique(unlist(lapply(splits, `[[`, "test"), use.names = FALSE)))
  cv <- list(zscore_var = NA_real_, zscore_mean = NA_real_, rmse = NA_real_,
             n_pred = 0L, k = length(splits), method = method)
  if (length(held) >= 10L && cv$k >= 2L) {
    .msg("kriging_adequacy(): cross-validating the kriging variance over ", cv$k, " folds ...")
    # The folds are run here rather than by gstat::krige.cv(), which subsets
    # the data per fold but not `weights`.  With the nugget zeroed in vk, a
    # held-out location's own error variance goes back into its variance.
    kcv <- tryCatch({
      pred <- pvar <- rep(NA_real_, nrow(kk))
      for (s in splits) {
        if (!length(s$train)) next
        kf <- gstat::krige(..z ~ 1, locations = kk[s$train, , drop = FALSE],
                           newdata = kk[s$test, , drop = FALSE], model = vk,
                           weights = w[s$train], nmax = as.integer(nmax), debug.level = 0)
        pred[s$test] <- suppressWarnings(as.numeric(kf$var1.pred))
        pvar[s$test] <- suppressWarnings(as.numeric(kf$var1.var)) + err_var[s$test]
      }
      data.frame(residual = kk$..z[held] - pred[held],
                 zscore = suppressWarnings((kk$..z[held] - pred[held]) / sqrt(pvar[held])))
    }, error = function(e) {
      .log_warn(paste0("kriging_adequacy(): gstat::krige() failed in cross-validation (%s); ",
                       "no cross-validation statistic."), conditionMessage(e))
      NULL
    })
    if (!is.null(kcv)) {
      zs <- suppressWarnings(as.numeric(kcv$zscore)); zs <- zs[is.finite(zs)]
      res <- suppressWarnings(as.numeric(kcv$residual)); res <- res[is.finite(res)]
      cv$zscore_var  <- if (length(zs) > 1L) stats::var(zs) else NA_real_
      cv$zscore_mean <- if (length(zs)) mean(zs) else NA_real_
      cv$rmse        <- if (length(res)) sqrt(mean(res^2)) else NA_real_
      cv$n_pred      <- length(zs)
      if (cv$n_pred < length(held))
        .warn_and_log(paste0("kriging_adequacy(): gstat::krige() returned no cross-validation ",
                             "prediction it could use for %d of %d held-out location(s); ",
                             "the statistic uses the other %d."),
                      length(held) - cv$n_pred, length(held), cv$n_pred)
    }
  } else {
    .log_warn("kriging_adequacy(): too few points or folds for the cross-validation statistic.")
  }

  structure(out,
            variogram = vm, sill = sill, nugget = nugget,
            range = suppressWarnings(as.numeric(sac)),
            range_identified = range_identified,
            rejected_reason = rejected_reason,
            cv = cv, nmax = as.integer(nmax), n_points = nrow(kp),
            n_locations = nrow(kk),
            cells_left_out = c(shape = sum(sliver), size = sum(nb$too_big)),
            response_var = response_var,
            class = c("kriging_adequacy", class(out)))
}


#' The whole feature's bounding-box area over its net area, per cell
#'
#' What gstat's polygon discretisation costs in proportion to:
#' spsample(type = "regular") lays its grid over the feature's bounding box
#' and keeps the points inside, as sp computes it (the box of all parts
#' together, holes subtracted from the area).  Inf for a cell with no area.
#' @keywords internal
#' @noRd
.cell_box_ratio <- function(cells) {
  g <- sf::st_geometry(cells)
  a <- suppressWarnings(as.numeric(sf::st_area(g)))
  box <- vapply(g, function(x) {
    b <- sf::st_bbox(x)
    as.numeric((b[["xmax"]] - b[["xmin"]]) * (b[["ymax"]] - b[["ymin"]]))
  }, numeric(1))
  ifelse(is.finite(a) & a > 0, box / a, Inf)
}


#' Which locations each cell is block-kriged from
#'
#' gstat takes the `nmax` locations nearest a block's centre (the polygon's
#' label point, sp::coordinates(), which is what predict.gstat() hands it).
#' A cell whose own locations all fall among those keeps that neighbourhood
#' and is kriged in the one batch call.  One that holds a location the
#' neighbourhood leaves out is marked `enlarged` and gets `sel`: all of its
#' own locations plus the `nmax` nearest outside it; `too_big` when that
#' exceeds `cap`.  `n_used` is the neighbourhood size per cell, NA for a cell
#' left out; `size` the system an enlarged or too-big cell needs.
#' @keywords internal
#' @noRd
.cell_neighbourhoods <- function(cells, kk, loc_cell, nmax, cap, skip) {
  nc <- nrow(cells); nl <- nrow(kk)
  n_used <- ifelse(skip, NA_integer_, min(nmax, nl))
  enlarged <- too_big <- rep(FALSE, nc)
  size <- rep(NA_integer_, nc)
  sel <- vector("list", nc)
  own_of <- split(seq_len(nl), factor(loc_cell, levels = seq_len(nc)))
  todo <- which(!skip & lengths(own_of) > 0L)
  if (!length(todo)) return(list(n_used = n_used, enlarged = enlarged,
                                 too_big = too_big, size = size, sel = sel))
  cen <- sp::coordinates(sf::as_Spatial(sf::st_set_crs(
    sf::st_geometry(cells)[todo], sf::NA_crs_)))
  xy <- sf::st_coordinates(kk)[, 1:2, drop = FALSE]
  for (j in seq_along(todo)) {
    ci <- todo[j]; own <- own_of[[ci]]
    d <- sqrt((xy[, 1] - cen[j, 1])^2 + (xy[, 2] - cen[j, 2])^2)
    # gstat's nmax nearest already hold every location of the cell.
    if (sum(d <= max(d[own])) <= nmax) next
    n_out <- min(nmax, nl - length(own))
    size[ci] <- length(own) + n_out
    if (size[ci] > cap) {
      too_big[ci] <- TRUE; n_used[ci] <- NA_integer_
      next
    }
    d[own] <- Inf
    sel[[ci]] <- c(own, order(d)[seq_len(n_out)])
    enlarged[ci] <- TRUE; n_used[ci] <- size[ci]
  }
  list(n_used = n_used, enlarged = enlarged, too_big = too_big, size = size, sel = sel)
}


#' Cross-validation splits in locations, from any fold shape the package accepts
#'
#' Row indices into `kk` (one row per distinct location, `first` marking the
#' first row of each location in `pts`).  A make_folds() result keeps its own
#' train sets, which for buffered_loo and nndm are not the complement of the
#' test set; fold labels become k-fold splits.  A location is held out under
#' its first row and never trains the split it is held out in.  Splits with
#' nothing held out are dropped.
#' @keywords internal
#' @noRd
.location_splits <- function(folds, pts, first) {
  rid <- pts$..row_id[first]
  is_split <- function(s) is.list(s) && !is.null(s$train) && !is.null(s$test)
  sp <- if (is.list(folds) && is.list(folds$folds) && length(folds$folds) &&
            all(vapply(folds$folds, is_split, logical(1)))) {
    lapply(folds$folds, function(s) {
      te <- which(rid %in% s$test)
      list(train = setdiff(which(rid %in% s$train), te), test = te)
    })
  } else {
    lab <- .fold_labels_for(folds, pts)[first]
    ok <- which(is.finite(lab))
    lapply(unique(lab[ok]), function(f)
      list(train = ok[lab[ok] != f], test = ok[lab[ok] == f]))
  }
  Filter(function(s) length(s$test) > 0L, sp)
}


#' The variance each cell's mean has with no data, \eqn{\bar C(B,B)}
#'
#' The covariance averaged over all pairs of points of the cell, on the
#' discretisation gstat's predict() gives a polygon (its default `sps.args`:
#' spsample(n = 500, type = "regular", offset = c(0.5, 0.5))), and without
#' the nugget, which gstat leaves out of block-to-block covariances.  This is
#' the block-kriging variance of a cell no datum informs, before ordinary
#' kriging adds the variance of the estimated mean.  0 for a pure nugget.
#' @keywords internal
#' @noRd
.block_prior_var <- function(cells, vm) {
  sv <- vm[as.character(vm$model) != "Nug", , drop = FALSE]
  c0 <- sum(sv$psill)
  if (!nrow(sv) || !(c0 > 0)) return(rep(0, nrow(cells)))
  geo <- sf::as_Spatial(sf::st_set_crs(sf::st_geometry(cells), sf::NA_crs_))
  vapply(seq_along(geo), function(i) {
    g <- sp::coordinates(sp::spsample(geo[i], n = 500, type = "regular",
                                      offset = c(0.5, 0.5)))
    m <- nrow(g)
    if (m < 2L) return(c0)
    # Each pair once from dist(), both orders, plus the m zero-distance pairs.
    cv <- gstat::variogramLine(sv, dist_vector = as.numeric(stats::dist(g)),
                               covariance = TRUE)$gamma
    (2 * sum(cv) + m * c0) / m^2
  }, numeric(1))
}


#' Fold labels for a layer, from any of the fold shapes the package accepts
#' @keywords internal
#' @noRd
.fold_labels_for <- function(folds, pts) {
  n <- nrow(pts)
  if (is.list(folds) && !is.null(folds$assignment)) {
    a <- folds$assignment
    if ("..row_id" %in% names(pts))
      return(a$fold[match(pts$..row_id, a$row_id)])
    return(a$fold[match(seq_len(n), a$row_id)])
  }
  if (is.atomic(folds) && length(folds) == n) return(as.integer(as.factor(folds)))
  stop("kriging_adequacy(): `folds` must be a make_folds() result on the same ",
       "layer, or a vector of fold labels with one entry per row.", call. = FALSE)
}


#' @export
print.kriging_adequacy <- function(x, ...) {
  # `[` on a data frame keeps the class and drops the attributes, so a column
  # or row subset of this object arrives here with no variogram, no `cv` and no
  # point count, and the summary below cannot be built from it.  Subsetting a
  # table is an ordinary thing to do, and knitr reaches print() through
  # knit_print.data.frame() without being asked, so the subset prints as what it
  # now is instead of failing on the first missing piece.
  if (is.null(attr(x, "variogram", exact = TRUE)) || is.null(attr(x, "n_points"))) {
    y <- x
    class(y) <- setdiff(class(y), "kriging_adequacy")
    cat("Block-kriging adequacy (subset; the fitted summary is not carried",
        "by a subset)\n")
    print(y)
    return(invisible(x))
  }
  df <- sf::st_drop_geometry(x)
  cv <- attr(x, "cv")
  n_loc <- attr(x, "n_locations") %||% attr(x, "n_points")
  pts_word <- if (n_loc < attr(x, "n_points")) "locations" else "points"
  cat(sprintf("Block-kriging adequacy over %d cells (%d points%s, nmax %d)\n",
              nrow(df), attr(x, "n_points"),
              if (pts_word == "locations") sprintf(" at %d distinct locations", n_loc) else "",
              attr(x, "nmax")))
  vm <- attr(x, "variogram", exact = TRUE)
  cat(sprintf("  variogram: %s; sill %.3g, nugget %.3g (%.0f%%), range %s\n",
              paste(sprintf("%s(%.3g, %.3g)", vm$model, vm$psill, vm$range), collapse = " + "),
              attr(x, "sill"), attr(x, "nugget"),
              100 * attr(x, "nugget") / attr(x, "sill"),
              if (isTRUE(attr(x, "range_identified")))
                sprintf("%.1f", attr(x, "range"))
              else if (is.character(attr(x, "rejected_reason")) &&
                       !is.na(attr(x, "rejected_reason")))
                sprintf("not identified (%s)", attr(x, "rejected_reason"))
              else "not identified"))
  # Every line below counts only the cells with an estimate, so say first
  # how many have none, and why: left out by shape or size, or unsolved.
  kr_ok <- is.finite(df$kr_pred) & is.finite(df$kr_var)
  if (!all(kr_ok)) {
    lo <- attr(x, "cells_left_out")
    n_shape <- if (length(lo) && !is.na(lo["shape"])) lo[["shape"]] else 0L
    n_size  <- if (length(lo) && !is.na(lo["size"])) lo[["size"]] else 0L
    n_unsolved <- sum(!kr_ok) - n_shape - n_size
    why <- c(if (n_shape) sprintf("%d left out as too thin or scattered for their bounding box", n_shape),
             if (n_size) sprintf("%d left out as holding more locations than max_neighbours", n_size),
             if (n_unsolved > 0L) {
               if (n_shape + n_size) sprintf("%d whose kriging system gstat could not solve",
                                             n_unsolved)
               else "gstat could not solve their kriging systems"
             })
    cat(sprintf("  no kriged estimate for %d of the %d cells (%s)\n",
                sum(!kr_ok), nrow(df), paste(why, collapse = "; ")))
  }
  r <- df$kr_ratio[is.finite(df$kr_ratio)]
  if (length(r))
    cat(sprintf(paste0("  kriging variance / no-data variance of the cell mean (1 = nothing ",
                       "from the data): median %.3f, range %.3f-%.3f; %d cell(s) above 0.5\n"),
                stats::median(r), min(r), max(r), sum(r > 0.5)))
  pop <- df[is.finite(df$kr_exceeds_design), , drop = FALSE]
  if (nrow(pop))
    cat(sprintf("  kriging variance exceeds s^2/n in %d of the %d cell(s) with two or more points\n",
                sum(pop$kr_exceeds_design), nrow(pop)))
  sh <- df$kr_shift[is.finite(df$kr_shift)]
  if (length(sh))
    cat(sprintf("  kriged minus plain mean: |shift| > 1 SE in %d of %d cell(s), > 2 SE in %d\n",
                sum(abs(sh) > 1), length(sh), sum(abs(sh) > 2)))
  n_empty <- sum(df$n == 0L)
  n_empty_kr <- sum(df$n == 0L & kr_ok)
  cat(sprintf("  empty cells: %d (kriged estimate and variance available for %s)\n",
              n_empty, if (n_empty_kr == n_empty) "each" else sprintf("%d of them", n_empty_kr)))
  # Named for the scheme that ran: every scheme used to print as "blocked
  # CV", random folds and leave-one-out included.
  cv_label <- switch(cv$method %||% "block_kfold",
                     block_kfold        = "blocked CV",
                     random_kfold       = "random k-fold CV",
                     buffered_loo       = "buffered leave-one-out CV",
                     nndm               = "NNDM leave-one-out CV",
                     leave_location_out = "leave-location-out CV",
                     "CV")
  if (is.finite(cv$zscore_var %||% NA_real_))
    cat(sprintf(paste0("  %s (%s, %d folds, %d %s): var of standardised ",
                       "error %.2f (1 = kriging variance correct; above 1 = ",
                       "understated), mean %.2f, RMSE %.3g\n"),
                cv_label, cv$method, cv$k, cv$n_pred, pts_word, cv$zscore_var,
                cv$zscore_mean, cv$rmse))
  else cat(sprintf("  %s: not computed\n", cv_label))
  invisible(x)
}
