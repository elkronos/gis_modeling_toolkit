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
#' no data at all, and, where the cell has points, the error variance of its
#' plain mean under the same variogram and the share of that variance kriging
#' removes; and it scores the variogram itself by blocked cross-validation.
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
#'   \item{\code{plain_var}}{The error variance of the cell's plain mean as
#'     an estimate of the cell's average, under the same variogram and the
#'     same discretisation of the cell as \code{kr_var}:
#'     \eqn{a' C a - 2 a' \bar c(B) + \bar C(B,B)}, with \eqn{a} the plain
#'     mean's weights on the cell's observations, \eqn{C} their covariances
#'     (the nugget on the diagonal), \eqn{\bar c(B)} their covariances with
#'     the cell and \eqn{\bar C(B,B)} the covariance within it.  It is the
#'     model's counterpart of \code{se^2}.  \eqn{s^2/n} is the variance the
#'     plain mean would have if the cell's points were independent draws;
#'     under spatial correlation it understates the error, most where a
#'     cell's points sit close together, so \code{se^2 / plain_var} below 1
#'     says the naive standard error is too small if the variogram is right.
#'     Defined from one point up.  \code{NA} for a cell with no point and for
#'     a cell left out of the kriging.}
#'   \item{\code{kr_gain}}{\code{1 - kr_var / plain_var}: the share of the
#'     plain mean's error variance that kriging removes in the cell.  Both
#'     variances come from the one variogram, so the gain does not move when
#'     sill and nugget are multiplied by a constant (unless locations are
#'     repeated, whose replicate variance is read off the data and does not
#'     scale with the model).  It is at least 0: the
#'     kriging neighbourhood always holds the cell's own points, the plain
#'     mean is one unbiased combination of them, and kriging is the best
#'     such combination.  It is what the variogram expects, not what
#'     happened: read it beside \code{attr(, "cv")$zscore_var}, which says
#'     whether the variogram's variances are the right size.  \code{NA}
#'     where \code{plain_var} is.  \code{kr_var} is not compared with
#'     \eqn{s^2/n}.  That comparison
#'     sets a variance from the model against one from the cell's own
#'     spread, so it moves with the fitted sill: the same data with the same
#'     kriged estimates put 4, 27 and 67 of 70 cells above as sill and
#'     nugget were multiplied by 1/3, 1 and 3.  And it does not mark cells
#'     where kriging fails to pay: against known cell averages on simulated
#'     fields (eight runs, clustered and uniform layouts, estimated and true
#'     variograms) kriging cut the plain mean's error by 20--53 percent in
#'     the cells above and by 27--60 percent in the rest.}
#'   \item{\code{kr_shift}}{The kriged estimate minus the plain mean, over
#'     the standard deviation the variogram gives that difference,
#'     \code{sqrt(plain_var - kr_var)} (the kriging error is uncorrelated
#'     with the difference, so the two variances subtract).  Under the right
#'     variogram it has variance 1: about 32 percent of cells beyond 1 and 5
#'     percent beyond 2.  More than that says the plain means and the
#'     variogram disagree, so the counts \code{print()} gives are a check of
#'     the variogram on the cells themselves, beside the cross-validation
#'     statistic.  In units of the naive standard error \eqn{s/\sqrt{n}}
#'     the same shift has a variance of 1.4--8.4 on clustered layouts, so
#'     counts beyond 1 and 2 in that unit overstate the disagreement.
#'     \code{NA} for a cell with no point, and where the two variances agree
#'     to within a ten-thousandth, which is where kriging returns the plain
#'     mean.}
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
#' On the same simulated fields, with 16, 36 and 64 square cells, under
#' uniform sampling and under clustered sampling (eight clusters of 60-unit
#' spread), eight draws each, with the true variogram and with the estimated
#' one (96 runs, 2,978 cells with points among them).  \code{kr_var} was
#' below \code{plain_var} in every such cell, and \code{kr_gain} ran from
#' 0.14 to 0.86 (median 0.48 uniform, 0.52 clustered).  Against the true
#' cell averages kriging cut the plain mean's squared error by 54--56
#' percent (its RMSE by 33 percent) on both layouts, where mean
#' \code{kr_var} over mean \code{plain_var} predicted 55--58 percent.  In
#' the third of cells with the lowest \code{kr_gain} the cut in squared
#' error was 39--49 percent and in the third with the highest 59--66
#' percent, in all four combinations of layout and variogram, so the column
#' ranks cells by what kriging buys.  \code{kr_shift} had a mean square of
#' 0.99--1.10, with 33--36 percent of cells beyond 1 and 4--6 percent beyond
#' 2.
#'
#' What the layout changed was the naive standard error and the empty cells.
#' Under uniform sampling \eqn{s^2/n} was below \code{plain_var} in 52--57
#' percent of cells with two or more points (median ratio 0.88--0.96), and at
#' most two cells were empty.  Under clustered sampling it was below in
#' 76--98 percent (median ratio 0.19--0.58), and 1--5 of 16 cells, 8--17 of
#' 36 and 23--36 of 64 were empty, each with a kriged estimate and variance
#' where the plain mean has nothing.  So on these fields kriging removed
#' about half of the plain mean's error variance on either layout;
#' clustering decides how far the naive standard error can be trusted and
#' how many cells have no plain mean at all.
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
#' extrapolation; a variogram that falls with distance, or a range below the
#' shortest lag, describes the data poorly at some lags; and a response in
#' which \code{estimate_sac_range()} found no spatial structure comes with a
#' pure-nugget model, which gives kriging no correlation to use.
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
#' \eqn{s_w^2} is the part that differs from visit to visit and averages
#' down, and what the nugget holds beyond it is micro-scale variation the
#' visits share, so the mean carries error variance
#' \eqn{\max(c_0 - s_w^2, 0) + s_w^2/m}.  That goes to \pkg{gstat} as a known
#' measurement error (its \code{weights}) on the model with the nugget set
#' to zero.  For replicates that differ only by measurement error this is
#' exactly the kriging of every observation, and for identical replicates it
#' is the kriging of one.  \code{plain_var} uses the same error variances:
#' the observations a cell has at one location share the micro-scale part
#' and average the rest.
#'
#' When \eqn{s_w^2} is above \eqn{c_0} the replicates show more variation
#' between visits than the model has room for (with a zero nugget, any at
#' all).  The replicate variance is then used as it is: a mean of \eqn{m} is
#' kriged with error variance \eqn{s_w^2/m}, and a location seen once with
#' \eqn{s_w^2} where the model gives it \eqn{c_0}.  Capping it at
#' \eqn{c_0} would krige the means as more exact than their replicates show
#' them to be: on ten simulated layers of 200 stations visited twice with a
#' measurement variance of 0.5, under the field's own variogram with no
#' nugget, the cross-validation statistic reads 2.96--4.42 with the cap and
#' 0.84--1.15 without it.  The warning then carries the class
#' \code{"spatialkit_replicate_var_above_nugget"} and both numbers, which
#' \code{attr(, "replicates")} records and \code{print()} repeats.  The
#' replicate variance is an estimate (its degrees of freedom are in the
#' warning and the record), and the variogram still says \eqn{c_0}: refit
#' it with a nugget of at least \eqn{s_w^2}.
#'
#' In the cross-validation a location is held out whole, under the fold of
#' its first row, and its error variance is part of its standardised error.
#' A cell or held-out location \pkg{gstat} still cannot krige is reported
#' \code{NA} with a warning, and \code{print()} says how many.
#'
#' @param assigned_points_sf Points with a cell identifier column, as
#'   \code{\link{assign_features_to_polygons}()} returns.
#' @param response_var The response column: numeric or logical. A logical
#'   column is kriged as 0/1, so its cell mean is a share. Any other class
#'   (a factor, an ordered factor, character, a date) is an error naming the
#'   column and its class, and is never converted: a factor's level codes are
#'   an arbitrary relabelling, and kriging them gave a variogram, kriged
#'   estimates and plain means that all changed with the order of the
#'   levels. For a category, pass an indicator (\code{x == "level"}).
#' @param cells_sf The cell polygons, with the matching ID column.  A layer
#'   with no CRS is taken to be in the points' CRS (and points with none in
#'   the cells'), with a warning.
#' @param id_col Preferred name of the ID column, found as
#'   \code{\link{summarize_by_cell}()} finds it: on the points the first of
#'   \code{id_col}, \code{"poly_id"}, \code{"polygon_id"} and
#'   \code{"cell_id"}; on the cells the first of that column,
#'   \code{"poly_id"}, \code{"polygon_id"}, \code{"id"}, \code{"cell_id"}
#'   and \code{"grid_id"}.  IDs are matched as text, whole numbers written
#'   out in full, so a double \code{1e5} matches an integer \code{100000};
#'   a point whose ID matches no cell is counted in no cell, with a
#'   warning.  Default \code{"poly_id"}.
#' @param sac Optional \code{sac_range} carrying a variogram model.
#' @param folds Optional \code{\link{make_folds}()} result on
#'   \code{assigned_points_sf} for the cross-validation statistic; built here
#'   with \code{block_kfold} when \code{NULL}.  As in \code{cv_*()}, folds
#'   whose recorded rows sit at other locations in
#'   \code{assigned_points_sf} (built on another layer, such as the points
#'   before \code{\link{assign_features_to_polygons}()} dropped some) are
#'   refused.
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
#'   error), \code{plain_var} (the plain mean's error variance under the
#'   variogram), \code{kr_pred}, \code{kr_var}, \code{kr_ratio},
#'   \code{kr_gain}, \code{kr_shift} and \code{kr_n_used} (how many
#'   locations the cell was kriged from; \code{NA} for a cell left out).
#'   Attributes:
#'   \code{variogram} (the model frame), \code{sill}, \code{nugget},
#'   \code{range}, \code{range_identified}, \code{rejected_reason} (why the
#'   range was refused, from \code{sac}; \code{NA} when it was identified),
#'   \code{cv} (a list:
#'   \code{zscore_var}, \code{zscore_mean}, \code{rmse}, \code{n_pred},
#'   \code{k}, \code{method}), \code{nmax}, \code{n_points} (the points
#'   used), \code{n_locations} (the distinct locations among them, which
#'   the kriging used), \code{cells_left_out} (counts of cells left out,
#'   named \code{shape} for \code{max_box_ratio} and \code{size} for
#'   \code{max_neighbours}) and \code{replicates} (\code{NULL} when no two
#'   points share a location; otherwise a list with \code{n_points}, the
#'   points that share one, \code{within_var}, their pooled within-location
#'   variance, \code{df}, its degrees of freedom, \code{nugget}, the
#'   variogram's, \code{nugget_used}, the larger of the two, which is the
#'   error variance a location seen once is kriged with, and
#'   \code{above_nugget}, \code{TRUE} when the replicate variance exceeds
#'   the variogram's nugget).
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
  # The response's type, before anything converts it.  as.numeric() under
  # suppressWarnings() used to be the first thing done to the column, so a
  # factor was kriged on its level codes (and its plain means were means of
  # codes): on one three-level factor, reordering the levels moved the nugget
  # 0.102 / 0.436 / 0.271 and changed which cells were flagged, with nothing
  # said.  A character column of numbers and a Date were converted the same
  # way, and a character column of labels ended in "fewer than 10 points with
  # a finite response", which blames the data.  estimate_sac_range() refuses
  # a factor, but it only ever saw the converted copy, and is not called at
  # all when `sac` is supplied.
  resp_raw <- sf::st_drop_geometry(assigned_points_sf)[[response_var]]
  if (!(is.numeric(resp_raw) || is.logical(resp_raw)))
    stop(sprintf(paste0(
      "kriging_adequacy(): response '%s' is %s, and kriging needs a numeric or ",
      "logical variable. %s For a category, pass an indicator column ",
      "(%s == \"<level>\"): it is kriged as 0/1 and its cell mean is that ",
      "level's share."),
      response_var,
      if (is.factor(resp_raw)) "a factor"
      else sprintf("of class %s", paste(class(resp_raw), collapse = "/")),
      if (is.factor(resp_raw))
        paste0("A factor's codes are an arbitrary relabelling of its levels: the ",
               "variogram, the kriged estimates and the plain cell means would ",
               "all change when the levels are reordered.")
      else
        paste0("It is not converted for you; convert it with as.numeric() ",
               "first if its numbers are what you mean."),
      response_var), call. = FALSE)
  if (!is.numeric(nmax) || length(nmax) != 1L || !is.finite(nmax) || nmax < 1)
    stop("kriging_adequacy(): `nmax` must be a single positive number.", call. = FALSE)

  # --- the ID column, as summarize_by_cell() finds it ---
  # The same candidates in the same order: on the points the columns
  # assign_features_to_polygons() writes, and on the cells, after the points'
  # own column, every one it reads the polygons' IDs from.  Cells keyed by
  # 'id' or 'grid_id' were refused here while summarize_by_cell() joined them.
  pts_candidates <- unique(c(id_col, "poly_id", "polygon_id", "cell_id"))
  id_pts <- pts_candidates[pts_candidates %in% names(assigned_points_sf)]
  if (!length(id_pts))
    stop("kriging_adequacy(): could not find a cell ID column in ",
         "`assigned_points_sf`. Looked for: ",
         paste(pts_candidates, collapse = ", "), ".", call. = FALSE)
  id_pts <- id_pts[[1L]]
  cells_candidates <- unique(c(id_pts, "poly_id", "polygon_id", "id", "cell_id",
                               "grid_id"))
  id_cells <- cells_candidates[cells_candidates %in% names(cells_sf)]
  if (!length(id_cells))
    stop("kriging_adequacy(): could not find a cell ID column in `cells_sf`. ",
         "Looked for: ", paste(cells_candidates, collapse = ", "), ".",
         call. = FALSE)
  id_cells <- id_cells[[1L]]

  if (!is.numeric(max_neighbours) || length(max_neighbours) != 1L ||
      !is.finite(max_neighbours) || max_neighbours < 1)
    stop("kriging_adequacy(): `max_neighbours` must be a single positive number.", call. = FALSE)
  if (!is.numeric(max_box_ratio) || length(max_box_ratio) != 1L ||
      is.na(max_box_ratio) || max_box_ratio < 1)
    stop("kriging_adequacy(): `max_box_ratio` must be a single number of at least 1.",
         call. = FALSE)

  # Folds are row IDs, matched by position on a layer without `..row_id`.
  # Built on prep_model_data()'s points and applied to the fewer that
  # assign_features_to_polygons() kept, they held out the wrong points (RMSE
  # 1.82 against 2.11 with folds built here) with nothing said, while cv_*()
  # refuse the same folds through the probe make_folds() records.  Checked on
  # the layer as passed, before anything is moved or pointized.
  .check_fold_provenance(folds, assigned_points_sf, "kriging_adequacy")

  # --- geometry: projected points, cells in the same CRS ---
  # A layer with no CRS is taken to be in the other's, as passed and before
  # anything is moved: the points were assigned to these cells, so they share
  # a coordinate space.  .align_crs() leaves a CRS-less layer as it is, and
  # gstat::krige() then failed on the mismatch with an internal assertion.
  pts <- assigned_points_sf
  cells_in <- cells_sf
  pts_na   <- is.na(sf::st_crs(pts))
  cells_na <- is.na(sf::st_crs(cells_in))
  if (cells_na && !pts_na)
    cells_in <- .transform_or_stamp(cells_in, sf::st_crs(pts), what = "cells_sf",
                                    caller = "kriging_adequacy")
  else if (pts_na && !cells_na)
    pts <- .transform_or_stamp(pts, sf::st_crs(cells_in),
                               what = "assigned_points_sf",
                               caller = "kriging_adequacy")
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
  # Both without a CRS, and the points put in the variogram's: the cells
  # follow them the same way.
  cells <- if (is.na(sf::st_crs(cells_in)) && !is.na(sf::st_crs(pts)))
    .transform_or_stamp(cells_in, sf::st_crs(pts), what = "cells_sf",
                        caller = "kriging_adequacy")
  else .align_crs(cells_in, pts)
  cells <- .safe_make_valid(cells)
  # Numeric or logical by the check above, so nothing is coerced to NA here.
  z <- as.numeric(sf::st_drop_geometry(pts)[[response_var]])
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
  sac_supplied <- !is.null(sac)
  if (!sac_supplied) {
    .msg("kriging_adequacy(): estimating the variogram of ", response_var, " ...")
    sac <- estimate_sac_range(kp, "..z", seed = seed)
  }
  vm <- attr(sac, "variogram_model")
  if (is.null(vm) || !is.data.frame(vm)) {
    # Why there is none: a classed refusal says (both fits singular), and a
    # bare NA is a run that could fit nothing.  Blaming `sac` when none was
    # passed sent the user looking at an argument they had not used.
    why <- attr(sac, "rejected_reason")
    why <- if (is.character(why) && length(why) == 1L && !is.na(why)) why
      else "estimate_sac_range() returned a bare NA"
    stop(sprintf(paste0("kriging_adequacy(): %s carries no variogram model, so ",
                        "there is nothing to krige with: %s."),
                 if (sac_supplied) "`sac`"
                 else sprintf("the variogram estimated here from '%s'", response_var),
                 why), call. = FALSE)
  }
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
  # said of every refusal, and a variogram that falls with distance, one that
  # runs past the fitted lags and one with no spatial structure are different
  # findings about the model the variances rest on.
  rejected_reason <- if (range_identified) NA_character_ else
    as.character(attr(sac, "rejected_reason") %||% "no effective range")[1L]
  if (!range_identified) {
    # A layer with no spatial structure comes with a pure nugget
    # (.sac_no_structure(), R/utils.R).  That model is what the data support,
    # so "one the data did not pin down" was the wrong thing to say of it.
    if (.sac_no_structure(sac) && all(fam == "Nug"))
      .warn_and_log(paste0("kriging_adequacy(): the variogram's range was not ",
                           "identified (%s), and the model that comes with ",
                           "that verdict is a pure nugget: no spatial ",
                           "correlation was found for kriging to use."),
                    rejected_reason)
    else
      .warn_and_log(paste0("kriging_adequacy(): the variogram's range was not ",
                           "identified (%s), so the model the kriging variances ",
                           "and ratios rest on is one the data did not pin down."),
                    rejected_reason)
  }
  # A variogram of residuals on predictors describes what the predictors left
  # over, and it is the response itself that is kriged here, with no
  # predictors: the variance they explained is missing from the model.
  # Measured with a spatially structured covariate: the variance of the
  # standardised CV errors went from 0.67-1.53 to 4.3-5.2.
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
  # between visits and averages down by the count; what the nugget holds
  # beyond it is micro-scale variation the visits share.  That goes to gstat
  # as a known measurement error (`weights` = 1 / variance) on the model with
  # its nugget zeroed: a location seen once gets the whole nugget (or the
  # replicate variance, when that is the larger); replicates differing only
  # by measurement error give the kriging of every observation; identical
  # replicates give the kriging of one.
  xy  <- sf::st_coordinates(kp)[, 1:2, drop = FALSE]
  key <- paste(xy[, 1], xy[, 2])
  loc <- match(key, unique(key))
  first <- !duplicated(loc)
  kk <- kp[first, , drop = FALSE]
  vk <- vm
  w <- NULL
  err_var <- rep(0, nrow(kk))    # already in vk's nugget when w is NULL
  replicates <- NULL             # the record of the step below, for the result
  # The two parts of the nugget, as the plain mean's variance below needs
  # them: `micro` is shared by every observation at a location, `me` differs
  # from one observation to the next.  With one observation per location the
  # two cannot be told apart and do not need to be: the nugget is then that
  # observation's own.
  micro <- nugget
  me    <- 0
  if (any(!first)) {
    m_loc <- tabulate(loc)
    kk$..z <- as.numeric(rowsum(kp$..z, loc)) / m_loc
    # The replicates' own estimate of what differs between visits.  It used to
    # be capped at the nugget, the most the model lets replicates differ by:
    # 40 stations visited twice with a replicate variance of 0.09, under a
    # model with no nugget, were reported as "0 of 0" and kriged as exact
    # (CV z-variance 1.61, where the uncapped variance below gives 1.32).  The
    # replicates measure the variance between visits directly, so it is used
    # as it is; what the nugget holds beyond it is the micro-scale part, and
    # when the replicates exceed the nugget that part is zero.  The rule is
    # continuous where the two meet: adding the replicate variance on top of
    # the whole nugget instead would take a mean of three from c0 / 3 to
    # c0 + c0 / 3 as the replicate variance crosses c0.
    df_w <- sum(m_loc - 1L)
    sw2  <- sum((kp$..z - kk$..z[loc])^2) / df_w
    me    <- sw2
    micro <- max(nugget - sw2, 0)
    above <- is.finite(sw2) && sw2 > nugget + 1e-8 * sill
    err_var <- micro + me / m_loc
    vk$psill[fam == "Nug"] <- 0
    w <- 1 / err_var
    n_rep <- sum(m_loc[loc] > 1L)
    replicates <- list(n_points = n_rep, within_var = sw2, df = as.integer(df_w),
                       nugget = nugget, nugget_used = max(nugget, sw2),
                       above_nugget = above)
    if (above)
      .warn_classed("spatialkit_replicate_var_above_nugget", paste0(
        "kriging_adequacy(): %d point(s) share a location with another; ",
        "kriging from the mean at each of the %d distinct locations. Their ",
        "pooled within-location variance is %.3g (%d degrees of freedom), ",
        "above the variogram's nugget of %.3g, so the model has less room than ",
        "the replicates show between visits. The replicate variance is used ",
        "as the measurement-error variance: a location seen m times is kriged ",
        "with error variance %.3g / m, a nugget of %.3g for a location seen ",
        "once, and kr_var, plain_var and the cross-validation statistic are ",
        "computed with it. Refit the variogram with a nugget of at least %.3g ",
        "so the model says the same; attr(, \"replicates\") records both ",
        "numbers."),
        n_rep, nrow(kk), sw2, as.integer(df_w), nugget, sw2, sw2, sw2,
        data = list(within_var = sw2, nugget = nugget))
    else
      .warn_and_log(paste0("kriging_adequacy(): %d point(s) share a location with another; ",
                           "kriging from the mean at each of the %d distinct locations, with ",
                           "the part of the nugget that differs between replicates (%.3g of ",
                           "%.3g) divided by their count."),
                    n_rep, nrow(kk), me, nugget)
  }

  # summarize_by_cell()'s conversion: as.character() writes a double from 1e5
  # on as "1e+05", so every such point missed its integer-keyed cell and the
  # cell came back with n = 0, no plain mean, and a kriging neighbourhood
  # without its own points.  Said, not silent, when an ID still matches no
  # cell, as summarize_by_cell() says it.
  ids_pts <- .id_as_character(sf::st_drop_geometry(kp)[[id_pts]])
  ids_cells <- .id_as_character(sf::st_drop_geometry(cells)[[id_cells]])
  unmatched <- !is.na(ids_pts) & !(ids_pts %in% ids_cells)
  if (any(unmatched))
    .warn_and_log(paste0("kriging_adequacy(): %d of the %d point(s) with a cell ID ",
                         "carry one that matches no `cells_sf$%s` (%s), so they ",
                         "are in no cell's n, mean or se, though the kriging still ",
                         "uses them. Check that `cells_sf` is the layer the points ",
                         "were assigned to."),
                  sum(unmatched), sum(!is.na(ids_pts)), id_cells,
                  paste(utils::head(unique(ids_pts[unmatched]), 5L), collapse = ", "))

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
  #
  # And the plain mean's own error variance under the same model, plain_var:
  # what kr_var has to be set against to say what kriging buys in a cell.
  # kr_var used to be compared with s^2/n, the variance the plain mean would
  # have if the cell's points were independent draws, as `kr_exceeds_design`.
  # That flag was read as marking cells where kriging does not pay, and it
  # does not: against known cell averages kriging cut the error by 20-53% in
  # flagged cells and 27-60% in the rest.  It was not on one footing either,
  # one variance from the model and one from the cell's own spread: the same
  # data with identical kriged estimates flagged 4, 27 and 67 of 70 cells as
  # the sill went x1/3, x1 and x3.  Both variances now come from the one
  # variogram and the one discretisation of the cell.  Computed for the cells
  # that are kriged and hold a point; the cells left out are left out of
  # this as well (a sliver for its discretisation, an oversized cell for the
  # pairs of its own locations).
  xy_loc  <- xy[first, , drop = FALSE]
  rows_of <- split(seq_len(nrow(kp)),
                   factor(match(ids_pts, ids_cells), levels = seq_len(nrow(cells))))
  own_pts <- vector("list", nrow(cells))
  for (ci in which(kriged & lengths(rows_of) > 0L)) {
    l   <- loc[rows_of[[ci]]]
    u   <- unique(l)
    cnt <- tabulate(match(l, u), nbins = length(u))
    # A location's `cnt` observations in the cell share its micro-scale part
    # and average their own measurement errors: the error variance the
    # kriging gives a mean of that many, and the whole nugget for one
    # observation at a location seen once.
    own_pts[[ci]] <- list(xy = xy_loc[u, , drop = FALSE], a = cnt / length(l),
                          e = micro + me / cnt)
  }
  prior_var <- plain_var <- rep(NA_real_, nrow(cells))
  if (any(!sliver)) {
    bv <- .block_model_vars(cells[!sliver, , drop = FALSE], vm, own_pts[!sliver])
    prior_var[!sliver] <- bv$prior
    plain_var[!sliver] <- bv$plain
  }

  # --- the plain means and their naive standard error, per cell ---
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
  out$plain_var <- plain_var
  out$kr_pred <- kr_pred
  out$kr_var <- kr_var
  out$kr_ratio <- ifelse(prior_var > 0, pmin(pmax(kr_var / prior_var, 0), 1), NA_real_)
  # The share of the plain mean's error variance that kriging removes.  gstat's
  # block variance and the same variance written out here agree to about
  # 2e-6 of their size, so a gain within a ten-thousandth below 0 is 0.
  gain <- ifelse(is.finite(plain_var) & plain_var > 0, 1 - kr_var / plain_var, NA_real_)
  gain[is.finite(gain) & gain < 0 & gain > -1e-4] <- 0
  out$kr_gain <- gain
  # The kriged estimate minus the plain mean, over the standard deviation the
  # model gives that difference.  The kriging error is uncorrelated with any
  # difference of two unbiased combinations of the neighbourhood's data, so
  # Var(plain - truth) = Var(kriged - truth) + Var(kriged - plain), and the
  # difference has variance plain_var - kr_var.  Over s/sqrt(n), as it was,
  # its variance on clustered layouts ran 1.4-8.4 against a nominal 1, so
  # the printed count of shifts beyond 2 reached 40% of cells where 3% is
  # the figure on the model's own scale.  Not defined where the two
  # variances agree to within a ten-thousandth, which is where kriging
  # returns the plain mean and both terms of the ratio are rounding error.
  shift_var <- plain_var - kr_var
  out$kr_shift <- ifelse(is.finite(shift_var) & shift_var > 1e-4 * plain_var,
                         (kr_pred - out$mean) / sqrt(pmax(shift_var, 0)), NA_real_)

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
            replicates = replicates,
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


#' Two variances of a cell's mean under the variogram, on one discretisation
#'
#' \code{prior} is the variance each cell's mean has with no data,
#' \eqn{\bar C(B,B)}: the covariance averaged over all pairs of points of the
#' cell, on the discretisation gstat's predict() gives a polygon (its default
#' `sps.args`: spsample(n = 500, type = "regular", offset = c(0.5, 0.5))),
#' and without the nugget, which gstat leaves out of block-to-block
#' covariances.  This is the block-kriging variance of a cell no datum
#' informs, before ordinary kriging adds the variance of the estimated mean.
#' 0 for a pure nugget.
#'
#' \code{plain} is the error variance of the plain mean of the cell's own
#' observations as an estimate of the cell's average: with the cell's
#' locations \eqn{p_l}, the share \eqn{a_l} of the cell's observations at
#' each, and the error variance \eqn{e_l} of the mean of those (the nugget
#' for one observation; see the caller for replicates),
#' \deqn{V = \sum_{l,l'} a_l a_{l'} C(p_l - p_{l'}) + \sum_l a_l^2 e_l
#'       - 2 \sum_l a_l \bar C(p_l, B) + \bar C(B,B),}
#' with \eqn{C} the covariance without the nugget, which sits on the data
#' diagonal only, as gstat block-kriges it.  The plain mean is one unbiased
#' linear combination of the data in the cell's kriging neighbourhood, and
#' ordinary kriging is the best of them, so the block-kriging variance of a
#' cell never exceeds this.  Both come from the same 500 points.
#'
#' @param cells The cells (none of them a sliver).
#' @param vm The variogram model frame.
#' @param own Optional list, one entry per cell: \code{NULL} for a cell
#'   whose plain variance is not wanted, else a list with \code{xy} (a
#'   two-column matrix of the cell's distinct locations), \code{a} (their
#'   weights in the plain mean, summing to 1) and \code{e} (their error
#'   variances).
#' @return A list with \code{prior} and \code{plain}, one value per cell;
#'   \code{plain} is \code{NA} where \code{own} gives nothing.
#' @keywords internal
#' @noRd
.block_model_vars <- function(cells, vm, own = NULL) {
  nc <- nrow(cells)
  sv <- vm[as.character(vm$model) != "Nug", , drop = FALSE]
  c0 <- sum(sv$psill)
  structured <- nrow(sv) > 0L && isTRUE(c0 > 0)
  has_own <- if (is.null(own)) rep(FALSE, nc) else !vapply(own, is.null, logical(1))
  prior <- rep(0, nc)
  plain <- rep(NA_real_, nc)
  # A pure nugget: no covariance between distinct places, so a cell mean has
  # no prior variance and the plain mean's error is its data's own.
  if (!structured) {
    for (i in which(has_own)) plain[i] <- sum(own[[i]]$a^2 * own[[i]]$e)
    return(list(prior = prior, plain = plain))
  }
  covf <- function(d) gstat::variogramLine(sv, dist_vector = d, covariance = TRUE)$gamma
  # sum_ij wa_i wb_j C(|A_i - B_j|), a block of rows at a time so a cell of
  # thousands of locations does not hold all its pairs at once.
  wsum <- function(A, wa, B, wb) {
    step <- max(1L, as.integer(2e6 %/% nrow(B)))
    tot <- 0
    for (s in seq(1L, nrow(A), by = step)) {
      i <- s:min(nrow(A), s + step - 1L)
      d <- sqrt(outer(A[i, 1L], B[, 1L], "-")^2 + outer(A[i, 2L], B[, 2L], "-")^2)
      cv <- matrix(covf(as.numeric(d)), nrow = length(i))
      tot <- tot + sum(wa[i] * as.numeric(cv %*% wb))
    }
    tot
  }
  geo <- sf::as_Spatial(sf::st_set_crs(sf::st_geometry(cells), sf::NA_crs_))
  for (i in seq_len(nc)) {
    g <- sp::coordinates(sp::spsample(geo[i], n = 500, type = "regular",
                                      offset = c(0.5, 0.5)))
    m <- nrow(g)
    prior[i] <- if (m < 2L) c0 else {
      # Each pair once from dist(), both orders, plus the m zero-distance pairs.
      cv <- covf(as.numeric(stats::dist(g)))
      (2 * sum(cv) + m * c0) / m^2
    }
    if (has_own[i] && m >= 1L) {
      o <- own[[i]]
      plain[i] <- wsum(o$xy, o$a, o$xy, o$a) + sum(o$a^2 * o$e) -
        2 * wsum(o$xy, o$a, g[, 1:2, drop = FALSE], rep(1 / m, m)) + prior[i]
    }
  }
  list(prior = prior, plain = pmax(plain, 0))
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
  # Kriging against the plain mean, both variances from the one variogram.
  # A subset or an object from before these columns existed has none.
  gn <- df$kr_gain[is.finite(df$kr_gain)]
  if (length(gn))
    cat(sprintf(paste0("  kr_gain, the share of the plain mean's error variance that kriging ",
                       "removes under this variogram (1 - kr_var / plain_var): median %.2f, ",
                       "range %.2f-%.2f over the %d cell(s) with points and a kriged ",
                       "estimate\n"),
                stats::median(gn), min(gn), max(gn), length(gn)))
  sh <- df$kr_shift[is.finite(df$kr_shift)]
  if (length(sh))
    cat(sprintf(paste0("  kriged minus plain mean, over the standard deviation the variogram ",
                       "gives that difference (sqrt(plain_var - kr_var)): |shift| > 1 in %d ",
                       "of %d cell(s), > 2 in %d (about %.0f and %.0f when the variogram is ",
                       "right)\n"),
                sum(abs(sh) > 1), length(sh), sum(abs(sh) > 2),
                0.317 * length(sh), 0.046 * length(sh)))
  # The naive standard error against the same model: s^2/n is the plain
  # mean's variance if the cell's points were independent draws.
  nv <- (df$se^2 / df$plain_var)[is.finite(df$se) & is.finite(df$plain_var) &
                                   df$plain_var > 0]
  if (length(nv))
    cat(sprintf(paste0("  naive s^2/n below plain_var (the naive standard error is too ",
                       "small if the variogram is right) in %d of the %d cell(s) with two ",
                       "or more points; median s^2/n over plain_var %.2f\n"),
                sum(nv < 1), length(nv), stats::median(nv)))
  rp <- attr(x, "replicates")
  if (is.list(rp) && isTRUE(rp$above_nugget))
    cat(sprintf(paste0("  replicates: pooled within-location variance %.3g is above the ",
                       "nugget %.3g and is used as the measurement-error variance, a nugget ",
                       "of %.3g for a location seen once\n"),
                rp$within_var, rp$nugget, rp$nugget_used))
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
                     knndm              = "kNNDM k-fold CV",
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
