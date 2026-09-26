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
#' \code{gstat::krige()} on the fold labels \code{\link{make_folds}()} built,
#' so the folds carry the same separation the package uses everywhere else.
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
#' model whose range was not identified (a bare \code{NA} estimate with the
#' model attached) is used with a warning: its sill was never reached by the
#' data, so the ratios rest on an extrapolation.  Requires \pkg{gstat}.
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
#' @param nmax The largest number of neighbours each kriging system uses
#'   (\code{gstat}'s \code{nmax}).  Default 50.
#' @param quiet Suppress progress messages.  Default \code{TRUE}.
#' @return An \code{sf} object of class \code{"kriging_adequacy"}, one row
#'   per cell with the cell geometry and: the ID column, \code{n} (points in
#'   the cell), \code{mean} (the plain mean), \code{se} (its naive standard
#'   error), \code{kr_pred}, \code{kr_var}, \code{kr_ratio},
#'   \code{kr_exceeds_design} and \code{kr_shift}.  Attributes:
#'   \code{variogram} (the model frame), \code{sill}, \code{nugget},
#'   \code{range}, \code{range_identified}, \code{cv} (a list:
#'   \code{zscore_var}, \code{zscore_mean}, \code{rmse}, \code{n_pred},
#'   \code{k}, \code{method}), \code{nmax}, \code{n_points} (the points
#'   used) and \code{n_locations} (the distinct locations among them, which
#'   the kriging used).
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
                             k = 5L, seed = 123L, nmax = 50L, quiet = TRUE) {
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

  # --- geometry: projected points, cells in the same CRS ---
  pts <- assigned_points_sf
  if (!all(sf::st_geometry_type(pts, by_geometry = TRUE) == "POINT"))
    pts <- coerce_to_points(pts, "auto")
  pts <- ensure_projected(pts)
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
  if (!range_identified)
    .warn_and_log(paste0("kriging_adequacy(): the variogram's range was not ",
                         "identified (%s), so its sill was never reached by the ",
                         "data; the kriging variances and ratios rest on an ",
                         "extrapolated model."),
                  attr(sac, "rejected_reason") %||% "no effective range")

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

  # --- block kriging onto the cells ---
  .msg("kriging_adequacy(): block kriging onto ", nrow(cells), " cells ...")
  bk <- tryCatch(
    gstat::krige(..z ~ 1, locations = kk, newdata = cells, model = vk,
                 weights = w, nmax = as.integer(nmax), debug.level = 0),
    error = function(e)
      stop("kriging_adequacy(): gstat::krige() failed: ", conditionMessage(e),
           call. = FALSE))
  kr_pred <- suppressWarnings(as.numeric(bk$var1.pred))
  kr_var  <- suppressWarnings(as.numeric(bk$var1.var))
  kr_var[is.finite(kr_var) & kr_var < 0] <- 0
  # gstat answers a system it cannot solve with NA and, at debug.level 0,
  # nothing else.
  n_na <- sum(!is.finite(kr_pred) | !is.finite(kr_var))
  if (n_na)
    .warn_and_log(paste0("kriging_adequacy(): gstat::krige() returned no estimate for %d ",
                         "of %d cell(s) (a kriging system it could not solve); their ",
                         "kr_ columns are NA."), n_na, length(kr_pred))
  # The scale for kr_ratio.  kr_var is the variance of a cell MEAN, and with
  # no data it levels off at the cell's own C(B,B), not at the point sill:
  # divided by the sill, an empty cell far from every datum read about 0.1,
  # and a large empty cell ranked below small populated ones.
  prior_var <- .block_prior_var(cells, vm)

  # --- the plain means and their design-based variance, per cell ---
  ids_pts <- as.character(sf::st_drop_geometry(kp)[[id_pts]])
  ids_cells <- as.character(sf::st_drop_geometry(cells)[[id_cells]])
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

  # --- the cross-validation statistic ---
  if (is.null(folds)) {
    folds <- make_folds(kp, k = k, method = "block_kfold", seed = seed)
    fold_vec <- folds$assignment$fold[match(kp$..row_id, folds$assignment$row_id)]
    method <- "block_kfold"
  } else {
    fold_vec <- .fold_labels_for(folds, kp)
    method <- if (is.list(folds) && !is.null(folds$method)) folds$method else "supplied labels"
  }
  # A location is held out whole, under the fold of its first row: a
  # replicate predicted from its twin in the training folds is no test.
  fold_vec <- fold_vec[first]
  cv <- list(zscore_var = NA_real_, zscore_mean = NA_real_, rmse = NA_real_,
             n_pred = 0L, k = length(unique(fold_vec[is.finite(fold_vec)])),
             method = method)
  keep <- is.finite(fold_vec)
  if (sum(keep) >= 10L && cv$k >= 2L) {
    .msg("kriging_adequacy(): cross-validating the kriging variance over ", cv$k, " folds ...")
    # The folds are run here rather than by gstat::krige.cv(), which subsets
    # the data per fold but not `weights`.  With the nugget zeroed in vk, a
    # held-out location's own error variance goes back into its variance.
    kc <- kk[keep, , drop = FALSE]; fk <- fold_vec[keep]
    wc <- w[keep]; ec <- err_var[keep]
    kcv <- tryCatch({
      pred <- pvar <- rep(NA_real_, nrow(kc))
      for (f in unique(fk)) {
        te <- fk == f
        kf <- gstat::krige(..z ~ 1, locations = kc[!te, , drop = FALSE],
                           newdata = kc[te, , drop = FALSE], model = vk, weights = wc[!te],
                           nmax = as.integer(nmax), debug.level = 0)
        pred[te] <- suppressWarnings(as.numeric(kf$var1.pred))
        pvar[te] <- suppressWarnings(as.numeric(kf$var1.var)) + ec[te]
      }
      data.frame(residual = kc$..z - pred,
                 zscore = suppressWarnings((kc$..z - pred) / sqrt(pvar)))
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
      if (cv$n_pred < nrow(kc))
        .warn_and_log(paste0("kriging_adequacy(): gstat::krige() returned no cross-validation ",
                             "prediction it could use for %d of %d held-out location(s); ",
                             "the statistic uses the other %d."),
                      nrow(kc) - cv$n_pred, nrow(kc), cv$n_pred)
    }
  } else {
    .log_warn("kriging_adequacy(): too few points or folds for the cross-validation statistic.")
  }

  structure(out,
            variogram = vm, sill = sill, nugget = nugget,
            range = suppressWarnings(as.numeric(sac)),
            range_identified = range_identified,
            cv = cv, nmax = as.integer(nmax), n_points = nrow(kp),
            n_locations = nrow(kk),
            response_var = response_var,
            class = c("kriging_adequacy", class(out)))
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
                sprintf("%.1f", attr(x, "range")) else "not identified"))
  # Every line below counts only the cells with an estimate, so say first
  # how many have none.
  kr_ok <- is.finite(df$kr_pred) & is.finite(df$kr_var)
  if (!all(kr_ok))
    cat(sprintf(paste0("  no kriged estimate for %d of the %d cells (gstat could not ",
                       "solve their kriging systems)\n"), sum(!kr_ok), nrow(df)))
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
  if (is.finite(cv$zscore_var %||% NA_real_))
    cat(sprintf(paste0("  blocked CV (%s, %d folds, %d %s): var of standardised ",
                       "error %.2f (1 = kriging variance correct; above 1 = ",
                       "understated), mean %.2f, RMSE %.3g\n"),
                cv$method, cv$k, cv$n_pred, pts_word, cv$zscore_var, cv$zscore_mean, cv$rmse))
  else cat("  blocked CV: not computed\n")
  invisible(x)
}
