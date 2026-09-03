#' Assign features to polygons and attach a polygon ID
#'
#' Joins an sf layer of input features to a polygon layer via spatial join.
#'
#' This is the second step of the package's pipeline: it labels every
#' observation with the cell it falls in, which is what
#' [summarize_by_cell()] then aggregates over. Reach for it directly (rather
#' than for [sf::st_join()]) when the join has to be *unambiguous* --- it
#' resolves features matching several polygons by an explicit `tie_break` rule
#' instead of silently duplicating rows, so the assigned layer keeps one row
#' per input feature and cell-level counts mean what they say.
#'
#' @param features_sf An sf object containing features to assign.
#' @param polygons_sf An sf or sfc polygonal layer.
#' @param polygon_id_col Name of the polygon identifier column. Default "poly_id".
#' @param keep_unassigned Logical; retain features that fall inside no
#'   polygon, carrying \code{NA} in the ID column. Default FALSE, which drops
#'   them.
#' @param predicate Binary spatial predicate function. Default sf::st_intersects.
#' @param largest Logical; when \code{features_sf} is itself polygonal, keep
#'   the polygon with the largest overlap. Default TRUE. Ignored for point and
#'   line features, and silently dropped if the \code{predicate} does not
#'   support it (\code{sf::st_intersects} does).
#' @param tie_break Strategy for resolving features that match multiple
#'   polygons: \code{"smallest_area"} (default) keeps the polygon with the
#'   smallest area, \code{"first"} keeps the first match (original order-dependent
#'   behavior).
#' @return An sf object with `polygon_id_col` attached, one row per input
#'   feature (fewer if `keep_unassigned = FALSE` dropped unmatched ones), in
#'   the CRS `features_sf` arrived in. Any column of `features_sf` whose name
#'   would collide with the polygon ID column is dropped before the spatial
#'   join (with a warning), so re-assigning an already-assigned layer replaces
#'   the old IDs rather than failing. If *no* feature falls inside any polygon
#'   the result is empty (or all-`NA` with `keep_unassigned = TRUE`) and a
#'   warning is raised, since the usual cause is two layers in different
#'   places --- a CRS that could only be stamped, not reprojected.
#' @examples
#' library(sf)
#' set.seed(1)
#' pts <- st_as_sf(
#'   data.frame(x = runif(20, 0, 100), y = runif(20, 0, 100), val = rnorm(20)),
#'   coords = c("x", "y"), crs = 32632
#' )
#' bnd <- st_sf(geometry = st_sfc(st_polygon(list(rbind(
#'   c(0, 0), c(100, 0), c(100, 100), c(0, 100), c(0, 0)
#' ))), crs = 32632))
#' grid <- create_grid_polygons(bnd, target_cells = 9, type = "square")
#' assigned <- assign_features_to_polygons(pts, grid)
#' table(assigned$poly_id)
#' @family aggregation
#' @seealso [build_tessellation()] to build the polygon layer;
#'   [summarize_by_cell()] for the aggregation step that consumes the result.
#' @export
assign_features_to_polygons <- function(
    features_sf, polygons_sf, polygon_id_col = "poly_id",
    keep_unassigned = FALSE, predicate = sf::st_intersects,
    largest = TRUE,
    tie_break = c("smallest_area", "first")
) {
  if (!inherits(features_sf, "sf"))
    stop("assign_features_to_polygons(): `features_sf` must be an sf object.")
  if (!inherits(polygons_sf, "sf")) {
    if (inherits(polygons_sf, "sfc")) polygons_sf <- sf::st_as_sf(polygons_sf)
    else stop("assign_features_to_polygons(): `polygons_sf` must be sf/sfc.")
  }
  if (nrow(polygons_sf) == 0L)
    stop("assign_features_to_polygons(): `polygons_sf` has zero rows.")
  tie_break <- match.arg(tie_break)

  orig_crs <- sf::st_crs(features_sf)
  hh <- harmonize_crs(features_sf, polygons_sf)
  f <- hh$a; p <- hh$b

  id_candidates <- c(polygon_id_col, "poly_id", "polygon_id", "id", "cell_id", "grid_id")
  id_col <- id_candidates[id_candidates %in% names(p)][1]
  if (length(id_col) == 0L || is.na(id_col)) {
    id_col <- polygon_id_col
    p[[id_col]] <- seq_len(nrow(p))
  }

  p_sel <- p[, id_col, drop = FALSE]

  # st_join() suffixes columns present on both sides (`poly_id.x` /
  # `poly_id.y`), which would defeat the rename below and leave
  # `polygon_id_col` absent -- silently returning zero rows.  This is reachable
  # simply by re-assigning already-assigned points, or points that came out of
  # summarize_by_cell().  Drop the colliding column(s) up front instead.
  collide <- intersect(unique(c(id_col, polygon_id_col)), names(f))
  if (length(collide)) {
    .warn_and_log(
      "assign_features_to_polygons(): `features_sf` already carries column(s) %s, which would collide with the polygon ID column; dropping them before the spatial join. Rename them first if you need to keep them.",
      paste(sprintf("'%s'", collide), collapse = ", ")
    )
    f <- f[, setdiff(names(f), collide), drop = FALSE]
  }

  f$`..pre_join_row_id` <- seq_len(nrow(f))
  
  f_gtypes <- unique(as.character(sf::st_geometry_type(f, by_geometry = TRUE)))
  use_largest <- isTRUE(largest) &&
    all(f_gtypes %in% c("POLYGON", "MULTIPOLYGON"))

  join_args <- list(x = f, y = p_sel, join = predicate, left = TRUE)
  join_ok <- tryCatch({
    if (use_largest) join_args$largest <- TRUE
    do.call(sf::st_join, join_args)
  }, error = function(e) {
    # Fall back without `largest` if the predicate doesn't support it
    join_args$largest <- NULL
    do.call(sf::st_join, join_args)
  })
  joined <- join_ok

  if (!identical(id_col, polygon_id_col)) {
    names(joined)[names(joined) == id_col] <- polygon_id_col
  }
  if (!polygon_id_col %in% names(joined)) {
    stop(sprintf(
      paste0("assign_features_to_polygons(): the spatial join did not produce ",
             "the expected '%s' column (join result has: %s). This usually ",
             "means a column-name collision between `features_sf` and ",
             "`polygons_sf`; rename the offending column in `features_sf`."),
      polygon_id_col, paste(names(joined), collapse = ", ")
    ), call. = FALSE)
  }


  # Deterministic tie-break for features matching multiple polygons.
  # The previous approach (keep first duplicate) was order-dependent.
  dup_mask <- duplicated(joined[["..pre_join_row_id"]])
  if (any(dup_mask) && identical(tie_break, "smallest_area")) {
    # For each duplicated feature, keep the polygon with the smallest area
    # (the most specific / tightest-fitting polygon).
    poly_areas <- suppressWarnings(as.numeric(sf::st_area(p)))
    poly_areas[!is.finite(poly_areas)] <- Inf
    names(poly_areas) <- as.character(p[[id_col]])
    
    joined$`..poly_area` <- poly_areas[as.character(joined[[polygon_id_col]])]
    joined$`..poly_area`[is.na(joined$`..poly_area`)] <- Inf
    
    # Within each group of duplicates, keep the row with smallest area
    joined <- joined[order(joined[["..pre_join_row_id"]], joined[["..poly_area"]]), ,
                     drop = FALSE]
    joined <- joined[!duplicated(joined[["..pre_join_row_id"]]), , drop = FALSE]
    joined[["..poly_area"]] <- NULL
  } else {
    joined <- joined[!duplicated(joined[["..pre_join_row_id"]]), , drop = FALSE]
  }
  joined[["..pre_join_row_id"]] <- NULL

  n_unassigned <- sum(is.na(joined[[polygon_id_col]]))
  # Zero matches is almost never what the caller meant, and it is silent: the
  # function returns an empty sf and every downstream step reports zero rows
  # without saying why.  The usual cause is a CRS mismatch that could only be
  # stamped rather than reprojected, which puts the two layers in different
  # places.
  if (n_unassigned == nrow(joined) && nrow(joined) > 0L)
    .warn_and_log(paste0("assign_features_to_polygons(): none of the %d ",
                         "feature(s) fall inside any polygon. Check that the ",
                         "two layers cover the same ground -- a CRS that had to ",
                         "be stamped rather than reprojected is the usual ",
                         "cause -- or pass a different `predicate`."),
                  nrow(joined))
  else if (n_unassigned > 0L)
    .log_info("assign_features_to_polygons(): %d of %d feature(s) fall inside no polygon.",
              n_unassigned, nrow(joined))

  if (!keep_unassigned) {
    joined <- joined[!is.na(joined[[polygon_id_col]]), , drop = FALSE]
  }

  if (!is.na(orig_crs)) joined <- sf::st_transform(joined, orig_crs)
  joined
}


#' Correlation implied by a fitted gstat variogram model
#'
#' Converts a fitted variogram to a correlation function of distance, giving
#' the correlation between two \strong{distinct} observations a distance
#' \code{h} apart: \code{partial_sill * f(h) / (nugget + partial_sill)}.
#'
#' There is deliberately no special case at \code{h = 0}.  The nugget captures
#' measurement error and variation below the sampling resolution, so two
#' distinct observations each carry independent nugget noise and are less than
#' perfectly correlated even when they coincide.  Self-correlation of 1 belongs
#' on the diagonal of the correlation matrix and is imposed by the caller.
#'
#' @param vgm_model A fitted \code{gstat} variogram model (a data frame with
#'   \code{model}, \code{psill} and \code{range} columns).
#' @return A function of distance returning correlation, or \code{NULL} if the
#'   model cannot be interpreted.
#' @keywords internal
#' @noRd
.vgm_correlation_fn <- function(vgm_model) {
  if (is.null(vgm_model) || !is.data.frame(vgm_model)) return(NULL)
  if (!all(c("model", "psill", "range") %in% names(vgm_model))) return(NULL)

  is_nug <- as.character(vgm_model$model) == "Nug"
  nugget <- sum(vgm_model$psill[is_nug], na.rm = TRUE)
  struct <- vgm_model[!is_nug, , drop = FALSE]
  if (nrow(struct) == 0L) return(NULL)

  psill <- sum(struct$psill, na.rm = TRUE)
  total <- nugget + psill
  if (!is.finite(total) || total <= 0 || !is.finite(psill) || psill <= 0)
    return(NULL)

  rng  <- struct$range[which.max(struct$psill)]
  type <- as.character(struct$model[which.max(struct$psill)])
  if (!is.finite(rng) || rng <= 0) return(NULL)

  ratio <- psill / total
  function(h) {
    f <- switch(
      type,
      Exp = exp(-h / rng),
      Sph = ifelse(h >= rng, 0, 1 - 1.5 * (h / rng) + 0.5 * (h / rng)^3),
      Gau = exp(-(h / rng)^2),
      exp(-h / rng)                     # sensible default for other families
    )
    # No special case at h = 0: this is the correlation between two DISTINCT
    # observations, which the nugget discounts even when they coincide.  An
    # observation with itself correlates at 1, and the caller imposes that on
    # the diagonal.
    pmin(pmax(ratio * f, 0), 1)
  }
}


#' Standard error of a cell mean under a design effect
#'
#' Two corrections, not one.  A design effect inflates the variance of the mean
#' to \eqn{\sigma^2 \mathrm{deff} / n}, which is what \code{s/sqrt(n/deff)}
#' applies --- but under the same within-cell correlation the ordinary sample
#' variance is \emph{also} biased downward.  For exchangeable correlation
#' \eqn{\rho} (Kish's own assumption), with \eqn{\mathrm{deff} = 1 + (n-1)\rho}:
#' \deqn{E[s^2] = \sigma^2 (n - \mathrm{deff}) / (n - 1)}
#' so \eqn{s^2} understates \eqn{\sigma^2} by very nearly the factor deff
#' inflates the mean's variance by, and the two errors compound rather than
#' cancel.  Measured 95\% CI coverage of the uncorrected form at \eqn{n = 30}
#' over 20,000 replicates: 0.921 at \eqn{\rho = 0.2}, 0.844 at \eqn{\rho = 0.5},
#' 0.628 at \eqn{\rho = 0.8}.  With the \eqn{\sqrt{(n-1)/(n-\mathrm{deff})}}
#' rescale below: 0.948, 0.950, 0.949.
#'
#' At \code{deff = 1} the factor is exactly 1, so the default path is
#' bit-identical to \code{s/sqrt(n)}.
#'
#' @param s Within-cell standard deviation.
#' @param n Number of non-missing observations in the cell.
#' @param deff Design effect for that cell, \eqn{\geq 1}.
#' @return The standard error, or \code{NA_real_} when the cell carries no
#'   information about \eqn{\sigma} (\code{deff >= n}: complete redundancy).
#' @keywords internal
#' @noRd
.se_with_deff <- function(s, n, deff, correct_s2 = TRUE) {
  if (!is.finite(s) || !is.finite(n) || n <= 1L) return(NA_real_)
  if (!is.finite(deff) || deff <= 1) return(s / sqrt(n))
  # correct_s2 = FALSE for a design effect the CALLER supplied as a constant.
  # The E[s^2] = sigma^2 (n - deff)/(n - 1) correction below is derived from a
  # within-cell correlation structure, so it is right for an ICC- or
  # variogram-derived deff and unjustified for an externally chosen number:
  # applied to deff = 2 it doubled a 3-point cell's SE instead of multiplying
  # it by sqrt(2), and returned NA for every cell with n <= deff -- neither of
  # which is the "uniform inflation" the argument is documented as.
  if (!isTRUE(correct_s2)) return(s * sqrt(deff / n))
  # deff is bounded above by n (every observation a copy of every other).  At
  # that bound the cell holds one observation's worth of information and s
  # carries none about sigma, so there is no standard error to report.
  if (deff >= n) return(NA_real_)
  s * sqrt(deff / n) * sqrt((n - 1) / (n - deff))
}


#' Per-cell design effect from a fitted variogram
#'
#' For \code{n} observations in a cell with correlation matrix \code{R}, the
#' effective sample size of the mean is \code{n^2 / sum(R)}, so the design
#' effect is \code{sum(R) / n}.  This generalises Kish's
#' \code{1 + (n - 1) * rho}, which is the special case of a constant
#' off-diagonal correlation: substituting \code{R = I + rho(J - I)} gives
#' \code{sum(R) / n = 1 + (n - 1) * rho} exactly.  Using the variogram instead
#' lets correlation decay with distance, which is the whole point of having
#' fitted one -- Kish assumes every pair in a cell is equally correlated
#' regardless of how far apart they are, and that assumption degrades as cells
#' get larger.
#'
#' @param coords Numeric matrix of coordinates.
#' @param cell_id Vector of cell identifiers, same length as \code{nrow(coords)}.
#' @param cor_fn Correlation function from \code{.vgm_correlation_fn()}.
#' @param max_n Cells larger than this are subsampled before forming the
#'   \code{n x n} correlation matrix.  Default 500.
#' @param seed RNG seed for that subsampling.
#' @return Named numeric vector of design effects, one per cell.
#' @keywords internal
#' @noRd
.cell_deff_variogram <- function(coords, cell_id, cor_fn, max_n = 500L,
                                 seed = 42L) {
  r_bar <- .cell_rbar_variogram(coords, cell_id, cor_fn, max_n = max_n,
                                seed = seed)
  n_i   <- table(factor(cell_id, levels = names(r_bar)))
  n_i   <- as.numeric(n_i[names(r_bar)])
  out   <- 1 + (n_i - 1) * r_bar
  # Bounded below by 1 (independence) and above by the cell size (complete
  # redundancy); a NULL cor_fn leaves NA throughout.
  out[!is.na(out)] <- pmin(pmax(out[!is.na(out)], 1), n_i[!is.na(out)])
  out[n_i <= 1L & !is.na(r_bar)] <- 1
  out
}


#' Mean off-diagonal correlation within each cell, from a fitted variogram
#'
#' The quantity a design effect is built from.  For a cell of \eqn{n} points
#' with correlation matrix \eqn{R}, \eqn{\mathrm{deff} = \sum R / n = 1 + (n-1)
#' \bar{r}} where \eqn{\bar{r}} is the mean off-diagonal correlation.  It is
#' \eqn{\bar{r}}, not deff, that is returned, because the design effect a
#' \emph{column} needs is the one for the observations that column actually
#' has: a cell of 10 points with 2 non-missing responses has a response design
#' effect of \eqn{1 + (2-1)\bar{r}}, not \eqn{1 + (10-1)\bar{r}}.  Measured:
#' adding 8 NA-response rows to a 2-observation cell moved its response SE from
#' 8.46 to 26.38 when the cell-level deff was applied to it; NA rows carry no
#' information and must not move anything.
#'
#' Subsampling estimates \eqn{\bar{r}} just as well as the full cell does --
#' the subsample's pairwise-distance distribution is the cell's -- which is why
#' the per-cell quantity kept is this and not \code{sum(R)/n_used}.
#'
#' @param coords Numeric matrix of \strong{projected} coordinates, in the CRS
#'   the variogram was fitted in.  A variogram range is a length in that CRS;
#'   evaluating it at degree distances saturates every pair.
#' @param cell_id Vector of cell identifiers, same length as \code{nrow(coords)}.
#' @param cor_fn Correlation function from \code{.vgm_correlation_fn()}.
#' @param max_n Cells larger than this are subsampled before forming the
#'   correlation matrix.  Default 500.
#' @param seed RNG seed for that subsampling.
#' @return Named numeric vector of mean off-diagonal correlations, one per
#'   cell; \code{NA} when \code{cor_fn} is \code{NULL}, \code{0} for a
#'   single-point cell (no pairs).
#' @keywords internal
#' @noRd
.cell_rbar_variogram <- function(coords, cell_id, cor_fn, max_n = 500L,
                                 seed = 42L) {
  ids <- unique(cell_id)
  ids <- ids[!is.na(ids)]
  out <- stats::setNames(rep(NA_real_, length(ids)), as.character(ids))
  if (is.null(cor_fn)) return(out)
  for (id in ids) {
    idx <- which(cell_id == id)
    out[[as.character(id)]] <- .rbar_from_coords(coords[idx, , drop = FALSE],
                                                 cor_fn, max_n = max_n,
                                                 seed = seed)
  }
  out
}


#' Mean off-diagonal correlation among a set of locations
#'
#' The single-cell core of \code{.cell_rbar_variogram()}, exposed so the
#' \code{..se_} closures can call it on exactly the rows a column is
#' non-missing for.  Subsampling (above \code{max_n}) is seeded and the RNG is
#' restored, so repeated calls from inside \code{dplyr::summarise()} neither
#' drift nor disturb the caller's stream.
#'
#' @return A scalar in \eqn{[0, 1]}; \code{0} when there is at most one point.
#' @keywords internal
#' @noRd
.rbar_from_coords <- function(xy, cor_fn, max_n = 500L, seed = 42L) {
  n_i <- nrow(xy)
  if (is.null(cor_fn) || is.null(n_i) || n_i <= 1L) return(0)
  if (n_i > max_n) {
    cleanup <- .with_seed(seed)
    on.exit(cleanup(), add = TRUE)
    xy <- xy[sample.int(n_i, max_n), , drop = FALSE]
  }
  d <- as.matrix(stats::dist(xy))
  # Rebuild explicitly rather than relying on cor_fn() to preserve `dim`.
  # A correlation function written as, say, rep(1, length(h)) returns a bare
  # vector, and diag()<- would then fail.
  R <- matrix(as.numeric(cor_fn(as.numeric(d))), nrow = nrow(d), ncol = ncol(d))
  diag(R) <- 1
  n_used <- nrow(R)
  r_bar  <- (sum(R) - n_used) / (n_used * (n_used - 1))
  min(max(r_bar, 0), 1)
}


#' Summarize features by polygon/cell ID
#'
#' Aggregates an sf point dataset into one row per cell. By default computes
#' counts and means, but the aggregation function is configurable.
#'
#' This is the third step of the package's pipeline, taking the labelled layer
#' from [assign_features_to_polygons()] down to cell level. What distinguishes
#' it from a plain `dplyr::group_by()` + `summarise()` is that it carries the
#' *uncertainty* of each aggregate with it: alongside every mean it returns a
#' within-cell standard deviation, a standard error and an observation count,
#' and it can correct that standard error for within-cell spatial
#' autocorrelation via `deff`. Reach for it whenever the cell-level values will
#' be modelled or mapped, because a cell mean over 2 observations and one over
#' 200 are not the same measurement and nothing downstream can tell them apart
#' otherwise.
#'
#' In addition to user-specified aggregation functions, this function always
#' computes within-cell standard deviation (`..sd_<var>`) and standard error
#' (`..se_<var>`) for every numeric response/predictor column, plus an `n`
#' column (rows falling in the cell) and a `cell_weight` column.
#' These columns let downstream models account for the fact that a cell with
#' 2 observations carries more aggregation uncertainty than one with 200.
#'
#' `cell_weight` is the *effective* sample size of the primary variable --- the
#' response when one was supplied, otherwise the first predictor. It counts
#' that variable's non-missing rows, not all rows (a cell of 10 rows with 3
#' finite responses carries 3 observations' worth of information about the
#' response, not 10), and it is divided by that cell's design effect when
#' `deff` applied one. With `deff = 1` and no missing values it equals `n`.
#' Pass it as the `weights` argument of a downstream regression.
#'
#' @section Spatial autocorrelation and standard-error bias:
#' **Important:** By default (`deff = 1`), the `..se_*` columns are computed as
#' `sd / sqrt(n)`, which assumes observations within each cell are independent.
#' When data are spatially autocorrelated — the common case for the spatial
#' workflows this package supports — within-cell observations are typically
#' positively correlated, so the effective sample size is smaller than `n`.
#' The naive SE is therefore **anticonservative** (too small), and downstream
#' weighted regressions using `cell_weight` or `..se_*` columns will produce
#' overconfident standard errors for cells with strong intra-cell correlation.
#'
#' Setting `deff = "kish"` applies an approximate correction using Kish's
#' design effect. Separate intra-class correlations (ICCs) are estimated for
#' response and predictor variables via a one-way random-effects decomposition
#' across all cells. Each variable type's ICC is used for its own SE
#' adjustment, and each cell's effective sample size is reduced to
#' `n_i / (1 + (n_i - 1) * rho)`. This is a first-order correction that
#' does not require a full spatial covariance model but does require enough
#' cells and observations for a stable ICC estimate.
#' You may also pass a fixed numeric design effect (e.g. `deff = 2`) to
#' uniformly inflate standard errors: an externally supplied constant is
#' applied as `sd * sqrt(deff / n)`, exactly `sqrt(deff)` times the naive SE in
#' every cell. (The `E[s^2]` correction that the estimated design effects also
#' apply is derived from within-cell correlation and would not be justified for
#' a number the caller chose.)
#'
#' Even with the Kish correction, the adjusted SE is an approximation.
#' For rigorous inference under spatial dependence, consider fitting an
#' explicit spatial covariance model (e.g. via \code{\link{fit_bayesian_spatial_model}}).
#'
#' @section What the standard error estimates:
#' The `..se_*` columns are the standard error of the cell mean **as an
#' estimate of the population (grand) mean** — the unconditional quantity, in
#' which the cell's own realised deviation is part of the error. That is the
#' right quantity when cells are treated as samples from a common population,
#' and the design-effect correction is calibrated for it: measured 95% interval
#' coverage of the grand mean is 0.95 with `deff = "kish"` (and 0.29 with the
#' naive SE) on exchangeable within-cell correlation, and 0.93 with
#' `deff = "variogram"` on a simulated Gaussian field.
#'
#' It is **not** the standard error of the cell's own mean (the block average
#' over that cell), which is what a cell-level map or a regression on cell
#' values usually wants. For that quantity the naive `sd / sqrt(n)` is the
#' better estimate — measured coverage 0.95, against essentially 1.00 for the
#' design-effect-corrected SE, which is about five times too wide. Use `deff`
#' when the cell means feed a population-level inference; leave it at 1 when
#' they are measurements of the cells themselves.
#'
#' @section Design effects and variable types:
#' `deff = "kish"` estimates a separate ICC for response and predictor
#' variables and applies each to its own columns. `deff = "variogram"` fits or
#' accepts **one** correlation function and applies it to every numeric column,
#' because a variogram is a property of the field being modelled rather than of
#' a variable type; a predictor whose spatial structure differs markedly from
#' the response's will have its SE corrected by the response's correlation.
#' Note also that supplying `predictor_vars` makes the internally estimated
#' variogram a **residual** variogram (that is what
#' [estimate_sac_range()] does with predictors), so the response's own design
#' effect depends on which predictors are listed. Pass `sac` explicitly when
#' you want control over which variogram is used.
#'
#' @param assigned_points_sf An sf object with a cell identifier column.
#' @param response_var Optional response column name for per-cell aggregation.
#' @param predictor_vars Optional predictor column names for per-cell aggregation.
#' @param id_col Preferred name of the polygon/cell ID column.
#' @param agg_funs Named list of aggregation functions. Default
#'   \code{list(mean = \(x) mean(x, na.rm = TRUE))}. Additional common options:
#'   \code{median}, \code{sum}, \code{sd}.
#' @param cells_sf Optional polygon sf layer to join cell geometries onto
#'   the output. When supplied, the return value is an sf object with
#'   the polygon geometry from cells_sf, with one row per cell in `cells_sf`
#'   --- cells that no feature fell in are kept, with `NA` summaries. Duplicate
#'   ID values in `cells_sf` would multiply those rows, so they are reported
#'   with a warning. When NULL (default), a plain data.frame/tibble is
#'   returned (previous behaviour).
#' @param deff Design-effect adjustment for standard errors. One of:
#'   \describe{
#'     \item{`1` (default)}{No adjustment; the classic IID standard error
#'       `sd / sqrt(n)`, which assumes the observations within a cell are
#'       independent. No `"deff_applied"` attribute is attached, so a result
#'       carrying none was computed this way.}
#'     \item{`"variogram"`}{Compute a per-cell design effect from a fitted
#'       variogram: for `n` points in a cell with correlation matrix `R`, the
#'       effective sample size of the mean is `n^2 / sum(R)`, so
#'       `deff = sum(R) / n`. This generalises Kish -- substituting a constant
#'       off-diagonal correlation recovers `1 + (n - 1) * rho` exactly -- but
#'       lets correlation decay with distance, which matters increasingly as
#'       cells get larger and Kish's single-`rho` assumption degrades. Supply
#'       the fit via `sac`, or it is estimated when `response_var` is given and
#'       'gstat' is available.}
#'     \item{`"kish"`}{Estimate per-variable-type intra-class correlations
#'       (ICCs) from the grouped data using a one-way random-effects ANOVA
#'       decomposition — one ICC for the response variable and a separate
#'       ICC for the predictor variables — then apply Kish's formula per
#'       cell: `deff_i = 1 + (n_i - 1) * rho`. When multiple columns are
#'       pooled for a single ICC estimate (e.g. several predictor variables),
#'       each column is z-scored before pooling so that variables with
#'       different scales contribute equally to the variance decomposition.
#'       The response-specific ICC
#'       is used for response SEs and the predictor-specific ICC for
#'       predictor SEs. The response's ICC is never applied to predictor
#'       columns or vice versa, and it is the response's ICC (not the
#'       predictors') that sets `cell_weight` whenever a response was given.
#'       Requires at least 2 cells with 2+ observations and at least 2 residual
#'       degrees of freedom (`N - k >= 2`); the ICC is taken as 0 --- no
#'       correction, no `"deff_applied"` attribute --- otherwise, and likewise
#'       when the estimate itself comes out at or below 0.}
#'     \item{A positive number}{Applied as a uniform design effect to every
#'       cell, as `sd * sqrt(deff / n)` --- exactly `sqrt(deff)` times the
#'       naive SE. Use when you have an external estimate of the design
#'       effect. Anything that is not a single number `>= 1` --- including a
#'       value below 1, which would *shrink* the standard errors --- is
#'       refused with a warning and replaced by 1.}
#'   }
#' @param sac Optional `sac_range` object from [estimate_sac_range()], used
#'   when `deff = "variogram"`. Supplying one avoids re-fitting the variogram
#'   and lets you inspect the fit the design effect is based on. A `sac_range`
#'   whose fit was *rejected* (its `status` is not `"ok"`) carries no usable
#'   correlation function, so `deff` falls back to 1 with a warning rather than
#'   correcting by a shape that was not trusted enough to report a range.
#' @param deff_max_n Cells with more than this many points are subsampled
#'   before forming the `n x n` correlation matrix used by
#'   `deff = "variogram"`. Default 500.
#' @param quiet Logical; suppress messages. Default TRUE.
#' @return A tibble/data.frame (or sf if cells_sf given) with per-cell
#'   summaries: the ID column, `n` (rows in the cell --- an input column also
#'   called `n` is not allowed to shadow it), one column per `agg_funs` entry
#'   per variable, `..sd_*` / `..se_*` for every numeric response and predictor,
#'   and `cell_weight`.
#'
#'   When a correction was actually applied, an attribute `"deff_applied"` is
#'   attached recording it: `method` plus `icc_resp`/`icc_pred` for `"kish"`,
#'   `deff`/`rbar`/`crs`/`max_n` for `"variogram"`, and `deff` alone for a
#'   fixed number. When `cells_sf` is supplied, *every* per-cell vector in that
#'   attribute (`deff` and `rbar` alike) is realigned to the joined row order,
#'   so `deff[i]` and `rbar[i]` still describe row `i`; cells with no
#'   observations carry `NA`. No attribute is attached when no correction was
#'   applied --- `deff = 1`, a `deff = "kish"` ICC of 0, or a `"variogram"`
#'   request that could not be fitted.
#'
#'   The ID column keeps its input type when `cells_sf`'s ID column and the
#'   summarised IDs already have the same class. When the classes differ, both
#'   are coerced to character in order to join (logged as a warning), and the
#'   returned ID column is therefore character.
#' @examples
#' library(sf)
#' set.seed(1)
#' n <- 200
#' x <- runif(n, 0, 100)
#' y <- runif(n, 0, 100)
#' # A response with spatial structure, so the within-cell ICC is not zero.
#' pts <- st_as_sf(
#'   data.frame(x = x, y = y, val = 0.05 * x + 0.05 * y + rnorm(n, sd = 0.5)),
#'   coords = c("x", "y"), crs = 32632
#' )
#' bnd <- st_sf(geometry = st_sfc(st_polygon(list(rbind(
#'   c(0, 0), c(100, 0), c(100, 100), c(0, 100), c(0, 0)
#' ))), crs = 32632))
#' grid <- create_grid_polygons(bnd, target_cells = 9, type = "square")
#' assigned <- assign_features_to_polygons(pts, grid)
#'
#' # IID standard errors (default) vs Kish design-effect adjustment
#' naive <- summarize_by_cell(assigned, response_var = "val")
#' kish  <- summarize_by_cell(assigned, response_var = "val", deff = "kish")
#' data.frame(n = naive$n,
#'            se_naive = naive$..se_resp_val,
#'            se_kish  = kish$..se_resp_val)
#' attr(kish, "deff_applied")   # method, icc_resp, icc_pred, per-cell deff
#' @family aggregation
#' @seealso [assign_features_to_polygons()], which produces the input layer;
#'   [build_tessellation()] for the cells themselves.
#' @export
summarize_by_cell <- function(assigned_points_sf,
                              response_var   = NULL,
                              predictor_vars = NULL,
                              id_col         = "poly_id",
                              agg_funs       = list(mean = function(x) mean(x, na.rm = TRUE)),
                              cells_sf       = NULL,
                              deff           = 1,
                              sac            = NULL,
                              deff_max_n     = 500L,
                              quiet          = TRUE) {
  .msg <- function(...) if (!quiet) message(...)
  df <- sf::st_drop_geometry(assigned_points_sf)

  # --- locate ID column ---
  id_candidates <- unique(c(id_col, "poly_id", "polygon_id", "cell_id"))
  id_found <- id_candidates[id_candidates %in% names(df)]
  if (length(id_found) == 0L) {
    stop("summarize_by_cell(): could not find an ID column. Looked for: ",
         paste(id_candidates, collapse = ", "), ".")
  }
  id_col <- id_found[[1]]
  .msg(sprintf("summarize_by_cell(): using id_col = '%s'", id_col))

  # --- validate agg_funs ---
  if (!is.list(agg_funs) || length(agg_funs) == 0L) {
    .log_warn("summarize_by_cell(): invalid agg_funs; falling back to mean.")
    agg_funs <- list(mean = function(x) mean(x, na.rm = TRUE))
  }
  if (is.null(names(agg_funs)) || any(!nzchar(names(agg_funs)))) {
    names(agg_funs) <- paste0("agg", seq_along(agg_funs))
  }

  # --- validate / resolve design effect ---
  use_kish <- identical(deff, "kish")
  use_vgm  <- identical(deff, "variogram")
  if (use_vgm) {
    # The SE closures below do arithmetic on `deff`; the per-cell variogram
    # values are applied after summarising, so neutralise it here.
    deff <- 1
  }
  if (!use_kish && !use_vgm) {
    if (!is.numeric(deff) || length(deff) != 1L || deff < 1) {
      # A real warning, not a log line: silently substituting 1 for a value the
      # caller chose means the standard errors are not the ones they asked for,
      # and nothing in the returned object records that (no `deff_applied` is
      # attached when deff == 1).
      .warn_and_log("summarize_by_cell(): `deff` must be a single number >= 1, \"kish\" or \"variogram\"; got %s. Falling back to deff = 1, so the standard errors are the uncorrected ones.",
                    paste(format(deff), collapse = ", "))
      deff <- 1
    }
  }

  # Helper: estimate intra-class correlation via one-way random-effects ANOVA.
  # Returns a single rho in [0, 1] or 0 if estimation fails.
  .estimate_icc <- function(df, id_col, var_cols) {
    # Pool across all supplied numeric columns for a single global ICC estimate.
    # Each variable is z-scored before pooling so that variables with larger
    # scales do not dominate the between/within variance decomposition.
    #
    # The pooled group is (VARIABLE, cell), not cell alone.  Pooling m columns
    # under the same cell label treats a cell's m z-scores -- one per variable
    # -- as m draws from one distribution, so independent per-variable cell
    # effects average away in the cell mean and the between-cell sum of
    # squares shrinks by about 1/m.  Measured (200 replicates, true rho = 0.5,
    # m = 4): pooled-by-cell ICC 0.121 against per-column ANOVA 0.495; the
    # predictor deff at n = 20 came out 3.3 where 10.5 was correct, and every
    # predictor SE was 44% too small.  Grouping by (variable, cell) is the
    # ordinary pooled one-way ANOVA across the m*k variable-cell groups, and
    # recovers 0.49 on the same data.
    vals <- numeric(0)
    grps <- character(0)
    for (v in var_cols) {
      # A row with no cell id belongs to no group: keep_unassigned = TRUE
      # output carries NA ids, and an NA level makes group_means[grps] NA
      # downstream ("missing value where TRUE/FALSE needed").
      ok <- !is.na(df[[v]]) & !is.na(df[[id_col]])
      raw <- df[[v]][ok]
      s <- sd(raw)
      if (is.na(s) || s < .Machine$double.eps) {
        # Constant column — skip, contributes no variance information
        next
      }
      vals <- c(vals, (raw - mean(raw)) / s)
      grps <- c(grps, paste(v, as.character(df[[id_col]][ok]), sep = "\r"))
    }
    grps <- factor(grps)
    k <- nlevels(grps)
    N <- length(vals)
    if (k < 2L || N < 4L) return(0)

    ni <- tabulate(grps)
    ni <- ni[ni > 0L]
    # The documentation promises "at least 2 cells with 2+ observations; falls
    # back to deff = 1 otherwise", and the k >= 2 / N >= 4 test above does not
    # deliver it.  All-singleton groups give SSW = 0, MSW = 0 and therefore an
    # ICC of exactly 1 -- a design effect of n from data carrying no
    # within-cell information at all.  A single group with two members gives an
    # ICC from one within-group degree of freedom, which is then applied to
    # every cell.  Both are noise dressed as a measurement.
    if (sum(ni >= 2L) < 2L || (N - k) < 2L) return(0)
    grand_mean <- mean(vals)
    group_means <- tapply(vals, grps, mean)

    # Between-group and within-group sum of squares
    SSB <- sum(ni * (group_means - grand_mean)^2)
    SSW <- sum((vals - group_means[grps])^2)

    MSB <- SSB / (k - 1)
    MSW <- if ((N - k) > 0) SSW / (N - k) else 0

    # Unbalanced n0 (Donner 1986)
    n0 <- (N - sum(ni^2) / N) / (k - 1)
    if (n0 < 1) return(0)

    rho <- (MSB - MSW) / (MSB + (n0 - 1) * MSW)
    # Clamp to [0, 1] — negative ICC implies no positive autocorrelation
    max(0, min(1, rho))
  }

  # --- resolve numeric columns for response and predictors ---
  .resolve_numeric <- function(df, cols, label) {
    keep <- cols[cols %in% names(df)]
    if (length(keep) == 0L) return(character(0))
    is_num <- vapply(keep, function(nm) is.numeric(df[[nm]]), logical(1))
    if (!all(is_num)) {
      .msg(sprintf("summarize_by_cell(): non-numeric %s columns skipped: %s",
                   label, paste(keep[!is_num], collapse = ", ")))
    }
    keep[is_num]
  }

  resp_num <- character(0)
  pred_num <- character(0)

  if (!is.null(response_var)) {
    if (response_var %in% names(df)) {
      resp_num <- .resolve_numeric(df, response_var, "response")
    } else {
      .msg(sprintf("summarize_by_cell(): response_var '%s' not found; skipping.",
                   response_var))
    }
  }
  if (!is.null(predictor_vars)) {
    present <- predictor_vars[predictor_vars %in% names(df)]
    if (length(present)) {
      pred_num <- .resolve_numeric(df, present, "predictor")
    } else {
      .msg("summarize_by_cell(): none of the requested predictor_vars are present; skipping.")
    }
  }

  # --- estimate ICC for Kish correction if requested ---
  resp_rho <- NULL
  pred_rho <- NULL
  if (use_kish) {
    if (length(resp_num) > 0L) {
      resp_rho <- .estimate_icc(df, id_col, resp_num)
      if (resp_rho > 0) {
        .msg(sprintf(
          "summarize_by_cell(): Kish deff using response ICC = %.4f", resp_rho
        ))
      } else {
        .msg("summarize_by_cell(): response ICC estimate <= 0; Kish deff not applied for response.")
      }
    }
    if (length(pred_num) > 0L) {
      pred_rho <- .estimate_icc(df, id_col, pred_num)
      if (pred_rho > 0) {
        .msg(sprintf(
          "summarize_by_cell(): Kish deff using predictor ICC = %.4f", pred_rho
        ))
      } else {
        .msg("summarize_by_cell(): predictor ICC estimate <= 0; Kish deff not applied for predictors.")
      }
    }
    if (is.null(resp_rho) && is.null(pred_rho)) {
      .msg("summarize_by_cell(): no numeric columns for ICC estimation; deff = 1.")
      resp_rho <- 0
      pred_rho <- 0
    }
    # Fill in zeros for the type that has no columns
    if (is.null(resp_rho)) resp_rho <- 0
    if (is.null(pred_rho)) pred_rho <- 0
  }

  # --- build combined function list for across() ---
  # Each prefix (resp_, pred_) gets user agg_funs + sd + se, with names
  # constructed so that .names = "{.fn}{.col}" produces the correct columns.
  # e.g. fn name "resp_mean_" + col "temp" -> "resp_mean_temp"
  #      fn name "..sd_resp_" + col "temp" -> "..sd_resp_temp"
  .build_funs <- function(prefix, rho = NULL, rbar = NULL, vgm = NULL) {
    fns <- list()
    for (nm in names(agg_funs)) {
      fn <- agg_funs[[nm]]
      # Force early binding of fn in the closure
      fns[[paste0(prefix, nm, "_")]] <- local({
        f <- fn; function(x) f(x)
      })
    }
    fns[[paste0("..sd_", prefix)]] <- function(x) {
      if (sum(!is.na(x)) > 1L) sd(x, na.rm = TRUE) else NA_real_
    }
    # SE with design-effect adjustment
    fns[[paste0("..se_", prefix)]] <- local({
      .deff <- deff
      .rho  <- rho
      .rbar <- rbar        # per-cell rbar over ALL rows (fast path, no NAs)
      .vgm  <- vgm         # list(coords, cor_fn, max_n): per-column exact path
      .kish <- use_kish
      .id   <- id_col
      function(x) {
        ok      <- !is.na(x)
        n_valid <- sum(ok)
        if (n_valid <= 1L) return(NA_real_)
        s <- sd(x, na.rm = TRUE)
        # Every path builds the design effect from n_valid -- the observations
        # THIS column has in THIS cell -- never from the cell's row count.  An
        # NA row contributes nothing to the mean, so it cannot contribute to
        # the mean's clustering either.
        if (.kish) {
          deff_i <- max(1, 1 + (n_valid - 1) * .rho)
        } else if (!is.null(.rbar)) {
          # Variogram path.  The mean pairwise correlation is over the
          # locations of the observations that ENTER this column's mean, so
          # when the column has NAs in this cell it is recomputed on those rows
          # alone; the precomputed all-rows value is reused when there are
          # none (the ordinary case, and the same number either way).
          if (all(ok)) {
            g  <- as.character(dplyr::cur_group()[[.id]])
            rb <- if (g %in% names(.rbar)) .rbar[[g]] else NA_real_
          } else {
            rows <- dplyr::cur_group_rows()[ok]
            rb   <- .rbar_from_coords(.vgm$coords[rows, , drop = FALSE],
                                      .vgm$cor_fn, max_n = .vgm$max_n)
          }
          deff_i <- if (is.finite(rb)) min(max(1, 1 + (n_valid - 1) * rb), n_valid)
                    else 1
        } else {
          deff_i <- .deff
        }
        # NOT s / sqrt(n_valid / deff_i): that corrects the mean's variance for
        # clustering but leaves s^2 biased low by the same clustering.  See
        # .se_with_deff() for the derivation and the measured coverage.  The
        # E[s^2] correction applies only to a design effect this function
        # DERIVED from the data's own clustering; a numeric `deff` supplied by
        # the caller is applied as the uniform inflation it is documented to be.
        .se_with_deff(s, n_valid, deff_i,
                      correct_s2 = .kish || !is.null(.rbar))
      }
    })
    fns
  }

  # --- per-cell design effect from the variogram ------------------------------
  # Computed separately from the aggregation because the across() closures see
  # only a column of values, never the geometry.  SE scales as sqrt(deff), so
  # the SEs are formed with deff = 1 below and rescaled per cell afterwards.
  vgm_deff <- NULL
  vgm_rbar <- NULL
  vgm_bits <- NULL
  if (use_vgm) {
    # A REJECTED sac is not a fit.  estimate_sac_range() returns NA with a
    # `rejected_reason` when the fitted range runs past the longest lag the
    # variogram covers, or when the optimiser stopped at its iteration limit --
    # and it says in as many words that such a range "is where the optimiser
    # halted rather than a fitted parameter" and must not be used.  It still
    # attaches `variogram_model` so the fit can be INSPECTED, and reading that
    # attribute without checking the value took the refused model as gospel:
    # the correlation function then saturated at every within-cell distance,
    # deff came out equal to n, cell_weight collapsed to 1 and standard errors
    # were inflated 40-50x.  Treat it like no fit at all.
    if (!is.null(sac) && is.na(suppressWarnings(as.numeric(sac))[1L]) &&
        !is.null(attr(sac, "rejected_reason"))) {
      .warn_and_log(paste0("summarize_by_cell(): the supplied `sac` reports no ",
                           "usable range (%s), so its fitted variogram model ",
                           "cannot size a design effect. Falling back to ",
                           "deff = 1."),
                    as.character(attr(sac, "rejected_reason")))
      sac <- NULL
    }
    vgm_model <- attr(sac, "variogram_model")
    vgm_crs   <- attr(sac, "crs")
    if (is.null(vgm_model)) {
      if (!is.null(response_var) && requireNamespace("gstat", quietly = TRUE)) {
        .msg("summarize_by_cell(): no fitted variogram supplied; estimating one.")
        est <- try(estimate_sac_range(assigned_points_sf, response_var,
                                      predictor_vars = predictor_vars),
                   silent = TRUE)
        # Same test on the internally estimated fit: a rejected range means the
        # model behind it is not usable either.
        if (!inherits(est, "try-error") &&
            !(is.na(suppressWarnings(as.numeric(est))[1L]) &&
              !is.null(attr(est, "rejected_reason")))) {
          vgm_model <- attr(est, "variogram_model")
          vgm_crs   <- attr(est, "crs")
        } else if (!inherits(est, "try-error")) {
          .warn_and_log(paste0("summarize_by_cell(): the estimated variogram ",
                               "reports no usable range (%s), so it cannot size ",
                               "a design effect. Falling back to deff = 1."),
                        as.character(attr(est, "rejected_reason")))
        }
      }
    }
    cor_fn <- .vgm_correlation_fn(vgm_model)
    if (is.null(cor_fn)) {
      .log_warn(paste0("summarize_by_cell(): deff = \"variogram\" requires a ",
                       "fitted variogram model; none was available (pass one ",
                       "via `sac = estimate_sac_range(...)`). Falling back to ",
                       "deff = 1."))
      use_vgm <- FALSE
      deff <- 1
    } else {
      # st_coordinates() returns one row per VERTEX, so any non-POINT geometry
      # (this function accepts POLYGON and MULTIPOINT features) would misalign
      # coords_mat with `df` and feed the wrong points into every cell.
      # coerce_to_points() takes representative points, matching the guard in
      # estimate_sac_range().
      pts_for_deff <- assigned_points_sf
      if (!all(sf::st_geometry_type(pts_for_deff, by_geometry = TRUE) == "POINT"))
        pts_for_deff <- coerce_to_points(pts_for_deff, "auto")

      # The variogram range is a LENGTH in the CRS the variogram was fitted in
      # -- estimate_sac_range() always fits in a projected CRS, in metres.
      # Evaluating the correlation function at distances taken from
      # unprojected input (degrees) saturates every pair: measured on the same
      # points, cells and variogram, median SE 0.657 with UTM input against
      # 235.9 with EPSG:4326 input, a factor of 354.  Put the points in the
      # variogram's own CRS when it is recorded, else project them the way
      # estimate_sac_range() would have.
      vgm_crs_ok <- !is.null(vgm_crs) && inherits(vgm_crs, "crs") && !is.na(vgm_crs)
      pts_for_deff <- if (vgm_crs_ok) {
        .transform_or_stamp(pts_for_deff, vgm_crs, what = "assigned_points_sf",
                            caller = "summarize_by_cell")
      } else {
        ensure_projected(pts_for_deff)
      }
      coords_mat <- sf::st_coordinates(pts_for_deff)[, 1:2, drop = FALSE]

      # Per-cell mean off-diagonal correlation, NOT a per-cell deff: the design
      # effect each column needs is 1 + (n_valid - 1) * rbar for ITS non-missing
      # count, and that is formed inside the ..se_ closures.
      vgm_rbar <- .cell_rbar_variogram(coords_mat, df[[id_col]], cor_fn,
                                       max_n = deff_max_n)
      vgm_deff <- .cell_deff_variogram(coords_mat, df[[id_col]], cor_fn,
                                       max_n = deff_max_n)
      vgm_bits <- list(coords = coords_mat, cor_fn = cor_fn, max_n = deff_max_n)
      .msg(sprintf(
        "summarize_by_cell(): variogram deff across cells (at the cell row count): median %.3f, max %.3f",
        stats::median(vgm_deff, na.rm = TRUE), max(vgm_deff, na.rm = TRUE)
      ))
    }
  }

  # --- single grouped summarise for all columns ---
  grouped <- df |> dplyr::group_by(.data[[id_col]])

  has_resp <- length(resp_num) > 0L
  has_pred <- length(pred_num) > 0L

  if (has_resp && has_pred) {
    out <- grouped |>
      dplyr::summarise(
        ..n_rows = dplyr::n(),
        dplyr::across(dplyr::all_of(resp_num), .build_funs("resp_", rho = resp_rho, rbar = vgm_rbar, vgm = vgm_bits),
                      .names = "{.fn}{.col}"),
        dplyr::across(dplyr::all_of(pred_num), .build_funs("pred_", rho = pred_rho, rbar = vgm_rbar, vgm = vgm_bits),
                      .names = "{.fn}{.col}"),
        .groups = "drop"
      )
  } else if (has_resp) {
    out <- grouped |>
      dplyr::summarise(
        ..n_rows = dplyr::n(),
        dplyr::across(dplyr::all_of(resp_num), .build_funs("resp_", rho = resp_rho, rbar = vgm_rbar, vgm = vgm_bits),
                      .names = "{.fn}{.col}"),
        .groups = "drop"
      )
  } else if (has_pred) {
    out <- grouped |>
      dplyr::summarise(
        ..n_rows = dplyr::n(),
        dplyr::across(dplyr::all_of(pred_num), .build_funs("pred_", rho = pred_rho, rbar = vgm_rbar, vgm = vgm_bits),
                      .names = "{.fn}{.col}"),
        .groups = "drop"
      )
  } else {
    out <- grouped |>
      dplyr::summarise(..n_rows = dplyr::n(), .groups = "drop")
  }
  # dplyr::summarise() makes each newly created column visible to the
  # expressions that follow it, so naming the row count `n` here put it in
  # scope for the across() calls below: a response or predictor column
  # literally named `n` was silently summarised as the ROW COUNT (resp_mean_n
  # equal to n, sd and se NA), with no message.  Build it under a reserved
  # name and rename once nothing can shadow it.
  names(out)[names(out) == "..n_rows"] <- "n"

  # cell_weight is the EFFECTIVE sample size of the primary variable -- the
  # response when there is one, else the first predictor -- so it is built from
  # that variable's non-missing count per cell, not from the row count: a cell
  # of 10 rows with 3 finite responses and deff = 2 has 1.5 observations' worth
  # of information about the response, not 5.  (`n` still reports rows.)
  primary_col <- if (has_resp) resp_num[[1L]] else if (has_pred) pred_num[[1L]] else NULL
  n_valid_primary <- if (is.null(primary_col)) out$n else {
    cnt <- tapply(!is.na(df[[primary_col]]), df[[id_col]], sum)
    v   <- as.numeric(cnt[as.character(out[[id_col]])])
    v[is.na(v)] <- 0
    v
  }
  out$cell_weight <- n_valid_primary

  # Record the design effect used and adjust cell_weight for effective n
  # cell_weight uses the response ICC (primary outcome); falls back to
  # predictor ICC when no response variable was supplied.
  if (use_kish) {
    # The fallback to the PREDICTOR's ICC is for the case the comment above
    # describes -- no response variable at all.  Applying it whenever the
    # response's own ICC came out <= 0 meant that adding a predictor to the
    # call changed the response's regression weight (measured: by a factor of
    # 20), deflating a cell whose response is uncorrelated on the strength of a
    # predictor that happens to be clustered.  The response SE is already left
    # at deff = 1 in that case, as documented; cell_weight has to agree with it.
    primary_rho <- if (has_resp) max(resp_rho, 0)
                   else if (has_pred) max(pred_rho, 0)
                   else 0
    if (primary_rho > 0) {
      deff_per_cell <- pmax(1, 1 + (n_valid_primary - 1) * primary_rho)
      out$cell_weight <- n_valid_primary / deff_per_cell
      attr(out, "deff_applied") <- list(
        method   = "kish",
        icc_resp = if (has_resp) resp_rho else NA_real_,
        icc_pred = if (has_pred) pred_rho else NA_real_,
        deff     = deff_per_cell
      )
    }
  } else if (use_vgm && !is.null(vgm_rbar)) {
    # The ..se_ closures have already applied the per-column design effect
    # (see .build_funs()); nothing is rescaled here.  What is recorded is the
    # design effect of the PRIMARY variable at its own non-missing count, which
    # is also what cell_weight is built from.
    # Same rule as the closures: rbar over the primary column's observed rows.
    rb <- vapply(seq_len(nrow(out)), function(i) {
      cell <- out[[id_col]][i]
      rows <- which(df[[id_col]] == cell & !is.na(df[[primary_col]]))
      if (length(rows) == n_valid_primary[i] && length(rows) == sum(df[[id_col]] == cell, na.rm = TRUE)) {
        v <- vgm_rbar[as.character(cell)]
        if (is.finite(v)) unname(v) else 0
      } else {
        .rbar_from_coords(coords_mat[rows, , drop = FALSE], cor_fn, max_n = deff_max_n)
      }
    }, numeric(1))
    d_i <- pmin(pmax(1, 1 + (n_valid_primary - 1) * rb), pmax(n_valid_primary, 1))

    out$cell_weight <- n_valid_primary / d_i
    attr(out, "deff_applied") <- list(
      method = "variogram",
      deff   = d_i,
      rbar   = rb,
      crs    = if (isTRUE(vgm_crs_ok)) vgm_crs else sf::st_crs(pts_for_deff),
      max_n  = deff_max_n
    )
  } else if (!use_kish && is.numeric(deff) && deff > 1) {
    out$cell_weight <- n_valid_primary / deff
    attr(out, "deff_applied") <- list(method = "fixed", deff = deff)
  }
  
  if (!is.null(cells_sf)) {
    if (!inherits(cells_sf, "sf")) {
      .log_warn("summarize_by_cell(): cells_sf is not an sf object; returning plain data.frame.")
    } else {
      # Locate the matching ID column in cells_sf
      cells_id_candidates <- unique(c(id_col, "poly_id", "polygon_id", "cell_id"))
      cells_id_found <- cells_id_candidates[cells_id_candidates %in% names(cells_sf)]
      if (length(cells_id_found) > 0L) {
        cells_id <- cells_id_found[[1]]
        # Keep only geometry + id from cells to avoid column collisions
        cells_slim <- cells_sf[, cells_id, drop = FALSE]
        # A duplicated cell ID makes left_join() emit one summary row per
        # matching cell, so sum(n) exceeded the number of points assigned (39
        # for 30 points; 46 for 23 after st_cast(., "POLYGON") split a
        # multi-part cell and each piece kept its poly_id).  Nothing checked
        # it, and the inflated total is easy to mistake for a real count.
        if (anyDuplicated(cells_slim[[cells_id]]))
          .warn_and_log(paste0("summarize_by_cell(): `cells_sf$%s` has %d ",
                               "duplicated value(s), so each summary row is ",
                               "repeated once per matching cell and totals such ",
                               "as sum(n) will exceed the number of points. Give ",
                               "every cell its own ID (see ",
                               "ensure_stable_poly_id())."),
                        cells_id, sum(duplicated(cells_slim[[cells_id]])))
        if (cells_id != id_col) {
          names(cells_slim)[names(cells_slim) == cells_id] <- id_col
        }
        # Join on the native type when both sides already agree.  Coercing
        # unconditionally made the returned ID type depend on an unrelated
        # argument and turned integer IDs into "1", "10", "2", ...
        if (!identical(class(out[[id_col]]), class(cells_slim[[id_col]]))) {
          .log_warn(
            "summarize_by_cell(): cells_sf$%s is <%s> but the summarised IDs are <%s>; coercing both to character to join, so the returned '%s' column is character.",
            id_col, paste(class(cells_slim[[id_col]]), collapse = "/"),
            paste(class(out[[id_col]]), collapse = "/"), id_col
          )
          out[[id_col]] <- as.character(out[[id_col]])
          cells_slim[[id_col]] <- as.character(cells_slim[[id_col]])
        }

        # dplyr::left_join() rebuilds attributes from the `x` template, which
        # silently drops "deff_applied".  Save it, then re-attach it, mapping
        # the per-cell vector onto the joined row order (cells with no
        # observations get NA).
        deff_attr   <- attr(out, "deff_applied")
        pre_join_id <- as.character(out[[id_col]])

        out <- dplyr::left_join(cells_slim, out, by = id_col)

        if (!is.null(deff_attr)) {
          # EVERY per-cell vector has to follow the join, not just $deff.
          # Realigning $deff alone left $rbar in pre-join order and pre-join
          # length, so deff[i] and rbar[i] described different cells and the
          # identity deff = 1 + (n-1) * rbar no longer held row-wise.
          for (fld in c("deff", "rbar")) {
            v <- deff_attr[[fld]]
            if (!is.null(v) && length(v) == length(pre_join_id)) {
              lookup <- stats::setNames(v, pre_join_id)
              deff_attr[[fld]] <- unname(lookup[as.character(out[[id_col]])])
            }
          }
          attr(out, "deff_applied") <- deff_attr
        }
      } else {
        .log_warn("summarize_by_cell(): cells_sf has no matching ID column; returning plain data.frame.")
      }
    }
  }

  out
}
