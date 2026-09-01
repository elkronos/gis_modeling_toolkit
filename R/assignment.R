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
#' @param keep_unassigned Logical; retain unmatched features. Default FALSE.
#' @param predicate Binary spatial predicate function. Default sf::st_intersects.
#' @param largest Logical; for polygon-on-polygon joins with overlapping
#'   polygons, keep the polygon with the largest overlap. Default TRUE. Only
#'   effective with predicates that support it (e.g. st_intersects).
#' @param tie_break Strategy for resolving features that match multiple
#'   polygons: \code{"smallest_area"} (default) keeps the polygon with the
#'   smallest area, \code{"first"} keeps the first match (original order-dependent
#'   behavior).
#' @return An sf object with polygon_id_col attached. Any column of
#'   `features_sf` whose name would collide with the polygon ID column is
#'   dropped before the spatial join (with a warning), so re-assigning an
#'   already-assigned layer replaces the old IDs rather than failing.
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
    .log_warn(
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
.se_with_deff <- function(s, n, deff) {
  if (!is.finite(s) || !is.finite(n) || n <= 1L) return(NA_real_)
  if (!is.finite(deff) || deff <= 1) return(s / sqrt(n))
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
  ids <- unique(cell_id)
  out <- stats::setNames(rep(NA_real_, length(ids)), as.character(ids))
  if (is.null(cor_fn)) return(out)

  cleanup <- .with_seed(seed)
  on.exit(cleanup(), add = TRUE)

  for (id in ids) {
    idx <- which(cell_id == id)
    n_i <- length(idx)
    if (n_i <= 1L) { out[[as.character(id)]] <- 1; next }
    if (n_i > max_n) idx <- sample(idx, max_n)

    d <- as.matrix(stats::dist(coords[idx, , drop = FALSE]))
    # Rebuild explicitly rather than relying on cor_fn() to preserve `dim`.
    # A correlation function written as, say, rep(1, length(h)) returns a bare
    # vector, and diag()<- would then fail.
    R <- matrix(as.numeric(cor_fn(as.numeric(d))), nrow = nrow(d), ncol = ncol(d))
    diag(R) <- 1
    n_used <- nrow(R)
    # deff = sum(R)/n_i for the WHOLE cell, which is 1 + (n_i - 1) * Rbar with
    # Rbar the mean off-diagonal correlation.  Subsampling estimates Rbar just
    # as well (the subsample's pairwise-distance distribution is the cell's),
    # but sum(R)/n_used answers for a cell of size n_used -- so a subsampled
    # cell used to be reported at the design effect of `max_n` points rather
    # than of its own n_i.  Measured: 4000 points, exponential correlation with
    # a 60-unit range, true deff 1821.8; sum(R)/n_used at max_n = 500 gave
    # 228.6, the rescale below gives 1825.2.
    #
    # When nothing was subsampled, n_used == n_i and this is algebraically
    # identical to the old sum(R)/n_i, so unsubsampled cells are unchanged.
    r_bar <- (sum(R) - n_used) / (n_used * (n_used - 1))
    out[[as.character(id)]] <- min(max(1 + (n_i - 1) * r_bar, 1), n_i)
  }
  out
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
#' (`..se_<var>`) for every numeric response/predictor column, plus a
#' `cell_weight` column equal to the observation count.
#' These columns let downstream models account for the fact that a cell with
#' 2 observations carries more aggregation uncertainty than one with 200.
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
#' uniformly inflate standard errors.
#'
#' Even with the Kish correction, the adjusted SE is an approximation.
#' For rigorous inference under spatial dependence, consider fitting an
#' explicit spatial covariance model (e.g. via \code{\link{fit_bayesian_spatial_model}}).
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
#'   the polygon geometry from cells_sf. When NULL (default), a plain
#'   data.frame/tibble is returned (previous behaviour).
#' @param deff Design-effect adjustment for standard errors. One of:
#'   \describe{
#'     \item{`1` (default)}{No adjustment; classic IID standard error.
#'       Equivalent to previous behaviour but now emits a message (when
#'       `quiet = FALSE`) reminding that SEs assume independence.}
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
#'       predictor SEs. Requires at least 2 cells with 2+ observations;
#'       falls back to `deff = 1` otherwise.}
#'     \item{A positive number}{Applied as a uniform design effect to every
#'       cell. Use when you have an external estimate of the design effect.}
#'   }
#' @param sac Optional `sac_range` object from [estimate_sac_range()], used
#'   when `deff = "variogram"`. Supplying one avoids re-fitting the variogram
#'   and lets you inspect the fit the design effect is based on.
#' @param deff_max_n Cells with more than this many points are subsampled
#'   before forming the `n x n` correlation matrix used by
#'   `deff = "variogram"`. Default 500.
#' @param quiet Logical; suppress messages. Default TRUE.
#' @return A tibble/data.frame (or sf if cells_sf given) with per-cell summaries
#'   including `n`, `cell_weight`, and `..sd_*` / `..se_*` columns.
#'   When `deff != 1`, an attribute `"deff_applied"` is attached to the result
#'   recording the design effect(s) used — including when `cells_sf` is
#'   supplied, in which case its per-cell `deff` vector is realigned to the
#'   joined row order and cells with no observations carry `NA`. The one
#'   exception is `deff = "kish"` with an estimated ICC of 0: no correction is
#'   applied, so no attribute is attached.
#'
#'   The ID column keeps its input type when `cells_sf`'s ID column and the
#'   summarised IDs already have the same class. When the classes differ, both
#'   are coerced to character in order to join (logged as a warning), and the
#'   returned ID column is therefore character.
#' @examples
#' library(sf)
#' set.seed(1)
#' pts <- st_as_sf(
#'   data.frame(x = runif(40, 0, 100), y = runif(40, 0, 100), val = rnorm(40)),
#'   coords = c("x", "y"), crs = 32632
#' )
#' bnd <- st_sf(geometry = st_sfc(st_polygon(list(rbind(
#'   c(0, 0), c(100, 0), c(100, 100), c(0, 100), c(0, 0)
#' ))), crs = 32632))
#' grid <- create_grid_polygons(bnd, target_cells = 9, type = "square")
#' assigned <- assign_features_to_polygons(pts, grid)
#'
#' # IID standard errors (default) vs Kish design-effect adjustment
#' cells <- summarize_by_cell(assigned, response_var = "val", deff = "kish")
#' cells
#' attr(cells, "deff_applied")
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
      .log_warn("summarize_by_cell(): deff must be >= 1, \"kish\" or \"variogram\"; falling back to 1.")
      deff <- 1
    }
  }

  # Helper: estimate intra-class correlation via one-way random-effects ANOVA.
  # Returns a single rho in [0, 1] or 0 if estimation fails.
  .estimate_icc <- function(df, id_col, var_cols) {
    # Pool across all supplied numeric columns for a single global ICC estimate.
    # Each variable is z-scored before pooling so that variables with larger
    # scales do not dominate the between/within variance decomposition.
    vals <- numeric(0)
    grps <- character(0)
    for (v in var_cols) {
      ok <- !is.na(df[[v]])
      raw <- df[[v]][ok]
      s <- sd(raw)
      if (is.na(s) || s < .Machine$double.eps) {
        # Constant column — skip, contributes no variance information
        next
      }
      vals <- c(vals, (raw - mean(raw)) / s)
      grps <- c(grps, as.character(df[[id_col]][ok]))
    }
    grps <- factor(grps)
    k <- nlevels(grps)
    N <- length(vals)
    if (k < 2L || N < 4L) return(0)

    ni <- tabulate(grps)
    ni <- ni[ni > 0L]
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
  .build_funs <- function(prefix, rho = NULL) {
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
      .kish <- use_kish
      function(x) {
        n_valid <- sum(!is.na(x))
        if (n_valid <= 1L) return(NA_real_)
        s <- sd(x, na.rm = TRUE)
        if (.kish) {
          deff_i <- max(1, 1 + (n_valid - 1) * .rho)
        } else {
          deff_i <- .deff
        }
        # NOT s / sqrt(n_valid / deff_i): that corrects the mean's variance for
        # clustering but leaves s^2 biased low by the same clustering.  See
        # .se_with_deff() for the derivation and the measured coverage.
        .se_with_deff(s, n_valid, deff_i)
      }
    })
    fns
  }

  # --- per-cell design effect from the variogram ------------------------------
  # Computed separately from the aggregation because the across() closures see
  # only a column of values, never the geometry.  SE scales as sqrt(deff), so
  # the SEs are formed with deff = 1 below and rescaled per cell afterwards.
  vgm_deff <- NULL
  if (use_vgm) {
    vgm_model <- attr(sac, "variogram_model")
    if (is.null(vgm_model)) {
      if (!is.null(response_var) && requireNamespace("gstat", quietly = TRUE)) {
        .msg("summarize_by_cell(): no fitted variogram supplied; estimating one.")
        est <- try(estimate_sac_range(assigned_points_sf, response_var,
                                      predictor_vars = predictor_vars),
                   silent = TRUE)
        if (!inherits(est, "try-error")) vgm_model <- attr(est, "variogram_model")
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
      coords_mat <- sf::st_coordinates(pts_for_deff)[, 1:2, drop = FALSE]
      vgm_deff <- .cell_deff_variogram(coords_mat, df[[id_col]], cor_fn,
                                       max_n = deff_max_n)
      .msg(sprintf(
        "summarize_by_cell(): variogram deff across cells: median %.3f, max %.3f",
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
        n = dplyr::n(),
        dplyr::across(dplyr::all_of(resp_num), .build_funs("resp_", rho = resp_rho),
                      .names = "{.fn}{.col}"),
        dplyr::across(dplyr::all_of(pred_num), .build_funs("pred_", rho = pred_rho),
                      .names = "{.fn}{.col}"),
        .groups = "drop"
      )
  } else if (has_resp) {
    out <- grouped |>
      dplyr::summarise(
        n = dplyr::n(),
        dplyr::across(dplyr::all_of(resp_num), .build_funs("resp_", rho = resp_rho),
                      .names = "{.fn}{.col}"),
        .groups = "drop"
      )
  } else if (has_pred) {
    out <- grouped |>
      dplyr::summarise(
        n = dplyr::n(),
        dplyr::across(dplyr::all_of(pred_num), .build_funs("pred_", rho = pred_rho),
                      .names = "{.fn}{.col}"),
        .groups = "drop"
      )
  } else {
    out <- grouped |>
      dplyr::summarise(n = dplyr::n(), .groups = "drop")
  }
  out$cell_weight <- out$n

  # Record the design effect used and adjust cell_weight for effective n
  # cell_weight uses the response ICC (primary outcome); falls back to
  # predictor ICC when no response variable was supplied.
  if (use_kish) {
    primary_rho <- if (has_resp && resp_rho > 0) resp_rho
                   else if (has_pred && pred_rho > 0) pred_rho
                   else 0
    if (primary_rho > 0) {
      deff_per_cell <- pmax(1, 1 + (out$n - 1) * primary_rho)
      out$cell_weight <- out$n / deff_per_cell
      attr(out, "deff_applied") <- list(
        method   = "kish",
        icc_resp = if (has_resp) resp_rho else NA_real_,
        icc_pred = if (has_pred) pred_rho else NA_real_,
        deff     = deff_per_cell
      )
    }
  } else if (use_vgm && !is.null(vgm_deff)) {
    # Match each cell's deff by id, then rescale.  The ..se_ closures ran with
    # deff = 1, so they hold s / sqrt(n); the factor below is exactly what
    # .se_with_deff() would have returned had the per-cell deff been available
    # inside the closure -- sqrt(deff) for the mean's variance AND
    # sqrt((n - 1)/(n - deff)) for the downward bias in s^2 itself.  Rescaling
    # by sqrt(deff) alone would reproduce the under-coverage documented there.
    d_i <- unname(vgm_deff[as.character(out[[id_col]])])
    d_i[!is.finite(d_i) | d_i < 1] <- 1

    n_i   <- out$n
    infl  <- sqrt(d_i) * sqrt((n_i - 1) / (n_i - d_i))
    # deff >= n is complete redundancy: no information about sigma is left.
    infl[!is.finite(infl) | d_i >= n_i] <- NA_real_
    infl[d_i <= 1] <- 1

    se_cols <- grep("^\\.\\.se_", names(out), value = TRUE)
    for (cn in se_cols) out[[cn]] <- out[[cn]] * infl

    out$cell_weight <- out$n / d_i
    attr(out, "deff_applied") <- list(
      method = "variogram",
      deff   = d_i,
      max_n  = deff_max_n
    )
  } else if (!use_kish && is.numeric(deff) && deff > 1) {
    out$cell_weight <- out$n / deff
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
          if (length(deff_attr$deff) == length(pre_join_id)) {
            lookup <- stats::setNames(deff_attr$deff, pre_join_id)
            deff_attr$deff <- unname(lookup[as.character(out[[id_col]])])
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
