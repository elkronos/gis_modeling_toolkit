#' Assign features to polygons and attach a polygon ID
#'
#' Joins an sf layer of input features to a polygon layer via spatial join.
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
#' @return An sf object with polygon_id_col attached.
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


#' Summarize features by polygon/cell ID
#'
#' Aggregates an sf point dataset into one row per cell. By default computes
#' counts and means, but the aggregation function is configurable.
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
#' explicit spatial covariance model (e.g. via \code{\link{fit_bayesian_model}}).
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
#' @param quiet Logical; suppress messages. Default TRUE.
#' @return A tibble/data.frame (or sf if cells_sf given) with per-cell summaries
#'   including `n`, `cell_weight`, and `..sd_*` / `..se_*` columns.
#'   When `deff != 1`, an attribute `"deff_applied"` is attached to the result
#'   recording the design effect(s) used.
#' @export
summarize_by_cell <- function(assigned_points_sf,
                              response_var   = NULL,
                              predictor_vars = NULL,
                              id_col         = "poly_id",
                              agg_funs       = list(mean = function(x) mean(x, na.rm = TRUE)),
                              cells_sf       = NULL,
                              deff           = 1,
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
  if (!use_kish) {
    if (!is.numeric(deff) || length(deff) != 1L || deff < 1) {
      .log_warn("summarize_by_cell(): deff must be >= 1 or \"kish\"; falling back to 1.")
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
      .msg(sprintf("summarize_by_cell(): non-numeric columns skipped: %s",
                   paste(keep[!is_num], collapse = ", ")))
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
        s / sqrt(n_valid / deff_i)
      }
    })
    fns
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
        out[[id_col]] <- as.character(out[[id_col]])
        cells_slim[[id_col]] <- as.character(cells_slim[[id_col]])
        out <- dplyr::left_join(cells_slim, out, by = id_col)
      } else {
        .log_warn("summarize_by_cell(): cells_sf has no matching ID column; returning plain data.frame.")
      }
    }
  }

  out
}
