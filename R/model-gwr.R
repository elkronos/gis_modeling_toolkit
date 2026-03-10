# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

#' Validate a GWR kernel name
#'
#' GWmodel accepts kernel names as character strings directly (unlike spgwr
#' which required function objects).
#'
#' @param kernel Character scalar.
#' @return The validated kernel string.
#' @keywords internal
#' @noRd
.validate_kernel <- function(kernel) {
  valid <- c("bisquare", "gaussian", "tricube", "boxcar", "exponential")
  kernel <- tolower(kernel)
  if (!kernel %in% valid) {
    .log_warn("fit_gwr_model(): unknown kernel '%s'; falling back to 'bisquare'.", kernel)
    kernel <- "bisquare"
  }
  kernel
}


#' Coerce an sf object to SpatialPointsDataFrame for GWmodel
#'
#' GWmodel's core functions (gwr.basic, bw.gwr, gwr.predict) currently
#' require Spatial* inputs.  Unlike the archived spgwr package, GWmodel is
#' actively maintained and may gain native sf support in the future;
#' centralizing the coercion here makes a future migration trivial.
#'
#' @param data_sf An sf object with POINT geometry.
#' @param keep_cols Character vector of column names to retain. NULL = all.
#' @return A SpatialPointsDataFrame.
#' @keywords internal
#' @noRd
.to_sp <- function(data_sf, keep_cols = NULL) {
  if (!requireNamespace("sp", quietly = TRUE))
    stop(".to_sp(): package 'sp' is required for GWmodel.", call. = FALSE)
  if (!is.null(keep_cols)) {
    keep_cols <- intersect(keep_cols, names(data_sf))
    data_sf <- data_sf[, keep_cols, drop = FALSE]
  }
  methods::as(data_sf, "Spatial")
}


#' Compute a sensible fallback bandwidth from spatial data
#'
#' For adaptive mode, returns an integer (number of nearest neighbours).
#' For fixed mode, returns a distance in CRS units (~1/3 of diagonal extent).
#'
#' @param sp_dat A Spatial* object with coordinates.
#' @param adaptive Logical; whether adaptive bandwidth is used.
#' @return Numeric scalar.
#' @keywords internal
#' @noRd
.fallback_bandwidth <- function(sp_dat, adaptive) {
  if (adaptive) {
    n <- nrow(sp_dat@coords)
    # Aim for ~50 nearest neighbours, clamped to [10, 0.9*n]
    return(as.integer(min(max(10L, 50L), floor(0.9 * n))))
  }
  bb <- sp::bbox(sp_dat)
  dx <- bb[1, 2] - bb[1, 1]
  dy <- bb[2, 2] - bb[2, 1]
  diag <- sqrt(dx^2 + dy^2)
  max(diag / 3, .Machine$double.eps)
}


#' Extract fitted or predicted values from a GWmodel GWR result
#'
#' **Unified extraction function** used by both in-sample evaluation
#' (10_evaluation.R) and cross-validation prediction (09_cross_validation.R).
#' Replaces the previously duplicated .extract_gwr_fitted() and
#' .extract_gwr_predictions() functions.
#'
#' Implements four strategies in order:
#'   1. Look for a direct prediction/fitted column in the SDF.
#'   2. Reconstruct from local coefficients × design matrix.
#'   3. Compute y − residual (works for in-sample gwr.basic results).
#'   4. Check if the response variable was placed in the SDF.
#'
#' @param gwr_obj The GWR result object (from gwr.basic or gwr.predict).
#'   Must have a \code{$SDF} component.
#' @param data_sf The sf data used to fit or predict (for design matrix
#'   reconstruction and observed values).
#' @param formula The regression formula.
#' @param n Expected number of observations/predictions.
#' @param response_var Name of the response variable. Defaults to the LHS
#'   of \code{formula}.
#' @return Numeric vector of fitted/predicted values (NA where extraction failed).
#' @keywords internal
#' @noRd
.extract_gwr_values <- function(gwr_obj, data_sf, formula, n,
                                response_var = all.vars(formula)[1]) {
  na_vec <- rep(NA_real_, n)
  
  sdf <- tryCatch(gwr_obj$SDF, error = function(e) NULL)
  if (is.null(sdf)) return(na_vec)
  
  sdf_data <- tryCatch({
    if (inherits(sdf, "sf")) {
      sf::st_drop_geometry(sdf)
    } else if (methods::is(sdf, "Spatial")) {
      sdf@data
    } else {
      as.data.frame(sdf)
    }
  }, error = function(e) NULL)
  if (is.null(sdf_data) || nrow(sdf_data) != n) return(na_vec)
  
  # Strategy 1: direct prediction/fitted column.
  # GWmodel::gwr.basic stores fitted values in "yhat";
  # gwr.predict may use "prediction" or the response name.
  pred_col_names <- c("yhat", "pred", "prediction", "fitted", "fit")
  hit <- names(sdf_data)[tolower(names(sdf_data)) %in% pred_col_names]
  if (length(hit) >= 1L) {
    vals <- suppressWarnings(as.numeric(sdf_data[[hit[1]]]))
    if (length(vals) == n && any(is.finite(vals))) return(vals)
  }
  
  # Strategy 2: reconstruct from local coefficients × design matrix
  mm <- tryCatch(
    stats::model.matrix(formula, data = sf::st_drop_geometry(data_sf)),
    error = function(e) NULL
  )
  if (!is.null(mm) && nrow(mm) == n) {
    # Normalize names so that e.g. "(Intercept)" and "Intercept" both
    # become "intercept", avoiding silent drops of the intercept term.
    .norm_names <- function(x) tolower(gsub("[^[:alnum:]_.]", "", x))
    cn_mm  <- .norm_names(colnames(mm))
    cn_sdf <- .norm_names(names(sdf_data))
    shared <- intersect(cn_mm, cn_sdf)
    if (length(shared) >= 1L) {
      mm_sub <- mm[, match(shared, cn_mm), drop = FALSE]
      cf_sub <- as.matrix(sdf_data[, match(shared, cn_sdf), drop = FALSE])
      if (nrow(cf_sub) == n) {
        y_hat <- rowSums(mm_sub * cf_sub)
        if (any(is.finite(y_hat))) return(y_hat)
      }
    }
  }
  
  # Strategy 3: y − residual (works for in-sample / gwr.basic)
  resid_cols <- c("residual", "gwr.e", "resid")
  resid_hit <- names(sdf_data)[tolower(names(sdf_data)) %in% resid_cols]
  if (length(resid_hit) >= 1L &&
      response_var %in% names(sf::st_drop_geometry(data_sf))) {
    y_obs <- sf::st_drop_geometry(data_sf)[[response_var]]
    resid <- suppressWarnings(as.numeric(sdf_data[[resid_hit[1]]]))
    if (length(resid) == n && length(y_obs) == n) return(y_obs - resid)
  }
  
  # Strategy 4: response variable placed in the SDF directly
  if (response_var %in% names(sdf_data)) {
    vals <- suppressWarnings(as.numeric(sdf_data[[response_var]]))
    if (length(vals) == n && any(is.finite(vals))) return(vals)
  }
  
  .log_warn(".extract_gwr_values(): all extraction strategies failed.")
  na_vec
}


# -----------------------------------------------------------------------------
# Main GWR fitting function
# -----------------------------------------------------------------------------

#' Fit a Geographically Weighted Regression (GWR) via GWmodel
#'
#' Fits a GWR using GWmodel on an sf dataset with either adaptive or fixed
#' bandwidth.
#'
#' @param data_sf An sf object with response, predictors, and geometries.
#' @param response_var Response column name.
#' @param predictor_vars Predictor column names.
#' @param adaptive Logical; use adaptive bandwidth. Default TRUE. When TRUE,
#'   bandwidth is an integer number of nearest neighbours. When FALSE,
#'   bandwidth is a fixed distance in CRS units.
#' @param bandwidth Optional numeric bandwidth value. For adaptive mode this
#'   is an integer (number of neighbours); for fixed mode a distance in CRS
#'   units. If NULL (default), bandwidth is selected automatically via
#'   \code{GWmodel::bw.gwr()}.
#' @param kernel Kernel function type. One of "bisquare" (default),
#'   "gaussian", "tricube", "boxcar", "exponential".
#' @param .already_prepped Logical (internal). If \code{TRUE}, skip the
#'   \code{prep_model_data()} call because the caller has already projected,
#'   coerced, and filtered the data.  Used by the CV internals to avoid a
#'   redundant second pass on every fold.  End users should leave this at the
#'   default \code{FALSE}.
#'
#' @section Collinearity diagnostics:
#' The function checks the condition number of the predictor matrix and warns
#' when it exceeds a threshold.
#' A **global** condition number is computed on the full predictor matrix.
#' In addition, a **local** spot-check is performed at a small random sample
#' of locations: for each sampled point, the nearest neighbours within the
#' bandwidth window are selected and the condition number of that local
#' (weighted) design sub-matrix is evaluated.  If the fraction of sampled
#' locations with an extreme local condition number (> 1e6) exceeds 25\%, a
#' separate warning is issued.
#'
#' Because the local spot-check examines only a subset of locations (up to 30
#' by default), it may not detect every problematic neighbourhood.  Users
#' working with highly clustered data or near-collinear predictors should
#' consider a full local-collinearity audit as a post-fit diagnostic.
#'
#' @return A \code{gwr_fit} object (inherits from \code{spatial_fit}).
#'   Supports \code{predict()}, \code{fitted()}, \code{residuals()},
#'   \code{coef()}, \code{summary()}, and \code{model_metrics()}.
#'   Model-specific metadata lives in \code{$info} (bandwidth, adaptive,
#'   kernel, AICc).  The raw GWmodel result is in \code{$engine}.
#' @export
fit_gwr_model <- function(data_sf, response_var, predictor_vars,
                          adaptive = TRUE,
                          bandwidth = NULL,
                          kernel = c("bisquare", "gaussian", "tricube",
                                     "boxcar", "exponential"),
                          .already_prepped = FALSE) {
  if (!inherits(data_sf, "sf"))
    stop("fit_gwr_model(): `data_sf` must be an sf object.")
  if (!requireNamespace("GWmodel", quietly = TRUE))
    stop("fit_gwr_model(): package 'GWmodel' is required. Install with install.packages('GWmodel').",
         call. = FALSE)
  if (!requireNamespace("sp", quietly = TRUE))
    stop("fit_gwr_model(): package 'sp' is required (for GWmodel interop).",
         call. = FALSE)
  kernel <- match.arg(kernel)
  kernel <- .validate_kernel(kernel)
  
  # prep_model_data() handles: point coercion, CRS projection, NA/non-finite
  # row removal — no need to duplicate that logic here.
  # When called from CV internals the data is already prepped; skip the
  
  # redundant pass to avoid re-projecting, re-coercing geometry, and
  # re-scanning for NAs on every fold.
  if (isTRUE(.already_prepped)) {
    dat <- data_sf
  } else {
    dat <- prep_model_data(
      data_sf = data_sf, response_var = response_var,
      predictor_vars = predictor_vars, pointize = "auto"
    )
  }
  
  # Warn or error if response looks non-continuous.
  # Gaussian GWR assumes a continuous response; binary data should error.
  resp_vals <- sf::st_drop_geometry(dat)[[response_var]]
  if (is.numeric(resp_vals)) {
    is_integer_like <- all(resp_vals == round(resp_vals), na.rm = TRUE)
    n_unique <- length(unique(resp_vals[is.finite(resp_vals)]))
    if (is_integer_like && n_unique <= 2L) {
      stop(
        sprintf("fit_gwr_model(): response '%s' is binary (%d unique values). Gaussian GWR is invalid for binary outcomes. Consider GWmodel::ggwr.basic() with family = 'binomial'.",
                response_var, n_unique),
        call. = FALSE
      )
    } else if (is_integer_like && n_unique <= 10L) {
      warning(
        sprintf("fit_gwr_model(): response '%s' is integer-valued with only %d unique values. Gaussian GWR assumes a continuous response; results may be unreliable for counts or ordinal outcomes.",
                response_var, n_unique),
        call. = FALSE
      )
    } else if (is_integer_like && n_unique <= 30L) {
      .log_warn(
        "fit_gwr_model(): response '%s' appears integer-valued (%d unique values). Verify that a Gaussian GWR is appropriate.",
        response_var, n_unique
      )
    }
  }
  
  n_obs <- nrow(dat)
  n_params <- length(predictor_vars) + 1L  # +1 for intercept
  
  pred_df <- sf::st_drop_geometry(dat)[, predictor_vars, drop = FALSE]
  for (pv in predictor_vars) {
    if (is.numeric(pred_df[[pv]]) && stats::sd(pred_df[[pv]], na.rm = TRUE) < .Machine$double.eps * 100) {
      .log_warn("fit_gwr_model(): predictor '%s' has near-zero variance; GWR may be unstable.", pv)
    }
  }
  if (length(predictor_vars) >= 2L) {
    num_preds <- predictor_vars[vapply(pred_df[predictor_vars], is.numeric, logical(1))]
    if (length(num_preds) >= 2L) {
      xmat <- as.matrix(pred_df[, num_preds, drop = FALSE])
      cn <- tryCatch(kappa(xmat, exact = FALSE), error = function(e) Inf)
      if (is.finite(cn) && cn > 1e6) {
        .log_warn(
          "fit_gwr_model(): global predictor matrix condition number = %.0f (collinearity risk). Note: local collinearity within bandwidth windows may be substantially worse than this global value.",
          cn
        )
      }
      # --- Local collinearity spot-check ---
      # Sample a few locations and check the condition number of the local
      # (nearest-neighbour) design sub-matrix.  This catches cases where the
      # global condition number looks benign but spatially clustered subsets
      # have near-zero predictor variance.
      coords <- sf::st_coordinates(dat)
      n_spot <- min(30L, n_obs)
      spot_idx <- if (n_obs <= 30L) seq_len(n_obs) else sample.int(n_obs, n_spot)
      n_local_extreme <- 0L
      for (si in spot_idx) {
        dists <- sqrt((coords[, 1] - coords[si, 1])^2 +
                        (coords[, 2] - coords[si, 2])^2)
        # Use the same neighbour count that will be used for fitting;
        # for adaptive, take bw-nearest; for fixed, take points within bw.
        if (adaptive) {
          local_bw <- if (!is.null(bandwidth)) as.integer(bandwidth) else min(50L, n_obs)
          nn_idx <- order(dists)[seq_len(min(local_bw, n_obs))]
        } else {
          local_bw_dist <- if (!is.null(bandwidth)) as.numeric(bandwidth) else Inf
          nn_idx <- which(dists <= local_bw_dist)
          if (length(nn_idx) < length(num_preds) + 1L) nn_idx <- order(dists)[seq_len(length(num_preds) + 1L)]
        }
        local_xmat <- xmat[nn_idx, , drop = FALSE]
        local_cn <- tryCatch(kappa(local_xmat, exact = FALSE), error = function(e) Inf)
        if (is.finite(local_cn) && local_cn > 1e6) n_local_extreme <- n_local_extreme + 1L
      }
      frac_extreme <- n_local_extreme / n_spot
      if (frac_extreme > 0.25) {
        .log_warn(
          "fit_gwr_model(): local collinearity spot-check: %.0f%% of %d sampled locations have condition number > 1e6. Local regressions may be unstable within bandwidth windows.",
          frac_extreme * 100, n_spot
        )
      } else if (n_local_extreme > 0L) {
        .log_warn(
          "fit_gwr_model(): local collinearity spot-check: %d of %d sampled locations have condition number > 1e6.",
          n_local_extreme, n_spot
        )
      }
    }
  }
  
  if (n_obs < n_params * 3L) {
    .log_warn("fit_gwr_model(): only %d observations for %d parameters; local regressions will be underdetermined.",
              n_obs, n_params)
  }
  
  # NOTE: dat is already projected by prep_model_data(), so no
  # additional .is_longlat() / ensure_projected() call is needed.
  
  fml <- stats::reformulate(termlabels = predictor_vars, response = response_var)
  
  needed_cols <- unique(c(response_var, predictor_vars))
  sp_dat <- .to_sp(dat, needed_cols)
  
  # --- Bandwidth selection ---
  bandwidth_is_fallback <- FALSE
  if (is.null(bandwidth)) {
    bw <- tryCatch(
      suppressWarnings(
        GWmodel::bw.gwr(fml, data = sp_dat, approach = "AICc",
                        kernel = kernel, adaptive = adaptive)
      ),
      error = function(e) {
        .log_warn("fit_gwr_model(): bw.gwr() failed: %s", conditionMessage(e))
        NA_real_
      }
    )
    
    if (!is.finite(bw) || is.na(bw) || bw <= 0) {
      bw <- .fallback_bandwidth(sp_dat, adaptive)
      bandwidth_is_fallback <- TRUE
      warning(
        sprintf(
          "fit_gwr_model(): automatic bandwidth selection failed; using arbitrary fallback bandwidth = %.4f. This fallback has no relationship to the data's spatial structure and may produce a poor fit. Consider supplying an explicit `bandwidth` argument.",
          bw
        ),
        call. = FALSE
      )
    }
  } else {
    bw <- as.numeric(bandwidth)
  }
  
  # Clamp bandwidth to safe range
  if (adaptive) {
    bw <- as.integer(round(bw))
    min_bw <- n_params + 1L
    max_bw <- n_obs
    if (bw < min_bw) {
      .log_warn("fit_gwr_model(): adaptive bandwidth %d too small for %d params; clamping to %d.",
                bw, n_params, min_bw)
      bw <- min_bw
    }
    if (bw > max_bw) bw <- max_bw
  }
  
  # --- Fit GWR ---
  fit <- tryCatch(
    GWmodel::gwr.basic(formula = fml, data = sp_dat, bw = bw,
                       kernel = kernel, adaptive = adaptive),
    error = function(e) {
      stop(sprintf("fit_gwr_model(): GWR fit failed: %s", conditionMessage(e)),
           call. = FALSE)
    }
  )
  
  # AICc extraction
  AICc_val <- NA_real_
  if (!is.null(fit$GW.diagnostic) && !is.null(fit$GW.diagnostic$AICc)) {
    AICc_val <- suppressWarnings(as.numeric(fit$GW.diagnostic$AICc))
  }
  
  new_spatial_fit(
    subclass       = "gwr_fit",
    engine         = fit,
    formula        = fml,
    response_var   = response_var,
    predictor_vars = predictor_vars,
    data_sf        = dat,
    info           = list(
      bandwidth             = as.numeric(bw),
      adaptive              = adaptive,
      kernel                = kernel,
      AICc                  = AICc_val,
      bandwidth_is_fallback = bandwidth_is_fallback
    )
  )
}