# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

#' Coerce an sf object to SpatialPointsDataFrame
#'
#' A fit's engine keeps its results as a Spatial object, the layout GWmodel's
#' functions return, and the GWR functions take the training rows and their
#' two coordinates from one.  Centralising the coercion here keeps that in one
#' place.
#'
#' @param data_sf An sf object with POINT geometry.
#' @param keep_cols Character vector of column names to retain. NULL = all.
#' @return A SpatialPointsDataFrame with two coordinate columns (any Z or M
#'   is dropped).
#' @keywords internal
#' @noRd
.to_sp <- function(data_sf, keep_cols = NULL) {
  if (!requireNamespace("sp", quietly = TRUE))
    stop(".to_sp(): package 'sp' is required: a GWR fit keeps its data and ",
         "its local estimates as Spatial objects.", call. = FALSE)
  if (!is.null(keep_cols)) {
    keep_cols <- intersect(keep_cols, names(data_sf))
    data_sf <- data_sf[, keep_cols, drop = FALSE]
  }
  # A GWR is strictly 2-D: distances are taken between the first two
  # coordinates.  prep_model_data() already drops Z/M; this covers the callers
  # that skip it (fit_gwr_model(.already_prepped = TRUE)) and a fit made
  # before it did.
  data_sf <- .drop_zm(data_sf)
  methods::as(data_sf, "Spatial")
}


#' Is the AICc outside its domain?
#'
#' The AICc of a GWR as GWmodel defines it (Hurvich et al. 1998, its C++
#' \code{gwr_diag1()} and \code{AICc_rss1()}; \code{.gwr_diagnostics()}
#' computes the same) is
#' \deqn{n \log(RSS/n) + n \log(2\pi) + n (n + tr S) / (n - 2 - tr S),}
#' defined only for \eqn{tr S < n - 2}.  Past that point the penalty changes
#' sign, and a (near-)interpolating fit scores a huge NEGATIVE AICc that ranks
#' it first everywhere it is compared.  GWmodel's AIC from the same fit is
#' \eqn{n \log(RSS/n) + n \log(2\pi) + n + tr S}, so
#' \eqn{AICc - AIC = (n + tr S)(2 + tr S) / (n - 2 - tr S)}: AICc exceeds AIC
#' exactly when the denominator is positive.  The domain can therefore be
#' read off the two numbers alone, which is what a table of them (a selection
#' sweep's, or one an injected engine returns) allows; \code{fit_gwr_model()}
#' has the trace itself and checks both.  An
#' infinite AICc (a residual sum of squares of exactly 0, or
#' \eqn{tr S = n - 2}) is outside the domain too.
#'
#' @param aic,aicc Numeric vectors of equal length, the AIC and the AICc.
#' @return Logical vector: \code{TRUE} where AICc is known to be undefined,
#'   \code{FALSE} where it is defined or cannot be checked (\code{NA} AIC).
#' @keywords internal
#' @noRd
.gwr_aicc_undefined <- function(aic, aicc) {
  aic  <- suppressWarnings(as.numeric(aic))
  aicc <- suppressWarnings(as.numeric(aicc))
  is.infinite(aicc) | (is.finite(aic) & is.finite(aicc) & aicc <= aic)
}


#' The hat-matrix trace, recovered from the AIC
#'
#' \code{AIC - (n log(RSS/n) + n log(2 pi) + n)}; see
#' \code{.gwr_aicc_undefined()}.  For a fit that carries the AIC and the
#' residual sum of squares but not the trace: it loses precision when
#' the residual sum of squares is tiny, and is \code{NA} when it is zero.
#' @keywords internal
#' @noRd
.gwr_trace_s <- function(aic, rss, n) {
  tr <- suppressWarnings(as.numeric(aic) -
                           (n * log(as.numeric(rss) / n) + n * log(2 * pi) + n))
  if (length(tr) == 1L && is.finite(tr)) tr else NA_real_
}


#' Smallest adaptive bandwidth that leaves every local fit a residual df
#'
#' GWmodel's bisquare and tricube kernels give the k-th nearest neighbour,
#' which sets the adaptive distance, a weight of exactly zero, so an adaptive
#' bandwidth of k fits each window on k - 1 points.  A floor of
#' \code{n_params + 1} therefore made every window an exact interpolator for
#' those two kernels (tr S = n, R2 = 1, AICc near -n^2/2).  The boxcar keeps
#' the k-th neighbour and the Gaussian and exponential kernels weight every
#' point, so \code{n_params + 1} leaves them one residual degree of freedom.
#' No floor keeps tr S below n - 2 at small n; \code{.gwr_aicc_undefined()}
#' is the guard for that.
#'
#' @param n_params Parameters in the local model, the intercept included.
#' @param kernel The kernel name.
#' @return Integer.
#' @keywords internal
#' @noRd
.gwr_min_adaptive_bw <- function(n_params, kernel) {
  as.integer(n_params) + if (kernel %in% c("bisquare", "tricube")) 2L else 1L
}


#' Warn when a fixed GWR bandwidth is tiny against the data's extent
#'
#' A fixed bandwidth is a length in the CRS the fit runs in, and that is not
#' necessarily the CRS the caller was looking at: prep_model_data() projects
#' geographic input, so 0.2 "degrees" becomes 0.2 metres and every local
#' regression holds its own observation alone.  Shared by
#' \code{fit_gwr_model()} and \code{gwr_model_selection()}.
#'
#' @param dat The prepared sf layer the fit runs on.
#' @param bandwidth The supplied fixed bandwidth.
#' @param fn Name of the calling function, for the message.
#' @return Invisibly \code{NULL}; called for the warning.
#' @keywords internal
#' @noRd
.gwr_warn_tiny_fixed_bw <- function(dat, bandwidth, fn) {
  fit_units <- tryCatch(sf::st_crs(dat)$units_gdal, error = function(e) NULL)
  extent_x  <- suppressWarnings(tryCatch({
    bb <- sf::st_bbox(dat); as.numeric(bb[["xmax"]] - bb[["xmin"]])
  }, error = function(e) NA_real_))
  if (length(extent_x) == 1L && is.finite(extent_x) && extent_x > 0 &&
      as.numeric(bandwidth) < extent_x / 1e4) {
    .warn_and_log(paste0("%s(): a fixed bandwidth of %s is less ",
                         "than a ten-thousandth of the data's extent (%s %s ",
                         "across). The bandwidth is a distance in the CRS the ",
                         "fit runs in (%s), which prep_model_data() may have ",
                         "chosen -- geographic input is projected first, so a ",
                         "value in degrees is read as metres. Every local ",
                         "window is then likely to hold its own observation ",
                         "alone."),
                  fn, format(as.numeric(bandwidth)), format(signif(extent_x, 4)),
                  if (is.null(fit_units) || is.na(fit_units)) "unit"
                  else fit_units,
                  tryCatch(sf::st_crs(dat)$input, error = function(e) "unknown"))
  }
  invisible(NULL)
}


#' What makes a local window rank deficient, and what to do about it
#'
#' The closing sentences of the warning \code{fit_gwr_model()} raises when
#' the local solve finds windows of deficient rank.
#'
#' @param selected The bandwidth, formatted, when it was selected by the
#'   search and not supplied; \code{NULL} otherwise.
#' @return A character string starting with a space.
#' @keywords internal
#' @noRd
.gwr_rank_hint <- function(selected = NULL) {
  paste0(" Typical causes are a predictor that is constant inside a window ",
         "(a regional indicator, say), a bandwidth that leaves a window ",
         "fewer observations than parameters (a fixed bandwidth in the wrong ",
         "units, for one), and neighbours tied at the kernel's edge (a ",
         "regular grid), which the bisquare and tricube kernels give weight ",
         "0. Use a larger bandwidth or drop the predictor; a predictor that ",
         "is constant over whole regions can only be given one global ",
         "coefficient (a mixed GWR, GWmodel::gwr.mixed()).",
         if (is.null(selected)) "" else sprintf(paste0(
           " No bandwidth was supplied: the search selected this one (%s) ",
           "on the fitted values, which those windows do estimate; supply a ",
           "larger one as `bandwidth` if every coefficient is wanted at ",
           "every location."), selected))
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
    # Aim for ~50 nearest neighbours, clamped to [10, 0.9*n].
    # (The previous expression min(max(10L, 50L), ...) collapsed the lower
    # clamp to the constant 50; for very small n the result could fall
    # below 10 and n never entered the lower bound.)
    return(as.integer(max(10L, min(50L, floor(0.9 * n)))))
  }
  # Guard locally rather than relying on the caller, matching .to_sp().  Both
  # current callers check first, so this only changes the message a future
  # caller would see.
  if (!requireNamespace("sp", quietly = TRUE))
    stop(".fallback_bandwidth(): package 'sp' is required for the fixed-",
         "bandwidth extent calculation.", call. = FALSE)
  bb <- sp::bbox(sp_dat)
  dx <- bb[1, 2] - bb[1, 1]
  dy <- bb[2, 2] - bb[2, 1]
  diag <- sqrt(dx^2 + dy^2)
  max(diag / 3, .Machine$double.eps)
}


#' Extract fitted or predicted values from a GWR engine object
#'
#' Used by \code{fitted.gwr_fit()} (in-sample evaluation,
#' R/model-classes.R).  \code{mode = "predict"} reads the table of a
#' \code{GWmodel::gwr.predict()} result, which \code{predict.gwr_fit()} does
#' not produce (see \code{.gwr_predict_at()}).
#'
#' A fit made by \code{fit_gwr_model()} carries its fitted values in
#' \code{$fitted}, and those are returned.  For any other engine object
#' (GWmodel's own, or one from an earlier version) four strategies follow, in
#' order:
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
                                response_var = all.vars(formula)[1],
                                mode = c("auto", "basic", "predict")) {
  mode <- match.arg(mode)
  na_vec <- rep(NA_real_, n)

  # A fit made by fit_gwr_model() carries its fitted values outright.  They
  # are not always the coefficients times the design: in a window of
  # deficient rank a coefficient is NA while the fitted value is estimable.
  # Matched exactly: `$fitted` would also take an element called
  # `fitted.values`.
  fv <- tryCatch(gwr_obj[["fitted", exact = TRUE]], error = function(e) NULL)
  if (mode != "predict" && is.numeric(fv) && length(fv) == n)
    return(as.numeric(fv))

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
  #
  # Two properties of this search are load-bearing.
  #
  # (a) It runs in PREFERENCE order, not in the SDF's column order.
  #     gwr.basic's SDF is c(colnames(betas), "y", "yhat", "residual", ...),
  #     so the local COEFFICIENTS come first.  Matching the whole name vector
  #     with %in% and taking hit[1] therefore picked whichever candidate name
  #     the DATA happened to list earliest: with a predictor called `fit`,
  #     `pred`, `prediction`, `fitted` or `yhat`, that is the coefficient
  #     surface for that predictor, not the fitted values.
  #
  # (b) Any candidate that is also a model term is excluded outright.  In
  #     gwr.basic's SDF such a column is certainly a coefficient surface; in
  #     gwr.predict's it is suffixed "_coef" and would not match anyway.
  #
  # The failure this prevents is silent: fitted() returned a coefficient
  # column, the in-sample R2 went negative, and residuals(), summary(),
  # model_metrics(), compare_models() and every cv_gwr() fold consumed it
  # without a warning.
  pred_col_names <- c("yhat", "pred", "prediction", "fitted", "fit")
  sdf_lower   <- tolower(names(sdf_data))
  model_terms <- tolower(c(response_var, all.vars(formula)[-1L],
                           "Intercept", "(Intercept)"))
  # The model-term exclusion is for gwr.basic()'s SDF, where the local
  # coefficient columns carry the BARE predictor names.  gwr.predict()'s SDF
  # suffixes them "_coef" (Intercept_coef, x_coef, ...) and names its
  # prediction column "prediction" -- so there, a predictor that happens to
  # be called `prediction` is NOT a reason to skip the column literally named
  # prediction; it IS the prediction, and excluding it made predict() return
  # all NA (0 finite of 100) while fitted() was fine.  Apply the exclusion only
  # when no "_coef" columns exist.
  # `mode` says which GWmodel call produced this SDF, because inferring it from
  # the column names cannot be made safe: a PREDICTOR whose name happens to end
  # in "_coef" made a gwr.basic() SDF look like a gwr.predict() one, which
  # switched off the model-term exclusion -- and with a second predictor named
  # `yhat`, fitted() then returned that predictor's local coefficient surface
  # instead of the fitted values (R2 -2.28 against 0.986), silently, for
  # residuals(), summary() and model_metrics() alike.  "auto" keeps the old
  # inference for any caller that has not been updated.
  has_coef_suffix <- switch(mode,
    basic   = FALSE,
    predict = TRUE,
    any(grepl("_coef$", sdf_lower))
  )
  hit <- NULL
  for (cand in pred_col_names) {
    idx <- which(sdf_lower == cand &
                   (has_coef_suffix | !(sdf_lower %in% model_terms)))
    if (length(idx) >= 1L) { hit <- idx[[1L]]; break }
  }
  if (!is.null(hit)) {
    # Index positionally: an SDF can carry two columns of the same name (a
    # predictor named "yhat" alongside gwr.basic's own), and [[<name>]] would
    # silently take the first.
    vals <- suppressWarnings(as.numeric(sdf_data[[hit]]))
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
    # EVERY model-matrix column must be matched.  With a partial match this
    # reconstructs a linear predictor missing one or more terms and returns it
    # as the fitted value: plausible-looking numbers that are simply wrong,
    # feeding fitted(), residuals(), summary() and every metric with no
    # warning.  A GWmodel column rename or a differently-spelled factor
    # contrast is enough to trigger it, so fall through to strategies 3 and 4
    # instead of guessing.
    if (length(shared) == ncol(mm)) {
      mm_sub <- mm[, match(shared, cn_mm), drop = FALSE]
      cf_df  <- sdf_data[, match(shared, cn_sdf), drop = FALSE]
      # A non-numeric coefficient column cannot be used either way: as.matrix()
      # would coerce the WHOLE selection to character and make `mm_sub *
      # cf_sub` throw an uncaught "non-numeric argument to binary operator",
      # while data.matrix() would quietly substitute factor codes for it.
      # Refuse and fall through instead.
      if (!all(vapply(cf_df, is.numeric, logical(1)))) {
        .log_warn(paste0(".extract_gwr_values(): the fit's SDF carries a ",
                         "non-numeric column for a model term, so the local ",
                         "coefficients cannot be multiplied through; trying ",
                         "the remaining strategies."))
      } else {
        cf_sub <- data.matrix(cf_df)
        if (nrow(cf_sub) == n) {
          y_hat <- rowSums(mm_sub * cf_sub)
          if (any(is.finite(y_hat))) return(y_hat)
        }
      }
    } else if (length(shared) >= 1L) {
      .log_warn(paste0(".extract_gwr_values(): only %d of %d model-matrix ",
                       "term(s) matched a column in the fit's SDF, so the ",
                       "local coefficients cannot reconstruct the full linear ",
                       "predictor; trying the remaining strategies."),
                length(shared), ncol(mm))
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

#' Fit a Geographically Weighted Regression (GWR)
#'
#' Fits a Gaussian GWR on an sf dataset with an adaptive or a fixed
#' bandwidth.  The kernels, the bandwidth conventions and the diagnostics are
#' those of the GWmodel package.  The local regressions are solved by this
#' package, each as a weighted least squares by pivoted QR, so a window whose
#' design is rank deficient neither stops the fit nor returns coefficients
#' made of rounding error.
#'
#' @param data_sf An sf object with response, predictors, and geometries.
#' @param response_var Response column name.
#' @param predictor_vars Predictor column names (numeric columns; a name given
#'   twice counts once).  \code{character(0)} fits the model with a local
#'   intercept alone: the fitted value at a location is the kernel-weighted
#'   mean of the response around it.  That is the baseline a GWR with
#'   predictors has to beat, since any GWR beats a global mean through its
#'   local intercept, and \code{\link{select_features_forward}()} fits it as
#'   the null model of a GWR sweep.
#' @param adaptive Logical; use adaptive bandwidth. Default TRUE. When TRUE,
#'   bandwidth is an integer number of nearest neighbours. When FALSE,
#'   bandwidth is a fixed distance in CRS units.
#' @param bandwidth Optional numeric bandwidth value. For adaptive mode this
#'   is an integer (number of neighbours); for fixed mode a distance in the
#'   units of the **projected** CRS the fit runs in. \code{prep_model_data()}
#'   projects geographic input before the bandwidth is used, so 0.2 supplied
#'   for lon/lat data is 0.2 metres, not 0.2 degrees. Read the CRS off
#'   \code{sf::st_crs(fit$data_sf)}, or pass \code{target_crs} to
#'   \code{\link{ensure_projected}} beforehand to fix the units yourself.
#'   If NULL (default), the bandwidth with the lowest AICc is searched for
#'   (the lowest leave-one-location-out error when observations share a
#'   location); see \strong{Bandwidth selection}.  With \code{adaptive = FALSE}, a
#'   bandwidth smaller than a ten-thousandth of the data's extent raises a
#'   warning naming the extent and the CRS the fit runs in: every local
#'   window then holds its own observation alone.  With \code{adaptive = TRUE}
#'   the count is rounded, and one too small for the model is raised, with a
#'   warning, to the smallest that gives every local regression more points
#'   of non-zero weight than parameters: the number of predictors plus 3 for
#'   the bisquare and tricube kernels (which give the farthest neighbour in a
#'   window weight 0), plus 2 for the others.  That floor is enough unless
#'   several neighbours tie at the kernel's edge (a regular grid), which
#'   leaves a window fewer weighted points; then use a larger bandwidth.
#'   A count above the number of
#'   observations is capped at it, with a warning (a distance meant for
#'   \code{adaptive = FALSE}, most often).
#' @param kernel Kernel function type. One of "bisquare" (default),
#'   "gaussian", "tricube", "boxcar", "exponential".
#' @param .already_prepped Logical (internal). If \code{TRUE}, skip the
#'   \code{prep_model_data()} call because the caller has already projected,
#'   coerced, and filtered the data.  Used by the CV internals to avoid a
#'   redundant second pass on every fold.  End users should leave this at the
#'   default \code{FALSE}.
#'
#' @section The local solve:
#' At every observation the model is a weighted least squares on the
#' observations whose kernel weight is above zero.  It is solved by a
#' Householder QR with column pivoting of the window's weighted design, never
#' by inverting \eqn{X^\top W X}{X'WX}.  Before the factorisation each row is
#' multiplied by the square root of its weight over the window's largest
#' weight, and each predictor is divided by its standard deviation over the
#' whole data and centred in the window, so nothing below depends on a
#' predictor's origin or units: an indicator coded 0/1 and one coded
#' 288.7/289.3 give the same fit.
#'
#' One tolerance, \code{1e-7}, decides what a window can estimate:
#' \itemize{
#'   \item \strong{Rank.}  A predictor is aliased in a window when the part
#'     of it that the terms entered before it do not explain varies there by
#'     less than \code{1e-7} of its standard deviation over the study area
#'     (the length of that part over the weighted rows).  A predictor that is
#'     constant among the weighted observations, two predictors that are
#'     exactly collinear there, and a window with fewer weighted observations
#'     than parameters are the cases.  No weight is cut off: under a Gaussian
#'     or exponential kernel every observation enters, and distant ones stop
#'     counting only through this rule.
#'   \item \strong{Coefficients.}  A coefficient is returned when it has one
#'     value over all least-squares solutions of the window, and is \code{NA}
#'     otherwise.  The slope of an aliased predictor is \code{NA}.  So is the
#'     intercept when a predictor is constant in the window at a value other
#'     than 0, because it would absorb that value times the slope; when the
#'     constant is 0 the intercept is identified and returned.  The slopes of
#'     the other predictors are returned.
#'   \item \strong{Fitted values.}  The fitted value at an observation is
#'     always estimable (its own row is in its window) and is returned, also
#'     where a coefficient is \code{NA}.
#' }
#' A window of deficient rank does not stop the fit.  A warning of class
#' \code{"spatialkit_gwr_rank_deficient"} says how many windows are affected
#' and how many coefficients of each term are \code{NA};
#' \code{info$n_rank_deficient}, \code{info$n_local_singular},
#' \code{info$nonfinite_coef} and the \code{rank} column of
#' \code{info$local_collinearity} locate them.  When many windows hold a
#' predictor constant, a local coefficient is the wrong model for it: its
#' effect is identified only as one global coefficient (a mixed GWR, as
#' \code{GWmodel::gwr.mixed()} fits).
#'
#' The leverage of each observation, the trace of the hat matrix, the
#' residual sum of squares, AICc and the coefficients' standard errors come
#' from the same decomposition.  On a window of full rank every one of them
#' equals what \code{GWmodel::gwr.basic()} returns, to rounding error: the
#' largest relative difference in any cell of the six fits the package's
#' tests compare is \code{3e-12}.  They are not bit-identical, since GWmodel
#' inverts the normal equations.
#'
#' @section Bandwidth selection:
#' With \code{bandwidth = NULL} the AICc is scored on 16 bandwidths spaced
#' evenly in the logarithm over the whole range, and then refined by a
#' golden-section search between the grid points either side of the best
#' (and, for a count of neighbours, at the counts next to the one it ends
#' on).  That is about 20 to 35 fits of the model in all.
#' The range runs, for an adaptive bandwidth, from the smallest count the
#' model's parameters allow to every observation; for a fixed one, from the
#' distance at which half of the windows hold one observation more than the
#' model has parameters (a third of it for the Gaussian and exponential
#' kernels) to the largest distance between two observations.  The bandwidth
#' returned is the best of all those scored.
#'
#' A golden-section search over the whole range, which
#' \code{GWmodel::bw.gwr()} runs, assumes the criterion has one minimum.  It
#' has several when the coefficients vary at a range that only a few dozen
#' observations resolve, and that search then returns a bandwidth near the
#' top of the range, a nearly global regression, as the selected one.
#'
#' A candidate bandwidth with windows of deficient rank is scored on its
#' fitted values, which those windows do estimate.  One at which AICc is
#' undefined (\eqn{tr(S)} not below \eqn{n - 2}) cannot be selected, and
#' neither can one at which a term of the model is estimable at no location:
#' the fit would have no coefficient for it anywhere.
#'
#' The criterion at a candidate comes from the local solve.  Under a Gaussian
#' or an exponential kernel, which weight every observation at every
#' location, the fitted value and leverage of a window that is far from rank
#' deficient are computed from its centred cross-products, many windows to
#' one matrix product, and the windows that are not are left to the QR; the
#' two routes agree to rounding error (\code{1e-9} relative in the package's
#' tests) and take the same rank decisions.  A default-bandwidth fit with one
#' of those kernels and two predictors then takes 1.1 to 1.6 times as long
#' as \code{GWmodel::bw.gwr()} followed by \code{gwr.basic()} (2.0 s
#' against 1.6 s at 1,000 observations and 7.4 s against 5.7 s at 2,000,
#' fixed Gaussian, on one machine), and 3 to 5 times as long by the QR
#' alone.  With several strongly correlated predictors more windows are
#' left to the QR and the fit is slower than that.  The fit at the selected
#' bandwidth is always the QR's.
#'
#' Everything scored is kept in \code{info$bandwidth_search}:
#' \code{approach} (the criterion minimised, \code{"AICc"} or \code{"CV"}),
#' \code{by_location} and \code{n_locations} (see \strong{Observations
#' that share a location}), \code{curve}
#' (bandwidth, criterion, trace of the hat matrix, windows of deficient rank,
#' and whether the point belongs to the grid or the refinement),
#' \code{undefined} (how many of the bandwidths scored have no criterion,
#' by reason),
#' \code{n_local_minima} and \code{multimodal} (local minima on the grid),
#' \code{within} and \code{flat} (the range of bandwidths whose AICc is
#' within 2 of the minimum, and whether half the grid is),
#' \code{near_global} (the widest bandwidth of the range is within 2 of the
#' minimum), \code{at_upper} and \code{at_lower} (the minimum is at an end
#' of the range).  When any of these holds, one logged line says so (on the
#' console, unless \code{\link{spatialkit_quiet}()} is on; in the trace only
#' for the folds of a cross-validation); it is a
#' message, not a warning, because a minimum at or near the top of the range
#' is the ordinary outcome for a process whose coefficients do not vary.
#' Read \code{near_global} as "AICc does not tell this fit from one global
#' regression".
#'
#' AICc is an in-sample criterion and rewards interpolation between
#' neighbours.  The bandwidth it selects can score worse than a wider one on
#' held-out spatial blocks; \code{\link{cv_gwr}} with a supplied
#' \code{bandwidth} is the check.
#'
#' @section Observations that share a location:
#' Rows with the same coordinates (several years at one station, replicate
#' samples of one plot, records geocoded to one address) get one kernel
#' weight at every regression point.  \code{info$n_locations} is the number
#' of distinct locations, and when it is below \code{n} a warning of class
#' \code{"spatialkit_gwr_repeated_locations"} gives the counts.  Three things
#' follow.
#'
#' \strong{The bandwidth is not selected by AICc.}  AICc counts the rows at
#' one location as independent observations.  When a location has an effect
#' of its own, which is the usual case, a window that holds little but the
#' location reproduces its mean, the residual sum of squares falls to the
#' spread within locations, and AICc is lowest where the windows
#' interpolate: on 100 simulated sites with 3 observations each it is 58 at
#' 5 to 10 neighbours, 90 at 13 and 221 at 58, while the error of prediction
#' at new sites is 46 at 10 neighbours, 1.35 at 13 and 0.15 at 150.  A floor on
#' the bandwidth does not repair that, because the criterion keeps falling
#' towards the floor.  With \code{bandwidth = NULL} the search therefore
#' minimises the
#' \strong{leave-one-location-out} sum of squares: each location's
#' observations are predicted by the local regression there without any of
#' them, the kernel keeping its width, and the squared errors are summed.
#' It is computed from the hat matrix of the fit at each bandwidth, without
#' refitting.  A bandwidth at which some location cannot be predicted from
#' the others cannot be selected.  \code{info$bandwidth_search$approach} is
#' then \code{"CV"} and \code{by_location} is \code{TRUE}.  On a layer
#' without repeat locations nothing changes: the criterion is AICc.
#'
#' \strong{The range of an adaptive search counts locations.}  It starts at
#' the smallest count that gives every window one location more than the
#' model has parameters, and for a fixed bandwidth at the distance at which
#' half the windows hold that many locations.
#'
#' \strong{AICc and the standard errors still assume independence.}
#' \code{info$AICc}, the comparison with the global regression in
#' \code{summary()} and the standard errors of \code{coef(se = TRUE)} and
#' \code{predict(se.fit = TRUE)} treat every row as an independent
#' observation.  With an effect shared by the rows of a location the standard
#' errors are too small, and AICc compared across bandwidths favours the
#' small ones.  \code{info$AICc} is \code{NA}, with a warning, when the fit
#' (nearly) reproduces the mean of the observations at each location:
#' \eqn{tr(S)} is then not below the number of rows that differ in location
#' or in a predictor, minus 2, and the value the formula gives measures only
#' the spread among observations that share a location.  If the repeats are
#' replicates of one measurement, fit their means, one row per location.
#'
#' @section Collinearity diagnostics:
#' The function computes **scaled condition indices** of the design and
#' warns when one exceeds 30, the conventional threshold, which Wheeler &
#' Tiefelsdorf (2005) carry over to the local designs of GWR.  An index is
#' the ratio of the largest to the smallest singular value (from an SVD) after
#' each column is scaled (Belsley, Kuh & Welsch 1980); an exactly singular
#' design gives `Inf`, which counts as the worst case, not an exempt one.
#' Scaling makes the index independent of the predictors' units;
#' \code{kappa()} on the raw matrix is not, and a threshold on it is a
#' threshold on nothing in particular.
#'
#' A **global** index is computed on the predictors centred at their means
#' (with the intercept, which centring makes orthogonal to them), and kept as
#' `info$condition_index`.  It measures how nearly the predictors are
#' collinear with one another over the whole study area; it is 1 for a single
#' predictor, and a change of origin (degrees C or kelvin, a year or years
#' since 2000) does not move it.
#'
#' **Local** indices are then computed at **every** location, on the weighted
#' design of the local regression there: each row weighted by the square root
#' of its kernel weight at the bandwidth the model is fitted with (supplied or
#' selected), with rows of negligible weight (below 1e-8) dropped.  A window
#' left with fewer rows than columns counts as singular.  Two indices are
#' kept for each window, as columns of `info$local_collinearity`:
#' \describe{
#'   \item{\code{cn}}{Belsley's index of the intercept plus the predictors,
#'     scaled to unit length but not centred.  A predictor whose values in the
#'     window are far from 0 against their spread (a year, a temperature in
#'     kelvin) is collinear with the intercept and raises it: the local
#'     intercept is then an extrapolation to 0 and is ill-determined, but the
#'     slopes are not.}
#'   \item{\code{cn_slopes}}{The index for the slopes: the predictors centred at
#'     their weighted mean in the window and each divided by its standard
#'     deviation over the whole study area.  It is 1 when the predictors vary
#'     as much, and as independently, inside the window as they do across the
#'     study area; it grows as a predictor becomes nearly constant inside the
#'     window (a regional covariate) or two predictors move together there.
#'     It does not depend on the predictors' origin or units.}
#' }
#' A window's slopes count as collinear when `cn_slopes` is above 30 or
#' singular.  Those windows are counted in `info$n_local_collinear`.  A
#' predictor that is constant, or nearly so, inside a window is caught this
#' way whether it is alone or has company, so a single predictor is surveyed
#' too.  `cn` does not enter: the local solve centres each predictor in its
#' window, so a predictor far from zero (a year, a temperature in kelvin, a
#' coordinate) has the slopes it would have at any other origin, and is not
#' flagged, however large `cn` is.  The one exception is a `cn` that is
#' singular to machine precision, where the predictor's spread is lost in the
#' rounding of its own values.
#'
#' A warning is issued whenever **any** location has collinear slopes; the
#' wording reports a percentage when more than 25% of locations are affected
#' and a count otherwise.  Both are real R warnings, not log lines.
#' Coefficients at a near-singular window are least-squares estimates with
#' little to rest on, and can be implausibly large.  An **exactly** singular
#' window (an indicator that is constant inside it, or fewer observations
#' than parameters) is of deficient rank: the coefficients it cannot estimate
#' are `NA` and the fit is returned (see \strong{The local solve}).
#' A window with only `cn` above 30 raises no warning:
#' `plot(fit, type = "coefficients", term = "Intercept")` masks it, and slope
#' maps do not.  Centre such a predictor if you want an interpretable local
#' intercept.  The survey is a caution about how well a coefficient is
#' determined, on rows whose weight is not negligible; the solve's rank is
#' the statement about whether it is determined at all.
#'
#' A further warning counts local regressions that have no fitted value.
#' That happens where the kernel weights are undefined: with an adaptive
#' bandwidth of `k`, a location where `k` or more observations share the same
#' coordinates has a kernel of zero width, and every kernel but the boxcar
#' divides 0 by 0 there.  The warning names that cause when it applies.
#' `fitted()`, `residuals()`, `summary()` and [model_metrics()] all drop those
#' rows, so when this warning fires the metrics describe only the part of the
#' study area that fitted.
#'
#' @return A \code{gwr_fit} object (inherits from \code{spatial_fit}).
#'   Supports \code{predict()}, \code{fitted()}, \code{residuals()},
#'   \code{coef()}, \code{summary()}, and \code{model_metrics()}.
#'   Model-specific metadata lives in \code{$info}: bandwidth, adaptive,
#'   kernel, AICc (\code{NA}, with a warning, where AICc is
#'   undefined: its effective number of parameters \eqn{tr(S)} is not below
#'   \eqn{n - 2}, the local regressions all but interpolate the data, and the
#'   large negative value the formula gives would rank the fit above any
#'   other; also where the fit reproduces the mean of the observations at
#'   each location, see \strong{Observations that share a location}),
#'   \code{bandwidth_is_fallback} (\code{TRUE} when no bandwidth could
#'   be selected: the search failed and the arbitrary fallback was used, or
#'   no bandwidth has a defined criterion and the widest was used),
#'   \code{condition_index} (the global index), \code{local_collinearity}
#'   (one row per observation: \code{row}, \code{x}, \code{y},
#'   \code{n_window}, \code{cn} and \code{cn_slopes}, see
#'   \strong{Collinearity diagnostics}, then \code{rank}, the rank the local
#'   solve found, and \code{leverage}, the observation's diagonal of the hat
#'   matrix), \code{n_local_collinear} (the
#'   locations whose slopes count as collinear), \code{n_local_singular}
#'   (the locations with a local coefficient that is \code{NA}),
#'   \code{n_rank_deficient} (the windows whose rank is below the number of
#'   parameters), \code{n_saturated} (the windows whose fitted value is their
#'   own observation: leverage 1),
#'   \code{nonfinite_coef} (a logical matrix, one row per observation and
#'   one column per term, with \code{Intercept} first, \code{TRUE} where the
#'   local coefficient is \code{NA}, so the count in
#'   \code{n_local_singular} can be placed), \code{n_dropped} (the rows
#'   \code{prep_model_data()} removed for missing or non-finite values or a
#'   bad geometry, so \code{$n} can be read against \code{nrow(data_sf)}),
#'   \code{n_locations} (the distinct locations of the rows fitted; see
#'   \strong{Observations that share a location}),
#'   \code{tr_S}, \code{enp}, \code{edf} and \code{sigma2} (the trace of the
#'   hat matrix, the effective number of parameters
#'   \eqn{2\,tr(S) - tr(S^\top S)}, the effective residual degrees of freedom
#'   and the residual variance, as GWmodel defines them),
#'   \code{global_ols} (\code{RSS}, \code{RMSE} and \code{AICc} of one
#'   least-squares regression on the same rows, the AICc by the same formula
#'   with the regression's rank as its trace; a GWR whose AICc is no lower
#'   has not shown that its coefficients vary, and \code{summary()} prints
#'   the two side by side),
#'   \code{rank_tolerance}, \code{max_kernel_width} (the widest kernel at any
#'   observation, which \code{predict()} measures new locations against) and
#'   \code{bandwidth_search} (see \strong{Bandwidth selection}; \code{NULL}
#'   when a bandwidth was supplied).
#'
#'   \code{$engine} holds the fit in the layout of the object
#'   \code{GWmodel::gwr.basic()} returns, under the same names: \code{SDF} (a
#'   \code{SpatialPointsDataFrame} of the local coefficients, \code{y},
#'   \code{yhat}, \code{residual}, \code{CV_Score}, \code{Stud_residual}, the
#'   standard errors \code{<term>_SE}, the t values \code{<term>_TV} and
#'   \code{Local_R2}), \code{GW.diagnostic}, \code{GW.arguments} and
#'   \code{lm} (the global least-squares fit).  It is a plain list, not a
#'   \code{gwrm} object.
#' @family model fitting
#' @examples
#' if (requireNamespace("sp", quietly = TRUE)) {
#'   library(sf)
#'   set.seed(1)
#'   n <- 60
#'   dat <- st_as_sf(
#'     data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
#'                elev = rnorm(n)),
#'     coords = c("x", "y"), crs = 32632
#'   )
#'   dat$price <- 10 + 0.01 * (st_coordinates(dat)[, 1] - 5e5) +
#'     2 * dat$elev + rnorm(n)
#'   fit <- fit_gwr_model(dat, "price", "elev", bandwidth = 30)
#'   print(summary(fit))     # print(): only a block's last value shows on its own
#'   head(predict(fit, newdata = dat))   # newdata is re-projected if needed
#'   head(coef(fit, se = TRUE)$se)       # local standard errors
#' }
#' @export
fit_gwr_model <- function(data_sf, response_var, predictor_vars,
                          adaptive = TRUE,
                          bandwidth = NULL,
                          kernel = c("bisquare", "gaussian", "tricube",
                                     "boxcar", "exponential"),
                          .already_prepped = FALSE) {
  if (!inherits(data_sf, "sf"))
    stop("fit_gwr_model(): `data_sf` must be an sf object.")
  # GWmodel is not needed: the local regressions, the bandwidth search and
  # the diagnostics are computed here.  sp is, for the Spatial table the
  # fit's engine keeps its results in.
  if (!requireNamespace("sp", quietly = TRUE))
    stop("fit_gwr_model(): package 'sp' is required (the fit keeps its ",
         "local estimates in a SpatialPointsDataFrame). Install with ",
         "install.packages('sp').", call. = FALSE)
  # match.arg() is the whole of kernel validation: an invalid value never gets
  # past it.
  kernel <- match.arg(kernel)

  # character(0) is the model without predictors: a local intercept alone,
  # which is the kernel-weighted mean of the response at each location.  It
  # is the baseline a GWR with predictors has to beat (a GWR beats a global
  # mean through its local intercept whatever its predictors are), and
  # select_features_forward() fits it as its null model.
  if (!is.character(predictor_vars))
    stop("fit_gwr_model(): `predictor_vars` must be a character vector of ",
         "column names (character(0) for the model with a local intercept ",
         "alone).", call. = FALSE)
  # A name given twice is one term: the formula collapses it, so the fit was
  # right, but the collinearity checks ran on the duplicated matrix and warned
  # "exactly singular" and "100% of locations collinear" (Inf condition
  # index), plot(type = "coefficients") then masked every location, and
  # n_params counted the name twice.  gwr_model_selection() already collapses
  # its candidates the same way.
  predictor_vars <- unique(predictor_vars)

  # `bandwidth` is used unchecked below (twice: in the local-collinearity
  # spot-check and as the fitted bandwidth), and the clamping block further
  # down runs only in the adaptive branch.  Unvalidated, NA gives "missing
  # value where TRUE/FALSE needed", a length-2 vector gives "the condition has
  # length > 1", and with adaptive = FALSE a zero or negative distance reaches
  # the kernel untouched.
  if (!is.null(bandwidth) &&
      (!is.numeric(bandwidth) || length(bandwidth) != 1L ||
       !is.finite(bandwidth) || bandwidth <= 0))
    stop(sprintf(paste0("fit_gwr_model(): `bandwidth` must be a single ",
                        "positive finite number (or NULL to select one ",
                        "automatically). With adaptive = %s it is %s."),
                 if (isTRUE(adaptive)) "TRUE" else "FALSE",
                 if (isTRUE(adaptive))
                   "the number of nearest neighbours in each local window"
                 else
                   "a distance in the CRS units of `data_sf`"),
         call. = FALSE)

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
  # Rows prep_model_data() removed on the way in, read from the record it
  # leaves -- which a layer prepared by the caller may still carry.  0 when
  # there is none, and when the record no longer describes this layer: a
  # fold's subset inside cv_*() has neither, and cv_*() reports its own.
  n_dropped <- as.integer(.get_row_record(dat, "dropped")$n %||% 0L)

  # Require plain POINT, mirroring fit_bayesian_spatial_model(), and for the
  # same reason: st_coordinates() on a multi-vertex MULTIPOINT or a POLYGON
  # returns one row per VERTEX, so the local-collinearity spot-check below
  # would index `coords[si, ]` against a matrix that no longer has one row per
  # observation.  prep_model_data() coerces to POINT, so this only fires for
  # `.already_prepped = TRUE` callers handing over unprepped geometry.
  gtypes <- as.character(sf::st_geometry_type(dat, by_geometry = TRUE))
  if (!all(gtypes == "POINT"))
    stop("fit_gwr_model(): geometry must be POINT after prep. ",
         "Run prep_model_data() (or coerce_to_points()) first.", call. = FALSE)

  # A factor or character predictor is refused.  A factor is several
  # model-matrix columns under one name, which the collinearity survey, the
  # coefficient table and gwr_model_selection() (one variable, one column)
  # are not written for; gwr_model_selection() refuses the same columns on
  # this function's account.  Reject the column here, where it can be named.
  non_num <- predictor_vars[
    !vapply(sf::st_drop_geometry(dat)[, predictor_vars, drop = FALSE],
            is.numeric, logical(1))]
  if (length(non_num) > 0L)
    stop(sprintf(paste0("fit_gwr_model(): predictor(s) %s are not numeric. ",
                        "GWmodel would fit a factor or character predictor as ",
                        "one local coefficient per contrast column; this ",
                        "package's GWR functions take numeric predictors only. ",
                        "Encode the column as numeric indicator ",
                        "column(s) yourself, and pass those as predictors."),
                 paste(sQuote(non_num), collapse = ", ")),
         call. = FALSE)


  # With n <= p + 1 no local regression can have a residual degree of freedom
  # -- every window is at best exactly determined -- so the fit that came
  # back was all-NA fitted values behind three warnings (n = 2 warned "only
  # 2 observations", "fallback bandwidth", then "2 of 2 local regressions
  # singular").  One error, up front, before the response is even looked at.
  if (nrow(dat) <= length(predictor_vars) + 1L)
    stop(sprintf(paste0("fit_gwr_model(): %d observation(s) for %d parameters ",
                        "(intercept + %d predictor(s)); a GWR needs at least ",
                        "%d, and far more to be useful."),
                 nrow(dat), length(predictor_vars) + 1L,
                 length(predictor_vars), length(predictor_vars) + 2L),
         call. = FALSE)

  # Warn or error if response looks non-continuous.
  # Gaussian GWR assumes a continuous response; binary data should error.
  #
  # The three degenerate cases are separated because folding them together
  # misdiagnoses two of them: all(numeric(0) == round(numeric(0)), na.rm =
  # TRUE) is TRUE, so an all-dropped dataset used to be reported as "binary
  # (0 unique values)" and a constant response as "binary (1 unique value)",
  # while a genuinely binary non-integer response (1.5 / 2.5) failed the
  # integer-like gate and passed unremarked.
  #
  # REGRESSION NOTE: the two-distinct-value error is gated on `is_integer_like`
  # and must stay that way.  "Exactly two distinct finite values" is NOT the
  # same thing as "binary": a left-censored or saturated measurement (every
  # observation at a detection limit, say 0.0031, or at a ceiling, 12.7401) has
  # two distinct values and is perfectly continuous -- Gaussian GWR on it is a
  # well-defined least-squares problem, and the error's advice to switch to
  # family = "binomial" is nonsense for such values.  The guard also runs once
  # per fold inside cv_gwr(), where a small training fold can legitimately hold
  # only two distinct values.  Non-integer 2-valued responses therefore take
  # the warning path below, not a hard stop.
  resp_vals <- sf::st_drop_geometry(dat)[[response_var]]
  # TRUE/FALSE is the 0/1 coding by another name, and must meet the same
  # guard: it used to bypass it (the guard sat behind is.numeric()) and fit a
  # Gaussian GWR to a binary outcome without a word.
  if (is.logical(resp_vals)) resp_vals <- as.numeric(resp_vals)
  # A factor or character response is refused here, as fit_rf_model() and
  # fit_bayesian_spatial_model() refuse it, with a message that names the
  # column; past this point it would fail in the arithmetic, which does not.
  if (!is.numeric(resp_vals))
    stop(sprintf(paste0("fit_gwr_model(): response '%s' is not numeric (it is ",
                        "%s). GWR here is a Gaussian regression: convert a ",
                        "number stored as text with as.numeric() first; for a ",
                        "categorical outcome use GWmodel::ggwr.basic()."),
                 response_var, class(resp_vals)[1L]),
         call. = FALSE)
  if (is.numeric(resp_vals)) {
    usable   <- resp_vals[is.finite(resp_vals)]
    n_usable <- length(usable)
    uniq     <- unique(usable)
    n_unique <- length(uniq)

    if (n_usable < 2L) {
      stop(
        sprintf("fit_gwr_model(): response '%s' has %d usable (finite, non-missing) value(s) after cleaning; a GWR needs at least two. Check for NA/Inf in the response and predictors.",
                response_var, n_usable),
        call. = FALSE
      )
    }
    if (n_unique < 2L) {
      stop(
        sprintf("fit_gwr_model(): response '%s' is constant (a single distinct finite value, %s). There is nothing to regress.",
                response_var, format(uniq[[1L]])),
        call. = FALSE
      )
    }

    is_integer_like <- all(usable == round(usable))
    if (n_unique == 2L && is_integer_like) {
      stop(
        sprintf("fit_gwr_model(): response '%s' is binary (2 distinct values: %s). Gaussian GWR is invalid for binary outcomes. Consider GWmodel::ggwr.basic() with family = 'binomial'.",
                response_var, paste(format(sort(uniq)), collapse = ", ")),
        call. = FALSE
      )
    }
    if (n_unique == 2L) {
      # Reached only when the response is NOT integer-like (the integer-like
      # case stopped above).  Two distinct non-integer values is the signature
      # of a censored or saturated measurement, which is continuous: fit it,
      # but say that the design is degenerate.
      warning(
        sprintf("fit_gwr_model(): response '%s' has only 2 distinct finite values (%s). The fit is a well-defined least-squares problem, but check that the response is genuinely continuous (e.g. censored at a detection limit) rather than a coded category; if it is categorical, use GWmodel::ggwr.basic() with family = 'binomial'.",
                response_var, paste(format(sort(uniq)), collapse = ", ")),
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
  local_collinearity <- NULL
  # One numeric predictor is enough to survey.  The design is cbind(1, xmat),
  # and a predictor that is (nearly) constant inside a window is collinear
  # with the INTERCEPT -- the case the local survey exists for.  Gating on two
  # predictors skipped it: a regional covariate nearly constant within each
  # cluster gave local slopes of -97 to 221 around a true 3 with no warning,
  # while adding a noise predictor to the same data warned at every location.
  if (length(predictor_vars) >= 1L) {
    num_preds <- predictor_vars[vapply(pred_df[predictor_vars], is.numeric, logical(1))]
    if (length(num_preds) >= 1L) {
      xmat <- as.matrix(pred_df[, num_preds, drop = FALSE])
      # The scaled condition INDEX (Belsley), thresholded at the literature's
      # 30 -- not kappa() on the raw matrix at 1e6, which depends on the
      # predictors' units and let a design with condition index 1322 through.
      # The predictors are CENTRED first.  Uncentred, a predictor far from 0
      # against its spread (kelvin, a year) is collinear with the intercept,
      # so the same field in kelvin warned "collinearity risk" (index 230)
      # where in degrees C it did not, although every slope is the same.  The
      # centred index is the standard one: 1 for a single predictor, and above
      # 30 only when predictors move together over the study area.  A
      # predictor constant over the whole data set is a zero column here
      # (Inf).  Collinearity with the intercept INSIDE a window is the local
      # survey's job.
      cn <- .condition_index(cbind(1, sweep(xmat, 2L, colMeans(xmat))))
      # A non-finite condition number is an EXACTLY singular design -- the worst
      # case there is -- and `is.finite(cn) && ...` silently let it through.
      if (!is.finite(cn) || cn > 30) {
        .warn_and_log(
          "fit_gwr_model(): global predictors (centred) have scaled condition index %s, above the conventional 30 (collinearity risk). Note: local collinearity within bandwidth windows may be substantially worse than this global value.",
          if (is.finite(cn)) sprintf("%.0f", cn) else "infinite (exactly singular)"
        )
      }
      # The LOCAL collinearity spot-check is deferred until the bandwidth is
      # known -- see .gwr_local_collinearity() below the bandwidth
      # selection.  Running it here meant that on the default path
      # (bandwidth = NULL, which is most calls) it used a stand-in window --
      # min(50, n) neighbours when adaptive, the whole data set when fixed --
      # instead of the window the model is actually fitted with, so the
      # documented warning simply never fired: identical data warned when the
      # bandwidth was passed explicitly and stayed silent when the same value
      # was selected automatically.
      local_collinearity <- list(xmat = xmat, num_preds = num_preds, cn = cn)
    }
  }
  
  if (n_obs < n_params * 3L) {
    .log_warn("fit_gwr_model(): only %d observations for %d parameters; local regressions will be underdetermined.",
              n_obs, n_params)
  }
  
  # NOTE: dat is already projected by prep_model_data(), so no
  # additional .is_longlat() / ensure_projected() call is needed.
  
  # reformulate() refuses an empty set of terms; "1" is the intercept alone.
  fml <- stats::reformulate(
    termlabels = if (length(predictor_vars)) predictor_vars else "1",
    response = response_var)

  needed_cols <- unique(c(response_var, predictor_vars))
  sp_dat <- .to_sp(dat, needed_cols)
  # The local regressions are solved here (.gwr_local_solve()), on the model
  # matrix and the two coordinates GWmodel's functions were handed.
  dsn    <- .gwr_design(fml, sp_dat@data)
  X      <- dsn$X
  y      <- dsn$y
  coords <- sp::coordinates(sp_dat)
  # The locations of the layer.  Rows that share one are weighted alike at
  # every regression point; see "Observations that share a location" in the
  # documentation for what follows from that.
  loc        <- .gwr_locations(coords)
  n_distinct <- if (loc$n < n_obs) .gwr_n_distinct(X, loc$id) else n_obs

  # --- Bandwidth selection ---
  bandwidth_is_fallback <- FALSE
  bw_fallback_raw       <- NA_real_
  bw_search             <- NULL
  bw_none_defined       <- FALSE
  # A fixed bandwidth is a length in the CRS the fit runs in, and that is not
  # necessarily the CRS the caller was looking at: prep_model_data() projects
  # geographic input, so 0.2 "degrees" becomes 0.2 metres and every local
  # regression holds its own observation alone.  Say so at the moment the
  # mismatch is visible.
  if (!is.null(bandwidth) && !adaptive)
    .gwr_warn_tiny_fixed_bw(dat, bandwidth, "fit_gwr_model")
  if (is.null(bandwidth)) {
    # A grid first, then a refinement round its minimum (.gwr_bw_search()).
    # A golden-section search over the whole range, GWmodel::bw.gwr()'s,
    # assumes one minimum and returns a near-global bandwidth when the
    # coefficients vary at short range.  cv_gwr(bandwidth = NULL) runs this
    # once per fold.
    bw_search <- tryCatch(
      .gwr_bw_search(X, y, coords, kernel, adaptive, approach = "AICc"),
      error = function(e) {
        .log_warn("fit_gwr_model(): the bandwidth search failed: %s",
                  conditionMessage(e))
        NULL
      }
    )
    bw <- if (is.null(bw_search)) NA_real_ else bw_search$bandwidth

    # Anything that is not one usable positive number takes the fallback.
    if (length(bw) != 1L || !is.numeric(bw) || !is.finite(bw) || bw <= 0) {
      bw <- .fallback_bandwidth(sp_dat, adaptive)
      bandwidth_is_fallback <- TRUE
      bw_fallback_raw       <- as.numeric(bw)
    } else if (isTRUE(bw_search$none_defined)) {
      # Nothing to select: the widest window is used, and it is not an
      # optimised bandwidth, which compare_models() says of a fallback too.
      bandwidth_is_fallback <- TRUE
      bw_none_defined       <- TRUE
      .warn_and_log(paste0("fit_gwr_model(): no bandwidth from %s to %s%s has ",
                           "a defined %s, so none could be selected: at ",
                           "each of the %d tried %s. The widest, %s, is ",
                           "used. Supply `bandwidth`, or use fewer ",
                           "predictors or more observations."),
                    format(signif(bw_search$range[1L], 6)),
                    format(signif(bw_search$range[2L], 6)),
                    if (adaptive) " neighbours" else "",
                    .gwr_bw_criterion_name(bw_search, score = TRUE),
                    bw_search$n_evaluations,
                    .gwr_bw_undefined_text(bw_search, n_obs),
                    format(signif(bw, 6)))
    } else {
      # What the grid showed, in one line: a message on the console for a
      # direct call, never a warning (a minimum at the top of the range is
      # the ordinary outcome when the coefficients do not vary).  A fold of a
      # cross-validation writes it to the trace only, or it would be printed
      # once per fold.
      note <- .gwr_bw_search_note(bw_search, adaptive)
      if (nzchar(note))
        (if (isTRUE(.already_prepped)) .log_info else .log_warn)(
          paste0("fit_gwr_model(): bandwidth %s%s selected by %s ",
                 "from %d scored between %s and %s: %s See ",
                 "info$bandwidth_search."),
          format(signif(bw, 6)), if (adaptive) " neighbours" else "",
          .gwr_bw_criterion_name(bw_search),
          bw_search$n_evaluations,
          format(signif(bw_search$range[1L], 6)),
          format(signif(bw_search$range[2L], 6)), note)
    }
  } else {
    bw <- as.numeric(bandwidth)
  }

  # Clamp bandwidth to safe range
  if (adaptive) {
    # The up-front check admits any positive finite number, which is right for
    # a fixed-distance bandwidth but not for an adaptive one: a distance in
    # metres passed while forgetting adaptive = FALSE can exceed
    # .Machine$integer.max, as.integer() makes it NA, and the clamp below --
    # which exists precisely to cap an over-large adaptive bandwidth at n_obs
    # -- aborts on "missing value where TRUE/FALSE needed" instead.
    .check_scalar(bw, "bandwidth", "fit_gwr_model", min = 1,
                  max = .Machine$integer.max,
                  what = "a single number of nearest neighbours when adaptive = TRUE")
    bw <- as.integer(round(bw))
    # n_params + 1 left bisquare and tricube windows with n_params points
    # (their bw-th neighbour has weight 0), so a small supplied bandwidth was
    # clamped onto an exact interpolator: R2 = 1 and an AICc of -15033 that
    # compare_models() ranked first.  See .gwr_min_adaptive_bw().  Capped at
    # n_obs: with n = p + 2 there is no larger window, and the AICc guard
    # below reports what is left.
    min_bw <- min(.gwr_min_adaptive_bw(n_params, kernel), n_obs)
    max_bw <- n_obs
    if (bw < min_bw) {
      # A warning, not a log line: the fit is not at the bandwidth the caller
      # asked for.  The fallback's own warning already names its clamp.
      if (bandwidth_is_fallback)
        .log_warn("fit_gwr_model(): adaptive bandwidth %d too small for %d params; clamping to %d.",
                  bw, n_params, min_bw)
      else
        .warn_and_log(paste0("fit_gwr_model(): an adaptive bandwidth of %d ",
                             "neighbours is too small for %d parameters with ",
                             "the %s kernel (a local regression needs more ",
                             "points with non-zero weight than parameters); ",
                             "using %d.%s"),
                      bw, n_params, kernel, min_bw,
                      if (kernel %in% c("bisquare", "tricube"))
                        paste0(" That is enough unless several neighbours tie ",
                               "at the kernel's edge (a regular grid), which ",
                               "gives them weight 0 too; then use a larger ",
                               "bandwidth.") else "")
      bw <- min_bw
    }
    # Capped at n, and said.  A supplied count above n is most often a
    # distance passed with adaptive left at TRUE: 1500 (metres) on 200 points
    # became a near-global fit whose local slopes varied half as much as the
    # intended one's.  The search never goes past n, and the fallback's own
    # warning names its clamp.
    if (bw > max_bw) {
      if (!is.null(bandwidth))
        .warn_and_log(paste0("fit_gwr_model(): an adaptive bandwidth of %d ",
                             "neighbours exceeds the %d observations; using ",
                             "%d. With adaptive = TRUE the bandwidth is a ",
                             "count of neighbours; for a distance in CRS ",
                             "units, set adaptive = FALSE."),
                      bw, n_obs, n_obs)
      bw <- max_bw
    }
  }

  # The fallback warning is issued AFTER the clamp, not before it.  The clamp
  # can move the value -- a 6-row dataset clamps the arbitrary 10 down to 6 --
  # and a warning naming a bandwidth the fit did not use sends the reader
  # looking for a number that appears nowhere in the result.  Both are named
  # when they differ.
  if (bandwidth_is_fallback && !bw_none_defined) {
    bw_used <- as.numeric(bw)
    warning(
      sprintf(
        paste0("fit_gwr_model(): automatic bandwidth selection failed; using ",
               "arbitrary fallback bandwidth = %.4f%s. This fallback has no ",
               "relationship to the data's spatial structure and may produce ",
               "a poor fit. Consider supplying an explicit `bandwidth` ",
               "argument."),
        bw_used,
        if (isTRUE(bw_used == bw_fallback_raw)) ""
        else sprintf(" (derived as %.4f, then clamped for %d observation(s) and %d parameter(s))",
                     bw_fallback_raw, n_obs, n_params)
      ),
      call. = FALSE
    )
  }

  # Now that the bandwidth is known, survey the local collinearity of every
  # window the model will actually be fitted with, and keep the per-location
  # values: the coefficient map masks on them, and the warning below is the
  # exact fraction rather than a sampled one.
  global_cn <- if (!is.null(local_collinearity)) local_collinearity$cn else NA_real_
  local_cn_df <- NULL
  if (!is.null(local_collinearity)) {
    local_cn_df <- .gwr_local_collinearity(sf::st_coordinates(dat),
                                           local_collinearity$xmat,
                                           adaptive, bw, kernel)
    .gwr_local_collinearity_warn(local_cn_df)
  }

  # --- Fit GWR ---
  # One pivoted QR per window.  A window of deficient rank does not stop the
  # fit and does not return rounding error: the coefficients it cannot
  # estimate are NA, and its fitted value is returned when it is estimable,
  # which at a training point it always is.
  sv <- .gwr_local_solve(X, y, coords, bw, kernel, adaptive, what = "fit")
  dg <- .gwr_diagnostics(y, sv$fitted, sv$leverage, sv$var_fit)
  fit <- .gwr_engine(sv, dg, fml, sp_dat, X, y, bw, kernel, adaptive)

  # The per-row, per-term mask of local coefficients that are not estimates
  # (`info$nonfinite_coef`), not only its row count: the count says how many
  # local regressions are affected, the mask says which, and which term.
  nonfinite_mask <- !is.finite(sv$coef)
  bad_rows    <- apply(nonfinite_mask, 1L, any)
  n_bad_local <- sum(bad_rows)
  no_fit      <- is.na(sv$fitted)
  deficient   <- sv$rank < n_params
  n_saturated <- sum(sv$leverage >= 1 - .gwr_rank_tol(), na.rm = TRUE)

  # A window with no fitted value.  Nothing downstream says so by itself:
  # fitted(), residuals(), summary() and model_metrics() all drop the rows
  # without one, so a fit in which 182 of 200 local regressions had none
  # reported n = 18 and R2 = 0.96 as though that were the whole model.  The
  # cause is undefined kernel weights: an adaptive bandwidth whose bw-th
  # nearest neighbour is at distance zero -- bw or more observations at one
  # location -- gives every kernel but the boxcar 0/0, and no row of such a
  # window carries weight.
  if (any(no_fit)) {
    n_zero_width <- if (isTRUE(adaptive)) sum(no_fit & sv$n_window == 0L) else 0L
    cause <- if (n_zero_width > 0L)
      sprintf(paste0(" -- at %d of them %d or more observations share one ",
                     "location, so the adaptive kernel has zero width and its ",
                     "weights are 0/0. Use a bandwidth larger than the most ",
                     "observations at any one location, or aggregate repeat ",
                     "observations first"), n_zero_width, as.integer(bw))
    else " (undefined kernel weights)"
    .warn_and_log(paste0("fit_gwr_model(): %d of %d local regression(s) ",
                         "returned non-finite coefficients%s. fitted(), ",
                         "residuals() and every metric are computed from the ",
                         "%d that succeeded, so they describe part of the ",
                         "study area only."),
                  sum(no_fit), n_obs, cause, n_obs - sum(no_fit))
  }

  # Windows of deficient rank that do have a fitted value: some coefficients
  # are NA there, and the caller has to be told which and how many.
  part <- deficient & !no_fit
  if (any(part)) {
    by_term <- colSums(nonfinite_mask[part, , drop = FALSE])
    sat_part <- sum(part & sv$leverage >= 1 - .gwr_rank_tol(), na.rm = TRUE)
    .gwr_warn(
      "spatialkit_gwr_rank_deficient",
      paste0("fit_gwr_model(): %d of %d local regressions are rank deficient ",
             "at the bandwidth in use (%s%s, %s kernel): among the ",
             "observations that carry weight there a predictor does not ",
             "vary, the predictors are exactly collinear, or fewer ",
             "observations than the %d parameters carry weight. The ",
             "coefficients those windows cannot estimate are NA (%s) and are ",
             "counted in info$n_local_singular; their fitted values are ",
             "estimable and are returned%s.%s"),
      sum(part), n_obs, format(bw), if (adaptive) " neighbours" else "",
      kernel, n_params,
      paste(sprintf("%s at %d", colnames(nonfinite_mask), by_term),
            collapse = ", "),
      if (sat_part > 0L)
        sprintf(paste0("; %d of those windows hold as many weighted ",
                       "observations as they have estimable coefficients, so ",
                       "the fitted value there is the observation itself"),
                sat_part) else "",
      .gwr_rank_hint(
        selected = if (is.null(bandwidth) && !bandwidth_is_fallback)
          format(bw) else NULL))
  }

  # AICc extraction
  AICc_val <- dg$AICc
  # Outside its domain (tr S >= n - 2) the AICc formula gives a large negative
  # number, -15033 for a 100-point fit whose windows interpolate, and
  # compare_models() and print() took it as the best fit there is.  Report
  # NA and say why; see .gwr_aicc_undefined().
  if (isTRUE(.gwr_aicc_undefined(dg$AIC, dg$AICc)) ||
      isTRUE(dg$tr_S >= n_obs - 2)) {
    trS <- dg$tr_S
    # "Use a larger bandwidth" cannot be followed once an adaptive one is
    # already every observation: the remedy is then fewer parameters, more
    # data, or a kernel that does not give the farthest point weight 0.
    remedy <- if (isTRUE(adaptive) && bw >= n_obs)
      sprintf(paste0("Even the widest adaptive window (all %d observations) ",
                     "leaves too few residual degrees of freedom for %d ",
                     "parameters: use fewer predictors, more observations, ",
                     "or a gaussian or exponential kernel."),
              n_obs, n_params)
    else
      sprintf("The bandwidth (%s%s) is too small for %d parameters; use a larger one.",
              format(bw), if (adaptive) " neighbours" else "", n_params)
    .warn_and_log(paste0("fit_gwr_model(): AICc is undefined for this fit ",
                         "and is reported as NA. Its effective number of ",
                         "parameters, tr(S) = %s, is not below n - 2 = %d, ",
                         "so the local regressions (nearly) interpolate the ",
                         "data, and there the AICc formula turns ",
                         "negative: its value (%s) would rank this fit above ",
                         "better ones. %s"),
                  if (is.na(trS)) "?" else format(round(trS, 2)),
                  n_obs - 2L, format(signif(AICc_val, 6)), remedy)
    AICc_val <- NA_real_
  } else if (n_distinct < n_obs && isTRUE(dg$tr_S >= n_distinct - 2)) {
    # Rows that share a location and every predictor value get one fitted
    # value, so this fit is a fit to the means of the groups they form, and
    # it has (all but) one parameter for each mean.  The formula still gives
    # a number, because it counts every row; that number measures the spread
    # inside the groups and nothing the bandwidth is responsible for.
    .warn_and_log(paste0("fit_gwr_model(): AICc is reported as NA for this ",
                         "fit. The %d observations are %d distinct ",
                         "combinations of location and predictor values, and ",
                         "the fit's effective number of parameters, tr(S) = ",
                         "%s, is not below that number minus 2 (%d): the ",
                         "local regressions (nearly) reproduce the mean of ",
                         "the observations at each location, so the value ",
                         "the formula gives (%s) measures only the spread ",
                         "among observations that share a location, and it ",
                         "can rank this fit above one that predicts at a ",
                         "new location. Use a larger bandwidth."),
                  n_obs, n_distinct, format(round(dg$tr_S, 2)),
                  n_distinct - 2L, format(signif(AICc_val, 6)))
    AICc_val <- NA_real_
  }
  if (!is.finite(AICc_val)) AICc_val <- NA_real_

  # Several observations at one location.  Said once, with the counts: the
  # numbers below that assume independent observations, and how the
  # bandwidth was chosen.  A fold of a cross-validation writes it to the
  # trace only, as it does the note on the search.
  if (loc$n < n_obs) {
    repeat_msg <- sprintf(
      paste0("fit_gwr_model(): the %d observations are at %d distinct ",
             "location(s) (%d hold more than one observation, the most %d). ",
             "AICc, the standard errors and the comparison with the global ",
             "regression count the observations at one location as ",
             "independent. They ",
             "are not when a location has an effect of its own: AICc then ",
             "falls as the bandwidth shrinks and the standard errors are ",
             "too small.%s If the repeats are replicates of one ",
             "measurement, fit their means, one row per location."),
      n_obs, loc$n, sum(loc$size > 1L), max(loc$size),
      if (is.null(bandwidth) && !bandwidth_is_fallback &&
          isTRUE(bw_search$by_location))
        sprintf(paste0(" No bandwidth was supplied: it was selected by %s ",
                       "(each location predicted from the others), not by ",
                       "AICc."), .gwr_bw_criterion_name(bw_search))
      else "")
    if (isTRUE(.already_prepped)) .log_info("%s", repeat_msg)
    else .gwr_warn("spatialkit_gwr_repeated_locations", "%s", repeat_msg)
  }

  # The benchmark a GWR is read against: one least-squares regression on the
  # same rows, scored by the same AICc formula (its hat-matrix trace is its
  # rank).  A GWR whose AICc is no lower has not shown that its coefficients
  # vary.
  ols_rss  <- sum(stats::residuals(fit$lm)^2)
  ols_rank <- fit$lm$rank
  global_ols <- list(
    RSS  = ols_rss,
    RMSE = sqrt(ols_rss / n_obs),
    AICc = if (isTRUE(ols_rss > 0) && n_obs - 2 - ols_rank > 0)
      n_obs * log(ols_rss / n_obs) + n_obs * log(2 * pi) +
        n_obs * ((n_obs + ols_rank) / (n_obs - 2 - ols_rank)) else NA_real_)

  # The survey and the solve, side by side for each window: the condition
  # indices, then the rank the solve found and the leverage of the window's
  # own observation.
  if (is.null(local_cn_df))
    local_cn_df <- data.frame(row = seq_len(n_obs), x = coords[, 1L],
                              y = coords[, 2L], n_window = sv$n_window,
                              cn = NA_real_, cn_slopes = NA_real_)
  n_local_collinear <- if (is.null(local_collinearity)) NA_integer_ else
    sum(.gwr_slopes_collinear(local_cn_df, own_solve = TRUE))
  local_cn_df$rank     <- sv$rank
  local_cn_df$leverage <- sv$leverage

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
      bandwidth_is_fallback = bandwidth_is_fallback,
      # The global scaled condition index of the centred numeric predictors
      # (with the intercept); the number to compare across candidate
      # predictor sets.  NA only when there is no numeric predictor, which
      # the non-numeric refusal above makes unreachable.
      condition_index       = global_cn,
      # One row per observation: the two condition indices (cn, with the
      # intercept and uncentred; cn_slopes, for the slopes) of the
      # kernel-weighted local design at that location, the rank the local
      # solve found there and the leverage of the location's own observation.
      local_collinearity    = local_cn_df,
      # The locations whose SLOPES are collinear; an ill-determined
      # intercept alone (cn > 30 only) is not counted.
      n_local_collinear     = n_local_collinear,
      # The locations with a coefficient that is not an estimate (NA), the
      # windows of deficient rank among them, and the windows whose fitted
      # value is their own observation.
      n_local_singular      = n_bad_local,
      n_rank_deficient      = sum(deficient),
      n_saturated           = n_saturated,
      # One row per observation, one column per term (Intercept first):
      # TRUE where the local coefficient is not an estimate.
      nonfinite_coef        = nonfinite_mask,
      n_dropped             = n_dropped,
      # The distinct locations of the training rows: fewer than `n` when
      # observations share a location.
      n_locations           = loc$n,
      # The global numbers of the fit, from the same decomposition as the
      # coefficients: trace of the hat matrix, effective number of parameters
      # and degrees of freedom, residual variance.
      tr_S                  = dg$tr_S,
      enp                   = dg$enp,
      edf                   = dg$edf,
      sigma2                = dg$sigma2,
      # One global least-squares fit on the same rows, by the same criterion.
      global_ols            = global_ols,
      rank_tolerance        = .gwr_rank_tol(),
      # The widest kernel at any training location: predict() compares a new
      # location's distance to its nearest training point with it.
      max_kernel_width      = suppressWarnings(max(sv$width, na.rm = TRUE)),
      # What the bandwidth search scored, NULL when a bandwidth was supplied.
      bandwidth_search      = bw_search
    )
  )
}


#' The distance parameter of the kernel at one regression point
#'
#' The bandwidth itself when it is a distance; for an adaptive bandwidth the
#' distance to the \code{bw}-th nearest point, and \code{bw / length(d)}
#' times the largest distance when \code{bw} exceeds the number of points,
#' as GWmodel has it.
#'
#' @param d Numeric vector of distances from the regression point.
#' @param bw Bandwidth: a distance, or a neighbour count when adaptive.
#' @param adaptive Logical.
#' @return One number.
#' @keywords internal
#' @noRd
.gw_kernel_width <- function(d, bw, adaptive) {
  if (!isTRUE(adaptive)) return(as.numeric(bw))
  k <- max(1L, as.integer(round(bw)))
  if (k <= length(d)) sort.int(d, partial = k)[k] else k / length(d) * max(d)
}


#' Kernel weights of GWmodel's five kernels
#'
#' The same definitions, written the same way, as \code{GWmodel::gw.weight()}
#' (its C++ \code{gw_weight_*()}), so a window of the local solve, and of the
#' collinearity survey, is the window \code{GWmodel::gwr.basic()} would
#' weight.  For an adaptive bandwidth \code{bw}
#' is a neighbour count and the kernel's distance parameter is the distance to
#' the \code{bw}-th nearest point; a count above \code{length(d)} widens the
#' kernel to \code{bw / length(d)} times the largest distance, as GWmodel does.
#'
#' Two details are GWmodel's and matter.  The truncated kernels test
#' \code{d > h}, so a point exactly at the kernel's edge is inside a boxcar
#' window (weight 1) and on the rim of a bisquare or tricube one (weight 0).
#' And a zero-width kernel (\code{h = 0}: the \code{bw} nearest points share
#' one location) is not an indicator of those points: every kernel but the
#' boxcar gives them \code{0/0 = NaN}, and a window whose weights are
#' undefined has no weighted row and no local regression.
#'
#' @param d Numeric vector of distances from the regression point.
#' @param bw Bandwidth: a distance, or a neighbour count when adaptive.
#' @param kernel One of the validated kernel names.
#' @param adaptive Logical.
#' @param h The kernel's distance parameter, when the caller has it already
#'   (\code{.gw_kernel_width()}).
#' @return Numeric weights, one per element of \code{d} (\code{NaN} where
#'   GWmodel's are).
#' @keywords internal
#' @noRd
.gw_kernel_weights <- function(d, bw, kernel, adaptive,
                               h = .gw_kernel_width(d, bw, adaptive)) {
  # The spelling below is GWmodel's, not a simplification of it: d / h would
  # give the same weights for h > 0 but not GWmodel's NaN for h = 0.
  switch(kernel,
         gaussian    = exp(d^2 / (-2 * h^2)),
         exponential = exp(-d / h),
         tricube     = ifelse(d > h, 0, (1 - d^3 / h^3)^3),
         boxcar      = ifelse(d > h, 0, 1),
         ifelse(d > h, 0, (1 - d^2 / h^2)^2))       # bisquare
}


#' The rows of a window and their kernel weights
#'
#' What \code{which(w > 0)} and \code{w[w > 0]} give for
#' \code{w <- .gw_kernel_weights(d, , kernel, , h = h)}, without computing a
#' weight for every point: the three kernels with an edge are evaluated only
#' at the points not beyond it.  The values are the same to the last bit (the
#' expressions are those of \code{.gw_kernel_weights()}), an undefined weight
#' (a zero-width kernel) is dropped as \code{w > 0} drops it, and the rows
#' come back in their order in \code{d}.
#'
#' @param d Numeric vector of distances from the regression point.
#' @param h The kernel's distance parameter (\code{.gw_kernel_width()}).
#' @param kernel One of the validated kernel names.
#' @return A list: \code{rows} (integer indices into \code{d}) and \code{w}
#'   (their weights, all above zero).
#' @keywords internal
#' @noRd
.gw_window <- function(d, h, kernel) {
  if (kernel == "gaussian" || kernel == "exponential") {
    w <- if (kernel == "gaussian") exp(d^2 / (-2 * h^2)) else exp(-d / h)
    rows <- which(w > 0)
    return(list(rows = rows, w = w[rows]))
  }
  inside <- which(!(d > h))
  di <- d[inside]
  w <- switch(kernel,
              tricube = (1 - di^3 / h^3)^3,
              boxcar  = rep(1, length(di)),
              (1 - di^2 / h^2)^2)                       # bisquare
  keep <- which(w > 0)
  list(rows = inside[keep], w = w[keep])
}


#' The one tolerance of the local solve
#'
#' Used twice in \code{.gwr_local_solve()}, on the same triangular factor: a
#' diagonal entry below it ends the window's rank, and a coefficient or a
#' fitted value is estimable when what the window cannot place of it is no
#' larger than this times its largest term.  \code{stats::lm()} uses the same
#' number for its own rank decision.
#'
#' @return \code{1e-7}.
#' @keywords internal
#' @noRd
.gwr_rank_tol <- function() 1e-7


#' One weighted least squares per regression point, by pivoted QR
#'
#' The local regression at a point \eqn{u} minimises
#' \eqn{\sum_r w_r(u)\,(y_r - x_r^\top\beta)^2} over the training rows whose
#' kernel weight is above zero.  It is solved here, for every row of
#' \code{rp}, by a Householder QR with column pivoting (LAPACK
#' \code{dgeqp3}, through \code{qr(LAPACK = TRUE)}) of the window's weighted
#' design, and never by inverting \eqn{X^\top W X}: the inverse squares the
#' design's condition number, and on a window of deficient rank it returns
#' rounding error or stops, according to how the predictor happens to be
#' coded.
#'
#' \strong{The design that is factored.}  Rows are multiplied by
#' \eqn{\sqrt{w / \max w}}, so the largest weight in a window is 1 whatever
#' the kernel's scale.  Each predictor is divided by its standard deviation
#' over the whole training data and centred, first on the regression point's
#' own value (an exact subtraction for nearby values, which leaves a predictor
#' that is constant in the window exactly zero) and then on its weighted mean
#' in the window.  The centred columns are orthogonal to the intercept, so
#' neither the rank nor any result depends on a predictor's origin or units:
#' a temperature in kelvin and in degrees C, an indicator coded 0/1 and
#' 288.7/289.3, give the same ranks, fitted values and slopes.
#'
#' \strong{The rank rule.}  With column pivoting the diagonal of the
#' triangular factor does not increase, and each entry is the length of what
#' is left of the next column once the columns already entered are projected
#' out.  The window's rank is the number of leading entries of at least
#' \code{tol}.  In words: a predictor is aliased in a window when the part of
#' it that the terms before it do not explain varies there, summed over the
#' weighted rows, by less than \code{tol} of its standard deviation over the
#' study area.  No weight is cut off: under a Gaussian or an exponential
#' kernel every row enters, and distant rows whose weights are tiny against
#' the window's largest stop counting only through this rule.
#'
#' \strong{Estimability.}  In a window of deficient rank a linear function
#' \eqn{c^\top\beta} of the coefficients has one value for every least-squares
#' solution exactly when \eqn{c} lies in the row space of the weighted
#' design.  That is tested on the same factor with the same \code{tol}: with
#' \eqn{R = [R_{11}\; R_{12}]} in pivot order, the part of \eqn{c} the window
#' cannot place is \eqn{c_2 - R_{12}^\top R_{11}^{-\top} c_1}, and \eqn{c} is
#' estimable when no entry of it exceeds \code{tol} times the largest entry
#' of \eqn{c}.  Every original coefficient is tested (its row of the
#' transformation back from the centred, scaled columns) and reported
#' \code{NA} when it fails, not only the column the pivoting left out: with a
#' predictor constant at \eqn{v} in a window its slope is aliased and so is
#' the intercept, which would otherwise absorb \eqn{v} times the slope, unless
#' \eqn{v = 0}.  The fitted or predicted value is tested the same way, with
#' \eqn{c} the regression point's own row: it is always estimable at a
#' training point (its row is in its own window), and at a new location it is
#' \code{NA} when a predictor there takes a value the window's rows cannot
#' place.
#'
#' \strong{What comes from the same decomposition.}  With
#' \eqn{t = R_{11}^{-\top} c_1} and \eqn{Q_1} the first \code{rank} columns of
#' \eqn{Q}, \eqn{c^\top\hat\beta = t^\top Q_1^\top W^{1/2} y}, so the weights
#' the estimate puts on the observations are \eqn{W^{1/2} Q_1 t}.  For the
#' regression point's own row those weights are its row of the hat matrix:
#' the leverage is its entry at the point itself, and its sum of squares is
#' the point's share of \eqn{tr(S^\top S)}.  For a coefficient the sum of
#' squares is the factor its variance multiplies \eqn{\sigma^2} by, the
#' diagonal of \eqn{C C^\top} with \eqn{C = (X^\top W X)^{-1} X^\top W} in
#' GWmodel's notation.  On a window of full rank all of these equal
#' \code{GWmodel::gwr.basic()}'s to rounding error.
#'
#' @param X Training design matrix, the intercept (a column of ones) first.
#'   A design without an intercept is scaled by each column's root mean
#'   square and not centred.
#' @param y Response, one value per row of \code{X}.
#' @param coords Two-column matrix of the training coordinates.
#' @param bw,kernel,adaptive The bandwidth settings; see
#'   \code{.gw_kernel_weights()}.  An adaptive bandwidth counts neighbours
#'   among the training points.
#' @param rp Two-column matrix of regression points, or \code{NULL} for the
#'   training points themselves.
#' @param X_rp Design matrix at \code{rp} (same columns as \code{X}).
#' @param what \code{"fit"}: coefficients, their variance factors, fitted
#'   values, hat rows.  \code{"score"}: fitted values and leverage only, for
#'   the bandwidth search and the selection sweep.  \code{"predict"}: fitted
#'   values at \code{rp}, with their variance factor when \code{se = TRUE}.
#' @param se For \code{what = "predict"}: also return the variance factor of
#'   each prediction.
#' @param tol The tolerance; see above.
#' @param loc The training layer's locations (\code{.gwr_locations()}), for
#'   \code{what = "score"} when the leave-one-location-out score is wanted;
#'   \code{NULL} otherwise.
#' @return A list.  Always \code{fitted} (\code{NA} where not estimable),
#'   \code{rank} (0 where no row carries weight), \code{n_window} (rows with
#'   weight above zero), \code{n_coef} and \code{width} (the kernel's distance
#'   parameter at each regression point).  When the regression points are the
#'   training points, \code{leverage} and \code{n_estimable} (for each term,
#'   the number of windows that estimate it); with \code{loc} on a layer whose
#'   locations repeat, \code{hat_loc} (for each row at a shared location, its
#'   hat-matrix entries on the rows of that location, in row order).  For
#'   \code{"fit"} also \code{coef} and
#'   \code{coef_var} (one row per regression point, one column per term,
#'   \code{NA} where not estimable; multiply \code{coef_var} by
#'   \eqn{\sigma^2} for the variance), \code{var_fit} (the fitted value's
#'   variance factor, its hat row's sum of squares), \code{q_diag} (the
#'   diagonal of \eqn{(I - S)^\top (I - S)}) and \code{local_r2} (GWmodel's
#'   local \eqn{R^2}).  For \code{"predict"} with \code{se = TRUE},
#'   \code{var_fit}.
#' @keywords internal
#' @noRd
.gwr_local_solve <- function(X, y, coords, bw, kernel, adaptive,
                             rp = NULL, X_rp = NULL,
                             what = c("fit", "score", "predict"),
                             se = FALSE, tol = .gwr_rank_tol(), loc = NULL) {
  what <- match.arg(what)
  self <- is.null(rp)
  if (self) { rp <- coords; X_rp <- X }
  n <- nrow(X); kk <- ncol(X); m <- nrow(rp)
  icpt <- kk >= 1L && all(X[, 1L] == 1)
  jp <- if (icpt) seq_len(kk)[-1L] else seq_len(kk)
  p  <- length(jp)
  s  <- rep(1, p)
  if (p) {
    s <- if (icpt) apply(X[, jp, drop = FALSE], 2L, stats::sd)
         else sqrt(colMeans(X[, jp, drop = FALSE]^2))
    # A predictor with no spread over the whole data has no scale to offer;
    # it is exactly constant in every window and aliased there whatever the
    # divisor.
    s[!is.finite(s) | s <= 0] <- 1
  }
  inv_s <- 1 / s
  Xs    <- sweep(X[, jp, drop = FALSE], 2L, s, "/")
  Xs_rp <- sweep(X_rp[, jp, drop = FALSE], 2L, s, "/")
  cx <- coords[, 1L]; cy <- coords[, 2L]; rx <- rp[, 1L]; ry <- rp[, 2L]
  full   <- what == "fit"
  need_q <- full || (what == "predict" && isTRUE(se))

  yhat  <- rep(NA_real_, m)
  rk    <- integer(m)
  nw    <- integer(m)
  width <- rep(NA_real_, m)
  near  <- if (self) NULL else rep(NA_real_, m)
  # A regression point without a location or with a missing predictor has
  # no local regression; its row stays NA.
  usable <- is.finite(rx) & is.finite(ry) &
    rowSums(!is.finite(X_rp)) == 0L
  lev   <- if (self) rep(NA_real_, m) else NULL
  vfit  <- if (need_q) rep(NA_real_, m) else NULL
  # The windows at which each term is estimable, counted for the training
  # points: a term that no window estimates has no coefficient anywhere.
  nest  <- if (self) integer(kk) else NULL
  # Each training row's hat entries on the rows that share its location, for
  # the leave-one-location-out score (see .gwr_lolo_residuals()).  Kept only
  # where a location holds more than one row.
  by_loc <- self && !need_q && !is.null(loc) && loc$n < n
  if (by_loc) {
    loc_rows <- split(seq_len(n), loc$id)
    shared   <- loc$size[loc$id] > 1L
    hat_loc  <- vector("list", n)
  }
  if (full) {
    B <- matrix(NA_real_, m, kk, dimnames = list(NULL, colnames(X)))
    V <- B
    qd <- numeric(n); tssw <- numeric(n); rssw <- numeric(n)
    dyb <- (y - mean(y))^2
  }
  # Each original coefficient as a function of the centred, scaled ones: a
  # slope is its scaled slope over the predictor's standard deviation, and
  # the intercept takes every centring back out.
  back <- function(xf, mu)
    if (!p) matrix(1, 1L, 1L)
    else if (icpt) rbind(c(1, -(xf + mu)), cbind(0, diag(inv_s, p, p)),
                         deparse.level = 0L)
    else diag(inv_s, p, p)

  todo <- which(usable)
  # A kernel that weights every row: the score of the windows that are far
  # from rank deficient comes from their cross-products, a block at a time
  # (.gwr_score_dense()); the loop below takes the others.
  if (self && what == "score" && icpt && n > 1L &&
      kernel %in% c("gaussian", "exponential")) {
    ds  <- .gwr_score_dense(X, y, coords, bw, kernel, adaptive,
                            loc = if (by_loc) loc else NULL)
    hit <- which(ds$done & usable)
    yhat[hit]  <- ds$fitted[hit]
    lev[hit]   <- ds$leverage[hit]
    rk[hit]    <- kk
    nw[hit]    <- ds$n_window[hit]
    width[hit] <- ds$width[hit]
    nest <- nest + length(hit)
    if (by_loc) hat_loc[hit] <- ds$hat_loc[hit]
    todo <- setdiff(todo, hit)
  }

  for (i in todo) {
    d <- sqrt((cx - rx[i])^2 + (cy - ry[i])^2)
    if (!self) near[i] <- min(d)
    h <- .gw_kernel_width(d, bw, adaptive)
    width[i] <- h
    # Rows with weight above zero.  An undefined weight (0/0 from a
    # zero-width adaptive kernel) is not above zero, so such a window has no
    # rows, rank 0 and nothing estimable.
    win  <- .gw_window(d, h, kernel)
    rows <- win$rows
    mk <- length(rows)
    nw[i] <- mk
    if (!mk) next
    wk <- win$w
    wt <- wk / max(wk)
    sw <- sqrt(wt)
    xf <- NULL; mu <- NULL
    if (!p) {
      Z <- matrix(sw, ncol = 1L)
      cfit <- 1
    } else if (icpt) {
      xf <- Xs_rp[i, ]
      D  <- Xs[rows, , drop = FALSE] - rep(xf, each = mk)
      mu <- colSums(wt * D) / sum(wt)
      Z  <- cbind(sw, sw * (D - rep(mu, each = mk)), deparse.level = 0L)
      # The regression point's own row in these columns: its predictors are
      # the centre of the first subtraction, so they are 0 - mu.
      cfit <- c(1, -mu)
    } else {
      xf <- Xs_rp[i, ]
      Z  <- sw * Xs[rows, , drop = FALSE]
      cfit <- xf
    }
    q   <- qr(Z, LAPACK = TRUE)
    qm  <- q$qr
    piv <- q$pivot
    nr  <- min(mk, kk)
    dg  <- abs(qm[cbind(seq_len(nr), seq_len(nr))])
    r <- 0L
    while (r < nr && is.finite(dg[r + 1L]) && dg[r + 1L] >= tol) r <- r + 1L
    rk[i] <- r
    if (!r) next
    i1  <- seq_len(r)
    R11 <- qm[i1, i1, drop = FALSE]
    qty <- qr.qty(q, sw * y[rows])[i1]
    # The fitted value first, by the same operations whatever is asked of the
    # window, so that predict() at a training point returns fitted()'s value
    # to the last bit.  `cfit` is the regression point's own row.
    cP <- cfit[piv]
    t1 <- backsolve(R11, cP[i1], transpose = TRUE)
    est1 <- TRUE
    if (r < kk) {
      i2  <- (r + 1L):kk
      R12 <- qm[i1, i2, drop = FALSE]
      est1 <- max(abs(cP[i2] - drop(crossprod(R12, t1)))) <=
        tol * max(abs(cfit))
    }
    fi <- if (self) match(i, rows) else NA_integer_
    # A training point's own row is a row of its window, at the window's
    # largest weight, so it is in the row space whatever the test's rounding.
    if (!is.na(fi)) est1 <- TRUE
    if (est1) yhat[i] <- sum(t1 * qty)
    if (!need_q) {
      if (self && est1 && !is.na(fi)) {
        # The point's own row of the weighted design is sqrt(wt) * cfit, so
        # its leverage needs no Q.
        lev[i] <- wt[fi] * sum(t1^2)
        # Nor do its hat entries on the rows at its own location: the first
        # `r` columns of Q are the pivoted design times the inverse of R11.
        if (by_loc && shared[i]) {
          fg <- match(loc_rows[[loc$id[i]]], rows)
          if (!anyNA(fg))
            hat_loc[[i]] <- sw[fg] * drop(Z[fg, piv[i1], drop = FALSE] %*%
                                            backsolve(R11, t1))
        }
      }
      if (self) {
        if (r == kk) nest <- nest + 1L
        else {
          Tm  <- back(xf, mu)
          TP  <- Tm[, piv, drop = FALSE]
          T1  <- backsolve(R11, t(TP[, i1, drop = FALSE]), transpose = TRUE)
          rho <- TP[, i2, drop = FALSE] - crossprod(T1, R12)
          nest <- nest + (apply(abs(rho), 1L, max) <=
                            tol * apply(abs(Tm), 1L, max))
        }
      }
      next
    }
    # Column j of H: the weights function j's estimate puts on the window's
    # observations.  Its first column is the regression point's hat row.
    A    <- sw * qr.qy(q, diag(1, mk, r))
    srow <- drop(A %*% t1)
    if (est1) vfit[i] <- sum(srow^2)
    if (!full) next
    if (self && est1 && !is.na(fi)) {
      lev[i] <- srow[fi]
      qd[rows] <- qd[rows] + srow^2
      qd[i]    <- qd[i] + 1 - 2 * srow[fi]
      tssw[rows] <- tssw[rows] + wk * dyb[i]
      rssw[rows] <- rssw[rows] + wk * (y[i] - yhat[i])^2
    }
    Tm <- back(xf, mu)
    TP <- Tm[, piv, drop = FALSE]
    T1 <- backsolve(R11, t(TP[, i1, drop = FALSE]), transpose = TRUE)
    ok <- rep(TRUE, kk)
    if (r < kk) {
      rho <- TP[, i2, drop = FALSE] - crossprod(T1, R12)
      ok  <- apply(abs(rho), 1L, max) <= tol * apply(abs(Tm), 1L, max)
    }
    if (self) nest <- nest + ok
    B[i, ok] <- drop(crossprod(T1, qty))[ok]
    V[i, ok] <- colSums((A %*% T1)^2)[ok]
  }

  out <- list(fitted = yhat, rank = rk, n_window = nw, n_coef = kk,
              width = width)
  if (self) {
    out$leverage <- lev
    out$n_estimable <- as.integer(nest)
  } else out$nearest <- near
  if (by_loc) out$hat_loc <- hat_loc
  if (need_q) out$var_fit <- vfit
  if (full) {
    out$coef <- B
    out$coef_var <- V
    out$q_diag <- qd
    out$local_r2 <- (tssw - rssw) / tssw
  }
  out
}


#' Fitted values and leverages for a kernel that weights every row
#'
#' A Gaussian or an exponential kernel gives every training row a weight at
#' every regression point, so each window of \code{.gwr_local_solve()} is a
#' QR of the whole design: \code{n} factorisations of \code{n} rows at each
#' bandwidth a search scores, which made a default-bandwidth fit three to
#' five times slower than one through \code{GWmodel::bw.gwr()}.  What a
#' score needs, the fitted value and the leverage at each training point, is
#' computed here for many windows at once.
#'
#' With the predictors centred on their means over the data and divided by
#' their standard deviations, the weighted sums a window's regression needs
#' (\eqn{\sum w}, \eqn{\sum w x}, \eqn{\sum w x x^\top}, \eqn{\sum w y},
#' \eqn{\sum w x y}) are one matrix product of the kernel weights of a block
#' of regression points with the columns they are sums of.  From them, for
#' each window, the cross-product of the predictors centred at their weighted
#' mean, \eqn{C = \sum w x x^\top - (\sum w x)(\sum w x)^\top / \sum w}; the
#' slopes solve \eqn{C b = c_y}, the fitted value is
#' \eqn{\bar y_w + a^\top b} with \eqn{a} the point's predictors less their
#' weighted mean, and the leverage is \eqn{1/\sum w + a^\top C^{-1} a} (the
#' point's own weight is the largest, 1).  \eqn{C} is the cross-product of
#' the columns \code{.gwr_local_solve()} factors, so these are the same
#' numbers by another route.  The systems are solved by a Cholesky
#' factorisation of each window's correlation matrix, written out over the
#' block.
#'
#' \strong{Which windows are answered here.}  Forming \eqn{C} squares the
#' design's condition number, so a window is answered only where that costs
#' nothing that shows, and is left to the QR otherwise (\code{done} is
#' \code{FALSE}):
#' \itemize{
#'   \item every weight is defined, and each predictor's centred sum of
#'     squares is at least \code{1e-8} and at least \code{1e-4} of its
#'     uncentred one, so the subtraction that centres it keeps twelve digits;
#'   \item the determinant of the window's correlation matrix is at least
#'     \eqn{e \times 10^{-4}}.  The other eigenvalues multiply to less than
#'     \eqn{e} (they average at most \eqn{p/(p-1)}), so the smallest
#'     eigenvalue is at least \code{1e-4} and the Cholesky solve is good to
#'     about eleven digits.
#' }
#' Under both conditions the smallest singular value of the window's design
#' is at least \code{1e-6}, ten times the solve's tolerance, so the QR finds
#' the window of full rank as well: no rank decision is taken here.
#'
#' @param X,y,coords,bw,kernel,adaptive As for \code{.gwr_local_solve()};
#'   \code{X} has the intercept first and \code{kernel} is \code{"gaussian"}
#'   or \code{"exponential"}.
#' @param loc The layer's locations when each row's hat entries on the rows
#'   of its own location are wanted (see \code{.gwr_local_solve()}), else
#'   \code{NULL}.
#' @param cells Largest number of cells of a block's weight matrix.
#' @return A list, one value per training row: \code{done}, \code{fitted},
#'   \code{leverage}, \code{width}, \code{n_window} and, with \code{loc},
#'   \code{hat_loc}.  The values of a row that is not \code{done} are not to
#'   be used.
#' @keywords internal
#' @noRd
.gwr_score_dense <- function(X, y, coords, bw, kernel, adaptive, loc = NULL,
                             cells = 2^21) {
  n <- nrow(X); p <- ncol(X) - 1L
  jp <- seq_len(p) + 1L
  yhat <- rep(NA_real_, n); lev <- rep(NA_real_, n)
  width <- rep(NA_real_, n); nw <- integer(n); done <- logical(n)
  by_loc <- !is.null(loc) && loc$n < n
  if (by_loc) {
    loc_rows <- split(seq_len(n), loc$id)
    shared   <- loc$size[loc$id] > 1L
    hat_loc  <- vector("list", n)
  }
  if (p) {
    s <- apply(X[, jp, drop = FALSE], 2L, stats::sd)
    s[!is.finite(s) | s <= 0] <- 1
    Xc <- sweep(sweep(X[, jp, drop = FALSE], 2L,
                      colMeans(X[, jp, drop = FALSE])), 2L, s, "/")
    # The pairs (a, b), a <= b, of predictors, and where each sits.
    ia <- rep(seq_len(p), times = rev(seq_len(p)))
    ib <- unlist(lapply(seq_len(p), function(a) a:p))
    q  <- length(ia)
    at <- matrix(0L, p, p)
    at[cbind(ia, ib)] <- seq_len(q); at[cbind(ib, ia)] <- seq_len(q)
    dg <- at[cbind(seq_len(p), seq_len(p))]
    M  <- cbind(1, y, Xc, Xc * y, Xc[, ia, drop = FALSE] * Xc[, ib, drop = FALSE],
                deparse.level = 0L)
  } else {
    M <- cbind(1, y, deparse.level = 0L)
  }
  cx <- coords[, 1L]; cy <- coords[, 2L]
  B  <- max(1L, min(n, as.integer(cells %/% n)))
  gauss <- kernel == "gaussian"
  for (start in seq.int(1L, n, by = B)) {
    I  <- start:min(n, start + B - 1L)
    nb <- length(I)
    # The distances and weights of .gw_window(), a block of regression
    # points at a time.
    D <- sqrt(outer(cx[I], cx, "-")^2 + outer(cy[I], cy, "-")^2)
    h <- if (isTRUE(adaptive))
      vapply(seq_len(nb), function(r) .gw_kernel_width(D[r, ], bw, TRUE),
             numeric(1))
      else rep(as.numeric(bw), nb)
    W <- if (gauss) exp(D^2 / (-2 * h^2)) else exp(-D / h)
    width[I] <- h
    nw[I] <- as.integer(rowSums(W > 0, na.rm = TRUE))
    Mo   <- W %*% M
    S0   <- Mo[, 1L]
    ybar <- Mo[, 2L] / S0
    ok   <- is.finite(S0) & is.finite(h) & h > 0
    if (!p) {
      ok <- ok & is.finite(ybar)
      ok[is.na(ok)] <- FALSE
      yhat[I[ok]] <- ybar[ok]; lev[I[ok]] <- 1 / S0[ok]; done[I] <- ok
      if (by_loc) for (r in which(ok & shared[I])) {
        i <- I[r]
        hat_loc[[i]] <- rep(1 / S0[r], length(loc_rows[[loc$id[i]]]))
      }
      next
    }
    S1   <- Mo[, 2L + seq_len(p), drop = FALSE]
    Sxy  <- Mo[, 2L + p + seq_len(p), drop = FALSE]
    S2   <- Mo[, 2L + 2L * p + seq_len(q), drop = FALSE]
    xbar <- S1 / S0
    Cm   <- S2 - S1[, ia, drop = FALSE] * xbar[, ib, drop = FALSE]
    cyv  <- Sxy - S1 * ybar
    a    <- Xc[I, , drop = FALSE] - xbar
    cjj  <- Cm[, dg, drop = FALSE]
    thin <- rowSums(!(cjj >= pmax(1e-8, 1e-4 * S2[, dg, drop = FALSE])))
    ok   <- ok & !is.na(thin) & thin == 0L
    sdw  <- suppressWarnings(sqrt(cjj))
    # Cholesky of each window's correlation matrix, the windows of the block
    # side by side: L[, i, j] holds entry (i, j) of every window's factor.
    L    <- array(0, c(nb, p, p))
    dets <- rep(1, nb)
    for (j in seq_len(p)) {
      v <- rep(1, nb)
      if (j > 1L) v <- v - rowSums(matrix(L[, j, seq_len(j - 1L)], nb)^2)
      dets <- dets * v
      ljj  <- suppressWarnings(sqrt(v))
      L[, j, j] <- ljj
      if (j < p) for (i in (j + 1L):p) {
        e <- Cm[, at[i, j]] / (sdw[, i] * sdw[, j])
        if (j > 1L)
          e <- e - rowSums(matrix(L[, i, seq_len(j - 1L)], nb) *
                             matrix(L[, j, seq_len(j - 1L)], nb))
        L[, i, j] <- e / ljj
      }
    }
    ok <- ok & is.finite(dets) & dets >= exp(1) * 1e-4
    # Solve the windows' systems for one right-hand side each.
    solve_block <- function(R) {
      Z <- R
      for (j in seq_len(p)) {
        e <- R[, j]
        if (j > 1L)
          e <- e - rowSums(matrix(L[, j, seq_len(j - 1L)], nb) *
                             Z[, seq_len(j - 1L), drop = FALSE])
        Z[, j] <- e / L[, j, j]
      }
      sol <- Z
      for (j in rev(seq_len(p))) {
        e <- Z[, j]
        if (j < p)
          e <- e - rowSums(matrix(L[, (j + 1L):p, j], nb) *
                             sol[, (j + 1L):p, drop = FALSE])
        sol[, j] <- e / L[, j, j]
      }
      sol / sdw
    }
    slope <- solve_block(cyv / sdw)
    u     <- solve_block(a / sdw)
    yh    <- ybar + rowSums(a * slope)
    lv    <- 1 / S0 + rowSums(a * u)
    ok    <- ok & is.finite(yh) & is.finite(lv)
    ok[is.na(ok)] <- FALSE
    yhat[I[ok]] <- yh[ok]; lev[I[ok]] <- lv[ok]; done[I] <- ok
    # The hat entries of a row on the rows at its own location, which share
    # its window and its weight of 1.
    if (by_loc) for (r in which(ok & shared[I])) {
      i <- I[r]
      g <- loc_rows[[loc$id[i]]]
      hat_loc[[i]] <- 1 / S0[r] +
        drop((Xc[g, , drop = FALSE] - rep(xbar[r, ], each = length(g))) %*% u[r, ])
    }
  }
  out <- list(done = done, fitted = yhat, leverage = lev, width = width,
              n_window = nw)
  if (by_loc) out$hat_loc <- hat_loc
  out
}


#' The design matrix and response of a GWR, and the design at new rows
#'
#' Built as GWmodel builds them (\code{model.frame()} and
#' \code{model.matrix()} on the formula), with the intercept's column named
#' \code{Intercept}.  The design at new rows comes from the fitted terms, so
#' its columns are the coefficients' columns in the same order (a predictor
#' named twice gives one column, not two); \code{na.pass} keeps its rows
#' aligned with the new locations.
#'
#' @param formula The model formula.
#' @param data data.frame of the training rows.
#' @param newdata Optional data.frame of new rows, with the predictors.
#' @return A list: \code{X}, \code{y} and, with \code{newdata}, \code{X_new}.
#' @keywords internal
#' @noRd
.gwr_design <- function(formula, data, newdata = NULL) {
  mf <- stats::model.frame(formula, data = data, drop.unused.levels = TRUE)
  tt <- stats::terms(mf)
  X  <- stats::model.matrix(tt, mf)
  y  <- as.numeric(stats::model.response(mf))
  name_int <- function(m) {
    i <- match("(Intercept)", colnames(m))
    if (!is.na(i)) colnames(m)[i] <- "Intercept"
    m
  }
  out <- list(X = name_int(X), y = y)
  if (!is.null(newdata)) {
    tt_new <- stats::delete.response(tt)
    out$X_new <- name_int(stats::model.matrix(tt_new, stats::model.frame(
      tt_new, newdata, na.action = stats::na.pass,
      xlev = stats::.getXlevels(tt, mf))))
  }
  out
}


#' Assemble the fit's engine object in the layout GWmodel uses
#'
#' \code{fit$engine} was the object \code{GWmodel::gwr.basic()} returned, and
#' code that reads it (the package's own \code{coef()} and \code{fitted()},
#' and users' scripts) reads \code{$SDF} and \code{$GW.diagnostic}.  The
#' local solve's results are therefore laid out the same way: \code{SDF} is a
#' \code{SpatialPointsDataFrame} with the coefficient columns first, then
#' \code{y}, \code{yhat}, \code{residual}, \code{CV_Score},
#' \code{Stud_residual}, the standard errors (\code{<term>_SE}), the t values
#' (\code{<term>_TV}) and \code{Local_R2}; \code{GW.diagnostic} and
#' \code{GW.arguments} carry GWmodel's fields under GWmodel's names, and
#' \code{lm} is the global least-squares fit.  \code{CV_Score} is 0, as
#' \code{gwr.basic()} leaves it unless asked for it.  A coefficient a window
#' cannot estimate is \code{NA}, with its standard error and t value.
#'
#' @param sv What \code{.gwr_local_solve(what = "fit")} returned.
#' @param dg What \code{.gwr_diagnostics()} returned for it.
#' @param formula,sp_dat,X,y The model formula, the training layer as a
#'   Spatial object, the design matrix and the response.
#' @param bw,kernel,adaptive The bandwidth settings.
#' @return A list: \code{GW.arguments}, \code{GW.diagnostic}, \code{lm},
#'   \code{SDF}, \code{fitted} (the fitted values, \code{NA} where not
#'   estimable) and \code{solver} (the method and its tolerance).
#' @keywords internal
#' @noRd
.gwr_engine <- function(sv, dg, formula, sp_dat, X, y, bw, kernel, adaptive) {
  n     <- length(y)
  betas <- sv$coef
  yhat  <- sv$fitted
  residual <- y - yhat
  stud  <- residual / sqrt(dg$sigma2 * sv$q_diag)
  se    <- sqrt(dg$sigma2 * sv$coef_var)
  tv    <- betas / se
  gwres <- data.frame(betas, y, yhat, residual, numeric(n), stud, se, tv,
                      sv$local_r2)
  nm <- colnames(betas)
  colnames(gwres) <- c(nm, "y", "yhat", "residual", "CV_Score",
                       "Stud_residual", paste(nm, "SE", sep = "_"),
                       paste(nm, "TV", sep = "_"), "Local_R2")
  # The training rows' own names, as GWmodel's table carries them.
  rownames(gwres) <- rownames(sp_dat@data)
  crd <- sp::coordinates(sp_dat)
  rownames(crd) <- rownames(gwres)
  sdf <- sp::SpatialPointsDataFrame(coords = crd, data = gwres,
                                    proj4string = sp_dat@proj4string,
                                    match.ID = FALSE)
  lms <- stats::lm(formula, data = sp_dat@data)
  lms$x <- X
  lms$y <- y
  list(
    GW.arguments = list(formula = formula, rp.given = FALSE, hatmatrix = TRUE,
                        bw = bw, kernel = kernel, adaptive = adaptive, p = 2,
                        theta = 0, longlat = FALSE, DM.given = FALSE,
                        F123.test = FALSE),
    GW.diagnostic = list(RSS.gw = dg$RSS, AIC = dg$AIC, AICc = dg$AICc,
                         enp = dg$enp, edf = dg$edf, gw.R2 = dg$R2,
                         gwR2.adj = dg$R2_adj, BIC = dg$BIC),
    lm = lms,
    SDF = sdf,
    fitted = yhat,
    solver = list(method = "pivoted Householder QR per window (LAPACK dgeqp3)",
                  tolerance = .gwr_rank_tol())
  )
}


#' GWmodel's global diagnostics, from the local solve
#'
#' The formulas of GWmodel's \code{gwr_diag1()}, on the trace of the hat
#' matrix and of \eqn{S^\top S} that \code{.gwr_local_solve()} returns:
#' \deqn{AIC = n \log(RSS/n) + n \log(2\pi) + n + tr S,}
#' \deqn{AICc = n \log(RSS/n) + n \log(2\pi) + n (n + tr S) / (n - 2 - tr S),}
#' \deqn{BIC = n \log(RSS/n) + n \log(2\pi) + \log(n)\, tr S,}
#' effective number of parameters \eqn{2\,tr S - tr S^\top S}, effective
#' degrees of freedom \eqn{n - 2\,tr S + tr S^\top S}, and the residual
#' variance \eqn{RSS / edf}.  AICc is returned as the formula gives it, also
#' outside its domain (\eqn{tr S \ge n - 2}); \code{.gwr_aicc_undefined()} is
#' the guard.  Every value is \code{NA} when a fitted value is: the fit then
#' covers part of the data only and has no residual sum of squares.
#'
#' @param y Response.
#' @param fitted,leverage Fitted values and leverages, one per observation.
#' @param var_fit Sum of squares of each hat row, or \code{NULL} when only
#'   AIC, AICc, BIC and the trace are wanted.
#' @return A named list: \code{RSS}, \code{tr_S}, \code{tr_StS}, \code{AIC},
#'   \code{AICc}, \code{BIC}, \code{enp}, \code{edf}, \code{sigma2},
#'   \code{R2}, \code{R2_adj}.
#' @keywords internal
#' @noRd
.gwr_diagnostics <- function(y, fitted, leverage, var_fit = NULL) {
  n <- length(y)
  na <- list(RSS = NA_real_, tr_S = NA_real_, tr_StS = NA_real_,
             AIC = NA_real_, AICc = NA_real_, BIC = NA_real_, enp = NA_real_,
             edf = NA_real_, sigma2 = NA_real_, R2 = NA_real_,
             R2_adj = NA_real_)
  if (!n || anyNA(fitted) || anyNA(leverage)) return(na)
  rss  <- sum((y - fitted)^2)
  tr_s <- sum(leverage)
  base <- n * log(rss / n) + n * log(2 * pi)
  out <- na
  out$RSS  <- rss
  out$tr_S <- tr_s
  out$AIC  <- base + n + tr_s
  out$AICc <- base + n * ((n + tr_s) / (n - 2 - tr_s))
  out$BIC  <- base + log(n) * tr_s
  if (!is.null(var_fit) && !anyNA(var_fit)) {
    tr_sts <- sum(var_fit)
    edf <- n - 2 * tr_s + tr_sts
    out$tr_StS <- tr_sts
    out$enp <- 2 * tr_s - tr_sts
    out$edf <- edf
    out$sigma2 <- if (isTRUE(edf > 0)) rss / edf else NA_real_
    yss <- sum((y - mean(y))^2)
    out$R2 <- 1 - rss / yss
    out$R2_adj <- 1 - (1 - out$R2) * (n - 1) / (edf - 1)
  }
  out
}


#' The distinct locations of a layer
#'
#' Two rows are at one location when both coordinates are equal, which is
#' when the kernel gives them one weight at every regression point.  Compared
#' as numbers, not as printed text.
#'
#' @param coords Two-column coordinate matrix.
#' @return A list: \code{id} (the location of each row, numbered from 1 in
#'   order of first appearance), \code{n} (the number of locations) and
#'   \code{size} (the rows at each).
#' @keywords internal
#' @noRd
.gwr_locations <- function(coords) {
  n <- nrow(coords)
  if (n < 2L) return(list(id = seq_len(n), n = n, size = rep(1L, n)))
  cx <- coords[, 1L]; cy <- coords[, 2L]
  o  <- order(cx, cy)
  sx <- cx[o]; sy <- cy[o]
  first <- c(TRUE, sx[-1L] != sx[-n] | sy[-1L] != sy[-n])
  first[is.na(first)] <- TRUE
  id <- integer(n)
  id[o] <- cumsum(first)
  id <- match(id, unique(id))
  list(id = id, n = max(id), size = tabulate(id))
}


#' The number of rows that differ in location or in a predictor
#'
#' Rows at one location with one set of predictor values receive one fitted
#' value at every bandwidth, so a GWR of \code{n} such rows is a GWR of the
#' means of the groups they form.  Its hat-matrix trace cannot exceed the
#' number of groups, and at that number the fit reproduces every group's mean.
#'
#' @param X Design matrix.
#' @param id Location of each row (\code{.gwr_locations()$id}).
#' @return One integer, \code{nrow(X)} when no two rows agree in both.
#' @keywords internal
#' @noRd
.gwr_n_distinct <- function(X, id) {
  n <- nrow(X)
  if (n < 2L) return(n)
  o  <- do.call(order, c(list(id), lapply(seq_len(ncol(X)), function(j) X[, j])))
  Xo <- X[o, , drop = FALSE]
  io <- id[o]
  same <- io[-1L] == io[-n] &
    rowSums(Xo[-1L, , drop = FALSE] != Xo[-n, , drop = FALSE]) == 0L
  as.integer(n - sum(same, na.rm = TRUE))
}


#' Leave-one-location-out residuals, from the hat matrix
#'
#' The residual of each observation when its location is predicted by the
#' local regression there without any observation of that location, the
#' kernel keeping the width it has in the fit.  The rows at one location
#' share one window and one weight, so the local regression there is one
#' weighted least squares and the usual deletion identity holds for the
#' group: with \eqn{H_{gg}} the block of the hat matrix on the location's
#' rows and \eqn{e_g} their residuals, the residuals without the location
#' are \eqn{(I - H_{gg})^{-1} e_g}.  For a location with one row that is
#' \eqn{e_i / (1 - h_{ii})}, the leave-one-out residual.
#'
#' \eqn{I - H_{gg}} is singular exactly when the other locations in the
#' window cannot place one of the location's rows: nothing predicts it but
#' its own observations.  The residuals of such a location are \code{NA}.
#' Singular is judged at the solve's tolerance, as a leverage of 1 is.
#'
#' @param resid Residuals of the fit, one per row.
#' @param leverage Leverages, one per row.
#' @param hat_loc What \code{.gwr_local_solve(loc = )} returned as
#'   \code{hat_loc}.
#' @param loc The layer's locations (\code{.gwr_locations()}).
#' @param tol The solve's tolerance.
#' @return Numeric vector, one per row, \code{NA} where no prediction without
#'   the location exists.
#' @keywords internal
#' @noRd
.gwr_lolo_residuals <- function(resid, leverage, hat_loc, loc,
                                tol = .gwr_rank_tol()) {
  n   <- length(resid)
  out <- rep(NA_real_, n)
  one <- loc$size[loc$id] == 1L
  ok  <- one & is.finite(resid) & is.finite(leverage) & leverage < 1 - tol
  out[ok] <- resid[ok] / (1 - leverage[ok])
  if (all(one)) return(out)
  rows <- split(seq_len(n), loc$id)
  for (g in which(loc$size > 1L)) {
    r  <- rows[[g]]
    hb <- hat_loc[r]
    if (anyNA(resid[r]) || any(vapply(hb, is.null, logical(1)))) next
    M  <- diag(1, length(r)) - do.call(rbind, hb)
    sv <- svd(M, nu = 0L, nv = 0L)$d
    if (!all(is.finite(sv)) || min(sv) < tol) next
    out[r] <- solve(M, resid[r])
  }
  out
}


#' The bandwidth criterion at one bandwidth
#'
#' One pass of \code{.gwr_local_solve(what = "score")} and the criterion from
#' its fitted values and leverages.  \code{"AICc"} is GWmodel's corrected AIC
#' (see \code{.gwr_diagnostics()}).  \code{"CV"} is the leave-one-location-out
#' sum of squares (\code{.gwr_lolo_residuals()}).  On a layer without repeat
#' locations that is the leave-one-out sum
#' \eqn{\sum (e_i / (1 - h_{ii}))^2}, which equals refitting each window
#' without its own observation, as \code{GWmodel::bw.gwr()} does.  Where
#' locations repeat, leaving one row out leaves its twins in, and the score
#' then rewards a window that holds little but the location itself.
#'
#' A bandwidth is scored \code{Inf}, so that no search settles on it, and
#' \code{why} names the first of these that applies (a fitted value, then
#' the residual sum of squares, then the criterion's own domain, then the
#' terms):
#' \describe{
#'   \item{\code{"fitted"}}{a fitted value is not estimable;}
#'   \item{\code{"rss"}}{the residual sum of squares is 0;}
#'   \item{\code{"term"}}{a term is estimable at no location, so the model
#'     has no coefficient for it anywhere;}
#'   \item{\code{"df"}}{AICc is outside its domain (\eqn{tr S \ge n - 2});}
#'   \item{\code{"df_distinct"}}{\eqn{tr S} is not below the number of rows
#'     that differ in location or in a predictor, minus 2
#'     (\code{.gwr_n_distinct()}): the fit all but reproduces the mean of
#'     each group of rows that share both, and the AICc of that fit to the
#'     means is outside its domain;}
#'   \item{\code{"leave_out"}}{for \code{"CV"}, a location rests on its own
#'     observations (for a single row, leverage 1 to the solve's tolerance),
#'     so its prediction without them does not exist.}
#' }
#'
#' @param X,y,coords,bw,kernel,adaptive As for \code{.gwr_local_solve()}.
#' @param approach \code{"AICc"} or \code{"CV"}.
#' @param loc The layer's locations (\code{.gwr_locations()}); worked out
#'   when \code{NULL}.
#' @param n_distinct \code{.gwr_n_distinct()} of the design; worked out when
#'   \code{NULL}.
#' @return A list: \code{value} (the criterion, possibly \code{Inf}),
#'   \code{why} (\code{""} when the criterion is defined), \code{RSS},
#'   \code{tr_S}, \code{AIC}, \code{AICc} (as the formulas give them,
#'   \code{NA} where a fitted value is not estimable), \code{n_deficient}
#'   (windows of deficient rank) and \code{n_distinct}.
#' @keywords internal
#' @noRd
.gwr_criterion <- function(X, y, coords, bw, kernel, adaptive,
                           approach = "AICc", loc = NULL, n_distinct = NULL) {
  n  <- length(y)
  cv <- identical(approach, "CV")
  if (is.null(loc)) loc <- .gwr_locations(coords)
  repeats <- loc$n < n
  if (is.null(n_distinct))
    n_distinct <- if (repeats) .gwr_n_distinct(X, loc$id) else n
  sv <- .gwr_local_solve(X, y, coords, bw, kernel, adaptive, what = "score",
                         loc = if (cv && repeats) loc else NULL)
  dg <- .gwr_diagnostics(y, sv$fitted, sv$leverage)
  value <- Inf
  why <- ""
  # The first reason that applies is the one reported; a term that no window
  # estimates comes last, since a fit that interpolates has none either.
  no_term <- any(sv$n_estimable == 0L)
  if (!is.finite(dg$RSS)) {
    why <- "fitted"
  } else if (!(dg$RSS > 0)) {
    why <- "rss"
  } else if (cv) {
    loo <- if (repeats)
      .gwr_lolo_residuals(y - sv$fitted, sv$leverage, sv$hat_loc, loc)
    else if (all(sv$leverage < 1 - .gwr_rank_tol()))
      (y - sv$fitted) / (1 - sv$leverage)
    else NA_real_
    if (anyNA(loo)) why <- "leave_out"
    else if (no_term) why <- "term"
    else value <- sum(loo^2)
  } else if (!(dg$tr_S < n - 2) || !is.finite(dg$AICc)) {
    why <- "df"
  } else if (!(dg$tr_S < n_distinct - 2)) {
    why <- "df_distinct"
  } else if (no_term) {
    why <- "term"
  } else {
    value <- dg$AICc
  }
  list(value = value, why = why, RSS = dg$RSS, tr_S = dg$tr_S, AIC = dg$AIC,
       AICc = dg$AICc, n_deficient = sum(sv$rank < sv$n_coef),
       n_distinct = n_distinct)
}


#' The range a bandwidth search covers
#'
#' Adaptive: from the smallest count that leaves every window a residual
#' degree of freedom (\code{.gwr_min_adaptive_bw()}) to every observation.
#' On a layer whose locations repeat the count has to buy locations, not
#' rows: the lower end is then the smallest count at which every window
#' holds one location more than the model has parameters (the rows of the
#' nearest such locations, and one more row for the bisquare and tricube
#' kernels, which give the farthest row weight 0).
#'
#' Fixed: up to the largest distance between two observations.  The lower end
#' is the distance at which half of the windows hold one location more than
#' the model has parameters (the median, over the observations, of the
#' distance to the location that completes such a window); a third of that
#' for the Gaussian and exponential kernels, which weight points beyond their
#' bandwidth.  It is the median and not the largest of those distances
#' because one isolated observation would otherwise set the floor for the
#' whole layer: its own window may be short of points at a bandwidth that
#' suits the rest, the local solve returns what that window can estimate, and
#' the criterion counts it.  Below the lower end most windows interpolate and
#' the criterion is undefined.
#'
#' \code{GWmodel::bw.gwr()} searches adaptive counts from 20 up, whatever the
#' model, and fixed distances from a 5000th of the largest, where most
#' windows hold their own point alone.
#'
#' @param coords Two-column coordinate matrix.
#' @param n_params Parameters of the local model, the intercept included.
#' @param kernel,adaptive Kernel name; whether the bandwidth counts
#'   neighbours.
#' @param loc The layer's locations (\code{.gwr_locations()}).
#' @return \code{c(lower, upper)}; equal when there is one candidate, and
#'   \code{NA} when every observation is at one location.
#' @keywords internal
#' @noRd
.gwr_bw_range <- function(coords, n_params, kernel, adaptive,
                          loc = .gwr_locations(coords)) {
  n  <- nrow(coords)
  nl <- loc$n
  repeats <- nl < n
  if (isTRUE(adaptive) && !repeats) {
    lo <- min(.gwr_min_adaptive_bw(n_params, kernel), n)
    return(as.numeric(c(lo, n)))
  }
  # One row per location, in the order the locations are numbered (which is
  # the rows' own order when none repeats).
  at <- !duplicated(loc$id)
  ux <- coords[at, 1L]; uy <- coords[at, 2L]
  if (isTRUE(adaptive)) {
    need <- min(as.integer(n_params) + 1L, nl)
    edge <- kernel %in% c("bisquare", "tricube")
    lo <- 0
    for (g in seq_len(nl)) {
      d  <- sqrt((ux - ux[g])^2 + (uy - uy[g])^2)
      dn <- sort.int(d, partial = need)[need]
      # A kernel with an edge weights the rows strictly inside the distance
      # of the count-th nearest row; the others weight that row too.
      lo <- max(lo, if (edge) sum(loc$size[d <= dn]) + 1L
                    else sum(loc$size[d < dn]) + 1L)
    }
    return(as.numeric(c(min(lo, n), n)))
  }
  # The (n_params + 2)-th smallest distance from a location, itself included
  # at 0: with a kernel that gives its edge weight 0 the window then holds
  # n_params + 1 weighted locations.
  kth <- min(as.integer(n_params) + 2L, nl)
  dk <- numeric(nl); hi <- 0
  for (g in seq_len(nl)) {
    d <- sqrt((ux - ux[g])^2 + (uy - uy[g])^2)
    dk[g] <- sort.int(d, partial = kth)[kth]
    hi <- max(hi, max(d))
  }
  if (!is.finite(hi) || hi <= 0) return(c(NA_real_, NA_real_))
  lo <- stats::median(dk[loc$id])
  if (kernel %in% c("gaussian", "exponential")) lo <- lo / 3
  # Every observation at one of very few locations: no distance to go by.
  if (!is.finite(lo) || lo <= 0) lo <- hi / 5000
  c(min(lo, hi), hi)
}


#' Search for the bandwidth that minimises AICc or the leave-one-out score
#'
#' \code{GWmodel::bw.gwr()} runs a golden-section search over the whole
#' range, which is right when the criterion has one minimum.  It has several
#' whenever the coefficients vary at a range only a few dozen points resolve:
#' the search then throws away the lower part of the range at its first
#' step and returns a near-global bandwidth labelled as selected.  This
#' search scores a coarse grid first and refines inside the bracket round the
#' grid's minimum:
#' \enumerate{
#'   \item \code{n_grid} bandwidths spaced evenly in the logarithm between the
#'     two ends of \code{.gwr_bw_range()} (whole numbers when adaptive);
#'   \item a golden-section search between the grid points either side of the
#'     best one, to a bracket of two neighbours (adaptive) or a thousandth of
#'     the bandwidth (fixed), at most \code{max_refine} steps; for an
#'     adaptive bandwidth the counts within two of the best are then scored
#'     too, while that finds a better one (three times at most), because the
#'     criterion is jagged from one count to the next;
#'   \item the best of every bandwidth scored, grid included, is returned.
#' }
#' Each bandwidth is scored once (\code{.gwr_criterion()}, through the local
#' QR solve), so a candidate with windows of deficient rank is scored on its
#' estimable fitted values and not as \code{Inf} or as rounding error.
#'
#' \strong{Layers whose locations repeat.}  When two or more rows share a
#' location, \code{"AICc"} is not used and the search minimises the
#' leave-one-location-out score instead (\code{approach} comes back as
#' \code{"CV"}, \code{requested} as \code{"AICc"}, \code{by_location} as
#' \code{TRUE}).  AICc counts the rows at one location as independent
#' observations.  When a location has an effect of its own, a window that
#' holds little but the location reproduces its mean, the residual sum of
#' squares falls to the spread within locations, and AICc is lowest at the
#' bandwidths where the windows interpolate: on 100 sites with 3 rows each
#' it was 58 at 5 to 10 neighbours, 90 at 13 and 221 at 58, while the error
#' of prediction at new sites was 46 at 10 neighbours, 1.35 at 13 and 0.15
#' at 150.  A floor on the bandwidth does not help, since the criterion goes
#' on falling towards it.  Predicting each location from the others asks the
#' question a bandwidth is chosen for.  \code{"CV"} leaves out whole locations for the same reason: with a
#' twin left in, the leave-one-out prediction of a row is its twin.
#'
#' What the grid shows is kept: the number of local minima on it, the range
#' of bandwidths whose criterion is within 2 of the minimum (1 percent of it
#' for \code{"CV"}), whether the widest bandwidth of the range is one of
#' them, and whether the minimum sits at either end of the range.
#'
#' @param X,y,coords,kernel,adaptive As for \code{.gwr_local_solve()}.
#' @param approach \code{"AICc"} or \code{"CV"}.
#' @param n_grid Number of grid bandwidths.
#' @param max_refine Largest number of golden-section steps.
#' @return A list: \code{bandwidth} (\code{NA} when the range is empty),
#'   \code{criterion}, \code{approach} (the criterion minimised),
#'   \code{requested} (the one asked for), \code{by_location} (\code{TRUE}
#'   when locations repeat, so that the score leaves out whole locations),
#'   \code{n_locations}, \code{range}, \code{curve} (a
#'   data.frame of every bandwidth scored: \code{bandwidth}, \code{criterion},
#'   \code{tr_S}, \code{n_deficient} and \code{stage}, \code{"grid"} or
#'   \code{"refine"}), \code{n_evaluations}, \code{n_local_minima},
#'   \code{multimodal}, \code{within} (the smallest and largest scored
#'   bandwidth within the margin of the minimum), \code{flat} (at least half
#'   of the grid is within the margin), \code{near_global} (the widest
#'   bandwidth of the range is within the margin, so the criterion does not
#'   tell the selected fit from a global one), \code{at_upper},
#'   \code{at_lower}, \code{none_defined} (no bandwidth has a finite
#'   criterion; the widest is then returned) and \code{undefined} (how many
#'   of the bandwidths scored have no criterion, by the reason
#'   \code{.gwr_criterion()} gives).
#' @keywords internal
#' @noRd
.gwr_bw_search <- function(X, y, coords, kernel, adaptive, approach = "AICc",
                           n_grid = 16L, max_refine = 12L) {
  adaptive <- isTRUE(adaptive)
  loc <- .gwr_locations(coords)
  by_location <- loc$n < nrow(coords)
  requested <- approach
  # AICc counts the rows at one location as independent; see above.
  if (by_location && identical(approach, "AICc")) approach <- "CV"
  n_distinct <- if (by_location) .gwr_n_distinct(X, loc$id) else nrow(X)
  rng <- .gwr_bw_range(coords, ncol(X), kernel, adaptive, loc)
  out <- list(bandwidth = NA_real_, criterion = NA_real_, approach = approach,
              requested = requested, by_location = by_location,
              n_locations = loc$n,
              range = rng, curve = NULL, n_evaluations = 0L,
              n_local_minima = NA_integer_, multimodal = NA, within = rng,
              flat = NA, near_global = NA, at_upper = NA, at_lower = NA,
              none_defined = NA, undefined = integer(0))
  if (anyNA(rng)) return(out)
  lo <- rng[1L]; hi <- rng[2L]

  bws <- numeric(0); val <- numeric(0); trs <- numeric(0)
  ndef <- integer(0); stage <- character(0); why <- character(0)
  score <- function(b, st) {
    hit <- match(b, bws)
    if (!is.na(hit)) return(val[hit])
    cr <- tryCatch(.gwr_criterion(X, y, coords, b, kernel, adaptive, approach,
                                  loc = loc, n_distinct = n_distinct),
                   error = function(e) NULL)
    v <- if (is.null(cr) || !is.finite(cr$value)) Inf else cr$value
    bws <<- c(bws, b); val <<- c(val, v)
    trs <<- c(trs, if (is.null(cr)) NA_real_ else cr$tr_S)
    ndef <<- c(ndef, if (is.null(cr)) NA_integer_ else cr$n_deficient)
    stage <<- c(stage, st)
    why <<- c(why, if (is.null(cr)) "error" else cr$why)
    v
  }

  grid <- if (hi > lo) exp(seq(log(lo), log(hi), length.out = max(2L, n_grid)))
          else lo
  if (adaptive) grid <- unique(pmin(pmax(round(grid), lo), hi))
  grid[1L] <- lo; grid[length(grid)] <- hi
  gval <- vapply(grid, score, numeric(1), st = "grid")

  if (!any(is.finite(gval))) {
    # Nothing to choose between: the widest window is the only one that can
    # have residual degrees of freedom, and the caller says so.
    out$bandwidth <- hi
    out$none_defined <- TRUE
  } else {
    out$none_defined <- FALSE
    j <- which.min(gval)
    a <- grid[max(j - 1L, 1L)]; b <- grid[min(j + 1L, length(grid))]
    gr <- (sqrt(5) - 1) / 2
    snap <- function(z) if (adaptive) round(z) else z
    wide <- function() if (adaptive) (b - a) > 2 else
      (b - a) > 1e-3 * max(a, .Machine$double.eps)
    if (wide()) {
      x1 <- snap(b - gr * (b - a)); x2 <- snap(a + gr * (b - a))
      f1 <- score(x1, "refine"); f2 <- score(x2, "refine")
      steps <- 0L
      while (wide() && x1 < x2 && steps < max_refine) {
        # Ties, two Inf among them, close in on the wider bandwidth.
        if (f1 < f2) {
          b <- x2; x2 <- x1; f2 <- f1
          x1 <- snap(b - gr * (b - a)); f1 <- score(x1, "refine")
        } else {
          a <- x1; x1 <- x2; f1 <- f2
          x2 <- snap(a + gr * (b - a)); f2 <- score(x2, "refine")
        }
        steps <- steps + 1L
      }
    }
    if (adaptive) {
      # Whole numbers: the few left inside the bracket are each scored.
      if ((b - a) >= 1 && (b - a) <= 4)
        for (z in seq(a, b)) score(z, "refine")
      # Over whole numbers the criterion is jagged from one count to the
      # next, so the counts within two of the best are scored as well, again
      # while that finds a better one, three times at most.
      for (pass in 1:3) {
        at   <- bws[which.min(val)]
        cand <- setdiff(seq(max(lo, at - 2), min(hi, at + 2)), bws)
        if (!length(cand)) break
        was <- min(val)
        for (z in cand) score(z, "refine")
        if (!(min(val) < was)) break
      }
    }
    best <- which.min(val)
    out$bandwidth <- bws[best]
    out$criterion <- val[best]
  }

  ord <- order(bws)
  out$curve <- data.frame(bandwidth = bws[ord], criterion = val[ord],
                          tr_S = trs[ord], n_deficient = ndef[ord],
                          stage = stage[ord], stringsAsFactors = FALSE)
  out$n_evaluations <- length(bws)
  # Why the bandwidths without a criterion have none, for the caller's
  # message when that is all of them.
  bad <- why[nzchar(why)]
  out$undefined <- vapply(split(bad, bad), length, integer(1))
  fin <- is.finite(gval)
  # Local minima on the grid: below both neighbours (an end point has one).
  ng <- length(grid)
  left  <- c(Inf, gval[-ng]); right <- c(gval[-1L], Inf)
  out$n_local_minima <- sum(fin & gval < left & gval < right)
  out$multimodal <- out$n_local_minima > 1L
  if (isFALSE(out$none_defined)) {
    margin <- if (identical(approach, "CV")) 0.01 * abs(out$criterion) else 2
    near <- is.finite(val) & val <= out$criterion + margin
    out$within <- range(bws[near])
    out$flat <- sum(fin & gval <= out$criterion + margin) >= ng / 2 && ng >= 4L
    out$near_global <- isTRUE(gval[ng] <= out$criterion + margin)
    out$at_upper <- isTRUE(out$bandwidth >= hi)
    out$at_lower <- isTRUE(out$bandwidth <= lo) && hi > lo
  }
  out
}


#' The name of the criterion a bandwidth search minimised
#'
#' @param bs What \code{.gwr_bw_search()} returned.
#' @param score \code{TRUE} for the noun a bandwidth "has" (a score),
#'   \code{FALSE} for the method it was selected by.
#' @return \code{"AICc"}, or the name of the cross-validation.
#' @keywords internal
#' @noRd
.gwr_bw_criterion_name <- function(bs, score = FALSE) {
  if (!identical(bs$approach, "CV")) return("AICc")
  paste0(if (isTRUE(bs$by_location)) "leave-one-location-out"
         else "leave-one-out",
         if (score) " score" else " cross-validation")
}


#' Why the bandwidths a search scored have no criterion
#'
#' For the warning raised when that is every one of them.
#'
#' @param bs What \code{.gwr_bw_search()} returned.
#' @param n Number of observations.
#' @return The reasons, joined by "or", as the clause that follows "at each
#'   of the bandwidths tried".
#' @keywords internal
#' @noRd
.gwr_bw_undefined_text <- function(bs, n) {
  text <- c(
    df          = sprintf(paste0("the local regressions leave too few ",
                                 "residual degrees of freedom (tr(S) is not ",
                                 "below n - 2 = %d)"), n - 2L),
    df_distinct = paste0("the local regressions all but reproduce the mean ",
                         "of each group of observations that share a ",
                         "location and predictor values"),
    fitted      = "a fitted value cannot be estimated",
    term        = "a term can be estimated at no location",
    leave_out   = if (isTRUE(bs$by_location))
      "a location cannot be predicted from the others"
      else "an observation cannot be predicted from the others",
    rss         = "the residual sum of squares is 0",
    error       = "the local regressions could not be solved")
  hit <- intersect(names(text), names(bs$undefined))
  if (!length(hit)) hit <- c("df", "fitted")
  paste(text[hit], collapse = ", or ")
}


#' One sentence on what a bandwidth search's grid showed
#'
#' For the line \code{fit_gwr_model()} and \code{gwr_model_selection()}
#' log once per search.  Empty when there is nothing to remark on.
#'
#' @param bs What \code{.gwr_bw_search()} returned.
#' @param adaptive Whether the bandwidth counts neighbours.
#' @return A character string, possibly \code{""}.
#' @keywords internal
#' @noRd
.gwr_bw_search_note <- function(bs, adaptive) {
  if (is.null(bs) || !isFALSE(bs$none_defined)) return("")
  unit <- if (isTRUE(adaptive)) " neighbours" else ""
  bits <- character(0)
  if (isTRUE(bs$at_upper))
    bits <- c(bits, sprintf(paste0(
      "the minimum is at the top of the range (%s%s), so the fit is close to ",
      "one global regression"), format(signif(bs$bandwidth, 6)), unit))
  else if (isTRUE(bs$near_global))
    bits <- c(bits, sprintf(paste0(
      "the widest bandwidth of the range (%s%s) scores within %s of the ",
      "minimum, so the criterion does not tell this fit from one global ",
      "regression"), format(signif(bs$range[2L], 6)), unit,
      if (identical(bs$approach, "CV")) "1 percent" else "2"))
  if (isTRUE(bs$at_lower))
    bits <- c(bits, sprintf(paste0(
      "the minimum is at the bottom of the range (%s%s), the smallest ",
      "bandwidth the model's parameters allow"),
      format(signif(bs$bandwidth, 6)), unit))
  if (isTRUE(bs$multimodal))
    bits <- c(bits, sprintf(paste0(
      "the criterion has %d local minima on the %d-point grid, so the ",
      "bandwidth is the best of those scored and not the one minimum of a ",
      "smooth curve"), bs$n_local_minima, sum(bs$curve$stage == "grid")))
  if (isTRUE(bs$flat))
    bits <- c(bits, sprintf(paste0(
      "the criterion is flat: bandwidths from %s to %s%s are within %s of ",
      "the minimum"), format(signif(bs$within[1L], 6)),
      format(signif(bs$within[2L], 6)), unit,
      if (identical(bs$approach, "CV")) "1 percent" else "2"))
  if (!length(bits)) return("")
  paste0(paste(bits, collapse = "; "), ".")
}


#' How far the widest kernel of a fit reaches
#'
#' What \code{predict.gwr_fit()} compares a new location's distance to its
#' nearest training point with.  The widest kernel is the bandwidth itself
#' when it is fixed; for an adaptive one it is the largest distance from a
#' training point to its \code{bw}-th nearest training point, which
#' \code{fit_gwr_model()} keeps as \code{info$max_kernel_width} and which is
#' worked out from the training coordinates for a fit that does not carry it
#' (one put together by hand).  A bisquare, tricube or boxcar kernel reaches
#' as far as its width and no farther.  A Gaussian or exponential kernel has
#' no edge; its reach is taken as three widths, where a point's weight is
#' 0.011 (Gaussian) or 0.05 (exponential) of the weight at the centre.
#'
#' @param object The \code{gwr_fit}.
#' @param coords Two-column matrix of the training coordinates.
#' @param bw,adaptive,kernel The fit's bandwidth settings.
#' @return A list: \code{width} (the widest kernel) and \code{reach}.
#' @keywords internal
#' @noRd
.gwr_kernel_reach <- function(object, coords, bw, adaptive, kernel) {
  width <- if (!isTRUE(adaptive)) as.numeric(bw) else {
    kept <- suppressWarnings(as.numeric(object$info$max_kernel_width %||% NA_real_))
    if (length(kept) == 1L && is.finite(kept)) kept else {
      cx <- coords[, 1L]; cy <- coords[, 2L]
      out <- 0
      for (i in seq_along(cx))
        out <- max(out, .gw_kernel_width(sqrt((cx - cx[i])^2 + (cy - cy[i])^2),
                                         bw, TRUE))
      out
    }
  }
  list(width = width,
       reach = if (kernel %in% c("gaussian", "exponential")) 3 * width else width)
}


#' Raise a warning with a condition class, and log it
#'
#' \code{.warn_and_log()} with a class, for the GWR conditions a caller may
#' want to catch by name (\code{tryCatch(spatialkit_gwr_rank_deficient = )},
#' \code{withCallingHandlers()}) and not by matching message text.
#'
#' @param class The condition class.
#' @param fmt,... As for \code{sprintf()}.
#' @param data Named list of fields to put on the condition, so that a
#'   handler reads the counts and does not parse the message.
#' @param tallied \code{TRUE} for a warning that
#'   \code{.gwr_tally_predict_warnings()} collects: while it is collecting,
#'   the line is not logged here, since the total is logged once at the end.
#' @return The message, invisibly.
#' @keywords internal
#' @noRd
.gwr_warn <- function(class, fmt, ..., data = NULL, tallied = FALSE) {
  # Logged and raised by the package's one classed-warning helper, except
  # while a collector is totalling: then it is raised for the collector only.
  if (!(isTRUE(tallied) && isTRUE(.gwr_predict_state$tally > 0L)))
    return(.warn_classed(class, fmt, ..., data = data %||% list()))
  msg <- sprintf(fmt, ...)
  warning(do.call(warningCondition,
                  c(list(message = msg, class = class, call = NULL), data)))
  invisible(msg)
}


# Whether .gwr_tally_predict_warnings() is collecting (its nesting depth).
.gwr_predict_state <- new.env(parent = emptyenv())

# The two warnings in which predict() on a GWR counts locations, and which
# .gwr_tally_predict_warnings() totals.
.gwr_counted_classes <- c("spatialkit_gwr_extrapolation",
                          "spatialkit_gwr_not_estimable")


#' The text of the two counted warnings of predict.gwr_fit()
#'
#' Shared by \code{predict.gwr_fit()}, which raises them for one call, and
#' \code{.gwr_tally_predict_warnings()}, which raises them once for many.
#'
#' @param caller Name of the function the message is from.
#' @param n_far,n Locations beyond the reach, and locations predicted.
#' @param reach_text The reach, described.
#' @param nearest Smallest and largest distance from a counted location to
#'   its nearest training point.
#' @param calls Number of calls the counts are totalled over, when more than
#'   one.
#' @param refits \code{TRUE} when the calls predicted the same rows from
#'   several fits (a sweep of cross-validations), so that a location is
#'   counted once for each fit that put it out of reach.
#' @return The message.
#' @keywords internal
#' @noRd
.gwr_extrapolation_msg <- function(caller, n_far, n, reach_text, nearest,
                                   calls = 1L, refits = FALSE) {
  sprintf(paste0("%s(): %d of %d location(s) are farther from every ",
                 "training point than the widest kernel of the fit reaches ",
                 "(%s): the nearest training point is %s away%s. No local ",
                 "regression of the fit rested on points that far off, so ",
                 "the prediction(s) there are extrapolations in space. They ",
                 "are returned as computed; keep them out of a map or a ",
                 "score unless that is intended."),
          caller, n_far, n, reach_text,
          if (nearest[2L] > nearest[1L])
            sprintf("%s to %s", format(signif(nearest[1L], 4)),
                    format(signif(nearest[2L], 4)))
          else format(signif(nearest[1L], 4)),
          .gwr_totals_clause(calls, refits))
}

# " (totals over the N prediction calls that found any)", for a message that
# totals more than one call; with `refits`, also that a location predicted by
# several fits is counted once for each.
.gwr_totals_clause <- function(calls, refits = FALSE) {
  if (!(calls > 1L)) return("")
  sprintf(" (totals over the %d prediction calls that found any%s)", calls,
          if (isTRUE(refits))
            "; a location predicted by several fits is counted once for each"
          else "")
}

#' The text of the warning that counts predictions that are NA
#'
#' See \code{.gwr_extrapolation_msg()}.
#'
#' @param caller Name of the function the message is from.
#' @param n_na,n,n_empty Locations without a prediction, locations
#'   predicted, and those without a prediction where no training point
#'   carries weight.
#' @param bw_text The bandwidth, described.
#' @param calls Number of calls the counts are totalled over.
#' @param refits As for \code{.gwr_extrapolation_msg()}.
#' @return The message.
#' @keywords internal
#' @noRd
.gwr_not_estimable_msg <- function(caller, n_na, n, n_empty, bw_text,
                                   calls = 1L, refits = FALSE) {
  sprintf(paste0("%s(): %d of %d location(s) have no estimable local ",
                 "regression, so their predictions are NA; the other %d are ",
                 "unaffected%s. At %d of them no training point carries ",
                 "weight within the %s. At %d the training points in the ",
                 "window cannot place the new row: a predictor is constant ",
                 "among them at a value other than the new row's, the ",
                 "predictors there are exactly collinear, or fewer points ",
                 "carry weight than the row needs."),
          caller, n_na, n, n - n_na,
          .gwr_totals_clause(calls, refits),
          n_empty, bw_text, n_na - n_empty)
}


#' Raise predict.gwr_fit()'s counted warnings once for a loop of calls
#'
#' \code{predict.gwr_fit()} counts the locations it cannot estimate
#' (\code{"spatialkit_gwr_not_estimable"}) and those beyond the reach of the
#' fit (\code{"spatialkit_gwr_extrapolation"}) for the call it is in.  A
#' caller that predicts in pieces, a surface in chunks or a cross-validation
#' in folds, would pass on one warning per piece, each with that piece's
#' counts.  Wrapped in this, the pieces' warnings are collected from the
#' fields the conditions carry and each class is raised once, with the totals
#' and under the caller's name.
#'
#' Collectors nest: one that runs inside another raises its totals for the
#' outer one alone, which adds them to its own, so a sweep of
#' cross-validations (\code{cv_block_size_sweep()},
#' \code{select_features_forward()}) raises each class once for the sweep.
#'
#' @param expr The code that calls \code{predict()} one or more times.
#' @param caller Name of the calling function, for the messages.
#' @param n_total Number of locations predicted over all calls, when the
#'   caller knows it; otherwise the totals are over the calls that warned.
#' @param refits \code{TRUE} for a caller whose calls predict the same rows
#'   from several fits: the message then says that a location is counted once
#'   for each.
#' @return The value of \code{expr}.
#' @keywords internal
#' @noRd
.gwr_tally_predict_warnings <- function(expr, caller, n_total = NULL,
                                        refits = FALSE) {
  far <- list(calls = 0L, n_far = 0L, n = 0L, nearest = numeric(0),
              reach = numeric(0), reach_text = "")
  na  <- list(calls = 0L, n_na = 0L, n = 0L, n_empty = 0L,
              bw_text = character(0))
  depth <- .gwr_predict_state$tally %||% 0L
  .gwr_predict_state$tally <- depth + 1L
  on.exit(.gwr_predict_state$tally <- depth, add = TRUE)
  value <- withCallingHandlers(
    expr,
    spatialkit_gwr_extrapolation = function(w) {
      far$calls   <<- far$calls + (w$calls %||% 1L)
      far$n_far   <<- far$n_far + w$n_far
      far$n       <<- far$n + w$n
      far$nearest <<- range(c(far$nearest, w$nearest))
      far$reach   <<- range(c(far$reach, w$reach))
      far$reach_text <<- w$reach_text
      invokeRestart("muffleWarning")
    },
    spatialkit_gwr_not_estimable = function(w) {
      na$calls   <<- na$calls + (w$calls %||% 1L)
      na$n_na    <<- na$n_na + w$n_na
      na$n       <<- na$n + w$n
      na$n_empty <<- na$n_empty + w$n_empty
      na$bw_text <<- union(na$bw_text, w$bw_text)
      invokeRestart("muffleWarning")
    })
  # Back at the caller's depth before the totals are raised, so that they
  # are logged here unless the caller is collected in turn.
  .gwr_predict_state$tally <- depth
  if (na$n_na > 0L) {
    n <- if (is.null(n_total)) na$n else as.integer(n_total)
    # One fit has one bandwidth.  Folds that each searched for their own do
    # not, and the message then names no single one.
    bw_text <- if (length(na$bw_text) == 1L) na$bw_text
      else "bandwidth of the fit that predicted it (the fits' bandwidths differ)"
    .gwr_warn("spatialkit_gwr_not_estimable", "%s",
              .gwr_not_estimable_msg(caller, na$n_na, n, na$n_empty,
                                     bw_text, na$calls, refits),
              data = list(n_na = na$n_na, n = n, n_empty = na$n_empty,
                          bw_text = bw_text, calls = na$calls),
              tallied = TRUE)
  }
  if (far$n_far > 0L) {
    n <- if (is.null(n_total)) far$n else as.integer(n_total)
    # One fit has one reach.  The folds of a cross-validation are different
    # fits, each with its own: the message then gives the smallest and the
    # largest, and `reach` carries both.
    one_reach <- !(far$reach[2L] > far$reach[1L])
    reach      <- if (one_reach) far$reach[1L] else far$reach
    reach_text <- if (one_reach) far$reach_text
      else sprintf(paste0("%s to %s in the units of the fits' CRS: each call ",
                          "predicted from a fit of its own, with its own reach"),
                   format(signif(far$reach[1L], 4)),
                   format(signif(far$reach[2L], 4)))
    .gwr_warn("spatialkit_gwr_extrapolation", "%s",
              .gwr_extrapolation_msg(caller, far$n_far, n, reach_text,
                                     far$nearest, far$calls, refits),
              data = list(n_far = far$n_far, n = n, reach = reach,
                          reach_text = reach_text, nearest = far$nearest,
                          calls = far$calls),
              tallied = TRUE)
  }
  value
}


#' Survey the local collinearity of every fitting window
#'
#' At each observation, forms the kernel-weighted local design that the local
#' regression there inverts, and computes two condition indices of it.
#' \describe{
#'   \item{\code{cn}}{Belsley's scaled condition index of \eqn{W^{1/2} X},
#'     intercept included (\code{cbind(1, xmat)}), columns scaled to unit
#'     length but not centred.  This is the weighted design as it stands, and
#'     the literature's local index (Wheeler 2007;
#'     \code{GWmodel::gwr.collin.diagno()} computes the same).  It is right
#'     about the INTERCEPT: a predictor whose values in the window are far
#'     from 0 against their spread makes the local intercept an extrapolation
#'     to 0.  It is wrong about the slopes: a slope's precision,
#'     \eqn{\sigma^2 / \sum w (x - \bar{x}_w)^2}, does not depend on where the
#'     predictor's origin is, yet the same field in kelvin scored a median 454
#'     where in degrees C it scored 27.}
#'   \item{\code{cn_slopes}}{The index for the slopes.  The predictors are
#'     centred at their weighted mean in the window, which makes them exactly
#'     orthogonal to the intercept (so no origin enters), and each is divided
#'     by its standard deviation over the whole data, so the reference scale
#'     is the study area's own spread: \eqn{Z = \sqrt{w / \sum w}\,
#'     (x - \bar{x}_w) / s}, row by row.  The index is
#'     \code{max(1, d_max) / min(1, d_min)} over the singular values \code{d}
#'     of \eqn{Z}: 1 when the window's predictors vary as much, and as
#'     independently, as across the study area.  It is large when a predictor
#'     is (nearly) constant inside the window -- a regional covariate,
#'     collinear with the intercept and nothing else -- and when two
#'     predictors move together there.  The centring has to be LOCAL: centred
#'     at the global mean, a cluster whose covariate sits at the global mean
#'     scored near 1 although its local slope was undetermined.  Over the
#'     whole data the index reduces to the standard centred scaled condition
#'     index.}
#' }
#' A non-finite index counts as extreme.  \code{kappa()} returns \code{Inf}
#' for an exactly singular matrix, and \code{is.finite(cn) && cn > 1e6}
#' discarded precisely the worst case.  See \code{.gwr_slopes_collinear()}
#' for how the two indices are combined.
#'
#' The weights are the kernel's own (Wheeler and Tiefelsdorf 2005 diagnose
#' GWR collinearity on the weighted design), so a bisquare window's edge
#' points, which contribute almost nothing to the fit, contribute almost
#' nothing here either; rows with negligible weight are dropped before the
#' SVD.  So are rows whose weight is undefined (\code{NaN}, a zero-width
#' adaptive kernel at co-located points), and a window left with fewer rows
#' than columns counts as singular (\code{Inf}) in both indices: the local
#' solve finds it of deficient rank.  The survey drops rows whose weight is
#' below 1e-8 and the solve drops none, so under a Gaussian or exponential
#' kernel the survey can call a window singular that the solve finds of full
#' rank on distant, lightly weighted rows; the survey's verdict is then the
#' caution that the coefficients rest on very little.  Every location is
#' surveyed, not a sample: the map of coefficients needs a value at each, and
#' the survey's cost is the same order as the fit's.
#'
#' @param coords Matrix of coordinates, one row per observation.
#' @param xmat Numeric matrix of the numeric predictors.
#' @param adaptive Logical; TRUE when the bandwidth counts neighbours.
#' @param bw The bandwidth the model is fitted with.
#' @param kernel The kernel name.
#' @return A data.frame with one row per observation: \code{row}, \code{x},
#'   \code{y}, \code{n_window} (points with non-negligible weight), \code{cn}
#'   (the uncentred scaled condition index with the intercept) and
#'   \code{cn_slopes} (the slope index), each \code{Inf} when singular.
#' @keywords internal
#' @noRd
.gwr_local_collinearity <- function(coords, xmat, adaptive, bw, kernel) {
  n_obs <- nrow(coords)
  out <- data.frame(row = seq_len(n_obs), x = coords[, 1], y = coords[, 2],
                    n_window = NA_integer_, cn = NA_real_, cn_slopes = NA_real_)
  # One predictor is enough: X below adds the intercept, and a predictor
  # constant inside a window is collinear with it.
  if (!is.matrix(xmat) || ncol(xmat) < 1L || n_obs < 1L) return(out)
  if (is.null(bw) || !is.finite(bw)) return(out)
  X <- cbind(1, xmat)
  # The slope index's reference scale: each predictor's spread over the
  # whole data.  A predictor with none has no slope to estimate anywhere.
  s_glob <- apply(xmat, 2L, stats::sd)
  s_ok   <- all(is.finite(s_glob) & s_glob > 0)
  for (i in seq_len(n_obs)) {
    d <- sqrt((coords[, 1] - coords[i, 1])^2 + (coords[, 2] - coords[i, 2])^2)
    w <- .gw_kernel_weights(d, bw, kernel, adaptive)
    keep <- which(is.finite(w) & w > 1e-8)
    out$n_window[i] <- length(keep)
    if (length(keep) < ncol(X)) {
      out$cn[i] <- Inf
      out$cn_slopes[i] <- Inf
      next
    }
    out$cn[i] <- .condition_index(sqrt(w[keep]) * X[keep, , drop = FALSE])
    out$cn_slopes[i] <- .gwr_slope_index(w[keep], xmat[keep, , drop = FALSE],
                                         s_glob, s_ok)
  }
  out
}


#' The origin-free condition index of a window's slopes
#'
#' See \code{.gwr_local_collinearity()}.  A predictor that is exactly
#' constant inside the window is singular (\code{Inf}) outright: after
#' centring its column is rounding noise, whose size depends on the value it
#' is constant at.
#'
#' @param w Positive kernel weights of the window's rows.
#' @param xk The window's rows of the predictor matrix.
#' @param s_glob Each predictor's standard deviation over the whole data.
#' @param s_ok Whether every \code{s_glob} is finite and positive.
#' @return A number of at least 1, or \code{Inf}.
#' @keywords internal
#' @noRd
.gwr_slope_index <- function(w, xk, s_glob,
                             s_ok = all(is.finite(s_glob) & s_glob > 0)) {
  if (!isTRUE(s_ok) || nrow(xk) < ncol(xk) + 1L) return(Inf)
  if (any(apply(xk, 2L, function(v) max(v) == min(v)))) return(Inf)
  ww <- w / sum(w)
  z  <- sqrt(ww) * sweep(sweep(xk, 2L, colSums(ww * xk)), 2L, s_glob, "/")
  sv <- tryCatch(svd(z, nu = 0, nv = 0)$d, error = function(e) NULL)
  if (is.null(sv) || !length(sv) || any(!is.finite(sv))) return(Inf)
  if (min(sv) <= .Machine$double.eps * max(1, sv)) return(Inf)
  max(1, max(sv)) / min(1, min(sv))
}


#' Which windows' local slopes count as collinear
#'
#' A window's slopes are collinear when their own index (\code{cn_slopes}) is
#' above 30 or singular.  The uncentred index with the intercept (\code{cn})
#' says nothing about the slopes of a fit this package solved: the local
#' solve centres each predictor in its window before it factors the design,
#' so a predictor's origin does not enter.  With a predictor shifted by 1e7
#' (\code{cn} about 3e7) the slopes of a 200-point fit moved by 1e-9, which
#' is the precision the shifted column itself has.  Only a \code{cn} that is
#' singular to machine precision counts: the predictor's spread is then lost
#' in the rounding of its own values.
#'
#' A fit whose engine is GWmodel's object (one made by an earlier version)
#' was solved from the uncentred \eqn{X^\top W X}, whose condition number is
#' about \code{cn^2}: shifting a predictor's origin until \code{cn} reached
#' 1.6e6 moved its local slopes by 0.4\% of their spread, 1.6e7 by 25\%, and
#' 1.6e8 made GWmodel stop with "matrix is singular".  For such a fit a
#' \code{cn} above 1e6 counts as collinear too.  A fit of this package is
#' told by the \code{rank} column its survey carries.
#'
#' A survey without \code{cn_slopes} (a fit made before it existed) falls
#' back to \code{cn > 30}.  \code{NA} (a survey that did not run) counts as
#' collinear, as it always has.
#'
#' @param lc A survey, as \code{.gwr_local_collinearity()} returns it.
#' @param own_solve Whether the fit was solved by this package's local
#'   solve.  By default, whether the survey carries the solve's \code{rank}.
#' @return Logical vector, one element per row of \code{lc}.
#' @keywords internal
#' @noRd
.gwr_slopes_collinear <- function(lc, own_solve = !is.null(lc$rank)) {
  cn <- lc$cn
  if (is.null(lc$cn_slopes)) return(!is.finite(cn) | cn > 30)
  cs <- lc$cn_slopes
  out <- !is.finite(cs) | cs > 30 | !is.finite(cn)
  if (!isTRUE(own_solve)) out <- out | cn > 1e6
  out
}


#' Warn about the local collinearity survey's findings
#'
#' The two messages the sampled spot-check used to raise, now on the exact
#' fraction of locations, and about the slopes (\code{.gwr_slopes_collinear()}):
#' a window whose only problem is an ill-determined intercept (\code{cn}
#' above 30, slopes fine) is logged, not warned about.
#' @keywords internal
#' @noRd
.gwr_local_collinearity_warn <- function(local_cn_df) {
  if (is.null(local_cn_df) || !nrow(local_cn_df) || all(is.na(local_cn_df$cn)))
    return(invisible(NULL))
  n <- nrow(local_cn_df)
  # Called by fit_gwr_model() before the solve's columns are added to the
  # survey, for a fit the package's own solve is about to make.
  slope_bad <- .gwr_slopes_collinear(local_cn_df, own_solve = TRUE)
  n_extreme <- sum(slope_bad)
  frac <- n_extreme / n
  # What happens at those windows: a near-singular window returns
  # coefficients that are least squares and may be implausibly large; an
  # exactly singular one returns NA for the coefficients it cannot estimate,
  # and the warning raised after the fit counts them.
  if (frac > 0.25) {
    .warn_and_log(
      "fit_gwr_model(): local collinearity: %.0f%% of %d locations have a collinear local design (the predictors, centred in the kernel-weighted window, have scaled condition index > 30 against their spread over the study area, or the window is singular) at the bandwidth in use. Local regressions there are unstable and their coefficients may be implausibly large, and a coefficient that an exactly singular window cannot estimate is NA; plot(fit, type = \"coefficients\") masks them.",
      frac * 100, n
    )
  } else if (n_extreme > 0L) {
    .warn_and_log(
      "fit_gwr_model(): local collinearity: %d of %d locations have a collinear local design (the predictors, centred in the kernel-weighted window, have scaled condition index > 30 against their spread over the study area, or the window is singular) at the bandwidth in use.",
      n_extreme, n
    )
  }
  # The intercept alone: a predictor far from 0 against its local spread.
  # That says nothing about the slopes, so it is not a warning; the Intercept
  # map masks those locations.
  n_int <- sum(!slope_bad & (!is.finite(local_cn_df$cn) | local_cn_df$cn > 30))
  if (n_int > 0L)
    .log_info(paste0("fit_gwr_model(): the local intercept is an extrapolation ",
                     "at %d of %d locations (a predictor's local values are far ",
                     "from 0 against their spread: condition index with the ",
                     "intercept > 30). The slopes there are not affected; centre ",
                     "the predictor to map an interpretable intercept."),
              n_int, n)
  invisible(NULL)
}


#' GWR predictions at new locations, one local regression each
#'
#' Solves the local regression at every row of \code{rp} with
#' \code{.gwr_local_solve()}, on the kernel weights of the training points
#' (an adaptive bandwidth counts neighbours among the TRAINING points), and
#' returns \eqn{x^\top\hat\beta(u)}.  Each location stands alone: a window
#' that carries no weight, or whose rows cannot place the new row's
#' predictors, gives \code{NA} for that location and leaves the others as
#' they are.  At a training location with that row's own predictors the
#' value is \code{fitted()}'s, to the last bit.
#'
#' @param formula The fitted model formula.
#' @param sp_train SpatialPointsDataFrame of the training data, 2-D.
#' @param newdata_df data.frame of the new rows, with the predictors.
#' @param rp Two-column numeric matrix of the new rows' coordinates, in the
#'   training CRS.
#' @param bw,kernel,adaptive The fit's bandwidth settings.
#' @param se Also return the variance factor of each prediction.
#' @return A list, each element with one value per row of \code{rp}:
#'   \code{fitted} (\code{NA} where the local regression cannot be
#'   estimated), \code{rank}, \code{n_window}, \code{nearest} (distance to
#'   the nearest training point), \code{width} (the kernel's distance
#'   parameter there) and, with \code{se = TRUE}, \code{var_fit} (multiply by
#'   the residual variance for the variance of the fitted value).
#' @keywords internal
#' @noRd
.gwr_predict_at <- function(formula, sp_train, newdata_df, rp, bw, kernel,
                            adaptive, se = FALSE) {
  dsn <- .gwr_design(formula, sp_train@data, newdata_df)
  k   <- ncol(dsn$X)
  if (nrow(dsn$X_new) != nrow(rp) || ncol(dsn$X_new) != k)
    stop(sprintf(paste0("the design matrix for newdata is %d x %d, but the ",
                        "model has %d term(s) for %d location(s)."),
                 nrow(dsn$X_new), ncol(dsn$X_new), k, nrow(rp)), call. = FALSE)
  sv <- .gwr_local_solve(dsn$X, dsn$y, sp::coordinates(sp_train), bw, kernel,
                         adaptive, rp = rp, X_rp = dsn$X_new,
                         what = "predict", se = se)
  sv$fitted <- unname(sv$fitted)
  sv$fitted[!is.finite(sv$fitted)] <- NA_real_
  sv
}


