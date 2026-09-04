# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

#' Validate a GWR kernel name
#'
#' GWmodel accepts kernel names as character strings directly (unlike spgwr
#' which required function objects).
#'
#' \strong{Currently unreachable.}  Every entry point that takes a kernel --
#' \code{fit_gwr_model()}, \code{gwr_model_selection()} and \code{cv_gwr()} --
#' declares it as a \code{c("bisquare", ...)} default and runs
#' \code{match.arg()} on it, which rejects any value this function would have
#' to repair.  \code{cv_gwr()} (R/cross-validation.R) is its only caller and
#' calls it on the line \emph{after} its own \code{match.arg()}, so the
#' fallback branch below cannot execute.  It is kept, rather than deleted,
#' only because that caller lives in another file; if the redundant call there
#' is removed, remove this too.  Do not add a comment anywhere claiming it
#' "earns its keep" in \code{cv_gwr()} -- it does not.
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


#' Extract fitted or predicted values from a GWmodel GWR result
#'
#' **Unified extraction function** used by \code{fitted.gwr_fit()} (in-sample
#' evaluation) and \code{predict.gwr_fit()} (prediction at new locations), both
#' in R/model-classes.R.  Replaces the previously duplicated
#' .extract_gwr_fitted() and .extract_gwr_predictions() functions.
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
                                response_var = all.vars(formula)[1],
                                mode = c("auto", "basic", "predict")) {
  mode <- match.arg(mode)
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
        .log_warn(paste0(".extract_gwr_values(): GWmodel's SDF carries a ",
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
                       "term(s) matched a column in GWmodel's SDF, so the ",
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
#'   is an integer (number of neighbours); for fixed mode a distance in the
#'   units of the **projected** CRS the fit runs in -- \code{prep_model_data()}
#'   projects geographic input before the bandwidth is used, so 0.2 supplied
#'   for lon/lat data is 0.2 metres, not 0.2 degrees. Read the CRS off
#'   \code{sf::st_crs(fit$data_sf)}, or pass \code{target_crs} to
#'   \code{\link{ensure_projected}} beforehand to fix the units yourself.
#'   If NULL (default), bandwidth is selected automatically via
#'   \code{GWmodel::bw.gwr()}.  With \code{adaptive = FALSE}, a bandwidth
#'   smaller than a ten-thousandth of the data's extent raises a warning naming
#'   the extent and the CRS the fit runs in: every local window is then likely
#'   to be empty, which used to produce a fit whose coefficients were all
#'   \code{NaN} with nothing raised anywhere.
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
#' In addition, a **local** spot-check is performed at up to 30 locations --
#' every location when there are 30 or fewer, otherwise a fixed sample of 30
#' drawn under a constant seed, so the diagnostic is reproducible and the count
#' is not configurable.  For each sampled point the nearest neighbours within
#' the bandwidth window -- the bandwidth the model is actually fitted with, not
#' a stand-in -- are selected and the condition number of that local design
#' sub-matrix is evaluated.  That sub-matrix is the predictors **plus an
#' intercept column**, matching the design GWmodel fits, and is unweighted; the
#' global condition number is computed on the predictors alone, so the two
#' numbers are not directly comparable.  An indicator that is constant inside a
#' window is collinear with the intercept and with nothing else, which is why
#' the intercept has to be there.  A non-finite condition number counts as
#' extreme: `kappa()` returns `Inf` for an exactly singular design, which is the
#' worst case, not an exempt one.
#'
#' A warning is issued whenever **any** sampled location has a singular or
#' near-singular local design; the wording reports a percentage when more than
#' 25\% of sampled locations are affected and a count otherwise.  Both are real
#' R warnings, not log lines.
#'
#' After the fit, the local coefficient surfaces are scanned and a further
#' warning counts local regressions that came back non-finite -- their windows
#' were singular.  `fitted()`, `residuals()`, `summary()` and
#' [model_metrics()] all drop those rows, so when this warning fires the
#' metrics describe only the part of the study area that fitted.
#'
#' Because the local spot-check examines only a subset of locations, it may not
#' detect every problematic neighbourhood.  Users working with highly clustered
#' data or near-collinear predictors should consider a full local-collinearity
#' audit as a post-fit diagnostic.
#'
#' @return A \code{gwr_fit} object (inherits from \code{spatial_fit}).
#'   Supports \code{predict()}, \code{fitted()}, \code{residuals()},
#'   \code{coef()}, \code{summary()}, and \code{model_metrics()}.
#'   Model-specific metadata lives in \code{$info} (bandwidth, adaptive,
#'   kernel, AICc, and \code{bandwidth_is_fallback} -- \code{TRUE} when
#'   automatic selection failed and the arbitrary fallback was used).  The raw GWmodel result is in \code{$engine}.
#' @family model fitting
#' @examples
#' \donttest{
#' if (requireNamespace("GWmodel", quietly = TRUE) &&
#'     requireNamespace("sp", quietly = TRUE)) {
#'   library(sf)
#'   set.seed(1)
#'   n <- 60
#'   dat <- st_as_sf(
#'     data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), elev = rnorm(n)),
#'     coords = c("x", "y"), crs = 32632
#'   )
#'   dat$price <- 10 + 0.01 * st_coordinates(dat)[, 1] + 2 * dat$elev + rnorm(n)
#'   fit <- fit_gwr_model(dat, "price", "elev", bandwidth = 30)
#'   summary(fit)
#'   head(predict(fit, newdata = dat))   # newdata is re-projected if needed
#' }
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
  if (!requireNamespace("GWmodel", quietly = TRUE))
    stop("fit_gwr_model(): package 'GWmodel' is required. Install with install.packages('GWmodel').",
         call. = FALSE)
  if (!requireNamespace("sp", quietly = TRUE))
    stop("fit_gwr_model(): package 'sp' is required (for GWmodel interop).",
         call. = FALSE)
  # match.arg() is the whole of kernel validation: an invalid value never gets
  # past it.  (.validate_kernel() below is called by cv_gwr() but only ever
  # after that function's own match.arg(), so it is unreachable there too --
  # see its @noRd block.)
  kernel <- match.arg(kernel)

  # prep_model_data() accepts character(0) so an intercept-only spatial GP can
  # be fitted; a GWR needs at least one predictor to vary locally, and
  # reformulate() would otherwise fail with a bare message naming neither
  # this function nor the argument.
  if (!is.character(predictor_vars) || length(predictor_vars) < 1L)
    stop("fit_gwr_model(): `predictor_vars` must name at least one predictor; ",
         "there are no local coefficients to estimate otherwise.",
         call. = FALSE)

  # `bandwidth` is used unchecked below (twice: in the local-collinearity
  # spot-check and as the fitted bandwidth), and the clamping block further
  # down runs only in the adaptive branch.  Unvalidated, NA gives "missing
  # value where TRUE/FALSE needed", a length-2 vector gives "the condition has
  # length > 1", and with adaptive = FALSE a zero or negative distance reaches
  # GWmodel::gwr.basic() untouched.
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

  # GWmodel fits and predicts through two different code paths, and only one
  # of them expands contrasts: gwr.basic() builds its design with
  # model.matrix(), so a factor predictor fits cleanly, while gwr.predict()
  # indexes the prediction frame by the raw variable names and multiplies the
  # result, which for a factor column is "non-numeric argument to binary
  # operator".  predict.gwr_fit() catches that and returns all NA, so the
  # model appears to fit and then silently predicts nothing.  Reject the
  # column here instead, where it can be named.
  non_num <- predictor_vars[
    !vapply(sf::st_drop_geometry(dat)[, predictor_vars, drop = FALSE],
            is.numeric, logical(1))]
  if (length(non_num) > 0L)
    stop(sprintf(paste0("fit_gwr_model(): predictor(s) %s are not numeric. ",
                        "GWmodel fits a factor or character predictor (via ",
                        "model.matrix contrasts) but cannot predict from it ",
                        "-- gwr.predict() does not expand contrasts and would ",
                        "return all NA. Encode the column as numeric indicator ",
                        "column(s) yourself, and pass those as predictors."),
                 paste(sQuote(non_num), collapse = ", ")),
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
  if (length(predictor_vars) >= 2L) {
    num_preds <- predictor_vars[vapply(pred_df[predictor_vars], is.numeric, logical(1))]
    if (length(num_preds) >= 2L) {
      xmat <- as.matrix(pred_df[, num_preds, drop = FALSE])
      cn <- tryCatch(kappa(xmat, exact = FALSE), error = function(e) Inf)
      # A non-finite condition number is an EXACTLY singular design -- the worst
      # case there is -- and `is.finite(cn) && ...` silently let it through.
      if (!is.finite(cn) || cn > 1e6) {
        .warn_and_log(
          "fit_gwr_model(): global predictor matrix condition number = %s (collinearity risk). Note: local collinearity within bandwidth windows may be substantially worse than this global value.",
          if (is.finite(cn)) sprintf("%.0f", cn) else "infinite (exactly singular)"
        )
      }
      # The LOCAL collinearity spot-check is deferred until the bandwidth is
      # known -- see .gwr_local_collinearity_check() below the bandwidth
      # selection.  Running it here meant that on the default path
      # (bandwidth = NULL, which is most calls) it used a stand-in window --
      # min(50, n) neighbours when adaptive, the whole data set when fixed --
      # instead of the window the model is actually fitted with, so the
      # documented warning simply never fired: identical data warned when the
      # bandwidth was passed explicitly and stayed silent when the same value
      # was selected automatically.
      local_collinearity <- list(xmat = xmat, num_preds = num_preds)
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
  bw_fallback_raw       <- NA_real_
  # A fixed bandwidth is a length in the CRS the fit runs in, and that is not
  # necessarily the CRS the caller was looking at: prep_model_data() projects
  # geographic input, so 0.2 "degrees" becomes 0.2 metres and every local
  # regression is empty -- all coefficients NaN, all fitted values NA,
  # summary() reporting n = 0 -- with nothing raised anywhere.  Say so at the
  # moment the mismatch is visible.
  if (!is.null(bandwidth) && !adaptive) {
    fit_units <- tryCatch(sf::st_crs(dat)$units_gdal, error = function(e) NULL)
    extent_x  <- suppressWarnings({
      bb <- sf::st_bbox(dat); as.numeric(bb[["xmax"]] - bb[["xmin"]])
    })
    if (is.finite(extent_x) && extent_x > 0 &&
        as.numeric(bandwidth) < extent_x / 1e4) {
      .warn_and_log(paste0("fit_gwr_model(): a fixed bandwidth of %s is less ",
                           "than a ten-thousandth of the data's extent (%s %s ",
                           "across). The bandwidth is a distance in the CRS the ",
                           "fit runs in (%s), which prep_model_data() may have ",
                           "chosen -- geographic input is projected first, so a ",
                           "value in degrees is read as metres. Every local ",
                           "window is likely to be empty."),
                    format(as.numeric(bandwidth)), format(signif(extent_x, 4)),
                    if (is.null(fit_units) || is.na(fit_units)) "unit"
                    else fit_units,
                    tryCatch(sf::st_crs(dat)$input, error = function(e) "unknown"))
    }
  }
  if (is.null(bandwidth)) {
    # .gwr_quietly(): bw.gwr() writes its golden-section search trace with bare
    # cat(), which neither suppressMessages() nor suppressWarnings() touches.
    # cv_gwr(bandwidth = NULL) calls this once per fold, so without it a
    # five-fold CV dumps five full traces.  gwr_model_selection() already does
    # the same thing and documents why.
    bw <- tryCatch(
      .gwr_quietly(suppressWarnings(
        GWmodel::bw.gwr(fml, data = sp_dat, approach = "AICc",
                        kernel = kernel, adaptive = adaptive)
      )),
      error = function(e) {
        .log_warn("fit_gwr_model(): bw.gwr() failed: %s", conditionMessage(e))
        NA_real_
      }
    )

    # Test the LENGTH first.  This condition sits outside the tryCatch above,
    # so a bw.gwr() that returns numeric(0) or a length-2 vector -- neither is
    # an error, so neither is caught -- reached `!is.finite(bw)` and raised a
    # bare "'length = 2' in coercion to 'logical(1)'" from a line that names
    # neither the function nor the cause.  Anything that is not one usable
    # positive number takes the fallback instead.
    if (length(bw) != 1L || !is.numeric(bw) || !is.finite(bw) || bw <= 0) {
      bw <- .fallback_bandwidth(sp_dat, adaptive)
      bandwidth_is_fallback <- TRUE
      bw_fallback_raw       <- as.numeric(bw)
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

  # The fallback warning is issued AFTER the clamp, not before it.  The clamp
  # can move the value -- a 6-row dataset clamps the arbitrary 10 down to 6 --
  # and a warning naming a bandwidth the fit did not use sends the reader
  # looking for a number that appears nowhere in the result.  Both are named
  # when they differ.
  if (bandwidth_is_fallback) {
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

  # Now that the bandwidth is known, run the local collinearity spot-check on
  # the window the model will actually be fitted with.
  if (!is.null(local_collinearity))
    .gwr_local_collinearity_check(sf::st_coordinates(dat),
                                  local_collinearity$xmat,
                                  adaptive, bw, n_obs)

  # --- Fit GWR ---
  fit <- tryCatch(
    GWmodel::gwr.basic(formula = fml, data = sp_dat, bw = bw,
                       kernel = kernel, adaptive = adaptive),
    error = function(e) {
      stop(sprintf("fit_gwr_model(): GWR fit failed: %s", conditionMessage(e)),
           call. = FALSE)
    }
  )
  
  # A local regression whose design is singular comes back as NaN rather than
  # as an error, and nothing downstream says so: fitted(), residuals(),
  # summary() and model_metrics() all drop the non-finite rows, so a fit in
  # which 182 of 200 local regressions failed reported n = 18 and R2 = 0.96 as
  # though that were the whole model.  Count them and say so.
  n_bad_local <- tryCatch({
    sdf <- fit$SDF
    if (is.null(sdf)) 0L else {
      b <- as.data.frame(sdf)
      cols <- intersect(c("Intercept", predictor_vars), names(b))
      if (!length(cols)) 0L
      else {
        ok <- vapply(b[cols],
                     function(z) is.finite(suppressWarnings(as.numeric(z))),
                     logical(nrow(b)))
        if (!is.matrix(ok)) ok <- matrix(ok, nrow = nrow(b))
        sum(!apply(ok, 1L, all))
      }
    }
  }, error = function(e) 0L)
  if (n_bad_local > 0L)
    .warn_and_log(paste0("fit_gwr_model(): %d of %d local regression(s) ",
                         "returned non-finite coefficients -- their windows are ",
                         "singular. fitted(), residuals() and every metric are ",
                         "computed from the %d that succeeded, so they describe ",
                         "part of the study area only."),
                  n_bad_local, n_obs, n_obs - n_bad_local)

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


#' Spot-check local collinearity inside the fitting window
#'
#' Samples locations, forms the LOCAL design matrix each of them will be fitted
#' with, and counts how many are numerically singular or nearly so.
#'
#' Two things this deliberately does differently from the global check:
#' \itemize{
#'   \item The intercept column is included (\code{cbind(1, xmat)}).  GWmodel
#'     fits an intercept, and the case this check exists for -- an indicator
#'     that is constant inside a window -- is collinear with the INTERCEPT and
#'     with nothing else, so a check on the predictors alone cannot see it.
#'   \item A non-finite condition number counts as extreme.  \code{kappa()}
#'     returns \code{Inf} for an exactly singular matrix, and
#'     \code{is.finite(cn) && cn > 1e6} discarded precisely the worst case.
#' }
#'
#' @param coords Matrix of coordinates, one row per observation.
#' @param xmat Numeric matrix of the numeric predictors.
#' @param adaptive Logical; TRUE when the bandwidth counts neighbours.
#' @param bw The bandwidth the model will actually be fitted with.
#' @param n_obs Number of observations.
#' @return Invisibly \code{NULL}; called for the warning.
#' @keywords internal
#' @noRd
.gwr_local_collinearity_check <- function(coords, xmat, adaptive, bw, n_obs) {
  if (!is.matrix(xmat) || ncol(xmat) < 2L || n_obs < 1L) return(invisible(NULL))
  if (is.null(bw) || !is.finite(bw)) return(invisible(NULL))
  n_spot <- min(30L, n_obs)
  spot_idx <- if (n_obs <= 30L) {
    seq_len(n_obs)
  } else {
    cleanup_spot <- .with_seed(42L)
    on.exit(cleanup_spot(), add = TRUE)
    sample.int(n_obs, n_spot)
  }
  n_extreme <- 0L
  for (si in spot_idx) {
    dists <- sqrt((coords[, 1] - coords[si, 1])^2 +
                    (coords[, 2] - coords[si, 2])^2)
    if (adaptive) {
      nn_idx <- order(dists)[seq_len(min(max(1L, as.integer(bw)), n_obs))]
    } else {
      nn_idx <- which(dists <= as.numeric(bw))
      if (length(nn_idx) < ncol(xmat) + 1L)
        nn_idx <- order(dists)[seq_len(min(ncol(xmat) + 1L, n_obs))]
    }
    local_xmat <- cbind(1, xmat[nn_idx, , drop = FALSE])
    local_cn <- tryCatch(kappa(local_xmat, exact = FALSE),
                         error = function(e) Inf)
    if (!is.finite(local_cn) || local_cn > 1e6) n_extreme <- n_extreme + 1L
  }
  frac <- n_extreme / n_spot
  if (frac > 0.25) {
    .warn_and_log(
      "fit_gwr_model(): local collinearity spot-check: %.0f%% of %d sampled locations have a singular or near-singular local design (condition number > 1e6) at the bandwidth in use. Local regressions there are unstable and their coefficients may come back non-finite.",
      frac * 100, n_spot
    )
  } else if (n_extreme > 0L) {
    .warn_and_log(
      "fit_gwr_model(): local collinearity spot-check: %d of %d sampled locations have a singular or near-singular local design (condition number > 1e6) at the bandwidth in use.",
      n_extreme, n_spot
    )
  }
  invisible(NULL)
}
