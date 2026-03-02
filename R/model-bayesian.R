#' Fit a Bayesian spatial regression with a 2D Gaussian Process (via brms)
#'
#' @param data_sf An sf object with response, predictors, and geometries.
#' @param response_var Response column name.
#' @param predictor_vars Predictor column names.
#' @param family A brms family object. Default NULL (resolved to
#'   \code{brms::gaussian()} lazily after the brms availability check).
#' @param gp_k Positive integer for GP rank, or NULL (default) for automatic
#'   selection based on dataset size: min(n/3, max(15, sqrt(n))).
#' @param gp_c Positive numeric for GP scale. Default 1.5.
#' @param prior Optional brms prior specification. When NULL and
#'   \code{standardize_predictors = TRUE}, weakly informative
#'   \code{normal(0, 5)} priors are set on regression coefficients.
#'   A data-informed GP length-scale prior is always appended
#'   automatically unless \code{prior} already contains an entry with
#'   \code{class = "lscale"}.
#' @param chains Number of MCMC chains. Default 4.
#' @param iter Total iterations per chain. Default 2000.
#' @param warmup Warmup iterations. Default floor(iter/2).
#' @param cores Number of parallel cores.
#' @param seed Integer seed. Default 123.
#' @param backend "auto", "cmdstanr", or "rstan".
#' @param control Named list of sampler controls.
#' @param compute_loo Logical; compute PSIS-LOO. Default TRUE.
#' @param standardize_predictors Logical; center and scale numeric predictors
#'   before fitting. Default FALSE. When TRUE, the scaling parameters are
#'   stored in the return value so predictions can be computed correctly.
#' @param check_convergence Logical; after fitting, check for divergences,
#'   low ESS, and high R-hat and issue warnings. Default TRUE.
#' @param pointize Strategy for non-point geometry coercion.
#' @param boundary Optional polygonal sf/sfc for CRS harmonization.
#' @param .already_prepped Logical (internal). If \code{TRUE}, skip the
#'   \code{prep_model_data()} call because the caller has already projected,
#'   coerced, and filtered the data.  Used by the CV internals to avoid a
#'   redundant second pass on every fold.  End users should leave this at the
#'   default \code{FALSE}.
#' @details
#' \strong{Coordinate scaling and anisotropy.}
#' Before fitting the GP, X and Y coordinates are each centered and divided by
#' their own standard deviation (lines 110–111).
#' Because the two axes are scaled independently, an isotropic
#' squared-exponential kernel in the \emph{scaled} space corresponds to an
#' \strong{anisotropic} kernel in the original CRS: the effective length-scale
#' in the X direction (in CRS units) differs from the Y direction whenever
#' \code{sd(X) != sd(Y)}.
#'
#' This per-axis standardization is deliberate — it stabilises the GP
#' numerically when the coordinate extents differ dramatically (common in
#' projected CRSs where easting and northing span very different ranges) — but
#' users who expect the GP to be isotropic in geographic distance should be
#' aware of this behaviour.
#'
#' If true isotropy in the original CRS is desired, one could use a single
#' scaling factor such as \code{max(sd(X), sd(Y))} for both axes.
#' The stored \code{$info$coord_scaling} list includes a
#' \code{scaling_type} element (\code{"anisotropic"}) so downstream code can
#' detect which strategy was used.
#'
#' @return A \code{bayesian_fit} object (inherits from \code{spatial_fit}).
#'   Supports \code{predict()}, \code{fitted()}, \code{residuals()},
#'   \code{coef()}, \code{summary()}, and \code{model_metrics()}.
#'   Model-specific metadata lives in \code{$info} (coord_scaling,
#'   predictor_scaling, gp_k, loo, looic, convergence_ok,
#'   convergence_diagnostics).  The raw brmsfit is in \code{$engine}.
#' @export
fit_bayesian_spatial_model <- function(
    data_sf, response_var, predictor_vars,
    family      = NULL,
    
    gp_k        = NULL,
    gp_c        = 1.5,
    prior       = NULL,
    chains      = 4,
    iter        = 2000,
    warmup      = floor(iter / 2),
    cores       = max(1L, parallel::detectCores() - 1L),
    seed        = 123,
    backend     = c("auto", "cmdstanr", "rstan"),
    control     = list(adapt_delta = 0.9, max_treedepth = 12),
    compute_loo = TRUE,
    standardize_predictors = FALSE,
    check_convergence = TRUE,
    pointize    = "auto",
    boundary    = NULL,
    .already_prepped = FALSE
) {
  if (!inherits(data_sf, "sf"))
    stop("fit_bayesian_spatial_model(): `data_sf` must be an sf object.")
  if (!requireNamespace("brms", quietly = TRUE))
    stop("fit_bayesian_spatial_model(): package 'brms' is required.")

  # Resolve family lazily — avoids evaluating brms::gaussian() at
  # source time, which would defeat the requireNamespace() guard above.
  if (is.null(family)) family <- brms::gaussian()

  backend <- match.arg(backend)
  if (identical(backend, "auto"))
    backend <- if (requireNamespace("cmdstanr", quietly = TRUE)) "cmdstanr" else "rstan"

  # When called from CV internals the data is already prepped; skip the
  # redundant pass to avoid re-projecting, re-coercing geometry, and
  # re-scanning for NAs on every fold.
  if (isTRUE(.already_prepped)) {
    dat_sf <- data_sf
  } else {
    dat_sf <- prep_model_data(
      data_sf = data_sf, response_var = response_var,
      predictor_vars = predictor_vars, boundary = boundary, pointize = pointize
    )
  }

  gtypes <- as.character(sf::st_geometry_type(dat_sf, by_geometry = TRUE))
  if (!all(gtypes %in% c("POINT", "MULTIPOINT")))
    stop("fit_bayesian_spatial_model(): geometry must be POINT/MULTIPOINT after prep.")

  # NOTE: prep_model_data() already guarantees a projected CRS, so no
  # additional .is_longlat() / ensure_projected() call is needed here.
  
  if (is.null(gp_k)) {
    n_dat <- nrow(dat_sf)
    gp_k <- as.integer(min(n_dat / 3, max(15L, floor(sqrt(n_dat)))))
    gp_k <- max(5L, gp_k)  # absolute floor of 5
    .log_info("fit_bayesian_spatial_model(): auto-selected gp_k = %d for n = %d.", gp_k, n_dat)
  }
  gp_k <- as.integer(gp_k)

  coords <- sf::st_coordinates(dat_sf)
  if (!all(c("X", "Y") %in% colnames(coords))) colnames(coords)[1:2] <- c("X", "Y")
  
  # Per-axis (anisotropic) standardization: each coordinate is centered and
  # divided by its own SD.  This means the GP kernel is isotropic in the
  # *scaled* space but anisotropic in the original CRS whenever sd(X) != sd(Y).
  # See @details in the roxygen block above.
  x_raw <- coords[, "X"]
  y_raw <- coords[, "Y"]
  x_center <- mean(x_raw);  x_scale <- max(sd(x_raw), .Machine$double.eps)
  y_center <- mean(y_raw);  y_scale <- max(sd(y_raw), .Machine$double.eps)
  coord_scaling <- list(
    x_center = x_center, x_scale = x_scale,
    y_center = y_center, y_scale = y_scale,
    scaling_type = "anisotropic"
  )
  dat_df <- cbind(
    sf::st_drop_geometry(dat_sf),
    `..x` = (x_raw - x_center) / x_scale,
    `..y` = (y_raw - y_center) / y_scale
  )

  # ---- Data-informed GP length-scale prior ----
  # gp_lengthscale_bounds() is calibrated for the squared-exponential
  # kernel that brms::gp() uses.  We derive a weakly informative

  # normal(0, ·) prior on lscale from the distance structure of the
  # (scaled) coordinates.
  ls_bounds <- gp_lengthscale_bounds(
    cbind(dat_df[["..x"]], dat_df[["..y"]])
  )
  ls_prior_sd <- max(ls_bounds[["upper"]], ls_bounds[["lower"]] * 1.5)
  .log_info(
    "fit_bayesian_spatial_model(): GP length-scale bounds [%.4f, %.4f] on scaled coords; using normal(0, %.4f) lscale prior.",
    ls_bounds[["lower"]], ls_bounds[["upper"]], ls_prior_sd
  )

  # Optional predictor standardization — store transform params so
  # predictions can be computed correctly on new data.
  predictor_scaling <- NULL
  if (isTRUE(standardize_predictors) && length(predictor_vars) > 0L) {
    predictor_scaling <- list()
    for (pv in predictor_vars) {
      if (is.numeric(dat_df[[pv]])) {
        pv_mean <- mean(dat_df[[pv]], na.rm = TRUE)
        pv_sd   <- max(sd(dat_df[[pv]], na.rm = TRUE), .Machine$double.eps)
        dat_df[[pv]] <- (dat_df[[pv]] - pv_mean) / pv_sd
        predictor_scaling[[pv]] <- list(center = pv_mean, scale = pv_sd)
      }
    }
    if (length(predictor_scaling) > 0L)
      .log_info("fit_bayesian_spatial_model(): standardized %d numeric predictor(s).",
                length(predictor_scaling))
  }

  if (length(predictor_vars) == 0L) {
    rhs_terms <- "1"
  } else {
    base_fml  <- stats::reformulate(termlabels = predictor_vars)
    rhs_terms <- as.character(base_fml)[2L]
  }
  gp_term <- sprintf("gp(..x, ..y, k = %s, c = %s)", gp_k, gp_c)
  fml <- stats::as.formula(sprintf("%s ~ %s + %s", response_var, rhs_terms, gp_term))

  # Build priors in two independent steps:
  # 1. Regression coefficient priors (only when no user prior and predictors
  #    are standardized).
  # 2. GP length-scale prior (always, unless the user's prior already
  #    includes an "lscale" class entry).

  if (is.null(prior)) {
    prior_parts <- list()
    if (isTRUE(standardize_predictors) && length(predictor_vars) > 0L) {
      prior_parts <- c(prior_parts, list(
        brms::set_prior("normal(0, 5)", class = "b")
      ))
      .log_info("fit_bayesian_spatial_model(): using weakly informative normal(0,5) priors on standardized coefficients.")
    }
    prior <- if (length(prior_parts) > 0L) Reduce(`+`, prior_parts) else NULL
  }

  # Always append the data-informed GP length-scale prior unless the

  # user explicitly supplied an lscale prior.
  user_has_lscale <- !is.null(prior) &&
    inherits(prior, "brmsprior") &&
    any(prior$class == "lscale")

  if (!user_has_lscale) {
    ls_prior <- brms::set_prior(
      sprintf("normal(0, %s)", ls_prior_sd), class = "lscale"
    )
    prior <- if (is.null(prior)) ls_prior else prior + ls_prior
    .log_info(
      "fit_bayesian_spatial_model(): appending data-informed GP length-scale prior normal(0, %.4f).",
      ls_prior_sd
    )
  } else {
    .log_info(
      "fit_bayesian_spatial_model(): user-supplied prior already includes lscale class; skipping automatic GP length-scale prior."
    )
  }

  brm_args <- list(
    formula = fml, data = dat_df, family = family, prior = prior,
    chains = chains, iter = iter, warmup = warmup, cores = cores,
    seed = seed, control = control
  )
  if (backend %in% c("cmdstanr", "rstan")) brm_args$backend <- backend

  fit <- try(do.call(brms::brm, brm_args), silent = TRUE)
  if (inherits(fit, "try-error"))
    stop(sprintf("fit_bayesian_spatial_model(): brms fit failed: %s", as.character(fit)))

  # ---- Convergence diagnostics ----
  convergence_ok <- TRUE
  convergence_diagnostics <- list()
  if (isTRUE(check_convergence) && inherits(fit, "brmsfit")) {
    # Divergent transitions
    np <- tryCatch(brms::nuts_params(fit), error = function(e) NULL)
    if (!is.null(np)) {
      n_divergent <- sum(np$Value[np$Parameter == "divergent__"], na.rm = TRUE)
      convergence_diagnostics$n_divergent <- n_divergent
      if (n_divergent > 0L) {
        convergence_ok <- FALSE
        .log_warn(
          "fit_bayesian_spatial_model(): %d divergent transition(s) detected. Consider increasing adapt_delta (currently %.2f) or reparameterizing the model.",
          n_divergent, control$adapt_delta %||% 0.8
        )
      }
    }
    
    # R-hat
    rhat_vals <- tryCatch({
      rh <- brms::rhat(fit)
      if (is.numeric(rh)) rh else NULL
    }, error = function(e) NULL)
    
    if (!is.null(rhat_vals)) {
      max_rhat <- max(rhat_vals, na.rm = TRUE)
      convergence_diagnostics$max_rhat <- max_rhat
      bad_rhat <- names(rhat_vals)[rhat_vals > 1.05]
      if (length(bad_rhat) > 0L) {
        convergence_ok <- FALSE
        .log_warn(
          "fit_bayesian_spatial_model(): %d parameter(s) have R-hat > 1.05 (max = %.3f): %s. The model may not have converged.",
          length(bad_rhat), max_rhat,
          paste(head(bad_rhat, 5), collapse = ", ")
        )
      }
    }
    
    # Effective sample size ratio
    neff_vals <- tryCatch({
      ne <- brms::neff_ratio(fit)
      if (is.numeric(ne)) ne else NULL
    }, error = function(e) NULL)
    
    if (!is.null(neff_vals)) {
      min_neff <- min(neff_vals, na.rm = TRUE)
      convergence_diagnostics$min_neff_ratio <- min_neff
      low_neff <- names(neff_vals)[neff_vals < 0.1]
      if (length(low_neff) > 0L) {
        convergence_ok <- FALSE
        .log_warn(
          "fit_bayesian_spatial_model(): %d parameter(s) have effective sample size ratio < 0.1 (min = %.3f): %s. Consider running longer chains.",
          length(low_neff), min_neff,
          paste(head(low_neff, 5), collapse = ", ")
        )
      }
    }
  }

  loo_obj <- NULL; looic <- NA_real_
  if (isTRUE(compute_loo)) {
    loo_try <- try(brms::loo(fit), silent = TRUE)
    if (!inherits(loo_try, "try-error")) {
      loo_obj <- loo_try
      est <- try(loo_obj$estimates, silent = TRUE)
      if (!inherits(est, "try-error") && is.matrix(est)) {
        rn <- rownames(est)
        if ("looic" %in% rn)
          looic <- suppressWarnings(as.numeric(est["looic", "Estimate"]))
        else if ("elpd_loo" %in% rn)
          looic <- suppressWarnings(as.numeric(-2 * est["elpd_loo", "Estimate"]))
      }
    } else {
      .log_warn("fit_bayesian_spatial_model(): LOO computation failed.")
    }
  }

  new_spatial_fit(
    subclass       = "bayesian_fit",
    engine         = fit,
    formula        = fml,
    response_var   = response_var,
    predictor_vars = predictor_vars,
    data_sf        = dat_sf,
    info           = list(
      loo                      = loo_obj,
      looic                    = looic,
      coords                   = c("..x", "..y"),
      coord_scaling            = coord_scaling,
      gp_k                     = gp_k,
      gp_lengthscale_bounds    = ls_bounds,
      convergence_ok           = convergence_ok,
      convergence_diagnostics  = convergence_diagnostics,
      predictor_scaling        = predictor_scaling
    )
  )
}
