# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

#' Filter a named argument list to only parameters accepted by a function
#' @keywords internal
#' @noRd
.filter_args <- function(fun, args_list) {
  fml <- try(formals(fun), silent = TRUE)
  if (inherits(fml, "try-error")) return(list())
  if ("..." %in% names(fml)) return(args_list)
  args_list[names(args_list) %in% names(fml)]
}


#' Check whether a model backend is available
#' @keywords internal
#' @noRd
.model_available <- function(model_name) {
  switch(model_name,
    "GWR" = requireNamespace("GWmodel", quietly = TRUE) &&
      requireNamespace("sp", quietly = TRUE),
    "Bayesian" = requireNamespace("brms", quietly = TRUE),
    FALSE
  )
}


# ---------------------------------------------------------------------------
# Post-fit residual spatial autocorrelation check (Moran's I)
# ---------------------------------------------------------------------------

#' Build a row-standardised k-nearest-neighbour sparse weight matrix
#'
#' For each observation the \code{k} closest neighbours receive weight 1;
#' all other pairs receive weight 0.
#' The resulting matrix is then row-standardised so that each row sums to 1.
#' This is the standard default in spatial statistics (Anselin, 1988) and is
#' far more robust to irregularly-spaced or clustered data than an
#' inverse-distance scheme, which gives enormous weight to very close pairs
#' and can inflate Moran's I significance.
#'
#' Uses \pkg{FNN} for O(n·k) kd-tree nearest-neighbour lookup when available,
#' avoiding the O(n²) full distance matrix.  Returns a
#' \code{Matrix::sparseMatrix} (dgCMatrix), which keeps memory proportional to
#' O(n·k) instead of O(n²).  Falls back to the dense brute-force path when
#' \pkg{FNN} or \pkg{Matrix} are not installed.
#'
#' @param coords Numeric matrix (n x 2) of projected coordinates.
#' @param k Integer number of nearest neighbours (default 8).
#' @return A row-standardised weight matrix W — sparse (dgCMatrix) when
#'   \pkg{Matrix} is available, otherwise a dense base matrix.
#' @keywords internal
#' @noRd
.build_knn_weights <- function(coords, k = 8L) {
  n <- nrow(coords)
  k <- min(as.integer(k), n - 1L)
  if (k < 1L) k <- 1L

  has_fnn    <- requireNamespace("FNN",    quietly = TRUE)
  has_matrix <- requireNamespace("Matrix", quietly = TRUE)

  if (has_fnn && has_matrix) {
    # --- Fast path: O(n*k) kd-tree lookup + sparse matrix ----
    nn_idx <- FNN::get.knn(coords, k = k)$nn.index          # n x k matrix
    row_i  <- rep(seq_len(n), each = k)
    col_j  <- as.integer(t(nn_idx))
    W <- Matrix::sparseMatrix(
      i = row_i, j = col_j, x = 1 / k,                      # row-standardised
      dims = c(n, n), repr = "C"
    )
  } else {
    # --- Fallback: dense O(n²) path when packages are missing ----
    if (!has_fnn && n > 5000L)
      stop("n = ", n, " requires FNN for k-NN weights (dense fallback would allocate an n*n matrix). Install FNN with install.packages(\"FNN\").", call. = FALSE)
    if (has_fnn) {
      nn_idx <- FNN::get.knn(coords, k = k)$nn.index
    } else {
      dmat <- as.matrix(stats::dist(coords))
      diag(dmat) <- Inf
      nn_idx <- t(apply(dmat, 1, function(row) order(row)[seq_len(k)]))
    }
    W <- matrix(0, n, n)
    for (i in seq_len(n)) W[i, nn_idx[i, ]] <- 1
    rs <- rowSums(W)
    rs[rs == 0] <- 1
    W <- W / rs
  }

  W
}


#' Compute Moran's I on the residuals of a fitted spatial model
#'
#' Given a \code{spatial_fit} object (GWR or Bayesian), extracts the
#' residuals and the observation coordinates, builds a spatial weight
#' matrix, and computes the Moran's I statistic together with
#' its analytical expectation and variance under the randomisation
#' assumption (Cliff & Ord).  A z-score and two-sided p-value are
#' provided so the caller can assess whether statistically significant
#' spatial autocorrelation remains after fitting.
#'
#' By default, weights are constructed as a k-nearest-neighbour (k = 8)
#' binary matrix, row-standardised.  Users may supply their own weight
#' matrix via the \code{weights} argument.
#'
#' @param fit A \code{spatial_fit} object (from \code{fit_gwr_model} or
#'   \code{fit_bayesian_spatial_model}).
#' @param alternative Character: \code{"two.sided"} (default),
#'   \code{"greater"} (positive autocorrelation), or \code{"less"}.
#' @param weights Optional user-supplied n x n weight matrix.  When
#'   \code{NULL} (the default), a k-nearest-neighbour binary weight matrix
#'   (k = 8, row-standardised) is built from the observation coordinates.
#'   If a non-row-standardised matrix is supplied (i.e. rows do not all sum
#'   to 1), the Cliff & Ord variance formula is still valid for general W
#'   and the computation proceeds, but a warning is emitted because the
#'   magnitude of I is not directly comparable to results obtained with
#'   row-standardised weights.
#' @param k Integer number of nearest neighbours used when building the
#'   default weight matrix (ignored when \code{weights} is supplied).
#'   Default 8.
#' @return A list with components:
#'   \describe{
#'     \item{observed}{Numeric scalar, Moran's I statistic.}
#'     \item{expected}{Expected I under the null of no spatial
#'       autocorrelation, \eqn{-1/(n-1)}.}
#'     \item{sd}{Standard deviation of I under the randomisation
#'       assumption.}
#'     \item{z}{Standardised z-score, \eqn{(I - E[I]) / sd(I)}.}
#'     \item{p_value}{Two-sided (or one-sided) p-value from the
#'       normal approximation.}
#'     \item{n}{Number of observations used.}
#'   }
#'   Returns \code{NULL} with a warning if computation fails (e.g. fewer
#'   than 4 valid residuals).
#' @export
residual_morans_i <- function(fit,
                              alternative = c("two.sided", "greater", "less"),
                              weights = NULL,
                              k = 8L) {
  alternative <- match.arg(alternative)

  if (!inherits(fit, "spatial_fit")) {
    .log_warn("residual_morans_i(): `fit` is not a spatial_fit object.")
    return(NULL)
  }

  # --- Extract residuals & coordinates ---
  resid <- tryCatch(residuals(fit), error = function(e) NULL)
  if (is.null(resid) || length(resid) < 4L) {
    .log_warn("residual_morans_i(): could not extract enough residuals (n < 4).")
    return(NULL)
  }

  coords <- tryCatch({
    pts <- ensure_projected(coerce_to_points(fit$data_sf, "auto"))
    sf::st_coordinates(pts)[, 1:2, drop = FALSE]
  }, error = function(e) NULL)

  if (is.null(coords) || nrow(coords) != length(resid)) {
    .log_warn("residual_morans_i(): coordinate extraction failed.")
    return(NULL)
  }

  # Drop any non-finite residuals
  ok <- is.finite(resid)
  if (sum(ok) < 4L) {
    .log_warn("residual_morans_i(): fewer than 4 finite residuals.")
    return(NULL)
  }
  resid  <- resid[ok]
  coords <- coords[ok, , drop = FALSE]
  n      <- length(resid)

  # --- Weight matrix ---
  if (!is.null(weights)) {
    if (!is.matrix(weights) || nrow(weights) != n || ncol(weights) != n) {
      .log_warn("residual_morans_i(): user-supplied `weights` has wrong dimensions; falling back to kNN.")
      W <- .build_knn_weights(coords, k = k)
    } else {
      W <- weights
      # Validate row-standardisation: each row should sum to 1 (or 0 for
      # isolates).  If not, the Moran's I statistic is still computed
      # correctly for a general W via the Cliff & Ord formula, but the
      # interpretation differs from the row-standardised case.
      rs <- rowSums(W)
      non_zero <- rs[rs != 0]
      if (length(non_zero) > 0 && !isTRUE(all.equal(non_zero, rep(1, length(non_zero)), tolerance = 1e-8))) {
        .log_warn(
          paste0("residual_morans_i(): user-supplied `weights` is not row-standardised ",
                 "(row sums range from %.4g to %.4g). The Moran's I statistic and ",
                 "variance are computed using the general Cliff & Ord formula, which ",
                 "remains valid, but the magnitude of I is not directly comparable to ",
                 "results from row-standardised weights."),
          min(non_zero), max(non_zero)
        )
      }
    }
  } else {
    W <- .build_knn_weights(coords, k = k)
  }
  S0 <- sum(W)
  if (S0 < .Machine$double.eps || sum(resid^2) < .Machine$double.eps) {
    .log_warn("residual_morans_i(): degenerate weights or zero residual variance.")
    return(NULL)
  }

  # --- Moran's I ---
  resid_c <- resid - mean(resid)
  I <- (n / S0) * as.numeric(crossprod(resid_c, W %*% resid_c)) / sum(resid_c^2)

  # --- Analytical expectation & variance (randomisation assumption) ---
  EI <- -1 / (n - 1)

  # Cliff & Ord variance under randomisation
  # S1 = 0.5 * sum((W + t(W))^2) rewritten via the identity
  #    = sum(W^2) + sum(W * t(W))
  # so that sparse W never materialises the denser (W + t(W)) intermediate.
  Wt <- t(W)
  S1 <- sum(W * W) + sum(W * Wt)
  S2 <- sum((rowSums(W) + colSums(W))^2)
  m2 <- sum(resid_c^2) / n
  m4 <- sum(resid_c^4) / n
  b2 <- m4 / (m2^2)                          # kurtosis

  A  <- n * ((n^2 - 3 * n + 3) * S1 - n * S2 + 3 * S0^2)
  D  <- n * (n - 1) * (n - 2) * (n - 3) * S0^2
  C  <- (n^2 - n) * S1 - 2 * n * S2 + 6 * S0^2
  VI <- (A - b2 * C) / D - EI^2

  sd_I <- if (VI > 0) sqrt(VI) else NA_real_
  z    <- if (is.finite(sd_I) && sd_I > 0) (I - EI) / sd_I else NA_real_

  p <- if (is.finite(z)) {
    switch(alternative,
      two.sided = 2 * stats::pnorm(abs(z), lower.tail = FALSE),
      greater   = stats::pnorm(z, lower.tail = FALSE),
      less      = stats::pnorm(z, lower.tail = TRUE)
    )
  } else {
    NA_real_
  }

  list(observed = I, expected = EI, sd = sd_I, z = z,
       p_value = p, n = n)
}


#' Compute residual Moran's I for every model in a named list of fits
#'
#' Convenience wrapper that calls \code{residual_morans_i()} on each
#' element and returns a data.frame suitable for joining to the metrics
#' table produced by \code{compare_models()}.
#'
#' @param fits Named list of \code{spatial_fit} objects.
#' @return A data.frame with columns \code{model}, \code{resid_morans_I},
#'   \code{resid_morans_z}, and \code{resid_morans_p}.
#' @keywords internal
#' @noRd
.residual_morans_table <- function(fits) {
  rows <- lapply(names(fits), function(nm) {
    mi <- residual_morans_i(fits[[nm]])
    if (is.null(mi)) {
      data.frame(model = nm, resid_morans_I = NA_real_,
                 resid_morans_z = NA_real_, resid_morans_p = NA_real_,
                 stringsAsFactors = FALSE)
    } else {
      data.frame(model = nm, resid_morans_I = mi$observed,
                 resid_morans_z = mi$z, resid_morans_p = mi$p_value,
                 stringsAsFactors = FALSE)
    }
  })
  do.call(rbind, rows)
}


# ---------------------------------------------------------------------------
# evaluate_insample: metrics from already-fit spatial_fit objects
# ---------------------------------------------------------------------------

#' Compute in-sample (or out-of-sample) metrics for fitted spatial models
#'
#' Accepts a single \code{spatial_fit} object or a named list of them.
#' Does NOT refit — uses \code{fitted()} for in-sample and
#' \code{predict()} for new data.
#'
#' @param fits A \code{spatial_fit} object, or a named list of them
#'   (e.g. \code{list(GWR = gwr_obj, Bayesian = bayes_obj)}).
#' @param newdata Optional sf object for out-of-sample evaluation.
#'   Must contain the response variable and all predictors.
#'   If NULL, in-sample metrics are computed.
#' @param ... Extra arguments passed to predict().
#' @return A data.frame with one row per model and columns for
#'   model name and all regression metrics.
#' @export
evaluate_insample <- function(fits, newdata = NULL, ...) {
  # Accept a single fit

  if (inherits(fits, "spatial_fit")) {
    fits <- stats::setNames(list(fits), class(fits)[1L])
  }
  if (!is.list(fits) || length(fits) == 0L)
    stop("evaluate_insample(): `fits` must be a spatial_fit or a named list of them.")

  rows <- lapply(names(fits), function(nm) {
    obj <- fits[[nm]]
    if (!inherits(obj, "spatial_fit")) {
      .log_warn("evaluate_insample(): '%s' is not a spatial_fit; skipping.", nm)
      return(NULL)
    }
    met <- model_metrics(obj, newdata = newdata, ...)
    cbind(data.frame(model = nm, stringsAsFactors = FALSE), met)
  })

  do.call(rbind, Filter(Negate(is.null), rows))
}


# ---------------------------------------------------------------------------
# compare_models: pure comparison from a list of fit objects
# ---------------------------------------------------------------------------

#' Side-by-side comparison of fitted spatial models
#'
#' Takes a named list of already-fit \code{spatial_fit} objects and produces
#' a tidy comparison table including in-sample metrics and model-specific
#' information criteria (AICc, LOOIC).
#'
#' @param fits A named list of \code{spatial_fit} objects.
#' @param newdata Optional sf for out-of-sample evaluation.
#' @param ... Extra arguments passed to predict().
#' @return A data.frame comparing all models.
#' @export
compare_models <- function(fits, newdata = NULL, ...) {
  if (!is.list(fits) || length(fits) == 0L)
    stop("compare_models(): `fits` must be a non-empty named list of spatial_fit objects.")

  met_df <- evaluate_insample(fits, newdata = newdata, ...)

  # Append model-specific information criteria
  met_df$AICc  <- NA_real_
  met_df$LOOIC <- NA_real_
  met_df$bandwidth_is_fallback <- NA
  for (i in seq_len(nrow(met_df))) {
    nm  <- met_df$model[i]
    obj <- fits[[nm]]
    if (inherits(obj, "gwr_fit")) {
      met_df$AICc[i] <- obj$info$AICc %||% NA_real_
      is_fb <- isTRUE(obj$info$bandwidth_is_fallback)
      met_df$bandwidth_is_fallback[i] <- is_fb
      if (is_fb) {
        warning(
          sprintf(
            "compare_models(): GWR model '%s' used a fallback bandwidth (%.4f). Metrics may be unreliable; consider re-fitting with an explicit bandwidth.",
            nm, obj$info$bandwidth
          ),
          call. = FALSE
        )
      }
    }
    if (inherits(obj, "bayesian_fit"))
      met_df$LOOIC[i] <- obj$info$looic %||% NA_real_
  }

  # --- Post-fit residual spatial autocorrelation check ---
  moran_df <- .residual_morans_table(fits)
  met_df   <- merge(met_df, moran_df, by = "model", all.x = TRUE, sort = FALSE)

  # Emit warnings for models whose residuals still show significant

  # spatial autocorrelation (alpha = 0.05)
  for (i in seq_len(nrow(met_df))) {
    p_val <- met_df$resid_morans_p[i]
    I_val <- met_df$resid_morans_I[i]
    if (is.finite(p_val) && p_val < 0.05) {
      .log_warn(
        paste0("compare_models(): residuals of '%s' show significant ",
               "spatial autocorrelation (Moran's I = %.4f, p = %.4g). ",
               "The model may not fully capture the spatial structure."),
        met_df$model[i], I_val, p_val
      )
    }
  }

  met_df
}


# ---------------------------------------------------------------------------
# compare_models_cv: cross-validated comparison
# ---------------------------------------------------------------------------

#' Cross-validated comparison of spatial models
#'
#' Fits and cross-validates one or more model types, returning a unified
#' comparison table.  Unlike \code{compare_models()}, this function does
#' perform fitting (inside CV folds), because CV inherently requires
#' repeated fitting.
#'
#' @param data_sf An sf object.
#' @param response_var Response column name.
#' @param predictor_vars Predictor column names.
#' @param models Character vector: subset of c("GWR", "Bayesian").
#' @param k Number of folds. Default 5.
#' @param seed RNG seed. Default 123.
#' @param folds Optional precomputed fold splits.
#' @param boundary Optional polygon sf/sfc.
#' @param pointize Geometry coercion strategy.
#' @param gwr_args Extra arguments for fit_gwr_model() / cv_gwr().
#' @param bayes_args Extra arguments for fit_bayesian_spatial_model().
#' @param summary "mean" or "median" for Bayesian predictions.
#' @param quiet Logical; suppress messages.
#' @return A list with overall, by_fold, and per-model cv_results.
#' @export
compare_models_cv <- function(
    data_sf, response_var, predictor_vars,
    models = c("GWR", "Bayesian"),
    k = 5, seed = 123, folds = NULL, boundary = NULL, pointize = "auto",
    gwr_args = list(), bayes_args = list(),
    summary = c("mean", "median"),
    quiet = FALSE
) {
  summary <- match.arg(summary)
  .msg <- function(...) if (!quiet) message(...)

  if (!inherits(data_sf, "sf"))
    stop("compare_models_cv(): `data_sf` must be an sf object.")
  models <- unique(intersect(models, c("GWR", "Bayesian")))
  if (length(models) == 0L) models <- "GWR"

  for (m in models) {
    if (!.model_available(m)) {
      .msg(sprintf("compare_models_cv(): dropping %s (package/function unavailable).", m))
      models <- setdiff(models, m)
    }
  }
  if (length(models) == 0L) stop("compare_models_cv(): no viable models.")

  if (!("..row_id" %in% names(data_sf)))
    data_sf$`..row_id` <- seq_len(nrow(data_sf))

  cleanup <- .with_seed(seed)
  on.exit(cleanup(), add = TRUE)

  comparison_rows <- list(); by_fold_rows <- list(); cv_results <- list()

  if ("GWR" %in% models) {
    .msg("compare_models_cv(): running CV for GWR ...")
    base <- c(
      list(data_sf = data_sf, response_var = response_var,
           predictor_vars = predictor_vars, folds = folds,
           k = k, seed = seed, boundary = boundary, pointize = pointize),
      gwr_args
    )
    gwr_cv <- do.call(cv_gwr, .filter_args(cv_gwr, base))
    cv_results$gwr_cv <- gwr_cv
    ov <- try(as.data.frame(gwr_cv$overall), silent = TRUE)
    if (inherits(ov, "try-error") || nrow(ov) == 0L)
      ov <- data.frame(RMSE = NA_real_, MAE = NA_real_, MAPE = NA_real_, SMAPE = NA_real_,
                       R2 = NA_real_, Adj_R2 = NA_real_)
    ov$model <- "GWR"
    comparison_rows[["GWR"]] <- ov
    bf <- try(as.data.frame(gwr_cv$fold_metrics), silent = TRUE)
    if (!inherits(bf, "try-error") && nrow(bf)) {
      bf$model <- "GWR"; by_fold_rows[["GWR"]] <- bf
    }
  }

  if ("Bayesian" %in% models) {
    .msg("compare_models_cv(): running CV for Bayesian ...")
    base <- list(data_sf = data_sf, response_var = response_var,
                 predictor_vars = predictor_vars, folds = folds,
                 k = k, seed = seed, boundary = boundary,
                 pointize = pointize, summary = summary,
                 fit_args = bayes_args)
    bayes_cv <- do.call(cv_bayes, .filter_args(cv_bayes, base))
    cv_results$bayes_cv <- bayes_cv
    ov <- try(as.data.frame(bayes_cv$overall), silent = TRUE)
    if (inherits(ov, "try-error") || nrow(ov) == 0L)
      ov <- data.frame(RMSE = NA_real_, MAE = NA_real_, MAPE = NA_real_, SMAPE = NA_real_,
                       R2 = NA_real_, Adj_R2 = NA_real_)
    ov$model <- "Bayesian"
    comparison_rows[["Bayesian"]] <- ov
    bf <- try(as.data.frame(bayes_cv$fold_metrics), silent = TRUE)
    if (!inherits(bf, "try-error") && nrow(bf)) {
      bf$model <- "Bayesian"; by_fold_rows[["Bayesian"]] <- bf
    }
  }

  c(list(overall = dplyr::bind_rows(comparison_rows),
         by_fold = dplyr::bind_rows(by_fold_rows)),
    cv_results)
}


# ===========================================================================
# Legacy wrappers  (backward compatibility)
# ===========================================================================

#' Evaluate spatial models (legacy interface)
#'
#' Thin wrapper that preserves the original \code{evaluate_models()} call
#' signature.  New code should use \code{compare_models()} (for already-fit
#' objects) or \code{compare_models_cv()} (for CV) instead.
#'
#' @inheritParams compare_models_cv
#' @param do_cv Logical; use cross-validation. Default TRUE.
#' @param gwr_args,bayes_args Extra arguments for model fitting.
#' @return A list with CV or in-sample results.
#' @export
evaluate_models <- function(
    data_sf, response_var, predictor_vars,
    do_cv = TRUE, folds = NULL, k = 5, seed = 123,
    boundary = NULL, pointize = "auto", gwr_args = list(), bayes_args = list(),
    summary = c("mean", "median"), models = c("GWR", "Bayesian"), quiet = FALSE
) {
  summary <- match.arg(summary)
  .msg <- function(...) if (!quiet) message(...)

  if (!inherits(data_sf, "sf"))
    stop("evaluate_models(): `data_sf` must be an sf object.")
  models <- unique(intersect(models, c("GWR", "Bayesian")))
  if (length(models) == 0L) models <- "GWR"

  for (m in models) {
    if (!.model_available(m)) {
      .msg(sprintf("evaluate_models(): dropping %s (package/function unavailable).", m))
      models <- setdiff(models, m)
    }
  }
  if (length(models) == 0L) stop("evaluate_models(): no viable models.")

  if (!("..row_id" %in% names(data_sf)))
    data_sf$`..row_id` <- seq_len(nrow(data_sf))

  # ==== CV path — delegate to compare_models_cv ====
  if (isTRUE(do_cv)) {
    result <- compare_models_cv(
      data_sf = data_sf, response_var = response_var,
      predictor_vars = predictor_vars, models = models,
      k = k, seed = seed, folds = folds, boundary = boundary,
      pointize = pointize, gwr_args = gwr_args, bayes_args = bayes_args,
      summary = summary, quiet = quiet
    )
    # Reshape to match legacy return structure
    out <- list()
    if (!is.null(result$gwr_cv))   out$gwr_cv   <- result$gwr_cv
    if (!is.null(result$bayes_cv)) out$bayes_cv <- result$bayes_cv
    out$comparison <- result$overall
    return(out)
  }

  # ==== In-sample path — fit then evaluate ====
  dat_sf <- prep_model_data(data_sf, response_var, predictor_vars, boundary, pointize)

  fit_list <- list()

  if ("GWR" %in% models) {
    .msg("evaluate_models(): fitting GWR ...")
    fit_list$GWR <- do.call(fit_gwr_model,
      c(list(data_sf = dat_sf, response_var = response_var,
             predictor_vars = predictor_vars), gwr_args))
  }

  if ("Bayesian" %in% models) {
    .msg("evaluate_models(): fitting Bayesian ...")
    fit_list$Bayesian <- do.call(fit_bayesian_spatial_model,
      c(list(data_sf = dat_sf, response_var = response_var,
             predictor_vars = predictor_vars), bayes_args))
  }

  # Build legacy-shaped return
  ret <- list(
    formula = deparse(stats::reformulate(predictor_vars, response_var)),
    data    = dat_sf
  )

  if (!is.null(fit_list$GWR)) {
    gwr_obj <- fit_list$GWR
    ret$gwr <- list(
      fit       = gwr_obj,           # the spatial_fit object (new)
      model     = gwr_obj$engine,    # raw GWmodel result (legacy compat)
      bandwidth = gwr_obj$info$bandwidth,
      adaptive  = gwr_obj$info$adaptive,
      kernel    = gwr_obj$info$kernel,
      AICc      = gwr_obj$info$AICc %||% NA_real_
    )
  }

  if (!is.null(fit_list$Bayesian)) {
    bay_obj <- fit_list$Bayesian
    ret$bayes <- list(
      fit            = bay_obj,           # the spatial_fit object (new)
      model          = bay_obj$engine,    # raw brmsfit (legacy compat)
      coord_scaling  = bay_obj$info$coord_scaling,
      loo            = bay_obj$info$loo,
      looic          = bay_obj$info$looic %||% NA_real_
    )
  }

  ret$metrics <- compare_models(fit_list)

  # Attach per-model residual Moran's I diagnostics for programmatic access
  ret$residual_morans <- lapply(fit_list, residual_morans_i)

  ret
}


#' Cross-validated comparison with optional tessellation (legacy interface)
#'
#' Thin wrapper preserving the original \code{evaluate_models_cv()} call
#' signature.  New code should use \code{compare_models_cv()} directly.
#'
#' @inheritParams compare_models_cv
#' @param tess_method Tessellation type for diagnostics.
#' @param tess_args Extra arguments for tessellation builders.
#' @return A list with overall, by_fold, tessellation.
#' @export
evaluate_models_cv <- function(
    data_sf, response_var, predictor_vars,
    k = 5, seed = 123, folds = NULL, boundary = NULL, pointize = "auto",
    tess_method = c("grid", "hex", "square", "voronoi", "triangles"),
    tess_args = list(), summary = c("mean", "median"),
    models = c("GWR", "Bayesian"),
    gwr_args = list(), bayes_args = list(),
    quiet = FALSE
) {
  summary     <- match.arg(summary)
  tess_method <- match.arg(tess_method)

  # --- Build tessellation for diagnostics ---
  proj_pts <- try(ensure_projected(coerce_to_points(data_sf, pointize)), silent = TRUE)
  if (inherits(proj_pts, "try-error")) proj_pts <- ensure_projected(data_sf)

  boundary_proj <- if (!is.null(boundary)) {
    ensure_projected(boundary, proj_pts)
  } else {
    g    <- sf::st_geometry(proj_pts)
    hull <- sf::st_convex_hull(sf::st_union(g))
    bb   <- sf::st_bbox(hull)
    diag <- sqrt((bb$xmax - bb$xmin)^2 + (bb$ymax - bb$ymin)^2)
    sf::st_buffer(hull, dist = 0.01 * diag)
  }

  cells <- NULL
  if (tess_method %in% c("grid", "hex", "square")) {
    grid_type <- if (tess_method %in% c("hex", "square")) tess_method else "square"
    args <- modifyList(list(boundary = boundary_proj, target_cells = 30,
                            type = grid_type), tess_args)
    cells <- do.call(create_grid_polygons, args)
  } else if (tess_method %in% c("voronoi", "triangles")) {
    args <- modifyList(list(points_sf = coerce_to_points(proj_pts, "auto"),
                            boundary = boundary_proj, method = tess_method),
                       tess_args)
    cells <- do.call(build_tessellation, args)$cells
  }
  tessellation <- list(method = tess_method, args = tess_args, cells = cells)

  # --- Run CV models via the new compare_models_cv ---
  result <- compare_models_cv(
    data_sf = data_sf, response_var = response_var,
    predictor_vars = predictor_vars, models = models,
    k = k, seed = seed, folds = folds, boundary = boundary,
    pointize = pointize, gwr_args = gwr_args, bayes_args = bayes_args,
    summary = summary, quiet = quiet
  )

  c(list(overall = result$overall, by_fold = result$by_fold,
         tessellation = tessellation),
    result[grep("_cv$", names(result))])
}
