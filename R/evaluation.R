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


#' Merge caller-supplied extras onto compare_models_cv()'s base arguments
#'
#' `c(base, extra)` would produce two entries of the same name whenever an
#' extra collides with a base argument, and `do.call()` then fails with
#' "formal argument ... matched by multiple actual arguments".  A user passing
#' `rf_args = list(seed = 3)` hits this immediately, so extras replace base
#' entries by name instead of being appended alongside them.
#'
#' The four arguments that define *what is being compared* are protected: a
#' per-model override of the data, the response, the predictors or the folds
#' would silently make the models incomparable, which is the whole point of the
#' function.
#'
#' @param base Named list of arguments built by `compare_models_cv()`.
#' @param extra Named list supplied by the caller (`gwr_args` etc.).
#' @param arg_label Name of the caller-facing argument, for messages.
#' @return A named list with no duplicate names.
#' @keywords internal
#' @noRd
.merge_args <- function(base, extra, arg_label) {
  if (is.null(extra) || length(extra) == 0L) return(base)
  if (is.null(names(extra)) || any(!nzchar(names(extra))))
    stop(sprintf("compare_models_cv(): every element of `%s` must be named.",
                 arg_label), call. = FALSE)

  protected <- c("data_sf", "response_var", "predictor_vars", "folds")
  clash <- intersect(names(extra), protected)
  if (length(clash) > 0L) {
    warning(sprintf(
      paste0("compare_models_cv(): ignoring %s in `%s` -- these are set by ",
             "compare_models_cv() itself so that every model is scored on the ",
             "same data and the same folds."),
      paste(sQuote(clash), collapse = ", "), arg_label), call. = FALSE)
    extra <- extra[setdiff(names(extra), protected)]
  }
  utils::modifyList(base, extra)
}


#' Check whether a model backend is available
#' @keywords internal
#' @noRd
.model_available <- function(model_name) {
  switch(model_name,
    "GWR" = requireNamespace("GWmodel", quietly = TRUE) &&
      requireNamespace("sp", quietly = TRUE),
    "Bayesian" = requireNamespace("brms", quietly = TRUE),
    "RF" = requireNamespace("ranger", quietly = TRUE),
    FALSE
  )
}


#' Warn about argument-list entries a CV wrapper will silently discard
#'
#' \code{.filter_args()} drops anything the target function has no formal for
#' (and no \code{...} to absorb it), which is correct but silent: a mistyped or
#' misplaced entry then simply has no effect.  This names the casualties.
#'
#' @param fun The function the arguments are destined for.
#' @param args_list The user-supplied argument list.
#' @param arg_name Name of the parameter it arrived as, for the message.
#' @param fun_name Name of \code{fun}, for the message.
#' @return \code{NULL}, invisibly.  Called for the warning.
#' @keywords internal
#' @noRd
.warn_dropped_args <- function(fun, args_list, arg_name, fun_name) {
  if (!length(args_list)) return(invisible(NULL))
  fml <- try(formals(fun), silent = TRUE)
  if (inherits(fml, "try-error") || "..." %in% names(fml))
    return(invisible(NULL))
  dropped <- setdiff(names(args_list), names(fml))
  if (length(dropped) == 0L) return(invisible(NULL))
  warning(sprintf(
    paste0("compare_models_cv(): ignoring %s entr%s %s -- %s() has no such ",
           "argument. Only arguments of %s() are forwarded from `%s`; pass ",
           "anything else by calling %s() yourself."),
    arg_name, if (length(dropped) == 1L) "y" else "ies",
    paste(sQuote(dropped), collapse = ", "), fun_name, fun_name,
    arg_name, fun_name), call. = FALSE)
  invisible(NULL)
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
#' @param use_fnn Logical; use \pkg{FNN}'s kd-tree for the neighbour lookup.
#'   Defaults to whether \pkg{FNN} is installed.  Exposed so the dense
#'   brute-force fallback is reachable on a machine that has FNN.
#' @param use_matrix Logical; return a \pkg{Matrix} sparse matrix rather than a
#'   dense base matrix.  Defaults to whether \pkg{Matrix} is installed.
#'   Exposed for the same reason as \code{use_fnn}.  The fast sparse path
#'   requires both backends; with either disabled the dense path runs.
#' @return A row-standardised weight matrix W — sparse (dgCMatrix) when
#'   \pkg{Matrix} is available, otherwise a dense base matrix.
#' @keywords internal
#' @noRd
.build_knn_weights <- function(coords, k = 8L,
                              use_fnn    = requireNamespace("FNN",    quietly = TRUE),
                              use_matrix = requireNamespace("Matrix", quietly = TRUE)) {
  n <- nrow(coords)
  k <- min(as.integer(k), n - 1L)
  if (k < 1L) k <- 1L

  # Backend availability is a parameter rather than a direct requireNamespace()
  # call so the dense fallback -- and its size guard -- can be exercised on a
  # machine that has FNN installed.  Defaults preserve the previous behaviour.
  has_fnn    <- isTRUE(use_fnn)
  has_matrix <- isTRUE(use_matrix)

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
    # The guard has to test BOTH backends, not just FNN: this branch is entered
    # whenever either is missing and always allocates W <- matrix(0, n, n), so
    # keying it on FNN alone let an unbounded n x n allocation through whenever
    # FNN was present but Matrix was not (or use_matrix = FALSE).
    if (!(has_fnn && has_matrix) && n > 5000L)
      stop("n = ", n, " requires FNN for k-NN weights, and Matrix to hold them sparsely (the dense fallback would allocate an n*n matrix). Install both with install.packages(c(\"FNN\", \"Matrix\")).", call. = FALSE)
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
#' @param weights Optional user-supplied n x n weight matrix — a base
#'   matrix or a \pkg{Matrix}-package matrix (e.g. a sparse dgCMatrix).
#'   When \code{NULL} (the default), a k-nearest-neighbour binary weight matrix
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
#' @family model evaluation
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
#'   residual_morans_i(fit)  # z near 0 / p large = no residual structure
#' }
#' }
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
    # Accept base matrices AND Matrix-package sparse/dense matrices
    # (is.matrix() is FALSE for e.g. dgCMatrix, which previously caused a
    # misleading "wrong dimensions" warning and a silent kNN fallback).
    is_mat_like <- is.matrix(weights) || inherits(weights, "Matrix")
    if (!is_mat_like || nrow(weights) != n || ncol(weights) != n) {
      .log_warn("residual_morans_i(): user-supplied `weights` is not an n x n matrix (n = %d); falling back to kNN.", n)
      W <- .build_knn_weights(coords, k = k)
    } else {
      W <- weights
      # Validate row-standardisation: each row should sum to 1 (or 0 for
      # isolates).  If not, the Moran's I statistic is still computed
      # correctly for a general W via the Cliff & Ord formula, but the
      # interpretation differs from the row-standardised case.
      # (rowSums() does not S4-dispatch for Matrix classes from package
      # code when Matrix is loaded but not attached — same as the
      # variance block below.)
      rs <- if (inherits(W, "Matrix")) Matrix::rowSums(W) else rowSums(W)
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

  # Moran's I is a function of the CENTRED residuals, so the degeneracy guard
  # has to test their sum of squares rather than the raw one.  Constant
  # non-zero residuals pass sum(resid^2) but give resid_c == 0, and the
  # variance block below then computes m2 = 0 -> b2 = NaN -> VI = NaN, so
  # `if (VI > 0)` raises "missing value where TRUE/FALSE needed" instead of
  # this function returning NULL.
  resid_c <- resid - mean(resid)
  ss_c    <- sum(resid_c^2)
  if (S0 < .Machine$double.eps || ss_c < .Machine$double.eps) {
    .log_warn("residual_morans_i(): degenerate weights or zero residual variance.")
    return(NULL)
  }

  # --- Moran's I ---
  # sum(resid_c * (W %*% resid_c)) rather than crossprod(): see the dispatch
  # note in the variance block below.  The two are numerically identical.
  I <- (n / S0) * sum(resid_c * (W %*% resid_c)) / ss_c

  # --- Analytical expectation & variance (randomisation assumption) ---
  EI <- -1 / (n - 1)

  # Cliff & Ord variance under randomisation
  # S1 = 0.5 * sum((W + t(W))^2) rewritten via the identity
  #    = sum(W^2) + sum(W * t(W))
  # so that sparse W never materialises the denser (W + t(W)) intermediate.
  #
  # NOTE: t(), rowSums(), colSums() and crossprod() are plain base functions
  # that do NOT dispatch to Matrix's S4 methods when called from package code
  # with Matrix loaded-but-not-attached, so sparse W needs the Matrix::
  # generics explicitly (or, for crossprod, the primitive-only rewrite used
  # for I above).  Only primitives -- *, %*%, and sum() -- dispatch on their
  # own; crossprod is neither a primitive nor an internal generic, which is
  # why it belongs on this list and not with them.
  is_sparse <- inherits(W, "Matrix")
  Wt <- if (is_sparse) Matrix::t(W) else t(W)
  S1 <- sum(W * W) + sum(W * Wt)
  rs <- if (is_sparse) Matrix::rowSums(W) else rowSums(W)
  cs <- if (is_sparse) Matrix::colSums(W) else colSums(W)
  S2 <- sum((rs + cs)^2)
  m2 <- ss_c / n
  m4 <- sum(resid_c^4) / n
  b2 <- m4 / (m2^2)                          # kurtosis

  A  <- n * ((n^2 - 3 * n + 3) * S1 - n * S2 + 3 * S0^2)
  D  <- (n - 1) * (n - 2) * (n - 3) * S0^2
  C  <- (n^2 - n) * S1 - 2 * n * S2 + 6 * S0^2
  VI <- (A - b2 * C) / D - EI^2

  # is.finite() as well as > 0: a degenerate kurtosis (b2) can make VI NaN,
  # and `if (NaN > 0)` is an error rather than FALSE.
  sd_I <- if (is.finite(VI) && VI > 0) sqrt(VI) else NA_real_
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
#'   (e.g. \code{list(GWR = gwr_obj, Bayesian = bayes_obj)}).  The names are
#'   used as the model labels and every element must have one; an unnamed
#'   list is an error.
#' @param newdata Optional sf object for out-of-sample evaluation.
#'   Must contain the response variable and all predictors.
#'   If NULL, in-sample metrics are computed.
#' @param ... Extra arguments passed to predict().
#' @return A data.frame with one row per model and columns for
#'   model name and all regression metrics.
#' @family model evaluation
#' @export
evaluate_insample <- function(fits, newdata = NULL, ...) {
  # Accept a single fit

  if (inherits(fits, "spatial_fit")) {
    fits <- stats::setNames(list(fits), class(fits)[1L])
  }
  if (!is.list(fits) || length(fits) == 0L)
    stop("evaluate_insample(): `fits` must be a spatial_fit or a named list of them.")

  # The loop below is over names(fits), which is NULL for an unnamed list --
  # so lapply() would iterate zero times and this function would return NULL,
  # silently, having passed the check above.  compare_models() then turns that
  # NULL into a plain list and dies in seq_len(nrow(met_df)) with "argument
  # must be coercible to non-negative integer", nowhere near the real cause.
  nms <- names(fits)
  if (is.null(nms) || anyNA(nms) || !all(nzchar(nms)))
    stop("evaluate_insample(): `fits` must be a NAMED list of spatial_fit ",
         "objects -- the names label the models in the output. Supply them, ",
         "e.g. list(GWR = gwr_fit, Bayesian = bayes_fit).", call. = FALSE)

  rows <- lapply(nms, function(nm) {
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
#' @family model evaluation
#' @export
compare_models <- function(fits, newdata = NULL, ...) {
  # Wrap a bare fit exactly as evaluate_insample() does.  Without this, a
  # single spatial_fit passes the is.list() check below (a fit *is* a list),
  # and every downstream loop then iterates the fit's own components instead
  # of a set of models -- yielding all-NA Moran's I columns and one
  # "'fit' is not a spatial_fit object" log line per component.
  if (inherits(fits, "spatial_fit"))
    fits <- stats::setNames(list(fits), class(fits)[1L])
  if (!is.list(fits) || length(fits) == 0L)
    stop("compare_models(): `fits` must be a spatial_fit or a non-empty named list of them.")

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
#' @param models Character vector: any subset of \code{c("GWR", "Bayesian",
#'   "RF")}, in any order.  Each is cross-validated on the same \code{folds}.
#'   Names outside that set raise a warning and are dropped; if nothing
#'   recognised remains, this is an error rather than a silent fallback.
#'   A recognised model whose backend package is not installed is dropped with
#'   a message.
#' @param k Number of folds. Default 5.
#' @param seed RNG seed. Default 123.
#' @param folds Optional precomputed fold splits.
#' @param boundary Optional polygon sf/sfc.
#' @param pointize Geometry coercion strategy.
#' @param gwr_args Extra arguments for \code{\link{cv_gwr}}.  Only names that
#'   are formal arguments of \code{cv_gwr()} are forwarded --- it has no
#'   \code{...} --- so entries meant for \code{fit_gwr_model()} alone (e.g.
#'   \code{longlat}) cannot be passed this way.  Anything dropped is named in a
#'   warning; call \code{cv_gwr()} directly if you need it.
#' @param bayes_args Extra arguments for \code{fit_bayesian_spatial_model()}.
#'   Forwarded whole as \code{cv_bayes(fit_args = )}, so an unrecognised name
#'   is an error from \code{fit_bayesian_spatial_model()} rather than a silent
#'   drop.  \code{compute_loo}, \code{boundary} and \code{pointize} are
#'   overridden by the CV internals.
#' @param rf_args Extra arguments for \code{\link{cv_rf}}, which passes
#'   anything it does not recognise on to \code{\link{fit_rf_model}} and thence
#'   to \code{ranger::ranger()}.
#' @param summary "mean" or "median" for Bayesian predictions.
#' @param quiet Logical; suppress messages.
#' @return A list with overall, by_fold, and per-model cv_results
#'   (\code{gwr_cv}, \code{bayes_cv}, \code{rf_cv} for the models that ran).
#' @family model evaluation
#' @examples
#' \donttest{
#' if (requireNamespace("ranger", quietly = TRUE)) {
#'   library(sf)
#'   set.seed(1)
#'   n <- 120
#'   dat <- st_as_sf(
#'     data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), elev = rnorm(n)),
#'     coords = c("x", "y"), crs = 32632
#'   )
#'   dat$price <- 10 + 0.01 * st_coordinates(dat)[, 1] + 2 * dat$elev + rnorm(n)
#'   cmp <- compare_models_cv(dat, "price", "elev", models = "RF", k = 3,
#'                            rf_args = list(num_trees = 100))
#'   cmp$overall
#' }
#' }
#' @export
compare_models_cv <- function(
    data_sf, response_var, predictor_vars,
    models = c("GWR", "Bayesian"),
    k = 5, seed = 123, folds = NULL, boundary = NULL, pointize = "auto",
    gwr_args = list(), bayes_args = list(), rf_args = list(),
    summary = c("mean", "median"),
    quiet = FALSE
) {
  summary <- match.arg(summary)
  .msg <- function(...) if (!quiet) message(...)

  if (!inherits(data_sf, "sf"))
    stop("compare_models_cv(): `data_sf` must be an sf object.")

  # Unrecognised model names used to be dropped by a bare intersect(), so
  # models = "RF" silently ran GWR instead of what was asked for -- a wrong
  # answer presented as the right one.
  known    <- c("GWR", "Bayesian", "RF")
  models   <- unique(as.character(models))
  unknown  <- setdiff(models, known)
  models   <- intersect(models, known)
  if (length(unknown) > 0L)
    warning(sprintf(
      "compare_models_cv(): ignoring unknown model(s) %s. Supported: %s.",
      paste(sQuote(unknown), collapse = ", "),
      paste(sQuote(known), collapse = ", ")), call. = FALSE)
  if (length(models) == 0L)
    stop(sprintf("compare_models_cv(): no recognised model requested. Supported: %s.",
                 paste(sQuote(known), collapse = ", ")), call. = FALSE)

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
    .warn_dropped_args(cv_gwr, gwr_args, "gwr_args", "cv_gwr")
    base <- .merge_args(
      list(data_sf = data_sf, response_var = response_var,
           predictor_vars = predictor_vars, folds = folds,
           k = k, seed = seed, boundary = boundary, pointize = pointize),
      gwr_args, "gwr_args"
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

  if ("RF" %in% models) {
    .msg("compare_models_cv(): running CV for RF ...")
    # cv_rf() has `...` (forwarded to fit_rf_model() and on to ranger), so
    # nothing in rf_args is filtered out; `pointize` is deliberately not
    # passed, because cv_rf() has no such formal and it would land in ranger.
    base <- .merge_args(
      list(data_sf = data_sf, response_var = response_var,
           predictor_vars = predictor_vars, folds = folds,
           k = k, seed = seed, boundary = boundary),
      rf_args, "rf_args"
    )
    rf_cv <- do.call(cv_rf, base)
    cv_results$rf_cv <- rf_cv
    ov <- try(as.data.frame(rf_cv$overall), silent = TRUE)
    if (inherits(ov, "try-error") || nrow(ov) == 0L)
      ov <- data.frame(RMSE = NA_real_, MAE = NA_real_, MAPE = NA_real_, SMAPE = NA_real_,
                       R2 = NA_real_, Adj_R2 = NA_real_)
    ov$model <- "RF"
    comparison_rows[["RF"]] <- ov
    bf <- try(as.data.frame(rf_cv$fold_metrics), silent = TRUE)
    if (!inherits(bf, "try-error") && nrow(bf)) {
      bf$model <- "RF"; by_fold_rows[["RF"]] <- bf
    }
  }

  c(list(overall = dplyr::bind_rows(comparison_rows),
         by_fold = dplyr::bind_rows(by_fold_rows)),
    cv_results)
}
