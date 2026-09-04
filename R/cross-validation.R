# -----------------------------------------------------------------------------
# Internal: shared CV helpers
# -----------------------------------------------------------------------------

#' Safe tibble/data.frame constructor
#' @keywords internal
#' @noRd
.safe_tibble <- function(...) {
  if (requireNamespace("tibble", quietly = TRUE)) {
    tibble::tibble(...)
  } else {
    data.frame(..., stringsAsFactors = FALSE)
  }
}


#' Re-map pre-built fold indices to the subset of rows that survive prep
#'
#' @param folds List of \code{list(train, test)} using original row IDs.  A
#'   fold whose train and test sets intersect is an error: the model would be
#'   fitted and scored on the same rows and the result reported as a CV score.
#'   IDs present in no fold-eligible row are dropped and their count logged.
#' @param keep_idx Integer vector of surviving row IDs.
#' @param k Fallback fold count if folds is NULL.
#' @param seed RNG seed for random fold creation.
#' @return List of \code{list(train, test, fold_id)} with IDs restricted to
#'   \code{keep_idx}.  \code{fold_id} is the fold's index in the ORIGINAL
#'   \code{folds} object, carried through so that dropping an unusable fold
#'   does not renumber the survivors: downstream \code{fold} columns then still
#'   agree with \code{make_folds()$assignment$fold}.
#' @keywords internal
#' @noRd
.remap_folds <- function(folds, keep_idx, k = 5L, seed = 123L) {
  # keep_idx is the data's ..row_id.  A duplicate there means match() will
  # send every fold entry for that ID to its first row: the others are never
  # scored, and can sit in train and test at once.
  if (anyDuplicated(keep_idx))
    stop(sprintf(paste0("cross-validation: `..row_id` has %d duplicated ",
                        "value(s). Row IDs are what fold splits are made of, so ",
                        "every row needs its own; drop the column to have the ",
                        "rows numbered, or make the IDs unique."),
                 sum(duplicated(keep_idx))), call. = FALSE)
  if (is.null(folds)) {
    .log_warn(
      ".remap_folds(): no fold specification provided; falling back to random k-fold CV (k=%d). Random folds leak spatial autocorrelation and overstate out-of-sample performance.",
      k
    )
    warning(
      ".remap_folds(): falling back to random k-fold CV. For spatial data, use make_folds(method='block_kfold') to avoid optimistic performance estimates.",
      call. = FALSE
    )
    cleanup <- .with_seed(seed)
    on.exit(cleanup(), add = TRUE)
    assign_vec <- sample(rep(seq_len(k), length.out = length(keep_idx)))
    out <- vector("list", k)
    for (i in seq_len(k))
      out[[i]] <- list(train = keep_idx[assign_vec != i],
                       test  = keep_idx[assign_vec == i],
                       fold_id = i)
    return(out)
  }

  if (!is.null(folds$folds) && is.list(folds$folds)) {
    folds <- folds$folds
  }

  # Carry the ORIGINAL fold index alongside each split.  Dropping an unusable
  # fold below removes an element from the list, and without this every later
  # fold would be silently renumbered, so fold_metrics$fold and
  # predictions$fold would no longer line up with make_folds()$assignment$fold.
  # A hand-built `folds` list is documented as accepted by every cv_*(), and
  # two ways of getting it wrong went entirely unremarked.  Train and test
  # overlapping is not cross-validation at all -- the model is fitted and
  # scored on the same rows and the result is reported as a CV score (measured:
  # RMSE 0.50 against 0.97 for the same data properly split) -- so refuse it.
  # IDs that name no row were dropped silently by na.omit(match()), so a
  # typo'd or stale fold list looked like a clean run; those are counted and
  # reported, because prep_model_data() legitimately drops rows and the two
  # cases have to be told apart by the user, not by us.
  unknown_total <- 0L
  for (j in seq_along(folds)) {
    f <- folds[[j]]
    ov <- intersect(f$train, f$test)
    if (length(ov))
      stop(sprintf(paste0("cross-validation: fold %d has %d row ID(s) in BOTH ",
                          "its train and test sets (e.g. %s). A fold that ",
                          "trains on its own test rows is not a ",
                          "cross-validation split; rebuild the folds with ",
                          "make_folds()."),
                   j, length(ov),
                   paste(utils::head(format(ov), 3L), collapse = ", ")),
           call. = FALSE)
    unknown_total <- unknown_total +
      sum(is.na(match(c(f$train, f$test), keep_idx)))
  }
  if (unknown_total > 0L)
    .log_info(paste0(".remap_folds(): %d fold entr(y/ies) name row IDs that are ",
                     "not in the data and were dropped. This is expected when ",
                     "rows were removed for missing values; if it is not, the ",
                     "folds were built on different data."),
              unknown_total)

  remapped <- lapply(seq_along(folds), function(j) {
    f <- folds[[j]]
    list(
      train = keep_idx[stats::na.omit(match(f$train, keep_idx))],
      test  = keep_idx[stats::na.omit(match(f$test, keep_idx))],
      fold_id = j
    )
  })

  # An empty TRAINING set is just as fatal as an empty test set and was not
  # checked at all: buffered_loo with a buffer spanning the data, or
  # random_kfold on a single row, produces folds nothing can be fitted on.
  # They used to sail through here and get dropped one at a time inside
  # .cv_fit_one_fold(), so the only symptom was a generic "all folds failed"
  # warning at the very end.
  no_test  <- vapply(remapped, function(f) length(f$test) == 0L, logical(1))
  no_train <- vapply(remapped, function(f) length(f$train) < 2L, logical(1))
  if (any(no_test))
    .log_warn(".remap_folds(): %d fold(s) have empty test sets after remapping.",
              sum(no_test))
  if (any(no_train))
    .log_warn(".remap_folds(): %d fold(s) have fewer than 2 training rows after remapping and cannot be fitted.",
              sum(no_train))
  # A real warning as well as the log line: a fold that vanishes here never
  # reaches the fitter, and used to be invisible in n_folds_attempted -- a
  # five-fold set reported attempted = succeeded = 4 with no condition.
  n_drop <- sum(no_test | no_train)
  if (n_drop > 0L)
    warning(sprintf(paste0("cross-validation: %d of %d fold(s) dropped before ",
                           "fitting (empty test set, or fewer than 2 training ",
                           "rows, after incomplete rows were removed). ",
                           "n_folds_attempted counts the folds supplied, so ",
                           "compare it with n_folds_succeeded."),
                    n_drop, length(remapped)), call. = FALSE)

  remapped[!(no_test | no_train)]
}


#' Compute overall metrics from a predictions data.frame
#'
#' When the predictions data.frame contains a \code{y_train_mean} column
#' (one value per observation, equal to the training-fold mean for the fold
#' that produced that prediction), R² is computed against the correct
#' out-of-sample baseline: each observation's contribution to TSS uses the
#' training mean from its fold, not the pooled test-set mean.
#'
#' Adjusted R² is intentionally reported as \code{NA} for pooled CV metrics
#' because there are no well-defined degrees of freedom when predictions are
#' aggregated across independently fitted folds.  Per-fold Adj_R² is also
#' set to \code{NA} for GWR (whose effective parameter count is not the
#' global predictor count p) and for the Bayesian GP model (complex
#' effective degrees of freedom).  Adj_R² is only meaningful for models
#' with a fixed, global number of coefficients.
#'
#' @keywords internal
#' @noRd
.cv_overall_metrics <- function(preds) {
  ok <- is.finite(preds$y) & is.finite(preds$yhat)
  if (!any(ok))
    return(data.frame(RMSE = NA_real_, MAE = NA_real_, MAPE = NA_real_,
                      SMAPE = NA_real_, R2 = NA_real_, Adj_R2 = NA_real_,
                      n_pred = 0L))

  # Pass per-observation training-fold means when available so that
  # .compute_reg_metrics() uses the correct out-of-sample R² baseline.
  ytm <- if ("y_train_mean" %in% names(preds)) preds$y_train_mean else NULL

  # Intentionally omit `p` here (pass NULL) so that pooled Adj_R² is NA.

  # Adjusted R² requires well-defined degrees of freedom (n, p) from a
  # single model fit.  When predictions are pooled across CV folds, `n` is
  # the total number of held-out observations while each fold was fit
  # independently with its own training set — the resulting Adj_R² has no
  # clean statistical interpretation.  Per-fold Adj_R² (computed in
  # .cv_run_folds()) is valid and should be used instead.
  met <- .compute_reg_metrics(preds$y, preds$yhat, p = NULL, y_train_mean = ytm)
  data.frame(RMSE = met$RMSE, MAE = met$MAE, MAPE = met$MAPE, SMAPE = met$SMAPE,
             R2 = met$R2, Adj_R2 = met$Adj_R2, n_pred = met$n,
             stringsAsFactors = FALSE)
}


# -----------------------------------------------------------------------------
# Shared fold runner
# -----------------------------------------------------------------------------

#' Recover the message text from a \code{try-error} object
#'
#' \code{conditionMessage()} has no method for a \code{"try-error"} object, so
#' it must be reached through the condition hanging off the object as an
#' attribute; the deparsed string is the fallback for a \code{try-error}
#' carrying none.  Shared by the per-fold catch and the parallel-worker catch.
#'
#' @param e A \code{try-error} object.
#' @return Character(1), whitespace-trimmed and collapsed to one line.
#' @keywords internal
#' @noRd
.try_error_message <- function(e) {
  cond <- attr(e, "condition")
  txt <- if (is.null(cond)) paste(as.character(e), collapse = " ")
         else conditionMessage(cond)
  txt <- gsub("[\r\n]+", " ", paste(txt, collapse = " "))
  trimws(txt)
}


#' Render the first fold error as a sentence to append to a warning
#'
#' Turns \code{.cv_run_folds()}'s \code{fit_errors} into \code{" First error:
#' <msg>"}, or \code{""} when no fold reported one.  Without this, "all 5 folds
#' failed" is the whole diagnosis a user gets when the backend package is not
#' installed -- the word "brms" (or "GWmodel") never appears, even though
#' calling the fitter directly says so plainly.
#'
#' @param res The list returned by \code{.cv_run_folds()}.
#' @return Character(1); empty when there is nothing to add.
#' @keywords internal
#' @noRd
.cv_first_error_suffix <- function(res) {
  msgs <- res$fit_errors
  if (is.null(msgs) || length(msgs) == 0L) return("")
  msg <- msgs[[1L]]
  if (!is.character(msg) || is.na(msg) || !nzchar(msg)) return("")
  if (nchar(msg) > 300L) msg <- paste0(substr(msg, 1L, 300L), "...")
  paste0(" First error: ", msg)
}


#' A small, projection-invariant fingerprint of the rows a fold set describes
#'
#' Fold splits are lists of \code{..row_id} values, and row IDs are just
#' \code{seq_len(nrow())} unless the caller supplied them.  Handing
#' \code{cv_gwr()} a \code{folds} object built from a \emph{different} dataset
#' of the same size therefore "worked": every ID matched, every fold was
#' populated, and the model was scored on splits that describe other
#' observations entirely.  Nothing in the result said so.
#'
#' The probe stores the location of up to \code{max_probe} rows, spread evenly
#' over the rows, so the check can ask whether row 37 is still the same
#' feature.  Four properties matter, and the first three were each broken by
#' the previous version of this function:
#' \itemize{
#'   \item \strong{Taken on the caller's input, on both sides.}
#'     \code{make_folds()} probes the geometry it was handed, before
#'     pointization; the \code{cv_*()} wrappers probe the \code{data_sf} they
#'     were handed, before \code{prep_model_data()}.  Both therefore see the
#'     same features, and a polygon layer reduced with a different
#'     \code{pointize} in the two calls no longer reads as different data.
#'     Non-point geometry is reduced to its centroid here, identically on both
#'     sides, purely as a stable representative location.
#'   \item \strong{Row IDs kept in their own type.}  A user-supplied character
#'     \code{..row_id} used to be coerced with \code{as.integer()}, becoming
#'     all-\code{NA} and matching row 1 on every probe.
#'   \item \strong{Compared numerically, with a tolerance.}  Formatting the
#'     coordinates with \code{"\%.7g"} and comparing strings flipped on the
#'     ~1 in 5000 coordinates that a reprojection round trip moved across a
#'     rounding boundary (by 2.6e-9 degrees).  Locations are now compared
#'     within \code{1e-6} degrees (about 10 cm) after transformation to
#'     EPSG:4326, or within a relative \code{1e-9} for CRS-less coordinates.
#'   \item \strong{Tolerant of dropped rows.}  Per-row values, so the
#'     comparison is made over whichever probe rows survived the complete-case
#'     filter.
#' }
#'
#' @param x An sf object carrying a \code{..row_id} column.
#' @param max_probe Maximum number of rows to record.  64 keeps the folds
#'   object small while making a same-size different-dataset collision
#'   effectively impossible.
#' @return A list with \code{row_id} (the IDs as supplied), numeric \code{x}
#'   and \code{y}, and \code{lonlat} (whether the coordinates are in
#'   EPSG:4326, i.e. the input carried a CRS), or \code{NULL} when no probe
#'   can be taken.
#' @keywords internal
#' @noRd
.fold_row_probe <- function(x, max_probe = 64L) {
  tryCatch({
    if (!inherits(x, "sf") || !("..row_id" %in% names(x)) || nrow(x) == 0L)
      return(NULL)
    ids  <- x[["..row_id"]]
    take <- unique(round(seq(1, nrow(x), length.out = min(nrow(x), max_probe))))
    g    <- sf::st_geometry(x)[take]
    if (!all(sf::st_geometry_type(g, by_geometry = TRUE) == "POINT"))
      g <- suppressWarnings(sf::st_centroid(g))
    cr     <- suppressWarnings(sf::st_crs(x))
    lonlat <- FALSE
    if (!is.na(cr)) {
      g      <- suppressWarnings(sf::st_transform(g, 4326))
      lonlat <- TRUE
    }
    xy <- suppressWarnings(sf::st_coordinates(g))
    if (is.null(xy) || nrow(xy) != length(take)) return(NULL)
    list(row_id = ids[take], x = as.numeric(xy[, 1L]), y = as.numeric(xy[, 2L]),
         lonlat = lonlat)
  }, error = function(e) NULL)
}


#' Refuse fold splits that describe a different dataset
#'
#' Compares \code{folds$params$row_probe} against the data being
#' cross-validated -- the caller's own \code{data_sf}, before preparation --
#' over whichever probe rows are present.  A fold set built by an older
#' version of this package carries no probe and is passed through unchecked,
#' as is one whose IDs cannot be matched (\code{NA} IDs) or whose coordinate
#' space cannot be compared (one side carried a CRS and the other did not).
#'
#' @param folds A \code{make_folds()} return value, or \code{NULL}.
#' @param data_sf The sf being cross-validated, as the caller supplied it,
#'   with \code{..row_id} stamped.
#' @param caller Name used in the error message.
#' @return \code{invisible(NULL)}; called for the error.
#' @keywords internal
#' @noRd
.check_fold_probe <- function(folds, data_sf, caller) {
  probe <- tryCatch(folds$params$row_probe, error = function(e) NULL)
  if (is.null(probe) || is.null(probe$row_id) || length(probe$row_id) == 0L ||
      is.null(probe$x) || anyNA(probe$row_id))
    return(invisible(NULL))
  now <- .fold_row_probe(data_sf, max_probe = nrow(data_sf))
  if (is.null(now)) return(invisible(NULL))
  if (!identical(isTRUE(probe$lonlat), isTRUE(now$lonlat))) {
    .log_info(paste0("%s(): the supplied `folds` were built on data %s a CRS ",
                     "and this data %s one, so their locations cannot be ",
                     "compared; skipping the provenance check."),
              caller, if (isTRUE(probe$lonlat)) "with" else "without",
              if (isTRUE(now$lonlat)) "has" else "lacks")
    return(invisible(NULL))
  }

  m  <- match(probe$row_id, now$row_id)
  ok <- !is.na(m)
  if (!any(ok)) {
    stop(sprintf(paste0("%s(): none of the row IDs in `folds` are present in ",
                        "the data. The folds were built from a different ",
                        "dataset; rebuild them with make_folds() on this one."),
                 caller), call. = FALSE)
  }
  tol <- if (isTRUE(probe$lonlat)) 1e-6
         else max(1e-9 * max(abs(c(now$x, now$y)), na.rm = TRUE), 1e-9)
  dx  <- abs(probe$x[ok] - now$x[m[ok]])
  dy  <- abs(probe$y[ok] - now$y[m[ok]])
  # An empty or non-finite geometry gives NaN on both sides, and NaN <= tol is
  # NA -- so `bad` was NA and `if (bad > 0L)` died with R's internal "missing
  # value where TRUE/FALSE needed", turning make_folds() on a layer with one
  # bad geometry (which make_folds itself drops and reports) into a hard error
  # on every subsequent cv_*() call.  Such a row cannot be compared, so it is
  # not evidence of a different dataset: exclude it, exactly as a row whose ID
  # is absent is excluded.
  cmp <- is.finite(dx) & is.finite(dy)
  bad <- sum(cmp & !(dx <= tol & dy <= tol))
  ok[ok] <- cmp
  if (bad > 0L)
    stop(sprintf(paste0("%s(): the supplied `folds` were built from different ",
                        "data -- %d of %d checked row IDs sit at a different ",
                        "location here. Fold splits are lists of row IDs, so ",
                        "folds from another dataset of the same size apply ",
                        "silently and score the model on splits that describe ",
                        "other observations. Rebuild them with make_folds() on ",
                        "this data."),
                 caller, bad, sum(ok)), call. = FALSE)
  invisible(NULL)
}


#' Fit-predict a single CV fold
#'
#' Encapsulates the per-fold work so it can be called sequentially or in
#' parallel.  Returns \code{NULL} when the fold is unusable before any work
#' starts, or \code{list(error = <message>)} when the fit or the prediction
#' threw, so the caller can report the cause rather than only the count.
#'
#' @keywords internal
#' @noRd
.cv_fit_one_fold <- function(i, dat_sf, response_var, remapped_fold,
                             keep_idx, fit_one, fold_info_fn, predict_args,
                             p) {
  # The label to report: the fold's index in the ORIGINAL fold object, so that
  # dropping an unusable fold in .remap_folds() does not renumber the rest.
  # Falls back to the loop position for hand-built folds carrying no fold_id.
  fold_lab <- remapped_fold$fold_id %||% i

  tr_pos <- stats::na.omit(match(remapped_fold$train, keep_idx))
  te_pos <- stats::na.omit(match(remapped_fold$test, keep_idx))
  if (length(tr_pos) < 2L || length(te_pos) < 1L) return(NULL)

  train_sf <- dat_sf[tr_pos, , drop = FALSE]
  test_sf  <- dat_sf[te_pos, , drop = FALSE]

  # Drop geometry once per fold to avoid redundant copies on wide data frames
  train_df <- sf::st_drop_geometry(train_sf)
  test_df  <- sf::st_drop_geometry(test_sf)

  # Fit model on training fold.
  #
  # The error TEXT is returned, not just logged: when every fold fails for the
  # same reason -- overwhelmingly the "package 'brms'/'GWmodel' is required"
  # case -- the caller's "all N folds failed" warning is the only thing the
  # user sees, and without the cause it never names the missing backend.
  fit_obj <- try(fit_one(train_sf), silent = TRUE)
  if (inherits(fit_obj, "try-error")) {
    msg <- .try_error_message(fit_obj)
    .log_warn(".cv_run_folds(): fold %d fit failed; skipping. Cause: %s",
              fold_lab, msg)
    return(list(error = msg))
  }
  if (!inherits(fit_obj, "spatial_fit")) {
    msg <- sprintf("fit_fn() returned a %s, not a spatial_fit",
                   paste(class(fit_obj), collapse = "/"))
    .log_warn(".cv_run_folds(): fold %d did not return a spatial_fit; skipping.", fold_lab)
    return(list(error = msg))
  }

  # Predict on test fold via the S3 generic
  y_true <- test_df[[response_var]]
  y_hat  <- try(
    do.call(predict, c(list(object = fit_obj, newdata = test_sf), predict_args)),
    silent = TRUE
  )
  if (inherits(y_hat, "try-error") || !is.numeric(y_hat)) {
    msg <- if (inherits(y_hat, "try-error")) .try_error_message(y_hat) else
      sprintf("predict() returned a %s, not a numeric vector",
              paste(class(y_hat), collapse = "/"))
    .log_warn(".cv_run_folds(): fold %d predict failed; skipping. Cause: %s",
              fold_lab, msg)
    return(list(error = msg))
  }
  # cv_spatial() is documented as the extensibility point for arbitrary
  # learners, so a fit_fn whose predict() returns the wrong length is
  # user-reachable.  Nothing downstream notices: .compute_reg_metrics() and
  # the data.frame() below both RECYCLE silently, so a length-2 y_hat against
  # 4 test rows yields a 4-row frame with the predictions repeated and metrics
  # computed against fabricated pairs.
  if (length(y_hat) != length(y_true)) {
    .log_warn(".cv_run_folds(): fold %d predicted %d value(s) for %d test row(s); skipping.",
              fold_lab, length(y_hat), length(y_true))
    return(NULL)
  }

  # Training-set mean: the correct null-model baseline for out-of-sample
  # R².  Using the test-set mean instead would give the null model credit
  # for knowing information that was not available at prediction time,
  # systematically inflating CV R².
  y_train <- train_df[[response_var]]
  y_train_mean <- mean(y_train[is.finite(y_train)], na.rm = TRUE)

  met <- .compute_reg_metrics(y_true, y_hat, p = p,
                              y_train_mean = y_train_mean)
  if (met$n == 0L) return(NULL)

  # Base fold stats
  fs <- data.frame(
    fold = fold_lab, n_train = length(tr_pos), n_test = length(y_true),
    n_pred = met$n,
    RMSE = met$RMSE, MAE = met$MAE, MAPE = met$MAPE, SMAPE = met$SMAPE,
    R2 = met$R2, Adj_R2 = met$Adj_R2, stringsAsFactors = FALSE
  )

  # Append model-specific per-fold info (bandwidth, gp_k, CRPS, coverage …)
  if (!is.null(fold_info_fn)) {
    extra <- try(fold_info_fn(fit_obj, test_sf, y_true, y_hat), silent = TRUE)
    if (!inherits(extra, "try-error") && is.list(extra)) {
      for (cn in names(extra)) fs[[cn]] <- extra[[cn]]
    }
  }

  # Prediction rows — include the training-fold mean so that
  # .cv_overall_metrics() can compute the pooled R² with the
  # correct per-observation baseline.
  pr <- data.frame(
    `..row_id` = test_sf$`..row_id`, fold = fold_lab,
    y = as.numeric(y_true), yhat = as.numeric(y_hat),
    y_train_mean = y_train_mean,
    stringsAsFactors = FALSE
  )

  list(pred_row = pr, fold_stat = fs)
}


#' Empirical CRPS via the energy (NRG) form for equally weighted draws
#'
#' Computes the continuous ranked probability score per observation:
#' CRPS(F, y) = E|X - y| - 0.5 * E|X - X'|, where X, X' are independent
#' draws from the predictive distribution.  The spread term uses the
#' Gini-mean-difference identity on sorted draws:
#' 0.5 * E|X - X'| = (1/m^2) * sum_i x_(i) * (2i - m - 1)
#' for m equally weighted draws (Gneiting & Raftery 2007, Eq. 21; see also
#' the scoringRules package documentation).
#'
#' @param draws Numeric matrix (m draws x n observations) of posterior
#'   predictive draws.
#' @param y Numeric vector of length n of observed values.
#' @return Numeric vector of length n with per-observation CRPS.
#' @keywords internal
#' @noRd
.crps_energy <- function(draws, y) {
  m <- nrow(draws)
  mae_terms <- colMeans(abs(sweep(draws, 2L, y)))
  # Vectorised within-column sort in base R: order by column, then by value.
  sorted_draws <- matrix(draws[order(col(draws), draws)], nrow = m)
  weights <- (2 * seq_len(m) - m - 1) / (m * m)
  spread_terms <- colSums(sorted_draws * weights)
  mae_terms - spread_terms
}


#' Resolve parallel settings into a usable core count
#'
#' Returns 1L for sequential execution or an integer >= 2 for parallel.
#' On Windows, \code{parallel::mclapply()} falls back to serial, so we
#' warn and return 1L.
#'
#' @param parallel Logical or positive integer.  \code{TRUE} auto-detects,
#'   an integer sets the core count explicitly, \code{FALSE} is sequential.
#' @param n_cores Deprecated alias kept for backwards compatibility.
#'   Overrides \code{parallel} when not \code{NULL}.
#' @return Integer number of worker cores to use (1 = sequential).
#' @keywords internal
#' @noRd
.resolve_n_cores <- function(parallel = FALSE, n_cores = NULL) {
  if (!is.null(n_cores)) {
    eff <- .sanitize_core_count(n_cores)
  } else if (isTRUE(parallel)) {
    # detectCores() can return NA on some platforms; sanitize before use.
    eff <- .sanitize_core_count(parallel::detectCores(logical = FALSE) - 1L)
  } else if (is.numeric(parallel) && length(parallel) == 1L &&
             !is.na(parallel) && parallel > 1) {
    eff <- .sanitize_core_count(parallel)
  } else {
    return(1L)
  }
  if (.Platform$OS.type == "windows" && eff > 1L) {
    message("cv parallel: forked parallelism (mclapply) is not available on Windows. ",
            "Falling back to sequential execution. For Windows parallelism, consider ",
            "the 'future' / 'future.apply' packages.")
    return(1L)
  }
  eff
}


#' Run a fit-predict loop across CV folds
#'
#' Fits a model on each training fold and predicts onto the test fold
#' using the standardised \code{predict()} S3 method.  This is the
#' single implementation that both \code{cv_gwr()} and \code{cv_bayes()}
#' delegate to.
#'
#' When \code{parallel = TRUE} (or an integer > 1), folds are fitted in
#' parallel using \code{parallel::mclapply()}, which yields near-linear
#' speedup on macOS and Linux.  On Windows, forked parallelism is not
#' available and execution falls back to sequential with a message.
#'
#' @param dat_sf Prepared sf data (projected, clean).
#' @param response_var Character(1).
#' @param predictor_vars Character vector.
#' @param remapped_folds List of list(train, test).
#' @param keep_idx Integer vector of surviving row IDs.
#' @param fit_one A function(train_sf, ...) that returns a \code{spatial_fit}.
#' @param fold_info_fn Optional function(spatial_fit) -> named list of extra
#'   per-fold metadata columns (e.g. bandwidth, gp_k).
#' @param predict_args Named list of extra arguments for predict().
#' @param p Integer number of predictors for Adj R² (NULL to skip).
#'   Only meaningful when the model uses a fixed, global set of p
#'   coefficients (e.g. a linear model).  For models with spatially
#'   varying coefficients (GWR) or complex effective degrees of freedom
#'   (GP-based models), pass NULL so that per-fold Adj_R² is reported as
#'   NA rather than a misleadingly favourable value.
#' @param parallel Logical or positive integer.  If \code{TRUE},
#'   auto-detect the number of cores; if an integer > 1, use that many
#'   cores; if \code{FALSE} (default), run sequentially.
#' @param n_cores \emph{Deprecated.}
#'   Explicit core count; overrides \code{parallel} when set.
#' @param seed Integer RNG seed, or \code{NULL} to leave fold RNG unseeded.
#'   When supplied, one seed per fold is drawn in the parent process and
#'   applied inside the fold worker, so results depend on (seed, fold index)
#'   alone and \code{parallel = TRUE} reproduces \code{parallel = FALSE}
#'   exactly.
#' @return List with pred_rows and fold_stats.
#' @keywords internal
#' @noRd
.cv_run_folds <- function(dat_sf, response_var, predictor_vars,
                          remapped_folds, keep_idx, fit_one,
                          fold_info_fn = NULL, predict_args = list(),
                          p = NULL, parallel = FALSE, n_cores = NULL,
                          seed = NULL) {
  cores   <- .resolve_n_cores(parallel, n_cores)
  n_folds <- length(remapped_folds)

  # Draw one seed per fold in the parent so that each fold's RNG stream is a
  # function of (seed, fold index) only -- never of the execution path.  This
  # is what makes parallel and sequential runs bit-identical.  Forked workers
  # cannot guarantee that on their own: mclapply() seeds each child from the
  # current time and process ID unless the L'Ecuyer-CMRG generator is in use,
  # so without this the parallel path is irreproducible for any fit_fn that
  # consumes RNG (see cv_spatial(), which accepts an arbitrary learner).
  fold_seeds <- if (is.null(seed)) {
    rep(NA_integer_, n_folds)
  } else {
    cleanup_draw <- .with_seed(seed)
    on.exit(cleanup_draw(), add = TRUE)
    sample.int(.Machine$integer.max, n_folds)
  }

  fold_worker <- function(i) {
    # .with_seed() saves and restores .Random.seed, so seeding a fold never
    # leaks into the caller's RNG state.
    if (!is.na(fold_seeds[i])) {
      cleanup_i <- .with_seed(fold_seeds[i])
      on.exit(cleanup_i(), add = TRUE)
    }
    .cv_fit_one_fold(
      i = i, dat_sf = dat_sf, response_var = response_var,
      remapped_fold = remapped_folds[[i]], keep_idx = keep_idx,
      fit_one = fit_one, fold_info_fn = fold_info_fn,
      predict_args = predict_args, p = p
    )
  }

  if (cores > 1L) {
    message(sprintf("cv: running %d folds in parallel on %d cores.",
                    n_folds, cores))
    results <- parallel::mclapply(
      seq_along(remapped_folds), fold_worker, mc.cores = cores
    )
  } else {
    results <- lapply(seq_along(remapped_folds), fold_worker)
  }

  # mclapply() hands back a try-error OBJECT (not NULL) when a child errors or
  # is killed, and Negate(is.null) keeps it -- the subsequent [[ "pred_row" ]]
  # then fails with "subscript out of bounds", destroying the real diagnosis.
  failed <- vapply(results, function(z) inherits(z, "try-error"), logical(1))
  worker_msgs <- character(0)
  if (any(failed)) {
    worker_msgs <- vapply(results[failed], .try_error_message, character(1))
    .log_warn("cv: %d fold(s) failed in a parallel worker: %s",
              sum(failed), paste(unique(worker_msgs), collapse = "; "))
  }
  results <- results[!failed]
  results <- Filter(Negate(is.null), results)

  # Folds that threw come back as list(error = <message>) rather than NULL, so
  # the caller's "all N folds failed" warning can name the cause.  `$` (not
  # `[[`) because a successful fold's list has no "error" element and `[[`
  # would abort with "subscript out of bounds".
  is_err <- vapply(results, function(z) !is.null(z$error), logical(1))
  fit_errors <- c(worker_msgs,
                  unlist(lapply(results[is_err], `[[`, "error"),
                         use.names = FALSE))
  results <- results[!is_err]

  pred_rows  <- lapply(results, `[[`, "pred_row")
  fold_stats <- lapply(results, `[[`, "fold_stat")

  list(pred_rows = pred_rows, fold_stats = fold_stats,
       fit_errors = as.character(fit_errors))
}


# -----------------------------------------------------------------------------
# Spatial autocorrelation range estimation
# -----------------------------------------------------------------------------

#' Nearest-neighbour distance from each feature of \code{query} to \code{data}
#'
#' Uses \code{sf::st_nearest_feature()} so the search is indexed rather than a
#' dense cross-distance matrix, which would be prohibitive at the sizes NNDM is
#' guarded to.
#'
#' @param query,data sf layers in a common CRS.
#' @return Numeric vector, one distance per row of \code{query}.
#' @keywords internal
#' @noRd
.nn_dist_to <- function(query, data) {
  idx <- sf::st_nearest_feature(query, data)
  as.numeric(sf::st_distance(sf::st_geometry(query),
                             sf::st_geometry(data)[idx],
                             by_element = TRUE))
}


#' Estimate the spatial autocorrelation range from data
#'
#' Fits exponential (or spherical) variogram models and returns the
#' \emph{effective range} — the distance at which the semivariance reaches
#' ~95 \% of the sill.
#'
#' To guard against anisotropy, the function first estimates directional
#' variograms at 0° (N–S), 45°, 90° (E–W) and 135° azimuths, each with a
#' ±22.5° tolerance.  Those four windows tile all 180 distinct azimuths
#' exactly once, with no gap and no double-counted pair.
#'
#' An \strong{omnidirectional} variogram is fitted as well, and it is the
#' default answer.  Splitting 180° four ways leaves each directional variogram
#' about a quarter of the point pairs, and the maximum of four noisy estimates
#' is biased upward: on simulated \emph{isotropic} fields the max-of-four came
#' back about 40% above the truth, so blocks were sized 40% too wide and every
#' fold lost training data for no reason.  The all-pairs fit is the
#' best-powered estimate available, so the directional maximum is returned only
#' when anisotropy is \emph{established}: \strong{all four} directions must have
#' produced a usable fit, the ratio of the largest to the smallest must exceed
#' 1.5, \strong{and} the widest must stand more than 1.5× above the all-pairs
#' estimate.  (One further case: when the omnidirectional fit is itself
#' unusable, the directional sweep is all there is and its maximum is returned
#' whatever the ratio.  \code{anisotropy_used} records which answer you got.)
#' That is the
#' conservative choice where anisotropy is real (blocks must be at least as
#' large as the longest autocorrelation range to avoid leakage) without paying
#' for it where it is not.
#'
#' Sweeping only 0° and 90° would leave the azimuths between 23° and 67°, and
#' between 113° and 157°, covered by neither window: on simulated fields with a
#' 3:1 anisotropy and a true major-axis range of 300, a two-direction sweep
#' recovered 255 and 249 for major axes at 0° and 90° but only 151 and 147 at
#' 45° and 135°.  Since \code{make_folds(auto_range = TRUE)} sizes its blocks
#' from this number, that halved the blocks for a diagonally oriented field.
#'
#' A direction whose fit fails, does not converge, or reports a range beyond
#' the longest fitted lag is excluded and recorded as \code{NA} in the
#' \code{directional} attribute.  Fewer than two usable directions leaves the
#' omnidirectional fit as the only estimate.
#'
#' Every variogram model is fitted \strong{with a nugget}.  A nugget-free model
#' forces the curve through the origin, and on any real measurement (which has
#' one) \pkg{gstat}'s default N/h² weights buy that constraint by collapsing
#' the range: with a 50% nugget the fitted range came back at about 0.45 of the
#' truth, so \code{make_folds(auto_range = TRUE)} built blocks less than half
#' the correlation length it reported.
#'
#' A log warning is emitted when the directional maximum is used; where the
#' all-pairs estimate is available it names both the ratio and that estimate.  A
#' log note is emitted instead when the directional ranges vary but the spread
#' is consistent with sampling noise.
#'
#' The returned range is in the coordinate units of the (projected) data and
#' can be passed directly to \code{make_folds(block_size = ...)} to ensure
#' that CV blocks are at least as wide as the autocorrelation range.
#'
#' @param points_sf An sf object with point geometries (will be projected
#'   automatically if in geographic CRS).  Non-POINT geometry is reduced to
#'   representative points; any Z or M dimension is dropped, because
#'   \code{gstat::variogram()} uses every coordinate dimension and an XYZ layer
#'   would otherwise return a range in 3-D while every consumer of it works in
#'   2-D map distance; and rows with empty or non-finite coordinates are dropped
#'   with a logged count.
#' @param response_var Character(1) name of the response column.
#' @param predictor_vars Optional character vector.  When supplied, an OLS
#'   residual variogram is fitted instead of a raw-response variogram, which
#'   better reflects the autocorrelation that the spatial model must handle.
#' @param n_max Maximum number of points to subsample before fitting.
#'   Variogram estimation is O(n²) so this keeps runtime bounded.
#' @param cutoff Fraction of the maximum inter-point distance to use as
#'   the variogram lag cutoff.  Default 0.5.
#' @param range_frac Positive numeric.  A fitted range exceeding
#'   \code{range_frac * cutoff * max_dist} -- that is, beyond the longest lag
#'   the empirical variogram was actually fitted over -- is treated as
#'   unidentified and \code{NA_real_} is returned.
#'   \code{gstat::fit.variogram()} yields a finite number even when the
#'   variogram never reaches a sill, and such a value is extrapolation past the
#'   observed lags rather than a long autocorrelation range.  Passing it to
#'   \code{make_folds(auto_range = TRUE)} would collapse the block grid to a
#'   single block.  Default 1.0; raise it to accept ranges extrapolated beyond
#'   the fitted lags.
#' @param seed RNG seed for the \code{n_max} subsample, restored afterwards so
#'   the caller's random stream is untouched.  Default \code{123L}: the
#'   subsample is an internal approximation rather than part of the answer, and
#'   leaving it unseeded made the returned range differ between runs on
#'   identical input (19531, 19589, 19605 on three calls) and silently advanced
#'   the caller's RNG.  Pass \code{NULL} for the old unseeded behaviour, or a
#'   different number to check how sensitive the estimate is to the subsample.
#'   Ignored when \code{nrow(points_sf) <= n_max}, where nothing is sampled.
#' @return A single value of class \code{sac_range}, which behaves as an
#'   ordinary number.  There are three shapes, and they carry different
#'   attributes:
#'   \describe{
#'     \item{Success}{A positive effective range in projected coordinate units,
#'       with the fit attached as attributes \code{directional} (the 0°, 45°,
#'       90° and 135° ranges, named by azimuth), \code{anisotropy} (largest
#'       over smallest), \code{anisotropy_used} (logical: whether the returned
#'       range is the directional maximum rather than the all-pairs estimate),
#'       \code{crs} (the projected CRS the variogram was
#'       fitted in -- the unit of the range), \code{max_dist},
#'       \code{cutoff_dist}, \code{variogram} (the empirical variogram) and
#'       \code{variogram_model} (the fitted \code{gstat} model), so the fit can
#'       be inspected rather than trusted.}
#'     \item{Rejected range}{\code{NA_real_} when a range was fitted but exceeds
#'       \code{range_frac * cutoff * max_dist} and is therefore unidentified
#'       (see \code{range_frac}).  It is classed \code{sac_range} as well, so it
#'       prints as a bare \code{NA} rather than dumping its attributes, and it
#'       carries \code{max_dist}, \code{cutoff_dist}, \code{variogram} and
#'       \code{variogram_model} --- the evidence for the rejection --- plus
#'       \code{rejected_range} (the value that was refused) and
#'       \code{rejected_reason}, plus \code{crs} — so the units the rejected
#'       number was in stay recoverable, which is what
#'       \code{plot(type = "variogram")} labels its axis from.  It does
#'       \strong{not} carry \code{directional} or \code{anisotropy}.}
#'     \item{No fit}{A bare, attribute-less \code{NA_real_} when estimation
#'       could not be attempted or produced nothing at all (\pkg{gstat}
#'       missing, fewer than 30 finite values, a variable with no variance, a
#'       degenerate extent, or a singular variogram fit).}
#'   }
#'   Attributes and the class do not affect \code{is.na()} or
#'   \code{is.finite()}, so every downstream guard treats all three the same
#'   way it always did.
#' @family cross-validation
#' @examples
#' if (requireNamespace("gstat", quietly = TRUE)) {
#'   library(sf)
#'   set.seed(1)
#'   n <- 80
#'   xy <- data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000))
#'   xy$z <- sin(xy$x / 200) + rnorm(n, sd = 0.2)
#'   pts <- st_as_sf(xy, coords = c("x", "y"), crs = 32632)
#'   estimate_sac_range(pts, response_var = "z")
#' }
#' @export
estimate_sac_range <- function(points_sf, response_var,
                               predictor_vars = NULL,
                               n_max = 5000L, cutoff = 0.5,
                               range_frac = 1.0, seed = 123L) {
  if (!requireNamespace("gstat", quietly = TRUE)) {
    .log_warn("estimate_sac_range(): package 'gstat' is required for variogram estimation; returning NA.")
    return(NA_real_)
  }
  if (!inherits(points_sf, "sf"))
    stop("estimate_sac_range(): `points_sf` must be an sf object.", call. = FALSE)
  # MULTIPOINT must be coerced too, not merely admitted: st_coordinates()
  # returns one row per VERTEX, so any multi-vertex feature makes xy[i, ] a
  # different feature than row i, and every fold below misaligns silently.
  # coerce_to_points() takes centroids, matching prep_model_data().
  if (!all(sf::st_geometry_type(points_sf, by_geometry = TRUE) == "POINT"))
    points_sf <- coerce_to_points(points_sf, "auto")

  # gstat::variogram() uses EVERY coordinate dimension, so an XYZ layer had its
  # elevation folded into each lag and the returned "range" was a length in
  # 3-D -- while the block grid, the buffered-LOO buffer, nndm's neighbour
  # distances and summarize_by_cell(deff = "variogram") all work in 2-D map
  # distance.  Measured: 413.6 against 136.8 for the same stations with a
  # 0-2000 m elevation column.  Nothing warned.
  points_sf <- sf::st_zm(points_sf, drop = TRUE, what = "ZM")

  # A row with an empty or non-finite geometry makes gstat fail in EVERY
  # direction, so one bad point among 200 turned the whole layer's estimate
  # into NA under the message "variogram model fit failed" -- blaming the fit
  # rather than the row, and disagreeing with make_folds(), which drops such
  # rows and says how many.
  bad_geom <- sf::st_is_empty(points_sf)
  if (!all(bad_geom)) {
    xy_chk <- suppressWarnings(sf::st_coordinates(points_sf[!bad_geom, ]))
    if (nrow(xy_chk) == sum(!bad_geom) && ncol(xy_chk) >= 2L)
      bad_geom[!bad_geom] <- !stats::complete.cases(xy_chk[, 1:2, drop = FALSE])
  }
  if (any(bad_geom)) {
    .log_warn("estimate_sac_range(): dropping %d point(s) with empty or non-finite coordinates.",
              sum(bad_geom))
    points_sf <- points_sf[!bad_geom, , drop = FALSE]
  }

  pts <- ensure_projected(points_sf)

  # Subsample if large.
  #
  # `seed` defaults to a constant, NOT to NULL.  .with_seed(NULL) deliberately
  # neither seeds nor restores -- the right behaviour for a function whose
  # randomness is part of its answer -- but this subsample is not that: it is
  # an internal approximation to keep an O(n^2) variogram tractable, and its
  # only effects are to make the returned range irreproducible above n_max
  # (measured: 19531, 19589, 19605 on three calls with identical input) and to
  # advance the caller's RNG stream as a side effect of a function that looks
  # like a summary statistic.  make_folds(auto_range = TRUE) then sizes its
  # blocks from a number that changes between runs.
  n <- nrow(pts)
  if (n > n_max) {
    cleanup <- .with_seed(seed)
    on.exit(cleanup(), add = TRUE)
    idx <- sample.int(n, n_max)
    pts <- pts[idx, , drop = FALSE]
    n <- n_max
  }

  if (n < 30L) {
    .log_warn("estimate_sac_range(): fewer than 30 points; variogram estimate unreliable. Returning NA.")
    return(NA_real_)
  }

  # Build the variable to model: raw response or OLS residuals.
  #
  # Validate the column first.  as.numeric() on a factor returns its LEVEL
  # CODES -- 1, 2, 3 in whatever order the levels happen to sit -- so a factor
  # response produced a variogram of an arbitrary integer relabelling of the
  # categories, and the estimated range changed when the levels were reordered
  # (measured: 3700 against 2497 on the same data).  A character column becomes
  # all-NA and is caught only by the "too few finite values" guard downstream,
  # which blames the data rather than the column type.
  if (!is.character(response_var) || length(response_var) != 1L || is.na(response_var))
    stop("estimate_sac_range(): `response_var` must be a single column name.",
         call. = FALSE)
  pts_df <- sf::st_drop_geometry(pts)
  if (!(response_var %in% names(pts_df)))
    stop(sprintf("estimate_sac_range(): column '%s' not found in `points_sf`.",
                 response_var), call. = FALSE)
  y <- pts_df[[response_var]]
  if (!is.numeric(y)) {
    if (is.logical(y)) {
      y <- as.numeric(y)                 # 0/1 is a well-defined variogram target
    } else {
      stop(sprintf(paste0("estimate_sac_range(): response '%s' is %s, and a ",
                          "variogram needs a numeric variable. A factor's ",
                          "codes are an arbitrary relabelling of its levels ",
                          "-- fitting a variogram to them yields a range that ",
                          "changes when the levels are reordered. Encode the ",
                          "column numerically first."),
                   response_var,
                   if (is.factor(y)) "a factor" else sprintf("of class %s",
                     paste(class(y), collapse = "/"))),
           call. = FALSE)
    }
  }
  if (!is.null(predictor_vars) && length(predictor_vars) > 0L) {
    df <- pts_df
    ok_preds <- intersect(predictor_vars, names(df))
    if (length(ok_preds) > 0L) {
      fml <- stats::reformulate(ok_preds, response_var)
      lm_fit <- try(stats::lm(fml, data = df, na.action = stats::na.exclude),
                     silent = TRUE)
      if (!inherits(lm_fit, "try-error")) {
        resid <- stats::residuals(lm_fit)
        if (length(resid) != nrow(pts)) {
          .log_warn("estimate_sac_range(): OLS residual length (%d) does not match data rows (%d); using raw response.",
                    length(resid), nrow(pts))
        } else {
          y <- resid
        }
      }
    }
  }

  pts$..sac_var <- as.numeric(y)
  pts <- pts[is.finite(pts$..sac_var), , drop = FALSE]
  if (nrow(pts) < 30L) {
    .log_warn("estimate_sac_range(): too few finite values after filtering; returning NA.")
    return(NA_real_)
  }
  # A variable with no variance has no autocorrelation structure to estimate:
  # every semivariance is 0, and gstat's fit returned a finite "range" (168
  # for an exactly-explained response, 673 for a constant one) that was then
  # used to size blocks.  There is no range; say so.
  if (stats::sd(pts$..sac_var) < sqrt(.Machine$double.eps) *
        max(1, abs(mean(pts$..sac_var)))) {
    .log_warn(paste0("estimate_sac_range(): the variable being modelled is ",
                     "constant (zero variance%s), so it has no autocorrelation ",
                     "range. Returning NA."),
              if (!is.null(predictor_vars) && length(predictor_vars) > 0L)
                " -- the OLS residuals are all zero, so the predictors explain the response exactly"
              else "")
    return(NA_real_)
  }

  # Empirical variogram
  max_dist <- try({
    bb <- sf::st_bbox(pts)
    sqrt((bb["xmax"] - bb["xmin"])^2 + (bb["ymax"] - bb["ymin"])^2)
  }, silent = TRUE)

  if (inherits(max_dist, "try-error") || !is.finite(max_dist) || max_dist <= 0)
    return(NA_real_)

  cutoff_dist <- as.numeric(cutoff * max_dist)

  # --- Helper: fit a variogram model and return the effective range --------
  # gstat::fit.variogram() signals a failed fit by setting attr(., "singular")
  # to TRUE and emitting a warning -- it returns normally rather than throwing.
  # Testing only for a try-error therefore (a) made the "Sph" fallback below
  # unreachable, and (b) let a singular fit's `range` flow straight out as the
  # estimated autocorrelation range, which make_folds(auto_range = TRUE) then
  # sizes spatial blocks from.  The `max_supported` guard further down only
  # catches ranges beyond the fitted lags, not fitting artefacts inside them.
  # A NON-CONVERGED fit is unusable for the same reason a singular one is: the
  # returned `range` is wherever the optimiser happened to stop, not a fitted
  # parameter.  gstat signals it with a real R warning ("No convergence after
  # 200 iterations") and then returns normally, so try() -- which catches only
  # errors -- let both the warning escape to the user and the untrustworthy
  # range flow onward.
  #
  # This matters more since the sweep went to four directions: each variogram
  # gets about half the point pairs a two-direction sweep gave it, so a
  # direction failing to converge is now routine rather than exceptional, and
  # `plot(fit, type = "variogram")` started emitting a bare gstat warning on
  # ordinary data.  The caller already handles an unusable direction -- it is
  # excluded from the maximum, and the isotropic fallback (which pools every
  # pair, and therefore converges far more readily) takes over when fewer than
  # two directions survive.  So the answer is to refuse the fit, not to pass a
  # warning up about a fit nothing was going to use.
  #
  # Other warnings are muffled but LOGGED rather than dropped: they say
  # something about the data even when the fit is usable.
  .fit_one_vgm <- function(vg, model_type) {
    converged <- TRUE
    # WITH a nugget.  gstat::vgm(model = "Exp") alone is a nugget-free model,
    # which forces the fitted curve through the origin; on any real
    # measurement (which has one) gstat's default N/h^2 weights then buy that
    # constraint by collapsing the range.  Measured on simulated fields with a
    # 50% nugget the nugget-free fit returned about 0.45 of the true range, so
    # make_folds(auto_range = TRUE) built blocks less than half the
    # correlation length it reported.  A nugget model on the same empirical
    # variogram recovered the truth (ratio 1.09).  The partial-sill range is
    # the one to keep, and fit.variogram() returns the nugget as row 1 and the
    # structured component as row 2, which .vgm_range_of() already reads.
    m <- withCallingHandlers(
      try(gstat::fit.variogram(
            vg, gstat::vgm(psill = NA, model = model_type, range = NA,
                           nugget = NA)),
          silent = TRUE),
      warning = function(w) {
        msg <- conditionMessage(w)
        if (grepl("convergence", msg, ignore.case = TRUE)) {
          converged <<- FALSE
        } else {
          .log_info("estimate_sac_range(): gstat::fit.variogram(%s) warned: %s",
                    model_type, msg)
        }
        invokeRestart("muffleWarning")
      }
    )
    if (inherits(m, "try-error") || !is.data.frame(m)) return(NULL)
    if (isTRUE(attr(m, "singular"))) return(NULL)
    # NOT `return(NULL)` on non-convergence.  The range it carries must never
    # size a block -- that is what the `converged` flag below is for -- but the
    # MODEL and its empirical variogram are still the most useful thing a user
    # can look at, and a sill-less variogram is precisely the case worth
    # looking at.  Discarding it here made plot(fit, type = "variogram")
    # error out with "the residual variogram could not be fitted" on exactly
    # that input.  Mark it and let the caller decide.
    attr(m, "converged") <- converged
    m
  }

  .fit_vgm_range <- function(vg) {
    if (inherits(vg, "try-error") || !inherits(vg, "data.frame") || NROW(vg) < 3L) return(NA_real_)
    vgm_model <- .fit_one_vgm(vg, "Exp")
    if (is.null(vgm_model)) vgm_model <- .fit_one_vgm(vg, "Sph")
    if (is.null(vgm_model)) {
      # Both the exponential and the spherical fit were singular or errored:
      # there is no identified range, so say so rather than returning one.
      .log_warn("estimate_sac_range(): both the exponential and the spherical fit to this variogram are singular; no range is identified from it.")
      return(NA_real_)
    }
    spatial_rows <- vgm_model[vgm_model$model != "Nug", , drop = FALSE]
    if (nrow(spatial_rows) == 0L || all(!is.finite(spatial_rows$range)))
      return(NA_real_)
    raw_range  <- max(spatial_rows$range, na.rm = TRUE)
    model_type <- spatial_rows$model[which.max(spatial_rows$range)]
    eff <- if (identical(as.character(model_type), "Exp")) 3 * raw_range else raw_range
    if (!is.finite(eff) || eff <= 0) return(NA_real_)
    # Carry the fitted model out so callers can inspect the fit rather than
    # trust a bare number, and the convergence flag with it: a range from an
    # optimiser that stopped at its iteration limit is where it happened to
    # stop, not a fitted parameter, so it must not size a block -- while the
    # variogram behind it stays available to plot.
    structure(eff, vgm_model = vgm_model,
              converged = !identical(attr(vgm_model, "converged"), FALSE))
  }

  # The largest lag the empirical variogram is fitted over.  A fitted range
  # beyond it is extrapolation past every observed lag rather than a long
  # correlation length -- see the rejection block further down, which applies
  # the same bound to the final answer.  Hoisted here so each DIRECTION can be
  # tested against it too.
  max_supported <- range_frac * cutoff_dist

  # --- Directional variograms (0/45/90/135 deg, tolerance 22.5 deg) --------
  # gstat uses azimuth in degrees clockwise from north.  0 = N-S, 90 = E-W.
  #
  # FOUR directions, not two.  A +/-22.5 window around 0 and 90 covers
  # [337.5, 22.5] and [67.5, 112.5] -- exactly 90 of the 180 distinct azimuths.
  # Every direction between 23 and 67 degrees, and between 113 and 157, fell
  # into NEITHER window, so a field whose major axis lay there was measured by
  # two variograms that both cut across it.  Measured on simulated anisotropic
  # fields (ratio 3:1, true major-axis range 300): the estimate came back at
  # 255 and 249 for major axes at 0 and 90 degrees, but 151 and 147 at 45 and
  # 135 -- half the true range.  make_folds(auto_range = TRUE) sizes its blocks
  # from this number, so a diagonally-oriented field silently got blocks half
  # as wide as the correlation it was meant to separate, which is the exact
  # leakage blocked CV exists to prevent.
  #
  # c(0, 45, 90, 135) at +/-22.5 tiles all 180 azimuths with no overlap and no
  # gap.  The cost is that each variogram uses about half as many pairs as a
  # 2-direction sweep would, which is why the isotropic fallback below matters:
  # on sparse data one of the four is more likely to fail to fit.
  dir_az   <- c(0, 45, 90, 135)
  dir_fits <- lapply(dir_az, function(az) {
    vg_dir <- try(
      gstat::variogram(..sac_var ~ 1, data = pts,
                       cutoff = cutoff_dist,
                       alpha = az, tol.hor = 22.5),
      silent = TRUE
    )
    list(vg = vg_dir, fit = .fit_vgm_range(vg_dir))
  })
  dir_ranges <- vapply(dir_fits, function(f) as.numeric(f$fit), numeric(1))

  # A direction counts only if its fit is BOTH finite and identified, i.e.
  # within the longest lag the variogram was fitted over.  The global check
  # further down applies exactly this test to the final answer; applying it per
  # direction as well is what keeps one starved direction from deciding the
  # result.  Splitting 180 degrees four ways leaves each variogram about half
  # the point pairs a two-direction sweep would give it, so on small or
  # irregular samples one direction can come back with a fit that never reaches
  # a sill -- observed at 15908 and 24982 against 82-196 for the other three.
  # Taking the max of that is not "conservative", it is reading a failed fit as
  # a long correlation length, and the global guard then discards the whole
  # estimate even though three directions agreed.
  dir_conv <- vapply(dir_fits,
                     function(f) !identical(attr(f$fit, "converged"), FALSE),
                     logical(1))
  dir_ok      <- is.finite(dir_ranges) & dir_ranges <= max_supported & dir_conv
  # Two usable directions are enough to say something about anisotropy.  Below
  # that, the isotropic variogram -- which pools every pair and is therefore
  # the stable estimate -- is the honest fallback.
  dir_success <- sum(dir_ok) >= 2L

  # Defined here, not inside the branch below: the success return reports it as
  # the `directional` attribute on BOTH paths, and the isotropic fallback never
  # enters that branch.
  usable          <- dir_ranges
  usable[!dir_ok] <- NA_real_
  anisotropy    <- NA_real_
  aniso_used    <- FALSE
  fit_converged <- TRUE
  vg_used       <- NULL
  vgm_used    <- NULL

  # --- Select the effective range ------------------------------------------
  # The ISOTROPIC variogram is fitted unconditionally, not only as a fallback.
  # It pools every point pair, so it is the best-powered estimate available,
  # and it is the yardstick the directional sweep has to beat before its
  # maximum is believed: splitting 180 degrees four ways leaves each
  # directional variogram about a quarter of the pairs, and the maximum of
  # four noisy estimates is biased upward.  Measured on ISOTROPIC simulated
  # fields, the max-of-four came in about 40% above the truth and the
  # "notable anisotropy" warning fired on the majority of them -- so the
  # warning was mostly false alarms and the number it reported was mostly
  # sampling error.  Blocks sized 40% too wide are not "conservative", they
  # throw away training data in every fold.
  vg_iso_always  <- try(
    gstat::variogram(..sac_var ~ 1, data = pts, cutoff = cutoff_dist),
    silent = TRUE
  )
  iso_fit_always <- .fit_vgm_range(vg_iso_always)
  iso_ok <- is.finite(iso_fit_always) &&
            as.numeric(iso_fit_always) <= max_supported &&
            !identical(attr(iso_fit_always, "converged"), FALSE)

  if (dir_success) {
    dir_max    <- max(usable, na.rm = TRUE)
    anisotropy <- dir_max / min(usable, na.rm = TRUE)
    winner     <- which.max(usable)
    if (any(!dir_ok))
      .log_info(paste0("estimate_sac_range(): %d of %d directional variograms ",
                       "did not yield an identified range (%s) and were ",
                       "excluded."),
                sum(!dir_ok), length(dir_az),
                paste(sprintf("%d\u00b0", dir_az[!dir_ok]), collapse = ", "))

    # Anisotropy has to clear TWO hurdles before the directional maximum is
    # used: a ratio the sweep itself considers notable, and a maximum that
    # stands clearly above the all-pairs estimate.  On an isotropic field the
    # second hurdle is what the noise cannot clear, because the isotropic fit
    # is centred on the truth while the directional maximum scatters around
    # it.
    # Three hurdles, and the first is the one the noise cannot fake: ALL FOUR
    # directions must have produced a usable fit.  Geometric anisotropy is a
    # smooth function of azimuth, so it shows up in every direction; a sweep in
    # which two directions failed to reach a sill has no azimuthal pattern to
    # report, only two noisy numbers, and taking their ratio as evidence is how
    # a field with a true range of 80 came back at 248 (the all-pairs fit said
    # 84).  Then the usual two: a ratio the sweep considers notable, and a
    # maximum that stands clearly above the all-pairs estimate.
    aniso_established <- sum(dir_ok) == length(dir_az) &&
      is.finite(anisotropy) && anisotropy > 1.5 &&
      (!iso_ok || dir_max > 1.5 * as.numeric(iso_fit_always))

    if (aniso_established) {
      aniso_used      <- TRUE
      effective_range <- dir_max
      vg_used  <- dir_fits[[winner]]$vg
      vgm_used <- attr(dir_fits[[winner]]$fit, "vgm_model")
      .log_warn(
        paste0("estimate_sac_range(): directional ranges vary by a factor of ",
               "%.1f and the widest, %.1f, stands above the all-directions ",
               "estimate of %s, so the maximum is used -- the conservative ",
               "choice for sizing blocks. Directional ranges: %s. Each ",
               "direction sees about a quarter of the point pairs, so at ",
               "modest sample sizes part of this spread is sampling noise; ",
               "pass an explicit block_size if you know the field's anisotropy."),
        anisotropy, dir_max,
        if (iso_ok) sprintf("%.1f", as.numeric(iso_fit_always)) else "n/a",
        paste(sprintf("%d\u00b0 = %.1f", dir_az[dir_ok], dir_ranges[dir_ok]),
              collapse = ", ")
      )
    } else if (iso_ok) {
      effective_range <- as.numeric(iso_fit_always)
      vg_used  <- vg_iso_always
      vgm_used <- attr(iso_fit_always, "vgm_model")
      if (is.finite(anisotropy) && anisotropy > 1.5)
        .log_info(
          paste0("estimate_sac_range(): the directional ranges vary by a factor ",
                 "of %.1f (%s), but the widest does not stand above the ",
                 "all-directions estimate (%.1f), so the spread is consistent ",
                 "with sampling noise -- each direction sees about a quarter of ",
                 "the point pairs. Using the all-directions estimate. Pass an ",
                 "explicit block_size if you know the field is anisotropic."),
          anisotropy,
          paste(sprintf("%d\u00b0 = %.1f", dir_az[dir_ok], dir_ranges[dir_ok]),
                collapse = ", "),
          as.numeric(iso_fit_always))
    } else {
      # The isotropic fit is unusable; the directional sweep is all there is.
      aniso_used      <- TRUE
      effective_range <- dir_max
      vg_used  <- dir_fits[[winner]]$vg
      vgm_used <- attr(dir_fits[[winner]]$fit, "vgm_model")
      if (is.finite(anisotropy) && anisotropy > 1.5)
        .log_warn(
          "estimate_sac_range(): notable anisotropy detected (range ratio %.1f). Directional ranges: %s. Using the maximum.",
          anisotropy,
          paste(sprintf("%d\u00b0 = %.1f", dir_az[dir_ok], dir_ranges[dir_ok]),
                collapse = ", ")
        )
    }
  } else {
    # --- Isotropic variogram (fallback when directional fits fail) ----------
    vg_iso <- try(
      gstat::variogram(..sac_var ~ 1, data = pts, cutoff = cutoff_dist),
      silent = TRUE
    )
    iso_range <- .fit_vgm_range(vg_iso)
    if (is.finite(iso_range)) {
      # A gstat fit that stopped at its iteration limit reports a range that is
      # wherever the optimiser happened to be, not a fitted parameter.  Record
      # it so the rejection block below refuses the VALUE while keeping the
      # variogram for inspection.
      fit_converged <- !identical(attr(iso_range, "converged"), FALSE)
      # as.numeric() strips the vgm_model attribute .fit_vgm_range() attaches;
      # it is re-attached under its documented name below, and leaving both
      # would ship the same object under two attribute names.
      effective_range <- as.numeric(iso_range)
      vg_used  <- vg_iso
      vgm_used <- attr(iso_range, "vgm_model")
    } else {
      # Neither directional nor isotropic succeeded
      .log_warn("estimate_sac_range(): variogram model fit failed; returning NA.")
      return(NA_real_)
    }
  }

  if (!is.finite(effective_range) || effective_range <= 0) {
    .log_warn("estimate_sac_range(): estimated range is non-positive or non-finite; returning NA.")
    return(NA_real_)
  }

  # --- Reject ranges the data cannot actually support -----------------------
  # gstat::fit.variogram() returns a finite number even when the empirical
  # variogram never reaches a sill.  The range is then unidentified and the
  # value is a fitting artefact, not an autocorrelation range.
  #
  # The test is against the CUTOFF distance, not the extent: the empirical
  # variogram is only computed out to `cutoff * max_dist`, so a fitted range
  # beyond that is extrapolation past the longest lag ever observed.  No amount
  # of data at shorter lags identifies it.  (With the default cutoff = 0.5 this
  # is half the extent, but the rationale is the fitted lag limit rather than an
  # arbitrary fraction of the study area.)
  #
  # Returning such a value silently is worse than returning NA, because
  # make_folds(auto_range = TRUE) sizes spatial blocks from it: a range
  # spanning the data yields one block covering everything, which silently
  # defeats blocked cross-validation.
  over_cutoff   <- is.finite(max_supported) && effective_range > max_supported
  # Non-convergence is refused on the same terms and for the same reason: the
  # number is not a fitted parameter.  gstat signals it with a warning and
  # returns anyway, which is why it needs its own test rather than riding on
  # the cutoff bound -- a non-converged range can land inside the bound and
  # would otherwise have sized a block.
  if (over_cutoff || !fit_converged) {
    if (over_cutoff) {
      .log_warn(
        paste0("estimate_sac_range(): fitted range (%.0f) exceeds the largest ",
               "lag the variogram was fitted over (%.4g = %s x cutoff %.0f); the ",
               "empirical variogram never reached a sill, so the range is ",
               "unidentified rather than long. Returning NA. Raise `cutoff` to ",
               "fit longer lags, supply `predictor_vars` to detrend, or set a ",
               "block size explicitly."),
        effective_range, max_supported, format(range_frac), cutoff_dist
      )
    } else {
      .log_warn(
        paste0("estimate_sac_range(): the variogram model did not converge ",
               "(gstat stopped at its iteration limit), so the range it ",
               "reports (%.0f) is where the optimiser halted rather than a ",
               "fitted parameter. Returning NA. Raise `cutoff` to fit longer ",
               "lags, supply `predictor_vars` to detrend, or set a block size ",
               "explicitly. The empirical variogram is attached for ",
               "inspection: plot(type = \"variogram\")."),
        effective_range
      )
    }
    # The VALUE is NA -- the range is genuinely unidentified and must not be
    # used to size blocks -- but the variogram that justified the rejection has
    # already been computed, and throwing it away leaves the user no way to see
    # WHY.  It is exactly the case where plot(type = "variogram") is most worth
    # looking at: a curve that never reaches a sill.  Attributes do not affect
    # is.finite()/is.na(), so every downstream guard behaves as before.
    #
    # The class is set here too, exactly as on the success return.  Without it
    # print.sac_range() never fired for a rejected range, so printing one
    # dumped the whole empirical variogram and the fitted gstat model to the
    # console as raw attributes instead of showing "NA".
    return(structure(
      NA_real_,
      class           = c("sac_range", "numeric"),
      max_dist        = as.numeric(max_dist),
      cutoff_dist     = as.numeric(cutoff_dist),
      crs             = sf::st_crs(pts),
      variogram       = vg_used,
      variogram_model = vgm_used,
      rejected_range  = as.numeric(effective_range),
      rejected_reason = if (over_cutoff)
        "fitted range exceeds the largest lag fitted"
      else "variogram model did not converge"
    ))
  }

  structure(
    effective_range,
    class           = c("sac_range", "numeric"),
    directional     = stats::setNames(usable, as.character(dir_az)),
    anisotropy      = anisotropy,
    # Whether the returned number is the directional maximum (anisotropy
    # established) or the all-pairs estimate.  Without it the `directional`
    # attribute alone cannot tell a caller which of the two was used, and
    # plot(type = "variogram") needs to know whose variogram it is drawing.
    anisotropy_used = isTRUE(aniso_used),
    max_dist        = as.numeric(max_dist),
    cutoff_dist     = as.numeric(cutoff_dist),
    # The CRS the variogram was fitted in.  Its range is a length in these
    # units; summarize_by_cell(deff = "variogram") transforms the points to it
    # before evaluating the correlation function at within-cell distances.
    crs             = sf::st_crs(pts),
    variogram       = vg_used,
    variogram_model = vgm_used
  )
}


#' Print a spatial autocorrelation range
#'
#' Prints the effective range as a plain number, with the directional fit
#' summarised beneath it when one is available.
#'
#' @param x An object of class \code{sac_range}.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.sac_range <- function(x, ...) {
  cat(format(as.numeric(x)), "\n")
  d <- attr(x, "directional")
  a <- attr(x, "anisotropy")
  # any(), not all(): a direction whose variogram never reached a sill is
  # recorded as NA and excluded from the maximum, and suppressing the whole
  # line in that case hides exactly the diagnostic worth seeing.
  if (!is.null(d) && length(d) > 0L && any(is.finite(d))) {
    # Names, not fixed positions: the azimuth sweep is c(0, 45, 90, 135) and
    # an older stored object may carry only c(0, 90).
    labs <- names(d)
    if (is.null(labs)) labs <- as.character(seq_along(d) - 1L)
    cat("  directional: ",
        paste(sprintf("%s deg = %s", labs,
                      ifelse(is.finite(d), format(unname(d)), "unidentified")),
              collapse = ", "),
        sep = "")
    if (is.finite(a)) cat(sprintf("  (ratio %.2f)", a))
    cat("\n")
  }
  invisible(x)
}


#' Compute grid dimensions that respect a minimum block size
#'
#' Given a bounding box and a minimum block edge length, returns nx/ny
#' values such that each cell is at least \code{block_size} wide and tall.
#'
#' @param bb An sf bbox.
#' @param block_size Positive numeric minimum block edge length (CRS units).
#' @return Named list with \code{nx} and \code{ny}.
#' @keywords internal
#' @noRd
# Ceiling on nx * ny for blocked CV.  Matches create_grid_polygons()'s
# max_cells default in spirit: large enough that no legitimate blocking scheme
# reaches it (a 1000 x 1000 grid is 1e6 blocks for a k-fold CV), small enough
# that a unit mistake is refused in milliseconds rather than after minutes of
# allocation.
.block_max_cells <- 1e6

.block_dims_from_size <- function(bb, block_size) {
  w <- as.numeric(bb["xmax"] - bb["xmin"])
  h <- as.numeric(bb["ymax"] - bb["ymin"])
  nx <- max(1L, floor(w / block_size))
  ny <- max(1L, floor(h / block_size))
  list(nx = nx, ny = ny)
}


#' Label the CRS the folds were actually built in
#'
#' \code{make_folds()} projects geographic input with \code{ensure_projected()},
#' so \code{block_size}, \code{sac_range} and \code{buffer} are lengths in a CRS
#' the caller may never have chosen.  Recording a short label in
#' \code{folds$params$crs} makes those units recoverable.  Prefers the
#' \code{AUTHORITY:CODE} form when there is one, falls back to the CRS's own
#' input string, and returns \code{NA_character_} for a missing CRS.
#'
#' @param x An sf/sfc object.
#' @return Character(1).
#' @keywords internal
#' @noRd
.fold_crs_label <- function(x) {
  cr <- tryCatch(sf::st_crs(x), error = function(e) NULL)
  if (is.null(cr) || is.na(cr)) return(NA_character_)
  epsg <- tryCatch(cr$epsg, error = function(e) NULL)
  if (!is.null(epsg) && length(epsg) == 1L && !is.na(epsg))
    return(paste0("EPSG:", epsg))
  inp <- tryCatch(cr$input, error = function(e) NULL)
  if (!is.null(inp) && length(inp) == 1L && !is.na(inp) && nzchar(inp))
    return(as.character(inp))
  wkt <- tryCatch(cr$wkt, error = function(e) NULL)
  if (!is.null(wkt) && length(wkt) == 1L && !is.na(wkt) && nzchar(wkt))
    return(as.character(wkt))
  NA_character_
}


# -----------------------------------------------------------------------------
# Fold Construction
# -----------------------------------------------------------------------------

#' Create spatial cross-validation folds
#'
#' Builds train/test splits using random K-fold, spatial block K-fold, or
#' buffered leave-one-out strategies.
#'
#' For \code{block_kfold}, the default grid sizing is purely geometric and
#' unrelated to the autocorrelation range of the data.  When blocks are
#' smaller than the autocorrelation range, spatially correlated observations
#' leak across folds and CV metrics become optimistic.  Use
#' \code{block_size} to set a minimum block edge length (in CRS units), or
#' set \code{auto_range = TRUE} to estimate the range from an empirical
#' variogram and enforce it automatically.
#'
#' @param points_sf An sf object.  Any Z or M dimension is dropped before
#'   folding: \code{sf::st_distance()} uses every coordinate dimension, so an
#'   XYZ layer would otherwise have elevation folded into every buffer, block
#'   and neighbour distance.  CRS-less points are aligned to a \code{boundary}
#'   or \code{prediction_points} that carries a CRS — reprojected when the
#'   coordinates look like lon/lat, otherwise stamped without reprojection,
#'   warning either way.
#' @param k Integer; number of folds.  Must be a single whole number >= 1 —
#'   a fraction, \code{NA} or a vector is an error, because a non-integer used
#'   to truncate silently and leave the last rows in no test set at all.
#'   Required for \code{"random_kfold"}, \code{"block_kfold"} and
#'   \code{"leave_location_out"}; \code{"buffered_loo"} and \code{"nndm"} are
#'   leave-one-out schemes and ignore it.  Not every method honours it: \code{"buffered_loo"} and \code{"nndm"} are
#'   leave-one-out schemes and always return \code{k = n} regardless of what
#'   was asked for, while \code{"block_kfold"} lowers it when the grid yields
#'   fewer than \code{k} non-empty blocks and \code{"leave_location_out"}
#'   lowers it when there are fewer than \code{k} distinct groups.  Read the
#'   \code{k} element of the returned list rather than assuming the requested
#'   value; a reduction is logged.
#' @param method One of \code{"random_kfold"}, \code{"block_kfold"},
#'   \code{"buffered_loo"}, \code{"leave_location_out"} or \code{"nndm"}.  See
#'   \strong{Details} for what each one does and when it is appropriate.
#' @param seed Optional integer RNG seed.
#' @param block_nx,block_ny Optional grid dimensions for block_kfold.
#'   Ignored when \code{block_size} or \code{auto_range} override them.
#' @param block_multiplier Numeric; target blocks multiplier. Default 3.
#' @param block_size Optional positive numeric minimum block edge length,
#'   \strong{in the units of the CRS the folds are built in}.  When supplied,
#'   grid dimensions are clamped so that every block is at least this wide and
#'   tall.  Takes precedence over \code{block_nx}/\code{block_ny} and
#'   \code{block_multiplier}.
#'
#'   Which CRS that is depends on the input.  Projected input is used as it
#'   stands, so \code{block_size} is in your own CRS's units.  Geographic
#'   (lon/lat) input is projected first by \code{\link{ensure_projected}()},
#'   which picks a local UTM zone or, at wide extents, an equal-area
#'   projection — a CRS you did not choose, whose units are metres but whose
#'   identity varies with the data.  \code{block_size} is then interpreted in
#'   \emph{that} CRS.  The CRS actually used is recorded in
#'   \code{params$crs} of the returned list; project the data yourself before
#'   calling if you want to fix the units in advance.
#'
#'   A \code{block_size} in the wrong unit asks for an enormous grid, so a
#'   request above 1,000,000 blocks is refused with an error naming the grid
#'   dimensions, the extent and the CRS's units, rather than being built.
#' @param auto_range Logical.  If \code{TRUE}, the spatial autocorrelation
#'   range is estimated via \code{estimate_sac_range()} — which fits
#'   directional variograms to account for anisotropy — and used as the
#'   minimum \code{block_size}.  Requires \code{response_var}.  An explicit
#'   \code{block_size} takes precedence.  Default \code{FALSE}.
#' @param range_frac Passed through to \code{estimate_sac_range()} when
#'   \code{auto_range = TRUE}.  A fitted range beyond the longest lag the
#'   empirical variogram was fitted over is rejected as unidentified, and block
#'   sizing falls back to geometry rather than collapsing to a single block.
#'   Default 1.0.
#' @param response_var Character(1) response column name.  Required when
#'   \code{auto_range = TRUE}.
#' @param predictor_vars Optional character vector of predictor column names.
#'   Passed to \code{estimate_sac_range()} for residual variogram estimation.
#' @param boundary Optional polygonal sf/sfc for block_kfold.
#' @param buffer Positive numeric distance for buffered_loo.
#' @param group_var Character(1) naming a column of \code{points_sf} that
#'   identifies the location each observation belongs to.  Required for
#'   \code{method = "leave_location_out"}, which keeps every observation from a
#'   location together in the same fold.  Repeated measurements at the same
#'   site otherwise get split across folds, and the model is scored partly on
#'   sites it has already seen -- which random k-fold reports as excellent
#'   performance.
#' @param prediction_points Optional \code{sf} layer of the locations you
#'   actually intend to predict onto.  Required for \code{method = "nndm"}.
#'   The grid from \code{\link{predict_surface}()} is the natural choice; a
#'   non-POINT layer (grid cells, polygons) is reduced to representative points
#'   first, so the target distances are point-to-point rather than
#'   point-to-polygon — the latter is zero for every cell that contains a
#'   training point, which pulls the target distribution towards zero and
#'   degenerates the CV towards plain leave-one-out.
#' @param min_train For \code{method = "nndm"}: the smallest fraction of the
#'   data any fold's training set may be reduced to by neighbour exclusion.
#'   Default \code{0.5}, as in \code{CAST::nndm()}.
#' @param phi For \code{method = "nndm"}: the largest nearest-neighbour
#'   distance the exclusion will push a held-out point to, in the CRS the
#'   folds are built in.  Default \code{NULL} = the largest
#'   prediction-to-training distance.
#' @details
#' \strong{Fold methods.}
#' \code{"random_kfold"} ignores geography entirely and will overstate
#' performance on autocorrelated data.  \code{"block_kfold"} separates folds
#' geographically.  \code{"buffered_loo"} holds out one point at a time and
#' excludes everything within a fixed \code{buffer}.
#'
#' \code{"leave_location_out"} groups by \code{group_var}, so all
#' observations from a location share a fold.
#'
#' \code{"nndm"} implements the distance-matching principle of Milà et al.
#' (2022): rather than choosing a buffer arbitrarily, it sizes the exclusion
#' around each held-out point so that the resulting training-to-test distance
#' distribution approaches the distribution of distances from your actual
#' prediction locations to the training data.
#'
#' The procedure is the paper's own (as in \code{CAST::nndm()}), and it is
#' deterministic.  Let \eqn{G_{ij}} be the empirical distribution of
#' prediction-to-nearest-training distances and \eqn{G_j^*} the distribution
#' of each held-out point's nearest remaining training point.  Starting from
#' plain leave-one-out, the point with the smallest \eqn{G_j^*} at which the
#' realised distribution exceeds the target -- \eqn{G_j^*(r) > G_{ij}(r)} --
#' has its nearest training neighbour removed, and this repeats until no such
#' point remains, subject to two limits: a point's nearest-neighbour distance
#' is never pushed beyond \code{phi} (default: the largest prediction distance,
#' since a training point already further than every prediction distance has
#' nothing to match), and no fold's training set is stripped below
#' \code{min_train} of the data.
#'
#' The realised distribution is then never \emph{more optimistic} than the
#' target: \eqn{G_j^*(r) \le G_{ij}(r)} up to the granularity of the
#' neighbour distances, which is the property the method exists to deliver.
#' An earlier version of this package drew one random radius per point from
#' \eqn{G_{ij}} and excluded up to the order statistic \emph{closest} to it,
#' which rounds down half the time: on a two-cluster layout the realised
#' distribution exceeded the target by up to 0.17 (13\% of folds had a
#' nearest training point within 50 m against a target of 9\%), an
#' \emph{optimistic} cross-validation.  \code{params$max_ecdf_excess} reports
#' the largest remaining excess; compare \code{params$target_median} with
#' \code{params$realised_median} as well.
#' @references
#' Mila, C., Mateu, J., Pebesma, E. and Meyer, H. (2022). Nearest neighbour
#' distance matching Leave-One-Out Cross-Validation for map validation.
#' \emph{Methods in Ecology and Evolution} \strong{13}, 1304-1316.
#' \doi{10.1111/2041-210X.13851}
#' @param drop_empty_blocks Logical. Default TRUE.
#' @return A list with method, k, folds, assignment, params.  The
#'   \code{train}/\code{test} elements of each fold contain \code{..row_id}
#'   values (equal to row positions when the input has no pre-existing
#'   \code{..row_id} column), consistent with the \code{assignment} tibble.
#'   The returned \code{k} is the number of folds actually built, which is not
#'   always the \code{k} that was requested (see the \code{k} argument above),
#'   and \code{length(folds)} always matches it.
#'
#'   For the methods that work in projected space — \code{"block_kfold"},
#'   \code{"buffered_loo"} and \code{"nndm"} — \code{params} carries a
#'   \code{crs} element naming the CRS the folds were built in (an
#'   \code{"EPSG:code"} string where there is one, otherwise the CRS's input
#'   definition).  Every length in \code{params} — \code{block_size},
#'   \code{sac_range}, \code{buffer}, \code{median_buffer} — is in that CRS's
#'   units, which for geographic input is a CRS
#'   \code{\link{ensure_projected}()} chose rather than one you passed.
#'
#'   Rows whose geometry is empty or has non-finite coordinates are dropped
#'   before folding, with a logged warning naming the count; they appear in no
#'   fold and in no \code{assignment} row.
#' @family cross-validation
#' @examples
#' library(sf)
#' set.seed(1)
#' pts <- st_as_sf(
#'   data.frame(x = runif(30, 0, 1000), y = runif(30, 0, 1000)),
#'   coords = c("x", "y"), crs = 32632
#' )
#' folds <- make_folds(pts, k = 3, method = "block_kfold", seed = 42)
#' folds$assignment          # fold membership per row
#' lengths(folds$folds[[1]]) # train/test row-ID splits
#'
#' # Buffered leave-one-out: neighbours within 100 units excluded from training
#' loo <- make_folds(pts, k = 1, method = "buffered_loo", buffer = 100)
#' @export
make_folds <- function(points_sf, k,
                       method = c("random_kfold", "block_kfold", "buffered_loo",
                                  "leave_location_out", "nndm"),
                       seed = NULL, block_nx = NULL, block_ny = NULL,
                       block_multiplier = 3, block_size = NULL,
                       auto_range = FALSE, range_frac = 1.0, response_var = NULL,
                       group_var = NULL, prediction_points = NULL,
                       predictor_vars = NULL, boundary = NULL,
                       buffer = NULL, min_train = 0.5, phi = NULL,
                       drop_empty_blocks = TRUE) {
  method <- match.arg(method)

  cleanup <- .with_seed(seed)
  on.exit(cleanup(), add = TRUE)

  if (!inherits(points_sf, "sf")) stop("make_folds(): `points_sf` must be an sf object.")
  # Guard zero rows here rather than letting each method fail in its own way:
  # block_kfold reaches sf::st_bbox(), which returns all-NA for an empty layer
  # and makes st_as_sfc() abort with an opaque "!anyNA(x) is not TRUE".
  if (nrow(points_sf) == 0L)
    stop("make_folds(): `points_sf` has no rows; there is nothing to split into folds.",
         call. = FALSE)
  # MULTIPOINT must be coerced too, not merely admitted: st_coordinates()
  # returns one row per VERTEX, so any multi-vertex feature makes xy[i, ] a
  # different feature than row i, and every fold below misaligns silently.
  # coerce_to_points() takes centroids, matching prep_model_data().
  if (!("..row_id" %in% names(points_sf)))
    points_sf$..row_id <- seq_len(nrow(points_sf))
  # A duplicated ID cannot be split on: match() resolves it to its first row,
  # so the other rows are never scored and can sit in train and test at once.
  if (anyDuplicated(points_sf$..row_id))
    stop(sprintf(paste0("make_folds(): `..row_id` has %d duplicated value(s). ",
                        "Row IDs are what fold splits are made of, so every row ",
                        "needs its own; drop the column to let make_folds() ",
                        "number the rows, or make the IDs unique."),
                 sum(duplicated(points_sf$..row_id))), call. = FALSE)
  if (anyNA(points_sf$..row_id))
    stop("make_folds(): `..row_id` contains NA; every row needs an ID.",
         call. = FALSE)
  # k was never validated.  A non-integer k truncated silently inside
  # rep(floor(n/k), k): only floor(k) folds were built, the size remainder was
  # written past the end of the vector, and the last rows of the permutation
  # landed in NO test set (assignment$fold == 0) with no condition raised --
  # while `k` was echoed back unchanged, so length(folds) != k, contradicting
  # @return.  cv_spatial() then scored 18 of 20 rows and reported
  # attempted == succeeded.
  if (!missing(k) && !is.null(k)) {
    if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k != round(k) ||
        k < 1)
      stop("make_folds(): `k` must be a single whole number >= 1; got ",
           paste(format(k), collapse = ", "), ".", call. = FALSE)
    k <- as.integer(k)
  }
  # The provenance probe is taken on the geometry AS SUPPLIED -- before
  # pointization -- because that is what the cv_*() wrappers will probe too.
  row_probe <- .fold_row_probe(points_sf)
  if (!all(sf::st_geometry_type(points_sf, by_geometry = TRUE) == "POINT"))
    points_sf <- coerce_to_points(points_sf, "auto")

  # Drop rows with no usable coordinates, AFTER ..row_id is stamped so the
  # survivors keep their original row identities and the dropped rows simply
  # never appear in a fold.  st_coordinates() yields one all-NA row per EMPTY
  # POINT rather than zero rows, so a row-count check alone lets them through:
  # block_kfold's st_intersects() then returns integer(0) for them, ..block_id
  # goes NA, st_distance() is all-NA and the nearest-block rescue used to die
  # with "replacement has length zero".  Matches voronoi_seeds_kmeans(), which
  # drops the same rows with the same warning.  Applied for every method --
  # random_kfold would happily put an unplottable point in a fold, and nndm
  # and buffered_loo both feed the coordinates to distance code.
  # Drop any Z/M dimension.  gstat::variogram() and sf::st_distance() use ALL
  # coordinate dimensions, so an XYZ layer (a sounding, a soil profile, a lidar
  # return) had its elevation folded into every lag: estimate_sac_range()
  # returned a length in 3-D while every consumer of it -- the block grid,
  # buffered_loo's buffer, nndm's neighbour distances, summarize_by_cell --
  # works in 2-D map distance.  Nothing warned.
  if (any(c("Z", "M") %in% sf::st_dimension(points_sf, NA_if_empty = FALSE)) ||
      !is.null(attr(sf::st_geometry(points_sf), "z_range")) ||
      !is.null(attr(sf::st_geometry(points_sf), "m_range"))) {
    points_sf <- sf::st_zm(points_sf, drop = TRUE, what = "ZM")
  }

  bad_geom <- sf::st_is_empty(points_sf)
  if (!all(bad_geom)) {
    xy_chk <- suppressWarnings(sf::st_coordinates(points_sf[!bad_geom, ]))
    if (nrow(xy_chk) == sum(!bad_geom) && ncol(xy_chk) >= 2L)
      bad_geom[!bad_geom] <- !stats::complete.cases(xy_chk[, 1:2, drop = FALSE])
  }
  if (any(bad_geom)) {
    .log_warn("make_folds(): dropping %d point(s) with empty or non-finite coordinates.",
              sum(bad_geom))
    points_sf <- points_sf[!bad_geom, , drop = FALSE]
    if (nrow(points_sf) == 0L)
      stop("make_folds(): `points_sf` has no usable coordinates; there is nothing to split into folds.",
           call. = FALSE)
  }

  # One probe for every return path: .ret() is the single exit, so the folds
  # object cannot ship without the fingerprint that lets a later cv_*() refuse
  # it if it is handed the wrong data.  (Taken above, on the input geometry.)
  .ret <- function(method, k, folds, assignment, params)
    list(method = method, k = k, folds = folds, assignment = assignment,
         params = c(params, list(row_probe = row_probe)))

  # ---- RANDOM K-FOLD ----
  if (method == "random_kfold") {
    n <- nrow(points_sf)
    if (k < 2) k <- 2L
    if (k > n) { .log_warn("make_folds(random_kfold): k > n; reducing."); k <- n }
    row_ids <- points_sf$..row_id
    idx <- sample.int(n, n)
    sizes <- rep(floor(n / k), k)
    remainder <- n - sum(sizes)
    if (remainder > 0) sizes[seq_len(remainder)] <- sizes[seq_len(remainder)] + 1L
    splits <- vector("list", k); start <- 1L; assign_vec <- integer(n)
    for (j in seq_len(k)) {
      stop_i <- start + sizes[j] - 1L
      test_idx <- idx[start:stop_i]; train_idx <- setdiff(idx, test_idx)
      # Express splits in ..row_id values (not positions) so that
      # .remap_folds() / .cv_fit_one_fold() interpret them correctly even
      # when the input rows carry non-sequential IDs (e.g. after
      # prep_model_data() has dropped rows inside cv_gwr()/cv_bayes()).
      splits[[j]] <- list(train = row_ids[train_idx], test = row_ids[test_idx])
      assign_vec[test_idx] <- j; start <- stop_i + 1L
    }
    return(.ret(method, k, splits,
                .safe_tibble(row_id = points_sf$..row_id, fold = assign_vec),
                list(seed = seed)))
  }

  # ---- BLOCK K-FOLD ----
  if (method == "block_kfold") {
    if (k < 2) k <- 2L
    pts <- ensure_projected(points_sf)
    # A CRS-less `points_sf` leaves .crs_or_null(pts) NULL, so the boundary
    # kept its own CRS and sf aborted the intersection with
    # "st_crs(x) == st_crs(y) is not TRUE" -- while prep_model_data(), and
    # therefore every cv_*() wrapper, handles the identical combination by
    # aligning the points to the boundary.  Do the same here.
    if (is.na(sf::st_crs(pts)) && !is.null(boundary) &&
        !is.na(sf::st_crs(boundary))) {
      pts <- .transform_or_stamp(pts, sf::st_crs(boundary),
                                 what = "points_sf", caller = "make_folds")
    }
    reg <- if (!is.null(boundary)) {
      b <- ensure_projected(boundary, .crs_or_null(pts))
      if (inherits(b, "sfc")) b <- sf::st_sf(geometry = b)
      b <- .safe_make_valid(sf::st_union(b))
      mat <- sf::st_intersects(pts, b, sparse = FALSE)
      inside_any <- apply(mat, 1L, any)
      if (!all(inside_any)) {
        bb_pts <- sf::st_as_sfc(sf::st_bbox(pts)) |> sf::st_set_crs(sf::st_crs(pts))
        b <- suppressWarnings(
          sf::st_union(.safe_make_valid(sf::st_sf(geometry = c(b, bb_pts)))))
      }
      b
    } else {
      sf::st_as_sfc(sf::st_bbox(pts)) |> sf::st_set_crs(sf::st_crs(pts)) |> sf::st_sf()
    }

    # --- Autocorrelation-aware block sizing ---
    sac_range <- NA_real_
    if (isTRUE(auto_range) && !is.null(response_var)) {
      # `seed` here is make_folds()'s own, whose default is NULL -- meaning
      # "do not seed the FOLD assignment".  Forwarding that NULL re-opened the
      # unseeded n_max subsample inside estimate_sac_range(): the range, and
      # with it the block size and every fold, changed on each default call
      # (measured 199.45 / 204.18 / 208.99 across three identical calls).  The
      # subsample is an internal approximation, not part of the answer, so it
      # keeps estimate_sac_range()'s own reproducible default when make_folds()
      # was not given a seed.
      sac_range <- estimate_sac_range(pts, response_var = response_var,
                                      predictor_vars = predictor_vars,
                                      range_frac = range_frac,
                                      seed = if (is.null(seed)) 123L else seed)
      if (is.finite(sac_range) && sac_range > 0) {
        message(sprintf(
          "make_folds(block_kfold): estimated spatial autocorrelation range = %.1f CRS units; using as minimum block size.",
          sac_range
        ))
        # auto_range sets block_size only if the caller didn't supply one
        if (is.null(block_size)) {
          block_size <- sac_range
        } else if (block_size < sac_range) {
          .log_warn(
            "make_folds(block_kfold): supplied block_size (%.1f) is smaller than the estimated autocorrelation range (%.1f). Spatial CV may still leak correlated information.",
            block_size, sac_range
          )
          warning(
            sprintf("make_folds(): block_size (%.1f) < estimated autocorrelation range (%.1f). Consider increasing block_size to reduce information leakage across folds.",
                    block_size, sac_range),
            call. = FALSE
          )
        }
      } else {
        .log_warn("make_folds(block_kfold): auto_range requested but estimation returned NA; falling back to geometric blocks.")
      }
    } else if (isTRUE(auto_range) && is.null(response_var)) {
      .log_warn("make_folds(block_kfold): auto_range = TRUE but response_var is NULL; cannot estimate range. Falling back to geometric blocks.")
      warning("make_folds(): auto_range requires response_var; ignoring.", call. = FALSE)
    }

    bb <- sf::st_bbox(reg)

    # Determine grid dimensions: block_size constrains nx/ny
    if (!is.null(block_size) && is.numeric(block_size) && block_size > 0) {
      size_dims <- .block_dims_from_size(bb, block_size)
      # If the caller also supplied explicit block_nx/block_ny, warn about override
      if (!is.null(block_nx) || !is.null(block_ny)) {
        .log_warn(
          "make_folds(block_kfold): block_size (%.1f) overrides explicit block_nx/block_ny.",
          block_size
        )
      }
      nx <- size_dims$nx
      ny <- size_dims$ny

      # Ensure at least k blocks so each fold can get one
      if (nx * ny < k) {
        .log_warn(
          "make_folds(block_kfold): block_size produces only %d blocks (< k = %d). Reducing k to match.",
          nx * ny, k
        )
        k <- max(2L, nx * ny)
      }
    } else if (is.null(block_nx) || is.null(block_ny)) {
      w  <- as.numeric(bb["xmax"] - bb["xmin"])
      h  <- as.numeric(bb["ymax"] - bb["ymin"])
      ratio <- if (h > 0) w / h else 1
      target_blocks <- max(1L, round(block_multiplier * k))
      nx <- max(1L, round(sqrt(target_blocks * ratio)))
      ny <- max(1L, round(max(1, target_blocks / nx)))

      # Diagnostic: warn if resulting block size is small relative to SAC range
      if (is.finite(sac_range) && sac_range > 0) {
        cell_w <- w / nx
        cell_h <- h / ny
        min_cell <- min(cell_w, cell_h)
        if (min_cell < sac_range) {
          .log_warn(
            "make_folds(block_kfold): geometric block size (%.1f) is smaller than the estimated autocorrelation range (%.1f). Consider setting block_size >= %.0f or auto_range = TRUE to avoid information leakage.",
            min_cell, sac_range, sac_range
          )
          warning(
            sprintf("make_folds(): block dimension (%.1f) < autocorrelation range (%.1f). Spatial CV may leak correlated information across folds. Pass block_size = %.0f or auto_range = TRUE.",
                    min_cell, sac_range, ceiling(sac_range)),
            call. = FALSE
          )
        }
      }
    } else {
      nx <- as.integer(block_nx); ny <- as.integer(block_ny)

      # Diagnostic: warn if user-supplied nx/ny yield blocks smaller than SAC
      if (is.finite(sac_range) && sac_range > 0) {
        w  <- as.numeric(bb["xmax"] - bb["xmin"])
        h  <- as.numeric(bb["ymax"] - bb["ymin"])
        cell_w <- w / nx; cell_h <- h / ny
        min_cell <- min(cell_w, cell_h)
        if (min_cell < sac_range) {
          .log_warn(
            "make_folds(block_kfold): user-supplied grid (%dx%d) yields blocks of ~%.1f units, smaller than estimated autocorrelation range (%.1f).",
            nx, ny, min_cell, sac_range
          )
          warning(
            sprintf("make_folds(): block_nx/block_ny yield blocks smaller than autocorrelation range (%.1f). Consider using block_size = %.0f.",
                    sac_range, ceiling(sac_range)),
            call. = FALSE
          )
        }
      }
    }

    # st_make_grid() builds every cell before anything downstream can look at
    # the count, so an unnoticed unit mistake -- block_size in kilometres on a
    # metre CRS, or a CRS-less unit square taken for lon/lat, where 0.25 means
    # a quarter-metre block on a 109 km extent -- asks for 1e8 to 1e11 cells
    # and exhausts memory instead of being refused.  create_grid_polygons()
    # guards exactly this mistake with max_cells and names the CRS units;
    # blocked CV is the sibling that did not.
    n_cells_est <- as.numeric(nx) * as.numeric(ny)
    if (is.finite(n_cells_est) && n_cells_est > .block_max_cells) {
      unit_lbl <- tryCatch({
        u <- sf::st_crs(reg)$units_gdal
        if (is.null(u) || is.na(u) || !nzchar(u)) "CRS units" else u
      }, error = function(e) "CRS units")
      stop(sprintf(paste0("make_folds(block_kfold): the requested grid is %d x ",
                          "%d = %s cells, above the %s this function will ",
                          "build. Check that `block_size` (%s) is expressed in ",
                          "the data's CRS units (%s) over an extent of %s x %s; ",
                          "a value in the wrong unit is the usual cause."),
                  nx, ny, format(n_cells_est, big.mark = ",", scientific = FALSE),
                  format(.block_max_cells, big.mark = ",", scientific = FALSE),
                  if (is.null(block_size)) "unset" else format(block_size),
                  unit_lbl,
                  format(signif(as.numeric(bb["xmax"] - bb["xmin"]), 4)),
                  format(signif(as.numeric(bb["ymax"] - bb["ymin"]), 4))),
           call. = FALSE)
    }

    grid <- sf::st_make_grid(reg, n = c(nx, ny), what = "polygons", square = TRUE)
    grid <- .safe_make_valid(grid)
    reg_union <- .safe_make_valid(sf::st_union(reg))
    grid <- suppressWarnings(sf::st_intersection(grid, reg_union))
    grid_sf <- sf::st_as_sf(grid)
    hits <- sf::st_intersects(pts, grid_sf)
    block_id <- vapply(hits, function(ix) if (length(ix)) ix[1] else NA_integer_, 1L)
    pts$..block_id <- block_id

    if (drop_empty_blocks) {
      used_blocks <- sort(unique(pts$..block_id[!is.na(pts$..block_id)]))
      grid_sf <- grid_sf[used_blocks, , drop = FALSE]
      pts$..block_id <- match(pts$..block_id, used_blocks)
    }
    if (anyNA(pts$..block_id)) {
      cent <- suppressWarnings(sf::st_centroid(sf::st_geometry(grid_sf)))
      na_idx <- which(is.na(pts$..block_id))
      dmat <- as.matrix(sf::st_distance(sf::st_geometry(pts[na_idx, ]), cent))
      # vapply(), not apply(): a point whose distances are all NA (an empty
      # geometry that survived upstream, or a grid with no finite centroid)
      # makes which.min() return integer(0), and apply() then simplifies the
      # whole result to a list() -- which assigns back as "replacement has
      # length zero".  Keep the NA instead and let the B/k guards below react.
      pts$..block_id[na_idx] <- vapply(
        seq_len(nrow(dmat)),
        function(i) {
          w <- which.min(dmat[i, ])
          if (length(w)) as.integer(w) else NA_integer_
        },
        integer(1)
      )
    }
    B <- max(pts$..block_id, na.rm = TRUE)
    # One block means one fold whose training set is empty -- blocked CV
    # silently degenerating into nothing at all.  It happens whenever the block
    # size exceeds half the extent, which an accepted autocorrelation range can
    # do (estimate_sac_range() only rejects ranges above half the DIAGONAL, a
    # factor of sqrt(2) looser than what .block_dims_from_size() needs).
    if (B < 2L) {
      # block_size is NULL when the caller set block_nx/block_ny directly, and
      # sprintf() on a zero-length argument yields character(0) -- an empty
      # error message.
      how <- if (is.null(block_size))
        "the requested block_nx/block_ny" else
        sprintf("the block size (%s)", format(block_size))
      stop(sprintf(paste0("make_folds(block_kfold): %s produces a single block ",
                          "covering the whole extent, so there is no spatial ",
                          "split to make and the one fold would have an empty ",
                          "training set. Pass a smaller `block_size` (or more ",
                          "blocks), or set auto_range = FALSE if the size came ",
                          "from the estimated autocorrelation range."), how),
           call. = FALSE)
    }
    if (B < k) { .log_warn("make_folds(block_kfold): blocks < k; reducing k."); k <- B }

    blk_sizes <- as.integer(table(factor(pts$..block_id, levels = seq_len(B))))
    order_blk <- order(blk_sizes, decreasing = TRUE)
    fold_loads <- integer(k); fold_blocks <- vector("list", k)
    for (i in seq_along(order_blk)) {
      j <- which(fold_loads == min(fold_loads))
      if (length(j) > 1) j <- sample(j, 1L)
      fold_blocks[[j]] <- c(fold_blocks[[j]], seq_len(B)[order_blk[i]])
      fold_loads[j] <- fold_loads[j] + blk_sizes[order_blk[i]]
    }

    if (min(fold_loads) > 0L && max(fold_loads) > 3 * min(fold_loads)) {
      .log_warn("make_folds(block_kfold): fold size imbalance \u2014 largest fold has %d obs vs %d in smallest.",
                max(fold_loads), min(fold_loads))
    }

    splits <- vector("list", k); assign_vec <- integer(nrow(pts))
    row_ids <- pts$..row_id
    for (j in seq_len(k)) {
      test_idx <- which(pts$..block_id %in% fold_blocks[[j]])
      # Row-ID (not positional) splits — see comment in the random_kfold branch.
      splits[[j]] <- list(train = row_ids[setdiff(seq_len(nrow(pts)), test_idx)],
                          test  = row_ids[test_idx])
      assign_vec[test_idx] <- j
    }
    return(.ret(method, k, splits,
                .safe_tibble(row_id = pts$..row_id, fold = assign_vec),
                list(seed = seed, grid_nx = nx, grid_ny = ny, blocks_used = B,
                     block_multiplier = block_multiplier,
                     block_size = block_size,
                     sac_range = sac_range,
                     auto_range = auto_range,
                     # `block_size` and `sac_range` are lengths in the CRS the
                     # folds were actually built in, which is NOT necessarily
                     # the CRS the caller passed: geographic input is projected
                     # by ensure_projected() to a CRS chosen here.  Record it
                     # so the units are recoverable.
                     crs = .fold_crs_label(pts),
                     boundary_supplied = !is.null(boundary))))
  }

  # ---- BUFFERED LOO ----
  if (method == "buffered_loo") {
    if (is.null(buffer) || !is.numeric(buffer) || buffer <= 0)
      stop("make_folds(buffered_loo): `buffer` (positive numeric) is required.")
    pts <- ensure_projected(points_sf)
    n <- nrow(pts)
    # The cost is QUADRATIC in n whatever the buffer: every one of the n
    # splits stores its own training-row vector of length about n, so the
    # fold object alone is ~4 n^2 bytes.  The old threshold of 20000 admitted
    # 1.6 GB of splits (measured 1.7 GB of R objects at n = 19999), before
    # .remap_folds() copied them.  The cap is now stated in bytes and the
    # message says what the request would cost.
    split_gb <- 4 * as.numeric(n)^2 / 1024^3
    if (split_gb > 0.6) {
      stop(sprintf(paste0(
        "make_folds(buffered_loo): n = %d would produce %d leave-one-out ",
        "splits holding about %.2f GB of row indices (the training vector of ",
        "every split is ~n long, so storage grows as n^2), before the model is ",
        "refitted %d times. Use 'block_kfold' instead, or subset the data."),
        n, n, split_gb, n), call. = FALSE)
    } else if (split_gb > 0.1) {
      .log_warn(
        "make_folds(buffered_loo): n = %d produces %d splits holding about %.2f GB of row indices, and the model will be refitted %d times; this may be slow.",
        n, n, split_gb, n
      )
    }
    nb <- sf::st_is_within_distance(sf::st_geometry(pts), sf::st_geometry(pts),
                                     dist = buffer)
    splits <- vector("list", n)
    row_ids <- pts$..row_id
    for (i in seq_len(n)) {
      excl <- sort(unique(nb[[i]]))
      # Row-ID (not positional) splits — see comment in the random_kfold branch.
      splits[[i]] <- list(train = row_ids[setdiff(seq_len(n), excl)],
                          test  = row_ids[i])
    }

    # A buffer that spans the data leaves EVERY fold with an empty training
    # set.  Those folds used to pass straight through .remap_folds() and get
    # dropped one at a time inside .cv_fit_one_fold(), so the only symptom was
    # a generic "all folds failed" warning at the very end of the CV run --
    # long after the cause could be acted on.  Say it here instead.
    n_train_each <- vapply(splits, function(s) length(s$train), integer(1))
    if (!any(n_train_each >= 2L))
      stop(sprintf(paste0("make_folds(buffered_loo): a buffer of %s excludes ",
                          "so much of the data that no fold retains 2 training ",
                          "points (largest training set: %d of %d). The buffer ",
                          "spans the data; use a smaller one, or 'block_kfold'."),
                   format(buffer), max(n_train_each), n), call. = FALSE)

    return(.ret(method, n, splits,
                .safe_tibble(row_id = pts$..row_id, fold = seq_len(n)),
                list(buffer = buffer, crs = .fold_crs_label(pts))))
  }

  # ---- LEAVE-LOCATION-OUT (grouped) ----
  if (method == "leave_location_out") {
    if (is.null(group_var) || !(group_var %in% names(points_sf)))
      stop("make_folds(leave_location_out): `group_var` must name a column of ",
           "`points_sf` identifying the location each observation belongs to.",
           call. = FALSE)

    grp <- as.character(sf::st_drop_geometry(points_sf)[[group_var]])
    # An empty string is refused on the same terms as NA: it is not a
    # location.  It also used to be worse than NA -- names<- and [ by name
    # treat "" as "no name", so grp_fold[""] returned NA and every row with a
    # blank label silently got fold NA, entered no test set, and put 36 NA
    # row-ids into the train splits, with no condition raised.
    if (anyNA(grp) || any(!nzchar(grp)))
      stop("make_folds(leave_location_out): `group_var` contains NA or empty ",
           "labels; every observation must belong to a named location.",
           call. = FALSE)

    row_ids <- points_sf$..row_id
    ug <- unique(grp)
    n_g <- length(ug)
    if (n_g < 2L)
      stop("make_folds(leave_location_out): need at least 2 distinct groups.",
           call. = FALSE)

    # k groups per fold; k >= n_g degenerates to leave-one-group-out.
    if (k > n_g) {
      .log_warn("make_folds(leave_location_out): k = %d exceeds the %d distinct groups; using leave-one-group-out.", k, n_g)
      k <- n_g
    }
    if (k < 2L) k <- 2L

    grp_fold <- sample(rep(seq_len(k), length.out = n_g))
    # Positional lookup via match(), never by name: name-indexing is what
    # turned an unusual label into an NA fold.
    fold_of_row <- grp_fold[match(grp, ug)]

    splits <- vector("list", k)
    for (j in seq_len(k)) {
      test_rows <- row_ids[fold_of_row == j]
      splits[[j]] <- list(train = row_ids[fold_of_row != j], test = test_rows)
    }

    return(.ret(method, k, splits,
                .safe_tibble(row_id = row_ids, fold = unname(fold_of_row)),
                list(seed = seed, group_var = group_var, n_groups = n_g)))
  }

  # ---- NNDM: nearest-neighbour distance matching ----
  if (method == "nndm") {
    pts <- ensure_projected(points_sf)
    n <- nrow(pts)
    if (n > 5000L)
      stop(sprintf("make_folds(nndm): n = %d exceeds the safety threshold of 5000. Fold construction sorts distances from every point to every other (O(n^2) time), and NNDM then produces n leave-one-out folds, so the model is refitted n times. Use 'block_kfold' instead, or subset your data.", n),
           call. = FALSE)
    if (n < 3L)
      stop("make_folds(nndm): need at least 3 points.", call. = FALSE)

    if (is.null(prediction_points))
      stop("make_folds(nndm): `prediction_points` is required -- NNDM matches ",
           "the training-to-test distance distribution to the distance ",
           "distribution from your actual prediction locations to the training ",
           "data. Use the grid you intend to predict onto (see ",
           "predict_surface()).", call. = FALSE)
    # Prediction locations are frequently the grid CELLS a surface will be
    # drawn on, not their centres.  Every other geometry input in this file is
    # pointized first; this one was not, so a POLYGON layer gave point-to-
    # POLYGON distances -- zero for every cell that contains a training point.
    # The target distance distribution was then pulled towards 0, far fewer
    # neighbours were excluded, and the distance-matched CV quietly degenerated
    # towards plain LOO: the optimistic direction.
    if (inherits(prediction_points, "sfc"))
      prediction_points <- sf::st_sf(geometry = prediction_points)
    if (!all(sf::st_geometry_type(prediction_points,
                                  by_geometry = TRUE) == "POINT")) {
      .log_info(paste0("make_folds(nndm): `prediction_points` are not POINTs; ",
                       "using representative points so the target distances ",
                       "are point-to-point."))
      prediction_points <- coerce_to_points(prediction_points, "auto")
    }
    if (is.na(sf::st_crs(pts)) && !is.na(sf::st_crs(prediction_points))) {
      pts <- .transform_or_stamp(pts, sf::st_crs(prediction_points),
                                 what = "points_sf", caller = "make_folds")
    }
    pred <- ensure_projected(prediction_points, target_crs = .crs_or_null(pts))
    pred <- sf::st_zm(pred, drop = TRUE, what = "ZM")

    # Target: distance from each prediction location to its nearest training
    # point.  This is the distance regime the model will actually face, and the
    # one an arbitrary fixed buffer has no reason to reproduce.
    g_target <- .nn_dist_to(pred, pts)
    g_target <- g_target[is.finite(g_target)]
    if (!length(g_target))
      stop("make_folds(nndm): could not compute prediction-to-training ",
           "distances.", call. = FALSE)

    if (!is.numeric(min_train) || length(min_train) != 1L ||
        !is.finite(min_train) || min_train <= 0 || min_train >= 1)
      stop("make_folds(nndm): `min_train` must be a single number in (0, 1).",
           call. = FALSE)
    if (is.null(phi)) phi <- max(g_target)
    if (!is.numeric(phi) || length(phi) != 1L || !is.finite(phi) || phi < 0)
      stop("make_folds(nndm): `phi` must be a single non-negative number.",
           call. = FALSE)

    # ---- The paper's procedure (Mila et al. 2022; CAST::nndm) ---------------
    # Deterministic: no radii are drawn.  Starting from plain LOO, the held-out
    # point with the SMALLEST current nearest-neighbour distance at which the
    # realised distribution exceeds the target has its nearest training
    # neighbour removed, and this repeats until no such point remains.
    #
    # Implemented as a single sweep over the points in increasing order of
    # their current nearest-neighbour distance.  Removing a neighbour only ever
    # INCREASES that point's distance, which moves it later in the order and
    # lowers the realised distribution at every smaller value -- so a position
    # that has passed the test can never fail it again, and the sweep never
    # has to look back.  The greedy choice ("smallest violator first") is
    # exactly the paper's.
    xy      <- sf::st_coordinates(pts)[, 1:2, drop = FALSE]
    row_ids <- pts$..row_id
    rmin    <- min_train * n
    # A point may lose at most n - 1 - rmin neighbours, so only that many
    # (plus one) of its sorted neighbour distances are ever consulted.  The
    # rule the sweep applies is (n - 1 - removed) > rmin, which permits
    # ceiling(n - 1 - rmin) removals -- floor() left one column short whenever
    # n * min_train is fractional (any odd n at the default 0.5), and the
    # ncol() guard in the loop then stopped the sweep one removal early, so
    # every fold that reached the floor kept one training point more than the
    # paper's procedure.
    k_need <- as.integer(min(n - 1L, ceiling(n - 1L - rmin) + 1L))
    k_need <- max(1L, k_need)
    if (requireNamespace("FNN", quietly = TRUE)) {
      kn <- FNN::get.knn(xy, k = min(k_need + 1L, n - 1L))
      nn_d <- kn$nn.dist; nn_i <- kn$nn.index
      # get.knn() means to exclude the query point, but an exact tie defeats
      # it: with co-located points it returns the query's OWN index in place
      # of one of its duplicates (verified: rbind(c(0,0), c(0,0), ...) gives
      # row 1 the neighbour list 1, 3, 4, ... -- row 2, its twin, is absent).
      # Simply dropping the self entry and padding with Inf left the twin
      # missing from the list altogether, so the sweep could never exclude it:
      # the fold holding out a repeated measurement trained on its
      # exact-location twin (leakage) while params$realised_distances reported
      # a large buffer for a fold that had none.  The displaced entry is
      # always a co-located point at distance 0, so put one back.
      key <- paste(xy[, 1L], xy[, 2L], sep = "\r")
      colocated <- if (anyDuplicated(key)) split(seq_len(n), key) else NULL
      for (i in seq_len(n)) {
        self <- which(nn_i[i, ] == i)
        if (!length(self)) next
        twins <- if (is.null(colocated)) integer(0)
                 else setdiff(colocated[[key[i]]], c(i, nn_i[i, -self]))
        for (s in seq_along(self)) {
          if (s <= length(twins)) {
            # Distance in that slot is already 0, which is correct for a twin.
            nn_i[i, self[s]] <- twins[s]
          } else {
            # No duplicate to restore (a self entry with no twin should not
            # happen, but do not let it become a phantom neighbour).
            nn_d[i, self[s]] <- Inf
            nn_i[i, self[s]] <- NA_integer_
          }
        }
        # Keep each row sorted by distance, as the sweep assumes.
        o <- order(nn_d[i, ])
        nn_d[i, ] <- nn_d[i, o]; nn_i[i, ] <- nn_i[i, o]
      }
      if (ncol(nn_d) < k_need) {
        pad  <- k_need - ncol(nn_d)
        nn_d <- cbind(nn_d, matrix(Inf, n, pad))
        nn_i <- cbind(nn_i, matrix(NA_integer_, n, pad))
      }
      nn_d <- nn_d[, seq_len(k_need), drop = FALSE]
      nn_i <- nn_i[, seq_len(k_need), drop = FALSE]
    } else {
      nn_d <- matrix(Inf, n, k_need); nn_i <- matrix(NA_integer_, n, k_need)
      for (i in seq_len(n)) {
        d <- sqrt((xy[, 1] - xy[i, 1])^2 + (xy[, 2] - xy[i, 2])^2)
        d[i] <- Inf
        o <- order(d)[seq_len(k_need)]
        nn_d[i, ] <- d[o]; nn_i[i, ] <- o
      }
    }

    g_sorted <- sort(g_target)
    n_g      <- length(g_sorted)
    G_target <- function(r) findInterval(r, g_sorted) / n_g   # right-continuous ECDF

    removed <- integer(n)
    Gjstar  <- nn_d[, 1L]
    o  <- order(Gjstar); sv <- Gjstar[o]; si <- o
    k  <- 1L
    n_iter <- 0L
    while (k <= n) {
      r   <- sv[k]
      j   <- si[k]
      cnt <- findInterval(r, sv)                    # realised count <= r
      violates <- is.finite(r) && (cnt / n) > G_target(r) + 1e-12 &&
                  r <= phi && (n - 1L - removed[j]) > rmin &&
                  removed[j] + 1L < ncol(nn_d)
      if (violates) {
        removed[j] <- removed[j] + 1L
        newv <- nn_d[j, removed[j] + 1L]
        sv <- sv[-k]; si <- si[-k]
        # Insert after every smaller value and, among equal values, after
        # those belonging to points with a smaller index -- the same tie
        # order as which.min() in the reference implementation, so that on
        # clustered data (where pushed points pile up at the same
        # cluster-to-cluster distance) the SAME points get pushed.
        pos <- findInterval(newv, sv)
        while (pos > 0L && sv[pos] == newv && si[pos] > j) pos <- pos - 1L
        sv <- append(sv, newv, after = pos)
        si <- append(si, j,    after = pos)
        n_iter <- n_iter + 1L
      } else {
        k <- k + 1L
      }
    }
    Gjstar[si] <- sv

    splits     <- vector("list", n)
    n_excluded <- removed
    realised   <- Gjstar
    for (i in seq_len(n)) {
      drop_idx <- if (removed[i] > 0L) nn_i[i, seq_len(removed[i])] else integer(0)
      drop_idx <- drop_idx[!is.na(drop_idx)]
      train_i  <- setdiff(seq_len(n), c(i, drop_idx))
      splits[[i]] <- list(train = row_ids[train_i], test = row_ids[i])
    }
    # The paper's own diagnostic: largest excess of realised over target.
    fin <- is.finite(realised)
    max_excess <- if (any(fin)) {
      rs <- sort(realised[fin])
      max(findInterval(rs, rs) / length(rs) - G_target(rs))
    } else NA_real_

    # Exclusion is not always needed.  When prediction locations sit no further
    # from the training data than training points sit from each other, plain
    # LOO already reproduces the target distribution and nothing is removed --
    # that is the correct answer, not a failure.  Conversely NNDM cannot pull
    # training points closer, so it cannot match a target that is shorter than
    # the training nearest-neighbour distances.
    .log_info(paste0("make_folds(nndm): %d folds. Target prediction-to-training ",
                     "distance: median %.1f. Realised test-to-training distance: ",
                     "median %.1f, excluding a median of %.1f training point(s) ",
                     "per fold (%d removals in total; largest remaining excess ",
                     "of realised over target ECDF %.3f)."),
              n, stats::median(g_target), stats::median(realised[fin]),
              stats::median(n_excluded), n_iter, max_excess)

    return(.ret(method, n, splits,
                .safe_tibble(row_id = row_ids, fold = seq_len(n)),
                list(seed = seed, n_prediction_points = nrow(pred),
                     crs = .fold_crs_label(pts),
                     min_train = min_train, phi = phi,
                     # The effective buffer each fold ended up with.
                     median_buffer   = stats::median(realised[fin]),
                     median_excluded = stats::median(n_excluded),
                     n_removed_total = n_iter,
                     max_ecdf_excess = max_excess,
                     target_median   = stats::median(g_target),
                     realised_median = stats::median(realised[fin]),
                     target_distances   = g_target,
                     realised_distances = realised)))
  }

  stop("make_folds(): unsupported method.")
}


# -----------------------------------------------------------------------------
# GWR Cross-Validation (now delegates to .cv_run_folds)
# -----------------------------------------------------------------------------

#' K-fold cross-validation for GWR
#'
#' Refits a geographically weighted regression from scratch on each training
#' fold and scores it on the held-out fold, so the reported error is what the
#' model achieves at locations it did not see.  Reach for it whenever you need
#' a defensible accuracy figure for a GWR: the in-sample \eqn{R^2} that
#' \code{\link{model_metrics}()} reports on a \code{gwr_fit} is close to
#' meaningless, because a local regression with a small bandwidth can track the
#' training points almost exactly.  Bandwidth is re-selected per fold unless you
#' fix it with \code{bandwidth}, which keeps the selection itself inside the
#' cross-validation rather than tuning on the full data first.
#'
#' Folds default to spatial blocks (\code{\link{make_folds}(method =
#' "block_kfold")}), not random ones -- with autocorrelated data a random
#' split leaves a held-out point's neighbours in the training set and the score
#' comes back flattering.  Use \code{\link{cv_bayes}()} for the same treatment
#' of a Bayesian GP model, \code{\link{cv_rf}()} for a forest, and
#' \code{\link{compare_models_cv}()} to score several backends on one set of
#' folds.
#'
#' @param data_sf An sf object.
#' @param response_var Response column name.
#' @param predictor_vars Predictor column names.
#' @param folds Optional fold definitions: a \code{\link{make_folds}()} return
#'   value, or a bare list of \code{list(train =, test =)} pairs of
#'   \code{..row_id} values.  Train and test must be disjoint — a fold that
#'   trains on its own test rows is not a cross-validation split and is refused
#'   with an error — and IDs naming no row in the prepared data are dropped with
#'   a logged count (expected when rows were removed for missing values; a sign
#'   the folds came from other data when they were not).
#' @param k Number of folds. Default 5.
#' @param seed RNG seed. Default 123.
#' @param adaptive Logical; use adaptive bandwidth. Default TRUE.
#' @param bandwidth Optional bandwidth value.
#' @param kernel Kernel function type.
#' @param boundary Optional polygonal sf/sfc for CRS alignment.
#' @param pointize Geometry coercion strategy.
#' @param block_size Optional minimum block edge length for spatial CV blocks
#'   (projected CRS units).  Ensures blocks are at least as large as the
#'   spatial autocorrelation range.
#' @param auto_range Logical.  If \code{TRUE} and \code{folds} is \code{NULL},
#'   estimate the autocorrelation range and use it as the minimum block size.
#'   Default \code{FALSE}.
#' @param parallel Logical or positive integer.  If \code{TRUE},
#'   auto-detect the number of cores and fit folds in parallel via
#'   \code{parallel::mclapply()} (macOS / Linux; falls back to sequential
#'   on Windows).  If an integer > 1, use that many cores.  Default
#'   \code{FALSE} (sequential).
#' @return A list with \code{overall}, \code{fold_metrics},
#'   \code{predictions}, \code{folds}, \code{n_folds_attempted},
#'   \code{n_folds_succeeded}, \code{formula} and \code{adaptive}.  The two
#'   fold counts make a run where every fold failed visible in the return
#'   value rather than only in a warning, since \code{overall} is a
#'   well-formed all-\code{NA} row either way.
#' @family cross-validation
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
#'   cv <- cv_gwr(dat, "price", "elev", k = 3, bandwidth = 30)
#'   cv$overall
#'   cv$fold_metrics
#' }
#' }
#' @export
cv_gwr <- function(data_sf, response_var, predictor_vars,
                   folds = NULL, k = 5, seed = 123,
                   adaptive = TRUE, bandwidth = NULL,
                   kernel = c("bisquare", "gaussian", "tricube",
                              "boxcar", "exponential"),
                   boundary = NULL, pointize = "auto",
                   block_size = NULL, auto_range = FALSE,
                   parallel = FALSE) {
  if (!inherits(data_sf, "sf")) stop("cv_gwr(): `data_sf` must be an sf object.")
  if (!requireNamespace("GWmodel", quietly = TRUE))
    stop("cv_gwr(): package 'GWmodel' is required.", call. = FALSE)
  if (!requireNamespace("sp", quietly = TRUE))
    stop("cv_gwr(): package 'sp' is required (for GWmodel interop).", call. = FALSE)

  kernel <- match.arg(kernel)
  kernel <- .validate_kernel(kernel)

  if (!("..row_id" %in% names(data_sf))) data_sf$`..row_id` <- seq_len(nrow(data_sf))
  dat_sf <- prep_model_data(data_sf, response_var, predictor_vars, boundary, pointize)
  if (!("..row_id" %in% names(dat_sf)))
    stop("cv_gwr(): `prep_model_data()` must preserve `..row_id`.")

  keep_idx <- dat_sf$`..row_id`

  # Refuse folds built from other data before anything is fitted: the splits
  # are row IDs, so a wrong `folds` of the right size applies silently.
  .check_fold_probe(folds, data_sf, "cv_gwr")

  if (is.null(folds)) {
    message("cv_gwr(): no folds supplied \u2014 using spatial block k-fold CV (k=", k, ").")
    folds <- make_folds(dat_sf, k = k, method = "block_kfold",
                        seed = seed, boundary = boundary,
                        block_size = block_size, auto_range = auto_range,
                        response_var = response_var,
                        predictor_vars = predictor_vars)
  }

  remapped_folds <- .remap_folds(folds, keep_idx, k, seed)

  cleanup_cv <- .with_seed(seed)
  on.exit(cleanup_cv(), add = TRUE)

  # Per-fold fitting function: fit_gwr_model returns a gwr_fit S3 object
  # Data is already prepped by the outer cv_gwr() call; skip the redundant
  # prep_model_data() pass inside fit_gwr_model() for each fold.
  fit_one <- function(train_sf) {
    fit_gwr_model(
      data_sf = train_sf, response_var = response_var,
      predictor_vars = predictor_vars,
      adaptive = adaptive, bandwidth = bandwidth, kernel = kernel,
      .already_prepped = TRUE
    )
  }

  # Extra per-fold info: bandwidth

  fold_info_fn <- function(fit_obj, test_sf, y_true, y_hat) {
    list(bandwidth = fit_obj$info$bandwidth %||% NA_real_)
  }

  # GWR fits (p+1) local coefficients at every training location, so the

  # effective number of parameters is >> the global predictor count p.
  # Passing p here would yield a per-fold Adj_R² that drastically overstates
  # parsimony.  We set p = NULL so that per-fold Adj_R² is reported as NA,
  # consistent with the pooled metric and with cv_bayes().
  res <- .cv_run_folds(
    dat_sf = dat_sf, response_var = response_var,
    predictor_vars = predictor_vars,
    remapped_folds = remapped_folds, keep_idx = keep_idx,
    fit_one = fit_one, fold_info_fn = fold_info_fn,
    p = NULL, parallel = parallel, seed = seed
  )

  preds <- if (length(res$pred_rows)) do.call(rbind, res$pred_rows) else
    data.frame(`..row_id` = integer(), fold = integer(),
               y = numeric(), yhat = numeric(), y_train_mean = numeric())
  folds_df <- if (length(res$fold_stats)) as.data.frame(dplyr::bind_rows(res$fold_stats)) else
    data.frame(fold = integer(), n_train = integer(), n_test = integer(),
               n_pred = integer(),
               RMSE = numeric(), MAE = numeric(), MAPE = numeric(),
               SMAPE = numeric(), R2 = numeric(), Adj_R2 = numeric(),
               bandwidth = numeric())

  # The folds SUPPLIED (or built), not the ones that survived .remap_folds():
  # a fold dropped there is exactly the kind of thing this count exists to
  # make visible against n_folds_succeeded.
  n_attempted <- length(if (!is.null(folds$folds)) folds$folds else folds)
  n_succeeded <- length(res$fold_stats)
  if (n_succeeded == 0L && n_attempted > 0L) {
    why <- .cv_first_error_suffix(res)
    .log_warn("cv_gwr(): all %d folds failed to produce predictions; results are empty.%s",
              n_attempted, why)
    warning("cv_gwr(): all folds failed; cross-validation results contain no predictions.",
            why, call. = FALSE)
  } else if (n_succeeded < n_attempted) {
    .log_warn("cv_gwr(): %d of %d folds produced predictions.", n_succeeded, n_attempted)
  }

  list(overall = .cv_overall_metrics(preds), fold_metrics = folds_df,
       predictions = preds, folds = remapped_folds,
       # Reported so that a run where every fold failed is visible in the
       # return value, not only in a warning: `overall` is a well-formed
       # all-NA row either way.  cv_spatial() and cv_rf() report the same two
       # fields, so all four CV entry points share one shape.
       n_folds_attempted = n_attempted, n_folds_succeeded = n_succeeded,
       formula = deparse(stats::reformulate(predictor_vars, response_var)),
       adaptive = adaptive)
}


# -----------------------------------------------------------------------------
# Bayesian Cross-Validation (now delegates to .cv_run_folds)
# -----------------------------------------------------------------------------

#' K-fold cross-validation for the Bayesian spatial model
#'
#' Refits the Gaussian-process model of
#' \code{\link{fit_bayesian_spatial_model}()} on each training fold and scores
#' it on the held-out fold.  Beyond the point-prediction metrics the other CV
#' wrappers report, this one scores the whole predictive \emph{distribution}:
#' \code{predictive_coverage} says what fraction of held-out observations fell
#' inside the 50/80/95\% intervals, and \code{mean_CRPS} rates sharpness and
#' calibration together.  That is the reason to reach for it -- a Bayesian model
#' is usually chosen for its uncertainty, and only held-out coverage shows
#' whether those intervals are honest at locations the model has not seen.
#'
#' It is the most expensive wrapper in the package by a wide margin: every fold
#' is a full MCMC run.  Use few folds, and \code{parallel = TRUE} if you have
#' the cores.  For a cheap first pass on the same question, cross-validate a
#' forest with \code{\link{cv_rf}()} and come back here once the predictor set
#' has settled.
#'
#' @param data_sf An sf object.
#' @param response_var Response column name.
#' @param predictor_vars Predictor column names.
#' @param folds Optional fold definitions: a \code{\link{make_folds}()} return
#'   value, or a bare list of \code{list(train =, test =)} pairs of
#'   \code{..row_id} values.  Train and test must be disjoint — a fold that
#'   trains on its own test rows is not a cross-validation split and is refused
#'   with an error — and IDs naming no row in the prepared data are dropped with
#'   a logged count (expected when rows were removed for missing values; a sign
#'   the folds came from other data when they were not).
#' @param k Number of folds. Default 5.
#' @param seed RNG seed. Default 123.
#' @param boundary Optional polygonal sf/sfc for CRS alignment.
#' @param pointize Geometry coercion strategy.
#' @param fit_args Named list of extra arguments for fit_bayesian_spatial_model().
#'   A user-supplied \code{gp_k} is respected in every fold; when omitted, the
#'   GP rank is auto-selected per training fold.  \code{compute_loo},
#'   \code{boundary}, and \code{pointize} are always overridden by the CV
#'   internals.
#' @param summary "mean" or "median" for posterior predictions.
#' @param compute_pred_intervals Logical; compute predictive intervals.
#' @param coverage_levels Numeric vector of coverage levels.
#' @param block_size Optional minimum block edge length for spatial CV blocks
#'   (projected CRS units).
#' @param auto_range Logical.  If \code{TRUE} and \code{folds} is \code{NULL},
#'   estimate the autocorrelation range and use it as the minimum block size.
#'   Default \code{FALSE}.
#' @param parallel Logical or positive integer.  If \code{TRUE},
#'   auto-detect the number of cores and fit folds in parallel via
#'   \code{parallel::mclapply()} (macOS / Linux; falls back to sequential
#'   on Windows).  If an integer > 1, use that many cores.  Default
#'   \code{FALSE} (sequential).  Bayesian folds with full MCMC runs
#'   are the primary beneficiary of this option.
#' @return A list with \code{overall}, \code{fold_metrics},
#'   \code{predictions}, \code{folds}, \code{n_folds_attempted},
#'   \code{n_folds_succeeded}, \code{formula} and \code{predictive_coverage}.
#'   The two fold counts make a run where every fold failed visible in the
#'   return value rather than only in a warning.
#'   The \code{predictive_coverage} entries (one per
#'   \code{coverage_levels} value, plus \code{mean_CRPS}) are averages across
#'   folds \strong{weighted by each fold's \code{n_pred}}, because the per-fold
#'   values in \code{fold_metrics} are themselves means over that fold's test
#'   rows; an unweighted average would not be the pooled quantity when fold
#'   sizes differ, which for spatially blocked folds they routinely do.
#' @family cross-validation
#' @examples
#' \dontrun{
#' if (requireNamespace("brms", quietly = TRUE)) {
#'   library(sf)
#'   set.seed(1)
#'   n <- 60
#'   dat <- st_as_sf(
#'     data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), elev = rnorm(n)),
#'     coords = c("x", "y"), crs = 32632
#'   )
#'   dat$price <- 10 + 0.01 * st_coordinates(dat)[, 1] + 2 * dat$elev + rnorm(n)
#'   cv <- cv_bayes(dat, "price", "elev", k = 2,
#'                  fit_args = list(chains = 2, iter = 500))
#'   cv$overall
#'   cv$predictive_coverage  # coverage at 50/80/95% plus mean CRPS
#' }
#' }
#' @export
cv_bayes <- function(data_sf, response_var, predictor_vars,
                     folds = NULL, k = 5, seed = 123, boundary = NULL,
                     pointize = "auto", fit_args = list(),
                     summary = c("mean", "median"),
                     compute_pred_intervals = TRUE,
                     coverage_levels = c(0.50, 0.80, 0.95),
                     block_size = NULL, auto_range = FALSE,
                     parallel = FALSE) {
  summary <- match.arg(summary)
  if (!inherits(data_sf, "sf")) stop("cv_bayes(): `data_sf` must be an sf object.")

  if (!("..row_id" %in% names(data_sf))) data_sf$`..row_id` <- seq_len(nrow(data_sf))
  dat_sf <- prep_model_data(data_sf, response_var, predictor_vars, boundary, pointize)
  if (!("..row_id" %in% names(dat_sf)))
    stop("cv_bayes(): `prep_model_data()` must preserve `..row_id`.")

  keep_idx <- dat_sf$`..row_id`
  n_pred   <- length(predictor_vars)

  # Refuse folds built from other data before anything is fitted: the splits
  # are row IDs, so a wrong `folds` of the right size applies silently.
  .check_fold_probe(folds, data_sf, "cv_bayes")

  if (is.null(folds)) {
    message("cv_bayes(): no folds supplied \u2014 using spatial block k-fold CV (k=", k, ").")
    folds <- make_folds(dat_sf, k = k, method = "block_kfold",
                        seed = seed, boundary = boundary,
                        block_size = block_size, auto_range = auto_range,
                        response_var = response_var,
                        predictor_vars = predictor_vars)
  }

  remapped_folds <- .remap_folds(folds, keep_idx, k, seed)

  cleanup_cv <- .with_seed(seed)
  on.exit(cleanup_cv(), add = TRUE)

  # Per-fold fitting function: returns a bayesian_fit S3 object
  # Data is already prepped by the outer cv_bayes() call; skip the redundant
  # prep_model_data() pass inside fit_bayesian_spatial_model() for each fold.
  fit_one <- function(train_sf) {
    # NOTE: a user-supplied fit_args$gp_k is respected as-is.  When absent,
    # fit_bayesian_spatial_model()'s default (gp_k = NULL) auto-selects a
    # rank per training fold.  (Previously `gp_k = NULL` was passed through
    # modifyList(), which *removes* NULL elements and therefore silently
    # deleted the user's gp_k.)
    fold_fit_args <- modifyList(fit_args, list(
      compute_loo = FALSE, boundary = boundary,
      pointize = pointize,
      .already_prepped = TRUE
    ))
    do.call(fit_bayesian_spatial_model,
            c(list(data_sf = train_sf, response_var = response_var,
                   predictor_vars = predictor_vars), fold_fit_args))
  }

  # Extra per-fold info: gp_k, CRPS, predictive coverage
  fold_info_fn <- function(fit_obj, test_sf, y_true, y_hat) {
    extras <- list(
      gp_k       = as.integer(fit_obj$info$gp_k %||% NA_integer_),
      gp_n_basis = as.integer(fit_obj$info$gp_n_basis %||% NA_integer_),
      n_draws    = NA_integer_,
      CRPS    = NA_real_
    )
    # Pre-initialise coverage columns so every fold emits the same schema
    # even when the posterior-draw step fails for some folds; heterogeneous
    # per-fold columns would otherwise break row-binding of fold_metrics.
    for (cl in coverage_levels) {
      extras[[sprintf("coverage_%.0f", cl * 100)]] <- NA_real_
    }

    # Full posterior predictive draws for intervals + CRPS
    if (isTRUE(compute_pred_intervals)) {
      ppred_draws <- try(
        predict(fit_obj, newdata = test_sf, type = "predict", draws = TRUE),
        silent = TRUE
      )
      if (!inherits(ppred_draws, "try-error") && is.matrix(ppred_draws)) {
        extras$n_draws <- nrow(ppred_draws)

        # Coverage at each level
        for (cl in coverage_levels) {
          alpha <- (1 - cl) / 2
          lwr <- apply(ppred_draws, 2L, stats::quantile, probs = alpha)
          upr <- apply(ppred_draws, 2L, stats::quantile, probs = 1 - alpha)
          in_interval <- y_true >= lwr & y_true <= upr
          extras[[sprintf("coverage_%.0f", cl * 100)]] <- mean(in_interval, na.rm = TRUE)
        }

        # Empirical CRPS via the NRG (energy) form, vectorised over
        # observations — see .crps_energy() for the formula and reference.
        crps_per_obs <- .crps_energy(ppred_draws, y_true)
        extras$CRPS  <- mean(crps_per_obs, na.rm = TRUE)
      }
    }

    extras
  }

  res <- .cv_run_folds(
    dat_sf = dat_sf, response_var = response_var,
    predictor_vars = predictor_vars,
    remapped_folds = remapped_folds, keep_idx = keep_idx,
    fit_one = fit_one, fold_info_fn = fold_info_fn,
    predict_args = list(summary = summary),
    p = NULL, parallel = parallel, seed = seed
  )

  # Assemble predictions — add yhat_sd column via a second pass
  pred_rows_with_sd <- lapply(res$pred_rows, function(pr) {
    pr$yhat_sd <- NA_real_  # placeholder; could be enriched
    pr
  })

  preds <- if (length(pred_rows_with_sd)) do.call(rbind, pred_rows_with_sd) else
    data.frame(`..row_id` = integer(), fold = integer(),
               y = numeric(), yhat = numeric(), y_train_mean = numeric(),
               yhat_sd = numeric())
  # bind_rows() (rather than rbind) tolerates folds whose extras differ,
  # filling missing columns with NA instead of erroring.
  folds_df <- if (length(res$fold_stats)) as.data.frame(dplyr::bind_rows(res$fold_stats)) else
    data.frame(fold = integer(), n_train = integer(), n_test = integer(),
               n_pred = integer(),
               RMSE = numeric(), MAE = numeric(), MAPE = numeric(),
               SMAPE = numeric(), R2 = numeric(), Adj_R2 = numeric(),
               CRPS = numeric(), gp_k = integer(),
               gp_n_basis = integer(), n_draws = integer())

  # The folds SUPPLIED (or built), not the ones that survived .remap_folds():
  # a fold dropped there is exactly the kind of thing this count exists to
  # make visible against n_folds_succeeded.
  n_attempted <- length(if (!is.null(folds$folds)) folds$folds else folds)
  n_succeeded <- length(res$fold_stats)
  if (n_succeeded == 0L && n_attempted > 0L) {
    why <- .cv_first_error_suffix(res)
    .log_warn("cv_bayes(): all %d folds failed to produce predictions; results are empty.%s",
              n_attempted, why)
    warning("cv_bayes(): all folds failed; cross-validation results contain no predictions.",
            why, call. = FALSE)
  } else if (n_succeeded < n_attempted) {
    .log_warn("cv_bayes(): %d of %d folds produced predictions.", n_succeeded, n_attempted)
  }

  list(overall = .cv_overall_metrics(preds), fold_metrics = folds_df,
       predictions = preds, folds = remapped_folds,
       # Reported so that a run where every fold failed is visible in the
       # return value, not only in a warning -- which matters most here,
       # where a fold dies on Stan sampling rather than on bad input.
       n_folds_attempted = n_attempted, n_folds_succeeded = n_succeeded,
       formula = deparse(stats::reformulate(predictor_vars, response_var)),
       # Weighted by each fold's n_pred.  Per-fold coverage and CRPS are
       # themselves means over that fold's test rows, so an unweighted mean
       # across folds is NOT the pooled quantity as soon as fold sizes differ
       # -- and block_kfold tolerates a 3:1 imbalance before it even warns, so
       # they routinely do.  (The `overall` metrics above are already pooled
       # correctly, from the stacked prediction rows; this was the odd one out.)
       predictive_coverage = if (nrow(folds_df) > 0L && "CRPS" %in% names(folds_df)) {
         w_all <- if (!is.null(folds_df$n_pred))
           suppressWarnings(as.numeric(folds_df$n_pred)) else rep(1, nrow(folds_df))
         wm <- function(v) {
           v  <- suppressWarnings(as.numeric(v))
           ok <- is.finite(v) & is.finite(w_all) & w_all > 0
           if (!any(ok)) return(NA_real_)
           sum(v[ok] * w_all[ok]) / sum(w_all[ok])
         }
         cov_cols  <- grep("^coverage_", names(folds_df), value = TRUE)
         cov_means <- vapply(cov_cols, function(cn) wm(folds_df[[cn]]), numeric(1))
         c(as.list(cov_means), mean_CRPS = wm(folds_df$CRPS))
       } else NULL)
}


# -----------------------------------------------------------------------------
# Unified CV (model-agnostic)
# -----------------------------------------------------------------------------

#' Model-agnostic spatial cross-validation
#'
#' Run K-fold CV for any model that returns a \code{spatial_fit} object.
#' This is the extensibility point: to plug in a new model type, supply
#' a \code{fit_fn(train_sf)} that returns a \code{spatial_fit}.
#'
#' @param data_sf An sf object.
#' @param response_var Response column name.
#' @param predictor_vars Predictor column names.
#' @param fit_fn A function of one argument, the training slice of
#'   \code{data_sf}, returning a \code{spatial_fit} built with
#'   \code{\link{new_spatial_fit}()}.  It is called once per fold on the
#'   training rows only, so anything done inside it -- scaling, tuning, an inner
#'   variable sweep -- is already nested and leak-free.  The \code{subclass} it
#'   stamps must have a \code{predict.<subclass>()} method registered, because
#'   that is how each fold is scored.
#' @param folds Optional fold definitions: a \code{\link{make_folds}()} return
#'   value, or a bare list of \code{list(train =, test =)} pairs of
#'   \code{..row_id} values.  Train and test must be disjoint — a fold that
#'   trains on its own test rows is not a cross-validation split and is refused
#'   with an error — and IDs naming no row in the prepared data are dropped with
#'   a logged count.
#'   Built via \code{block_kfold} when \code{NULL}.
#' @param k Number of folds.
#' @param seed RNG seed.
#' @param boundary Optional boundary for fold construction.
#' @param pointize Geometry coercion strategy.
#' @param predict_args Extra arguments for predict().
#' @param fold_info_fn Optional function for per-fold extras.
#' @param p Number of predictors for Adj R² (NULL to skip).  Only meaningful
#'   for models with a fixed global parameter count; pass NULL for models
#'   with spatially varying coefficients (e.g. GWR).
#' @param block_size Optional minimum block edge length for spatial CV blocks
#'   (projected CRS units).
#' @param auto_range Logical.  If \code{TRUE} and \code{folds} is \code{NULL},
#'   estimate the autocorrelation range and use it as the minimum block size.
#'   Default \code{FALSE}.
#' @param parallel Logical or positive integer.  If \code{TRUE},
#'   auto-detect the number of cores and fit folds in parallel via
#'   \code{parallel::mclapply()} (macOS / Linux; falls back to sequential
#'   on Windows).  If an integer > 1, use that many cores.  Default
#'   \code{FALSE} (sequential).
#' @return A list with \code{overall}, \code{fold_metrics}, \code{predictions},
#'   \code{folds}, and the two fold counts \code{n_folds_attempted} and
#'   \code{n_folds_succeeded}.  The counts are reported deliberately: a
#'   \code{fit_fn} that fails on every fold otherwise looks like a successful
#'   run that happened to score \code{NA}, so compare them before trusting
#'   \code{overall}.  The \code{fold} column of \code{fold_metrics} and
#'   \code{predictions} carries the fold's index in the \code{folds} object
#'   that was supplied, so it lines up with
#'   \code{make_folds()$assignment$fold} even when some folds were unusable
#'   and dropped.
#' @family cross-validation
#' @seealso \code{\link{new_spatial_fit}()} for the constructor a \code{fit_fn}
#'   must use; \code{\link{cv_gwr}()}, \code{\link{cv_bayes}()} and
#'   \code{\link{cv_rf}()} for the built-in backends, which are thin wrappers
#'   over this function.
#' @examples
#' library(sf)
#' set.seed(1)
#' n <- 80
#' site <- st_as_sf(
#'   data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), elev = rnorm(n)),
#'   coords = c("x", "y"), crs = 32632
#' )
#' site$price <- 10 + 0.01 * st_coordinates(site)[, 1] + 2 * site$elev + rnorm(n)
#'
#' # 1. A fit_fn returning a spatial_fit of your own subclass.
#' lm_fit <- function(train_sf) {
#'   new_spatial_fit(
#'     subclass       = "lm_fit",
#'     engine         = lm(price ~ elev, st_drop_geometry(train_sf)),
#'     formula        = price ~ elev,
#'     response_var   = "price",
#'     predictor_vars = "elev",
#'     data_sf        = train_sf
#'   )
#' }
#'
#' # 2. The predict() method cv_spatial() scores each fold with. Without it
#' #    every fold fails and `overall` comes back all-NA.
#' predict.lm_fit <- function(object, newdata = NULL, ...) {
#'   if (is.null(newdata)) newdata <- object$data_sf
#'   as.numeric(stats::predict(object$engine, st_drop_geometry(newdata)))
#' }
#' registerS3method("predict", "lm_fit", predict.lm_fit)
#'
#' cv <- cv_spatial(site, "price", "elev", fit_fn = lm_fit, k = 3, seed = 1)
#' cv$overall
#' # Compare these before trusting the metrics above.
#' c(attempted = cv$n_folds_attempted, succeeded = cv$n_folds_succeeded)
#' @export
cv_spatial <- function(data_sf, response_var, predictor_vars,
                       fit_fn, folds = NULL, k = 5, seed = 123,
                       boundary = NULL, pointize = "auto",
                       predict_args = list(), fold_info_fn = NULL,
                       p = NULL, block_size = NULL,
                       auto_range = FALSE, parallel = FALSE) {
  if (!inherits(data_sf, "sf")) stop("cv_spatial(): `data_sf` must be an sf object.")
  if (!is.function(fit_fn)) stop("cv_spatial(): `fit_fn` must be a function.")

  if (!("..row_id" %in% names(data_sf))) data_sf$`..row_id` <- seq_len(nrow(data_sf))
  dat_sf <- prep_model_data(data_sf, response_var, predictor_vars, boundary, pointize)
  keep_idx <- dat_sf$`..row_id`

  # Refuse folds built from other data before anything is fitted: the splits
  # are row IDs, so a wrong `folds` of the right size applies silently.
  .check_fold_probe(folds, data_sf, "cv_spatial")

  if (is.null(folds)) {
    message("cv_spatial(): no folds supplied \u2014 using spatial block k-fold CV (k=", k, ").")
    folds <- make_folds(dat_sf, k = k, method = "block_kfold",
                        seed = seed, boundary = boundary,
                        block_size = block_size, auto_range = auto_range,
                        response_var = response_var,
                        predictor_vars = predictor_vars)
  }

  remapped_folds <- .remap_folds(folds, keep_idx, k, seed)

  cleanup_cv <- .with_seed(seed)
  on.exit(cleanup_cv(), add = TRUE)

  res <- .cv_run_folds(
    dat_sf = dat_sf, response_var = response_var,
    predictor_vars = predictor_vars,
    remapped_folds = remapped_folds, keep_idx = keep_idx,
    fit_one = fit_fn, fold_info_fn = fold_info_fn,
    predict_args = predict_args, p = p,
    parallel = parallel, seed = seed
  )

  preds <- if (length(res$pred_rows)) do.call(rbind, res$pred_rows) else
    data.frame(`..row_id` = integer(), fold = integer(),
               y = numeric(), yhat = numeric(), y_train_mean = numeric())
  folds_df <- if (length(res$fold_stats)) as.data.frame(dplyr::bind_rows(res$fold_stats)) else
    data.frame()

  # cv_gwr() and cv_bayes() both raise a real condition here; cv_spatial() used
  # to return an all-NA `overall` and an empty data.frame with nothing at R
  # condition level, so a fit_fn that failed on every fold looked like a
  # successful run that happened to score NA.
  # The folds SUPPLIED (or built), not the ones that survived .remap_folds():
  # a fold dropped there is exactly the kind of thing this count exists to
  # make visible against n_folds_succeeded.
  n_attempted <- length(if (!is.null(folds$folds)) folds$folds else folds)
  n_succeeded <- length(res$fold_stats)
  if (n_succeeded == 0L && n_attempted > 0L) {
    why <- .cv_first_error_suffix(res)
    .log_warn("cv_spatial(): all %d folds failed to produce predictions; results are empty.%s",
              n_attempted, why)
    warning("cv_spatial(): all folds failed; cross-validation results contain ",
            "no predictions.", why, call. = FALSE)
  } else if (n_succeeded < n_attempted) {
    .log_warn("cv_spatial(): %d of %d folds produced predictions.",
              n_succeeded, n_attempted)
  }

  list(overall = .cv_overall_metrics(preds),
       fold_metrics = folds_df, predictions = preds,
       folds = remapped_folds,
       n_folds_attempted = n_attempted, n_folds_succeeded = n_succeeded)
}
