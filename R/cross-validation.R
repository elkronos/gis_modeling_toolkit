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
#'   agree with \code{make_folds()$assignment$fold}.  Splits that already
#'   carry a \code{fold_id} (the \code{$folds} of a \code{cv_*()} result)
#'   keep it, when every one is a distinct whole number of at least 1.
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
    .warn_and_log(
      "cross-validation: falling back to random k-fold CV. For spatial data, use make_folds(method='block_kfold') to avoid optimistic performance estimates."
    )
    cleanup <- .with_seed(seed)
    on.exit(cleanup(), add = TRUE)
    assign_vec <- sample(rep(seq_len(k), length.out = length(keep_idx)))
    out <- vector("list", k)
    for (i in seq_len(k))
      out[[i]] <- list(train = keep_idx[assign_vec != i],
                       test  = keep_idx[assign_vec == i],
                       fold_id = i)
    attr(out, "dropped") <- data.frame(fold = integer(0), reason = character(0),
                                       stringsAsFactors = FALSE)
    attr(out, "orphans") <- keep_idx[0]
    attr(out, "n_unknown_ids") <- 0L
    return(out)
  }

  if (!is.null(folds$folds) && is.list(folds$folds)) {
    folds <- folds$folds
  }

  # Every split must carry `train` and `test` BY NAME.  area_of_applicability()
  # has always refused anything else and said which three shapes it takes; this
  # path read f$train / f$test straight, got NULL for both, and produced empty
  # folds -- so the row-coverage warning below fired and blamed folds "built on
  # a different or subsetted layer", which is not the cause, and the run
  # returned an all-NA `overall` with n_folds_succeeded = 0.  The realistic way
  # in is a fold list from another package: blockCV::cv_spatial()$folds_list
  # holds two UNNAMED vectors per fold (train first).  Its $folds_ids is a
  # label vector, which .folds_from_labels() accepts as it stands.
  if (!is.list(folds) || length(folds) == 0L)
    stop("cross-validation: `folds` is not a recognised fold object. Pass a ",
         "make_folds() result, a list of `list(train =, test =)` splits, or a ",
         "vector of fold labels, one per row.", call. = FALSE)
  bad <- which(!vapply(folds, function(z)
    is.list(z) && !is.null(z$train) && !is.null(z$test), logical(1)))
  if (length(bad)) {
    positional <- all(vapply(folds, function(z)
      is.list(z) && length(z) == 2L && is.null(names(z)), logical(1)))
    stop(sprintf(paste0(
      "cross-validation: %d of %d fold(s) carry no `train`/`test` element ",
      "(e.g. fold %s). `folds` must be a make_folds() result, a list of ",
      "`list(train =, test =)` splits, or a vector of fold labels, one per ",
      "row.%s"),
      length(bad), length(folds), paste(utils::head(bad, 3L), collapse = ", "),
      if (positional) paste0(
        " These splits look positional -- two unnamed vectors each -- so name ",
        "them `train` and `test`, or pass a label vector instead: ",
        "blockCV::cv_spatial() returns one as `$folds_ids`.") else ""),
      call. = FALSE)
  }

  # Carry the ORIGINAL fold index alongside each split.  Dropping an unusable
  # fold below removes an element from the list, and without this every later
  # fold would be silently renumbered, so fold_metrics$fold and
  # predictions$fold would no longer line up with make_folds()$assignment$fold.
  # A split that already carries a fold_id keeps it.  A cv_*() result's $folds
  # holds the splits that survived, each with the fold_id it was reported
  # under, so a dropped fold leaves a gap (1, 2, 4, 5).  Numbering those by
  # position relabelled every fold after the gap when the same splits were
  # handed to a second cv_*() -- its fold 3 was the first run's fold 4 --
  # while fold_separation() labels them by the fold_id they carry.  Carried
  # ids are kept when every split has a usable, distinct one; anything else
  # is numbered by position.
  fold_ids <- vapply(folds, function(f) {
    v <- f$fold_id
    if (is.numeric(v) && length(v) == 1L && is.finite(v) && v >= 1 &&
        v == round(v) && v <= .Machine$integer.max) as.integer(v)
    else NA_integer_
  }, integer(1))
  if (anyNA(fold_ids) || anyDuplicated(fold_ids))
    fold_ids <- seq_along(folds)

  # A hand-built `folds` list is documented as accepted by every cv_*(), and
  # two ways of getting it wrong went entirely unremarked.  Train and test
  # overlapping is not cross-validation at all -- the model is fitted and
  # scored on the same rows and the result is reported as a CV score (measured:
  # RMSE 0.50 against 0.97 for the same data properly split) -- so refuse it.
  # IDs that name no row were dropped silently by na.omit(match()), so a
  # typo'd or stale fold list looked like a clean run; those are counted and
  # reported, because prep_model_data() legitimately drops rows and the two
  # cases have to be told apart by the user, not by us.
  # Collected rather than counted: a row the folds name but the data does not
  # have appears once in each fold's test set and once in every OTHER fold's
  # training set, so summing per-fold counts reported it k times over -- five
  # missing rows read as 15, 20 or 50 at k = 3, 4 or 10, which is not a
  # number a user can compare with n_dropped or with the size of their data.
  unknown_ids <- vector("list", length(folds))
  for (j in seq_along(folds)) {
    f <- folds[[j]]
    ov <- intersect(f$train, f$test)
    if (length(ov))
      stop(sprintf(paste0("cross-validation: fold %d has %d row ID(s) in BOTH ",
                          "its train and test sets (e.g. %s). A fold that ",
                          "trains on its own test rows is not a ",
                          "cross-validation split; rebuild the folds with ",
                          "make_folds()."),
                   fold_ids[j], length(ov),
                   paste(utils::head(format(ov), 3L), collapse = ", ")),
           call. = FALSE)
    ent <- c(f$train, f$test)
    unknown_ids[[j]] <- ent[is.na(match(ent, keep_idx))]
  }
  unknown_total <- length(unique(unlist(unknown_ids, use.names = FALSE)))
  if (unknown_total > 0L)
    .log_info(paste0("cross-validation: the folds name %d row ID(s) that are ",
                     "not in the data, which were dropped. This is expected when ",
                     "rows were removed for missing values; if it is not, the ",
                     "folds were built on different data."),
              unknown_total)

  remapped <- lapply(seq_along(folds), function(j) {
    f <- folds[[j]]
    list(
      train = keep_idx[stats::na.omit(match(f$train, keep_idx))],
      test  = keep_idx[stats::na.omit(match(f$test, keep_idx))],
      fold_id = fold_ids[j]
    )
  })

  # The mirror of the check above.  Fold IDs naming no row are reported; rows
  # that no fold names were not -- a folds object built on site[1:45, ] and
  # applied to all 90 rows scored 45 of them and reported attempted =
  # succeeded = 3 with nothing said.  Those rows enter no fit and no score,
  # and a caller who did not intend that needs to hear it as an R condition.
  covered <- unique(unlist(lapply(remapped, function(f) c(f$train, f$test)),
                           use.names = FALSE))
  orphan  <- setdiff(keep_idx, covered)
  if (length(orphan))
    .warn_and_log(paste0("cross-validation: %d of %d rows in the data are ",
                         "named by no fold (e.g. row ID %s). They enter no ",
                         "training set and are never scored. If the folds ",
                         "were built on a different or subsetted layer, ",
                         "rebuild them with make_folds() on this one."),
                  length(orphan), length(keep_idx),
                  paste(utils::head(format(orphan), 3L), collapse = ", "))

  # An empty TRAINING set is just as fatal as an empty test set and was not
  # checked at all: buffered_loo with a buffer spanning the data, or
  # random_kfold on a single row, produces folds nothing can be fitted on.
  # They used to sail through here and get dropped one at a time inside
  # .cv_fit_one_fold(), so the only symptom was a generic "all folds failed"
  # warning at the very end.
  no_test  <- vapply(remapped, function(f) length(f$test) == 0L, logical(1))
  no_train <- vapply(remapped, function(f) length(f$train) < 2L, logical(1))
  if (any(no_test))
    .log_warn("cross-validation: %d fold(s) have empty test sets after remapping.",
              sum(no_test))
  if (any(no_train))
    .log_warn("cross-validation: %d fold(s) have fewer than 2 training rows after remapping and cannot be fitted.",
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

  # What was dropped, and which rows no fold named, travel with the result
  # rather than only through the warnings above: the cv_*() functions return
  # them as fold_status and orphan_rows.
  out <- remapped[!(no_test | no_train)]
  drop_idx <- which(no_test | no_train)
  attr(out, "dropped") <- data.frame(
    fold   = vapply(drop_idx, function(j) as.integer(remapped[[j]]$fold_id), integer(1)),
    # as.character(): ifelse() on a zero-length logical returns logical(0),
    # so with no fold dropped the column came back a different type from the
    # one the populated frame has.
    reason = as.character(ifelse(no_test[drop_idx],
                                 "empty test set after remapping",
                                 "fewer than 2 training rows after remapping")),
    stringsAsFactors = FALSE)
  attr(out, "orphans") <- orphan
  attr(out, "n_unknown_ids") <- as.integer(unknown_total)
  out
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
#' @param preds The stacked prediction rows.
#' @param metrics Optional user scoring function (see \code{cv_spatial()}),
#'   applied to the pooled finite pairs; its values are appended as columns.
#' @keywords internal
#' @noRd
.cv_overall_metrics <- function(preds, metrics = NULL) {
  ok <- is.finite(preds$y) & is.finite(preds$yhat)
  if (!any(ok)) {
    out <- data.frame(RMSE = NA_real_, MAE = NA_real_, MAPE = NA_real_,
                      SMAPE = NA_real_, R2 = NA_real_, Adj_R2 = NA_real_,
                      n_pred = 0L, n_MAPE = 0L, n_SMAPE = 0L)
    for (cn in .user_metric_names(metrics)) out[[cn]] <- NA_real_
    return(out)
  }

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
  out <- data.frame(RMSE = met$RMSE, MAE = met$MAE, MAPE = met$MAPE, SMAPE = met$SMAPE,
                    R2 = met$R2, Adj_R2 = met$Adj_R2, n_pred = met$n,
                    n_MAPE = met$n_MAPE, n_SMAPE = met$n_SMAPE,
                    stringsAsFactors = FALSE)
  if (!is.null(metrics)) {
    um <- .apply_user_metrics(metrics, preds$y, preds$yhat,
                              where = "the pooled predictions")
    if (is.null(um)) um <- .na_user_metrics(metrics)
    for (cn in names(um)) out[[cn]] <- um[[cn]]
  }
  out
}


# -----------------------------------------------------------------------------
# User-supplied scoring functions (the `metrics` argument of the cv_*())
# -----------------------------------------------------------------------------

#' The columns every fold-metrics and overall frame already has
#'
#' A user metric may not reuse one of these names: it would silently
#' overwrite the built-in value, and \code{compare_models_cv()} reads several
#' of them by name.  \code{mean_CRPS} is here because
#' \code{compare_models_cv()} writes it into the Bayesian row of
#' \code{overall} over a user column of that name.  The per-fold extras
#' (\code{CRPS}, \code{coverage_*}, \code{bandwidth}, ...) depend on the
#' backend and its arguments, so \code{.cv_fit_one_fold()} checks those as it
#' writes them.
#' @keywords internal
#' @noRd
.cv_reserved_metric_cols <- c("fold", "n_train", "n_test", "n_pred", "RMSE",
                              "MAE", "MAPE", "SMAPE", "R2", "Adj_R2", "n_MAPE",
                              "n_SMAPE", "model", "mean_CRPS")

#' Validate the `metrics` argument of the cv_*() functions
#'
#' @param metrics \code{NULL}, or a function of at least two arguments.
#' @param caller Name for the error message.
#' @return \code{metrics}, invisibly validated.
#' @keywords internal
#' @noRd
.check_metrics_fn <- function(metrics, caller) {
  if (is.null(metrics)) return(NULL)
  if (!is.function(metrics) || length(formals(metrics)) < 2L)
    stop(sprintf(paste0("%s(): `metrics` must be a function of two arguments, ",
                        "(y, yhat), returning a named numeric vector; got %s."),
                 caller,
                 if (is.function(metrics)) "a function of fewer than two arguments"
                 else paste0("an object of class ", class(metrics)[1L])),
         call. = FALSE)
  metrics
}

#' Coerce what a user metric returned into a named list of scalars
#'
#' Accepts a named numeric vector, a named list of length-one numerics, or a
#' one-row data frame.  Anything else, a missing or duplicated name, a
#' non-scalar element, or a name that collides with a built-in column is an
#' error: a scoring function that returns the wrong shape is a programming
#' mistake to surface, not a fold to skip.
#'
#' @param out The value returned by the user's function.
#' @return Named list of length-one numerics (\code{NA} where not finite).
#' @keywords internal
#' @noRd
.as_user_metric_list <- function(out) {
  if (is.data.frame(out)) {
    if (nrow(out) != 1L)
      stop("`metrics` returned a data frame with ", nrow(out),
           " rows; one row of named scalars is required.", call. = FALSE)
    out <- as.list(out)
  }
  if (is.atomic(out)) out <- as.list(out)
  if (!is.list(out) || length(out) == 0L)
    stop("`metrics` must return a named numeric vector (or a named list of ",
         "scalars); got an object of class ", class(out)[1L], ".", call. = FALSE)
  nm <- names(out)
  if (is.null(nm) || any(is.na(nm)) || any(!nzchar(nm)))
    stop("`metrics` must return a vector whose every element is named.",
         call. = FALSE)
  if (anyDuplicated(nm))
    stop("`metrics` returned duplicated names: ",
         paste(unique(nm[duplicated(nm)]), collapse = ", "), ".", call. = FALSE)
  clash <- intersect(nm, .cv_reserved_metric_cols)
  if (length(clash))
    stop("`metrics` returned names that are already columns of the metrics ",
         "frames: ", paste(clash, collapse = ", "), ". Use other names.",
         call. = FALSE)
  vals <- vector("list", length(out))
  for (i in seq_along(out)) {
    v <- out[[i]]
    if (length(v) != 1L || !(is.numeric(v) || is.logical(v)))
      stop("`metrics` must return one number per name; element '", nm[i],
           "' is not a single numeric value.", call. = FALSE)
    v <- as.numeric(v)
    vals[[i]] <- if (is.finite(v)) v else NA_real_
  }
  names(vals) <- nm
  vals
}

#' Apply the user's scoring function to one set of (y, yhat) pairs
#'
#' Only the pairs the built-in metrics use (both finite) reach the function,
#' so its columns are averages over the same rows as \code{RMSE}.
#' A function that throws is logged and contributes nothing there (the
#' column is then \code{NA} in that frame); one that returns the wrong shape
#' is an error, see \code{.as_user_metric_list()}.
#'
#' @param metrics The validated function.
#' @param y,yhat Observed and predicted values.
#' @param where Text for the log line ("fold 3", "the pooled predictions").
#' @return Named list of scalars, or \code{NULL} when the function threw.
#' @keywords internal
#' @noRd
.apply_user_metrics <- function(metrics, y, yhat, where) {
  ok  <- is.finite(y) & is.finite(yhat)
  out <- try(metrics(as.numeric(y[ok]), as.numeric(yhat[ok])), silent = TRUE)
  if (inherits(out, "try-error")) {
    .log_warn("cross-validation: `metrics` failed on %s: %s. Its columns are NA there.",
              where, .try_error_message(out))
    return(NULL)
  }
  .as_user_metric_list(out)
}

#' The names a user metric will produce, for typing an empty frame
#'
#' Asks the function on zero-length input; a function that cannot answer that
#' (throws, or returns the wrong shape) contributes no columns to an empty
#' frame, which is the one place its columns can be absent.
#' @keywords internal
#' @noRd
.user_metric_names <- function(metrics) {
  if (is.null(metrics)) return(character(0))
  out <- try(.as_user_metric_list(metrics(numeric(0), numeric(0))), silent = TRUE)
  if (inherits(out, "try-error")) character(0) else names(out)
}

#' An all-NA stand-in for a user metric that threw
#' @keywords internal
#' @noRd
.na_user_metrics <- function(metrics) {
  nm <- .user_metric_names(metrics)
  stats::setNames(as.list(rep(NA_real_, length(nm))), nm)
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


#' The per-fold metrics frame with no rows, typed
#'
#' Every \code{cv_*()} returns this shape when no fold produced a prediction,
#' so downstream code that subsets on a metric column works the same whether
#' the run succeeded or not.  \code{extra} names backend-specific columns
#' (\code{bandwidth} for GWR).
#'
#' @keywords internal
#' @noRd
.empty_fold_metrics <- function(extra = character(0), metrics = NULL) {
  base <- data.frame(fold = integer(), n_train = integer(), n_test = integer(),
                     n_pred = integer(),
                     RMSE = numeric(), MAE = numeric(), MAPE = numeric(),
                     SMAPE = numeric(), R2 = numeric(), Adj_R2 = numeric(),
                     n_MAPE = integer(), n_SMAPE = integer())
  for (e in c(extra, .user_metric_names(metrics))) base[[e]] <- numeric()
  base
}


#' Render the first fold error as a sentence to append to a warning
#'
#' Turns \code{.cv_run_folds()}'s \code{fit_errors} into \code{" First error:
#' <msg>"}, or \code{""} when no fold reported one.  Without this, "all 5 folds
#' failed" is the whole diagnosis a user gets when the backend package is not
#' installed.  The word "brms" (or "sp", or "ranger") never appears, even
#' though calling the fitter directly says so plainly.
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
#'   \item \strong{Independent of \code{sf_use_s2()}.}  On a lon/lat layer
#'     \code{st_centroid()} is spherical with s2 on and planar with it off,
#'     and the two differ by up to 5e-4 degrees on county polygons, far past
#'     the tolerance: folds built before \code{sf_use_s2(FALSE)} (a common
#'     workaround for invalid polygons), or saved and read in a session set
#'     the other way, were refused as "built from different data".  The
#'     centroid is now always the planar one, which is also the one that
#'     does not fail on invalid geometry.  Probes of this kind carry
#'     \code{kind = "planar_centroid"}; one without \code{kind} was taken
#'     under whatever setting was current, and is checked the same way.
#' }
#'
#' @param x An sf object carrying a \code{..row_id} column.
#' @param max_probe Maximum number of rows to record.  64 keeps the folds
#'   object small while making a same-size different-dataset collision
#'   effectively impossible.
#' @param legacy Take the centroid under the current \code{sf_use_s2()}
#'   setting, as probes without a \code{kind} were taken, to check one.
#' @return A list with \code{row_id} (the IDs as supplied), numeric \code{x}
#'   and \code{y}, \code{lonlat} (whether the coordinates are in EPSG:4326,
#'   i.e. the input carried a CRS), \code{kind}, and \code{points} (whether
#'   every geometry of \code{x} is a POINT, so that a refusal can say when
#'   folds built on a pointized copy meet the polygons they came from), or
#'   \code{NULL} when no probe can be taken.
#' @keywords internal
#' @noRd
.fold_row_probe <- function(x, max_probe = 64L, legacy = FALSE) {
  tryCatch({
    if (!inherits(x, "sf") || !("..row_id" %in% names(x)) || nrow(x) == 0L)
      return(NULL)
    ids  <- x[["..row_id"]]
    take <- unique(round(seq(1, nrow(x), length.out = min(nrow(x), max_probe))))
    all_points <- all(sf::st_geometry_type(x, by_geometry = TRUE) == "POINT")
    g    <- sf::st_geometry(x)[take]
    if (!all(sf::st_geometry_type(g, by_geometry = TRUE) == "POINT"))
      g <- if (isTRUE(legacy)) suppressWarnings(sf::st_centroid(g))
           else .planar_centroid(g)
    cr     <- suppressWarnings(sf::st_crs(x))
    lonlat <- FALSE
    if (!is.na(cr)) {
      g      <- suppressWarnings(sf::st_transform(g, 4326))
      lonlat <- TRUE
    }
    xy <- suppressWarnings(sf::st_coordinates(g))
    if (is.null(xy) || nrow(xy) != length(take)) return(NULL)
    list(row_id = ids[take], x = as.numeric(xy[, 1L]), y = as.numeric(xy[, 2L]),
         lonlat = lonlat,
         kind = if (isTRUE(legacy)) "session_centroid" else "planar_centroid",
         points = all_points)
  }, error = function(e) NULL)
}


#' The planar centroid of a geometry set, whatever \code{sf_use_s2()} says
#'
#' For \code{.fold_row_probe()}: s2 is switched off for the one call on a
#' lon/lat set and restored on exit.
#' @keywords internal
#' @noRd
.planar_centroid <- function(g) {
  if (isTRUE(sf::st_is_longlat(g)) && isTRUE(sf::sf_use_s2())) {
    suppressMessages(sf::sf_use_s2(FALSE))
    on.exit(suppressMessages(sf::sf_use_s2(TRUE)), add = TRUE)
  }
  suppressWarnings(sf::st_centroid(g))
}


#' Turn a vector of fold labels into train/test splits keyed by `..row_id`
#'
#' \code{area_of_applicability()} accepts a plain label vector (one per row);
#' every \code{cv_*()} documents the same three shapes but a vector reached
#' \code{.remap_folds()} and died with R's bare "$ operator is invalid for
#' atomic vectors".  Called right after \code{..row_id} is assigned, so the
#' i-th label maps to the i-th row's ID whatever rows are dropped later.
#'
#' @keywords internal
#' @noRd
.folds_from_labels <- function(folds, data_sf, caller) {
  if (is.null(folds) || is.list(folds) || !is.atomic(folds)) return(folds)
  n <- nrow(data_sf)
  if (length(folds) != n)
    stop(sprintf("%s(): `folds` has %d labels but the data has %d rows.",
                 caller, length(folds), n), call. = FALSE)
  # as.factor() orders its levels with sort(), which uses LC_COLLATE: labels
  # differing only in case or punctuation ("north"/"North") land in a
  # different order under C than under en_US, and since the fold NUMBER is the
  # level's position, fold_metrics$fold, predictions$fold and
  # fold_status$fold then name different groups on different machines.  The
  # partition is unaffected; the numbering was not reproducible.  sort() with
  # method = "radix" is always C-collation, so the numbering is now a property
  # of the labels alone.
  # That order is for character labels only.  Numbers sorted as strings run
  # "1", "10", "11", "12", "2", ..., so with ten or more numeric labels (site
  # IDs for leave-location-out) output fold 2 was the user's label 10, and
  # area_of_applicability(), which numbers the same vector numerically,
  # disagreed.  Numbers are ordered as numbers -- unique() after
  # as.character() because two doubles can print alike -- and a factor keeps
  # the level order its owner gave it.
  lv <- if (is.factor(folds)) levels(droplevels(folds))
        else if (is.numeric(folds) || is.logical(folds))
          unique(as.character(sort(unique(folds))))
        else sort(unique(as.character(folds)), method = "radix")
  f <- factor(as.character(folds), levels = lv)
  f <- droplevels(f)
  if (anyNA(f))
    stop(caller, "(): `folds` contains missing labels.", call. = FALSE)
  if (nlevels(f) < 2L)
    stop(caller, "(): `folds` must define at least two non-empty folds.",
         call. = FALSE)
  ids <- data_sf[["..row_id"]]
  lapply(levels(f), function(lv) {
    te <- ids[f == lv]
    list(train = setdiff(ids, te), test = te)
  })
}


#' Refuse fold splits that describe a different dataset
#'
#' Compares \code{folds$params$row_probe} against the data being
#' cross-validated (the caller's own \code{data_sf}, before preparation) over
#' whichever probe rows are present.  A fold set built by an older
#' version of this package carries no probe and is passed through unchecked,
#' as is one whose IDs cannot be matched (\code{NA} IDs) or whose coordinate
#' space cannot be compared (one side carried a CRS and the other did not).
#' When the probe was taken on POINT geometry and the data is not, or the
#' reverse, a mismatch is refused as that, not as folds from different data.
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
  # A probe with no `kind` predates the planar centroid and was taken under
  # the sf_use_s2() setting of its day; recompute the same way, so folds
  # saved by an older version are not refused for the change of method.
  now <- .fold_row_probe(data_sf, max_probe = nrow(data_sf),
                         legacy = is.null(probe$kind))
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
  # Folds built on coerce_to_points() of this very layer (or on the polygons,
  # handed a pointized copy) name the same rows, but each polygon's probe
  # point is its centroid while the pointized copy's is whatever point
  # `pointize` chose -- st_point_on_surface() under "auto" -- so every
  # non-convex feature "moved", and the refusal said the folds came from a
  # different dataset.  The locations cannot be compared across that change;
  # say that, and what to do.  A probe without `points` predates the field.
  if (bad > 0L && is.logical(probe$points) && length(probe$points) == 1L &&
      !is.na(probe$points) && !identical(probe$points, now$points))
    stop(sprintf(paste0("%s(): the supplied `folds` were built on %s geometry ",
                        "and this data has %s geometry, so their row locations ",
                        "cannot be matched (%d of %d checked row IDs sit at a ",
                        "different point): the point coerce_to_points() gives ",
                        "a polygon is in general not the centroid compared ",
                        "here. Build the folds with make_folds() on the layer ",
                        "passed here; it reduces polygons to points itself."),
                 caller,
                 if (probe$points) "POINT" else "non-POINT",
                 if (isTRUE(now$points)) "POINT" else "non-POINT",
                 bad, sum(ok)), call. = FALSE)
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


#' Give a cv_*() boundary that has no CRS one, once, naming the caller
#'
#' A \code{cv_*()} call hands \code{boundary} to \code{prep_model_data()},
#' which takes the data's CRS from it, and to \code{make_folds()}, which clips
#' the blocks to it.  Each interpreted a CRS-less boundary on its own: a
#' lon/lat one raised two R warnings, one from each, naming
#' \code{ensure_projected()} rather than the caller or the argument.  Called
#' with \code{to = NULL} before preparation, a boundary whose coordinates
#' look like lon/lat is taken as EPSG:4326 -- what \code{prep_model_data()}
#' would assume -- with one warning.  Called with the prepared data as
#' \code{to} before \code{make_folds()}, a boundary still without a CRS is
#' stamped with the data's, as \code{make_folds()} would stamp it, with one
#' warning.  Either way the next consumer finds a CRS and says nothing.
#'
#' @param boundary The caller's \code{boundary}, or \code{NULL}.
#' @param caller Calling function name, for the warning.
#' @param to \code{NULL}, or the prepared data whose CRS to stamp.
#' @return \code{boundary}, with a CRS where one was assumed.
#' @keywords internal
#' @noRd
.cv_boundary_crs <- function(boundary, caller, to = NULL) {
  if (!(inherits(boundary, "sf") || inherits(boundary, "sfc")) ||
      !is.na(sf::st_crs(boundary)))
    return(boundary)
  if (is.null(to)) {
    ll <- .looks_like_lonlat(boundary)
    if (!isTRUE(ll$lonlat)) return(boundary)
    .warn_and_log(
      paste0("%s(): `boundary` has no CRS; its coordinates look like lon/lat ",
             "(xmin=%.2f, xmax=%.2f, ymin=%.2f, ymax=%.2f), so it is taken as ",
             "EPSG:4326. Set the CRS explicitly with sf::st_crs() to suppress ",
             "this."),
      caller, ll$bb[["xmin"]], ll$bb[["xmax"]], ll$bb[["ymin"]], ll$bb[["ymax"]])
    return(sf::st_set_crs(boundary, 4326))
  }
  crs <- .crs_or_null(to)
  if (is.null(crs)) return(boundary)
  .transform_or_stamp(boundary, crs, what = "boundary", caller = caller)
}


#' Fit-predict a single CV fold
#'
#' Encapsulates the per-fold work so it can be called sequentially or in
#' parallel.  Returns \code{list(skip = <reason>)} when the fold is unusable
#' before any work starts or produces nothing scorable, or
#' \code{list(error = <message>)} when the fit or the prediction threw, so the
#' caller can report the cause as well as the count.
#'
#' @keywords internal
#' @noRd
.cv_fit_one_fold <- function(i, dat_sf, response_var, remapped_fold,
                             keep_idx, fit_one, fold_info_fn, predict_args,
                             p, metrics = NULL) {
  # The label to report: the fold's index in the ORIGINAL fold object, so that
  # dropping an unusable fold in .remap_folds() does not renumber the rest.
  # Falls back to the loop position for hand-built folds carrying no fold_id.
  fold_lab <- remapped_fold$fold_id %||% i

  tr_pos <- stats::na.omit(match(remapped_fold$train, keep_idx))
  te_pos <- stats::na.omit(match(remapped_fold$test, keep_idx))
  # A fold with nothing to fit on or nothing to score is skipped, and the
  # reason travels back with it so the caller's fold_status can say which.
  if (length(tr_pos) < 2L || length(te_pos) < 1L)
    return(list(skip = sprintf("%d training row(s) and %d test row(s) matched the data",
                               length(tr_pos), length(te_pos))))

  train_sf <- dat_sf[tr_pos, , drop = FALSE]
  test_sf  <- dat_sf[te_pos, , drop = FALSE]

  # Drop geometry once per fold to avoid redundant copies on wide data frames
  train_df <- sf::st_drop_geometry(train_sf)
  test_df  <- sf::st_drop_geometry(test_sf)

  # Fit model on training fold.
  #
  # The error TEXT is returned, not just logged: when every fold fails for the
  # same reason -- overwhelmingly the "package 'brms' is required" kind of
  # case -- the caller's "all N folds failed" warning is the only thing the
  # user sees, and without the cause it never names the missing backend.
  fit_obj <- try(fit_one(train_sf), silent = TRUE)
  if (inherits(fit_obj, "try-error")) {
    msg <- .try_error_message(fit_obj)
    .log_warn("cross-validation: fold %d fit failed; skipping. Cause: %s",
              fold_lab, msg)
    return(list(error = msg))
  }
  if (!inherits(fit_obj, "spatial_fit")) {
    msg <- sprintf("fit_fn() returned a %s, not a spatial_fit",
                   paste(class(fit_obj), collapse = "/"))
    .log_warn("cross-validation: fold %d did not return a spatial_fit; skipping.", fold_lab)
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
    .log_warn("cross-validation: fold %d predict failed; skipping. Cause: %s",
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
    .log_warn("cross-validation: fold %d predicted %d value(s) for %d test row(s); skipping.",
              fold_lab, length(y_hat), length(y_true))
    return(list(skip = sprintf("predict() returned %d value(s) for %d test row(s)",
                               length(y_hat), length(y_true))))
  }

  # Training-set mean: the baseline out-of-sample R² is measured against,
  # because it is the only null prediction available at prediction time.
  # (The held-out rows' own mean would be a null model that knows the test
  # data.  It fits them at least as well as any other constant, so it gives
  # the LOWER R², not a higher one; the choice is about using training
  # information only.)  model_metrics(newdata =) uses the same baseline.
  y_train <- train_df[[response_var]]
  y_train_mean <- mean(y_train[is.finite(y_train)], na.rm = TRUE)

  met <- .compute_reg_metrics(y_true, y_hat, p = p,
                              y_train_mean = y_train_mean)
  if (met$n == 0L)
    return(list(skip = "no finite (observed, predicted) pair among the test rows"))

  # Base fold stats
  fs <- data.frame(
    fold = fold_lab, n_train = length(tr_pos), n_test = length(y_true),
    n_pred = met$n,
    RMSE = met$RMSE, MAE = met$MAE, MAPE = met$MAPE, SMAPE = met$SMAPE,
    R2 = met$R2, Adj_R2 = met$Adj_R2,
    n_MAPE = met$n_MAPE, n_SMAPE = met$n_SMAPE, stringsAsFactors = FALSE
  )

  # Append model-specific per-fold info (bandwidth, gp_k, CRPS, coverage …).
  # A `..per_row` element is the exception: a data frame with one row per
  # test observation (cv_bayes() puts the posterior predictive SD there),
  # which goes into the prediction rows below rather than the fold stats.
  # A fold_info_fn that THROWS gets what a throwing user `metrics` function
  # gets: logged, its columns NA for this fold, the fold kept.  Its extras
  # used to be dropped without a word -- cv_bayes() given coverage_levels =
  # c(50, 80, 95) lost gp_k, n_draws, CRPS and every coverage column on every
  # fold with fold_status "ok" -- so the cause now also rides back as a note
  # for fold_status$message.  One that returns the wrong SHAPE is an error,
  # again as for `metrics`: a vector where a scalar belongs died later in
  # `fs[[cn]] <-` with R's "replacement has 0 rows".
  per_row <- NULL
  note    <- NULL
  if (!is.null(fold_info_fn)) {
    extra <- try(fold_info_fn(fit_obj, test_sf, y_true, y_hat), silent = TRUE)
    # A named vector -- the shape `metrics` accepts -- is taken as the list it
    # stands for.  It used to fall through every branch below: its columns
    # never appeared and fold_status said "ok".  as.list(NULL) is list(), so
    # a NULL return stays a no-op; an unnamed vector fails the naming check.
    if (!inherits(extra, "try-error") && (is.null(extra) || is.atomic(extra)))
      extra <- as.list(extra)
    if (inherits(extra, "try-error")) {
      note <- sprintf("fold_info_fn failed: %s", .try_error_message(extra))
      .log_warn("cross-validation: fold %s: %s. Its columns are NA there.",
                format(fold_lab), note)
    } else if (is.list(extra)) {
      if (is.data.frame(extra$..per_row)) {
        .check_fold_per_row(extra$..per_row)
        if (nrow(extra$..per_row) == length(y_true))
          per_row <- extra$..per_row
        else
          .log_warn(paste0("cross-validation: fold %s: fold_info_fn's `..per_row` ",
                           "has %d rows for %d test rows and was dropped."),
                    format(fold_lab), nrow(extra$..per_row), length(y_true))
      }
      extra$..per_row <- NULL
      .check_fold_extras(extra, names(fs))
      for (cn in names(extra)) fs[[cn]] <- extra[[cn]]
    } else {
      stop("cross-validation: `fold_info_fn` must return a named list (or a ",
           "named vector) of per-fold values; it returned an object of class ",
           class(extra)[1L], ".", call. = FALSE)
    }
  }

  # The user's scoring function, on the same finite pairs the built-in
  # metrics used.  Its pooled counterpart is applied in .cv_overall_metrics().
  # Its names are checked against the columns this fold already has, not only
  # against the static built-in list: the extras (CRPS, coverage_*, gp_k,
  # n_draws, bandwidth, or whatever a fold_info_fn returns) are only known
  # here, and a user metric of the same name overwrote them silently --
  # cv_bayes()'s predictive_coverage then reported the user's number.
  if (!is.null(metrics)) {
    um <- .apply_user_metrics(metrics, y_true, y_hat,
                              where = sprintf("fold %s", format(fold_lab)))
    if (is.null(um)) um <- .na_user_metrics(metrics)
    clash <- intersect(names(um), names(fs))
    if (length(clash))
      stop("`metrics` returned names that are already columns of the metrics ",
           "frames: ", paste(clash, collapse = ", "), ". Use other names.",
           call. = FALSE)
    for (cn in names(um)) fs[[cn]] <- um[[cn]]
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
  if (!is.null(per_row)) pr <- cbind(pr, per_row)

  list(pred_row = pr, fold_stat = fs, note = note)
}


#' Validate the \code{..per_row} a \code{fold_info_fn} returned
#'
#' Its columns are bound onto the fold's prediction rows, so each must be
#' named, once, with a name \code{predictions} does not already have.  A
#' clash was not caught: \code{dplyr::bind_rows()} renamed both copies
#' (\code{yhat...4}, \code{yhat...6}) with only a message, and
#' \code{overall} then found no \code{yhat} and came back all \code{NA}
#' with \code{n_pred = 0}, beside finite per-fold scores.
#'
#' @param per_row The \code{..per_row} data frame.
#' @return \code{invisible(NULL)}; called for the error.
#' @keywords internal
#' @noRd
.check_fold_per_row <- function(per_row) {
  nm <- names(per_row)
  if (!length(nm)) return(invisible(NULL))
  if (anyNA(nm) || !all(nzchar(nm)))
    stop("cross-validation: `fold_info_fn`'s `..per_row` must have a name ",
         "for every column.", call. = FALSE)
  if (anyDuplicated(nm))
    stop("cross-validation: `fold_info_fn`'s `..per_row` has duplicated ",
         "column names: ", paste(unique(nm[duplicated(nm)]), collapse = ", "),
         ".", call. = FALSE)
  clash <- intersect(nm, .cv_prediction_cols)
  if (length(clash))
    stop("cross-validation: `fold_info_fn`'s `..per_row` has columns that ",
         "`predictions` already has: ", paste(clash, collapse = ", "),
         ". Use other names.", call. = FALSE)
  invisible(NULL)
}

#' The columns every fold's prediction rows carry before any \code{..per_row}
#' @keywords internal
#' @noRd
.cv_prediction_cols <- c("..row_id", "fold", "y", "yhat", "y_train_mean")


#' Validate what a \code{fold_info_fn} returned, before it is written
#'
#' Every element (\code{..per_row} already removed) becomes one cell of the
#' fold's row of \code{fold_metrics}, so each must be named, once, with a
#' name that is not already a column, and hold one value.
#'
#' @param extra The list, without \code{..per_row}.
#' @param taken The column names the fold's row already has.
#' @return \code{invisible(NULL)}; called for the error.
#' @keywords internal
#' @noRd
.check_fold_extras <- function(extra, taken) {
  if (!length(extra)) return(invisible(NULL))
  nm <- names(extra)
  if (is.null(nm) || anyNA(nm) || !all(nzchar(nm)))
    stop("cross-validation: `fold_info_fn` must return a list whose every ",
         "element is named.", call. = FALSE)
  if (anyDuplicated(nm))
    stop("cross-validation: `fold_info_fn` returned duplicated names: ",
         paste(unique(nm[duplicated(nm)]), collapse = ", "), ".", call. = FALSE)
  clash <- intersect(nm, union(taken, .cv_reserved_metric_cols))
  if (length(clash))
    stop("cross-validation: `fold_info_fn` returned names that are already ",
         "columns of fold_metrics: ", paste(clash, collapse = ", "),
         ". Use other names.", call. = FALSE)
  bad <- nm[!vapply(extra, function(v) is.atomic(v) && length(v) == 1L,
                    logical(1))]
  if (length(bad))
    stop("cross-validation: `fold_info_fn` must return one value per name ",
         "(a data frame with one row per test row goes in `..per_row`); ",
         "element '", bad[1L], "' is not a single value.", call. = FALSE)
  invisible(NULL)
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
  # `m` is nrow(), an INTEGER, so m * m overflows above 46,340 draws and the
  # whole weight vector goes NA -- which made every CRPS, and with it
  # `mean_CRPS` and the calibration summary, silently NA behind one
  # "NAs produced by integer overflow" warning.  48,000 draws is an ordinary
  # cv_bayes(fit_args = list(chains = 4, iter = 13000)) run.  Do the
  # arithmetic in double.
  md <- as.numeric(m)
  weights <- (2 * seq_len(m) - md - 1) / (md * md)
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
    # A session-wide mc.cores opt-in caps the auto-detected count, the
    # convention brms and parallel::mclapply() itself follow.
    auto <- .sanitize_core_count(parallel::detectCores(logical = FALSE) - 1L)
    eff  <- min(auto, .sanitize_core_count(getOption("mc.cores", auto), auto))
  } else if (is.numeric(parallel) && length(parallel) == 1L &&
             !is.na(parallel) && parallel > 1) {
    eff <- .sanitize_core_count(parallel)
  } else {
    return(1L)
  }
  # More workers than cores is not parallelism, it is a fork bomb with a
  # memory bill: cap at the machine, and at two under R CMD check, which is
  # the limit CRAN's check farm enforces.
  n_machine <- .sanitize_core_count(parallel::detectCores(logical = TRUE),
                                    fallback = eff)
  if (eff > n_machine) {
    message("cv parallel: ", eff, " workers requested on a machine with ",
            n_machine, " cores; using ", n_machine, ".")
    eff <- n_machine
  }
  if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_")) && eff > 2L) eff <- 2L
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
#' \code{cv_gwr()} never asks for it (see the comment there).
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
#'   NA in place of a misleadingly favourable value.
#' @param parallel Logical or positive integer.  If \code{TRUE},
#'   auto-detect the number of cores; if an integer > 1, use that many
#'   cores; if \code{FALSE} (default), run sequentially.
#' @param n_cores \emph{Deprecated.}
#'   Explicit core count; overrides \code{parallel} when set.
#' @param seed Integer RNG seed, or \code{NULL} to draw the per-fold seeds
#'   from the session's RNG stream.  Either way one seed per fold is drawn in
#'   the parent process and applied inside the fold worker, so results depend
#'   on (seed, fold index) alone and \code{parallel = TRUE} reproduces
#'   \code{parallel = FALSE} exactly; with \code{NULL}, a \code{set.seed()}
#'   before the call makes both reproducible.
#' @param metrics Optional user scoring function, already validated by
#'   \code{.check_metrics_fn()}; applied per fold here.
#' @return List with pred_rows and fold_stats.
#' @keywords internal
#' @noRd
.cv_run_folds <- function(dat_sf, response_var, predictor_vars,
                          remapped_folds, keep_idx, fit_one,
                          fold_info_fn = NULL, predict_args = list(),
                          p = NULL, parallel = FALSE, n_cores = NULL,
                          seed = NULL, metrics = NULL) {
  cores   <- .resolve_n_cores(parallel, n_cores)
  n_folds <- length(remapped_folds)

  # Draw one seed per fold in the parent so that each fold's RNG stream is a
  # function of (seed, fold index) only -- never of the execution path.  This
  # is what makes parallel and sequential runs bit-identical.  Forked workers
  # cannot guarantee that on their own: mclapply() seeds each child from the
  # current time and process ID unless the L'Ecuyer-CMRG generator is in use,
  # so without this the parallel path is irreproducible for any fit_fn that
  # consumes RNG (see cv_spatial(), which accepts an arbitrary learner).
  # seed = NULL draws the per-fold seeds from the session's stream instead of
  # leaving the folds unseeded.  Unseeded forked workers are seeded by
  # mclapply() from the clock and the process ID, so `set.seed(1); cv_*(seed =
  # NULL, parallel = 2)` was irreproducible while the sequential call was --
  # against the README's unqualified promise.  Drawing from the caller's
  # stream (and advancing it, as any RNG-consuming call would) makes both
  # paths a function of the state set.seed() left, and identical to each other.
  fold_seeds <- if (is.null(seed)) {
    sample.int(.Machine$integer.max, n_folds)
  } else {
    cleanup_draw <- .with_seed(seed)
    on.exit(cleanup_draw(), add = TRUE)
    sample.int(.Machine$integer.max, n_folds)
  }

  # The label each fold reports: its index in the original fold object (see
  # .cv_fit_one_fold()), or its position here when it carries none.
  labels <- vapply(remapped_folds, function(f) as.integer(f$fold_id %||% NA), integer(1))
  labels[is.na(labels)] <- seq_along(remapped_folds)[is.na(labels)]

  fold_worker <- function(i) {
    # .with_seed() saves and restores .Random.seed, so seeding a fold never
    # leaks into the caller's RNG state.
    if (!is.na(fold_seeds[i])) {
      cleanup_i <- .with_seed(fold_seeds[i])
      on.exit(cleanup_i(), add = TRUE)
    }
    # Which fold is running, for a caller that collects the folds' warnings
    # and names the folds they came from (.cv_collect_fold_warnings()).
    outer_fold <- .cv_fold_state$fold
    .cv_fold_state$fold <- labels[i]
    on.exit(.cv_fold_state$fold <- outer_fold, add = TRUE)
    .cv_fit_one_fold(
      i = i, dat_sf = dat_sf, response_var = response_var,
      remapped_fold = remapped_folds[[i]], keep_idx = keep_idx,
      fit_one = fit_one, fold_info_fn = fold_info_fn,
      predict_args = predict_args, p = p, metrics = metrics
    )
  }

  if (cores > 1L) {
    message(sprintf("cv: running %d folds in parallel on %d cores.",
                    n_folds, cores))
    # A warning raised inside a forked child dies with the child: R
    # conditions do not cross the fork, so fit_gwr_model()'s documented
    # integer-response warning, raised in every fold, reached nobody under
    # parallel = 2 while the sequential run showed all four.  Collect them in
    # the worker and re-raise in the parent, once per distinct message.  A
    # warning of one of the package's own classes comes back as the condition
    # it was, class and fields kept, so that a handler written for the class
    # sees it in a parallel run as in a sequential one (.cv_relay_warning()).
    # An error that escapes .cv_fit_one_fold() -- the documented shape errors
    # of `metrics` and `fold_info_fn` -- aborts a sequential run.  In a worker
    # it became that fold's try-error, and the run carried on without it.  It
    # is caught here, carried back, and raised again in the parent below.
    caught_worker <- function(i) {
      caught <- list()
      res  <- tryCatch(
        withCallingHandlers(
          fold_worker(i),
          warning = function(w) {
            caught[[length(caught) + 1L]] <<- .cv_relay_warning(w)
            invokeRestart("muffleWarning")
          }),
        error = function(e) list(escaped = conditionMessage(e)))
      # .cv_fit_one_fold() returns NULL for an unusable fold; NULL cannot
      # carry an attribute, so wrap the pair instead.
      list(res = res, fold_warnings = caught)
    }
    # mc.preschedule = FALSE: one fork per fold.  Prescheduled, mclapply()
    # hands each core a chunk of folds, copies one fold's try-error to every
    # fold of its chunk and returns NULL for every fold of a core that dies,
    # so a fold that succeeded was discarded, or reported with another
    # fold's error, for sharing a core with a failure.  k forks cost little
    # next to k model fits, and the per-fold seeds are drawn above, so the
    # results do not change.
    results <- parallel::mclapply(
      seq_along(remapped_folds), caught_worker, mc.cores = cores,
      mc.preschedule = FALSE
    )
    relayed <- do.call(c, lapply(results, function(z)
      if (!inherits(z, "try-error")) as.list(z$fold_warnings)))
    said <- character(0)
    for (w in relayed) {
      classed <- inherits(w, "condition")
      # The two warnings in which predict() on a GWR counts locations are
      # raised every time: the caller adds their counts up, and two folds can
      # count the same.
      if (!(classed && inherits(w, .gwr_counted_classes))) {
        key <- if (classed) paste(class(w)[1L], conditionMessage(w), sep = "\r") else w
        if (key %in% said) next
        said <- c(said, key)
      }
      if (classed) warning(w) else warning(w, call. = FALSE)
    }
    results <- lapply(results, function(z) if (inherits(z, "try-error")) z else z$res)
    esc <- which(vapply(results, function(z)
      is.list(z) && !inherits(z, "try-error") && !is.null(z$escaped), logical(1)))
    if (length(esc)) {
      j <- esc[1L]
      stop(sprintf("%s (fold %s, raised in a parallel worker)",
                   results[[j]]$escaped,
                   format(remapped_folds[[j]]$fold_id %||% j)), call. = FALSE)
    }
  } else {
    results <- lapply(seq_along(remapped_folds), fold_worker)
  }

  # One status per fold, in the folds' own order, before anything is filtered:
  # what each fold did is the diagnosis a caller needs when "3 of 5 folds
  # produced predictions" is all the console kept.  mclapply() hands back a
  # try-error OBJECT when a child errors, and NULL for a child that was killed
  # (out of memory, a segfault) -- a worker error, not a skip, and it was
  # reported as "skipped" and left out of fit_errors; a fold that threw comes
  # back as list(error = <message>); one skipped before fitting or after
  # predicting as list(skip = <reason>); a successful one carries pred_row and
  # fold_stat, and a note when its fold_info_fn failed.  `$` (not `[[`)
  # throughout, because a successful fold's list has no "error" element and
  # `[[` would abort.
  status <- character(n_folds); msg <- character(n_folds)
  for (i in seq_len(n_folds)) {
    z <- results[[i]]
    if (inherits(z, "try-error")) {
      status[i] <- "worker_error"; msg[i] <- .try_error_message(z)
    } else if (is.null(z)) {
      status[i] <- "worker_error"
      msg[i] <- paste0("the parallel worker delivered no result (it was ",
                       "killed or crashed, for example for lack of memory)")
    } else if (!is.null(z$error)) {
      status[i] <- "error"; msg[i] <- z$error
    } else if (!is.null(z$skip)) {
      status[i] <- "skipped"; msg[i] <- z$skip
    } else {
      status[i] <- "ok"; msg[i] <- z$note %||% ""
    }
  }
  fold_status <- data.frame(fold = labels, status = status, message = msg,
                            stringsAsFactors = FALSE)
  if (any(status == "worker_error"))
    .log_warn("cv: %d fold(s) failed in a parallel worker: %s",
              sum(status == "worker_error"),
              paste(unique(msg[status == "worker_error"]), collapse = "; "))
  fit_errors <- msg[status %in% c("worker_error", "error")]
  results <- results[status == "ok"]

  pred_rows  <- lapply(results, `[[`, "pred_row")
  fold_stats <- lapply(results, `[[`, "fold_stat")

  list(pred_rows = pred_rows, fold_stats = fold_stats,
       fit_errors = as.character(fit_errors), fold_status = fold_status)
}


# The fold a cross-validation is fitting and predicting at this moment (its
# label), while .cv_run_folds() runs it; NULL otherwise.
.cv_fold_state <- new.env(parent = emptyenv())

#' A fold's warning in the form a parallel worker hands back
#'
#' R conditions do not cross a fork, so a worker collects what its fold
#' raised and the parent raises it again.  A warning of one of the package's
#' own classes (\code{spatialkit_*}) goes back as the condition itself, with
#' its fields and the fold's label as \code{fold}; its call is dropped, as
#' the package's own carry none.  Any other warning goes back as its text:
#' a condition of another package can hold a call or a backtrace with the
#' fold's data in it.
#'
#' @param w The warning condition.
#' @return The condition, or its message.
#' @keywords internal
#' @noRd
.cv_relay_warning <- function(w) {
  if (!any(startsWith(class(w), "spatialkit_"))) return(conditionMessage(w))
  w$call <- NULL
  if (is.null(w$fold)) w$fold <- .cv_fold_state$fold
  w
}


#' Collect the warnings a Bayesian fit raises in each fold
#'
#' \code{fit_bayesian_spatial_model()} raises up to four classed warnings
#' about a fit: more basis functions than the rows support
#' (\code{"spatialkit_gp_basis_size"}), a variogram range below what the
#' basis resolves (\code{"spatialkit_gp_resolution"}), a sampler that did not
#' converge (\code{"spatialkit_convergence"}) and a basis the fitted model
#' shows to be too coarse (\code{"spatialkit_gp_adequacy"}).  None says which
#' fold it came from, and a cross-validation passed each on once per fold.
#' Wrapped in this they are muffled as they arrive and returned with the
#' label of the fold that raised them, for
#' \code{.cv_raise_fold_warnings()} to raise once per class.  Each fold's
#' own line stays in the log.
#'
#' @param expr The code that runs the folds.
#' @return A list: \code{value}, the value of \code{expr}, and \code{said},
#'   one entry per warning with \code{class}, \code{fold}, \code{message}
#'   and, for a convergence warning, \code{failed} (the checks that failed).
#' @keywords internal
#' @noRd
.cv_collect_fold_warnings <- function(expr) {
  said <- list()
  note <- function(w) {
    # The fold a relayed warning was marked with in its worker, or the one
    # running now.  A warning from no fold is not this function's to hold
    # back, and is passed on as it is.
    fold <- as.integer(w$fold %||% .cv_fold_state$fold %||% NA_integer_)[1L]
    if (is.na(fold)) return(invisible(NULL))
    said[[length(said) + 1L]] <<- list(
      class   = intersect(.cv_fold_warning_classes, class(w))[1L],
      fold    = fold,
      message = conditionMessage(w),
      failed  = as.character(w$failed %||% character(0)))
    invokeRestart("muffleWarning")
  }
  value <- withCallingHandlers(expr,
                               spatialkit_gp_basis_size = note,
                               spatialkit_gp_resolution = note,
                               spatialkit_convergence   = note,
                               spatialkit_gp_adequacy   = note)
  list(value = value, said = said)
}

.cv_fold_warning_classes <- c("spatialkit_gp_basis_size", "spatialkit_gp_resolution",
                              "spatialkit_convergence", "spatialkit_gp_adequacy")

# The checks of the sampler, as fit_bayesian_spatial_model() names them on a
# convergence warning (`failed`) and as a message names them.
.cv_convergence_checks <- c(divergent = "divergent transitions", rhat = "R-hat",
                            ess = "effective sample size", ebfmi = "E-BFMI")


#' Raise the folds' Bayesian warnings once per class, naming the folds
#'
#' The other half of \code{.cv_collect_fold_warnings()}.  One warning per
#' class, of that class, under the caller's name.  Each names the folds,
#' carries them as \code{folds} and every fold's own text as
#' \code{fold_messages} (named by fold), and quotes the first fold's text.
#'
#' The convergence warning covers the folds that were scored: a fold whose
#' sampler did not converge is scored like any other, and that is what the
#' warning is about.  \code{cv_bayes()} knows them from
#' \code{fold_metrics$convergence_ok}; a caller without that column knows
#' them from the warnings collected.
#'
#' @param caller Function name the messages carry.
#' @param said What \code{.cv_collect_fold_warnings()} collected.
#' @param scored Labels of the folds that produced scores.
#' @param n_run Number of folds run.
#' @param not_converged Labels of the scored folds whose fit reports
#'   \code{convergence_ok = FALSE}, or \code{NULL} when the caller has no
#'   such record.
#' @return \code{invisible(NULL)}; called for the warnings.
#' @keywords internal
#' @noRd
.cv_raise_fold_warnings <- function(caller, said, scored, n_run,
                                    not_converged = NULL) {
  of_class <- function(cl) {
    x <- Filter(function(z) identical(z$class, cl), said)
    x <- x[!duplicated(vapply(x, function(z) z$fold, integer(1)))]
    x[order(vapply(x, function(z) z$fold, integer(1)))]
  }
  own_text <- function(z) sub("^fit_bayesian_spatial_model\\(\\): ", "", z$message)
  texts_of <- function(x) stats::setNames(
    vapply(x, function(z) z$message, character(1)),
    vapply(x, function(z) as.character(z$fold), character(1)))
  others <- function(x) if (length(x) > 1L)
    " The log has the other folds' text, and the condition carries every fold's as `fold_messages`."
  else ""
  in_cv_bayes <- identical(caller, "cv_bayes")

  # ---- The sampler ----
  conv <- of_class("spatialkit_convergence")
  conv_folds <- vapply(conv, function(z) z$fold, integer(1))
  bad <- sort(unique(as.integer(
    if (is.null(not_converged)) intersect(conv_folds, scored) else not_converged)))
  if (length(bad)) {
    conv <- conv[conv_folds %in% bad]
    failed <- stats::setNames(lapply(conv, function(z) z$failed),
                              vapply(conv, function(z) as.character(z$fold), character(1)))
    failed <- failed[lengths(failed) > 0L]
    any_failed <- unique(unlist(failed, use.names = FALSE))
    where <- if (in_cv_bayes) "`fit_args`" else "the learner's call"
    detail <- if (length(failed) == length(bad)) {
      # Every such fold said which checks failed: name them, and give the
      # remedy of each check that failed somewhere.
      remedy <- c(
        if ("divergent" %in% any_failed)
          sprintf("raise `control = list(adapt_delta = ...)` in %s for the divergent transitions", where),
        if (any(c("rhat", "ess") %in% any_failed))
          sprintf("raise `iter` in %s for %s", where,
                  paste(c(if ("rhat" %in% any_failed) "R-hat",
                          if ("ess" %in% any_failed) "the effective sample size"),
                        collapse = " and ")),
        if ("ebfmi" %in% any_failed)
          paste0("a low E-BFMI comes from the shape of the posterior, so a longer ",
                 "warmup can help and otherwise the model or its priors need changing"))
      sprintf(" The checks that failed: %s. The log has each fold's values.%s",
              paste(sprintf("fold %s: %s", names(failed), vapply(failed, function(f)
                paste(.cv_convergence_checks[intersect(names(.cv_convergence_checks), f)],
                      collapse = ", "), character(1))), collapse = "; "),
              if (length(remedy)) {
                r <- paste(remedy, collapse = "; ")
                sprintf(" %s%s.", toupper(substr(r, 1L, 1L)), substr(r, 2L, nchar(r)))
              } else "")
    } else {
      sprintf(paste0(" The log says which checks failed (divergent transitions, ",
                     "R-hat, effective sample size, E-BFMI): raise ",
                     "`control = list(adapt_delta = ...)` in %s for divergent ",
                     "transitions, `iter` for R-hat and the effective sample size."),
              where)
    }
    .warn_classed(
      "spatialkit_convergence",
      paste0("%s(): the sampler did not converge in %d of %d fold(s) (fold %s); ",
             "their scores enter `overall` like the others but come from ",
             "unreliable posteriors.%s%s"),
      caller, length(bad), length(scored), paste(bad, collapse = ", "),
      if (in_cv_bayes) " fold_metrics$convergence_ok marks them." else "",
      detail,
      data = list(folds = bad, failed = failed, fold_messages = texts_of(conv)))
  }

  # ---- The basis against the rows ----
  big <- of_class("spatialkit_gp_basis_size")
  if (length(big))
    .warn_classed(
      "spatialkit_gp_basis_size",
      paste0("%s(): in %d of %d fold(s) (fold %s) the training rows are fewer ",
             "than one per ten basis functions. Fold %d: %s%s"),
      caller, length(big), n_run,
      paste(vapply(big, function(z) z$fold, integer(1)), collapse = ", "),
      big[[1L]]$fold, own_text(big[[1L]]), others(big),
      data = list(folds = vapply(big, function(z) z$fold, integer(1)),
                  fold_messages = texts_of(big)))

  # ---- The basis against a variogram fitted before sampling ----
  res <- of_class("spatialkit_gp_resolution")
  if (length(res))
    .warn_classed(
      "spatialkit_gp_resolution",
      paste0("%s(): in %d of %d fold(s) (fold %s) the Gaussian-process basis ",
             "may not resolve the autocorrelation range. Fold %d: %s%s"),
      caller, length(res), n_run,
      paste(vapply(res, function(z) z$fold, integer(1)), collapse = ", "),
      res[[1L]]$fold, own_text(res[[1L]]), others(res),
      data = list(folds = vapply(res, function(z) z$fold, integer(1)),
                  fold_messages = texts_of(res)))

  # ---- The basis against the fitted model ----
  adq <- of_class("spatialkit_gp_adequacy")
  if (length(adq))
    .warn_classed(
      "spatialkit_gp_adequacy",
      paste0("%s(): in %d of %d fold(s) (fold %s) the fitted model shows a ",
             "Gaussian-process basis too coarse for the structure in the ",
             "data%s. Fold %d: %s%s"),
      caller, length(adq), n_run,
      paste(vapply(adq, function(z) z$fold, integer(1)), collapse = ", "),
      if (in_cv_bayes) "; fold_metrics$gp_adequacy has each fold's status" else "",
      adq[[1L]]$fold, own_text(adq[[1L]]), others(adq),
      data = list(folds = vapply(adq, function(z) z$fold, integer(1)),
                  fold_messages = texts_of(adq)))
  invisible(NULL)
}


#' The per-fold status of a cross-validation, dropped folds included
#'
#' \code{.remap_folds()} drops folds with an empty test set or fewer than two
#' training rows before the runner sees them and records which on an
#' attribute; the runner reports every fold it ran.  This merges the two into
#' one frame in fold order, so \code{n_folds_attempted - n_folds_succeeded}
#' always has names and reasons beside it.
#'
#' @param remapped_folds The list \code{.remap_folds()} returned.
#' @param res The list \code{.cv_run_folds()} returned.
#' @return A data.frame with \code{fold}, \code{status} (\code{"ok"},
#'   \code{"error"}, \code{"skipped"}, \code{"worker_error"} or
#'   \code{"dropped"}) and \code{message}.
#' @keywords internal
#' @noRd
.cv_fold_status <- function(remapped_folds, res) {
  ran <- res$fold_status
  if (is.null(ran))
    ran <- data.frame(fold = integer(0), status = character(0), message = character(0),
                      stringsAsFactors = FALSE)
  dropped <- attr(remapped_folds, "dropped")
  if (is.data.frame(dropped) && nrow(dropped))
    ran <- rbind(ran, data.frame(fold = as.integer(dropped$fold), status = "dropped",
                                 message = as.character(dropped$reason),
                                 stringsAsFactors = FALSE))
  ran <- ran[order(ran$fold), , drop = FALSE]
  rownames(ran) <- NULL
  ran
}


#' Warn when some folds, but not all, produced no predictions
#'
#' \code{overall} is pooled over the folds that succeeded, so a run that lost
#' a fold reports a score over fewer rows -- and not a random few: in spatial
#' block CV the fold that fails is usually the hardest extrapolation block
#' (the only rows of a factor level, a region no training fold covers).  This
#' used to be a log line only, which \code{tryCatch()},
#' \code{expect_warning()} and \code{options(warn = 2)} never see.  Folds
#' dropped before fitting are left out of the warning: \code{.remap_folds()}
#' has already raised one for them.  The all-folds-failed case is the
#' callers' own, louder warning.
#'
#' @param caller Function name the message carries.
#' @param res The list \code{.cv_run_folds()} returned.
#' @param preds The stacked prediction rows \code{overall} is pooled from.
#' @param n_rows Rows in the prepared data.
#' @param n_attempted,n_succeeded The fold counts the caller returns.
#' @keywords internal
#' @noRd
.cv_warn_failed_folds <- function(caller, res, preds, n_rows,
                                  n_attempted, n_succeeded) {
  st  <- res$fold_status
  bad <- if (is.data.frame(st)) st[st$status != "ok", , drop = FALSE] else NULL
  if (is.null(bad) || nrow(bad) == 0L) {
    if (n_succeeded < n_attempted)
      .log_warn("%s(): %d of %d folds produced predictions.",
                caller, n_succeeded, n_attempted)
    return(invisible(NULL))
  }
  .warn_and_log(paste0(
    "%s(): %d of %d fold(s) failed (%s), so `overall` pools the other folds ",
    "only and covers %d of the %d rows. A fold that fails is often the ",
    "hardest to predict (a region or a factor level no training fold ",
    "covers), so `overall` may flatter the model; `fold_status` gives each ",
    "fold's cause."),
    caller, nrow(bad), n_attempted,
    paste(sprintf("fold %s: %s", bad$fold, bad$status), collapse = ", "),
    sum(is.finite(preds$y) & is.finite(preds$yhat)), n_rows)
}


#' Is this the warning \code{.cv_warn_failed_folds()} raises for \code{caller}?
#'
#' For a caller that runs \code{cv_spatial()} many times and reports the
#' consequence in its own terms (\code{select_features_forward()}), so the
#' two cannot drift apart.
#' @keywords internal
#' @noRd
.is_failed_folds_warning <- function(w, caller = "cv_spatial") {
  startsWith(conditionMessage(w), paste0(caller, "(): ")) &&
    grepl("^[^:]+: [0-9]+ of [0-9]+ fold\\(s\\) failed \\(",
          conditionMessage(w))
}


#' The fold list without \code{.remap_folds()}'s bookkeeping attributes
#'
#' The dropped-fold frame, the orphan IDs and the unknown-ID count ride on the
#' remapped list so the runner and \code{.cv_fold_status()} can read them; the
#' \code{folds} element the cv_*() functions return carries them as their own
#' elements instead, so the list itself is handed back plain.
#'
#' @keywords internal
#' @noRd
.bare_folds <- function(x) {
  attr(x, "dropped") <- NULL
  attr(x, "orphans") <- NULL
  attr(x, "n_unknown_ids") <- NULL
  x
}


# -----------------------------------------------------------------------------
# Spatial autocorrelation range estimation
# -----------------------------------------------------------------------------

#' Nearest-neighbour distance from each feature of \code{query} to \code{data}
#'
#' Uses \code{sf::st_nearest_feature()} so the search is indexed and avoids a
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


#' Fit a spatial trend and an exponential-plus-nugget covariance by REML
#'
#' The detrending step of \code{estimate_sac_range(detrend = "reml")}.
#' \code{nlme::gls()} estimates the trend coefficients and the covariance
#' parameters (range, nugget proportion, variance) jointly by residual maximum
#' likelihood, so the trend is fitted by generalised least squares under the
#' fitted correlation and the range is the REML estimate.  No variogram is
#' fitted to residuals, which is what removes the residual-variogram bias
#' (Lark, Cullis and Welham 2006).  The fit is
#' \eqn{O(n^3)}, so it runs on at most \code{max_n} rows (a seeded random
#' subsample when there are more) and the trend coefficients are then applied
#' to every row.
#'
#' A spatial correlation structure cannot take a zero distance between
#' distinct observations, so the fit uses one row a location.  Where a
#' location has several rows (repeat visits), the one used is drawn at
#' random under \code{seed}.  It is not "the first": the rows arrive sorted
#' by coordinates and then by the response, so the first row of a location is
#' its smallest response, and a fit to those is a fit to the minimum of each
#' location's visits (with four visits of unit noise the intercept came out
#' 1.04 low and the nugget at 0.41 for a truth of 1).  It is not the mean of
#' the visits either: the mean of \eqn{m} visits carries a \eqn{1/m} share of
#' the variance between visits, and the nugget fitted to means is short by
#' that factor.  The rows are put in a fixed order first (coordinates, then
#' the other columns), so neither draw depends on the order they came in.
#'
#' The REML likelihood of a short range on a large layer has a valley between
#' the optimum and the plateau where the fit is the same as with no
#' correlation, and a start in that valley can climb either side.  So the
#' start is chosen from the likelihood, and the fit kept is the best of the
#' starts tried, by its REML log-likelihood, never the first that converges.
#' The profile REML log-likelihood is evaluated on a grid (range parameter:
#' the median nearest-neighbour distance of the rows used, times 1/4 to 64 in
#' powers of two, capped at a third of the extent; nugget proportion 0.05,
#' 0.3, 0.6 and 0.9), which costs one Cholesky factorisation a point, and
#' \code{gls()} starts from the best grid point.  When that fit fails, or
#' ends below the grid's best value, the three extent-based starts (a tenth,
#' a thirtieth and a third of the extent) are fitted as well and the fit
#' with the largest log-likelihood is kept.  A fit is not discarded for a
#' nugget near 0 or near 1: the likelihood decides.
#'
#' @param mf A model frame (trend formula already applied, incomplete rows
#'   excluded) with the projected coordinates in columns \code{.sac_x} and
#'   \code{.sac_y}.
#' @param fml The trend formula.
#' @param extent A length scale of the layer, used to bound starting ranges.
#' @return \code{NULL} when no fit converged, else a list with \code{beta}
#'   (named trend coefficients), \code{range} (the exponential range
#'   parameter), \code{nugget_prop}, \code{sigma2}, \code{n_used},
#'   \code{subsampled}, \code{n_colocated} (the rows left out because another
#'   row of their location was drawn), \code{pair_floor} (the 30th-shortest
#'   distance between the rows the fit used, after the draw of one row a
#'   location and the subsample: the distance a range has to reach for 30
#'   pairs of those points to lie inside it), \code{logLik} (the REML
#'   log-likelihood of the fit kept), \code{logLik_null} (that of the same
#'   trend with independent errors, \code{NA} when it could not be fitted)
#'   and \code{n_fits} (how many starts were fitted).
#' @keywords internal
#' @noRd
.reml_trend <- function(mf, fml, extent, max_n = 400L, seed = 123L) {
  if (!requireNamespace("nlme", quietly = TRUE)) return(NULL)
  d <- mf
  # A fixed order of the rows: coordinates, then every other column.  Both
  # draws below index rows, and gls() sums over them in order.
  sort_cols <- c(".sac_x", ".sac_y", setdiff(names(d), c(".sac_x", ".sac_y")))
  ord <- tryCatch(
    do.call(order, c(unname(lapply(d[sort_cols], function(x)
      if (is.numeric(x) || is.logical(x)) as.numeric(x) else as.character(x))),
      list(method = "radix"))),
    error = function(e) NULL)
  if (length(ord) == nrow(d)) d <- d[ord, , drop = FALSE]
  dup <- duplicated(d[, c(".sac_x", ".sac_y")])
  if (any(dup) || nrow(d) > max_n) {
    cleanup <- .with_seed(seed)
    on.exit(cleanup(), add = TRUE)
  }
  # One row a location, drawn at random: the first row of each location in a
  # random order of the rows.  Nothing is drawn when no location repeats, so
  # a layer of distinct locations gets the subsample it always got.
  n_colocated <- 0L
  if (any(dup)) {
    n_repeated <- sum(!duplicated(d[dup, c(".sac_x", ".sac_y"), drop = FALSE]))
    perm <- sample.int(nrow(d))
    keep <- sort(perm[!duplicated(d[perm, c(".sac_x", ".sac_y"), drop = FALSE])])
    n_colocated <- nrow(d) - length(keep)
    d <- d[keep, , drop = FALSE]
    .log_info(paste0("estimate_sac_range(): %d location(s) have more than one ",
                     "row. The REML fit uses one row a location, drawn at ",
                     "random under `seed` (%d of %d row(s)); the trend it fits ",
                     "is then removed from every row."),
              n_repeated, nrow(d), nrow(d) + n_colocated)
  }
  n <- nrow(d)
  subsampled <- FALSE
  if (n > max_n) {
    d <- d[sample.int(n, max_n), , drop = FALSE]
    subsampled <- TRUE
  }
  n <- nrow(d)
  if (n < 30L) return(NULL)
  # The REML range is fitted to these point pairs, not to a binned variogram,
  # so the shortest range it can identify is set by how many of them are
  # close: 30 pairs is the usual minimum per lag (Journel and Huijbregts).
  D  <- as.matrix(stats::dist(d[, c(".sac_x", ".sac_y")]))
  pd <- D[upper.tri(D)]
  pair_floor <- if (length(pd) >= 30L) sort(pd, partial = 30L)[30L] else NA_real_

  # One gls() fit from a starting range and nugget proportion; NULL when it
  # fails or its parameters are not usable.
  fit_from <- function(r0, nug0) {
    fit <- tryCatch(
      withCallingHandlers(
        nlme::gls(fml, data = d,
                  correlation = nlme::corExp(value = c(r0, nug0),
                                             form = ~ .sac_x + .sac_y,
                                             nugget = TRUE),
                  method = "REML"),
        warning = function(w) {
          .log_info("estimate_sac_range(): nlme::gls() warned during REML detrending: %s",
                    conditionMessage(w))
          invokeRestart("muffleWarning")
        }),
      error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    cs <- try(stats::coef(fit$modelStruct$corStruct, unconstrained = FALSE),
              silent = TRUE)
    if (inherits(cs, "try-error") || !all(is.finite(cs))) return(NULL)
    rng <- as.numeric(cs[["range"]]); nug <- as.numeric(cs[["nugget"]])
    s2  <- as.numeric(fit$sigma)^2
    ll  <- suppressWarnings(as.numeric(stats::logLik(fit)))
    if (!is.finite(rng) || rng <= 0 || !is.finite(s2) || s2 <= 0) return(NULL)
    list(beta = stats::coef(fit), range = rng, nugget_prop = nug, sigma2 = s2,
         logLik = if (length(ll) == 1L && is.finite(ll)) ll else NA_real_)
  }

  # The profile REML log-likelihood of (range, nugget proportion), up to a
  # constant: trend and variance profiled out.  NULL when the design cannot
  # be built, and then every start is fitted and compared by gls()'s own
  # log-likelihood.
  X <- try(stats::model.matrix(fml, d), silent = TRUE)
  y <- try(stats::model.response(stats::model.frame(fml, d)), silent = TRUE)
  profile_ll <- NULL
  if (!inherits(X, "try-error") && !inherits(y, "try-error") &&
      is.numeric(y) && nrow(X) == n && length(y) == n) {
    profile_ll <- function(r, nug) {
      R <- (1 - nug) * exp(-D / r)
      diag(R) <- 1
      U <- tryCatch(chol(R), error = function(e) NULL)
      if (is.null(U)) return(-Inf)
      Xs  <- backsolve(U, X, transpose = TRUE)
      ys  <- backsolve(U, y, transpose = TRUE)
      qrx <- qr(Xs)
      df  <- n - qrx$rank
      rss <- sum(qr.resid(qrx, ys)^2)
      if (df < 1L || !is.finite(rss) || rss <= 0) return(-Inf)
      rd  <- abs(diag(qr.R(qrx)))[seq_len(qrx$rank)]
      -0.5 * (df * log(rss / df) + 2 * sum(log(diag(U))) + 2 * sum(log(rd)))
    }
  }

  floor_r <- sqrt(.Machine$double.eps)
  legacy  <- unique(pmax(extent * c(1 / 10, 1 / 30, 1 / 3), floor_r))
  starts  <- lapply(legacy, function(r) c(r, 0.1))
  grid_best <- -Inf
  if (!is.null(profile_ll)) {
    nn <- D; diag(nn) <- Inf
    nn_med <- stats::median(apply(nn, 1L, min))
    r_grid <- unique(pmax(pmin(nn_med * 2^(-2:6), extent / 3), floor_r))
    g <- expand.grid(r = r_grid, nug = c(0.05, 0.3, 0.6, 0.9))
    g$ll <- mapply(profile_ll, g$r, g$nug)
    if (any(is.finite(g$ll))) {
      j <- which.max(g$ll)
      grid_best <- g$ll[j]
      starts <- c(list(c(g$r[j], g$nug[j])), starts)
    }
  }

  best <- NULL
  n_fits <- 0L
  for (s in starts) {
    f <- fit_from(s[1L], s[2L])
    n_fits <- n_fits + 1L
    if (!is.null(f)) {
      # Fits are compared on one scale: the profile value at the fitted
      # parameters where it is available, gls()'s log-likelihood otherwise.
      f$score <- if (!is.null(profile_ll)) profile_ll(f$range, f$nugget_prop)
                 else f$logLik
      if (is.null(best) || isTRUE(f$score > best$score)) best <- f
    }
    # The grid start is fitted first.  When its fit ends at least as high as
    # the grid's best value there is nothing a start in the valley could add.
    if (n_fits == 1L && is.finite(grid_best) && !is.null(best) &&
        isTRUE(best$score >= grid_best - 0.01)) break
  }
  if (is.null(best)) return(NULL)

  ll_null <- tryCatch(
    suppressWarnings(as.numeric(stats::logLik(nlme::gls(fml, data = d, method = "REML")))),
    error = function(e) NA_real_)
  list(beta = best$beta, range = best$range, nugget_prop = best$nugget_prop,
       sigma2 = best$sigma2, n_used = n, subsampled = subsampled,
       n_colocated = as.integer(n_colocated),
       pair_floor = pair_floor, logLik = best$logLik,
       logLik_null = if (length(ll_null) == 1L && is.finite(ll_null)) ll_null else NA_real_,
       n_fits = n_fits)
}


#' Does the empirical variogram decrease with distance over its shorter lags?
#'
#' A semivariance that falls as distance grows leaves the range unidentified:
#' a model with a monotone rise to a sill is being fitted to a curve that has
#' none.  The net change in \code{gamma} across the lags in the shorter half
#' of the variogram, weighted by the pairs supporting each step, is compared
#' with the mean semivariance over those lags; a net fall of more than
#' \code{tol} of that mean is a decrease.
#'
#' The steps are taken between lag classes of equal width, as many as there
#' are bins in the shorter half.  Bins that are already equally spaced are
#' their own classes.  Geometrically spaced bins (\code{lags = "adaptive"}
#' in \code{estimate_sac_range()}) crowd the short lags, and their pair
#' counts grow with the square of the lag, so the pair weights would let the
#' last few steps speak for all of them; pooling the bins that share a class
#' (pair-weighted) puts the steps back on an even footing in distance, and a
#' variogram that rises steeply over its first bins and then levels off is
#' not called decreasing.
#'
#' The tolerance is fixed while the noise in the short-lag bins grows as the
#' sample shrinks, so small samples trip it on ordinary fields.  The refusal
#' is the conservative outcome (geometric blocks, with a warning), so it is
#' left as it is.  An unremoved trend is \emph{not} what produces this
#' shape.  A trend makes the variogram rise without reaching a sill, which
#' the over-cutoff rejection catches, so the message aimed at this case must
#' not say "trend".
#'
#' @param vg An empirical variogram from \code{gstat::variogram()}.
#' @return \code{TRUE} for a net decrease, \code{FALSE} otherwise (including
#'   when fewer than \code{min_bins} lag classes fall in the shorter half).
#' @keywords internal
#' @noRd
.variogram_decreasing <- function(vg, frac = 0.5, min_bins = 4L, tol = 0.15) {
  if (!is.data.frame(vg) || !all(c("dist", "gamma", "np") %in% names(vg)))
    return(FALSE)
  vg <- vg[is.finite(vg$gamma) & is.finite(vg$dist) & is.finite(vg$np) & vg$np > 0, ,
           drop = FALSE]
  if (nrow(vg) < min_bins) return(FALSE)
  vg <- vg[order(vg$dist), , drop = FALSE]
  upper <- frac * max(vg$dist)
  short <- vg[vg$dist <= upper, , drop = FALSE]
  if (nrow(short) < min_bins || !(upper > 0)) return(FALSE)
  # Equal-width lag classes over the shorter half; a class takes the
  # pair-weighted mean of the bins in it.
  k   <- nrow(short)
  cls <- pmin(k, pmax(1L, as.integer(ceiling(short$dist / (upper / k) - 1e-9))))
  np    <- as.numeric(tapply(short$np, cls, sum))
  gamma <- as.numeric(tapply(short$np * short$gamma, cls, sum)) / np
  if (length(np) < min_bins) return(FALSE)
  w   <- pmin(utils::head(np, -1L), utils::tail(np, -1L))
  if (!any(w > 0)) return(FALSE)
  net <- sum(w * diff(gamma)) / sum(w) * (length(np) - 1L)
  scale <- stats::weighted.mean(gamma, np)
  is.finite(net) && is.finite(scale) && scale > 0 && net < -tol * scale
}


#' The k nearest other rows of a coordinate matrix
#'
#' \pkg{FNN}'s kd-tree when it is installed, a blocked brute-force search
#' otherwise.  Meant for distinct locations: on exact duplicates the kd-tree
#' can drop a twin.
#'
#' @param xy Numeric matrix, two columns.
#' @param k Number of neighbours; capped at \code{nrow(xy) - 1}.
#' @return A list with \code{index} and \code{dist}, each \code{nrow(xy)} by
#'   \code{k}, neighbours in order of distance.
#' @keywords internal
#' @noRd
.sac_knn <- function(xy, k) {
  n <- nrow(xy)
  k <- min(as.integer(k), n - 1L)
  if (k < 1L)
    return(list(index = matrix(integer(0), n, 0L), dist = matrix(numeric(0), n, 0L)))
  if (requireNamespace("FNN", quietly = TRUE)) {
    nn <- FNN::get.knn(xy, k = k)
    return(list(index = matrix(nn$nn.index, n, k), dist = matrix(nn$nn.dist, n, k)))
  }
  idx <- matrix(0L, n, k); dst <- matrix(0, n, k)
  for (s in seq(1L, n, by = 500L)) {
    r <- s:min(n, s + 499L)
    d <- sqrt(outer(xy[r, 1L], xy[, 1L], "-")^2 + outer(xy[r, 2L], xy[, 2L], "-")^2)
    d[cbind(seq_along(r), r)] <- Inf
    for (i in seq_along(r)) {
      o <- order(d[i, ])[seq_len(k)]
      idx[r[i], ] <- o
      dst[r[i], ] <- d[i, o]
    }
  }
  list(index = idx, dist = dst)
}


#' Lag boundaries spaced geometrically from the scale of the closest pairs
#'
#' The first bin runs from zero to the distance within which
#' \code{max(50, n / 4)} pairs of distinct locations lie (read from
#' nearest-neighbour distances, and capped at a quarter of the cutoff); the
#' rest are spaced geometrically from there to the cutoff.  Equal bins out to
#' half the extent are far wider than the spacing of close pairs on a dense
#' layer, so structure shorter than one bin is invisible to them and is
#' counted as nugget.
#'
#' The size of the first bin is a trade.  A larger one, \code{max(150, n /
#' 2)} pairs, is steadier on a layer with no structure, but the test for
#' spatial structure already refuses those (white noise at n = 300, 400
#' draws: a finite range in 8 with this rule and in 10 with the larger bin,
#' a nugget under a tenth of the sill in 3 and in 2), and on structured
#' layers the larger bin loses the short lags the nugget is read from: with
#' an effective range of 60 on 300 points over a 1000 m square the nugget
#' came out under 0.05 (truth 0.2) in 24 of 34 identified fits against 16,
#' and at 80 points in 19 of 24 against 14 of 25.
#'
#' The ratio between boundaries is at least \code{sqrt(1 + min_np / m0)},
#' which is what gives a ring \code{min_np} pairs when pairs grow with the
#' square of the distance, so the number of bins past the first is
#' \code{n_lags} or fewer, and at least 5.
#'
#' @param xy Coordinates of the points the variogram is computed on.
#' @param cutoff_dist The largest lag.
#' @param n_lags The most bins wanted past the first.
#' @return Increasing upper boundaries; the last is \code{cutoff_dist}.
#' @keywords internal
#' @noRd
.sac_lag_boundaries <- function(xy, cutoff_dist, n_lags = 15L, min_np = 30L) {
  u  <- unique(xy)
  nu <- nrow(u)
  m0 <- max(50, nrow(xy) / 4)
  lo <- cutoff_dist / 4
  if (nu >= 3L) {
    k  <- min(max(10L, as.integer(ceiling(4 * m0 / nu)) + 2L), nu - 1L)
    kd <- sort(as.numeric(.sac_knn(u, k)$dist))
    kd <- kd[is.finite(kd) & kd > 0]
    # Each pair is listed from both of its ends.  The margin keeps the pair
    # the boundary was read from inside the first bin whatever the rounding
    # of its distance, so the bins hold the same pairs in any unit.
    if (length(kd))
      lo <- min(kd[min(length(kd), as.integer(ceiling(2 * m0)))] * (1 + 1e-9),
                cutoff_dist / 4)
  }
  nb <- floor(log(cutoff_dist / lo) / log(sqrt(1 + min_np / m0)))
  nb <- max(5L, min(as.integer(n_lags), as.integer(nb)))
  exp(seq(log(lo), log(cutoff_dist), length.out = nb + 1L))
}


#' Merge neighbouring variogram bins until each holds enough pairs
#'
#' From the shortest lag up, bins are pooled until the pool holds
#' \code{min_np} pairs; a short tail joins the bin before it.  The pooled
#' semivariance and distance are the pair-weighted means, which is what the
#' classical estimator gives for the union of the bins, so no pair is lost.
#'
#' @param vg An empirical variogram (\code{np}, \code{dist}, \code{gamma}).
#' @return The variogram with merged rows, ordered by distance.
#' @keywords internal
#' @noRd
.sac_merge_bins <- function(vg, min_np = 30L) {
  if (!is.data.frame(vg) || nrow(vg) < 2L) return(vg)
  vg <- vg[order(vg$dist), , drop = FALSE]
  if (all(vg$np >= min_np)) return(vg)
  grp <- integer(nrow(vg)); g <- 1L; acc <- 0
  for (i in seq_len(nrow(vg))) {
    grp[i] <- g
    acc <- acc + vg$np[i]
    if (acc >= min_np) { g <- g + 1L; acc <- 0 }
  }
  if (acc > 0 && g > 1L) grp[grp == g] <- g - 1L
  np  <- as.numeric(tapply(vg$np, grp, sum))
  out <- vg[!duplicated(grp), , drop = FALSE]
  out$dist  <- as.numeric(tapply(vg$np * vg$dist,  grp, sum)) / np
  out$gamma <- as.numeric(tapply(vg$np * vg$gamma, grp, sum)) / np
  out$np    <- np
  rownames(out) <- NULL
  out
}


#' Take the pairs at zero separation out of the first bin of a variogram
#'
#' \pkg{gstat} counts two rows that share a location as a pair in the first
#' lag bin, of the all-pairs variogram and of each directional one.  Such a
#' pair is not at a lag.  The semivariance jumps at zero (between visits to
#' one location it is the variance between visits; at any positive
#' separation the variation between locations is added), and a bin that
#' mixes the two sides of the jump gets a mean distance and a semivariance
#' in between, which a model reads as a rise from a low nugget.  On 100
#' locations visited 5 times, with independent location effects of variance
#' 1 and a variance between visits of 0.3, the fitted nugget was 0 and the
#' range twice the spacing of the locations, where the variogram at positive
#' lags is flat at 1.3.  So the first bin keeps its pairs at positive
#' separations only, and the semivariance at zero separation is reported on
#' its own (\code{zero_lag_semivariance}).
#'
#' @param vg An empirical variogram from \code{gstat::variogram()}, rows by
#'   increasing distance.
#' @param zero_pairs,zero_ss The number of pairs at zero separation and the
#'   sum of their half squared differences, in the units of \code{vg$gamma}.
#' @param first_bound The upper edge of the first bin.
#' @return \code{vg} with the pairs taken out of its first row, or without
#'   that row when it held no others.  Unchanged when there are no such
#'   pairs, or when the first row cannot hold them (fewer pairs than that,
#'   or a remainder that would lie outside the bin or below zero), which
#'   would mean \pkg{gstat} had not put them there.
#' @keywords internal
#' @noRd
.sac_drop_zero_lag <- function(vg, zero_pairs, zero_ss, first_bound) {
  if (!is.data.frame(vg) || !nrow(vg) || !isTRUE(zero_pairs > 0)) return(vg)
  np1  <- as.numeric(vg$np[1L])
  rest <- np1 - as.numeric(zero_pairs)
  if (!is.finite(rest) || rest < 0) return(vg)
  if (rest == 0) {
    # Nothing but pairs at zero separation: there is no bin left.
    if (isTRUE(vg$dist[1L] <= 1e-12 * first_bound)) vg <- vg[-1L, , drop = FALSE]
    return(vg)
  }
  d <- np1 * vg$dist[1L] / rest
  g <- (np1 * vg$gamma[1L] - zero_ss) / rest
  if (!is.finite(d) || !is.finite(g) || d > first_bound * (1 + 1e-8) ||
      g < -1e-8 * abs(vg$gamma[1L])) return(vg)
  vg$np[1L]    <- rest
  vg$dist[1L]  <- d
  vg$gamma[1L] <- max(g, 0)
  vg
}


#' Fit a nugget-plus-structure variogram model by profile weighted least squares
#'
#' The fit behind \code{estimate_sac_range()}.  The criterion is
#' \pkg{gstat}'s default for \code{fit.variogram()}: squared error weighted
#' by \eqn{N / h^2}.  The lags are divided by the cutoff and the
#' semivariances by their pair-weighted mean first, so nothing in the fit
#' depends on a unit.  For a trial range the nugget and the partial sill are
#' the two-parameter least-squares solution with both kept non-negative,
#' which is exact (the unconstrained solution, or one of the two edges).
#' The only search is therefore over the range: a logarithmic grid of
#' effective ranges from a quarter of the first lag to \code{upper} times
#' the cutoff, then a one-dimensional refinement around the best grid point.
#' There is no start value, no iteration limit and no singular fit, so the
#' same variogram gives the same model whatever the units, and a change in
#' the last bits of the data moves the result by about as much.
#'
#' The nugget can be given a floor above zero (\code{nugget_min}): the
#' semivariance of pairs at zero separation, where a layer has repeat visits
#' to its locations.  That is the variance between visits, and the limit of
#' the semivariance from positive lags cannot lie under it.  The solution
#' stays exact: the unconstrained one where it clears the floor, else the
#' better of the fit with the nugget on its floor and the constant (the
#' weighted mean semivariance, or the floor when that is higher).
#'
#' @param vg An empirical variogram (\code{np}, \code{dist}, \code{gamma}).
#' @param cutoff_dist The largest lag, in the units of \code{vg$dist}.
#' @param family \code{"Exp"} or \code{"Sph"}.
#' @param upper Top of the range grid, as a multiple of the cutoff.
#' @param nugget_min The least nugget allowed, in the units of
#'   \code{vg$gamma}; 0 for none.
#' @return \code{NULL} when fewer than three bins hold pairs.  Otherwise a
#'   \pkg{gstat} variogram model.  When no model with a range fits better
#'   than a constant it is a pure nugget (one row); else a nugget row and a
#'   structured row.  Attributes: \code{converged} (\code{FALSE} when the
#'   criterion is still falling at either end of the grid, so that the range,
#'   the nugget and the sill are the grid's edge and not a minimum),
#'   \code{SSErr} (the weighted squared error on the scaled axes) and
#'   \code{gain} (the share of the constant model's weighted squared error
#'   that the model removes).
#' @keywords internal
#' @noRd
.sac_fit_model <- function(vg, cutoff_dist, family = "Exp", upper = 4,
                           n_grid = 121L, nugget_min = 0) {
  if (!is.data.frame(vg) || !all(c("np", "dist", "gamma") %in% names(vg)))
    return(NULL)
  ok <- is.finite(vg$gamma) & is.finite(vg$dist) & is.finite(vg$np) &
    vg$np > 0 & vg$dist > 0
  if (sum(ok) < 3L || !is.finite(cutoff_dist) || cutoff_dist <= 0) return(NULL)
  h  <- vg$dist[ok] / cutoff_dist
  np <- as.numeric(vg$np[ok])
  gbar <- sum(np * vg$gamma[ok]) / sum(np)
  if (!is.finite(gbar) || gbar <= 0) return(NULL)
  g <- vg$gamma[ok] / gbar
  w <- np / h^2
  w <- w / sum(w)
  shape <- if (identical(family, "Sph"))
    function(r) ifelse(h < r, 1.5 * h / r - 0.5 * (h / r)^3, 1)
  else function(r) 1 - exp(-3 * h / r)          # r is the EFFECTIVE range
  mg   <- sum(w * g)
  # The nugget's floor on the scaled axis: 0 unless one was given.
  n0 <- if (is.finite(nugget_min) && nugget_min > 0) nugget_min / gbar else 0
  # The best constant: the weighted mean, or the floor when that is higher.
  c0   <- max(mg, n0)
  sse0 <- sum(w * (g - c0)^2)
  # Nugget, partial sill and weighted squared error at effective range r.
  # Two terms, each with a lower bound: the least-squares solution is the
  # unconstrained one where it is allowed, and otherwise lies on an edge,
  # the constant (no partial sill) or the nugget on its floor.
  solve_at <- function(r) {
    f   <- shape(r)
    mf  <- sum(w * f)
    df  <- f - mf
    sff <- sum(w * df^2)
    best <- c(c0, 0, sse0)
    if (sff > 1e-12) {
      b1 <- sum(w * df * (g - mg)) / sff
      b0 <- mg - b1 * mf
      if (b1 > 0 && b0 >= n0) {
        s <- sum(w * (g - b0 - b1 * f)^2)
        if (s < best[3L]) best <- c(b0, b1, s)
      } else if (b1 > 0 || n0 > mg) {
        # The nugget on its floor, and the partial sill solved alone.
        c1 <- sum(w * f * (g - n0)) / sum(w * f^2)
        s  <- sum(w * (g - n0 - c1 * f)^2)
        if (c1 > 0 && s < best[3L]) best <- c(n0, c1, s)
      }
    }
    best
  }
  lr  <- seq(log(min(h) / 4), log(upper), length.out = n_grid)
  val <- vapply(lr, function(x) solve_at(exp(x))[3L], numeric(1))
  j   <- which.min(val)
  r_hat <- exp(lr[j])
  at_edge <- j == 1L || j == n_grid
  if (!at_edge) {
    o <- stats::optimize(function(x) solve_at(exp(x))[3L],
                         c(lr[j - 1L], lr[j + 1L]), tol = 1e-10)
    if (isTRUE(o$objective <= val[j])) r_hat <- exp(o$minimum)
  }
  b <- solve_at(r_hat)
  if (!(b[2L] > 0) || !(b[3L] < sse0 * (1 - 1e-8))) {
    m <- gstat::vgm(psill = if (n0 > mg) nugget_min else mg * gbar,
                    model = "Nug", range = 0)
    attr(m, "converged") <- TRUE
    attr(m, "SSErr") <- sse0
    attr(m, "gain") <- 0
    return(m)
  }
  # A nugget on its floor is the floor, to the last bit (dividing by the
  # scale and multiplying back does not promise that).
  nug <- if (n0 > 0 && b[1L] <= n0) nugget_min else b[1L] * gbar
  m <- gstat::vgm(psill = b[2L] * gbar, model = family,
                  range = r_hat * cutoff_dist / (if (identical(family, "Sph")) 1 else 3),
                  nugget = nug)
  attr(m, "singular")  <- FALSE
  attr(m, "converged") <- !at_edge
  attr(m, "SSErr")     <- b[3L]
  attr(m, "gain")      <- 1 - b[3L] / sse0
  m
}


#' Add a second, longer structure of a given range to a fitted variogram model
#'
#' One exponential structure fitted to lags that resolve the spacing of the
#' points follows the shortest structure a layer has.  When the layer also
#' has a weaker structure many times longer, the range of that fit is the
#' short one, and blocks sized from it leave the long structure across
#' folds.  Equal lags, blind to the short structure, fit the long one.  This
#' puts the two together on one variogram: a nugget, a short exponential
#' structure and a long one, both ranges given, with the nugget and the two
#' partial sills solved by non-negative least squares (exactly, by trying
#' each set of free terms) under the weights of \code{.sac_fit_model()}.
#' The short range is not searched again.  Left free beside a long
#' structure it shrinks towards the first lag and takes the nugget with it:
#' on a layer with a nugget of 0.4 that put the nugget at 0 in 2 of 5 fits.
#'
#' @param vg An empirical variogram (\code{np}, \code{dist}, \code{gamma}).
#' @param cutoff_dist The largest lag.
#' @param r_short,r_long Effective ranges of the two structures.
#' @param min_share The least share of the model's semivariance at the
#'   largest lag that the long structure has to carry.
#' @param converged Whether the fit that gave \code{r_long} found a minimum
#'   inside the ranges it searched.  When it did not, \code{r_long} is the
#'   edge of that search, and the long structure's partial sill, solved at
#'   that range, is no more a fitted value than the range is.
#' @param nugget_min The least nugget allowed, as in \code{.sac_fit_model()}.
#' @return \code{NULL} when the variogram does not support the long
#'   structure: its share is under \code{min_share}, or the short structure
#'   gets no partial sill.  Otherwise a \pkg{gstat} model with a nugget and
#'   two exponential rows, short then long, with attributes \code{converged}
#'   (as given: the model's sill is a fitted value only when the long range
#'   is), \code{SSErr} and \code{long_share}.
#' @keywords internal
#' @noRd
.sac_fit_nested <- function(vg, cutoff_dist, r_short, r_long, min_share = 0.10,
                            converged = TRUE, nugget_min = 0) {
  if (!is.data.frame(vg) || !all(c("np", "dist", "gamma") %in% names(vg)))
    return(NULL)
  ok <- is.finite(vg$gamma) & is.finite(vg$dist) & is.finite(vg$np) &
    vg$np > 0 & vg$dist > 0
  if (sum(ok) < 5L || !is.finite(r_short) || !is.finite(r_long) ||
      r_short <= 0 || r_long <= r_short) return(NULL)
  h  <- vg$dist[ok] / cutoff_dist
  np <- as.numeric(vg$np[ok])
  gbar <- sum(np * vg$gamma[ok]) / sum(np)
  if (!is.finite(gbar) || gbar <= 0) return(NULL)
  g  <- vg$gamma[ok] / gbar
  w  <- np / h^2
  w  <- w / sum(w)
  rs <- r_short / cutoff_dist
  rl <- r_long / cutoff_dist
  # The nugget's floor is taken off first and put back after, so the three
  # terms solved are all at or above zero.
  n0 <- if (is.finite(nugget_min) && nugget_min > 0) nugget_min / gbar else 0
  g  <- g - n0
  X <- cbind(1, 1 - exp(-3 * h / rs), 1 - exp(-3 * h / rl))
  M <- crossprod(X * w, X)
  v <- as.numeric(crossprod(X * w, g))
  gwg <- sum(w * g * g)
  best <- c(0, 0, 0, gwg)
  for (S in list(1L, 2L, 3L, c(1L, 2L), c(1L, 3L), c(2L, 3L), 1:3)) {
    b <- tryCatch(solve(M[S, S, drop = FALSE], v[S]), error = function(e) NULL)
    if (is.null(b) || any(!is.finite(b)) || any(b < 0)) next
    sse <- gwg - 2 * sum(b * v[S]) + sum(b * as.numeric(M[S, S, drop = FALSE] %*% b))
    if (sse < best[4L]) { full <- numeric(3); full[S] <- b; best <- c(full, sse) }
  }
  on_floor <- n0 > 0 && !(best[1L] > 0)
  best[1L] <- best[1L] + n0
  # The long structure's part of the model at the largest lag.
  long_at  <- best[3L] * (1 - exp(-3 / rl))
  total_at <- best[1L] + best[2L] * (1 - exp(-3 / rs)) + long_at
  share <- if (total_at > 0) long_at / total_at else 0
  if (!(best[2L] > 0) || !(best[3L] > 0) || !(share >= min_share)) return(NULL)
  m <- gstat::vgm(psill = best[3L] * gbar, model = "Exp", range = r_long / 3,
                  add.to = gstat::vgm(psill = best[2L] * gbar, model = "Exp",
                                      range = r_short / 3,
                                      nugget = if (on_floor) nugget_min
                                               else best[1L] * gbar))
  attr(m, "singular")   <- FALSE
  attr(m, "converged")  <- isTRUE(converged)
  attr(m, "SSErr")      <- best[4L]
  attr(m, "long_share") <- share
  m
}


#' Test for spatial structure at the nearest-neighbour scale
#'
#' Moran's I of the modelled variable on row-standardised k-nearest-neighbour
#' weights, for k = 2, 4, 8 and 16, each standardised with the Cliff and Ord
#' moments for regression residuals (\code{.morans_residual_moments()}), so
#' the same statistic serves a raw response (an intercept-only design) and
#' the residuals of a trend on predictors (that trend's design).  The p-value
#' is the smallest one-sided normal p over the four graphs, multiplied by
#' their number (Bonferroni).  It needs no random numbers.  Several graphs
#' are used because one neighbourhood size misses structure at other scales.
#' Measured on a 1000 m square: on white noise the adjusted p-value is below
#' 0.05 in 0.041 of 1000 layers at n = 300 and 0.053 at n = 100; with a
#' nugget of a fifth of the sill it is below 0.05 for 88 percent of fields
#' of effective range 60 at n = 300, 75 percent of range 150 at n = 100, and
#' 21 percent of range 30 at n = 300, where the range is about the spacing
#' of the points and there is little for any test to find.
#'
#' Rows that share a location (repeat visits) are not neighbours in space:
#' they are the same place.  Left as rows, each visit's nearest neighbours
#' are the other visits to its own location, and the statistic measures how
#' alike those are, which a location's own effect produces with no structure
#' between locations at all: on 100 locations visited 5 times, with
#' independent location effects, it found "structure" in 40 of 40 layers.
#' So the test is made on locations: the response and each column of the
#' design are averaged over the rows of a location, and the neighbours are
#' the nearest other locations.  The moments are those of independent,
#' equal-variance errors; means over unequal numbers of visits differ in
#' variance, which they do not allow for, so with very unequal visit counts
#' the p-value is approximate.
#'
#' @param xy Coordinates (n by 2).
#' @param y The response.
#' @param X Design matrix of the trend removed from it, or \code{NULL} for an
#'   intercept alone.  The statistic is computed on the OLS residuals of
#'   \code{y} on \code{X}.
#' @return \code{NULL} when it cannot be computed, else a list with \code{p}
#'   (the adjusted p-value), \code{morans_i} (Moran's I on the 8 nearest
#'   neighbours, or on the graph nearest to 8 when the layer is too small),
#'   \code{table} (\code{k}, \code{morans_i}, \code{z}, \code{p} for
#'   each graph) and \code{n_locations} (the distinct locations the test was
#'   made on; fewer than the rows when locations repeat).
#' @keywords internal
#' @noRd
.sac_structure_test <- function(xy, y, X = NULL, ks = c(2L, 4L, 8L, 16L)) {
  n <- length(y)
  if (is.null(X)) X <- matrix(1, n, 1L)
  if (!is.matrix(X) || nrow(X) != n || nrow(xy) != n) return(NULL)
  y <- as.numeric(y)
  loc <- .sac_locations(xy)
  if (length(loc$first) < n) {
    xy <- xy[loc$first, , drop = FALSE]
    y  <- as.numeric(rowsum(y, loc$id)) / loc$size
    X  <- rowsum(X, loc$id) / loc$size
    n  <- length(y)
  }
  e <- tryCatch(qr.resid(qr(X), y), error = function(e) NULL)
  if (is.null(e)) return(NULL)
  den <- sum(e^2)
  if (!is.finite(den) || den <= 0) return(NULL)
  ks <- ks[ks <= n - 2L]
  if (!length(ks)) return(NULL)
  # One kd-tree query serves every graph when no distance ties with the next
  # at a graph's edge: the k nearest then get 1 / k each, which is what
  # .build_knn_weights() returns.  Ties (a regular grid) go to
  # .build_knn_weights(), whose tie rule does not depend on the order of the
  # rows.  The locations are distinct by now.
  nn <- NULL
  if (requireNamespace("FNN", quietly = TRUE) &&
      requireNamespace("Matrix", quietly = TRUE)) {
    kq <- min(max(ks) + 1L, n - 1L)
    nn <- FNN::get.knn(xy, k = kq)
    nn$tol <- 1e-9 * max(1, max(abs(xy)))
  }
  weights_for <- function(k) {
    if (!is.null(nn) && k < ncol(nn$nn.dist) &&
        !any(nn$nn.dist[, k + 1L] - nn$nn.dist[, k] <=
               nn$nn.dist[, k] * 1e-9 + nn$tol))
      return(Matrix::sparseMatrix(i = rep(seq_len(n), k),
                                  j = as.integer(nn$nn.index[, seq_len(k)]),
                                  x = 1 / k, dims = c(n, n), repr = "C"))
    .build_knn_weights(xy, k = k)
  }
  rows <- lapply(ks, function(k) {
    W  <- weights_for(k)
    sp <- inherits(W, "Matrix")
    S0 <- sum(W)
    I  <- (n / S0) * sum(e * as.numeric(W %*% e)) / den
    mom <- .morans_residual_moments(W, X, S0, sp)
    if (is.null(mom) || !isTRUE(mom$VI > 0)) return(c(k, I, NA_real_, NA_real_))
    z <- (I - mom$EI) / sqrt(mom$VI)
    c(k, I, z, stats::pnorm(z, lower.tail = FALSE))
  })
  tab <- as.data.frame(do.call(rbind, rows))
  names(tab) <- c("k", "morans_i", "z", "p")
  tab$k <- as.integer(tab$k)
  okp <- is.finite(tab$p)
  list(p = if (any(okp)) min(1, sum(okp) * min(tab$p[okp])) else NA_real_,
       morans_i = tab$morans_i[which.min(abs(tab$k - 8L))],
       table = tab, n_locations = as.integer(n))
}


#' Group the rows of a coordinate matrix by location
#'
#' Rows are at one location when both coordinates are equal to the last bit.
#'
#' @param xy Coordinates (n by 2).
#' @return A list with \code{id} (the location of each row, numbered in order
#'   of first appearance), \code{first} (the first row of each location) and
#'   \code{size} (the rows at each).
#' @keywords internal
#' @noRd
.sac_locations <- function(xy) {
  key <- paste(sprintf("%.17g", xy[, 1L] + 0), sprintf("%.17g", xy[, 2L] + 0), sep = "\r")
  id  <- match(key, unique(key))
  list(id = id, first = which(!duplicated(id)), size = tabulate(id, max(id, 0L)))
}


#' Semivariance at the two shortest separations a layer has, with no model
#'
#' Two numbers that need no variogram model and bracket the nugget, which is
#' the limit of the semivariance as the separation shrinks to zero.
#'
#' \strong{Nearest-neighbour pairs.}  Each location is paired with the
#' nearest other location, each such pair of locations counted once.  Half
#' the mean squared difference between a row of one and a row of the other,
#' averaged over the pairs of locations, is the semivariance at about the
#' nearest-neighbour distance.  The semivariance does not fall as the lag
#' shrinks, so this bounds the nugget from above, up to its sampling error.
#' One pair a location, and not the closest pairs of the layer, because
#' those come from its densest part alone.  Where a location has several
#' rows, every pairing of its rows with those of its neighbour counts, with
#' equal weight, and the pair of locations still counts once.
#'
#' \strong{Pairs at zero separation.}  Rows that share a location (repeat
#' visits).  Half the mean squared difference over every pair of rows at one
#' location is the variance between visits: the part of the nugget that is
#' measurement, without the variation at scales shorter than the spacing of
#' the locations, and so a bound on the nugget from below.  Every pair is
#' used, not one partner a row: the rows arrive sorted by their value, and
#' pairing each with the next pairs neighbouring order statistics, which
#' reads low (0.205 for a variance of 0.300 at 8 visits).
#'
#' @param xy Coordinates (n by 2).
#' @param v The modelled variable.
#' @return A list with \code{gamma}, \code{pairs} and \code{dist} (the
#'   semivariance of the nearest-neighbour pairs of locations, their number
#'   and their median separation; \code{NA}, 0 and \code{NA} when there is
#'   one location), and \code{zero_gamma}, \code{zero_pairs},
#'   \code{zero_ss} and \code{zero_df} (the semivariance of the pairs of
#'   rows at one location, their number, the sum of their half squared
#'   differences, and the degrees of freedom behind it: the rows less the
#'   locations; \code{NA}, 0, 0 and 0 when no location repeats).
#' @keywords internal
#' @noRd
.sac_short_lag <- function(xy, v) {
  n <- nrow(xy)
  out <- list(gamma = NA_real_, pairs = 0L, dist = NA_real_,
              zero_gamma = NA_real_, zero_pairs = 0L, zero_ss = 0, zero_df = 0L)
  if (n < 2L) return(out)
  v   <- as.numeric(v)
  loc <- .sac_locations(xy)
  G   <- length(loc$first)
  m   <- loc$size
  mu  <- v[loc$first]
  ss  <- numeric(G)                    # sum of squares about the location mean
  if (G < n) {
    mu <- as.numeric(rowsum(v, loc$id)) / m
    ss <- as.numeric(rowsum((v - mu[loc$id])^2, loc$id))
    # Over the m (m - 1) / 2 pairs of rows at a location, the half squared
    # differences sum to m / 2 times its sum of squares about its mean.
    zp <- sum(m * (m - 1) / 2)
    zs <- sum(m / 2 * ss)
    out$zero_pairs <- if (zp <= .Machine$integer.max) as.integer(zp) else zp
    out$zero_ss    <- zs
    out$zero_gamma <- zs / zp
    out$zero_df    <- n - G
  }
  if (G < 2L) return(out)
  nn <- .sac_knn(xy[loc$first, , drop = FALSE], 1L)
  i <- seq_len(G); j <- nn$index[, 1L]
  a <- pmin(i, j); b <- pmax(i, j)
  keep <- !duplicated(cbind(a, b))
  a <- a[keep]; b <- b[keep]
  # The mean squared difference between a row of location a and a row of
  # location b: the squared difference of the two means plus each
  # location's mean square about its own.
  out$gamma <- mean(0.5 * ((mu[a] - mu[b])^2 + ss[a] / m[a] + ss[b] / m[b]))
  out$pairs <- length(a)
  out$dist  <- stats::median(nn$dist[keep, 1L])
  out
}


#' The nugget of an estimated autocorrelation range
#'
#' The nugget variance of the variogram model behind a
#' \code{\link{estimate_sac_range}()} result: the model's semivariance as
#' the separation shrinks to zero, i.e. measurement error plus variation at
#' scales shorter than the first lag of the empirical variogram.  It is
#' extrapolated to zero from that lag, not observed, so it is only as good
#' as the first lag is short: \code{attr(x, "first_lag")} says where the
#' extrapolation starts.
#' With \code{lags = "adaptive"} (the default) that is the separation of the
#' closest pairs of points; with \code{lags = "equal"} it is the middle of a
#' bin \code{cutoff * max_dist / n_lags} wide, which on a dense layer can be
#' many times the spacing of close pairs, and structure shorter than that
#' bin is then counted as nugget.
#'
#' Three more attributes say how far to trust it.
#' \code{attr(x, "nugget_at_bound")} is \code{TRUE} when the fit put the
#' nugget at its lower bound (0, or within 1e-4 of the sill): that is the
#' bound and not an estimate, and a criterion that uses it as a noise
#' variance has no penalty left.  \code{attr(x, "short_lag_semivariance")}
#' is the semivariance of nearest-neighbour pairs, computed with no model:
#' each location is paired with the nearest other location
#' (\code{short_lag_pairs} and \code{short_lag_dist} give the number of
#' such pairs and their median separation).  A semivariance does not fall
#' as the lag shrinks, so this value bounds the nugget from above, up to its
#' sampling error; a fitted nugget well above it is contradicted by the
#' closest pairs.  \code{attr(x, "zero_lag_semivariance")} exists where
#' locations were visited more than once (\code{zero_lag_pairs} is not 0):
#' it is the semivariance of every pair of rows that share a location, that
#' is, the variance between visits to one place.  That is the measurement
#' part of the nugget, observed and not extrapolated, without the variation
#' at scales below the spacing of the locations, so it bounds the nugget
#' from below.  Such pairs are at no lag and are in no bin of the
#' variogram: the model, and so the nugget, describe positive separations.
#' The variogram fit uses the bound, though: where it rests on 30 degrees
#' of freedom or more (rows less locations) the nugget is not let under it,
#' and \code{attr(x, "nugget_floor")} holds the value (0 when there is no
#' floor).  With \code{detrend = "reml"} the model is built from the REML
#' parameters and takes no floor: its nugget can come out below the
#' variance between visits, and \code{nugget_floor} is 0.
#' A nugget equal to its floor is the variance between visits, a
#' measurement, and \code{nugget_at_bound} stays \code{FALSE} for it.
#'
#' When \code{attr(x, "lag_fallback")} is \code{TRUE} the result is the
#' equal-lag one although adaptive lags were asked for (see "The lags" in
#' \code{\link{estimate_sac_range}}), and its nugget is extrapolated from
#' the first equal bin.
#'
#' The nugget is carried as the \code{nugget} attribute of every classed
#' result, identified or rejected, because it is the number a resolution
#' criterion for a tessellation needs (the short-lag variance that no cell
#' can average away).  On a result refused for want of spatial structure it
#' is the variance of the modelled variable, which is the nugget of the
#' model that refusal supports.
#'
#' @param x A \code{sac_range} object, or anything else.
#' @return A single number: the nugget in the units of the response's
#'   variance; \code{NA_real_} when \code{x} carries no fitted model (a bare
#'   \code{NA} from a run that could not fit anything, a rejected result with
#'   too few lag bins to fit, or an object that is not a \code{sac_range}).
#' @seealso \code{\link{estimate_sac_range}}, which produces the object.
#' @family cross-validation
#' @examples
#' if (requireNamespace("gstat", quietly = TRUE)) {
#'   library(sf)
#'   # A field with a real nugget: half a unit of white noise on a unit sill.
#'   set.seed(3)
#'   n <- 250
#'   xy <- data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000))
#'   D  <- as.matrix(dist(xy))
#'   xy$z <- as.numeric(t(chol(exp(-D / 150) + diag(0.5, n))) %*% rnorm(n))
#'   r <- estimate_sac_range(st_as_sf(xy, coords = c("x", "y"), crs = 32632), "z")
#'   print(sac_nugget(r))  # the fitted nugget variance, on the sill's scale
#'   attr(r, "nugget_at_bound")         # FALSE: an estimate, not the bound
#'   attr(r, "short_lag_semivariance")  # what the closest pairs say, no model
#'   attr(r, "zero_lag_pairs")          # 0: no location was visited twice
#'   sac_nugget(NA)        # nothing fitted: NA
#' }
#' @export
sac_nugget <- function(x) {
  v <- attr(x, "nugget")
  if (is.null(v)) return(NA_real_)
  v <- suppressWarnings(as.numeric(v))
  if (length(v) != 1L || !is.finite(v)) NA_real_ else v
}


#' Estimate the spatial autocorrelation range from data
#'
#' Fits an exponential variogram model with a nugget and returns the
#' \emph{effective range}: three times the fitted range parameter, which is
#' where the semivariance reaches about 95 percent of the sill.  It is the
#' distance beyond which two observations are (near) uncorrelated, which is
#' what a block or a buffer has to exceed.
#'
#' The result is a function of the data alone.  The same layer gives the
#' same answer whatever the unit of the response or of the coordinates and
#' whatever the order of its rows, and a change in the last bits of the
#' values changes the answer by about as much.  Three things make it so.
#' The rows are sorted by their coordinates as supplied, then by the
#' response and the predictors, before anything else, so the projection
#' chosen for lon/lat input, the \code{n_max} subsample and the pair sums
#' do not depend on their order.  (Rows that share a location are therefore
#' in the order of their values, and no step treats the first of them as
#' special.)  The variogram is computed on the standardised variable, and a
#' response is called constant only when its standard deviation is under
#' 1e-12 of its largest absolute value, which no change of unit alters.
#' And the model is fitted without iteration: the criterion is \pkg{gstat}'s weighted least
#' squares (weights \eqn{N / h^2}) on lags divided by the cutoff; for each
#' trial range the nugget and the partial sill are solved exactly, with both
#' kept at or above zero; and the range is found by a search over a
#' logarithmic grid from a quarter of the first lag to four cutoffs, refined
#' in one dimension.  There is no start value, no iteration limit and no
#' fit that fails to converge.
#'
#' Only the exponential family is fitted.  On fields smoother than
#' exponential that makes the range long, which is the safe side (blocks too
#' large, cross-validation pessimistic).  Choosing the family by the smaller
#' weighted sum of squares would correct a spherical field but send
#' exponential fields low, the direction that leaks.
#'
#' \strong{Spatial structure is tested first.}  A variogram model can be
#' fitted to anything, and on a layer with no structure the fit returns a
#' finite range for a large share of draws.  So Moran's I of the modelled
#' variable is computed on its 2, 4, 8 and 16 nearest neighbours, each
#' standardised with the moments that fit the variable (those of regression
#' residuals when \code{predictor_vars} are given), and the smallest
#' one-sided p-value is multiplied by four.  When that adjusted p-value is
#' not below \code{structure_alpha}, the result is \code{NA} with
#' \code{rejected_reason} \code{"no spatial structure detected at the
#' nearest-neighbour scale"}, a pure-nugget model and no fitted range.  The
#' test uses no random numbers.  \code{structure_p}, \code{morans_i} and
#' \code{structure_test} are attached to every classed result.  Measured on
#' a 1000 m square with 300 points, the test refuses 388 of 400 white-noise
#' layers and 8 get a finite range, against 196 with the test off.
#' Its power is limited where few neighbours lie inside the range.  With a
#' nugget of a fifth of the sill it finds structure in 88 percent of fields
#' of effective range 60 on 300 points, and in 21 percent at range 30, which
#' is about the spacing of those points.  A short range on a sparse layer
#' can therefore be refused; \code{structure_alpha = 1} switches the
#' refusal off.  Rows that share a location (repeat visits) are one place
#' and not each other's neighbours, so on such a layer the test is made on
#' the means of the locations: what it asks is whether locations near each
#' other are alike, which a location's own effect, however large, does not
#' answer.  With very unequal numbers of visits the p-value is approximate.
#'
#' \strong{The lags.}  With \code{lags = "adaptive"} (the default) the first
#' bin runs from zero to the distance within which the closest
#' \code{max(50, n / 4)} pairs of distinct locations lie, and up to
#' \code{n_lags} further bins are spaced geometrically from there to the
#' cutoff; a bin holding fewer than 30 pairs is merged with the next.  The
#' short lags are then resolved at the spacing of the points, which is what
#' the nugget and a short range are read from.  With \code{lags = "equal"}
#' there are \code{n_lags} bins of equal width out to the cutoff, as
#' \pkg{gstat} bins by default; each is \code{cutoff * max_dist / n_lags}
#' wide, a structure shorter than the first of them is invisible and is
#' counted as nugget, and the range fitted is that of whatever longer
#' structure the layer has.  Below about a hundred points neither binning
#' resolves much and the adaptive bins bring no gain.
#'
#' Two rows at one location are a pair at zero separation, and are in no
#' bin, whichever the binning: the variogram and the model describe positive
#' separations.  Their semivariance, the variance between visits to a
#' location, is returned on its own as \code{zero_lag_semivariance}, and
#' the nugget fitted to the variogram is not let under it when the layer has
#' 30 or more rows beyond one a location to estimate it from (a model from
#' \code{detrend = "reml"} takes no such floor); see
#' \code{\link{sac_nugget}}.
#'
#' The adaptive lags put little weight on the long ones, and a weak, long
#' structure under much noise can be clear at equal lags while the adaptive
#' fit runs past the largest lag.  So when the adaptive range is refused as
#' past the fitted lags and the equal-lag fit has a range it would itself
#' return, the result is the equal-lag one, whole, exactly as
#' \code{lags = "equal"} gives it, with \code{lags} saying \code{"equal"}
#' and \code{lag_fallback} set to \code{TRUE}.  Its nugget is then the
#' equal-lag nugget.  No other refusal falls back.
#'
#' \strong{Structure at two scales.}  One structure fitted to lags that
#' resolve the spacing of the points follows the shortest structure of the
#' layer.  A layer can have a weaker structure several times longer as well,
#' and blocks sized from the short one would leave the long one across
#' folds.  So with \code{lags = "adaptive"} the equal-lag fit is made too,
#' and when its range is at least twice the adaptive one, a model with a
#' nugget and both structures, at those two ranges, is fitted to the
#' adaptive variogram.  That model is returned (\code{structures} is 2 and
#' \code{short_range} holds the range of the short structure) and the range
#' is the long one, unless the long structure carries under 3 percent of
#' the model's semivariance at the largest lag: then the lags that resolve
#' the layer give the equal-lag range nothing, and the single structure
#' stands.  The bar is low because the two mistakes are not alike.  A long
#' structure that is missed sizes blocks too small, which leaks.  One that
#' is believed in error returns the equal-lag range, which is too long and
#' makes a validation pessimistic, and \code{short_range} still holds the
#' short one.  When the equal-lag range is one that \code{lags = "equal"}
#' would itself refuse (past the fitted lags, as under an unremoved trend),
#' the model with it is refused as well, the answer is taken away and not
#' lengthened, and the long structure has to carry 10 percent.
#'
#' Measured on a clustered layer of 1200 points over a 1000 m square, 30
#' draws each.  With structures of range 20 and 200 both are kept in 20
#' draws, and the range returned is under half the equal-lag one in 1; with
#' 40 and 400, in 18 and 1; with 30 and 300, in 21 and 0; with 15 and 450,
#' in 28 and 1.  With one structure of range 40 a second is kept in 4, and
#' the range is within a factor 2 of the truth in 27 draws (19 with equal
#' lags); with one of range 30, in 7 and 25 (14).  The nugget comes from the
#' adaptive lags either way: where it is 0.3, the median is 0.32 to 0.35,
#' and 0.49 to 0.65 with equal lags.  A long structure whose equal-lag
#' range is under twice the adaptive one is not separated: the range is
#' then that of the single fit, which lies between the two.
#'
#' The estimate is the \strong{omnidirectional} (all-pairs) fit.  Directional
#' variograms are fitted as well, at 0° (N–S), 45°, 90° (E–W) and 135°
#' azimuths with a ±22.5° tolerance.  Those four windows tile all 180 distinct
#' azimuths exactly once, and their ranges are returned in the
#' \code{directional} attribute, with their largest-over-smallest ratio in
#' \code{anisotropy}.  They are a diagnostic and never the answer, for two
#' reasons.  Each direction sees about a quarter of the point pairs, and the
#' maximum of four quarter-sample fits is biased upward.  And the windows
#' are fixed to the coordinate axes, so any answer built from them changes
#' when the layer is rotated, which a property of the field must not do.
#' The all-pairs fit is the best-powered estimate available and is invariant
#' to rotation.  When it identifies no range, no range is returned: the
#' largest directional range does not stand in for it.
#'
#' Where a field is \emph{known} to be anisotropic, blocks must be at least
#' as large as the longest autocorrelation range to avoid leakage, and the
#' conservative choice is to size them from the longest directional range
#' explicitly.  Read it from \code{directional_fitted}, not
#' \code{directional}: on a strongly anisotropic field the major axis is the
#' direction most likely to run past the fitted lags, which leaves it
#' \code{NA} in \code{directional}, so \code{max()} of that is \code{NA}, or
#' with \code{na.rm = TRUE} the second-longest range.  Check
#' \code{directional_status} first: a major axis marked \code{"over_cutoff"}
#' has no identified range at all, and a longer \code{cutoff} or
#' \code{\link{make_folds}(method = "nndm")} is the way on.  A ratio above
#' 1.5 is written to the package log at INFO level with that advice, which
#' reaches the session log file but not the console (the ratio passes 1.5 on
#' most isotropic fields too); \code{print()} shows the directional ranges
#' and the ratio, and \code{attr(range, "anisotropy")} holds it.
#'
#' A direction whose variogram is flat or has too few bins, whose range runs
#' beyond the longest fitted lag, or whose fit has no minimum inside the
#' ranges searched is excluded and recorded as \code{NA} in the
#' \code{directional} attribute.
#'
#' The model is fitted \strong{with a nugget}.  A nugget-free model forces
#' the curve through the origin, and on any real measurement (which has one)
#' the \eqn{N / h^2} weights buy that constraint by collapsing the range.
#' The fitted nugget can sit at its lower bound of zero; the result then
#' carries \code{nugget_at_bound = TRUE}, and \code{\link{sac_nugget}} says
#' what else is attached to judge a nugget by.
#'
#' The returned range is in the coordinate units of the (projected) data and
#' can be passed directly to \code{make_folds(block_size = ...)} so that CV
#' blocks are at least as wide as the autocorrelation range.
#'
#' @param points_sf An sf object with point geometries (will be projected
#'   automatically if in geographic CRS).  Non-POINT geometry is reduced to
#'   representative points; any Z or M dimension is dropped, because
#'   \code{gstat::variogram()} uses every coordinate dimension and an XYZ layer
#'   would otherwise return a range in 3-D while every consumer of it works in
#'   2-D map distance; and rows with empty or non-finite coordinates are dropped
#'   with a logged count.
#' @param response_var Character(1) name of the response column.  Rows where
#'   it is missing or infinite are left out before anything else is done, so
#'   they neither take part in the choice of projection for lon/lat input
#'   nor use up the \code{n_max} subsample.  For a
#'   count or other response whose variance tracks its mean, see the section
#'   on non-Gaussian responses: the range is still estimated, but it is a
#'   less reliable number than for a Gaussian response.
#' @param predictor_vars Optional character vector.  When supplied, the
#'   trend on these predictors is removed first and the variogram describes
#'   the residual autocorrelation, the part a spatial model has to handle
#'   once the covariates have done their work.  How the trend is removed is
#'   set by \code{detrend}, and it matters: see "Detrending and the
#'   residual-variogram bias".  Rows whose response or predictor is missing
#'   or infinite are left out of that fit, and so of the variogram, with a
#'   logged count.  When that leaves fewer than 30 rows the result is the
#'   unclassed \code{NA} that says how many were left (see Value).  One
#'   case differs, as it always has: when no row at all has every predictor
#'   (one is missing on every row) no trend can be fitted, and the variogram
#'   is of the raw response, with an R warning and \code{detrended} set to
#'   \code{FALSE}.
#' @param n_max Maximum number of points the variogram is computed on, at
#'   least 30; \code{Inf} for no limit.  Variogram estimation is O(n²) so
#'   this keeps runtime bounded.  A layer with more rows that carry a value
#'   is cut to a seeded random subsample of this size, drawn after rows with
#'   nothing to model have been left out and with the rows in a fixed order,
#'   so the subsample does not depend on how the layer is sorted.  The
#'   counts are returned as the attributes \code{n_used}, \code{n_missing}
#'   and \code{n_subsampled}.
#' @param cutoff Fraction of the maximum inter-point distance to use as the
#'   variogram lag cutoff, above 0 and at most 1.  That distance is the
#'   farthest pair, found on the convex hull, and not the bounding-box
#'   diagonal, which depends on how the axes are oriented.  Default 0.5.
#' @param range_frac Positive numeric.  A fitted range exceeding
#'   \code{range_frac * cutoff * max_dist} (that is, beyond the longest lag
#'   the empirical variogram was actually fitted over) is treated as
#'   unidentified and \code{NA_real_} is returned.
#'   A model fit yields a finite number even when the
#'   variogram never reaches a sill, and such a value extrapolates past the
#'   observed lags instead of measuring a long autocorrelation range.  Passing
#'   it to \code{make_folds(auto_range = TRUE)} would collapse the block grid to
#'   a single block.  Default 1.0; raise it to accept ranges extrapolated
#'   beyond the fitted lags.  The bound does not guarantee room for two
#'   blocks: at the defaults it is half the farthest-pair distance, about
#'   0.71 of the side of a square layer, and a block grid needs a range below
#'   half the width of the bounding box in one direction or the other.  An
#'   accepted range between the two leaves
#'   \code{make_folds(auto_range = TRUE)} room for a single block of that
#'   size (see its \code{auto_range} argument for what it does then).
#'   Lowering \code{range_frac} to fit the
#'   grid would turn those estimates into \code{NA} and the blocks into
#'   geometric ones smaller than the range.
#' @param seed RNG seed for the \code{n_max} subsample, restored afterwards so
#'   the caller's random stream is untouched.  Default \code{123L}: the
#'   subsample is an internal approximation and no part of the answer, and
#'   leaving it unseeded would make the returned range differ between runs on
#'   identical input and silently advance the caller's RNG.  Pass \code{NULL}
#'   for an unseeded subsample, or a
#'   different number to check how sensitive the estimate is to the subsample.
#'   Ignored when no more than \code{n_max} rows carry a value, where nothing
#'   is sampled.  The \code{reml_max_n} subsample uses the same seed, and so
#'   does the draw of one row a location that \code{detrend = "reml"} makes
#'   on a layer with repeat visits.
#' @param detrend How the trend on \code{predictor_vars} is removed;
#'   ignored when there are none.  \code{"ols"} (default) fits it by ordinary
#'   least squares and fits the variogram to the residuals, which
#'   underestimates the range (see the section
#'   below).  \code{"reml"} fits the trend and an exponential-plus-nugget
#'   covariance together by residual maximum likelihood with
#'   \code{nlme::gls()}, returns the REML range, and attaches the empirical
#'   variogram of the REML residuals for inspection.  It needs \pkg{nlme},
#'   costs \eqn{O(n^3)}, and so runs on at most \code{reml_max_n} points;
#'   when it does not converge the \code{"ols"} path runs instead with an R
#'   warning saying so.  The start of the fit is chosen from the profile
#'   REML likelihood on a grid of ranges and nugget proportions, and when
#'   more than one start is fitted the fit with the largest log-likelihood
#'   is kept.
#' @param reml_max_n Positive integer, at least 30.  With \code{detrend =
#'   "reml"}, the trend and covariance are fitted on a seeded random subsample
#'   of this many points when the layer has more; the fitted trend is then
#'   removed from every point.  The fit takes one row a location: where a
#'   location has several rows (repeat visits), one of them is drawn at
#'   random under \code{seed}, so that the nugget fitted is the variance
#'   between visits and not that of their minimum or of their mean.
#'   Default 400.  Raise it for a better-determined fit at a cubic cost.
#' @param keep_directional_fits \code{TRUE} or \code{FALSE}.  Attach each
#'   direction's empirical
#'   variogram and fitted model as \code{directional_fits}?  Defaults to
#'   \code{FALSE}: the four variograms are most of the object's size, while
#'   the numbers
#'   read from them (\code{directional}, \code{directional_fitted},
#'   \code{directional_status}, \code{anisotropy}) are attached either
#'   way, and \code{plot()} draws the effective variogram from its own
#'   attribute.  Set \code{TRUE} to inspect the directional curves.
#' @param lags How the lags are binned: \code{"adaptive"} (default) or
#'   \code{"equal"}; see "The lags" above.
#' @param n_lags A whole number from 5 to 100, default 15.  With
#'   \code{lags = "equal"}, the number of bins.  With
#'   \code{lags = "adaptive"}, the most bins spaced geometrically past the
#'   first; fewer are used when that many would hold under 30 pairs each.
#' @param structure_alpha The level of the test for spatial structure, above
#'   0 and at most 1; default 0.05.  A layer whose adjusted p-value is not
#'   below it gets no range.  \code{1} switches the refusal off; the test
#'   is still computed and attached.
#' @section Detrending and the residual-variogram bias:
#' Fitting a variogram to the residuals of a least-squares trend
#' underestimates both the sill and the range, because the trend fit absorbs
#' part of the long-wavelength spatial variation (Lark, Cullis and Welham
#' 2006).  Blocks sized from that range are then too small and a blocked
#' validation is less conservative than it claims.  How large the effect is
#' depends on how smooth the trend terms are in space, and on how many there
#' are: for ordinary covariates it is a few percent, and for trend surfaces
#' in the coordinates it is large.  Iterating between a GLS trend
#' fit and a variogram refit (Neuman and Jacobson 1984) recovers only part of
#' it, because the variogram of GLS residuals
#' is biased too; REML does not fit a variogram to residuals at all, which is
#' why \code{detrend = "reml"} returns its own range estimate.  Its price is
#' a single family (exponential with nugget), a cubic cost in \code{n}, and
#' a somewhat wider sampling spread.  The default stays \code{"ols"} so that
#' existing scripts keep their method; a script that detrends on smooth
#' or coordinate-based terms should pass \code{detrend = "reml"}, or size its
#' blocks with a margin.  Detrending removes a trend the predictors explain.
#' It does not make an unidentified range identified.
#'
#' @section Count and other non-Gaussian responses:
#' The empirical variogram assumes second-order stationarity: a variance that
#' is the same everywhere, so that semivariance depends on separation alone.
#' A count response breaks that assumption by construction, because its
#' variance tracks its mean, and so does any response with a mean-variance
#' relationship (rates, proportions, skewed amounts).  On such data the
#' variogram mixes distance-dependence with mean-dependence, and the fitted
#' range is not the thing it claims to be: where the mean is high the
#' semivariance is inflated whatever the distance, which flattens the curve
#' and can lengthen or shorten the apparent range depending on where the
#' high-mean regions sit.
#'
#' Supplying \code{predictor_vars} to detrend helps, because it removes the
#' part of the mean the covariates explain, but it does not fix the variance
#' structure: the residuals of an OLS fit to a count still have a variance
#' that tracks the fitted mean.  The principled remedy is a variogram of
#' Pearson residuals from a model in the right family, which this function
#' does not compute.  Treat the range from a count response as no
#' more than an order of magnitude, size blocks conservatively from it, and
#' prefer \code{\link{make_folds}(method = "nndm")}, which does not depend on a
#' fitted range at all.
#'
#' @return A single number, of class \code{sac_range} in the first two of the
#'   three shapes below and an unclassed \code{NA} in the third; all three
#'   behave as an ordinary number.  The shapes carry different attributes:
#'   \describe{
#'     \item{Success}{A positive effective range in projected coordinate units,
#'       with the fit attached as attributes \code{directional} (the 0°, 45°,
#'       90° and 135° ranges, named by azimuth; \code{NA} where that
#'       direction's fit was unusable), \code{anisotropy} (largest
#'       over smallest), \code{anisotropy_used} (always \code{FALSE}: the
#'       all-pairs fit is the estimate; kept for code that reads it),
#'       \code{directional_status}
#'       (per azimuth, why a direction is
#'       \code{NA} in \code{directional}: \code{"ok"}, \code{"over_cutoff"}
#'       (its range ran past the largest lag fitted), \code{"not_converged"}
#'       (its fit has no minimum inside the ranges searched)
#'       or \code{"no_fit"} (a flat variogram, or too few bins)),
#'       \code{directional_fitted} (the range each
#'       direction's fit reported whether or not it was usable, so a refused
#'       directional range stays recoverable) and \code{directional_fits} (a
#'       list by azimuth of each direction's empirical \code{variogram} and
#'       fitted \code{model}, \code{NULL} where there is none, and
#'       \code{NULL} altogether unless \code{keep_directional_fits = TRUE}),
#'       \code{detrended} (logical: whether the variogram is of the
#'       residuals on \code{predictor_vars} or of the raw response.  A
#'       missing predictor is an error, and a failed detrending fit
#'       warns and falls back to the raw response with this set to
#'       \code{FALSE}), \code{detrend_method} (\code{"ols"} or \code{"reml"}
#'       when detrended, \code{NA} otherwise), \code{reml} (with
#'       \code{detrend = "reml"}: a list with \code{n_used},
#'       \code{subsampled}, \code{nugget_prop}, \code{sigma2},
#'       \code{logLik} (the REML log-likelihood of the fit kept),
#'       \code{logLik_gain} (its gain over the same trend with independent
#'       errors) and \code{n_fits} (the starts fitted); \code{NULL}
#'       otherwise),
#'       \code{crs} (the projected CRS the variogram was
#'       fitted in: the unit of the range), \code{max_dist},
#'       \code{cutoff_dist}, \code{variogram} (the empirical variogram, in
#'       the units of the response's variance),
#'       \code{variogram_model} (the fitted model, in \pkg{gstat}'s layout,
#'       or with \code{detrend = "reml"} one built from the REML
#'       parameters), \code{nugget} (that model's nugget variance; see
#'       \code{\link{sac_nugget}}), \code{nugget_at_bound} (\code{TRUE} when
#'       the nugget is at its lower bound, within 1e-4 of the sill, and so
#'       is not an estimate), \code{first_lag} (the mean separation in the
#'       first bin: the shortest lag resolved, from which the nugget is
#'       extrapolated), \code{bins_below_range} (how many bins lie inside
#'       the fitted range; with fewer than three the split between nugget
#'       and structure rests on very little),
#'       \code{short_lag_semivariance}, \code{short_lag_pairs} and
#'       \code{short_lag_dist} (the semivariance of nearest-neighbour pairs
#'       of locations with no model, their number and their median
#'       separation: an upper bound on the nugget),
#'       \code{zero_lag_semivariance} and \code{zero_lag_pairs} (the
#'       semivariance of the pairs of rows that share a location, which is
#'       the variance between visits to one place and a lower bound on the
#'       nugget, and their number; \code{NA} and 0 when no location has two
#'       rows), \code{nugget_floor} (the value the fitted nugget was not let
#'       under: \code{zero_lag_semivariance} when the rows outnumber the
#'       locations by 30 or more, else 0, and always 0 with
#'       \code{detrend = "reml"}), \code{structure_p},
#'       \code{morans_i} and
#'       \code{structure_test} (the adjusted p-value of the test for
#'       spatial structure, Moran's I on the 8 nearest neighbours, and the
#'       statistic, z and p for each neighbourhood size),
#'       \code{structures} and \code{short_range} (how many structures the
#'       model has besides the nugget, and the effective range of the short
#'       one when there are two, \code{NA} otherwise), \code{lags},
#'       \code{lag_boundaries} and \code{lag_fallback} (the binning behind
#'       the result, its upper bin edges, and \code{TRUE} when adaptive lags
#'       were asked for and the equal-lag result is returned; see "The
#'       lags"), and
#'       \code{n_used}, \code{n_missing} and \code{n_subsampled} (the points
#'       the variogram was computed on, the rows left out for a missing
#'       response or predictor, and the rows left out by the \code{n_max}
#'       subsample), so the fit can be inspected and need not
#'       be taken on trust.}
#'     \item{Rejected range}{\code{NA_real_} when no range is identified.
#'       It is classed \code{sac_range} as well, so it prints as \code{NA}
#'       without dumping its attributes, and it carries every attribute of
#'       a success (the evidence for the rejection), plus
#'       \code{rejected_range} (the value that was refused, \code{NA} when
#'       none was fitted) and \code{rejected_reason}, one of:
#'       \itemize{
#'         \item \code{"no spatial structure detected at the
#'           nearest-neighbour scale"}: the test above.  The model attached
#'           is a pure nugget at the variance of the modelled variable, and
#'           nothing fitted is handed on.
#'         \item \code{"the empirical variogram is flat: no spatial
#'           structure at the lags resolved"}: no model with a range fits
#'           the variogram better than a constant.  The model attached is
#'           that constant, a pure nugget.
#'         \item \code{"empirical variogram decreases with distance"}: the
#'           semivariance falls over the shorter lags (a net fall of more
#'           than 15 percent of the mean semivariance there, weighted by
#'           pairs), which is the shape of a periodic, hole-effect
#'           structure or of a variance that differs between a dense
#'           cluster and the rest of the layer.  Sampling noise in the
#'           short-lag bins of a small sample can make that fall too.
#'         \item \code{"fitted range exceeds the largest lag fitted"}: the
#'           range is past \code{range_frac * cutoff * max_dist} (see
#'           \code{range_frac}).  An unremoved trend ends here when the
#'           variogram rises past the fitted lags; a milder trend only
#'           lengthens the fitted range and passes, which is what
#'           \code{predictor_vars} is for.
#'         \item \code{"fitted range is below the shortest lag fitted"}:
#'           too few pairs of points lie inside the range to identify it,
#'           because it is shorter than the shortest lag the
#'           empirical variogram resolves (\code{first_lag}) or, with
#'           \code{detrend = "reml"}, than the distance within
#'           which 30 pairs of the points the REML fit used lie, when that is
#'           shorter (the REML range is fitted to the point pairs, not to the
#'           bins).  A structure that short cannot be told from a nugget.
#'           \code{range_floor} holds the distance the refused range fell
#'           short of.
#'         \item \code{"fitted range is non-positive or non-finite"}.
#'         \item \code{"no variogram model could be fitted (singular
#'           fits)"}: fewer than three lag bins hold point pairs.
#'           \code{variogram_model} is \code{NULL} and \code{nugget} is
#'           \code{NA}.
#'       }
#'       The reason names one ground for the refusal, in the order above.
#'       \code{attr(attr(range, "variogram_model"), "converged")} says one
#'       thing more about the model: it is \code{FALSE} when the
#'       least-squares criterion is still falling at an edge of the ranges
#'       searched (a variogram with no sill inside four cutoffs), and the
#'       nugget and sill of such a model are that edge and not fitted
#'       values.  For a model with two structures the flag is that of the
#'       equal-lag fit the long range came from: when that fit has no
#'       minimum, the long structure's range and partial sill are the edge
#'       of its search and the flag is \code{FALSE}, although the nugget
#'       and the short structure were fitted at the adaptive lags.  A model
#'       built from REML parameters carries no such
#'       attribute.  \code{\link{resolution_profile}()} reads it, and takes
#'       no nugget from a model marked \code{FALSE}.  Results
#'       stored by earlier versions can carry
#'       \code{"variogram model did not converge"} as their reason and
#'       \code{anisotropy_used = TRUE}; neither is produced any more.}
#'     \item{No fit}{An unclassed \code{NA_real_} when estimation could not
#'       be attempted at all, whose one attribute, \code{rejected_reason},
#'       says why: \code{"package 'gstat', which the variogram needs, is not
#'       installed"}, \code{"<n> points, fewer than the 30 a variogram range
#'       is estimated from"}, \code{"<n> point(s) with a finite value to
#'       model, fewer than the 30 a variogram range is estimated from"}
#'       (the count is of the whole layer's rows with a value),
#'       \code{"the response is constant"}, \code{"the residuals on
#'       predictor_vars are constant: the predictors explain the response
#'       exactly"} or \code{"the points have no extent (the largest distance
#'       between them is zero or could not be computed)"}.  Nothing is
#'       fitted, so none of the other attributes above exist: no
#'       \code{variogram}, \code{variogram_model} or \code{rejected_range},
#'       which is what tells it from a range that was fitted and refused.}
#'   }
#'   Attributes and the class do not affect \code{is.na()} or
#'   \code{is.finite()}, so every downstream guard treats all three the same
#'   way it always did.
#' @references
#' Cliff, A. D. and Ord, J. K. (1981). \emph{Spatial Processes: Models and
#' Applications}. Pion, London.
#'
#' Lark, R. M., Cullis, B. R. and Welham, S. J. (2006). On spatial prediction
#' of soil properties in the presence of a spatial trend: the empirical best
#' linear unbiased predictor (E-BLUP) with REML. \emph{European Journal of
#' Soil Science}, 57(6), 787--799. \doi{10.1111/j.1365-2389.2005.00768.x}
#'
#' Neuman, S. P. and Jacobson, E. A. (1984). Analysis of nonintrinsic spatial
#' variability by residual kriging with application to regional groundwater
#' levels. \emph{Mathematical Geology}, 16(5), 499--521.
#' \doi{10.1007/BF01886329}
#' @seealso \code{\link{sac_nugget}} for the nugget behind the estimate,
#'   \code{\link{plot.sac_range}} to see the variogram the estimate rests on.
#' @family cross-validation
#' @examples
#' if (requireNamespace("gstat", quietly = TRUE)) {
#'   library(sf)
#'   # A Gaussian random field with an exponential covariance: range
#'   # parameter 100, so the true effective range is 3 x 100 = 300, plus a
#'   # small nugget.
#'   set.seed(9)
#'   n <- 150
#'   xy <- data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000))
#'   D  <- as.matrix(dist(xy))
#'   xy$z <- as.numeric(t(chol(exp(-D / 100) + diag(0.1, n))) %*% rnorm(n))
#'   pts <- st_as_sf(xy, coords = c("x", "y"), crs = 32632)
#'   r <- estimate_sac_range(pts, response_var = "z")
#'   # print() throughout, because only the last value of a braced block is
#'   # shown on its own, and every one of these is worth reading.
#'   print(r)                              # the effective range, in metres
#'   print(attr(r, "directional"))         # the four directional ranges
#'   print(attr(r, "variogram_model"))     # the fitted model behind it
#'   print(attr(r, "structure_p"))         # the test for spatial structure
#'
#'   # A field whose range the data cannot pin down: the variogram never
#'   # reaches a sill within the lags fitted, so the answer is NA with the
#'   # refused value attached rather than a long range asserted.
#'   xy$trend <- sin(xy$x / 400) + rnorm(n, sd = 0.2)
#'   r2 <- estimate_sac_range(st_as_sf(xy, coords = c("x", "y"), crs = 32632),
#'                            response_var = "trend")
#'   print(r2)
#'   attr(r2, "rejected_reason")
#'   attr(r2, "rejected_range")
#'
#'   # White noise: no structure, so no range.
#'   xy$noise <- rnorm(n)
#'   r3 <- estimate_sac_range(st_as_sf(xy, coords = c("x", "y"), crs = 32632),
#'                            response_var = "noise")
#'   attr(r3, "rejected_reason")
#' }
#' @export
estimate_sac_range <- function(points_sf, response_var,
                               predictor_vars = NULL,
                               n_max = 5000L, cutoff = 0.5,
                               range_frac = 1.0, seed = 123L,
                               detrend = c("ols", "reml"),
                               reml_max_n = 400L,
                               keep_directional_fits = FALSE,
                               lags = c("adaptive", "equal"),
                               n_lags = 15L,
                               structure_alpha = 0.05) {
  detrend <- match.arg(detrend)
  lags    <- match.arg(lags)
  if (!is.numeric(reml_max_n) || length(reml_max_n) != 1L ||
      !is.finite(reml_max_n) || reml_max_n < 30)
    stop("estimate_sac_range(): `reml_max_n` must be a single number of at ",
         "least 30.", call. = FALSE)
  # A bad value used to come back as an unrelated refusal ("no variogram
  # model could be fitted" for cutoff = 0) or as a bare R error that named no
  # argument ("missing value where TRUE/FALSE needed").
  bad_arg <- function(arg, what, x)
    stop(sprintf("estimate_sac_range(): `%s` must be %s; got %s.", arg, what,
                 if (is.null(x)) "NULL"
                 else if (length(x) != 1L) sprintf("%s of length %d", class(x)[1L], length(x))
                 else if (is.numeric(x) || is.logical(x)) format(x)
                 else sprintf("%s (%s)", class(x)[1L], format(x))),
         call. = FALSE)
  one_number <- function(x) is.numeric(x) && length(x) == 1L && !is.na(x)
  if (!one_number(cutoff) || !(cutoff > 0) || cutoff > 1)
    bad_arg("cutoff", paste0("a single number above 0 and at most 1 (the share ",
                             "of the largest distance between two points that ",
                             "the lags run to)"), cutoff)
  if (!one_number(range_frac) || !is.finite(range_frac) || !(range_frac > 0))
    bad_arg("range_frac", "a single finite number above 0", range_frac)
  if (!one_number(n_max) || n_max < 30)
    bad_arg("n_max", paste0("a single number of at least 30 (the fewest points ",
                            "a variogram range is estimated from; Inf for no ",
                            "subsample)"), n_max)
  if (!is.logical(keep_directional_fits) || length(keep_directional_fits) != 1L ||
      is.na(keep_directional_fits))
    bad_arg("keep_directional_fits", "TRUE or FALSE", keep_directional_fits)
  if (!is.null(seed) && (!one_number(seed) || !is.finite(seed) ||
                         abs(seed) > .Machine$integer.max))
    bad_arg("seed", "NULL or a single whole number", seed)
  if (!one_number(n_lags) || !is.finite(n_lags) || n_lags < 5 || n_lags > 100 ||
      n_lags != floor(n_lags))
    bad_arg("n_lags", "a single whole number from 5 to 100", n_lags)
  if (!one_number(structure_alpha) || !(structure_alpha > 0) || structure_alpha > 1)
    bad_arg("structure_alpha", paste0("a single number above 0 and at most 1 ",
                                      "(1 switches the refusal off)"), structure_alpha)
  n_max  <- if (is.finite(n_max)) floor(n_max) else Inf
  n_lags <- as.integer(n_lags)

  # The early NA returns carry their reason as `rejected_reason`, so a caller
  # (make_folds(auto_range = TRUE), kriging_adequacy(), summarize_by_cell())
  # can say why no range came back: the log line alone is invisible under
  # spatialkit_quiet(), knitr and tryCatch().  They carry nothing else -- no
  # variogram, no rejected_range -- which is what tells them from a range that
  # was fitted and refused.
  if (!requireNamespace("gstat", quietly = TRUE)) {
    .log_warn("estimate_sac_range(): package 'gstat' is required for variogram estimation; returning NA.")
    return(structure(NA_real_, rejected_reason =
      "package 'gstat', which the variogram needs, is not installed"))
  }
  if (!inherits(points_sf, "sf"))
    stop("estimate_sac_range(): `points_sf` must be an sf object.", call. = FALSE)
  # MULTIPOINT must be coerced too, not merely admitted: st_coordinates()
  # returns one row per VERTEX, so any multi-vertex feature makes xy[i, ] a
  # different feature than row i, and every fold below misaligns silently.
  # coerce_to_points() takes centroids, matching prep_model_data().
  if (!all(sf::st_geometry_type(points_sf, by_geometry = TRUE) == "POINT"))
    points_sf <- coerce_to_points(points_sf, "auto")

  # gstat::variogram() uses EVERY coordinate dimension, so an XYZ layer has
  # its elevation folded into each lag and the "range" is a length in 3-D,
  # while the block grid, the buffered-LOO buffer, nndm's neighbour distances
  # and summarize_by_cell(deff = "variogram") all work in 2-D map distance.
  points_sf <- sf::st_zm(points_sf, drop = TRUE, what = "ZM")

  # A row with an empty or non-finite geometry makes gstat fail in every
  # direction; it is dropped, with a count, as make_folds() drops it.
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

  if (nrow(points_sf) < 30L) {
    .log_warn("estimate_sac_range(): fewer than 30 points; variogram estimate unreliable. Returning NA.")
    return(structure(NA_real_, rejected_reason = sprintf(
      "%d points, fewer than the 30 a variogram range is estimated from",
      nrow(points_sf))))
  }

  # The variable to model: validate the column first.  as.numeric() on a
  # factor returns its LEVEL CODES, an arbitrary relabelling of the
  # categories, so the range of a factor response changes when the levels are
  # reordered; a character column becomes all-NA and would be blamed on the
  # data rather than on the column type.
  if (!is.character(response_var) || length(response_var) != 1L || is.na(response_var))
    stop("estimate_sac_range(): `response_var` must be a single column name.",
         call. = FALSE)
  df_all <- sf::st_drop_geometry(points_sf)
  if (!(response_var %in% names(df_all)))
    stop(sprintf("estimate_sac_range(): column '%s' not found in `points_sf`.",
                 response_var), call. = FALSE)
  y_all <- df_all[[response_var]]
  if (!is.numeric(y_all)) {
    if (is.logical(y_all)) {
      y_all <- as.numeric(y_all)         # 0/1 is a well-defined variogram target
    } else {
      stop(sprintf(paste0("estimate_sac_range(): response '%s' is %s, and a ",
                          "variogram needs a numeric variable. A factor's ",
                          "codes are an arbitrary relabelling of its levels ",
                          "-- fitting a variogram to them yields a range that ",
                          "changes when the levels are reordered. Encode the ",
                          "column numerically first."),
                   response_var,
                   if (is.factor(y_all)) "a factor" else sprintf("of class %s",
                     paste(class(y_all), collapse = "/"))),
           call. = FALSE)
    }
  }
  has_pred <- !is.null(predictor_vars) && length(predictor_vars) > 0L
  if (has_pred) {
    # A predictor that is not a column is an error, as it is everywhere else
    # in the package: dropped silently, the RAW response would be modelled.
    missing_preds <- setdiff(predictor_vars, names(df_all))
    if (length(missing_preds))
      stop("estimate_sac_range(): predictor_vars ",
           paste(sQuote(missing_preds), collapse = ", "),
           " not found in the data.", call. = FALSE)
  }

  # Rows with nothing to model go first, before the projection is chosen and
  # before the n_max subsample is drawn.  Drawn over every row, the subsample
  # spent its budget on rows that never enter the variogram: 300 valued rows
  # among 20,300 were estimated from the 74 that fell in the draw, and 100
  # among 20,100 came back as "20 point(s) with a finite value".  lm()'s
  # na.exclude drops NA but not Inf, so the rows are chosen here: a finite
  # response and every predictor present and, if numeric, finite.
  fin <- is.finite(y_all)
  if (has_pred) {
    fin_pred <- Reduce(`&`, lapply(predictor_vars, function(v) {
      x <- df_all[[v]]
      if (is.numeric(x)) is.finite(x) else !is.na(x)
    }), TRUE)
    if (!any(fin & fin_pred) && sum(fin) >= 30L) {
      # No row at all has every predictor (one is missing on every row), so
      # no trend can be fitted: the detrending fit fails, as it always did
      # here, and the raw response is modelled.  That is a different
      # estimand, so it is an R warning, not a log line.  When some rows are
      # complete but fewer than 30, there is a trend to fit and too little
      # to fit a variogram to: that is the NA below, with their count.
      .warn_and_log(paste0("estimate_sac_range(): the %s detrending on %s ",
                           "failed (none of the %d row(s) with a response has ",
                           "every predictor present and finite); the ",
                           "variogram is fitted to the RAW response instead, ",
                           "which includes the trend."),
                    if (identical(detrend, "reml")) "REML" else "OLS",
                    paste(predictor_vars, collapse = " + "), sum(fin))
      has_pred <- FALSE
    } else {
      fin <- fin & fin_pred
    }
  }
  n_missing <- sum(!fin)
  if (n_missing > 0L)
    .log_info(if (has_pred)
                paste0("estimate_sac_range(): %d of %d row(s) have a missing or ",
                       "non-finite response or predictor and are left out of ",
                       "the detrending fit and the variogram.")
              else
                paste0("estimate_sac_range(): %d of %d row(s) have a missing or ",
                       "non-finite response and are left out of the variogram."),
              n_missing, length(fin))
  if (sum(fin) < 30L) {
    .log_warn("estimate_sac_range(): too few finite values after filtering; returning NA.")
    return(structure(NA_real_, rejected_reason = sprintf(paste0(
      "%d point(s) with a finite value to model, fewer than the 30 a ",
      "variogram range is estimated from"), sum(fin))))
  }
  if (n_missing > 0L) points_sf <- points_sf[fin, , drop = FALSE]

  # Rows in a canonical order: the coordinates as supplied first, then the
  # response and the predictors to break ties between repeat visits to one
  # location.  Three steps below read the rows in order.  The projection of a
  # lon/lat layer is chosen by scoring candidates on rows taken at even
  # spacing down the layer; the subsample indexes rows; and gstat sums the
  # pairs of a bin in row order.  Unsorted, the same layer with its rows
  # permuted got another subsample and a variogram that differed in its last
  # bits; sorted after the projection, a layer over a continental extent
  # still got an Albers CRS in one row order and a Lambert azimuthal one in
  # another, and a range of 506,838 against 505,530.  So the sort comes
  # before all three.  Nothing below may read meaning into the order of rows
  # that share a location: there it is the order of their values.
  xy_in <- sf::st_coordinates(points_sf)[, 1:2, drop = FALSE]
  df_in <- sf::st_drop_geometry(points_sf)
  ord <- do.call(order, c(
    list(xy_in[, 1L], xy_in[, 2L], as.numeric(df_in[[response_var]])),
    if (has_pred) lapply(predictor_vars, function(v) {
      x <- df_in[[v]]
      if (is.numeric(x) || is.logical(x)) as.numeric(x) else as.character(x)
    }),
    list(method = "radix")))
  points_sf <- points_sf[ord, , drop = FALSE]

  pts <- ensure_projected(points_sf)

  # Subsample if large.  `seed` defaults to a constant, NOT to NULL: the
  # subsample is an internal approximation that keeps an O(n^2) variogram
  # tractable, not part of the answer, so it must not make the range differ
  # between runs on identical input or advance the caller's RNG stream.
  n <- nrow(pts)
  n_subsampled <- 0L
  if (n > n_max) {
    cleanup <- .with_seed(seed)
    on.exit(cleanup(), add = TRUE)
    idx <- sort(sample.int(n, n_max))
    pts <- pts[idx, , drop = FALSE]
    n_subsampled <- as.integer(n - n_max)
    n <- nrow(pts)
    .log_info("estimate_sac_range(): the variogram is computed on a random subsample of %d of the %d point(s) with a value (n_max).",
              n, n + n_subsampled)
  }

  pts_df <- sf::st_drop_geometry(pts)
  y_raw  <- as.numeric(pts_df[[response_var]])
  y <- y_raw
  xy <- sf::st_coordinates(pts)[, 1:2, drop = FALSE]
  detrended <- FALSE
  detrend_method <- NA_character_
  reml_fit <- NULL
  # The design of the trend removed from the response, for the structure
  # test: an intercept alone unless a detrending fit succeeds.
  X_design <- matrix(1, n, 1L)
  if (has_pred) {
    df <- pts_df
    fml <- stats::reformulate(predictor_vars, response_var)
    # REML: trend and covariance fitted together, so the trend is a GLS fit
    # under the fitted correlation and the range is the REML estimate, not a
    # variogram of residuals.  On failure the OLS path below runs instead.
    if (identical(detrend, "reml")) {
      df$.sac_x <- xy[, 1L]; df$.sac_y <- xy[, 2L]
      mf <- df[, c(response_var, predictor_vars, ".sac_x", ".sac_y"), drop = FALSE]
      extent <- max(diff(range(xy[, 1L])), diff(range(xy[, 2L])), 1)
      reml_fit <- .reml_trend(mf, fml, extent = extent, max_n = reml_max_n,
                              seed = seed)
      if (is.null(reml_fit)) {
        .warn_and_log(paste0("estimate_sac_range(): the REML detrending on %s ",
                             "did not converge%s; falling back to OLS ",
                             "detrending, whose residual variogram is biased ",
                             "toward a shorter range (see ?estimate_sac_range)."),
                      paste(predictor_vars, collapse = " + "),
                      if (!requireNamespace("nlme", quietly = TRUE))
                        " (package 'nlme' is not installed)" else "")
      } else {
        X <- try(stats::model.matrix(fml, mf), silent = TRUE)
        b <- if (inherits(X, "try-error")) NA_real_ else reml_fit$beta[colnames(X)]
        if (inherits(X, "try-error") || anyNA(b) || nrow(X) != n) {
          .warn_and_log("estimate_sac_range(): the REML trend coefficients do not match the design; falling back to OLS detrending.")
          reml_fit <- NULL
        } else {
          y <- as.numeric(mf[[response_var]] - X %*% b)
          X_design <- X
          detrended <- TRUE
          detrend_method <- "reml"
          .log_info(paste0("estimate_sac_range(): REML detrending on %s fitted ",
                           "on %d point(s)%s: range parameter %.1f, nugget ",
                           "proportion %.2f."),
                    paste(predictor_vars, collapse = " + "), reml_fit$n_used,
                    if (isTRUE(reml_fit$subsampled)) " (random subsample)" else "",
                    reml_fit$range, reml_fit$nugget_prop)
        }
      }
    }
    lm_fit <- if (is.null(reml_fit)) try(stats::lm(fml, data = df), silent = TRUE) else NULL
    if (is.null(lm_fit)) {
      # REML handled the trend above; nothing to do here.
    } else if (inherits(lm_fit, "try-error")) {
      # Falling through to the raw response is a different estimand, so it
      # is an R warning, not a log line.
      .warn_and_log(paste0("estimate_sac_range(): the OLS detrending on %s ",
                           "failed (%s); the variogram is fitted to the RAW ",
                           "response instead, which includes the trend."),
                    paste(predictor_vars, collapse = " + "),
                    .try_error_message(lm_fit))
    } else {
      resid <- stats::residuals(lm_fit)
      if (length(resid) != n) {
        .warn_and_log("estimate_sac_range(): OLS residual length (%d) does not match data rows (%d); the variogram is fitted to the RAW response instead.",
                      length(resid), n)
      } else {
        y <- as.numeric(resid)
        Xl <- try(stats::model.matrix(lm_fit), silent = TRUE)
        if (!inherits(Xl, "try-error") && is.matrix(Xl) && nrow(Xl) == n) X_design <- Xl
        detrended <- TRUE
        detrend_method <- "ols"
      }
    }
  }

  ok_y <- is.finite(y)
  if (sum(ok_y) < 30L) {
    .log_warn("estimate_sac_range(): too few finite values after filtering; returning NA.")
    return(structure(NA_real_, rejected_reason = sprintf(paste0(
      "%d point(s) with a finite value to model, fewer than the 30 a ",
      "variogram range is estimated from"), sum(ok_y))))
  }
  if (!all(ok_y)) {
    pts <- pts[ok_y, , drop = FALSE]; xy <- xy[ok_y, , drop = FALSE]
    y <- y[ok_y]; y_raw <- y_raw[ok_y]; X_design <- X_design[ok_y, , drop = FALSE]
    n <- length(y)
  }
  # A variable with no variance has no autocorrelation structure to estimate:
  # every semivariance is 0.  There is no range; say so.  "No variance" is
  # judged against the size of the response's own values and against nothing
  # absolute: a spread under 1e-12 of the largest of them is rounding (the
  # residuals of an exact fit come out near 1e-16 of it), and anything wider
  # is data.  The floor used to be sqrt(.Machine$double.eps) times max(1,
  # |mean|), which called a response constant once its unit made its
  # standard deviation smaller than 1.5e-8, and one of standard deviation 1
  # on a mean of 1e8 as well.  The size is divided out before the squares
  # are taken, so neither end of the floating-point range is in the way.
  y_size <- max(abs(y_raw))
  y_mean <- mean(y)
  y_sd   <- if (is.finite(y_size) && y_size > 0) y_size * stats::sd(y / y_size) else 0
  if (!is.finite(y_sd) || !(y_sd > 1e-12 * y_size)) {
    .log_warn(paste0("estimate_sac_range(): the variable being modelled is ",
                     "constant (zero variance%s), so it has no autocorrelation ",
                     "range. Returning NA."),
              if (detrended)
                " -- the residuals are all zero, so the predictors explain the response exactly"
              else "")
    return(structure(NA_real_, rejected_reason = if (detrended)
      "the residuals on predictor_vars are constant: the predictors explain the response exactly"
      else "the response is constant"))
  }
  y_var <- y_sd^2
  # The variogram is computed on the standardised variable and put back on
  # the response's scale afterwards, and the model is fitted on scaled axes
  # (see .sac_fit_model()), so nothing below depends on the unit of the
  # response.  Only the geometry and this one column go to gstat.
  pts <- sf::st_sf(..sac_var = (y - y_mean) / y_sd, geometry = sf::st_geometry(pts))

  # The lag cutoff is a fraction of the maximum inter-point distance, as
  # documented, NOT of the bounding-box diagonal.  The diagonal is a property
  # of the axes, not of the points: it grows by up to sqrt(2) when the same
  # layer is rotated 45 degrees.  The farthest pair lies on the convex hull,
  # whose vertex count is tiny, so this is cheap at any n.
  max_dist <- try({
    hull <- sf::st_convex_hull(sf::st_union(sf::st_geometry(pts)))
    hv   <- unique(sf::st_coordinates(hull)[, 1:2, drop = FALSE])
    if (nrow(hv) < 2L) 0 else max(stats::dist(hv))
  }, silent = TRUE)

  if (inherits(max_dist, "try-error") || !is.finite(max_dist) || max_dist <= 0) {
    .log_warn("estimate_sac_range(): the points have no extent to fit a variogram over. Returning NA.")
    return(structure(NA_real_, rejected_reason = paste0(
      "the points have no extent (the largest distance between them is zero ",
      "or could not be computed)")))
  }

  cutoff_dist <- as.numeric(cutoff * max_dist)

  # --- Is there spatial structure at all? ----------------------------------
  # A variogram model can be fitted to anything, and on white noise (300
  # points) it returned a finite range for 24 of 60 layers.  The
  # test is made before any range is believed; see .sac_structure_test().
  # The response and the trend's design go in, so that residuals are judged
  # against the null that fits them.
  structure_test <- tryCatch(.sac_structure_test(xy, y_raw, X_design),
                             error = function(e) {
    .log_info("estimate_sac_range(): the test for spatial structure could not be computed (%s); no range is refused on its account.",
              conditionMessage(e))
    NULL
  })
  structure_p <- if (is.null(structure_test)) NA_real_ else structure_test$p
  no_structure <- structure_alpha < 1 && is.finite(structure_p) &&
    structure_p >= structure_alpha
  # The semivariance of nearest-neighbour pairs and of pairs at zero
  # separation.  Computed on the standardised variable, as gstat sees it
  # (`zero_ss` stays on that scale: emp_vgm() takes those pairs out of the
  # first bin with it), and reported in the units of the response's variance.
  short_lag <- .sac_short_lag(xy, pts$..sac_var)
  short_lag$gamma      <- short_lag$gamma * y_var
  short_lag$zero_gamma <- short_lag$zero_gamma * y_var
  # The variance between visits to a location is a floor under the nugget:
  # the semivariance at a positive lag has that variance in it, and the
  # variation between locations on top.  It is used when it rests on 30
  # degrees of freedom or more (rows less locations).  Without it, on 100
  # locations visited 5 times with a range of 300 between them, the nugget
  # extrapolated from the positive lags was under the variance between
  # visits in 22 of 27 fits and exactly 0 in 7, while 1000 pairs of visits
  # put that variance at 0.30 (0.28 to 0.32) for a truth of 0.3.
  nugget_floor <- if (isTRUE(short_lag$zero_df >= 30L) &&
                      is.finite(short_lag$zero_gamma)) short_lag$zero_gamma else 0

  # --- Lag bins -------------------------------------------------------------
  # "adaptive": the first bin ends at the distance within which the closest
  # max(50, n / 4) pairs lie and the rest are spaced geometrically out to
  # the cutoff, with bins of fewer than 30 pairs merged into the next.
  # "equal": n_lags bins of equal width, gstat's default at n_lags = 15.
  # Both sets of boundaries exist whichever is asked for: the default looks
  # at the equal bins as well (a second, longer structure; a range the
  # adaptive bins leave past the lags).
  bounds_of <- list(
    adaptive = if (identical(lags, "adaptive"))
      .sac_lag_boundaries(xy, cutoff_dist, n_lags = n_lags) else NULL,
    equal    = cutoff_dist * seq_len(n_lags) / n_lags)
  emp_vgm <- function(alpha = NULL, binning = lags) {
    adaptive <- identical(binning, "adaptive")
    vg <- try(
      if (is.null(alpha)) {
        if (adaptive)
          gstat::variogram(..sac_var ~ 1, data = pts, cutoff = cutoff_dist,
                           boundaries = bounds_of$adaptive)
        else
          gstat::variogram(..sac_var ~ 1, data = pts, cutoff = cutoff_dist,
                           width = cutoff_dist / n_lags)
      } else {
        if (adaptive)
          gstat::variogram(..sac_var ~ 1, data = pts, cutoff = cutoff_dist,
                           boundaries = bounds_of$adaptive, alpha = alpha,
                           tol.hor = 22.5)
        else
          gstat::variogram(..sac_var ~ 1, data = pts, cutoff = cutoff_dist,
                           width = cutoff_dist / n_lags, alpha = alpha, tol.hor = 22.5)
      },
      silent = TRUE)
    if (inherits(vg, "try-error") || !inherits(vg, "data.frame") || !NROW(vg))
      return(NULL)
    # Rows that share a location are pairs at zero separation, which gstat
    # puts in the first bin of every variogram; they are not at a lag.
    vg <- .sac_drop_zero_lag(vg, short_lag$zero_pairs, short_lag$zero_ss,
                             first_bound = bounds_of[[binning]][1L])
    if (!NROW(vg)) return(NULL)
    if (adaptive) vg <- .sac_merge_bins(vg)
    vg$gamma <- vg$gamma * y_var
    vg
  }

  # The largest lag the empirical variogram is fitted over.  A fitted range
  # beyond it is extrapolation past every observed lag rather than a long
  # correlation length; each DIRECTION is tested against it too.
  max_supported <- range_frac * cutoff_dist
  # Top of the range grid, in cutoffs: past any range that can be accepted.
  grid_upper <- max(4, 2 * range_frac)

  # --- Fit a variogram model and return its effective range ----------------
  # An exponential model WITH a nugget.  A nugget-free model forces the curve
  # through the origin, and on any real measurement (which has one) the
  # N/h^2 weights buy that constraint by collapsing the range.  NA when
  # fewer than three bins hold pairs, or when no model with a range fits
  # better than a constant; the model, and whether its minimum lies inside
  # the ranges searched, travel as attributes.
  .fit_vgm_range <- function(vg) {
    m <- .sac_fit_model(vg, cutoff_dist, "Exp", upper = grid_upper,
                        nugget_min = nugget_floor)
    if (!is.data.frame(m)) return(NA_real_)
    sp <- m[as.character(m$model) != "Nug", , drop = FALSE]
    if (nrow(sp) == 0L || !all(is.finite(sp$range)))
      return(structure(NA_real_, flat_model = m))
    eff <- 3 * max(sp$range)
    if (!is.finite(eff) || eff <= 0) return(NA_real_)
    structure(eff, vgm_model = m,
              converged = !identical(attr(m, "converged"), FALSE))
  }

  # --- One binning: the four directional variograms and the all-pairs one ---
  # Directions at 0/45/90/135 degrees with a tolerance of 22.5: gstat counts
  # azimuth clockwise from north (0 = N-S, 90 = E-W), and these four windows
  # tile all 180 azimuths with no overlap and no gap.  They are a diagnostic
  # and never the answer: each sees about a quarter of the pairs, and the
  # windows are fixed to the coordinate axes.
  #
  # The all-pairs variogram is the estimate.  It pools every point pair, so
  # it is the best-powered fit available and the only one that does not
  # change when the layer is rotated.  The directional maximum never stands
  # in for it.  It used to, when this fit was singular or had not converged:
  # on isotropic fields (150 sites, effective range 30 km) that was 6 of 86
  # identified ranges, at a median 3.9 times the truth against 1.4 for the
  # all-pairs answers.
  dir_az <- c(0, 45, 90, 135)
  first_lag_of <- function(vg) {
    if (!is.data.frame(vg)) return(NA_real_)
    d1 <- vg$dist[is.finite(vg$dist) & is.finite(vg$np) & vg$np > 0]
    if (length(d1)) min(d1) else NA_real_
  }
  sweep_lags <- function(binning, vg_iso = emp_vgm(binning = binning),
                         iso_fit = .fit_vgm_range(vg_iso)) {
    dir_fits <- lapply(dir_az, function(az) {
      vg_dir <- emp_vgm(az, binning = binning)
      list(vg = vg_dir, fit = .fit_vgm_range(vg_dir))
    })
    dir_ranges <- vapply(dir_fits, function(f) as.numeric(f$fit), numeric(1))
    dir_conv <- vapply(dir_fits,
                       function(f) !identical(attr(f$fit, "converged"), FALSE),
                       logical(1))
    # A direction counts only if its fit is BOTH finite and identified, i.e.
    # within the longest lag the variogram was fitted over.
    dir_ok <- is.finite(dir_ranges) & dir_ranges <= max_supported & dir_conv
    # The per-azimuth outcome, kept rather than collapsed: no model, a range
    # past the fitted lags and a criterion with no minimum inside the ranges
    # searched are three different findings.
    dir_status <- ifelse(!is.finite(dir_ranges), "no_fit",
                  ifelse(dir_ranges > max_supported, "over_cutoff",
                  ifelse(!dir_conv, "not_converged", "ok")))
    usable <- dir_ranges
    usable[!dir_ok] <- NA_real_
    list(lags = binning, lag_bounds = bounds_of[[binning]],
         vg_iso = vg_iso, iso_fit = iso_fit, first_lag = first_lag_of(vg_iso),
         dir_ranges = dir_ranges, dir_ok = dir_ok,
         dir_success = sum(dir_ok) >= 2L,
         dir_status = stats::setNames(dir_status, as.character(dir_az)),
         dir_fitted = stats::setNames(dir_ranges, as.character(dir_az)),
         dir_detail = stats::setNames(lapply(dir_fits, function(f)
           list(variogram = if (inherits(f$vg, "data.frame")) f$vg else NULL,
                model     = attr(f$fit, "vgm_model"))), as.character(dir_az)),
         usable = usable,
         anisotropy = if (sum(dir_ok) >= 2L)
           max(usable, na.rm = TRUE) / min(usable, na.rm = TRUE) else NA_real_)
  }
  # What the refusals further down will say of a range fitted to a
  # variogram, in their order: "none" (no model with a range), "decreasing",
  # "over" (past the fitted lags; a fit with no minimum inside the ranges
  # searched ends here or at "under"), "under" (below the first lag), "ok".
  range_verdict <- function(fit, vg, first_lag) {
    r <- as.numeric(fit)
    if (!is.finite(r) || r <= 0) return("none")
    if (.variogram_decreasing(vg)) return("decreasing")
    if (is.finite(max_supported) && r > max_supported) return("over")
    if (is.finite(first_lag) && r < first_lag) return("under")
    "ok"
  }

  S <- sweep_lags(lags)

  # --- What the equal lags add to the adaptive ones -------------------------
  # One structure fitted to lags that resolve the closest pairs follows the
  # SHORTEST structure of the layer, and those lags leave little weight on
  # the long ones.  Two things follow, and the equal-lag fit of the all-pairs
  # variogram answers both.
  #
  # A second, longer structure.  A layer that also has a weaker structure
  # several times longer would have its blocks sized from the short one,
  # which is the direction that leaks.  Equal lags do not see the short
  # structure and fit the long one.  When the equal-lag range is at least
  # twice the adaptive one, a model with both structures, at those two
  # ranges, is fitted to the adaptive variogram (.sac_fit_nested()).  It is
  # kept, and the range is the long one, unless the long structure gets
  # under 3 percent of that model's semivariance at the largest lag: then
  # the lags that resolve the layer give the equal-lag range nothing, and
  # the single structure stands.  The bar is low on purpose.  With a ratio
  # of 3 and a share of 10 percent the long structure was dropped where the
  # field has one: on a clustered layer of 1200 points over a 1000 m square
  # (30 draws each) with structures of range 20 and 200, the range returned
  # was under half the equal-lag one in 8 draws, and is in 1 now; with 40
  # and 400, in 8 and 1; with 30 and 300, in 8 and 0.  The price is paid on
  # layers with one short structure, where a noisy equal-lag range is
  # believed more often: with one structure of range 40 a second is kept in
  # 4 of 30 draws (0 before) and the range is within a factor 2 of the
  # truth in 27 (29 before, 19 with equal lags); with range 30, 7 and 25 (0
  # and 29, and 14 with equal lags).  That error is on the long side, where
  # blocks are too large and a validation pessimistic, and `short_range`
  # still holds the short structure.  When the equal-lag range is one that
  # lags = "equal" would itself refuse (past the fitted lags, most often),
  # keeping it takes the answer away instead of lengthening it, and the bar
  # stays at 10 percent: an unremoved trend clears it (of 12 draws of a
  # short structure under a linear trend the equal-lag range ran past the
  # lags in 11, and 9 of those cleared it), a noisy long lag mostly does not.
  #
  # A range the adaptive lags leave past the largest lag.  A weak, long
  # structure under much noise (a nugget of half the sill or more) can show
  # at equal lags and not at adaptive ones, whose long lags carry little
  # weight.  When the adaptive range is refused as past the fitted lags and
  # the equal-lag fit has a range it would itself return, the equal-lag
  # result is the answer, whole: its variogram, its model, its directional
  # sweep and its first lag, as lags = "equal" gives them, with
  # `lag_fallback` set.  Of 120 such layers (four designs) that was 14.
  # Nothing falls back from the other refusals: a range under the first
  # adaptive lag or a flat adaptive variogram is no structure the wider
  # equal bins can resolve, and on white noise that passed the structure
  # test their equal-lag fits are what came back finite.
  short_range  <- NA_real_
  lag_fallback <- FALSE
  if (identical(lags, "adaptive") && is.null(reml_fit)) {
    ada <- range_verdict(S$iso_fit, S$vg_iso, S$first_lag)
    if (ada %in% c("ok", "over")) {
      eq_vg  <- emp_vgm(binning = "equal")
      eq_fit <- .fit_vgm_range(eq_vg)
      eq_conv <- !identical(attr(eq_fit, "converged"), FALSE)
      # Whether lags = "equal" would return the equal-lag range.
      eq_ok <- eq_conv &&
        identical(range_verdict(eq_fit, eq_vg, first_lag_of(eq_vg)), "ok")
      if (identical(ada, "ok") && is.finite(eq_fit) &&
          as.numeric(eq_fit) >= 2 * as.numeric(S$iso_fit)) {
        two <- .sac_fit_nested(S$vg_iso, cutoff_dist,
                               r_short = as.numeric(S$iso_fit),
                               r_long = as.numeric(eq_fit),
                               min_share = if (eq_ok) 0.03 else 0.10,
                               converged = eq_conv, nugget_min = nugget_floor)
        if (is.data.frame(two)) {
          short_range <- 3 * two$range[2L]
          .log_info(paste0("estimate_sac_range(): structure at two scales. The ",
                           "equal-lag fit has an effective range of %.4g, %.1f ",
                           "times the %.4g fitted on the adaptive lags, and a ",
                           "structure of that range carries %.0f%% of the ",
                           "semivariance at the largest lag. The model has both ",
                           "(short range %.4g) and the longer range is the ",
                           "estimate.%s"),
                    as.numeric(eq_fit), as.numeric(eq_fit) / as.numeric(S$iso_fit),
                    as.numeric(S$iso_fit), 100 * attr(two, "long_share"),
                    short_range,
                    if (eq_conv) "" else paste0(
                      " The equal-lag fit has no minimum inside the ranges ",
                      "searched, so the long range and its partial sill are the ",
                      "edge of the search and the model is marked converged = ",
                      "FALSE."))
          S$iso_fit <- structure(as.numeric(eq_fit), vgm_model = two,
                                 converged = eq_conv)
        }
      } else if (identical(ada, "over") && eq_ok) {
        .log_info(paste0("estimate_sac_range(): the range fitted on the ",
                         "adaptive lags (%.4g) is past the largest lag fitted ",
                         "(%.4g), and the equal-lag fit has a range inside the ",
                         "lags (%.4g). The equal-lag result is returned, as ",
                         "lags = \"equal\" gives it (attr(range, ",
                         "\"lag_fallback\") is TRUE): its nugget is ",
                         "extrapolated from a first bin %.4g wide."),
                  as.numeric(S$iso_fit), max_supported, as.numeric(eq_fit),
                  cutoff_dist / n_lags)
        S <- sweep_lags("equal", vg_iso = eq_vg, iso_fit = eq_fit)
        lag_fallback <- TRUE
      }
    }
  }
  # The binning behind the result: the one asked for unless it fell back.
  lags_used   <- S$lags
  lag_bounds  <- S$lag_bounds
  vg_iso      <- S$vg_iso
  iso_fit     <- S$iso_fit
  first_lag   <- S$first_lag
  dir_ranges  <- S$dir_ranges
  dir_ok      <- S$dir_ok
  dir_success <- S$dir_success
  dir_status  <- S$dir_status
  dir_fitted  <- S$dir_fitted
  usable      <- S$usable
  anisotropy  <- S$anisotropy
  # Four empirical variograms are most of what the object weighs and nothing
  # in the package reads them, so they are kept on request only.
  dir_detail_out <- if (isTRUE(keep_directional_fits)) S$dir_detail else NULL

  vg_used  <- vg_iso
  vgm_used <- NULL
  effective_range <- NA_real_
  if (!is.null(reml_fit)) {
    # --- REML detrending: the range is the REML estimate -------------------
    # Not a variogram fitted to the REML residuals: a residual variogram is
    # biased toward a shorter range whatever fitted the trend, and the REML
    # parameters are the estimate that does not carry that bias.  The
    # empirical variogram of the residuals is attached for inspection, and
    # the model line drawn through it is the REML model.
    effective_range <- 3 * reml_fit$range
    vgm_used <- gstat::vgm(psill  = reml_fit$sigma2 * (1 - reml_fit$nugget_prop),
                           model  = "Exp", range = reml_fit$range,
                           nugget = reml_fit$sigma2 * reml_fit$nugget_prop)
  } else if (is.finite(iso_fit)) {
    effective_range <- as.numeric(iso_fit)
    vgm_used <- attr(iso_fit, "vgm_model")
  } else {
    # No structured model: a flat variogram (a pure nugget is the best fit),
    # or too few bins to fit one to.
    vgm_used <- attr(iso_fit, "flat_model")
  }

  # The nugget of whatever model stands behind the answer, identified or not:
  # its semivariance as the separation shrinks to zero.  A fit that runs
  # into the nugget's lower bound reports 0 (or, by REML, a value within
  # rounding of it), which is the bound and not an estimate; it is marked,
  # never handed on bare.
  nugget_val <- .vgm_nugget_of(vgm_used)
  sill_val   <- if (is.data.frame(vgm_used)) sum(as.numeric(vgm_used$psill)) else NA_real_
  nugget_at_bound <- if (is.finite(nugget_val) && is.finite(sill_val) && sill_val > 0)
    nugget_val <= 1e-4 * sill_val else NA
  # Bins inside the fitted range; with two structures, inside the short one,
  # which is the one the nugget is told apart from.
  bins_below <- if (is.finite(effective_range) && is.data.frame(vg_used))
    sum(vg_used$dist < (if (is.finite(short_range)) short_range else effective_range))
  else NA_integer_
  # The REML fit's own numbers, on the refusals as on the success.
  reml_attr <- if (is.null(reml_fit)) NULL else
    list(n_used = reml_fit$n_used, subsampled = isTRUE(reml_fit$subsampled),
         nugget_prop = reml_fit$nugget_prop, sigma2 = reml_fit$sigma2,
         logLik = reml_fit$logLik,
         logLik_gain = reml_fit$logLik - reml_fit$logLik_null,
         n_fits = reml_fit$n_fits)

  # Every classed return is built here, so each carries the same evidence:
  # the directional sweep, the variogram, the model, the nugget and how far
  # it can be trusted, the structure test, and what the lags and the row
  # counts were.
  finish <- function(value, ...) {
    structure(
      value,
      class           = c("sac_range", "numeric"),
      directional     = stats::setNames(usable, as.character(dir_az)),
      anisotropy      = anisotropy,
      # Always FALSE: the all-pairs fit is the estimate or there is none.
      # Kept so that code reading it, and objects stored with it, still work.
      anisotropy_used = FALSE,
      directional_status = dir_status,
      directional_fitted = dir_fitted,
      directional_fits   = dir_detail_out,
      detrended       = isTRUE(detrended),
      detrend_method  = detrend_method,
      reml            = reml_attr,
      max_dist        = as.numeric(max_dist),
      cutoff_dist     = as.numeric(cutoff_dist),
      # The CRS the variogram was fitted in: the unit of the range.
      crs             = sf::st_crs(pts),
      variogram       = vg_used,
      variogram_model = vgm_used,
      nugget          = nugget_val,
      nugget_at_bound = nugget_at_bound,
      first_lag       = as.numeric(first_lag),
      bins_below_range = bins_below,
      short_lag_semivariance = short_lag$gamma,
      short_lag_pairs = short_lag$pairs,
      short_lag_dist  = short_lag$dist,
      zero_lag_semivariance = short_lag$zero_gamma,
      zero_lag_pairs  = short_lag$zero_pairs,
      # The floor the variogram fit held the nugget to; the REML model is
      # not fitted to the variogram and has none.
      nugget_floor    = if (is.null(reml_fit)) nugget_floor else 0,
      structure_p     = structure_p,
      morans_i        = if (is.null(structure_test)) NA_real_ else structure_test$morans_i,
      structure_test  = if (is.null(structure_test)) NULL else structure_test$table,
      structures      = if (is.data.frame(vgm_used))
        sum(as.character(vgm_used$model) != "Nug") else NA_integer_,
      short_range     = as.numeric(short_range),
      lags            = lags_used,
      lag_fallback    = lag_fallback,
      lag_boundaries  = as.numeric(lag_bounds),
      n_used          = as.integer(n),
      n_missing       = as.integer(n_missing),
      n_subsampled    = as.integer(n_subsampled),
      ...
    )
  }

  # --- No spatial structure: no range, whatever was fitted ------------------
  if (no_structure) {
    best <- structure_test$table[which.min(structure_test$table$p), , drop = FALSE]
    .log_warn(
      paste0("estimate_sac_range(): no spatial structure detected at the ",
             "nearest-neighbour scale: Moran's I of %s on the %d nearest ",
             "neighbours is %.3f (z = %.2f), the strongest of the %d ",
             "neighbourhoods tested, and the adjusted p-value is %.3f, not ",
             "below structure_alpha = %s.%s Neighbouring values are no more ",
             "alike than distant ones, so there is no range to estimate. ",
             "Returning NA. plot() the returned value to see the variogram; ",
             "set structure_alpha = 1 to fit a range regardless."),
      if (detrended) "the residuals on predictor_vars" else "the response",
      best$k, best$morans_i, best$z, nrow(structure_test$table), structure_p,
      format(structure_alpha),
      if (isTRUE(structure_test$n_locations < n))
        sprintf(paste0(" The %d rows lie at %d locations, and the test is ",
                       "made on the location means: rows that share a ",
                       "location are not neighbours."),
                n, structure_test$n_locations)
      else "")
    # The model the test supports is a pure nugget at the sample variance.
    # Whatever was fitted is not attached: a spurious fit's nugget was under
    # half the sill in 7 of 60 white-noise layers at n = 300.
    vgm_used <- gstat::vgm(psill = y_var, model = "Nug", range = 0)
    attr(vgm_used, "converged") <- TRUE
    nugget_val <- y_var
    nugget_at_bound <- FALSE
    bins_below <- NA_integer_
    short_range <- NA_real_
    return(finish(NA_real_, rejected_range = NA_real_,
                  rejected_reason = .sac_reason_no_structure))
  }

  if (is.null(reml_fit) && !is.finite(iso_fit)) {
    if (is.data.frame(vgm_used)) {
      .log_warn(paste0("estimate_sac_range(): the empirical variogram is flat: ",
                       "no model with a range fits it better than a constant ",
                       "(%.3g), so there is no range to estimate at the lags ",
                       "resolved (%.3g to %.3g). Returning NA. The empirical ",
                       "variogram is attached for inspection: call plot() on ",
                       "the returned value."),
                nugget_val, first_lag, cutoff_dist)
      return(finish(NA_real_, rejected_range = NA_real_,
                    rejected_reason = .sac_reason_flat))
    }
    .log_warn(paste0("estimate_sac_range(): no variogram model could be fitted: ",
                     "fewer than three lag bins hold point pairs; returning NA. ",
                     "The empirical variogram is attached for inspection: call ",
                     "plot() on the returned value."))
    return(finish(NA_real_, rejected_range = NA_real_,
                  rejected_reason = "no variogram model could be fitted (singular fits)"))
  }

  if (any(!dir_ok) && dir_success)
    .log_info(paste0("estimate_sac_range(): %d of %d directional variograms ",
                     "did not yield an identified range (%s) and were ",
                     "excluded from the directional diagnostic."),
              sum(!dir_ok), length(dir_az),
              paste(sprintf("%d\u00b0", dir_az[!dir_ok]), collapse = ", "))

  if (!is.finite(effective_range) || effective_range <= 0) {
    .log_warn("estimate_sac_range(): estimated range is non-positive or non-finite; returning NA.")
    return(finish(NA_real_, rejected_range = as.numeric(effective_range),
                  rejected_reason = "fitted range is non-positive or non-finite"))
  }

  # --- Reject ranges the data cannot actually support -----------------------
  # A model fit returns a finite number even when the empirical variogram
  # never reaches a sill.  The range is then unidentified and the value is a
  # fitting artefact, not an autocorrelation range.  The test is against the
  # CUTOFF distance, not the extent: a fitted range beyond the longest lag
  # observed is extrapolation.  Returning such a value silently is worse than
  # returning NA, because make_folds(auto_range = TRUE) sizes spatial blocks
  # from it: a range spanning the data yields one block covering everything.
  over_cutoff <- is.finite(max_supported) && effective_range > max_supported
  # A semivariance that FALLS with distance over the shorter lags has no sill
  # to reach, so the range fitted through it is not identified either.  This
  # is not the trend signature (a trend makes the curve rise without a sill,
  # which `over_cutoff` catches) but a periodic (hole-effect) structure, or
  # a variance that differs between a dense cluster and the rest of the layer.
  decreasing  <- .variogram_decreasing(vg_used)
  # The mirror of `over_cutoff` at the other end of the lags: a range too
  # short for enough pairs of points to lie inside it is not identified, and
  # the data cannot tell it from a nugget.  For the variogram fit the bound
  # is the shortest lag the variogram resolves (the mean separation in its
  # first bin).  The REML range is fitted to the point pairs themselves, not
  # to the bins, so its bound is the distance within which the points it
  # used have 30 pairs (the usual minimum per lag), capped at that first lag
  # so that it never refuses what the bin bound accepts.
  range_floor <- first_lag
  floor_is_pairs <- FALSE
  if (!is.null(reml_fit) && is.finite(reml_fit$pair_floor %||% NA_real_) &&
      !isTRUE(reml_fit$pair_floor >= first_lag)) {
    range_floor    <- reml_fit$pair_floor
    floor_is_pairs <- TRUE
  }
  under_lag <- is.finite(range_floor) && effective_range < range_floor

  if (over_cutoff || decreasing || under_lag) {
    # What to try next.  "Supply `predictor_vars`" must not be said of a
    # range that is already of the residuals on them.  REML is suggested only
    # where it was not asked for (a REML fit that failed has fallen back to
    # OLS).
    detrend_advice <- if (isTRUE(detrended))
      paste0("try predictors that carry the trend (coordinate terms, for ",
             "instance)",
             if (identical(detrend_method, "ols") && !identical(detrend, "reml"))
               " or `detrend = \"reml\"`" else "")
    else "supply `predictor_vars` to detrend"
    if (decreasing) {
      .log_warn(
        paste0("estimate_sac_range(): the empirical variogram decreases with ",
               "distance over its shorter lags, so no range is identified ",
               "from it (the fit reported %.0f). That shape is what a periodic ",
               "(hole-effect) structure produces, or a variance that differs ",
               "between a dense cluster and the rest of the layer; an ",
               "unremoved trend makes a variogram rise without a sill, which ",
               "is a different signal.%s Returning NA. Inspect it with plot() on ",
               "the returned value, and set a block size explicitly."),
        effective_range,
        # The fixed 15% tolerance does not widen with the sampling noise of
        # the short-lag bins, so small samples trip it on ordinary fields.
        if (n < 100L)
          sprintf(paste0(" With %d points the short-lag bins are noisy enough ",
                         "for an ordinary field to show this shape as well."), n)
        else ""
      )
    } else if (over_cutoff) {
      .log_warn(
        paste0("estimate_sac_range(): fitted range (%.0f) exceeds the largest ",
               "lag the variogram was fitted over (%.4g = %s x cutoff %.0f); the ",
               "empirical variogram never reached a sill, so the range is ",
               "unidentified rather than long. Returning NA. Raise `cutoff` to ",
               "fit longer lags, %s, or set a block size explicitly. (A ",
               "variogram that is flat from the first lag, with no spatial ",
               "structure to find, can end here too: plot() the returned value ",
               "to tell the two apart.)"),
        effective_range, max_supported, format(range_frac), cutoff_dist,
        detrend_advice
      )
    } else if (!is.null(reml_fit)) {
      .log_warn(
        paste0("estimate_sac_range(): the REML range (%.3g) is shorter than ",
               "%s (%.3g): fewer than 30 pairs of the %d points the REML fit ",
               "used lie inside it, too few to identify a range, and a ",
               "structure that short cannot be told from a nugget. Returning ",
               "NA."),
        effective_range,
        if (floor_is_pairs) "the distance within which 30 of its point pairs lie"
        else paste0("the shortest lag the variogram of the residuals resolves ",
                    "(the mean separation in its first bin)"),
        range_floor, as.integer(reml_fit$n_used)
      )
    } else if (identical(lags_used, "equal")) {
      .log_warn(
        paste0("estimate_sac_range(): the fitted range (%.3g) is shorter than ",
               "the shortest lag the variogram resolves (%.3g, the mean ",
               "separation in its first bin), so no lag inside it was fitted ",
               "and it cannot be told from a nugget. Returning NA. If the ",
               "empirical variogram is at its sill from the first bin, re-run ",
               "with a smaller `cutoff` (the bins are cutoff * max_dist / ",
               "n_lags wide) or with lags = \"adaptive\", which starts the ",
               "bins at the closest pairs; plot() the returned value to see it."),
        effective_range, range_floor
      )
    } else {
      .log_warn(
        paste0("estimate_sac_range(): the fitted range (%.3g) is shorter than ",
               "the shortest lag the variogram resolves (%.3g, the mean ",
               "separation of the closest pairs of points, which fill its ",
               "first bin), so no lag inside it was fitted and it cannot be ",
               "told from a nugget. Returning NA. The points are too far apart ",
               "to measure a range that short; plot() the returned value to ",
               "see the variogram."),
        effective_range, range_floor
      )
    }
    # The VALUE is NA: the range is unidentified and must not size blocks.
    # The variogram that justified the rejection travels with it, so the
    # refusal can be looked at: attributes do not affect is.finite()/is.na().
    return(finish(
      NA_real_,
      rejected_range  = as.numeric(effective_range),
      rejected_reason = if (decreasing)
        "empirical variogram decreases with distance"
      else if (over_cutoff)
        "fitted range exceeds the largest lag fitted"
      else "fitted range is below the shortest lag fitted",
      # The distance the refused range fell short of, for that refusal only.
      range_floor     = if (!decreasing && !over_cutoff)
        as.numeric(range_floor) else NULL
    ))
  }

  if (is.finite(anisotropy) && anisotropy > 1.5 && is.null(reml_fit))
    .log_info(
      paste0("estimate_sac_range(): the directional ranges vary by a factor ",
             "of %.1f (%s). Each direction sees about a quarter of the ",
             "point pairs and the windows are fixed to the coordinate axes, ",
             "so this spread is expected on an isotropic field too; the ",
             "all-directions estimate (%.1f) is used. If the field is known ",
             "to be anisotropic, size blocks from the longest directional ",
             "range instead: max(attr(range, \"directional_fitted\"), ",
             "na.rm = TRUE), after checking attr(range, ",
             "\"directional_status\"), since a direction that is not \"ok\" ",
             "has no identified range."),
      anisotropy,
      paste(sprintf("%d\u00b0 = %.1f", dir_az[dir_ok], dir_ranges[dir_ok]),
            collapse = ", "),
      effective_range)
  if (isTRUE(nugget_at_bound))
    .log_info(paste0("estimate_sac_range(): the fitted nugget (%.3g on a sill ",
                     "of %.3g) is at its lower bound, so it is a bound and not ",
                     "an estimate (attr(range, \"nugget_at_bound\") is TRUE). ",
                     "The semivariance of nearest-neighbour pairs, %.3g at a ",
                     "median separation of %.3g, bounds the nugget from above."),
              nugget_val, sill_val, short_lag$gamma, short_lag$dist)

  finish(effective_range)
}


#' The nugget psill of a gstat variogram model, or NA
#' @keywords internal
#' @noRd
.vgm_nugget_of <- function(vm) {
  if (!is.data.frame(vm) || !all(c("model", "psill") %in% names(vm))) return(NA_real_)
  v <- sum(as.numeric(vm$psill[as.character(vm$model) == "Nug"]))
  if (length(v) != 1L || !is.finite(v)) NA_real_ else v
}


#' Whether the fit behind a variogram model found a minimum, or NA
#'
#' estimate_sac_range() marks every model it fits with a `converged`
#' attribute (see .sac_fit_model()): FALSE when the weighted least-squares
#' criterion is still falling at an edge of the ranges searched, so that the
#' range, the sill and the nugget are that edge and not a minimum.  That is
#' a variogram with no sill inside four cutoffs (or a structure far below
#' the first lag), and its nugget is not a fitted value either.
#' `rejected_reason` names one ground for a refusal and does not say this:
#' such a fit is "fitted range exceeds the largest lag fitted" (or "... is
#' below the shortest lag fitted"), as is a fit with a proper minimum just
#' past the lags, whose nugget is sound.  The flag is what tells the two
#' apart, and resolution_profile() reads it.  A model with two structures
#' carries the flag of the equal-lag fit that gave the long range (see
#' .sac_fit_nested()): FALSE there means the long range and its partial sill
#' are the edge of that search, so the model's sill is not a fitted value.
#' Its nugget and short structure were fitted at the adaptive lags and are
#' sound, but the flag is one bit for the whole model and is read as it
#' always was: nothing is taken from a model marked FALSE.  NA when the
#' model says nothing
#' either way: one built from REML parameters, one made by hand, or no model
#' at all.  Objects stored by earlier versions carry gstat's verdict here
#' (FALSE when its optimiser stopped at the iteration limit), which is read
#' the same way.
#' @keywords internal
#' @noRd
.vgm_converged <- function(vm) {
  v <- attr(vm, "converged", exact = TRUE)
  if (is.logical(v) && length(v) == 1L && !is.na(v)) v else NA
}


#' Print a spatial autocorrelation range
#'
#' Prints the effective range as a plain number, with the directional fit
#' summarised beneath it when one is available.  A direction whose fit was
#' unusable is labelled with why (\code{directional_status}) and the range
#' its fit reported (\code{directional_fitted}) when the object carries
#' them, and \code{unidentified} otherwise.  A last line names the unit and
#' the CRS the range is a length in (\code{attr(x, "crs")}, which for
#' lon/lat input is the projected CRS the estimate chose; for a layer with
#' no CRS, it says the range is in that layer's own coordinate units), and
#' whether the variogram is of the response or of its residuals on
#' \code{predictor_vars} (\code{detrended}, \code{detrend_method}).
#'
#' @param x An object of class \code{sac_range}.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @family print methods
#' @export
print.sac_range <- function(x, ...) {
  cat(format(as.numeric(x)), "\n")
  d <- attr(x, "directional")
  a <- attr(x, "anisotropy")
  st <- attr(x, "directional_status")
  ft <- attr(x, "directional_fitted")
  # any(), not all(): a direction whose variogram never reached a sill is
  # recorded as NA and excluded from the maximum, and suppressing the whole
  # line in that case hides exactly the diagnostic worth seeing.  With the
  # per-azimuth status attached the line is worth printing even when no
  # direction was usable: four ranges past the fitted lags say something.
  if (!is.null(d) && length(d) > 0L &&
      (any(is.finite(d)) || (!is.null(st) && length(st) == length(d)))) {
    # Names, not fixed positions: the azimuth sweep is c(0, 45, 90, 135) and
    # an older stored object may carry only c(0, 90).
    labs <- names(d)
    if (is.null(labs)) labs <- as.character(seq_along(d) - 1L)
    failed <- vapply(seq_along(d), function(i) {
      if (is.null(st) || length(st) != length(d) || is.null(ft) ||
          length(ft) != length(d)) return("unidentified")
      why <- switch(as.character(st[[i]]),
                    over_cutoff   = "past the fitted lags",
                    not_converged = "not converged",
                    no_fit        = "no fit",
                    "unidentified")
      if (is.finite(ft[[i]])) sprintf("%s (%s)", format(unname(ft[[i]])), why)
      else why
    }, character(1))
    cat("  directional: ",
        paste(sprintf("%s deg = %s", labs,
                      ifelse(is.finite(d), format(unname(d)), failed)),
              collapse = ", "),
        sep = "")
    if (is.finite(a)) cat(sprintf("  (ratio %.2f)", a))
    cat("\n")
  }
  # What the number is a length in, and of what.  Lon/lat input is fitted in a
  # UTM or equal-area CRS picked for it, and a detrended estimate is the range
  # of the residuals rather than of the response.  Both ride as attributes,
  # and not every function that takes the object reconciles them with its
  # own data, so they are shown where a mismatch can be seen before the
  # object is passed on.
  cr <- attr(x, "crs")
  dm <- attr(x, "detrend_method")
  what <- if (isTRUE(attr(x, "detrended")))
    sprintf("; variogram of the residuals on predictor_vars (%s)",
            if (is.character(dm) && length(dm) == 1L && !is.na(dm)) dm else "detrended")
  else if (isFALSE(attr(x, "detrended"))) "; variogram of the response itself"
  else ""
  if (inherits(cr, "crs") && !is.na(cr)) {
    u <- tryCatch(cr$units_gdal, error = function(e) NULL)
    unit <- if (is.character(u) && length(u) == 1L && !is.na(u) && nzchar(u))
      switch(u, metre = "metres", kilometre = "kilometres", foot = "feet",
             "US survey foot" = "US survey feet", degree = "degrees", u)
    else "CRS units"
    cat(sprintf("  in %s of %s%s\n", unit, .fold_crs_label(cr), what))
  } else if (inherits(cr, "crs")) {
    # A layer with no CRS: the range is in its own, unnamed coordinate units,
    # and whether it is of residuals matters as much as with one.
    cat(sprintf("  in the coordinate units of a layer with no CRS%s\n", what))
  }
  invisible(x)
}


#' Compute grid dimensions that respect a minimum block size
#'
#' Given a bounding box and a minimum block edge length, returns nx/ny
#' values such that each cell is at least \code{block_size} long on every
#' axis of the extent that is at least that long; a shorter axis gets one
#' row or column, whose cells are as long as the extent there.
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

#' Estimate an autocorrelation range for a diagnostic, quietly
#'
#' \code{make_folds()} calls this when \code{auto_range} is off but a response
#' column is available, so that the leakage warning on geometric blocks can
#' fire.  The range returned here sizes nothing.
#'
#' Skips the estimate (returning \code{NA}) when \code{gstat} is not
#' installed or there are fewer than 30 points (the estimator's own floor),
#' logging why at INFO level.  Otherwise runs \code{estimate_sac_range()} with
#' its console echo silenced and any R warning it raises muffled, because the
#' caller did not ask for a range and an unidentified one only means the
#' diagnostic cannot be given.  A failure of any other kind also yields
#' \code{NA}.
#'
#' @return A single number (possibly \code{NA}), stripped of the
#'   \code{sac_range} class so that nothing downstream mistakes it for a
#'   user-requested estimate.
#' @keywords internal
#' @noRd
.sac_range_for_diagnostic <- function(pts, response_var, predictor_vars,
                                      range_frac, seed) {
  if (!requireNamespace("gstat", quietly = TRUE)) {
    .log_info("make_folds(block_kfold): leakage check skipped: package 'gstat' is not installed, so no autocorrelation range could be estimated to compare the blocks against.")
    return(NA_real_)
  }
  if (nrow(pts) < 30L) {
    .log_info("make_folds(block_kfold): leakage check skipped: fewer than 30 points, below the floor at which a variogram range is estimated.")
    return(NA_real_)
  }
  r <- tryCatch(
    logger::with_log_threshold(
      withCallingHandlers(
        estimate_sac_range(pts, response_var = response_var,
                           predictor_vars = predictor_vars,
                           range_frac = range_frac, seed = seed),
        warning = function(w) invokeRestart("muffleWarning")),
      threshold = logger::FATAL, namespace = "spatialkit", index = 2),
    error = function(e) NA_real_)
  r <- suppressWarnings(as.numeric(r))
  if (length(r) != 1L || !is.finite(r) || r <= 0) NA_real_ else r
}

#' Create spatial cross-validation folds
#'
#' Builds train/test splits by random k-fold, spatial block k-fold,
#' leave-location-out, buffered leave-one-out, nearest-neighbour distance
#' matching (NNDM) leave-one-out or its k-fold version (kNNDM).
#'
#' For \code{block_kfold}, the default grid sizing is purely geometric and
#' unrelated to the autocorrelation range of the data.  When blocks are
#' smaller than the autocorrelation range, spatially correlated observations
#' leak across folds and CV metrics become optimistic.  Use
#' \code{block_size} to set a minimum block edge length (in CRS units), or
#' set \code{auto_range = TRUE} to estimate the range from an empirical
#' variogram and enforce it automatically.  When a response is given and the
#' blocks come out smaller than the estimated range, the warning says what
#' range-sized blocks would yield on this layout before recommending them:
#' the setting and the grid it builds when that gives \code{k} folds, and
#' otherwise how many folds it gives (or that it leaves a single block), the
#' largest block size that still gives \code{k} folds, and
#' \code{method = "nndm"} as the design that needs no block size.
#'
#' @param points_sf An sf object.  Any Z or M dimension is dropped before
#'   folding: \code{sf::st_distance()} uses every coordinate dimension, so an
#'   XYZ layer would otherwise have elevation folded into every buffer, block
#'   and neighbour distance.  CRS-less points are aligned to a \code{boundary}
#'   or \code{prediction_points} that carries a CRS.  They are reprojected when
#'   the coordinates look like lon/lat and otherwise stamped without
#'   reprojection, with a warning either way.
#' @param k Integer; number of folds.  Must be a single whole number >= 1,
#'   and is required except for the two leave-one-out methods.
#'   A fraction, \code{NA} or a vector is an error, because a non-integer used
#'   to truncate silently and leave the last rows in no test set at all.
#'   Not every method honours it.  \code{"buffered_loo"} and \code{"nndm"} are
#'   leave-one-out schemes and always return one fold per held-out point
#'   (\code{k = n}, or \code{max_test} when that is given) regardless of what
#'   was asked for; \code{"block_kfold"} lowers it when the grid yields fewer
#'   than \code{k} non-empty blocks, and \code{"leave_location_out"} lowers it
#'   when there are fewer than \code{k} distinct groups.  \code{"random_kfold"}
#'   and \code{"knndm"} lower it only when the layer has fewer than \code{k}
#'   rows.  \code{k = 1} is
#'   raised to 2 by \code{"random_kfold"}, \code{"block_kfold"},
#'   \code{"leave_location_out"} and \code{"knndm"}, since one fold has no
#'   training set.  Read
#'   the \code{k} element of the returned list, and do not assume the
#'   requested value; \code{params$k_requested} keeps the value asked for.
#'   Returning fewer folds than asked raises a warning of class
#'   \code{"spatialkit_k_reduced"}, whose message gives both counts, the
#'   reason and a setting that works; the condition carries them as
#'   \code{k_requested} and \code{k}, so
#'   \code{tryCatch(spatialkit_k_reduced = )} or
#'   \code{withCallingHandlers()} can act on it without reading the text.
#'   (Before this version a reduction was a log line only.)
#' @param method One of \code{"random_kfold"}, \code{"block_kfold"},
#'   \code{"buffered_loo"}, \code{"leave_location_out"}, \code{"nndm"} or
#'   \code{"knndm"}.  See \strong{Details} for what each one does and when it
#'   is appropriate.
#' @param seed Optional integer RNG seed.  \code{"nndm"} is deterministic and
#'   reads it only to draw the \code{max_test} subset; \code{"buffered_loo"}
#'   likewise.  \code{"knndm"} uses it for the random partition among its
#'   candidates, for ties in the packing, and above 4,000 points for the
#'   k-means starts, so its folds are reproducible only with a seed.
#' @param block_nx,block_ny Optional grid dimensions for block_kfold, each a
#'   single whole number >= 1.  Give both, or give one and the other is
#'   derived from the extent's aspect ratio so that the blocks are roughly
#'   square.  Ignored when \code{block_size} or \code{auto_range} override
#'   them.
#' @param block_multiplier A single positive number, default 3.  When
#'   neither \code{block_size} nor \code{block_nx}/\code{block_ny} is given,
#'   the automatic grid aims for
#'   \code{block_multiplier * k} blocks over the extent (aspect-preserving),
#'   so each fold holds out about \code{block_multiplier} blocks.  An extent
#'   more than about \code{block_multiplier * k} times as wide as it is tall
#'   (or as tall as it is wide) gets a single row (or column) of that many
#'   blocks.  With 1,
#'   every fold is one contiguous region and the score depends heavily on
#'   which region each fold happened to get; with many, the blocks shrink
#'   towards single points and the scheme drifts back towards random k-fold.
#'   3 is a compromise between those two, not a published constant.  The
#'   block size that matters for leakage is the autocorrelation range, which
#'   is what \code{block_size} and \code{auto_range} control.
#' @param block_size Optional positive numeric minimum block edge length,
#'   \strong{in the units of the CRS the folds are built in}.  Takes
#'   precedence over \code{block_nx}/\code{block_ny} and
#'   \code{block_multiplier}.  When supplied, a whole number of blocks is
#'   fitted across each axis of the extent: as many as leave every block at
#'   least this long on that axis.  On an axis at least as long as
#'   \code{block_size} the blocks built are therefore from one to just under
#'   two times this size.  An axis shorter than \code{block_size} cannot
#'   hold a block that long: a single row or column spans it, and its blocks
#'   are as long as the extent on that axis, which is less than the size
#'   asked (a 39.8 km by 5 km extent with \code{block_size = 10000} gets
#'   three blocks 13.3 km wide and 5 km tall).  \code{params$block_width} and
#'   \code{params$block_height} record what was built.  A size that leaves
#'   fewer than \code{k} blocks holding points returns that many folds,
#'   with a warning (see \code{k}); one that leaves a single block is an
#'   error.  Both name the largest size that still gives \code{k} folds.
#'
#'   Which CRS that is depends on the input.  Projected input is used as it
#'   stands, so \code{block_size} is in your own CRS's units.  Geographic
#'   (lon/lat) input is projected first by \code{\link{ensure_projected}()},
#'   which picks a local UTM zone or, at wide extents, an equal-area
#'   projection.  That is a CRS you did not choose, whose units are metres but
#'   whose identity varies with the data.  \code{block_size} is then
#'   interpreted in \emph{that} CRS.  The CRS actually used is recorded in
#'   \code{params$crs} of the returned list; project the data yourself before
#'   calling if you want to fix the units in advance.
#'
#'   A \code{block_size} in the wrong unit asks for an enormous grid, so a
#'   request above 1,000,000 blocks is refused with an error naming the grid
#'   dimensions, the extent and the CRS's units.
#' @param auto_range Logical.  If \code{TRUE}, the spatial autocorrelation
#'   range is estimated via \code{estimate_sac_range()} (the omnidirectional,
#'   all-pairs range; the directional ranges it also fits are a diagnostic
#'   only, so on a field known to be anisotropic pass
#'   \code{block_size = max(attr(r, "directional_fitted"), na.rm = TRUE)}
#'   yourself after checking \code{attr(r, "directional_status")}, see
#'   \code{\link{estimate_sac_range}}) and used as the
#'   minimum \code{block_size}.  Requires \code{response_var}.  An explicit
#'   \code{block_size} takes precedence.  Default \code{FALSE}.  Sizing
#'   blocks from the autocorrelation range is the recommendation of Roberts
#'   et al. (2017) and what \pkg{blockCV} (Valavi et al. 2019) automates.
#'   \pkg{blockCV} takes the fitted variogram's range
#'   \emph{parameter} as the block size, whereas this uses the
#'   \emph{effective} range \code{estimate_sac_range()} returns (three
#'   times that parameter for an exponential fit), so its blocks are larger
#'   than \pkg{blockCV}'s from the same variogram.  When no range is
#'   identified, the geometric grid is used instead, with a warning that
#'   gives the reason.  Range-sized blocks do not fit every layout.  On a
#'   square extent at \code{k = 5} a range up to a third of the side gives
#'   five folds; up to half, a 2 x 2 grid and four folds, returned with the
#'   \code{"spatialkit_k_reduced"} warning; and above half the extent in
#'   both directions there is room for a single block, and
#'   \code{make_folds()} stops with an error naming the size, rather than
#'   return one fold with an empty training set.  The warning and the error
#'   name the largest block size that still gives \code{k} folds on this
#'   layout and how far below the range it is; \code{method = "nndm"} is the
#'   design that does not depend on a block size.  The message printed when
#'   the range is used states the grid and the block size built.
#' @param range_frac Passed through to \code{estimate_sac_range()} when
#'   \code{auto_range = TRUE}.  A fitted range beyond the longest lag the
#'   empirical variogram was fitted over is rejected as unidentified, and block
#'   sizing falls back to geometry (with a warning).  A range within that
#'   bound can still be too wide for two blocks; see \code{auto_range}.
#'   Default 1.0.
#' @param response_var Character(1) response column name.  Required when
#'   \code{auto_range = TRUE}.  For \code{"block_kfold"} a name that is not a
#'   column of \code{points_sf} is an error, whether or not
#'   \code{auto_range} is set: with it off the response still feeds the
#'   leakage warning, which a misspelt name would silently disable.
#' @param predictor_vars Optional character vector of predictor column names.
#'   Passed to \code{estimate_sac_range()} for residual variogram estimation.
#' @param boundary Optional polygonal sf/sfc for block_kfold.  The grid is
#'   clipped to it; a cell that the boundary only touches at a corner or
#'   along an edge is not a block.  A boundary without a CRS is brought into
#'   the points' CRS: reprojected from EPSG:4326 when its coordinates look
#'   like lon/lat, otherwise stamped, with a warning either way.  The same
#'   goes for \code{prediction_points} and \code{blocks}.
#' @param buffer For \code{"buffered_loo"}: a single positive number, the
#'   distance within which the held-out point's neighbours are excluded from
#'   its training set.  Like \code{block_size} it is in the units of the CRS
#'   the folds are built in (\code{params$crs}), which for geographic
#'   (lon/lat) input is the metre CRS \code{\link{ensure_projected}()}
#'   chooses: 0.1 means 0.1 m, not 0.1 degrees.  A buffer that excludes no
#'   neighbour from any fold makes the scheme plain leave-one-out, and is
#'   warned about.
#' @param group_var Character(1) naming a column of \code{points_sf} that
#'   identifies the location each observation belongs to.  Required for
#'   \code{method = "leave_location_out"}, which keeps every observation from a
#'   location together in the same fold.  Repeated measurements at the same
#'   site otherwise get split across folds, and the model is scored partly on
#'   sites it has already seen, which random k-fold reports as excellent
#'   performance.
#' @param prediction_points Optional \code{sf} layer of the locations you
#'   actually intend to predict onto.  Required for \code{method = "nndm"} and
#'   \code{method = "knndm"}.
#'   The grid from \code{\link{predict_surface}()} is the natural choice; a
#'   non-POINT layer (grid cells, polygons) is reduced to representative points
#'   first, so the target distances are point-to-point.  A point-to-polygon
#'   distance is zero for every cell that contains a training point, which
#'   pulls the target distribution towards zero and degenerates the CV towards
#'   plain leave-one-out.
#' @param min_train For \code{method = "nndm"}: the smallest fraction of the
#'   data any fold's training set may be reduced to by neighbour exclusion.
#'   Default \code{0.5}, as in \code{CAST::nndm()}.  Where it binds, the
#'   distance matching stops short and the cross-validation stays optimistic;
#'   see \strong{Details}.  For \code{method = "knndm"}: the smallest fraction
#'   of the data a fold may leave for training, so a candidate whose largest
#'   fold holds more than \code{1 - min_train} of the rows is set aside (the
#'   \code{maxp = 0.5} of \code{CAST::knndm()} is \code{min_train = 0.5}).
#'   Equal folds leave \code{1 - 1/k} of the data, so a value above that is
#'   an error.
#' @param phi For \code{method = "nndm"}: the distance up to which the two
#'   nearest-neighbour distance distributions are matched, in the CRS the
#'   folds are built in.  Matching is attempted only while a held-out
#'   point's nearest-neighbour distance is at most \code{phi}; the exclusion
#'   that takes it past \code{phi} is the last one, so a realised distance
#'   can exceed \code{phi} by up to one neighbour step, as in
#'   \code{CAST::nndm()}.  In Mila et al. (2022), and in
#'   \code{CAST::nndm()}, \eqn{\phi} is the autocorrelation range of the
#'   outcome: beyond it observations are effectively independent, so
#'   matching is unnecessary.  \code{\link{estimate_sac_range}()} gives such
#'   a value.  Default \code{NULL} = the largest prediction-to-training
#'   distance, which matches everywhere (\code{CAST}'s \code{phi = "max"}).
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
#' (2022): it sizes the exclusion around each held-out point, with no arbitrary
#' buffer, so that the resulting training-to-test distance distribution
#' approaches the distribution of distances from your actual prediction
#' locations to the training data.
#'
#' The procedure is the paper's, and it is deterministic.  Let \eqn{G_{ij}}
#' be the empirical distribution of
#' prediction-to-nearest-training distances and \eqn{G_j^*} the distribution
#' of each held-out point's nearest remaining training point.  Starting from
#' plain leave-one-out, the point with the smallest \eqn{G_j^*} at which the
#' realised distribution exceeds the target (\eqn{G_j^*(r) > G_{ij}(r)}) has
#' its nearest training neighbour removed, and this repeats until no such
#' point remains, subject to two limits: a point is matched only while its
#' nearest-neighbour distance is at most \code{phi} (default: the largest
#' prediction distance, since a training point already further than every
#' prediction distance has nothing to match), so the exclusion that takes it
#' past \code{phi} is its last and a realised distance can exceed \code{phi}
#' by up to one neighbour step; and no fold's training set is stripped below
#' \code{min_train} of the data.
#'
#' It differs from \code{CAST::nndm()} in two details, so the folds agree
#' closely with CAST's but not exactly.  The removal rule is strict: a
#' neighbour is removed whenever the realised distribution exceeds the target,
#' whereas CAST removes one only while the realised distribution, less the
#' point about to move, is still at or above the target.  The rule here
#' therefore removes up to one point more per distance value (a few percent
#' more removals in all on clustered layouts), erring on the pessimistic side.
#' And ties in \eqn{G_j^*} are broken by the points' coordinates, not by
#' their row index as in CAST, so the folds do not depend on the order of the
#' rows.
#'
#' Where neither limit binds, the realised distribution of
#' \emph{nearest-neighbour distances} is then never shorter than the target:
#' \eqn{G_j^*(r) \le G_{ij}(r)} up to the granularity of the neighbour
#' distances, which is the property the method exists to deliver.  It is a
#' statement about the distance to the one nearest training point.  A learner
#' that draws on several neighbours (a forest, a k-nearest-neighbour mean, a
#' kriging neighbourhood) also depends on how far the second and later
#' neighbours are, which NNDM does not match, so on clustered samples its
#' cross-validated error can still read low with the nearest-neighbour
#' distances matched exactly.  Beyond \code{phi} no matching is attempted, by
#' design.  \code{min_train} is different: when the prediction locations lie
#' further from the samples than a fold can be made to hold out (clustered
#' samples and a prediction domain well beyond them, the layout NNDM is meant
#' for), it stops the matching early and the realised distances stay
#' too short.  \code{params$n_at_min_train} counts the folds that were held
#' at the floor while still closer than the target allows, and
#' \code{make_folds()} warns when that leaves more than one point's worth of
#' excess at or below \code{phi}.  Lower \code{min_train} to match further, or
#' read the cross-validated score as an upper bound on performance at the
#' prediction locations.
#'
#' The match has a second side.  \code{params$max_ecdf_deficit} is the
#' largest amount by which the target distribution runs above the realised
#' one, the mirror of \code{params$max_ecdf_excess}: the share of prediction
#' points within some distance of the data that the folds do not reproduce.
#' Part of it is the granularity of the matching.  Excluding a neighbour
#' moves a held-out point out to its next remaining neighbour, and on
#' clustered samples that can be the next cluster, well past the prediction
#' distances.
#'
#' The other part no setting reaches.  Excluding neighbours can only lengthen
#' a held-out distance, so when the prediction points sit \emph{nearer} to
#' the data than the samples sit to each other (prediction points among the
#' samples, or a regular sample predicted between its nodes) the folds
#' returned are plain leave-one-out and are still farther from their
#' training data than the prediction points are.  The cross-validated error
#' is then pessimistic for those locations: an upper bound on the error
#' there.  \code{params$loo_ecdf_deficit} is the largest amount by which
#' the target distribution runs above that of plain leave-one-out.  When it
#' exceeds \code{params$ecdf_deficit_tol}, three times
#' \eqn{\sqrt{1/n + 1/m}} for \eqn{n} folds and \eqn{m} prediction points
#' (about three times the chance difference between two samples of one
#' distribution), \code{make_folds()} raises a warning of class
#' \code{"spatialkit_nndm_out_of_reach"} and sets
#' \code{params$target_out_of_reach} to \code{TRUE}.
#'
#' On small samples no deficit can exceed that tolerance.  A deficit is a
#' difference of two shares, so it is at most 1, and the tolerance is 1 or
#' more whenever \eqn{nm/(n+m) \le 9}: for nine folds or fewer whatever the
#' number of prediction points, and for up to eighteen folds with as many
#' prediction points.  There the same warning is raised in one case only,
#' complete separation: every prediction point nearer to the data than any
#' sample is to its nearest other sample (a deficit of exactly 1), provided
#' the two sets are large enough for that ordering to mean something, which
#' is taken as \code{choose(n + m, n) >= 1000} ways of interleaving them (17
#' prediction points for 3 folds, 7 for 6 folds, 4 for 12 folds).
#' \code{params$target_separated} records the ordering whether or not it
#' was enough to warn, and the warning's condition carries it as
#' \code{separated}.  With very few samples it can hold for a grid over the
#' samples' own extent, because a few points hug the box they span; the
#' folds are then farther than the grid all the same.  A partial shortfall
#' on such a sample (most prediction points nearer than leave-one-out
#' reaches, not all) is not warned about: the tolerance is 1 or more at
#' these sizes, which no difference of shares exceeds, so only complete
#' separation can be reported.  Compare
#' \code{params$target_median} with \code{params$realised_median}.
#'
#' \code{params$wasserstein} is the Wasserstein distance between the target
#' and realised distributions, in distance units, and
#' \code{params$wasserstein_loo} the same for plain leave-one-out: equal
#' when nothing was excluded, smaller when the matching helped.
#'
#' \strong{kNNDM.}  \code{"knndm"} is the k-fold version of that idea
#' (Linnenbrink et al. 2024): \code{k} folds, and so \code{k} model fits,
#' whose held-out distances are as close to the same target as a set of
#' candidates allows.  It has no row limit, and is the method to use above
#' the 5,000 rows \code{"nndm"} takes or when a model is too costly to refit
#' once per point.  It is not an equal of NNDM below that: a fold design has
#' far fewer degrees of freedom than one exclusion set per point, so the
#' match is coarser.
#'
#' The candidates are a random k-fold partition (the one
#' \code{"random_kfold"} builds with the same \code{seed}) and, for each of
#' up to 50 cluster counts spaced evenly on the log scale from \code{k} to
#' one fewer than the number of distinct locations, the samples clustered on
#' their coordinates and the clusters merged into \code{k} folds.  Few large
#' clusters give folds that are whole regions, far from their training data;
#' many small ones approach the random partition.  Up to 4,000 points the
#' clustering is Ward's (\code{stats::hclust(method = "ward.D2")}), cut at
#' each count.  Above that, where the n by n distance matrix it needs is no
#' longer affordable, k-means first reduces the samples to 4,000 groups
#' (Lloyd's algorithm from random sample locations, at most ten passes) and
#' Ward's method is continued from those groups; a count above 4,000 is a
#' k-means run of its own.  \code{params$clustering} says which
#' (\code{"ward"} or \code{"kmeans+ward"}).  Clusters are merged by the packing
#' \code{"block_kfold"} uses for its blocks, largest first and each to the
#' fold with the fewest points so far (see \strong{Fold balance}).
#' \code{CAST::knndm()} merges clusters by a rule of its own, which this was
#' \emph{not} checked against, so the folds are not expected to equal CAST's.
#'
#' Each candidate is scored by the Wasserstein distance \eqn{W} between its
#' held-out distances (each point's distance to the nearest point of another
#' fold) and the target distances, and the smallest \eqn{W} is returned.  A
#' candidate is set aside without being scored when its largest fold leaves
#' less than \code{min_train} of the rows for training or holds more than
#' \code{balance_tol} times the points of its smallest; the random partition
#' is always kept.  So \code{params$wasserstein} is never above
#' \code{params$wasserstein_random}, the \eqn{W} of that random partition,
#' and \code{params$candidates} lists every candidate with its cluster count,
#' fold sizes and \eqn{W}.  \code{params$chosen} is \code{"random"} when no
#' clustering came closer, which is the right answer for samples spread
#' evenly over the prediction area.
#'
#' Read \eqn{W} against the mean target distance.  On 400 points in 12
#' clusters predicted onto a grid over their extent (12 draws, 5 and 10
#' folds) the folds returned had \eqn{W} at 0.32 of the mean target distance,
#' against 0.70 to 0.72 for a random partition and 0.01 for \code{"nndm"}.
#' A 5-nearest-neighbour learner scored 9 to 11 percent high on them, where
#' random folds scored 23 to 25 percent low and \code{"nndm"} 6 percent low.
#' On 2,000 points in 40 clusters of unequal size the ratio was 0.07 to
#' 0.10, against 0.56 to 0.59 for a random partition and 0.003 to 0.006 for
#' \code{"nndm"}; on 30,000 points of that layout it was 0.11 to 0.15, and
#' on 100,000 points 0.27.  So the method removes most of the mismatch of a
#' random partition and leaves more of it than NNDM does, in either
#' direction: \eqn{W} has no sign.  On 400 evenly spread points the random
#' partition or one close to it came back and the two scored alike (3.3 and
#' 3.5 percent high).  With few folds a random partition is itself farther
#' than the target, since each fold trains on \code{1 - 1/k} of the points:
#' at \code{k = 2} on evenly spread points the mean held-out distance was 9
#' to 53 percent above the target mean (60 draws of 100 to 1,000 points),
#' and no candidate is nearer.
#'
#' Construction took 2 s at 8,000 points and 10 s at 30,000 with ten folds
#' (37 s at 100,000) with \pkg{FNN} installed, and three to four times that
#' without it; it grows with the number of folds.
#'
#' Two cases are warned about, each with a condition class.  When the
#' prediction points sit nearer to the data than a sample is to its nearest
#' other sample, no fold design reaches them, since a held-out point is never
#' nearer to its training data than that; the test is the one
#' \code{"nndm"} applies, the warning has the same class,
#' \code{"spatialkit_nndm_out_of_reach"}, and what comes back is the random
#' partition or a candidate as fine as one.  In the other direction, the
#' closest candidate can still be
#' short of the target: prediction points well beyond the samples lie farther
#' from the data than \code{k} folds of those samples can be held apart.
#' When the mean held-out distance of the folds returned is below the mean
#' target distance by more than a quarter of it and by at least half of
#' \eqn{W} (the test \code{\link{fold_separation}()} words as
#' \code{"closer"}), \code{make_folds()} raises a warning of class
#' \code{"spatialkit_knndm_short_of_target"} and sets
#' \code{params$short_of_target}: the cross-validated error then tends to
#' read low for those locations.  Fewer folds hold out larger regions.
#' Short of either warning the folds can still differ from the target at
#' some distances; \code{params$max_ecdf_excess} and
#' \code{params$max_ecdf_deficit} give the largest difference each way, and
#' \code{\link{fold_separation}()} reads the folds against the same
#' prediction points fold by fold.
#'
#' An earlier version of this package drew one random radius per point from
#' \eqn{G_{ij}} and excluded up to the order statistic \emph{closest} to it,
#' which rounds down half the time: on a two-cluster layout the realised
#' distribution exceeded the target by up to 0.17 (13% of folds had a
#' nearest training point within 50 m against a target of 9%), an
#' \emph{optimistic} cross-validation.  \code{params$max_ecdf_excess} reports
#' the largest remaining excess; compare \code{params$target_median} with
#' \code{params$realised_median} as well.
#' @references
#' Linnenbrink, J., Mila, C., Ludwig, M. and Meyer, H. (2024). kNNDM CV:
#' k-fold nearest-neighbour distance matching cross-validation for map
#' accuracy estimation. \emph{Geoscientific Model Development} \strong{17},
#' 5897-5912. \doi{10.5194/gmd-17-5897-2024}
#'
#' Mila, C., Mateu, J., Pebesma, E. and Meyer, H. (2022). Nearest neighbour
#' distance matching Leave-One-Out Cross-Validation for map validation.
#' \emph{Methods in Ecology and Evolution} \strong{13}, 1304-1316.
#' \doi{10.1111/2041-210X.13851}
#'
#' Roberts, D. R., Bahn, V., Ciuti, S., Boyce, M. S., Elith, J.,
#' Guillera-Arroita, G., Hauenstein, S., Lahoz-Monfort, J. J., Schroder, B.,
#' Thuiller, W., Warton, D. I., Wintle, B. A., Hartig, F. and Dormann, C. F.
#' (2017). Cross-validation strategies for data with temporal, spatial,
#' hierarchical, or phylogenetic structure. \emph{Ecography} \strong{40},
#' 913-929. \doi{10.1111/ecog.02881}
#'
#' Valavi, R., Elith, J., Lahoz-Monfort, J. J. and Guillera-Arroita, G. (2019).
#' blockCV: An R package for generating spatially or environmentally separated
#' folds for k-fold cross-validation of species distribution models.
#' \emph{Methods in Ecology and Evolution} \strong{10}, 225-232.
#' \doi{10.1111/2041-210X.13107}
#' @param drop_empty_blocks Logical. Default TRUE.  With \code{FALSE} the
#'   blocks that hold no point are kept and packed into folds too, but
#'   \code{k} is still lowered to the number of blocks that hold points, so
#'   no fold is left without test points.
#' @param blocks Optional polygon layer (\code{sf} or \code{sfc}, POLYGON or
#'   MULTIPOLYGON, at least two features) to use as the blocks of
#'   \code{"block_kfold"} in place of the grid this function would otherwise
#'   build: the \code{$cells} of a \code{\link{build_tessellation}()} result,
#'   hexagons, watersheds, administrative units, the \code{$blocks} of
#'   \code{blockCV::cv_spatial()}.  Each point takes the block that contains
#'   it, and the blocks are then assigned to folds exactly as grid cells are.
#'   \code{block_size}, \code{block_nx}/\code{block_ny},
#'   \code{block_multiplier} and \code{boundary} have nothing to act on and
#'   are ignored (logged); \code{auto_range} only compares the estimated range
#'   against the blocks.  See \strong{Supplied blocks} below for the CRS,
#'   overlap and coverage rules.  An error for any other \code{method}.
#' @param balance_tol A number of at least 1, default 3.  For \code{"block_kfold"}:
#'   the ratio of the largest fold's point count to the smallest's above
#'   which the folds are reported as imbalanced, with a warning (an R
#'   condition, also logged) that names both counts.  \code{Inf} disables the
#'   check.  The value is the tolerance of a check, not a target the packing
#'   aims for: see \strong{Fold balance} below for what the packing can and
#'   cannot do.  \code{params$balance_ratio} carries the ratio achieved.
#'   For \code{"knndm"} it is a limit and not a check: a candidate whose
#'   largest fold holds more than \code{balance_tol} times the points of its
#'   smallest is set aside, so the folds returned are always within it.
#' @param max_test For \code{"buffered_loo"} and \code{"nndm"}: optional
#'   whole number of at least 2.  Folds are built for a random subset of that
#'   many held-out points (drawn under \code{seed}), and every row stays
#'   available for training, so each fold returned is the fold the full
#'   design has for that point.  The model is then refitted \code{max_test}
#'   times and not \code{n} times, and the folds hold
#'   \code{4 * n * max_test} bytes of row indices in place of
#'   \code{4 * n^2}, which lets \code{"buffered_loo"} run on layers its
#'   storage limit otherwise refuses.  It does not lift the 5,000-row limit of
#'   \code{"nndm"}, whose cost is the matching itself; for that there is
#'   \code{"knndm"}.  Rows that are not held out have \code{NA} in
#'   \code{assignment$fold}.  For \code{"nndm"} the matching still runs over
#'   all \code{n} points and \code{params} describes all of them.
#'   \code{\link{area_of_applicability}()} needs every training row in a
#'   test fold and refuses a subset.  Default
#'   \code{NULL}: every point is held out once.  An error for the k-fold
#'   methods, whose number of folds is \code{k}.
#' @section Supplied blocks:
#' A polygon layer passed as \code{blocks} is aligned to the points the way
#' \code{boundary} is: CRS-less points are aligned to blocks that carry a CRS
#' (reprojected if they look like lon/lat, otherwise stamped, warning either
#' way), and the blocks are then brought into the CRS the folds are built in.
#' A point inside more than one block is given the first (lowest row) that
#' contains it, as for a point on the shared edge of two grid cells.  When the
#' blocks that caught such a point share area instead of an edge, the layer
#' overlaps and is not a partition, and this is warned about.  A point
#' inside no block is assigned to the nearest one, by distance to the polygon
#' itself, and the count of such points is warned about, unless they sit
#' within a millionth of the extent of a block, which is an edge that
#' reprojection or clipping moved by a rounding error.  Blocks that hold no
#' point are dropped when \code{drop_empty_blocks = TRUE}, and \code{k} is
#' lowered to the number of blocks that hold points when that is smaller.
#' \code{params$n_blocks} is the number of blocks before empties were dropped,
#' \code{params$blocks_used} the number after (so with
#' \code{drop_empty_blocks = FALSE} the two are equal, and
#' \code{sum(params$block_sizes > 0)} is how many of them hold points),
#' \code{params$grid_nx} and
#' \code{params$grid_ny} are \code{NA}, and \code{params$block_scale} is the
#' median over blocks that hold points of the side of the square with the
#' block's area.  That is the length compared against the autocorrelation range
#' for the leakage warning, since a polygon has no single edge length.
#'
#' The connection to the rest of the package is
#' \code{\link{build_tessellation}()}: every shape it builds can be a block
#' design here, including Voronoi cells around
#' \code{\link{get_voronoi_seeds}(method = "kmeans")} seeds, which adapt to the
#' density of the points.
#'
#' @section Fold balance:
#' \code{"block_kfold"} and \code{"knndm"} balance the number of points per
#' fold, and they do so by packing: blocks (for \code{"knndm"}, clusters of
#' the samples) are taken largest first and each goes to the
#' fold with the fewest points so far, ties broken at random.  This is the
#' longest-processing-time rule for multiway partitioning.  Measured against
#' the optimum by enumeration (two folds, up to ten blocks, heavy-tailed block
#' sizes) it is optimal in 72 percent of cases and within two points of
#' optimal on average; a local search with 30 random restarts moved the
#' largest-to-smallest ratio by 0.004 on average and never brought a packing
#' above the 3:1 tolerance below it.  An imbalance past the tolerance is
#' therefore in the points per block, which no assignment of blocks to folds
#' can even out, and this function offers no search over packings.  The
#' remedy is the block design: smaller blocks, or blocks that adapt to the
#' density of the points passed through \code{blocks}.  On clustered layouts
#' where the geometric grid exceeded 3:1 in 22 percent of cases (median ratio
#' 1.8, worst 6.3), Voronoi cells around 15 k-means seeds never exceeded 1.3
#' (median 1.09).  Density-adaptive blocks are smaller where points are
#' dense, so check \code{params$block_scale} against the autocorrelation
#' range as you would a grid.
#'
#' \code{"knndm"} chooses its clusters, so for it \code{balance_tol} is a
#' limit: a clustering whose packing exceeds it is not a candidate, and the
#' folds returned are within it (\code{params$fold_sizes},
#' \code{params$balance_ratio}).
#'
#' The other methods do not balance point counts.  \code{"random_kfold"} is
#' balanced by construction (fold sizes differ by at most one);
#' \code{"buffered_loo"} and \code{"nndm"} hold out one point per fold;
#' \code{"leave_location_out"} gives each fold the same number of
#' \emph{locations} (to within one), so folds differ by as much as the
#' locations' sizes do.
#'
#' @return A list with method, k, folds, assignment, params.  The
#'   \code{train}/\code{test} elements of each fold contain \code{..row_id}
#'   values (equal to row positions when the input has no pre-existing
#'   \code{..row_id} column), consistent with the \code{assignment} tibble.
#'   The returned \code{k} is the number of folds actually built, which is not
#'   always the \code{k} that was requested (see the \code{k} argument above),
#'   and \code{length(folds)} always matches it.
#'
#'   For \code{"block_kfold"} the block design is returned with the folds.
#'   \code{assignment} has a third column, \code{block_id}: the block each
#'   point fell in, numbered as the rows of \code{params$blocks}, an
#'   \code{sf} layer of the block polygons in the CRS the folds were built
#'   in.  Those rows are \strong{not} the rows of the layer the blocks came
#'   from: \code{drop_empty_blocks = TRUE} (the default) removes the blocks
#'   that hold no point and renumbers the rest, so with supplied
#'   \code{blocks} a nine-polygon layer of which three are empty comes back
#'   as six rows numbered 1 to 6.  \code{params$blocks$source_row} is the
#'   row each one came from in that original layer (a grid cell's index in
#'   the full \code{grid_nx} by \code{grid_ny} grid, or the row of the
#'   \code{blocks} argument), so
#'   \code{blocks[params$blocks$source_row, ]} recovers them with their own
#'   columns and in their own order.  It runs from 1 to
#'   \code{params$n_blocks} and is the identity when nothing was dropped.
#'   For a grid \code{params$n_blocks} is \code{grid_nx * grid_ny}, the
#'   cells a \code{boundary} clips away included, and
#'   \code{params$blocks_used} is the number of blocks returned.
#'   \code{params$block_sizes} is the number of points in each block, indexed
#'   by \code{block_id} (zeros are empty blocks that
#'   \code{drop_empty_blocks = FALSE} kept), and \code{params$fold_blocks}
#'   is a list with one integer vector per fold naming the blocks packed into
#'   it.  Between them the folds account for every block exactly once,
#'   empty ones included, so a fold's territory on the map is all of its
#'   blocks and not merely the ones that happen to hold points.  So
#'   \code{table(assignment$fold)} can be traced back to the blocks
#'   it is made of, a fold can be seen to be one contiguous region or
#'   several, and the blocks can be drawn over the data
#'   (\code{\link{plot_folds}()} does so).
#'
#'   For \code{"block_kfold"}, \code{params} carries a
#'   \code{params$blocks_supplied} that says whether the blocks came from
#'   \code{blocks} or from a grid built here, and
#'   \code{params$boundary_supplied} whether a \code{boundary} was given.
#'   \code{params$block_width} and \code{params$block_height} are the width
#'   and height of the grid cells built, before any clipping to a
#'   \code{boundary} (\code{NA} for supplied \code{blocks}); they are what
#'   \code{block_size} turned into: at least that size on an axis long
#'   enough to hold it, and the whole extent on an axis that is not (see
#'   \code{block_size}).
#'   Every method's \code{params} carries \code{params$k_requested}, the
#'   \code{k} that was asked for (\code{NA} when none was given), beside
#'   the \code{k} returned, and
#'   \code{params$row_probe}, a small sample of row IDs and coordinates
#'   that every \code{cv_*()} compares against the data it is handed, so
#'   folds built from a different layer of the same size are refused, never
#'   applied silently.
#'
#'   For \code{"knndm"}, \code{params} carries \code{chosen}
#'   (\code{"random"} or \code{"clusters"}), \code{n_clusters} (the number of
#'   clusters packed into the folds; \code{NA} for the random partition),
#'   \code{clustering} (\code{"ward"} or \code{"kmeans+ward"}; \code{NA} when
#'   no cluster count was tried), \code{candidates} (a data frame with one row
#'   per candidate, the random partition first: \code{candidate},
#'   \code{n_clusters}, \code{largest_fold}, \code{balance_ratio},
#'   \code{admissible} and \code{W}, which is \code{NA} for a candidate set
#'   aside), \code{wasserstein}, \code{wasserstein_random} and
#'   \code{wasserstein_loo} (the Wasserstein distance to the target of the
#'   folds returned, of the random partition among the candidates and of
#'   plain leave-one-out, in distance units), \code{fold_sizes},
#'   \code{balance_ratio}, \code{target_distances} and
#'   \code{realised_distances} (one per prediction point, and one per row:
#'   its distance to the nearest point of another fold) with their medians
#'   and means, and the diagnostics \code{"nndm"} reports under the same
#'   names: \code{max_ecdf_excess}, \code{max_ecdf_deficit},
#'   \code{loo_ecdf_deficit}, \code{ecdf_deficit_tol},
#'   \code{target_separated} and \code{target_out_of_reach}, with
#'   \code{short_of_target} for the opposite case (see \strong{Details}).
#'
#'   With \code{max_test}, \code{"buffered_loo"} and \code{"nndm"} return
#'   that many folds, \code{assignment$fold} is \code{NA} for the rows not
#'   held out, and \code{params$max_test} records the value.
#'
#'   For the methods that work in projected space, \code{params} also carries a
#'   \code{crs} element naming the CRS the folds were built in (an
#'   \code{"EPSG:code"} string where there is one, otherwise the CRS's input
#'   definition).  Every length in \code{params} (\code{block_size},
#'   \code{sac_range}, \code{buffer}, \code{median_buffer}) is in that CRS's
#'   units, which for geographic input is a CRS
#'   \code{\link{ensure_projected}()} chose for you, and not one you passed.
#'
#'   Rows whose geometry is empty or has non-finite coordinates are dropped
#'   before folding, with a logged warning naming the count; they appear in no
#'   fold and in no \code{assignment} row.
#' @family cross-validation
#' @examples
#' library(sf)
#' set.seed(1)
#' pts <- st_as_sf(
#'   data.frame(x = 5e5 + runif(30, 0, 1000), y = 5e6 + runif(30, 0, 1000)),
#'   coords = c("x", "y"), crs = 32632
#' )
#' folds <- make_folds(pts, k = 3, method = "block_kfold", seed = 42)
#' folds$assignment          # fold and block membership per row
#' lengths(folds$folds[[1]]) # train/test row-ID splits
#' folds$params$block_sizes  # points per block; params$fold_blocks packs them
#'
#' # Buffered leave-one-out: neighbours within 100 units excluded from training
#' loo <- make_folds(pts, k = 1, method = "buffered_loo", buffer = 100)
#'
#' # Three folds matched to the distances from a prediction grid to the data
#' grid <- st_make_grid(pts, n = c(10, 10), what = "centers")
#' km <- make_folds(pts, k = 3, method = "knndm", prediction_points = grid,
#'                  seed = 42)
#' km$params[c("chosen", "n_clusters", "wasserstein", "wasserstein_random")]
#' @export
make_folds <- function(points_sf, k,
                       method = c("random_kfold", "block_kfold", "buffered_loo",
                                  "leave_location_out", "nndm", "knndm"),
                       seed = NULL, block_nx = NULL, block_ny = NULL,
                       block_multiplier = 3, block_size = NULL,
                       auto_range = FALSE, range_frac = 1.0, response_var = NULL,
                       group_var = NULL, prediction_points = NULL,
                       predictor_vars = NULL, boundary = NULL,
                       buffer = NULL, min_train = 0.5, phi = NULL,
                       drop_empty_blocks = TRUE, blocks = NULL,
                       balance_tol = 3, max_test = NULL) {
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
  } else if (method %in% c("random_kfold", "block_kfold", "leave_location_out",
                           "knndm")) {
    # A missing or NULL k skipped the check above and died at the first
    # `if (k < 2)` with R's "argument is of length zero".  The two
    # leave-one-out methods never read k, so only the k-fold ones need it.
    stop(sprintf(paste0("make_folds(): `k` (the number of folds) is required ",
                        "for method = \"%s\"."), method), call. = FALSE)
  }
  # The k that was asked for travels with the k that is returned.  Three
  # methods can return fewer folds than asked, and a result whose `k` is 4
  # could not say whether 4 or 5 had been requested.
  k_requested <- if (!missing(k) && !is.null(k)) k else NA_integer_
  # block_size was tested with `is.numeric(block_size) && block_size > 0` and
  # anything failing that was silently ignored -- yet echoed back unchanged in
  # params$block_size, so a negative, zero or character value looked honoured.
  # NA and a length-2 vector reached the grid arithmetic and died as internal
  # R errors.  Validate it once, the way `k` is.  A `units` object passes
  # is.numeric() but fails `block_size <= 0` inside the units package with a
  # message that never names the argument, so it is refused here by name.
  if (!is.null(block_size) &&
      (inherits(block_size, "units") || !is.numeric(block_size) ||
       length(block_size) != 1L || !is.finite(block_size) || block_size <= 0))
    stop("make_folds(): `block_size` must be a single positive number in the ",
         "units of the data's CRS; got ",
         paste(format(block_size), collapse = ", "), ".", call. = FALSE)
  # block_nx/block_ny were never validated: 0, a negative, NA or a vector
  # reached st_make_grid() and failed there with sf's or base R's own
  # errors, and 2.7 was truncated to 2 without notice.
  for (nm in c("block_nx", "block_ny")) {
    v <- get(nm)
    if (!is.null(v) &&
        (inherits(v, "units") || !is.numeric(v) || length(v) != 1L ||
         !is.finite(v) || v != round(v) || v < 1))
      stop(sprintf("make_folds(): `%s` must be a single whole number >= 1; got %s.",
                   nm, paste(format(v), collapse = ", ")), call. = FALSE)
  }
  # block_multiplier was never validated: NA died inside st_make_grid() on
  # 'length.out', a vector silently used its largest element, and a `units`
  # value was silently taken as a plain count.
  if (inherits(block_multiplier, "units") || !is.numeric(block_multiplier) ||
      length(block_multiplier) != 1L || !is.finite(block_multiplier) ||
      block_multiplier <= 0)
    stop("make_folds(): `block_multiplier` must be a single positive number; got ",
         if (length(block_multiplier)) paste(format(block_multiplier), collapse = ", ")
         else "a value of length 0", ".", call. = FALSE)
  # A supplied block design is refused for the other methods rather than
  # ignored: hexagons that silently became random folds would be the worst
  # outcome.  The polygon check reuses .assert_sf(), which already recognises
  # a whole build_tessellation() result and says to pass its `$cells`.
  if (!is.null(blocks)) {
    if (method != "block_kfold")
      stop(sprintf(paste0("make_folds(): `blocks` is only used by method = ",
                          "\"block_kfold\"; got method = \"%s\"."), method),
           call. = FALSE)
    if (inherits(blocks, "sfc")) blocks <- sf::st_sf(geometry = blocks)
    .assert_sf(blocks, c("POLYGON", "MULTIPOLYGON"), "blocks", caller = "make_folds")
    if (nrow(blocks) < 2L)
      stop("make_folds(): `blocks` must hold at least 2 polygons; got ",
           nrow(blocks), ".", call. = FALSE)
  }
  if (!is.numeric(balance_tol) || length(balance_tol) != 1L || is.na(balance_tol) ||
      balance_tol < 1)
    stop("make_folds(): `balance_tol` must be a single number >= 1 (Inf disables ",
         "the imbalance warning); got ",
         paste(format(balance_tol), collapse = ", "), ".", call. = FALSE)
  # `max_test` belongs to the two leave-one-out methods.  It is refused for
  # the others, as `blocks` is: a k-fold design that quietly ignored it would
  # score every row while the caller believed a subset had been asked for.
  if (!is.null(max_test)) {
    if (!(method %in% c("buffered_loo", "nndm")))
      stop(sprintf(paste0(
        "make_folds(): `max_test` is only used by the leave-one-out methods, ",
        "\"buffered_loo\" and \"nndm\", where it limits how many held-out ",
        "points are built and scored; got method = \"%s\", whose number of ",
        "folds is `k`."), method), call. = FALSE)
    if (inherits(max_test, "units") || !is.numeric(max_test) ||
        length(max_test) != 1L || !is.finite(max_test) ||
        max_test != round(max_test) || max_test < 2 ||
        max_test > .Machine$integer.max)
      stop("make_folds(): `max_test` must be a single whole number >= 2 (the ",
           "number of held-out points to build folds for); got ",
           if (length(max_test)) paste(format(max_test), collapse = ", ")
           else "a value of length 0", ".", call. = FALSE)
    max_test <- as.integer(max_test)
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
         params = c(params, list(k_requested = k_requested,
                                 row_probe = row_probe)))

  # ---- RANDOM K-FOLD ----
  if (method == "random_kfold") {
    n <- nrow(points_sf)
    if (k < 2) k <- 2L
    if (k > n) k <- n
    # Fewer folds than asked is a different cross-validation from the one the
    # caller wrote down.  It was a log line only, which tryCatch(), knitr and
    # spatialkit_quiet() never show.
    if (!is.na(k_requested) && k < k_requested)
      .warn_k_reduced(k_requested, k, paste0(
        "make_folds(random_kfold): returning %d folds, not the k = %d asked ",
        "for: the layer has %d rows with usable coordinates and every fold ",
        "holds out at least one. This is leave-one-out; pass k <= %d for ",
        "fewer, larger folds."), k, k_requested, n, n)
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
    blocks_supplied <- !is.null(blocks)
    if (blocks_supplied) {
      # The blocks ARE the design: the region a grid would be built over and
      # the arguments that size a grid have nothing to act on.  Say so rather
      # than let a caller believe a block_size or boundary was honoured.
      if (!is.null(boundary))
        .log_warn("make_folds(block_kfold): `boundary` is ignored when `blocks` are supplied; the blocks define the region.")
      if (!is.null(block_size) || !is.null(block_nx) || !is.null(block_ny))
        .log_warn("make_folds(block_kfold): `block_size`, `block_nx` and `block_ny` are ignored when `blocks` are supplied; the blocks are used as given.")
      boundary <- NULL; block_size <- NULL; block_nx <- NULL; block_ny <- NULL
      # Same CRS alignment as `boundary` below: CRS-less points are aligned to
      # the blocks, and the blocks are then brought into the points' CRS.
      if (is.na(sf::st_crs(pts)) && !is.na(sf::st_crs(blocks)))
        pts <- .transform_or_stamp(pts, sf::st_crs(blocks),
                                   what = "points_sf", caller = "make_folds")
      if (!is.null(.crs_or_null(pts)))
        blocks <- .transform_or_stamp(blocks, sf::st_crs(pts),
                                      what = "blocks", caller = "make_folds")
      blocks <- .safe_make_valid(blocks)
    }
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
      # A CRS-less boundary is aligned to the points as every other function
      # given one aligns it: an R warning naming this function and argument.
      # ensure_projected()'s stamp was a log line naming neither, invisible
      # to tryCatch(), knitr and spatialkit_quiet().
      b <- if (is.null(.crs_or_null(pts))) ensure_projected(boundary)
           else .transform_or_stamp(boundary, sf::st_crs(pts),
                                    what = "boundary", caller = "make_folds")
      if (inherits(b, "sfc")) b <- sf::st_sf(geometry = b)
      b <- .safe_make_valid(sf::st_union(b))
      mat <- sf::st_intersects(pts, b, sparse = FALSE)
      inside_any <- apply(mat, 1L, any)
      if (!any(inside_any))
        # No point inside the boundary at all is almost certainly two layers
        # in different places -- a CRS that could only be stamped -- and
        # extending the region to cover them would hide it.
        stop("make_folds(block_kfold): none of the ", nrow(pts), " points fall ",
             "inside `boundary`. Check that the two layers cover the same ",
             "ground; a CRS that had to be stamped rather than reprojected is ",
             "the usual cause.", call. = FALSE)
      if (!all(inside_any)) {
        # Points outside the boundary were silently absorbed by extending
        # the region to the points' bounding box.  Say so: those points get
        # blocks, but the boundary the caller drew is not the one used.
        .warn_and_log(paste0("make_folds(block_kfold): %d of %d points fall ",
                             "outside `boundary`; the block region has been ",
                             "extended to the points' bounding box so every ",
                             "point receives a block."),
                      sum(!inside_any), nrow(pts))
        bb_pts <- sf::st_as_sfc(sf::st_bbox(pts)) |> sf::st_set_crs(sf::st_crs(pts))
        b <- suppressWarnings(
          sf::st_union(.safe_make_valid(sf::st_sf(geometry = c(b, bb_pts)))))
      }
      b
    } else {
      sf::st_as_sfc(sf::st_bbox(pts)) |> sf::st_set_crs(sf::st_crs(pts)) |> sf::st_sf()
    }

    # --- Autocorrelation-aware block sizing ---
    # A misspelt response_var was an error under auto_range = TRUE (from
    # estimate_sac_range()), but with auto_range off the diagnostic-only
    # estimate below swallowed it and the leakage check silently never ran.
    if (!is.null(response_var)) {
      if (!is.character(response_var) || length(response_var) != 1L ||
          is.na(response_var))
        stop("make_folds(block_kfold): `response_var` must be a single column ",
             "name; got ", paste(format(response_var), collapse = ", "), ".",
             call. = FALSE)
      if (!(response_var %in% names(points_sf)))
        stop(sprintf("make_folds(block_kfold): `response_var` '%s' is not a column of `points_sf`.",
                     response_var), call. = FALSE)
    }
    sac_range <- NA_real_
    size_from_range <- FALSE
    # The extent the grid is laid over.  Read before the range is, because
    # what a range-sized block yields depends on it.
    bb <- sf::st_bbox(reg)
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
      if (is.finite(sac_range) && sac_range > 0 && blocks_supplied) {
        # Supplied blocks cannot be resized; the range still feeds the
        # leakage diagnostic made on the blocks further down.
        message(sprintf(
          "make_folds(block_kfold): estimated spatial autocorrelation range = %.1f CRS units; the supplied `blocks` are used as given, so the range is only compared against them.",
          sac_range
        ))
      } else if (is.finite(sac_range) && sac_range > 0) {
        # auto_range sets block_size only if the caller didn't supply one
        if (is.null(block_size)) {
          block_size <- sac_range
          size_from_range <- TRUE
          # The range is a MINIMUM edge: a whole number of blocks is fitted
          # across each axis, so on an axis at least as long as the range a
          # block is 1 to under 2 times the range, and an axis shorter than
          # the range is spanned by one row or column, shorter than the
          # range.  Say what is built, since that, not the range, is the
          # design.
          built <- .block_dims_from_size(bb, block_size)
          message(sprintf(
            "make_folds(block_kfold): estimated spatial autocorrelation range = %.1f CRS units; using as minimum block size. Over the %s extent that is a %s x %s grid of blocks %s wide and %s tall.",
            sac_range, .block_extent_chr(bb),
            format(built$nx, scientific = FALSE), format(built$ny, scientific = FALSE),
            format(signif(as.numeric(bb["xmax"] - bb["xmin"]) / built$nx, 4)),
            format(signif(as.numeric(bb["ymax"] - bb["ymin"]) / built$ny, 4))
          ))
        } else {
          message(sprintf(
            "make_folds(block_kfold): estimated spatial autocorrelation range = %.1f CRS units; using as minimum block size.",
            sac_range
          ))
          if (block_size < sac_range)
            .warn_and_log(
              "make_folds(): block_size (%.1f) < estimated autocorrelation range (%.1f). Spatial CV may leak correlated information across folds. %s",
              block_size, sac_range,
              .block_range_advice(pts, reg, bb, sac_range, k, mention_auto = FALSE)
            )
        }
      } else {
        # The caller asked for range-sized blocks and is not getting them.
        # This was a log line only, the quietest of the three outcomes (the
        # success path messages, a missing response_var warns), so under
        # knitr, spatialkit_quiet or tryCatch() a CV result could not show
        # that its blocks were never sized from the data.
        why <- attr(sac_range, "rejected_reason")
        # estimate_sac_range() attaches a rejected_reason to every NA it
        # returns, including the ones it gives before fitting anything.  A
        # `sac_range` without one (hand-built, or saved by an older version)
        # still gets a reason: the two floors that are cheap to test here, and
        # the rest of the short list.
        if (!(is.character(why) && length(why) == 1L && !is.na(why)))
          why <- if (!requireNamespace("gstat", quietly = TRUE))
            "package 'gstat', which the variogram needs, is not installed"
          else if (nrow(pts) < 30L)
            sprintf("%d points, fewer than the 30 a variogram range is estimated from",
                    nrow(pts))
          else paste0("estimate_sac_range() returned NA before fitting: fewer ",
                      "than 30 finite values to model, a variable with no ",
                      "variance, or points with no extent; its log line says ",
                      "which")
        .warn_and_log(paste0("make_folds(block_kfold): auto_range = TRUE, but ",
                             "no autocorrelation range was identified (%s); ",
                             "falling back to geometric blocks, which are not ",
                             "sized from the data. See ?estimate_sac_range, or ",
                             "pass `block_size`."),
                      why)
      }
    } else if (isTRUE(auto_range) && is.null(response_var)) {
      .warn_and_log("make_folds(): auto_range requires response_var; ignoring.")
    } else if (!is.null(response_var)) {
      # auto_range is off, but a response is to hand -- which it always is when
      # a cv_*() function built the folds -- so estimate the range for the
      # LEAKAGE DIAGNOSTIC alone.  It sizes nothing: the blocks below are the
      # same geometric blocks as before, and the folds do not change.  Before
      # this branch existed, `sac_range` stayed NA on every default call, so
      # the two "block smaller than the autocorrelation range" warnings below
      # could fire only in the one configuration (auto_range = TRUE) that had
      # already sized the blocks from the range and therefore never needed
      # them.  The estimate's own log lines are silenced on the console
      # (index 2): the caller did not ask for a range, and a variogram that
      # never reached a sill simply means no diagnostic can be given here.
      # The file trace (index 1) keeps them.
      sac_range <- .sac_range_for_diagnostic(
        pts, response_var, predictor_vars, range_frac,
        seed = if (is.null(seed)) 123L else seed)
    }

    if (!blocks_supplied) {
      # Determine grid dimensions: block_size constrains nx/ny
      if (!is.null(block_size) && is.numeric(block_size) && block_size > 0) {
        size_dims <- .block_dims_from_size(bb, block_size)
        # The auto_range branch above already compared a supplied block_size to
        # the range it estimated.  When the range came from the diagnostic-only
        # branch instead, make the same comparison here: a hand-set block_size
        # below the range leaks exactly as a geometric one does.
        if (!isTRUE(auto_range) && is.finite(sac_range) && sac_range > 0 &&
            block_size < sac_range) {
          .warn_and_log(
            "make_folds(): block_size (%.1f) < estimated autocorrelation range (%.1f). Spatial CV may leak correlated information across folds. %s",
            block_size, sac_range,
            .block_range_advice(pts, reg, bb, sac_range, k, mention_auto = FALSE)
          )
        }
        # If the caller also supplied explicit block_nx/block_ny, warn about override
        if (!is.null(block_nx) || !is.null(block_ny)) {
          .log_warn(
            "make_folds(block_kfold): block_size (%.1f) overrides explicit block_nx/block_ny.",
            block_size
          )
        }
        nx <- size_dims$nx
        ny <- size_dims$ny

        # A grid of fewer than k cells lowers k.  That is decided once, further
        # down, from the blocks that hold points: a cell count announced here
        # could be lowered again there (a 3 x 3 grid with four occupied cells
        # is four folds, not nine), and a single block is an error.
      } else if (is.null(block_nx) && is.null(block_ny)) {
        w  <- as.numeric(bb["xmax"] - bb["xmin"])
        h  <- as.numeric(bb["ymax"] - bb["ymin"])
        # Points on one horizontal line have h == 0.  Treating that as a
        # square (ratio 1) gave a 4 x 4 grid whose rows all collapse onto the
        # line, so only 4 blocks existed and k = 5 was lowered to 4.
        ratio <- if (h > 0) w / h else if (w > 0) Inf else 1
        target_blocks <- max(1L, round(block_multiplier * k))
        # nx is capped at the target just as ny is floored at 1.  Without the
        # cap the rule was not symmetric: a tall extent got 1 x 15 blocks,
        # but the same layer turned on its side got round(sqrt(15 * w/h))
        # columns -- 39 x 1 on a 10 km x 100 m corridor -- blocks less than
        # half as long, and a scheme drifting towards random k-fold.
        nx <- min(target_blocks, max(1L, round(sqrt(target_blocks * ratio))))
        ny <- max(1L, round(max(1, target_blocks / nx)))

        # Diagnostic: warn if resulting block size is small relative to SAC range.
        # Only a dimension that is split has blocks bordering each other
        # across it: a single row spans the whole height and leaks to no
        # other block that way.  Taking the minimum over both compared the
        # range with 0 on points on a line (and with the width of a
        # corridor), and the advice to pass block_size = <range> then gave
        # SHORTER blocks along it.
        split_cells <- c(if (nx > 1) w / nx, if (ny > 1) h / ny)
        if (is.finite(sac_range) && sac_range > 0 && length(split_cells)) {
          min_cell <- min(split_cells)
          if (min_cell < sac_range) {
            # The advice used to be "Pass block_size = <range> or auto_range
            # = TRUE" whatever the layout.  On a square extent at k = 5 that
            # gives five folds only for a range between a quarter and a third
            # of the side: up to a half it returned four folds, and above it
            # stopped with the single-block error, whose own advice led back
            # here.  .block_range_advice() builds the grid the advice asks
            # for and says what comes back.
            .warn_and_log(
              "make_folds(): block dimension (%.1f) < autocorrelation range (%.1f). Spatial CV may leak correlated information across folds. %s",
              min_cell, sac_range,
              .block_range_advice(pts, reg, bb, sac_range, k)
            )
          }
        }
      } else {
        w  <- as.numeric(bb["xmax"] - bb["xmin"])
        h  <- as.numeric(bb["ymax"] - bb["ymin"])
        # Giving only one of the two used to send the call to the automatic
        # grid, silently: block_nx = 10 alone gave a 3 x 4 grid.  Honour the
        # one given and derive the other so the blocks are roughly square.
        if (is.null(block_ny)) {
          nx <- as.integer(block_nx)
          ny <- if (w > 0) max(1, round(nx * h / w)) else 1
          .log_info("make_folds(block_kfold): only block_nx given; block_ny = %s derived from the extent's aspect ratio.", format(ny, scientific = FALSE))
        } else if (is.null(block_nx)) {
          ny <- as.integer(block_ny)
          nx <- if (h > 0) max(1, round(ny * w / h)) else 1
          .log_info("make_folds(block_kfold): only block_ny given; block_nx = %s derived from the extent's aspect ratio.", format(nx, scientific = FALSE))
        } else {
          nx <- as.integer(block_nx); ny <- as.integer(block_ny)
        }

        # Diagnostic: warn if user-supplied nx/ny yield blocks smaller than SAC
        # (over the dimensions that are split; see the automatic grid above).
        split_cells <- c(if (nx > 1) w / nx, if (ny > 1) h / ny)
        if (is.finite(sac_range) && sac_range > 0 && length(split_cells)) {
          min_cell <- min(split_cells)
          if (min_cell < sac_range) {
            .warn_and_log(
              "make_folds(): block_nx/block_ny yield blocks smaller than autocorrelation range (%.1f). Spatial CV may leak correlated information across folds. %s",
              sac_range,
              .block_range_advice(pts, reg, bb, sac_range, k)
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
      # A block_size hundreds of orders of magnitude too small overflows the
      # product to Inf, which the old `is.finite() &&` let past the guard and
      # into st_make_grid()'s "result would be too long a vector".
      if (!is.finite(n_cells_est) || n_cells_est > .block_max_cells) {
        unit_lbl <- tryCatch({
          u <- sf::st_crs(reg)$units_gdal
          if (is.null(u) || is.na(u) || !nzchar(u)) "CRS units" else u
        }, error = function(e) "CRS units")
        # %s, not %d: nx and ny are doubles from floor(), and a block_size in
        # the wrong unit -- the very mistake this guard exists to explain --
        # gives counts past 2^31 that %d refuses with "invalid format".
        fmt_count <- function(v, big = "") {
          if (!is.finite(v)) "more than 1e308"
          else if (v >= 1e15) format(v, digits = 3, scientific = TRUE)
          else format(v, big.mark = big, scientific = FALSE)
        }
        extent_lbl <- sprintf("%s x %s",
                              format(signif(as.numeric(bb["xmax"] - bb["xmin"]), 4)),
                              format(signif(as.numeric(bb["ymax"] - bb["ymin"]), 4)))
        # Name what produced the grid.  This always blamed `block_size`,
        # printing "unset" when block_nx/block_ny or the automatic grid had
        # sized it -- an argument the caller never passed.
        check <- if (!is.null(block_size) && size_from_range)
          sprintf(paste0("The block size is the estimated autocorrelation range ",
                         "(%s %s) over an extent of %s; pass `block_size` ",
                         "yourself, or set auto_range = FALSE."),
                  format(as.numeric(block_size)), unit_lbl, extent_lbl)
        else if (!is.null(block_size))
          sprintf(paste0("Check that `block_size` (%s) is expressed in the ",
                         "data's CRS units (%s) over an extent of %s; a value ",
                         "in the wrong unit is the usual cause."),
                  format(as.numeric(block_size)), unit_lbl, extent_lbl)
        else if (!is.null(block_nx) || !is.null(block_ny))
          sprintf(paste0("The grid is the one `block_nx`/`block_ny` ask for over ",
                         "an extent of %s (%s); ask for fewer blocks."),
                  extent_lbl, unit_lbl)
        else
          sprintf(paste0("The automatic grid aims for `block_multiplier` (%s) ",
                         "x k (%d) blocks over an extent of %s (%s); pass a ",
                         "smaller `block_multiplier`."),
                  format(block_multiplier), k, extent_lbl, unit_lbl)
        stop(sprintf(paste0("make_folds(block_kfold): the requested grid is %s x ",
                            "%s = %s cells, above the %s this function will ",
                            "build. %s"),
                    fmt_count(nx), fmt_count(ny), fmt_count(n_cells_est, ","),
                    format(.block_max_cells, big.mark = ",", scientific = FALSE),
                    check),
             call. = FALSE)
      }

      # .block_grid() lays the cells and clips them to the region.  The
      # leakage advice and the k warning build their trial grids with the same
      # function, so what they say a block size yields is what it yields.
      built_grid <- .block_grid(reg, nx, ny)
      grid_sf  <- built_grid$grid_sf
      cell_idx <- built_grid$cell_idx
      n_blocks <- as.integer(nx * ny)
    } else {
      # The caller's polygons are the blocks.  Their row order is the block
      # id, so `assignment` can be joined back to the layer that was passed.
      grid_sf <- sf::st_as_sf(blocks)
      nx <- NA_integer_; ny <- NA_integer_
      cell_idx <- seq_len(nrow(grid_sf))
      n_blocks <- nrow(grid_sf)
    }
    hits <- sf::st_intersects(pts, grid_sf)
    block_id <- vapply(hits, function(ix) if (length(ix)) ix[1] else NA_integer_, 1L)
    pts$..block_id <- block_id

    if (blocks_supplied) {
      # A point in two grid cells sits on their shared edge and either cell
      # will do.  A point in two SUPPLIED blocks may mean the same -- adjacent
      # polygons share edges too -- or that the design overlaps and is not a
      # partition at all.  Only the pairs that actually caught a point are
      # tested, by whether they share area rather than a line.
      multi <- which(lengths(hits) > 1L)
      if (length(multi)) {
        pairs <- unique(t(vapply(unclass(hits)[multi],
                                 function(ix) sort(ix[1:2]), integer(2))))
        g <- sf::st_geometry(grid_sf)
        overlaps <- vapply(seq_len(nrow(pairs)), function(r) {
          a <- suppressWarnings(sf::st_intersection(g[pairs[r, 1L]], g[pairs[r, 2L]]))
          length(a) > 0L && sum(as.numeric(sf::st_area(a))) > 0
        }, logical(1))
        if (any(overlaps))
          .warn_and_log(paste0("make_folds(block_kfold): %d point(s) fall inside ",
                               "more than one of the supplied `blocks`, which ",
                               "overlap; each is assigned to the first block ",
                               "(lowest row) that contains it."),
                        length(multi))
      }
      na_idx <- which(is.na(pts$..block_id))
      if (length(na_idx) == nrow(pts))
        stop("make_folds(block_kfold): none of the ", nrow(pts), " points fall ",
             "inside any of the supplied `blocks`. Check that the two layers ",
             "cover the same ground; a CRS that had to be stamped rather than ",
             "reprojected is the usual cause.", call. = FALSE)
      if (length(na_idx)) {
        # A point outside every block goes to the nearest one, by distance to
        # the polygon itself (a supplied polygon need not be compact, so its
        # centroid is no guide).  Points within a millionth of the extent of
        # a block are on an edge that reprojection or clipping moved by a
        # rounding error, not outside the design; the rest are reported.
        dmat <- as.matrix(sf::st_distance(sf::st_geometry(pts[na_idx, ]),
                                          sf::st_geometry(grid_sf)))
        nearest <- vapply(seq_len(nrow(dmat)), function(i) {
          w <- which.min(dmat[i, ])
          if (length(w)) as.integer(w) else NA_integer_
        }, integer(1))
        d_min <- as.numeric(dmat[cbind(seq_along(nearest), nearest)])
        tol <- 1e-6 * max(as.numeric(bb["xmax"] - bb["xmin"]),
                          as.numeric(bb["ymax"] - bb["ymin"]), 0)
        n_far <- sum(is.finite(d_min) & d_min > tol)
        if (n_far > 0L)
          .warn_and_log(paste0("make_folds(block_kfold): %d of %d points fall ",
                               "outside every supplied block (the farthest by ",
                               "%.1f units); each has been assigned to the ",
                               "nearest block."),
                        n_far, nrow(pts), max(d_min[is.finite(d_min)]))
        pts$..block_id[na_idx] <- nearest
      }
    }

    # Which row of the ORIGINAL block layer each surviving block came from.
    # Dropping the empty ones renumbers `..block_id`, so without this the
    # returned design cannot be tied back to the layer the caller supplied --
    # and a join by row position silently mis-attributes every block after
    # the first gap.
    block_source_row <- cell_idx
    if (drop_empty_blocks) {
      used_blocks <- sort(unique(pts$..block_id[!is.na(pts$..block_id)]))
      grid_sf <- grid_sf[used_blocks, , drop = FALSE]
      pts$..block_id <- match(pts$..block_id, used_blocks)
      block_source_row <- cell_idx[used_blocks]
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
    # The number of blocks that hold a point, which is what k is limited by.
    # This was the highest block id a point fell in: the same number once
    # empty blocks are dropped and renumbered, but under
    # drop_empty_blocks = FALSE only an artifact of the numbering.  With two
    # clusters on a 4 x 4 grid it was 16, so k = 5 was kept for 2 occupied
    # blocks and three folds came back with no test points at all -- no
    # warning, because the balance check skips an empty fold.
    B <- length(unique(pts$..block_id[!is.na(pts$..block_id)]))
    # One block means one fold whose training set is empty -- blocked CV
    # silently degenerating into nothing at all.  It happens whenever the block
    # size exceeds half the extent, which an accepted autocorrelation range can
    # do (estimate_sac_range() only rejects ranges above half the DIAGONAL, a
    # factor of sqrt(2) looser than what .block_dims_from_size() needs).
    if (B < 2L) {
      # block_size is NULL when the caller set block_nx/block_ny directly, and
      # sprintf() on a zero-length argument yields character(0) -- an empty
      # error message.
      # Three ways to arrive here; name the one the caller actually used, not
      # an argument they never passed.
      if (blocks_supplied)
        stop(sprintf(paste0("make_folds(block_kfold): all %d points fall in the ",
                            "same one of the supplied `blocks`, so there is no ",
                            "spatial split to make and the one fold would have ",
                            "an empty training set. Supply smaller blocks."),
                     nrow(pts)), call. = FALSE)
      how <- if (size_from_range)
        sprintf("the block size (%s, the estimated autocorrelation range)",
                format(as.numeric(block_size)))
      else if (!is.null(block_size))
        sprintf("the block size (%s)", format(as.numeric(block_size)))
      else if (!is.null(block_nx) || !is.null(block_ny))
        "the requested block_nx/block_ny"
      else
        "the automatic grid (block_multiplier x k blocks over the extent)"
      # What to do instead, worked out for this layout.  The old text said
      # "set auto_range = FALSE if the size came from the estimated
      # autocorrelation range", which returns the default grid and its leakage
      # warning, and that warning said to set auto_range = TRUE.
      alt <- if (!is.null(block_size))
        .block_largest_size_for_k(pts, reg, bb, k, below = block_size) else NULL
      todo <- if (size_from_range)
        paste0("Block folds cannot be separated by this range here. ",
               if (is.null(alt))
                 sprintf("No block size down to 1/200 of the extent gives %d blocks that hold points. ", k)
               else sprintf(paste0("The largest block_size that gives %d ",
                                   "folds is %s, %.0f%% of the range. "),
                            k, alt$chr, 100 * alt$size / as.numeric(sac_range)),
               .block_other_method(nrow(pts)))
      else if (!is.null(block_size)) {
        if (is.null(alt))
          sprintf("No block size down to 1/200 of the extent gives %d blocks that hold points; pass a smaller `k`, or blocks that follow the points (`blocks`).", k)
        else sprintf("Pass a smaller `block_size`: the largest that gives %d folds is %s.",
                     k, alt$chr)
      } else if (!is.null(block_nx) || !is.null(block_ny))
        "Ask for more blocks."
      else
        "Every point lies in one block of the grid; block folds need points in at least two."
      stop(sprintf(paste0("make_folds(block_kfold): %s produces a single block ",
                          "covering the whole extent (%s), so there is no ",
                          "spatial split to make and the one fold would have ",
                          "an empty training set. %s"),
                   how, .block_extent_chr(bb), todo),
           call. = FALSE)
    }
    if (B < k) {
      # Raised as a condition with its own class, not only logged: a result
      # with fewer folds than asked is a different cross-validation, and
      # under knitr, tryCatch() or spatialkit_quiet() the log line never
      # reached the caller.  Raised once, here, from the blocks that hold
      # points, whichever of the grid's cell count or its empty cells is the
      # cause.
      k_asked <- k
      k <- B
      why <- if (blocks_supplied)
        sprintf("%d of the %d supplied blocks hold points", B, n_blocks)
      else if (size_from_range)
        sprintf(paste0("the estimated autocorrelation range (%s), used as the ",
                       "block size, builds a %s x %s grid over the %s extent ",
                       "and %d of its blocks hold points"),
                format(as.numeric(block_size)), format(nx, scientific = FALSE),
                format(ny, scientific = FALSE), .block_extent_chr(bb), B)
      else if (!is.null(block_size))
        sprintf(paste0("block_size = %s builds a %s x %s grid over the %s ",
                       "extent and %d of its blocks hold points"),
                format(as.numeric(block_size)), format(nx, scientific = FALSE),
                format(ny, scientific = FALSE), .block_extent_chr(bb), B)
      else if (!is.null(block_nx) || !is.null(block_ny))
        sprintf("the requested %s x %s grid has %d blocks that hold points",
                format(nx, scientific = FALSE), format(ny, scientific = FALSE), B)
      else
        sprintf("the automatic %s x %s grid has %d blocks that hold points",
                format(nx, scientific = FALSE), format(ny, scientific = FALSE), B)
      alt <- if (!blocks_supplied && !is.null(block_size))
        .block_largest_size_for_k(pts, reg, bb, k_asked, below = block_size) else NULL
      todo <- if (blocks_supplied)
        sprintf("Supply blocks that divide the points further, or pass k = %d.", B)
      else if (size_from_range)
        paste0(if (is.null(alt)) ""
               else sprintf(paste0("The largest block_size that gives %d folds ",
                                   "is %s, %.0f%% of the range, so its folds ",
                                   "are separated by less than the range. "),
                            k_asked, alt$chr, 100 * alt$size / as.numeric(sac_range)),
               sprintf("Pass k = %d to keep range-sized blocks. ", B),
               .block_other_method(nrow(pts)))
      else if (!is.null(block_size)) {
        if (is.null(alt)) sprintf("Pass k = %d, or blocks that follow the points (`blocks`).", B)
        else sprintf(paste0("For %d folds pass block_size = %s (the largest ",
                            "that gives %d), or pass k = %d to keep these blocks."),
                     k_asked, alt$chr, k_asked, B)
      } else if (!is.null(block_nx) || !is.null(block_ny))
        sprintf("Ask for a finer grid, or pass k = %d.", B)
      else
        sprintf("Pass a larger `block_multiplier`, blocks that follow the points (`blocks`), or k = %d.", B)
      .warn_k_reduced(k_asked, k, paste0(
        "make_folds(block_kfold): returning %d folds, not the k = %d asked ",
        "for: %s, and every fold needs a block of its own. %s"),
        k, k_asked, why, todo)
    }

    # Counted over every block in the grid, not over 1..B.  B used to be the
    # highest block id a point fell in, which under drop_empty_blocks = FALSE
    # says nothing about how many blocks there are: on a layer whose points
    # sit in one quadrant of an 8x8 grid, it was 27, so blocks 28-64 were
    # packed into no fold at all while the 12 equally empty blocks below 27
    # were -- the
    # same kind of block treated two ways depending on where it fell in an
    # arbitrary numbering.  Counting over the grid makes `fold_blocks` cover
    # every block `blocks` and `block_sizes` describe.  It cannot move a
    # point: order() is stable, so the populated blocks are still placed
    # first, in the same order, against the same fold loads; the extra
    # zero-size blocks land at the tail and add no observations to any fold.
    blk_sizes <- as.integer(table(factor(pts$..block_id,
                                         levels = seq_len(nrow(grid_sf)))))

    # Leakage diagnostic for supplied blocks, the counterpart of the grid
    # branches above: the scale of a polygon is the side of the square with
    # its area, taken over the blocks that hold points.
    block_scale <- NA_real_
    if (blocks_supplied) {
      block_scale <- stats::median(
        sqrt(as.numeric(sf::st_area(sf::st_geometry(grid_sf)[which(blk_sizes > 0L)]))),
        na.rm = TRUE)
      if (is.finite(sac_range) && sac_range > 0 && is.finite(block_scale) &&
          block_scale < sac_range) {
        .warn_and_log(
          "make_folds(): the supplied blocks (median scale %.1f) are smaller than the autocorrelation range (%.1f). Spatial CV may leak correlated information across folds; supply blocks at least %.0f units across.",
          block_scale, sac_range, ceiling(sac_range)
        )
      }
    }

    # Largest block first, each to the fold with the fewest points so far:
    # the longest-processing-time rule for multiway partitioning.  Measured
    # against the optimum by enumeration (k = 2, up to 10 blocks, heavy-tailed
    # sizes) it is optimal in 72% of cases and within 2 points of optimal on
    # average, and a local search plus 30 random restarts moved the max/min
    # ratio by 0.004 on average and never brought an over-tolerance packing
    # under 3:1 -- imbalance past the tolerance lives in the block sizes, not
    # in the packing, so no search is offered.  The remedy is `blocks`.
    # The rule is .pack_largest_first(), which method = "knndm" uses for its
    # clusters too.  Each fold's blocks are listed in the order they were
    # placed.
    packed     <- .pack_largest_first(blk_sizes, k)
    fold_loads <- packed$loads
    fold_blocks <- lapply(seq_len(k), function(j)
      packed$order[packed$fold_of[packed$order] == j])

    # The residual imbalance is checked against the tolerance and, past it,
    # raised as a warning a pipeline can catch -- not only logged.  No fold
    # can be empty here: k is at most the number of blocks holding points,
    # and the packing gives each of the k largest blocks a fold of its own
    # (every fold is at load 0 until it has one), so the Inf branch is only
    # a guard against dividing by zero.
    balance_ratio <- if (min(fold_loads) > 0L) max(fold_loads) / min(fold_loads) else Inf
    if (min(fold_loads) > 0L && balance_ratio > balance_tol) {
      .warn_and_log(paste0("make_folds(block_kfold): fold size imbalance -- ",
                           "largest fold has %d obs vs %d in smallest (ratio ",
                           "%.2f, tolerance %s). The imbalance is in the ",
                           "points per block, which no assignment of blocks ",
                           "to folds can even out; use smaller blocks, or ",
                           "supply count-adaptive ones through `blocks`."),
                    max(fold_loads), min(fold_loads), balance_ratio,
                    format(balance_tol))
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
    # The block design itself travels with the folds: which block each point
    # fell in (a third `assignment` column), the points per block indexed by
    # block_id -- over every block still in the grid, so empties are counted
    # when drop_empty_blocks = FALSE kept them -- which blocks each fold was
    # packed from, and the block polygons in the CRS the folds were built in,
    # so the design can be drawn over the data and a fold seen to be one
    # region or several.
    block_polys <- sf::st_sf(block_id = seq_len(nrow(grid_sf)),
                             source_row = as.integer(block_source_row),
                             geometry = sf::st_geometry(grid_sf))
    return(.ret(method, k, splits,
                .safe_tibble(row_id = pts$..row_id, fold = assign_vec,
                             block_id = as.integer(pts$..block_id)),
                # The blocks the design actually has, which is what `blocks`
                # and `block_sizes` describe and what the folds partition.
                # This was `B`, the highest block id a point fell in -- equal
                # to the block count under the default, but an artifact of the
                # numbering under drop_empty_blocks = FALSE, where it reported
                # 27 for a 64-block grid that had dropped nothing.
                list(seed = seed, grid_nx = nx, grid_ny = ny,
                     blocks_used = nrow(grid_sf),
                     n_blocks = n_blocks,
                     blocks_supplied = blocks_supplied,
                     block_scale = block_scale,
                     block_sizes = blk_sizes,
                     fold_blocks = fold_blocks,
                     blocks = block_polys,
                     block_multiplier = block_multiplier,
                     block_size = block_size,
                     # What was built.  `block_size` is a minimum edge and a
                     # whole number of cells is fitted across each axis, so
                     # on an axis at least as long as `block_size` a cell is
                     # 1 to under 2 times as long (asked 4000, built
                     # 4981.6).  An axis shorter than `block_size` is spanned
                     # by one row or column, and the cell is then as long as
                     # the extent there: shorter than asked (asked 10000 on
                     # a 39840 x 4969 extent, built 13280 x 4969).  NA for
                     # supplied blocks, which have no single width; see
                     # `block_scale`.
                     block_width  = if (blocks_supplied) NA_real_ else
                       as.numeric(bb["xmax"] - bb["xmin"]) / nx,
                     block_height = if (blocks_supplied) NA_real_ else
                       as.numeric(bb["ymax"] - bb["ymin"]) / ny,
                     sac_range = sac_range,
                     auto_range = auto_range,
                     balance_ratio = balance_ratio,
                     balance_tol = balance_tol,
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
    # NA, numeric(0) and a vector used to reach `buffer <= 0` and die with
    # base R's "missing value where TRUE/FALSE needed" or "length = 2 in
    # coercion" -- and an NA is exactly what estimate_sac_range() returns
    # when it identifies no range.  A `units` object failed inside the units
    # package without naming `buffer`.  Inf is let through: the check below
    # reports that it spans the data, which is the more useful message.
    if (is.null(buffer) || inherits(buffer, "units") || !is.numeric(buffer) ||
        length(buffer) != 1L || is.na(buffer) || buffer <= 0)
      stop(sprintf(paste0(
        "make_folds(buffered_loo): `buffer` must be a single positive number, ",
        "in the units of the CRS the folds are built in (metres for lon/lat ",
        "input, which is projected first); got %s.%s"),
        if (is.null(buffer)) "NULL"
        else if (!length(buffer)) "a value of length 0"
        else paste(format(buffer), collapse = ", "),
        if (length(buffer) == 1L && is.na(buffer))
          " An NA from estimate_sac_range() means no range was identified; choose a buffer yourself."
        else ""), call. = FALSE)
    pts <- ensure_projected(points_sf)
    n <- nrow(pts)
    # With `max_test`, folds are built for a seeded random subset of the
    # held-out points only.  Every row stays available for training, so the
    # folds built are the ones the full design would have for those points.
    scored <- .loo_scored_rows(n, max_test)
    n_s    <- length(scored)
    # The cost is QUADRATIC in n whatever the buffer: every one of the n
    # splits stores its own training-row vector of length about n, so the
    # fold object alone is ~4 n^2 bytes.  The old threshold of 20000 admitted
    # 1.6 GB of splits (measured 1.7 GB of R objects at n = 19999), before
    # .remap_folds() copied them.  The cap is now stated in bytes and the
    # message says what the request would cost.  `max_test` splits cost
    # 4 * n * max_test bytes, which is what makes the method usable past the
    # cap.
    split_gb <- 4 * as.numeric(n) * as.numeric(n_s) / 1024^3
    if (split_gb > 0.6) {
      stop(sprintf(paste0(
        "make_folds(buffered_loo): n = %d would produce %d leave-one-out ",
        "splits holding about %.2f GB of row indices (the training vector of ",
        "every split is ~n long, so storage grows as n^2), before the model is ",
        "refitted %d times. Use 'block_kfold' instead, or subset the data.%s"),
        n, n_s, split_gb, n_s,
        # The largest subset that stays under the cap, as a round number.
        sprintf(paste0(" Or pass `max_test` to build and score a random ",
                       "subset of the held-out points while every row stays ",
                       "available for training: max_test = %s keeps the ",
                       "splits under 0.6 GB here."),
                format(max(2, floor(.signif_down(
                  0.6 * 1024^3 / (4 * as.numeric(n)), 2L))),
                  scientific = FALSE, trim = TRUE))), call. = FALSE)
    } else if (split_gb > 0.1) {
      .log_warn(
        "make_folds(buffered_loo): n = %d produces %d splits holding about %.2f GB of row indices, and the model will be refitted %d times; this may be slow.",
        n, n_s, split_gb, n_s
      )
    }
    nb <- sf::st_is_within_distance(sf::st_geometry(pts)[scored],
                                     sf::st_geometry(pts), dist = buffer)
    splits <- vector("list", n_s)
    row_ids <- pts$..row_id
    for (i in seq_len(n_s)) {
      # The held-out point is within any positive buffer of itself, so it is
      # among the excluded rows.
      excl <- sort(unique(c(scored[i], nb[[i]])))
      # Row-ID (not positional) splits — see comment in the random_kfold branch.
      splits[[i]] <- list(train = row_ids[setdiff(seq_len(n), excl)],
                          test  = row_ids[scored[i]])
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
    # The opposite failure: a buffer below every nearest-neighbour distance
    # excludes nothing, and buffered LOO is then plain LOO -- the optimistic
    # scheme it exists to replace.  The usual cause is a buffer in degrees on
    # lon/lat input (0.1 read as 0.1 m once the data are projected), which
    # used to pass without a word.  A co-located duplicate is caught by any
    # positive buffer, so every nb[[i]] of length 1 means nothing was dropped.
    if (all(lengths(nb) <= 1L)) {
      crs_lbl <- .fold_crs_label(pts)
      .warn_and_log(paste0("make_folds(buffered_loo): a buffer of %s excludes ",
                           "no neighbour from any fold, so this is plain ",
                           "leave-one-out. `buffer` is in the units of %s (the ",
                           "CRS the folds are built in; metres for lon/lat ",
                           "input), and every point is further than that from ",
                           "its nearest neighbour."),
                    format(buffer),
                    if (is.na(crs_lbl)) "the data's CRS" else crs_lbl)
    }

    # A row that is not held out has no fold: NA in `assignment`.
    prm <- list(buffer = buffer, crs = .fold_crs_label(pts))
    if (!is.null(max_test)) prm$max_test <- max_test
    return(.ret(method, n_s, splits,
                .safe_tibble(row_id = pts$..row_id,
                             fold = match(seq_len(n), scored)),
                prm))
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
    if (k > n_g) k <- n_g
    if (k < 2L) k <- 2L
    if (!is.na(k_requested) && k < k_requested)
      .warn_k_reduced(k_requested, k, paste0(
        "make_folds(leave_location_out): returning %d folds, not the k = %d ",
        "asked for: `%s` has %d distinct groups and a group is never split ",
        "across folds. This is leave-one-group-out; pass k <= %d to hold out ",
        "several groups per fold."), k, k_requested, group_var, n_g, n_g)

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
      stop(sprintf(paste0(
        "make_folds(nndm): n = %d exceeds the safety threshold of 5000. Fold ",
        "construction keeps a sorted table of up to n/2 neighbours per point ",
        "(memory grows as n^2), and where min_train binds for most points it ",
        "makes about n^2/2 removals at O(n) each, so the worst case is O(n^3) ",
        "time (about nine minutes at n = 3000). NNDM then produces n ",
        "leave-one-out folds, so the model is refitted n times. Use ",
        "method = \"knndm\", which matches k folds to the same prediction ",
        "points and has no row limit, or 'block_kfold', or subset your data."),
        n), call. = FALSE)
    if (n < 3L)
      stop("make_folds(nndm): need at least 3 points.", call. = FALSE)

    # The prediction points, aligned to the samples, and the distance from
    # each to its nearest sample: the target of the matching.
    tg <- .nndm_target(pts, prediction_points, "nndm")
    pts <- tg$pts; pred <- tg$pred; g_target <- tg$g_target

    # A `units` object passes is.numeric() and then fails the comparison
    # inside the units package with a message that names no argument.
    if (inherits(min_train, "units") || !is.numeric(min_train) ||
        length(min_train) != 1L || !is.finite(min_train) || min_train <= 0 ||
        min_train >= 1)
      stop("make_folds(nndm): `min_train` must be a single number in (0, 1).",
           call. = FALSE)
    if (is.null(phi)) phi <- max(g_target)
    if (inherits(phi, "units") || !is.numeric(phi) || length(phi) != 1L ||
        !is.finite(phi) || phi < 0)
      stop("make_folds(nndm): `phi` must be a single non-negative number, in ",
           "the units of the CRS the folds are built in (a plain number, not a ",
           "units object).", call. = FALSE)

    # ---- The paper's procedure (Mila et al. 2022) ----------------------------
    # Deterministic: no radii are drawn.  Starting from plain LOO, the held-out
    # point with the SMALLEST current nearest-neighbour distance at which the
    # realised distribution exceeds the target has its nearest training
    # neighbour removed, and this repeats until no such point remains.
    # CAST::nndm() tests (cnt - 1)/n >= G instead of cnt/n > G, and breaks
    # ties by row index; the strict test here removes up to one point more
    # per distance value, the pessimistic side (see @details).
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
      # Release get.knn()'s own copy now.  While it is referenced, the
      # repair below and the column subset after it work on duplicates, and
      # the originals -- two n x n/2 tables, a few hundred MB at the 5000
      # cap -- stayed alive beside them through the sweep and the splits.
      kn <- NULL
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
    at_floor <- logical(n)
    Gjstar  <- nn_d[, 1L]
    # Plain leave-one-out's distances, kept for scale: the Wasserstein
    # distance of the folds to the target is reported beside this one's.
    loo_d   <- Gjstar
    # Ties in Gjstar -- every mutual-nearest-neighbour pair, all of a regular
    # grid -- are broken by a key that depends on the GEOMETRY, not on the
    # row order: the point's rank in (x, y) lexicographic order.  CAST breaks
    # them by row index (which.min), so shuffling the rows of one layer gave
    # 6-18 different folds out of 300 and moved 72 of 100 points on a 10 x 10
    # grid.  Identical data must give identical folds.
    tie_key <- order(order(xy[, 1L], xy[, 2L]))
    o  <- order(Gjstar, tie_key); sv <- Gjstar[o]; si <- o
    k  <- 1L
    n_iter <- 0L
    while (k <= n) {
      r   <- sv[k]
      j   <- si[k]
      cnt <- findInterval(r, sv)                    # realised count <= r
      over <- is.finite(r) && (cnt / n) > G_target(r) + 1e-12 && r <= phi
      violates <- over && (n - 1L - removed[j]) > rmin &&
                  removed[j] + 1L < ncol(nn_d)
      if (violates) {
        removed[j] <- removed[j] + 1L
        newv <- nn_d[j, removed[j] + 1L]
        sv <- sv[-k]; si <- si[-k]
        # Insert after every smaller value and, among equal values, after
        # those belonging to points with a smaller tie_key, so that on
        # clustered data (where pushed points pile up at the same
        # cluster-to-cluster distance) the SAME points get pushed whatever
        # the row order.
        pos <- findInterval(newv, sv)
        while (pos > 0L && sv[pos] == newv && tie_key[si[pos]] > tie_key[j])
          pos <- pos - 1L
        sv <- append(sv, newv, after = pos)
        si <- append(si, j,    after = pos)
        n_iter <- n_iter + 1L
      } else {
        # Still more optimistic than the target here, but min_train forbids
        # stripping this fold's training set any further: the matching stops
        # short for this point.  Counted so it can be reported (see below).
        if (over) at_floor[j] <- TRUE
        k <- k + 1L
      }
    }
    Gjstar[si] <- sv
    nn_d <- NULL                     # only the neighbour ids are needed now

    # The matching above runs over every point, with or without `max_test`:
    # what a fold excludes depends on all the others.  `max_test` only limits
    # which of the n folds are built and returned, a seeded random subset, so
    # the folds returned are the ones the full design has for those points
    # and every diagnostic below still describes all n.
    scored     <- .loo_scored_rows(n, max_test)
    splits     <- vector("list", length(scored))
    n_excluded <- removed
    realised   <- Gjstar
    for (s in seq_along(scored)) {
      i <- scored[s]
      drop_idx <- if (removed[i] > 0L) nn_i[i, seq_len(removed[i])] else integer(0)
      drop_idx <- drop_idx[!is.na(drop_idx)]
      train_i  <- setdiff(seq_len(n), c(i, drop_idx))
      splits[[s]] <- list(train = row_ids[train_i], test = row_ids[i])
    }
    # The paper's own diagnostic: largest excess of realised over target.
    fin <- is.finite(realised)
    max_excess <- if (any(fin)) {
      rs <- sort(realised[fin])
      max(findInterval(rs, rs) / length(rs) - G_target(rs))
    } else NA_real_

    # min_train binds where the prediction locations lie further from the
    # samples than a fold can be made to hold out -- clustered samples and a
    # prediction grid well beyond them, the layout NNDM is meant for.  The
    # realised distances then stay optimistic (one cluster predicted onto a
    # 20 km grid: median 1171 m against a target of 6704 m, 96 of 100 folds
    # at the floor), and the only signal was an INFO line in a log file.
    # Warn when the floor left more than one point's worth of excess at or
    # below phi; beyond phi no matching is attempted by design, so an excess
    # there is not the floor's doing.
    n_at_floor <- sum(at_floor)
    excess_phi <- if (any(fin & realised <= phi)) {
      rs <- sort(realised[fin])
      rp <- rs[rs <= phi]
      max(findInterval(rp, rs) / length(rs) - G_target(rp))
    } else NA_real_
    # The message states the statistic that fired it.  It used to lead with
    # the two medians, which describe the whole distributions, not the short
    # distances the test is about: the realised median could sit above the
    # target's while the warning said the folds were "more optimistic".
    if (n_at_floor > 0L && is.finite(excess_phi) && excess_phi > 1 / n + 1e-9)
      .warn_and_log(paste0(
        "make_folds(nndm): min_train = %s stopped the distance matching in %d ",
        "of %d folds, which keep a training point closer than the prediction ",
        "distances allow. At distances up to phi (%.4g) the share of folds ",
        "whose nearest training point is within a distance runs above the ",
        "target's share by as much as %.3f (more than one fold's worth, 1/%d), ",
        "so short distances are over-represented and the CV score is ",
        "optimistic for the prediction locations. (The medians, %.4g realised ",
        "and %.4g target, summarise all distances and need not differ the same ",
        "way.) The prediction locations lie further from the samples than a ",
        "fold can be made to hold out: read the CV score as an upper bound on ",
        "performance there, or lower min_train."),
        format(min_train), n_at_floor, n, phi, excess_phi, n,
        stats::median(realised[fin]), stats::median(g_target))

    # The mirror of the excess: the largest amount by which the TARGET
    # distribution runs above the realised one, the share of prediction points
    # within some distance of the data that the folds do not reproduce.  It
    # has two sources.  Excluding a neighbour moves a held-out point out to
    # its next remaining neighbour, which on clustered samples can be far
    # beyond the prediction distances: that part is the granularity of the
    # matching.  The other part was there before anything was excluded.
    max_deficit <- .ecdf_deficit(realised, g_sorted)
    # Exclusion can only lengthen a held-out distance, so where the target
    # runs above plain leave-one-out's own distribution no setting can reach
    # it.  When prediction locations sit about as far from the data as the
    # samples sit from each other, plain leave-one-out already matches and
    # nothing is removed, which is the right answer.  When they sit NEARER
    # (prediction points among the samples, a regular sample predicted
    # between its nodes), the folds returned are plain leave-one-out and
    # still too far, and that was an INFO line only: prediction points 5 m
    # from 200 samples gave a target median of 6.6 against a realised 355.7,
    # with no removal, an excess of 0 and no warning.  The rule, and what it
    # was calibrated on, is .nndm_reach(); method = "knndm" applies the same
    # one.
    reach        <- .nndm_reach(loo_d, g_target)
    loo_deficit  <- reach$deficit
    deficit_tol  <- reach$tolerance
    past_tol     <- reach$past_tol
    loo_fin      <- reach$loo_fin
    separated    <- reach$separated
    out_of_reach <- reach$out_of_reach
    # The Wasserstein distance between the two distributions (the area
    # between their ECDFs, in distance units), with plain leave-one-out's for
    # scale: equal when nothing was excluded, smaller when matching helped.
    w_real <- if (any(fin)) .fold_sep_compare(realised[fin], g_target)$w1 else NA_real_
    w_loo  <- if (length(loo_fin)) .fold_sep_compare(loo_fin, g_target)$w1 else NA_real_
    if (out_of_reach)
      .warn_nndm_out_of_reach(loo_deficit, deficit_tol, separated, paste0(
        "make_folds(nndm): the prediction points sit nearer to the data than ",
        "leave-one-out can reach. Their distance to the nearest sample has ",
        "median %.4g; a sample's nearest other sample is a median of %.4g ",
        "away, and the folds hold each point out at a median of %.4g from ",
        "its training data (%s). %s ",
        "Excluding neighbours can only lengthen a held-out distance, so no ",
        "setting closes this gap: the cross-validated error is pessimistic ",
        "for these prediction locations, an upper bound on the error there. ",
        "If the prediction points were meant to cover a wider area than the ",
        "samples, check that `prediction_points` is that layer."),
        stats::median(g_target), stats::median(loo_fin),
        stats::median(realised[fin]),
        if (n_iter == 0L) "plain leave-one-out, no neighbour excluded"
        else sprintf("%d neighbour exclusion(s) in all", n_iter),
        # The evidence, by the rule that fired.
        if (past_tol) sprintf(paste0(
          "The share of prediction points within a distance of the data runs ",
          "above leave-one-out's share by as much as %.3f, past the %.3f ",
          "allowed for %d folds and %d prediction points."),
          loo_deficit, deficit_tol, n, n_g)
        else sprintf(paste0(
          "Every one of the %d prediction points is nearer to the data (%s ",
          "at most) than any sample is to its nearest other sample (%s at ",
          "least). With %d folds and %d prediction points the allowance for a ",
          "difference between the two shares is %.3f, more than a difference ",
          "of shares can be, so the warning rests on that complete separation."),
          n_g, format(signif(max(g_target), 4)), format(signif(min(loo_fin), 4)),
          n, n_g, deficit_tol))

    .log_info(paste0("make_folds(nndm): %d folds%s. Target prediction-to-training ",
                     "distance: median %.1f. Realised test-to-training distance: ",
                     "median %.1f, excluding a median of %.1f training point(s) ",
                     "per fold (%d removals in total; largest remaining excess ",
                     "of realised over target ECDF %.3f, largest deficit %.3f)."),
              n,
              if (length(scored) < n)
                sprintf(" matched, %d of them built (max_test)", length(scored))
              else "",
              stats::median(g_target), stats::median(realised[fin]),
              stats::median(n_excluded), n_iter, max_excess, max_deficit)

    return(.ret(method, length(scored), splits,
                # A row that is not held out has no fold: NA.
                .safe_tibble(row_id = row_ids, fold = match(seq_len(n), scored)),
                c(list(seed = seed, n_prediction_points = nrow(pred),
                     crs = .fold_crs_label(pts),
                     min_train = min_train, phi = phi,
                     # The effective buffer each fold ended up with.
                     median_buffer   = stats::median(realised[fin]),
                     median_excluded = stats::median(n_excluded),
                     n_removed_total = n_iter,
                     max_ecdf_excess = max_excess,
                     # The other side of the match: how far the target runs
                     # above the folds' distribution, and above plain
                     # leave-one-out's, which no exclusion can reach.
                     # target_out_of_reach is TRUE when the warning above
                     # was raised: the latter is past its tolerance, or the
                     # two sets of distances are completely separated
                     # (target_separated) and large enough to say so.
                     max_ecdf_deficit    = max_deficit,
                     loo_ecdf_deficit    = loo_deficit,
                     ecdf_deficit_tol    = deficit_tol,
                     target_separated    = separated,
                     target_out_of_reach = out_of_reach,
                     wasserstein         = w_real,
                     wasserstein_loo     = w_loo,
                     n_at_min_train  = n_at_floor,
                     target_median   = stats::median(g_target),
                     realised_median = stats::median(realised[fin]),
                     target_distances   = g_target,
                     realised_distances = realised),
                  # Present only when a subset of the folds was asked for.
                  if (!is.null(max_test)) list(max_test = max_test))))
  }

  # ---- kNNDM: k-fold nearest-neighbour distance matching ----
  if (method == "knndm") {
    pts <- ensure_projected(points_sf)
    n <- nrow(pts)
    if (n < 3L)
      stop("make_folds(knndm): need at least 3 points.", call. = FALSE)

    # The same target as method = "nndm": the distance from each prediction
    # point to its nearest sample.
    tg <- .nndm_target(pts, prediction_points, "knndm")
    pts <- tg$pts; pred <- tg$pred; g_target <- tg$g_target

    if (inherits(min_train, "units") || !is.numeric(min_train) ||
        length(min_train) != 1L || !is.finite(min_train) || min_train <= 0 ||
        min_train >= 1)
      stop("make_folds(knndm): `min_train` must be a single number in (0, 1).",
           call. = FALSE)

    if (k < 2L) k <- 2L
    if (k > n) k <- n
    if (!is.na(k_requested) && k < k_requested)
      .warn_k_reduced(k_requested, k, paste0(
        "make_folds(knndm): returning %d folds, not the k = %d asked for: ",
        "the layer has %d rows with usable coordinates and every fold holds ",
        "out at least one. This is leave-one-out; pass k <= %d for fewer, ",
        "larger folds."), k, k_requested, n, n)
    # k folds hold 1/k of the rows each when they are equal, so no k-fold
    # design leaves a training set above 1 - 1/k of the data.  A min_train
    # past that asks for something no candidate can give.
    if ((1 - min_train) < 1 / k - 1e-9)
      stop(sprintf(paste0(
        "make_folds(knndm): min_train = %s cannot be met with k = %d folds: ",
        "equal folds hold 1/%d of the rows each, which leaves a training set ",
        "of %s of the data at most. Use min_train <= %s, or more folds."),
        format(min_train), k, k, format(signif(1 - 1 / k, 3)),
        format(.signif_down(1 - 1 / k, 3L))), call. = FALSE)

    xy      <- sf::st_coordinates(pts)[, 1:2, drop = FALSE]
    row_ids <- pts$..row_id
    res <- .knndm_folds(xy, g_target, k, min_train = min_train,
                        balance_tol = balance_tol)
    lab      <- res$fold
    realised <- res$realised
    fin      <- is.finite(realised)
    cand     <- res$candidates
    chosen   <- cand[res$chosen, , drop = FALSE]
    is_random <- identical(chosen$candidate, "random")
    fold_sizes <- tabulate(lab, nbins = k)

    splits <- lapply(seq_len(k), function(j)
      # Row-ID (not positional) splits -- see the random_kfold branch.
      list(train = row_ids[lab != j], test = row_ids[lab == j]))

    # Both sides of the match, as method = "nndm" reports them.
    g_sorted    <- sort(g_target)
    max_excess  <- if (any(fin)) .ecdf_deficit(g_target, sort(realised[fin]))
                   else NA_real_
    max_deficit <- .ecdf_deficit(realised, g_sorted)
    # No fold design holds a point out nearer than its nearest other sample,
    # so a target that plain leave-one-out cannot reach is out of reach of
    # every candidate here too.  The rule is the one "nndm" applies, on the
    # same two sets of distances.
    loo_d  <- res$loo_distances
    reach  <- .nndm_reach(loo_d, g_target)
    w_loo  <- if (length(reach$loo_fin))
      .fold_sep_compare(reach$loo_fin, g_target)$w1 else NA_real_
    n_g    <- length(g_target)
    design <- if (is_random) "the random partition"
              else sprintf("%d clusters of the samples packed into %d folds",
                           chosen$n_clusters, k)
    if (reach$out_of_reach)
      .warn_nndm_out_of_reach(reach$deficit, reach$tolerance, reach$separated,
        paste0(
        "make_folds(knndm): the prediction points sit nearer to the data than ",
        "any fold can hold a point out. Their distance to the nearest sample ",
        "has median %.4g; a sample's nearest other sample is a median of %.4g ",
        "away, and the %d folds returned (%s) hold each point out at a median ",
        "of %.4g from its training data. %s ",
        "A held-out point is never nearer to its training data than to its ",
        "nearest other sample, so no number of folds and no clustering ",
        "closes this gap: the cross-validated error is pessimistic for these ",
        "prediction locations, an upper bound on the error there. If the ",
        "prediction points were meant to cover a wider area than the samples, ",
        "check that `prediction_points` is that layer."),
        stats::median(g_target), stats::median(reach$loo_fin), k, design,
        stats::median(realised[fin]),
        # The evidence, by the rule that fired.
        if (reach$past_tol) sprintf(paste0(
          "The share of prediction points within a distance of the data runs ",
          "above the share of samples within that distance of another sample ",
          "by as much as %.3f, past the %.3f allowed for %d samples and %d ",
          "prediction points."),
          reach$deficit, reach$tolerance, n, n_g)
        else sprintf(paste0(
          "Every one of the %d prediction points is nearer to the data (%s ",
          "at most) than any sample is to its nearest other sample (%s at ",
          "least). With %d samples and %d prediction points the allowance for ",
          "a difference between the two shares is %.3f, which no ",
          "difference of shares can exceed, so the warning rests on that complete ",
          "separation."),
          n_g, format(signif(max(g_target), 4)),
          format(signif(min(reach$loo_fin), 4)), n, n_g, reach$tolerance))

    n_set_aside <- sum(cand$candidate == "clusters" & !cand$admissible)
    # The other direction is the optimistic one.  The candidate kept is the
    # closest of those tried, which is not the same as close: prediction
    # points well beyond the samples lie farther from the data than k folds
    # of these samples can be held apart (12 clusters over a 10 km box,
    # predicted onto a 30 km one: held-out median 1427 against a target of
    # 7220, and nothing said).  The test is the one fold_separation() words
    # as "closer" on the distance scale: the mean held-out distance below the
    # target's by more than a quarter of it, and by at least half of W, so
    # that the gap is a shift to one side and not two distributions crossing.
    # Over 20 draws each of uniform, 12-cluster and 4-cluster samples of 100,
    # 300 and 1,000 points at 2, 5 and 10 folds: with the prediction grid
    # over the samples' own extent it fired in 0 of 180 uniform layouts, 4 of
    # 180 with 12 clusters and 52 of 180 with 4 clusters (two or ten folds
    # cannot follow four clusters, and the grid fills the gaps between them);
    # with the grid over a box 1.5, 2 and 3 times as wide, in 315, 529 and
    # 540 of 540.  The out-of-reach warning above fired in none of those
    # 2,160 and in 540 of 540 with prediction points 5 m from the samples.
    mean_real   <- mean(realised[fin])
    mean_target <- mean(g_target)
    short <- identical(.fold_sep_reading(mean_real - mean_target, chosen$W,
                                         .fold_sep_tol_dist * mean_target),
                       "closer")
    if (short) {
      is_cl <- cand$candidate == "clusters"
      coarsest <- !is_random &&
        chosen$n_clusters == min(cand$n_clusters[is_cl & cand$admissible])
      # The clusterings with fewer clusters than the one kept, all of them
      # set aside when it is the coarsest allowed.
      n_coarser <- if (is_random) 0L else
        sum(is_cl & !cand$admissible & cand$n_clusters < chosen$n_clusters)
      .warn_knndm_short(mean_real, mean_target, paste0(
        "make_folds(knndm): the %d folds returned (%s) hold each point out ",
        "at a mean of %.4g from its training data, and the prediction points ",
        "sit a mean of %.4g from the data (medians %.4g and %.4g). This is ",
        "the closest of the %d candidates scored (Wasserstein distance %.4g; ",
        "the random partition %.4g) and it is still short of the target, so ",
        "the cross-validated error tends to read low for these prediction ",
        "locations: read the score as an upper bound on performance there.%s"),
        k, design, mean_real, mean_target, stats::median(realised[fin]),
        stats::median(g_target), sum(cand$admissible), chosen$W, res$w_random,
        if (coarsest && n_coarser > 0L) sprintf(paste0(
          " It is the coarsest clustering allowed: the %d with fewer clusters ",
          "left a fold with more than %s of the rows (min_train = %s) or more ",
          "than %s times the points of another (balance_tol). Fewer folds ",
          "hold out larger regions."),
          n_coarser, format(signif(1 - min_train, 3)), format(min_train),
          format(balance_tol))
        else if (coarsest) sprintf(paste0(
          " It is the coarsest clustering there is, one cluster per fold, so ",
          "these samples cannot be held further apart in %d folds. Fewer ",
          "folds hold out larger regions."), k)
        else "")
    }
    .log_info(paste0("make_folds(knndm): %d folds from %s. Target ",
                     "prediction-to-data distance: median %.1f. Held-out-to-",
                     "training distance: median %.1f. Wasserstein distance to ",
                     "the target %.4g (random partition %.4g, leave-one-out ",
                     "%.4g). %d cluster count(s) tried by %s clustering, %d ",
                     "set aside by min_train or balance_tol."),
              k, design, stats::median(g_target), stats::median(realised[fin]),
              chosen$W, res$w_random, w_loo,
              sum(cand$candidate == "clusters"),
              if (is.na(res$clustering)) "no" else res$clustering, n_set_aside)

    return(.ret(method, k, splits,
                .safe_tibble(row_id = row_ids, fold = lab),
                list(seed = seed, n_prediction_points = nrow(pred),
                     crs = .fold_crs_label(pts),
                     min_train = min_train,
                     # What was chosen: "random", or the clustering and the
                     # number of clusters packed into the k folds.
                     clustering   = res$clustering,
                     chosen       = chosen$candidate,
                     n_clusters   = chosen$n_clusters,
                     # Every candidate tried, the random partition first.
                     candidates   = cand,
                     # The distance between the folds' held-out distances and
                     # the target's, beside the same for the random partition
                     # among the candidates and for plain leave-one-out.
                     wasserstein        = chosen$W,
                     wasserstein_random = res$w_random,
                     wasserstein_loo    = w_loo,
                     max_ecdf_excess    = max_excess,
                     max_ecdf_deficit   = max_deficit,
                     loo_ecdf_deficit    = reach$deficit,
                     ecdf_deficit_tol    = reach$tolerance,
                     target_separated    = reach$separated,
                     target_out_of_reach = reach$out_of_reach,
                     # TRUE when the warning above was raised: the folds are
                     # still closer to their training data than the
                     # prediction points are to the samples.
                     short_of_target     = short,
                     fold_sizes    = fold_sizes,
                     balance_ratio = max(fold_sizes) / min(fold_sizes),
                     balance_tol   = balance_tol,
                     target_median   = stats::median(g_target),
                     realised_median = stats::median(realised[fin]),
                     target_mean     = mean_target,
                     realised_mean   = mean_real,
                     target_distances   = g_target,
                     realised_distances = realised)))
  }

  stop("make_folds(): unsupported method.")
}


# -----------------------------------------------------------------------------
# make_folds() helpers: conditions a caller can catch, and what a block size
# yields on a given layout
# -----------------------------------------------------------------------------

#' Warn that fewer folds are returned than were asked for
#'
#' \code{.warn_and_log()} with a condition class and two fields, so that a
#' caller running many fold constructions can catch exactly this case with
#' \code{tryCatch(spatialkit_k_reduced = )} or
#' \code{withCallingHandlers()} and read \code{k_requested} and \code{k}
#' from the condition, without matching message text.
#'
#' @param k_requested,k The number of folds asked for and returned.
#' @param fmt,... The message, as for \code{sprintf()}.
#' @keywords internal
#' @noRd
.warn_k_reduced <- function(k_requested, k, fmt, ...) {
  .warn_classed("spatialkit_k_reduced", fmt, ...,
                data = list(k_requested = as.integer(k_requested),
                            k = as.integer(k)))
}


#' Warn that NNDM's prediction points are nearer than leave-one-out reaches
#'
#' Classed as \code{spatialkit_nndm_out_of_reach}, with the statistic and its
#' tolerance as fields \code{deficit} and \code{tolerance}, and
#' \code{separated}: whether every prediction distance is shorter than every
#' leave-one-out distance, which is what the warning rests on where the
#' tolerance is 1 or more.
#'
#' @param deficit,tolerance The largest amount by which the target ECDF runs
#'   above leave-one-out's, and the amount allowed.
#' @param separated Logical: complete separation of the two sets of
#'   distances.
#' @param fmt,... The message, as for \code{sprintf()}.
#' @keywords internal
#' @noRd
.warn_nndm_out_of_reach <- function(deficit, tolerance, separated, fmt, ...) {
  .warn_classed("spatialkit_nndm_out_of_reach", fmt, ...,
                data = list(deficit = deficit, tolerance = tolerance,
                            separated = isTRUE(separated)))
}


#' The prediction points of the distance-matched methods, and their distances
#'
#' What \code{method = "nndm"} and \code{method = "knndm"} share: the
#' prediction points reduced to POINTs and brought into the samples' CRS, and
#' the distance from each to its nearest sample, which is the target both
#' methods match.
#'
#' Prediction locations are frequently the grid CELLS a surface will be drawn
#' on, not their centres.  A POLYGON layer gives point-to-POLYGON distances,
#' zero for every cell that contains a training point: the target
#' distribution is pulled towards 0 and the matching degenerates towards
#' plain leave-one-out, the optimistic direction.  So the layer is pointized
#' first, as every other geometry input of \code{make_folds()} is.
#'
#' @param pts The projected samples.
#' @param prediction_points The argument as passed.
#' @param label \code{"nndm"} or \code{"knndm"}, for the messages.
#' @return A list: \code{pts} (stamped with the prediction points' CRS when
#'   they had none), \code{pred} and \code{g_target}, the finite distances.
#' @keywords internal
#' @noRd
.nndm_target <- function(pts, prediction_points, label) {
  if (is.null(prediction_points))
    stop("make_folds(", label, "): `prediction_points` is required -- ",
         if (identical(label, "knndm")) "kNNDM" else "NNDM", " matches ",
         "the training-to-test distance distribution to the distance ",
         "distribution from your actual prediction locations to the training ",
         "data. Use the grid you intend to predict onto (see ",
         "predict_surface()).", call. = FALSE)
  if (inherits(prediction_points, "sfc"))
    prediction_points <- sf::st_sf(geometry = prediction_points)
  if (!all(sf::st_geometry_type(prediction_points,
                                by_geometry = TRUE) == "POINT")) {
    .log_info(paste0("make_folds(%s): `prediction_points` are not POINTs; ",
                     "using representative points so the target distances ",
                     "are point-to-point."), label)
    prediction_points <- coerce_to_points(prediction_points, "auto")
  }
  if (is.na(sf::st_crs(pts)) && !is.na(sf::st_crs(prediction_points))) {
    pts <- .transform_or_stamp(pts, sf::st_crs(prediction_points),
                               what = "points_sf", caller = "make_folds")
  }
  # Aligned as `boundary` is: a CRS-less layer gets the warning, naming it.
  pred <- if (is.null(.crs_or_null(pts))) ensure_projected(prediction_points)
          else .transform_or_stamp(prediction_points, sf::st_crs(pts),
                                   what = "prediction_points",
                                   caller = "make_folds")
  pred <- sf::st_zm(pred, drop = TRUE, what = "ZM")

  # Target: distance from each prediction location to its nearest training
  # point.  This is the distance regime the model will actually face, and the
  # one an arbitrary fixed buffer has no reason to reproduce.
  g_target <- .nn_dist_to(pred, pts)
  g_target <- g_target[is.finite(g_target)]
  if (!length(g_target))
    stop("make_folds(", label, "): could not compute prediction-to-training ",
         "distances.", call. = FALSE)
  list(pts = pts, pred = pred, g_target = g_target)
}


#' Largest amount by which a target ECDF runs above that of a set of distances
#'
#' \code{max(G(r) - F(r))} over every distance \code{r} in either set, with
#' \code{G} the ECDF of the target and \code{F} that of \code{d}: the share of
#' prediction points within some distance of the data that \code{d} does not
#' reproduce.  Swapping the two arguments gives the excess.
#'
#' @param d Distances; non-finite values are dropped.
#' @param g_sorted The target distances, sorted.
#' @return A number in \code{[0, 1]}, or \code{NA} when \code{d} has no finite
#'   value.
#' @keywords internal
#' @noRd
.ecdf_deficit <- function(d, g_sorted) {
  d <- sort(d[is.finite(d)])
  if (!length(d)) return(NA_real_)
  at <- sort(unique(c(d, g_sorted)))
  max(findInterval(at, g_sorted) / length(g_sorted) -
        findInterval(at, d) / length(d))
}


#' Are the prediction points nearer to the data than any fold can hold out?
#'
#' A held-out point is never nearer to its training data than to its nearest
#' other sample, the distance plain leave-one-out gives.  So where the target
#' distribution runs above leave-one-out's, no fold design reaches it: not
#' NNDM, which only excludes neighbours, and not a k-fold design, whose
#' training sets are smaller still.  This is the rule both methods apply.
#'
#' Two ECDFs of one distribution differ by chance by about
#' \code{sqrt(1/n + 1/m)}.  Three times that was exceeded in none of 556
#' layouts whose folds match their target (uniform, clustered and graded
#' samples of 100 to 1000 points against grids and random prediction points,
#' 30 to 4000 of them; the largest ratio was 2.4) and in all 136 layouts
#' built out of reach (the smallest was 3.6).  One fold's worth, \code{1/n},
#' fires on every matching layout.
#'
#' That tolerance is 1 or more whenever \code{n * m / (n + m) <= 9}: for nine
#' samples or fewer whatever the target, and for up to eighteen with as many
#' prediction points.  A deficit is a difference of two shares and cannot
#' exceed 1, so there the rule above cannot fire and plain leave-one-out came
#' back without a word however near the prediction points were (12 samples,
#' prediction points 5 m from each: target median 5, realised 1097, no
#' warning).  The one case that needs no tolerance is complete separation:
#' every prediction distance shorter than every leave-one-out distance, a
#' deficit of exactly 1.  It is reported when the two sets are large enough
#' for that ordering to mean something, taken as
#' \code{choose(n + m, n) >= 1000} ways of interleaving them (17 prediction
#' points for 3 samples, 7 for 6, 4 for 12).  The count is a floor on the
#' sizes, not an error rate, and was set on layouts run through
#' \code{make_folds(method = "nndm")}: uniform, clustered and graded samples
#' of 3 to 18 points with prediction points over the samples' own extent, 300
#' layouts per cell.  With 2 to 5 random prediction points the two sets were
#' separated in 1055 of 8100 layouts, with medians as close as 1.06 to 1; the
#' floor passes over all but 12 of them.  With grids of 9 to 400 points it
#' warned in 262 of 900 layouts of 3 samples (the grids above the floor), 46
#' of 900 of 4, 18 of 1200 of 5, 6 of 1200 of 6 and none of 6000 of 8 to 18.
#' In every layout warned about, leave-one-out sat at least 1.38 times as far
#' as the target by the median and 1.8 times by the mean, so the statement
#' was true of those folds: a few samples hug the box they span.  Of layouts
#' built out of reach (prediction points 5 m from each of 9, 12 or 16
#' samples), it warned in 300 of 300, none of them by the share rule.
#'
#' @param loo_d Each sample's distance to its nearest other sample.
#' @param g_target The finite prediction-to-data distances.
#' @return A list: \code{deficit} (the largest amount by which the target
#'   ECDF runs above leave-one-out's), \code{tolerance}, \code{past_tol},
#'   \code{loo_fin} (the finite leave-one-out distances), \code{separated}
#'   and \code{out_of_reach}.
#' @keywords internal
#' @noRd
.nndm_reach <- function(loo_d, g_target) {
  n <- length(loo_d); n_g <- length(g_target)
  deficit   <- .ecdf_deficit(loo_d, sort(g_target))
  tolerance <- 3 * sqrt(1 / n + 1 / n_g)
  past_tol  <- is.finite(deficit) && deficit > tolerance
  loo_fin   <- loo_d[is.finite(loo_d)]
  separated <- length(loo_fin) > 0L && max(g_target) < min(loo_fin)
  sep_counts <- separated && lchoose(n + n_g, n) >= log(1000)
  list(deficit = deficit, tolerance = tolerance, past_tol = past_tol,
       loo_fin = loo_fin, separated = separated,
       out_of_reach = past_tol || sep_counts)
}


#' The held-out rows of a leave-one-out design
#'
#' Every row, or with \code{max_test} a random subset of that size, sorted.
#' The draw uses the session RNG, which \code{make_folds()} has seeded when
#' it was given a \code{seed}.
#'
#' @param n Number of rows.
#' @param max_test \code{NULL}, or the number of held-out points to keep.
#' @return Sorted row positions.
#' @keywords internal
#' @noRd
.loo_scored_rows <- function(n, max_test = NULL) {
  if (is.null(max_test) || max_test >= n) return(seq_len(n))
  sort(sample.int(n, max_test))
}


#' Pack blocks into k folds, largest first
#'
#' The longest-processing-time rule for multiway partitioning: the blocks
#' are taken in decreasing order of size and each goes to the fold with the
#' fewest points so far, ties broken at random.  \code{"block_kfold"} packs
#' its grid cells this way and \code{"knndm"} its clusters.
#'
#' @param sizes Integer points per block.
#' @param k Number of folds.
#' @return A list: \code{fold_of} (the fold of each block), \code{loads}
#'   (points per fold) and \code{order} (the order the blocks were placed
#'   in).
#' @keywords internal
#' @noRd
.pack_largest_first <- function(sizes, k) {
  ord <- order(sizes, decreasing = TRUE)
  loads <- integer(k); fold_of <- integer(length(sizes))
  for (b in ord) {
    j <- which(loads == min(loads))
    if (length(j) > 1) j <- sample(j, 1L)
    fold_of[b] <- j
    loads[j] <- loads[j] + sizes[b]
  }
  list(fold_of = fold_of, loads = loads, order = ord)
}


#' The fold of each row under a random k-fold partition
#'
#' The partition \code{method = "random_kfold"} builds: one permutation of
#' the rows cut into \code{k} runs whose sizes differ by at most one.  Drawn
#' from the same RNG state it gives the same folds, which is how the random
#' candidate of \code{"knndm"} is the partition \code{"random_kfold"} returns
#' for the same \code{seed}.
#'
#' @param n Number of rows.
#' @param k Number of folds, at most \code{n}.
#' @return Integer fold labels, one per row.
#' @keywords internal
#' @noRd
.random_fold_labels <- function(n, k) {
  idx <- sample.int(n, n)
  sizes <- rep(floor(n / k), k)
  remainder <- n - sum(sizes)
  if (remainder > 0) sizes[seq_len(remainder)] <- sizes[seq_len(remainder)] + 1L
  lab <- integer(n)
  lab[idx] <- rep.int(seq_len(k), sizes)
  lab
}


#' Nearest row of \code{data} to each row of \code{query}, in the plane
#'
#' \pkg{FNN}'s kd-tree when it is installed, and otherwise
#' \code{sf::st_nearest_feature()}, which is indexed too.  The distance is
#' Euclidean in both.
#'
#' @param query,data Two-column coordinate matrices.
#' @param use_fnn Whether to use \pkg{FNN}.
#' @return A list: \code{index} (row of \code{data}) and \code{dist}.
#' @keywords internal
#' @noRd
.nn_xy <- function(query, data,
                   use_fnn = requireNamespace("FNN", quietly = TRUE)) {
  if (isTRUE(use_fnn)) {
    kn <- FNN::get.knnx(data, query, k = 1L)
    return(list(index = as.integer(kn$nn.index[, 1L]),
                dist  = as.numeric(kn$nn.dist[, 1L])))
  }
  idx <- sf::st_nearest_feature(.xy_points(query), .xy_points(data))
  list(index = idx,
       dist  = sqrt((query[, 1L] - data[idx, 1L])^2 +
                    (query[, 2L] - data[idx, 2L])^2))
}


#' A coordinate matrix as planar points with no CRS
#' @keywords internal
#' @noRd
.xy_points <- function(m)
  sf::st_geometry(sf::st_as_sf(data.frame(x = m[, 1L], y = m[, 2L]),
                               coords = c("x", "y")))


#' Each point's distance to the nearest other point
#'
#' The held-out distance of plain leave-one-out.  A co-located duplicate is
#' at distance 0.
#'
#' @param xy Two-column coordinate matrix with at least two rows.
#' @param use_fnn Whether to use \pkg{FNN}.
#' @return A numeric vector, one distance per row.
#' @keywords internal
#' @noRd
.loo_dist_xy <- function(xy, use_fnn = requireNamespace("FNN", quietly = TRUE)) {
  if (isTRUE(use_fnn))
    return(as.numeric(FNN::get.knn(xy, k = 1L)$nn.dist[, 1L]))
  # sf::st_nearest_feature(x) with `y` missing gives each feature's nearest
  # among the others, in n^2 time: 68 s at 8,000 points.  The indexed search
  # needs two sets.  Two different rows differ in at least one bit of their
  # (zero-based) row number, so a row's nearest other point lies in the
  # opposite half for one of those bits, and no point of an opposite half is
  # nearer than it.  The smallest opposite-half distance over the bits is
  # therefore the distance wanted, from log2(n) pairs of indexed searches
  # (0.5 s at 8,000 points).
  n <- nrow(xy)
  row0 <- seq_len(n) - 1L
  d <- rep(Inf, n)
  for (bit in seq_len(max(1L, ceiling(log2(n)))) - 1L) {
    half <- bitwAnd(bitwShiftR(row0, bit), 1L) + 1L
    d <- pmin(d, .fold_nn_dist(xy, half, 2L, FALSE), na.rm = TRUE)
  }
  d
}


#' Each point's distance to the nearest point of another fold
#'
#' @param xy Two-column coordinate matrix.
#' @param lab Integer fold label per row, in \code{1..k}.
#' @param k Number of folds.
#' @param use_fnn Whether to use \pkg{FNN}.
#' @return A numeric vector, one distance per row; \code{NA} for a point
#'   whose fold holds every row.
#' @keywords internal
#' @noRd
.knndm_heldout_dist <- function(xy, lab, k, use_fnn) {
  # One point per fold is leave-one-out, and one nearest-neighbour search
  # then stands for the n searches of the loop over folds.
  if (k == nrow(xy) && !anyDuplicated(lab)) return(.loo_dist_xy(xy, use_fnn))
  .fold_nn_dist(xy, lab, k, use_fnn)
}


#' The loop behind \code{.knndm_heldout_dist()}: one search per fold
#' @keywords internal
#' @noRd
.fold_nn_dist <- function(xy, lab, k, use_fnn) {
  d <- rep(NA_real_, nrow(xy))
  for (f in seq_len(k)) {
    te <- lab == f
    if (!any(te) || all(te)) next
    d[te] <- .nn_xy(xy[te, , drop = FALSE], xy[!te, , drop = FALSE],
                    use_fnn)$dist
  }
  d
}


#' Cluster counts tried by kNNDM
#'
#' From \code{k} clusters (one per fold) to one fewer than the number of
#' distinct locations (nearly every point a cluster of its own), evenly
#' spaced on the log scale.
#'
#' @param k Number of folds.
#' @param n_distinct Number of distinct locations.
#' @param len Number of counts before duplicates are dropped.
#' @return Increasing integers, or \code{integer(0)} when there are too few
#'   distinct locations to form more than \code{k} clusters.
#' @keywords internal
#' @noRd
.knndm_cluster_counts <- function(k, n_distinct, len = 50L) {
  top <- n_distinct - 1L
  if (top < k) return(integer(0))
  unique(as.integer(round(exp(seq(log(k), log(top), length.out = len)))))
}


#' k-means clusters of the coordinates, for kNNDM on large samples
#'
#' Used to reduce a large sample to groups Ward's method can start from, and
#' for cluster counts above the number of those groups.  Lloyd's algorithm from \code{q} distinct sample locations drawn at random,
#' stopped after \code{iter_max} passes or when no point changes cluster.
#' The nearest-centre step is an indexed search, so a pass costs about
#' \code{n log q}, where \code{stats::kmeans()} compares every point with
#' every centre (\code{n q}; at 20,000 points and 2,000 centres, 0.8 s
#' against 0.09 s for ten passes).  The clusters are candidates to be scored
#' by their distance to the target, so a run that has not converged costs
#' nothing but a slightly different candidate.  A centre that attracts no
#' point is dropped, so fewer than \code{q} clusters can come back.
#'
#' @param xy Two-column coordinate matrix.
#' @param xy_u Its distinct rows.
#' @param q Number of clusters asked for, at most \code{nrow(xy_u)}.
#' @param use_fnn Whether to use \pkg{FNN}.
#' @param iter_max Passes.
#' @return Integer cluster labels \code{1..q'}, one per row of \code{xy}.
#' @keywords internal
#' @noRd
.knndm_kmeans <- function(xy, xy_u, q, use_fnn, iter_max = 10L) {
  cen <- xy_u[sample.int(nrow(xy_u), q), , drop = FALSE]
  cl <- integer(0)
  for (it in seq_len(iter_max)) {
    new <- .nn_xy(xy, cen, use_fnn)$index
    new <- match(new, sort(unique(new)))
    if (identical(new, cl)) break
    cl <- new
    cen <- unname(rowsum(xy, cl) / tabulate(cl))
  }
  cl
}


#' Ward's clustering continued from groups of points
#'
#' \code{stats::hclust(method = "ward.D2")} started in the middle of the
#' dendrogram: the groups stand for merges already made, and the tree above
#' them is the one Ward's method builds from there.  The dissimilarity of two
#' groups of sizes \eqn{a} and \eqn{b} whose centroids are \eqn{d} apart is
#' \eqn{\sqrt{2ab/(a+b)}\,d}, which is \eqn{d} for two single points and is
#' what the Lance-Williams update of \code{"ward.D2"} carries, and the sizes
#' go in as \code{members}.  Started from an 80-cluster cut of a full Ward
#' tree on 600 points it reproduces every coarser cut of that tree.
#'
#' @param xy Two-column coordinate matrix.
#' @param grp Integer group label per row, \code{1..m} with \code{m >= 2}.
#' @return An \code{hclust} tree over the \code{m} groups.
#' @keywords internal
#' @noRd
.ward_from_groups <- function(xy, grp) {
  sz  <- tabulate(grp)
  m   <- length(sz)
  dc  <- stats::dist(rowsum(xy, grp) / sz)
  # A `dist` object runs down the columns of the lower triangle.
  col <- rep(seq_len(m - 1L), times = (m - 1L):1L)
  a   <- as.numeric(sz[col])
  b   <- as.numeric(sz[sequence((m - 1L):1L) + col])
  col <- NULL
  dc[] <- sqrt(2 * a * b / (a + b)) * as.numeric(dc)
  stats::hclust(dc, method = "ward.D2", members = sz)
}


#' Choose k folds whose held-out distances are closest to a target
#'
#' The core of \code{make_folds(method = "knndm")}, after Linnenbrink et al.
#' (2024).  The candidates are a random k-fold partition and, for each of a
#' ladder of cluster counts, the samples clustered on their coordinates and
#' the clusters packed into \code{k} folds largest first.  Each admissible
#' candidate is scored by the Wasserstein distance \eqn{W} between its
#' held-out distances (each point's distance to the nearest point of another
#' fold) and the target distances, and the smallest \eqn{W} wins; a tie goes
#' to the earlier candidate, the random partition first.
#'
#' A cluster candidate is admissible when its largest fold leaves at least
#' \code{min_train} of the rows for training (or holds no more than an equal
#' share, \code{ceiling(n / k)}, where \code{min_train} allows less than
#' that) and its largest fold holds at most \code{balance_tol} times as many
#' points as its smallest.  The random partition is always admissible.
#'
#' Clustering is Ward's (\code{stats::hclust(method = "ward.D2")}), cut at
#' each count, up to \code{ward_max} points, which needs the n by n distance
#' matrix (64 MB at 4,000 points, 1.6 GB at 20,000).  Above that the samples
#' are first reduced to \code{ward_max} groups by k-means
#' (\code{.knndm_kmeans()}) and Ward's method is continued from those groups
#' (\code{.ward_from_groups()}); a count above the number of groups is a
#' k-means run of its own.  On a 40-cluster layout of unequal cluster sizes
#' with a grid target (10 folds, 3 draws each) the folds chosen had \eqn{W}
#' at 0.09 to 0.12 of the mean target distance with Ward on all 8,000 points
#' and 0.09 to 0.11 from 4,000 groups, where k-means at every count gave 0.16
#' to 0.19; at 30,000 points, 0.11 to 0.15 from 4,000 groups against 0.18 to
#' 0.25; at 100,000 points 0.27 (one draw), the count chosen (12,648) being
#' above the groups.  A random partition stood at 0.58 to 0.62 throughout.
#' With and without \pkg{FNN} the folds were identical in all six draws at
#' 8,000 and 30,000 points.
#'
#' @param xy Two-column coordinate matrix of the samples.
#' @param g_target Finite prediction-to-data distances.
#' @param k Number of folds, between 2 and \code{nrow(xy)}.
#' @param min_train,balance_tol As in \code{make_folds()}.
#' @param use_fnn Whether to use \pkg{FNN} for the nearest-neighbour searches.
#' @param ward_max Largest sample clustered hierarchically.
#' @param n_ladder Number of cluster counts tried, before duplicates.
#' @return A list: \code{fold} (label per row), \code{realised} (held-out
#'   distance per row), \code{chosen} (row of \code{candidates}),
#'   \code{candidates} (a data frame: \code{candidate}, \code{n_clusters},
#'   \code{largest_fold}, \code{balance_ratio}, \code{admissible},
#'   \code{W}), \code{w_random}, \code{clustering} (\code{"ward"},
#'   \code{"kmeans+ward"}, or \code{NA} when no cluster count was tried) and
#'   \code{loo_distances}.
#' @keywords internal
#' @noRd
.knndm_folds <- function(xy, g_target, k, min_train = 0.5, balance_tol = 3,
                         use_fnn = requireNamespace("FNN", quietly = TRUE),
                         ward_max = 4000L, n_ladder = 50L) {
  n <- nrow(xy)
  w_of <- function(d) {
    d <- d[is.finite(d)]
    if (length(d)) .fold_sep_compare(d, g_target)$w1 else NA_real_
  }
  one <- function(candidate, n_clusters, loads, admissible, W)
    data.frame(candidate = candidate, n_clusters = as.integer(n_clusters),
               largest_fold = if (length(loads)) max(loads) else NA_integer_,
               balance_ratio = if (length(loads) && min(loads) > 0)
                 max(loads) / min(loads) else NA_real_,
               admissible = admissible, W = W, stringsAsFactors = FALSE)
  # The largest fold a cluster candidate may have.
  cap <- max(ceiling(n / k), floor((1 - min_train) * n + 1e-9))

  # The random partition is drawn before anything else touches the RNG, so
  # that it is the one "random_kfold" builds from the same seed.
  lab_best <- .random_fold_labels(n, k)
  d_best   <- .knndm_heldout_dist(xy, lab_best, k, use_fnn)
  w_random <- w_of(d_best)
  cand <- list(one("random", NA_integer_, tabulate(lab_best, nbins = k), TRUE,
                   w_random))
  w_best <- if (is.finite(w_random)) w_random else Inf
  best <- 1L

  distinct <- !duplicated(xy)
  qs <- .knndm_cluster_counts(k, sum(distinct), n_ladder)
  clustering <- NA_character_
  if (length(qs)) {
    use_ward <- n <= ward_max
    clustering <- if (use_ward) "ward" else "kmeans+ward"
    if (use_ward) {
      tree <- stats::hclust(stats::dist(xy), method = "ward.D2")
      groups <- NULL; n_groups <- n
    } else {
      # Too many points for the n by n distance matrix: k-means reduces them
      # to at most `ward_max` groups and Ward's method carries on from those.
      xy_u   <- xy[distinct, , drop = FALSE]
      groups <- .knndm_kmeans(xy, xy_u, min(ward_max, nrow(xy_u) - 1L), use_fnn)
      n_groups <- max(groups)
      tree <- if (n_groups >= 2L) .ward_from_groups(xy, groups) else NULL
    }
    for (q in qs) {
      cl <- if (use_ward) stats::cutree(tree, k = q)
            else if (!is.null(tree) && q <= n_groups)
              stats::cutree(tree, k = q)[groups]
            # More clusters than groups: k-means at that count.
            else .knndm_kmeans(xy, xy_u, q, use_fnn)
      sizes <- tabulate(cl)
      if (length(sizes) < k) {
        cand[[length(cand) + 1L]] <- one("clusters", length(sizes), integer(0),
                                         FALSE, NA_real_)
        next
      }
      pk  <- .pack_largest_first(sizes, k)
      ok  <- min(pk$loads) > 0L && max(pk$loads) <= cap &&
             max(pk$loads) / min(pk$loads) <= balance_tol
      w_q <- NA_real_
      if (ok) {
        lab <- pk$fold_of[cl]
        d   <- .knndm_heldout_dist(xy, lab, k, use_fnn)
        w_q <- w_of(d)
        if (is.finite(w_q) && w_q < w_best) {
          w_best <- w_q; lab_best <- lab; d_best <- d
          best <- length(cand) + 1L
        }
      }
      cand[[length(cand) + 1L]] <- one("clusters", length(sizes), pk$loads, ok,
                                       w_q)
    }
  }
  list(fold = lab_best, realised = d_best, chosen = best,
       candidates = do.call(rbind, cand), w_random = w_random,
       clustering = clustering,
       loo_distances = .loo_dist_xy(xy, use_fnn))
}


#' Warn that kNNDM's best candidate is still closer than the target
#'
#' Classed as \code{spatialkit_knndm_short_of_target}, with the two mean
#' distances as fields \code{held_out_mean} and \code{target_mean}.
#'
#' @param held_out_mean,target_mean Mean held-out-to-training distance of the
#'   folds returned, and mean prediction-to-data distance.
#' @param fmt,... The message, as for \code{sprintf()}.
#' @keywords internal
#' @noRd
.warn_knndm_short <- function(held_out_mean, target_mean, fmt, ...) {
  .warn_classed("spatialkit_knndm_short_of_target", fmt, ...,
                data = list(held_out_mean = held_out_mean,
                            target_mean = target_mean))
}


#' The extent of a bbox, for a message
#' @keywords internal
#' @noRd
.block_extent_chr <- function(bb) {
  sprintf("%s x %s",
          format(signif(as.numeric(bb["xmax"] - bb["xmin"]), 4)),
          format(signif(as.numeric(bb["ymax"] - bb["ymin"]), 4)))
}


#' Lay the block grid of block_kfold over a region
#'
#' The cells of an \code{nx} by \code{ny} grid over the bounding box of
#' \code{reg}, clipped to it.  Clipping to a boundary drops the cells outside
#' it and renumbers the rest, so a cell's index in the full grid is carried
#' from the "idx" attribute: that is what \code{params$blocks$source_row}
#' promises.  Where the boundary only touches a cell at a corner or along an
#' edge the intersection is a zero-area POINT or LINESTRING, which used to be
#' kept as a block of its own, packed into a fold when empty blocks were kept
#' and able to catch a point lying exactly on the boundary as a one-point
#' block.  Only the areal pieces are kept, unless there are none (points on
#' one straight line have a region of zero area, and every piece of it is a
#' line).
#'
#' @param reg The region: an sf/sfc polygon layer (or the points' bbox).
#' @param nx,ny Grid dimensions.
#' @return A list with \code{grid_sf} (the blocks) and \code{cell_idx} (each
#'   block's index in the full grid).
#' @keywords internal
#' @noRd
.block_grid <- function(reg, nx, ny) {
  grid <- sf::st_make_grid(reg, n = c(nx, ny), what = "polygons", square = TRUE)
  grid <- .safe_make_valid(grid)
  reg_union <- .safe_make_valid(sf::st_union(reg))
  grid <- suppressWarnings(sf::st_intersection(grid, reg_union))
  cell_idx <- attr(grid, "idx")
  cell_idx <- if (is.matrix(cell_idx) && nrow(cell_idx) == length(grid))
    as.integer(cell_idx[, 1L]) else seq_along(grid)
  areal <- sf::st_dimension(grid) %in% 2L
  if (any(areal) && !all(areal)) {
    grid <- grid[areal]; cell_idx <- cell_idx[areal]
  }
  list(grid_sf = sf::st_as_sf(grid), cell_idx = cell_idx)
}


#' How many blocks of an nx by ny grid hold points
#'
#' Builds the grid as \code{make_folds()} does and assigns the points as it
#' does (the first cell that contains each), so the count is the number of
#' folds that grid can support.  \code{NA} for a grid too large to build.
#'
#' @param pts Projected points.
#' @param reg The block region.
#' @param nx,ny Grid dimensions.
#' @return An integer count, or \code{NA}.
#' @keywords internal
#' @noRd
.block_occupied <- function(pts, reg, nx, ny) {
  n_cells <- as.numeric(nx) * as.numeric(ny)
  if (!is.finite(n_cells) || n_cells > 1e5) return(NA_integer_)
  g <- tryCatch(.block_grid(reg, nx, ny), error = function(e) NULL)
  if (is.null(g)) return(NA_integer_)
  hits <- sf::st_intersects(pts, g$grid_sf)
  id <- vapply(hits, function(ix) if (length(ix)) ix[1] else NA_integer_, 1L)
  length(unique(id[!is.na(id)]))
}


#' The same count from the coordinates alone, for searching
#'
#' Bins the points into the cells of the grid over \code{bb} by arithmetic.
#' A point on a shared edge goes to the lower cell, as the first-cell rule of
#' \code{make_folds()} sends it.  Rounding can still differ from the polygon
#' test for a point within an ulp of an edge, so a size found this way is
#' confirmed with \code{.block_occupied()} before it is quoted.
#'
#' @param xy Coordinate matrix.
#' @param bb The region's bbox.
#' @param nx,ny Grid dimensions.
#' @return An integer count.
#' @keywords internal
#' @noRd
.block_occupied_fast <- function(xy, bb, nx, ny) {
  w <- as.numeric(bb["xmax"] - bb["xmin"]); h <- as.numeric(bb["ymax"] - bb["ymin"])
  ix <- if (w > 0) pmin(nx, pmax(1, ceiling((xy[, 1L] - as.numeric(bb["xmin"])) / (w / nx))))
        else rep(1, nrow(xy))
  iy <- if (h > 0) pmin(ny, pmax(1, ceiling((xy[, 2L] - as.numeric(bb["ymin"])) / (h / ny))))
        else rep(1, nrow(xy))
  length(unique((iy - 1) * nx + ix))
}


#' A length no smaller than \code{x}, rounded for a message
#'
#' The size a leakage warning names must be at least the range, and must read
#' back as the number that was tried.  Rounded up to a whole unit from 10
#' upward, and to three significant figures below it (a range of 0.3712 km
#' used to be quoted as a block size of 1).
#'
#' @param x A positive number.
#' @return A number no smaller than \code{x}.
#' @keywords internal
#' @noRd
.block_size_at_least <- function(x) {
  if (!is.finite(x) || x <= 0) return(x)
  if (x >= 10) return(ceiling(x))
  p <- 10^(floor(log10(x)) - 2)
  ceiling(x / p) * p
}


#' The largest block size, below a bound, that gives k blocks holding points
#'
#' The grid \code{.block_dims_from_size()} builds changes only at sizes
#' \code{w / i} and \code{h / j}, so those (a hair under each, down to a
#' two-hundredth of the extent) are the candidates, largest first.  A
#' candidate whose grid has at least \code{k} occupied cells is rounded down
#' to three significant figures, as it will be quoted, and the grid at that
#' quoted number is then built and counted for real: the size returned is one
#' that \code{make_folds(block_size = )} has in effect already been run with.
#'
#' @param pts Projected points.
#' @param reg The block region.
#' @param bb Its bbox.
#' @param k Blocks needed.
#' @param below Sizes at or above this are not considered.
#' @return \code{NULL} when none is found, else a list with \code{size}
#'   (numeric), \code{chr} (as quoted), \code{nx} and \code{ny}.
#' @keywords internal
#' @noRd
.block_largest_size_for_k <- function(pts, reg, bb, k, below = Inf) {
  below <- as.numeric(below)[1L]
  w <- as.numeric(bb["xmax"] - bb["xmin"]); h <- as.numeric(bb["ymax"] - bb["ymin"])
  cand <- c(if (is.finite(w) && w > 0) w / seq_len(200L),
            if (is.finite(h) && h > 0) h / seq_len(200L)) * (1 - 1e-9)
  cand <- sort(unique(cand[is.finite(cand) & cand > 0 & cand < below]),
               decreasing = TRUE)
  if (!length(cand)) return(NULL)
  xy <- sf::st_coordinates(pts)[, 1:2, drop = FALSE]
  confirmed <- 0L
  for (b in cand) {
    d <- .block_dims_from_size(bb, b)
    if (as.numeric(d$nx) * as.numeric(d$ny) < k) next
    if (.block_occupied_fast(xy, bb, d$nx, d$ny) < k) next
    chr  <- format(.signif_down(b), scientific = FALSE, trim = TRUE)
    size <- suppressWarnings(as.numeric(chr))
    if (!is.finite(size) || size <= 0 || size >= below) next
    ds  <- .block_dims_from_size(bb, size)
    occ <- .block_occupied(pts, reg, ds$nx, ds$ny)
    if (!is.na(occ) && occ >= k)
      return(list(size = size, chr = chr, nx = ds$nx, ny = ds$ny))
    confirmed <- confirmed + 1L
    if (confirmed >= 8L) break
  }
  NULL
}


#' Name the fold method that does not depend on a block size
#'
#' @param n Number of points.
#' @return Character(1).
#' @keywords internal
#' @noRd
.block_other_method <- function(n) {
  if (n <= 5000L)
    return(paste0("method = \"nndm\" does not depend on a block size: it ",
                  "builds leave-one-out folds matched to the distances from ",
                  "`prediction_points` to the data."))
  # Above the leave-one-out method's row limit the k-fold version is the one
  # that runs.  The sentence used to end at the limit, naming a method and
  # then saying it could not be used.
  sprintf(paste0("method = \"knndm\" does not depend on a block size: it ",
                 "builds k folds matched to the distances from ",
                 "`prediction_points` to the data. (method = \"nndm\", its ",
                 "leave-one-out version, takes at most 5000 rows, and this ",
                 "layer has %d.)"), n)
}


#' Warn that block-building arguments arrived beside supplied folds
#'
#' \code{block_size} and \code{auto_range} are read only when a
#' \code{cv_*()} function builds its own folds.  Beside supplied
#' \code{folds} they used to be dropped without a word, so
#' \code{cv_spatial(folds = random_folds, auto_range = TRUE)} ran on the
#' random folds and looked like a range-sized blocked run.  \code{k} is not
#' checked: the package's own callers always pass it with their folds, and
#' it is documented as used only when \code{folds} is \code{NULL}.  The
#' warning carries the class \code{spatialkit_fold_args_ignored}.
#'
#' @param folds The \code{folds} argument as the caller received it.
#' @param block_size,auto_range The two arguments.
#' @param caller Name for the message.
#' @return \code{TRUE}, invisibly, when it warned.
#' @keywords internal
#' @noRd
.warn_fold_args_ignored <- function(folds, block_size, auto_range, caller) {
  if (is.null(folds)) return(invisible(FALSE))
  given <- c(if (!is.null(block_size)) "`block_size`",
             if (isTRUE(auto_range)) "`auto_range = TRUE`")
  if (!length(given)) return(invisible(FALSE))
  several <- length(given) > 1L
  msg <- sprintf(paste0(
    "%s(): %s %s no effect because `folds` is supplied: the folds are used ",
    "as given and no blocks are built. Pass %s to make_folds() when building ",
    "the folds, or leave `folds` NULL to have them built here."),
    caller, paste(given, collapse = " and "), if (several) "have" else "has",
    if (several) "them" else "it")
  .warn_classed("spatialkit_fold_args_ignored", "%s", msg)
  invisible(TRUE)
}


#' What following a leakage warning's advice yields on this layout
#'
#' The leakage warnings of \code{make_folds(method = "block_kfold")} used to
#' end "Pass block_size = <range> or auto_range = TRUE" whatever the extent.
#' This builds the grid that size asks for, counts the blocks that hold
#' points, and words the advice by the outcome: the setting and what it
#' builds when that gives \code{k} folds; otherwise how many folds it gives
#' (or that it stops with the single-block error), the largest block size
#' that still gives \code{k} folds and how far below the range that is, and
#' the method that does not depend on a block size.
#'
#' @param pts Projected points.
#' @param reg The block region.
#' @param bb Its bbox.
#' @param sac_range The estimated range.
#' @param k Folds asked for.
#' @param mention_auto Whether \code{auto_range = TRUE} alone would size the
#'   blocks from the range (it would not beside an explicit \code{block_size}).
#' @return Character(1), one or more sentences.
#' @keywords internal
#' @noRd
.block_range_advice <- function(pts, reg, bb, sac_range, k, mention_auto = TRUE) {
  # An estimate_sac_range() result carries a class and a dozen attributes,
  # which arithmetic and format() would drag into the text.
  sac_range <- as.numeric(sac_range)[1L]
  ext <- .block_extent_chr(bb)
  size_chr <- format(.block_size_at_least(sac_range), scientific = FALSE, trim = TRUE)
  size <- as.numeric(size_chr)
  d   <- .block_dims_from_size(bb, size)
  occ <- .block_occupied(pts, reg, d$nx, d$ny)
  # auto_range = TRUE sizes the grid from the range itself, not from the
  # rounded figure quoted here; it is named only when both build one grid.
  da <- .block_dims_from_size(bb, sac_range)
  auto <- isTRUE(mention_auto) && identical(c(da$nx, da$ny), c(d$nx, d$ny))
  setting <- if (auto) sprintf("block_size = %s or auto_range = TRUE", size_chr)
             else sprintf("block_size = %s", size_chr)
  grid_chr <- sprintf("%s x %s", format(d$nx, scientific = FALSE),
                      format(d$ny, scientific = FALSE))
  if (is.na(occ)) return(sprintf("Pass %s.", setting))
  if (occ >= k)
    return(sprintf(paste0("Pass %s: over this %s extent that builds a %s grid ",
                          "in which %d blocks hold points, enough for k = %d ",
                          "folds."), setting, ext, grid_chr, occ, k))
  alt <- .block_largest_size_for_k(pts, reg, bb, k, below = size)
  alt_txt <- if (is.null(alt))
    sprintf("No block size down to 1/200 of the extent gives %d blocks that hold points.", k)
  else sprintf(paste0("The largest block_size that gives %d folds is %s, ",
                      "%.0f%% of the range, so its folds are separated by ",
                      "less than the range."), k, alt$chr, 100 * alt$size / as.numeric(sac_range))
  other <- .block_other_method(nrow(pts))
  if (occ >= 2L)
    return(sprintf(paste0("Blocks of the range leave %d blocks that hold ",
                          "points over this %s extent (a %s grid), so %s ",
                          "returns %d folds, not %d. %s %s"),
                   occ, ext, grid_chr, setting, occ, k, alt_txt, other))
  sprintf(paste0("A block of the range leaves a single block that holds ",
                 "points over this %s extent, so block folds cannot be ",
                 "separated by the range here: %s stops with an error. %s %s"),
          ext, setting, alt_txt, other)
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
#' cross-validation, with no tuning on the full data first.
#'
#' Folds default to spatial blocks (\code{\link{make_folds}(method =
#' "block_kfold")}), not random ones.  With autocorrelated data a random split
#' leaves a held-out point's neighbours in the training set and the score comes
#' back flattering.  Use \code{\link{cv_bayes}()} for the same treatment
#' of a Bayesian GP model, \code{\link{cv_rf}()} for a forest, and
#' \code{\link{compare_models_cv}()} to score several backends on one set of
#' folds.
#'
#' @param data_sf An sf object.
#' @param response_var Response column name.
#' @param predictor_vars Predictor column names.
#' @param folds Optional fold definitions, in any of three shapes: a
#'   \code{\link{make_folds}()} return value; a bare list of
#'   \code{list(train =, test =)} pairs of \code{..row_id} values; or a vector
#'   of fold labels, one per row, which becomes leave-that-label-out splits.
#'   The label vector is how folds built by another package are used here
#'   (\code{blockCV::cv_spatial()} returns one as \code{$folds_ids}), since its
#'   \code{$folds_list} holds two \emph{unnamed} vectors per fold and is
#'   refused by name.  Train and test must be disjoint: a fold that trains on
#'   its own test rows is not a cross-validation split and is refused with an
#'   error.  IDs naming no row in the prepared data are dropped with
#'   a logged count (expected when rows were removed for missing values; a sign
#'   the folds came from other data when they were not).
#' @param k Number of folds to build when \code{folds} is \code{NULL}; not
#'   read when \code{folds} is supplied.  Default 5.
#' @param seed RNG seed. Default 123.
#' @param adaptive Logical; use adaptive bandwidth. Default TRUE.
#' @param bandwidth Optional bandwidth, applied to every fold. For
#'   \code{adaptive = TRUE} an integer number of neighbours; for
#'   \code{adaptive = FALSE} a distance in the units of the \strong{projected}
#'   CRS the folds are fitted in (geographic input is projected first, so a
#'   value in degrees is read as metres). \code{NULL} (default) selects a
#'   bandwidth per fold from that fold's training rows, by the search
#'   \code{\link{fit_gwr_model}()} runs (its section Bandwidth selection:
#'   AICc, or the leave-one-location-out sum of squares when observations
#'   share a location), and reports it in \code{fold_metrics$bandwidth}.
#' @param kernel Kernel function type.
#' @param boundary Optional polygonal sf/sfc for CRS alignment.
#' @param pointize Geometry coercion strategy.
#' @param block_size Optional minimum block edge length (projected CRS
#'   units) for the blocks built when \code{folds} is \code{NULL}; see
#'   \code{\link{make_folds}()}.  It is a size you choose and is not
#'   compared with the autocorrelation range unless one can be estimated, in
#'   which case a smaller size is warned about.  Beside supplied
#'   \code{folds} it has nothing to size and raises a warning (class
#'   \code{spatialkit_fold_args_ignored}).
#' @param auto_range Logical.  If \code{TRUE} and \code{folds} is \code{NULL},
#'   estimate the autocorrelation range and use it as the minimum block size.
#'   Default \code{FALSE}.  \code{TRUE} beside supplied \code{folds} raises
#'   the same warning as \code{block_size}.
#' @param parallel Accepted so that every \code{cv_*()} function takes the
#'   same arguments, but the GWR folds always run one after another in this
#'   R process, where the warnings \code{predict()} raises for each fold are
#'   collected into one.  A value asking for more than one core
#'   raises a warning saying so.  Default \code{FALSE}.
#' @param metrics Optional scoring function of your own, a
#'   \code{function(y, yhat)} returning a named numeric vector; its names
#'   become columns of \code{fold_metrics} (per fold) and \code{overall}
#'   (pooled) beside the built-in ones.  See \strong{Your own metrics} on
#'   \code{\link{cv_spatial}()} for the contract.
#' @inheritSection model_metrics Percentage errors on responses with zeros
#' @return A list with \code{overall}, \code{fold_metrics},
#'   \code{predictions}, \code{folds}, \code{n_folds_attempted},
#'   \code{n_folds_succeeded}, \code{fold_status}, \code{orphan_rows},
#'   \code{n_unknown_ids}, \code{n_dropped}, \code{formula} and
#'   \code{adaptive}.  The two
#'   fold counts make a run where every fold failed visible in the return
#'   value itself, beyond the warning, since \code{overall} is a
#'   well-formed all-\code{NA} row either way, and \code{fold_status} (one
#'   row per fold: \code{fold}, \code{status}, \code{message}) says why each
#'   missing fold is missing.  See \code{\link{cv_spatial}} for the five
#'   statuses and for \code{orphan_rows}.  \code{Adj_R2} is \code{NA}
#'   in both \code{overall} and \code{fold_metrics}: the pooled predictions
#'   have no single parameter count, and a GWR's effective parameter count is
#'   not its predictor count (see \code{\link{cv_spatial}}).  A held-out
#'   location whose local regression cannot be estimated has \code{yhat}
#'   \code{NA}, and one beyond the reach of its fold's fit is predicted and
#'   scored; each case raises one warning for the whole run (classes
#'   \code{"spatialkit_gwr_not_estimable"} and
#'   \code{"spatialkit_gwr_extrapolation"}), with the totals over the folds
#'   that had any.  See \code{\link{predict.gwr_fit}}.
#' @family cross-validation
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
#'   cv <- cv_gwr(dat, "price", "elev", k = 3, bandwidth = 30)
#'   print(cv$overall)       # pooled over the held-out rows of every fold
#'   cv$fold_metrics         # and fold by fold
#' }
#' @export
cv_gwr <- function(data_sf, response_var, predictor_vars,
                   folds = NULL, k = 5, seed = 123,
                   adaptive = TRUE, bandwidth = NULL,
                   kernel = c("bisquare", "gaussian", "tricube",
                              "boxcar", "exponential"),
                   boundary = NULL, pointize = "auto",
                   block_size = NULL, auto_range = FALSE,
                   parallel = FALSE, metrics = NULL) {
  if (!inherits(data_sf, "sf")) stop("cv_gwr(): `data_sf` must be an sf object.")
  metrics <- .check_metrics_fn(metrics, "cv_gwr")
  if (!requireNamespace("sp", quietly = TRUE))
    stop("cv_gwr(): package 'sp' is required (a GWR fit keeps its local ",
         "estimates in a SpatialPointsDataFrame).", call. = FALSE)

  # match.arg() is the whole of kernel validation: it refuses a wrong case,
  # NA or a vector, so a separate repair step after it could never run.
  kernel <- match.arg(kernel)
  # An adaptive count fit_gwr_model() would refuse is refused once, here,
  # instead of in every fold (which returned "all folds failed").
  if (!is.null(bandwidth) && isTRUE(adaptive))
    .check_scalar(bandwidth, "bandwidth", "cv_gwr", min = 1,
                  max = .Machine$integer.max,
                  what = "a single number of nearest neighbours when adaptive = TRUE")

  if (!("..row_id" %in% names(data_sf))) data_sf$`..row_id` <- seq_len(nrow(data_sf))
  folds <- .folds_from_labels(folds, data_sf, "cv_gwr")
  boundary <- .cv_boundary_crs(boundary, "cv_gwr")
  dat_sf <- prep_model_data(data_sf, response_var, predictor_vars, boundary, pointize)
  if (!("..row_id" %in% names(dat_sf)))
    stop("cv_gwr(): `prep_model_data()` must preserve `..row_id`.")

  keep_idx <- dat_sf$`..row_id`

  # Refuse folds built from other data before anything is fitted: the splits
  # are row IDs, so a wrong `folds` of the right size applies silently.
  .check_fold_probe(folds, data_sf, "cv_gwr")
  .warn_fold_args_ignored(folds, block_size, auto_range, "cv_gwr")

  if (is.null(folds)) {
    message("cv_gwr(): no folds supplied -- using spatial block k-fold CV (k=", k, ").")
    folds <- make_folds(dat_sf, k = k, method = "block_kfold",
                        seed = seed,
                        boundary = .cv_boundary_crs(boundary, "cv_gwr", to = dat_sf),
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
  #
  # GWR folds never fork.  The rule was made for the engine this function
  # had: GWmodel is built with OpenMP, GNU libgomp is not fork-safe, and once
  # any GWR had run in the session a forked child's first parallel region
  # waited on a futex for ever, so parallel::mclapply() never returned.  The
  # local regressions are now solved by this package and no GWmodel code
  # runs in a fold.  The folds still run here, one after another: forked
  # folds have not been tested with the package's own solve (its batched
  # scoring is BLAS matrix products, and an OpenMP-threaded BLAS is the same
  # hazard).  A PSOCK cluster would avoid the fork, but its workers need this
  # same spatialkit installed, which pkgload::load_all() or any development
  # copy does not give them, and they cannot see what a `metrics` closure
  # reads from the global environment.
  if (.resolve_n_cores(parallel) > 1L) {
    .warn_and_log(paste0(
      "cv_gwr(): `parallel` is ignored and the %d folds run one after ",
      "another: GWR folds are not run in forked workers."),
      length(remapped_folds))
    parallel <- FALSE
  }
  # predict() on a GWR warns per call about held-out locations it cannot
  # estimate and about those beyond the reach of the fold's fit.  The folds
  # run in this process, so the per-fold warnings are collected here and each
  # class is raised once, with the totals over the folds that had any.  No
  # n_total: fold schemes differ in how many rows they predict.
  res <- .gwr_tally_predict_warnings(
    .cv_run_folds(
      dat_sf = dat_sf, response_var = response_var,
      predictor_vars = predictor_vars,
      remapped_folds = remapped_folds, keep_idx = keep_idx,
      fit_one = fit_one, fold_info_fn = fold_info_fn,
      p = NULL, parallel = parallel, seed = seed, metrics = metrics
    ),
    caller = "cv_gwr")

  preds <- if (length(res$pred_rows)) do.call(rbind, res$pred_rows) else
    data.frame(`..row_id` = integer(), fold = integer(),
               y = numeric(), yhat = numeric(), y_train_mean = numeric())
  folds_df <- if (length(res$fold_stats)) as.data.frame(dplyr::bind_rows(res$fold_stats)) else
    .empty_fold_metrics("bandwidth", metrics)

  # The folds SUPPLIED (or built), not the ones that survived .remap_folds():
  # a fold dropped there is exactly the kind of thing this count exists to
  # make visible against n_folds_succeeded.
  n_attempted <- length(if (!is.null(folds$folds)) folds$folds else folds)
  n_succeeded <- length(res$fold_stats)
  if (n_succeeded == 0L && n_attempted > 0L) {
    why <- .cv_first_error_suffix(res)
    .warn_and_log(paste0("cv_gwr(): all folds failed (all %d folds failed to ",
                         "produce predictions); cross-validation results ",
                         "contain no predictions.%s"),
                  n_attempted, why)
  } else {
    .cv_warn_failed_folds("cv_gwr", res, preds, length(keep_idx),
                          n_attempted, n_succeeded)
  }

  list(overall = .cv_overall_metrics(preds, metrics), fold_metrics = folds_df,
       predictions = preds, folds = .bare_folds(remapped_folds),
       # Reported so that a run where every fold failed is visible in the
       # return value, not only in a warning: `overall` is a well-formed
       # all-NA row either way.  cv_spatial() and cv_rf() report the same two
       # fields, so all four CV entry points share one shape.
       n_folds_attempted = n_attempted, n_folds_succeeded = n_succeeded,
       fold_status = .cv_fold_status(remapped_folds, res),
       orphan_rows = attr(remapped_folds, "orphans"),
       n_unknown_ids = attr(remapped_folds, "n_unknown_ids"),
       n_dropped = as.integer(.get_row_record(dat_sf, "dropped")$n %||% 0L),
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
#' inside the 50/80/95% intervals, and \code{mean_CRPS} rates sharpness and
#' calibration together.  That is the reason to reach for it.  A Bayesian model
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
#' @param folds Optional fold definitions, in any of three shapes: a
#'   \code{\link{make_folds}()} return value; a bare list of
#'   \code{list(train =, test =)} pairs of \code{..row_id} values; or a vector
#'   of fold labels, one per row, which becomes leave-that-label-out splits.
#'   The label vector is how folds built by another package are used here
#'   (\code{blockCV::cv_spatial()} returns one as \code{$folds_ids}), since its
#'   \code{$folds_list} holds two \emph{unnamed} vectors per fold and is
#'   refused by name.  Train and test must be disjoint: a fold that trains on
#'   its own test rows is not a cross-validation split and is refused with an
#'   error.  IDs naming no row in the prepared data are dropped with
#'   a logged count (expected when rows were removed for missing values; a sign
#'   the folds came from other data when they were not).
#' @param k Number of folds to build when \code{folds} is \code{NULL}; not
#'   read when \code{folds} is supplied.  Default 5.
#' @param seed RNG seed. Default 123.  It seeds fold construction \strong{and},
#'   through a per-fold draw, each fold's Stan sampler, so two seeds give
#'   different posteriors even on identical \code{folds}. A \code{seed} in
#'   \code{fit_args} overrides the per-fold draw with one fixed sampler seed.
#' @param boundary Optional polygonal sf/sfc for CRS alignment.
#' @param pointize Geometry coercion strategy.
#' @param fit_args Named list of extra arguments for fit_bayesian_spatial_model().
#'   A user-supplied \code{gp_k} is respected in every fold; when omitted, the
#'   GP rank is auto-selected per training fold.  \code{compute_loo},
#'   \code{boundary}, and \code{pointize} are always overridden by the CV
#'   internals.  A categorical or ordinal \code{family}
#'   (\code{brms::categorical()}, \code{cumulative}, \code{sratio},
#'   \code{cratio}, \code{acat}) is refused before anything is fitted: its
#'   prediction is a probability per response category, and every score
#'   here needs one number per row.
#' @param summary "mean" or "median" for posterior predictions.
#' @param compute_pred_intervals Logical; compute predictive intervals.
#' @param coverage_levels Numeric vector of the nominal coverage levels to
#'   score, as proportions strictly between 0 and 1 (\code{0.95}, not
#'   \code{95}), each given once; anything else is an error.  Each becomes a
#'   column \code{coverage_<percent>} of \code{fold_metrics}, named at full
#'   precision (\code{0.975} gives \code{coverage_97.5}).
#' @param block_size Optional minimum block edge length (projected CRS
#'   units) for the blocks built when \code{folds} is \code{NULL}; see
#'   \code{\link{make_folds}()}.  Beside supplied \code{folds} it has
#'   nothing to size and raises a warning (class
#'   \code{spatialkit_fold_args_ignored}).
#' @param auto_range Logical.  If \code{TRUE} and \code{folds} is \code{NULL},
#'   estimate the autocorrelation range and use it as the minimum block size.
#'   Default \code{FALSE}.  \code{TRUE} beside supplied \code{folds} raises
#'   the same warning as \code{block_size}.
#' @param parallel Logical or positive integer.  If \code{TRUE},
#'   auto-detect the number of cores and fit folds in parallel via
#'   \code{parallel::mclapply()} (macOS / Linux; falls back to sequential
#'   on Windows).  If an integer > 1, use that many cores.  Default
#'   \code{FALSE} (sequential).  Bayesian folds with full MCMC runs
#'   are the primary beneficiary of this option.  Mind the memory: every
#'   fold compiles its own Stan model, and one compilation can take several
#'   GB (3.6 GB was measured), so \code{parallel = n} runs \code{n} of them
#'   at once.  A compiler killed for lack of memory fails its fold with
#'   rstan's \code{"invalid connection"} error, which \code{fold_status}
#'   records; use fewer cores if you see it.
#' @param metrics Optional scoring function of your own, a
#'   \code{function(y, yhat)} returning a named numeric vector; its names
#'   become columns of \code{fold_metrics} (per fold) and \code{overall}
#'   (pooled) beside the built-in ones.  See \strong{Your own metrics} on
#'   \code{\link{cv_spatial}()} for the contract.  \code{yhat} is the
#'   posterior predictive \code{summary} (mean or median) of each held-out
#'   row; a score that needs the draws belongs in \code{predictive_coverage}
#'   and \code{CRPS}, which this function computes itself.
#' @section Warnings from the folds:
#' \code{\link{fit_bayesian_spatial_model}()} warns about a fit under four
#' condition classes: \code{spatialkit_gp_basis_size},
#' \code{spatialkit_gp_resolution}, \code{spatialkit_convergence} and
#' \code{spatialkit_gp_adequacy} (see its \strong{Warnings}).  Each fold is
#' a fit, and none of those warnings says which fold it came from, so they
#' are collected as the folds run and each class is raised once after the
#' last fold, under this function's name and with the same class.  Each
#' warning names the folds and carries \code{folds} (their labels) and
#' \code{fold_messages} (every fold's own text, named by fold) as fields of
#' the condition.  The convergence warning covers the folds that were
#' scored and gives the checks each failed, also as the field
#' \code{failed}; the other three quote the first fold's own text.  Each
#' fold's own line is in the log.  Warnings that Stan or \pkg{brms} raise
#' themselves are not collected and are passed on.
#' @inheritSection model_metrics Percentage errors on responses with zeros
#' @inheritSection model_metrics Which metrics survive a non-Gaussian response
#' @return A list with \code{overall}, \code{fold_metrics},
#'   \code{predictions}, \code{folds}, \code{n_folds_attempted},
#'   \code{n_folds_succeeded}, \code{fold_status}, \code{orphan_rows},
#'   \code{n_unknown_ids}, \code{n_dropped}, \code{formula},
#'   \code{predictive_coverage} and \code{coverage_levels} (the nominal
#'   levels, named by their \code{coverage_*} column).
#'   The two fold counts make a run where every fold failed visible in the
#'   return value itself, beyond the warning, and \code{fold_status} (one
#'   row per fold: \code{fold}, \code{status}, \code{message}) keeps the
#'   reason each missing fold is missing, including the error text of a fold
#'   whose sampler failed, where a long run's console output would not.
#'   See \code{\link{cv_spatial}} for the five statuses and for
#'   \code{orphan_rows}.  \code{predictions} carries,
#'   beyond the columns its siblings share, \code{yhat_sd}: the posterior
#'   predictive standard deviation of each held-out row, from the same draws
#'   that give the coverage below (\code{NA} when
#'   \code{compute_pred_intervals = FALSE} or the draws failed for that fold).
#'   \code{overall$Adj_R2} is always \code{NA}, as for every \code{cv_*()}:
#'   see \code{\link{cv_spatial}}.
#'   \code{fold_metrics} carries, beyond the columns its siblings share,
#'   \code{gp_k}, \code{gp_n_basis}, \code{n_draws}, \code{CRPS}, the
#'   \code{coverage_*} columns, \code{convergence_ok} and
#'   \code{gp_adequacy}.  \code{convergence_ok} is \code{TRUE} or
#'   \code{FALSE} as \code{\link{fit_bayesian_spatial_model}()} judged that
#'   fold's sampler (divergent transitions, R-hat, effective sample size and
#'   E-BFMI; see its \code{check_convergence}), \code{NA}
#'   when \code{fit_args} sets \code{check_convergence = FALSE}.  A fold that
#'   did not converge is scored like the others, so a run with any
#'   \code{FALSE} raises one warning of class \code{spatialkit_convergence}
#'   naming those folds and the checks each failed.
#'   \code{gp_adequacy} is what that fold's posterior said about its
#'   Gaussian-process basis (\code{"ok"}, \code{"below_resolution"},
#'   \code{"unresolved_short_range"}, \code{"no_signal_above_floor"} or
#'   \code{"not_checked"}; see "What the basis resolves" in
#'   \code{\link{fit_bayesian_spatial_model}}).  A fold whose basis was
#'   found too coarse is scored like the others too.  See \strong{Warnings
#'   from the folds}.
#'   The \code{predictive_coverage} entries (one per
#'   \code{coverage_levels} value, plus \code{mean_CRPS}) are averages across
#'   folds \strong{weighted by each fold's \code{n_pred}}, because the per-fold
#'   values in \code{fold_metrics} are themselves means over that fold's test
#'   rows; an unweighted average would not be the pooled quantity when fold
#'   sizes differ, which for spatially blocked folds they routinely do.
#' @family cross-validation
#' @examples
#' \dontrun{
#' # Not run: fits with Stan, which needs a working C++ toolchain and takes
#' # minutes of MCMC -- both outside what an example may assume.
#' if (requireNamespace("brms", quietly = TRUE)) {
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
#'   # Two short chains per fold keep this to a few minutes and leave the
#'   # posterior rough, so the intervals below will run narrow; a run to
#'   # report uses brms's defaults (chains = 4, iter = 2000).
#'   cv <- cv_bayes(dat, "price", "elev", k = 2,
#'                  fit_args = list(chains = 2, iter = 1000))
#'   print(cv$overall)
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
                     parallel = FALSE, metrics = NULL) {
  summary <- match.arg(summary)
  if (!inherits(data_sf, "sf")) stop("cv_bayes(): `data_sf` must be an sf object.")
  # A categorical or ordinal family predicts a probability per response
  # category, and every score here needs one number per row: predict(type =
  # "epred") stops for these families, so each fold compiled and sampled a
  # full model and was then discarded (k = 2 on 50 rows: 4.35 minutes for an
  # empty result).  Refused before anything is fitted, as
  # fit_bayesian_spatial_model() refuses a factor under bernoulli().
  fam <- .brms_family_name(if (is.list(fit_args)) fit_args$family)
  if (!is.na(fam) && fam %in% .brms_category_families)
    stop(sprintf(paste0("cv_bayes(): the %s family gives a probability per ",
                        "response category, not one predicted number per row, ",
                        "and cross-validation here scores one number per row ",
                        "(RMSE, CRPS, interval coverage). Every fold would run ",
                        "its MCMC and then fail to be scored, so nothing is ",
                        "fitted. For a two-level outcome, code it as 0/1 and use ",
                        "brms::bernoulli()."), sQuote(fam)), call. = FALSE)
  metrics <- .check_metrics_fn(metrics, "cv_bayes")
  # Named by column: the nominal level travels with the result from here on,
  # instead of being read back from a column name.
  cov_lv  <- .check_coverage_levels(coverage_levels)

  if (!("..row_id" %in% names(data_sf))) data_sf$`..row_id` <- seq_len(nrow(data_sf))
  folds <- .folds_from_labels(folds, data_sf, "cv_bayes")
  boundary <- .cv_boundary_crs(boundary, "cv_bayes")
  dat_sf <- prep_model_data(data_sf, response_var, predictor_vars, boundary, pointize)
  if (!("..row_id" %in% names(dat_sf)))
    stop("cv_bayes(): `prep_model_data()` must preserve `..row_id`.")

  keep_idx <- dat_sf$`..row_id`
  n_pred   <- length(predictor_vars)

  # Refuse folds built from other data before anything is fitted: the splits
  # are row IDs, so a wrong `folds` of the right size applies silently.
  .check_fold_probe(folds, data_sf, "cv_bayes")
  .warn_fold_args_ignored(folds, block_size, auto_range, "cv_bayes")

  if (is.null(folds)) {
    message("cv_bayes(): no folds supplied -- using spatial block k-fold CV (k=", k, ").")
    folds <- make_folds(dat_sf, k = k, method = "block_kfold",
                        seed = seed,
                        boundary = .cv_boundary_crs(boundary, "cv_bayes", to = dat_sf),
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
    # `seed` reaches the SAMPLER, as cv_rf()'s reaches the forest.  It did
    # not: fit_bayesian_spatial_model() carries seed = 123 and fold_fit_args
    # never set it, so every fold of every run sampled from Stan seed 123 and
    # cv_bayes(seed = ) changed nothing on fixed folds.  .cv_run_folds() runs
    # each fold under its own seeded stream, so a draw here is distinct and
    # reproducible per (seed, fold).  A seed in fit_args still wins.
    fold_fit_args <- modifyList(fit_args, list(
      compute_loo = FALSE, boundary = boundary,
      pointize = pointize,
      .already_prepped = TRUE
    ))
    if (is.null(fit_args$seed))
      fold_fit_args$seed <- sample.int(.Machine$integer.max, 1L)
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
      CRPS    = NA_real_,
      # TRUE / FALSE as fit_bayesian_spatial_model() judged the fold's
      # sampler (divergences, R-hat, ESS, E-BFMI), NA when nothing was
      # checked.  A fold whose posterior did not converge scores like any
      # other, so the table has to say which ones those are.
      convergence_ok = {
        v <- fit_obj$info$convergence_ok
        if (is.logical(v) && length(v) == 1L) v else NA
      },
      # What the fold's posterior said about its GP basis
      # ($info$gp_adequacy$status): a fold whose basis was too coarse scores
      # like any other.
      gp_adequacy = {
        v <- fit_obj$info$gp_adequacy$status
        if (is.character(v) && length(v) == 1L) v else NA_character_
      }
    )
    # Pre-initialise coverage columns so every fold emits the same schema
    # even when the posterior-draw step fails for some folds; heterogeneous
    # per-fold columns would otherwise break row-binding of fold_metrics.
    for (cn in names(cov_lv)) extras[[cn]] <- NA_real_

    # Full posterior predictive draws for intervals + CRPS
    if (isTRUE(compute_pred_intervals)) {
      ppred_draws <- try(
        predict(fit_obj, newdata = test_sf, type = "predict", draws = TRUE),
        silent = TRUE
      )
      if (inherits(ppred_draws, "try-error"))
        .log_warn(paste0("cv_bayes(): the posterior predictive draws failed on a ",
                         "fold (%s); its coverage and CRPS are NA."),
                  .try_error_message(ppred_draws))
      if (!inherits(ppred_draws, "try-error") && is.matrix(ppred_draws)) {
        extras$n_draws <- nrow(ppred_draws)
        # Per-row posterior predictive SD, for $predictions$yhat_sd.  The
        # column existed as an unconditional NA placeholder; the draws that
        # give the coverage below give this for free.
        if (nrow(ppred_draws) >= 2L && ncol(ppred_draws) == length(y_true))
          extras$..per_row <- data.frame(
            yhat_sd = apply(ppred_draws, 2L, stats::sd),
            stringsAsFactors = FALSE)

        # Coverage at each level, and the empirical CRPS via the NRG (energy)
        # form (see .crps_energy()).  In a try of their own, so that one
        # failure here -- quantile() refuses a draw matrix holding an NA --
        # costs these columns and not gp_k, n_draws and yhat_sd with them.
        scores <- try({
          cov <- vapply(cov_lv, function(cl) {
            alpha <- (1 - cl) / 2
            lwr <- apply(ppred_draws, 2L, stats::quantile, probs = alpha)
            upr <- apply(ppred_draws, 2L, stats::quantile, probs = 1 - alpha)
            mean(y_true >= lwr & y_true <= upr, na.rm = TRUE)
          }, numeric(1))
          list(cov = cov,
               crps = mean(.crps_energy(ppred_draws, y_true), na.rm = TRUE))
        }, silent = TRUE)
        if (inherits(scores, "try-error")) {
          .log_warn(paste0("cv_bayes(): coverage and CRPS could not be computed ",
                           "from a fold's draws (%s); they are NA there."),
                    .try_error_message(scores))
        } else {
          for (cn in names(cov_lv)) extras[[cn]] <- scores$cov[[cn]]
          extras$CRPS <- scores$crps
        }
      }
    }

    extras
  }

  # Each fold's fit can raise three classed warnings of its own: a sampler
  # that did not converge, a variogram range below what the basis resolves,
  # and a basis the fitted model shows to be too coarse.  None says which
  # fold it is, and passed on they arrived once per fold beside this
  # function's own summary of the same folds.  They are collected here with
  # their fold and raised once per class after the folds, naming the folds.
  collected <- .cv_collect_fold_warnings(.cv_run_folds(
    dat_sf = dat_sf, response_var = response_var,
    predictor_vars = predictor_vars,
    remapped_folds = remapped_folds, keep_idx = keep_idx,
    fit_one = fit_one, fold_info_fn = fold_info_fn,
    predict_args = list(summary = summary),
    p = NULL, parallel = parallel, seed = seed, metrics = metrics
  ))
  res <- collected$value

  # Every fold carries a yhat_sd column: the posterior predictive SD where
  # the draws succeeded, NA where compute_pred_intervals = FALSE or the
  # draw step failed, so the rows bind whatever happened per fold.
  pred_rows_with_sd <- lapply(res$pred_rows, function(pr) {
    if (is.null(pr$yhat_sd)) pr$yhat_sd <- NA_real_
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
               n_MAPE = integer(), n_SMAPE = integer(),
               CRPS = numeric(), gp_k = integer(),
               gp_n_basis = integer(), n_draws = integer(),
               convergence_ok = logical(), gp_adequacy = character())
  if (!nrow(folds_df))
    for (cn in .user_metric_names(metrics)) folds_df[[cn]] <- numeric()

  # The folds SUPPLIED (or built), not the ones that survived .remap_folds():
  # a fold dropped there is exactly the kind of thing this count exists to
  # make visible against n_folds_succeeded.
  n_attempted <- length(if (!is.null(folds$folds)) folds$folds else folds)
  n_succeeded <- length(res$fold_stats)
  if (n_succeeded == 0L && n_attempted > 0L) {
    why <- .cv_first_error_suffix(res)
    .warn_and_log(paste0("cv_bayes(): all folds failed (all %d folds failed to ",
                         "produce predictions); cross-validation results ",
                         "contain no predictions.%s"),
                  n_attempted, why)
  } else {
    .cv_warn_failed_folds("cv_bayes", res, preds, length(keep_idx),
                          n_attempted, n_succeeded)
  }

  # A fold whose sampler did not converge is scored like any other: its
  # RMSE, CRPS and coverage enter `overall`.  One R warning per run names
  # those folds and the checks each failed; fold_metrics$convergence_ok marks
  # them.  The two warnings about the basis are raised the same way.
  .cv_raise_fold_warnings(
    "cv_bayes", collected$said, scored = folds_df$fold,
    n_run = length(remapped_folds),
    not_converged = if (nrow(folds_df) && "convergence_ok" %in% names(folds_df))
      folds_df$fold[folds_df$convergence_ok %in% FALSE] else integer(0))

  list(overall = .cv_overall_metrics(preds, metrics), fold_metrics = folds_df,
       predictions = preds, folds = .bare_folds(remapped_folds),
       # Reported so that a run where every fold failed is visible in the
       # return value, not only in a warning -- which matters most here,
       # where a fold dies on Stan sampling rather than on bad input.
       n_folds_attempted = n_attempted, n_folds_succeeded = n_succeeded,
       fold_status = .cv_fold_status(remapped_folds, res),
       orphan_rows = attr(remapped_folds, "orphans"),
       n_unknown_ids = attr(remapped_folds, "n_unknown_ids"),
       n_dropped = as.integer(.get_row_record(dat_sf, "dropped")$n %||% 0L),
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
         cov_cols  <- intersect(names(cov_lv), names(folds_df))
         cov_means <- vapply(cov_cols, function(cn) wm(folds_df[[cn]]), numeric(1))
         c(as.list(cov_means), mean_CRPS = wm(folds_df$CRPS))
       } else NULL,
       # The nominal level behind each coverage column, named by column.  The
       # names used to be the only record, rounded to a whole percent: 0.995
       # was drawn at 1.00 by plot_calibration(), and 0.975 and 0.985 were
       # both "coverage_98", one overwriting the other.
       coverage_levels = cov_lv)
}


#' Validate cv_bayes()'s coverage_levels and name their columns
#'
#' Every level must lie strictly between 0 and 1.  The column is named from
#' the level in percent at full precision (\code{coverage_95},
#' \code{coverage_97.5}), so no two distinct levels share a column; a level
#' given twice is an error.
#'
#' @param x The \code{coverage_levels} argument.
#' @return A numeric vector of the levels, named by column.
#' @keywords internal
#' @noRd
.check_coverage_levels <- function(x) {
  if (is.null(x) || (is.numeric(x) && !length(x)))
    return(stats::setNames(numeric(0), character(0)))
  if (!is.numeric(x) || any(!is.finite(x)) || any(x <= 0 | x >= 1)) {
    hint <- if (is.numeric(x) && all(is.finite(x)) && all(x > 1 & x < 100))
      sprintf(" For percentages, divide by 100: c(%s).",
              paste(as.character(x / 100), collapse = ", ")) else ""
    stop(sprintf(paste0("cv_bayes(): `coverage_levels` must be numbers strictly ",
                        "between 0 and 1 (0.95 for a 95%% interval); got %s.%s"),
                 paste(as.character(x), collapse = ", "), hint),
         call. = FALSE)
  }
  cols <- paste0("coverage_", vapply(as.numeric(x), function(cl)
    format(round(cl * 100, 6), scientific = FALSE, trim = TRUE, digits = 15),
    character(1)))
  if (anyDuplicated(cols))
    stop(sprintf("cv_bayes(): `coverage_levels` gives the level %s more than once.",
                 paste(unique(as.character(x[duplicated(cols)])),
                       collapse = ", ")), call. = FALSE)
  stats::setNames(as.numeric(x), cols)
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
#' @section Name collision with blockCV:
#' The \pkg{blockCV} package exports a function of the same name that does
#' the opposite job: \code{blockCV::cv_spatial()} \emph{builds} spatial
#' folds, where this function \emph{runs} a cross-validation over folds it is
#' given.  With both packages attached, whichever was attached last masks the
#' other; \code{spatialkit::cv_spatial()} always resolves to this one.  The
#' two cooperate: \code{blockCV::cv_spatial()} returns
#' its fold assignment as \code{$folds_ids}, a vector of fold labels, and
#' that vector is accepted directly as the \code{folds} argument here and in
#' every other \code{cv_*()} function.  Fold construction is
#' \pkg{blockCV}'s home ground and this package does not try to match its
#' breadth there; what this package adds is the path from irregular points
#' through data-drawn regions to a cross-validated, compared model.
#'
#' @param data_sf An sf object.
#' @param response_var Response column name.
#' @param predictor_vars Predictor column names.
#' @param fit_fn A function of one argument, the training slice of
#'   \code{data_sf}, returning a \code{spatial_fit} built with
#'   \code{\link{new_spatial_fit}()}.  It is called once per fold on the
#'   training rows only, so anything done inside it (scaling, tuning, an inner
#'   variable sweep) is already nested and leak-free.  The \code{subclass} it
#'   stamps must have a \code{predict.<subclass>()} method registered, because
#'   that is how each fold is scored.
#' @param folds Optional fold definitions, in any of three shapes: a
#'   \code{\link{make_folds}()} return value; a bare list of
#'   \code{list(train =, test =)} pairs of \code{..row_id} values; or a vector
#'   of fold labels, one per row, which becomes leave-that-label-out splits.
#'   The label vector is how folds built by another package are used here
#'   (\code{blockCV::cv_spatial()} returns one as \code{$folds_ids}), since its
#'   \code{$folds_list} holds two \emph{unnamed} vectors per fold and is
#'   refused by name.  Train and test must be disjoint: a fold that trains on
#'   its own test rows is not a cross-validation split and is refused with an
#'   error.  IDs naming no row in the prepared data are dropped with
#'   a logged count.
#'   Built via \code{block_kfold} when \code{NULL}.
#' @param k Number of folds to build when \code{folds} is \code{NULL}; not
#'   read when \code{folds} is supplied.
#' @param seed RNG seed.
#' @param boundary Optional boundary for fold construction.
#' @param pointize Geometry coercion strategy.
#' @param predict_args Extra arguments for predict().
#' @param fold_info_fn Optional \code{function(fit, test_sf, y, yhat)}
#'   returning a named list (or a named vector) of per-fold extras (a
#'   bandwidth, a tuning value, anything read off the fitted object), added
#'   as columns of \code{fold_metrics}.  It sees the fit and the held-out
#'   layer, which
#'   \code{metrics} does not; it is applied per fold only, and its values
#'   are not pooled.  An element \code{..per_row} that is a data frame with
#'   one row per held-out observation is spliced into \code{predictions}
#'   instead (\code{NA} in the rows of a fold that returned none; one of the
#'   wrong length is dropped and logged); its columns must be named, once,
#'   and not reuse a column \code{predictions} already has (\code{..row_id},
#'   \code{fold}, \code{y}, \code{yhat}, \code{y_train_mean}).  Every
#'   other element must be named, hold one value, and not reuse a column
#'   \code{fold_metrics} already has;
#'   anything else is an error.  A \code{fold_info_fn} that throws on a fold
#'   is logged, its columns are \code{NA} for that fold, and
#'   \code{fold_status$message} says so; the fold is kept.
#' @param p Number of predictors for Adj R² (NULL to skip).  Only meaningful
#'   for models with a fixed global parameter count; pass NULL for models
#'   with spatially varying coefficients (e.g. GWR).
#' @param block_size Optional minimum block edge length (projected CRS
#'   units) for the blocks built when \code{folds} is \code{NULL}; see
#'   \code{\link{make_folds}()}.  Beside supplied \code{folds} it has
#'   nothing to size and raises a warning (class
#'   \code{spatialkit_fold_args_ignored}).
#' @param auto_range Logical.  If \code{TRUE} and \code{folds} is \code{NULL},
#'   estimate the autocorrelation range and use it as the minimum block size.
#'   Default \code{FALSE}.  \code{TRUE} beside supplied \code{folds} raises
#'   the same warning as \code{block_size}.
#' @param parallel Logical or positive integer.  If \code{TRUE},
#'   auto-detect the number of cores and fit folds in parallel via
#'   \code{parallel::mclapply()} (macOS / Linux; falls back to sequential
#'   on Windows).  If an integer > 1, use that many cores.  Default
#'   \code{FALSE} (sequential).  A learner that runs OpenMP code (a
#'   function of GWmodel, or an xgboost built with GNU libgomp) can hang the
#'   forked workers once it has run in the session; keep such a
#'   \code{fit_fn} sequential.
#' @param metrics Optional scoring function of your own; see \strong{Your
#'   own metrics} below.  Default \code{NULL}: the built-in metrics only.
#' @param .caller Internal. The name the messages carry, so a wrapper such as
#'   \code{\link{cv_rf}} reports itself in place of \code{cv_spatial()}.
#' @section Your own metrics:
#' The built-in columns (\code{RMSE}, \code{MAE}, \code{MAPE},
#' \code{SMAPE}, \code{R2}, \code{Adj_R2}) are the Gaussian regression
#' set, and the section below says which of them survive a count or a
#' bounded response.  \code{metrics} is the way to score what they cannot: a
#' \code{function(y, yhat)} that returns a named numeric vector (a named
#' list of scalars, or a one-row data frame, also serve), for example a
#' Poisson deviance, a log score on a probability, or a loss with your own
#' weights.  It is applied twice, in the same way the built-in metrics are:
#' once per fold, to that fold's held-out rows, so each name becomes a column
#' of \code{fold_metrics}; and once to the pooled out-of-sample predictions
#' of every fold, so each name becomes a column of \code{overall}.  Only the
#' pairs the built-in metrics use reach the function (both \code{y} and
#' \code{yhat} finite), so its columns describe the same rows as
#' \code{RMSE}, and \code{n_pred} counts them.
#'
#' The contract: every element named, names unique and not one of the
#' built-in column names, one number per name.  The built-in names include
#' the per-fold extras of the backend or of your \code{fold_info_fn}
#' (\code{bandwidth} for \code{cv_gwr()}; \code{CRPS}, \code{coverage_*},
#' \code{gp_k}, \code{gp_n_basis} and \code{n_draws} for \code{cv_bayes()}),
#' and \code{mean_CRPS}, which \code{compare_models_cv()} writes.  Anything
#' else is an error,
#' because a scoring function that returns the wrong shape is a mistake to
#' surface instead of a fold to skip.  A function that \emph{throws} on a
#' fold is logged and its columns are \code{NA} for that fold (and for
#' \code{overall}, if it throws on the pooled predictions); a fold is never
#' dropped for it.  When no fold produced a prediction the empty
#' \code{fold_metrics} frame still carries the function's columns, typed,
#' provided the function can be called on zero-length input.
#'
#' \code{fold_info_fn} is the per-fold half of the same mechanism, with
#' access to the fitted object and the held-out layer; \code{metrics} sees
#' only the two vectors but is also pooled.  \code{compare_models_cv()}
#' hands one \code{metrics} to every backend, so the columns are comparable
#' across the rows of its \code{overall}.
#' @section Warnings a learner raises in every fold:
#' A warning from a fold's fit or from its \code{predict()} reaches the
#' caller as it is raised, so a learner that warns in every fold warns once
#' per fold.  Two groups of the package's own warnings are collected over the
#' folds instead and raised once per call, under this function's name and
#' with the class they had:
#' \itemize{
#'   \item the two in which \code{predict()} on a GWR counts locations
#'     (\code{"spatialkit_gwr_not_estimable"} and
#'     \code{"spatialkit_gwr_extrapolation"}), with the totals over the folds
#'     that had any, as \code{\link{cv_gwr}()} raises them (see
#'     \code{\link{predict.gwr_fit}});
#'   \item the four \code{\link{fit_bayesian_spatial_model}()} raises about
#'     a fit, naming the folds, as \code{\link{cv_bayes}()} raises them (see
#'     \strong{Warnings from the folds} there).
#' }
#' With \code{parallel} above one core a fold's warnings are raised after
#' the folds have run, once per distinct text; one of the package's own
#' classes keeps its class and its fields.
#' @inheritSection model_metrics Percentage errors on responses with zeros
#' @return A list with \code{overall}, \code{fold_metrics}, \code{predictions},
#'   \code{folds}, the two fold counts \code{n_folds_attempted} and
#'   \code{n_folds_succeeded}, and four elements that say what happened to
#'   the difference between them and to the rows: \code{fold_status},
#'   \code{orphan_rows}, \code{n_unknown_ids} and \code{n_dropped}.  The
#'   counts are reported deliberately: a
#'   \code{fit_fn} that fails on every fold otherwise looks like a successful
#'   run that happened to score \code{NA}, so compare them before trusting
#'   \code{overall}.  \code{fold_status} is a data.frame with one row per
#'   fold supplied (\code{fold}, \code{status} and \code{message}), where
#'   \code{status} is \code{"ok"} (\code{message} is empty unless
#'   \code{fold_info_fn} failed on the fold); \code{"error"} (the fit or its
#'   \code{predict()} threw; \code{message} is the error text);
#'   \code{"skipped"} (nothing scorable: too few matched rows, a prediction
#'   of the wrong length, or no finite observed/predicted pair);
#'   \code{"dropped"} (an empty test set, or fewer than two training rows,
#'   once incomplete rows were removed, so the fold never reached the fitter);
#'   or \code{"worker_error"} (a parallel worker died, for example killed for
#'   lack of memory).  Each fold runs in its own worker, so a failure costs
#'   that fold only, and an error that stops a sequential run (a
#'   \code{metrics} or \code{fold_info_fn} return value of the wrong shape)
#'   stops a parallel one too, naming the fold.  Every fold missing
#'   from \code{fold_metrics} has its reason there, which matters most when
#'   the console output of a long run is gone.  When some folds, but not
#'   all, end as \code{"error"}, \code{"skipped"} or \code{"worker_error"},
#'   the function warns, naming them and how many rows \code{overall}
#'   covers: it is pooled over the folds that produced predictions, and the
#'   fold that fails is often the hardest to predict, so it may flatter the
#'   model.  (A \code{"dropped"} fold has its own warning.)
#'   \code{orphan_rows} holds the
#'   \code{..row_id}s of rows in the data that no fold names (they enter no
#'   training set and are never scored; non-empty only when the folds were
#'   built on a different or subsetted layer), and \code{n_unknown_ids}
#'   counts the distinct row IDs the folds name that the data does not have.
#'   Each such row is named by every fold, once as a test row and once in
#'   each other fold's training set, and this counts the row, not the
#'   mentions (expected when rows were removed for missing values).
#'   \code{n_dropped} is how many rows \code{prep_model_data()} removed for
#'   missing or non-finite values or a bad geometry before any fold was
#'   fitted.  The \code{fold} column of
#'   \code{fold_metrics}, \code{predictions} and \code{fold_status} carries
#'   the fold's index in the \code{folds} object that was supplied, so it
#'   lines up with \code{make_folds()$assignment$fold} even when some folds
#'   were unusable and dropped.  Splits that already carry a \code{fold_id},
#'   as the \code{folds} of a \code{cv_*()} result do, keep it: handing
#'   one run's \code{folds} to another labels each fold as the first run
#'   and \code{\link{fold_separation}()} do, a dropped fold's gap
#'   included.  \code{R2} is out-of-sample \eqn{R^2}: the
#'   total sum of squares is taken about the mean of the \emph{training} rows
#'   (in \code{overall}, each held-out row about its own fold's training
#'   mean), the null prediction available when the fold is predicted, and
#'   not about the held-out rows' own mean, a null model that would know the
#'   test data.  On spatial blocks the two can differ widely; \code{R2} is
#'   below 0 when the model predicts worse than the training mean.
#'   \code{\link{model_metrics}(newdata = )} uses the same baseline.
#'   \code{overall$Adj_R2} is always \code{NA}: the
#'   pooled out-of-sample predictions come from \code{k} separately fitted
#'   models and have no single parameter count to adjust for.  The per-fold
#'   \code{fold_metrics$Adj_R2} carries the adjusted value when \code{p} is
#'   supplied, and is \code{NA} otherwise.
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
#'   data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
#'              elev = rnorm(n)),
#'   coords = c("x", "y"), crs = 32632
#' )
#' site$price <- 10 + 0.01 * (st_coordinates(site)[, 1] - 5e5) +
#'   2 * site$elev + rnorm(n)
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
#' # Compare these before trusting the metrics above; fold_status says why
#' # any fold is missing from fold_metrics.
#' c(attempted = cv$n_folds_attempted, succeeded = cv$n_folds_succeeded)
#' cv$fold_status
#'
#' # 3. A metric of your own beside the built-in ones: per fold and pooled.
#' med_ae <- function(y, yhat) c(MedAE = stats::median(abs(y - yhat)))
#' cv2 <- cv_spatial(site, "price", "elev", fit_fn = lm_fit, k = 3, seed = 1,
#'                   metrics = med_ae)
#' cv2$overall$MedAE
#' cv2$fold_metrics[, c("fold", "RMSE", "MedAE")]
#' @export
cv_spatial <- function(data_sf, response_var, predictor_vars,
                       fit_fn, folds = NULL, k = 5, seed = 123,
                       boundary = NULL, pointize = "auto",
                       predict_args = list(), fold_info_fn = NULL,
                       p = NULL, block_size = NULL,
                       auto_range = FALSE, parallel = FALSE,
                       metrics = NULL, .caller = "cv_spatial") {
  if (!inherits(data_sf, "sf")) stop(.caller, "(): `data_sf` must be an sf object.", call. = FALSE)
  if (!is.function(fit_fn)) stop(.caller, "(): `fit_fn` must be a function.", call. = FALSE)
  metrics <- .check_metrics_fn(metrics, .caller)

  if (!("..row_id" %in% names(data_sf))) data_sf$`..row_id` <- seq_len(nrow(data_sf))
  folds <- .folds_from_labels(folds, data_sf, .caller)
  boundary <- .cv_boundary_crs(boundary, .caller)
  dat_sf <- prep_model_data(data_sf, response_var, predictor_vars, boundary, pointize)
  keep_idx <- dat_sf$`..row_id`

  # Refuse folds built from other data before anything is fitted: the splits
  # are row IDs, so a wrong `folds` of the right size applies silently.
  .check_fold_probe(folds, data_sf, .caller)
  .warn_fold_args_ignored(folds, block_size, auto_range, .caller)

  if (is.null(folds)) {
    message(.caller, "(): no folds supplied -- using spatial block k-fold CV (k=", k, ").")
    folds <- make_folds(dat_sf, k = k, method = "block_kfold",
                        seed = seed,
                        boundary = .cv_boundary_crs(boundary, .caller, to = dat_sf),
                        block_size = block_size, auto_range = auto_range,
                        response_var = response_var,
                        predictor_vars = predictor_vars)
  }

  remapped_folds <- .remap_folds(folds, keep_idx, k, seed)

  cleanup_cv <- .with_seed(seed)
  on.exit(cleanup_cv(), add = TRUE)

  # A `fit_fn` that returns a GWR is predicted fold by fold like cv_gwr()'s,
  # and predict() on a GWR counts per call the held-out locations it cannot
  # estimate and those beyond the reach of the fold's fit.  They are
  # collected here as cv_gwr() collects them, and each class is raised once
  # with the totals.  No other backend signals either class, so for every
  # other learner this does nothing.
  # A `fit_fn` that returns a Bayesian fit raises that fit's classed warnings
  # in every fold (see cv_bayes()); they are collected with their fold and
  # raised once per class below.
  collected <- .cv_collect_fold_warnings(.gwr_tally_predict_warnings(
    .cv_run_folds(
      dat_sf = dat_sf, response_var = response_var,
      predictor_vars = predictor_vars,
      remapped_folds = remapped_folds, keep_idx = keep_idx,
      fit_one = fit_fn, fold_info_fn = fold_info_fn,
      predict_args = predict_args, p = p,
      parallel = parallel, seed = seed, metrics = metrics
    ),
    caller = .caller))
  res <- collected$value

  # bind_rows(), as for fold_stats below: a fold whose fold_info_fn returned
  # no `..per_row` (conditionally, with the wrong row count, or because it
  # threw) has fewer columns, and rbind() died on that with "numbers of
  # columns of arguments do not match" after every fold had been fitted.
  # The missing cells are NA.
  preds <- if (length(res$pred_rows)) as.data.frame(dplyr::bind_rows(res$pred_rows)) else
    data.frame(`..row_id` = integer(), fold = integer(),
               y = numeric(), yhat = numeric(), y_train_mean = numeric())
  # Typed even when empty: cv_gwr() and cv_bayes() return a 0-row frame with
  # the metric columns, and a bare data.frame() here made
  # subset(fold_metrics, RMSE < 5) error for two of the four cv_*().
  folds_df <- if (length(res$fold_stats)) as.data.frame(dplyr::bind_rows(res$fold_stats)) else
    .empty_fold_metrics(metrics = metrics)

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
    .warn_and_log(paste0("%s(): all folds failed (all %d folds failed to ",
                         "produce predictions); cross-validation results ",
                         "contain no predictions.%s"),
                  .caller, n_attempted, why)
  } else {
    .cv_warn_failed_folds(.caller, res, preds, length(keep_idx),
                          n_attempted, n_succeeded)
  }
  if (length(collected$said))
    .cv_raise_fold_warnings(.caller, collected$said, scored = folds_df$fold,
                            n_run = length(remapped_folds))

  list(overall = .cv_overall_metrics(preds, metrics),
       fold_metrics = folds_df, predictions = preds,
       folds = .bare_folds(remapped_folds),
       n_folds_attempted = n_attempted, n_folds_succeeded = n_succeeded,
       fold_status = .cv_fold_status(remapped_folds, res),
       orphan_rows = attr(remapped_folds, "orphans"),
       n_unknown_ids = attr(remapped_folds, "n_unknown_ids"),
       n_dropped = as.integer(.get_row_record(dat_sf, "dropped")$n %||% 0L))
}
