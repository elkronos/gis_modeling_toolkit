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
      "cross-validation: no fold specification provided; falling back to random k-fold CV (k=%d). Random folds leak spatial autocorrelation and overstate out-of-sample performance.",
      k
    )
    warning(
      "cross-validation: falling back to random k-fold CV. For spatial data, use make_folds(method='block_kfold') to avoid optimistic performance estimates.",
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
                   j, length(ov),
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
      fold_id = j
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
#' of them by name.
#' @keywords internal
#' @noRd
.cv_reserved_metric_cols <- c("fold", "n_train", "n_test", "n_pred", "RMSE",
                              "MAE", "MAPE", "SMAPE", "R2", "Adj_R2", "n_MAPE",
                              "n_SMAPE", "model")

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
#' installed.  The word "brms" (or "GWmodel") never appears, even though
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
  f <- factor(as.character(folds),
              levels = sort(unique(as.character(folds)), method = "radix"))
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
  # same reason -- overwhelmingly the "package 'brms'/'GWmodel' is required"
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

  # Training-set mean: the correct null-model baseline for out-of-sample
  # R².  Using the test-set mean instead would give the null model credit
  # for knowing information that was not available at prediction time,
  # systematically inflating CV R².
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
  per_row <- NULL
  if (!is.null(fold_info_fn)) {
    extra <- try(fold_info_fn(fit_obj, test_sf, y_true, y_hat), silent = TRUE)
    if (!inherits(extra, "try-error") && is.list(extra)) {
      if (is.data.frame(extra$..per_row) &&
          nrow(extra$..per_row) == length(y_true))
        per_row <- extra$..per_row
      extra$..per_row <- NULL
      for (cn in names(extra)) fs[[cn]] <- extra[[cn]]
    }
  }

  # The user's scoring function, on the same finite pairs the built-in
  # metrics used.  Its pooled counterpart is applied in .cv_overall_metrics().
  if (!is.null(metrics)) {
    um <- .apply_user_metrics(metrics, y_true, y_hat,
                              where = sprintf("fold %s", format(fold_lab)))
    if (is.null(um)) um <- .na_user_metrics(metrics)
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
#' \code{cv_gwr()} never asks for it: GWmodel's OpenMP code deadlocks in a
#' forked worker (see the comment there).
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
    # the worker and re-raise in the parent, once per distinct message.
    caught_worker <- function(i) {
      msgs <- character(0)
      res  <- withCallingHandlers(
        fold_worker(i),
        warning = function(w) {
          msgs <<- c(msgs, conditionMessage(w))
          invokeRestart("muffleWarning")
        })
      # .cv_fit_one_fold() returns NULL for an unusable fold; NULL cannot
      # carry an attribute, so wrap the pair instead.
      list(res = res, fold_warnings = msgs)
    }
    results <- parallel::mclapply(
      seq_along(remapped_folds), caught_worker, mc.cores = cores
    )
    relayed <- unique(unlist(lapply(results, function(z)
      if (!inherits(z, "try-error")) z$fold_warnings), use.names = FALSE))
    for (m in relayed) warning(m, call. = FALSE)
    results <- lapply(results, function(z) if (inherits(z, "try-error")) z else z$res)
  } else {
    results <- lapply(seq_along(remapped_folds), fold_worker)
  }

  # One status per fold, in the folds' own order, before anything is filtered:
  # what each fold did is the diagnosis a caller needs when "3 of 5 folds
  # produced predictions" is all the console kept.  mclapply() hands back a
  # try-error OBJECT (not NULL) when a child errors or is killed; a fold that
  # threw comes back as list(error = <message>); one skipped before fitting
  # or after predicting as list(skip = <reason>); a successful one carries
  # pred_row and fold_stat.  `$` (not `[[`) throughout, because a successful
  # fold's list has no "error" element and `[[` would abort.
  labels <- vapply(remapped_folds, function(f) as.integer(f$fold_id %||% NA), integer(1))
  labels[is.na(labels)] <- seq_along(remapped_folds)[is.na(labels)]
  status <- character(n_folds); msg <- character(n_folds)
  for (i in seq_len(n_folds)) {
    z <- results[[i]]
    if (inherits(z, "try-error")) {
      status[i] <- "worker_error"; msg[i] <- .try_error_message(z)
    } else if (is.null(z)) {
      status[i] <- "skipped"; msg[i] <- "no result returned"
    } else if (!is.null(z$error)) {
      status[i] <- "error"; msg[i] <- z$error
    } else if (!is.null(z$skip)) {
      status[i] <- "skipped"; msg[i] <- z$skip
    } else {
      status[i] <- "ok"; msg[i] <- ""
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
#' subsample when there are more; exact duplicate locations are dropped first,
#' because a spatial correlation structure cannot take a zero distance between
#' distinct observations) and the trend coefficients are then applied to every
#' row.  Several starting ranges are tried; the first fit that converges is
#' kept.
#'
#' @param mf A model frame (trend formula already applied, incomplete rows
#'   excluded) with the projected coordinates in columns \code{.sac_x} and
#'   \code{.sac_y}.
#' @param fml The trend formula.
#' @param extent A length scale of the layer, used to pick starting ranges.
#' @return \code{NULL} when no fit converged, else a list with \code{beta}
#'   (named trend coefficients), \code{range} (the exponential range
#'   parameter), \code{nugget_prop}, \code{sigma2}, \code{n_used} and
#'   \code{subsampled}.
#' @keywords internal
#' @noRd
.reml_trend <- function(mf, fml, extent, max_n = 400L, seed = 123L) {
  if (!requireNamespace("nlme", quietly = TRUE)) return(NULL)
  d <- mf
  # Duplicate locations: keep the first of each.
  dup <- duplicated(d[, c(".sac_x", ".sac_y")])
  if (any(dup)) d <- d[!dup, , drop = FALSE]
  n <- nrow(d)
  subsampled <- FALSE
  if (n > max_n) {
    cleanup <- .with_seed(seed)
    on.exit(cleanup(), add = TRUE)
    d <- d[sample.int(n, max_n), , drop = FALSE]
    subsampled <- TRUE
  }
  if (nrow(d) < 30L) return(NULL)
  starts <- unique(pmax(extent * c(1 / 10, 1 / 30, 1 / 3), sqrt(.Machine$double.eps)))
  for (r0 in starts) {
    fit <- tryCatch(
      withCallingHandlers(
        nlme::gls(fml, data = d,
                  correlation = nlme::corExp(value = c(r0, 0.1),
                                             form = ~ .sac_x + .sac_y,
                                             nugget = TRUE),
                  method = "REML"),
        warning = function(w) {
          .log_info("estimate_sac_range(): nlme::gls() warned during REML detrending: %s",
                    conditionMessage(w))
          invokeRestart("muffleWarning")
        }),
      error = function(e) NULL)
    if (is.null(fit)) next
    cs <- try(stats::coef(fit$modelStruct$corStruct, unconstrained = FALSE),
              silent = TRUE)
    if (inherits(cs, "try-error") || !all(is.finite(cs))) next
    rng <- as.numeric(cs[["range"]]); nug <- as.numeric(cs[["nugget"]])
    s2  <- as.numeric(fit$sigma)^2
    if (!is.finite(rng) || rng <= 0 || !is.finite(s2) || s2 <= 0) next
    return(list(beta = stats::coef(fit), range = rng, nugget_prop = nug,
                sigma2 = s2, n_used = nrow(d), subsampled = subsampled))
  }
  NULL
}


#' Does the empirical variogram decrease with distance over its shorter lags?
#'
#' A semivariance that falls as distance grows leaves the range unidentified:
#' a model with a monotone rise to a sill is being fitted to a curve that has
#' none.  The net change in \code{gamma} across the bins in the shorter half
#' of the lags, weighted by the pairs supporting each step, is compared with
#' the mean semivariance over those bins; a net fall of more than \code{tol}
#' of that mean is a decrease.  Measured on 60 draws each (n = 250, tol =
#' 0.15): 0 of an exponential field, 2 percent of white noise, 98 percent of a
#' field with a periodic (hole-effect) component, 100 percent of a layer whose
#' variance differs between a dense cluster and the rest.  The tolerance is
#' fixed while the noise in the short-lag bins grows as the sample shrinks, so
#' small samples trip it on ordinary fields: an exponential field with
#' effective range 300 and nugget 0.2 was refused in 7--9 of 60 draws at
#' n = 30, 3--6 at n = 50 and 0--1 at n = 100, and with range 150 in 15--16
#' of 60 at n = 30.  The refusal is the conservative outcome (geometric
#' blocks, with a warning), so it is left as it is.  An unremoved trend
#' is \emph{not} what produces this shape.  A trend makes the variogram rise
#' without reaching a sill, which the over-cutoff rejection catches, so the
#' message aimed at this case must not say "trend".
#'
#' @param vg An empirical variogram from \code{gstat::variogram()}.
#' @return \code{TRUE} for a net decrease, \code{FALSE} otherwise (including
#'   when fewer than \code{min_bins} bins fall in the shorter half).
#' @keywords internal
#' @noRd
.variogram_decreasing <- function(vg, frac = 0.5, min_bins = 4L, tol = 0.15) {
  if (!is.data.frame(vg) || !all(c("dist", "gamma", "np") %in% names(vg)))
    return(FALSE)
  vg <- vg[is.finite(vg$gamma) & is.finite(vg$dist) & is.finite(vg$np) & vg$np > 0, ,
           drop = FALSE]
  if (nrow(vg) < min_bins) return(FALSE)
  vg <- vg[order(vg$dist), , drop = FALSE]
  short <- vg[vg$dist <= frac * max(vg$dist), , drop = FALSE]
  if (nrow(short) < min_bins) return(FALSE)
  w   <- pmin(utils::head(short$np, -1L), utils::tail(short$np, -1L))
  if (!any(w > 0)) return(FALSE)
  net <- sum(w * diff(short$gamma)) / sum(w) * (nrow(short) - 1L)
  scale <- stats::weighted.mean(short$gamma, short$np)
  is.finite(net) && is.finite(scale) && scale > 0 && net < -tol * scale
}


#' The nugget of an estimated autocorrelation range
#'
#' The nugget variance of the variogram model behind a
#' \code{\link{estimate_sac_range}()} result: the semivariance at zero
#' separation, i.e. measurement error plus variation at scales shorter than
#' the closest pair.  It is carried as the \code{nugget} attribute of every
#' classed result, identified or rejected, because it is the number a
#' resolution criterion for a tessellation needs (the short-lag variance that
#' no cell can average away).
#'
#' @param x A \code{sac_range} object, or anything else.
#' @return A single number: the nugget in the units of the response's
#'   variance; \code{NA_real_} when \code{x} carries no fitted model (a bare
#'   \code{NA} from a run that could not fit anything, a rejected result whose
#'   fits were all singular, or an object that is not a \code{sac_range}).
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
#' Fits exponential (or spherical) variogram models and returns the
#' \emph{effective range}: for the exponential model, three times the fitted
#' range parameter, which is where the semivariance reaches ~95 % of the
#' sill; for the spherical model (fitted only when the exponential fit is
#' singular or does not converge) the fitted range itself, which is where the
#' spherical semivariance reaches its sill exactly.  Both are the distance
#' beyond which two observations are (near) uncorrelated, which is what a
#' block or a buffer has to exceed.
#'
#' The exponential model is kept whenever it converges, without comparing it
#' with the spherical fit, and on fields smoother than exponential that makes
#' the range long.  Measured on simulated fields (n = 300 on a 1000 m square,
#' 30 draws each): about 1.8--2.1 times the practical range of a Gaussian
#' covariance, and 1.3--1.4 times the range of a spherical one, while an
#' exponential field came back at 0.97 of its effective range.  The error is
#' on the safe side (blocks too large, cross-validation pessimistic), and it
#' is kept on purpose: choosing the family by the smaller weighted sum of
#' squares corrects the spherical case but sends exponential fields low, to
#' about 0.82 of the truth, which is the direction that leaks.
#'
#' The estimate is the \strong{omnidirectional} (all-pairs) fit.  Directional
#' variograms are fitted as well, at 0° (N–S), 45°, 90° (E–W) and 135°
#' azimuths with a ±22.5° tolerance.  Those four windows tile all 180 distinct
#' azimuths exactly once, and their ranges are returned in the
#' \code{directional} attribute, with their largest-over-smallest ratio in
#' \code{anisotropy}.  They are a diagnostic, not the answer, for two reasons.
#' Each direction sees about a quarter of the point pairs, and the maximum of
#' four quarter-sample fits is biased upward: on simulated \emph{isotropic}
#' fields it came in about 40% above the truth, and no hurdle placed in
#' front of it (all four directions fitted, ratio above 1.5, maximum above
#' 1.5× the all-pairs fit) kept it out.  One isotropic field rotated in 10°
#' steps "established" anisotropy in 14 of 18 orientations.  And the windows
#' are fixed to the coordinate axes, so any answer built from them changes
#' when the layer is rotated, which a property of the field must not do.  The
#' all-pairs fit is the best-powered estimate available and is invariant to
#' rotation.
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
#' and the ratio, and \code{attr(range, "anisotropy")} holds it.  Only when
#' the omnidirectional fit is singular or did not converge is the directional
#' maximum returned in its place, and \code{anisotropy_used} is \code{TRUE}
#' in that case alone.  An omnidirectional fit that converged to a range past
#' the fitted lags is refused (see the Value section) whatever the directions
#' found: the directions that reached a sill are the shorter ones, so their
#' maximum is a lower bound, not an estimate.
#'
#' A direction whose fit fails, does not converge, or reports a range beyond
#' the longest fitted lag is excluded and recorded as \code{NA} in the
#' \code{directional} attribute.
#'
#' Every variogram model is fitted \strong{with a nugget}.  A nugget-free model
#' forces the curve through the origin, and on any real measurement (which has
#' one) \pkg{gstat}'s default N/h² weights buy that constraint by collapsing
#' the range: with a 50% nugget the fitted range came back at about 0.45 of the
#' truth, so \code{make_folds(auto_range = TRUE)} built blocks less than half
#' the correlation length it reported.
#'
#' The lags are binned the way \pkg{gstat} bins them by default, 15 bins out
#' to the cutoff, each \code{cutoff * max_dist / 15} wide (about 47 m on a
#' 1000 m square at the defaults).  A range spanning only one or two bins is
#' resolved coarsely and comes out long: exponential fields with an effective
#' range of 60 m (n = 300 on a 1000 m square, 30 draws) returned a median of
#' 89--102 m, where the same fields binned over a 200 m cutoff gave 65--68,
#' and at a range of 300 m there was no bias.  When the estimate is within a
#' few bin widths of zero, run it again with a smaller \code{cutoff}.
#'
#' Nothing tests whether the layer has spatial structure at all.  On white
#' noise (n = 300 on a 1000 m square, 30 draws) the estimate was a finite,
#' spurious range (57--533 m) in 8 draws and a refusal in the rest, mostly
#' as past the fitted lags or not converged, and only once as no model
#' fitted; with \code{detrend = "reml"} it was finite in 10 of 30.  A
#' spurious range errs towards larger blocks, so the harm is mostly lost
#' training data, but a caller who needs to know whether there is any
#' structure should look at the variogram (\code{plot()} on the result)
#' rather than at whether the answer is \code{NA}.
#'
#' A log warning is emitted when the directional maximum is used; where the
#' all-pairs estimate is available it names both the ratio and that estimate.  A
#' log note is emitted instead when the directional ranges vary but the spread
#' is consistent with sampling noise.
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
#' @param response_var Character(1) name of the response column.  For a
#'   count or other response whose variance tracks its mean, see the section
#'   on non-Gaussian responses: the range is still estimated, but it is a
#'   less reliable number than for a Gaussian response.
#' @param predictor_vars Optional character vector.  When supplied, the
#'   trend on these predictors is removed first and the variogram describes
#'   the residual autocorrelation, the part a spatial model has to handle
#'   once the covariates have done their work.  How the trend is removed is
#'   set by \code{detrend}, and it matters: see "Detrending and the
#'   residual-variogram bias".
#' @param n_max Maximum number of points to subsample before fitting.
#'   Variogram estimation is O(n²) so this keeps runtime bounded.
#' @param cutoff Fraction of the maximum inter-point distance to use as the
#'   variogram lag cutoff.  That distance is the farthest pair, found on the
#'   convex hull, and not the bounding-box diagonal, which depends on how the
#'   axes are oriented.  Default 0.5.
#' @param range_frac Positive numeric.  A fitted range exceeding
#'   \code{range_frac * cutoff * max_dist} (that is, beyond the longest lag
#'   the empirical variogram was actually fitted over) is treated as
#'   unidentified and \code{NA_real_} is returned.
#'   \code{gstat::fit.variogram()} yields a finite number even when the
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
#'   size (see its \code{auto_range} argument for what it does then).  On a
#'   1000 m square with an exponential field of effective range 570 (n = 300),
#'   9 of 30 draws were accepted in that band.  Lowering \code{range_frac} to fit the
#'   grid would turn those estimates into \code{NA} and the blocks into
#'   geometric ones smaller than the range.
#' @param seed RNG seed for the \code{n_max} subsample, restored afterwards so
#'   the caller's random stream is untouched.  Default \code{123L}: the
#'   subsample is an internal approximation and no part of the answer, and
#'   leaving it unseeded made the returned range differ between runs on
#'   identical input (19531, 19589, 19605 on three calls) and silently advanced
#'   the caller's RNG.  Pass \code{NULL} for the old unseeded behaviour, or a
#'   different number to check how sensitive the estimate is to the subsample.
#'   Ignored when \code{nrow(points_sf) <= n_max}, where nothing is sampled.
#'   The \code{reml_max_n} subsample uses the same seed.
#' @param detrend How the trend on \code{predictor_vars} is removed;
#'   ignored when there are none.  \code{"ols"} (default) fits it by ordinary
#'   least squares and fits the variogram to the residuals.  This is the
#'   long-standing behaviour, which underestimates the range (see the section
#'   below).  \code{"reml"} fits the trend and an exponential-plus-nugget
#'   covariance together by residual maximum likelihood with
#'   \code{nlme::gls()}, returns the REML range, and attaches the empirical
#'   variogram of the REML residuals for inspection.  It needs \pkg{nlme},
#'   costs \eqn{O(n^3)} (about 2 s at 300 points, 12 s at 500, 45 s at 800),
#'   and so runs on at most \code{reml_max_n} points; when it does not
#'   converge the \code{"ols"} path runs instead with an R warning saying so.
#' @param reml_max_n Positive integer, at least 30.  With \code{detrend =
#'   "reml"}, the trend and covariance are fitted on a seeded random subsample
#'   of this many points when the layer has more (exact duplicate locations
#'   are dropped first); the fitted trend is then removed from every point.
#'   Default 400.  Raise it for a better-determined fit at the cost above.
#' @param keep_directional_fits Logical.  Attach each direction's empirical
#'   variogram and fitted model as \code{directional_fits}?  Defaults to
#'   \code{FALSE}: the four variograms are most of the object's size (42.1 KB
#'   of 59.3 KB at \eqn{n = 400}, and the difference between a 52.4 KB and a
#'   94.6 KB \code{make_folds(auto_range = TRUE)} result), while the numbers
#'   read from them (\code{directional}, \code{directional_fitted},
#'   \code{directional_status}, \code{anisotropy}) are attached either
#'   way, and \code{plot()} draws the effective variogram from its own
#'   attribute.  Set \code{TRUE} to inspect the directional curves.
#' @section Detrending and the residual-variogram bias:
#' Fitting a variogram to the residuals of a least-squares trend
#' underestimates both the sill and the range, because the trend fit absorbs
#' part of the long-wavelength spatial variation (Lark, Cullis and Welham
#' 2006).  Blocks sized from that range are then too small and a blocked
#' validation is less conservative than it claims.  How large the effect is
#' depends on how smooth the trend terms are in space, and on how many there
#' are.  Measured for this estimator on simulated exponential fields (n =
#' 300, true effective range 300, nugget 0.2, 40--60 draws), as the median
#' ratio of the estimate from the trend-removed data to the estimate from
#' the true field:
#' \itemize{
#'   \item a white-noise covariate: OLS 1.00, REML 0.99 (no bias to speak of);
#'   \item a spatially smooth covariate (a random field with range 300 or
#'     1000): OLS 0.97, REML 0.95--0.98;
#'   \item a linear trend in the coordinates: OLS 0.92, REML 1.04;
#'   \item a quadratic trend in the coordinates (five terms): OLS 0.75,
#'     REML 1.06.
#' }
#' So for ordinary covariates the OLS bias is a few percent, and for trend
#' surfaces in the coordinates it is large.  Iterating between a GLS trend
#' fit and a variogram refit (Neuman and Jacobson 1984) recovers only part of
#' it (0.80 in the quadratic case), because the variogram of GLS residuals
#' is biased too; REML does not fit a variogram to residuals at all, which is
#' why \code{detrend = "reml"} returns its own range estimate.  Its price is
#' a single family (exponential with nugget), a cubic cost in \code{n}, and
#' a somewhat wider sampling spread.  The default stays \code{"ols"} so that
#' existing scripts return what they did; a script that detrends on smooth
#' or coordinate-based terms should pass \code{detrend = "reml"}, or size its
#' blocks with a margin.
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
#' does not compute.  Until it does, treat the range from a count response as no
#' more than an order of magnitude, size blocks conservatively from it, and
#' prefer \code{\link{make_folds}(method = "nndm")}, which does not depend on a
#' fitted range at all.
#'
#' @return A single number, of class \code{sac_range} in the first two of the
#'   three shapes below and a bare \code{NA} in the third; all three behave
#'   as an ordinary number.  The shapes carry different attributes:
#'   \describe{
#'     \item{Success}{A positive effective range in projected coordinate units,
#'       with the fit attached as attributes \code{directional} (the 0°, 45°,
#'       90° and 135° ranges, named by azimuth; \code{NA} where that
#'       direction's fit was unusable), \code{anisotropy} (largest
#'       over smallest), \code{anisotropy_used} (logical: \code{TRUE} only when
#'       the all-pairs fit was singular or did not converge and the
#'       directional maximum stands in for it), \code{directional_status}
#'       (per azimuth, why a direction is
#'       \code{NA} in \code{directional}: \code{"ok"}, \code{"over_cutoff"}
#'       (its range ran past the largest lag fitted), \code{"not_converged"}
#'       or \code{"no_fit"}), \code{directional_fitted} (the range each
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
#'       \code{subsampled}, \code{nugget_prop} and \code{sigma2} from the
#'       REML fit; \code{NULL} otherwise),
#'       \code{crs} (the projected CRS the variogram was
#'       fitted in: the unit of the range), \code{max_dist},
#'       \code{cutoff_dist}, \code{variogram} (the empirical variogram),
#'       \code{variogram_model} (the fitted \code{gstat} model, or with
#'       \code{detrend = "reml"} a \code{gstat} model built from the REML
#'       parameters) and \code{nugget} (that model's nugget variance; see
#'       \code{\link{sac_nugget}}), so the fit can be inspected and need not
#'       be taken on trust.}
#'     \item{Rejected range}{\code{NA_real_} when a range was fitted but is
#'       not identified: it exceeds \code{range_frac * cutoff * max_dist} (see
#'       \code{range_frac}), which applies to the all-pairs fit even when some
#'       directions reached a sill; or it is shorter than the shortest lag the
#'       empirical variogram resolves (the mean separation in its first
#'       bin), below which a structure cannot be told from a nugget; or the
#'       model did not converge; or the empirical
#'       variogram \emph{decreases} with distance over its shorter lags (a
#'       net fall of more than 15 percent of the mean semivariance there,
#'       weighted by pairs), which is the shape of a periodic, hole-effect
#'       structure or of a variance that differs between a dense cluster and
#'       the rest of the layer.  Sampling noise in the short-lag bins of a
#'       small sample can make that fall too: on exponential fields it
#'       refused 7--9 of 60 draws at n = 30, 3--6 at n = 50 and 0--1 at
#'       n = 100 (effective range 300 on a 1000 m square), and 15--16 of 60
#'       at n = 30 with a range of 150.  An unremoved trend makes the
#'       variogram rise instead; when it rises past the fitted lags the first
#'       test catches it, but a milder trend only lengthens the fitted range
#'       and passes, which is what \code{predictor_vars} is for.  Last, the
#'       fitted range can be non-positive.  It is classed \code{sac_range} as
#'       well, so it prints as \code{NA} without dumping its
#'       attributes, and it carries \code{max_dist}, \code{cutoff_dist},
#'       \code{variogram}, \code{variogram_model} and \code{nugget} (the
#'       evidence for the rejection), plus \code{rejected_range} (the value
#'       that was refused), \code{rejected_reason} (one of
#'       \code{"fitted range exceeds the largest lag fitted"},
#'       \code{"fitted range is below the shortest lag fitted"},
#'       \code{"variogram model did not converge"},
#'       \code{"empirical variogram decreases with distance"},
#'       \code{"fitted range is non-positive or non-finite"},
#'       \code{"no variogram model could be fitted (singular fits)"}), \code{crs}
#'       (so the units the rejected number was in stay recoverable, which is
#'       what \code{plot()} labels its axis from) and
#'       \code{detrend_method}.  It carries \code{directional},
#'       \code{anisotropy}, \code{anisotropy_used}, \code{directional_status},
#'       \code{directional_fitted} and, with
#'       \code{keep_directional_fits = TRUE}, \code{directional_fits} as well:
#'       the directional sweep runs whatever becomes of the all-pairs fit, and
#'       its per-azimuth outcome is what says whether any direction reached
#'       a sill the pooled variogram did not, or whether every direction ran
#'       past the fitted lags alike.  The same shape, with
#'       \code{rejected_range = NA}, \code{variogram_model = NULL} and
#'       \code{nugget = NA}, is returned when no variogram model could be
#'       fitted at all (both the exponential and the spherical fit singular,
#'       which a flat, nugget-only variogram can produce, though on white
#'       noise it was the outcome in only 1 of 30 draws: see above);
#'       \code{rejected_reason} says so and the empirical variogram is still
#'       attached.}
#'     \item{No fit}{A bare, attribute-less \code{NA_real_} when estimation
#'       could not be attempted at all: \pkg{gstat} missing, fewer than 30
#'       finite values, a variable with no variance, or a degenerate extent.
#'       Without \pkg{gstat} nothing is fitted, so none of the attributes
#'       above exist either.}
#'   }
#'   Attributes and the class do not affect \code{is.na()} or
#'   \code{is.finite()}, so every downstream guard treats all three the same
#'   way it always did.
#' @references
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
#'   print(attr(r, "variogram_model"))     # the fitted gstat model behind it
#'
#'   # A field whose range the data cannot pin down: the variogram never
#'   # reaches a sill within the lags fitted, so the answer is NA with the
#'   # refused value attached rather than a long range asserted.
#'   xy$trend <- sin(xy$x / 400) + rnorm(n, sd = 0.2)
#'   r2 <- estimate_sac_range(st_as_sf(xy, coords = c("x", "y"), crs = 32632),
#'                            response_var = "trend")
#'   print(r2)
#'   attr(r2, "rejected_range")
#' }
#' @export
estimate_sac_range <- function(points_sf, response_var,
                               predictor_vars = NULL,
                               n_max = 5000L, cutoff = 0.5,
                               range_frac = 1.0, seed = 123L,
                               detrend = c("ols", "reml"),
                               reml_max_n = 400L,
                               keep_directional_fits = FALSE) {
  detrend <- match.arg(detrend)
  if (!is.numeric(reml_max_n) || length(reml_max_n) != 1L ||
      !is.finite(reml_max_n) || reml_max_n < 30)
    stop("estimate_sac_range(): `reml_max_n` must be a single number of at ",
         "least 30.", call. = FALSE)
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
  # into NA under a "no variogram model could be fitted" message -- blaming
  # the fit rather than the row, and disagreeing with make_folds(), which drops such
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
  detrended <- FALSE
  detrend_method <- NA_character_
  reml_fit <- NULL
  if (!is.null(predictor_vars) && length(predictor_vars) > 0L) {
    df <- pts_df
    # A predictor that is not a column is an error, as it is everywhere else
    # in the package: intersect() used to drop it silently, and with every
    # name unknown the RAW response was modelled -- so make_folds(auto_range =
    # TRUE) sized blocks from range 88.5 instead of the residual range 362.2
    # with nothing said.
    missing_preds <- setdiff(predictor_vars, names(df))
    if (length(missing_preds))
      stop("estimate_sac_range(): predictor_vars ",
           paste(sQuote(missing_preds), collapse = ", "),
           " not found in the data.", call. = FALSE)
    fml <- stats::reformulate(predictor_vars, response_var)
    # REML: trend and covariance fitted together, so the trend is a GLS fit
    # under the fitted correlation and the range is the REML estimate, not a
    # variogram of residuals.  Measured on simulated fields (n = 300, true
    # effective range 300, 40-60 draws): the OLS residual variogram returned
    # a median 0.97 of the oracle range for a spatially smooth covariate, 0.92
    # for a linear trend in the coordinates and 0.75 for a quadratic one; REML
    # 0.95-1.06 throughout.  On failure the OLS path below runs instead.
    if (identical(detrend, "reml")) {
      xy_tr <- sf::st_coordinates(pts)[, 1:2, drop = FALSE]
      df$.sac_x <- xy_tr[, 1]; df$.sac_y <- xy_tr[, 2]
      # Complete rows only, chosen here rather than by na.action so the same
      # row set can carry the coordinates alongside the model frame and the
      # residuals can be put back at their positions afterwards.
      keep <- stats::complete.cases(df[, c(response_var, predictor_vars), drop = FALSE]) &
        is.finite(df$.sac_x) & is.finite(df$.sac_y)
      mf <- df[keep, c(response_var, predictor_vars, ".sac_x", ".sac_y"), drop = FALSE]
      if (sum(keep) >= 30L) {
        extent <- max(diff(range(xy_tr[keep, 1])), diff(range(xy_tr[keep, 2])), 1)
        reml_fit <- .reml_trend(mf, fml, extent = extent, max_n = reml_max_n,
                                seed = seed)
      }
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
        if (inherits(X, "try-error") || anyNA(b) || nrow(X) != sum(keep)) {
          .warn_and_log("estimate_sac_range(): the REML trend coefficients do not match the design; falling back to OLS detrending.")
          reml_fit <- NULL
        } else {
          r_all <- rep(NA_real_, nrow(df))
          r_all[keep] <- as.numeric(mf[[response_var]] - X %*% b)
          y <- r_all
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
    lm_fit <- if (is.null(reml_fit))
      try(stats::lm(fml, data = df, na.action = stats::na.exclude), silent = TRUE)
    else NULL
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
      if (length(resid) != nrow(pts)) {
        .warn_and_log("estimate_sac_range(): OLS residual length (%d) does not match data rows (%d); the variogram is fitted to the RAW response instead.",
                      length(resid), nrow(pts))
      } else {
        y <- resid
        detrended <- TRUE
        detrend_method <- "ols"
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

  # Empirical variogram.  The lag cutoff is a fraction of the maximum
  # inter-point distance, as documented -- NOT of the bounding-box diagonal.
  # The diagonal is a property of the axes, not of the points: it grows by up
  # to sqrt(2) when the same layer is rotated 45 degrees, so the variogram was
  # binned differently and the fitted range moved (257 vs 394 on one field;
  # pinning the cutoff restored 257 exactly).  The farthest pair lies on the
  # convex hull, whose vertex count is tiny, so this is cheap at any n.
  max_dist <- try({
    hull <- sf::st_convex_hull(sf::st_union(sf::st_geometry(pts)))
    hv   <- unique(sf::st_coordinates(hull)[, 1:2, drop = FALSE])
    if (nrow(hv) < 2L) 0 else max(stats::dist(hv))
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
    fit_from <- function(range0) {
      converged <- TRUE
      m <- withCallingHandlers(
        try(gstat::fit.variogram(
              vg, gstat::vgm(psill = NA, model = model_type, range = range0,
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
      attr(m, "converged") <- converged
      m
    }
    # SEVERAL starting ranges, not gstat's one.  With `range = NA` gstat
    # starts the optimiser at a third of the longest lag; for a field whose
    # range is a small fraction of the extent that start is ten times too
    # long, the nugget and the partial sill trade off along the way, and
    # whether the Gauss-Newton iteration lands or collapses to a singular
    # model depends on floating-point details -- the same 250-point field fitted
    # on one machine and came back singular on another.  Shorter and longer
    # starts are tried as well, and the winner is the converged, non-singular
    # fit with the smallest weighted sum of squares (gstat's own criterion,
    # attr "SSErr"), which is a property of the data rather than of the path
    # the optimiser took.  A non-converged fit is kept only when nothing
    # converged: its range must never size a block -- that is what the
    # `converged` flag is for -- but the model and its empirical variogram are
    # still the most useful thing a user can look at, and a sill-less
    # variogram is precisely the case worth looking at.
    dmax <- max(vg$dist, na.rm = TRUE)
    starts <- if (is.finite(dmax) && dmax > 0)
      c(NA_real_, dmax / 10, dmax / 30, dmax / 2) else NA_real_
    fits <- Filter(Negate(is.null), lapply(starts, fit_from))
    if (!length(fits)) return(NULL)
    conv <- vapply(fits, function(m) !identical(attr(m, "converged"), FALSE),
                   logical(1))
    if (any(conv)) fits <- fits[conv]
    sse <- vapply(fits, function(m) {
      v <- attr(m, "SSErr")
      if (is.null(v) || !is.finite(v)) Inf else as.numeric(v)
    }, numeric(1))
    fits[[which.min(sse)]]
  }

  .fit_vgm_range <- function(vg) {
    if (inherits(vg, "try-error") || !inherits(vg, "data.frame") || NROW(vg) < 3L) return(NA_real_)
    # Exponential first, spherical as the fallback -- but a CONVERGED fit of
    # either beats a non-converged fit of the preferred one.  On a field with
    # no nugget the exponential fit runs into the nugget's zero bound and
    # gstat reports non-convergence from every start, while the spherical fit
    # converges cleanly; taking the exponential's non-converged range would
    # then throw away a usable estimate (and the caller must refuse a
    # non-converged range, so the omnidirectional fit would count as failed).
    fits <- list(Exp = .fit_one_vgm(vg, "Exp"))
    conv <- function(m) !is.null(m) && !identical(attr(m, "converged"), FALSE)
    if (!conv(fits$Exp)) fits$Sph <- .fit_one_vgm(vg, "Sph")
    vgm_model <- if (conv(fits$Exp)) fits$Exp
      else if (conv(fits$Sph)) fits$Sph
      else if (!is.null(fits$Exp)) fits$Exp
      else fits$Sph
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

  # The per-azimuth outcome, kept rather than collapsed: `usable` below holds
  # NA for every direction that failed, but a fit that never converged, one
  # whose range ran past the fitted lags and one that could not be fitted at
  # all are three different findings, and the refused range itself -- the
  # most informative number when a field is anisotropic and the sweep could
  # not use it -- was unrecoverable.  All three ride on every classed return,
  # the rejected ones included: the four fits ran regardless, and the
  # rejected path is precisely where a user needs to know whether the field
  # is anisotropic.  The empirical variogram and fitted model of each
  # direction travel too, so the sweep can be drawn.
  dir_status <- ifelse(!is.finite(dir_ranges), "no_fit",
                ifelse(!dir_conv, "not_converged",
                ifelse(dir_ranges > max_supported, "over_cutoff", "ok")))
  dir_status <- stats::setNames(dir_status, as.character(dir_az))
  dir_fitted <- stats::setNames(dir_ranges, as.character(dir_az))
  dir_detail <- stats::setNames(lapply(dir_fits, function(f)
    list(variogram = if (inherits(f$vg, "data.frame")) f$vg else NULL,
         model     = attr(f$fit, "vgm_model"))), as.character(dir_az))
  # Four empirical variograms, one per azimuth, are most of what this object
  # weighs: at n = 400 they are 42.1 KB of a 59.3 KB estimate, and
  # make_folds(auto_range = TRUE) parks one in `params`, which took that
  # folds object from 52.4 KB to 94.6 KB -- for something nothing in the
  # package reads.  plot() draws the effective variogram from the
  # `variogram` attribute, and the per-direction ranges and outcomes are in
  # `directional`, `directional_fitted` and `directional_status` either way.
  # Kept on request, for looking at the directional curves themselves.
  dir_detail_out <- if (isTRUE(keep_directional_fits)) dir_detail else NULL

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
  # A converged all-pairs fit whose range runs past the fitted lags is not an
  # unusable fit but a finding: the pooled variogram, which sees every pair,
  # reached no sill.  The directional maximum must not stand in for it.  The
  # directions that did reach a sill are by construction the SHORTER ones (a
  # trend's cross-slope directions, an anisotropic field's minor axes), so
  # their maximum is a lower bound on the range, not an estimate of it, and
  # whether two of them happened to fit flipped the answer between NA and a
  # finite range from one draw to the next: with an east-west trend on an
  # exponential field of range 150, 12 of 30 draws returned 131-596 this way
  # with anisotropy_used = TRUE, and 16 others NA.  It goes to the rejection
  # below, as the documentation always said a trend would.
  iso_over <- is.finite(iso_fit_always) &&
              !identical(attr(iso_fit_always, "converged"), FALSE) &&
              as.numeric(iso_fit_always) > max_supported

  if (!is.null(reml_fit)) {
    # --- REML detrending: the range is the REML estimate ---------------------
    # Not a variogram fitted to the REML residuals: a residual variogram is
    # biased toward a shorter range whatever fitted the trend (measured 0.80
    # of the oracle for GLS residuals under a quadratic trend), and the REML
    # parameters are the estimate that does not carry that bias.  The
    # empirical variogram of the residuals is attached for inspection -- the
    # model line drawn through it is the REML model, and a curve sitting
    # below that line at long lags is the bias made visible -- and the
    # directional sweep on the residuals stays the diagnostic it always is.
    effective_range <- 3 * reml_fit$range
    vgm_used <- gstat::vgm(psill  = reml_fit$sigma2 * (1 - reml_fit$nugget_prop),
                           model  = "Exp", range = reml_fit$range,
                           nugget = reml_fit$sigma2 * reml_fit$nugget_prop)
    vg_used  <- if (inherits(vg_iso_always, "data.frame")) vg_iso_always else NULL
    if (dir_success)
      anisotropy <- max(usable, na.rm = TRUE) / min(usable, na.rm = TRUE)
  } else if (dir_success && !iso_over) {
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
    # The directional maximum is NEVER preferred over a usable all-pairs fit.
    # Three hurdles were tried (all four directions fitted, ratio > 1.5,
    # maximum > 1.5 x the all-pairs estimate) and still let the noise through:
    # on ONE isotropic exponential field (true range 150, n = 300) they
    # declared anisotropy in 14 of 18 axis orientations and returned ranges
    # from 225 to 529 -- a 2.35x spread produced by nothing but the direction
    # the axes happened to point.  Four windows fixed in CRS azimuth cannot
    # give a rotation-invariant answer, and the maximum of four quarter-sample
    # fits is biased upward whatever hurdle is put in front of it.  So the
    # all-pairs range is the estimate; the directional ranges are reported as
    # a diagnostic, and a caller who KNOWS the field is anisotropic can size
    # blocks from the longest directional range explicitly.  Not from
    # max(attr(x, "directional")): on exactly such a field the major axis is
    # the direction most likely to run past the fitted lags, which makes it NA
    # there and the maximum NA (or, with na.rm = TRUE, the second-longest).
    if (iso_ok) {
      effective_range <- as.numeric(iso_fit_always)
      vg_used  <- vg_iso_always
      vgm_used <- attr(iso_fit_always, "vgm_model")
      if (is.finite(anisotropy) && anisotropy > 1.5)
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
          as.numeric(iso_fit_always))
    } else {
      # The isotropic fit is singular or did not converge (one that converged
      # past the fitted lags goes to the refusal instead, see `iso_over`); the
      # directional sweep is all there is.
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
    # Also reached when the all-pairs fit converged past the fitted lags
    # (`iso_over`), whatever the directions did, so the refusal below sees it.
    # The same all-pairs variogram and fit as above; it was recomputed here,
    # which cost a second fit and logged its failure twice.
    vg_iso    <- vg_iso_always
    iso_range <- iso_fit_always
    if (dir_success)
      anisotropy <- max(usable, na.rm = TRUE) / min(usable, na.rm = TRUE)
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
      # Neither directional nor isotropic succeeded.  The VALUE is NA, but the
      # empirical variogram is still the thing to look at: both fits being
      # singular is one thing a flat, nugget-only variogram produces --
      # residuals with no spatial structure at the lags resolved -- and
      # returning a bare NA left plot(type = "variogram") unable to draw
      # exactly that picture ("could not be fitted; there may be too few
      # finite residuals").  Only one: on 30 draws of white noise this branch
      # was reached once, and eight came back with a finite (spurious) range.
      .log_warn(paste0("estimate_sac_range(): no variogram model could be fitted ",
                       "(the exponential and spherical fits are both singular, ",
                       "which a flat, nugget-only variogram can produce); ",
                       "returning NA. The empirical variogram is attached for ",
                       "inspection: call plot() on the returned value."))
      return(structure(
        NA_real_,
        class           = c("sac_range", "numeric"),
        directional     = stats::setNames(usable, as.character(dir_az)),
        anisotropy      = anisotropy,
        anisotropy_used = FALSE,
        directional_status = dir_status,
        directional_fitted = dir_fitted,
        directional_fits   = dir_detail_out,
        detrended       = isTRUE(detrended),
        detrend_method  = detrend_method,
        max_dist        = as.numeric(max_dist),
        cutoff_dist     = as.numeric(cutoff_dist),
        crs             = sf::st_crs(pts),
        variogram       = if (inherits(vg_iso, "data.frame")) vg_iso else NULL,
        variogram_model = NULL,
        nugget          = NA_real_,
        rejected_range  = NA_real_,
        rejected_reason = "no variogram model could be fitted (singular fits)"
      ))
    }
  }

  # The nugget of whatever model stands behind the answer, identified or not:
  # the semivariance at zero separation, which a resolution criterion needs
  # and which used to be reachable only by reading gstat's row layout.
  nugget_val <- .vgm_nugget_of(vgm_used)

  if (!is.finite(effective_range) || effective_range <= 0) {
    # Classed like every other refusal, so the evidence travels with the NA;
    # this used to be the one path that returned a bare, attribute-less NA
    # from inside a completed fit.
    .log_warn("estimate_sac_range(): estimated range is non-positive or non-finite; returning NA.")
    return(structure(
      NA_real_,
      class           = c("sac_range", "numeric"),
      directional     = stats::setNames(usable, as.character(dir_az)),
      anisotropy      = anisotropy,
      anisotropy_used = isTRUE(aniso_used),
      directional_status = dir_status,
      directional_fitted = dir_fitted,
      directional_fits   = dir_detail_out,
      detrended       = isTRUE(detrended),
      detrend_method  = detrend_method,
      max_dist        = as.numeric(max_dist),
      cutoff_dist     = as.numeric(cutoff_dist),
      crs             = sf::st_crs(pts),
      variogram       = vg_used,
      variogram_model = vgm_used,
      nugget          = nugget_val,
      rejected_range  = as.numeric(effective_range),
      rejected_reason = "fitted range is non-positive or non-finite"
    ))
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
  # A semivariance that FALLS with distance over the shorter lags has no sill
  # to reach, so the range fitted through it is not identified either.  This
  # is not the trend signature -- a trend makes the curve rise without a sill,
  # which `over_cutoff` catches -- but a periodic (hole-effect) structure, or
  # a variance that differs between a dense cluster and the rest of the layer
  # (see .variogram_decreasing() for the measured rates).
  decreasing    <- .variogram_decreasing(vg_used)
  # The mirror of `over_cutoff` at the other end of the lags.  A range shorter
  # than the shortest lag the variogram resolves (the mean separation in its
  # first non-empty bin) describes a structure that has died out before the
  # closest pairs of points, which the data cannot tell from a nugget.  It is
  # what white noise gives the REML fit, which is not fitted to the binned
  # variogram at all: on 30 draws of iid noise (n = 300, 1000 m square) it
  # returned ranges of 0.18-23.6 m in 19, against a first lag of about 30 m,
  # and one of 0.27 m sized a 3642 x 3676 block grid.  The weighted fit to
  # the bins did not go below the first lag on those draws, so nothing it
  # identified before is refused.
  first_lag <- if (is.data.frame(vg_used) && all(c("dist", "np") %in% names(vg_used))) {
    d1 <- vg_used$dist[is.finite(vg_used$dist) & is.finite(vg_used$np) & vg_used$np > 0]
    if (length(d1)) min(d1) else NA_real_
  } else NA_real_
  under_lag <- is.finite(first_lag) && effective_range < first_lag
  # Non-convergence is refused on the same terms and for the same reason: the
  # number is not a fitted parameter.  gstat signals it with a warning and
  # returns anyway, which is why it needs its own test rather than riding on
  # the cutoff bound -- a non-converged range can land inside the bound and
  # would otherwise have sized a block.
  if (over_cutoff || !fit_converged || decreasing || under_lag) {
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
        # the short-lag bins: 7-9 of 60 ordinary exponential fields were
        # refused at n = 30, 0-1 at n = 100 (see .variogram_decreasing()).
        if (nrow(pts) < 100L)
          sprintf(paste0(" With %d points the short-lag bins are noisy enough ",
                         "for an ordinary field to show this shape as well."),
                  nrow(pts))
        else ""
      )
    } else if (over_cutoff) {
      .log_warn(
        paste0("estimate_sac_range(): fitted range (%.0f) exceeds the largest ",
               "lag the variogram was fitted over (%.4g = %s x cutoff %.0f); the ",
               "empirical variogram never reached a sill, so the range is ",
               "unidentified rather than long. Returning NA. Raise `cutoff` to ",
               "fit longer lags, supply `predictor_vars` to detrend, or set a ",
               "block size explicitly. (A variogram that is flat from the ",
               "first lag, with no spatial structure to find, can end here ",
               "too: plot() the returned value to tell the two apart.)"),
        effective_range, max_supported, format(range_frac), cutoff_dist
      )
    } else if (fit_converged) {
      .log_warn(
        paste0("estimate_sac_range(): the fitted range (%.3g) is shorter than ",
               "the shortest lag the variogram resolves (%.3g, the mean ",
               "separation in its first bin), so the structure it describes ",
               "dies out before the closest pairs of points and cannot be told ",
               "from a nugget. That is what a layer with no spatial structure ",
               "at the lags resolved gives. Returning NA."),
        effective_range, first_lag
      )
    } else {
      .log_warn(
        paste0("estimate_sac_range(): the variogram model did not converge ",
               "(gstat stopped at its iteration limit), so the range it ",
               "reports (%.0f) is where the optimiser halted rather than a ",
               "fitted parameter. Returning NA. Raise `cutoff` to fit longer ",
               "lags, supply `predictor_vars` to detrend, or set a block size ",
               "explicitly. The empirical variogram is attached for ",
               "inspection: call plot() on the returned value."),
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
      directional     = stats::setNames(usable, as.character(dir_az)),
      anisotropy      = anisotropy,
      anisotropy_used = isTRUE(aniso_used),
      directional_status = dir_status,
      directional_fitted = dir_fitted,
      directional_fits   = dir_detail_out,
      detrended       = isTRUE(detrended),
      detrend_method  = detrend_method,
      max_dist        = as.numeric(max_dist),
      cutoff_dist     = as.numeric(cutoff_dist),
      crs             = sf::st_crs(pts),
      variogram       = vg_used,
      variogram_model = vgm_used,
      nugget          = nugget_val,
      rejected_range  = as.numeric(effective_range),
      rejected_reason = if (decreasing)
        "empirical variogram decreases with distance"
      else if (over_cutoff)
        "fitted range exceeds the largest lag fitted"
      else if (fit_converged)
        "fitted range is below the shortest lag fitted"
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
    # Per azimuth: why a direction is NA in `directional`, the range its fit
    # reported whether or not it was usable, and the empirical variogram and
    # model behind it.
    directional_status = dir_status,
    directional_fitted = dir_fitted,
    directional_fits   = dir_detail_out,
    # Whether the variogram is of OLS residuals on `predictor_vars` (TRUE) or
    # of the raw response.  make_folds(auto_range = TRUE) and
    # summarize_by_cell(deff = "variogram") both need to know which.
    detrended       = isTRUE(detrended),
    # "ols" or "reml" when detrended, NA otherwise; and the REML fit's own
    # numbers when it was used, so the caller can see how much of the layer
    # the trend was estimated on.
    detrend_method  = detrend_method,
    reml            = if (is.null(reml_fit)) NULL else
      list(n_used = reml_fit$n_used, subsampled = isTRUE(reml_fit$subsampled),
           nugget_prop = reml_fit$nugget_prop, sigma2 = reml_fit$sigma2),
    max_dist        = as.numeric(max_dist),
    cutoff_dist     = as.numeric(cutoff_dist),
    # The CRS the variogram was fitted in.  Its range is a length in these
    # units; summarize_by_cell(deff = "variogram") transforms the points to it
    # before evaluating the correlation function at within-cell distances.
    crs             = sf::st_crs(pts),
    variogram       = vg_used,
    variogram_model = vgm_used,
    nugget          = nugget_val
  )
}


#' The nugget psill of a gstat variogram model, or NA
#' @keywords internal
#' @noRd
.vgm_nugget_of <- function(vm) {
  if (!is.data.frame(vm) || !all(c("model", "psill") %in% names(vm))) return(NA_real_)
  v <- sum(as.numeric(vm$psill[as.character(vm$model) == "Nug"]))
  if (length(v) != 1L || !is.finite(v)) NA_real_ else v
}


#' Print a spatial autocorrelation range
#'
#' Prints the effective range as a plain number, with the directional fit
#' summarised beneath it when one is available.  A direction whose fit was
#' unusable is labelled with why (\code{directional_status}) and the range
#' its fit reported (\code{directional_fitted}) when the object carries
#' them, and \code{unidentified} otherwise.  A last line names the unit and
#' the CRS the range is a length in (\code{attr(x, "crs")}, which for
#' lon/lat input is the projected CRS the estimate chose), and whether the
#' variogram is of the response or of its residuals on
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
  if (inherits(cr, "crs") && !is.na(cr)) {
    u <- tryCatch(cr$units_gdal, error = function(e) NULL)
    unit <- if (is.character(u) && length(u) == 1L && !is.na(u) && nzchar(u))
      switch(u, metre = "metres", kilometre = "kilometres", foot = "feet",
             "US survey foot" = "US survey feet", degree = "degrees", u)
    else "CRS units"
    dm <- attr(x, "detrend_method")
    what <- if (isTRUE(attr(x, "detrended")))
      sprintf("; variogram of the residuals on predictor_vars (%s)",
              if (is.character(dm) && length(dm) == 1L && !is.na(dm)) dm else "detrended")
    else if (isFALSE(attr(x, "detrended"))) "; variogram of the response itself"
    else ""
    cat(sprintf("  in %s of %s%s\n", unit, .fold_crs_label(cr), what))
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
#'   or \code{prediction_points} that carries a CRS.  They are reprojected when
#'   the coordinates look like lon/lat and otherwise stamped without
#'   reprojection, with a warning either way.
#' @param k Integer; number of folds.  Must be a single whole number >= 1.
#'   A fraction, \code{NA} or a vector is an error, because a non-integer used
#'   to truncate silently and leave the last rows in no test set at all.
#'   Not every method honours it.  \code{"buffered_loo"} and \code{"nndm"} are
#'   leave-one-out schemes and always return \code{k = n} regardless of what
#'   was asked for; \code{"block_kfold"} lowers it when the grid yields fewer
#'   than \code{k} non-empty blocks, and \code{"leave_location_out"} lowers it
#'   when there are fewer than \code{k} distinct groups.  Read the \code{k}
#'   element of the returned list, and do not assume the requested value.  A
#'   reduction is written to the package log and raises no R warning, so
#'   \code{tryCatch(warning = )} will not see it and \code{suppressWarnings()}
#'   will not hide it.
#' @param method One of \code{"random_kfold"}, \code{"block_kfold"},
#'   \code{"buffered_loo"}, \code{"leave_location_out"} or \code{"nndm"}.  See
#'   \strong{Details} for what each one does and when it is appropriate.
#' @param seed Optional integer RNG seed.
#' @param block_nx,block_ny Optional grid dimensions for block_kfold.
#'   Ignored when \code{block_size} or \code{auto_range} override them.
#' @param block_multiplier Numeric, default 3.  When neither \code{block_size}
#'   nor \code{block_nx}/\code{block_ny} is given, the automatic grid aims for
#'   \code{block_multiplier * k} blocks over the extent (aspect-preserving),
#'   so each fold holds out about \code{block_multiplier} blocks.  With 1,
#'   every fold is one contiguous region and the score depends heavily on
#'   which region each fold happened to get; with many, the blocks shrink
#'   towards single points and the scheme drifts back towards random k-fold.
#'   3 is a compromise between those two, not a published constant.  The
#'   block size that matters for leakage is the autocorrelation range, which
#'   is what \code{block_size} and \code{auto_range} control.
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
#'   range is estimated via \code{estimate_sac_range()} (which fits
#'   directional variograms to account for anisotropy) and used as the
#'   minimum \code{block_size}.  Requires \code{response_var}.  An explicit
#'   \code{block_size} takes precedence.  Default \code{FALSE}.  Sizing
#'   blocks from the autocorrelation range is the recommendation of Roberts
#'   et al. (2017) and what \pkg{blockCV} (Valavi et al. 2019) automates.
#'   \pkg{blockCV} takes the fitted variogram's range
#'   \emph{parameter} as the block size, whereas this uses the
#'   \emph{effective} range \code{estimate_sac_range()} returns (three
#'   times that parameter for an exponential fit), so its blocks are larger
#'   than \pkg{blockCV}'s from the same variogram.
#' @param range_frac Passed through to \code{estimate_sac_range()} when
#'   \code{auto_range = TRUE}.  A fitted range beyond the longest lag the
#'   empirical variogram was fitted over is rejected as unidentified, and block
#'   sizing falls back to geometry, so the grid does not collapse to a single
#'   block.
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
#'   sites it has already seen, which random k-fold reports as excellent
#'   performance.
#' @param prediction_points Optional \code{sf} layer of the locations you
#'   actually intend to predict onto.  Required for \code{method = "nndm"}.
#'   The grid from \code{\link{predict_surface}()} is the natural choice; a
#'   non-POINT layer (grid cells, polygons) is reduced to representative points
#'   first, so the target distances are point-to-point.  A point-to-polygon
#'   distance is zero for every cell that contains a training point, which
#'   pulls the target distribution towards zero and degenerates the CV towards
#'   plain leave-one-out.
#' @param min_train For \code{method = "nndm"}: the smallest fraction of the
#'   data any fold's training set may be reduced to by neighbour exclusion.
#'   Default \code{0.5}, as in \code{CAST::nndm()}.
#' @param phi For \code{method = "nndm"}: the distance up to which the two
#'   nearest-neighbour distance distributions are matched, in the CRS the
#'   folds are built in; the exclusion never pushes a held-out point's
#'   nearest neighbour beyond it.  In Mila et al. (2022), and in
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
#' The procedure is the paper's own (as in \code{CAST::nndm()}), and it is
#' deterministic.  Let \eqn{G_{ij}} be the empirical distribution of
#' prediction-to-nearest-training distances and \eqn{G_j^*} the distribution
#' of each held-out point's nearest remaining training point.  Starting from
#' plain leave-one-out, the point with the smallest \eqn{G_j^*} at which the
#' realised distribution exceeds the target (\eqn{G_j^*(r) > G_{ij}(r)}) has
#' its nearest training neighbour removed, and this repeats until no such
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
#' distribution exceeded the target by up to 0.17 (13% of folds had a
#' nearest training point within 50 m against a target of 9%), an
#' \emph{optimistic} cross-validation.  \code{params$max_ecdf_excess} reports
#' the largest remaining excess; compare \code{params$target_median} with
#' \code{params$realised_median} as well.
#' @references
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
#' @param drop_empty_blocks Logical. Default TRUE.
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
#' Only \code{"block_kfold"} balances the number of points per fold, and it
#' does so by packing: blocks are taken largest first and each goes to the
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
#'   For the methods that work in projected space (\code{"block_kfold"},
#'   \code{"buffered_loo"} and \code{"nndm"}), \code{params} carries a
#'   \code{params$blocks_supplied} that says whether the blocks came from
#'   \code{blocks} or from a grid built here, and
#'   \code{params$boundary_supplied} whether a \code{boundary} was given;
#'   \code{params$row_probe} is a small sample of row IDs and coordinates
#'   that every \code{cv_*()} compares against the data it is handed, so
#'   folds built from a different layer of the same size are refused, never
#'   applied silently.
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
                       drop_empty_blocks = TRUE, blocks = NULL,
                       balance_tol = 3) {
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
  # block_size was tested with `is.numeric(block_size) && block_size > 0` and
  # anything failing that was silently ignored -- yet echoed back unchanged in
  # params$block_size, so a negative, zero or character value looked honoured.
  # NA and a length-2 vector reached the grid arithmetic and died as internal
  # R errors.  Validate it once, the way `k` is.
  if (!is.null(block_size) &&
      (!is.numeric(block_size) || length(block_size) != 1L ||
       !is.finite(block_size) || block_size <= 0))
    stop("make_folds(): `block_size` must be a single positive number in the ",
         "units of the data's CRS; got ",
         paste(format(block_size), collapse = ", "), ".", call. = FALSE)
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
        blocks <- ensure_projected(blocks, .crs_or_null(pts))
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
      b <- ensure_projected(boundary, .crs_or_null(pts))
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
      if (is.finite(sac_range) && sac_range > 0 && blocks_supplied) {
        # Supplied blocks cannot be resized; the range still feeds the
        # leakage diagnostic made on the blocks further down.
        message(sprintf(
          "make_folds(block_kfold): estimated spatial autocorrelation range = %.1f CRS units; the supplied `blocks` are used as given, so the range is only compared against them.",
          sac_range
        ))
      } else if (is.finite(sac_range) && sac_range > 0) {
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

    bb <- sf::st_bbox(reg)

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
          .log_warn(
            "make_folds(block_kfold): supplied block_size (%.1f) is smaller than the estimated autocorrelation range (%.1f). Spatial CV may leak correlated information.",
            block_size, sac_range
          )
          warning(
            sprintf("make_folds(): block_size (%.1f) < estimated autocorrelation range (%.1f). Consider increasing block_size to reduce information leakage across folds.",
                    block_size, sac_range),
            call. = FALSE
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
        # %s, not %d: nx and ny are doubles from floor(), and a block_size in
        # the wrong unit -- the very mistake this guard exists to explain --
        # gives counts past 2^31 that %d refuses with "invalid format".
        stop(sprintf(paste0("make_folds(block_kfold): the requested grid is %s x ",
                            "%s = %s cells, above the %s this function will ",
                            "build. Check that `block_size` (%s) is expressed in ",
                            "the data's CRS units (%s) over an extent of %s x %s; ",
                            "a value in the wrong unit is the usual cause."),
                    format(nx, scientific = FALSE), format(ny, scientific = FALSE),
                    format(n_cells_est, big.mark = ",", scientific = FALSE),
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
    } else {
      # The caller's polygons are the blocks.  Their row order is the block
      # id, so `assignment` can be joined back to the layer that was passed.
      grid_sf <- sf::st_as_sf(blocks)
      nx <- NA_integer_; ny <- NA_integer_
    }
    n_blocks <- nrow(grid_sf)
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
    block_source_row <- seq_len(nrow(grid_sf))
    if (drop_empty_blocks) {
      used_blocks <- sort(unique(pts$..block_id[!is.na(pts$..block_id)]))
      grid_sf <- grid_sf[used_blocks, , drop = FALSE]
      pts$..block_id <- match(pts$..block_id, used_blocks)
      block_source_row <- used_blocks
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
      # Three ways to arrive here; name the one the caller actually used, not
      # an argument they never passed.
      if (blocks_supplied)
        stop(sprintf(paste0("make_folds(block_kfold): all %d points fall in the ",
                            "same one of the supplied `blocks`, so there is no ",
                            "spatial split to make and the one fold would have ",
                            "an empty training set. Supply smaller blocks."),
                     nrow(pts)), call. = FALSE)
      how <- if (!is.null(block_size))
        sprintf("the block size (%s)", format(block_size))
      else if (!is.null(block_nx) || !is.null(block_ny))
        "the requested block_nx/block_ny"
      else
        "the automatic grid (block_multiplier x k blocks over the extent)"
      stop(sprintf(paste0("make_folds(block_kfold): %s produces a single block ",
                          "covering the whole extent, so there is no spatial ",
                          "split to make and the one fold would have an empty ",
                          "training set. Pass a smaller `block_size` (or more ",
                          "blocks), or set auto_range = FALSE if the size came ",
                          "from the estimated autocorrelation range."), how),
           call. = FALSE)
    }
    if (B < k) { .log_warn("make_folds(block_kfold): blocks < k; reducing k."); k <- B }

    # Counted over every block in the grid, not over 1..B.  B is the highest
    # block id a point fell in, which under drop_empty_blocks = FALSE says
    # nothing about how many blocks there are: on a layer whose points sit in
    # one quadrant of an 8x8 grid, B was 27, so blocks 28-64 were packed into
    # no fold at all while the 12 equally empty blocks below 27 were -- the
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
        .log_warn(
          "make_folds(block_kfold): the supplied blocks have a median scale (sqrt of area) of %.1f units, smaller than the estimated autocorrelation range (%.1f).",
          block_scale, sac_range
        )
        warning(
          sprintf("make_folds(): the supplied blocks (median scale %.1f) are smaller than the autocorrelation range (%.1f). Spatial CV may leak correlated information across folds; supply blocks at least %.0f units across.",
                  block_scale, sac_range, ceiling(sac_range)),
          call. = FALSE
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
    order_blk <- order(blk_sizes, decreasing = TRUE)
    fold_loads <- integer(k); fold_blocks <- vector("list", k)
    for (i in seq_along(order_blk)) {
      j <- which(fold_loads == min(fold_loads))
      if (length(j) > 1) j <- sample(j, 1L)
      fold_blocks[[j]] <- c(fold_blocks[[j]], order_blk[i])
      fold_loads[j] <- fold_loads[j] + blk_sizes[order_blk[i]]
    }

    # The residual imbalance is checked against the tolerance and, past it,
    # raised as a warning a pipeline can catch -- not only logged.
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
      violates <- is.finite(r) && (cnt / n) > G_target(r) + 1e-12 &&
                  r <= phi && (n - 1L - removed[j]) > rmin &&
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
#' @param k Number of folds. Default 5.
#' @param seed RNG seed. Default 123.
#' @param adaptive Logical; use adaptive bandwidth. Default TRUE.
#' @param bandwidth Optional bandwidth, applied to every fold. For
#'   \code{adaptive = TRUE} an integer number of neighbours; for
#'   \code{adaptive = FALSE} a distance in the units of the \strong{projected}
#'   CRS the folds are fitted in (geographic input is projected first, so a
#'   value in degrees is read as metres). \code{NULL} (default) selects a
#'   bandwidth per fold with \code{GWmodel::bw.gwr()}, which is reported in
#'   \code{fold_metrics$bandwidth}. See \code{\link{fit_gwr_model}}.
#' @param kernel Kernel function type.
#' @param boundary Optional polygonal sf/sfc for CRS alignment.
#' @param pointize Geometry coercion strategy.
#' @param block_size Optional minimum block edge length for spatial CV blocks
#'   (projected CRS units).  Blocks are then at least as large as the
#'   spatial autocorrelation range.
#' @param auto_range Logical.  If \code{TRUE} and \code{folds} is \code{NULL},
#'   estimate the autocorrelation range and use it as the minimum block size.
#'   Default \code{FALSE}.
#' @param parallel Accepted so that every \code{cv_*()} function takes the
#'   same arguments, but the GWR folds always run one after another in this
#'   R process.  GWmodel is built with OpenMP, and OpenMP (GNU libgomp)
#'   deadlocks \code{parallel::mclapply()}'s forked workers once a GWR has
#'   been fitted in the session, for example by \code{fit_gwr_model()}, so
#'   forking would hang the call.  A value asking for more than one core
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
#'   not its predictor count (see \code{\link{cv_spatial}}).
#' @family cross-validation
#' @examples
#' if (requireNamespace("GWmodel", quietly = TRUE) &&
#'     requireNamespace("sp", quietly = TRUE)) {
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
  if (!requireNamespace("GWmodel", quietly = TRUE))
    stop("cv_gwr(): package 'GWmodel' is required.", call. = FALSE)
  if (!requireNamespace("sp", quietly = TRUE))
    stop("cv_gwr(): package 'sp' is required (for GWmodel interop).", call. = FALSE)

  kernel <- match.arg(kernel)
  kernel <- .validate_kernel(kernel)

  if (!("..row_id" %in% names(data_sf))) data_sf$`..row_id` <- seq_len(nrow(data_sf))
  folds <- .folds_from_labels(folds, data_sf, "cv_gwr")
  dat_sf <- prep_model_data(data_sf, response_var, predictor_vars, boundary, pointize)
  if (!("..row_id" %in% names(dat_sf)))
    stop("cv_gwr(): `prep_model_data()` must preserve `..row_id`.")

  keep_idx <- dat_sf$`..row_id`

  # Refuse folds built from other data before anything is fitted: the splits
  # are row IDs, so a wrong `folds` of the right size applies silently.
  .check_fold_probe(folds, data_sf, "cv_gwr")

  if (is.null(folds)) {
    message("cv_gwr(): no folds supplied -- using spatial block k-fold CV (k=", k, ").")
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
  #
  # GWR folds never fork.  GWmodel is built with OpenMP, and GNU libgomp is
  # not fork-safe: once any GWR has run in this session (fit_gwr_model(), a
  # sequential cv_gwr(), select_gwr_variables(), a bare GWmodel::bw.gwr()) the
  # parent holds an OpenMP thread pool that a forked child inherits without
  # its threads, so the child's first parallel region waits on a futex for
  # ever and parallel::mclapply() never returns.  That was the ordinary
  # fit-then-cross-validate order, and it hung with no timeout.  Nothing
  # visible from R says whether the pool exists, so the folds run here, one
  # after another.  A PSOCK cluster would avoid the fork, but its workers
  # need this same spatialkit installed, which pkgload::load_all() or any
  # development copy does not give them, and they cannot see what a
  # `metrics` closure reads from the global environment.  The other cv_*()
  # keep mclapply(): ranger runs its own threads, not libgomp's.
  if (.resolve_n_cores(parallel) > 1L) {
    .warn_and_log(paste0(
      "cv_gwr(): `parallel` is ignored and the %d folds run one after ",
      "another. GWmodel runs OpenMP code, which deadlocks forked (mclapply) ",
      "workers once a GWR has been fitted in the session."),
      length(remapped_folds))
    parallel <- FALSE
  }
  res <- .cv_run_folds(
    dat_sf = dat_sf, response_var = response_var,
    predictor_vars = predictor_vars,
    remapped_folds = remapped_folds, keep_idx = keep_idx,
    fit_one = fit_one, fold_info_fn = fold_info_fn,
    p = NULL, parallel = parallel, seed = seed, metrics = metrics
  )

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
    .log_warn("cv_gwr(): all %d folds failed to produce predictions; results are empty.%s",
              n_attempted, why)
    warning("cv_gwr(): all folds failed; cross-validation results contain no predictions.",
            why, call. = FALSE)
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
#' @param k Number of folds. Default 5.
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
#' @param metrics Optional scoring function of your own, a
#'   \code{function(y, yhat)} returning a named numeric vector; its names
#'   become columns of \code{fold_metrics} (per fold) and \code{overall}
#'   (pooled) beside the built-in ones.  See \strong{Your own metrics} on
#'   \code{\link{cv_spatial}()} for the contract.  \code{yhat} is the
#'   posterior predictive \code{summary} (mean or median) of each held-out
#'   row; a score that needs the draws belongs in \code{predictive_coverage}
#'   and \code{CRPS}, which this function computes itself.
#' @inheritSection model_metrics Percentage errors on responses with zeros
#' @inheritSection model_metrics Which metrics survive a non-Gaussian response
#' @return A list with \code{overall}, \code{fold_metrics},
#'   \code{predictions}, \code{folds}, \code{n_folds_attempted},
#'   \code{n_folds_succeeded}, \code{fold_status}, \code{orphan_rows},
#'   \code{n_unknown_ids}, \code{n_dropped}, \code{formula} and
#'   \code{predictive_coverage}.
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
  metrics <- .check_metrics_fn(metrics, "cv_bayes")

  if (!("..row_id" %in% names(data_sf))) data_sf$`..row_id` <- seq_len(nrow(data_sf))
  folds <- .folds_from_labels(folds, data_sf, "cv_bayes")
  dat_sf <- prep_model_data(data_sf, response_var, predictor_vars, boundary, pointize)
  if (!("..row_id" %in% names(dat_sf)))
    stop("cv_bayes(): `prep_model_data()` must preserve `..row_id`.")

  keep_idx <- dat_sf$`..row_id`
  n_pred   <- length(predictor_vars)

  # Refuse folds built from other data before anything is fitted: the splits
  # are row IDs, so a wrong `folds` of the right size applies silently.
  .check_fold_probe(folds, data_sf, "cv_bayes")

  if (is.null(folds)) {
    message("cv_bayes(): no folds supplied -- using spatial block k-fold CV (k=", k, ").")
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
        # Per-row posterior predictive SD, for $predictions$yhat_sd.  The
        # column existed as an unconditional NA placeholder; the draws that
        # give the coverage below give this for free.
        if (nrow(ppred_draws) >= 2L && ncol(ppred_draws) == length(y_true))
          extras$..per_row <- data.frame(
            yhat_sd = apply(ppred_draws, 2L, stats::sd),
            stringsAsFactors = FALSE)

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
    p = NULL, parallel = parallel, seed = seed, metrics = metrics
  )

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
               gp_n_basis = integer(), n_draws = integer())
  if (!nrow(folds_df))
    for (cn in .user_metric_names(metrics)) folds_df[[cn]] <- numeric()

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
  } else {
    .cv_warn_failed_folds("cv_bayes", res, preds, length(keep_idx),
                          n_attempted, n_succeeded)
  }

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
#' @param k Number of folds.
#' @param seed RNG seed.
#' @param boundary Optional boundary for fold construction.
#' @param pointize Geometry coercion strategy.
#' @param predict_args Extra arguments for predict().
#' @param fold_info_fn Optional \code{function(fit, test_sf, y, yhat)}
#'   returning a named list of per-fold extras (a bandwidth, a tuning value,
#'   anything read off the fitted object), added as columns of
#'   \code{fold_metrics}.  It sees the fit and the held-out layer, which
#'   \code{metrics} does not; it is applied per fold only, and its values
#'   are not pooled.  An element \code{..per_row} that is a data frame with
#'   one row per held-out observation is spliced into \code{predictions}
#'   instead.
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
#'   \code{FALSE} (sequential).  A learner that runs OpenMP code (GWmodel,
#'   or an xgboost built with GNU libgomp) can hang the forked workers once
#'   it has run in the session; keep such a \code{fit_fn} sequential, as
#'   \code{\link{cv_gwr}()} does.
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
#' built-in column names, one number per name.  Anything else is an error,
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
#'   \code{status} is \code{"ok"}; \code{"error"} (the fit or its
#'   \code{predict()} threw; \code{message} is the error text);
#'   \code{"skipped"} (nothing scorable: too few matched rows, a prediction
#'   of the wrong length, or no finite observed/predicted pair);
#'   \code{"dropped"} (an empty test set, or fewer than two training rows,
#'   once incomplete rows were removed, so the fold never reached the fitter);
#'   or \code{"worker_error"} (a parallel worker died).  Every fold missing
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
#'   were unusable and dropped.  \code{overall$Adj_R2} is always \code{NA}: the
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
  dat_sf <- prep_model_data(data_sf, response_var, predictor_vars, boundary, pointize)
  keep_idx <- dat_sf$`..row_id`

  # Refuse folds built from other data before anything is fitted: the splits
  # are row IDs, so a wrong `folds` of the right size applies silently.
  .check_fold_probe(folds, data_sf, .caller)

  if (is.null(folds)) {
    message(.caller, "(): no folds supplied -- using spatial block k-fold CV (k=", k, ").")
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
    parallel = parallel, seed = seed, metrics = metrics
  )

  preds <- if (length(res$pred_rows)) do.call(rbind, res$pred_rows) else
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
    .log_warn("%s(): all %d folds failed to produce predictions; results are empty.%s",
              .caller, n_attempted, why)
    warning(.caller, "(): all folds failed; cross-validation results contain ",
            "no predictions.", why, call. = FALSE)
  } else {
    .cv_warn_failed_folds(.caller, res, preds, length(keep_idx),
                          n_attempted, n_succeeded)
  }

  list(overall = .cv_overall_metrics(preds, metrics),
       fold_metrics = folds_df, predictions = preds,
       folds = .bare_folds(remapped_folds),
       n_folds_attempted = n_attempted, n_folds_succeeded = n_succeeded,
       fold_status = .cv_fold_status(remapped_folds, res),
       orphan_rows = attr(remapped_folds, "orphans"),
       n_unknown_ids = attr(remapped_folds, "n_unknown_ids"),
       n_dropped = as.integer(.get_row_record(dat_sf, "dropped")$n %||% 0L))
}
