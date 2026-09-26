# =============================================================================
# Forward feature selection with spatially blocked inner folds
# =============================================================================

#' Greedy forward feature selection with spatially blocked inner folds
#'
#' Selects predictors by repeatedly adding whichever candidate most improves a
#' cross-validated score, stopping when no candidate improves it by more than
#' \code{tol}.
#'
#' \strong{The inner folds must be spatial, and that is the entire point.}
#' Nested selection is only worth doing if the inner loop is blocked the same
#' way the outer one is.  Random inner folds inside blocked outer folds select
#' variables that look predictive only because nearby points leak between
#' train and test.  The outer loop then reports honest-looking numbers for a
#' dishonestly chosen feature set, which is worse than not selecting at all,
#' because the dishonesty is now hidden behind a defensible-looking validation.
#' \code{method} therefore defaults to \code{"block_kfold"} and logs a loud
#' caution if set to \code{"random_kfold"} (a deliberate choice, so it is not
#' raised as an R warning).
#'
#' Call this \emph{inside} the \code{fit_fn} you pass to \code{cv_spatial()}.
#' \code{.cv_fit_one_fold()} calls \code{fit_fn(train_sf)} on the training
#' slice only, so anything done inside it is automatically nested and
#' leak-free; no extra plumbing is needed.  The cost grows fast: a sweep over
#' \code{p} candidates costs roughly \code{p^2 / 2 * k} model fits, and nesting
#' that inside \code{n} outer leave-one-out folds multiplies it by \code{n}.
#' \code{max_fits} guards against that.
#'
#' @section The score is not a performance estimate:
#' \code{$score} is the cross-validated \code{metric} of the winning set at
#' the final step: the best of every candidate set the sweep scored.  That
#' is the number the selection optimised, and a number optimised over many
#' candidates is optimistically biased by construction: Cawley and Talbot
#' (2010) show the bias can exceed the genuine differences between the models
#' being compared.  Quote it as the selection criterion, not as the
#' performance of the selected model.  An honest performance estimate needs
#' data the selection never saw.  Two ways to get one: run this function
#' inside the \code{fit_fn} of \code{\link{cv_spatial}()}, so the outer
#' folds score a model whose predictors were chosen on the inner ones alone;
#' or pass \code{select_on = "split"}, which selects on one spatial half of
#' \code{train_sf} and returns the selected set's score on the other as
#' \code{score_holdout}.
#'
#' @param train_sf Training data (\code{sf}).
#' @param response_var Character(1).
#' @param candidate_vars Character vector of predictors to choose among.
#' @param fit_fn A function \code{(train_sf, predictor_vars)} returning a
#'   \code{spatial_fit}.  The signature takes two arguments because selection
#'   has to refit with different predictor sets.
#' @param k Inner fold count. Default 5.
#' @param method Inner fold method. Default \code{"block_kfold"}.
#' @param block_size Passed to \code{make_folds()}; inherit the outer block
#'   size so inner and outer blocks are on the same spatial scale.
#' @param metric Score to optimise: \code{"RMSE"}, \code{"MAE"} (minimised) or
#'   \code{"R2"} (maximised). Default \code{"RMSE"}.
#' @param tol Minimum improvement required to accept a variable.  Default 0,
#'   meaning any improvement is accepted. The first variable is judged against
#'   the null (intercept-only) model, so \code{tol} bites from step 1, but
#'   only when that null model can be scored. Backends that refuse a
#'   zero-length \code{predictor_vars} (\code{\link{fit_rf_model}} and
#'   \code{\link{fit_gwr_model}} both do) have no null score, and there the
#'   first variable is accepted unconditionally.
#' @param max_vars Optional cap on how many predictors to select.
#' @param max_fits Abort if the sweep would exceed this many model fits.
#'   Default 5000.
#' @param seed RNG seed.  It governs both the inner fold construction and the
#'   cross-validation itself: it is forwarded to \code{cv_spatial(seed = )},
#'   which draws one RNG stream per fold from it, so it also seeds the
#'   \emph{learner} inside every fold.  A stochastic \code{fit_fn} is therefore
#'   reproducible from this one value.
#' @param quiet Logical; suppress this function's progress \code{message()}s.
#'   It does not silence R warnings, nor the package's console log echo
#'   (see \code{\link{spatialkit_quiet}} for that). Default \code{FALSE}.
#' @param auto_range Logical.  If \code{TRUE} and \code{method} is
#'   \code{"block_kfold"}, the autocorrelation range is estimated from the
#'   response (detrended on the candidates) and used as the minimum block
#'   size of the inner folds, exactly as in \code{\link{make_folds}()};
#'   \code{block_size} still applies as a floor.  Default \code{FALSE}, which
#'   keeps the geometric blocks.  Inner blocks smaller than the range let a
#'   candidate be selected for spatial proximity to the response rather than
#'   for predicting it, which is the same failure the \code{random_kfold}
#'   caution above exists to prevent, so the leakage warning
#'   \code{make_folds()} raises applies here with more force than usual.
#' @param select_on \code{"all"} (default) runs the sweep on every row of
#'   \code{train_sf}.  \code{"split"} runs it on one spatially blocked half,
#'   then fits the selected set on that half and scores it on the other:
#'   \code{score_holdout} is then the selected model's \code{metric} on rows
#'   whose response the sweep never read (the sweep's own \code{score} is
#'   not; see "The score is not a performance estimate").  It is one
#'   estimate from one region: the halves share a border with no buffer, so
#'   rows near it are still correlated with the selection half.
#'   Both halves come back in \code{$split}.  See the "Post-selection
#'   inference" section of \code{\link{determine_optimal_levels}} for the
#'   trade: coverage for half the sample.
#' @return A list of class \code{"feature_selection"} (so that
#'   \code{\link{plot.feature_selection}()} draws the selection path) with
#'   \code{selected} (the chosen predictors, in the order they were added),
#'   \code{score}, \code{score_holdout}, \code{history}, \code{params} and
#'   \code{split}.
#'   \code{score} is the winning set's cross-validated \code{metric} at the
#'   final step: the \strong{selection-internal} optimum, optimistically
#'   biased because it was chosen as the best of many (see the section above),
#'   and \code{NA} when nothing was selected.  \code{history} is a data.frame
#'   with \code{step}, \code{variable}, \code{score} and \code{n_pred},
#'   holding every candidate evaluated at every step; when the null model
#'   could be scored it also carries a \code{step = 0} row named
#'   \code{"<none>"} giving that baseline, so the first variable's gain can be
#'   read off directly.  Every set is scored on the same rows: those the null
#'   model's cross-validation predicted or, when there is no null model,
#'   those any step-1 set predicted; \code{params$n_scored} counts them (a
#'   warning says so when that is fewer than all).  \code{n_pred} is how many
#'   rows the set's cross-validation predicted.  A set that left some of the
#'   scored rows unpredicted, because a fold failed for it, has \code{score}
#'   \code{NA}, with a warning naming it: scored on the rows it did predict
#'   it would be compared on fewer, usually easier, rows than its rivals.  A
#'   factor with a level found in one spatial block only is the usual case,
#'   and cannot be selected.
#'   \code{score_holdout} is \code{NA} unless \code{select_on = "split"}, and
#'   then the selected set's \code{metric} when fitted on the selection half
#'   and predicted on the estimation half (\eqn{R^2} against the selection
#'   half's mean, the out-of-sample convention); \code{NA} when nothing was
#'   selected or the prediction failed.  \code{split} is \code{NULL} or a
#'   list with \code{selection} and \code{estimation}, integer row positions
#'   in \code{train_sf} as passed; rows the completeness filter above dropped
#'   are in neither.
#'   \code{params} records \code{metric}, \code{method}, \code{k},
#'   \code{tol}, \code{seed}, \code{auto_range}, \code{select_on},
#'   \code{n_candidates}, \code{estimated_fits} and \code{n_scored}.
#' @references
#' Cawley, G. C. and Talbot, N. L. C. (2010). On over-fitting in model
#' selection and subsequent selection bias in performance evaluation.
#' \emph{Journal of Machine Learning Research}, 11, 2079--2107.
#' \url{https://jmlr.org/papers/v11/cawley10a.html}
#' @family cross-validation
#' @examples
#' if (requireNamespace("GWmodel", quietly = TRUE) &&
#'     requireNamespace("sp", quietly = TRUE)) {
#'   library(sf)
#'   set.seed(1)
#'   n <- 120
#'   pts <- st_as_sf(
#'     data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
#'                a = rnorm(n), b = rnorm(n), noise = rnorm(n)),
#'     coords = c("x", "y"), crs = 32632
#'   )
#'   pts$resp <- 3 * pts$a + 2 * pts$b + rnorm(n, 0, 0.3)
#'
#'   fit_fn <- function(tr, vars) fit_gwr_model(tr, "resp", vars, bandwidth = 30)
#'   sel <- select_features_forward(pts, "resp", c("a", "b", "noise"), fit_fn,
#'                                  k = 3, quiet = TRUE)
#'   sel$selected
#' }
#' @export
select_features_forward <- function(train_sf, response_var, candidate_vars,
                                    fit_fn, k = 5,
                                    method = c("block_kfold", "random_kfold"),
                                    block_size = NULL,
                                    metric = c("RMSE", "MAE", "R2"),
                                    tol = 0, max_vars = NULL,
                                    max_fits = 5000L, seed = 123,
                                    quiet = FALSE, auto_range = FALSE,
                                    select_on = c("all", "split")) {
  method <- match.arg(method)
  metric <- match.arg(metric)
  select_on <- match.arg(select_on)
  .msg <- function(...) if (!quiet) message(...)

  if (!inherits(train_sf, "sf"))
    stop("select_features_forward(): `train_sf` must be an sf object.", call. = FALSE)
  if (!is.function(fit_fn))
    stop("select_features_forward(): `fit_fn` must be a function of ",
         "(train_sf, predictor_vars).", call. = FALSE)
  candidate_vars <- unique(candidate_vars)
  missing_v <- setdiff(candidate_vars, names(train_sf))
  if (length(missing_v) > 0L)
    stop("select_features_forward(): candidate(s) absent from `train_sf`: ",
         paste(sQuote(missing_v), collapse = ", "), call. = FALSE)
  if (length(candidate_vars) == 0L)
    stop("select_features_forward(): no candidate variables supplied.", call. = FALSE)

  # Every candidate set must be scored on the SAME observations.  cv_spatial()
  # re-runs prep_model_data() per candidate, which drops rows that are NA in
  # THAT candidate's columns -- so without this, a step compares an RMSE over
  # 150 rows against one over 120, and can prefer a variable for having an
  # easier surviving subset rather than for predicting better.
  #
  # The test must match prep_model_data()'s exactly.  complete.cases() alone
  # lets Inf through, while prep_model_data() also drops on all(is.finite(r))
  # over the numeric columns -- so a candidate carrying a single Inf would
  # still reach cv_spatial() with a smaller row set than its rivals, which is
  # precisely the comparison this filter exists to prevent.
  fs_df <- sf::st_drop_geometry(train_sf)[, c(response_var, candidate_vars),
                                          drop = FALSE]
  keep <- stats::complete.cases(fs_df)
  num_mask <- vapply(fs_df, is.numeric, logical(1))
  if (any(num_mask))
    keep <- keep & apply(as.matrix(fs_df[, num_mask, drop = FALSE]), 1L,
                         function(r) all(is.finite(r)))
  if (!all(keep)) {
    if (sum(keep) < 2L)
      stop("select_features_forward(): fewer than two rows are complete and ",
           "finite across the response and every candidate.", call. = FALSE)
    .log_warn(paste0("select_features_forward(): dropping %d row(s) incomplete ",
                     "or non-finite across the response and candidates, so ",
                     "every candidate set is scored on the same %d ",
                     "observations."),
              sum(!keep), sum(keep))
    train_sf <- train_sf[keep, , drop = FALSE]
  }

  # Sample splitting: the sweep sees the selection half only; the estimation
  # half scores the chosen set afterwards.  .spatial_half_split()'s positions
  # index train_sf as it stands here, after the completeness filter.
  split <- NULL
  holdout_sf <- NULL
  if (identical(select_on, "split")) {
    if (!all(sf::st_geometry_type(train_sf, by_geometry = TRUE) == "POINT"))
      train_sf <- coerce_to_points(train_sf, "auto")
    split <- .spatial_half_split(train_sf, seed = seed,
                                 caller = "select_features_forward")
    holdout_sf <- train_sf[split$estimation, , drop = FALSE]
    train_sf   <- train_sf[split$selection, , drop = FALSE]
    # Returned as positions in the layer the CALLER passed, as print() says
    # and as determine_optimal_levels() returns them.  They were positions
    # after the completeness filter, so with 7 incomplete rows
    # pts[fs$split$estimation, ] held 47 selection-half rows and 4 of the
    # dropped ones: the one use of the split, estimating on rows the
    # selection never saw, got rows it had seen.
    kept <- which(keep)
    split$selection  <- kept[split$selection]
    split$estimation <- kept[split$estimation]
    .msg(sprintf("select_features_forward(): selecting on %d points, scoring the result on the other %d.",
                 nrow(train_sf), nrow(holdout_sf)))
  }

  if (identical(method, "random_kfold"))
    .log_warn(paste0("select_features_forward(): random_kfold inner folds ",
                     "leak spatial autocorrelation between train and test. ",
                     "Variables can then be selected for being spatially ",
                     "close to the response rather than predictive of it, and ",
                     "an outer blocked loop will report honest-looking numbers ",
                     "for a dishonestly chosen feature set. Use ",
                     "'block_kfold' unless you know why you are not."))

  p <- length(candidate_vars)
  max_steps <- if (is.null(max_vars)) p else min(p, as.integer(max_vars))
  # sweep cost: p + (p-1) + ... over at most max_steps rounds, times k folds
  est_fits <- sum(p - seq_len(max_steps) + 1L) * k
  if (est_fits > max_fits)
    stop(sprintf(paste0("select_features_forward(): this sweep would fit about ",
                        "%d models (%d candidates, %d folds), above max_fits = %d. ",
                        "Reduce candidates, lower k, set max_vars, or raise ",
                        "max_fits deliberately. Nesting a sweep inside ",
                        "leave-one-out outer folds multiplies this by n."),
                 est_fits, p, k, max_fits), call. = FALSE)

  better <- if (metric == "R2") function(a, b) a > b else function(a, b) a < b
  worst  <- if (metric == "R2") -Inf else Inf

  # Build the inner folds ONCE, before the sweep.
  #
  # Fold construction does not depend on the candidate set: make_folds() reads
  # `predictor_vars` only to detrend the range estimate (under auto_range, and
  # for the leakage diagnostic), which uses the full candidate list whatever
  # the sweep has selected so far, and every other argument is fixed across
  # the sweep.  Rebuilding them inside score_set() therefore produced the same
  # splits p^2/2 times -- repeating every block-size warning as many times --
  # and, worse, put make_folds() OUTSIDE the try() below, so a single
  # fold-construction failure (block_kfold raising when the geometry collapses
  # to one block, say) killed the whole sweep instead of the candidate being
  # scored NA.  It cannot be a per-candidate NA in any case: if the folds
  # cannot be built, no candidate is scorable, so this is one informative
  # error instead of p^2/2 silent ones.
  folds <- try(make_folds(train_sf, k = k, method = method, seed = seed,
                          block_size = block_size, auto_range = auto_range,
                          response_var = response_var,
                          predictor_vars = candidate_vars), silent = TRUE)
  if (inherits(folds, "try-error"))
    stop("select_features_forward(): could not build the inner CV folds, so no ",
         "candidate can be scored: ",
         trimws(conditionMessage(attr(folds, "condition"))),
         " Adjust `method`, `k` or `block_size`.", call. = FALSE)

  # Every set is scored on ONE fixed row set, `ref_ids`: the rows the null
  # model predicted, or, when there is no null model (RF and GWR refuse an
  # empty predictor set), the rows any step-1 candidate predicted.
  # cv_spatial()'s `overall` pools whichever folds survived, so a set whose
  # fit or predict failed on a fold -- a factor level found in one block
  # only, the ordinary case -- was scored on fewer rows, usually the easier
  # ones, and could win for that alone: a noise factor was chosen over the
  # true driver (RF RMSE 2.35 on 192 rows against 2.63 on 250).  A set that
  # leaves a reference row unpredicted is scored NA instead, and said so.
  # Rows the reference runs did not predict (a fold that fails for every set,
  # from geometry rather than predictors) drop out of every score alike.
  ref_ids <- NULL
  run_set <- function(vars) {
    inner_fit <- function(tr) fit_fn(tr, vars)
    # cv_spatial()'s own partial-failure warning is muffled: the consequence
    # for the sweep is reported below, once per step, in the sweep's terms.
    res <- try(withCallingHandlers(
      suppressMessages(
        cv_spatial(train_sf, response_var, vars, fit_fn = inner_fit,
                   folds = folds, seed = seed)),
      warning = function(w)
        if (.is_failed_folds_warning(w)) invokeRestart("muffleWarning")),
      silent = TRUE)
    if (inherits(res, "try-error") || !is.data.frame(res$predictions))
      return(NULL)
    pr <- res$predictions
    pr[is.finite(pr$y) & is.finite(pr$yhat), , drop = FALSE]
  }
  n_pred_of <- function(pr) if (is.null(pr)) 0L else nrow(pr)
  covers_ref <- function(pr) !is.null(pr) && all(ref_ids %in% pr$`..row_id`)
  score_on <- function(pr) {
    if (!length(ref_ids) || !covers_ref(pr)) return(NA_real_)
    val <- .cv_overall_metrics(pr[pr$`..row_id` %in% ref_ids, , drop = FALSE])[[metric]]
    if (is.null(val) || !is.finite(val)) NA_real_ else as.numeric(val)
  }
  set_ref <- function(ids, from) {
    ref_ids <<- ids
    if (length(ids) && length(ids) < nrow(train_sf))
      .warn_and_log(paste0(
        "select_features_forward(): %s predicted only %d of the %d rows, so ",
        "every candidate set is scored on those %d; the rest drop out of ",
        "every score alike."),
        from, length(ids), nrow(train_sf), length(ids))
  }

  selected  <- character(0)
  remaining <- candidate_vars
  history   <- list()

  # Score the EMPTY set first, so step 1 has something to beat.  `best` used to
  # start at Inf/-Inf with the stopping test gated on is.finite(best), which
  # made the first step unconditional: whatever the winning candidate's score,
  # it was accepted.  On a pure-noise response a "predictive" feature was
  # therefore selected in 100% of runs, and in 23% of them the returned set was
  # WORSE in CV RMSE than fitting nothing at all -- while the documentation
  # says the sweep stops "when no candidate improves it by more than `tol`".
  # An intercept-only fit is not something every backend can do, so fall back
  # to the old unconditional first step when it fails, and say so.
  # A backend that cannot fit an intercept-only model -- fit_rf_model() and
  # fit_gwr_model() both refuse a zero-length predictor_vars -- makes
  # cv_spatial() report "all folds failed", which is expected here and handled
  # by the fallback below.  suppressWarnings() silences the R condition, but
  # the k "fold i fit failed" lines and the "all k folds failed" line are
  # LOGGER records, which no condition handler touches: every successful
  # RF/GWR run printed them, identical to a genuinely failed run's, even with
  # quiet = TRUE.  Raise the console threshold for the probe alone; the file
  # trace (index 1) keeps the lines, where a diagnostic belongs.  The rows
  # the null model predicts, when it can be fitted, are the rows every set
  # is scored on.
  null_run <- logger::with_log_threshold(
    suppressWarnings(run_set(character(0))),
    threshold = logger::FATAL, namespace = "spatialkit", index = 2)
  if (n_pred_of(null_run) > 0L)
    set_ref(sort(unique(null_run$`..row_id`)), "the null (intercept-only) model")
  null_score <- score_on(null_run)
  best <- if (is.finite(null_score)) null_score else worst
  if (!is.finite(null_score))
    .msg("select_features_forward(): the null (intercept-only) model could not ",
         "be scored, so the first variable is accepted unconditionally.")
  else
    history[[length(history) + 1L]] <- data.frame(
      step = 0L, variable = "<none>", score = unname(null_score),
      n_pred = n_pred_of(null_run), stringsAsFactors = FALSE
    )

  repeat {
    if (length(remaining) == 0L || length(selected) >= max_steps) break

    runs <- lapply(remaining, function(v) run_set(c(selected, v)))
    if (is.null(ref_ids))
      set_ref(sort(unique(unlist(lapply(runs, `[[`, "..row_id")))),
              "the step-1 candidate sets together")
    step_scores <- vapply(runs, score_on, numeric(1))
    names(step_scores) <- remaining
    n_pred <- vapply(runs, n_pred_of, integer(1))
    history[[length(history) + 1L]] <- data.frame(
      step = length(selected) + 1L, variable = remaining,
      score = unname(step_scores), n_pred = n_pred, stringsAsFactors = FALSE
    )
    # Sets that predicted some reference rows but not all.  One predicting
    # none has already raised cv_spatial()'s "all folds failed" warning.
    short <- which(n_pred > 0L & !vapply(runs, covers_ref, logical(1)))
    if (length(short))
      .warn_and_log(paste0(
        "select_features_forward(): step %d: %s left some of the %d scored ",
        "rows unpredicted (a fold failed: a factor level found in one block ",
        "only, say), so %s scored NA rather than on fewer, easier rows."),
        length(selected) + 1L,
        paste(vapply(short, function(j) sprintf(
          "{%s} (%d of them predicted)",
          paste(c(selected, remaining[j]), collapse = ", "),
          sum(ref_ids %in% runs[[j]]$`..row_id`)), character(1)),
          collapse = ", "),
        length(ref_ids), if (length(short) == 1L) "it is" else "they are")

    if (all(is.na(step_scores))) {
      .msg("select_features_forward(): every candidate failed to score at step ",
           length(selected) + 1L, "; stopping.")
      break
    }

    # The commonest way to call this wrongly is with a learner written for
    # cv_spatial(), which takes ONE argument: `function(train_sf, ...)` swallows
    # `vars` into the dots and fits the same model whatever it is handed.  The
    # training layer still carries every column, so nothing errors; every
    # candidate simply scores bit-identically to the null model, nothing
    # improves on it, and the function returns an empty selection and an NA
    # score with -- until this -- no word about why.  Two different predictor
    # sets do not produce the same cross-validated metric to the last digit,
    # so the signature is unambiguous.
    if (length(selected) == 0L && length(step_scores) >= 2L &&
        is.finite(null_score) && all(is.finite(step_scores)) &&
        isTRUE(all.equal(unname(step_scores), rep(null_score, length(step_scores)),
                         tolerance = 1e-12)))
      .warn_and_log(paste0(
        "select_features_forward(): every candidate scored exactly what the ",
        "intercept-only model scored (%s = %.6g), so `fit_fn` appears to be ",
        "ignoring its second argument. It must be a function of ",
        "(train_sf, predictor_vars) that fits only the variables it is given; ",
        "a learner written for cv_spatial() takes one argument and needs a ",
        "wrapper."), metric, null_score)

    idx  <- if (metric == "R2") which.max(step_scores) else which.min(step_scores)
    cand <- remaining[idx]
    cand_score <- step_scores[[idx]]

    gain <- if (metric == "R2") cand_score - best else best - cand_score
    if (is.finite(best) && !(gain > tol)) {
      .msg(sprintf("select_features_forward(): no candidate improves %s by more than %g; stopping.",
                   metric, tol))
      break
    }

    selected  <- c(selected, cand)
    remaining <- setdiff(remaining, cand)
    best      <- cand_score
    .msg(sprintf("select_features_forward(): + %s  (%s = %.4f)",
                 cand, metric, cand_score))
  }

  # The honest score: the selected set fitted on the selection half and
  # predicted on the estimation half, which the sweep never saw.
  score_holdout <- NA_real_
  if (!is.null(holdout_sf) && length(selected) > 0L) {
    score_holdout <- tryCatch({
      fit <- fit_fn(train_sf, selected)
      yh  <- as.numeric(do.call(stats::predict, list(object = fit, newdata = holdout_sf)))
      y   <- as.numeric(sf::st_drop_geometry(holdout_sf)[[response_var]])
      if (length(yh) != length(y))
        stop(sprintf("predict() returned %d values for %d rows", length(yh), length(y)))
      y_tr <- as.numeric(sf::st_drop_geometry(train_sf)[[response_var]])
      met <- .compute_reg_metrics(y, yh, p = NULL,
                                  y_train_mean = mean(y_tr[is.finite(y_tr)]))
      v <- met[[metric]]
      if (is.null(v) || !is.finite(v)) NA_real_ else as.numeric(v)
    }, error = function(e) {
      .warn_and_log("select_features_forward(): the hold-out score could not be computed (%s); score_holdout is NA.",
                    conditionMessage(e))
      NA_real_
    })
  }

  # Classed so that plot() finds plot.feature_selection(); the object is still
  # the same list, and `$`, `[[`, names() and printing all behave as before.
  structure(list(
    selected = selected,
    # Never hand back the `worst` sentinel as if it were a score: when nothing
    # was selected there is no score, and Inf / -Inf reads as a real number to
    # any caller that compares it.
    score    = if (length(selected) == 0L) NA_real_ else best,
    score_holdout = score_holdout,
    history  = if (length(history)) do.call(rbind, history) else
      data.frame(step = integer(0), variable = character(0), score = numeric(0),
                 n_pred = integer(0)),
    params   = list(metric = metric, method = method, k = k, tol = tol,
                    seed = seed, auto_range = isTRUE(auto_range),
                    select_on = select_on,
                    n_candidates = p, estimated_fits = est_fits,
                    n_scored = length(ref_ids)),
    split    = split
  ), class = c("feature_selection", "list"))
}
