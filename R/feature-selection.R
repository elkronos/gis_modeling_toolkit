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
#' train and test -- and the outer loop then reports honest-looking numbers for
#' a dishonestly chosen feature set, which is worse than not selecting at all,
#' because the dishonesty is now hidden behind a defensible-looking validation.
#' \code{method} therefore defaults to \code{"block_kfold"} and logs a loud
#' caution if set to \code{"random_kfold"} (a deliberate choice, so it is not
#' raised as an R warning).
#'
#' Call this \emph{inside} the \code{fit_fn} you pass to \code{cv_spatial()}.
#' \code{.cv_fit_one_fold()} calls \code{fit_fn(train_sf)} on the training
#' slice only, so anything done inside it is automatically nested and
#' leak-free; no extra plumbing is needed.  Note the cost: a sweep over
#' \code{p} candidates costs roughly \code{p^2 / 2 * k} model fits, and nesting
#' that inside \code{n} outer leave-one-out folds multiplies it by \code{n}.
#' \code{max_fits} guards against that.
#'
#' @param train_sf Training data (\code{sf}).
#' @param response_var Character(1).
#' @param candidate_vars Character vector of predictors to choose among.
#' @param fit_fn A function \code{(train_sf, predictor_vars)} returning a
#'   \code{spatial_fit}.  Note the two-argument signature: selection needs to
#'   refit with different predictor sets.
#' @param k Inner fold count. Default 5.
#' @param method Inner fold method. Default \code{"block_kfold"}.
#' @param block_size Passed to \code{make_folds()}; inherit the outer block
#'   size so inner and outer blocks are on the same spatial scale.
#' @param metric Score to optimise: \code{"RMSE"}, \code{"MAE"} (minimised) or
#'   \code{"R2"} (maximised). Default \code{"RMSE"}.
#' @param tol Minimum improvement required to accept a variable.  Default 0,
#'   meaning any improvement is accepted. The first variable is judged against
#'   the null (intercept-only) model, so \code{tol} bites from step 1 --- but
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
#' @param quiet Suppress progress messages. Default FALSE.
#' @return A list with \code{selected} (the chosen predictors, in the order
#'   they were added), \code{score} (their cross-validated \code{metric}),
#'   \code{history} and \code{params}. \code{history} is a data.frame with
#'   \code{step}, \code{variable} and \code{score}, holding every candidate
#'   evaluated at every step; when the null model could be scored it also
#'   carries a \code{step = 0} row named \code{"<none>"} giving that
#'   baseline, so the first variable's gain can be read off directly.
#' @family cross-validation
#' @examples
#' \donttest{
#' if (requireNamespace("GWmodel", quietly = TRUE) &&
#'     requireNamespace("sp", quietly = TRUE)) {
#'   library(sf)
#'   set.seed(1)
#'   n <- 120
#'   pts <- st_as_sf(
#'     data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000),
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
#' }
#' @export
select_features_forward <- function(train_sf, response_var, candidate_vars,
                                    fit_fn, k = 5,
                                    method = c("block_kfold", "random_kfold"),
                                    block_size = NULL,
                                    metric = c("RMSE", "MAE", "R2"),
                                    tol = 0, max_vars = NULL,
                                    max_fits = 5000L, seed = 123,
                                    quiet = FALSE) {
  method <- match.arg(method)
  metric <- match.arg(metric)
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
  # `predictor_vars` only when auto_range = TRUE, which this function never
  # enables, and every other argument is fixed across the sweep.  Rebuilding
  # them inside score_set() therefore produced the same splits p^2/2 times --
  # repeating every block-size warning as many times -- and, worse, put
  # make_folds() OUTSIDE the try() below, so a single fold-construction failure
  # (block_kfold raising when the geometry collapses to one block, say) killed
  # the whole sweep instead of the candidate being scored NA.  It cannot be a
  # per-candidate NA in any case: if the folds cannot be built, no candidate is
  # scorable, so this is one informative error instead of p^2/2 silent ones.
  folds <- try(make_folds(train_sf, k = k, method = method, seed = seed,
                          block_size = block_size, response_var = response_var,
                          predictor_vars = candidate_vars), silent = TRUE)
  if (inherits(folds, "try-error"))
    stop("select_features_forward(): could not build the inner CV folds, so no ",
         "candidate can be scored: ",
         trimws(conditionMessage(attr(folds, "condition"))),
         " Adjust `method`, `k` or `block_size`.", call. = FALSE)

  score_set <- function(vars) {
    inner_fit <- function(tr) fit_fn(tr, vars)
    res <- try(suppressMessages(
      cv_spatial(train_sf, response_var, vars, fit_fn = inner_fit,
                 folds = folds, seed = seed)), silent = TRUE)
    if (inherits(res, "try-error") || is.null(res$overall)) return(NA_real_)
    val <- res$overall[[metric]]
    if (is.null(val) || !is.finite(val)) NA_real_ else as.numeric(val)
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
  # suppressWarnings(): a backend that cannot fit an intercept-only model --
  # fit_rf_model() and fit_gwr_model() both refuse a zero-length
  # predictor_vars -- makes cv_spatial() report "all folds failed", which is
  # expected here and handled by the fallback below.  Letting it escape turned
  # a silent, working call into one that warns about a fit the user never
  # asked for.
  null_score <- suppressWarnings(score_set(character(0)))
  best <- if (is.finite(null_score)) null_score else worst
  if (!is.finite(null_score))
    .msg("select_features_forward(): the null (intercept-only) model could not ",
         "be scored, so the first variable is accepted unconditionally.")
  else
    history[[length(history) + 1L]] <- data.frame(
      step = 0L, variable = "<none>", score = unname(null_score),
      stringsAsFactors = FALSE
    )

  repeat {
    if (length(remaining) == 0L || length(selected) >= max_steps) break

    step_scores <- vapply(remaining, function(v) score_set(c(selected, v)),
                          numeric(1))
    history[[length(history) + 1L]] <- data.frame(
      step = length(selected) + 1L, variable = remaining,
      score = unname(step_scores), stringsAsFactors = FALSE
    )

    if (all(is.na(step_scores))) {
      .msg("select_features_forward(): every candidate failed to score at step ",
           length(selected) + 1L, "; stopping.")
      break
    }

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

  list(
    selected = selected,
    # Never hand back the `worst` sentinel as if it were a score: when nothing
    # was selected there is no score, and Inf / -Inf reads as a real number to
    # any caller that compares it.
    score    = if (length(selected) == 0L) NA_real_ else best,
    history  = if (length(history)) do.call(rbind, history) else
      data.frame(step = integer(0), variable = character(0), score = numeric(0)),
    params   = list(metric = metric, method = method, k = k, tol = tol,
                    seed = seed, n_candidates = p, estimated_fits = est_fits)
  )
}
