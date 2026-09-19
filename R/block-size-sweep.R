# =============================================================================
# Cross-validation error against block size: the inflation of a random-CV
# number as something you see rather than something you argue.
# =============================================================================

#' Cross-validate at a ladder of block sizes
#'
#' A blocked cross-validation with blocks smaller than the autocorrelation
#' range leaks: every held-out point has a near-identical neighbour in the
#' training set, and the score is optimistic in proportion.  The single
#' number a \code{cv_*()} call returns cannot show this.  This runs the same
#' cross-validation at a ladder of block sizes, and by default once with
#' random folds as the fully leaky reference.  It returns the metric at each,
#' with the estimated autocorrelation range alongside so that the curve can be
#' read against it: it rises as the blocks pass the range and then plateaus,
#' and the height of the rise is how much the random-fold number overstated
#' the model.
#'
#' @section The fit budget:
#' Each block size is a full cross-validation, so the cost is
#' \code{length(block_sizes) * k} fits, plus \code{k} for the random
#' reference.  \code{max_fits} caps that (default 60: six sizes at
#' \code{k = 5}, plus the reference).  A sweep that would run past the cap
#' refuses to start, naming the number of fits it would have needed.  Raise
#' \code{max_fits} deliberately; the RF example below takes seconds, a Bayesian
#' \code{fit_fn} takes minutes per fit.
#'
#' @section The ladder:
#' When \code{block_sizes} is \code{NULL}, \code{n_sizes} values are
#' log-spaced from a twenty-fifth to a half of the shorter side of the
#' data's extent, and any size at which the grid would hold fewer than
#' \code{k} blocks is dropped, so every point on the curve is a \code{k}-fold
#' cross-validation of the same shape.  Sizes are in the units of the CRS the
#' folds are built in (\code{make_folds()}'s \code{params$crs}, metres for
#' geographic input), and the returned table records that CRS.
#'
#' @param data_sf An sf object with the response and predictors.
#' @param response_var,predictor_vars Column names.
#' @param fit_fn A \code{function(train_sf)} returning a \code{spatial_fit},
#'   as for \code{\link{cv_spatial}()}; see the example for wrapping a
#'   built-in backend.
#' @param block_sizes Optional numeric vector of block edge lengths to sweep,
#'   in the CRS units the folds are built in.  Default \code{NULL}: the ladder
#'   described above.
#' @param n_sizes Number of sizes in the default ladder.  Default 6.
#' @param k Folds per cross-validation.  Default 5.
#' @param metric Which column of \code{overall} to read.  Default
#'   \code{"RMSE"}.
#' @param include_random Also cross-validate with random folds, as the leaky
#'   reference.  Default \code{TRUE}.
#' @param max_fits The fit budget; see above.  Default 60.
#' @param sac Optional \code{sac_range} from \code{\link{estimate_sac_range}()}
#'   to mark on the curve.  Default \code{NULL}: estimated here from the
#'   response, detrended on \code{predictor_vars}, when \pkg{gstat} is
#'   installed.
#' @param seed Seed for the fold construction at every size.
#' @param quiet Suppress the progress messages.  Default \code{FALSE}.
#' @param ... Passed to \code{\link{cv_spatial}()} at every size
#'   (\code{predict_args}, \code{p}, \code{parallel}, \code{metrics}, ...).
#'   \code{block_size}, \code{folds}, \code{k}, \code{seed} and
#'   \code{auto_range} are set here and cannot be passed.
#' @return A data.frame of class \code{"block_size_sweep"} with one row per
#'   cross-validation: \code{block_size} (\code{NA} for the random
#'   reference), \code{method}, \code{blocks_used}, \code{k} (the folds
#'   actually built), \code{n_folds_succeeded}, \code{value} (the pooled
#'   metric), \code{fold_min}, \code{fold_max} and \code{fold_sd} (its spread
#'   across folds).  Attributes: \code{metric}, \code{sac_range} (the
#'   effective range, or \code{NA}), \code{crs}, \code{n_fits}, and
#'   \code{results}, the full \code{cv_spatial()} result at every size.
#'   \code{plot()} draws it.
#' @family cross-validation
#' @seealso \code{\link{plot.block_size_sweep}()}.
#' @examples
#' if (requireNamespace("ranger", quietly = TRUE) &&
#'     requireNamespace("gstat", quietly = TRUE)) {
#'   library(sf)
#'   set.seed(1)
#'   n <- 200
#'   x <- runif(n, 0, 1000); y <- runif(n, 0, 1000)
#'   d <- as.matrix(dist(cbind(x, y)))
#'   field <- as.numeric(t(chol(exp(-d / 100) + diag(1e-6, n))) %*% rnorm(n))
#'   dat <- st_as_sf(data.frame(x = x, y = y, a = rnorm(n)), coords = c("x", "y"),
#'                   crs = 32632)
#'   dat$z <- field + 0.5 * dat$a + rnorm(n, 0, 0.2)
#'   rf_fn <- function(train_sf)
#'     fit_rf_model(train_sf, "z", "a", include_coords = TRUE, num_trees = 100, seed = 1)
#'   sw <- cv_block_size_sweep(dat, "z", "a", fit_fn = rf_fn, k = 4, n_sizes = 4,
#'                             quiet = TRUE)
#'   sw
#'   if (requireNamespace("ggplot2", quietly = TRUE)) plot(sw)
#' }
#' @export
cv_block_size_sweep <- function(data_sf, response_var, predictor_vars, fit_fn,
                                block_sizes = NULL, n_sizes = 6L, k = 5L,
                                metric = "RMSE", include_random = TRUE,
                                max_fits = 60L, sac = NULL, seed = 123L,
                                quiet = FALSE, ...) {
  .msg <- function(...) if (!quiet) message(...)
  if (!inherits(data_sf, "sf"))
    stop("cv_block_size_sweep(): `data_sf` must be an sf object.", call. = FALSE)
  if (!is.function(fit_fn))
    stop("cv_block_size_sweep(): `fit_fn` must be a function(train_sf) returning a ",
         "spatial_fit, as for cv_spatial().", call. = FALSE)
  if (!is.character(metric) || length(metric) != 1L)
    stop("cv_block_size_sweep(): `metric` must be a single column name of `overall`.",
         call. = FALSE)
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k < 2 || k != round(k))
    stop("cv_block_size_sweep(): `k` must be a single whole number >= 2.", call. = FALSE)
  k <- as.integer(k)
  dots <- list(...)
  fixed <- intersect(names(dots), c("block_size", "folds", "k", "seed", "auto_range"))
  if (length(fixed))
    stop("cv_block_size_sweep(): ", paste(sQuote(fixed), collapse = ", "),
         " cannot be passed through `...`; the sweep sets them.", call. = FALSE)

  # The extent the ladder is built over, in the CRS the folds will use.
  pts <- prep_model_data(data_sf, response_var, predictor_vars)
  bb  <- sf::st_bbox(pts)
  w <- as.numeric(bb["xmax"] - bb["xmin"]); h <- as.numeric(bb["ymax"] - bb["ymin"])
  side <- min(w, h)
  if (!is.finite(side) || side <= 0)
    stop("cv_block_size_sweep(): the data have no extent to build blocks over.",
         call. = FALSE)

  if (is.null(block_sizes)) {
    if (!is.numeric(n_sizes) || length(n_sizes) != 1L || n_sizes < 1)
      stop("cv_block_size_sweep(): `n_sizes` must be a positive number.", call. = FALSE)
    # The top sits a hair under half the side so that floor(side / bs) is 2
    # rather than a floating-point 1, which would drop the largest size.
    block_sizes <- exp(seq(log(side / 25), log(side / 2 * (1 - 1e-9)),
                           length.out = as.integer(n_sizes)))
  } else {
    if (!is.numeric(block_sizes) || !length(block_sizes) || any(!is.finite(block_sizes)) ||
        any(block_sizes <= 0))
      stop("cv_block_size_sweep(): `block_sizes` must be positive numbers.", call. = FALSE)
    block_sizes <- sort(unique(as.numeric(block_sizes)))
  }
  # Every point on the curve is a k-fold CV of the same shape: drop the sizes
  # whose grid holds fewer than k blocks (make_folds() would lower k there).
  n_blocks <- vapply(block_sizes, function(bs) {
    d <- .block_dims_from_size(bb, bs); as.numeric(d$nx) * as.numeric(d$ny)
  }, numeric(1))
  dropped <- block_sizes[n_blocks < k]
  block_sizes <- block_sizes[n_blocks >= k]
  if (length(dropped))
    .log_info("cv_block_size_sweep(): dropping %d block size(s) whose grid holds fewer than k = %d blocks: %s.",
              length(dropped), k, paste(signif(dropped, 3), collapse = ", "))
  if (!length(block_sizes))
    stop("cv_block_size_sweep(): no block size leaves at least k = ", k,
         " blocks over the extent (", signif(w, 3), " x ", signif(h, 3),
         "); pass smaller `block_sizes` or a smaller `k`.", call. = FALSE)

  n_fits <- length(block_sizes) * k + if (isTRUE(include_random)) k else 0L
  if (n_fits > max_fits)
    stop(sprintf(paste0("cv_block_size_sweep(): %d block sizes x %d folds%s = %d ",
                        "model fits, above the budget max_fits = %d. Raise ",
                        "`max_fits` if that cost is intended, or lower ",
                        "`n_sizes` or `k`."),
                 length(block_sizes), k,
                 if (isTRUE(include_random)) sprintf(" + %d for the random reference", k) else "",
                 n_fits, as.integer(max_fits)), call. = FALSE)

  # The autocorrelation range to mark: supplied, or estimated once here.
  sac_range <- NA_real_
  if (!is.null(sac)) {
    sac_range <- suppressWarnings(as.numeric(sac))
  } else if (requireNamespace("gstat", quietly = TRUE)) {
    est <- try(logger::with_log_threshold(
      estimate_sac_range(pts, response_var, predictor_vars = predictor_vars,
                         seed = seed),
      threshold = logger::FATAL, namespace = "spatialkit", index = 2), silent = TRUE)
    if (!inherits(est, "try-error")) sac_range <- suppressWarnings(as.numeric(est))
  }

  # The folds are built here, on the layer as passed (as compare_models_cv()
  # does), so that make_folds()'s own record of the grid -- blocks_used, the
  # k actually built -- is to hand; cv_spatial() then runs on them.  No
  # response is handed to make_folds(): its leakage warning would fire at
  # every size below the range, and this sweep is that diagnostic drawn out.
  if (!("..row_id" %in% names(data_sf))) data_sf$..row_id <- seq_len(nrow(data_sf))
  run_one <- function(bs, method) {
    f <- if (identical(method, "random_kfold"))
      make_folds(data_sf, k = k, method = "random_kfold", seed = seed)
    else
      make_folds(data_sf, k = k, method = "block_kfold", block_size = bs, seed = seed)
    args <- c(list(data_sf = data_sf, response_var = response_var,
                   predictor_vars = predictor_vars, fit_fn = fit_fn,
                   folds = f, k = k, seed = seed), dots)
    cv <- do.call(cv_spatial, args)
    ov <- as.data.frame(cv$overall)
    if (!(metric %in% names(ov)))
      stop(sprintf("cv_block_size_sweep(): '%s' is not a column of cv_spatial()$overall (available: %s).",
                   metric, paste(names(ov), collapse = ", ")), call. = FALSE)
    fm <- as.data.frame(cv$fold_metrics)
    fv <- if (metric %in% names(fm)) suppressWarnings(as.numeric(fm[[metric]])) else numeric(0)
    fv <- fv[is.finite(fv)]
    list(cv = cv, row = data.frame(
      block_size = bs, method = method,
      blocks_used = if (!is.null(f$params$blocks_used)) as.integer(f$params$blocks_used) else NA_integer_,
      k = as.integer(f$k),
      n_folds_succeeded = as.integer(cv$n_folds_succeeded),
      value = suppressWarnings(as.numeric(ov[[metric]][1L])),
      fold_min = if (length(fv)) min(fv) else NA_real_,
      fold_max = if (length(fv)) max(fv) else NA_real_,
      fold_sd  = if (length(fv) > 1L) stats::sd(fv) else NA_real_,
      stringsAsFactors = FALSE))
  }

  rows <- list(); results <- list()
  if (isTRUE(include_random)) {
    .msg("cv_block_size_sweep(): random folds (the leaky reference) ...")
    one <- suppressMessages(run_one(NA_real_, "random_kfold"))
    results[["random_kfold"]] <- one$cv
    rows[[length(rows) + 1L]] <- one$row
  }
  for (bs in block_sizes) {
    .msg(sprintf("cv_block_size_sweep(): block size %s ...", format(signif(bs, 4))))
    one <- suppressMessages(run_one(bs, "block_kfold"))
    results[[as.character(signif(bs, 6))]] <- one$cv
    rows[[length(rows) + 1L]] <- one$row
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  structure(out,
            metric = metric, sac_range = sac_range,
            crs = .fold_crs_label(pts), k = k, n_fits = n_fits,
            response_var = response_var, results = results,
            class = c("block_size_sweep", "data.frame"))
}


#' @export
print.block_size_sweep <- function(x, ...) {
  metric <- attr(x, "metric"); r <- attr(x, "sac_range")
  # `[` keeps the class.  A column subset such as `x[, 1:3]` loses the
  # attributes this header is built from and the columns the table below
  # selects, so it used to fail on `is.finite(NULL)`.  Print what is left as a
  # plain table instead; knitr calls print() on a data frame unasked, so this
  # path is reachable from a document.  (A row subset keeps everything and
  # prints normally.)
  needed <- c("block_size", "method", "blocks_used", "k", "n_folds_succeeded",
              "value", "fold_min", "fold_max")
  if (is.null(metric) || !all(needed %in% names(x))) {
    cat("Block-size sweep (subset; the run summary is not carried by a subset)\n\n")
    print(as.data.frame(unclass(x)), row.names = FALSE)
    return(invisible(x))
  }
  cat(sprintf("Cross-validation %s against block size (%d fits, k = %d, sizes in %s)\n",
              metric, attr(x, "n_fits"), attr(x, "k"), attr(x, "crs") %||% "CRS units"))
  if (is.finite(r)) cat(sprintf("  estimated autocorrelation range: %.1f\n", r))
  else cat("  autocorrelation range: not identified\n")
  df <- as.data.frame(x)
  df$block_size <- ifelse(is.na(df$block_size), "random", format(signif(df$block_size, 4)))
  print(df[, c("block_size", "method", "blocks_used", "k", "n_folds_succeeded",
               "value", "fold_min", "fold_max")], row.names = FALSE, digits = 4)
  invisible(x)
}


#' Plot cross-validation error against block size
#'
#' The curve from \code{\link{cv_block_size_sweep}()}: the pooled metric at
#' each block size, the fold-to-fold range as a band, the random-fold
#' reference as a dashed line, and the estimated autocorrelation range as a
#' vertical marker.  Blocks smaller than the range leak, so the curve rises
#' from the reference towards the range and plateaus beyond it; the height
#' of the rise is what the random-fold number overstated.
#'
#' @param x A \code{block_size_sweep}.
#' @param ... Ignored.
#' @return A \code{ggplot} object.
#' @family plotting
#' @examples
#' if (requireNamespace("ranger", quietly = TRUE) &&
#'     requireNamespace("gstat", quietly = TRUE) &&
#'     requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(sf)
#'   set.seed(1)
#'   n <- 200
#'   x <- runif(n, 0, 1000); y <- runif(n, 0, 1000)
#'   d <- as.matrix(dist(cbind(x, y)))
#'   field <- as.numeric(t(chol(exp(-d / 100) + diag(1e-6, n))) %*% rnorm(n))
#'   dat <- st_as_sf(data.frame(x = x, y = y, a = rnorm(n)), coords = c("x", "y"),
#'                   crs = 32632)
#'   dat$z <- field + 0.5 * dat$a + rnorm(n, 0, 0.2)
#'   rf_fn <- function(train_sf)
#'     fit_rf_model(train_sf, "z", "a", include_coords = TRUE, num_trees = 100,
#'                  seed = 1)
#'   sw <- cv_block_size_sweep(dat, "z", "a", fit_fn = rf_fn, k = 4, n_sizes = 4,
#'                             quiet = TRUE)
#'   # Error rises from the random-fold reference (dashed) towards the estimated
#'   # autocorrelation range (vertical marker) and plateaus past it.  The height
#'   # of that rise is what random folds were hiding.
#'   plot(sw)
#' }
#' @export
plot.block_size_sweep <- function(x, ...) {
  .need_ggplot("plot.block_size_sweep()")
  df <- as.data.frame(x)
  metric <- attr(x, "metric"); r <- attr(x, "sac_range")
  blk <- df[df$method == "block_kfold" & is.finite(df$value), , drop = FALSE]
  if (!nrow(blk))
    stop("plot.block_size_sweep(): no block size produced a finite metric.", call. = FALSE)
  rnd <- df[df$method == "random_kfold" & is.finite(df$value), , drop = FALSE]
  unit <- attr(x, "crs") %||% "CRS units"

  p <- ggplot2::ggplot(blk, ggplot2::aes(x = .data$block_size, y = .data$value))
  if (any(is.finite(blk$fold_min) & is.finite(blk$fold_max)))
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$fold_min, ymax = .data$fold_max),
                                  fill = "#2166AC", alpha = 0.12)
  p <- p + ggplot2::geom_line(colour = "#2166AC") +
    ggplot2::geom_point(colour = "#2166AC", size = 2.5)
  if (nrow(rnd))
    p <- p + ggplot2::geom_hline(yintercept = rnd$value[1L], linetype = "dashed",
                                 colour = "grey40")
  if (is.finite(r))
    p <- p + ggplot2::geom_vline(xintercept = r, linetype = "dotted", colour = "#B2182B")
  minimise <- !(metric %in% c("R2", "Adj_R2", "coverage_50", "coverage_80", "coverage_95"))
  p + ggplot2::scale_x_log10() +
    ggplot2::labs(
      title = sprintf("Cross-validated %s against block size", metric),
      subtitle = if (is.finite(r))
        sprintf("Estimated autocorrelation range %.0f (dotted); blocks below it leak", r)
      else "Autocorrelation range not identified; no marker drawn",
      caption = paste(c(
        if (nrow(rnd)) sprintf("Dashed: random folds, %s = %.3g (the leaky reference)", metric, rnd$value[1L]),
        sprintf("Band: fold-to-fold range; %d folds per size, %d fits in total; %s is better",
                attr(x, "k"), attr(x, "n_fits"), if (minimise) "lower" else "higher")),
        collapse = "\n"),
      x = sprintf("Block edge length (%s, log scale)", unit), y = metric) +
    ggplot2::theme_minimal()
}
