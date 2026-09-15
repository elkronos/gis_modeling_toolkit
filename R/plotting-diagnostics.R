# =============================================================================
# Diagnostic plots: the curve behind a chosen point, the folds behind a pooled
# number, the distribution behind a count.  Every function here draws data the
# package already returns; nothing is recomputed.
# =============================================================================

# -----------------------------------------------------------------------------
# 9.1  Per-fold metric strip
# -----------------------------------------------------------------------------

#' Plot one cross-validation metric fold by fold
#'
#' A pooled RMSE of 3.2 can come from 3.2 in every fold or from 1.1 in eight
#' folds and 14 in one --- a model that works, and a model that fails in one
#' region --- and the pooled number cannot tell the two apart.  This draws the
#' metric of each fold as a point, sized by the number of held-out
#' predictions the fold contributed, with the pooled value from
#' \code{overall} as a horizontal line, so the spread behind the number is
#' visible.  A \code{\link{compare_models_cv}()} result draws one panel per
#' model on a shared scale, folds aligned, which is the comparison the shared
#' fold set was built for.
#'
#' Any column of \code{fold_metrics} can be drawn, including a backend's
#' extras (\code{bandwidth}, \code{CRPS}, \code{coverage_95}) and columns a
#' \code{metrics} function added.  A column that is \code{NA} in every fold is
#' refused with a message saying why rather than drawn as an empty panel:
#' \code{Adj_R2} is \code{NA} for every backend unless \code{p} was passed to
#' \code{cv_spatial()}, by design.  The pooled line is drawn from
#' \code{overall} when it carries the metric, from
#' \code{predictive_coverage} for \code{cv_bayes()}'s coverage and CRPS
#' columns, and not at all for a per-fold extra that has no pooled
#' counterpart (a bandwidth), in which case the caption says so.
#'
#' @param cv The list returned by \code{\link{cv_spatial}()},
#'   \code{\link{cv_gwr}()}, \code{\link{cv_bayes}()}, \code{\link{cv_rf}()}
#'   or \code{\link{compare_models_cv}()}.
#' @param metric Character(1) naming a column of \code{fold_metrics} (or of
#'   \code{by_fold}).  Default \code{"RMSE"}.
#' @param ... Ignored.
#' @return A \code{ggplot} object.
#' @family plotting
#' @examples
#' if (requireNamespace("ranger", quietly = TRUE) &&
#'     requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(sf)
#'   set.seed(1)
#'   n <- 150
#'   dat <- st_as_sf(
#'     data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), a = rnorm(n)),
#'     coords = c("x", "y"), crs = 32632
#'   )
#'   dat$z <- 2 * dat$a + 0.003 * st_coordinates(dat)[, 1] + rnorm(n, 0, 0.5)
#'   cv <- cv_rf(dat, "z", "a", k = 5, num_trees = 100)
#'   plot_cv_metrics(cv, "RMSE")
#' }
#' @export
plot_cv_metrics <- function(cv, metric = "RMSE", ...) {
  .need_ggplot("plot_cv_metrics()")
  if (!is.character(metric) || length(metric) != 1L || !nzchar(metric))
    stop("plot_cv_metrics(): `metric` must be a single column name.", call. = FALSE)

  # One shape for the two inputs: a per-fold frame with a `model` column, and
  # a one-row-per-model overall frame.
  if (is.list(cv) && !is.null(cv$by_fold) && !is.null(cv$overall) &&
      "model" %in% names(cv$overall)) {
    fm <- as.data.frame(cv$by_fold)
    ov <- as.data.frame(cv$overall)
    pooled <- .pooled_metric_by_model(cv, ov, metric)
    what <- "compare_models_cv()"
  } else if (is.list(cv) && !is.null(cv$fold_metrics) && !is.null(cv$overall)) {
    fm <- as.data.frame(cv$fold_metrics)
    fm$model <- rep("model", nrow(fm))
    ov <- as.data.frame(cv$overall)
    pooled <- data.frame(model = "model",
                         value = .pooled_metric_one(cv, ov, metric),
                         stringsAsFactors = FALSE)
    what <- "cv_*()"
  } else {
    stop("plot_cv_metrics(): `cv` must be the list returned by cv_spatial(), ",
         "cv_gwr(), cv_bayes(), cv_rf() or compare_models_cv().", call. = FALSE)
  }

  if (!nrow(fm))
    stop("plot_cv_metrics(): the result has no per-fold metrics; every fold ",
         "failed, so there is nothing to draw (see n_folds_succeeded).",
         call. = FALSE)
  if (!(metric %in% names(fm)))
    stop(sprintf(paste0("plot_cv_metrics(): '%s' is not a column of the per-fold ",
                        "metrics. Available: %s."),
                 metric, paste(setdiff(names(fm), c("fold", "model")), collapse = ", ")),
         call. = FALSE)
  v <- suppressWarnings(as.numeric(fm[[metric]]))
  if (!any(is.finite(v))) {
    why <- if (identical(metric, "Adj_R2"))
      " It is NA for every backend unless `p` was passed to cv_spatial(), by design."
    else if (grepl("^coverage_|^CRPS$", metric))
      " The posterior predictive draws that give it may have failed in every fold (compute_pred_intervals = FALSE gives the same)."
    else ""
    stop(sprintf("plot_cv_metrics(): '%s' is NA in every fold, so there is nothing to draw.%s",
                 metric, why), call. = FALSE)
  }

  df <- data.frame(model = as.character(fm$model), fold = fm$fold, value = v,
                   n_pred = if ("n_pred" %in% names(fm)) as.numeric(fm$n_pred) else NA_real_,
                   stringsAsFactors = FALSE)
  df <- df[is.finite(df$value), , drop = FALSE]
  df$fold <- factor(df$fold, levels = sort(unique(df$fold)))
  df$model <- factor(df$model, levels = unique(df$model))
  pooled$model <- factor(pooled$model, levels = levels(df$model))
  pooled <- pooled[is.finite(pooled$value), , drop = FALSE]

  n_models <- nlevels(df$model)
  title <- sprintf("%s by fold", metric)
  caption <- if (nrow(pooled))
    "Dashed line: the pooled value from `overall`"
  else
    sprintf("No pooled value: `%s` is a per-fold quantity with no counterpart in `overall`", metric)

  size_aes <- if (any(is.finite(df$n_pred))) ggplot2::aes(size = .data$n_pred) else NULL
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$fold, y = .data$value))
  if (nrow(pooled))
    p <- p + ggplot2::geom_hline(data = pooled, ggplot2::aes(yintercept = .data$value),
                                 linetype = "dashed", colour = "#B2182B")
  p <- p + ggplot2::geom_point(mapping = size_aes, alpha = 0.8) +
    ggplot2::scale_size_continuous(name = "Held-out rows", range = c(1.5, 5)) +
    ggplot2::labs(title = title, x = "Fold", y = metric, caption = caption) +
    ggplot2::theme_minimal()
  if (n_models > 1L)
    p <- p + ggplot2::facet_wrap(~ model, ncol = 1L, scales = "fixed")
  p
}

#' The pooled value of a metric for a single cv_*() result
#' @keywords internal
#' @noRd
.pooled_metric_one <- function(cv, ov, metric) {
  if (nrow(ov) >= 1L && metric %in% names(ov)) {
    val <- suppressWarnings(as.numeric(ov[[metric]][1L]))
    if (is.finite(val)) return(val)
  }
  pc <- cv$predictive_coverage
  if (is.list(pc)) {
    key <- if (identical(metric, "CRPS")) "mean_CRPS" else metric
    if (!is.null(pc[[key]])) return(suppressWarnings(as.numeric(pc[[key]])))
  }
  NA_real_
}

#' The pooled values of a metric per model for a compare_models_cv() result
#' @keywords internal
#' @noRd
.pooled_metric_by_model <- function(cmp, ov, metric) {
  models <- as.character(ov$model)
  vals <- vapply(seq_along(models), function(i) {
    if (metric %in% names(ov)) {
      v <- suppressWarnings(as.numeric(ov[[metric]][i]))
      if (is.finite(v)) return(v)
    }
    # The Bayesian row's coverage and CRPS are in `overall` since 3.3; an
    # older result, or a per-fold extra, may still be reachable through the
    # backend's own list.
    sub <- switch(models[i], GWR = cmp$gwr_cv, Bayesian = cmp$bayes_cv,
                  RF = cmp$rf_cv, NULL)
    if (is.list(sub) && !is.null(sub$overall))
      return(.pooled_metric_one(sub, as.data.frame(sub$overall), metric))
    NA_real_
  }, numeric(1))
  data.frame(model = models, value = vals, stringsAsFactors = FALSE)
}


# -----------------------------------------------------------------------------
# 9.2  Area-of-applicability dissimilarity distribution
# -----------------------------------------------------------------------------

#' Plot the dissimilarity distribution behind an area of applicability
#'
#' \code{n_outside} says how many prediction locations fall outside the area
#' of applicability; it does not say whether the rest sit comfortably inside
#' or crowd against the threshold, nor how far outside the outsiders are.
#' This draws the dissimilarity index of the prediction locations against
#' that of the cross-validated training data, with the threshold marked, so
#' the prediction set can be read as mostly inside, marginal or largely
#' outside.  The training curve is the reference the threshold was derived
#' from: its upper tail ends where the threshold is (or below it, when the
#' outlier fence removed the tail).
#'
#' @param x An \code{aoa} object from \code{\link{area_of_applicability}()}.
#' @param type \code{"ecdf"} (default), the two empirical distribution
#'   functions on one axis, or \code{"histogram"}, the prediction DI as bars
#'   with the training DI as an outline.
#' @param ... Ignored.
#' @return A \code{ggplot} object.
#' @family plotting
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(sf)
#'   set.seed(2)
#'   n <- 200
#'   train <- st_as_sf(
#'     data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000),
#'                a = rnorm(n), b = rnorm(n)),
#'     coords = c("x", "y"), crs = 32632)
#'   train$z <- train$a - train$b + rnorm(n, 0, 0.3)
#'   # Prediction locations whose predictor `a` drifts beyond the training range.
#'   new <- st_as_sf(
#'     data.frame(x = runif(100, 0, 1000), y = runif(100, 0, 1000),
#'                a = rnorm(100, mean = 2), b = rnorm(100)),
#'     coords = c("x", "y"), crs = 32632)
#'   aoa <- area_of_applicability(new, train_sf = train, predictor_vars = c("a", "b"))
#'   plot(aoa)
#'   plot(aoa, type = "histogram")
#' }
#' @export
plot.aoa <- function(x, type = c("ecdf", "histogram"), ...) {
  .need_ggplot("plot.aoa()")
  type <- match.arg(type)
  if (!inherits(x, "aoa") || is.null(x$aoa) || !("DI" %in% names(x$aoa)))
    stop("plot.aoa(): `x` must be the object returned by area_of_applicability().",
         call. = FALSE)
  di_new <- suppressWarnings(as.numeric(sf::st_drop_geometry(x$aoa)$DI))
  di_tr  <- suppressWarnings(as.numeric(x$train_DI))
  thr    <- as.numeric(x$threshold)
  di_new <- di_new[is.finite(di_new)]
  di_tr  <- di_tr[is.finite(di_tr)]
  if (!length(di_new))
    stop("plot.aoa(): no finite dissimilarity index at any prediction location ",
         "(every row had a missing or non-finite predictor).", call. = FALSE)

  df <- rbind(
    data.frame(set = "Prediction locations", DI = di_new, stringsAsFactors = FALSE),
    if (length(di_tr)) data.frame(set = "Training (cross-validated)", DI = di_tr,
                                  stringsAsFactors = FALSE)
  )
  df$set <- factor(df$set, levels = c("Training (cross-validated)", "Prediction locations"))

  n_all <- x$n_new %||% length(di_new)
  share_out <- if (length(di_new)) mean(di_new > thr) else NA_real_
  # Where the prediction set sits relative to the threshold: the fraction
  # inside, and how close the inside ones run to the edge.
  q_in <- if (any(di_new <= thr)) stats::quantile(di_new[di_new <= thr], 0.9) / thr else NA_real_
  subtitle <- sprintf(
    "%d of %d prediction locations outside (DI > %.3g)%s",
    x$n_outside %||% sum(di_new > thr), n_all, thr,
    if (is.finite(q_in))
      sprintf("; 90%% of those inside sit below %.0f%% of the threshold", 100 * q_in)
    else "")
  caption <- sprintf("Threshold %s: %s",
                     if (isTRUE(x$params$threshold_supplied)) "supplied" else "from the training DI",
                     if (isTRUE(x$params$folds_supplied))
                       sprintf("training DI is cross-validated over %s folds",
                               .aoa_folds_label(x$params$folds_method))
                     else "training DI is not cross-validated (no folds), so the threshold is optimistic")

  if (type == "ecdf") {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$DI, colour = .data$set)) +
      ggplot2::stat_ecdf(geom = "step", linewidth = 0.8) +
      ggplot2::geom_vline(xintercept = thr, linetype = "dashed", colour = "#B2182B") +
      ggplot2::scale_colour_manual(values = c("Training (cross-validated)" = "grey45",
                                              "Prediction locations" = "#2166AC"),
                                   name = NULL) +
      ggplot2::labs(title = "Dissimilarity index: prediction locations against training",
                    subtitle = subtitle, caption = caption,
                    x = "Dissimilarity index (DI)", y = "Cumulative share") +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom")
    return(p)
  }

  new_df <- df[df$set == "Prediction locations", , drop = FALSE]
  tr_df  <- df[df$set != "Prediction locations", , drop = FALSE]
  p <- ggplot2::ggplot() +
    ggplot2::geom_histogram(data = new_df,
                            ggplot2::aes(x = .data$DI, y = ggplot2::after_stat(.data$density)),
                            bins = 30, fill = "#2166AC", alpha = 0.6)
  if (nrow(tr_df))
    p <- p + ggplot2::geom_freqpoly(data = tr_df,
                                    ggplot2::aes(x = .data$DI, y = ggplot2::after_stat(.data$density)),
                                    bins = 30, colour = "grey30", linewidth = 0.7)
  p + ggplot2::geom_vline(xintercept = thr, linetype = "dashed", colour = "#B2182B") +
    ggplot2::labs(title = "Dissimilarity index: prediction locations (bars) against training (line)",
                  subtitle = subtitle, caption = caption,
                  x = "Dissimilarity index (DI)", y = "Density") +
    ggplot2::theme_minimal()
}


# -----------------------------------------------------------------------------
# 9.4  Interval calibration of cv_bayes()
# -----------------------------------------------------------------------------

#' Plot the interval calibration of a Bayesian cross-validation
#'
#' \code{\link{cv_bayes}()} scores every fold's posterior predictive intervals
#' at each of \code{coverage_levels}: the share of held-out observations the
#' interval of that nominal level contained.  One coverage number cannot show
#' a pattern; the pairs can.  This draws observed coverage against nominal
#' with the diagonal, one point per level for the fold-weighted pooled value
#' (from \code{predictive_coverage}) and one faint point per fold, so
#' systematic over-confidence (points below the line) or intervals wider than
#' they need to be (above it) are read at a glance.  Three levels is a thin
#' curve; pass \code{coverage_levels = seq(0.1, 0.9, by = 0.1)} to
#' \code{cv_bayes()} for a full one --- the levels are read off the column
#' names, so whatever was computed is drawn.
#'
#' @param cv The list returned by \code{\link{cv_bayes}()}, or a
#'   \code{\link{compare_models_cv}()} result that ran the Bayesian backend
#'   (its \code{$bayes_cv} is used).
#' @param ... Ignored.
#' @return A \code{ggplot} object.
#' @family plotting
#' @examples
#' \donttest{
#' if (requireNamespace("brms", quietly = TRUE) &&
#'     requireNamespace("ggplot2", quietly = TRUE)) {
#'   # cv <- cv_bayes(dat, "z", "a", k = 3, coverage_levels = seq(0.1, 0.9, 0.2))
#'   # plot_calibration(cv)
#' }
#' }
#' @export
plot_calibration <- function(cv, ...) {
  .need_ggplot("plot_calibration()")
  if (is.list(cv) && !is.null(cv$bayes_cv)) cv <- cv$bayes_cv
  if (!is.list(cv) || is.null(cv$fold_metrics))
    stop("plot_calibration(): `cv` must be the list returned by cv_bayes(), or a ",
         "compare_models_cv() result whose Bayesian backend ran.", call. = FALSE)
  fm <- as.data.frame(cv$fold_metrics)
  cov_cols <- grep("^coverage_[0-9]+$", names(fm), value = TRUE)
  if (!length(cov_cols))
    stop("plot_calibration(): the result carries no coverage_* columns; ",
         "cv_bayes() computes them when compute_pred_intervals = TRUE and at ",
         "least one fold produced posterior predictive draws.", call. = FALSE)
  nominal <- as.numeric(sub("^coverage_", "", cov_cols)) / 100

  per_fold <- do.call(rbind, lapply(seq_along(cov_cols), function(j) {
    data.frame(fold = fm$fold, nominal = nominal[j],
               observed = suppressWarnings(as.numeric(fm[[cov_cols[j]]])),
               n_pred = if ("n_pred" %in% names(fm)) as.numeric(fm$n_pred) else NA_real_,
               stringsAsFactors = FALSE)
  }))
  per_fold <- per_fold[is.finite(per_fold$observed), , drop = FALSE]
  if (!nrow(per_fold))
    stop("plot_calibration(): coverage is NA in every fold (the posterior ",
         "predictive draws failed everywhere), so there is nothing to draw.",
         call. = FALSE)

  pc <- cv$predictive_coverage
  pooled <- data.frame(nominal = nominal, observed = vapply(cov_cols, function(cn) {
    if (is.list(pc) && !is.null(pc[[cn]])) suppressWarnings(as.numeric(pc[[cn]]))
    else {
      # Fold-weighted mean, the same pooling cv_bayes() does.
      sub <- per_fold[per_fold$nominal == nominal[match(cn, cov_cols)], , drop = FALSE]
      w <- if (all(is.finite(sub$n_pred))) sub$n_pred else rep(1, nrow(sub))
      if (nrow(sub)) sum(sub$observed * w) / sum(w) else NA_real_
    }
  }, numeric(1)))
  pooled <- pooled[is.finite(pooled$observed), , drop = FALSE]

  gap <- pooled$observed - pooled$nominal
  verdict <- if (!length(gap)) ""
    else if (all(gap < -0.05)) "intervals are too narrow at every level (over-confident)"
    else if (all(gap > 0.05)) "intervals are wider than they need to be at every level"
    else if (max(abs(gap)) <= 0.05) "coverage is within 5 points of nominal at every level"
    else "coverage departs from nominal by more than 5 points at some level"
  n_levels <- length(nominal)

  p <- ggplot2::ggplot() +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_point(data = per_fold,
                        ggplot2::aes(x = .data$nominal, y = .data$observed),
                        colour = "grey60", alpha = 0.5, size = 1.8) +
    ggplot2::geom_line(data = pooled, ggplot2::aes(x = .data$nominal, y = .data$observed),
                       colour = "#2166AC") +
    ggplot2::geom_point(data = pooled, ggplot2::aes(x = .data$nominal, y = .data$observed),
                        colour = "#2166AC", size = 3) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(title = "Interval calibration of the posterior predictive intervals",
                  subtitle = verdict,
                  caption = sprintf(paste0("Blue: pooled across folds (weighted by held-out rows); ",
                                           "grey: individual folds. %d nominal level%s%s."),
                                    n_levels, if (n_levels == 1L) "" else "s",
                                    if (n_levels < 5L)
                                      "; pass coverage_levels = seq(0.1, 0.9, by = 0.1) to cv_bayes() for a full curve"
                                    else ""),
                  x = "Nominal coverage", y = "Observed coverage") +
    ggplot2::theme_minimal()
  p
}


# -----------------------------------------------------------------------------
# 9.5  One sweep plot, and the callers that can feed it
# -----------------------------------------------------------------------------

#' Draw a criterion against a swept parameter, with the chosen point marked
#'
#' The one drawing routine behind \code{plot.resolution_profile()},
#' \code{plot.feature_selection()} and \code{plot.gwr_model_selection()}.  A
#' sharp optimum means the data chose the point; a flat curve means the
#' rule did, and the chosen value then deserves less weight than it reads
#' with --- which is what the picture is for.
#'
#' @param df Data frame with columns \code{x}, \code{y}, \code{panel}
#'   (facet; one level for a single panel) and optionally \code{label}
#'   (text at the point) and \code{role} (\code{"candidate"} points are drawn
#'   faint, everything else full).
#' @param chosen Data frame with \code{panel}, \code{x}, \code{y}: the point
#'   the rule picked in each panel (may have zero rows).
#' @param flat Optional data frame with \code{panel}, \code{xmin},
#'   \code{xmax}: the region over which the criterion is within tolerance of
#'   its optimum, shaded.
#' @param title,subtitle,caption,x_lab,y_lab Labels.
#' @param log_x Log-scale the x axis.
#' @param connect Draw a line through the non-candidate points of each panel.
#' @return A \code{ggplot} object.
#' @keywords internal
#' @noRd
.plot_sweep <- function(df, chosen, flat = NULL, title = NULL, subtitle = NULL,
                        caption = NULL, x_lab = "x", y_lab = "criterion",
                        log_x = FALSE, connect = TRUE) {
  df$panel <- factor(df$panel, levels = unique(df$panel))
  if (is.null(df$role)) df$role <- "value"
  main <- df[df$role != "candidate", , drop = FALSE]
  cand <- df[df$role == "candidate", , drop = FALSE]
  p <- ggplot2::ggplot()
  if (!is.null(flat) && nrow(flat)) {
    flat$panel <- factor(flat$panel, levels = levels(df$panel))
    p <- p + ggplot2::geom_rect(data = flat,
                                ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                                             ymin = -Inf, ymax = Inf),
                                fill = "#2166AC", alpha = 0.10)
  }
  if (nrow(cand))
    p <- p + ggplot2::geom_point(data = cand, ggplot2::aes(x = .data$x, y = .data$y),
                                 colour = "grey65", alpha = 0.6, size = 1.6)
  if (connect && nrow(main))
    p <- p + ggplot2::geom_line(data = main, ggplot2::aes(x = .data$x, y = .data$y),
                                colour = "grey30")
  if (nrow(main))
    p <- p + ggplot2::geom_point(data = main, ggplot2::aes(x = .data$x, y = .data$y),
                                 colour = "grey30", size = 2)
  if (!is.null(main$label) && any(nzchar(main$label)))
    p <- p + ggplot2::geom_text(data = main[nzchar(main$label), , drop = FALSE],
                                ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
                                vjust = -0.8, size = 3, colour = "grey20")
  if (!is.null(chosen) && nrow(chosen)) {
    chosen$panel <- factor(chosen$panel, levels = levels(df$panel))
    p <- p + ggplot2::geom_vline(data = chosen, ggplot2::aes(xintercept = .data$x),
                                 linetype = "dotted", colour = "#B2182B") +
      ggplot2::geom_point(data = chosen, ggplot2::aes(x = .data$x, y = .data$y),
                          colour = "#B2182B", size = 3.5)
  }
  p <- p + ggplot2::labs(title = title, subtitle = subtitle, caption = caption,
                         x = x_lab, y = y_lab) +
    ggplot2::theme_minimal()
  if (log_x) p <- p + ggplot2::scale_x_log10()
  if (nlevels(df$panel) > 1L)
    p <- p + ggplot2::facet_wrap(~ panel, ncol = 1L, scales = "free_y")
  p
}


#' Plot a resolution profile
#'
#' The criteria of a \code{\link{resolution_profile}()} against the number of
#' cells, one panel per criterion on a shared x axis, with the level each
#' criterion selects marked and the region over which it is within
#' \code{tol} of its optimum shaded --- the flat region
#' \code{\link{select_resolution}()} reports, drawn.  A criterion whose
#' optimum sits at the support ceiling or the range floor is captioned as
#' such, because there the bound is choosing, not the criterion.
#'
#' @param x A \code{resolution_profile}.
#' @param criteria Character vector of criteria to draw, any of
#'   \code{"cp"}, \code{"reliability"}, \code{"elbow"}, \code{"moran_z"},
#'   \code{"wss"}.  Default: every one of the four selectable criteria that
#'   is finite at some level.  \code{"wss"} is the raw curve behind
#'   \code{"elbow"} and has no selected point.
#' @param tol Passed to \code{\link{select_resolution}()} for the flat region.
#' @param ... Ignored.
#' @return A \code{ggplot} object.
#' @family plotting
#' @examples
#' if (requireNamespace("gstat", quietly = TRUE) &&
#'     requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(sf)
#'   set.seed(3)
#'   n <- 400
#'   xy <- data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000))
#'   xy$z <- sin(xy$x / 200) + cos(xy$y / 250) + rnorm(n, sd = 0.3)
#'   pts <- st_as_sf(xy, coords = c("x", "y"), crs = 32632)
#'   prof <- resolution_profile(pts, response_var = "z", n_levels = 10)
#'   plot(prof)
#'   plot(prof, criteria = c("cp", "wss"))
#' }
#' @export
plot.resolution_profile <- function(x, criteria = NULL, tol = 0.02, ...) {
  .need_ggplot("plot.resolution_profile()")
  if (!inherits(x, "resolution_profile"))
    stop("plot.resolution_profile(): `x` must come from resolution_profile().",
         call. = FALSE)
  selectable <- c("cp", "reliability", "elbow", "moran_z")
  allowed <- c(selectable, "wss")
  if (is.null(criteria)) {
    criteria <- selectable[vapply(selectable, function(cn)
      any(is.finite(suppressWarnings(as.numeric(x[[cn]])))), logical(1))]
    if (!length(criteria)) criteria <- "wss"
  }
  bad <- setdiff(criteria, allowed)
  if (length(bad))
    stop("plot.resolution_profile(): unknown criteria: ", paste(bad, collapse = ", "),
         ". Choose from ", paste(allowed, collapse = ", "), ".", call. = FALSE)

  lv <- as.numeric(x$levels)
  labels <- c(cp = "Mallows' Cp (lower is better)",
              reliability = "Reliability of cell means (higher is better)",
              elbow = "WSS elbow statistic (higher is better)",
              moran_z = "|Moran's z| of cell means (lower is better)",
              wss = "Within-cluster sum of squares")
  rows <- list(); chosen <- list(); flat <- list(); notes <- character(0)
  for (cn in criteria) {
    v <- suppressWarnings(as.numeric(x[[cn]]))
    if (cn == "moran_z") v <- abs(v)
    ok <- is.finite(v)
    rows[[cn]] <- data.frame(panel = labels[[cn]], x = lv[ok], y = v[ok],
                             stringsAsFactors = FALSE)
    if (cn %in% selectable && any(ok)) {
      sel <- try(select_resolution(x, criterion = cn, tol = tol), silent = TRUE)
      if (!inherits(sel, "try-error")) {
        chosen[[cn]] <- data.frame(panel = labels[[cn]], x = sel$best,
                                   y = v[match(sel$best, lv)], stringsAsFactors = FALSE)
        if (length(sel$flat))
          flat[[cn]] <- data.frame(panel = labels[[cn]], xmin = min(sel$flat),
                                   xmax = max(sel$flat), stringsAsFactors = FALSE)
        if (isTRUE(sel$at_ceiling))
          notes <- c(notes, sprintf("%s: optimum at the support ceiling (the bound is choosing)", cn))
        else if (isTRUE(sel$at_floor))
          notes <- c(notes, sprintf("%s: optimum at the range floor (the bound is choosing)", cn))
      }
    }
  }
  df <- do.call(rbind, rows)
  if (!nrow(df))
    stop("plot.resolution_profile(): none of the requested criteria is finite at ",
         "any level.", call. = FALSE)
  bounds <- attr(x, "bounds")
  subtitle <- if (is.list(bounds))
    sprintf("%d levels from %d to %d cells (floor %s, ceiling %s)%s",
            nrow(x), min(lv), max(lv),
            if (is.finite(bounds$floor %||% NA)) format(bounds$floor) else "none",
            if (is.finite(bounds$ceiling %||% NA)) format(bounds$ceiling) else "none",
            if (isTRUE(bounds$supported)) "" else "; the floor exceeds the ceiling")
  else NULL
  caption <- paste(c("Red: the level each criterion selects; shaded: within tolerance of its optimum",
                     notes), collapse = "\n")
  .plot_sweep(df, chosen = do.call(rbind, chosen), flat = do.call(rbind, flat),
              title = "Resolution profile", subtitle = subtitle, caption = caption,
              x_lab = "Number of cells", y_lab = NULL,
              log_x = length(lv) > 2L && max(lv) / min(lv) > 8)
}


#' Plot the path of a forward feature selection
#'
#' \code{\link{select_features_forward}()} scores every candidate at every
#' step and keeps the best; its \code{history} holds all of them.  This draws
#' the accepted variable's score at each step as the path, every other
#' candidate's score at that step as a faint point, and the step at which the
#' selection stopped in red --- so the picture says whether the last variable
#' was a clear gain or the first that happened to clear \code{tol}, and
#' whether the runner-up would have done as well.  The scores are the
#' selection's own cross-validated criterion, optimistically biased by the
#' selection (see the help page's section on that); when a hold-out score
#' was computed (\code{select_on = "split"}) it is drawn as a separate mark
#' at the final step and named in the caption.
#'
#' @param x The list returned by \code{\link{select_features_forward}()}.
#' @param ... Ignored.
#' @return A \code{ggplot} object.
#' @family plotting
#' @examples
#' if (requireNamespace("ranger", quietly = TRUE) &&
#'     requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(sf)
#'   set.seed(4)
#'   n <- 150
#'   dat <- st_as_sf(
#'     data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000),
#'                a = rnorm(n), b = rnorm(n), c = rnorm(n)),
#'     coords = c("x", "y"), crs = 32632)
#'   dat$z <- 2 * dat$a - dat$b + rnorm(n, 0, 0.5)
#'   fit_fn <- function(train_sf, vars)
#'     fit_rf_model(train_sf, "z", vars, num_trees = 80, seed = 1)
#'   sel <- select_features_forward(dat, "z", c("a", "b", "c"), fit_fn = fit_fn,
#'                                  k = 3, quiet = TRUE)
#'   plot(sel)
#' }
#' @export
plot.feature_selection <- function(x, ...) {
  .need_ggplot("plot.feature_selection()")
  h <- x$history
  if (!is.data.frame(h) || !all(c("step", "variable", "score") %in% names(h)))
    stop("plot.feature_selection(): `x` must be the list returned by ",
         "select_features_forward(), with its `history` frame.", call. = FALSE)
  h <- h[is.finite(suppressWarnings(as.numeric(h$score))), , drop = FALSE]
  if (!nrow(h))
    stop("plot.feature_selection(): no candidate produced a finite score at any ",
         "step, so there is nothing to draw.", call. = FALSE)
  metric   <- x$params$metric %||% "score"
  minimise <- !(metric %in% c("R2", "Adj_R2"))
  selected <- x$selected
  n_sel    <- length(selected)

  # The accepted variable at step s is selected[s]; the path runs through
  # its score.  Steps past n_sel were scored and rejected (nothing cleared
  # tol), and their best candidate is drawn as part of the path too, so the
  # stop is visible as a flattening rather than as a cut.
  h$role  <- "candidate"
  h$label <- ""
  path_rows <- integer(0)
  for (s in sort(unique(h$step))) {
    idx <- which(h$step == s)
    if (s == 0L) { pick <- idx[1L] }
    else if (s <= n_sel) { pick <- idx[match(selected[s], h$variable[idx])] }
    else { pick <- idx[if (minimise) which.min(h$score[idx]) else which.max(h$score[idx])] }
    if (length(pick) == 1L && !is.na(pick)) path_rows <- c(path_rows, pick)
  }
  h$role[path_rows]  <- "path"
  h$label[path_rows] <- h$variable[path_rows]
  df <- data.frame(panel = metric, x = h$step, y = as.numeric(h$score),
                   role = h$role, label = h$label, stringsAsFactors = FALSE)
  chosen <- if (n_sel > 0L && n_sel %in% h$step) {
    i <- path_rows[h$step[path_rows] == n_sel]
    data.frame(panel = metric, x = n_sel, y = as.numeric(h$score[i[1L]]))
  } else data.frame(panel = character(0), x = numeric(0), y = numeric(0))

  has_null <- any(h$step == 0L)
  holdout  <- suppressWarnings(as.numeric(x$score_holdout %||% NA_real_))
  caption <- paste(c(
    sprintf("%s (%s is better); grey points are the other candidates at each step",
            metric, if (minimise) "lower" else "higher"),
    if (has_null) "Step 0 is the intercept-only model"
    else "The intercept-only model could not be scored, so the path starts at the first variable",
    if (n_sel == 0L) "Nothing was selected"
    else sprintf("Selected: %s", paste(selected, collapse = ", ")),
    if (is.finite(holdout))
      sprintf("Blue: hold-out score of the selected set on the estimation half (%.3g)", holdout)
  ), collapse = "\n")

  p <- .plot_sweep(df, chosen = chosen, title = "Forward selection path",
                   subtitle = sprintf("Selection criterion, %s over %s folds",
                                      x$params$method %||% "cross-validation",
                                      format(x$params$k %||% "k")),
                   caption = caption, x_lab = "Step (variables in the model)",
                   y_lab = metric) +
    ggplot2::scale_x_continuous(breaks = sort(unique(df$x)))
  if (is.finite(holdout) && n_sel > 0L)
    p <- p + ggplot2::geom_point(data = data.frame(x = n_sel, y = holdout),
                                 ggplot2::aes(x = .data$x, y = .data$y),
                                 colour = "#2166AC", shape = 17, size = 3.5)
  p
}


#' Plot a GWR model selection
#'
#' \code{\link{gwr_model_selection}()} ranks every model it evaluated on
#' AICc.  This draws each model's criterion against its number of
#' predictors, the winner in red, so the gap between the best model and the
#' runners-up --- which the ranked table shows only as numbers --- is read
#' as a shape: a winner well below the rest was chosen by the data, a
#' winner a fraction of an AICc unit ahead of three others was chosen by the
#' tie-break.  The criterion is in-sample and the caption carries the label
#' \code{gwr_model_selection()} attached to it, including any note that it
#' was read positionally.
#'
#' @param x A \code{gwr_model_selection} object.
#' @param ... Ignored.
#' @return A \code{ggplot} object.
#' @family plotting
#' @examples
#' \donttest{
#' if (requireNamespace("GWmodel", quietly = TRUE) &&
#'     requireNamespace("ggplot2", quietly = TRUE)) {
#'   # sel <- gwr_model_selection(dat, "z", c("a", "b", "c"))
#'   # plot(sel)
#' }
#' }
#' @export
plot.gwr_model_selection <- function(x, ...) {
  .need_ggplot("plot.gwr_model_selection()")
  tab <- x$table
  if (!is.data.frame(tab) || !all(c("rank", "n_vars", "criterion") %in% names(tab)))
    stop("plot.gwr_model_selection(): `x` must be the object returned by ",
         "gwr_model_selection().", call. = FALSE)
  v  <- suppressWarnings(as.numeric(tab$criterion))
  ok <- is.finite(v)
  if (!any(ok))
    stop("plot.gwr_model_selection(): no model has a finite criterion.", call. = FALSE)
  crit_label <- as.character(x$criterion %||% "AICc")
  models <- data.frame(x = as.numeric(tab$n_vars[ok]), y = v[ok],
                       rank = as.integer(tab$rank[ok]),
                       vars = as.character(tab$variables[ok]),
                       stringsAsFactors = FALSE)
  # Every model as a faint candidate; the best model of each size as the
  # path (what one more variable buys); the winner labelled with its terms.
  cand <- data.frame(panel = crit_label, x = models$x, y = models$y,
                     role = "candidate", label = "", stringsAsFactors = FALSE)
  best_i <- vapply(split(seq_len(nrow(models)), models$x),
                   function(i) i[which.min(models$y[i])], integer(1))
  path <- data.frame(panel = crit_label, x = models$x[best_i], y = models$y[best_i],
                     role = "path",
                     label = ifelse(models$rank[best_i] == 1L, models$vars[best_i], ""),
                     stringsAsFactors = FALSE)
  df <- rbind(cand, path)
  win <- which(models$rank == 1L)[1L]
  chosen <- data.frame(panel = crit_label, x = models$x[win], y = models$y[win])
  n_finite <- sum(ok)
  gap <- if (n_finite >= 2L) sort(v[ok])[2L] - min(v[ok]) else NA_real_
  .plot_sweep(df, chosen = chosen,
              title = "GWR model selection",
              subtitle = sprintf("%d models evaluated at bandwidth %s (%s)%s",
                                 nrow(tab), format(signif(x$bandwidth, 4)),
                                 if (isTRUE(x$adaptive)) "adaptive" else "fixed",
                                 if (is.finite(gap)) sprintf("; the winner leads the runner-up by %.2f", gap) else ""),
              caption = paste0("Criterion: ", crit_label, " (in-sample; lower is better). ",
                               "Line: the best model of each size."),
              x_lab = "Number of predictors", y_lab = crit_label) +
    ggplot2::scale_x_continuous(breaks = sort(unique(df$x)))
}
