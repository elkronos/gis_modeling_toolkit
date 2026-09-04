# =============================================================================
# Plot methods for fitted models and fold schemes
# =============================================================================

.need_ggplot <- function(what) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop(what, " requires the 'ggplot2' package. Install it with ",
         "install.packages(\"ggplot2\").", call. = FALSE)
}


#' Plot a fitted spatial model
#'
#' Diagnostic plots for a \code{spatial_fit}.  The package previously shipped
#' \code{print()} and \code{summary()} methods but no \code{plot()}, so the
#' checks most likely to reveal a problem -- is there structure left in the
#' residuals, and where is it -- had to be written by hand each time.
#'
#' @param x A \code{spatial_fit}.
#' @param type One of:
#'   \describe{
#'     \item{\code{"residuals"}}{Residuals mapped at the training locations.
#'       Spatial structure here is the signal that the model has not captured
#'       the autocorrelation.}
#'     \item{\code{"observed_predicted"}}{Observed against fitted, with a 1:1
#'       reference line.}
#'     \item{\code{"variogram"}}{Empirical variogram of the residuals with the
#'       fitted model overlaid, so the fit can be judged rather than trusted.
#'       The distance axis is labelled in the units of the CRS the variogram
#'       was actually fitted in, which is not necessarily the fit's own CRS
#'       (lon/lat data are projected first). A single-direction fit names its
#'       azimuth in the subtitle; a fit that did not converge says so in the
#'       caption, since the overlaid model line is then not a fit to believe.
#'       Requires 'gstat'.}
#'   }
#' @param ... Ignored.
#' @return A \code{ggplot} object.
#' @family plotting
#' @examples
#' \donttest{
#' # Works on any spatial_fit; a forest keeps the example free of the optional
#' # GWR/Stan backends.
#' if (requireNamespace("ranger", quietly = TRUE) &&
#'     requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(sf)
#'   set.seed(1)
#'   n <- 120
#'   pts <- st_as_sf(
#'     data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), elev = rnorm(n)),
#'     coords = c("x", "y"), crs = 32632
#'   )
#'   pts$price <- 10 + 0.01 * st_coordinates(pts)[, 1] + 2 * pts$elev + rnorm(n)
#'   fit <- fit_rf_model(pts, "price", "elev", num_trees = 100, seed = 1)
#'   plot(fit, type = "residuals")
#'   plot(fit, type = "observed_predicted")
#'   if (requireNamespace("gstat", quietly = TRUE))
#'     plot(fit, type = "variogram")
#' }
#' }
#' @export
plot.spatial_fit <- function(x, type = c("residuals", "observed_predicted",
                                         "variogram"), ...) {
  type <- match.arg(type)
  .need_ggplot("plot.spatial_fit()")

  res <- try(stats::residuals(x), silent = TRUE)
  if (inherits(res, "try-error") || !is.numeric(res))
    stop("plot.spatial_fit(): could not extract residuals from this fit.",
         call. = FALSE)
  fit_vals <- try(stats::fitted(x), silent = TRUE)

  dat <- x$data_sf
  if (!inherits(dat, "sf"))
    stop("plot.spatial_fit(): the fit carries no training geometry.",
         call. = FALSE)
  if (length(res) != nrow(dat))
    stop("plot.spatial_fit(): residual length (", length(res),
         ") does not match the training data (", nrow(dat), ").", call. = FALSE)
  # An all-NA residual vector (e.g. a bayesian_fit whose posterior_epred
  # failed) otherwise yields limits = c(Inf, -Inf), which ggplot builds without
  # complaint and paints entirely in na.value -- a uniformly grey map that
  # looks like a result.
  if (!any(is.finite(as.numeric(res))))
    stop("plot.spatial_fit(): the model produced no finite residuals, so there ",
         "is nothing to plot.", call. = FALSE)

  if (type == "residuals") {
    dat$.resid <- as.numeric(res)
    lim <- max(abs(dat$.resid), na.rm = TRUE)
    # A perfectly-fitting model gives lim = 0, so limits = c(0, 0): a
    # degenerate diverging scale whose breaks collapse onto a single value.
    # The all-NA case a few lines above is handled; this one fell through.
    if (!is.finite(lim) || lim <= 0) lim <- 1
    return(
      ggplot2::ggplot(dat) +
        ggplot2::geom_sf(ggplot2::aes(colour = .data$.resid)) +
        ggplot2::scale_colour_gradient2(
          low = "#2166AC", mid = "grey92", high = "#B2182B",
          midpoint = 0, limits = c(-lim, lim), name = "Residual"
        ) +
        ggplot2::labs(
          title = "Model residuals",
          subtitle = "Visible spatial structure means unmodelled autocorrelation"
        ) +
        ggplot2::theme_minimal()
    )
  }

  if (type == "observed_predicted") {
    if (inherits(fit_vals, "try-error") || !is.numeric(fit_vals))
      stop("plot.spatial_fit(): could not extract fitted values.", call. = FALSE)
    obs <- as.numeric(fit_vals) + as.numeric(res)
    df  <- data.frame(observed = obs, predicted = as.numeric(fit_vals))
    return(
      ggplot2::ggplot(df, ggplot2::aes(x = .data$observed, y = .data$predicted)) +
        ggplot2::geom_abline(slope = 1, intercept = 0,
                             linetype = "dashed", colour = "grey50") +
        ggplot2::geom_point(alpha = 0.6) +
        ggplot2::labs(title = "Observed vs predicted",
                      x = "Observed", y = "Predicted") +
        ggplot2::theme_minimal()
    )
  }

  # ---- variogram of residuals ----------------------------------------------
  if (!requireNamespace("gstat", quietly = TRUE))
    stop("plot.spatial_fit(type = \"variogram\") requires the 'gstat' package.",
         call. = FALSE)

  dat$.resid <- as.numeric(res)
  sac <- estimate_sac_range(dat, ".resid")
  vg  <- attr(sac, "variogram")
  vm  <- attr(sac, "variogram_model")
  if (is.null(vg))
    stop("plot.spatial_fit(): the residual variogram could not be fitted; ",
         "there may be too few finite residuals to model.", call. = FALSE)

  # The axis is in the units of the CRS the VARIOGRAM was fitted in, which
  # estimate_sac_range() chose with ensure_projected() -- not necessarily the
  # fit's own CRS.  Labelling it "CRS units" for a lon/lat fit named degrees
  # while the numbers were metres of an auto-chosen UTM zone that appeared
  # nowhere on the plot.  estimate_sac_range() returns the CRS for exactly this.
  vg_units <- tryCatch({
    u <- sf::st_crs(attr(sac, "crs"))$units_gdal
    if (is.null(u) || is.na(u) || !nzchar(u)) "CRS units" else u
  }, error = function(e) "CRS units")

  # Which variogram is this?  When anisotropy was established,
  # estimate_sac_range() returns the single azimuth with the widest range --
  # about a quarter of the point pairs -- and calling that plainly "Residual
  # variogram" overstates the leftover structure for the plot whose whole
  # purpose is to show how much there is.
  vg_title <- "Residual variogram"
  if (isTRUE(attr(sac, "anisotropy_used"))) {
    dir_r <- attr(sac, "directional")
    az <- if (!is.null(dir_r) && length(dir_r))
      names(dir_r)[which.max(replace(dir_r, is.na(dir_r), -Inf))] else NA
    vg_title <- if (is.na(az)) "Residual variogram (widest direction only)"
      else sprintf("Residual variogram, %s\u00b0 \u00b1 22.5\u00b0 (the widest of four directions)", az)
  }

  p <- ggplot2::ggplot(vg, ggplot2::aes(x = .data$dist, y = .data$gamma)) +
    ggplot2::geom_point(ggplot2::aes(size = .data$np), alpha = 0.7) +
    ggplot2::scale_size_continuous(name = "Pairs") +
    ggplot2::labs(title = vg_title,
                  x = sprintf("Distance (%s)", vg_units), y = "Semivariance") +
    ggplot2::theme_minimal()

  if (!is.null(vm)) {
    line <- try(gstat::variogramLine(vm, maxdist = max(vg$dist, na.rm = TRUE),
                                     n = 200), silent = TRUE)
    if (!inherits(line, "try-error"))
      p <- p + ggplot2::geom_line(data = line,
                                  ggplot2::aes(x = .data$dist, y = .data$gamma),
                                  colour = "#B2182B", linewidth = 0.8)
  }
  if (is.finite(sac)) {
    p <- p + ggplot2::geom_vline(xintercept = as.numeric(sac),
                                 linetype = "dotted", colour = "#B2182B") +
      ggplot2::labs(subtitle = sprintf("Effective range = %.1f %s",
                                       as.numeric(sac), vg_units))
  } else if (!is.null(attr(sac, "rejected_range"))) {
    # A variogram that never reaches a sill is the single most useful thing to
    # SEE, so draw it and label why no range is marked, rather than refusing.
    # estimate_sac_range() records WHICH of two reasons applied; hard-coding
    # the sill wording for both captioned a non-converged fit with a sentence
    # its own two numbers refute (range 11 "exceeds" a cutoff of 682).
    reason <- attr(sac, "rejected_reason")
    p <- p + ggplot2::labs(
      subtitle = if (identical(reason, "variogram model did not converge"))
        sprintf(paste0("No effective range: the fit stopped at gstat's ",
                       "iteration limit, so the range it reports (%.0f) is ",
                       "where the optimiser halted, not a fitted parameter."),
                attr(sac, "rejected_range"))
      else
        sprintf(paste0("No effective range: the fitted range (%.0f) ",
                       "exceeds the largest lag fitted (%.0f), so the ",
                       "variogram never reached a sill."),
                attr(sac, "rejected_range"),
                attr(sac, "cutoff_dist")))
  }
  p
}


#' Map a cross-validation fold scheme
#'
#' Shows which fold each observation belongs to.  This is the fastest way to
#' see whether spatial blocks are actually separating the data, or whether the
#' blocks are smaller than the autocorrelation range and therefore leaking.
#'
#' @param folds A list returned by \code{make_folds()}.
#' @param points_sf The \code{sf} layer the folds were built from.
#' @param boundary Optional polygonal \code{sf}/\code{sfc} to draw underneath.
#' @return A \code{ggplot} object.
#' @family plotting
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(sf)
#'   set.seed(1)
#'   n <- 80
#'   pts <- st_as_sf(
#'     data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000)),
#'     coords = c("x", "y"), crs = 32632
#'   )
#'   f <- make_folds(pts, k = 5, method = "block_kfold", block_size = 300)
#'   plot_folds(f, pts)
#' }
#' @export
plot_folds <- function(folds, points_sf, boundary = NULL) {
  .need_ggplot("plot_folds()")
  # `method` and `k` are read straight into the title; a folds list missing
  # either collapses sprintf() to character(0), which ggplot2 accepts and
  # renders as a plot with no title at all rather than erroring.
  if (!is.list(folds) || is.null(folds$assignment) ||
      is.null(folds$method) || is.null(folds$k))
    stop("plot_folds(): `folds` must be the list returned by make_folds() ",
         "(with `assignment`, `method` and `k`).", call. = FALSE)
  if (!inherits(points_sf, "sf"))
    stop("plot_folds(): `points_sf` must be an sf object.", call. = FALSE)

  asg <- folds$assignment
  if (!all(c("row_id", "fold") %in% names(asg)))
    stop("plot_folds(): the fold assignment table lacks `row_id`/`fold`.",
         call. = FALSE)

  dat <- points_sf
  if (!("..row_id" %in% names(dat))) dat$..row_id <- seq_len(nrow(dat))
  dat$.fold <- factor(asg$fold[match(dat$..row_id, asg$row_id)])

  n_assigned <- sum(!is.na(dat$.fold))
  if (n_assigned == 0L)
    stop("plot_folds(): no points matched the fold assignment; were `folds` ",
         "built from a different layer?", call. = FALSE)

  p <- ggplot2::ggplot()
  if (!is.null(boundary))
    p <- p + ggplot2::geom_sf(data = sf::st_geometry(boundary),
                              fill = NA, colour = "grey60")
  p +
    ggplot2::geom_sf(data = dat, ggplot2::aes(colour = .data$.fold)) +
    ggplot2::scale_colour_viridis_d(name = "Fold", na.value = "grey80") +
    ggplot2::labs(
      title = sprintf("%s folds (k = %s)", folds$method, folds$k),
      subtitle = "Blocks smaller than the autocorrelation range still leak"
    ) +
    ggplot2::theme_minimal()
}
