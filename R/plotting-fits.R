# =============================================================================
# Plot methods for fitted models and fold schemes
# =============================================================================

.need_ggplot <- function(what) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop(what, " requires the 'ggplot2' package. Install it with ",
         "install.packages(\"ggplot2\").", call. = FALSE)
}


# ggplot2 clips a title, subtitle or caption that is wider than the figure
# instead of wrapping it, and the loss is silent: on a help page a sentence
# ending "exceeds the largest lag fitted (682)" was drawn as "exceeds the
# largest lag fit".  These labels are built from numbers whose width is not
# known until the fit is in hand, so the break has to be found at draw time.
# The widths passed in were set from the six-inch figure the reference pages
# and the articles are drawn at, measuring the drawn text rather than
# counting characters.  The subtitle face runs 0.074 to 0.078 inches a
# character, so about 68 fill the 5.3 inches left of it.  The caption face
# runs 0.058 to 0.060, against a width the panel sets rather than the figure:
# a right-hand legend pulls the caption's right edge in, from 5.9 inches to
# 5.1 on the variogram plot, and the text runs leftwards from there.  Both
# callers pass a count short of the limit, which is the room a run of wide
# characters needs.  Breaks already in the string are kept.
.wrap_label <- function(x, width) {
  if (is.null(x) || length(x) != 1L || is.na(x) || !nzchar(x)) return(x)
  parts <- strsplit(as.character(x), "\n", fixed = TRUE)[[1]]
  wrapped <- vapply(parts,
                    function(s) paste(strwrap(s, width = width), collapse = "\n"),
                    character(1), USE.NAMES = FALSE)
  paste(wrapped, collapse = "\n")
}


#' Plot a fitted spatial model
#'
#' Diagnostic plots for a \code{spatial_fit}.  The package previously shipped
#' \code{print()} and \code{summary()} methods but no \code{plot()}, so the
#' checks most likely to reveal a problem (is there structure left in the
#' residuals, and where is it) had to be written by hand each time.
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
#'       fitted model overlaid, so the fit can be judged rather than trusted,
#'       and (unless \code{response = FALSE}) the variogram of the
#'       response itself on the same points and lags, drawn hollow with a
#'       dashed fit.  The gap between the two curves is the spatial structure
#'       the model absorbed: a residual sill well below the response sill
#'       means most of it, two curves that coincide mean none.  When both
#'       effective ranges were identified the caption gives the residual sill
#'       as a share of the response sill and the two ranges; when either
#'       variogram reached no sill the caption says so and compares nothing,
#'       because a sill the data never reached is not a number to divide by.
#'       The residual range is expected to come out shorter and the residual
#'       sill lower even when the model is right, because residuals of a
#'       fitted trend understate the variogram (see
#'       \code{\link{estimate_sac_range}()}, "Detrending and the
#'       residual-variogram bias").
#'       The distance axis is labelled in the units of the CRS the variogram
#'       was actually fitted in, which is not necessarily the fit's own CRS
#'       (lon/lat data are projected first). A single-direction fit names its
#'       azimuth in the title; a fit that identified no range says why in the
#'       subtitle, since the overlaid model line is then not a fit to believe.
#'       Requires 'gstat'.}
#'     \item{\code{"coefficients"}}{For a GWR fit only: the local coefficient
#'       of one \code{term} mapped at the training locations, which is the
#'       reason to fit GWR at all.  Locations where the local design is
#'       collinear (the kernel-weighted window's scaled condition index is
#'       above 30, or the window is singular) are drawn hollow and grey
#'       (\code{mask = TRUE}), because the smooth surface a naive map draws
#'       over them is the picture of an unstable estimate, not of a
#'       relationship; the subtitle counts them.  The condition indices are
#'       the fit's \code{info$local_collinearity}, computed for every
#'       location when the model was fitted.  A diverging scale centred on
#'       zero is used when the coefficient changes sign, otherwise a
#'       sequential one.}
#'   }
#' @param response Logical, default \code{TRUE}: for \code{type =
#'   "variogram"}, overlay the response's own variogram.  Ignored by the
#'   other types.
#' @param term For \code{type = "coefficients"}: which local coefficient to
#'   map, one of the names \code{coef(x)} returns.  Default \code{NULL}: the
#'   first predictor.  Ignored by the other types.
#' @param mask For \code{type = "coefficients"}: whether to draw locations
#'   whose local design is collinear (scaled condition index of the
#'   kernel-weighted window above 30, or singular) as hollow grey points
#'   instead of colouring them by a coefficient that is not to be believed
#'   there.  Default \code{TRUE}.  Locations whose coefficient is non-finite
#'   are masked either way.
#' @param ... Ignored.
#' @return A \code{ggplot} object.
#' @family plotting
#' @examples
#' # Works on any spatial_fit; a forest keeps the example free of the optional
#' # GWR/Stan backends.
#' if (requireNamespace("ranger", quietly = TRUE) &&
#'     requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(sf)
#'   # price depends on elevation, which the forest sees, and on a spatially
#'   # correlated field it does not: that field is what the residual plots
#'   # are there to find.
#'   set.seed(2)
#'   n <- 120
#'   xy <- data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
#'                    elev = rnorm(n))
#'   D  <- as.matrix(dist(xy[, c("x", "y")]))
#'   xy$price <- 10 + 2 * xy$elev +
#'     as.numeric(t(chol(exp(-D / 100) + diag(0.3, n))) %*% rnorm(n))
#'   pts <- st_as_sf(xy, coords = c("x", "y"), crs = 32632)
#'   fit <- fit_rf_model(pts, "price", "elev", num_trees = 100, seed = 1)
#'   # print() each one: inside a braced block only the last value is drawn.
#'   print(plot(fit, type = "residuals"))          # structure left in the residuals
#'   print(plot(fit, type = "observed_predicted"))
#'   if (requireNamespace("gstat", quietly = TRUE))
#'     plot(fit, type = "variogram")   # residual and response variograms compared
#' }
#' @seealso \code{\link{coef.gwr_fit}()} for the coefficients themselves.
#' @export
plot.spatial_fit <- function(x, type = c("residuals", "observed_predicted",
                                         "variogram", "coefficients"),
                             response = TRUE, term = NULL, mask = TRUE, ...) {
  type <- match.arg(type)
  .need_ggplot("plot.spatial_fit()")
  if (type == "coefficients")
    return(.plot_gwr_coefficients(x, term = term, mask = mask))

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
  if (is.null(attr(sac, "variogram", exact = TRUE)))
    stop("plot.spatial_fit(): the residual variogram could not be computed; ",
         "there may be too few finite residuals, or the fit's data_sf may ",
         "carry no usable geometry.", call. = FALSE)
  # The response's own variogram, on the same points and lags: the gap
  # between the two curves is the spatial structure the model absorbed.  Its
  # estimate is diagnostic-only, so its log lines stay off the console.
  resp <- NULL
  if (isTRUE(response) && is.character(x$response_var) &&
      x$response_var %in% names(dat)) {
    resp <- try(logger::with_log_threshold(
      estimate_sac_range(dat, x$response_var),
      threshold = logger::FATAL, namespace = "spatialkit", index = 2),
      silent = TRUE)
    if (inherits(resp, "try-error") || is.null(attr(resp, "variogram", exact = TRUE))) resp <- NULL
  }
  .draw_sac_variogram(sac, what = "Residual variogram", overlay = resp,
                      overlay_label = sprintf("Response (%s)", x$response_var))
}


#' Map one local GWR coefficient, masking collinear windows
#'
#' @param x A \code{gwr_fit}.
#' @param term Coefficient name, or \code{NULL} for the first predictor.
#' @param mask Mask collinear windows.
#' @return A \code{ggplot} object.
#' @keywords internal
#' @noRd
.plot_gwr_coefficients <- function(x, term = NULL, mask = TRUE) {
  if (!inherits(x, "gwr_fit"))
    stop("plot.spatial_fit(type = \"coefficients\"): local coefficients exist ",
         "for GWR fits only (class gwr_fit); this fit is a ",
         paste(setdiff(class(x), "spatial_fit"), collapse = "/"), ".", call. = FALSE)
  dat <- x$data_sf
  if (!inherits(dat, "sf"))
    stop("plot.spatial_fit(): the fit carries no training geometry.", call. = FALSE)
  cf <- stats::coef(x)
  if (!is.data.frame(cf) || !nrow(cf))
    stop("plot.spatial_fit(type = \"coefficients\"): coef() returned no local ",
         "coefficients.", call. = FALSE)
  if (nrow(cf) != nrow(dat))
    stop(sprintf(paste0("plot.spatial_fit(type = \"coefficients\"): %d rows of ",
                        "local coefficients for %d training locations; the two ",
                        "cannot be aligned."), nrow(cf), nrow(dat)), call. = FALSE)
  if (is.null(term)) {
    term <- intersect(x$predictor_vars, names(cf))[1L]
    if (is.na(term)) term <- names(cf)[1L]
  }
  if (!is.character(term) || length(term) != 1L || !(term %in% names(cf)))
    stop(sprintf("plot.spatial_fit(type = \"coefficients\"): `term` must be one of %s.",
                 paste(sQuote(names(cf)), collapse = ", ")), call. = FALSE)
  vals <- suppressWarnings(as.numeric(cf[[term]]))
  lc <- x$info$local_collinearity
  cn_bad <- if (is.data.frame(lc) && nrow(lc) == nrow(dat))
    (!is.finite(lc$cn) | lc$cn > 30) else rep(FALSE, nrow(dat))
  non_finite <- !is.finite(vals)
  masked <- non_finite | (isTRUE(mask) & cn_bad)
  if (all(masked))
    stop("plot.spatial_fit(type = \"coefficients\"): every location is masked ",
         "(non-finite coefficient, or a collinear local design at all of them); ",
         "there is no surface to draw.", call. = FALSE)

  dat$.coef   <- vals
  dat$.masked <- masked
  shown <- dat[!masked, , drop = FALSE]
  hidden <- dat[masked, , drop = FALSE]
  rng <- range(shown$.coef, na.rm = TRUE)
  diverging <- rng[1] < 0 && rng[2] > 0

  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = shown, ggplot2::aes(colour = .data$.coef), size = 1.8)
  if (nrow(hidden))
    p <- p + ggplot2::geom_sf(data = hidden, shape = 1, colour = "grey55", size = 1.8)
  p <- p + if (diverging) {
    lim <- max(abs(rng))
    ggplot2::scale_colour_gradient2(low = "#2166AC", mid = "grey92", high = "#B2182B",
                                    midpoint = 0, limits = c(-lim, lim), name = term)
  } else {
    ggplot2::scale_colour_viridis_c(name = term)
  }
  n_cn <- sum(cn_bad & !non_finite); n_nf <- sum(non_finite)
  subtitle <- if (any(masked))
    sprintf("%d of %d locations masked (hollow):\n%s", sum(masked), nrow(dat),
            paste(c(if (isTRUE(mask) && n_cn) sprintf("%d with a collinear local design (condition index > 30)", n_cn),
                    if (n_nf) sprintf("%d with a non-finite coefficient", n_nf)),
                  collapse = "; "))
  else if (isTRUE(mask) && is.data.frame(lc))
    "No location masked: every local design is well conditioned"
  else if (!isTRUE(mask) && any(cn_bad))
    sprintf("mask = FALSE: %d location(s) with a collinear local design are drawn as if reliable", sum(cn_bad))
  else "Collinearity not surveyed (fewer than two numeric predictors)"
  p + ggplot2::labs(
    title = sprintf("Local coefficient of %s", term),
    subtitle = subtitle,
    caption = sprintf("%s bandwidth %s, %s kernel%s",
                      if (isTRUE(x$info$adaptive)) "Adaptive" else "Fixed",
                      format(signif(x$info$bandwidth, 4)), x$info$kernel,
                      if (isTRUE(x$info$bandwidth_is_fallback)) " (fallback bandwidth)" else "")) +
    ggplot2::theme_minimal()
}


#' Plot an estimated spatial autocorrelation range
#'
#' Draws the empirical variogram that \code{\link{estimate_sac_range}()}
#' attaches to its result, with the fitted model and the effective range
#' overlaid where a range was identified, and a subtitle saying why not where
#' it was not.  A variogram that never reaches a sill, or that has no spatial
#' structure at the lags resolved, is the single most useful thing to
#' \emph{see} when a range comes back \code{NA}, so those cases are drawn
#' rather than refused.
#'
#' Nothing is recomputed: the plot reads the \code{variogram},
#' \code{variogram_model}, \code{crs}, \code{directional} and
#' \code{anisotropy_used} attributes the estimate already carries.  The
#' distance axis is in the units of the CRS the variogram was actually fitted
#' in, which \code{estimate_sac_range()} may have chosen itself for lon/lat
#' input.
#'
#' @param x An object of class \code{sac_range}, as returned by
#'   \code{\link{estimate_sac_range}()}.
#' @param ... Ignored.
#' @return A \code{ggplot} object.
#' @seealso \code{\link{plot.spatial_fit}(type = "variogram")}, which draws
#'   the same picture for a fitted model's residuals.
#' @family plotting
#' @examples
#' if (requireNamespace("gstat", quietly = TRUE) &&
#'     requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(sf)
#'   # An exponential field with range parameter 150 (effective range about
#'   # 450 m) and a nugget of half the sill, so the fitted model, the nugget
#'   # and the range line all have something to show.
#'   set.seed(3)
#'   n <- 250
#'   xy <- data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000))
#'   D  <- as.matrix(dist(xy))
#'   xy$z <- as.numeric(t(chol(exp(-D / 150) + diag(0.5, n))) %*% rnorm(n))
#'   pts <- st_as_sf(xy, coords = c("x", "y"), crs = 32632)
#'   r <- estimate_sac_range(pts, response_var = "z")
#'   plot(r)
#'
#'   # A field whose variogram never reaches a sill: the plot still draws it,
#'   # and the subtitle says why no range is marked.
#'   xy$trend <- sin(xy$x / 400) + rnorm(n, sd = 0.2)
#'   r2 <- estimate_sac_range(st_as_sf(xy, coords = c("x", "y"), crs = 32632),
#'                            response_var = "trend")
#'   plot(r2)
#' }
#' @export
plot.sac_range <- function(x, ...) {
  .need_ggplot("plot.sac_range()")
  if (!inherits(x, "sac_range"))
    stop("plot.sac_range(): `x` must be an object returned by ",
         "estimate_sac_range().", call. = FALSE)
  if (is.null(attr(x, "variogram", exact = TRUE)))
    stop("plot.sac_range(): this estimate carries no empirical variogram to ",
         "draw. estimate_sac_range() returns a bare NA, with nothing attached, ",
         "when the input has too few usable points or the variogram could not ",
         "be computed at all.", call. = FALSE)
  .draw_sac_variogram(x, what = "Empirical variogram")
}


#' Draw a sac_range object's variogram
#'
#' The one drawing routine behind \code{plot.sac_range()} and
#' \code{plot.spatial_fit(type = "variogram")}.  \code{what} is the plain
#' title used when the all-pairs variogram is drawn; when anisotropy was
#' established and a single direction is being shown, the title says so and
#' names the azimuth, because calling a quarter of the pairs by the plain
#' name overstates the structure the plot exists to show.
#'
#' @param sac A \code{sac_range} with a non-\code{NULL} \code{variogram}
#'   attribute.
#' @param what Character(1) title stem.
#' @param overlay Optional second \code{sac_range} (the response's, when
#'   \code{sac} is the residuals') drawn hollow with a dashed model line, so
#'   the structure the model absorbed is the gap between the two.
#' @param overlay_label Legend label for the overlay.
#' @return A \code{ggplot} object.
#' @keywords internal
#' @noRd
.draw_sac_variogram <- function(sac, what = "Empirical variogram",
                                overlay = NULL, overlay_label = "Response") {
  vg  <- attr(sac, "variogram", exact = TRUE)
  vm  <- attr(sac, "variogram_model")
  if (!is.null(overlay) && is.null(attr(overlay, "variogram", exact = TRUE))) overlay <- NULL

  # The axis is in the units of the CRS the VARIOGRAM was fitted in, which
  # estimate_sac_range() chose with ensure_projected() -- not necessarily the
  # caller's own CRS.  Labelling it "CRS units" for a lon/lat input named
  # degrees while the numbers were metres of an auto-chosen UTM zone that
  # appeared nowhere on the plot.  estimate_sac_range() returns the CRS for
  # exactly this.
  vg_units <- tryCatch({
    u <- sf::st_crs(attr(sac, "crs"))$units_gdal
    if (is.null(u) || is.na(u) || !nzchar(u)) "CRS units" else u
  }, error = function(e) "CRS units")

  # Which variogram is this?  When anisotropy was established,
  # estimate_sac_range() returns the single azimuth with the widest range --
  # about a quarter of the point pairs -- and calling that plainly by `what`
  # overstates the structure for the plot whose whole purpose is to show how
  # much there is.
  vg_title <- what
  if (isTRUE(attr(sac, "anisotropy_used"))) {
    dir_r <- attr(sac, "directional")
    az <- if (!is.null(dir_r) && length(dir_r))
      names(dir_r)[which.max(replace(dir_r, is.na(dir_r), -Inf))] else NA
    vg_title <- if (is.na(az)) sprintf("%s (widest direction only)", what)
      else sprintf("%s, %s\u00b0 \u00b1 22.5\u00b0 (the widest of four directions)",
                   what, az)
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

  # The overlay: hollow points and a dashed model line for the second
  # variogram, sharing the axes.  Nothing about the main curve changes.
  overlay_caption <- NULL
  if (!is.null(overlay)) {
    ovg <- attr(overlay, "variogram", exact = TRUE)
    ovm <- attr(overlay, "variogram_model")
    p <- p + ggplot2::geom_point(data = ovg,
                                 ggplot2::aes(x = .data$dist, y = .data$gamma,
                                              size = .data$np),
                                 shape = 1, colour = "grey35", alpha = 0.8)
    if (!is.null(ovm)) {
      oline <- try(gstat::variogramLine(ovm, maxdist = max(vg$dist, ovg$dist, na.rm = TRUE),
                                        n = 200), silent = TRUE)
      if (!inherits(oline, "try-error"))
        p <- p + ggplot2::geom_line(data = oline,
                                    ggplot2::aes(x = .data$dist, y = .data$gamma),
                                    colour = "grey35", linetype = "dashed",
                                    linewidth = 0.8)
    }
    # The sills are compared only when BOTH ranges were identified: a model
    # whose range ran past the lags fitted has a sill the data never reached,
    # and a ratio of two such numbers reads as a finding while being noise.
    sill_of <- function(m) if (is.data.frame(m) && "psill" %in% names(m))
      sum(as.numeric(m$psill), na.rm = TRUE) else NA_real_
    s_res <- sill_of(vm); s_resp <- sill_of(ovm)
    both_ok <- is.finite(sac) && is.finite(overlay) &&
      is.finite(s_res) && is.finite(s_resp) && s_resp > 0
    overlay_caption <- paste(c(
      sprintf("Hollow points, dashed line: %s. Filled points, solid line: %s.",
              overlay_label, tolower(what)),
      if (both_ok)
        sprintf("Residual sill is %.0f%% of the response sill; effective ranges %.0f (residuals) and %.0f (response).",
                100 * s_res / s_resp, as.numeric(sac), as.numeric(overlay))
      else sprintf("Sills not compared: %s.",
                   if (!is.finite(sac) && !is.finite(overlay))
                     "neither variogram reached an identified sill"
                   else if (!is.finite(sac))
                     "the residual variogram reached no identified sill"
                   else if (!is.finite(overlay))
                     "the response variogram reached no identified sill"
                   else "a variogram model could not be fitted")
    ), collapse = "\n")
    p <- p + ggplot2::labs(caption = .wrap_label(overlay_caption, 72))
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
      subtitle = .wrap_label(
        if (identical(reason, "variogram model did not converge"))
          sprintf(paste0("No effective range: the fit stopped at gstat's ",
                         "iteration limit, so the range it reports (%.0f) is ",
                         "where the optimiser halted, not a fitted parameter."),
                  attr(sac, "rejected_range"))
        else if (isTRUE(grepl("^no variogram model", reason)))
          # Both model fits singular: the picture of residuals with no spatial
          # structure at the lags resolved.  There is no model line to draw.
          paste0("No effective range: no variogram model could be fitted ",
                 "(a flat, nugget-only variogram -- no spatial structure at ",
                 "these lags).")
        else if (identical(reason, "empirical variogram decreases with distance"))
          # The shape of a periodic structure or of a variance that differs
          # between a dense cluster and the rest -- not of a trend, which rises
          # without a sill.
          sprintf(paste0("No effective range: the semivariance falls with ",
                         "distance over the shorter lags, so the fitted range ",
                         "(%.0f) is not identified (periodic structure, or a ",
                         "variance that differs across the layer)."),
                  attr(sac, "rejected_range"))
        else if (identical(reason, "fitted range exceeds the largest lag fitted"))
          sprintf(paste0("No effective range: the fitted range (%.0f) ",
                         "exceeds the largest lag fitted (%.0f), so the ",
                         "variogram never reached a sill."),
                  attr(sac, "rejected_range"),
                  attr(sac, "cutoff_dist"))
        else
          # A reason this method does not know by name: say it verbatim rather
          # than caption it with another case's sentence.
          sprintf("No effective range: %s.", as.character(reason)),
      65))
  }
  p
}


#' Map a cross-validation fold scheme
#'
#' Shows which fold each observation belongs to.  This is the fastest way to
#' see whether spatial blocks are actually separating the data, or whether the
#' blocks are smaller than the autocorrelation range and therefore leaking.
#' For \code{"block_kfold"} folds the block outlines are drawn too, from
#' \code{folds$params$blocks}, so a fold can be seen to be one region or
#' several and an empty block can be seen to be empty.  The subtitle states
#' the parameter that decides whether the scheme leaks, read from
#' \code{folds$params}: the block size (and the estimated range when
#' \code{auto_range} found one), the buffer, the number of location groups,
#' or the median NNDM exclusion.
#'
#' @param folds A list returned by \code{make_folds()}.
#' @param points_sf The \code{sf} layer the folds were built from.
#' @param boundary Optional polygonal \code{sf}/\code{sfc} to draw underneath.
#' @param blocks Logical; draw the block polygons when \code{folds} carries
#'   them.  Default \code{TRUE}.
#' @return A \code{ggplot} object.
#' @family plotting
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(sf)
#'   set.seed(1)
#'   n <- 80
#'   pts <- st_as_sf(
#'     data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000)),
#'     coords = c("x", "y"), crs = 32632
#'   )
#'   f <- make_folds(pts, k = 5, method = "block_kfold", block_size = 300)
#'   plot_folds(f, pts)
#' }
#' @export
plot_folds <- function(folds, points_sf, boundary = NULL, blocks = TRUE) {
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
  # The block design, when the folds carry it (block_kfold): outlines under
  # the points, in the CRS the folds were built in -- coord_sf() brings the
  # layers to one CRS.
  blk <- folds$params$blocks
  if (isTRUE(blocks) && inherits(blk, "sf") && nrow(blk) > 0L)
    p <- p + ggplot2::geom_sf(data = sf::st_geometry(blk),
                              fill = NA, colour = "grey45", linewidth = 0.25)
  p +
    ggplot2::geom_sf(data = dat, ggplot2::aes(colour = .data$.fold)) +
    ggplot2::scale_colour_viridis_d(name = "Fold", na.value = "grey80") +
    ggplot2::labs(
      title = sprintf("%s folds (k = %s)", folds$method, folds$k),
      subtitle = .fold_subtitle(folds)
    ) +
    ggplot2::theme_minimal()
}

# The subtitle states the one parameter that decides whether the scheme
# leaks, read from what make_folds() recorded.  It used to be a fixed line
# about block size, which was wrong for every method that has no blocks.
.fold_subtitle <- function(folds) {
  prm <- folds$params
  fin <- function(v) length(v) == 1L && is.finite(suppressWarnings(as.numeric(v)))
  num <- function(v) format(signif(as.numeric(v), 3), big.mark = ",")
  units <- if (is.character(prm$crs) && length(prm$crs) == 1L && nzchar(prm$crs))
    sprintf(" (%s units)", prm$crs) else ""
  switch(as.character(folds$method),
    random_kfold =
      "Random folds: a held-out point's neighbours stay in the training set",
    block_kfold = {
      size <- if (isTRUE(prm$blocks_supplied))
        sprintf("%s supplied blocks", format(prm$n_blocks %||% nrow(prm$blocks)))
      else if (fin(prm$block_size))
        sprintf("Block size %s%s", num(prm$block_size), units)
      else if (fin(prm$grid_nx) && fin(prm$grid_ny))
        sprintf("%s x %s block grid", format(prm$grid_nx), format(prm$grid_ny))
      else "Blocked folds"
      # Two lines: one long one is clipped at the width a help page or a
      # vignette draws these at.
      range <- if (!fin(prm$sac_range))
        "Blocks smaller than the autocorrelation range leak"
      else if (fin(prm$block_size) &&
               isTRUE(all.equal(as.numeric(prm$block_size),
                                as.numeric(prm$sac_range))))
        "Sized from the estimated autocorrelation range"
      else
        sprintf("Estimated range %s; blocks below it leak", num(prm$sac_range))
      paste(size, range, sep = "\n")
    },
    buffered_loo = if (fin(prm$buffer))
      sprintf("Leave-one-out with a %s buffer%s", num(prm$buffer), units)
    else "Leave-one-out with an exclusion buffer",
    leave_location_out = if (fin(prm$n_groups) && is.character(prm$group_var))
      sprintf("%s groups of `%s` across the folds", format(prm$n_groups),
              prm$group_var[1L])
    else "One location group at a time",
    nndm = if (fin(prm$median_buffer))
      sprintf("Distance-matched exclusion, median %s%s", num(prm$median_buffer),
              units)
    else "Distance-matched exclusion (NNDM)",
    NULL)
}
