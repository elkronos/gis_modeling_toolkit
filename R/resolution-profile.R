# =============================================================================
# Resolution profile: every candidate level count, scored on several criteria
# at once, with the bounds the data impose and the region over which the
# criteria are flat.  determine_optimal_levels() answers "which k"; this
# answers "what does each k cost and buy", which is the question a
# tessellation's resolution has to be defended on.
# =============================================================================

#' Mean correlation between two random points in a rectangle
#'
#' Quadrature on a fixed \eqn{q \times q} grid over the unit square, scaled to
#' a rectangle of sides \code{sx} and \code{sy}: the mean of the correlation
#' function over the off-diagonal pairs of grid points.  For a square cell
#' this is the quantity Krige's additivity relation needs -- the variance of
#' the average of a stationary field over a block is the partial sill times
#' this mean correlation.
#'
#' @param cor_fn Correlation function of distance, from
#'   \code{.vgm_correlation_fn()} (it already excludes the nugget).
#' @param sx,sy Rectangle sides in the coordinate units.
#' @param q Grid points per side.  400 pairs of pairs at \code{q = 20}; the
#'   quadrature error is well below the sampling noise it feeds into.
#' @return A scalar in \eqn{[0, 1]}.
#' @keywords internal
#' @noRd
.rbar_rect <- function(cor_fn, sx, sy, q = 20L) {
  if (!is.finite(sx) || !is.finite(sy) || sx <= 0 || sy <= 0) return(NA_real_)
  u  <- (seq_len(q) - 0.5) / q
  g  <- as.matrix(expand.grid(x = u * sx, y = u * sy))
  d  <- as.numeric(stats::dist(g))            # off-diagonal pairs, once each
  r  <- as.numeric(cor_fn(d))
  mean(r[is.finite(r)])
}


#' Analytic reliability of cell means at a candidate resolution
#'
#' Krige's additivity relation splits the point-support variance of a field
#' over a domain \eqn{V} into a between-block and a within-block part,
#' \eqn{D^2(0|V) = D^2(v|V) + D^2(0|v)}, and both parts are computable from a
#' fitted variogram for any block size before anything is aggregated
#' (Cressie 1996).  With partial sill \eqn{c}, nugget \eqn{c_0}, mean
#' within-block correlation \eqn{\bar r(v)} and mean within-domain correlation
#' \eqn{\bar r(V)}:
#' \itemize{
#'   \item between-cell variance of the true block means:
#'     \eqn{c\,(\bar r(v) - \bar r(V))};
#'   \item sampling variance of a cell mean around its block mean, from
#'     \eqn{n} points treated as exchangeable within the block:
#'     \eqn{(c_0 + c\,(1 - \bar r(v))) / n}.
#' }
#' The reliability is the first over their sum: the share of the spread in
#' the cell means that is signal rather than sampling noise, the same
#' quantity as the shrinkage factor of Fay and Herriot (1979).  It falls
#' toward zero as cells grow to the domain (nothing left between them) and
#' as cells shrink to a handful of points (all noise), so it has an interior
#' optimum -- unlike "minimise the SE of the cell means", which bigger cells
#' always win.
#'
#' Validated by simulation before it went on the profile (exponential fields,
#' 600 points on a 1000-unit extent, 5-8 replicates each, empirical
#' reliability from the true block means of a fine grid): analytic and
#' empirical argmax 6 vs 5 at effective range 300 with nugget 0.3, 5 vs 5
#' with nugget 2, and at effective range 90 an analytic optimum of 13 inside
#' a band flat to within 2 percent from 6 to 24, against an empirical argmax
#' of 5 that varied from 2 to 6 across replicates.  The curve is broad: read
#' the flat region, not the argmax.
#'
#' @param L Number of cells.
#' @param area Domain area in coordinate units squared.
#' @param n_total Points in the layer.
#' @param nugget,psill Nugget and total partial sill of the variogram model.
#' @param cor_fn Correlation function of distance.
#' @param rbar_V Mean within-domain correlation (precomputed).
#' @return A scalar in \eqn{[0, 1]}, or \code{NA_real_}.
#' @keywords internal
#' @noRd
.reliability_at <- function(L, area, n_total, nugget, psill, cor_fn, rbar_V) {
  if (!is.finite(L) || L < 1 || !is.finite(area) || area <= 0) return(NA_real_)
  s  <- sqrt(area / L)
  n  <- n_total / L
  rv <- .rbar_rect(cor_fn, s, s)
  if (!is.finite(rv) || !is.finite(rbar_V)) return(NA_real_)
  between <- psill * (rv - rbar_V)
  noise   <- (nugget + psill * (1 - rv)) / n
  if (!is.finite(between) || !is.finite(noise) || between <= 0) return(0)
  between / (between + noise)
}


#' Score every candidate number of cells on several criteria at once
#'
#' A tessellation's resolution is the number of cells the points are cut
#' into, and no single number settles it: the within-cluster sum of squares
#' of the coordinates has an elbow, a piecewise-constant approximation of the
#' response has a Mallows \eqn{C_p}, the residual autocorrelation of the cell
#' means has a \eqn{z}, and the cell means have a reliability.  This function
#' computes all of them at every level of a ladder the data bound, and
#' returns the table, so that the level chosen --- by
#' \code{\link{select_resolution}()} or by eye --- can be defended with the
#' whole profile rather than one criterion's argmin.  Mallows (1973) presents
#' \eqn{C_p} itself as a display of the bias-variance trade-off across
#' candidates rather than a rule that picks one; this is that display, with
#' the other criteria alongside.
#'
#' @section The ladder and its bounds:
#' Cells are k-means clusters of the (projected) coordinates, fitted at each
#' level as the best of \code{nstart} k-means++ restarts (see
#' \code{\link{determine_optimal_levels}} for why).  Levels are spaced
#' logarithmically, because cell diameter scales as \eqn{L^{-1/2}}: a unit
#' step wastes fits at large \eqn{L} and starves resolution at small.  The
#' ladder runs from a floor to a ceiling the data impose.  The ceiling is
#' \code{floor(n / min_cell_n)}: cells with fewer than \code{min_cell_n} points
#' on average have too little support, and the model-aware criteria are not
#' computable below nine cells in any case.  The floor is
#' \code{ceiling(area / range^2)} when an autocorrelation range is available:
#' cells wider than the range average over more than one patch of the field.
#' When the floor exceeds the ceiling the data cannot support a tessellation
#' that respects their own correlation structure; that is reported as a
#' finding (a logged warning, and \code{attr(x, "bounds")$supported} is
#' \code{FALSE}) and the ladder runs from 2 to the ceiling anyway, so the
#' profile still shows what each level costs.
#'
#' @section The criteria, and how each behaved when measured:
#' \describe{
#'   \item{\code{elbow}}{The signed distance of the WSS curve below the chord
#'     from its first to its last level, the classical elbow statistic
#'     (larger is better).  Geometry only; it knows nothing of the response.}
#'   \item{\code{cp}}{Mallows' \eqn{C_p} of the piecewise-constant
#'     approximation of the response (or of its OLS residuals on
#'     \code{predictor_vars}) by cell means: \eqn{RSS(L)/n + 2 \tau^2 L / n},
#'     with \eqn{\tau^2} the nugget of the fitted variogram (lower is better).
#'     \strong{Measured on simulated exponential fields (600 points on a
#'     1000-unit extent, sill 1, 20 replicates): with a nugget of 0.3 its
#'     minimum sat at the support ceiling in every replicate at effective
#'     ranges of 90, 300 and 900; with a nugget of 2 it was interior (median
#'     44, range 8--66).}  On a smooth field with little noise the
#'     approximation keeps improving as cells shrink and the penalty is too
#'     small to stop it, so \eqn{C_p} says "as fine as the support allows"
#'     and \code{min_cell_n} is what is choosing; \code{select_resolution()}
#'     says so when that happens.  It becomes a genuine interior criterion
#'     only when the nugget is a large share of the sill.}
#'   \item{\code{moran_z}}{The standardised deviate of Moran's I on the
#'     residuals of the cell means regressed on the cell-mean predictors (an
#'     intercept alone when there are none): how much spatial structure the
#'     tessellation has left unexplained (\eqn{|z|} smaller is better).
#'     Calibrated and flat in \eqn{L} on a response with no structure --- see
#'     \code{\link{determine_optimal_levels}} --- so it separates levels only
#'     where structure remains.  \code{NA} at nine cells or fewer.}
#'   \item{\code{reliability}}{The share of the spread in the cell means that
#'     is between-cell signal rather than sampling noise, from the fitted
#'     variogram alone via Krige's additivity relation (Cressie 1996) ---
#'     the shrinkage factor of Fay and Herriot (1979) --- for square cells of
#'     the level's average area with the level's average point count (larger
#'     is better).  It has an interior optimum, and a broad one: validated
#'     against the empirical reliability of true block means on simulated
#'     fields, the analytic and empirical optima agreed to within a level or
#'     two where the empirical estimate was stable, and the band within 2
#'     percent of the maximum spanned a factor of 3--6 in \eqn{L}.  Read the
#'     flat region, not the argmax.  \code{NA} without a usable variogram.}
#' }
#' \code{cp} and \code{reliability} answer different questions --- how well
#' the cells represent the field, and whether the cell values are
#' distinguishable from noise --- and can disagree; both are shown so the
#' choice between them is made knowingly.
#'
#' @param data_sf An sf object of points (other geometries are reduced to
#'   representative points).
#' @param response_var Optional response column name (numeric or logical).
#'   Enables \code{cp}, \code{moran_z}, and --- via a variogram estimated
#'   from it --- the floor of the ladder and \code{reliability}.
#' @param predictor_vars Optional predictor column names (numeric or
#'   logical).  With them, \code{cp} scores the OLS residuals of the response
#'   on the predictors, the variogram is estimated from those residuals, and
#'   \code{moran_z} regresses the cell means on the cell-mean predictors.
#' @param levels Optional integer vector of level counts to score, replacing
#'   the ladder; values outside \code{[2, n - 1]} are dropped.
#' @param n_levels Number of levels on the ladder.  Default 20.
#' @param min_cell_n Minimum average number of points per cell that a level
#'   must keep; sets the ceiling.  Default 9.
#' @param sample_n Points are subsampled to this many before anything is
#'   fitted, as in \code{determine_optimal_levels()}.  Default 1500.  The
#'   support columns describe the subsample.
#' @param nstart k-means++ restarts per level.  Default 25.
#' @param seed RNG seed for the subsample and the restarts; restored
#'   afterwards.  Default 123.
#' @param sac Optional \code{sac_range} object from
#'   \code{\link{estimate_sac_range}()} to take the range, nugget and
#'   correlation function from --- pass one fitted with \code{detrend =
#'   "reml"}, say, or on a residual field of your choosing.  When
#'   \code{NULL} and a response is given, one is estimated on the subsample
#'   with the same \code{predictor_vars}.
#' @param quiet Logical; suppress this function's progress \code{message()}s.
#'   Default \code{TRUE}.
#' @param select_on \code{"all"} (default) profiles every point;
#'   \code{"split"} profiles one spatially blocked half and returns the other
#'   half as the set to estimate on, in the \code{"split"} attribute.  See
#'   the "Post-selection inference" section of
#'   \code{\link{determine_optimal_levels}}; the profile reads the response
#'   whenever \code{response_var} is given.
#' @return A data.frame of class \code{resolution_profile} with one row per
#'   level and columns \code{levels}, \code{wss}, \code{wss_spread} (relative
#'   spread of WSS across the restarts), \code{elbow}, \code{cell_n_min},
#'   \code{cell_n_median}, \code{cell_diam_median} (twice the median RMS
#'   radius of the cells, in coordinate units), \code{rss}, \code{cp},
#'   \code{moran_i}, \code{moran_z} and \code{reliability}; columns a missing
#'   input leaves undefined are \code{NA}.  Attributes: \code{bounds} (a list
#'   with \code{floor}, \code{ceiling}, \code{supported}, \code{area},
#'   \code{range}, \code{n}, \code{min_cell_n}), \code{variogram} (a list with
#'   \code{nugget}, \code{psill}, \code{range}, \code{model}; \code{NULL}
#'   when none was usable), \code{variable} (\code{"response"},
#'   \code{"residuals"} or \code{NA}), \code{wss_bumps}, \code{nstart},
#'   \code{sac} (the range object used) and, with \code{select_on =
#'   "split"}, \code{split} (a list with \code{selection} and
#'   \code{estimation}, integer row positions in \code{data_sf}).
#' @references
#' Cressie, N. (1996). Change of support and the modifiable areal unit
#' problem. \emph{Geographical Systems}, 3(2--3), 159--180.
#'
#' Fay, R. E. and Herriot, R. A. (1979). Estimates of income for small places:
#' an application of James-Stein procedures to census data. \emph{Journal of
#' the American Statistical Association}, 74(366), 269--277.
#' \doi{10.1080/01621459.1979.10482505}
#'
#' Mallows, C. L. (1973). Some comments on \eqn{C_p}. \emph{Technometrics},
#' 15(4), 661--675. \doi{10.1080/00401706.1973.10489103}
#' @family aggregation
#' @seealso \code{\link{select_resolution}()} to read a level and its flat
#'   region off the profile; \code{\link{determine_optimal_levels}()} for the
#'   integer-vector interface; \code{\link{build_tessellation}()} and
#'   \code{\link{get_voronoi_seeds}()}, which accept the profile or a
#'   \code{select_resolution()} result directly as \code{approx_n_cells}
#'   and \code{n}.
#' @examples
#' if (requireNamespace("gstat", quietly = TRUE)) {
#'   library(sf)
#'   set.seed(2)
#'   n <- 400
#'   xy <- data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000))
#'   D  <- as.matrix(dist(xy))
#'   xy$z <- as.numeric(t(chol(exp(-D / 100) + diag(0.3, n))) %*% rnorm(n))
#'   pts <- st_as_sf(xy, coords = c("x", "y"), crs = 32632)
#'   prof <- resolution_profile(pts, response_var = "z", n_levels = 12)
#'   prof
#'   select_resolution(prof, criterion = "reliability")
#' }
#' @export
resolution_profile <- function(data_sf, response_var = NULL, predictor_vars = NULL,
                               levels = NULL, n_levels = 20L, min_cell_n = 9L,
                               sample_n = 1500L, nstart = 25L, seed = 123L,
                               sac = NULL, quiet = TRUE,
                               select_on = c("all", "split")) {
  .msg <- function(...) if (!quiet) message(...)
  select_on <- match.arg(select_on)
  if (!inherits(data_sf, "sf"))
    stop("resolution_profile(): `data_sf` must be an sf object.", call. = FALSE)
  if (!is.numeric(min_cell_n) || length(min_cell_n) != 1L || !is.finite(min_cell_n) ||
      min_cell_n < 1)
    stop("resolution_profile(): `min_cell_n` must be a single number >= 1.", call. = FALSE)
  if (!is.numeric(nstart) || length(nstart) != 1L || !is.finite(nstart) || nstart < 1)
    stop("resolution_profile(): `nstart` must be a single number >= 1.", call. = FALSE)
  nstart <- as.integer(nstart)
  has_resp <- !is.null(response_var)
  has_pred <- !is.null(predictor_vars) && length(predictor_vars) > 0L
  if (has_pred && !has_resp)
    stop("resolution_profile(): `predictor_vars` needs a `response_var`.", call. = FALSE)

  if (!all(sf::st_geometry_type(data_sf, by_geometry = TRUE) == "POINT"))
    data_sf <- coerce_to_points(data_sf, "auto")
  data_sf <- ensure_projected(data_sf)
  xy <- sf::st_coordinates(data_sf)[, 1:2, drop = FALSE]
  ok_xy <- is.finite(xy[, 1]) & is.finite(xy[, 2])
  if (!all(ok_xy)) {
    .log_warn("resolution_profile(): dropping %d point(s) with empty or non-finite coordinates.",
              sum(!ok_xy))
    data_sf <- data_sf[ok_xy, , drop = FALSE]
    xy <- xy[ok_xy, , drop = FALSE]
  }
  if (nrow(xy) < 3L)
    stop("resolution_profile(): fewer than three points.", call. = FALSE)

  # Sample splitting, on the full layer before the subsample, so the
  # positions index `data_sf` as passed (rows dropped above excepted).
  split <- NULL
  if (identical(select_on, "split")) {
    split <- .spatial_half_split(data_sf, seed = seed, caller = "resolution_profile")
    if (!all(ok_xy)) {
      # Positions refer to the layer after the drop; map them back.
      kept <- which(ok_xy)
      split$selection  <- kept[split$selection]
      split$estimation <- kept[split$estimation]
      keep_sel <- match(split$selection, kept)
    } else keep_sel <- split$selection
    data_sf <- data_sf[keep_sel, , drop = FALSE]
    xy <- xy[keep_sel, , drop = FALSE]
  }

  resp <- NULL; pred <- NULL
  if (has_resp) {
    df <- sf::st_drop_geometry(data_sf)
    if (!(response_var %in% names(df)))
      stop(sprintf("resolution_profile(): column '%s' not found.", response_var), call. = FALSE)
    if (!(is.numeric(df[[response_var]]) || is.logical(df[[response_var]])))
      stop(sprintf("resolution_profile(): `response_var` must be numeric or logical; '%s' is %s.",
                   response_var, class(df[[response_var]])[1L]), call. = FALSE)
    resp <- as.numeric(df[[response_var]])
    if (has_pred) {
      missing_p <- setdiff(predictor_vars, names(df))
      if (length(missing_p))
        stop("resolution_profile(): predictor_vars ",
             paste(sQuote(missing_p), collapse = ", "), " not found.", call. = FALSE)
      non_num <- predictor_vars[!vapply(predictor_vars, function(v)
        is.numeric(df[[v]]) || is.logical(df[[v]]), logical(1))]
      if (length(non_num))
        stop("resolution_profile(): `predictor_vars` must be numeric or logical; ",
             paste(sQuote(non_num), collapse = ", "), " not.", call. = FALSE)
      pred <- as.matrix(df[, predictor_vars, drop = FALSE])
      storage.mode(pred) <- "double"
    }
  }

  cleanup <- .with_seed(seed)
  on.exit(cleanup(), add = TRUE)

  # Subsample, keeping everything aligned.
  n_all <- nrow(xy)
  if (n_all > sample_n) {
    idx <- sample(seq_len(n_all), sample_n)
    xy <- xy[idx, , drop = FALSE]
    data_sf <- data_sf[idx, , drop = FALSE]
    if (has_resp) resp <- resp[idx]
    if (has_pred) pred <- pred[idx, , drop = FALSE]
  }
  n <- nrow(xy)

  # The variable the cells have to represent: the response, or what the
  # predictors leave of it.  Rows with a non-finite value drop out of the RSS
  # and the Moran statistic but stay in the geometry.
  variable <- NA_character_
  y <- NULL
  if (has_resp) {
    y <- resp
    variable <- "response"
    if (has_pred) {
      fit <- try(stats::lm.fit(x = cbind(1, pred), y = resp), silent = TRUE)
      ok_rows <- is.finite(resp) & apply(is.finite(pred), 1L, all)
      if (!inherits(fit, "try-error") && sum(ok_rows) > ncol(pred) + 1L) {
        fit <- stats::lm.fit(x = cbind(1, pred[ok_rows, , drop = FALSE]), y = resp[ok_rows])
        y <- rep(NA_real_, n); y[ok_rows] <- fit$residuals
        variable <- "residuals"
      } else {
        .log_warn("resolution_profile(): the OLS fit on `predictor_vars` failed; scoring the raw response.")
      }
    }
  }

  # Variogram: supplied, or estimated on the subsample.
  vg <- NULL
  if (is.null(sac) && has_resp && requireNamespace("gstat", quietly = TRUE)) {
    sac <- tryCatch(
      suppressWarnings(estimate_sac_range(data_sf, response_var,
                                          predictor_vars = if (has_pred) predictor_vars else NULL,
                                          seed = seed)),
      error = function(e) {
        .log_warn("resolution_profile(): estimate_sac_range() failed (%s); no variogram-based columns.",
                  conditionMessage(e))
        NULL
      })
  }
  range_eff <- NA_real_
  if (!is.null(sac)) {
    vm <- attr(sac, "variogram_model")
    r  <- suppressWarnings(as.numeric(sac))
    if (length(r) == 1L && is.finite(r) && r > 0) range_eff <- r
    cor_fn <- if (is.data.frame(vm)) .vgm_correlation_fn(vm) else NULL
    if (!is.null(cor_fn)) {
      vg <- list(nugget = .vgm_nugget_of(vm),
                 psill  = sum(as.numeric(vm$psill[as.character(vm$model) != "Nug"])),
                 range  = range_eff, model = vm, cor_fn = cor_fn)
      if (!is.finite(vg$nugget)) vg$nugget <- 0
      if (!is.finite(vg$psill) || vg$psill <= 0) vg <- NULL
    }
    if (is.null(vg))
      .log_info("resolution_profile(): the variogram carries no usable model; `cp` and `reliability` are NA.")
  }

  # Bounds.
  hull <- sf::st_convex_hull(sf::st_union(sf::st_geometry(data_sf)))
  area <- suppressWarnings(as.numeric(sf::st_area(hull)))
  bb   <- sf::st_bbox(data_sf)
  bbw  <- as.numeric(bb["xmax"] - bb["xmin"]); bbh <- as.numeric(bb["ymax"] - bb["ymin"])
  if (!is.finite(area) || area <= 0) area <- bbw * bbh
  n_uniq  <- nrow(unique(round(xy, 8)))
  ceiling_L <- max(2L, min(as.integer(floor(n / min_cell_n)), n_uniq - 1L))
  floor_L   <- if (is.finite(range_eff) && range_eff > 0)
    max(2L, as.integer(ceiling(area / range_eff^2))) else 2L
  supported <- floor_L <= ceiling_L
  if (!supported)
    .log_warn(paste0("resolution_profile(): cells no wider than the autocorrelation ",
                     "range (%.0f) would need at least %d of them, but %d points ",
                     "at min_cell_n = %d support at most %d. The data cannot ",
                     "support a tessellation that respects their own ",
                     "correlation structure; the profile runs from 2 to %d so ",
                     "the cost of each level is still visible."),
              range_eff, floor_L, n, as.integer(min_cell_n), ceiling_L, ceiling_L)
  lo <- if (supported) floor_L else 2L
  if (is.null(levels)) {
    levels <- unique(as.integer(round(exp(seq(log(lo), log(ceiling_L),
                                              length.out = max(2L, as.integer(n_levels)))))))
    levels <- sort(unique(c(lo, levels, ceiling_L)))
  } else {
    levels <- sort(unique(as.integer(levels)))
    levels <- levels[levels >= 2L & levels <= n_uniq - 1L]
    if (!length(levels))
      stop("resolution_profile(): no usable value in `levels`.", call. = FALSE)
  }
  bounds <- list(floor = floor_L, ceiling = ceiling_L, supported = supported,
                 area = area, range = range_eff, n = n,
                 min_cell_n = as.integer(min_cell_n))

  rbar_V <- if (!is.null(vg)) .rbar_rect(vg$cor_fn, bbw, bbh) else NA_real_
  pred_for_moran <- if (has_pred) pred else matrix(numeric(0), nrow = n, ncol = 0L)

  # The sweep.
  m <- length(levels)
  out <- data.frame(levels = levels, wss = NA_real_, wss_spread = NA_real_,
                    elbow = NA_real_, cell_n_min = NA_integer_,
                    cell_n_median = NA_real_, cell_diam_median = NA_real_,
                    rss = NA_real_, cp = NA_real_, moran_i = NA_real_,
                    moran_z = NA_real_, reliability = NA_real_)
  for (i in seq_len(m)) {
    L  <- levels[i]
    kb <- .kmeans_best(xy, L, nstart = nstart)
    if (is.null(kb)) {
      .log_warn("resolution_profile(): k-means failed at %d cells; that level is NA.", L)
      next
    }
    km <- kb$km
    out$wss[i]        <- kb$wss
    out$wss_spread[i] <- kb$spread
    sizes <- tabulate(km$cluster, nbins = L)
    out$cell_n_min[i]    <- as.integer(min(sizes))
    out$cell_n_median[i] <- stats::median(sizes)
    rad <- sqrt(km$withinss / pmax(sizes, 1L))
    out$cell_diam_median[i] <- 2 * stats::median(rad[sizes >= 2L])
    if (!is.null(y)) {
      okr <- is.finite(y)
      if (sum(okr) > L) {
        cm  <- stats::ave(y[okr], km$cluster[okr])
        out$rss[i] <- sum((y[okr] - cm)^2)
        if (!is.null(vg))
          out$cp[i] <- out$rss[i] / sum(okr) + 2 * vg$nugget * L / sum(okr)
      }
      mi <- .morans_i_for_k(xy[okr, , drop = FALSE], y[okr],
                            pred_for_moran[okr, , drop = FALSE], km$cluster[okr])
      out$moran_i[i] <- mi[["I"]]
      out$moran_z[i] <- mi[["z"]]
    }
    if (!is.null(vg))
      out$reliability[i] <- .reliability_at(L, area, n, vg$nugget, vg$psill,
                                            vg$cor_fn, rbar_V)
  }

  # Elbow distance over the ladder, and the bumps on it.
  fin <- is.finite(out$wss)
  if (sum(fin) >= 3L) {
    k_norm   <- (levels[fin] - min(levels[fin])) / max(1, diff(range(levels[fin])))
    w        <- out$wss[fin]
    wss_norm <- (w - min(w)) / max(.Machine$double.eps, max(w) - min(w))
    x1 <- k_norm[1L]; y1 <- wss_norm[1L]
    x2 <- k_norm[length(k_norm)]; y2 <- wss_norm[length(wss_norm)]
    line_len <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
    out$elbow[fin] <- if (line_len < .Machine$double.eps) 0 else
      .below_chord(k_norm, wss_norm, x1, y1, x2, y2, line_len)
  }
  wss_bumps <- .wss_bumps(out$wss[fin])
  if (wss_bumps > 0L)
    .log_warn(paste0("resolution_profile(): the WSS curve rises at %d step(s) of ",
                     "the ladder even with %d k-means++ restarts; the elbow and ",
                     "the neighbouring levels' spread are partly optimisation noise."),
              wss_bumps, nstart)

  structure(out,
            class      = c("resolution_profile", "data.frame"),
            bounds     = bounds,
            variogram  = if (is.null(vg)) NULL else vg[c("nugget", "psill", "range", "model")],
            variable   = variable,
            wss_bumps  = wss_bumps,
            nstart     = nstart,
            sac        = sac,
            split      = split)
}


#' Print a resolution profile
#'
#' @param x A \code{resolution_profile}.
#' @param digits Significant digits for the table.  Default 3.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.resolution_profile <- function(x, digits = 3L, ...) {
  b <- attr(x, "bounds")
  # A subset keeps the class and loses the attributes, and `x[, 1:3]` also
  # loses the `levels` column the ladder line is built from.  Print the table
  # as what it now is instead of erroring on the missing pieces; knitr calls
  # print() on a data frame without being asked, so this path is reachable
  # from a document as well as from the console.
  if (is.null(b) || !("levels" %in% names(x))) {
    cat("Resolution profile (subset; the ladder summary is not carried by a",
        "subset)\n\n")
    print(as.data.frame(unclass(x)), row.names = FALSE)
    return(invisible(x))
  }
  cat("Resolution profile:", nrow(x), "levels on", b$n, "points\n")
  cat(sprintf("  ladder      : %d to %d cells (floor %s, ceiling %d at min_cell_n = %d)%s\n",
              min(x$levels), max(x$levels),
              if (is.finite(b$range)) sprintf("%d from range %.0f", b$floor, b$range)
              else "2 (no range)",
              b$ceiling, b$min_cell_n,
              if (isTRUE(b$supported)) "" else "  -- floor above ceiling: not supported"))
  vg <- attr(x, "variogram")
  cat(sprintf("  variogram   : %s\n",
              if (is.null(vg)) "none usable (cp and reliability are NA)" else
                sprintf("nugget %.3g, partial sill %.3g, range %s", vg$nugget, vg$psill,
                        if (is.finite(vg$range)) sprintf("%.0f", vg$range) else "unidentified")))
  cat(sprintf("  scored on   : %s; %d k-means++ restarts per level; WSS rises at %d step(s)\n",
              if (is.na(attr(x, "variable"))) "geometry only" else attr(x, "variable"),
              attr(x, "nstart"), attr(x, "wss_bumps")))
  sp <- attr(x, "split")
  if (!is.null(sp))
    cat(sprintf("  split       : selected on %d points; estimate on the other %d (attr \"split\")\n",
                length(sp$selection), length(sp$estimation)))
  cat("\n")
  tab <- as.data.frame(unclass(x))
  attributes(tab)[setdiff(names(attributes(tab)), c("names", "row.names", "class"))] <- NULL
  num <- vapply(tab, is.double, logical(1))
  tab[num] <- lapply(tab[num], signif, digits = digits)
  print(tab, row.names = FALSE)
  invisible(x)
}


#' Read a level, and the region over which it is not distinguishable, off a profile
#'
#' Picks the level a criterion prefers, together with the \emph{flat region}:
#' every level whose criterion value is within \code{tol} of the optimum.  On
#' the criteria this package computes the flat region is routinely wide ---
#' the reliability curve is flat to within 2 percent over a factor of 3--6 in
#' the number of cells, and \eqn{C_p} on a smooth field descends to the
#' support ceiling --- so the region is the answer, and the argmin only a
#' point in it.  When the optimum sits at the ladder's ceiling or floor the
#' result says so, because a bound is then doing the choosing rather than the
#' criterion (see \code{\link{resolution_profile}} for what each criterion
#' measures and how it behaved on simulated fields).
#'
#' @param profile A \code{resolution_profile}.
#' @param criterion Which column decides: \code{"cp"} (minimised),
#'   \code{"reliability"} (maximised), \code{"elbow"} (maximised) or
#'   \code{"moran_z"} (\eqn{|z|} minimised).
#' @param tol Width of the flat region.  For \code{cp} and \code{reliability}
#'   it is relative to the optimum's value (\code{0.02} keeps levels within 2
#'   percent of it); for \code{elbow} and \code{moran_z}, whose optimum can
#'   be zero, it is relative to the criterion's range over the ladder.
#' @return A list of class \code{resolution_selection} with \code{best} (the
#'   level), \code{flat} (the levels in the flat region, ascending),
#'   \code{criterion}, \code{value} (the optimum), \code{at_ceiling} and
#'   \code{at_floor} (logical: the optimum is the last or first level of the
#'   ladder), \code{n_levels} and \code{values} (the criterion at every
#'   level).
#' @family aggregation
#' @examples
#' if (requireNamespace("gstat", quietly = TRUE)) {
#'   library(sf)
#'   set.seed(2)
#'   n <- 400
#'   xy <- data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000))
#'   D  <- as.matrix(dist(xy))
#'   xy$z <- as.numeric(t(chol(exp(-D / 100) + diag(0.3, n))) %*% rnorm(n))
#'   pts <- st_as_sf(xy, coords = c("x", "y"), crs = 32632)
#'   prof <- resolution_profile(pts, response_var = "z", n_levels = 12)
#'
#'   sel <- select_resolution(prof, criterion = "reliability")
#'   sel                      # the level, and the flat region around it
#'   sel$flat                 # every level within `tol` of the optimum
#'   sel$at_ceiling           # TRUE would mean the ladder, not the criterion, chose
#'
#'   # A different criterion can prefer a different level while agreeing on the
#'   # region: the flat region is the answer, the argmin a point in it.
#'   select_resolution(prof, criterion = "cp")$flat
#' }
#' @export
select_resolution <- function(profile,
                              criterion = c("cp", "reliability", "elbow", "moran_z"),
                              tol = 0.02) {
  if (!inherits(profile, "resolution_profile"))
    stop("select_resolution(): `profile` must come from resolution_profile().",
         call. = FALSE)
  criterion <- match.arg(criterion)
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0)
    stop("select_resolution(): `tol` must be a single non-negative number.", call. = FALSE)
  v <- as.numeric(profile[[criterion]])
  if (criterion == "moran_z") v <- abs(v)
  ok <- is.finite(v)
  if (!any(ok))
    stop(sprintf(paste0("select_resolution(): `%s` is NA at every level%s."),
                 criterion,
                 switch(criterion,
                        cp = " (it needs a response and a usable variogram)",
                        reliability = " (it needs a usable variogram)",
                        moran_z = " (it needs a response and more than nine cells)",
                        "")), call. = FALSE)
  maximise <- criterion %in% c("reliability", "elbow")
  lv <- profile$levels
  opt <- if (maximise) max(v[ok]) else min(v[ok])
  best <- lv[ok][if (maximise) which.max(v[ok]) else which.min(v[ok])]
  flat <- if (criterion %in% c("cp", "reliability")) {
    if (maximise) lv[ok][v[ok] >= opt * (1 - tol)] else lv[ok][v[ok] <= opt * (1 + tol)]
  } else {
    band <- tol * diff(range(v[ok]))
    if (maximise) lv[ok][v[ok] >= opt - band] else lv[ok][v[ok] <= opt + band]
  }
  structure(list(best = as.integer(best), flat = sort(as.integer(flat)),
                 criterion = criterion, value = opt,
                 at_ceiling = best == max(lv), at_floor = best == min(lv),
                 n_levels = length(lv),
                 values = stats::setNames(v, lv)),
            class = "resolution_selection")
}


#' @export
print.resolution_selection <- function(x, ...) {
  cat(sprintf("Resolution by %s: %d cells\n", x$criterion, x$best))
  cat(sprintf("  flat region : %s (%d of %d levels)\n",
              if (length(x$flat) > 1L) sprintf("%d to %d", min(x$flat), max(x$flat))
              else as.character(x$flat),
              length(x$flat), x$n_levels))
  if (isTRUE(x$at_ceiling))
    cat("  note        : the optimum is the support ceiling (n / min_cell_n); the\n",
        "               bound is choosing, not the criterion. Lower min_cell_n to\n",
        "               see whether the criterion keeps descending.\n", sep = "")
  if (isTRUE(x$at_floor))
    cat("  note        : the optimum is the first level of the ladder; the floor\n",
        "               is choosing, not the criterion.\n", sep = "")
  invisible(x)
}


#' Resolve a cell count from a number or from a level-selection result
#'
#' The functions that need a cell count -- \code{build_tessellation()}'s
#' \code{approx_n_cells}, \code{get_voronoi_seeds()}'s \code{n} -- accept,
#' besides a number, the object the level-selection step produced:
#' \code{determine_optimal_levels()}'s integer vector of ranked candidates
#' (the first is used), a \code{resolution_selection} (its \code{$best}), or a
#' \code{resolution_profile} (read with \code{select_resolution()} at its
#' default criterion).  The count and where it came from are returned so the
#' caller can record both on its output.
#'
#' @param x \code{NULL}, a number, a numeric vector, a
#'   \code{resolution_selection} or a \code{resolution_profile}.
#' @param arg,caller Names for the messages.
#' @return \code{NULL} when \code{x} is \code{NULL}; otherwise a list with
#'   \code{n} (a single number) and \code{from} (a character description, or
#'   \code{NULL} when \code{x} was a plain number).
#' @keywords internal
#' @noRd
.resolve_cell_count <- function(x, arg, caller) {
  if (is.null(x)) return(NULL)
  from <- NULL
  if (inherits(x, "resolution_selection")) {
    n <- x$best
    from <- sprintf("select_resolution(criterion = \"%s\")%s", x$criterion,
                    if (isTRUE(x$at_ceiling)) ", an optimum at the support ceiling"
                    else if (isTRUE(x$at_floor)) ", an optimum at the range floor"
                    else "")
  } else if (inherits(x, "resolution_profile")) {
    sel <- select_resolution(x)
    n <- sel$best
    from <- sprintf("resolution_profile() read with select_resolution(criterion = \"%s\")%s",
                    sel$criterion,
                    if (isTRUE(sel$at_ceiling)) ", an optimum at the support ceiling"
                    else if (isTRUE(sel$at_floor)) ", an optimum at the range floor"
                    else "")
    .log_info("%s(): `%s` is a resolution profile; read with select_resolution()'s default criterion (%s): %d cells.",
              caller, arg, sel$criterion, n)
  } else if (is.numeric(x) && length(x) >= 1L && !is.list(x)) {
    n <- x[[1L]]
    if (length(x) > 1L) {
      from <- sprintf("the first of %d ranked candidates (%s)", length(x),
                      paste(format(x, trim = TRUE), collapse = ", "))
      .log_info("%s(): `%s` holds %d candidates, as determine_optimal_levels() returns them; using the first (%s).",
                caller, arg, length(x), format(n))
    }
  } else {
    stop(sprintf(paste0("%s(): `%s` must be a number, the integer vector ",
                        "determine_optimal_levels() returns, a select_resolution() ",
                        "result or a resolution_profile(); got an object of class %s."),
                 caller, arg, paste(class(x), collapse = "/")), call. = FALSE)
  }
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1)
    stop(sprintf("%s(): `%s` must resolve to a positive number of cells; got %s.",
                 caller, arg, paste(format(n), collapse = ", ")), call. = FALSE)
  list(n = as.numeric(n), from = from)
}
