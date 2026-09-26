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
#' this is the quantity Krige's additivity relation needs.  The variance of
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
#' optimum, unlike "minimise the SE of the cell means", which bigger cells
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
#' returns the table, so that the level chosen, by
#' \code{\link{select_resolution}()} or by eye, can be defended with the
#' whole profile in place of one criterion's argmin.  Mallows (1973) presents
#' \eqn{C_p} itself as a display of the bias-variance trade-off across
#' candidates; he does not offer it as a rule that picks one.  This is that
#' display, with the other criteria alongside.
#'
#' @section The ladder and its bounds:
#' Cells are k-means clusters of the (projected) coordinates, fitted at each
#' level as the best of \code{nstart} k-means++ restarts (see
#' \code{\link{determine_optimal_levels}} for why).  Levels are spaced
#' logarithmically, because cell diameter scales as \eqn{L^{-1/2}}: a unit
#' step wastes fits at large \eqn{L} and starves resolution at small.  The
#' ladder runs from a floor to a ceiling the data impose.  The ceiling is
#' \code{floor(n / min_cell_n)}, with \eqn{n} every point of the layer (not
#' the subsample): cells with fewer than \code{min_cell_n} points
#' on average have too little support, and Moran's z is not computable at
#' nine cells or fewer in any case.  The floor is
#' \code{ceiling(area / range^2)} when an autocorrelation range is available:
#' cells wider than the range average over more than one patch of the field.
#' When the floor exceeds the ceiling the data cannot support a tessellation
#' that respects their own correlation structure; that is reported as a
#' finding (a logged warning, and \code{attr(x, "bounds")$supported} is
#' \code{FALSE}) and the ladder runs from 2 to the ceiling anyway, so the
#' profile still shows what each level costs.  On a layer larger than
#' \code{sample_n} the ceiling is also held to half the subsample (two
#' subsample points per cell, the least a fitted cell can be scored on);
#' \code{ceiling_from} is then \code{"sample_n"}, a floor above that is logged
#' with a request to raise \code{sample_n}, and it does not make
#' \code{supported} \code{FALSE}.
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
#'     It estimates the error of predicting a new observation by the mean of
#'     its cell.  When the cells are fitted to a subsample of \eqn{m} of the
#'     \eqn{N} points with a response, the penalty is split between the two:
#'     \eqn{RSS(L)/m + \tau^2 L_m / m + \tau^2 L / N}, where the first two
#'     terms estimate the approximation error from the subsample (adding back
#'     the optimism of its own cell means, over the \eqn{L_m} cells its scored
#'     points fall in) and the last is the variance of cell means built from
#'     all \eqn{N}, which is what the tessellation will carry.  With no
#'     subsample it is the formula above.
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
#'     Calibrated and flat in \eqn{L} on a response with no structure (see
#'     \code{\link{determine_optimal_levels}}), so it separates levels only
#'     where structure remains.  \code{NA} at nine cells or fewer.}
#'   \item{\code{reliability}}{The between-cell signal's share of the spread
#'     in the cell means, from the fitted variogram alone via Krige's
#'     additivity relation (Cressie 1996), for square cells of the level's
#'     average area holding the level's average share of the layer's points
#'     with a response (larger is better).
#'     This is the shrinkage factor of Fay and Herriot (1979).  It has an
#'     interior optimum, and a broad one: validated against the empirical
#'     reliability of true block means on simulated fields, the analytic and
#'     empirical optima agreed to within a level or two where the empirical
#'     estimate was stable, and the band within 2 percent of the maximum
#'     spanned a factor of 3--6 in \eqn{L}.  Read the flat region, not the
#'     argmax.  \code{NA} without a usable variogram.}
#' }
#' \code{cp} and \code{reliability} answer different questions: how well the
#' cells represent the field, and whether the cell values are distinguishable
#' from noise.  They can disagree, and both are shown so the choice between
#' them is made knowingly.
#'
#' @param data_sf An sf object of points (other geometries are reduced to
#'   representative points).
#' @param response_var Optional response column name (numeric or logical).
#'   Enables \code{cp} and \code{moran_z}.  A variogram estimated from it also
#'   sets the floor of the ladder and \code{reliability}.  Rows where it, or a
#'   predictor, is missing or non-finite stay in the geometry and are left out
#'   of the OLS fit, the RSS, \code{cp} and \code{moran_z}; a logged warning
#'   gives their number.
#' @param predictor_vars Optional predictor column names (numeric or
#'   logical).  With them, \code{cp} scores the OLS residuals of the response
#'   on the predictors, the variogram is estimated from those residuals, and
#'   \code{moran_z} regresses the cell means on the cell-mean predictors.
#' @param levels Optional integer vector of level counts to score, replacing
#'   the ladder; values below 2, or at or above the number of distinct
#'   locations, are dropped (k-means cannot place more centres than there are
#'   distinct points).
#' @param n_levels Number of levels on the ladder.  Default 20.
#' @param min_cell_n Minimum average number of points per cell that a level
#'   must keep; sets the ceiling.  Default 9.
#' @param sample_n Points are subsampled to this many before the k-means
#'   fits, as in \code{determine_optimal_levels()}.  Default 1500.  The
#'   columns read off the fitted cells (\code{wss}, the \code{cell_} columns,
#'   \code{rss}, \code{moran_i}, \code{moran_z}) describe the subsample; the
#'   bounds, \code{supported}, the variance term of \code{cp} and
#'   \code{reliability} describe every point of the layer, so the answer does
#'   not change with \code{sample_n} except through the fits.
#' @param nstart k-means++ restarts per level.  Default 25.
#' @param seed RNG seed for the subsample and the restarts; restored
#'   afterwards.  Default 123.
#' @param sac Optional \code{sac_range} object from
#'   \code{\link{estimate_sac_range}()} to take the range, nugget and
#'   correlation function from.  Pass one fitted with \code{detrend =
#'   "reml"}, say, or on a residual field of your choosing.  When
#'   \code{NULL} and a response is given, one is estimated on the subsample
#'   with the same \code{predictor_vars}.
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
#'   with \code{floor}, \code{ceiling}, \code{ceiling_from} (\code{"min_cell_n"},
#'   \code{"distinct locations"} or \code{"sample_n"}, whichever bound it),
#'   \code{supported}, \code{area}, \code{range}, \code{n} (the points in the
#'   layer), \code{n_sample} (the points the k-means fits ran on),
#'   \code{n_distinct}, \code{min_cell_n}), \code{variogram} (a list with
#'   \code{nugget}, \code{psill}, \code{range}, \code{model}; \code{NULL}
#'   when none was usable), \code{variable} (\code{"response"},
#'   \code{"residuals"} or \code{NA}), \code{wss_bumps}, \code{nstart},
#'   \code{sac} (the range object used) and, with \code{select_on =
#'   "split"}, \code{split} (a \code{spatialkit_split}: \code{selection}
#'   and \code{estimation}, integer row positions in \code{data_sf}, with
#'   the \code{method} and \code{seed} that made them).
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
#'   # An exponential field with range parameter 200 (true effective range
#'   # 600 m) on a 1 km square, with a nugget of 0.6 on a unit sill: enough
#'   # noise for Mallows' Cp to have an interior optimum rather than descend
#'   # to the ceiling.
#'   set.seed(2)
#'   n <- 400
#'   xy <- data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000))
#'   D  <- as.matrix(dist(xy))
#'   xy$z <- as.numeric(t(chol(exp(-D / 200) + diag(0.6, n))) %*% rnorm(n))
#'   pts <- st_as_sf(xy, coords = c("x", "y"), crs = 32632)
#'   prof <- resolution_profile(pts, response_var = "z", n_levels = 12)
#'   print(prof)               # one row per level; print() because only the
#'                             # last value of a braced block is shown
#'   select_resolution(prof, criterion = "cp")
#' }
#' @export
resolution_profile <- function(data_sf, response_var = NULL, predictor_vars = NULL,
                               levels = NULL, n_levels = 20L, min_cell_n = 9L,
                               sample_n = 1500L, nstart = 25L, seed = 123L,
                               sac = NULL, select_on = c("all", "split")) {
  select_on <- match.arg(select_on)
  if (!inherits(data_sf, "sf"))
    stop("resolution_profile(): `data_sf` must be an sf object.", call. = FALSE)
  if (!is.numeric(min_cell_n) || length(min_cell_n) != 1L || !is.finite(min_cell_n) ||
      min_cell_n < 1)
    stop("resolution_profile(): `min_cell_n` must be a single number >= 1.", call. = FALSE)
  if (!is.numeric(nstart) || length(nstart) != 1L || !is.finite(nstart) || nstart < 1)
    stop("resolution_profile(): `nstart` must be a single number >= 1.", call. = FALSE)
  nstart <- as.integer(nstart)
  # These three used to fail deep inside seq() or an if() with R's own message,
  # and a fractional min_cell_n printed a ceiling that did not match it.
  if (!is.numeric(n_levels) || length(n_levels) != 1L || !is.finite(n_levels) || n_levels < 2)
    stop("resolution_profile(): `n_levels` must be a single number >= 2.", call. = FALSE)
  if (!is.numeric(sample_n) || length(sample_n) != 1L || !is.finite(sample_n) || sample_n < 3)
    stop("resolution_profile(): `sample_n` must be a single number >= 3.", call. = FALSE)
  min_cell_n <- as.integer(floor(min_cell_n))
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

  # Rows the response criteria can read: a finite response and, with
  # predictors, finite predictors.  The rest stay in the geometry.
  resp_ok <- NULL
  if (has_resp) {
    resp_ok <- is.finite(resp)
    if (has_pred) resp_ok <- resp_ok & apply(is.finite(pred), 1L, all)
    if (!all(resp_ok))
      .log_warn(paste0("resolution_profile(): %d of %d row(s) have a missing or ",
                       "non-finite %s; they stay in the geometry but are left out ",
                       "of the OLS fit, the RSS, Cp and Moran's z."),
                sum(!resp_ok), length(resp_ok),
                if (has_pred) "response or predictor" else "response")
  }

  # The bounds describe the layer the cells will be built on and aggregated
  # over, so they are read off every point: the extent, the distinct
  # locations and the number of points that carry a response.  Only the
  # k-means work below runs on the subsample.  The support ceiling, the
  # verdict, Cp's penalty and the reliability used to take n from the
  # subsample, which capped every layer larger than sample_n at
  # floor(sample_n / min_cell_n) cells and made the answer depend on sample_n.
  n_all  <- nrow(xy)
  n_resp <- if (has_resp) sum(resp_ok) else n_all
  hull <- sf::st_convex_hull(sf::st_union(sf::st_geometry(data_sf)))
  area <- suppressWarnings(as.numeric(sf::st_area(hull)))
  bb   <- sf::st_bbox(data_sf)
  bbw  <- as.numeric(bb["xmax"] - bb["xmin"]); bbh <- as.numeric(bb["ymax"] - bb["ymin"])
  if (!is.finite(area) || area <= 0) area <- bbw * bbh
  # duplicated() on a complex vector, because unique() on an n x 2 matrix
  # pastes every row into a string: 0.01 s against 5 s at a million points.
  n_uniq_all <- sum(!duplicated(complex(real = round(xy[, 1], 8),
                                        imaginary = round(xy[, 2], 8))))

  # Subsample, keeping everything aligned.
  if (n_all > sample_n) {
    idx <- sample(seq_len(n_all), sample_n)
    xy <- xy[idx, , drop = FALSE]
    data_sf <- data_sf[idx, , drop = FALSE]
    if (has_resp) { resp <- resp[idx]; resp_ok <- resp_ok[idx] }
    if (has_pred) pred <- pred[idx, , drop = FALSE]
  }
  n <- nrow(xy)

  # The variable the cells have to represent: the response, or what the
  # predictors leave of it.  Rows with a non-finite value drop out of the RSS
  # and the Moran statistic but stay in the geometry.
  variable <- NA_character_
  y <- NULL
  if (has_resp) {
    y <- rep(NA_real_, n)
    y[resp_ok] <- resp[resp_ok]
    variable <- "response"
    if (has_pred) {
      # Fitted on the complete rows only.  lm.fit() refuses any NA, and a
      # first fit on every row turned one missing value into "the OLS fit
      # failed": the raw response, trend and all, was then scored against a
      # variogram estimate_sac_range() had fitted to the residuals.
      fit <- if (sum(resp_ok) > ncol(pred) + 1L)
        try(stats::lm.fit(x = cbind(1, pred[resp_ok, , drop = FALSE]), y = resp[resp_ok]),
            silent = TRUE)
      if (!is.null(fit) && !inherits(fit, "try-error")) {
        y[resp_ok] <- fit$residuals
        variable <- "residuals"
      } else {
        .log_warn(paste0("resolution_profile(): the OLS fit on `predictor_vars` failed ",
                         "on the %d complete row(s); scoring the raw response."),
                  sum(resp_ok))
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

  # Bounds.  `area`, `n_uniq_all` and `n_all` are the layer's (see above);
  # `n_uniq` is the subsample's, the locations k-means can place centres on.
  n_uniq  <- nrow(unique(round(xy, 8)))
  if (n == n_all) n_uniq_all <- n_uniq
  by_support <- floor(n_all / min_cell_n)
  # The ceiling is the smallest of what the layer's point count supports,
  # what its distinct locations allow (k-means cannot place more centres than
  # there are distinct points) and, when the fits run on a subsample, what the
  # subsample can fit: two of its points per cell, and one short of its
  # distinct locations.  Which one bound it is recorded, because the print and
  # the "bound is choosing" notes name it.
  data_ceiling <- min(by_support, n_uniq_all - 1L)
  fit_cap      <- if (n < n_all) min(floor(n / 2), n_uniq - 1L) else Inf
  ceiling_L    <- max(2L, as.integer(min(data_ceiling, fit_cap)))
  ceiling_from <- if (fit_cap < data_ceiling) "sample_n"
                  else if (by_support <= n_uniq_all - 1L) "min_cell_n"
                  else "distinct locations"
  # Kept as a double until the comparison: a short range on a continental
  # extent puts area / range^2 past .Machine$integer.max, and as.integer() of
  # that is NA, which turned the "not supported" branch into an abort.
  floor_raw <- if (is.finite(range_eff) && range_eff > 0)
    max(2, ceiling(area / range_eff^2)) else 2
  # The verdict is about the data, so it is taken against the layer's
  # ceiling; a subsample too small to reach the floor is a setting to change,
  # and is said separately.
  supported <- floor_raw <= max(2, data_ceiling)
  floor_L   <- if (floor_raw <= .Machine$integer.max) as.integer(floor_raw) else NA_integer_
  if (!supported)
    .log_warn(paste0("resolution_profile(): cells no wider than the autocorrelation ",
                     "range (%.0f) would need at least %.0f of them, but %d points ",
                     "at min_cell_n = %d support at most %d. The data cannot ",
                     "support a tessellation that respects their own ",
                     "correlation structure; the profile runs from 2 to %d so ",
                     "the cost of each level is still visible."),
              range_eff, floor_raw, n_all, min_cell_n,
              max(2L, as.integer(data_ceiling)), ceiling_L)
  else if (floor_raw > ceiling_L)
    .log_warn(paste0("resolution_profile(): cells no wider than the autocorrelation ",
                     "range (%.0f) need at least %.0f of them, which the %d points ",
                     "support, but k-means on the %d-point subsample can fit at ",
                     "most %d; raise `sample_n`. The profile runs from 2 to %d."),
              range_eff, floor_raw, n_all, n, ceiling_L, ceiling_L)
  lo <- if (supported && floor_raw <= ceiling_L) floor_L else 2L
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
  bounds <- list(floor = floor_L, ceiling = ceiling_L, ceiling_from = ceiling_from,
                 supported = supported, area = area, range = range_eff, n = n_all,
                 n_sample = n, n_distinct = n_uniq_all, min_cell_n = min_cell_n)

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
      m_ok <- sum(okr)
      L_ok <- length(unique(km$cluster[okr]))     # cells holding a scored row
      if (m_ok > L_ok) {
        cm  <- stats::ave(y[okr], km$cluster[okr])
        out$rss[i] <- sum((y[okr] - cm)^2)
        # Cp for the cells of the whole layer, estimated from the subsample.
        # RSS / m + tau^2 L_ok / m estimates the approximation error plus
        # tau^2 (it adds back the optimism of the subsample's own cell means,
        # tau^2 per cell); tau^2 L / N is the sampling variance of cell means
        # built from the layer's N rows with a response.  With no subsample
        # (m = N, L_ok = L) it is Mallows' RSS / n + 2 tau^2 L / n.
        if (!is.null(vg))
          out$cp[i] <- out$rss[i] / m_ok + vg$nugget * (L_ok / m_ok + L / n_resp)
      }
      mi <- .morans_i_for_k(xy[okr, , drop = FALSE], y[okr],
                            pred_for_moran[okr, , drop = FALSE], km$cluster[okr])
      out$moran_i[i] <- mi[["I"]]
      out$moran_z[i] <- mi[["z"]]
    }
    # The layer's point count: the cell means are built from every point,
    # not from the subsample the partition was fitted on.
    if (!is.null(vg))
      out$reliability[i] <- .reliability_at(L, area, n_resp, vg$nugget, vg$psill,
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
#' @family print methods
#' @export
print.resolution_profile <- function(x, digits = 3L, ...) {
  b <- attr(x, "bounds")
  # `[` keeps the class here.  A column subset such as `x[, 1:3]` loses the
  # attributes and can lose the `levels` column the ladder line is built from;
  # a row subset keeps both, and prints honestly below (its ladder line is the
  # rows it has, its floor and ceiling are the run's).  The one row subset that
  # cannot print is the empty one: min() and max() of no levels are Inf and
  # -Inf, and `%d` refuses those.  knitr calls print() on a data frame without
  # being asked, so every one of these paths is reachable from a document.
  if (nrow(x) == 0L) {
    cat("Resolution profile: 0 levels (an empty subset)\n")
    return(invisible(x))
  }
  if (is.null(b) || !("levels" %in% names(x))) {
    cat("Resolution profile (subset; the ladder summary is not carried by a",
        "subset)\n\n")
    print(as.data.frame(unclass(x)), row.names = FALSE)
    return(invisible(x))
  }
  # `n` is the layer; the k-means fits may have run on a subsample of it.
  n_fit <- b$n_sample %||% b$n
  cat("Resolution profile:", nrow(x), "levels on", b$n,
      if (isTRUE(n_fit < b$n)) sprintf("points (k-means fitted to a subsample of %d)\n", n_fit)
      else "points\n")
  cat(sprintf("  ladder      : %d to %d cells (floor %s, ceiling %d from %s)%s\n",
              min(x$levels), max(x$levels),
              if (!is.finite(b$range)) "2 (no range)"
              else if (is.na(b$floor)) sprintf("beyond integer range from range %.0f", b$range)
              else sprintf("%d from range %.0f", b$floor, b$range),
              b$ceiling,
              if (identical(b$ceiling_from, "distinct locations"))
                sprintf("%d distinct locations", b$n_distinct %||% NA_integer_)
              else if (identical(b$ceiling_from, "sample_n"))
                sprintf("the %d-point subsample", n_fit)
              else sprintf("min_cell_n = %d", b$min_cell_n),
              if (!isTRUE(b$supported)) "  -- floor above ceiling: not supported"
              else if (isTRUE(b$floor > b$ceiling)) "  -- floor above ceiling: raise sample_n"
              else ""))
  vg <- attr(x, "variogram", exact = TRUE)
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


# A flat region is every level within `tol` of the optimum, and the criterion
# curves are not monotone, so the region is a SET and routinely has holes in
# it: over 25 simulated fields, 9 of 100 criterion bands skipped at least one
# rung, one of them printing as "12 to 19" while rejecting five of the seven
# levels inside that range.  Adjacency is a fact about the ladder rather than
# about the numbers (19 and 21 are neighbours when nothing was scored between
# them), so the ladder has to come in with the band.
# The levels a criterion actually scored, read off the named value vector
# select_resolution() carries.
.scored_levels <- function(values) {
  if (is.null(values)) return(integer(0))
  as.integer(names(values))[is.finite(values)]
}

.band_runs <- function(flat, ladder) {
  flat <- sort(unique(as.integer(flat)))
  if (!length(flat)) return(NULL)
  # `ladder` is the levels the criterion was SCORED at, not every level on the
  # profile: a level k-means failed at is NA for every criterion, and treating
  # it as a rejection would print a hole the criterion never made.
  # Sorted here, not trusted: a profile re-sorted by a criterion column
  # (`prof[order(prof$cp), ]`) keeps its class and hands its levels over in
  # that order, and adjacency read off an unsorted ladder split every band.
  ladder <- sort(unique(as.integer(ladder)))
  pos <- match(flat, ladder)
  if (anyNA(pos)) return(data.frame(from = min(flat), to = max(flat)))
  cuts <- c(0L, which(diff(pos) != 1L), length(pos))
  do.call(rbind, lapply(seq_len(length(cuts) - 1L), function(i) {
    r <- flat[(cuts[i] + 1L):cuts[i + 1L]]
    data.frame(from = min(r), to = max(r))
  }))
}

# A band is a set of levels on a discrete ladder, drawn on a continuous axis.
# A run of a single level has zero width and vanishes, so each run is widened
# to the half-gap either side: the shading then covers the rungs it names and
# stops halfway to the ones it does not.  The end rungs mirror their one
# neighbouring gap so they get a band the same shape as the rest.
.band_rects <- function(flat, scored, ladder = scored, log_x = FALSE) {
  # Two ladders: adjacency (whether a gap is a rejection) is read off the
  # levels the criterion was SCORED at, while the half-gaps are measured on
  # the full ladder the axis shows, so an unscored rung between two accepted
  # ones is not shaded over and an end rung mirrors its real neighbour.
  runs <- .band_runs(flat, scored)
  if (is.null(runs)) return(NULL)
  lv <- sort(unique(as.numeric(ladder)))
  plain <- data.frame(xmin = runs$from, xmax = runs$to)
  if (length(lv) < 2L) return(plain)
  mid <- if (isTRUE(log_x)) function(a, b) sqrt(a * b) else function(a, b) (a + b) / 2
  inner <- mid(lv[-length(lv)], lv[-1L])
  edges <- c(if (isTRUE(log_x)) lv[1L]^2 / inner[1L] else 2 * lv[1L] - inner[1L],
             inner,
             if (isTRUE(log_x)) lv[length(lv)]^2 / inner[length(inner)]
             else 2 * lv[length(lv)] - inner[length(inner)])
  i <- match(runs$from, lv); j <- match(runs$to, lv)
  if (anyNA(i) || anyNA(j) || !all(is.finite(edges))) return(plain)
  data.frame(xmin = edges[i], xmax = edges[j + 1L])
}


.band_label <- function(flat, ladder) {
  runs <- .band_runs(flat, ladder)
  if (is.null(runs)) return("none")
  paste(ifelse(runs$from == runs$to, as.character(runs$from),
               sprintf("%d to %d", runs$from, runs$to)), collapse = ", ")
}


#' Read a level, and the region over which it is not distinguishable, off a profile
#'
#' Picks the level a criterion prefers, together with the \emph{flat region}:
#' every level whose criterion value is within \code{tol} of the optimum.  On
#' the criteria this package computes the flat region is routinely wide: the
#' reliability curve is flat to within 2 percent over a factor of 3--6 in the
#' number of cells, and \eqn{C_p} on a smooth field descends to the support
#' ceiling.  The region is the answer, and the argmin only a point in it.
#' When the optimum sits at an end of the levels the criterion was scored at,
#' the result says so and names the bound, because a bound is then doing the
#' choosing rather than the criterion (see
#' \code{\link{resolution_profile}} for what each criterion measures and how
#' it behaved on simulated fields).
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
#'   level), \code{flat} (the levels in the flat region, ascending; a set,
#'   which can skip a rung), \code{criterion}, \code{value} (the optimum),
#'   \code{at_ceiling} and \code{at_floor} (logical: the optimum is the last
#'   or first of the levels this criterion was scored at, which for
#'   \code{moran_z} starts above nine cells), \code{edge} (which bound that
#'   is, in words: the support ceiling, the subsample's ceiling, the range
#'   floor, the ladder's own end,
#'   or the first or last level the criterion is computable at; \code{NA} for
#'   an interior optimum), \code{n_levels} and \code{values} (the criterion
#'   at every level, \code{NA} where it could not be computed).
#' @family aggregation
#' @examples
#' if (requireNamespace("gstat", quietly = TRUE)) {
#'   library(sf)
#'   # An exponential field with range parameter 200 (true effective range
#'   # 600 m) on a 1 km square, with a nugget of 0.6 on a unit sill: enough
#'   # noise for Mallows' Cp to have an interior optimum rather than descend
#'   # to the ceiling.
#'   set.seed(2)
#'   n <- 400
#'   xy <- data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000))
#'   D  <- as.matrix(dist(xy))
#'   xy$z <- as.numeric(t(chol(exp(-D / 200) + diag(0.6, n))) %*% rnorm(n))
#'   pts <- st_as_sf(xy, coords = c("x", "y"), crs = 32632)
#'   prof <- resolution_profile(pts, response_var = "z", n_levels = 12)
#'
#'   sel <- select_resolution(prof, criterion = "cp")
#'   print(sel)               # the level, and the flat region around it
#'   print(sel$flat)          # every level within `tol` of the optimum
#'   print(sel$edge)          # NA here: the optimum is interior
#'
#'   # Reliability prefers coarse cells and here runs into the floor the
#'   # autocorrelation range sets, which the result says in words.
#'   rel <- select_resolution(prof, criterion = "reliability")
#'   print(rel)
#'   c(at_floor = rel$at_floor, edge = rel$edge)
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
  # `[` keeps the class, so a subset arrives here looking like a profile.  An
  # empty one or a column subset missing the criterion used to fall through to
  # the "NA at every level" error below, which sends the user off to add a
  # response or fix a variogram when the real problem is the object.
  if (nrow(profile) == 0L)
    stop("select_resolution(): `profile` has no levels; it is an empty subset ",
         "of a resolution profile.", call. = FALSE)
  if (!criterion %in% names(profile) || !"levels" %in% names(profile))
    stop(sprintf(paste0("select_resolution(): `profile` has no `%s` column%s; ",
                        "it is a column subset of a resolution profile. Pass ",
                        "the profile resolution_profile() returned."),
                 if (criterion %in% names(profile)) "levels" else criterion,
                 if (!criterion %in% names(profile) && !"levels" %in% names(profile))
                   " (nor `levels`)" else ""),
         call. = FALSE)
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
  # A negative optimum sends the multiplicative threshold PAST the optimum, so
  # not even the argmin qualifies and the band comes back empty.  The optimum
  # is within any tolerance of itself by definition, so that is the floor.
  if (!length(flat)) flat <- best
  # The edge flags are relative to the levels this criterion was SCORED at,
  # not to the ladder: Moran's z is NA at nine cells or fewer, and an optimum
  # at the first level it could be computed at is a bound choosing just as
  # much as the ladder's end is.  `edge` then says WHICH bound, read off the
  # profile's own record of why the ladder stops where it does, because the
  # notes used to name "the support ceiling (n / min_cell_n)" for a ladder
  # the caller had set with `levels =`.
  scored <- lv[ok]
  at_ceiling <- best == max(scored); at_floor <- best == min(scored)
  edge <- .ladder_edge(best, scored, lv, attr(profile, "bounds"))
  structure(list(best = as.integer(best), flat = sort(as.integer(flat)),
                 criterion = criterion, value = opt,
                 at_ceiling = at_ceiling, at_floor = at_floor, edge = edge,
                 n_levels = length(lv),
                 values = stats::setNames(v, lv)),
            class = "resolution_selection")
}


# Which bound an optimum at the end of the scored levels is sitting on.  NA
# when it is interior.  The wording is shared by the print methods, the
# profile plot's caption and the tessellation builders' provenance note.
.ladder_edge <- function(best, scored, ladder, bounds) {
  if (best == max(scored) && best == min(scored)) return("the only level scored")
  if (best == max(scored)) {
    if (max(scored) < max(ladder))
      return("the last level the criterion is computable at")
    if (is.list(bounds) && identical(as.integer(max(ladder)), as.integer(bounds$ceiling)))
      return(if (identical(bounds$ceiling_from, "distinct locations"))
               "the support ceiling (one short of the distinct locations)"
             else if (identical(bounds$ceiling_from, "sample_n"))
               "the subsample's ceiling (two subsample points per cell; raise sample_n)"
             else "the support ceiling (n / min_cell_n)")
    return("the last level of the ladder")
  }
  if (best == min(scored)) {
    if (min(scored) > min(ladder))
      return("the first level the criterion is computable at")
    if (is.list(bounds) && isTRUE(bounds$supported) && is.finite(bounds$range %||% NA) &&
        identical(as.integer(min(ladder)), as.integer(bounds$floor)))
      return("the range floor (area / range^2)")
    return("the first level of the ladder")
  }
  NA_character_
}


#' @export
print.resolution_selection <- function(x, ...) {
  cat(sprintf("Resolution by %s: %d cells\n", x$criterion, x$best))
  cat(sprintf("  flat region : %s (%d of %d levels)\n",
              .band_label(x$flat, .scored_levels(x$values)),
              length(x$flat), x$n_levels))
  edge <- x$edge %||% NA_character_
  if (!is.na(edge)) {
    hint <- if (grepl("min_cell_n", edge, fixed = TRUE))
      " Lower min_cell_n to see whether the criterion keeps going."
    else if (grepl("raise sample_n", edge, fixed = TRUE))
      " Raise sample_n to see whether the criterion keeps going."
    else if (grepl("range floor", edge, fixed = TRUE))
      " Fewer cells would be wider than the range and average over more than one patch of the field."
    else ""
    # strwrap() collapses runs of spaces, so the aligned label is put back
    # after wrapping rather than wrapped with the sentence.
    body <- strwrap(sprintf("the optimum is %s; the bound is choosing, not the criterion.%s",
                            edge, hint), width = 62)
    cat(paste0(c("  note        : ", rep("                ", length(body) - 1L)), body),
        sep = "\n")
  }
  invisible(x)
}


#' Every criterion's pick, side by side
#'
#' \code{\link{select_resolution}()} reads one criterion at a time.  This puts
#' all of them in one table: the level each prefers, the flat region around
#' it, and whether a ladder bound is doing the choosing rather than the
#' criterion.  The closing line gives the levels that lie in \emph{every} flat
#' region, the cell counts no criterion objects to.
#'
#' A flat region is a set, not an interval.  The criterion curves are not
#' monotone, so a region can skip a rung of the ladder, and the table prints
#' what the criterion actually accepts (\code{"26, 31"}, not \code{"26 to
#' 31"}) rather than a range that would quietly include the levels it
#' rejected.
#'
#' That set is often empty, and an empty one is a result rather than a
#' failure.  The criteria answer different questions: how well the cells
#' represent the field (\eqn{C_p}), whether the cell values are
#' distinguishable from noise (\code{reliability}), where the within-cluster
#' sum of squares bends (\code{elbow}), and whether the cell means still carry
#' autocorrelation (\code{moran_z}).  A field with no single right resolution
#' shows up here as disjoint bands, and the spread between the picks is
#' printed for the same reason.
#'
#' Nothing in the table is a decision procedure.  Each flat region is
#' routinely wide, the choice within it belongs to the analyst, and
#' \code{\link{plot.resolution_profile}()} draws the curves the bands were
#' read from.
#'
#' @param object A \code{\link{resolution_profile}()}.
#' @param criteria Character vector, any of \code{"cp"},
#'   \code{"reliability"}, \code{"elbow"}, \code{"moran_z"}.  Default: every
#'   one of the four that is finite at some level.
#' @param tol Passed to \code{\link{select_resolution}()} for the flat
#'   region.  Default 0.02.
#' @param ... Ignored.
#' @return A data.frame of class \code{resolution_summary}, one row per
#'   criterion in the order given, with columns \code{criterion},
#'   \code{best}, \code{flat_min}, \code{flat_max}, \code{n_flat},
#'   \code{value} (the optimum itself, on that criterion's own scale and so
#'   not comparable across rows, which is why the print method leaves it out),
#'   \code{at_floor}, \code{at_ceiling} and \code{edge} (which bound an
#'   edge optimum sits on, as \code{\link{select_resolution}()} reports it;
#'   \code{NA} when interior).  Attributes: \code{bands} (a named list holding each
#'   criterion's flat region in full, since \code{flat_min} and
#'   \code{flat_max} are only its ends and the region can have holes in it),
#'   \code{common} (the levels in every flat region, an integer vector that is
#'   empty when the regions do not overlap), \code{scored} (the levels each
#'   criterion returned a finite value at, which is what decides whether a gap
#'   in a band is a rejection or a level nothing was computed at),
#'   \code{ladder} (every level on the profile), \code{n_levels}, \code{tol}
#'   and \code{variable} (what the criteria were scored on).  The print method recomputes the closing
#'   comparison from \code{bands}, so a row subset of the result reads
#'   honestly.
#' @family aggregation
#' @seealso \code{\link{select_resolution}()} for one criterion, with the
#'   per-level values attached; \code{\link{plot.resolution_profile}()} for
#'   the curves behind these bands.
#' @examples
#' if (requireNamespace("gstat", quietly = TRUE)) {
#'   library(sf)
#'   # An exponential field with range parameter 200 (true effective range
#'   # 600 m) on a 1 km square, with a nugget of 0.6 on a unit sill: enough
#'   # noise for Mallows' Cp to have an interior optimum rather than descend
#'   # to the ceiling.
#'   set.seed(2)
#'   n <- 400
#'   xy <- data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000))
#'   D  <- as.matrix(dist(xy))
#'   xy$z <- as.numeric(t(chol(exp(-D / 200) + diag(0.6, n))) %*% rnorm(n))
#'   pts <- st_as_sf(xy, coords = c("x", "y"), crs = 32632)
#'   prof <- resolution_profile(pts, response_var = "z", n_levels = 12)
#'
#'   # print() is explicit because only the last value of a braced block is
#'   # shown, and the table is the thing worth seeing here.
#'   print(summary(prof))
#'   print(attr(summary(prof), "common")) # the levels all of them accept, if any
#'   attr(summary(prof), "bands")         # each criterion's region in full
#' }
#' @export
summary.resolution_profile <- function(object, criteria = NULL, tol = 0.02, ...) {
  if (!inherits(object, "resolution_profile"))
    stop("summary.resolution_profile(): `object` must come from resolution_profile().",
         call. = FALSE)
  # `[` keeps the class, so a subset arrives here looking like a profile.
  if (nrow(object) == 0L || !("levels" %in% names(object)))
    stop("summary.resolution_profile(): `object` has no levels to summarise; it is an ",
         if (nrow(object) == 0L) "empty" else "incomplete",
         " subset of a resolution profile.", call. = FALSE)
  # `tol` is checked here as well as in select_resolution(), because the try()
  # below would otherwise swallow that check and report a bad tolerance as a
  # fault in the profile's data.
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0)
    stop("summary.resolution_profile(): `tol` must be a single non-negative number.",
         call. = FALSE)
  selectable <- c("cp", "reliability", "elbow", "moran_z")
  # A criterion the caller named is reported on by select_resolution() itself,
  # whose message says which input it is missing.  Only the defaulted set is
  # filtered here, and then a failure is unexpected rather than informative.
  explicit <- !is.null(criteria)
  if (explicit) {
    # setdiff() coerces, so a factor used to pass the check below and then die
    # inside match.arg() with no mention of this function.
    if (!is.character(criteria))
      stop("summary.resolution_profile(): `criteria` must be a character vector; got ",
           class(criteria)[1L], ".", call. = FALSE)
    if (!length(criteria))
      stop("summary.resolution_profile(): `criteria` is empty. Name at least one of ",
           paste(selectable, collapse = ", "), ", or leave it NULL for every one ",
           "the profile can score.", call. = FALSE)
    if (anyDuplicated(criteria))
      stop("summary.resolution_profile(): `criteria` repeats ",
           paste(unique(criteria[duplicated(criteria)]), collapse = ", "),
           ". A criterion counted twice agrees with itself, which would make the ",
           "closing comparison meaningless.", call. = FALSE)
  } else {
    criteria <- selectable[vapply(selectable, function(cn)
      cn %in% names(object) &&
        any(is.finite(suppressWarnings(as.numeric(object[[cn]])))), logical(1))]
  }
  bad <- setdiff(criteria, selectable)
  if (length(bad))
    stop("summary.resolution_profile(): unknown criteria: ", paste(bad, collapse = ", "),
         ". Choose from ", paste(selectable, collapse = ", "), ".", call. = FALSE)
  if (!length(criteria)) {
    # Two different faults reach here, and the wrong message sends the caller
    # off to fix the data when the object is what is wrong.
    if (!any(selectable %in% names(object)))
      stop("summary.resolution_profile(): `object` carries none of the ",
           "criterion columns (", paste(selectable, collapse = ", "),
           "); it is a column subset of a resolution profile. Pass the ",
           "profile resolution_profile() returned.", call. = FALSE)
    stop("summary.resolution_profile(): no criterion is finite at any level. ",
         "cp and reliability need a usable variogram, moran_z a response and ",
         "more than nine cells; a geometry-only profile carries elbow alone.",
         call. = FALSE)
  }

  sels <- if (explicit)
    lapply(criteria, function(cn) select_resolution(object, criterion = cn, tol = tol))
  else {
    got <- lapply(criteria, function(cn)
      try(select_resolution(object, criterion = cn, tol = tol), silent = TRUE))
    ok <- !vapply(got, inherits, logical(1), "try-error")
    if (!any(ok))
      stop("summary.resolution_profile(): no criterion could be read off this profile.",
           call. = FALSE)
    got[ok]
  }

  out <- data.frame(
    criterion  = vapply(sels, function(s) s$criterion, character(1)),
    best       = vapply(sels, function(s) as.integer(s$best), integer(1)),
    flat_min   = vapply(sels, function(s) min(s$flat), integer(1)),
    flat_max   = vapply(sels, function(s) max(s$flat), integer(1)),
    n_flat     = vapply(sels, function(s) length(s$flat), integer(1)),
    value      = vapply(sels, function(s) as.numeric(s$value), numeric(1)),
    at_floor   = vapply(sels, function(s) isTRUE(s$at_floor), logical(1)),
    at_ceiling = vapply(sels, function(s) isTRUE(s$at_ceiling), logical(1)),
    edge       = vapply(sels, function(s) s$edge %||% NA_character_, character(1)),
    stringsAsFactors = FALSE)

  # The levels no criterion objects to.  Intersecting the regions rather than
  # the picks is the point: a single level chosen by one criterion says
  # nothing about whether the others can live with it.
  cnames <- vapply(sels, function(s) s$criterion, character(1))
  bands  <- stats::setNames(lapply(sels, function(s) s$flat), cnames)
  scored <- stats::setNames(lapply(sels, function(s) .scored_levels(s$values)), cnames)
  common <- Reduce(intersect, bands)
  lv <- as.integer(object$levels)
  # The whole ladder, not its ends: the print method needs it to tell a band
  # with a hole in it from a solid run, and `[` on the result keeps it, so a
  # subset can still be read honestly.
  structure(out, class = c("resolution_summary", "data.frame"),
            common = sort(unique(as.integer(common))), bands = bands,
            scored = scored, ladder = lv,
            n_levels = nrow(object), tol = tol, variable = attr(object, "variable"))
}


#' @export
print.resolution_summary <- function(x, ...) {
  need <- c("criterion", "best", "flat_min", "flat_max", "n_flat",
            "at_floor", "at_ceiling")
  # `[` keeps the class, so a column subset arrives here looking like a
  # summary.  It drops every custom attribute while keeping the columns, so
  # BOTH have to be checked: a subset that merely drops `value` keeps all of
  # the columns below and then hands min()/max() no levels, which are Inf and
  # -Inf, which "%d" refuses.  knitr prints a data frame without being asked,
  # so this path is reachable from a document.
  if (!all(need %in% names(x)) || !length(attr(x, "ladder"))) {
    cat("Resolution picks (subset; the comparison is not carried by a column subset)\n\n")
    print(as.data.frame(unclass(x)), row.names = FALSE)
    return(invisible(x))
  }
  lad <- attr(x, "ladder"); bands <- attr(x, "bands"); scored <- attr(x, "scored")
  n_lv <- attr(x, "n_levels") %||% length(lad)
  rungs_of <- function(cn) if (!is.null(scored[[cn]])) scored[[cn]] else lad
  band_of <- function(i) {
    cn <- x$criterion[i]
    b <- if (!is.null(bands[[cn]])) bands[[cn]] else seq(x$flat_min[i], x$flat_max[i])
    .band_label(b, rungs_of(cn))
  }
  cat(sprintf("Resolution picks: %d criteri%s over %d level%s (%d to %d cells)\n\n",
              nrow(x), if (nrow(x) == 1L) "on" else "a",
              n_lv, if (n_lv == 1L) "" else "s",
              min(lad), max(lad)))
  tab <- data.frame(
    criterion = x$criterion,
    best = x$best,
    `flat region` = vapply(seq_len(nrow(x)), band_of, character(1)),
    `levels in band` = x$n_flat,
    check.names = FALSE, stringsAsFactors = FALSE)
  print(tab, row.names = FALSE)
  # The notes go under the table rather than in a column of it: a listed band
  # such as "19 to 21, 26, 31 to 33" is wide, and a fifth column of sentences
  # pushed the whole thing past eighty characters and wrapped it.
  # One line per distinct bound, naming the criteria sitting on it.  `edge`
  # is absent from an object built before it existed; the flags still say
  # that a bound chose, just not which.
  edge <- if ("edge" %in% names(x)) x$edge
          else ifelse(x$at_ceiling, "the last level of the ladder",
                      ifelse(x$at_floor, "the first level of the ladder", NA_character_))
  on_edge <- !is.na(edge)
  if (any(on_edge)) {
    # Wrapped, because with four criteria named one line reached 87
    # columns, which is the wrapping this block was moved here to avoid.
    cat("\n")
    for (e in unique(edge[on_edge]))
      cat(strwrap(sprintf("%s: the optimum is %s.",
                          paste(x$criterion[on_edge & edge == e], collapse = ", "), e),
                  width = 78, prefix = "  ", exdent = 2), sep = "\n")
    cat("  There the bound is choosing, not the criterion.\n")
  }
  # With one criterion there is nothing to compare: the spread is zero and the
  # intersection is that criterion's own band, so the closing block is skipped.
  if (nrow(x) > 1L) {
    # Recomputed from the rows in hand rather than read off the attribute.
    # `[` carries `common` through a row subset unchanged, so a two-criterion
    # slice of a four-criterion table used to print the parent's verdict: on
    # 25 simulated fields that was wrong for 39 of 250 subsets, denying an
    # overlap the two rows above it plainly showed.
    cm <- if (!is.null(bands) && all(x$criterion %in% names(bands)))
      sort(unique(Reduce(intersect, bands[x$criterion])))
    else sort(unique(attr(x, "common")))
    cat("\n")
    rng <- range(x$best)
    cat(sprintf("  picks span %d to %d cells (%.1fx)\n", rng[1L], rng[2L],
                rng[2L] / max(rng[1L], 1L)))
    if (length(cm))
      cat(sprintf("  in every flat region: %s (%d level%s)\n",
                  .band_label(cm, Reduce(intersect, lapply(x$criterion, rungs_of))),
                  length(cm), if (length(cm) == 1L) "" else "s"))
    else
      cat("  no level is in every flat region: the criteria disagree over the\n",
          "  whole ladder. plot() draws the curves they were read from.\n", sep = "")
  }
  invisible(x)
}


#' Resolve a cell count from a number or from a level-selection result
#'
#' The functions that need a cell count (\code{build_tessellation()}'s
#' \code{approx_n_cells}, \code{get_voronoi_seeds()}'s \code{n}) accept,
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
  edge_note <- function(sel) {
    e <- sel$edge %||% NA_character_
    if (is.na(e)) "" else sprintf(", an optimum at %s", e)
  }
  if (inherits(x, "resolution_selection")) {
    n <- x$best
    from <- sprintf("select_resolution(criterion = \"%s\")%s", x$criterion, edge_note(x))
  } else if (inherits(x, "resolution_profile")) {
    # The default criterion is cp, which a geometry-only profile cannot score;
    # fall through to the first criterion that is finite somewhere rather
    # than fail with cp's own message, which names neither this function
    # nor the way round it.
    usable <- Filter(function(cn) cn %in% names(x) &&
                       any(is.finite(suppressWarnings(as.numeric(x[[cn]])))),
                     c("cp", "reliability", "elbow", "moran_z"))
    if (!length(usable))
      stop(sprintf(paste0("%s(): `%s` is a resolution profile with no criterion finite at ",
                          "any level, so no cell count can be read off it."),
                   caller, arg), call. = FALSE)
    sel <- select_resolution(x, criterion = usable[[1L]])
    n <- sel$best
    from <- sprintf("resolution_profile() read with select_resolution(criterion = \"%s\")%s",
                    sel$criterion, edge_note(sel))
    .log_info("%s(): `%s` is a resolution profile; read with select_resolution(criterion = \"%s\"): %d cells.",
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
  # The integer ceiling matters as much as the floor: every consumer of this
  # count reaches as.integer() eventually, which turns anything above
  # .Machine$integer.max into NA and then aborts inside a clamp that names
  # neither the argument nor the cause.
  if (is.numeric(n) && length(n) == 1L && is.finite(n) && n > .Machine$integer.max)
    stop(sprintf(paste0("%s(): `%s` must resolve to at most %s cells (R's ",
                        "largest integer); got %s."),
                 caller, arg, format(.Machine$integer.max), format(n)),
         call. = FALSE)
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1)
    stop(sprintf("%s(): `%s` must resolve to a positive number of cells; got %s.",
                 caller, arg, paste(format(n), collapse = ", ")), call. = FALSE)
  list(n = as.numeric(n), from = from)
}
