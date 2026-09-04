#' Signed perpendicular distance of points below the chord of a curve
#'
#' Positive when \code{(x, y)} lies below the straight line from
#' \code{(x1, y1)} to \code{(x2, y2)}, negative above it, scaled to a
#' perpendicular distance.  The knee of a decreasing WSS curve is the point of
#' greatest sag \emph{below} that chord; an unsigned distance would let a
#' concave bump above it win instead.
#'
#' @keywords internal
#' @noRd
.below_chord <- function(x, y, x1, y1, x2, y2, line_len) {
  y_chord <- y1 + (y2 - y1) * (x - x1) / (x2 - x1)
  (y_chord - y) * (x2 - x1) / line_len
}


#' Select an elbow (knee) from a WSS curve
#'
#' Heuristically selects the "elbow" from a vector of within-cluster sum of
#' squares (WSS) values as a function of cluster count k.
#'
#' Replaced with the standard perpendicular-distance-to-line method: draw a
#' line from (k_min, WSS_min) to (k_max, WSS_max), and pick the k whose
#' WSS deviates most from that line. This is robust to smooth curves and
#' matches the widely-used "kneedle" approach.
#'
#' @param wss Numeric vector of WSS indexed by k.
#' @param max_k Integer upper bound on k.
#' @param min_k Integer lower bound on k.
#' @param return_neighbors Logical; return neighboring k values.
#' @return A list with knee_k, candidates, diagnostics.
#' @keywords internal
#' @noRd
.elbow_from_wss <- function(wss, max_k = length(wss), min_k = 1L,
                            return_neighbors = TRUE) {
  if (!is.numeric(wss) || length(wss) < 2L)
    stop(".elbow_from_wss(): `wss` must be numeric length >= 2.")
  min_k <- as.integer(max(1L, min_k))
  max_k <- as.integer(min(length(wss), max_k))
  if (min_k >= max_k) stop(".elbow_from_wss(): need at least two k values.")

  k_idx <- seq.int(min_k, max_k)
  wss_k <- as.numeric(wss[k_idx])

  .make_candidates <- function(knee) {
    if (!return_neighbors) return(knee)
    sort(unique(pmin(max_k, pmax(min_k, c(knee - 1L, knee, knee + 1L)))))
  }

  if (length(wss_k) < 3L) {
    knee_k <- floor((min_k + max_k) / 2)
    return(list(
      knee_k = knee_k, candidates = .make_candidates(knee_k),
      diagnostics = list(wss = wss_k, d1 = diff(wss_k), d2 = numeric(0))
    ))
  }
  
  k_norm   <- (k_idx - min(k_idx)) / max(1, max(k_idx) - min(k_idx))
  wss_norm <- (wss_k - min(wss_k)) / max(.Machine$double.eps, max(wss_k) - min(wss_k))

  # Line from first point to last point
  x1 <- k_norm[1];   y1 <- wss_norm[1]
  x2 <- k_norm[length(k_norm)]; y2 <- wss_norm[length(wss_norm)]
  # Perpendicular distance from each point to the line
  line_len <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
  if (line_len < .Machine$double.eps) {
    # Degenerate: constant WSS
    knee_k <- floor((min_k + max_k) / 2)
  } else {
    # SIGNED deviation, positive BELOW the chord.  A WSS curve is decreasing
    # and (nearly) convex, so the knee is the point that sags furthest under
    # the line from first to last -- but abs() let a point ABOVE the chord
    # (a concave bump: a k where k-means fell into a worse local optimum than
    # its neighbours) win with the same magnitude, and it did.  Measured on a
    # curve with one such bump, the "knee" was the bump.
    perp_dist <- .below_chord(k_norm, wss_norm, x1, y1, x2, y2, line_len)
    knee_k <- k_idx[which.max(perp_dist)]
  }

  d1 <- diff(wss_k)
  d2 <- diff(d1)

  list(knee_k = knee_k, candidates = .make_candidates(knee_k),
       diagnostics = list(wss = wss_k, d1 = d1, d2 = d2))
}


#' The "not computable" return of .morans_i_for_k()
#'
#' Kept as a named constant so every early exit has the same shape as the
#' success path; a bare \code{NA_real_} from one of them would silently make
#' \code{moran_z[k]} the \emph{I} of the next candidate.
#'
#' @return \code{c(I = NA_real_, z = NA_real_)}.
#' @keywords internal
#' @noRd
.morans_na <- function() c(I = NA_real_, z = NA_real_)


#' Compute Moran's I for residuals at a given tessellation resolution
#'
#' For a given k-means cluster assignment, fits OLS on cell-level means and
#' computes Moran's I on the residuals using a k-nearest-neighbour (k = 8)
#' binary weight matrix, row-standardised, together with its standardised
#' deviate under the Cliff & Ord regression-residual moments.
#'
#' \strong{Rank on \code{z}, not on \code{I}.}  \eqn{E[I]} and \eqn{Var[I]}
#' both depend on the number of cells, so \eqn{|I|} shrinks as \code{k} grows
#' whether or not the tessellation is capturing anything.  Over 300 replicates
#' of a response with no spatial structure, mean \eqn{|I|} fell monotonically
#' from 0.114 at \code{k = 10} to 0.050 at \code{k = 60}; mean \eqn{|z|} over
#' the same runs was 0.77, 0.78, 0.76, 0.80 against the theoretical
#' \eqn{E|N(0,1)| = 0.798}.  Ranking on \eqn{|I|} therefore prefers the finest
#' tessellation for arithmetic reasons rather than statistical ones.
#'
#' @param xy Numeric matrix of coordinates.
#' @param response Numeric vector of response values.
#' @param predictors Numeric matrix of predictor values.
#' @param cluster_ids Integer vector of cluster assignments. The cluster count
#'   is derived from this vector, so it is not passed separately.
#' @return A named numeric vector \code{c(I = , z = )}: Moran's I on the
#'   cell-level OLS residuals, and its standardised deviate.  Values of
#'   \code{z} near 0 indicate the resolution captures the spatial pattern;
#'   positive values indicate residual spatial autocorrelation remains.  Both
#'   are \code{NA_real_} when the statistic is not computable at this \code{k}.
#' @keywords internal
#' @noRd
.morans_i_for_k <- function(xy, response, predictors, cluster_ids) {
  # Aggregate to cell-level means
  cell_ids <- sort(unique(cluster_ids))
  n_cells <- length(cell_ids)
  if (n_cells < 4L) return(.morans_na())

  cell_resp <- numeric(n_cells)
  cell_xy   <- matrix(0, n_cells, 2)
  cell_pred <- matrix(0, n_cells, ncol(predictors))

  for (j in seq_along(cell_ids)) {
    mask <- cluster_ids == cell_ids[j]
    cell_resp[j]    <- mean(response[mask], na.rm = TRUE)
    cell_xy[j, ]    <- colMeans(xy[mask, , drop = FALSE])
    cell_pred[j, ]  <- colMeans(predictors[mask, , drop = FALSE])
  }

  # Fit OLS on cell means
  ok <- is.finite(cell_resp) & apply(is.finite(cell_pred), 1, all)
  if (sum(ok) < 4L) return(.morans_na())

  fit <- try(stats::lm.fit(x = cbind(1, cell_pred[ok, , drop = FALSE]),
                            y = cell_resp[ok]),
             silent = TRUE)
  if (inherits(fit, "try-error")) return(.morans_na())
  resid <- fit$residuals
  n <- length(resid)

  # k-nearest-neighbour weight matrix via shared helper (sparse when possible)
  n_neighbors <- min(8L, n - 1L)
  if (n_neighbors < 1L) return(.morans_na())

  # Refuse to report a number that carries no information.  When every cell is
  # a neighbour of every other (n <= n_neighbors + 1, i.e. n <= 9 at the
  # default of 8), W is the complete row-standardised matrix W_ij = 1/(n-1),
  # so W %*% e = -e/(n-1) for ANY mean-zero residual vector, S0 = n, and
  # Moran's I collapses to exactly -1/(n - 1) whatever the data are.  It is not
  # merely uninformative but biased for level selection: |I| = 1/(n-1) falls
  # monotonically in the number of cells for arithmetic reasons alone, so
  # criterion = "morans_i" would rank the largest candidate k first every time.
  # NA excludes these candidates instead; determine_optimal_levels() falls back
  # to the geometric ranking when none of them clears the floor.
  if (n <= n_neighbors + 1L) return(.morans_na())

  W <- .build_knn_weights(cell_xy[ok, , drop = FALSE], k = n_neighbors)

  # Moran's I = (n / S0) * (e' W e) / (e' e)
  S0 <- sum(W)
  if (S0 < .Machine$double.eps || sum(resid^2) < .Machine$double.eps)
    return(.morans_na())
  resid_c <- resid - mean(resid)
  # sum(resid_c * (W %*% resid_c)) rather than crossprod(): .build_knn_weights()
  # returns a sparse Matrix when FNN and Matrix are installed, and
  # base::crossprod() does not dispatch on the dgeMatrix that W %*% resid_c
  # produces ("requires numeric/complex matrix/vector arguments").  The two
  # forms are numerically identical.  Matches residual_morans_i().
  I <- as.numeric((n / S0) * sum(resid_c * (W %*% resid_c)) / sum(resid_c^2))

  # The ranking needs a STANDARDISED deviate, not |I|.  E[I] and Var(I) both
  # depend on k, so |I| shrinks with the number of cells for arithmetic reasons
  # that have nothing to do with the data.  Measured over 300 replicates of a
  # response with NO spatial structure (n = 1200, two noise predictors), mean
  # |I| fell monotonically from 0.1136 at k = 10 to 0.0502 at k = 60 -- -55.9%
  # -- so ranking on |I| prefers the finest candidate whatever the data say.
  # Over the same runs mean |z| was 0.769, 0.778, 0.755, 0.801 against the
  # theoretical E|N(0,1)| = 0.798, with sd(z) 0.96-1.02 and a two-sided 5%
  # rejection rate of 0.040-0.057.  It is calibrated, and flat in k.
  #
  # These are Cliff & Ord's regression-residual moments, and they are EXACT
  # here: `resid` is by construction the OLS residual of the cell means on
  # cbind(1, cell_pred), which is the one case the formula is derived for.
  mom <- .morans_residual_moments(W = W, X = cbind(1, cell_pred[ok, , drop = FALSE]),
                                  S0 = S0, is_sparse = inherits(W, "Matrix"))
  # No usable moments means no usable z -- and z is what the ranking reads.
  # Returning a finite I beside an NA z would be the documented "not
  # computable" shape in one element and a number in the other.
  if (is.null(mom) || !is.finite(mom$VI) || mom$VI <= 0) return(.morans_na())
  c(I = I, z = (I - mom$EI) / sqrt(mom$VI))
}


#' Determine an optimal number of spatial levels via an elbow heuristic
#'
#' Computes a WSS curve over k=1..K_max using k-means on projected feature
#' coordinates and selects candidate k values around the elbow.
#'
#' When \code{response_var} and \code{predictor_vars} are provided, the
#' geometric WSS elbow is supplemented with Moran's I computed on OLS
#' residuals at each candidate k.  The Moran's I profile measures how much
#' spatial autocorrelation in the response remains *unexplained* at a given
#' tessellation resolution — a direct reflection of the spatial process being
#' modeled, rather than mere geometric compactness of coordinates.  The
#' combined criterion selects the k that best balances geometric parsimony
#' and residual spatial independence.
#'
#' To keep memory use and runtime bounded for large \code{max_levels}, the
#' initial k-means sweep records only within-cluster sum-of-squares (WSS)
#' without retaining cluster assignments.  Moran's I is then evaluated
#' lazily: k-means is re-run only for a focused neighbourhood around the
#' elbow (±4 by default, or ±\code{top_n} if larger), so that only the most
#' promising candidate k values incur the cost of the full Moran's I
#' computation.
#'
#' \strong{The model-aware criteria rank on the standardised deviate, not on
#' |Moran's I|.}  Both \eqn{E[I]} and \eqn{Var[I]} depend on the number of
#' cells, so \eqn{|I|} falls as \code{k} grows whether or not the finer
#' tessellation is capturing anything.  Measured over 300 replicates of a
#' response with \emph{no} spatial structure, mean \eqn{|I|} fell monotonically
#' from 0.114 at \code{k = 10} to 0.050 at \code{k = 60} (\eqn{-56\%}), which
#' made an \eqn{|I|} ranking prefer the largest candidate for arithmetic
#' reasons alone.  Candidates are therefore ordered by
#' \eqn{|z| = |I - E[I]| / \mathrm{sd}(I)} using the Cliff & Ord regression
#' residual moments --- exact here, because the cell-level residuals are OLS
#' residuals by construction.  Over the same runs \eqn{z} had mean \eqn{\approx
#' 0}, \eqn{\mathrm{sd} \approx 1} and a two-sided 5\% rejection rate of
#' 0.040--0.057 at every \code{k}.  Both quantities are reported in the
#' \code{"diagnostics"} attribute, as \code{moran_i} and \code{moran_z}.
#'
#' \strong{Resolution floor on the model-aware criteria.}  Moran's I is
#' computed on cell-level residuals with an 8-nearest-neighbour weight matrix,
#' so it only carries information once there are more than nine cells.  At nine
#' or fewer, every cell is a neighbour of every other, the row-standardised
#' weight matrix is complete, and Moran's I collapses to exactly
#' \eqn{-1/(k - 1)} for \emph{any} residual vector — a function of \code{k}
#' alone, and one whose magnitude shrinks monotonically with \code{k}, which
#' would make \code{criterion = "morans_i"} prefer the largest candidate every
#' time.  Those candidates therefore return \code{NA} and are excluded from the
#' model-aware ranking.  When no candidate in the elbow neighbourhood clears
#' the floor — which is the usual outcome for small \code{max_levels} — the
#' whole call falls back to the geometric ranking and logs a warning; raise
#' \code{max_levels} above roughly 10 if you want the model-aware criteria to
#' contribute.  Under \code{criterion = "combined"}, a candidate below the
#' floor that sits alongside candidates above it is ranked last on the Moran's
#' I axis while still competing on the geometric axis.
#'
#' @param data_sf An sf object.
#' @param max_levels Integer upper bound on levels. Default 12.
#' @param top_n Integer; how many candidates to return. Default 3. Under
#'   \code{criterion = "geometric"} the candidate set is the elbow and its two
#'   immediate neighbours, so at most 3 values are ever returned no matter how
#'   large \code{top_n} is; only the model-aware criteria can return more.
#' @param sample_n Integer; subsample size for speed. Default 1500.
#' @param set_seed Integer RNG seed. Default 123.
#' @param response_var Optional response column name. When provided alongside
#'   \code{predictor_vars}, enables model-aware level selection via Moran's I
#'   on OLS residuals. Must be numeric or logical (logicals are read as 0/1);
#'   a factor or character response raises an error rather than being coerced,
#'   because the residuals of an OLS fit to arbitrary level codes carry no
#'   meaning to test for autocorrelation.
#' @param predictor_vars Optional predictor column names. Must be numeric or
#'   logical (logicals are read as 0/1); factor/character columns raise an
#'   error.
#' @param criterion One of \code{"geometric"} (default when no response given),
#'   \code{"morans_i"} (select the k whose residual Moran's I is least
#'   \emph{significant}), or \code{"combined"} (rank-average of WSS elbow
#'   distance and that same quantity).  Falls back to \code{"geometric"} if
#'   response/predictors are unavailable, and also when no candidate clears the
#'   nine-cell resolution floor described in \strong{Details}.
#' @return An integer vector of candidate level counts. When
#'   \code{criterion != "geometric"}, an attribute \code{"diagnostics"} is
#'   attached with per-k Moran's I values (\code{moran_i}) and their
#'   standardised deviates (\code{moran_z}) — except when the model-aware path
#'   itself falls back to the geometric result (no viable k in the elbow
#'   neighbourhood, or Moran's I could not be computed for any candidate), in
#'   which case no diagnostics are available and the attribute is absent. Both
#'   fallbacks are logged as warnings.
#' @examples
#' library(sf)
#' set.seed(1)
#' # Two clearly separated clusters: the elbow should sit near k = 2
#' pts <- st_as_sf(
#'   data.frame(x = c(runif(25, 0, 10), runif(25, 90, 100)),
#'              y = c(runif(25, 0, 10), runif(25, 90, 100))),
#'   coords = c("x", "y"), crs = 32632
#' )
#' determine_optimal_levels(pts, max_levels = 6)
#' @family aggregation
#' @seealso [build_tessellation()], which takes the chosen level count as
#'   `approx_n_cells`; [assign_features_to_polygons()] and
#'   [summarize_by_cell()] for the steps that follow.
#' @export
determine_optimal_levels <- function(data_sf, max_levels = 12L, top_n = 3L,
                                     sample_n = 1500L, set_seed = 123L,
                                     response_var = NULL,
                                     predictor_vars = NULL,
                                     criterion = c("geometric", "morans_i",
                                                    "combined")) {
  if (!inherits(data_sf, "sf"))
    stop("determine_optimal_levels(): `data_sf` must be an sf object.")

  criterion <- match.arg(criterion)
  has_model_vars <- !is.null(response_var) && !is.null(predictor_vars) &&
    response_var %in% names(data_sf) &&
    all(predictor_vars %in% names(data_sf))

  # Auto-upgrade to combined when model variables are available
  if (has_model_vars && criterion == "geometric") {
    criterion <- "combined"
    .log_info("determine_optimal_levels(): response_var and predictor_vars supplied; using combined criterion (geometric + Moran's I).")
  }
  # Fall back if model variables not available for model-aware criteria
  if (!has_model_vars && criterion != "geometric") {
    .log_warn("determine_optimal_levels(): criterion='%s' requires response_var and predictor_vars; falling back to geometric.", criterion)
    criterion <- "geometric"
  }

  # MULTIPOINT must be coerced too, not merely admitted: st_coordinates()
  # returns one row per VERTEX, so any multi-vertex feature makes xy[i, ] a
  # different feature than row i of resp_vec/pred_mat, and every index below
  # reads the wrong rows silently.  Matches the guard in estimate_sac_range().
  if (!all(sf::st_geometry_type(data_sf, by_geometry = TRUE) == "POINT")) {
    data_sf <- coerce_to_points(data_sf, "auto")
  }
  data_sf <- ensure_projected(data_sf)

  xy <- sf::st_coordinates(data_sf)[, 1:2, drop = FALSE]
  n  <- nrow(xy)
  if (n < 3L) return(1L)

  cleanup <- .with_seed(set_seed)
  on.exit(cleanup(), add = TRUE)

  # Extract model variables before subsampling so indices stay aligned
  resp_vec <- pred_mat <- NULL
  if (has_model_vars) {
    df <- sf::st_drop_geometry(data_sf)
    # A factor or character predictor makes as.matrix() return a CHARACTER
    # matrix, which dies deep inside colMeans() with "'x' must be numeric".
    # Name the offending columns here instead.
    #
    # Logicals are NOT offending: as.matrix() on a logical column gives a
    # logical matrix, which colMeans() handles, and storage.mode() below makes
    # the 0/1 coding explicit.  fit_rf_model()/cv_rf()/predict() all accept
    # logical predictors, so rejecting them here would be inconsistent.
    non_num <- predictor_vars[!vapply(predictor_vars,
                                      function(v) is.numeric(df[[v]]) ||
                                        is.logical(df[[v]]),
                                      logical(1))]
    if (length(non_num)) {
      stop(sprintf(
        paste0("determine_optimal_levels(): `predictor_vars` must be numeric ",
               "or logical; %s %s not. Encode factor/character predictors ",
               "numerically (e.g. with model.matrix()) before calling."),
        paste(sprintf("'%s'", non_num), collapse = ", "),
        if (length(non_num) == 1L) "is" else "are"
      ), call. = FALSE)
    }
    # The same check the predictors get, for the same reason.  as.numeric() on
    # a factor returns its LEVEL CODES, so a factor response was silently
    # turned into an arbitrary integer relabelling of the categories and the
    # model-aware criteria ran an OLS on it: re-ordering the levels of the same
    # factor changed the chosen k (11 -> 12) and every moran_z.  A character
    # response becomes all-NA and is caught only downstream, where the message
    # blames the data rather than the column type.
    if (!(is.numeric(df[[response_var]]) || is.logical(df[[response_var]])))
      stop(sprintf(
        paste0("determine_optimal_levels(): `response_var` must be numeric or ",
               "logical; '%s' is %s. as.numeric() on a factor returns its level ",
               "CODES, so the model-aware criteria would be fitted to an ",
               "arbitrary relabelling of the categories."),
        response_var, class(df[[response_var]])[1L]), call. = FALSE)
    resp_vec <- as.numeric(df[[response_var]])
    pred_mat <- as.matrix(df[, predictor_vars, drop = FALSE])
    storage.mode(pred_mat) <- "double"
  }

  if (n > sample_n) {
    idx <- sample(seq_len(n), sample_n)
    xy <- xy[idx, , drop = FALSE]
    if (has_model_vars) {
      resp_vec <- resp_vec[idx]
      pred_mat <- pred_mat[idx, , drop = FALSE]
    }
  }

  k_max <- max(2L, min(as.integer(max_levels), nrow(xy) - 1L))
  
  n_uniq <- nrow(unique(round(xy, 8)))
  k_max <- min(k_max, n_uniq - 1L)
  if (k_max < 2L) return(1L)

  wss <- numeric(k_max)
  failed_k <- integer(0)

  for (k in seq_len(k_max)) {
    if (k == 1L) {
      ctr <- colMeans(xy)
      wss[k] <- sum(rowSums((xy - matrix(ctr, nrow(xy), 2, byrow = TRUE))^2))
    } else {
      km <- try(stats::kmeans(xy, centers = k, iter.max = 50, nstart = 5), silent = TRUE)
      if (inherits(km, "try-error")) {
        wss[k] <- wss[k - 1L]
        failed_k <- c(failed_k, k)
      } else {
        wss[k] <- km$tot.withinss
      }
    }
  }

  if (length(failed_k) > 0L) {
    .log_warn("determine_optimal_levels(): kmeans failed for k = %s; interpolating WSS.",
              paste(failed_k, collapse = ", "))
    for (fk in sort(failed_k)) {
      lo <- max(1L, fk - 1L)
      hi <- min(k_max, fk + 1L)
      while (hi %in% failed_k && hi < k_max) hi <- hi + 1L
      if (hi %in% failed_k) {
        k_max <- lo
        break
      }
      wss[fk] <- (wss[lo] + wss[hi]) / 2
    }
    if (k_max < 2L) return(1L)
  }

  elbow <- .elbow_from_wss(wss, max_k = k_max, min_k = 1L, return_neighbors = TRUE)

  if (criterion == "geometric") {
    out <- as.integer(head(elbow$candidates, max(1L, as.integer(top_n))))
    out[out < 1L]    <- 1L
    out[out > k_max] <- k_max
    return(unique(out))
  }

  # --- Model-aware criteria: Moran's I only for elbow neighbourhood ---
  # Rather than running k-means for every k (wasteful when k_max is large),

  # we evaluate a focused neighbourhood around the elbow.  A window of ±4
  # around the knee is wide enough to capture the Moran's I minimum near the
  # geometric elbow while avoiding O(k_max) redundant k-means fits.
  knee_k   <- elbow$knee_k
  margin   <- max(4L, as.integer(top_n))
  eval_ks  <- seq.int(max(2L, knee_k - margin), min(k_max, knee_k + margin))
  eval_ks  <- setdiff(eval_ks, failed_k)

  if (length(eval_ks) == 0L) {
    .log_warn("determine_optimal_levels(): no viable k values in elbow neighbourhood; falling back to geometric.")
    out <- as.integer(head(elbow$candidates, max(1L, as.integer(top_n))))
    out[out < 1L] <- 1L; out[out > k_max] <- k_max
    return(unique(out))
  }

  # Run k-means only for the candidate k values and compute Moran's I.
  # The WSS of the re-run clustering is recorded (wss_eval) so that the
  # combined ranking below compares elbow distance and Moran's I computed
  # on the *same* clustering — the sweep's RNG state differs, so its WSS
  # can come from a different local optimum than the Moran evaluation.
  moran_vals <- rep(NA_real_, k_max)   # raw I, reported in $diagnostics
  moran_z    <- rep(NA_real_, k_max)   # standardised deviate, used for ranking
  wss_eval   <- wss
  for (k in eval_ks) {
    km <- try(stats::kmeans(xy, centers = k, iter.max = 50, nstart = 5),
              silent = TRUE)
    if (inherits(km, "try-error")) next
    wss_eval[k]   <- km$tot.withinss
    mi            <- .morans_i_for_k(xy, resp_vec, pred_mat, km$cluster)
    moran_vals[k] <- mi[["I"]]
    moran_z[k]    <- mi[["z"]]
  }

  # Ranking is on |z|.  |I| is not comparable across k -- see .morans_i_for_k().
  valid_moran <- is.finite(moran_z[eval_ks])

  if (!any(valid_moran)) {
    .log_warn("determine_optimal_levels(): Moran's I could not be computed; falling back to geometric.")
    out <- as.integer(head(elbow$candidates, max(1L, as.integer(top_n))))
    out[out < 1L] <- 1L; out[out > k_max] <- k_max
    return(unique(out))
  }

  if (criterion == "morans_i") {
    # Select k that minimizes |Moran's I| among evaluated candidates.
    #
    # Rank ONLY the candidates that actually produced a finite Moran's I.
    # Ranking all of 1:k_max and truncating to top_n padded the answer with k
    # values that were never evaluated: the unevaluated entries all sit at Inf,
    # order() breaks those ties by index, and head() then appended 1, 2, 3, ...
    # whenever top_n exceeded the number of finite candidates -- including k
    # below the resolution floor, and k = 1, which is not a tessellation.
    finite_ks <- eval_ks[is.finite(moran_z[eval_ks])]
    if (length(finite_ks) == 0L) {
      .log_warn(paste0("determine_optimal_levels(): no candidate produced a ",
                       "finite Moran's I (every candidate is at or below the ",
                       "resolution floor); falling back to geometric."))
      out <- as.integer(head(elbow$candidates, max(1L, as.integer(top_n))))
      out[out < 1L] <- 1L; out[out > k_max] <- k_max
      return(unique(out))
    }
    ranked <- finite_ks[order(abs(moran_z[finite_ks]))]
    out <- as.integer(head(ranked, max(1L, as.integer(top_n))))
    out[out < 1L] <- 1L; out[out > k_max] <- k_max
    out <- unique(out)
    attr(out, "diagnostics") <- list(moran_i = moran_vals, moran_z = moran_z,
                                      wss = wss[1:k_max], eval_ks = eval_ks)
    return(out)
  }

  # --- Combined: rank-average of WSS elbow distance and |Moran's I| ---
  # Only rank over the evaluated neighbourhood to keep dimensions aligned;
  # use wss_eval so both criteria reflect the same clustering per k.
  k_norm   <- (eval_ks - min(eval_ks)) / max(1, max(eval_ks) - min(eval_ks))
  wss_sub  <- wss_eval[eval_ks]
  wss_norm <- (wss_sub - min(wss_sub)) / max(.Machine$double.eps, max(wss_sub) - min(wss_sub))
  x1 <- k_norm[1]; y1 <- wss_norm[1]
  x2 <- k_norm[length(k_norm)]; y2 <- wss_norm[length(wss_norm)]
  line_len <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
  if (line_len < .Machine$double.eps) {
    perp_dist <- rep(0, length(eval_ks))
  } else {
    # Signed (positive below the chord), as in .select_elbow(); see there.
    perp_dist <- .below_chord(k_norm, wss_norm, x1, y1, x2, y2, line_len)
  }

  # Rank both criteria (lower rank = better)
  rank_elbow <- rank(-perp_dist, ties.method = "average")  # higher distance = better
  # |z|, not |I|: the two rank candidates differently and only |z| is
  # comparable across k.  See .morans_i_for_k().
  abs_moran_sub <- abs(moran_z[eval_ks])
  abs_moran_sub[!is.finite(abs_moran_sub)] <- max(abs_moran_sub[is.finite(abs_moran_sub)], 1) + 1
  rank_moran <- rank(abs_moran_sub, ties.method = "average")  # lower |z| = better

  combined_rank <- (rank_elbow + rank_moran) / 2
  best_idx <- order(combined_rank)
  out <- as.integer(eval_ks[head(best_idx, max(1L, as.integer(top_n)))])
  out[out < 1L] <- 1L; out[out > k_max] <- k_max
  out <- unique(out)
  attr(out, "diagnostics") <- list(
    moran_i = moran_vals, moran_z = moran_z, wss = wss[1:k_max],
    wss_eval = wss_eval[1:k_max],
    combined_rank = stats::setNames(combined_rank, eval_ks),
    eval_ks = eval_ks,
    criterion = "combined"
  )
  out
}
