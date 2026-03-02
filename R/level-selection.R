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
    perp_dist <- abs((y2 - y1) * k_norm - (x2 - x1) * wss_norm +
                       x2 * y1 - y2 * x1) / line_len
    knee_k <- k_idx[which.max(perp_dist)]
  }

  d1 <- diff(wss_k)
  d2 <- diff(d1)

  list(knee_k = knee_k, candidates = .make_candidates(knee_k),
       diagnostics = list(wss = wss_k, d1 = d1, d2 = d2))
}


#' Compute Moran's I for residuals at a given tessellation resolution
#'
#' For a given k-means cluster assignment, fits OLS on cell-level means and
#' computes Moran's I on the residuals using a k-nearest-neighbour (k = 8)
#' binary weight matrix, row-standardised.
#' Lower absolute Moran's I suggests the tessellation resolution adequately
#' captures the spatial structure in the data.
#'
#' @param xy Numeric matrix of coordinates.
#' @param response Numeric vector of response values.
#' @param predictors Numeric matrix of predictor values.
#' @param cluster_ids Integer vector of cluster assignments.
#' @param k Number of clusters.
#' @return Numeric scalar: Moran's I statistic (values near 0 indicate the
#'   tessellation resolution captures the spatial pattern; positive values
#'   indicate residual spatial autocorrelation remains).
#' @keywords internal
#' @noRd
.morans_i_for_k <- function(xy, response, predictors, cluster_ids, k) {
  # Aggregate to cell-level means
  cell_ids <- sort(unique(cluster_ids))
  n_cells <- length(cell_ids)
  if (n_cells < 4L) return(NA_real_)

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
  if (sum(ok) < 4L) return(NA_real_)

  fit <- try(stats::lm.fit(x = cbind(1, cell_pred[ok, , drop = FALSE]),
                            y = cell_resp[ok]),
             silent = TRUE)
  if (inherits(fit, "try-error")) return(NA_real_)
  resid <- fit$residuals
  n <- length(resid)

  # k-nearest-neighbour weight matrix via shared helper (sparse when possible)
  n_neighbors <- min(8L, n - 1L)
  if (n_neighbors < 1L) return(NA_real_)

  W <- .build_knn_weights(cell_xy[ok, , drop = FALSE], k = n_neighbors)

  # Moran's I = (n / S0) * (e' W e) / (e' e)
  S0 <- sum(W)
  if (S0 < .Machine$double.eps || sum(resid^2) < .Machine$double.eps)
    return(NA_real_)
  resid_c <- resid - mean(resid)
  I <- (n / S0) * as.numeric(crossprod(resid_c, W %*% resid_c)) / sum(resid_c^2)
  I
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
#' @param data_sf An sf object.
#' @param max_levels Integer upper bound on levels. Default 12.
#' @param top_n Integer; how many candidates to return. Default 3.
#' @param sample_n Integer; subsample size for speed. Default 1500.
#' @param set_seed Integer RNG seed. Default 123.
#' @param response_var Optional response column name. When provided alongside
#'   \code{predictor_vars}, enables model-aware level selection via Moran's I
#'   on OLS residuals.
#' @param predictor_vars Optional predictor column names.
#' @param criterion One of \code{"geometric"} (default when no response given),
#'   \code{"morans_i"} (select k that minimizes |Moran's I|), or
#'   \code{"combined"} (rank-average of WSS elbow distance and |Moran's I|).
#'   Falls back to \code{"geometric"} if response/predictors are unavailable.
#' @return An integer vector of candidate level counts. When
#'   \code{criterion != "geometric"}, an attribute \code{"diagnostics"} is
#'   attached with per-k Moran's I values.
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

  if (!all(sf::st_geometry_type(data_sf, by_geometry = TRUE) %in%
           c("POINT", "MULTIPOINT"))) {
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
    resp_vec <- as.numeric(df[[response_var]])
    pred_mat <- as.matrix(df[, predictor_vars, drop = FALSE])
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

  # Run k-means only for the candidate k values and compute Moran's I
  moran_vals <- rep(NA_real_, k_max)
  for (k in eval_ks) {
    km <- try(stats::kmeans(xy, centers = k, iter.max = 50, nstart = 5),
              silent = TRUE)
    if (inherits(km, "try-error")) next
    moran_vals[k] <- .morans_i_for_k(xy, resp_vec, pred_mat,
                                       km$cluster, k)
  }

  valid_moran <- is.finite(moran_vals[eval_ks])

  if (!any(valid_moran)) {
    .log_warn("determine_optimal_levels(): Moran's I could not be computed; falling back to geometric.")
    out <- as.integer(head(elbow$candidates, max(1L, as.integer(top_n))))
    out[out < 1L] <- 1L; out[out > k_max] <- k_max
    return(unique(out))
  }

  if (criterion == "morans_i") {
    # Select k that minimizes |Moran's I| among evaluated candidates.
    abs_moran <- rep(Inf, k_max)
    abs_moran[eval_ks] <- abs(moran_vals[eval_ks])
    abs_moran[!is.finite(abs_moran)] <- Inf
    ranked <- order(abs_moran)
    out <- as.integer(head(ranked, max(1L, as.integer(top_n))))
    out[out < 1L] <- 1L; out[out > k_max] <- k_max
    out <- unique(out)
    attr(out, "diagnostics") <- list(moran_i = moran_vals, wss = wss[1:k_max],
                                      eval_ks = eval_ks)
    return(out)
  }

  # --- Combined: rank-average of WSS elbow distance and |Moran's I| ---
  # Only rank over the evaluated neighbourhood to keep dimensions aligned.
  k_norm   <- (eval_ks - min(eval_ks)) / max(1, max(eval_ks) - min(eval_ks))
  wss_sub  <- wss[eval_ks]
  wss_norm <- (wss_sub - min(wss_sub)) / max(.Machine$double.eps, max(wss_sub) - min(wss_sub))
  x1 <- k_norm[1]; y1 <- wss_norm[1]
  x2 <- k_norm[length(k_norm)]; y2 <- wss_norm[length(wss_norm)]
  line_len <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
  if (line_len < .Machine$double.eps) {
    perp_dist <- rep(0, length(eval_ks))
  } else {
    perp_dist <- abs((y2 - y1) * k_norm - (x2 - x1) * wss_norm +
                       x2 * y1 - y2 * x1) / line_len
  }

  # Rank both criteria (lower rank = better)
  rank_elbow <- rank(-perp_dist, ties.method = "average")  # higher distance = better
  abs_moran_sub <- abs(moran_vals[eval_ks])
  abs_moran_sub[!is.finite(abs_moran_sub)] <- max(abs_moran_sub[is.finite(abs_moran_sub)], 1) + 1
  rank_moran <- rank(abs_moran_sub, ties.method = "average")  # lower |I| = better

  combined_rank <- (rank_elbow + rank_moran) / 2
  best_idx <- order(combined_rank)
  out <- as.integer(eval_ks[head(best_idx, max(1L, as.integer(top_n)))])
  out[out < 1L] <- 1L; out[out > k_max] <- k_max
  out <- unique(out)
  attr(out, "diagnostics") <- list(
    moran_i = moran_vals, wss = wss[1:k_max],
    combined_rank = stats::setNames(combined_rank, eval_ks),
    eval_ks = eval_ks,
    criterion = "combined"
  )
  out
}
