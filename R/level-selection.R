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
#' The maximum-distance-to-chord rule: draw the chord from (k_min, WSS_min)
#' to (k_max, WSS_max) and take the k whose WSS sags furthest below it (see
#' \code{.below_chord()}).  This is the classical elbow construction.  It is
#' \emph{not} Kneedle (Satopaa et al. 2011), which takes the first local
#' maximum of the normalised difference curve that clears a sensitivity
#' threshold: the two agree on smooth single-knee curves and can disagree on
#' shouldered ones (for the 20-value curve 1000, 600, 560, 555, 552, 550,
#' 300, 100, 60, 50, 45, 42, 40, 39, 38, 37, 36, 35, 34, 33 this rule answers
#' k = 8 where Kneedle answers k = 2).
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

  # Knee FIRST.  The candidates were sorted ascending, which put the knee in
  # the middle: `k[1]` and `top_n = 1` -- documented as "the top-ranked
  # candidate" -- returned knee - 1 on every geometric call (the default, and
  # the fallback every model-aware call takes below the nine-cell floor).  The
  # help example asked for "two clearly separated clusters" and answered 1.
  # Same code in 1.0.0.  The model-aware path already ranks best-first, so
  # position 1 now means the same thing on both.
  .make_candidates <- function(knee) {
    if (!return_neighbors) return(knee)
    unique(pmin(max_k, pmax(min_k, c(knee, knee - 1L, knee + 1L))))
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


#' k-means++ seeding
#'
#' The first centre is a point drawn at random; each later centre is a point
#' drawn with probability proportional to its squared distance from the
#' nearest centre already chosen (Arthur and Vassilvitskii 2007).  Draws from
#' the current RNG stream, so the caller seeds it.
#'
#' @param xy Numeric matrix of coordinates.
#' @param k Number of centres, at most \code{nrow(xy)}.
#' @return A \code{k x ncol(xy)} matrix of starting centres.
#' @keywords internal
#' @noRd
.kmeanspp_centers <- function(xy, k) {
  n <- nrow(xy)
  idx <- integer(k)
  idx[1L] <- sample.int(n, 1L)
  if (k > 1L) {
    d2 <- rowSums((xy - matrix(xy[idx[1L], ], n, ncol(xy), byrow = TRUE))^2)
    for (j in 2:k) {
      tot <- sum(d2)
      # Every remaining point coincides with a centre: fall back to a uniform
      # draw among the points not yet chosen.
      idx[j] <- if (is.finite(tot) && tot > 0)
        sample.int(n, 1L, prob = d2 / tot)
      else sample(setdiff(seq_len(n), idx[seq_len(j - 1L)]), 1L)
      d2 <- pmin(d2, rowSums((xy - matrix(xy[idx[j], ], n, ncol(xy), byrow = TRUE))^2))
    }
  }
  xy[idx, , drop = FALSE]
}


#' The best of several k-means++ restarts, with their spread
#'
#' \code{stats::kmeans(nstart = )} restarts from uniform-random centres and
#' returns only the best solution.  Every point on a WSS curve built that way
#' is a local optimum, and with few restarts the curve mixes genuine level
#' effects with optimisation noise: on clustered layouts a sweep over
#' \code{k = 1..30} at \code{nstart = 5} had one or two \emph{increases} in
#' WSS in three of five draws, and the elbow rule once selected such a bump
#' (see \code{.elbow_from_wss()}).  Steinley (2003) shows far more restarts
#' are needed than practitioners use; Fränti and Sieranoja (2019) show that
#' k-means++ seeding cuts how many are needed and that the gain saturates,
#' which is the basis for a fixed budget.  With 25 k-means++ restarts the same
#' sweeps had no increase at all.
#'
#' Each restart is seeded by \code{.kmeanspp_centers()} and run as a single
#' start, so the spread of \code{tot.withinss} across restarts is available:
#' it says how rough the objective is at this \code{k}, which a profile can
#' report so a flat region is honestly wide.
#'
#' @param xy Numeric matrix of coordinates.
#' @param k Number of clusters (\code{>= 2}).
#' @param nstart Number of restarts.
#' @param iter.max Passed to \code{stats::kmeans()}.
#' @return \code{NULL} when every restart failed; otherwise a list with
#'   \code{km} (the best \code{kmeans} fit), \code{wss} (its
#'   \code{tot.withinss}), \code{spread} (\code{(max - min) / min} of
#'   \code{tot.withinss} across the restarts that ran) and \code{n_ok}.
#' @keywords internal
#' @noRd
.kmeans_best <- function(xy, k, nstart = 25L, iter.max = 50L) {
  best <- NULL; w <- rep(NA_real_, nstart)
  for (r in seq_len(nstart)) {
    km <- try(suppressWarnings(stats::kmeans(xy, centers = .kmeanspp_centers(xy, k),
                                             iter.max = iter.max, nstart = 1L)),
              silent = TRUE)
    if (inherits(km, "try-error")) next
    w[r] <- km$tot.withinss
    if (is.null(best) || km$tot.withinss < best$tot.withinss) best <- km
  }
  if (is.null(best)) return(NULL)
  ok <- w[is.finite(w)]
  list(km = best, wss = best$tot.withinss,
       spread = if (length(ok) > 1L && min(ok) > 0) (max(ok) - min(ok)) / min(ok) else 0,
       n_ok = length(ok))
}


#' Count the increases in a WSS curve
#'
#' A k-means WSS curve over increasing \code{k} can only rise where some
#' \code{k} landed in a worse local optimum than its neighbour, so
#' \code{sum(diff(wss) > 0)} is a direct count of optimisation bumps: zero
#' means the curve is at least monotone.  \code{.elbow_from_wss()} already
#' computes the differences; this is the check nobody read off them.
#'
#' @keywords internal
#' @noRd
.wss_bumps <- function(wss) {
  wss <- as.numeric(wss)
  if (length(wss) < 2L) return(0L)
  sum(diff(wss) > 0, na.rm = TRUE)
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
#' \strong{The WSS curve is read for its shape, so it is fitted to be
#' smooth.}  Every point on it is a k-means local optimum, and with a few
#' random restarts the curve mixes level effects with optimisation noise:
#' measured on clustered layouts, a sweep over \code{k = 1..30} at
#' \code{stats::kmeans(nstart = 5)} \emph{rose} at one or two steps in three
#' of five draws, and an earlier form of the elbow rule once selected such a
#' bump.  Each \code{k} is therefore fitted as the best of 25 restarts seeded
#' by k-means++ (Arthur and Vassilvitskii 2007), the budget at which the gain
#' from further restarts saturates (Fränti and Sieranoja 2019; Steinley 2003
#' on why the usual handful is not enough).  The same sweeps then had no
#' increase at all.  A curve that still rises somewhere is reported --- a
#' logged warning names the number of rising steps, and the model-aware
#' paths return it as \code{wss_bumps} in the \code{"diagnostics"} attribute
#' beside \code{wss_spread}, the relative spread of WSS across the restarts
#' at each \code{k} --- but not refused: a bumpy curve is uncertain, not
#' unidentified.  Because the optimiser changed, a selection made by an
#' earlier version on a curve that had such a bump can differ from the one
#' made now; where the earlier curve was clean, the answer is the same.
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
#' alone.  The criterion ranks on \eqn{|z|}, not on \eqn{|I|}, and at the
#' floor the residual moments give \eqn{E[I] = I} and \eqn{\mathrm{Var}[I] =
#' 0} identically (the algebra holds to \eqn{10^{-16}}), so the standardised
#' deviate is \eqn{0/0}: it carries no information about the tessellation, and
#' whichever way rounding noise resolves it those candidates would rank first
#' or last on nothing.  They therefore return \code{NA} and are excluded from
#' the model-aware ranking.  When no candidate in the elbow neighbourhood clears
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
#'   nine-cell resolution floor described in \strong{Details}.  Note that
#'   supplying both \code{response_var} and \code{predictor_vars} upgrades
#'   \code{"geometric"} to \code{"combined"}: the selection then depends on
#'   the response (see "Post-selection inference").
#' @param select_on \code{"all"} (default) selects on every point;
#'   \code{"split"} selects on one spatially blocked half of the points and
#'   returns the other half as the set to estimate on, so that the standard
#'   errors computed downstream on the chosen cells are not post-selection.
#'   See "Post-selection inference".
#' @section Post-selection inference:
#' When the selection reads the response --- here, whenever both
#' \code{response_var} and \code{predictor_vars} are supplied --- everything
#' estimated afterwards on the chosen cells is estimated on data that already
#' influenced the choice, and its standard errors are post-selection ones:
#' descriptive, not at nominal coverage (Gao, Bien and Witten 2022; Chen and
#' Witten 2023 give the exact selective test for two k-means clusters, the
#' first link of this chain).  The exposure is narrower than "the partition
#' was chosen on the response": every partition here is k-means on the
#' coordinates alone, and the response only decides which \emph{count} is
#' ranked first.  It is not zero, because the count determines every cell the
#' downstream standard errors are computed over.
#'
#' \code{select_on = "split"} is sample splitting: the layer is cut into two
#' spatially blocked halves (\code{\link{make_folds}(k = 2, method =
#' "block_kfold")}), the selection runs on the first half only, and the row
#' positions of both halves come back in the \code{"split"} attribute
#' (\code{selection} and \code{estimation}).  Build the tessellation on
#' every point --- cells are geometry --- but aggregate and fit on
#' \code{data_sf[attr(x, "split")$estimation, ]}, which the selection never
#' saw; that restores nominal coverage with no new theory.  The price is
#' precision: half the points estimate, and García Rasines and Young (2023)
#' show a \emph{contiguous} spatial half is less efficient than the
#' exchangeable split the i.i.d. theory assumes, because the two halves are
#' not interchangeable.  Two alternatives keep the whole sample --- data
#' thinning for count responses (Neufeld et al. 2024) and data fission for
#' Gaussian-like ones (Leiner et al. 2023) --- and are not implemented here;
#' the split needs no distributional assumption, which is why it comes first.
#' Selection on coordinates alone (\code{"geometric"} with no response) is
#' not exposed in this way, and \code{"split"} then changes nothing but the
#' attribute.
#' @return An integer vector of candidate level counts, \strong{best first}:
#'   under the geometric criterion the elbow, then its lower and upper
#'   neighbours; under the model-aware criteria the candidates in rank order.
#'   \code{k[1]} is therefore the top-ranked count on every path, and
#'   \code{top_n = 1} returns it alone. When
#'   \code{criterion != "geometric"}, an attribute \code{"diagnostics"} is
#'   attached with per-k Moran's I values (\code{moran_i}) and their
#'   standardised deviates (\code{moran_z}), the WSS curve (\code{wss}) with
#'   the relative between-restart spread at each \code{k} (\code{wss_spread}),
#'   the number of rising steps on it (\code{wss_bumps}) and the restart
#'   budget (\code{nstart}) — except when the model-aware path
#'   itself falls back to the geometric result (no viable k in the elbow
#'   neighbourhood, or Moran's I could not be computed for any candidate), in
#'   which case no diagnostics are available and the attribute is absent. Both
#'   fallbacks are logged as warnings.  The geometric path returns a plain
#'   integer vector; a rising WSS curve is still logged there.  With
#'   \code{select_on = "split"} every path adds a \code{"split"} attribute:
#'   a list with \code{selection} and \code{estimation} (integer row
#'   positions in \code{data_sf}), \code{method} and \code{seed}.  For a full
#'   per-level table --- criteria, cell support, restart spread, the flat
#'   region --- see \code{\link{resolution_profile}()}.
#' @references
#' Arthur, D. and Vassilvitskii, S. (2007). k-means++: the advantages of
#' careful seeding. \emph{Proceedings of the 18th Annual ACM-SIAM Symposium
#' on Discrete Algorithms}, 1027--1035.
#'
#' Fränti, P. and Sieranoja, S. (2019). How much can k-means be improved by
#' using better initialization and repeats? \emph{Pattern Recognition}, 93,
#' 95--112. \doi{10.1016/j.patcog.2019.04.014}
#'
#' Steinley, D. (2003). Local optima in K-means clustering: what you don't
#' know may hurt you. \emph{Psychological Methods}, 8(3), 294--304.
#' \doi{10.1037/1082-989X.8.3.294}
#'
#' Gao, L. L., Bien, J. and Witten, D. (2022). Selective inference for
#' hierarchical clustering. \emph{Journal of the American Statistical
#' Association}, 119, 332--342. \doi{10.1080/01621459.2022.2116331}
#'
#' Chen, Y. T. and Witten, D. M. (2023). Selective inference for k-means
#' clustering. \emph{Journal of Machine Learning Research}, 24(152), 1--41.
#' \url{https://jmlr.org/papers/v24/22-0371.html}
#'
#' García Rasines, D. and Young, G. A. (2023). Splitting strategies for
#' post-selection inference. \emph{Biometrika}, 110(3), 597--614.
#' \doi{10.1093/biomet/asac070}
#'
#' Leiner, J., Duan, B., Wasserman, L. and Ramdas, A. (2023). Data fission:
#' splitting a single data point. \emph{Journal of the American Statistical
#' Association}, 120(549), 135--146. \doi{10.1080/01621459.2023.2270748}
#'
#' Neufeld, A., Dharamshi, A., Gao, L. L. and Witten, D. (2024). Data thinning
#' for convolution-closed distributions. \emph{Journal of Machine Learning
#' Research}, 25(57), 1--35. \url{https://jmlr.org/papers/v25/23-0446.html}
#' @examples
#' library(sf)
#' set.seed(1)
#' # Two clearly separated clusters: the elbow should sit near k = 2
#' pts <- st_as_sf(
#'   data.frame(x = c(runif(25, 0, 10), runif(25, 90, 100)),
#'              y = c(runif(25, 0, 10), runif(25, 90, 100))),
#'   coords = c("x", "y"), crs = 32632
#' )
#' determine_optimal_levels(pts, max_levels = 6)   # 2 1 3: the elbow first
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
                                                    "combined"),
                                     select_on = c("all", "split")) {
  if (!inherits(data_sf, "sf"))
    stop("determine_optimal_levels(): `data_sf` must be an sf object.")

  criterion <- match.arg(criterion)
  select_on <- match.arg(select_on)
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

  # Sample splitting: select on one spatial half, hand the other back for
  # estimation.  Done on the full layer, before the subsample, so the
  # positions returned index `data_sf` as the caller passed it.
  split <- NULL
  if (identical(select_on, "split")) {
    split <- .spatial_half_split(data_sf, seed = set_seed,
                                 caller = "determine_optimal_levels")
    data_sf <- data_sf[split$selection, , drop = FALSE]
  }
  .with_split <- function(out) {
    if (!is.null(split)) attr(out, "split") <- split
    out
  }

  xy <- sf::st_coordinates(data_sf)[, 1:2, drop = FALSE]
  n  <- nrow(xy)
  if (n < 3L) return(.with_split(1L))

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
  if (k_max < 2L) return(.with_split(1L))

  # Say it BEFORE the sweep: the model-aware criteria carry no information at
  # nine cells or fewer (see the resolution floor in Details), so with
  # k_max <= 9 no candidate can ever clear it and the call is going to fall
  # back to the geometric ranking whatever the data say.  The sweep used to
  # run first and the fallback was logged afterwards.
  if (criterion != "geometric" && k_max <= 9L)
    .log_warn(paste0("determine_optimal_levels(): max_levels leaves k_max = %d, ",
                     "and the model-aware criteria carry no information at ",
                     "nine cells or fewer; criterion = '%s' will fall back to ",
                     "the geometric ranking. Raise max_levels above 9 (well ",
                     "above, so the elbow neighbourhood reaches past the ",
                     "floor) for the criteria to contribute."),
              k_max, criterion)

  # 25 k-means++ restarts per k rather than stats::kmeans(nstart = 5): the
  # WSS curve is read for its shape, and optimisation noise on it is what
  # produced the concave bump the elbow rule once selected.  See
  # .kmeans_best() for the measured effect and the citations.
  nstart <- 25L
  wss <- numeric(k_max)
  wss_spread <- numeric(k_max)
  failed_k <- integer(0)

  for (k in seq_len(k_max)) {
    if (k == 1L) {
      ctr <- colMeans(xy)
      wss[k] <- sum(rowSums((xy - matrix(ctr, nrow(xy), 2, byrow = TRUE))^2))
    } else {
      # Assignments are discarded on purpose: only the WSS is kept from this
      # sweep, which bounds memory for large max_levels (see Details).
      kb <- .kmeans_best(xy, k, nstart = nstart)
      if (is.null(kb)) {
        wss[k] <- wss[k - 1L]
        failed_k <- c(failed_k, k)
      } else {
        wss[k] <- kb$wss
        wss_spread[k] <- kb$spread
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
    if (k_max < 2L) return(.with_split(1L))
  }

  elbow <- .elbow_from_wss(wss, max_k = k_max, min_k = 1L, return_neighbors = TRUE)

  # A WSS curve that rises anywhere is one where some k landed in a worse
  # optimum than its neighbour, and the elbow read from it is partly noise.
  # Reported, not refused: a bumpy curve is uncertain, not unidentified.
  wss_bumps <- .wss_bumps(wss[seq_len(k_max)])
  if (wss_bumps > 0L)
    .log_warn(paste0("determine_optimal_levels(): the WSS curve rises at %d ",
                     "step(s) between k = 1 and %d, so some k landed in a ",
                     "worse local optimum than its neighbour even with %d ",
                     "k-means++ restarts; the elbow read from it is partly ",
                     "optimisation noise. Treat the neighbouring candidates ",
                     "as equivalent."),
              wss_bumps, k_max, nstart)

  if (criterion == "geometric") {
    out <- as.integer(head(elbow$candidates, max(1L, as.integer(top_n))))
    out[out < 1L]    <- 1L
    out[out > k_max] <- k_max
    return(.with_split(unique(out)))
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
    return(.with_split(unique(out)))
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
    kb <- .kmeans_best(xy, k, nstart = nstart)
    if (is.null(kb)) next
    km <- kb$km
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
    return(.with_split(unique(out)))
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
      return(.with_split(unique(out)))
    }
    ranked <- finite_ks[order(abs(moran_z[finite_ks]))]
    out <- as.integer(head(ranked, max(1L, as.integer(top_n))))
    out[out < 1L] <- 1L; out[out > k_max] <- k_max
    out <- unique(out)
    attr(out, "diagnostics") <- list(moran_i = moran_vals, moran_z = moran_z,
                                      wss = wss[1:k_max],
                                      wss_spread = wss_spread[1:k_max],
                                      wss_bumps = wss_bumps, nstart = nstart,
                                      eval_ks = eval_ks)
    return(.with_split(out))
  }

  # --- Combined: rank-average of WSS elbow distance and |z| of Moran's I ---
  # |z|, not |I|: I's attainable range is set by the eigenvalues of the weights
  # matrix, which is rebuilt at every k, so |I| is not comparable across k.
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
    wss_spread = wss_spread[1:k_max],
    wss_bumps = wss_bumps, nstart = nstart,
    combined_rank = stats::setNames(combined_rank, eval_ks),
    eval_ks = eval_ks,
    criterion = "combined"
  )
  .with_split(out)
}
