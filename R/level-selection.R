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


#' How far a WSS curve sags below a power law, and whether that is an elbow
#'
#' Points with no cluster structure have a WSS close to \eqn{c/k}: each of
#' \eqn{k} cells covers about \eqn{1/k} of the extent, and a cell's mean
#' squared distance to its centre scales with its area.  On linear axes that
#' curve is convex everywhere, and the classical chord rule lands where
#' \eqn{c/k} sags furthest below its own chord, \eqn{k = \sqrt{k_{min}
#' k_{max}}}: the ladder chooses, whatever the data.  On log-log axes
#' \eqn{c/k} is a straight line, and so is any power law.  The sag is
#' therefore measured there: \eqn{\log} WSS against \eqn{\log k}, as the
#' vertical distance below the straight line joining the first and last
#' \eqn{k}, in natural-log units (0.22 means the WSS is a fifth below the
#' power law through the ends).  Separated clusters fall faster than a power
#' law until there is one cell per cluster and like \eqn{c/k} after, which is
#' a bend at the cluster count.
#'
#' \code{.ELBOW_MIN_SAG} is the least sag reported as an elbow:
#' \eqn{\log 1.25 \approx 0.22}.  Measured with 25 k-means++ restarts, 60 to
#' 1500 points and ladders to 3--40: uniform layouts over squares, discs,
#' triangles, a 1.5:1 rectangle, an L-shape, density gradients, jittered
#' lattices and Gaussian blobs sagged at most 0.16; two to ten separated
#' clusters at least 0.6 once the ladder passed the cluster count; four
#' touching clusters 0.12--0.7, so a bend that weak can go unreported on a
#' small sample.  An elongated extent also bends, at about its aspect ratio,
#' because the first cuts go across its long axis (0.1--0.2 for a 2:1
#' rectangle, 0.2--0.45 for 4:1, 0.19 for the county centroids of North
#' Carolina): that is its shape, not clusters, and past the threshold it is
#' reported as an elbow.
#' @param k,wss Level counts (increasing, positive) and their WSS.
#' @return The sag at each \code{k}; \code{NA} throughout when fewer than
#'   three levels or a non-positive WSS leave no line to measure against.
#' @keywords internal
#' @noRd
.ELBOW_MIN_SAG <- log(1.25)
.elbow_sag <- function(k, wss) {
  m <- length(k)
  y <- suppressWarnings(log(as.numeric(wss)))
  if (m < 3L || !all(is.finite(y))) return(rep(NA_real_, m))
  x <- log(as.numeric(k))
  y[1L] + (y[m] - y[1L]) * (x - x[1L]) / (x[m] - x[1L]) - y
}


#' The elbow of a WSS ladder, including a fall to zero
#'
#' \code{.elbow_sag()} on the levels whose WSS is positive, and the reading
#' both callers take from it.  A level whose WSS is zero has a cell on every
#' distinct location, which a ladder reaches when locations repeat (stations
#' visited many times).  It has no place on log axes, and log(0) would take
#' the line with it, so it is left out of the line.  But the fall to zero is
#' the sharpest bend a WSS curve can make: when the positive levels have no
#' elbow, the first zero level is the elbow, one cell per location.  Its sag
#' is measured with its WSS floored at \code{.ELBOW_ZERO_TOL} of the first
#' level's, below the line through the first level and the last positive one
#' (the power law \eqn{c/k} through the first level when there is no other
#' positive one).  Leaving the zero out and reading the rest would say
#' "no cluster structure" and answer 3 on five stations visited thirty times
#' each, where one metre of jitter answers 5.
#'
#' Zero is relative: at most \code{.ELBOW_ZERO_TOL} (\eqn{10^{-12}}) of the
#' first level's WSS.  k-means leaves floating-point residue where the exact
#' value is 0 (7.8e-17 at the station count of one layer of repeat visits),
#' and the log of that dragged the line down as far as log(0) would.
#'
#' @param k Levels, increasing, with \code{k[1]} the anchor of the line (1 in
#'   both callers).
#' @param wss Their WSS, \code{wss[1]} the anchor's (the total sum of
#'   squares at \code{k = 1}).
#' @return A list with \code{sag} (at every level; \code{NA} at a zero level
#'   unless it is the elbow), \code{structured}, \code{knee} (the level of
#'   greatest sag, \code{NA} when not structured) and \code{zero} (which
#'   levels count as zero).
#' @keywords internal
#' @noRd
.ELBOW_ZERO_TOL <- 1e-12
.elbow_read <- function(k, wss) {
  m  <- length(k)
  wss <- as.numeric(wss)
  w1 <- wss[1L]
  zero <- if (is.finite(w1) && w1 > 0)
    is.finite(wss) & wss >= 0 & wss <= .ELBOW_ZERO_TOL * w1
  else rep(FALSE, m)
  sag <- rep(NA_real_, m)
  sag[!zero] <- .elbow_sag(k[!zero], wss[!zero])
  structured <- any(is.finite(sag)) && max(sag, na.rm = TRUE) >= .ELBOW_MIN_SAG
  if (!structured && any(zero)) {
    z  <- which(zero)[1L]
    p  <- which(!zero)
    lw <- suppressWarnings(log(wss[p]))
    lk <- log(as.numeric(k[p]))
    if (all(is.finite(lw))) {
      np <- length(p)
      slope <- if (np >= 2L) (lw[np] - lw[1L]) / (lk[np] - lk[1L]) else -1
      sag[z] <- lw[1L] + slope * (log(as.numeric(k[z])) - lk[1L]) -
        log(.ELBOW_ZERO_TOL * w1)
      structured <- sag[z] >= .ELBOW_MIN_SAG
    }
  }
  list(sag = sag, structured = structured,
       knee = if (structured) k[which.max(sag)] else NA, zero = zero)
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
#' The chord is drawn on log-log axes (see \code{.elbow_sag()}), where a
#' curve with no cluster structure is straight, and a fall to a WSS of zero
#' counts as a bend (see \code{.elbow_read()}).  When the sag there does not
#' reach \code{.ELBOW_MIN_SAG} the curve has no elbow: \code{structured} is
#' \code{FALSE} and \code{knee_k} is the linear-axis chord rule's answer,
#' which on such a curve is set by \code{min_k} and \code{max_k}, not by the
#' data.  The caller says so.
#'
#' @param wss Numeric vector of WSS indexed by k.
#' @param max_k Integer upper bound on k.
#' @param min_k Integer lower bound on k.
#' @param return_neighbors Logical; return neighbouring k values.
#' @return A list with knee_k, candidates, structured (whether the curve has
#'   an elbow) and diagnostics (with the log-log \code{sag} at each k).
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

  # Read on log-log axes, with a fall to zero counted as a bend (see
  # .elbow_read()).  Two levels are too few for a line, but not for a fall
  # to zero: two stations visited many times have their elbow at 2.
  rd <- .elbow_read(k_idx, wss_k)
  if (length(wss_k) < 3L) {
    knee_k <- if (rd$structured) rd$knee else floor((min_k + max_k) / 2)
    return(list(
      knee_k = knee_k, candidates = .make_candidates(knee_k),
      structured = rd$structured,
      diagnostics = list(wss = wss_k, d1 = diff(wss_k), d2 = numeric(0),
                         sag = rd$sag)
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
  # The same rule on log-log axes decides.  On linear axes a curve with no
  # cluster structure, WSS ~ c / k, still has a point furthest below its
  # chord, near sqrt(min_k * max_k), and that is what used to be returned as
  # the elbow: 4 at the default max_levels = 12, 13 at 160, on any uniform
  # layer.  On log-log axes that curve is straight (see .elbow_sag()).  The
  # linear answer is kept only as the fallback, flagged, when there is no
  # bend there.
  # A k with a WSS of 0 (to within 1e-12 of the total) has a cell on every
  # distinct location, which the sweep reaches when locations repeat.  It is
  # left out of the line, since log(0) would take the whole line with it,
  # and when the rest of the curve has no elbow it is the elbow: the WSS
  # falling to zero is the sharpest bend there is.  It used to be only left
  # out, and five stations visited thirty times each then had "no cluster
  # structure" and got 3.  Anything else that is not positive still leaves
  # no line at all.
  sag <- rd$sag
  structured <- rd$structured
  if (structured) knee_k <- rd$knee

  d1 <- diff(wss_k)
  d2 <- diff(d1)

  list(knee_k = knee_k, candidates = .make_candidates(knee_k), structured = structured,
       diagnostics = list(wss = wss_k, d1 = d1, d2 = d2, sag = sag))
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
    # One vector per coordinate, and each draw by inverting the cumulative
    # squared distance: sample.int(prob = ) sorts all n weights to draw one
    # point, and the n x 2 matrices built for every distance update were the
    # rest of the cost.  Seeding took 73 percent of a resolution profile's
    # time; measured on 5000 points it fell from 0.63 s to 0.08 s at
    # k = 1000, and a profile of 3000 points from 35 s to 10 s.  The draw
    # has the same law, P(i) = d2[i] / sum(d2), but maps
    # the random number to a different point, so it seeds different centres
    # for the same seed than it did up to 2.0.0.
    cols <- lapply(seq_len(ncol(xy)), function(j) xy[, j])
    dist2 <- function(i) {
      s <- 0
      for (v in cols) s <- s + (v - v[i])^2
      s
    }
    d2 <- dist2(idx[1L])
    for (j in 2:k) {
      cs  <- cumsum(d2)
      tot <- cs[n]
      # Every remaining point coincides with a centre: fall back to a uniform
      # draw among the points not yet chosen.  A point at distance 0 adds
      # nothing to the cumulative sum, so the inversion never lands on one.
      idx[j] <- if (is.finite(tot) && tot > 0)
        min(n, findInterval(stats::runif(1L) * tot, cs) + 1L)
      else sample(setdiff(seq_len(n), idx[seq_len(j - 1L)]), 1L)
      d2 <- pmin(d2, dist2(idx[j]))
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
#' report so that a flat region is shown at its real width.
#'
#' With \code{init = "random"} each restart starts from \code{k} distinct
#' locations drawn at random (\code{.random_centers()}), the draws
#' \code{stats::kmeans(centers = k, nstart = nstart)} makes, so the best of
#' them is the fit that call returns, and their spread is kept as well.
#'
#' @param xy Numeric matrix of coordinates, in the caller's unit.  The
#'   restarts see them in the frame of \code{.kmeans_frame()}.
#' @param k Number of clusters (\code{>= 2}).
#' @param nstart Number of restarts.
#' @param iter.max Passed to \code{stats::kmeans()}.
#' @param init \code{"kmeans++"} (default) or \code{"random"}: how each
#'   restart's starting centres are drawn.
#' @return \code{NULL} when every restart failed; otherwise a list with
#'   \code{km} (the best \code{kmeans} fit), \code{wss} (its
#'   \code{tot.withinss}, both in the caller's coordinates), \code{spread}
#'   (\code{(max - min) / min} of \code{tot.withinss} across the restarts
#'   that ran, in the frame) and \code{n_ok}.
#' @keywords internal
#' @noRd
.kmeans_best <- function(xy, k, nstart = 25L, iter.max = 50L, init = "kmeans++") {
  best <- NULL; w <- rep(NA_real_, nstart)
  random <- identical(init, "random")
  # The starts are drawn from, run on and ranked by the coordinates in the
  # frame of .kmeans_frame(), which does not move with the unit or the
  # origin; the kept fit is taken back to the caller's coordinates.
  z <- .kmeans_frame(xy)$z
  # The distinct rows random starts are drawn from, found once.  Fewer of
  # them than centres is what stats::kmeans() refuses as "more cluster
  # centers than distinct data points".
  cn <- if (random) unique(z) else NULL
  if (random && nrow(cn) < k) return(NULL)
  for (r in seq_len(nstart)) {
    km <- if (random)
      try(suppressWarnings(stats::kmeans(z, centers = .random_centers(z, cn, k, nstart),
                                         iter.max = iter.max, nstart = 1L)),
          silent = TRUE)
    else
      try(suppressWarnings(stats::kmeans(z, centers = .kmeanspp_centers(z, k),
                                         iter.max = iter.max, nstart = 1L)),
          silent = TRUE)
    if (inherits(km, "try-error")) next
    w[r] <- km$tot.withinss
    if (is.null(best) || km$tot.withinss < best$tot.withinss) best <- km
  }
  if (is.null(best)) return(NULL)
  ok <- w[is.finite(w)]
  best <- .kmeans_back(best, xy)
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
#' tessellation on arithmetic grounds alone.
#'
#' @param xy Numeric matrix of coordinates.
#' @param response Numeric vector of response values.
#' @param predictors Numeric matrix of predictor values.
#' @param cluster_ids Integer vector of cluster assignments. The cluster count
#'   is derived from this vector, so it is not passed separately.
#' @param grouped Whether the cell means are formed by grouped sums, one pass
#'   over the rows for all the cells.  \code{FALSE} (default) makes a pass
#'   over all the rows for each cell, which costs rows times cells and is
#'   quadratic in the rows when the cells are a fixed size.  The two give the
#'   same means to rounding, so a caller whose numbers are held to the last
#'   bit keeps the default.
#' @param use_fnn,use_matrix Whether \pkg{FNN} and \pkg{Matrix} may be used
#'   for the neighbour weights of the cell means; handed to
#'   \code{.build_knn_weights()}, which stops above \code{.KNN_DENSE_MAX}
#'   cells without both.
#' @return A named numeric vector \code{c(I = , z = )}: Moran's I on the
#'   cell-level OLS residuals, and its standardised deviate.  With no
#'   predictors those residuals are the centred cell means, so this is the
#'   autocorrelation among the cell means: a \code{z} near 0 says that
#'   neighbouring cell means are unrelated, which happens when the cells are
#'   as wide as the range of the field or wider (and at about 10 to 20 cells
#'   whatever the data, where the 8-neighbour graph is close to complete),
#'   and a positive \code{z} that neighbouring cells resemble each other,
#'   which they do more and more as the cells get finer than the range.  It
#'   is not a measure of the structure the cells leave unexplained.  Both
#'   are \code{NA_real_} when the statistic is not computable at this \code{k}.
#' @keywords internal
#' @noRd
.morans_i_for_k <- function(xy, response, predictors, cluster_ids, grouped = FALSE,
                            use_fnn = .has_backend("FNN"),
                            use_matrix = .has_backend("Matrix")) {
  # Aggregate to cell-level means
  cell_ids <- sort(unique(cluster_ids))
  n_cells <- length(cell_ids)
  if (n_cells < 4L) return(.morans_na())

  # I and z do not change under a + b * y, so the response is taken where
  # the arithmetic reads its spread.  A response whose raw values lose their
  # spread to rounding (1e15 plus a spread of 1) is taken about its centre
  # (.regression_centre()), or the cell means and their regression work on
  # the rounding of the offset; the centre is 0 for every other response,
  # and so for the column resolution_profile() passes, which it has already
  # centred.  A response whose largest absolute value lies outside 2^-255
  # to 2^255 is divided by the power of two at or below it, which is exact,
  # or the squares below overflow or lose their precision.
  raw <- response
  centre <- .regression_centre(raw, bounded = FALSE)
  shift_scale <- function(r) {
    a <- suppressWarnings(max(abs(r[is.finite(r)])))
    p <- 1
    if (is.finite(a) && a > 0 && (a < 2^-255 || a > 2^255)) {
      p <- 2^min(floor(log2(a)), 1023)
      if (p > a) p <- p / 2
    }
    list(r = r / p, p = p)
  }
  rs <- shift_scale(raw - centre)
  response <- rs$r

  if (isTRUE(grouped)) {
    # The means of the loop below from one pass over the rows: sums by cell
    # over counts by cell.  A missing response is left out of its cell's mean,
    # as mean(na.rm = TRUE) leaves it out, and a cell with none gets NaN.
    g    <- match(cluster_ids, cell_ids)
    cnt  <- tabulate(g, nbins = n_cells)
    cell_means <- function(r) {
      seen <- !is.na(r)
      as.numeric(rowsum(ifelse(seen, r, 0), g)) / as.numeric(rowsum(as.numeric(seen), g))
    }
    cell_rows <- as.numeric(rowsum(as.numeric(!is.na(response)), g))
    cell_all  <- as.numeric(cnt)
    cell_xy   <- unname(rowsum(xy[, 1:2, drop = FALSE], g)) / cnt
    pred_means <- function(P) if (ncol(P)) unname(rowsum(P, g)) / cnt
                              else matrix(0, n_cells, 0L)
  } else {
    cell_xy   <- matrix(0, n_cells, 2)
    masks <- lapply(cell_ids, function(id) cluster_ids == id)
    cell_means <- function(r)
      vapply(masks, function(mask) mean(r[mask], na.rm = TRUE), numeric(1))
    cell_rows <- vapply(masks, function(mask) sum(!is.na(response[mask])), numeric(1))
    cell_all  <- vapply(masks, function(mask) as.numeric(sum(mask)), numeric(1))
    pred_means <- function(P) {
      m <- matrix(0, n_cells, ncol(P))
      for (j in seq_along(cell_ids)) m[j, ] <- colMeans(P[masks[[j]], , drop = FALSE])
      m
    }

    for (j in seq_along(cell_ids))
      cell_xy[j, ] <- colMeans(xy[masks[[j]], , drop = FALSE])
  }
  # A predictor at an offset (1e8 plus a spread of 1), which the fit below
  # drops as collinear with the intercept or fits on the rounding of the
  # offset, is taken about its centre (.centre_columns()) before its cell
  # means are formed; the regression does not change under a shift of a
  # predictor, and at 1e15 the cell means of the raw column would carry
  # rounding of the order of its spread.  The retry and the moments below
  # read the same design; the retry's test for an exact fit reads it with
  # the centres added back (.rounding_design()).
  predictors <- .centre_columns(predictors)
  centres <- attr(predictors, "centres")
  if (!is.null(centres)) centres <- c(0, centres)

  # Residuals of the rows within the rounding of a row: the predictors fit
  # the response exactly, and its cell means differ by rounding alone.  The
  # rows are fitted by OLS and read as resolution_profile() reads its own
  # residuals (.residuals_rounding() on the rows), about the response's
  # centre where the residuals of the raw values are under 1e-12 of them,
  # as the cell means are below: the rounding of a fit on 1e15 + 1e5 w + f
  # at its raw values holds the spread of f.  The
  # bounds on the cell means below take the rounding of each row as its
  # own, and so shrink with the rows of a cell; where the rounding is the
  # same over the rows of a cell, as on a response at 1e15 that changes by
  # a unit in the last place over the layer, they call it spread.
  rows_rounding <- function() {
    P <- as.matrix(predictors)
    okr <- is.finite(response)
    if (ncol(P)) okr <- okr & rowSums(!is.finite(P)) == 0
    if (sum(okr) <= ncol(P) + 1L) return(FALSE)
    X <- cbind(1, P[okr, , drop = FALSE])
    yd <- response[okr]
    p <- rs$p
    fr <- try(stats::lm.fit(x = X, y = yd), silent = TRUE)
    if (inherits(fr, "try-error")) return(FALSE)
    e <- fr$residuals
    if (centre == 0 && sum((e - mean(e))^2) <= 1e-12 * sum(yd^2)) {
      c0 <- .residual_centre(raw, bounded = FALSE)
      if (c0 != 0) {
        rs0 <- shift_scale(raw - c0)
        f0 <- try(stats::lm.fit(x = X, y = rs0$r[okr]), silent = TRUE)
        if (!inherits(f0, "try-error")) {
          fr <- f0
          yd <- rs0$r[okr]
          p <- rs0$p
        }
      }
    }
    rd <- .rounding_design(X, fr$coefficients, centres)
    isTRUE(.residuals_rounding(fr$residuals, raw[okr] / p, yd, rd$X, rd$b))
  }
  if (rows_rounding()) return(.morans_na())

  cell_pred <- pred_means(predictors)
  cell_resp <- cell_means(response)

  # Fit OLS on cell means
  ok_pred <- apply(is.finite(cell_pred), 1, all)
  ok <- is.finite(cell_resp) & ok_pred
  if (sum(ok) < 4L) return(.morans_na())

  fit <- try(stats::lm.fit(x = cbind(1, cell_pred[ok, , drop = FALSE]),
                            y = cell_resp[ok]),
             silent = TRUE)
  if (inherits(fit, "try-error")) return(.morans_na())
  resid <- fit$residuals

  # Residuals within the rounding of the cell means they are regressed from
  # (.no_spread()) leave nothing to correlate: a constant response, or one
  # the predictors fit exactly.  The bound is relative, so a response in
  # small units is not refused for its unit.
  no_spread <- .no_spread(sum(resid^2), cell_resp[ok])
  # The rounding of cell means: .POOL_ULPS units in the last place of the
  # raw response over the square root of each cell's rows, plus the
  # rounding of the sums the means are formed from, summed over the cells
  # on the residual degrees of freedom, or the rounding of the regression,
  # taken on the size of the terms of the cell design
  # (.residuals_rounding() with the rows of each cell).  The mean of n_j
  # rows carries the rounding of its values over sqrt(n_j), and so does its
  # spread, so the bound of a single row would refuse cell means with real
  # spread where the cells hold many rows.  The sum of n_j terms carries a
  # rounding that grows with them: the bound takes one
  # .Machine$double.eps times the mean absolute term over the cell's rows
  # (of the response as its means are formed, and of each predictor times
  # its coefficient) times sqrt(n_j), on the mean.  A response taken about
  # its centre above lost its spread to the rounding of its raw values, and
  # the residuals of an exact fit are then that rounding, which the bound on
  # the centred cell means alone (.no_spread()) calls spread.
  mean_rounding <- function(fit, cell_resp, ok, p, r) {
    rd <- .rounding_design(cbind(1, cell_pred[ok, , drop = FALSE]), fit$coefficients,
                           centres)
    b <- abs(as.numeric(fit$coefficients[-1L]))
    b[!is.finite(b)] <- 0
    m <- sqrt(cell_rows[ok]) * cell_means(abs(r))[ok]
    if (length(b))
      m <- m + sqrt(cell_all[ok]) * drop(pred_means(abs(predictors))[ok, , drop = FALSE] %*% b)
    .residuals_rounding(fit$residuals, raw[is.finite(raw)] / p, cell_resp[ok],
                        rd$X, rd$b, n = cell_rows[ok], m = .Machine$double.eps * m)
  }
  if (centre != 0 && !no_spread)
    no_spread <- mean_rounding(fit, cell_resp, ok, rs$p, response)
  # A response whose raw values show their spread is not centred above, but
  # its residuals can still be small beside the cell means (1e15 plus 1e5
  # times a predictor plus a spread of 1, or 1e11 plus a spread of 1), and
  # the regression's rounding is then of the order of the residuals times
  # .Machine$double.eps over their ratio to the cell means.  Where their
  # sum of squares is under 1e-12 of the cell means', so that the rounding
  # can reach 1e-10 of them, or where they show no spread at all, the cell
  # means are formed and regressed again about the response's centre
  # (.residual_centre(), as resolution_profile() takes its fit since
  # BUG-36), and refused only if those residuals are within the rounding of
  # the regression, taken on the size of the fit's terms, or of the cell
  # means (mean_rounding(), above).  Residuals do not change under a shift,
  # so this changes nothing but the rounding.
  if (centre == 0 && (no_spread || sum(resid^2) <= 1e-12 * sum(cell_resp[ok]^2))) {
    c0 <- .residual_centre(raw, bounded = FALSE)
    if (c0 != 0) {
      rs0 <- shift_scale(raw - c0)
      cr0 <- cell_means(rs0$r)
      ok0 <- is.finite(cr0) & ok_pred
      fit0 <- if (sum(ok0) >= 4L)
        try(stats::lm.fit(x = cbind(1, cell_pred[ok0, , drop = FALSE]), y = cr0[ok0]),
            silent = TRUE)
      if (!is.null(fit0) && !inherits(fit0, "try-error")) {
        cell_resp <- cr0; ok <- ok0; resid <- fit0$residuals
        no_spread <- mean_rounding(fit0, cell_resp, ok, rs0$p, rs0$r)
      }
    }
  }
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

  W <- .build_knn_weights(cell_xy[ok, , drop = FALSE], k = n_neighbors,
                          use_fnn = use_fnn, use_matrix = use_matrix)

  # Moran's I = (n / S0) * (e' W e) / (e' e)
  S0 <- sum(W)
  if (S0 < .Machine$double.eps || no_spread)
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
  # These are Cliff & Ord's regression-residual moments.  `resid` is by
  # construction the OLS residual of the cell means on cbind(1, cell_pred),
  # the case the formula is derived for, but the derivation also assumes
  # errors of equal variance, and a mean over n_j points has variance
  # sigma^2 / n_j.  Measured under the null at 20 cells, z stayed calibrated
  # on uniform, gradient and moderately clustered layouts (uniform: mean
  # 0.02); next to single-point cells beside cells of 70 or more points its
  # mean was 0.2 to 0.34 and it rejected 7 to 8 percent at the 5 percent
  # level.  Weighting the regression by n_j would remove that; it is not done.
  mom <- .morans_residual_moments(W = W, X = cbind(1, cell_pred[ok, , drop = FALSE]),
                                  S0 = S0, is_sparse = inherits(W, "Matrix"))
  # No usable moments means no usable z -- and z is what the ranking reads.
  # Returning a finite I beside an NA z would be the documented "not
  # computable" shape in one element and a number in the other.
  if (is.null(mom) || !is.finite(mom$VI) || mom$VI <= 0) return(.morans_na())
  c(I = I, z = (I - mom$EI) / sqrt(mom$VI))
}


#' Choose a cell count from the elbow of the k-means WSS curve
#'
#' Fits k-means to the coordinates at \eqn{k = 1} to \code{max_levels} cells
#' and returns the count at the elbow of the within-cluster sum of squares
#' (WSS), read on log-log axes, with its two neighbours.  Points with no
#' cluster structure have no elbow, and the result is then \code{NA} with a
#' warning.  For a count chosen from the response, use
#' \code{\link{resolution_profile}()}.
#'
#' @section How the elbow is read:
#' Points with no cluster structure have a WSS curve close to \eqn{c/k}.  The
#' classical rule, the point furthest below the chord from the first to the
#' last \eqn{k}, still finds a "knee" on such a curve on linear axes, at about
#' \eqn{\sqrt{K_{max}}}: 4 at the default \code{max_levels = 12} and 13 at
#' 160, whatever the data.  On \eqn{\log k} against \eqn{\log} WSS that curve
#' is a straight line, while separated clusters fall faster than it until
#' there is one cell per cluster and like it after, a bend at the cluster
#' count.  The elbow is therefore the \eqn{k} whose \eqn{\log} WSS sags
#' furthest below the straight line joining \eqn{k = 1} and \eqn{K_{max}} on
#' those axes, and it counts as one only when the sag is at least
#' \eqn{\log 1.25} (the WSS a fifth below the power law through the ends).
#'
#' An elongated extent also bends, at about its aspect ratio, because the
#' first cuts go across its long axis; past the threshold that bend is
#' reported as an elbow, and it describes the shape of the extent rather
#' than clusters in it.  When locations repeat (stations visited many
#' times), the sweep can reach one cell per distinct location, where the WSS
#' is zero (to within \eqn{10^{-12}} of the total).  That \eqn{k} has no
#' place on log axes and is left out of the line, but the fall to zero is
#' the sharpest bend there is: when the rest of the curve has no elbow, the
#' number of distinct locations is the elbow.
#'
#' @section When there is no elbow:
#' With no elbow the function returns \code{NA}, so that a count which is
#' not a finding is not passed on as one.  The result carries
#' \code{structured = FALSE} and the largest sag as \code{sag}, and a
#' warning of class \code{spatialkit_no_elbow} says why.  The count the
#' linear-axis rule gives there is the ladder's and not the data's: about
#' \eqn{\sqrt{K_{max}}}, 3 or 4 cells at the default, close to one mean for
#' the whole layer.  The warning carries it as \code{k} and the largest sag
#' as \code{sag}.  \code{on_no_elbow = "warn"} returns that count all the
#' same, with the same warning, and \code{"error"} stops.
#'
#' @section Model-aware criteria:
#' Under \code{criterion = "morans_i"} or \code{"combined"}, with a
#' \code{response_var}, the elbow is supplemented with Moran's I of the cell
#' means at each candidate \eqn{k}: the cell means are regressed on the
#' cell-mean predictors (on an intercept alone when no
#' \code{predictor_vars} are given) and Moran's I is computed on the
#' residuals over each cell's 8 nearest cells.  Under \code{"geometric"},
#' the default, \code{response_var} and \code{predictor_vars} are not read,
#' whether one or both are given.
#'
#' Moran's I of the cell means tests whether neighbouring cells are
#' independent.  It does not measure the autocorrelation the cells leave
#' unexplained in the response, and it does not say how well they represent
#' it: on an autocorrelated response it grows as the cells get finer than
#' the range, because neighbouring cells then resemble each other, and it is
#' close to 0 at about 10 to 20 cells whatever the data, where the neighbour
#' graph is nearly complete.
#'
#' \code{"morans_i"} prefers the candidate with the smallest \eqn{|z|}: the
#' count at which neighbouring cell means stop resembling each other, cells
#' about as wide as the range.  Only the elbow's neighbourhood is scored, so
#' the response moves the answer by at most a few cells either side of the
#' elbow.  That ranking is by a diagnostic of the cell means and not by a
#' count rule, and it is returned with a warning of class
#' \code{spatialkit_diagnostic_ranking} that says so and gives the count
#' ranked first, the first count the statistic could be computed at and the
#' counts scored.  \code{\link{select_resolution}()} refuses the same
#' statistic unless asked with \code{diagnostic = TRUE}; here the criterion
#' is asked for by name, so the count is returned and the warning carries
#' the caveat.
#'
#' \code{"morans_i"} ranks on the standardised deviate
#' \eqn{|z| = |I - E[I]| / \mathrm{sd}(I)}, using the Cliff and Ord
#' regression residual moments, and not on \eqn{|I|}: both \eqn{E[I]} and
#' \eqn{Var[I]} depend on the number of cells, so \eqn{|I|} falls as
#' \eqn{k} grows whether or not the finer cells capture anything, and a
#' ranking by \eqn{|I|} would prefer the largest candidate for arithmetic
#' reasons alone.  The cell-level residuals are OLS residuals by
#' construction, which is the case those moments are derived for, but the
#' derivation also assumes errors of equal variance, and a cell mean over
#' \eqn{n_j} points has a variance proportional to \eqn{1/n_j}, so
#' \eqn{z} is less reliable where single-point cells sit next to cells of
#' many points.  Where structure remains, \eqn{|z|} mixes its size with the
#' number of cells it is measured on, since \eqn{\mathrm{sd}(I)} shrinks as
#' cells are added.  Both quantities are reported in the
#' \code{"diagnostics"} attribute, as \code{moran_i} and \code{moran_z}.
#'
#' Moran's I is computed with an 8-nearest-neighbour weight matrix, so it
#' carries information only once there are more than nine cells.  At nine or
#' fewer, every cell is a neighbour of every other, the row-standardised
#' weight matrix is complete, and Moran's I is exactly \eqn{-1/(k - 1)} for
#' any residual vector.  The residual moments then give \eqn{E[I] = I} and
#' \eqn{\mathrm{Var}[I] = 0}, so the standardised deviate is \eqn{0/0}: those
#' candidates return \code{NA} and are left out of the ranking of
#' \code{"morans_i"}.  When no candidate in the elbow's neighbourhood clears
#' the floor, the call falls back to the geometric ranking and logs a
#' warning that says so.  On points with no cluster structure the elbow
#' sits near \eqn{\sqrt{K_{max}}} and the neighbourhood is the elbow plus or
#' minus \code{max(4, top_n)}, so the neighbourhood reaches ten cells only
#' from about \code{max_levels = 40}.  On such a layer
#' \code{\link{resolution_profile}()}, which scores Moran's z at every level
#' of its ladder, is the model-aware view.
#'
#' \code{"combined"} is the elbow's ranking of the counts around it, best
#' sag first, with Moran's I of the cell means and its \eqn{z} beside it in
#' the diagnostics; Moran's I does not take part in the ranking, and the
#' result is not an estimate of the number of clusters.  The sag is read off
#' the sweep's WSS curve, the one the elbow is read from, so with an elbow
#' the first count is the elbow and the counts after it come in order of
#' sag; they can differ from the two neighbours \code{"geometric"} returns.
#' An elbow below ten cells is a count Moran's I cannot score: the geometric
#' ranking is returned, the elbow and its two neighbours, with a logged
#' warning, and the diagnostics record it (see Value).  When the WSS curve
#' has no elbow (under \code{on_no_elbow = "warn"}) the sag is flat, every
#' candidate tied, and the window comes out in order of distance from the
#' count the chord rule lands on, which is the geometric order.  The k-means
#' re-run at each candidate, whose cells Moran's I is computed on, can land
#' in another optimum than the sweep; its WSS is reported as
#' \code{wss_eval} and does not rank.  Where Moran's I can be computed at no
#' candidate of a window that reaches ten cells (a constant response, too
#' few cells holding a row with one, or too few in the half that
#' \code{select_on = "split"} reads), \code{moran_i} and \code{moran_z} are
#' \code{NA}, the log gives the reason at the WARN level, and the counts are
#' ranked as they would be otherwise.
#'
#' Use \code{"geometric"} when the cell count should follow the clustering
#' of the points, \code{"morans_i"} for a ranking by the diagnostic, and
#' \code{\link{resolution_profile}()} when the response should drive the
#' count.
#'
#' @section The k-means sweep:
#' The WSS curve is read for its shape, so it is fitted to be smooth.  Every
#' point on it is a k-means local optimum, and with a few random restarts
#' the curve mixes level effects with optimisation noise and can rise from
#' one \eqn{k} to the next.  Each \eqn{k} is therefore fitted as the best of
#' 25 restarts seeded by k-means++ (Arthur and Vassilvitskii 2007).  A fixed
#' budget of restarts follows Fränti and Sieranoja (2019), and Steinley (2003)
#' shows why the usual handful is not enough; 25 is the package's choice, not
#' a published figure.  A
#' curve that still rises somewhere is reported and not refused.  A logged
#' warning names the number of rising steps, and the model-aware paths
#' return it as \code{wss_bumps} in the \code{"diagnostics"} attribute beside
#' \code{wss_spread}, the relative spread of WSS across the restarts at each
#' \eqn{k}.
#'
#' The sweep is plain k-means, whatever the cells are built with
#' afterwards.  This function has no \code{partition} argument: the elbow is
#' a reading of the sum of squares that plain k-means minimises, and the
#' cells behind Moran's z are that partition's.  A count returned here says
#' nothing about density-weighted cells; \code{\link{resolution_profile}()},
#' whose cells are of that kind by default, scores counts of those.
#'
#' To keep memory use and run time bounded for large \code{max_levels}, the
#' sweep records only the WSS at each \eqn{k}, without the cluster
#' assignments.  Moran's I is then evaluated on demand: k-means is re-run
#' only for a neighbourhood around the elbow (plus or minus 4 by default, or
#' plus or minus \code{top_n} if larger).  The re-run at each count draws
#' its starts from a stream of its own, seeded before the sweep, so for one
#' \code{set_seed} a count gets the same cells, \code{wss_eval} and Moran's
#' I in any neighbourhood that holds it, whichever \code{top_n} and
#' \code{max_levels} set the neighbourhood: a larger \code{top_n} adds
#' counts to the neighbourhood and leaves the scores of the others as they
#' were.
#'
#' @param data_sf An sf object.  Non-point geometries are reduced to points
#'   with \code{\link{coerce_to_points}()}.  Features without a location are
#'   dropped with a warning, and a layer that mixes XY and XYZ geometries or
#'   carries M is read in x and y (see \link{spatialkit-geometry}).
#' @param max_levels Integer upper bound on levels. Default 12.  The sweep
#'   also stops at the number of distinct locations (k-means cannot place
#'   more centres) and one short of the number of points.
#' @param top_n Integer; how many candidates to return. Default 3. Under
#'   \code{criterion = "geometric"} the candidate set is the elbow and its two
#'   immediate neighbours, so at most 3 values are returned however large
#'   \code{top_n} is; only the model-aware criteria can return more.  Under
#'   those it also sets the counts scored, the elbow plus or minus
#'   \code{max(4, top_n)}: a \code{top_n} above 4 scores more counts, and
#'   under \code{"morans_i"} one of them can rank first.  A count scores the
#'   same whichever \code{top_n} puts it in the window.
#' @param sample_n Integer; subsample size for speed. Default 1500.
#' @param set_seed Integer RNG seed for the subsample and the k-means++
#'   restarts; restored afterwards. Default 123.  The rows are put in
#'   coordinate order before either, so the answer does not depend on the
#'   order they come in.
#' @param response_var Name of the response column (numeric or logical), or
#'   \code{NULL}.  Read under \code{criterion = "morans_i"}, which ranks the
#'   counts around the elbow by Moran's I of the cell means, and
#'   \code{"combined"}, which reports that statistic beside the elbow's
#'   ranking; under \code{"geometric"}, the default, it is not read.
#'   Logicals are read as 0/1; a factor or character response raises an
#'   error and is never coerced, because the residuals of an OLS fit to
#'   arbitrary level codes carry no meaning to test for autocorrelation.
#'   Rows where it, or a predictor, is missing or non-finite stay in the WSS
#'   sweep and the cells and are left out of Moran's I; a logged warning
#'   gives their number.
#' @param predictor_vars Optional predictor column names, read where
#'   \code{response_var} is. Must be numeric or
#'   logical (logicals are read as 0/1); factor/character columns raise an
#'   error.  Without them the model-aware criteria regress the cell means on
#'   an intercept alone, as \code{\link{resolution_profile}()} does.
#' @param criterion One of \code{"geometric"} (default: the elbow and its
#'   two neighbours, read off the coordinates alone), \code{"morans_i"}
#'   (among the candidates around the elbow, the \eqn{k} at which Moran's I
#'   of the cell means, or of their regression residuals, is least
#'   significant: the count at which neighbouring cells stop resembling each
#'   other, which need not be the count that represents the response best),
#'   or \code{"combined"} (the candidates around the elbow in the order of
#'   the WSS curve's log-log sag, the quantity the elbow is read from, with
#'   Moran's I of the cell means and its \eqn{z} in the diagnostics; the
#'   statistic does not take part in the ranking).  Under
#'   \code{"morans_i"} a warning of class
#'   \code{spatialkit_diagnostic_ranking} says that the ranking is by a
#'   diagnostic of the cell means; its condition carries \code{k} (the
#'   count ranked first), \code{criterion}, \code{first_computable} (the
#'   first count at which the statistic is finite) and \code{window}
#'   (the first and last count tried).  \code{"morans_i"} and
#'   \code{"combined"} fall back to \code{"geometric"}, with a warning, when
#'   \code{response_var} is not given, and with a logged warning when no
#'   candidate clears the nine-cell floor described under "Model-aware
#'   criteria".  A \code{response_var} or \code{predictor_vars} naming a
#'   column that is not in \code{data_sf} is an error under every criterion.
#'   Supplying both \code{response_var} and \code{predictor_vars} does not
#'   change the criterion: under \code{"geometric"} neither is read, and a
#'   line logged at the INFO level says that \code{"combined"} reports
#'   Moran's I.
#' @param select_on \code{"all"} (default) selects on every point;
#'   \code{"split"} reads the response on one spatially blocked half of the
#'   points only and returns the other half as the set to estimate on, so
#'   that the standard errors computed downstream on the chosen cells are not
#'   post-selection.  The count is still chosen for the whole layer.  See
#'   "Post-selection inference".
#' @param on_no_elbow What to do when the count is not a finding about the
#'   data: the WSS curve has no elbow, the ladder is too short to show one,
#'   or there are too few points or locations for a ladder at all.  Each
#'   choice raises a condition of class \code{spatialkit_no_elbow} there,
#'   which carries \code{k} (the count the chord rule or the ladder gives),
#'   \code{sag} (the largest log-log sag of the curve; an elbow needs
#'   \eqn{\log 1.25}) and \code{structured = FALSE}, so a caller can catch
#'   it without matching the message.
#'   \code{"na"} (default) returns \code{NA_integer_} there, with a warning
#'   of that class, and every result under \code{"na"} carries the
#'   attributes \code{structured} (\code{TRUE} when the WSS curve has an
#'   elbow, \code{FALSE} with the \code{NA}) and \code{sag}, so one test of
#'   \code{attr(k, "structured")} covers both outcomes.  \code{"warn"}
#'   returns the counts the chord rule or the ladder gives (see Value), with
#'   the same warning; its results carry neither attribute.  \code{"error"}
#'   stops with an error of that class.  Under the model-aware criteria the
#'   candidates are scored around the elbow: with no elbow, \code{"na"} and
#'   \code{"error"} score none, and under \code{"warn"} the ladder places
#'   that window.
#' @section Post-selection inference:
#' When the selection reads the response (here, under
#' \code{criterion = "morans_i"} with a response), everything
#' estimated afterwards on the chosen cells is estimated on data that already
#' influenced the choice, and its standard errors are post-selection ones:
#' descriptive, not at nominal coverage (Gao, Bien and Witten 2022; Chen and
#' Witten 2023 give the exact selective test for two k-means clusters, the
#' first link of this chain).  The exposure is narrower than "the partition
#' was chosen on the response": every partition here is k-means on the
#' coordinates alone, and the response only decides which count is ranked
#' first.  It is not zero, because the count determines every cell the
#' downstream standard errors are computed over.  Under \code{"combined"}
#' the response reaches the diagnostics and not the ranking, which is the
#' same whether or not Moran's I can be computed.
#'
#' \code{select_on = "split"} is sample splitting: the layer is cut into two
#' spatially blocked halves (\code{\link{make_folds}(k = 2, method =
#' "block_kfold")}), the criteria that read the response (Moran's I on the
#' cell means) read the first half only, and the row positions of both
#' halves come back in the \code{"split"} attribute (\code{selection} and
#' \code{estimation}).  The WSS curve and the k-means cells still use every
#' point: they read coordinates alone, and the count is for a tessellation
#' of every point, so it is chosen on that layer's extent and clusters
#' rather than on half of them.  Build the tessellation on every point
#' (cells are geometry), but aggregate and fit on
#' \code{data_sf[attr(x, "split")$estimation, ]}, whose response the
#' selection never saw.  That keeps the selection's use of the response out
#' of the estimates, but only as far as the two halves are independent: the
#' split has no buffer, so points near the border between them are
#' correlated with the selection half over the autocorrelation range, and
#' coverage is nominal only when that range is short against the blocks.
#' The price is precision: half the points estimate, and García Rasines and
#' Young (2023) show a contiguous spatial half is less efficient than the
#' exchangeable split the i.i.d. theory assumes, because the two halves are
#' not interchangeable.  Two alternatives keep the whole sample: data
#' thinning for count responses (Neufeld et al. 2024) and data fission for
#' Gaussian-like ones (Leiner et al. 2023).  Neither is implemented here;
#' the split needs no distributional assumption.  Under \code{"combined"}
#' the split changes only the rows the diagnostics read.  Selection on
#' coordinates alone (\code{"geometric"}, whatever variables are given) is
#' not exposed in this way, and \code{"split"} then changes nothing but the
#' attribute.
#'
#' The estimation rows cover only the estimation half's blocks, while the
#' cells cover the whole layer, so the estimation rows do not fill every
#' cell.  A cell inside the selection half gets none and comes back
#' \code{NA}; a cell across the border between the halves is estimated from
#' its estimation-half points alone, and its standard error describes that
#' part, not the cell.  Count, per cell, how many of its points are
#' estimation rows, and read inferential results only off cells whose
#' points all are.
#' @return An integer vector of candidate counts, best first:
#'   \describe{
#'     \item{\code{"geometric"}}{The elbow, then its lower and upper
#'       neighbours.}
#'     \item{\code{"morans_i"}}{The candidates in the order of \eqn{|z|},
#'       with a warning of class \code{spatialkit_diagnostic_ranking}.}
#'     \item{\code{"combined"}}{The candidates around the elbow in the order
#'       of the sag.}
#'     \item{No elbow}{\code{NA_integer_}, with a warning of class
#'       \code{spatialkit_no_elbow}; under \code{on_no_elbow = "warn"}, the
#'       ladder's choice and its neighbours, with the same warning.  So it
#'       is when the ladder is too short to show an elbow (\code{max_levels}
#'       of 1 or 2), and for a layer with fewer than three points or every
#'       point on one location, which has no ladder (1 under
#'       \code{"warn"}).  A layer with no points is an error.}
#'   }
#'   \code{k[1]} is the top-ranked count, or \code{NA}, on every path, and
#'   \code{top_n = 1} returns it alone.  The vector can carry these
#'   attributes:
#'   \describe{
#'     \item{\code{structured}, \code{sag}}{Under the default
#'       \code{on_no_elbow = "na"}, on every result: \code{structured} is
#'       \code{TRUE} when the WSS curve has an elbow and \code{FALSE} with the
#'       \code{NA}, and \code{sag} is the curve's largest log-log sag.  None
#'       under \code{"warn"} or \code{"error"}.  The geometric path carries no
#'       other attribute.}
#'     \item{\code{diagnostics}}{Under \code{"morans_i"} and
#'       \code{"combined"}: a list with per-\eqn{k} Moran's I
#'       (\code{moran_i}) and standardised deviates (\code{moran_z}), the WSS
#'       curve (\code{wss}) with the relative between-restart spread at each
#'       \eqn{k} (\code{wss_spread}), the number of rising steps on it
#'       (\code{wss_bumps}), the restart budget (\code{nstart}), the elbow
#'       the scored neighbourhood was drawn around (\code{knee_k}), whether
#'       the WSS curve has an elbow (\code{structured}) and its largest
#'       log-log sag (\code{sag}), the \eqn{k} at which k-means failed
#'       (\code{failed_k}; their \code{wss} entries are interpolated from the
#'       neighbours, not measured) and the \eqn{k} the model-aware pass
#'       scored (\code{eval_ks}, the elbow's neighbourhood).  Under
#'       \code{"combined"} it also carries the WSS of the re-run clustering at
#'       those \eqn{k} (\code{wss_eval}, the clustering Moran's I is computed
#'       on; reported, not ranked on), the rank of the sag of \code{wss}
#'       that ordered them (\code{combined_rank}, named by \eqn{k}; with no
#'       elbow every candidate ties and the order is by distance from
#'       \code{knee_k}) and \code{criterion = "combined"}.  When
#'       \code{"combined"} returned the geometric ranking because the elbow is
#'       below ten cells, it carries \code{criterion = "geometric"} and
#'       \code{fallback} (the reason) in place of \code{combined_rank}.
#'       The \code{NA} the default returns on a model-aware call whose WSS
#'       curve has no elbow (or a ladder too short to show one) carries the
#'       attribute too, with what the sweep measured: \code{moran_i} and
#'       \code{moran_z} all \code{NA} and \code{eval_ks} empty, since no
#'       count was scored, \code{knee_k} the count the chord rule lands on,
#'       and \code{structured = FALSE}.  Under \code{"combined"}, Moran's I
#'       computable at no candidate of a neighbourhood that reaches ten cells
#'       leaves the attribute in place, with \code{moran_i} and
#'       \code{moran_z} all \code{NA}, and a logged warning says so.
#'       The attribute is absent when the model-aware path falls back to the
#'       geometric result (no viable \eqn{k} in the elbow's neighbourhood, a
#'       neighbourhood that ends at nine cells or fewer, or, under
#'       \code{"morans_i"}, Moran's I computable at no candidate); these
#'       fallbacks are logged as warnings.  It is absent too on the \code{NA}
#'       of a layer with no ladder (fewer than three points, one location, or
#'       k-means failing at every \eqn{k} from 2).}
#'     \item{\code{split}}{With \code{select_on = "split"}, on every path: a
#'       list with \code{selection} and \code{estimation} (integer row
#'       positions in \code{data_sf}), \code{method} and \code{seed}.}
#'   }
#'   A rising WSS curve is logged on every path.  For a per-level table of
#'   criteria, cell support, restart spread and the flat region, see
#'   \code{\link{resolution_profile}()}.
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
#'   data.frame(x = 5e5 + c(runif(25, 0, 10), runif(25, 90, 100)),
#'              y = 5e6 + c(runif(25, 0, 10), runif(25, 90, 100))),
#'   coords = c("x", "y"), crs = 32632
#' )
#' determine_optimal_levels(pts, max_levels = 6)   # 2 1 3: the elbow first
#'
#' # Points spread evenly have no elbow: NA, with structured = FALSE and a
#' # warning of class spatialkit_no_elbow that says why (printed here
#' # wrapped, in place of the warning)
#' even <- st_as_sf(data.frame(x = 5e5 + runif(200, 0, 100),
#'                             y = 5e6 + runif(200, 0, 100)),
#'                  coords = c("x", "y"), crs = 32632)
#' withCallingHandlers(
#'   determine_optimal_levels(even, max_levels = 6),
#'   spatialkit_no_elbow = function(w) {
#'     writeLines(strwrap(conditionMessage(w), 76))
#'     invokeRestart("muffleWarning")
#'   }
#' )
#' @family resolution
#' @seealso [build_tessellation()] and [get_voronoi_seeds()], which accept
#'   this function's result directly as `approx_n_cells` and `n` (the first
#'   candidate is used, and the output records that it came from here; an
#'   `NA` from it stops them with a note on where it came from);
#'   [assign_features_to_polygons()] and [summarize_by_cell()] for the steps
#'   that follow.
#' @export
determine_optimal_levels <- function(data_sf, max_levels = 12L, top_n = 3L,
                                     sample_n = 1500L, set_seed = 123L,
                                     response_var = NULL,
                                     predictor_vars = NULL,
                                     criterion = c("geometric", "morans_i",
                                                    "combined"),
                                     select_on = c("all", "split"),
                                     on_no_elbow = c("na", "warn", "error")) {
  # A tie of a neighbour search that the rounding of the coordinates can
  # have parted is said once, when the call returns (.tie_watch()).
  tie_end <- .tie_watch("determine_optimal_levels", "data_sf")
  on.exit(tie_end(returnValue(tie_end), data_sf), add = TRUE)
  if (!inherits(data_sf, "sf"))
    stop("determine_optimal_levels(): `data_sf` must be an sf object.")

  criterion <- match.arg(criterion)
  select_on <- match.arg(select_on)
  on_no_elbow <- match.arg(on_no_elbow)
  # A layer with no rows used to return 1, a cell count for nothing.
  if (nrow(data_sf) == 0L)
    stop("determine_optimal_levels(): `data_sf` has no rows, so there is no layer to ",
         "choose a cell count for.", call. = FALSE)
  # A named column that is not there is an error, as it is in
  # resolution_profile().  It used to count as "no model variables": a typo
  # in response_var silently kept the geometric criterion that supplying both
  # variables then upgraded to "combined" (7 6 8 against 5 6 4 on one layer,
  # with no warning and no diagnostics), and under "morans_i" or "combined"
  # the log line blamed variables that had been supplied.
  if (!is.null(response_var)) {
    if (!is.character(response_var) || length(response_var) != 1L || is.na(response_var))
      stop("determine_optimal_levels(): `response_var` must be a single column name.",
           call. = FALSE)
    if (!(response_var %in% names(data_sf)))
      stop(sprintf("determine_optimal_levels(): column '%s' not found.", response_var),
           call. = FALSE)
  }
  missing_p <- setdiff(predictor_vars, names(data_sf))
  if (length(missing_p))
    stop("determine_optimal_levels(): predictor_vars ",
         paste(sQuote(missing_p, FALSE), collapse = ", "), " not found.", call. = FALSE)
  # The criterion the caller names, or the default, stands.  Both model
  # variables used to switch "geometric" to "combined", also when
  # "geometric" was named (match.arg() cannot tell a named default from an
  # omitted one), so the count read the response unasked.
  if (criterion == "geometric" && !is.null(response_var) && !is.null(predictor_vars))
    .log_info(paste0("determine_optimal_levels(): criterion = 'geometric' reads the ",
                     "coordinates alone, so `response_var` and `predictor_vars` are not ",
                     "read; criterion = 'combined' adds Moran's I of the cell means at ",
                     "the counts around the elbow to the result's diagnostics."))
  # Fall back when a model-aware criterion has no response to read.  A
  # warning, not a log line: the caller asked for a criterion and gets another.
  # A response without predictors is enough: .morans_i_for_k() regresses the
  # cell means on an intercept alone, as resolution_profile() has it do, and
  # a caller who asked for "morans_i" with a response used to be told it
  # "needs both" and handed the geometric ranking.
  if (is.null(response_var) && criterion != "geometric") {
    .warn_and_log(paste0("determine_optimal_levels(): criterion = '%s' needs a ",
                         "response_var, and none was given%s; falling back to geometric."),
                  criterion,
                  if (is.null(predictor_vars)) "" else " (predictor_vars alone score nothing)")
    criterion <- "geometric"
  }
  # The model-aware pass runs whenever a model-aware criterion stands, with
  # the predictors given or with none.
  has_model_vars <- criterion != "geometric"
  if (has_model_vars && is.null(predictor_vars)) predictor_vars <- character(0)

  # MULTIPOINT must be coerced too, not merely admitted: st_coordinates()
  # returns one row per VERTEX, so any multi-vertex feature makes xy[i, ] a
  # different feature than row i of resp_vec/pred_mat, and every index below
  # reads the wrong rows silently.  Matches the guard in estimate_sac_range().
  # A layer that mixes XY and XYZ (or XYM) geometries is read in x and y,
  # before the coerce and again after it: sf reads a layer in the dimension
  # of its first geometry.  So is one whose geometries carry M, which GEOS
  # refuses in the coerce (.xy_if_mixed_or_m()).  The layer is not
  # returned, so no warning.  A polygon s2 cannot take a centroid of in
  # degrees is refused in this function's words, naming `data_sf`
  # (.pointize_auto()).
  # A row without a location (an empty point, or a coordinate that is not
  # finite or lies beyond 1e150) is read on the coordinates as given: it is
  # emptied here, before the coerce and the projection, which could
  # otherwise bring it back finite, and the layer is then read as it is read
  # with that row empty.  Points in an sfc_GEOMETRY column are read as the
  # same points in sfc_POINT (.typed_points()): sf::st_coordinates() refuses
  # that class.
  data_sf <- .empty_unlocated(.xy_if_mixed_or_m(.typed_points(data_sf)))
  if (!all(sf::st_geometry_type(data_sf, by_geometry = TRUE) == "POINT")) {
    data_sf <- .xy_if_mixed(.pointize_auto(data_sf, "determine_optimal_levels",
                                           "`data_sf`"))
  }
  data_sf <- ensure_projected(data_sf)

  # A point with no location has no place in a partition of the plane.
  # One POINT EMPTY used to make the k = 1 WSS NA and k-means fail at every
  # k, and the failure handler then returned 1 where the clean layer gave
  # 2 1 3, with nothing but a log line about interpolating the WSS; a point
  # at 1e200 stopped the call with base R's "missing value where TRUE/FALSE
  # needed".  The bound is read again after the projection, which can take a
  # location away (a latitude beyond 90 comes back NaN).
  ok_xy <- .located_rows(data_sf)
  if (!all(ok_xy)) {
    .warn_and_log("determine_optimal_levels(): dropping %d point(s) without a location (%s).",
                  sum(!ok_xy), .UNLOCATED_WHAT)
    data_sf <- data_sf[ok_xy, , drop = FALSE]
  }
  xy <- sf::st_coordinates(data_sf)[, 1:2, drop = FALSE]

  # Sample splitting: read the response on one spatial half, hand the other
  # back for estimation.  Done on the full layer, before the subsample, so
  # the positions returned index `data_sf` as the caller passed it (rows
  # dropped above excepted: they are mapped back).
  split <- NULL
  if (identical(select_on, "split")) {
    split <- .spatial_half_split(data_sf, seed = set_seed,
                                 caller = "determine_optimal_levels")
    sel_rows <- seq_len(nrow(data_sf)) %in% split$selection
    if (!all(ok_xy)) {
      kept <- which(ok_xy)
      split$selection  <- kept[split$selection]
      split$estimation <- kept[split$estimation]
    }
    # The split exists so that a selection made on the response can be
    # estimated on points it never saw, so only the Moran's I pass, which
    # reads the response, is held to the selection half (`moran_rows` below).
    # The WSS sweep and the k-means partitions read coordinates alone and
    # keep every point: the count is for a tessellation of every point.  The
    # whole selection used to run on the half, and a count chosen for the
    # half's extent and clusters (2 where the layer has 4) was then applied
    # to the whole layer.
  }
  # `found` is set once the elbow has been read, under on_no_elbow = "na":
  # the result then carries `structured` and `sag` whichever way it went, so
  # one test of attr(, "structured") covers both.
  found <- NULL
  .with_split <- function(out) {
    if (!is.null(split)) attr(out, "split") <- split
    if (!is.null(found)) {
      attr(out, "structured") <- found$structured
      attr(out, "sag") <- found$sag
    }
    out
  }

  # A count that is not a finding about the data: what the caller is told,
  # and what is returned, under `on_no_elbow`.  The warning carries a class
  # and the numbers: code that runs this over many layers could otherwise
  # tell such a count from an elbow only by matching the message.  It is
  # raised beside the NA too: the NA used to come with a log line alone,
  # which tryCatch(), withCallingHandlers() and options(warn = 2) do not
  # see.  `k` on the condition is the count the chord rule or the ladder
  # gives, the one "warn" returns.  Returns NULL when the caller goes on to
  # return the count ("warn"), or the value to return ("na"), with
  # `diagnostics` attached when there are any.
  .not_a_finding <- function(k, sag, fmt, ..., diagnostics = NULL) {
    msg <- sprintf(fmt, ...)
    fields <- list(k = as.integer(k), sag = as.numeric(sag), structured = FALSE)
    if (identical(on_no_elbow, "error"))
      stop(do.call(errorCondition,
                   c(list(message = msg, class = "spatialkit_no_elbow", call = NULL), fields)))
    if (identical(on_no_elbow, "na")) {
      .warn_classed("spatialkit_no_elbow", "%s Returning NA (on_no_elbow = \"na\").", msg,
                    data = fields)
      out <- structure(NA_integer_, structured = FALSE, sag = as.numeric(sag))
      if (!is.null(diagnostics)) attr(out, "diagnostics") <- diagnostics
      return(.with_split(out))
    }
    .warn_classed("spatialkit_no_elbow", "%s", msg, data = fields)
    NULL
  }

  n  <- nrow(xy)
  if (n < 1L)
    stop("determine_optimal_levels(): no point of `data_sf` has coordinates, so there ",
         "is no layer to choose a cell count for.", call. = FALSE)
  # One or two points used to return 1 with nothing said.
  if (n < 3L) {
    na <- .not_a_finding(1L, NA_real_,
                         paste0("determine_optimal_levels(): %d point%s too few for a ",
                                "ladder of cell counts (it takes three). k = 1 is not a ",
                                "finding about the data."),
                         n, if (n == 1L) " is" else "s are")
    return(if (is.null(na)) .with_split(1L) else na)
  }

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

  # The rows whose response the model-aware criteria may read: the selection
  # half under select_on = "split", and only rows with a finite response and
  # finite predictors.  A cell mean over a row with one missing predictor was
  # NA and dropped the whole cell, and how many cells that removed depended on
  # the points per cell, so z was computed on a different subset of cells at
  # each k: three missing values in 400 rows made every z in the window NA
  # and the call fell back to the geometric ranking.
  moran_rows <- if (!is.null(split) && has_model_vars) sel_rows else rep(TRUE, n)
  if (has_model_vars) {
    complete <- is.finite(resp_vec) & rowSums(!is.finite(pred_mat)) == 0
    if (!all(complete))
      .log_warn(paste0("determine_optimal_levels(): %d of %d row(s) have a missing or ",
                       "non-finite response or predictor; they stay in the WSS sweep ",
                       "and the cells but are left out of Moran's I."),
                sum(!complete), n)
    moran_rows <- moran_rows & complete
  }

  # Subsample from the rows in a canonical order, and hand k-means++ the
  # points in that order even with no subsample to draw: both index rows, so
  # the same layer with its rows permuted used to get different cells and a
  # different count.  Coordinates first; the response and the predictors
  # break ties between repeat visits to one location.
  # First coordinates equal but for their last bits are one column
  # (.near_columns()), so the order does not move with the unit or the origin.
  ord <- do.call(order, c(list(.near_columns(xy[, 1], xy[, 2]), xy[, 2], xy[, 1]),
                          if (has_model_vars)
                            c(list(resp_vec),
                              lapply(seq_len(ncol(pred_mat)), function(j) pred_mat[, j]))))
  idx <- if (n > sample_n) ord[sample(seq_len(n), sample_n)] else ord
  # The distinct locations of the whole layer, for the messages below: on a
  # subsample the ladder is bounded by the subsample's points and locations,
  # and a larger sample_n raises those bounds up to the layer's.
  n_uniq_layer <- if (length(idx) < n) .n_locations(xy) else NA_integer_
  xy <- xy[idx, , drop = FALSE]
  moran_rows <- moran_rows[idx]
  if (has_model_vars) {
    resp_vec <- resp_vec[idx]
    pred_mat <- pred_mat[idx, , drop = FALSE]
  }

  ml    <- as.integer(max_levels)
  k_max <- max(2L, min(ml, nrow(xy) - 1L))

  # One centre per distinct location at most, which with repeat visits is
  # below nrow(xy) - 1 (stats::kmeans() refuses as many centres as points).
  # It was one short of the distinct locations, so five stations visited
  # thirty times each could not reach k = 5.  Counted as stats::kmeans() is
  # handed them (.n_locations()), so the count does not move with the unit.
  n_uniq <- .n_locations(xy)
  k_max <- min(k_max, n_uniq)
  if (k_max < 2L) {
    na <- .not_a_finding(1L, NA_real_,
                         paste0("determine_optimal_levels(): the %d points sit on %d ",
                                "distinct location, so there is no ladder of cell counts. ",
                                "k = 1 is not a finding about the data."),
                         nrow(xy), n_uniq)
    return(if (is.null(na)) .with_split(1L) else na)
  }
  # Which bound ended the ladder, for the messages below: they used to name
  # max_levels whatever did, and on five stations at max_levels = 12 it was
  # the stations.
  # On a subsample the bound is named as the subsample's, with sample_n and
  # the rows of the layer: it used to read "one short of the 9 points" on a
  # layer of 400 rows at sample_n = 9.
  drawn <- nrow(xy) < n
  k_max_from <- if (k_max < max(2L, ml) && k_max == n_uniq)
    (if (drawn && n_uniq < n_uniq_layer)
       sprintf("the %d distinct locations in the sample_n = %d rows drawn from %d",
               n_uniq, nrow(xy), n)
     else sprintf("the %d distinct locations", n_uniq))
  else if (k_max < max(2L, ml))
    (if (drawn) sprintf("one short of the sample_n = %d points drawn from %d", nrow(xy), n)
     else sprintf("one short of the %d points", nrow(xy)))
  else if (ml < 2L) sprintf("max_levels = %d, raised to 2", ml)
  else "max_levels"
  # How far k_max can go with max_levels left out (cap_drawn), and with the
  # whole layer drawn as well (cap_layer).  The advice in the messages below
  # is read from these and not from k_max_from: where max_levels ties the
  # other bound (8 locations at max_levels = 8), k_max_from is "max_levels"
  # and no max_levels takes k_max higher.
  cap_drawn <- min(nrow(xy) - 1L, n_uniq)
  cap_layer <- if (drawn) min(n - 1L, n_uniq_layer) else cap_drawn
  kmeans_cut <- FALSE

  # Say it BEFORE the sweep: the model-aware criteria carry no information at
  # nine cells or fewer (see the resolution floor in Details), so with
  # k_max <= 9 no candidate can ever clear it and the call is going to fall
  # back to the geometric ranking whatever the data say.  The sweep used to
  # run first and the fallback was logged afterwards.
  # Under "na" and "error" a curve with no elbow ends the call before any
  # ranking, so the fallback is promised only for a curve that has one.
  # The line names the bound that ended the ladder (k_max_from), and advises
  # raising max_levels only where a larger max_levels takes k_max past nine
  # cells: on eight distinct locations no max_levels takes k_max above 8.
  # On a subsample a larger sample_n can.
  # "combined" ranks by the sag and reads no response at any count, so it
  # is not told that a larger max_levels lets the criteria contribute: in a
  # window that ends at nine cells or fewer it returns the geometric ranking
  # by a fixed rule.
  if (criterion != "geometric" && k_max <= 9L) {
    bound <- if (startsWith(k_max_from, "the ")) paste(k_max_from, "leave")
             else if (ml < 2L) paste0(k_max_from, ", leaves")
             else paste(k_max_from, "leaves")
    if_elbow <- switch(on_no_elbow,
                       warn = "",
                       na = " if the WSS curve has an elbow (with none the call returns NA)",
                       error = " if the WSS curve has an elbow (with none the call stops)")
    if (criterion == "morans_i")
      .log_warn(paste0("determine_optimal_levels(): %s k_max = %d, ",
                       "and the model-aware criteria carry no information at ",
                       "nine cells or fewer; criterion = '%s' will fall back to ",
                       "the geometric ranking%s. %s"),
                bound, k_max, criterion, if_elbow,
                if (cap_drawn >= 10L)
                  paste0("Raise max_levels above 9 (well above, so the elbow ",
                         "neighbourhood reaches past the floor) for the criteria ",
                         "to contribute.")
                else if (cap_layer >= 10L)
                  sprintf(paste0("Raise sample_n, which draws %d of the %d rows, %sfor ",
                                 "the criteria to contribute."),
                          nrow(xy), n,
                          if (ml <= 9L)
                            paste0("and max_levels above 9 (well above, so the elbow ",
                                   "neighbourhood reaches past the floor) ")
                          else "")
                else sprintf(paste0("No max_levels%s takes k_max above %d on these ",
                                    "points, so the criteria cannot contribute."),
                             if (drawn) " or sample_n" else "", cap_layer))
    else
      .log_warn(paste0("determine_optimal_levels(): %s k_max = %d, so the window ",
                       "around the elbow ends at nine cells or fewer, and there ",
                       "criterion = 'combined' returns the geometric ranking%s. ",
                       "Moran's I, which it reports beside its ranking, is ",
                       "computed from ten cells."),
                bound, k_max, if_elbow)
  }

  # 25 k-means++ restarts per k rather than stats::kmeans(nstart = 5): the
  # WSS curve is read for its shape, and optimisation noise on it is what
  # produced the concave bump the elbow rule once selected.  See
  # .kmeans_best() for the measured effect and the citations.
  # The sweep is plain k-means, the kind "kmeans" of .partition_centres(),
  # from k-means++ starts: the elbow is a reading of the sum of squares that
  # kind minimises.  resolution_profile(partition = "density") is where a
  # density-weighted partition is scored.
  nstart <- 25L
  wss <- numeric(k_max)
  wss_spread <- numeric(k_max)
  failed_k <- integer(0)

  # The integer that seeds the model-aware re-runs (see the loop over
  # `eval_ks` below), read here, before the sweep.  Read after it, it came
  # after as many draws as the sweep took, which grow with k_max, so
  # max_levels 30 and 31 gave one knee and one window and other cells,
  # Moran's I and ranking at the same counts.  Read here it depends on the
  # stream and the distinct locations alone.  The stream is put back after
  # the read, so the sweep starts where it did, and the sweep, the elbow and
  # the geometric ranking are the same under every criterion.  The bound
  # keeps the integer plus any count of the ladder (at most n_uniq) an
  # integer.  With no stream yet (set_seed = NULL on a fresh session) the
  # read starts one, and the sweep goes on from it.
  rerun_seed <- NULL
  if (has_model_vars) {
    stream_before <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    rerun_seed <- sample.int(.Machine$integer.max - n_uniq, 1L)
    if (!is.null(stream_before))
      assign(".Random.seed", stream_before, envir = .GlobalEnv)
  }

  for (k in seq_len(k_max)) {
    if (k == 1L) {
      ctr <- colMeans(xy)
      wss[k] <- sum(rowSums((xy - matrix(ctr, nrow(xy), 2, byrow = TRUE))^2))
    } else {
      # Assignments are discarded on purpose: only the WSS is kept from this
      # sweep, which bounds memory for large max_levels (see Details).
      kb <- .partition_centres(xy, k, partition = "kmeans", init = "kmeans++",
                               nstart = nstart)
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
        k_max_from <- sprintf("k-means failing from k = %d", lo + 1L)
        kmeans_cut <- TRUE
        break
      }
      wss[fk] <- (wss[lo] + wss[hi]) / 2
    }
    if (k_max < 2L) {
      na <- .not_a_finding(1L, NA_real_,
                           paste0("determine_optimal_levels(): k-means failed at every k ",
                                  "from 2, so there is no ladder of cell counts. k = 1 is ",
                                  "not a finding about the data."))
      return(if (is.null(na)) .with_split(1L) else na)
    }
  }

  elbow <- .elbow_from_wss(wss, max_k = k_max, min_k = 1L, return_neighbors = TRUE)
  # No bend on log-log axes: the points have no cluster structure the
  # ladder can see, and the chord rule's answer is set by max_levels.  It
  # used to be handed on as if the data had chosen it, and then with a
  # warning; by default (on_no_elbow = "na") NA is returned in its place, and
  # "warn" returns the count; both raise the classed warning.
  # A ladder of two levels has no line to test, and used to be told it fell
  # in one.
  # The largest sag the curve reaches: what an elbow is measured by, and
  # reported whether or not it is one.
  max_sag <- suppressWarnings(max(elbow$diagnostics$sag, na.rm = TRUE))
  if (!is.finite(max_sag)) max_sag <- NA_real_
  # What the sweep measured, for the NA of a model-aware call: the curve, the
  # chord rule's count and the sag, with no count scored (every Moran's I
  # NA, `eval_ks` empty), in the shape of the diagnostics of a ranked call.
  na_diagnostics <- if (has_model_vars && identical(on_no_elbow, "na"))
    list(moran_i = rep(NA_real_, k_max), moran_z = rep(NA_real_, k_max),
         wss = wss[seq_len(k_max)], wss_spread = wss_spread[seq_len(k_max)],
         wss_bumps = .wss_bumps(wss[seq_len(k_max)]), nstart = nstart,
         knee_k = elbow$knee_k, failed_k = failed_k, eval_ks = integer(0),
         structured = FALSE, sag = max_sag)
  not_found <- NULL
  if (!isTRUE(elbow$structured) && k_max < 3L)
    not_found <- .not_a_finding(elbow$knee_k, max_sag,
                  paste0("determine_optimal_levels(): a ladder of k = 1 to %d (ended ",
                         "by %s) is too short to read an elbow from: it takes three ",
                         "levels to see a bend. k = %d is not a finding about the ",
                         "data; %s"),
                  k_max, k_max_from, elbow$knee_k,
                  # max_levels is the way on only where a larger one takes
                  # k_max to 3; at 2 locations and max_levels = 2 none does.
                  if (kmeans_cut) "choose the count on other grounds."
                  else if (cap_drawn >= 3L) "raise max_levels to 3 or more."
                  else if (cap_layer >= 3L)
                    sprintf("raise sample_n, which draws %d of the %d rows%s.",
                            nrow(xy), n,
                            if (ml < 3L) ", and max_levels to 3 or more" else "")
                  else "choose the count on other grounds.",
                  diagnostics = na_diagnostics)
  else if (!isTRUE(elbow$structured))
    not_found <- .not_a_finding(elbow$knee_k, max_sag,
                  paste0("determine_optimal_levels(): the WSS curve has no elbow: on ",
                         "log-log axes it falls in a straight line, as it does for ",
                         "points with no cluster structure. k = %d is where the chord ",
                         "rule lands on such a curve, set by the ladder (k = 1 to %d, ",
                         "%s) rather than by the data%s. Choose the count on ",
                         "other grounds, e.g. resolution_profile() with a response."),
                  elbow$knee_k, k_max, k_max_from,
                  # Only "warn" goes on to score the window around it; "na"
                  # and "error" end here.
                  if (criterion == "geometric") ""
                  else if (identical(on_no_elbow, "warn"))
                    "; the model-aware criteria are evaluated around it"
                  else "; with no elbow, the model-aware criteria have no window to score",
                  diagnostics = na_diagnostics)
  # on_no_elbow = "na": NA, with the flag and the sag, in place of the count;
  # and the same two attributes on a count the elbow did choose.
  if (!is.null(not_found)) return(not_found)
  if (identical(on_no_elbow, "na")) found <- list(structured = TRUE, sag = max_sag)

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
  # The WSS of the re-run clustering is recorded (wss_eval) beside the
  # Moran's I computed on it: the re-run draws its starts from another
  # stream than the sweep, so it can land in a different local optimum than
  # the sweep did.  No ranking reads it; "combined" ranks by the sweep's
  # sag, the one the elbow is read from.
  #
  # Each re-run has a stream of its own: the re-run at k is seeded with the
  # integer read before the sweep plus k, then the stream is put back.  So a
  # count gets the same cells, wss_eval and Moran's I in every window that
  # holds it, whichever top_n sets the window and whatever max_levels (or
  # anything else that lengthens the sweep) is, and a window that starts
  # lower does not move the counts the narrower window holds.  The integer
  # depends on set_seed, or with set_seed = NULL on the caller's stream,
  # which the read does not move.
  moran_vals <- rep(NA_real_, k_max)   # raw I, reported in $diagnostics
  moran_z    <- rep(NA_real_, k_max)   # standardised deviate, ranks under "morans_i"
  wss_eval   <- wss
  for (k in eval_ks) {
    restore_stream <- .with_seed(rerun_seed + k)
    kb <- .partition_centres(xy, k, partition = "kmeans", init = "kmeans++",
                             nstart = nstart)
    restore_stream()
    if (is.null(kb)) next
    km <- kb$km
    wss_eval[k]   <- km$tot.withinss
    # Cells of the whole layer, means of the rows the selection may read.
    mi            <- .morans_i_for_k(xy[moran_rows, , drop = FALSE], resp_vec[moran_rows],
                                     pred_mat[moran_rows, , drop = FALSE],
                                     km$cluster[moran_rows])
    moran_vals[k] <- mi[["I"]]
    moran_z[k]    <- mi[["z"]]
  }

  # Ranking is on |z|.  |I| is not comparable across k -- see .morans_i_for_k().
  valid_moran <- is.finite(moran_z[eval_ks])

  # Said, with a class, whenever Moran's z orders the candidates returned,
  # which is under criterion = "morans_i" alone.  It is a diagnostic of the
  # cell means and not a count rule (see Details), and select_resolution()
  # refuses the same statistic unless asked with `diagnostic = TRUE`.  Here
  # the caller named the criterion, so the count is returned, and what
  # ranked it is not left to the help page: with a response alone,
  # criterion = "morans_i" returned 10 to 13 cells on 16 of 16 white-noise
  # layers and 16 of 16 exponential fields of 500 points, 10 being the first
  # count the statistic can be computed at.  The way on it names runs on the
  # default profile, where select_resolution(criterion = "cp"), which it used
  # to name, is an error.
  .ranked_by_diagnostic <- function(k_first) {
    scored <- eval_ks[valid_moran]
    .warn_classed(
      "spatialkit_diagnostic_ranking",
      paste0("determine_optimal_levels(): criterion = '%s' ranks the candidates by ",
             "|Moran's z| of %s over each cell's 8 nearest cells. That is a diagnostic ",
             "of whether neighbouring cells are independent, not a count rule: |z| is ",
             "near 0 at about 10 to 20 cells whatever the response and grows as the ",
             "cells get finer than the range, so the smallest |z| does not say how well ",
             "the cells represent the response. k = %d is ranked first; the statistic is ",
             "computable from k = %d, and the window tried was k = %d to %d. Use ",
             "criterion = \"geometric\" for the elbow's ranking, or resolution_profile() ",
             "and select_resolution() for a count chosen for summaries of the response."),
      criterion,
      if (ncol(pred_mat) > 0L)
        "the residuals of the cell means regressed on the cell-mean predictors"
      else "the cell means",
      as.integer(k_first), as.integer(min(scored)),
      as.integer(min(eval_ks)), as.integer(max(eval_ks)),
      data = list(k = as.integer(k_first), criterion = criterion,
                  first_computable = as.integer(min(scored)),
                  window = as.integer(range(eval_ks))))
  }

  if (!any(valid_moran)) {
    # A window that ends at nine cells or fewer cannot score anything (see
    # the resolution floor in Details), and on points with no cluster
    # structure the elbow sits near sqrt(max_levels), so that is the usual
    # case below max_levels = 40: measured on 1000 uniform points,
    # max_levels = 12, 20 and 30 all fell back, and 40 scored k = 10 and 11
    # alone.  The warning used to say only that Moran's I "could not be
    # computed", and the one before the sweep fires only at k_max <= 9.
    # Such a window falls back under either criterion, set by the ladder and
    # not by the response.  "combined" reads no response at any count, so
    # its line says the rule and not a want of information (as before the
    # sweep).
    if (max(eval_ks) <= 9L) {
      if (k_max > 9L && criterion == "combined")
        .log_warn(paste0("determine_optimal_levels(): the window around the elbow, ",
                         "k = %d to %d, ends at nine cells or fewer, and there ",
                         "criterion = 'combined' returns the geometric ranking. ",
                         "Moran's I, which it reports beside its ranking, is computed ",
                         "from ten cells; resolution_profile() scores every level of ",
                         "its ladder."),
                  min(eval_ks), max(eval_ks))
      else if (k_max > 9L)
        .log_warn(paste0("determine_optimal_levels(): the model-aware criteria score ",
                         "only the elbow's neighbourhood, k = %d to %d, and carry no ",
                         "information at nine cells or fewer, so criterion = '%s' ",
                         "falls back to the geometric ranking. On points with no ",
                         "cluster structure the elbow sits near sqrt(max_levels), and ",
                         "the neighbourhood reaches ten cells from about max_levels = ",
                         "40; resolution_profile() scores every level of its ladder."),
                  min(eval_ks), max(eval_ks), criterion)
    } else if (criterion == "morans_i") {
      .log_warn(paste0("determine_optimal_levels(): Moran's I could not be computed at ",
                       "any k from %d to %d in the elbow's neighbourhood (too few of ",
                       "the cells hold a row with a response, the cell-mean ",
                       "regression is singular, or the cell means differ by rounding ",
                       "alone); falling back to geometric."),
                max(10L, min(eval_ks)), max(eval_ks))
    } else {
      # "combined" ranks by the sag, which reads no response, so the ranking
      # goes on with moran_i and moran_z NA.  It fell back to the geometric
      # ranking here, without diagnostics, as it had to while Moran's z took
      # part in its ranking, and the response, or the half select_on =
      # "split" reads, then decided which of the two orders came back:
      # 11 10 12 on a split of eleven clusters whose whole layer gives
      # 11 12 10.
      .log_warn(paste0("determine_optimal_levels(): Moran's I could not be computed at ",
                       "any k from %d to %d in the elbow's neighbourhood (too few of ",
                       "the cells hold a row with a response, the cell-mean ",
                       "regression is singular, or the cell means differ by rounding ",
                       "alone); criterion = 'combined' ranks the ",
                       "counts by the sag as usual, with moran_i and moran_z NA in ",
                       "the diagnostics."),
                max(10L, min(eval_ks)), max(eval_ks))
    }
    if (max(eval_ks) <= 9L || criterion == "morans_i") {
      out <- as.integer(head(elbow$candidates, max(1L, as.integer(top_n))))
      out[out < 1L] <- 1L; out[out > k_max] <- k_max
      return(.with_split(unique(out)))
    }
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
    .ranked_by_diagnostic(out[1L])
    attr(out, "diagnostics") <- list(moran_i = moran_vals, moran_z = moran_z,
                                      wss = wss[1:k_max],
                                      wss_spread = wss_spread[1:k_max],
                                      wss_bumps = wss_bumps, nstart = nstart,
                                      knee_k = knee_k, failed_k = failed_k,
                                      eval_ks = eval_ks,
                                      structured = isTRUE(elbow$structured),
                                      sag = max_sag)
    return(.with_split(out))
  }

  # --- Combined: the elbow's sag ranks; Moran's I is reported beside it ---
  # Only the evaluated neighbourhood is ranked, to keep dimensions aligned.
  #
  # The axis is the log-log sag the elbow itself is read from
  # (.elbow_read()), below the line from k = 1 to k_max, on the sweep's WSS
  # curve, so the elbow ranks first.  It was the chord on linear axes across
  # the window, the rule the elbow stopped using because on a curve like
  # c / k it lands near sqrt(first x last) whatever the data: on eight
  # separated clusters it ranked 6 or 7 first where the elbow is 8.  With no
  # elbow (under on_no_elbow = "warn") the axis is flat, every candidate
  # tied, and the window comes out in order of distance from the chord
  # rule's count, the geometric order.
  #
  # Moran's z of the cell means took part in the ranking, averaged with the
  # sag, and no longer does.  On a curve with no elbow the sag was flat and
  # the counts Moran's z cannot score (nine cells or fewer) took the last
  # place on its axis, so the lowest count it scores, ten, or the count
  # after it came first: on 600 uniform points with a smooth field, 10 6 5
  # where the geometric call returns 6 5 7.  With an elbow, its |z| is near
  # 0 at the coarsest counts it scores whatever the response, and it pulled
  # ten ahead: on eleven clusters with a smooth field, 10 11 12 where the
  # elbow is 11.  moran_i and moran_z stay in the diagnostics.
  #
  # An elbow below ten cells is a count Moran's I cannot score (the nine-cell
  # floor).  There the geometric ranking is returned, the elbow and its two
  # neighbours, and the diagnostics say so.  (When Moran's z still took part,
  # ranking the window put ten first whatever the response did.)
  if (isTRUE(elbow$structured) && knee_k <= 9L) {
    .log_warn(paste0("determine_optimal_levels(): the WSS elbow is at k = %d, ",
                     "below the ten cells Moran's I needs, so criterion = '%s' ",
                     "returns the geometric ranking, with Moran's I in the ",
                     "diagnostics where it could be computed."), knee_k, criterion)
    out <- as.integer(head(elbow$candidates, max(1L, as.integer(top_n))))
    out[out < 1L] <- 1L; out[out > k_max] <- k_max
    out <- unique(out)
    attr(out, "diagnostics") <- list(
      moran_i = moran_vals, moran_z = moran_z, wss = wss[1:k_max],
      wss_eval = wss_eval[1:k_max],
      wss_spread = wss_spread[1:k_max],
      wss_bumps = wss_bumps, nstart = nstart,
      knee_k = knee_k, failed_k = failed_k,
      eval_ks = eval_ks,
      criterion = "geometric",
      fallback = sprintf(paste0("the elbow (k = %d) is below the ten cells ",
                                "Moran's I needs"), knee_k),
      structured = isTRUE(elbow$structured), sag = max_sag
    )
    return(.with_split(out))
  }

  n_eval <- length(eval_ks)
  if (isTRUE(elbow$structured)) {
    # The sweep's sag, the one the elbow was read from, so the elbow (its
    # largest) ranks first.  The re-run's sag ranked here before: where the
    # re-run landed in a worse optimum at the elbow, another count went
    # first and the elbow could drop out of the counts returned (16 17 18
    # on eighteen loose clusters whose elbow is 15).
    sag_w <- elbow$diagnostics$sag[eval_ks]
    rank_elbow <- rank(-sag_w, na.last = TRUE, ties.method = "average")  # more sag = better
  } else {
    rank_elbow <- rep((n_eval + 1) / 2, n_eval)
  }
  # The rank of the sag is the ranking.  Moran's z does not take part (see
  # above); it is still computed and kept in the diagnostics.
  combined_rank <- rank_elbow
  # Exact ties (with a flat axis, every candidate) go to the k nearest the
  # elbow the window was drawn around; order() alone hands them over
  # smallest first, 2 and 3 on uniform points.
  best_idx <- order(combined_rank, abs(eval_ks - knee_k))
  out <- as.integer(eval_ks[head(best_idx, max(1L, as.integer(top_n)))])
  out[out < 1L] <- 1L; out[out > k_max] <- k_max
  out <- unique(out)
  attr(out, "diagnostics") <- list(
    moran_i = moran_vals, moran_z = moran_z, wss = wss[1:k_max],
    wss_eval = wss_eval[1:k_max],
    wss_spread = wss_spread[1:k_max],
    wss_bumps = wss_bumps, nstart = nstart,
    # The geometric elbow the neighbourhood was drawn around, and the k at
    # which k-means failed and whose WSS entries are interpolated -- those
    # values on the curve are not measurements.
    knee_k = knee_k, failed_k = failed_k,
    combined_rank = stats::setNames(combined_rank, eval_ks),
    eval_ks = eval_ks,
    criterion = "combined",
    structured = isTRUE(elbow$structured), sag = max_sag
  )
  .with_split(out)
}
