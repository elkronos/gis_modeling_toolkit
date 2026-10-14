# tests/testthat/helper-pool-reference.R
# ---------------------------------------------------------------------------
# Reference implementations of the region values, their leave-one-out errors
# and the size rule, for the tests of R/pool.R.
#
# Everything above the line "Layouts" is copied unchanged, comments and all,
# from the experiments that tested the estimates, with the constants the
# functions read.  It is what the package's pooled value, its own-data
# weight, its two leave-one-out errors and its count at a size have to
# reproduce.  The comments are the experiments' own: a "tract" there is a
# row, and a "training" row is a row with a value.
#
#   local-outputs/core-factorial/20_estimates.R
#     EST, NEIGHBOURS, PAR_DIV, PAR_NSTART, PAR_SEED, oneway(), reml_oneway(),
#     weight(), reldiff(), parent_of(), region_estimates(),
#     predict_partition()
#   local-outputs/core-factorial/30_count_rules.R
#     loo_plain(), pick_rung()
#   local-outputs/core-factorial/10_partitions.R
#     rung_k()
#   local-outputs/estimate-design/20_lib.R
#     wfun(), reml_means(), fit_regions(), loo_int()
#
# region_estimates(want = "L") is the pooled value and want = "P" the plain
# mean; fit_regions() returns the parts of the pooled value (sigma2, tau2, w,
# t, Q, nb) and the leave-one-out residuals of both estimates (eL, eP).
#
# parent_of() calls kpp_centres(), which is not copied: it runs only for the
# estimate "PAR", and fit_regions() asks for "P", "R" and "L".
#
# The functions call FNN: a test that uses them skips without it.
# ---------------------------------------------------------------------------

# ---- core-factorial/20_estimates.R ---------------------------------------------

EST        <- c("P", "R", "L", "PAR")
NEIGHBOURS <- 8L
PAR_DIV    <- 6
PAR_NSTART <- 5L
PAR_SEED   <- 123L

# ---- copied from the prototype -----------------------------------------------------

oneway <- function(y, gi, K) {
  n  <- tabulate(gi, K)
  stopifnot(all(n > 0L))
  S  <- as.numeric(rowsum(y, gi))
  m  <- S / n
  N  <- length(y)
  mu <- sum(y) / N
  if (N <= K) stop("no group holds two training tracts: no within-group variance")
  sigma2 <- sum((y - m[gi])^2) / (N - K)
  raw <- if (K > 1L) {
    (sum(n * (m - mu)^2) / (K - 1L) - sigma2) / ((N - sum(n^2) / N) / (K - 1L))
  } else NA_real_
  tau2 <- if (K > 1L) max(0, raw) else 0
  list(n = n, S = S, m = m, N = N, K = K, mu = mu, sigma2 = sigma2, tau2 = tau2, tau2_raw = raw)
}

reml_oneway <- function(a) {
  ssw <- a$sigma2 * (a$N - a$K)
  at <- function(r) {
    w  <- a$n / (1 + a$n * r)
    mu <- sum(w * a$m) / sum(w)
    s2 <- (ssw + sum(w * (a$m - mu)^2)) / (a$N - 1)
    list(dev = (a$N - 1) * log(s2) + sum(log1p(a$n * r)) + log(sum(w)), mu = mu, sigma2 = s2,
         tau2 = r * s2)
  }
  o  <- stats::optimize(function(lr) at(exp(lr))$dev, c(-25, 15), tol = 1e-10)
  b  <- at(exp(o$minimum))
  z  <- at(0)
  if (z$dev <= b$dev) z else b
}

weight <- function(tau2, sigma2, n) {
  v <- sigma2 / n
  ifelse(tau2 + v > 0, tau2 / (tau2 + v), 1)
}

reldiff <- function(a, b) ifelse(a == b, 0, abs(a - b) / pmax(abs(a), abs(b)))

# ---- the four estimates of every region --------------------------------------------

# The parent of each region: k-means of the regions' centroids at Kp parents.
parent_of <- function(cen, Kp, seed = PAR_SEED) {
  K <- nrow(cen)
  if (Kp < 2L || K < 3L) return(rep(1L, K))
  Kp <- min(Kp, K - 1L)
  pc <- kpp_centres(cen, Kp, seed, nstart = PAR_NSTART)
  p  <- FNN::get.knnx(pc, cen, k = 1L)$nn.index[, 1]
  match(p, sort(unique(p)))
}

# Training values `y` in regions `gi` (1..K, none empty), training
# coordinates `xy`.  Returns a K x length(want) matrix of estimates.
region_estimates <- function(y, gi, K, xy, want = EST) {
  a <- oneway(y, gi, K)
  out  <- matrix(NA_real_, K, length(want), dimnames = list(NULL, want))
  secs <- stats::setNames(numeric(length(want)), want)
  cen  <- if (any(c("L", "PAR") %in% want)) rowsum(xy, gi) / a$n else NULL
  for (e in want) {
    t0 <- proc.time()[[3]]
    if (e == "P") {
      out[, e] <- a$m
    } else if (e == "R") {
      r  <- reml_oneway(a)
      wr <- weight(r$tau2, r$sigma2, a$n)
      out[, e] <- wr * a$m + (1 - wr) * r$mu
    } else if (e == "L") {
      if (K > 1L) {
        nb   <- FNN::get.knn(cen, k = min(NEIGHBOURS, K - 1L))$nn.index
        t3   <- rowSums(matrix(a$S[nb], K)) / rowSums(matrix(a$n[nb], K))
        tau3 <- max(0, (sum(a$n * (a$m - t3)^2) - K * a$sigma2) / a$N)
        w3   <- weight(tau3, a$sigma2, a$n)
        out[, e] <- w3 * a$m + (1 - w3) * t3
      } else out[, e] <- a$m
    } else if (e == "PAR") {
      parent <- parent_of(cen, as.integer(round(K / PAR_DIV)))
      P   <- max(parent)
      np  <- as.numeric(rowsum(a$n, parent))
      pm  <- as.numeric(rowsum(a$S, parent)) / np
      ssb <- sum(a$n * (a$m - pm[parent])^2)
      dfb <- K - P
      s2b <- 0
      if (dfb > 0L) {
        n0  <- (a$N - sum(as.numeric(rowsum(a$n^2, parent)) / np)) / dfb
        s2b <- max(0, (ssb / dfb - a$sigma2) / n0)
      }
      B <- weight(s2b, a$sigma2, a$n)
      if (s2b > 0) {
        w  <- 1 / (s2b + a$sigma2 / a$n)
        pm <- as.numeric(rowsum(w * a$m, parent)) / as.numeric(rowsum(w, parent))
      } else B <- rep(0, K)
      out[, e] <- pm[parent] + B * (a$m - pm[parent])
    } else stop("unknown estimate ", e)
    secs[e] <- proc.time()[[3]] - t0
  }
  list(est = out, a = a, secs = secs)
}

# One fold, one grouping: g_tr and g_te are the group ids of the training and
# the held-out tracts (any type).  `p_empty`: the study's saved prediction,
# used under P for a held-out tract of an empty group (NULL: the rule itself).
predict_partition <- function(g_tr, g_te, y_tr, xy_tr, xy_te, want = EST, p_empty = NULL) {
  u  <- sort(unique(g_tr))
  K  <- length(u)
  gi <- match(g_tr, u)
  ti <- match(g_te, u)
  empty <- is.na(ti)
  fi <- ti
  if (any(empty)) {
    fi[empty] <- gi[FNN::get.knnx(xy_tr, xy_te[empty, , drop = FALSE], k = 1L)$nn.index[, 1]]
  }
  r <- region_estimates(y_tr, gi, K, xy_tr, want)
  pred <- r$est[fi, , drop = FALSE]
  fill_gap <- 0
  if (any(empty) && !is.null(p_empty) && "P" %in% want) {
    fill_gap <- max(reldiff(pred[empty, "P"], p_empty[empty]))
    pred[empty, "P"] <- p_empty[empty]
  }
  list(pred = pred, n_own = ifelse(empty, 0L, r$a$n[ti]), K = K, a = r$a, gi = gi, secs = r$secs,
       fill_gap = fill_gap, n_empty = sum(empty))
}

# ---- core-factorial/30_count_rules.R -------------------------------------------

# Mean squared closed-form leave-one-out error of plain averages.  `a` is
# oneway() of the training values in the regions `gi`; `nn1` the index of
# each training tract's nearest other training tract.
loo_plain <- function(y_tr, gi, a, nn1) {
  nj <- a$n[gi]
  lr <- (y_tr - a$m[gi]) * nj / pmax(nj - 1L, 1L)
  s  <- nj == 1L
  if (any(s)) lr[s] <- y_tr[s] - a$m[gi[nn1[s]]]
  mean(lr^2)
}

# The rung a criterion chooses (lower is better; the first of equals).
pick_rung <- function(crit) which.min(crit)

# ---- core-factorial/10_partitions.R --------------------------------------------

# The count of a rung: m training tracts per region.
rung_k <- function(n_train, m) max(2L, as.integer(round(n_train / m)))

# ---- estimate-design/20_lib.R --------------------------------------------------

wfun <- function(tau2, v) ifelse(tau2 + v > 0, tau2 / (tau2 + v), 1)

# REML of tau2 for K means with known variances v (see the header, Rv).
reml_means <- function(m, v) {
  at <- function(t2) {
    w  <- 1 / (t2 + v)
    mu <- sum(w * m) / sum(w)
    list(dev = sum(log(t2 + v)) + log(sum(w)) + sum(w * (m - mu)^2), mu = mu, tau2 = t2)
  }
  s <- stats::var(m)
  if (!is.finite(s) || s <= 0) return(at(0))
  o <- stats::optimize(function(lt) at(exp(lt))$dev, c(log(s) - 25, log(s) + 5), tol = 1e-10)
  b <- at(exp(o$minimum))
  z <- at(0)
  if (z$dev <= b$dev) z else b
}

# Everything about one fold and one grouping.  g_tr, g_te: the group of every
# training and held-out tract.  nn1: each training tract's nearest other
# training tract.  s2_tr: the training tracts' sampling variances (NULL: the
# dataset has none).  p_empty: the study's saved prediction, used under P for
# a held-out tract of an empty group (NULL: the fill rule itself).
fit_regions <- function(g_tr, g_te, y_tr, xy_tr, xy_te, nn1, s2_tr = NULL, p_empty = NULL) {
  stopifnot(all(is.finite(y_tr)), is.null(s2_tr) || all(is.finite(s2_tr) & s2_tr > 0))
  u  <- sort(unique(g_tr))
  K  <- length(u)
  gi <- match(g_tr, u)
  ti <- match(g_te, u)
  empty <- is.na(ti)
  fi <- ti
  if (any(empty)) fi[empty] <- gi[FNN::get.knnx(xy_tr, xy_te[empty, , drop = FALSE], k = 1L)$nn.index[, 1]]
  r  <- region_estimates(y_tr, gi, K, xy_tr, want = c("P", "R", "L"))      # the factorial's P, R, L
  a  <- r$a
  n  <- a$n; N <- a$N
  # L's parts once more (the expressions of region_estimates()), kept for what follows
  cen <- rowsum(xy_tr, gi) / n
  nb  <- FNN::get.knn(cen, k = min(NEIGHBOURS, K - 1L))$nn.index
  Sn  <- rowSums(matrix(a$S[nb], K))
  Nn  <- rowSums(matrix(n[nb], K))
  t3  <- Sn / Nn
  Q   <- sum(n * (a$m - t3)^2)
  tau2 <- max(0, (Q - K * a$sigma2) / N)
  w   <- weight(tau2, a$sigma2, n)
  L   <- w * a$m + (1 - w) * t3
  stopifnot(identical(L, unname(r$est[, "L"])))
  # leave one training tract out
  nj  <- n[gi]
  n1  <- nj - 1L
  d1  <- pmax(n1, 1L)
  res <- y_tr - a$m[gi]
  ssw <- as.numeric(rowsum(res^2, gi))
  sng <- nj == 1L
  kP  <- gi
  kP[sng] <- gi[nn1[sng]]                          # the region whose mean predicts a tract alone in its own
  eP  <- res * nj / d1
  eP[sng] <- y_tr[sng] - a$m[kP[sng]]
  m1  <- (a$S[gi] - y_tr) / d1
  w1  <- ifelse(n1 >= 1L, weight(tau2, a$sigma2, d1), 0)
  eL  <- y_tr - (w1 * m1 + (1 - w1) * t3[gi])
  est <- cbind(P = a$m, R = unname(r$est[, "R"]), L = L, Pv = NA_real_, Rv = NA_real_, Lv = NA_real_)
  out <- list(K = K, gi = gi, ti = ti, fi = fi, empty = empty, n = n, S = a$S, m = a$m, N = N, sigma2 = a$sigma2,
              ssw = ssw, nb = nb, Sn = Sn, Nn = Nn, t = t3, Q = Q, tau2 = tau2, w = w, nj = nj, n1 = n1, kP = kP,
              eP = eP, eL = eL, w1 = w1, loo_P = loo_plain(y_tr, gi, a, nn1), loo_L = mean(eL^2),
              loo_L_se = stats::sd(eL^2) / sqrt(N), has_v = !is.null(s2_tr))
  if (!is.null(s2_tr)) {
    A     <- sum(s2_tr * (1 - 1 / nj))
    su2   <- max(0, (sum(ssw) - A) / (N - K))
    om    <- 1 / (su2 + s2_tr)
    s2sum <- as.numeric(rowsum(s2_tr, gi))
    sbar  <- s2sum / n
    v     <- (su2 + sbar) / n
    tauv  <- max(0, (Q - sum(n * v)) / N)
    wv    <- wfun(tauv, v)
    rm    <- reml_means(a$m, v)
    Bv    <- wfun(rm$tau2, v)
    est[, "Pv"] <- as.numeric(rowsum(om * y_tr, gi)) / as.numeric(rowsum(om, gi))
    est[, "Lv"] <- wv * a$m + (1 - wv) * t3
    est[, "Rv"] <- Bv * a$m + (1 - Bv) * rm$mu
    sb1  <- (s2sum[gi] - s2_tr) / d1
    wv1  <- ifelse(n1 >= 1L, wfun(tauv, (su2 + sb1) / d1), 0)
    out  <- c(out, list(su2 = su2, sbar = sbar, v = v, tauv = tauv, wv = wv, taur = rm$tau2, wv1 = wv1,
                        v1 = (su2 + sb1) / d1,                 # the plain mean's variance with tract i left out (repair)
                        eLv = y_tr - (wv1 * m1 + (1 - wv1) * t3[gi])))
  }
  out$est <- est
  pred <- est[fi, , drop = FALSE]
  if (any(empty) && !is.null(p_empty)) pred[empty, "P"] <- p_empty[empty]
  out$pred  <- pred
  out$n_own <- ifelse(empty, 0L, n[fi])
  out
}

loo_int <- function(centre, r_tr, c_tr, c_te) {
  q <- stats::quantile(r_tr / c_tr, c(0.05, 0.95), names = FALSE, type = 7)
  cbind(centre + c_te * q[1], centre + c_te * q[2])
}

# ---- Layouts -----------------------------------------------------------------------

# The layouts the engine is checked on.  Ten by size, as (rows, rows per
# region, clustered), and three with the number of regions set: 2 and 3
# regions on 60 rows, 4 on 90 clustered rows.
POOL_LAYOUTS <- list(
  list(N = 600L,  per = 18,  clustered = FALSE), list(N = 600L,  per = 4.5, clustered = FALSE),
  list(N = 600L,  per = 2,   clustered = TRUE),  list(N = 3000L, per = 4.5, clustered = TRUE),
  list(N = 3000L, per = 9,   clustered = FALSE), list(N = 3000L, per = 3,   clustered = TRUE),
  list(N = 3000L, per = 2,   clustered = FALSE), list(N = 40L,   per = 4.5, clustered = FALSE),
  list(N = 12L,   per = 4.5, clustered = FALSE), list(N = 9L,    per = 4.5, clustered = TRUE),
  list(N = 60L, clustered = FALSE, K = 2L, seed = 5L),
  list(N = 60L, clustered = FALSE, K = 3L, seed = 6L),
  list(N = 90L, clustered = TRUE,  K = 4L, seed = 7L))

# Seed the generator R starts with, whatever kind the session has set, and
# return the function that puts the session's kind and stream back.
pool_seed <- function(seed) {
  had  <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old  <- if (had) get(".Random.seed", envir = globalenv())
  kind <- RNGkind()
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  function() {
    suppressWarnings(do.call(RNGkind, as.list(kind)))
    if (had) assign(".Random.seed", old, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
      rm(".Random.seed", envir = globalenv())
  }
}

# One layout: `N` locations on a 100 km square, spread evenly or in twelve
# Gaussian clusters, a response with a smooth trend and noise, and the
# regions of a k-means fit at rung_k(N, per) centres (or `K`): each row is in
# the region of its nearest centre, and the regions that hold a row are
# numbered 1 to K.  The seed is the one the layouts were first run with,
# 100 + N + 7 * per.  The draws are R's default generator's whatever kind
# the session has set, and the caller's random-number stream is left as it
# was.
pool_layout <- function(N, per = NA, clustered = FALSE, K = NULL,
                        seed = 100 + N + 7 * per) {
  cleanup <- pool_seed(seed)
  on.exit(cleanup(), add = TRUE)
  xy <- if (clustered) {
    c0 <- cbind(stats::runif(12, 0, 1e5), stats::runif(12, 0, 1e5))
    i  <- sample(12, N, replace = TRUE)
    cbind(c0[i, 1] + stats::rnorm(N, 0, 4000), c0[i, 2] + stats::rnorm(N, 0, 4000))
  } else cbind(stats::runif(N, 0, 1e5), stats::runif(N, 0, 1e5))
  y <- sin(xy[, 1] / 2e4) + cos(xy[, 2] / 1.5e4) + stats::rnorm(N, 0, 0.6)
  if (is.null(K)) K <- rung_k(N, per)
  km <- suppressWarnings(stats::kmeans(xy, centers = K, nstart = 2, iter.max = 100))
  a  <- spatialkit:::.nearest_centre(xy[, 1], xy[, 2], km$centers[, 1], km$centers[, 2])
  gi <- match(a, sort(unique(a)))
  list(xy = xy, y = y, gi = gi, K = max(gi))
}

# Each row's nearest other row, as the reference functions take it: the
# experiments' own expression.  It serves layouts without repeated locations
# only.  Where a location repeats, FNN can return a row as its own
# neighbour, and which of several rows at one distance it returns is not
# fixed; such a layout is refused.
pool_nearest_other <- function(xy) {
  if (anyDuplicated(xy))
    stop("pool_nearest_other(): two rows share a location; the reference serves ",
         "layouts without repeated locations.")
  nn <- FNN::get.knn(xy, k = 1L)$nn.index[, 1]
  if (any(nn == seq_len(nrow(xy))))
    stop("pool_nearest_other(): a row came back as its own neighbour.")
  nn
}
