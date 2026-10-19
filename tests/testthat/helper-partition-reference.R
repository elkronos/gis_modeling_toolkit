# tests/testthat/helper-partition-reference.R
# ---------------------------------------------------------------------------
# Reference implementations of two partition kinds, for
# test-core-partition-reference.R.
#
# local_density(), wkmeanspp(), den_centres() and rnd_centres() are copied
# unchanged from the experiment that tested the kinds
# (local-outputs/core-factorial/10_partitions.R), with the constants they
# read.  They are what the package's density kind and its random starts have
# to reproduce.
#
# The density functions call FNN: a test that uses them skips without it.
# ---------------------------------------------------------------------------

RND_NSTART <- 5L
RND_ITER   <- 100L
DEN_NSTART <- 3L
DEN_ITER   <- 40L
DEN_NN     <- 8L
DEN_POWER  <- 1

# ---- RND ------------------------------------------------------------------

rnd_centres <- function(xy_tr, k, seed) {
  warned <- 0L
  one <- function(alg) {
    set.seed(seed)
    withCallingHandlers(
      tryCatch(stats::kmeans(xy_tr, centers = k, nstart = RND_NSTART, iter.max = RND_ITER, algorithm = alg),
               error = function(e) e),
      warning = function(w) { warned <<- warned + 1L; invokeRestart("muffleWarning") })
  }
  km <- one("Hartigan-Wong")
  fell <- inherits(km, "error")
  if (fell) km <- one("Lloyd")          # only if Hartigan-Wong stops; noted
  if (inherits(km, "error")) stop("k-means failed under both algorithms: ", conditionMessage(km))
  cen <- unname(km$centers)
  attr(cen, "warned") <- warned
  attr(cen, "fell") <- fell
  cen
}

# ---- DEN ------------------------------------------------------------------

# Points per unit area around each point, from the distance to its 8th
# nearest neighbour (not below 1 m).
local_density <- function(xy, k = DEN_NN) {
  k <- min(k, nrow(xy) - 1L)
  d <- FNN::get.knn(xy, k = k)$nn.dist[, k]
  k / (pi * pmax(d, 1)^2)
}

wkmeanspp <- function(xy, k, w) {
  n <- nrow(xy)
  idx <- integer(k)
  idx[1L] <- sample.int(n, 1L, prob = w)
  x1 <- xy[, 1]; x2 <- xy[, 2]
  d2 <- (x1 - x1[idx[1L]])^2 + (x2 - x2[idx[1L]])^2
  if (k > 1L) for (j in 2:k) {
    cs  <- cumsum(w * d2)
    tot <- cs[n]
    idx[j] <- if (is.finite(tot) && tot > 0) min(n, findInterval(stats::runif(1L) * tot, cs) + 1L)
              else sample.int(n, 1L)
    d2 <- pmin(d2, (x1 - x1[idx[j]])^2 + (x2 - x2[idx[j]])^2)
  }
  xy[idx, , drop = FALSE]
}

# Weighted Lloyd: lowers sum_i w_i |x_i - c_a(i)|^2; a centre is the weighted
# mean of its points, and a centre left without a point stays where it is.
den_centres <- function(xy_tr, k, seed) {
  set.seed(seed)
  w <- local_density(xy_tr)^DEN_POWER
  w <- w / mean(w)
  best <- NULL
  for (r in seq_len(DEN_NSTART)) {
    cen <- wkmeanspp(xy_tr, k, w)
    for (it in seq_len(DEN_ITER)) {
      a  <- FNN::get.knnx(cen, xy_tr, k = 1L)$nn.index[, 1]
      sw <- rowsum(w, a, reorder = TRUE)
      sx <- rowsum(xy_tr * w, a, reorder = TRUE)
      ids <- as.integer(rownames(sw))
      new <- cen
      new[ids, ] <- sx / as.numeric(sw)
      moved <- max(abs(new - cen))
      cen <- new
      if (moved < 1e-6) break
    }
    nn  <- FNN::get.knnx(cen, xy_tr, k = 1L)
    obj <- sum(w * nn$nn.dist[, 1]^2)
    if (is.null(best) || obj < best$obj) best <- list(cen = cen, obj = obj)
  }
  unname(best$cen)
}

# ---- Point sets and measures ----------------------------------------------

# A clustered layout: `n_bg` points spread over a square of side `side` and
# `n_cl` more in Gaussian clusters, so that the point density varies by a
# factor of a hundred or so across the layer.  The rows come out shuffled.
# The caller's random-number stream is left as it was.
partition_clustered_xy <- function(n_bg = 300L, n_cl = 900L, clusters = 4L, side = 1000,
                                   origin = c(0, 0), seed = 11L) {
  cleanup <- spatialkit:::.with_seed(seed)
  on.exit(cleanup(), add = TRUE)
  bg  <- cbind(stats::runif(n_bg, 0, side), stats::runif(n_bg, 0, side))
  ctr <- cbind(stats::runif(clusters, 0.2, 0.8) * side, stats::runif(clusters, 0.2, 0.8) * side)
  who <- sample.int(clusters, n_cl, replace = TRUE)
  cl  <- ctr[who, , drop = FALSE] + matrix(stats::rnorm(2L * n_cl, sd = side / 40), n_cl, 2L)
  xy  <- rbind(bg, cl)
  xy  <- xy[sample.int(nrow(xy)), , drop = FALSE]
  unname(sweep(xy, 2L, origin, "+"))
}

# A layout with more dense places than a coarse partition has cells, of two
# sizes: `cities` large Gaussian clusters hold half of the points (in
# decreasing shares), `towns` small tight ones hold 30 percent, and the rest
# are spread over the square.  The rows come out shuffled, and the caller's
# random-number stream is left as it was.
partition_towns_xy <- function(n = 2000L, cities = 3L, towns = 25L, side = 1000, seed = 1L) {
  cleanup <- spatialkit:::.with_seed(seed)
  on.exit(cleanup(), add = TRUE)
  n_city <- round(0.5 * n); n_town <- round(0.3 * n); n_bg <- n - n_city - n_town
  cc   <- cbind(stats::runif(cities, 0.2, 0.8), stats::runif(cities, 0.2, 0.8)) * side
  wc   <- sample.int(cities, n_city, replace = TRUE, prob = cities:1)
  city <- cc[wc, , drop = FALSE] + matrix(stats::rnorm(2L * n_city, sd = side / 40), n_city, 2L)
  tc   <- cbind(stats::runif(towns), stats::runif(towns)) * side
  wt   <- sample.int(towns, n_town, replace = TRUE)
  town <- tc[wt, , drop = FALSE] + matrix(stats::rnorm(2L * n_town, sd = side / 250), n_town, 2L)
  xy   <- rbind(city, town, cbind(stats::runif(n_bg), stats::runif(n_bg)) * side)
  unname(xy[sample.int(nrow(xy)), , drop = FALSE])
}

# `n` points spread evenly over a square of side `side`.  The caller's
# random-number stream is left as it was.
partition_uniform_xy <- function(n = 1200L, side = 1000, seed = 41L) {
  cleanup <- spatialkit:::.with_seed(seed)
  on.exit(cleanup(), add = TRUE)
  cbind(stats::runif(n), stats::runif(n)) * side
}

# A coordinate matrix as an sf POINT layer, with any further columns.
partition_points_sf <- function(xy, crs = 32632, ...) {
  d <- data.frame(x = xy[, 1L], y = xy[, 2L], ...)
  sf::st_as_sf(d, coords = c("x", "y"), crs = crs)
}

# The coordinates of an sf POINT layer as a bare two-column matrix.
partition_xy <- function(x) unname(sf::st_coordinates(x)[, 1:2, drop = FALSE])

# The rows in coordinate order: the order the package's density kind works in.
partition_sorted <- function(xy) xy[order(xy[, 1L], xy[, 2L]), , drop = FALSE]

# The nearest centre of every point, by brute force; the lower centre number
# wins a tie.
partition_nearest <- function(xy, cen) {
  d2 <- outer(xy[, 1L], cen[, 1L], "-")^2 + outer(xy[, 2L], cen[, 2L], "-")^2
  apply(d2, 1L, which.min)
}

# The coefficient of variation of the cell counts: 0 when every cell holds
# the same number of points.
partition_count_cv <- function(sizes) stats::sd(sizes) / mean(sizes)

# Run a call whose warnings are not what the test is about.
partition_quietly <- function(expr) suppressWarnings(suppressMessages(expr))

# Run a call and keep the warnings it raises: list(value, warnings), the
# warnings as condition objects.
partition_warnings <- function(expr) {
  got <- list()
  value <- withCallingHandlers(expr, warning = function(w) {
    got[[length(got) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = got)
}

# Those of a partition_warnings() result that have a class.
partition_of_class <- function(res, class)
  Filter(function(w) inherits(w, class), res$warnings)
