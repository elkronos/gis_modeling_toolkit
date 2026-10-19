# =============================================================================
# Region values pooled toward their neighbours: the engine.
#
# For one column, with its rows in regions 1..K.  A region's plain mean is
# m_j = S_j / n_j.  Its pooled value is
#
#   L_j = w_j m_j + (1 - w_j) t_j
#
#   t_j     the target: the mean of the rows of the 8 regions whose centroids
#           are nearest region j's (all the other regions when there are
#           fewer than 9), region j itself never among them
#   sigma2  the spread of the rows about their region's mean,
#           sum (y_i - m_g(i))^2 / (N - K)
#   tau2    the mean squared distance between a region's level and its
#           target, max(0, (Q - K sigma2) / N) with Q = sum n_j (m_j - t_j)^2;
#           the target's own error is part of it
#   w_j     the own weight, tau2 / (tau2 + sigma2 / n_j): it grows with the
#           region's count
#
# Both estimates have a leave-one-out error: each row predicted by the value
# of its region with that row taken out, sigma2, tau2, the centroids, the
# neighbours and the targets held.
#
# .pool_fit() is the one place this arithmetic is written.  Its expressions
# are held, bit for bit, to the reference functions in
# tests/testthat/helper-pool-reference.R, so an edit that reorders a sum is a
# change of result.  .pool_column() prepares one column of an assigned layer
# for it.  Everything in this file is internal.
#
# Two things the result must not depend on, and how each is kept out:
#
#   the order of the rows   The centroids are summed with the rows in
#                           coordinate order, and a row's nearest other row
#                           is found through the distinct locations in that
#                           order.  Sums of the response are taken in the
#                           rows' own order, which moves a result in its last
#                           bits and no further.
#   whether FNN is there    Every search has an exact path in base R, and the
#                           kd-tree path returns the same answer to the last
#                           bit: a row whose candidates tie is searched again
#                           by the exact rule.
# =============================================================================

# Regions in a target.
.POOL_NEIGHBOURS <- 8L

# Rows per region: the size a count is set by when none is given, and the
# sizes of a ladder of sizes, coarsest first.
.OBS_PER_CELL_DEFAULT <- 4.5
.OBS_PER_CELL_LADDER  <- c(18, 12, 9, 6, 4.5, 3, 2)

# Relative difference up to which two leave-one-out errors, computed at two
# ends of the pipeline on the same rows and regions, are the same number.
.RECORD_TOL <- 1e-8

# An interval's ends are called unstable when fewer residuals than this are
# expected beyond each of them.
.FEW_RESIDUALS <- 5

# sigma2 is called unstable when it rests on fewer values beyond one per
# region than this.
.POOL_MIN_DF <- 30L

# Regions beat one mean for all rows when the one mean's leave-one-out error
# is above theirs by more than this fraction.
.STRUCTURE_TOL <- 0.02

# The most pairs an exact search compares, unless
# options(spatialkit.exact_pairs_max = ) sets another number
# (.exact_pairs_max()), and the most point-to-centre comparisons of a
# partition fit at a fixed size.
.EXACT_PAIRS_MAX <- 2e9
.FIT_WORK_MAX    <- 1e10


# -----------------------------------------------------------------------------
# Backends and the bound on exact searches
# -----------------------------------------------------------------------------

#' Whether an optional package can be loaded
#'
#' How the code in this file asks whether \pkg{FNN} is there.  The answer is
#' passed down as \code{use_fnn}, so replacing this one function switches the
#' backend off for a whole call.
#'
#' @param pkg Package name.
#' @return \code{TRUE} or \code{FALSE}.
#' @keywords internal
#' @noRd
.has_backend <- function(pkg) requireNamespace(pkg, quietly = TRUE)


#' The most pairs an exact search may compare
#'
#' @return \code{getOption("spatialkit.exact_pairs_max")}, or
#'   \code{.EXACT_PAIRS_MAX} when the option is not set.  A value that is not
#'   one number of pairs is an error and not a silent return to the default.
#' @keywords internal
#' @noRd
.exact_pairs_max <- function() {
  lim <- getOption("spatialkit.exact_pairs_max", .EXACT_PAIRS_MAX)
  if (!is.numeric(lim) || length(lim) != 1L || is.na(lim) || lim < 0)
    stop("options(spatialkit.exact_pairs_max = ) must be a single number of pairs, ",
         "0 or more; got ",
         if (length(lim) != 1L) sprintf("%s of length %d", class(lim)[1L], length(lim))
         else if (is.numeric(lim)) format(lim) else class(lim)[1L],
         ".", call. = FALSE)
  as.numeric(lim)
}


#' Stop before an exact search that would compare too many pairs
#'
#' An exact search compares every row asked for with every candidate, so its
#' cost is a product the caller knows before it starts.  Above the limit the
#' call stops, with an error of class \code{spatialkit_work_limit} that names
#' the search, the product, the limit and the way on.
#'
#' @param work The number of pairs the search would compare.
#' @param limit The most it may compare: \code{.exact_pairs_max()}.
#' @param what A short name of the search, kept in the condition's
#'   \code{what} field.
#' @param caller Name of the exported function, for the message.
#' @param remedy The message's last sentence.
#' @param search The search in the words of the message; \code{what} when not
#'   given.
#' @param data Further fields for the condition, as a named list.
#' @return \code{NULL}, invisibly, when \code{work} is within \code{limit}.
#' @keywords internal
#' @noRd
.work_guard <- function(work, limit, what, caller,
                        remedy = "Install FNN: the result is the same.",
                        search = what, data = list()) {
  if (!(work > limit)) return(invisible(NULL))
  msg <- sprintf(paste0("%s(): %s would compare %.2g pairs without the FNN package, above ",
                        "the limit of %.2g (options(spatialkit.exact_pairs_max = )). %s"),
                 caller, search, work, limit, remedy)
  stop(do.call(errorCondition,
               c(list(message = msg, class = "spatialkit_work_limit", call = NULL,
                      work = as.numeric(work), limit = as.numeric(limit), what = what),
                 data)))
}


# -----------------------------------------------------------------------------
# Small pieces of the arithmetic
# -----------------------------------------------------------------------------

#' The weight a pooled value gives the region's own mean
#'
#' @param tau2 Mean squared distance between a region's level and its target.
#' @param v Variance of the region's own mean, \code{sigma2 / n}.
#' @return \code{tau2 / (tau2 + v)}; 1 where both are 0.
#' @keywords internal
#' @noRd
.pool_weight <- function(tau2, v) ifelse(tau2 + v > 0, tau2 / (tau2 + v), 1)


#' Whether a sum of squares is rounding error on the scale of the values
#'
#' The test \code{resolution_profile()} applies to a constant response.  A
#' column with no spread still has a sum of squares above zero when its value
#' is not exact in binary (0.1, 1/3), and dividing one such sum by another
#' gives weights that mean nothing.
#'
#' @param ss A sum of squared deviations of \code{y}.
#' @param y The values.
#' @return \code{TRUE} when \code{ss} is at most
#'   \code{100 * .Machine$double.eps^2 * sum(y^2)}.
#' @keywords internal
#' @noRd
.no_spread <- function(ss, y) ss <= 100 * .Machine$double.eps^2 * sum(y^2)


#' The sum of squares inside the regions, about means without the rounding of
#' their sums
#'
#' A region's mean is its sum over its count, and a sum of values that are
#' not exact in binary is off by a rounding error that grows with the count.
#' In a region of 500 rows that all hold 0.1 the mean is not 0.1, and the sum
#' of squares about it is that error and not zero: above 100 rows or so it
#' passes the bound of \code{.no_spread()}, and a column with no spread
#' inside any region would be fitted as if it had some.  Here each mean is
#' corrected by the mean of the residuals about it, which is how
#' \code{mean()} corrects its own sum, and the squares are taken about the
#' corrected means.  The result is never above the sum about the uncorrected
#' means by more than rounding.
#'
#' @param res The rows' residuals about their region's mean, \code{y - m[gi]}.
#' @param gi Region of each row, \code{1..K}.
#' @param n Rows per region.
#' @return One number.
#' @keywords internal
#' @noRd
.ss_within <- function(res, gi, n) sum((res - (as.numeric(rowsum(res, gi)) / n)[gi])^2)


#' A coordinate matrix as the searches need it
#'
#' The kd-tree is handed every column of the matrix it is given, and the
#' exact search and the tie rule read the first two.  A third column (a
#' height, a part number) would therefore give one answer with \pkg{FNN} and
#' another without, so it is refused and not dropped.  Whole-number storage
#' is made double: \code{rowsum()} of integers returns \code{NA} without a
#' warning once a sum passes \code{.Machine$integer.max}.
#'
#' @param xy The matrix.
#' @param rows The number of rows it must have, or \code{NULL} for any.
#' @param who Name of the calling function, for the message.
#' @param arg Name of the argument, for the message.
#' @return \code{xy} as a matrix of doubles.
#' @keywords internal
#' @noRd
.pool_xy <- function(xy, rows, who, arg = "xy") {
  if (!is.matrix(xy) || !is.numeric(xy) || ncol(xy) != 2L ||
      (!is.null(rows) && nrow(xy) != rows))
    stop(sprintf("%s(): `%s` must be a numeric matrix with two columns%s.", who, arg,
                 if (is.null(rows)) "" else " and one row per value"), call. = FALSE)
  if (!is.double(xy)) storage.mode(xy) <- "double"
  xy
}


#' The number of regions a size gives
#'
#' @param n Number of rows the size divides.
#' @param obs_per_cell Rows per region.
#' @return \code{n / obs_per_cell} rounded as R rounds (a half goes to the
#'   even number), and at least 2.
#' @keywords internal
#' @noRd
.count_at <- function(n, obs_per_cell) pmax(2L, as.integer(round(n / obs_per_cell)))


#' The two offsets of the interval for one observation
#'
#' Quantiles of the leave-one-out residuals at \code{(1 - level) / 2} and its
#' complement.  The probability is rounded to 12 places: \code{(1 - 0.9) / 2}
#' is \code{0.04999999999999999} in floating point, and a quantile taken
#' there differs in its last bits from the one at 0.05.
#'
#' @param e Leave-one-out residuals.
#' @param level The level of the interval, strictly between 0 and 1.
#' @return \code{c(lo, hi)}, to be added to a region's value; two \code{NA}
#'   when there are fewer than two residuals or one is missing.
#' @keywords internal
#' @noRd
.obs_offsets <- function(e, level) {
  if (length(e) < 2L || anyNA(e)) return(c(NA_real_, NA_real_))
  p <- round((1 - level) / 2, 12)
  stats::quantile(e, c(p, 1 - p), names = FALSE, type = 7)
}


#' Whether regions predict a left-out row better than one mean for all rows
#'
#' Two things have to hold.  The leave-one-out error of one mean is above
#' the regions' by more than \code{tol}; and the region means differ from one
#' mean by more than chance, by an F test at 5 percent after a Bonferroni
#' adjustment for the levels the count was searched over.  The first alone
#' has no sampling error in it.  This is the rule \code{select_resolution()}
#' applies to \code{cp}, with the leave-one-out error in its place.
#'
#' @param loo Leave-one-out error of the regions' values.
#' @param loo_one The same error with one mean for all rows.
#' @param p The p-value of the F test of the region means against one mean.
#' @param tol The margin, as a fraction of \code{loo}.
#' @param n_searched Number of levels the count was searched over.
#' @return A list: \code{structured} (\code{TRUE}, \code{FALSE}, or \code{NA}
#'   when an error is not finite, or when the regions beat one mean and there
#'   is no p-value), \code{structure_p} (the adjusted p-value) and
#'   \code{beats} (the first condition alone).
#' @keywords internal
#' @noRd
.no_structure <- function(loo, loo_one, p, tol, n_searched = 1L) {
  if (!length(p)) p <- NA_real_
  p_adj <- min(1, p * n_searched)
  beats <- !(loo_one <= loo * (1 + tol))
  structured <-
    if (!is.finite(loo) || !is.finite(loo_one)) NA
    else if (!beats || (is.finite(p_adj) && p_adj >= 0.05)) FALSE
    else if (is.finite(p_adj)) TRUE
    else NA
  list(structured = structured, structure_p = p_adj, beats = beats)
}


#' Column names in the words of a warning
#'
#' @param cols Character vector of names.
#' @return One string: \code{'a'}; \code{'a' and 'b'}; \code{'a', 'b' and
#'   'c'}; beyond six names the first six and the number of the rest,
#'   \code{'a', 'b', 'c', 'd', 'e', 'f' and 4 more}.
#' @keywords internal
#' @noRd
.cols_text <- function(cols) {
  q <- sprintf("'%s'", cols)
  k <- length(q)
  if (k == 0L) return("")
  if (k == 1L) return(q)
  if (k <= 6L) return(paste(paste(q[-k], collapse = ", "), "and", q[k]))
  sprintf("%s and %d more", paste(q[1:6], collapse = ", "), k - 6L)
}


# -----------------------------------------------------------------------------
# Neighbours of the regions
# -----------------------------------------------------------------------------

#' The regions whose centroids are nearest each region's
#'
#' The centroids are put in coordinate order (first coordinate, then second),
#' searched, and the result mapped back to the regions' numbers.  A tie
#' between two candidates therefore goes to the lower first coordinate, then
#' the lower second, whatever the regions are numbered.  Regions whose
#' centroids coincide exactly are told apart by their number, the one place
#' the numbering decides anything.  A region is never its own neighbour.
#'
#' A centroid that shares its point with \code{k} others or more is not
#' searched for at all: its neighbours are the first \code{k} of them in that
#' order, all at distance 0.  That is what a search returns for it, and it
#' keeps a layer whose regions lie on one point by the thousand (rows that
#' carry the coordinates of a coarser unit, or a placeholder) from costing
#' the square of their number: a kd-tree asked about such a centroid has to
#' list every other one on its point before it can tell which come first.
#'
#' For the rest, with \pkg{FNN} the kd-tree returns one neighbour more than
#' is wanted, where there are that many other regions (ten regions or more
#' for eight neighbours), and R recomputes the squared distances.  A row is
#' settled when it is not among its own neighbours and its distances rise by
#' more than one part in a thousand million from each neighbour to the next;
#' it keeps the tree's first \code{k}.  Every other row is searched again by
#' \code{.knn_tied_rows()}, by distance and then by row, which is the rule of
#' the exact search.  Without \pkg{FNN} each of them is searched exactly,
#' which compares it with every centroid; the pairs that makes, the rows
#' searched for times the centroids, are what \code{.work_guard()} stops
#' above the limit.  The two paths return the same matrix.
#'
#' @param cen Two-column matrix of centroids, one row per region.
#' @param k Number of neighbours wanted.
#' @param use_fnn Whether to use \pkg{FNN}.
#' @param caller Name of the exported function, for the work guard's message.
#' @return An integer matrix with a row per region and \code{min(k, K - 1)}
#'   columns, nearest first.
#' @keywords internal
#' @noRd
.region_neighbours <- function(cen, k, use_fnn, caller = "resolution_profile") {
  cen <- .pool_xy(cen, NULL, ".region_neighbours", "cen")
  K <- nrow(cen)
  k <- min(as.integer(k), K - 1L)
  if (k < 1L) return(matrix(0L, K, 0L))
  o <- order(cen[, 1L], cen[, 2L])
  u <- unname(cen[o, , drop = FALSE])
  x <- u[, 1L]; y <- u[, 2L]
  # Centroids on one point stand together in this order, lowest number first.
  first <- c(TRUE, x[-1L] != x[-K] | y[-1L] != y[-K])
  run   <- cumsum(first)
  from  <- which(first)[run]                # where each centroid's point starts
  many  <- tabulate(run)[run] > k           # k others or more on its own point
  idx   <- matrix(0L, K, k)
  if (any(many)) {
    # The first k + 1 of the point, less the centroid itself when it is one.
    at <- which(many)
    j  <- matrix(rep(seq_len(k) - 1L, each = length(at)), length(at), k)
    idx[at, ] <- from[at] + j + (j >= at - from[at])
  }
  rest <- which(!many)
  if (!isTRUE(use_fnn)) {
    .work_guard(as.numeric(length(rest)) * K, .exact_pairs_max(), what = "cell neighbours",
                caller = caller,
                search = sprintf("finding each cell's neighbours among %d cells", K))
    idx[rest, ] <- .knn_exact_rows(x, y, rest, k)
  } else if (length(rest)) {
    kq  <- min(k + 1L, K - 1L)
    nn  <- matrix(FNN::get.knn(u, k = kq)$nn.index, K, kq)[rest, , drop = FALSE]
    d2  <- matrix((x[nn] - x[rest])^2 + (y[nn] - y[rest])^2, length(rest), kq)
    own <- rowSums(nn == rest) > 0L
    tied <- if (kq > 1L)
      rowSums(d2[, -1L, drop = FALSE] <= d2[, -kq, drop = FALSE] * (1 + 1e-9)) > 0L
    else rep(FALSE, length(rest))
    idx[rest, ] <- nn[, seq_len(k), drop = FALSE]
    redo <- rest[own | tied]
    if (length(redo)) idx[redo, ] <- .knn_tied_rows(u, redo, k)
  }
  out <- matrix(0L, K, k)
  out[o, ] <- o[idx]
  out
}


#' Centroids of the regions and each region's neighbours
#'
#' The centroid of a region is the mean of the coordinates of its rows.  The
#' rows are summed in coordinate order, so the centroids are the same numbers
#' to the last bit whatever order the rows come in.  On gridded locations
#' that is what keeps the neighbours fixed: summed in row order the centroids
#' of a lattice differ in their last bits from one row order to the next, and
#' exact ties between them fall differently.
#'
#' @param gi Region of each row, \code{1..K}.
#' @param K Number of regions.
#' @param xy Two-column matrix of the rows' coordinates; any other number of
#'   columns is an error.
#' @param n Rows per region.
#' @param use_fnn Whether to use \pkg{FNN}.
#' @param caller Name of the exported function, for the work guard's message.
#' @return A list: \code{cen} (a \code{K x 2} matrix) and \code{nb} (the
#'   matrix of \code{.region_neighbours()} for \code{.POOL_NEIGHBOURS}
#'   neighbours).  It depends on the rows and their regions and not on the
#'   values, so one serves every column that has the same rows.
#' @keywords internal
#' @noRd
.pool_geometry <- function(gi, K, xy, n, use_fnn, caller = "resolution_profile") {
  xy  <- .pool_xy(xy, length(gi), ".pool_geometry")
  o   <- order(xy[, 1L], xy[, 2L])
  cen <- rowsum(xy[o, , drop = FALSE], gi[o]) / n
  list(cen = cen, nb = .region_neighbours(cen, .POOL_NEIGHBOURS, use_fnn, caller))
}


# -----------------------------------------------------------------------------
# The nearest other row
# -----------------------------------------------------------------------------

#' For chosen locations, the nearest other location, searched
#'
#' With \pkg{FNN} the kd-tree returns the three nearest locations of each one
#' asked for, itself among them, and R recomputes the squared distances.  A
#' location is settled when it is one of its own three and its two nearest
#' others differ by more than one part in a thousand million.  Every other
#' location is searched again by \code{.knn_tied_rows()}, which settles a tie
#' by row, so a lattice costs no comparison of every pair.  Without \pkg{FNN}
#' the search is exact and is stopped by \code{.work_guard()} above the
#' limit.  The two paths return the same locations.
#'
#' @param u Two-column matrix of distinct locations in coordinate order, at
#'   least two.
#' @param q The locations to search for, as rows of \code{u}.
#' @param use_fnn Whether to use \pkg{FNN}.
#' @param caller Name of the exported function, for the work guard's message.
#' @return For each of \code{q} the row of \code{u} nearest it and not it; a
#'   tie goes to the lower row, which is the lower first coordinate and then
#'   the lower second.
#' @keywords internal
#' @noRd
.nearest_location_search <- function(u, q, use_fnn, caller = "resolution_profile") {
  U  <- nrow(u)
  nq <- length(q)
  if (!isTRUE(use_fnn)) {
    .work_guard(as.numeric(nq) * U, .exact_pairs_max(), what = "nearest other location",
                caller = caller,
                search = sprintf("finding the nearest other location of %d of %d locations",
                                 nq, U))
    return(.knn_exact_rows(u[, 1L], u[, 2L], q, 1L)[, 1L])
  }
  kq <- min(3L, U)
  nn <- matrix(FNN::get.knnx(u, u[q, , drop = FALSE], k = kq)$nn.index, nq, kq)
  d2 <- matrix((u[nn, 1L] - u[q, 1L])^2 + (u[nn, 2L] - u[q, 2L])^2, nq, kq)
  own <- nn == q
  d2[own] <- Inf
  at   <- cbind(seq_len(nq), max.col(-d2, ties.method = "first"))
  best <- d2[at]
  out  <- nn[at]
  d2[at] <- Inf
  # The nearest of what is left: the second nearest other location, when the
  # three hold one.
  second <- if (kq > 2L) pmin(d2[, 1L], d2[, 2L], d2[, 3L]) else rep(Inf, nq)
  redo <- which(rowSums(own) != 1L | !(second > best * (1 + 1e-9)))
  if (length(redo)) out[redo] <- .knn_tied_rows(u, q[redo], 1L)[, 1L]
  out
}


#' For chosen locations, the nearest other location, searched once
#'
#' The nearest other location of a location does not depend on the regions.
#' With \code{cache} given, each one found is kept there and no location is
#' searched twice, so fits of the same rows in one set of regions after
#' another share their searches.  The cache holds the locations it was filled
#' for and starts again when it is handed others.
#'
#' @param u,q,use_fnn,caller As for \code{.nearest_location_search()}.
#' @param cache \code{NULL}, or an environment.
#' @return As \code{.nearest_location_search()}.
#' @keywords internal
#' @noRd
.nearest_location <- function(u, q, use_fnn, cache = NULL, caller = "resolution_profile") {
  if (is.null(cache)) return(.nearest_location_search(u, q, use_fnn, caller))
  if (!identical(cache$u, u)) {
    cache$u <- u
    cache$nearest <- rep(NA_integer_, nrow(u))
  }
  todo <- q[is.na(cache$nearest[q])]
  if (length(todo))
    cache$nearest[todo] <- .nearest_location_search(u, todo, use_fnn, caller)
  cache$nearest[q]
}


#' For chosen rows, a nearest other row
#'
#' Chosen so that the region of the row returned does not depend on the order
#' of the rows.
#' \enumerate{
#'   \item The rows are put in the order of their first coordinate, their
#'     second, then their region.  A run of equal coordinates is a location,
#'     and the first row of a location is the one in the lowest-numbered
#'     region.
#'   \item A row that shares its location with other rows takes the first
#'     other row there: the first row of the location, or the second when it
#'     is the first itself.
#'   \item A row alone at its location takes the first row of its nearest
#'     other location; a tie between locations goes to the lower first
#'     coordinate, then the lower second.
#' }
#' An exact search on the rows as given breaks a tie by row number, so on
#' gridded locations the neighbour would change with the order of the rows.
#' Where no two locations tie and none repeats, the row returned is the
#' nearest other row.
#'
#' @param xy Two-column matrix of finite coordinates, one row per row; any
#'   other number of columns is an error.
#' @param gi Region of each row.
#' @param rows The rows to find a neighbour for.
#' @param use_fnn Whether to use \pkg{FNN}.  The answer is the same without.
#' @param cache \code{NULL}, or an environment in which the nearest other
#'   location of each location searched is kept (\code{.nearest_location()}).
#' @param caller Name of the exported function, for the work guard's message.
#' @return For each of \code{rows}, a row of \code{xy} that is not it;
#'   \code{NA} when \code{xy} has one row.
#' @keywords internal
#' @noRd
.nearest_other <- function(xy, gi, rows, use_fnn, cache = NULL,
                           caller = "resolution_profile") {
  xy <- .pool_xy(xy, length(gi), ".nearest_other")
  if (!length(rows)) return(integer(0))
  n  <- nrow(xy)
  o  <- order(xy[, 1L], xy[, 2L], gi)
  xs <- xy[o, , drop = FALSE]
  first <- c(TRUE, xs[-1L, 1L] != xs[-n, 1L] | xs[-1L, 2L] != xs[-n, 2L])
  loc   <- cumsum(first)                 # the location of each sorted row
  head  <- which(first)                  # the first sorted row of each location
  m     <- tabulate(loc)                 # rows at each location
  u     <- unname(xs[first, , drop = FALSE])
  at    <- integer(n)
  at[o] <- seq_len(n)
  pos   <- at[rows]                      # where each row asked for stands
  L     <- loc[pos]
  ans   <- rep(NA_integer_, length(rows))
  shared <- m[L] > 1L
  if (any(shared)) {
    h <- head[L[shared]]
    ans[shared] <- ifelse(pos[shared] == h, h + 1L, h)
  }
  if (!all(shared) && nrow(u) > 1L) {
    q  <- sort(unique(L[!shared]))
    nl <- .nearest_location(u, q, use_fnn, cache, caller)
    ans[!shared] <- head[nl[match(L[!shared], q)]]
  }
  o[ans]
}


# -----------------------------------------------------------------------------
# The fit
# -----------------------------------------------------------------------------

#' Refuse centroids and neighbours that cannot be those of \code{K} regions
#'
#' A fit reads the \code{geo} it is handed in place of computing its own.
#' One of another number of regions would not stop it: \code{S[nb]} is
#' reshaped to \code{K} rows whatever \code{nb} holds, and the targets come
#' out as sums over other regions, with no \code{NA} and no error.  This
#' checks what can be checked without computing the geometry again: the two
#' shapes, that every neighbour is a region, and that none is the region
#' itself.  It cannot tell the geometry of these regions from that of other
#' regions of the same number; the caller answers for that.
#'
#' @param geo The list to check: \code{cen} and \code{nb}.
#' @param K Number of regions.
#' @return \code{NULL}, invisibly; an error when \code{geo} does not fit.
#' @keywords internal
#' @noRd
.pool_geo_check <- function(geo, K) {
  k  <- min(.POOL_NEIGHBOURS, K - 1L)
  ok <- is.list(geo) &&
    is.matrix(geo$cen) && is.numeric(geo$cen) && identical(dim(geo$cen), c(K, 2L)) &&
    !anyNA(geo$cen) &&
    is.matrix(geo$nb) && is.numeric(geo$nb) && identical(dim(geo$nb), c(K, k)) &&
    !anyNA(geo$nb) && all(geo$nb >= 1L & geo$nb <= K) && !any(geo$nb == seq_len(K))
  if (!ok)
    stop(sprintf(paste0(".pool_fit(): `geo` must hold the centroids and neighbours of these ",
                        "%d regions: `cen` %d x 2, and `nb` %d x %d with every entry another ",
                        "region."), K, K, K, k), call. = FALSE)
  invisible(NULL)
}


#' Pooled values, plain means and their leave-one-out errors for one column
#'
#' The first of these that holds decides what is returned, and \code{status}
#' is its name.
#' \describe{
#'   \item{\code{"one cell"}}{\code{K} is 1.  The value is the plain mean and
#'     its weight 1; there is no target and no \code{tau2}.  Each residual is
#'     the plain one, \code{(y - m) N / (N - 1)}, and every leave-one-out
#'     error their mean square.  For a constant column \code{sigma2} and the
#'     residuals are exactly 0 and the errors \code{NA}.  With one row,
#'     \code{sigma2}, the residuals and the errors are \code{NA}.}
#'   \item{\code{"no cell holds two values"}}{\code{N <= K}: \code{sigma2}
#'     has no estimate and nothing is fitted.  \code{n}, \code{S}, \code{m},
#'     \code{N}, \code{K}, \code{constant} and \code{flat} are returned and
#'     everything else is \code{NA} (\code{NULL} for the residuals).}
#'   \item{\code{"no spread inside any cell"}}{\code{flat}: the sum of
#'     squares inside the regions is rounding error, by \code{.no_spread()}
#'     on the sum the fit uses or on the same sum about means corrected for
#'     the rounding of their own sums (\code{.ss_within()}).  \code{sigma2}
#'     is 0, every weight 1 and the value the plain mean; the residual of a
#'     row that shares its region is 0.  When the column is constant over all
#'     its rows as well (\code{constant}), \code{tau2} is 0, every residual
#'     is 0 and every leave-one-out error is \code{NA}: no level fits a
#'     constant better than another.  These are set, not computed, because
#'     the formulas would divide one rounding error by another.}
#'   \item{\code{"no signal between cells"}}{\code{tau2} is 0: every weight
#'     is 0 and each value is its target.}
#'   \item{\code{"ok"}}{None of the above.}
#' }
#'
#' Leaving one row out: under plain means a row in a region of two or more is
#' predicted by the mean of the others, and a row alone in its region by the
#' mean of the region of its nearest other row (\code{.nearest_other()}).
#' Under pooled values the region's mean and own weight are computed again
#' without the row and \code{sigma2}, \code{tau2} and the target are held; a
#' region left empty takes its target.  The target never holds the row, since
#' a region is not its own neighbour.  Rows at one location are left out one
#' at a time.
#'
#' @param y Finite values, at least one.  They are read as doubles, so a
#'   column of whole numbers gives the fit of \code{as.numeric()} of it.
#' @param gi Region of each, integers \code{1..K}, every region present.
#' @param K Number of regions.
#' @param xy Two-column matrix of finite coordinates in a projected CRS, one
#'   row per value.  A matrix with a third column is an error: the caller
#'   takes the two it means.
#' @param want Which estimates' leave-one-out residuals to compute: a
#'   non-empty subset of \code{c("pooled", "mean")}.  Everything else is
#'   computed whatever it holds.
#' @param nn1 For each row its nearest other row, or \code{NULL} to have it
#'   found for the rows that need one: those alone in their region, when
#'   \code{"mean"} is wanted.
#' @param geo \code{NULL}, or \code{fit[c("cen", "nb")]} of an earlier fit on
#'   the same rows and regions; neither is then computed again.  One whose
#'   shape is not that of \code{K} regions is an error
#'   (\code{.pool_geo_check()}); one of other regions of the same number
#'   cannot be told from the right one and gives wrong targets.
#' @param use_fnn Whether to use \pkg{FNN}.  The answer is the same without.
#' @param cache \code{NULL}, or an environment for \code{.nearest_other()}.
#' @param caller Name of the exported function, for the message of the error
#'   that stops an exact search above the limit
#'   (\code{spatialkit_work_limit}).
#' @return A list: \code{n}, \code{S}, \code{m} (rows, sum and plain mean of
#'   each region), \code{N}, \code{K}, \code{sigma2}, \code{cen}, \code{nb},
#'   \code{target}, \code{Q}, \code{tau2}, \code{w}, \code{value} (the pooled
#'   value), \code{resid_pooled}, \code{resid_mean} (leave-one-out residuals,
#'   one per row; \code{NULL} for an estimate not in \code{want}),
#'   \code{loo_pooled}, \code{loo_pooled_se}, \code{loo_mean},
#'   \code{loo_mean_se} (their mean square and its standard error over rows;
#'   \code{NA} for an estimate not in \code{want}), \code{loo_one},
#'   \code{loo_one_se} (the same for one mean for all rows), \code{f_p} (the
#'   p-value of the F test of the region means against one mean; \code{NA}
#'   when \code{flat}, with one region, or when \code{N - K < 1}),
#'   \code{constant}, \code{flat} and \code{status}.
#' @keywords internal
#' @noRd
.pool_fit <- function(y, gi, K, xy, want = c("pooled", "mean"), nn1 = NULL, geo = NULL,
                      use_fnn = .has_backend("FNN"), cache = NULL,
                      caller = "resolution_profile") {
  want <- unique(want)
  if (!length(want) || !all(want %in% c("pooled", "mean")))
    stop(".pool_fit(): `want` must be \"pooled\", \"mean\" or both.", call. = FALSE)
  # Whole numbers are summed as doubles: rowsum() of integers returns NA
  # without a warning once a region's sum passes .Machine$integer.max.
  y <- as.numeric(y)
  K <- as.integer(K)
  N <- length(y)
  if (N < 1L || length(gi) != N)
    stop(".pool_fit(): `y` must hold at least one value and `gi` the region of each.",
         call. = FALSE)
  # tabulate() leaves out what is not among 1..K, so the first comparison
  # fails for a number outside them and the second for one with no row.
  n <- tabulate(gi, K)
  if (sum(n) != N || !all(n > 0L))
    stop(".pool_fit(): the regions must be numbered 1 to K with a row in each.",
         call. = FALSE)
  S <- as.numeric(rowsum(y, gi))
  xy <- .pool_xy(xy, N, ".pool_fit")
  if (!is.null(geo)) .pool_geo_check(geo, K)
  m <- S / n
  res   <- y - m[gi]
  ss_w  <- sum(res^2)
  tss   <- sum((y - mean(y))^2)
  const <- .no_spread(tss, y)                       # constant over the rows
  # No spread inside any region.  The second test is on the sum of squares the
  # fit uses; the third is on the same sum about means that carry no rounding
  # of their own, without which a column constant inside large regions fails
  # the second (.ss_within()).
  flat  <- const || .no_spread(ss_w, y) || .no_spread(.ss_within(res, gi, n), y)
  # Mean square of a set of residuals, and its standard error over rows.
  loo_of <- function(e) c(mean(e^2), stats::sd(e^2) / sqrt(N))
  none   <- c(NA_real_, NA_real_)

  out <- list(n = n, S = S, m = m, N = N, K = K, sigma2 = NA_real_,
              cen = NA_real_, nb = NA_integer_, target = rep(NA_real_, K),
              Q = NA_real_, tau2 = NA_real_, w = rep(NA_real_, K),
              value = rep(NA_real_, K), resid_pooled = NULL, resid_mean = NULL,
              loo_pooled = NA_real_, loo_pooled_se = NA_real_,
              loo_mean = NA_real_, loo_mean_se = NA_real_,
              loo_one = NA_real_, loo_one_se = NA_real_, f_p = NA_real_,
              constant = const, flat = flat, status = NA_character_)

  # ---- 1. One region: nothing to pool toward -------------------------------
  if (K == 1L) {
    if (is.null(geo)) geo <- .pool_geometry(gi, K, xy, n, use_fnn, caller)
    e   <- if (N == 1L) NA_real_ else if (const) numeric(N) else res * N / (N - 1)
    loo <- if (N == 1L || const) none else loo_of(e)
    out$sigma2 <- if (N == 1L) NA_real_ else if (const) 0 else ss_w / (N - 1)
    out$cen <- geo$cen
    out$nb  <- geo$nb
    out$w   <- 1
    out$value <- m
    if ("pooled" %in% want) {
      out$resid_pooled <- e; out$loo_pooled <- loo[1L]; out$loo_pooled_se <- loo[2L]
    }
    if ("mean" %in% want) {
      out$resid_mean <- e; out$loo_mean <- loo[1L]; out$loo_mean_se <- loo[2L]
    }
    out$loo_one <- loo[1L]; out$loo_one_se <- loo[2L]
    out$flat    <- const
    out$status  <- "one cell"
    return(out)
  }

  # ---- 2. No region holds two rows: sigma2 has no estimate -----------------
  if (N <= K) {
    out$status <- "no cell holds two values"
    return(out)
  }

  # ---- 3 to 5. The fit ------------------------------------------------------
  sigma2 <- if (flat) 0 else ss_w / (N - K)
  if (is.null(geo)) geo <- .pool_geometry(gi, K, xy, n, use_fnn, caller)
  cen    <- geo$cen                                 # centroids, summed in coordinate order
  nb     <- geo$nb                                  # K x min(8, K - 1)
  Sn     <- rowSums(matrix(S[nb], K)); Nn <- rowSums(matrix(n[nb], K))
  target <- Sn / Nn
  Q      <- sum(n * (m - target)^2)
  tau2   <- if (const) 0 else max(0, (Q - K * sigma2) / N)
  w      <- if (flat) rep(1, K) else .pool_weight(tau2, sigma2 / n)
  value  <- if (flat) m else w * m + (1 - w) * target

  # Leave one row out.
  nj <- n[gi]; n1 <- nj - 1L; d1 <- pmax(n1, 1L)
  if ("mean" %in% want) {
    eP  <- if (flat) numeric(N) else res * nj / d1
    sng <- nj == 1L
    if (!const && any(sng)) {
      if (is.null(nn1)) {
        nn1 <- integer(N)
        nn1[sng] <- .nearest_other(xy, gi, which(sng), use_fnn, cache, caller)
      }
      eP[sng] <- y[sng] - m[gi[nn1[sng]]]
    }
    loo <- if (const) none else loo_of(eP)
    out$resid_mean <- eP; out$loo_mean <- loo[1L]; out$loo_mean_se <- loo[2L]
  }
  if ("pooled" %in% want) {
    m1 <- (S[gi] - y) / d1
    w1 <- ifelse(n1 >= 1L, if (flat) 1 else .pool_weight(tau2, sigma2 / d1), 0)
    eL <- if (const) numeric(N) else y - (w1 * m1 + (1 - w1) * target[gi])
    if (flat && !const) eL[n1 >= 1L] <- 0
    loo <- if (const) none else loo_of(eL)
    out$resid_pooled <- eL; out$loo_pooled <- loo[1L]; out$loo_pooled_se <- loo[2L]
  }
  # One mean for all rows: the reference every count has to beat.
  e1  <- (y - mean(y)) * N / (N - 1)
  loo <- if (const) none else loo_of(e1)

  out$sigma2 <- sigma2
  out$cen    <- cen
  out$nb     <- nb
  out$target <- target
  out$Q      <- Q
  out$tau2   <- tau2
  out$w      <- w
  out$value  <- value
  out$loo_one <- loo[1L]; out$loo_one_se <- loo[2L]
  out$f_p <- if (flat) NA_real_ else
    stats::pf(max(((tss - ss_w) / (K - 1)) / (ss_w / (N - K)), 0), K - 1, N - K,
              lower.tail = FALSE)
  out$status <- if (flat) "no spread inside any cell"
                else if (tau2 == 0) "no signal between cells"
                else "ok"
  out
}


# -----------------------------------------------------------------------------
# One column of an assigned layer
# -----------------------------------------------------------------------------

#' Region values of one column, from its rows
#'
#' A row is scored when its value is finite, its region's ID is not missing
#' and its coordinates are finite; every other row is, for this column, as if
#' absent.  The regions are the distinct IDs among the scored rows, in the
#' order \code{u[order(u, method = "radix")]}: ascending for numbers, the
#' C-locale order for strings on every platform, the order of the levels for
#' a factor.
#'
#' \code{status} is that of \code{.pool_fit()}, or one of three set here,
#' under each of which the values are \code{NA}: \code{"no coordinates"}
#' (\code{xy} is \code{NULL}) and \code{"no located value"} (no scored row),
#' both before any fit, and \code{"work limit"} when the fit is stopped
#' before an exact search above the limit.  That stop is an error of class
#' \code{spatialkit_work_limit} inside the fit and a status here: the caller
#' still has its plain columns to return.  No warning is raised here.
#'
#' @param y The column's values, one per row.
#' @param id The region of each row, of any type; \code{NA} for a row in
#'   none.
#' @param xy Two-column matrix of the rows' coordinates in a projected CRS,
#'   or \code{NULL} when the rows have none.
#' @param obs_level \code{NULL}, or the level of the interval for one
#'   observation.
#' @param want,use_fnn As for \code{.pool_fit()}.
#' @param geo \code{NULL}, or the \code{geo} this function returned for an
#'   earlier column with the same \code{mask} and the same \code{id}.  The
#'   fit refuses one of another number of regions and cannot tell one of
#'   other regions of the same number, so the caller answers for passing it
#'   only between columns whose scored rows are the same.
#' @return A list.  Per region, in the order of \code{ids}: \code{ids},
#'   \code{n} (scored rows), \code{value}, \code{w}, \code{target}, and
#'   \code{obs_lo}, \code{obs_hi} (\code{NULL} without \code{obs_level}).
#'   \code{record}: the column's entry in the record of a summary, with
#'   \code{column} and \code{role} left \code{NA} for the caller.
#'   \code{mask}: which rows are scored.  \code{geo}: the centroids and
#'   neighbours the fit used, its own or the ones handed in, for a later
#'   column with the same mask; \code{NULL} when no fit used any.
#'   \code{fit}: the list of \code{.pool_fit()}, or \code{NULL} when no fit
#'   ran or it was stopped.
#'   \code{work}, \code{limit}: the pairs and the limit of the search that
#'   stopped the fit, \code{NA} otherwise.
#' @keywords internal
#' @noRd
.pool_column <- function(y, id, xy, obs_level, want, use_fnn, geo = NULL) {
  rows   <- length(y)
  valued <- is.finite(y) & !is.na(id)
  offsets <- if (is.null(obs_level)) NULL else c(NA_real_, NA_real_)
  record <- list(column = NA_character_, role = NA_character_,
                 n = 0L, cells = 0L, df = NA_integer_, neighbours = NA_integer_,
                 sigma2 = NA_real_, tau2 = NA_real_, own_weight = NA_real_,
                 share_alone = NA_real_, loo = NA_real_, loo_se = NA_real_,
                 loo_one_cell = NA_real_, structured = NA, structure_p = NA_real_,
                 offsets = offsets, repeated = 0L, skipped = 0L, constant = NA,
                 status = NA_character_)
  out <- list(ids = id[0L], n = integer(0), value = numeric(0), w = numeric(0),
              target = numeric(0), obs_lo = NULL, obs_hi = NULL, record = record,
              mask = logical(rows), geo = NULL, fit = NULL,
              work = NA_real_, limit = NA_real_)
  if (!is.null(obs_level)) out$obs_lo <- out$obs_hi <- numeric(0)

  if (is.null(xy)) {
    out$record$skipped <- sum(valued)
    out$record$status  <- "no coordinates"
    return(out)
  }
  located <- is.finite(xy[, 1L]) & is.finite(xy[, 2L])
  mask    <- valued & located
  out$mask <- mask
  out$record$skipped <- sum(valued & !located)
  if (!any(mask)) {
    out$record$status <- "no located value"
    return(out)
  }

  ys  <- y[mask]                            # the fit reads it as doubles
  ids <- id[mask]
  xs  <- cbind(xy[mask, 1L], xy[mask, 2L])
  u   <- unique(ids)
  u   <- u[order(u, method = "radix")]
  gi  <- match(ids, u)
  K   <- length(u)
  N   <- length(ys)
  n   <- tabulate(gi, K)
  na_K <- rep(NA_real_, K)
  out$ids <- u
  out$n   <- n
  out$value <- out$w <- out$target <- na_K
  if (!is.null(obs_level)) out$obs_lo <- out$obs_hi <- na_K
  out$record$n           <- N
  out$record$cells       <- K
  out$record$df          <- N - K
  out$record$neighbours  <- min(.POOL_NEIGHBOURS, K - 1L)
  out$record$share_alone <- sum(n == 1L) / N
  out$record$repeated    <- .repeated_rows(xs)

  fit <- tryCatch(
    .pool_fit(ys, gi, K, xs, want = want, geo = geo, use_fnn = use_fnn,
              caller = "summarize_by_cell"),
    spatialkit_work_limit = function(e) e)
  if (inherits(fit, "spatialkit_work_limit")) {
    out$work  <- fit$work
    out$limit <- fit$limit
    out$record$status <- "work limit"
    return(out)
  }

  out$fit    <- fit
  out$value  <- fit$value
  out$w      <- fit$w
  out$target <- fit$target
  if (is.matrix(fit$cen)) out$geo <- fit[c("cen", "nb")]
  out$record$sigma2       <- fit$sigma2
  out$record$tau2         <- fit$tau2
  out$record$own_weight   <- stats::median(fit$w)
  out$record$loo          <- fit$loo_pooled
  out$record$loo_se       <- fit$loo_pooled_se
  out$record$loo_one_cell <- fit$loo_one
  out$record$constant     <- fit$constant
  out$record$status       <- fit$status
  if (K >= 2L && !fit$flat && is.finite(fit$loo_pooled)) {
    s <- .no_structure(fit$loo_pooled, fit$loo_one, fit$f_p, tol = .STRUCTURE_TOL)
    out$record$structured  <- s$structured
    out$record$structure_p <- s$structure_p
  }
  if (!is.null(obs_level)) {
    off <- .obs_offsets(fit$resid_pooled, obs_level)
    out$record$offsets <- off
    out$obs_lo <- fit$value + off[1L]
    out$obs_hi <- fit$value + off[2L]
  }
  out
}


#' How many rows share their location with another row
#'
#' @param xy Two-column matrix of finite coordinates.
#' @return The number of rows whose coordinates are those of another row.
#' @keywords internal
#' @noRd
.repeated_rows <- function(xy) {
  n <- nrow(xy)
  if (n < 2L) return(0L)
  o <- order(xy[, 1L], xy[, 2L])
  x <- xy[o, 1L]; y <- xy[o, 2L]
  as_before <- c(FALSE, x[-1L] == x[-n] & y[-1L] == y[-n])
  sum(as_before | c(as_before[-1L], FALSE))
}


# -----------------------------------------------------------------------------
# The intervals a summary's frame holds
# -----------------------------------------------------------------------------

#' The intervals in a summary's frame, and what each covers
#'
#' One row per family of interval columns in the frame.  The interval for one
#' observation is centred on the pooled value; the t interval is centred on
#' the plain mean, and what it covers follows the design effect: the cell's
#' own mean at a design effect of 1, the cell mean as an estimate of the
#' population mean under any other.  Neither row is written when its level
#' was not asked for.
#'
#' @param obs_level Level of the interval for one observation, or
#'   \code{NULL}.
#' @param conf_level Level of the plain mean's t interval, or \code{NULL}.
#' @param deff_label The design effect the t interval was built with, as it
#'   is written in a call: \code{1}, \code{"kish"} (with its quotes),
#'   \code{"variogram"} or a number.  It reads \code{1} exactly when the
#'   interval carries none, which is what \code{.deff_label()} writes.
#' @return A data frame with the columns \code{columns} (the stems of the two
#'   interval columns), \code{covers}, \code{centre}, \code{level} and
#'   \code{method}; no row when both levels are \code{NULL}.
#' @keywords internal
#' @noRd
.intervals_table <- function(obs_level, conf_level, deff_label) {
  obs <- !is.null(obs_level)
  ci  <- !is.null(conf_level)
  own <- identical(as.character(deff_label), "1")
  data.frame(
    columns = c(character(0),
                if (obs) "..obs_lo_, ..obs_hi_",
                if (ci)  "..ci_lo_, ..ci_hi_"),
    covers  = c(character(0),
                if (obs) paste0("one further observation in the cell, located like the data; ",
                                "on average over observations, not cell by cell"),
                if (ci && own) paste0("the cell's own mean, when its observations are spread ",
                                      "through the cell; not a range for one observation"),
                if (ci && !own) paste0("the cell mean as an estimate of the population mean; ",
                                       "not a range for one observation")),
    centre  = c(character(0), if (obs) "pooled value", if (ci) "plain mean"),
    level   = c(numeric(0), if (obs) as.numeric(obs_level), if (ci) as.numeric(conf_level)),
    method  = c(character(0),
                if (obs) "quantiles of leave-one-out residuals",
                if (ci)  sprintf("t, deff = %s", deff_label)),
    stringsAsFactors = FALSE)
}


# -----------------------------------------------------------------------------
# Moran's z of the cell means above the dense limit
# -----------------------------------------------------------------------------

#' Whether Moran's z is left out, on every ladder, where it cannot be computed
#'
#' Above \code{.KNN_DENSE_MAX} cells the neighbour weights of the cell means
#' need \pkg{FNN} and \pkg{Matrix}.  When either is missing a profile
#' leaves \code{moran_i} and \code{moran_z} \code{NA} at such a level, with
#' a warning, whatever its ladder.  Were this \code{FALSE}, only a ladder
#' set by \code{obs_per_cell} would do so, and every other ladder would stop
#' there with the error of \code{.build_knn_weights()}.
#'
#' @return \code{TRUE}.
#' @keywords internal
#' @noRd
.moran_skip_any_ladder <- function() TRUE
