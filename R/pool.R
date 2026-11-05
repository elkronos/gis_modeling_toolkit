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
# Three things the result must not depend on, and how each is kept out:
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
#   the unit and origin     The searches call two squared distances within
#                           .TIE_TOL of each other a tie, so on gridded
#                           locations the last bits of the centroids, which
#                           move with the unit and the origin of the
#                           coordinates, do not decide between neighbours.
#                           Of centroids still tied, the region of the lower
#                           number is taken, so regions on one point stay
#                           told apart by number when rounding parts them.
#                           Centroids within the reach of a tie of one
#                           another count as one point, so a point with more
#                           regions than a region has neighbours stays one.
# =============================================================================

# Regions in a target.
.POOL_NEIGHBOURS <- 8L

# Relative difference up to which two squared distances are a tie in the
# searches for neighbours, here and in the density partition (its neighbour
# and nearest-centre searches, R/partition.R), and in the nearest-centre
# searches of the scoring, the assignment and cell_values_at().  Centroids of
# a lattice that tie exactly in one unit differ in their last bits in
# another, or after a shift of the origin, and each unit in the last place a
# centroid is off moves its squared distances by about 5e-16 times the
# coordinates over the spacing.  A centroid summed from many rows can be off
# by several.  So a tie holds while the spacing is above about a
# ten-millionth of the coordinates, which leaves room for that.  The same
# fraction of the largest coordinate is the reach within which the pooling
# counts centroids as one point (.region_neighbours()).
.TIE_TOL <- 1e-7

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


#' Counts, or a count and a limit on it, written in full
#'
#' Thousands separators and no scientific notation, so that two counts that
#' differ print differently: 1.0004e10 is "10,004,000,000", where \code{"\%.2g"}
#' wrote both it and 1e10 as "1e+10".  Each element is written on its own,
#' with 7 significant digits, or with more, up to 15, where 7 would write two
#' different numbers of \code{x} alike: a limit of 1e10 - 0.4 is
#' "9,999,999,999.6" beside a work of 1e10.  So the numbers of one message
#' are passed together.
#'
#' @param x A numeric vector.
#' @return A character vector of the length of \code{x}.
#' @keywords internal
#' @noRd
.count_text <- function(x) {
  x <- as.numeric(x)
  write <- function(d)
    vapply(x, function(v) format(v, big.mark = ",", scientific = FALSE, trim = TRUE, digits = d),
           character(1))
  distinct <- !duplicated(x)
  for (d in 7:15) {
    out <- write(d)
    if (!anyDuplicated(out[distinct])) break
  }
  out
}


# What the texts that fire when FNN cannot be loaded say to do.  FNN is in
# Imports and installs with spatialkit, so asking the user to install it
# misstates the case: it is there and could not be loaded.
.FNN_REMEDY <- paste0("FNN, which installs with spatialkit, could not be loaded: reinstall ",
                      "it; the result is the same.")


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
#' @param remedy The message's last sentence: \code{.FNN_REMEDY} by default.
#' @param search The search in the words of the message; \code{what} when not
#'   given.
#' @param data Further fields for the condition, as a named list.
#' @return \code{NULL}, invisibly, when \code{work} is within \code{limit}.
#' @keywords internal
#' @noRd
.work_guard <- function(work, limit, what, caller,
                        remedy = .FNN_REMEDY,
                        search = what, data = list()) {
  if (!(work > limit)) return(invisible(NULL))
  nums <- .count_text(c(work, limit))
  msg <- sprintf(paste0("%s(): %s would compare %s pairs without the FNN package, above ",
                        "the limit of %s (options(spatialkit.exact_pairs_max = )). %s"),
                 caller, search, nums[1L], nums[2L], remedy)
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
#' A column with no spread still has a sum of squares above zero when its
#' value is not exact in binary (0.1, 1/3), and dividing one such sum by
#' another gives weights that mean nothing.  The scale is that of \code{y}.
#' On raw values far from 0 it calls a column constant that is not (1e15
#' plus a spread of 1, whose values are stored to an eighth), and on values
#' whose squares overflow or underflow (1e200 or 1e-170 times a spread of 1)
#' it holds whatever the spread.  So \code{.pool_fit()} applies it both to
#' the raw values and to the column centred and scaled
#' (\code{.spread_readings()}), where a second bound keeps the rounding of
#' the values from counting as spread (\code{.scaled_rounding()}), and the
#' test of a constant column of \code{resolution_profile()}'s estimate reads
#' both (\code{.constant_column()}).  The profile's test of a constant
#' response for \code{cp} applies it to the raw values, and where it finds
#' no spread there reads the column again in the same way.  Residuals it
#' reads with \code{.residuals_rounding()} alone.
#'
#' @param ss A sum of squared deviations of \code{y}.
#' @param y The values.
#' @return \code{TRUE} when \code{ss} is at most
#'   \code{100 * .Machine$double.eps^2 * sum(y^2)}.
#' @keywords internal
#' @noRd
.no_spread <- function(ss, y) ss <= 100 * .Machine$double.eps^2 * sum(y^2)


# The rounding a column's values carry, in units in the last place of its
# largest absolute value: a standard deviation of at most this many is no
# spread on the column centred and scaled.  A value stored once is off by
# half a unit at most, and the same number computed two ways (0.1 and
# 0.3 / 3, 2 and sqrt(2)^2, a unit conversion done in two orders) by a unit
# or so: a column of the two gives a standard deviation of 0.5.  A standard
# deviation of 1 at 1e15, where a unit is 0.125, is 8 units, and one of
# 1e-9 at 1e6 is 8.6.
.POOL_ULPS <- 3


# Beyond these a column's squares are not numbers a double holds to its
# precision.  Above 1e150 the sum of the squares of 180 million values
# passes the largest double, and below a spread of 1e-150 the squares of the
# spread fall near the smallest normal double and lose their bits.  A column
# of fewer rows can be fitted a little beyond them; the bounds are set where
# the squares of many values can overflow or lose their precision, not where
# every column does.  The bound above is .POOL_XY_MAX's (R/assignment.R).
.POOL_Y_MAX <- 1e150
.POOL_Y_MIN <- 1e-150

# One of those bounds as a message writes it: "1e150", "1e-150".
.bound_text <- function(x) sub("e+", "e", format(x), fixed = TRUE)


#' A column centred on its median and scaled by its largest distance from it
#'
#' The column the second reading of the test for no spread is taken on, and
#' the one \code{.pool_fit()} fits where the raw values lose their spread to
#' rounding.  The values are divided by the power of two at or below their
#' largest absolute value, the median is subtracted, and the result is
#' divided by the power of two at or below its own largest absolute value,
#' so that the largest lies in [1, 2).  Division by a power of two is exact,
#' and so is the subtraction wherever a value lies within a factor of two of
#' the median, which is where raw values lose their spread.  Equal values
#' stay equal, and a column of one value gives zeros.
#'
#' @param y Finite values, at least one.
#' @return A list: \code{y}, the scaled column; \code{centre}, the median in
#'   the column's units, and \code{scale}, the product of the two divisors,
#'   so that \code{centre + scale * y} gives the values back to rounding;
#'   \code{spread}, the largest distance of a value from the median in the
#'   column's units (0 for a column of one value); and \code{unit}, one unit
#'   in the last place of the largest absolute value, in the units of the
#'   scaled column (\code{.Machine$double.eps} divided by the second divisor).
#' @keywords internal
#' @noRd
.spread_scaled <- function(y) {
  eps <- .Machine$double.eps
  a <- max(abs(y))
  if (!(a > 0)) return(list(y = y, centre = 0, scale = 1, spread = 0, unit = eps))
  # log2() of the largest doubles rounds up to 1024, whose power is Inf, and
  # log2() of a value just below a power of two can round up to it.
  p  <- 2^min(floor(log2(a)), 1023)
  if (p > a) p <- p / 2
  u  <- y / p
  md <- stats::median(u)
  u  <- u - md
  b  <- max(abs(u))
  if (!(b > 0)) return(list(y = u, centre = md * p, scale = p, spread = 0, unit = eps))
  q  <- 2^floor(log2(b))
  if (q > b) q <- q / 2
  list(y = u / q, centre = md * p, scale = p * q, spread = b * p, unit = eps / q)
}


#' Whether a sum of squares of the scaled column is no more than rounding
#'
#' The test of \code{.no_spread()} on the column of \code{.spread_scaled()},
#' and a second bound: the standard deviation the sum gives, over \code{df}
#' degrees of freedom, is at most \code{.POOL_ULPS} units in the last place
#' of the column's largest absolute value.  Centred and scaled, two values
#' one unit in the last place apart are 0 and 1, and without the second
#' bound any two values that differ at all would count as spread: 0.1 beside
#' 0.3 / 3 would be fitted on its last bit, which is what \code{.no_spread()}
#' is there to prevent.
#'
#' @param ss A sum of squared deviations of \code{sc$y}.
#' @param df Its degrees of freedom: the rows less one for the sum about the
#'   overall mean, the rows less the regions for the sum inside them.
#' @param sc \code{.spread_scaled()} of the column.
#' @return \code{TRUE} or \code{FALSE}.
#' @keywords internal
#' @noRd
.scaled_rounding <- function(ss, df, sc)
  .no_spread(ss, sc$y) || ss <= max(df, 0) * (.POOL_ULPS * sc$unit)^2


#' Whether a column has spread, constant or inside its regions, on the
#' column centred and scaled
#'
#' The readings \code{.pool_fit()} takes on its raw values, taken on the
#' column of \code{.spread_scaled()} with the bound of
#' \code{.scaled_rounding()}: no spread over all the rows, and none inside
#' any region, by the sum about the region means or by the same sum about
#' means corrected for the rounding of their own sums (\code{.ss_within()}).
#'
#' @param sc \code{.spread_scaled()} of the column.
#' @param gi,n Region of each row, and rows per region.
#' @return A list: \code{constant} (no spread over all the rows) and
#'   \code{flat} (none inside any region).
#' @keywords internal
#' @noRd
.spread_readings <- function(sc, gi, n) {
  y <- sc$y
  N <- length(y)
  K <- length(n)
  res   <- y - (as.numeric(rowsum(y, gi)) / n)[gi]
  const <- .scaled_rounding(sum((y - mean(y))^2), N - 1L, sc)
  list(constant = const,
       flat = const || .scaled_rounding(sum(res^2), N - K, sc) ||
         .scaled_rounding(.ss_within(res, gi, n), N - K, sc))
}


#' Whether a column's values lie beyond what the arithmetic of the fit holds
#'
#' Read for a column that is not constant: a column of one value needs no
#' squares, at any size.
#'
#' @param y The values.
#' @param sc \code{.spread_scaled(y)}.
#' @return \code{TRUE} when an absolute value is above \code{.POOL_Y_MAX} or
#'   the largest distance from the median is below \code{.POOL_Y_MIN}.
#' @keywords internal
#' @noRd
.pool_out_of_range <- function(y, sc) max(abs(y)) > .POOL_Y_MAX || sc$spread < .POOL_Y_MIN


#' A fit of the scaled column, put back in the column's units
#'
#' The estimator is unchanged by a shift and a scale of the column, so the
#' weights, the p-value of the F test, the centroids, the neighbours, the
#' counts, \code{constant}, \code{flat} and the status of the fit of
#' \code{.spread_scaled()$y} are the column's.  The numbers in squared units
#' are multiplied by the square of the scale and the residuals by the scale,
#' which, both being powers of two, is exact; the targets and the pooled
#' values are \code{centre + scale * x}; the sums and plain means are those
#' of the column itself, and where the value is the plain mean (one region,
#' or no spread inside any region) it is that mean.
#'
#' @param f The list \code{.pool_fit()} returned for the scaled column.
#' @param S,m The column's sums and plain means by region.
#' @param sc \code{.spread_scaled()} of the column.
#' @return \code{f} in the column's units.
#' @keywords internal
#' @noRd
.pool_unscale <- function(f, S, m, sc) {
  s <- sc$scale; s2 <- s * s
  for (v in c("sigma2", "Q", "tau2", "loo_pooled", "loo_pooled_se", "loo_mean",
              "loo_mean_se", "loo_one", "loo_one_se"))
    f[[v]] <- f[[v]] * s2
  f$S <- S
  f$m <- m
  f$target <- sc$centre + s * f$target
  f$value  <- if (f$status %in% c("one cell", "no spread inside any cell")) m
              else sc$centre + s * f$value
  if (!is.null(f$resid_pooled)) f$resid_pooled <- f$resid_pooled * s
  if (!is.null(f$resid_mean))   f$resid_mean   <- f$resid_mean * s
  f
}


#' Whether a column is constant, by the test of the fit
#'
#' What \code{.pool_fit()} returns as \code{constant} for the same values:
#' no spread on the raw values, and none on the column centred and scaled
#' beyond the rounding of its values (\code{.scaled_rounding()}).
#' \code{resolution_profile()} reads it for its estimate's columns and, where
#' the raw values show no spread, for the test of a constant response for
#' \code{cp}; \code{estimate_sac_range()} where its floor finds no spread on
#' the raw response (\code{.regression_centre()}); and the Kish ICC of
#' \code{summarize_by_cell()} for a column whose standard deviation is below
#' \code{.Machine$double.eps} or not finite.
#'
#' @param y The values, at least one.
#' @param sc \code{.spread_scaled(y)}.
#' @return \code{TRUE} or \code{FALSE}.
#' @keywords internal
#' @noRd
.constant_column <- function(y, sc = .spread_scaled(y))
  .no_spread(sum((y - mean(y))^2), y) &&
    .scaled_rounding(sum((sc$y - mean(sc$y))^2), length(y) - 1L, sc)


#' One unit in the last place of the largest absolute value of a column
#'
#' @param y The values.
#' @return The power of two at or below \code{max(abs(y))} times
#'   \code{.Machine$double.eps}; 0 for a column of zeros.
#' @keywords internal
#' @noRd
.ulp_of <- function(y) {
  a <- max(abs(y))
  if (!(a > 0)) return(0)
  # As in .spread_scaled(): log2() of the largest doubles rounds up to 1024,
  # and of a value just below a power of two can round up to it.
  p <- 2^min(floor(log2(a)), 1023)
  if (p > a) p <- p / 2
  p * .Machine$double.eps
}


#' The value a response is regressed about where its raw values lose their
#' spread
#'
#' A response at a large offset (1e15 plus a spread of 1) is stored to an
#' eighth, and a regression on an intercept of the raw values works on the
#' rounding of the offset, which there is larger than the spread: OLS
#' residuals with a sum of squares of 90,000 where the spread gives 614.
#' Residuals, Moran's I of cell means and the test for spatial structure do
#' not change under a shift, so such a response is regressed about the
#' centre of \code{.spread_scaled()}.  That is where its raw values show no
#' spread (a standard deviation of at most 1e-12 of the largest absolute
#' value, \code{estimate_sac_range()}'s floor) and the column centred and
#' scaled does (\code{.constant_column()}).  Everywhere else the centre is
#' 0, so a response whose raw values show their spread is computed as
#' before, and so is a constant one.
#'
#' @param y The values; non-finite ones are left out.
#' @param bounded \code{TRUE} to return 0 also for values beyond what the
#'   pooling's arithmetic holds (\code{.pool_out_of_range()}), whose columns
#'   \code{resolution_profile()} leaves \code{NA} for that reason.
#' @return The centre, in the units of \code{y}, or 0.
#' @keywords internal
#' @noRd
.regression_centre <- function(y, bounded = TRUE) {
  y <- y[is.finite(y)]
  if (length(y) < 2L) return(0)
  a <- max(abs(y))
  if (!(a > 0) || isTRUE(stats::sd(y / a) > 1e-12)) return(0)
  sc <- .spread_scaled(y)
  if (.constant_column(y, sc) || (bounded && .pool_out_of_range(y, sc))) return(0)
  sc$centre
}


#' The value a response is regressed about again where its residuals lose
#' their spread
#'
#' A response whose raw values show their spread through a trend (1e15 plus
#' 1e5 times a predictor plus a spread of 1) is regressed on its raw values,
#' and the regression works on the rounding of the offset: its residuals,
#' which should hold the spread of 1, hold that rounding as well, and the
#' bound of the regression's rounding on the raw values holds them both.
#' Where the residuals of the raw values show no spread, the response is
#' regressed again about the centre of \code{.spread_scaled()}, and the
#' residuals of that fit are read.  Residuals do not change under a shift,
#' so this changes nothing but the rounding.
#'
#' @param y The values; non-finite ones are left out.
#' @param bounded \code{TRUE} to return 0 for values beyond what the
#'   pooling's arithmetic holds (\code{.pool_out_of_range()}), as
#'   \code{.regression_centre()} does.
#' @return The centre, in the units of \code{y}, or 0.
#' @keywords internal
#' @noRd
.residual_centre <- function(y, bounded = TRUE) {
  y <- y[is.finite(y)]
  if (length(y) < 2L || !(max(abs(y)) > 0)) return(0)
  sc <- .spread_scaled(y)
  if (bounded && .pool_out_of_range(y, sc)) return(0)
  sc$centre
}


# A predictor whose centre is more than this many times its spread is
# taken about its centre before it is regressed on.  The design cbind(1, x)
# has a condition number of about |centre| / spread, and the rounding of
# the fit grows with it: at 1e7 a predictor that explains most of the
# response moves the residuals by 1e-4 of their size, and from about 1e7
# lm.fit() drops it.  Below the bound the raw column is fitted as it is.
.PRED_CENTRE_RATIO <- 1e3


#' Predictor columns at an offset taken about their centres
#'
#' The QR decomposition behind \code{lm.fit()} and \code{lm()} drops a
#' column whose part orthogonal to the earlier columns is under 1e-7 of its
#' norm.  A predictor whose spread is under about 1e-7 of its size (a
#' timestamp in seconds, a northing, 1e8 plus a spread of 1) is then taken
#' as collinear with the intercept and left out, and the residuals are those
#' of the intercept alone; a little below that the column is kept, but the
#' fit's rounding grows with the ratio of its size to its spread.  Taken
#' about its centre the same column keeps its spread, and the OLS fit of a
#' model with an intercept does not change under a shift of a predictor.
#' A column is taken less the centre of \code{.spread_scaled()} of its
#' finite values where that centre is more than \code{.PRED_CENTRE_RATIO}
#' times the column's spread; every other column is returned as it is, so a
#' fit on such columns gives what it gave to the last bit.  A column of
#' times, dates or time differences (\code{POSIXt}, \code{Date},
#' \code{difftime}), which \code{lm()} reads as numbers, is read as the
#' same numbers.  A column of one value is returned as it is, and is
#' dropped by the fit as before.
#'
#' @param P A numeric matrix or a data frame of predictor columns.
#' @return \code{P} with each column at an offset less its centre, as a
#'   plain numeric column.  Where a column is taken less its centre,
#'   \code{attr(, "centres")} holds the centre of each column (0 for one
#'   left as it is), for \code{.rounding_design()}.
#' @keywords internal
#' @noRd
.centre_columns <- function(P) {
  centres <- numeric(ncol(P))
  for (j in seq_len(ncol(P))) {
    v <- if (is.data.frame(P)) P[[j]] else P[, j]
    if (inherits(v, c("POSIXt", "Date", "difftime"))) v <- as.numeric(v)
    if (!is.numeric(v) || is.matrix(v)) next
    f <- v[is.finite(v)]
    if (length(f) < 2L || !(max(abs(f)) > 0)) next
    sc <- .spread_scaled(f)
    if (!(sc$spread > 0) || !(abs(sc$centre) > .PRED_CENTRE_RATIO * sc$spread)) next
    if (is.data.frame(P)) P[[j]] <- v - sc$centre else P[, j] <- v - sc$centre
    centres[j] <- sc$centre
  }
  if (any(centres != 0)) attr(P, "centres") <- centres
  P
}


#' The design and coefficients the test for an exact fit reads
#'
#' \code{.residuals_rounding()} takes the rounding of an exact fit on the
#' size of the fit's terms.  A response computed from a predictor as it is
#' given carries rounding of the size of that predictor's terms (5 times a
#' northing of 4.5e6 is 2.25e7, on a layer 4,000 m across), and the terms
#' of the design taken about its centres (\code{.centre_columns()}) are
#' thousands of times smaller.  So the test reads the design with the
#' centres added back, and the intercept less the centres times their
#' coefficients, which give the same fitted values.  Where that design does
#' not have full rank (\code{qr()}, the tolerance \code{lm.fit()} reads),
#' the offset drops a column from it, its terms are not those of a fit, and
#' the design as fitted is read.
#'
#' @param X The design as fitted, the intercept's column of ones first.
#' @param b Its coefficients.
#' @param centres The centre taken off each column of \code{X}, 0 for the
#'   intercept and for a column left as it is; \code{NULL} for none.
#' @return A list: \code{X} and \code{b}.
#' @keywords internal
#' @noRd
.rounding_design <- function(X, b, centres) {
  if (is.null(centres) || !any(centres != 0)) return(list(X = X, b = b))
  Xr <- X + rep(centres, each = nrow(X))
  if (qr(Xr)$rank < ncol(Xr)) return(list(X = X, b = b))
  bf <- as.numeric(b)
  bf[!is.finite(bf)] <- 0
  b[1L] <- b[1L] - sum(centres * bf)
  list(X = Xr, b = b)
}


#' An OLS fit on the predictors, those at an offset taken about their
#' centres
#'
#' \code{lm.fit(X, y)} on the design with its predictor columns passed
#' through \code{.centre_columns()}.  A design with no column at an offset
#' is fitted as it is, to the last bit.
#'
#' @param X The design, the intercept's column of ones first.
#' @param y The response.
#' @return A list: \code{fit}, the \code{lm.fit()} result or the
#'   \code{try-error}; \code{X}, the design it was fitted on; and
#'   \code{centres}, the centre taken off each of its columns, or
#'   \code{NULL} where none was (for \code{.rounding_design()}).
#' @keywords internal
#' @noRd
.lm_fit_centred <- function(X, y) {
  centres <- NULL
  if (ncol(X) >= 2L) {
    P <- .centre_columns(X[, -1L, drop = FALSE])
    if (!is.null(attr(P, "centres"))) centres <- c(0, attr(P, "centres"))
    X[, -1L] <- P
  }
  list(fit = try(stats::lm.fit(x = X, y = y), silent = TRUE), X = X, centres = centres)
}


# The rounding of an exact fit's residuals, in .Machine$double.eps times the
# size of the fit's terms (see .residuals_rounding()): the multiple of
# .no_spread(), whose bound is 100 eps^2 times the sum of squares.
.FIT_EPS <- 10


#' Whether residuals are no more than the rounding of the response they are
#' of
#'
#' The residuals of a response the predictors fit exactly are rounding: of
#' the values, a few units in the last place of the raw response, and of
#' the regression.  Either bound holds them.  The first is
#' \code{.POOL_ULPS} units in the last place of the raw response's largest
#' absolute value, as a standard deviation over the residuals' degrees of
#' freedom.  The second is taken on the size of the fit's terms: each row's
#' absolute intercept plus the absolute products of its predictors and
#' their coefficients.  An exact fit's response is those terms summed, and
#' its rounding is at most a few \code{.Machine$double.eps} times their size,
#' whatever the number of rows.  The residuals an OLS fit returns carry
#' more, and more with more rows: \code{lm.fit()} leaves those of 7 w - 2 w2
#' on 800 rows at 13 \code{.Machine$double.eps} times the root mean square
#' of the terms, and those of 2.5 x on 30,000 rows of an integer x at 600
#' times.  So the residuals are taken again from the coefficients,
#' \code{yd - X b}, less their own OLS fit on \code{X}, which leaves the
#' rounding of the response and of that sum.  On exact fits with one to five
#' predictors, 30 to 100,000 rows, at offsets 0 and 1e6, raw and about the
#' response's centre, those residuals that the first bound does not hold
#' come out at most 0.37 \code{.Machine$double.eps} times the root mean
#' square of the terms.  The bound is a root mean square of
#' \code{.FIT_EPS} (10) \code{.Machine$double.eps} times that of the terms,
#' the multiple \code{.no_spread()} takes on a column, and residuals above
#' it count as spread: those of 1e3 w + 1e-11 f on w, a field f, are 45 times.
#' The sums are taken on the residuals and terms divided by the power of two
#' at or below the largest term or value, which is exact, so that they
#' neither overflow nor underflow.
#'
#' @param e The residuals the fit returned.
#' @param y The raw response.
#' @param yd The response as regressed (centred by
#'   \code{.regression_centre()} or \code{.residual_centre()}, or raw).
#' @param X The fit's design, with the intercept's column of ones; its rows
#'   are those of \code{yd}.  \code{NULL} where the design is not known,
#'   and the regression's bound is then that of \code{.no_spread()} on
#'   \code{yd}, read on \code{e}.
#' @param b The fit's coefficients; \code{NA} (an aliased column) counts as
#'   0.
#' @return \code{TRUE} or \code{FALSE}.
#' @keywords internal
#' @noRd
.residuals_rounding <- function(e, y, yd, X = NULL, b = NULL) {
  ss <- function(r) sum((r - mean(r))^2)
  if (is.null(X) || is.null(b))
    return(ss(e) <= max(100 * .Machine$double.eps^2 * sum(yd^2),
                        (length(e) - 1L) * (.POOL_ULPS * .ulp_of(y))^2))
  X <- as.matrix(X)
  b <- as.numeric(b)
  b[!is.finite(b)] <- 0
  size <- drop(abs(X) %*% abs(b))
  a <- max(size, abs(yd))
  # The power of two at or below it, so that the division is exact.
  a <- if (is.finite(a) && a > 0) 2^min(floor(log2(a)), 1023) else 1
  if (a > max(size, abs(yd))) a <- a / 2
  # .POOL_ULPS units in the last place of the raw response, as a standard
  # deviation over n - 1 degrees of freedom, on the scale of `a`.
  ulps <- function(r) (length(r) - 1L) * (.POOL_ULPS * .ulp_of(y) / a)^2
  if (ss(e / a) <= ulps(e)) return(TRUE)
  r <- yd / a - drop(X %*% (b / a))
  fit <- try(stats::lm.fit(X, r), silent = TRUE)
  if (!inherits(fit, "try-error")) {
    b2 <- fit$coefficients
    b2[!is.finite(b2)] <- 0
    r <- r - drop(X %*% b2)
  }
  ss(r) <= max(ulps(r), (.FIT_EPS * .Machine$double.eps)^2 * sum((size / a)^2))
}


#' The mean square of a set of residuals and its standard error over rows
#'
#' \code{c(mean(e^2), sd(e^2) / sqrt(length(e)))}.  The standard error
#' squares the squares, which pass the largest double for residuals above
#' about 1e77, where it came out \code{Inf}, and fall below the smallest
#' normal double for residuals under about 1e-77, where it lost its
#' precision and further down came out 0.  Residuals whose largest lies
#' outside 2^-255 to 2^255 are therefore divided by the power of two at or
#' below it first, and both numbers multiplied back by its square, which is
#' exact; inside those bounds the numbers are computed as they stand.
#'
#' @param e Finite residuals, at least two.
#' @return Two numbers.
#' @keywords internal
#' @noRd
.loo_mean_square <- function(e) {
  a <- max(abs(e))
  if (!isTRUE(a > 0) || (a >= 2^-255 && a <= 2^255))
    return(c(mean(e^2), stats::sd(e^2) / sqrt(length(e))))
  s  <- 2^floor(log2(a))
  e2 <- (e / s)^2
  c(mean(e2), stats::sd(e2) / sqrt(length(e))) * s^2
}


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
#' searched, and the result mapped back to the regions' numbers.  Two
#' candidates whose squared distances differ by no more than \code{.TIE_TOL}
#' (one part in ten million) are a tie, and a tie goes to the lower first
#' coordinate, then the lower second, first coordinates within that part of
#' the distance counting as equal (\code{.knn_pick_rows()}), whatever the
#' regions are numbered.  So the last bits of the centroids, which move with
#' the unit and the origin of the coordinates, do not decide.  Of candidates
#' still tied, second coordinates within that part of the distance counting
#' as equal too, the region of the lower number is taken (\code{num = o}),
#' not the lower place in coordinate order.  So regions whose centroids
#' coincide are told apart by their number whether they coincide exactly or
#' only to rounding, as another unit or origin can leave them; that is the
#' one place the numbering decides anything.  A region is never its own
#' neighbour.
#'
#' A centroid that shares its point with \code{k} others or more is not
#' searched for at all: its neighbours are the first \code{k} of them by
#' region number.  It keeps a layer whose regions lie on one point by the
#' thousand (rows that carry the coordinates of a coarser unit, or a
#' placeholder) from costing the square of their number: a kd-tree asked
#' about such a centroid has to list every other one on its point before it
#' can tell which come first.  A point is found with the reach of a tie,
#' \code{.TIE_TOL} times the largest absolute coordinate of the centroids
#' (\code{.shared_points()}): centroids that chains of pairs within it in x
#' and in y link are on one point while they lie within it of one another,
#' and so are equal centroids; a chain that spans further is linked again at
#' half the reach, and so on.  Another unit, or an origin moved by little
#' against the coordinates, parts the centroids of a point by rounding, far
#' less than the reach, so they stay on it and each keeps the same
#' \code{k}.  The reach follows the largest coordinate, so an origin moved
#' far, to the corner of the centroids say, can part centroids that lie
#' within it without coinciding.  Where more than \code{k} are on a point,
#' centroids on it that do not coincide take each other by number too, not
#' by distance.  With \code{k} or fewer on a point they are searched for;
#' where the point's centroids coincide each takes all the others, and only
#' the order of its neighbours can follow the rounding.
#'
#' For the rest, with \pkg{FNN} the kd-tree returns one neighbour more than
#' is wanted, where there are that many other regions (ten regions or more
#' for eight neighbours), and R recomputes the squared distances.  A row is
#' settled when it is not among its own neighbours and its squared distances
#' rise by more than \code{.TIE_TOL} and one part in a thousand million from
#' each neighbour to the next; it keeps the tree's first \code{k}, which is
#' what the rule picks.  Every other row is searched again by
#' \code{.knn_tied_rows()} with the rule of the exact search.  Without
#' \pkg{FNN} each of them is searched exactly, which compares it with every
#' centroid; the pairs that makes, the rows searched for times the centroids,
#' are what \code{.work_guard()} stops above the limit.  The two paths return
#' the same matrix.
#'
#' @param cen Two-column matrix of centroids, one row per region.
#' @param k Number of neighbours wanted.
#' @param use_fnn Whether to use \pkg{FNN}.
#' @param caller Name of the exported function, for the work guard's message.
#' @return An integer matrix with a row per region and \code{min(k, K - 1)}
#'   columns, nearest first to within \code{.TIE_TOL}.
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
  # Centroids within the reach of a tie of one another are on one point.
  h <- .TIE_TOL * max(abs(u))
  if (!is.finite(h)) h <- 0
  pt    <- .shared_points(x, y, h)
  many  <- tabulate(pt)[pt] > k             # k others or more on its own point
  idx   <- matrix(0L, K, k)
  if (any(many)) {
    # The point's centroids by number; each takes the first k + 1 of them,
    # less itself when it is one.
    at   <- which(many)
    m    <- at[order(pt[at], o[at])]
    from <- match(pt[m], pt[m])             # where each point starts in m
    pos  <- integer(K); pos[m] <- seq_along(m)
    f    <- from[pos[at]]
    j    <- matrix(rep(seq_len(k) - 1L, each = length(at)), length(at), k)
    idx[at, ] <- m[f + j + (j >= pos[at] - f)]
  }
  rest <- which(!many)
  if (!isTRUE(use_fnn)) {
    .work_guard(as.numeric(length(rest)) * K, .exact_pairs_max(), what = "cell neighbours",
                caller = caller,
                search = sprintf("finding each cell's neighbours among %d cells", K))
    idx[rest, ] <- .knn_exact_rows(x, y, rest, k, tol = .TIE_TOL, num = o)
  } else if (length(rest)) {
    kq  <- min(k + 1L, K - 1L)
    nn  <- matrix(FNN::get.knn(u, k = kq)$nn.index, K, kq)[rest, , drop = FALSE]
    d2  <- matrix((x[nn] - x[rest])^2 + (y[nn] - y[rest])^2, length(rest), kq)
    own <- rowSums(nn == rest) > 0L
    tied <- if (kq > 1L)
      rowSums(d2[, -1L, drop = FALSE] <= d2[, -kq, drop = FALSE] * (1 + .TIE_TOL + 1e-9)) > 0L
    else rep(FALSE, length(rest))
    idx[rest, ] <- nn[, seq_len(k), drop = FALSE]
    redo <- rest[own | tied]
    if (length(redo)) idx[redo, ] <- .knn_tied_rows(u, redo, k, tol = .TIE_TOL, num = o)
  }
  out <- matrix(0L, K, k)
  out[o, ] <- o[idx]
  out
}


#' Centroids on one point, to within a reach
#'
#' Centroids that chains of pairs within \code{h} of each other in x and in
#' y link are one point, as long as they all lie within \code{h} of one
#' another in x and in y.  A chain that spans further is not a point: its
#' centroids are linked again at \code{h / 2}, the chains of that wider than
#' \code{h / 2} again at \code{h / 4}, and so on, until no chain is wider
#' than the reach it was linked at.  So rounding that moves the centroids of
#' a point by far less than \code{h} keeps them on it where other centroids
#' that close chain on beyond \code{h}, and a point never spans more than
#' \code{h}.  Each halving parts at least one pair of the chain's centroids
#' once the reach is below their distance, so it ends; centroids it leaves
#' apart to the last are on the point of the centroids equal to them.  With
#' \code{h = 0} a point is the centroids that are equal.
#'
#' @param x,y Coordinates of the centroids.
#' @param h The reach, 0 or more.
#' @return An integer vector: the point of each centroid, numbered from 1.
#' @keywords internal
#' @noRd
.shared_points <- function(x, y, h) {
  n <- length(x)
  if (n == 0L) return(integer(0))
  o  <- order(x, y)
  xs <- x[o]; ys <- y[o]
  first <- c(TRUE, xs[-1L] != xs[-n] | ys[-1L] != ys[-n])
  d  <- cumsum(first)                       # distinct centroid of each, in order
  nd <- d[n]
  xd <- xs[first]; yd <- ys[first]
  pt <- seq_len(nd)                         # equal centroids, until a reach joins more
  at <- if (isTRUE(h > 0) && nd >= 2L) seq_len(nd) else integer(0)
  ch <- rep(1L, length(at))                 # the chain each was on at the last reach
  top <- nd
  while (length(at) > 1L && h > 0) {
    ch <- .reach_chains(xd[at], yd[at], ch, h)
    wide <- .chain_span(xd[at], ch) > h | .chain_span(yd[at], ch) > h
    pt[at[!wide]] <- top + ch[!wide]
    top <- top + length(at)
    at <- at[wide]; ch <- ch[wide]
    h  <- h / 2
  }
  out <- integer(n); out[o] <- pt[d]
  match(out, unique(out))
}


#' How far each centroid's chain spans along one coordinate
#'
#' @param v One coordinate of the centroids.
#' @param g The chain of each, as positive integers.
#' @return For each centroid, the greatest less the least \code{v} of its
#'   chain.
#' @keywords internal
#' @noRd
.chain_span <- function(v, g) {
  lo <- hi <- numeric(max(g))
  ov <- order(g, v)
  lo[rev(g[ov])] <- rev(v[ov]); hi[g[ov]] <- v[ov]
  hi[g] - lo[g]
}


#' Chains of distinct centroids within a reach of each other
#'
#' Centroids within \code{h} of each other in x and in y are linked, and a
#' chain is the centroids that links join, directly or through others.  Only
#' centroids of one group (\code{g}, the chain each was on at a wider reach)
#' are linked.  Each group is put in squares of side \code{h} counted from
#' its own least x and y, so the squares' numbers stay small however small
#' \code{h} gets: a group from a wider reach is a chain at twice \code{h},
#' and spans at most twice \code{h} times its centroids.  Two in one square are
#' within \code{h} of each other, and two within \code{h} lie in one square
#' or in two that touch.  Squares side by side are linked when the nearest of
#' their centroids along that side are within \code{h}; squares corner to
#' corner are linked when a pair of their centroids is.
#'
#' @param x,y Coordinates of distinct centroids.
#' @param g The group of each.
#' @param h The reach, above 0.
#' @return An integer vector: the chain of each centroid, as the place in
#'   \code{x} of one centroid of it.
#' @keywords internal
#' @noRd
.reach_chains <- function(x, y, g, h) {
  nd <- length(x)
  g  <- match(g, unique(g))
  lo_x <- numeric(max(g)); lo_y <- numeric(max(g))
  ox <- order(g, x); oy <- order(g, y)
  lo_x[rev(g[ox])] <- rev(x[ox]); lo_y[rev(g[oy])] <- rev(y[oy])
  gx <- floor((x - lo_x[g]) / h) + 1; gy <- floor((y - lo_y[g]) / h) + 1
  nx <- max(gx) + 2; ny <- max(gy) + 2
  key <- if ((max(g) + 1) * nx * ny < 2^52)
    function(i, j) (g * nx + i) * ny + j
  else function(i, j) paste(g, i, j)
  gk   <- key(gx, gy)
  cell <- match(gk, gk)                     # the first centroid of each square
  # The least and the greatest x and y of each square.
  lo_x <- hi_x <- lo_y <- hi_y <- numeric(nd)
  ox <- order(cell, x); oy <- order(cell, y)
  lo_x[rev(cell[ox])] <- rev(x[ox]); hi_x[cell[ox]] <- x[ox]
  lo_y[rev(cell[oy])] <- rev(y[oy]); hi_y[cell[oy]] <- y[oy]
  sq  <- which(cell == seq_len(nd))
  mem <- split(seq_len(nd), cell)
  e1 <- e2 <- integer(0)
  for (dir in list(c(1, 0), c(0, 1), c(1, 1), c(1, -1))) {
    nb <- cell[match(key(gx + dir[1L], gy + dir[2L])[sq], gk)]
    a  <- sq[!is.na(nb)]; b <- nb[!is.na(nb)]
    if (!length(a)) next
    linked <- if (dir[2L] == 0) lo_x[b] - hi_x[a] <= h
    else if (dir[1L] == 0) lo_y[b] - hi_y[a] <= h
    else {
      ma <- mem[as.character(a)]; mb <- mem[as.character(b)]
      la <- lengths(ma); lb <- lengths(mb)
      pa <- unlist(mapply(function(u, v) rep(u, length(v)), ma, mb, SIMPLIFY = FALSE),
                   use.names = FALSE)
      pb <- unlist(mapply(function(u, v) rep(v, each = length(u)), ma, mb, SIMPLIFY = FALSE),
                   use.names = FALSE)
      ok <- abs(x[pb] - x[pa]) <= h & abs(y[pb] - y[pa]) <= h
      seq_along(a) %in% rep(seq_along(a), la * lb)[ok]
    }
    e1 <- c(e1, a[linked]); e2 <- c(e2, b[linked])
  }
  # Join the linked squares: each points to the lowest square of its chain.
  lab <- cell
  repeat {
    r1 <- lab[e1]; r2 <- lab[e2]
    hi <- pmax(r1, r2); lo <- pmin(r1, r2)
    if (!any(hi != lo)) break
    k  <- order(hi, lo)
    f  <- !duplicated(hi[k])
    lab[hi[k][f]] <- pmin(lab[hi[k][f]], lo[k][f])
    repeat {
      l2 <- lab[lab]
      if (identical(l2, lab)) break
      lab <- l2
    }
  }
  lab
}


#' Centroids of the regions and each region's neighbours
#'
#' The centroid of a region is the mean of the coordinates of its rows.  The
#' rows are summed in coordinate order, so the centroids are the same numbers
#' to the last bit whatever order the rows come in.  On gridded locations
#' that keeps the neighbours fixed under another order of the rows: summed in
#' row order the centroids of a lattice differ in their last bits from one
#' row order to the next.  Under another unit or origin the centroids differ
#' in their last bits all the same, and there the tie tolerance of
#' \code{.region_neighbours()} keeps the neighbours fixed.
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
#' Two locations whose squared distances differ by no more than
#' \code{.TIE_TOL} are a tie, as in \code{.region_neighbours()}.  With
#' \pkg{FNN} the kd-tree returns the three nearest locations of each one
#' asked for, itself among them, and R recomputes the squared distances.  A
#' location is settled when it is one of its own three and the squared
#' distances of its two nearest others differ by more than \code{.TIE_TOL}
#' and one part in a thousand million.  Every other location is searched
#' again by \code{.knn_tied_rows()}, which settles a tie by the rule, so a
#' lattice costs no comparison of every pair.  Without \pkg{FNN} the search
#' is exact and is stopped by \code{.work_guard()} above the limit.  The two
#' paths return the same locations.
#'
#' @param u Two-column matrix of distinct locations in coordinate order, at
#'   least two.
#' @param q The locations to search for, as rows of \code{u}.
#' @param use_fnn Whether to use \pkg{FNN}.
#' @param caller Name of the exported function, for the work guard's message.
#' @return For each of \code{q} the row of \code{u} nearest it and not it; a
#'   tie goes to the lower first coordinate, then the lower second (first
#'   coordinates within \code{.TIE_TOL} of the distance counting as equal),
#'   which is the lower row.
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
    return(.knn_exact_rows(u[, 1L], u[, 2L], q, 1L, tol = .TIE_TOL)[, 1L])
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
  redo <- which(rowSums(own) != 1L | !(second > best * (1 + .TIE_TOL + 1e-9)))
  if (length(redo)) out[redo] <- .knn_tied_rows(u, q[redo], 1L, tol = .TIE_TOL)[, 1L]
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
#'     other location; a tie between locations, squared distances within
#'     \code{.TIE_TOL} of each other, goes to the lower first coordinate,
#'     then the lower second (\code{.nearest_location_search()}).
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
#'   \item{\code{"values out of range"}}{The column is not constant, and an
#'     absolute value is above \code{.POOL_Y_MAX} or the largest distance of
#'     a value from the median is below \code{.POOL_Y_MIN}: the squares the
#'     fit is made of are not numbers a double holds to its precision.
#'     Nothing is fitted, and what is returned is as under \code{"no cell
#'     holds two values"}.}
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
#'     the rounding of their own sums (\code{.ss_within()}), on the raw values
#'     and on the column centred and scaled alike, where a standard deviation
#'     of at most \code{.POOL_ULPS} units in the last place is rounding too
#'     (below).  \code{sigma2} is 0, every weight 1 and the value the plain
#'     mean; the residual of a row that shares its region is 0.  When the column is constant over all
#'     its rows as well (\code{constant}), \code{tau2} is 0, every residual
#'     is 0 and every leave-one-out error is \code{NA}: no level fits a
#'     constant better than another.  These are set, not computed, because
#'     the formulas would divide one rounding error by another.}
#'   \item{\code{"no signal between cells"}}{\code{tau2} is 0: every weight
#'     is 0 and each value is its target.}
#'   \item{\code{"ok"}}{None of the above.}
#' }
#'
#' No spread is read twice: on the raw values, and on the column centred on
#' its median and scaled by its largest distance from it
#' (\code{.spread_scaled()}, \code{.spread_readings()}).  \code{constant}
#' and \code{flat} are \code{TRUE} where neither reading finds spread.  The
#' raw reading compares the sums of squares with the size of the values, so
#' 1e15 plus a spread of 1, whose values are stored to an eighth, has none by
#' it, and values whose squares overflow or underflow have none at any
#' spread; the scaled reading compares them with the size of the
#' differences, which no shift or scale of the column changes, and with the
#' rounding of the values: a standard deviation, over all the rows or inside
#' the regions, of at most \code{.POOL_ULPS} (3) units in the last place of
#' the largest absolute value is no spread by it either
#' (\code{.scaled_rounding()}).  So a column that varies only in its last
#' bits, 0.1 beside 0.3 / 3, has no spread by either reading, while 1e15 plus
#' a spread of 1 (8 such units) has spread by the scaled one.  Where the
#' scaled reading finds spread the raw one does not, the raw arithmetic has
#' lost it as well, and the fit is made on the scaled column, where the
#' differences are exact, with the column's readings, and put back in the
#' column's units (\code{.pool_unscale()}).  Every other column is fitted on
#' its raw values, so where the raw values show their spread the scaled
#' reading changes no number.
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
#' @param readings \code{NULL}, or for the fit of the scaled column the
#'   function makes of itself, \code{list(constant, flat)} as it read them on
#'   the column: the values are then neither read, refused nor scaled again.
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
                      caller = "resolution_profile", readings = NULL) {
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
  sc <- NULL
  rescale <- FALSE
  if (is.null(readings)) {
    const <- .no_spread(tss, y)                     # constant over the rows
    # No spread inside any region.  The second test is on the sum of squares
    # the fit uses; the third is on the same sum about means that carry no
    # rounding of their own, without which a column constant inside large
    # regions fails the second (.ss_within()).
    flat  <- const || .no_spread(ss_w, y) || .no_spread(.ss_within(res, gi, n), y)
    # The same readings on the column centred and scaled, where the raw ones
    # find no spread: they compare the sums with the size of the values, so
    # 1e15 plus a spread of 1 has none by them, and values whose squares
    # overflow or underflow have none at any spread.  There a standard
    # deviation of .POOL_ULPS units in the last place or less is rounding,
    # so a column that varies only in its last bits (0.1 beside 0.3 / 3)
    # stays without spread.  A column has no spread where neither reading
    # finds any.  Where the scaled reading finds spread the raw one misses,
    # the raw arithmetic has lost it too, and the fit is made on the scaled
    # column (below); every other column is fitted on its raw values.
    sc <- .spread_scaled(y)
    if (flat) {
      on <- .spread_readings(sc, gi, n)
      rescale <- (const && !on$constant) || !on$flat
      const <- const && on$constant
      flat  <- on$flat
    }
  } else {
    # The fit of the scaled column: the readings are those of the column.
    const <- readings$constant
    flat  <- readings$flat
  }
  # Mean square of a set of residuals, and its standard error over rows.
  loo_of <- .loo_mean_square
  none   <- c(NA_real_, NA_real_)

  out <- list(n = n, S = S, m = m, N = N, K = K, sigma2 = NA_real_,
              cen = NA_real_, nb = NA_integer_, target = rep(NA_real_, K),
              Q = NA_real_, tau2 = NA_real_, w = rep(NA_real_, K),
              value = rep(NA_real_, K), resid_pooled = NULL, resid_mean = NULL,
              loo_pooled = NA_real_, loo_pooled_se = NA_real_,
              loo_mean = NA_real_, loo_mean_se = NA_real_,
              loo_one = NA_real_, loo_one_se = NA_real_, f_p = NA_real_,
              constant = const, flat = flat, status = NA_character_)

  # ---- 0. Beyond the arithmetic, or spread lost to it ----------------------
  if (!is.null(sc) && !const && .pool_out_of_range(y, sc)) {
    out$status <- "values out of range"
    return(out)
  }
  # The fit of the scaled column takes the column's readings, and is not
  # read or scaled again.
  if (rescale)
    return(.pool_unscale(.pool_fit(sc$y, gi, K, xy, want = want, nn1 = nn1, geo = geo,
                                   use_fnn = use_fnn, cache = cache, caller = caller,
                                   readings = list(constant = const, flat = flat)),
                         S, m, sc))

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
#'   \code{"variogram"} or a number.  It reads \code{1}, alone or before a
#'   semicolon that names the columns whose intervals are \code{NA}, exactly
#'   when the interval carries none, which is what \code{.deff_label()}
#'   writes.
#' @return A data frame with the columns \code{columns} (the stems of the two
#'   interval columns), \code{covers}, \code{centre}, \code{level} and
#'   \code{method}; no row when both levels are \code{NULL}.
#' @keywords internal
#' @noRd
.intervals_table <- function(obs_level, conf_level, deff_label) {
  obs <- !is.null(obs_level)
  ci  <- !is.null(conf_level)
  # A 1 that goes on to name the columns left NA ("1; NA for 'z', ...")
  # describes the default's interval for every other column.
  own <- grepl("^1(;|$)", as.character(deff_label))
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
