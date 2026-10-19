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


#' Mean correlation between two random points in a polygon
#'
#' The domain term of Krige's relation, \eqn{\bar r(V)}, over the study area
#' itself rather than its bounding box: quadrature on a regular grid of about
#' \eqn{q^2} points clipped to \code{domain}.  The profile measures the
#' domain's area on the convex hull, and a bounding box that the hull fills
#' only partly (a diagonal strip, a triangle) spreads the pairs further apart
#' than the domain does, which made \eqn{\bar r(V)} too small, the
#' between-cell variance too large and reliability too high at coarse levels.
#' Measured on 600 points in a 3000 x 120 strip, reliability picked 8 cells
#' axis-aligned and 2 rotated by 45 degrees (10 and 8 on a square), while
#' the hull's area did not move; over the hull both orientations pick 6
#' (10 on the square).  The grid is isotropic, where \code{.rbar_rect()}'s
#' \eqn{q \times q} grid on a thin rectangle is 25 times coarser along it
#' than across.
#'
#' @param cor_fn Correlation function of distance.
#' @param domain A polygon (sfc) in the coordinate units of \code{cor_fn}.
#' @param area Its area.
#' @param q Grid points per side of a square of the same area.
#' @return A scalar in \eqn{[0, 1]}; the bounding box's value when too few
#'   grid points fall inside (a domain with next to no area).
#' @keywords internal
#' @noRd
.rbar_domain <- function(cor_fn, domain, area, q = 20L) {
  bb <- sf::st_bbox(domain)
  w  <- as.numeric(bb[["xmax"]] - bb[["xmin"]]); h <- as.numeric(bb[["ymax"]] - bb[["ymin"]])
  if (!is.finite(area) || area <= 0 || !is.finite(w) || !is.finite(h) || w <= 0 || h <= 0)
    return(.rbar_rect(cor_fn, w, h, q))
  # Spacing for about q^2 points inside the domain, but never more than 40
  # q^2 grid nodes over the box, so a sliver of a hull costs no more than a
  # point-in-polygon test on 16,000 points.
  s  <- max(sqrt(area) / q, sqrt(w * h / (40 * q^2)))
  gx <- seq(bb[["xmin"]] + s / 2, bb[["xmax"]], by = s)
  gy <- seq(bb[["ymin"]] + s / 2, bb[["ymax"]], by = s)
  g  <- as.matrix(expand.grid(x = gx, y = gy))
  inside <- lengths(sf::st_intersects(
    sf::st_as_sf(as.data.frame(g), coords = c("x", "y"), crs = sf::st_crs(domain)),
    domain)) > 0L
  if (sum(inside) < 10L) return(.rbar_rect(cor_fn, w, h, q))
  r <- as.numeric(cor_fn(as.numeric(stats::dist(g[inside, , drop = FALSE]))))
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


# -----------------------------------------------------------------------------
# A ladder set by a size, the estimate's columns, and the two bounds on cost
# -----------------------------------------------------------------------------

#' A refused value in the words of an error message
#'
#' @param x The value the caller passed.
#' @return One string: a string in its quotes, numbers as they print (the
#'   first six of them), \code{numeric(0)} for a vector of no length, and the
#'   class of anything else.
#' @keywords internal
#' @noRd
.got_text <- function(x) {
  if (is.null(x)) return("NULL")
  if (!(is.character(x) || is.numeric(x) || is.logical(x)) || is.object(x))
    return(sprintf("an object of class %s", class(x)[1L]))
  if (!length(x)) return(sprintf("%s(0)", class(x)[1L]))
  # Each number as it prints alone: format() of the vector would pad them to
  # a common number of decimals.
  v <- if (is.character(x)) ifelse(is.na(x), "NA", sprintf("\"%s\"", x))
       else vapply(x, function(one) format(one), "")
  txt <- paste(utils::head(v, 6L), collapse = ", ")
  if (length(v) > 6L) txt <- paste0(txt, ", ...")
  if (length(v) > 1L) sprintf("c(%s)", txt) else txt
}


#' Items joined as a sentence lists them
#'
#' @param x Character vector, at least one item.
#' @return \code{a}; \code{a and b}; \code{a, b and c}.
#' @keywords internal
#' @noRd
.and_text <- function(x) {
  k <- length(x)
  if (k <= 1L) return(paste(x))
  paste(paste(x[-k], collapse = ", "), "and", x[k])
}


#' Whether a value is the string "auto"
#'
#' Read by its value alone.  A string taken from a named vector or list
#' (\code{cfg["sac"]}) carries a name, and one built by another package can
#' carry a class; \code{identical(x, "auto")} is \code{FALSE} of both.
#'
#' @param x The value the caller passed.
#' @return \code{TRUE} or \code{FALSE}.
#' @keywords internal
#' @noRd
.is_auto <- function(x)
  is.character(x) && length(x) == 1L && identical(as.character(unclass(x)), "auto")


#' Rows per cell, in the words of a message
#'
#' @param x The number of rows per cell.
#' @return \code{"3 rows per cell"}, \code{"1.03 rows per cell"}; \code{"1 row
#'   per cell"} when the number prints as 1.
#' @keywords internal
#' @noRd
.per_cell_text <- function(x) {
  v <- sprintf("%.3g", x)
  paste(v, if (identical(v, "1")) "row" else "rows", "per cell")
}


#' The finest level of a ladder of sizes
#'
#' The count the finest size gives, held to the most centres k-means can
#' place.  It needs the rows and the sizes alone, so the guards on the cost
#' of a ladder of sizes read it before anything is estimated;
#' \code{.size_ladder()} returns the same number as its \code{ceiling}.
#'
#' @param n The rows a size divides.
#' @param sizes Rows per cell.
#' @param cap The most centres k-means can place.
#' @return One integer, 2 or more.
#' @keywords internal
#' @noRd
.finest_level <- function(n, sizes, cap)
  as.integer(min(.count_at(n, min(sizes)), max(2L, as.integer(cap))))


#' The levels of a ladder of sizes, its verdict, and what is said of it
#'
#' The level of a size \code{s} is the count \code{.count_at()} gives for it,
#' held to the most centres k-means can place.  The finest size takes the
#' part \code{min_cell_n} has on the automatic ladder: it sets the ceiling,
#' \code{supported} and the warning.  One event raises one warning.
#' \describe{
#'   \item{Fewer than two cells of the finest size, or one distinct
#'     location}{\code{spatialkit_not_supported}; the ladder is the one level
#'     of 2 cells and \code{supported} is \code{FALSE}.  Nothing else is
#'     raised.}
#'   \item{A coarser size gives fewer than two cells}{Given as a number: its
#'     level is 2, with \code{spatialkit_size_not_reached}.  From
#'     \code{"auto"}: the size is dropped and nothing is said, since the
#'     caller did not ask for it.}
#'   \item{A size asks for more cells than can be placed}{Its level is the
#'     most that can, with \code{spatialkit_size_not_reached}.}
#'   \item{The range floor is above the finest level}{
#'     \code{spatialkit_not_supported}, as on the automatic ladder.}
#' }
#'
#' @param sizes Rows per cell, largest first, without repeats.
#' @param sizes_auto Whether they came from \code{obs_per_cell = "auto"}.
#' @param n The rows a size divides: those with a response, or every point.
#' @param n_points,n_distinct The layer's points and its distinct locations.
#' @param cap The most centres k-means can place on them.
#' @param floor_raw The range floor, a double (2 without a range).
#' @param range_eff The range it came from, for the message.
#' @param support_what,rows_carry What \code{n} counts and what a counted row
#'   carries, in the words of the messages.
#' @return A list: \code{levels} (sorted, distinct), \code{size_levels} (the
#'   level of each size, named by size; \code{NA} for a size dropped),
#'   \code{ceiling} (the finest level), \code{ceiling_from}, \code{supported}
#'   and \code{too_few}.
#' @keywords internal
#' @noRd
.size_ladder <- function(sizes, sizes_auto, n, n_points, n_distinct, cap, floor_raw,
                         range_eff, support_what, rows_carry) {
  cap    <- as.integer(cap)
  wanted <- round(n / sizes)                 # the cells each size asks for
  s_min  <- sizes[length(sizes)]
  want   <- wanted[length(sizes)]            # of the finest size
  too_few   <- min(want, cap) < 2
  supported <- !too_few && floor_raw <= max(2, min(want, cap))
  ceiling_from <- if (cap < 2L) "fewer than two distinct locations"
                  else if (want < 2) "fewer than two cells of the size"
                  else if (want > cap) "distinct locations"
                  else "obs_per_cell"
  size_levels <- pmin(.count_at(n, sizes), max(2L, cap))
  if (!too_few && isTRUE(sizes_auto)) size_levels[wanted < 2] <- NA_integer_
  names(size_levels) <- as.character(sizes)
  levels <- sort(unique(size_levels[!is.na(size_levels)]))
  finest <- max(levels)
  size_text <- function(s) format(s)

  if (too_few && cap < 2L) {
    .warn_classed("spatialkit_not_supported",
                  paste0("resolution_profile(): the %d points sit on %d distinct location%s, ",
                         "so there is no ladder of cell counts to choose from. The profile ",
                         "scores 2 cells so that level is still visible, and `supported` is ",
                         "FALSE. A tessellation needs points at two locations or more."),
                  n_points, n_distinct, if (n_distinct == 1L) "" else "s",
                  data = list(reason = ceiling_from, n = n, min_cell_n = NA_integer_,
                              obs_per_cell = s_min))
  } else if (too_few) {
    # The largest size that still gives two cells: round(n / s) is 2 from
    # n / s = 1.5 up.  Written to two decimals, rounded down.
    most <- floor(n / 1.5 * 100) / 100
    .warn_classed("spatialkit_not_supported",
                  paste0("resolution_profile(): %d %s at obs_per_cell = %s make%s %d cell%s, and ",
                         "it takes two, so there is no cell count of that size. The profile ",
                         "scores 2 cells (%s) so that level is still visible, ",
                         "and `supported` is FALSE. %s"),
                  n,
                  if (n == 1L) sub("^rows", "row", sub("^points$", "point", support_what))
                  else support_what,
                  size_text(s_min), if (n == 1L) "s" else "", as.integer(want),
                  if (want == 1) "" else "s", .per_cell_text(n / 2),
                  if (most < 1)
                    sprintf("Two cells take %s on at least two rows, whatever `obs_per_cell` is.",
                            rows_carry)
                  else sprintf("An `obs_per_cell` of at most %s gives two cells here.",
                               format(most)),
                  data = list(reason = ceiling_from, n = n, min_cell_n = NA_integer_,
                              obs_per_cell = s_min))
  } else {
    for (j in seq_along(sizes)) {
      if (wanted[j] < 2 && !isTRUE(sizes_auto))
        .warn_classed("spatialkit_size_not_reached",
                      paste0("resolution_profile(): obs_per_cell = %s on %d rows is fewer ",
                             "than two cells; two cells are used (%s)."),
                      size_text(sizes[j]), n, .per_cell_text(n / 2),
                      data = list(obs_per_cell = sizes[j], wanted = as.integer(wanted[j]),
                                  used = 2L, reason = "fewer than two cells"))
      else if (wanted[j] > cap)
        .warn_classed("spatialkit_size_not_reached",
                      paste0("resolution_profile(): obs_per_cell = %s asks for %d cells, more ",
                             "than k-means can place on %s; %d cells are used (%s)."),
                      size_text(sizes[j]), as.integer(wanted[j]),
                      # With no location repeated the bound is one short of
                      # the points, which stats::kmeans() needs.
                      if (n_distinct < n_points)
                        sprintf("the layer's %d distinct locations", n_distinct)
                      else sprintf("%d points", n_points),
                      cap, .per_cell_text(n / cap),
                      data = list(obs_per_cell = sizes[j], wanted = as.integer(wanted[j]),
                                  used = cap, reason = "more cells than can be placed"))
    }
    if (!supported)
      .warn_classed("spatialkit_not_supported",
                    paste0("resolution_profile(): cells no wider than the autocorrelation ",
                           "range (%.0f) would need at least %.0f of them, but %d %s ",
                           "at obs_per_cell = %s support at most %d. The data cannot ",
                           "support a tessellation that respects their own ",
                           "correlation structure; the profile scores %s so ",
                           "the cost of each level is still visible."),
                    range_eff, floor_raw, n, support_what, size_text(s_min), finest,
                    if (length(levels) == 1L) sprintf("the %d cells that size gives", levels)
                    else sprintf("the %d to %d cells the sizes give", min(levels), max(levels)),
                    data = list(reason = "floor above ceiling", n = n,
                                min_cell_n = NA_integer_, obs_per_cell = s_min))
  }
  list(levels = levels, size_levels = size_levels, ceiling = finest,
       ceiling_from = ceiling_from, supported = supported, too_few = too_few)
}


#' The most point-to-centre comparisons a fit at a fixed size may make
#'
#' @return \code{getOption("spatialkit.fit_work_max")}, or
#'   \code{.FIT_WORK_MAX} when the option is not set.  A value that is not
#'   one number is an error and not a silent return to the default.
#' @keywords internal
#' @noRd
.fit_work_max <- function() {
  lim <- getOption("spatialkit.fit_work_max", .FIT_WORK_MAX)
  if (!is.numeric(lim) || length(lim) != 1L || is.na(lim) || lim < 0)
    stop("options(spatialkit.fit_work_max = ) must be a single number of comparisons, ",
         "0 or more; got ", .got_text(lim), ".", call. = FALSE)
  as.numeric(lim)
}


#' Stop before partition fits that compare too many points with centres
#'
#' On a ladder of sizes the count of cells grows with the rows, so a step
#' that compares every point with every centre is quadratic in the rows.
#' Which fits have such a step:
#' \describe{
#'   \item{kind \code{"kmeans"}, either starting rule}{\code{stats::kmeans()}
#'     compares points with centres in every pass, and a k-means++ start
#'     makes one pass over the points per centre.}
#'   \item{kind \code{"density"} from k-means++ starts, with \pkg{FNN}}{The
#'     passes use a kd-tree; the start makes one draw per centre over all the
#'     locations.}
#'   \item{kind \code{"density"} without \pkg{FNN}}{Every pass compares every
#'     location with every centre, and a start takes up to \code{iter_max}
#'     passes: the count is that worst case.}
#'   \item{kind \code{"density"} from random starts, with \pkg{FNN}}{None,
#'     and no work is counted.}
#' }
#' The work is the points times the cells of the finest level times the
#' starts, times \code{iter_max} in the third case.
#'
#' @param n_points The rows the centres are fitted to.
#' @param cells The finest level.
#' @param nstart Number of starts, as resolved.
#' @param part The partition record: \code{kind} and \code{init}.
#' @param iter_max The iteration limit of a start.
#' @param use_fnn Whether \pkg{FNN} is used.  Evaluated for the kind
#'   \code{"density"}, and for the kind \code{"kmeans"} only when the call
#'   stops.
#' @return \code{NULL}, invisibly; an error of class
#'   \code{spatialkit_work_limit} when the work is above
#'   \code{.fit_work_max()}, with \code{work}, \code{limit}, \code{what}
#'   (\code{"partition fit"}), \code{kind}, \code{init}, \code{nstart} and
#'   \code{fnn}.
#' @keywords internal
#' @noRd
.fit_guard <- function(n_points, cells, nstart, part, iter_max, use_fnn) {
  density <- identical(part$kind, "density")
  # `use_fnn` is read here for the density kind alone, whose count depends on
  # it.  For the kind "kmeans" it is read only when the call stops, for the
  # condition's `fnn`, so a k-means profile that passes has not asked.
  slow <- density && !isTRUE(use_fnn)
  if (density && !slow && identical(part$init, "random")) return(invisible(NULL))
  passes <- if (slow) as.numeric(iter_max) else 1
  work   <- as.numeric(n_points) * as.numeric(cells) * as.numeric(nstart) * passes
  limit  <- .fit_work_max()
  if (!(work > limit)) return(invisible(NULL))
  fnn <- isTRUE(use_fnn)
  one <- identical(as.integer(nstart), 1L)
  with <- if (density && !fnn)
    sprintf("%d start%s of up to %d pass%s each", as.integer(nstart), if (one) "" else "s",
            as.integer(iter_max), if (identical(as.integer(iter_max), 1L)) "" else "es")
  else sprintf("%d %s%s start%s", as.integer(nstart), if (density) "density-weighted " else "",
               if (identical(part$init, "random")) "random" else "k-means++",
               if (one) "" else "s")
  msg <- sprintf(paste0("resolution_profile(): fitting %d cells to %d rows with %s is %.2g ",
                        "point-to-centre comparisons, above the limit of %.2g ",
                        "(options(spatialkit.fit_work_max = )). Ways on: fewer starts ",
                        "(`nstart`), which keep the best of fewer; a coarser `obs_per_cell`; ",
                        "partition = \"density\" with init = \"random\" and the FNN package ",
                        "installed, which has no step that compares every point with every ",
                        "centre and leaves a cell without a point more often; or a higher ",
                        "limit."),
                 as.integer(cells), as.integer(n_points), with, work, limit)
  stop(errorCondition(msg, class = "spatialkit_work_limit", call = NULL,
                      work = work, limit = limit, what = "partition fit",
                      kind = part$kind, init = part$init, nstart = as.integer(nstart),
                      fnn = fnn))
}


#' The estimate's numbers at one level of a profile
#'
#' Each scored row of the layer goes to its nearest centre, which is the
#' Voronoi cell the seeding functions hand on and \code{build_tessellation()}
#' draws.  Centres that hold no scored row are dropped and the rest numbered
#' in the order of their centres.  Then the engine gives the numbers
#' (\code{.pool_fit()}).  The cluster vector of the fit is not read: on a
#' subsample it covers the subsample only, and a start cut off at its
#' iteration limit can leave a row outside its nearest centre's cell.
#'
#' @param xy Two-column matrix of the scored rows' coordinates.
#' @param y Their values.
#' @param centres Two-column matrix of the level's centres.
#' @param estimate \code{"mean"} or \code{"pooled"}.
#' @param use_fnn Whether \pkg{FNN} is used.  The numbers are the same
#'   without it.
#' @param cache The environment the levels of one profile share for the
#'   nearest other location of a location.
#' @return A list: \code{cells} (the cells that hold a scored row),
#'   \code{scored_per_cell} (rows over those cells), \code{loo} and
#'   \code{loo_se} (the leave-one-out error of the estimate and its standard
#'   error over rows), \code{share_alone} (the share of rows that are the
#'   only scored row of their cell), \code{own_weight} (the median own weight
#'   of a pooled value; \code{NA} under plain means and for a constant
#'   column), \code{sigma2}, \code{tau2} (\code{NA} under plain means) and
#'   \code{f_p}.
#' @keywords internal
#' @noRd
.estimate_at_level <- function(xy, y, centres, estimate, use_fnn, cache) {
  N  <- length(y)
  g  <- .nearest_centre(xy[, 1L], xy[, 2L], centres[, 1L], centres[, 2L], use_fnn)
  u  <- sort(unique(g))
  gi <- match(g, u)
  K  <- length(u)
  fit <- .pool_fit(y, gi, K, xy, want = estimate, use_fnn = use_fnn, cache = cache,
                   caller = "resolution_profile")
  pooled <- identical(estimate, "pooled")
  list(cells = K, scored_per_cell = N / K,
       loo    = if (pooled) fit$loo_pooled else fit$loo_mean,
       loo_se = if (pooled) fit$loo_pooled_se else fit$loo_mean_se,
       share_alone = sum(fit$n == 1L) / N,
       own_weight  = if (pooled && !isTRUE(fit$constant)) stats::median(fit$w) else NA_real_,
       sigma2 = fit$sigma2, tau2 = if (pooled) fit$tau2 else NA_real_, f_p = fit$f_p)
}


#' Score every candidate number of cells on several criteria at once
#'
#' A tessellation's resolution is the number of cells the points are cut
#' into, and no single number settles it: the within-cluster sum of squares
#' of the coordinates has an elbow, a piecewise-constant approximation of the
#' response has a Mallows \eqn{C_p}, the residual autocorrelation of the cell
#' means has a \eqn{z}, the cell means have a reliability, and a cell value
#' has a leave-one-out error.  On a ladder of counts this function computes
#' all of them at every level the data bound, and returns the table, so that
#' the level chosen, by
#' \code{\link{select_resolution}()} or by eye, can be defended with the
#' whole profile in place of one criterion's argmin.  Mallows (1973) presents
#' \eqn{C_p} itself as a display of the bias-variance trade-off across
#' candidates; he does not offer it as a rule that picks one.  This is that
#' display, with the other criteria alongside.
#'
#' Called with a numeric response and no predictors and none of
#' \code{levels}, \code{n_levels}, \code{min_cell_n}, \code{sample_n} or
#' \code{range_floor}, the profile is one level: 4.5 rows per cell
#' (\code{obs_per_cell = "auto"}), scored for pooled values
#' (\code{estimate = "auto"}) by the leave-one-out error \code{loo}, with
#' no variogram fitted (\code{sac = "auto"}), in regions of the density
#' kind (\code{partition = "density"}) kept from 3 k-means++ starts
#' (\code{nstart = NULL}).  \code{cp} and \code{reliability} need the
#' variogram and \code{elbow} needs a curve, so the three are \code{NA}
#' there; \code{moran_z} is computed.  \code{\link{select_resolution}()}
#' reads the count off it by the size.  Any of those five ladder arguments
#' gives the
#' ladder of counts described below, scored for pooled values as well.
#' \code{estimate = NULL, obs_per_cell = NULL, sac = NULL, partition =
#' "kmeans"} gives the profile of earlier versions: the ladder, its columns
#' and its k-means regions from 25 starts.
#'
#' @section The partition:
#' The regions at every level are the Voronoi cells of a set of centres: a
#' point belongs to the centre nearest to it.  \code{partition} says what the
#' centres minimise, and every level of the ladder is built with the same
#' kind.
#' \describe{
#'   \item{\code{"kmeans"}}{The sum, over the points, of the
#'     squared distance to the nearest centre.  The centres follow the square
#'     root of the point density.  The cells are compact, smaller where the
#'     points are dense and larger where they are sparse, and a cell in a
#'     dense place holds several times the points of a cell in a sparse
#'     place.}
#'   \item{\code{"density"} (default)}{The same sum with each point weighted
#'     by its local density to the power \code{density_power}.  The local
#'     density at a point is the number of points in the disc that reaches
#'     its \code{density_k}-th nearest neighbour, over the area of that disc.
#'     With \code{density_power = 1} the centres follow the point density
#'     itself.  The cells are compact: many small ones where the points are
#'     dense, a few large ones where they are sparse.  They hold about the
#'     same number of points each once there are more cells than dense
#'     places (clusters of points).  With fewer cells than dense places a
#'     cell takes in whole clusters, large and small, and the counts are
#'     about as uneven as under \code{"kmeans"}.  \code{density_power = 0}
#'     weights every point alike, which is the objective of
#'     \code{"kmeans"}; a power above 1 goes past equal counts and gives the
#'     cells of dense places fewer points than those of sparse places.}
#' }
#' Two simulations show both ends.  The measure is the coefficient of
#' variation of the cell counts (their standard deviation over their mean,
#' 0 when every cell holds the same number), as medians over 20 layers.  On
#' 1,200 points, 900 of them in four clusters over a sparse background, cut
#' into 60 cells, it was 1.2 under \code{"kmeans"}, where the fullest cell
#' held 42 times the points of the smallest, and 0.29 under
#' \code{"density"}, where it held 4 times.  On 2,000 points with 28 dense
#' places (half of the points in 3 large clusters, 30 percent in 25 small
#' ones, the rest spread out) it was 0.90 under \code{"kmeans"} and 0.97
#' under \code{"density"} at 10 cells and 1.08 and 1.01 at 20 cells; at 60
#' cells it was 0.95 and 0.31, and at 130 cells 1.05 and 0.43.  A ladder
#' starts at 2 cells, so its low levels are usually below the number of
#' dense places, and the counts in \code{cell_n_min} and
#' \code{cell_n_median} show where the two kinds part.
#'
#' The density kind is the default.  The default size, 4.5 rows per cell,
#' puts the count above the number of dense places on a layer whose
#' clusters hold more than a few observations each, so the default kind and
#' the default count go together.  A coarse ladder the caller sets (with
#' \code{levels}, large sizes, or the low levels of the automatic ladder)
#' may not, and at those levels the density kind does not even the counts.
#' \code{partition = "kmeans"} gives the k-means regions, with 25
#' k-means++ starts, as earlier versions built them.
#'
#' \code{init} says how each of the \code{nstart} starts is drawn before it
#' is refined; the best start is kept.  \code{"kmeans++"} draws the centres
#' one at a time, each a point taken in proportion to its squared distance
#' from the centres already drawn (times its weight, under
#' \code{"density"}), which spreads the start over the points.
#' \code{"random"} draws \code{k} distinct locations at random, as
#' \code{\link[stats]{kmeans}()} does.
#'
#' A start of the kind \code{"kmeans"} is refined by
#' \code{\link[stats]{kmeans}()} (Hartigan and Wong), for at most 50
#' iterations from a k-means++ start and 100 from a random one.  A start of
#' the kind \code{"density"} is refined by Lloyd's iteration, for at most 100
#' steps: each point goes to its nearest centre, each centre moves to the
#' weighted mean of its points, and a centre left without a point stays
#' where it is.  \code{iter_max} sets another limit.  The start with the
#' lowest sum is kept, and a kept start that was cut off at the limit is
#' reported (see "Conditions"): its centres are those of the last step.
#' Large layers can need more than 100 steps: on a simulated layer of 30,000
#' points no start of six at 100 and 300 cells took more than 68, and on one
#' of 60,000 points three of six took 102 to 152.
#'
#' The density kind depends on the points alone.  It works on the distinct
#' locations in coordinate order, so the order of the rows does not matter.
#' Its two tolerances are fractions of the extent of the points (the diagonal
#' of their bounding box), so the unit of the coordinates does not matter:
#' the neighbour distance is not taken below a millionth of the extent, and
#' the iteration stops when no centre moves by \eqn{10^{-12}} of it.  Points
#' at one location count as that many points at it and share a cell; the
#' neighbour distance is measured to the \code{density_k}-th nearest
#' \emph{other} location, so it is not zero where locations repeat.  The
#' neighbour searches use \pkg{FNN} when it is installed.  The result is the
#' same without it; the search then takes time that grows with the square of
#' the number of points, about thirty times as long for 10,000 points and
#' 2,000 centres.
#'
#' The columns keep their meaning under \code{"density"}, because each is
#' read off the cells and not off the way they were drawn.  \code{wss} is the
#' within-cell sum of squares of the coordinates with every point counted
#' once, taken about each cell's mean, so \code{wss} and \code{elbow}
#' compare between kinds; the weighted sum the centres minimise is kept in
#' \code{attr(, "partition")$objective}, and \code{wss_spread} is its spread
#' over the starts.  \code{cp}, \code{moran_z} and the test against one cell
#' need cells drawn without the response, which both kinds give.
#' \code{reliability} is computed for square cells of the level's average
#' area and count and never reads the cells, so it is the same number under
#' either kind.  The \code{elbow} is found more readily under
#' \code{"density"} where the points cluster: once each cluster has its
#' centres, further centres go into the clusters, the sparse ground between
#' them keeps its sum of squares, and the curve flattens.  On simulated
#' layers of 800 points in 3 to 8 clusters over a sparse background an elbow
#' was read on 24 of 30 under \code{"density"} and on 10 under
#' \code{"kmeans"}; on 100 uniform layers, on none under either.
#'
#' A profile of the density kind, or of either kind from random starts,
#' carries \code{attr(, "partition")}, and
#' \code{\link{select_resolution}()}, the seeding functions and
#' \code{\link{build_tessellation}()} hand that record on with the regions
#' (see "Value").  A result with no such record was built by plain k-means:
#' from k-means++ starts when it comes from a profile, and from random
#' starts when a seeding function fitted it.
#'
#' @section The ladder and its bounds:
#' This is the ladder of a call that sets \code{levels}, \code{n_levels},
#' \code{min_cell_n}, \code{sample_n} or \code{range_floor}, or passes
#' \code{obs_per_cell = NULL}, and of a call with no estimate: one without
#' a response, with \code{predictor_vars}, or with a logical or wholly
#' missing response.  Cells are clusters of the (projected)
#' coordinates of the kind \code{partition} names, fitted at each level as
#' the best of \code{nstart} k-means++ starts (see
#' \code{\link{determine_optimal_levels}} for why) unless \code{init} says
#' otherwise (see "The partition").  Levels are spaced
#' logarithmically, because cell diameter scales as \eqn{L^{-1/2}}: a unit
#' step wastes fits at large \eqn{L} and starves resolution at small.  The
#' ladder runs from 2 to a ceiling the data impose.  The ceiling is
#' \code{floor(n / min_cell_n)}: cells with fewer than \code{min_cell_n}
#' rows on average have too little support.  With a \code{response_var},
#' \eqn{n} is the number of rows of the layer (not of the subsample) that
#' are scored: those with a response and, with \code{predictor_vars}, every
#' predictor, since the cell means the criteria read are built from those
#' rows alone.  \code{ceiling_from} is then \code{"min_cell_n (rows with a
#' response)"} when some rows are not scored (the label is the same when it
#' is a predictor they lack), and a warning of class
#' \code{spatialkit_response_ceiling} gives both counts and names the bound
#' in force.  Without a
#' \code{response_var}, or when no row has a response, \eqn{n} is every
#' point.  The ceiling is also held to the number of distinct locations, one
#' cell on each, and to one short of the number of points, which
#' \code{stats::kmeans()} needs.  The count \eqn{n / L} is the rows per cell
#' only when the scored rows are spread over the whole extent.
#' Where they sit in one part of it, the cells elsewhere hold none and the
#' occupied cells hold more than \eqn{n / L}, so the ceiling is lower than
#' those cells need: \code{cells_resp} and \code{cell_resp_mean} give the
#' cells that hold a scored row and the scored rows per such cell at each
#' level (of the subsample, when the cells are fitted to one), and
#' profiling the scored rows alone bounds the ladder by the area they
#' cover.
#'
#' When fewer than two cells meet \code{min_cell_n}
#' (\code{floor(n / min_cell_n) < 2}), or the layer has one distinct
#' location, there is no ladder: the profile scores two cells so that level
#' is still visible (or the \code{levels} asked for, as given),
#' \code{supported} is \code{FALSE}, \code{ceiling_from}
#' is \code{"fewer than two supported cells"} (or \code{"fewer than two
#' distinct locations"}) and a warning of class
#' \code{spatialkit_not_supported} names both numbers.
#'
#' \code{min_cell_n} sets the ceiling of this ladder.  It is a mean number
#' of rows per cell and not a minimum: a cell of a level at the ceiling can
#' hold fewer (\code{cell_n_min}).  It is floored to a whole number, so 4.5
#' is read as 4.  With \code{levels} it does not choose the levels scored,
#' and it still sets \code{bounds$ceiling}, \code{supported} and the
#' \code{spatialkit_not_supported} warning.
#'
#' The range floor is \code{ceiling(area / range^2)} when an autocorrelation
#' range is available: cells wider than the range average over more than one
#' patch of the field.  It is reported in the bounds and the print, and by
#' default it is not applied: the ladder starts at 2.  When the layer
#' supports the floor and it lies between the ladder's ends, the floor is
#' one of the ladder's levels all the same (so the ladder can have one more
#' level than \code{n_levels}): \code{select_resolution(criterion =
#' "reliability")} reads from the floor up and returns the floor itself, not
#' the next level the spacing happened to put above it.  The range is an
#' estimate, and a response with no spatial structure can still be given
#' one: over 100 simulated white-noise layers of 900 points the variogram
#' returned a range on 29, and started at the floor, 11 of those ladders
#' began at 11 to 98 cells of a ceiling of 100 and \eqn{C_p} took the floor
#' or a level beside it, where the best count was 2 to 8.
#' \code{range_floor = TRUE} starts the ladder at the floor, as the default
#' did up to version 2.0.0.  When
#' the floor exceeds the ceiling the data cannot support a tessellation
#' that respects their own correlation structure; that is reported as a
#' finding (a warning of class \code{spatialkit_not_supported}, and
#' \code{attr(x, "bounds")$supported} is \code{FALSE}) and the ladder runs
#' from 2 to the ceiling, so the profile still shows what each level costs.
#' A floor equal to the ceiling would leave one level and nothing to compare
#' it with, so that ladder starts at 2 as well.  On a layer larger than
#' \code{sample_n} the ceiling is also held to half the subsample's scored
#' rows (two per cell, the least a fitted cell can be scored on; half the
#' subsample without a response); \code{ceiling_from} is then
#' \code{"sample_n"}, a floor above that is logged with a request to raise
#' \code{sample_n}, and it does not make \code{supported} \code{FALSE}.
#'
#' The floor moves with the range estimate, and the ceiling with \eqn{n}.
#' To compare profiles of similar data (the folds of a cross-validation,
#' say), pass the same \code{levels} to each.  \code{reliability} is the one
#' criterion that still reads the floor: \code{\link{select_resolution}()}
#' holds it to the levels at or above the floor by default, so its pick can
#' differ between two such profiles by a factor of 10 or more when the floor
#' is supported on one and not on the other; pass \code{range_floor = FALSE}
#' there to read the whole ladder.
#'
#' A ladder that \code{obs_per_cell} sets is a ladder of sizes (see
#' "Observations per cell").  It has one level per size \eqn{s},
#' \code{max(2, round(n / s))} cells, held to the most centres k-means can
#' place: the distinct locations, and one short of the points.  Its centres
#' are fitted to every row.  No subsample is drawn, and \code{n_levels},
#' \code{min_cell_n}, \code{sample_n} and \code{range_floor} are not read.
#' The finest size takes the part \code{min_cell_n} has on the ladder above:
#' \code{bounds$ceiling} is the finest level, and \code{supported} is
#' \code{FALSE} when the finest size gives fewer than two cells or the range
#' floor is above the finest level.  The range floor is reported and not
#' applied.
#'
#' A ladder of sizes under \code{sac = "auto"} fits no variogram, so
#' \code{cp} and \code{reliability} are \code{NA} at every level there, and
#' \code{attr(, "cp_noise")} says why.  \code{sac = NULL}, or a \code{sac} of
#' the caller's own, gives them.  \code{sac = "auto"} and
#' \code{obs_per_cell = "auto"} are the defaults, so the default profile
#' has neither column.
#'
#' @section Observations per cell:
#' \code{obs_per_cell} sets the ladder by a size: the mean number of rows a
#' cell is to hold.  The level of a size \eqn{s} is
#' \code{max(2, round(n / s))} cells, with \eqn{n} the rows the size
#' divides (see the argument).  A size says how many observations a cell
#' value is to stand on, and that is the caller's choice: no criterion makes
#' it.  Under a size neither the cells, which are fitted to the coordinates,
#' nor their number, which is a count of rows, depend on the values of the
#' response.  Under \code{select_resolution(criterion = "loo")} the count
#' does.
#'
#' For pooled values the leave-one-out error tends to keep falling as the
#' cells shrink, so its minimum tends to be the finest size offered, and the
#' end of the ladder would then choose the count.  In a simulation
#' of 200 layers, each 1,500 clustered observations of a field with
#' covariance \code{exp(-d / 100)} plus noise on a square of side 1,000
#' (\code{dev/sim-pooled-values.R} in the package's repository;
#' \code{\link{summarize_by_cell}} describes it), profiled at the sizes 18,
#' 12, 9, 6, 4.5, 3 and 2 with \code{partition = "density"}, the lowest
#' \code{loo} of the pooled value was at the size 2 on 176 layers, at 3 on
#' 22 and at 4.5 on 2.  That of the plain mean was at 2 on none: it was at 6
#' on 86 layers, at 4.5 on 51, at 9 on 49, at 12 on 12, and at 18 and at 3
#' on one each.  \code{select_resolution(criterion = "size")} therefore sets
#' the count of pooled values from a size, 4.5 unless another is given, and
#' \code{criterion = "loo"} reads the lowest error for plain means.  These
#' figures are of those seven sizes and say nothing of another ladder.
#'
#' Two columns say what a size leaves a cell value to stand on.
#' \code{share_alone} is the share of the scored rows that are the only
#' scored row of their cell: the plain mean of such a cell is one
#' observation, and its pooled value has the smallest own weight a cell can
#' have.  \code{own_weight} is the median weight a pooled value gives the
#' cell's own mean: at 1 the pooled values are the plain means, and at 0
#' each is the mean of the cell's neighbours.
#'
#' \code{obs_per_cell = "auto"}, the default, is resolved against the other
#' ladder arguments, in this order.  Whether the caller supplied an argument
#' is what counts, so a value a wrapper passes on counts as supplied.
#' \enumerate{
#'   \item With \code{levels}, \code{"auto"} is \code{NULL}: the levels are
#'     scored as given.  A numeric \code{obs_per_cell} beside \code{levels}
#'     is an error, since each sets the ladder.
#'   \item With \code{n_levels}, \code{min_cell_n}, \code{sample_n} or
#'     \code{range_floor} supplied, \code{"auto"} is \code{NULL}: those
#'     arguments describe the ladder of "The ladder and its bounds", and the
#'     caller gets it, with the estimate's columns when there is an
#'     estimate.
#'   \item Otherwise \code{"auto"} is 4.5 when the estimate is
#'     \code{"pooled"}, the seven sizes 18, 12, 9, 6, 4.5, 3 and 2 when it
#'     is \code{"mean"}, and \code{NULL} without an estimate.
#'   \item A numeric \code{obs_per_cell} beside \code{n_levels},
#'     \code{min_cell_n}, \code{sample_n} or \code{range_floor = TRUE} sets
#'     the ladder, with a warning of class \code{spatialkit_ladder_unused}
#'     that names the arguments that are not read.
#' }
#'
#' On a ladder of sizes the centres are fitted to every row and the count of
#' cells grows with the rows, so a step that compares every point with every
#' centre grows with the square of the rows.  \code{partition = "kmeans"}
#' has such a step in every pass and in a k-means++ start.  \code{partition
#' = "density"} has it in a k-means++ start and, without \pkg{FNN}, in every
#' pass; with \pkg{FNN} and \code{init = "random"} it has none.  The profile
#' stops before any fit, with an error of class
#' \code{spatialkit_work_limit}, when the points times the cells of the
#' finest level times \code{nstart} (times \code{iter_max} for the density
#' kind without \pkg{FNN}) is above
#' \code{getOption("spatialkit.fit_work_max", 1e10)}.  The message names the
#' ways on: fewer starts, which keep the best of fewer; a coarser size;
#' \code{partition = "density"} with \code{init = "random"} and \pkg{FNN}
#' installed, which leaves a cell without a point more often (see
#' \code{spatialkit_empty_cells} under "Conditions"); or a higher limit.
#' The limit counts the finest level, and a ladder of several sizes fits
#' every one of its levels.  At the default size the work grows with the
#' square of the rows, and the limit is set where the default call's fit
#' takes about a minute: a call just below it runs that long, and prints
#' nothing while it runs.
#'
#' Without \pkg{FNN} the density kind also finds its weights by comparing
#' every pair of distinct locations.  A profile on a ladder of sizes stops,
#' with the same class, when the square of their number is above
#' \code{getOption("spatialkit.exact_pairs_max", 2e9)}, which is about
#' 44,700 locations.  The result of every search is the same with and
#' without the package.  \code{moran_i} and \code{moran_z} need \pkg{FNN}
#' and \pkg{Matrix} at a level where more than 5,000 cells hold a scored
#' row.  Without them a ladder of sizes has \code{NA} in both columns
#' there, with a warning of class \code{spatialkit_moran_skipped}, and so
#' does a ladder set any other way (see "Conditions").
#' At a level of up to 5,000 such cells, when either package is missing, the
#' neighbour weights of the cell means are held as a dense matrix of cells
#' by cells.
#'
#' @section The criteria, and what each returns in practice:
#' The four columns are not four interchangeable count rules.
#' \code{select_resolution()} reads \code{cp} by default on a profile
#' without an estimate; on one with an estimate it reads the size (pooled
#' values) or \code{loo} (plain means).  \code{reliability}
#' needs a variogram range and is \code{NA} without one; \code{elbow} is
#' \code{NA} unless the points cluster; and \code{moran_z} is a diagnostic
#' and not a count rule.  With an \code{estimate} the profile has a fifth
#' column, \code{loo}, which is the last entry.  Each entry below says what
#' the column measures and what it returns in practice.
#' \describe{
#'   \item{\code{elbow}}{How far the WSS curve sags below a power law: on
#'     log-log axes, \eqn{\log} WSS below the straight line from \eqn{k = 1}
#'     (the total sum of squares) to the last level, in natural-log units
#'     (larger is better).  Points with no cluster structure have a WSS close
#'     to \eqn{c/k}, which is straight on those axes, so the column is
#'     \code{NA} at every level unless the largest sag reaches
#'     \eqn{\log 1.25}, and the print says there is no elbow (the rule and its
#'     calibration are in \code{\link{determine_optimal_levels}}).  The
#'     classical chord on linear axes found a "knee" on such a layer anyway,
#'     at about \eqn{\sqrt{L_{first} L_{last}}}, where the ladder's ends put
#'     it.  A level whose WSS is 0 (to within \eqn{10^{-12}} of the total),
#'     one cell on every distinct location, which the ladder reaches when
#'     locations repeat, is left out of the line; when the other levels have
#'     no elbow, that fall to zero is the elbow, and its sag is measured with
#'     the WSS floored at \eqn{10^{-12}} of the total, which puts it far
#'     above any other level's.  Geometry only; it knows nothing of the
#'     response.  \strong{In practice it is \code{NA} unless the points
#'     cluster}, and where they do it names the number of clusters, which
#'     is usually the first level or two of the ladder.}
#'   \item{\code{cp}}{Mallows' \eqn{C_p} of the piecewise-constant
#'     approximation of the response (or of its OLS residuals on
#'     \code{predictor_vars}) by cell means: \eqn{RSS(L)/n + 2 \tau^2 L / n},
#'     with \eqn{\tau^2} the nugget of the fitted variogram (lower is better).
#'     It estimates the error of re-measuring the sampled sites: each scored
#'     row observed again at its own location, with fresh noise of variance
#'     \eqn{\tau^2}, and predicted by the mean of its cell.  For an
#'     observation at a new location the penalty needs the spread of the
#'     field inside a cell as well, which the nugget leaves out.  This
#'     \eqn{\tau^2} is the nugget, the noise of one measurement; it is not
#'     the \code{tau2} of a pooled value, which is the mean squared distance
#'     between a cell's level and its target.  When the
#'     cells are fitted to a subsample of \eqn{m} of the
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
#'     only when the nugget is a large share of the sill.  With a nugget of 0
#'     (under \eqn{10^{-4}} of the sill; usually a fit clipped at its lower
#'     bound) the penalty is 0 and \eqn{C_p} descends to the ceiling whatever
#'     the field; the profile warns.  When the variogram gives no nugget (none
#'     was fitted, its fit has no minimum inside the ranges searched, its
#'     range is below the shortest lag, it found no spatial structure, or
#'     \code{sac} is a range alone), \eqn{\tau^2} is Mallows' own
#'     choice instead: the residual mean square \eqn{RSS / (m - L_m)} of the
#'     finest level whose cells hold at least two scored rows on average.  It
#'     counts the structure within those cells as noise too, so on average it
#'     is no smaller than the nugget and errs towards fewer cells; the profile warns,
#'     and \code{attr(, "cp_noise")} records which was used.  When no level's
#'     cells hold that many, \code{cp} is \code{NA} at every level, as it is
#'     for a constant response, which the profile warns about.
#'     \strong{In practice its minimum is often a bound}, the support
#'     ceiling when the nugget is a small share of the sill, and
#'     \code{select_resolution()} then names the bound in \code{edge}.  The
#'     same formula at one cell,
#'     with the total sum of squares as its RSS, is kept in
#'     \code{attr(, "one_cell")}, with the p-value at every level of the F
#'     test of that level's cell means against one mean.  A level has to
#'     beat one cell: \code{select_resolution()} returns
#'     \code{structured = FALSE} when \eqn{C_p} at one cell is within
#'     \code{tol} of the optimum, or when the cell means at the optimum
#'     differ from one mean by no more than chance allows, because one mean
#'     for the whole layer then predicts as well as the cells and the count
#'     is not identified (see \code{\link{select_resolution}} for the rule
#'     and its measured error rates).}
#'   \item{\code{moran_z}}{\strong{A diagnostic of the cell means, not a
#'     count rule.}  The standardised deviate of Moran's I on the cell
#'     means, or with \code{predictor_vars} on the residuals of the cell
#'     means regressed on the cell-mean predictors, over the 8 nearest
#'     neighbours of each cell.  It tests whether neighbouring cell means
#'     are independent.  It does not measure what the cells leave
#'     unexplained, and it does not say how well they represent the
#'     response.  On any autocorrelated response it grows as the cells get
#'     finer than the range, because neighbouring cells then resemble each
#'     other, and it has almost no
#'     power below about 20 cells, where the 8-neighbour graph between cells
#'     is close to complete and Moran's I sits near \eqn{-1/(L - 1)}
#'     whatever the data.  The level with the smallest \eqn{|z|} is
#'     therefore among the coarsest it is computed at on a structured
#'     response, whatever the layer and its ceiling, and arbitrary on a
#'     response with no structure, where
#'     \eqn{z} is calibrated and flat in \eqn{L}.
#'     \code{select_resolution()} refuses to read it as a criterion unless
#'     asked with \code{diagnostic = TRUE}, and neither \code{summary()} nor
#'     the functions that take a profile as a cell count read it by
#'     default.  \code{NA} at nine cells or fewer.}
#'   \item{\code{reliability}}{The between-cell signal's share of the spread
#'     in the cell means, from the fitted variogram alone via Krige's
#'     additivity relation (Cressie 1996), for square cells of the level's
#'     average area holding the level's average share of the layer's points
#'     with a response (larger is better).
#'     This is the shrinkage factor of Fay and Herriot (1979).  The curve
#'     has an interior optimum, and a broad one: validated against the
#'     empirical reliability of true block means on simulated fields, the
#'     analytic and empirical optima agreed to within a level or two where
#'     the empirical estimate was stable, and the band within 2 percent of
#'     the maximum spanned a factor of 3--6 in \eqn{L}.  \strong{That
#'     optimum usually lies below the range floor}: on 20 simulated fields
#'     of 900 points the curve read from the floor up was highest at the
#'     floor on 19.
#'     \code{select_resolution()} therefore reads it at the levels at or
#'     above the floor (its \code{range_floor} argument), and what it
#'     returns is then the floor, which is a level of every ladder the
#'     profile builds: the coarsest
#'     count the range allows, with the range choosing and not the curve.
#'     Read the flat region, not the argmax.  The domain term is taken over
#'     the convex hull the area is measured on, so rotating the layer does
#'     not move it.  \code{NA} without a usable variogram, and when the
#'     variogram's range was rejected (see \code{sac}).}
#'   \item{\code{loo}}{With an \code{estimate}: the leave-one-out mean
#'     squared error of that cell value (lower is better).  It estimates the
#'     squared error of predicting one further observation, at a location
#'     drawn the way the data's locations were, by the value of the cell it
#'     falls in.  That error depends on the cell value, so the column is
#'     computed for the estimate the profile is given and for no other.
#'     Each scored row of the layer is predicted by the value of its cell
#'     computed without it, in the cell of its nearest centre, which is the
#'     cell the seeds will build.  Under \code{"mean"} a row in a cell of
#'     two or more is predicted by the mean of the others, and a row alone
#'     in its cell by the mean of the cell of its nearest other scored row.
#'     Under \code{"pooled"} the cell's mean and its own weight are computed
#'     again without the row; the spread inside the cells, the spread
#'     between a cell and its target and the targets are held; and a cell
#'     left empty takes its target.  It needs no variogram.  It describes
#'     locations among the observations and says nothing of a location away
#'     from them: it is a criterion for comparing levels where data are
#'     near, not an estimate of the error of a map.  \code{loo_se} is the
#'     standard error, over rows, of that level's mean squared error.  It is
#'     not the standard error of a difference between two levels, whose
#'     errors come from the same rows.  \code{scored_per_cell} counts the
#'     rows \code{loo} read, every scored row of the layer, over the cells
#'     that hold one; \code{cell_resp_mean} counts the rows the other
#'     response criteria read, which are the subsample's when the cells are
#'     fitted to one.  \strong{In the simulation described under
#'     "Observations per cell", on a ladder of the seven sizes 18 to 2, its
#'     minimum was at the finest size on 176 of 200 layers for pooled values
#'     and on none for plain means}: for pooled values a size sets the count,
#'     and \code{select_resolution(criterion = "size")} reads it.}
#' }
#' \code{cp} and \code{reliability} answer different questions: how well the
#' cells represent the field, and whether the cell values are distinguishable
#' from noise.  They can disagree, and both are shown so the choice between
#' them is made knowingly.
#'
#' @section Conditions:
#' What changes the meaning of a result is raised as an R warning with a
#' class, so \code{tryCatch()}, \code{withCallingHandlers()} and
#' \code{options(warn = 2)} see it without matching message text; the rest
#' is logged (see \code{\link{spatialkit_quiet}}).
#' \describe{
#'   \item{\code{spatialkit_not_supported}}{\code{resolution_profile()}:
#'     \code{supported} is \code{FALSE}, because the range floor exceeds the
#'     ceiling or because fewer than two cells meet \code{min_cell_n}.}
#'   \item{\code{spatialkit_response_ceiling}}{\code{resolution_profile()}:
#'     rows that are not scored (no response, or with \code{predictor_vars}
#'     a predictor missing) lowered the ceiling.  The message names the
#'     bound in force, \code{min_cell_n} or, on a subsample, \code{sample_n};
#'     the condition carries \code{n}, \code{n_response}, \code{ceiling},
#'     \code{ceiling_from}, \code{n_sample} and \code{n_scored}.}
#'   \item{\code{spatialkit_constant_response}}{\code{resolution_profile()}:
#'     the scored variable is constant to within rounding (a response the
#'     predictors fit exactly leaves such residuals), so \code{cp},
#'     \code{moran_i} and \code{moran_z} are \code{NA} at every level.}
#'   \item{\code{spatialkit_empty_cells}}{\code{resolution_profile()},
#'     \code{\link{get_voronoi_seeds}()} and
#'     \code{\link{voronoi_seeds_kmeans}()} under \code{partition =
#'     "density"}: a centre holds no point, so there are fewer occupied cells
#'     than the count.  A centre that loses its points stays where it is, and
#'     the start that is kept is the one with the lowest sum, whether or not
#'     all its centres hold points.  It happens at high counts, mostly from
#'     \code{init = "random"}.  On simulated layers of 800 points (in 4
#'     clusters, in 28, and uniform; 20 of each) no partition into 80 cells
#'     or fewer had such a centre.  At 400 cells, from random starts with
#'     the best of 10 kept, 2, 5 and 1 of the 20 layers had one; from
#'     k-means++ starts none of 840 partitions at 40 to 600 cells had one,
#'     though it can happen.  The condition carries \code{level}, the count,
#'     and \code{empty}.}
#'   \item{\code{spatialkit_not_converged}}{\code{resolution_profile()},
#'     \code{\link{get_voronoi_seeds}()} and
#'     \code{\link{voronoi_seeds_kmeans}()}, for every call except the
#'     plain k-means each has always run: the start that was kept was cut
#'     off at the iteration limit, so its centres are those of the last step
#'     and a longer run can still move points between cells.  The condition
#'     carries \code{level} (the count; for a profile every level it
#'     happened at), \code{iter} and \code{iter_max}.  Raise
#'     \code{iter_max} (\code{kmeans_iter} in
#'     \code{get_voronoi_seeds()}).}
#'   \item{\code{spatialkit_partition_unused}}{\code{\link{get_voronoi_seeds}()}
#'     and \code{\link{voronoi_seeds_kmeans}()}, given a profile or a
#'     \code{select_resolution()} result together with a \code{partition},
#'     an \code{init} or density settings other than those the profile was
#'     built with: the scored centres are returned, and what the call asked
#'     for is not used.  The
#'     condition carries \code{kind} and \code{init}, those of the centres
#'     returned.}
#'   \item{\code{spatialkit_bound_chosen}}{\code{\link{get_voronoi_seeds}()},
#'     \code{\link{voronoi_seeds_kmeans}()} and
#'     \code{\link{build_tessellation}()}, given a profile or a
#'     \code{select_resolution()} result whose optimum sits on a bound of
#'     the ladder: the bound chose the count.  The condition carries
#'     \code{edge}, \code{criterion} and \code{n}.}
#'   \item{\code{spatialkit_no_structure}}{The same functions, given a count
#'     whose \code{structured} is \code{FALSE}: one cell does as well, so
#'     the count is not identified.  A count set by a size stands, since no
#'     criterion chose it, and the warning then says that the cell values
#'     will say little that one mean for the layer does not.
#'     \code{\link{summarize_by_cell}(estimate = "pooled")} raises the class
#'     for its own cells.}
#'   \item{\code{spatialkit_no_elbow}}{\code{\link{determine_optimal_levels}()}:
#'     the WSS curve has no elbow, or the ladder is too short to show one,
#'     so the count returned is the ladder's and not a finding.  The
#'     condition carries \code{k}, \code{sag} and \code{structured}.}
#'   \item{\code{spatialkit_diagnostic_ranking}}{\code{\link{determine_optimal_levels}()}
#'     under \code{criterion = "morans_i"} or \code{"combined"}: Moran's z
#'     of the cell means, a diagnostic and not a count rule, ordered the
#'     candidates returned.  The condition carries \code{k},
#'     \code{criterion}, \code{first_computable}, \code{window} and
#'     \code{upgraded}.}
#'   \item{\code{spatialkit_ladder_unused}}{\code{resolution_profile()}: a
#'     numeric \code{obs_per_cell} sets the ladder, and the caller also
#'     supplied \code{n_levels}, \code{min_cell_n}, \code{sample_n} or
#'     \code{range_floor = TRUE}, which are then not read.  The condition
#'     carries \code{unused} (the argument names) and \code{obs_per_cell}.}
#'   \item{\code{spatialkit_size_not_reached}}{\code{resolution_profile()}
#'     on a ladder of sizes: a size does not give the count it asks for.
#'     Either a size given as a number gives fewer than two cells while a
#'     finer one does not, and its level is 2 (\code{reason} is
#'     \code{"fewer than two cells"}); or a size asks for more cells than
#'     k-means can place, and its level is the most that can (\code{reason}
#'     is \code{"more cells than can be placed"}).  The condition carries
#'     \code{obs_per_cell},
#'     \code{wanted}, \code{used} and \code{reason}.  When the finest size
#'     gives fewer than two cells the warning is
#'     \code{spatialkit_not_supported}, with \code{obs_per_cell} among its
#'     fields and \code{min_cell_n} \code{NA}; and
#'     \code{spatialkit_response_ceiling} is not raised on such a ladder,
#'     since a size already counts the rows with a response.}
#'   \item{\code{spatialkit_repeated_locations}}{\code{resolution_profile()}
#'     with an \code{estimate}: scored rows share a location.  Each row is
#'     left out alone and the others at its location stay, so \code{loo}
#'     then describes a further observation at a place already observed.
#'     Raised once, after the levels are scored, and not when \code{loo} is
#'     \code{NA} at every level.  The condition carries \code{n} (the rows
#'     that share a location), \code{total} (the rows \code{loo} read) and
#'     \code{column}.}
#'   \item{\code{spatialkit_moran_skipped}}{\code{resolution_profile()}
#'     with a level of more than 5,000 cells, on any ladder, when
#'     \pkg{FNN} or \pkg{Matrix} is not installed: \code{moran_i} and
#'     \code{moran_z} are \code{NA} at a level where more than 5,000 cells
#'     hold a scored row, and every other column is computed.  Raised once,
#'     before the first fit, when such a level can occur: the cells that
#'     hold a scored row are counted after a level's fit, so where rows
#'     without a response leave fewer of them the two columns are computed
#'     after all.  The condition carries \code{limit}, \code{levels} (the
#'     levels above the limit) and \code{missing}.}
#' }
#' Two refusals are errors with a class: \code{spatialkit_diagnostic_criterion}
#' (\code{select_resolution(criterion = "moran_z")} without
#' \code{diagnostic = TRUE}) and \code{spatialkit_no_default_criterion} (a
#' profile passed as a cell count whose \code{cp} is \code{NA} at every
#' level while another criterion is finite).  The remaining warnings of
#' these functions (points dropped, a variogram that gives \code{cp} no
#' nugget, a zero nugget, a variogram of another variable) are R warnings
#' without a class; rows left out for a missing response, a rising WSS
#' curve and a subsample too small to reach the floor are logged only.
#'
#' A third refusal is an error of class \code{spatialkit_work_limit}:
#' \code{resolution_profile()} stops before a step whose count of
#' comparisons is above a limit.  The condition carries \code{work} (the
#' count), \code{limit} and \code{what} (the step).  On a ladder that
#' \code{obs_per_cell} sets, the count of cells grows with the rows, so two
#' costs are bounded before the variogram is estimated and before any fit.
#' The fits: the points times the cells of
#' the finest level times \code{nstart} may not pass
#' \code{getOption("spatialkit.fit_work_max", 1e10)}
#' (\code{what = "partition fit"}; the condition also carries \code{kind},
#' \code{init}, \code{nstart} and \code{fnn}).  That product counts the
#' step in which every point is compared with every centre, which
#' \code{partition = "kmeans"} has in every pass and in a k-means++ start,
#' and \code{partition = "density"} has in a k-means++ start.  Without
#' \pkg{FNN} the density kind compares every point with every centre in
#' every pass, and the product is multiplied by \code{iter_max}; with
#' \pkg{FNN} and \code{init = "random"} it has no such step and nothing is
#' counted.  The density kind's weights: without \pkg{FNN} they come from a
#' search that compares every pair of distinct locations, and the square of
#' their number may not pass \code{getOption("spatialkit.exact_pairs_max",
#' 2e9)} (\code{what = "density weights"}); \code{density_power = 0} has no
#' weights and no such search.  With an \code{estimate} and without
#' \pkg{FNN}, on any ladder, the same limit bounds the three exact searches
#' the estimate's columns make: each scored row's nearest centre, each cell's
#' nearest cells, and, under \code{"mean"}, the nearest other location of a
#' row alone in its cell.  The result of every search is the same with and
#' without \pkg{FNN}.
#'
#' The warnings of the pooled values and of the count record, which
#' \code{\link{summarize_by_cell}()} raises, are listed under "Conditions" on
#' its page.
#'
#' @param data_sf An sf object of points (other geometries are reduced to
#'   representative points).  Features with empty or non-finite coordinates
#'   are dropped with a warning.
#' @param response_var Optional response column name (numeric or logical).
#'   Enables \code{cp} and \code{moran_z}, and, when it is numeric, the
#'   pooled estimate the default profile scores by \code{loo}.  A variogram
#'   estimated from it also gives the range floor the bounds report and
#'   \code{reliability}; it is estimated on a ladder of counts or when
#'   \code{sac} asks for it, and not on the default ladder of one size, so
#'   \code{cp} and \code{reliability} are \code{NA} there.  Rows
#'   where it, or a predictor, is missing or non-finite stay in the geometry
#'   and are left out of the OLS fit, the RSS, \code{cp}, \code{moran_z} and
#'   the count the support ceiling divides by \code{min_cell_n}; a logged
#'   warning gives their number, and an R warning says so when they lower
#'   the ceiling.
#' @param predictor_vars Optional predictor column names (numeric or
#'   logical).  With them, \code{cp} scores the OLS residuals of the response
#'   on the predictors, the variogram is estimated from those residuals, and
#'   \code{moran_z} regresses the cell means on the cell-mean predictors.
#' @param levels Optional integer vector of level counts to score, replacing
#'   the ladder; values below 2, above the number of distinct locations, or
#'   at or above the number of points, are dropped (k-means cannot place more
#'   centres than there are distinct points, and \code{stats::kmeans()}
#'   refuses as many centres as points).
#' @param n_levels Number of levels on the ladder.  Default 20.  The ladder
#'   has fewer when rounding puts two of them on one count, and one more
#'   when the range floor is added as a level (see "The ladder and its
#'   bounds").  Supplying it, or \code{min_cell_n}, \code{sample_n},
#'   \code{levels} or \code{range_floor}, gives that ladder in place of the
#'   ladder of sizes \code{obs_per_cell = "auto"} would set.
#' @param min_cell_n Minimum average number of rows per cell that a level
#'   must keep (scored rows when \code{response_var} is given: those with a
#'   response and every predictor; points otherwise); sets the ceiling.
#'   Default 9.  A layer that cannot give two cells that many is profiled at
#'   two cells, or at the \code{levels} asked for, with
#'   \code{supported = FALSE} and a warning.  It is floored to a whole
#'   number, and it is not read on a ladder \code{obs_per_cell} sets (see
#'   "The ladder and its bounds").
#' @param sample_n Points are subsampled to this many before the k-means
#'   fits, as in \code{determine_optimal_levels()}, on the ladder of "The
#'   ladder and its bounds"; on a ladder of sizes, the default, the centres
#'   are fitted to every row.  Default 1500.  The
#'   columns read off the fitted cells (\code{wss}, the \code{cell_} and
#'   \code{cells_} columns, \code{rss}, \code{moran_i}, \code{moran_z})
#'   describe the subsample; the bounds, \code{supported}, the variance term
#'   of \code{cp} and \code{reliability} describe the whole layer, so the
#'   answer does not change with \code{sample_n} except through the fits.
#' @param nstart Starts per level; the best is kept.  They are k-means++
#'   starts unless \code{init} says otherwise.  \code{NULL} (default) is the
#'   kind's own number: 3 for the kind \code{"density"} and 25 for the kind
#'   \code{"kmeans"}, the number earlier versions ran.  A number is run
#'   under either \code{partition}: \code{nstart = 25} gives the density
#'   kind the starts it had when 25 was the default.
#'   \code{attr(, "nstart")} holds the number that was run.
#' @param seed RNG seed for the subsample and the restarts; restored
#'   afterwards.  Default 123.  The rows are put in coordinate order before
#'   either, so the profile does not depend on the order they come in.
#' @param sac Optional \code{sac_range} object from
#'   \code{\link{estimate_sac_range}()} to take the range, nugget and
#'   correlation function from.  Pass one fitted with \code{detrend =
#'   "reml"}, say.  It must describe the variable the profile scores: the
#'   raw response without \code{predictor_vars}, the residuals on them with;
#'   a sac whose \code{detrended} attribute says otherwise is used with a
#'   warning.  Its range is read in its own CRS (\code{attr(sac, "crs")}),
#'   to which the points are transformed first, so the area, the floor and
#'   \code{cell_diam_median} are then in that CRS's units.  A sac whose
#'   range was rejected (\code{NA} with a \code{rejected_reason}) gives
#'   \code{cp} its nugget, with a warning, and leaves \code{reliability}
#'   \code{NA}, since that needs the range; one whose model did not converge
#'   (its fit has no minimum inside the ranges searched, so its nugget and
#'   sill are not fitted values), or whose range is below the shortest lag
#'   fitted (a structure that cannot be told from a nugget, so the nugget is
#'   not identified either),
#'   gives neither, and \code{cp} takes its noise variance from the finest
#'   level (see \code{cp} above).  Whether the model converged is read from
#'   the model (\code{attr(attr(sac, "variogram_model"), "converged")}, which
#'   \code{\link{estimate_sac_range}()} sets and describes), not from the
#'   reason: the reason names one ground for the refusal, and such a fit is
#'   refused as past the largest lag.  The warning then names both.  A
#'   sac that carries no such flag (one made by hand, or built from REML
#'   parameters) is judged by its
#'   reason.  Under \code{select_on = "split"} the sac must come from
#'   the selection half alone: run the profile once without it, fit the sac
#'   on \code{data_sf[attr(p, "split")$selection, ]} and pass it to a second
#'   call with the same \code{seed}, which makes the same split.  When
#'   \code{NULL} and a response is given, one is estimated on the subsample
#'   (its selection half under \code{"split"}) with the same
#'   \code{predictor_vars}.  A plain number is taken as the range alone, in
#'   the units of the CRS the profile is computed in (metres for lon/lat
#'   input): it gives the range floor, \code{reliability}, which needs a
#'   fitted model, is \code{NA}, and \code{cp} takes its noise variance from
#'   the finest level.  A sac refused because the response shows no spatial
#'   structure (either of the two \code{rejected_reason} values of
#'   \code{\link{estimate_sac_range}()} that say so; the model it attaches
#'   is a pure nugget) gives \code{cp} no
#'   nugget either, and \code{select_resolution()} then returns
#'   \code{structured = FALSE}.  A \code{units} object is refused
#'   rather than read as a number in whatever unit it was written in.
#'   The string \code{"auto"}, the default, is \code{NULL} on every ladder
#'   but one that \code{obs_per_cell} sets.  On that ladder no variogram is
#'   estimated or
#'   read: \code{cp} and \code{reliability} are \code{NA} at every level,
#'   \code{cp} takes no noise variance from the finest level in its place,
#'   \code{attr(, "cp_noise")} is a list whose \code{source} is
#'   \code{"none"} and whose \code{reason} says why, and the
#'   \code{variogram} and \code{sac} attributes are \code{NULL}.  The
#'   other columns do not read the variogram and are the numbers
#'   \code{sac = NULL} gives.  \code{sac = NULL}, the default of earlier
#'   versions, estimates the variogram on a ladder of sizes too, and gives
#'   \code{cp} and \code{reliability} there.
#' @param select_on \code{"all"} (default) profiles every point;
#'   \code{"split"} reads the response on one spatially blocked half only
#'   (the OLS fit, the variogram, \code{rss}, \code{cp} and \code{moran_z})
#'   and returns the other half as the set to estimate on, in the
#'   \code{"split"} attribute.  The cells, \code{wss}, \code{elbow}, and the
#'   extent and point counts behind the bounds, \code{cp}'s variance term and
#'   \code{reliability} still come from every point, because the
#'   tessellation the count is for is built on every point: the levels are
#'   cell counts for the whole layer, and the estimation half's response
#'   never touches them.  See the "Post-selection inference" section of
#'   \code{\link{determine_optimal_levels}}; the profile reads the response
#'   whenever \code{response_var} is given.
#' @param range_floor \code{FALSE} (default) starts the ladder at 2 whatever
#'   the range, and the range floor is reported in the bounds and the print
#'   without being applied.  \code{TRUE} starts the ladder at the floor when
#'   the data support it and the floor is below the ceiling, as the default
#'   did up to version 2.0.0; an estimated range then sets the fewest cells
#'   any criterion can choose (see "The ladder and its bounds").
#' @param partition What the centres of the regions minimise, at every level
#'   of the ladder.  \code{"density"} (default): the sum of squared
#'   distances from the points to their nearest centre, each point weighted
#'   by its local density to the power \code{density_power}; compact cells
#'   that, at the default power, hold about the same number of points each
#'   once there are more cells than dense places, and counts about as uneven
#'   as \code{"kmeans"} gives below that.  \code{"kmeans"}: the same sum
#'   unweighted; compact cells whose counts are several times higher in
#'   dense places than in sparse ones.  \code{"kmeans"} was the default of
#'   earlier versions, and with \code{nstart = NULL} it runs their 25
#'   starts.  See "The partition".
#' @param init How each start is drawn.  \code{"kmeans++"} (default): the
#'   centres one at a time, each in proportion to its squared distance from
#'   those already drawn, times its weight under \code{partition =
#'   "density"}.  \code{"random"}: \code{k} distinct locations at random, as
#'   \code{\link[stats]{kmeans}()} draws them.
#' @param density_power For \code{partition = "density"}: the power the local
#'   density is raised to for a point's weight, a number from 0 to 2.
#'   Default 1, which evens out the cell counts once the cells outnumber the
#'   dense places (see "The partition").  \code{0} weights every
#'   point alike, the objective of \code{"kmeans"} (reached by Lloyd's
#'   iteration, so the centres are close to that kind's and not equal to
#'   them); values between give counts between the two; values above 1 give
#'   the cells of dense places fewer points than those of sparse places.
#'   Setting it under \code{partition = "kmeans"}, which does not read it,
#'   is a warning.
#' @param density_k For \code{partition = "density"}: the number of nearest
#'   neighbours the local density is read from.  Default 8.  It is held to
#'   one fewer than the distinct locations.  A small number follows the
#'   density closely and noisily; a large one smooths it.  Setting it under
#'   \code{partition = "kmeans"} is a warning.
#' @param iter_max The most iterations a start runs for: a single positive
#'   whole number, or \code{NULL} (default) for the kind's own limit, 50 for
#'   \code{partition = "kmeans"} from k-means++ starts and 100 for random
#'   starts and for \code{partition = "density"}.  A kept start that was cut
#'   off at the limit is a warning of class
#'   \code{"spatialkit_not_converged"}, except for the kind \code{"kmeans"}
#'   from k-means++ starts under \code{iter_max = NULL}, which reports
#'   nothing, as it always has.
#' @param estimate The cell value the profile scores, which is the one
#'   \code{\link{summarize_by_cell}()} is to report: \code{"auto"}
#'   (default; see below), \code{NULL} for none, \code{"mean"} for the
#'   plain mean of a cell's rows, \code{"pooled"} for the mean pulled toward
#'   the mean of the 8 nearest cells.  With an estimate the profile gains
#'   five columns and \code{attr(, "estimate")} (see "Value"), and every
#'   column and attribute it has without one is unchanged.  \code{loo} is
#'   the leave-one-out mean squared error of that estimate: each row predicted
#'   by the value of its cell computed without it.  Under \code{"mean"} a
#'   row in a cell of two or more is predicted by the mean of the others,
#'   and a row alone in its cell by the mean of the cell of its nearest
#'   other row.  Under \code{"pooled"} the cell's mean and its own weight
#'   are computed again without the row, and the spread inside the cells,
#'   the spread between a cell and its target, and the targets are held.
#'   It is computed on every row of the layer with a finite response (of the
#'   selection half under \code{select_on = "split"}), not on the subsample
#'   the centres were fitted to, and each row is scored in the cell of its
#'   nearest centre, the cell \code{\link{get_voronoi_seeds}()} hands on.
#'   It needs no variogram.  Rows at one location are left out one at a
#'   time.  \code{"auto"} is \code{"pooled"} when there is a numeric
#'   \code{response_var} with a finite value and no \code{predictor_vars},
#'   and \code{NULL} otherwise.  \code{NULL}, the default of earlier
#'   versions, scores none and gives the profile without the five columns.
#'   \code{"mean"} and \code{"pooled"} are
#'   errors without a \code{response_var}, with \code{predictor_vars} (the
#'   profile would then score residuals, and the summary reports cell values
#'   of the response), with a logical response, and with a response that has
#'   no finite value.
#' @param obs_per_cell Rows per cell: \code{"auto"} (default; see below);
#'   \code{NULL} for the ladder described under "The ladder and its
#'   bounds", the default of earlier versions, whatever the other arguments
#'   are; or one number or several, each at least 1 (fractions are
#'   allowed), for a ladder of sizes.  The level of a size \eqn{s} is
#'   \code{max(2, round(n / s))} cells, with \eqn{n} the rows that have a
#'   response (and, with \code{predictor_vars}, every predictor), or every
#'   point when there is no \code{response_var} or no row has a response.
#'   Under \code{select_on = "split"} \eqn{n} counts the rows of both
#'   halves, since the cells are built on every point.  On a ladder of sizes
#'   the centres are fitted to every row, no subsample is drawn, and
#'   \code{n_levels}, \code{min_cell_n}, \code{sample_n} and
#'   \code{range_floor} are not read: supplying \code{n_levels},
#'   \code{min_cell_n}, \code{sample_n} or \code{range_floor = TRUE} beside
#'   a numeric \code{obs_per_cell} is a warning of class
#'   \code{spatialkit_ladder_unused}, and \code{levels} beside it is an
#'   error.  A level is held to the most centres k-means can place (the
#'   distinct locations, and one short of the points).  \code{"auto"} is
#'   \code{NULL} when \code{levels} is given or the caller supplies any of
#'   those four arguments; otherwise 4.5 when the estimate is
#'   \code{"pooled"}, the seven sizes 18, 12, 9, 6, 4.5, 3 and 2 when it is
#'   \code{"mean"}, and \code{NULL} without an estimate.  The cost of a
#'   ladder of sizes is bounded before any fit (see "Conditions").
#' @return A data.frame of class \code{resolution_profile} with one row per
#'   level and columns \code{levels}, \code{wss} (the within-cell sum of
#'   squares of the coordinates, every point counted once, under either
#'   \code{partition}), \code{wss_spread} (relative
#'   spread of WSS across the restarts; under \code{partition = "density"},
#'   of the weighted sum the starts minimise), \code{elbow}, \code{cell_n_min},
#'   \code{cell_n_median}, \code{cell_diam_median} (twice the median RMS
#'   radius of the cells, in coordinate units: about 0.8 of the side of a
#'   square cell of the same area, and less for finer cells, so a size to
#'   compare levels by rather than a width), \code{rss}, \code{cp},
#'   \code{moran_i}, \code{moran_z}, \code{reliability}, \code{cells_resp}
#'   (the cells that hold at least one scored row) and \code{cell_resp_mean}
#'   (the scored rows per such cell; compare it with \code{min_cell_n}).
#'   Like the other columns read off the fitted cells, the last two count
#'   rows of the subsample when the layer is larger than \code{sample_n}:
#'   the cells of the layer then hold about \code{n_response / n_scored}
#'   times as many scored rows as \code{cell_resp_mean} shows.
#'   Columns a missing input leaves undefined are \code{NA}.  Attributes:
#'   \code{bounds} (a list with \code{floor}, \code{ceiling},
#'   \code{ceiling_from} (\code{"min_cell_n"}, \code{"min_cell_n (rows with a
#'   response)"}, \code{"distinct locations"}, \code{"sample_n"},
#'   \code{"fewer than two supported cells"} or \code{"fewer than two
#'   distinct locations"}, whichever bound it), \code{supported}, \code{area},
#'   \code{range}, \code{n} (the points in the layer), \code{n_sample} (the
#'   points the k-means fits ran on), \code{n_distinct}, \code{min_cell_n},
#'   \code{range_floor}, \code{n_response} (the rows of the layer that carry
#'   a response and, with \code{predictor_vars}, every predictor, which the
#'   support ceiling counts; \code{NA} without a
#'   \code{response_var}), \code{n_scored} (the rows of the subsample the
#'   response criteria read; \code{NA} without a \code{response_var}),
#'   \code{levels_given} (whether \code{levels} set the ladder) and
#'   \code{predictors} (whether \code{predictor_vars} were given)),
#'   \code{variogram} (a list with \code{nugget}, \code{psill}, \code{range}
#'   (\code{NA} when rejected), \code{model} and \code{detrended} (whether
#'   the sac says it is a variogram of residuals; \code{NA} when it does not
#'   say); \code{NULL} when none was usable), \code{variable} (\code{"response"},
#'   \code{"residuals"} or \code{NA}), \code{wss_bumps}, \code{nstart},
#'   \code{sac} (the range object used) and, with \code{select_on =
#'   "split"}, \code{split} (a \code{spatialkit_split}: \code{selection}
#'   and \code{estimation}, integer row positions in \code{data_sf}, with
#'   the \code{method} and \code{seed} that made them), \code{cp_noise} (a
#'   list with \code{source}, \code{"variogram nugget"} or
#'   \code{"finest-level residual mean square"}, \code{value}, the
#'   \eqn{\tau^2} used, and \code{level}, the level it was read at or
#'   \code{NA}; \code{NULL} when \code{cp} has neither), \code{one_cell}
#'   (what one cell for the whole layer scores, the reference a level has to
#'   beat: a list with \code{n}, the scored rows, \code{rss}, their total
#'   sum of squares, \code{cp}, the same formula at one cell or \code{NA}
#'   when \code{cp} has no noise variance, \code{constant}, whether the
#'   scored variable is constant to within rounding, \code{cp_levels}, the
#'   \code{cp} column it belongs to, and \code{f_p}, named by level, the
#'   p-value of the F test of that level's cell means against one mean,
#'   \eqn{(TSS - RSS) / (L - 1)} over \eqn{RSS / (m - L)} on the \eqn{L}
#'   cells that hold a scored row and the \eqn{m} scored rows; \code{NULL}
#'   without a response), \code{centres}
#'   (a list named by level of the two-column matrices of the scored
#'   partitions' centres), \code{centres_crs} (the CRS they are in) and
#'   \code{hull} (the convex hull of the points profiled, an \code{sfc} in
#'   that CRS).  A
#'   k-means partition is the Voronoi partition of its centres, so these
#'   seed the very cells each level was scored on: \code{select_resolution()}
#'   returns the chosen level's as \code{$seeds}, and
#'   \code{\link{get_voronoi_seeds}()} and \code{\link{voronoi_seeds_kmeans}()}
#'   use them when given the selection or the profile.  The centres carry the
#'   hull as their \code{"cover"} attribute, which those two functions hand on
#'   so that Voronoi cells built on the seeds reach the edge of the points.
#'   A profile built with \code{partition = "density"} or \code{init =
#'   "random"} also carries \code{partition}: a list with \code{kind},
#'   \code{init}, \code{density_power} and \code{density_k} (both \code{NA}
#'   for \code{kind = "kmeans"}); for \code{kind = "density"},
#'   \code{objective}, the weighted sum the centres minimise, named by
#'   level; \code{iter_max}, the iteration limit; and \code{iter} and
#'   \code{converged}, the iterations each level's kept start took and
#'   whether it stopped before the limit, named by level.  The selection
#'   carries the first four as \code{$partition}, its
#'   seeds and the seeds the seeding functions return as
#'   \code{attr(, "partition")}, and a Voronoi tessellation built on them as
#'   \code{params$partition}.  A profile of the kind \code{"kmeans"} from
#'   k-means++ starts has no such attribute, and neither has anything made
#'   from it.
#'
#'   With an \code{estimate} the frame has five more columns, after
#'   \code{cell_resp_mean}: \code{scored_per_cell} (the rows \code{loo}
#'   read over the cells that hold one; the layer's counterpart of
#'   \code{cell_resp_mean}, which counts the subsample's rows), \code{loo}
#'   (the leave-one-out mean squared error of the estimate; lower is
#'   better), \code{loo_se} (the standard error over rows of that level's
#'   mean squared error; it is not the standard error of a difference
#'   between two levels, whose errors come from the same rows),
#'   \code{share_alone} (the share of the scored rows that are the only
#'   scored row of their cell) and \code{own_weight} (the median over cells
#'   of the weight a pooled value gives the cell's own mean; \code{NA}
#'   under \code{"mean"}).  \code{loo}, \code{loo_se} and
#'   \code{own_weight} are \code{NA} at a level where no cell holds two
#'   scored rows, and at every level when the response is constant over the
#'   rows \code{loo} reads.  It also carries \code{estimate}: a list with
#'   \code{value} (\code{"mean"} or \code{"pooled"}), \code{neighbours}
#'   (8, or \code{NA} under \code{"mean"}), \code{response} (the column),
#'   \code{n} (the rows \code{loo} read), \code{rows} (\code{"all"} or
#'   \code{"selection half"}), \code{features} (\code{"points"}, or
#'   \code{"other"} when the layer held features that were reduced to
#'   points), \code{seed}, \code{constant} (whether the response is constant
#'   over those rows; \code{FALSE} for fewer than two), \code{loo_one_cell}
#'   and \code{loo_one_cell_se} (the same error with one mean for all those
#'   rows, the reference a level has to beat), \code{loo_levels} (the
#'   \code{loo} column those references belong to, named by level), and,
#'   named by level, \code{sigma2} (the spread of the rows about their
#'   cell's mean), \code{tau2} (the mean squared distance between a cell's
#'   level and its target; \code{NA} under \code{"mean"}), \code{cells}
#'   (the cells that hold a scored row) and \code{f_p} (the p-value of the F
#'   test of those cells' means against one mean).
#'
#'   On a ladder that \code{obs_per_cell} sets, \code{bounds} describes that
#'   ladder: \code{ceiling} is the finest level; \code{ceiling_from} is
#'   \code{"obs_per_cell"}, \code{"distinct locations"},
#'   \code{"fewer than two cells of the size"} or
#'   \code{"fewer than two distinct locations"};
#'   \code{supported} is \code{FALSE} when the finest size gives fewer than
#'   two cells or the range floor is above the finest level;
#'   \code{levels_given} is \code{TRUE}; \code{n_sample} is every point;
#'   \code{min_cell_n} holds the argument, which is not read; and two fields
#'   are added, \code{obs_per_cell} (the sizes, largest first) and
#'   \code{size_levels} (the level of each, named by size; \code{NA} for a
#'   size of \code{"auto"} that gives fewer than two cells and was left
#'   out).
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
#'   # 600 m) on a 1 km square, with a nugget of 0.6 on a unit sill.
#'   set.seed(4)
#'   n <- 400
#'   xy <- data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000))
#'   D  <- as.matrix(dist(xy))
#'   xy$z <- as.numeric(t(chol(exp(-D / 200) + diag(0.6, n))) %*% rnorm(n))
#'   pts <- st_as_sf(xy, coords = c("x", "y"), crs = 32632)
#'   # n_levels asks for the ladder of counts.
#'   prof <- resolution_profile(pts, response_var = "z", n_levels = 12)
#'   print(prof)               # one row per level; print() because only the
#'                             # last value of a braced block is shown
#'   attr(prof, "bounds")$floor   # the range floor: reported, not applied
#'   print(select_resolution(prof, criterion = "cp"))
#'
#'   # The default profile: one level at 4.5 rows per cell, scored for
#'   # pooled values, in density-weighted cells, with no variogram.
#'   prof1 <- resolution_profile(pts, response_var = "z")
#'   print(prof1)
#'   select_resolution(prof1)   # the count the size gives
#' }
#' @export
resolution_profile <- function(data_sf, response_var = NULL, predictor_vars = NULL,
                               levels = NULL, n_levels = 20L, min_cell_n = 9L,
                               sample_n = 1500L, nstart = NULL, seed = 123L,
                               sac = "auto", select_on = c("all", "split"),
                               range_floor = FALSE,
                               partition = c("density", "kmeans"),
                               init = c("kmeans++", "random"),
                               density_power = 1, density_k = 8L,
                               iter_max = NULL,
                               estimate = "auto", obs_per_cell = "auto") {
  # Which of the arguments that describe the automatic ladder the caller
  # supplied.  Read first, before any of them is altered: `obs_per_cell` is
  # resolved against them below.  A value a wrapper passes on counts as
  # supplied.
  ladder_given <- c(n_levels = !missing(n_levels), min_cell_n = !missing(min_cell_n),
                    sample_n = !missing(sample_n), range_floor = !missing(range_floor))
  select_on <- match.arg(select_on)
  partition <- match.arg(partition)
  # As match.arg() reads it, with a message that names this function when
  # the value is not one of the two.
  init <- if (identical(init, .partition_inits())) init[1L]
          else .match_init(init, "resolution_profile")
  # Whether the caller set the ladder: select_resolution() holds `reliability`
  # to the range floor only on a ladder this function built.
  levels_given <- !is.null(levels)
  if (!inherits(data_sf, "sf"))
    stop("resolution_profile(): `data_sf` must be an sf object.", call. = FALSE)
  if (!is.logical(range_floor) || length(range_floor) != 1L || is.na(range_floor))
    stop("resolution_profile(): `range_floor` must be TRUE or FALSE.", call. = FALSE)
  if (!is.numeric(min_cell_n) || length(min_cell_n) != 1L || !is.finite(min_cell_n) ||
      min_cell_n < 1)
    stop("resolution_profile(): `min_cell_n` must be a single number >= 1.", call. = FALSE)
  # NULL is the kind's own number of starts: the 25 this function has always
  # run for the kind "kmeans", .NSTART_DENSITY for the kind "density".  The
  # number it resolves to is the one the result records.
  if (is.null(nstart))
    nstart <- if (identical(partition, "density")) .NSTART_DENSITY else 25L
  if (!is.numeric(nstart) || length(nstart) != 1L || !is.finite(nstart) || nstart < 1)
    stop("resolution_profile(): `nstart` must be NULL or a single number >= 1.", call. = FALSE)
  nstart <- as.integer(nstart)
  # The kind of partition every level of the ladder is built with, and how its
  # starts are drawn (R/partition.R).  `part` is the record the profile, its
  # selections and the seeds made from them carry.
  part <- .check_partition_args(partition, init, density_power, density_k,
                                "resolution_profile")
  # The present kind, the best of `nstart` k-means++ starts, leaves no record:
  # a profile made with it is the object it always was.
  part_default <- identical(part$kind, "kmeans") && identical(part$init, "kmeans++")
  # Iterations per start, unless `iter_max` sets them.  50 for k-means from
  # k-means++ starts, as always.  Random starts lie further from an optimum,
  # and Lloyd's iteration of the density kind takes more steps than Hartigan
  # and Wong's: 100 for both.  Over 192 density runs on 300 to 6,000 points
  # (uniform and clustered, 2 to 667 cells) the longest took 66 steps to stop
  # moving, and 5 percent took more than 40.  Larger layers take more: a
  # start that is cut off is reported after the sweep.
  if (!is.null(iter_max))
    iter_max <- .check_whole_count(iter_max, "iter_max", "resolution_profile",
                                   "iterations", null_ok = TRUE)
  part_iter <- if (!is.null(iter_max)) iter_max else if (part_default) 50L else 100L
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
  # as.numeric() strips a `units` object to its number in whatever unit it
  # was written in, so set_units(1.5, "km") became a range of 1.5 CRS units
  # (metres): a floor of 41 million cells, with nothing said.  Refused by
  # name, as fold_separation() and cv_block_size_sweep() refuse it; and a
  # character string, which as.numeric() also read, with it.
  # "auto" is resolved once the ladder is known: no variogram on a ladder
  # `obs_per_cell` sets, and NULL on every other.  Until then it is NULL, so
  # that nothing below reads it as a variogram the caller supplied.
  sac_auto <- .is_auto(sac)
  if (!is.null(sac) && !sac_auto && (inherits(sac, "units") || !is.numeric(sac)))
    stop("resolution_profile(): `sac` must be NULL, \"auto\", an estimate_sac_range() ",
         "result, or a plain number in the units of the CRS the profile is computed ",
         "in (metres for lon/lat input); got ",
         if (inherits(sac, "units")) format(sac)
         else sprintf("an object of class %s", class(sac)[1L]),
         ".", call. = FALSE)
  if (sac_auto) sac <- NULL
  # `estimate` and `obs_per_cell`: what the arguments alone decide is checked
  # here, after every check above; what needs the response, once it is read.
  # A string is read by its value alone, as match.arg() reads `partition`:
  # one taken from a named vector (cfg["estimate"]) or carrying a class is
  # the string it holds.  Every comparison below is identical(), which a
  # name or a class defeats, so both are dropped here, once.
  if (is.character(estimate)) estimate <- as.character(unclass(estimate))
  if (!is.null(estimate) &&
      !(is.character(estimate) && length(estimate) == 1L && !is.na(estimate) &&
          estimate %in% c("mean", "pooled", "auto")))
    stop(sprintf(paste0("resolution_profile(): `estimate` must be NULL, \"mean\", ",
                        "\"pooled\" or \"auto\"; got %s."), .got_text(estimate)),
         call. = FALSE)
  est_named <- !is.null(estimate) && !identical(estimate, "auto")
  if (est_named && !has_resp)
    stop("resolution_profile(): `estimate` names the value of a response, so it needs ",
         "a `response_var`.", call. = FALSE)
  if (est_named && has_pred)
    stop("resolution_profile(): `estimate` cannot be combined with `predictor_vars`. ",
         "With predictors the profile scores the residuals of the response on them, ",
         "and the summary reports region values of the response itself, so the ",
         "criterion would be computed for a number nobody reports. Drop ",
         "`predictor_vars`, or leave `estimate` NULL.", call. = FALSE)
  obs_auto <- .is_auto(obs_per_cell)
  # Plain numbers: a `units` object is numeric too, and comparing one with 1
  # stops in that package's words, which name neither argument nor function.
  if (!is.null(obs_per_cell) && !obs_auto &&
      !(is.numeric(obs_per_cell) && !is.object(obs_per_cell) &&
          length(obs_per_cell) >= 1L &&
          all(is.finite(obs_per_cell)) && all(obs_per_cell >= 1)))
    stop(sprintf(paste0("resolution_profile(): `obs_per_cell` must be NULL, \"auto\", or ",
                        "numbers of at least 1 (observations per cell; fractions are ",
                        "allowed); got %s."), .got_text(obs_per_cell)), call. = FALSE)
  if (!is.null(obs_per_cell) && !obs_auto && levels_given)
    stop("resolution_profile(): pass one of `obs_per_cell` and `levels`: each sets the ",
         "ladder.", call. = FALSE)

  # Whether the layer came as points: a feature that is not one is scored at
  # the point it is reduced to, and the record of an estimate says so.
  was_points <- all(sf::st_geometry_type(data_sf, by_geometry = TRUE) == "POINT")
  if (!was_points)
    data_sf <- coerce_to_points(data_sf, "auto")
  data_sf <- ensure_projected(data_sf)
  xy <- sf::st_coordinates(data_sf)[, 1:2, drop = FALSE]
  ok_xy <- is.finite(xy[, 1]) & is.finite(xy[, 2])
  if (!all(ok_xy)) {
    # An R warning, as in determine_optimal_levels(): a log line alone could
    # not be caught, and spatialkit_quiet() hid it.
    .warn_and_log("resolution_profile(): dropping %d point(s) with empty or non-finite coordinates.",
                  sum(!ok_xy))
    data_sf <- data_sf[ok_xy, , drop = FALSE]
    xy <- xy[ok_xy, , drop = FALSE]
  }
  if (nrow(xy) < 3L)
    stop("resolution_profile(): fewer than three points.", call. = FALSE)

  # Sample splitting, on the full layer before the subsample, so the
  # positions index `data_sf` as passed (rows dropped above excepted).
  # The cells are geometry, and the tessellation the count is for is built on
  # every point, so the partitions, the WSS and the bounds keep the whole
  # layer; only what reads the response is held to the selection half
  # (`in_sel`, below).  The profile used to run on the half alone, and the
  # count it chose for the half's extent, point count and clusters was then
  # applied to the whole layer.
  split <- NULL
  in_sel <- NULL
  if (identical(select_on, "split")) {
    split <- .spatial_half_split(data_sf, seed = seed, caller = "resolution_profile")
    in_sel <- seq_len(nrow(xy)) %in% split$selection
    if (!all(ok_xy)) {
      # Positions refer to the layer after the drop; map them back.
      kept <- which(ok_xy)
      split$selection  <- kept[split$selection]
      split$estimation <- kept[split$estimation]
    }
    # The split cannot tell which rows a supplied variogram was fitted on, and
    # the floor, Cp's penalty and reliability all come from it: one fitted on
    # the whole layer carries the estimation half's response into the choice.
    if (!is.null(sac))
      .log_warn(paste0("resolution_profile(): select_on = \"split\" with a supplied ",
                       "`sac`: it sets the floor, Cp's penalty and reliability, so it ",
                       "must be fitted on the selection half alone, ",
                       "data_sf[attr(x, \"split\")$selection, ] (the split is the same ",
                       "for the same layer and seed). One fitted on the whole layer lets ",
                       "the estimation half's response steer the count."))
  }

  # A supplied variogram's range is a length in the CRS it was fitted in
  # (estimate_sac_range() records it as attr(sac, "crs")), and the area, the
  # floor and every distance below are measured in data_sf's.  Read in another
  # projection it moved the floor, the verdict and the reliability optimum
  # without a word: a range in US feet put the floor of a metre layer at 2
  # where the same range in metres put it at 8.  Work in the variogram's CRS,
  # as summarize_by_cell() does.  After the split, so the halves are the ones
  # a first call without the sac returned; a sac with no recorded, or a
  # geographic, CRS is taken as it is.
  sac_crs <- if (!is.null(sac)) attr(sac, "crs") else NULL
  if (inherits(sac_crs, "crs") && !is.na(sac_crs) && !isTRUE(sf::st_is_longlat(sac_crs)) &&
      !is.na(sf::st_crs(data_sf)) && !isTRUE(sf::st_crs(data_sf) == sac_crs)) {
    .log_info("resolution_profile(): working in the CRS of `sac` (%s), in which its range is a length.",
              sac_crs$input %||% "unnamed")
    data_sf <- .transform_or_stamp(data_sf, sac_crs, what = "data_sf",
                                   caller = "resolution_profile")
    xy <- sf::st_coordinates(data_sf)[, 1:2, drop = FALSE]
    if (!all(is.finite(xy)))
      stop("resolution_profile(): some points have no coordinates in the CRS of `sac`.",
           call. = FALSE)
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

  # The estimate the profile scores, now that the response is read.  "auto"
  # is the pooled value wherever naming it would not be an error, and no
  # estimate otherwise, so a call that runs without the argument runs with
  # "auto".
  est <- estimate
  if (!is.null(est)) {
    resp_logical <- has_resp && is.logical(df[[response_var]])
    resp_any     <- has_resp && any(is.finite(resp))
    if (identical(est, "auto")) {
      est <- if (has_resp && !has_pred && !resp_logical && resp_any) "pooled" else NULL
    } else {
      if (resp_logical)
        stop(sprintf(paste0("resolution_profile(): `estimate` scores the cell value ",
                            "summarize_by_cell() reports, and that function summarises ",
                            "numeric columns only; '%s' is logical. Convert it with ",
                            "as.numeric() first, or leave `estimate` NULL."), response_var),
             call. = FALSE)
      if (!resp_any)
        stop(sprintf(paste0("resolution_profile(): `estimate` needs rows with a value of ",
                            "'%s', and none has a finite one."), response_var), call. = FALSE)
    }
  }
  # The ladder `obs_per_cell` sets: sizes, largest first, or none.  "auto"
  # gives way to a ladder the caller describes (`levels`, or any of the four
  # arguments of the automatic ladder), and is otherwise the size that goes
  # with the estimate: one for pooled values, seven for plain means, none
  # without an estimate.
  # (`ladder_sizes`, not `sizes`: the sweep below keeps the points per cell
  # of a level under that name.)
  ladder_sizes <- NULL
  sizes_auto <- FALSE
  if (obs_auto) {
    if (!levels_given && !any(ladder_given)) {
      ladder_sizes <- if (identical(est, "pooled")) .OBS_PER_CELL_DEFAULT
                      else if (identical(est, "mean")) .OBS_PER_CELL_LADDER
      sizes_auto <- !is.null(ladder_sizes)
    }
  } else if (!is.null(obs_per_cell)) {
    ladder_sizes <- sort(unique(as.numeric(obs_per_cell)), decreasing = TRUE)
    # `range_floor = FALSE` asks for what a ladder of sizes does anyway.
    unused <- names(ladder_given)[ladder_given &
                                    (names(ladder_given) != "range_floor" | range_floor)]
    if (length(unused))
      .warn_classed("spatialkit_ladder_unused",
                    paste0("resolution_profile(): %s %s not read when `obs_per_cell` sets ",
                           "the ladder: the levels are the sizes asked for, and the centres ",
                           "are fitted to every row. Leave `obs_per_cell` NULL for the ",
                           "ladder %s."),
                    .and_text(sprintf("`%s`", unused)),
                    if (length(unused) > 1L) "are" else "is",
                    if (length(unused) > 1L) "they describe" else "it describes",
                    data = list(unused = unused, obs_per_cell = ladder_sizes))
  }
  size_ladder <- !is.null(ladder_sizes)
  # No variogram is estimated or read on a ladder of sizes under sac = "auto".
  no_vg <- sac_auto && size_ladder

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
                       "of the OLS fit, the RSS, Cp, Moran's z and the rows the ",
                       "support ceiling counts."),
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
  # The support ceiling counts the rows that carry a response (`n_support`),
  # as Cp's variance term and the reliability do: a cell's mean is built from
  # those rows.  It used to count every point, so 900 points with a response
  # on 225 got a ceiling of 100 "from min_cell_n = 9" at 2.25 scored rows per
  # cell, and at a Cp pick of 84 cells summarize_by_cell() returned 4 cells
  # with no mean and 21 with no interval.  A geometry-only profile counts
  # every point, and so does one in which no row has a response: nothing can
  # be scored there, and the geometry columns keep their ladder.
  n_all  <- nrow(xy)
  n_resp <- if (has_resp) sum(resp_ok) else n_all
  n_support <- if (has_resp && n_resp > 0L) n_resp else n_all
  hull <- sf::st_convex_hull(sf::st_union(sf::st_geometry(data_sf)))
  area <- suppressWarnings(as.numeric(sf::st_area(hull)))
  bb   <- sf::st_bbox(data_sf)
  bbw  <- as.numeric(bb["xmax"] - bb["xmin"]); bbh <- as.numeric(bb["ymax"] - bb["ymin"])
  if (!is.finite(area) || area <= 0) area <- bbw * bbh
  # duplicated() on a complex vector, because unique() on an n x 2 matrix
  # pastes every row into a string: 0.01 s against 5 s at a million points.
  n_uniq_all <- sum(!duplicated(complex(real = round(xy[, 1], 8),
                                        imaginary = round(xy[, 2], 8))))

  # Subsample, keeping everything aligned, from the rows in a canonical order:
  # the subsample and every k-means++ draw index rows, so the same layer with
  # its rows permuted used to get different cells and a different count (a
  # WSS 2.6 percent apart and Cp picks of 222, 173 and 135 cells over
  # permutations of one 2000-point layer).  Coordinates first; the response
  # and the predictors break ties between repeat visits to one location.
  # Every point is reordered even with no subsample to draw.
  ord <- do.call(order, c(list(xy[, 1], xy[, 2]),
                          if (has_resp) list(resp),
                          if (has_pred) lapply(seq_len(ncol(pred)), function(j) pred[, j])))
  # What the estimate's columns read: every row of the layer with a finite
  # response (of the selection half under select_on = "split"), in that
  # order.  Kept here, before the subsample overwrites the layer's
  # coordinates, response and flags with its own.
  lay <- NULL
  if (!is.null(est)) {
    scored <- resp_ok[ord]
    if (!is.null(in_sel)) scored <- scored & in_sel[ord]
    lay <- list(xy = unname(xy[ord[scored], , drop = FALSE]), y = resp[ord[scored]])
  }
  # A size is a number of rows of the layer, so the centres of a ladder of
  # sizes are fitted to every row: no subsample is drawn, and `sample_n` is
  # not read.
  idx <- if (!size_ladder && n_all > sample_n) ord[sample(seq_len(n_all), sample_n)] else ord
  xy <- xy[idx, , drop = FALSE]
  data_sf <- data_sf[idx, , drop = FALSE]
  if (has_resp) { resp <- resp[idx]; resp_ok <- resp_ok[idx] }
  if (has_pred) pred <- pred[idx, , drop = FALSE]
  if (!is.null(in_sel)) in_sel <- in_sel[idx]
  n <- nrow(xy)

  # The variable the cells have to represent: the response, or what the
  # predictors leave of it.  Rows with a non-finite value drop out of the RSS
  # and the Moran statistic but stay in the geometry, and so, under
  # select_on = "split", do the rows of the estimation half.
  variable <- NA_character_
  y <- NULL
  if (has_resp) {
    use <- if (is.null(in_sel)) resp_ok else resp_ok & in_sel
    y <- rep(NA_real_, n)
    y[use] <- resp[use]
    variable <- "response"
    if (has_pred) {
      # Fitted on the complete rows only.  lm.fit() refuses any NA, and a
      # first fit on every row turned one missing value into "the OLS fit
      # failed": the raw response, trend and all, was then scored against a
      # variogram estimate_sac_range() had fitted to the residuals.
      fit <- if (sum(use) > ncol(pred) + 1L)
        try(stats::lm.fit(x = cbind(1, pred[use, , drop = FALSE]), y = resp[use]),
            silent = TRUE)
      if (!is.null(fit) && !inherits(fit, "try-error")) {
        y[use] <- fit$residuals
        variable <- "residuals"
      } else {
        .log_warn(paste0("resolution_profile(): the OLS fit on `predictor_vars` failed ",
                         "on the %d complete row(s); scoring the raw response."),
                  sum(use))
      }
    }
  }

  # The distinct locations.  `n_uniq_all` is the layer's (see above);
  # `n_uniq` is the subsample's, the locations k-means can place centres on.
  n_uniq  <- nrow(unique(round(xy, 8)))
  if (n == n_all) n_uniq_all <- n_uniq
  # k-means can put a centre on every distinct location but no more, and
  # stats::kmeans() refuses as many centres as points (see the present
  # ladder's ceiling, below).
  k_cap_all    <- min(n_uniq_all, n_all - 1L)

  # The two optional backends.  Each is asked for once, when a step first
  # needs the answer, and every later step is handed that answer: the fits,
  # each row's nearest centre, the estimate's fit, Moran's z and the guards.
  # These are the answers each of those steps reads by default, and a step
  # that never reads one never asks: a geometry-only profile of the kind
  # "kmeans" loads neither package, as it never has.
  backend <- local({
    known <- c(FNN = NA, Matrix = NA)
    function(pkg) {
      if (is.na(known[[pkg]])) known[[pkg]] <<- isTRUE(.has_backend(pkg))
      known[[pkg]]
    }
  })
  # Without FNN each scored row is compared with every centre of a level to
  # find its own.  The finest level is the largest such search.
  centre_guard <- function(finest)
    if (!is.null(est) && !backend("FNN"))
      .work_guard(as.numeric(length(lay$y)) * finest, .exact_pairs_max(),
                  what = "nearest centre", caller = "resolution_profile",
                  search = sprintf("finding the nearest of %d centres for each of %d rows",
                                   finest, length(lay$y)))
  # On a ladder of sizes the count grows with the rows, so two costs that the
  # subsample bounds on every other ladder are bounded here: the search
  # behind the density kind's weights without FNN, which compares every pair
  # of distinct locations, and the fits themselves.  The finest level follows
  # from the rows and the finest size alone (.size_ladder() gives the same
  # number below), so the guards run before the variogram is estimated: a
  # call one of them stops has fitted nothing and raised no warning about a
  # column it never returns.
  if (size_ladder) {
    finest <- .finest_level(n_support, ladder_sizes, k_cap_all)
    if (identical(part$kind, "density") && density_power > 0 && !backend("FNN"))
      .work_guard(as.numeric(n_uniq_all)^2, .exact_pairs_max(), what = "density weights",
                  caller = "resolution_profile",
                  search = sprintf(paste0("the neighbour search behind the weights of ",
                                          "partition = \"density\", over %d distinct locations,"),
                                   n_uniq_all))
    .fit_guard(n_all, finest, nstart, part, part_iter, backend("FNN"))
    centre_guard(finest)
  }

  # Variogram: supplied, or estimated on the subsample (its selection half
  # under select_on = "split": the range reads the response).
  vg <- NULL
  # The messages below name the variogram as the caller knows it: their own
  # `sac`, or the one estimated here.  They said "`sac`" either way, and
  # told a caller who had passed none to pass a different one.
  sac_given <- !is.null(sac)
  sac_what  <- if (sac_given) "`sac`"
               else "the variogram estimated from `response_var` (attr(<profile>, \"sac\"))"
  if (!no_vg && is.null(sac) && has_resp && requireNamespace("gstat", quietly = TRUE)) {
    # Whatever the estimate prints to standard output is held here and logged
    # once: neither suppressWarnings() nor spatialkit_quiet() reaches a
    # print().  estimate_sac_range() fits its model without gstat's iterative
    # fit, whose "a possible solution MIGHT be to scale semivariances and/or
    # distances" this was written for (nine times on one profile of 2,634
    # points), so nothing is expected here now.
    printed <- character(0)
    sac <- tryCatch({
      sac_est <- NULL
      printed <- utils::capture.output(
        sac_est <- suppressWarnings(estimate_sac_range(
          if (is.null(in_sel)) data_sf else data_sf[in_sel, , drop = FALSE],
          response_var,
          predictor_vars = if (has_pred) predictor_vars else NULL,
          seed = seed)),
        type = "output")
      sac_est
    }, error = function(e) {
      .log_warn("resolution_profile(): estimate_sac_range() failed (%s); no variogram-based columns.",
                conditionMessage(e))
      NULL
    })
    if (length(printed))
      .log_info("resolution_profile(): the variogram fit printed %d line(s) to the console, held back here: %s",
                length(printed), paste(unique(printed), collapse = " | "))
  }
  range_eff <- NA_real_
  if (!is.null(sac)) {
    vm <- attr(sac, "variogram_model")
    r  <- suppressWarnings(as.numeric(sac))[1L]
    if (is.finite(r) && r > 0) range_eff <- r
    # A REJECTED range (NA with a `rejected_reason`) still carries a model,
    # as evidence for the rejection, and that model used to be read as if it
    # had been accepted: reliability then took a correlation function whose
    # range could be many times the extent, which pins its optimum to the
    # first level, and nothing said so.  Here the nugget, fitted at the short
    # lags, still serves Cp; the correlation function, which is the
    # unidentified range, does not.  Three refusals are whole.  A model
    # marked as not converged: its nugget and sill are where the fit stopped
    # without a minimum.  A range below the shortest lag: the structure
    # cannot be told from a nugget, so the split between the two is exactly
    # what is not identified.  On white noise detrended by REML that model's
    # nugget was 6e-7 on a sill of 0.99, and Cp, with no penalty, ran to the
    # support ceiling (7 cells with the nugget at 0.98).  And a layer with no
    # spatial structure, further down.
    rejected <- attr(sac, "rejected_reason")
    rejected <- if (!is.finite(r) && length(rejected)) as.character(rejected)[1L] else NULL
    # Whether the model converged is read off the model, not off the reason.
    # estimate_sac_range() gives a refusal one reason and marks the model
    # besides: attr(model, "converged") is FALSE when the least-squares
    # criterion is still falling at an edge of the ranges searched (a
    # variogram with no sill inside four cutoffs), and such a fit is refused
    # as "fitted range exceeds the largest lag fitted".  Its nugget is then
    # that edge and not a fitted value: with an east-west trend on an
    # exponential field (300 points, 25 draws) Cp was given nuggets of
    # 1.18-1.62 where the field's was 0.1.  See .vgm_converged().  A sac that
    # carries no flag (made by hand, built from REML parameters, or saved
    # before the flag existed) is judged by its reason; and a reason that
    # says the fit did not converge, which only a stored object carries, is
    # believed whatever the flag says.
    not_converged <- !is.null(rejected) &&
      (isFALSE(.vgm_converged(vm)) || grepl("converge", rejected, fixed = TRUE))
    cor_fn <- if (is.data.frame(vm)) .vgm_correlation_fn(vm) else NULL
    if (!is.null(cor_fn)) {
      vg <- list(nugget = .vgm_nugget_of(vm),
                 psill  = sum(as.numeric(vm$psill[as.character(vm$model) != "Nug"])),
                 range  = range_eff, model = vm, cor_fn = cor_fn,
                 detrended = if (length(attr(sac, "detrended")) == 1L)
                   as.logical(attr(sac, "detrended")) else NA)
      if (!is.finite(vg$nugget)) vg$nugget <- 0
      if (!is.finite(vg$psill) || vg$psill <= 0) vg <- NULL
    }
    if (!is.null(vg) && !is.null(rejected)) {
      if (not_converged) {
        # Both facts are named: why the range was refused, and that the fit
        # behind it did not converge, which the reason leaves out whenever it
        # names the other ground.  The close says what `cp` is left with in
        # every case.  The fallback further down reads the residual mean
        # square of the finest level whose cells hold at least two scored
        # rows on average, and where no level does (only finer `levels` were
        # asked for, or the response is missing on most rows) `cp` is NA at
        # every level.  "From the finest level" alone was untrue there, and
        # more sacs arrive here now that the model's own verdict is read: one
        # refused as past the lags gave `cp` 0.82 and 0.87 from its nugget of
        # 0.5 at 200 and 250 cells of 300 points, and gives NA.
        .warn_and_log(paste0("resolution_profile(): %s reports no usable range (%s)%s; ",
                             "its nugget and sill are where the fit stopped without ",
                             "reaching a minimum, not ",
                             "fitted values, so neither is used: `reliability` is NA and ",
                             "`cp` takes its noise variance from the finest level whose ",
                             "cells hold at least two scored rows on average, and is NA ",
                             "at every level when none does."),
                      sac_what, rejected,
                      if (grepl("converge", rejected, fixed = TRUE)) ""
                      else ", and its variogram model did not converge")
        vg <- NULL
      } else if (grepl("shortest lag", rejected, fixed = TRUE)) {
        .warn_and_log(paste0("resolution_profile(): %s reports no usable range (%s): ",
                             "a structure that dies out before the closest pairs of ",
                             "points cannot be told from a nugget, so its nugget (%.3g) ",
                             "is not identified either: `reliability` is NA and `cp` ",
                             "takes its noise variance from the finest level."),
                      sac_what, rejected, vg$nugget)
        vg <- NULL
      } else if (.sac_no_structure(sac)) {
        # A variogram refused because the response shows no spatial structure
        # (.sac_no_structure(), R/utils.R).  estimate_sac_range() attaches a
        # pure nugget to that refusal, which has no correlation function and
        # does not reach this branch: it takes the log line below.  This is
        # for a sac that says so and still carries a model with a range (one
        # built by hand, or stored by a version that attached the fitted
        # model).  On white noise such a model's nugget is whatever a
        # spurious short range left of the variance: 0 to 0.65 against a true
        # 1 on simulated layers, with which Cp ran to 24 to 100 cells of 100.
        # Its nugget is no more identified than its range, so Cp takes its
        # noise variance from the finest level, as for a range below the
        # shortest lag.
        .warn_and_log(paste0("resolution_profile(): %s reports no usable range (%s), ",
                             "so its nugget (%.3g) is not a fitted value either: ",
                             "`reliability` is NA and `cp` takes its noise variance ",
                             "from the finest level."),
                      sac_what, rejected, vg$nugget)
        vg <- NULL
      } else {
        .warn_and_log(paste0("resolution_profile(): %s reports no usable range (%s). ",
                             "`cp` still takes its nugget (%.3g), fitted at the short lags; ",
                             "`reliability` needs the range and is NA."),
                      sac_what, rejected, vg$nugget)
        vg$cor_fn <- NULL
      }
    } else if (is.null(vg)) {
      if (.sac_no_structure(sac))
        .log_info(paste0("resolution_profile(): %s found no spatial structure (%s), ",
                         "and its model is a pure nugget; `reliability` is NA and ",
                         "`cp` takes its noise variance from the finest level."),
                  sac_what, rejected %||% "no reason recorded")
      else
        .log_info("resolution_profile(): the variogram carries no usable model; `reliability` is NA and `cp` takes its noise variance from the finest level.")
    }
    # A nugget of 0 is usually the fit's lower bound, not an estimate
    # (estimate_sac_range() marks it as `nugget_at_bound`), and with it Cp's
    # penalty is 0: Cp
    # is the RSS alone and falls with every finer level to the ceiling,
    # whatever the field.  Zero to within 1e-4 of the sill: a REML fit stops
    # short of its bound, and a nugget of 6e-7 on a sill of 0.99 slipped past
    # a test for exactly 0.
    if (has_resp && !is.null(vg) && vg$nugget <= 1e-4 * (vg$nugget + vg$psill))
      .warn_and_log(paste0("resolution_profile(): the variogram's nugget is 0 (%.2g ",
                           "against a sill of %.3g), so Cp has ",
                           "no penalty and falls as the cells get finer: its minimum ",
                           "will sit at or near the finest level whatever the field. ",
                           "A zero nugget is usually the fit's lower bound; %s"),
                    vg$nugget, vg$nugget + vg$psill,
                    if (sac_given) "check it with plot() on the sac, or pass one whose nugget you trust."
                    else paste0("check it with plot(attr(<profile>, \"sac\")), or pass a ",
                                "`sac` whose nugget you trust."))
    # The variogram has to describe the variable the profile scores.  A
    # residual variogram's nugget put into Cp for the raw response (or the
    # reverse) moved the Cp pick from 2-4 cells to the ceiling of 44 in five
    # of five simulated fields, silently.  A hand-built sac that does not say
    # which it is cannot be checked.
    if (!is.null(vg) && !is.na(vg$detrended) && !is.na(variable) &&
        !identical(vg$detrended, identical(variable, "residuals")))
      .warn_and_log(paste0("resolution_profile(): %s is a variogram of %s, but the ",
                           "profile scores %s, so `cp` and `reliability` take their ",
                           "nugget, sill and range from a different variable. Pass a ",
                           "`sac` fitted to %s."),
                    sac_what,
                    if (vg$detrended) "detrended residuals" else "the raw response",
                    if (identical(variable, "residuals"))
                      "the OLS residuals on `predictor_vars`" else "the raw response",
                    if (identical(variable, "residuals"))
                      "the same residuals (predictor_vars as here)"
                    else "the raw response (no predictor_vars)")
  }

  # Bounds.  `area`, `n_uniq_all`, `n_all` and `k_cap_all` are the layer's
  # and `n_uniq` the subsample's (see above).
  # The rows of the subsample the response criteria read (under select_on =
  # "split", those of the selection half): what a fitted cell is scored on.
  n_scored <- if (has_resp) sum(is.finite(y)) else NA_integer_
  # What a row has to carry to be scored, in the words of the messages: with
  # predictors a row with a response and no predictor is left out as well,
  # and "450 of the 900 points carry a response" was said of a layer on
  # which every point had one.
  rows_carry <- if (has_pred) "a response and every predictor" else "a response"
  # Kept as a double until the comparison: a short range on a continental
  # extent puts area / range^2 past .Machine$integer.max, and as.integer() of
  # that is NA, which turned the "not supported" branch into an abort.
  floor_raw <- if (is.finite(range_eff) && range_eff > 0)
    max(2, ceiling(area / range_eff^2)) else 2
  floor_L   <- if (floor_raw <= .Machine$integer.max) as.integer(floor_raw) else NA_integer_
  support_what <- if (n_support < n_all) paste("rows with", rows_carry) else "points"
  if (size_ladder) {
    # A ladder of sizes: the levels are the counts the sizes give, scored as
    # given, and the finest size takes the part `min_cell_n` has in the
    # ceiling, the verdict and the warning (.size_ladder()).
    sl <- .size_ladder(ladder_sizes, sizes_auto, n_support, n_all, n_uniq_all, k_cap_all,
                       floor_raw, range_eff, support_what, rows_carry)
    levels       <- sl$levels
    supported    <- sl$supported
    ceiling_L    <- sl$ceiling
    ceiling_from <- sl$ceiling_from
  } else {
    by_support <- floor(n_support / min_cell_n)
    # The ceiling is the smallest of what the layer's point count supports,
    # what its distinct locations allow and, when the fits run on a subsample,
    # what the subsample can fit: two of its points per cell.  k-means can put
    # a centre on every distinct location but no more, and stats::kmeans()
    # refuses as many centres as points, so that bound is the distinct
    # locations when some repeat and one short of the points when none do.  It
    # used to be one short of the distinct locations in every case, which
    # dropped k = 5 on five stations visited thirty times each.  Which bound
    # holds is recorded, because the print and the "bound is choosing" notes
    # name it.
    k_cap_fit    <- min(n_uniq, n - 1L)
    data_ceiling <- min(by_support, k_cap_all)
    # Two scored subsample rows per cell, not two subsample points: with 20,000
    # points, a response on 5,000 and sample_n = 1500, half the subsample left
    # about 0.5 scored rows per cell at a ceiling of 750.
    fit_rows     <- if (has_resp && n_resp > 0L) n_scored else n
    fit_cap      <- if (n < n_all) min(floor(fit_rows / 2), k_cap_fit) else Inf
    # A ladder the caller set is scored as given.  Settled here, before the
    # messages below, because they name it: with `levels = 2:4` on ten points
    # the warning said "the profile scores 2 cells" over a table of three
    # levels.
    if (levels_given) {
      levels <- sort(unique(as.integer(levels)))
      levels <- levels[!is.na(levels) & levels >= 2L & levels <= k_cap_fit]
      if (!length(levels))
        stop("resolution_profile(): no usable value in `levels`.", call. = FALSE)
    }
    # Fewer than two cells that meet min_cell_n, or one distinct location, is
    # no ladder at all.  The ceiling used to be raised to 2 without a word:
    # ten points at min_cell_n = 9 came back as "ceiling 2 from min_cell_n",
    # supported.
    too_few <- data_ceiling < 2
    ceiling_L    <- max(2L, as.integer(min(data_ceiling, fit_cap)))
    ceiling_from <- if (k_cap_all < 2L) "fewer than two distinct locations"
                    else if (by_support < 2) "fewer than two supported cells"
                    else if (fit_cap < data_ceiling) "sample_n"
                    else if (by_support <= k_cap_all)
                      (if (n_support < n_all) "min_cell_n (rows with a response)" else "min_cell_n")
                    else "distinct locations"
    # Said when the rows that are not scored moved the ceiling, in terms of
    # the bound in force.  It used to name min_cell_n whichever bound held:
    # with 20,000 points, 5,000 scored and sample_n = 1500 it read "185 cells,
    # not the 750 that every point at min_cell_n = 9 would give", where both
    # numbers are sample_n's (half the 370 scored subsample rows, and half the
    # subsample) and min_cell_n allows 555 and 2,222.
    if (n_support < n_all && !too_few) {
      all_fit_cap <- if (n < n_all) min(floor(n / 2), k_cap_fit) else Inf
      all_by      <- min(floor(n_all / min_cell_n), k_cap_all)
      all_ceiling <- max(2L, as.integer(min(all_by, all_fit_cap)))
      all_what    <- if (all_fit_cap < all_by) "half the subsample"
                     else sprintf("every point at min_cell_n = %d", min_cell_n)
      # On a subsample `cell_resp_mean` counts subsample rows, and the cells of
      # the layer hold more by this factor (12.3 against 167 in the case above).
      per_layer <- if (n < n_all && isTRUE(n_scored > 0L)) n_support / n_scored else NA_real_
      if (all_ceiling != ceiling_L) {
        msg <- if (identical(ceiling_from, "sample_n"))
          sprintf(paste0("resolution_profile(): %d of the %d points of the subsample carry %s ",
                         "(%d of the layer's %d), so the ceiling is %d cells, two scored ",
                         "subsample rows per cell, not the %d that %s would give. `sample_n` ",
                         "sets this ceiling, not `min_cell_n`: the layer's %d rows with %s ",
                         "support %d cells at min_cell_n = %d, and a larger `sample_n` raises ",
                         "the ceiling towards that. `cell_resp_mean` counts the scored rows ",
                         "of the subsample per occupied cell; the cells of the layer hold ",
                         "about %.3g times as many."),
                  n_scored, n, rows_carry, n_support, n_all, ceiling_L, all_ceiling, all_what,
                  n_support, rows_carry, as.integer(min(by_support, k_cap_all)), min_cell_n,
                  per_layer)
        else
          sprintf(paste0("resolution_profile(): %d of the %d points carry %s, so the ceiling ",
                         "is %d cells, not the %d that %s would give: the cell means the ",
                         "criteria score are built from the rows with %s. Where those rows ",
                         "sit in one part of the extent the cells elsewhere hold none, and ",
                         "this ceiling is lower than the occupied cells need; ",
                         "`cell_resp_mean` gives the scored rows per occupied cell at each ",
                         "level%s, and a profile of the rows with %s alone bounds the ladder ",
                         "by the area they cover."),
                  n_support, n_all, rows_carry, ceiling_L, all_ceiling, all_what, rows_carry,
                  if (is.finite(per_layer))
                    sprintf(" (rows of the subsample; the cells of the layer hold about %.3g times as many)",
                            per_layer)
                  else "",
                  rows_carry)
        .warn_classed("spatialkit_response_ceiling", "%s", msg,
                      data = list(n = n_all, n_response = n_support, ceiling = ceiling_L,
                                  ceiling_from = ceiling_from, n_sample = n,
                                  n_scored = as.integer(n_scored)))
      }
    }
    # The verdict is about the data, so it is taken against the layer's
    # ceiling; a subsample too small to reach the floor is a setting to change,
    # and is said separately.
    supported <- !too_few && floor_raw <= max(2, data_ceiling)
    # R warnings with a class, not log lines: `supported = FALSE` changes what
    # every count read off the profile means, and a log line reaches neither
    # tryCatch() nor options(warn = 2).
    if (too_few) {
      # The largest min_cell_n that leaves two cells.  Below one there is none:
      # a single scored row was told to "lower `min_cell_n` (at most 0 here)".
      half <- as.integer(floor(n_support / 2))
      .warn_classed("spatialkit_not_supported",
                    paste0("resolution_profile(): %s, so %s. %s, and `supported` is FALSE. %s"),
                    if (k_cap_all < 2L)
                      sprintf("the %d points sit on %d distinct location%s", n_all, n_uniq_all,
                              if (n_uniq_all == 1L) "" else "s")
                    else sprintf("%d %s at min_cell_n = %d support%s %d cell%s, and it takes two",
                                 n_support,
                                 if (n_support == 1L) sub("^rows", "row", sub("^points$", "point", support_what))
                                 else support_what,
                                 min_cell_n, if (n_support == 1L) "s" else "",
                                 as.integer(by_support), if (by_support == 1) "" else "s"),
                    if (levels_given) "no cell count meets `min_cell_n`"
                    else "there is no ladder of cell counts to choose from",
                    # A ladder the caller set is scored as given, and said to be.
                    if (levels_given)
                      sprintf("The %d level%s asked for in `levels` (%s) %s scored as given",
                              length(levels), if (length(levels) == 1L) "" else "s",
                              if (length(levels) == 1L) sprintf("%d cells", levels)
                              else sprintf("%d to %d cells", min(levels), max(levels)),
                              if (length(levels) == 1L) "is" else "are")
                    else "The profile scores 2 cells so that level is still visible",
                    if (k_cap_all < 2L) "A tessellation needs points at two locations or more."
                    else if (half < 1L)
                      sprintf("Two cells take %s on at least two rows, whatever `min_cell_n` is.",
                              rows_carry)
                    else sprintf("Lower `min_cell_n` (at most %d here), or add %s.", half,
                                 support_what),
                    data = list(reason = ceiling_from, n = n_support, min_cell_n = min_cell_n))
    } else if (!supported)
      .warn_classed("spatialkit_not_supported",
                    paste0("resolution_profile(): cells no wider than the autocorrelation ",
                           "range (%.0f) would need at least %.0f of them, but %d %s ",
                           "at min_cell_n = %d support at most %d. The data cannot ",
                           "support a tessellation that respects their own ",
                           "correlation structure; the profile runs from 2 to %d so ",
                           "the cost of each level is still visible."),
                    range_eff, floor_raw, n_support, support_what, min_cell_n,
                    max(2L, as.integer(data_ceiling)), ceiling_L,
                    data = list(reason = "floor above ceiling", n = n_support,
                                min_cell_n = min_cell_n))
    else if (floor_raw > ceiling_L)
      .log_warn(paste0("resolution_profile(): cells no wider than the autocorrelation ",
                       "range (%.0f) need at least %.0f of them, which the %d %s ",
                       "support, but k-means on the %d-point subsample can fit at ",
                       "most %d; raise `sample_n`. The profile runs from 2 to %d."),
                range_eff, floor_raw, n_support, support_what, n, ceiling_L, ceiling_L)
    # With range_floor = FALSE, the default, the ladder starts at 2 whatever the
    # range, and the floor is only reported.  The range is an estimate.  It is
    # re-estimated on every subset, and a floor that applies on one fold and
    # not on the next starts their ladders at different levels: reliability and
    # the elbow, which sit near the bottom of a ladder, then move between folds
    # with the ladder's start.  And on a response with no structure a
    # spurious range of 40 to 122 km started 11 of 100 white-noise ladders at
    # 11 to 98 cells of 100, where Cp then took the floor or a level beside it.
    # A floor equal to the ceiling is not applied either: it left a ladder of
    # one level (3 cells at 30 points, 5 at 45), "the only level scored".
    lo <- if (range_floor && supported && floor_raw < ceiling_L) floor_L else 2L
    if (range_floor && supported && floor_raw == ceiling_L && ceiling_L > 2L && !levels_given)
      .log_info(paste0("resolution_profile(): the range floor (%d cells) equals the ",
                       "ceiling, which would leave one level and nothing to compare ",
                       "it with; the ladder runs from 2 to %d."), ceiling_L, ceiling_L)
    if (!levels_given) {
      levels <- unique(as.integer(round(exp(seq(log(lo), log(ceiling_L),
                                                length.out = max(2L, as.integer(n_levels)))))))
      # The floor is a rung of the ladder whenever the layer supports it and it
      # lies between the ladder's ends, applied or not.  select_resolution()
      # holds `reliability` to the levels at or above the floor, and without a
      # rung on it the pick was the next rung up: 54 cells for a floor of 45,
      # 81 for 68 at the default 20 levels, and a step of the ladder is a
      # factor of 4.5 at 5 levels.  One level more for the other criteria.
      floor_rung <- if (supported && is.finite(range_eff) && !is.na(floor_L) &&
                        floor_L > lo && floor_L < ceiling_L) floor_L else NULL
      levels <- sort(unique(c(lo, levels, floor_rung, ceiling_L)))
    }
  }
  bounds <- list(floor = floor_L, ceiling = ceiling_L, ceiling_from = ceiling_from,
                 supported = supported, area = area, range = range_eff, n = n_all,
                 n_sample = n, n_distinct = n_uniq_all, min_cell_n = min_cell_n,
                 range_floor = range_floor,
                 n_response = if (has_resp) as.integer(n_resp) else NA_integer_,
                 n_scored = as.integer(n_scored),
                 levels_given = levels_given || size_ladder,
                 predictors = has_pred)
  # The sizes, largest first, and the level each gives: NA for a size that
  # was dropped.
  if (size_ladder) {
    bounds$obs_per_cell <- ladder_sizes
    bounds$size_levels  <- sl$size_levels
  }

  # A ladder of sizes passed its guards before the variogram (above).  On
  # every other ladder the levels are known only now, and the one search to
  # bound is the estimate's, of each scored row's nearest centre.
  if (!size_ladder) centre_guard(max(levels))

  # Over the hull the area is measured on, not its bounding box (see
  # .rbar_domain()); reliability needs an identified range.
  rbar_V <- if (!is.null(vg$cor_fn)) .rbar_domain(vg$cor_fn, hull, area) else NA_real_
  pred_for_moran <- if (has_pred) pred else matrix(numeric(0), nrow = n, ncol = 0L)

  # The sweep.
  m <- length(levels)
  out <- data.frame(levels = levels, wss = NA_real_, wss_spread = NA_real_,
                    elbow = NA_real_, cell_n_min = NA_integer_,
                    cell_n_median = NA_real_, cell_diam_median = NA_real_,
                    rss = NA_real_, cp = NA_real_, moran_i = NA_real_,
                    moran_z = NA_real_, reliability = NA_real_,
                    cells_resp = NA_integer_, cell_resp_mean = NA_real_)
  # One cell for the whole layer: the level every other has to beat.  Its RSS
  # is the total sum of squares of the scored rows, exact and free; its Cp is
  # filled in below, once the noise variance is known.  The ladder has no
  # such level, so on a response with no spatial structure the argmin over
  # it, which is noise, used to be reported with nothing to say that one
  # mean does as well (picks of 2 to 750 cells over 29 white-noise layers).
  one_cell <- NULL
  is_const <- FALSE
  if (!is.null(y) && sum(is.finite(y)) >= 2L) {
    y1 <- y[is.finite(y)]
    one_cell <- list(n = length(y1), rss = sum((y1 - mean(y1))^2), cp = NA_real_)
    # Constant to within rounding: no level fits it better than another, the
    # variogram has nothing to fit and the fallback noise variance is 0.
    # Rounding is measured on the scale of the response: OLS residuals of a
    # response the predictors fit exactly are all rounding error, of size
    # 1e-15, and against their own sum of squares they never looked constant
    # (Cp was then computed with a noise variance of 7.9e-31).
    one_cell$constant <- one_cell$rss <=
      100 * .Machine$double.eps^2 * max(sum(y1^2), sum(resp[is.finite(y)]^2))
    is_const <- isTRUE(one_cell$constant)
    # Nothing is scored on a constant: `cp`, `moran_i` and `moran_z` stay NA
    # (see the sweep), which is what the message says.  A response of
    # 1e6 + 1e-9 x noise was called constant here and then given a finite
    # `cp` computed from its rounding error, and a count of 2.
    if (is_const)
      .warn_classed("spatialkit_constant_response",
                    paste0("resolution_profile(): the %s is constant over the %d ",
                           "scored row(s) (its total sum of squares, %.3g, is ",
                           "rounding error at most), so no cell count fits it better ",
                           "than another: `cp`, `moran_i` and `moran_z` are NA at ",
                           "every level."),
                    if (identical(variable, "residuals")) "OLS residual of the response on `predictor_vars`"
                    else "response", one_cell$n, one_cell$rss,
                    data = list(n = one_cell$n))
  }
  # Moran's z of the cell means needs FNN and Matrix above .KNN_DENSE_MAX
  # cells.  On a ladder of sizes that is 22,500 rows away at 4.5 per cell, so
  # there the step is left out at such a level, with one warning here, before
  # any fit.  Every other ladder does the same while .moran_skip_any_ladder()
  # is TRUE; with it FALSE such a ladder stops with the error of
  # .build_knn_weights().
  # The step is left out by the cells that hold a scored row, which are known
  # only after a level's fit, so the warning says what is known here: where
  # that happens.  It is raised when it can happen: the step runs at all (a
  # response that is not constant), and both a level and the scored rows are
  # above the limit.  Only then are the two backends asked for here.
  moran_off <- (size_ladder || .moran_skip_any_ladder()) && !is.null(y) && !is_const &&
    max(levels) > .KNN_DENSE_MAX && isTRUE(n_scored > .KNN_DENSE_MAX) &&
    !(backend("FNN") && backend("Matrix"))
  if (moran_off) {
    lacking <- c(if (!backend("FNN")) "FNN", if (!backend("Matrix")) "Matrix")
    .warn_classed("spatialkit_moran_skipped",
                  paste0("resolution_profile(): `moran_i` and `moran_z` will be NA at a level ",
                         "where more than %d cells hold a scored row (this ladder reaches %d ",
                         "cells), because %s not ",
                         "installed: Moran's z of the cell means needs FNN to find each ",
                         "cell's neighbours and Matrix to hold them. Every other column is ",
                         "computed."),
                  .KNN_DENSE_MAX, max(levels),
                  if (length(lacking) > 1L) "FNN and Matrix are" else paste(lacking, "is"),
                  data = list(limit = .KNN_DENSE_MAX,
                              levels = as.integer(levels[levels > .KNN_DENSE_MAX]),
                              missing = lacking))
  }
  # The estimate's columns, one value per level, and what every level is
  # compared with: one mean for all the rows `loo` reads.
  est_cols <- NULL
  est_const <- FALSE
  est_one <- c(NA_real_, NA_real_)
  est_cache <- NULL
  if (!is.null(est)) {
    na_m <- rep(NA_real_, m)
    est_cols <- list(scored_per_cell = na_m, loo = na_m, loo_se = na_m, share_alone = na_m,
                     own_weight = na_m, sigma2 = na_m, tau2 = na_m,
                     cells = rep(NA_integer_, m), f_p = na_m)
    n_lay <- length(lay$y)
    # The engine's test of a constant column, and its error of one mean.  It
    # takes two rows to be constant, as it does for `one_cell`: one scored
    # row is too few rows, and the constant-response warning is not raised
    # for it.
    est_const <- n_lay >= 2L && .no_spread(sum((lay$y - mean(lay$y))^2), lay$y)
    if (n_lay >= 2L && !est_const) {
      e1 <- (lay$y - mean(lay$y)) * n_lay / (n_lay - 1)
      est_one <- c(mean(e1^2), stats::sd(e1^2) / sqrt(n_lay))
    }
    # A row alone in its region is predicted, under plain means, from the
    # region of its nearest other row.  The nearest other location of a
    # location is the same at every level, so every level's fit shares one
    # cache: filled here in one search with FNN, and on first need without.
    est_cache <- new.env(parent = emptyenv())
    if (identical(est, "mean") && n_lay >= 2L && !est_const && backend("FNN"))
      .nearest_other(lay$xy, rep.int(1L, n_lay), seq_len(n_lay), backend("FNN"), est_cache,
                     "resolution_profile")
  }
  # The centres of the partition each level was scored on.  A k-means
  # partition is the Voronoi partition of its centres, so these are the cells
  # the criteria judged; get_voronoi_seeds() and voronoi_seeds_kmeans() hand
  # them on (via select_resolution()) instead of running a k-means of their
  # own, whose cells were not the ones scored.
  centres <- vector("list", m)
  names(centres) <- as.character(levels)
  # Rows and occupied cells behind each level's RSS, for Cp's fallback noise
  # estimate below.
  m_ok_v <- rep(NA_real_, m); L_ok_v <- rep(NA_real_, m)
  # What the columns mean under partition = "density".  The cells are then
  # the Voronoi cells of centres that minimise a density-weighted sum, and
  # each column is still read off those cells:
  #   wss, elbow    The within-cell sum of squares of the coordinates, every
  #                 point counted once, about each cell's own mean.  That is
  #                 what the kind "kmeans" reports (its centres are the cell
  #                 means), so the column and the elbow read from it compare
  #                 between kinds, and the curve still starts from the total
  #                 sum of squares at one cell.  The weighted sum the centres
  #                 minimise is not a sum of squares of the layer: it weights
  #                 a point by its density, and its one-cell value is not the
  #                 total sum of squares.  It is kept in attr(, "partition").
  #                 The elbow rule's threshold holds on points with no
  #                 clusters: over 100 uniform layers of 400 to 1,500 points
  #                 the largest sag was 0.006 under "kmeans" and 0.048 under
  #                 "density", against the log(1.25) = 0.22 an elbow takes,
  #                 and on 60 layers of 800 points whose density varies
  #                 without clusters (a gradient, two plateaus, a radial
  #                 fall) an elbow was read on 0 and on 2.  Read off the
  #                 weighted sum it was read on 10 of those 60, which is the
  #                 other reason to read the unweighted one.  Where the
  #                 points do cluster the density kind's curve bends more
  #                 sharply: past the number of clusters its centres go into
  #                 the clusters, and the sparse ground keeps its sum of
  #                 squares.  On 30 layers of 3 to 8 clusters over a sparse
  #                 background an elbow was read on 24 against 10.
  #   wss_spread    The spread, over the starts, of what the starts minimise
  #                 and are ranked by: the weighted sum.
  #   cell_n_*, cell_diam_median, cells_resp, cell_resp_mean
  #                 Counts and sizes of the cells, whatever drew them.
  #   rss, cp       The cells are drawn from the coordinates alone under
  #                 either kind (the weights are a function of the
  #                 coordinates), so the cell means are a linear fit with one
  #                 parameter per occupied cell, which is all Mallows' Cp
  #                 assumes.  Its variance term, tau^2 times the cells over
  #                 the rows, is the sum over the cells of (rows in the cell
  #                 / rows) x tau^2 / (rows in the cell), and holds however
  #                 the rows are shared out.
  #   moran_i, moran_z, the F test against one cell
  #                 Statistics of the cell means under the hypothesis of no
  #                 structure, valid for any cells drawn without the
  #                 response.  Moran's z assumes cell means of equal
  #                 variance, which cells of about equal counts come closer
  #                 to.
  #   reliability   Not read off the cells at all: it is computed for square
  #                 cells of the level's average area that hold the level's
  #                 average count, so it is the same number under either
  #                 kind.  The cells of the density kind match its counts
  #                 more closely and its areas less.
  # Nothing here is undefined under the density kind, so no column is NA for
  # the kind's sake.  A level at which no start gives a partition is NA under
  # either kind.
  objective_v <- rep(NA_real_, m)
  # Iterations the kept start of each level took, and whether it stopped
  # before the limit.
  iter_v <- rep(NA_integer_, m); converged_v <- rep(NA, m)
  # Occupied cells per level.  Every cell of the kind "kmeans" holds a point;
  # a centre of the density kind that loses its points stays where it is,
  # and its cell then adds no mean to Cp's variance term.
  L_fit_v <- as.numeric(levels)
  for (i in seq_len(m)) {
    L  <- levels[i]
    kb <- .partition_centres(xy, L, partition = part$kind, init = part$init,
                             nstart = nstart, iter.max = part_iter,
                             density_power = density_power, density_k = density_k,
                             use_fnn = backend("FNN"))
    if (is.null(kb)) {
      .log_warn("resolution_profile(): %s failed at %d cells; that level is NA.",
                if (identical(part$kind, "density")) "the density-weighted k-means"
                else "k-means", L)
      next
    }
    cl <- kb$cluster
    centres[[i]] <- kb$centers
    out$wss[i]        <- kb$wss
    out$wss_spread[i] <- kb$spread
    objective_v[i]    <- kb$objective
    iter_v[i]         <- as.integer(kb$iter)
    converged_v[i]    <- isTRUE(kb$converged)
    sizes <- tabulate(cl, nbins = L)
    L_fit <- sum(sizes > 0L)
    L_fit_v[i] <- L_fit
    if (L_fit < L)
      .warn_classed("spatialkit_empty_cells",
                    paste0("resolution_profile(): at %d cells, %d of the centres hold no ",
                           "point (a centre that loses its points stays where it is), so ",
                           "the level has %d occupied cells: `cell_n_min` is 0 there, and ",
                           "`cp` counts the occupied cells."),
                    L, L - L_fit, L_fit,
                    data = list(level = as.integer(L), empty = as.integer(L - L_fit)))
    out$cell_n_min[i]    <- as.integer(min(sizes))
    out$cell_n_median[i] <- stats::median(sizes)
    rad <- sqrt(kb$withinss / pmax(sizes, 1L))
    out$cell_diam_median[i] <- 2 * stats::median(rad[sizes >= 2L])
    if (!is.null(y)) {
      okr <- is.finite(y)
      m_ok <- sum(okr)
      L_ok <- length(unique(cl[okr]))             # cells holding a scored row
      out$cells_resp[i] <- L_ok
      if (L_ok > 0L) out$cell_resp_mean[i] <- m_ok / L_ok
      if (m_ok > L_ok) {
        cm  <- stats::ave(y[okr], cl[okr])
        out$rss[i] <- sum((y[okr] - cm)^2)
        m_ok_v[i] <- m_ok; L_ok_v[i] <- L_ok
        # Cp for the cells of the whole layer, estimated from the subsample.
        # RSS / m + tau^2 L_ok / m estimates the approximation error plus
        # tau^2 (it adds back the optimism of the subsample's own cell means,
        # tau^2 per cell); tau^2 L / N is the sampling variance of cell means
        # built from the layer's N rows with a response.  With no subsample
        # (m = N, L_ok = L) it is Mallows' RSS / n + 2 tau^2 L / n.
        if (!is.null(vg) && !is_const)
          out$cp[i] <- out$rss[i] / m_ok + vg$nugget * (L_ok / m_ok + L_fit / n_resp)
      }
      # Moran's I of cell means that differ by rounding error alone is
      # rounding error too.
      # Left NA as well where the neighbour weights of more than
      # .KNN_DENSE_MAX cell means cannot be built (`moran_off`, above).  On a
      # ladder of sizes the cell means come from grouped sums: the loop makes
      # a pass over all the rows for every cell, which is quadratic in the
      # rows at a fixed size.
      if (!is_const && !(moran_off && L_ok > .KNN_DENSE_MAX)) {
        mi <- .morans_i_for_k(xy[okr, , drop = FALSE], y[okr],
                              pred_for_moran[okr, , drop = FALSE], cl[okr],
                              grouped = size_ladder, use_fnn = backend("FNN"),
                              use_matrix = backend("Matrix"))
        out$moran_i[i] <- mi[["I"]]
        out$moran_z[i] <- mi[["z"]]
      }
    }
    # The layer's point count: the cell means are built from every point,
    # not from the subsample the partition was fitted on.
    if (!is.null(vg$cor_fn))
      out$reliability[i] <- .reliability_at(L, area, n_resp, vg$nugget, vg$psill,
                                            vg$cor_fn, rbar_V)
    # The estimate, last: on every scored row of the layer, each in the cell
    # of its nearest centre.  It draws no random number and changes nothing
    # above, so the columns a profile has without an estimate are the same
    # numbers with one.
    if (!is.null(est) && length(lay$y)) {
      lv <- .estimate_at_level(lay$xy, lay$y, kb$centers, est, backend("FNN"), est_cache)
      for (nm in names(est_cols)) est_cols[[nm]][i] <- lv[[nm]]
    }
  }
  # Each row is left out alone, and the other rows at its location stay.
  # Said once, after the levels are scored and where one of them has a
  # `loo`: on a constant response, and where no level could be scored, the
  # sentence would describe a number the profile does not hold.
  if (!is.null(est) && any(is.finite(est_cols$loo))) {
    repeated <- .repeated_rows(lay$xy)
    if (repeated > 0L)
      .warn_classed("spatialkit_repeated_locations",
                    paste0("resolution_profile(): %d of the %d rows with a value of '%s' ",
                           "share a location with another row. The leave-one-out error ",
                           "leaves one row out and keeps the others at its location, so it, ",
                           "and an interval built from it, describes a further observation ",
                           "at a place already observed, not an observation at a new place."),
                    repeated, length(lay$y), response_var,
                    data = list(n = as.integer(repeated), total = length(lay$y),
                                column = response_var))
  }

  # Cp's noise variance.  From the variogram's nugget when it gives one; when
  # it gives none (no variogram, a fit that did not converge, a range refused
  # as below the shortest lag, no model at all) Cp used to be NA at every
  # level, and a benchmark then scored whatever rule stood in for it, one
  # training set in six on real data.  The fallback is Mallows' own choice:
  # the residual mean square of the most flexible model, here the finest
  # level whose cells hold at least two scored rows on average.  It counts
  # the structure within those cells as noise as well, so on average it is
  # no smaller than the nugget and errs towards fewer cells.
  cp_noise <- if (!is.null(vg) && !is_const)
    list(source = "variogram nugget", value = vg$nugget, level = NA_integer_)
  else NULL
  # A ladder of sizes under sac = "auto" fits no variogram, and `cp` gets no
  # noise variance from the finest level in its place: on such a ladder that
  # level is one the caller chose, so the number would be set by the ladder
  # asked for.  The reason is recorded where the noise variance would be.
  if (no_vg)
    cp_noise <- list(source = "none", value = NA_real_, level = NA_integer_,
                     reason = paste0("no variogram: a size ladder fits one only when ",
                                     "`sac` is given"))
  if (is.null(vg) && !is.null(y) && !is_const && !no_vg) {
    dfree <- m_ok_v - L_ok_v
    cand  <- which(is.finite(out$rss) & is.finite(dfree) & dfree >= L_ok_v)
    if (length(cand)) {
      j  <- cand[which.max(levels[cand])]
      s2 <- out$rss[j] / dfree[j]
      if (is.finite(s2) && s2 > 0) {
        fin_cp <- is.finite(out$rss) & is.finite(m_ok_v) & m_ok_v > 0
        out$cp[fin_cp] <- out$rss[fin_cp] / m_ok_v[fin_cp] +
          s2 * (L_ok_v[fin_cp] / m_ok_v[fin_cp] + L_fit_v[fin_cp] / n_resp)
        cp_noise <- list(source = "finest-level residual mean square", value = s2,
                         level = as.integer(levels[j]))
        # A variogram that found no spatial structure does carry a model, the
        # pure nugget summarize_by_cell() and kriging_adequacy() use, so "no
        # usable variogram model" was untrue of it and contradicted the line
        # logged above.  What it lacks is a nugget apart from the sill.
        why <- if (is.null(sac))
          "no variogram was fitted (gstat is not installed, or estimate_sac_range() failed; see the log)"
        else if (!is.data.frame(attr(sac, "variogram_model")))
          sprintf("%s carries no variogram model and so no nugget", sac_what)
        else if (.sac_no_structure(sac))
          sprintf("%s found no spatial structure and gives no nugget that can be told from the sill",
                  sac_what)
        else sprintf("%s gives no usable variogram model", sac_what)
        .warn_and_log(paste0("resolution_profile(): %s, so `cp` takes its noise ",
                             "variance from the residual mean square of the finest ",
                             "level (%d cells: %.3g), Mallows' own choice. It counts ",
                             "the structure within those cells as noise too, so it ",
                             "errs towards fewer cells; attr(<profile>, \"cp_noise\") ",
                             "records it. `reliability` needs the variogram's range ",
                             "and is NA."),
                      why, levels[j], s2)
      }
    }
  }

  # Cp at one cell, by the formula of the sweep (L = L_ok = 1) and with the
  # noise variance the levels were given.  `cp_levels` ties it to the column
  # it belongs to: select_resolution() compares the two only while that
  # column is the one computed here.
  if (!is.null(one_cell)) {
    if (!is.null(cp_noise) && is.finite(cp_noise$value))
      one_cell$cp <- one_cell$rss / one_cell$n +
        cp_noise$value * (1 / one_cell$n + 1 / n_resp)
    one_cell$cp_levels <- stats::setNames(out$cp, levels)
    # Whether a level's cells differ from one mean by more than chance: the
    # F test of the cell means, (TSS - RSS) / (cells - 1) over RSS / (rows -
    # cells), on the cells that hold a scored row.  The cells are fitted to
    # the coordinates alone, so with no spatial structure in the scored
    # variable the statistic has its F distribution at every level, and the
    # test needs no noise variance.  select_resolution() reads it at the Cp
    # optimum: a comparison of Cp at one cell with Cp there has no sampling
    # error in it, and none at all to work with when Cp's noise variance is
    # 0 (Cp is then RSS / n, which falls with every finer level).
    df1 <- L_ok_v - 1; df2 <- m_ok_v - L_ok_v
    f_p <- rep(NA_real_, m)
    f_ok <- !is_const & is.finite(out$rss) & is.finite(df1) & df1 >= 1 & is.finite(df2) & df2 >= 1
    if (any(f_ok)) {
      f_stat <- ((one_cell$rss - out$rss[f_ok]) / df1[f_ok]) / (out$rss[f_ok] / df2[f_ok])
      f_p[f_ok] <- stats::pf(pmax(f_stat, 0), df1[f_ok], df2[f_ok], lower.tail = FALSE)
    }
    one_cell$f_p <- stats::setNames(f_p, levels)
  }
  attr(out, "one_cell") <- one_cell

  # The elbow over the ladder, and the bumps on it.  The sag of log WSS below
  # the straight line from k = 1 to the last level on log-log axes (see
  # .elbow_sag()); NA at every level when the curve has no elbow.  On linear
  # axes the chord from the first level to the last always has a point
  # furthest below it, at about sqrt(first x last level) on a layer with no
  # cluster structure, so min_cell_n and the range floor chose the "elbow",
  # and a geometry-only profile handed that count to build_tessellation().
  # Anchored at k = 1, whose WSS is the total sum of squares and exact,
  # rather than at the first level: the ladder often starts at the range
  # floor, and a uniform square's WSS sits above c / k at k = 2 (two halves
  # keep 5/8 of the total, not 1/2), so a chord from 2 finds a "knee" at its
  # four quadrants.
  fin <- is.finite(out$wss)
  # A level with a WSS of 0 (to within 1e-12 of the total) has a cell on
  # every distinct location.  It is left out of the line, which log(0) would
  # take with it, and when the rest of the curve has no elbow it is the
  # elbow, the fall to zero being the sharpest bend there is (see
  # .elbow_read()).  Leaving it out alone made five stations visited thirty
  # times each a layer with no elbow.
  if (any(fin)) {
    tss <- sum(sweep(xy, 2L, colMeans(xy))^2)
    rd  <- .elbow_read(c(1L, levels[fin]), c(tss, out$wss[fin]))
    if (rd$structured)
      out$elbow[fin] <- rd$sag[-1L]
    else if (sum(!rd$zero[-1L]) >= 2L)
      .log_info(paste0("resolution_profile(): the WSS curve has no elbow (on log-log ",
                       "axes it falls in a straight line, as it does for points with ",
                       "no cluster structure); `elbow` is NA."))
  }
  wss_bumps <- .wss_bumps(out$wss[fin])
  if (wss_bumps > 0L && part_default)
    .log_warn(paste0("resolution_profile(): the WSS curve rises at %d step(s) of ",
                     "the ladder even with %d k-means++ restarts; the elbow and ",
                     "the neighbouring levels' spread are partly optimisation noise."),
              wss_bumps, nstart)
  else if (wss_bumps > 0L && identical(part$kind, "kmeans"))
    .log_warn(paste0("resolution_profile(): the WSS curve rises at %d step(s) of ",
                     "the ladder even with %s per level; the elbow and the ",
                     "neighbouring levels' spread are partly optimisation noise."),
              wss_bumps, .partition_label(part, nstart))
  else if (wss_bumps > 0L)
    # The density kind minimises the weighted sum, so a rise of the
    # unweighted one is not by itself a start that went wrong; a rise of the
    # weighted sum is.
    .log_warn(paste0("resolution_profile(): the WSS curve rises at %d step(s) of ",
                     "the ladder. The centres of partition = \"density\" minimise the ",
                     "density-weighted sum (attr(<profile>, \"partition\")$objective), ",
                     "which rises at %d step(s) with %s per level; the elbow is read ",
                     "off the unweighted curve and is that much less sure."),
              wss_bumps, .wss_bumps(objective_v[is.finite(objective_v)]),
              .partition_label(part, nstart))

  # A level whose kept start was cut off at the iteration limit: its centres
  # are those of the last step, and a longer run can still move points
  # between cells.  Said for every call but the one this function has always
  # made (k-means from k-means++ starts with its own limit), whose result and
  # conditions are unchanged.
  cut_off <- which(converged_v %in% FALSE)
  if (length(cut_off) && (!part_default || !is.null(iter_max)))
    .warn_not_converged("resolution_profile", levels[cut_off], iter_v[cut_off], part_iter,
                        part$kind, "iter_max", of = m)

  # The record of the kind that built the regions.  Absent for the present
  # kind, so that a profile made with it is unchanged; .partition_record()
  # reads either.  The weighted sum each level's centres minimise goes with
  # it under the density kind, since `wss` there is the unweighted one, and
  # so do the iterations each level's kept start took and whether it stopped
  # before the limit.
  part_attr <- NULL
  if (!part_default) {
    part_attr <- part[.partition_fields()]
    if (identical(part$kind, "density"))
      part_attr$objective <- stats::setNames(objective_v, levels)
    part_attr$iter_max  <- as.integer(part_iter)
    part_attr$iter      <- stats::setNames(iter_v, levels)
    part_attr$converged <- stats::setNames(converged_v, levels)
  }

  # The estimate's five columns, after the columns a profile has without one.
  if (!is.null(est))
    for (nm in c("scored_per_cell", "loo", "loo_se", "share_alone", "own_weight"))
      out[[nm]] <- est_cols[[nm]]

  out <- structure(out,
            class      = c("resolution_profile", "data.frame"),
            bounds     = bounds,
            variogram  = if (is.null(vg)) NULL
                         else vg[c("nugget", "psill", "range", "model", "detrended")],
            variable   = variable,
            wss_bumps  = wss_bumps,
            nstart     = nstart,
            sac        = sac,
            split      = split,
            cp_noise   = cp_noise,
            centres    = centres,
            centres_crs = sf::st_crs(data_sf),
            # The convex hull of every point profiled, in that CRS: the ground
            # the scored centres stand for.  get_voronoi_seeds() hands it on as
            # the seeds' cover when it is given the profile and neither points
            # nor a boundary, so their Voronoi cells reach the points' hull and
            # not only the centres'.
            hull       = hull)
  if (!is.null(part_attr)) attr(out, "partition") <- part_attr
  # What the estimate's columns were computed for and against.
  if (!is.null(est)) {
    attr(out, "estimate") <- list(
      value        = est,
      neighbours   = if (identical(est, "pooled")) .POOL_NEIGHBOURS else NA_integer_,
      response     = response_var,
      n            = length(lay$y),
      rows         = if (is.null(split)) "all" else "selection half",
      features     = if (was_points) "points" else "other",
      seed         = seed,
      constant     = est_const,
      loo_one_cell = est_one[1L], loo_one_cell_se = est_one[2L],
      loo_levels   = stats::setNames(est_cols$loo, levels),
      sigma2       = stats::setNames(est_cols$sigma2, levels),
      tau2         = stats::setNames(est_cols$tau2, levels),
      cells        = stats::setNames(est_cols$cells, levels),
      f_p          = stats::setNames(est_cols$f_p, levels))
  }
  out
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
  # The rows the support ceiling counts, when they are fewer than the points.
  n_rsp <- b$n_response %||% NA_integer_
  some_missing <- isTRUE(n_rsp > 0L && n_rsp < b$n)
  too_few <- b$ceiling_from %in% c("fewer than two supported cells",
                                   "fewer than two distinct locations")
  # What a row carries to be counted: with predictors, a predictor missing
  # leaves it out as much as a response missing does.
  carry <- if (isTRUE(b$predictors)) "with a response and every predictor" else "with a response"
  cat("Resolution profile:", nrow(x), "levels on", b$n,
      if (isTRUE(n_fit < b$n) && some_missing)
        sprintf("points (%d %s; k-means fitted to a subsample of %d)\n", n_rsp, carry, n_fit)
      else if (isTRUE(n_fit < b$n)) sprintf("points (k-means fitted to a subsample of %d)\n", n_fit)
      else if (some_missing) sprintf("points (%d %s)\n", n_rsp, carry)
      else "points\n")
  # A ladder `obs_per_cell` set: its levels are the counts the sizes give,
  # so the line names the sizes and the rows they divide, and no ceiling.
  if (!is.null(b$obs_per_cell)) {
    n_div <- if (isTRUE(n_rsp > 0L)) n_rsp else b$n
    # The sizes that have a level: a size of "auto" that gave fewer than two
    # cells was left out of the ladder, and is left out of its line.
    built <- if (length(b$size_levels) == length(b$obs_per_cell))
      b$obs_per_cell[!is.na(b$size_levels)] else b$obs_per_cell
    cat(sprintf("  ladder      : %s cells (obs_per_cell = %s on the %d %s; centres fitted to every row; floor %s)%s\n",
                if (min(x$levels) == max(x$levels)) sprintf("%d", max(x$levels))
                else sprintf("%d to %d", min(x$levels), max(x$levels)),
                paste(vapply(built, function(s) format(s), ""), collapse = ", "),
                n_div,
                if (isTRUE(n_rsp > 0L)) paste(if (n_div == 1L) "row" else "rows", carry)
                else "points",
                if (!is.finite(b$range)) "2 (no range)"
                else paste0(if (is.na(b$floor)) sprintf("beyond integer range from range %.0f", b$range)
                            else sprintf("%d from range %.0f", b$floor, b$range),
                            ", not applied"),
                if (!isTRUE(b$supported) && identical(b$ceiling_from, "fewer than two distinct locations"))
                  "  -- the points sit on one location: not supported"
                else if (!isTRUE(b$supported) && identical(b$ceiling_from, "fewer than two cells of the size"))
                  "  -- fewer than two cells of the size: not supported"
                else if (!isTRUE(b$supported)) "  -- floor above the finest level: not supported"
                else ""))
  } else {
    cat(sprintf("  ladder      : %d to %d cells (floor %s, ceiling %d%s)%s\n",
                min(x$levels), max(x$levels),
                if (!is.finite(b$range)) "2 (no range)"
                else paste0(if (is.na(b$floor)) sprintf("beyond integer range from range %.0f", b$range)
                            else sprintf("%d from range %.0f", b$floor, b$range),
                            if (isFALSE(b$range_floor)) ", not applied" else ""),
                b$ceiling,
                # With no location repeated the cap is one short of the points
                # (stats::kmeans() needs fewer centres than points), and "from
                # 30 distinct locations" named a number that did not bind.
                if (identical(b$ceiling_from, "distinct locations") &&
                    isTRUE(b$ceiling < (b$n_distinct %||% NA_integer_)))
                  sprintf(", one short of the %d points", b$n)
                else if (identical(b$ceiling_from, "distinct locations"))
                  sprintf(" from %d distinct locations", b$n_distinct %||% NA_integer_)
                else if (identical(b$ceiling_from, "sample_n") &&
                         isTRUE((b$n_scored %||% NA_integer_) < n_fit))
                  sprintf(" from the %d scored rows of the %d-point subsample", b$n_scored, n_fit)
                else if (identical(b$ceiling_from, "sample_n"))
                  sprintf(" from the %d-point subsample", n_fit)
                else if (identical(b$ceiling_from, "min_cell_n (rows with a response)"))
                  sprintf(" from min_cell_n = %d on the %d rows %s",
                          b$min_cell_n, n_rsp, carry)
                else if (identical(b$ceiling_from, "fewer than two supported cells"))
                  sprintf(": fewer than two cells meet min_cell_n = %d", b$min_cell_n)
                else if (identical(b$ceiling_from, "fewer than two distinct locations"))
                  ": the points sit on one location"
                else sprintf(" from min_cell_n = %d", b$min_cell_n),
                # "no ladder" is said of a ladder the profile built, one level
                # of 2 cells; levels the caller asked for are in the table.
                if (!isTRUE(b$supported) && isTRUE(too_few) && isTRUE(b$levels_given))
                  "  -- levels as given: not supported"
                else if (!isTRUE(b$supported) && isTRUE(too_few)) "  -- no ladder: not supported"
                else if (!isTRUE(b$supported)) "  -- floor above ceiling: not supported"
                else if (isTRUE(b$floor > b$ceiling)) "  -- floor above ceiling: raise sample_n"
                else ""))
  }
  vg <- attr(x, "variogram", exact = TRUE)
  cn <- attr(x, "cp_noise", exact = TRUE)
  cat(sprintf("  variogram   : %s\n",
              if (is.null(vg))
                (if (identical(cn$source, "finest-level residual mean square"))
                   sprintf(paste0("none usable (reliability is NA; cp's noise ",
                                  "variance %.3g is the residual mean square at %d cells)"),
                           cn$value, cn$level)
                 # A ladder of sizes that was given no `sac` fits none.
                 else if (identical(cn$source, "none"))
                   paste0("none fitted: a size ladder fits one only when `sac` is given ",
                          "(cp and reliability are NA)")
                 else "none usable (cp and reliability are NA)")
              else
                sprintf("nugget %.3g, partial sill %.3g, range %s", vg$nugget, vg$psill,
                        if (is.finite(vg$range)) sprintf("%.0f", vg$range) else "unidentified")))
  # The kind that built the cells: unrecorded for k-means from k-means++
  # starts, whose line is the one it always was.
  part <- attr(x, "partition", exact = TRUE)
  cat(sprintf("  scored on   : %s; %s per level; WSS rises at %d step(s)\n",
              if (is.na(attr(x, "variable"))) "geometry only" else attr(x, "variable"),
              if (is.list(part) && !is.null(part$kind)) .partition_label(part, attr(x, "nstart"))
              else sprintf("%d k-means++ restarts", attr(x, "nstart")),
              attr(x, "wss_bumps")))
  # What one cell for the whole layer scores: the ladder starts at 2, so the
  # level every other has to beat is not a row of the table.
  one <- attr(x, "one_cell", exact = TRUE)
  if (is.list(one) && is.finite(one$rss %||% NA_real_))
    cat(sprintf("  one cell    : rss %s%s (one mean for the layer; a level has to beat it)\n",
                format(signif(one$rss, digits)),
                if (isTRUE(one$constant)) ", a constant"
                else if (is.finite(one$cp %||% NA_real_))
                  sprintf(", cp %s", format(signif(one$cp, digits)))
                else ""))
  # The estimate the profile scores, the rows its leave-one-out error reads,
  # and that error with one mean for all of them.
  est <- attr(x, "estimate", exact = TRUE)
  if (is.list(est) && !is.null(est$value))
    cat(sprintf("  estimate    : %s; loo on %d row%s of the %s; one cell: loo %s\n",
                if (identical(est$value, "pooled"))
                  sprintf("pooled toward the %d nearest cells", est$neighbours)
                else "plain means",
                est$n, if (isTRUE(est$n == 1L)) "" else "s",
                if (identical(est$rows, "selection half")) "selection half" else "layer",
                if (isTRUE(est$constant)) "NA (the response is constant)"
                else if (isTRUE(est$n < 2L)) "NA (fewer than two rows)"
                else format(signif(est$loo_one_cell, digits))))
  if (all(c("elbow", "wss") %in% names(x)) && !any(is.finite(x$elbow)) &&
      sum(is.finite(x$wss)) >= 2L)
    cat("  elbow       : none; the WSS curve falls as it does with no cluster structure\n")
  sp <- attr(x, "split")
  if (!is.null(sp))
    cat(sprintf("  split       : response read on %d points; estimate on the other %d (attr \"split\")\n",
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
#' every level whose criterion value is within \code{tol} of the optimum.
#' On the default profile, one size scored for pooled values, the default
#' criterion is that size: it returns the count the size gives, with no
#' optimum and no flat region (see "A profile built with an estimate"); what
#' follows is about the criteria read on a ladder of counts.  On
#' the criteria this package computes the flat region is routinely wide: the
#' reliability curve is flat to within 2 percent over a factor of 3--6 in the
#' number of cells, and \eqn{C_p} on a smooth field descends to the support
#' ceiling.  The region is the answer, and the argmin only a point in it.
#' When the optimum sits at an end of the levels the criterion was scored at,
#' the result says so and names the bound, because a bound is then doing the
#' choosing rather than the criterion (see
#' \code{\link{resolution_profile}} for what each criterion measures and
#' what it returns in practice).  The functions that take the result as a
#' cell count raise that as a warning of class
#' \code{spatialkit_bound_chosen}.
#'
#' @section What each criterion returns:
#' \describe{
#'   \item{\code{cp}}{The level that minimises Mallows' \eqn{C_p}: the count
#'     at which re-measuring the sampled sites, each predicted by its cell's
#'     mean, has the least estimated error (\code{\link{resolution_profile}}
#'     says what that leaves out for a new location).  On a smooth response
#'     with little noise that is the support ceiling.  The result also says whether the optimum beats one
#'     cell for the whole layer (\code{attr(profile, "one_cell")}), in
#'     \code{structured}.  It is \code{TRUE} when two things hold:
#'     \eqn{C_p} at one cell is above the optimum by more than \code{tol},
#'     and the cell means at the optimum differ from one mean by more than
#'     chance allows.  The second is the F test of the cells against one
#'     mean, \eqn{(TSS - RSS) / (L - 1)} over \eqn{RSS / (m - L)}, at 5
#'     percent after a Bonferroni adjustment for the levels the optimum was
#'     searched over (\code{structure_p} is the adjusted p-value).  The
#'     cells are fitted to the coordinates alone, so the test holds its
#'     level on a response with no spatial structure whatever \eqn{C_p}'s
#'     noise variance is, and it is what decides when that variance is 0:
#'     \eqn{C_p} is then \eqn{RSS / n}, which falls with every finer level
#'     by chance alone, and the first comparison says nothing.
#'     When either fails, \code{structured} is \code{FALSE}: one mean for
#'     the whole layer predicts as well as any level scored, so the
#'     response shows no spatial structure the cells can use and the count
#'     is not identified.  \code{structured} is also \code{FALSE} when the
#'     profile's variogram was refused for showing no spatial structure.
#'     \code{best} is still the argmin.
#'     \strong{Measured on white noise} (uniform points, the automatic
#'     ladder of \code{obs_per_cell = NULL}), \code{structured} was
#'     \code{TRUE} on 3 of 200 layers of 60 points, 7 of 200 of 120, 5 of
#'     150 of 300 and none of 80 of 600 or of
#'     40 of 1,500 with the noise variance taken from the finest level; and
#'     with the package's own variogram, whose nugget was 0 on 27 of the
#'     180 layers, on none of 120 layers of 120 points and 3 of 60 of 300.
#'     The first comparison alone gave 27, 20, 8, 0 and 0, and then 38 and
#'     11.  The F distribution is exact for a Gaussian response, and the
#'     test held its level on others: with no structure in a lognormal
#'     response, a 0/1 response with mean 0.1 and a \eqn{t_3} response,
#'     \code{structured} was \code{TRUE} on 2 to 10 of 200 layers at each
#'     of 60, 120 and 300 points.  A \code{FALSE} on a small layer is weak
#'     evidence that there is
#'     no structure: on exponential fields (range parameter 100 on a
#'     1000-unit extent, nugget 0.3 on a unit sill) \code{structured} was
#'     \code{TRUE} on 12 of 40 layers of 60 points, 36 of 40 of 120 and 40
#'     of 40 of 300.}
#'   \item{\code{reliability}}{The coarsest count the range allows.  The
#'     reliability curve peaks below the range floor, so it is read at the
#'     levels at or above the floor the profile reports (when the layer
#'     supports the floor and the profile built its own ladder).  The floor
#'     is a level of every such ladder, so the result is the floor, with
#'     \code{edge} naming it (on a profile whose k-means failed at that
#'     level, the first level above it).  \code{range_floor = FALSE} reads
#'     every level and returns the curve's own maximum, which is usually a
#'     count below the floor.}
#'   \item{\code{elbow}}{The number of clusters of the points, when they
#'     cluster; an error otherwise, because the column is \code{NA} at every
#'     level unless the WSS curve bends.}
#'   \item{\code{moran_z}}{Refused unless \code{diagnostic = TRUE}.  It is
#'     the level at which neighbouring cell means (or, with predictors,
#'     their regression residuals) are least autocorrelated.  \eqn{|z|} is
#'     near 0 at 10 to 20 cells whatever the data and grows as cells get
#'     finer than the range, so its minimum is usually the first level it
#'     can be computed at, or one near it, and it says nothing about how
#'     well the cells represent the response.}
#'   \item{\code{loo}}{The level with the lowest leave-one-out error of the
#'     estimate the profile was built for
#'     (\code{\link{resolution_profile}(estimate = )}); of levels with the
#'     same error, the one with the fewest cells.  The flat region is every
#'     level whose error is within \code{tol} of the lowest.  The column
#'     needs no variogram.  \code{structured} is read as for \code{cp}, with
#'     the leave-one-out error in place of \eqn{C_p}: it is \code{TRUE} when
#'     the error of one mean for all scored rows is above the optimum by
#'     more than \code{tol} and the F test of the cells at the optimum
#'     against one mean, adjusted for the levels searched, has
#'     \eqn{p < 0.05}; \code{one_cell} then holds the error of one mean.
#'     Under this criterion the count depends on the values of the response.
#'     When the profile's estimate is the pooled value and the lowest error
#'     is at the finest level scored, a note says that the end of the ladder
#'     chose the count.  For pooled values the error tends to keep falling to
#'     the finest level: in the simulation described under "Observations per
#'     cell" in \code{\link{resolution_profile}}, the lowest error of the
#'     pooled value was at the finest of seven sizes on 176 of 200 layers,
#'     and that of the plain mean on none.  On a ladder \code{obs_per_cell}
#'     set, an optimum at
#'     an end has the \code{edge} \code{"the finest size offered
#'     (obs_per_cell)"} or \code{"the coarsest size offered
#'     (obs_per_cell)"}; a finest level the profile held to the most cells
#'     k-means can place keeps the label of the support ceiling.}
#'   \item{\code{size}}{The count a size gives,
#'     \code{max(2, round(n / obs_per_cell))}, with \code{n} the rows the
#'     size divides: the rows with a response, or every point for a profile
#'     without one.  No optimum is searched for.  The level returned is the
#'     first of three that applies: the level the profile built for that
#'     size; the count itself, when the profile has a partition at it;
#'     otherwise the built level nearest to the count on the log scale (of
#'     two as near, the one with fewer cells), with a warning of class
#'     \code{spatialkit_size_not_reached}, which says so when the count is
#'     more than k-means can place on the layer.  It reads any profile, one
#'     built without a response or without an estimate included, but not a
#'     column subset, which has lost the row count.  \code{flat} is the level,
#'     \code{value} the rows per cell it leaves (\code{n / best}), and
#'     \code{at_ceiling}, \code{at_floor} and \code{edge} are \code{NA}: no
#'     bound chose anything, so the functions that take the selection as a
#'     count raise no \code{spatialkit_bound_chosen}.  Under this rule
#'     neither the cells, which are fitted to the coordinates, nor their
#'     number, which is a count of rows, depend on the values of the
#'     response.
#'     A size does not ask whether the response has structure the cells can
#'     use, so \code{structured} asks it separately, on a profile built with
#'     an \code{estimate}.  It is \code{FALSE} when the leave-one-out error
#'     of one mean for all scored rows is within \code{tol} of the error at
#'     the level, or below it, or when the F test of the cells at the level
#'     against one mean has \eqn{p \ge 0.05}.  The count stands either way,
#'     since the caller set it: \code{structured = FALSE} says that the cell
#'     values will say little that one mean for the layer does not.  The
#'     functions that take the selection as a count then raise
#'     \code{spatialkit_no_structure}.  On a profile without an estimate
#'     \code{structured} is \code{NA}.}
#'   \item{\code{auto}}{The rule that goes with the profile's estimate:
#'     \code{"size"} for a profile built with \code{estimate = "pooled"},
#'     \code{"loo"} for \code{estimate = "mean"}, and \code{"cp"} for a
#'     profile without an estimate.  A column subset of a profile has lost
#'     the record of its estimate and is read by \code{"cp"}.
#'     \code{criterion} in the result is the criterion applied, never
#'     \code{"auto"}.}
#' }
#'
#' @section A profile built with an estimate:
#' This function's default criterion is \code{"auto"}: \code{"size"} for a
#' profile of pooled values, \code{"loo"} for one of plain means and
#' \code{"cp"} for a profile without an estimate.  A profile built with
#' \code{estimate} that is passed straight to
#' \code{\link{get_voronoi_seeds}()}, \code{\link{voronoi_seeds_kmeans}()} or
#' \code{\link{build_tessellation}()} as the count is read the same way, so
#' \code{select_resolution(profile)} and the profile itself give the same
#' count.  A profile of pooled values built with a ladder argument
#' (\code{n_levels} and the others) is read by the size as well: the
#' default returns the level nearest the count the size gives, on the
#' automatic ladder its top level, with \code{spatialkit_size_not_reached}.
#' The default was \code{"cp"} in earlier versions, whatever the profile;
#' \code{criterion = "cp"} reads it as before on a profile that has it.  \code{cp} is computed on such a profile as on any other and
#' scores plain means; when the profile's estimate is the pooled value the
#' result's note says so.
#'
#' \code{"cp"} and \code{"reliability"} cannot be read off a profile whose
#' ladder \code{obs_per_cell} set under \code{sac = "auto"}.  Such a profile
#' fits no variogram, both columns are \code{NA} at every level, and the
#' call is an error that names the ways to get them: build the profile with
#' \code{sac = NULL}, with a \code{sac} of your own, or with
#' \code{obs_per_cell = NULL}.  The default profile is such a profile.  For
#' the count \code{"cp"} gave in earlier versions, build the profile with
#' \code{estimate = NULL, obs_per_cell = NULL}, and with \code{partition =
#' "kmeans"} for the same regions.  On a ladder of one size that has a
#' variogram either criterion returns the one level, with \code{edge}
#' \code{"the only level scored"}.
#'
#' For a ladder of one size, \code{summary()} of the profile lists each
#' criterion that is finite at the level and says that the size set the
#' count, and \code{plot()} draws each as one point, with no selected level
#' and no band.  When none is finite, as on a constant response, where
#' \code{loo} is \code{NA} too, \code{summary()} stops with an error that
#' says so.
#'
#' @param profile A \code{resolution_profile}.
#' @param criterion Which rule sets the count: \code{"cp"} (minimised),
#'   \code{"reliability"} (maximised at or above the range floor),
#'   \code{"elbow"} (maximised), \code{"loo"} (the leave-one-out error of
#'   the profile's estimate, minimised), \code{"size"} (the count
#'   \code{obs_per_cell} gives), \code{"auto"} (\code{"size"} for a profile
#'   of pooled values, \code{"loo"} for one of plain means, \code{"cp"} for
#'   a profile without an estimate) or, with \code{diagnostic = TRUE},
#'   \code{"moran_z"} (\eqn{|z|} minimised).  Default \code{"auto"}.
#'   \code{"cp"} was the default of earlier versions; name it to read
#'   \eqn{C_p} on a profile that has it.
#' @param tol Width of the flat region.  For \code{cp}, \code{reliability}
#'   and \code{loo} it is relative to the optimum's value (\code{0.02} keeps
#'   levels within 2 percent of it); for \code{elbow} and \code{moran_z},
#'   whose optimum can be zero, it is relative to the criterion's range over
#'   the ladder.  Under \code{"size"} there is no flat region, and
#'   \code{tol} is the margin by which the cells have to beat one mean for
#'   \code{structured} and nothing else.
#' @param range_floor For \code{criterion = "reliability"} only.
#'   \code{TRUE} (default) reads the criterion at the levels at or above
#'   the range floor in \code{attr(profile, "bounds")}, which is where
#'   \code{resolution_profile()} started the ladder by default up to version
#'   2.0.0.  It applies when the layer supports the floor, some scored level
#'   reaches it, and the profile built its own ladder (a ladder set with
#'   \code{levels} is read as given).  \code{FALSE} reads every level.
#' @param diagnostic \code{criterion = "moran_z"} is refused, with an error
#'   of class \code{spatialkit_diagnostic_criterion} that says what the
#'   column measures, unless this is \code{TRUE}: \code{moran_z} is a
#'   diagnostic of the cell means, not a count rule.  Default \code{FALSE}.
#' @param obs_per_cell For \code{criterion = "size"}: the mean number of
#'   rows per cell wanted, a single number of at least 1 (fractions are
#'   allowed).  \code{NULL} (default) is the profile's own size when it was
#'   built at exactly one (\code{\link{resolution_profile}(obs_per_cell =
#'   )}), and 4.5 otherwise.  Given with a criterion that is not
#'   \code{"size"} it is an error.
#' @return A list of class \code{resolution_selection} with \code{best} (the
#'   level), \code{flat} (the levels in the flat region, ascending; a set,
#'   which can skip a rung), \code{criterion}, \code{value} (the optimum),
#'   \code{at_ceiling} and \code{at_floor} (logical: the optimum is the last
#'   or first of the levels this criterion was read at, which for
#'   \code{moran_z} starts above nine cells and for \code{reliability} at
#'   the range floor; both \code{NA} under \code{"size"}, which searches for
#'   no optimum), \code{edge} (which bound that
#'   is, in words: the support ceiling, the subsample's ceiling, the range
#'   floor, the ladder's own end, the finest or the coarsest size offered
#'   on a ladder \code{obs_per_cell} set, the only level scored,
#'   or the first or last level the criterion is computable at; \code{NA} for
#'   an interior optimum and under \code{"size"}), \code{n_levels},
#'   \code{values} (the criterion
#'   at every level, \code{NA} where it could not be computed),
#'   \code{structured} (for \code{cp}: \code{FALSE} when \code{cp} at one
#'   cell is within \code{tol} of the optimum, when the cell means at the
#'   optimum differ from one mean by no more than chance allows, or when
#'   the variogram found no spatial structure; \code{TRUE} when the optimum
#'   beats one cell on both counts; see "What each criterion returns".
#'   For \code{loo} and \code{size} on a profile built with
#'   \code{estimate}: \code{FALSE} when the leave-one-out error of one mean
#'   for all scored rows is within \code{tol} of the error at \code{best},
#'   or below it, or when the cell means at \code{best} differ from one
#'   mean by no more than chance allows; \code{TRUE} when the cells beat
#'   one mean on both counts.
#'   \code{TRUE} for an \code{elbow}; \code{NA} for \code{reliability}
#'   and \code{moran_z}, and for a profile that carries no one-cell value),
#'   \code{one_cell} (\code{cp} at one cell under \code{cp}; the
#'   leave-one-out error of one mean for all scored rows under \code{loo}
#'   and \code{size}; \code{NA} otherwise),
#'   \code{structure_p} (the p-value of the F test of the cells at
#'   \code{best} against one mean, multiplied by the number of levels
#'   searched and capped at 1; a size searches one level; \code{NA} where
#'   \code{structured} is not read off \code{cp} or the leave-one-out
#'   error),
#'   \code{range_floor} (the floor \code{reliability} was held to, or
#'   \code{NA}), \code{note} (what the flags mean, in words: a character
#'   vector, one note an element; \code{NULL} when there is nothing to
#'   say.  The notes written under \code{"loo"} and \code{"size"}, and for
#'   a profile built with \code{estimate} or with one size, carry names:
#'   \code{size} (the size that set the count), \code{held} (a size the
#'   profile held to the most cells k-means can place), \code{structure}
#'   (why the count is not identified), \code{finest} (the lowest error of
#'   pooled values is at the finest level scored), \code{estimate}
#'   (\code{cp} scores plain means, on a profile of pooled values) and
#'   \code{ladder} (the ladder is one size).  A note of \code{cp} or
#'   \code{moran_z} about the criterion itself has no name.  The warning
#'   of class \code{spatialkit_no_structure} that the functions taking the
#'   selection as a count raise carries every note but those named
#'   \code{size}, \code{held}, \code{finest}, \code{estimate} and
#'   \code{ladder}) and
#'   \code{seeds} (an sf POINT layer, with \code{seed_id}, of the centres of
#'   the partition the profile scored at \code{best}, in the CRS it was
#'   computed in; \code{NULL} for a profile that does not carry them).  Their
#'   Voronoi cells are the cells that were scored;
#'   \code{\link{get_voronoi_seeds}(method = "kmeans", n = <this>)} and
#'   \code{\link{voronoi_seeds_kmeans}(k = <this>)} return them.  For a
#'   profile built with \code{partition = "density"} or \code{init =
#'   "random"} the list also has \code{partition} (\code{kind},
#'   \code{init}, \code{density_power}, \code{density_k}), and the seeds
#'   carry it as \code{attr(, "partition")}; the regions are those the
#'   profile built, whatever kind built them.
#'
#'   For \code{"loo"} and \code{"size"}, and for any criterion on a profile
#'   built with \code{estimate}, the list has these elements after
#'   \code{seeds} (\code{partition}, when there is one, stays last):
#'   \code{estimate} (the profile's, \code{"mean"} or \code{"pooled"};
#'   \code{NA} for a profile without one);
#'   \code{obs_per_cell} (rows per cell at \code{best}: the rows a size
#'   divides over \code{best}; \code{NA} for an object that keeps no row
#'   count);
#'   \code{obs_per_cell_asked} (the size asked for; \code{NA} unless the
#'   criterion is \code{"size"});
#'   \code{scored_per_cell} (the profile's column at \code{best}: the rows
#'   the leave-one-out error read over the cells that hold one.  It is not
#'   \code{obs_per_cell}: it leaves out rows the error did not read and
#'   cells that hold none of them, as under \code{select_on = "split"});
#'   \code{loo} and \code{loo_one_cell} (the leave-one-out error at
#'   \code{best} and that of one mean for all scored rows);
#'   \code{share_alone} and \code{own_weight} (the profile's columns at
#'   \code{best}); and \code{count_rule}, a list of class
#'   \code{spatialkit_count_rule} that records how the count was set:
#'   \code{rule} (the criterion applied), \code{estimate} (the estimate the
#'   rule was computed for: the profile's for \code{"size"} and
#'   \code{"loo"}, \code{"mean"} for \code{"cp"}, \code{NA} otherwise),
#'   \code{cells}, \code{obs_per_cell}, \code{obs_per_cell_asked},
#'   \code{response}, \code{n} (the rows the size divides), \code{n_points},
#'   \code{n_fit} (the rows the centres were fitted to; fewer than
#'   \code{n_points} for a subsample), \code{seed}, \code{features},
#'   \code{scored} (what the profile scored at \code{best} for its
#'   estimate: \code{estimate}, \code{n}, \code{cells}, \code{loo},
#'   \code{loo_se}, \code{loo_one_cell}, \code{share_alone},
#'   \code{own_weight}, \code{structured}, \code{rows} and \code{centres},
#'   the matrix of the centres the rows were scored in, in \code{crs} and
#'   in the profile's order; \code{NULL} without an estimate),
#'   \code{partition}, \code{crs} and
#'   \code{scored_cells} (\code{NA}: the record stands for a count).
#'   \code{seeds} carries the same record as \code{attr(, "count_rule")}
#'   with \code{scored_cells = TRUE}, since they are the centres that were
#'   scored.  The record travels on with the seeds, the cells and the
#'   assigned layer to \code{\link{summarize_by_cell}()}, which warns when
#'   it reports another cell value than the count was chosen for (class
#'   \code{spatialkit_estimate_mismatch}) and when the cells or rows it is
#'   given are not the ones the count was scored on
#'   (\code{spatialkit_cells_differ}).  Under \code{"loo"} and \code{"size"} \code{one_cell} holds
#'   the leave-one-out error of one mean and \code{range_floor} is
#'   \code{NA}.  \code{criterion} is the criterion applied, never
#'   \code{"auto"}.
#' @family aggregation
#' @examples
#' if (requireNamespace("gstat", quietly = TRUE)) {
#'   library(sf)
#'   # An exponential field with range parameter 200 (true effective range
#'   # 600 m) on a 1 km square, with a nugget of 0.6 on a unit sill.
#'   set.seed(4)
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
#'   print(sel$edge)          # the bound the optimum sits on; NA when interior
#'   print(sel$structured)    # FALSE when one cell for the layer does as well
#'
#'   # Reliability prefers coarse cells: it is read at or above the floor the
#'   # autocorrelation range sets, and the result says in words that the
#'   # floor is choosing.
#'   rel <- select_resolution(prof, criterion = "reliability")
#'   print(rel)
#'   print(c(at_floor = rel$at_floor, edge = rel$edge))
#'
#'   # The default profile is one size scored for pooled values, and the
#'   # default criterion, "auto", reads it by that size.
#'   select_resolution(resolution_profile(pts, response_var = "z"))
#' }
#' @export
select_resolution <- function(profile,
                              criterion = c("auto", "cp", "reliability", "elbow",
                                            "moran_z", "loo", "size"),
                              tol = 0.02, range_floor = TRUE, diagnostic = FALSE,
                              obs_per_cell = NULL) {
  if (!inherits(profile, "resolution_profile"))
    stop("select_resolution(): `profile` must come from resolution_profile().",
         call. = FALSE)
  criterion <- match.arg(criterion)
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0)
    stop("select_resolution(): `tol` must be a single non-negative number.", call. = FALSE)
  if (!is.logical(range_floor) || length(range_floor) != 1L || is.na(range_floor))
    stop("select_resolution(): `range_floor` must be TRUE or FALSE.", call. = FALSE)
  if (!is.logical(diagnostic) || length(diagnostic) != 1L || is.na(diagnostic))
    stop("select_resolution(): `diagnostic` must be TRUE or FALSE.", call. = FALSE)
  if (!is.null(obs_per_cell) &&
      !(is.numeric(obs_per_cell) && !is.object(obs_per_cell) && length(obs_per_cell) == 1L &&
        is.finite(obs_per_cell) && obs_per_cell >= 1))
    stop(sprintf(paste0("select_resolution(): `obs_per_cell` must be NULL or a single ",
                        "number of at least 1; got %s."), .got_text(obs_per_cell)),
         call. = FALSE)
  # The estimate the profile scores, when it names one.  "auto" is the rule
  # that goes with it: a size for pooled values, the leave-one-out minimum
  # for plain means, and Cp for a profile that names none (one built without
  # `estimate`, and a column subset, which has lost its attributes).
  est <- attr(profile, "estimate", exact = TRUE)
  est_value <- .profile_estimate(profile)
  asked <- criterion
  if (criterion == "auto")
    criterion <- if (identical(est_value, "pooled")) "size"
                 else if (identical(est_value, "mean")) "loo"
                 else "cp"
  if (!is.null(obs_per_cell) && criterion != "size")
    stop(sprintf(paste0("select_resolution(): `obs_per_cell` is read by criterion = \"size\" ",
                        "only, and criterion = \"%s\" was asked for%s."),
                 asked,
                 if (asked == "auto")
                   sprintf(", which reads this profile by \"%s\"", criterion)
                 else ""),
         call. = FALSE)
  # `[` keeps the class, so a subset arrives here looking like a profile.  An
  # empty one or a column subset missing the criterion used to fall through to
  # the "NA at every level" error below, which sends the user off to add a
  # response or fix a variogram when the real problem is the object.
  if (nrow(profile) == 0L)
    stop("select_resolution(): `profile` has no levels; it is an empty subset ",
         "of a resolution profile.", call. = FALSE)
  if (criterion == "size")
    return(.select_size(profile, obs_per_cell, tol, est_value))
  if (criterion == "loo" && !"loo" %in% names(profile))
    stop(paste0("select_resolution(): `profile` has no `loo` column; build it with ",
                "resolution_profile(estimate = \"mean\") or estimate = \"pooled\"."),
         call. = FALSE)
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
  one <- attr(profile, "one_cell", exact = TRUE)
  # The sizes of a ladder `obs_per_cell` set, when the profile is on one.
  ladder_sizes <- if (is.list(attr(profile, "bounds", exact = TRUE)))
    attr(profile, "bounds", exact = TRUE)$obs_per_cell else NULL
  # A constant response is its own case: Cp is NA because no level fits a
  # constant better than another, and the message used to ask for "a level
  # whose cells hold two scored rows", which every level of that profile had.
  if (!any(ok) && criterion %in% c("cp", "moran_z") && is.list(one) && isTRUE(one$constant))
    stop(sprintf(paste0("select_resolution(): `%s` is NA at every level because the %s is ",
                        "constant: its total sum of squares is rounding error at most, ",
                        "so no cell count fits it better than another."),
                 criterion,
                 if (identical(attr(profile, "variable"), "residuals"))
                   "OLS residual of the response on the predictors" else "response"),
         call. = FALSE)
  # The same for the leave-one-out error, whose reference is the estimate's:
  # the rows it read do not vary.
  if (!any(ok) && criterion == "loo" && is.list(est) && isTRUE(est$constant))
    stop(sprintf(paste0("select_resolution(): `%s` is NA at every level because the %s is ",
                        "constant: its total sum of squares is rounding error at most, ",
                        "so no cell count fits it better than another."),
                 criterion, "response"),
         call. = FALSE)
  # A ladder of sizes that was given no `sac` fitted no variogram, so Cp has
  # no noise variance and the reliability no range.  Said here, with the
  # ways to get the column, and not as a profile that lacks a response.
  if (!any(ok) && criterion %in% c("cp", "reliability") &&
      identical(attr(profile, "cp_noise", exact = TRUE)$source, "none"))
    stop(sprintf(paste0("select_resolution(): `%s` is NA at every level: this profile is on ",
                        "a size ladder (obs_per_cell = %s), which fits no variogram unless ",
                        "`sac` is given, and `%s` takes its %s from one. Build the profile ",
                        "with `sac = NULL` to have the variogram estimated, with a `sac` of ",
                        "your own, or with obs_per_cell = NULL for the automatic ladder; ",
                        "criterion = \"size\" reads this profile."),
                 criterion, .got_text(ladder_sizes), criterion,
                 if (criterion == "cp") "noise variance" else "range"),
         call. = FALSE)
  if (!any(ok))
    stop(sprintf(paste0("select_resolution(): `%s` is NA at every level%s."),
                 criterion,
                 switch(criterion,
                        cp = paste0(" (it needs a response, and a level whose cells hold ",
                                    "two scored rows each on average for its noise ",
                                    "variance when the variogram gives no nugget)"),
                        reliability = " (it needs a usable variogram with an identified range)",
                        moran_z = " (it needs a response and more than nine cells)",
                        elbow = paste0(" (the WSS curve has no elbow: on log-log axes it ",
                                       "falls in a straight line, as it does for points ",
                                       "with no cluster structure, so the geometry names ",
                                       "no cell count)"),
                        loo = paste0(" (it needs a level at which some cell holds two ",
                                     "rows with a response)"),
                        "")), call. = FALSE)
  maximise <- criterion %in% c("reliability", "elbow")
  lv <- profile$levels
  bounds <- attr(profile, "bounds")
  with_pred <- identical(attr(profile, "variable"), "residuals")
  # moran_z is Moran's z of the cell means.  Its smallest |z| is where that
  # test has least power, about 10 to 20 cells, whatever the layer, and it
  # used to be handed over as a cell count like any other criterion's.
  what_z <- if (with_pred)
    "the residuals of the cell means regressed on the cell-mean predictors"
  else "the cell means"
  if (criterion == "moran_z" && !diagnostic)
    stop(errorCondition(
      sprintf(paste0("select_resolution(): `moran_z` is a diagnostic, not a count rule, ",
                     "and is read only when asked for with `diagnostic = TRUE`. It is ",
                     "Moran's z of %s over each cell's 8 nearest cells: a test of whether ",
                     "neighbouring cells are independent. It is near 0 at about 10 to 20 ",
                     "cells whatever the data, where that neighbour graph is close to ",
                     "complete, and grows as the cells get finer than the range, so its ",
                     "smallest |z| (here at %d cells; it is computable from %d and the ",
                     "ladder runs to %d) does not say how well the cells represent the ",
                     "response. `cp` does."),
              what_z, as.integer(lv[ok][which.min(v[ok])]), as.integer(min(lv[ok])),
              as.integer(max(lv))),
      class = "spatialkit_diagnostic_criterion", call = NULL))
  # Reliability tends to peak at a few cells, below the range floor, so on a
  # ladder from 2 its maximum is a count whose cells are wider than the
  # range.  It is read at or above the floor the profile reports, which is
  # where the ladder itself used to start.  A ladder the caller set with
  # `levels` is read as given, as it always was.
  held_floor <- NA_integer_
  if (criterion == "reliability" && range_floor && is.list(bounds) &&
      isTRUE(bounds$supported) && isFALSE(bounds$levels_given) &&
      is.finite(bounds$range %||% NA_real_) && is.finite(bounds$floor %||% NA_integer_) &&
      any(ok & lv >= bounds$floor)) {
    ok <- ok & lv >= bounds$floor
    held_floor <- as.integer(bounds$floor)
  }
  opt <- if (maximise) max(v[ok]) else min(v[ok])
  best <- lv[ok][if (maximise) which.max(v[ok]) else which.min(v[ok])]
  # Of levels with the same leave-one-out error, the one with the fewest
  # cells, whatever order the rows of the profile are in.
  if (criterion == "loo") best <- min(lv[ok][v[ok] == opt])
  flat <- if (criterion %in% c("cp", "reliability", "loo")) {
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
  edge <- .ladder_edge(best, scored, lv, bounds, held_floor = held_floor)
  # Whether the criterion found anything to choose on.  For Cp: does any
  # level beat one cell for the whole layer by more than `tol`?  When one
  # cell is inside the flat region, one mean predicts as well as the cells
  # and the argmin is a draw from the noise.  Compared only while the `cp`
  # column is the one the one-cell value was computed with (a profile whose
  # column was edited has no such reference).
  structured <- NA
  one_cp <- NA_real_
  structure_p <- NA_real_
  note <- NULL
  if (criterion == "cp") {
    ref <- if (is.list(one)) one$cp_levels else NULL
    same <- is.list(one) && is.finite(one$cp %||% NA_real_) && !is.null(ref) &&
      all(as.character(lv) %in% names(ref)) &&
      isTRUE(all.equal(unname(ref[as.character(lv)]), as.numeric(profile[["cp"]])))
    if (same) {
      one_cp <- one$cp
      # Two things have to hold for the count to be called identified.  Cp at
      # the optimum is below Cp at one cell by more than `tol`; and the cell
      # means there differ from one mean by more than chance, by the F test
      # the profile kept at every level, at 5 percent after a Bonferroni
      # adjustment for the levels the optimum was searched over.  The first
      # alone has no sampling error in it: on white noise it called the count
      # identified on 14 percent of layers of 60 points, 10 of 120 and 5 of
      # 300 (both together: 2, 4 and 3), and on every layer whose Cp had a
      # noise variance of 0, where Cp is RSS / n and a finer level beats one
      # cell by more than `tol` by chance alone.
      beats <- !(one_cp <= opt * (1 + tol))
      n_searched <- sum(ok)
      p_raw <- if (is.null(one$f_p)) NA_real_ else unname(one$f_p[as.character(best)])[1L]
      if (is.finite(p_raw)) structure_p <- min(1, p_raw * n_searched)
      by_chance <- is.finite(structure_p) && structure_p >= 0.05
      structured <- if (!beats || by_chance) FALSE
                    else if (is.finite(structure_p)) TRUE
                    else NA
      p_text <- if (is.finite(structure_p))
        sprintf("p = %s after adjusting for the %d level%s searched",
                format(signif(structure_p, 2)), n_searched, if (n_searched == 1L) "" else "s")
      else ""
      # Cp's noise variance, when it is no penalty at all.
      cn <- attr(profile, "cp_noise", exact = TRUE)
      vgm <- attr(profile, "variogram", exact = TRUE)
      no_penalty <- is.list(cn) && is.finite(cn$value %||% NA_real_) &&
        (cn$value <= 0 ||
           (is.list(vgm) && is.finite(vgm$psill %||% NA_real_) &&
              cn$value <= 1e-4 * (cn$value + vgm$psill)))
      if (!beats)
        note <- sprintf(paste0("Cp at one cell (%s) is within %s%% of the optimum (%s): one ",
                               "mean for the whole layer predicts as well as any level ",
                               "scored, so the response shows no spatial structure the ",
                               "cells can use and the count is not identified."),
                        format(signif(one_cp, 3)), format(100 * tol), format(signif(opt, 3)))
      else if (by_chance)
        note <- sprintf(paste0("the cell means at %d cells differ from one mean for the ",
                               "whole layer by no more than chance allows (F test of the ",
                               "cells against one mean, %s), so the response shows no ",
                               "spatial structure the cells can use and the count is not ",
                               "identified. %s"),
                        as.integer(best), p_text,
                        if (no_penalty)
                          paste0("Cp's noise variance is 0: Cp has no penalty and falls ",
                                 "with every finer level whatever the response.")
                        else sprintf(paste0("Cp at one cell is %s against %s at the optimum, ",
                                            "a gap that chance produces on a layer of this size."),
                                     format(signif(one_cp, 3)), format(signif(opt, 3))))
      else if (isTRUE(structured) && no_penalty)
        note <- sprintf(paste0("Cp's noise variance is 0, so Cp is RSS / n and falls with ",
                               "every finer level by chance alone: its comparison with one ",
                               "cell says nothing, and `structured` rests on the F test of ",
                               "the cells against one mean alone (%s)."), p_text)
      else if (is.na(structured))
        note <- paste0("Cp at one cell is above the optimum by more than `tol`, but the ",
                       "profile carries no F test of the cells against one mean at this ",
                       "level, so whether that is more than chance is not known.")
    }
    sac_prof <- attr(profile, "sac", exact = TRUE)
    why_sac  <- attr(sac_prof, "rejected_reason")
    if (.sac_no_structure(sac_prof)) {
      structured <- FALSE
      note <- c(note, sprintf(paste0("the profile's variogram found no spatial structure in ",
                                     "the response (%s), so the count is not identified."),
                              as.character(why_sac)[1L]))
    }
  } else if (criterion == "elbow") {
    structured <- TRUE
  } else if (criterion == "moran_z") {
    note <- sprintf(paste0("moran_z is a diagnostic, not a count rule: this is the level ",
                           "at which %s are least autocorrelated across neighbouring ",
                           "cells, not the level that represents the response best."),
                    what_z)
  } else if (criterion == "loo") {
    # The same two questions as for Cp, asked of the leave-one-out error of
    # the profile's estimate: is it below that of one mean for every scored
    # row by more than `tol`, and do the cell means at the optimum differ
    # from one mean by more than chance, after the adjustment for the levels
    # searched.  Read only while the `loo` column is the one the references
    # were recorded with.
    n_searched <- sum(ok)
    st <- .estimate_structure(profile, best, tol, n_searched = n_searched)
    structured  <- st$structured
    one_cp      <- st$one_cell
    structure_p <- st$structure_p
    p_text <- if (is.finite(structure_p))
      sprintf("p = %s after adjusting for the %d level%s searched",
              format(signif(structure_p, 2)), n_searched, if (n_searched == 1L) "" else "s")
    else ""
    if (isFALSE(st$beats))
      note <- c(structure = sprintf(
        paste0("the leave-one-out error at one cell (%s) is within %s%% of the optimum ",
               "(%s): one mean for the whole layer predicts as well as any level scored, ",
               "so the response shows no spatial structure the cells can use and the ",
               "count is not identified."),
        format(signif(one_cp, 3)), format(100 * tol), format(signif(opt, 3))))
    else if (isFALSE(structured))
      note <- c(structure = sprintf(
        paste0("the cell means at %d cells differ from one mean for the whole layer by ",
               "no more than chance allows (F test of the cells against one mean, %s), so ",
               "the response shows no spatial structure the cells can use and the count ",
               "is not identified. The leave-one-out error at one cell is %s against %s ",
               "at the optimum, a gap that chance produces on a layer of this size."),
        as.integer(best), p_text, format(signif(one_cp, 3)), format(signif(opt, 3))))
    else if (isTRUE(st$beats) && is.na(structured))
      note <- c(structure = paste0(
        "the leave-one-out error at one cell is above the optimum by more than `tol`, ",
        "but the profile carries no F test of the cells against one mean at this level, ",
        "so whether that is more than chance is not known."))
    # For pooled values the error can keep falling to the end of the ladder.
    # When it did, the ladder's end set the count.
    if (identical(est_value, "pooled") && best == max(lv[ok]))
      note <- c(note, finest = paste0(
        "the lowest leave-one-out error of the pooled value is at the finest level ",
        "scored, so the end of the ladder chose the count. criterion = \"size\" sets ",
        "it by a size."))
  }
  # Cp scores plain means whatever the profile's estimate is.
  if (criterion == "cp" && identical(est_value, "pooled"))
    note <- c(note, estimate = paste0(
      "cp scores plain means, and this profile's estimate is the pooled value; ",
      "criterion = \"size\" is the rule for it."))
  # A ladder of one size is one level, so the criterion had nothing to choose
  # among.
  if (criterion %in% c("cp", "reliability") && length(ladder_sizes) == 1L)
    note <- c(note, ladder = sprintf(
      paste0("the ladder is one size (obs_per_cell = %s), so %s had one level to choose ",
             "from. Build the profile with obs_per_cell = NULL, or with several sizes, ",
             "for a ladder %s can choose on."),
      .got_text(ladder_sizes), criterion, criterion))
  sel <- structure(list(best = as.integer(best), flat = sort(as.integer(flat)),
                 criterion = criterion, value = opt,
                 at_ceiling = at_ceiling, at_floor = at_floor, edge = edge,
                 n_levels = length(lv),
                 values = stats::setNames(v, lv),
                 structured = structured, one_cell = one_cp,
                 structure_p = structure_p,
                 range_floor = held_floor, note = note,
                 # The centres of the partition this level was scored on, for
                 # get_voronoi_seeds() and voronoi_seeds_kmeans() to seed the
                 # same cells from.  NULL for a profile made before they were
                 # kept.
                 seeds = .profile_seeds(profile, best)),
            class = "resolution_selection")
  .selection_close(sel, profile, est_value, tol)
}


#' The estimate a profile names, by value
#'
#' @param profile A \code{resolution_profile}, or a subset of one.
#' @return \code{"mean"} or \code{"pooled"}; \code{NA} for a profile built
#'   without an estimate, for a column subset, which has lost the attribute,
#'   and for an attribute that names neither.
#' @keywords internal
#' @noRd
.profile_estimate <- function(profile) {
  est <- attr(profile, "estimate", exact = TRUE)
  if (is.list(est) && is.character(est$value) && length(est$value) == 1L &&
      est$value %in% c("mean", "pooled"))
    as.character(unclass(est$value)) else NA_character_
}


#' Whether `loo` is read where no criterion is named
#'
#' \code{summary()} and \code{plot()} of a profile take \code{loo} into
#' their default set where it is the leave-one-out error
#' \code{resolution_profile()} wrote: the profile names an estimate, or the
#' object is a column subset that kept \code{loo} and none of \code{cp},
#' \code{reliability} and \code{elbow}.  A column called \code{loo} that a
#' caller added to a profile built without an estimate is read only when
#' \code{criteria} names it, so the table and the plot of such a profile are
#' the ones they were before the column had a meaning here.
#'
#' @param profile A \code{resolution_profile}, or a subset of one.
#' @return \code{TRUE} or \code{FALSE}.
#' @keywords internal
#' @noRd
.loo_by_default <- function(profile)
  "loo" %in% names(profile) &&
    (is.list(attr(profile, "estimate", exact = TRUE)) ||
       !any(c("cp", "reliability", "elbow") %in% names(profile)))


# The names select_resolution() gives the notes it writes about something
# other than whether the count is identified: the size that set the count, a
# size the profile held, an optimum at the finest level of a pooled profile,
# the estimate `cp` scores, and a ladder of one size.  The note that says why
# a count is not identified is named `structure` under "size" and "loo" and
# has no name under the criteria that wrote one before.  The warning raised
# when the count becomes cells pastes every note but these five.
.NOTES_NOT_STRUCTURE <- c("size", "held", "finest", "estimate", "ladder")


#' What a selection carries beyond the criterion's own answer
#'
#' On the path that names an estimate (the profile carries one, or the
#' criterion is \code{"loo"} or \code{"size"}) the selection gains what the
#' profile scored at the chosen level and the record of how the count was
#' set, and its seeds carry that record.  The kind that built the profile's
#' cells, when the profile records one, is the last element either way.  A
#' selection from a profile without an estimate by \code{"cp"},
#' \code{"reliability"}, \code{"elbow"} or \code{"moran_z"} gains nothing
#' but that kind, as before.
#'
#' @param sel The selection, with the elements every criterion returns.
#' @param profile The \code{resolution_profile} it was read off.
#' @param est_value The profile's estimate, \code{"mean"} or
#'   \code{"pooled"}; \code{NA} for a profile that names none.
#' @param tol The tolerance of the call.
#' @param asked The size asked for; \code{NA} unless the criterion is
#'   \code{"size"}.
#' @return The selection.
#' @keywords internal
#' @noRd
.selection_close <- function(sel, profile, est_value, tol, asked = NA_real_) {
  if (!is.na(est_value) || sel$criterion %in% c("loo", "size")) {
    best <- sel$best
    n    <- .count_rule_rows(attr(profile, "bounds", exact = TRUE))
    est  <- attr(profile, "estimate", exact = TRUE)
    sel$estimate           <- est_value
    # Rows per cell at the count chosen: the rows a size divides over the
    # centres asked for.  `scored_per_cell` is another number: the rows the
    # leave-one-out error read, over the cells that hold one.
    sel$obs_per_cell       <- if (is.na(n)) NA_real_ else n / best
    sel$obs_per_cell_asked <- as.numeric(asked)
    sel$scored_per_cell    <- .profile_at(profile, best, "scored_per_cell")
    sel$loo                <- .profile_at(profile, best, "loo")
    sel$loo_one_cell       <- if (!is.na(est_value) && "loo" %in% names(profile))
      as.numeric(est$loo_one_cell %||% NA_real_)[1L] else NA_real_
    sel$share_alone        <- .profile_at(profile, best, "share_alone")
    sel$own_weight         <- .profile_at(profile, best, "own_weight")
    # Whether the cells at this count beat one mean for the profile's
    # estimate.  "size" and "loo" asked that themselves; under another
    # criterion it is asked here of the one level chosen.
    scored_structured <- if (sel$criterion %in% c("loo", "size")) sel$structured
                         else .estimate_structure(profile, best, tol, n_searched = 1L)$structured
    rec <- .count_rule_new(profile, rule = sel$criterion, cells = best,
                           obs_per_cell_asked = asked, structured = scored_structured,
                           scored_cells = NA)
    sel$count_rule <- rec
    # The seeds are the centres the profile scored at this count.
    if (inherits(sel$seeds, "sf")) {
      rec$scored_cells <- TRUE
      attr(sel$seeds, "count_rule") <- rec
    }
  }
  # The kind that built the profile's cells, when it is not the present one
  # (see resolution_profile()): the selection names it, and its seeds carry
  # it on.
  part <- attr(profile, "partition", exact = TRUE)
  if (is.list(part) && !is.null(part$kind))
    sel$partition <- part[.partition_fields()]
  sel
}


#' The count a size gives, read off a profile
#'
#' \code{select_resolution(criterion = "size")}.  The count is
#' \code{.count_at(n, obs_per_cell)}, with \code{n} the rows the size
#' divides, and the level returned is the first of three that applies: the
#' level the profile built for that size; the count itself, when the profile
#' has a partition at it; the built level nearest to the count on the log
#' scale, with a warning.  No optimum is searched for, so no bound can
#' choose the count and the edge flags are \code{NA}.
#'
#' A size does not ask whether the response has structure the cells can
#' use, so that is asked separately, of the profile's estimate at the level
#' returned: \code{.estimate_structure()} with one level searched.
#'
#' @param profile A \code{resolution_profile}.
#' @param obs_per_cell The size asked for, or \code{NULL} for the profile's
#'   own size when it was built at exactly one, and
#'   \code{.OBS_PER_CELL_DEFAULT} otherwise.
#' @param tol The margin of the structure test.
#' @param est_value The profile's estimate, or \code{NA}.
#' @return A \code{resolution_selection}.
#' @keywords internal
#' @noRd
.select_size <- function(profile, obs_per_cell, tol, est_value) {
  bounds <- attr(profile, "bounds", exact = TRUE)
  n <- .count_rule_rows(bounds)
  if (is.na(n))
    stop(paste0("select_resolution(): criterion = \"size\" needs the number of rows the ",
                "size divides, which a profile keeps in attr(, \"bounds\"); this object ",
                "has none (it is a column subset of a resolution profile). Pass the ",
                "profile resolution_profile() returned."), call. = FALSE)
  if (!"levels" %in% names(profile))
    stop(paste0("select_resolution(): `profile` has no `levels` column; it is a column ",
                "subset of a resolution profile. Pass the profile resolution_profile() ",
                "returned."), call. = FALSE)
  sizes <- bounds$obs_per_cell
  size  <- if (!is.null(obs_per_cell)) as.numeric(obs_per_cell)
           else if (length(sizes) == 1L) as.numeric(sizes)
           else .OBS_PER_CELL_DEFAULT
  k_star <- .count_at(n, size)
  lv <- profile$levels
  # The levels at which a partition was built.
  built <- if ("wss" %in% names(profile)) is.finite(suppressWarnings(as.numeric(profile$wss)))
           else rep(TRUE, length(lv))
  if (!any(built))
    stop(paste0("select_resolution(): no level of `profile` has a partition (`wss` is NA ",
                "at every level), so no level can stand for a size."), call. = FALSE)
  # (i) The level the profile built for this size.  It is the count the size
  # gives unless the profile held the size to the most cells k-means can
  # place; the profile said so when it was built, so nothing is raised here.
  at <- match(size, as.numeric(sizes))
  own <- if (!is.na(at) && length(bounds$size_levels) == length(sizes))
    as.integer(bounds$size_levels[[at]]) else NA_integer_
  held <- FALSE
  if (!is.na(own) && own %in% lv[built]) {
    best <- own
    held <- best != k_star
  } else if (k_star %in% lv[built]) {
    # (ii) The count itself is a level.
    best <- k_star
  } else {
    # (iii) The built level nearest to the count on the log scale; of two as
    # near, the one with fewer cells.  Ratios, not differences of logarithms:
    # 2 and 8 are as near to 4 as each other, and their logarithms differ
    # from log(4) by amounts that need not round alike.
    cand <- sort(unique(lv[built]))
    best <- cand[which.min(pmax(cand, k_star) / pmin(cand, k_star))]
    # The most cells k-means can place on the layer: its distinct locations
    # where some repeat, one short of the points where none does.  A profile
    # fitted to every point counted the locations; one fitted to a subsample
    # counted the subsample's, and one short of the points is then still a
    # bound.  A count above it is one no profile scores, and the advice to
    # build one at that size is left out.
    one <- function(v) is.numeric(v) && length(v) == 1L && is.finite(v)
    repeats <- one(bounds$n) && one(bounds$n_distinct) && one(bounds$n_sample) &&
      bounds$n_sample >= bounds$n && bounds$n_distinct < bounds$n
    cap <- if (!one(bounds$n)) NA_integer_
           else if (repeats) as.integer(bounds$n_distinct)
           else as.integer(bounds$n) - 1L
    .warn_classed("spatialkit_size_not_reached",
                  paste0("select_resolution(): obs_per_cell = %s asks for %d cells, which ",
                         "is not a level of this profile; the nearest level, %d cells (%s), ",
                         "is used. %s"),
                  format(size), k_star, as.integer(best), .per_cell_text(n / best),
                  if (!is.na(cap) && k_star > cap)
                    sprintf(paste0("That count is more than k-means can place on %s (%d ",
                                   "cells), so no profile scores it."),
                            if (repeats)
                              sprintf("the layer's %d distinct locations",
                                      as.integer(bounds$n_distinct))
                            else sprintf("%d points", as.integer(bounds$n)), cap)
                  else sprintf("Build the profile with obs_per_cell = %s to score that count.",
                               format(size)),
                  data = list(obs_per_cell = size, wanted = k_star, used = as.integer(best),
                              reason = "not a level"))
  }
  best <- as.integer(best)
  note <- c(size = sprintf(paste0("the count is set by obs_per_cell = %s, the caller's ",
                                  "choice; no criterion chose it."), format(size)))
  if (held)
    note <- c(note, held = sprintf(
      paste0("the profile holds this size to %d cells, the most k-means can place on ",
             "this layer; obs_per_cell = %s asks for %d."), best, format(size), k_star))
  # Whether cells of this size beat one mean for the whole layer.  The count
  # stands either way: the caller set it.
  st <- .estimate_structure(profile, best, tol, n_searched = 1L)
  if (isFALSE(st$structured))
    note <- c(note, structure = paste0(
      if (isFALSE(st$beats))
        sprintf(paste0("at %d cells the leave-one-out error of the %s (%s) is not below ",
                       "that of one mean for the whole layer (%s) by more than %s%%"),
                best, if (identical(est_value, "pooled")) "pooled value" else "plain mean",
                format(signif(st$loo, 3)), format(signif(st$one_cell, 3)), format(100 * tol))
      else
        sprintf(paste0("the cell means at %d cells differ from one mean for the whole ",
                       "layer by no more than chance allows (F test, p = %s)"),
                best, format(signif(st$structure_p, 2))),
      ": the response shows no spatial structure that cells of this size can use."))
  sel <- structure(list(best = best, flat = best,
                 criterion = "size", value = n / best,
                 at_ceiling = NA, at_floor = NA, edge = NA_character_,
                 n_levels = length(lv),
                 values = stats::setNames(n / lv, lv),
                 structured = st$structured, one_cell = st$one_cell,
                 structure_p = st$structure_p,
                 range_floor = NA_integer_, note = note,
                 seeds = .profile_seeds(profile, best)),
            class = "resolution_selection")
  .selection_close(sel, profile, est_value, tol, asked = size)
}


#' The scored partition's centres at one level, as an sf POINT layer
#'
#' @param profile A \code{resolution_profile}.
#' @param L The level (cell count).
#' @return An sf POINT layer of \code{L} centres in the CRS the profile was
#'   computed in, or \code{NULL} when the profile does not carry them.
#' @keywords internal
#' @noRd
.profile_seeds <- function(profile, L) {
  cen <- attr(profile, "centres", exact = TRUE)
  xy  <- if (is.list(cen)) cen[[as.character(L)]] else NULL
  if (!is.matrix(xy) || nrow(xy) < 1L || ncol(xy) < 2L) return(NULL)
  crs <- attr(profile, "centres_crs", exact = TRUE)
  s <- sf::st_as_sf(data.frame(seed_id = seq_len(nrow(xy)), x = xy[, 1L], y = xy[, 2L]),
                    coords = c("x", "y"),
                    crs = if (inherits(crs, "crs")) crs else sf::NA_crs_)
  attr(s, "nstart") <- attr(profile, "nstart", exact = TRUE)
  # The ground the centres stand for: the hull of the points profiled.  A
  # select_resolution() result holds these centres and not the profile, so
  # without it a selection given to get_voronoi_seeds() with neither points
  # nor a boundary had no cover to hand on.
  attr(s, "cover")  <- attr(profile, "hull", exact = TRUE)
  # The kind that built the cells these centres seed, when the profile
  # records one.
  part <- attr(profile, "partition", exact = TRUE)
  if (is.list(part) && !is.null(part$kind))
    attr(s, "partition") <- part[.partition_fields()]
  s
}

# Which bound an optimum at the end of the scored levels is sitting on.  NA
# when it is interior.  The wording is shared by the print methods, the
# profile plot's caption and the tessellation builders' provenance note.
.ladder_edge <- function(best, scored, ladder, bounds, held_floor = NA_integer_) {
  if (best == max(scored) && best == min(scored) &&
      (!is.finite(held_floor) || length(unique(ladder)) < 2L))
    return("the only level scored")
  # The first level at or above the range floor, for a criterion read from
  # the floor up (reliability): the floor is the bound there, whether or not
  # the ladder has a rung exactly on it.
  if (is.finite(held_floor) && best == min(scored))
    return(if (best == held_floor) "the range floor (area / range^2)"
           else "the first level at or above the range floor (area / range^2)")
  # A ladder `obs_per_cell` set ends where the sizes asked for end: no
  # ceiling of support and no range floor bounds it.  One case apart: a
  # finest size the profile held to the most cells k-means can place.  That
  # level is the layer's bound and not the caller's, and it keeps the label
  # of the support ceiling.
  size_ladder <- is.list(bounds) && !is.null(bounds$obs_per_cell)
  size_held <- size_ladder && identical(bounds$ceiling_from, "distinct locations")
  if (best == max(scored)) {
    if (max(scored) < max(ladder))
      return("the last level the criterion is computable at")
    if (size_ladder && !size_held) return("the finest size offered (obs_per_cell)")
    if (is.list(bounds) && identical(as.integer(max(ladder)), as.integer(bounds$ceiling)))
      return(if (identical(bounds$ceiling_from, "distinct locations") &&
                 isTRUE(bounds$ceiling < (bounds$n_distinct %||% NA_integer_)))
               "the support ceiling (one short of the points, the most k-means can fit)"
             else if (identical(bounds$ceiling_from, "distinct locations"))
               "the support ceiling (as many cells as the distinct locations allow)"
             else if (identical(bounds$ceiling_from, "sample_n") &&
                      isTRUE((bounds$n_scored %||% NA_integer_) < (bounds$n_sample %||% NA_integer_)))
               "the subsample's ceiling (two scored subsample rows per cell; raise sample_n)"
             else if (identical(bounds$ceiling_from, "sample_n"))
               "the subsample's ceiling (two subsample points per cell; raise sample_n)"
             else if (identical(bounds$ceiling_from, "min_cell_n (rows with a response)"))
               (if (isTRUE(bounds$predictors))
                  "the support ceiling (rows with a response and every predictor / min_cell_n)"
                else "the support ceiling (rows with a response / min_cell_n)")
             else "the support ceiling (n / min_cell_n)")
    return("the last level of the ladder")
  }
  if (best == min(scored)) {
    if (min(scored) > min(ladder))
      return("the first level the criterion is computable at")
    if (size_ladder) return("the coarsest size offered (obs_per_cell)")
    if (is.list(bounds) && isTRUE(bounds$supported) && is.finite(bounds$range %||% NA) &&
        !isFALSE(bounds$range_floor) &&
        identical(as.integer(min(ladder)), as.integer(bounds$floor)))
      return("the range floor (area / range^2)")
    return("the first level of the ladder")
  }
  NA_character_
}


# What to do about an optimum on a bound, by the bound.  Shared by the print
# method and by the warning the functions that take the count raise.
.edge_hint <- function(edge) {
  if (is.null(edge) || is.na(edge)) return("")
  if (grepl("min_cell_n", edge, fixed = TRUE))
    " Lower min_cell_n to see whether the criterion keeps going."
  else if (grepl("raise sample_n", edge, fixed = TRUE))
    " Raise sample_n to see whether the criterion keeps going."
  else if (grepl("range floor", edge, fixed = TRUE))
    " Fewer cells would be wider than the range and average over more than one patch of the field."
  else if (grepl("only level", edge, fixed = TRUE))
    " One level leaves nothing to compare it with."
  else if (grepl("(obs_per_cell)", edge, fixed = TRUE))
    " `obs_per_cell` sets the sizes offered; how small a cell may be is the caller's choice."
  else ""
}


#' @export
print.resolution_selection <- function(x, ...) {
  cat(sprintf("Resolution by %s: %d cells\n", x$criterion, x$best))
  by_size <- identical(x$criterion, "size")
  rec <- if (inherits(x$count_rule, "spatialkit_count_rule")) x$count_rule else NULL
  # A size has no flat region: its line gives the size asked for and the
  # rows per cell the count leaves.
  if (by_size) {
    asked <- x$obs_per_cell_asked %||% NA_real_
    n_div <- if (is.null(rec)) NA_integer_ else rec$n
    cat(sprintf("  size        : %s asked (obs_per_cell); %s rows / %d cells = %.2f\n",
                paste(format(asked), if (isTRUE(asked == 1)) "row per cell" else "rows per cell"),
                format(n_div), x$best, x$value))
  } else {
    cat(sprintf("  flat region : %s (%d of %d levels)\n",
                .band_label(x$flat, .scored_levels(x$values)),
                length(x$flat), x$n_levels))
  }
  # The kind that built the cells, when the profile recorded one; a selection
  # from a k-means profile prints as it always did.
  if (is.list(x$partition) && !is.null(x$partition$kind))
    cat(sprintf("  cells       : %s\n", .partition_label(x$partition)))
  # strwrap() collapses runs of spaces, so the aligned label is put back
  # after wrapping rather than wrapped with the sentence.
  say <- function(text, label = "note", width = 62) {
    body <- strwrap(text, width = width)
    cat(paste0(c(sprintf("  %-12s: ", label), rep("                ", length(body) - 1L)), body),
        sep = "\n")
  }
  # What the profile scored at this count for its estimate, under the two
  # criteria that read it.
  est <- x$estimate %||% NA_character_
  if (x$criterion %in% c("size", "loo") && length(est) == 1L && !is.na(est)) {
    cat(sprintf("  estimate    : %s\n",
                if (identical(est, "pooled"))
                  sprintf("pooled toward the %d nearest cells", .POOL_NEIGHBOURS)
                else "plain means"))
    alone <- x$share_alone %||% NA_real_
    # The rows the leave-one-out error read and the cells that hold one.  The
    # line says so where they are not every row and every cell: under the
    # split the error read one half of the rows, and a centre with no scored
    # row near it (rows without a response) is a cell the share says nothing
    # of.
    scored <- if (!is.null(rec) && is.list(rec$scored)) rec$scored else NULL
    half <- identical(scored$rows, "selection half")
    occupied <- scored$cells %||% NA_integer_
    fewer <- length(occupied) == 1L && is.finite(occupied) && occupied < x$best
    per <- x$scored_per_cell %||% NA_real_
    if (is.finite(alone)) {
      w <- x$own_weight %||% NA_real_
      say(paste0(if (half)
                   sprintf("in the selection half (%s rows%s, %.2f per cell): ",
                           format(scored$n),
                           if (fewer) sprintf(" in %d of the %d cells", occupied, x$best) else "",
                           per)
                 else if (fewer)
                   sprintf("in the %d of the %d cells that hold a scored row (%.2f per cell): ",
                           occupied, x$best, per)
                 else "",
                 sprintf("%.1f%% of the rows are alone in their cell", 100 * alone),
                 if (is.finite(w)) sprintf("; median own weight %.2f", w) else ""),
          label = "at this size", width = 76)
    }
    # Left out when the profile has neither number: a constant response.
    loo <- x$loo %||% NA_real_
    one <- x$loo_one_cell %||% NA_real_
    if (isTRUE(is.finite(loo)) || isTRUE(is.finite(one)))
      cat(sprintf("  loo         : %s (one cell for %s: %s)\n",
                  format(signif(loo, 3)), if (half) "that half" else "the layer",
                  format(signif(one, 3))))
  }
  edge <- x$edge %||% NA_character_
  if (!is.na(edge))
    say(sprintf("the optimum is %s; the bound is choosing, not the criterion.%s",
                edge, .edge_hint(edge)))
  # What `structured = FALSE` and a diagnostic criterion mean, in words.
  for (nt in x$note) say(nt)
  invisible(x)
}


#' Every criterion's pick, side by side
#'
#' \code{\link{select_resolution}()} reads one criterion at a time.  This puts
#' the count criteria in one table: the level each prefers, the flat region
#' around it, and whether a ladder bound is doing the choosing rather than
#' the criterion.  The closing line gives the levels that lie in \emph{every}
#' flat region, the cell counts no criterion objects to.
#'
#' The criteria are not peers, and the table should be read with what each
#' returns in practice in mind (see \code{\link{resolution_profile}()}):
#' \code{reliability} is read at or above the range floor and usually
#' returns it, \code{elbow} is absent unless the points cluster, and
#' \code{cp}'s optimum is often the support ceiling.  \code{moran_z} is a
#' diagnostic of the cell means and not a count rule, so it is not in the
#' default set; named in \code{criteria} it gets a row, and it is left out
#' of the closing comparison.
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
#' distinguishable from noise (\code{reliability}), and where the
#' within-cluster sum of squares bends (\code{elbow}).  A field with no
#' single right resolution shows up here as disjoint bands, and the spread
#' between the picks is printed for the same reason.
#'
#' Nothing in the table is a decision procedure.  Each flat region is
#' routinely wide, the choice within it belongs to the analyst, and
#' \code{\link{plot.resolution_profile}()} draws the curves the bands were
#' read from.
#'
#' A profile whose ladder is one size (\code{obs_per_cell} set to one
#' number, or the default profile) has one level, and the size set it.  The
#' table then has one row per criterion that is finite at that level, and
#' the line under it says
#' that the count is set by the size and that no criterion chooses among
#' levels.  The result carries the size as \code{attr(, "obs_per_cell")}.
#'
#' @param object A \code{\link{resolution_profile}()}.
#' @param criteria Character vector, any of \code{"cp"},
#'   \code{"reliability"}, \code{"elbow"}, \code{"moran_z"}, \code{"loo"}.
#'   Default: every one of \code{"cp"}, \code{"reliability"},
#'   \code{"elbow"} and \code{"loo"} that is finite at some level
#'   (\code{loo} is a column of a profile built with \code{estimate}; a
#'   column of that name on a profile built without one is read only when
#'   named here).
#'   \code{"moran_z"}, a diagnostic, is read only when named here.
#'   \code{"size"} and \code{"auto"} are not criteria of the table: a size
#'   sets the count and has no flat region to compare.
#' @param tol Passed to \code{\link{select_resolution}()} for the flat
#'   region.  Default 0.02.
#' @param ... Ignored.
#' @return A data.frame of class \code{resolution_summary}, one row per
#'   criterion in the order given, with columns \code{criterion},
#'   \code{best}, \code{flat_min}, \code{flat_max}, \code{n_flat},
#'   \code{value} (the optimum itself, on that criterion's own scale and so
#'   not comparable across rows, which is why the print method leaves it out),
#'   \code{at_floor}, \code{at_ceiling}, \code{edge} (which bound an
#'   edge optimum sits on, as \code{\link{select_resolution}()} reports it;
#'   \code{NA} when interior) and \code{structured} (as
#'   \code{select_resolution()} reports it: \code{FALSE} for a \code{cp}
#'   that one cell matches).  Attributes: \code{bands} (a named list holding each
#'   criterion's flat region in full, since \code{flat_min} and
#'   \code{flat_max} are only its ends and the region can have holes in it),
#'   \code{common} (the levels in every flat region of the count criteria,
#'   \code{moran_z} left out; an integer vector that is
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
#'   # 600 m) on a 1 km square, with a nugget of 0.6 on a unit sill.
#'   set.seed(4)
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
  # `loo` is the leave-one-out error of the estimate a profile was built for.
  # A profile without one has no such column, and a column of that name a
  # caller added to it is read only when named (.loo_by_default()), so its
  # table is unchanged.
  selectable <- c("cp", "reliability", "elbow", "moran_z", "loo")
  # moran_z is a diagnostic of the cell means, not a count rule: its smallest
  # |z| usually sits below Cp's flat region, which left the default table's
  # "in every flat region" line empty.  It is read when named, and then kept
  # out of the closing comparison.
  count_rules <- setdiff(selectable, "moran_z")
  # A ladder of one size is one level: the size set the count, and no
  # criterion has levels to choose among.
  bounds <- attr(object, "bounds", exact = TRUE)
  one_size <- if (is.list(bounds) && length(bounds$obs_per_cell) == 1L)
    as.numeric(bounds$obs_per_cell) else NULL
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
    read <- if (.loo_by_default(object)) count_rules else setdiff(count_rules, "loo")
    criteria <- read[vapply(read, function(cn)
      cn %in% names(object) &&
        any(is.finite(suppressWarnings(as.numeric(object[[cn]])))), logical(1))]
  }
  bad <- setdiff(criteria, selectable)
  # What select_resolution(criterion = "auto") reads this profile by.
  est_value <- .profile_estimate(object)
  auto_is <- if (identical(est_value, "pooled")) "size"
             else if (identical(est_value, "mean")) "loo"
             else "cp"
  if (length(bad))
    stop("summary.resolution_profile(): unknown criteria: ", paste(bad, collapse = ", "),
         ". Choose from ", paste(selectable, collapse = ", "), ".",
         # "auto" names the rule that goes with the profile's estimate; said
         # with what it is on this profile.
         if ("auto" %in% bad)
           sprintf(paste0(" \"auto\" names no criterion of the table; on this profile ",
                          "select_resolution() takes it as \"%s\"."), auto_is)
         else "",
         # A size sets the count; it has no optimum and no flat region to put
         # beside the others'.
         if ("size" %in% bad || ("auto" %in% bad && auto_is == "size"))
           paste0(" A size sets the count and has no flat region to compare; ",
                  "select_resolution(criterion = \"size\") reads it.")
         else "",
         call. = FALSE)
  if (!length(criteria)) {
    # Two different faults reach here, and the wrong message sends the caller
    # off to fix the data when the object is what is wrong.
    if (!any(setdiff(selectable, "loo") %in% names(object)))
      stop("summary.resolution_profile(): `object` carries none of the ",
           "criterion columns (", paste(setdiff(selectable, "loo"), collapse = ", "),
           "); it is a column subset of a resolution profile. Pass the ",
           "profile resolution_profile() returned.", call. = FALSE)
    stop("summary.resolution_profile(): no count criterion is finite at any level. ",
         "cp needs a response, reliability a usable variogram; a geometry-only ",
         "profile carries elbow alone, ",
         "and elbow is NA when the WSS curve has no elbow.",
         if ("moran_z" %in% names(object) &&
             any(is.finite(suppressWarnings(as.numeric(object[["moran_z"]])))))
           " moran_z is finite, and is a diagnostic: name it in `criteria` to see it."
         else "",
         if (!is.null(one_size))
           sprintf(paste0(" The ladder is one size (obs_per_cell = %s); the count is set ",
                          "by the size, and select_resolution(criterion = \"size\") reads it."),
                   format(one_size))
         else if (is.list(bounds) && length(bounds$obs_per_cell))
           # Several sizes: the three criteria can all be NA for a reason the
           # sentences above do not give.
           sprintf(paste0(" The profile is on a ladder of sizes (obs_per_cell = %s)%s; ",
                          "select_resolution(criterion = \"size\") reads it."),
                   .got_text(bounds$obs_per_cell),
                   if (identical(attr(object, "cp_noise", exact = TRUE)$source, "none"))
                     paste0(", which fits no variogram unless `sac` is given, so it has no ",
                            "cp and no reliability")
                   else "")
         else "",
         call. = FALSE)
  }

  sels <- if (explicit)
    # Naming moran_z in `criteria` is the request select_resolution() asks for.
    lapply(criteria, function(cn) select_resolution(object, criterion = cn, tol = tol,
                                                     diagnostic = TRUE))
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
    structured = vapply(sels, function(s) as.logical(s$structured %||% NA)[1L], logical(1)),
    stringsAsFactors = FALSE)

  # The levels no criterion objects to.  Intersecting the regions rather than
  # the picks is the point: a single level chosen by one criterion says
  # nothing about whether the others can live with it.
  cnames <- vapply(sels, function(s) s$criterion, character(1))
  bands  <- stats::setNames(lapply(sels, function(s) s$flat), cnames)
  scored <- stats::setNames(lapply(sels, function(s) .scored_levels(s$values)), cnames)
  # Over the count criteria: a diagnostic has no say in which counts they
  # agree on (its own band when it is the only row).
  common <- Reduce(intersect, bands[if (any(cnames %in% count_rules)) cnames %in% count_rules
                                    else TRUE])
  lv <- as.integer(object$levels)
  # The whole ladder, not its ends: the print method needs it to tell a band
  # with a hole in it from a solid run, and `[` on the result keeps it, so a
  # subset can still be read honestly.
  out <- structure(out, class = c("resolution_summary", "data.frame"),
            common = sort(unique(as.integer(common))), bands = bands,
            scored = scored, ladder = lv,
            n_levels = nrow(object), tol = tol, variable = attr(object, "variable"))
  # The kind that built the cells, when the profile recorded one; the summary
  # of a k-means profile is the object it always was.
  part <- attr(object, "partition", exact = TRUE)
  if (is.list(part) && !is.null(part$kind))
    attr(out, "partition") <- part[.partition_fields()]
  # The one size of the ladder, for the print method's line; absent on every
  # other ladder.
  if (!is.null(one_size)) attr(out, "obs_per_cell") <- one_size
  out
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
  rungs_of <- function(cn) if (!is.na(cn) && !is.null(scored[[cn]])) scored[[cn]] else lad
  band_of <- function(i) {
    cn <- x$criterion[i]
    # A row index past the last row (`s[2, ]` of a one-row table) is a row of
    # NA, which seq() refuses.
    if (is.na(cn) || !is.finite(x$flat_min[i]) || !is.finite(x$flat_max[i])) return("none")
    b <- if (!is.null(bands[[cn]])) bands[[cn]] else seq(x$flat_min[i], x$flat_max[i])
    .band_label(b, rungs_of(cn))
  }
  cat(sprintf("Resolution picks: %d criteri%s over %d level%s (%d to %d cells)\n\n",
              nrow(x), if (nrow(x) == 1L) "on" else "a",
              n_lv, if (n_lv == 1L) "" else "s",
              min(lad), max(lad)))
  # The kind that built the cells, when the profile recorded one.
  part <- attr(x, "partition", exact = TRUE)
  if (is.list(part) && !is.null(part$kind))
    cat(sprintf("  cells: %s\n\n", .partition_label(part)))
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
  # One blank line between the table and the first note, whichever it is.
  gap <- FALSE
  open_notes <- function() if (!gap) { cat("\n"); gap <<- TRUE }
  # A ladder of one size has one level, and it is the size that set it: said
  # in place of "the optimum is the only level scored".
  one_size <- attr(x, "obs_per_cell", exact = TRUE)
  if (is.numeric(one_size) && length(one_size) == 1L) {
    open_notes()
    cat(strwrap(sprintf(paste0("the ladder is one size (obs_per_cell = %s); the count is ",
                               "set by the size, and no criterion chooses among levels."),
                        format(one_size)),
                width = 78, prefix = "  ", exdent = 2), sep = "\n")
  } else if (any(on_edge)) {
    # Wrapped, because with four criteria named one line reached 87
    # columns, which is the wrapping this block was moved here to avoid.
    open_notes()
    for (e in unique(edge[on_edge]))
      cat(strwrap(sprintf("%s: the optimum is %s.",
                          paste(x$criterion[on_edge & edge == e], collapse = ", "), e),
                  width = 78, prefix = "  ", exdent = 2), sep = "\n")
    cat("  There the bound is choosing, not the criterion.\n")
  }
  # A count one cell matches is not identified, whatever the table shows.
  if ("structured" %in% names(x) && any(x$structured %in% FALSE)) {
    open_notes()
    cat(strwrap(sprintf(paste0("%s: one cell for the whole layer does as well, so the ",
                               "response shows no spatial structure the cells can use ",
                               "and this count is not identified."),
                        paste(x$criterion[x$structured %in% FALSE], collapse = ", ")),
                width = 78, prefix = "  ", exdent = 2), sep = "\n")
  }
  # moran_z gets its row when asked for, and no part in the comparison: it is
  # a diagnostic of the cell means, and its band says nothing about which
  # counts represent the response.
  rule <- !is.na(x$criterion) & x$criterion != "moran_z"
  if (any(x$criterion %in% "moran_z")) {
    open_notes()
    cat(strwrap(paste0("moran_z is a diagnostic of the cell means, not a count rule",
                       if (sum(rule) > 1L) "; it is left out of the comparison below." else "."),
                width = 78, prefix = "  ", exdent = 2), sep = "\n")
  }
  # With one criterion there is nothing to compare: the spread is zero and the
  # intersection is that criterion's own band, so the closing block is skipped.
  if (sum(rule) > 1L) {
    crit <- x$criterion[rule]
    # Recomputed from the rows in hand rather than read off the attribute.
    # `[` carries `common` through a row subset unchanged, so a two-criterion
    # slice of a four-criterion table used to print the parent's verdict: on
    # 25 simulated fields that was wrong for 39 of 250 subsets, denying an
    # overlap the two rows above it plainly showed.
    cm <- if (!is.null(bands) && all(crit %in% names(bands)))
      sort(unique(Reduce(intersect, bands[crit])))
    else sort(unique(attr(x, "common")))
    cat("\n")
    rng <- range(x$best[rule])
    cat(sprintf("  picks span %d to %d cells (%.1fx)\n", rng[1L], rng[2L],
                rng[2L] / max(rng[1L], 1L)))
    if (length(cm))
      cat(sprintf("  in every flat region: %s (%d level%s)\n",
                  .band_label(cm, Reduce(intersect, lapply(crit, rungs_of))),
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
#' \code{resolution_profile} (read with \code{select_resolution()}: by
#' \code{cp} when it names no estimate, by \code{criterion = "auto"} when it
#' names one).  The count and where it came from are
#' returned so the caller can record both on its output.
#'
#' Whether the profile carries an \code{"estimate"} attribute is tested
#' first, and a profile that has one goes straight to
#' \code{select_resolution(x, criterion = "auto")}: a size for pooled
#' values, the leave-one-out minimum for plain means.  The test that some
#' criterion is finite, below, is not applied to it: a profile of one size
#' that fitted no variogram has no \code{cp}, \code{reliability} or
#' \code{elbow} and still has a count.  An error \code{select_resolution()}
#' raises on the way is passed on as it is.
#'
#' A profile without an estimate is read by \code{cp} or not at all, with
#' one exception: a
#' geometry-only profile whose one finite criterion is \code{elbow} is read
#' by that, since it is the only thing such a profile can mean.  A profile
#' whose \code{cp} is \code{NA} at every level while \code{reliability} (or
#' \code{elbow} beside a response) is finite stops with an error of class
#' \code{spatialkit_no_default_criterion} naming what is available, as
#' \code{select_resolution()} on the same object does.  \code{moran_z} is a
#' diagnostic and is never read here.
#'
#' A count that a bound of the ladder chose is raised as a warning of class
#' \code{spatialkit_bound_chosen} (fields \code{edge}, \code{criterion},
#' \code{n}), and a count that one cell matches as one of class
#' \code{spatialkit_no_structure}: the caller of the tessellation functions
#' is then told by an R condition, where there used to be a log line and a
#' free-text attribute.
#'
#' @param x \code{NULL}, a number, a numeric vector, a
#'   \code{resolution_selection} or a \code{resolution_profile}.
#' @param arg,caller Names for the messages.
#' @return \code{NULL} when \code{x} is \code{NULL}; otherwise a list with
#'   \code{n} (a single number), \code{from} (a character description, or
#'   \code{NULL} when \code{x} was a plain number), \code{seeds} (the
#'   centres of the partition the profile scored at \code{n} cells, from
#'   \code{.profile_seeds()}, or \code{NULL} when \code{x} carries none) and,
#'   for a selection or a profile, \code{criterion}, \code{edge},
#'   \code{at_floor}, \code{at_ceiling} and \code{structured} as
#'   \code{select_resolution()} reports them (\code{NULL} otherwise), and
#'   \code{count_rule}, the selection's record of how the count was set
#'   (\code{NULL} for a number, and for a selection that carries none).
#' @keywords internal
#' @noRd
.resolve_cell_count <- function(x, arg, caller) {
  if (is.null(x)) return(NULL)
  from <- NULL
  seeds <- NULL
  sel <- NULL
  edge_note <- function(sel) {
    e <- sel$edge %||% NA_character_
    if (is.na(e)) "" else sprintf(", an optimum at %s", e)
  }
  # Whether a selection's count was set by a size.  The size is then part of
  # the rule, and the provenance note names it with the criterion.
  by_size <- function(sel) {
    asked <- sel$obs_per_cell_asked %||% NA_real_
    identical(sel$criterion, "size") && is.numeric(asked) && length(asked) == 1L &&
      is.finite(asked)
  }
  # The call that made a selection, as the provenance note writes it.
  call_text <- function(sel) {
    if (by_size(sel))
      sprintf("select_resolution(criterion = \"size\", obs_per_cell = %s)",
              format(sel$obs_per_cell_asked))
    else sprintf("select_resolution(criterion = \"%s\")", sel$criterion)
  }
  if (inherits(x, "resolution_selection")) {
    n <- x$best
    sel <- x
    # For every criterion but a size this is the line it always was, whatever
    # the selection holds as its `criterion`.
    from <- if (by_size(x)) paste0(call_text(x), edge_note(x))
            else sprintf("select_resolution(criterion = \"%s\")%s", x$criterion, edge_note(x))
    seeds <- x$seeds
  } else if (inherits(x, "resolution_profile") &&
             !is.null(attr(x, "estimate", exact = TRUE))) {
    # A profile that names an estimate is read by the rule that goes with
    # it: a size for pooled values, the leave-one-out minimum for plain
    # means (criterion = "auto").  Tested first, before any criterion is
    # asked to be finite: a profile of one size that fitted no variogram has
    # no `cp`, `reliability` or `elbow` and still has a count, and so does a
    # pooled profile of a constant response.  An error select_resolution()
    # raises on the way is passed on as it is.
    sel <- select_resolution(x, criterion = "auto")
    n <- sel$best
    seeds <- sel$seeds
    from <- sprintf("resolution_profile() read with %s%s", call_text(sel), edge_note(sel))
    .log_info("%s(): `%s` is a resolution profile; read with %s: %d cells.",
              caller, arg, call_text(sel), n)
  } else if (inherits(x, "resolution_profile")) {
    # A profile without an estimate is read by cp, as the default criterion
    # of select_resolution() reads one.  It
    # used to fall through cp, reliability, elbow, moran_z to the first that
    # was finite somewhere, with an INFO log line the console does not show:
    # a geometry-only profile given a variogram was read by reliability and
    # handed on 13 cells where the response profile of the same layer gives
    # 100, and a response profile with no cp was read by moran_z.
    # select_resolution() on either object is an error, and so is this.
    finite_in <- function(cn) cn %in% names(x) &&
      any(is.finite(suppressWarnings(as.numeric(x[[cn]]))))
    usable <- Filter(finite_in, c("cp", "reliability", "elbow"))
    geometry_only <- is.na(attr(x, "variable") %||% NA_character_)
    # The sizes of a ladder `obs_per_cell` set.  Such a profile has a count,
    # the one a size gives, though none of the three criteria may be finite
    # on it: a ladder of sizes fits no variogram unless `sac` is given, and
    # one size is one level, which has no WSS curve.
    ladder_sizes <- if (is.list(attr(x, "bounds", exact = TRUE)))
      attr(x, "bounds", exact = TRUE)$obs_per_cell else NULL
    # A geometry-only profile of a layer with no cluster structure lands here
    # too: its elbow is NA rather than the sqrt(first x last level) the
    # linear chord rule used to hand on as if the data had chosen it.
    if (!length(usable))
      stop(sprintf(paste0("%s(): `%s` is a resolution profile with no criterion finite at ",
                          "any level, so no cell count can be read off it.%s%s"),
                   caller, arg,
                   if (length(ladder_sizes))
                     sprintf(paste0(" It is on a ladder of sizes (obs_per_cell = %s) and was ",
                                    "built without an `estimate`, so it is read by `cp`, or by ",
                                    "`elbow` when it has no `response_var`, and neither is ",
                                    "finite on it%s. A size reads it: pass ",
                                    "select_resolution(<profile>, criterion = \"size\") as ",
                                    "`%s`, or pass a number."),
                             .got_text(ladder_sizes),
                             if (!geometry_only &&
                                 identical(attr(x, "cp_noise", exact = TRUE)$source, "none"))
                               paste0(" (a ladder of sizes fits no variogram unless `sac` is ",
                                      "given, and `cp` takes its noise variance from one)")
                             else "",
                             arg)
                   else if (geometry_only)
                     paste0(" It is geometry-only, and its WSS curve has no elbow: the ",
                            "points have no cluster structure to choose a count. Pass a ",
                            "number, or profile with a `response_var`.")
                   else "",
                   if (finite_in("moran_z"))
                     paste0(" (`moran_z` is finite on it, but it is a diagnostic of the ",
                            "cell means and is not read as a count.)")
                   else ""), call. = FALSE)
    # The one fall-through kept: a geometry-only profile whose only finite
    # criterion is the elbow can mean nothing else.
    crit <- if ("cp" %in% usable) "cp"
            else if (geometry_only && identical(usable, "elbow")) "elbow"
            else NULL
    if (is.null(crit))
      stop(errorCondition(
        sprintf(paste0("%s(): `%s` is a resolution profile whose `cp`, the criterion a ",
                       "profile is read by, is NA at every level%s, so no cell count is ",
                       "read off it. Finite on it: %s. To use one, pass ",
                       "select_resolution(<profile>, criterion = \"%s\") as `%s`; or ",
                       "pass a number.%s"),
                caller, arg,
                if (geometry_only) " (the profile has no `response_var`)" else "",
                paste(usable, collapse = ", "), usable[[1L]], arg,
                # A ladder of sizes has one more reader.
                if (length(ladder_sizes))
                  sprintf(paste0(" The profile is on a ladder of sizes (obs_per_cell = %s): ",
                                 "criterion = \"size\" reads it by a size."),
                          .got_text(ladder_sizes))
                else ""),
        class = "spatialkit_no_default_criterion", call = NULL))
    sel <- select_resolution(x, criterion = crit)
    n <- sel$best
    seeds <- sel$seeds
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
    stop(sprintf("%s(): `%s` must resolve to a positive number of cells; got %s.%s",
                 caller, arg, paste(format(n), collapse = ", "),
                 # determine_optimal_levels(on_no_elbow = "na") hands this on.
                 if (isFALSE(attr(x, "structured", exact = TRUE)))
                   paste0(" It is what determine_optimal_levels(on_no_elbow = \"na\") ",
                          "returns when the WSS curve has no elbow: the points name no ",
                          "count. Pass a number.")
                 else ""), call. = FALSE)
  # Seeds only when they are the partition at exactly this count; a selection
  # edited by hand, or one whose level had no successful k-means, has none.
  if (!inherits(seeds, "sf") || nrow(seeds) != n) seeds <- NULL
  # The count is about to become cells.  When a bound of the ladder chose
  # it, or one cell does as well, the caller is told by an R warning with a
  # class: a count that sat at the ceiling used to be marked by nothing but
  # the text of attr(, "n_from").
  if (!is.null(sel)) {
    e <- sel$edge %||% NA_character_
    if (length(e) == 1L && !is.na(e))
      .warn_classed("spatialkit_bound_chosen",
                    paste0("%s(): `%s` is %d cells, the optimum of `%s`, and that optimum is ",
                           "%s: the bound chose the count, not the criterion.%s Pass a ",
                           "number as `%s` to set the count yourself."),
                    caller, arg, as.integer(n), sel$criterion %||% "the criterion", e,
                    .edge_hint(e), arg,
                    data = list(edge = e, criterion = sel$criterion, n = as.integer(n)))
    # The notes that say why the count is not identified: every note of the
    # selection, as before, but the ones select_resolution() writes about
    # something else (which rule set the count, a size the profile held, the
    # end of the ladder, which estimate the criterion scores, a ladder of one
    # size).  Those carry one of five names and are left out here.  A note
    # without a name, or with any other, is pasted as it always was, so a
    # selection whose notes the caller wrote is told in the caller's words.
    why <- sel$note
    other <- names(why) %in% .NOTES_NOT_STRUCTURE
    if (any(other)) {
      why <- why[!other]
      if (!length(why)) why <- NULL
    }
    if (isFALSE(sel$structured) && by_size(sel))
      # A size set this count and no criterion, so the count stands; what
      # the caller is told is that the cells will say little.  A selection
      # that names "size" and no size asked for is not one this function
      # wrote, and is told in the present words.
      .warn_classed("spatialkit_no_structure",
                    paste0("%s(): `%s` is %d cells, set by obs_per_cell = %s, and %s The ",
                           "count stands, since a size set it and no criterion; the cell ",
                           "values will say little that one mean for the layer does not."),
                    caller, arg, as.integer(n), format(sel$obs_per_cell_asked %||% NA_real_),
                    paste(why %||% paste0("one mean for the whole layer predicts a left-out ",
                                          "observation as well as cells of this size do."),
                          collapse = " "),
                    data = list(criterion = sel$criterion, n = as.integer(n),
                                one_cell = sel$one_cell))
    else if (isFALSE(sel$structured))
      .warn_classed("spatialkit_no_structure",
                    paste0("%s(): `%s` is %d cells, the optimum of `%s`, but %s Pass a ",
                           "number as `%s` if a count is wanted on other grounds."),
                    caller, arg, as.integer(n), sel$criterion %||% "the criterion",
                    paste(why %||% paste0("one cell for the whole layer does as well, ",
                                          "so the count is not identified."),
                          collapse = " "),
                    arg,
                    data = list(criterion = sel$criterion, n = as.integer(n),
                                one_cell = sel$one_cell))
  }
  list(n = as.numeric(n), from = from, seeds = seeds,
       criterion = sel$criterion, edge = sel$edge, at_floor = sel$at_floor,
       at_ceiling = sel$at_ceiling, structured = sel$structured,
       # How the count was set, for the functions that turn it into cells:
       # NULL for a number, and for a count a present criterion read off a
       # profile without an estimate.
       count_rule = sel$count_rule)
}
