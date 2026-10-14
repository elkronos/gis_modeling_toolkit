# =============================================================================
# The partition engine: the centres whose Voronoi cells are the regions.
#
# Two kinds of partition, each from one of two starting rules:
#
#   kind "kmeans"   minimises the sum over the points of the squared distance
#                   to the nearest centre.  Its centres follow the square root
#                   of the point density, so a cell in a dense place holds
#                   several times the points of a cell in a sparse place.
#   kind "density"  minimises the same sum with each point weighted by
#                   (local density)^power.  With power 1 the centres follow
#                   the point density itself.  Once there are more cells than
#                   dense places (clusters of points) the cells hold about
#                   the same number of points each; with fewer, a cell takes
#                   in whole dense places, large and small, and the counts
#                   are about as uneven as those of the kind "kmeans".
#
#   start "kmeans++"  the first centre is a point drawn at random, each later
#                     one a point drawn in proportion to its (weighted)
#                     squared distance from the nearest centre so far
#   start "random"    k distinct locations drawn at random
#
# .partition_centres() is the one entry point.  resolution_profile(),
# determine_optimal_levels(), get_voronoi_seeds() and voronoi_seeds_kmeans()
# all take their centres from it.
# =============================================================================

# The starting rules, in the order the arguments list them.
.partition_inits <- function() c("kmeans++", "random")

# Two constants of the density kind, both fractions of the extent of the
# points (the diagonal of their bounding box), so that neither depends on the
# unit the coordinates are in.
#
# The neighbour distance a local density is read from is not taken below
# 1e-6 of the extent.  Repeated locations never reach it: the distance is
# measured between distinct locations (see .density_weights()).  It is for
# locations that are distinct only in their last digits.
#
# Lloyd's iteration stops when no coordinate of a centre changes by 1e-12 of
# the extent or more.  A step in which a point changes cell moves a centre by
# far more than that, so in practice it stops when the cells stop changing.
.DENSITY_FLOOR_REL <- 1e-6
.DENSITY_TOL_REL   <- 1e-12

# Starts of the density kind when the caller sets no number: the best of
# three is kept.  The kind "kmeans" keeps the number each function has always
# run.
.NSTART_DENSITY <- 3L


# The four fields of a partition record that say what built the regions.
.partition_fields <- function() c("kind", "init", "density_power", "density_k")


#' Match `init` with a message that names the function
#'
#' \code{match.arg()} on a value that is not a string stops with "'arg' must
#' be NULL or a character vector", which names neither the argument nor the
#' function.
#'
#' @param init The value the caller passed.
#' @param caller Function name for the message.
#' @param null_ok Whether \code{NULL} is a value of its own (the seeding
#'   functions' default).
#' @return \code{"kmeans++"} or \code{"random"}, or \code{NULL} when
#'   \code{init} is \code{NULL} and \code{null_ok}.
#' @keywords internal
#' @noRd
.match_init <- function(init, caller, null_ok = FALSE) {
  if (is.null(init) && isTRUE(null_ok)) return(NULL)
  choices <- .partition_inits()
  hit <- if (is.character(init) && length(init) == 1L && !is.na(init))
    pmatch(init, choices) else NA_integer_
  if (is.na(hit))
    stop(sprintf("%s(): `init` should be one of %s%s; got %s.", caller,
                 paste(sprintf("\"%s\"", choices), collapse = ", "),
                 if (isTRUE(null_ok)) ", or NULL" else "",
                 if (is.character(init) && length(init) == 1L && !is.na(init))
                   sprintf("\"%s\"", init)
                 else if (length(init) != 1L)
                   sprintf("%s of length %d", class(init)[1L], length(init))
                 else if (is.character(init)) "NA"
                 else class(init)[1L]),
         call. = FALSE)
  choices[hit]
}


#' Validate the partition arguments of a public function
#'
#' \code{partition} and \code{init} arrive already matched (\code{partition}
#' by \code{match.arg()} in the caller, as the neighbouring arguments are;
#' \code{init} by \code{.match_init()}).
#' \code{density_power} and \code{density_k} are checked here whatever the
#' kind, so that a value that makes no sense is refused and not ignored, and
#' a value set under \code{partition = "kmeans"}, where neither is read, is
#' reported: here, or by the caller through \code{.warn_density_unused()}
#' once it knows that a k-means run is what the call gets (a seeding function
#' handed a selection returns scored centres, which may be density ones).
#'
#' @param partition,init The kind and the starting rule in force.
#' @param density_power,density_k As the caller received them.
#' @param caller Function name for the messages.
#' @param warn_unused Whether to report here a density setting that
#'   \code{partition = "kmeans"} does not read.
#' @return A list with \code{kind}, \code{init}, \code{density_power} and
#'   \code{density_k} (the last two \code{NA} for \code{kind = "kmeans"}),
#'   the record the result objects carry, and two elements for the caller:
#'   \code{asked} (the two density settings as given) and \code{unused} (the
#'   names of those set to something other than the default under
#'   \code{kind = "kmeans"}).
#' @keywords internal
#' @noRd
.check_partition_args <- function(partition, init, density_power, density_k, caller,
                                  warn_unused = TRUE) {
  .check_scalar(density_power, "density_power", caller, min = 0, max = 2,
                what = "a single number from 0 to 2")
  .check_scalar(density_k, "density_k", caller, min = 1,
                max = .Machine$integer.max,
                what = "a single positive whole number of neighbours")
  if (density_k != floor(density_k))
    stop(sprintf("%s(): `density_k` must be a single positive whole number of neighbours; got %s.",
                 caller, format(density_k)), call. = FALSE)
  asked <- list(density_power = as.numeric(density_power), density_k = as.integer(density_k))
  if (identical(partition, "kmeans")) {
    part <- list(kind = "kmeans", init = init,
                 density_power = NA_real_, density_k = NA_integer_, asked = asked,
                 unused = c(if (!isTRUE(all.equal(density_power, 1))) "density_power",
                            if (density_k != 8) "density_k"))
    if (isTRUE(warn_unused)) .warn_density_unused(part, caller)
    return(part)
  }
  list(kind = "density", init = init,
       density_power = asked$density_power, density_k = asked$density_k,
       asked = asked, unused = character(0))
}


#' Say that a density setting was given where no density partition is built
#'
#' @param part A \code{.check_partition_args()} result.
#' @param caller Function name for the message.
#' @return \code{NULL}, invisibly; a warning when \code{part$unused} names a
#'   setting.
#' @keywords internal
#' @noRd
.warn_density_unused <- function(part, caller) {
  set <- part$unused
  if (length(set))
    .warn_and_log(paste0("%s(): %s %s read by partition = \"density\" only, and ",
                         "partition = \"kmeans\" ignores %s."),
                  caller, paste(sprintf("`%s`", set), collapse = " and "),
                  if (length(set) > 1L) "are" else "is",
                  if (length(set) > 1L) "them" else "it")
  invisible(NULL)
}


#' Check a count of starts or an iteration limit before the engine reads it
#'
#' @param x The value.
#' @param arg,caller Names for the message.
#' @param what What the number counts, for the message.
#' @param null_ok Whether the argument also takes \code{NULL}, which the
#'   caller has dealt with; the message then names it.
#' @return \code{x} as an integer.
#' @keywords internal
#' @noRd
.check_whole_count <- function(x, arg, caller, what, null_ok = FALSE) {
  what <- sprintf("%sa single positive whole number of %s",
                  if (isTRUE(null_ok)) "NULL or " else "", what)
  .check_scalar(x, arg, caller, min = 1, max = .Machine$integer.max, what = what)
  if (x != floor(x))
    stop(sprintf("%s(): `%s` must be %s; got %s.", caller, arg, what, format(x)),
         call. = FALSE)
  as.integer(x)
}


#' The partition record a result object carries, or the default it stands for
#'
#' A result built by a kind or a starting rule other than its function's own
#' default carries \code{attr(, "partition")} (a tessellation carries it as
#' \code{params$partition}, a \code{select_resolution()} result as
#' \code{$partition}).  A profile, a selection or a layer of k-means seeds
#' with no record was built by plain k-means.
#'
#' An object that is not a k-means result has no kind: seeds of
#' \code{method = "random"} or \code{"provided"}, a grid or a triangulation,
#' and anything else.  So has a Voronoi tessellation with no record, because
#' its cells do not say whether the points they were built on were k-means
#' seeds or the observations themselves.
#'
#' @param x A resolution profile, a selection, a layer of seeds or a
#'   tessellation.
#' @return A list with \code{kind}, \code{init}, \code{density_power} and
#'   \code{density_k}.  For a k-means result with no record \code{kind} is
#'   \code{"kmeans"} and \code{init} is \code{NA}: the starting rule is then
#'   the default of the function that made it.  For an object that is not a
#'   k-means result, or does not say, every element is \code{NA}.
#' @keywords internal
#' @noRd
.partition_record <- function(x) {
  rec <- attr(x, "partition", exact = TRUE)
  if (is.null(rec) && is.list(x) && !inherits(x, "data.frame"))
    rec <- if (is.list(x[["params"]])) x[["params"]][["partition"]] else x[["partition"]]
  if (is.list(rec) && !is.null(rec$kind))
    return(rec[.partition_fields()])
  kmeans_made <- inherits(x, c("resolution_profile", "resolution_selection")) ||
    (inherits(x, "sf") && "method" %in% names(x) && nrow(x) > 0L &&
       all(x[["method"]] %in% "kmeans"))
  list(kind = if (kmeans_made) "kmeans" else NA_character_, init = NA_character_,
       density_power = NA_real_, density_k = NA_integer_)
}


#' A partition's kind and starting rule in words
#'
#' For the print method of a profile and the log lines of the seeding
#' functions.
#'
#' @param rec A partition record.
#' @param nstart Number of starts, or \code{NULL} to leave it out.
#' @return One string, such as \code{"25 density-weighted k-means++ starts
#'   (power 1, 8 neighbours)"}.
#' @keywords internal
#' @noRd
.partition_label <- function(rec, nstart = NULL) {
  rule <- if (identical(rec$init, "random")) "random" else "k-means++"
  paste0(if (is.null(nstart)) "" else paste0(format(nstart), " "),
         if (identical(rec$kind, "density")) "density-weighted " else "",
         rule, " start", if (is.null(nstart) || !identical(as.numeric(nstart), 1)) "s" else "",
         if (identical(rec$kind, "density"))
           sprintf(" (power %s, %d neighbour%s)", format(rec$density_power),
                   as.integer(rec$density_k),
                   if (identical(as.integer(rec$density_k), 1L)) "" else "s")
         else "")
}


#' Say that some centres of a density partition hold no point
#'
#' A centre of the density kind that loses its points stays where it is
#' (see \code{.density_kmeans()}), so a partition can come back with fewer
#' occupied cells than centres, and the start that is kept can be such a
#' one: the starts are ranked by the weighted sum alone.  It happens at high
#' counts, from about a fifth of the points up, and mostly from random
#' starts (the figures are in the "Conditions" list of
#' \code{resolution_profile()}).  It is said when it happens: the count
#' asked for is then not the count of regions that hold data.  The kind
#' \code{"kmeans"} never returns such a centre.
#'
#' @param sizes Points per centre.
#' @param caller Function name for the message.
#' @param what What the centres are called there (\code{"seeds"}).
#' @return \code{NULL}, invisibly; a warning of class
#'   \code{spatialkit_empty_cells} when a size is 0.
#' @keywords internal
#' @noRd
.warn_empty_seeds <- function(sizes, caller, what = "seeds") {
  empty <- sum(sizes == 0L)
  if (empty > 0L)
    .warn_classed("spatialkit_empty_cells",
                  paste0("%s(): %d of the %d %s hold no point (a centre of ",
                         "partition = \"density\" that loses its points stays where it ",
                         "is), so %d of their Voronoi cells will hold a point."),
                  caller, empty, length(sizes), what, length(sizes) - empty,
                  data = list(level = length(sizes), empty = as.integer(empty)))
  invisible(NULL)
}


#' Say that the start that was kept had not converged at the iteration limit
#'
#' The centres returned are then those of the last step.  A longer run moves
#' them further and can move points between cells, so the partition is close
#' to, and not at, the local optimum the start leads to.
#'
#' @param caller Function name for the message.
#' @param level The count, or the counts, of cells at which it happened.
#' @param iter Iterations the kept start took at each of them.
#' @param iter_max The limit.
#' @param kind The kind of partition.
#' @param arg The argument of \code{caller} that raises the limit.
#' @param of Number of levels the caller built, when it built more than one.
#' @return \code{NULL}, invisibly, after a warning of class
#'   \code{spatialkit_not_converged} that carries \code{level}, \code{iter}
#'   and \code{iter_max}.
#' @keywords internal
#' @noRd
.warn_not_converged <- function(caller, level, iter, iter_max, kind, arg, of = NULL) {
  level <- as.integer(level)
  one   <- identical(as.integer(iter_max), 1L)
  steps <- if (identical(kind, "density"))
    paste(if (one) "step" else "steps", "of Lloyd's iteration")
  else paste(if (one) "iteration" else "iterations", "of k-means")
  shown <- if (length(level) > 6L)
    paste0(paste(utils::head(level, 6L), collapse = ", "), ", ...")
  else paste(level, collapse = ", ")
  where <- if (is.null(of)) sprintf("at %s cells", shown)
           else sprintf("at %d of the %d levels (%s cells)", length(level), as.integer(of), shown)
  .warn_classed("spatialkit_not_converged",
                paste0("%s(): %s the start that was kept had not converged after %d %s ",
                       "(partition = \"%s\"). The centres are those of the last step; a ",
                       "longer run can still move points between cells. Raise `%s` to ",
                       "let it finish."),
                caller, where, as.integer(iter_max), steps, kind, arg,
                data = list(level = level, iter = as.integer(iter),
                            iter_max = as.integer(iter_max)))
  invisible(NULL)
}


#' The centres of a partition of the points into k regions
#'
#' The one place the package computes region centres.
#'
#' \strong{Kind \code{"kmeans"}.}  With \code{init = "kmeans++"} this is
#' \code{.kmeans_best()} as it always was: the best of \code{nstart}
#' k-means++ starts, each one run of \code{stats::kmeans()} from its seeds.
#' With \code{init = "random"} each start is \code{k} distinct locations
#' drawn at random, as \code{stats::kmeans(centers = k)} draws them;
#' \code{one_call = TRUE} hands the whole job to one
#' \code{stats::kmeans(centers = k, nstart = nstart)} call, which is what the
#' seeding functions have always run, warnings and errors included.
#'
#' \strong{Kind \code{"density"}.}  See \code{.density_kmeans()}.
#'
#' @param xy Numeric matrix of finite coordinates; the first two columns are
#'   used.
#' @param k Number of centres.
#' @param partition \code{"kmeans"} or \code{"density"}.
#' @param init \code{"kmeans++"} or \code{"random"}.
#' @param nstart Number of starts; the best is kept.
#' @param iter.max Iteration limit of each start: \code{stats::kmeans()}'s
#'   for the kind \code{"kmeans"}, Lloyd steps for the kind \code{"density"}.
#' @param density_power,density_k The weight's power and the number of
#'   neighbours the local density is read from.  Kind \code{"density"} only.
#' @param sort_rows Kind \code{"kmeans"} only: put the rows in coordinate
#'   order before the starts are drawn, so that the result does not depend on
#'   the order they come in.  \code{FALSE} where the caller has already put
#'   them in an order of its own, and where the rows are clustered as given.
#'   The kind \code{"density"} always works in coordinate order.
#' @param one_call See above.
#' @param use_fnn Whether the kind \code{"density"} may use \pkg{FNN} to
#'   shortlist neighbours.  The answer is the same without it.
#' @return \code{NULL} when no start produced a partition; otherwise a list
#'   with \code{centers} (a \code{k x 2} matrix), \code{cluster} (the centre
#'   of each row of \code{xy}, in the order given), \code{size} (points per
#'   centre), \code{withinss} (per cell, the unweighted sum of squared
#'   distances of its points from their mean), \code{wss} (their total),
#'   \code{objective} (what the starts minimised and were ranked by: the same
#'   total for the kind \code{"kmeans"}, the weighted sum for the kind
#'   \code{"density"}), \code{spread} (\code{(max - min) / min} of the
#'   objective over the starts that ran; \code{NA} under \code{one_call}),
#'   \code{n_ok} (starts that ran), \code{iter} (iterations of the best
#'   start; \code{stats::kmeans()} reports one more than the limit for a fit
#'   it cut off, which is held to the limit here except under
#'   \code{one_call}, whose record is the one the seeding functions have
#'   always returned), \code{converged} (whether that start stopped before
#'   the limit:
#'   for the kind \code{"kmeans"}, whether \code{stats::kmeans()} returned
#'   it with no fault) and \code{km} (the \code{stats::kmeans()} fit;
#'   \code{NULL} for the kind \code{"density"}).
#' @keywords internal
#' @noRd
.partition_centres <- function(xy, k, partition = "kmeans", init = "kmeans++",
                               nstart = 25L, iter.max = 50L,
                               density_power = 1, density_k = 8L,
                               sort_rows = FALSE, one_call = FALSE,
                               use_fnn = requireNamespace("FNN", quietly = TRUE)) {
  if (identical(partition, "density"))
    return(.density_kmeans(xy, k, power = density_power, k_nn = density_k,
                           init = init, nstart = nstart, iter.max = iter.max,
                           use_fnn = use_fnn))
  o <- if (isTRUE(sort_rows)) order(xy[, 1L], xy[, 2L]) else NULL
  xs <- if (is.null(o)) xy else xy[o, , drop = FALSE]
  if (identical(init, "random") && isTRUE(one_call)) {
    km <- stats::kmeans(xs, centers = k, iter.max = iter.max, nstart = nstart)
    spread <- NA_real_; n_ok <- NA_integer_
  } else {
    kb <- .kmeans_best(xs, k, nstart = nstart, iter.max = iter.max, init = init)
    if (is.null(kb)) return(NULL)
    km <- kb$km; spread <- kb$spread; n_ok <- kb$n_ok
  }
  cluster <- km$cluster
  if (!is.null(o)) { cluster <- integer(length(o)); cluster[o] <- km$cluster }
  list(centers = unname(km$centers[, 1:2, drop = FALSE]), cluster = cluster,
       size = as.integer(km$size), withinss = as.numeric(km$withinss),
       wss = km$tot.withinss, objective = km$tot.withinss,
       spread = spread, n_ok = n_ok,
       iter = if (identical(init, "random") && isTRUE(one_call)) as.integer(km$iter)
              else min(as.integer(km$iter), as.integer(iter.max)),
       # `ifault` is 2 for a fit that ran out of iterations and 4 for one
       # whose quick-transfer stage did; stats::kmeans() warns of both.
       converged = is.null(km$ifault) || identical(as.integer(km$ifault), 0L),
       km = km)
}


#' Starting centres drawn as stats::kmeans(centers = k) draws them
#'
#' \code{k} distinct locations: rows of the distinct rows \code{cn} drawn
#' without replacement.  With a single start \code{stats::kmeans()} first
#' draws among all the rows and falls back to the distinct ones only when two
#' of the rows drawn coincide; that is followed here, so the best of
#' \code{nstart} starts made this way is the fit
#' \code{stats::kmeans(xy, centers = k, nstart = nstart)} returns.
#'
#' @param xy Numeric matrix of coordinates.
#' @param cn Its distinct rows, \code{unique(xy)}.
#' @param k Number of centres.
#' @param nstart Number of starts the caller runs.
#' @return A \code{k}-row matrix of starting centres.
#' @keywords internal
#' @noRd
.random_centers <- function(xy, cn, k, nstart) {
  if (nstart == 1L) {
    cen <- xy[sample.int(nrow(xy), k), , drop = FALSE]
    if (!any(duplicated(cen))) return(cen)
  }
  cn[sample.int(nrow(cn), k), , drop = FALSE]
}


# -----------------------------------------------------------------------------
# The density kind
# -----------------------------------------------------------------------------

#' The distinct locations of a set of points, in coordinate order
#'
#' The density kind works on locations: a place visited thirty times is one
#' location with thirty points.  The order (first coordinate, then second) is
#' the same whatever order the rows come in, which is what makes the kind
#' independent of it; rows at one location are interchangeable.
#'
#' @param xy Two-column matrix of finite coordinates.
#' @return A list: \code{order} (the permutation that sorts the rows),
#'   \code{u} (the distinct locations, sorted), \code{m} (points at each) and
#'   \code{gid} (the location of each sorted row).
#' @keywords internal
#' @noRd
.partition_locations <- function(xy) {
  n  <- nrow(xy)
  o  <- order(xy[, 1L], xy[, 2L])
  xs <- xy[o, , drop = FALSE]
  first <- c(TRUE, xs[-1L, 1L] != xs[-n, 1L] | xs[-1L, 2L] != xs[-n, 2L])
  gid <- cumsum(first)
  list(order = o, u = unname(xs[first, , drop = FALSE]), m = tabulate(gid), gid = gid)
}


#' For chosen rows, the k nearest other locations by exact search
#'
#' Every squared distance is \code{(x - x')^2 + (y - y')^2} in R's own
#' arithmetic, and ties go to the lower row number.  This is the definition
#' the kd-tree path of \code{.knn_locations()} has to agree with.
#'
#' @param x,y Coordinates of the locations.
#' @param rows The rows to search for.
#' @param k Number of neighbours, below \code{length(x)}.
#' @return An integer matrix, \code{length(rows) x k}.
#' @keywords internal
#' @noRd
.knn_exact_rows <- function(x, y, rows, k) {
  n <- length(x)
  out <- matrix(0L, length(rows), k)
  if (!length(rows)) return(out)
  # Blocks of rows, so that the block-by-n matrix of distances stays small.
  per <- max(1L, as.integer(2e6 %/% n))
  for (from in seq.int(1L, length(rows), by = per)) {
    at <- from:min(from + per - 1L, length(rows))
    r  <- rows[at]
    d2 <- outer(x[r], x, "-")^2 + outer(y[r], y, "-")^2
    d2[cbind(seq_along(r), r)] <- Inf            # a location is not its own neighbour
    # order() leaves ties in their original order: the lower row number first.
    out[at, ] <- t(apply(d2, 1L, function(v) order(v)[seq_len(k)]))
  }
  out
}


#' For rows whose neighbours tie, the k nearest other locations from a kd-tree
#'
#' The answer of \code{.knn_exact_rows()} without its comparison of each row
#' with every location.  \pkg{FNN} shortlists more neighbours than are
#' wanted, R recomputes their distances, and the \code{k} nearest are taken
#' with the same rule: by distance, ties to the lower row number.  A
#' shortlist is long enough when its furthest member is further than the
#' \code{k}-th, by more than one part in a thousand million: every location
#' as near as the \code{k}-th is then on it.  A row whose shortlist ends
#' inside the tie gets one twice as long, up to all the other locations.  On
#' a lattice, or on coordinates rounded to the metre, most rows tie and a
#' short list settles each, so the time grows with the number of locations
#' and not with its square.
#'
#' @param u Two-column matrix of distinct locations.
#' @param rows The rows to search for.
#' @param k Number of neighbours, below \code{nrow(u)}.
#' @return An integer matrix, \code{length(rows) x k}.
#' @keywords internal
#' @noRd
.knn_tied_rows <- function(u, rows, k) {
  n <- nrow(u); x <- u[, 1L]; y <- u[, 2L]
  out  <- matrix(0L, length(rows), k)
  todo <- seq_along(rows)
  q <- min(n - 1L, 2L * (k + 1L))
  while (length(todo)) {
    r  <- rows[todo]; nr <- length(r)
    # The location itself comes back as its own nearest neighbour.
    nn <- matrix(FNN::get.knnx(u, u[r, , drop = FALSE], k = q + 1L)$nn.index, nr, q + 1L)
    self <- nn == r
    none <- which(rowSums(self) != 1L)
    if (length(none)) { self[none, ] <- FALSE; self[none, q + 1L] <- TRUE }
    cand <- matrix(t(nn)[!t(self)], nr, q, byrow = TRUE)
    d2   <- matrix((x[cand] - x[r])^2 + (y[cand] - y[r])^2, nr, q)
    # Each row's candidates by distance, then by row number: column i of
    # `ord` is row i's, nearest first.
    ord  <- matrix(order(rep.int(seq_len(nr), q), d2, cand), nrow = q)
    pick <- ord[seq_len(k), , drop = FALSE]
    dk   <- d2[pick[k, ]]
    far  <- d2[ord[q, ]]
    done <- q >= n - 1L | far > dk * (1 + 1e-9)
    # as.vector(): a two-column matrix of positions would be read as pairs
    # of (row, column).
    if (any(done))
      out[todo[done], ] <- matrix(cand[as.vector(pick[, done, drop = FALSE])],
                                  sum(done), k, byrow = TRUE)
    todo <- todo[!done]
    q <- min(n - 1L, 2L * q)
  }
  out
}


#' The k nearest other locations of every location
#'
#' With \pkg{FNN} the kd-tree shortlists \code{k + 1} neighbours and R
#' recomputes their distances.  A row whose next neighbour is no further than
#' its \code{k}-th, to within one part in a thousand million, is settled by
#' \code{.knn_tied_rows()} with the rule of the exact search, so the set
#' returned for every row is the one the exact search gives and the two paths
#' agree to the last bit.  Without \pkg{FNN} every row is searched exactly,
#' which compares each location with every other.
#'
#' @param u Two-column matrix of distinct locations.
#' @param k Number of neighbours, below \code{nrow(u)}.
#' @param use_fnn Whether to use \pkg{FNN}.
#' @return A list: \code{idx} (an \code{nrow(u) x k} matrix of row numbers;
#'   the set matters, the order within a row does not) and \code{d} (the
#'   distance to the furthest of them, the \code{k}-th nearest).
#' @keywords internal
#' @noRd
.knn_locations <- function(u, k, use_fnn = requireNamespace("FNN", quietly = TRUE)) {
  n <- nrow(u); x <- u[, 1L]; y <- u[, 2L]
  d2_of <- function(idx) matrix((x[idx] - x)^2 + (y[idx] - y)^2, n, ncol(idx))
  if (isTRUE(use_fnn)) {
    kq  <- min(k + 1L, n - 1L)
    idx <- matrix(FNN::get.knn(u, k = kq)$nn.index, n, kq)
    d2  <- d2_of(idx)
    dk  <- do.call(pmax, lapply(seq_len(k), function(j) d2[, j]))
    close <- if (kq > k) which(d2[, kq] <= dk * (1 + 1e-9)) else integer(0)
    idx <- idx[, seq_len(k), drop = FALSE]
    if (length(close)) {
      idx[close, ] <- .knn_tied_rows(u, close, k)
      d2 <- d2_of(idx)
      dk <- do.call(pmax, lapply(seq_len(k), function(j) d2[, j]))
    }
    return(list(idx = idx, d = sqrt(dk)))
  }
  idx <- .knn_exact_rows(x, y, seq_len(n), k)
  d2  <- d2_of(idx)
  list(idx = idx, d = sqrt(do.call(pmax, lapply(seq_len(k), function(j) d2[, j]))))
}


#' For chosen points, the nearest centre by exact search
#'
#' One pass over the centres, each compared with every point at once.  A
#' centre replaces the best so far only when it is strictly nearer, so a tie
#' goes to the lower centre number.
#'
#' @param x,y Coordinates of the points.
#' @param cx,cy Coordinates of the centres.
#' @return Integer vector: the centre of each point.
#' @keywords internal
#' @noRd
.nearest_exact <- function(x, y, cx, cy) {
  best <- (x - cx[1L])^2 + (y - cy[1L])^2
  idx  <- rep.int(1L, length(x))
  for (j in seq_along(cx)[-1L]) {
    d2 <- (x - cx[j])^2 + (y - cy[j])^2
    nearer <- d2 < best
    best[nearer] <- d2[nearer]
    idx[nearer]  <- j
  }
  idx
}


#' The nearest centre of every point
#'
#' With \pkg{FNN} the kd-tree shortlists the two nearest centres and R
#' recomputes both distances.  A point whose two candidates are within one
#' part in a thousand million of each other is searched exactly, so the
#' answer is the one \code{.nearest_exact()} gives for every point, with or
#' without \pkg{FNN}: on a lattice, where a point can be the same distance
#' from two centres, the cells do not depend on which is installed.  Up to 16
#' centres the exact search is the faster of the two and is used either way
#' (measured at 1,500 to 100,000 points; the kd-tree wins from about 32).
#'
#' @param x,y Coordinates of the points.
#' @param cx,cy Coordinates of the centres.
#' @param use_fnn Whether to use \pkg{FNN}.
#' @return Integer vector: the centre of each point.
#' @keywords internal
#' @noRd
.nearest_centre <- function(x, y, cx, cy,
                            use_fnn = requireNamespace("FNN", quietly = TRUE)) {
  K <- length(cx)
  if (K == 1L) return(rep.int(1L, length(x)))
  if (!isTRUE(use_fnn) || K <= 16L) return(.nearest_exact(x, y, cx, cy))
  nn <- FNN::get.knnx(cbind(cx, cy), cbind(x, y), k = 2L)$nn.index
  a <- nn[, 1L]; b <- nn[, 2L]
  da <- (x - cx[a])^2 + (y - cy[a])^2
  db <- (x - cx[b])^2 + (y - cy[b])^2
  out <- a
  swap <- db < da
  out[swap] <- b[swap]
  close <- which(pmax(da, db) <= pmin(da, db) * (1 + 1e-9))
  if (length(close)) out[close] <- .nearest_exact(x[close], y[close], cx, cy)
  out
}


#' Density weights of the distinct locations
#'
#' The local density at a location is the number of other points at it and
#' at its \code{k_nn} nearest other locations, over the area of the disc that
#' reaches the furthest of them: \eqn{count / (\pi d^2)}.  Where no location
#' repeats this is \eqn{k_{nn} / (\pi d^2)} with \eqn{d} the distance to the
#' \code{k_nn}-th nearest neighbour, the usual nearest-neighbour estimate.
#' Measuring \eqn{d} between distinct locations keeps it above zero where
#' locations repeat: a point with \code{k_nn} others on top of it has a
#' \code{k_nn}-th neighbour at distance zero and no density by the usual
#' estimate, and a floor on that distance gives every such point the same
#' enormous weight, which hands the whole partition to the repeated
#' locations.  Counted this way a place visited thirty times is thirty points
#' in a disc of ordinary size.
#'
#' The weight of a point is its density to the power \code{power}, scaled so
#' that the weights of the points average one.
#'
#' @param u Distinct locations, already divided by the extent.
#' @param m Points at each location.
#' @param power The power, from 0 to 2.
#' @param k_nn Number of neighbours.
#' @param use_fnn Whether to use \pkg{FNN}.
#' @return The weight of one point at each location.  All one when
#'   \code{power} is 0 or there is one location.
#' @keywords internal
#' @noRd
.density_weights <- function(u, m, power, k_nn, use_fnn) {
  n_u <- nrow(u)
  if (power == 0 || n_u < 2L) return(rep(1, n_u))
  kk  <- min(as.integer(k_nn), n_u - 1L)
  nn  <- .knn_locations(u, kk, use_fnn)
  d   <- pmax(nn$d, .DENSITY_FLOOR_REL)
  cnt <- (m - 1) + rowSums(matrix(m[nn$idx], n_u, kk))
  w   <- (cnt / (pi * d^2))^power
  w / mean(rep.int(w, m))
}


#' Weighted k-means++ starting centres
#'
#' The first centre is a location drawn with probability proportional to its
#' weight; each later one is drawn with probability proportional to its
#' weight times its squared distance from the nearest centre so far (Arthur
#' and Vassilvitskii 2007, with weights).  A location that already holds a
#' centre has distance zero and is not drawn again.
#'
#' @param x,y Coordinates of the locations.
#' @param W Weight of each location (the weights of its points, summed).
#' @param k Number of centres, at most the number of locations.
#' @return Row numbers of the \code{k} starting centres.
#' @keywords internal
#' @noRd
.wkmeanspp <- function(x, y, W, k) {
  n <- length(x)
  idx <- integer(k)
  idx[1L] <- sample.int(n, 1L, prob = W)
  d2 <- (x - x[idx[1L]])^2 + (y - y[idx[1L]])^2
  if (k > 1L) for (j in 2:k) {
    wd  <- W * d2
    cs  <- cumsum(wd)
    tot <- cs[n]
    # The first location whose cumulative weight passes the draw: what
    # findInterval(draw, cs) + 1 returns, without its checks that `cs` is
    # sorted and has no NA, which read the whole vector twice more.
    pick <- if (is.finite(tot) && tot > 0)
      min(n, sum(cs <= stats::runif(1L) * tot) + 1L) else NA_integer_
    # A draw that rounds up to the total lands on the last location, which
    # may already hold a centre; the last location that can be drawn is the
    # one meant.
    if (!is.na(pick) && !(wd[pick] > 0)) pick <- max(which(wd > 0))
    if (is.na(pick)) {
      free <- setdiff(seq_len(n), idx[seq_len(j - 1L)])
      pick <- free[sample.int(length(free), 1L)]
    }
    idx[j] <- pick
    d2 <- pmin.int(d2, (x - x[pick])^2 + (y - y[pick])^2)
  }
  idx
}


#' Density-weighted k-means
#'
#' Each point gets the weight (local density)^\code{power}, and the centres
#' minimise the weighted sum of squared distances from each point to its
#' nearest centre.  A start is \code{k} locations (weighted k-means++, or
#' drawn at random), refined by Lloyd's iteration: each location goes to its
#' nearest centre, each centre moves to the weighted mean of its points, and
#' a centre left without a point stays where it is.  The iteration stops when
#' no coordinate of a centre changes by \code{.DENSITY_TOL_REL} of the extent
#' or more, or after \code{iter.max} steps.  The start with the lowest
#' weighted sum is kept.
#'
#' What makes the result depend on nothing but the points:
#' \itemize{
#'   \item \emph{Row order.}  The work is done on the distinct locations in
#'     coordinate order, so the random draws index the same locations
#'     whatever order the rows come in.
#'   \item \emph{Units.}  The coordinates are centred on their bounding box
#'     and divided by its diagonal first.  The floor on the neighbour
#'     distance and the stopping tolerance are fractions of that diagonal,
#'     and the centres are put back in the caller's units at the end.
#'   \item \emph{Repeated locations.}  See \code{.density_weights()}.  Points
#'     at one location always share a cell.
#'   \item \emph{\pkg{FNN}.}  It shortlists neighbours; every distance that
#'     decides anything is computed in R, and near-ties are settled by an
#'     exact search with one tie rule (see \code{.nearest_centre()}).
#' }
#'
#' With \code{power = 0} every weight is one and this is plain k-means by
#' Lloyd's iteration: the objective of the kind \code{"kmeans"}, reached by a
#' different algorithm, so the centres are close to that kind's and not equal
#' to them.
#'
#' @param xy Numeric matrix of finite coordinates; the first two columns are
#'   used.
#' @param k Number of centres, from 1 to the number of distinct locations.
#' @param power,k_nn The weight's power and the number of neighbours.
#' @param init \code{"kmeans++"} or \code{"random"}.
#' @param nstart,iter.max Starts, and Lloyd steps per start.
#' @param use_fnn Whether to use \pkg{FNN}.
#' @return As \code{.partition_centres()}, with \code{km} \code{NULL} and
#'   two further elements: \code{weights} (the weight of each row of
#'   \code{xy}) and \code{converged} (whether the best start stopped before
#'   \code{iter.max}).  \code{NULL} when \code{k} is below 1 or above the
#'   number of distinct locations.
#' @keywords internal
#' @noRd
.density_kmeans <- function(xy, k, power = 1, k_nn = 8L, init = "kmeans++",
                            nstart = 3L, iter.max = 50L,
                            use_fnn = requireNamespace("FNN", quietly = TRUE)) {
  xy <- unname(as.matrix(xy)[, 1:2, drop = FALSE])
  storage.mode(xy) <- "double"
  if (!all(is.finite(xy)))
    stop(".density_kmeans(): non-finite coordinates.", call. = FALSE)
  n <- nrow(xy); k <- as.integer(k)
  loc <- .partition_locations(xy)
  n_u <- nrow(loc$u)
  if (n < 1L || is.na(k) || k < 1L || k > n_u) return(NULL)

  # Centred on the bounding box and divided by its diagonal.  The midpoint
  # and the diagonal are read off minima and maxima, which do not depend on
  # row order.
  lo <- c(min(loc$u[, 1L]), min(loc$u[, 2L])); hi <- c(max(loc$u[, 1L]), max(loc$u[, 2L]))
  mid <- (lo + hi) / 2
  extent <- sqrt(sum((hi - lo)^2))
  if (!is.finite(extent)) return(NULL)
  if (!(extent > 0)) extent <- 1                  # one location
  zx <- (loc$u[, 1L] - mid[1L]) / extent
  zy <- (loc$u[, 2L] - mid[2L]) / extent
  m  <- loc$m

  # Without FNN every search is exact and compares each point with every
  # other, or with every centre.  Said once, where it starts to be felt.
  if (!isTRUE(use_fnn) && n_u > 5000L)
    .log_warn_once("density_no_fnn",
                   paste0("partition = \"density\" on %d locations without the FNN ",
                          "package: the neighbour searches compare every pair of ",
                          "locations, and every location with every centre at each ",
                          "step. The result is the same; install FNN to get it sooner."),
                   n_u)
  w  <- .density_weights(cbind(zx, zy), m, power, k_nn, use_fnn)
  W  <- m * w
  zw <- cbind(W, zx * W, zy * W)

  nstart <- max(1L, as.integer(nstart)); iter.max <- max(1L, as.integer(iter.max))
  best <- NULL
  objs <- rep(NA_real_, nstart)
  for (r in seq_len(nstart)) {
    start <- if (identical(init, "random")) sample.int(n_u, k) else .wkmeanspp(zx, zy, W, k)
    cx <- zx[start]; cy <- zy[start]
    it <- 0L; moved <- Inf
    while (it < iter.max) {
      it <- it + 1L
      a  <- .nearest_centre(zx, zy, cx, cy, use_fnn)
      s  <- rowsum(zw, a, reorder = TRUE)
      id <- as.integer(rownames(s))
      nx <- cx; ny <- cy
      nx[id] <- s[, 2L] / s[, 1L]
      ny[id] <- s[, 3L] / s[, 1L]
      moved <- max(abs(nx - cx), abs(ny - cy))
      cx <- nx; cy <- ny
      if (moved < .DENSITY_TOL_REL) break
    }
    # The cells are those of the centres as they now stand.
    a   <- .nearest_centre(zx, zy, cx, cy, use_fnn)
    obj <- sum(W * ((zx - cx[a])^2 + (zy - cy[a])^2))
    if (!is.finite(obj) || !all(is.finite(cx)) || !all(is.finite(cy))) next
    objs[r] <- obj
    if (is.null(best) || obj < best$obj)
      best <- list(cx = cx, cy = cy, a = a, obj = obj, iter = it,
                   converged = moved < .DENSITY_TOL_REL)
  }
  if (is.null(best)) return(NULL)

  a <- best$a
  size <- as.integer(tabulate(rep.int(a, m), nbins = k))
  # Unweighted, about each cell's own mean: the within-cell sum of squares of
  # the coordinates, which depends on the cells alone and is the quantity the
  # kind "kmeans" reports.
  sm  <- rowsum(cbind(m, zx * m, zy * m), a, reorder = TRUE)
  id  <- as.integer(rownames(sm))
  mx  <- numeric(k); my <- numeric(k)
  mx[id] <- sm[, 2L] / sm[, 1L]; my[id] <- sm[, 3L] / sm[, 1L]
  ss  <- rowsum(m * ((zx - mx[a])^2 + (zy - my[a])^2), a, reorder = TRUE)
  withinss <- numeric(k)
  withinss[as.integer(rownames(ss))] <- ss[, 1L] * extent^2

  cluster <- integer(n); cluster[loc$order] <- a[loc$gid]
  weights <- numeric(n); weights[loc$order] <- w[loc$gid]
  ok <- objs[is.finite(objs)]
  list(centers = cbind(best$cx * extent + mid[1L], best$cy * extent + mid[2L]),
       cluster = cluster, size = size, withinss = withinss, wss = sum(withinss),
       objective = best$obj * extent^2,
       spread = if (length(ok) > 1L && min(ok) > 0) (max(ok) - min(ok)) / min(ok) else 0,
       n_ok = length(ok), iter = best$iter, km = NULL,
       weights = weights, converged = best$converged)
}
