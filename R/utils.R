# ===========================================================================
# utils.R — Shared utilities, internal helpers, and core metric functions
#
# Migrated from 00_setup.R. All library() calls removed; dependencies
# declared in DESCRIPTION Imports.
# ===========================================================================

#' Null-coalescing operator
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

#' Resolve an object's CRS, or NULL when it has none
#'
#' `ensure_projected(target_crs = )` deliberately errors on an unusable
#' `target_crs`, so that a typo cannot silently turn it into a no-op.  Internal
#' callers, though, routinely derive the target from another object
#' (`sf::st_crs(training_data)`), and that object is legitimately allowed to
#' carry no CRS.  Passing `NA_crs_` straight through would turn a supported
#' CRS-less workflow into a hard error, so those call sites funnel the value
#' through here and get `NULL` ("no target, pick one automatically") instead.
#'
#' @param x An sf/sfc object or anything `sf::st_crs()` accepts.
#' @return An `sf::crs` object, or NULL when the CRS is missing.
#' @keywords internal
#' @noRd
.crs_or_null <- function(x) {
  cr <- tryCatch(sf::st_crs(x), error = function(e) NA)
  if (length(cr) == 0L || all(is.na(cr))) NULL else cr
}

#' Reproject to a target CRS, or stamp it when the source has none
#'
#' `sf::st_transform()` refuses a CRS-less object ("cannot transform sfc object
#' with missing crs"), so any `if (!is.null(crs)) st_transform(x, crs)` turns a
#' `crs =` argument into a hard error for exactly the users most likely to pass
#' it: those whose data carries no CRS. Reprojection is impossible there, but
#' assumption is not. The lon/lat heuristic runs first, exactly as
#' `ensure_projected()` does: coordinates that look like degrees are taken as
#' EPSG:4326 and REPROJECTED; only when they do not is the target stamped. Both
#' branches warn.  A target in longitude and latitude is not stamped on
#' coordinates that cannot be in it: that is refused.
#'
#' @param x An sf or sfc object.
#' @param crs Target CRS (anything `sf::st_crs()` accepts).
#' @param what Label naming `x` in the warning.
#' @param caller Calling function name, for the warning.
#' @param geo_remedy Words appended to the refusal's remedy (after "Set the
#'   CRS of `what` with sf::st_crs()"), or NULL.
#' @param lat_max The latitude, north or south, beyond which a CRS-less `x`
#'   is refused a target in longitude and latitude.  90 by default.  A
#'   caller that tests the rows of the stamped layer for a latitude beyond
#'   90 and names them (create_grid_polygons() and build_tessellation(),
#'   for `boundary`) passes 180, so that its refusal names the row.
#' @return `x` reprojected to, or stamped with, `crs`.  A CRS-less `x` with
#'   an x outside -180 to 360 or a y beyond `lat_max` north or south is
#'   refused when `crs` is in longitude and latitude: its coordinates cannot
#'   be in that CRS.
#' @keywords internal
#' @noRd
.transform_or_stamp <- function(x, crs, what = "input", caller = "spatialkit",
                                geo_remedy = NULL, lat_max = 90) {
  if (is.null(crs)) return(x)
  # An NA target means "leave the coordinates in whatever space they are in":
  # build_tessellation() passes it for CRS-less points that the lon/lat
  # heuristic declined to interpret, so that the grid is built in the same
  # unnamed space as the points rather than being projected on its own.
  if (inherits(crs, "crs") && is.na(crs)) return(x)
  if (is.na(sf::st_crs(x))) {
    # Stamping is the LAST resort, not the first.  ensure_projected() -- which
    # this helper's own documentation says it matches -- runs the lon/lat
    # heuristic first and REPROJECTS when the coordinates look like degrees.
    # Stamping them instead put identical input 5,400 km apart depending on
    # which of the two paths a caller happened to take, and it is how
    # summarize_by_cell(deff = "variogram") came to compare degree separations
    # against a metre range (deff = n for every cell) and how
    # assign_features_to_polygons() came to return zero rows for points that
    # build_tessellation() had just tessellated correctly.
    ll <- .looks_like_lonlat(x)
    if (isTRUE(ll$lonlat)) {
      .warn_and_log(
        paste0("%s(): `%s` has no CRS; its coordinates look like lon/lat ",
               "(xmin=%.2f, xmax=%.2f, ymin=%.2f, ymax=%.2f), so they are ",
               "taken as EPSG:4326 and reprojected to the target CRS. Set the ",
               "CRS explicitly with sf::st_crs() to suppress this."),
        caller, what, ll$bb[["xmin"]], ll$bb[["xmax"]],
        ll$bb[["ymin"]], ll$bb[["ymax"]])
      x0 <- sf::st_set_crs(x, 4326)
      out <- .transform_keep_unlocated(x0, crs)
      .warn_rows_unplaced(x0, out, caller = caller, arg = what, remedy = .UNPLACED_REMEDY)
      return(out)
    }
    # Named by its label: most callers have no `crs` argument, so "the
    # supplied `crs`" sent the reader looking for one they had not passed.
    lbl <- tryCatch(sf::st_crs(crs)$input, error = function(e) NULL)
    lbl <- if (is.character(lbl) && length(lbl) == 1L && !is.na(lbl) && nzchar(lbl))
      sprintf("the target CRS ('%s')", lbl) else "the target CRS"
    # A target in longitude and latitude is not stamped on coordinates that
    # cannot be in it: an x outside longitude -180 to 360, or a y beyond 90
    # north or south.  The stamp drew metres at "latitude" 4e6, and sf added
    # its bare "bounding box has potentially an invalid value range for
    # longlat data" after the warning above.  A caller that tests the rows
    # of the stamped layer for a latitude beyond 90 and refuses them by row
    # (the builders, for `boundary`, with class spatialkit_cells_degenerate)
    # passes lat_max = 180, so that a y between 90 and 180 reaches that
    # refusal.  The box is the one .looks_like_lonlat() read, of the
    # coordinates that have a location.
    bb <- ll$bb
    if (.is_longlat(sf::st_crs(crs)) && !is.null(bb) &&
        all(is.finite(bb)) &&
        (bb[["xmin"]] < -180 || bb[["xmax"]] > 360 ||
         bb[["ymin"]] < -lat_max || bb[["ymax"]] > lat_max))
      stop(sprintf(paste0(
        "%s(): `%s` has no CRS and its coordinates lie outside longitude -180 ",
        "to 360 and latitude -90 to 90 (xmin=%.6g, xmax=%.6g, ymin=%.6g, ",
        "ymax=%.6g), so they cannot be in %s, which is in longitude and ",
        "latitude. Set the CRS of `%s` with sf::st_crs()%s."),
        caller, what, bb[["xmin"]], bb[["xmax"]], bb[["ymin"]], bb[["ymax"]],
        lbl, what, if (is.null(geo_remedy)) "" else paste0(", ", geo_remedy)),
        call. = FALSE)
    .warn_and_log(
      "%s(): `%s` has no CRS and its coordinates do not look like lon/lat, so it cannot be reprojected; stamping %s WITHOUT reprojection. Verify the coordinates are already expressed in that CRS, or set the input CRS with sf::st_crs().",
      caller, what, lbl)
    return(sf::st_set_crs(x, crs))
  }
  # A row with a coordinate that is not finite or lies beyond 1e150 comes
  # back empty (.transform_keep_unlocated()), not as the rest of the feature.
  # A row of a projected `x` with a coordinate beyond .LONLAT_MOVE_MAX comes
  # out empty when `crs` is in degrees (.lonlat_of_placed(), which it moves
  # with): moved from Web Mercator at x = 1e20 it never came back from PROJ.
  # A row with a location that `crs` cannot place in full comes back empty,
  # or without the vertices it cannot place, and a warning of class
  # "spatialkit_rows_unplaced" names it in the caller's words: a ring with
  # a vertex at latitude 95 lost that vertex in silence.
  out <- .transform_keep_unlocated(x, crs)
  .warn_rows_unplaced(x, out, caller = caller, arg = what, remedy = .UNPLACED_REMEDY)
  out
}

#' Check whether an sf/sfc object has a geographic (lon/lat) CRS
#'
#' @param x An sf or sfc object.
#' @return Logical; TRUE if x has a geographic CRS, FALSE otherwise.
#' @keywords internal
#' @noRd
.is_longlat <- function(x) {
  cr <- sf::st_crs(x)
  if (is.na(cr)) return(FALSE)
  isTRUE(tryCatch(sf::st_is_longlat(cr), error = function(e) FALSE))
}

#' Areas in the squared unit of the layer's CRS
#'
#' \code{sf::st_area()} gives a projected layer's area in its CRS's unit
#' squared when sf can name that unit (\code{st_crs()$ud_unit}), and in
#' square metres when it cannot and the CRS carries a \code{+to_meter}
#' (a unit of 1000 km, or of 3.28 m): there it multiplies GEOS's area by
#' \code{to_meter^2}.  The package reads areas beside lengths in the CRS's
#' unit, so such an area is taken as GEOS gives it, in that unit squared.
#' Every other layer gets \code{sf::st_area()}'s numbers as they are: lon/lat
#' layers their geodesic square metres, and layers in metres, in a named
#' unit or with no CRS their planar areas.
#'
#' @param x An sf or sfc object.
#' @return A plain numeric vector, one area per geometry.
#' @keywords internal
#' @noRd
.planar_area <- function(x) {
  g  <- sf::st_geometry(x)
  cr <- sf::st_crs(g)
  if (!is.na(cr) && !.is_longlat(g) && is.null(cr$ud_unit)) {
    to_m <- cr$to_meter
    if (!is.null(to_m) && length(to_m) == 1L && !is.na(to_m) && to_m != 1)
      return(as.numeric(sf::st_area(sf::st_set_crs(g, NA))))
  }
  as.numeric(sf::st_area(g))
}

#' A point layer read in its typed class
#'
#' A geometry column of the general class \code{sfc_GEOMETRY} can hold
#' points only: \code{sf::st_cast(x, "GEOMETRY")} builds one, and so do
#' \code{c()} of such a column with points and a point set with
#' \code{[[<-} in place of the last other geometry of a mixed column.
#' \code{sf::st_coordinates()} refuses that class ("not implemented for
#' objects of class sfc_GEOMETRY"), and \code{sf::st_dimension()} and
#' \code{sf::st_transform()} stop with "Not compatible with STRSXP" on the
#' column \code{sf::st_cast()} returns.  When every geometry is a POINT, empty ones included, the
#' column is cast to \code{sfc_POINT}, which holds the same geometries; the
#' layer is then read as the same points in that class are.  When every
#' geometry is a MULTIPOINT, the column is cast to \code{sfc_MULTIPOINT} in
#' the same way.  Any other layer, one that mixes points with other types or
#' POINT with MULTIPOINT among them, is returned as it came.
#'
#' @param x An sf or sfc object.
#' @return \code{x}, with its geometry column in \code{sfc_POINT} when it
#'   was an \code{sfc_GEOMETRY} of points, and in \code{sfc_MULTIPOINT}
#'   when it was one of multipoints.
#' @keywords internal
#' @noRd
.typed_points <- function(x) {
  is_sf <- inherits(x, "sf")
  if (!is_sf && !inherits(x, "sfc")) return(x)
  g <- if (is_sf) sf::st_geometry(x) else x
  if (!inherits(g, "sfc_GEOMETRY") || length(g) == 0L) return(x)
  types <- unique(as.character(sf::st_geometry_type(g, by_geometry = TRUE)))
  if (length(types) != 1L || !types %in% c("POINT", "MULTIPOINT")) return(x)
  g <- sf::st_cast(g, types)
  if (!is_sf) return(g)
  sf::st_geometry(x) <- g
  x
}


# The one place the package hands a line to logger.  Every helper formats its
# message with sprintf() first, so it is marked skip_formatter(): no formatter
# on any index -- including one logger copied into this namespace from a
# user's global configuration -- gets to read a `%` or a `{` in it as syntax.
# And a log line is never worth the computation that produced it.  An
# appender that fails (a session temp directory deleted under the file trace,
# a user's own appender that throws) is swallowed here rather than aborting
# the caller, so the R warning .warn_and_log() raises next still arrives.
# `raising` tells the console appender that the line is about to be raised as
# an R warning too; see .sk_console_appender() in zzz.R.
.sk_log_state <- new.env(parent = emptyenv())
.sk_log <- function(level, msg, raising = FALSE) {
  .sk_log_state$raising <- raising
  on.exit(.sk_log_state$raising <- FALSE, add = TRUE)
  tryCatch(logger::log_level(level, logger::skip_formatter(msg),
                             namespace = "spatialkit"),
           error = function(e) NULL)
  invisible(msg)
}


#' Structured warning via logger
#' @keywords internal
#' @noRd
.log_warn <- function(fmt, ...) {
  # A caution that is also raised as an R warning goes through
  # .warn_and_log(), which marks its line `raising`, so a knitted document
  # shows it once, as the warning, rather than as a message and a warning.
  # A line logged here is a caution only, and reaches the document as a
  # message (see .sk_console_appender()).
  .sk_log(logger::WARN, sprintf(fmt, ...))
}

# A warning about a CHOICE, logged once per session under `key`.  A choice
# made once is warned about once: fit_rf_model(include_coords = TRUE) inside
# a five-fold cv_rf() or a twenty-fit cv_block_size_sweep() used to print the
# same paragraph on every fit, which reads as twenty problems rather than one
# decision.  The first line says it will not repeat.
.warned_once <- new.env(parent = emptyenv())
.log_warn_once <- function(key, fmt, ...) {
  if (isTRUE(.warned_once[[key]])) return(invisible(FALSE))
  .warned_once[[key]] <- TRUE
  .log_warn(paste0(fmt, " (This is logged once per session.)"), ...)
  invisible(TRUE)
}


#' Validate a scalar argument before it reaches \code{as.integer()} or \code{if ()}
#'
#' Seven exported functions took a count or a distance and passed it straight
#' into \code{as.integer()} or an \code{if ()} test.  Unvalidated, a length-2
#' vector aborts with "'length = 2' in coercion to 'logical(1)'", \code{NA}
#' and \code{Inf} abort with "missing value where TRUE/FALSE needed" (because
#' \code{as.integer()} of either is \code{NA}), and so does any value above
#' \code{.Machine$integer.max} -- none of which names the argument the caller
#' passed.  \code{resolution_profile()} already guards that last trap
#' explicitly; this is the same guard, shared.
#'
#' @param x The value to check.
#' @param arg,caller Names used in the message.
#' @param min,max Inclusive bounds.  Pass \code{max = .Machine$integer.max}
#'   wherever the value later reaches \code{as.integer()}.
#' @param what A noun phrase describing what the argument is, for the message.
#' @return \code{x}, invisibly, when it is valid; otherwise an error.
#' @keywords internal
#' @noRd
.check_scalar <- function(x, arg, caller, min = -Inf, max = Inf,
                          what = "a single finite number") {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x))
    stop(sprintf("%s(): `%s` must be %s; got %s.", caller, arg, what,
                 if (length(x) != 1L)
                   sprintf("%s of length %d", class(x)[1L], length(x))
                 else if (is.numeric(x)) format(x)
                 else class(x)[1L]),
         call. = FALSE)
  if (x < min)
    stop(sprintf("%s(): `%s` must be %s and at least %s; got %s.",
                 caller, arg, what, format(min), format(x)), call. = FALSE)
  if (x > max)
    stop(sprintf("%s(): `%s` must be %s and at most %s%s; got %s.",
                 caller, arg, what, format(max),
                 if (identical(max, .Machine$integer.max))
                   " (R's largest integer)" else "",
                 format(x)),
         call. = FALSE)
  invisible(x)
}


#' Log a warning AND raise it as an R condition
#'
#' A logger line is invisible to \code{tryCatch(warning = )},
#' \code{withCallingHandlers()}, \code{testthat::expect_warning()} and
#' \code{options(warn = 2)}.  For a situation the caller should be able to
#' catch or escalate (data dropped, a CRS assumed, a result degraded), the
#' log line is not enough on its own.  Used wherever the reference manual says
#' the function \emph{warns}; purely methodological cautions stay
#' \code{.log_warn()} and their documentation says "logged".
#'
#' @keywords internal
#' @noRd
.warn_and_log <- function(fmt, ...) {
  msg <- sprintf(fmt, ...)
  # Logged first, so the trace keeps the line even when the caller catches
  # the warning with tryCatch(); .sk_log() cannot fail, so the warning
  # always follows.
  .sk_log(logger::WARN, msg, raising = TRUE)
  warning(msg, call. = FALSE)
  invisible(msg)
}


#' Structured info message via logger
#' @keywords internal
#' @noRd
.log_info <- function(fmt, ...) {
  .sk_log(logger::INFO, sprintf(fmt, ...))
}


#' Make geometries valid
#'
#' sf >= 1.0 (a hard dependency of this package) provides st_make_valid()
#' natively, so no fallback path is needed.  (An earlier lwgeom fallback was
#' removed: modern lwgeom no longer exports st_make_valid, and the branch
#' was unreachable anyway with sf >= 1.0 installed.)
#'
#' @param g An sf, sfc, or sfg object.
#' @return The same object with repaired geometries.
#' @keywords internal
#' @noRd
.safe_make_valid <- function(g) {
  suppressWarnings(sf::st_make_valid(g))
}


#' Align the CRS of `a` to match `b` (transform if both defined & differ)
#'
#' A row of `a` with a coordinate that is not finite or lies beyond 1e150
#' comes back empty, and so does a row of a projected `a` with a coordinate
#' beyond `.LONLAT_MOVE_MAX` when `b` is in degrees
#' (`.transform_keep_unlocated()`, `.lonlat_of_placed()`).
#'
#' A row of `a` with a location that the CRS of `b` cannot place in full
#' comes back empty, or without the vertices it cannot place, and a warning
#' of class `spatialkit_rows_unplaced` names it in the words of `caller` and
#' `arg` (`.warn_rows_unplaced()`).
#'
#' @param a sf or sfc object to re-project.
#' @param b sf or sfc reference whose CRS takes precedence.
#' @param caller,arg The function the user called and the argument `a` came
#'   from, for the warning.
#' @param label `a` as the warning names it; the argument in backticks by
#'   default.
#' @return `a`, possibly transformed.
#' @keywords internal
#' @noRd
.align_crs <- function(a, b, caller = "spatialkit", arg = "input",
                       label = sprintf("`%s`", arg)) {
  if (is.null(a) || is.null(b)) return(a)
  if (is.na(sf::st_crs(a)) || is.na(sf::st_crs(b))) return(a)
  if (sf::st_crs(a) == sf::st_crs(b)) return(a)
  # A row with a coordinate that is not finite or lies beyond 1e150 comes
  # back empty (.transform_keep_unlocated()), not as the rest of the feature.
  out <- .transform_keep_unlocated(a, sf::st_crs(b))
  .warn_rows_unplaced(a, out, caller = caller, arg = arg, label = label,
                      remedy = .UNPLACED_REMEDY)
  out
}


#' Remove duplicate features by a coordinate key that does not move with
#' the unit or the origin
#'
#' Works on sf objects and bare sfc vectors with POINT geometry, and also
#' handles MULTIPOINT (or mixed) input by keying each feature on all of its
#' vertices, keeping the deduplication mask aligned per feature.
#'
#' The key is a vertex's place in the frame of \code{.kmeans_frame()}: the
#' coordinates centred on the middle of the bounding box of every vertex with
#' a location, divided by its diagonal and rounded to \code{.KMEANS_STEP}
#' (about a billionth of it).  Ten decimals of the layer's unit merged copies
#' 3e-8 m apart in kilometres and not in metres, and copies 1e-9 m apart with
#' one false origin and not with another; points within about a billionth of
#' the extent of each other now key alike in any unit or origin, but for a
#' point that lies within its own rounding noise of the middle of a step.  A
#' vertex without a location (not finite, or beyond 1e150) is keyed by its
#' coordinates as they are.
#'
#' @param g An sf object or sfc vector with POINT/MULTIPOINT geometry.
#' @return `g` with duplicates removed (same class as input).
#' @keywords internal
#' @noRd
.dedup_points <- function(g) {
  # Points in an sfc_GEOMETRY column are keyed as the same points in
  # sfc_POINT (.typed_points()): sf::st_coordinates() refuses that class.
  geom <- .typed_points(if (inherits(g, "sfc")) g else sf::st_geometry(g))
  n <- length(geom)
  if (n == 0L) {
    if (inherits(g, "sfc")) return(g[integer(0)])
    return(g[FALSE, , drop = FALSE])
  }

  gtypes <- as.character(sf::st_geometry_type(geom, by_geometry = TRUE))
  fast <- all(gtypes == "POINT")
  # The vertices of each feature: one row per feature on the fast path;
  # MULTIPOINT (or mixed) input has one row per *vertex* in
  # st_coordinates(), which would misalign a per-feature mask.
  verts <- if (fast) {
    m <- sf::st_coordinates(geom)[, 1:2, drop = FALSE]
    list(xy = m, of = seq_len(n))
  } else {
    ms <- lapply(seq_len(n), function(i) sf::st_coordinates(geom[i])[, 1:2, drop = FALSE])
    list(xy = do.call(rbind, ms), of = rep(seq_len(n), vapply(ms, nrow, integer(1))))
  }
  xy <- verts$xy
  ok <- is.finite(xy[, 1L]) & is.finite(xy[, 2L]) &
    abs(xy[, 1L]) <= .POOL_XY_MAX & abs(xy[, 2L]) <= .POOL_XY_MAX
  vkey <- character(nrow(xy))
  if (any(ok)) {
    z <- .kmeans_frame(xy[ok, , drop = FALSE])$to(xy[ok, , drop = FALSE])
    # + 0 reads -0 as 0, so a point that rounds to the middle from below
    # keys as one that rounds to it from above.
    vkey[ok] <- sprintf("%a %a", z[, 1L] + 0, z[, 2L] + 0)
  }
  vkey[!ok] <- sprintf("raw %a %a", xy[!ok, 1L], xy[!ok, 2L])
  key <- if (fast) vkey else {
    k <- rep("<empty>", n)
    have <- unique(verts$of)
    k[have] <- vapply(split(vkey, factor(verts$of, levels = have)),
                      paste, character(1), collapse = ";")
    k
  }

  mask <- !duplicated(key)
  if (inherits(g, "sfc")) g[mask] else g[mask, , drop = FALSE]
}


#' Assert that an object is sf with one of the expected geometry types
#'
#' @param x An object.
#' @param what Character vector of acceptable geometry type names.
#' @param label Label used in error messages.
#' @param caller Name of the user-facing function, for the message prefix.
#'   Defaults to the name of the function that called this one.
#' @keywords internal
#' @noRd
.assert_sf <- function(x, what = c("POINT", "POLYGON", "MULTIPOLYGON"),
                       label = deparse(substitute(x)),
                       caller = NULL) {
  if (is.null(caller)) {
    cl <- sys.call(-1L)
    caller <- if (is.null(cl)) "spatialkit" else
      sub("^.*:::?", "", deparse(cl[[1L]])[1L])   # drop a pkg:: prefix
  }
  if (!inherits(x, "sf"))
    stop(sprintf("%s(): `%s` must be an sf object%s.", caller, label,
                 if (is.list(x) && !is.null(x$cells))
                   " (this looks like a build_tessellation() result; pass its `$cells`)"
                 else ""),
         call. = FALSE)
  gcls <- unique(as.character(sf::st_geometry_type(x, by_geometry = TRUE)))
  # `all`, not `any`: a mixed-geometry layer with one acceptable type used to
  # pass, contradicting the error text below and letting e.g. a POINT/POLYGON
  # mix through a POINT-only check.  The explicit length check keeps a zero-row
  # layer failing as it did under `any()` -- all() is vacuously TRUE on the
  # empty set.
  if (length(gcls) == 0L || !all(gcls %in% what))
    stop(sprintf("%s(): `%s` geometry must be one of: %s (found: %s).",
                 caller, label, paste(what, collapse = ", "),
                 paste(gcls, collapse = ", ")),
         call. = FALSE)
}


#' Temporarily set the RNG seed and restore it on exit
#'
#' Returns an on.exit-compatible cleanup expression. Call inside a function:
#'   cleanup <- .with_seed(seed); on.exit(cleanup(), add = TRUE)
#'
#' @param seed Integer seed, or NULL to skip.
#' @return A zero-argument function that restores the previous RNG state.
#' @keywords internal
#' @noRd
.with_seed <- function(seed) {
  if (is.null(seed)) return(function() invisible(NULL))
  old_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (old_exists) get(".Random.seed", envir = .GlobalEnv) else NULL
  set.seed(seed)
  function() {
    if (old_exists)
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv))
      rm(".Random.seed", envir = .GlobalEnv)
  }
}


#' The power of two at or below a positive value
#'
#' The first divisor of \code{.spread_scaled()}: \code{log2()} of the
#' largest doubles rounds up to 1024, whose power is \code{Inf}, and of a
#' value just below a power of two can round up to it.
#'
#' @param a A positive finite value.
#' @return A power of two in \code{(a / 2, a]}.
#' @keywords internal
#' @noRd
.pow2_floor <- function(a) {
  s <- 2^min(floor(log2(a)), 1023)
  if (s > a) s / 2 else s
}


#' The mean of a column at any size
#'
#' \code{mean()} sums in a double where a long double is not wider, and a
#' column near the largest double then has a mean of \code{Inf}.  The
#' column is divided by the power of two at or below its largest absolute
#' value, which is exact, and the mean multiplied back.
#'
#' @param y Finite values, at least one.
#' @return The mean of \code{y}.
#' @keywords internal
#' @noRd
.mean_pow2 <- function(y) {
  a <- max(abs(y))
  if (!(a > 0)) return(0)
  s <- .pow2_floor(a)
  mean(y / s) * s
}


#' The training mean an out-of-sample R-squared is scored against
#'
#' \code{mean()} of the finite training responses, less \code{centre} where
#' the fit took the response about one.  Where that mean is finite it is
#' returned as it stands, so the baseline is the one \code{mean()} gives.
#' Where it is not, the values are near the largest double and
#' \code{mean()}, summing in a double where a long double is not wider, has
#' overflowed: the mean is then \code{.mean_pow2()} of the values, and
#' where a value less the centre overflows, twice \code{.mean_pow2()} of the
#' halves.
#'
#' @param y The training responses; values that are not finite are dropped.
#' @param centre One finite value, 0 where the fit took no centre.
#' @return The mean, \code{NaN} where no value is finite (as \code{mean()}).
#' @keywords internal
#' @noRd
.train_mean <- function(y, centre = 0) {
  y <- y[is.finite(y)]
  d <- if (centre != 0) y - centre else y
  m <- mean(d)
  if (is.finite(m) || !length(y)) return(m)
  if (all(is.finite(d))) return(.mean_pow2(d))
  .mean_pow2(y / 2 - centre / 2) * 2
}


#' A sum of squared differences at any size
#'
#' The differences \code{u - v} are divided by the power of two at or below
#' their largest absolute value, so that the largest lies in [1, 2) and no
#' square underflows or overflows where the sum is not itself beyond a
#' double.  Division by a power of two is exact.  Where a difference
#' overflows, the differences are taken on the halves, \code{f = 2}.
#'
#' @param u,v The values; \code{v} one value or one per value of \code{u}.
#' @return A list with \code{m}, \code{s} and \code{f}, the sum being
#'   \code{m * (s * f)^2}, \code{m} between 1 and \code{4 * length(u)} (0
#'   where every difference is 0); \code{NULL} where a difference is not
#'   finite on the halves either (\code{v} not finite).
#' @keywords internal
#' @noRd
.ss_pow2 <- function(u, v) {
  d <- u - v
  f <- 1
  if (!all(is.finite(d))) {
    d <- u / 2 - v / 2
    f <- 2
  }
  a <- max(abs(d))
  if (!is.finite(a)) return(NULL)
  if (!(a > 0)) return(list(m = 0, s = 1, f = 1))
  s <- .pow2_floor(a)
  list(m = sum((d / s)^2), s = s, f = f)
}


#' Compute basic regression metrics (RMSE, MAE, MAPE, SMAPE, R-squared, Adjusted R-squared)
#'
#' Shared across model-prep, cross-validation, and evaluation modules.
#'
#' @param y Numeric vector of observed values.
#' @param yhat Numeric vector of predicted values.
#' @param p Integer number of predictors (for Adjusted R-squared). Default NULL
#'   (Adjusted R-squared omitted).
#' @param y_train_mean Baseline mean for R-squared computation. When NULL (default),
#'   the mean of `y` is used. A scalar value (e.g., per-fold training mean) is
#'   used directly as the baseline. A per-observation vector (length matching `y`)
#'   is filtered in parallel with `y` and `yhat` to remove non-finite cases.
#'   Any other length is an error. Recycling it would silently produce a wrong
#'   R-squared.
#' @param centre The value `y`, `yhat` and `y_train_mean` are given less of:
#'   0 (the default) for values in the response's units, or the centre a fit
#'   regressed its response about (one value, or one per observation for
#'   pooled folds, each the centre of its own fold's fit). RMSE, MAE and R2
#'   do not change under a shift, so they read the centred values, which
#'   keep the spread a double rounds away at the response's size; MAPE and
#'   SMAPE, which read the size of the values, read the centred error over
#'   the values with the centre added back. With 0 the arithmetic is that of
#'   values in the response's units, to the last bit.
#' @return A data.frame with n, RMSE, MAE, MAPE, SMAPE, R2, Adj_R2, n_MAPE and
#'   n_SMAPE. `Adj_R2` is always present, and is `NA` when `p` is NULL or
#'   `n <= p + 1`.
#'
#'   `n` counts the finite (y, yhat) pairs. It is **not** the number of rows
#'   `MAPE` and `SMAPE` were averaged over: both have a denominator that can be
#'   zero, and each drops the rows where its own denominator vanishes (`MAPE`
#'   where `y == 0`, `SMAPE` where `|y| + |yhat| == 0`), returning `NA` only
#'   when no row qualifies. "Zero" is relative to the data, as for every
#'   metric here: a denominator no larger than `100 * .Machine$double.eps`
#'   times the largest one. `R2` is `NA` where the total sum of squares has
#'   an RMS deviation no larger than that fraction of the RMS of `y` and `y`,
#'   read again centred and scaled (`.constant_column()`), shows no spread
#'   beyond its rounding. The result therefore does not depend on the
#'   response's units or origin.
#'   `n_MAPE` and `n_SMAPE` are the row counts each was
#'   actually averaged over, so that a percentage error over a subset is
#'   labelled as one; they equal `n` whenever no row was dropped, and are `0`
#'   in the empty frame. They sit last so that code addressing the first seven
#'   columns by position is unaffected. The user-facing help
#'   (see `model_metrics()`'s "Percentage errors on responses with zeros")
#'   still tells callers to prefer RMSE/MAE/R2 on a response that can be zero:
#'   the count makes the subsetting visible, not harmless.
#' @keywords internal
#' @noRd
.compute_reg_metrics <- function(y, yhat, p = NULL, y_train_mean = NULL,
                                 centre = 0) {
  ok <- is.finite(y) & is.finite(yhat)

  # A scalar baseline (per-fold training mean) is used as-is; a per-observation
  # vector (pooled) is filtered in parallel with `y`/`yhat`.  Anything else
  # would recycle against the filtered `y` and silently distort R-squared.
  if (!is.null(y_train_mean) && length(y_train_mean) != 1L) {
    if (length(y_train_mean) != length(ok))
      stop(sprintf(
        ".compute_reg_metrics(): `y_train_mean` must be length 1 (a scalar baseline) or length %d (one per observation); got %d.",
        length(ok), length(y_train_mean)
      ), call. = FALSE)
    y_train_mean <- y_train_mean[ok]
  }
  if (length(centre) != 1L) {
    if (length(centre) != length(ok))
      stop(sprintf(
        ".compute_reg_metrics(): `centre` must be length 1 or length %d (one per observation); got %d.",
        length(ok), length(centre)
      ), call. = FALSE)
    centre <- centre[ok]
  }
  if (!all(is.finite(centre)))
    stop(".compute_reg_metrics(): `centre` must be finite.", call. = FALSE)
  centred <- any(centre != 0)

  y <- y[ok]; yhat <- yhat[ok]
  n <- length(y)
  if (n == 0L) return(data.frame(n = 0L, RMSE = NA_real_, MAE = NA_real_,
                                 MAPE = NA_real_, SMAPE = NA_real_,
                                 R2 = NA_real_, Adj_R2 = NA_real_,
                                 n_MAPE = 0L, n_SMAPE = 0L))
  rss  <- sum((y - yhat)^2)

  baseline <- if (!is.null(y_train_mean)) y_train_mean else mean(y)
  tss  <- sum((y - baseline)^2)
  # Squares of errors below about 1e-162 fall below the smallest double, and
  # squares of large errors sum past the largest one: on the values as they
  # stand RSS is 0 (RMSE 0, R2 1 whatever the fit) or Inf (RMSE Inf, R2 NA
  # or -Inf), and TSS the same.  Where a sum is not finite, or is below
  # n * .POOL_Y_MIN^2 while its terms are not all 0, that sum alone is read
  # again on its own terms divided by a power of two (.ss_pow2()), which is
  # exact: the sum is m * (s * f)^2.  The other sum keeps its value, as m
  # with s = f = 1, so a far prediction or baseline does not move the scale
  # of the sum it is not part of.  RMSE is multiplied back, and R2 takes
  # the ratio of the sums as m_rss / m_tss times the square of the ratio of
  # their scales.  With both sums as they stand, every product below is by
  # 1 and the arithmetic is the old one to the bit.
  lim <- n * .POOL_Y_MIN^2
  sr <- if (!is.finite(rss) || (rss < lim && any(y != yhat))) .ss_pow2(y, yhat)
  if (is.null(sr)) sr <- list(m = rss, s = 1, f = 1)
  st <- if (!is.finite(tss) || (tss < lim && any(y != baseline)))
    .ss_pow2(y, if (!is.null(y_train_mean)) y_train_mean else .mean_pow2(y))
  if (is.null(st)) st <- list(m = tss, s = 1, f = 1)
  rmse <- sqrt(sr$m / n) * sr$s * sr$f
  mae  <- mean(abs(y - yhat))

  # "Zero" means zero at the scale of the data, 100 machine epsilons of its
  # magnitude, for every metric.  The thresholds were absolute (1e-14 for a
  # denominator, var(y) > 2.2e-16 for R2), so a response in small units --
  # sd below about 1.5e-8 -- lost its R2 while RMSE and MAE were fine, and
  # select_features_forward(metric = "R2") then selected nothing.  A constant
  # response still has a TSS of exactly 0 (or of rounding noise, below the
  # threshold), so its R2 stays NA.
  tol <- 100 * .Machine$double.eps
  # The percentage errors read the size of the values, so centred values
  # are read with their centre added back; the error itself is that of the
  # centred values.
  ya  <- if (centred) y + centre else y
  yha <- if (centred) yhat + centre else yhat
  nz <- abs(ya) > tol * max(abs(ya))
  mape <- if (any(nz)) mean(abs((y[nz] - yhat[nz]) / ya[nz])) * 100 else NA_real_

  denom <- abs(ya) + abs(yha)
  smape_ok <- denom > tol * max(denom)
  smape <- if (any(smape_ok)) {
    mean(2 * abs(y[smape_ok] - yhat[smape_ok]) / denom[smape_ok]) * 100
  } else NA_real_

  # TSS against the squared magnitude of y: R2 needs the spread about the
  # baseline to exceed rounding, i.e. an RMS deviation above tol times the
  # RMS of y.  (Relative to sum(y^2) itself, not squared tol, it would turn
  # R2 NA for an ordinary response on a large offset, 1e8 +- 0.1.)  At a
  # larger offset the bound reads the size of the values, not their
  # spread: 1e15 plus a spread of 2.3 has an RMS deviation under 100 eps
  # of its RMS.  Where the bound finds no spread, y is read again centred
  # and scaled (.constant_column()), and R2 is NA only where that reading
  # finds none either.  TSS and RSS are sums of differences, which a double
  # holds at the offset, so only the test changes.  The test reads y's own
  # spread whatever the baseline.
  # On a TSS read again the test reads y in the units of its scale.  An
  # exact fit has RSS 0 at any scale (0 times an Inf ratio would be NaN).
  # .constant_column() takes y as it stands: its own second reading is on
  # the column scaled, and where its raw sums underflow or overflow they do
  # so together and leave the answer to that reading.
  ss_y <- sum(((y / st$f) / st$s)^2)
  g  <- ((sr$s / st$f) / st$s) * sr$f
  rr <- if (sr$m == 0) 0 else sr$m / st$m * g * g
  r2 <- if (st$m > tol^2 * ss_y) 1 - rr
        else if (is.finite(st$m) && st$m > 0 && !.constant_column(y)) 1 - rr
        else NA_real_

  adj_r2 <- NA_real_
  if (!is.null(p) && is.finite(r2) && n > (p + 1L)) {
    adj_r2 <- 1 - (1 - r2) * (n - 1) / (n - p - 1)
  }

  # The rows each percentage error was averaged over.  Appended after the
  # seven original columns on purpose: the metric block's positions are
  # stable, and every frame built from this one adds these two at the end of
  # its own metric block (.cv_overall_metrics(), .cv_fit_one_fold()).
  data.frame(n = n, RMSE = rmse, MAE = mae, MAPE = mape, SMAPE = smape,
             R2 = r2, Adj_R2 = adj_r2,
             n_MAPE = sum(nz), n_SMAPE = sum(smape_ok))
}


#' Sanitize a worker-core count to a safe positive integer
#'
#' \code{parallel::detectCores()} is documented to return \code{NA} on
#' platforms where the core count cannot be determined; arithmetic on that
#' NA (e.g. \code{max(1L, detectCores() - 1L)}) propagates it and later
#' crashes \code{if (cores > 1L)} checks or gets passed to backends.
#' This helper collapses NA/non-finite/invalid input to \code{fallback}.
#'
#' @param cores Candidate core count (any type).
#' @param fallback Integer used when `cores` is unusable. Default 1L.
#' @return A positive integer scalar.
#' @keywords internal
#' @noRd
.sanitize_core_count <- function(cores, fallback = 1L) {
  cores <- suppressWarnings(as.integer(cores[1L]))
  if (length(cores) != 1L || is.na(cores) || cores < 1L) fallback else cores
}


#' Compute a distance matrix or vector with optional subsampling
#'
#' @param xy Numeric matrix of coordinates.
#' @param max_n Maximum number of rows to use. Default 1000.
#' @param seed RNG seed for reproducible subsampling. Default 42.
#' @return Numeric vector of pairwise distances.
#' @keywords internal
#' @noRd
.safe_dist <- function(xy, max_n = 1000L, seed = 42L) {
  n <- nrow(xy)
  if (n <= 1L) return(numeric(0))
  if (n > max_n) {
    cleanup <- .with_seed(seed)
    on.exit(cleanup(), add = TRUE)
    xy <- xy[sample.int(n, max_n), , drop = FALSE]
  }
  as.numeric(stats::dist(xy))
}


#' Refuse arguments in `...` that nothing downstream will read
#'
#' A function whose `...` is forwarded to another, or documented as
#' "Ignored", has no terminal check anywhere: a name that matches no formal
#' vanishes without a condition.  For the evaluation functions the one
#' argument they exist to take is `newdata`, and a one-character slip
#' (`newdta = hold`) silently turned an out-of-sample RMSE of 25.24 into the
#' in-sample 1.086 with the same return shape.  `create_grid_polygons()`, which
#' has no `...`, rejects the same typo.
#'
#' @param dots `list(...)` from the caller.
#' @param fn Name of the user-facing function, for the message.
#' @param allowed Character vector of names that ARE consumed downstream.
#' @return `dots`, invisibly, when every name is allowed.
#' @keywords internal
#' @noRd
.check_dots <- function(dots, fn, allowed = character(0)) {
  if (!length(dots)) return(invisible(dots))
  nm <- names(dots)
  if (is.null(nm)) nm <- rep("", length(dots))
  unnamed <- !nzchar(nm)
  bad     <- !unnamed & !(nm %in% allowed)
  if (any(unnamed))
    stop(fn, "(): ", sum(unnamed), " unnamed argument(s) in `...` would be ",
         "ignored. Name every argument.", call. = FALSE)
  if (any(bad))
    stop(fn, "(): unused argument(s) ",
         paste(sprintf("`%s`", nm[bad]), collapse = ", "),
         if (length(allowed)) paste0(" (this method accepts ",
                                     paste(sprintf("`%s`", allowed), collapse = ", "),
                                     " through `...`)")
         else " (this method takes nothing through `...`)",
         ". A misspelt `newdata` here would silently return in-sample values.",
         call. = FALSE)
  invisible(dots)
}


#' Catch a misspelt `newdata` on the evaluation functions
#'
#' `model_metrics()`, `evaluate_insample()` and `compare_models()` forward
#' `...` to `predict()`, which checks it, but only on the out-of-sample
#' branch.  A misspelt `newdata` leaves `newdata` NULL and the typo in `...`,
#' the in-sample branch never calls `predict()`, and the in-sample metrics
#' come back with the out-of-sample return shape.  Arguments in `...` with no
#' `newdata` are therefore an error here.
#'
#' @param dots `list(...)`.
#' @param newdata The caller's `newdata`.
#' @param fn User-facing function name.
#' @keywords internal
#' @noRd
.check_dots_newdata <- function(dots, newdata, fn) {
  if (is.null(newdata) && length(dots)) {
    nm <- names(dots)
    if (is.null(nm)) nm <- rep("", length(dots))
    nm <- ifelse(nzchar(nm), sprintf("`%s`", nm), "<unnamed>")
    stop(fn, "(): argument(s) ", paste(nm, collapse = ", "),
         " were supplied but `newdata` was not; nothing reads them in the ",
         "in-sample case, and the in-sample metrics would be returned as if ",
         "they were out-of-sample. Did you mean `newdata = `?", call. = FALSE)
  }
  invisible(dots)
}


#' Belsley condition index of a design matrix
#'
#' The collinearity diagnostic the literature actually thresholds: the ratio
#' of the largest to the smallest singular value of the design matrix after
#' each column has been scaled to unit Euclidean length (Belsley, Kuh & Welsch
#' 1980; for GWR, Wheeler & Tiefelsdorf 2005).  Columns are scaled but NOT
#' centred, so an intercept column stays in and a column that is constant
#' inside a window shows up as collinear with it.  The conventional threshold
#' is 30.  \code{kappa()} on the raw matrix depends on the predictors' units.
#' A design with condition index 1322, whose local coefficients ran from -86
#' to +150 around a true value of 2, had a raw kappa under 1e6 and raised
#' nothing.
#'
#' @param X Numeric matrix.
#' @return A non-negative number; \code{Inf} for an exactly singular design or
#'   a zero column.
#' @keywords internal
#' @noRd
.condition_index <- function(X) {
  X <- as.matrix(X)
  if (!nrow(X) || !ncol(X)) return(NA_real_)
  if (nrow(X) < ncol(X)) return(Inf)
  nrm <- sqrt(colSums(X^2))
  if (any(!is.finite(nrm)) || any(nrm == 0)) return(Inf)
  sv <- tryCatch(svd(sweep(X, 2L, nrm, "/"), nu = 0, nv = 0)$d,
                 error = function(e) NULL)
  if (is.null(sv) || !length(sv)) return(Inf)
  if (min(sv) <= .Machine$double.eps * max(sv)) return(Inf)
  max(sv) / min(sv)
}


#' Ellipsoidal (WGS84) geodesic distance between lon/lat pairs
#'
#' Vincenty's inverse formula, vectorised over pairs.  \code{sf::st_distance()}
#' on lon/lat geometry uses s2's SPHERE (R = 6371 km), and the sphere-to-WGS84
#' gap of 0.24-0.56% is the same size as the projection distortions
#' \code{.crs_distance_error()} compares.  The "measured error X% vs Y%"
#' figures were therefore off by up to half a percentage point and the
#' least-distorting candidate was mis-ranked in 16 of 40 random wide extents.
#' The ellipsoidal distance needs no Suggests package (\pkg{lwgeom} is not a
#' dependency).
#'
#' @param lon1,lat1,lon2,lat2 Numeric vectors in decimal degrees, recycled.
#' @return Distances in metres, \code{NA} for a pair with a coordinate that
#'   is not finite.  Nearly antipodal pairs, where the iteration does not
#'   converge, fall back to the spherical distance.
#' @keywords internal
#' @noRd
.geod_distance <- function(lon1, lat1, lon2, lat2) {
  a <- 6378137; f <- 1 / 298.257223563; b <- a * (1 - f)
  rad <- pi / 180
  # A pair with a coordinate that is not finite is left out of the
  # iteration and given NA.  Its test below was NA, and the update of
  # `active` stopped with "NAs are not allowed in subscripted assignments",
  # so one such row made every candidate CRS score NA.  Its latitudes are
  # set to 0 so that tan() gives no warning on them.
  n <- max(length(lon1), length(lat1), length(lon2), length(lat2))
  lon1 <- rep_len(lon1, n); lat1 <- rep_len(lat1, n)
  lon2 <- rep_len(lon2, n); lat2 <- rep_len(lat2, n)
  fin <- is.finite(lon1 + lat1 + lon2 + lat2)
  lat1[!fin] <- 0; lat2[!fin] <- 0
  U1 <- atan((1 - f) * tan(lat1 * rad)); U2 <- atan((1 - f) * tan(lat2 * rad))
  # A pair whose longitudes differ by more than 180 degrees is measured the
  # short way round, so that the stopping test below reads its true
  # longitude difference: a pair a few millimetres across longitude 180 has
  # a difference close to 360 degrees as written.  Each part is taken from
  # its own side of 180, where the subtraction is exact, so the difference
  # keeps the precision of a pair away from 180.
  dlon <- lon2 - lon1
  w <- which(fin & dlon < -180)
  dlon[w] <- (lon2[w] + 180) + (180 - lon1[w])
  e <- which(fin & dlon > 180)
  dlon[e] <- (lon2[e] - 180) - (lon1[e] + 180)
  L  <- dlon * rad
  sU1 <- sin(U1); cU1 <- cos(U1); sU2 <- sin(U2); cU2 <- cos(U2)
  lambda <- L
  sinSigma <- cosSigma <- sigma <- cosSqAlpha <- cos2SigmaM <- numeric(n)
  active <- fin
  for (iter in seq_len(200L)) {
    sl <- sin(lambda[active]); cl <- cos(lambda[active])
    ss <- sqrt((cU2[active] * sl)^2 +
               (cU1[active] * sU2[active] - sU1[active] * cU2[active] * cl)^2)
    cs <- sU1[active] * sU2[active] + cU1[active] * cU2[active] * cl
    sg <- atan2(ss, cs)
    sinAlpha <- ifelse(ss == 0, 0, cU1[active] * cU2[active] * sl / ss)
    csa <- 1 - sinAlpha^2
    c2sm <- ifelse(csa == 0, 0, cs - 2 * sU1[active] * sU2[active] / csa)
    C  <- f / 16 * csa * (4 + f * (4 - 3 * csa))
    lambda_new <- L[active] + (1 - C) * f * sinAlpha *
      (sg + C * ss * (c2sm + C * cs * (-1 + 2 * c2sm^2)))
    sinSigma[active] <- ss; cosSigma[active] <- cs; sigma[active] <- sg
    cosSqAlpha[active] <- csa; cos2SigmaM[active] <- c2sm
    # The change is measured against L, not in radians: each pass shrinks
    # it by about f cos^2 U, so on an east-west pair under about
    # 1e-12 / (f cos^2 U) rad the first pass already moved lambda by less
    # than 1e-12, and sigma, taken from the starting lambda, was short by
    # f cos^2 U (0.34% at the equator, on pairs of a few millimetres).  A
    # pair 1 rad or more apart in longitude keeps the test in radians, so
    # nearly antipodal pairs fall back to the sphere as before.  A pair
    # along a meridian (L = 0) stops on the first pass, as does a pair at
    # longitudes 180 and -180 on one parallel.
    done <- abs(lambda_new - lambda[active]) <= 1e-12 * pmin(1, abs(L[active]))
    lambda[active] <- lambda_new
    active[active] <- !done
    if (!any(active)) break
  }
  uSq <- cosSqAlpha * (a^2 - b^2) / b^2
  A <- 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)))
  B <- uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)))
  dSigma <- B * sinSigma * (cos2SigmaM + B / 4 *
    (cosSigma * (-1 + 2 * cos2SigmaM^2) -
     B / 6 * cos2SigmaM * (-3 + 4 * sinSigma^2) * (-3 + 4 * cos2SigmaM^2)))
  d <- b * A * (sigma - dSigma)
  # Non-converged (near-antipodal) pairs: spherical fallback.
  if (any(active)) {
    la1 <- lat1[active] * rad; la2 <- lat2[active] * rad
    dl  <- L[active]
    d[active] <- 6371008.8 * acos(pmin(1, pmax(-1,
      sin(la1) * sin(la2) + cos(la1) * cos(la2) * cos(dl))))
  }
  d[fin & lon1 == lon2 & lat1 == lat2] <- 0
  d[!fin] <- NA_real_
  d
}


# ---------------------------------------------------------------------------
# Row-indexed record attributes
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# estimate_sac_range()'s verdict "no spatial structure"
# ---------------------------------------------------------------------------
# estimate_sac_range() refuses a range in two ways that both say the modelled
# variable has no spatial structure: its test on the nearest neighbours finds
# none, or the variogram is flat at the lags it resolves.  Either way the
# model it attaches is a pure nugget.  resolution_profile(),
# select_resolution(), summarize_by_cell(deff = "variogram"),
# kriging_adequacy() and plot.sac_range() each act on that verdict, so the
# wording lives here and they ask .sac_no_structure() instead of matching a
# string of their own.
.sac_reason_no_structure <-
  "no spatial structure detected at the nearest-neighbour scale"
.sac_reason_flat <-
  "the empirical variogram is flat: no spatial structure at the lags resolved"

# TRUE when `sac` is a range estimate refused because the variable shows no
# spatial structure.  A reason that was not written by estimate_sac_range()
# (a stored object, a `sac` built by hand) counts when it says "no spatial
# structure", which both of the reasons above do.
.sac_no_structure <- function(sac) {
  if (is.null(sac)) return(FALSE)
  if (is.finite(suppressWarnings(as.numeric(sac))[1L])) return(FALSE)
  why <- attr(sac, "rejected_reason", exact = TRUE)
  is.character(why) && length(why) >= 1L && !is.na(why[1L]) &&
    grepl("no spatial structure", why[1L], fixed = TRUE)
}

# estimate_sac_range()'s floor of 30 counts distinct locations as well as
# rows: rows that share a location are one place to the variogram and to its
# test for spatial structure.  This is the reason of the NA it returns below
# that floor.  It starts with the count that fell short, as the floor on rows
# does, and is not a verdict on structure (it does not say "no spatial
# structure", so .sac_no_structure() is FALSE for it).  `subsampled` is TRUE
# when the layer has 30 or more locations and the n_max subsample kept fewer.
.sac_reason_few_locations <- function(n_locations, n_points, subsampled = FALSE)
  sprintf(paste0("%d distinct location(s) among the %d point(s) %s, fewer than ",
                 "the 30 a variogram range is estimated from"),
          as.integer(n_locations), as.integer(n_points),
          if (isTRUE(subsampled)) "of the n_max subsample"
          else "with a finite value to model")


# Two functions hand back a layer with a record of what happened to its rows:
# prep_model_data() records the rows it dropped, assign_features_to_polygons()
# records the features that matched more than one polygon and had the tie-break
# rule decide for them.  Each record describes the rows the layer had when it
# was built -- and `[.sf` copies attributes through verbatim, so a subset of
# such a layer used to carry, and report, its parent's numbers: after
# prep_model_data() dropped 3 of 40 rows, `attr(pre[1:20, ], "dropped")$n` was
# still 3 with `which` pointing at rows the subset does not contain, and a fit
# built on that subset with `.already_prepped = TRUE` reported
# `info$n_dropped = 3` for 20 rows nothing had been dropped from.
#
# A layer carrying a record is therefore classed, and the `[` method below
# removes the record rather than let a stale one be read as current.  The
# record is stamped with the row count it was computed for as well, because a
# few paths bypass `[`: sf::st_transform() rebuilds the class vector with "sf"
# first, which takes dispatch away from the method.  Everything that reads a
# record goes through .get_row_record(), which refuses one whose stamp no
# longer matches the layer, so a number this package reports is never taken
# from a record that has stopped describing the rows in front of it.

# Registered with S4 so that a classed layer still satisfies dispatch written
# for "sf": methods::as(x, "Spatial") -- which .to_sp() calls on its way into
# GWmodel -- and terra::vect() and friends look the class up in the S4 table,
# and an unregistered class ahead of "sf" fails them with "no method or
# default for coercing". The class now goes after "sf" (see below), but a
# layer built by an earlier version, read back from an .rds, has it first.
# Registering only the chain up to "sf" is deliberate: sf itself registers
# c("sf", "data.frame"), and naming "data.frame" here as well is rejected as
# inconsistent with that.
setOldClass(c("spatialkit_rows", "sf"))

.row_record_attrs <- c("dropped", "ties", "unassigned")

# Attach `value` as the `which` record of `x`, stamped and classed.
#
# The class goes right after "sf", not ahead of it.  vctrs -- behind
# dplyr::bind_rows(), vctrs::vec_rbind() and dplyr::union_all() -- reads the
# first class, and an unknown one ahead of "sf" sent two layers with the same
# record (bind_rows(a, a), or two equal-sized batches) down its same-type
# path into its sf restore method, which failed with 'attr(obj, "sf_column")
# does not point to a geometry column'.  After "sf", vctrs binds an sf, and
# the result carries no record: it describes neither input's rows.  `[`
# still removes the record, reaching `[.spatialkit_rows` through the
# NextMethod() in sf's own `[` method, as it always did for the output of
# sf::st_transform(), which puts "sf" first.
.set_row_record <- function(x, which, value) {
  value$n_rows <- as.integer(nrow(x))
  attr(x, which) <- value
  if (!inherits(x, "spatialkit_rows")) {
    cl <- class(x)
    at <- match("sf", cl)
    class(x) <- if (is.na(at)) c("spatialkit_rows", cl)
                else append(cl, "spatialkit_rows", after = at)
  }
  x
}

# Read the `which` record of `x`, or NULL when there is none and when the one
# there is describes a different set of rows.  An unstamped record is one a
# caller built by hand or one from an object made by an older version; it is
# returned as given rather than second-guessed.
.get_row_record <- function(x, which) {
  rec <- attr(x, which, exact = TRUE)
  if (is.null(rec) || !is.list(rec)) return(NULL)
  if (is.null(rec$n_rows)) return(rec)
  if (!identical(as.integer(rec$n_rows), as.integer(nrow(x)))) return(NULL)
  rec
}


#' Subset a layer that carries a row record
#'
#' \code{\link{prep_model_data}()} and
#' \code{\link{assign_features_to_polygons}()} return a layer with an
#' attribute recording what happened to its rows (\code{"dropped"} for the
#' first; \code{"ties"} and \code{"unassigned"} for the second).  Those
#' records describe the rows the layer was
#' built with, and \code{[} on an \code{sf} object copies attributes through
#' unchanged, which would leave a subset reporting its parent's numbers for a
#' different set of rows.  Subsetting therefore returns a plain layer with the
#' record removed, and so do the \pkg{dplyr} verbs that select or reorder
#' rows (\code{filter()}, \code{slice()}, \code{arrange()},
#' \code{distinct()}); read the record from the layer the function returned,
#' before subsetting it.  Binding such layers (\code{rbind()},
#' \code{dplyr::bind_rows()}) likewise returns a plain \code{sf} layer.
#'
#' Each record carries \code{n_rows}, the number of rows it was made for.
#' \code{sf::st_drop_geometry()} keeps the rows, and with them the record: it
#' returns a data frame of class \code{c("spatialkit_rows", "data.frame")}.
#' Binding such data frames with
#' \code{rbind()} or \code{dplyr::bind_rows()} keeps the first one's record
#' and class, so the record then describes only the first input's rows: its
#' \code{n_rows} no longer equals \code{nrow()} of the result.  The
#' package's own readers ignore a record whose \code{n_rows} does not match;
#' when reading \code{attr(x, "dropped")}, \code{attr(x, "ties")} or
#' \code{attr(x, "unassigned")} yourself from a layer that has been through
#' such steps, check it the same way.
#'
#' @param x A layer returned by \code{\link{prep_model_data}()} or
#'   \code{\link{assign_features_to_polygons}()}.
#' @param ... Passed to the underlying \code{sf} or data frame method.
#' @return The subset, without the row records and without this class.  A
#'   subset that is not a data frame (a single column taken with
#'   \code{drop = TRUE}) is returned unchanged.
#' @examples
#' library(sf)
#' dat <- st_as_sf(
#'   data.frame(x = 1:5, y = 5:1,
#'              resp = c(1, 2, NA, 4, 5), pred = c(1, 2, 3, 4, Inf)),
#'   coords = c("x", "y"), crs = 32632
#' )
#' clean <- prep_model_data(dat, "resp", "pred")   # warns: 2 of 5 rows dropped
#' attr(clean, "dropped")$n          # 2
#' attr(clean[1:2, ], "dropped")     # NULL -- the record does not follow
#' @seealso \code{\link{prep_model_data}()} and
#'   \code{\link{assign_features_to_polygons}()}, the two functions that
#'   attach the records this method removes.
#' @export
`[.spatialkit_rows` <- function(x, ...) {
  y <- NextMethod()
  if (!is.data.frame(y)) return(y)
  for (nm in .row_record_attrs) attr(y, nm) <- NULL
  class(y) <- setdiff(class(y), "spatialkit_rows")
  y
}

# dplyr's row verbs -- filter(), slice(), arrange(), distinct() -- reach the
# data through dplyr_row_slice(), not `[`, and kept the record: filter(a,
# v >= 3) on 5 rows returned 3 rows still reporting ties$n = 3.  They drop it
# as `[` does.  sf's own dplyr methods strip "sf" from the class and call
# NextMethod(), so this runs for an sf layer too.  Registered in .onLoad()
# on dplyr's generic.
.dplyr_row_slice_spatialkit_rows <- function(data, i, ...) {
  y <- NextMethod()
  if (!is.data.frame(y)) return(y)
  for (nm in .row_record_attrs) attr(y, nm) <- NULL
  class(y) <- setdiff(class(y), "spatialkit_rows")
  y
}


# ---------------------------------------------------------------------------
# Diagnoses written to stderr by compiled code
# ---------------------------------------------------------------------------

# ranger's C++ layer writes its real diagnosis straight to stderr and then
# throws a generic condition.  `ranger(mtry = 99)` prints
#   Error: mtry can not be larger than number of variables in data.
#          Ranger will EXIT now.
# while the R error a caller can catch reads "User interrupt or internal
# error." -- which names neither the argument nor the problem.  That line is
# not an R condition, so suppressMessages(), withCallingHandlers() and
# tryCatch() all miss it: it escapes every handler to the console, and into a
# CI log where it reads as an error from a test that is passing, while the
# error the package reports says nothing actionable and fold_status records
# the same nothing.  Diverting the message stream for the duration of the call
# gets both right: the useful text becomes the reported reason, and nothing is
# printed behind the caller's back.
#
# Only one message sink may be active at a time, and testthat, knitr and
# capture.output(type = "message") all use it, so the diversion happens only
# when the stream is undiverted (sink number 2 is stderr itself).  Otherwise
# the call runs exactly as it did before and the caller falls back to R's own
# message.

# Run `fun()`, returning its value or the error it threw, plus whatever it
# wrote to stderr.  Anything written by a call that SUCCEEDED is passed
# straight through to stderr afterwards, so ordinary progress and warning
# output from compiled code is not swallowed.
.call_capturing_stderr <- function(fun, path = tempfile("spatialkit-stderr-")) {
  err <- NULL
  if (!identical(as.integer(sink.number(type = "message")), 2L)) {
    val <- tryCatch(fun(), error = function(e) { err <<- e; NULL })
    return(list(value = val, error = err, stderr = character(0)))
  }
  # A session temp directory deleted under a running session leaves nowhere
  # to divert to, and file() then failed the whole call with "cannot open the
  # connection".  Run it undiverted instead, as when a sink is active.
  con <- tryCatch(suppressWarnings(file(path, open = "wt")),
                  error = function(e) NULL)
  if (is.null(con)) {
    val <- tryCatch(fun(), error = function(e) { err <<- e; NULL })
    return(list(value = val, error = err, stderr = character(0)))
  }
  open <- TRUE
  # Restore the stream whatever happens, an interrupt included: leaving a
  # session with its messages diverted to a deleted temp file would silence
  # every warning from here on.
  on.exit({
    if (!identical(as.integer(sink.number(type = "message")), 2L))
      sink(type = "message")
    if (open) close(con)
    unlink(path)
  }, add = TRUE)
  sink(con, type = "message")
  val <- tryCatch(fun(), error = function(e) { err <<- e; NULL })
  sink(type = "message")
  close(con); open <- FALSE
  txt <- tryCatch(readLines(path, warn = FALSE), error = function(e) character(0))
  unlink(path)
  if (is.null(err) && length(txt))
    cat(txt, file = stderr(), sep = "\n", fill = TRUE)
  list(value = val, error = err, stderr = txt)
}

# The reason to report for a failed call: what the compiled code wrote, when
# R's own message is the placeholder it throws alongside it, and both when the
# R message says something of its own.
.stderr_reason <- function(e, stderr_lines) {
  msg <- conditionMessage(e)
  lines <- trimws(stderr_lines)
  lines <- lines[nzchar(lines)]
  if (!length(lines)) return(msg)
  hit <- grep("^Error", lines)
  diag <- lines[if (length(hit)) hit[1L] else length(lines)]
  diag <- sub("^Error\\s*:?\\s*", "", diag)
  diag <- sub("\\s*Ranger will EXIT now\\.?$", "", diag)
  diag <- trimws(diag)
  if (!nzchar(diag)) return(msg)
  # "User interrupt or internal error." is what ranger throws for every
  # failure its C++ layer diagnoses, so it carries no information of its own.
  if (grepl("User interrupt or internal error", msg, fixed = TRUE)) diag
  else if (identical(msg, diag)) msg
  else paste0(msg, " (", diag, ")")
}


#' Log a warning and raise it with a condition class
#'
#' \code{.warn_and_log()} with a condition class, so that a caller can catch
#' exactly this case with \code{tryCatch(<class> = )} or
#' \code{withCallingHandlers()}, rather than by matching message text.  Used
#' where code running many calls has to tell one outcome from another, and
#' for a warning that changes what a result means.  Every classed warning of
#' the package is raised here: the helpers that fix a class
#' (\code{.warn_deff_fallback()}, \code{.warn_k_reduced()},
#' \code{.aoa_warn_zero_threshold()} and the like) and \code{.gwr_warn()}
#' call it, so each is logged once and raised the same way.  The classes are
#' named \code{spatialkit_<what>}.
#'
#' @param class Character vector of condition classes, most specific first
#'   (\code{"spatialkit_..."}).
#' @param fmt,... Passed to \code{sprintf()}.
#' @param data Named list of fields to carry on the condition (counts a
#'   handler can read without parsing the message).
#' @keywords internal
#' @noRd
.warn_classed <- function(class, fmt, ..., data = list()) {
  msg <- sprintf(fmt, ...)
  # Taken over by the caller for the length of this call (see
  # .answer_for_classed()): neither logged nor raised.
  if (any(class %in% .answered_state$classes)) return(invisible(msg))
  # Already said, in these words, by an earlier call of the same run (see
  # .once_per_run()): neither logged nor raised again.
  if (.said_in_run(class, msg)) return(invisible(msg))
  # raising = TRUE: the warning below shows it, so a knitted document must
  # not repeat the line as a message (see .sk_console_appender()).
  .sk_log(logger::WARN, msg, raising = TRUE)
  warning(do.call(warningCondition,
                  c(list(message = msg, class = class, call = NULL), data)))
  invisible(msg)
}


# ---------------------------------------------------------------------------
# A classed warning the caller answers for itself
# ---------------------------------------------------------------------------
# cell_values_at() places its locations with assign_features_to_polygons()
# and reports a location in no polygon row by row, so that function's
# warning that no feature falls in any polygon says nothing new there.  A
# calling handler can muffle the warning, but not the line .warn_classed()
# logs before it raises it, which reaches the console.  A function that
# answers for a class itself opens a span here around the one call.  While
# the span is open a warning of that class is neither logged nor raised;
# every other warning is logged and raised as usual.

# The classes taken over by the open spans.
.answered_state <- new.env(parent = emptyenv())

#' Open a span in which a classed warning is neither logged nor raised
#'
#' @param classes Condition classes the caller answers for itself.
#' @return A function that closes the span; call it when the call the span
#'   was opened for returns or stops (\code{tryCatch(, finally = )}).  A span
#'   opened inside another adds its classes to the outer one's, and closing
#'   it leaves the outer one's.
#' @keywords internal
#' @noRd
.answer_for_classed <- function(classes) {
  outer <- .answered_state$classes
  .answered_state$classes <- union(outer, classes)
  function() {
    .answered_state$classes <- outer
    invisible(NULL)
  }
}


# ---------------------------------------------------------------------------
# One warning for a run of calls on the same rows
# ---------------------------------------------------------------------------
# compare_models_cv() runs one cv_*() per model and cv_block_size_sweep() one
# cv_spatial() per block size.  Each of those prepares the layer for itself,
# so each found the same incomplete rows and raised the same
# "spatialkit_rows_dropped" warning: two for two models, eight for the default
# ladder.  A function that makes such a run opens it here with the classes it
# wants said once.  While the run is open, a warning of one of those classes
# whose text has been raised before is neither logged nor raised again; one
# with other text (other rows, another layer) is raised as usual.

# The open run: `classes`, and `said`, the class and text of what was raised.
.once_state <- new.env(parent = emptyenv())

#' Open a run in which a classed warning is raised once per text
#'
#' @param classes Condition classes to raise once per distinct message.
#' @return A function that closes the run; register it with \code{on.exit()}.
#'   A run opened inside another adds its classes to the outer one's and
#'   shares its record, so the text is still raised once.
#' @keywords internal
#' @noRd
.once_per_run <- function(classes) {
  outer <- .once_state$run
  .once_state$run <- list(classes = union(outer$classes, classes),
                          said = outer$said %||% character(0))
  function() {
    if (is.null(outer)) {
      .once_state$run <- NULL
    } else {
      outer$said <- .once_state$run$said %||% outer$said
      .once_state$run <- outer
    }
    invisible(NULL)
  }
}

# TRUE when a warning of this class and text was already raised in the open
# run; records it otherwise.  FALSE with no run open or a class it does not
# cover.
.said_in_run <- function(class, msg) {
  run <- .once_state$run
  if (is.null(run) || !any(class %in% run$classes)) return(FALSE)
  key <- paste(class[1L], msg, sep = "\r")
  if (key %in% run$said) return(TRUE)
  .once_state$run$said <- c(run$said, key)
  FALSE
}

# The three classed warnings prep_model_data() raises about a layer, and
# the one make_folds() raises about coordinates far from their origin
# (.warn_xy_imprecise()).  A function that prepares one layer several times,
# or makes folds on it several times, raises each once (.once_per_run()).
.prep_warning_classes <- c("spatialkit_rows_dropped",
                           "spatialkit_response_imprecise",
                           "spatialkit_constant_predictor",
                           "spatialkit_xy_imprecise")


#' Warn where a response lies so far from zero that doubles there do not
#' hold a millionth of its spread
#'
#' At 1e15 doubles are 0.125 apart, so a response of standard deviation 1
#' stored there is rounded to eighths, and so is every mean, prediction and
#' fold score computed at that size: the range, the folds and the selection
#' read off it are those of the rounded values.  The values are what the
#' caller gave and are not refused; this says once per call that the
#' results read that rounding, where one unit in the last place of the
#' largest absolute value (\code{.ulp_of()}) is more than a millionth of
#' the standard deviation (\code{.cell_sd()}, which holds at any size).  A
#' response with no spread beyond a few units in the last place
#' (\code{.constant_column()}), or of fewer than two finite values, is left
#' to the tests for a constant.  In a run of calls (\code{.once_per_run()}) it
#' is said once for the column, by whichever function reads it first.
#'
#' @param y The response.
#' @param caller The function that reads it, named in the message.
#' @param arg The name of the response column, named in the message.
#' @return \code{TRUE} (invisibly) when it warned, else \code{FALSE}.
#' @keywords internal
#' @noRd
.warn_response_imprecise <- function(y, caller, arg) {
  if (!(is.numeric(y) || is.logical(y))) return(invisible(FALSE))
  y <- as.numeric(y)
  y <- y[is.finite(y)]
  if (length(y) < 2L) return(invisible(FALSE))
  u <- .ulp_of(y)
  s <- .cell_sd(y)
  if (!(is.finite(s) && s > 0 && u > 1e-6 * s)) return(invisible(FALSE))
  # A spread of a few units in the last place is the package's reading of a
  # constant (.constant_column()), which the tests for one report.
  if (.constant_column(y)) return(invisible(FALSE))
  # Taken over by the caller (.answer_for_classed()): not said, and not
  # recorded for the run either, so that a later call of the run that the
  # user sees still says it (make_folds()'s range diagnostic muffles every
  # warning of its estimate_sac_range() call).
  if ("spatialkit_response_imprecise" %in% .answered_state$classes)
    return(invisible(FALSE))
  # Keyed on the column, not on the caller's words, so that a run says it
  # once whichever of its functions reads the response first.
  if (.said_in_run("spatialkit_response_imprecise", paste("column", arg)))
    return(invisible(FALSE))
  .warn_classed(
    "spatialkit_response_imprecise",
    paste0("%s(): `%s` lies so far from zero that its values lose their ",
           "precision: doubles there are %s apart, more than a millionth of ",
           "their standard deviation (%s), and the results read that ",
           "rounding. Subtract a constant from it first."),
    caller, arg, format(signif(u, 3)), format(signif(s, 3)),
    data = list(spacing = u, spread = s))
  invisible(TRUE)
}


# The least length whose square is a normal double (about 1.49e-154).  A
# buffer by a shorter distance, or a box with a shorter side, has an area
# that underflows.
.SQRT_XMIN <- sqrt(.Machine$double.xmin)


#' Warn where coordinates lie so far from their origin that doubles there do
#' not hold a millionth of the layer's extent
#'
#' At 1e15 doubles are 0.125 apart, so a layer 4 km across stored there is
#' rounded to eighths of a metre, and so is every edge computed at that size:
#' the edges of blocks (a corner plus a multiple of a block's side), the
#' lines of a grid and the buffer drawn around points.  A point that lies
#' within that rounding of an edge near the origin can then fall on the
#' other side of it.  The coordinates are what the caller gave and are not
#' refused; this says once per call that the result reads that rounding,
#' where one unit in the last place of the largest coordinate
#' (\code{.ulp_of()}) is more than a millionth of the extent, the bound of
#' the package's precision rule.  A layer with no extent (every location at
#' one point) has no spread to resolve and is left to the caller's own
#' tests.  Rows without a location (\code{.located_rows()}) are not read.
#' In a run of calls (\code{.once_per_run()}) the same text is said once.
#'
#' @param x The layer (sf or sfc), or a bounding box of its located
#'   coordinates in a form \code{.edge_tol()} reads.  A list of layers is
#'   read as one, all of them in one CRS.
#' @param caller The function that reads it, named in the message.
#' @param arg The names of the arguments \code{x} holds, named in the
#'   message.
#' @param extent The extent the spacing is compared with; by default the
#'   diagonal of the located coordinates' bounding box.
#' @param of What \code{extent} is, in the message.
#' @param reads What reads the rounding, in the message.
#' @param remedy The remedy, in the message.
#' @return The spacing of doubles at the coordinates and the extent, in a
#'   list, with \code{warned} \code{TRUE} when it warned (invisibly).
#' @keywords internal
#' @noRd
.warn_xy_imprecise <- function(x, caller, arg, extent = NULL,
                               of = if (length(arg) > 1L) "their extent" else "its extent",
                               reads = "so edges and ties are read at that rounding",
                               remedy = "Shift the layer nearer the origin, or project it.") {
  out <- list(warned = FALSE, spacing = NA_real_, extent = NA_real_)
  layers <- if (inherits(x, c("sf", "sfc")) || !is.list(x)) list(x) else x
  boxes <- lapply(layers, function(l) {
    if (!inherits(l, c("sf", "sfc"))) return(l)
    g <- sf::st_geometry(l)
    if (!length(g)) return(NULL)
    # The box of the whole layer leaves out empty rows; rows with a
    # coordinate that is not finite or lies beyond the bound are left out
    # only where the box shows one, so a layer of a million cells is not
    # read row by row.
    bb <- unclass(sf::st_bbox(g))
    if (all(is.finite(bb)) && max(abs(bb)) <= .POOL_XY_MAX) return(bb)
    g <- g[.located_rows(g)]
    if (!length(g)) return(NULL)
    unclass(sf::st_bbox(g))
  })
  boxes <- Filter(Negate(is.null), boxes)
  if (!length(boxes)) return(invisible(out))
  bb <- unlist(boxes)
  v <- suppressWarnings(as.numeric(bb))
  v <- v[is.finite(v)]
  if (!length(v)) return(invisible(out))
  u <- .ulp_of(v)
  d <- if (is.null(extent)) .bb_diagonal(bb) else extent
  out$spacing <- u
  out$extent  <- d
  if (!(is.finite(d) && d > 0 && u > 1e-6 * d)) return(invisible(out))
  crs_of <- vapply(layers, function(l)
    inherits(l, c("sf", "sfc", "bbox")) && !is.na(sf::st_crs(l)), logical(1))
  who <- if (length(arg) > 1L)
    paste0(paste0("`", arg[-length(arg)], "`", collapse = ", "), " and `",
           arg[length(arg)], "` lie")
  else paste0("`", arg, "` lies")
  .warn_classed(
    "spatialkit_xy_imprecise",
    paste0("%s(): %s so far from the origin%s that %s coordinates lose ",
           "their precision: they are rounded to %s, more than a millionth of ",
           "%s (%s), %s. %s"),
    caller, who,
    if (!any(crs_of)) "" else if (length(arg) > 1L) " of their CRS" else " of its CRS",
    if (length(arg) > 1L) "their" else "its",
    format(signif(u, 3)), of, format(signif(d, 3)), reads, remedy,
    data = list(spacing = u, extent = d))
  out$warned <- TRUE
  invisible(out)
}
