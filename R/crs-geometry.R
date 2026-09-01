# -----------------------------------------------------------------------------
# CRS Selection
# -----------------------------------------------------------------------------

#' Pick a sensible local projected CRS for an sf/sfc object
#'
#' Chooses an appropriate projected coordinate reference system for spatial
#' data, favoring a UTM zone based on the dataset's geographic centroid when
#' possible, and falling back to Web Mercator (EPSG:3857) when longitude/latitude
#' information cannot be reliably determined.
#'
#' Extents reaching more than 5 degrees from the central meridian of the
#' candidate UTM zone get an equal-area projection centred on the data instead:
#' Albers for mid-latitudes, Lambert Azimuthal Equal Area for equatorial and
#' high-latitude extents.  Forcing continental data into one UTM zone produces
#' percent-scale distance errors that propagate silently into variogram ranges,
#' block sizes, GWR bandwidth and GP length-scales.  Note that only longitude
#' offset matters: transverse Mercator error scales with distance from the
#' central meridian, so tall narrow north-south extents keep their UTM zone.
#'
#' @param x An sf or sfc object.
#' @return An sf::crs object.
#' @keywords internal
#' @noRd
.pick_local_projected_crs <- function(x) {
  if (!inherits(x, c("sf", "sfc"))) return(sf::NA_crs_)
  crs <- sf::st_crs(x)
  if (!is.na(crs) && !.is_longlat(x)) return(crs)

  x_ll <- tryCatch({
    if (is.na(crs)) stop("No CRS set.")
    sf::st_transform(x, 4326)
  }, error = function(e) x)

  if (is.na(sf::st_crs(x_ll)) || !.is_longlat(x_ll)) return(sf::st_crs(3857))

  ctr <- sf::st_coordinates(sf::st_centroid(sf::st_union(sf::st_geometry(x_ll))))
  if (!is.numeric(ctr) || length(ctr) < 2) return(sf::st_crs(3857))
  lon <- ctr[1]; lat <- ctr[2]

  # st_centroid() on empty or degenerate geometry can return NA/NaN, which
  # passes the is.numeric() check above.  A non-finite centroid cannot place
  # either a UTM zone or an equal-area projection, so fall back explicitly
  # rather than letting NA propagate into the comparisons below.
  if (!is.finite(lon) || !is.finite(lat)) {
    .log_warn(
      ".pick_local_projected_crs(): could not compute a finite centroid (lon = %s, lat = %s); falling back to EPSG:3857.",
      format(lon), format(lat)
    )
    return(sf::st_crs(3857))
  }

  # ---- Reject extents too wide for a single UTM zone ----
  # Transverse Mercator scale error grows with distance from the *central
  # meridian*: k = k0 * (1 + x^2 / (2 R^2)) where x ~ R * dlon * cos(lat).
  # Two consequences drive the logic below:
  #
  #   * Only the LONGITUDE offset matters.  Latitude span does not inflate the
  #     error -- cos(lat) shrinks x, so a tall narrow north-south extent is
  #     precisely what UTM is designed for (a 4-deg-wide, 60-deg-tall strip
  #     peaks at about +0.02%).  Switching such data to an equal-area
  #     projection would make it worse, not better: LAEA distorts distance by
  #     about -1.1% at 30 deg from its centre.
  #   * What matters is the offset from the CENTRAL MERIDIAN OF THE SELECTED
  #     ZONE, not the raw bbox span, because the centroid-derived zone is not
  #     generally centred on the data.
  #
  # Beyond about 5 deg from the central meridian (~0.3% error) an equal-area
  # projection centred on the data is the better choice.  Forcing CONUS into a
  # single zone puts the extent edge 31 deg off the meridian: about +7.5%.
  # That error is silent and it propagates -- estimate_sac_range() returns a
  # range in CRS units, make_folds(block_kfold) sizes blocks in CRS units, and
  # both the GWR bandwidth and the GP length-scale read projected coordinates.
  bb       <- sf::st_bbox(x_ll)
  lon_min  <- as.numeric(bb["xmin"]); lon_max <- as.numeric(bb["xmax"])
  span_lon <- lon_max - lon_min
  span_lat <- as.numeric(bb["ymax"] - bb["ymin"])

  # A bbox wider than a hemisphere cannot be told apart from one that merely
  # straddles the antimeridian, and in the latter case the centroid lands on
  # the opposite side of the planet -- which would centre an equal-area
  # projection 180 deg from the data.  Fall back to the documented
  # can't-determine-reliably behaviour instead of guessing.
  if (is.finite(span_lon) && span_lon > 180) {
    .log_warn(
      paste0(".pick_local_projected_crs(): longitude extent spans %.1f deg. ",
             "This is either global coverage or data straddling the ",
             "antimeridian, which a bounding box cannot distinguish; the ",
             "centroid is unreliable either way. Falling back to EPSG:3857. ",
             "Pass target_crs to ensure_projected() to choose a projection ",
             "suited to your extent."),
      span_lon
    )
    return(sf::st_crs(3857))
  }

  # as.integer() is belt-and-braces: floor() already yields an integral double,
  # which sprintf("%d", ...) accepts, but making the type explicit removes the
  # dependency on that coercion in the log messages below.
  cand_zone <- as.integer(max(1, min(60, floor((lon + 180) / 6) + 1)))
  cand_cm   <- 6 * cand_zone - 183           # central meridian of that zone
  lon_off   <- max(abs(lon_min - cand_cm), abs(lon_max - cand_cm))

  if (is.finite(lon_off) && lon_off > 5) {
    # Albers suits mid-latitude extents; its standard parallels are
    # conventionally placed at the 1/6 and 5/6 points of the latitude span.
    # Lambert Azimuthal Equal Area covers equatorial and high-latitude cases,
    # where Albers standard parallels degenerate.
    if (abs(lat) > 20 && abs(lat) < 70) {
      lat1 <- as.numeric(bb["ymin"]) + span_lat / 6
      lat2 <- as.numeric(bb["ymax"]) - span_lat / 6
      if (!is.finite(lat1) || !is.finite(lat2) || isTRUE(all.equal(lat1, lat2))) {
        lat1 <- lat - 5; lat2 <- lat + 5
      }
      .log_warn(
        paste0(".pick_local_projected_crs(): extent reaches %.1f deg from the ",
               "central meridian of UTM zone %d, well beyond the 3 deg the ",
               "zone is designed for (%.1f deg longitude span in total). Using ",
               "Albers equal-area (lat_1=%.1f, lat_2=%.1f, lon_0=%.1f) instead; ",
               "a single UTM zone would distort distances by several percent, ",
               "which propagates into variogram ranges, block sizes, GWR ",
               "bandwidth and GP length-scales. Pass target_crs to ",
               "ensure_projected() to override."),
        lon_off, cand_zone, span_lon, lat1, lat2, lon
      )
      return(sf::st_crs(sprintf(
        "+proj=aea +lat_1=%f +lat_2=%f +lat_0=%f +lon_0=%f +datum=WGS84 +units=m +no_defs",
        lat1, lat2, lat, lon)))
    }

    .log_warn(
      paste0(".pick_local_projected_crs(): extent reaches %.1f deg from the ",
             "central meridian of UTM zone %d, well beyond the 3 deg the zone ",
             "is designed for (%.1f deg longitude span in total). Using Lambert ",
             "Azimuthal Equal Area centred on (%.1f, %.1f) instead; a single ",
             "UTM zone would distort distances by several percent, which ",
             "propagates into variogram ranges, block sizes, GWR bandwidth and ",
             "GP length-scales. Pass target_crs to ensure_projected() to ",
             "override."),
      lon_off, cand_zone, span_lon, lon, lat
    )
    return(sf::st_crs(sprintf(
      "+proj=laea +lat_0=%f +lon_0=%f +datum=WGS84 +units=m +no_defs",
      lat, lon)))
  }

  # Reuse the candidate zone computed above so the zone named in the messages
  # and the EPSG code returned here can never diverge.
  epsg <- if (!is.na(lat) && lat < 0) 32700 + cand_zone else 32600 + cand_zone
  sf::st_crs(epsg)
}

# -----------------------------------------------------------------------------
# Projection Enforcement
# -----------------------------------------------------------------------------

#' Ensure an object has a projected CRS (with sensible defaults)
#'
#' Coerces spatial objects to a projected coordinate reference system suitable
#' for distance/area calculations.
#'
#' @details
#' An object that already has a projected CRS is returned untouched. Only
#' geographic (lon/lat) input is transformed, and the CRS chosen depends on the
#' extent of the data — it is **not** always UTM:
#'
#' \describe{
#'   \item{Local extents}{The UTM zone containing the data's centre
#'     (EPSG:326xx north of the equator, EPSG:327xx south). Distances and areas
#'     are close to true over a few degrees of longitude, which is the case
#'     this package is usually in.}
#'   \item{Wide extents}{Once the data reach well beyond the roughly 3 degrees
#'     a UTM zone is designed for, a single zone would distort distances by
#'     several percent — and that error propagates straight into variogram
#'     ranges, block sizes, GWR bandwidths and GP length-scales. An equal-area
#'     projection centred on the data is used instead: Albers conic
#'     (`+proj=aea`) for extents wider than tall, Lambert azimuthal
#'     (`+proj=laea`) otherwise. A warning names the projection, the span that
#'     triggered it, and this argument.}
#'   \item{Missing CRS}{If `x` has no CRS at all but its bounding box looks
#'     like lon/lat, EPSG:4326 is assumed with a warning, then the rules above
#'     apply. Set the CRS explicitly to suppress it.}
#' }
#'
#' `target_crs` overrides all of this: pass it whenever you need a specific,
#' reproducible projection — comparing runs, matching an existing layer, or
#' fixing the units that [make_folds()]'s `block_size` will be interpreted in.
#'
#' @param x An sf or sfc object (other objects returned unchanged).
#' @param target_crs Optional target CRS (sf object, integer EPSG, or crs).
#'   Must resolve to a usable CRS via [sf::st_crs()]; an unusable value (one
#'   that resolves to `NA_crs_`) raises an error rather than silently leaving
#'   `x` unprojected.
#' @return x, potentially with a new projected CRS.
#' @examples
#' library(sf)
#' pts_ll <- st_as_sf(
#'   data.frame(lon = c(9.1, 9.2), lat = c(48.7, 48.8)),
#'   coords = c("lon", "lat"), crs = 4326
#' )
#' # A local extent gets the containing UTM zone.
#' st_crs(ensure_projected(pts_ll))$epsg  # 32632
#'
#' # A continental extent gets an equal-area projection instead, with a
#' # warning naming it -- see Details.
#' wide <- st_as_sf(
#'   data.frame(lon = c(-120, -70), lat = c(30, 48)),
#'   coords = c("lon", "lat"), crs = 4326
#' )
#' st_crs(ensure_projected(wide))$proj4string
#'
#' # target_crs overrides the choice entirely.
#' st_crs(ensure_projected(pts_ll, target_crs = 3035))$epsg  # 3035
#' @export
ensure_projected <- function(x, target_crs = NULL) {
  if (!(inherits(x, "sf") || inherits(x, "sfc"))) return(x)

  if (!is.null(target_crs)) {
    tcrs <- sf::st_crs(target_crs)
    if (is.na(tcrs)) {
      # Silently returning `x` here would make ensure_projected() a no-op and
      # let unprojected coordinates flow into distance/area computations.
      stop(sprintf(
        paste0("ensure_projected(): `target_crs` does not resolve to a usable ",
               "CRS (sf::st_crs() returned NA for an object of class %s). Pass ",
               "an EPSG code, a proj string, an sf::crs object, or an sf/sfc ",
               "object that carries a CRS -- or omit `target_crs` to let a ",
               "local projected CRS be chosen automatically."),
        paste(class(target_crs), collapse = "/")
      ), call. = FALSE)
    }
    xcrs <- sf::st_crs(x)
    if (is.na(xcrs)) {
      # No source CRS: we cannot reproject, only assume.  Make the
      # assumption loud so silently misplaced data is traceable.
      .log_warn(
        "ensure_projected(): input has no CRS; stamping the supplied target CRS ('%s') WITHOUT reprojection. Verify the coordinates are already expressed in that CRS, or set the input CRS explicitly with sf::st_crs().",
        tcrs$input %||% "unknown"
      )
      sf::st_crs(x) <- tcrs
      return(x)
    }
    if (!identical(xcrs, tcrs)) x <- sf::st_transform(x, tcrs)
    return(x)
  }

  xcrs <- sf::st_crs(x)
  if (is.na(xcrs)) {
    bb <- try(sf::st_bbox(x), silent = TRUE)
    if (!inherits(bb, "try-error") &&
        all(is.finite(bb)) &&
        bb["xmin"] >= -180 && bb["xmax"] <= 180 &&
        bb["ymin"] >= -90  && bb["ymax"] <= 90) {

      ## Secondary heuristic: guard against small projected coordinates
      ## (e.g. a local site survey in metres with coords in [0, 50]) that
      ## happen to fall inside the lon/lat envelope.  We require that the
      ## data look *positively* like geographic coordinates:
      ##
      ##  (a) at least one coordinate has meaningful fractional degrees
      ##      (|frac| > 0.001 and < 0.999), OR
      ##  (b) the bounding-box extent exceeds 1 degree in at least one
      ##      axis — too large to be a typical small-site survey.
      x_range <- bb["xmax"] - bb["xmin"]
      y_range <- bb["ymax"] - bb["ymin"]
      has_large_extent <- (x_range > 1) || (y_range > 1)

      has_geo_precision <- FALSE
      if (!has_large_extent) {
        coords_mat <- try(sf::st_coordinates(x), silent = TRUE)
        if (!inherits(coords_mat, "try-error") && nrow(coords_mat) > 0) {
          frac <- abs(c(coords_mat[, 1], coords_mat[, 2]) %% 1)
          has_geo_precision <- any(frac > 0.001 & frac < 0.999)
        }
      }

      if (has_large_extent || has_geo_precision) {
        .log_warn(
          "ensure_projected(): CRS is missing; assuming EPSG:4326 because bbox looks like lon/lat (xmin=%.2f, xmax=%.2f, ymin=%.2f, ymax=%.2f). Set CRS explicitly to suppress this warning.",
          bb["xmin"], bb["xmax"], bb["ymin"], bb["ymax"]
        )
        sf::st_crs(x) <- sf::st_crs(4326)
        tr <- .pick_local_projected_crs(x)
        x <- sf::st_transform(x, tr)
      } else {
        .log_warn(
          "ensure_projected(): CRS is missing and coordinates fall within the lon/lat envelope, but they lack the decimal precision or extent typical of geographic coordinates (xmin=%.2f, xmax=%.2f, ymin=%.2f, ymax=%.2f). Not assuming EPSG:4326. Set CRS explicitly with sf::st_crs(x) <- sf::st_crs(...).",
          bb["xmin"], bb["xmax"], bb["ymin"], bb["ymax"]
        )
      }
    }
    return(x)
  }

  if (.is_longlat(x)) {
    tr <- .pick_local_projected_crs(x)
    x <- sf::st_transform(x, tr)
  }
  x
}

# -----------------------------------------------------------------------------
# CRS Harmonization
# -----------------------------------------------------------------------------

#' Harmonize CRS between two spatial objects
#'
#' Aligns two sf objects to a common CRS.
#'
#' When one input carries no CRS, the other object's CRS (or `target_crs`) is
#' *stamped* onto it with [sf::st_set_crs()]: the coordinates are **not**
#' reprojected, only relabelled.  That assumption is announced with a warning,
#' matching what [ensure_projected()] does in the same situation.
#'
#' @param a,b Objects of class sf or sfc.
#' @param prefer Which object's CRS to keep ("a" or "b").
#' @param target_crs Optional target CRS to apply to both.
#' @param on_transform_error What to do when st_transform() fails:
#'   \code{"stop"} (default) raises an error immediately;
#'   \code{"set_crs"} falls back to st_set_crs() (UNSAFE — coordinates are
#'   NOT reprojected, only the CRS label is overwritten). The \code{"set_crs"}
#'   option exists only for rare edge cases where you are certain the
#'   coordinates already match the target CRS definition.
#' @return A named list with components a and b.
#' @export
harmonize_crs <- function(a, b, prefer = c("a", "b"), target_crs = NULL,
                          on_transform_error = c("stop", "set_crs")) {
  if (!inherits(a, c("sf", "sfc"))) stop("harmonize_crs(): `a` must be sf or sfc.")
  if (!inherits(b, c("sf", "sfc"))) stop("harmonize_crs(): `b` must be sf or sfc.")
  prefer <- match.arg(prefer)
  on_transform_error <- match.arg(on_transform_error)

  crs_a <- sf::st_crs(a)
  crs_b <- sf::st_crs(b)

  .safe_transform <- function(x, to) {
    tryCatch(
      sf::st_transform(x, to),
      error = function(e) {
        if (identical(on_transform_error, "set_crs")) {
          .log_warn(
            "harmonize_crs(): st_transform() failed (%s); falling back to st_set_crs(). WARNING: coordinates are NOT reprojected \u2014 downstream distances, joins, and areas may be wrong. Set on_transform_error='stop' (the default) to surface this error instead.",
            conditionMessage(e)
          )
          sf::st_set_crs(x, to)
        } else {
          stop(sprintf(
            "harmonize_crs(): st_transform() failed: %s. If you are certain the coordinates already match the target CRS, pass on_transform_error='set_crs' to override (not recommended).",
            conditionMessage(e)
          ), call. = FALSE)
        }
      }
    )
  }

  # st_set_crs() only stamps a label: the coordinates are NOT moved.  Whenever
  # that happens on CRS-less input, say so, matching the loud assumption
  # ensure_projected() makes in the same situation.
  .warn_stamped <- function(which_obj, to) {
    lbl <- tryCatch({
      inp <- sf::st_crs(to)$input
      if (is.null(inp) || length(inp) != 1L || is.na(inp) ||
          !nzchar(as.character(inp))) "unknown" else as.character(inp)
    }, error = function(e) "unknown")
    .log_warn(
      "harmonize_crs(): `%s` has no CRS; stamping '%s' onto it WITHOUT reprojection. Verify the coordinates are already expressed in that CRS, or set it explicitly with sf::st_crs().",
      which_obj, lbl
    )
  }

  if (!is.null(target_crs)) {
    if (is.na(crs_a)) {
      .warn_stamped("a", target_crs)
      a <- sf::st_set_crs(a, target_crs)
    } else {
      a <- .safe_transform(a, target_crs)
    }
    if (is.na(crs_b)) {
      .warn_stamped("b", target_crs)
      b <- sf::st_set_crs(b, target_crs)
    } else {
      b <- .safe_transform(b, target_crs)
    }
    return(list(a = a, b = b))
  }

  if (is.na(crs_a) && is.na(crs_b)) return(list(a = a, b = b))

  if (is.na(crs_a) && !is.na(crs_b)) {
    .warn_stamped("a", crs_b)
    a <- sf::st_set_crs(a, crs_b)
    return(list(a = a, b = b))
  }
  if (!is.na(crs_a) && is.na(crs_b)) {
    .warn_stamped("b", crs_a)
    b <- sf::st_set_crs(b, crs_a)
    return(list(a = a, b = b))
  }

  if (identical(sf::st_crs(a), sf::st_crs(b))) return(list(a = a, b = b))

  if (prefer == "a") {
    b <- .safe_transform(b, crs_a)
  } else {
    a <- .safe_transform(a, crs_b)
  }
  list(a = a, b = b)
}

# -----------------------------------------------------------------------------
# Geometry Helpers
# -----------------------------------------------------------------------------

#' Fast center points of per-feature bounding boxes
#'
#' @param x An sf object.
#' @return An sfc (POINT) vector.
#' @keywords internal
#' @noRd
.bbox_center_sfc <- function(x) {
  stopifnot(inherits(x, "sf"))
  pts <- lapply(sf::st_geometry(x), function(g) {
    bb <- sf::st_bbox(g)
    sf::st_point(c((bb["xmin"] + bb["xmax"]) / 2, (bb["ymin"] + bb["ymax"]) / 2))
  })
  sf::st_sfc(pts, crs = sf::st_crs(x))
}

# -----------------------------------------------------------------------------
# Point Coercion
# -----------------------------------------------------------------------------

#' Coerce arbitrary geometries to representative points
#'
#' Converts the geometry column of an sf object to POINTs using one of several
#' strategies.
#'
#' LINESTRING midpoints are sampled with [sf::st_line_sample()], which yields
#' no point for an EMPTY LINESTRING.  Rather than silently misaligning the
#' result (or letting sf crash), such input raises an error; drop empty
#' geometries first with `x <- x[!sf::st_is_empty(x), ]`.
#'
#' @param x An sf object.
#' @param mode One of "auto", "centroid", "point_on_surface", "surface",
#'   "line_midpoint", "bbox_center".
#' @param tmp_project Logical; temporarily project for line-based midpoints.
#' @return An sf object with geometry coerced to POINTs.
#' @examples
#' library(sf)
#' poly <- st_sf(
#'   id = 1,
#'   geometry = st_sfc(st_polygon(list(rbind(
#'     c(0, 0), c(2, 0), c(2, 2), c(0, 2), c(0, 0)
#'   ))), crs = 32632)
#' )
#' coerce_to_points(poly, "auto")  # interior representative point
#' @export
coerce_to_points <- function(
    x,
    mode = c("auto", "centroid", "point_on_surface", "surface",
             "line_midpoint", "bbox_center"),
    tmp_project = TRUE
) {
  stopifnot(inherits(x, "sf"))
  mode <- match.arg(mode)
  if (identical(mode, "surface")) mode <- "point_on_surface"
  if (nrow(x) == 0L) return(x)

  g   <- sf::st_geometry(x)
  crs <- sf::st_crs(x)

  # -- bbox_center ---
  if (mode == "bbox_center") {
    return(sf::st_set_geometry(x, .bbox_center_sfc(x)))
  }

  # -- direct spherical-safe ops ---
  if (mode == "centroid") {
    return(sf::st_set_geometry(x, suppressWarnings(sf::st_centroid(g))))
  }
  if (mode == "point_on_surface") {
    return(sf::st_set_geometry(x, sf::st_point_on_surface(g)))
  }

  is_ll <- .is_longlat(x)

  # An EMPTY LINESTRING has no midpoint: st_line_sample() yields an empty
  # MULTIPOINT that st_cast(, "POINT") silently drops (and in sf 1.0.x the
  # call segfaults outright), so the sampled midpoints would no longer align
  # 1:1 with the rows they are scattered back into.  Reject before sampling.
  .guard_empty_lines <- function(geom, idx_ls) {
    empty <- which(sf::st_is_empty(geom))
    if (!length(empty)) return(invisible(NULL))
    shown <- idx_ls[empty][seq_len(min(5L, length(empty)))]
    stop(sprintf(
      paste0("coerce_to_points(): %d of %d LINESTRING feature(s) are EMPTY ",
             "(row(s) %s%s); st_line_sample() yields no midpoint for them, ",
             "which would misalign the result. Drop them first, e.g. ",
             "x <- x[!sf::st_is_empty(x), ]."),
      length(empty), length(idx_ls),
      paste(shown, collapse = ", "),
      if (length(empty) > 5L) ", ..." else ""
    ), call. = FALSE)
  }

  # Backstop for any other way the sampled count could diverge from the number
  # of LINESTRING rows being filled.
  .check_midpoint_alignment <- function(midps, idx_ls) {
    if (length(midps) == length(idx_ls)) return(invisible(NULL))
    stop(sprintf(
      paste0("coerce_to_points(): st_line_sample() returned %d midpoint(s) for ",
             "%d LINESTRING feature(s); the result would be misaligned."),
      length(midps), length(idx_ls)
    ), call. = FALSE)
  }

  if (mode == "line_midpoint") {
    gtypes <- as.character(sf::st_geometry_type(g, by_geometry = TRUE))
    if (any(gtypes %in% c("MULTILINESTRING", "GEOMETRYCOLLECTION"))) {
      stop("line_midpoint only supports LINESTRING; cast MULTILINESTRING first.")
    }
    idx_ls <- which(gtypes == "LINESTRING")
    idx_other <- which(gtypes != "LINESTRING")
    out <- vector("list", length(g))

    # Batch all LINESTRINGs: project once, sample all, back-transform once
    if (length(idx_ls)) {
      if (is_ll && !tmp_project) {
        ctr <- suppressWarnings(sf::st_centroid(g[idx_ls]))
        out[idx_ls] <- as.list(ctr)
      } else {
        g_ls      <- g[idx_ls]
        .guard_empty_lines(g_ls, idx_ls)
        g_ls_sf   <- sf::st_sf(geometry = g_ls)
        g_ls_proj <- if (tmp_project) ensure_projected(g_ls_sf) else g_ls_sf
        midps     <- sf::st_line_sample(sf::st_geometry(g_ls_proj), sample = 0.5)
        midps     <- sf::st_cast(midps, "POINT")
        if (!identical(sf::st_crs(g_ls_proj), crs))
          midps <- sf::st_transform(midps, crs)
        .check_midpoint_alignment(midps, idx_ls)
        out[idx_ls] <- as.list(midps)
      }
    }
    # Non-LINESTRING fallback to centroid
    if (length(idx_other)) {
      ctr <- suppressWarnings(sf::st_centroid(g[idx_other]))
      out[idx_other] <- as.list(ctr)
    }
    return(sf::st_set_geometry(x, sf::st_sfc(out, crs = crs)))
  }

  # ---- mode == "auto" ----
  # Batch by geometry type for performance instead of per-feature loop
  g   <- sf::st_geometry(x)
  out <- vector("list", length(g))
  gtypes <- as.character(sf::st_geometry_type(g, by_geometry = TRUE))

  # --- POINT: passthrough ---
  idx_pt <- which(gtypes == "POINT")
  if (length(idx_pt)) out[idx_pt] <- as.list(g[idx_pt])

  # --- MULTIPOINT: centroid of sub-points ---
  idx_mpt <- which(gtypes == "MULTIPOINT")
  if (length(idx_mpt)) {
    for (i in idx_mpt) {
      out[[i]] <- suppressWarnings(sf::st_centroid(g[i]))[[1L]]
    }
  }

  # --- POLYGON / MULTIPOLYGON: vectorized point_on_surface ---
  idx_poly <- which(gtypes %in% c("POLYGON", "MULTIPOLYGON"))
  if (length(idx_poly)) {
    pos <- sf::st_point_on_surface(g[idx_poly])
    out[idx_poly] <- as.list(pos)
  }

  # --- LINESTRING: midpoint via line_sample (batched) ---
  idx_ls <- which(gtypes == "LINESTRING")
  if (length(idx_ls)) {
    if (is_ll && !tmp_project) {
      ctr <- suppressWarnings(sf::st_centroid(g[idx_ls]))
      out[idx_ls] <- as.list(ctr)
    } else {
      g_ls     <- g[idx_ls]
      .guard_empty_lines(g_ls, idx_ls)
      g_ls_sf  <- sf::st_sf(geometry = g_ls)
      g_ls_proj <- if (tmp_project) ensure_projected(g_ls_sf) else g_ls_sf
      midps    <- sf::st_line_sample(sf::st_geometry(g_ls_proj), sample = 0.5)
      midps    <- sf::st_cast(midps, "POINT")
      if (!identical(sf::st_crs(g_ls_proj), crs))
        midps <- sf::st_transform(midps, crs)
      .check_midpoint_alignment(midps, idx_ls)
      out[idx_ls] <- as.list(midps)
    }
  }

  # --- MULTILINESTRING: longest part midpoint (batched projection) ---
  idx_mls <- which(gtypes == "MULTILINESTRING")
  if (length(idx_mls)) {
    if (is_ll && !tmp_project) {
      ctr <- suppressWarnings(sf::st_centroid(g[idx_mls]))
      out[idx_mls] <- as.list(ctr)
    } else {
      g_mls     <- g[idx_mls]
      g_mls_sf  <- sf::st_sf(geometry = g_mls)
      g_mls_proj <- if (tmp_project) ensure_projected(g_mls_sf) else g_mls_sf
      proj_geom  <- sf::st_geometry(g_mls_proj)
      proj_crs   <- sf::st_crs(g_mls_proj)
      for (j in seq_along(idx_mls)) {
        parts <- suppressWarnings(sf::st_cast(proj_geom[j], "LINESTRING"))
        if (length(parts) == 0L) {
          out[[idx_mls[j]]] <- suppressWarnings(sf::st_centroid(g[idx_mls[j]]))[[1]]
        } else {
          lens <- as.numeric(sf::st_length(parts))
          k    <- if (length(lens)) which.max(lens) else 1L
          mp   <- sf::st_line_sample(parts[k], sample = 0.5)
          mp   <- sf::st_cast(mp, "POINT")
          if (!identical(proj_crs, crs)) mp <- sf::st_transform(mp, crs)
          out[[idx_mls[j]]] <- mp[[1]]
        }
      }
    }
  }

  # --- Anything else: centroid fallback ---
  idx_other <- which(!gtypes %in% c("POINT", "MULTIPOINT", "POLYGON",
                                     "MULTIPOLYGON", "LINESTRING", "MULTILINESTRING"))
  if (length(idx_other)) {
    ctr <- suppressWarnings(sf::st_centroid(g[idx_other]))
    out[idx_other] <- as.list(ctr)
  }

  sf::st_set_geometry(x, sf::st_sfc(out, crs = crs))
}
