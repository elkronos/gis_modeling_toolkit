# -----------------------------------------------------------------------------
# CRS Selection
# -----------------------------------------------------------------------------

#' Evaluate an expression with sf's spherical engine (s2) switched on
#'
#' On lon/lat data sf hands st_area(), st_centroid(), st_union() and
#' st_sample() to s2 when \code{sf::sf_use_s2()} is TRUE.  With it FALSE,
#' areas and sampling need lwgeom, which this package does not depend on
#' ("package lwgeom required" from every Voronoi tessellation, even a
#' projected one, because the stable-ID sort key is measured in lon/lat), and
#' centroids and unions become planar arithmetic on degrees, which moves the
#' centre that picks a UTM zone and prints sf's warnings past \code{quiet}.
#' The package's own measurements on the sphere therefore run with s2 on
#' whatever the session says, and the session's setting is restored on the
#' way out, error or not.  Nothing is toggled when s2 is already on.
#'
#' @param expr Expression to evaluate.
#' @return The value of \code{expr}.
#' @keywords internal
#' @noRd
.with_s2 <- function(expr) {
  if (isTRUE(sf::sf_use_s2())) return(expr)
  suppressMessages(sf::sf_use_s2(TRUE))
  on.exit(suppressMessages(sf::sf_use_s2(FALSE)), add = TRUE)
  expr
}


#' Centre of a lon/lat layer on the sphere
#'
#' The centroid s2 gives for the union of the layer, which is what
#' \code{.pick_local_projected_crs()} has always used with s2 on (the
#' default).  Three things went wrong when it was taken with a bare
#' \code{st_centroid(st_union())}: with \code{sf_use_s2(FALSE)} the centre
#' was planar in degrees, so the UTM zone chosen for data near a zone edge
#' depended on a session option; sf's warning and message about that got past
#' \code{quiet}; and with s2 on, a polygon GEOS accepts but s2 rejects (a
#' repeated vertex, common in real shapefiles) stopped
#' \code{ensure_projected()} with "Edge 1 is degenerate" although a plain
#' \code{st_transform()} of it works.  Now s2 is always used; a geometry it
#' rejects is repaired and tried again; and if that fails too, the centre is
#' the normalised mean of the vertices' unit vectors (for a point layer
#' exactly what s2 returns), which needs no valid geometry at all.
#'
#' @param x_ll An sf/sfc object in a geographic CRS.
#' @return A 1 x 2 matrix (\code{X}, \code{Y}) in degrees, possibly
#'   non-finite when no centre exists (the caller falls back then).
#' @keywords internal
#' @noRd
.lonlat_centre <- function(x_ll) {
  g <- sf::st_geometry(x_ll)
  centre_of <- function(geom) {
    ctr <- sf::st_coordinates(sf::st_centroid(sf::st_union(geom)))
    if (!is.numeric(ctr) || length(ctr) < 2L) stop("no centroid")
    ctr[1L, 1:2, drop = FALSE]
  }
  ctr <- tryCatch(.with_s2(centre_of(g)), error = function(e) NULL)
  if (is.null(ctr))
    ctr <- tryCatch(.with_s2(centre_of(.safe_make_valid(g))), error = function(e) NULL)
  if (is.null(ctr)) {
    xy <- tryCatch(sf::st_coordinates(g)[, 1:2, drop = FALSE],
                   error = function(e) matrix(numeric(0), 0L, 2L))
    xy <- xy[is.finite(xy[, 1L]) & is.finite(xy[, 2L]), , drop = FALSE]
    lam <- xy[, 1L] * pi / 180; phi <- xy[, 2L] * pi / 180
    v   <- c(sum(cos(phi) * cos(lam)), sum(cos(phi) * sin(lam)), sum(sin(phi)))
    ctr <- matrix(if (nrow(xy) && sqrt(sum(v^2)) > 1e-12)
                    c(atan2(v[2L], v[1L]), atan2(v[3L], sqrt(v[1L]^2 + v[2L]^2))) * 180 / pi
                  else c(NA_real_, NA_real_),
                  1L, 2L, dimnames = list(NULL, c("X", "Y")))
  }
  ctr
}


#' Pick a sensible local projected CRS for an sf/sfc object
#'
#' Chooses an appropriate projected coordinate reference system for spatial
#' data, favouring a UTM zone based on the dataset's geographic centroid when
#' the extent is narrow enough for one.
#'
#' Beyond about 5 degrees from the central meridian of the candidate zone, the
#' zone, a Lambert azimuthal equal-area centred on the data and (unless its
#' cone constant degenerates) an Albers conic are each scored by
#' \code{.crs_distance_error()}, which projects representative points of the
#' data and compares planar with geodesic pairwise distances.  The least
#' distorting is returned, which may still be the zone.  Forcing continental
#' data into one UTM zone produces percent-scale distance errors that propagate
#' silently into variogram ranges, block sizes, GWR bandwidth and GP
#' length-scales; a badly chosen equal-area projection does the same, which is
#' why the choice is measured rather than assumed.  Only longitude offset
#' triggers the comparison: transverse Mercator error scales with distance from
#' the central meridian, so tall narrow north-south extents keep their zone
#' without being scored.
#'
#' Data straddling the antimeridian are detected from the one very large gap in
#' the sorted longitudes and given an equal-area projection centred on the true
#' extent, since Web Mercator (EPSG:3857) would SPLIT a wrapped layer.  Data
#' that span more than 180 degrees of longitude with no such gap surround a
#' pole; when every point also lies on one side of the equator, the layer
#' circles that pole (Antarctic stations, a pan-Arctic network) and gets a
#' Lambert azimuthal equal-area centred on it, provided that measures a
#' smaller distance error than the global fallback.  Only coverage that is
#' left -- spanning both hemispheres, or a low-latitude belt the polar
#' projection fits worse -- falls back to Web Mercator (Equal Earth for
#' \code{purpose = "area"}).
#'
#' @param x An sf or sfc object.
#' @return A list with \code{crs} (the chosen \code{sf::crs}) and
#'   \code{candidates}, a data.frame with one row per projection considered:
#'   \code{name}, \code{crs} (its definition as a string),
#'   \code{distance_error} (the measured worst-case relative distance error
#'   over sampled pairs, \code{NA} where it could not be measured) and
#'   \code{chosen}.  Where only one projection was in play (a zone kept on
#'   a local extent, the equal-area projection for a wrapped layer), that
#'   one is measured and reported alone; a layer circling a pole reports the
#'   polar projection and the global one it was measured against.
#'   \code{candidates} is \code{NULL} only where no local projection was
#'   chosen: non-geographic input, no finite centroid, or an extent that
#'   falls back to the global projection.
#' @keywords internal
#' @noRd
.pick_local_projected_crs <- function(x, purpose = c("distance", "area")) {
  purpose <- match.arg(purpose)
  if (!inherits(x, c("sf", "sfc"))) return(list(crs = sf::NA_crs_, candidates = NULL))
  crs <- sf::st_crs(x)
  if (!is.na(crs) && !.is_longlat(x)) return(list(crs = crs, candidates = NULL))
  # The measured candidates travel back with the choice.  On the comparison
  # paths the warning quotes two of the figures and drops the rest; here
  # every candidate's error is kept, since the number that propagates into
  # ranges, block sizes, bandwidths and length-scales is worth having whether
  # or not it was large enough to change the choice.

  x_ll <- tryCatch({
    if (is.na(crs)) stop("No CRS set.")
    sf::st_transform(x, 4326)
  }, error = function(e) x)

  scored <- function(crs_obj, names, crs_list, err, best) {
    list(crs = crs_obj, candidates = data.frame(
      name = names,
      crs = vapply(crs_list, function(cc) cc$input %||% NA_character_, character(1)),
      distance_error = as.numeric(err),
      chosen = seq_along(names) == best,
      stringsAsFactors = FALSE))
  }

  # The global fallback: Web Mercator when distances are what matter (and
  # nothing local fits), Equal Earth when areas are -- Mercator's area
  # distortion is unbounded, while Equal Earth is equal-area by construction.
  global_crs <- function() {
    if (purpose == "area") {
      ee <- tryCatch(sf::st_crs("+proj=eqearth +datum=WGS84 +units=m +no_defs"),
                     error = function(e) sf::NA_crs_)
      if (!is.na(ee)) return(ee)
      return(sf::st_crs("+proj=moll +datum=WGS84 +units=m +no_defs"))
    }
    sf::st_crs(3857)
  }
  if (is.na(sf::st_crs(x_ll)) || !.is_longlat(x_ll))
    return(list(crs = global_crs(), candidates = NULL))

  # On the sphere whatever sf_use_s2() says, and without failing on a polygon
  # s2 rejects: see .lonlat_centre().
  ctr <- .lonlat_centre(x_ll)
  if (!is.numeric(ctr) || length(ctr) < 2) return(list(crs = global_crs(), candidates = NULL))
  lon <- ctr[1]; lat <- ctr[2]

  # st_centroid() on empty or degenerate geometry can return NA/NaN, which
  # passes the is.numeric() check above.  A non-finite centroid cannot place
  # either a UTM zone or an equal-area projection, so fall back explicitly
  # rather than letting NA propagate into the comparisons below.
  if (!is.finite(lon) || !is.finite(lat)) {
    .log_warn(
      ".pick_local_projected_crs(): could not compute a finite centroid (lon = %s, lat = %s); falling back to %s.",
      format(lon), format(lat), if (purpose == "area") "Equal Earth" else "EPSG:3857"
    )
    return(list(crs = global_crs(), candidates = NULL))
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

  # A bbox wider than a hemisphere is either global coverage or data merely
  # straddling the antimeridian.  The BOX cannot tell them apart -- but the
  # coordinates can: if the data wrap, the sorted longitudes contain one gap
  # far larger than any other (the empty stretch the box spans the wrong way
  # round), and the true extent is 360 minus that gap.  This matters because
  # EPSG:3857 is the one answer that is worse than doing nothing for wrapped
  # data: it SPLITS the layer, so two stations 41 km apart come out 40,068 km
  # apart and every distance downstream -- variogram range, block size, GWR
  # bandwidth, GP length-scale -- is destroyed.
  if (is.finite(span_lon) && span_lon > 180) {
    lons  <- sf::st_coordinates(sf::st_geometry(x_ll))[, 1L]
    lons  <- sort(lons[is.finite(lons)])
    wrapped <- FALSE
    if (length(lons) > 1L) {
      gaps     <- diff(lons)
      max_gap  <- max(gaps)
      # True span once the wrap is undone.  Treat it as wrapped only when
      # that span is comfortably smaller than a hemisphere, so genuinely
      # global coverage (many small gaps) still falls through.
      true_span <- 360 - max_gap
      wrapped   <- is.finite(max_gap) && true_span <= 180
    }

    if (wrapped) {
      # Re-express longitudes on [0, 360) and recurse on the shifted layer, so
      # the zone/equal-area logic below sees a contiguous extent.  The CRS it
      # picks is expressed with an explicit +lon_0/+pm, so transforming the
      # ORIGINAL (unshifted) coordinates into it is correct: PROJ wraps
      # longitudes itself.
      lon_shift <- lons[lons < 0] + 360
      lon_all   <- c(lons[lons >= 0], lon_shift)
      lon_ctr   <- mean(range(lon_all))
      lon_ctr   <- ((lon_ctr + 180) %% 360) - 180   # back onto [-180, 180)
      .log_warn(
        paste0(".pick_local_projected_crs(): the bounding box spans %.1f deg of ",
               "longitude, but the coordinates straddle the antimeridian: the ",
               "true extent is %.1f deg. Using an equal-area projection centred ",
               "on lon_0=%.1f. (EPSG:3857 would have split the layer in two.)"),
        span_lon, 360 - max(diff(lons)), lon_ctr
      )
      wrap_crs <- sf::st_crs(sprintf(
        "+proj=laea +lat_0=%f +lon_0=%f +datum=WGS84 +units=m +no_defs",
        lat, lon_ctr))
      return(scored(wrap_crs,
                    sprintf("Lambert azimuthal equal-area centred on (%.1f, %.1f)",
                            lon_ctr, lat),
                    list(wrap_crs), .crs_distance_error(x_ll, wrap_crs), 1L))
    }

    # No gap of 180 deg or more means the longitudes surround a pole.  With
    # every point on one side of the equator the layer circles THAT pole --
    # Antarctic stations, a pan-Arctic network -- and is not global coverage.
    # Web Mercator splits such a layer at +/-180 and stretches it towards the
    # pole: rings of Antarctic stations measured worst-case distance errors
    # of 15,000-20,000% in it (a 111 km pair came out 445 km, the South Pole
    # at y = -2.4e8 m), where a Lambert azimuthal centred on the pole gave
    # about 2%.  That projection is equal-area, so it serves purpose = "area"
    # too.  Its distortion grows away from the pole (about 40% for a belt
    # reaching the equator), so it is measured against the global fallback
    # and used only when it does better; data spanning both hemispheres never
    # get here and keep the global fallback as before.
    lat_min <- as.numeric(bb["ymin"]); lat_max <- as.numeric(bb["ymax"])
    if (is.finite(lat_min) && is.finite(lat_max) && (lat_min >= 0 || lat_max <= 0)) {
      north     <- lat_min >= 0
      polar_crs <- sf::st_crs(sprintf(
        "+proj=laea +lat_0=%d +lon_0=0 +datum=WGS84 +units=m +no_defs",
        if (north) 90L else -90L))
      glob      <- global_crs()
      glob_name <- if (purpose != "area") "Web Mercator (EPSG:3857)" else
        if (grepl("eqearth", glob$input, fixed = TRUE)) "Equal Earth" else "Mollweide"
      cands     <- list(
        list(name = sprintf("Lambert azimuthal equal-area centred on the %s Pole",
                            if (north) "North" else "South"),
             crs = polar_crs),
        list(name = glob_name, crs = glob))
      err <- vapply(cands, function(cd) .crs_distance_error(x_ll, cd$crs), numeric(1))
      if (is.finite(err[1L]) && (!is.finite(err[2L]) || err[1L] < err[2L])) {
        .log_warn(
          paste0(".pick_local_projected_crs(): longitude extent spans %.1f deg ",
                 "without straddling the antimeridian, and every point lies %s ",
                 "of the equator (latitude %.1f to %.1f): the layer circles the ",
                 "%s Pole. Using %s: measured worst-case distance error %.2f%% ",
                 "against %s for %s. Pass target_crs to ensure_projected() to ",
                 "override."),
          span_lon, if (north) "north" else "south", lat_min, lat_max,
          if (north) "North" else "South", cands[[1L]]$name, 100 * err[1L],
          if (is.finite(err[2L])) sprintf("%.2f%%", 100 * err[2L]) else "not measurable",
          cands[[2L]]$name)
        return(scored(polar_crs, vapply(cands, `[[`, character(1), "name"),
                      lapply(cands, `[[`, "crs"), err, 1L))
      }
    }

    .log_warn(
      paste0(".pick_local_projected_crs(): longitude extent spans %.1f deg and ",
             "the coordinates do not straddle the antimeridian, so this is ",
             "global coverage; no local projection fits it. Falling back to ",
             "%s. Pass target_crs to ensure_projected() to choose a ",
             "projection suited to your extent."),
      span_lon, if (purpose == "area") "Equal Earth (equal-area)" else "EPSG:3857"
    )
    return(list(crs = global_crs(), candidates = NULL))
  }

  # as.integer() is belt-and-braces: floor() already yields an integral double,
  # which sprintf("%d", ...) accepts, but making the type explicit removes the
  # dependency on that coercion in the log messages below.
  cand_zone <- as.integer(max(1, min(60, floor((lon + 180) / 6) + 1)))
  cand_cm   <- 6 * cand_zone - 183           # central meridian of that zone
  lon_off   <- max(abs(lon_min - cand_cm), abs(lon_max - cand_cm))

  utm_epsg <- if (!is.na(lat) && lat < 0) 32700 + cand_zone else 32600 + cand_zone

  if (purpose == "area") {
    # Areas are what matter, so the zone -- conformal, not equal-area -- is
    # not a candidate at any extent: a UTM zone's area scale runs from 0.9992
    # at the central meridian to about 1.002 at the zone edge, and far worse
    # beyond it.  Both equal-area candidates preserve area exactly; the one
    # kept is the one that distorts distances least, because everything else
    # in the package -- ranges, block sizes, bandwidths -- still reads
    # distances off the same coordinates.
    lat1 <- as.numeric(bb["ymin"]) + span_lat / 6
    lat2 <- as.numeric(bb["ymax"]) - span_lat / 6
    if (!is.finite(lat1) || !is.finite(lat2) || isTRUE(all.equal(lat1, lat2))) {
      lat1 <- lat - 5; lat2 <- lat + 5
    }
    cone_n <- (sin(lat1 * pi / 180) + sin(lat2 * pi / 180)) / 2
    cands <- list(list(
      name = sprintf("Lambert azimuthal equal-area centred on (%.1f, %.1f)", lon, lat),
      crs = sf::st_crs(sprintf(
        "+proj=laea +lat_0=%f +lon_0=%f +datum=WGS84 +units=m +no_defs", lat, lon))))
    if (is.finite(cone_n) && abs(cone_n) >= 0.05) {
      cands <- c(cands, list(list(
        name = sprintf("Albers equal-area (lat_1=%.1f, lat_2=%.1f, lon_0=%.1f)",
                       lat1, lat2, lon),
        crs = sf::st_crs(sprintf(
          paste0("+proj=aea +lat_1=%f +lat_2=%f +lat_0=%f +lon_0=%f ",
                 "+datum=WGS84 +units=m +no_defs"),
          lat1, lat2, lat, lon)))))
    }
    err  <- vapply(cands, function(cd) .crs_distance_error(x_ll, cd$crs), numeric(1))
    best <- which.min(ifelse(is.finite(err), err, Inf))
    if (!length(best) || !is.finite(err[best])) best <- 1L
    .log_info(
      paste0(".pick_local_projected_crs(): purpose = \"area\": using %s ",
             "(equal-area; worst-case distance error %s over sampled pairs) ",
             "rather than UTM zone %d."),
      cands[[best]]$name,
      if (is.finite(err[best])) sprintf("%.2f%%", 100 * err[best]) else "not measurable",
      cand_zone
    )
    return(scored(cands[[best]]$crs, vapply(cands, `[[`, character(1), "name"),
                  lapply(cands, `[[`, "crs"), err, best))
  }

  if (is.finite(lon_off) && lon_off > 5) {
    # Which projection is actually best is not something a rule of thumb gets
    # right.  The previous heuristic -- centroid latitude picks Albers or
    # LAEA, Albers' standard parallels come from the bbox -- could return a
    # projection an ORDER OF MAGNITUDE worse than the zone it was rejecting:
    # a trans-equatorial extent collapses the conic's cone constant (15.6%
    # distance error against UTM's 1.7%, and at lat_1 = -lat_2 exactly PROJ
    # refuses the string and every caller aborts with "invalid crs"), and a
    # very tall extent defeats any conic whatever its parallels.
    #
    # So measure it instead of guessing.  Project a deterministic sample of
    # the data's own points into each candidate, compare the planar pairwise
    # distances with the geodesic ones, and keep whichever candidate distorts
    # least.  That is exactly the quantity the choice exists to protect --
    # estimate_sac_range() returns a range in CRS units, make_folds() sizes
    # blocks in CRS units, and both the GWR bandwidth and the GP length-scale
    # read projected coordinates.
    lat1 <- as.numeric(bb["ymin"]) + span_lat / 6
    lat2 <- as.numeric(bb["ymax"]) - span_lat / 6
    if (!is.finite(lat1) || !is.finite(lat2) || isTRUE(all.equal(lat1, lat2))) {
      lat1 <- lat - 5; lat2 <- lat + 5
    }
    cone_n <- (sin(lat1 * pi / 180) + sin(lat2 * pi / 180)) / 2

    cands <- list(
      list(name = sprintf("UTM zone %d", cand_zone), crs = sf::st_crs(utm_epsg)),
      list(name = sprintf("Lambert azimuthal equal-area centred on (%.1f, %.1f)",
                          lon, lat),
           crs = sf::st_crs(sprintf(
             "+proj=laea +lat_0=%f +lon_0=%f +datum=WGS84 +units=m +no_defs",
             lat, lon)))
    )
    # A conic whose parallels straddle the equator is not merely poor, it is
    # degenerate: |lat_1 + lat_2| = 0 is rejected by PROJ outright.  Leave it
    # out of the comparison rather than letting st_crs() raise.
    if (is.finite(cone_n) && abs(cone_n) >= 0.05) {
      cands <- c(cands, list(list(
        name = sprintf("Albers equal-area (lat_1=%.1f, lat_2=%.1f, lon_0=%.1f)",
                       lat1, lat2, lon),
        crs = sf::st_crs(sprintf(
          paste0("+proj=aea +lat_1=%f +lat_2=%f +lat_0=%f +lon_0=%f ",
                 "+datum=WGS84 +units=m +no_defs"),
          lat1, lat2, lat, lon)))))
    }

    err <- vapply(cands, function(cd) .crs_distance_error(x_ll, cd$crs),
                  numeric(1))
    best <- which.min(ifelse(is.finite(err), err, Inf))
    if (!length(best) || !is.finite(err[best])) best <- 1L   # UTM fallback

    if (identical(cands[[best]]$name, cands[[1L]]$name)) {
      # The figures are only quotable when the comparison actually ran; a
      # candidate that could not be scored comes back NA, and a message
      # reporting "NA% vs NA%" asserts a comparison that never happened.
      measured <- all(is.finite(err))
      .log_warn(
        paste0(".pick_local_projected_crs(): extent reaches %.1f deg from the ",
               "central meridian of UTM zone %d (%.1f deg of longitude in ",
               "total); the zone is kept%s. Pass target_crs to ",
               "ensure_projected() to override."),
        lon_off, cand_zone, span_lon,
        if (measured)
          sprintf(", because on this extent every equal-area candidate distorts distances more than it does (%.2f%% vs %.2f%% worst-case over sampled pairs)",
                  100 * min(err[-1L]), 100 * err[1L])
        else " because the candidate projections could not be scored on this layer"
      )
      return(scored(sf::st_crs(utm_epsg), vapply(cands, `[[`, character(1), "name"),
                    lapply(cands, `[[`, "crs"), err, 1L))
    }

    .log_warn(
      paste0(".pick_local_projected_crs(): extent reaches %.1f deg from the ",
             "central meridian of UTM zone %d (%.1f deg of longitude in ",
             "total). Using %s instead: measured worst-case distance error ",
             "%.2f%% against the zone's %.2f%%, and that error propagates into ",
             "variogram ranges, block sizes, GWR bandwidth and GP ",
             "length-scales. Pass target_crs to ensure_projected() to ",
             "override."),
      lon_off, cand_zone, span_lon, cands[[best]]$name,
      100 * err[best], 100 * err[1L]
    )
    return(scored(cands[[best]]$crs, vapply(cands, `[[`, character(1), "name"),
                  lapply(cands, `[[`, "crs"), err, best))
  }

  # Reuse the candidate zone computed above so the zone named in the messages
  # and the EPSG code returned here can never diverge.  The zone's own error
  # is measured here as well -- about 20 ms, bounded by the 40-point sample --
  # so the ordinary path reports the figure the comparison paths quote.
  zone <- sf::st_crs(utm_epsg)
  scored(zone, sprintf("UTM zone %d", cand_zone), list(zone),
         .crs_distance_error(x_ll, zone), 1L)
}


#' Worst-case relative distance error of a projection on a point layer
#'
#' Projects a deterministic sample of \code{x_ll} (lon/lat) into \code{crs} and
#' compares the planar pairwise distances with the geodesic ones, returning the
#' largest relative discrepancy.  Used by \code{.pick_local_projected_crs()} to
#' choose between candidate projections by measurement rather than by rule of
#' thumb.
#'
#' The sample is taken by evenly spaced index (no RNG), so the result is
#' reproducible; \code{max_n} keeps the pairwise work bounded.
#'
#' @param x_ll An sf/sfc object in a geographic CRS.  Non-POINT geometry is
#'   reduced to representative points first, so the two distance vectors are
#'   the same length (\code{st_coordinates()} yields one row per vertex).
#'   With fewer than \code{max_n} features the outline's vertices, densified
#'   along its edges, are added to those points, so that a single study-area
#'   polygon is measured across its extent rather than not at all.
#' @param crs Candidate \code{sf::crs}.
#' @param max_n Maximum number of points to sample.  Default 40 (780 pairs).
#' @return Numeric worst-case \code{|d_planar / d_geodesic - 1|}, or \code{NA}
#'   when it cannot be computed (which the caller treats as "unusable").
#' @keywords internal
#' @noRd
.crs_distance_error <- function(x_ll, crs, max_n = 40L) {
  tryCatch({
    g <- sf::st_geometry(x_ll)
    g <- g[!sf::st_is_empty(g)]
    # Representative POINTS, not the geometries themselves.  st_coordinates()
    # returns one row per VERTEX, so for a polygon or line layer the projected
    # distance vector was a different length from the geodesic one: the
    # comparison recycled, R raised "longer object length is not a multiple of
    # shorter object length" at the caller, every candidate scored NA, and the
    # selection silently fell back to the UTM zone while the log line reported
    # "NA% vs NA%".  Every county-polygon layer took that path.  The empty
    # parts go first: one inside a line feature segfaults GEOS here (see
    # .drop_empty_parts()), and tryCatch() cannot catch a crash, so any
    # lon/lat line layer with a null part took ensure_projected() down.
    if (!all(sf::st_geometry_type(g, by_geometry = TRUE) == "POINT")) {
      g_full <- .drop_empty_parts(g)
      g <- suppressWarnings(sf::st_point_on_surface(g_full))
      # One point per feature is nothing to measure on a study-area outline:
      # a single polygon gave one point, every candidate scored NA, and the
      # selector kept the UTM zone at any extent -- a CONUS outline got zone
      # 15 (13.7% worst-case error) where its vertices score Albers at 2.4%,
      # and prep_model_data(boundary =) moved a whole analysis into it.  A
      # handful of features gives a handful of pairs, none near the edges
      # where a zone distorts most.  So below `max_n` points, add the
      # outline's own vertices, densified along each edge until there are
      # about `max_n` of them (a box has only its four corners).  Layers with
      # `max_n` features or more are measured as before, and so is a
      # GEOMETRYCOLLECTION layer, whose vertices st_coordinates() refuses.
      xy <- if (length(g) < max_n)
        tryCatch(sf::st_coordinates(g_full), error = function(e) NULL)
      if (!is.null(xy)) {
        ring <- if (ncol(xy) > 2L)
          do.call(paste, as.data.frame(xy[, -(1:2), drop = FALSE])) else rep("1", nrow(xy))
        xy <- xy[, 1:2, drop = FALSE]
        ok <- is.finite(xy[, 1L]) & is.finite(xy[, 2L])
        xy <- xy[ok, , drop = FALSE]; ring <- ring[ok]
        if (nrow(xy) > 0L) {
          per_edge <- max(0L, ceiling(max_n / max(1L, nrow(unique(xy)))) - 1L)
          if (per_edge > 0L && nrow(xy) > 1L) {
            same <- ring[-1L] == ring[-length(ring)]
            a <- xy[-nrow(xy), , drop = FALSE][same, , drop = FALSE]
            b <- xy[-1L, , drop = FALSE][same, , drop = FALSE]
            f <- rep(seq_len(per_edge) / (per_edge + 1), each = nrow(a))
            a <- a[rep(seq_len(nrow(a)), per_edge), , drop = FALSE]
            b <- b[rep(seq_len(nrow(b)), per_edge), , drop = FALSE]
            xy <- rbind(xy, a + f * (b - a))
          }
          xy <- unique(xy)
          g <- c(sf::st_geometry(g), sf::st_geometry(sf::st_as_sf(
            data.frame(x = xy[, 1L], y = xy[, 2L]), coords = c("x", "y"),
            crs = sf::st_crs(g))))
        }
      }
    }
    n <- length(g)
    if (n < 2L) return(NA_real_)
    if (n > max_n) g <- g[unique(round(seq(1, n, length.out = max_n)))]
    # Ellipsoidal geodesics, not sf's s2 sphere: the sphere is itself
    # 0.24-0.56% off WGS84, which is the size of the errors being ranked.
    ll  <- sf::st_coordinates(sf::st_transform(g, 4326))[, 1:2, drop = FALSE]
    ij  <- which(lower.tri(matrix(0, nrow(ll), nrow(ll))), arr.ind = TRUE)
    d_geo <- .geod_distance(ll[ij[, 2], 1], ll[ij[, 2], 2],
                            ll[ij[, 1], 1], ll[ij[, 1], 2])
    d_prj <- as.numeric(stats::dist(
      sf::st_coordinates(sf::st_transform(g, crs))[, 1:2, drop = FALSE]))
    keep  <- is.finite(d_geo) & is.finite(d_prj) & d_geo > 0
    if (!any(keep)) return(NA_real_)
    max(abs(d_prj[keep] / d_geo[keep] - 1))
  }, error = function(e) NA_real_)
}


#' Worst-case relative area distortion of a projection over a layer
#'
#' The counterpart of \code{.crs_distance_error()} for the property a density
#' or a rate depends on.  Probe polygons are measured twice: their planar
#' area in \code{crs} and their geodesic area on the globe.  The probes are the
#' layer's own polygons when it has them (up to \code{max_n}, by evenly spaced
#' index), otherwise an \code{n x n} grid over its bounding box.  An equal-area
#' projection makes the ratio of the two the same for every probe (up to the
#' sphere-versus-ellipsoid factor of the s2 areas, which drifts slowly with
#' latitude and stays within a few tenths of a percent), so the figure
#' returned is the spread of that ratio: the largest relative departure from
#' its median.  A conformal projection such as a UTM zone gives a ratio that
#' varies as the square of its scale factor: 0.25 percent across a zone,
#' 14 percent when the conterminous United States is forced into one.
#'
#' @param x An sf/sfc object with a CRS.
#' @param crs The projection to score; default the CRS of \code{x}.
#' @param n Probe grid size when \code{x} has no polygons.
#' @param max_n Largest number of the layer's own polygons to measure.
#' @param grid Logical; probe with the \code{n x n} grid over the bounding box
#'   even when \code{x} has polygons.  A single study-area polygon is one
#'   probe, and one ratio has no spread to measure.
#' @return Numeric worst-case \code{|ratio / median(ratio) - 1|}, or
#'   \code{NA} when it cannot be computed (no CRS, no area, geodesic areas
#'   unavailable).
#' @keywords internal
#' @noRd
.crs_area_error <- function(x, crs = NULL, n = 6L, max_n = 200L, grid = FALSE) {
  tryCatch({
    if (is.null(crs)) crs <- sf::st_crs(x)
    if (is.na(crs) || is.na(sf::st_crs(x))) return(NA_real_)
    g <- sf::st_geometry(x)
    g <- g[!sf::st_is_empty(g)]
    if (!length(g)) return(NA_real_)
    types <- as.character(sf::st_geometry_type(g, by_geometry = TRUE))
    probe <- if (!isTRUE(grid) && all(types %in% c("POLYGON", "MULTIPOLYGON"))) {
      if (length(g) > max_n) g[unique(round(seq(1, length(g), length.out = max_n)))] else g
    } else {
      bb <- sf::st_bbox(sf::st_transform(g, crs))
      sf::st_make_grid(sf::st_as_sfc(bb), n = c(n, n), what = "polygons")
    }
    probe  <- sf::st_transform(probe, crs)
    # Densify before going to lon/lat: s2 reads each edge as a great circle,
    # and over a continental probe that is not the straight edge the planar
    # area was measured on -- an Equal Earth grid over a near-global extent
    # measured 10% "distortion".  A projected CRS only; densifying lon/lat
    # needs lwgeom.
    if (!isTRUE(sf::st_is_longlat(crs))) {
      pb  <- sf::st_bbox(probe)
      ext <- max(as.numeric(pb["xmax"] - pb["xmin"]), as.numeric(pb["ymax"] - pb["ymin"]))
      if (is.finite(ext) && ext > 0) probe <- sf::st_segmentize(probe, dfMaxLength = ext / 100)
    }
    planar <- as.numeric(sf::st_area(probe))
    # Geodesic areas on the sphere whatever sf_use_s2() says: with it off, sf
    # asks lwgeom for them, which is not a dependency, the error became NA
    # here, and summarize_by_cell(area = TRUE) then refused every grid while
    # ensure_projected(purpose = "area") skipped its distortion check.
    geod   <- .with_s2(as.numeric(sf::st_area(sf::st_transform(probe, 4326))))
    ratio  <- planar / geod
    ok <- is.finite(ratio) & is.finite(geod) & geod > 0
    if (sum(ok) < 2L) return(NA_real_)
    max(abs(ratio[ok] / stats::median(ratio[ok]) - 1))
  }, error = function(e) NA_real_)
}

#' The area-distortion tolerance below which a CRS is taken as equal-area
#'
#' One percent.  Measured with .crs_area_error(): a UTM zone edge to edge
#' 0.25%, a 2.5-degree extent inside one 0.04%, an equal-area projection a
#' few tenths of a percent (the s2 sphere against the WGS84 ellipsoid), the
#' conterminous US in one zone 13.9%, Web Mercator over 2.5 degrees of
#' latitude at 48N 4.1%.  The tolerance sits in the gap.
#' @keywords internal
#' @noRd
.area_error_tol <- 0.01

# -----------------------------------------------------------------------------
# Projection Enforcement
# -----------------------------------------------------------------------------

#' Return sampled points to the coordinate space of the input
#'
#' The line-midpoint branches project temporarily so that "halfway along the
#' line" is measured in a length unit, then bring the midpoints back.  When
#' the input had NO CRS and \code{ensure_projected()} interpreted it as
#' lon/lat, "back" is EPSG:4326 (the space the input's numbers were in), with
#' the CRS then stripped again so the output matches the input.
#' \code{sf::st_transform(x, NA_crs_)} is an error ("crs not found"), which
#' is what every CRS-less LINESTRING layer inside the lon/lat envelope used
#' to die with.
#'
#' @param pts sfc of points in the projected CRS.
#' @param proj_obj The object \code{ensure_projected()} returned (carries
#'   \code{attr(, "crs_assumed")} when an assumption was made).
#' @param crs The input's \code{sf::st_crs()}, possibly \code{NA}.
#' @return \code{pts} in the input's coordinate space.
#' @keywords internal
#' @noRd
.back_to_input_crs <- function(pts, proj_obj, crs) {
  pcrs <- sf::st_crs(proj_obj)
  if (identical(pcrs, crs)) return(pts)
  if (is.na(crs)) {
    assumed <- attr(proj_obj, "crs_assumed")
    if (!is.null(assumed) && !is.na(pcrs))
      pts <- sf::st_transform(pts, sf::st_crs(assumed))
    return(sf::st_set_crs(pts, NA))
  }
  sf::st_transform(pts, crs)
}


#' Do CRS-less coordinates look like longitude/latitude?
#'
#' The single heuristic behind every "assume EPSG:4326" decision in the
#' package, so that the decision is the SAME wherever it is taken.  It used
#' to live only in the no-target branch of \code{ensure_projected()}: fitting
#' assumed lon/lat and projected, while every \code{predict()} method (which
#' passes a target) stamped the fit's projected CRS onto the raw numbers.
#' The same CRS-less rows then sat in two different places at fit and at
#' predict time, and \code{predict(fit, newdata = training_rows)} disagreed
#' with \code{fitted(fit)} by up to the response's standard deviation.
#'
#' Coordinates are taken to be lon/lat when the bounding box fits the
#' [-180, 180] x [-90, 90] envelope AND the data look \emph{positively}
#' geographic: an extent above one degree on some axis, OR fractional
#' coordinates with the precision of decimal degrees.
#'
#' The two tests are a disjunction, and the extent test is evaluated first, so
#' any CRS-less planar layer that fits inside the envelope and is more than one
#' unit across is treated as lon/lat, a 50 m site survey included.  That is a
#' deliberate trade: a CRS-less layer is ambiguous by construction, and the
#' failure modes are not symmetric.  Reading true degrees as planar metres
#' makes every distance in the package meaningless with no way to notice;
#' reading a small planar survey as degrees produces coordinates that are
#' obviously wrong and a warning that names the assumption.  Requiring BOTH
#' tests would not help: a genuine study area 0.01 degrees across passes only
#' the precision test, and a [0, 1]-normalised planar layer passes it too.
#' When the data are planar, set the CRS explicitly, as the warning says.
#'
#' @param x An sf/sfc object with no CRS.
#' @return A list: \code{lonlat} (logical) and \code{bb} (the bbox, or
#'   \code{NULL} if it could not be taken).
#' @keywords internal
#' @noRd
.looks_like_lonlat <- function(x) {
  bb <- try(sf::st_bbox(x), silent = TRUE)
  if (inherits(bb, "try-error") || !all(is.finite(bb)))
    return(list(lonlat = FALSE, bb = NULL))
  in_env <- bb["xmin"] >= -180 && bb["xmax"] <= 180 &&
            bb["ymin"] >= -90  && bb["ymax"] <= 90
  if (!in_env) return(list(lonlat = FALSE, bb = bb))

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
  list(lonlat = isTRUE(has_large_extent || has_geo_precision), bb = bb)
}


#' Re-apply the CRS assumption a fit was built under to CRS-less new data
#'
#' \code{ensure_projected()} records \code{attr(, "crs_assumed")} on data it
#' interpreted as lon/lat.  A prediction frame with no CRS is given that same
#' interpretation before being aligned to the fit, so that newdata drawn from
#' the training rows (even a subset whose own bounding box would not have
#' triggered the heuristic) lands where the training rows did.
#'
#' @param newdata sf, possibly CRS-less.
#' @param training_sf The fit's \code{data_sf}.
#' @param caller Name for the log line.
#' @param what Name of the argument \code{newdata} came from, for the warning.
#'   \code{predict_surface()} replays the assumption on \code{grid} and
#'   \code{boundary} too, and a message naming \code{newdata} there sends the
#'   reader looking for an argument they did not pass.
#' @return \code{newdata}: with a CRS set when a POSITIVE assumption was
#'   replayed; marked \code{crs_assumed = "none"} when the fit itself had no
#'   CRS, so that \code{ensure_projected()} leaves it in the space the fit used;
#'   unchanged otherwise.
#' @keywords internal
#' @noRd
.replay_crs_assumption <- function(newdata, training_sf, caller = "predict",
                                   what = "newdata") {
  if (!inherits(newdata, "sf") && !inherits(newdata, "sfc")) return(newdata)
  if (!is.na(sf::st_crs(newdata))) return(newdata)
  assumed <- attr(training_sf, "crs_assumed")

  # A NEGATIVE decision has to be replayed too, and it is the harder half.
  # When the training data had no CRS and the heuristic declined to call them
  # lon/lat, the fit works in the raw input space.  Nothing was recorded, so
  # every predict() re-ran the heuristic on `newdata` ALONE -- and a subset of
  # those same training rows whose own bounding box happens to sit inside the
  # lon/lat envelope is judged differently from the whole: the subset is taken
  # as degrees, reprojected, and predicted about 1e6 m from where it was
  # fitted.  predict(fit, training_subset) then disagreed with
  # fitted(fit)[subset] by more than the response's standard deviation while
  # predict(fit, full_training_set) agreed exactly.  Mark the newdata so
  # ensure_projected() leaves it in the same untouched space the fit used.
  if ((is.null(assumed) || identical(assumed, "none")) &&
      is.na(sf::st_crs(training_sf))) {
    attr(newdata, "crs_assumed") <- "none"
    return(newdata)
  }

  if (is.null(assumed) || identical(assumed, "none")) return(newdata)
  .warn_and_log(paste0("%s(): `%s` has no CRS; interpreting it as %s, the ",
                   "assumption the model was fitted under, so that it is placed ",
                   "where the training data were. Set the CRS explicitly to ",
                   "suppress this."), caller, what, assumed)
  sf::st_set_crs(newdata, sf::st_crs(assumed))
}


#' Ensure an object has a projected CRS (with sensible defaults)
#'
#' Coerces spatial objects to a projected coordinate reference system suitable
#' for distance/area calculations.
#'
#' @details
#' An object that already has a projected CRS is returned untouched. Only
#' geographic (lon/lat) input is transformed, and the CRS chosen depends on the
#' extent of the data. It is **not** always UTM:
#'
#' \describe{
#'   \item{Local extents}{The UTM zone containing the data's centre
#'     (EPSG:326xx north of the equator, EPSG:327xx south). Distances and areas
#'     are close to true over a few degrees of longitude, which is the case
#'     this package is usually in. The centre is the centroid on the sphere,
#'     computed with s2 whatever [sf::sf_use_s2()] is set to, so data near a
#'     zone edge get the same zone in every session.}
#'   \item{Wide extents}{Once the data reach well beyond the roughly 3 degrees
#'     a UTM zone is designed for, a single zone can distort distances by
#'     several percent, and that error propagates straight into variogram
#'     ranges, block sizes, GWR bandwidths and GP length-scales. Which
#'     projection is actually best is then **measured, not assumed**: the zone,
#'     a Lambert azimuthal equal-area centred on the data and (where its
#'     standard parallels do not degenerate) an Albers conic are each scored by
#'     projecting representative points of the data (a non-POINT layer is
#'     reduced to one point per feature, plus the vertices of its outline when
#'     it has fewer than 40 features, so a single study-area polygon is scored
#'     too) and comparing planar with geodesic pairwise distances, and the one
#'     that distorts least is used.
#'     The choice, both error figures and this argument are **logged** (see the
#'     logging note under [spatialkit_quiet()]); they are not R warnings, so
#'     `tryCatch(warning = )` does not see them.}
#'   \item{Antimeridian}{Data straddling ±180° have a bounding box wider than a
#'     hemisphere. The wrap is detected from the coordinates (one very large
#'     gap in the sorted longitudes) and an equal-area projection centred on
#'     the true extent is used.}
#'   \item{Around a pole}{Data spanning more than 180 degrees of longitude
#'     with no such gap surround a pole. When every point also lies on one
#'     side of the equator (Antarctic stations, a pan-Arctic network), a
#'     Lambert azimuthal equal-area centred on that pole is used, provided it
#'     measures a smaller distance error than the global fallback. Web
#'     Mercator splits such a layer at +/-180 degrees and stretches it
#'     towards the pole: a ring of Antarctic stations measured a worst-case
#'     distance error near 20,000 percent in it, against about 2 percent in
#'     the polar projection. Only the coverage left over, spanning both
#'     hemispheres or a low-latitude belt the polar projection fits worse,
#'     falls back to EPSG:3857 (Equal Earth for `purpose = "area"`).}
#'   \item{Missing CRS}{With no `target_crs`, a bounding box that looks like
#'     lon/lat means EPSG:4326 is assumed (a real warning) and the rules above
#'     then apply; coordinates the heuristic declines are left exactly as they
#'     are. With `target_crs` supplied there is no source CRS to reproject
#'     from, so the same heuristic decides between two outcomes: lon/lat-looking
#'     coordinates are read as EPSG:4326 and reprojected to the target (a real
#'     warning), and anything else has the target **stamped on without
#'     reprojection**. That is a relabel, logged only, so verify the
#'     coordinates really are in that CRS. Set the CRS explicitly to suppress
#'     either.}
#' }
#'
#' `target_crs` overrides all of this. Pass it whenever you need a specific,
#' reproducible projection: comparing runs, matching an existing layer, or
#' fixing the units that [make_folds()]'s `block_size` will be interpreted in.
#'
#' @param x An sf or sfc object (other objects returned unchanged).
#' @param target_crs Optional target CRS (sf object, integer EPSG, or crs).
#'   Must resolve to a usable CRS via [sf::st_crs()]; an unusable value (one
#'   that resolves to `NA_crs_`) raises an error, so `x` is never left
#'   silently unprojected.
#' @param purpose Which property the projection is for.  `"distance"` (the
#'   default, and everything above): the candidate that distorts pairwise
#'   distances least, which is what ranges, block sizes, bandwidths and
#'   length-scales read off the coordinates.  `"area"`: densities or rates
#'   per cell are going to be computed, so the CRS must be equal-area.  For
#'   lon/lat input the choice is then made among equal-area projections only:
#'   a Lambert azimuthal centred on the data, or an Albers conic where its
#'   parallels do not degenerate, whichever distorts distances less.  A UTM
#'   zone (conformal, not equal-area) never enters that comparison, and global
#'   coverage gets Equal Earth in place of Web Mercator.  Already-projected
#'   input is still returned untouched, but its area distortion over the
#'   extent is measured (the spread of planar-to-geodesic area ratios over
#'   probe polygons) and logged as a warning when it exceeds 1 percent.
#'   Measured: a zone's own width edge to edge, 0.25 percent; the
#'   conterminous United States forced into one zone, 14 percent; Web
#'   Mercator over 2.5 degrees of latitude at 48N, 4 percent; an equal-area
#'   projection, a few tenths of a percent, which is the sphere the
#'   geodesic areas are computed on against the ellipsoid the projection
#'   uses.  [summarize_by_cell()] applies the same measurement before it
#'   computes a density.  Ignored when `target_crs` is given.
#' @return x, potentially with a new projected CRS.  When a projection was
#'   chosen here (lon/lat input, no \code{target_crs}) the result carries
#'   \code{attr(x, "crs_choice")}: a data.frame with one row per projection
#'   considered, holding \code{name}, \code{crs} (its definition), the measured
#'   worst-case \code{distance_error} (relative, over sampled pairs;
#'   \code{NA} where it could not be measured) and \code{chosen}.  The
#'   figure the log line quotes for the winner is therefore recoverable for
#'   every candidate, and is measured for the single candidate on the paths
#'   where no comparison runs (a UTM zone on a local extent, the equal-area
#'   projection chosen for a layer straddling the antimeridian).  It is
#'   \code{NULL} exactly when no local projection was chosen here: input
#'   that already carried a projected CRS, a \code{target_crs} you supplied,
#'   or an extent no local projection fits, which falls back to Web Mercator
#'   or Equal Earth.  CRS-less input
#'   additionally carries \code{attr(x, "crs_assumed")}: \code{"EPSG:4326"} when
#'   the lon/lat heuristic fired, \code{"none"} when it declined.  That
#'   attribute is also read on the way IN.  An object already carrying
#'   \code{"none"} is returned untouched, with the heuristic skipped, which is
#'   how a \code{predict()} method replays a fit's negative decision so that a
#'   subset of the training rows is not judged differently from the whole.
#' @family spatial data preparation
#' @examples
#' library(sf)
#' pts_ll <- st_as_sf(
#'   data.frame(lon = c(9.1, 9.2), lat = c(48.7, 48.8)),
#'   coords = c("lon", "lat"), crs = 4326
#' )
#' # A local extent gets the containing UTM zone; the zone's measured
#' # distance error over the extent travels with the result.
#' st_crs(ensure_projected(pts_ll))$epsg  # 32632
#' attr(ensure_projected(pts_ll), "crs_choice")
#'
#' # A continental extent is scored against the zone and may get an equal-area
#' # projection instead; the choice and both error figures are LOGGED, not
#' # warned -- see Details and ?spatialkit_quiet.
#' wide <- st_as_sf(
#'   data.frame(lon = c(-120, -70), lat = c(30, 48)),
#'   coords = c("lon", "lat"), crs = 4326
#' )
#' st_crs(ensure_projected(wide))$proj4string
#'
#' # target_crs overrides the choice entirely.
#' st_crs(ensure_projected(pts_ll, target_crs = 3035))$epsg  # 3035
#'
#' # For densities per cell the CRS has to be equal-area.  Two candidates are
#' # scored, an Albers conic and a Lambert azimuthal, both centred on the data,
#' # and whichever distorts distance less over the extent is used: here Albers.
#' st_crs(ensure_projected(pts_ll, purpose = "area"))$proj4string
#' @export
ensure_projected <- function(x, target_crs = NULL, purpose = c("distance", "area")) {
  purpose <- match.arg(purpose)
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
      # No source CRS: we cannot reproject, only assume -- and the assumption
      # must be the SAME one the no-target branch below makes, or the same
      # rows land in different places at fit time (no target) and at predict
      # time (target = the fit's CRS).  Measured: predict(fit, newdata =
      # training rows) differed from fitted(fit) by up to one response SD.
      ll <- .looks_like_lonlat(x)
      if (isTRUE(ll$lonlat) && !isTRUE(sf::st_is_longlat(tcrs))) {
        .warn_and_log(
          "ensure_projected(): input has no CRS; its coordinates look like lon/lat (xmin=%.2f, xmax=%.2f, ymin=%.2f, ymax=%.2f), so they are taken as EPSG:4326 and reprojected to the target CRS ('%s'). Set the CRS explicitly to suppress this.",
          ll$bb["xmin"], ll$bb["xmax"], ll$bb["ymin"], ll$bb["ymax"],
          tcrs$input %||% "unknown"
        )
        sf::st_crs(x) <- sf::st_crs(4326)
        x <- sf::st_transform(x, tcrs)
        attr(x, "crs_assumed") <- "EPSG:4326"
        return(x)
      }
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
    # A predict() method has already replayed the fit's decision: the training
    # data had no CRS and were NOT taken for lon/lat, so these coordinates
    # belong in the same untouched space.  Re-running the heuristic here would
    # judge a subset differently from the whole (see .replay_crs_assumption()).
    if (identical(attr(x, "crs_assumed"), "none")) return(x)

    ll <- .looks_like_lonlat(x)
    if (isTRUE(ll$lonlat)) {
      .warn_and_log(
        "ensure_projected(): CRS is missing; assuming EPSG:4326 because bbox looks like lon/lat (xmin=%.2f, xmax=%.2f, ymin=%.2f, ymax=%.2f). Set CRS explicitly to suppress this warning.",
        ll$bb["xmin"], ll$bb["xmax"], ll$bb["ymin"], ll$bb["ymax"]
      )
      sf::st_crs(x) <- sf::st_crs(4326)
      tr <- .pick_local_projected_crs(x, purpose = purpose)
      x <- sf::st_transform(x, tr$crs)
      attr(x, "crs_choice") <- tr$candidates
      # Recorded so a predict() on CRS-less newdata can replay the same
      # interpretation (see .replay_crs_assumption()).
      attr(x, "crs_assumed") <- "EPSG:4326"
    } else if (!is.null(ll$bb) &&
               ll$bb["xmin"] >= -180 && ll$bb["xmax"] <= 180 &&
               ll$bb["ymin"] >= -90  && ll$bb["ymax"] <= 90) {
      .log_warn(
        "ensure_projected(): CRS is missing and coordinates fall within the lon/lat envelope, but they lack the decimal precision or extent typical of geographic coordinates (xmin=%.2f, xmax=%.2f, ymin=%.2f, ymax=%.2f). Not assuming EPSG:4326. Set CRS explicitly with sf::st_crs(x) <- sf::st_crs(...).",
        ll$bb["xmin"], ll$bb["xmax"], ll$bb["ymin"], ll$bb["ymax"]
      )
    }
    # Record the NEGATIVE decision as well, so a predict() on CRS-less newdata
    # replays "these are planar, leave them alone" instead of re-deciding from
    # newdata's own bounding box.
    if (!isTRUE(attr(x, "crs_assumed") == "EPSG:4326"))
      attr(x, "crs_assumed") <- "none"
    return(x)
  }

  if (.is_longlat(x)) {
    # The projections considered and their measured distance errors ride on
    # the result as `crs_choice`: the number the warning above quotes for
    # the winner and the zone is kept for every candidate, and on the
    # ordinary path -- a zone kept without a comparison -- for the zone.
    tr <- .pick_local_projected_crs(x, purpose = purpose)
    x <- sf::st_transform(x, tr$crs)
    attr(x, "crs_choice") <- tr$candidates
  } else if (purpose == "area") {
    # Projected input is never reprojected, but when areas are what matter
    # the caller should know whether this CRS can deliver them.  Measured, as
    # the distance choice is: the spread of planar-to-geodesic area ratios
    # over probe polygons on the extent.
    err <- .crs_area_error(x)
    if (is.finite(err) && err > .area_error_tol)
      .log_warn(
        paste0("ensure_projected(purpose = \"area\"): the input's CRS (%s) ",
               "distorts areas across this extent by up to %.1f%% (planar ",
               "against geodesic, spread over probe polygons); densities and ",
               "rates computed in it are not comparable between cells. It is ",
               "returned as it is. Pass an equal-area target_crs, or start ",
               "from lon/lat input and let purpose = \"area\" choose one."),
        .fold_crs_label(x), 100 * err)
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
#' When one input carries no CRS, the same lon/lat heuristic
#' [ensure_projected()] uses decides what happens, so both entry points place
#' identical data in the same place: coordinates that look like degrees are
#' taken as EPSG:4326 and **reprojected** to the other object's CRS (or
#' `target_crs`); coordinates that do not are *stamped* with
#' [sf::st_set_crs()], which relabels without moving them.  Either way the
#' assumption is announced with a warning.
#'
#' @param a,b Objects of class sf or sfc.
#' @param prefer Which object's CRS to keep ("a" or "b").
#' @param target_crs Optional target CRS to apply to both: anything
#'   [sf::st_crs()] accepts, including an sf or sfc object, whose CRS is used.
#' @param on_transform_error What to do when st_transform() fails:
#'   \code{"stop"} (default) raises an error immediately;
#'   \code{"set_crs"} falls back to st_set_crs() (UNSAFE: coordinates are
#'   NOT reprojected, only the CRS label is overwritten). The \code{"set_crs"}
#'   option exists only for rare edge cases where you are certain the
#'   coordinates already match the target CRS definition.
#' @return A named list with components a and b.
#' @family spatial data preparation
#' @examples
#' library(sf)
#' a <- st_as_sf(data.frame(x = c(500000, 500100), y = c(4000000, 4000100)),
#'               coords = c("x", "y"), crs = 32632)
#' b <- st_transform(a, 4326)                 # same points, lon/lat
#' h <- harmonize_crs(a, b)                    # b is brought into a's CRS
#' c(a = st_crs(h$a)$epsg, b = st_crs(h$b)$epsg)
#' st_crs(h$a) == st_crs(h$b)
#' @export
harmonize_crs <- function(a, b, prefer = c("a", "b"), target_crs = NULL,
                          on_transform_error = c("stop", "set_crs")) {
  if (!inherits(a, c("sf", "sfc"))) stop("harmonize_crs(): `a` must be sf or sfc.")
  if (!inherits(b, c("sf", "sfc"))) stop("harmonize_crs(): `b` must be sf or sfc.")
  prefer <- match.arg(prefer)
  on_transform_error <- match.arg(on_transform_error)
  # A layer as the target means its CRS, as it does for ensure_projected().
  # Passed through as it was, st_transform() read a multi-row sf as a list of
  # candidate CRSs ("the condition has length > 1") and refused a one-row one.
  # Only sf/sfc are converted here: st_crs() on a string it cannot parse
  # throws before st_transform() runs, which would bypass on_transform_error.
  if (inherits(target_crs, c("sf", "sfc"))) target_crs <- sf::st_crs(target_crs)

  crs_a <- sf::st_crs(a)
  crs_b <- sf::st_crs(b)

  .safe_transform <- function(x, to) {
    tryCatch(
      sf::st_transform(x, to),
      error = function(e) {
        if (identical(on_transform_error, "set_crs")) {
          .warn_and_log(
            "harmonize_crs(): st_transform() failed (%s); falling back to st_set_crs(). WARNING: coordinates are NOT reprojected -- downstream distances, joins, and areas may be wrong. Set on_transform_error='stop' (the default) to surface this error instead.",
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

  # st_set_crs() only stamps a label: the coordinates are NOT moved.  Stamping
  # is therefore the LAST resort: ensure_projected() runs the lon/lat
  # heuristic first and reprojects when the coordinates look like degrees, and
  # this function's Details said it matched that behaviour while in fact
  # always stamping -- placing identical input 5,400 km apart depending on
  # which entry point a caller used.  .resolve_crsless() closes the gap: it
  # reprojects lon/lat-looking coordinates and only stamps otherwise, warning
  # either way.
  .resolve_crsless <- function(x, which_obj, to) {
    lbl <- tryCatch({
      inp <- sf::st_crs(to)$input
      if (is.null(inp) || length(inp) != 1L || is.na(inp) ||
          !nzchar(as.character(inp))) "unknown" else as.character(inp)
    }, error = function(e) "unknown")
    ll <- .looks_like_lonlat(x)
    if (isTRUE(ll$lonlat)) {
      .warn_and_log(
        paste0("harmonize_crs(): `%s` has no CRS; its coordinates look like ",
               "lon/lat (xmin=%.2f, xmax=%.2f, ymin=%.2f, ymax=%.2f), so they ",
               "are taken as EPSG:4326 and reprojected to '%s'. Set the CRS ",
               "explicitly with sf::st_crs() to suppress this."),
        which_obj, ll$bb[["xmin"]], ll$bb[["xmax"]], ll$bb[["ymin"]],
        ll$bb[["ymax"]], lbl)
      return(.safe_transform(sf::st_set_crs(x, 4326), to))
    }
    .warn_and_log(
      "harmonize_crs(): `%s` has no CRS and its coordinates do not look like lon/lat; stamping '%s' onto it WITHOUT reprojection. Verify the coordinates are already expressed in that CRS, or set it explicitly with sf::st_crs().",
      which_obj, lbl
    )
    sf::st_set_crs(x, to)
  }

  if (!is.null(target_crs)) {
    if (is.na(crs_a)) {
      a <- .resolve_crsless(a, "a", target_crs)
    } else {
      a <- .safe_transform(a, target_crs)
    }
    if (is.na(crs_b)) {
      b <- .resolve_crsless(b, "b", target_crs)
    } else {
      b <- .safe_transform(b, target_crs)
    }
    return(list(a = a, b = b))
  }

  if (is.na(crs_a) && is.na(crs_b)) return(list(a = a, b = b))

  if (is.na(crs_a) && !is.na(crs_b)) {
    a <- .resolve_crsless(a, "a", crs_b)
    return(list(a = a, b = b))
  }
  if (!is.na(crs_a) && is.na(crs_b)) {
    b <- .resolve_crsless(b, "b", crs_a)
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
    sf::st_point(unname(c((bb[["xmin"]] + bb[["xmax"]]) / 2,
                          (bb[["ymin"]] + bb[["ymax"]]) / 2)))
  })
  sf::st_sfc(pts, crs = sf::st_crs(x))
}


#' Drop the EMPTY parts of multi-part geometries
#'
#' Every call to [sf::st_point_on_surface()] in the package goes through this
#' first.  GEOS (3.12.1, which sf 1.0.x links) SEGFAULTS computing the
#' interior point of a non-empty geometry that holds an EMPTY line: a
#' MULTILINESTRING with an empty part beside a real one, or a
#' GEOMETRYCOLLECTION with an empty LINESTRING among its members.  The R
#' session is lost, not merely the call.  An empty POLYGON member does not
#' crash but is worse in its way: GEOS takes the interior point from the
#' highest dimension present, finds that dimension empty, and returns
#' POINT EMPTY for a geometry that has a line in it.
#'
#' An empty part adds no points to the geometry, so dropping it changes
#' nothing but those two failures.  A feature left with no parts is EMPTY as
#' a whole, which GEOS handles (it gives POINT EMPTY).  Features are never
#' removed, so the result stays aligned row for row with the input, and a
#' feature with no empty part is returned exactly as it was.
#'
#' @param x An sf or sfc object.
#' @return \code{x}, with the empty parts removed from MULTILINESTRING,
#'   MULTIPOLYGON and GEOMETRYCOLLECTION features (recursively for a
#'   collection's members).
#' @keywords internal
#' @noRd
.drop_empty_parts <- function(x) {
  if (inherits(x, "sf")) {
    g  <- sf::st_geometry(x)
    g2 <- .drop_empty_parts(g)
    return(if (identical(g2, g)) x else sf::st_set_geometry(x, g2))
  }
  multi <- c("MULTILINESTRING", "MULTIPOLYGON", "GEOMETRYCOLLECTION")
  idx   <- which(as.character(sf::st_geometry_type(x, by_geometry = TRUE)) %in% multi)
  if (!length(idx)) return(x)

  # Emptiness read off the structure, without a GEOS call per part: a matrix
  # with no rows (a LINESTRING or a ring; sf refuses NA in one), a POINT with
  # no coordinates, or a list with no non-empty element (a POLYGON's rings,
  # the parts of a multi-geometry, a collection's members).
  is_empty <- function(s) {
    if (is.list(s)) return(all(vapply(s, is_empty, logical(1))))
    length(s) == 0L || (!is.matrix(s) && all(is.na(s)))
  }
  has_empty <- function(s) {
    if (!inherits(s, multi)) return(FALSE)
    for (k in unclass(s)) if (is_empty(k) || has_empty(k)) return(TRUE)
    FALSE
  }
  strip <- function(s) {
    if (!inherits(s, multi)) return(s)
    kids <- unclass(s)
    if (inherits(s, "GEOMETRYCOLLECTION")) kids <- lapply(kids, strip)
    structure(kids[!vapply(kids, is_empty, logical(1))], class = class(s))
  }

  # Only features that hold an empty part are rebuilt; the rest, which is
  # nearly always all of them, are left exactly as they were.
  for (i in idx[vapply(unclass(x)[idx], has_empty, logical(1))])
    x[[i]] <- strip(x[[i]])
  x
}

# -----------------------------------------------------------------------------
# Point Coercion
# -----------------------------------------------------------------------------

#' Coerce arbitrary geometries to representative points
#'
#' Converts the geometry column of an sf object to POINTs using one of several
#' strategies.
#'
#' The result has one row per row of `x`, in the same order.  An EMPTY
#' geometry of any type, lines included, becomes an EMPTY POINT in its own
#' row; [prep_model_data()] and [make_folds()] then drop such rows, as they
#' drop any other empty geometry.  Empty lines are never handed to
#' [sf::st_line_sample()]: it yields no midpoint for them, which would
#' misalign the result, and with sf 1.0.x an empty MULTILINESTRING (or an
#' empty part of one) crashed the R session.  An empty part inside a
#' non-empty feature is ignored, so the feature gets the point its other
#' parts give; GEOS's interior point, used by `"point_on_surface"` and by
#' the temporary projection's choice of CRS, segfaulted on an empty line
#' part too.
#'
#' @param x An sf object.
#' @param mode One of "auto", "centroid", "point_on_surface", "surface",
#'   "line_midpoint", "bbox_center".
#' @param tmp_project Logical; temporarily project for line-based midpoints.
#'   When \code{x} has no CRS and its coordinates fall inside the lon/lat
#'   envelope, that temporary projection interprets them as EPSG:4326 (with a
#'   warning) and the midpoints returned are geodesic ones brought back to the
#'   input's numbers, not planar midpoints.  Set the CRS, or pass
#'   \code{tmp_project = FALSE}, for planar data.
#' @return An sf object with geometry coerced to POINTs, row for row with
#'   `x`; an empty input geometry gives an empty POINT.
#' @family spatial data preparation
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

  # -- direct ops ---
  # Not spherical-safe, whatever this heading used to say.  On lon/lat input
  # st_centroid() is spherical only while sf_use_s2() is TRUE: with it FALSE
  # it is planar in degrees (a box -120..-60 x 50..75 got a centre 155 km
  # from the s2 one) and its warning is suppressed here.  st_point_on_surface()
  # is GEOS, planar in degrees under either setting.
  if (mode == "centroid") {
    return(sf::st_set_geometry(x, suppressWarnings(sf::st_centroid(g))))
  }
  if (mode == "point_on_surface") {
    # Empty parts dropped first: GEOS segfaults on an empty line inside a
    # non-empty feature (see .drop_empty_parts()).
    return(sf::st_set_geometry(x, sf::st_point_on_surface(.drop_empty_parts(g))))
  }

  is_ll <- .is_longlat(x)

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

  # Midpoints of the LINESTRING rows `idx_ls`, one per row, batched: project
  # once, sample all, back-transform once.  An EMPTY line has no midpoint:
  # st_line_sample() yields an empty MULTIPOINT that st_cast(, "POINT")
  # silently drops (and in sf 1.0.x the call can segfault outright), so the
  # samples would no longer align 1:1 with the rows they are scattered back
  # into.  Empty rows never reach the sampler.  They get an EMPTY POINT, as an
  # empty polygon, point or collection already does, so the rows stay aligned
  # and prep_model_data() and make_folds() drop them like any empty geometry.
  # This used to be an error, which made a line layer with one null geometry
  # the only kind of layer those "drop empty rows" paths could not clean.
  .line_midpoints <- function(idx_ls) {
    res  <- rep(list(sf::st_point()), length(idx_ls))
    full <- which(!sf::st_is_empty(g[idx_ls]))
    if (!length(full)) return(res)
    g_ls_sf   <- sf::st_sf(geometry = g[idx_ls[full]])
    g_ls_proj <- if (tmp_project) ensure_projected(g_ls_sf) else g_ls_sf
    midps     <- sf::st_line_sample(sf::st_geometry(g_ls_proj), sample = 0.5)
    midps     <- sf::st_cast(midps, "POINT")
    midps     <- .back_to_input_crs(midps, g_ls_proj, crs)
    .check_midpoint_alignment(midps, full)
    res[full] <- as.list(midps)
    res
  }

  if (mode == "line_midpoint") {
    gtypes <- as.character(sf::st_geometry_type(g, by_geometry = TRUE))
    if (any(gtypes %in% c("MULTILINESTRING", "GEOMETRYCOLLECTION"))) {
      stop("coerce_to_points(): method \"line_midpoint\" only supports LINESTRING; cast MULTILINESTRING first.", call. = FALSE)
    }
    idx_ls <- which(gtypes == "LINESTRING")
    idx_other <- which(gtypes != "LINESTRING")
    out <- vector("list", length(g))

    if (length(idx_ls)) {
      if (is_ll && !tmp_project) {
        ctr <- suppressWarnings(sf::st_centroid(g[idx_ls]))
        out[idx_ls] <- as.list(ctr)
      } else {
        out[idx_ls] <- .line_midpoints(idx_ls)
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
    pos <- sf::st_point_on_surface(.drop_empty_parts(g[idx_poly]))
    out[idx_poly] <- as.list(pos)
  }

  # --- LINESTRING: midpoint via line_sample (batched) ---
  idx_ls <- which(gtypes == "LINESTRING")
  if (length(idx_ls)) {
    if (is_ll && !tmp_project) {
      ctr <- suppressWarnings(sf::st_centroid(g[idx_ls]))
      out[idx_ls] <- as.list(ctr)
    } else {
      out[idx_ls] <- .line_midpoints(idx_ls)
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
        # st_cast() turns an EMPTY MULTILINESTRING into ONE empty LINESTRING,
        # not zero parts, and a MULTILINESTRING can also carry an empty part
        # beside real ones.  Either reached st_line_sample() below, which
        # segfaults on an empty line in sf 1.0.x and took the R session with
        # it (the usual source: a null geometry in a line layer, which
        # GeoPackage and st_read()'s promote_to_multi return as
        # MULTILINESTRING EMPTY).  Sample only parts that have a midpoint; a
        # feature with none gets an EMPTY POINT, as an empty LINESTRING does.
        parts <- suppressWarnings(sf::st_cast(proj_geom[j], "LINESTRING"))
        parts <- parts[!sf::st_is_empty(parts)]
        if (length(parts) == 0L) {
          out[[idx_mls[j]]] <- sf::st_point()
        } else {
          lens <- as.numeric(sf::st_length(parts))
          k    <- if (length(lens)) which.max(lens) else 1L
          mp   <- sf::st_line_sample(parts[k], sample = 0.5)
          mp   <- sf::st_cast(mp, "POINT")
          mp   <- .back_to_input_crs(mp, g_mls_proj, crs)
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
