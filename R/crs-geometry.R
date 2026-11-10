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
#' \code{st_transform()} of it works.  Now s2 is always used, and a layer s2
#' rejects is repaired and tried again.  If s2 rejects the repaired layer
#' too, the features it refuses are repaired one by one (\code{.s2_mend()}:
#' without their repeated vertices, then with crossing edges split), and the
#' centre is the centroid of the union of the features s2 then accepts; a
#' log line names the features s2 still refuses, which are left out.  Only
#' when s2 accepts none of the features is the centre the normalised mean of
#' the vertices' unit vectors (for a point layer exactly what s2 returns),
#' which needs no valid geometry; each ring's closing vertex counts once
#' there, so the mean does not depend on where a ring starts, and a log
#' line says it was used.  A row with a coordinate that is not finite or lies
#' beyond \code{.POOL_XY_MAX} is read as an empty row
#' (\code{.empty_unlocated()}).
#'
#' @param x_ll An sf/sfc object in a geographic CRS.
#' @return A 1 x 2 matrix (\code{X}, \code{Y}) in degrees, possibly
#'   non-finite when no centre exists (the caller falls back then).
#' @keywords internal
#' @noRd
.lonlat_centre <- function(x_ll) {
  # A row without a location places nothing: it is read as an empty row.
  g <- sf::st_geometry(.empty_unlocated(x_ll))
  centre_of <- function(geom) {
    ctr <- sf::st_coordinates(sf::st_centroid(sf::st_union(geom)))
    if (!is.numeric(ctr) || length(ctr) < 2L) stop("no centroid")
    ctr[1L, 1:2, drop = FALSE]
  }
  ctr <- tryCatch(.with_s2(centre_of(g)), error = function(e) NULL)
  if (!is.null(ctr)) return(ctr)
  ctr <- tryCatch(.with_s2(centre_of(.safe_make_valid(g))), error = function(e) NULL)
  if (!is.null(ctr)) return(ctr)

  # One feature s2 refuses even repaired (a bowtie: s2's repair leaves it as
  # it is) sent the whole layer to the mean of its vertices, in silence: on
  # 100 one-degree squares, one of them drawn with 4,000 vertices, the
  # centre moved 4 degrees to that square, and ensure_projected() chose an
  # Albers CRS centred there in place of UTM zone 14.  The features s2
  # refuses are repaired one by one, and those it refuses even then are
  # left out of the centre, with a log line.
  rows_text <- function(r) {
    k <- length(r)
    if (k == 1L) sprintf("row %d", r)
    else if (k <= 6L) sprintf("rows %s and %d", paste(r[-k], collapse = ", "), r[k])
    else sprintf("rows %s and %d more", paste(r[1:6], collapse = ", "), k - 6L)
  }
  full <- which(!sf::st_is_empty(g))
  m    <- tryCatch(.with_s2(.s2_mend(g, split = TRUE)), error = function(e) NULL)
  keep <- if (is.null(m)) integer(0) else setdiff(full, m$refused)
  if (length(keep)) {
    ctr <- tryCatch(.with_s2(centre_of(m$g[keep])), error = function(e) NULL)
    # Counted over the features that are not empty: an empty row holds no
    # location, and counted in, the figures did not add up.
    if (!is.null(ctr) && length(m$refused))
      .log_info(paste0(".pick_local_projected_crs(): s2 refuses %d of the %s ",
                       "in longitude and latitude even repaired (%s); ",
                       "the centre that places the projection is taken over the ",
                       "other %d."),
                length(m$refused),
                if (length(full) < length(g))
                  sprintf("%d features that are not empty", length(full))
                else sprintf("%d features", length(g)),
                rows_text(m$refused), length(keep))
    if (!is.null(ctr)) return(ctr)
  }

  # The mean of the vertices' unit vectors, read feature by feature on a
  # layer of mixed types or of collections, which sf::st_coordinates()
  # refuses.  Each ring's closing vertex counts once: counted twice, the
  # mean moved with the vertex a ring starts at (a bowtie listed from each of
  # its four corners gave four centres and four Albers CRSs).
  use <- if (length(keep)) keep else full
  xy  <- tryCatch(.open_ring_xy(g[use]), error = function(e) matrix(numeric(0), 0L, 2L))
  lam <- xy[, 1L] * pi / 180; phi <- xy[, 2L] * pi / 180
  v   <- c(sum(cos(phi) * cos(lam)), sum(cos(phi) * sin(lam)), sum(sin(phi)))
  ctr <- matrix(if (nrow(xy) && sqrt(sum(v^2)) > 1e-12)
                  c(atan2(v[2L], v[1L]), atan2(v[3L], sqrt(v[1L]^2 + v[2L]^2))) * 180 / pi
                else c(NA_real_, NA_real_),
                1L, 2L, dimnames = list(NULL, c("X", "Y")))
  if (all(is.finite(ctr)))
    .log_info(paste0(".pick_local_projected_crs(): %s, so the centre that places ",
                     "the projection is the mean of the vertices of %s."),
              if (is.null(m)) "s2 stopped while the features it refuses were repaired"
              else if (length(keep))
                sprintf("s2 takes no union of the %d features it accepts", length(keep))
              else "s2 refuses every feature in longitude and latitude, even repaired",
              if (length(keep)) "those features" else "the layer")
  ctr
}


#' Whether an Albers conic with these standard parallels is a candidate
#'
#' The cone constant must be at least 0.05 in absolute value (parallels that
#' straddle the equator degenerate, and PROJ refuses lat_1 = -lat_2), and
#' both parallels must lie between -90 and 90: the fallback parallels
#' \code{lat - 5} and \code{lat + 5} pass a pole when the centre is within 5
#' degrees of it, and PROJ refuses them.  A conic at a pole is the polar
#' azimuthal case, which the Lambert azimuthal candidate centred there
#' covers.
#'
#' @param lat1,lat2 Standard parallels in degrees.
#' @return \code{TRUE} or \code{FALSE}.
#' @keywords internal
#' @noRd
.albers_parallels_ok <- function(lat1, lat2) {
  cone_n <- (sin(lat1 * pi / 180) + sin(lat2 * pi / 180)) / 2
  isTRUE(is.finite(cone_n) && abs(cone_n) >= 0.05 &&
           abs(lat1) <= 90 && abs(lat2) <= 90)
}

#' The largest offset of a longitude box from a meridian
#'
#' The largest distance in longitude, on the circle, from \code{cm} to any
#' longitude between \code{lon_min} and \code{lon_max}: 180 when the meridian
#' opposite \code{cm} lies in the box, otherwise the larger of the two ends'
#' offsets, each wrapped to 0 to 180.  The zone of the CRS chooser comes from
#' the centre on the sphere, which for a layer across a pole, or one whose
#' longitudes are given outside -180 to 180, can lie far outside the box, and
#' the unwrapped difference of the box's ends then came out above 180.
#'
#' @param lon_min,lon_max The box's ends in degrees, \code{lon_min <= lon_max}.
#' @param cm The meridian in degrees.
#' @return A number from 0 to 180, or a non-finite number when the box or
#'   \code{cm} is not finite.
#' @keywords internal
#' @noRd
.lon_offset_from <- function(lon_min, lon_max, cm) {
  wrap <- function(x) { d <- abs(x - cm) %% 360; pmin(d, 360 - d) }
  off <- max(wrap(lon_min), wrap(lon_max))
  if (!is.finite(off)) return(off)
  # The meridian opposite cm lies in the box when, measured eastwards from
  # lon_min, it comes no later than lon_max.
  if ((cm + 180 - lon_min) %% 360 <= lon_max - lon_min) 180 else off
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
#' Rows with a coordinate that is not finite or lies beyond
#' \code{.POOL_XY_MAX} take no part in the choice: it is made on the layer
#' with those rows emptied (\code{.empty_unlocated()}).  So do rows with a
#' vertex at a latitude beyond 90 (\code{.beyond_pole_features()}), read on
#' the copy in WGS84.
#'
#' On the comparison path and the \code{purpose = "area"} path a candidate
#' is used only if it places every located vertex of the layer
#' (\code{.places_vertices()}; a vertex at a latitude beyond 90 takes no
#' part, \code{.located_vertices()}).  One that does not gets
#' \code{distance_error} \code{NA}, is not chosen, and the log line names
#' it.  When no candidate places every vertex, the global fallback is used.
#' The Albers candidate is built only where \code{.albers_parallels_ok()}
#' holds for its standard parallels.
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
  # Rows without a location (a coordinate that is not finite or lies beyond
  # .POOL_XY_MAX) take no part in the choice: the box, the longitudes of the
  # wrap test, the centre and the scores are read as with those rows empty.
  # Read, a longitude of 1e200 made a box "across the antimeridian" and a
  # lon_0 PROJ refuses, and a latitude of 1e200 moved the standard parallels
  # of every other row's projection.
  x <- .empty_unlocated(x)
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
  # A row with a vertex at a latitude beyond 90 has no place on the Earth,
  # and no CRS places it: it takes no part in the choice either, the whole
  # row, as a row with a coordinate beyond 1e150.  Read, one station at
  # (-96, 95) among 39 over the conterminous United States moved the centre
  # and the Albers parallels (41.2 and 84.2 in place of 32.9 and 43.0), and
  # a swapped row at (40, -96) put the box 139 degrees from the zone's
  # central meridian.  Tested in degrees on the WGS84 copy; the layer
  # returned to the caller is not changed.
  x_ll <- .empty_unlocated(x_ll, which(.beyond_pole_features(sf::st_geometry(x_ll))))
  # No row has a location: nothing places a projection.  Read on, the centre
  # came out as NaN and the log line below quoted "lon = NaN, lat = NaN".
  # Empty parts are dropped first: sf::st_is_empty() calls a feature empty
  # when its first part is empty (.drop_empty_parts()).
  if (all(sf::st_is_empty(.drop_empty_parts(sf::st_geometry(x_ll))))) {
    .log_warn(paste0(
      ".pick_local_projected_crs(): no row has a location (each geometry is ",
      "empty, or has a coordinate that is not finite or lies beyond 1e150, ",
      "or a latitude beyond 90); falling back to %s."),
      if (purpose == "area") "Equal Earth" else "EPSG:3857")
    return(list(crs = global_crs(), candidates = NULL))
  }

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
    # Read feature by feature on a layer of mixed types or of collections,
    # which sf::st_coordinates() refuses: the refusal stopped
    # ensure_projected() on a box across 180 listed as a POLYGON beside a
    # MULTIPOLYGON.
    lons  <- .vertex_coords(sf::st_geometry(x_ll))[, 1L]
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
  # The largest offset of the box from that meridian, on the circle.
  lon_off   <- .lon_offset_from(lon_min, lon_max, cand_cm)

  utm_epsg <- if (!is.na(lat) && lat < 0) 32700 + cand_zone else 32600 + cand_zone

  # A candidate is used only if it places every located vertex of the layer.
  # PROJ gives no coordinates for a point a UTM zone cannot place (about 90
  # degrees or more from its central meridian on the equator), and the
  # score left such points out: LINESTRING (0 0, 180 0) kept UTM zone 46,
  # scored without both of its ends, and ensure_projected() returned the
  # line empty, without a word.  The vertices are read once, and only where
  # candidates are compared.
  vtx <- NULL
  places <- function(cd) {
    if (is.null(vtx)) vtx <<- .located_vertices(x_ll)
    .places_vertices(vtx, sf::st_crs(x_ll), cd$crs)
  }
  # Named in the log line: the candidates that cannot place the layer.
  not_placed <- function(cands, placed, verb = "cannot")
    if (all(placed)) "" else
      sprintf(" %s %s place every location of the layer, so %s not used.",
              paste(vapply(cands[!placed], `[[`, character(1), "name"),
                    collapse = " and "),
              verb, if (sum(!placed) == 1L) "it is" else "they are")
  # No candidate places every vertex: the global fallback, as for global
  # coverage.
  none_placed <- function(cands) {
    .log_warn(
      paste0(".pick_local_projected_crs(): no candidate projection can place ",
             "every location of the layer (%s cannot). Falling back to %s. ",
             "Pass target_crs to ensure_projected() to choose a projection ",
             "suited to your extent."),
      paste(vapply(cands, `[[`, character(1), "name"), collapse = " and "),
      if (purpose == "area") "Equal Earth (equal-area)" else "EPSG:3857")
    list(crs = global_crs(), candidates = NULL)
  }

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
    cands <- list(list(
      name = sprintf("Lambert azimuthal equal-area centred on (%.1f, %.1f)", lon, lat),
      crs = sf::st_crs(sprintf(
        "+proj=laea +lat_0=%f +lon_0=%f +datum=WGS84 +units=m +no_defs", lat, lon))))
    if (.albers_parallels_ok(lat1, lat2)) {
      cands <- c(cands, list(list(
        name = sprintf("Albers equal-area (lat_1=%.1f, lat_2=%.1f, lon_0=%.1f)",
                       lat1, lat2, lon),
        crs = sf::st_crs(sprintf(
          paste0("+proj=aea +lat_1=%f +lat_2=%f +lat_0=%f +lon_0=%f ",
                 "+datum=WGS84 +units=m +no_defs"),
          lat1, lat2, lat, lon)))))
    }
    placed <- vapply(cands, places, logical(1))
    if (!any(placed)) return(none_placed(cands))
    err  <- vapply(seq_along(cands), function(i)
      if (placed[i]) .crs_distance_error(x_ll, cands[[i]]$crs) else NA_real_, numeric(1))
    best <- which.min(ifelse(is.finite(err), err, Inf))
    if (!length(best) || !is.finite(err[best])) best <- which(placed)[1L]
    .log_info(
      paste0(".pick_local_projected_crs(): purpose = \"area\": using %s ",
             "(equal-area; worst-case distance error %s over sampled pairs) ",
             "rather than UTM zone %d.%s"),
      cands[[best]]$name,
      if (is.finite(err[best])) sprintf("%.2f%%", 100 * err[best]) else "not measurable",
      cand_zone, not_placed(cands, placed)
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

    cands <- list(
      list(name = sprintf("UTM zone %d", cand_zone), crs = sf::st_crs(utm_epsg)),
      list(name = sprintf("Lambert azimuthal equal-area centred on (%.1f, %.1f)",
                          lon, lat),
           crs = sf::st_crs(sprintf(
             "+proj=laea +lat_0=%f +lon_0=%f +datum=WGS84 +units=m +no_defs",
             lat, lon)))
    )
    # A conic whose parallels straddle the equator is not merely poor, it is
    # degenerate: |lat_1 + lat_2| = 0 is rejected by PROJ outright.  Fallback
    # parallels that pass a pole (a centre within 5 degrees of it) are
    # rejected too.  Leave such a conic out of the comparison rather than
    # letting st_crs() raise.
    if (.albers_parallels_ok(lat1, lat2)) {
      cands <- c(cands, list(list(
        name = sprintf("Albers equal-area (lat_1=%.1f, lat_2=%.1f, lon_0=%.1f)",
                       lat1, lat2, lon),
        crs = sf::st_crs(sprintf(
          paste0("+proj=aea +lat_1=%f +lat_2=%f +lat_0=%f +lon_0=%f ",
                 "+datum=WGS84 +units=m +no_defs"),
          lat1, lat2, lat, lon)))))
    }

    # A candidate that cannot place every located vertex is scored NA and
    # not chosen; the zone is no longer kept by default when it is such a
    # candidate.
    placed <- vapply(cands, places, logical(1))
    if (!any(placed)) return(none_placed(cands))
    err <- vapply(seq_along(cands), function(i)
      if (placed[i]) .crs_distance_error(x_ll, cands[[i]]$crs) else NA_real_,
      numeric(1))
    best <- which.min(ifelse(is.finite(err), err, Inf))
    # Nothing measured: the zone when it places the layer, else the first
    # candidate that does.
    if (!length(best) || !is.finite(err[best])) best <- which(placed)[1L]

    if (identical(cands[[best]]$name, cands[[1L]]$name)) {
      # The figures are only quotable when the comparison actually ran; a
      # candidate that could not be scored comes back NA, and a message
      # reporting "NA% vs NA%" asserts a comparison that never happened.
      measured <- all(is.finite(err[placed]))
      .log_warn(
        paste0(".pick_local_projected_crs(): extent reaches %.1f deg from the ",
               "central meridian of UTM zone %d (%.1f deg of longitude in ",
               "total); the zone is kept%s.%s Pass target_crs to ",
               "ensure_projected() to override."),
        lon_off, cand_zone, span_lon,
        if (!any(placed[-1L]))
          " because no other candidate can place every location of the layer"
        else if (measured)
          sprintf(", because on this extent every equal-area candidate distorts distances more than it does (%.2f%% vs %.2f%% worst-case over sampled pairs)",
                  100 * min(err[-1L][placed[-1L]]), 100 * err[1L])
        else " because the candidate projections could not be scored on this layer",
        not_placed(cands, placed)
      )
      return(scored(sf::st_crs(utm_epsg), vapply(cands, `[[`, character(1), "name"),
                    lapply(cands, `[[`, "crs"), err, 1L))
    }

    .log_warn(
      paste0(".pick_local_projected_crs(): extent reaches %.1f deg from the ",
             "central meridian of UTM zone %d (%.1f deg of longitude in ",
             "total). Using %s instead: %s%s%s%s Pass target_crs to ",
             "ensure_projected() to override."),
      lon_off, cand_zone, span_lon, cands[[best]]$name,
      if (is.finite(err[best]))
        sprintf("measured worst-case distance error %.2f%%", 100 * err[best])
      else "its worst-case distance error could not be measured",
      if (placed[1L]) sprintf(" against the zone's %.2f%%", 100 * err[1L]) else "",
      if (is.finite(err[best]))
        paste0(", and that error propagates into variogram ranges, block sizes, ",
               "GWR bandwidth and GP length-scales.")
      else ".",
      not_placed(cands, placed)
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
#' The sample is taken at even spacing along the points sorted by coordinate
#' (no RNG, and no dependence on the order of the rows), so the result is
#' reproducible; \code{max_n} keeps the pairwise work bounded.  Copies of one
#' location count once, and rows without a finite location are left out with
#' the empty ones, so the sample does not depend on how many rows list a
#' place.  A detailed outline, more than \code{4 * max_n} distinct vertices,
#' is thinned the same way first, along its distinct vertices sorted by
#' coordinate.
#'
#' @param x_ll An sf/sfc object in a geographic CRS.  Non-POINT geometry is
#'   reduced to representative points first, so the two distance vectors are
#'   the same length (\code{st_coordinates()} yields one row per vertex).
#'   With fewer than \code{max_n} distinct points the outline's vertices,
#'   densified along its edges, are added to those points, so that a single
#'   study-area polygon is measured across its extent rather than not at
#'   all, however many rows repeat it.  Each
#'   edge is filled from the same end whichever way its ring runs, and an
#'   added point within about 1e-9 degrees of a point listed before it is
#'   left out, so an edge two features share, or a vertex on another
#'   feature's edge, adds no pair of points a few bits apart.  A MULTIPOINT
#'   adds its member points and nothing between them.  The layer is read in
#'   x and y, and a layer of mixed types or of collections is read feature
#'   by feature (\code{.vertex_coords()}).
#' @param crs Candidate \code{sf::crs}.
#' @param max_n Maximum number of points to sample.  Default 40 (780 pairs).
#' @return Numeric worst-case \code{|d_planar / d_geodesic - 1|} over the
#'   sampled pairs 1 mm or more apart on the ellipsoid, or \code{NA} when it
#'   cannot be computed, no pair is that long, or \code{crs} cannot place a
#'   sampled point that has a location in degrees (which the caller treats
#'   as "unusable").
#' @keywords internal
#' @noRd
.crs_distance_error <- function(x_ll, crs, max_n = 40L) {
  tryCatch({
    # Empty parts are dropped before empty features: sf::st_is_empty() calls
    # a GEOMETRYCOLLECTION, MULTIPOLYGON or MULTILINESTRING empty when its
    # first part is empty, and a box listed after an empty line in one
    # GEOMETRYCOLLECTION was left out, so every candidate scored NA.  An
    # empty part inside a line feature also segfaults GEOS below (see
    # .drop_empty_parts()), and tryCatch() cannot catch a crash, so any
    # lon/lat line layer with a null part took ensure_projected() down.
    g <- .drop_empty_parts(sf::st_geometry(x_ll))
    g <- g[!sf::st_is_empty(g)]
    # Scored in x and y.  With Z or M kept, the ring key read Z and M, so the
    # edges of a POLYGON Z whose vertices differ in Z were not filled (a box
    # over the conterminous United States chose another CRS than its XY
    # copy); the interior point of a LINESTRING Z kept its Z beside XY
    # vertices, so sf::st_transform() stopped, and GEOS refuses M: every
    # candidate scored NA and the UTM zone was kept.
    if (any(.geom_dims(g) != "XY")) g <- sf::st_zm(g, drop = TRUE, what = "ZM")
    # Representative POINTS, not the geometries themselves.  st_coordinates()
    # returns one row per VERTEX, so for a polygon or line layer the projected
    # distance vector was a different length from the geodesic one: the
    # comparison recycled, R raised "longer object length is not a multiple of
    # shorter object length" at the caller, every candidate scored NA, and the
    # selection silently fell back to the UTM zone while the log line reported
    # "NA% vs NA%".  Every county-polygon layer took that path.
    if (!all(sf::st_geometry_type(g, by_geometry = TRUE) == "POINT")) {
      g_full <- g
      g <- suppressWarnings(sf::st_point_on_surface(g_full))
      # GEOS reads a feature that crosses +-180 the long way round, so its
      # interior point landed on the far side of the globe ((-0.5, -17) for a
      # box around Fiji) and the error reported for a projection accurate to
      # 0.02% on the box was 25%.  Such a feature spans more than 180 degrees
      # of raw longitude; take its spherical centroid instead.  When s2
      # refuses one of these features, the centroids are taken feature by
      # feature (.s2_centroids()): a feature s2 refuses for a repeated vertex
      # gets the centroid of the same ring, one with crossing edges that of
      # its rebuild, and only one s2 refuses after both keeps GEOS's point.
      # Taken in one call, one refusal kept GEOS's point for every wide
      # feature, and a box around Fiji with one vertex listed twice scored
      # 25% where the same box without the repeat scores 0.02%.
      bb_all <- sf::st_bbox(g_full)
      if (isTRUE(sf::st_is_longlat(g_full)) && all(is.finite(bb_all)) &&
          bb_all[["xmax"]] - bb_all[["xmin"]] > 180) {
        wide <- which(vapply(g_full, function(s) {
          b <- sf::st_bbox(s)
          isTRUE(b[["xmax"]] - b[["xmin"]] > 180)
        }, logical(1)))
        ctr <- if (length(wide))
          tryCatch(.with_s2(.s2_centroids(g_full[wide], split = TRUE)),
                   error = function(e) NULL)
        if (!is.null(ctr)) {
          took <- setdiff(seq_along(wide), ctr$refused)
          if (length(took)) g[wide[took]] <- ctr$points[took]
        }
      }
      # One point per distinct location.  A long-format layer lists one
      # polygon over many rows, and each row gave a copy of one point: from
      # 40 rows no outline was added and every pair had zero length, so every
      # candidate scored NA and the UTM zone was kept (a box over the
      # conterminous United States, 1.76% in the Lambert azimuthal CRS on one
      # row); below that the copies took places in the sample.
      g <- .distinct_located_points(g)
      # One point per feature is nothing to measure on a study-area outline:
      # a single polygon gave one point, every candidate scored NA, and the
      # selector kept the UTM zone at any extent -- a CONUS outline got zone
      # 15 (13.7% worst-case error) where its vertices score Albers at 2.4%,
      # and prep_model_data(boundary =) moved a whole analysis into it.  A
      # handful of features gives a handful of pairs, none near the edges
      # where a zone distorts most.  So below `max_n` points, add the
      # outline's own vertices, densified along each edge until there are
      # about `max_n` of them (a box has only its four corners).  Layers with
      # `max_n` distinct points or more are measured on those points.  A
      # layer of mixed types or of collections is read feature by feature
      # (.vertex_coords()): sf::st_coordinates() refuses it, and a box listed
      # as one GEOMETRYCOLLECTION scored NA for every candidate.  Exact
      # copies of a feature are read once: they add no
      # vertex, and each copy's ring was densified again, which cost 0.9 s
      # per candidate on 10,000 copies of a box.
      xy <- if (length(g) < max_n)
        tryCatch(.vertex_coords(g_full[!duplicated(unclass(g_full))]),
                 error = function(e) NULL)
      # Only about `max_n` of these points are measured (the evenly spaced
      # subsample below), so thin a detailed outline to a few times that
      # first, at even spacing along its distinct vertices sorted by
      # coordinate.  The rows of st_coordinates() run feature by feature in
      # the order of the layer, and thinned along them the vertices kept
      # changed with that order: five detailed outlines over the conterminous
      # United States got an Albers CRS in three of eight feature orders and
      # a Lambert azimuthal one in the other five.  A vertex listed twice (a
      # ring's closing vertex, a border two features share) counts once, so
      # a ring listed from another vertex is thinned the same.  unique() on
      # the whole vertex matrix, twice, and a ring key pasted for every
      # vertex made the score cost 2.8 s per candidate on a 300,000-vertex
      # outline (0.04 s before the outline was scored), and
      # ensure_projected() 14 s; one sort costs a few hundredths of a second
      # there.  A thinned outline is not densified: the rows kept are not
      # neighbours along a ring, so there is no edge between them to fill.
      # The outline is thinned when more than `4 * max_n` distinct vertices
      # remain.  Counted by row, 33 copies of a box listed 165 rows of 4
      # vertices, and the box was thinned to its corners and not filled.
      # `%in% TRUE` drops the comparisons that are NA, at non-finite
      # vertices, which the `ok` filter below would drop anyway.
      thinned <- FALSE
      if (!is.null(xy) && nrow(xy) > 0L) {
        o  <- order(xy[, 1L], xy[, 2L])
        s  <- xy[o, 1:2, drop = FALSE]
        o  <- o[c(TRUE, s[-1L, 1L] != s[-nrow(s), 1L] |
                        s[-1L, 2L] != s[-nrow(s), 2L]) %in% TRUE]
        thinned <- length(o) > 4L * max_n
        if (thinned)
          xy <- xy[sort(o[unique(round(seq(1, length(o), length.out = 4L * max_n)))]), ,
                   drop = FALSE]
      }
      if (!is.null(xy)) {
        # The ring key is read from the L columns alone, so that no other
        # coordinate can end an edge.
        lc   <- grepl("^L[0-9]+$", colnames(xy))
        ring <- if (any(lc))
          do.call(paste, as.data.frame(xy[, lc, drop = FALSE])) else rep("1", nrow(xy))
        xy <- xy[, 1:2, drop = FALSE]
        ok <- is.finite(xy[, 1L]) & is.finite(xy[, 2L])
        xy <- xy[ok, , drop = FALSE]; ring <- ring[ok]
        if (nrow(xy) > 0L) {
          # A MULTIPOINT has no edges: its member points are added as they
          # are, and nothing is filled between them.  st_coordinates() gives
          # the members of a feature one L1 value, which the ring key read as
          # a ring, and points were added on the segments between consecutive
          # members, at places the layer does not have, moving with the order
          # of the members: two MULTIPOINTs of five points scored 0.75% in
          # UTM zone 14 in one order and 0.78% in another.  The points of a
          # POINT or MULTIPOINT inside a layer of mixed types each get a ring
          # key of their own (.vertex_coords()).
          per_edge <- if (thinned || inherits(g_full, "sfc_MULTIPOINT")) 0L else
            max(0L, ceiling(max_n / max(1L, nrow(unique(xy)))) - 1L)
          if (per_edge > 0L && nrow(xy) > 1L) {
            same <- ring[-1L] == ring[-length(ring)]
            a <- xy[-nrow(xy), , drop = FALSE][same, , drop = FALSE]
            b <- xy[-1L, , drop = FALSE][same, , drop = FALSE]
            # Each edge is filled from the end with the smaller longitude (the
            # smaller latitude when the two are equal), whichever way its ring
            # runs.  An edge two features share is listed once by each, from
            # opposite ends.  Filled from each end in turn, some copies of an
            # added point differed in their last bits; unique() kept both,
            # and the pair, about 1e-9 m apart, set every candidate's score.
            # Two boxes of 12 by 15 degrees, one above the other, scored 76.9%
            # in UTM zone 13 and 86.8% in EPSG:5070, where their outline as
            # one box scores 0.5% and 10.9%, and ensure_projected() kept the
            # zone with "76.92% vs 76.92%".  Filled from the same end, the two
            # copies are the same point.
            sw <- a[, 1L] > b[, 1L] | (a[, 1L] == b[, 1L] & a[, 2L] > b[, 2L])
            ab <- a[sw, , drop = FALSE]
            a[sw, ] <- b[sw, ]
            b[sw, ] <- ab
            f <- rep(seq_len(per_edge) / (per_edge + 1), each = nrow(a))
            a <- a[rep(seq_len(nrow(a)), per_edge), , drop = FALSE]
            b <- b[rep(seq_len(nrow(b)), per_edge), , drop = FALSE]
            # Along the edge as it runs on the globe, the short way round.
            # Interpolated in raw degrees, the edge of an outline from 177 to
            # -178 was filled with points near longitude 0, and a box around
            # Fiji reported a 164% distance error for a projection accurate to
            # 0.02% on it.
            dl <- b[, 1L] - a[, 1L]
            dl <- ((dl + 180) %% 360) - 180
            lon <- ((a[, 1L] + f * dl + 180) %% 360) - 180
            lat <- a[, 2L] + f * (b[, 2L] - a[, 2L])
            # An added point can still fall on a point listed before it, with
            # other last bits: a vertex of another feature that lies on this
            # edge (a T-junction, or a corner touching the edge), or a point
            # added along an edge that overlaps this one in part.  Rounding
            # the added points alone moved them off such vertices.  An added
            # point is left out when it shares a cell of 1e-9 degrees (about
            # 0.1 mm) with an earlier point, the cells taken in each
            # coordinate from either of two grids half a cell apart.  Two
            # points less than half a cell apart in both coordinates share a
            # cell under one of the four pairings of the grids, so copies a
            # few bits apart are never kept either side of a cell's border;
            # points 1e-9 degrees or more apart in one coordinate are both
            # kept.  The vertices stay as they are: two that differ in their
            # last bits are left to the floor on a pair's length below.  The
            # added points go in
            # sorted by coordinate, so the copy kept does not depend on the
            # order of the features.
            add <- cbind(lon, lat)
            add <- add[order(add[, 1L], add[, 2L]), , drop = FALSE]
            p <- rbind(xy, add)
            near <- logical(nrow(p))
            for (kx in list(round(p[, 1L] * 1e9), floor(p[, 1L] * 1e9)))
              for (ky in list(round(p[, 2L] * 1e9), floor(p[, 2L] * 1e9)))
                near <- near | duplicated(paste(kx, ky))
            xy <- rbind(xy, add[!near[-seq_len(nrow(xy))], , drop = FALSE])
          }
          xy <- unique(xy)
          g <- c(sf::st_geometry(g), sf::st_geometry(sf::st_as_sf(
            data.frame(x = xy[, 1L], y = xy[, 2L]), coords = c("x", "y"),
            crs = sf::st_crs(g))))
        }
      }
    } else {
      # One point per distinct location, as for the points of other layers.
      # The sample of `max_n` was taken over rows: 39 stations with station
      # 1 listed 1,001 times (repeat visits) sampled 4 of the 39, and the
      # zone's error read 2.32% where the stations score 5.07%.  A row at a
      # coordinate that is not finite made .geod_distance() stop, and every
      # candidate scored NA.
      g <- .distinct_located_points(g)
    }
    n <- length(g)
    if (n < 2L) return(NA_real_)
    if (n > max_n) {
      # Evenly spaced along the points sorted by coordinate, not along the
      # rows as they come.  By row, the 40 points measured changed with the
      # order of the layer, the candidates' errors with them, and near a tie
      # the projection chosen: 400 points over the conterminous United States
      # got an Albers CRS in some row orders and a Lambert azimuthal one in
      # others, and everything measured in that CRS afterwards (a variogram
      # range, the blocks sized from it) moved with the order of the rows.
      cc  <- sf::st_coordinates(g)[, 1:2, drop = FALSE]
      ord <- order(cc[, 1L], cc[, 2L])
      g   <- g[ord[unique(round(seq(1, n, length.out = max_n)))]]
    }
    # Ellipsoidal geodesics, not sf's s2 sphere: the sphere is itself
    # 0.24-0.56% off WGS84, which is the size of the errors being ranked.
    ll  <- sf::st_coordinates(sf::st_transform(g, 4326))[, 1:2, drop = FALSE]
    ij  <- which(lower.tri(matrix(0, nrow(ll), nrow(ll))), arr.ind = TRUE)
    d_geo <- .geod_distance(ll[ij[, 2], 1], ll[ij[, 2], 2],
                            ll[ij[, 1], 1], ll[ij[, 1], 2])
    # A point placed in degrees that `crs` cannot place makes the score NA.
    # Its pairs were left out, so a UTM zone was scored on the points it
    # can place: on LINESTRING (0 0, 180 0) zone 46 scored without both
    # ends of the line, and was kept.
    # A point at a latitude beyond 90 has no place on the Earth, which no
    # CRS places, and its pairs are left out as before
    # (.located_vertices()).
    pp <- sf::st_coordinates(sf::st_transform(g, crs))[, 1:2, drop = FALSE]
    if (nrow(pp) != nrow(ll) ||
        any(is.finite(ll[, 1L]) & is.finite(ll[, 2L]) & abs(ll[, 2L]) <= 90 &
            !(is.finite(pp[, 1L]) & is.finite(pp[, 2L]))))
      return(NA_real_)
    d_prj <- as.numeric(stats::dist(pp))
    # Pairs under 1 mm are left out.  Two listings of one place that are not
    # the same double (a corner two features share, written -107.2 in one
    # and -110.6 + 3.4 in the other; a point listed twice a few bits apart;
    # (180, y) and (-180, y)) made a pair about 1e-9 m long whose planar
    # length is rounding noise, or 40,075 km in Web Mercator, and that one
    # ratio set the score of every candidate: two boxes side by side, their
    # shared corners typed that way, scored 13.7% in UTM zone 13, where
    # their outline as one box scores 0.26%.  Planar coordinates up to 2e7 m
    # from the origin carry rounding of a few 1e-9 m, a few parts in 1e6 of
    # 1 mm.  A layer with no pair 1 mm or longer scores NA, as a layer of
    # one location does.
    keep  <- is.finite(d_geo) & is.finite(d_prj) & d_geo >= 1e-3
    if (!any(keep)) return(NA_real_)
    max(abs(d_prj[keep] / d_geo[keep] - 1))
  }, error = function(e) NA_real_)
}


#' The vertices of a layer of any geometry class, with their ring keys
#'
#' \code{sf::st_coordinates()} for a layer of one geometry type.  A layer
#' that mixes types (class \code{sfc_GEOMETRY}) or holds collections
#' (\code{sfc_GEOMETRYCOLLECTION}), which \code{sf::st_coordinates()}
#' refuses, is read part by part: each collection is unpacked into its
#' members, and every part that is not a collection is a leaf, numbered in
#' the order of the layer.  The leaves are read one geometry type at a time
#' and put back in that order.  \code{L1} is the leaf, \code{L2} numbers the
#' points of a POINT or MULTIPOINT leaf (0 for other leaves), so that no two
#' of them share a ring key and nothing is filled between them, and
#' \code{L3} and \code{L4} are the leaf's own ring and part (0 where it has
#' none).  Empty parts and empty leaves are left out, and Z and M are
#' dropped from such a layer.
#'
#' @param g An sfc.
#' @return A numeric matrix with columns \code{X} and \code{Y} first, then
#'   the columns of \code{sf::st_coordinates()} or \code{L1} to \code{L4};
#'   one row per vertex.
#' @keywords internal
#' @noRd
.vertex_coords <- function(g) {
  if (!inherits(g, c("sfc_GEOMETRY", "sfc_GEOMETRYCOLLECTION")))
    return(sf::st_coordinates(g))
  if (any(.geom_dims(g) != "XY")) g <- sf::st_zm(g, drop = TRUE, what = "ZM")
  # Empty parts first: sf::st_is_empty() calls a MULTIPOLYGON or
  # MULTILINESTRING leaf empty when its first part is empty.
  g <- .drop_empty_parts(g)
  leaves_of <- function(s)
    if (inherits(s, "GEOMETRYCOLLECTION"))
      unlist(lapply(unclass(s), leaves_of), recursive = FALSE)
    else list(s)
  leaves <- unlist(lapply(unclass(g), leaves_of), recursive = FALSE)
  leaves <- leaves[!vapply(leaves, sf::st_is_empty, logical(1))]
  cols   <- c("X", "Y", "L1", "L2", "L3", "L4")
  if (!length(leaves)) return(matrix(numeric(0), 0L, length(cols), dimnames = list(NULL, cols)))
  ty   <- vapply(leaves, function(s) class(s)[2L], character(1))
  out  <- lapply(unique(ty), function(t) {
    idx <- which(ty == t)
    cc  <- sf::st_coordinates(sf::st_sfc(leaves[idx]))
    lc  <- which(grepl("^L[0-9]+$", colnames(cc)))
    # The last L column numbers the leaves of this type; those before it
    # are the leaf's own rings and parts.
    leaf <- if (length(lc)) idx[cc[, lc[length(lc)]]] else idx
    own  <- cbind(cc[, lc[-length(lc)], drop = FALSE], matrix(0, nrow(cc), 2L))[, 1:2, drop = FALSE]
    pt   <- if (t %in% c("POINT", "MULTIPOINT")) seq_len(nrow(cc)) else rep(0L, nrow(cc))
    cbind(cc[, 1L], cc[, 2L], leaf, pt, own)
  })
  out <- do.call(rbind, out)
  out <- out[order(out[, 3L]), , drop = FALSE]
  dimnames(out) <- list(NULL, cols)
  out
}


#' One point per distinct location, rows without a finite one left out
#'
#' Keyed on exact coordinates through a sort on x and then y, not
#' \code{unique()} on the coordinate matrix.  The rows kept keep their
#' order, and of the copies of one location the first is kept.  An object
#' whose coordinates are not one row per point is returned as it came.
#'
#' @param g An sfc of non-empty points.
#' @return \code{g} without copies and without rows at a coordinate that is
#'   not finite.
#' @keywords internal
#' @noRd
.distinct_located_points <- function(g) {
  cc <- tryCatch(sf::st_coordinates(g), error = function(e) NULL)
  if (is.null(cc) || nrow(cc) != length(g)) return(g)
  ok <- which(is.finite(cc[, 1L]) & is.finite(cc[, 2L]))
  if (!length(ok)) return(g[0L])
  o <- ok[order(cc[ok, 1L], cc[ok, 2L])]
  s <- cc[o, 1:2, drop = FALSE]
  keep <- o[c(TRUE, s[-1L, 1L] != s[-nrow(s), 1L] | s[-1L, 2L] != s[-nrow(s), 2L])]
  if (length(keep) == length(g)) g else g[sort(keep)]
}


#' The distinct located vertices of a layer in longitude and latitude
#'
#' Read with \code{.vertex_coords()}, so a layer of mixed types or of
#' collections is read feature by feature.  A vertex with a coordinate that
#' is not finite or lies beyond \code{.POOL_XY_MAX} is left out, and so is
#' one at a latitude beyond 90, which no CRS places; copies of one vertex
#' count once.  Left in, a swapped row at (40, -96) among 39 stations over
#' the conterminous United States failed every candidate, and the layer got
#' Web Mercator, which cannot place that row either.
#'
#' @param x An sf or sfc object in longitude and latitude.
#' @return A two-column matrix of x and y; \code{NULL} when the vertices
#'   cannot be read.
#' @keywords internal
#' @noRd
.located_vertices <- function(x) {
  xy <- tryCatch(.vertex_coords(sf::st_geometry(x)), error = function(e) NULL)
  if (is.null(xy)) return(NULL)
  xy <- xy[, 1:2, drop = FALSE]
  ok <- is.finite(xy[, 1L]) & is.finite(xy[, 2L]) &
    abs(xy[, 1L]) <= .POOL_XY_MAX & abs(xy[, 2L]) <= 90
  xy <- xy[ok, , drop = FALSE]
  if (nrow(xy) < 2L) return(unname(xy))
  o <- order(xy[, 1L], xy[, 2L])
  s <- xy[o, , drop = FALSE]
  unname(s[c(TRUE, s[-1L, 1L] != s[-nrow(s), 1L] | s[-1L, 2L] != s[-nrow(s), 2L]), ,
           drop = FALSE])
}


#' Which features have a vertex at a latitude beyond 90
#'
#' Such a vertex has no place on the Earth, and no CRS places it.  Read on
#' the second coordinate of each vertex, in degrees; a coordinate that is
#' not finite does not count here (\code{.unlocated_features()} reads it).
#'
#' @param g An sfc in longitude and latitude, in degrees.
#' @return A logical vector, one element per feature.
#' @keywords internal
#' @noRd
.beyond_pole_features <- function(g) {
  n   <- length(g)
  out <- logical(n)
  if (!n) return(out)
  bb <- sf::st_bbox(g)
  if (all(is.finite(bb[c("ymin", "ymax")])) &&
      bb[["ymin"]] >= -90 && bb[["ymax"]] <= 90) return(out)
  far <- function(y) any(is.finite(y) & abs(y) > 90)
  # A POINT is a vector, a LINESTRING or MULTIPOINT a matrix, and the other
  # types (collections included) lists of those.
  bad <- function(s) {
    if (is.list(s)) return(any(vapply(s, bad, logical(1))))
    if (is.matrix(s)) return(nrow(s) > 0L && far(s[, 2L]))
    length(s) >= 2L && far(s[2L])
  }
  vapply(unclass(g), bad, logical(1))
}


#' Whether a CRS places every one of a set of points
#'
#' The points are moved from \code{from} to \code{crs} with
#' \code{sf::st_transform()} as one MULTIPOINT.  GDAL returns a MULTIPOINT
#' with a member it cannot move as an empty geometry, and drops such a
#' vertex from a line or a ring; either way the count or a coordinate tells.
#'
#' @param xy A two-column matrix (\code{.located_vertices()}); \code{NULL}
#'   when the vertices could not be read, which places nothing that can be
#'   tested and gives \code{TRUE}.
#' @param from The CRS of \code{xy}.
#' @param crs The candidate CRS.
#' @return \code{TRUE} when every point comes back with finite coordinates.
#' @keywords internal
#' @noRd
.places_vertices <- function(xy, from, crs) {
  if (is.null(xy) || !nrow(xy)) return(TRUE)
  pp <- tryCatch(
    sf::st_coordinates(sf::st_transform(sf::st_sfc(sf::st_multipoint(xy), crs = from), crs)),
    error = function(e) NULL)
  !is.null(pp) && nrow(pp) == nrow(xy) && all(is.finite(pp[, 1:2]))
}


#' The rows of a layer that a move left without some of their location
#'
#' Compares a layer before a move between CRSs with the layer after it.  A
#' row that had a location (not empty, and no coordinate that is not finite
#' or lies beyond \code{.POOL_XY_MAX}; \code{.unlocated_features()}) and
#' comes back empty, with fewer coordinates, or with a coordinate that is no
#' longer finite was not placed by the move: \code{sf::st_transform()}
#' returns a point it cannot move as an empty point and drops such a vertex
#' from a line or a ring.
#'
#' @param before,after The layer before and after the move, row for row.
#' @return The positions of those rows, possibly none.
#' @keywords internal
#' @noRd
.rows_unplaced <- function(before, after) {
  gb  <- sf::st_geometry(before)
  ga  <- sf::st_geometry(after)
  if (length(gb) != length(ga) || !length(gb)) return(integer(0))
  had <- !sf::st_is_empty(.drop_empty_parts(gb)) & !.unlocated_features(gb)
  if (!any(had)) return(integer(0))
  gone <- sf::st_is_empty(.drop_empty_parts(ga)) & had
  lost <- function(b, a) length(a) != length(b) || any(is.finite(b) & !is.finite(a))
  # One pass over the whole layer when nothing is lost, which is nearly
  # always; row by row otherwise.
  if (!any(gone) && !lost(unlist(unclass(gb)[had]), unlist(unclass(ga)[had])))
    return(integer(0))
  rows <- which(had)
  rows[gone[rows] | vapply(rows, function(i)
    lost(unlist(unclass(gb[[i]])), unlist(unclass(ga[[i]]))), logical(1))]
}


#' Warn that a move out of longitude and latitude left rows unplaced
#'
#' @param before,after As for \code{.rows_unplaced()}.
#' @param caller,arg The function and the argument the warning names.
#' @return \code{NULL}, invisibly.  The warning, of class
#'   \code{"spatialkit_rows_unplaced"}, carries \code{rows} and \code{n}.
#' @keywords internal
#' @noRd
.warn_rows_unplaced <- function(before, after, caller = "ensure_projected", arg = "x") {
  rows <- .rows_unplaced(before, after)
  if (!length(rows)) return(invisible(NULL))
  n <- length(sf::st_geometry(before))
  .warn_classed("spatialkit_rows_unplaced",
                paste0("%s(): the move of `%s` to %s cannot place every vertex of %s ",
                       "of its %d, which %s a location in longitude and latitude; ",
                       "%s back empty or without the vertices it cannot place. ",
                       "Pass a target_crs that holds every location of the layer."),
                caller, arg, .fold_crs_label(after), .rows_text(rows), n,
                if (length(rows) == 1L) "has" else "have",
                if (length(rows) == 1L) "it comes" else "they come",
                data = list(rows = rows, n = n))
  invisible(NULL)
}


#' Rows named in a message
#'
#' @param r Row numbers.
#' @return "row 3", "rows 2, 5 and 7", or the first six and how many more.
#' @keywords internal
#' @noRd
.rows_text <- function(r) {
  k <- length(r)
  if (k == 1L) sprintf("row %d", r)
  else if (k <= 6L) sprintf("rows %s and %d", paste(r[-k], collapse = ", "), r[k])
  else sprintf("rows %s and %d more", paste(r[1:6], collapse = ", "), k - 6L)
}


#' Sort keys of polygons, read off their coordinates
#'
#' One row per polygon of \code{g}: the bounding box (\code{xmin},
#' \code{ymin}, \code{xmax}, \code{ymax}), then the means of the vertices'
#' x, y, x^2, y^2 and x * y.  The means are taken over each ring's vertices
#' without its closing vertex, so they do not depend on where a ring starts
#' or which way it runs, and they tell apart polygons that share a bounding
#' box: the two triangles of a split square, or a square and the diamond
#' inscribed in it.  Two polygons agree on all nine when they are copies of
#' one polygon, which measure the same, and otherwise only in constructed
#' cases, such as a polygon whose vertices average to the centre of its box
#' and the same polygon turned half a circle about that centre; such a pair
#' keeps its row order.  The coordinates are read off the geometry lists in
#' one pass per geometry type, about 2 microseconds a polygon.
#'
#' @param g An sfc of POLYGON and MULTIPOLYGON geometries, none empty.
#' @return A numeric matrix, \code{length(g)} rows and nine columns.
#' @keywords internal
#' @noRd
.polygon_sort_keys <- function(g) {
  nm    <- c("xmin", "ymin", "xmax", "ymax", "mx", "my", "mxx", "myy", "mxy")
  out   <- matrix(NA_real_, length(g), length(nm), dimnames = list(NULL, nm))
  types <- as.character(sf::st_geometry_type(g, by_geometry = TRUE))
  for (ty in unique(types)) {
    w <- which(types == ty)
    # The rings as one list of coordinate matrices, and the feature each
    # belongs to.  Read off the lists directly: sf::st_coordinates() alone
    # costs about 7 microseconds a polygon.
    gw    <- unclass(g)[w]
    polys <- if (ty == "MULTIPOLYGON") unlist(gw, recursive = FALSE) else gw
    feat  <- if (ty == "MULTIPOLYGON") rep(seq_along(w), lengths(gw)) else seq_along(w)
    rings <- unlist(polys, recursive = FALSE)
    f_r   <- rep(feat, lengths(polys))
    d     <- vapply(rings, dim, integer(2))
    from  <- cumsum(c(1L, d[1L, ] * d[2L, ]))[seq_len(ncol(d))]
    # A matrix unlists column by column: a ring's x, then its y.  The last
    # row of each ring is its closing vertex, a repeat of its first.
    nk <- pmax(d[1L, ] - 1L, 0L)
    v  <- unlist(rings, use.names = FALSE)
    x  <- v[sequence(nk, from = from)]
    y  <- v[sequence(nk, from = from + d[1L, ])]
    f  <- rep(f_r, nk)
    s  <- rowsum(cbind(x, y, x * x, y * y, x * y), f)
    k  <- as.integer(rownames(s))
    out[w[k], 5:9] <- s / tabulate(f, nbins = length(w))[k]
    # The box: each feature's first and last vertex in the order of x, then
    # of y.
    for (j in 1:2) {
      z  <- if (j == 1L) x else y
      o  <- order(f, z)
      fo <- f[o]
      lo <- !duplicated(fo); hi <- !duplicated(fo, fromLast = TRUE)
      out[w[fo[lo]], j]      <- z[o][lo]
      out[w[fo[hi]], j + 2L] <- z[o][hi]
    }
  }
  out
}


#' Worst-case relative area distortion of a projection over a layer
#'
#' The counterpart of \code{.crs_distance_error()} for the property a density
#' or a rate depends on.  Probe polygons are measured twice: their planar
#' area in \code{crs} and their geodesic area on the globe.  The probes are the
#' layer's own polygons when it has them, otherwise an \code{n x n} grid over
#' its bounding box.  Copies of a polygon are one probe: the probes are the
#' distinct polygons, the first copy of each kept, so a layer that lists a
#' polygon once per period, or one cell twice, gives the figure of its
#' distinct polygons.  Above \code{max_n} polygons, \code{max_n} are taken at
#' even spacing along the polygons sorted by \code{.polygon_sort_keys()} (the
#' bounding box in the layer's CRS, then the means of the vertices), and the
#' polygons with the smallest \code{xmin} and \code{ymin} and the largest
#' \code{xmax} and \code{ymax} are added, so the figure does not depend on the
#' order of the rows.  An equal-area projection makes the ratio of the two the
#' same for every probe (up to the
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
#' @param max_n Number of the layer's own polygons to measure at even spacing
#'   when it has more; up to three more are added at the extremes.
#' @param grid Logical; probe with the \code{n x n} grid over the bounding box
#'   even when \code{x} has polygons.  A single study-area polygon is one
#'   probe, and one ratio has no spread to measure.
#' @param one_by_grid Logical; a layer of one distinct polygon whose planar
#'   and geodesic areas are measured is measured on the \code{n x n} grid
#'   over its bounding box instead of giving \code{NA}.
#' @return Numeric worst-case \code{|ratio / median(ratio) - 1|}, or
#'   \code{NA} when it cannot be computed (no CRS, fewer than two probes
#'   with an area, a probe whose geodesic area is unavailable, a probe
#'   polygon \code{sf::st_is_valid()} rejects in a projected \code{crs}, a
#'   probe with a coordinate beyond \code{.LONLAT_MOVE_MAX} in a projected
#'   \code{crs}).  A probe
#'   with two vertices that fall on one point in degrees is measured there
#'   without the repeat; a probe s2 refuses even then makes the figure
#'   \code{NA}.  With \code{one_by_grid}, a layer of one distinct polygon
#'   that is measured gets the grid's figure.
#' @keywords internal
#' @noRd
.crs_area_error <- function(x, crs = NULL, n = 6L, max_n = 200L, grid = FALSE,
                            one_by_grid = FALSE) {
  tryCatch({
    if (is.null(crs)) crs <- sf::st_crs(x)
    if (is.na(crs) || is.na(sf::st_crs(x))) return(NA_real_)
    g <- sf::st_geometry(x)
    g <- g[!sf::st_is_empty(g)]
    if (!length(g)) return(NA_real_)
    types <- as.character(sf::st_geometry_type(g, by_geometry = TRUE))
    own   <- !isTRUE(grid) && all(types %in% c("POLYGON", "MULTIPOLYGON"))
    probe <- if (own) {
      # One probe per distinct polygon.  Counted by row, copies of a polygon
      # weighted the median by their number, took places in the sample above
      # `max_n`, and two copies of one polygon passed as two probes: three
      # 50 km squares of a UTM zone measured 0.93% once each and 1.04% with
      # the left one listed three times (a warning, and a refused density),
      # and one Web Mercator box listed twice measured 0 where alone it is
      # unmeasured.  Copies agree on every sort key, so only polygons whose
      # keys tie with a neighbour's are compared, exactly, and the first
      # copy of each is kept; polygons that tie on the keys without being
      # copies are both kept.  Below `max_n` polygons a failure to read the
      # keys leaves the layer as it is, as before the keys were read there.
      k <- if (length(g) > max_n) .polygon_sort_keys(g)
           else tryCatch(.polygon_sort_keys(g), error = function(e) NULL)
      tie <- FALSE
      if (!is.null(k)) {
        rk <- order(k[, 1L], k[, 2L], k[, 3L], k[, 4L], k[, 5L], k[, 6L], k[, 7L],
                    k[, 8L], k[, 9L])
        ks <- k[rk, , drop = FALSE]
        tie <- c(FALSE, rowSums(ks[-1L, , drop = FALSE] != ks[-nrow(ks), , drop = FALSE]) == 0) %in% TRUE
      }
      if (any(tie)) {
        cand <- which(tie | c(tie[-1L], FALSE))
        dup  <- logical(length(rk))
        dup[cand] <- duplicated(unclass(g)[rk[cand]])
        if (any(dup)) {
          keep <- sort(rk[!dup])
          g  <- g[keep]
          k  <- k[keep, , drop = FALSE]
          rk <- order(k[, 1L], k[, 2L], k[, 3L], k[, 4L], k[, 5L], k[, 6L], k[, 7L],
                      k[, 8L], k[, 9L])
        }
      }
      if (length(g) > max_n) {
        # Evenly spaced along the polygons sorted by their bounding boxes, not
        # along the rows: by row, which polygons were measured changed with
        # the order of the layer, and with it the figure and whether
        # summarize_by_cell(area = TRUE) refused the cells (300 cells of a
        # UTM zone, one of them 700 km east of its central meridian, were
        # refused in five of eight row orders).  The keys are read off the
        # coordinates and never fail, where st_point_on_surface() or
        # st_centroid() can fail on a polygon GEOS or s2 rejects, and a
        # failure here would make summarize_by_cell(area = TRUE) refuse the
        # cells as unmeasurable.  Polygons that share a bounding box are
        # ordered by the means of their vertices: on the box alone they kept
        # their row order, and 300 Delaunay triangles of gridded points, each
        # sharing its box with the other half of its square, were refused in
        # 2 of 21 row orders and given a density in the other 19.  The
        # polygons with the smallest xmin and ymin and the largest xmax and
        # ymax are always measured: even spacing along the sorted polygons
        # can pass over the one that reaches furthest in y, and under a conic
        # CRS 300 cells whose ratios spread by 7.5 percent, one of them far
        # north of the rest, got a density.  The keys cost about 2
        # microseconds a polygon; sf::st_bbox() on each polygon cost 8.
        ks <- k[rk, , drop = FALSE]
        at <- c(round(seq(1, length(g), length.out = max_n)),
                which.min(ks[, "xmin"]), which.min(ks[, "ymin"]),
                which.max(ks[, "xmax"]), which.max(ks[, "ymax"]))
        g[rk[sort(unique(at))]]
      } else g
    } else {
      bb <- sf::st_bbox(sf::st_transform(g, crs))
      sf::st_make_grid(sf::st_as_sfc(bb), n = c(n, n), what = "polygons")
    }
    probe  <- sf::st_transform(probe, crs)
    # A polygon GEOS calls invalid in a projected CRS has a planar area that
    # is not the area it encloses (a hole outside its shell is subtracted,
    # the lobes of a bowtie cancel), while s2 may measure it in degrees: the
    # ratio then measures the polygon, not the CRS.  Three 10 km squares of
    # UTM zone 15 beside a square with a hole outside its shell measured 40
    # percent, and ensure_projected(purpose = "area") logged it.  Such a
    # layer is left unmeasured here; ensure_projected() then measures the
    # grid over the extent, and summarize_by_cell(area = TRUE) refuses the
    # invalid cells before it gets here.
    if (own && !isTRUE(sf::st_is_longlat(crs)) &&
        !all(sf::st_is_valid(probe) %in% TRUE)) return(NA_real_)
    # A probe of a projected CRS with a coordinate beyond .LONLAT_MOVE_MAX is
    # not on the Earth and has no geodesic area, so the figure is NA, as for
    # a probe s2 refuses below.  Moved to EPSG:4326 from Web Mercator at
    # x = 1e20, it never came back from PROJ.
    if (!isTRUE(sf::st_is_longlat(crs)) && any(.beyond_lonlat_move(probe)))
      return(NA_real_)
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
    # When the one call fails, each probe is measured on its own, and a probe
    # s2 refuses is measured again without the vertices that repeat the one
    # before them in degrees.  Two corners a fraction of a nanometre apart in
    # the layer's CRS can fall on one point in degrees, and s2 refuses the
    # ring ("duplicate vertex"); an exact repeated vertex, which GEOS
    # accepts, does the same.  4 of 133 k-means cells kept in EPSG:32615 for
    # rows in degrees did, the one call made the figure NA, and
    # summarize_by_cell(area = TRUE) refused cells that are valid in their
    # own CRS.  Dropping a repeat leaves the ring's points in degrees as
    # they were, so no area moves.  st_make_valid() is not used: it snaps
    # the vertices to 1e-7 degrees, which moves an area by up to about the
    # perimeter times a centimetre, and on cells a few metres across moves
    # the ratio by more than the tolerance.
    p_ll <- sf::st_transform(probe, 4326)
    geod <- .with_s2(tryCatch(as.numeric(sf::st_area(p_ll)), error = function(e)
      vapply(seq_along(p_ll), function(i)
        tryCatch(as.numeric(sf::st_area(p_ll[i])), error = function(e)
          tryCatch(as.numeric(sf::st_area(.drop_repeated_vertices(p_ll[i]))),
                   error = function(e) NA_real_)),
        numeric(1))))
    # A probe s2 refuses even then (crossing edges, a spike) leaves the
    # figure unmeasured.  Leaving it out let the probes that carry the
    # distortion drop out of the spread: 20 Web Mercator cells from 30N to
    # 67N, 18 of them with a repeated vertex, measured 7e-15 over the two
    # left, where all twenty give 1.87.
    if (anyNA(geod)) return(NA_real_)
    ratio  <- planar / geod
    ok <- is.finite(ratio) & is.finite(geod) & geod > 0
    if (sum(ok) < 2L) {
      # One polygon is one ratio, with no spread to measure.  With
      # one_by_grid, a layer of one distinct polygon whose areas are both
      # measured is measured on the grid over its bounding box, as a layer
      # without polygons is; a polygon s2 refuses stays unmeasured above.
      if (isTRUE(one_by_grid) && own && length(g) == 1L && sum(ok) == 1L)
        return(.crs_area_error(x, crs = crs, n = n, max_n = max_n, grid = TRUE))
      return(NA_real_)
    }
    max(abs(ratio[ok] / stats::median(ratio[ok]) - 1))
  }, error = function(e) NA_real_)
}

#' A polygon without the vertices that repeat the one before them
#'
#' Each ring keeps its first vertex and every vertex whose x or y differs from
#' the one before it, so the ring's points and the area it encloses are as
#' they were, and the ring stays closed.  Other coordinates (z, m) are kept
#' and not compared.  In a GEOMETRYCOLLECTION the rings of each polygon
#' member are treated so, and the other members are kept as they were.  Used
#' by \code{.crs_area_error()} on a probe s2 refuses in degrees, and by
#' \code{.s2_mend()}.
#'
#' With \code{tol} above 0 a vertex is also dropped when neither its x nor
#' its y differs by more than \code{tol} from the last vertex kept, and a
#' vertex that close to the ring's last vertex is dropped before it; the
#' last vertex is kept, so the ring stays closed.  No vertex moves, and the
#' area changes by less than the ring's length times \code{tol}.  A ring
#' that would keep fewer than four vertices is kept as it was.  Used by
#' create_grid_polygons() on clipped cells s2 refuses in degrees.
#'
#' @param g An sfc of one POLYGON, MULTIPOLYGON or GEOMETRYCOLLECTION.
#' @param tol How far apart, in the units of the coordinates, two vertices
#'   are read as one; 0 for exact repeats only.
#' @return The same sfc with each ring's repeats dropped.
#' @keywords internal
#' @noRd
.drop_repeated_vertices <- function(g, tol = 0) {
  ring <- function(r) {
    if (nrow(r) < 2L) return(r)
    if (!(tol > 0)) {
      moved <- r[-1L, 1L] != r[-nrow(r), 1L] | r[-1L, 2L] != r[-nrow(r), 2L]
      return(r[c(TRUE, moved), , drop = FALSE])
    }
    n <- nrow(r)
    apart <- function(i, j) max(abs(r[i, 1:2] - r[j, 1:2])) > tol
    keep <- 1L
    for (k in seq_len(max(0L, n - 2L)) + 1L)
      if (apart(k, keep[length(keep)])) keep <- c(keep, k)
    while (length(keep) > 1L && !apart(keep[length(keep)], n))
      keep <- keep[-length(keep)]
    if (length(keep) + 1L < 4L) return(r)
    r[c(keep, n), , drop = FALSE]
  }
  one <- function(geom) {
    cls <- class(geom)
    out <- if (inherits(geom, "GEOMETRYCOLLECTION")) lapply(unclass(geom), one)
    else if (inherits(geom, "MULTIPOLYGON"))
      lapply(unclass(geom), function(p) lapply(p, ring))
    else if (inherits(geom, "POLYGON")) lapply(unclass(geom), ring)
    else return(geom)
    class(out) <- cls
    out
  }
  sf::st_sfc(one(g[[1L]]), crs = sf::st_crs(g))
}

#' Features of a lon/lat layer that s2 refuses, repaired one by one
#'
#' The features s2 refuses (\code{.s2_refused()}) are tried again, each on
#' its own: first without the vertices that repeat the one before them
#' (\code{.drop_repeated_vertices()}, also for the polygons inside a
#' collection), which leaves the ring's points in degrees as they were, and
#' then, with \code{split}, rebuilt with crossing
#' edges split (\code{.split_crossing_edges()}).  Every feature s2 accepts is
#' returned as it was.  Run with s2 on (\code{.with_s2()}): the rebuild is
#' s2's.
#'
#' @param g An sfc in a geographic CRS.
#' @param split Logical; try the rebuild with crossing edges split for a
#'   feature s2 still refuses without its repeated vertices.
#' @return A list: \code{g}, the layer with each feature repaired that a
#'   retry made s2 accept; \code{refused}, the positions s2 refuses after
#'   every retry (left as they were); \code{mended}, the positions repaired.
#' @keywords internal
#' @noRd
.s2_mend <- function(g, split = TRUE) {
  bad <- .s2_refused(g)
  refused <- integer(0)
  for (i in bad) {
    polygonal <- inherits(g[[i]], c("POLYGON", "MULTIPOLYGON", "GEOMETRYCOLLECTION"))
    if (polygonal) {
      d <- tryCatch(.drop_repeated_vertices(g[i]), error = function(e) NULL)
      if (!is.null(d) && !length(.s2_refused(d))) {
        g[i] <- d
        next
      }
    }
    if (isTRUE(split)) {
      f <- .split_crossing_edges(g[i])
      if (!is.null(f)) {
        g[i] <- f
        next
      }
    }
    refused <- c(refused, i)
  }
  list(g = g, refused = refused, mended = setdiff(bad, refused))
}

#' Spherical centroids of a lon/lat layer, feature by feature where s2
#' refuses one
#'
#' \code{sf::st_centroid()} of the whole layer, under the engine sf uses in
#' the session.  When that stops because s2 refuses a feature, the features
#' are repaired one by one (\code{.s2_mend()}) and the centroids taken again:
#' a feature s2 accepts as it is keeps the centroid it has, a feature s2
#' refuses only for a repeated vertex gets the centroid of the same ring, and
#' with \code{split} a feature with crossing edges gets the centroid of its
#' rebuild.  A feature s2 refuses after every retry gets an empty point and
#' its position in \code{refused}; one refusal no longer takes the centroids
#' of the other features with it.
#'
#' @param g An sfc.
#' @param split As for \code{.s2_mend()}.
#' @return A list: \code{points}, an sfc of points row for row with \code{g};
#'   \code{refused}, the positions s2 refuses after every retry; \code{error},
#'   s2's error for the features still refused (for the whole layer when
#'   none is), or \code{NULL}.  An error that is not s2's refusal of a
#'   feature of \code{g} is raised as it came.
#' @keywords internal
#' @noRd
.s2_centroids <- function(g, split = TRUE) {
  pts <- tryCatch(suppressWarnings(sf::st_centroid(g)), error = function(e) e)
  if (!inherits(pts, "error"))
    return(list(points = pts, refused = integer(0), error = NULL))
  if (!isTRUE(sf::st_is_longlat(g)) || !isTRUE(sf::sf_use_s2()) ||
      !.is_s2_ring_error(pts))
    stop(pts)
  m <- .s2_mend(g, split)
  if (!length(m$mended) && !length(m$refused)) stop(pts)
  keep <- setdiff(seq_along(g), m$refused)
  out  <- rep(list(sf::st_point()), length(g))
  if (length(keep))
    out[keep] <- as.list(suppressWarnings(sf::st_centroid(m$g[keep])))
  # The error quoted is that of the features still refused: s2's error for
  # the whole layer named the first feature it refused, which may be one a
  # retry mended (a repeated vertex beside a bowtie).
  err <- pts
  if (length(m$refused)) {
    e2 <- tryCatch(suppressWarnings(sf::st_centroid(m$g[m$refused])),
                   error = function(e) e)
    if (inherits(e2, "error")) err <- e2
  }
  list(points = sf::st_sfc(out, crs = sf::st_crs(g)), refused = m$refused,
       error = err)
}

#' The features whose planar centroid lies outside their bounding box
#'
#' GEOS weights each part of a ring by its signed area.  Where a ring's
#' edges cross (a bowtie), or a hole lies outside its shell, parts of the
#' feature count with opposite signs, and the point \code{sf::st_centroid()}
#' gives is not a centroid of the feature: a bowtie 10 km across in a
#' projected CRS, one corner 1 m off symmetric, got a point over 33,000 km
#' away.
#' The centroid of a feature lies in its bounding box, so a point outside
#' the box of a feature GEOS reports invalid is such a point.  Only the
#' centroids GEOS takes are read (a projected CRS, no CRS, or degrees with
#' s2 off): a spherical centroid can lie outside the box of its vertices in
#' degrees.
#'
#' @param g An sfc.
#' @param pts The centroids \code{sf::st_centroid()} gave, row for row with
#'   \code{g}.
#' @return A list: \code{at}, the positions in \code{g} of the features
#'   whose point lies outside their box, and \code{reasons}, GEOS's reason
#'   for each (\code{sf::st_is_valid(reason = TRUE)}).
#' @keywords internal
#' @noRd
.centroids_off_box <- function(g, pts) {
  none <- list(at = integer(0), reasons = character(0))
  if (!length(g) || (isTRUE(sf::st_is_longlat(g)) && isTRUE(sf::sf_use_s2())))
    return(none)
  ok <- tryCatch(suppressWarnings(sf::st_is_valid(g)), error = function(e) NULL)
  if (length(ok) != length(g)) return(none)
  bad <- which(!(ok %in% TRUE))
  if (!length(bad)) return(none)
  xy <- suppressWarnings(sf::st_coordinates(pts[bad]))
  if (nrow(xy) != length(bad)) {
    # An empty point has no coordinates; read the rows one by one.
    xy <- t(vapply(bad, function(i) {
      p <- unclass(pts[[i]])
      if (length(p) >= 2L) as.numeric(p[1:2]) else c(NA_real_, NA_real_)
    }, numeric(2)))
  }
  off <- vapply(seq_along(bad), function(j) {
    x <- xy[j, 1L]; y <- xy[j, 2L]
    if (!is.finite(x) || !is.finite(y)) return(FALSE)
    b <- as.numeric(sf::st_bbox(g[[bad[j]]]))
    if (!all(is.finite(b))) return(FALSE)
    tol <- 1e-9 * max(b[3L] - b[1L], b[4L] - b[2L])
    x < b[1L] - tol || x > b[3L] + tol || y < b[2L] - tol || y > b[4L] + tol
  }, logical(1))
  at <- bad[off]
  if (!length(at)) return(none)
  why <- tryCatch(suppressWarnings(sf::st_is_valid(g[at], reason = TRUE)),
                  error = function(e) NULL)
  if (length(why) != length(at)) why <- rep("invalid geometry", length(at))
  why[is.na(why)] <- "invalid geometry"
  list(at = at, reasons = unname(why))
}

#' Refuse the centroids GEOS puts outside their features
#'
#' Raises the error of class \code{"spatialkit_features_degenerate"} for
#' the features \code{.centroids_off_box()} found, with \code{site =
#' "geos_centroid"}, so a caller that raises it again in its own words
#' (\code{.pointize_refused()}) can tell it from s2's refusal.
#'
#' @param rows The rows of \code{arg}.
#' @param reasons GEOS's reason for each.
#' @param n The number of rows of \code{arg}.
#' @param caller The function the error names.
#' @param arg The argument the error names, with its backquotes.
#' @param remedy The closing sentence.
#' @return Does not return.  The condition carries \code{rows}, \code{n},
#'   \code{reasons}, \code{parent} (\code{NULL}) and \code{site}.
#' @keywords internal
#' @noRd
.raise_centroids_off <- function(rows, reasons, n, caller, arg, remedy) {
  k <- length(rows)
  shown <- seq_len(min(k, 6L))
  named <- paste0(paste(sprintf("row %d: \"%s\"", rows[shown], reasons[shown]),
                        collapse = "; "),
                  if (k > 6L) sprintf("; and %d more", k - 6L) else "")
  msg <- sprintf(paste0(
    "%s(): the %s GEOS takes of %d of the %d features of %s %s outside %s ",
    "bounding box (%s). GEOS weights each part of a ring by its signed area, so where ",
    "the edges of a ring cross, or a hole lies outside its shell, parts of the ",
    "feature count with opposite signs and the point is not a centroid of it. %s"),
    caller, if (k == 1L) "centroid" else "centroids", k, n, arg,
    if (k == 1L) "lies" else "lie", if (k == 1L) "the feature's" else "each feature's",
    named, remedy)
  stop(errorCondition(msg, class = "spatialkit_features_degenerate", call = NULL,
                      rows = rows, n = n, reasons = reasons, parent = NULL,
                      repaired = FALSE, s2_forced = FALSE, site = "geos_centroid"))
}

#' Refuse the lines the temporary projection cannot place
#'
#' Raises the error of class \code{"spatialkit_features_degenerate"} for
#' the rows of \code{arg} whose line the temporary projection of
#' \code{coerce_to_points()} returned empty, without some of its vertices,
#' or with a coordinate that is not finite.  \code{sf::st_line_sample()}
#' segfaults on an empty line, and a crash cannot be caught, so such a line
#' is never sampled.  \code{site = "projection"} lets a caller that raises
#' it again in its own words (\code{.pointize_refused()}) tell it from the
#' other refusals.
#'
#' @param rows The rows of \code{arg}.
#' @param n The number of rows of \code{arg}.
#' @param caller The function the error names.
#' @param arg The argument the error names, with its backquotes.
#' @param crs_label The CRS of the temporary projection, as
#'   \code{.fold_crs_label()} gives it.
#' @param remedy The closing sentence.
#' @return Does not return.  The condition carries \code{rows}, \code{n},
#'   \code{reasons}, \code{parent} (\code{NULL}), \code{site} and
#'   \code{crs} (\code{crs_label}).
#' @keywords internal
#' @noRd
.raise_unplaced_lines <- function(rows, n, caller, arg, crs_label, remedy) {
  k <- length(rows)
  msg <- sprintf(paste0(
    "%s(): the temporary projection to %s, which takes the midpoints of lines, ",
    "cannot place every vertex of %s of the %d features of %s: %s back empty or ",
    "without the vertices it cannot place, so no midpoint is taken. %s"),
    caller, crs_label, .rows_text(rows), n, arg,
    if (k == 1L) "it comes" else "they come", remedy)
  stop(errorCondition(msg, class = "spatialkit_features_degenerate", call = NULL,
                      rows = rows, n = n,
                      reasons = rep(sprintf("not placed in %s", crs_label), k),
                      parent = NULL, repaired = FALSE, s2_forced = FALSE,
                      site = "projection", crs = crs_label))
}

#' The vertices of a layer with each ring's closing vertex counted once
#'
#' Read feature by feature off the geometry lists, for a layer of any type,
#' mixed types and collections included.  A ring, or a line, whose last
#' vertex repeats its first gives that vertex once, so a mean of the
#' vertices does not depend on where a ring starts.  Coordinates beyond x
#' and y are not read; rows that are not finite are left out.
#'
#' @param g An sfc.
#' @return A two-column matrix of x and y.
#' @keywords internal
#' @noRd
.open_ring_xy <- function(g) {
  none <- matrix(numeric(0), 0L, 2L)
  stack <- function(l) do.call(rbind, c(list(none), l))
  xy2 <- function(m) {
    m <- unclass(m)
    if (!is.matrix(m) || ncol(m) < 2L) return(none)
    m <- m[, 1:2, drop = FALSE]
    dimnames(m) <- NULL
    m
  }
  open <- function(m) {
    m <- xy2(m)
    k <- nrow(m)
    if (k >= 2L && isTRUE(m[1L, 1L] == m[k, 1L] && m[1L, 2L] == m[k, 2L]))
      m[-k, , drop = FALSE] else m
  }
  leaf <- function(s) {
    if (inherits(s, "GEOMETRYCOLLECTION")) return(stack(lapply(unclass(s), leaf)))
    if (inherits(s, "MULTIPOLYGON"))
      return(stack(lapply(unclass(s), function(p) stack(lapply(p, open)))))
    if (inherits(s, c("POLYGON", "MULTILINESTRING"))) return(stack(lapply(unclass(s), open)))
    if (inherits(s, "LINESTRING")) return(open(s))
    if (inherits(s, "MULTIPOINT")) return(xy2(s))
    v <- as.numeric(unclass(s))
    if (length(v) >= 2L) matrix(v[1:2], 1L) else none
  }
  xy <- stack(lapply(g, leaf))
  xy[is.finite(xy[, 1L]) & is.finite(xy[, 2L]), , drop = FALSE]
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
#' \code{[-180, 180] x [-90, 90]} envelope AND the data look \emph{positively}
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
#' the precision test, and a \code{[0, 1]}-normalised planar layer passes it too.
#' When the data are planar, set the CRS explicitly, as the warning says.
#'
#' The precision test reads the coordinates of a layer that mixes XY and XYZ
#' (or XYM) geometries in XY: \code{sf::st_coordinates()} reads such a layer
#' in the dimension of its first geometry and stops on it, which made a
#' layer in degrees look planar.  It leaves empty geometries out: an empty
#' polygon stops \code{sf::st_coordinates()} too.  Empty parts go first, as
#' \code{sf::st_is_empty()} calls a GEOMETRYCOLLECTION, MULTIPOLYGON or
#' MULTILINESTRING empty when its first part is empty.  A layer that mixes
#' geometry types or holds collections, which \code{sf::st_coordinates()}
#' refuses, is read feature by feature (\code{.vertex_coords()}): the
#' refusal made a layer in degrees that mixes types look planar.
#'
#' Both tests read the coordinates that have a location: a coordinate that is
#' not finite or lies beyond \code{.POOL_XY_MAX} (the test
#' \code{.values_at_xy()} applies) is left out.  One such row made the
#' layer's box infinite or put it off the envelope, and the whole layer was
#' then read as planar.  The coordinates are read for the box only when the
#' layer's own box is not finite or reaches beyond the bound.
#'
#' @param x An sf/sfc object with no CRS.
#' @return A list: \code{lonlat} (logical) and \code{bb} (the bounding box
#'   of the coordinates that have a location, or \code{NULL} when it could
#'   not be taken or no coordinate has one).
#' @keywords internal
#' @noRd
.looks_like_lonlat <- function(x) {
  bb <- try(sf::st_bbox(x), silent = TRUE)
  if (inherits(bb, "try-error")) return(list(lonlat = FALSE, bb = NULL))
  # The coordinates of the layer in XY, without its empty geometries, the
  # ones without a location left out; NULL when they cannot be read.
  located_xy <- function() {
    xc <- x
    if (inherits(x, c("sf", "sfc"))) {
      xc <- sf::st_geometry(x)
      xc <- .drop_empty_parts(xc)
      xc <- .xy_if_mixed(xc[!sf::st_is_empty(xc)])
    }
    cc <- try(.vertex_coords(xc), silent = TRUE)
    if (inherits(cc, "try-error")) return(NULL)
    keep <- is.finite(cc[, 1L]) & is.finite(cc[, 2L]) &
            abs(cc[, 1L]) <= .POOL_XY_MAX & abs(cc[, 2L]) <= .POOL_XY_MAX
    cc[keep, 1:2, drop = FALSE]
  }
  xy <- NULL
  if (!all(is.finite(bb)) || any(abs(bb) > .POOL_XY_MAX)) {
    xy <- located_xy()
    if (is.null(xy) || !nrow(xy)) return(list(lonlat = FALSE, bb = NULL))
    bb <- sf::st_bbox(c(xmin = min(xy[, 1L]), ymin = min(xy[, 2L]),
                        xmax = max(xy[, 1L]), ymax = max(xy[, 2L])),
                      crs = sf::st_crs(bb))
  }
  in_env <- bb["xmin"] >= -180 && bb["xmax"] <= 180 &&
            bb["ymin"] >= -90  && bb["ymax"] <= 90
  if (!in_env) return(list(lonlat = FALSE, bb = bb))

  x_range <- bb["xmax"] - bb["xmin"]
  y_range <- bb["ymax"] - bb["ymin"]
  has_large_extent <- (x_range > 1) || (y_range > 1)
  has_geo_precision <- FALSE
  if (!has_large_extent) {
    coords_mat <- if (is.null(xy)) located_xy() else xy
    if (!is.null(coords_mat) && nrow(coords_mat) > 0) {
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


#' Move a layer to another CRS with its rows without a location left empty
#'
#' A row with a coordinate that is not finite or lies beyond
#' \code{.POOL_XY_MAX} (\code{.unlocated_features()}) is read as given,
#' before the move, and comes back empty, in its own type and dimension
#' (\code{.empty_unlocated()}).  Moved as it is, a polygon or line with one
#' such vertex lost that vertex in \code{sf::st_transform()} and came back
#' as a smaller finite shape: a square with a vertex at longitude 1e200 came
#' back as a triangle.  From a projected CRS such a vertex can also come
#' back finite: a northing of 1e200 in EPSG:32616 lands at a finite point in
#' EPSG:32617.  Every reader that moves a caller's layer moves it here:
#' \code{ensure_projected()}, \code{harmonize_crs()},
#' \code{.transform_or_stamp()}, \code{.align_crs()}, \code{.cover_follow()}
#' (the builders' `cover`), the move of the polygons into a scored record's
#' CRS in \code{assign_features_to_polygons()}, and
#' \code{.values_at_free_cells()}.  The move itself is
#' \code{.lonlat_of_placed()}: a row of a projected layer with a coordinate
#' beyond \code{.LONLAT_MOVE_MAX} is emptied before a move to a CRS in
#' degrees, which from Web Mercator at x = 1e20 never came back from PROJ.
#' Any other layer whose rows all have a location is moved by
#' \code{sf::st_transform()} alone.  Points held in a column of the general
#' class \code{sfc_GEOMETRY} are moved in \code{sfc_POINT}
#' (\code{.typed_points()}): \code{sf::st_transform()} stops with "Not
#' compatible with STRSXP" on the column \code{sf::st_cast(x, "GEOMETRY")}
#' returns, and such a layer comes back in \code{sfc_POINT}.
#'
#' @param x An sf or sfc object with a CRS.
#' @param crs The CRS to move it to.
#' @return \code{x} in \code{crs}, with those rows empty.
#' @keywords internal
#' @noRd
.transform_keep_unlocated <- function(x, crs) {
  x <- .typed_points(x)
  rows <- which(.unlocated_features(sf::st_geometry(x)))
  .empty_unlocated(.lonlat_of_placed(x, crs), rows)
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
#'     standard parallels do not degenerate or pass a pole) an Albers conic are
#'     each scored by projecting representative points of the data (a non-POINT layer is
#'     reduced to one point per feature, plus the vertices of its outline when
#'     fewer than 40 distinct points remain, so a single study-area polygon is
#'     scored too, however many rows repeat it; copies of one location count
#'     once, and rows without a finite location are left out) and comparing
#'     planar with geodesic distances over the pairs at least 1 mm apart, and
#'     the one that distorts least is used. A candidate must place every
#'     location of the layer to be used: one under which a vertex of the
#'     layer has no coordinates (a UTM zone about 90 degrees or more from its
#'     central meridian) gets `distance_error` `NA` and is not chosen, and
#'     the log line names it. A vertex at a latitude beyond 90, which no CRS
#'     places, takes no part in this test. When no candidate places every
#'     vertex, Web Mercator (Equal Earth for `purpose = "area"`) is used.
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
#'     are. The box is the one of the coordinates that have a location (see
#'     the next item). With `target_crs` supplied there is no source CRS to reproject
#'     from, so the same heuristic decides between two outcomes: lon/lat-looking
#'     coordinates are read as EPSG:4326 and reprojected to the target (a real
#'     warning), and anything else has the target **stamped on without
#'     reprojection**. That is a relabel, logged only, so verify the
#'     coordinates really are in that CRS. Set the CRS explicitly to suppress
#'     either.}
#'   \item{Rows without a location}{A coordinate that is not finite or lies
#'     beyond 1e150 has no location, and nor does a row whose geometry is
#'     empty. The lon/lat test of a layer without a CRS reads the
#'     coordinates that have a location, vertex by vertex: a polygon with one
#'     such vertex is read on its other vertices. The choice of a projection
#'     reads the rows that have a location: a row with one such coordinate
#'     takes no part in it, as an empty row would, and a layer in which no
#'     row has a location gets Web Mercator (Equal Earth for
#'     `purpose = "area"`). Moved to another CRS, these rows come back
#'     empty, in their own geometry type;
#'     [sf::st_transform()] would drop the vertex it cannot move, or put a
#'     coordinate beyond 1e150 at a finite place, and return the rest of the
#'     feature. A row that has a location and that the move out of lon/lat
#'     cannot place in full comes back empty, or without the vertices the
#'     CRS cannot place, and a warning of class `spatialkit_rows_unplaced`
#'     names such rows (`rows` on the condition). With the projection chosen
#'     here that happens only at a coordinate no CRS places, such as a
#'     latitude beyond 90; with a `target_crs`, wherever that CRS cannot
#'     place the layer. The choice of a projection does not read a row with
#'     a vertex at a latitude beyond 90 either: it is made on the other
#'     rows, and a layer each row of which is empty, has no location or has
#'     such a vertex gets Web Mercator (Equal Earth for `purpose = "area"`).}
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
#'   parallels do not degenerate or pass a pole, whichever distorts distances
#'   less.  A UTM zone (conformal, not equal-area) never enters that
#'   comparison, and global
#'   coverage gets Equal Earth in place of Web Mercator.  Already-projected
#'   input is still returned untouched, but its area distortion over the
#'   extent is measured (the spread of planar-to-geodesic area ratios over
#'   probe polygons) and logged as a warning when it exceeds 1 percent.  The
#'   probes are the layer's distinct polygons (exact copies of a polygon
#'   count once), or a 6 by 6 grid over its bounding box for a layer without
#'   polygons, a layer of one distinct polygon, which gives one ratio and no
#'   spread, and a layer with a polygon s2 cannot measure in degrees (edges
#'   that cross, a spike) or that [sf::st_is_valid()] rejects in its CRS (a
#'   hole outside its shell), whose planar area is not the area it encloses.
#'   When the grid gives no figure either (a single
#'   point, a line along one axis), a warning says that the distortion could
#'   not be measured.
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
#'   A layer that mixes XY and XYZ (or XYM) geometries is returned in XY
#'   when it is transformed, with a warning of class
#'   \code{"spatialkit_zm_dropped"}: [sf::st_transform()] reads a layer in
#'   the dimension of its first geometry and misplaces the others.  The
#'   projection is then chosen from the layer in XY.  A layer in one
#'   dimension is transformed with its Z and M, and so is a layer whose
#'   only geometries of another dimension than the first are empty and not
#'   points.  A layer that is not transformed is returned with its own geometry.
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
        x0 <- .xy_if_mixed_warn(x, "ensure_projected", "x")
        x  <- .transform_keep_unlocated(x0, tcrs)
        .warn_rows_unplaced(x0, x)
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
    if (!identical(xcrs, tcrs)) {
      x0 <- .xy_if_mixed_warn(x, "ensure_projected", "x")
      x  <- .transform_keep_unlocated(x0, tcrs)
      # A move out of lon/lat that leaves a located row empty, or without
      # some of its vertices, is said.
      if (isTRUE(sf::st_is_longlat(xcrs))) .warn_rows_unplaced(x0, x)
    }
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
      # A layer that mixes dimensions is read in XY before the CRS is chosen
      # from it, as below.
      x <- .xy_if_mixed_warn(x, "ensure_projected", "x")
      tr <- .pick_local_projected_crs(x, purpose = purpose)
      x0 <- x
      x  <- .transform_keep_unlocated(x0, tr$crs)
      .warn_rows_unplaced(x0, x)
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
    # A layer that mixes XY and XYZ (or XYM) geometries is read in XY before
    # the CRS is chosen from it: sf reads the layer in the dimension of its
    # first geometry, so the centre and the scores read a mangled ring (a
    # layer in Missouri got UTM zone 25 South), and the transform returns it.
    x <- .xy_if_mixed_warn(x, "ensure_projected", "x")
    tr <- .pick_local_projected_crs(x, purpose = purpose)
    # Where candidates are compared, one is chosen only if it places every
    # located vertex.  A located row can still come back empty, or without
    # some of its vertices, from the paths that compare nothing (a zone on a
    # local extent, the global fallback) at a coordinate no CRS places, such
    # as a latitude beyond 90.  sf::st_transform() does that without a
    # word; it is said here.
    x0 <- x
    x  <- .transform_keep_unlocated(x0, tr$crs)
    .warn_rows_unplaced(x0, x)
    attr(x, "crs_choice") <- tr$candidates
  } else if (purpose == "area") {
    # Projected input is never reprojected, but when areas are what matter
    # the caller should know whether this CRS can deliver them.  Measured, as
    # the distance choice is: the spread of planar-to-geodesic area ratios
    # over probe polygons on the extent.
    # Measured on the layer in XY when it mixes dimensions; the layer is
    # returned as it came.  A layer of one distinct polygon gives one ratio
    # and no spread, and a polygon s2 refuses in degrees (edges that cross,
    # a spike) leaves the polygons unmeasured; the CRS is then measured on
    # the grid over the extent, as for a layer without polygons.  So is a
    # layer with a polygon GEOS calls invalid in its CRS, whose planar area
    # is not the area it encloses.  The first two were
    # passed over without a word: a box from 100W to 90W and 30N to 60N in
    # Web Mercator logged nothing, where split in two it logs 26.2%.
    xy  <- .xy_if_mixed(x)
    err <- .crs_area_error(xy)
    if (!is.finite(err)) err <- .crs_area_error(xy, grid = TRUE)
    if (!is.finite(err))
      .log_warn(
        paste0("ensure_projected(purpose = \"area\"): the area distortion of the ",
               "input's CRS (%s) across this extent could not be measured: no two ",
               "probes, among its polygons or the cells of a 6 by 6 grid over its ",
               "bounding box, have both a planar and a geodesic area. It is ",
               "returned as it is."),
        .fold_crs_label(x))
    else if (err > .area_error_tol)
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
#' @return A named list with components a and b.  A layer that is
#'   transformed and mixes XY and XYZ (or XYM) geometries is returned in XY,
#'   with a warning of class \code{"spatialkit_zm_dropped"}:
#'   [sf::st_transform()] reads a layer in the dimension of its first
#'   geometry and misplaces the others.  A layer in one dimension is
#'   transformed with its Z and M, and so is a layer whose only geometries
#'   of another dimension than the first are empty and not points.  A row
#'   of a transformed layer with a coordinate that is not finite or lies
#'   beyond 1e150 comes back empty, in its own geometry type: it has no
#'   location to move, and [sf::st_transform()] would drop that vertex and
#'   return the rest of the feature, or put the coordinate at a finite
#'   place.  So does a row of a projected layer moved to a CRS in degrees
#'   with a coordinate beyond 1e10 in its CRS's own unit, further out than
#'   any projected CRS places a point of the Earth.  A layer that is not
#'   transformed is returned as it came.
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

  # A layer that mixes XY and XYZ (or XYM) geometries is transformed, and
  # returned, in XY, with a warning: sf reads the layer in the dimension of
  # its first geometry, so its transform returns an XYZ ring among XY rows
  # as a ring of other points, or stops with "OGR error".
  .safe_transform <- function(x, to, which_obj) {
    x <- .xy_if_mixed_warn(x, "harmonize_crs", which_obj)
    # A row with a coordinate that is not finite or lies beyond 1e150 comes
    # back empty (.transform_keep_unlocated()): moved as it is, a polygon
    # lost that vertex and came back as a smaller finite shape.  A row of a
    # projected layer with a coordinate beyond .LONLAT_MOVE_MAX comes out
    # empty when the target is in degrees (.lonlat_of_placed(), which that
    # helper moves with): moved from Web Mercator at x = 1e20 it never came
    # back from PROJ.
    tryCatch(
      .transform_keep_unlocated(x, to),
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
      return(.safe_transform(sf::st_set_crs(x, 4326), to, which_obj))
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
      a <- .safe_transform(a, target_crs, "a")
    }
    if (is.na(crs_b)) {
      b <- .resolve_crsless(b, "b", target_crs)
    } else {
      b <- .safe_transform(b, target_crs, "b")
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
    b <- .safe_transform(b, crs_a, "b")
  } else {
    a <- .safe_transform(a, crs_b, "a")
  }
  list(a = a, b = b)
}

# -----------------------------------------------------------------------------
# Geometry Helpers
# -----------------------------------------------------------------------------

#' The coordinate dimensions of a layer's geometries
#'
#' Read from each geometry: sf takes the dimension of a layer from its first
#' geometry, so its own test cannot tell a layer that mixes XY and XYZ
#' geometries from one that does not.  An empty geometry has a dimension
#' too.  \code{.drop_zm()} counts the values of a point layer instead, which
#' cannot tell a mixed layer from a uniform one either.
#'
#' @param x An sf or sfc object.
#' @return The distinct dimensions, among \code{"XY"}, \code{"XYZ"},
#'   \code{"XYM"} and \code{"XYZM"}, in the order they first occur.
#' @keywords internal
#' @noRd
.geom_dims <- function(x) {
  g <- if (inherits(x, "sfc")) x else sf::st_geometry(x)
  unique(vapply(g, function(s) class(s)[1L], character(1), USE.NAMES = FALSE))
}

#' Whether a layer mixes geometries of different dimensions
#' @param x An sf or sfc object.
#' @return \code{TRUE} when its geometries do not all have one dimension.
#' @keywords internal
#' @noRd
.mixes_zm <- function(x) length(.geom_dims(x)) > 1L

#' A layer read in x and y when it mixes dimensions
#'
#' sf reads a layer that mixes XY and XYZ (or XYM) geometries in the
#' dimension of its first: \code{sf::st_coordinates()} of such a point layer
#' lays the values out of step or stops on its \code{dimnames}, and
#' \code{sf::st_transform()} of one with an XYZ point first stops with "OGR
#' error".  Z and M are then dropped from every geometry.  A layer in one
#' dimension is returned as it came.
#'
#' @param x An sf or sfc object.
#' @return \code{x}, or its copy in XY when it mixes dimensions.
#' @keywords internal
#' @noRd
.xy_if_mixed <- function(x) {
  if (.mixes_zm(x)) sf::st_zm(x, drop = TRUE, what = "ZM") else x
}

#' A layer only read for its x and y, in XY when it mixes dimensions or
#' carries M
#'
#' For the layers a builder or a seeder reads and does not return (the
#' tessellation builders' \code{boundary}, points and cover, the seeders'
#' \code{boundary}): the cells and seeds made from them are defined in x and
#' y.  A layer that mixes dimensions is read in XY, as by
#' \code{.xy_if_mixed()}, and so is one whose geometries carry M (XYM or
#' XYZM), which GEOS refuses ("GEOS does not support XYM or XYZM
#' geometries").  A layer XY or XYZ throughout is returned as it came.
#'
#' @param x An sf or sfc object.
#' @return \code{x}, or its copy in XY.
#' @keywords internal
#' @noRd
.xy_if_mixed_or_m <- function(x) {
  d <- .geom_dims(x)
  if (length(d) > 1L || any(d %in% c("XYM", "XYZM")))
    sf::st_zm(x, drop = TRUE, what = "ZM")
  else x
}

#' A layer about to be transformed and returned, read in x and y when it
#' mixes dimensions, with a warning
#'
#' \code{sf::st_transform()} reads a layer in the dimension of its first
#' geometry: an XYZ ring among XY rows comes back as a ring of other points,
#' and an XY geometry among XYZ-first rows stops it with "OGR error".  A
#' function that returns its transformed copy of the caller's layer hands
#' it over in XY when the layer mixes dimensions, and says so with a warning
#' of class \code{"spatialkit_zm_dropped"}, which carries the rows that held
#' Z or M as \code{rows}.  A layer in one dimension is returned as it came,
#' and is transformed with its Z.  So is a layer whose only geometries of
#' another dimension than the first are empty and not points: they have no
#' coordinates to misplace, and sf transforms the rest with their Z.  An
#' empty point holds \code{NA} coordinates, which sf misreads as any other.
#'
#' @param x An sf or sfc object.
#' @param caller The function whose result it is, for the message.
#' @param arg The argument that holds the layer, for the message.
#' @param why The end of the message: what reads the layer wrongly and what
#'   comes back.
#' @return \code{x}, or its copy in XY when it mixes dimensions.
#' @keywords internal
#' @noRd
.xy_if_mixed_warn <- function(x, caller, arg,
                              why = paste0("which sf::st_transform() reads ",
                                           "wrongly, so Z and M are dropped and ",
                                           "the layer is returned in XY")) {
  g <- if (inherits(x, "sfc")) x else sf::st_geometry(x)
  d <- vapply(g, function(s) class(s)[1L], character(1), USE.NAMES = FALSE)
  if (length(unique(d)) < 2L) return(x)
  # sf writes every geometry in the dimension of the first; an empty one
  # after it has no coordinates to be written wrongly, save an empty point,
  # whose NA coordinates are written as values.
  bare <- sf::st_is_empty(g) & !vapply(g, inherits, logical(1), "POINT")
  if (!any(d[-1L] != d[1L] & !bare[-1L])) return(x)
  zm <- which(d != "XY")
  .warn_classed("spatialkit_zm_dropped",
                paste0("%s(): `%s` mixes geometries with and without Z or M ",
                       "coordinates (%d of its %d geometries carry them), %s. ",
                       "Drop them with sf::st_zm() first to silence this."),
                caller, arg, length(zm), length(d), why, data = list(rows = zm))
  sf::st_zm(x, drop = TRUE, what = "ZM")
}

#' A layer about to be reduced to points, in XY when a geometry the mode
#' reduces carries M
#'
#' For \code{coerce_to_points()}, after its flatten of a layer that mixes
#' dimensions.  GEOS does not read M, so its centroid and interior point
#' stop on a geometry that carries it.  When the mode reduces such a
#' geometry to a point (any under \code{"centroid"},
#' \code{"point_on_surface"} and \code{"line_midpoint"}; any that is not a
#' POINT under \code{"auto"}), Z and M are dropped from the whole layer
#' with a warning of class \code{"spatialkit_zm_dropped"}, which carries the
#' rows whose geometry carries M as \code{rows}.  Empty geometries are not
#' counted.  When the only geometries the mode reduces that carry M are
#' empty, Z and M are dropped from those alone, without a warning: they hold
#' no coordinates, so nothing is lost, and GEOS refuses them all the same.
#' Any other layer is returned as it came.
#'
#' @param x An sf object.
#' @param mode The mode of \code{coerce_to_points()}, after
#'   \code{"surface"} is read as \code{"point_on_surface"}.
#' @return \code{x}, or its copy in XY.
#' @keywords internal
#' @noRd
.xy_if_reduced_m <- function(x, mode) {
  g <- sf::st_geometry(x)
  d <- vapply(g, function(s) class(s)[1L], character(1), USE.NAMES = FALSE)
  has_m <- d %in% c("XYM", "XYZM")
  if (!any(has_m)) return(x)
  reduced <- if (identical(mode, "auto"))
    has_m & !vapply(g, inherits, logical(1), "POINT")
  else has_m
  if (!any(reduced)) return(x)
  has_m <- has_m & !sf::st_is_empty(g)
  if (!any(reduced & has_m)) {
    lst <- as.list(g)
    lst[reduced] <- lapply(lst[reduced], sf::st_zm, drop = TRUE, what = "ZM")
    return(sf::st_set_geometry(x, sf::st_sfc(lst, crs = sf::st_crs(g),
                                             precision = sf::st_precision(g))))
  }
  rows <- which(has_m)
  .warn_classed("spatialkit_zm_dropped",
                paste0("coerce_to_points(): %d of the %d geometries of `x` ",
                       "carry M coordinates, which GEOS does not read, so Z ",
                       "and M are dropped and the points are returned in XY. ",
                       "Drop them with sf::st_zm() first to silence this."),
                length(rows), length(d), data = list(rows = rows))
  sf::st_zm(x, drop = TRUE, what = "ZM")
}

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

#' Planar lengths and midpoints of lines, in the unit of their CRS
#'
#' \code{sf::st_length()} gives a projected line's length in its CRS's unit
#' when sf can name that unit (\code{st_crs()$ud_unit}), and in metres when
#' it cannot and the CRS carries a \code{+to_meter} (a unit of 1000 km, or
#' of 3.28 m).  \code{sf::st_line_sample()} hands that length to GDAL as a
#' distance along the line in the CRS's unit, so a line in a unit of 2 m got
#' its "midpoint" at its end, and one in a unit of 0.5 m a quarter of the way
#' along.  Both are taken here with the CRS set aside, so GEOS measures and
#' GDAL samples the planar line in the CRS's own unit, and the CRS is put
#' back on the points.  Lines in metres, in a named unit or with no CRS get
#' the numbers \code{sf::st_length()} and \code{sf::st_line_sample()} give
#' them, to the bit.
#'
#' @param x An sfc of non-empty LINESTRINGs in a projected CRS or none.
#' @return \code{.planar_line_lengths()}: a plain numeric vector, one length
#'   per line.  \code{.planar_line_midpoints()}: an sfc of POINTs in the CRS
#'   of \code{x}, one per line.
#' @keywords internal
#' @noRd
.planar_line_lengths <- function(x) {
  if (.is_longlat(x)) stop(".planar_line_lengths(): lines in degrees.", call. = FALSE)
  as.numeric(sf::st_length(sf::st_set_crs(x, NA)))
}

# The midpoints of .planar_line_lengths(): see there.
.planar_line_midpoints <- function(x) {
  if (.is_longlat(x)) stop(".planar_line_midpoints(): lines in degrees.", call. = FALSE)
  mp <- sf::st_line_sample(sf::st_set_crs(x, NA), sample = 0.5)
  sf::st_set_crs(sf::st_cast(mp, "POINT"), sf::st_crs(x))
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
#' In mode `"auto"` a MULTILINESTRING gets the midpoint of its longest part.
#' When several parts share the greatest length (to one part in ten
#' million, measured in the projection used), their midpoints are taken
#' and the one lowest in x, then in y, is kept, so the point does not depend
#' on the order in which the feature lists its parts.  x values within one
#' part in ten million of the diagonal of the midpoints' bounding box count
#' as one, so y decides between midpoints whose x differ only in their last
#' bits, and the point does not move with the origin of the coordinates.
#'
#' A layer whose geometries do not all have the same dimensions (XY beside
#' XYZ or XYM, as `rbind()` of an XY and an XYZ layer gives) is read in x
#' and y: sf hands GEOS and s2 every geometry in the dimension of the first,
#' and the points of such a layer are returned in one dimension.
#' Z and M are dropped, the points come back in XY, and a warning of class
#' `spatialkit_zm_dropped` says so, with the rows that carried Z or M as its
#' `rows`.  Mode `"bbox_center"` reads the box of each geometry and gives XY
#' points for any layer, with no warning.  A layer in one dimension, XY or
#' XYZ throughout, is read as it is, and so is a layer whose only
#' geometries of another dimension than the first are empty and not points.
#'
#' GEOS does not read M coordinates.  When the mode reduces a geometry that
#' carries M (XYM or XYZM) to a point (any such geometry under
#' `"centroid"`, `"point_on_surface"` and `"line_midpoint"`, any that is
#' not a POINT under `"auto"`), the layer is reduced in x and y: Z and M
#' are dropped from every geometry, the points come back in XY, and a
#' warning of class `spatialkit_zm_dropped` says so, with the rows whose
#' geometry carries M as its `rows`.  Empty geometries are not counted.
#' A POINT layer that carries M is returned as it is under `"auto"`, and
#' `"bbox_center"` gives XY points with no warning, as for any layer.
#'
#' On a layer in longitude and latitude with [sf::sf_use_s2()] on, sf takes
#' centroids with s2 (under `"centroid"`, for the rows other than
#' LINESTRINGs under `"line_midpoint"`, and for collections under
#' `"auto"`), and s2 refuses some polygons GEOS accepts.  A polygon s2
#' refuses because a vertex repeats the one before it, on its own or inside
#' a collection, gets the centroid of the same ring without the repeat; the
#' ring's points are unchanged.  A
#' polygon s2 refuses for another reason (edges that cross, a figure eight)
#' stops the call with an error of class `spatialkit_features_degenerate`
#' that gives its rows of `x` (`rows` on the condition), quotes s2, and
#' says to transform `x` to a projected CRS or to use another `mode`.  With
#' s2 off, the centroids are GEOS's, planar in degrees, as before.
#'
#' GEOS, which takes the centroids in a projected CRS, without a CRS, and in
#' degrees with s2 off, weights each part of a ring by its signed area.
#' Where the edges of a ring cross (a bowtie), or a hole lies outside its
#' shell, parts of the feature count with opposite signs, and the point can
#' lie far outside the feature: a bowtie 10 km across, one corner 1 m off
#' symmetric, got a point over 33,000 km away.  A centroid outside the bounding
#' box of a feature GEOS reports invalid stops the call with an error of
#' the same class, `spatialkit_features_degenerate`, that gives its rows of
#' `x` (`rows`), quotes GEOS's reason for each (`reasons`), carries
#' `site = "geos_centroid"`, and says to pass `mode = "point_on_surface"`
#' or to repair the features with [sf::st_make_valid()].  A centroid
#' inside the box is returned as GEOS gives it.
#'
#' Modes `"auto"` and `"line_midpoint"` take the midpoint of a line in a
#' temporary projection (`tmp_project`), the one [ensure_projected()]
#' chooses for the lines.  A line that projection returns empty or without
#' some of its vertices is never sampled, since [sf::st_line_sample()] ends
#' the R session on an empty line: the call stops with an error of the
#' same class, `spatialkit_features_degenerate`, that gives its rows of `x`
#' (`rows`), carries `site = "projection"`, and says to pass
#' `mode = "centroid"` or `tmp_project = FALSE`, which take the point in
#' longitude and latitude, without the projection.  A line with a
#' coordinate that is not finite or lies beyond 1e150 has no location, and
#' gets an EMPTY POINT under these modes, as an empty line does.
#'
#' The functions that reduce a layer they are given to points this way
#' ([prep_model_data()], [make_folds()], [fold_separation()],
#' [estimate_sac_range()], [resolution_profile()],
#' [determine_optimal_levels()], [cell_values_at()], [summarize_by_cell()],
#' [select_features_forward()], [compare_models_cv()],
#' [area_of_applicability()], [predict_surface()]) raise each of these
#' errors in their own name, with the same class and fields, naming their
#' argument and its rows.
#'
#' @param x An sf object.
#' @param mode One of "auto", "centroid", "point_on_surface", "surface",
#'   "line_midpoint", "bbox_center".
#' @param tmp_project Logical; temporarily project for line-based midpoints.
#'   When \code{x} has no CRS and the lon/lat heuristic of
#'   \code{\link{ensure_projected}()} takes its coordinates for degrees
#'   (inside the lon/lat envelope and more than one unit across, or with
#'   decimal-degree precision), that temporary projection interprets them as
#'   EPSG:4326 (with a warning) and the midpoints returned are geodesic ones
#'   brought back to the input's numbers, not planar midpoints.  Set the CRS,
#'   or pass \code{tmp_project = FALSE}, for planar data.
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

  crs <- sf::st_crs(x)

  # -- bbox_center ---
  # Read geometry by geometry, so a layer that mixes dimensions gives the
  # boxes of its rows; the points are XY whatever the input.
  if (mode == "bbox_center") {
    return(sf::st_set_geometry(x, .bbox_center_sfc(x)))
  }

  # A layer that mixes XY and XYZ (or XYM) geometries is read in x and y:
  # sf hands GEOS and s2 every geometry in the dimension of the first, so an
  # XYZ ring among XY rows (or the reverse) stopped GEOS or was read as a
  # ring of other points, and st_transform() of the lines stopped with "OGR
  # error".  The points come back in XY, with the warning of class
  # spatialkit_zm_dropped.  A layer whose geometry types each keep one
  # dimension (XY polygons beside a POINT Z) is read in x and y too: the
  # types are read apart below, but the points would come back mixed.  The
  # warning therefore says what is done, not that the layer was misread.  A
  # layer in one dimension takes the path below as it did, and so does one
  # whose only geometries of another dimension than the first are empty and
  # not points (see .xy_if_mixed_warn()).
  x <- .xy_if_mixed_warn(x, "coerce_to_points", "x",
                         why = paste0("so it is read in x and y: Z and M are ",
                                      "dropped and the points are returned ",
                                      "in XY"))
  # GEOS does not read M: its centroid, interior point and the rest stopped
  # with sf's "GEOS does not support XYM or XYZM geometries" on a layer whose
  # geometries carry M.  When the mode reduces such a geometry to a point
  # (every one under "centroid", "point_on_surface" and "line_midpoint",
  # every one that is not a POINT under "auto"), Z and M are dropped from
  # the whole layer and the points come back in XY, with the warning of
  # class spatialkit_zm_dropped.  GEOS's points of a POLYGON Z are XY too,
  # so no measure is lost that a point could carry.  POINT rows under
  # "auto" are returned as they are, with their M.  Empty geometries hold no
  # coordinates and are not counted; when they are the only ones the mode
  # reduces that carry M, Z and M are dropped from them alone, without a
  # word, since GEOS refuses an empty XYM polygon too.
  x <- .xy_if_reduced_m(x, mode)
  g <- sf::st_geometry(x)

  # Centroids of the rows `rows` (positions in g; sub is g[rows]).  On
  # lon/lat input with s2 on, s2 takes them, and a polygon s2 refuses
  # stopped the call with s2's bare message ("Loop 0 is not valid: Edge 1 is
  # degenerate (duplicate vertex)"), which names neither the row nor a way
  # on, even for a polygon GEOS accepts.  A polygon s2 refuses for a vertex
  # that repeats the one before it, also inside a collection, gets the
  # centroid of the same ring without the repeat (.s2_centroids()); a
  # collection was not retried, so under "auto" and "centroid" it stopped
  # the call.  One s2 refuses even then (crossing
  # edges, a figure eight) stops the call with the features' class, naming
  # its rows of `x`.  It is not repaired: its centroid would be that of
  # another shape.  The message offers mode "auto" only when every row
  # refused is a POLYGON or MULTIPOLYGON, which "auto" reduces with a point
  # on the surface; "auto" takes the centroid of a collection, through this
  # function, and refuses the same row again.  The case is kept on the
  # condition (auto_helps) for the callers that raise it again.
  centroids_of <- function(sub, rows) {
    sc <- .s2_centroids(sub, split = FALSE)
    if (length(sc$refused)) {
      bad  <- rows[sc$refused]
      why  <- .s2_refused_rows(sub[sc$refused])
      auto_helps <- all(as.character(sf::st_geometry_type(
        sub[sc$refused], by_geometry = TRUE)) %in% c("POLYGON", "MULTIPOLYGON"))
      .raise_cells_degenerate(
        sc$error, bad,
        if (!is.null(why) && length(why$reasons) == length(bad)) why$reasons
        else rep(conditionMessage(sc$error), length(bad)),
        length(g), "coerce_to_points", "`x`", "sf::st_centroid()",
        "call coerce_to_points() again; the centroids are then taken in that CRS",
        of = sprintf("the %d features of `x`", length(g)), kind = "features",
        more = if (mode == "auto" || !auto_helps)
          "Or pass mode = \"point_on_surface\", which takes a point on each geometry with GEOS."
        else
          "Or pass mode = \"auto\", which takes a point on the surface of each polygon with GEOS.",
        auto_helps = auto_helps)
    }
    # A centroid GEOS puts outside its feature (a ring whose edges cross,
    # a hole outside its shell) is refused, naming its rows of `x`; it was
    # returned without a word, for a bowtie over 33,000 km away.
    off <- .centroids_off_box(sub, sc$points)
    if (length(off$at))
      .raise_centroids_off(
        rows[off$at], off$reasons, length(g), "coerce_to_points", "`x`",
        paste0("Pass mode = \"point_on_surface\", which takes a point on each geometry ",
               "with GEOS, or repair the features with sf::st_make_valid() and call ",
               "coerce_to_points() again."))
    sc$points
  }

  # -- direct ops ---
  # Not spherical-safe, whatever this heading used to say.  On lon/lat input
  # st_centroid() is spherical only while sf_use_s2() is TRUE: with it FALSE
  # it is planar in degrees (a box -120..-60 x 50..75 got a centre 155 km
  # from the s2 one) and its warning is suppressed here.  st_point_on_surface()
  # is GEOS, planar in degrees under either setting.
  if (mode == "centroid") {
    return(sf::st_set_geometry(x, centroids_of(g, seq_along(g))))
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
  # The temporary projection of the lines whose midpoints are taken, rows
  # `rows` of `x`.  A line it returns empty, without some of its vertices,
  # or with a coordinate that is not finite is refused, naming its rows of
  # `x`: st_line_sample() segfaults on an empty line, and a crash cannot be
  # caught.  LINESTRING (0 0, 180 0) in UTM zone 46 ended the R session.
  # ensure_projected()'s warning about such rows would name the rows of this
  # subset, so it is not given; the refusal names them.
  projected_lines <- function(lines_sf, rows) {
    done <- .answer_for_classed("spatialkit_rows_unplaced")
    proj <- tryCatch(ensure_projected(lines_sf), finally = done())
    bad  <- .rows_unplaced(lines_sf, proj)
    if (length(bad))
      .raise_unplaced_lines(
        rows[bad], length(g), "coerce_to_points", "`x`", .fold_crs_label(proj),
        paste0("Pass mode = \"centroid\", or tmp_project = FALSE, which take the ",
               "point in longitude and latitude, without the projection."))
    proj
  }

  # A line with a coordinate that is not finite or lies beyond .POOL_XY_MAX
  # has no location: it is read as an empty line, as elsewhere.  In lon/lat
  # the temporary projection emptied it, and st_line_sample() segfaulted on
  # the empty copy.
  .line_midpoints <- function(idx_ls) {
    res  <- rep(list(sf::st_point()), length(idx_ls))
    full <- which(!sf::st_is_empty(g[idx_ls]) & !.unlocated_features(g[idx_ls]))
    if (!length(full)) return(res)
    g_ls_sf   <- sf::st_sf(geometry = g[idx_ls[full]])
    g_ls_proj <- if (tmp_project) projected_lines(g_ls_sf, idx_ls[full]) else g_ls_sf
    midps     <- .planar_line_midpoints(sf::st_geometry(g_ls_proj))
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
      ctr <- centroids_of(g[idx_other], idx_other)
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
      # Empty parts are stripped before anything else sees them (see
      # .drop_empty_parts()): how st_transform() and st_cast() treat an empty
      # part varies across sf, GDAL and GEOS versions -- on macOS builds a
      # feature holding one came back EMPTY as a whole, losing its real part.
      g_mls     <- .drop_empty_parts(g[idx_mls])
      g_mls_sf  <- sf::st_sf(geometry = g_mls)
      # A located row the projection cannot place is refused, as for a
      # LINESTRING; it came back as an empty point.  A row without a
      # location comes back empty from the projection and gets an empty
      # point below.
      g_mls_proj <- if (tmp_project) projected_lines(g_mls_sf, idx_mls) else g_mls_sf
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
        # The parts are read off the structure rather than st_cast(), for the
        # same reason as above.
        mats  <- Filter(function(m) is.matrix(m) && nrow(m) > 0L,
                        unclass(proj_geom[[j]]))
        if (length(mats) == 0L) {
          out[[idx_mls[j]]] <- sf::st_point()
        } else {
          parts <- sf::st_sfc(lapply(mats, sf::st_linestring), crs = proj_crs)
          lens <- .planar_line_lengths(parts)
          k    <- if (length(lens)) which.max(lens) else 1L
          # Parts of equal greatest length are all sampled, and the midpoint
          # lowest in x, then in y, is kept: which.max() alone took the first
          # in the feature's part order, so the same geometry with its parts
          # listed the other way round got another point.  Lengths within
          # .TIE_TOL of the greatest are equal: congruent parts differ in
          # their last bits after a shift of the origin.  So do midpoints of
          # one x reached from other coordinates, and x values within
          # .TIE_TOL of the diagonal of the midpoints' box are one column
          # (.near_columns()), so y decides between them.
          k    <- if (length(k)) which(lens >= lens[k] * (1 - .TIE_TOL)) else k
          mp   <- .planar_line_midpoints(parts[k])
          mp   <- .back_to_input_crs(mp, g_mls_proj, crs)
          if (length(mp) > 1L) {
            # The column of x, then y, then Z or M where the parts carry
            # them, then x.
            xy <- unname(as.data.frame(sf::st_coordinates(mp)))
            mp <- mp[do.call(order, c(list(.near_columns(xy[[1L]], xy[[2L]])),
                                      xy[-1L], xy[1L]))[1L]]
          }
          out[[idx_mls[j]]] <- mp[[1]]
        }
      }
    }
  }

  # --- Anything else: centroid fallback ---
  idx_other <- which(!gtypes %in% c("POINT", "MULTIPOINT", "POLYGON",
                                     "MULTIPOLYGON", "LINESTRING", "MULTILINESTRING"))
  if (length(idx_other)) {
    ctr <- centroids_of(g[idx_other], idx_other)
    out[idx_other] <- as.list(ctr)
  }

  sf::st_set_geometry(x, sf::st_sfc(out, crs = crs))
}
