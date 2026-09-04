# -----------------------------------------------------------------------------
# CRS Selection
# -----------------------------------------------------------------------------

#' Pick a sensible local projected CRS for an sf/sfc object
#'
#' Chooses an appropriate projected coordinate reference system for spatial
#' data, favouring a UTM zone based on the dataset's geographic centroid when
#' the extent is narrow enough for one.
#'
#' Beyond about 5 degrees from the central meridian of the candidate zone, the
#' zone, a Lambert azimuthal equal-area centred on the data and (unless its
#' cone constant degenerates) an Albers conic are each scored by
#' \code{.crs_distance_error()} -- projecting representative points of the data
#' and comparing planar with geodesic pairwise distances -- and the least
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
#' extent; only genuinely global coverage falls back to Web Mercator
#' (EPSG:3857), which would otherwise SPLIT a wrapped layer.
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
      return(sf::st_crs(sprintf(
        "+proj=laea +lat_0=%f +lon_0=%f +datum=WGS84 +units=m +no_defs",
        lat, lon_ctr)))
    }

    .log_warn(
      paste0(".pick_local_projected_crs(): longitude extent spans %.1f deg and ",
             "the coordinates do not straddle the antimeridian, so this is ",
             "global coverage; no local projection fits it. Falling back to ",
             "EPSG:3857. Pass target_crs to ensure_projected() to choose a ",
             "projection suited to your extent."),
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

  utm_epsg <- if (!is.na(lat) && lat < 0) 32700 + cand_zone else 32600 + cand_zone

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
      return(sf::st_crs(utm_epsg))
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
    return(cands[[best]]$crs)
  }

  # Reuse the candidate zone computed above so the zone named in the messages
  # and the EPSG code returned here can never diverge.
  sf::st_crs(utm_epsg)
}


#' Worst-case relative distance error of a projection on a point layer
#'
#' Projects a deterministic sample of \code{x_ll} (lon/lat) into \code{crs} and
#' compares the planar pairwise distances with the geodesic ones, returning the
#' largest relative discrepancy.  Used by \code{.pick_local_projected_crs()} to
#' choose between candidate projections by measurement rather than by rule of
#' thumb.
#'
#' The sample is taken by evenly spaced index (no RNG), so the answer is
#' reproducible; \code{max_n} keeps the pairwise work bounded.
#'
#' @param x_ll An sf/sfc object in a geographic CRS.  Non-POINT geometry is
#'   reduced to representative points first, so the two distance vectors are
#'   the same length (\code{st_coordinates()} yields one row per vertex).
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
    # "NA% vs NA%".  Every county-polygon layer took that path.
    if (!all(sf::st_geometry_type(g, by_geometry = TRUE) == "POINT"))
      g <- suppressWarnings(sf::st_point_on_surface(g))
    n <- length(g)
    if (n < 2L) return(NA_real_)
    if (n > max_n) g <- g[unique(round(seq(1, n, length.out = max_n)))]
    d_geo <- suppressMessages(as.numeric(sf::st_distance(g)))
    d_prj <- as.numeric(stats::dist(
      sf::st_coordinates(sf::st_transform(g, crs))[, 1:2, drop = FALSE]))
    # st_distance() returns the full matrix, dist() the lower triangle.
    d_geo <- as.numeric(stats::as.dist(matrix(d_geo, nrow = length(g))))
    keep  <- is.finite(d_geo) & is.finite(d_prj) & d_geo > 0
    if (!any(keep)) return(NA_real_)
    max(abs(d_prj[keep] / d_geo[keep] - 1))
  }, error = function(e) NA_real_)
}

# -----------------------------------------------------------------------------
# Projection Enforcement
# -----------------------------------------------------------------------------

#' Return sampled points to the coordinate space of the input
#'
#' The line-midpoint branches project temporarily so that "halfway along the
#' line" is measured in a length unit, then bring the midpoints back.  When
#' the input had NO CRS and \code{ensure_projected()} interpreted it as
#' lon/lat, "back" is EPSG:4326 -- the space the input's numbers were in --
#' with the CRS then stripped again so the output matches the input.
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
#' assumed lon/lat and projected, while every \code{predict()} method -- which
#' passes a target -- stamped the fit's projected CRS onto the raw numbers.
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
#' unit across is treated as lon/lat -- a 50 m site survey included.  That is a
#' deliberate trade: a CRS-less layer is ambiguous by construction, and the
#' failure modes are not symmetric.  Reading true degrees as planar metres
#' makes every distance in the package meaningless with no way to notice;
#' reading a small planar survey as degrees produces coordinates that are
#' obviously wrong and a warning that names the assumption.  Requiring BOTH
#' tests would not help: a genuine study area 0.01 degrees across passes only
#' the precision test, and a [0, 1]-normalised planar layer passes it too.
#' Set the CRS explicitly -- the warning says so -- when the data are planar.
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
#' the training rows -- even a subset whose own bounding box would not have
#' triggered the heuristic -- lands where the training rows did.
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
#' extent of the data — it is **not** always UTM:
#'
#' \describe{
#'   \item{Local extents}{The UTM zone containing the data's centre
#'     (EPSG:326xx north of the equator, EPSG:327xx south). Distances and areas
#'     are close to true over a few degrees of longitude, which is the case
#'     this package is usually in.}
#'   \item{Wide extents}{Once the data reach well beyond the roughly 3 degrees
#'     a UTM zone is designed for, a single zone can distort distances by
#'     several percent — and that error propagates straight into variogram
#'     ranges, block sizes, GWR bandwidths and GP length-scales. Which
#'     projection is actually best is then **measured, not assumed**: the zone,
#'     a Lambert azimuthal equal-area centred on the data and (where its
#'     standard parallels do not degenerate) an Albers conic are each scored by
#'     projecting representative points of the data — a non-POINT layer is
#'     reduced to points first — and comparing planar with geodesic pairwise
#'     distances, and the one that distorts least is used.
#'     The choice, both error figures and this argument are **logged** (see the
#'     logging note under [spatialkit_quiet()]); they are not R warnings, so
#'     `tryCatch(warning = )` does not see them.}
#'   \item{Antimeridian}{Data straddling ±180° have a bounding box wider than a
#'     hemisphere. The wrap is detected from the coordinates (one very large
#'     gap in the sorted longitudes) and an equal-area projection centred on
#'     the true extent is used. Only genuinely global coverage falls back to
#'     EPSG:3857.}
#'   \item{Missing CRS}{With no `target_crs`, a bounding box that looks like
#'     lon/lat means EPSG:4326 is assumed (a real warning) and the rules above
#'     then apply; coordinates the heuristic declines are left exactly as they
#'     are. With `target_crs` supplied there is no source CRS to reproject
#'     from, so the same heuristic decides between two outcomes: lon/lat-looking
#'     coordinates are read as EPSG:4326 and reprojected to the target (a real
#'     warning), and anything else has the target **stamped on without
#'     reprojection** — a relabel, logged only, so verify the coordinates really
#'     are in that CRS. Set the CRS explicitly to suppress either.}
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
#' @return x, potentially with a new projected CRS.  CRS-less input
#'   additionally carries \code{attr(x, "crs_assumed")}: \code{"EPSG:4326"} when
#'   the lon/lat heuristic fired, \code{"none"} when it declined.  That
#'   attribute is also read on the way IN — an object already carrying
#'   \code{"none"} is returned untouched, with the heuristic skipped, which is
#'   how a \code{predict()} method replays a fit's negative decision so that a
#'   subset of the training rows is not judged differently from the whole.
#' @examples
#' library(sf)
#' pts_ll <- st_as_sf(
#'   data.frame(lon = c(9.1, 9.2), lat = c(48.7, 48.8)),
#'   coords = c("lon", "lat"), crs = 4326
#' )
#' # A local extent gets the containing UTM zone.
#' st_crs(ensure_projected(pts_ll))$epsg  # 32632
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
      tr <- .pick_local_projected_crs(x)
      x <- sf::st_transform(x, tr)
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
          .warn_and_log(
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
#'   When \code{x} has no CRS and its coordinates fall inside the lon/lat
#'   envelope, that temporary projection interprets them as EPSG:4326 (with a
#'   warning) and the midpoints returned are geodesic ones brought back to the
#'   input's numbers, not planar midpoints.  Set the CRS, or pass
#'   \code{tmp_project = FALSE}, for planar data.
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
        midps     <- .back_to_input_crs(midps, g_ls_proj, crs)
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
      midps    <- .back_to_input_crs(midps, g_ls_proj, crs)
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
