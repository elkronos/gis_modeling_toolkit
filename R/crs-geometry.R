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

  zone <- floor((lon + 180) / 6) + 1
  zone <- max(1, min(60, zone))
  epsg <- if (!is.na(lat) && lat < 0) 32700 + zone else 32600 + zone
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
#' @param x An sf or sfc object (other objects returned unchanged).
#' @param target_crs Optional target CRS (sf object, integer EPSG, or crs).
#' @return x, potentially with a new projected CRS.
#' @export
ensure_projected <- function(x, target_crs = NULL) {
  if (!(inherits(x, "sf") || inherits(x, "sfc"))) return(x)

  if (!is.null(target_crs)) {
    tcrs <- sf::st_crs(target_crs)
    if (is.na(tcrs)) return(x)
    xcrs <- sf::st_crs(x)
    if (is.na(xcrs)) {
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
            "harmonize_crs(): st_transform() failed (%s); falling back to st_set_crs(). WARNING: coordinates are NOT reprojected — downstream distances, joins, and areas may be wrong. Set on_transform_error='stop' (the default) to surface this error instead.",
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

  if (!is.null(target_crs)) {
    a <- if (!is.na(crs_a)) .safe_transform(a, target_crs) else sf::st_set_crs(a, target_crs)
    b <- if (!is.na(crs_b)) .safe_transform(b, target_crs) else sf::st_set_crs(b, target_crs)
    return(list(a = a, b = b))
  }

  if (is.na(crs_a) && is.na(crs_b)) return(list(a = a, b = b))

  if (is.na(crs_a) && !is.na(crs_b)) {
    a <- sf::st_set_crs(a, crs_b)
    return(list(a = a, b = b))
  }
  if (!is.na(crs_a) && is.na(crs_b)) {
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
#' @param x An sf object.
#' @param mode One of "auto", "centroid", "point_on_surface", "surface",
#'   "line_midpoint", "bbox_center".
#' @param tmp_project Logical; temporarily project for line-based midpoints.
#' @return An sf object with geometry coerced to POINTs.
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

  .as_sf_single <- function(geom, crs) {
    sfc1 <- if (inherits(geom, "sfc")) {
      if (is.na(sf::st_crs(geom)) && !is.na(crs)) sf::st_crs(geom) <- crs
      geom
    } else {
      sf::st_sfc(geom, crs = crs)
    }
    sf::st_sf(geometry = sfc1)
  }

  midpoint_linestring <- function(geom) {
    if (is_ll && !tmp_project) {
      return(suppressWarnings(sf::st_centroid(geom)))
    }
    x1     <- .as_sf_single(geom, crs)
    x_proj <- if (tmp_project) ensure_projected(x1) else x1
    midp   <- sf::st_line_sample(sf::st_geometry(x_proj), sample = 0.5)
    midp   <- sf::st_cast(midp, "POINT")
    if (!identical(sf::st_crs(x_proj), crs)) midp <- sf::st_transform(midp, crs)
    sf::st_geometry(midp)
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
        g_ls_sf   <- sf::st_sf(geometry = g_ls)
        g_ls_proj <- if (tmp_project) ensure_projected(g_ls_sf) else g_ls_sf
        midps     <- sf::st_line_sample(sf::st_geometry(g_ls_proj), sample = 0.5)
        midps     <- sf::st_cast(midps, "POINT")
        if (!identical(sf::st_crs(g_ls_proj), crs))
          midps <- sf::st_transform(midps, crs)
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
      out[[i]] <- sf::st_centroid(g[i])[[1L]]
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
      g_ls_sf  <- sf::st_sf(geometry = g_ls)
      g_ls_proj <- if (tmp_project) ensure_projected(g_ls_sf) else g_ls_sf
      midps    <- sf::st_line_sample(sf::st_geometry(g_ls_proj), sample = 0.5)
      midps    <- sf::st_cast(midps, "POINT")
      if (!identical(sf::st_crs(g_ls_proj), crs))
        midps <- sf::st_transform(midps, crs)
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
