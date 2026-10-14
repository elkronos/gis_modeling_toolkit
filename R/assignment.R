#' Assign features to polygons and attach a polygon ID
#'
#' Joins an sf layer of input features to a polygon layer via spatial join.
#'
#' This is the second step of the package's pipeline: it labels every
#' observation with the cell it falls in, which is what
#' [summarize_by_cell()] then aggregates over. Prefer it to [sf::st_join()]
#' when the join has to be *unambiguous*. It resolves features matching
#' several polygons by an explicit `tie_break` rule instead of silently
#' duplicating rows, so the assigned layer keeps one row per input feature
#' and cell-level counts mean what they say.
#'
#' The join runs in the CRS of `polygons_sf` whenever that CRS is projected,
#' so cell edges are the straight lines the cells were drawn with and overlap
#' areas are planar. A copy of `features_sf` is transformed for it, and the
#' features come back with the coordinates they arrived with. Otherwise (the
#' polygons are in lon/lat, or carry no CRS) the join runs in the CRS of
#' `features_sf`.
#'
#' @param features_sf An sf object containing features to assign.
#' @param polygons_sf An sf or sfc polygonal layer.
#' @param polygon_id_col Name of the polygon identifier column: the column of
#'   `polygons_sf` the IDs are read from, and the name they are returned
#'   under. Default "poly_id". When `polygons_sf` has no column of that name,
#'   the IDs are read from the first of `"poly_id"`, `"polygon_id"`, `"id"`,
#'   `"cell_id"` and `"grid_id"` it carries, and a message names the column
#'   that was used; when it has none of them, the polygons are numbered
#'   `1, 2, ...` in row order, and a message says so. The column the IDs are
#'   read from must hold one distinct, non-missing value per polygon: a
#'   duplicated or missing ID is an error, because features in different
#'   polygons would come back under one ID. That includes a unit stored as
#'   several rows that share its ID (a region split into `POLYGON` parts):
#'   merge the parts into one row per ID first, for example with
#'   `dplyr::summarise(dplyr::group_by(polygons_sf, <ID column>))`, which
#'   keeps the unit one cell. Numbering the rows would make each part a cell
#'   of its own.
#' @param keep_unassigned Logical; retain features that fall inside no
#'   polygon, carrying \code{NA} in the ID column. Default FALSE, which drops
#'   them with a warning of class \code{"spatialkit_unassigned"} that gives
#'   their number.
#' @param predicate Binary spatial predicate function. Default sf::st_intersects.
#'   Not used when \code{largest} applies: sf then assigns polygon features by
#'   overlap area and never calls the predicate.
#' @param largest Logical; when \code{features_sf} is itself polygonal, keep
#'   the polygon with the largest overlap. Default TRUE. Ignored for point and
#'   line features. A feature that only touches the polygon layer (shares an
#'   edge or a corner with it, with no overlap area) has no largest overlap
#'   and is unassigned, in any CRS; with \code{largest = FALSE} the default
#'   \code{st_intersects} counts touching, so such a feature is assigned.
#'   A feature that overlaps two or more polygons by exactly the same area
#'   (to 9 significant digits; a square split evenly across a cell edge) is
#'   given to one of them by \code{tie_break} and counted in \code{"ties"},
#'   so the choice does not depend on the order of the polygon rows.
#'   Invalid geometries (usually a self-intersecting ring),
#'   whose overlap is undefined, are repaired with \code{sf::st_make_valid()}
#'   for the join, with a warning, and returned as they arrived. If the
#'   overlap still cannot be computed the function stops: falling back to
#'   \code{predicate} and \code{tie_break} would change the rule for every
#'   feature in the layer, so pass \code{largest = FALSE} to ask for that.
#' @param tie_break Strategy for resolving features that match multiple
#'   polygons (with \code{largest}, that overlap several polygons equally):
#'   \code{"smallest_area"} (default) keeps the polygon with the
#'   smallest area and, among polygons of equal area (a point on the shared
#'   edge of two grid cells), the one whose bounding-box centre is lowest,
#'   then leftmost, so the choice does not depend on the order of the rows;
#'   \code{"first"} keeps the first match (original order-dependent
#'   behavior).
#' @param quiet Logical; suppress this function's \code{message()}s (the one
#'   naming a fallback ID column). It does not silence R warnings. Default
#'   \code{FALSE}.
#' @return An sf object with `polygon_id_col` attached, one row per input
#'   feature (fewer if `keep_unassigned = FALSE` dropped unmatched ones), in
#'   the CRS `features_sf` arrived in. A column of `features_sf` already
#'   called `polygon_id_col` is dropped before the spatial join (with a
#'   warning), so re-assigning an already-assigned layer replaces the old IDs
#'   and does not fail. Every other column is kept, including one named like
#'   the polygons' own ID column when that is read from a fallback such as
#'   `"id"` (a site `id` joined to cells keyed by `id`). If *no* feature
#'   falls inside any polygon
#'   the result is empty (or all-`NA` with `keep_unassigned = TRUE`) and a
#'   warning is raised, since the usual cause is two layers in different
#'   places (a CRS that could only be stamped, not reprojected).
#'
#'   Features that fall inside no polygon are counted on the attribute
#'   `"unassigned"`: a list with `n`, `which` (their row positions in
#'   `features_sf`) and `n_rows`. With `keep_unassigned = FALSE` they are not
#'   in the result, and a warning of class `"spatialkit_unassigned"` gives
#'   their number (`tryCatch(spatialkit_unassigned = )` catches it without
#'   matching the message; the condition carries `n` and `total`). The usual
#'   causes are a polygon layer that stops short of the data (a coastline
#'   drawn at another scale, a tessellation of seeds built without a
#'   `boundary` or `cover`; see [build_tessellation()]) and points outside
#'   the study area. With `keep_unassigned = TRUE` they are kept with an
#'   `NA` ID and no warning is raised, unless every feature is unassigned.
#'
#'   The
#'   attribute `"ties"` records how many features matched more than one
#'   polygon (with `largest`, overlapped several by exactly the same area)
#'   and had the `tie_break` rule decide for them: a list with `n`,
#'   `which` (their row positions in `features_sf`), `rule` and `n_rows`
#'   (the number of rows returned, which the record was made for). A large
#'   `n` means the polygon layer overlaps, and per-cell counts built from the
#'   result depend on the rule. The record describes the rows this call
#'   returned and does not survive subsetting: `joined[i, ]`, like
#'   `dplyr::filter()`, `slice()` or `arrange()` of it, is a plain layer
#'   with no `"ties"` attribute, so nothing reports the parent's count
#'   for a different set of rows. `sf::st_drop_geometry()` keeps the record,
#'   since the rows are the same; see \code{\link{[.spatialkit_rows}} for
#'   what binding such data frames does. The `"unassigned"` record is
#'   stamped and removed the same way.
#' @examples
#' library(sf)
#' set.seed(1)
#' pts <- st_as_sf(
#'   data.frame(x = 5e5 + runif(20, 0, 100), y = 5e6 + runif(20, 0, 100),
#'              val = rnorm(20)),
#'   coords = c("x", "y"), crs = 32632
#' )
#' bnd <- st_sf(geometry = st_sfc(st_polygon(list(rbind(
#'   c(5e5, 5e6), c(5e5 + 100, 5e6), c(5e5 + 100, 5e6 + 100),
#'   c(5e5, 5e6 + 100), c(5e5, 5e6)
#' ))), crs = 32632))
#' grid <- create_grid_polygons(bnd, target_cells = 9, type = "square")
#' assigned <- assign_features_to_polygons(pts, grid)
#' table(assigned$poly_id)
#' @family aggregation
#' @seealso [build_tessellation()] to build the polygon layer;
#'   [summarize_by_cell()] for the aggregation step that consumes the result.
#' @export
assign_features_to_polygons <- function(
    features_sf, polygons_sf, polygon_id_col = "poly_id",
    keep_unassigned = FALSE, predicate = sf::st_intersects,
    largest = TRUE,
    tie_break = c("smallest_area", "first"),
    quiet = FALSE
) {
  .msg <- function(...) if (!isTRUE(quiet)) message(...)
  if (!inherits(features_sf, "sf"))
    stop("assign_features_to_polygons(): `features_sf` must be an sf object.")
  if (!inherits(polygons_sf, "sf")) {
    if (inherits(polygons_sf, "sfc")) polygons_sf <- sf::st_as_sf(polygons_sf)
    else stop("assign_features_to_polygons(): `polygons_sf` must be sf/sfc.")
  }
  if (nrow(polygons_sf) == 0L)
    stop("assign_features_to_polygons(): `polygons_sf` has zero rows.")
  tie_break <- match.arg(tie_break)

  orig_crs <- sf::st_crs(features_sf)
  # The join runs in the polygons' CRS whenever that is projected: it is the
  # CRS the cells are drawn in, so their edges are straight and overlap areas
  # planar there.  Keeping the features' CRS instead pulled the package's own
  # projected cells into lon/lat, bending every edge into a great-circle arc,
  # and with s2 the largest-overlap join failed on a degenerate intersection
  # piece -- on sf's nc counties against a 36-cell grid, 55 of 100 counties
  # ended up outside their largest-overlap cell.
  crs_p <- sf::st_crs(polygons_sf)
  join_in_p <- !is.na(crs_p) && !isTRUE(sf::st_is_longlat(crs_p))
  hh <- harmonize_crs(features_sf, polygons_sf,
                      prefer = if (join_in_p) "b" else "a")
  f <- hh$a; p <- hh$b
  # `f` is only the join's copy.  The rows returned keep the geometry the
  # features arrived with; a CRS-less layer is resolved into the polygons'
  # CRS and returned there, as before.
  f_geom <- if (is.na(orig_crs)) sf::st_geometry(f) else sf::st_geometry(features_sf)

  id_candidates <- unique(c(polygon_id_col, "poly_id", "polygon_id", "id", "cell_id", "grid_id"))
  id_col <- id_candidates[id_candidates %in% names(p)][1]
  # Which column the IDs came from is said whenever it is not the one asked
  # for.  A fallback used to be taken without a word: nine polygons carrying a
  # zone attribute called `id` came back as three "polygons", and nine keyed
  # by `GEOID` came back numbered 1 to 9 with the GEOIDs unused.
  id_source <- "the column named by `polygon_id_col`"
  if (length(id_col) == 0L || is.na(id_col)) {
    id_col <- polygon_id_col
    p[[id_col]] <- seq_len(nrow(p))
    id_source <- "row numbers"
    .log_info(paste0("assign_features_to_polygons(): `polygons_sf` has no '%s' ",
                     "column and none of the fallback ID columns; polygons numbered ",
                     "1 to %d in row order."), polygon_id_col, nrow(p))
    .msg(sprintf(paste0(
      "assign_features_to_polygons(): `polygons_sf` has no '%s' column, and none ",
      "of the fallback ID columns (%s); the polygons are numbered 1 to %d in row ",
      "order. To use an ID the polygons already carry, name its column with ",
      "`polygon_id_col`."),
      polygon_id_col,
      paste(sprintf("'%s'", setdiff(id_candidates, polygon_id_col)), collapse = ", "),
      nrow(p)))
  } else if (!identical(id_col, polygon_id_col)) {
    id_source <- sprintf("read because `polygons_sf` has no '%s' column", polygon_id_col)
    .log_info(paste0("assign_features_to_polygons(): `polygons_sf` has no '%s' ",
                     "column; polygon IDs read from its '%s' column."),
              polygon_id_col, id_col)
    .msg(sprintf(paste0(
      "assign_features_to_polygons(): `polygons_sf` has no '%s' column; the ",
      "polygon IDs are read from its '%s' column and returned as `%s`."),
      polygon_id_col, id_col, polygon_id_col))
  }
  # One distinct, non-missing ID per polygon.  A duplicated ID merged the
  # features of two polygons under one label (a `poly_id` repeated once gave 8
  # IDs for 9 polygons), and a missing one made the features inside that
  # polygon look unassigned, so they were dropped: neither raised a condition.
  # Rows that share an ID on purpose (one unit stored as several POLYGON
  # rows) are refused with the rest, since the two cases look the same from
  # here; the message says how to keep such a unit whole, because numbering
  # the rows would split it into one cell per part.
  ids_p  <- p[[id_col]]
  id_na  <- sum(is.na(ids_p))
  id_dup <- sum(duplicated(ids_p[!is.na(ids_p)]))
  if (id_na > 0L || id_dup > 0L)
    stop(sprintf(paste0(
      "assign_features_to_polygons(): the polygon ID column '%s' (%s) holds %d ",
      "duplicated and %d missing value(s) among %d polygon(s), so features in ",
      "different polygons would share an ID or lose theirs. Give every polygon ",
      "its own ID (polygons_sf$%s <- seq_len(nrow(polygons_sf)), or ",
      "ensure_stable_poly_id()), or name the column that holds one with ",
      "`polygon_id_col`.%s"),
      id_col, id_source, id_dup, id_na, nrow(p), polygon_id_col,
      if (id_dup > 0L) sprintf(paste0(
        " If the rows that share an ID are the parts of one unit (a region ",
        "stored as several polygons), merge them into one row per ID first, so ",
        "the unit stays one cell: dplyr::summarise(dplyr::group_by(polygons_sf, ",
        "%s)). Numbering the rows would make each part a cell of its own."),
        if (identical(make.names(id_col), id_col)) id_col
        else sprintf("`%s`", id_col))
      else ""), call. = FALSE)

  # The polygons' ID travels through the join under a reserved name and is
  # renamed to `polygon_id_col` afterwards.  Under its source name, a
  # features column of that name -- a site `id` joined to cells keyed by
  # `id` -- collided inside st_join() and was dropped, although the result
  # only ever gains `polygon_id_col`.
  join_id <- "..poly_id"
  p_sel <- p[, id_col, drop = FALSE]
  names(p_sel)[names(p_sel) == id_col] <- join_id

  # st_join() suffixes columns present on both sides (`poly_id.x` /
  # `poly_id.y`), which would defeat the rename below and leave
  # `polygon_id_col` absent -- silently returning zero rows.  This is reachable
  # simply by re-assigning already-assigned points, or points that came out of
  # summarize_by_cell().  Drop the colliding column(s) up front instead.
  collide <- intersect(unique(c(polygon_id_col, join_id)), names(f))
  if (length(collide)) {
    .warn_and_log(
      "assign_features_to_polygons(): `features_sf` already carries column(s) %s, which would collide with the polygon ID column; dropping them before the spatial join. Rename them first if you need to keep them.",
      paste(sprintf("'%s'", collide), collapse = ", ")
    )
    f <- f[, setdiff(names(f), collide), drop = FALSE]
  }

  f$`..pre_join_row_id` <- seq_len(nrow(f))

  f_gtypes <- unique(as.character(sf::st_geometry_type(f, by_geometry = TRUE)))
  use_largest <- isTRUE(largest) &&
    all(f_gtypes %in% c("POLYGON", "MULTIPOLYGON"))

  if (use_largest) {
    # sf measures the overlap with st_intersection(), which GEOS cannot do on
    # an invalid ring: one bow-tie parcel threw a TopologyException, and a
    # catch-all retry without `largest` used to swallow it and reassign EVERY
    # straddling feature in the layer by `tie_break` (95 of 200 changed cell),
    # with no R warning.  An invalid ring has no well-defined overlap anyway
    # (a bow-tie's two lobes cancel to zero area), so the copies used for the
    # join are repaired, and the caller is told how many.
    bad_f <- !(sf::st_is_valid(f) %in% TRUE)
    bad_p <- !(sf::st_is_valid(p_sel) %in% TRUE)
    if (any(bad_f) || any(bad_p)) {
      .warn_and_log(paste0(
        "assign_features_to_polygons(): %d of %d feature(s) and %d of %d ",
        "polygon(s) have invalid geometry (usually a self-intersecting ring), ",
        "so their overlap areas are undefined; they were repaired with ",
        "sf::st_make_valid() for the largest-overlap join. The features are ",
        "returned with the geometry they arrived with; repair them yourself ",
        "to silence this."),
        sum(bad_f), nrow(f), sum(bad_p), nrow(p_sel))
      if (any(bad_f)) {
        g <- sf::st_geometry(f)
        g[bad_f] <- .safe_make_valid(g[bad_f])
        f <- sf::st_set_geometry(f, g)
      }
      if (any(bad_p)) {
        g <- sf::st_geometry(p_sel)
        g[bad_p] <- .safe_make_valid(g[bad_p])
        p_sel <- sf::st_set_geometry(p_sel, g)
      }
    }
  }

  join_args <- list(x = f, y = p_sel, join = predicate, left = TRUE)
  largest_ties <- integer(0)
  if (use_largest) {
    # With `largest = TRUE` sf never calls `predicate` (it intersects the
    # layers and keeps the biggest piece), so no predicate can reject it: an
    # error here is a geometry failure.  Falling back to `tie_break` would
    # change the rule for every feature in the layer, so stop instead.
    #
    # sf's st_intersection() warns "attribute variables are assumed to be
    # spatially constant throughout all geometries" unless every attribute
    # is flagged constant, and st_join() adds its own unflagged grouping
    # columns before calling it, so every ordinary call raised it (setting
    # st_agr() on the inputs does not help).  Keeping one cell per feature is
    # exactly that assumption, so the warning says nothing about the data;
    # it alone is muffled.
    joined <- tryCatch(
      withCallingHandlers(
        do.call(sf::st_join, c(join_args, largest = TRUE)),
        warning = function(w) {
          if (grepl("attribute variables are assumed to be spatially constant",
                    conditionMessage(w), fixed = TRUE))
            invokeRestart("muffleWarning")
        }),
      error = function(e) stop(sprintf(paste0(
        "assign_features_to_polygons(): the largest-overlap join failed (%s). ",
        "Pass largest = FALSE to assign every feature by `predicate` and ",
        "`tie_break` instead%s."),
        conditionMessage(e),
        if (isTRUE(sf::st_is_longlat(f)))
          ", or project both layers first (see ensure_projected())" else ""),
        call. = FALSE))
    # sf keeps the largest intersection piece without asking whether it has
    # any area, so under GEOS a feature that only shares an edge or a corner
    # with the cells was assigned to one of them with zero overlap, while
    # under s2 (lon/lat, s2 on) the same feature had no piece and came back
    # unassigned: sf's 100 nc counties against 50 of them as cells gave 70
    # rows projected and 50 in lon/lat.  A feature with no overlap has no largest
    # overlap, so it is unassigned whatever the CRS.  s2 already does this.
    if (!(isTRUE(sf::st_is_longlat(f)) && isTRUE(sf::sf_use_s2()))) {
      hit <- which(!is.na(joined[[join_id]]))
      if (length(hit)) {
        # "2********": the interiors meet in an area.
        overlaps <- suppressMessages(sf::st_relate(f, p_sel, pattern = "2********"))
        ids_p    <- p_sel[[join_id]]
        fi       <- joined[["..pre_join_row_id"]][hit]
        touch    <- !vapply(seq_along(hit), function(k)
          joined[[join_id]][hit[k]] %in% ids_p[overlaps[[fi[k]]]], logical(1))
        if (any(touch)) {
          .log_info(paste0("assign_features_to_polygons(): %d polygon feature(s) ",
                           "only touch the polygon layer (no overlap area) and ",
                           "are left unassigned."), sum(touch))
          joined[[join_id]][hit[touch]] <- NA
        }
      }
    }
    # sf keeps which.max() of the overlap areas, so a feature split evenly
    # between two cells went to whichever comes first in `polygons_sf`, and
    # no tie was recorded: a 40 m square across the edge of cells 1 and 2
    # went to cell 1, or to cell 2 with the polygon rows reversed.  Features
    # overlapping two or more cells are measured again, and an exact tie
    # (areas equal to 9 significant digits, as for the cell areas below) is
    # decided by `tie_break` among the equally largest cells.
    lt <- .largest_overlap_ties(f, p_sel, joined, join_id, tie_break)
    if (length(lt$rows)) {
      joined[[join_id]][lt$rows] <- lt$ids
      largest_ties <- joined[["..pre_join_row_id"]][lt$rows]
      .log_info(paste0("assign_features_to_polygons(): %d of %d polygon ",
                       "feature(s) overlap two or more polygons by exactly the ",
                       "same area; the '%s' rule chose one for each."),
                length(lt$rows), nrow(f), tie_break)
    }
  } else {
    joined <- do.call(sf::st_join, join_args)
  }
  # The join ran on copies (moved into the polygons' CRS, repaired); hand
  # back the geometry the caller passed, in the CRS it arrived in, rather
  # than a transform round trip of it.
  joined <- sf::st_set_geometry(joined, f_geom[joined[["..pre_join_row_id"]]])

  if (!identical(join_id, polygon_id_col)) {
    names(joined)[names(joined) == join_id] <- polygon_id_col
  }
  if (!polygon_id_col %in% names(joined)) {
    stop(sprintf(
      paste0("assign_features_to_polygons(): the spatial join did not produce ",
             "the expected '%s' column (join result has: %s). This usually ",
             "means a column-name collision between `features_sf` and ",
             "`polygons_sf`; rename the offending column in `features_sf`."),
      polygon_id_col, paste(names(joined), collapse = ", ")
    ), call. = FALSE)
  }


  # Deterministic tie-break for features matching multiple polygons.
  # The previous approach (keep first duplicate) was order-dependent.
  dup_mask <- duplicated(joined[["..pre_join_row_id"]])
  # How many features the tie-break decided, and which: this used to be a
  # branch condition and nothing else, yet a tie-break firing on a third of
  # the features means the polygon layer overlaps and every cell count
  # downstream is suspect.  Reported on the result as `ties`, and logged.
  tie_rows <- sort(unique(joined[["..pre_join_row_id"]][dup_mask]))
  if (length(tie_rows))
    .log_info(paste0("assign_features_to_polygons(): %d of %d feature(s) fall ",
                     "inside more than one polygon; the '%s' rule chose one for ",
                     "each. A large share means the polygon layer overlaps."),
              length(tie_rows), nrow(f), tie_break)
  if (any(dup_mask) && identical(tie_break, "smallest_area")) {
    # For each duplicated feature, keep the polygon with the smallest area
    # (the most specific / tightest-fitting polygon).  Areas are compared to
    # 9 significant digits, so equal cells whose computed areas differ in the
    # last bits count as equal.
    poly_areas <- suppressWarnings(as.numeric(sf::st_area(p)))
    poly_areas[!is.finite(poly_areas)] <- Inf
    names(poly_areas) <- as.character(p[[id_col]])

    joined$`..poly_area` <- signif(poly_areas[as.character(joined[[polygon_id_col]])], 9L)
    joined$`..poly_area`[is.na(joined$`..poly_area`)] <- Inf

    # Equal areas -- every cell of a regular grid, for a point on a shared
    # edge -- used to fall through to row order, so the answer depended on
    # the order of the polygon layer after all: edge points at x = 100 went
    # to cells 1, 4, 7, or to 2, 5, 8 with the rows reversed.  The candidate
    # whose bounding-box centre is lowest, then leftmost, wins instead.  On
    # a create_grid_polygons() grid, which is numbered from the lower left a
    # row at a time, that is the cell the row order picked.
    tie_ids <- unique(as.character(
      joined[[polygon_id_col]][joined[["..pre_join_row_id"]] %in% tie_rows]))
    tie_ids <- tie_ids[!is.na(tie_ids)]
    p_tie   <- sf::st_geometry(p)[match(tie_ids, as.character(p[[id_col]]))]
    ctr     <- vapply(p_tie, function(g) {
      bb <- sf::st_bbox(g)
      c((bb[["xmin"]] + bb[["xmax"]]) / 2, (bb[["ymin"]] + bb[["ymax"]]) / 2)
    }, numeric(2))
    ctr_x <- stats::setNames(ctr[1, ], tie_ids)
    ctr_y <- stats::setNames(ctr[2, ], tie_ids)
    ids_j <- as.character(joined[[polygon_id_col]])
    key_y <- unname(ctr_y[ids_j]); key_y[!is.finite(key_y)] <- Inf
    key_x <- unname(ctr_x[ids_j]); key_x[!is.finite(key_x)] <- Inf

    # Within each group of duplicates, keep the row with smallest area
    joined <- joined[order(joined[["..pre_join_row_id"]], joined[["..poly_area"]],
                           key_y, key_x), , drop = FALSE]
    joined <- joined[!duplicated(joined[["..pre_join_row_id"]]), , drop = FALSE]
    joined[["..poly_area"]] <- NULL
  } else {
    joined <- joined[!duplicated(joined[["..pre_join_row_id"]]), , drop = FALSE]
  }
  # Which features matched no polygon, by their row in `features_sf`, read
  # before the row key goes.
  unassigned_rows <- sort(joined[["..pre_join_row_id"]][is.na(joined[[polygon_id_col]])])
  joined[["..pre_join_row_id"]] <- NULL

  n_unassigned <- length(unassigned_rows)
  n_features   <- nrow(joined)
  # Zero matches is almost never what the caller meant, and it is silent: the
  # function returns an empty sf and every downstream step reports zero rows
  # without saying why.  The usual cause is a CRS mismatch that could only be
  # stamped rather than reprojected, which puts the two layers in different
  # places.
  if (n_unassigned == n_features && n_features > 0L)
    .warn_classed("spatialkit_unassigned",
                  paste0("assign_features_to_polygons(): none of the %d ",
                         "feature(s) fall inside any polygon. Check that the ",
                         "two layers cover the same ground -- a CRS that had to ",
                         "be stamped rather than reprojected is the usual ",
                         "cause -- or pass a different `predicate`."),
                  n_features, data = list(n = n_unassigned, total = n_features))
  else if (n_unassigned > 0L && !keep_unassigned)
    # Rows removed from the result are data dropped, which the caller has to
    # be able to catch: this was one INFO line in the temp log, and a
    # tessellation of 25 k-means seeds built without a boundary lost 764 of
    # 3,000 observations (the edge of the study area, so the cell means moved
    # too) with no condition raised and nothing on the result.
    .warn_classed("spatialkit_unassigned",
                  paste0("assign_features_to_polygons(): %d of %d feature(s) ",
                         "(%.1f%%) fall inside no polygon and are left out of ",
                         "the result; attr(, \"unassigned\")$which gives their ",
                         "rows in `features_sf`. Check that the polygons cover ",
                         "the features (cells built from seeds need a ",
                         "`boundary` or a `cover`; see build_tessellation()), or ",
                         "pass keep_unassigned = TRUE to keep them with an NA ID."),
                  n_unassigned, n_features, 100 * n_unassigned / n_features,
                  data = list(n = n_unassigned, total = n_features))
  else if (n_unassigned > 0L)
    .log_info("assign_features_to_polygons(): %d of %d feature(s) fall inside no polygon and are kept with an NA ID.",
              n_unassigned, n_features)

  if (!keep_unassigned) {
    joined <- joined[!is.na(joined[[polygon_id_col]]), , drop = FALSE]
  }

  # Stamped and classed: the record names row positions, so it must not
  # survive a subset that renumbers or removes them.
  tie_rows <- sort(unique(c(tie_rows, largest_ties)))
  joined <- .set_row_record(joined, "ties",
                            list(n = length(tie_rows),
                                 which = as.integer(tie_rows),
                                 rule = tie_break))
  # The features no polygon holds, beside the ties: the count used to exist
  # only in a log line.  Stamped the same way, and removed by a subset with
  # the other records named in .row_record_attrs (R/utils.R).
  joined <- .set_row_record(joined, "unassigned",
                            list(n = n_unassigned,
                                 which = as.integer(unassigned_rows)))
  joined
}


# Exact ties in a largest-overlap join.  `joined` is st_join(largest = TRUE)
# of `f` onto `p_sel`, the polygons' ID in column `join_id`.  The features
# assigned there that overlap two or more polygons are intersected with them
# again; where the largest overlap is shared (areas equal to 9 significant
# digits) the polygon is picked among the equally largest by `tie_break`:
# "smallest_area" takes the smallest polygon, then the one whose bounding-box
# centre is lowest, then leftmost, as for the predicate join; "first" takes
# the first by polygon row.  Returns the rows of `joined` that had a tie and
# the ID each gets.
.largest_overlap_ties <- function(f, p_sel, joined, join_id, tie_break) {
  none <- list(rows = integer(0), ids = NULL)
  rows <- which(!is.na(joined[[join_id]]))
  if (!length(rows)) return(none)
  gf <- sf::st_geometry(f)[joined[["..pre_join_row_id"]][rows]]
  gp <- sf::st_geometry(p_sel)
  cand  <- suppressMessages(sf::st_intersects(gf, gp))
  multi <- which(lengths(cand) >= 2L)
  if (!length(multi)) return(none)
  pieces <- tryCatch(suppressMessages(suppressWarnings(
    sf::st_intersection(gf[multi], gp))), error = function(e) NULL)
  if (is.null(pieces) || !length(pieces)) return(none)
  idx  <- attr(pieces, "idx")
  area <- signif(suppressWarnings(as.numeric(sf::st_area(pieces))), 9L)
  area[!is.finite(area)] <- 0
  out_rows <- integer(0)
  out_ids  <- p_sel[[join_id]][integer(0)]
  for (s in split(seq_along(area), idx[, 1L])) {
    top <- max(area[s])
    if (!(top > 0)) next
    tied <- sort(unique(idx[s[area[s] == top], 2L]))
    if (length(tied) < 2L) next
    pick <- if (identical(tie_break, "first")) tied[1L] else {
      ga <- signif(suppressWarnings(as.numeric(sf::st_area(gp[tied]))), 9L)
      ga[!is.finite(ga)] <- Inf
      ctr <- vapply(gp[tied], function(g) {
        bb <- sf::st_bbox(g)
        c((bb[["xmin"]] + bb[["xmax"]]) / 2, (bb[["ymin"]] + bb[["ymax"]]) / 2)
      }, numeric(2))
      tied[order(ga, ctr[2, ], ctr[1, ])[1L]]
    }
    out_rows <- c(out_rows, rows[multi[idx[s[1L], 1L]]])
    out_ids  <- c(out_ids, p_sel[[join_id]][pick])
  }
  list(rows = out_rows, ids = out_ids)
}


#' Correlation implied by a fitted gstat variogram model
#'
#' Converts a fitted variogram to a correlation function of distance, giving
#' the correlation between two \strong{distinct} observations a distance
#' \code{h} apart.  For a nested model with nugget \eqn{c_0} and structured
#' components \eqn{(c_i, a_i)} the semivariance is
#' \eqn{\gamma(h) = c_0 + \sum_i c_i (1 - f_i(h))}, so the correlation is
#' \eqn{1 - \gamma(h) / (c_0 + \sum_i c_i) = \sum_i c_i f_i(h) / (c_0 + \sum_i c_i)}:
#' every structured component contributes, weighted by its partial sill.  The
#' function used to read the single largest component only, which for a
#' user-built \code{Nug + Exp + Sph} model gave 0.108 at \code{h = 200} where
#' \code{gstat::variogramLine()} implies 0.197.
#'
#' Only the exponential, spherical and Gaussian families are implemented.  A
#' model carrying any other family (Matern, power, circular, ...) returns
#' \code{NULL} and is never silently read as exponential, so the caller
#' falls back to \code{deff = 1} and says so.  \code{\link{estimate_sac_range}}
#' only ever produces single-component \code{Exp} or \code{Sph} models; other
#' shapes reach this function through a user-built \code{sac}.
#'
#' A component with a 2-D geometric anisotropy (gstat's \code{ang1} and
#' \code{anis1}, from \code{vgm(..., anis = c(angle, ratio))}) cannot be
#' evaluated at a scalar distance, so the function of distance reads it at
#' its major range, and the returned function carries an attribute
#' \code{"cor_xy"}: a function of a coordinate matrix returning the full
#' correlation matrix, with each component's separations rotated and its
#' minor axis stretched by \code{1 / anis1} as gstat does.  Callers that hold
#' coordinates use that.  The 3-D terms (\code{ang2}, \code{ang3},
#' \code{anis2}) do not act on 2-D coordinates and are ignored.
#'
#' There is deliberately no special case at \code{h = 0}.  The nugget captures
#' measurement error and variation below the sampling resolution, so two
#' distinct observations each carry independent nugget noise and are less than
#' perfectly correlated even when they coincide.  Self-correlation of 1 belongs
#' on the diagonal of the correlation matrix and is imposed by the caller.
#'
#' @param vgm_model A fitted \code{gstat} variogram model (a data frame with
#'   \code{model}, \code{psill} and \code{range} columns).
#' @return A function of distance returning correlation, or \code{NULL} if the
#'   model cannot be interpreted or uses a family this function does not
#'   implement.
#' @keywords internal
#' @noRd
.vgm_correlation_fn <- function(vgm_model) {
  if (is.null(vgm_model) || !is.data.frame(vgm_model)) return(NULL)
  if (!all(c("model", "psill", "range") %in% names(vgm_model))) return(NULL)

  fam    <- as.character(vgm_model$model)
  is_nug <- fam == "Nug"
  nugget <- sum(vgm_model$psill[is_nug], na.rm = TRUE)
  struct <- vgm_model[!is_nug, , drop = FALSE]
  if (nrow(struct) == 0L) return(NULL)
  if (!all(as.character(struct$model) %in% names(.vgm_shape_fns))) return(NULL)

  psill <- sum(struct$psill, na.rm = TRUE)
  total <- nugget + psill
  if (!is.finite(total) || total <= 0 || !is.finite(psill) || psill <= 0)
    return(NULL)
  c_i <- as.numeric(struct$psill)
  a_i <- as.numeric(struct$range)
  if (any(!is.finite(c_i) | c_i < 0) || any(!is.finite(a_i) | a_i <= 0))
    return(NULL)
  f_i <- .vgm_shape_fns[as.character(struct$model)]

  fn <- function(h) {
    num <- 0
    for (i in seq_along(c_i)) num <- num + c_i[i] * f_i[[i]](h, a_i[i])
    # No special case at h = 0: this is the correlation between two DISTINCT
    # observations, which the nugget discounts even when they coincide.  An
    # observation with itself correlates at 1, and the caller imposes that on
    # the diagonal.
    pmin(pmax(num / total, 0), 1)
  }

  # A 2-D anisotropy was read as isotropic at the major range: for
  # vgm(0.8, "Exp", 300, 0.2, anis = c(0, 0.2)) the correlation at 50 m
  # east-west came out 0.677 where gstat's is 0.348, and the median design
  # effect 14.1 against 7.7.  Distance alone cannot carry a direction, so
  # coordinate-holding callers get a function of the coordinates.
  ang  <- if ("ang1" %in% names(struct)) as.numeric(struct$ang1) else rep(0, nrow(struct))
  anis <- if ("anis1" %in% names(struct)) as.numeric(struct$anis1) else rep(1, nrow(struct))
  ang[!is.finite(ang)] <- 0
  anis[!is.finite(anis) | anis <= 0 | anis > 1] <- 1
  if (any(anis < 1)) {
    attr(fn, "cor_xy") <- function(xy) {
      num <- 0
      for (i in seq_along(c_i)) {
        # gstat: ang1 is the major axis's direction, clockwise from north;
        # the separation along the minor axis is stretched by 1 / anis1.
        th <- ang[i] * pi / 180
        u  <- xy[, 1] * sin(th) + xy[, 2] * cos(th)
        v  <- (xy[, 1] * cos(th) - xy[, 2] * sin(th)) / anis[i]
        d  <- as.matrix(stats::dist(cbind(u, v)))
        num <- num + c_i[i] * f_i[[i]](d, a_i[i])
      }
      pmin(pmax(num / total, 0), 1)
    }
  }
  fn
}

# Correlation shape f(h; a) of each supported gstat family, so that the
# semivariance of a component with partial sill c is c * (1 - f(h; a)).
# These are gstat's own parametrisations (variogramLine() reproduces them):
# the exponential's `a` is the range PARAMETER, reaching 95% of the sill at
# 3a; the spherical reaches its sill exactly at a; the Gaussian at ~1.73a.
.vgm_shape_fns <- list(
  Exp = function(h, a) exp(-h / a),
  Sph = function(h, a) ifelse(h >= a, 0, 1 - 1.5 * (h / a) + 0.5 * (h / a)^3),
  Gau = function(h, a) exp(-(h / a)^2)
)

# TRUE for a variogram model with a positive nugget and no structured sill
# (no structured component, or only ones of partial sill 0, of the families
# above): it implies no correlation between distinct observations.
# .vgm_correlation_fn() returns NULL for it, deliberately, since it has no
# correlation function to offer; callers that can use "uncorrelated" ask this.
.vgm_is_pure_nugget <- function(vgm_model) {
  if (!is.data.frame(vgm_model) || !all(c("model", "psill") %in% names(vgm_model)))
    return(FALSE)
  fam <- as.character(vgm_model$model)
  ps  <- suppressWarnings(as.numeric(vgm_model$psill))
  if (!length(fam) || !all(fam %in% c("Nug", names(.vgm_shape_fns))) ||
      !all(is.finite(ps)))
    return(FALSE)
  sum(ps[fam == "Nug"]) > 0 && all(ps[fam != "Nug"] == 0)
}


#' Standard error of a cell mean under a design effect
#'
#' Two corrections, not one.  A design effect inflates the variance of the mean
#' to \eqn{\sigma^2 \mathrm{deff} / n}, which is what \code{s/sqrt(n/deff)}
#' applies.  But under the same within-cell correlation the ordinary sample
#' variance is \emph{also} biased downward.  For exchangeable correlation
#' \eqn{\rho} (Kish's own assumption), with \eqn{\mathrm{deff} = 1 + (n-1)\rho}:
#' \deqn{E[s^2] = \sigma^2 (n - \mathrm{deff}) / (n - 1)}
#' so \eqn{s^2} understates \eqn{\sigma^2} by very nearly the factor deff
#' inflates the mean's variance by, and the two errors compound rather than
#' cancel.  Measured 95% CI coverage of the uncorrected form at \eqn{n = 30}
#' over 20,000 replicates: 0.921 at \eqn{\rho = 0.2}, 0.844 at \eqn{\rho = 0.5},
#' 0.628 at \eqn{\rho = 0.8}.  With the \eqn{\sqrt{(n-1)/(n-\mathrm{deff})}}
#' rescale below: 0.948, 0.950, 0.949.
#'
#' At \code{deff = 1} the factor is exactly 1, so the default path is
#' bit-identical to \code{s/sqrt(n)}.
#'
#' @param s Within-cell standard deviation.
#' @param n Number of non-missing observations in the cell.
#' @param deff Design effect for that cell, \eqn{\geq 1}.
#' @return The standard error, or \code{NA_real_} when the cell carries no
#'   information about \eqn{\sigma} (\code{deff >= n}: complete redundancy).
#' @keywords internal
#' @noRd
.se_with_deff <- function(s, n, deff, correct_s2 = TRUE) {
  if (!is.finite(s) || !is.finite(n) || n <= 1L) return(NA_real_)
  if (!is.finite(deff) || deff <= 1) return(s / sqrt(n))
  # correct_s2 = FALSE for a design effect the CALLER supplied as a constant.
  # The E[s^2] = sigma^2 (n - deff)/(n - 1) correction below is derived from a
  # within-cell correlation structure, so it is right for an ICC- or
  # variogram-derived deff and unjustified for an externally chosen number:
  # applied to deff = 2 it doubled a 3-point cell's SE instead of multiplying
  # it by sqrt(2), and returned NA for every cell with n <= deff -- neither of
  # which is the "uniform inflation" the argument is documented as.
  if (!isTRUE(correct_s2)) return(s * sqrt(deff / n))
  # deff is bounded above by n (every observation a copy of every other).  At
  # that bound the cell holds one observation's worth of information and s
  # carries none about sigma, so there is no standard error to report.
  if (deff >= n) return(NA_real_)
  s * sqrt(deff / n) * sqrt((n - 1) / (n - deff))
}


#' Per-cell design effect from a fitted variogram
#'
#' For \code{n} observations in a cell with correlation matrix \code{R}, the
#' effective sample size of the mean is \code{n^2 / sum(R)}, so the design
#' effect is \code{sum(R) / n}.  This generalises Kish's
#' \code{1 + (n - 1) * rho}, which is the special case of a constant
#' off-diagonal correlation: substituting \code{R = I + rho(J - I)} gives
#' \code{sum(R) / n = 1 + (n - 1) * rho} exactly.  Using the variogram instead
#' lets correlation decay with distance, which is the whole point of having
#' fitted one.  Kish assumes every pair in a cell is equally correlated
#' regardless of how far apart they are, and that assumption degrades as cells
#' get larger.
#'
#' @param coords Numeric matrix of coordinates.
#' @param cell_id Vector of cell identifiers, same length as \code{nrow(coords)}.
#' @param cor_fn Correlation function from \code{.vgm_correlation_fn()}.
#' @param max_n Cells larger than this are subsampled before forming the
#'   \code{n x n} correlation matrix.  Default 500.
#' @param seed RNG seed for that subsampling.
#' @param r_bar The per-cell mean correlations, when the caller already has
#'   them from \code{.cell_cor_stats_variogram()}: building every cell's
#'   correlation matrix a second time only to recompute them doubled the
#'   cost of \code{summarize_by_cell(deff = "variogram")}.
#' @return Named numeric vector of design effects, one per cell.
#' @keywords internal
#' @noRd
.cell_deff_variogram <- function(coords, cell_id, cor_fn, max_n = 500L,
                                 seed = 42L, r_bar = NULL) {
  if (is.null(r_bar))
    r_bar <- .cell_rbar_variogram(coords, cell_id, cor_fn, max_n = max_n,
                                  seed = seed)
  n_i   <- table(factor(cell_id, levels = names(r_bar)))
  n_i   <- as.numeric(n_i[names(r_bar)])
  out   <- 1 + (n_i - 1) * r_bar
  # Bounded below by 1 (independence) and above by the cell size (complete
  # redundancy); a NULL cor_fn leaves NA throughout.
  out[!is.na(out)] <- pmin(pmax(out[!is.na(out)], 1), n_i[!is.na(out)])
  out[n_i <= 1L & !is.na(r_bar)] <- 1
  out
}


#' Mean off-diagonal correlation within each cell, from a fitted variogram
#'
#' The quantity a design effect is built from.  For a cell of \eqn{n} points
#' with correlation matrix \eqn{R}, \eqn{\mathrm{deff} = \sum R / n = 1 + (n-1)
#' \bar{r}} where \eqn{\bar{r}} is the mean off-diagonal correlation.  It is
#' \eqn{\bar{r}}, not deff, that is returned, because the design effect a
#' \emph{column} needs is the one for the observations that column actually
#' has: a cell of 10 points with 2 non-missing responses has a response design
#' effect of \eqn{1 + (2-1)\bar{r}}, not \eqn{1 + (10-1)\bar{r}}.  Measured:
#' adding 8 NA-response rows to a 2-observation cell moved its response SE from
#' 8.46 to 26.38 when the cell-level deff was applied to it; NA rows carry no
#' information and must not move anything.
#'
#' Subsampling estimates \eqn{\bar{r}} just as well as the full cell does
#' (the subsample's pairwise-distance distribution is the cell's), which is
#' why the per-cell quantity kept is this and not \code{sum(R)/n_used}.
#'
#' @param coords Numeric matrix of \strong{projected} coordinates, in the CRS
#'   the variogram was fitted in.  A variogram range is a length in that CRS;
#'   evaluating it at degree distances saturates every pair.
#' @param cell_id Vector of cell identifiers, same length as \code{nrow(coords)}.
#' @param cor_fn Correlation function from \code{.vgm_correlation_fn()}.
#' @param max_n Cells larger than this are subsampled before forming the
#'   correlation matrix.  Default 500.
#' @param seed RNG seed for that subsampling.
#' @return Named numeric vector of mean off-diagonal correlations, one per
#'   cell; \code{NA} when \code{cor_fn} is \code{NULL}, \code{0} for a
#'   single-point cell (no pairs).
#' @keywords internal
#' @noRd
.cell_rbar_variogram <- function(coords, cell_id, cor_fn, max_n = 500L,
                                 seed = 42L) {
  .cell_cor_stats_variogram(coords, cell_id, cor_fn, max_n = max_n,
                            seed = seed)$rbar
}


#' Per-cell correlation statistics from a fitted variogram
#'
#' One pass over the cells building each correlation matrix once, returning
#' both quantities \code{summarize_by_cell()} takes from it: the mean
#' off-diagonal correlation (\code{rbar}, see \code{.cell_rbar_variogram()})
#' and the degrees-of-freedom ratio of the within-cell variance
#' (\code{df_ratio}, see \code{.cor_stats_from_coords()}).
#'
#' @return A list of two named numeric vectors, \code{rbar} and
#'   \code{df_ratio}, one entry per cell; both \code{NA} when \code{cor_fn}
#'   is \code{NULL}.
#' @keywords internal
#' @noRd
.cell_cor_stats_variogram <- function(coords, cell_id, cor_fn, max_n = 500L,
                                      seed = 42L) {
  ids <- unique(cell_id)
  ids <- ids[!is.na(ids)]
  rbar <- stats::setNames(rep(NA_real_, length(ids)), as.character(ids))
  dfr  <- rbar
  if (is.null(cor_fn)) return(list(rbar = rbar, df_ratio = dfr))
  for (id in ids) {
    idx <- which(cell_id == id)
    st  <- .cor_stats_from_coords(coords[idx, , drop = FALSE], cor_fn,
                                  max_n = max_n, seed = seed)
    rbar[[as.character(id)]] <- st[["rbar"]]
    dfr[[as.character(id)]]  <- st[["df_ratio"]]
  }
  list(rbar = rbar, df_ratio = dfr)
}


#' Mean off-diagonal correlation among a set of locations
#'
#' The single-cell core of \code{.cell_rbar_variogram()}, exposed so the
#' \code{..se_} closures can call it on exactly the rows a column is
#' non-missing for.  Subsampling (above \code{max_n}) is seeded and the RNG is
#' restored, so repeated calls from inside \code{dplyr::summarise()} neither
#' drift nor disturb the caller's stream.
#'
#' @return A scalar in \eqn{[0, 1]}; \code{0} when there is at most one point.
#' @keywords internal
#' @noRd
.rbar_from_coords <- function(xy, cor_fn, max_n = 500L, seed = 42L) {
  unname(.cor_stats_from_coords(xy, cor_fn, max_n = max_n, seed = seed)[["rbar"]])
}


#' Mean off-diagonal correlation and variance degrees of freedom for one cell
#'
#' Builds the correlation matrix \eqn{R} of a set of locations once and takes
#' two numbers from it.  \code{rbar} is the mean off-diagonal correlation that
#' \code{.rbar_from_coords()} returns.  \code{df_ratio} is the Satterthwaite
#' degrees of freedom of the within-cell sample variance \eqn{s^2} under
#' \eqn{R}, divided by \eqn{n - 1}: with \eqn{A = I - J/n} the centring
#' projection, \eqn{s^2 = y'Ay/(n-1)} is a weighted sum of chi-squares whose
#' moment-matched df is \eqn{[\mathrm{tr}(AR)]^2 / \mathrm{tr}(ARAR)}.  Under
#' exchangeable correlation the ratio is exactly 1 (\eqn{AR} has \eqn{n-1}
#' equal non-zero eigenvalues), which is why a t interval on the Kish path
#' keeps \eqn{n - 1} df; under correlation that decays with distance it is
#' below 1, and a t interval that ignores that under-covers (measured 0.918 at
#' a range of 400 on a 1000-unit domain against 0.958 with this df).
#' \eqn{\mathrm{tr}(AR) = n - \mathrm{deff}}, so the same \eqn{R} gives the
#' design effect and the df with no extra pass.
#'
#' Returning the ratio instead of the df lets a subsampled cell
#' (\code{n > max_n}) scale it back to its own \eqn{n - 1}: the subsample's
#' pairwise-distance distribution is the cell's, so the ratio transfers where
#' the raw df would not.
#'
#' Rows with a non-finite coordinate (an empty point) have no location to
#' correlate, so they are left out: the statistics describe the located
#' points.  An anisotropic model is evaluated through its \code{"cor_xy"}
#' attribute (see \code{.vgm_correlation_fn()}).
#'
#' @return Named numeric vector \code{c(rbar =, df_ratio =)}; \code{c(0, 1)}
#'   when there is at most one located point or no correlation function.
#' @keywords internal
#' @noRd
.cor_stats_from_coords <- function(xy, cor_fn, max_n = 500L, seed = 42L) {
  # An empty point's NA coordinates put NA in R and stopped the df test
  # below with "missing value where TRUE/FALSE needed".
  if (!is.null(xy) && length(dim(xy)) == 2L && nrow(xy) > 0L) {
    located <- is.finite(xy[, 1]) & is.finite(xy[, 2])
    if (!all(located)) xy <- xy[located, , drop = FALSE]
  }
  n_i <- nrow(xy)
  if (is.null(cor_fn) || is.null(n_i) || n_i <= 1L)
    return(c(rbar = 0, df_ratio = 1))
  if (n_i > max_n) {
    cleanup <- .with_seed(seed)
    on.exit(cleanup(), add = TRUE)
    xy <- xy[sample.int(n_i, max_n), , drop = FALSE]
  }
  # Rebuild explicitly rather than relying on cor_fn() to preserve `dim`.
  # A correlation function written as, say, rep(1, length(h)) returns a bare
  # vector, and diag()<- would then fail.
  cor_xy <- attr(cor_fn, "cor_xy", exact = TRUE)
  R <- if (is.function(cor_xy)) {
    matrix(as.numeric(cor_xy(xy)), nrow = nrow(xy), ncol = nrow(xy))
  } else {
    d <- as.matrix(stats::dist(xy))
    matrix(as.numeric(cor_fn(as.numeric(d))), nrow = nrow(d), ncol = ncol(d))
  }
  diag(R) <- 1
  n_used <- nrow(R)
  r_bar  <- (sum(R) - n_used) / (n_used * (n_used - 1))
  # tr(AR) = n - sum(R)/n and, with c = colMeans(R) (R is symmetric, so row
  # means are the same), tr(ARAR) = sum(R^2) - 2 n sum(c^2) + (sum(c))^2.
  cm     <- colMeans(R)
  tr_ar  <- n_used - sum(R) / n_used
  tr_ar2 <- sum(R * R) - 2 * n_used * sum(cm^2) + sum(cm)^2
  df_r   <- if (tr_ar2 > 0) (tr_ar^2 / tr_ar2) / (n_used - 1) else 0
  c(rbar = min(max(r_bar, 0), 1), df_ratio = min(max(df_r, 0), 1))
}


#' Warn that a requested design effect was not the one applied
#'
#' The classed warning for a design effect that fell back: to 1 for a refused
#' number, and to \code{"kish"} (or, when its ICCs are 0, to nothing) for a
#' \code{"variogram"} request with no usable model.  The class is
#' \code{"spatialkit_deff_fallback"} in both cases.
#'
#' @keywords internal
#' @noRd
.warn_deff_fallback <- function(fmt, ..., data = list()) {
  .warn_classed("spatialkit_deff_fallback", fmt, ..., data = data)
}


#' Summarize features by polygon/cell ID
#'
#' Aggregates an sf point dataset into one row per cell. By default computes
#' counts and means, but the aggregation function is configurable.
#'
#' This is the third step of the package's pipeline, taking the labelled layer
#' from [assign_features_to_polygons()] down to cell level. What distinguishes
#' it from a plain `dplyr::group_by()` + `summarise()` is that it carries the
#' *uncertainty* of each aggregate with it: alongside every mean it returns a
#' within-cell standard deviation, a standard error and an observation count,
#' and it can correct that standard error for within-cell spatial
#' autocorrelation via `deff`. Reach for it whenever the cell-level values will
#' be modelled or mapped, because a cell mean over 2 observations and one over
#' 200 are not the same measurement and nothing downstream can tell them apart
#' otherwise.
#'
#' In addition to user-specified aggregation functions, this function always
#' computes within-cell standard deviation (`..sd_<var>`) and standard error
#' (`..se_<var>`) for every numeric response/predictor column, plus an `n`
#' column (rows falling in the cell) and a `cell_weight` column.
#' These columns let downstream models account for the fact that a cell with
#' 2 observations carries more aggregation uncertainty than one with 200.
#'
#' `cell_weight` is the *effective* sample size of the primary variable: the
#' response when one was supplied, otherwise the first predictor. It counts
#' that variable's non-missing rows, not all rows (a cell of 10 rows with 3
#' finite responses carries 3 observations' worth of information about the
#' response, not 10), and it is divided by that cell's design effect when
#' `deff` applied one. With `deff = 1` and no missing values it equals `n`.
#' Pass it as the `weights` argument of a downstream regression.
#'
#' An infinite value (`Inf`, `-Inf`) in a response or predictor column is
#' treated as missing, as `NA` and `NaN` are, in every statistic of that
#' column: the `agg_funs` results, `..sd_*`, `..se_*`, the interval columns,
#' `..neff_*`, `cell_weight` and the design effect. `n` still counts the row.
#' A warning of class `"spatialkit_nonfinite"` names the columns and gives the
#' number of rows and cells affected (the condition carries `n` and
#' `columns`). [kriging_adequacy()], [resolution_profile()],
#' [estimate_sac_range()] and [prep_model_data()] leave such rows out too, so
#' the steps of one pipeline agree on a cell's mean. A function in `agg_funs`
#' therefore never sees an infinite value: `max` returns the largest finite
#' value, and infinities cannot be counted through `agg_funs`.
#'
#' @section Spatial autocorrelation and standard-error bias:
#' By default (`deff = 1`), the `..se_*` columns are computed as
#' `sd / sqrt(n)`, which treats the observations within each cell as
#' independent. When data are spatially autocorrelated (the common case for the
#' spatial workflows this package supports), within-cell observations are
#' typically positively correlated: they share the cell's departure from the
#' population mean, so as an estimate of the **population (grand) mean** a cell
#' mean has an effective sample size smaller than `n`. For that estimand the
#' naive SE is **anticonservative** (too small), and a downstream weighted
#' regression that uses `cell_weight` or the `..se_*` columns for population-level
#' inference will produce overconfident standard errors for cells with strong
#' intra-cell correlation. For the cell's **own** mean the naive SE is the right
#' one when the cell's points are spread through it, and the corrected SE is too
#' wide; see "What the standard error estimates" before setting `deff`.
#'
#' Setting `deff = "kish"` applies an approximate correction using Kish's
#' design effect. An intra-class correlation (ICC) is estimated for every
#' summarised column on its own, the response and each predictor, via a
#' one-way random-effects decomposition across all cells. Each column's ICC
#' is used for that column's SE adjustment, and each cell's effective sample
#' size is reduced to `n_i / (1 + (n_i - 1) * rho)`. This is a first-order
#' correction that does not require a full spatial covariance model but does
#' require enough cells and observations for a stable ICC estimate.
#'
#' A design effect estimated from the data (`"kish"` or `"variogram"`) comes
#' with a second, small-sample correction. The within-cell correlation that
#' inflates the variance of the mean to `sigma^2 * deff / n` also biases the
#' within-cell sample variance downward: under exchangeable correlation `rho`
#' (Kish's own assumption), with `deff = 1 + (n - 1) * rho`,
#' `E[s^2] = sigma^2 * (n - deff) / (n - 1)`, so `s^2` understates `sigma^2`
#' by very nearly the factor by which `deff` inflates the mean's variance, and
#' the two errors compound rather than cancel. The standard error is therefore
#' `s * sqrt(deff / n) * sqrt((n - 1) / (n - deff))`, and `NA` where
#' `deff >= n` (the cell then holds one observation's worth of information
#' and `s` carries none about `sigma`). This is the package's own derivation,
#' not taken from a reference. Measured 95% interval coverage at `n = 30`
#' over 20,000 replicates: 0.921, 0.844 and 0.628 at `rho` = 0.2, 0.5 and 0.8
#' with `s * sqrt(deff / n)` alone, against 0.948, 0.950 and 0.949 with the
#' rescaling.
#'
#' You may also pass a fixed numeric design effect (e.g. `deff = 2`) to
#' uniformly inflate standard errors: an externally supplied constant is
#' applied as `sd * sqrt(deff / n)`, exactly `sqrt(deff)` times the naive SE in
#' every cell. (The `E[s^2]` correction that the estimated design effects also
#' apply is derived from within-cell correlation and would not be justified for
#' a number the caller chose.)
#'
#' Even with the Kish correction, the adjusted SE is an approximation.
#' For rigorous inference under spatial dependence, consider fitting an
#' explicit spatial covariance model (e.g. via \code{\link{fit_bayesian_spatial_model}}).
#'
#' @section What the standard error estimates:
#' The `..se_*` columns are the standard error of the cell mean **as an
#' estimate of the population (grand) mean**: the unconditional quantity, in
#' which the cell's own realised deviation is part of the error. That is the
#' right quantity when cells are treated as samples from a common population,
#' and the design-effect correction is calibrated for it: measured 95% interval
#' coverage of the grand mean is 0.95 with `deff = "kish"` (and 0.29 with the
#' naive SE) on exchangeable within-cell correlation, and 0.93 with
#' `deff = "variogram"` on a simulated Gaussian field.
#'
#' It is **not** the standard error of the cell's own mean (the block average
#' over that cell), which is what a cell-level map or a regression on cell
#' values usually wants. For that quantity the naive `sd / sqrt(n)` is the
#' better of the two on offer: measured coverage 0.95 under exchangeable
#' within-cell correlation, against very nearly 1.00 for the
#' design-effect-corrected SE, which is too wide by the factor
#' `sqrt(deff / (1 - rho))` (4.6 at 20 points a cell and `rho = 0.5`). That 0.95
#' is exact under the exchangeable model and holds under a spatial covariance
#' model only when the cell's points are spread through the cell; with
#' *clustered* sampling inside a cell it is anticonservative for the block
#' average too (measured 0.58), and the right quantity there is a
#' block-kriging variance, which this function does not compute. Use `deff`
#' when the cell means feed a population-level inference; leave it at 1 when
#' they are measurements of the cells themselves and the sampling within cells
#' is reasonably uniform.
#'
#' @section Design effects and variable types:
#' `deff = "kish"` estimates a separate ICC for each summarised column, the
#' response and every predictor, and applies each to that column alone:
#' listing or removing a predictor leaves every other column's standard error,
#' effective n and interval as they were. (One ICC used to be pooled over all
#' the predictors, which was right for a column only when they all shared it:
#' beside a predictor with an ICC of 0.87, one with an ICC of 0 had its
#' standard error multiplied by 3.1.) `deff = "variogram"` fits or
#' accepts **one** correlation function and applies it to every numeric column,
#' because a variogram is a property of the field being modelled rather than of
#' a variable type; a predictor whose spatial structure differs markedly from
#' the response's will have its SE corrected by the response's correlation.
#' The internally estimated variogram is fitted to the **response itself**,
#' never to OLS residuals, whatever `predictor_vars` holds: the `..se_resp_*`
#' columns estimate the grand mean of the response, so the correlation to
#' correct for is the response's own. (A residual variogram, whose correlation
#' is that of the part the predictors do not explain, is weaker; using it here
#' dropped grand-mean coverage from 0.93 to 0.51 the moment a predictor was
#' listed.) Pass `sac` explicitly when you want a different variogram, such as
#' a residual one from `estimate_sac_range(..., predictor_vars = )`. A `sac`
#' whose `attr(sac, "detrended")` is `TRUE` is used as given, but when it
#' corrects response columns a warning says that their standard errors are
#' understated, as [kriging_adequacy()] warns about the same mismatch.
#'
#' @section Confidence intervals:
#' With `conf_level` set, every numeric response and predictor column gains
#' four more columns: `..neff_*`, that column's effective sample size in the
#' cell (its non-missing count over its design effect, the per-column version
#' of `cell_weight`); `..df_*`, the degrees of freedom the interval uses; and
#' `..ci_lo_*` / `..ci_hi_*`, a t interval for the **cell mean as an estimate
#' of the grand mean**. The estimand is the same as `..se_*`'s, so everything
#' in "What the standard error estimates" applies to it, including that it is
#' not an interval for the cell's own block average. The interval is
#' `mean +/- qt((1 + conf_level) / 2, df) * se`, centred on the plain mean of
#' the column's non-missing values whatever `agg_funs` computes, and is `NA`
#' wherever the standard error is (a single observation; complete redundancy
#' under `deff`). `..neff_*` and `..df_*` are `NA` where the column has one
#' non-missing value or none, like the standard error; `cell_weight` still
#' counts such a cell (1, or `1 / deff` for a numeric `deff`), so the two
#' differ there.
#'
#' The degrees of freedom are **not** `neff - 1`. The interval's spread comes
#' from the within-cell sample variance, and under exchangeable correlation
#' (`deff = "kish"`) that variance keeps its `n - 1` degrees of freedom
#' whatever the design effect. The design effect inflates the mean's
#' variance and biases `s^2`, both of which the standard error already
#' corrects, and the resulting pivot is exactly t on `n - 1` df. Measured 95%
#' coverage of the grand mean on the Kish path, 20 cells of 20 at an ICC of
#' 0.2 / 0.6 / 0.9: 0.954 / 0.953 / 0.952 with `n - 1`, against
#' 0.992 / 1.000 / 1.000 with `neff - 1`, which is not an interval so much as
#' a refusal to say anything. The effective-sample-size degrees of freedom of
#' Faes et al. (2009) belong to a mean estimated across many correlated units
#' whose variance is estimated from the spread between them; they do not
#' transfer to a single cell's mean with a variance estimated from inside it.
#' `..df_*` is therefore `n - 1` at `deff = 1`, for a numeric `deff` and for
#' `"kish"`. Under `deff = "variogram"` correlation decays with distance and
#' the within-cell variance loses degrees of freedom to it; `..df_*` is then
#' the Satterthwaite (1946) moment-matched df of `s^2` under the fitted
#' correlation matrix, a fraction of `n - 1` that shrinks as the range grows.
#' Measured on a Gaussian field with an exponential range of 150 and 400 on a
#' 1000-unit domain, 16 cells of about 25 points: 0.960 and 0.958 with that
#' df, against 0.931 and 0.918 with `n - 1`, and 0.34 and 0.19 for the naive
#' `deff = 1` interval. The interval is only as good as the design effect
#' under it: a mis-specified variogram, or an ICC estimated from too few
#' cells, moves the coverage with it.
#'
#' @param assigned_points_sf An sf object with a cell identifier column.
#' @param response_var Optional response column name for per-cell
#'   aggregation: a single character string. A name that is not a column, or
#'   a column that is not numeric, is skipped with a warning.
#' @param predictor_vars Optional predictor column names for per-cell
#'   aggregation. Names that are not columns, and columns that are not
#'   numeric, are skipped with a warning.
#' @param id_col Preferred name of the polygon/cell ID column.
#' @param agg_funs Named list of aggregation functions. Default
#'   \code{list(mean = \(x) mean(x, na.rm = TRUE))}. Additional common options:
#'   \code{median}, \code{sum}, \code{sd}. A single function
#'   (\code{agg_funs = median}) or a character vector of function names
#'   (\code{c("median", "sum")}) is also accepted. A single function is named
#'   after the expression passed: a name gives that name (\code{median}, or
#'   \code{f} for a variable \code{f} holding a function),
#'   \code{stats::median} gives \code{median}, and any other expression gives
#'   \code{agg1}; pass a named list to choose the name. Anything else falls
#'   back to the default mean with a warning. The functions receive each
#'   column with its infinite values already set to \code{NA} (see
#'   "Details"), so they need \code{na.rm = TRUE} as they do for missing
#'   values.
#' @param cells_sf Optional polygon sf layer to join cell geometries onto
#'   the output. When supplied, the return value is an sf object with
#'   the polygon geometry from cells_sf, with one row per cell in `cells_sf`.
#'   Cells that no feature fell in are kept, with `NA` summaries. Duplicate
#'   ID values in `cells_sf` would multiply those rows, so they are reported
#'   with a warning. Its ID column is the first of `id_col`, `"poly_id"`,
#'   `"polygon_id"`, `"id"`, `"cell_id"` and `"grid_id"` it carries, the
#'   list and order [assign_features_to_polygons()] reads the polygons' IDs
#'   from. A summarised ID that matches no cell is reported with a warning,
#'   since the join drops it with its points; a `cells_sf` with none of
#'   those columns, or one that is not an sf object, gives a warning and a
#'   plain data frame (an error with `area = TRUE`). When NULL (default), a
#'   plain data.frame/tibble is returned (previous behaviour).
#' @param deff Design-effect adjustment for standard errors. One of:
#'   \describe{
#'     \item{`1` (default)}{No adjustment; the classic IID standard error
#'       `sd / sqrt(n)`, which assumes the observations within a cell are
#'       independent. No `"deff_applied"` attribute is attached, so a result
#'       carrying none was computed this way.}
#'     \item{`"variogram"`}{Compute a per-cell design effect from a fitted
#'       variogram: for `n` points in a cell with correlation matrix `R`, the
#'       effective sample size of the mean is `n^2 / sum(R)`, so
#'       `deff = sum(R) / n`. This generalises Kish (substituting a constant
#'       off-diagonal correlation recovers `1 + (n - 1) * rho` exactly) but
#'       lets correlation decay with distance, which matters increasingly as
#'       cells get larger and Kish's single-`rho` assumption degrades. Supply
#'       the fit via `sac`, or it is estimated when `response_var` is given and
#'       'gstat' is available. Exponential, spherical and Gaussian models are
#'       supported, with a nugget and with several structured components
#'       (each weighted by its partial sill) and with gstat's 2-D geometric
#'       anisotropy (`vgm(..., anis = c(angle, ratio))`), applied as gstat
#'       applies it. A request that cannot be met falls back to `"kish"`,
#'       with a warning of class `"spatialkit_deff_fallback"` that gives the
#'       reason and the ICCs: a model of any other family (named in the
#'       warning), or no usable model (none supplied and none could be
#'       estimated, or a rejected fit that could not be replaced by an
#'       estimate). The standard errors, intervals and `cell_weight` are then
#'       Kish's, and the uncorrected ones only where Kish's own ICC is 0; the
#'       `deff_method` column says which on every row (see "Value"). The
#'       variogram route needs an identified range, which a layer does not
#'       always have, so read `deff_method` for the fallback. A model with
#'       no structured component (a pure
#'       nugget) implies that distinct observations are uncorrelated, so it
#'       is applied as a design effect of 1 in every cell, not treated as a
#'       fallback. That includes the model [estimate_sac_range()] returns
#'       when it finds no spatial structure (an `NA` range whose
#'       `rejected_reason` says so): the answer is then a design effect of
#'       1, `deff_method` reads `"variogram"`, and no fallback warning is
#'       raised. Points with empty
#'       geometry count towards their cells' values but not towards the
#'       correlation, with a warning.}
#'     \item{`"kish"`}{Estimate an intra-class correlation (ICC) for each
#'       summarised column from the grouped data using a one-way
#'       random-effects ANOVA decomposition (one ICC for the response and one
#'       for each predictor), then apply Kish's formula per
#'       cell: `deff_i = 1 + (n_i - 1) * rho`. The ICC is the ANOVA
#'       (method-of-moments) estimator with Donner's `n0` for unequal cell
#'       sizes, not the REML estimate a mixed model returns: on a single
#'       unbalanced draw the two can differ by 0.1--0.2 (one check with cell
#'       sizes 3 to 77 and a true ICC of 0.5 gave 0.33 against REML's 0.51),
#'       while balanced designs agree to about 0.01. Neither is wrong, so do
#'       not read the difference as a defect. Each column's ICC is used
#'       for that column's SEs and for no other: the response's ICC is never
#'       applied to a predictor, nor one predictor's to another. It is the
#'       response's ICC that sets `cell_weight` whenever a response was
#'       given, and the first predictor's own ICC otherwise.
#'       Requires at least 2 cells with 2+ observations and at least 2 residual
#'       degrees of freedom (`N - k >= 2`); the ICC is taken as 0 (no
#'       correction for that column) otherwise, and likewise when the
#'       estimate itself comes out at or below 0. The `"deff_applied"`
#'       attribute is attached when any ICC is positive, so when all are
#'       0 there is none.}
#'     \item{A positive number}{Applied as a uniform design effect to every
#'       cell, as `sd * sqrt(deff / n)`, exactly `sqrt(deff)` times the
#'       naive SE. Use when you have an external estimate of the design
#'       effect. Anything that is not a single finite number `>= 1`
#'       (including a value below 1, which would *shrink* the standard
#'       errors, and `NA` or `Inf`) is refused with a warning and replaced
#'       by 1.}
#'   }
#' @param sac Optional `sac_range` object from [estimate_sac_range()], used
#'   when `deff = "variogram"`. Supplying one avoids re-fitting the variogram
#'   and lets you inspect the fit the design effect is based on. A `sac_range`
#'   whose fit was *rejected* (its value is `NA` and it carries a
#'   `rejected_reason` attribute) carries no usable correlation function, so
#'   it is set aside rather than correcting by a shape that was not trusted
#'   enough to report a range: the variogram is then estimated as if no `sac`
#'   had been given, with a plain warning saying so, or, where that is not
#'   possible, `deff` falls back to `"kish"` with the fallback warning, which
#'   names the rejection. The exception is a `sac` refused because no spatial
#'   structure was found, whose model is a pure nugget: it is used, and gives
#'   a design effect of 1. A `sac` with no `variogram_model` attribute -- a plain
#'   number or a `units` object, say -- is a range without a correlation
#'   function, and is set aside the same way. A `sac` fitted to residuals
#'   (`attr(sac, "detrended")` `TRUE`) is used as given, with a warning when
#'   it corrects response columns (see "Design effects and variable types").
#' @param deff_max_n Cells with more than this many points are subsampled
#'   before forming the `n x n` correlation matrix used by
#'   `deff = "variogram"`. Default 500. It must be a single number of at
#'   least 2 when `deff = "variogram"`; anything else is an error.
#' @param quiet Logical; suppress this function's progress \code{message()}s.
#'   It does not silence R warnings, nor the package's console log echo
#'   (see \code{\link{spatialkit_quiet}} for that). Default \code{TRUE},
#'   unlike the tessellation functions, whose default is \code{FALSE}.
#' @param area Logical, default `FALSE`.  With `TRUE`, and `cells_sf`
#'   supplied, the result gains `cell_area` (each cell's planar area in the
#'   squared units of `cells_sf`'s CRS) and `n_per_area` (the count of rows
#'   in the cell over that area, a point density; a rate of anything else
#'   is that thing's `agg_funs` sum over `cell_area`).  The request is
#'   **refused** with an error, not answered with a number, when the cells'
#'   CRS distorts areas across them by more than 1 percent, measured as the
#'   spread of planar-to-geodesic area ratios over the cells: a density is a
#'   comparison between cells and is meaningless where the map scale differs
#'   from one cell to the next.  Inside a UTM zone the spread is under 0.3
#'   percent and the request goes through; the conterminous United States
#'   forced into one zone (14 percent), or a few degrees of latitude in Web
#'   Mercator (4 percent at 48N), does not.  [ensure_projected()] with
#'   `purpose = "area"` chooses an equal-area CRS for lon/lat input; build the
#'   cells in it.  Lon/lat cells are measured geodesically instead: their
#'   `cell_area` is `sf::st_area()`'s area in square metres (on the sphere
#'   with s2, sf's default), not a planar area in squared degrees, and the
#'   distortion check passes by construction; with s2 switched off, sf needs
#'   the lwgeom package for that area and the request is refused without
#'   it.  The measured spread is attached as `attr(, "area_error")`.
#' @param conf_level Optional confidence level in (0, 1), such as `0.95`.
#'   When given, every numeric response and predictor column also gets
#'   `..neff_*`, `..df_*`, `..ci_lo_*` and `..ci_hi_*` (see "Confidence
#'   intervals"). Default `NULL`: no interval columns, and the frame is
#'   exactly what it was before this argument existed.
#' @param estimate Which value of a cell the frame reports: `"mean"`
#'   (default) or `"pooled"`. With `"mean"` the frame is the one described
#'   above. With `"pooled"` every numeric response and predictor column also
#'   gets a pooled value, `w * m + (1 - w) * t`: the cell's own mean `m` and a
#'   target `t`, the mean of the rows of the 8 cells whose observations'
#'   centroids are nearest (all the other cells when there are fewer than 10;
#'   the polygons are not read). The own weight is
#'   `w = tau2 / (tau2 + sigma2 / n)`, where `n` is the cell's count,
#'   `sigma2` the spread of the observations about their cell's mean and
#'   `tau2` the mean squared distance between a cell's level and its target.
#'   Each column is pooled on its own rows: those with a value, a cell and a
#'   location. `agg_funs`, `deff`, `sac` and `conf_level` do not enter the
#'   pooled columns, and every other column is computed as with `"mean"`.
#'   Lon/lat coordinates, and coordinates without a CRS that look like
#'   degrees, are projected for the distances between centroids, as
#'   [resolution_profile()] projects them. The projection is chosen from
#'   every row of the layer, those without a value or a cell included, so on
#'   such input the rows a column does not use can still decide the CRS its
#'   distances are measured in; a layer projected beforehand is used in its
#'   own CRS. See "Value" for the columns, the attributes and the warnings.
#' @param obs_level Optional level in (0, 1), such as `0.90`, of an interval
#'   for one further observation in a cell, located like the data. It needs
#'   `estimate = "pooled"`: the interval is the pooled value plus the
#'   `(1 - obs_level) / 2` and `(1 + obs_level) / 2` quantiles of the
#'   leave-one-out residuals of the pooled values, so its two offsets are the
#'   same in every cell. It is a statement about observations on average,
#'   not about each cell. It is a separate argument from `conf_level`, whose
#'   interval is for a mean and is not a range for an observation. Default
#'   `NULL`: no such interval.
#' @return A tibble/data.frame (or sf if cells_sf given) with per-cell
#'   summaries: the ID column, `n` (rows in the cell), one column per
#'   `agg_funs` entry per variable, `..sd_*` / `..se_*` for every numeric
#'   response and predictor, `..neff_*` / `..df_*` / `..ci_lo_*` /
#'   `..ci_hi_*` for the same columns when `conf_level` is given,
#'   `cell_weight`, `deff_applied` and `deff_method` when a design effect
#'   was requested (any `deff` other than 1), and `cell_area` / `n_per_area`
#'   when `area = TRUE`. An input column also called `n` is not allowed to
#'   shadow the count.
#'
#'   `deff_applied` is `TRUE` on every row when a design-effect correction
#'   was applied (exactly when the `"deff_applied"` attribute below is
#'   attached; under `"kish"`, when any standard-error column was corrected)
#'   and `FALSE` when the standard errors are the uncorrected ones: a
#'   refused `deff`, or `"kish"` ICCs of 0 for every column (and `NA` on a
#'   `cells_sf` row no point fell in). `deff_method` names the correction
#'   the standard errors carry, since their column names do not change with
#'   it: `"variogram"`, `"kish"` or `"fixed"` when the one requested was
#'   applied, `"none"` when none was, and `"kish (variogram fallback)"` or
#'   `"none (variogram fallback)"` when a `"variogram"` request could not be
#'   met and Kish stood in for it (the second when Kish's ICCs were 0 too).
#'   So `deff_applied` is `TRUE` on a Kish fallback row: test
#'   `deff_method == "variogram"` to keep only rows the variogram corrected.
#'   Unlike the attribute, both columns survive `rbind()` and
#'   `dplyr::bind_rows()` of many results. A fallback is also signalled by a
#'   warning of class `"spatialkit_deff_fallback"` (a Kish ICC of 0 under a
#'   `"kish"` request is reported on `attr(, "icc")` instead), which
#'   `tryCatch(spatialkit_deff_fallback = )` catches without matching the
#'   message. It is raised when the correction requested was not the one
#'   applied: a refused number (the standard errors are then the uncorrected
#'   ones) or a `"variogram"` request that fell back to Kish; the condition
#'   of the second carries `requested` and `used`.
#'
#'   When a correction was actually applied, an attribute `"deff_applied"` is
#'   attached recording it: `method` plus `icc_resp`/`icc_pred` and `deff`
#'   for `"kish"` (`icc_pred` is a vector named by predictor column, one ICC
#'   each; `deff` is the primary variable's per-cell design effect,
#'   all 1 when only a predictor's ICC was positive; `requested` and
#'   `fallback_reason` are added when a `"variogram"` request fell back to
#'   it),
#'   `deff`/`deff_rows`/`rbar`/`crs`/`max_n` for `"variogram"` (`deff` is
#'   the design effect at the primary variable's non-missing count per cell,
#'   `deff_rows` at the cell's row count, which is the vector the log line
#'   summarises as a median and a max), and `deff` alone for a fixed number.
#'   When `cells_sf` is supplied, *every* per-cell vector in that attribute
#'   (`deff`, `deff_rows` and `rbar` alike) is realigned to the joined row
#'   order, so `deff[i]` and `rbar[i]` still describe row `i`; cells with no
#'   observations carry `NA`. No attribute is attached when no correction was
#'   applied: `deff = 1`, or a `deff = "kish"` request (or Kish fallback)
#'   whose ICCs are all 0. A `deff = "kish"` request, and a `"variogram"`
#'   request that fell back to it, always
#'   records the ICCs it estimated on an attribute `"icc"` (`resp`, and
#'   `pred`, a vector named by predictor column; `NA` for a variable type
#'   with no numeric column), whether or not
#'   they were positive enough to apply, so a result with no `"deff_applied"`
#'   still says what the ICCs came out as.
#'
#'   The ID column keeps its input type when `cells_sf`'s ID column and the
#'   summarised IDs already have the same class. When the classes differ, both
#'   are coerced to character in order to join (logged as a warning), and the
#'   returned ID column is therefore character. Whole numbers are written out
#'   in full for that (`"100000"`, never `"1e+05"`), so an integer and a
#'   double ID of the same cell still match.
#'
#'   With `estimate = "pooled"`, each numeric column `v` with its prefix
#'   (`resp_` or `pred_`) gains `..pooled_<prefix><v>` (the pooled value),
#'   `..pooled_w_<prefix><v>` (the own weight, 0 to 1; it is not a regression
#'   weight) and `..pooled_target_<prefix><v>` (the target), and with
#'   `obs_level` also `..obs_lo_<prefix><v>` and `..obs_hi_<prefix><v>`.
#'   They follow every column listed above and precede `cell_area` and
#'   `n_per_area`, the response's first, then each predictor's. A cell with
#'   no row that has a value, a cell and a location gets `NA`, never a value
#'   made from its neighbours. Two attributes are attached.
#'   `attr(, "pooled")`, of class `"spatialkit_pooled"`, records what was
#'   fitted: `toward`, `neighbours`, `crs` (the CRS the centroid distances
#'   were measured in), `obs_level`, `conf_level`, `intervals` and `columns`,
#'   a list with one entry per pooled column holding `column`, `role`
#'   (`"response"` or `"predictor"`), `n` (rows used), `cells`, `df` (values
#'   beyond one per cell, which `sigma2` rests on), `neighbours` (the cells
#'   in each target), `sigma2`, `tau2`, `own_weight` (the median weight),
#'   `share_alone` (the share of rows that are the only one in their cell),
#'   `loo` and `loo_se` (the leave-one-out mean squared error of the pooled
#'   value, and its standard error over rows), `loo_one_cell` (the same error
#'   with one mean for all rows), `structured` and `structure_p` (whether the
#'   cells predict a left-out observation better than one mean, and the
#'   p-value of the F test of the cell means against one mean, which has to
#'   be below 0.05 as well), `offsets` (the two offsets of the interval),
#'   `repeated`, `skipped`, `constant` and `status`. `attr(, "intervals")` is
#'   the same data frame as `attr(, "pooled")$intervals`: one row per family
#'   of interval columns in the frame, saying what it covers, what it is
#'   centred on, its level and its method. Away from the observations no
#'   interval here is calibrated.
#'
#'   Where the pooled columns of a column cannot be computed they are `NA`,
#'   with a warning of class `"spatialkit_pool_unavailable"` and not an
#'   error: a layer without geometry, no row with a value, a cell and a
#'   location, no cell with two values, or more cells than the neighbour
#'   search may compare without the 'FNN' package
#'   (`options(spatialkit.exact_pairs_max = )`, 2e9 pairs by default). The
#'   other warnings of the pooled columns are of class
#'   `"spatialkit_pool_degenerate"` (one cell, no spread inside any cell, or
#'   every own weight 0), `"spatialkit_pool_few_cells"` (2 to 9 cells, where
#'   each target is the mean of all the other cells),
#'   `"spatialkit_pool_few_df"` (`sigma2` rests on fewer than 30 values
#'   beyond one per cell), `"spatialkit_no_structure"` (`structured` is
#'   `FALSE`: the pooled values predict a left-out observation no better
#'   than one mean for all rows, or the F test does not tell the cell means
#'   from one mean; the text says which, and `beats` in the condition is
#'   `TRUE` for the second), `"spatialkit_pool_rows_skipped"` (rows with a
#'   value and a cell and no location) and, with `obs_level`,
#'   `"spatialkit_repeated_locations"` and `"spatialkit_few_residuals"`.
#'   Each is raised once per call and reason, for all the columns it applies
#'   to, and carries them in `columns`; the numbers of each column are in
#'   `attr(, "pooled")`. The two about few cells and few values are raised
#'   for columns whose weights were estimated, and the two of `obs_level` for
#'   columns whose interval is read from residuals, so a constant column
#'   raises `"spatialkit_pool_degenerate"` and no other.
#' @examples
#' library(sf)
#' set.seed(1)
#' n <- 200
#' east  <- runif(n, 0, 100)
#' north <- runif(n, 0, 100)
#' # A response with spatial structure, so the within-cell ICC is not zero.
#' pts <- st_as_sf(
#'   data.frame(x = 5e5 + east, y = 5e6 + north,
#'              val = 5 + 0.02 * east + 0.02 * north + rnorm(n, sd = 0.5)),
#'   coords = c("x", "y"), crs = 32632
#' )
#' bnd <- st_sf(geometry = st_sfc(st_polygon(list(rbind(
#'   c(5e5, 5e6), c(5e5 + 100, 5e6), c(5e5 + 100, 5e6 + 100),
#'   c(5e5, 5e6 + 100), c(5e5, 5e6)
#' ))), crs = 32632))
#' grid <- create_grid_polygons(bnd, target_cells = 9, type = "square")
#' assigned <- assign_features_to_polygons(pts, grid)
#'
#' # IID standard errors (default) vs Kish design-effect adjustment.  The
#' # correction is large here, and meant to be: about 22 points per cell that
#' # all share the cell's part of the trend carry far fewer than 22
#' # independent pieces of information about the cell mean.
#' naive <- summarize_by_cell(assigned, response_var = "val")
#' kish  <- summarize_by_cell(assigned, response_var = "val", deff = "kish")
#' data.frame(n = naive$n,
#'            se_naive = naive$..se_resp_val,
#'            se_kish  = kish$..se_resp_val)
#' attr(kish, "deff_applied")   # method, icc_resp, icc_pred, per-cell deff
#'
#' # A 95% interval for each cell mean as an estimate of the grand mean, on
#' # the Kish-corrected standard error and n - 1 degrees of freedom.
#' ci <- summarize_by_cell(assigned, response_var = "val", deff = "kish",
#'                         conf_level = 0.95)
#' ci[, c("poly_id", "n", "resp_mean_val", "..neff_resp_val", "..df_resp_val",
#'        "..ci_lo_resp_val", "..ci_hi_resp_val")]
#' @references
#' Faes, C., Molenberghs, G., Aerts, M., Verbeke, G. and Kenward, M. G. (2009).
#' The effective sample size and an alternative small-sample
#' degrees-of-freedom method. \emph{The American Statistician}, 63(4),
#' 389--399. \doi{10.1198/tast.2009.08196}
#'
#' Satterthwaite, F. E. (1946). An approximate distribution of estimates of
#' variance components. \emph{Biometrics Bulletin}, 2(6), 110--114.
#' \doi{10.2307/3002019}
#' @family aggregation
#' @seealso [assign_features_to_polygons()], which produces the input layer;
#'   [build_tessellation()] for the cells themselves.
#' @export
summarize_by_cell <- function(assigned_points_sf,
                              response_var   = NULL,
                              predictor_vars = NULL,
                              id_col         = "poly_id",
                              agg_funs       = list(mean = function(x) mean(x, na.rm = TRUE)),
                              cells_sf       = NULL,
                              deff           = 1,
                              sac            = NULL,
                              deff_max_n     = 500L,
                              quiet          = TRUE,
                              conf_level     = NULL,
                              area           = FALSE,
                              estimate       = c("mean", "pooled"),
                              obs_level      = NULL) {
  .msg <- function(...) if (!quiet) message(...)
  estimate <- match.arg(estimate)
  if (!is.null(conf_level) &&
      (!is.numeric(conf_level) || length(conf_level) != 1L ||
       !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1))
    stop("summarize_by_cell(): `conf_level` must be a single number strictly ",
         "between 0 and 1 (0.95 for a 95% interval), or NULL for no intervals.",
         call. = FALSE)
  # Strictly between 0 and 1 as the interval reads the level: through its tail
  # probability rounded to 12 decimals (.obs_offsets()).  That is 0.5 for a
  # level of 0 and anything below, 0 for a level of 1 and anything above, and
  # also for a level within 1e-12 of either, which would give the median
  # twice or the range of the residuals.
  if (!is.null(obs_level) &&
      (!is.numeric(obs_level) || length(obs_level) != 1L || !is.finite(obs_level) ||
       round((1 - obs_level) / 2, 12) <= 0 || round((1 - obs_level) / 2, 12) >= 0.5))
    stop("summarize_by_cell(): `obs_level` must be a single number strictly ",
         "between 0 and 1 (0.90 for a 90% interval), or NULL for no interval.",
         call. = FALSE)
  # The interval for one observation has one centre, the pooled value.  With
  # plain means it is refused, not centred on the plain mean without a word.
  if (!is.null(obs_level) && identical(estimate, "mean"))
    stop("summarize_by_cell(): `obs_level` needs estimate = \"pooled\": the interval ",
         "for an observation is built from the leave-one-out residuals of the pooled ",
         "value and is centred on it. `conf_level` gives the plain mean's t interval, ",
         "which is not a range for an observation.", call. = FALSE)
  if (!is.logical(area) || length(area) != 1L || is.na(area))
    stop("summarize_by_cell(): `area` must be TRUE or FALSE.", call. = FALSE)
  # c("val", "flag") used to stop with "the condition has length > 1", and a
  # non-character name was looked up and skipped without a word.
  if (!is.null(response_var) &&
      (!is.character(response_var) || length(response_var) != 1L ||
       is.na(response_var) || !nzchar(response_var)))
    stop("summarize_by_cell(): `response_var` must be a single column name ",
         "(a character string), or NULL; list further numeric columns in ",
         "`predictor_vars`.", call. = FALSE)
  if (!is.null(predictor_vars) && !is.character(predictor_vars))
    stop("summarize_by_cell(): `predictor_vars` must be a character vector of ",
         "column names, or NULL.", call. = FALSE)
  # Whether the caller asked for a design effect at all.  Only then does the
  # result carry the per-row `deff_applied` column, so a default call returns
  # exactly the frame it always has.
  deff_requested <- !(is.numeric(deff) && length(deff) == 1L &&
                        isTRUE(deff == 1))
  # A density needs the cells, and it needs them in a CRS whose areas are
  # comparable.  Both are checked before any summary is computed, so a
  # request that cannot be honoured fails at once rather than after the
  # aggregation.
  if (isTRUE(area)) {
    if (!inherits(cells_sf, "sf"))
      stop("summarize_by_cell(): `area = TRUE` needs `cells_sf`, the polygon layer ",
           "the cells' areas are read from.", call. = FALSE)
    if (is.na(sf::st_crs(cells_sf)))
      stop("summarize_by_cell(): `area = TRUE` needs `cells_sf` to carry a CRS: ",
           "without one there is no way to tell whether its areas are ",
           "comparable between cells. Set it with sf::st_crs().", call. = FALSE)
    area_err <- .crs_area_error(cells_sf)
    if (!is.finite(area_err))
      stop("summarize_by_cell(): `area = TRUE`: the area distortion of `cells_sf`'s ",
           "CRS could not be measured (no polygon areas, or geodesic areas ",
           "unavailable), so a density cannot be vouched for.", call. = FALSE)
    if (area_err > .area_error_tol)
      stop(sprintf(paste0("summarize_by_cell(): `area = TRUE` refused: the CRS of ",
                          "`cells_sf` (%s) distorts areas across the cells by up ",
                          "to %.1f%% (planar against geodesic), so a density or ",
                          "rate computed in it would not be comparable between ",
                          "cells. Build the cells in an equal-area CRS -- ",
                          "ensure_projected(purpose = \"area\") picks one for ",
                          "lon/lat input -- and pass those."),
                   .fold_crs_label(cells_sf), 100 * area_err), call. = FALSE)
  }
  df <- sf::st_drop_geometry(assigned_points_sf)

  # --- locate ID column ---
  id_candidates <- unique(c(id_col, "poly_id", "polygon_id", "cell_id"))
  id_found <- id_candidates[id_candidates %in% names(df)]
  if (length(id_found) == 0L) {
    stop("summarize_by_cell(): could not find an ID column. Looked for: ",
         paste(id_candidates, collapse = ", "), ".")
  }
  id_col <- id_found[[1]]
  .msg(sprintf("summarize_by_cell(): using id_col = '%s'", id_col))

  # --- validate agg_funs ---
  # A bare function (`agg_funs = median`) or a vector of function names
  # (`"median"`) is what the documented options invite, and both used to be
  # replaced by the mean with only a log line: resp_mean_v = 1.228 came back
  # for a median of 0.927.
  if (is.function(agg_funs)) {
    # Named after the expression passed; `stats::median` used to give
    # resp_agg1_v where `median` gave resp_median_v.
    fn_expr  <- substitute(agg_funs)
    fn_name  <- if (is.name(fn_expr)) as.character(fn_expr)
                else if (is.call(fn_expr) && length(fn_expr) == 3L &&
                         (identical(fn_expr[[1L]], as.name("::")) ||
                          identical(fn_expr[[1L]], as.name(":::"))))
                  as.character(fn_expr[[3L]])
                else "agg1"
    agg_funs <- stats::setNames(list(agg_funs), fn_name)
  } else if (is.character(agg_funs) && length(agg_funs) > 0L && !anyNA(agg_funs)) {
    fns <- lapply(agg_funs, function(nm)
      tryCatch(match.fun(nm), error = function(e) NULL))
    if (all(vapply(fns, is.function, logical(1))))
      agg_funs <- stats::setNames(fns, agg_funs)
  }
  if (!is.list(agg_funs) || length(agg_funs) == 0L) {
    .warn_and_log(paste0("summarize_by_cell(): `agg_funs` must be a named list ",
                         "of functions, a function, or the names of functions; ",
                         "got %s. Falling back to the mean."),
                  if (is.character(agg_funs))
                    paste(sprintf("'%s'", agg_funs), collapse = ", ")
                  else paste0("an object of class ", class(agg_funs)[1L]))
    agg_funs <- list(mean = function(x) mean(x, na.rm = TRUE))
  }
  if (is.null(names(agg_funs)) || any(!nzchar(names(agg_funs)))) {
    names(agg_funs) <- paste0("agg", seq_along(agg_funs))
  }

  # --- validate / resolve design effect ---
  use_kish <- identical(deff, "kish")
  use_vgm  <- identical(deff, "variogram")
  if (use_vgm) {
    # A subsample of 0 or 1 points has no pairs: the mean correlation came out
    # NaN, the design effect 1 and the SEs uncorrected, on rows still marked
    # deff_applied = TRUE; NA stopped with "missing value where TRUE/FALSE
    # needed".  Checked only here, where the argument is used.
    .check_scalar(deff_max_n, "deff_max_n", "summarize_by_cell", min = 2,
                  max = .Machine$integer.max)
    # The SE closures below do arithmetic on `deff`; the per-cell variogram
    # values are applied after summarising, so neutralise it here.
    deff <- 1
  }
  if (!use_kish && !use_vgm) {
    # is.finite(): NA and NaN made the test itself NA ("missing value where
    # TRUE/FALSE needed"), and Inf passed it, gave uncorrected SEs (see
    # .se_with_deff()) beside a cell_weight of 0 and a record of deff = Inf.
    if (!is.numeric(deff) || length(deff) != 1L || !is.finite(deff) || deff < 1) {
      # A real warning, not a log line: silently substituting 1 for a value the
      # caller chose means the standard errors are not the ones they asked for,
      # and nothing in the returned object records that (no `deff_applied` is
      # attached when deff == 1).
      .warn_deff_fallback("summarize_by_cell(): `deff` must be a single number >= 1 (and finite), \"kish\" or \"variogram\"; got %s. Falling back to deff = 1, so the standard errors are the uncorrected ones.",
                          paste(format(deff), collapse = ", "))
      deff <- 1
    }
  }

  # Helper: estimate intra-class correlation via one-way random-effects ANOVA.
  # Returns a single rho in [0, 1] or 0 if estimation fails.
  #
  # Called with ONE column at a time (see .kish_iccs() below): an ICC is a
  # property of a variable, and one pooled value for several predictors was
  # right for none of them unless they all shared it.  With an ICC-0 and an
  # ICC-0.87 predictor on 60 cells of 12, the pooled 0.42 multiplied both
  # standard errors by 3.13 where their own ICCs give 1.00 and 8.96, so the
  # grand-mean interval covered 1.000 and 0.492 against 0.958 and 0.950.
  # The body still accepts several columns, and what follows describes that.
  .estimate_icc <- function(df, id_col, var_cols) {
    # Each variable is z-scored before pooling so that variables with larger
    # scales do not dominate the between/within variance decomposition.
    #
    # The pooled group is (VARIABLE, cell), not cell alone.  Pooling m columns
    # under the same cell label treats a cell's m z-scores -- one per variable
    # -- as m draws from one distribution, so independent per-variable cell
    # effects average away in the cell mean and the between-cell sum of
    # squares shrinks by about 1/m.  Measured (200 replicates, true rho = 0.5,
    # m = 4): pooled-by-cell ICC 0.121 against per-column ANOVA 0.495; the
    # predictor deff at n = 20 came out 3.3 where 10.5 was correct, and every
    # predictor SE was 44% too small.  Grouping by (variable, cell) is the
    # ordinary pooled one-way ANOVA across the m*k variable-cell groups, and
    # recovers 0.49 on the same data.
    vals <- numeric(0)
    grps <- character(0)
    for (v in var_cols) {
      # A row with no cell id belongs to no group: keep_unassigned = TRUE
      # output carries NA ids, and an NA level makes group_means[grps] NA
      # downstream ("missing value where TRUE/FALSE needed").
      ok <- !is.na(df[[v]]) & !is.na(df[[id_col]])
      raw <- df[[v]][ok]
      s <- sd(raw)
      if (is.na(s) || s < .Machine$double.eps) {
        # Constant column — skip, contributes no variance information
        next
      }
      vals <- c(vals, (raw - mean(raw)) / s)
      grps <- c(grps, paste(v, as.character(df[[id_col]][ok]), sep = "\r"))
    }
    grps <- factor(grps)
    k <- nlevels(grps)
    N <- length(vals)
    if (k < 2L || N < 4L) return(0)

    ni <- tabulate(grps)
    ni <- ni[ni > 0L]
    # The documentation promises "at least 2 cells with 2+ observations; falls
    # back to deff = 1 otherwise", and the k >= 2 / N >= 4 test above does not
    # deliver it.  All-singleton groups give SSW = 0, MSW = 0 and therefore an
    # ICC of exactly 1 -- a design effect of n from data carrying no
    # within-cell information at all.  A single group with two members gives an
    # ICC from one within-group degree of freedom, which is then applied to
    # every cell.  Both are noise dressed as a measurement.
    if (sum(ni >= 2L) < 2L || (N - k) < 2L) return(0)
    grand_mean <- mean(vals)
    group_means <- tapply(vals, grps, mean)

    # Between-group and within-group sum of squares
    SSB <- sum(ni * (group_means - grand_mean)^2)
    SSW <- sum((vals - group_means[grps])^2)

    MSB <- SSB / (k - 1)
    MSW <- if ((N - k) > 0) SSW / (N - k) else 0

    # Unbalanced n0 (Donner 1986)
    n0 <- (N - sum(ni^2) / N) / (k - 1)
    if (n0 < 1) return(0)

    rho <- (MSB - MSW) / (MSB + (n0 - 1) * MSW)
    # Clamp to [0, 1] — negative ICC implies no positive autocorrelation
    max(0, min(1, rho))
  }

  # --- resolve numeric columns for response and predictors ---
  # A column that was asked for and cannot be summarised is a warning, not a
  # progress message: .msg() prints nothing under the default quiet = TRUE,
  # so a misspelt response_var returned a frame with no resp_* column, a
  # cell_weight equal to n and (under "kish") no design effect, without a
  # word.
  .resolve_numeric <- function(df, cols, label) {
    keep <- cols[cols %in% names(df)]
    if (length(keep) == 0L) return(character(0))
    is_num <- vapply(keep, function(nm) is.numeric(df[[nm]]), logical(1))
    if (!all(is_num)) {
      .warn_and_log("summarize_by_cell(): %s column(s) %s are not numeric and are not summarised.",
                    label, paste(sprintf("'%s'", keep[!is_num]), collapse = ", "))
    }
    keep[is_num]
  }

  resp_num <- character(0)
  pred_num <- character(0)

  if (!is.null(response_var)) {
    if (response_var %in% names(df)) {
      resp_num <- .resolve_numeric(df, response_var, "response")
    } else {
      .warn_and_log("summarize_by_cell(): response_var '%s' is not a column of `assigned_points_sf`; no resp_* columns are returned.",
                    response_var)
    }
  }
  if (!is.null(predictor_vars)) {
    absent <- setdiff(predictor_vars, names(df))
    if (length(absent))
      .warn_and_log("summarize_by_cell(): predictor_vars %s are not columns of `assigned_points_sf` and are not summarised.",
                    paste(sprintf("'%s'", absent), collapse = ", "))
    present <- predictor_vars[predictor_vars %in% names(df)]
    if (length(present)) pred_num <- .resolve_numeric(df, present, "predictor")
  }

  # --- infinite values are missing values ---
  # One rule, applied once, before anything reads the columns: an infinite
  # value in a summarised column is treated as NA is.  Left in, one Inf made
  # its cell's mean Inf and its SE NA while ..neff_ and cell_weight still
  # counted the row, and under deff = "kish" it set the ICC to 0 (sd() of the
  # column is NaN), so EVERY cell's standard error came back uncorrected (0.169
  # against 1.365 on the clean cells) with no condition raised.
  # kriging_adequacy(), resolution_profile(), estimate_sac_range() and
  # prep_model_data() already drop such rows, so one pipeline reported two
  # means for the same cell.  Masking `df` here, and not at each site that
  # counts or averages, is what keeps the mean, sd, se, interval, effective n,
  # cell_weight, the Kish ICC and the variogram's row selection on the same
  # rows.  NaN is missing already (is.na(NaN) is TRUE), so it is not counted.
  inf_cols <- character(0)
  inf_note <- character(0)
  inf_rows <- 0L
  for (v in unique(c(resp_num, pred_num))) {
    is_inf <- is.infinite(df[[v]])
    if (!any(is_inf)) next
    inf_cols <- c(inf_cols, v)
    inf_rows <- inf_rows + sum(is_inf)
    inf_note <- c(inf_note, sprintf(
      "%s '%s' %d row(s) in %d cell(s)",
      if (v %in% resp_num) "response" else "predictor", v, sum(is_inf),
      length(unique(df[[id_col]][is_inf]))))
    df[[v]][is_inf] <- NA
  }
  if (length(inf_cols))
    .warn_classed("spatialkit_nonfinite", paste0(
      "summarize_by_cell(): %d infinite value(s) are treated as missing: %s. ",
      "Every statistic of %s (`agg_funs`, sd, se, the interval, the effective ",
      "n, cell_weight and the design effect) leaves them out, as it leaves out ",
      "NA; `n` still counts the rows. Recode or remove them to silence this."),
      inf_rows, paste(inf_note, collapse = "; "),
      if (length(inf_cols) > 1L) "those columns" else "that column",
      data = list(n = inf_rows, columns = inf_cols))

  # --- estimate ICC for Kish correction if requested ---
  # One ICC per summarised column, the response's and each predictor's own.
  # Called here for deff = "kish", and from the variogram block below when a
  # "variogram" request cannot be met and falls back to Kish.
  resp_rho <- NULL
  pred_rho <- NULL
  .kish_iccs <- function() {
    resp <- NULL
    pred <- NULL
    if (length(resp_num) > 0L) {
      resp <- .estimate_icc(df, id_col, resp_num)
      if (resp > 0) {
        .msg(sprintf(
          "summarize_by_cell(): Kish deff using response ICC = %.4f", resp
        ))
      } else {
        .msg("summarize_by_cell(): response ICC estimate <= 0; Kish deff not applied for response.")
      }
    }
    if (length(pred_num) > 0L) {
      # Named by column; the SE closures look their own column up.
      pred <- vapply(pred_num, function(v) .estimate_icc(df, id_col, v),
                     numeric(1))
      names(pred) <- pred_num
      for (v in pred_num) {
        if (pred[[v]] > 0) {
          .msg(sprintf(
            "summarize_by_cell(): Kish deff using predictor ICC = %.4f for '%s'",
            pred[[v]], v
          ))
        } else {
          .msg(sprintf(paste0("summarize_by_cell(): predictor ICC estimate <= 0 for ",
                              "'%s'; Kish deff not applied for it."), v))
        }
      }
    }
    if (is.null(resp) && is.null(pred))
      .msg("summarize_by_cell(): no numeric columns for ICC estimation; deff = 1.")
    # Zero for the type that has no columns
    list(resp = if (is.null(resp)) 0 else resp,
         pred = if (is.null(pred)) 0 else pred)
  }
  if (use_kish) {
    icc_est  <- .kish_iccs()
    resp_rho <- icc_est$resp
    pred_rho <- icc_est$pred
  }
  # Set when a "variogram" request fell back to Kish: why, for the record.
  vgm_fallback <- NULL

  # --- build combined function list for across() ---
  # Each prefix (resp_, pred_) gets user agg_funs + sd + se, with names
  # constructed so that .names = "{.fn}{.col}" produces the correct columns.
  # e.g. fn name "resp_mean_" + col "temp" -> "resp_mean_temp"
  #      fn name "..sd_resp_" + col "temp" -> "..sd_resp_temp"
  .build_funs <- function(prefix, rho = NULL, rbar = NULL, vgm = NULL) {
    fns <- list()
    for (nm in names(agg_funs)) {
      fn <- agg_funs[[nm]]
      # Force early binding of fn in the closure
      fns[[paste0(prefix, nm, "_")]] <- local({
        f <- fn; function(x) f(x)
      })
    }
    fns[[paste0("..sd_", prefix)]] <- function(x) {
      if (sum(!is.na(x)) > 1L) sd(x, na.rm = TRUE) else NA_real_
    }
    # The per-column, per-cell quantities behind the SE and (when asked for)
    # the interval: the column's own non-missing count, sd, design effect and
    # the degrees of freedom of its variance.  One function so that every
    # column built from these numbers is built from the SAME numbers.
    .col_stats <- local({
      .deff <- deff
      .rho  <- rho
      .rbar <- rbar        # per-cell rbar over ALL rows (fast path, no NAs)
      .dfr  <- vgm$df_ratio  # per-cell df ratio over ALL rows, same path
      .vgm  <- vgm         # list(coords, cor_fn, max_n): per-column exact path
      .kish <- use_kish
      .id   <- id_col
      # The ..se_, ..neff_, ..df_ and ..ci_ closures each call this for the
      # same column in the same cell (summarise() runs each of them over
      # every cell in turn), and on the per-column variogram path every call
      # rebuilt the cell's correlation matrix: five builds per column and
      # cell with conf_level set.  The statistics depend only on which rows
      # enter, so they are kept against those rows.
      .seen <- new.env(parent = emptyenv())
      function(x) {
        ok      <- !is.na(x)
        n_valid <- sum(ok)
        if (n_valid <= 1L)
          return(list(n_valid = n_valid, s = NA_real_, deff = NA_real_,
                      df = NA_real_))
        s <- sd(x, na.rm = TRUE)
        # The df of s^2.  Exactly n_valid - 1 under independence, under a
        # caller-supplied constant deff (a uniform inflation of the SE says
        # nothing about the variance estimate) and under Kish's exchangeable
        # correlation; only correlation that decays with distance changes it
        # (see .cor_stats_from_coords()).
        df_i <- n_valid - 1
        # Every path builds the design effect from n_valid -- the observations
        # THIS column has in THIS cell -- never from the cell's row count.  An
        # NA row contributes nothing to the mean, so it cannot contribute to
        # the mean's clustering either.
        if (.kish) {
          # A named `.rho` holds one ICC per column (the predictors'); each
          # column takes its own.  across() runs one closure list over all
          # the predictor columns, so the column is asked for here.
          rho_i  <- if (!is.null(names(.rho))) .rho[[dplyr::cur_column()]] else .rho
          deff_i <- max(1, 1 + (n_valid - 1) * rho_i)
        } else if (!is.null(.rbar)) {
          # Variogram path.  The mean pairwise correlation is over the
          # locations of the observations that ENTER this column's mean, so
          # when the column has NAs in this cell it is recomputed on those rows
          # alone; the precomputed all-rows value is reused when there are
          # none (the ordinary case, and the same number either way).
          if (all(ok)) {
            g  <- as.character(dplyr::cur_group()[[.id]])
            rb <- if (g %in% names(.rbar)) .rbar[[g]] else NA_real_
            dr <- if (!is.null(.dfr) && g %in% names(.dfr)) .dfr[[g]] else NA_real_
          } else {
            rows <- dplyr::cur_group_rows()[ok]
            key  <- paste(rows, collapse = ",")
            st   <- .seen[[key]]
            if (is.null(st)) {
              st <- .cor_stats_from_coords(.vgm$coords[rows, , drop = FALSE],
                                           .vgm$cor_fn, max_n = .vgm$max_n)
              assign(key, st, envir = .seen)
            }
            rb <- st[["rbar"]]; dr <- st[["df_ratio"]]
          }
          deff_i <- if (is.finite(rb)) min(max(1, 1 + (n_valid - 1) * rb), n_valid)
                    else 1
          if (is.finite(dr)) df_i <- dr * (n_valid - 1)
        } else {
          deff_i <- .deff
        }
        list(n_valid = n_valid, s = s, deff = deff_i, df = df_i)
      }
    })
    .se_of <- local({
      .correct <- use_kish || !is.null(rbar)
      function(st) {
        if (st$n_valid <= 1L) return(NA_real_)
        # NOT s / sqrt(n_valid / deff_i): that corrects the mean's variance for
        # clustering but leaves s^2 biased low by the same clustering.  See
        # .se_with_deff() for the derivation and the measured coverage.  The
        # E[s^2] correction applies only to a design effect this function
        # DERIVED from the data's own clustering; a numeric `deff` supplied by
        # the caller is applied as the uniform inflation it is documented to be.
        .se_with_deff(st$s, st$n_valid, st$deff, correct_s2 = .correct)
      }
    })
    # SE with design-effect adjustment
    fns[[paste0("..se_", prefix)]] <- function(x) .se_of(.col_stats(x))

    if (!is.null(conf_level)) {
      # Effective sample size of THIS column in this cell: its non-missing
      # count over its design effect.  cell_weight is the same number for the
      # primary variable only.
      fns[[paste0("..neff_", prefix)]] <- function(x) {
        st <- .col_stats(x)
        if (st$n_valid <= 1L) return(NA_real_)
        st$n_valid / st$deff
      }
      # The degrees of freedom the interval below uses.  Not neff - 1: the
      # interval's spread is estimated from the within-cell s^2, whose df is
      # n - 1 under exchangeable correlation whatever the design effect (the
      # pivot is exactly t on n - 1 df there), and a Satterthwaite fraction
      # of n - 1 under a distance-decaying correlation.  Measured 95% coverage
      # on the Kish path at rho = 0.2 / 0.6 / 0.9: 0.954 / 0.953 / 0.952 with
      # n - 1, against 0.992 / 1.000 / 1.000 with neff - 1.
      fns[[paste0("..df_", prefix)]] <- function(x) {
        st <- .col_stats(x)
        if (st$n_valid <= 1L) return(NA_real_)
        st$df
      }
      .ci_bound <- local({
        .p <- 1 - (1 - conf_level) / 2
        function(x, side) {
          st <- .col_stats(x)
          se <- .se_of(st)
          if (!is.finite(se) || !is.finite(st$df) || st$df <= 0) return(NA_real_)
          mean(x, na.rm = TRUE) + side * stats::qt(.p, df = st$df) * se
        }
      })
      fns[[paste0("..ci_lo_", prefix)]] <- function(x) .ci_bound(x, -1)
      fns[[paste0("..ci_hi_", prefix)]] <- function(x) .ci_bound(x, +1)
    }
    fns
  }

  # --- per-cell design effect from the variogram ------------------------------
  # Computed separately from the aggregation because the across() closures see
  # only a column of values, never the geometry.  SE scales as sqrt(deff), so
  # the SEs are formed with deff = 1 below and rescaled per cell afterwards.
  vgm_deff <- NULL
  vgm_rbar <- NULL
  vgm_bits <- NULL
  if (use_vgm) {
    # A REJECTED sac is not a fit.  estimate_sac_range() returns NA with a
    # `rejected_reason` when the fitted range runs past the longest lag the
    # variogram covers, or when the optimiser stopped at its iteration limit --
    # and it says in as many words that such a range "is where the optimiser
    # halted rather than a fitted parameter" and must not be used.  It still
    # attaches `variogram_model` so the fit can be INSPECTED, and reading that
    # attribute without checking the value took the refused model as gospel:
    # the correlation function then saturated at every within-cell distance,
    # deff came out equal to n, cell_weight collapsed to 1 and standard errors
    # were inflated 40-50x.  Treat it like no fit at all.
    #
    # What replaces it is said once it is known.  This used to raise the
    # classed fallback warning ("Falling back to deff = 1") right here, and
    # then a variogram estimated from `response_var` was applied after all
    # (deff_applied TRUE, median deff 5.2), so tryCatch() on the class threw
    # a corrected result away; with the estimate rejected too, one fallback
    # gave two R warnings.  Now a replaced `sac` gets a plain warning below,
    # and an unreplaced one is named in the single fallback warning.
    no_fit    <- NULL       # why no model is available, for the warning below
    set_aside <- NULL       # why a supplied `sac` was not used
    set_aside_so <- ""      # and what that means, for the replacement warning
    # One refusal is an answer and not a failed fit.  When
    # estimate_sac_range() finds no spatial structure it returns NA with a
    # pure-nugget model at the variance of the variable (.sac_no_structure(),
    # R/utils.R).  Under that model distinct observations are uncorrelated and
    # the design effect is 1 exactly, which is the pure-nugget route below.
    # Set aside with the other refusals, it raised the fallback warning and
    # handed the request to Kish for a layer the variogram had just answered.
    no_structure <- function(s)
      .sac_no_structure(s) && .vgm_is_pure_nugget(attr(s, "variogram_model"))
    structure_note <- NULL  # the refusal that gave the pure nugget, for the log
    if (no_structure(sac)) {
      structure_note <- as.character(attr(sac, "rejected_reason"))[1L]
    } else if (!is.null(sac) && is.na(suppressWarnings(as.numeric(sac))[1L]) &&
        !is.null(attr(sac, "rejected_reason"))) {
      set_aside <- sprintf("the supplied `sac` reports no usable range (%s)",
                           as.character(attr(sac, "rejected_reason"))[1L])
      set_aside_so <- ", so its fitted variogram model cannot size a design effect"
      sac <- NULL
    } else if (!is.null(sac) && is.null(attr(sac, "variogram_model"))) {
      # A bare range -- a number, a units object -- carries no correlation
      # function, and it was passed over without a word: the design effect
      # came from a variogram estimated here, and nothing said that the value
      # given had not been used.
      set_aside <- sprintf(paste0("the supplied `sac` (%s) carries no fitted ",
                                  "variogram model, and a range alone cannot ",
                                  "size a design effect"),
                           if (is.numeric(sac)) paste(format(sac), collapse = ", ")
                           else sprintf("an object of class %s", class(sac)[1L]))
      sac <- NULL
    }
    # A residual variogram used for the response's own standard errors is
    # used as documented, but said so: kriging_adequacy() warns about the
    # same attribute on the same object.
    detrended_sac <- !is.null(sac) && isTRUE(attr(sac, "detrended"))
    vgm_model <- attr(sac, "variogram_model")
    vgm_crs   <- attr(sac, "crs")
    if (is.null(vgm_model)) {
      # Why there is no model, for the fallback warning below: the only
      # signal used to be a log line, which spatialkit_quiet() hides and no
      # tryCatch() sees, while the standard errors came back uncorrected.
      if (is.null(response_var)) {
        no_fit <- "there is no `response_var` to fit one to"
      } else if (!requireNamespace("gstat", quietly = TRUE)) {
        no_fit <- "estimating one needs the 'gstat' package, which is not installed"
      } else {
        .msg("summarize_by_cell(): no fitted variogram supplied; estimating one.")
        # On the RESPONSE, not on OLS residuals.  The ..se_resp_* columns are
        # the SE of the cell mean as an estimate of the grand mean of the
        # response, so the correlation that has to be corrected for is the
        # response's own.  Passing predictor_vars here fitted a RESIDUAL
        # variogram -- the correlation of the part the predictors do not
        # explain, which is weaker -- and grand-mean coverage fell from 0.927
        # to 0.507 the moment a predictor was listed (which is the normal
        # thing to do, since that is also how per-cell predictor summaries are
        # requested).  A caller who wants the residual field's correlation
        # can pass it through `sac`.
        est <- try(estimate_sac_range(assigned_points_sf, response_var),
                   silent = TRUE)
        # Same test on the internally estimated fit: a rejected range means the
        # model behind it is not usable either.
        if (inherits(est, "try-error")) {
          no_fit <- sprintf("estimate_sac_range() failed: %s",
                            conditionMessage(attr(est, "condition")))
        } else if (no_structure(est)) {
          # No structure in the response: the pure nugget is the model.
          structure_note <- as.character(attr(est, "rejected_reason"))[1L]
          vgm_model <- attr(est, "variogram_model")
          vgm_crs   <- attr(est, "crs")
        } else if (is.na(suppressWarnings(as.numeric(est))[1L]) &&
                   !is.null(attr(est, "rejected_reason"))) {
          no_fit <- sprintf(paste0("the variogram estimated from `response_var` ",
                                   "reports no usable range (%s)"),
                            as.character(attr(est, "rejected_reason"))[1L])
        } else {
          vgm_model <- attr(est, "variogram_model")
          vgm_crs   <- attr(est, "crs")
          if (is.null(vgm_model))
            no_fit <- "estimate_sac_range() returned no fitted model"
        }
      }
    }
    cor_fn <- .vgm_correlation_fn(vgm_model)
    if (is.null(cor_fn) && .vgm_is_pure_nugget(vgm_model)) {
      # No structured component: distinct observations are uncorrelated and
      # a design effect of 1 is the exact answer, not a fallback.  It used to
      # be reported as "the supplied model could not be read", with the
      # fallback warning and deff_applied FALSE.
      .log_info(paste0("summarize_by_cell(): the variogram model is a pure ",
                       "nugget%s, so distinct observations are uncorrelated and ",
                       "the design effect is 1 in every cell."),
                if (is.null(structure_note)) ""
                else sprintf(" (the variogram's verdict: %s)", structure_note))
      cor_fn <- function(h) 0 * h
    }
    if (is.null(cor_fn)) {
      # Name the family when that is the reason: a user-built Matern or power
      # model used to be read as exponential without a word.
      other <- if (is.data.frame(vgm_model) && "model" %in% names(vgm_model))
        setdiff(unique(as.character(vgm_model$model)),
                c("Nug", names(.vgm_shape_fns))) else character(0)
      why <- if (length(other))
        sprintf(paste0("deff = \"variogram\" supports exponential, spherical and ",
                       "Gaussian variogram models (plus a nugget); the supplied ",
                       "model uses %s."),
                paste(sQuote(other, FALSE), collapse = ", "))
      else
        sprintf(paste0("deff = \"variogram\" requires a fitted variogram model; ",
                       "none was available (%s). Pass one via ",
                       "`sac = estimate_sac_range(...)`."),
                paste(c(set_aside,
                        if (is.null(no_fit)) "the supplied model could not be read"
                        else no_fit),
                      collapse = "; "))
      # The request cannot be met.  It used to fall back to deff = 1, the
      # uncorrected interval, under the same column names: the interval with
      # no correction at all for a caller who asked for one.  A range goes
      # unidentified often enough that this is an ordinary case, not a
      # corner.  Kish needs no variogram, so it is the fallback; the
      # uncorrected SE is returned only when Kish's own ICC is 0, and either
      # way the rows say which.
      use_vgm  <- FALSE
      use_kish <- TRUE
      deff     <- 1
      vgm_fallback <- why
      icc_est  <- .kish_iccs()
      resp_rho <- icc_est$resp
      pred_rho <- icc_est$pred
      icc_txt <- paste(c(
        if (length(resp_num)) sprintf("response '%s' %.3f", resp_num[[1L]], resp_rho),
        if (length(pred_num)) sprintf("predictor '%s' %.3f", pred_num, pred_rho)),
        collapse = ", ")
      if (max(resp_rho, pred_rho, 0) > 0)
        .warn_deff_fallback(paste0(
          "summarize_by_cell(): %s Falling back to deff = \"kish\" (ICC: %s): ",
          "the standard errors carry Kish's correction in place of the ",
          "variogram's, and `deff_method` reads \"kish (variogram fallback)\" ",
          "on every row."),
          why, icc_txt, data = list(requested = "variogram", used = "kish"))
      else
        .warn_deff_fallback(paste0(
          "summarize_by_cell(): %s Falling back to deff = \"kish\", whose ICC ",
          "is 0 (%s), so the standard errors are the uncorrected ones; ",
          "`deff_method` reads \"none (variogram fallback)\" on every row."),
          why,
          if (nzchar(icc_txt)) icc_txt else "no numeric column to estimate it from",
          data = list(requested = "variogram", used = "none"))
    } else {
      if (!is.null(set_aside))
        .warn_and_log(paste0("summarize_by_cell(): %s%s. It was set aside, ",
                             "and the design effect uses the variogram estimated ",
                             "from `response_var` instead."),
                      set_aside, set_aside_so)
      if (detrended_sac && length(resp_num))
        .warn_and_log(paste0("summarize_by_cell(): `sac` is the variogram of the ",
                             "residuals on predictors (detrend = \"%s\"), but the ",
                             "..se_resp_* columns estimate the grand mean of the ",
                             "response, whose own correlation is the one to correct ",
                             "for; a residual variogram is weaker, so those standard ",
                             "errors are understated. Pass a variogram of the ",
                             "response (estimate_sac_range() without predictor_vars) ",
                             "or sac = NULL."),
                      as.character(attr(sac, "detrend_method") %||% "ols")[1L])
      # st_coordinates() returns one row per VERTEX, so any non-POINT geometry
      # (this function accepts POLYGON and MULTIPOINT features) would misalign
      # coords_mat with `df` and feed the wrong points into every cell.
      # coerce_to_points() takes representative points, matching the guard in
      # estimate_sac_range().
      pts_for_deff <- assigned_points_sf
      if (!all(sf::st_geometry_type(pts_for_deff, by_geometry = TRUE) == "POINT"))
        pts_for_deff <- coerce_to_points(pts_for_deff, "auto")

      # The variogram range is a LENGTH in the CRS the variogram was fitted in
      # -- estimate_sac_range() always fits in a projected CRS, in metres.
      # Evaluating the correlation function at distances taken from
      # unprojected input (degrees) saturates every pair: measured on the same
      # points, cells and variogram, median SE 0.657 with UTM input against
      # 235.9 with EPSG:4326 input, a factor of 354.  Put the points in the
      # variogram's own CRS when it is recorded, else project them the way
      # estimate_sac_range() would have.
      vgm_crs_ok <- !is.null(vgm_crs) && inherits(vgm_crs, "crs") && !is.na(vgm_crs)
      pts_for_deff <- if (vgm_crs_ok) {
        .transform_or_stamp(pts_for_deff, vgm_crs, what = "assigned_points_sf",
                            caller = "summarize_by_cell")
      } else {
        ensure_projected(pts_for_deff)
      }
      coords_mat <- sf::st_coordinates(pts_for_deff)[, 1:2, drop = FALSE]
      # An empty point has no location to correlate.  It still counts towards
      # its cell's values, and .cor_stats_from_coords() leaves it out of the
      # correlation; it used to stop the call with "missing value where
      # TRUE/FALSE needed", while deff = "kish" took the same layer.
      unlocated <- !(is.finite(coords_mat[, 1]) & is.finite(coords_mat[, 2])) &
        !is.na(df[[id_col]])
      if (any(unlocated))
        .warn_and_log(paste0("summarize_by_cell(): %d point(s) have empty or ",
                             "non-finite coordinates; their values are ",
                             "summarised, but the variogram design effect of ",
                             "their cells is computed from the located points ",
                             "only."),
                      sum(unlocated))

      # Per-cell mean off-diagonal correlation, NOT a per-cell deff: the design
      # effect each column needs is 1 + (n_valid - 1) * rbar for ITS non-missing
      # count, and that is formed inside the ..se_ closures.
      vgm_stats <- .cell_cor_stats_variogram(coords_mat, df[[id_col]], cor_fn,
                                             max_n = deff_max_n)
      vgm_rbar <- vgm_stats$rbar
      vgm_deff <- .cell_deff_variogram(coords_mat, df[[id_col]], cor_fn,
                                       max_n = deff_max_n, r_bar = vgm_rbar)
      vgm_bits <- list(coords = coords_mat, cor_fn = cor_fn, max_n = deff_max_n,
                       df_ratio = vgm_stats$df_ratio)
      .msg(sprintf(
        "summarize_by_cell(): variogram deff across cells (at the cell row count): median %.3f, max %.3f",
        stats::median(vgm_deff, na.rm = TRUE), max(vgm_deff, na.rm = TRUE)
      ))
    }
  }

  # --- single grouped summarise for all columns ---
  grouped <- df |> dplyr::group_by(.data[[id_col]])

  has_resp <- length(resp_num) > 0L
  has_pred <- length(pred_num) > 0L

  if (has_resp && has_pred) {
    out <- grouped |>
      dplyr::summarise(
        ..n_rows = dplyr::n(),
        dplyr::across(dplyr::all_of(resp_num), .build_funs("resp_", rho = resp_rho, rbar = vgm_rbar, vgm = vgm_bits),
                      .names = "{.fn}{.col}"),
        dplyr::across(dplyr::all_of(pred_num), .build_funs("pred_", rho = pred_rho, rbar = vgm_rbar, vgm = vgm_bits),
                      .names = "{.fn}{.col}"),
        .groups = "drop"
      )
  } else if (has_resp) {
    out <- grouped |>
      dplyr::summarise(
        ..n_rows = dplyr::n(),
        dplyr::across(dplyr::all_of(resp_num), .build_funs("resp_", rho = resp_rho, rbar = vgm_rbar, vgm = vgm_bits),
                      .names = "{.fn}{.col}"),
        .groups = "drop"
      )
  } else if (has_pred) {
    out <- grouped |>
      dplyr::summarise(
        ..n_rows = dplyr::n(),
        dplyr::across(dplyr::all_of(pred_num), .build_funs("pred_", rho = pred_rho, rbar = vgm_rbar, vgm = vgm_bits),
                      .names = "{.fn}{.col}"),
        .groups = "drop"
      )
  } else {
    out <- grouped |>
      dplyr::summarise(..n_rows = dplyr::n(), .groups = "drop")
  }
  # dplyr::summarise() makes each newly created column visible to the
  # expressions that follow it, so naming the row count `n` here put it in
  # scope for the across() calls below: a response or predictor column
  # literally named `n` was silently summarised as the ROW COUNT (resp_mean_n
  # equal to n, sd and se NA), with no message.  Build it under a reserved
  # name and rename once nothing can shadow it.
  names(out)[names(out) == "..n_rows"] <- "n"

  # cell_weight is the EFFECTIVE sample size of the primary variable -- the
  # response when there is one, else the first predictor -- so it is built from
  # that variable's non-missing count per cell, not from the row count: a cell
  # of 10 rows with 3 finite responses and deff = 2 has 1.5 observations' worth
  # of information about the response, not 5.  (`n` still reports rows.)
  primary_col <- if (has_resp) resp_num[[1L]] else if (has_pred) pred_num[[1L]] else NULL
  n_valid_primary <- if (is.null(primary_col)) out$n else {
    # exclude = NULL keeps the group of rows with no cell ID (from
    # keep_unassigned = TRUE), which is a summary row too: tapply() dropped
    # it, and its cell_weight came out 0 beside n = 5 and a finite SE.
    cnt <- tapply(!is.na(df[[primary_col]]), factor(df[[id_col]], exclude = NULL), sum)
    v   <- as.numeric(cnt[match(as.character(out[[id_col]]), names(cnt))])
    v[is.na(v)] <- 0
    v
  }
  out$cell_weight <- n_valid_primary

  # Record the design effect used and adjust cell_weight for effective n
  # cell_weight uses the response ICC (primary outcome); falls back to
  # predictor ICC when no response variable was supplied.
  if (use_kish) {
    # The fallback to the PREDICTOR's ICC is for the case the comment above
    # describes -- no response variable at all.  Applying it whenever the
    # response's own ICC came out <= 0 meant that adding a predictor to the
    # call changed the response's regression weight (measured: by a factor of
    # 20), deflating a cell whose response is uncorrelated on the strength of a
    # predictor that happens to be clustered.  The response SE is already left
    # at deff = 1 in that case, as documented; cell_weight has to agree with it.
    # With no response the primary variable is the FIRST predictor, and the
    # weight is built from that predictor's own ICC.  From an ICC pooled over
    # every predictor listed, an unclustered primary predictor's weight fell
    # from 12 to 2.12 when a clustered one was listed after it.
    primary_rho <- if (has_resp) max(resp_rho, 0)
                   else if (has_pred) max(pred_rho[[1L]], 0)
                   else 0
    deff_per_cell <- pmax(1, 1 + (n_valid_primary - 1) * primary_rho)
    if (primary_rho > 0)
      out$cell_weight <- n_valid_primary / deff_per_cell
    # Recorded whenever ANY ICC corrected its column.  Keyed on the
    # primary variable alone, a clustered predictor beside an unclustered
    # response had its SEs inflated 11x on rows marked deff_applied = FALSE,
    # with no attribute.  `deff` stays the primary variable's (all 1 then).
    if (max(resp_rho, pred_rho, 0) > 0) {
      attr(out, "deff_applied") <- c(list(
        method   = "kish",
        icc_resp = if (has_resp) resp_rho else NA_real_,
        # One ICC per predictor column, named by column.
        icc_pred = if (has_pred) pred_rho else NA_real_,
        deff     = deff_per_cell
      ),
      # A "variogram" request that fell back here says so, and why.
      if (!is.null(vgm_fallback))
        list(requested = "variogram", fallback_reason = vgm_fallback))
    }
  } else if (use_vgm && !is.null(vgm_rbar)) {
    # The ..se_ closures have already applied the per-column design effect
    # (see .build_funs()); nothing is rescaled here.  What is recorded is the
    # design effect of the PRIMARY variable at its own non-missing count, which
    # is also what cell_weight is built from.
    # Same rule as the closures: rbar over the primary column's observed rows.
    rb <- vapply(seq_len(nrow(out)), function(i) {
      cell <- out[[id_col]][i]
      # `primary_col` is NULL when neither response_var nor predictor_vars was
      # supplied (see its definition above), and df[[NULL]] aborts with
      # "attempt to select less than one element in get1index" -- nowhere near
      # the cause.  The variogram design effect is a function of the cell's
      # COORDINATES and the fitted correlation, not of any column's values, so
      # with no value column there is nothing to filter on and every point in
      # the cell counts.
      rows <- if (is.null(primary_col))
        which(df[[id_col]] == cell)
      else
        which(df[[id_col]] == cell & !is.na(df[[primary_col]]))
      if (length(rows) == n_valid_primary[i] && length(rows) == sum(df[[id_col]] == cell, na.rm = TRUE)) {
        v <- vgm_rbar[as.character(cell)]
        if (is.finite(v)) unname(v) else 0
      } else {
        .rbar_from_coords(coords_mat[rows, , drop = FALSE], cor_fn, max_n = deff_max_n)
      }
    }, numeric(1))
    d_i <- pmin(pmax(1, 1 + (n_valid_primary - 1) * rb), pmax(n_valid_primary, 1))

    out$cell_weight <- n_valid_primary / d_i
    # `deff_rows` is the design effect at the cell's ROW count -- the vector
    # the log line above reduces to a median and a max -- beside `deff`,
    # which is at the primary variable's non-missing count.  The two agree
    # wherever that variable is complete.
    deff_rows <- if (is.null(vgm_deff)) NULL else
      unname(vgm_deff[as.character(out[[id_col]])])
    attr(out, "deff_applied") <- list(
      method    = "variogram",
      deff      = d_i,
      deff_rows = deff_rows,
      rbar      = rb,
      crs       = if (isTRUE(vgm_crs_ok)) vgm_crs else sf::st_crs(pts_for_deff),
      max_n     = deff_max_n
    )
  } else if (!use_kish && is.numeric(deff) && deff > 1) {
    out$cell_weight <- n_valid_primary / deff
    attr(out, "deff_applied") <- list(method = "fixed", deff = deff)
  }
  # The ICCs a Kish request estimated, whether or not either was positive
  # enough to apply: `deff_applied` is absent exactly when nothing was
  # applied, which is when a user most wants to see what the ICC came out
  # as.  NA for a variable type that had no numeric column.
  if (use_kish)
    attr(out, "icc") <- list(
      resp = if (has_resp) as.numeric(resp_rho) else NA_real_,
      # Named by predictor column: each has its own.
      pred = if (has_pred) pred_rho else NA_real_)
  # The same fact per row.  The attribute does not survive rbind(),
  # dplyr::bind_rows() or most dplyr verbs, so results combined across calls
  # (a simulation, a benchmark over many layers) could not tell a row whose
  # design effect fell back to 1 from one that was corrected, and a fallback
  # fires when the variogram is rejected -- often when correlation is
  # strongest.  Only when a design effect was requested, so the default
  # frame is unchanged.
  if (deff_requested) {
    applied <- !is.null(attr(out, "deff_applied"))
    out$deff_applied <- applied
    # Which correction the standard errors carry, by name, since the columns
    # keep their names whatever was applied: "variogram", "kish" or "fixed"
    # when the one asked for was applied, "none" when it was not (a refused
    # number, Kish ICCs of 0), and "kish (variogram fallback)" or
    # "none (variogram fallback)" when a "variogram" request could not be
    # met and Kish stood in for it.
    out$deff_method <- paste0(
      if (applied) attr(out, "deff_applied")$method else "none",
      if (!is.null(vgm_fallback)) " (variogram fallback)" else "")
  }

  # --- pooled values ---
  # Appended after every column above, so each of those keeps its name, its
  # values and its position, and before the join below adds the geometry.
  # Nothing here reads `agg_funs` or `deff`: a column is pooled on its own
  # scored rows, from their values, their cells and their coordinates.
  if (identical(estimate, "pooled") && (has_resp || has_pred)) {
    # The design effect the t interval carries is the one recorded above as
    # applied, not the one asked for: a refused number, and Kish ICCs of 0 for
    # every column, leave the default's numbers, and the table of intervals
    # then describes them as the default's.
    deff_carried <- attr(out, "deff_applied")
    pooled <- .pooled_columns(
      assigned_points_sf, df, id_col, out[[id_col]], resp_num, pred_num,
      obs_level = obs_level, conf_level = conf_level,
      deff_label = .deff_label(
        deff_carried,
        # Under Kish each column has its own ICC, and one of 0 leaves that
        # column's standard error as it is.
        kish_zero = if (use_kish)
          unique(c(if (has_resp && resp_rho <= 0) resp_num,
                   if (has_pred) pred_num[pred_rho <= 0]))))
    for (nm in names(pooled$columns)) out[[nm]] <- pooled$columns[[nm]]
    attr(out, "pooled")    <- pooled$record
    attr(out, "intervals") <- pooled$record$intervals
    .pooled_warnings(pooled$events, obs_level)
  }

  if (!is.null(cells_sf)) {
    if (!inherits(cells_sf, "sf")) {
      # A warning, not a log line: the documented return is an sf layer.
      .warn_and_log("summarize_by_cell(): `cells_sf` is not an sf object; returning a plain data frame without cell geometry.")
    } else {
      # Locate the matching ID column in cells_sf: the candidates, in the
      # order, that assign_features_to_polygons() reads the polygons' IDs
      # from, so the cells are read the way the points were labelled.  'id'
      # and 'grid_id' were missing here, so cells keyed by 'id' (common in a
      # shapefile) were assigned cleanly and then summarised to a plain table
      # with no geometry, and cells with 'id' and a differently numbered
      # 'cell_id' were joined on 'cell_id', putting most summaries on the
      # wrong polygons.
      cells_id_candidates <- unique(c(id_col, "poly_id", "polygon_id", "id",
                                      "cell_id", "grid_id"))
      cells_id_found <- cells_id_candidates[cells_id_candidates %in% names(cells_sf)]
      if (length(cells_id_found) == 0L) {
        msg <- sprintf(paste0("`cells_sf` has none of the ID columns %s, so the ",
                              "summaries cannot be joined to it"),
                       paste(sprintf("'%s'", cells_id_candidates), collapse = ", "))
        if (isTRUE(area))
          stop("summarize_by_cell(): `area = TRUE` needs the cells' areas, but ",
               msg, ". Rename the cells' ID column to '", id_col, "'.",
               call. = FALSE)
        .warn_and_log("summarize_by_cell(): %s; returning a plain data frame without cell geometry. Rename the cells' ID column to '%s'.",
                      msg, id_col)
      } else {
        cells_id <- cells_id_found[[1]]
        if (cells_id != id_col)
          .log_info(paste0("summarize_by_cell(): `cells_sf` has no '%s' column; ",
                           "joining the summaries on its '%s' column."),
                    id_col, cells_id)
        # Keep only geometry + id from cells to avoid column collisions
        cells_slim <- cells_sf[, cells_id, drop = FALSE]
        # A duplicated cell ID makes left_join() emit one summary row per
        # matching cell, so sum(n) exceeded the number of points assigned (39
        # for 30 points; 46 for 23 after st_cast(., "POLYGON") split a
        # multi-part cell and each piece kept its poly_id).  Nothing checked
        # it, and the inflated total is easy to mistake for a real count.
        if (anyDuplicated(cells_slim[[cells_id]]))
          .warn_and_log(paste0("summarize_by_cell(): `cells_sf$%s` has %d ",
                               "duplicated value(s), so each summary row is ",
                               "repeated once per matching cell and totals such ",
                               "as sum(n) will exceed the number of points. Give ",
                               "every cell its own ID (see ",
                               "ensure_stable_poly_id())."),
                        cells_id, sum(duplicated(cells_slim[[cells_id]])))
        if (cells_id != id_col) {
          names(cells_slim)[names(cells_slim) == cells_id] <- id_col
        }
        # Join on the native type when both sides already agree.  Coercing
        # unconditionally made the returned ID type depend on an unrelated
        # argument and turned integer IDs into "1", "10", "2", ...
        if (!identical(class(out[[id_col]]), class(cells_slim[[id_col]]))) {
          .log_warn(
            "summarize_by_cell(): cells_sf$%s is <%s> but the summarised IDs are <%s>; coercing both to character to join, so the returned '%s' column is character.",
            id_col, paste(class(cells_slim[[id_col]]), collapse = "/"),
            paste(class(out[[id_col]]), collapse = "/"), id_col
          )
          out[[id_col]] <- .id_as_character(out[[id_col]])
          cells_slim[[id_col]] <- .id_as_character(cells_slim[[id_col]])
        }
        # A summarised cell that matches no polygon is dropped by the join
        # below, with its points.  That was silent -- the coercion above lost
        # every cell whose double ID R prints in scientific notation, 11 of 30
        # points in one check -- and it is also what a `cells_sf` other than
        # the layer the points were assigned to looks like.
        ids_out   <- out[[id_col]]
        unmatched <- !is.na(ids_out) & !(ids_out %in% cells_slim[[id_col]])
        if (any(unmatched))
          .warn_and_log(paste0("summarize_by_cell(): %d of the %d summarised cell ",
                               "ID(s) (%s) match no `cells_sf$%s`, so they and the ",
                               "%d point(s) in them are not in the result. Check ",
                               "that `cells_sf` is the layer the points were ",
                               "assigned to."),
                        sum(unmatched), sum(!is.na(ids_out)),
                        paste(utils::head(.id_as_character(ids_out[unmatched]), 5L),
                              collapse = ", "),
                        cells_id, sum(out$n[unmatched]))

        # dplyr::left_join() rebuilds attributes from the `x` template, which
        # silently drops "deff_applied".  Save it, then re-attach it, mapping
        # the per-cell vector onto the joined row order (cells with no
        # observations get NA).
        deff_attr   <- attr(out, "deff_applied")
        icc_attr    <- attr(out, "icc")
        pooled_attr <- attr(out, "pooled")
        pre_join_id <- as.character(out[[id_col]])

        out <- dplyr::left_join(cells_slim, out, by = id_col)
        if (!is.null(icc_attr)) attr(out, "icc") <- icc_attr
        # The record of the pooled columns and its table of intervals, which
        # the join drops as it drops "icc".
        if (!is.null(pooled_attr)) {
          attr(out, "pooled")    <- pooled_attr
          attr(out, "intervals") <- pooled_attr$intervals
        }

        if (isTRUE(area)) {
          # Planar areas in the cells' own (verified equal-area) CRS, as
          # plain numbers in its squared units; the density is points per
          # unit area.  Cells with no observations have n = NA after the
          # join and get a density of NA, not 0: no point was counted there
          # because the layer had none, which is not the same finding.
          cell_area <- suppressWarnings(as.numeric(sf::st_area(sf::st_geometry(out))))
          cell_area[!is.finite(cell_area) | cell_area <= 0] <- NA_real_
          out$cell_area  <- cell_area
          out$n_per_area <- as.numeric(out$n) / cell_area
          attr(out, "area_error") <- area_err
        }

        if (!is.null(deff_attr)) {
          # EVERY per-cell vector has to follow the join, not just $deff.
          # Realigning $deff alone left $rbar in pre-join order and pre-join
          # length, so deff[i] and rbar[i] described different cells and the
          # identity deff = 1 + (n-1) * rbar no longer held row-wise.
          # A fixed deff is one number, not a per-cell vector: with exactly
          # one cell summarised its length matched, and deff = 2 came back
          # as c(2, NA, NA, ...).
          flds <- if (identical(deff_attr$method, "fixed")) character(0)
                  else c("deff", "deff_rows", "rbar")
          for (fld in flds) {
            v <- deff_attr[[fld]]
            if (!is.null(v) && length(v) == length(pre_join_id)) {
              lookup <- stats::setNames(v, pre_join_id)
              deff_attr[[fld]] <- unname(lookup[as.character(out[[id_col]])])
            }
          }
          attr(out, "deff_applied") <- deff_attr
        }
      }
    }
  }

  out
}


# as.character() writes a double in scientific notation from 1e5 on
# ("1e+05") but an integer or a string in full ("100000"), so a join across
# ID classes lost every cell whose double ID R prints that way.  Whole
# numbers are written out digit by digit, and anything else as
# as.character() writes it (format(scientific = FALSE) would pad every
# element to a common width).
.id_as_character <- function(x) {
  if (!is.double(x) || is.object(x)) return(as.character(x))
  out   <- as.character(x)
  whole <- is.finite(x) & x == trunc(x) & abs(x) < 2^53
  out[whole] <- sprintf("%.0f", x[whole] + 0)   # + 0 turns -0 into 0
  out
}


# -----------------------------------------------------------------------------
# Pooled values in summarize_by_cell()
# -----------------------------------------------------------------------------

# Beyond this a coordinate is not read as a location: the squared distance
# between two such points is not a finite number, and the two neighbour
# searches then stop or disagree.
.POOL_XY_MAX <- 1e150

#' The coordinates the pooled values are computed from
#'
#' One row per row of the layer, in a projected CRS.  Features that are not
#' points are reduced to one point each first, and the points are then
#' projected, in that order, which is the order \code{resolution_profile()}
#' takes: lon/lat coordinates, and CRS-less ones that look like degrees, go
#' through \code{ensure_projected()}, so degrees are never read as distances.
#' The projection is chosen from every row of the layer, as the profile
#' chooses it.
#'
#' A feature with an empty geometry has \code{NA} coordinates, and so has one
#' with a coordinate beyond \code{.POOL_XY_MAX}: both are rows without a
#' location to the caller.  The warning sf raises when it takes a point on
#' the surface of a lon/lat polygon is not passed on: the point is one of the
#' polygon's, and it is projected before any distance is taken.
#'
#' @param layer The assigned layer.
#' @return \code{NULL} when the layer has no geometry (a plain data frame);
#'   otherwise a list with \code{xy}, a two-column matrix, and \code{crs},
#'   the CRS its distances are measured in.
#' @keywords internal
#' @noRd
.pooled_xy <- function(layer) {
  if (!inherits(layer, "sf")) return(NULL)
  if (nrow(layer) == 0L)
    return(list(xy = matrix(numeric(0), 0L, 2L), crs = sf::st_crs(layer)))
  if (!all(sf::st_geometry_type(layer, by_geometry = TRUE) == "POINT"))
    layer <- withCallingHandlers(
      coerce_to_points(layer, "auto"),
      warning = function(w)
        if (grepl("st_point_on_surface may not give correct results",
                  conditionMessage(w), fixed = TRUE))
          invokeRestart("muffleWarning"))
  layer <- ensure_projected(layer)
  xy <- unname(sf::st_coordinates(layer)[, 1:2, drop = FALSE])
  far <- abs(xy[, 1L]) > .POOL_XY_MAX | abs(xy[, 2L]) > .POOL_XY_MAX
  xy[which(far), ] <- NA_real_
  list(xy = xy, crs = sf::st_crs(layer))
}


#' The design effect the t interval carries, as a call writes it
#'
#' Read from what \code{summarize_by_cell()} recorded as applied, so that the
#' table of intervals describes the standard errors the frame holds: a
#' request that was refused, or that corrected no column, reads \code{1}.
#'
#' @param applied The \code{"deff_applied"} attribute of the summary, or
#'   \code{NULL} when no correction was applied.
#' @param kish_zero Under Kish's correction, the columns whose ICC is 0 and
#'   whose standard errors are therefore the uncorrected ones.
#' @return One string for \code{.intervals_table()}: \code{1},
#'   \code{"kish"} or \code{"variogram"} with their quotes, or the number.
#'   It reads \code{1} only when nothing was applied: the table takes the
#'   words of its row from that.
#' @keywords internal
#' @noRd
.deff_label <- function(applied, kish_zero = character(0)) {
  if (is.null(applied)) return("1")
  switch(applied$method,
    kish = paste0(
      "\"kish\"",
      if (length(kish_zero))
        sprintf("; 1 for %s, whose ICC is 0", .cols_text(kish_zero))),
    variogram = "\"variogram\"",
    {
      # A number applied is above 1, and is written so that it does not read 1.
      lab <- format(applied$deff, digits = 15)
      if (identical(lab, "1")) sprintf("%.17g", applied$deff) else lab
    })
}


#' Pooled columns, their record and the events to warn about
#'
#' What \code{summarize_by_cell(estimate = "pooled")} adds to its grouped
#' summary.  Every numeric column is fitted once, on its own scored rows
#' (\code{.pool_column()}), the response first and then the predictors in the
#' order given; a column named as both is fitted once and returned under both
#' prefixes.  The centroids and neighbours of the cells depend on which rows
#' are scored and not on their values, so a column whose scored rows are
#' those of an earlier column is handed that column's.
#'
#' @param layer The assigned layer as the caller passed it.
#' @param df Its attribute table, infinite values already set to \code{NA}.
#' @param id_col Name of the ID column of \code{df}.
#' @param out_ids The IDs of the summary's rows, in its order.
#' @param resp_num,pred_num The numeric response and predictor columns.
#' @param obs_level,conf_level The two levels of the call, or \code{NULL}.
#' @param deff_label The design effect of the t interval, as
#'   \code{.intervals_table()} takes it.
#' @return A list: \code{columns}, the new columns by name, in the order they
#'   are appended; \code{record}, the list of class \code{spatialkit_pooled};
#'   and \code{events}, what \code{.pooled_warnings()} reads.
#' @keywords internal
#' @noRd
.pooled_columns <- function(layer, df, id_col, out_ids, resp_num, pred_num,
                            obs_level, conf_level, deff_label) {
  use_fnn <- .has_backend("FNN")
  id   <- df[[id_col]]
  loc  <- .pooled_xy(layer)
  xy   <- loc$xy
  cols <- unique(c(resp_num, pred_num))
  located <- if (is.null(xy)) NULL else is.finite(xy[, 1L]) & is.finite(xy[, 2L])

  fits <- stats::setNames(vector("list", length(cols)), cols)
  kept <- list()                      # the scored rows and geometry of earlier columns
  unlocated <- logical(nrow(df))      # rows with a value, a cell and no location
  for (v in cols) {
    geo <- NULL
    if (!is.null(xy)) {
      valued    <- is.finite(df[[v]]) & !is.na(id)
      scored    <- valued & located
      unlocated <- unlocated | (valued & !located)
      for (k in kept) if (identical(k$mask, scored)) { geo <- k$geo; break }
    }
    res <- .pool_column(df[[v]], id, xy, obs_level, want = "pooled",
                        use_fnn = use_fnn, geo = geo)
    if (is.null(geo) && !is.null(res$geo))
      kept[[length(kept) + 1L]] <- list(mask = res$mask, geo = res$geo)
    # Kept of a column: what the frame, the record and the warnings read.  Its
    # fit, with the residuals, and its own copy of the mask and the geometry
    # are let go here, so that many columns do not hold them all at once.
    fits[[v]] <- res[c("ids", "value", "w", "target", "obs_lo", "obs_hi", "record",
                       "work", "limit")]
  }

  columns <- list()
  records <- list()
  roles <- rbind(
    data.frame(prefix = rep("resp_", length(resp_num)), column = resp_num,
               role = rep("response", length(resp_num)), stringsAsFactors = FALSE),
    data.frame(prefix = rep("pred_", length(unique(pred_num))), column = unique(pred_num),
               role = rep("predictor", length(unique(pred_num))), stringsAsFactors = FALSE))
  for (i in seq_len(nrow(roles))) {
    v   <- roles$column[i]
    nm  <- paste0(roles$prefix[i], v)
    res <- fits[[v]]
    at  <- match(out_ids, res$ids)
    columns[[paste0("..pooled_", nm)]]        <- res$value[at]
    columns[[paste0("..pooled_w_", nm)]]      <- res$w[at]
    columns[[paste0("..pooled_target_", nm)]] <- res$target[at]
    if (!is.null(obs_level)) {
      columns[[paste0("..obs_lo_", nm)]] <- res$obs_lo[at]
      columns[[paste0("..obs_hi_", nm)]] <- res$obs_hi[at]
    }
    rec <- res$record
    rec$column <- v
    rec$role   <- roles$role[i]
    records[[nm]] <- rec
  }

  record <- structure(
    list(toward     = "neighbours",
         neighbours = .POOL_NEIGHBOURS,
         crs        = if (is.null(loc)) sf::NA_crs_ else loc$crs,
         obs_level  = obs_level,
         conf_level = conf_level,
         intervals  = .intervals_table(obs_level, conf_level, deff_label),
         columns    = records),
    class = "spatialkit_pooled")

  events <- list(
    records = lapply(fits, function(r) r$record),
    w_min   = vapply(fits, function(r)
      if (length(r$w) && !all(is.na(r$w))) min(r$w, na.rm = TRUE) else NA_real_, numeric(1)),
    work    = vapply(fits, function(r) as.numeric(r$work), numeric(1)),
    limit   = vapply(fits, function(r) as.numeric(r$limit), numeric(1)),
    skipped_rows  = sum(unlocated),
    skipped_cells = length(unique(id[unlocated])))
  list(columns = columns, record = record, events = events)
}


#' Raise what happened to the pooled columns, once per call and reason
#'
#' Each event is raised once, for all the columns it applies to: the text
#' names them (\code{.cols_text()}) and the condition carries them in
#' \code{columns}, in the order of the frame.  A field "by column" is a
#' vector named by those columns.  The numbers of each column are in
#' \code{attr(, "pooled")}.
#'
#' @param ev The \code{events} of \code{.pooled_columns()}.
#' @param obs_level The level of the interval for one observation, or
#'   \code{NULL}.
#' @return \code{NULL}, invisibly.
#' @keywords internal
#' @noRd
.pooled_warnings <- function(ev, obs_level) {
  rec  <- ev$records
  cols <- names(rec)
  if (!length(cols)) return(invisible(NULL))
  of <- function(field, as = as.numeric)
    stats::setNames(as(unlist(lapply(rec, function(r) r[[field]]))), cols)
  status <- of("status", as.character)
  n      <- of("n", as.integer)
  cells  <- of("cells", as.integer)
  df     <- of("df", as.integer)
  one    <- function(w) length(w) == 1L
  listed <- function(x) paste(format(x, trim = TRUE), collapse = ", ")

  # ---- the pooled columns are NA --------------------------------------------
  w <- cols[status == "no coordinates"]
  if (length(w))
    .warn_classed("spatialkit_pool_unavailable", paste0(
      "summarize_by_cell(): estimate = \"pooled\" needs the observations' locations, ",
      "and `assigned_points_sf` carries no geometry: the pooled columns are NA. The ",
      "plain columns are unaffected."),
      data = list(reason = "no coordinates", columns = w))
  w <- cols[status == "no located value"]
  if (length(w))
    .warn_classed("spatialkit_pool_unavailable", paste0(
      "summarize_by_cell(): the pooled columns of %s are NA: no row has a value of ",
      "%s, a cell and a location."),
      .cols_text(w), if (one(w)) "it" else "any of them",
      data = list(reason = "no located value", columns = w))
  w <- cols[status == "no cell holds two values"]
  if (length(w))
    .warn_classed("spatialkit_pool_unavailable", paste0(
      "summarize_by_cell(): the pooled columns of %s are NA: %s, so the spread of ",
      "values inside a cell cannot be estimated and nothing can set a weight. The ",
      "plain columns are unaffected."),
      .cols_text(w),
      if (one(w)) "no cell holds two values of it"
      else "for none of them does any cell hold two values",
      data = list(reason = "no cell holds two values", columns = w,
                  n = n[w], cells = cells[w]))
  w <- cols[status == "work limit"]
  if (length(w)) {
    top <- w[which.max(cells[w])]             # the column with the most cells
    .warn_classed("spatialkit_pool_unavailable", paste0(
      "summarize_by_cell(): the pooled columns of %s are NA: finding each cell's ",
      "neighbours among %d cells without the FNN package would compare %.2g pairs, ",
      "above the limit of %.2g (options(spatialkit.exact_pairs_max = )). Install FNN: ",
      "the result is the same. The plain columns are unaffected."),
      .cols_text(w), cells[[top]], ev$work[[top]], ev$limit[[top]],
      data = list(reason = "work limit", columns = w, work = ev$work[[top]],
                  limit = ev$limit[[top]], cells = cells[w]))
  }

  # ---- the pooled values carry no pooling ------------------------------------
  w <- cols[status == "one cell"]
  if (length(w))
    .warn_classed("spatialkit_pool_degenerate", paste0(
      "summarize_by_cell(): the values of %s lie in one cell, so there is nothing to ",
      "pool toward: the pooled value is the plain mean and its weight is 1."),
      .cols_text(w),
      data = list(reason = "one cell", columns = w, n = n[w]))
  w <- cols[status == "no spread inside any cell"]
  if (length(w))
    .warn_classed("spatialkit_pool_degenerate", paste0(
      "summarize_by_cell(): the values of %s do not vary inside any cell, so every ",
      "own weight is 1 and the pooled value is the plain mean. attr(, \"pooled\") says ",
      "which columns are constant."),
      .cols_text(w),
      data = list(reason = "no spread inside any cell", columns = w,
                  constant = of("constant", as.logical)[w]))
  w <- cols[status == "no signal between cells"]
  if (length(w))
    .warn_classed("spatialkit_pool_degenerate", paste0(
      "summarize_by_cell(): every own weight of %s is 0: the cell means differ from ",
      "their neighbours' means by no more than the spread inside the cells explains. ",
      "Each pooled value is then the mean of the cell's nearest cells and holds none ",
      "of the cell's own values."),
      .cols_text(w),
      data = list(reason = "no signal between cells", columns = w,
                  sigma2 = of("sigma2")[w]))

  # ---- values were computed, and rest on little -------------------------------
  # The columns whose weights were estimated.  Where no value varies inside a
  # cell the weights are set to 1 and `sigma2` to 0: nothing rests on the
  # count of cells or of values there, and nothing is said of either.
  fitted <- status %in% c("ok", "no signal between cells")
  w <- cols[fitted & cells >= 2L & cells <= 9L]
  if (length(w))
    .warn_classed("spatialkit_pool_few_cells", paste0(
      "summarize_by_cell(): %s %s values in fewer than 10 cells (%s), so each cell's ",
      "target is the mean of all the other cells, not of cells near it. A cell with a ",
      "small own weight then gets a pooled value on the far side of the overall mean ",
      "from its plain mean, and two cells can change order. The `..pooled_w_*` columns ",
      "hold the weights; estimate = \"mean\" gives the plain means."),
      .cols_text(w), if (one(w)) "has" else "have", listed(cells[w]),
      data = list(columns = w, cells = cells[w], own_weight = ev$w_min[w]))
  w <- cols[fitted & cells >= 2L & df >= 1L & df < .POOL_MIN_DF]
  if (length(w))
    .warn_classed("spatialkit_pool_few_df", paste0(
      "summarize_by_cell(): the spread of %s inside a cell is estimated from fewer ",
      "than %d values beyond one per cell (%s), and every own weight rests on it: the ",
      "pooled values and their weights are unstable."),
      .cols_text(w), .POOL_MIN_DF, listed(df[w]),
      data = list(columns = w, needed = .POOL_MIN_DF, df = df[w], n = n[w],
                  cells = cells[w]))

  # ---- the cells do no better than one mean -----------------------------------
  # The rule has two tests (.no_structure()), and the text says which failed.
  # `worse`: one mean predicts a left-out row as well as the pooled values do.
  # `untold`: the pooled values predict it better, and the F test does not
  # tell the cell means from one mean, which is what it returns where few
  # values beyond one per cell stand behind it.
  structured <- of("structured", as.logical)
  w <- cols[!is.na(structured) & !structured]
  if (length(w)) {
    loo <- of("loo"); loo_one <- of("loo_one_cell"); p <- of("structure_p")
    own <- of("own_weight")
    beats  <- .beats_one_mean(loo[w], loo_one[w])
    worse  <- w[!beats]
    untold <- w[beats]
    says <- sprintf(paste0("leave-one-out error %s against %s; F test of the cells ",
                           "against one mean, p = %s; median own weight %.2f"),
                    vapply(loo[w], function(v) format(signif(v, 4)), character(1)),
                    vapply(loo_one[w], function(v) format(signif(v, 4)), character(1)),
                    vapply(p[w], function(v) format(signif(v, 2)), character(1)),
                    own[w])
    # One column: its numbers stand in the sentence.  Several: they follow the
    # sentences, which are what the reader has to see and R prints a warning's
    # first getOption("warning.length") characters only.
    inside <- if (one(w)) sprintf(" (%s)", says) else ""
    text <- c(
      if (length(worse)) sprintf(paste0(
        "the pooled values of %s predict a left-out observation no better than one ",
        "mean for all rows%s, so %s no spatial structure these cells can use. The ",
        "values are returned. A map of them looks smooth because each is mostly its ",
        "neighbours' mean, not because the response varies smoothly."),
        .cols_text(worse), inside, if (one(worse)) "it shows" else "they show"),
      if (length(untold)) sprintf(paste0(
        "%s pooled values of %s predict a left-out observation better than one mean ",
        "for all rows, but the F test does not tell the cell means from one mean%s, so ",
        "no spatial structure these cells can use is shown%s"),
        if (length(worse)) "The" else "the", .cols_text(untold), inside,
        if (length(worse)) sprintf(" for %s either.", if (one(untold)) "it" else "them")
        else ". The values are returned."),
      if (!one(w)) .each_text(sprintf("'%s': %s.", w, says)))
    .warn_classed("spatialkit_no_structure", "summarize_by_cell(): %s",
      paste(text, collapse = " "),
      data = list(criterion = "loo", n = cells[w], one_cell = loo_one[w], columns = w,
                  estimate = "pooled", loo = loo[w], structure_p = p[w],
                  beats = stats::setNames(beats, w)))
  }

  # ---- rows left out, and what the interval rests on --------------------------
  if (ev$skipped_rows > 0L)
    .warn_classed("spatialkit_pool_rows_skipped", paste0(
      "summarize_by_cell(): %d row(s) in %d cell(s) have empty or non-finite ",
      "coordinates. They count in `n` and the plain columns and take no part in the ",
      "pooled columns."),
      ev$skipped_rows, ev$skipped_cells,
      data = list(n = ev$skipped_rows, cells = ev$skipped_cells,
                  columns = cols[of("skipped", as.integer) > 0L]))
  if (!is.null(obs_level)) {
    # The columns whose interval is read from residuals: both offsets are
    # numbers, and the column is not constant.  The interval of a constant
    # column has no width because its residuals are set to 0; it rests on no
    # count of them and says nothing of where a further observation lies.
    has <- vapply(rec, function(r) !is.null(r$offsets) && all(is.finite(r$offsets)) &&
                    !isTRUE(r$constant), logical(1))
    repeated <- of("repeated", as.integer)
    w <- cols[has & repeated > 0L]
    if (length(w))
      .warn_classed("spatialkit_repeated_locations", paste0(
        "summarize_by_cell(): rows with a value share a location with another row ",
        "(%s). The leave-one-out error leaves one row out and keeps the others at its ",
        "location, so it, and an interval built from it, describes a further ",
        "observation at a place already observed, not an observation at a new place."),
        .each_text(sprintf("%d of the %d rows of '%s'", repeated[w], n[w], w), sep = "; ",
                   more = "and %d more"),
        data = list(columns = w, n = repeated[w], total = n[w]))
    # The count at which five residuals are expected beyond each end.  It is a
    # whole number wherever it fits one, and a double above that.
    p_lo   <- round((1 - obs_level) / 2, 12)
    needed <- ceiling(round(.FEW_RESIDUALS / p_lo, 6))
    if (needed <= .Machine$integer.max) needed <- as.integer(needed)
    w <- cols[has & n * p_lo < .FEW_RESIDUALS]
    if (length(w))
      .warn_classed("spatialkit_few_residuals", paste0(
        "summarize_by_cell(): the %s%% interval of %s rests on fewer than %s ",
        "leave-one-out residuals (%s), the number at which five are expected beyond ",
        "each end; its ends are unstable."),
        .level_text(obs_level), .cols_text(w), format(needed, scientific = FALSE),
        listed(n[w]),
        data = list(columns = w, needed = needed, level = obs_level, n = n[w]))
  }
  invisible(NULL)
}


#' Whether the pooled values predict a left-out row better than one mean
#'
#' The first of the two tests of \code{.no_structure()}, with the margin the
#' summary hands that rule.  The warning and the print read it to say which
#' test a column failed.
#'
#' @param loo,loo_one The leave-one-out errors of the pooled values and of one
#'   mean for all rows.
#' @return \code{TRUE} where one mean's error is above the pooled values' by
#'   more than \code{.STRUCTURE_TOL}.
#' @keywords internal
#' @noRd
.beats_one_mean <- function(loo, loo_one) !(loo_one <= loo * (1 + .STRUCTURE_TOL))


#' The numbers of each column, for at most six of them
#'
#' A warning that gives numbers column by column gives them for the first six
#' columns, as \code{.cols_text()} names six: the text of a warning is printed
#' to \code{getOption("warning.length")} characters, and the numbers of every
#' column are in \code{attr(, "pooled")}.
#'
#' @param pieces One string per column.
#' @param sep What joins them.
#' @param more What stands for the columns left out, with \code{\%d} for
#'   their number.
#' @return One string.
#' @keywords internal
#' @noRd
.each_text <- function(pieces, sep = " ", more = "And %d more, in attr(, \"pooled\").") {
  k <- length(pieces)
  if (k <= 6L) return(paste(pieces, collapse = sep))
  paste(c(pieces[1:6], sprintf(more, k - 6L)), collapse = sep)
}


#' A level as a percentage, with the digits it was given with
#'
#' @param level A level strictly between 0 and 1.
#' @return \code{"90"} for 0.9 and \code{"99.9999999"} for 0.999999999;
#'   \code{format()} alone writes seven digits and reads the second as 100.
#' @keywords internal
#' @noRd
.level_text <- function(level) format(100 * level, digits = 15)


#' @export
print.spatialkit_pooled <- function(x, ...) {
  label  <- function(name) sprintf("  %-12s: ", name)
  indent <- strrep(" ", 16L)
  sci    <- function(v) sprintf("%.2e", v)
  pct    <- function(level) paste0(.level_text(level), "%")
  count  <- function(k, what) sprintf("%d %s%s", as.integer(k), what, if (k == 1L) "" else "s")
  # The words beside an interval's columns, read from what the table says the
  # interval covers.
  target <- function(covers)
    if (startsWith(covers, "one further observation")) "ONE OBSERVATION"
    else if (grepl("population mean", covers, fixed = TRUE)) "THE POPULATION MEAN"
    else "THE CELL'S MEAN"
  offset <- function(o) sprintf("value %s %s", if (o < 0) "-" else "+",
                                format(signif(abs(o), 3)))
  iv    <- x$intervals
  shown <- rep(FALSE, NROW(iv))
  obs   <- which(iv$centre == "pooled value")

  for (nm in names(x$columns)) {
    r <- x$columns[[nm]]
    cat(sprintf("Cell values: %s, %s in %s\n", r$column, count(r$n, "row"),
                count(r$cells, "cell")))
    if (!identical(r$status, "ok")) {
      cat(label("status"), r$status, "\n", sep = "")
      next
    }
    cat(label("value"), "..pooled_", nm, ", pooled toward ",
        if (r$cells <= 9L) sprintf("all the other cells (%d of them)", r$neighbours)
        else sprintf("the %d nearest cells (centroids of the rows)", r$neighbours),
        "\n", sep = "")
    cat(label("variances"), "inside cells ", sci(r$sigma2), " (from ",
        count(r$df, "value"), " beyond one per cell);\n", indent,
        "between a cell and its target ", sci(r$tau2), "\n", sep = "")
    cat(label("own weight"), sprintf("median %.2f (..pooled_w_%s)\n", r$own_weight, nm),
        sep = "")
    cat(label("alone"), sprintf("%.1f%% of the rows are the only one in their cell\n",
                                100 * r$share_alone), sep = "")
    cat(label("loo"), sci(r$loo),
        " (leave-one-out mean squared error of the pooled value);\n", indent,
        "one mean for all rows: ", sci(r$loo_one_cell),
        # Which of the rule's two tests failed (.no_structure()).
        if (isFALSE(r$structured)) {
          if (.beats_one_mean(r$loo, r$loo_one_cell))
            sprintf(": the F test does not tell these cells from one mean (p = %s)",
                    format(signif(r$structure_p, 2)))
          else ": these cells do no better than one mean"
        },
        "\n", sep = "")
    if (length(obs) && !is.null(r$offsets)) {
      cat(label("interval"), pct(iv$level[obs]), " for ", target(iv$covers[obs]),
          sprintf(" (..obs_lo_%s, ..obs_hi_%s):\n", nm, nm), indent,
          if (all(is.finite(r$offsets)))
            sprintf("%s to %s, quantiles of %s", offset(r$offsets[1L]),
                    offset(r$offsets[2L]), count(r$n, "leave-one-out residual"))
          else "not available",
          "\n", sep = "")
      shown[obs] <- TRUE
    }
  }

  rest <- which(!shown)
  if (length(rest)) {
    cat("Also in the frame\n")
    for (i in rest) {
      covers <- strsplit(iv$covers[i], "; ", fixed = TRUE)[[1L]]
      head   <- sprintf("  %s : ", iv$columns[i])
      cat(head, pct(iv$level[i]), " for ", target(iv$covers[i]), " (", iv$method[i],
          "), centred on the ", iv$centre[i],
          if (length(covers) > 1L)
            paste0(";\n", strrep(" ", nchar(head)), paste(covers[-1L], collapse = "; ")),
          "\n", sep = "")
    }
  }
  cat("Away from the observations no interval here is calibrated.\n")
  invisible(x)
}
