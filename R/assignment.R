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
#' areas are planar. Cells in lon/lat whose record of their count says they
#' are the cells a profile scored, and holds its centres in a projected CRS
#' (the Voronoi cells [build_tessellation()] builds on the centres a profile
#' scored and returns in a geographic `crs`; see Value), were drawn in the
#' record's CRS, and the join runs there. A copy of `features_sf` is
#' transformed for it, and the features come back with the coordinates they
#' arrived with. Otherwise (other polygons in lon/lat, such as a hex or
#' square lattice sized by a selection or the cells of centres refitted to
#' every point, or polygons without a CRS) the join runs in the CRS of
#' `features_sf`. In longitude and latitude, with s2 on (`sf::sf_use_s2()`,
#' sf's default), sf joins on the sphere with s2, which refuses a polygon
#' that is not valid there, such as one with two corners it reads as one
#' point or with edges that cross; cells returned in lon/lat can hold two
#' corners that differ only by rounding and fall on one point in degrees
#' (see `crs` in [create_voronoi_polygons()]). A join by `predicate` then
#' stops with an error of class `"spatialkit_cells_degenerate"` that gives
#' the rows of the polygons s2 refuses (the condition carries them as
#' `rows`), names the cause and says to transform `polygons_sf` to a
#' projected CRS, where the join runs. The largest-overlap join repairs such
#' polygons first (see `largest`), and stops with that error when s2 still
#' refuses a polygon it repaired, such as one whose hole crosses its outer
#' ring.
#' A polygon of `features_sf` that s2 refuses stops the join by `predicate`
#' with an error of class `"spatialkit_features_degenerate"`, with the same
#' fields, that gives its rows, quotes s2's reason and says to transform
#' `features_sf` to a projected CRS, where the join runs; so does one s2
#' refuses even after the repair of the largest-overlap join. When that
#' repair makes s2 accept the polygons of both layers and s2 then refuses a
#' piece of an overlap of the two layers (a figure eight, one vertex listed
#' twice, is repaired into two parts that touch at a point, and the piece can
#' hold that point twice), the join stops with the error for the layer whose
#' polygons s2 refused as they arrived; with s2 on (sf's default), the join
#' by `predicate` meets those polygons too. When s2 refused no polygon of
#' either layer as it arrived, the error says to pass `largest = FALSE` or to
#' project both layers. Whether s2 refuses such a piece follows the last bits
#' of the vertices s2 computes for it, which differ between platforms and s2
#' versions: a figure eight pinched on the edge of a polygon of the other
#' layer can give a piece that s2 refuses on macOS on arm64 and accepts on
#' Linux on x86-64. Where s2 accepts the piece, the join runs and gives each
#' feature the polygon of its largest overlap.
#' In longitude and latitude the largest-overlap join
#' runs on s2, and the areas `tie_break` compares are measured on the sphere
#' with s2, whatever [sf::sf_use_s2()] says, as the package's other areas in
#' degrees are; the join by `predicate` runs with the engine the session
#' names. With s2 off that engine is GEOS, which reads these polygons in
#' degrees, so the largest-overlap join's error for a polygon s2 refuses also
#' offers `largest = FALSE`.
#'
#' sf reads a layer in the dimension of its first geometry, so s2, GEOS and
#' [sf::st_transform()] read a layer that mixes XY and XYZ (or XYM)
#' geometries wrongly; `rbind()` of an XY and an XYZ layer gives one. When
#' either layer mixes dimensions, the join runs on a copy of it without Z and
#' M, and the rows come back with the geometry they arrived with, Z and M
#' included. The same holds, with no warning, for a layer whose geometries
#' carry M (XYM or XYZM), which GEOS does not read: the join is defined in
#' x and y. Features without a CRS that mix dimensions and are moved from
#' lon/lat into the polygons' CRS come back moved and in XY, with a warning
#' of class `"spatialkit_zm_dropped"`, unless their only geometries of
#' another dimension than the first are empty and not points. Features
#' without a CRS in one dimension that carry M come back with their own
#' geometry: with the CRS they were given, or moved from lon/lat with
#' [sf::st_transform()], which keeps their M.
#'
#' Layers so far from the origin of their coordinates that doubles there are
#' spaced more than a millionth of the extent of the two layers apart (0.125
#' at 1e15, beside a layer a few kilometres across) are joined with their
#' edges and points held to that spacing, so a point that lies within it of
#' an edge near the origin can fall in the polygon beside it. The join runs,
#' and a warning of class `"spatialkit_xy_imprecise"` gives the spacing and
#' the extent (the condition carries them as `spacing` and `extent`).
#'
#' @param features_sf An sf object containing features to assign. A layer
#'   whose geometries mix dimensions or carry M (XYM or XYZM) is joined in x
#'   and y, and the rows returned keep their Z and M (see Details).
#' @param polygons_sf An sf or sfc polygonal layer. Polygons whose
#'   geometries mix dimensions or carry M (XYM or XYZM) are read in x and y.
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
#'   their number. Features that \code{tie_break} gives no polygon (see
#'   there) are kept or dropped with them.
#' @param predicate Binary spatial predicate function. Default sf::st_intersects.
#'   Not used when \code{largest} applies: sf then assigns polygon features by
#'   overlap area and never calls the predicate. With \code{sf::st_intersects}
#'   a POINT or MULTIPOINT feature with a point within a margin of a
#'   polygon is in it: one part in 1e10
#'   of the largest coordinate of the polygons, 0.44 mm at a UTM northing of
#'   4,400 km, at most a millionth of the diagonal of their bounding box
#'   and at least four times the spacing of doubles at that coordinate (in
#'   longitude and latitude, that many degrees of arc on the
#'   sphere while \code{sf::sf_use_s2()} is on, and no margin while it is
#'   off). A point on the edge two polygons share is then in both in any
#'   unit and origin of the coordinates, and \code{tie_break} decides it,
#'   where the last bits of its coordinates put it inside one polygon or the
#'   other.
#' @param largest Logical; when \code{features_sf} is itself polygonal, keep
#'   the polygon with the largest overlap. Default TRUE. Ignored for point and
#'   line features. A feature that only touches the polygon layer (shares an
#'   edge or a corner with it, with no overlap area) has no largest overlap
#'   and is unassigned, in any CRS; with \code{largest = FALSE} the default
#'   \code{st_intersects} counts touching, so such a feature is assigned.
#'   A feature that overlaps two or more polygons by exactly the same area
#'   (to 9 significant digits; a square split evenly across a cell edge) is
#'   counted in \code{"ties"} and given to the one \code{tie_break} names,
#'   so the choice does not depend on the order of the polygon rows; where
#'   \code{"smallest_area"} names none, the feature is given no polygon (see
#'   \code{tie_break}).
#'   Invalid geometries (usually a self-intersecting ring),
#'   whose overlap is undefined, are repaired with \code{sf::st_make_valid()}
#'   for the join, with a warning, and returned as they arrived. If the
#'   overlap still cannot be computed the function stops: falling back to
#'   \code{predicate} and \code{tie_break} would change the rule for every
#'   feature in the layer, so pass \code{largest = FALSE} to ask for that.
#'   In longitude and latitude with s2 on (sf's default), where s2 refuses
#'   a polygon of either layer as it arrived, the join by \code{predicate}
#'   meets that refusal too, and the error names the polygons and the
#'   transform instead (see Details).
#'   In longitude and latitude the validity check, the repair and the
#'   overlaps run on s2 whatever \code{sf::sf_use_s2()} says; with s2 off
#'   the join by \code{predicate} runs with GEOS, and the error that names
#'   the polygons also offers \code{largest = FALSE}.
#' @param tie_break Strategy for resolving features that match multiple
#'   polygons (with \code{largest}, that overlap several polygons equally):
#'   \code{"smallest_area"} (default) keeps the polygon with the
#'   smallest area and, among polygons of equal area (a point on the shared
#'   edge of two grid cells), the one whose bounding-box centre is lowest,
#'   then leftmost, then the one whose vertices' mean is lowest, then
#'   leftmost (the two triangles of a split square share a bounding box),
#'   so the choice does not depend on the order of the rows. Areas are
#'   compared to 9 significant digits, and centres and means within one part
#'   in 1e10 of the largest coordinate of the polygons, at most a
#'   millionth of the diagonal of their bounding box and at least four
#'   times the spacing of doubles at that coordinate (with any
#'   \code{predicate}, and in longitude and latitude whatever
#'   \code{sf::sf_use_s2()} says), so the choice does not depend on the unit
#'   or the origin of the coordinates either. Where polygons equal on every
#'   key are left (one polygon under two IDs, or a square and a rhombus of
#'   its area about its centre), the rule names no polygon: the feature is
#'   given none, \code{NA} in the ID column, counted in \code{"ties"} and
#'   listed in its \code{unresolved}, and a warning of class
#'   \code{"spatialkit_ties_unresolved"} gives how many such features there
#'   are (the condition carries \code{n} and \code{total}); it is kept or
#'   dropped as \code{keep_unassigned} says.
#'   Areas in longitude and latitude are measured on the sphere with s2,
#'   whatever \code{sf::sf_use_s2()} says. A polygon among those compared
#'   whose geometry is invalid (a ring whose edges cross, whose lobes GEOS
#'   measures as cancelling) is measured on a copy repaired with
#'   \code{sf::st_make_valid()}, with a warning, and returned as it arrived;
#'   in degrees, one s2 refuses even repaired (a bowtie) stops the call with
#'   the error of class \code{"spatialkit_cells_degenerate"} (see Details).
#'   \code{"first"} keeps the first match, the polygon that comes first
#'   in the rows of \code{polygons_sf} (original order-dependent
#'   behavior). On cells built from the centres a profile scored, a point on
#'   the edge between two of them is placed by the profile's rule whichever
#'   is given (see Value).
#' @param quiet Logical; suppress this function's \code{message()}s (the one
#'   naming a fallback ID column). It does not silence R warnings. Default
#'   \code{FALSE}.
#' @return An sf object with `polygon_id_col` attached, one row per input
#'   feature (fewer if `keep_unassigned = FALSE` dropped unmatched ones, and
#'   features `tie_break` gives no polygon), in
#'   the CRS `features_sf` arrived in. A column of `features_sf` already
#'   called `polygon_id_col` is dropped before the spatial join (with a
#'   warning), so re-assigning an already-assigned layer replaces the old IDs
#'   and does not fail. Every other column is kept, including one named like
#'   the polygons' own ID column when that is read from a fallback such as
#'   `"id"` (a site `id` joined to cells keyed by `id`). If *no* feature
#'   falls inside any polygon
#'   the result is empty (or all-`NA` with `keep_unassigned = TRUE`) and a
#'   warning is raised, since the usual cause is two layers in different
#'   places (a CRS that could only be stamped, not reprojected). A layer
#'   without a CRS is brought into the other layer's CRS as
#'   [harmonize_crs()] does; beside a layer in longitude and latitude, one
#'   with an x outside -180 to 360 or a y beyond 90 north or south cannot be
#'   in degrees, and the call stops with an error that names `features_sf`
#'   or `polygons_sf` and gives its box.
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
#'   the study area. A feature with empty geometry, or with a coordinate
#'   that is not finite or lies beyond 1e150, or in a geographic CRS a
#'   vertex at a latitude beyond 90 degrees (read as given, before any
#'   move into the polygons' CRS), has no location to assign: it matches no
#'   polygon, and the warning says how many of them there are. With
#'   `keep_unassigned = TRUE` they are kept with an
#'   `NA` ID and no warning is raised, unless every feature is unassigned.
#'   A polygon of `polygons_sf` with a coordinate that is not finite or lies
#'   beyond 1e150 (read as given) has no location either: it takes no
#'   feature, and the features inside the rest of its ring are unassigned.
#'
#'   The
#'   attribute `"ties"` records how many features matched more than one
#'   polygon (with `largest`, overlapped several by exactly the same area)
#'   and had the `tie_break` rule decide for them (on cells built from the
#'   centres a profile scored, the profile's rule then places the points
#'   among them, as the last paragraph says): a list with `n`,
#'   `which` (their row positions in `features_sf`), `rule`, `n_rows`
#'   (the number of rows returned, which the record was made for) and, when
#'   the rule gave some of them no polygon (see `tie_break`), `unresolved`,
#'   the row positions of those among `which`. A large
#'   `n` means the polygon layer overlaps, and per-cell counts built from the
#'   result depend on the rule. The record describes the rows this call
#'   returned and does not survive subsetting: `joined[i, ]`, like
#'   `dplyr::filter()`, `slice()` or `arrange()` of it, is a plain layer
#'   with no `"ties"` attribute, so nothing reports the parent's count
#'   for a different set of rows. `sf::st_drop_geometry()` keeps the record,
#'   since the rows are the same; see \code{\link{[.spatialkit_rows}} for
#'   what binding such data frames does. The `"unassigned"` record is
#'   stamped and removed the same way.
#'
#'   When `polygons_sf` carries `attr(, "count_rule")`, the record of how
#'   its cell count was set (cells from [build_tessellation()] on seeds
#'   placed from a [select_resolution()] result or from a profile built with
#'   an `estimate`), the result carries the same record.
#'   [summarize_by_cell()] reads it and checks that the rows and cells it is
#'   given are the ones the count was scored on. It is a plain attribute: a
#'   row subset keeps it and `dplyr::mutate()`, `dplyr::select()` and
#'   `dplyr::arrange()` drop it, which is why the summary also looks for it
#'   on its `cells_sf`. The profile scores a feature that is not a point at
#'   one point of it (`coerce_to_points(, "auto")`), in the cell of the
#'   nearest centre; this function gives a polygon feature to the cell it
#'   overlaps most. A feature that straddles a cell boundary can therefore
#'   be assigned to another cell than the one it was scored in, and the
#'   summary does not compare the two for such a layer. To assign the
#'   features as they were scored, assign their points:
#'   `assign_features_to_polygons(coerce_to_points(features_sf, "auto"), polygons_sf)`.
#'
#'   When that record says the cells are the ones the profile scored
#'   (`scored_cells` is `TRUE`: the cells [build_tessellation()] builds on a
#'   selection's seeds, or on the seeds a seeding function returns for one)
#'   and holds the centres it scored (`count_rule$scored$centres`), a point
#'   on the edge between two or more of those cells is given the cell
#'   [resolution_profile()] scored it in: that of its nearest centre,
#'   measured in the CRS the record names, and of centres the same distance
#'   away the lower-numbered one (a distance counts as the same as the least
#'   when its square is no more than the square of the least distance
#'   lengthened by a margin, plus one part in ten million of the least
#'   squared distance; the margin is one part in 1e10 of the largest
#'   coordinate of the centres, at most a millionth of the diagonal of their
#'   bounding box and at least four times the spacing of doubles at that
#'   coordinate). Such points occur on
#'   gridded locations (raster cell centres, a regular sampling grid), where
#'   a row can lie the same distance from two centres and the polygons hold
#'   it in either cell or in both. The rows each cell holds are then the rows
#'   its count was scored on, and [summarize_by_cell()]'s check of them reads
#'   `"match"`. A point counts as on the edge when it lies within a margin of
#'   both cells: one part in 1e10 of the largest coordinate of the cells,
#'   0.44 mm at a UTM northing of 4,400 km, at most a millionth of the
#'   diagonal of their bounding box and at least four times the spacing of
#'   doubles at that coordinate. For a point that two centres tie
#'   for while the cell holding it is that of the strictly nearer one, the
#'   margin also takes in the point's distance from the boundary between the
#'   two, measured in the CRS the record names: midway between centres far
#'   apart the tie reaches about 1.25e-8 of their distance from the boundary,
#'   which can exceed that margin.
#'   [build_tessellation()] draws such cells in the CRS the record names,
#'   around the centres as the record holds them, whatever CRS the seeds came
#'   back in (unless its `crs` names another projected CRS), so their corners
#'   carry only the rounding of that computation, far less than the margin.
#'   Cells it returns in a
#'   geographic `crs` are moved back into the record's CRS for the join (see
#'   Details). In a UTM zone the round trip through degrees moves a corner by
#'   about 2e-9 m; in the equal-area CRS the profile works in for a layer
#'   wider than a zone it can move one by more than the margin (how far
#'   depends on the PROJ build), so a point on such an edge can keep the
#'   polygon that holds it, and the summary's check then reads `"differ"`. A
#'   point farther from its nearest centre's cell keeps the polygon that
#'   holds it, and so do features that are not points, points without a CRS
#'   beside a record with one or with a CRS beside a record without one, and
#'   a join in longitude and latitude (polygons
#'   in degrees beside a record in degrees). A centre's cell
#'   is the one cell that holds it and no other centre; a centre outside a
#'   boundary that clipped its cell away from it has the cell, among those
#'   that hold no centre, whose inside is nearest it. Cells built from other
#'   seeds leave most centres without a cell, and their points to the
#'   polygons. Polygons without the record, such as polygons
#'   the caller brings, and cells whose record says `scored_cells` is
#'   `FALSE` or `NA` (other cells built from the count, such as a lattice
#'   sized by it) are placed by the polygons alone.
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
  # The record of how the cells' count was set, when they carry one: copied
  # onto the assigned rows, where summarize_by_cell() looks for it first.
  count_rule <- .count_rule_get(polygons_sf)

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
  # A layer that mixes XY and XYZ (or XYM) geometries is joined in XY.  sf
  # reads a layer in the dimension of its first geometry: s2 drops Z and M
  # only when the first geometry carries them, and GEOS and the transform
  # read an XYZ ring among XY rows as XY pairs.  With a clean XYZ square in
  # row 3 of six XY squares, the point in it fell in no polygon and three
  # others in two; in a projected CRS GEOS stopped on the ring.  So is a
  # layer whose geometries carry M (XYM or XYZM): GEOS does not read M and
  # stopped the join, the overlap and the validity check with "GEOS does not
  # support XYM or XYZM geometries".  The join is defined in x and y and
  # reads M nowhere.  Every geometric step of the join below reads these
  # copies; a layer XY or XYZ throughout is read as it came.
  features_dims <- .geom_dims(features_sf)
  features_flat <- length(features_dims) > 1L || any(features_dims %in% c("XYM", "XYZM"))
  features_xy <- if (features_flat) sf::st_zm(features_sf, drop = TRUE, what = "ZM")
                 else features_sf
  # Points in an sfc_GEOMETRY column are joined as the same points in
  # sfc_POINT (.typed_points()): on the column sf::st_cast() returns, the
  # join stopped in sf::st_dimension() with "Not compatible with STRSXP".
  features_xy <- .typed_points(features_xy)
  # A feature with a coordinate that is not finite or lies beyond 1e150 has
  # no location.  It is read as given, before any move, and the join's copy
  # of it is emptied below, so it matches no polygon and is counted with the
  # empty features: moved into the polygons' CRS, a northing of 1e200 in
  # EPSG:32616 came back finite in EPSG:32617 and was joined there.
  # A feature of a layer in a geographic CRS with a vertex at a latitude
  # beyond 90 degrees has no location either (lat = TRUE).  Its copy is
  # emptied here, before the move and the join read it: it fell in no
  # polygon, but sf's bare "bounding box has potentially an invalid value
  # range for longlat data" came with it, and the texts below called it a
  # feature the polygons miss.
  features_ll <- .is_longlat(features_xy)
  unlocated_f <- .unlocated_features(sf::st_geometry(features_xy), lat = TRUE)
  features_xy <- .empty_unlocated(features_xy, .beyond_pole_rows(features_xy), lat = TRUE)
  # A polygon with a vertex that is not finite or lies beyond 1e150, read as
  # given, has no location either, and its copy is emptied before the join:
  # read as it was, one with an NA or Inf vertex caught the points inside the
  # rest of its ring (and some of its neighbours'), and one at 1e200 stopped
  # GEOS ("orientationIndex encountered NaN/Inf numbers").  It takes no
  # point, and the points in it are counted with the other unassigned ones.
  polygons_xy <- .empty_unlocated(.xy_if_mixed_or_m(polygons_sf))
  # A row with a location that the move cannot place in full comes back
  # empty, or without the vertices it cannot place, and is said in this
  # function's words, naming `features_sf` or `polygons_sf`, in place of
  # harmonize_crs()'s, which names `a` or `b`.
  # A layer without a CRS beside one in longitude and latitude is not given
  # that CRS when its coordinates cannot be in it (x outside -180 to 360, y
  # beyond 90 north or south), and the refusal names `features_sf` or
  # `polygons_sf`, not harmonize_crs()'s `a` or `b`.  Stamped, features in
  # metres fell in no polygon, and polygons in metres stopped on an s2
  # refusal that blamed the lon/lat reading the stamp had made.
  if (is.na(sf::st_crs(features_xy)) && .is_longlat(polygons_xy))
    .refuse_lonlat_stamp(.looks_like_lonlat(features_xy)$bb, sf::st_crs(polygons_xy),
                         caller = "assign_features_to_polygons", what = "features_sf")
  if (is.na(sf::st_crs(polygons_xy)) && .is_longlat(features_xy))
    .refuse_lonlat_stamp(.looks_like_lonlat(polygons_xy)$bb, sf::st_crs(features_xy),
                         caller = "assign_features_to_polygons", what = "polygons_sf")
  done <- .answer_for_classed("spatialkit_rows_unplaced")
  hh <- tryCatch(harmonize_crs(features_xy, polygons_xy,
                               prefer = if (join_in_p) "b" else "a"),
                 finally = done())
  f <- hh$a; p <- hh$b
  # A layer without a CRS that harmonize_crs() read as lon/lat and moved is
  # compared as EPSG:4326, so the warning names longitude and latitude.
  .warn_rows_unplaced(.as_harmonize_read(features_xy, f), f,
                      caller = "assign_features_to_polygons",
                      arg = "features_sf", remedy = .UNPLACED_REMEDY)
  # `f` is only the join's copy.  The rows returned keep the geometry the
  # features arrived with, Z and M included; a CRS-less layer is resolved
  # into the polygons' CRS and returned there, as before.  A CRS-less layer
  # read in XY for the join keeps its Z and M when it was only given a CRS.
  # When it was moved from lon/lat, a layer in one dimension (M throughout)
  # comes back moved with its M, since sf transforms such a layer right; one
  # that mixes dimensions comes back moved, in XY, and a warning says so,
  # unless sf transforms the layer right (below).
  f_geom <- if (is.na(orig_crs)) sf::st_geometry(f) else sf::st_geometry(features_sf)
  if (is.na(orig_crs) && features_flat) {
    if (identical(f_geom, sf::st_set_crs(sf::st_geometry(features_xy), sf::st_crs(f))))
      f_geom <- sf::st_set_crs(sf::st_geometry(features_sf), sf::st_crs(f))
    else if (identical(.xy_if_mixed_warn(features_sf, "assign_features_to_polygons",
                                         "features_sf"), features_sf))
      # A layer in one dimension (M throughout), or one whose only geometries
      # of another dimension than the first are empty and not points, which
      # sf transforms right: the rows are moved with their Z and M, from
      # EPSG:4326 as harmonize_crs() took them.
      f_geom <- .transform_keep_unlocated(
        sf::st_set_crs(sf::st_geometry(features_sf), 4326), sf::st_crs(f))
  }
  # The join's copy of a feature without a location is an empty geometry of
  # its own type and dimension; the row returned keeps `f_geom`.
  f <- .empty_unlocated(f, which(unlocated_f))
  # Cells in longitude and latitude whose record says a profile scored them
  # (the Voronoi cells of the scored centres that build_tessellation()
  # returns in a geographic `crs`) were drawn in the record's projected CRS,
  # and are joined there: their edges are straight in it, and the profile
  # measured the rows in it.  The join's copy of the features is moved there
  # from whatever CRS it is in; the rows returned keep `f_geom`, set above.
  # On the sphere s2 refused some of these cells: where four centres stand
  # the same distance from one point, two cells meet at two corners a
  # fraction of a nanometre apart, which the transform to degrees puts at
  # one point ("duplicate vertex"), and the join of rows in degrees stopped.
  # The polygons are moved from `polygons_xy`, not from `p`, which holds
  # them in the features' CRS when that is projected.  Projected cells are
  # joined in their own CRS whatever CRS their record names, and cells
  # without a CRS have none to be moved from.  A polygon with a vertex that
  # is not finite or lies beyond 1e150, read as given, comes back empty
  # (.transform_keep_unlocated()), as harmonize_crs() returns it: moved as
  # it was, it came back as a smaller shape that took points.
  # The features are moved with .transform_keep_unlocated() as well, and a
  # row either move cannot place in full is said, as above.  The polygons
  # are said for the move whose copy is joined.
  crs_s <- if (!join_in_p && !is.na(crs_p)) .scored_crs(count_rule) else NULL
  if (is.null(crs_s)) {
    .warn_rows_unplaced(.as_harmonize_read(polygons_xy, p), p,
                        caller = "assign_features_to_polygons",
                        arg = "polygons_sf", remedy = .UNPLACED_REMEDY)
  } else {
    if (!isTRUE(sf::st_crs(f) == crs_s)) {
      f0 <- f
      f  <- .transform_keep_unlocated(f0, crs_s)
      .warn_rows_unplaced(f0, f, caller = "assign_features_to_polygons",
                          arg = "features_sf", remedy = .UNPLACED_REMEDY)
    }
    p <- .transform_keep_unlocated(polygons_xy, crs_s)
    .warn_rows_unplaced(polygons_xy, p, caller = "assign_features_to_polygons",
                        arg = "polygons_sf", remedy = .UNPLACED_REMEDY)
  }

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

  # In longitude and latitude the overlaps and the areas the join compares
  # are measured on the sphere with s2, whatever sf::sf_use_s2() says, as the
  # package's other areas in degrees are (.with_s2()).  With s2 off sf hands
  # st_area() to lwgeom, which the package does not depend on: the
  # tie-break and the largest-overlap join stopped with "package lwgeom
  # required, please install it first", and the latter advised `largest =
  # FALSE`, which stopped the same way.  The join by `predicate` keeps the
  # session's engine.
  on_sphere <- isTRUE(sf::st_is_longlat(sf::st_crs(p_sel)))
  s2_if_sphere <- function(expr) if (on_sphere) .with_s2(expr) else expr
  # The session's engine, read before any step is put on s2: with s2 off the
  # join by `predicate` runs with GEOS, which reads in degrees polygons s2
  # refuses, so the largest-overlap join's error still offers it.
  s2_session <- isTRUE(sf::sf_use_s2())

  # The layers as they arrived at the join, before the largest-overlap
  # join's repair: s2's refusal of a piece of an overlap is traced to the
  # polygons it refused as they arrived (.stop_join_degenerate()).
  f_in <- f
  p_in <- p_sel
  # The edges of the polygons, and the points read against them, are held
  # to the spacing of doubles at the coordinates' size; where that is more
  # than a millionth of the layers' extent it is said (.warn_xy_imprecise()).
  .warn_xy_imprecise(list(f, p_sel), "assign_features_to_polygons",
                     c("features_sf", "polygons_sf"))
  # The margin within which a point lies on a polygon's edge, and within
  # which the tie-break's bounding-box centres and vertex means count as
  # equal: one part in 1e10 of the largest coordinate of the polygons, at
  # most a millionth of their diagonal (.edge_tol()), the margin of the
  # edge rule of cells a profile scored.
  bb_p    <- sf::st_bbox(p_sel)
  tie_tol <- .edge_tol(bb_p)
  if (use_largest) s2_if_sphere({
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
  })

  join_args <- list(x = f, y = p_sel, join = predicate, left = TRUE)
  largest_ties <- integer(0)
  # Features the smallest-area rule leaves between polygons it cannot tell
  # apart, by their row in `features_sf`: given no polygon.
  unresolved_rows <- integer(0)
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
    # In degrees the join runs on s2 whatever the session says (above), and
    # so does its handler, which asks s2 which polygons it refuses.
    joined <- s2_if_sphere(tryCatch(
      withCallingHandlers(
        do.call(sf::st_join, c(join_args, largest = TRUE)),
        warning = function(w) {
          if (grepl("attribute variables are assumed to be spatially constant",
                    conditionMessage(w), fixed = TRUE))
            invokeRestart("muffleWarning")
        }),
      error = function(e) {
        # The polygons s2 refuses were repaired above, and s2 accepts most
        # of them then.  It still refuses some (a hole that crosses its
        # outer ring stays two rings whose edges cross), and `largest =
        # FALSE` would meet the same refusal, so the error names the cause
        # and the transform that helps instead of advising it.  So it does
        # where the repair made s2 accept a polygon and s2 then refused a
        # piece of its overlap with the other layer (a figure eight,
        # repaired into two parts that touch at a point): the join by
        # `predicate` meets that polygon as it arrived, where the session
        # has s2 on; with s2 off it runs with GEOS, and the error offers it.
        .stop_join_degenerate(e, f_in, p_in, f, p_sel, is.na(orig_crs), is.na(crs_p),
                              predicate_on_s2 = s2_session)
        # `largest = FALSE` is advised only where that join does not meet
        # the same error: an engine sf lacks (lwgeom) stops it too.
        in_degrees <- isTRUE(sf::st_is_longlat(f))
        advice <- if (!grepl("lwgeom", conditionMessage(e), fixed = TRUE))
          paste0(" Pass largest = FALSE to assign every feature by `predicate` and ",
                 "`tie_break` instead",
                 if (in_degrees) ", or project both layers first (see ensure_projected())",
                 ".")
        else if (in_degrees) " Project both layers first (see ensure_projected())."
        else ""
        stop(sprintf("assign_features_to_polygons(): the largest-overlap join failed (%s).%s",
                     conditionMessage(e), advice),
             call. = FALSE)
      }))
    # sf keeps the largest intersection piece without asking whether it has
    # any area, so under GEOS a feature that only shares an edge or a corner
    # with the cells was assigned to one of them with zero overlap, while
    # under s2 (lon/lat, s2 on) the same feature had no piece and came back
    # unassigned: sf's 100 nc counties against 50 of them as cells gave 70
    # rows projected and 50 in lon/lat.  A feature with no overlap has no largest
    # overlap, so it is unassigned whatever the CRS.  s2 already does this,
    # and joins every layer in degrees here.
    if (!on_sphere) {
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
    lt <- s2_if_sphere(.largest_overlap_ties(f, p_sel, joined, join_id, tie_break,
                                             tie_tol))
    if (length(lt$rows)) {
      joined[[join_id]][lt$rows] <- lt$ids
      largest_ties <- joined[["..pre_join_row_id"]][lt$rows]
      unresolved_rows <- largest_ties[is.na(lt$ids)]
      .log_info(paste0("assign_features_to_polygons(): %d of %d polygon ",
                       "feature(s) overlap two or more polygons by exactly the ",
                       "same area; the '%s' rule chose one for each."),
                length(lt$rows), nrow(f), tie_break)
    }
  } else {
    # In longitude and latitude sf joins with s2, which refuses a polygon
    # with two corners it reads as one point.  Cells returned in degrees hold
    # such corners on gridded locations, and the join stopped with s2's
    # "Loop 0 is not valid: Edge 20 is degenerate (duplicate vertex)", which
    # names neither the polygon nor a way on: the k-means cells of a 30 by 20
    # lattice at 30 m returned in EPSG:4326 stopped it once they had lost
    # their record or their CRS.  Such a refusal, and every other refusal
    # of a polygon of either layer by s2, now stops with the cause and the
    # remedy; every other error passes as it was raised.
    # A point on the edge two polygons share lies in both, and the default
    # predicate found both only while the edge passed through it exactly: a
    # change of unit or a shift of the origin moves the last bits of the
    # coordinates, which put a lattice point a hair inside one cell or the
    # other (25 of 925 points of a 37 x 25 lattice changed cell in US feet,
    # 11 with the false origin moved by 0.1 m).  With sf::st_intersects() a
    # point within `tie_tol` of a polygon is in it, so an edge point is a tie
    # in every frame and `tie_break` decides it (.intersects_margin()).
    if (identical(predicate, sf::st_intersects))
      join_args$join <- function(x, y, ...) .intersects_margin(x, y, bb_p)
    joined <- withCallingHandlers(
      do.call(sf::st_join, join_args),
      error = function(e)
        .stop_join_degenerate(e, f, p_sel, f, p_sel, is.na(orig_crs), is.na(crs_p)))
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
    # last bits count as equal.  Only the polygons the tied features fall in
    # are measured, those whose areas the rule compares.
    tie_ids <- unique(as.character(
      joined[[polygon_id_col]][joined[["..pre_join_row_id"]] %in% tie_rows]))
    tie_ids <- tie_ids[!is.na(tie_ids)]
    tie_at  <- match(tie_ids, as.character(p[[id_col]]))
    p_tie   <- sf::st_geometry(p)[tie_at]
    poly_areas <- s2_if_sphere(.tie_break_areas(p_tie, tie_at, nrow(p), is.na(crs_p),
                                                s2_session))
    poly_areas[!is.finite(poly_areas)] <- Inf
    names(poly_areas) <- tie_ids

    joined$`..poly_area` <- signif(poly_areas[as.character(joined[[polygon_id_col]])], 9L)
    joined$`..poly_area`[is.na(joined$`..poly_area`)] <- Inf

    # Equal areas -- every cell of a regular grid, for a point on a shared
    # edge -- used to fall through to row order, so the answer depended on
    # the order of the polygon layer after all: edge points at x = 100 went
    # to cells 1, 4, 7, or to 2, 5, 8 with the rows reversed.  The candidate
    # whose bounding-box centre is lowest, then leftmost, wins instead.  On
    # a create_grid_polygons() grid, which is numbered from the lower left a
    # row at a time, that is the cell the row order picked.  The centres are
    # compared within `tie_tol`: in another unit or origin the centres of a
    # row of cells differ in their last bits, which chose the cell of an edge
    # point.  Candidates that share a bounding box (the two triangles of a
    # split square) are told apart by the mean of their vertices, lowest,
    # then leftmost; candidates equal on every key (one polygon under two
    # IDs, a square and a rhombus that share its centre) went to the first
    # in the row order of `polygons_sf`, and are now given no polygon
    # (.tie_rule_rows()).
    ctr     <- vapply(p_tie, function(g) {
      bb <- sf::st_bbox(g)
      c((bb[["xmin"]] + bb[["xmax"]]) / 2, (bb[["ymin"]] + bb[["ymax"]]) / 2)
    }, numeric(2))
    vm    <- .vertex_means(p_tie)
    at_j  <- match(as.character(joined[[polygon_id_col]]), tie_ids)
    key   <- cbind(joined[["..poly_area"]], ctr[2, at_j], ctr[1, at_j],
                   vm[at_j, 2L], vm[at_j, 1L])
    grp   <- joined[["..pre_join_row_id"]]
    keep  <- .tie_rule_rows(key, grp, tie_tol)
    joined <- joined[keep, , drop = FALSE]
    many  <- duplicated(joined[["..pre_join_row_id"]])
    open  <- unique(joined[["..pre_join_row_id"]][many])
    joined <- joined[!many, , drop = FALSE]
    joined[[polygon_id_col]][joined[["..pre_join_row_id"]] %in% open] <- NA
    unresolved_rows <- c(unresolved_rows, open)
    joined[["..poly_area"]] <- NULL
  } else {
    joined <- joined[!duplicated(joined[["..pre_join_row_id"]]), , drop = FALSE]
  }
  # On cells built from the centres a profile scored, a point on the edge
  # between two of them goes to the cell the profile scored it in, whatever
  # `tie_break` picked or which side of the edge rounding put it on: its
  # nearest centre, a tie to the lower-numbered one.  Then the rows each cell
  # holds are the rows its count was scored on.
  if (!is.null(count_rule)) {
    on_edge <- .scored_edge_cells(f, p, id_col, joined, polygon_id_col, count_rule,
                                  features_xy)
    if (length(on_edge$rows)) {
      joined[[polygon_id_col]][on_edge$rows] <- on_edge$ids
      .log_info(paste0("assign_features_to_polygons(): %d point(s) on the edge between two ",
                       "cells resolution_profile() scored were given the cell it scored ",
                       "them in, that of their nearest centre."),
                length(on_edge$rows))
    }
  }
  # Which features matched no polygon, by their row in `features_sf`, read
  # before the row key goes.  A feature the tie-break left between polygons
  # it cannot tell apart matched some, and is reported with the ties.
  unresolved_rows <- sort(unique(unresolved_rows))
  unassigned_rows <- sort(setdiff(
    joined[["..pre_join_row_id"]][is.na(joined[[polygon_id_col]])], unresolved_rows))
  joined[["..pre_join_row_id"]] <- NULL
  if (length(unresolved_rows))
    .warn_classed("spatialkit_ties_unresolved",
                  paste0("assign_features_to_polygons(): %d of %d feature(s) match two or ",
                         "more polygons that the smallest-area rule of `tie_break` cannot ",
                         "tell apart (equal in area, in bounding-box centre and in the mean ",
                         "of their vertices, such as one polygon under two IDs), so they ",
                         "are given no polygon and %s; attr(, \"ties\")$unresolved gives ",
                         "their rows in `features_sf`."),
                  length(unresolved_rows), nrow(joined),
                  if (keep_unassigned) "kept with an NA ID" else "left out of the result",
                  data = list(n = length(unresolved_rows), total = nrow(joined)))

  n_unassigned <- length(unassigned_rows)
  n_features   <- nrow(joined)
  # A feature without a location (an empty geometry, or a coordinate that is
  # not finite or lies beyond 1e150, read as given) is in no polygon too, but
  # for want of a location, not because the polygons miss it: the texts say
  # so, as build_tessellation() does, and keep their advice for the rest.
  no_loc  <- sf::st_is_empty(f_geom) | unlocated_f
  n_empty <- if (n_unassigned) sum(no_loc[unassigned_rows]) else 0L
  # Zero matches is almost never what the caller meant, and it is silent: the
  # function returns an empty sf and every downstream step reports zero rows
  # without saying why.  The usual cause is a CRS mismatch that could only be
  # stamped rather than reprojected, which puts the two layers in different
  # places.
  if (n_unassigned == n_features && n_features > 0L)
    .warn_classed("spatialkit_unassigned",
                  paste0("assign_features_to_polygons(): none of the %d ",
                         "feature(s) fall inside any polygon. %s"),
                  n_features,
                  if (n_empty == n_features)
                    sprintf("They have no location (%s), so there is nothing to assign.",
                            .unlocated_what(features_ll))
                  else paste0(
                    if (n_empty > 0L)
                      sprintf("%d of them %s no location (%s); for the rest, check", n_empty,
                              if (n_empty == 1L) "has" else "have", .unlocated_what(features_ll))
                    else "Check",
                    " that the two layers cover the same ground -- a CRS that had to ",
                    "be stamped rather than reprojected is the usual ",
                    "cause -- or pass a different `predicate`."),
                  data = list(n = n_unassigned, total = n_features))
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
                         "rows in `features_sf`. %s"),
                  n_unassigned, n_features, 100 * n_unassigned / n_features,
                  if (n_empty == n_unassigned)
                    sprintf(paste0("They have no location (%s), so there is nothing ",
                                   "to assign; pass keep_unassigned = TRUE to keep them ",
                                   "with an NA ID."), .unlocated_what(features_ll))
                  else paste0(
                    if (n_empty > 0L)
                      sprintf("%d of them %s no location (%s); for the rest, check", n_empty,
                              if (n_empty == 1L) "has" else "have", .unlocated_what(features_ll))
                    else "Check",
                    " that the polygons cover ",
                    "the features (cells built from seeds need a ",
                    "`boundary` or a `cover`; see build_tessellation()), or ",
                    "pass keep_unassigned = TRUE to keep them with an NA ID."),
                  data = list(n = n_unassigned, total = n_features))
  else if (n_unassigned > 0L)
    .log_info("assign_features_to_polygons(): %d of %d feature(s) fall inside no polygon and are kept with an NA ID%s.",
              n_unassigned, n_features,
              if (n_empty > 0L) sprintf(" (%d of them %s no location: %s)", n_empty,
                                        if (n_empty == 1L) "has" else "have",
                                        .unlocated_what(features_ll)) else "")

  if (!keep_unassigned) {
    joined <- joined[!is.na(joined[[polygon_id_col]]), , drop = FALSE]
  }

  # Stamped and classed: the record names row positions, so it must not
  # survive a subset that renumbers or removes them.
  tie_rows <- sort(unique(c(tie_rows, largest_ties)))
  joined <- .set_row_record(joined, "ties",
                            c(list(n = length(tie_rows),
                                   which = as.integer(tie_rows),
                                   rule = tie_break),
                              # Only when the rule gave some feature no polygon.
                              if (length(unresolved_rows))
                                list(unresolved = as.integer(unresolved_rows))))
  # The features no polygon holds, beside the ties: the count used to exist
  # only in a log line.  Stamped the same way, and removed by a subset with
  # the other records named in .row_record_attrs (R/utils.R).
  joined <- .set_row_record(joined, "unassigned",
                            list(n = n_unassigned,
                                 which = as.integer(unassigned_rows)))
  .count_rule_put(joined, count_rule)
}


# What the caller of a join does once the polygons are transformed.
.JOIN_AGAIN <- "assign again; the join then runs in that CRS"

#' Stop with the cause when s2 refuses a polygon of a layer in degrees
#'
#' sf measures geometry in longitude and latitude with s2, which refuses a
#' polygon that is not valid on the sphere: a ring with two corners it reads
#' as one point ("Loop 0 is not valid: Edge 20 is degenerate (duplicate
#' vertex)") or with edges that cross; a hole that crosses its outer ring
#' ("Loop 0 edge 1 crosses loop 1 edge 0"); two rings that share an edge
#' ("Loop 1 edge 2 has duplicate near loop 2 edge 3").  Its error names
#' neither the polygon nor a way on.  Cells built in a projected CRS and
#' returned in degrees hold such corners where four seeds stand the same
#' distance from one point: two cells meet there at two corners that differ
#' only by rounding, which the transform puts at one point.  Called from an
#' error handler around an sf call on \code{polygons}.
#'
#' @param e The error caught.  s2's refusals of a polygon all name a ring
#'   ("Loop N ...", "... crosses loop N ...", "Inconsistent loop orientations
#'   detected" for oriented polygons), so an error that names none is not
#'   one; sf's own check of the polygons below then tells s2's refusal of
#'   these polygons from that of another layer's.
#' @param polygons The polygons as the failed call had them (sf or sfc).
#' @param caller Name of the exported function, for the message.
#' @param layer The argument the polygons came in, as the message names it.
#' @param task What cannot run, as in "so \code{task} cannot run in degrees".
#' @param again What to do once the polygons are transformed.
#' @param no_crs Whether the polygons came without a CRS and were read as
#'   longitude and latitude.
#' @param pick The projected CRS the remedy names after the one the cells
#'   were built in, as in "or \code{pick}".
#' @param ... Passed to \code{.raise_cells_degenerate()} (\code{of},
#'   \code{source}, \code{remedy}, \code{site}, \code{s2_forced}).
#' @return \code{NULL}, invisibly, when \code{e} is not s2's refusal of a
#'   polygon or s2 refuses no polygon of \code{polygons} (the refusal is
#'   another layer's, or of a piece the call made): the caller then raises
#'   \code{e} as it would.  Otherwise it stops with the error
#'   \code{.raise_cells_degenerate()} builds.
#' @keywords internal
#' @noRd
.stop_cells_degenerate <- function(e, polygons, caller, layer, task, again,
                                   no_crs = FALSE,
                                   pick = "the one ensure_projected() picks for it",
                                   ...) {
  if (!.is_s2_ring_error(e)) return(invisible(NULL))
  g <- sf::st_geometry(polygons)
  if (!isTRUE(sf::st_is_longlat(g)) || !isTRUE(sf::sf_use_s2())) return(invisible(NULL))
  bad <- .s2_refused_rows(g)
  if (is.null(bad) || !length(bad$rows)) return(invisible(NULL))
  .raise_cells_degenerate(e, bad$rows, bad$reasons, length(g), caller, layer, task,
                          again, no_crs, pick = pick, ...)
}

# The ID columns summarize_by_cell() looks for, in order: among the rows'
# columns (`.sbc_id_candidates()`), and, with the one found there as
# `id_col`, among the cells' columns (`.cells_id_candidates()`), in the order
# assign_features_to_polygons() reads the polygons' IDs from.
.sbc_id_candidates <- function(id_col)
  unique(c(id_col, "poly_id", "polygon_id", "cell_id"))
.cells_id_candidates <- function(id_col)
  unique(c(id_col, "poly_id", "polygon_id", "id", "cell_id", "grid_id"))

#' Refuse projected cells GEOS calls invalid, before a density is taken
#'
#' \code{summarize_by_cell(area = TRUE)} divides by each cell's planar area,
#' which GEOS takes on the cell's rings as they are: a ring whose edges cross
#' gets the signed sum of its lobes, and parts that overlap are counted
#' twice.  Every non-empty cell is checked with \code{sf::st_is_valid()} in
#' its CRS (\code{NA} counts as invalid), and the request is refused when
#' one is invalid, naming the first six such cells by the ID column the
#' summaries are joined on (by row when there is none), GEOS's reason for
#' the first, and the remedy, \code{sf::st_make_valid()}.
#'
#' @param cells The cells as they are measured (in XY when they mix
#'   dimensions or carry M), row for row with \code{cells_sf}.
#' @param cells_sf The cells as passed, for their ID column.
#' @param id_col,assigned_points_sf As passed to \code{summarize_by_cell()}:
#'   the ID column is found as the join finds it.
#' @return \code{NULL}, invisibly, when every non-empty cell is valid;
#'   otherwise it stops.
#' @keywords internal
#' @noRd
.stop_cells_invalid <- function(cells, cells_sf, id_col, assigned_points_sf) {
  g    <- sf::st_geometry(cells)
  keep <- which(!sf::st_is_empty(g))
  if (!length(keep)) return(invisible(NULL))
  bad <- keep[!(suppressWarnings(sf::st_is_valid(g[keep])) %in% TRUE)]
  if (!length(bad)) return(invisible(NULL))
  why <- tryCatch(suppressWarnings(sf::st_is_valid(g[bad[1L]], reason = TRUE)),
                  error = function(e) NA_character_)
  # The ID column of the join, found as the join finds it; by row when the
  # cells have none (the join then refuses the request with its own text).
  idc <- tryCatch({
    pid <- .sbc_id_candidates(id_col)
    pid <- pid[pid %in% names(assigned_points_sf)]
    cid <- .cells_id_candidates(if (length(pid)) pid[[1L]] else id_col)
    cid <- cid[cid %in% names(cells_sf)]
    if (length(cid)) cid[[1L]] else NULL
  }, error = function(e) NULL)
  label <- if (is.null(idc)) as.character(bad)
           else as.character(sf::st_drop_geometry(cells_sf)[[idc]][bad])
  k     <- length(bad)
  what  <- if (is.null(idc)) "row" else idc
  whats <- if (is.null(idc) && k > 1L) "rows" else what
  named <- if (k == 1L) sprintf("%s %s", whats, label)
           else if (k <= 6L) sprintf("%s %s and %s", whats,
                                     paste(label[-k], collapse = ", "), label[k])
           else sprintf("%s %s and %d more", whats, paste(label[1:6], collapse = ", "), k - 6L)
  on   <- if (k == 1L) "" else sprintf(" for %s %s", what, label[1L])
  geos <- if (length(why) == 1L && !is.na(why) && nzchar(why))
    sprintf("GEOS's reason%s: \"%s\"", on, why)
  else sprintf("GEOS gives no reason%s", on)
  of <- if (length(g) == 1L) "the one cell of `cells_sf`"
        else sprintf("%d of the %d cells of `cells_sf`", k, length(g))
  stop(sprintf(paste0(
    "summarize_by_cell(): `area = TRUE` refused: %s %s in its CRS (%s; %s), so %s ",
    "planar area, and a density over it, cannot be vouched for. Repair the cells ",
    "with sf::st_make_valid() and pass them again."),
    of, if (k == 1L) "is not a valid polygon" else "are not valid polygons",
    named, geos, if (k == 1L) "its" else "each one's"), call. = FALSE)
}

# s2's refusals of a polygon all name a ring ("Loop N ...", "... crosses
# loop N ...", "Inconsistent loop orientations detected").
.is_s2_ring_error <- function(e)
  grepl("[Ll]oop [0-9]+|loop orientations", conditionMessage(e))

# The rows of `x` (sf or sfc) that sf's validity check refuses, with its
# reason for each, as a list of `rows` and `reasons`; NULL when the check
# cannot say.  With `polygonal`, only polygons count: a ring error is not a
# point's or a line's.  Run it under the engine whose refusal is asked about.
.s2_refused_rows <- function(x, polygonal = FALSE) {
  g <- sf::st_geometry(x)
  check <- function(...) tryCatch(suppressWarnings(sf::st_is_valid(g, ...)),
                                  error = function(err) NULL)
  ok  <- check()
  why <- check(reason = TRUE)
  if (length(ok) != length(g) || length(why) != length(g)) return(NULL)
  bad <- ok %in% FALSE
  if (polygonal) bad <- bad & (sf::st_dimension(g) %in% 2L)
  list(rows = which(bad), reasons = unname(why[bad]))
}

#' Stop with the cause when s2 refuses a polygon of either layer of the join
#'
#' The handler of both joins of \code{assign_features_to_polygons()}.  s2
#' can refuse a polygon of \code{polygons_sf} or of \code{features_sf}, and
#' on the largest-overlap join, which repairs the polygons it refuses first,
#' a piece of an overlap of two polygons it accepts: a figure eight, one
#' vertex listed twice, is repaired into two parts that touch at a point, and
#' a piece of its overlap with the other layer can hold that point twice.
#' The layers the join used are examined first, the polygons before the
#' features, then the layers as they arrived, in the same order; the first
#' layer s2 refuses a polygon of is named.  So a polygon s2 refuses even
#' repaired is named before one the repair made s2 accept.  When s2 refuses
#' no polygon of either layer as it arrived, the join by \code{predicate}
#' runs, and the caller's message, which advises it, stands.
#'
#' @param e The error caught.
#' @param f_in,p_in The features and polygons as they arrived at the join.
#' @param f_used,p_used The features and polygons the join used (repaired on
#'   the largest-overlap join; the same as \code{f_in} and \code{p_in} on the
#'   join by \code{predicate}).
#' @param no_crs_f,no_crs_p Whether the features or the polygons came without
#'   a CRS and were read as longitude and latitude.
#' @param predicate_on_s2 Whether the join by \code{predicate} runs on s2:
#'   the session's \code{sf::sf_use_s2()}, read before the largest-overlap
#'   join was put on s2.  When it is \code{FALSE}, that join runs with GEOS,
#'   so the error names the largest-overlap join and offers
#'   \code{largest = FALSE}.
#' @return \code{NULL}, invisibly, when \code{e} is not s2's refusal of a
#'   ring or s2 refuses no polygon of either layer; otherwise it stops with
#'   the error \code{.raise_cells_degenerate()} builds, of class
#'   \code{"spatialkit_cells_degenerate"} for the polygons and
#'   \code{"spatialkit_features_degenerate"} for the features.
#' @keywords internal
#' @noRd
.stop_join_degenerate <- function(e, f_in, p_in, f_used, p_used, no_crs_f = FALSE,
                                  no_crs_p = FALSE, predicate_on_s2 = TRUE) {
  if (!.is_s2_ring_error(e)) return(invisible(NULL))
  if (!isTRUE(sf::st_is_longlat(sf::st_geometry(p_used))) || !isTRUE(sf::sf_use_s2()))
    return(invisible(NULL))
  same <- identical(f_in, f_used) && identical(p_in, p_used)
  pa <- .s2_refused_rows(p_used)
  pb <- if (same) pa else .s2_refused_rows(p_in)
  fa <- .s2_refused_rows(f_used, polygonal = TRUE)
  fb <- if (same) fa else .s2_refused_rows(f_in, polygonal = TRUE)
  if (is.null(pa) || is.null(pb) || is.null(fa) || is.null(fb)) return(invisible(NULL))
  n_p <- nrow(p_used); n_f <- nrow(f_used)
  task <- if (predicate_on_s2) "the join" else "the largest-overlap join"
  forced <- !isTRUE(predicate_on_s2)
  cells <- function(bad, repaired)
    .raise_cells_degenerate(e, bad$rows, bad$reasons, n_p, "assign_features_to_polygons",
                            "`polygons_sf`", task, .JOIN_AGAIN, no_crs_p,
                            repaired = repaired,
                            also = list(rows = fb$rows, reasons = fb$reasons, n = n_f),
                            s2_forced = forced, advise_predicate = forced)
  features <- function(bad, repaired) {
    # The join by `predicate` repairs nothing: the caller's repair helps
    # where it makes s2 accept every row named, which depends on the ring
    # (s2 accepts a repeated vertex and a figure eight repaired, not a
    # bowtie or a hole that crosses its outer ring).
    fixable <- !repaired && same && isTRUE(tryCatch(
      all(sf::st_is_valid(.safe_make_valid(sf::st_geometry(f_in)[bad$rows])) %in% TRUE),
      error = function(err) FALSE))
    .raise_cells_degenerate(e, bad$rows, bad$reasons, n_f, "assign_features_to_polygons",
                            "`features_sf`", task, .JOIN_AGAIN, no_crs_f,
                            of = sprintf("the %d features of `features_sf`", n_f),
                            kind = "features", repaired = repaired, fixable = fixable,
                            s2_forced = forced, advise_predicate = forced)
  }
  # A polygon s2 refuses as the join had it is what the join failed on; the
  # layers as they arrived are asked about only when s2 accepts both
  # repaired.  Asked the other way round, a figure-eight cell beside a
  # feature s2 refuses even repaired (a bowtie) was named as the cause, with
  # the feature's reason quoted as that of a piece of an overlap.
  if (length(pa$rows)) cells(pa, FALSE)
  if (length(fa$rows)) features(fa, FALSE)
  if (length(pb$rows)) cells(pb, TRUE)
  if (length(fb$rows)) features(fb, TRUE)
  invisible(NULL)
}

#' The errors of class "spatialkit_cells_degenerate" and
#' "spatialkit_features_degenerate"
#'
#' Builds and raises the error \code{.stop_cells_degenerate()} or
#' \code{.stop_join_degenerate()} found the cause of.  A function that hands
#' its polygons to another that raised this error calls it again with the
#' fields of that condition and its own words, so the message names the
#' argument the caller passed.
#'
#' @param parent s2's error, quoted in the message.
#' @param rows The rows of \code{layer} s2 refuses.
#' @param reasons s2's reason for each of \code{rows}
#'   (\code{sf::st_is_valid(reason = TRUE)}).
#' @param n The number of polygons s2 was handed.
#' @param caller,layer,task,again,no_crs As for \code{.stop_cells_degenerate()}.
#' @param of The polygons counted, as in "refuses 2 of \code{of}"; \code{NULL}
#'   for "the \code{n} polygons of \code{layer}".
#' @param kind \code{"cells"} for the polygons a call places or measures
#'   things in (the cells' class, with the sentence on corners that differ by
#'   rounding where s2 reads two corners as one point), or \code{"features"}
#'   for the features a join places (the features' class, without it).
#' @param repaired Whether the largest-overlap join repaired these rows into
#'   polygons s2 accepts, and s2 then refused a piece of an overlap.
#' @param also The rows of \code{features_sf} s2 refuses too, as a list of
#'   \code{rows}, \code{reasons} and \code{n}, or \code{NULL}.
#' @param fixable Whether \code{sf::st_make_valid()} makes s2 accept every
#'   row named (features only), so the message offers it.
#' @param s2_forced Whether the call ran on s2 although the session has
#'   \code{sf::sf_use_s2(FALSE)}: the message then does not call s2 the
#'   engine sf uses in degrees.
#' @param advise_predicate Whether the join by \code{predicate} runs on these
#'   layers (with GEOS, s2 off): the message then offers
#'   \code{largest = FALSE} and does not say that join meets the refusal.
#' @param more A sentence added after the remedy, or \code{NULL}: another way
#'   on that the caller offers.
#' @param pick As for \code{.stop_cells_degenerate()}; read for the cells of a
#'   layer with a CRS.
#' @param source Where the cells' corners that s2 reads as one point come
#'   from, the sentence after the one that names them; \code{NULL} for cells
#'   built in a projected CRS and returned in degrees.  \code{""} leaves out
#'   both sentences, for polygons that are not cells (a boundary).
#' @param remedy The remedy sentence, in place of the one built from
#'   \code{layer}, \code{again} and \code{pick}; \code{NULL} for that one.
#'   \code{more} is still added after it.
#' @param site A word stored on the condition as \code{site}, which tells a
#'   caller that raises it again in its own words which polygons were
#'   refused; \code{NULL} for none.
#' @param auto_helps Whether every row refused is a POLYGON or MULTIPOLYGON,
#'   which mode \code{"auto"} of \code{coerce_to_points()} reduces with a
#'   point on the surface, so a caller that raises the error again offers
#'   \code{"auto"}; \code{FALSE} when a row is another type (a collection),
#'   which \code{"auto"} takes the centroid of and refuses again.  Stored on
#'   the condition as \code{auto_helps}; \code{NULL} stores nothing.
#' @return Does not return.  The condition carries \code{rows}, \code{n},
#'   \code{reasons}, \code{parent}, \code{repaired} and \code{s2_forced},
#'   \code{features_rows} when \code{also} is given, \code{site} when
#'   \code{site} is given, and \code{auto_helps} when \code{auto_helps} is
#'   given.
#' @keywords internal
#' @noRd
.raise_cells_degenerate <- function(parent, rows, reasons, n, caller, layer, task, again,
                                    no_crs = FALSE, of = NULL,
                                    kind = c("cells", "features"), repaired = FALSE,
                                    also = NULL, fixable = FALSE, s2_forced = FALSE,
                                    advise_predicate = FALSE, more = NULL,
                                    pick = "the one ensure_projected() picks for it",
                                    source = NULL, remedy = NULL, site = NULL,
                                    auto_helps = NULL) {
  kind <- match.arg(kind)
  said <- trimws(strsplit(conditionMessage(parent), "\n", fixed = TRUE)[[1L]][1L])
  rows_text <- function(r) {
    k <- length(r)
    if (k == 1L) sprintf("row %d", r)
    else if (k <= 6L) sprintf("rows %s and %d", paste(r[-k], collapse = ", "), r[k])
    else sprintf("rows %s and %d more", paste(r[1:6], collapse = ", "), k - 6L)
  }
  # Each row with its own reason, for the rows rows_text() names.
  rows_why <- function(r, why) {
    k <- length(r); shown <- seq_len(min(k, 6L))
    paste0(paste(sprintf("row %d: \"%s\"", r[shown], why[shown]), collapse = "; "),
           if (k > 6L) sprintf("; and %d more", k - 6L) else "")
  }
  k <- length(rows)
  # s2's message for the call is quoted as the rows' own when it is one of
  # their reasons.  It can be another layer's (a feature s2 met first) or a
  # piece of an overlap the call made; each row's own reason is quoted then.
  own <- !isTRUE(repaired) && any(vapply(sub("^Loop [0-9]+: ", "", reasons), function(r)
    nzchar(r) && grepl(r, said, fixed = TRUE), logical(1)))
  named <- if (own) sprintf("%s; \"%s\"", rows_text(rows), said) else rows_why(rows, reasons)
  them <- if (k == 1L) c("it", "it", "it") else c("them", "them", "they")
  told <- if (isTRUE(repaired))
    paste0(sprintf(paste0("The largest-overlap join repaired %s with sf::st_make_valid(), ",
                          "which s2 accepts, and s2 then refused a piece of an overlap of the ",
                          "two layers (\"%s\")"), them[1L], said),
           if (isTRUE(advise_predicate)) ". "
           else sprintf("; the join by `predicate`, which repairs nothing, meets %s as %s arrived. ",
                        them[2L], them[3L]))
  else ""
  other_said <- if (!isTRUE(repaired) && !own) sprintf("s2's error reads \"%s\". ", said)
                else ""
  # s2 gives "duplicate vertex" for two corners it reads as one point; the
  # cause and its usual source are named for those polygons only.
  dup <- kind == "cells" & !identical(source, "") &
    grepl("duplicate vertex", reasons, fixed = TRUE)
  other <- which(!dup)
  cause <- if (!any(dup)) "" else paste0(
    if (!length(other)) {
      if (k == 1L) "It has two corners that s2 reads as one point. "
      else "Each has two corners that s2 reads as one point. "
    } else sprintf(paste0("%d of them %s two corners that s2 reads as one point, and ",
                          "s2 refuses the other %d for another reason (%s: \"%s\"). "),
                   sum(dup), if (sum(dup) == 1L) "has" else "have", length(other),
                   rows_text(rows[other[1L]]), reasons[other[1L]]),
    if (!is.null(source)) source else paste0(
    "Cells built in a projected CRS and returned in longitude and latitude (by ",
    "build_tessellation() or create_voronoi_polygons() with a geographic `crs`) can ",
    "hold such corners: where four seeds stand the same distance from one point, two ",
    "cells meet there at two corners that differ only by rounding, which the ",
    "transform to degrees puts at one point. "))
  also_text <- if (length(also$rows))
    sprintf("s2 also refuses %d of the %d features of `features_sf` (%s). ",
            length(also$rows), also$n, rows_why(also$rows, also$reasons))
  else ""
  fix <- if (!is.null(remedy)) remedy
  else if (kind == "features") paste0(
    if (isTRUE(no_crs))
      sprintf(paste0("%s has no CRS and was read as longitude and latitude: give it its ",
                     "CRS with sf::st_set_crs(), transform it to a projected CRS with ",
                     "sf::st_transform(), and %s."), layer, again)
    else
      sprintf(paste0("Transform %s to a projected CRS, such as the one ensure_projected() ",
                     "picks for it (sf::st_transform()), and %s."), layer, again),
    if (isTRUE(fixable) && !isTRUE(no_crs))
      sprintf(" Or repair %s with sf::st_make_valid(), which s2 accepts, and assign again.",
              them[1L]))
  else if (isTRUE(no_crs))
    sprintf(paste0("%s has no CRS and was read as longitude and latitude: give it its ",
                   "CRS with sf::st_set_crs(), transform it to a projected CRS, such as ",
                   "the one its cells were built in, with sf::st_transform(), and %s."),
            layer, again)
  else
    sprintf(paste0("Transform %s to a projected CRS, such as the one its cells were ",
                   "built in (sf::st_transform()) or %s, and %s."), layer, pick, again)
  # With s2 off the join by `predicate` runs with GEOS, which reads these
  # polygons in degrees.
  if (isTRUE(advise_predicate))
    fix <- paste0(fix, " Or pass largest = FALSE to assign every feature by ",
                  "`predicate` and `tie_break` instead; with sf::sf_use_s2(FALSE), as in ",
                  "this session, that join runs with GEOS.")
  if (length(more)) fix <- paste(fix, more)
  engine <- if (isTRUE(s2_forced))
    "s2, which the package uses here in longitude and latitude whatever sf::sf_use_s2() says,"
  else "s2, the engine sf uses for geometry in longitude and latitude,"
  if (is.null(of)) of <- sprintf("the %d polygons of %s", n, layer)
  msg <- sprintf(paste0(
    "%s(): %s refuses %d ",
    "of %s (%s), so %s cannot run in degrees. %s%s%s%s%s"),
    caller, engine, k, of, named, task, told, cause, other_said, also_text, fix)
  fields <- list(rows = rows, n = n, reasons = reasons, parent = parent,
                 repaired = isTRUE(repaired), s2_forced = isTRUE(s2_forced))
  if (!is.null(also)) fields$features_rows <- also$rows
  if (!is.null(site)) fields$site <- site
  if (!is.null(auto_helps)) fields$auto_helps <- isTRUE(auto_helps)
  stop(do.call(errorCondition, c(list(msg, class = paste0("spatialkit_", kind, "_degenerate"),
                                      call = NULL), fields)))
}


# The areas the smallest-area tie-break compares: `g`, the polygons the tied
# features fall in, rows `rows` of the `n` polygons of `polygons_sf`.  GEOS
# gives a ring whose edges cross the signed area of its lobes, which cancel
# (a bowtie of two 90,000 m^2 lobes measured 0 m^2), so a point in a valid
# cell and in a bowtie went to the bowtie as the smaller, without a word.
# Polygons the validity check refuses are measured on repaired copies, with
# the largest-overlap join's warning, and returned as they arrived.  Called
# under the engine of the join (s2 in degrees); a polygon s2 refuses even
# repaired (a bowtie in degrees) has no area there, and stops the call with
# the error that names it and the transform.  `s2_session` is the session's
# sf::sf_use_s2(), which the error reads to name the engine.
.tie_break_areas <- function(g, rows, n, no_crs = FALSE, s2_session = TRUE) {
  if (!length(g)) return(numeric(0))
  ok  <- tryCatch(suppressWarnings(sf::st_is_valid(g)), error = function(e) NULL)
  bad <- if (length(ok) == length(g)) which(!(ok %in% TRUE)) else seq_along(g)
  if (length(bad)) {
    .warn_and_log(paste0(
      "assign_features_to_polygons(): %d of the %d polygon(s) whose areas the ",
      "smallest-area tie-break compares have invalid geometry (usually a ",
      "self-intersecting ring), so their areas are undefined; they were repaired ",
      "with sf::st_make_valid() to measure them for the tie-break. The polygons are ",
      "not changed and the features are returned with the geometry they arrived ",
      "with; repair them yourself to silence this."),
      length(bad), length(g))
    g[bad] <- .safe_make_valid(g[bad])
  }
  withCallingHandlers(
    suppressWarnings(as.numeric(sf::st_area(g))),
    error = function(e) {
      if (!.is_s2_ring_error(e) || !isTRUE(sf::st_is_longlat(g)) ||
          !isTRUE(sf::sf_use_s2())) return()
      still <- .s2_refused_rows(g)
      if (is.null(still) || !length(still$rows)) return()
      .raise_cells_degenerate(e, rows[still$rows], still$reasons, n,
                              "assign_features_to_polygons", "`polygons_sf`",
                              "the smallest-area tie-break", .JOIN_AGAIN, no_crs,
                              s2_forced = !isTRUE(s2_session))
    })
}


# Exact ties in a largest-overlap join.  `joined` is st_join(largest = TRUE)
# of `f` onto `p_sel`, the polygons' ID in column `join_id`.  The features
# assigned there that overlap two or more polygons are intersected with them
# again; where the largest overlap is shared (areas equal to 9 significant
# digits) the polygon is picked among the equally largest by `tie_break`:
# "smallest_area" takes the smallest polygon, then the one whose bounding-box
# centre is lowest, then leftmost, then the one whose vertices' mean is
# lowest, then leftmost, centres and means within `tol` equal, as for the
# predicate join (.tie_rule_rows()), and gives NA where polygons equal on
# every key are left; "first" takes the first by polygon row.  Returns the
# rows of `joined` that had a tie and the ID each gets (NA for none).
.largest_overlap_ties <- function(f, p_sel, joined, join_id, tie_break, tol = 0) {
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
      vm <- .vertex_means(gp[tied])
      won <- tied[.tie_rule_rows(cbind(ga, ctr[2, ], ctr[1, ], vm[, 2L], vm[, 1L]),
                                 rep(1L, length(tied)), tol)]
      if (length(won) == 1L) won else NA_integer_
    }
    out_rows <- c(out_rows, rows[multi[idx[s[1L], 1L]]])
    out_ids  <- c(out_ids, p_sel[[join_id]][pick])
  }
  list(rows = out_rows, ids = out_ids)
}


#' The candidates the smallest-area rule keeps
#'
#' Each feature's candidates are a group of rows of \code{key}: the area of
#' the polygon (rounded to 9 significant digits by the caller), the
#' bounding-box centre's y and x, and the mean of the vertices' y and x
#' (\code{.vertex_means()}).  The candidates of least area are kept; then,
#' key by key, those within \code{tol} of the least among them.  A key that
#' is not finite counts as \code{Inf}.  A group left with more than one row
#' holds polygons the rule cannot tell apart: copies of one polygon, or
#' polygons that share an area, a bounding-box centre and a vertex mean (a
#' square and a rhombus about its centre).
#'
#' @param key Numeric matrix, one row per candidate, the five keys in order.
#' @param grp The feature of each row.
#' @param tol The margin for the centres and the means (\code{.edge_tol()}).
#' @return Logical vector, one per row: whether the rule keeps it.
#' @keywords internal
#' @noRd
.tie_rule_rows <- function(key, grp, tol = 0) {
  keep <- rep(TRUE, nrow(key))
  for (j in seq_len(ncol(key))) {
    k <- key[, j]
    k[!is.finite(k)] <- Inf
    lo <- stats::ave(ifelse(keep, k, Inf), grp, FUN = min)
    keep <- keep & if (j == 1L) k == lo else k <= lo + tol
  }
  keep
}


#' The mean of the vertices of polygons, as \code{.polygon_sort_keys()} takes it
#'
#' @param g An sfc.
#' @return A matrix, one row per geometry: the mean x and the mean y of its
#'   vertices, each ring's closing vertex left out; \code{Inf} for a
#'   geometry that is empty or not a POLYGON or MULTIPOLYGON.
#' @keywords internal
#' @noRd
.vertex_means <- function(g) {
  out <- matrix(Inf, length(g), 2L)
  ok  <- as.character(sf::st_geometry_type(g, by_geometry = TRUE)) %in%
    c("POLYGON", "MULTIPOLYGON") & !sf::st_is_empty(g)
  if (any(ok)) out[ok, ] <- .polygon_sort_keys(g[ok])[, c("mx", "my"), drop = FALSE]
  out
}


#' The default predicate of the join, with a point on an edge in both polygons
#'
#' \code{sf::st_intersects()}, and for each point feature (POINT or
#' MULTIPOINT) also the polygons it lies within \code{.edge_tol(bb)} of
#' (\code{.edge_join()}; a MULTIPOINT is within the margin of a polygon when
#' one of its points is), so a point on the edge two polygons share is in
#' both whichever side of it rounding puts the point.  Other geometries are
#' joined by \code{sf::st_intersects()} alone.
#'
#' @param x The features (sf).
#' @param y The polygons (sf), in the CRS of \code{x}.
#' @param bb The bounding box the margin is taken from: the polygons'.
#' @return A list, one integer vector of polygon rows per feature, as
#'   \code{sf::st_join()} takes from its \code{join}.
#' @keywords internal
#' @noRd
.intersects_margin <- function(x, y, bb) {
  hits <- unclass(sf::st_intersects(x, y))
  gx   <- sf::st_geometry(x)
  pt   <- which(as.character(sf::st_geometry_type(gx, by_geometry = TRUE)) %in%
                  c("POINT", "MULTIPOINT") & !sf::st_is_empty(gx))
  if (length(pt)) {
    near <- attr(.edge_join(gx[pt], y, bb), "hits")
    hits[pt] <- lapply(seq_along(pt), function(k)
      sort(unique(c(hits[[pt[k]]], as.integer(near[[k]])))))
  }
  lapply(hits, as.integer)
}


#' Points on the edge between cells a profile scored, placed as it placed them
#'
#' On cells whose record of their count sits beside the centres the profile
#' scored (\code{.scored_centres()}), \code{resolution_profile()} gave each
#' row to its nearest centre in the record's CRS, squared distances within
#' \code{.TIE_TOL} of the least, plus the margin \code{.centres_tol()} on
#' the distance (\code{.tie_bound()}), a tie and a tie to the lower-numbered
#' centre (\code{.nearest_centre()}), and the polygons hold a row that is the
#' same distance from two centres in either cell, or in both.  A point the
#' join placed in a cell, that lies within the margin of the edge of that
#' cell and of the cell of the centre the profile's rule gives it, is given
#' that cell.  The margin is \code{.edge_tol()}, plus, for a point the rule
#' gives to another centre than that of the cell it is in, the point's
#' distance from the bisector of the two centres in the join's CRS: a point
#' whose squared distances to two centres tie by that rule can lie
#' farther from their shared edge than \code{.edge_tol()} when the centres
#' are far apart (by about \code{1.25e-8} of their distance midway between
#' them).  A point farther from that cell than the margin keeps the cell it
#' has: it is not on their shared edge, so the cells are not the ones scored
#' there (built in another CRS, say), and the summary's check says so.
#'
#' The cell of a centre is the one cell that holds it, when it holds no
#' other centre: a Voronoi cell holds its own centre and no other, so cells
#' built from other seeds (a subset of the centres) leave most centres
#' without one, and the points tied to those keep their cells.  A centre in
#' no cell stands outside a boundary that clipped its cell away from it;
#' its cell is then the one cell holding no centre whose point on the
#' surface has it as its nearest centre, when no other such cell has.
#'
#' Edges are measured in the CRS of the join: a projected one, or none on
#' either side when the record names none.  Cells in longitude and latitude
#' beside a record in a projected CRS arrive here in that CRS, where
#' \code{assign_features_to_polygons()} joins them.  Polygons in longitude
#' and latitude (beside a record in degrees), features without a CRS beside
#' a record with one, or the reverse, are left as the join placed them.
#'
#' @param f The join's copy of the features, in the join's CRS, one row per
#'   row of \code{features}.
#' @param p The polygons, in the join's CRS, with their IDs in \code{id_col}.
#' @param id_col The ID column of \code{p}.
#' @param joined The join, one row per feature: \code{..pre_join_row_id}
#'   (the row of \code{f}) and the polygon's ID in \code{id}.
#' @param id The name of the ID column of \code{joined}.
#' @param rec The record of the count the polygons carry.
#' @param features The features as the caller passed them.
#' @return A list: \code{rows}, the rows of \code{joined} whose cell
#'   changes, and \code{ids}, the ID each gets.
#' @keywords internal
#' @noRd
.scored_edge_cells <- function(f, p, id_col, joined, id, rec, features) {
  none <- list(rows = integer(0), ids = p[[id_col]][integer(0)])
  cen <- .scored_centres(rec)
  if (is.null(cen)) return(none)
  crs_j <- sf::st_crs(p)
  crs_r <- rec$crs
  crs_f <- sf::st_crs(features)
  if (!is.na(crs_j) && isTRUE(sf::st_is_longlat(crs_j))) return(none)
  if (!is.na(crs_r) && isTRUE(sf::st_is_longlat(crs_r))) return(none)
  if (is.na(crs_r) != is.na(crs_j) || is.na(crs_r) != is.na(crs_f)) return(none)

  # The points the join placed in a cell.
  ids_p <- p[[id_col]]
  at    <- match(joined[[id]], ids_p)
  src   <- joined[["..pre_join_row_id"]]
  gf    <- sf::st_geometry(f)
  point <- sf::st_geometry_type(gf, by_geometry = TRUE) == "POINT" & !sf::st_is_empty(gf)
  rows  <- which(!is.na(at) & point[src])
  if (!length(rows)) return(none)

  # Their distance from their own cell's edge.
  gp  <- sf::st_geometry(p)
  tol <- .edge_tol(sf::st_bbox(gp))
  edge <- suppressWarnings(sf::st_boundary(gp))
  d    <- as.numeric(sf::st_distance(gf[src[rows]], edge[at[rows]], by_element = TRUE))

  # Their nearest centre, from their coordinates in the record's CRS, as the
  # profile measured them.  Every point is searched: one inside the band
  # where two centres tie can lie farther from the edge than `tol`.
  pts <- features[src[rows], , drop = FALSE]
  if (!is.na(crs_r) && !isTRUE(crs_f == crs_r)) pts <- sf::st_transform(pts, crs_r)
  xy <- unname(sf::st_coordinates(pts)[, 1:2, drop = FALSE])
  ok <- is.finite(xy[, 1L]) & is.finite(xy[, 2L]) & is.finite(d)
  rows <- rows[ok]; d <- d[ok]
  if (!length(rows)) return(none)
  g <- .nearest_centre(xy[ok, 1L], xy[ok, 2L], cen[, 1L], cen[, 2L],
                       use_fnn = .has_backend("FNN"))

  # The cell of each centre: the one cell that holds it and no other centre.
  where <- sf::st_as_sf(data.frame(x = cen[, 1L], y = cen[, 2L]), coords = c("x", "y"),
                        crs = crs_r)
  if (!is.na(crs_r) && !isTRUE(crs_r == crs_j)) where <- sf::st_transform(where, crs_j)
  hit  <- sf::st_intersects(where, gp)
  cell <- rep(NA_integer_, nrow(cen))
  one  <- lengths(hit) == 1L
  cell[one] <- vapply(hit[one], function(h) h[1L], integer(1))
  cell[cell %in% cell[duplicated(cell) & !is.na(cell)]] <- NA_integer_
  # A centre in no cell lies outside a boundary that clipped its cell away
  # from it.  Its cell is the one, among the cells that hold no centre, whose
  # inside is nearer it than any other centre: a point on that cell's surface
  # has it as its nearest centre, and no other such cell has.
  lost <- which(lengths(hit) == 0L)
  if (length(lost)) {
    empty <- setdiff(seq_along(gp), unlist(hit))
    if (length(empty)) {
      inside <- suppressWarnings(sf::st_point_on_surface(gp[empty]))
      if (!is.na(crs_r) && !isTRUE(crs_r == crs_j)) inside <- sf::st_transform(inside, crs_r)
      ixy <- unname(sf::st_coordinates(inside)[, 1:2, drop = FALSE])
      fin <- is.finite(ixy[, 1L]) & is.finite(ixy[, 2L])
      k <- rep(NA_integer_, length(empty))
      k[fin] <- .nearest_centre(ixy[fin, 1L], ixy[fin, 2L], cen[, 1L], cen[, 2L],
                                use_fnn = .has_backend("FNN"))
      own <- !is.na(k) & k %in% lost & !(k %in% k[duplicated(k) & !is.na(k)])
      cell[k[own]] <- empty[own]
    }
  }

  # A point goes to its centre's cell when it lies on the edge of both
  # cells.  The margin is `tol`, and for a point the rule gives to a centre
  # that ties with a strictly nearer one, the centre of the cell it is in,
  # also the point's distance from the bisector of the two: how far inside
  # that cell the tie reaches.  It is measured in the record's CRS, where
  # the two tie, and put in the units of the join's CRS by the ratio of the
  # centres' distance apart in the two.  For any other point it is `tol`
  # alone, so a point placed in another cell by cells drawn in another CRS
  # is not moved.
  to   <- cell[g]
  move <- which(!is.na(to) & to != at[rows])
  if (!length(move)) return(none)
  h    <- match(at[rows[move]], cell)
  gm   <- g[move]
  pxy  <- xy[ok, , drop = FALSE][move, , drop = FALSE]
  d2g  <- (pxy[, 1L] - cen[gm, 1L])^2 + (pxy[, 2L] - cen[gm, 2L])^2
  d2h  <- (pxy[, 1L] - cen[h, 1L])^2 + (pxy[, 2L] - cen[h, 2L])^2
  c_r  <- sqrt((cen[gm, 1L] - cen[h, 1L])^2 + (cen[gm, 2L] - cen[h, 2L])^2)
  wxy  <- unname(sf::st_coordinates(where)[, 1:2, drop = FALSE])
  c_j  <- sqrt((wxy[gm, 1L] - wxy[h, 1L])^2 + (wxy[gm, 2L] - wxy[h, 2L])^2)
  band <- (d2g - d2h) / (2 * c_r) * (c_j / c_r)
  band[!(!is.na(h) & is.finite(band) & d2h < d2g)] <- 0
  margin <- tol + band
  keep <- which(d[move] <= margin)
  move <- move[keep]; margin <- margin[keep]
  if (!length(move)) return(none)
  d_to <- as.numeric(sf::st_distance(gf[src[rows[move]]], gp[to[move]], by_element = TRUE))
  move <- move[which(d_to <= margin)]
  list(rows = rows[move], ids = ids_p[to[move]])
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
#' @param deff Design effect for that cell, \eqn{\geq 1}, or \code{NA} when
#'   it is unknown (a Kish ICC that cannot be estimated).
#' @return The standard error, or \code{NA_real_} when the cell carries no
#'   information about \eqn{\sigma} (\code{deff >= n}: complete redundancy)
#'   or the design effect is unknown.
#' @keywords internal
#' @noRd
.se_with_deff <- function(s, n, deff, correct_s2 = TRUE) {
  if (!is.finite(s) || !is.finite(n) || n <= 1L) return(NA_real_)
  # An unknown design effect is not a design effect of 1: the uncorrected SE
  # returned here would read as an estimate the data do not support.
  if (is.na(deff)) return(NA_real_)
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


#' Standard deviation of a cell's values, in any unit
#'
#' \code{sd()} squares the deviations from the mean.  The squares overflow
#' for values above about 1e154, where \code{sd()} is \code{Inf} (or
#' \code{NaN} where the mean overflows too), and fall below the smallest
#' normal double for a spread below about 1e-154, where they lose their bits
#' and \code{sd()} is 0 or off.  Where the raw standard deviation is not
#' finite or is below \code{.POOL_Y_MIN}, it is taken on
#' \code{.spread_scaled()} of the values and multiplied by its scale, a
#' power of two; elsewhere it is \code{sd()} of the values as they stand.
#' It is \code{Inf} only where the standard deviation is beyond the largest
#' double.  The raw reading is \code{sd(x, na.rm = TRUE)} even where
#' \code{x} has no \code{NA}: \code{sd(x)} takes another path through
#' \code{cov()} and can differ from it in the last bit.
#'
#' @param x Finite values or \code{NA}, at least two that are not \code{NA}.
#' @return The standard deviation of the values that are not \code{NA}.
#' @keywords internal
#' @noRd
.cell_sd <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.finite(s) && s >= .POOL_Y_MIN) return(s)
  sc <- .spread_scaled(x[!is.na(x)])
  sd(sc$y) * sc$scale
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
#' Rows with a coordinate that is not finite or lies beyond
#' \code{.POOL_XY_MAX} (an empty point among them) have no location to
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
    located <- is.finite(xy[, 1]) & is.finite(xy[, 2]) &
      abs(xy[, 1]) <= .POOL_XY_MAX & abs(xy[, 2]) <= .POOL_XY_MAX
    if (!all(located)) xy <- xy[located, , drop = FALSE]
  }
  n_i <- nrow(xy)
  if (is.null(cor_fn) || is.null(n_i) || n_i <= 1L)
    return(c(rbar = 0, df_ratio = 1))
  if (n_i > max_n) {
    # Drawn from the rows in coordinate order (.canonical_order()), so the
    # subsample is the same for the cell's rows in any order.  Drawn from the
    # rows as given, the same 1500-point cell got a standard error of 0.261
    # with its rows as given and 0.255 with them reversed.  Rows at one
    # location are interchangeable for the correlation, so the coordinates
    # are the only keys.
    ord <- .canonical_order(xy[, 1], xy[, 2], xy[, 1], xy[, 2])
    cleanup <- .with_seed(seed)
    on.exit(cleanup(), add = TRUE)
    xy <- xy[ord[sample.int(n_i, max_n)], , drop = FALSE]
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
#' number, and to \code{"kish"} (or, when no ICC is positive, to nothing) for a
#' \code{"variogram"} request with no usable model.  The class is
#' \code{"spatialkit_deff_fallback"} in both cases.
#'
#' @keywords internal
#' @noRd
.warn_deff_fallback <- function(fmt, ..., data = list()) {
  .warn_classed("spatialkit_deff_fallback", fmt, ..., data = data)
}


# The reason .estimate_icc() gives for a column that does not vary.  Its ICC
# is NA (the between and within variances are both 0), and its standard
# errors are 0 whatever the correlation, so no number is withheld for it.
.icc_reason_constant <- "the column does not vary"


#' Summarize features by polygon/cell ID
#'
#' Aggregates an sf point dataset into one row per cell. By default computes
#' counts and means, but the aggregation function is configurable. It also
#' returns a pooled value for each cell, its mean pulled toward the mean of
#' the cells around it (`estimate = "pooled"`, the default; `estimate =
#' "mean"` leaves it out), and can return a range for one further
#' observation.
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
#' `..sd_*` is the standard deviation of the cell's values in the units of
#' the column, at any size. Where the raw squares overflow or lose their bits
#' (a standard deviation that is not finite or is below 1e-150), it is taken
#' on the values centred on their median and divided by powers of two,
#' which keeps their spread. A cell whose
#' standard deviation is beyond the largest double (about 1.8e308) has
#' `..sd_*` of `Inf` and standard errors and interval of `NA`, with a warning
#' of class `"spatialkit_sd_overflow"` (the condition carries `columns` and
#' `cells`, the number of such cells of each column).
#'
#' # Spatial autocorrelation and standard-error bias
#'
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
#' require enough cells and observations for a stable ICC estimate. A column
#' whose ICC cannot be estimated at all gets `NA` standard errors (see
#' `deff`).
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
#' # What the standard error estimates
#'
#' What the `..se_*` columns estimate follows `deff`. At the default,
#' `deff = 1`, `sd / sqrt(n)` is the standard error of **the cell's own
#' mean**, the mean of the units the cell holds, and it is the right one when
#' the cell's observations are spread through the cell. Under
#' `deff = "kish"`, `"variogram"` or a number other than 1 the columns are the
#' standard error of the cell mean **as an estimate of the population (grand)
#' mean**: the unconditional quantity, in which the cell's own realised
#' deviation is part of the error, widened for the correlation inside the
#' cell. That is the right quantity when cells are treated as samples from a
#' common population: measured 95% interval coverage of the grand mean is
#' 0.95 with `deff = "kish"` (and 0.29 with the naive SE) on exchangeable
#' within-cell correlation, and 0.93 with `deff = "variogram"` on a simulated
#' Gaussian field.
#'
#' The corrected standard error is **not** that of the cell's own mean (the
#' block average over that cell), which is what a cell-level map or a
#' regression on cell values usually wants. For that quantity the naive
#' `sd / sqrt(n)` is the
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
#' # Design effects and variable types
#'
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
#' When that variogram finds no spatial structure (the pure nugget
#' [estimate_sac_range()] returns with that verdict), every column gets a
#' design effect of 1, the predictors included, and a warning of class
#' `"spatialkit_deff_no_structure"` names the predictor columns: the verdict
#' is about the variable the variogram was fitted to, not about them.
#' `deff = "kish"` gives each predictor its own ICC.
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
#' # Confidence intervals
#'
#' With `conf_level` set, every numeric response and predictor column gains
#' four more columns: `..neff_*`, that column's effective sample size in the
#' cell (its non-missing count over its design effect, the per-column version
#' of `cell_weight`); `..df_*`, the degrees of freedom the interval uses; and
#' `..ci_lo_*` / `..ci_hi_*`, a t interval for a mean. Which mean follows
#' `deff`, as for `..se_*`: at the default, `deff = 1`, **the cell's own
#' mean**, when its observations are spread through the cell; under
#' `deff = "kish"`, `"variogram"` or a number other than 1, the **cell mean as
#' an estimate of the population (grand) mean**. Everything in "What the
#' standard error estimates" applies to it, including that the corrected
#' interval is not one for the cell's own block average. The interval is
#' `mean +/- qt((1 + conf_level) / 2, df) * se`, centred on the plain mean of
#' the column's non-missing values whatever `agg_funs` computes, and is `NA`
#' wherever the standard error is (a single observation; complete redundancy
#' under `deff`). `..neff_*` and `..df_*` are `NA` where the column has one
#' non-missing value or none, like the standard error; `cell_weight` still
#' counts such a cell (1, or `1 / deff` for a numeric `deff`), so the two
#' differ there.
#'
#' At no setting is `..ci_*` a range for one observation: it is an interval
#' for a mean, and it narrows as a cell's count grows where the spread of
#' single observations does not. The interval for one further observation in
#' a cell has its own argument and columns, `obs_level` with
#' `estimate = "pooled"` and `..obs_lo_*` / `..obs_hi_*` (see "Intervals:
#' which one answers which question"). In the simulation described there,
#' `..ci_*` at `conf_level = 0.90` held 58.9 percent of the observations held
#' out among the others, and the interval for one observation at
#' `obs_level = 0.90` held 89.6 percent. Away from the observations no
#' interval here is calibrated.
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
#' # Pooled values
#'
#' With `estimate = "pooled"` every numeric column also gets a pooled value.
#' That is the default; `estimate = "mean"` leaves the pooled columns out
#' and returns the frame of earlier versions.
#' A cell's value, plain or pooled, estimates the cell's *level*: the
#' expected value of one further observation located in the cell the way its
#' observations were. It is not the average of a surface over the cell's
#' area, which would give weight to ground that holds no observation.
#'
#' For one column, write `n`, `S` and `m = S / n` for the count, the sum and
#' the mean of a cell's rows, `N` for the rows in all and `K` for the cells
#' that hold one. Then
#'
#' * `sigma2 = sum((y - m)^2) / (N - K)` is the spread of the observed values
#'   about their cell's mean;
#' * the target `t` of a cell is the mean of the rows of its neighbours,
#'   `sum(S) / sum(n)` over them. The neighbours are the 8 cells whose
#'   centroids are nearest, the centroid of a cell being the mean of the
#'   coordinates of its rows. The polygons are not read, so two cells are
#'   neighbours whether or not they touch, and a cell is never its own
#'   neighbour;
#' * `tau2 = max(0, (sum(n * (m - t)^2) - K * sigma2) / N)` is the mean
#'   squared distance between a cell's level and its target, the target's own
#'   error included;
#' * the own weight is `w = tau2 / (tau2 + sigma2 / n)`;
#' * the pooled value is `w * m + (1 - w) * t`.
#'
#' The own weight is the share of the pooled value that is the cell's own
#' mean. It grows with the cell's count and with `tau2`. At 1 the pooled
#' value is the plain mean; at 0 it is the target and holds none of the
#' cell's own values. Pooled values are for maps and for predicting
#' observations. `..pooled_w_*` is not a regression weight. The plain
#' columns are the same under either `estimate`, and they are the ones to
#' use, with `cell_weight`, in a regression on cell values: each column is
#' pulled toward its neighbours by its own amount, so a regression of pooled
#' responses on pooled predictors is not a regression of the levels.
#'
#' The cells can be any polygons the rows were assigned to, the package's
#' own or ones the caller brings, with IDs of any type: only each row's cell
#' and location are read. Each column is pooled on its own rows, those with
#' a value, a cell and a location, with its own centroids, neighbours,
#' `sigma2` and `tau2`; a
#' missing value in one column changes nothing in another. A cell with one
#' such row gets a value. A cell with none gets `NA`, never a value made from
#' its neighbours. The pooled columns do not read `agg_funs` or `deff`: the
#' mean that is pooled is the arithmetic mean whatever `agg_funs` computes,
#' and no design effect enters `sigma2`, `tau2`, the weights or the interval
#' of `obs_level`.
#'
#' Three cases change what the numbers rest on, and each is a warning (see
#' "Conditions"). With 2 to 9 cells the neighbours of a cell are all the
#' others, so its target is the mean of every row outside it, and the overall
#' mean lies between that target and the cell's own mean. A cell whose own
#' weight is below its share of the rows, `n / N`, then gets a pooled value
#' on the far side of the overall mean from its plain mean, and two cells can
#' change order. When few cells hold two rows, `sigma2` is estimated from few
#' values and every weight rests on it: `df` in `attr(, "pooled")` is the
#' number of values beyond one per cell. And when `tau2` comes out at 0,
#' every own weight is 0 and each pooled value is its target; with two
#' cells that is the other cell's mean, so the two values change places.
#'
#' Two candidate neighbours are a tie when their squared distances differ by
#' no more than one part in ten million. A tie goes to the centroid with the
#' lower x, then the lower y, coordinates within that part of the distance
#' counting as equal, and then to the cell whose ID sorts first. The
#' centroids are summed with the rows in coordinate order, so a cell's
#' neighbours do not depend on the order of the rows, and on gridded
#' locations they do not depend on the unit or the origin of the coordinates
#' either. The values agree to rounding under another order of the rows,
#' since a cell's sums are taken in the order of its rows, and under another
#' unit or origin where three to eight cells' centroids coincide, since the
#' rounding then sets the order in which a cell's neighbours on that point
#' are summed. A cell whose centroid shares its point with as many other
#' cells as it has neighbours, or more, takes its neighbours from them, the
#' first by ID, without a search. Centroids count as one point there when
#' they are equal, or when they lie within a ten-millionth of the largest
#' coordinate of one another, in x and in y. Where other centroids that close
#' to them chain on further than that, the chain is taken again at half that
#' reach, and so on, until its parts span no more than the reach they were
#' taken at. So those cells keep their neighbours and their values in another
#' unit, and with the origin moved by little against the coordinates. On such
#' a point, cells whose centroids lie that close without coinciding take each
#' other by ID as well, not by distance. The reach follows the largest
#' coordinate, so with the origin moved far, to the corner of the data say,
#' such centroids can lie beyond it, and their cells then take the nearest.
#' Distances are measured in the CRS named under `estimate`, so the same
#' layer in another projection can still give a cell other neighbours, and
#' so can a lattice whose spacing is below about a ten-millionth of its
#' coordinates.
#'
#' The figures on this page that describe pooled values and the two
#' intervals come from one simulation, `dev/sim-pooled-values.R` in the
#' package's repository, whose header holds its settings and seeds. It draws 200
#' layers on a square of side 1,000. Each has 2,100 clustered locations
#' (four fifths of them round 25 parent points with a standard deviation of
#' 40, the rest uniform) and a response that is a field with covariance
#' `exp(-d / 100)` plus noise of variance 0.3. The 300 locations nearest a
#' random centre are held out as a block, 300 of the others are held out at
#' random, and the 1,500 left are the layer: 333 cells at 4.5 observations a
#' cell from [resolution_profile()] with `partition = "density"`,
#' `select_resolution(criterion = "size")`, [get_voronoi_seeds()],
#' [build_tessellation()] and [assign_features_to_polygons()], summarised
#' with `estimate = "pooled"`. A held-out observation is given the value of
#' its cell by [cell_values_at()]. At the observations held out at random
#' the root mean squared error of the pooled value was 0.96 of the plain
#' mean's: 0.84 in cells that hold one of the layer's observations, 0.95 in
#' cells of 2 to 4, 0.98 in cells of 5 to 9 and 0.99 in cells of 10 or
#' more. In the held-out block it was 0.95.
#'
#' `structured` in `attr(, "pooled")` says whether these cells predict a
#' left-out observation better than one mean for all rows. It is `FALSE`
#' when the leave-one-out error of one mean (`loo_one_cell`) is within 2
#' percent of that of the pooled values (`loo`) or below it, or when the F
#' test of the cell means against one mean has `p >= 0.05`
#' (`structure_p`). `FALSE` raises a warning of class
#' `"spatialkit_no_structure"`, and the values are returned. On a response
#' without spatial structure they still make a smooth map, because each is
#' mostly the mean of its neighbours. With white noise as the response on
#' the simulation's 200 layers the median own weight was 0.11, `tau2` was 0
#' on 20 layers, the warning was raised on all 200, and a cell's pooled
#' value had a correlation of 0.79 with the mean of the pooled values of its
#' 8 nearest cells, where the plain means had one of -0.02.
#'
#' The pooling model implies an error variance for a pooled value as an
#' estimate of the cell's level, `(1 - w) * tau2`. It is not returned, and
#' neither is an interval for the level, because no check has confirmed an
#' interval built from it; `sigma2`, `tau2` and `..pooled_w_*` are in the
#' result for a reader who wants to form it. In the simulation, with a
#' cell's level taken as the mean of the noise-free field over further
#' locations drawn in it like the data's, the squared error of the pooled
#' value was 0.98 of `(1 - w) * tau2` among the observations and 6.9 times
#' it in the held-out block: away from the observations it is too small.
#'
#' Beside a sharp boundary the pooled value can be further from a cell's
#' level than the plain mean, because its target takes in cells across the
#' boundary. A second simulation measured this on fields whose level jumps
#' at county boundaries: one level for each county, alone and with a smooth
#' field added, on two layouts of real locations (892 and 2,634 points, each
#' in a county), with noise of 5, 20 or 50 percent of an observation's
#' variance, 200 layers of each, five random folds, and on each fold the
#' cells the package's defaults build from the training observations. A
#' held-out observation is a boundary row when every training observation
#' of its cell lies in its own county and more than half of those of the
#' cell's neighbours lie in other counties. Boundary rows were 21.75 to
#' 28.04 percent of the held-out observations. Over them the root mean
#' squared error of the pooled value against the noise-free field, as a
#' geometric mean over the layers of its ratio to the plain mean's, was
#' from 9.16 percent lower to 56.36 percent higher. On county levels alone
#' it was 43.21 and 56.36 percent higher with 5 percent noise, 11.23 and
#' 20.58 percent higher with 20 percent, and from 3.79 percent lower to 1.47
#' percent higher with 50 percent; with the smooth field added it was from
#' 9.16 percent lower to 1.59 percent higher. Over all held-out observations
#' the pooled value's error was the lower of the two in every setting, by
#' 1.75 to 18.73 percent.
#'
#' # Intervals: which one answers which question
#'
#' Two families of interval columns can be in the frame. Each has its own
#' argument, and neither changes the other.
#'
#' | Family | Columns | Target | Centre | Asked for by |
#' |---|---|---|---|---|
#' | observation | `..obs_lo_*`, `..obs_hi_*` | one further observation in the cell | pooled value | `obs_level` |
#' | mean | `..ci_lo_*`, `..ci_hi_*` | a mean, never one observation | plain mean | `conf_level` |
#'
#' The observation of the first is one located in the cell the way the
#' data's were. `obs_level` needs `estimate = "pooled"`: without it the call
#' is an error, not an interval round the plain mean. The mean the second
#' covers follows `deff`: at `deff = 1` the cell's own mean, when its
#' observations are spread through the cell; under `deff = "kish"`,
#' `"variogram"` or a number other than 1 the cell mean as an estimate of
#' the population mean (see "Confidence intervals", which also describes
#' `..neff_*` and `..df_*`). No interval and no standard error is returned
#' for a cell's level as the pooled value estimates it.
#' `attr(, "intervals")` says the same of the frame in hand, one row per
#' family present.
#'
#' The interval for one observation is the pooled value plus two quantiles
#' of the leave-one-out residuals of the pooled values. The residual `e` of
#' a row is its value minus the pooled value of its cell computed without
#' it: the cell's mean and own weight are computed again without the row,
#' and `sigma2`, `tau2` and the targets are held. With
#' `p = (1 - obs_level) / 2`,
#' `..obs_lo = pooled + quantile(e, p)` and
#' `..obs_hi = pooled + quantile(e, 1 - p)`.
#'
#' * It is a statement about observations on average: of further
#'   observations located like the data, about `obs_level` of them fall
#'   inside the interval of their cell. It does not hold cell by cell (the
#'   second simulation below gives the coverage in cells of one
#'   observation).
#' * The two offsets are the same in every cell. The interval's centre
#'   varies over the map and its width does not. It does not narrow where a
#'   cell holds more rows: the offsets are quantiles of the errors of
#'   predicting single observations, and the spread of observations inside a
#'   cell is in every one of them.
#' * It is wide. One observation is uncertain.
#' * Only the level 0.90 has been checked. Another level is computed the
#'   same way.
#' * Each residual is the error of predicting an observation from the other
#'   rows of its cell and from the cells around it, so the interval
#'   describes locations among the observations. A cell's value and interval
#'   describe the part of the cell where its observations lie, and a Voronoi
#'   cell can reach far past them toward the edge of the study area. The
#'   interval is too narrow for a cell with no observation, a gap in the
#'   layer or ground beyond the layer's edge.
#'
#' In the simulation described under "Pooled values", the interval at
#' `obs_level = 0.90` held 89.6 percent of the observations held out at
#' random among the others (83.7 to 94.0 percent in single layers) and 75.1
#' percent of those in the held-out block. `..ci_*` at `conf_level = 0.90`,
#' read as a range for one observation, held 58.9 and 49.2 percent of them,
#' over the held-out observations whose cell has that interval.
#' [cell_values_at()] returns each location's distance from the nearest
#' observation, and its help gives the coverage by that distance.
#'
#' The second simulation described under "Pooled values" read the
#' interval at `obs_level = 0.90` on smooth fields of two ranges and on
#' county levels with a smooth field, with 20 or 50 percent noise, on both
#' layouts: twelve settings. Setting by setting it held 89.38 to 89.79
#' percent of the held-out observations. Counted by the training
#' observations in an observation's cell (1, 2 to 4, 5 to 9, 10 or more), it
#' held at least 85.65 percent in every class but one: on county levels with
#' a smooth field, on the layout of 892 points with 20 percent noise, it
#' held 83.74 percent of the 4,878 observations whose cell held one training
#' observation. On fields of that kind it is not a 90 percent interval in
#' such cells: there it held 83.74 to 87.83 percent with 20 or 50 percent
#' noise, and 82.64 and 82.94 percent with 5 percent.
#'
#' Away from the observations no interval here is calibrated.
#'
#' # Conditions
#'
#' What the pooled columns and the count record have to say is raised as R
#' warnings with a class, so `tryCatch(<class> = )` and
#' `withCallingHandlers()` catch each without matching its text. A warning
#' of the pooled columns is raised once per call and reason, for all the
#' columns it applies to: its text names them and the condition carries them
#' in `columns`, in the order of the pooled columns. A field said to be by
#' column is a vector named by those columns, and the numbers of each column
#' are in `attr(, "pooled")`.
#'
#' * `"spatialkit_pool_unavailable"`: the pooled columns of a column are
#'   `NA`. The plain columns are returned, and this is never an error.
#'   `reason` is `"no coordinates"` (a layer without geometry),
#'   `"no located value"` (no row with a value, a cell and a location),
#'   `"no cell holds two values"` (the condition also carries `n` and
#'   `cells` by column), `"values out of range"` (values so large, or so
#'   close to their median, that the squares of many such values can
#'   overflow or lose their precision; the message gives the bounds, and the
#'   column rescaled by a power of ten is pooled) or `"work limit"`. The last
#'   is the case of many cells without the 'FNN' package: the neighbour
#'   search then compares every pair of cells and is not made above
#'   `getOption("spatialkit.exact_pairs_max", 2e9)` pairs, which is about
#'   44,700 cells (the condition also carries `work`, `limit`, and `cells`
#'   by column). With 'FNN' there is no such limit, and the result is the
#'   same.
#' * `"spatialkit_pool_degenerate"`: the pooled values were computed and
#'   carry no pooling. `reason` is `"one cell"` (the pooled value is the
#'   plain mean; also `n` by column), `"no spread inside any cell"` (the
#'   values do not vary inside any cell to their precision, so every own
#'   weight is 1; also `constant` by column) or `"no signal between
#'   cells"` (`tau2` is 0 and every own weight is 0; also `sigma2` by
#'   column).
#' * `"spatialkit_pool_few_cells"`: 2 to 9 cells, so each target is the mean
#'   of all the other cells, and a cell whose own weight is below its share
#'   of the column's rows gets a pooled value on the far side of the overall
#'   mean from its plain mean (`cells`, and `own_weight`, the smallest, by
#'   column). The text also names each column with two cells and `tau2` at
#'   0, which has `cells` 2 and `own_weight` 0: each cell's pooled value is
#'   then the other cell's mean, so the two values change places.
#' * `"spatialkit_pool_few_df"`: `sigma2` rests on fewer than 30 values
#'   beyond one per cell (`needed`, and `df`, `n` and `cells` by column).
#' * `"spatialkit_no_structure"`: `structured` is `FALSE` (see "Pooled
#'   values"). The text says which of the two tests failed, and `beats` is
#'   `TRUE` for a column whose pooled values beat one mean while the F test
#'   does not tell the cell means from one mean (also `criterion = "loo"`,
#'   `estimate = "pooled"`, and by column `n`, the cells, `one_cell`, `loo`
#'   and `structure_p`).
#' * `"spatialkit_pool_rows_skipped"`: rows with a value and a cell have no
#'   location: an empty point, a coordinate that is not finite or lies
#'   beyond 1e150, as given or in the CRS the distances are measured in, or
#'   in a geographic CRS a vertex at a latitude beyond 90 degrees.
#'   They count in `n` and the plain columns and take no part in the pooled
#'   ones (`n`, `cells`).
#' * `"spatialkit_repeated_locations"`, with `obs_level`: rows with a value
#'   share a location. A residual leaves one row out and keeps the others at
#'   its location, so the interval then describes a further observation at a
#'   place already observed (`n` and `total` by column).
#' * `"spatialkit_few_residuals"`, with `obs_level`: fewer than five
#'   residuals are expected beyond each end of the interval, which at 0.90
#'   is fewer than 100 residuals (`needed`, `level`, and `n` by column).
#' * `"spatialkit_estimate_mismatch"` and `"spatialkit_cells_differ"`: the
#'   two checks of the count record, described under "Value".
#'
#' The two about few cells and few values are raised for columns whose
#' weights were estimated, and the two of `obs_level` for columns whose
#' interval is read from residuals, so a constant column raises
#' `"spatialkit_pool_degenerate"` and no other.
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
#'       'gstat' is available. The distances are taken in x and y: a layer
#'       whose features, or the points taken from them, do not all have the
#'       same dimensions (some with a Z or M coordinate) is read without its
#'       Z and M. Exponential, spherical and Gaussian models are
#'       supported, with a nugget and with several structured components
#'       (each weighted by its partial sill) and with gstat's 2-D geometric
#'       anisotropy (`vgm(..., anis = c(angle, ratio))`), applied as gstat
#'       applies it. A request that cannot be met falls back to `"kish"`,
#'       with a warning of class `"spatialkit_deff_fallback"` that gives the
#'       reason and the ICCs: a model of any other family (named in the
#'       warning), or no usable model (none supplied and none could be
#'       estimated, or a rejected fit that could not be replaced by an
#'       estimate). The standard errors, intervals and `cell_weight` are then
#'       Kish's, the uncorrected ones where Kish's own ICC is 0 and `NA`
#'       where it cannot be estimated; the
#'       `deff_method` column says which on every row (see "Value"). The
#'       variogram route needs an identified range, which a layer does not
#'       always have, so read `deff_method` for the fallback. A layer whose
#'       geometries carry M (XYM or XYZM), which GEOS does not read, is
#'       read in x and y too. A model with
#'       no structured component (a pure
#'       nugget) implies that distinct observations are uncorrelated, so it
#'       is applied as a design effect of 1 in every cell, not treated as a
#'       fallback. That includes the model [estimate_sac_range()] returns
#'       when it finds no spatial structure (an `NA` range whose
#'       `rejected_reason` says so): the answer is then a design effect of
#'       1, `deff_method` reads `"variogram"`, and no fallback warning is
#'       raised. Predictor columns summarised under that verdict get the same
#'       design effect of 1, and a warning of class
#'       `"spatialkit_deff_no_structure"` names them (fields `columns` and
#'       `reason`). Points with empty
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
#'       An ICC needs at least 2 cells with 2+ observations of the column and
#'       at least 2 residual degrees of freedom (`N - k >= 2`). Where a column
#'       has fewer (all its values in one cell, say), its ICC cannot be
#'       estimated: it is `NA`, `attr(, "icc")$reason` says why, and the
#'       column's standard errors, `..neff_*` and `..ci_*` are `NA` in every
#'       cell of two or more observations (`..df_*` stays `n - 1`), as is
#'       `cell_weight` when that column sets it, with a warning of class
#'       `"spatialkit_icc_unavailable"` (fields `columns` and `reason`). Rows
#'       with no cell id (kept by `keep_unassigned = TRUE` in
#'       [assign_features_to_polygons()]) give a summary row of their own but
#'       no cell to the ICC, so a column whose values all fall in them, or
#'       differ only there, has an ICC that cannot be estimated. An
#'       estimate at or below 0 is another case: the ICC is 0 and the
#'       column's standard errors are the uncorrected ones. A column that
#'       does not vary (all its values equal, those of rows with no cell id
#'       included) has an ICC of `NA` and standard errors of 0, with no
#'       warning: its mean has no error whatever the correlation. A column
#'       whose values differ only by their rounding (0.1 beside 0.3 / 3) is
#'       read the same way, as one that does not vary, at any size of its
#'       values; its standard errors are 0 in a cell whose values are equal
#'       and of the size of that rounding in one whose values differ. The
#'       `"deff_applied"` attribute is attached when any ICC is positive, so
#'       when all are 0 or `NA` there is none.}
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
#'   least 2 when `deff = "variogram"`; anything else is an error. The
#'   subsample is drawn with a fixed seed from the cell's points in
#'   coordinate order, so it does not depend on the row order of
#'   `assigned_points_sf`, and the caller's random-number stream is left
#'   where it was.
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
#'   spread of planar-to-geodesic area ratios over the cells (over 200 of
#'   them when there are more, taken at even spacing along the cells sorted
#'   by their bounding boxes in that CRS and then by the means of their
#'   vertices, together with the cells that reach furthest along x and y
#'   either way, so the order of the rows does not matter): a density is a
#'   comparison between cells and is meaningless where the map scale differs
#'   from one cell to the next.  Inside a UTM zone the spread is under 0.3
#'   percent and the request goes through; the conterminous United States
#'   forced into one zone (14 percent), or a few degrees of latitude in Web
#'   Mercator (4 percent at 48N), does not.  Exact copies of one cell count
#'   once.  A layer of one distinct cell gives one ratio and no spread, so its
#'   distortion is measured on a 6 by 6 grid over the cell's bounding box.
#'   Projected cells are first checked with `sf::st_is_valid()` in their
#'   CRS, every non-empty cell before any summary: a cell GEOS calls invalid
#'   (edges that cross, a spike, parts that overlap) has no planar area to
#'   vouch for, and the request is refused with an error that names the
#'   first such cells and GEOS's reason and says to repair the cells with
#'   `sf::st_make_valid()`.  A valid projected cell with two vertices that
#'   fall on one point once in degrees, where s2 refuses it, is measured
#'   there without the repeat, which moves no area; when s2 refuses a valid
#'   cell it measures for another reason, the spread cannot be measured and
#'   the request is refused.
#'   [ensure_projected()] with
#'   `purpose = "area"` chooses an equal-area CRS for lon/lat input; build the
#'   cells in it.  Lon/lat cells are measured geodesically instead: their
#'   `cell_area` is `sf::st_area()`'s area in square metres (on the sphere
#'   with s2, sf's default), not a planar area in squared degrees, and on
#'   the cells s2 measures the distortion check passes by construction.
#'   Every cell's area is taken before any summary, and when s2 refuses a
#'   cell that is not valid on the sphere, such as one with two corners it
#'   reads as one point (cells returned in lon/lat can hold them; see `crs`
#'   in [create_voronoi_polygons()]), the request stops with an error of
#'   class `"spatialkit_cells_degenerate"` that gives the rows of `cells_sf`
#'   s2 refuses (`rows` on the condition) and says to transform the cells to
#'   a projected CRS, such as the one they were built in or the equal-area
#'   one `ensure_projected(purpose = "area")` picks.  With s2 switched off,
#'   sf needs the lwgeom package for that area and the request is refused
#'   without it.  The measured spread is attached as `attr(, "area_error")`.
#'   Cells that mix XY and XYZ (or XYM) geometries are measured, and their
#'   spread taken, without Z and M, since sf reads a layer in the dimension
#'   of its first geometry; so are cells that carry M throughout, which GEOS
#'   does not read.  The cells are returned with their own geometry.  In a
#'   projected CRS whose unit sf has no name for (a `+to_meter` such as
#'   1e6), `sf::st_area()` gives square metres; `cell_area` is the planar
#'   area in that CRS's unit squared all the same.
#' @param conf_level Optional confidence level in (0, 1), such as `0.95`.
#'   When given, every numeric response and predictor column also gets
#'   `..neff_*`, `..df_*`, `..ci_lo_*` and `..ci_hi_*` (see "Confidence
#'   intervals"). Default `NULL`: no such columns.
#' @param estimate Which value of a cell the frame reports: `"pooled"`
#'   (default) or `"mean"`. With `"mean"` the frame is the one described
#'   above, the frame of earlier versions, whose default this was. With
#'   `"pooled"` every numeric response and predictor column also
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
#'   own CRS. When the input carries the record of its cells' count (see
#'   "Value") and that record names a projected CRS, the one the profile
#'   worked in, the distances are measured in it, whatever CRS the layer is
#'   in. The distances are taken in x and y: a layer whose features, or the
#'   points taken from them, do not all have the same dimensions (some with a
#'   Z or M coordinate) is read without its Z and M. So is a layer whose
#'   geometries carry M (XYM or XYZM), which GEOS does not read. The value is matched as [match.arg()] matches a choice, so an
#'   abbreviation (`"m"`) is accepted and `NULL` is the default. See "Pooled
#'   values" for what the pooled value estimates and what its numbers rest
#'   on, "Value" for the columns and the attributes, and "Conditions" for
#'   the warnings.
#' @param obs_level Optional level in (0, 1), such as `0.90`, of an interval
#'   for one further observation in a cell, located like the data. It needs
#'   `estimate = "pooled"`: the interval is the pooled value plus the
#'   `(1 - obs_level) / 2` and `(1 + obs_level) / 2` quantiles of the
#'   leave-one-out residuals of the pooled values, so its two offsets are the
#'   same in every cell. It is a statement about observations on average,
#'   not about each cell. It is a separate argument from `conf_level`, whose
#'   interval is for a mean and is not a range for an observation. Any level
#'   strictly between 0 and 1 is computed the same way; only 0.90 has been
#'   checked (see "Intervals: which one answers which question"). Default
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
#'   and `FALSE` when none was: a refused `deff`, or `"kish"` ICCs of 0 or
#'   `NA` for every column (the standard errors are then the uncorrected
#'   ones, or `NA` where an ICC cannot be estimated). It is `NA` on a
#'   `cells_sf` row no point fell in. `deff_method` names the correction
#'   the standard errors carry, since their column names do not change with
#'   it: `"variogram"`, `"kish"` or `"fixed"` when the one requested was
#'   applied, `"none"` when none was, and `"kish (variogram fallback)"` or
#'   `"none (variogram fallback)"` when a `"variogram"` request could not be
#'   met and Kish stood in for it (the second when Kish's ICCs were 0 or
#'   `NA` too).
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
#'   of the second carries `requested` and `used`. Two more classes say what
#'   a design effect could not do: `"spatialkit_icc_unavailable"` names the
#'   columns whose Kish ICC cannot be estimated, and
#'   `"spatialkit_deff_no_structure"` the predictor columns given a design
#'   effect of 1 by a variogram that found no spatial structure.
#'
#'   When a correction was actually applied, an attribute `"deff_applied"` is
#'   attached recording it: `method` plus `icc_resp`/`icc_pred` and `deff`
#'   for `"kish"` (`icc_pred` is a vector named by predictor column, one ICC
#'   each; `deff` is the primary variable's per-cell design effect,
#'   all 1 when only a predictor's ICC was positive, except that it is `NA`
#'   in each cell of two or more of the primary variable's values when that
#'   variable's ICC cannot be estimated, and its entry in `icc_resp` or
#'   `icc_pred` is then `NA` too; `requested` and
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
#'   whose ICCs are all 0 or `NA`. A `deff = "kish"` request, and a `"variogram"`
#'   request that fell back to it, always
#'   records the ICCs it estimated on an attribute `"icc"` (`resp`, and
#'   `pred`, a vector named by predictor column; `NA` for a variable type
#'   with no numeric column and for a column whose ICC cannot be estimated;
#'   `reason`, a character vector named by column, is then added and says
#'   why), whether or not
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
#'   error, and the plain columns are returned. "Conditions" lists that
#'   class and the other warnings of the pooled columns, with what each
#'   condition carries.
#'
#'   When `assigned_points_sf` carries `attr(, "count_rule")`, or `cells_sf`
#'   does (the assigned layer's is read first), the result carries that
#'   record with one more element, `check`, whatever `estimate` is. The
#'   record is the one [select_resolution()] wrote: `rule`, `estimate` (the
#'   cell value the count was chosen for), `cells`, `obs_per_cell`,
#'   `response`, `scored` (what the profile scored at the count: the rows,
#'   the cells that hold one, the leave-one-out error among them, and the
#'   centres the rows were given to),
#'   `partition`, `crs` and `scored_cells`. It has class
#'   `"spatialkit_count_rule"` and a print method. `check$estimate` is
#'   `"match"` or `"mismatch"`: whether the summary reports the estimate the
#'   count was chosen for (`NA` for a count chosen for neither). A mismatch
#'   is a warning of class `"spatialkit_estimate_mismatch"` (fields
#'   `chosen_for`, `reported`, `rule`, `cells`), and the summary is returned
#'   as asked. `check$cells` is `"match"`, `"differ"` or `"not checked:
#'   <reason>"`: the response is fitted again on the rows and cells given,
#'   with the estimate the profile scored, and its scored rows, its occupied
#'   cells and its leave-one-out error are compared with the record's, the
#'   error to a relative tolerance of 1e-8. The first that differs is a
#'   warning of class `"spatialkit_cells_differ"` (fields `what`: `"n"`,
#'   `"cells"` or `"loo"`; `scored`; `found`): rows were dropped or moved
#'   between cells, or the layer is another one. The same class with `what =
#'   "record"` (fields `field`, `scored`, `found`) is raised when the
#'   assigned layer and `cells_sf` carry different records, so the rows were
#'   assigned to other cells than the `cells_sf` given: the cells of another
#'   selection, or the scored cells of a selection beside other cells built
#'   from its count, such as a lattice sized by it. `check$cells` is then
#'   `"differ"`, whatever the comparison of the rows finds. The comparison
#'   is not made, and `check$cells` says why, when the response is another
#'   column or was not summarised, the profile read the selection half
#'   (`select_on = "split"`), the cells are not the ones the profile scored
#'   (`scored_cells` is `FALSE` or `NA`), the profile scored no estimate, the
#'   input has no coordinates, the features are not points (see
#'   [assign_features_to_polygons()]), or the fit was stopped by the limit on
#'   an exact search without the 'FNN' package. Cells and rows without a
#'   record, such as polygons the caller brings, are not checked and raise
#'   neither warning.
#'
#'   A row that is the same distance from two centres lies on the boundary
#'   between their cells, which happens on gridded locations.
#'   [resolution_profile()] scores such a row in the cell of the lower-numbered
#'   centre (a distance counts as the same as the least when its square is no
#'   more than the square of the least distance lengthened by a margin, plus
#'   one part in ten million of the least squared distance; the margin is one
#'   part in 1e10 of the largest coordinate of the centres, at most a
#'   millionth of the diagonal of their bounding box and at least four times
#'   the spacing of doubles at that coordinate), and [assign_features_to_polygons()], on cells
#'   built from the centres the profile scored, gives a point on that boundary
#'   the same cell, measured in the CRS the record names, so the comparison
#'   reads `"match"` on gridded locations too. It does so on cells that carry
#'   the record, for rows that have a CRS exactly when the record names one.
#'   Rows without a CRS beside a record that names one, rows with a CRS beside
#'   a record that names none, rows assigned to the cells after they lost the
#'   record (written to a file and read back, say) and summarised with
#'   `cells_sf` set to the cells that carry it, and rows given their cells by
#'   another join are placed by the polygons alone, so a row on a boundary can
#'   sit in either cell.
#'   [build_tessellation()] draws those cells in the CRS the record names,
#'   whatever CRS the seeds arrive in, unless its `crs` names another
#'   projected CRS. Cells built in another projected CRS are the
#'   nearest-centre cells of that CRS, and cells transformed to another
#'   projected CRS are joined there, so a row near a boundary, on it or not,
#'   can sit in another cell than the profile put it in. So can a row joined
#'   in longitude and latitude, as rows and cells in degrees are when the
#'   cells do not carry the record or carry one in degrees (the Voronoi
#'   cells of the scored centres that [build_tessellation()] returns in a
#'   geographic `crs` are joined in the record's CRS while they carry it), a
#'   row in cells that carry the record but are not the Voronoi cells of its
#'   centres (seeds moved, or cells edited, after the selection), and a row
#'   that reached the record's CRS by another transformation than the
#'   profile's, which can settle a tie the other way. Cells returned in a
#'   geographic `crs` are moved back to the record's CRS for that join, and
#'   for a layer wider than a UTM zone the round trip can move a corner by
#'   more than the margin of [assign_features_to_polygons()], so a row on an
#'   edge there can sit in another cell too. The comparison then reports
#'   `"differ"`, with `what = "loo"` when it finds the same rows in the same
#'   number of cells with another leave-one-out error, and the warning names
#'   these cases. The summary's own numbers describe the cells as assigned.
#'
#'   The comparison uses the fit of the pooled columns when that is a fit of
#'   the estimate the profile scored, and fits the response once more
#'   otherwise. Without the 'FNN' package that fit compares every pair of
#'   cell centroids (pooled values) or, for each row alone in its cell,
#'   every other location (plain means), and above
#'   `options(spatialkit.exact_pairs_max = )` pairs it is not made
#'   (`"not checked: work limit"`). A summary is not compared when neither
#'   layer carries the record: `attr(x, "count_rule") <- NULL` on the
#'   assigned layer and on `cells_sf` removes it.
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
#' # Nine cells are too few for pooled values (each would be pulled toward
#' # all the other cells), so these summaries report plain means alone:
#' # estimate = "mean".
#' # IID standard errors (default) vs Kish design-effect adjustment.  The
#' # correction is large here, and meant to be: about 22 points per cell that
#' # all share the cell's part of the trend carry far fewer than 22
#' # independent pieces of information about the cell mean.
#' naive <- summarize_by_cell(assigned, response_var = "val", estimate = "mean")
#' kish  <- summarize_by_cell(assigned, response_var = "val", deff = "kish",
#'                            estimate = "mean")
#' data.frame(n = naive$n,
#'            se_naive = naive$..se_resp_val,
#'            se_kish  = kish$..se_resp_val)
#' attr(kish, "deff_applied")   # method, icc_resp, icc_pred, per-cell deff
#'
#' # A 95% interval for each cell mean as an estimate of the grand mean, on
#' # the Kish-corrected standard error and n - 1 degrees of freedom.
#' ci <- summarize_by_cell(assigned, response_var = "val", deff = "kish",
#'                         conf_level = 0.95, estimate = "mean")
#' ci[, c("poly_id", "n", "resp_mean_val", "..neff_resp_val", "..df_resp_val",
#'        "..ci_lo_resp_val", "..ci_hi_resp_val")]
#'
#' # Pooled values, and a range for one observation: a smooth field plus
#' # noise at 400 points, on an 8 x 8 grid, so that a cell holds about six.
#' set.seed(2)
#' east  <- runif(400, 0, 8000)
#' north <- runif(400, 0, 8000)
#' pts2 <- st_as_sf(
#'   data.frame(x = 5e5 + east, y = 5e6 + north,
#'              val = 10 + 3 * sin(east / 1500) + 3 * cos(north / 1200) +
#'                rnorm(400)),
#'   coords = c("x", "y"), crs = 32632
#' )
#' bnd2 <- st_sf(geometry = st_sfc(st_polygon(list(rbind(
#'   c(5e5, 5e6), c(5e5 + 8000, 5e6), c(5e5 + 8000, 5e6 + 8000),
#'   c(5e5, 5e6 + 8000), c(5e5, 5e6)
#' ))), crs = 32632))
#' grid2 <- create_grid_polygons(bnd2, target_cells = 64, type = "square")
#' assigned2 <- assign_features_to_polygons(pts2, grid2)
#' # estimate = "pooled" is the default.
#' pooled <- summarize_by_cell(assigned2, response_var = "val", obs_level = 0.90)
#' # The plain mean, the pooled value and the weight on the cell's own mean.
#' head(pooled[, c("poly_id", "n", "resp_mean_val", "..pooled_resp_val",
#'                 "..pooled_w_resp_val")])
#' # The range for one further observation in a cell.
#' head(pooled[, c("poly_id", "..obs_lo_resp_val", "..obs_hi_resp_val")])
#' attr(pooled, "pooled")   # what was fitted, and what the interval covers
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
                              estimate       = c("pooled", "mean"),
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
  # result carry the per-row `deff_applied` column, so a call at the default
  # deff = 1 gains no column for it.
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
    # Measured on the cells in XY when they mix XY and XYZ (or XYM)
    # geometries: GEOS and s2 read such a layer in the dimension of its first
    # geometry, and the area of an XYZ cell among XY ones came out wrong or
    # could not be taken.  So are cells that carry M throughout, which GEOS
    # does not read: the validity check stopped on them with GEOS's message.
    cells_xy <- .xy_if_mixed_or_m(cells_sf)
    if (isTRUE(sf::st_is_longlat(cells_sf))) {
      # Lon/lat cells: their cell_area below is s2's area of every cell, and
      # the check that follows takes s2's area for the planar one too.  s2
      # refuses a cell that is not valid on the sphere, such as one with two
      # corners it reads as one point, which cells returned in lon/lat hold
      # on gridded locations.  The check measures up to 200 of the cells: a
      # refused one among them made the request refused as unmeasurable, and
      # one outside them stopped the request after the summary with s2's
      # bare message.  Every cell is measured here, before any summary.  With
      # s2 off the areas are lwgeom's, and the check below refuses the
      # request without it, as the help says.
      if (isTRUE(sf::sf_use_s2()))
        tryCatch(invisible(sf::st_area(sf::st_geometry(cells_xy))), error = function(e)
          .stop_cells_degenerate(
            e, cells_xy, "summarize_by_cell", "`cells_sf`",
            "the measurement of the cells' areas (`area = TRUE`)",
            "call summarize_by_cell() again; the cells' areas are then planar areas in that CRS",
            pick = "the equal-area one ensure_projected(purpose = \"area\") picks for it"))
    } else {
      # A cell GEOS calls invalid has no planar area to vouch for: a ring
      # whose edges cross gets the signed sum of its lobes (a bowtie of
      # 4.25e7 m^2 got 2e7), and parts that overlap are counted twice.  The
      # distortion check measures up to 200 cells, so such a cell among them
      # made the request refused and one outside them got a density over
      # that area.  Every non-empty cell is checked, so the outcome does not
      # turn on the sample.  Lon/lat cells are s2's to judge, above.
      .stop_cells_invalid(cells_xy, cells_sf, id_col, assigned_points_sf)
    }
    # One distinct cell gives one ratio and no spread: it is measured on the
    # grid over its bounding box.  A single cell was refused as unmeasured,
    # though it has both areas.
    area_err <- .crs_area_error(cells_xy, one_by_grid = TRUE)
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
  # The record of how the cells' count was set, when the input carries one:
  # on the assigned rows first, on the cells second.  It is read from the
  # cells here, before that layer is cut down to its ID column.
  rule_rows  <- .count_rule_get(assigned_points_sf)
  rule_cells <- if (is.null(cells_sf)) NULL else .count_rule_get(cells_sf)
  count_rule <- if (!is.null(rule_rows)) rule_rows else rule_cells
  # The CRS the profile worked in, in which the pooled values' distances are
  # then measured too.
  count_crs <- if (!is.null(count_rule) && inherits(count_rule$crs, "crs") &&
                   !is.na(count_rule$crs) && !isTRUE(sf::st_is_longlat(count_rule$crs)))
    count_rule$crs else NULL
  df <- sf::st_drop_geometry(assigned_points_sf)

  # --- locate ID column ---
  id_candidates <- .sbc_id_candidates(id_col)
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
  # Returns a single rho in [0, 1], 0 for an estimate at or below 0, or NA
  # with a `reason` attribute when the data cannot give an estimate at all
  # (the values in one cell or in none, too few cells with two or more
  # values, or a column that does not vary).  The two used to be one 0, so
  # a column whose ICC was unidentified got the uncorrected standard error
  # under the label of an ICC of 0.
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
    na_icc <- function(why) structure(NA_real_, reason = why)
    vals <- numeric(0)
    grps <- character(0)
    n_const <- 0L
    n_all <- 0L  # non-missing values, with a cell id or without one
    n_id  <- 0L  # those of them with a cell id
    for (v in var_cols) {
      # Whether the column varies is read from all its values, those of rows
      # with no cell id included: they form a summary row of their own, and
      # its standard error is 0 whatever the correlation only when they are
      # equal too.  Read from the rows with a cell id alone, a column equal
      # on them and varying on the others was taken for one that does not
      # vary, and its no-id row kept the uncorrected standard error.
      # A standard deviation below .Machine$double.eps is not by itself no
      # spread: a column in small units has one, and a column near 1e200 has
      # one that overflows.  One above it is not by itself spread either: a
      # column of one value up to rounding has a standard deviation of a unit
      # in the last place of its values, above eps for values above about 4.
      # So the column is read again by the pool's test (.constant_column())
      # where its standard deviation is below the floor or not finite, or
      # its sum of squares shows no spread on the scale of its values
      # (.no_spread()), and it does not vary only where neither its raw
      # values nor its values centred and scaled show spread beyond their
      # rounding.  A column with spread fails .no_spread() and is read as
      # before, to the last bit.
      x_all <- df[[v]][!is.na(df[[v]])]
      n_all <- n_all + length(x_all)
      s_all <- sd(x_all)
      if (!is.na(s_all) &&
          (!(s_all >= .Machine$double.eps && is.finite(s_all)) ||
           .no_spread(sum((x_all - mean(x_all))^2), x_all)) &&
          .constant_column(x_all)) {
        # No spread: the column adds no variance information.
        n_const <- n_const + 1L
        next
      }
      # A row with no cell id belongs to no group: keep_unassigned = TRUE
      # output carries NA ids, and an NA level makes group_means[grps] NA
      # downstream ("missing value where TRUE/FALSE needed").
      ok <- !is.na(df[[v]]) & !is.na(df[[id_col]])
      n_id <- n_id + sum(ok)
      raw <- df[[v]][ok]
      s <- sd(raw)
      if (is.na(s)) next
      # The same reading on the rows with a cell id.  A column with spread
      # there is z-scored from .spread_scaled(), which divides by powers of
      # two, since its raw squares can underflow or overflow, or its raw
      # values lose their spread to rounding; the ICC does not change under
      # a shift and a scale.  Elsewhere the raw values are z-scored as they
      # stand.
      if (!(s >= .Machine$double.eps && is.finite(s)) ||
          .no_spread(sum((raw - mean(raw))^2), raw)) {
        if (.constant_column(raw)) next
        raw <- .spread_scaled(raw)$y
        s <- sd(raw)
      }
      vals <- c(vals, (raw - mean(raw)) / s)
      grps <- c(grps, paste(v, as.character(df[[id_col]][ok]), sep = "\r"))
    }
    grps <- factor(grps)
    k <- nlevels(grps)
    N <- length(vals)
    if (k == 0L)
      return(na_icc(
        if (n_const > 0L) .icc_reason_constant
        else if (n_all > 0L && n_id == 0L) "its values fall in no cell"
        else if (n_id >= 2L) "its values differ only on rows that fall in no cell"
        else if (n_all >= 4L) "only one of its values falls in a cell"
        else "it has fewer than four values"))
    if (k < 2L) return(na_icc("its values fall in one cell"))
    if (N < 4L) return(na_icc("it has fewer than four values"))

    ni <- tabulate(grps)
    ni <- ni[ni > 0L]
    # An ICC needs at least 2 cells with 2+ observations, and the k >= 2 /
    # N >= 4 test above does not deliver it.  All-singleton groups give SSW =
    # 0, MSW = 0 and therefore an ICC of exactly 1 -- a design effect of n from
    # data carrying no within-cell information at all.  A single group with
    # two members gives an ICC from one within-group degree of freedom, which
    # is then applied to every cell.  Both are noise dressed as a measurement,
    # so the ICC is NA: not estimated, which is not an estimate of 0.
    if (sum(ni >= 2L) < 2L || (N - k) < 2L)
      return(na_icc("fewer than two cells hold two or more of its values"))
    grand_mean <- mean(vals)
    group_means <- tapply(vals, grps, mean)

    # Between-group and within-group sum of squares
    SSB <- sum(ni * (group_means - grand_mean)^2)
    SSW <- sum((vals - group_means[grps])^2)

    MSB <- SSB / (k - 1)
    MSW <- if ((N - k) > 0) SSW / (N - k) else 0

    # Unbalanced n0 (Donner 1986)
    n0 <- (N - sum(ni^2) / N) / (k - 1)
    if (n0 < 1) return(na_icc("its cell sizes give an n0 below 1"))

    rho <- (MSB - MSW) / (MSB + (n0 - 1) * MSW)
    # Not reached after the guards above; kept so that a 0/0 can never pass
    # for an estimate.
    if (!is.finite(rho)) return(na_icc("the variance decomposition is undefined"))
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
  # An ICC that cannot be estimated is NA, and `icc_reason` (a character
  # vector named by column) says why for each such column.
  resp_rho <- NULL
  pred_rho <- NULL
  icc_reason <- character(0)
  .kish_iccs <- function() {
    resp <- NULL
    pred <- NULL
    reason <- character(0)
    # vapply() drops the `reason` attribute, so it is collected here.
    one <- function(v) {
      r <- .estimate_icc(df, id_col, v)
      if (is.na(r)) reason[[v]] <<- attr(r, "reason")
      as.numeric(r)
    }
    if (length(resp_num) > 0L) {
      resp <- one(resp_num)
      if (is.na(resp)) {
        .msg(sprintf("summarize_by_cell(): the response ICC cannot be estimated (%s).",
                     reason[[resp_num]]))
      } else if (resp > 0) {
        .msg(sprintf(
          "summarize_by_cell(): Kish deff using response ICC = %.4f", resp
        ))
      } else {
        .msg("summarize_by_cell(): response ICC estimate <= 0; Kish deff not applied for response.")
      }
    }
    if (length(pred_num) > 0L) {
      # Named by column; the SE closures look their own column up.
      pred <- vapply(pred_num, one, numeric(1))
      names(pred) <- pred_num
      for (v in pred_num) {
        if (is.na(pred[[v]])) {
          .msg(sprintf(paste0("summarize_by_cell(): the ICC of predictor '%s' ",
                              "cannot be estimated (%s)."), v, reason[[v]]))
        } else if (pred[[v]] > 0) {
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
         pred = if (is.null(pred)) 0 else pred,
         reason = reason)
  }
  # The columns whose ICC is NA because they do not vary: their standard
  # errors are 0 whatever the correlation, so they keep a design effect of 1.
  icc_const_cols <- function() names(icc_reason)[icc_reason == .icc_reason_constant]
  # The columns whose ICC cannot be estimated and whose numbers that leaves
  # unknown: not a column that does not vary, and only one with a cell of
  # two or more values (a cell of one value has no standard error to
  # withhold).  The group of rows with no cell id is a summary row too.
  icc_withheld <- function() {
    v <- setdiff(names(icc_reason), icc_const_cols())
    v[vapply(v, function(col)
      any(tapply(!is.na(df[[col]]), factor(df[[id_col]], exclude = NULL), sum) >= 2L),
      logical(1))]
  }
  if (use_kish) {
    icc_est  <- .kish_iccs()
    resp_rho <- icc_est$resp
    pred_rho <- icc_est$pred
    icc_reason <- icc_est$reason
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
    # .cell_sd(), not sd(): the raw squares overflow for values near 1e200
    # and underflow for a spread near 1e-170, which gave standard errors of
    # NA or 0 in every cell, with no warning.
    fns[[paste0("..sd_", prefix)]] <- function(x) {
      if (sum(!is.na(x)) > 1L) .cell_sd(x) else NA_real_
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
      .const <- icc_const_cols()
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
        s <- .cell_sd(x)
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
          # An ICC that cannot be estimated leaves the design effect unknown,
          # so the SE, effective n and interval are NA (the df of s^2 does
          # not depend on rho and stays n - 1).  A column that does not vary
          # has an SE of 0 whatever the correlation, so it keeps 1.
          deff_i <- if (!is.na(rho_i)) max(1, 1 + (n_valid - 1) * rho_i)
                    else if (dplyr::cur_column() %in% .const) 1
                    else NA_real_
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
    structure_from <- NULL  # whose verdict it is, for the predictors' warning
    if (no_structure(sac)) {
      structure_note <- as.character(attr(sac, "rejected_reason"))[1L]
      structure_from <- "of the supplied `sac`"
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
          structure_from <- "estimated from `response_var`"
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
      # That verdict is about the variable the variogram was fitted to, and
      # the design effect of 1 it gives is the answer for that variable.  A
      # predictor summarised beside it gets the same 1 whatever structure it
      # has of its own: beside a white-noise response, a predictor with a
      # range of its own kept the uncorrected standard errors without a word.
      # The predictors are named; nothing is changed.  A pure nugget the
      # caller built, with no verdict behind it, is applied as given.
      under <- setdiff(pred_num, resp_num)
      if (!is.null(structure_note) && length(under)) {
        several <- length(under) > 1L
        named <- sprintf("the predictor column%s %s", if (several) "s" else "",
                         .cols_text(under))
        .warn_classed("spatialkit_deff_no_structure", paste0(
          "summarize_by_cell(): the variogram %s finds no spatial structure (%s), ",
          "so %s. That verdict is about the variable the variogram was fitted ",
          "to: the standard errors%s of %s are the uncorrected ones, whatever ",
          "spatial structure %s of %s own. deff = \"kish\" estimates an ICC for ",
          "each column on its own."),
          structure_from, structure_note,
          if (length(resp_num)) sprintf("every column gets a design effect of 1, %s included", named)
          else sprintf("%s get%s a design effect of 1", named, if (several) "" else "s"),
          if (is.null(conf_level)) "" else ", effective n and intervals",
          if (several) "those columns" else "that column",
          if (several) "they have" else "it has", if (several) "their" else "its",
          data = list(columns = under, reason = structure_note))
      }
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
      # uncorrected SE is returned only when Kish's own ICC is 0 (NA when it
      # cannot be estimated), and either way the rows say which.
      use_vgm  <- FALSE
      use_kish <- TRUE
      deff     <- 1
      vgm_fallback <- why
      icc_est  <- .kish_iccs()
      resp_rho <- icc_est$resp
      pred_rho <- icc_est$pred
      icc_reason <- icc_est$reason
      icc_txt <- paste(c(
        if (length(resp_num)) sprintf("response '%s' %.3f", resp_num[[1L]], resp_rho),
        if (length(pred_num)) sprintf("predictor '%s' %.3f", pred_num, pred_rho)),
        collapse = ", ")
      # The columns whose standard errors are NA for want of an ICC.  The
      # texts of the two cases with none are as they were.
      withheld <- icc_withheld()
      withheld_txt <- if (length(withheld))
        sprintf("NA for %s, whose ICC cannot be estimated", .cols_text(withheld))
      else ""
      if (max(resp_rho, pred_rho, 0, na.rm = TRUE) > 0)
        .warn_deff_fallback(paste0(
          "summarize_by_cell(): %s Falling back to deff = \"kish\" (ICC: %s): ",
          "the standard errors carry Kish's correction in place of the ",
          "variogram's%s, and `deff_method` reads \"kish (variogram fallback)\" ",
          "on every row."),
          why, icc_txt, if (nzchar(withheld_txt)) sprintf(" (%s)", withheld_txt) else "",
          data = list(requested = "variogram", used = "kish"))
      else if (anyNA(c(resp_rho, pred_rho))) {
        # "The uncorrected ones" is said only of columns that have them: one
        # whose ICC is 0, or one that does not vary.  When every other column
        # is withheld, the text names the NA alone.
        kept <- setdiff(unique(c(resp_num, pred_num)),
                        setdiff(names(icc_reason), icc_const_cols()))
        .warn_deff_fallback(paste0(
          "summarize_by_cell(): %s Falling back to deff = \"kish\", whose ICC ",
          "is 0 or cannot be estimated (%s), so no column carries a correction: ",
          "the standard errors are %s; `deff_method` reads ",
          "\"none (variogram fallback)\" on every row."),
          why, icc_txt,
          if (!nzchar(withheld_txt)) "the uncorrected ones"
          else if (length(kept)) paste0("the uncorrected ones, and ", withheld_txt)
          else withheld_txt,
          data = list(requested = "variogram", used = "none"))
      }
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
      # A layer that mixes XY and XYZ (or XYM) geometries is read in x and y
      # (.xy_if_mixed()), before the coerce and the move into the variogram's
      # CRS and again after the coerce, whose points can mix dimensions: sf
      # reads such a layer in the dimension of its first geometry, so the
      # move stopped with "OGR error" or st_coordinates() below read other
      # coordinates.  So is a layer whose geometries carry M, which GEOS
      # does not read (.xy_if_mixed_or_m()): the coerce stopped on its
      # polygons.  These points are not returned, so nothing is said.  A
      # polygon s2 cannot take a centroid of in degrees is refused in this
      # function's words, naming `assigned_points_sf` (.pointize_auto()).
      pts_for_deff <- .xy_if_mixed_or_m(assigned_points_sf)
      # A row with a coordinate that is not finite or lies beyond 1e150 has
      # no location.  It is read as given, before the reduction to points
      # and the move, and emptied (.empty_unlocated()), as .pooled_xy()
      # reads the pooled path's rows: moved into the variogram's CRS, a
      # northing of 1e200 in EPSG:32616 came back finite and was correlated,
      # and in that CRS it was correlated at 1e200, with nothing said.
      pts_for_deff <- .empty_unlocated(pts_for_deff)
      if (!all(sf::st_geometry_type(pts_for_deff, by_geometry = TRUE) == "POINT"))
        pts_for_deff <- .xy_if_mixed(.pointize_auto(pts_for_deff, "summarize_by_cell",
                                                    "`assigned_points_sf`"))

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
      # The bound is read again after the move (.values_at_xy()): NA for an
      # empty point and for a coordinate that is not finite or lies beyond
      # 1e150.
      coords_mat <- .values_at_xy(pts_for_deff)
      # A point without a location has nothing to correlate.  It still counts
      # towards its cell's values, and .cor_stats_from_coords() leaves it out
      # of the correlation; it used to stop the call with "missing value
      # where TRUE/FALSE needed", while deff = "kish" took the same layer.
      unlocated <- is.na(coords_mat[, 1]) & !is.na(df[[id_col]])
      if (any(unlocated))
        .warn_and_log(paste0("summarize_by_cell(): %d point(s) have no location ",
                             "(%s); their values are summarised, but the ",
                             "variogram design effect of their cells is ",
                             "computed from the located points only."),
                      sum(unlocated), .UNLOCATED_WHAT)

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

  has_resp <- length(resp_num) > 0L
  has_pred <- length(pred_num) > 0L

  # --- an ICC that cannot be estimated withholds that column's numbers ---
  # Raised only where a number is withheld: a column that does not vary keeps
  # its standard errors of 0, and a column with no cell of two or more values
  # had no standard error to give.  Both still carry their reason on
  # attr(, "icc").
  if (use_kish && length(icc_reason)) {
    withheld <- icc_withheld()
    if (length(withheld)) {
      several <- length(withheld) > 1L
      primary <- if (has_resp) resp_num[[1L]] else if (has_pred) pred_num[[1L]]
      .warn_classed("spatialkit_icc_unavailable", paste0(
        "summarize_by_cell(): Kish's ICC cannot be estimated for %s, so the design ",
        "effect of %s is unknown: %s standard errors%s are NA in every cell with ",
        "two or more values%s. deff = 1 gives the uncorrected standard errors, ",
        "those of each cell's own mean."),
        paste(sprintf("'%s' (%s)", withheld, icc_reason[withheld]), collapse = ", "),
        if (several) "those columns" else "that column", if (several) "their" else "its",
        if (is.null(conf_level)) "" else ", effective n and intervals",
        if (!is.null(primary) && primary %in% withheld) ", and so is cell_weight" else "",
        data = list(columns = withheld, reason = icc_reason[withheld]))
    }
  }

  # --- single grouped summarise for all columns ---
  grouped <- df |> dplyr::group_by(.data[[id_col]])

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

  # --- a standard deviation beyond the largest double ---
  # .cell_sd() gives Inf only where the spread inside the cell is of order
  # 1e308; .se_with_deff() then returns NA, and that is said here.
  sd_over <- integer(0)
  for (sd_col in grep("^\\.\\.sd_(resp|pred)_", names(out), value = TRUE)) {
    k_inf <- sum(is.infinite(out[[sd_col]]))
    v <- sub("^\\.\\.sd_(resp|pred)_", "", sd_col)
    if (k_inf > 0L && !(v %in% names(sd_over))) sd_over[[v]] <- k_inf
  }
  if (length(sd_over)) {
    several <- length(sd_over) > 1L
    .warn_classed("spatialkit_sd_overflow", paste0(
      "summarize_by_cell(): the standard deviation of %s is beyond the largest ",
      "double (about 1.8e308), so `..sd_*` is Inf in %s and the standard errors%s ",
      "are NA there. Divide %s by a power of ten to summarise %s."),
      paste(sprintf("'%s' (%d %s)", names(sd_over), sd_over,
                    ifelse(sd_over == 1L, "cell", "cells")), collapse = ", "),
      if (sum(sd_over) == 1L) "that cell" else "those cells",
      if (is.null(conf_level)) "" else " and intervals",
      if (several) "those columns" else "the column",
      if (several) "them" else "it",
      data = list(columns = names(sd_over), cells = unname(sd_over)))
  }

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
    primary_rho <- if (has_resp) resp_rho
                   else if (has_pred) pred_rho[[1L]]
                   else 0
    if (is.na(primary_rho) && !(primary_col %in% icc_const_cols())) {
      # The primary variable's ICC cannot be estimated, so its effective n is
      # unknown in every cell of two or more values, as its ..neff_ is.  A
      # cell of one value has a design effect of 1 whatever the ICC, and a
      # cell of none a weight of 0.
      deff_per_cell <- ifelse(n_valid_primary <= 1, 1, NA_real_)
      out$cell_weight <- ifelse(n_valid_primary <= 1, n_valid_primary, NA_real_)
    } else {
      # A column that does not vary has an NA ICC and a design effect of 1.
      primary_rho <- if (is.na(primary_rho)) 0 else max(primary_rho, 0)
      deff_per_cell <- pmax(1, 1 + (n_valid_primary - 1) * primary_rho)
      if (primary_rho > 0)
        out$cell_weight <- n_valid_primary / deff_per_cell
    }
    # Recorded whenever ANY ICC corrected its column.  Keyed on the
    # primary variable alone, a clustered predictor beside an unclustered
    # response had its SEs inflated 11x on rows marked deff_applied = FALSE,
    # with no attribute.  `deff` stays the primary variable's (all 1 then,
    # or NA in its cells of two or more values when its ICC cannot be
    # estimated).
    if (max(resp_rho, pred_rho, 0, na.rm = TRUE) > 0) {
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
  # An ICC that cannot be estimated is NA, and `reason`, named by column,
  # says why; the field is there only then, so every other result keeps the
  # two fields it had.
  if (use_kish)
    attr(out, "icc") <- c(list(
      resp = if (has_resp) as.numeric(resp_rho) else NA_real_,
      # Named by predictor column: each has its own.
      pred = if (has_pred) pred_rho else NA_real_),
      if (length(icc_reason)) list(reason = icc_reason))
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
  pooled <- NULL
  if (identical(estimate, "pooled") && (has_resp || has_pred)) {
    # The design effect the t interval carries is the one recorded above as
    # applied, not the one asked for: a refused number, and Kish ICCs of 0 for
    # every column, leave the default's numbers, and the table of intervals
    # then describes them as the default's.  A column whose ICC cannot be
    # estimated is named in either case: its intervals are NA.
    deff_carried <- attr(out, "deff_applied")
    pooled <- .pooled_columns(
      assigned_points_sf, df, id_col, out[[id_col]], resp_num, pred_num,
      obs_level = obs_level, conf_level = conf_level,
      deff_label = .deff_label(
        deff_carried,
        # Under Kish each column has its own ICC, and one of 0 leaves that
        # column's standard error as it is.
        kish_zero = if (use_kish)
          unique(c(if (has_resp && isTRUE(resp_rho <= 0)) resp_num,
                   if (has_pred) pred_num[!is.na(pred_rho) & pred_rho <= 0])),
        # and one that cannot be estimated leaves it NA.
        kish_na = if (use_kish) icc_withheld() else character(0)),
      crs = count_crs)
    for (nm in names(pooled$columns)) out[[nm]] <- pooled$columns[[nm]]
    attr(out, "pooled")    <- pooled$record
    attr(out, "intervals") <- pooled$record$intervals
    .pooled_warnings(pooled$events, obs_level)
  }

  # --- the record of how the cells' count was set ---
  # Read and checked whatever `estimate` is, and attached with the result of
  # its checks: whether the summary reports the estimate the count was
  # chosen for, and whether the cells and rows it was given are the ones the
  # profile scored.  Nothing above depends on it but the CRS of the pooled
  # values' distances.
  if (!is.null(count_rule)) {
    # What the engine finds for the response on the rows summarised.  Under
    # pooled values with a count scored for them it is the fit just made;
    # otherwise one more fit of the response.
    found <- function(est) {
      if (identical(est, "pooled") && !is.null(pooled)) {
        r <- pooled$record$columns[[paste0("resp_", response_var)]]
        return(list(n = r$n, cells = r$cells, loo = r$loo, status = r$status))
      }
      loc <- if (!is.null(pooled)) pooled$loc
             else .pooled_xy(assigned_points_sf, count_crs, "summarize_by_cell",
                             "`assigned_points_sf`")
      res <- .pool_column(df[[response_var]], df[[id_col]], loc$xy, NULL, want = est,
                          use_fnn = .has_backend("FNN"))
      list(n = res$record$n, cells = res$record$cells,
           loo = if (is.null(res$fit)) NA_real_
                 else if (identical(est, "pooled")) res$fit$loo_pooled
                 else res$fit$loo_mean,
           status = res$record$status)
    }
    attr(out, "count_rule") <- .count_rule_check(
      count_rule, other = if (!is.null(rule_rows)) rule_cells,
      estimate = estimate, response_var = response_var, summarised = resp_num,
      layer = assigned_points_sf, found = found)
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
      cells_id_candidates <- .cells_id_candidates(id_col)
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
        rule_attr   <- attr(out, "count_rule")
        pre_join_id <- as.character(out[[id_col]])

        out <- dplyr::left_join(cells_slim, out, by = id_col)
        if (!is.null(icc_attr)) attr(out, "icc") <- icc_attr
        # The record of the count, with the result of its checks.  It is the
        # one read above, not whatever the cells' own attribute holds.
        attr(out, "count_rule") <- rule_attr
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
          # Taken on the cells in XY when they mix dimensions (an XYZ cell
          # among XY ones measured 1.6e14 square metres for 95 square km)
          # or carry M, which GEOS does not read and stopped on; `out` keeps
          # the cells' own geometry.  .planar_area() keeps the CRS's unit
          # where sf has no name for it: sf::st_area() gave square metres
          # there (1e6 for a square of 1000 m, which is 1e-6 in units of
          # 1000 km).
          cell_area <- suppressWarnings(
            .planar_area(.xy_if_mixed_or_m(sf::st_geometry(out))))
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

#' The features that have a coordinate without a location
#'
#' The test \code{.values_at_xy()} applies to a point, read on the x and y of
#' every vertex: a coordinate that is not finite or lies beyond
#' \code{.POOL_XY_MAX} is not a location.  A feature with one such vertex has
#' no location as a whole, since no centre, area or transform of it can be
#' taken.  Empty features are not counted: they are empty already.  The
#' features are read one by one only when some coordinate of the layer is
#' \code{NA}, not finite or beyond the bound; otherwise no feature can have
#' such a vertex.  That test reads the coordinates themselves, not
#' \code{sf::st_bbox()}: the box skips \code{NA}, so a ring with a vertex at
#' \code{(NA, y)} passed as located, and sf keeps the box it stored on the
#' column when a row is replaced with \code{[[<-}, so a row edited to
#' \code{Inf} or \code{1e200} passed too.
#'
#' With \code{lat = TRUE}, a feature of a layer in a geographic CRS with a
#' vertex at a latitude beyond 90 degrees (\code{.beyond_pole_features()})
#' has no location either: no place on the Earth is there.  The readers of a
#' caller's layer that empty its rows without a location before anything
#' else reads them pass it (\code{.empty_unlocated(lat = TRUE)}).  The
#' default reads the coordinates as numbers only, so the movers still
#' compare such a row before and after a move and say the move cannot place
#' it (\code{.rows_unplaced()}).
#'
#' @param g An sfc.
#' @param lat \code{TRUE} to read a latitude beyond 90 degrees, in a layer
#'   whose CRS is geographic, as no location.
#' @return A logical vector, one element per feature.
#' @keywords internal
#' @noRd
.unlocated_features <- function(g, lat = FALSE) {
  n   <- length(g)
  out <- logical(n)
  if (!n) return(out)
  if (isTRUE(lat) && .is_longlat(g))
    return(.unlocated_features(g) | .beyond_pole_features(g))
  # Empty polygons and lines hold no coordinate at all: v is then NULL.
  v <- unlist(unclass(g), use.names = FALSE)
  if (!length(v) || (!anyNA(v) && all(abs(v) <= .POOL_XY_MAX))) return(out)
  off <- function(v) any(!is.finite(v) | abs(v) > .POOL_XY_MAX)
  if (inherits(g, "sfc_POINT")) {
    # x and y of each point, read point by point: sf::st_coordinates() of a
    # layer of XY and XYZ points stops ("length of 'dimnames' [2] not equal
    # to array extent").
    xy <- matrix(unlist(lapply(unclass(g), `[`, 1:2), use.names = FALSE),
                 ncol = 2L, byrow = TRUE)
    empty <- is.na(xy[, 1L]) & is.na(xy[, 2L])
    return(!empty & (!is.finite(xy[, 1L]) | !is.finite(xy[, 2L]) |
                     abs(xy[, 1L]) > .POOL_XY_MAX | abs(xy[, 2L]) > .POOL_XY_MAX))
  }
  # A POINT is a vector, a LINESTRING or MULTIPOINT a matrix, and the other
  # types (collections included) lists of those.
  bad <- function(s) {
    if (is.list(s)) return(any(vapply(s, bad, logical(1))))
    if (is.matrix(s)) return(nrow(s) > 0L && off(s[, 1:2]))
    length(s) >= 2L && !all(is.na(s)) && off(s[1:2])
  }
  vapply(unclass(g), bad, logical(1))
}

#' Give rows an empty geometry of their own type and dimension
#'
#' A row without a location (\code{.unlocated_features()}) is emptied before a
#' reading it would otherwise shape: the bounding box a projection or a
#' lon/lat test is read from, a centre, a distance score.  The empty geometry
#' keeps the type and the dimensions of the feature, since sf needs one
#' dimension per layer; an empty POINT holds \code{NA} coordinates.  The
#' layer read with these rows emptied is read as the same layer with those
#' rows empty to begin with.
#'
#' sf warns "bounding box has potentially an invalid value range for longlat
#' data" whenever it reads the box of a lon/lat layer that holds a latitude
#' beyond 90; with \code{lat = TRUE} that text is not passed on while the
#' rows are emptied, and the layer returned no longer holds such a vertex,
#' so sf has no reason to say it again.
#'
#' @param x An sf or sfc object.
#' @param rows The rows to empty; \code{NULL} for the rows
#'   \code{.unlocated_features()} finds.
#' @param lat Passed to \code{.unlocated_features()}: \code{TRUE} empties
#'   the rows of a lon/lat layer with a latitude beyond 90 degrees too.
#' @return \code{x}, with those rows emptied; \code{x} itself when there are
#'   none.
#' @keywords internal
#' @noRd
.empty_unlocated <- function(x, rows = NULL, lat = FALSE) {
  if (isTRUE(lat))
    return(.muffle_lonlat_box(.empty_unlocated(
      x, rows %||% which(.unlocated_features(sf::st_geometry(x), lat = TRUE)))))
  g <- sf::st_geometry(x)
  if (is.null(rows)) rows <- which(.unlocated_features(g))
  if (!length(rows)) return(x)
  g[rows] <- sf::st_sfc(lapply(rows, function(i) {
    s   <- g[[i]]
    cls <- class(s)
    if (inherits(s, "POINT"))
      sf::st_point(rep(NA_real_, nchar(cls[1L])), dim = cls[1L])
    else if (is.matrix(s))
      structure(unclass(s)[0L, , drop = FALSE], class = cls)
    else structure(list(), class = cls)
  }), crs = sf::st_crs(g))
  if (inherits(x, "sf")) sf::st_set_geometry(x, g) else g
}

# sf's text for a lon/lat box that reaches beyond the poles or past
# longitude 360 (sf::st_is_longlat()).
.LONLAT_BOX_TEXT <- "bounding box has potentially an invalid value range for longlat data"

#' Evaluate an expression without passing on sf's lon/lat box warning
#'
#' Only that text is muffled; every other condition passes as it comes.
#' Used where the rows sf warns about are read as rows without a location
#' (\code{.empty_unlocated(lat = TRUE)}), so the caller says it in its own
#' words.
#'
#' @param expr The expression.
#' @return The value of \code{expr}.
#' @keywords internal
#' @noRd
.muffle_lonlat_box <- function(expr) {
  withCallingHandlers(expr, warning = function(w)
    if (identical(conditionMessage(w), .LONLAT_BOX_TEXT)) invokeRestart("muffleWarning"))
}

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
#' location to the caller.  The bound is read on the coordinates as given,
#' before the projection, and again after it: a move between CRSs can bring
#' a coordinate beyond the bound back finite (a northing of 1e200 in
#' EPSG:32616 lands at a finite point in EPSG:32617), and a transform of it
#' is not a place the caller gave.  The rows beyond it as given are emptied
#' before the projection, so they take no part in choosing it: the layer is
#' read as it is read with those rows empty.  The points on the surface of
#' lon/lat polygons are \code{coerce_to_points()}'s, read as every other
#' caller reads them: sf's warning about them is not passed on, and a
#' polygon written across longitude 180 gets a point on it (or the warning
#' of class \code{spatialkit_surface_point_off}).
#' A layer that mixes XY and XYZ (or XYM) geometries is read in XY, since sf
#' reads a layer in the dimension of its first geometry, and so is a layer
#' whose geometries carry M, which GEOS does not read
#' (\code{.xy_if_mixed_or_m()}); so are the points reduced from a layer whose
#' points do not all have one dimension.  No geometry is returned, so the
#' caller's keep their Z and M.
#'
#' @param layer The assigned layer.
#' @param crs \code{NULL}, or the projected CRS to measure in: the one the
#'   record of the cells' count names, in which the profile measured.  The
#'   points are then put in that CRS, and no other is chosen.
#' @param caller,arg The function and its argument that \code{layer} is,
#'   which a refusal of a polygon s2 cannot take a centroid of in degrees,
#'   or of a centroid GEOS puts outside its feature, names
#'   (\code{.pointize_auto()}); \code{NULL} raises the refusal of
#'   \code{coerce_to_points()} as it comes.
#' @return \code{NULL} when the layer has no geometry (a plain data frame);
#'   otherwise a list with \code{xy}, a two-column matrix, and \code{crs},
#'   the CRS its distances are measured in.
#' @keywords internal
#' @noRd
.pooled_xy <- function(layer, crs = NULL, caller = NULL, arg = NULL) {
  if (!inherits(layer, "sf")) return(NULL)
  if (nrow(layer) == 0L)
    return(list(xy = matrix(numeric(0), 0L, 2L), crs = sf::st_crs(layer)))
  # A layer that mixes XY and XYZ (or XYM) geometries is read in XY: sf reads
  # a layer in the dimension of its first geometry.  So is a layer whose
  # geometries carry M, which GEOS does not read (.xy_if_mixed_or_m()).  The
  # points reduced from a layer in one dimension can mix them (the point on
  # the surface of a POLYGON Z is XY beside the POINT Z rows kept), so they
  # are tested again.  Points in an sfc_GEOMETRY column are read as the same
  # points in sfc_POINT (.typed_points()): sf::st_coordinates() refuses that
  # class.
  layer <- .xy_if_mixed_or_m(.typed_points(layer))
  # A row of a layer in a geographic CRS with a vertex at a latitude beyond
  # 90 degrees has no location (.unlocated_features(lat = TRUE)), and is
  # emptied before the reduction to points: a polygon with such a vertex got
  # a point from GEOS, with sf's bare "bounding box has potentially an
  # invalid value range for longlat data", and was pooled there.
  layer <- .empty_unlocated(layer, .beyond_pole_rows(layer), lat = TRUE)
  if (!all(sf::st_geometry_type(layer, by_geometry = TRUE) == "POINT"))
    layer <- .xy_if_mixed(
      if (is.null(caller)) coerce_to_points(layer, "auto")
      else .pointize_auto(layer, caller, arg))
  given <- unname(sf::st_coordinates(layer)[, 1:2, drop = FALSE])
  far_given <- abs(given[, 1L]) > .POOL_XY_MAX | abs(given[, 2L]) > .POOL_XY_MAX
  # Those rows are emptied before the projection: the bounding box it is
  # chosen from, and the test of a layer without a CRS for lon/lat, then read
  # what they read with the rows emptied.  The empty point keeps the
  # dimensions of the feature, as sf needs one dimension per layer.
  layer <- .empty_unlocated(layer, which(far_given))
  # A row the move cannot place comes back empty and gets NA coordinates
  # below, and both callers (summarize_by_cell(), cell_values_at()) report
  # it in their own words, so ensure_projected()'s warning about it, which
  # names `x`, is not given as well.
  done  <- .answer_for_classed("spatialkit_rows_unplaced")
  layer <- tryCatch(if (is.null(crs)) ensure_projected(layer)
                    else ensure_projected(layer, target_crs = crs),
                    finally = done())
  xy <- unname(sf::st_coordinates(layer)[, 1:2, drop = FALSE])
  far <- abs(xy[, 1L]) > .POOL_XY_MAX | abs(xy[, 2L]) > .POOL_XY_MAX
  xy[which(far | far_given), ] <- NA_real_
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
#' @param kish_na Under Kish's correction, the columns whose ICC cannot be
#'   estimated and whose standard errors are therefore \code{NA}.
#' @return One string for \code{.intervals_table()}: \code{1},
#'   \code{"kish"} or \code{"variogram"} with their quotes, or the number.
#'   It reads \code{1}, alone or before a semicolon, only when nothing was
#'   applied: the table takes the words of its row from that.  The columns
#'   in \code{kish_na} are named after a semicolon, after \code{1} as after
#'   \code{"kish"}: a request that corrected no column can still leave one
#'   \code{NA}.
#' @keywords internal
#' @noRd
.deff_label <- function(applied, kish_zero = character(0), kish_na = character(0)) {
  na_txt <- if (length(kish_na))
    sprintf("; NA for %s, whose ICC cannot be estimated", .cols_text(kish_na))
  if (is.null(applied)) return(paste0("1", na_txt))
  switch(applied$method,
    kish = paste0(
      "\"kish\"",
      if (length(kish_zero))
        sprintf("; 1 for %s, whose ICC is 0", .cols_text(kish_zero)),
      na_txt),
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
#' @param crs \code{NULL}, or the CRS the distances are measured in, as
#'   \code{.pooled_xy()} takes it.
#' @return A list: \code{columns}, the new columns by name, in the order they
#'   are appended; \code{record}, the list of class \code{spatialkit_pooled};
#'   \code{events}, what \code{.pooled_warnings()} reads (with
#'   \code{response}, the response columns, which its texts tell from the
#'   predictors); and \code{loc}, the
#'   coordinates the columns were fitted on (\code{.pooled_xy()}'s list).
#' @keywords internal
#' @noRd
.pooled_columns <- function(layer, df, id_col, out_ids, resp_num, pred_num,
                            obs_level, conf_level, deff_label, crs = NULL) {
  use_fnn <- .has_backend("FNN")
  id   <- df[[id_col]]
  loc  <- .pooled_xy(layer, crs, "summarize_by_cell", "`assigned_points_sf`")
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
    skipped_cells = length(unique(id[unlocated])),
    skipped_ll    = inherits(layer, "sf") && .is_longlat(layer),
    response      = resp_num)
  list(columns = columns, record = record, events = events, loc = loc)
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
    nums <- .count_text(c(ev$work[[top]], ev$limit[[top]]))
    .warn_classed("spatialkit_pool_unavailable", paste0(
      "summarize_by_cell(): the pooled columns of %s are NA: finding each cell's ",
      "neighbours among %d cells without the FNN package would compare %s pairs, ",
      "above the limit of %s (options(spatialkit.exact_pairs_max = )). %s The plain ",
      "columns are unaffected."),
      .cols_text(w), cells[[top]], nums[1L], nums[2L], .FNN_REMEDY,
      data = list(reason = "work limit", columns = w, work = ev$work[[top]],
                  limit = ev$limit[[top]], cells = cells[w]))
  }
  w <- cols[status == "values out of range"]
  if (length(w))
    .warn_classed("spatialkit_pool_unavailable", paste0(
      "summarize_by_cell(): the pooled columns of %s are NA: %s values lie beyond what the ",
      "arithmetic of pooled values holds, which is values up to %s in absolute value whose ",
      "largest distance from their median is at least %s; beyond these bounds the squares ",
      "of many values can overflow or lose their precision. Rescaled by a power of ten to ",
      "within that, %s pooled. The plain columns are unaffected."),
      .cols_text(w), if (one(w)) "its" else "their", .bound_text(.POOL_Y_MAX),
      .bound_text(.POOL_Y_MIN), if (one(w)) "the column is" else "the columns are",
      data = list(reason = "values out of range", columns = w))

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
      "summarize_by_cell(): the values of %s do not vary inside any cell, to the ",
      "precision of the values, so every own weight is 1 and the pooled value is the ",
      "plain mean. attr(, \"pooled\") says which columns are constant."),
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
  if (length(w)) {
    # Two cells and tau2 at 0: every own weight is 0 and each cell's target
    # is the other cell's mean, so the pooled values are the plain means
    # exchanged.  The text names the columns this happened to.
    swap <- w[cells[w] == 2L & status[w] == "no signal between cells"]
    swapped <- if (!length(swap)) "" else sprintf(paste0(
      " For %s, with 2 cells and `tau2` at 0, each cell's pooled value is the other ",
      "cell's mean: the two cells' values change places."), .cols_text(swap))
    .warn_classed("spatialkit_pool_few_cells", paste0(
      "summarize_by_cell(): %s %s values in fewer than 10 cells (%s), so each cell's ",
      "target is the mean of all the other cells, not of cells near it. A cell whose ",
      "own weight is below its share of the column's rows then gets a pooled value on ",
      "the far side of the overall mean from its plain mean, and two cells can change ",
      "order.%s ",
      "The `..pooled_w_*` columns hold the weights; estimate = \"mean\" gives the plain ",
      "means."),
      .cols_text(w), if (one(w)) "has" else "have", listed(cells[w]), swapped,
      data = list(columns = w, cells = cells[w], own_weight = ev$w_min[w]))
  }
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
    # What a smooth-looking map does not show: the response, when the
    # response alone is named; a predictor is not the response.
    smooth <- if (all(worse %in% ev$response) && one(worse)) "the response varies"
              else if (one(worse)) "the column varies"
              else "the columns vary"
    text <- c(
      if (length(worse)) sprintf(paste0(
        "the pooled values of %s predict a left-out observation no better than one ",
        "mean for all rows%s, so %s no spatial structure these cells can use. The ",
        "values are returned. A map of them looks smooth because each is mostly its ",
        "neighbours' mean, not because %s smoothly."),
        .cols_text(worse), inside, if (one(worse)) "it shows" else "they show", smooth),
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
      "summarize_by_cell(): %d row(s) in %d cell(s) have no location (%s). They count ",
      "in `n` and the plain columns and take no part in the pooled columns."),
      ev$skipped_rows, ev$skipped_cells, .unlocated_what(ev$skipped_ll),
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
    # A column with tau2 = 0 is fitted, and each of its values is its target.
    # Its block is printed in full, with a line that says so in place of the
    # status: its cells can still predict a left-out row better than one
    # mean, and the loo line says whether they do.
    no_tau2 <- identical(r$status, "no signal between cells")
    if (!identical(r$status, "ok") && !no_tau2) {
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
    if (no_tau2)
      cat(label("note"),
          paste(strwrap(paste0("tau2 is 0: every own weight is 0, and each value is its ",
                               "neighbours' mean, with none of the cell's own rows in it"),
                        width = 62),
                collapse = paste0("\n", indent)),
          "\n", sep = "")
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
