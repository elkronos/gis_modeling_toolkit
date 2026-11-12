# Help topics shared by many functions: how input geometry is read, the
# condition classes the package raises, and the options it reads.  None
# documents an object, so each is a NULL topic with a name of its own.

#' How input geometry is read
#'
#' @description
#' The rules the package applies to the spatial layers it is given: how two
#' layers are brought into one CRS, how a layer without a CRS is read, how
#' mixed dimensions and rows without a location are handled, what a
#' boundary needs, and what \pkg{s2} refuses in longitude and latitude. A
#' function's own page says how it applies them: which of its arguments
#' they cover, and whether a row it leaves out is reported by a warning or
#' in the log.
#'
#' @section Coordinate reference systems:
#' A function that takes two layers, such as points and a boundary or
#' features and cells, brings them into one CRS before it reads either.
#' [harmonize_crs()] is the rule in its plain form: when both layers have a
#' CRS, one is transformed into the other's; a layer without one is read by
#' the test in "Layers without a CRS" and then reprojected or stamped, with
#' a warning.
#'
#' A function that measures distances moves a lon/lat layer into the
#' projected CRS [ensure_projected()] chooses from the layer's extent: the
#' UTM zone for a local extent, otherwise the candidate that distorts
#' distances least (its page gives the cases). A layer that already has a
#' projected CRS is used as it is. Distances, and arguments given as
#' distances, are in the units of the CRS the function works in, which its
#' page names. For a lon/lat layer the function projects, that is the CRS
#' [ensure_projected()] chose, whose unit is the metre.
#'
#' @section Layers without a CRS:
#' A layer without a CRS is read by a test of its coordinates. When their
#' bounding box lies within longitude -180 to 180 and latitude -90 to 90,
#' and they either span more than one unit on some axis or carry the
#' decimal places of degrees, they are taken as EPSG:4326, with a warning.
#' Otherwise they are left as they are. A small planar survey inside that
#' box is read as degrees, so set the CRS of such a layer with
#' `sf::st_crs(x) <- ...` before passing it.
#'
#' Given a CRS to read the layer in (`target_crs`, or the CRS of a second
#' layer), coordinates the test takes for degrees are reprojected from
#' EPSG:4326, and the others are stamped with that CRS and not moved. A
#' stamp is a relabel: check that the coordinates are in that CRS. A CRS in
#' longitude and latitude is never stamped on coordinates that cannot be in
#' it. An x outside -180 to 360, or a y beyond 90 north or south, stops the
#' call with an error that gives the bounding box of the coordinates.
#'
#' @section Mixed dimensions, Z and M:
#' \pkg{sf} reads a layer in the dimension of its first geometry, and GEOS
#' refuses M. A function that moves a layer into another CRS and returns
#' it, or reduces it to points, reads a layer that mixes XY and XYZ (or
#' XYM) geometries, as `rbind()` of an XY and an XYZ layer gives, in x and
#' y. It returns that layer in XY, with a warning of class
#' `spatialkit_zm_dropped` (its page names any exception). A layer in one
#' dimension, XYM and XYZM included, is moved with its Z and M, and so is a
#' layer whose only geometries of another dimension than the first are
#' empty and not points.
#'
#' A layer a function only reads for its x and y, such as a boundary or
#' the points and cover of a builder, is read in x and y when it mixes
#' dimensions or its geometries carry M (XYM or XYZM). [coerce_to_points()]
#' reduces a geometry that carries M in x and y as well, with the same
#' warning; its page gives the exceptions.
#'
#' @section Rows without a location:
#' A row has no location when its geometry is empty, or one of its
#' coordinates is not finite or lies beyond 1e150. Such a row takes no part
#' in the choice of a projection, and its coordinates none in the lon/lat
#' test. Moved into another CRS, it comes back empty, in its own geometry
#' type. [sf::st_transform()] alone would drop the vertex it cannot move,
#' or put a coordinate beyond 1e150 at a finite place, and return the rest
#' of the feature.
#'
#' In a geographic CRS a vertex at a latitude beyond 90 is no place on the
#' Earth. A row with one takes no part in the choice of a projection either.
#' [coerce_to_points()], [prep_model_data()] and the fits, the points of the
#' k-means seeders, the features of [assign_features_to_polygons()], the
#' pooled columns of [summarize_by_cell()] and [cell_values_at()], and
#' [ensure_stable_poly_id()] read such a row as a row without a location.
#' [coerce_to_points()] and [ensure_stable_poly_id()] say so with a warning
#' of class `spatialkit_lat_beyond_90` that names the rows.
#'
#' A row that has a location and that the target CRS cannot place in full
#' comes back empty, or without the vertices that CRS cannot place, with a
#' warning of class `spatialkit_rows_unplaced` that names the rows.
#'
#' What a function does with such a row (drops it, keeps it with `NA`, or
#' leaves it out of one step) is on its page.
#'
#' @section Boundaries:
#' A boundary is a polygonal sf or sfc layer: the study area that cells are
#' clipped to, that seeds are drawn in, or that fold blocks are cut from.
#' \itemize{
#'   \item When exactly one of the data layer and the boundary has a CRS,
#'     the other is read in it, with a warning. A boundary without a CRS
#'     given with lon/lat points is read as lon/lat when its coordinates
#'     fit longitude and latitude, and refused with an error otherwise.
#'   \item A row that is empty, or has a coordinate that is not finite or
#'     lies beyond 1e150, read as given and again after any move into
#'     another CRS, covers no ground and is left out. A boundary with no
#'     such row left is refused with an error that says so.
#'   \item A boundary in a geographic CRS with a vertex at a latitude beyond
#'     90 is refused with an error of class `spatialkit_cells_degenerate`
#'     that names its rows and says to check its CRS.
#'   \item A row the repair of the boundary leaves without a polygon (a ring
#'     with no height: every vertex on one line, or in degrees on one
#'     parallel) is left out as well, and a boundary with no polygon left is
#'     refused.
#'   \item A boundary that mixes dimensions or carries M is read in x and y.
#' }
#' [get_voronoi_seeds()] adds the rules for drawing points inside a
#' boundary, and [make_folds()] the rule for a block the boundary only
#' touches.
#'
#' @section Polygons in degrees:
#' With s2 on (`sf::sf_use_s2()`, sf's default), sf reads polygons in
#' longitude and latitude on the sphere, and s2 refuses a polygon that is
#' not valid there: one whose edges cross, or with two corners it reads as
#' one point. Cells returned in lon/lat can hold two such corners. This
#' happens where four Voronoi seeds stand the same distance from one point,
#' and the transform to degrees puts two corners that differ only by
#' rounding at one point (see `crs` in [create_voronoi_polygons()]).
#'
#' A function that joins, measures or unions such a polygon in degrees
#' stops with an error that names its rows, quotes s2's reason and says to
#' transform the layer to a projected CRS, such as the one
#' [ensure_projected()] picks for it. The class is
#' `spatialkit_cells_degenerate` for cells and boundaries, and
#' `spatialkit_features_degenerate` for features being assigned or reduced
#' to points. Some steps repair the polygons first, such as the
#' largest-overlap join of [assign_features_to_polygons()] and the union of
#' a boundary in [get_voronoi_seeds()], and stop only when s2 refuses a
#' repaired polygon.
#'
#' A hex or square grid built with a geographic `crs` is laid in degrees
#' ([create_grid_polygons()]). Voronoi cells asked for in a geographic `crs`
#' are built in a projected CRS and then transformed.
#'
#' @section Coordinates far from the origin:
#' Far from the origin of a CRS, doubles are spaced more widely: 0.125 apart
#' at 1e15. Where they are spaced more than a millionth of a layer's extent
#' apart, edges and points are held to that spacing, and a point within it
#' of an edge can fall on the other side. The functions that build cells,
#' blocks and joins still run, with a warning of class
#' `spatialkit_xy_imprecise` that gives the spacing and the extent;
#' [clip_target_for()] stops with an error of class
#' `spatialkit_xy_imprecise_error` where its buffer comes back empty or
#' with no area. Those that search for neighbours ([estimate_sac_range()],
#' [fit_gwr_model()], [residual_morans_i()]) give the same warning where a
#' distance lies past a tie by less than the rounding. Subtract a constant
#' from the coordinates to bring them nearer the origin.
#'
#' @seealso [spatialkit-conditions] for every condition class named here,
#'   [ensure_projected()] for the choice of a projected CRS, and
#'   [harmonize_crs()] for aligning two layers.
#' @family data preparation
#' @name spatialkit-geometry
NULL


#' Warnings and errors the package raises
#'
#' @description
#' Every condition class the package raises, with the functions that raise
#' it, what it means, what to do and the fields it carries. Each class can
#' be caught by name.
#'
#' @section Catching a condition by class:
#' Each condition below is an R warning or error whose class vector starts
#' with the class named, followed by `"warning"` or `"error"` and
#' `"condition"`. A handler named after the class catches it without
#' matching the message, and reads its fields with `$`:
#' \preformatted{
#' w <- tryCatch(determine_optimal_levels(pts),
#'               spatialkit_no_elbow = function(w) w)
#' w$k
#' }
#' `withCallingHandlers()` with `invokeRestart("muffleWarning")`, or
#' `suppressWarnings(expr, classes = )`, lets the call finish and drops the
#' warning. `options(warn = 2)` turns every warning into an error.
#'
#' Some warnings and errors carry no class of the package's own, such as
#' the refusal of an argument of the wrong type. Catch those with
#' `warning =` or `error =`.
#'
#' @section Warnings and log lines:
#' The package also writes a log through \pkg{logger}: a trace in a
#' temporary file, and a console echo of every line at the WARN level and
#' above. A log line is not an R condition. `tryCatch(warning = )` does not
#' catch it, and `suppressWarnings()` does not silence it. Where a help page
#' says a function raises a warning, it means an R warning; where it says
#' the function logs one, it means a log line. [spatialkit_quiet()] turns
#' the console echo off.
#'
#' The package writes a classed warning to the log as well as raising it,
#' so a warning muffled by a handler still prints its log line until the
#' echo is turned off.
#'
#' @section Condition classes:
#' Each entry gives the class, the functions that raise it, what it means
#' and what to do, and then the fields a handler reads with `$`. Each is a
#' warning unless marked as an error. An entry without fields is a
#' condition that carries only its message.
#'
#' * `spatialkit_aoa_zero_threshold`. Raised by [area_of_applicability()]
#'   with `folds`. The dissimilarity threshold is 0: held-out rows have
#'   exact copies on the training side of their fold. Use folds that keep
#'   copies together, such as `make_folds(method = "leave_location_out",
#'   group_var = )`.
#' * `spatialkit_bound_chosen`. Raised by [get_voronoi_seeds()],
#'   [voronoi_seeds_kmeans()] and [build_tessellation()], given a profile or
#'   a selection as the count. The optimum sits on a bound of the ladder, so
#'   a bound of the ladder set the count. Read `edge` in the
#'   [select_resolution()] result; [resolution_profile()] describes the
#'   bounds. Fields: `edge`, `criterion` and `n`.
#' * `spatialkit_cells_degenerate` (error; a warning from
#'   [ensure_stable_poly_id()]). Raised by [build_tessellation()],
#'   [create_voronoi_polygons()], [create_grid_polygons()],
#'   [clip_target_for()], [get_voronoi_seeds()], [voronoi_seeds_random()],
#'   [make_folds()], [assign_features_to_polygons()], [summarize_by_cell()],
#'   [cell_values_at()] and [ensure_stable_poly_id()]. \pkg{s2} refuses a
#'   cell or a boundary row in lon/lat even after repair, or a boundary in a
#'   geographic CRS has a vertex beyond latitude 90. From
#'   [ensure_stable_poly_id()]: the repair left a row with no polygon.
#'   Transform the layer to a projected CRS; for a latitude beyond 90, check
#'   its CRS. Fields: `rows`, `n` and `reasons`; the error also carries
#'   `parent`, `repaired` and `s2_forced`, and some steps add `site`,
#'   `cause` or `features_rows`.
#' * `spatialkit_cells_differ`. Raised by [summarize_by_cell()]. The rows or
#'   cells given are not those the count was scored on; `what` says which
#'   check differs. Summarise the layer the count was chosen on, assigned to
#'   the cells built from that selection. Fields: `what`, `scored` and
#'   `found`, and `field` for the count record.
#' * `spatialkit_constant_predictor`. Raised by [prep_model_data()], and so
#'   by every fitter and `cv_*()` function. A predictor takes one value in
#'   every row kept, so its coefficient cannot be told from the intercept.
#'   Drop it from `predictor_vars`.
#' * `spatialkit_constant_response`. Raised by [resolution_profile()]. The
#'   scored variable is constant to within rounding, so `cp`, `moran_i` and
#'   `moran_z` are `NA` at every level. Check the response. With
#'   `predictor_vars`, a response the predictors fit exactly leaves constant
#'   residuals. Fields: `n`.
#' * `spatialkit_convergence`. Raised by [fit_bayesian_spatial_model()],
#'   [cv_bayes()] and [cv_spatial()]. The sampler's diagnostics fail
#'   (divergent transitions, R-hat, effective sample size or E-BFMI), and
#'   `convergence_ok` is `FALSE`. Raise `adapt_delta` for divergent
#'   transitions; run more iterations for R-hat and effective sample size.
#'   Fields: `failed`, and from cross-validation also `folds` and
#'   `fold_messages`.
#' * `spatialkit_cp_noise`. Raised by [resolution_profile()]. `cp` does not
#'   take the variogram's nugget as it stands, or has no noise variance and
#'   is `NA` at every level. Read `attr(, "cp_noise")` for the variance `cp`
#'   used and its source. Fields: `source`, `value`, `level`, `cap`,
#'   `uncapped_source` and `uncapped_value`.
#' * `spatialkit_deff_fallback`. Raised by [summarize_by_cell()]. The design
#'   effect requested was not the one applied: a `"variogram"` request fell
#'   back to Kish's, or a number given was refused. Read the `deff_method`
#'   column; `deff_method == "variogram"` keeps the rows the variogram
#'   corrected. Fields, for a `"variogram"` request: `requested` and `used`.
#' * `spatialkit_deff_no_structure`. Raised by [summarize_by_cell()]. The
#'   variogram of the response found no spatial structure, so the predictor
#'   columns named get a design effect of 1 from a verdict about the
#'   response. Use `deff = "kish"` to give each predictor its own ICC.
#'   Fields: `columns` and `reason`.
#' * `spatialkit_diagnostic_criterion` (error). Raised by
#'   [select_resolution()]. `criterion = "moran_z"` was asked for without
#'   `diagnostic = TRUE`. `moran_z` is a diagnostic of the cell means, which
#'   [select_resolution()] reads only under `diagnostic = TRUE`. Pass
#'   `diagnostic = TRUE`, or choose another criterion.
#' * `spatialkit_diagnostic_ranking`. Raised by [determine_optimal_levels()]
#'   with `criterion = "morans_i"`. The candidates are ranked by Moran's z
#'   of the cell means, a diagnostic of the cells. Use `criterion =
#'   "geometric"`, or [resolution_profile()] when the response should drive
#'   the count. Fields: `k`, `criterion`, `first_computable` and `window`.
#' * `spatialkit_empty_cells`. Raised by [resolution_profile()],
#'   [get_voronoi_seeds()] and [voronoi_seeds_kmeans()] with `partition =
#'   "density"`. A centre holds no point, so fewer cells are occupied than
#'   the count. It happens at high counts, mostly from random starts.
#'   k-means++ starts (`init = "kmeans++"`) leave a cell empty less often.
#'   Fields: `level` and `empty`.
#' * `spatialkit_estimate_mismatch`. Raised by [summarize_by_cell()]. The
#'   summary reports another cell value than the count was chosen for: plain
#'   means against pooled values, or the reverse. Report the value the count
#'   was chosen for, or choose the count for the value you report. Fields:
#'   `chosen_for`, `reported`, `rule` and `cells`.
#' * `spatialkit_features_degenerate` (error). Raised by
#'   [coerce_to_points()], [prep_model_data()],
#'   [assign_features_to_polygons()] and [cell_values_at()]. \pkg{s2}
#'   refuses a polygon of the features in lon/lat, GEOS puts the centroid of
#'   an invalid feature outside it, or a temporary projection leaves a line
#'   empty. Transform the features to a projected CRS, choose another `mode`
#'   in [coerce_to_points()], or repair them with [sf::st_make_valid()].
#'   Fields: `rows`, `n`, `reasons`, `parent`, `repaired` and `s2_forced`,
#'   and some steps add `site`, `crs` or `auto_helps`.
#' * `spatialkit_few_residuals`. Raised by [summarize_by_cell()] with
#'   `obs_level`. Fewer than five residuals are expected beyond each end of
#'   the interval (at 0.90, fewer than 100 residuals), so its ends are
#'   unstable. Ask for a lower `obs_level`, which needs fewer residuals.
#'   Fields: `columns`, `needed`, `level` and `n`.
#' * `spatialkit_fold_args_ignored`. Raised by [cv_spatial()], [cv_rf()],
#'   [cv_gwr()], [cv_bayes()] and [compare_models_cv()]. `block_size` or
#'   `auto_range = TRUE` was given beside `folds`, which leaves them nothing
#'   to size. Drop them, or leave out `folds` so the function builds the
#'   blocks.
#' * `spatialkit_gp_adequacy`. Raised by [fit_bayesian_spatial_model()],
#'   [cv_bayes()] and [cv_spatial()]. The basis of the Bayesian GP is too
#'   coarse: the posterior asks for a shorter length-scale than it resolves,
#'   or the residuals are autocorrelated at short range. Follow the message:
#'   a finer basis (`gp_k`, `gp_resolution`) or, at its cap, kriging or a
#'   model of cell means. Fields, from cross-validation: `folds` and
#'   `fold_messages`.
#' * `spatialkit_gp_basis_size`. Raised by [fit_bayesian_spatial_model()],
#'   [cv_bayes()] and [cv_spatial()]. Fewer than one row per ten basis
#'   functions, so most basis weights are set by their prior. A smaller
#'   `gp_k` is cheaper and resolves less; with this few rows kriging fits
#'   the same kind of model. Fields, from cross-validation: `folds` and
#'   `fold_messages`.
#' * `spatialkit_gp_resolution`. Raised by [fit_bayesian_spatial_model()],
#'   [cv_bayes()] and [cv_spatial()]. A variogram fitted before compiling
#'   puts the response's range below what the basis resolves, and the
#'   response's autocorrelation agrees. The check after sampling tests this
#'   estimate. The message gives the basis size that would resolve the
#'   range. Fields, from cross-validation: `folds` and `fold_messages`.
#' * `spatialkit_gwr_extrapolation`. Raised by [predict.gwr_fit()],
#'   [predict_surface()], [cv_gwr()], [cv_spatial()],
#'   [cv_block_size_sweep()] and [select_features_forward()]. New locations
#'   lie beyond the reach of the widest kernel; they are predicted and
#'   counted. Where extrapolation is intended, use `suppressWarnings(expr,
#'   classes = "spatialkit_gwr_extrapolation")`. Fields: `n_far`, `n`,
#'   `reach`, `reach_text`, `nearest` and `calls`.
#' * `spatialkit_gwr_not_estimable`. Raised by [predict.gwr_fit()],
#'   [predict_surface()], [cv_gwr()], [cv_spatial()],
#'   [cv_block_size_sweep()] and [select_features_forward()]. The local
#'   regression at a new location cannot be estimated, so its prediction is
#'   `NA`. [predict.gwr_fit()] says when a window cannot be estimated.
#'   Fields: `n_na`, `n`, `n_empty`, `bw_text` and `calls`.
#' * `spatialkit_gwr_rank_deficient`. Raised by [fit_gwr_model()]. Some
#'   local windows have deficient rank, and their coefficients are `NA`; the
#'   fitted values are returned. Where many windows hold a predictor
#'   constant, fit it as one global coefficient (a mixed GWR, as
#'   `GWmodel::gwr.mixed()` fits).
#' * `spatialkit_gwr_repeated_locations`. Raised by [fit_gwr_model()] and
#'   [gwr_model_selection()]. Rows share coordinates. AICc counts them as
#'   independent, so a bandwidth left to the search (`bandwidth = NULL`) is
#'   chosen by leave-one-location-out error, and the AICc gaps between
#'   models make what a predictor adds look larger than it is. Read
#'   "Observations that share a location" in [fit_gwr_model()].
#' * `spatialkit_icc_unavailable`. Raised by [summarize_by_cell()]. A
#'   column's Kish ICC cannot be estimated, so its standard errors and
#'   intervals are `NA` in cells of two or more observations. `attr(,
#'   "icc")$reason` says why. Fields: `columns` and `reason`.
#' * `spatialkit_k_reduced`. Raised by [make_folds()], and the functions
#'   that build folds with it. Fewer folds were made than asked for. Read
#'   `k` in the result; the message gives the reason and a setting that
#'   works. Fields: `k_requested` and `k`.
#' * `spatialkit_knndm_short_of_target`. Raised by [make_folds()] with
#'   `method = "knndm"`. The held-out points sit closer to their training
#'   data than the prediction points do, so the cross-validated error tends
#'   to read low. Use fewer folds; each then holds out a larger part of the
#'   layer. Fields: `held_out_mean` and `target_mean`.
#' * `spatialkit_ladder_unused`. Raised by [resolution_profile()]. A numeric
#'   `obs_per_cell` sets the ladder, so the `n_levels`, `min_cell_n`,
#'   `sample_n` or `range_floor = TRUE` also given are not read. Pass
#'   `obs_per_cell` or those arguments, not both. Fields: `unused` and
#'   `obs_per_cell`.
#' * `spatialkit_lat_beyond_90`. Raised by [coerce_to_points()] and
#'   [ensure_stable_poly_id()]. Rows of a layer in a geographic CRS have a
#'   vertex at a latitude beyond 90, which is no place on the Earth.
#'   [coerce_to_points()] gives them an empty point; [ensure_stable_poly_id()]
#'   keys them as empty polygons and returns them as they came. Check the
#'   layer's CRS and the order of its coordinates. Fields: `rows` and `n`.
#' * `spatialkit_moran_skipped`. Raised by [resolution_profile()]. A level
#'   has more than 5,000 cells holding a scored row while \pkg{FNN} cannot
#'   be loaded or \pkg{Matrix} is not installed, so `moran_i` and `moran_z`
#'   are `NA` there. Install \pkg{Matrix}, or reinstall \pkg{FNN}. Fields:
#'   `limit`, `levels` and `missing`.
#' * `spatialkit_new_levels`. Raised by [predict.bayesian_fit()]. Rows carry
#'   a level of a categorical predictor that the model was not fitted with,
#'   so their predictions are `NA`. The warning names the predictor, the
#'   levels and the number of rows.
#' * `spatialkit_nndm_out_of_reach`. Raised by [make_folds()] with `method =
#'   "nndm"` or `"knndm"`. The prediction points sit nearer to the data than
#'   any fold can hold a point out, so the cross-validated error is
#'   pessimistic for them. No fold setting closes the gap; read the score as
#'   an upper bound on the error at those points. Fields: `deficit`,
#'   `tolerance` and `separated`.
#' * `spatialkit_no_cover`. Raised by [build_tessellation()] and
#'   [create_voronoi_polygons()]. Seeds came with neither a boundary nor a
#'   cover, so their cells cover only the seeds' own hull. Pass the study
#'   area as `boundary`, or the observations as `cover`.
#' * `spatialkit_no_default_criterion` (error). Raised by
#'   [get_voronoi_seeds()], [voronoi_seeds_kmeans()] and
#'   [build_tessellation()]. A profile passed as the count has `cp`, the
#'   criterion a profile is read by, `NA` at every level while another
#'   criterion is finite. Pass `select_resolution(profile, criterion = )`
#'   with a finite criterion, or a number.
#' * `spatialkit_no_elbow` (an error under `on_no_elbow = "error"`). Raised
#'   by [determine_optimal_levels()]. The WSS curve has no elbow, or the
#'   ladder is too short to show one, so the count it gives is set by the
#'   ladder. Choose the count on other grounds, such as
#'   [resolution_profile()] with a response. Fields: `k`, `sag` and
#'   `structured`.
#' * `spatialkit_no_scored_rows`. Raised by [resolution_profile()]. No row
#'   of the selection half has a response (`select_on = "split"`), or no row
#'   has a response and every predictor, so the columns that read those rows
#'   are `NA` at every level. `reason` says which; check where the rows with
#'   a response and every predictor lie. Fields: `reason`, `n`, `selection`,
#'   `n_response` and `column`.
#' * `spatialkit_no_structure`. Raised by [get_voronoi_seeds()],
#'   [voronoi_seeds_kmeans()], [build_tessellation()] and
#'   [summarize_by_cell()] with `estimate = "pooled"`. One cell for the
#'   whole layer does as well (`structured` is `FALSE`), so the cell values
#'   say little that one mean does not. The call goes on with the count.
#'   The message says which test failed, and whether a size set the count;
#'   see "Pooled values" in [summarize_by_cell()]. Fields: `criterion`, `n`
#'   and `one_cell`, and from [summarize_by_cell()] also `columns`,
#'   `estimate`, `loo`, `structure_p` and `beats`.
#' * `spatialkit_nonfinite`. Raised by [summarize_by_cell()]. `Inf` or
#'   `-Inf` in a response or predictor column is treated as missing in every
#'   statistic of that column; `n` still counts the row. The warning names
#'   the columns and gives the rows and cells affected. Fields: `n` and
#'   `columns`.
#' * `spatialkit_not_converged`. Raised by [resolution_profile()],
#'   [get_voronoi_seeds()] and [voronoi_seeds_kmeans()]. The start kept was
#'   cut off at the iteration limit, so a longer run could still move points
#'   between cells. Raise `iter_max` (`kmeans_iter` in
#'   [get_voronoi_seeds()]). Fields: `level`, `iter` and `iter_max`.
#' * `spatialkit_not_supported`. Raised by [resolution_profile()].
#'   `supported` is `FALSE`: the range floor exceeds the ceiling, fewer than
#'   two cells meet `min_cell_n`, or the points sit on one location.
#'   `attr(x, "bounds")` gives the floor, the ceiling and `ceiling_from`.
#'   Fields: `reason`, `n` and `min_cell_n`, and `obs_per_cell` for a ladder
#'   of sizes.
#' * `spatialkit_partition_unused`. Raised by [get_voronoi_seeds()] and
#'   [voronoi_seeds_kmeans()]. Given a profile or a selection, the call also
#'   named a `partition`, `init` or density setting other than the
#'   profile's; the scored centres are returned. Pass the count as a number
#'   (`n = sel$best`) for a fresh k-means. Fields: `kind` and `init`.
#' * `spatialkit_pool_degenerate`. Raised by [summarize_by_cell()] with
#'   `estimate = "pooled"`. The pooled values carry no pooling: one cell, no
#'   spread inside any cell, or no signal between cells. `reason` says which
#'   case applies. Fields: `reason` and `columns`, with `n`, `constant` or
#'   `sigma2` by case.
#' * `spatialkit_pool_few_cells`. Raised by [summarize_by_cell()] with
#'   `estimate = "pooled"`. The column has 2 to 9 cells, so each target is
#'   the mean of all the other cells, and a pooled value can land on the
#'   far side of the overall mean from its plain mean. With 10 cells or
#'   more each target is the mean of the observations in the 8 nearest
#'   cells. Fields: `columns`, `cells` and `own_weight`.
#' * `spatialkit_pool_few_df`. Raised by [summarize_by_cell()] with
#'   `estimate = "pooled"`. The spread inside the cells (`sigma2`), which
#'   every own weight rests on, comes from fewer than 30 values beyond one
#'   per cell. Fewer cells leave more values beyond one per cell. Fields:
#'   `columns`, `needed`, `df`, `n` and `cells`.
#' * `spatialkit_pool_rows_skipped`. Raised by [summarize_by_cell()] with
#'   `estimate = "pooled"`. Rows with a value and a cell have no location;
#'   they count in `n` and the plain columns and take no part in the pooled
#'   ones. See "Rows without a location" in [spatialkit-geometry]. Fields:
#'   `n`, `cells` and `columns`.
#' * `spatialkit_pool_unavailable`. Raised by [summarize_by_cell()] with
#'   `estimate = "pooled"` and [resolution_profile()] with an `estimate`.
#'   The pooled columns (in a profile, `loo`, `loo_se` and `own_weight`) are
#'   `NA`; `reason` says why. The plain columns are returned. For `"values
#'   out of range"`, divide the column by a power of ten; for `"work
#'   limit"`, see [spatialkit-options]. Fields: `reason` and `columns`, with
#'   `n` and `cells`, or `work`, `limit` and `cells`, by reason.
#' * `spatialkit_predict_failed`. Raised by [predict.bayesian_fit()]. The
#'   posterior draw failed, and every prediction is `NA`. The message gives
#'   the cause.
#' * `spatialkit_repeated_locations`. Raised by [summarize_by_cell()] with
#'   `obs_level` and [resolution_profile()] with an `estimate`. Rows share a
#'   location, so the interval or `loo` describes a further observation at a
#'   place already observed. Fields: `n`, the rows that share a location,
#'   and `total`, the rows read, with `columns` from [summarize_by_cell()]
#'   or `column` from [resolution_profile()].
#' * `spatialkit_replicate_var_above_nugget`. Raised by
#'   [kriging_adequacy()]. The replicates at one location vary between
#'   visits more than the variogram's nugget allows; their variance is used
#'   as it is. Refit the variogram with a nugget of at least the replicate
#'   variance. Fields: `within_var` and `nugget`.
#' * `spatialkit_response_ceiling`. Raised by [resolution_profile()]. Rows
#'   that are not scored (no response, or a predictor missing) lowered the
#'   ceiling on the count. The message names the bound in force. Fields:
#'   `n`, `n_response`, `ceiling`, `ceiling_from`, `n_sample` and
#'   `n_scored`.
#' * `spatialkit_response_imprecise`. Raised by [prep_model_data()], and so
#'   by every fitter and `cv_*()` function, and by [estimate_sac_range()],
#'   [kriging_adequacy()] and [select_features_forward()]. The response lies
#'   so far from zero that doubles there are spaced more than a millionth of
#'   its standard deviation apart, and every mean is rounded to that
#'   spacing. Subtract a constant from the response first. Fields: `spacing`
#'   and `spread`.
#' * `spatialkit_response_scale`. Raised by [fit_bayesian_spatial_model()].
#'   `prior` sets a prior on the response's scale, so a response far from
#'   zero, or with a standard deviation far from 1, is handed to \pkg{brms} as
#'   it is. See "The response's centre and scale" in
#'   [fit_bayesian_spatial_model()].
#' * `spatialkit_rows_dropped`. Raised by [prep_model_data()], and so by
#'   every fitter and `cv_*()` function, and by [gp_lengthscale_bounds()].
#'   Rows were dropped for a missing or non-finite value or for having no
#'   location. The message gives the counts and how far the dropped rows'
#'   mean response sits from the kept rows'. When that difference is large,
#'   the fit is biased towards the complete rows. Fields, from
#'   [gp_lengthscale_bounds()]: `n` and `which`.
#' * `spatialkit_rows_unplaced`. Raised by [ensure_projected()],
#'   [harmonize_crs()], [plot_tessellation_map()], and the other functions
#'   that move a layer into another CRS. The target CRS cannot place a row
#'   in full, so the row comes back empty or without those vertices. Correct
#'   or drop those rows. Fields: `rows` and `n`.
#' * `spatialkit_sd_overflow`. Raised by [summarize_by_cell()]. A cell's
#'   standard deviation is beyond the largest double, so its `..sd_*` is
#'   `Inf` and its standard errors and interval are `NA`. Divide the column
#'   by a power of ten. Fields: `columns` and `cells`.
#' * `spatialkit_size_not_reached`. Raised by [resolution_profile()],
#'   [select_resolution()] and [get_voronoi_seeds()]. A size does not give
#'   the count it asks for (fewer than two cells, or more than k-means can
#'   place), and the nearest level is used. Fields: `obs_per_cell`, `wanted`
#'   and `used`, the two counts, and `reason`.
#' * `spatialkit_surface_point_off`. Raised by [coerce_to_points()]. The
#'   point GEOS takes on a lon/lat polygon spanning more than 180 degrees of
#'   longitude does not lie on it read on the sphere, and that point is
#'   kept. Transform the layer to a projected CRS that holds the polygon,
#'   and call again. Fields: `rows`.
#' * `spatialkit_ties_unresolved`. Raised by
#'   [assign_features_to_polygons()]. The `tie_break` rule names no polygon
#'   for some features, since the polygons are equal on every key; those
#'   features get `NA`. They are kept or dropped as `keep_unassigned` says.
#'   Fields: `n` and `total`.
#' * `spatialkit_unassigned`. Raised by [assign_features_to_polygons()],
#'   [build_tessellation()] and [create_voronoi_polygons()]. Features fall
#'   inside no polygon, or points or seeds lie outside every cell. Give the
#'   cells a `boundary` or `cover` that reaches the data; `keep_unassigned =
#'   TRUE` keeps such features with `NA`. Fields: `n` and `total`.
#' * `spatialkit_variogram_unusable`. Raised by [resolution_profile()]. The
#'   variogram gives no usable range, so `reliability` is `NA`. `reason`
#'   gives the `rejected_reason` of [estimate_sac_range()]. Fields: `reason`,
#'   `converged`, `nugget` and `nugget_used`.
#' * `spatialkit_work_limit` (error). Raised by [resolution_profile()] and
#'   [cell_values_at()]. A search or a fit would make more comparisons than
#'   an option allows. The message says what to change; see
#'   [spatialkit-options]. Fields: `work`, `limit` and `what`, and for a fit
#'   also `kind`, `init`, `nstart` and `fnn`.
#' * `spatialkit_xy_imprecise`. Raised by [assign_features_to_polygons()],
#'   [build_tessellation()], [create_voronoi_polygons()],
#'   [create_grid_polygons()], [clip_target_for()], [make_folds()],
#'   [estimate_sac_range()], [fit_gwr_model()] and [residual_morans_i()].
#'   The coordinates lie so far from the origin of their CRS that doubles
#'   there are spaced more than a millionth of the extent apart, or a
#'   distance lies past a tie by less than that rounding. The call runs on
#'   the rounded coordinates. Subtract a constant from the coordinates; see
#'   "Coordinates far from the origin" in [spatialkit-geometry]. Fields:
#'   `spacing`, and `extent` or `distance`.
#' * `spatialkit_xy_imprecise_error` (error). Raised by [clip_target_for()].
#'   At that rounding the buffer comes back empty or with no area. Shift the
#'   points nearer the origin, or pass `boundary`.
#' * `spatialkit_zm_dropped`. Raised by [ensure_projected()],
#'   [harmonize_crs()], [coerce_to_points()],
#'   [assign_features_to_polygons()], [get_voronoi_seeds()],
#'   [ensure_stable_poly_id()], [predict_surface()] and
#'   [area_of_applicability()]. A layer that mixes XY and XYZ (or XYM)
#'   geometries came back in XY, or [coerce_to_points()] reduced geometries
#'   that carry M in x and y. See "Mixed dimensions, Z and M" in
#'   [spatialkit-geometry]. Fields: `rows`, the rows that carried Z or M.
#'
#' @seealso [spatialkit_quiet()] for the console log, [spatialkit-options]
#'   for the work limits, [spatialkit-geometry] for the rules behind the
#'   geometry classes, and `vignette("troubleshooting", package =
#'   "spatialkit")` for the messages met most often.
#' @family options and caches
#' @examples
#' library(sf)
#' set.seed(1)
#' # Points spread evenly have no elbow.
#' pts <- st_as_sf(data.frame(x = 5e5 + runif(200, 0, 100),
#'                            y = 5e6 + runif(200, 0, 100)),
#'                 coords = c("x", "y"), crs = 32632)
#'
#' # Catch the warning by its class and read its fields.
#' w <- tryCatch(determine_optimal_levels(pts, max_levels = 6),
#'               spatialkit_no_elbow = function(w) w)
#' class(w)
#' w$k
#'
#' # Let the call finish without the warning; the log line still prints
#' # until the console echo is off.
#' old <- spatialkit_quiet()
#' k <- suppressWarnings(determine_optimal_levels(pts, max_levels = 6),
#'                       classes = "spatialkit_no_elbow")
#' spatialkit_quiet(old)
#' k
#' @name spatialkit-conditions
NULL


#' Options that bound the work of a search or a fit
#'
#' @description
#' The package reads two options. Each bounds a step whose cost grows with
#' the square of the rows or of the cells, so that such a step stops, or is
#' skipped, before it starts. Set them with `options()`.
#'
#' @section The options:
#' \describe{
#'   \item{`spatialkit.exact_pairs_max`}{Default `2e9` pairs. The most pairs
#'     a search made without the \pkg{FNN} package may compare. Without
#'     \pkg{FNN}, a nearest-neighbour search (each cell's neighbours, each
#'     row's nearest centre, each location's nearest other location)
#'     compares every pair, and a search over the pairs of about 44,700
#'     cells or locations reaches the default. With \pkg{FNN}, which
#'     installs with spatialkit, these searches use a kd-tree and the limit
#'     is not read. The result of every search is the same with and
#'     without \pkg{FNN}. Above the limit:
#'     \itemize{
#'       \item [resolution_profile()] stops with an error of class
#'         `spatialkit_work_limit` before the search: on a ladder of sizes,
#'         for the weights of `partition = "density"`; with an `estimate`,
#'         on any ladder, for the estimate's searches.
#'       \item [cell_values_at()] stops with the same error before the
#'         search.
#'       \item [summarize_by_cell()] returns `NA` pooled columns with a
#'         warning of class `spatialkit_pool_unavailable` (`reason = "work
#'         limit"`), and does not run the check of the count record
#'         (`"not checked: work limit"`).
#'       \item [select_resolution()] hands on the centres fitted to the
#'         subsample without refitting them to every point, with a logged
#'         warning.
#'     }}
#'   \item{`spatialkit.fit_work_max`}{Default `1e10` point-to-centre
#'     comparisons. The most work the partition fits of a ladder of sizes
#'     may do: the points times the cells of the finest level times
#'     `nstart`, times `iter_max` for `partition = "density"` without
#'     \pkg{FNN}. With \pkg{FNN} and `init = "random"` the density kind has
#'     no step that compares every point with every centre, and nothing is
#'     counted. Above the limit, [resolution_profile()] stops with an error
#'     of class `spatialkit_work_limit` before any fit. At the default size
#'     the work grows with the square of the rows. [select_resolution()]
#'     refits the centres of a profile fitted to a subsample to every point,
#'     and counts the points times the cells of the count it selects, times
#'     `iter_max` for the density kind without \pkg{FNN}, whatever the
#'     start. Above the limit it hands on the centres fitted to the
#'     subsample, with a logged warning.}
#' }
#'
#' @section Setting a limit:
#' Each option takes a single number, 0 or more; `Inf` removes the limit.
#' Any other value is an error when the option is read. The error of class
#' `spatialkit_work_limit` carries `work` (the count), `limit` and `what`
#' (the step), and its message says what to change: for a fit, fewer
#' starts, a coarser `obs_per_cell`, random starts of the density kind or a
#' higher limit; for a search, reinstalling \pkg{FNN}. A
#' profile just below `spatialkit.fit_work_max` prints nothing while its
#' fits run.
#'
#' The console log is set with [spatialkit_quiet()], not with an option.
#'
#' @seealso [spatialkit-conditions] for the classes named here, and
#'   [resolution_profile()] for the work its fits do.
#' @family options and caches
#' @examples
#' getOption("spatialkit.exact_pairs_max", 2e9)
#' getOption("spatialkit.fit_work_max", 1e10)
#'
#' # Raise one limit for a session, then put it back.
#' old <- options(spatialkit.fit_work_max = 2e10)
#' getOption("spatialkit.fit_work_max")
#' options(old)
#' @name spatialkit-options
NULL
