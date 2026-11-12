# Warnings and errors the package raises

Every condition class the package raises, with the functions that raise
it, what it means, what to do and the fields it carries. Each class can
be caught by name.

## Catching a condition by class

Each condition below is an R warning or error whose class vector starts
with the class named, followed by `"warning"` or `"error"` and
`"condition"`. A handler named after the class catches it without
matching the message, and reads its fields with `$`:


    w <- tryCatch(determine_optimal_levels(pts),
                  spatialkit_no_elbow = function(w) w)
    w$k

[`withCallingHandlers()`](https://rdrr.io/r/base/conditions.html) with
`invokeRestart("muffleWarning")`, or
`suppressWarnings(expr, classes = )`, lets the call finish and drops the
warning. `options(warn = 2)` turns every warning into an error.

Some warnings and errors carry no class of the package's own, such as
the refusal of an argument of the wrong type. Catch those with
`warning =` or `error =`.

## Warnings and log lines

The package also writes a log through logger: a trace in a temporary
file, and a console echo of every line at the WARN level and above. A
log line is not an R condition. `tryCatch(warning = )` does not catch
it, and [`suppressWarnings()`](https://rdrr.io/r/base/warning.html) does
not silence it. Where a help page says a function raises a warning, it
means an R warning; where it says the function logs one, it means a log
line.
[`spatialkit_quiet()`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit_quiet.md)
turns the console echo off.

The package writes a classed warning to the log as well as raising it,
so a warning muffled by a handler still prints its log line until the
echo is turned off.

## Condition classes

Each entry gives the class, the functions that raise it, what it means
and what to do, and then the fields a handler reads with `$`. Each is a
warning unless marked as an error. An entry without fields is a
condition that carries only its message.

- `spatialkit_aoa_zero_threshold`. Raised by
  [`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md)
  with `folds`. The dissimilarity threshold is 0: held-out rows have
  exact copies on the training side of their fold. Use folds that keep
  copies together, such as
  `make_folds(method = "leave_location_out", group_var = )`.

- `spatialkit_bound_chosen`. Raised by
  [`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md),
  [`voronoi_seeds_kmeans()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md)
  and
  [`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md),
  given a profile or a selection as the count. The optimum sits on a
  bound of the ladder, so a bound of the ladder set the count. Read
  `edge` in the
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  result;
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
  describes the bounds. Fields: `edge`, `criterion` and `n`.

- `spatialkit_cells_degenerate` (error; a warning from
  [`ensure_stable_poly_id()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_stable_poly_id.md)).
  Raised by
  [`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md),
  [`create_voronoi_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_voronoi_polygons.md),
  [`create_grid_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons.md),
  [`clip_target_for()`](https://elkronos.github.io/gis_modeling_toolkit/reference/clip_target_for.md),
  [`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md),
  [`voronoi_seeds_random()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_random.md),
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md),
  [`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md),
  [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md),
  [`cell_values_at()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cell_values_at.md)
  and
  [`ensure_stable_poly_id()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_stable_poly_id.md).
  s2 refuses a cell or a boundary row in lon/lat even after repair, or a
  boundary in a geographic CRS has a vertex beyond latitude 90. From
  [`ensure_stable_poly_id()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_stable_poly_id.md):
  the repair left a row with no polygon. Transform the layer to a
  projected CRS; for a latitude beyond 90, check its CRS. Fields:
  `rows`, `n` and `reasons`; the error also carries `parent`, `repaired`
  and `s2_forced`, and some steps add `site`, `cause` or
  `features_rows`.

- `spatialkit_cells_differ`. Raised by
  [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md).
  The rows or cells given are not those the count was scored on; `what`
  says which check differs. Summarise the layer the count was chosen on,
  assigned to the cells built from that selection. Fields: `what`,
  `scored` and `found`, and `field` for the count record.

- `spatialkit_constant_predictor`. Raised by
  [`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md),
  and so by every fitter and `cv_*()` function. A predictor takes one
  value in every row kept, so its coefficient cannot be told from the
  intercept. Drop it from `predictor_vars`.

- `spatialkit_constant_response`. Raised by
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md).
  The scored variable is constant to within rounding, so `cp`, `moran_i`
  and `moran_z` are `NA` at every level. Check the response. With
  `predictor_vars`, a response the predictors fit exactly leaves
  constant residuals. Fields: `n`.

- `spatialkit_convergence`. Raised by
  [`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md),
  [`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
  and
  [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md).
  The sampler's diagnostics fail (divergent transitions, R-hat,
  effective sample size or E-BFMI), and `convergence_ok` is `FALSE`.
  Raise `adapt_delta` for divergent transitions; run more iterations for
  R-hat and effective sample size. Fields: `failed`, and from
  cross-validation also `folds` and `fold_messages`.

- `spatialkit_cp_noise`. Raised by
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md).
  `cp` does not take the variogram's nugget as it stands, or has no
  noise variance and is `NA` at every level. Read `attr(, "cp_noise")`
  for the variance `cp` used and its source. Fields: `source`, `value`,
  `level`, `cap`, `uncapped_source` and `uncapped_value`.

- `spatialkit_deff_fallback`. Raised by
  [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md).
  The design effect requested was not the one applied: a `"variogram"`
  request fell back to Kish's, or a number given was refused. Read the
  `deff_method` column; `deff_method == "variogram"` keeps the rows the
  variogram corrected. Fields, for a `"variogram"` request: `requested`
  and `used`.

- `spatialkit_deff_no_structure`. Raised by
  [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md).
  The variogram of the response found no spatial structure, so the
  predictor columns named get a design effect of 1 from a verdict about
  the response. Use `deff = "kish"` to give each predictor its own ICC.
  Fields: `columns` and `reason`.

- `spatialkit_diagnostic_criterion` (error). Raised by
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md).
  `criterion = "moran_z"` was asked for without `diagnostic = TRUE`.
  `moran_z` is a diagnostic of the cell means, which
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  reads only under `diagnostic = TRUE`. Pass `diagnostic = TRUE`, or
  choose another criterion.

- `spatialkit_diagnostic_ranking`. Raised by
  [`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
  with `criterion = "morans_i"`. The candidates are ranked by Moran's z
  of the cell means, a diagnostic of the cells. Use
  `criterion = "geometric"`, or
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
  when the response should drive the count. Fields: `k`, `criterion`,
  `first_computable` and `window`.

- `spatialkit_empty_cells`. Raised by
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md),
  [`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)
  and
  [`voronoi_seeds_kmeans()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md)
  with `partition = "density"`. A centre holds no point, so fewer cells
  are occupied than the count. It happens at high counts, mostly from
  random starts. k-means++ starts (`init = "kmeans++"`) leave a cell
  empty less often. Fields: `level` and `empty`.

- `spatialkit_estimate_mismatch`. Raised by
  [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md).
  The summary reports another cell value than the count was chosen for:
  plain means against pooled values, or the reverse. Report the value
  the count was chosen for, or choose the count for the value you
  report. Fields: `chosen_for`, `reported`, `rule` and `cells`.

- `spatialkit_features_degenerate` (error). Raised by
  [`coerce_to_points()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md),
  [`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md),
  [`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md)
  and
  [`cell_values_at()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cell_values_at.md).
  s2 refuses a polygon of the features in lon/lat, GEOS puts the
  centroid of an invalid feature outside it, or a temporary projection
  leaves a line empty. Transform the features to a projected CRS, choose
  another `mode` in
  [`coerce_to_points()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md),
  or repair them with
  [`sf::st_make_valid()`](https://r-spatial.github.io/sf/reference/valid.html).
  Fields: `rows`, `n`, `reasons`, `parent`, `repaired` and `s2_forced`,
  and some steps add `site`, `crs` or `auto_helps`.

- `spatialkit_few_residuals`. Raised by
  [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
  with `obs_level`. Fewer than five residuals are expected beyond each
  end of the interval (at 0.90, fewer than 100 residuals), so its ends
  are unstable. Ask for a lower `obs_level`, which needs fewer
  residuals. Fields: `columns`, `needed`, `level` and `n`.

- `spatialkit_fold_args_ignored`. Raised by
  [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md),
  [`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md),
  [`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md),
  [`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
  and
  [`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md).
  `block_size` or `auto_range = TRUE` was given beside `folds`, which
  leaves them nothing to size. Drop them, or leave out `folds` so the
  function builds the blocks.

- `spatialkit_gp_adequacy`. Raised by
  [`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md),
  [`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
  and
  [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md).
  The basis of the Bayesian GP is too coarse: the posterior asks for a
  shorter length-scale than it resolves, or the residuals are
  autocorrelated at short range. Follow the message: a finer basis
  (`gp_k`, `gp_resolution`) or, at its cap, kriging or a model of cell
  means. Fields, from cross-validation: `folds` and `fold_messages`.

- `spatialkit_gp_basis_size`. Raised by
  [`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md),
  [`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
  and
  [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md).
  Fewer than one row per ten basis functions, so most basis weights are
  set by their prior. A smaller `gp_k` is cheaper and resolves less;
  with this few rows kriging fits the same kind of model. Fields, from
  cross-validation: `folds` and `fold_messages`.

- `spatialkit_gp_resolution`. Raised by
  [`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md),
  [`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
  and
  [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md).
  A variogram fitted before compiling puts the response's range below
  what the basis resolves, and the response's autocorrelation agrees.
  The check after sampling tests this estimate. The message gives the
  basis size that would resolve the range. Fields, from
  cross-validation: `folds` and `fold_messages`.

- `spatialkit_gwr_extrapolation`. Raised by
  [`predict.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-predict.md),
  [`predict_surface()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict_surface.md),
  [`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md),
  [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md),
  [`cv_block_size_sweep()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_block_size_sweep.md)
  and
  [`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md).
  New locations lie beyond the reach of the widest kernel; they are
  predicted and counted. Where extrapolation is intended, use
  `suppressWarnings(expr, classes = "spatialkit_gwr_extrapolation")`.
  Fields: `n_far`, `n`, `reach`, `reach_text`, `nearest` and `calls`.

- `spatialkit_gwr_not_estimable`. Raised by
  [`predict.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-predict.md),
  [`predict_surface()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict_surface.md),
  [`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md),
  [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md),
  [`cv_block_size_sweep()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_block_size_sweep.md)
  and
  [`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md).
  The local regression at a new location cannot be estimated, so its
  prediction is `NA`.
  [`predict.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-predict.md)
  says when a window cannot be estimated. Fields: `n_na`, `n`,
  `n_empty`, `bw_text` and `calls`.

- `spatialkit_gwr_rank_deficient`. Raised by
  [`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md).
  Some local windows have deficient rank, and their coefficients are
  `NA`; the fitted values are returned. Where many windows hold a
  predictor constant, fit it as one global coefficient (a mixed GWR, as
  [`GWmodel::gwr.mixed()`](https://rdrr.io/pkg/GWmodel/man/gwr.mixed.html)
  fits).

- `spatialkit_gwr_repeated_locations`. Raised by
  [`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
  and
  [`gwr_model_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gwr_model_selection.md).
  Rows share coordinates. AICc counts them as independent, so a
  bandwidth left to the search (`bandwidth = NULL`) is chosen by
  leave-one-location-out error, and the AICc gaps between models make
  what a predictor adds look larger than it is. Read "Observations that
  share a location" in
  [`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md).

- `spatialkit_icc_unavailable`. Raised by
  [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md).
  A column's Kish ICC cannot be estimated, so its standard errors and
  intervals are `NA` in cells of two or more observations.
  `attr(, "icc")$reason` says why. Fields: `columns` and `reason`.

- `spatialkit_k_reduced`. Raised by
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md),
  and the functions that build folds with it. Fewer folds were made than
  asked for. Read `k` in the result; the message gives the reason and a
  setting that works. Fields: `k_requested` and `k`.

- `spatialkit_knndm_short_of_target`. Raised by
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  with `method = "knndm"`. The held-out points sit closer to their
  training data than the prediction points do, so the cross-validated
  error tends to read low. Use fewer folds; each then holds out a larger
  part of the layer. Fields: `held_out_mean` and `target_mean`.

- `spatialkit_ladder_unused`. Raised by
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md).
  A numeric `obs_per_cell` sets the ladder, so the `n_levels`,
  `min_cell_n`, `sample_n` or `range_floor = TRUE` also given are not
  read. Pass `obs_per_cell` or those arguments, not both. Fields:
  `unused` and `obs_per_cell`.

- `spatialkit_lat_beyond_90`. Raised by
  [`coerce_to_points()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md)
  and
  [`ensure_stable_poly_id()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_stable_poly_id.md).
  Rows of a layer in a geographic CRS have a vertex at a latitude beyond
  90, which is no place on the Earth.
  [`coerce_to_points()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md)
  gives them an empty point;
  [`ensure_stable_poly_id()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_stable_poly_id.md)
  keys them as empty polygons and returns them as they came. Check the
  layer's CRS and the order of its coordinates. Fields: `rows` and `n`.

- `spatialkit_moran_skipped`. Raised by
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md).
  A level has more than 5,000 cells holding a scored row while FNN
  cannot be loaded or Matrix is not installed, so `moran_i` and
  `moran_z` are `NA` there. Install Matrix, or reinstall FNN. Fields:
  `limit`, `levels` and `missing`.

- `spatialkit_new_levels`. Raised by
  [`predict.bayesian_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-predict.md).
  Rows carry a level of a categorical predictor that the model was not
  fitted with, so their predictions are `NA`. The warning names the
  predictor, the levels and the number of rows.

- `spatialkit_nndm_out_of_reach`. Raised by
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  with `method = "nndm"` or `"knndm"`. The prediction points sit nearer
  to the data than any fold can hold a point out, so the cross-validated
  error is pessimistic for them. No fold setting closes the gap; read
  the score as an upper bound on the error at those points. Fields:
  `deficit`, `tolerance` and `separated`.

- `spatialkit_no_cover`. Raised by
  [`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
  and
  [`create_voronoi_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_voronoi_polygons.md).
  Seeds came with neither a boundary nor a cover, so their cells cover
  only the seeds' own hull. Pass the study area as `boundary`, or the
  observations as `cover`.

- `spatialkit_no_default_criterion` (error). Raised by
  [`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md),
  [`voronoi_seeds_kmeans()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md)
  and
  [`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md).
  A profile passed as the count has `cp`, the criterion a profile is
  read by, `NA` at every level while another criterion is finite. Pass
  `select_resolution(profile, criterion = )` with a finite criterion, or
  a number.

- `spatialkit_no_elbow` (an error under `on_no_elbow = "error"`). Raised
  by
  [`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md).
  The WSS curve has no elbow, or the ladder is too short to show one, so
  the count it gives is set by the ladder. Choose the count on other
  grounds, such as
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
  with a response. Fields: `k`, `sag` and `structured`.

- `spatialkit_no_scored_rows`. Raised by
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md).
  No row of the selection half has a response (`select_on = "split"`),
  or no row has a response and every predictor, so the columns that read
  those rows are `NA` at every level. `reason` says which; check where
  the rows with a response and every predictor lie. Fields: `reason`,
  `n`, `selection`, `n_response` and `column`.

- `spatialkit_no_structure`. Raised by
  [`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md),
  [`voronoi_seeds_kmeans()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md),
  [`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
  and
  [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
  with `estimate = "pooled"`. One cell for the whole layer does as well
  (`structured` is `FALSE`), so the cell values say little that one mean
  does not. The call goes on with the count. The message says which test
  failed, and whether a size set the count; see "Pooled values" in
  [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md).
  Fields: `criterion`, `n` and `one_cell`, and from
  [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
  also `columns`, `estimate`, `loo`, `structure_p` and `beats`.

- `spatialkit_nonfinite`. Raised by
  [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md).
  `Inf` or `-Inf` in a response or predictor column is treated as
  missing in every statistic of that column; `n` still counts the row.
  The warning names the columns and gives the rows and cells affected.
  Fields: `n` and `columns`.

- `spatialkit_not_converged`. Raised by
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md),
  [`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)
  and
  [`voronoi_seeds_kmeans()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md).
  The start kept was cut off at the iteration limit, so a longer run
  could still move points between cells. Raise `iter_max` (`kmeans_iter`
  in
  [`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)).
  Fields: `level`, `iter` and `iter_max`.

- `spatialkit_not_supported`. Raised by
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md).
  `supported` is `FALSE`: the range floor exceeds the ceiling, fewer
  than two cells meet `min_cell_n`, or the points sit on one location.
  `attr(x, "bounds")` gives the floor, the ceiling and `ceiling_from`.
  Fields: `reason`, `n` and `min_cell_n`, and `obs_per_cell` for a
  ladder of sizes.

- `spatialkit_partition_unused`. Raised by
  [`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)
  and
  [`voronoi_seeds_kmeans()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md).
  Given a profile or a selection, the call also named a `partition`,
  `init` or density setting other than the profile's; the scored centres
  are returned. Pass the count as a number (`n = sel$best`) for a fresh
  k-means. Fields: `kind` and `init`.

- `spatialkit_pool_degenerate`. Raised by
  [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
  with `estimate = "pooled"`. The pooled values carry no pooling: one
  cell, no spread inside any cell, or no signal between cells. `reason`
  says which case applies. Fields: `reason` and `columns`, with `n`,
  `constant` or `sigma2` by case.

- `spatialkit_pool_few_cells`. Raised by
  [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
  with `estimate = "pooled"`. The column has 2 to 9 cells, so each
  target is the mean of all the other cells, and a pooled value can land
  on the far side of the overall mean from its plain mean. With 10 cells
  or more each target is the mean of the observations in the 8 nearest
  cells. Fields: `columns`, `cells` and `own_weight`.

- `spatialkit_pool_few_df`. Raised by
  [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
  with `estimate = "pooled"`. The spread inside the cells (`sigma2`),
  which every own weight rests on, comes from fewer than 30 values
  beyond one per cell. Fewer cells leave more values beyond one per
  cell. Fields: `columns`, `needed`, `df`, `n` and `cells`.

- `spatialkit_pool_rows_skipped`. Raised by
  [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
  with `estimate = "pooled"`. Rows with a value and a cell have no
  location; they count in `n` and the plain columns and take no part in
  the pooled ones. See "Rows without a location" in
  [spatialkit-geometry](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md).
  Fields: `n`, `cells` and `columns`.

- `spatialkit_pool_unavailable`. Raised by
  [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
  with `estimate = "pooled"` and
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
  with an `estimate`. The pooled columns (in a profile, `loo`, `loo_se`
  and `own_weight`) are `NA`; `reason` says why. The plain columns are
  returned. For `"values out of range"`, divide the column by a power of
  ten; for `"work limit"`, see
  [spatialkit-options](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-options.md).
  Fields: `reason` and `columns`, with `n` and `cells`, or `work`,
  `limit` and `cells`, by reason.

- `spatialkit_predict_failed`. Raised by
  [`predict.bayesian_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-predict.md).
  The posterior draw failed, and every prediction is `NA`. The message
  gives the cause.

- `spatialkit_repeated_locations`. Raised by
  [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
  with `obs_level` and
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
  with an `estimate`. Rows share a location, so the interval or `loo`
  describes a further observation at a place already observed. Fields:
  `n`, the rows that share a location, and `total`, the rows read, with
  `columns` from
  [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
  or `column` from
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md).

- `spatialkit_replicate_var_above_nugget`. Raised by
  [`kriging_adequacy()`](https://elkronos.github.io/gis_modeling_toolkit/reference/kriging_adequacy.md).
  The replicates at one location vary between visits more than the
  variogram's nugget allows; their variance is used as it is. Refit the
  variogram with a nugget of at least the replicate variance. Fields:
  `within_var` and `nugget`.

- `spatialkit_response_ceiling`. Raised by
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md).
  Rows that are not scored (no response, or a predictor missing) lowered
  the ceiling on the count. The message names the bound in force.
  Fields: `n`, `n_response`, `ceiling`, `ceiling_from`, `n_sample` and
  `n_scored`.

- `spatialkit_response_imprecise`. Raised by
  [`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md),
  and so by every fitter and `cv_*()` function, and by
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md),
  [`kriging_adequacy()`](https://elkronos.github.io/gis_modeling_toolkit/reference/kriging_adequacy.md)
  and
  [`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md).
  The response lies so far from zero that doubles there are spaced more
  than a millionth of its standard deviation apart, and every mean is
  rounded to that spacing. Subtract a constant from the response first.
  Fields: `spacing` and `spread`.

- `spatialkit_response_scale`. Raised by
  [`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md).
  `prior` sets a prior on the response's scale, so a response far from
  zero, or with a standard deviation far from 1, is handed to brms as it
  is. See "The response's centre and scale" in
  [`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md).

- `spatialkit_rows_dropped`. Raised by
  [`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md),
  and so by every fitter and `cv_*()` function, and by
  [`gp_lengthscale_bounds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gp_lengthscale_bounds.md).
  Rows were dropped for a missing or non-finite value or for having no
  location. The message gives the counts and how far the dropped rows'
  mean response sits from the kept rows'. When that difference is large,
  the fit is biased towards the complete rows. Fields, from
  [`gp_lengthscale_bounds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gp_lengthscale_bounds.md):
  `n` and `which`.

- `spatialkit_rows_unplaced`. Raised by
  [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md),
  [`harmonize_crs()`](https://elkronos.github.io/gis_modeling_toolkit/reference/harmonize_crs.md),
  [`plot_tessellation_map()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_tessellation_map.md),
  and the other functions that move a layer into another CRS. The target
  CRS cannot place a row in full, so the row comes back empty or without
  those vertices. Correct or drop those rows. Fields: `rows` and `n`.

- `spatialkit_sd_overflow`. Raised by
  [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md).
  A cell's standard deviation is beyond the largest double, so its
  `..sd_*` is `Inf` and its standard errors and interval are `NA`.
  Divide the column by a power of ten. Fields: `columns` and `cells`.

- `spatialkit_size_not_reached`. Raised by
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md),
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  and
  [`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md).
  A size does not give the count it asks for (fewer than two cells, or
  more than k-means can place), and the nearest level is used. Fields:
  `obs_per_cell`, `wanted` and `used`, the two counts, and `reason`.

- `spatialkit_surface_point_off`. Raised by
  [`coerce_to_points()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md).
  The point GEOS takes on a lon/lat polygon spanning more than 180
  degrees of longitude does not lie on it read on the sphere, and that
  point is kept. Transform the layer to a projected CRS that holds the
  polygon, and call again. Fields: `rows`.

- `spatialkit_ties_unresolved`. Raised by
  [`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md).
  The `tie_break` rule names no polygon for some features, since the
  polygons are equal on every key; those features get `NA`. They are
  kept or dropped as `keep_unassigned` says. Fields: `n` and `total`.

- `spatialkit_unassigned`. Raised by
  [`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md),
  [`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
  and
  [`create_voronoi_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_voronoi_polygons.md).
  Features fall inside no polygon, or points or seeds lie outside every
  cell. Give the cells a `boundary` or `cover` that reaches the data;
  `keep_unassigned = TRUE` keeps such features with `NA`. Fields: `n`
  and `total`.

- `spatialkit_variogram_unusable`. Raised by
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md).
  The variogram gives no usable range, so `reliability` is `NA`.
  `reason` gives the `rejected_reason` of
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md).
  Fields: `reason`, `converged`, `nugget` and `nugget_used`.

- `spatialkit_work_limit` (error). Raised by
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
  and
  [`cell_values_at()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cell_values_at.md).
  A search or a fit would make more comparisons than an option allows.
  The message says what to change; see
  [spatialkit-options](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-options.md).
  Fields: `work`, `limit` and `what`, and for a fit also `kind`, `init`,
  `nstart` and `fnn`.

- `spatialkit_xy_imprecise`. Raised by
  [`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md),
  [`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md),
  [`create_voronoi_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_voronoi_polygons.md),
  [`create_grid_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons.md),
  [`clip_target_for()`](https://elkronos.github.io/gis_modeling_toolkit/reference/clip_target_for.md),
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md),
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md),
  [`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
  and
  [`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md).
  The coordinates lie so far from the origin of their CRS that doubles
  there are spaced more than a millionth of the extent apart, or a
  distance lies past a tie by less than that rounding. The call runs on
  the rounded coordinates. Subtract a constant from the coordinates; see
  "Coordinates far from the origin" in
  [spatialkit-geometry](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md).
  Fields: `spacing`, and `extent` or `distance`.

- `spatialkit_xy_imprecise_error` (error). Raised by
  [`clip_target_for()`](https://elkronos.github.io/gis_modeling_toolkit/reference/clip_target_for.md).
  At that rounding the buffer comes back empty or with no area. Shift
  the points nearer the origin, or pass `boundary`.

- `spatialkit_zm_dropped`. Raised by
  [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md),
  [`harmonize_crs()`](https://elkronos.github.io/gis_modeling_toolkit/reference/harmonize_crs.md),
  [`coerce_to_points()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md),
  [`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md),
  [`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md),
  [`ensure_stable_poly_id()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_stable_poly_id.md),
  [`predict_surface()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict_surface.md)
  and
  [`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md).
  A layer that mixes XY and XYZ (or XYM) geometries came back in XY, or
  [`coerce_to_points()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md)
  reduced geometries that carry M in x and y. See "Mixed dimensions, Z
  and M" in
  [spatialkit-geometry](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md).
  Fields: `rows`, the rows that carried Z or M.

## See also

[`spatialkit_quiet()`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit_quiet.md)
for the console log,
[spatialkit-options](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-options.md)
for the work limits,
[spatialkit-geometry](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md)
for the rules behind the geometry classes, and
[`vignette("troubleshooting", package = "spatialkit")`](https://elkronos.github.io/gis_modeling_toolkit/articles/troubleshooting.md)
for the messages met most often.

Other options and caches:
[`clear_fitted_cache()`](https://elkronos.github.io/gis_modeling_toolkit/reference/clear_fitted_cache.md),
[`clear_grid_cache()`](https://elkronos.github.io/gis_modeling_toolkit/reference/clear_grid_cache.md),
[`spatialkit-options`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-options.md),
[`spatialkit_quiet()`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit_quiet.md)

## Examples

``` r
library(sf)
set.seed(1)
# Points spread evenly have no elbow.
pts <- st_as_sf(data.frame(x = 5e5 + runif(200, 0, 100),
                           y = 5e6 + runif(200, 0, 100)),
                coords = c("x", "y"), crs = 32632)

# Catch the warning by its class and read its fields.
w <- tryCatch(determine_optimal_levels(pts, max_levels = 6),
              spatialkit_no_elbow = function(w) w)
class(w)
#> [1] "spatialkit_no_elbow" "warning"             "condition"          
w$k
#> [1] 3

# Let the call finish without the warning; the log line still prints
# until the console echo is off.
old <- spatialkit_quiet()
k <- suppressWarnings(determine_optimal_levels(pts, max_levels = 6),
                      classes = "spatialkit_no_elbow")
spatialkit_quiet(old)
k
#> [1] NA
#> attr(,"structured")
#> [1] FALSE
#> attr(,"sag")
#> [1] 0
```
