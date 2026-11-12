# Package index

## Start here

The pipeline in order, the function for each step, and where each
default comes from.

- [`spatialkit`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-package.md)
  [`spatialkit-package`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-package.md)
  : spatialkit: Spatial Tessellation, Modelling, and Cross-Validation
  Toolkit

## Prepare the data

Put the observations in a projected CRS as points, and see how the
package reads input geometry.

- [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
  : Project an sf object to a CRS suited to its extent
- [`harmonize_crs()`](https://elkronos.github.io/gis_modeling_toolkit/reference/harmonize_crs.md)
  : Bring two spatial objects to one CRS
- [`coerce_to_points()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md)
  : Reduce each geometry to one representative point
- [`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)
  : Prepare an sf layer for model fitting
- [`` `[`( ``*`<spatialkit_rows>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/sub-.spatialkit_rows.md)
  : Subset a layer that carries a row record
- [`spatialkit-geometry`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md)
  : How input geometry is read

## Measure spatial autocorrelation

The distance over which observations stay correlated. It can size the
cross-validation blocks (`make_folds(auto_range = TRUE)`) and sets a
floor on the number of cells, which
`select_resolution(criterion = "reliability")` and
`resolution_profile(range_floor = TRUE)` apply.

- [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  [`print(`*`<sac_range>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  : Estimate the spatial autocorrelation range from data
- [`sac_nugget()`](https://elkronos.github.io/gis_modeling_toolkit/reference/sac_nugget.md)
  : Extract the nugget of an estimated autocorrelation range
- [`plot(`*`<sac_range>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.sac_range.md)
  : Plot an estimated spatial autocorrelation range

## Choose a resolution

How many cells the data support, and the counts the chosen criterion
scores within `tol` of its best (the flat region).

- [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
  [`print(`*`<resolution_profile>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
  : Score candidate cell counts on several criteria
- [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  : Choose a cell count from a resolution profile
- [`summary(`*`<resolution_profile>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/summary.resolution_profile.md)
  : Compare every criterion's pick on a resolution profile
- [`plot(`*`<resolution_profile>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.resolution_profile.md)
  : Plot a resolution profile
- [`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
  : Choose a cell count from the elbow of the k-means WSS curve

## Build the cells

Voronoi, hex, square and Delaunay cells. Start with
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md).
Hex and square cells need a boundary;
[`clip_target_for()`](https://elkronos.github.io/gis_modeling_toolkit/reference/clip_target_for.md)
builds one from the points. The seed functions place the points Voronoi
cells are drawn around, and the lower-level builders draw cells of one
kind or number them.

- [`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
  : Build Voronoi, triangle, hex or square cells from points
- [`clip_target_for()`](https://elkronos.github.io/gis_modeling_toolkit/reference/clip_target_for.md)
  : Build the area to tessellate within, from points or a boundary
- [`plot_tessellation_map()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_tessellation_map.md)
  : Plot a tessellation map with optional boundary, seeds, and features

### Place Voronoi seeds

- [`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)
  : Place seeds for a Voronoi tessellation
- [`voronoi_seeds_kmeans()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md)
  : Place seeds at k-means cluster centres
- [`voronoi_seeds_random()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_random.md)
  : Place seeds at random inside a boundary

### Lower-level builders and IDs

- [`create_voronoi_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_voronoi_polygons.md)
  : Build Voronoi cells around points
- [`create_grid_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons.md)
  : Create square or hexagonal grid polygons over a boundary
- [`create_grid_polygons_cached()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons_cached.md)
  : Create and cache grid polygons over a boundary
- [`ensure_stable_poly_id()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_stable_poly_id.md)
  : Give cells IDs that follow their position

## Assign and summarise

Put each observation in a cell, reduce to one row per cell with a count,
a standard error and a pooled value, and read cell values at other
locations.
[`kriging_adequacy()`](https://elkronos.github.io/gis_modeling_toolkit/reference/kriging_adequacy.md)
checks the cell mean against block kriging.

- [`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md)
  : Assign features to cells and attach a cell ID
- [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
  : Summarise assigned points by cell, with standard errors and pooled
  values
- [`cell_values_at()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cell_values_at.md)
  : Read cell values at new locations
- [`kriging_adequacy()`](https://elkronos.github.io/gis_modeling_toolkit/reference/kriging_adequacy.md)
  : Check whether block kriging would improve cell means

## Make folds

Spatial cross-validation folds, the block size, and how far the held-out
points sit from the training data, beside the map’s locations when you
pass them.

- [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  : Create spatial cross-validation folds
- [`plot_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_folds.md)
  : Map a cross-validation fold scheme
- [`fold_separation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fold_separation.md)
  : Measure how far the held-out points sit from their training data
- [`cv_block_size_sweep()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_block_size_sweep.md)
  [`plot(`*`<block_size_sweep>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_block_size_sweep.md)
  : Cross-validate at a ladder of block sizes

## Fit a model

Random forest, GWR and Bayesian GP backends behind one `spatial_fit`
class, the length-scale bounds the GP uses, and the constructor for a
backend of your own.

- [`fit_rf_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md)
  : Fit a random forest with ranger
- [`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
  : Fit a geographically weighted regression (GWR)
- [`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)
  : Fit a Bayesian Gaussian-process regression with brms
- [`gp_lengthscale_bounds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gp_lengthscale_bounds.md)
  : Compute length-scale bounds for a squared-exponential GP
- [`new_spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/new_spatial_fit.md)
  : Build a spatial_fit object for a backend of your own

## Use a fitted model

[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) for any fit;
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
[`residuals()`](https://rdrr.io/r/stats/residuals.html) for each
backend; and [`coef()`](https://rdrr.io/r/stats/coef.html) for GWR and
the Bayesian GP, since a forest has no coefficients. A forest’s fitted
values and residuals are out-of-bag.

- [`print(`*`<spatial_fit>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.spatial_fit.md)
  [`print(`*`<rf_fit>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.spatial_fit.md)
  : Print a fitted model
- [`summary(`*`<spatial_fit>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/summary.spatial_fit.md)
  : Summarise a fitted model
- [`plot(`*`<spatial_fit>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.spatial_fit.md)
  : Plot a fitted spatial model
- [`predict(`*`<gwr_fit>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-predict.md)
  [`predict(`*`<bayesian_fit>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-predict.md)
  [`predict(`*`<rf_fit>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-predict.md)
  : Predict from a fitted model
- [`fitted(`*`<gwr_fit>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-fitted.md)
  [`fitted(`*`<bayesian_fit>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-fitted.md)
  [`fitted(`*`<rf_fit>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-fitted.md)
  : Extract the fitted values of a fitted model
- [`residuals(`*`<gwr_fit>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-residuals.md)
  [`residuals(`*`<bayesian_fit>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-residuals.md)
  [`residuals(`*`<rf_fit>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-residuals.md)
  : Extract the residuals of a fitted model
- [`coef(`*`<gwr_fit>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-coef.md)
  [`coef(`*`<bayesian_fit>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-coef.md)
  [`coef(`*`<rf_fit>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-coef.md)
  : Extract the coefficients of a fitted model

## Check a fit in sample

Metrics on the training data (out-of-bag for a forest), and residual
Moran’s I, which tests for spatial structure the model left behind.
In-sample metrics do not compare across backends; compare backends with
[`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md).
[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md),
[`evaluate_insample()`](https://elkronos.github.io/gis_modeling_toolkit/reference/evaluate_insample.md)
and
[`compare_models()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models.md)
also score rows you held out, passed as `newdata`.

- [`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md)
  : Compute goodness-of-fit metrics for a fitted model
- [`evaluate_insample()`](https://elkronos.github.io/gis_modeling_toolkit/reference/evaluate_insample.md)
  : Compute the metrics of several fitted models
- [`compare_models()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models.md)
  : Compare fitted models in sample
- [`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md)
  : Test the residuals of a fitted model for spatial autocorrelation
- [`print(`*`<morans_i>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.morans_i.md)
  : Print a residual Moran's I result

## Cross-validate and compare

Score a model on held-out blocks, compare backends on the same folds,
and plot the per-fold spread and interval calibration.

- [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)
  : Cross-validate any learner with spatial folds
- [`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md)
  : Cross-validate a random forest with spatial folds
- [`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md)
  : Cross-validate a GWR with spatial folds
- [`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
  : Cross-validate a Bayesian GP with spatial folds
- [`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md)
  : Compare models on the same cross-validation folds
- [`plot_cv_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_cv_metrics.md)
  : Plot one cross-validation metric fold by fold
- [`plot_calibration()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_calibration.md)
  : Plot the interval calibration of a Bayesian cross-validation

## Choose predictors

Forward selection scored on spatial folds for any backend, and an AICc
search over predictor sets for GWR.

- [`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md)
  : Select predictors by forward selection on spatial folds
- [`plot(`*`<feature_selection>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.feature_selection.md)
  : Plot the path of a forward feature selection
- [`gwr_model_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gwr_model_selection.md)
  [`print(`*`<gwr_model_selection>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/gwr_model_selection.md)
  : Select GWR predictors by forward selection on AICc
- [`plot(`*`<gwr_model_selection>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.gwr_model_selection.md)
  : Plot a GWR model selection

## Predict a surface

Predict onto a regular grid, and mark the locations whose predictor
values are too unlike the training data for the cross-validated score to
hold there.

- [`predict_surface()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict_surface.md)
  : Predict a fitted model onto a regular grid
- [`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md)
  : Mark where a model's predictions apply
- [`print(`*`<aoa>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.aoa.md)
  : Print an area-of-applicability result
- [`plot(`*`<aoa>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.aoa.md)
  : Plot the dissimilarity distribution behind an area of applicability

## Set options, catch conditions and clear caches

The console log, the two work limits, every warning and error class the
package raises, and two in-session caches. To free memory, clear the
grid cache that
[`create_grid_polygons_cached()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons_cached.md)
fills; clear a Bayesian fit’s fitted-value cache after editing its
`brmsfit` by hand.

- [`spatialkit_quiet()`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit_quiet.md)
  : Turn the console log off or on
- [`spatialkit-options`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-options.md)
  : Options that bound the work of a search or a fit
- [`spatialkit-conditions`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-conditions.md)
  : Warnings and errors the package raises
- [`clear_grid_cache()`](https://elkronos.github.io/gis_modeling_toolkit/reference/clear_grid_cache.md)
  : Clear the in-session grid cache
- [`clear_fitted_cache()`](https://elkronos.github.io/gis_modeling_toolkit/reference/clear_fitted_cache.md)
  : Clear the cached fitted values of a Bayesian GP fit
