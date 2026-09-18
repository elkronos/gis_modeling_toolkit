# Package index

## Overview

- [`spatialkit`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-package.md)
  [`spatialkit-package`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-package.md)
  : spatialkit: Spatial Tessellation, Modeling, and Cross-Validation
  Toolkit

## Prepare the data

Coordinate reference systems, geometry coercion, and the clip target.

- [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
  : Ensure an object has a projected CRS (with sensible defaults)
- [`harmonize_crs()`](https://elkronos.github.io/gis_modeling_toolkit/reference/harmonize_crs.md)
  : Harmonize CRS between two spatial objects
- [`coerce_to_points()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md)
  : Coerce arbitrary geometries to representative points
- [`clip_target_for()`](https://elkronos.github.io/gis_modeling_toolkit/reference/clip_target_for.md)
  : Build a polygonal clip target from points and/or a boundary
- [`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)
  : Prepare and sanitize an sf dataset for spatial modeling

## Measure the spatial structure

The distance over which observations stay correlated. It sizes both the
cells and the cross-validation blocks.

- [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  : Estimate the spatial autocorrelation range from data
- [`sac_nugget()`](https://elkronos.github.io/gis_modeling_toolkit/reference/sac_nugget.md)
  : The nugget of an estimated autocorrelation range
- [`print(`*`<sac_range>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.sac_range.md)
  : Print a spatial autocorrelation range
- [`plot(`*`<sac_range>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.sac_range.md)
  : Plot an estimated spatial autocorrelation range

## Choose a resolution

How many cells, and how sure that number is.

- [`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
  : Determine an optimal number of spatial levels via an elbow heuristic
- [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
  : Score every candidate number of cells on several criteria at once
- [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  : Read a level, and the region over which it is not distinguishable,
  off a profile
- [`print(`*`<resolution_profile>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.resolution_profile.md)
  : Print a resolution profile
- [`plot(`*`<resolution_profile>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.resolution_profile.md)
  : Plot a resolution profile

## Tessellate

Build the regions, with reproducible cell identifiers.

- [`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
  : Build a tessellation (Voronoi, Delaunay triangles, hex grid, or
  square grid)
- [`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)
  : Generate seed points for Voronoi tessellation
- [`voronoi_seeds_kmeans()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md)
  : K-means seed generation from point coordinates
- [`voronoi_seeds_random()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_random.md)
  : Random seed generation within a polygonal boundary
- [`create_voronoi_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_voronoi_polygons.md)
  : Create Voronoi polygons from points with robust CRS and optional
  clipping
- [`create_grid_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons.md)
  : Create square or hexagonal grid polygons over a boundary
- [`create_grid_polygons_cached()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons_cached.md)
  : Create and cache grid polygons over a boundary
- [`clear_grid_cache()`](https://elkronos.github.io/gis_modeling_toolkit/reference/clear_grid_cache.md)
  : Clear the in-session grid cache
- [`ensure_stable_poly_id()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_stable_poly_id.md)
  : Create deterministic, stable polygon IDs based on spatial sort keys
- [`plot_tessellation_map()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_tessellation_map.md)
  : Plot a tessellation map with optional boundary, seeds, and features

## Assign and aggregate

One row per cell, with a count and a standard error for every aggregate,
and the record of which rows were dropped or tied on the way.

- [`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md)
  : Assign features to polygons and attach a polygon ID
- [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
  : Summarize features by polygon/cell ID
- [`` `[`( ``*`<spatialkit_rows>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/sub-.spatialkit_rows.md)
  : Subset a layer that carries a row record
- [`kriging_adequacy()`](https://elkronos.github.io/gis_modeling_toolkit/reference/kriging_adequacy.md)
  : Block-kriging adequacy diagnostics for a set of cells

## Fold

Spatial cross-validation folds, and how large the blocks should be.

- [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  : Create spatial cross-validation folds
- [`plot_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_folds.md)
  : Map a cross-validation fold scheme
- [`cv_block_size_sweep()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_block_size_sweep.md)
  : Cross-validate at a ladder of block sizes
- [`plot(`*`<block_size_sweep>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.block_size_sweep.md)
  : Plot cross-validation error against block size

## Fit

Three backends behind one `spatial_fit` class, and the constructor for a
backend of your own.

- [`new_spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/new_spatial_fit.md)
  : Build a spatial_fit S3 object
- [`fit_rf_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md)
  : Fit a random forest via ranger
- [`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
  : Fit a Geographically Weighted Regression (GWR) via GWmodel
- [`gwr_model_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gwr_model_selection.md)
  : Forward model selection for geographically weighted regression
- [`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)
  : Fit a Bayesian spatial regression with a 2D Gaussian Process (via
  brms)
- [`gp_lengthscale_bounds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gp_lengthscale_bounds.md)
  : Heuristic length-scale bounds for a squared-exponential GP
- [`clear_fitted_cache()`](https://elkronos.github.io/gis_modeling_toolkit/reference/clear_fitted_cache.md)
  : Clear cached fitted values for a Bayesian spatial model

## Methods on a fit

- [`print(`*`<spatial_fit>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.spatial_fit.md)
  : Print a fitted spatial model
- [`summary(`*`<spatial_fit>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/summary.spatial_fit.md)
  : Summarise a fitted spatial model
- [`plot(`*`<spatial_fit>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.spatial_fit.md)
  : Plot a fitted spatial model
- [`predict(`*`<rf_fit>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict.rf_fit.md)
  : Predict from a random forest fit
- [`predict(`*`<gwr_fit>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict.gwr_fit.md)
  : Predict from a GWR spatial model
- [`predict(`*`<bayesian_fit>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict.bayesian_fit.md)
  : Predict from a Bayesian spatial GP model
- [`fitted(`*`<rf_fit>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/fitted.rf_fit.md)
  : Out-of-bag predictions from a random forest fit
- [`fitted(`*`<gwr_fit>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/fitted.gwr_fit.md)
  : In-sample fitted values from a GWR fit
- [`fitted(`*`<bayesian_fit>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/fitted.bayesian_fit.md)
  : In-sample fitted values from a Bayesian spatial GP fit
- [`residuals(`*`<rf_fit>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/residuals.rf_fit.md)
  : Out-of-bag residuals from a random forest fit
- [`residuals(`*`<gwr_fit>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/residuals.gwr_fit.md)
  : In-sample residuals from a GWR fit
- [`residuals(`*`<bayesian_fit>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/residuals.bayesian_fit.md)
  : In-sample residuals from a Bayesian spatial GP fit
- [`coef(`*`<rf_fit>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/coef.rf_fit.md)
  : Coefficients are undefined for a random forest
- [`coef(`*`<gwr_fit>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/coef.gwr_fit.md)
  : Extract GWR local coefficients
- [`coef(`*`<bayesian_fit>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/coef.bayesian_fit.md)
  : Extract Bayesian model fixed-effect summaries
- [`print(`*`<rf_fit>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.rf_fit.md)
  : Print a random forest fit
- [`print(`*`<gwr_model_selection>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.gwr_model_selection.md)
  : Print a GWR model selection result
- [`plot(`*`<gwr_model_selection>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.gwr_model_selection.md)
  : Plot a GWR model selection

## Validate

Score on held-out blocks, compare backends on the same folds, and test
the residuals.

- [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)
  : Model-agnostic spatial cross-validation
- [`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md)
  : Cross-validate a random forest with spatial folds
- [`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md)
  : K-fold cross-validation for GWR
- [`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
  : K-fold cross-validation for the Bayesian spatial model
- [`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md)
  : Cross-validated comparison of spatial models
- [`plot_cv_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_cv_metrics.md)
  : Plot one cross-validation metric fold by fold
- [`plot_calibration()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_calibration.md)
  : Plot the interval calibration of a Bayesian cross-validation
- [`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md)
  : Greedy forward feature selection with spatially blocked inner folds
- [`plot(`*`<feature_selection>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.feature_selection.md)
  : Plot the path of a forward feature selection
- [`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md)
  : Compute goodness-of-fit metrics for a spatial model
- [`evaluate_insample()`](https://elkronos.github.io/gis_modeling_toolkit/reference/evaluate_insample.md)
  : Compute in-sample (or out-of-sample) metrics for fitted spatial
  models
- [`compare_models()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models.md)
  : Side-by-side comparison of fitted spatial models
- [`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md)
  : Compute Moran's I on the residuals of a fitted spatial model
- [`print(`*`<morans_i>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.morans_i.md)
  : Print a residual Moran's I result

## Predict, and check where the prediction applies

- [`predict_surface()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict_surface.md)
  : Predict a fitted spatial model onto a regular grid
- [`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md)
  : Area of applicability of a spatial prediction model
- [`print(`*`<aoa>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.aoa.md)
  : Print an area-of-applicability result
- [`plot(`*`<aoa>`*`)`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.aoa.md)
  : Plot the dissimilarity distribution behind an area of applicability

## Package options

- [`spatialkit_quiet()`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit_quiet.md)
  : Quieten (or restore) spatialkit's console log
