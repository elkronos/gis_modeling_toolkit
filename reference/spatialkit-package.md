# spatialkit: Spatial Tessellation, Modeling, and Cross-Validation Toolkit

Constructs analysis regions from the distribution of the data itself, as
an alternative to aggregating onto administrative boundaries that were
drawn for unrelated purposes. Seeds and builds Voronoi, Delaunay,
hexagonal and square tessellations with reproducible identifiers,
selects a cell count from the spatial structure of the observations,
assigns features to cells, and aggregates to cell level with optional
design-effect corrections so that standard errors account for
within-cell autocorrelation. Also manages coordinate reference systems.
Fits geographically weighted regression (via 'GWmodel'; Lu et al. (2014)
[doi:10.1080/10095020.2014.917453](https://doi.org/10.1080/10095020.2014.917453)
), Bayesian spatial Gaussian process regression (via 'brms', using the
Hilbert space approximation of Riutort-Mayol et al. (2023)
[doi:10.1007/s11222-022-10167-2](https://doi.org/10.1007/s11222-022-10167-2)
) and random forests (via 'ranger', with the permutation importance of
Strobl et al. (2007)
[doi:10.1186/1471-2105-8-25](https://doi.org/10.1186/1471-2105-8-25) ),
each behind one S3 class with consistent predict, fitted, residuals and
plot methods. Provides spatial cross-validation with random, block,
buffered, leave-location-out and nearest-neighbour distance-matched
folds (Mila et al. (2022)
[doi:10.1111/2041-210X.13851](https://doi.org/10.1111/2041-210X.13851)
), forward variable selection, model comparison, prediction onto a
regular surface, and the area of applicability of Meyer and Pebesma
(2021)
[doi:10.1111/2041-210X.13650](https://doi.org/10.1111/2041-210X.13650)
to flag where a fitted model extrapolates beyond its training data.

## The pipeline, in order

The package is built around one workflow. Each step names the function
that performs it; every step is optional except the ones your question
needs.

Read the order as an argument, not a menu. Steps 1 to 4 are the claim:
that regions drawn from the data's own spatial structure are a better
basis for aggregation and modelling than boundaries drawn for other
purposes. Step 6 is what you do with the regions. Steps 5, 7 and 9 are
the evidence layer: a claim that data-drawn regions beat borrowed ones
is only meaningful if it can be checked, and random folds over
autocorrelated data cannot check it. Spatial cross-validation and the
area of applicability are here because they are what stop the central
claim from being unfalsifiable, not because the package is a
cross-validation package.

1.  **Choose a resolution.**
    [`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
    reads a cell count out of the spatial structure of the observations,
    so you do not have to guess one.

2.  **Tessellate.**
    [`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
    turns the point pattern into analysis regions (Voronoi, Delaunay
    triangles, or a hex/square grid) with reproducible cell identifiers.
    [`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)
    controls where Voronoi seeds go.

3.  **Assign.**
    [`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md)
    labels every observation with the cell it falls in, resolving
    multi-match ties explicitly instead of duplicating rows.

4.  **Aggregate.**
    [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
    reduces to one row per cell, carrying a standard error and
    observation count with every aggregate, and can correct those errors
    for within-cell autocorrelation.

5.  **Fold.**
    [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
    builds spatial cross-validation folds: blocked, buffered,
    leave-location-out or nearest-neighbour distance-matched. Random
    folds flatter autocorrelated data; these do not.

6.  **Fit.**
    [`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
    for coefficients that vary across the map,
    [`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)
    for an explicit spatial Gaussian process with calibrated
    uncertainty, or
    [`fit_rf_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md)
    for predictive accuracy. All three return a `spatial_fit` with
    common [`predict()`](https://rdrr.io/r/stats/predict.html),
    [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
    [`residuals()`](https://rdrr.io/r/stats/residuals.html),
    [`summary()`](https://rdrr.io/r/base/summary.html) and
    [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods,
    and [`coef()`](https://rdrr.io/r/stats/coef.html) on the two that
    have coefficients (a forest has none, so
    [`coef()`](https://rdrr.io/r/stats/coef.html) on an `rf_fit` errors
    by design; use `$info$importance`); write your own backend with
    [`new_spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/new_spatial_fit.md).

7.  **Validate.**
    [`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md),
    [`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md),
    [`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md)
    or the model-agnostic
    [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)
    score a model on held-out blocks;
    [`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md)
    scores several backends on one set of folds.
    [`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md)
    tests whether spatial structure survives in the residuals, and
    [`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md)
    chooses predictors inside the cross-validation.

8.  **Predict.**
    [`predict_surface()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict_surface.md)
    projects a fit onto a regular grid;
    [`plot_tessellation_map()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_tessellation_map.md)
    and
    [`plot.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.spatial_fit.md)
    draw the results.

9.  **Check applicability.**
    [`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md)
    flags where that surface extrapolates beyond the training data. A
    cross-validation score says nothing about ground the model has never
    seen; this is what tells you where the map should not be believed.

Supporting these throughout,
[`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
and
[`coerce_to_points()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md)
handle coordinate reference systems and geometry coercion, and
[`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
estimates the distance over which observations remain correlated. That
is the number that should be setting your block size.

## Defaults and their sources

A stated design principle of this package is that defaults follow
current research rather than convention. That is a claim with a
maintenance cost, so this is the list it applies to, in two parts.

Defaults that cite a reference, each traceable to the help page of the
function named:

- `fit_rf_model(include_coords = FALSE)`, and
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) on a forest
  returning out-of-bag predictions: Meyer et al. (2019).

- Permutation importance over impurity importance in
  [`fit_rf_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md):
  Strobl et al. (2007).

- The Gaussian-process basis count and boundary factor derived from the
  length-scale-to-domain ratio in
  [`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md):
  Riutort-Mayol et al. (2023).

- The nearest-neighbour distance-matching folds of
  `make_folds(method = "nndm")` and their `min_train = 0.5`: Mila et al.
  (2022).

- The area-of-applicability threshold as the outlier-removed maximum of
  the training dissimilarity, with importance weights applied directly,
  without taking their square root, matching the reference
  implementation: Meyer and Pebesma (2021).

- The effective range of an exponential variogram as three times its
  range parameter, and the identifiability guard against ranges beyond
  half the maximum separation, in
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md).

- Cliff and Ord moments for the residual Moran's I in
  [`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md),
  with `null = "auto"`.

- The small-sample rescaling applied with every data-derived design
  effect in
  [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md),
  whose derivation and measured coverage are on that help page.

Defaults that were chosen, and are defensible, but do not rest on a
citation. They are listed so that they are not mistaken for the first
kind:

- The 25 k-means++ restarts per candidate level in the sweep of
  [`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md):
  k-means++ seeding follows Arthur and Vassilvitskii (2007) and a fixed
  budget follows Franti and Sieranoja (2019), but 25 is the budget at
  which the measured WSS curve stopped rising between levels, not a
  published figure.

- `max_levels = 12` and the unit-step ladder of candidate levels in
  [`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md).

- The 3:1 fold-imbalance tolerance (`balance_tol = 3`) in
  `make_folds(method = "block_kfold")`.

- The 1 percent area-distortion tolerance below which a CRS is taken as
  equal-area by `ensure_projected(purpose = "area")` and
  `summarize_by_cell(area = TRUE)`: measured, a UTM zone edge to edge is
  0.25 percent and a continent forced into one zone 14 percent, and the
  figure sits in the gap.

- `deff_max_n = 500`, `sample_n = 1500` and `top_n = 3`, the subsampling
  and shortlist sizes.

- `coverage_levels = c(0.50, 0.80, 0.95)` in
  [`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md).

- The condition-index cut of 30 in the collinearity checks of
  [`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md),
  a conventional rule of thumb.

- The 10 percent posterior-mass warning threshold on the GP length-scale
  in
  [`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md),
  which is this package's own operationalisation of a check the
  reference recommends, not a figure from the paper.

## Where to start

If you are reading a single page, read
[`vignette("getting-started", package = "spatialkit")`](https://elkronos.github.io/gis_modeling_toolkit/articles/getting-started.md):
it takes an `sf` layer of points through every step above, from a cell
count to a cross-validated model and a map with its area of
applicability, on North Carolina data, needing only ranger for the model
and ggplot2 for the maps beyond the hard dependencies.
[`vignette("spatialkit_nc_demo", package = "spatialkit")`](https://elkronos.github.io/gis_modeling_toolkit/articles/spatialkit_nc_demo.md)
is the longer worked example, with four tessellations, both fold schemes
side by side and a GWR fit; it takes its cell counts as given, and
[`vignette("resolution", package = "spatialkit")`](https://elkronos.github.io/gis_modeling_toolkit/articles/resolution.md)
is where those are argued for.

If you would rather run something, ten numbered scripts are installed
with the package. Each prints what it is doing and says what to look for
in a figure before drawing it:


    dir <- system.file("scripts", package = "spatialkit")
    list.files(dir)
    source(file.path(dir, "03-folds.R"))    # one topic
    source(file.path(dir, "00-run-all.R"))  # all ten

Set `SPATIALKIT_TOUR_OUTPUT` to a folder to write the figures there
instead of drawing them, and `SPATIALKIT_TOUR_PAUSE` to `"no"` to skip
the pause between figures that an interactive session gets.

## See also

[`vignette("getting-started", package = "spatialkit")`](https://elkronos.github.io/gis_modeling_toolkit/articles/getting-started.md)
for the pipeline on one page, and
[`vignette("spatialkit_nc_demo", package = "spatialkit")`](https://elkronos.github.io/gis_modeling_toolkit/articles/spatialkit_nc_demo.md)
for the longer worked example.

Useful entry points by task:
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
(build regions),
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
(aggregate to them),
[`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
(split them without flattering the model),
[`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md)
(score several models at once),
[`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md)
(find where not to trust the result).

## Author

**Maintainer**: Justin Chase <jchase.msu@gmail.com> \[copyright holder\]
