# spatialkit: Spatial Tessellation, Modelling, and Cross-Validation Toolkit

Constructs analysis regions from the distribution of the data itself, as
an alternative to aggregating onto administrative boundaries that were
drawn for unrelated purposes. Seeds and builds Voronoi, Delaunay,
hexagonal and square tessellations with reproducible identifiers.
Chooses a cell count, by default from a target number of observations
per cell, or by scoring a ladder of counts against the spatial structure
of the observations. Assigns features to cells and reports for each cell
a count, a mean with its standard error and a pooled value pulled toward
the mean of the observations in the 8 nearest cells. A design effect can
correct the standard errors for autocorrelation within the cells.
Projects lon/lat input to a projected coordinate reference system chosen
for the extent of the data. Fits geographically weighted regression (its
kernels, diagnostics and forward selection follow 'GWmodel'; Lu et al.
(2014)
[doi:10.1080/10095020.2014.917453](https://doi.org/10.1080/10095020.2014.917453)
), Bayesian Gaussian process regression (via 'brms', using the Hilbert
space approximation of Riutort-Mayol et al. (2023)
[doi:10.1007/s11222-022-10167-2](https://doi.org/10.1007/s11222-022-10167-2)
) and random forests (via 'ranger', with the permutation importance of
Strobl et al. (2007)
[doi:10.1186/1471-2105-8-25](https://doi.org/10.1186/1471-2105-8-25) ),
each behind one S3 class with consistent predict, fitted, residuals and
plot methods. Provides spatial cross-validation with random, block,
buffered, leave-location-out and nearest-neighbour distance-matched
folds (Milà et al. (2022)
[doi:10.1111/2041-210X.13851](https://doi.org/10.1111/2041-210X.13851) ;
Linnenbrink et al. (2024)
[doi:10.5194/gmd-17-5897-2024](https://doi.org/10.5194/gmd-17-5897-2024)
), forward variable selection, model comparison, prediction onto a
regular surface, and the area of applicability of Meyer and Pebesma
(2021)
[doi:10.1111/2041-210X.13650](https://doi.org/10.1111/2041-210X.13650)
to flag where a fitted model extrapolates beyond its training data.

## Where to start

[`vignette("getting-started", package = "spatialkit")`](https://elkronos.github.io/gis_modeling_toolkit/articles/getting-started.md)
takes a layer of points through the pipeline below, from a cell count to
a cross-validated model and a map that marks where its predictions
apply. Each step of the pipeline names the article that covers it.

Ten scripts, 01 to 10, each on one topic, are installed with the
package; `00-run-all.R` runs them in order. Each prints what it is doing
and says what to look for in a figure before drawing it:


    dir <- system.file("scripts", package = "spatialkit")
    list.files(dir)
    source(file.path(dir, "03-folds.R"))    # one topic
    source(file.path(dir, "00-run-all.R"))  # 01 to 10 in order

Set `SPATIALKIT_TOUR_OUTPUT` to a folder to write the figures there
instead of drawing them, and `SPATIALKIT_TOUR_PAUSE` to `"no"` to skip
the pause between figures that an interactive session gets.

Three help topics hold what many functions share:
[spatialkit-geometry](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md)
(how input geometry and CRSs are read),
[spatialkit-conditions](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-conditions.md)
(every warning and error class the package raises, and how to catch one)
and
[spatialkit-options](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-options.md)
(the two work limits).

## The pipeline, in order

Steps 1 to 5 build cells from the spatial structure of the observations
and summarise the observations onto them. Steps 6 to 10 fit a model to
the points and check it: spatial cross-validation scores the model on
blocks held away from its training data, and the area of applicability
marks where its predictions leave the range of that data. Use the steps
your question needs, in this order.

1.  **Prepare the data.**
    [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md),
    [`harmonize_crs()`](https://elkronos.github.io/gis_modeling_toolkit/reference/harmonize_crs.md),
    [`coerce_to_points()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md)
    and
    [`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)
    put the observations in a projected CRS as points with complete
    values. Article:
    [`vignette("getting-started")`](https://elkronos.github.io/gis_modeling_toolkit/articles/getting-started.md).

2.  **Measure spatial autocorrelation.**
    [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
    estimates the distance over which observations stay correlated, and
    [`sac_nugget()`](https://elkronos.github.io/gis_modeling_toolkit/reference/sac_nugget.md)
    returns the nugget of its variogram. The range sets a floor on the
    number of cells, which
    [`select_resolution`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)`(criterion = "reliability")`
    and
    [`resolution_profile`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)`(range_floor = TRUE)`
    apply. Sizing blocks from it is one way to set how far the held-out
    points sit from their training data, and suits a map that has to
    reach that far from the data. Article:
    [`vignette("spatial-cross-validation")`](https://elkronos.github.io/gis_modeling_toolkit/articles/spatial-cross-validation.md).

3.  **Choose a resolution.**
    [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
    scores candidate cell counts and
    [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
    reads the count off the profile; by default the count gives 4.5
    observations per cell.
    [`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
    reads a count off the elbow of the k-means WSS curve of the
    coordinates alone, and returns `NA` when the curve has no elbow.
    Article:
    [`vignette("resolution")`](https://elkronos.github.io/gis_modeling_toolkit/articles/resolution.md).

4.  **Build the cells.**
    [`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
    draws Voronoi, hex, square or Delaunay cells with reproducible
    identifiers;
    [`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)
    places the Voronoi seeds and
    [`clip_target_for()`](https://elkronos.github.io/gis_modeling_toolkit/reference/clip_target_for.md)
    builds the polygon that hex and square cells are clipped to.
    Article:
    [`vignette("cells")`](https://elkronos.github.io/gis_modeling_toolkit/articles/cells.md).

5.  **Assign and summarise.**
    [`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md)
    puts each observation in one cell, and
    [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
    reduces to one row per cell with a count, a standard error and a
    pooled value.
    [`cell_values_at()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cell_values_at.md)
    reads cell values at other locations, and
    [`kriging_adequacy()`](https://elkronos.github.io/gis_modeling_toolkit/reference/kriging_adequacy.md)
    checks the cell mean against block kriging. Article:
    [`vignette("cells")`](https://elkronos.github.io/gis_modeling_toolkit/articles/cells.md).

6.  **Make folds.**
    [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
    builds spatial cross-validation folds: blocked, buffered,
    leave-location-out or nearest-neighbour distance-matched. A
    cross-validated error is the error at the distance the folds put
    between the held-out points and their training data: it reads low
    for predictions made farther from the data and high for predictions
    made closer, so the folds to use are the ones matched to where the
    model will predict.
    [`fold_separation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fold_separation.md)
    with `prediction_points` compares the two distances.
    [`cv_block_size_sweep()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_block_size_sweep.md)
    scores a range of block sizes. Article:
    [`vignette("spatial-cross-validation")`](https://elkronos.github.io/gis_modeling_toolkit/articles/spatial-cross-validation.md).

7.  **Fit a model.**
    [`fit_rf_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md)
    (random forest),
    [`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
    (GWR, coefficients that vary across the map) and
    [`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)
    (Bayesian GP) return a `spatial_fit` with
    [`predict()`](https://rdrr.io/r/stats/predict.html),
    [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
    [`residuals()`](https://rdrr.io/r/stats/residuals.html),
    [`summary()`](https://rdrr.io/r/base/summary.html) and
    [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods,
    and [`coef()`](https://rdrr.io/r/stats/coef.html) for GWR and the
    Bayesian GP.
    [`new_spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/new_spatial_fit.md)
    builds a backend of your own.
    [`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md),
    [`evaluate_insample()`](https://elkronos.github.io/gis_modeling_toolkit/reference/evaluate_insample.md),
    [`compare_models()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models.md)
    and
    [`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md)
    check a fit in sample. Article:
    [`vignette("models")`](https://elkronos.github.io/gis_modeling_toolkit/articles/models.md).

8.  **Cross-validate and compare.**
    [`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md),
    [`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md),
    [`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
    and
    [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md),
    which takes any learner, score a model on held-out blocks;
    [`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md)
    scores several backends on one set of folds. Articles:
    [`vignette("spatial-cross-validation")`](https://elkronos.github.io/gis_modeling_toolkit/articles/spatial-cross-validation.md)
    and
    [`vignette("models")`](https://elkronos.github.io/gis_modeling_toolkit/articles/models.md).

9.  **Choose predictors.**
    [`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md)
    chooses predictors inside the cross-validation, and
    [`gwr_model_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gwr_model_selection.md)
    chooses GWR predictors by forward selection. Article:
    [`vignette("models")`](https://elkronos.github.io/gis_modeling_toolkit/articles/models.md).

10. **Predict a surface and check where it applies.**
    [`predict_surface()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict_surface.md)
    predicts onto a regular grid, and
    [`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md)
    marks where the predictors leave the range the model was trained and
    scored on. Article:
    [`vignette("models")`](https://elkronos.github.io/gis_modeling_toolkit/articles/models.md).

11. **Report.**
    [`ensure_stable_poly_id()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_stable_poly_id.md)
    gives the cells identifiers that do not depend on the order of the
    rows. Article:
    [`vignette("reporting")`](https://elkronos.github.io/gis_modeling_toolkit/articles/reporting.md).

[`vignette("spatialkit_nc_demo")`](https://elkronos.github.io/gis_modeling_toolkit/articles/spatialkit_nc_demo.md)
runs the pipeline on one simulated layer.

Two choices are the user's, and each has an argument with a default.
`estimate` is the cell value to report: the plain mean of a cell's
observations (`"mean"`), or the pooled value (`"pooled"`, the default),
which is that mean pulled toward the mean of the observations in the 8
nearest cells (nearest by the centroid of their observations).
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
scores a count for the value it is given, and
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
reports the value it is given and warns when the count was chosen for
the other one. `obs_per_cell` is the size: how many observations a cell
is to hold on average, 4.5 by default. For pooled values the size sets
the count
([`select_resolution`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)`(criterion = "size")`);
for plain means the lowest leave-one-out error can set it
(`criterion = "loo"`).
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)'s
default criterion, `"auto"`, takes the rule that goes with the profile's
estimate. The defaults that change results from 2.0.0, and the argument
that restores each, are listed under Breaking changes in
`news(package = "spatialkit")`.

## Defaults and their sources

Defaults that follow a reference, each cited on the help page of the
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
  `make_folds(method = "nndm")` and their `min_train = 0.5`: Milà et al.
  (2022).

- The area-of-applicability threshold as the outlier-removed maximum of
  the training dissimilarity, as the paper defines it (the reference
  implementation, CAST, uses the outlier fence itself, which is larger
  whenever a training value lies above it), with importance weights
  applied directly, without taking their square root, as CAST does:
  Meyer and Pebesma (2021).

- Cliff and Ord moments for the residual Moran's I in
  [`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md),
  with `null = "auto"`: Cliff and Ord (1981).

Defaults chosen without a citation:

- The effective range of an exponential variogram as three times its
  range parameter, and the identifiability guard against ranges beyond
  half the maximum separation, in
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md).

- The 25 k-means++ restarts per candidate level in the sweep of
  [`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md):
  k-means++ seeding follows Arthur and Vassilvitskii (2007) and a fixed
  budget follows Fränti and Sieranoja (2019); 25 is the package's
  choice, not a published figure.

- `max_levels = 12` and the unit-step ladder of candidate levels in
  [`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md).

- The 3:1 fold-imbalance tolerance (`balance_tol = 3`) in
  `make_folds(method = "block_kfold")`.

- The 1 percent area-distortion tolerance below which a CRS is taken as
  equal-area by `ensure_projected(purpose = "area")` and
  `summarize_by_cell(area = TRUE)`.

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

- The small-sample rescaling applied with every data-derived design
  effect in
  [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md):
  the package's own derivation from Kish's exchangeable-correlation
  model, given on that help page.

- The size of 4.5 observations per cell that
  `select_resolution(criterion = "size")` and
  `resolution_profile(obs_per_cell = "auto")` take when none is given.

- The density kind of partition, kept from the best of 3 k-means++
  starts, in
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md),
  [`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)
  and
  [`voronoi_seeds_kmeans()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md):
  a convention, a kind that puts small cells where observations are
  dense and holds the counts per cell close once the cells outnumber the
  dense places, from a number of starts that keeps the cost of a fit to
  every row down. `partition = "kmeans"` gives plain k-means.

- The 8 nearest cells a pooled value is pulled toward in
  `summarize_by_cell(estimate = "pooled")`: a convention, a count large
  enough to give a cell a stable target. It is not an argument.

- The 5 residuals expected beyond each end of the interval of
  `summarize_by_cell(obs_level = )`, below which the warning
  `spatialkit_few_residuals` is raised: a convention for when the ends
  of an interval read from quantiles are unstable.

- The 30 values beyond one per cell below which `spatialkit_pool_few_df`
  is raised: a convention for when the spread inside the cells, which
  every own weight rests on, is estimated from little.

- The 2 percent by which the cells have to beat one mean for the whole
  layer, in leave-one-out error, for a response to be called
  `structured` by `summarize_by_cell(estimate = "pooled")`: a
  convention, and the default of `tol` in
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md),
  which sets the same margin there.

- The two work limits, `getOption("spatialkit.exact_pairs_max", 2e9)`
  pairs for a search made without the FNN package and
  `getOption("spatialkit.fit_work_max", 1e10)` point-to-centre
  comparisons for the fits of a ladder that `obs_per_cell` sets:
  conventions that keep a search or a fit whose cost grows with the
  square of the rows or of the cells from starting. Above one,
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
  and
  [`cell_values_at()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cell_values_at.md)
  stop with an error, and
  [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
  returns `NA` pooled columns with a warning (see
  [spatialkit-options](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-options.md)).

## See also

Entry points by task:
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
(build cells),
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
(summarise observations onto them),
[`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
(split the data into cross-validation folds),
[`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md)
(score several models on the same folds),
[`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md)
(find where predictions leave the range of the training data).

## Author

**Maintainer**: Justin Chase <jchase.msu@gmail.com> \[copyright holder\]

Authors:

- Justin Chase <jchase.msu@gmail.com> \[copyright holder\]
