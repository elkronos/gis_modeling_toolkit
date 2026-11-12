#' @section Where to start:
#' \code{vignette("getting-started", package = "spatialkit")} takes a layer of
#' points through the pipeline below, from a cell count to a cross-validated
#' model and a map that marks where its predictions apply.  Each step of the
#' pipeline names the article that covers it.
#'
#' Ten scripts, 01 to 10, each on one topic, are installed with the package;
#' \code{00-run-all.R} runs them in order.  Each prints what it is doing and
#' says what to look for in a figure before drawing it:
#' \preformatted{
#' dir <- system.file("scripts", package = "spatialkit")
#' list.files(dir)
#' source(file.path(dir, "03-folds.R"))    # one topic
#' source(file.path(dir, "00-run-all.R"))  # 01 to 10 in order
#' }
#' Set \code{SPATIALKIT_TOUR_OUTPUT} to a folder to write the figures there
#' instead of drawing them, and \code{SPATIALKIT_TOUR_PAUSE} to \code{"no"}
#' to skip the pause between figures that an interactive session gets.
#'
#' Three help topics hold what many functions share:
#' \link{spatialkit-geometry} (how input geometry and CRSs are read),
#' \link{spatialkit-conditions} (every warning and error class the package
#' raises, and how to catch one) and \link{spatialkit-options} (the two work
#' limits).
#'
#' @section The pipeline, in order:
#' Steps 1 to 5 build cells from the spatial structure of the observations
#' and summarise the observations onto them.  Steps 6 to 10 fit a model to
#' the points and check it: spatial cross-validation scores the model on
#' blocks held away from its training data, and the area of applicability
#' marks where its predictions leave the range of that data.  Use the steps
#' your question needs, in this order.
#'
#' \enumerate{
#'   \item \strong{Prepare the data.} \code{\link{ensure_projected}()},
#'     \code{\link{harmonize_crs}()}, \code{\link{coerce_to_points}()} and
#'     \code{\link{prep_model_data}()} put the observations in a projected
#'     CRS as points with complete values.  Article:
#'     \code{vignette("getting-started")}.
#'   \item \strong{Measure spatial autocorrelation.}
#'     \code{\link{estimate_sac_range}()} estimates the distance over which
#'     observations stay correlated, and \code{\link{sac_nugget}()} returns
#'     the nugget of its variogram.  The range sets a floor on the number of
#'     cells, which \code{\link{select_resolution}(criterion = "reliability")}
#'     and \code{\link{resolution_profile}(range_floor = TRUE)} apply.  Sizing
#'     blocks from it is one way to set how far the held-out
#'     points sit from their training data, and suits a map that has to
#'     reach that far from the data.  Article:
#'     \code{vignette("spatial-cross-validation")}.
#'   \item \strong{Choose a resolution.} \code{\link{resolution_profile}()}
#'     scores candidate cell counts and \code{\link{select_resolution}()}
#'     reads the count off the profile; by default the count gives 4.5
#'     observations per cell.  \code{\link{determine_optimal_levels}()} reads
#'     a count off the elbow of the k-means WSS curve of the coordinates
#'     alone, and returns \code{NA} when the curve has no elbow.  Article:
#'     \code{vignette("resolution")}.
#'   \item \strong{Build the cells.} \code{\link{build_tessellation}()} draws
#'     Voronoi, hex, square or Delaunay cells with reproducible identifiers;
#'     \code{\link{get_voronoi_seeds}()} places the Voronoi seeds and
#'     \code{\link{clip_target_for}()} builds the polygon that hex and
#'     square cells are clipped to.  Article: \code{vignette("cells")}.
#'   \item \strong{Assign and summarise.}
#'     \code{\link{assign_features_to_polygons}()} puts each observation in
#'     one cell, and \code{\link{summarize_by_cell}()} reduces to one row per
#'     cell with a count, a standard error and a pooled value.
#'     \code{\link{cell_values_at}()} reads cell values at other locations,
#'     and \code{\link{kriging_adequacy}()} checks the cell mean against
#'     block kriging.  Article: \code{vignette("cells")}.
#'   \item \strong{Make folds.} \code{\link{make_folds}()} builds spatial
#'     cross-validation folds: blocked, buffered, leave-location-out or
#'     nearest-neighbour distance-matched.  A cross-validated error is the
#'     error at the distance the folds put between the held-out points and
#'     their training data: it reads low for predictions made farther from
#'     the data and high for predictions made closer, so the folds to use
#'     are the ones matched to where the model will predict.
#'     \code{\link{fold_separation}()} with \code{prediction_points}
#'     compares the two distances.  \code{\link{cv_block_size_sweep}()}
#'     scores a range of block sizes.
#'     Article: \code{vignette("spatial-cross-validation")}.
#'   \item \strong{Fit a model.} \code{\link{fit_rf_model}()} (random
#'     forest), \code{\link{fit_gwr_model}()} (GWR, coefficients that vary
#'     across the map) and \code{\link{fit_bayesian_spatial_model}()}
#'     (Bayesian GP) return a \code{spatial_fit} with \code{predict()},
#'     \code{fitted()}, \code{residuals()}, \code{summary()} and
#'     \code{plot()} methods, and \code{coef()} for GWR and the Bayesian GP.
#'     \code{\link{new_spatial_fit}()} builds a backend of your own.
#'     \code{\link{model_metrics}()}, \code{\link{evaluate_insample}()},
#'     \code{\link{compare_models}()} and \code{\link{residual_morans_i}()}
#'     check a fit in sample.  Article: \code{vignette("models")}.
#'   \item \strong{Cross-validate and compare.} \code{\link{cv_rf}()},
#'     \code{\link{cv_gwr}()}, \code{\link{cv_bayes}()} and
#'     \code{\link{cv_spatial}()}, which takes any learner, score a model on
#'     held-out blocks; \code{\link{compare_models_cv}()} scores several
#'     backends on one set of folds.  Articles:
#'     \code{vignette("spatial-cross-validation")} and
#'     \code{vignette("models")}.
#'   \item \strong{Choose predictors.}
#'     \code{\link{select_features_forward}()} chooses predictors inside the
#'     cross-validation, and \code{\link{gwr_model_selection}()} chooses GWR
#'     predictors by forward selection.  Article: \code{vignette("models")}.
#'   \item \strong{Predict a surface and check where it applies.}
#'     \code{\link{predict_surface}()} predicts onto a regular grid, and
#'     \code{\link{area_of_applicability}()} marks where the predictors leave
#'     the range the model was trained and scored on.  Article:
#'     \code{vignette("models")}.
#'   \item \strong{Report.} \code{\link{ensure_stable_poly_id}()} gives the
#'     cells identifiers that do not depend on the order of the rows.  Article:
#'     \code{vignette("reporting")}.
#' }
#'
#' \code{vignette("spatialkit_nc_demo")} runs the pipeline on one simulated
#' layer.
#'
#' Two choices are the user's, and each has an argument with a default.
#' \code{estimate} is the cell value to report: the plain mean of a cell's
#' observations (\code{"mean"}), or the pooled value (\code{"pooled"}, the
#' default), which is that mean pulled toward the mean of the observations
#' in the 8 nearest cells (nearest by the centroid of their observations).
#' \code{\link{resolution_profile}()} scores a
#' count for the value it is given, and \code{\link{summarize_by_cell}()}
#' reports the value it is given and warns when the count was chosen for
#' the other one.  \code{obs_per_cell} is the size: how many observations a
#' cell is to hold on average, 4.5 by default.  For pooled values the size
#' sets the count (\code{\link{select_resolution}(criterion = "size")}); for
#' plain means the lowest leave-one-out error can set it
#' (\code{criterion = "loo"}).  \code{select_resolution()}'s default
#' criterion, \code{"auto"}, takes the rule that goes with the profile's
#' estimate.  The defaults that change results from 2.0.0, and the argument
#' that restores each, are listed under Breaking changes in
#' \code{news(package = "spatialkit")}.
#'
#' @section Defaults and their sources:
#' Defaults that follow a reference, each cited on the help page of the
#' function named:
#' \itemize{
#'   \item \code{fit_rf_model(include_coords = FALSE)}, and \code{fitted()} on a
#'     forest returning out-of-bag predictions: Meyer et al. (2019).
#'   \item Permutation importance over impurity importance in
#'     \code{fit_rf_model()}: Strobl et al. (2007).
#'   \item The Gaussian-process basis count and boundary factor derived from
#'     the length-scale-to-domain ratio in
#'     \code{fit_bayesian_spatial_model()}: Riutort-Mayol et al. (2023).
#'   \item The nearest-neighbour distance-matching folds of
#'     \code{make_folds(method = "nndm")} and their \code{min_train = 0.5}:
#'     Milà et al. (2022).
#'   \item The area-of-applicability threshold as the outlier-removed maximum
#'     of the training dissimilarity, as the paper defines it (the reference
#'     implementation, CAST, uses the outlier fence itself, which is larger
#'     whenever a training value lies above it), with importance weights
#'     applied directly, without taking their square root, as CAST does:
#'     Meyer and Pebesma (2021).
#'   \item Cliff and Ord moments for the residual Moran's I in
#'     \code{residual_morans_i()}, with \code{null = "auto"}: Cliff and Ord
#'     (1981).
#' }
#'
#' Defaults chosen without a citation:
#' \itemize{
#'   \item The effective range of an exponential variogram as three times its
#'     range parameter, and the identifiability guard against ranges beyond
#'     half the maximum separation, in \code{estimate_sac_range()}.
#'   \item The 25 k-means++ restarts per candidate level in the sweep of
#'     \code{determine_optimal_levels()}: k-means++ seeding follows Arthur
#'     and Vassilvitskii (2007) and a fixed budget follows Fränti and
#'     Sieranoja (2019); 25 is the package's choice, not a published
#'     figure.
#'   \item \code{max_levels = 12} and the unit-step ladder of candidate
#'     levels in \code{determine_optimal_levels()}.
#'   \item The 3:1 fold-imbalance tolerance (\code{balance_tol = 3}) in
#'     \code{make_folds(method = "block_kfold")}.
#'   \item The 1 percent area-distortion tolerance below which a CRS is
#'     taken as equal-area by \code{ensure_projected(purpose = "area")} and
#'     \code{summarize_by_cell(area = TRUE)}.
#'   \item \code{deff_max_n = 500}, \code{sample_n = 1500} and
#'     \code{top_n = 3}, the subsampling and shortlist sizes.
#'   \item \code{coverage_levels = c(0.50, 0.80, 0.95)} in \code{cv_bayes()}.
#'   \item The condition-index cut of 30 in the collinearity checks of
#'     \code{fit_gwr_model()}, a conventional rule of thumb.
#'   \item The 10 percent posterior-mass warning threshold on the GP
#'     length-scale in \code{fit_bayesian_spatial_model()}, which is this
#'     package's own operationalisation of a check the reference recommends,
#'     not a figure from the paper.
#'   \item The small-sample rescaling applied with every data-derived design
#'     effect in \code{summarize_by_cell()}: the package's own derivation from
#'     Kish's exchangeable-correlation model, given on that help page.
#'   \item The size of 4.5 observations per cell that
#'     \code{select_resolution(criterion = "size")} and
#'     \code{resolution_profile(obs_per_cell = "auto")} take when none is
#'     given.
#'   \item The density kind of partition, kept from the best of 3
#'     k-means++ starts, in \code{resolution_profile()},
#'     \code{get_voronoi_seeds()} and \code{voronoi_seeds_kmeans()}: a
#'     convention, a kind that puts small cells where observations are dense
#'     and holds the counts per cell close once the cells outnumber the
#'     dense places, from a number of starts that keeps the cost of a fit
#'     to every row down.  \code{partition = "kmeans"} gives plain k-means.
#'   \item The 8 nearest cells a pooled value is pulled toward in
#'     \code{summarize_by_cell(estimate = "pooled")}: a convention, a count
#'     large enough to give a cell a stable target.  It is not an argument.
#'   \item The 5 residuals expected beyond each end of the interval of
#'     \code{summarize_by_cell(obs_level = )}, below which the warning
#'     \code{spatialkit_few_residuals} is raised: a convention for when the
#'     ends of an interval read from quantiles are unstable.
#'   \item The 30 values beyond one per cell below which
#'     \code{spatialkit_pool_few_df} is raised: a convention for when the
#'     spread inside the cells, which every own weight rests on, is estimated
#'     from little.
#'   \item The 2 percent by which the cells have to beat one mean for the
#'     whole layer, in leave-one-out error, for a response to be called
#'     \code{structured} by \code{summarize_by_cell(estimate = "pooled")}:
#'     a convention, and the default of \code{tol} in
#'     \code{select_resolution()}, which sets the same margin there.
#'   \item The two work limits,
#'     \code{getOption("spatialkit.exact_pairs_max", 2e9)} pairs for a search
#'     made without the \pkg{FNN} package and
#'     \code{getOption("spatialkit.fit_work_max", 1e10)} point-to-centre
#'     comparisons for the fits of a ladder that \code{obs_per_cell} sets:
#'     conventions that keep a search or a fit whose cost grows with the
#'     square of the rows or of the cells from starting.  Above one,
#'     \code{resolution_profile()} and \code{cell_values_at()} stop with an
#'     error, and \code{summarize_by_cell()} returns \code{NA} pooled columns
#'     with a warning (see \link{spatialkit-options}).
#' }
#'
#' @seealso
#' Entry points by task:
#' \code{\link{build_tessellation}()} (build cells),
#' \code{\link{summarize_by_cell}()} (summarise observations onto them),
#' \code{\link{make_folds}()} (split the data into cross-validation folds),
#' \code{\link{compare_models_cv}()} (score several models on the same folds),
#' \code{\link{area_of_applicability}()} (find where predictions leave the
#' range of the training data).
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom stats complete.cases dist fitted kmeans lm.fit median pnorm
#' @importFrom stats predict quantile reformulate residuals sd setNames
#' @importFrom methods as is setOldClass
#' @importFrom utils head modifyList
#' @importFrom dplyr .data
NULL
