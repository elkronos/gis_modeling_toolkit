#' @section The pipeline, in order:
#' The package is built around one workflow. Each step names the function that
#' performs it; every step is optional except the ones your question needs.
#'
#' \enumerate{
#'   \item \strong{Choose a resolution.} \code{\link{determine_optimal_levels}()}
#'     reads a cell count out of the spatial structure of the observations,
#'     rather than making you guess one.
#'   \item \strong{Tessellate.} \code{\link{build_tessellation}()} turns the
#'     point pattern into analysis regions --- Voronoi, Delaunay triangles, or a
#'     hex/square grid --- with reproducible cell identifiers.
#'     \code{\link{get_voronoi_seeds}()} controls where Voronoi seeds go.
#'   \item \strong{Assign.} \code{\link{assign_features_to_polygons}()} labels
#'     every observation with the cell it falls in, resolving multi-match ties
#'     explicitly rather than duplicating rows.
#'   \item \strong{Aggregate.} \code{\link{summarize_by_cell}()} reduces to one
#'     row per cell, carrying a standard error and observation count with every
#'     aggregate, and can correct those errors for within-cell autocorrelation.
#'   \item \strong{Fold.} \code{\link{make_folds}()} builds spatial
#'     cross-validation folds --- blocked, buffered, leave-location-out or
#'     nearest-neighbour distance-matched. Random folds flatter autocorrelated
#'     data; these do not.
#'   \item \strong{Fit.} \code{\link{fit_gwr_model}()} for coefficients that
#'     vary across the map, \code{\link{fit_bayesian_spatial_model}()} for an
#'     explicit spatial Gaussian process with calibrated uncertainty, or
#'     \code{\link{fit_rf_model}()} for predictive accuracy. All three return a
#'     \code{spatial_fit} with common \code{predict()}, \code{fitted()},
#'     \code{residuals()}, \code{coef()} and \code{plot()} methods; write your
#'     own backend with \code{\link{new_spatial_fit}()}.
#'   \item \strong{Validate.} \code{\link{cv_gwr}()},
#'     \code{\link{cv_bayes}()}, \code{\link{cv_rf}()} or the model-agnostic
#'     \code{\link{cv_spatial}()} score a model on held-out blocks;
#'     \code{\link{compare_models_cv}()} scores several backends on one set of
#'     folds. \code{\link{residual_morans_i}()} tests whether spatial structure
#'     survives in the residuals, and \code{\link{select_features_forward}()}
#'     chooses predictors inside the cross-validation.
#'   \item \strong{Predict.} \code{\link{predict_surface}()} projects a fit onto
#'     a regular grid; \code{\link{plot_tessellation_map}()} and
#'     \code{\link{plot.spatial_fit}()} draw the results.
#'   \item \strong{Check applicability.} \code{\link{area_of_applicability}()}
#'     flags where that surface extrapolates beyond the training data. A
#'     cross-validation score says nothing about ground the model has never
#'     seen; this is what tells you where the map should not be believed.
#' }
#'
#' Supporting these throughout, \code{\link{ensure_projected}()} and
#' \code{\link{coerce_to_points}()} handle coordinate reference systems and
#' geometry coercion, and \code{\link{estimate_sac_range}()} estimates the
#' distance over which observations remain correlated --- the number that
#' should be setting your block size.
#'
#' @section Where to start:
#' If you are reading a single page, read
#' \code{vignette("spatialkit_nc_demo", package = "spatialkit")}: it runs the
#' whole pipeline above on North Carolina data, with maps at each step.
#'
#' @seealso
#' \code{vignette("spatialkit_nc_demo", package = "spatialkit")} for the worked
#' end-to-end example.
#'
#' Useful entry points by task:
#' \code{\link{build_tessellation}()} (build regions),
#' \code{\link{summarize_by_cell}()} (aggregate to them),
#' \code{\link{make_folds}()} (split them honestly),
#' \code{\link{compare_models_cv}()} (score several models at once),
#' \code{\link{area_of_applicability}()} (find where not to trust the result).
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom stats complete.cases dist fitted kmeans lm.fit median pnorm
#' @importFrom stats predict quantile reformulate residuals sd setNames
#' @importFrom methods as is
#' @importFrom utils head modifyList
#' @importFrom dplyr .data
NULL
