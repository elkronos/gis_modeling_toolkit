#' Calibrate an inverse-gamma prior for a GP length-scale
#'
#' Returns shape and rate for an inverse-gamma placing \code{tail} probability
#' below \code{lower} and \code{tail} above \code{upper}.
#'
#' A half-normal \code{normal(0, sd)} -- the previous choice -- puts its mode at
#' zero, so most of its mass sits at length-scales shorter than the data can
#' identify.  That is precisely where the Hilbert-space GP degenerates: as the
#' length-scale shrinks below what the basis resolves, the marginal-variance /
#' length-scale posterior develops a funnel and the sampler divergences.  An
#' inverse-gamma with both tails pinned is the standard remedy (Betancourt,
#' "Robust Gaussian Process Modeling"), and the bounds needed to calibrate it
#' are already computed by \code{gp_lengthscale_bounds()}.
#'
#' Solved numerically: if \eqn{X \sim InvGamma(a, b)} then
#' \eqn{1/X \sim Gamma(a, rate = b)}, so both tail conditions can be written
#' against \code{pgamma()} and minimised on the log scale.
#'
#' @param lower,upper Positive numerics, \code{lower < upper}.
#' @param tail Target probability in each tail. Default 0.01.
#' @return A list with \code{shape}, \code{scale} and \code{ok}.  These map
#'   directly onto Stan's \code{inv_gamma(shape, scale)}.  \code{ok} is
#'   \code{FALSE} when the solve did not converge to the requested tails, in
#'   which case the caller should fall back.
#' @keywords internal
#' @noRd
.lscale_invgamma <- function(lower, upper, tail = 0.01) {
  bad <- list(shape = NA_real_, scale = NA_real_, ok = FALSE)
  if (!is.finite(lower) || !is.finite(upper) || lower <= 0 || upper <= lower)
    return(bad)

  # P(X < lower) = 1 - pgamma(1/lower, a, rate = b)
  # P(X > upper) =     pgamma(1/upper, a, rate = b)
  obj <- function(par) {
    a <- exp(par[[1]]); b <- exp(par[[2]])
    lo <- 1 - stats::pgamma(1 / lower, shape = a, rate = b)
    hi <- stats::pgamma(1 / upper, shape = a, rate = b)
    (lo - tail)^2 + (hi - tail)^2
  }

  # Start from an inverse-gamma whose mode is the geometric mean of the bounds:
  # mode = b / (a + 1), so with a = 3, b = 4 * sqrt(lower * upper).
  init <- c(log(3), log(4 * sqrt(lower * upper)))
  fit  <- tryCatch(
    stats::optim(init, obj, method = "Nelder-Mead",
                 control = list(maxit = 2000, reltol = 1e-12)),
    error = function(e) NULL
  )
  if (is.null(fit)) return(bad)

  a <- exp(fit$par[[1]]); b <- exp(fit$par[[2]])
  if (!is.finite(a) || !is.finite(b) || a <= 0 || b <= 0) return(bad)

  lo <- 1 - stats::pgamma(1 / lower, shape = a, rate = b)
  hi <- stats::pgamma(1 / upper, shape = a, rate = b)
  # Accept only if both tails landed near target; a sloppy fit would be a
  # worse prior than the half-normal it replaces.
  ok <- is.finite(lo) && is.finite(hi) &&
    abs(lo - tail) < 0.5 * tail && abs(hi - tail) < 0.5 * tail

  list(shape = a, scale = b, ok = ok)
}


#' Read the family name out of whatever brm() would accept
#'
#' \code{brms::brm()} takes a \code{stats} family function, an evaluated family
#' object, a \code{brmsfamily}, or a bare string.  All that is needed here is
#' the name, to decide whether a numeric response is required.
#'
#' @param family Any of the above, or NULL.
#' @return A length-one character, or \code{NA_character_} if the name cannot
#'   be determined (in which case no response-type check is made).
#' @keywords internal
#' @noRd
.brms_family_name <- function(family) {
  tryCatch({
    f <- family
    if (is.function(f)) f <- f()
    if (is.character(f)) return(as.character(f)[[1L]])
    if (is.list(f) && !is.null(f$family)) return(as.character(f$family)[[1L]])
    NA_character_
  }, error = function(e) NA_character_)
}


#' Build the brms gp() term for the spatial GP
#'
#' Kept separate from \code{fit_bayesian_spatial_model()} so the formula text
#' can be tested without compiling a Stan model.
#'
#' \code{scale = FALSE} is deliberate and load-bearing: \code{brms::gp()}
#' otherwise renormalises its covariates so the maximum pairwise distance is 1
#' and reports \code{lscale} in that space, while this package standardises the
#' coordinates itself and derives the length-scale prior, \code{c} and the
#' basis-adequacy threshold in those units.  Two normalisations would put every
#' length-scale quantity in the wrong space.
#'
#' @param gp_k Basis functions per dimension.
#' @param gp_c Boundary factor.
#' @param gp_iso Logical; single shared length-scale (TRUE) or one per axis.
#' @return A length-one character string.
#' @keywords internal
#' @noRd
.gp_formula_term <- function(gp_k, gp_c, gp_iso = FALSE) {
  sprintf("gp(..x, ..y, k = %s, c = %s, scale = FALSE, iso = %s)",
          gp_k, gp_c, if (isTRUE(gp_iso)) "TRUE" else "FALSE")
}


#' Fit a Bayesian spatial regression with a 2D Gaussian Process (via brms)
#'
#' Fits a regression whose residual spatial structure is modelled explicitly, as
#' a Gaussian process over the coordinates, rather than left in the errors.  Two
#' things follow, and they are the reasons to reach for this backend.  First,
#' every quantity comes with a posterior, so predictions carry calibrated
#' intervals instead of point estimates -- score them with
#' \code{\link{cv_bayes}()}, which reports held-out interval coverage and CRPS.
#' Second, the fitted length-scale is itself an estimate of how far the spatial
#' dependence reaches, a number you can read and report.
#'
#' Choose it over \code{\link{fit_gwr_model}()} when you want one global
#' relationship plus an explicit spatial random field, and uncertainty you can
#' defend; choose GWR instead when the question is how a coefficient
#' \emph{varies} across the map.  Choose \code{\link{fit_rf_model}()} when
#' predictive accuracy matters more than an interpretable model and the
#' response is non-linear in the predictors.  The cost here is time: this is
#' full MCMC via 'brms' and Stan, so it is minutes rather than seconds, and the
#' GP is fitted through a reduced-rank basis approximation whose size
#' (\code{gp_k}) trades fidelity against runtime.
#'
#' @param data_sf An sf object with response, predictors, and geometries.
#' @param response_var Response column name.
#' @param predictor_vars Predictor column names.  May be \code{character(0)}
#'   for an intercept-only model: the response is then explained by the
#'   intercept and the spatial Gaussian process alone, which is the right
#'   baseline for asking how much of the surface is spatial structure rather
#'   than covariate effect (and the natural null model for comparing against a
#'   covariate model with \code{\link{compare_models}}).
#' @param family A model family accepted by \code{brms::brm()} (a stats
#'   family function or a brms family object). Default NULL (resolved to
#'   \code{stats::gaussian()}).
#' @param gp_k Positive integer giving the number of GP basis functions
#'   \emph{per dimension}, or NULL (default) to derive it from the
#'   length-scale/domain ratio.  Note that the fitted model carries
#'   \code{gp_k^2} basis functions, not \code{gp_k} (see Details).
#' @param gp_iso Logical; passed to \code{brms::gp(iso = )}.  \code{FALSE}
#'   (the default) fits a separate length-scale per coordinate axis, letting the
#'   model learn any directional structure from the data.  \code{TRUE} fits a
#'   single shared length-scale, which -- because the coordinates are
#'   standardised per axis beforehand -- makes the kernel anisotropic in the
#'   original CRS by whatever ratio \code{sd(X)/sd(Y)} happens to take.  See
#'   Details.
#' @param gp_c Positive numeric boundary factor for the approximate GP, or
#'   NULL (default) to derive it alongside \code{gp_k}.  The boundary must be
#'   wide enough to contain the longest plausible correlation range; a value
#'   that is too small truncates the domain and degrades the approximation for
#'   smooth, long-range surfaces.
#' @param prior Optional brms prior specification. When NULL and
#'   \code{standardize_predictors = TRUE}, weakly informative
#'   \code{normal(0, 5)} priors are set on regression coefficients.
#'   A data-informed GP length-scale prior is always appended
#'   automatically unless \code{prior} already contains an entry with
#'   \code{class = "lscale"}.
#' @param chains Number of MCMC chains. Default 4.
#' @param iter Total iterations per chain. Default 2000.
#' @param warmup Warmup iterations. Default floor(iter/2).
#' @param cores Number of parallel cores.
#' @param seed Integer seed. Default 123.
#' @param backend "auto", "cmdstanr", or "rstan".
#' @param control Named list of sampler controls, \emph{merged} over the
#'   package defaults \code{list(adapt_delta = 0.9, max_treedepth = 12)} rather
#'   than replacing them.  Passing \code{list(max_treedepth = 15)} therefore
#'   keeps \code{adapt_delta = 0.9} -- which matters, because that is exactly
#'   the setting the divergence warning tells you to raise.
#' @param compute_loo Logical; compute PSIS-LOO. Default TRUE.
#' @param standardize_predictors Logical; center and scale numeric predictors
#'   before fitting. Default FALSE. When TRUE, the scaling parameters are
#'   stored in the return value so predictions can be computed correctly.
#' @param check_convergence Logical; after fitting, check for divergences,
#'   low ESS, and high R-hat and issue warnings. Default TRUE.
#' @param pointize Strategy for non-point geometry coercion.
#' @param boundary Optional polygonal sf/sfc for CRS harmonization.
#' @param .already_prepped Logical (internal). If \code{TRUE}, skip the
#'   \code{prep_model_data()} call because the caller has already projected,
#'   coerced, and filtered the data.  The data must then have plain POINT
#'   geometry (an error is raised otherwise).  Used by the CV internals to
#'   avoid a redundant second pass on every fold.  End users should leave
#'   this at the default \code{FALSE}.
#' @details
#' \strong{GP basis count and boundary factor.}
#' \code{brms::gp()} builds a full tensor grid over its covariates, so a term
#' \code{gp(..x, ..y, k = gp_k)} carries \code{gp_k^2} basis functions -- the
#' \code{gp_k} argument is the count \emph{per dimension}, not the total rank.
#' Both \code{gp_k} and \code{gp_c} are therefore chosen from the ratio of the
#' estimated length-scale to the domain extent, following
#' Riutort-Mayol et al. (2023), rather than from the number of observations:
#' \code{gp_c} is set large enough to contain the upper length-scale bound,
#' and \code{gp_k} large enough to resolve the lower one.  The derived value is
#' typically 21-25 per dimension and is largely independent of \code{n}.
#'
#' The domain extent used is the one \code{brms::gp(c = )} itself multiplies:
#' the full pooled range of the column-centred coordinates
#' (\code{brms:::choose_L()}), not the per-axis half-range in which
#' Riutort-Mayol et al. state their inequalities.  Both constraints are really
#' constraints on the boundary \eqn{L = c \times S}, so expressing them in
#' brms's units is what keeps \code{gp_c}, \code{gp_k} and
#' \code{$info$gp_ell_min} describing the basis brms actually builds.  A
#' \code{gp_c} derived on the half-range convention and handed to
#' \code{brms::gp()} produces a boundary twice as wide as intended, against
#' which \code{gp_k} under-resolves by a factor of two.
#'
#' The GP term is built with \code{scale = FALSE}.  \code{brms::gp()} otherwise
#' rescales its covariates so the maximum Euclidean distance between two points
#' is 1, and reports \code{lscale} in that space; since this function already
#' standardises the coordinates, and the length-scale prior, \code{gp_c} and the
#' adequacy check below are all expressed in those standardised units, a second
#' normalisation would leave every length-scale quantity in the wrong units.
#'
#' After fitting, the posterior length-scale is compared against the smallest
#' scale the chosen basis can resolve
#' (\code{1.75 * gp_c * S / gp_k}, stored as \code{$info$gp_ell_min}); a
#' warning is issued when more than 10\% of the posterior mass falls below it,
#' which is the signal that \code{gp_k} should be raised.
#'
#' \strong{Coordinate scaling and anisotropy.}
#' Before fitting the GP, X and Y coordinates are each centred and divided by
#' their own standard deviation.  This is a conditioning step: easting and
#' northing frequently span very different ranges in a projected CRS, and
#' handing Stan raw metres samples poorly.
#'
#' Because the axes are scaled independently, a \emph{single} shared
#' length-scale in the scaled space corresponds to an anisotropic kernel in the
#' original CRS, stretched by whatever ratio \code{sd(X)/sd(Y)} happens to
#' take.  That ratio is a property of how the sampling locations are laid out,
#' not of the process being modelled, so it is not a defensible source of
#' anisotropy.
#'
#' \code{gp_iso = FALSE} (the default) therefore fits one length-scale per
#' axis, letting the model estimate directional structure from the data instead
#' of inheriting it from the standardisation.  Set \code{gp_iso = TRUE} to
#' recover the previous single-length-scale behaviour.
#'
#' Note that \code{gp_iso} does not affect cost: \code{brms::gp()} builds a
#' tensor grid either way, so the model carries \code{gp_k^2} basis functions
#' regardless.  The stored \code{$info$coord_scaling} list records the scaling
#' strategy, and \code{$info$gp_iso} records which kernel was used.
#'
#' @return A \code{bayesian_fit} object (inherits from \code{spatial_fit}).
#'   Supports \code{predict()}, \code{fitted()}, \code{residuals()},
#'   \code{coef()}, \code{summary()}, and \code{model_metrics()}.
#'   Model-specific metadata lives in \code{$info} (coords -- the names of the
#'   scaled coordinate columns handed to \code{brms::gp()}; coord_scaling,
#'   predictor_scaling, gp_k, gp_c, gp_iso, gp_n_basis, gp_ell_min,
#'   gp_lengthscale_bounds -- the \code{c(lower, upper)} the length-scale prior
#'   was calibrated over; gp_lscale_prior -- the length-scale prior
#'   \code{brms::validate_prior()} reports the model will \emph{actually} use,
#'   which is not necessarily the one this function requested (several entries,
#'   semicolon-separated, if brms resolved the axes differently); loo, looic,
#'   convergence_ok,
#'   convergence_diagnostics).  The raw brmsfit is in \code{$engine}.
#' @family model fitting
#' @examples
#' \dontrun{
#' if (requireNamespace("brms", quietly = TRUE)) {
#'   library(sf)
#'   set.seed(1)
#'   n <- 60
#'   dat <- st_as_sf(
#'     data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), elev = rnorm(n)),
#'     coords = c("x", "y"), crs = 32632
#'   )
#'   dat$price <- 10 + 0.01 * st_coordinates(dat)[, 1] + 2 * dat$elev + rnorm(n)
#'   fit <- fit_bayesian_spatial_model(dat, "price", "elev",
#'                                     chains = 2, iter = 500,
#'                                     compute_loo = FALSE)
#'   summary(fit)
#'   head(predict(fit, newdata = dat))
#' }
#' }
#' @references
#' Riutort-Mayol, G., Burkner, P.-C., Andersen, M. R., Solin, A. and
#' Vehtari, A. (2023). Practical Hilbert space approximate Bayesian
#' Gaussian processes for probabilistic programming.
#' \emph{Statistics and Computing} \strong{33}, 1.
#' \doi{10.1007/s11222-022-10167-2}
#' @export
fit_bayesian_spatial_model <- function(
    data_sf, response_var, predictor_vars,
    family      = NULL,
    
    gp_k        = NULL,
    gp_c        = NULL,
    gp_iso      = FALSE,
    prior       = NULL,
    chains      = 4,
    iter        = 2000,
    warmup      = floor(iter / 2),
    cores       = max(1L, parallel::detectCores() - 1L),
    seed        = 123,
    backend     = c("auto", "cmdstanr", "rstan"),
    control     = list(),
    compute_loo = TRUE,
    standardize_predictors = FALSE,
    check_convergence = TRUE,
    pointize    = "auto",
    boundary    = NULL,
    .already_prepped = FALSE
) {
  if (!inherits(data_sf, "sf"))
    stop("fit_bayesian_spatial_model(): `data_sf` must be an sf object.")
  if (!requireNamespace("brms", quietly = TRUE))
    stop("fit_bayesian_spatial_model(): package 'brms' is required.")
  if (!is.character(response_var) || length(response_var) != 1L)
    stop("fit_bayesian_spatial_model(): `response_var` must be a single ",
         "column name.", call. = FALSE)
  if (!is.character(predictor_vars))
    stop("fit_bayesian_spatial_model(): `predictor_vars` must be a character ",
         "vector (character(0) for an intercept-only spatial GP).",
         call. = FALSE)

  # Merge rather than replace.  A user passing control = list(max_treedepth =
  # 15) would otherwise silently lose adapt_delta = 0.9 -- the very setting
  # the divergence warning further down tells them to raise.
  if (!is.list(control))
    stop("fit_bayesian_spatial_model(): `control` must be a named list of ",
         "sampler control arguments.", call. = FALSE)
  control <- utils::modifyList(list(adapt_delta = 0.9, max_treedepth = 12),
                               control)

  # Resolve family lazily.  NOTE: gaussian() comes from stats — brms does
  # NOT export a `gaussian` object (brms::gaussian errors); brm() accepts
  # stats family functions directly.
  if (is.null(family)) family <- stats::gaussian()

  backend <- match.arg(backend)
  if (identical(backend, "auto"))
    backend <- if (requireNamespace("cmdstanr", quietly = TRUE)) "cmdstanr" else "rstan"

  # detectCores() (used in the `cores` default) can return NA on some
  # platforms; collapse NA/invalid values to 1 before handing to brms.
  cores <- .sanitize_core_count(cores)

  # When called from CV internals the data is already prepped; skip the
  # redundant pass to avoid re-projecting, re-coercing geometry, and
  # re-scanning for NAs on every fold.
  if (isTRUE(.already_prepped)) {
    dat_sf <- data_sf
  } else {
    dat_sf <- prep_model_data(
      data_sf = data_sf, response_var = response_var,
      predictor_vars = predictor_vars, boundary = boundary, pointize = pointize
    )
  }

  # Require plain POINT: st_coordinates() on a multi-vertex MULTIPOINT
  # returns one row per vertex, which would silently misalign coordinates
  # with data rows in the cbind() below.  prep_model_data() now coerces
  # MULTIPOINT to POINT, so this only fires for `.already_prepped = TRUE`
  # callers passing unprepped geometry.
  gtypes <- as.character(sf::st_geometry_type(dat_sf, by_geometry = TRUE))
  if (!all(gtypes == "POINT"))
    stop("fit_bayesian_spatial_model(): geometry must be POINT after prep. ",
         "Run prep_model_data() (or coerce_to_points()) first.")

  # ---- Response validation ----
  # Both sibling backends validate the response (GWR rejects binary outcomes,
  # fit_rf_model() rejects non-numeric ones and fewer than two rows); without
  # this a character response reached brms::brm() and failed inside Stan's
  # data block, where nothing points back at the column.
  dat_cols <- sf::st_drop_geometry(dat_sf)
  if (!(response_var %in% names(dat_cols)))
    stop(sprintf("fit_bayesian_spatial_model(): response '%s' is absent from the data.",
                 response_var), call. = FALSE)
  if (nrow(dat_sf) < 2L)
    stop(sprintf("fit_bayesian_spatial_model(): %d usable row(s) after cleaning; at least two are needed to fit a spatial GP.",
                 nrow(dat_sf)), call. = FALSE)
  if (identical(.brms_family_name(family), "gaussian") &&
      !is.numeric(dat_cols[[response_var]]))
    stop(sprintf(paste0("fit_bayesian_spatial_model(): response '%s' is %s, ",
                        "not numeric, but the family is gaussian. Convert it ",
                        "to numeric, or pass a `family` that matches it (e.g. ",
                        "brms::bernoulli() for a 0/1 outcome)."),
                 response_var, class(dat_cols[[response_var]])[1L]),
         call. = FALSE)

  # NOTE: prep_model_data() coerces geometry to points and drops incomplete
  # rows, and projects the CRS *whenever one can be established* --
  # ensure_projected() passes a CRS-less dataset through unchanged when its
  # coordinates do not look like lon/lat, so this is not an unconditional
  # guarantee.  Set the CRS on `data_sf` if a projected fit matters.

  # `..x` and `..y` are this function's own reserved names: they are the
  # columns handed to brms::gp() and the ones predict()/fitted() rebuild from
  # $info$coord_scaling.  A data column already using either name survives the
  # cbind() below as a DUPLICATE, and both dat_df[["..x"]] and brms's
  # gp(..x, ..y) then resolve to the user's column -- so the GP is fitted over
  # arbitrary data while $info$coord_scaling records the transform that was
  # never applied, and predict time writes a different `..x` than was fitted.
  # Nothing downstream can detect this, so refuse it here.
  reserved <- intersect(c("..x", "..y"), names(sf::st_drop_geometry(dat_sf)))
  if (length(reserved) > 0L)
    stop(sprintf(paste0("fit_bayesian_spatial_model(): column(s) %s in ",
                        "`data_sf` collide with the reserved names this ",
                        "function gives the scaled coordinates it hands to ",
                        "brms::gp(). Rename them (they would otherwise be ",
                        "silently used as the GP's coordinates in place of the ",
                        "geometry)."),
                 paste(sQuote(reserved), collapse = ", ")),
         call. = FALSE)

  coords <- sf::st_coordinates(dat_sf)
  if (!all(c("X", "Y") %in% colnames(coords))) colnames(coords)[1:2] <- c("X", "Y")


  # Per-axis (anisotropic) standardization: each coordinate is centered and
  # divided by its own SD.  This means the GP kernel is isotropic in the
  # *scaled* space but anisotropic in the original CRS whenever sd(X) != sd(Y).
  # See @details in the roxygen block above.
  x_raw <- coords[, "X"]
  y_raw <- coords[, "Y"]
  x_center <- mean(x_raw);  x_scale <- max(sd(x_raw), .Machine$double.eps)
  y_center <- mean(y_raw);  y_scale <- max(sd(y_raw), .Machine$double.eps)
  coord_scaling <- list(
    x_center = x_center, x_scale = x_scale,
    y_center = y_center, y_scale = y_scale,
    scaling_type = "anisotropic"
  )
  dat_df <- cbind(
    sf::st_drop_geometry(dat_sf),
    `..x` = (x_raw - x_center) / x_scale,
    `..y` = (y_raw - y_center) / y_scale
  )

  # ---- Data-informed GP length-scale bounds ----
  # gp_lengthscale_bounds() is calibrated for the squared-exponential kernel
  # that brms::gp() uses.  These bounds drive three separate decisions further
  # down: the basis count and boundary factor (.gp_basis_spec()), and the
  # length-scale prior itself -- which is an inv_gamma() with both tails pinned
  # to these bounds whenever .lscale_invgamma() can calibrate one, and a
  # half-normal normal(0, ls_prior_sd) only as the fallback when it cannot.
  # The prior is CHOSEN ~120 lines below and logs its own message there;
  # naming one here would contradict it on a normal run.
  ls_bounds <- gp_lengthscale_bounds(
    cbind(dat_df[["..x"]], dat_df[["..y"]])
  )
  ls_prior_sd <- max(ls_bounds[["upper"]], ls_bounds[["lower"]] * 1.5)
  .log_info(
    "fit_bayesian_spatial_model(): GP length-scale bounds [%.4f, %.4f] on scaled coords.",
    ls_bounds[["lower"]], ls_bounds[["upper"]]
  )

  # ---- GP basis count and boundary factor ----
  # Both are derived from the length-scale/domain ratio (Riutort-Mayol et al.
  # 2023), NOT from n.  brms carries gp_k^2 basis functions here, so a rule
  # that scaled gp_k with sqrt(n) made the approximation cost grow as n while
  # adding no resolution the data supported.  NULL means "derive"; an explicit
  # value passes through untouched, which is the contract the CV internals
  # rely on when a user supplies gp_k per fold.
  gp_spec   <- .gp_basis_spec(cbind(dat_df[["..x"]], dat_df[["..y"]]), ls_bounds)
  gp_k_auto <- is.null(gp_k)
  if (gp_k_auto)      gp_k <- gp_spec$k
  if (is.null(gp_c))  gp_c <- gp_spec$c

  # Validate before the arithmetic below divides by gp_k.  Previously an
  # invalid user value only surfaced as a brms formula error; now it would also
  # produce a nonsensical resolvable-length-scale and a spurious diagnostic.
  if (!is.numeric(gp_k) || length(gp_k) != 1L || !is.finite(gp_k) || gp_k < 2)
    stop("fit_bayesian_spatial_model(): `gp_k` must be a single finite number >= 2.",
         call. = FALSE)
  if (!is.numeric(gp_c) || length(gp_c) != 1L || !is.finite(gp_c) || gp_c <= 1)
    stop("fit_bayesian_spatial_model(): `gp_c` must be a single finite number > 1.",
         call. = FALSE)
  gp_k <- as.integer(gp_k)
  if (gp_c < 1.25)
    .log_warn(
      paste0("fit_bayesian_spatial_model(): gp_c = %.2f is below 1.25, the ",
             "boundary factor brms itself defaults to (c = 5/4) and the floor ",
             "this package derives; the GP boundary may truncate the domain. ",
             "Note brms multiplies c by the full range of the centred ",
             "coordinates, not the half-range."),
      gp_c
    )

  # Smallest length-scale this (k, c) pair can resolve -- the inversion of
  # m >= 1.75 * L / ell with L = c * S.  gp_spec$S is brms's own domain
  # measure (the pooled range of the centred coordinates), so this is the
  # resolution of the basis brms will actually build, not of a notional one.
  # Reported here and re-checked against the posterior after fitting.
  gp_ell_min <- 1.75 * gp_c * gp_spec$S / gp_k
  .log_info(
    paste0("fit_bayesian_spatial_model(): GP basis k = %d per dimension, ",
           "c = %.2f; TOTAL basis functions = %d (n = %d). ",
           "Smallest resolvable length-scale = %.4f on scaled coords."),
    gp_k, gp_c, gp_k^2, nrow(dat_sf), gp_ell_min
  )
  if (isTRUE(gp_spec$capped) && isTRUE(gp_k_auto))
    .log_warn(
      paste0("fit_bayesian_spatial_model(): derived GP basis count was capped ",
             "at %d. The approximation may under-resolve short-range structure; ",
             "pass a larger gp_k explicitly if short-range spatial variation ",
             "matters."),
      gp_k^2
    )

  # Optional predictor standardization — store transform params so
  # predictions can be computed correctly on new data.
  predictor_scaling <- NULL
  if (isTRUE(standardize_predictors) && length(predictor_vars) > 0L) {
    predictor_scaling <- list()
    for (pv in predictor_vars) {
      if (is.numeric(dat_df[[pv]])) {
        pv_mean <- mean(dat_df[[pv]], na.rm = TRUE)
        pv_sd   <- max(sd(dat_df[[pv]], na.rm = TRUE), .Machine$double.eps)
        dat_df[[pv]] <- (dat_df[[pv]] - pv_mean) / pv_sd
        predictor_scaling[[pv]] <- list(center = pv_mean, scale = pv_sd)
      }
    }
    if (length(predictor_scaling) > 0L)
      .log_info("fit_bayesian_spatial_model(): standardized %d numeric predictor(s).",
                length(predictor_scaling))
  }

  # An intercept-only spatial GP: the response is explained by the intercept
  # and the Gaussian process alone.  prep_model_data() accepts character(0)
  # precisely so this is reachable from the public entry point, not only via
  # .already_prepped = TRUE.
  if (length(predictor_vars) == 0L) {
    rhs_terms <- "1"
    .log_info(paste0("fit_bayesian_spatial_model(): no predictors supplied; ",
                     "fitting an intercept-only model with the spatial GP as ",
                     "the only structure."))
  } else {
    base_fml  <- stats::reformulate(termlabels = predictor_vars)
    rhs_terms <- as.character(base_fml)[2L]
  }
  # scale = FALSE is essential, not cosmetic.  brms::gp() defaults to
  # scale = TRUE, which rescales the covariates so the maximum Euclidean
  # distance between any two points is 1 -- and the posterior `lscale` is then
  # reported in THAT space.  This function has already standardised the
  # coordinates itself (lines above), and gp_lengthscale_bounds(), gp_c and
  # gp_ell_min are all derived in those standardised units.  Leaving brms to
  # apply a second, different normalisation puts the length-scale prior and the
  # basis-adequacy check in units the model does not use -- roughly a factor of
  # the max pairwise distance (~4.9 for standardised 2D coordinates), which
  # makes the prior far too diffuse and the adequacy check fire on every fit.
  # With scale = FALSE there is exactly one coordinate scaling, ours.
  gp_term <- .gp_formula_term(gp_k, gp_c, gp_iso)
  fml <- stats::as.formula(sprintf("%s ~ %s + %s", response_var, rhs_terms, gp_term))

  # Build priors in two independent steps:
  # 1. Regression coefficient priors (only when no user prior and predictors
  #    are standardized).
  # 2. GP length-scale prior (always, unless the user's prior already
  #    includes an "lscale" class entry).

  if (is.null(prior)) {
    prior_parts <- list()
    if (isTRUE(standardize_predictors) && length(predictor_vars) > 0L) {
      prior_parts <- c(prior_parts, list(
        brms::set_prior("normal(0, 5)", class = "b")
      ))
      .log_info("fit_bayesian_spatial_model(): using weakly informative normal(0,5) priors on standardized coefficients.")
    }
    prior <- if (length(prior_parts) > 0L) Reduce(`+`, prior_parts) else NULL
  }

  # Always append the data-informed GP length-scale prior unless the

  # user explicitly supplied an lscale prior.
  user_has_lscale <- !is.null(prior) &&
    inherits(prior, "brmsprior") &&
    any(prior$class == "lscale")

  lscale_prior_spec <- NULL
  if (!user_has_lscale) {
    # Prefer an inverse-gamma with both tails pinned to the estimated bounds.
    # A half-normal puts its mode at zero, so most of its mass sits at
    # length-scales shorter than the basis can resolve -- exactly where the
    # Hilbert-space GP funnels and the sampler divergences.
    ig <- .lscale_invgamma(ls_bounds[["lower"]], ls_bounds[["upper"]])

    # %.10g, not %.6f.  %.6f is fixed-point, so a small scale rounds to the
    # literal "0.000000" (3e-7 does) and Stan rejects inv_gamma(a, 0) with an
    # opaque error from deep inside the model block.  Tightly clustered
    # coordinates get there, and .lscale_invgamma()'s `ok` test only checks
    # tail probabilities, never parameter magnitudes -- so guard those here.
    ig_usable <- isTRUE(ig$ok) &&
      is.finite(ig$shape) && is.finite(ig$scale) &&
      ig$shape > 0 && ig$scale > 0
    if (ig_usable) {
      lscale_prior_spec <- sprintf("inv_gamma(%.10g, %.10g)", ig$shape, ig$scale)
      .log_info(
        paste0("fit_bayesian_spatial_model(): GP length-scale prior ",
               "%s (1%% tails at %.4f and %.4f on scaled coords)."),
        lscale_prior_spec, ls_bounds[["lower"]], ls_bounds[["upper"]]
      )
    } else {
      # Bounds too wide (or too degenerate) to pin both tails; a badly fitted
      # inverse-gamma would be worse than the half-normal it replaces.
      # Guard the fallback's own scale for the same reason: normal(0, 0) is
      # not a prior Stan will accept either.
      sd_spec <- if (is.finite(ls_prior_sd) && ls_prior_sd > 0)
        ls_prior_sd else 1
      lscale_prior_spec <- sprintf("normal(0, %.10g)", sd_spec)
      .log_info(
        paste0("fit_bayesian_spatial_model(): could not calibrate an ",
               "inverse-gamma length-scale prior over [%.4f, %.4f]; falling ",
               "back to %s."),
        ls_bounds[["lower"]], ls_bounds[["upper"]], lscale_prior_spec
      )
    }

    # Attach the prior at COEFFICIENT level, not globally.
    #
    # brms::set_prior(spec, class = "lscale") with no `coef` is a *global*
    # prior, and brms only applies a global prior to coefficients that have no
    # individual prior of their own.  Every lscale coefficient always does --
    # brms assigns each one a default inv_gamma() -- so the calibrated prior
    # was silently dropped with the note "The global prior ... of class
    # 'lscale' will not be used in the model as all related coefficients have
    # individual priors already", and Stan received brms's defaults.  That made
    # gp_lengthscale_bounds(), .lscale_invgamma(), the tail calibration and the
    # %.10g guard all dead weight, and $info$gp_lscale_prior asserted a prior
    # the model did not have.  (A global class = "b" prior does stick, because
    # brms's `b` defaults are flat; the asymmetry is specific to lscale.)
    #
    # The coefficient names are read back from brms rather than hard-coded:
    # they embed the covariate names ("gp..x..y..x", "gp..x..y..y") and the
    # count depends on `iso`, so deriving them is the only version-safe way.
    ls_coefs <- tryCatch({
      gp_def <- brms::get_prior(fml, data = dat_df, family = family)
      gp_def$coef[gp_def$class == "lscale" & nzchar(gp_def$coef)]
    }, error = function(e) character(0))

    ls_prior <- if (length(ls_coefs) > 0L) {
      Reduce(`+`, lapply(ls_coefs, function(k)
        brms::set_prior(lscale_prior_spec, class = "lscale", coef = k)))
    } else {
      # No lscale coefficient to attach to (a brms that names them
      # differently, or a formula with no GP term).  Fall back to the global
      # form: it may be ignored, but it is the only thing left to say.
      .log_warn(paste0("fit_bayesian_spatial_model(): brms reported no ",
                       "coefficient-level 'lscale' priors for this formula, so ",
                       "the calibrated length-scale prior is attached globally ",
                       "and brms may not use it. Check $info$gp_lscale_prior ",
                       "against brms::make_stancode()."))
      brms::set_prior(lscale_prior_spec, class = "lscale")
    }
    prior <- if (is.null(prior)) ls_prior else prior + ls_prior
  } else {
    .log_info(
      "fit_bayesian_spatial_model(): user-supplied prior already includes lscale class; skipping automatic GP length-scale prior."
    )
  }

  # Record the prior brms will actually use, not the one that was requested.
  # validate_prior() resolves globals against coefficient-level defaults and
  # returns the full table brms hands to Stan, so reading the lscale rows back
  # out of it is the only way $info$gp_lscale_prior can be trusted.  On any
  # failure the requested spec is kept and said to be unverified.
  gp_lscale_prior_used <- lscale_prior_spec
  vp <- tryCatch(
    suppressWarnings(brms::validate_prior(prior, formula = fml, data = dat_df,
                                          family = family)),
    error = function(e) NULL
  )
  if (!is.null(vp) && is.data.frame(vp) && all(c("class", "coef", "prior") %in% names(vp))) {
    used <- unique(vp$prior[vp$class == "lscale" & nzchar(vp$coef) &
                              nzchar(vp$prior)])
    if (length(used) > 0L) {
      gp_lscale_prior_used <- paste(used, collapse = "; ")
      if (!is.null(lscale_prior_spec) && !all(used == lscale_prior_spec))
        .log_warn(paste0("fit_bayesian_spatial_model(): brms resolved the GP ",
                         "length-scale prior to %s, not the requested %s. ",
                         "$info$gp_lscale_prior records what brms will use."),
                  gp_lscale_prior_used, lscale_prior_spec)
    }
  }

  brm_args <- list(
    formula = fml, data = dat_df, family = family, prior = prior,
    chains = chains, iter = iter, warmup = warmup, cores = cores,
    seed = seed, control = control
  )
  if (backend %in% c("cmdstanr", "rstan")) brm_args$backend <- backend

  fit <- try(do.call(brms::brm, brm_args), silent = TRUE)
  if (inherits(fit, "try-error"))
    stop(sprintf("fit_bayesian_spatial_model(): brms fit failed: %s", as.character(fit)))

  # ---- Convergence diagnostics ----
  convergence_ok <- TRUE
  convergence_diagnostics <- list()
  if (isTRUE(check_convergence) && inherits(fit, "brmsfit")) {
    # Divergent transitions
    np <- tryCatch(brms::nuts_params(fit), error = function(e) NULL)
    if (!is.null(np)) {
      n_divergent <- sum(np$Value[np$Parameter == "divergent__"], na.rm = TRUE)
      convergence_diagnostics$n_divergent <- n_divergent
      if (n_divergent > 0L) {
        convergence_ok <- FALSE
        .log_warn(
          "fit_bayesian_spatial_model(): %d divergent transition(s) detected. Consider increasing adapt_delta (currently %.2f) or reparameterizing the model.",
          n_divergent, control$adapt_delta %||% 0.8
        )
      }
    }
    
    # R-hat
    rhat_vals <- tryCatch({
      rh <- brms::rhat(fit)
      if (is.numeric(rh)) rh else NULL
    }, error = function(e) NULL)
    
    if (!is.null(rhat_vals)) {
      # brms::rhat() returns NaN for a parameter that is constant across
      # draws, which `lprior` frequently is.  Two consequences had to be
      # guarded.  `names(rhat_vals)[rhat_vals > 1.05]` keeps one NA element per
      # NaN entry -- so convergence_ok was set FALSE and the warning named
      # "NA" -- while which() drops them.  And max(..., na.rm = TRUE) over an
      # all-NA vector is -Inf with a warning, not a diagnostic.
      max_rhat <- if (any(is.finite(rhat_vals)))
        max(rhat_vals[is.finite(rhat_vals)]) else NA_real_
      convergence_diagnostics$max_rhat <- max_rhat
      bad_rhat <- names(rhat_vals)[which(rhat_vals > 1.05)]
      if (length(bad_rhat) > 0L) {
        convergence_ok <- FALSE
        .log_warn(
          "fit_bayesian_spatial_model(): %d parameter(s) have R-hat > 1.05 (max = %.3f): %s. The model may not have converged.",
          length(bad_rhat), max_rhat,
          paste(utils::head(bad_rhat, 5), collapse = ", ")
        )
      }
    }
    
    # Effective sample size ratio
    neff_vals <- tryCatch({
      ne <- brms::neff_ratio(fit)
      if (is.numeric(ne)) ne else NULL
    }, error = function(e) NULL)
    
    if (!is.null(neff_vals)) {
      # Same NaN hazard as R-hat above: a constant parameter has no effective
      # sample size.  which() drops the NAs a logical subscript would keep as
      # NA names, and min() is guarded so an all-NaN vector reports NA rather
      # than +Inf.
      min_neff <- if (any(is.finite(neff_vals)))
        min(neff_vals[is.finite(neff_vals)]) else NA_real_
      convergence_diagnostics$min_neff_ratio <- min_neff
      low_neff <- names(neff_vals)[which(neff_vals < 0.1)]
      if (length(low_neff) > 0L) {
        convergence_ok <- FALSE
        .log_warn(
          "fit_bayesian_spatial_model(): %d parameter(s) have effective sample size ratio < 0.1 (min = %.3f): %s. Consider running longer chains.",
          length(low_neff), min_neff,
          paste(utils::head(low_neff, 5), collapse = ", ")
        )
      }
    }

    # ---- GP basis adequacy ----
    # gp_k is derived from a PRIOR length-scale bound computed from inter-point
    # spacing.  If the posterior length-scale lands below what (gp_k, gp_c) can
    # resolve, the Hilbert-space approximation is inadequate and nothing else
    # here would say so.  This is the diagnostic Riutort-Mayol et al. (2023)
    # recommend, and it is what makes a smaller default gp_k safe rather than
    # merely cheaper.
    #
    # brms names GP parameters lscale(gp...) / sdgp(gp...), embedding the
    # covariate names, so match on the prefix rather than a literal term name.
    # The tryCatch keeps a future brms naming change from erroring the fit.
    ls_draws <- tryCatch(
      as.matrix(fit, variable = "^lscale", regex = TRUE),
      error = function(e) NULL
    )
    if (!is.null(ls_draws) && length(ls_draws) > 0L) {
      frac_below <- mean(as.numeric(ls_draws) < gp_ell_min, na.rm = TRUE)
      convergence_diagnostics$gp_lscale_below_resolution <- frac_below
      if (is.finite(frac_below) && frac_below > 0.10)
        .log_warn(
          paste0("fit_bayesian_spatial_model(): %.0f%% of posterior ",
                 "length-scale draws fall below the smallest scale this basis ",
                 "can resolve (%.4f). Consider raising gp_k (currently %d, ",
                 "%d total basis functions)."),
          100 * frac_below, gp_ell_min, gp_k, gp_k^2
        )
    }
  }

  loo_obj <- NULL; looic <- NA_real_
  if (isTRUE(compute_loo)) {
    loo_try <- try(brms::loo(fit), silent = TRUE)
    if (!inherits(loo_try, "try-error")) {
      loo_obj <- loo_try
      est <- try(loo_obj$estimates, silent = TRUE)
      if (!inherits(est, "try-error") && is.matrix(est)) {
        rn <- rownames(est)
        if ("looic" %in% rn)
          looic <- suppressWarnings(as.numeric(est["looic", "Estimate"]))
        else if ("elpd_loo" %in% rn)
          looic <- suppressWarnings(as.numeric(-2 * est["elpd_loo", "Estimate"]))
      }
    } else {
      .log_warn("fit_bayesian_spatial_model(): LOO computation failed.")
    }
  }

  new_spatial_fit(
    subclass       = "bayesian_fit",
    engine         = fit,
    formula        = fml,
    response_var   = response_var,
    predictor_vars = predictor_vars,
    data_sf        = dat_sf,
    info           = list(
      loo                      = loo_obj,
      looic                    = looic,
      coords                   = c("..x", "..y"),
      coord_scaling            = coord_scaling,
      gp_k                     = gp_k,
      gp_c                     = gp_c,
      gp_iso                   = gp_iso,
      gp_lscale_prior          = gp_lscale_prior_used,
      gp_n_basis               = gp_k^2,
      gp_ell_min               = gp_ell_min,
      gp_lengthscale_bounds    = ls_bounds,
      convergence_ok           = convergence_ok,
      convergence_diagnostics  = convergence_diagnostics,
      predictor_scaling        = predictor_scaling
    )
  )
}
