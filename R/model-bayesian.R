#' Calibrate an inverse-gamma prior for a GP length-scale
#'
#' Returns shape and rate for an inverse-gamma placing \code{tail} probability
#' below \code{lower} and \code{tail} above \code{upper}.
#'
#' A half-normal \code{normal(0, sd)} (the previous choice) puts its mode at
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
  #
  # The two tail conditions have one exact solution, and it is a 1-D root.
  # With X = 1/Y, Y ~ Gamma(a, rate = b), and pgamma(q, a, rate = b) =
  # pgamma(q * b, a): the upper-tail condition gives b = upper * qgamma(tail,
  # a), and substituting into the lower-tail one leaves
  #     qgamma(1 - tail, a) / qgamma(tail, a) = upper / lower,
  # whose left side is monotone decreasing in a.  A 2-D Nelder-Mead on
  # (log a, log b) from a fixed start failed to reach the tails whenever
  # upper/lower exceeded ~65 -- every clustered layout -- and was sensitive
  # to the absolute scale of the bounds, so the fit silently fell back to the
  # half-normal this very function calls the worse prior.
  ratio <- upper / lower
  f <- function(loga) {
    a <- exp(loga)
    stats::qgamma(1 - tail, shape = a) / stats::qgamma(tail, shape = a) - ratio
  }
  root <- tryCatch({
    # a in [1e-4, 1e4]: ratio(1e-4) is astronomically large, ratio(1e4) is
    # ~1.05, so every ratio > 1.05 is bracketed.
    lo_a <- log(1e-4); hi_a <- log(1e4)
    if (f(hi_a) > 0) hi_a <- log(1e6)
    stats::uniroot(f, c(lo_a, hi_a), tol = 1e-12, maxiter = 1000L)$root
  }, error = function(e) NULL)
  if (is.null(root)) return(bad)
  a <- exp(root)
  b <- upper * stats::qgamma(tail, shape = a)
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


#' Evaluate a convergence accessor without posterior's capped-ESS warning
#'
#' \pkg{posterior} warns "The ESS has been capped to avoid unstable estimates."
#' once for every parameter whose effective sample size exceeds its draw
#' count.  In a Hilbert-space GP that is most of the \code{zgp} basis weights,
#' which mix well precisely because their chains are antithetic: an
#' \code{n = 80} fit raised 41 of them from a single
#' \code{brms::neff_ratio()} call.  R keeps only the first 50 warnings, so a
#' warning that mattered and came later -- loo's Pareto-k, raised after this
#' check -- was dropped.  The capped value is still the one to read; only that
#' message is muffled, and every other warning passes through.  Restricting
#' the check to non-\code{zgp} parameters instead would hide real flags: in the
#' fits the review ran, the only ratios below 0.1 were \code{zgp} weights.
#'
#' @param expr The call to evaluate.
#' @return The value of \code{expr}.
#' @keywords internal
#' @noRd
.muffle_ess_cap <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("The ESS has been capped", conditionMessage(w), fixed = TRUE))
      invokeRestart("muffleWarning")
  })
}


# The brms families whose response is a category rather than a number, and so
# the only ones fit_bayesian_spatial_model() lets a factor or character
# response through for.
.brms_category_families <- c("categorical", "cumulative", "sratio", "cratio",
                             "acat")


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


#' The coefficient-level length-scale parameters brms creates for a model
#'
#' One row per \code{lscale} coefficient, with everything that addresses it:
#' \code{coef}, \code{resp}, \code{dpar} and \code{nlpar}.  The address is more
#' than the \code{coef} name.  A family with several distributional parameters
#' gives each its own GP -- \code{brms::categorical()} one per non-reference
#' category (\code{muhi}, \code{mumid}), a \code{brms::mixture()} one per
#' component -- under the SAME \code{coef} names, so a prior built from
#' \code{coef} alone lands every row on \code{dpar = ""} twice, and brms
#' refuses the whole model with "Duplicated prior specifications are not
#' allowed" before sampling.
#'
#' @param fml,data,family As handed to \code{brms::brm()}.
#' @return A data.frame with character columns \code{coef}, \code{resp},
#'   \code{dpar}, \code{nlpar}, or \code{NULL} when brms reports none (or
#'   \code{brms::get_prior()} fails).
#' @keywords internal
#' @noRd
.lscale_coef_rows <- function(fml, data, family) {
  gp_def <- tryCatch(brms::get_prior(fml, data = data, family = family),
                     error = function(e) NULL)
  if (is.null(gp_def) || !all(c("class", "coef") %in% names(gp_def)))
    return(NULL)
  keep <- gp_def$class == "lscale" & nzchar(gp_def$coef)
  if (!any(keep)) return(NULL)
  field <- function(nm) {
    v <- if (nm %in% names(gp_def)) as.character(gp_def[[nm]][keep])
         else rep("", sum(keep))
    v[is.na(v)] <- ""
    v
  }
  unique(data.frame(coef = field("coef"), resp = field("resp"),
                    dpar = field("dpar"), nlpar = field("nlpar"),
                    stringsAsFactors = FALSE))
}


#' One coefficient-level lscale prior per row of .lscale_coef_rows()
#'
#' @param spec The prior, as a Stan distribution string.
#' @param rows A (non-empty) return value of \code{.lscale_coef_rows()}.
#' @return A \code{brmsprior}.
#' @keywords internal
#' @noRd
.lscale_prior_rows <- function(spec, rows) {
  Reduce(`+`, lapply(seq_len(nrow(rows)), function(i)
    brms::set_prior(spec, class = "lscale", coef = rows$coef[[i]],
                    resp = rows$resp[[i]], dpar = rows$dpar[[i]],
                    nlpar = rows$nlpar[[i]])))
}


#' The weakly informative slope prior, one row per distributional parameter
#'
#' \code{set_prior(spec, class = "b")} with no \code{dpar} addresses only the
#' slopes of the main parameter.  A \code{brms::categorical()} or
#' \code{brms::mixture()} family has none of those -- its slopes belong to
#' \code{muhi}/\code{mumid} or \code{mu1}/\code{mu2} -- so brms refused the
#' whole model with "The following priors do not correspond to any model
#' parameter: b ~ normal(0, 5)", a prior the user never wrote.  The
#' class-level \code{b} rows are therefore read back from
#' \code{brms::get_prior()}, as for the length-scale prior, and the prior is
#' set on each.  A family with one \code{mu} gets the single dpar-less row it
#' always had.
#'
#' @param spec The prior, as a Stan distribution string.
#' @param fml,data,family As handed to \code{brms::brm()}.
#' @return A \code{brmsprior}; the global \code{class = "b"} row when
#'   \code{brms::get_prior()} fails or reports no class-level \code{b} row.
#' @keywords internal
#' @noRd
.b_prior_rows <- function(spec, fml, data, family) {
  gp_def <- tryCatch(brms::get_prior(fml, data = data, family = family),
                     error = function(e) NULL)
  keep <- if (!is.null(gp_def) && all(c("class", "coef") %in% names(gp_def)))
    gp_def$class == "b" & !nzchar(gp_def$coef) else FALSE
  if (!any(keep)) return(brms::set_prior(spec, class = "b"))
  field <- function(nm) {
    v <- if (nm %in% names(gp_def)) as.character(gp_def[[nm]][keep])
         else rep("", sum(keep))
    v[is.na(v)] <- ""
    v
  }
  rows <- unique(data.frame(resp = field("resp"), dpar = field("dpar"),
                            nlpar = field("nlpar"), stringsAsFactors = FALSE))
  Reduce(`+`, lapply(seq_len(nrow(rows)), function(i)
    brms::set_prior(spec, class = "b", resp = rows$resp[[i]],
                    dpar = rows$dpar[[i]], nlpar = rows$nlpar[[i]])))
}


#' Label lscale coefficients for a log line, with their dpar/nlpar/resp
#' @keywords internal
#' @noRd
.lscale_row_labels <- function(rows) {
  where <- paste(rows$resp, rows$dpar, rows$nlpar, sep = "/")
  where <- gsub("^/+|/+$", "", gsub("/{2,}", "/", where))
  paste0(sQuote(rows$coef), ifelse(nzchar(where), paste0(" (", where, ")"), ""))
}


# ---------------------------------------------------------------------------
# The automatic slope prior
# ---------------------------------------------------------------------------

# brms families whose location parameter is on the scale of log(response)
# although they report an identity link.
.brms_log_location_families <- c("lognormal", "shifted_lognormal",
                                 "hurdle_lognormal")

#' The automatic slope prior for standardised predictors
#'
#' A slope per standard deviation of a predictor is on the scale of the linear
#' predictor.  Under an identity link that is the scale of the response, so a
#' fixed \code{normal(0, 5)} is a different prior for every choice of units:
#' it left 98.6 percent of a slope alone with income in thousands of dollars
#' and 0.007 percent of the same slope with income in dollars.  The prior is
#' therefore \code{normal(0, 2.5 * s)}, where \code{s} is the standard
#' deviation of whatever the linear predictor is a model of:
#' \itemize{
#'   \item \code{sd(y)} for an identity link (gaussian, student, skew_normal
#'     and the like, and mixtures of them);
#'   \item \code{sd(log(y))} for the lognormal families, which report an
#'     identity link while their location is on the log scale;
#'   \item nothing for any other link (log, logit, ...), where the linear
#'     predictor has a scale of its own and the fixed \code{normal(0, 5)}
#'     stays.
#' }
#' A predictor standardised to unit SD cannot move the response by much more
#' than one SD of the response per SD of its own, so 2.5 of them is weak
#' without being flat.  This is the default scaling of \pkg{rstanarm}
#' (Gabry and Goodrich, "Prior Distributions for rstanarm Models").
#'
#' @param y The response on the training rows.
#' @param family As handed to \code{brms::brm()}.
#' @return A list with \code{spec} (the Stan distribution string),
#'   \code{scale} (\code{s}; 1 on the unit scale), \code{basis}
#'   (\code{"sd"}, \code{"sd_log"} or \code{"unit"}) and \code{label} (how to
#'   say it).
#' @keywords internal
#' @noRd
.b_prior_auto <- function(y, family) {
  unit <- list(spec = "normal(0, 5)", scale = 1, basis = "unit",
               label = "fixed, on the scale of the link")
  fam <- tryCatch({
    f <- family
    if (is.function(f)) f <- f()
    if (is.character(f)) f <- brms::brmsfamily(f[[1L]])
    f
  }, error = function(e) NULL)
  if (!is.list(fam)) return(unit)
  parts <- if (is.list(fam$mix) && length(fam$mix) > 0L) fam$mix else list(fam)
  one <- function(p, nm) {
    v <- tryCatch(as.character(p[[nm]])[1L], error = function(e) NA_character_)
    if (length(v) != 1L) NA_character_ else v
  }
  links <- vapply(parts, one, character(1), nm = "link")
  fams  <- vapply(parts, one, character(1), nm = "family")
  if (anyNA(links) || !all(links == "identity")) return(unit)
  if (!(is.numeric(y) || is.logical(y))) return(unit)
  y <- as.numeric(y)
  on_log <- fams %in% .brms_log_location_families
  if (all(on_log)) {
    v <- y[is.finite(y) & y > 0]
    s <- if (length(v) > 1L) stats::sd(log(v)) else NA_real_
    basis <- "sd_log"; label <- "2.5 x SD of the log response"
  } else if (!any(on_log)) {
    v <- y[is.finite(y)]
    s <- if (length(v) > 1L) stats::sd(v) else NA_real_
    basis <- "sd"; label <- "2.5 x SD of the response"
  } else {
    return(unit)
  }
  # A constant response has no scale to give.  normal(0, 0) passes
  # brms::validate_prior() and Stan cannot use it, so fall back to one unit.
  if (!is.finite(s) || s <= 0) {
    s <- 1; basis <- "unit"; label <- "2.5; the response has no spread"
  }
  # %.10g, as for the length-scale prior: fixed-point formatting rounds a
  # small scale to the literal 0.
  list(spec = sprintf("normal(0, %.10g)", 2.5 * s), scale = s, basis = basis,
       label = label)
}


#' The slope prior brms will use, read back from validate_prior()
#'
#' @param vp The table \code{brms::validate_prior()} returned, or \code{NULL}.
#' @param requested The automatic prior's spec, or \code{NULL} when none was
#'   set.
#' @return One string: the class-level \code{b} priors (\code{"(flat)"} for
#'   brms's default), then any coefficient-level ones as \code{name: prior}.
#'   The requested spec, or \code{"(flat)"}, when \code{vp} cannot be read.
#' @keywords internal
#' @noRd
.b_prior_used <- function(vp, requested = NULL) {
  fallback <- if (is.null(requested)) "(flat)" else requested
  if (is.null(vp) || !is.data.frame(vp) ||
      !all(c("class", "coef", "prior") %in% names(vp)))
    return(fallback)
  is_b <- vp$class == "b"
  if (!any(is_b)) return(fallback)
  cls <- unique(vp$prior[is_b & !nzchar(vp$coef)])
  cls[!nzchar(cls)] <- "(flat)"
  if (!length(cls)) cls <- "(flat)"
  own <- is_b & nzchar(vp$coef) & nzchar(vp$prior)
  paste(c(unique(cls),
          if (any(own)) unique(sprintf("%s: %s", vp$coef[own], vp$prior[own]))),
        collapse = "; ")
}


# ---------------------------------------------------------------------------
# What the basis resolves, in map units
# ---------------------------------------------------------------------------

#' The linear unit of a layer's CRS, for a message
#' @keywords internal
#' @noRd
.gp_unit_label <- function(x) {
  u <- tryCatch(sf::st_crs(x)$units_gdal, error = function(e) NULL)
  if (length(u) != 1L || is.na(u) || !nzchar(u)) "CRS units" else u
}

#' A distance in map units, for a message
#' @keywords internal
#' @noRd
.gp_fmt_dist <- function(v) {
  ifelse(is.finite(v),
         format(signif(v, 3), big.mark = ",", scientific = FALSE, trim = TRUE),
         "NA")
}

#' Basis functions per dimension needed to resolve an effective range
#'
#' The inversion of \code{m >= 1.75 * c * S / ell}, with the effective range
#' turned into a length-scale on the axis where it is shortest in scaled
#' units (the one with the larger coordinate SD).
#'
#' @param range Effective range in map units.
#' @param gp_c,S Boundary factor and pooled centred range.
#' @param xy_scale The two coordinate SDs the axes were divided by.
#' @return An integer-valued numeric (not capped).
#' @keywords internal
#' @noRd
.gp_k_needed <- function(range, gp_c, S, xy_scale) {
  ell <- (range / .gp_se_range_factor) / max(xy_scale)
  ceiling(1.75 * gp_c * S / ell)
}


#' Parse the length-scale prior this package sets
#'
#' @param spec A prior string: \code{"inv_gamma(a, b)"} or
#'   \code{"normal(0, s)"} (a half-normal, since \code{lscale} has a lower
#'   bound of 0).  Anything else, or several different entries, is not parsed.
#' @return A list with \code{cdf} (a function) and \code{var}, or \code{NULL}.
#' @keywords internal
#' @noRd
.lscale_prior_parse <- function(spec) {
  if (!is.character(spec) || length(spec) != 1L || is.na(spec)) return(NULL)
  num <- "\\s*([-+0-9.eE]+)\\s*"
  m <- regmatches(spec, regexec(paste0("^inv_gamma\\(", num, ",", num, "\\)$"),
                                spec))[[1L]]
  if (length(m) == 3L) {
    a <- suppressWarnings(as.numeric(m[2L])); b <- suppressWarnings(as.numeric(m[3L]))
    if (!is.finite(a) || !is.finite(b) || a <= 0 || b <= 0) return(NULL)
    return(list(
      cdf = function(x) 1 - stats::pgamma(1 / x, shape = a, rate = b),
      var = if (a > 2) b^2 / ((a - 1)^2 * (a - 2)) else NA_real_))
  }
  m <- regmatches(spec, regexec(paste0("^normal\\(\\s*0\\s*,", num, "\\)$"),
                                spec))[[1L]]
  if (length(m) == 2L) {
    s <- suppressWarnings(as.numeric(m[2L]))
    if (!is.finite(s) || s <= 0) return(NULL)
    return(list(cdf = function(x) 2 * stats::pnorm(x / s) - 1,
                var = s^2 * (1 - 2 / pi)))
  }
  NULL
}


#' Moran's I of a variable's OLS residuals over a weight matrix
#'
#' One-sided (greater), with the Cliff and Ord regression-residual moments.
#'
#' @param y Numeric response.
#' @param X Design matrix with its intercept column.
#' @param W Weight matrix, base or \pkg{Matrix}.
#' @return A list with \code{I}, \code{z} and \code{p}, or \code{NULL} when
#'   it cannot be computed.
#' @keywords internal
#' @noRd
.gp_resid_moran <- function(y, X, W) {
  n  <- length(y)
  e  <- stats::lm.fit(X, y)$residuals
  e  <- e - mean(e)
  S0 <- sum(W)
  if (!is.finite(S0) || S0 <= 0 || sum(e^2) <= 0) return(NULL)
  I   <- as.numeric((n / S0) * sum(e * (W %*% e)) / sum(e^2))
  mom <- .morans_residual_moments(W = W, X = X, S0 = S0,
                                  is_sparse = inherits(W, "Matrix"))
  if (is.null(mom) || !is.finite(mom$VI) || mom$VI <= 0) return(NULL)
  z <- (I - mom$EI) / sqrt(mom$VI)
  list(I = I, z = z, p = stats::pnorm(z, lower.tail = FALSE))
}

#' Short-lag Moran's I of a variable's OLS residuals
#'
#' Over the \code{k} nearest neighbours.  Used before fitting to ask whether
#' the response has any spatial structure at short range before a variogram
#' range below the basis floor is taken at its word.
#' \code{estimate_sac_range()} tests for structure itself and refuses most
#' white-noise layers, and still returned a finite range below the floor in
#' 3 of 60 draws.
#'
#' @param xy Coordinate matrix.
#' @param y Numeric response.
#' @param X Design matrix without the intercept, or \code{NULL}.
#' @param k Neighbours.  4: the short lag is where structure the basis cannot
#'   resolve shows.
#' @return A list with \code{I}, \code{z} and \code{p} (one-sided, greater),
#'   or \code{NULL} when it cannot be computed.
#' @keywords internal
#' @noRd
.gp_short_lag_moran <- function(xy, y, X = NULL, k = 4L) {
  tryCatch({
    y  <- as.numeric(y)
    xy <- as.matrix(xy)[, 1:2, drop = FALSE]
    n  <- length(y)
    if (n < 10L || nrow(xy) != n || !all(is.finite(y)) || !all(is.finite(xy)))
      return(NULL)
    X <- if (is.null(X) || !length(X)) matrix(1, n, 1L) else cbind(1, as.matrix(X))
    .gp_resid_moran(y, X, .build_knn_weights(xy, k = k))
  }, error = function(e) NULL)
}

# The correlation that defines an effective range: the distance at which it
# has fallen to this (.gp_se_range_factor is the same 0.05 for the kernel).
.gp_range_correlation <- 0.05

#' Moran's I of a variable's OLS residuals at the basis floor's distance
#'
#' Over the pairs of points between 0.75 and 1.25 floors apart.  An effective
#' range is the distance at which the correlation has fallen to 0.05, so a
#' response whose range is below the floor is uncorrelated at the floor's
#' distance, and one that is still correlated there has structure the floor
#' does not bound.  A variogram can return a short range for such a response:
#' on squared-exponential fields with ranges of 60 and 90 km it returned 5 to
#' 7 km in 1 to 9 of 60 draws, each with significant short-lag
#' autocorrelation.  Asking that this statistic not be both above 0.05 and
#' significant took the warnings on those fields from 15 of 180 draws to 2,
#' on a sine with a 96 km wavelength from 8 of 60 to 0, and on fields with
#' ranges of 8 to 15 km from 178 of 180 to 176.
#'
#' The band needs every pairwise distance, so more than \code{max_n} rows are
#' thinned to \code{max_n} at even steps through the rows (no random draw:
#' the fit must leave the caller's random-number stream where it was).
#'
#' @param xy Coordinate matrix.
#' @param y Numeric response.
#' @param X Design matrix without the intercept, or \code{NULL}.
#' @param floor The floor, in map units.
#' @param max_n Most rows used.
#' @return A list with \code{I}, \code{z}, \code{p} (one-sided, greater),
#'   \code{band} (the two distances) and \code{n_pairs}, or \code{NULL} when
#'   it cannot be computed (no pair in the band, among other reasons).
#' @keywords internal
#' @noRd
.gp_floor_moran <- function(xy, y, X = NULL, floor, max_n = 1500L) {
  tryCatch({
    y  <- as.numeric(y)
    xy <- as.matrix(xy)[, 1:2, drop = FALSE]
    n  <- length(y)
    if (n < 10L || nrow(xy) != n || !all(is.finite(y)) || !all(is.finite(xy)) ||
        !is.finite(floor) || floor <= 0)
      return(NULL)
    X <- if (is.null(X) || !length(X)) matrix(1, n, 1L) else cbind(1, as.matrix(X))
    if (n > max_n) {
      keep <- unique(round(seq(1, n, length.out = max_n)))
      xy <- xy[keep, , drop = FALSE]; y <- y[keep]; X <- X[keep, , drop = FALSE]
    }
    band <- c(0.75, 1.25) * floor
    D <- as.matrix(stats::dist(xy))
    W <- (D >= band[1L] & D < band[2L]) * 1
    diag(W) <- 0
    out <- .gp_resid_moran(y, X, W)
    if (is.null(out)) return(NULL)
    out$band    <- band
    out$n_pairs <- sum(W) / 2
    out
  }, error = function(e) NULL)
}


#' Compare the response's estimated range with what the basis resolves
#'
#' Run before anything is compiled.  The floor of the basis is known from the
#' coordinates alone, and a range below it means minutes of sampling for a
#' surface that cannot hold the structure.
#'
#' A variogram range below the floor is taken as the reach of the structure
#' only when two things agree with it: the variable is autocorrelated at
#' short range, and it is no longer autocorrelated at the floor's distance.
#' Without the first the range is noise (a variogram fitted to white noise
#' can still return a short range: 3 of 60 draws).  Without the second the
#' range is contradicted: structure reaches at least to the floor, and the
#' variogram either fitted a smooth surface badly or found a short range
#' beside a longer one.
#'
#' @param dat_sf The prepared training layer.
#' @param response_var,predictor_vars Column names.
#' @param resolution Named \code{c(x, y)}: the smallest effective range the
#'   basis resolves along each axis, in map units.
#' @param gp_c,S,xy_scale As in \code{.gp_k_needed()}.
#' @return A list with \code{status} (\code{"ok"}, \code{"below_floor"},
#'   \code{"below_floor_no_autocorrelation"},
#'   \code{"below_floor_long_range"}, \code{"not_estimated"} or
#'   \code{"not_run"}), \code{range} (map units, or \code{NA}), \code{floor}
#'   (the coarser axis's), \code{k_needed}, \code{moran},
#'   \code{moran_floor} and \code{reason}.
#' @keywords internal
#' @noRd
.gp_prefit_range_check <- function(dat_sf, response_var, predictor_vars,
                                   resolution, gp_c, S, xy_scale) {
  out <- list(status = "not_run", range = NA_real_,
              floor = max(resolution), k_needed = NA_real_, moran = NULL,
              moran_floor = NULL, reason = NA_character_)
  df <- sf::st_drop_geometry(dat_sf)
  y  <- df[[response_var]]
  if (!(is.numeric(y) || is.logical(y))) {
    out$reason <- "the response is not numeric"
    return(out)
  }
  # The range of what the GP is asked to model: the response net of the
  # predictors when there are any.  A predictor estimate_sac_range() cannot
  # detrend on (it stops) falls back to the raw response.
  est <- function(pv) tryCatch(
    suppressWarnings(suppressMessages(
      estimate_sac_range(dat_sf, response_var,
                         predictor_vars = if (length(pv)) pv else NULL))),
    error = function(e) structure(NA_real_, rejected_reason = conditionMessage(e),
                                  failed = TRUE))
  r <- est(predictor_vars)
  if (isTRUE(attr(r, "failed")) && length(predictor_vars) > 0L)
    r <- est(character(0))
  rng <- suppressWarnings(as.numeric(r))[1L]

  design <- function() {
    X <- NULL
    if (length(predictor_vars) > 0L)
      X <- tryCatch(
        stats::model.matrix(stats::reformulate(predictor_vars), data = df)[, -1L, drop = FALSE],
        error = function(e) NULL)
    xy <- sf::st_coordinates(dat_sf)[, 1:2, drop = FALSE]
    if (!is.null(X) && nrow(X) != nrow(xy)) X <- NULL
    list(xy = xy, X = X)
  }
  moran <- function() {
    d <- design()
    .gp_short_lag_moran(d$xy, as.numeric(y), d$X)
  }

  if (!is.finite(rng) || rng <= 0) {
    out$status <- "not_estimated"
    out$reason <- attr(r, "rejected_reason") %||% "no range was returned"
    out$moran  <- moran()
    # Both tests, so the log can say which way an unestimated range leans.
    if (!is.null(out$moran) && isTRUE(out$moran$p < 0.05)) {
      d <- design()
      out$moran_floor <- .gp_floor_moran(d$xy, as.numeric(y), d$X,
                                         floor = out$floor)
    }
    return(out)
  }
  out$range    <- rng
  out$k_needed <- .gp_k_needed(rng, gp_c, S, xy_scale)
  if (rng >= out$floor) {
    out$status <- "ok"
    return(out)
  }
  # Below the floor.  A variogram fitted to a response with no spatial
  # structure can still return a short range (estimate_sac_range() refuses
  # most such layers itself; 3 of 60 white-noise draws got one), so the
  # range is taken at its word only when the short-lag autocorrelation is
  # there (or cannot be tested).
  out$moran <- moran()
  if (!(is.null(out$moran) || isTRUE(out$moran$p < 0.05))) {
    out$status <- "below_floor_no_autocorrelation"
    return(out)
  }
  # And only when the autocorrelation is gone by the floor's distance.  Any
  # real structure makes the short-lag test significant, long-range structure
  # included, so that test alone let a spuriously short range through.
  d <- design()
  out$moran_floor <- .gp_floor_moran(d$xy, as.numeric(y), d$X, floor = out$floor)
  mf <- out$moran_floor
  out$status <- if (!is.null(mf) && isTRUE(mf$I > .gp_range_correlation) &&
                    isTRUE(mf$p < 0.05))
    "below_floor_long_range" else "below_floor"
  out
}


# ---------------------------------------------------------------------------
# Sampler diagnostics
# ---------------------------------------------------------------------------

# The thresholds convergence_ok is judged by; see @param check_convergence
# for where each comes from.  R-hat has two: above .bayes_rhat_max fails
# whatever else holds, above .bayes_rhat_caution fails only beside too few
# effective draws (.bayes_rhat_tiers()).
.bayes_rhat_max      <- 1.05
.bayes_rhat_caution  <- 1.01
.bayes_ess_per_chain <- 100
.bayes_ebfmi_min     <- 0.3

#' Sort R-hat values into the ones that fail and the ones to caution on
#'
#' One cut at 1.01 over every parameter fails fits that have converged.  This
#' model has several hundred basis weights, and the largest of several hundred
#' R-hat values is near 1.01 by chance at the default run length.  Four
#' default fits (four chains of 1,000 post-warmup draws, about 540
#' parameters) had maxima of 1.005 to 1.011 with every effective sample size
#' at least 1,086 against 400 required, and one of the four was failed on a
#' single basis weight at 1.011.  Chains that have converged by construction
#' (540 independent AR(1) parameters, four chains of 1,000 draws, eight
#' replicates each) give the same picture: at an effective sample size near
#' 1,300 the largest R-hat was 1.008 to 1.011 and one replicate in eight had
#' a parameter above 1.01; near 450 it was 1.027 to 1.032, with 130 to 160
#' parameters above 1.01 and none above 1.05.
#'
#' The same simulation ties the two recommendations of Vehtari et al. (2021)
#' together: parameters with an effective sample size of 350 to 450 had a
#' mean R-hat of 1.009.  So a value above 1.01 is a failure when the
#' parameter's effective sample size is also short, a caution when it is not,
#' and a failure above 1.05 either way.
#'
#' @param rhat_vals Named R-hat values (\code{NaN} for a constant parameter).
#' @param neff_vals Named ratios from \code{brms::neff_ratio()}, or
#'   \code{NULL}.
#' @param n_draws,n_chains Post-warmup draws over all chains, and chains.
#' @return A list with \code{failed} and \code{caution}, the parameter names
#'   in each.  A parameter above 1.01 whose effective sample size cannot be
#'   read fails.
#' @keywords internal
#' @noRd
.bayes_rhat_tiers <- function(rhat_vals, neff_vals = NULL, n_draws = NA_real_,
                              n_chains = NA_real_) {
  nm    <- names(rhat_vals)
  above <- which(rhat_vals > .bayes_rhat_caution)
  if (!length(above)) return(list(failed = character(0), caution = character(0)))
  ess <- rep(NA_real_, length(rhat_vals))
  counts_known <- is.numeric(n_draws) && length(n_draws) == 1L &&
    is.finite(n_draws) && n_draws > 0 && is.numeric(n_chains) &&
    length(n_chains) == 1L && is.finite(n_chains) && n_chains > 0
  if (counts_known && is.numeric(neff_vals) && !is.null(names(neff_vals)) &&
      !is.null(nm))
    ess <- unname(neff_vals[nm]) * n_draws
  enough <- is.finite(ess) & ess >= .bayes_ess_per_chain * n_chains
  fails  <- above[rhat_vals[above] > .bayes_rhat_max | !enough[above]]
  list(failed  = nm[fails], caution = nm[setdiff(above, fails)])
}

#' Parameters whose bulk or tail effective sample size is too small
#'
#' \code{brms::neff_ratio()} is \code{min(bulk ESS, tail ESS) / draws}, so the
#' smaller of the two effective sample sizes is the ratio times the number of
#' post-warmup draws.  The rule is absolute, 100 per chain, as Stan's own
#' warning is: a share of the draws asks for 100 in all at
#' \code{chains = 2, iter = 1000} and for 2,000 at four chains of 10,000.
#'
#' @param neff_vals Named ratios from \code{brms::neff_ratio()}.
#' @param n_draws Post-warmup draws over all chains.
#' @param n_chains Number of chains.
#' @return The names of the parameters below
#'   \code{.bayes_ess_per_chain * n_chains}; \code{character(0)} when the
#'   draw or chain count is unknown.
#' @keywords internal
#' @noRd
.bayes_low_ess <- function(neff_vals, n_draws, n_chains) {
  if (!is.numeric(n_draws) || length(n_draws) != 1L || !is.finite(n_draws) ||
      !is.numeric(n_chains) || length(n_chains) != 1L || !is.finite(n_chains) ||
      n_draws <= 0 || n_chains <= 0)
    return(character(0))
  names(neff_vals)[which(neff_vals * n_draws < .bayes_ess_per_chain * n_chains)]
}

#' E-BFMI per chain from brms::nuts_params()
#'
#' The energy Bayesian fraction of missing information of Betancourt (2016):
#' the mean squared change in energy between successive transitions over the
#' variance of the energy.
#'
#' @param np The data.frame \code{brms::nuts_params()} returns.
#' @return A named numeric, one per chain; length 0 when the energies are
#'   absent.
#' @keywords internal
#' @noRd
.bayes_ebfmi <- function(np) {
  if (!is.data.frame(np) || !all(c("Parameter", "Value", "Chain") %in% names(np)))
    return(numeric(0))
  en <- np[np$Parameter == "energy__", , drop = FALSE]
  if (!nrow(en)) return(numeric(0))
  if ("Iteration" %in% names(en)) en <- en[order(en$Chain, en$Iteration), ]
  vapply(split(en$Value, en$Chain), function(e) {
    if (length(e) < 3L || !is.finite(stats::var(e)) || stats::var(e) <= 0)
      return(NA_real_)
    (sum(diff(e)^2) / length(e)) / stats::var(e)
  }, numeric(1))
}

#' What to do about each sampler check that failed
#'
#' The remedy differs by check.  More iterations bring R-hat down and the
#' effective sample size up.  They do not remove divergent transitions, which
#' need a smaller step size (a higher \code{adapt_delta}), and a low E-BFMI
#' comes from the shape of the posterior.  A fit whose only failure was two
#' divergent transitions used to be told to run more iterations.
#'
#' @param failed Character: any of \code{"divergent"}, \code{"rhat"},
#'   \code{"ess"} and \code{"ebfmi"}.
#' @param iter,chains,adapt_delta The settings the fit ran with.
#' @return One string, a sentence per remedy, in the order above.
#' @keywords internal
#' @noRd
.bayes_convergence_advice <- function(failed, iter, chains, adapt_delta = NULL) {
  out <- character(0)
  if ("divergent" %in% failed) {
    ad <- suppressWarnings(as.numeric(adapt_delta %||% NA_real_))[1L]
    out <- c(out, sprintf(
      paste0("For the divergent transitions, raise adapt_delta in `control`",
             "%s; more iterations do not remove them."),
      if (is.finite(ad)) sprintf(" (now %.2f)", ad) else ""))
  }
  what <- c(if ("rhat" %in% failed) "bring R-hat down",
            if ("ess" %in% failed) "raise the effective sample size")
  if (length(what) > 0L)
    out <- c(out, sprintf(
      "Run more iterations (iter = %s, chains = %s now) to %s.",
      format(iter), format(chains), paste(what, collapse = " and ")))
  if ("ebfmi" %in% failed)
    out <- c(out, paste0(
      "A low E-BFMI comes from the shape of the posterior, not the run ",
      "length: a longer warmup can help, and otherwise the model or its ",
      "priors need changing."))
  paste(out, collapse = " ")
}


# ---------------------------------------------------------------------------
# GP adequacy
# ---------------------------------------------------------------------------

#' Largest gap between a prior's distribution function and a sample's
#'
#' The Kolmogorov-Smirnov distance between the draws of a parameter and its
#' prior: 0 when the posterior is the prior, 1 when they share no mass.  It
#' sees a posterior that moved as well as one that narrowed, which the ratio
#' of variances does not: a length-scale posterior pushed to the long side of
#' its prior is wider than the prior and has moved almost all of its mass.
#'
#' @param x Numeric draws.
#' @param cdf The prior's distribution function.
#' @return A number in 0 to 1, or \code{NA} with fewer than two finite draws.
#' @keywords internal
#' @noRd
.gp_prior_gap <- function(x, cdf) {
  x <- sort(x[is.finite(x)])
  n <- length(x)
  if (n < 2L) return(NA_real_)
  Fx <- cdf(x)
  if (!all(is.finite(Fx))) return(NA_real_)
  max(seq_len(n) / n - Fx, Fx - (seq_len(n) - 1) / n)
}

# The thresholds of "not identified"; see .gp_adequacy().
.gp_share_floor       <- 0.01
.gp_contraction_floor <- 0.25
.gp_prior_gap_floor   <- 0.30

#' What the posterior says about the Gaussian process, against its basis
#'
#' Two questions the single "share of draws below the floor" could not
#' separate.  Did the data give the Gaussian process anything to carry?  If
#' so, does the posterior want a shorter length-scale than the basis
#' resolves?  The share below the floor cannot answer the first: with a
#' derived basis the prior puts under 1 percent of its mass below the floor,
#' so a posterior equal to the prior gives 1 percent against a 10 percent
#' trigger, and with a small user \code{gp_k} the share restates the prior
#' (0.984 prior, 0.986 posterior).
#'
#' "Not identified" is decided on one of two things, and \code{identified_by}
#' records which.
#' \itemize{
#'   \item With a residual standard deviation on the scale of the Gaussian
#'     process (a gaussian or lognormal family), on the share of the variance
#'     the process carries, \code{sdgp^2 / (sdgp^2 + sigma^2)}: not identified
#'     when its 2.5 percent quantile is under 0.01.  On nine fits that
#'     quantile was 0.000 to 0.003 for the three with nothing to find (white
#'     noise, a 30 km field under a 97 km floor, \code{gp_k = 6}) and 0.17 to
#'     0.94 for the six with a surface.
#'   \item Without one, on the length-scales: not identified when on every
#'     axis the posterior neither narrowed (it kept more than three quarters
#'     of its prior's variance) nor moved (the largest gap between its
#'     distribution function and the prior's is under 0.30).  On the
#'     length-scale draws of the same nine fits the gap was 0.07 to 0.15 for
#'     the three with nothing to find and 0.60 to 0.95 for the five whose
#'     posterior moved; the sixth sat at its prior's median and is identified
#'     by narrowing (0.63 and 0.81 of the variance gone, gaps 0.15 and 0.25).
#'     Draws from the prior itself give a gap under 0.20 in 999 of 1,000
#'     samples of 100.  This rule agrees with the share on eight of the nine
#'     and calls the \code{gp_k = 6} fit identified (one axis lost 0.32 of
#'     its prior's variance), as the rule it replaces did; it has been run
#'     on no fit of a family without a residual standard deviation.
#' }
#' The variance ratio alone is not a test of identification, and was used as
#' one: a surface smoother than the prior expects (a plane, a range longer
#' than the extent) moves the length-scale to the long side of its prior,
#' where the posterior is wider than the prior and the ratio's "contraction"
#' is negative.  Three fits whose Gaussian process carried 0.57 to 0.94 of
#' the variance at the 2.5 percent quantile had contractions of -0.34 to
#' -2.27 and were reported as having found no spatial structure.
#'
#' @param ls_draws Matrix of length-scale draws, one column per length-scale.
#' @param sdgp_draws,sigma_draws Numeric draws, or \code{NULL}.
#' @param gp_ell_min The smallest length-scale the basis resolves.
#' @param prior_spec The length-scale prior string, for the contraction and
#'   the gap.
#' @param xy_scale Coordinate SDs, to report ranges in map units.
#' @return A list; see \code{$info$gp_adequacy} in
#'   \code{?fit_bayesian_spatial_model}.  \code{status} is \code{"ok"},
#'   \code{"below_resolution"} or \code{"not_identified"}.  It is what the
#'   posterior says and no more: the caller reads the residuals next
#'   (\code{.gp_fit_resid_moran()}) and splits both \code{"not_identified"}
#'   and \code{"ok"} on what is left in them.
#' @keywords internal
#' @noRd
.gp_adequacy <- function(ls_draws, sdgp_draws = NULL, sigma_draws = NULL,
                         gp_ell_min, prior_spec = NULL, xy_scale = c(1, 1)) {
  qs <- c(0.025, 0.5, 0.975)
  ls_draws   <- as.matrix(ls_draws)
  frac_below <- mean(as.numeric(ls_draws) < gp_ell_min, na.rm = TRUE)

  pr <- .lscale_prior_parse(prior_spec)
  contraction <- prior_gap <- rep(NA_real_, ncol(ls_draws))
  names(contraction) <- names(prior_gap) <- colnames(ls_draws)
  if (!is.null(pr)) {
    if (is.finite(pr$var) && pr$var > 0)
      contraction[] <- 1 - apply(ls_draws, 2L, stats::var) / pr$var
    prior_gap[] <- apply(ls_draws, 2L, .gp_prior_gap, cdf = pr$cdf)
  }

  share <- NULL
  if (is.numeric(sdgp_draws) && is.numeric(sigma_draws) &&
      length(sdgp_draws) == length(sigma_draws) && length(sdgp_draws) > 0L) {
    sh <- sdgp_draws^2 / (sdgp_draws^2 + sigma_draws^2)
    share <- stats::quantile(sh, qs, names = FALSE, na.rm = TRUE)
    names(share) <- c("q2.5", "q50", "q97.5")
  }
  sdgp <- NULL
  if (is.numeric(sdgp_draws) && length(sdgp_draws) > 0L) {
    sdgp <- stats::quantile(sdgp_draws, qs, names = FALSE, na.rm = TRUE)
    names(sdgp) <- c("q2.5", "q50", "q97.5")
  }

  # The share decides wherever there is one.  The length-scales decide only
  # where there is not, and there a posterior that moved counts as much as
  # one that narrowed.
  identified <- NA; identified_by <- NA_character_
  if (!is.null(share) && is.finite(share[["q2.5"]])) {
    identified    <- share[["q2.5"]] >= .gp_share_floor
    identified_by <- "spatial_share"
  } else if (all(is.finite(prior_gap))) {
    narrowed <- is.finite(contraction) & contraction >= .gp_contraction_floor
    moved    <- prior_gap >= .gp_prior_gap_floor
    identified    <- any(narrowed | moved)
    identified_by <- "length_scale"
  }

  # The posterior effective range per length-scale, in map units.  One
  # length-scale per axis under gp_iso = FALSE; a shared one is reported on
  # both axes, each with its own scaling.
  # A family with several distributional parameters has one GP each, and so
  # more columns than axes; no single range describes that fit.
  eff <- NULL
  if (ncol(ls_draws) %in% c(1L, 2L)) {
    axes <- if (ncol(ls_draws) == 2L) list(x = 1L, y = 2L) else list(x = 1L, y = 1L)
    eff <- t(vapply(seq_along(axes), function(i)
      .gp_se_range_factor * xy_scale[[i]] *
        stats::quantile(ls_draws[, axes[[i]]], qs, names = FALSE, na.rm = TRUE),
      numeric(3)))
    dimnames(eff) <- list(names(axes), c("q2.5", "q50", "q97.5"))
  }

  status <- if (identical(identified, FALSE)) "not_identified"
            else if (is.finite(frac_below) && frac_below > 0.10) "below_resolution"
            else "ok"
  list(status = status, frac_below = frac_below, identified = identified,
       identified_by = identified_by, contraction = contraction,
       prior_gap = prior_gap, spatial_share = share, sdgp = sdgp,
       effective_range = eff)
}

# The test on the residuals that the status is split on: neighbours and level.
.gp_resid_k     <- 4L
.gp_resid_level <- 0.05

#' Short-lag Moran's I of a fitted model's residuals
#'
#' What the posterior of the Gaussian process cannot say: whether the surface
#' it fitted left structure behind.  A process that carries a share of the
#' variance at a range clear of the floor is all \code{.gp_adequacy()} sees
#' in a response with a regional trend plus structure finer than the basis,
#' and in one with the trend alone.  The residuals tell them apart, over the
#' few nearest neighbours and one-sided: ten fits of 400 points on a 100 km
#' square (floor 32 to 33 km) with a component at 8 or 12 km gave p of
#' 5.3e-10 to 0.014, five of them under a posterior that read "ok".  See
#' "What the basis resolves" in \code{?fit_bayesian_spatial_model} for the
#' fits with no such component.
#'
#' The fitted values come from \code{fitted()}, which caches them, so the
#' pass over the posterior is not repeated by \code{summary()} or
#' \code{residuals()}.
#'
#' @param obj A \code{bayesian_fit}.
#' @param fam_name The family's name, or \code{NA}.
#' @param k Neighbours.
#' @return A list with \code{test} (a list with \code{I}, \code{z}, the
#'   one-sided \code{p} and \code{k}, or \code{NULL}) and \code{reason} (why
#'   there is no test, or \code{NA}).
#' @keywords internal
#' @noRd
.gp_fit_resid_moran <- function(obj, fam_name = NA_character_, k = .gp_resid_k) {
  none <- function(why) list(test = NULL, reason = why)
  if (isTRUE(fam_name %in% .brms_category_families))
    return(none("an ordinal or categorical family has no residual per row"))
  res <- tryCatch(suppressWarnings(suppressMessages(stats::residuals(obj))),
                  error = function(e) e)
  if (inherits(res, "error"))
    return(none(sprintf("residuals() failed: %s", conditionMessage(res))))
  if (!is.numeric(res) || length(res) != obj$n || sum(is.finite(res)) < 10L)
    return(none("residuals() did not return a finite number for at least 10 rows"))
  mor <- tryCatch(suppressWarnings(suppressMessages(
    residual_morans_i(obj, alternative = "greater", k = k))),
    error = function(e) NULL)
  p <- suppressWarnings(as.numeric(mor$p_value %||% NA_real_))[1L]
  if (is.null(mor) || !is.finite(p))
    return(none("Moran's I could not be computed on the residuals"))
  list(test = list(I = mor$observed, z = mor$z, p = p, k = as.integer(k)),
       reason = NA_character_)
}


#' Fit a Bayesian spatial regression with a 2D Gaussian Process (via brms)
#'
#' Fits a regression whose residual spatial structure is modelled explicitly, as
#' a Gaussian process over the coordinates, and so kept out of the errors.  Two
#' things follow, and they are the reasons to reach for this backend.  First,
#' every quantity comes with a posterior, so predictions carry intervals
#' instead of point estimates.  Whether those intervals are calibrated is
#' something to measure, not assume: score them with
#' \code{\link{cv_bayes}()}, which reports held-out interval coverage and CRPS.
#' Second, the spatial term is a smooth regional surface at a stated
#' resolution.  The Gaussian process is fitted through a reduced-rank basis,
#' and that basis has a floor: an effective range below which it cannot
#' represent structure, reported in map units as \code{$info$gp_resolution}
#' and by \code{print()}.  With the defaults the floor is 14 to 33 percent of
#' the extent along each axis (about 32 km on a uniform 100 km square, 155 km
#' east-west and 63 km north-south on the census tracts of North Carolina).
#' Structure shorter than the floor is counted as noise: the fitted range comes
#' out too long, the residual standard deviation too large, and the intervals
#' too narrow where the structure is.  See "What the basis resolves" below.
#'
#' The fitted length-scale is an estimate of how far the spatial dependence
#' reaches only when its posterior sits clear of both that floor and its
#' prior.  \code{$info$gp_adequacy} says whether it does, and when it does
#' \code{$info$gp_adequacy$effective_range} holds the posterior effective range
#' per axis in map units (the raw \code{lscale} draws are in per-axis
#' standardised units and are not a distance).
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
#'   baseline for asking how much of the surface is spatial structure and how
#'   much covariate effect (and the natural null model for comparing against a
#'   covariate model with \code{\link{compare_models}}).
#' @param family A model family accepted by \code{brms::brm()}: a stats
#'   family function such as \code{poisson()}, or a brms family object such
#'   as \code{brms::zero_inflated_poisson()}, \code{brms::negbinomial()},
#'   \code{brms::hurdle_poisson()}, \code{brms::bernoulli()} or
#'   \code{brms::Beta()}.  Default \code{NULL}, resolved to
#'   \code{stats::gaussian()}.  The family reaches \code{brms::brm()}
#'   unchanged with the spatial GP term still in the formula.  A factor
#'   response is accepted only under \code{brms::categorical()} or an ordinal
#'   family, and those fits have no single expected value per row, so several
#'   methods cannot use them; see the section on non-Gaussian responses for
#'   what each family supports, and the count example below.
#' @param gp_k Positive integer giving the number of GP basis functions
#'   \emph{per dimension}, or NULL (default) to derive it from the
#'   length-scale/domain ratio.  The fitted model carries
#'   \code{gp_k^2} basis functions, not \code{gp_k} (see Details).  A value
#'   above the derived one makes the basis finer, and the length-scale
#'   prior's lower tail moves down to the length-scale that basis resolves,
#'   because a finer basis under the default prior gains almost nothing.  A
#'   value below the derived one leaves the prior where it was, mostly below
#'   what the basis resolves; that share is stored as
#'   \code{$info$gp_adequacy$prior_mass_below}, logged when it passes 10
#'   percent and flagged by \code{print()}.  To ask for a resolution in map
#'   units, use \code{gp_resolution} and leave this \code{NULL}.
#' @param gp_resolution The shortest effective range the basis should resolve,
#'   as a single positive number in the units of the data's CRS (metres for a
#'   UTM layer), or NULL (default) for the percentile rule described under
#'   Details.  \code{gp_k} is derived from it, and the length-scale prior's
#'   lower tail moves with the basis.  It can only refine: a value the default
#'   basis already resolves leaves the default basis in place.  It is an error
#'   when the basis it needs exceeds the cap of 50 functions per dimension
#'   (2,500 in all); the message gives the basis size needed and the shortest
#'   range the cap reaches on this layout.  Cost grows with
#'   \code{gp_k^2}.  Give \code{gp_resolution} or \code{gp_k}, not both.
#' @param gp_iso Logical; passed to \code{brms::gp(iso = )}.  \code{FALSE}
#'   (the default) fits a separate length-scale per coordinate axis, letting the
#'   model learn any directional structure from the data.  \code{TRUE} fits a
#'   single shared length-scale, which makes the kernel anisotropic in the
#'   original CRS by whatever ratio \code{sd(X)/sd(Y)} happens to take, because
#'   the coordinates are standardised per axis beforehand.  See Details.
#' @param gp_c Positive numeric boundary factor for the approximate GP, or
#'   NULL (default) to derive it alongside \code{gp_k}.  The boundary must be
#'   wide enough to contain the longest plausible correlation range; a value
#'   that is too small truncates the domain and degrades the approximation for
#'   smooth, long-range surfaces.  When you set \code{gp_c} and leave
#'   \code{gp_k = NULL}, \code{gp_k} is derived for \emph{your} boundary: a
#'   wider boundary needs more basis functions to resolve the same
#'   length-scale, so raising \code{gp_c} raises the derived \code{gp_k} with
#'   it, up to the cap of 50 per dimension (a capped value is logged).
#' @param prior Optional brms prior specification. When NULL and
#'   \code{standardize_predictors = TRUE}, a weakly informative prior is set on
#'   the regression coefficients, on the scale of the linear predictor:
#'   \code{normal(0, 2.5 * sd(y))} under an identity link (gaussian, student,
#'   skew-normal and mixtures of them), \code{normal(0, 2.5 * sd(log(y)))} for
#'   the lognormal families, whose location is on the log scale, and the fixed
#'   \code{normal(0, 5)} under any other link (log, logit and so on), where the
#'   linear predictor has a scale of its own.  \code{sd()} is taken over the
#'   training rows.  The result therefore does not depend on the units of the
#'   response: a fixed \code{normal(0, 5)} under an identity link, which is
#'   what earlier versions set, left a slope alone when income was in
#'   thousands of dollars and shrank the same slope to nothing when it was in
#'   dollars.  The prior brms used is stored as \code{$info$b_prior} and shown
#'   by \code{print()}.  When \code{prior} is given, no slope prior is added.
#'   A data-informed GP length-scale prior is always appended
#'   automatically unless \code{prior} already contains an entry with
#'   \code{class = "lscale"}.
#' @param chains Number of MCMC chains. Default 4.
#' @param iter Total iterations per chain. Default 2000.
#' @param warmup Warmup iterations. Default floor(iter/2).
#' @param cores Number of cores for the sampler, one chain per core. Default
#'   \code{getOption("mc.cores", 1L)}, the same convention \pkg{brms} uses
#'   itself, so \code{options(mc.cores = 4)} once per session runs the four
#'   default chains in parallel everywhere. The previous default of
#'   \code{detectCores() - 1} took every core but one on any machine, which
#'   is not what a shared server or a check farm wants.
#' @param seed Integer seed. Default 123.
#' @param backend "auto" (default), "cmdstanr", or "rstan". "auto" uses
#'   \pkg{cmdstanr} only when a CmdStan build is actually available, and
#'   \pkg{rstan} otherwise, which \pkg{brms} always brings. The \pkg{cmdstanr}
#'   package is a thin interface and can be installed without a CmdStan build
#'   (\code{cmdstanr::install_cmdstan()} builds one). An explicit "cmdstanr"
#'   with no usable CmdStan raises an error that says how to install it, so the
#'   failure does not come from inside the sampler.
#' @param control Named list of sampler controls, \emph{merged} over the
#'   package defaults \code{list(adapt_delta = 0.9, max_treedepth = 12)}
#'   instead of replacing them.  Passing \code{list(max_treedepth = 15)}
#'   therefore keeps \code{adapt_delta = 0.9}, which matters, because that is
#'   exactly the setting the divergence warning tells you to raise.
#' @param compute_loo Logical; compute PSIS-LOO. Default TRUE.
#' @param standardize_predictors Logical; center and scale numeric predictors
#'   before fitting. Default FALSE. When TRUE, the scaling parameters are
#'   stored in the return value (\code{$info$predictor_scaling}, a
#'   \code{center} and \code{scale} per predictor) so predictions can be
#'   computed correctly.  The model is then fitted on the standardised
#'   predictors, so \code{coef()} reports a slope per standard deviation of
#'   each predictor and an intercept at the predictor means, not the raw-unit
#'   values \code{stats::lm()} would give; see \code{\link{coef.bayesian_fit}}.
#' @param check_convergence Logical; after fitting, read the sampler's
#'   diagnostics and set \code{$info$convergence_ok}.  The verdict follows
#'   Stan's own diagnostics, with these thresholds:
#'   \itemize{
#'     \item any divergent transition after warmup;
#'     \item bulk or tail effective sample size below 100 per chain for any
#'       parameter (Vehtari et al. 2021; the rule \pkg{rstan} warns on).  It
#'       is an absolute count: a share of the draws, which earlier versions
#'       used, asks for 100 in all at \code{chains = 2, iter = 1000};
#'     \item rank-normalised split R-hat, in two tiers.  Above 1.05, the
#'       value \pkg{rstan}'s own sampler warning uses, for any parameter.
#'       Above 1.01, the value Vehtari et al. (2021) recommend, for a
#'       parameter whose effective sample size is also below the count in
#'       the previous item, or cannot be read;
#'     \item E-BFMI below 0.3 in any chain (Betancourt 2016).
#'   }
#'   An R-hat between 1.01 and 1.05 on a parameter with enough effective
#'   draws is a caution, not a failure: it is stored by name in
#'   \code{$info$convergence_diagnostics$rhat_caution}, logged, and counted
#'   by \code{print()}.  The reason is the number of parameters.  This model
#'   has one weight per basis function, about 540 parameters at the usual
#'   basis, and the largest of that many R-hat values sits near 1.01 in a
#'   fit that has converged: four fits at the default four chains of 2,000
#'   iterations had largest values of 1.005, 1.009, 1.009 and 1.011, with
#'   every effective sample size at least 1,086 against the 400 required.
#'   Chains that have converged by construction (540 independent
#'   parameters, four chains of 1,000 draws) had a largest value of 1.008 to
#'   1.011 at an effective sample size near 1,300, and 130 to 160 values
#'   above 1.01, none above 1.05, at an effective sample size near 450.
#'   Vehtari et al. give 1.01 beside an effective sample size above 400, and
#'   an R-hat of 1.01 is about what an effective sample size of 400 gives
#'   (1.009 on average in the same chains), so the two are read together
#'   here.
#'
#'   Transitions that stop at \code{max_treedepth} are counted and logged but
#'   do not change the verdict: Stan treats them as a cost, not as invalid
#'   draws.  A failed verdict sets \code{convergence_ok} to \code{FALSE},
#'   raises one R warning of class \code{spatialkit_convergence} naming each
#'   failed check with its count and its first five parameters, and
#'   \code{print()} flags it; \code{$info$convergence_diagnostics} holds every
#'   value.  The warning closes with the remedy for the checks that failed: a
#'   higher \code{adapt_delta} for divergent transitions, which more
#'   iterations do not remove, and more iterations for R-hat and effective
#'   sample size.  \code{\link{cv_bayes}()} marks such folds in
#'   \code{fold_metrics$convergence_ok} and \code{\link{compare_models}()}
#'   names such a model.  A short run usually fails, on effective sample
#'   size: two fits of two chains and 300 post-warmup draws each had two or
#'   three parameters below 100 effective draws per chain, and 66 and 76 of
#'   their 537 parameters above R-hat 1.01.
#'
#'   The same step checks the GP basis against the posterior and the
#'   residuals (see "What the basis resolves") and stores the result as
#'   \code{$info$gp_adequacy}.  The residuals need the fitted values, one
#'   pass over the posterior that \code{fitted()} caches for later calls.
#'   That check is about the approximation, not the sampler, so it does not
#'   change \code{convergence_ok}; \code{print()} and \code{summary()} show
#'   its status, and two of its states raise a warning of class
#'   \code{spatialkit_gp_adequacy}.  \code{FALSE} skips
#'   all of these and leaves \code{convergence_ok} \code{NA} (not checked)
#'   and the adequacy status \code{"not_checked"}.
#'   Default TRUE.
#' @param pointize Strategy for non-point geometry coercion.
#' @param boundary Optional polygonal sf/sfc whose CRS the data are projected
#'   to.  It is not a clip: rows outside it are kept.
#' @param .already_prepped Logical (internal). If \code{TRUE}, skip the
#'   \code{prep_model_data()} call because the caller has already projected,
#'   coerced, and filtered the data.  The data must then have plain POINT
#'   geometry (an error is raised otherwise).  Used by the CV internals to
#'   avoid a redundant second pass on every fold.  End users should leave
#'   this at the default \code{FALSE}.
#' @details
#' \strong{GP basis count and boundary factor.}
#' \code{brms::gp()} builds a full tensor grid over its covariates, so a term
#' \code{gp(..x, ..y, k = gp_k)} carries \code{gp_k^2} basis functions: the
#' \code{gp_k} argument is the count \emph{per dimension}, not the total rank.
#' Both \code{gp_k} and \code{gp_c} are therefore chosen from the ratio of the
#' estimated length-scale to the domain extent, following
#' Riutort-Mayol et al. (2023): \code{gp_c} is set large enough to contain the
#' upper length-scale bound, and \code{gp_k} large enough to resolve the lower
#' one.  The derived value depends on the layout, not on \code{n}: 23 on a
#' uniform square, 27 to 38 on the census tracts of Iowa, North Carolina and
#' Texas, and the cap of 50 on California's, whose tracts are clustered.
#'
#' The domain extent used is the one \code{brms::gp(c = )} itself multiplies:
#' the full pooled range of the column-centred coordinates
#' (\code{brms:::choose_L()}, taken over the \strong{unique} coordinate rows,
#' because \code{brms:::.data_gp()} reduces the covariates to unique rows first
#' under the default \code{gr = TRUE}, so repeat visits to one location do not
#' widen the domain), not the per-axis half-range in which
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
#' \strong{What the basis resolves.}
#' The smallest length-scale the chosen basis can resolve is
#' \code{1.75 * gp_c * S / gp_k}, stored in standardised units as
#' \code{$info$gp_ell_min}.  As an effective range (the distance at which the
#' squared-exponential correlation has fallen to 5 percent, 2.45 length-scales)
#' and in map units it is \code{$info$gp_resolution}, one value per axis
#' because each coordinate is divided by its own standard deviation.  The
#' floor is soft: the basis keeps 99 percent of a process's variance at the
#' floor, 69 percent at half of it and 41 percent at a third.  The default
#' length-scale prior is what makes it bind, with 1 percent of its mass below
#' the floor and about one part in a hundred million below half of it.
#'
#' Three checks hold the basis against the data.
#' \enumerate{
#'   \item Before anything is compiled, the autocorrelation range of the
#'     response (of its OLS residuals when there are predictors) is estimated
#'     with \code{\link{estimate_sac_range}()} and compared with the coarser
#'     axis's floor.  That range is a variogram estimate and can be wrong in
#'     both directions, so a range below the floor raises a warning only when
#'     two other things agree with it.  The same variable is autocorrelated
#'     at short range (Moran's I over 4 nearest neighbours, one-sided, at 5
#'     percent): a variogram fitted to a response with no spatial structure
#'     can still return a short range (in 3 of 60 white-noise draws;
#'     \code{estimate_sac_range()} refuses most such layers itself).  And it is
#'     no longer autocorrelated at the floor's distance (Moran's I over the
#'     pairs 0.75 to 1.25 floors apart is not both above 0.05, the
#'     correlation that defines an effective range, and significant at 5
#'     percent): a variogram can also return a short range for a surface
#'     that is smooth over long distances.  The warning has class
#'     \code{spatialkit_gp_resolution}, gives both numbers and the basis size
#'     that would resolve the range, says that it is an estimate which the
#'     check after sampling tests against the residuals of the fitted model,
#'     and says when the basis size is past the
#'     cap or past the number of rows, where kriging is the better tool if
#'     the range is right.  A short range the floor-distance test contradicts
#'     is logged and shown by \code{print()} without a warning: it can still
#'     mean short-range structure beside a longer one, which the basis counts
#'     as noise.  Nothing before sampling separates the two; the residuals
#'     after it do (the third check).  A range that cannot be
#'     estimated is logged with the reason and shown by \code{print()}; the
#'     result is in \code{$info$gp_prefit}.  In 60 simulated draws each
#'     (400 points on a 100 km square, floor 33 km) the warning was raised
#'     for 59 fields with a range of 8 to 15 km; for 0 or 1 fields
#'     that were squared-exponential with a range of 60 or 90 km, a sine or
#'     a plane; for 2 of white noise; and for 5 exponential fields with a
#'     60 km range.  An exponential field has structure at every scale, and
#'     a fit of one with a 50 km range had 62 percent of its length-scale
#'     draws below the floor, so those 5 are not plainly wrong.
#'   \item Also before sampling, the share of the length-scale prior below
#'     the floor is computed (\code{$info$gp_adequacy$prior_mass_below}).  A
#'     derived basis keeps it near 0.01; a small \code{gp_k} does not (0.98
#'     at \code{gp_k = 6} where the derived value is 23), and a share above
#'     0.10 is logged and flagged by \code{print()}.
#'   \item After sampling, under \code{check_convergence = TRUE}, the
#'     posterior and the residuals are read.  The posterior says whether the
#'     Gaussian process carries a share of the variance and where its
#'     length-scale sits against the floor.  It cannot say what the surface
#'     left behind, so the residuals are tested for autocorrelation at short
#'     range in every fit (Moran's I over 4 nearest neighbours, one-sided,
#'     at 5 percent; \code{$info$gp_adequacy$residual_moran}).
#'     \code{$info$gp_adequacy$status} is one of:
#'     \describe{
#'       \item{\code{"ok"}}{The Gaussian process carries a share of the
#'         variance, its length-scale is clear of the floor, and the
#'         residuals show no autocorrelation at short range.  A surface
#'         smoother than the prior expects (a plane, a range longer than the
#'         extent) is \code{"ok"}: its length-scale posterior sits on the
#'         long side of the prior and is wider than it.}
#'       \item{\code{"below_resolution"}}{The Gaussian process carries a
#'         share of the variance and more than 10 percent of the
#'         length-scale draws fall below the floor.  The 10 percent trigger
#'         is this package's operationalisation of the check Riutort-Mayol
#'         et al. recommend, not a figure from the paper.}
#'       \item{\code{"unresolved_short_range"}}{The residuals are
#'         autocorrelated at short range: the signature of structure finer
#'         than the basis, counted as noise.  It is reached two ways, and
#'         \code{$info$gp_adequacy$identified} says which.  The Gaussian
#'         process carries no variance that can be told from zero
#'         (\code{FALSE}): all of the structure is below the floor.  Or it
#'         carries a share and its length-scale is clear of the floor
#'         (\code{TRUE}): it fitted a surface, and finer structure sits
#'         beside it, as with a regional trend plus short-range structure.}
#'       \item{\code{"no_signal_above_floor"}}{The Gaussian process carries
#'         no variance that can be told from zero, and the residuals show
#'         nothing left over.  A response with no spatial structure gives
#'         this, so it is not a failure by itself.}
#'       \item{\code{"not_checked"}}{\code{check_convergence = FALSE}.}
#'     }
#'     The second and third raise a warning of class
#'     \code{spatialkit_gp_adequacy}; the fourth is logged.  When the
#'     residuals cannot be tested (an ordinal or categorical family has no
#'     residual per row; a posterior that cannot be drawn), the status rests
#'     on the posterior alone, a line is logged, the reason is in
#'     \code{$info$gp_adequacy$residual_untested}, and \code{print()} keeps
#'     its flag on a variogram range that was below the floor.
#'
#'     The residual test rests on fits of points on a 100 km square whose
#'     floor is 32 to 33 km.  Sampled fits: ten of a response with a
#'     component below the floor (a field with an 8 or 12 km range, alone or
#'     beside a plane or a field with a 90 km range) had p of 5.3e-10 to
#'     0.014, and five of the ten had a posterior that reads \code{"ok"};
#'     twelve of a response with none (planes, bowls, a sine, fields with
#'     ranges of 45 to 200 km) had p of 0.66 to 0.97.  At the posterior mode
#'     of the same model, computed without sampling, the test was
#'     significant in 0 of 160 fits of fields with ranges of 45 to 200 km
#'     and in 3 of 120 fits of a plane, a bowl and a sine (one noise draw of
#'     the 40, shared by the three, whose own Moran's I has p = 0.003);
#'     and in 16 to 18 of 20 fits of each of four responses with an 8 km
#'     field in them, where the length-scale at the mode was below the floor
#'     in 0 or 1 of the 20.  On a fit with nothing left over the test is
#'     conservative, because the surface takes up part of the noise and
#'     leaves the residuals slightly negatively autocorrelated (median I of
#'     -0.03 to -0.07).  It is not sharp just below the floor: a field with
#'     a 20 km range was significant in 8 of 20 fits, and there it is the
#'     length-scale draws that report it (\code{"below_resolution"}).  An
#'     exponential field has structure at every scale and can be reported
#'     with a range above the floor (8 of 20 fits at 60 km).  The 5 percent
#'     level is the usual one and is not tuned; with many rows a small
#'     leftover autocorrelation is significant, so read the size of I
#'     beside its p.
#'
#'     "No variance
#'     that can be told from zero" is decided on the Gaussian process's share
#'     of the variance, \code{sdgp^2 / (sdgp^2 + sigma^2)}: its 2.5 percent
#'     quantile is under 0.01.  On nine fits that quantile was 0.000 to 0.003
#'     for the three with nothing for the process to find and 0.17 to 0.94
#'     for the six with a surface; the threshold rests on those nine.  It
#'     does not mark off a response whose structure is all below the floor:
#'     on five fits of a field with an 8 km range under a 33 km floor the
#'     quantile was 0.000, 0.000, 0.001, 0.015 and 0.019, so two of the five
#'     count as carrying a share, and it is the residual test that reports
#'     them.  The
#'     share exists for a gaussian or lognormal family.  For any other
#'     family the length-scales decide instead: the process is not identified
#'     when on every axis the posterior kept more than three quarters of its
#'     prior's variance and the largest gap between its distribution function
#'     and the prior's is under 0.30.  That second rule has been checked on
#'     the length-scale draws of the same nine gaussian fits, where it agrees
#'     with the share on eight, and on no fit of another family.
#'     \code{$info$gp_adequacy$identified_by} says which rule was used.  The width of the length-scale posterior does not
#'     decide by itself, because a posterior that moved to longer
#'     length-scales is wider than its prior; and the share of draws below
#'     the floor cannot detect an unidentified process, because a posterior
#'     the data did not move puts 1 percent below the floor, the same as the
#'     prior.
#' }
#' A response whose range is below the floor needs a finer basis and a prior
#' that allows the shorter length-scale, both of which \code{gp_resolution}
#' sets, up to the cap.  Past the cap this model is the wrong tool: resolving
#' a 30 km range on the North Carolina tracts would take 156 functions per
#' dimension, 24,336 in all for 2,672 points.  Kriging has no such floor; see
#' \code{\link{kriging_adequacy}()} for cell means.
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
#' \code{gp_iso} does not affect cost: \code{brms::gp()} builds a
#' tensor grid either way, so the model carries \code{gp_k^2} basis functions
#' regardless.  The stored \code{$info$coord_scaling} list records the scaling
#' strategy, and \code{$info$gp_iso} records which kernel was used.
#'
#' @return A \code{bayesian_fit} object (inherits from \code{spatial_fit}).
#'   Supports \code{predict()}, \code{fitted()}, \code{residuals()},
#'   \code{coef()}, \code{summary()}, and \code{model_metrics()}.
#'   Model-specific metadata lives in \code{$info} (coords: the names of the
#'   scaled coordinate columns handed to \code{brms::gp()}; coord_scaling,
#'   predictor_scaling, gp_k, gp_c, gp_iso, gp_n_basis, gp_ell_min,
#'   gp_resolution: the smallest effective range the basis resolves, as
#'   \code{c(x, y)} in map units; gp_resolution_requested: the
#'   \code{gp_resolution} argument, or \code{NULL};
#'   gp_prefit: the check made before sampling, a list with \code{status}
#'   (\code{"ok"}; \code{"below_floor"}, the one that warns;
#'   \code{"below_floor_no_autocorrelation"} and
#'   \code{"below_floor_long_range"}, a range below the floor that the
#'   short-lag or the floor-distance test did not bear out;
#'   \code{"not_estimated"} or \code{"not_run"}), \code{range} (the
#'   estimated range in map units), \code{floor}, \code{k_needed},
#'   \code{moran} and \code{moran_floor} (the two Moran's I tests, each with
#'   \code{I}, \code{z} and a one-sided \code{p}), \code{reason} and
#'   \code{unit};
#'   gp_adequacy: a list with \code{status} (see Details), \code{frac_below}
#'   (the share of length-scale draws below \code{gp_ell_min}),
#'   \code{identified} and \code{identified_by} (\code{"spatial_share"} or
#'   \code{"length_scale"}, the rule that decided it; see Details),
#'   \code{contraction} (one minus the ratio of posterior to prior variance,
#'   per length-scale, negative when the posterior is the wider of the two)
#'   and \code{prior_gap} (the largest gap between the posterior's and the
#'   prior's distribution functions, per length-scale; both \code{NA} for a
#'   length-scale prior this function did not set), \code{spatial_share} and
#'   \code{sdgp} (2.5,
#'   50 and 97.5 percent quantiles; the share only for a gaussian or lognormal
#'   family), \code{effective_range} (the posterior effective range per axis
#'   in map units, same quantiles), \code{residual_moran} (the test on the
#'   residuals the status rests on: \code{I}, \code{z}, a one-sided \code{p}
#'   and \code{k}, the number of neighbours; absent when it could not be
#'   run, and \code{residual_untested} then says why), \code{resolution},
#'   \code{unit} and \code{prior_mass_below};
#'   b_prior: the slope prior brms used (\code{"(flat)"} for its default,
#'   \code{NA} for a model with no slopes), and b_prior_source: where it came
#'   from; predictor_levels: the levels each factor or character predictor
#'   has in the training rows, which \code{predict()} holds new data to;
#'   gp_S: the pooled centred range \code{brms::gp(c = )} multiplies;
#'   gp_cmeans: the column means brms centred the scaled coordinates on;
#'   gp_xy_range: the training extrema of the scaled coordinates, which
#'   \code{predict()} uses, with gp_S and gp_cmeans, to hold the GP boundary
#'   at its fitted value;
#'   gp_lengthscale_bounds: the \code{c(lower, upper)} the length-scale prior
#'   was calibrated over; gp_lscale_prior: the length-scale prior
#'   \code{brms::validate_prior()} reports the model will \emph{actually} use,
#'   which is not necessarily the one this function requested (several entries,
#'   semicolon-separated, if brms resolved the axes differently); loo, looic,
#'   convergence_ok (\code{TRUE} or \code{FALSE}, and \code{NA} when nothing
#'   was checked, as under \code{check_convergence = FALSE}),
#'   convergence_diagnostics: \code{n_divergent}, \code{n_max_treedepth},
#'   \code{ebfmi} (per chain), \code{max_rhat}, \code{min_neff_ratio},
#'   \code{min_ess} (the smallest bulk or tail effective sample size) and
#'   \code{ess_required} (100 per chain), and \code{rhat_failed} /
#'   \code{neff_failed} / \code{ess_failed}, the parameters that failed each
#'   check by name with their values (R-hat, ESS ratio and ESS; empty when
#'   none failed), which is what makes a failed check actionable,
#'   \code{rhat_caution}, the parameters whose R-hat is between 1.01 and 1.05
#'   with enough effective draws, which do not fail the verdict, plus
#'   \code{gp_lscale_below_resolution} and \code{gp_adequacy}, the share and
#'   status described under Details; and n_dropped: the rows
#'   \code{prep_model_data()} removed for missing or non-finite values or a
#'   bad geometry, so \code{$n} can be read against \code{nrow(data_sf)}).
#'   The raw brmsfit is in \code{$engine}.
#' @family model fitting
#' @section Warnings:
#' Every warning this function raises carries a condition class, so a caller
#' can catch or muffle exactly one of them by name:
#' \describe{
#'   \item{\code{spatialkit_rows_dropped},
#'     \code{spatialkit_constant_predictor}}{From
#'     \code{\link{prep_model_data}()}: incomplete rows were dropped, with the
#'     counts; a predictor takes one value in every row kept.}
#'   \item{\code{spatialkit_gp_basis_size}}{Fewer than one row per ten basis
#'     functions.}
#'   \item{\code{spatialkit_gp_resolution}}{A variogram fitted before
#'     compiling puts the response's range below what the basis resolves, and
#'     the autocorrelation of the response agrees.  An estimate; the check
#'     after sampling tests it against the residuals of the fitted model.}
#'   \item{\code{spatialkit_convergence}}{The sampler's diagnostics fail.}
#'   \item{\code{spatialkit_gp_adequacy}}{The basis is too coarse for the
#'     data: the posterior asks for a shorter length-scale than it resolves
#'     (\code{"below_resolution"}), or the residuals are autocorrelated at
#'     short range (\code{"unresolved_short_range"}).}
#' }
#' A response that takes a single value is an error, raised before anything
#' is compiled.  The convergence warning carries the checks that failed as
#' \code{failed} (\code{"divergent"}, \code{"ebfmi"}, \code{"rhat"},
#' \code{"ess"}).  \code{\link{cv_bayes}()} fits once per fold; it collects
#' the last four classes from its folds and raises each once for the run,
#' naming the folds, and \code{\link{cv_spatial}()} does the same for a
#' \code{fit_fn} that calls this function.
#' @section Non-Gaussian responses:
#' Nothing in this function is Gaussian-specific except its default.  A
#' numeric count, binary (0/1 or logical) or bounded response passes straight
#' through to brms under the family you name; zero-inflated and hurdle
#' counts, negative binomial, Bernoulli, beta, ordinal, categorical and
#' mixture families have all been verified to reach \code{brms::brm()} with
#' the GP term intact, the length-scale prior attached to each
#' distributional parameter's GP.
#'
#' The response check is family-aware.  A factor or character response is
#' accepted only under \code{brms::categorical()} or an ordinal family
#' (\code{cumulative}, \code{sratio}, \code{cratio}, \code{acat}), and refused
#' under every other family before anything is compiled; under gaussian a
#' logical response is refused too.  The case to watch is a two-level factor
#' under \code{brms::bernoulli()}: brms would fit it, but
#' \code{residuals()}, \code{summary()}, \code{model_metrics()} and
#' \code{\link{cv_bayes}()} could not score a factor, so convert it to 0/1
#' first.
#'
#' An ordinal or categorical fit has a probability per response category, not
#' one expected value per row.  \code{predict()} with its default
#' \code{type = "epred"}, \code{fitted()}, \code{residuals()},
#' \code{summary()} and \code{model_metrics()} therefore stop with a message
#' saying so, and \code{\link{cv_bayes}()} refuses the family before fitting
#' anything.  \code{predict(type = "predict", draws = TRUE)} returns the
#' posterior predicted categories, as category indices, for new rows as well
#' as the training ones (the share of draws in each category estimates its
#' probability), and \code{brms::posterior_epred(fit$engine)} the
#' probabilities for the training rows.
#'
#' For the numeric families two things follow.  First, the metrics that come
#' back from \code{\link{model_metrics}()} and the \code{cv_*()} functions are
#' not all meaningful for such a response: RMSE and MAE are, MAPE, SMAPE and R-squared
#' are Gaussian-shaped, and for this backend \code{\link{cv_bayes}()}'s CRPS
#' and interval coverage are the proper scores to read.  See
#' \code{\link{model_metrics}()}, section "Which metrics survive a
#' non-Gaussian response".  Second, \code{\link{estimate_sac_range}()} fits
#' a variogram to the raw response (or OLS residuals), whose variance tracks
#' its mean for a count; the range it reports is then less trustworthy than
#' for a Gaussian response, and its help page says how.
#'
#' One trap.  The response check reads the family's name through
#' \code{brms}'s own accessor; a family object it cannot name is treated as
#' unknown and the check is skipped entirely, without falling back to either
#' rule above.  A malformed \code{family} therefore buys less validation, not
#' more, and a wrong response type will surface as a Stan error, with no
#' message from this function.
#'
#' @section Spatial confounding:
#' A fixed-effect coefficient estimated alongside a spatial random effect is a
#' different quantity from the same coefficient in a non-spatial fit.  When a
#' covariate is itself spatially smooth, the GP term absorbs part of its
#' effect and the two estimates can disagree sharply; when the response is
#' smoother than the covariate, the estimate can shrink toward zero regardless
#' of the true effect (Bolin and Wallin 2025).  Under a correctly specified
#' spatial model this is not a bias but a change of estimand (Zimmerman and
#' Ver Hoef 2022): the spatial coefficient is the effect \emph{net of} whatever
#' the spatial field can explain, and the non-spatial one is not.  Which of the
#' two a user wants depends on the question, so the honest diagnostic is to
#' report both side by side and leave them unadjusted: fit the same formula
#' with \code{stats::lm()} or \code{stats::glm()} and compare.  Compare like
#' with like: under \code{standardize_predictors = TRUE} the coefficients here
#' are per standard deviation of each predictor, so either fit the
#' non-spatial model on the same standardised columns or divide these slopes
#' by \code{$info$predictor_scaling[[name]]$scale} first (see
#' \code{\link{coef.bayesian_fit}}).
#'
#' The literature on remedies is unsettled and this function takes no side.
#' Restricted spatial regression (Hughes and Haran 2013) projects the spatial
#' effect off the covariate space; later work argues it does not deliver what
#' it promises (Hanks et al. 2015; Khan and Calder 2022).  Within this backend,
#' modelling the covariate and the spatial field as explicitly correlated
#' (Marques, Kneib and Klein 2022) is expressible in a brms formula with no
#' new code, and a spectral adjustment (Guan et al. 2023) is the principled
#' route a Hilbert-space basis is well placed to support later.
#'
#' @examples
#' \dontrun{
#' # Not run: fits with Stan, which needs a working C++ toolchain and takes
#' # minutes of MCMC -- both outside what an example may assume.
#' if (requireNamespace("brms", quietly = TRUE)) {
#'   library(sf)
#'   set.seed(1)
#'   n <- 60
#'   dat <- st_as_sf(
#'     data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
#'                elev = rnorm(n)),
#'     coords = c("x", "y"), crs = 32632
#'   )
#'   dat$price <- 10 + 0.01 * (st_coordinates(dat)[, 1] - 5e5) +
#'     2 * dat$elev + rnorm(n)
#'   # Two short chains keep this to a few minutes; expect Stan to warn about
#'   # effective sample size, which is the cost of that.  Use its defaults
#'   # (chains = 4, iter = 2000) for a fit to report.
#'   fit <- fit_bayesian_spatial_model(dat, "price", "elev",
#'                                     chains = 2, iter = 1000,
#'                                     compute_loo = FALSE)
#'   print(summary(fit))
#'   head(predict(fit, newdata = dat))
#'
#'   # A zero-inflated count.  The family is the only thing that changes.
#'   # summary()'s R2 and MAPE assume a Gaussian response, so for a count
#'   # read the coefficients here and score the model with cv_bayes(), whose
#'   # CRPS and interval coverage are defined for any response.
#'   dat$count <- rpois(n, exp(0.5 + 0.8 * dat$elev)) * rbinom(n, 1, 0.7)
#'   fit_zip <- fit_bayesian_spatial_model(dat, "count", "elev",
#'                                         family = brms::zero_inflated_poisson(),
#'                                         chains = 2, iter = 1000,
#'                                         compute_loo = FALSE)
#'   coef(fit_zip)
#' }
#' }
#' @references
#' Betancourt, M. (2016). Diagnosing suboptimal cotangent disintegrations in
#' Hamiltonian Monte Carlo. arXiv:1604.00695.
#'
#' Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and Burkner, P.-C.
#' (2021). Rank-normalization, folding, and localization: an improved R-hat
#' for assessing convergence of MCMC (with discussion). \emph{Bayesian
#' Analysis} \strong{16}, 667--718. \doi{10.1214/20-BA1221}
#'
#' Riutort-Mayol, G., Burkner, P.-C., Andersen, M. R., Solin, A. and
#' Vehtari, A. (2023). Practical Hilbert space approximate Bayesian
#' Gaussian processes for probabilistic programming.
#' \emph{Statistics and Computing} \strong{33}, 17.
#' \doi{10.1007/s11222-022-10167-2}
#'
#' Bolin, D. and Wallin, J. (2025). Spatial self-confounding: smoothness-related
#' estimation bias in spatial regression models. \emph{Biometrika} \strong{113}.
#' \doi{10.1093/biomet/asaf076}
#'
#' Guan, Y., Page, G. L., Reich, B. J., Ventrucci, M. and Yang, S. (2023).
#' Spectral adjustment for spatial confounding. \emph{Biometrika}
#' \strong{110}, 699--719. \doi{10.1093/biomet/asac069}
#'
#' Hanks, E. M., Schliep, E. M., Hooten, M. B. and Hoeting, J. A. (2015).
#' Restricted spatial regression in practice: geostatistical models,
#' confounding, and robustness under model misspecification.
#' \emph{Environmetrics} \strong{26}, 243--254. \doi{10.1002/env.2331}
#'
#' Hughes, J. and Haran, M. (2013). Dimension reduction and alleviation of
#' confounding for spatial generalized linear mixed models. \emph{Journal of
#' the Royal Statistical Society: Series B} \strong{75}, 139--159.
#' \doi{10.1111/j.1467-9868.2012.01041.x}
#'
#' Khan, K. and Calder, C. A. (2022). Restricted spatial regression methods:
#' implications for inference. \emph{Journal of the American Statistical
#' Association} \strong{117}, 482--494. \doi{10.1080/01621459.2020.1788949}
#'
#' Marques, I., Kneib, T. and Klein, N. (2022). Mitigating spatial confounding
#' by explicitly correlating Gaussian random fields. \emph{Environmetrics}
#' \strong{33}, e2727. \doi{10.1002/env.2727}
#'
#' Zimmerman, D. L. and Ver Hoef, J. M. (2022). On deconfounding spatial
#' confounding in linear models. \emph{The American Statistician}
#' \strong{76}, 159--167. \doi{10.1080/00031305.2021.1946149}
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
    cores       = getOption("mc.cores", 1L),
    seed        = 123,
    backend     = c("auto", "cmdstanr", "rstan"),
    control     = list(),
    compute_loo = TRUE,
    standardize_predictors = FALSE,
    check_convergence = TRUE,
    pointize    = "auto",
    boundary    = NULL,
    gp_resolution = NULL,
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

  backend <- .resolve_stan_backend(match.arg(backend))

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
  # Rows prep_model_data() removed on the way in, read from the record it
  # leaves -- which a layer prepared by the caller may still carry.  0 when
  # there is none, and when the record no longer describes this layer: a
  # fold's subset inside cv_*() has neither, and cv_*() reports its own.
  n_dropped <- as.integer(.get_row_record(dat_sf, "dropped")$n %||% 0L)

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
  fam_name <- .brms_family_name(family)
  y_resp   <- dat_cols[[response_var]]
  if (identical(fam_name, "gaussian") && !is.numeric(y_resp))
    stop(sprintf(paste0("fit_bayesian_spatial_model(): response '%s' is %s, ",
                        "not numeric, but the family is gaussian. Convert it ",
                        "to numeric, or pass a `family` that matches it (e.g. ",
                        "brms::bernoulli() for an outcome coded 0/1)."),
                 response_var, class(y_resp)[1L]),
         call. = FALSE)
  # A factor or character response under any other family that is not
  # categorical or ordinal.  brms fits a two-level factor under bernoulli(), so
  # this used to pass, and then residuals() came back all NA, summary(),
  # model_metrics() and compare_models() stopped on "response is factor", and
  # cv_bayes() ran a whole fold of MCMC before aborting in the fold scoring.
  # Refused here, before anything is compiled.  A family whose name cannot be
  # read (NA) is not checked, as documented.
  if (!is.na(fam_name) && !identical(fam_name, "gaussian") &&
      !(is.numeric(y_resp) || is.logical(y_resp)) &&
      !(fam_name %in% .brms_category_families))
    stop(sprintf(paste0("fit_bayesian_spatial_model(): response '%s' is %s, ",
                        "not numeric, and the family is %s. Only ",
                        "brms::categorical() and the ordinal families (%s) ",
                        "take a categorical response here; convert it first, ",
                        "e.g. to 0/1 for brms::bernoulli() with as.integer(",
                        "%s == \"<the level to model>\"). residuals(), ",
                        "summary(), model_metrics() and cv_bayes() all need a ",
                        "numeric response to score the fit."),
                 response_var, class(y_resp)[1L], fam_name,
                 paste(setdiff(.brms_category_families, "categorical"),
                       collapse = ", "),
                 response_var),
         call. = FALSE)

  # A response with one value leaves nothing for the regression or the GP to
  # explain.  It used to reach the compiler: a minute of C++ and then a
  # sampler with sigma running into zero, for a model that could have been
  # refused from the column alone.
  if (length(unique(y_resp)) <= 1L)
    stop(sprintf(paste0("fit_bayesian_spatial_model(): response '%s' takes a ",
                        "single value (%s) in all %d row(s), so there is no ",
                        "variation for the regression or the Gaussian process ",
                        "to explain. Check the column and the filter that ",
                        "produced these rows."),
                 response_var, format(y_resp[[1L]]), nrow(dat_sf)),
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
  # rely on when a user supplies gp_k per fold.  A user's gp_c is handed to the
  # rule, so a derived gp_k is sized for the boundary actually used: k grows
  # with c, and a gp_c raised for a long-range surface (the advice under
  # @param gp_c) used to keep the k derived for the default c, a basis coarser
  # than the lower length-scale bound it is sized to resolve.  gp_c is
  # therefore validated first, before the rule reads it.
  if (!is.null(gp_c) &&
      (!is.numeric(gp_c) || length(gp_c) != 1L || !is.finite(gp_c) || gp_c <= 1))
    stop("fit_bayesian_spatial_model(): `gp_c` must be a single finite number > 1.",
         call. = FALSE)
  gp_spec   <- .gp_basis_spec(cbind(dat_df[["..x"]], dat_df[["..y"]]), ls_bounds,
                              c = gp_c)
  gp_k_auto <- is.null(gp_k)
  k_cap     <- as.integer(floor(sqrt(.gp_max_basis)))
  xy_scale  <- c(x = x_scale, y = y_scale)
  unit_lbl  <- .gp_unit_label(dat_sf)

  # ---- A requested resolution, in map units ----
  # The percentile rule sizes the basis from the layout alone; this sizes it
  # from the shortest effective range the caller wants represented.  It can
  # only refine: a request the derived basis already meets keeps the derived
  # basis.  Past the cap it is an error, not a silent cap, because the caller
  # asked for a range the model then cannot hold.
  gp_resolution_requested <- NULL
  if (!is.null(gp_resolution)) {
    if (!is.numeric(gp_resolution) || length(gp_resolution) != 1L ||
        !is.finite(gp_resolution) || gp_resolution <= 0)
      stop(sprintf(paste0("fit_bayesian_spatial_model(): `gp_resolution` must ",
                          "be a single positive number: the shortest effective ",
                          "range the basis should resolve, in the units of the ",
                          "data's CRS (%s)."), unit_lbl), call. = FALSE)
    if (!gp_k_auto)
      stop(paste0("fit_bayesian_spatial_model(): pass `gp_resolution` or ",
                  "`gp_k`, not both. `gp_resolution` derives `gp_k`."),
           call. = FALSE)
    c_use <- if (is.null(gp_c)) gp_spec$c else gp_c
    k_req <- .gp_k_needed(gp_resolution, c_use, gp_spec$S, xy_scale)
    if (!is.finite(k_req) || k_req > k_cap) {
      finest <- .gp_se_range_factor * (1.75 * c_use * gp_spec$S / k_cap) * xy_scale
      stop(sprintf(paste0("fit_bayesian_spatial_model(): `gp_resolution` = %s ",
                          "%s needs gp_k = %s per dimension (%s basis ",
                          "functions), above the cap of %d (%d). At the cap ",
                          "this layout resolves effective ranges from %s %s ",
                          "along x and %s %s along y. A shorter range is ",
                          "outside what the reduced-rank GP can represent ",
                          "here: use kriging, or model cell means (see ",
                          "kriging_adequacy())."),
                   .gp_fmt_dist(gp_resolution), unit_lbl,
                   format(k_req, big.mark = ",", scientific = FALSE),
                   format(k_req^2, big.mark = ",", scientific = FALSE),
                   k_cap, k_cap^2,
                   .gp_fmt_dist(finest[["x"]]), unit_lbl,
                   .gp_fmt_dist(finest[["y"]]), unit_lbl),
           call. = FALSE)
    }
    gp_resolution_requested <- as.numeric(gp_resolution)
    gp_k <- max(as.integer(k_req), gp_spec$k)
  }
  if (is.null(gp_k))  gp_k <- gp_spec$k
  if (is.null(gp_c))  gp_c <- gp_spec$c

  # Validate before the arithmetic below divides by gp_k.  Previously an
  # invalid user value only surfaced as a brms formula error; now it would also
  # produce a nonsensical resolvable-length-scale and a spurious diagnostic.
  if (!is.numeric(gp_k) || length(gp_k) != 1L || !is.finite(gp_k) || gp_k < 2)
    stop("fit_bayesian_spatial_model(): `gp_k` must be a single finite number >= 2.",
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
  # The same floor as an effective range in map units, per axis: what a user
  # can hold against the distances in their data.  The axes differ because
  # each coordinate was divided by its own SD.
  gp_resolution_xy <- .gp_se_range_factor * gp_ell_min * xy_scale
  .log_info(
    paste0("fit_bayesian_spatial_model(): GP basis k = %d per dimension, ",
           "c = %.2f; TOTAL basis functions = %d (n = %d). ",
           "Smallest resolvable length-scale = %.4f on scaled coords: ",
           "effective ranges from %s %s along x and %s %s along y."),
    gp_k, gp_c, gp_k^2, nrow(dat_sf), gp_ell_min,
    .gp_fmt_dist(gp_resolution_xy[["x"]]), unit_lbl,
    .gp_fmt_dist(gp_resolution_xy[["y"]]), unit_lbl
  )
  if (isTRUE(gp_spec$capped) && isTRUE(gp_k_auto) &&
      is.null(gp_resolution_requested))
    .log_warn(
      paste0("fit_bayesian_spatial_model(): the derived GP basis was capped ",
             "at %d per dimension (%d basis functions), the most this function ",
             "fits. It resolves effective ranges from %s %s along x and %s %s ",
             "along y; shorter structure is counted as noise."),
      gp_k, gp_k^2,
      .gp_fmt_dist(gp_resolution_xy[["x"]]), unit_lbl,
      .gp_fmt_dist(gp_resolution_xy[["y"]]), unit_lbl
    )

  # ---- The length-scale prior follows a basis the caller made finer ----
  # The prior's lower 1% tail sits at the percentile bound, and the derived
  # basis is sized to that bound.  A finer basis under the same prior gains
  # almost nothing, because the prior, not the basis, then holds the
  # length-scale up: on fields with a 15 km range, k = 50 with the percentile
  # prior fitted a 34 km range, the same as an exact GP under that prior, and
  # 17.6 km once the tail was moved to what k = 50 resolves.  So the tail
  # moves down with the basis.  It never moves up: a basis coarser than the
  # bound would pin the prior into a sliver (1.76 to 1.89 at gp_k = 6), and
  # that case is reported below as prior mass under the floor instead.  Only
  # a basis finer than the derived one (an explicit gp_k above it, or a
  # gp_resolution that needed more) moves the tail, so a default fit, and a
  # request the derived basis already meets, have the prior they always had.
  if (gp_k > gp_spec$k && gp_ell_min < ls_bounds[["lower"]]) {
    .log_info(
      paste0("fit_bayesian_spatial_model(): gp_k = %d resolves a shorter ",
             "length-scale (%.4f) than the prior's lower bound (%.4f); the ",
             "prior's lower 1%% tail moves to %.4f so the finer basis can be ",
             "used."),
      gp_k, gp_ell_min, ls_bounds[["lower"]], gp_ell_min
    )
    ls_bounds[["lower"]] <- gp_ell_min
  }

  # Fewer than one row per ten basis functions.  The approximation is still a
  # valid GP prior, but most weights are then set by that prior alone and the
  # sampler pays for all of them.
  if (nrow(dat_sf) < gp_k^2 / 10)
    .warn_classed(
      "spatialkit_gp_basis_size",
      paste0("fit_bayesian_spatial_model(): %d row(s) against %d basis ",
             "functions (gp_k = %d per dimension), fewer than one row per ten ",
             "basis functions. Most basis weights are then set by their prior, ",
             "not the data, and sampling is slow. A smaller `gp_k` is cheaper ",
             "but resolves less (now effective ranges from %s %s along x and ",
             "%s %s along y); with this few rows kriging fits the same kind ",
             "of model exactly."),
      nrow(dat_sf), gp_k^2, gp_k,
      .gp_fmt_dist(gp_resolution_xy[["x"]]), unit_lbl,
      .gp_fmt_dist(gp_resolution_xy[["y"]]), unit_lbl)

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
  # env = globalenv(), not this frame.  A formula carries its environment, and
  # brms stores the formula inside the brmsfit, so the default (the frame this
  # line runs in) made the fit and its $engine capture THIS frame -- which by
  # the end holds the brmsfit itself (`fit`), the data twice and every
  # intermediate.  An environment is serialised once but a list every time it
  # is met, so saveRDS() wrote the whole brmsfit twice.  The formula names
  # columns of the data and brms's own gp(); nothing in it needs this frame.
  fml <- stats::as.formula(sprintf("%s ~ %s + %s", response_var, rhs_terms, gp_term),
                           env = globalenv())

  # Build priors in two independent steps:
  # 1. Regression coefficient priors (only when no user prior and predictors
  #    are standardized).
  # 2. GP length-scale prior (always, unless the user's prior already
  #    includes an "lscale" class entry).

  prior_from_user <- !is.null(prior)
  b_prior_auto    <- NULL
  if (is.null(prior)) {
    prior_parts <- list()
    if (isTRUE(standardize_predictors) && length(predictor_vars) > 0L) {
      # The scale follows the response wherever the linear predictor is on
      # the response's scale (.b_prior_auto()).  A fixed normal(0, 5) is a
      # different prior for every choice of units: with income in dollars it
      # shrank a slope of 10,236 per SD to 0.3, with the same income in
      # thousands it left 98.8% of it, and convergence_ok was TRUE both times.
      b_prior_auto <- .b_prior_auto(y_resp, family)
      # One row per distributional parameter's slopes (.b_prior_rows()): a
      # dpar-less row matches no slope of a categorical or mixture model.
      prior_parts <- c(prior_parts, list(
        .b_prior_rows(b_prior_auto$spec, fml, dat_df, family)
      ))
      .log_info(paste0("fit_bayesian_spatial_model(): slope prior %s on the ",
                       "standardised coefficients (%s)."),
                b_prior_auto$spec, b_prior_auto$label)
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
    # So is each one's dpar/nlpar/resp: a categorical or mixture family repeats
    # the same names once per distributional parameter (see
    # .lscale_coef_rows()), and dropping that address made every such fit fail.
    ls_rows <- .lscale_coef_rows(fml, dat_df, family)

    ls_prior <- if (!is.null(ls_rows)) {
      .lscale_prior_rows(lscale_prior_spec, ls_rows)
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
    # The user's prior must actually REACH Stan.  A global lscale prior
    # (class = "lscale", no coef) is silently discarded by brms because every
    # lscale coefficient already carries an individual default -- the very
    # trap the package's own prior fell into.  Expand any global lscale row to
    # the coefficient level, exactly as is done for the automatic prior above,
    # so that supplying `prior = set_prior("normal(0, 1)", class = "lscale")`
    # does what it says.
    is_global_ls <- prior$class == "lscale" & !nzchar(prior$coef)
    if (any(is_global_ls)) {
      ls_rows <- .lscale_coef_rows(fml, dat_df, family)
      if (!is.null(ls_rows)) {
        # Each global row is expanded only onto the coefficients it addresses
        # -- those with its own resp/dpar/nlpar, which is how brms resolves a
        # global prior itself -- and never onto one the user already gave a
        # coefficient-level prior.  Expanding every row onto every coef name
        # turned a dpar-level prior (dpar = "mumid" on a categorical fit) into
        # rows for every category at once, duplicates brms refuses.
        pcol <- function(nm) {
          v <- if (nm %in% names(prior)) as.character(prior[[nm]])
               else rep("", nrow(prior))
          v[is.na(v)] <- ""
          v
        }
        p_resp <- pcol("resp"); p_dpar <- pcol("dpar"); p_nlpar <- pcol("nlpar")
        row_keys  <- paste(ls_rows$coef, ls_rows$resp, ls_rows$dpar,
                           ls_rows$nlpar, sep = "\r")
        user_keys <- paste(prior$coef, p_resp, p_dpar, p_nlpar, sep = "\r")[
          prior$class == "lscale" & nzchar(prior$coef)]
        expanded <- NULL
        done     <- rep(FALSE, nrow(prior))
        targets  <- character(0)
        for (i in which(is_global_ls)) {
          hit <- ls_rows$resp == p_resp[[i]] & ls_rows$dpar == p_dpar[[i]] &
            ls_rows$nlpar == p_nlpar[[i]] & !(row_keys %in% user_keys)
          if (!any(hit)) next
          tgt  <- ls_rows[hit, , drop = FALSE]
          part <- .lscale_prior_rows(prior$prior[[i]], tgt)
          expanded <- if (is.null(expanded)) part else expanded + part
          done[i]  <- TRUE
          targets  <- c(targets, .lscale_row_labels(tgt))
        }
        if (any(done)) {
          kept  <- prior[!done, , drop = FALSE]
          prior <- if (nrow(kept) > 0L) kept + expanded else expanded
          .log_info(paste0("fit_bayesian_spatial_model(): the supplied global ",
                           "'lscale' prior was attached to %s so that brms uses ",
                           "it (a global lscale prior is otherwise discarded)."),
                    paste(targets, collapse = ", "))
        }
      } else {
        .log_warn(paste0("fit_bayesian_spatial_model(): the supplied prior has a ",
                         "GLOBAL 'lscale' entry, which brms discards because ",
                         "every length-scale coefficient already has an ",
                         "individual prior, and brms reported no coefficient ",
                         "names to attach it to. Check $info$gp_lscale_prior ",
                         "against brms::make_stancode()."))
      }
    }
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

  # The slope prior, read from the same table: what brms will use, not what
  # was asked for.  NA when the model has no slopes.
  b_prior_used <- if (length(predictor_vars) > 0L)
    .b_prior_used(if (is.data.frame(vp)) vp else NULL, b_prior_auto$spec)
    else NA_character_
  b_prior_source <- if (length(predictor_vars) == 0L) NA_character_
    else if (!is.null(b_prior_auto)) sprintf("automatic: %s", b_prior_auto$label)
    else if (prior_from_user) "from `prior`"
    else "brms default"

  # ---- Prior mass below the basis floor ----
  # Known before sampling, and it bounds what the post-fit share can say: a
  # posterior that the data did not move leaves this much below the floor.
  # A derived basis keeps it near 0.01.  A small user gp_k does not: 0.984 at
  # gp_k = 6 on a layout whose derived value is 23, where the post-fit share
  # (0.986) only restated it.
  ls_prior_parsed  <- .lscale_prior_parse(gp_lscale_prior_used)
  prior_mass_below <- if (is.null(ls_prior_parsed)) NA_real_
                      else as.numeric(ls_prior_parsed$cdf(gp_ell_min))
  if (is.finite(prior_mass_below) && prior_mass_below > 0.10)
    .log_warn(
      paste0("fit_bayesian_spatial_model(): %.0f%% of the length-scale prior ",
             "(%s) lies below the smallest length-scale gp_k = %d resolves ",
             "(%.4f on scaled coords; effective ranges from %s %s along x and ",
             "%s %s along y). The basis is too coarse for the prior, so the ",
             "share of posterior draws below that floor will be high whatever ",
             "the data say.%s"),
      100 * prior_mass_below, gp_lscale_prior_used, gp_k, gp_ell_min,
      .gp_fmt_dist(gp_resolution_xy[["x"]]), unit_lbl,
      .gp_fmt_dist(gp_resolution_xy[["y"]]), unit_lbl,
      if (!gp_k_auto && gp_k < gp_spec$k)
        sprintf(" Leave `gp_k` NULL (the derived value is %d) or raise it.",
                gp_spec$k)
      else "")

  # ---- The response's own range against the floor ----
  gp_prefit <- .gp_prefit_range_check(
    dat_sf, response_var, predictor_vars, resolution = gp_resolution_xy,
    gp_c = gp_c, S = gp_spec$S, xy_scale = xy_scale)
  gp_prefit$unit <- unit_lbl
  range_of <- if (length(predictor_vars) > 0L)
    "the response net of the predictors" else "the response"
  if (identical(gp_prefit$status, "below_floor")) {
    # A variogram estimate, made before anything is fitted.  It is said as
    # one: the range can be wrong, and the check after sampling reads the
    # fitted model and its residuals.  The message used to state the range
    # as fact and tell the caller the model could not represent it.  It does
    # not say the later check "supersedes" it either: that check answers on
    # the residuals, and a status of "ok" over residuals nobody had read left
    # a warning that was right withdrawn by a check that had not looked.
    k_need <- gp_prefit$k_needed
    k_txt  <- format(k_need, big.mark = ",", scientific = FALSE)
    k2_txt <- format(k_need^2, big.mark = ",", scientific = FALSE)
    advice <- if (k_need > k_cap)
      sprintf(paste0("Resolving it would take gp_k = %s per dimension (%s ",
                     "basis functions), above the cap of %d; if the fitted ",
                     "model bears the range out, use kriging, or model cell ",
                     "means (see kriging_adequacy())."), k_txt, k2_txt, k_cap)
    else
      sprintf(paste0("Resolving it would take gp_k = %s per dimension (%s ",
                     "basis functions): gp_resolution = %s fits that basis.%s"),
              k_txt, k2_txt, format(signif(gp_prefit$range, 3),
                                    scientific = FALSE, trim = TRUE),
              if (k_need^2 > nrow(dat_sf))
                sprintf(paste0(" That is more basis functions than the %d ",
                               "row(s); if the range is right, kriging (see ",
                               "kriging_adequacy() for cell means) is usually ",
                               "the better tool."),
                        nrow(dat_sf))
              else "")
    mf <- gp_prefit$moran_floor
    .warn_classed(
      "spatialkit_gp_resolution",
      paste0("fit_bayesian_spatial_model(): a variogram fitted before ",
             "sampling puts the autocorrelation range of %s at %s %s, shorter ",
             "than the %s %s this basis resolves (gp_k = %d: effective ranges ",
             "from %s %s along x and %s %s along y). This is an estimate and ",
             "can be wrong%s; the check after sampling ($info$gp_adequacy, ",
             "shown by print()) reads the fitted model and tests its ",
             "residuals for what the surface left behind. If ",
             "the range is right, structure at that range is counted as ",
             "noise: the surface is smoothed over it and the intervals are ",
             "too narrow. %s"),
      range_of, .gp_fmt_dist(gp_prefit$range), unit_lbl,
      .gp_fmt_dist(gp_prefit$floor), unit_lbl, gp_k,
      .gp_fmt_dist(gp_resolution_xy[["x"]]), unit_lbl,
      .gp_fmt_dist(gp_resolution_xy[["y"]]), unit_lbl,
      if (!is.null(mf))
        sprintf(paste0(", though the autocorrelation at the floor's distance ",
                       "agrees with it (Moran's I = %.3f over pairs %s to %s ",
                       "%s apart)"), mf$I, .gp_fmt_dist(mf$band[1L]),
                .gp_fmt_dist(mf$band[2L]), unit_lbl)
      else "",
      advice)
  } else if (identical(gp_prefit$status, "below_floor_long_range")) {
    mf <- gp_prefit$moran_floor
    .log_warn(
      paste0("fit_bayesian_spatial_model(): a variogram fitted before ",
             "sampling puts the autocorrelation range of %s at %s %s, below ",
             "the %s %s this basis resolves, but %s is still autocorrelated ",
             "at that distance (Moran's I = %.3f over pairs %s to %s %s ",
             "apart, one-sided p = %.2g; an effective range ends where the ",
             "correlation falls to 0.05). The variogram's range is therefore ",
             "not the reach of the structure, and no warning is raised on ",
             "it. Either the variogram fitted a smooth surface badly, or ",
             "there is short-range structure beside a longer one, which this ",
             "basis counts as noise; the check after sampling ",
             "($info$gp_adequacy) tells the two apart on the residuals of ",
             "the fitted model."),
      range_of, .gp_fmt_dist(gp_prefit$range), unit_lbl,
      .gp_fmt_dist(gp_prefit$floor), unit_lbl, range_of, mf$I,
      .gp_fmt_dist(mf$band[1L]), .gp_fmt_dist(mf$band[2L]), unit_lbl, mf$p)
  } else if (identical(gp_prefit$status, "not_estimated")) {
    mor <- gp_prefit$moran
    .log_warn(
      paste0("fit_bayesian_spatial_model(): the autocorrelation range of %s ",
             "could not be estimated before fitting (%s), so the basis ",
             "(effective ranges from %s %s along x and %s %s along y) was not ",
             "checked against it.%s"),
      range_of, gp_prefit$reason,
      .gp_fmt_dist(gp_resolution_xy[["x"]]), unit_lbl,
      .gp_fmt_dist(gp_resolution_xy[["y"]]), unit_lbl,
      if (!is.null(mor) && isTRUE(mor$p < 0.05)) {
        # Any real structure is autocorrelated at short range, long-range
        # structure included, so that alone does not point below the floor.
        mf <- gp_prefit$moran_floor
        if (!is.null(mf) && isTRUE(mf$I > .gp_range_correlation) &&
            isTRUE(mf$p < 0.05))
          sprintf(paste0(" It is still autocorrelated at the floor's ",
                         "distance (Moran's I = %.3f over pairs %s to %s %s ",
                         "apart, one-sided p = %.2g), so its structure ",
                         "reaches at least as far as the basis resolves."),
                  mf$I, .gp_fmt_dist(mf$band[1L]), .gp_fmt_dist(mf$band[2L]),
                  unit_lbl, mf$p)
        else
          sprintf(paste0(" It is autocorrelated at short range (Moran's I = ",
                         "%.3f over 4 nearest neighbours, one-sided p = %.2g)",
                         "%s, so a range below that floor is not ruled out."),
                  mor$I, mor$p,
                  if (!is.null(mf))
                    sprintf(" and not at the floor's distance (Moran's I = %.3f)",
                            mf$I)
                  else "")
      } else "")
  } else if (identical(gp_prefit$status, "below_floor_no_autocorrelation")) {
    .log_info(
      paste0("fit_bayesian_spatial_model(): a variogram fitted to %s gave a ",
             "range of %s %s, below the basis floor, but it shows no ",
             "short-range autocorrelation (Moran's I = %.3f over 4 nearest ",
             "neighbours, one-sided p = %.2g), so that range is not read as ",
             "structure."),
      range_of, .gp_fmt_dist(gp_prefit$range), unit_lbl,
      gp_prefit$moran$I, gp_prefit$moran$p)
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
  # NA until a check has actually run.  It was TRUE from the start, so
  # check_convergence = FALSE returned convergence_ok = TRUE over an empty
  # diagnostics list -- a verdict on checks that never happened, for a fit
  # whose max R-hat was 1.28 -- and print() had nothing to caveat.
  #
  # The verdict follows Stan's own diagnostics (see @param check_convergence
  # for the source of each threshold).  It used to cut R-hat at 1.05 and to
  # judge effective sample size as a share of the draws, which equals 100
  # per chain only at the default 4 x 2000: a fit of 2 x 600 with bulk and
  # tail ESS of 130 and 137 passed while rstan raised both of its ESS
  # warnings on it.
  convergence_ok <- NA
  convergence_diagnostics <- list()
  conv_problems <- character(0)
  # Which checks failed, so the advice that closes the warning fits them.
  conv_failed   <- character(0)
  gp_adequacy <- list(status = "not_checked")
  if (isTRUE(check_convergence) && inherits(fit, "brmsfit")) {
    convergence_ok <- TRUE
    # Chains and post-warmup draws, read from the fit; the call's own
    # arguments when the fit cannot say.
    n_chains <- tryCatch(as.numeric(brms::nchains(fit)), error = function(e) NA_real_)
    if (length(n_chains) != 1L || !is.finite(n_chains) || n_chains < 1)
      n_chains <- suppressWarnings(as.numeric(chains))[1L]
    n_draws <- tryCatch(as.numeric(brms::ndraws(fit)), error = function(e) NA_real_)
    if (length(n_draws) != 1L || !is.finite(n_draws) || n_draws < 1)
      n_draws <- suppressWarnings(as.numeric(chains) *
                                    (as.numeric(iter) - as.numeric(warmup)))[1L]

    # Divergent transitions, tree-depth saturation and E-BFMI
    np <- tryCatch(brms::nuts_params(fit), error = function(e) NULL)
    if (!is.null(np)) {
      n_divergent <- sum(np$Value[np$Parameter == "divergent__"], na.rm = TRUE)
      convergence_diagnostics$n_divergent <- n_divergent
      if (n_divergent > 0L) {
        convergence_ok <- FALSE
        conv_failed   <- c(conv_failed, "divergent")
        conv_problems <- c(conv_problems, sprintf(
          "%d divergent transition(s) after warmup", n_divergent))
      }
      # Transitions that stopped at the tree-depth limit.  Stan calls this an
      # efficiency concern, not a validity one, so it is recorded and logged
      # and does not change the verdict.
      td <- np$Value[np$Parameter == "treedepth__"]
      if (length(td) > 0L) {
        max_td <- suppressWarnings(as.numeric(control$max_treedepth %||% 10))[1L]
        n_max_td <- sum(td >= max_td, na.rm = TRUE)
        convergence_diagnostics$n_max_treedepth <- n_max_td
        if (n_max_td > 0L)
          .log_warn(
            paste0("fit_bayesian_spatial_model(): %d transition(s) stopped at ",
                   "max_treedepth = %s. That slows sampling and does not by ",
                   "itself invalidate the draws; raise max_treedepth in ",
                   "`control` if the effective sample size is low."),
            n_max_td, format(max_td))
      }
      ebfmi <- .bayes_ebfmi(np)
      if (length(ebfmi) > 0L) {
        convergence_diagnostics$ebfmi <- ebfmi
        low_e <- which(ebfmi < .bayes_ebfmi_min)
        if (length(low_e) > 0L) {
          convergence_ok <- FALSE
          conv_failed   <- c(conv_failed, "ebfmi")
          conv_problems <- c(conv_problems, sprintf(
            "E-BFMI below %s in chain(s) %s (%s)", format(.bayes_ebfmi_min),
            paste(names(ebfmi)[low_e], collapse = ", "),
            paste(sprintf("%.2f", ebfmi[low_e]), collapse = ", ")))
        }
      }
    }
    
    # R-hat and effective sample size.  Both accessors run under
    # .muffle_ess_cap(): posterior's "The ESS has been capped" arrives once
    # per well-mixed basis weight.  Both are read before either is judged,
    # because an R-hat between 1.01 and 1.05 is judged beside the effective
    # sample size of the same parameter.
    rhat_vals <- tryCatch({
      rh <- .muffle_ess_cap(brms::rhat(fit))
      if (is.numeric(rh)) rh else NULL
    }, error = function(e) NULL)
    neff_vals <- tryCatch({
      ne <- .muffle_ess_cap(brms::neff_ratio(fit))
      if (is.numeric(ne)) ne else NULL
    }, error = function(e) NULL)

    if (!is.null(rhat_vals)) {
      if (is.null(names(rhat_vals)))
        names(rhat_vals) <- sprintf("parameter %d", seq_along(rhat_vals))
      # brms::rhat() returns NaN for a parameter that is constant across
      # draws, which `lprior` frequently is.  Two consequences had to be
      # guarded.  A logical subscript keeps one NA element per NaN entry --
      # so convergence_ok was set FALSE and the warning named "NA" -- while
      # which() drops them.  And max(..., na.rm = TRUE) over an all-NA
      # vector is -Inf with a warning, not a diagnostic.
      max_rhat <- if (any(is.finite(rhat_vals)))
        max(rhat_vals[is.finite(rhat_vals)]) else NA_real_
      convergence_diagnostics$max_rhat <- max_rhat
      # Two tiers (.bayes_rhat_tiers()): above 1.05 fails; above 1.01 fails
      # when that parameter's effective sample size is short as well, and is
      # a caution when it is not.  One cut at 1.01 over about 540 parameters
      # failed a default fit whose largest R-hat was 1.011 and whose smallest
      # effective sample size was 1,278 against 400 required.
      tiers    <- .bayes_rhat_tiers(rhat_vals, neff_vals, n_draws, n_chains)
      bad_rhat <- tiers$failed
      # Every failing parameter by name, with its value -- the message
      # names five, and "max R-hat 1.09" is not actionable where
      # "sdgp_gp..x..y has R-hat 1.09" is.  An empty named numeric when none.
      convergence_diagnostics$rhat_failed  <- rhat_vals[bad_rhat]
      convergence_diagnostics$rhat_caution <- rhat_vals[tiers$caution]
      n_caution <- length(tiers$caution)
      if (length(bad_rhat) > 0L) {
        convergence_ok <- FALSE
        conv_failed <- c(conv_failed, "rhat")
        n_hi <- sum(rhat_vals[bad_rhat] > .bayes_rhat_max)
        conv_problems <- c(conv_problems, sprintf(
          paste0("%d parameter(s) fail on R-hat (max %.3f): %d above %s, %d ",
                 "above %s with %s%s: %s"),
          length(bad_rhat), max_rhat, n_hi, format(.bayes_rhat_max),
          length(bad_rhat) - n_hi, format(.bayes_rhat_caution),
          if (is.null(neff_vals))
            "an effective sample size that could not be read"
          else "too few effective draws",
          if (n_caution > 0L)
            sprintf(" (%d more above %s have enough and are not counted)",
                    n_caution, format(.bayes_rhat_caution)) else "",
          paste(utils::head(bad_rhat, 5), collapse = ", ")))
      } else if (n_caution > 0L) {
        .log_warn(
          paste0("fit_bayesian_spatial_model(): %d of %d parameter(s) have ",
                 "R-hat between %s and %s (max %.3f), each with at least the ",
                 "%s effective draws required. Among this many parameters ",
                 "that happens by chance in a converged fit, so ",
                 "$info$convergence_ok does not fail on it; ",
                 "$info$convergence_diagnostics$rhat_caution names them."),
          n_caution, sum(is.finite(rhat_vals)), format(.bayes_rhat_caution),
          format(.bayes_rhat_max), max_rhat,
          format(.bayes_ess_per_chain * n_chains))
      }
    }

    # Effective sample size
    if (!is.null(neff_vals)) {
      # Same NaN hazard as R-hat above: a constant parameter has no effective
      # sample size.  which() drops the NAs a logical subscript would keep as
      # NA names, and min() is guarded so an all-NaN vector reports NA rather
      # than +Inf.
      min_neff <- if (any(is.finite(neff_vals)))
        min(neff_vals[is.finite(neff_vals)]) else NA_real_
      convergence_diagnostics$min_neff_ratio <- min_neff
      # neff_ratio is min(bulk, tail) ESS over the draws, so the smaller of
      # the two effective sample sizes is the ratio times the draws.
      convergence_diagnostics$min_ess     <- min_neff * n_draws
      convergence_diagnostics$ess_required <- .bayes_ess_per_chain * n_chains
      low_neff <- .bayes_low_ess(neff_vals, n_draws, n_chains)
      convergence_diagnostics$neff_failed <- neff_vals[low_neff]
      convergence_diagnostics$ess_failed  <- neff_vals[low_neff] * n_draws
      if (length(low_neff) > 0L) {
        convergence_ok <- FALSE
        conv_failed   <- c(conv_failed, "ess")
        conv_problems <- c(conv_problems, sprintf(
          paste0("%d parameter(s) with bulk or tail effective sample size ",
                 "under %s (%s per chain x %s chain(s); min %.0f): %s"),
          length(low_neff), format(.bayes_ess_per_chain * n_chains),
          format(.bayes_ess_per_chain), format(n_chains), min_neff * n_draws,
          paste(utils::head(low_neff, 5), collapse = ", ")))
      }
    }

    # ---- GP basis adequacy ----
    # gp_k is derived from a PRIOR length-scale bound computed from inter-point
    # spacing.  If the posterior length-scale lands below what (gp_k, gp_c) can
    # resolve, the Hilbert-space approximation is inadequate and nothing else
    # here would say so.  Riutort-Mayol et al. (2023) recommend checking that
    # the posterior length-scale lies inside the region the basis resolves and
    # raising m when it does not; the 10% posterior-mass trigger is this
    # package's operationalisation of that check, not a figure from the paper.
    #
    # That share alone is blind when the data do not inform the length-scale:
    # the posterior then equals the prior, which puts about 1% below the floor
    # by construction.  .gp_adequacy() therefore asks first whether the GP is
    # identified at all; the residual check that splits "no structure" from
    # "structure the basis cannot hold" needs the fitted values and runs once
    # the fit object exists, below.
    #
    # brms names GP parameters lscale(gp...) / sdgp(gp...), embedding the
    # covariate names, so match on the prefix rather than a literal term name.
    # The tryCatch keeps a future brms naming change from erroring the fit.
    ls_draws <- tryCatch(
      as.matrix(fit, variable = "^lscale", regex = TRUE),
      error = function(e) NULL
    )
    if (!is.null(ls_draws) && length(ls_draws) > 0L) {
      one_col <- function(pattern, regex) {
        m <- tryCatch(as.matrix(fit, variable = pattern, regex = regex),
                      error = function(e) NULL)
        if (is.matrix(m) && ncol(m) == 1L) as.numeric(m[, 1L]) else NULL
      }
      sdgp_draws  <- one_col("^sdgp", TRUE)
      # The share of variance the GP carries needs a residual SD on the same
      # scale as the GP: a Gaussian response, or one Gaussian on the log scale.
      sigma_draws <- if (fam_name %in% c("gaussian", "lognormal"))
        one_col("sigma", FALSE) else NULL
      gp_adequacy <- .gp_adequacy(ls_draws, sdgp_draws, sigma_draws,
                                  gp_ell_min = gp_ell_min,
                                  prior_spec = gp_lscale_prior_used,
                                  xy_scale   = xy_scale)
      convergence_diagnostics$gp_lscale_below_resolution <- gp_adequacy$frac_below
    }
    # Every accessor failed, so nothing was checked after all.
    if (!length(convergence_diagnostics)) convergence_ok <- NA
  }
  gp_adequacy$resolution       <- gp_resolution_xy
  gp_adequacy$unit             <- unit_lbl
  gp_adequacy$prior_mass_below <- prior_mass_below

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
      # With the cause: "LOO computation failed." alone left compare_models()
      # showing LOOIC NA with nothing anywhere saying why.
      .log_warn(paste0("fit_bayesian_spatial_model(): LOO computation failed, ",
                       "so $info$looic is NA. Cause: %s"),
                .try_error_message(loo_try))
    }
  }

  obj <- new_spatial_fit(
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
      # The same floor as an effective range in map units, per axis, and the
      # request that set it when there was one.
      gp_resolution            = gp_resolution_xy,
      gp_resolution_requested  = gp_resolution_requested,
      gp_prefit                = gp_prefit,
      gp_adequacy              = gp_adequacy,
      # The slope prior brms used and where it came from.
      b_prior                  = b_prior_used,
      b_prior_source           = b_prior_source,
      # The levels each categorical predictor has in the training rows, which
      # are the levels brms kept; predict() and fitted() hand brms these and
      # no others.
      predictor_levels         = .bayes_predictor_levels(dat_cols, predictor_vars),
      # The scaled-coordinate extrema of the training data, the column means
      # brms centred them on, and the pooled centred range it derived the
      # boundary L = c * S from.  predict() needs them: brms 2.17 to 2.22
      # store only Xgp/dmax/cmeans in the fit's GP basis, NOT the
      # Hilbert-space boundary L, so brms:::.data_gp() RECOMPUTES
      # L = c * range(centred newdata) from whatever rows predict() is handed.
      # Unless L is held, a point's prediction depends on which other points
      # share the predict() call -- predict_surface()'s chunk_size changed the
      # surface, and every cv_bayes() fold was scored against a basis the model
      # was not fitted with.  brms >= 2.23.0 stores L and reuses it.  See
      # .pin_gp_boundary_rows().
      gp_xy_range              = list(
        x = range(dat_df[["..x"]], na.rm = TRUE),
        y = range(dat_df[["..y"]], na.rm = TRUE)
      ),
      gp_S                     = gp_spec$S,
      gp_cmeans                = gp_spec$cmeans,
      gp_lengthscale_bounds    = ls_bounds,
      convergence_ok           = convergence_ok,
      n_dropped                = n_dropped,
      convergence_diagnostics  = convergence_diagnostics,
      predictor_scaling        = predictor_scaling
    )
  )

  # ---- What the surface left behind ----
  # The posterior of the GP cannot say whether structure was left in the
  # residuals, and that is the question for both of its verdicts.
  #
  # A GP the data did not identify looks the same in its hyperparameters
  # whether the response is white noise or has structure finer than the
  # basis: sdgp near zero and a length-scale posterior equal to its prior.
  # What separates them is what is left in the residuals at short range, over
  # the few nearest neighbours and one-sided: on a field with a 30 km range
  # under a 97 to 251 km floor the residual Moran's I had p = 0.003 at 4
  # neighbours and 0.065 at the default 8 two-sided, and white noise on the
  # same layout 0.36.
  #
  # A GP the data did identify looks the same whether the surface it fitted
  # is all there is or sits beside structure finer than the basis.  A weak
  # plane plus a field with an 8 km range under a 33 km floor gave a share of
  # the variance of 0.09 to 0.54 at ranges of 55 and 72 km, 2 percent of the
  # length-scale draws below the floor, and so "ok", over residuals with a
  # Moran's I of 0.143 (p = 5.6e-06) and a fitted surface that explained 19
  # percent of the response.  The residuals were read only for an
  # unidentified GP, so nothing was said.  They are read for every fit whose
  # posterior was, and "ok" is not returned over autocorrelated ones.
  if (!identical(gp_adequacy$status, "not_checked")) {
    rm_res <- .gp_fit_resid_moran(obj, fam_name)
    gp_adequacy$residual_moran <- rm_res$test
    if (is.null(rm_res$test)) gp_adequacy$residual_untested <- rm_res$reason
    left_over <- !is.null(rm_res$test) && rm_res$test$p < .gp_resid_level
    gp_adequacy$status <- switch(
      gp_adequacy$status,
      not_identified = if (left_over) "unresolved_short_range"
                       else "no_signal_above_floor",
      ok             = if (left_over) "unresolved_short_range" else "ok",
      gp_adequacy$status)
    obj$info$gp_adequacy <- gp_adequacy
  }
  # The adequacy status rides with the other diagnostics when any were read;
  # an unchecked fit keeps its empty list.
  if (length(obj$info$convergence_diagnostics) > 0L)
    obj$info$convergence_diagnostics$gp_adequacy <- gp_adequacy$status

  # ---- Say what was found ----
  # One classed R warning for a failed sampler verdict.  The checks used to
  # be log lines only, which spatialkit_quiet(), knitr and tryCatch() hide,
  # and under cmdstanr nothing else warns.
  if (identical(convergence_ok, FALSE) && length(conv_problems) > 0L)
    .warn_classed(
      "spatialkit_convergence",
      paste0("fit_bayesian_spatial_model(): the sampler's diagnostics fail, ",
             "so $info$convergence_ok is FALSE: %s. %s ",
             "$info$convergence_diagnostics has every value."),
      paste(conv_problems, collapse = "; "),
      .bayes_convergence_advice(conv_failed, iter = iter, chains = chains,
                                adapt_delta = control$adapt_delta),
      # Which checks failed ("divergent", "ebfmi", "rhat", "ess"), for a
      # caller that fits fold by fold and reports them by fold.
      data = list(failed = conv_failed))

  floor_txt <- sprintf("%s %s along x and %s %s along y",
                       .gp_fmt_dist(gp_resolution_xy[["x"]]), unit_lbl,
                       .gp_fmt_dist(gp_resolution_xy[["y"]]), unit_lbl)
  # What "not identified" rested on, as the clause that opens both messages:
  # the share of the variance where there is one, the length-scales where
  # there is not.
  id_txt <- function() {
    sh <- gp_adequacy$spatial_share
    if (identical(gp_adequacy$identified_by, "spatial_share") && !is.null(sh))
      return(sprintf(
        paste0("the Gaussian process carries no variance that can be told ",
               "from zero (its share of the variance has a 95%% interval of ",
               "%.3f to %.3f)"), sh[["q2.5"]], sh[["q97.5"]]))
    ct <- gp_adequacy$contraction
    pg <- gp_adequacy$prior_gap
    sprintf(
      paste0("the data did not move the Gaussian process's length-scale from ",
             "its prior (largest gap between the two distribution functions ",
             "%s%s)"),
      paste(sprintf("%.2f", pg), collapse = " and "),
      if (length(ct) && all(is.finite(ct)))
        sprintf("; %s of the prior's variance kept",
                paste(sprintf("%.0f%%", 100 * (1 - ct)), collapse = " and "))
      else "")
  }
  # What can still be done about a basis found short: a finer one if the cap
  # allows, and otherwise a different model.
  finest_cap <- .gp_se_range_factor * (1.75 * gp_c * gp_spec$S / k_cap) * xy_scale
  finer_txt <- if (!gp_k_auto && gp_k < gp_spec$k)
    sprintf(paste0("`gp_k` = %d is below the value derived for this layout ",
                   "(%d): leave `gp_k` NULL, or raise it."), gp_k, gp_spec$k)
  else if (gp_k >= k_cap)
    sprintf(paste0("The basis is at the cap of %d per dimension, so this ",
                   "model cannot hold a shorter range: use kriging, or model ",
                   "cell means (see kriging_adequacy())."), k_cap)
  else sprintf(
    paste0("A `gp_resolution` below %s %s gives a finer basis; the cap of %d ",
           "per dimension reaches %s %s along x and %s %s along y.%s"),
    .gp_fmt_dist(max(gp_resolution_xy)), unit_lbl, k_cap,
    .gp_fmt_dist(finest_cap[["x"]]), unit_lbl,
    .gp_fmt_dist(finest_cap[["y"]]), unit_lbl,
    # The variogram's short range, whichever way the check before fitting
    # read it: this text is used only where the fitted model agrees with it.
    if (gp_prefit$status %in% c("below_floor", "below_floor_long_range") &&
        is.finite(gp_prefit$range) && gp_prefit$range < max(finest_cap))
      sprintf(paste0(" The range a variogram gave before fitting (%s %s) is ",
                     "below that; if it is right, kriging, or a model of cell ",
                     "means (see kriging_adequacy()), suits these data better."),
              .gp_fmt_dist(gp_prefit$range), unit_lbl)
    else "")
  rm_test <- gp_adequacy$residual_moran
  # What an identified GP found, as the clause that opens its message: the
  # share of the variance where there is one, and the range it sits at.
  found_txt <- function() {
    sh  <- gp_adequacy$spatial_share
    eff <- gp_adequacy$effective_range
    at  <- if (is.matrix(eff) && all(is.finite(eff[, "q50"])))
      sprintf(" at an effective range of %s %s along x and %s %s along y",
              .gp_fmt_dist(eff[["x", "q50"]]), unit_lbl,
              .gp_fmt_dist(eff[["y", "q50"]]), unit_lbl)
      else ""
    if (!is.null(sh) && all(is.finite(sh)))
      sprintf(paste0("the Gaussian process carries %.2f to %.2f of the ",
                     "variance (95%% interval)%s"), sh[["q2.5"]], sh[["q97.5"]], at)
    else sprintf("the Gaussian process fitted a surface%s", at)
  }
  if (identical(gp_adequacy$status, "below_resolution")) {
    .warn_classed(
      "spatialkit_gp_adequacy",
      paste0("fit_bayesian_spatial_model(): %.0f%% of the posterior ",
             "length-scale draws fall below the smallest length-scale this ",
             "basis resolves (gp_k = %d, %d basis functions: effective ranges ",
             "from %s). The data ask for a shorter range than the basis can ",
             "hold, so the fitted range is too long and part of the structure ",
             "is counted as noise.%s %s"),
      100 * gp_adequacy$frac_below, gp_k, gp_k^2, floor_txt,
      if (!is.null(rm_test) && rm_test$p < .gp_resid_level)
        sprintf(paste0(" The residuals are autocorrelated at short range ",
                       "(Moran's I = %.3f over %d nearest neighbours, ",
                       "one-sided p = %.2g)."), rm_test$I, rm_test$k, rm_test$p)
      else "",
      finer_txt)
  } else if (identical(gp_adequacy$status, "unresolved_short_range") &&
             identical(gp_adequacy$identified, FALSE)) {
    .warn_classed(
      "spatialkit_gp_adequacy",
      paste0("fit_bayesian_spatial_model(): %s, but the residuals are ",
             "autocorrelated at short range (Moran's I = %.3f over %d nearest ",
             "neighbours, one-sided p = %.2g). Structure with an effective ",
             "range under %s cannot be held by this basis (gp_k = %d) and is ",
             "counted as noise, so the surface misses it and the intervals do ",
             "not cover it. %s"),
      id_txt(), rm_test$I, rm_test$k, rm_test$p, floor_txt, gp_k, finer_txt)
  } else if (identical(gp_adequacy$status, "unresolved_short_range")) {
    # An identified GP over residuals that are still autocorrelated: a
    # surface was found, and structure finer than the basis sits beside it.
    .warn_classed(
      "spatialkit_gp_adequacy",
      paste0("fit_bayesian_spatial_model(): %s, and the residuals are still ",
             "autocorrelated at short range (Moran's I = %.3f over %d nearest ",
             "neighbours, one-sided p = %.2g). The surface left structure ",
             "behind. Structure with an effective range under %s cannot be ",
             "held by this basis (gp_k = %d) and is counted as noise, so the ",
             "surface misses it and the intervals do not cover it. %s"),
      found_txt(), rm_test$I, rm_test$k, rm_test$p, floor_txt, gp_k, finer_txt)
  } else if (identical(gp_adequacy$status, "ok") && is.null(rm_test)) {
    # "ok" on the posterior alone.  Said, because the residuals are what a
    # surface beside finer structure is caught by.
    .log_warn(
      paste0("fit_bayesian_spatial_model(): the residuals could not be tested ",
             "for autocorrelation at short range (%s), so the adequacy status ",
             "\"ok\" rests on the posterior alone: %s, with %.0f%% of the ",
             "length-scale draws below what the basis resolves. Structure ",
             "finer than the basis beside that surface would not show there."),
      gp_adequacy$residual_untested %||% "no reason recorded", found_txt(),
      100 * gp_adequacy$frac_below)
  } else if (identical(gp_adequacy$status, "no_signal_above_floor")) {
    # Not a failure by itself: a response with no spatial structure gives the
    # same picture, and the residuals show nothing left over.
    .log_warn(
      paste0("fit_bayesian_spatial_model(): %s, and the residuals %s. A ",
             "response with no spatial structure gives this. Structure with ",
             "an effective range under %s could not be held by this basis in ",
             "any case and would be counted as noise."),
      id_txt(),
      if (!is.null(rm_test))
        sprintf(paste0("show no autocorrelation at short range (Moran's I = ",
                       "%.3f over %d nearest neighbours, one-sided p = %.2g)"),
                rm_test$I, rm_test$k, rm_test$p)
      else sprintf("could not be tested for autocorrelation (%s)",
                   gp_adequacy$residual_untested %||% "no reason recorded"),
      floor_txt)
  }

  obj
}


#' The levels each categorical predictor takes in the training rows
#'
#' brms drops unused levels when it builds the model frame, and at predict
#' time refuses any \code{newdata} factor whose \emph{level set} is not a
#' subset of the ones it kept, whatever values the rows carry.  These are
#' those levels: the values that occur, in the factor's own order.
#'
#' @param df The training data, without geometry.
#' @param predictor_vars Predictor names.
#' @return A named list, one element per factor or character predictor, each
#'   a list with \code{levels} and \code{ordered}; empty when there is none.
#' @keywords internal
#' @noRd
.bayes_predictor_levels <- function(df, predictor_vars) {
  out <- list()
  for (cn in intersect(predictor_vars, names(df))) {
    v <- df[[cn]]
    if (!(is.factor(v) || is.character(v))) next
    f <- droplevels(as.factor(v))
    out[[cn]] <- list(levels = levels(f), ordered = is.ordered(v))
  }
  out
}


#' Is a CmdStan build available to cmdstanr?
#'
#' The \pkg{cmdstanr} package is an interface; the CmdStan toolchain it drives
#' is installed separately (\code{cmdstanr::install_cmdstan()}) and is absent
#' on many machines that have the package: a runner that installed it to
#' satisfy Suggests, a laptop where only the package was installed.  Choosing
#' the backend on \code{requireNamespace("cmdstanr")} alone therefore sent
#' every fit into "CmdStan path has not been set yet. See ?set_cmdstan_path",
#' raised from inside the sampler, on such machines.
#'
#' @param has_pkg Whether the \pkg{cmdstanr} package can be loaded.
#' @param path_fn A function returning the CmdStan path, which errors when
#'   none is set.  Both arguments exist so the decision can be tested without
#'   a CmdStan installation.
#' @return \code{TRUE} when a CmdStan directory is set and exists.
#' @keywords internal
#' @noRd
.cmdstan_usable <- function(has_pkg = requireNamespace("cmdstanr", quietly = TRUE),
                            path_fn = function() cmdstanr::cmdstan_path()) {
  if (!isTRUE(has_pkg)) return(FALSE)
  p <- tryCatch(path_fn(), error = function(e) NULL)
  is.character(p) && length(p) == 1L && !is.na(p) && nzchar(p) && dir.exists(p)
}


#' Resolve the Stan backend for brms
#'
#' @param backend One of "auto", "cmdstanr", "rstan" (already matched).
#' @param cmdstan_ok Whether a CmdStan build is usable; see
#'   \code{.cmdstan_usable()}.
#' @return "cmdstanr" or "rstan".
#' @keywords internal
#' @noRd
.resolve_stan_backend <- function(backend, cmdstan_ok = .cmdstan_usable()) {
  if (identical(backend, "auto")) {
    chosen <- if (isTRUE(cmdstan_ok)) "cmdstanr" else "rstan"
    .log_info("fit_bayesian_spatial_model(): backend = \"auto\" resolved to %s%s.",
              chosen,
              if (identical(chosen, "rstan") &&
                  requireNamespace("cmdstanr", quietly = TRUE))
                " (the cmdstanr package is installed but no CmdStan build is; see cmdstanr::install_cmdstan())"
              else "")
    return(chosen)
  }
  if (identical(backend, "cmdstanr") && !isTRUE(cmdstan_ok))
    stop(paste0("fit_bayesian_spatial_model(): backend = \"cmdstanr\" needs the ",
                "cmdstanr package AND a CmdStan build, and ",
                if (requireNamespace("cmdstanr", quietly = TRUE))
                  "no CmdStan path is set (run cmdstanr::install_cmdstan(), or cmdstanr::set_cmdstan_path() if it is already built)"
                else "the cmdstanr package is not installed",
                ". Use backend = \"rstan\" for the backend brms always provides."),
         call. = FALSE)
  backend
}
