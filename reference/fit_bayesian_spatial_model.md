# Fit a Bayesian spatial regression with a 2D Gaussian Process (via brms)

Fits a regression whose residual spatial structure is modelled
explicitly, as a Gaussian process over the coordinates, and so kept out
of the errors. Two things follow, and they are the reasons to reach for
this backend. First, every quantity comes with a posterior, so
predictions carry calibrated intervals instead of point estimates. Score
them with
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md),
which reports held-out interval coverage and CRPS. Second, the fitted
length-scale is itself an estimate of how far the spatial dependence
reaches, a number you can read and report.

## Usage

``` r
fit_bayesian_spatial_model(
  data_sf,
  response_var,
  predictor_vars,
  family = NULL,
  gp_k = NULL,
  gp_c = NULL,
  gp_iso = FALSE,
  prior = NULL,
  chains = 4,
  iter = 2000,
  warmup = floor(iter/2),
  cores = getOption("mc.cores", 1L),
  seed = 123,
  backend = c("auto", "cmdstanr", "rstan"),
  control = list(),
  compute_loo = TRUE,
  standardize_predictors = FALSE,
  check_convergence = TRUE,
  pointize = "auto",
  boundary = NULL,
  .already_prepped = FALSE
)
```

## Arguments

- data_sf:

  An sf object with response, predictors, and geometries.

- response_var:

  Response column name.

- predictor_vars:

  Predictor column names. May be `character(0)` for an intercept-only
  model: the response is then explained by the intercept and the spatial
  Gaussian process alone, which is the right baseline for asking how
  much of the surface is spatial structure and how much covariate effect
  (and the natural null model for comparing against a covariate model
  with
  [`compare_models`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models.md)).

- family:

  A model family accepted by
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html): a
  stats family function such as
  [`poisson()`](https://rdrr.io/r/stats/family.html), or a brms family
  object such as
  [`brms::zero_inflated_poisson()`](https://paulbuerkner.com/brms/reference/brmsfamily.html),
  [`brms::negbinomial()`](https://paulbuerkner.com/brms/reference/brmsfamily.html),
  [`brms::hurdle_poisson()`](https://paulbuerkner.com/brms/reference/brmsfamily.html),
  [`brms::bernoulli()`](https://paulbuerkner.com/brms/reference/brmsfamily.html)
  or
  [`brms::Beta()`](https://paulbuerkner.com/brms/reference/brmsfamily.html).
  Default `NULL`, resolved to
  [`stats::gaussian()`](https://rdrr.io/r/stats/family.html). The family
  reaches
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html)
  unchanged with the spatial GP term still in the formula, so any
  response type brms can fit, this function can fit; see the section on
  non-Gaussian responses and the count example below.

- gp_k:

  Positive integer giving the number of GP basis functions *per
  dimension*, or NULL (default) to derive it from the
  length-scale/domain ratio. The fitted model carries `gp_k^2` basis
  functions, not `gp_k` (see Details).

- gp_c:

  Positive numeric boundary factor for the approximate GP, or NULL
  (default) to derive it alongside `gp_k`. The boundary must be wide
  enough to contain the longest plausible correlation range; a value
  that is too small truncates the domain and degrades the approximation
  for smooth, long-range surfaces.

- gp_iso:

  Logical; passed to `brms::gp(iso = )`. `FALSE` (the default) fits a
  separate length-scale per coordinate axis, letting the model learn any
  directional structure from the data. `TRUE` fits a single shared
  length-scale, which makes the kernel anisotropic in the original CRS
  by whatever ratio `sd(X)/sd(Y)` happens to take, because the
  coordinates are standardised per axis beforehand. See Details.

- prior:

  Optional brms prior specification. When NULL and
  `standardize_predictors = TRUE`, weakly informative `normal(0, 5)`
  priors are set on regression coefficients. A data-informed GP
  length-scale prior is always appended automatically unless `prior`
  already contains an entry with `class = "lscale"`.

- chains:

  Number of MCMC chains. Default 4.

- iter:

  Total iterations per chain. Default 2000.

- warmup:

  Warmup iterations. Default floor(iter/2).

- cores:

  Number of cores for the sampler, one chain per core. Default
  `getOption("mc.cores", 1L)`, the same convention brms uses itself, so
  `options(mc.cores = 4)` once per session runs the four default chains
  in parallel everywhere. The previous default of `detectCores() - 1`
  took every core but one on any machine, which is not what a shared
  server or a check farm wants.

- seed:

  Integer seed. Default 123.

- backend:

  "auto" (default), "cmdstanr", or "rstan". "auto" uses cmdstanr only
  when a CmdStan build is actually available, and rstan otherwise, which
  brms always brings. The cmdstanr package is a thin interface and can
  be installed without a CmdStan build
  ([`cmdstanr::install_cmdstan()`](https://mc-stan.org/cmdstanr/reference/install_cmdstan.html)
  builds one). An explicit "cmdstanr" with no usable CmdStan raises an
  error that says how to install it, so the failure does not come from
  inside the sampler.

- control:

  Named list of sampler controls, *merged* over the package defaults
  `list(adapt_delta = 0.9, max_treedepth = 12)` instead of replacing
  them. Passing `list(max_treedepth = 15)` therefore keeps
  `adapt_delta = 0.9`, which matters, because that is exactly the
  setting the divergence warning tells you to raise.

- compute_loo:

  Logical; compute PSIS-LOO. Default TRUE.

- standardize_predictors:

  Logical; center and scale numeric predictors before fitting. Default
  FALSE. When TRUE, the scaling parameters are stored in the return
  value so predictions can be computed correctly.

- check_convergence:

  Logical; after fitting, check for divergences, low ESS, and high R-hat
  and issue warnings. Default TRUE.

- pointize:

  Strategy for non-point geometry coercion.

- boundary:

  Optional polygonal sf/sfc for CRS harmonization.

- .already_prepped:

  Logical (internal). If `TRUE`, skip the
  [`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)
  call because the caller has already projected, coerced, and filtered
  the data. The data must then have plain POINT geometry (an error is
  raised otherwise). Used by the CV internals to avoid a redundant
  second pass on every fold. End users should leave this at the default
  `FALSE`.

## Value

A `bayesian_fit` object (inherits from `spatial_fit`). Supports
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html),
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`summary()`](https://rdrr.io/r/base/summary.html), and
[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md).
Model-specific metadata lives in `$info` (coords: the names of the
scaled coordinate columns handed to
[`brms::gp()`](https://paulbuerkner.com/brms/reference/gp.html);
coord_scaling, predictor_scaling, gp_k, gp_c, gp_iso, gp_n_basis,
gp_ell_min, gp_S: the pooled centred range `brms::gp(c = )` multiplies;
gp_xy_range: the training extrema of the scaled coordinates, which
[`predict()`](https://rdrr.io/r/stats/predict.html) uses to pin the GP
boundary; gp_lengthscale_bounds: the `c(lower, upper)` the length-scale
prior was calibrated over; gp_lscale_prior: the length-scale prior
[`brms::validate_prior()`](https://paulbuerkner.com/brms/reference/validate_prior.html)
reports the model will *actually* use, which is not necessarily the one
this function requested (several entries, semicolon-separated, if brms
resolved the axes differently); loo, looic, convergence_ok,
convergence_diagnostics: `n_divergent`, `max_rhat`, `min_neff_ratio`,
and `rhat_failed` / `neff_failed`, the parameters that failed each check
by name with their values (empty when none failed), which is what makes
a failed check actionable; and n_dropped: the rows
[`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)
removed for missing or non-finite values or a bad geometry, so `$n` can
be read against `nrow(data_sf)`). The raw brmsfit is in `$engine`.

## Details

Choose it over
[`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
when you want one global relationship plus an explicit spatial random
field, and uncertainty you can defend; choose GWR instead when the
question is how a coefficient *varies* across the map. Choose
[`fit_rf_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md)
when predictive accuracy matters more than an interpretable model and
the response is non-linear in the predictors. The cost here is time:
this is full MCMC via 'brms' and Stan, so it is minutes rather than
seconds, and the GP is fitted through a reduced-rank basis approximation
whose size (`gp_k`) trades fidelity against runtime.

**GP basis count and boundary factor.**
[`brms::gp()`](https://paulbuerkner.com/brms/reference/gp.html) builds a
full tensor grid over its covariates, so a term `gp(..x, ..y, k = gp_k)`
carries `gp_k^2` basis functions: the `gp_k` argument is the count *per
dimension*, not the total rank. Both `gp_k` and `gp_c` are therefore
chosen from the ratio of the estimated length-scale to the domain
extent, following Riutort-Mayol et al. (2023): `gp_c` is set large
enough to contain the upper length-scale bound, and `gp_k` large enough
to resolve the lower one. The derived value is typically 21-25 per
dimension and is largely independent of `n`.

The domain extent used is the one `brms::gp(c = )` itself multiplies:
the full pooled range of the column-centred coordinates
(`brms:::choose_L()`, taken over the **unique** coordinate rows, because
`brms:::.data_gp()` reduces the covariates to unique rows first under
the default `gr = TRUE`, so repeat visits to one location do not widen
the domain), not the per-axis half-range in which Riutort-Mayol et al.
state their inequalities. Both constraints are really constraints on the
boundary \\L = c \times S\\, so expressing them in brms's units is what
keeps `gp_c`, `gp_k` and `$info$gp_ell_min` describing the basis brms
actually builds. A `gp_c` derived on the half-range convention and
handed to
[`brms::gp()`](https://paulbuerkner.com/brms/reference/gp.html) produces
a boundary twice as wide as intended, against which `gp_k`
under-resolves by a factor of two.

The GP term is built with `scale = FALSE`.
[`brms::gp()`](https://paulbuerkner.com/brms/reference/gp.html)
otherwise rescales its covariates so the maximum Euclidean distance
between two points is 1, and reports `lscale` in that space; since this
function already standardises the coordinates, and the length-scale
prior, `gp_c` and the adequacy check below are all expressed in those
standardised units, a second normalisation would leave every
length-scale quantity in the wrong units.

After fitting, the posterior length-scale is compared against the
smallest scale the chosen basis can resolve (`1.75 * gp_c * S / gp_k`,
stored as `$info$gp_ell_min`); a warning is issued when more than 10% of
the posterior mass falls below it, which is the signal that `gp_k`
should be raised.

**Coordinate scaling and anisotropy.** Before fitting the GP, X and Y
coordinates are each centred and divided by their own standard
deviation. This is a conditioning step: easting and northing frequently
span very different ranges in a projected CRS, and handing Stan raw
metres samples poorly.

Because the axes are scaled independently, a *single* shared
length-scale in the scaled space corresponds to an anisotropic kernel in
the original CRS, stretched by whatever ratio `sd(X)/sd(Y)` happens to
take. That ratio is a property of how the sampling locations are laid
out, not of the process being modelled, so it is not a defensible source
of anisotropy.

`gp_iso = FALSE` (the default) therefore fits one length-scale per axis,
letting the model estimate directional structure from the data instead
of inheriting it from the standardisation. Set `gp_iso = TRUE` to
recover the previous single-length-scale behaviour.

`gp_iso` does not affect cost:
[`brms::gp()`](https://paulbuerkner.com/brms/reference/gp.html) builds a
tensor grid either way, so the model carries `gp_k^2` basis functions
regardless. The stored `$info$coord_scaling` list records the scaling
strategy, and `$info$gp_iso` records which kernel was used.

## Non-Gaussian responses

Nothing in this function is Gaussian-specific except its default. The
response check is family-aware: a non-numeric response is refused only
when the family resolves to gaussian, so a count, binary or bounded
response passes straight through to brms under the family you name.
Zero-inflated and hurdle counts, negative binomial, Bernoulli, beta and
ordinal families have all been verified to reach
[`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html) with
the GP term intact.

Two things follow. First, the metrics that come back from
[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md)
and the `cv_*()` functions are not all meaningful for such a response:
RMSE and MAE are, MAPE, SMAPE and R-squared are Gaussian-shaped, and for
this backend
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)'s
CRPS and interval coverage are the proper scores to read. See
[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md),
section "Which metrics survive a non-Gaussian response". Second,
[`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
fits a variogram to the raw response (or OLS residuals), whose variance
tracks its mean for a count; the range it reports is then less
trustworthy than for a Gaussian response, and its help page says how.

One trap. The response check reads the family's name through `brms`'s
own accessor; a family object it cannot name is treated as "not
gaussian" and the check is skipped entirely, without falling back to the
gaussian rule. A malformed `family` therefore buys less validation, not
more, and a wrong response type will surface as a Stan error, with no
message from this function.

## Spatial confounding

A fixed-effect coefficient estimated alongside a spatial random effect
is a different quantity from the same coefficient in a non-spatial fit.
When a covariate is itself spatially smooth, the GP term absorbs part of
its effect and the two estimates can disagree sharply; when the response
is smoother than the covariate, the estimate can shrink toward zero
regardless of the true effect (Bolin and Wallin 2025). Under a correctly
specified spatial model this is not a bias but a change of estimand
(Zimmerman and Ver Hoef 2022): the spatial coefficient is the effect
*net of* whatever the spatial field can explain, and the non-spatial one
is not. Which of the two a user wants depends on the question, so the
honest diagnostic is to report both side by side and leave them
unadjusted: fit the same formula with
[`stats::lm()`](https://rdrr.io/r/stats/lm.html) or
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) and compare.

The literature on remedies is unsettled and this function takes no side.
Restricted spatial regression (Hughes and Haran 2013) projects the
spatial effect off the covariate space; later work argues it does not
deliver what it promises (Hanks et al. 2015; Khan and Calder 2022).
Within this backend, modelling the covariate and the spatial field as
explicitly correlated (Marques, Kneib and Klein 2022) is expressible in
a brms formula with no new code, and a spectral adjustment (Guan et al.
2023) is the principled route a Hilbert-space basis is well placed to
support later.

## References

Riutort-Mayol, G., Burkner, P.-C., Andersen, M. R., Solin, A. and
Vehtari, A. (2023). Practical Hilbert space approximate Bayesian
Gaussian processes for probabilistic programming. *Statistics and
Computing* **33**, 17.
[doi:10.1007/s11222-022-10167-2](https://doi.org/10.1007/s11222-022-10167-2)

Bolin, D. and Wallin, J. (2025). Spatial self-confounding:
smoothness-related estimation bias in spatial regression models.
*Biometrika* **113**.
[doi:10.1093/biomet/asaf076](https://doi.org/10.1093/biomet/asaf076)

Guan, Y., Page, G. L., Reich, B. J., Ventrucci, M. and Yang, S. (2023).
Spectral adjustment for spatial confounding. *Biometrika* **110**,
699–719.
[doi:10.1093/biomet/asac069](https://doi.org/10.1093/biomet/asac069)

Hanks, E. M., Schliep, E. M., Hooten, M. B. and Hoeting, J. A. (2015).
Restricted spatial regression in practice: geostatistical models,
confounding, and robustness under model misspecification.
*Environmetrics* **26**, 243–254.
[doi:10.1002/env.2331](https://doi.org/10.1002/env.2331)

Hughes, J. and Haran, M. (2013). Dimension reduction and alleviation of
confounding for spatial generalized linear mixed models. *Journal of the
Royal Statistical Society: Series B* **75**, 139–159.
[doi:10.1111/j.1467-9868.2012.01041.x](https://doi.org/10.1111/j.1467-9868.2012.01041.x)

Khan, K. and Calder, C. A. (2022). Restricted spatial regression
methods: implications for inference. *Journal of the American
Statistical Association* **117**, 482–494.
[doi:10.1080/01621459.2020.1788949](https://doi.org/10.1080/01621459.2020.1788949)

Marques, I., Kneib, T. and Klein, N. (2022). Mitigating spatial
confounding by explicitly correlating Gaussian random fields.
*Environmetrics* **33**, e2727.
[doi:10.1002/env.2727](https://doi.org/10.1002/env.2727)

Zimmerman, D. L. and Ver Hoef, J. M. (2022). On deconfounding spatial
confounding in linear models. *The American Statistician* **76**,
159–167.
[doi:10.1080/00031305.2021.1946149](https://doi.org/10.1080/00031305.2021.1946149)

## See also

Other model fitting:
[`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md),
[`fit_rf_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md),
[`gp_lengthscale_bounds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gp_lengthscale_bounds.md),
[`new_spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/new_spatial_fit.md),
[`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Not run: fits with Stan, which needs a working C++ toolchain and takes
# minutes of MCMC -- both outside what an example may assume.
if (requireNamespace("brms", quietly = TRUE)) {
  library(sf)
  set.seed(1)
  n <- 60
  dat <- st_as_sf(
    data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
               elev = rnorm(n)),
    coords = c("x", "y"), crs = 32632
  )
  dat$price <- 10 + 0.01 * (st_coordinates(dat)[, 1] - 5e5) +
    2 * dat$elev + rnorm(n)
  # Two short chains keep this to a few minutes; expect Stan to warn about
  # effective sample size, which is the cost of that.  Use its defaults
  # (chains = 4, iter = 2000) for a fit to report.
  fit <- fit_bayesian_spatial_model(dat, "price", "elev",
                                    chains = 2, iter = 1000,
                                    compute_loo = FALSE)
  print(summary(fit))
  head(predict(fit, newdata = dat))

  # A zero-inflated count.  The family is the only thing that changes.
  # summary()'s R2 and MAPE assume a Gaussian response, so for a count
  # read the coefficients here and score the model with cv_bayes(), whose
  # CRPS and interval coverage are defined for any response.
  dat$count <- rpois(n, exp(0.5 + 0.8 * dat$elev)) * rbinom(n, 1, 0.7)
  fit_zip <- fit_bayesian_spatial_model(dat, "count", "elev",
                                        family = brms::zero_inflated_poisson(),
                                        chains = 2, iter = 1000,
                                        compute_loo = FALSE)
  coef(fit_zip)
}
} # }
```
