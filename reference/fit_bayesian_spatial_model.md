# Fit a Bayesian Gaussian-process regression with brms

Fits a regression whose residual spatial structure is modelled as a
Gaussian process over the coordinates, and so kept out of the errors.
Every quantity has a posterior, so predictions come with intervals;
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
reports their held-out coverage and CRPS. The spatial term is a smooth
surface at a stated resolution. The Gaussian process is fitted through a
reduced-rank basis, and that basis has a floor: an effective range below
which it cannot represent structure, reported in CRS units as
`$info$gp_resolution` and by
[`print()`](https://rdrr.io/r/base/print.html). With the defaults the
floor is about 32 km along each axis on points spread evenly over a 100
km square. Structure shorter than the floor is counted as noise: the
fitted range comes out too long, the residual standard deviation too
large, and the intervals too narrow where the structure is. See "What
the basis resolves" below.

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
  gp_resolution = NULL,
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
  Gaussian process alone, the baseline for asking how much of the
  surface is spatial structure and how much covariate effect, and the
  null model to compare a covariate model with in
  [`compare_models()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models.md).

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
  unchanged with the spatial GP term still in the formula. A factor
  response is accepted only under
  [`brms::categorical()`](https://paulbuerkner.com/brms/reference/brmsfamily.html)
  or an ordinal family, and those fits have no single expected value per
  row, so several methods cannot use them; see "Non-Gaussian responses"
  for what each family supports, and the count example below.

- gp_k:

  Positive integer giving the number of GP basis functions *per
  dimension*, or NULL (default) to derive it from the
  length-scale/domain ratio. The fitted model carries `gp_k^2` basis
  functions, not `gp_k` (see "GP basis count and boundary factor"). A
  value above the derived one makes the basis finer, and the
  length-scale prior's lower tail moves down to the length-scale that
  basis resolves, because a finer basis under the default prior gains
  almost nothing. A value below the derived one leaves the prior where
  it was, mostly below what the basis resolves; that share is stored as
  `$info$gp_adequacy$prior_mass_below`, logged when it passes 10 percent
  and flagged by [`print()`](https://rdrr.io/r/base/print.html). To ask
  for a resolution in CRS units, use `gp_resolution` and leave this
  `NULL`.

- gp_c:

  Positive numeric boundary factor for the approximate GP, or NULL
  (default) to derive it alongside `gp_k`. The boundary must be wide
  enough to contain the longest plausible effective range; a value that
  is too small truncates the domain and degrades the approximation for
  smooth, long-range surfaces. When you set `gp_c` and leave
  `gp_k = NULL`, `gp_k` is derived for *your* boundary: a wider boundary
  needs more basis functions to resolve the same length-scale, so
  raising `gp_c` raises the derived `gp_k` with it, up to the cap of 50
  per dimension (a capped value is logged).

- gp_iso:

  Logical; passed to `brms::gp(iso = )`. `FALSE` (the default) fits a
  separate length-scale per coordinate axis, so the model learns any
  directional structure from the data. `TRUE` fits a single length-scale
  shared by both standardised axes, which makes the kernel anisotropic
  in the original CRS by the ratio `sd(X)/sd(Y)`, because the
  coordinates are standardised per axis beforehand. See "Coordinate
  scaling and anisotropy".

- prior:

  Optional brms prior specification. When NULL and
  `standardize_predictors = TRUE`, a weakly informative prior is set on
  the regression coefficients, on the scale of the linear predictor:
  `normal(0, 2.5 * sd(y))` under an identity link (gaussian, student,
  skew-normal and mixtures of them), `normal(0, 2.5 * sd(log(y)))` for
  the lognormal families, whose location is on the log scale, and the
  fixed `normal(0, 5)` under any other link (log, logit and so on),
  where the linear predictor has a scale of its own.
  [`sd()`](https://rdrr.io/r/stats/sd.html) is taken over the training
  rows. The result therefore does not depend on the units of the
  response. The prior brms used is stored as `$info$b_prior` and shown
  by [`print()`](https://rdrr.io/r/base/print.html). When `prior` is
  given, no slope prior is added. A prior given here on a class that
  lives on the response's scale (`b`, `Intercept`, `sigma`, `sdgp`, and
  `beta` of `exgaussian()`) keeps the response as it is; see "The
  response's centre and scale". A GP length-scale prior set from the
  data is added unless `prior` already contains an entry with
  `class = "lscale"`.

- chains:

  Number of MCMC chains. Default 4.

- iter:

  Total iterations per chain. Default 2000.

- warmup:

  Warmup iterations. Default `floor(iter / 2)`.

- cores:

  Number of cores for the sampler, one chain per core. Default
  `getOption("mc.cores", 1L)`, the convention brms uses itself, so
  `options(mc.cores = 4)` once per session runs the four default chains
  in parallel everywhere.

- seed:

  Integer seed. Default 123.

- backend:

  `"auto"` (default), `"cmdstanr"` or `"rstan"`. `"auto"` uses cmdstanr
  only when a CmdStan build is available, and rstan otherwise, which
  brms always brings. The cmdstanr package is a thin interface and can
  be installed without a CmdStan build
  ([`cmdstanr::install_cmdstan()`](https://mc-stan.org/cmdstanr/reference/install_cmdstan.html)
  builds one). An explicit `"cmdstanr"` with no usable CmdStan raises an
  error that says how to install it, before the sampler starts.

- control:

  Named list of sampler controls, merged over the package defaults
  `list(adapt_delta = 0.99, max_treedepth = 12)` rather than replacing
  them, so `list(max_treedepth = 15)` keeps `adapt_delta = 0.99`. See
  "Sampler settings".

- compute_loo:

  Logical; compute PSIS-LOO. Default TRUE.

- standardize_predictors:

  Logical; centre and scale numeric predictors before fitting. Default
  FALSE. When TRUE, the centre and scale of each predictor are stored in
  `$info$predictor_scaling` and applied again by
  [`predict()`](https://rdrr.io/r/stats/predict.html). The model is then
  fitted on the standardised predictors, so
  [`coef()`](https://rdrr.io/r/stats/coef.html) reports a slope per
  standard deviation of each predictor and an intercept at the predictor
  means, not the raw-unit values
  [`stats::lm()`](https://rdrr.io/r/stats/lm.html) would give (see
  [spatial_fit-coef](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-coef.md)).

- check_convergence:

  Logical, default TRUE. After fitting, read the sampler's diagnostics
  and the GP basis against the posterior and the residuals, and set
  `$info$convergence_ok` and `$info$gp_adequacy`. See "Convergence
  checks" and "What the basis resolves". `FALSE` skips both, leaves
  `convergence_ok` `NA` (not checked) and sets the adequacy status to
  `"not_checked"`.

- pointize:

  How a feature that is not a point is reduced to one: one of `"auto"`
  (default), `"centroid"`, `"point_on_surface"`, `"surface"`,
  `"line_midpoint"` or `"bbox_center"`; see
  [`coerce_to_points()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md).

- boundary:

  Optional polygonal sf/sfc whose CRS the data are projected to. It is
  not a clip: rows outside it are kept.

- gp_resolution:

  The shortest effective range the basis should resolve, as a single
  positive number in the units of the data's CRS (metres for a UTM
  layer), or NULL (default) for the rule described under "GP basis count
  and boundary factor". `gp_k` is derived from it, and the length-scale
  prior's lower tail moves with the basis. It can only refine: a value
  the default basis already resolves leaves the default basis in place.
  It is an error when the basis it needs exceeds the cap of 50 functions
  per dimension (2,500 in all); the message gives the basis size needed
  and the shortest range the cap reaches on this layout. Cost grows with
  `gp_k^2`. Give `gp_resolution` or `gp_k`, not both.

- .already_prepped:

  Internal; leave at the default.

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
gp_ell_min, gp_resolution: the smallest effective range the basis
resolves, as `c(x, y)` in CRS units; gp_resolution_requested: the
`gp_resolution` argument, or `NULL`; gp_prefit: the check made before
sampling, a list with `status` (`"ok"`; `"below_floor"`, the one that
warns; `"below_floor_no_autocorrelation"` and
`"below_floor_long_range"`, a range below the floor that the short-lag
or the floor-distance test did not bear out; `"not_estimated"` or
`"not_run"`), `range` (the estimated range in CRS units), `floor`,
`k_needed`, `moran` and `moran_floor` (the two Moran's I tests, each
with `I`, `z` and a one-sided `p`; `NULL` where the test is not run, or
where the residuals are the rounding of an exact fit), `reason` and
`unit`; gp_adequacy: a list with `status` (see "What the basis
resolves"), `frac_below` (the share of length-scale draws below
`gp_ell_min`), `identified` and `identified_by` (`"spatial_share"` or
`"length_scale"`, the rule that decided it; see "What the basis
resolves"), `contraction` (one minus the ratio of posterior to prior
variance, per length-scale, negative when the posterior is the wider of
the two) and `prior_gap` (the largest gap between the posterior's and
the prior's distribution functions, per length-scale; both `NA` for a
length-scale prior this function did not set), `spatial_share` and
`sdgp` (2.5, 50 and 97.5 percent quantiles; the share only for a
gaussian or lognormal family), `effective_range` (the posterior
effective range per axis in CRS units, same quantiles), `residual_moran`
(the test on the residuals the status rests on: `I`, `z`, a one-sided
`p` and `k`, the number of neighbours; absent when it could not be run,
and `residual_untested` then says why), `resolution`, `unit` and
`prior_mass_below`; b_prior: the slope prior brms used (`"(flat)"` for
its default, `NA` for a model with no slopes), and b_prior_source: where
it came from; predictor_levels: the levels each factor or character
predictor has in the training rows, which
[`predict()`](https://rdrr.io/r/stats/predict.html) holds new data to;
gp_S: the pooled centred range `brms::gp(c = )` multiplies; gp_cmeans:
the column means brms centred the scaled coordinates on; gp_xy_range:
the training extrema of the scaled coordinates, which
[`predict()`](https://rdrr.io/r/stats/predict.html) uses, with gp_S and
gp_cmeans, to hold the GP boundary at its fitted value;
gp_lengthscale_bounds: the `c(lower, upper)` the length-scale prior was
calibrated over; gp_lscale_prior: the length-scale prior
[`brms::validate_prior()`](https://paulbuerkner.com/brms/reference/validate_prior.html)
reports the model will *actually* use, which is not necessarily the one
this function requested (several entries, semicolon-separated, if brms
resolved the axes differently); loo, looic, convergence_ok (`TRUE` or
`FALSE`, and `NA` when nothing was checked, as under
`check_convergence = FALSE`), convergence_diagnostics: `n_divergent`,
`n_max_treedepth`, `ebfmi` (per chain), `max_rhat`, `min_neff_ratio`,
`min_ess` (the smallest bulk or tail effective sample size) and
`ess_required` (100 per chain), and `rhat_failed` / `neff_failed` /
`ess_failed`, the parameters that failed each check by name with their
values (R-hat, ESS ratio and ESS; empty when none failed), which is what
makes a failed check actionable, `rhat_caution`, the parameters whose
R-hat is between 1.01 and 1.05 with enough effective draws, which do not
fail the verdict, plus `gp_lscale_below_resolution` and `gp_adequacy`,
the share and status described under "What the basis resolves"; and
n_dropped: the rows
[`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)
removed for missing or non-finite values or a bad geometry, so `$n` can
be read against `nrow(data_sf)`); response_scaling: the `centre` and
`scale` the response was handed to brms with, absent when it was handed
over as it is (see "The response's centre and scale"). The raw brmsfit
is in `$engine`.

## Details

The fitted length-scale is an estimate of how far the spatial dependence
reaches only when its posterior sits clear of both that floor and its
prior. `$info$gp_adequacy` says whether it does, and when it does
`$info$gp_adequacy$effective_range` holds the posterior effective range
per axis in CRS units (the raw `lscale` draws are in per-axis
standardised units and are not a distance).

Choose it over
[`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
when you want one global relationship plus a spatial random field, with
posterior intervals; choose GWR when the question is how a coefficient
varies across the map. Choose
[`fit_rf_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md)
when predictive accuracy matters more than an interpretable model and
the response is non-linear in the predictors. The cost is time: this is
full MCMC through brms and Stan, so a fit takes minutes rather than
seconds, and the basis size (`gp_k`) trades fidelity against run time.

## Sampler settings

`adapt_delta` defaults to 0.99, not Stan's 0.8, because one divergent
transition sets `convergence_ok` to `FALSE` (see "Convergence checks")
and lower values can leave a few on ordinary data. A lower value, such
as `list(adapt_delta = 0.9)`, runs faster, since the smaller step of
0.99 costs sampling time. Some divergent transitions can remain at 0.99
too, for example on a layer whose spatial field is weak against the
noise; then raise it to 0.995 or 0.999, which costs more time again.

## Convergence checks

Under `check_convergence = TRUE` the verdict in `$info$convergence_ok`
follows Stan's own diagnostics, with these thresholds:

- any divergent transition after warmup;

- bulk or tail effective sample size below 100 per chain for any
  parameter (Vehtari et al. 2021; the rule rstan warns on). It is an
  absolute count, not a share of the draws;

- rank-normalised split R-hat, in two tiers. Above 1.05, the value
  rstan's own sampler warning uses, for any parameter. Above 1.01, the
  value Vehtari et al. (2021) recommend, for a parameter whose effective
  sample size is also below the count in the previous item, or cannot be
  read;

- E-BFMI below 0.3 in any chain (Betancourt 2016).

An R-hat between 1.01 and 1.05 on a parameter with enough effective
draws is a caution, not a failure: it is stored by name in
`$info$convergence_diagnostics$rhat_caution`, logged, and counted by
[`print()`](https://rdrr.io/r/base/print.html). The model has one weight
per basis function, about 540 parameters at the basis derived for points
spread evenly over a square, and among that many the largest R-hat sits
near 1.01 in a fit that has converged. Vehtari et al. give 1.01 beside
an effective sample size above 400, so the two are read together here.

Transitions that stop at `max_treedepth` are counted and logged but do
not change the verdict: Stan treats them as a cost, not as invalid
draws. A failed verdict sets `convergence_ok` to `FALSE`, raises one R
warning of class `spatialkit_convergence` naming each failed check with
its count and its first five parameters, and
[`print()`](https://rdrr.io/r/base/print.html) flags it;
`$info$convergence_diagnostics` holds every value. The warning closes
with the remedy for the checks that failed: a higher `adapt_delta` for
divergent transitions, which more iterations do not remove, and more
iterations for R-hat and effective sample size. A short run usually
fails on effective sample size.
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
marks such folds in `fold_metrics$convergence_ok` and
[`compare_models()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models.md)
names such a model.

The same step checks the GP basis against the posterior and the
residuals (see "What the basis resolves") and stores the result as
`$info$gp_adequacy`. The residuals need the fitted values, one pass over
the posterior that
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) caches for
later calls. That check is about the approximation, not the sampler, so
it does not change `convergence_ok`;
[`print()`](https://rdrr.io/r/base/print.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) show its status, and
two of its states raise a warning of class `spatialkit_gp_adequacy`.

## GP basis count and boundary factor

[`brms::gp()`](https://paulbuerkner.com/brms/reference/gp.html) builds a
full tensor grid over its covariates, so a term `gp(..x, ..y, k = gp_k)`
carries `gp_k^2` basis functions: the `gp_k` argument is the count *per
dimension*, not the total rank. Both `gp_k` and `gp_c` are therefore
chosen from the ratio of the estimated length-scale to the domain
extent, following Riutort-Mayol et al. (2023): `gp_c` is set large
enough to contain the upper length-scale bound, and `gp_k` large enough
to resolve the lower one. The derived value depends on the layout, not
on `n`: about 23 for points spread evenly over a square; other layouts
give other values, up to the cap of 50.

The domain extent used is the one `brms::gp(c = )` itself multiplies:
the full pooled range of the column-centred coordinates
(`brms:::choose_L()`, taken over the unique coordinate rows, because
`brms:::.data_gp()` reduces the covariates to unique rows first under
the default `gr = TRUE`, so repeat visits to one location do not widen
the domain), not the per-axis half-range in which Riutort-Mayol et al.
state their inequalities. Both constraints are on the boundary \\L = c
\times S\\, and stating them in brms's units keeps `gp_c`, `gp_k` and
`$info$gp_ell_min` describing the basis brms builds. A `gp_c` derived on
the half-range convention and handed to
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

## What the basis resolves

The smallest length-scale the chosen basis can resolve is
`1.75 * gp_c * S / gp_k`, stored in standardised units as
`$info$gp_ell_min`. As an effective range (the distance at which the
squared-exponential correlation has fallen to 5 percent, 2.45
length-scales) and in CRS units it is `$info$gp_resolution`, one value
per axis because each coordinate is divided by its own standard
deviation. The floor is soft: the basis keeps 99 percent of a process's
variance at the floor, 69 percent at half of it and 41 percent at a
third. The default length-scale prior is what makes it bind, with 1
percent of its mass below the floor and about one part in a hundred
million below half of it.

Three checks hold the basis against the data.

1.  Before anything is compiled, the autocorrelation range of the
    response (of its OLS residuals when there are predictors) is
    estimated with
    [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
    and compared with the coarser axis's floor. That range is a
    variogram estimate and can be wrong in both directions, so a range
    below the floor raises a warning only when two other things agree
    with it. The same variable is autocorrelated at short range (Moran's
    I over 4 nearest neighbours, one-sided, at 5 percent): a variogram
    fitted to a response with no spatial structure can still return a
    short range, though
    [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
    refuses most such layers itself. And it is no longer autocorrelated
    at the floor's distance (Moran's I over the pairs 0.75 to 1.25
    floors apart is not both above 0.05, the correlation that defines an
    effective range, and significant at 5 percent): a variogram can also
    return a short range for a surface that is smooth over long
    distances. The warning has class `spatialkit_gp_resolution`, gives
    both numbers and the basis size that would resolve the range, says
    that it is an estimate which the check after sampling tests against
    the residuals of the fitted model, and says when the basis size is
    past the cap or past the number of rows, where kriging is the better
    tool if the range is right. A short range the floor-distance test
    contradicts is logged and shown by
    [`print()`](https://rdrr.io/r/base/print.html) without a warning: it
    can still mean short-range structure beside a longer one, which the
    basis counts as noise. Nothing before sampling separates the two;
    the residuals after it do (the third check). A range that cannot be
    estimated is logged with the reason and shown by
    [`print()`](https://rdrr.io/r/base/print.html); the result is in
    `$info$gp_prefit`. An exponential field has structure at every
    scale, so the warning on one is not plainly wrong.

2.  Also before sampling, the share of the length-scale prior below the
    floor is computed (`$info$gp_adequacy$prior_mass_below`). A derived
    basis keeps it near 0.01; a small `gp_k` does not, and a share above
    0.10 is logged and flagged by
    [`print()`](https://rdrr.io/r/base/print.html).

3.  After sampling, under `check_convergence = TRUE`, the posterior and
    the residuals are read. The posterior says whether the Gaussian
    process carries a share of the variance and where its length-scale
    sits against the floor. It cannot say what the surface left behind,
    so the residuals are tested for autocorrelation at short range in
    every fit (Moran's I over 4 nearest neighbours, one-sided, at 5
    percent; `$info$gp_adequacy$residual_moran`).
    `$info$gp_adequacy$status` is one of:

    - `"ok"`:

      The Gaussian process carries a share of the variance, its
      length-scale is clear of the floor, and the residuals show no
      autocorrelation at short range. A surface smoother than the prior
      expects (a plane, a range longer than the extent) is `"ok"`: its
      length-scale posterior sits on the long side of the prior and is
      wider than it.

    - `"below_resolution"`:

      The Gaussian process carries a share of the variance and more than
      10 percent of the length-scale draws fall below the floor. The 10
      percent trigger is this package's operationalisation of the check
      Riutort-Mayol et al. recommend, not a figure from the paper.

    - `"unresolved_short_range"`:

      The residuals are autocorrelated at short range: the signature of
      structure finer than the basis, counted as noise. It is reached
      two ways, and `$info$gp_adequacy$identified` says which. The
      Gaussian process carries no variance that can be told from zero
      (`FALSE`): all of the structure is below the floor. Or it carries
      a share and its length-scale is clear of the floor (`TRUE`): it
      fitted a surface, and finer structure sits beside it, as with a
      regional trend plus short-range structure.

    - `"no_signal_above_floor"`:

      The Gaussian process carries no variance that can be told from
      zero, and the residuals show nothing left over. A response with no
      spatial structure gives this, so it is not a failure by itself.

    - `"not_checked"`:

      `check_convergence = FALSE`.

    The second and third raise a warning of class
    `spatialkit_gp_adequacy`; the fourth is logged. When the residuals
    cannot be tested (an ordinal or categorical family has no residual
    per row; a posterior that cannot be drawn), the status rests on the
    posterior alone, a line is logged, the reason is in
    `$info$gp_adequacy$residual_untested`, and
    [`print()`](https://rdrr.io/r/base/print.html) keeps its flag on a
    variogram range that was below the floor.

    On a fit with nothing left over the test is conservative, because
    the surface takes up part of the noise and leaves the residuals
    slightly negatively autocorrelated. It is not sharp just below the
    floor, where the length-scale draws report the structure
    (`"below_resolution"`). An exponential field has structure at every
    scale and can be reported with a range above the floor. The 5
    percent level is the usual one and is not tuned; with many rows a
    small leftover autocorrelation is significant, so read the size of I
    beside its p.

    "No variance that can be told from zero" is decided on the Gaussian
    process's share of the variance, `sdgp^2 / (sdgp^2 + sigma^2)`: its
    2.5 percent quantile is under 0.01. It does not mark off a response
    whose structure is all below the floor: such a process can still
    carry a share, and the residual test reports it. The share exists
    for a gaussian or lognormal family. For any other family the
    length-scales decide instead: the process is not identified when on
    every axis the posterior kept more than three quarters of its
    prior's variance and the largest gap between its distribution
    function and the prior's is under 0.30.
    `$info$gp_adequacy$identified_by` says which rule was used. The
    width of the length-scale posterior does not decide by itself,
    because a posterior that moved to longer length-scales is wider than
    its prior; and the share of draws below the floor cannot detect an
    unidentified process, because a posterior the data did not move puts
    1 percent below the floor, the same as the prior.

A response whose range is below the floor needs a finer basis and a
prior that allows the shorter length-scale, both of which
`gp_resolution` sets, up to the cap. Past the cap this model cannot
resolve the range, and kriging, which has no such floor, is the better
tool; see
[`kriging_adequacy()`](https://elkronos.github.io/gis_modeling_toolkit/reference/kriging_adequacy.md)
for cell means.

## Coordinate scaling and anisotropy

Before fitting the GP, X and Y coordinates are each centred and divided
by their own standard deviation. This is a conditioning step: easting
and northing frequently span very different ranges in a projected CRS,
and handing Stan raw metres samples poorly.

Because the axes are scaled independently, a *single* shared
length-scale in the scaled space corresponds to an anisotropic kernel in
the original CRS, stretched by the ratio `sd(X)/sd(Y)`. That ratio is a
property of how the sampling locations are laid out, not of the process
being modelled.

`gp_iso = FALSE` (the default) therefore fits one length-scale per axis,
letting the model estimate directional structure from the data instead
of inheriting it from the standardisation. Set `gp_iso = TRUE` for one
length-scale shared by both standardised axes.

`gp_iso` does not affect cost:
[`brms::gp()`](https://paulbuerkner.com/brms/reference/gp.html) builds a
tensor grid either way, so the model carries `gp_k^2` basis functions
regardless. The stored `$info$coord_scaling` list records the scaling
strategy, and `$info$gp_iso` records which kernel was used.

## Warnings

Every warning this function raises carries a condition class, so a
caller can catch or muffle exactly one of them by name:

- `spatialkit_rows_dropped`, `spatialkit_constant_predictor`:

  From
  [`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md):
  incomplete rows were dropped, with the counts; a predictor takes one
  value in every row kept.

- `spatialkit_response_scale`:

  The response lies far from zero or its standard deviation far from 1,
  and `prior` sets a prior on a class on the response's scale, so the
  response is handed to brms as it is; see "The response's centre and
  scale".

- `spatialkit_gp_basis_size`:

  Fewer than one row per ten basis functions.

- `spatialkit_gp_resolution`:

  A variogram fitted before compiling puts the response's range below
  what the basis resolves, and the autocorrelation of the response
  agrees. An estimate; the check after sampling tests it against the
  residuals of the fitted model.

- `spatialkit_convergence`:

  The sampler's diagnostics fail.

- `spatialkit_gp_adequacy`:

  The basis is too coarse for the data: the posterior asks for a shorter
  length-scale than it resolves (`"below_resolution"`), or the residuals
  are autocorrelated at short range (`"unresolved_short_range"`).

A response that takes a single value is an error, raised before anything
is compiled, and so is one whose values differ only by rounding (1e15
and the next double). A fit that comes back from brms with no posterior
draws, as when Stan rejects every starting value, is an error that names
the response's size as brms was given it. The convergence warning
carries the checks that failed as `failed` (`"divergent"`, `"ebfmi"`,
`"rhat"`, `"ess"`).
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
fits once per fold; it collects the last four classes from its folds and
raises each once for the run, naming the folds, and
[`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)
does the same for a `fit_fn` that calls this function.

## The response's centre and scale

Stan starts every parameter between -2 and 2, and brms sets its default
priors for the intercept, the residual standard deviation and the
Gaussian process's standard deviation from the response it is given,
with scales of at least 2.5 in its units. A response at 1e8 plus a
spread of 1, or in units of 1e-20, lies far from both, and brms can miss
the data by orders of magnitude. So under the gaussian, student,
skew_normal, exgaussian and asym_laplace families with an identity link,
a response whose centre is more than 1,000 times its spread from zero is
handed to brms less that centre, and one whose standard deviation lies
outside \\2^{-20}\\ to \\2^{20}\\ is divided by the power of two nearest
that standard deviation. Under exgaussian and asym_laplace the rule is
narrower: a centre more than 4 times the spread from zero is taken off,
and a standard deviation outside \\2^{-2}\\ to \\2^{2}\\ is divided out,
because brms's default prior on exgaussian's `beta` and Stan's starting
values do not follow the response. Both steps are exact, and the model
is the same model of the response: the shift moves the location only,
and the scale multiplies the location, `sigma` and the mean `beta` of
exgaussian's exponential part, and leaves the degrees of freedom of
student, the skewness of skew_normal and the quantile of asym_laplace as
they are. Every other response, and every response under another family,
is handed over as it is.

[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`predict()`](https://rdrr.io/r/stats/predict.html) and its draws,
[`residuals()`](https://rdrr.io/r/stats/residuals.html),
[`coef()`](https://rdrr.io/r/stats/coef.html), `$info$loo`,
`$info$looic` and `$info$gp_adequacy$sdgp` are put back in the
response's units, and `$info$b_prior` gives the automatic slope prior in
those units. The same data from another origin or in other units then
give the same fit, up to Monte Carlo error. The centre and scale are in
`$info$response_scaling`. `$engine` is the fit of the response brms was
given, so brms's own functions on it
([`brms::fixef()`](https://rdrr.io/pkg/nlme/man/fixed.effects.html),
[`brms::posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html))
answer on that scale.

A prior given in `prior` on a class on the response's scale (`b`,
`Intercept`, `sigma`, `sdgp`, and `beta` of exgaussian) is read by brms
on the response it is given. With such a prior the response is handed
over as it is, and the warning `spatialkit_response_scale` says so.

## Non-Gaussian responses

Nothing in this function is Gaussian-specific except its default. A
numeric count, binary (0/1 or logical) or bounded response passes
through to brms under the family you name, with the GP term in the
formula and the length-scale prior attached to each distributional
parameter's GP.

The response check is family-aware. A factor or character response is
accepted only under
[`brms::categorical()`](https://paulbuerkner.com/brms/reference/brmsfamily.html)
or an ordinal family (`cumulative`, `sratio`, `cratio`, `acat`), and
refused under every other family before anything is compiled; under
gaussian a logical response is refused too. The case to watch is a
two-level factor under
[`brms::bernoulli()`](https://paulbuerkner.com/brms/reference/brmsfamily.html):
brms would fit it, but
[`residuals()`](https://rdrr.io/r/stats/residuals.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md)
and
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
could not score a factor, so convert it to 0/1 first.

An ordinal or categorical fit has a probability per response category,
not one expected value per row.
[`predict()`](https://rdrr.io/r/stats/predict.html) with its default
`type = "epred"`,
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html),
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md)
therefore stop with a message saying so, and
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
refuses the family before fitting anything.
`predict(type = "predict", draws = TRUE)` returns the posterior
predicted categories, as category indices, for new rows as well as the
training ones (the share of draws in each category estimates its
probability), and `brms::posterior_epred(fit$engine)` the probabilities
for the training rows.

For the numeric families two things follow. First, the metrics that come
back from
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

The response check reads the family's name through brms's own accessor.
If brms cannot name the family object, the check is skipped, and a wrong
response type surfaces as a Stan error with no message from this
function.

## Spatial confounding

A fixed-effect coefficient estimated alongside a spatial random effect
is a different quantity from the same coefficient in a non-spatial fit.
When a covariate is itself spatially smooth, the GP term absorbs part of
its effect and the two estimates can disagree sharply; when the response
is smoother than the covariate, the estimate can shrink toward zero
regardless of the true effect (Bolin and Wallin 2026). Under a correctly
specified spatial model this is not a bias but a change of estimand
(Zimmerman and Ver Hoef 2022): the spatial coefficient is the effect
*net of* whatever the spatial field can explain, and the non-spatial one
is not. The right one depends on the question, so report both side by
side, unadjusted: fit the same formula with
[`stats::lm()`](https://rdrr.io/r/stats/lm.html) or
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) and compare. Compare
like with like: under `standardize_predictors = TRUE` the coefficients
here are per standard deviation of each predictor, so either fit the
non-spatial model on the same standardised columns or divide these
slopes by `$info$predictor_scaling[[name]]$scale` first (see
[spatial_fit-coef](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-coef.md)).

Published work on remedies is unsettled, and this function takes no
side. Restricted spatial regression (Hughes and Haran 2013) projects the
spatial effect off the covariate space; later work argues it does not
deliver what it promises (Hanks et al. 2015; Khan and Calder 2022).
Within this backend, modelling the covariate and the spatial field as
explicitly correlated (Marques, Kneib and Klein 2022) is expressible in
a brms formula with no new code, and a spectral adjustment (Guan et al.
2023) is the principled route that suits a Hilbert-space basis; this
function does not implement it.

## References

Betancourt, M. (2016). Diagnosing suboptimal cotangent disintegrations
in Hamiltonian Monte Carlo. arXiv:1604.00695.

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and Burkner, P.-C.
(2021). Rank-normalization, folding, and localization: an improved R-hat
for assessing convergence of MCMC (with discussion). *Bayesian Analysis*
**16**, 667–718.
[doi:10.1214/20-BA1221](https://doi.org/10.1214/20-BA1221)

Riutort-Mayol, G., Burkner, P.-C., Andersen, M. R., Solin, A. and
Vehtari, A. (2023). Practical Hilbert space approximate Bayesian
Gaussian processes for probabilistic programming. *Statistics and
Computing* **33**, 17.
[doi:10.1007/s11222-022-10167-2](https://doi.org/10.1007/s11222-022-10167-2)

Bolin, D. and Wallin, J. (2026). Spatial self-confounding:
smoothness-related estimation bias in spatial regression models.
*Biometrika* **113**(1), asaf076.
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

[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
to score the fit on spatial folds;
[`predict()`](https://rdrr.io/r/stats/predict.html)
([spatial_fit-predict](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-predict.md))
and [`coef()`](https://rdrr.io/r/stats/coef.html)
([spatial_fit-coef](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-coef.md));
[`gp_lengthscale_bounds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gp_lengthscale_bounds.md)
for the length-scale range the prior is set over.

Other model fitting:
[`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md),
[`fit_rf_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md),
[`gp_lengthscale_bounds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gp_lengthscale_bounds.md),
[`new_spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/new_spatial_fit.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Not run: fits with Stan, which needs a working C++ toolchain and takes
# minutes of sampling.
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
  # Two short chains keep this to a few minutes; Stan may warn about
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
