# Additions to the annotated bibliography for `spatialkit`

Thirty peer-reviewed works **not among the 418 entries** of
`spatialkit_annotated_bibliography.md`, found by searching each of that
document's twelve sections for the top journals in the sub-field and for
2023–2026 work published after the review's own sources. Every candidate was
checked against the 418 by author, year and title before being kept.

**Verification.** Twenty-six entries were resolved individually through the
OpenAlex works API against their DOI, and title, journal, volume, year and
authors were compared to what is written here. Four have no DOI or resisted the
API and were confirmed against the publisher's or repository's page instead;
those are marked *(publisher page)*. Nothing here is cited from memory alone.
Six further works surfaced in the search are preprints and are listed at the
end under "Watch list" rather than in the numbered set.

**Reading the entries.** Numbering continues the bibliography's section
scheme (`§n.m+`), so an entry can be merged into the existing document under
the same section. Each carries the same relation tags the bibliography uses —
**supports / corrects / extends / informs** — and names the package function
or `dev/BACKLOG.md` item it bears on. The paragraph says *what changes* if the
work is taken seriously, not merely what it is about.

**Where the additions land, in one line each.**

* The `C_p` criterion in stage 5 now has its source paper; the k-means restart
  budget in 5.1 has a citable figure; and there is a *stability-based* number-
  of-clusters criterion that fits the "return a profile" design.
* Stage 6 gains three alternatives to sample splitting that do not throw away
  half the data, and an exact selective-inference method for k-means.
* Stage 4 and 10.3 gain formal isotropy tests to replace the ratio-of-ranges
  heuristic.
* Item 5.6 and REG-K8 gain their change-of-support citations.
* The declined item REG-K4 gains a paper that sharpens exactly when the
  `deff >= 1` clamp is right and when it is not.
* The Bayesian backend gains the spatial-confounding literature the review
  named but the bibliography did not carry.
* Two reviews and a Nature Communications perspective document the state of
  spatial-ML validation practice against which the package's evidence layer is
  positioned.

---

## §1+ · Change of support

**1.36+ Cressie (1996).** Change of support and the modifiable areal unit
problem. *Geographical Systems*, 3(2–3), 159–180. *(publisher page:
ro.uow.edu.au/infopapers/2392)*
  The paper that puts the MAUP on geostatistical footing: aggregating a point
  process to blocks is a change of support, and the variance of a block
  average is a computable function of the point variogram and the block
  geometry (Krige's additivity relation, the dispersion variance). This is the
  citation **5.6** currently lacks for its central claim that the between- and
  within-block variance terms are computable *before* anything is aggregated.
  It also reframes `summarize_by_cell()`'s `deff` correction as one term of a
  decomposition the package could report in full.
  **→ supports/extends** · `summarize_by_cell()`, BACKLOG 5.6, REG-K8

**1.37+ Bradley, Wikle & Holan (2016).** Bayesian spatial change of support
for count-valued survey data with application to the American Community
Survey. *Journal of the American Statistical Association*, 111(514), 472–487.
DOI: 10.1080/01621459.2015.1117471
  The modern model-based treatment of moving count data between incompatible
  supports, with the uncertainty carried through rather than dropped. Two
  things for the package. First, it is the principled answer to what a
  block-kriging aggregator (REG-K8) is for a *count* response — the case
  1.3 already flags as breaking the variogram's stationarity assumption.
  Second, the bibliography carries the same authors' 2017 JRSS-B regionalization
  paper (CAGE) but not this one, so the change-of-support half of that line of
  work was missing.
  **→ extends** · `summarize_by_cell()`, `estimate_sac_range()`, BACKLOG 1.3,
  REG-K8

**1.38+ Young & Gotway (2007).** Linking spatial data from different sources:
the effects of change of support. *Stochastic Environmental Research and Risk
Assessment*, 21(5), 589–600. DOI: 10.1007/s00477-007-0136-z
  A worked demonstration that fitting a model to data on one support and
  applying it on another biases both coefficients and their standard errors
  in predictable directions. Directly relevant to `predict_surface()`, which
  fits on points (or cell means) and predicts onto a regular grid: the grid's
  support is neither. The paper supplies the language a help-page caveat needs
  and the direction of the bias.
  **→ corrects/informs** · `predict_surface()`, `summarize_by_cell()`,
  BACKLOG 5.6

## §2+ · Tessellation

**2.34+ Carr, Olsen & White (1992).** Hexagon mosaic maps for display of
univariate and bivariate geographical data. *Cartography and Geographic
Information Systems*, 19(4), 228–236. DOI: 10.1559/152304092783721231
*(publisher page)*
  The origin of hexagonal binning as a *display* device, and an honest one:
  the argument for hexagons here is visual — no preferred axis, uniform
  neighbour distance on the page — not statistical. Worth citing beside the
  contested-advantage entries (Birch 2007; Khalfa & Hardyns 2025) because it
  shows where the hexagon preference actually comes from. Supports the
  backlog's demotion of hexagonal blocks from an item to a note under 7.1.
  **→ informs** · `create_grid_polygons()`, `plot_tessellation_map()`,
  BACKLOG 7.1

## §3+ · Choosing the number of cells

**3.37+ Mallows (1973).** Some comments on *C_p*. *Technometrics*, 15(4),
661–675. DOI: 10.1080/00401706.1973.10489103
  The source of the criterion stage 5.4 is built on, and it was not in the
  bibliography. Two things worth taking from the original rather than from
  secondary accounts. Mallows is explicit that *C_p* is a **tool for
  displaying** the bias–variance trade-off across candidate models, not a
  rule that picks one — which is exactly the "return a profile, not an
  integer" argument of 5.5, made by the criterion's author. And the paper
  discusses the case where the noise variance must itself be estimated, which
  is the role the nugget (4.2) plays in the resolution problem; the criterion's
  sharpness depends on that estimate, and Mallows says so.
  **→ supports** · `determine_optimal_levels()`, BACKLOG 4.2, 5.4, 5.5

**3.38+ Fränti & Sieranoja (2019).** How much can k-means be improved by using
better initialization and repeats? *Pattern Recognition*, 93, 95–112.
DOI: 10.1016/j.patcog.2019.04.014
  A large empirical study of exactly the two levers 5.1 proposes: k-means++
  initialisation and the number of restarts. The findings that matter here:
  repeats help most when clusters overlap, which for a spatial point pattern
  is the normal case rather than the exception; initialisation quality and
  repeat count trade off, so k-means++ reduces how many repeats are needed but
  does not remove the need; and the improvement from repeats saturates, which
  is the empirical basis for choosing a fixed budget rather than "as many as
  possible". Replace the bare Steinley (2003) citation in 5.1 with this as the
  primary source for the number.
  **→ supports** · `determine_optimal_levels()`, `get_voronoi_seeds()`,
  BACKLOG 5.1, 1.5

**3.39+ Tibshirani & Walther (2005).** Cluster validation by prediction
strength. *Journal of Computational and Graphical Statistics*, 14(3),
511–528. DOI: 10.1198/106186005X59243
  A number-of-clusters criterion built on **stability under resampling**: split
  the data, cluster each half, and ask how well one half's clusters predict
  the other's. It answers a different question from the elbow or from *C_p* —
  not "which k best represents the field" but "which k would I get again on a
  fresh sample" — and that is precisely the zonation-instability concern the
  bibliography's Openshaw & Taylor (1979) entry raises about the package.
  Belongs on the 5.5 profile as a column beside the between-restart spread,
  which measures the same thing more cheaply and less rigorously.
  **→ extends** · `determine_optimal_levels()`, BACKLOG 5.5, 6.1

**3.40+ Steinley (2006).** K-means clustering: a half-century synthesis.
*British Journal of Mathematical and Statistical Psychology*, 59(1), 1–34.
DOI: 10.1348/000711005X48266
  The review that the bibliography's Steinley (2003) entry is the empirical
  half of. Useful chiefly for its treatment of local optima and of the
  standardisation question: k-means on coordinates is scale-dependent, and a
  projected CRS makes the two axes commensurable in a way lon/lat does not —
  which is a reason `determine_optimal_levels()` must run in projected units
  that the package enforces but does not currently cite.
  **→ supports** · `determine_optimal_levels()`, `ensure_projected()`,
  BACKLOG 5.1

**3.41+ Feng, Barcelos, Gaboardi, Knaap, Wei, Wolf, Zhao & Rey (2022).**
spopt: a Python package for solving spatial optimization problems in PySAL.
*Journal of Open Source Software*, 7(74), 3330. DOI: 10.21105/joss.03330
  The reference implementation of the regionalization algorithms the
  bibliography discusses (max-p, SKATER, REDCAP, AZP, Ward with contiguity).
  Peer-reviewed software rather than a method paper, but it is what a user who
  wants constrained regionalization will actually reach for, and the package
  should name it as the neighbour it is not competing with — the same division
  of labour 1.4 draws with `blockCV`.
  **→ informs** · `build_tessellation()`, BACKLOG 1.4, 7.1

## §4+ · Moran's I

**4.44+ de Jong, Sprenger & van Veen (1984).** On extreme values of Moran's I
and Geary's c. *Geographical Analysis*, 16(1), 17–24.
DOI: 10.1111/j.1538-4632.1984.tb00797.x *(publisher page)*
  The attainable range of Moran's I is set by the eigenvalues of the weights
  matrix, so it changes whenever the neighbour structure changes. That is the
  sharpest available reason for something the package already does but does
  not explain: `determine_optimal_levels()` ranks on |z|, not |I|, because the
  weights matrix is rebuilt at every k and the scale of I moves with it. Cite
  it in the C2 row of the backlog's "already implemented" table and at 5.4,
  and consider reporting the bounds alongside I in `residual_morans_i()` so a
  value of 0.3 can be read against its own maximum.
  **→ supports/extends** · `residual_morans_i()`,
  `determine_optimal_levels()`, BACKLOG 5.4, REG-C2

## §5+ · Variograms and anisotropy

**5.41+ Guan, Sherman & Calvin (2004).** A nonparametric test for spatial
isotropy using subsampling. *Journal of the American Statistical
Association*, 99(467), 810–821. DOI: 10.1198/016214504000001150
  A formal test, with a p-value, for the thing `estimate_sac_range()` decides
  by a ratio-of-directional-ranges heuristic. The package's current rule — call
  it anisotropic when the range ratio exceeds a threshold and take the maximum
  — has no error rate. This gives one, and does so without assuming a
  parametric variogram, so it can be applied to the empirical variograms 10.3
  proposes to retain. The cost is subsampling, which on the point counts the
  package targets is affordable.
  **→ corrects/extends** · `estimate_sac_range()`, BACKLOG 4.3, 10.3

**5.42+ Weller & Hoeting (2016).** A review of nonparametric hypothesis tests
of isotropy properties in spatial data. *Statistical Science*, 31(3),
305–324. DOI: 10.1214/16-STS547
  The survey that places Guan et al. among its alternatives and, usefully,
  tabulates which tests need what — lattice versus irregular sampling, sample
  sizes at which each has power, and the `spTest` implementation. The point
  for the package is a caution the review states plainly: at a few hundred
  irregular points, most isotropy tests have little power, so a failure to
  reject is not evidence of isotropy. That should temper the message 4.3 and
  10.3 attach to a directional result.
  **→ informs** · `estimate_sac_range()`, BACKLOG 4.3, 10.3

## §6+ · Geographically weighted regression

**6.39+ Yu, Fotheringham, Li, Oshan, Kang & Wolf (2020).** Inference in
multiscale geographically weighted regression. *Geographical Analysis*,
52(1), 87–106. DOI: 10.1111/gean.12189
  Derives the hat matrix, effective number of parameters and a corrected
  critical value for MGWR — the inferential machinery that turns local
  coefficient surfaces from pictures into tests. Two consequences. The
  bandwidth-specific effective parameter count is what a per-location
  condition-number diagnostic (10.2) should be normalised against; and the
  corrected critical value is the multiple-testing adjustment the bibliography's
  da Silva & Fotheringham (2015) entry asks for, now derived for the
  multiscale case. The GWR backend is single-scale today (`gwr.basic`); this is
  the paper to read before wiring `gwr.multiscale`.
  **→ extends** · `fit_gwr_model()`, BACKLOG 9.6, 10.2, REG-F1..F7

**6.40+ Kao & Fotheringham (2026).** Software for MGWR: a comparison of six
algorithms. *Transactions in GIS*, 30(2). DOI: 10.1111/tgis.70211
  A head-to-head of the MGWR implementations, `GWmodel::gwr.multiscale` among
  them, on speed and on whether they reach the same answer. Directly relevant
  because the package has a hard dependency on one of the compared
  implementations and REG-F1..F7 contemplates exposing its multiscale routine:
  before that, know how it compares. Also a reminder that the bibliography's
  MGWR entries predate this comparison.
  **→ informs** · `fit_gwr_model()`, `gwr_model_selection()`, REG-F1..F7

**6.41+ Jiao & Tao (2025).** Geographical Gaussian process regression: a
spatial machine-learning model based on spatial similarity. *Geographical
Analysis*, 57(3). DOI: 10.1111/gean.12423
  A hybrid that puts a Gaussian process inside a geographically weighted
  frame. It sits between two of the package's backends and is worth knowing as
  a neighbour: it makes explicit that "coefficients vary in space" (GWR) and
  "residual structure varies in space" (GP) are different claims that one model
  can carry at once. Not a candidate backend — the package deliberately keeps
  three — but the right citation when a user asks why the two existing ones
  disagree.
  **→ informs** · `fit_gwr_model()`, `fit_bayesian_spatial_model()`,
  `compare_models_cv()`

## §7+ · The Bayesian backend and spatial confounding

The review's Tier 3 item 19 asks for a spatial-confounding caveat on the GP
backend's coefficients and names Paciorek (2010), Hanks et al. (2015), Khan &
Calder (2022) and Bolin & Wallin (2025), all of which the bibliography carries.
The four below are the ones it does not, and together they are the current
state of the argument.

**7.41+ Hughes & Haran (2013).** Dimension reduction and alleviation of
confounding for spatial generalized linear mixed models. *Journal of the
Royal Statistical Society: Series B*, 75(1), 139–159.
DOI: 10.1111/j.1467-9868.2012.01041.x
  The paper that made "restricted spatial regression" mainstream by projecting
  the spatial random effect off the covariate space, and a fast basis for it.
  The bibliography argues (via Hanks 2015 and Khan & Calder) that RSR should
  *not* be added; this is the paper being argued against, and a caveat that
  names only one side of a live dispute is weaker than one that names both.
  **→ informs** · `fit_bayesian_spatial_model()`, BACKLOG 1.1

**7.42+ Zimmerman & Ver Hoef (2022).** On deconfounding spatial confounding in
linear models. *The American Statistician*, 76(2), 159–167.
DOI: 10.1080/00031305.2021.1946149
  A short, clarifying paper: under the standard linear spatial model,
  "spatial confounding" is not a bias in the fixed-effect estimate at all when
  the model is correctly specified — the apparent conflict is between two
  different estimands. This is the cleanest statement of *why* the package
  should report spatial and non-spatial coefficients side by side (the
  review's recommendation) rather than adjust one toward the other: they
  estimate different things, and the user has to choose which they mean.
  **→ corrects/informs** · `fit_bayesian_spatial_model()`, `summary()` on a
  `bayesian_fit`, BACKLOG 1.1

**7.43+ Guan, Page, Reich, Ventrucci & Yang (2023).** Spectral adjustment for
spatial confounding. *Biometrika*, 110(3), 699–719.
DOI: 10.1093/biomet/asac069
  The current constructive proposal: treat confounding as a scale-dependent
  phenomenon and adjust in the spectral domain, so that fine-scale variation
  in the covariate identifies the effect while coarse-scale variation is
  discounted. Relevant to an HSGP backend in particular, because the
  Hilbert-space basis *is* a spectral representation — the machinery to do
  this is closer to hand in `fit_bayesian_spatial_model()` than in most
  implementations. Not a near-term item; the right citation for a "what would
  a principled fix look like" sentence in the caveat.
  **→ extends** · `fit_bayesian_spatial_model()`, `gp_lengthscale_bounds()`

**7.44+ Marques, Kneib & Klein (2022).** Mitigating spatial confounding by
explicitly correlating Gaussian random fields. *Environmetrics*, 33(5),
e2727. DOI: 10.1002/env.2727
  Models the covariate and the spatial effect as correlated random fields
  rather than orthogonalising them. Cited here for one reason: it is the
  approach that maps most naturally onto a `brms` formula, so if a user asks
  what they can *do* about confounding within the package's existing backend,
  this is the honest answer, and it needs no new code.
  **→ informs** · `fit_bayesian_spatial_model()`, BACKLOG 1.1

## §8+ · Random forests and spatial machine learning

**8.46+ Talebi, Peeters, Otto & Tolosana-Delgado (2022).** A truly spatial
random forests algorithm for geoscience data analysis and modelling.
*Mathematical Geosciences*, 54, 31–60. DOI: 10.1007/s11004-021-09946-w
  Builds spatial context into the forest by feeding each tree local
  neighbourhood statistics rather than raw coordinates. It is the direct
  alternative to the package's `include_coords = FALSE` default: instead of
  refusing location to the forest, give it location *summarised* so it cannot
  memorise. Belongs beside RF-GLS (REG-H4) as the other way to make a forest
  spatial, and it is the one that stays inside `ranger`'s interface.
  **→ extends** · `fit_rf_model()`, REG-H4

**8.47+ Meyer & Pebesma (2022).** Machine learning-based global maps of
ecological variables and the challenge of assessing them. *Nature
Communications*, 13, 2208. DOI: 10.1038/s41467-022-29838-9
  The most-cited statement of the package's evidence-layer thesis: a global
  map from a machine-learning model is only as trustworthy as its validation,
  random cross-validation of clustered samples overstates accuracy, and the
  area of applicability is what stops a map being believed where it should
  not be. The bibliography has the authors' 2021 AOA paper but not this, and
  this is the one to cite in `DESCRIPTION`-level positioning. Also the clearest
  motivation for making `check-brms` and blocked folds the default evidence.
  **→ supports** · `area_of_applicability()`, `make_folds()`, BACKLOG 1.4,
  2.6, 3.5

**8.48+ Koldasbayeva, Tregubova, Gasanov, Zaytsev, Petrovskaia & Burnaev
(2024).** Challenges in data-driven geospatial modeling for environmental
research and practice. *Nature Communications*, 15.
DOI: 10.1038/s41467-024-55240-8
  A perspective from outside the geography-and-ecology group that produced most
  of the bibliography's validation entries, reaching the same conclusions —
  spatial dependence breaks random splits, uncertainty must be reported, and
  generalisation to unsampled regions is the open problem. Independent
  convergence is worth more than another paper from the same school.
  **→ supports** · BACKLOG 1.4

**8.49+ Kopczewska (2022).** Spatial machine learning: new opportunities for
regional science. *The Annals of Regional Science*, 68(3), 713–755.
DOI: 10.1007/s00168-021-01101-x
  A survey of how spatial structure enters machine-learning models — as
  features, as sample weights, as validation design, as post-hoc diagnostics —
  with a taxonomy that the package's three backends map onto cleanly. Useful
  for 1.4 because it positions "cut the data into regions first, then model"
  as one recognised strategy among several, with named alternatives.
  **→ informs** · BACKLOG 1.4

**8.50+ Nikparvar & Thill (2021).** Machine learning of spatial data. *ISPRS
International Journal of Geo-Information*, 10(9), 600.
DOI: 10.3390/ijgi10090600
  The other general review, organised by data type rather than by method.
  Cited for completeness beside Kopczewska; between them they cover the
  literature a reader of the package's README would expect it to have engaged
  with.
  **→ informs** · BACKLOG 1.4

## §10+ · Area of applicability and uncertainty

**10.18+ Lou, Luo & Meng (2025).** GeoConformal prediction: a model-agnostic
framework for measuring the uncertainty of spatial prediction. *Annals of the
American Association of Geographers*, 115(8), 1971–.
DOI: 10.1080/24694452.2025.2516091
  Conformal prediction made spatial: prediction intervals with a finite-sample
  coverage guarantee, computed from the residuals of held-out data weighted by
  spatial proximity, for any backend. This is the model-agnostic complement to
  the AOA: the AOA says *where* a prediction should not be trusted, this says
  *how much* it should be, and both come from the same held-out folds. The
  bibliography has Mao, Martin & Reich (2023) for the statistical foundation;
  this is the applied, geographer-facing version, and the one a `spatialkit`
  user is likelier to have read.
  **→ extends** · `area_of_applicability()`, `predict_surface(se = )`,
  `cv_spatial()`, REG-J3/J4

## §11+ · Effective sample size and post-selection inference

**11.32+ Ferrer & Vallejos (2025).** Is the effective sample size always less
than *n*? A spatial regression approach. *Statistics & Probability Letters*,
218, 110309. DOI: 10.1016/j.spl.2024.110309
  Establishes the conditions under which the spatial effective sample size is
  bounded by *n* — and, by stating them, shows the bound is not universal:
  under negative spatial dependence the effective sample size can exceed the
  observed one. This bears directly on the declined item REG-K4, which keeps
  the `deff >= 1` clamp. The resolution is sharper than "keep" or "drop". On
  the `"variogram"` path the correlation function is non-negative by
  construction (exponential, spherical, Gaussian), so `deff >= 1` is a theorem
  there and the clamp never binds. On the `"kish"` path a negative ICC is a
  legitimate estimate under repulsive sampling, and clamping it discards real
  information. Keep the clamp on one path, document why on the other.
  **→ corrects** · `summarize_by_cell()`, BACKLOG 2.2, REG-K4

**11.33+ Chen & Witten (2023).** Selective inference for k-means clustering.
*Journal of Machine Learning Research*, 24(152), 1–41.
jmlr.org/papers/v24/22-0371.html *(publisher page; JMLR issues no DOIs)*
  An exact test for a difference in means between two clusters found by
  k-means on the same data, conditioning on the clustering event. This is the
  post-selection problem of stage 6 stated for exactly the algorithm
  `determine_optimal_levels()` and `get_voronoi_seeds()` use. It does not
  cover the downstream regression, but it settles the first link in the
  chain: whether two adjacent cells genuinely differ, after having chosen the
  cells by looking. The bibliography's Gao, Bien & Witten (2022) is the
  hierarchical-clustering predecessor.
  **→ extends** · `determine_optimal_levels()`, `summarize_by_cell()`,
  BACKLOG 6.1

**11.34+ Leiner, Duan, Wasserman & Ramdas (2023).** Data fission: splitting a
single data point. *Journal of the American Statistical Association*,
120(549), 135–146. DOI: 10.1080/01621459.2023.2270748
  Splits each observation into two dependent pieces by adding and subtracting
  noise, so that selection can be done on one piece and inference on the other
  **without halving the sample**. That is the answer to 6.1's stated cost —
  "`"split"` gives coverage at some loss of precision" — for Gaussian-like
  responses. It also has a known failure mode the discussion papers document:
  the noise scale is a tuning parameter with no data-driven default. The
  backlog should present it as the second option, not the first.
  **→ extends** · BACKLOG 6.1

**11.35+ Neufeld, Dharamshi, Gao & Witten (2024).** Data thinning for
convolution-closed distributions. *Journal of Machine Learning Research*,
25, paper 23-0446. jmlr.org/papers/v25/23-0446.html *(publisher page)*
  The count-data counterpart: a Poisson or negative-binomial observation is
  thinned into independent pieces exactly, with no tuning parameter. For the
  zero-inflated and count responses the package now supports via
  `family =` (1.1), this is the sample-splitting alternative that costs no
  precision at all, and it composes with 6.1's spatial half-split rather than
  replacing it.
  **→ extends** · BACKLOG 1.1, 6.1

**11.36+ García Rasines & Young (2023).** Splitting strategies for
post-selection inference. *Biometrika*, 110(3), 597–614.
DOI: 10.1093/biomet/asac070
  Compares randomised, deterministic and data-fission-style splits for
  post-selection inference and characterises when each is efficient. This is
  the paper that lets 6.1 choose between its own spatial half-split and the
  two entries above on stated grounds rather than by taste — in particular it
  shows the spatial case is unusual because a *contiguous* half-split is not
  an exchangeable one, which is both why it is needed and why its efficiency
  is lower than the i.i.d. theory suggests.
  **→ informs** · BACKLOG 6.1

---

## Watch list — preprints and unverified

Surfaced by the search, relevant, but either not yet peer-reviewed or not
verifiable through the metadata services during this pass. Not counted in the
thirty above.

* **Spatial conformal inference through localized quantile regression**
  (arXiv:2412.01098, Dec 2024). Spatial conformal intervals with local
  calibration; a candidate methods paper for the uncertainty side of 10.18+.
  I could not confirm a journal version — an *Annals of Statistics* 53(4)
  record that a search suggested turned out, on checking the DOI, to be a
  different paper (Liang & Foygel Barber 2025). Treat as preprint.
* **Assessing the performance of spatial cross-validation** (arXiv:2303.07334).
  Simulation comparison of spatial CV schemes behind the `spatialsample`
  package; tidymodels-aligned, so its findings on buffered blocks transfer to
  `make_folds(buffer = )`. No journal version found.
* **Demystifying spatial confounding** (arXiv:2309.16861) and **Robust spatial
  confounding adjustment via basis voting** (arXiv:2510.22464). Both continue
  the §7+ argument; both preprints.
* **A robust nonparametric test for spatial isotropy in lattice data**
  (arXiv:2605.18030, May 2026). Lattice only; noted for 4.3.
* **Effective sample size for functional spatial data** (arXiv:2601.20812).
  Out of scope; noted because it extends Ferrer & Vallejos.

## Method

Twelve topic searches mirroring the bibliography's sections, each targeting
the leading journals for that sub-literature and work from 2023 onward.
Candidates were matched against the 418 existing entries by first-author
surname and year and, where that failed, by title words, before any were kept;
several I expected to be absent were present under an online-first year (Lark,
Cullis & Welham as 2005; Khan & Calder as 2020; da Silva & Fotheringham as
2015; Wolf, Oshan & Fotheringham as 2017; Bolin & Wallin as 2025), which is
worth knowing when merging. DOIs were resolved individually through OpenAlex;
where the service rate-limited, the publisher or repository page was used and
the entry is marked. One search-suggested DOI was found to point at the wrong
paper and is recorded on the watch list rather than silently dropped.
