# convergence_ok follows Stan's own diagnostics.
#
# The verdict cut R-hat at 1.05 and judged effective sample size as a share
# of the draws (below 0.1 failed).  That share is 100 per chain only at the
# default 4 x 2000: a fit of two chains and 300 post-warmup draws each, with
# bulk and tail ESS of 130 and 137 and 66 parameters above R-hat 1.01, passed
# while rstan raised both of its ESS warnings on it.  Tree-depth saturation
# and E-BFMI were not read, and a failed verdict was a log line only.
#
# R-hat has two tiers.  One cut at 1.01 over every parameter fails fits that
# have converged: this model has about 540 parameters, and a default fit of
# four chains and 2,000 iterations whose smallest effective sample size was
# 1,278 (400 required) had one basis weight at 1.011.  So above 1.05 fails,
# above 1.01 fails beside too few effective draws for that parameter, and
# above 1.01 with enough is a caution.
#
# No Stan: brms::brm() returns a stand-in and the brms accessors the check
# reads are replaced, so the rules run on vectors chosen for each case.

.fbc_pts <- function(n = 120, seed = 1) {
  set.seed(seed)
  d <- sf::st_as_sf(
    data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
               a = stats::rnorm(n)),
    coords = c("x", "y"), crs = 32632)
  d$z <- 2 * d$a + stats::rnorm(n)
  d
}

.fbc_np <- function(n_div = 0, treedepth = rep(5, 40), energy = NULL,
                    chains = 2L) {
  per <- length(treedepth) / chains
  it  <- rep(seq_len(per), times = chains)
  ch  <- rep(seq_len(chains), each = per)
  if (is.null(energy)) {
    set.seed(99); energy <- stats::rnorm(length(treedepth))
  }
  rbind(
    data.frame(Chain = ch, Iteration = it, Parameter = "divergent__",
               Value = c(rep(1, n_div), rep(0, length(treedepth) - n_div))),
    data.frame(Chain = ch, Iteration = it, Parameter = "treedepth__",
               Value = treedepth),
    data.frame(Chain = ch, Iteration = it, Parameter = "energy__",
               Value = energy))
}

# Fit with the accessors returning what each case needs.  `...` reaches the
# fitter (chains, iter, control).
.fbc_fit <- function(rhat, neff, np = .fbc_np(), ...) {
  d <- .fbc_pts()
  local_mocked_bindings(
    brm = function(...) structure(list(), class = "brmsfit"),
    nuts_params = function(...) np,
    rhat = function(...) rhat,
    neff_ratio = function(...) neff,
    .package = "brms")
  local_mocked_bindings(estimate_sac_range = function(...) 5000,
                        .package = "spatialkit")
  w <- list()
  fit <- withCallingHandlers(
    suppressMessages(fit_bayesian_spatial_model(d, "z", "a", gp_k = 5,
                                                compute_loo = FALSE, ...)),
    warning = function(c) { w[[length(w) + 1L]] <<- c; invokeRestart("muffleWarning") })
  conv <- Filter(function(c) inherits(c, "spatialkit_convergence"), w)
  list(fit = fit, conv = conv)
}

.fbc_pars <- c("b_Intercept", "b_a", "sdgp_gp..x..y", "lscale_gp..x..y..x",
               "sigma", "lprior")
.fbc_vec  <- function(v) stats::setNames(rep_len(v, length(.fbc_pars)), .fbc_pars)


test_that("the effective-sample-size rule is 100 per chain, not a share of the draws", {
  low <- spatialkit:::.bayes_low_ess
  r <- c(a = 0.15, b = 0.5, c = NaN)
  # Two chains of 500: 0.15 of 1000 draws is 150, under the 200 required.
  # The share rule (below 0.1 fails) passed it.
  expect_identical(low(r, n_draws = 1000, n_chains = 2), "a")
  # Four chains of 10,000: 0.05 of 40,000 draws is 2,000, far above 400.
  # The share rule failed it.
  expect_identical(low(c(a = 0.05, b = 0.5), n_draws = 40000, n_chains = 4),
                   character(0))
  # At the default 4 x 1000 post-warmup draws the two rules agree.
  expect_identical(low(c(a = 0.09, b = 0.11), n_draws = 4000, n_chains = 4), "a")
  # Unknown draw or chain counts flag nothing; they do not error.
  expect_identical(low(r, n_draws = NA_real_, n_chains = 2), character(0))
  expect_identical(spatialkit:::.bayes_ess_per_chain, 100)
  expect_identical(spatialkit:::.bayes_rhat_max, 1.05)
  expect_identical(spatialkit:::.bayes_rhat_caution, 1.01)
  expect_identical(spatialkit:::.bayes_ebfmi_min, 0.3)
})

test_that(".bayes_rhat_tiers() fails above 1.05, and above 1.01 only beside too few effective draws", {
  tiers <- spatialkit:::.bayes_rhat_tiers
  rh <- c(a = 1.002, b = 1.011, c = 1.03, d = 1.06, e = NaN, f = 1.02)
  ne <- c(a = 0.5,   b = 0.32,  c = 0.05, d = 0.9,  e = NaN, f = 0.4)
  # Four chains of 1000: 400 effective draws required.  b has 1280 (caution),
  # c has 200 (fails), d is above 1.05 whatever it has, f has 1600 (caution).
  got <- tiers(rh, ne, n_draws = 4000, n_chains = 4)
  expect_identical(got$failed, c("c", "d"))
  expect_identical(got$caution, c("b", "f"))
  # A constant parameter (NaN) is in neither.
  expect_false("e" %in% c(got$failed, got$caution))
  # Nothing above 1.01: both empty.
  expect_identical(tiers(c(a = 1.001, b = 1.009), ne, 4000, 4),
                   list(failed = character(0), caution = character(0)))
  # An effective sample size that cannot be read does not excuse the R-hat:
  # no ratios at all, a parameter missing from them, unknown draw counts.
  expect_identical(tiers(rh, NULL, 4000, 4)$failed, c("b", "c", "d", "f"))
  expect_identical(tiers(rh, ne[c("a", "c", "d")], 4000, 4)$failed,
                   c("b", "c", "d", "f"))
  expect_identical(tiers(rh, ne, NA_real_, 4)$failed, c("b", "c", "d", "f"))
  expect_length(tiers(rh, NULL, 4000, 4)$caution, 0L)
  # The requirement follows the chains: at two chains of 500, c's 50 draws
  # are short and b's 320 are enough.
  got2 <- tiers(rh, ne, n_draws = 1000, n_chains = 2)
  expect_identical(got2$failed, c("c", "d"))
  expect_identical(got2$caution, c("b", "f"))
})

test_that(".bayes_convergence_advice() gives the remedy for the checks that failed", {
  adv <- spatialkit:::.bayes_convergence_advice
  # Divergences alone: adapt_delta, and no advice to run longer.
  d <- adv("divergent", iter = 2000, chains = 4, adapt_delta = 0.9)
  expect_match(d, "raise adapt_delta in `control` (now 0.90)", fixed = TRUE)
  expect_match(d, "more iterations do not remove them", fixed = TRUE)
  expect_false(grepl("Run more iterations", d, fixed = TRUE))
  # R-hat and effective sample size: more iterations, and nothing on
  # adapt_delta.
  r <- adv(c("rhat", "ess"), iter = 1000, chains = 2, adapt_delta = 0.9)
  expect_identical(r, paste0("Run more iterations (iter = 1000, chains = 2 now) ",
                             "to bring R-hat down and raise the effective sample size."))
  expect_identical(adv("ess", iter = 1000, chains = 2),
                   "Run more iterations (iter = 1000, chains = 2 now) to raise the effective sample size.")
  # E-BFMI is neither.
  e <- adv("ebfmi", iter = 2000, chains = 4)
  expect_match(e, "A low E-BFMI comes from the shape of the posterior", fixed = TRUE)
  expect_false(grepl("Run more iterations|adapt_delta", e))
  # All of them, in one string, each once.
  all4 <- adv(c("divergent", "ebfmi", "rhat", "ess"), 2000, 4, 0.95)
  expect_match(all4, "(now 0.95)", fixed = TRUE)
  expect_match(all4, "Run more iterations (iter = 2000, chains = 4 now)", fixed = TRUE)
  expect_match(all4, "E-BFMI", fixed = TRUE)
})

test_that(".bayes_ebfmi() is the energy fraction per chain", {
  set.seed(3)
  e1 <- stats::rnorm(200)                       # independent energies: about 2
  e2 <- cumsum(stats::rnorm(200, sd = 0.05))    # a slow random walk: near 0
  np <- data.frame(Chain = rep(1:2, each = 200), Iteration = rep(1:200, 2),
                   Parameter = "energy__", Value = c(e1, e2))
  got <- spatialkit:::.bayes_ebfmi(np)
  expect_named(got, c("1", "2"))
  expect_equal(got[["1"]], (sum(diff(e1)^2) / 200) / stats::var(e1))
  expect_gt(got[["1"]], 1)
  expect_lt(got[["2"]], 0.3)
  # No energies recorded: nothing to report.
  expect_length(spatialkit:::.bayes_ebfmi(np[0, ]), 0L)
  expect_length(spatialkit:::.bayes_ebfmi(data.frame(Parameter = "divergent__",
                                                     Value = 0)), 0L)
})

test_that("a clean fit passes and raises nothing", {
  skip_if_not_installed("brms")
  r <- .fbc_fit(rhat = .fbc_vec(1.002), neff = .fbc_vec(0.7))
  expect_true(isTRUE(r$fit$info$convergence_ok))
  expect_length(r$conv, 0L)
  cd <- r$fit$info$convergence_diagnostics
  expect_identical(cd$n_divergent, 0)
  expect_identical(cd$n_max_treedepth, 0L)
  expect_length(cd$rhat_failed, 0L)
  expect_length(cd$neff_failed, 0L)
  expect_length(cd$ess_failed, 0L)
  # 4 chains x (2000 - 1000) draws by default.
  expect_identical(cd$ess_required, 400)
  expect_equal(cd$min_ess, 0.7 * 4000)
  expect_named(cd$ebfmi, c("1", "2"))
})

test_that("R-hat above 1.05 fails the verdict and is raised with the parameters named", {
  skip_if_not_installed("brms")
  rh <- .fbc_vec(1.002)
  rh[["sdgp_gp..x..y"]] <- 1.07
  rh[["lprior"]] <- NaN                         # constant parameter: not a failure
  r <- .fbc_fit(rhat = rh, neff = .fbc_vec(0.7))
  expect_identical(r$fit$info$convergence_ok, FALSE)
  cd <- r$fit$info$convergence_diagnostics
  expect_identical(names(cd$rhat_failed), "sdgp_gp..x..y")
  expect_length(cd$rhat_caution, 0L)
  expect_equal(cd$max_rhat, 1.07)
  # One R warning, classed, where there was a log line.
  expect_length(r$conv, 1L)
  msg <- conditionMessage(r$conv[[1]])
  expect_match(msg, paste0("1 parameter(s) fail on R-hat (max 1.070): 1 above 1.05, ",
                           "0 above 1.01 with too few effective draws: sdgp_gp..x..y"),
               fixed = TRUE)
  expect_match(msg, "$info$convergence_ok is FALSE", fixed = TRUE)
  # R-hat is answered with more iterations, not with adapt_delta.
  expect_match(msg, "Run more iterations (iter = 2000, chains = 4 now) to bring R-hat down.",
               fixed = TRUE)
  expect_false(grepl("adapt_delta", msg, fixed = TRUE))
  txt <- paste(utils::capture.output(print(r$fit)), collapse = "\n")
  expect_match(txt, "Convergence warnings present", fixed = TRUE)
})

test_that("R-hat between 1.01 and 1.05 with enough effective draws is a caution, not a failure", {
  skip_if_not_installed("brms")
  # The fit the rule was changed for: about 540 parameters at the default
  # four chains of 1000 draws, one basis weight at 1.011, and a smallest
  # effective sample size of 1,278 against 400 required.  It was failed and
  # told to run more iterations.
  pars <- c(.fbc_pars, sprintf("zgp_gp..x..y[%d]", 1:529))
  rh <- stats::setNames(rep(1.003, length(pars)), pars)
  rh[["zgp_gp..x..y[63]"]] <- 1.011
  rh[["lprior"]] <- NaN
  ne <- stats::setNames(rep(0.6, length(pars)), pars)
  ne[["zgp_gp..x..y[63]"]] <- 1278 / 4000
  logged <- capture_spatialkit_log(r <- .fbc_fit(rhat = rh, neff = ne),
                                   level = logger::WARN)
  expect_true(isTRUE(r$fit$info$convergence_ok))
  expect_length(r$conv, 0L)
  cd <- r$fit$info$convergence_diagnostics
  expect_length(cd$rhat_failed, 0L)
  expect_identical(names(cd$rhat_caution), "zgp_gp..x..y[63]")
  expect_equal(unname(cd$rhat_caution), 1.011)
  expect_equal(cd$max_rhat, 1.011)
  # Recorded, logged and printed; not flagged.
  expect_true(log_has(logged, "1 of 534 parameter\\(s\\) have R-hat between 1.01 and 1.05 \\(max 1.011\\)"))
  expect_true(log_has(logged, "at least the 400 effective draws required"))
  txt <- paste(utils::capture.output(print(r$fit)), collapse = "\n")
  expect_match(txt, "Convergence: passed; 1 parameter(s) with R-hat between 1.01 and 1.05",
               fixed = TRUE)
  expect_false(grepl("Convergence warnings present", txt, fixed = TRUE))
  # print() on a summary says the same.  (Built by hand: the stand-in engine
  # has no fitted values for summary() to score.)
  smry <- structure(list(class = "bayesian_fit", formula = r$fit$formula, n = r$fit$n,
                         response_var = "z", predictor_vars = "a",
                         info = r$fit$info[setdiff(names(r$fit$info), ".cache")],
                         in_sample = data.frame(n = r$fit$n, RMSE = 1, MAE = 1, R2 = 0.5)),
                    class = "summary.spatial_fit")
  stx <- paste(utils::capture.output(print(smry)), collapse = "\n")
  expect_match(stx, "Convergence: passed; 1 parameter(s) with R-hat between 1.01 and 1.05",
               fixed = TRUE)
  expect_false(grepl("Convergence warnings present", stx, fixed = TRUE))

  # The same R-hat beside too few effective draws for that parameter fails,
  # on both counts.
  ne2 <- ne; ne2[["zgp_gp..x..y[63]"]] <- 300 / 4000
  r2 <- .fbc_fit(rhat = rh, neff = ne2)
  expect_identical(r2$fit$info$convergence_ok, FALSE)
  expect_identical(names(r2$fit$info$convergence_diagnostics$rhat_failed),
                   "zgp_gp..x..y[63]")
  expect_length(r2$fit$info$convergence_diagnostics$rhat_caution, 0L)
  m2 <- conditionMessage(r2$conv[[1]])
  expect_match(m2, paste0("1 parameter(s) fail on R-hat (max 1.011): 0 above 1.05, ",
                          "1 above 1.01 with too few effective draws: zgp_gp..x..y[63]"),
               fixed = TRUE)
  expect_match(m2, "effective sample size under 400", fixed = TRUE)
  expect_match(m2, "to bring R-hat down and raise the effective sample size", fixed = TRUE)

  # And it fails when the effective sample size cannot be read at all.
  r3 <- .fbc_fit(rhat = rh, neff = NULL)
  expect_identical(r3$fit$info$convergence_ok, FALSE)
  expect_identical(names(r3$fit$info$convergence_diagnostics$rhat_failed),
                   "zgp_gp..x..y[63]")
  expect_match(conditionMessage(r3$conv[[1]]),
               "1 above 1.01 with an effective sample size that could not be read: zgp_gp..x..y[63]",
               fixed = TRUE)

  # A clean fit prints no convergence line at all.
  r4 <- .fbc_fit(rhat = .fbc_vec(1.002), neff = .fbc_vec(0.7))
  expect_false(grepl("Convergence", paste(utils::capture.output(print(r4$fit)),
                                          collapse = "\n"), fixed = TRUE))
})

test_that("effective sample size is judged per chain in the fit", {
  skip_if_not_installed("brms")
  ne <- .fbc_vec(0.7)
  ne[["lscale_gp..x..y..x"]] <- 0.15
  # chains = 2, iter = 1000: 1000 post-warmup draws, so 150 effective draws
  # against 200 required.  convergence_ok was TRUE.
  r <- .fbc_fit(rhat = .fbc_vec(1.002), neff = ne, chains = 2, iter = 1000)
  expect_identical(r$fit$info$convergence_ok, FALSE)
  cd <- r$fit$info$convergence_diagnostics
  expect_identical(names(cd$neff_failed), "lscale_gp..x..y..x")
  expect_equal(unname(cd$ess_failed), 150)
  expect_identical(cd$ess_required, 200)
  expect_equal(cd$min_ess, 150)
  expect_length(r$conv, 1L)
  expect_match(conditionMessage(r$conv[[1]]),
               "1 parameter(s) with bulk or tail effective sample size under 200 (100 per chain x 2 chain(s); min 150): lscale_gp..x..y..x",
               fixed = TRUE)

  # The same ratio is enough at the default four chains of 1000 draws: 600.
  r4 <- .fbc_fit(rhat = .fbc_vec(1.002), neff = ne)
  expect_true(isTRUE(r4$fit$info$convergence_ok))
  # And a ratio the share rule failed passes when the chains are long.
  ne2 <- .fbc_vec(0.7); ne2[["sigma"]] <- 0.05
  r5 <- .fbc_fit(rhat = .fbc_vec(1.002), neff = ne2, chains = 4, iter = 20000)
  expect_true(isTRUE(r5$fit$info$convergence_ok))
  expect_equal(r5$fit$info$convergence_diagnostics$min_ess, 0.05 * 40000)
})

test_that("divergences and low E-BFMI fail; tree-depth saturation is recorded only", {
  skip_if_not_installed("brms")
  r <- .fbc_fit(rhat = .fbc_vec(1.002), neff = .fbc_vec(0.7),
                np = .fbc_np(n_div = 3))
  expect_identical(r$fit$info$convergence_ok, FALSE)
  msg <- conditionMessage(r$conv[[1]])
  expect_match(msg, "3 divergent transition(s) after warmup", fixed = TRUE)
  # The remedy is adapt_delta.  The warning used to close with "Run more
  # iterations" whatever had failed, and more iterations do not remove
  # divergent transitions.
  expect_match(msg, "raise adapt_delta in `control` (now 0.90); more iterations do not remove them",
               fixed = TRUE)
  expect_false(grepl("Run more iterations", msg, fixed = TRUE))
  # A user's adapt_delta is the one reported.
  r <- .fbc_fit(rhat = .fbc_vec(1.002), neff = .fbc_vec(0.7),
                np = .fbc_np(n_div = 2), control = list(adapt_delta = 0.95))
  expect_match(conditionMessage(r$conv[[1]]), "(now 0.95)", fixed = TRUE)

  set.seed(5)
  walk <- c(stats::rnorm(20), cumsum(stats::rnorm(20, sd = 0.02)))
  r <- .fbc_fit(rhat = .fbc_vec(1.002), neff = .fbc_vec(0.7),
                np = .fbc_np(energy = walk))
  expect_identical(r$fit$info$convergence_ok, FALSE)
  expect_lt(r$fit$info$convergence_diagnostics$ebfmi[["2"]], 0.3)
  expect_gt(r$fit$info$convergence_diagnostics$ebfmi[["1"]], 0.3)
  expect_match(conditionMessage(r$conv[[1]]), "E-BFMI below 0.3 in chain(s) 2",
               fixed = TRUE)
  expect_match(conditionMessage(r$conv[[1]]),
               "A low E-BFMI comes from the shape of the posterior", fixed = TRUE)
  expect_false(grepl("Run more iterations", conditionMessage(r$conv[[1]]), fixed = TRUE))

  # Transitions at the tree-depth limit (12, the package default): counted
  # and logged, verdict unchanged.
  td <- rep(5, 40); td[1:7] <- 12
  logged <- capture_spatialkit_log(
    r <- .fbc_fit(rhat = .fbc_vec(1.002), neff = .fbc_vec(0.7),
                  np = .fbc_np(treedepth = td)),
    level = logger::WARN)
  expect_true(isTRUE(r$fit$info$convergence_ok))
  expect_identical(r$fit$info$convergence_diagnostics$n_max_treedepth, 7L)
  expect_length(r$conv, 0L)
  expect_true(log_has(logged, "7 transition\\(s\\) stopped at max_treedepth = 12"))
})

test_that("several failures arrive as one warning", {
  skip_if_not_installed("brms")
  rh <- .fbc_vec(1.002); rh[c("sigma", "b_a")] <- c(1.08, 1.02)
  ne <- .fbc_vec(0.7);   ne[["sigma"]] <- 0.01
  r <- .fbc_fit(rhat = rh, neff = ne, np = .fbc_np(n_div = 1))
  expect_length(r$conv, 1L)
  msg <- conditionMessage(r$conv[[1]])
  expect_match(msg, "1 divergent transition(s)", fixed = TRUE)
  # sigma is above 1.05.  b_a is at 1.02 with 2,800 effective draws: a
  # caution, counted in the message and not in the failures.
  expect_match(msg, paste0("1 parameter(s) fail on R-hat (max 1.080): 1 above 1.05, ",
                           "0 above 1.01 with too few effective draws (1 more above ",
                           "1.01 have enough and are not counted): sigma"),
               fixed = TRUE)
  expect_identical(names(r$fit$info$convergence_diagnostics$rhat_caution), "b_a")
  expect_match(msg, "effective sample size under 400", fixed = TRUE)
  # Each failed check gets its own remedy, once.
  expect_match(msg, "raise adapt_delta in `control` (now 0.90)", fixed = TRUE)
  expect_match(msg, "Run more iterations (iter = 2000, chains = 4 now) to bring R-hat down and raise the effective sample size.",
               fixed = TRUE)
  # It can be caught by class.
  caught <- tryCatch(
    {
      d <- .fbc_pts()
      local_mocked_bindings(
        brm = function(...) structure(list(), class = "brmsfit"),
        nuts_params = function(...) .fbc_np(n_div = 2),
        rhat = function(...) .fbc_vec(1.002),
        neff_ratio = function(...) .fbc_vec(0.7),
        .package = "brms")
      local_mocked_bindings(estimate_sac_range = function(...) 5000,
                            .package = "spatialkit")
      suppressMessages(fit_bayesian_spatial_model(d, "z", "a", gp_k = 5,
                                                  compute_loo = FALSE))
      "no condition"
    },
    spatialkit_convergence = function(c) "caught")
  expect_identical(caught, "caught")
})
