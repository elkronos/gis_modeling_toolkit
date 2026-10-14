# tests/testthat/test-core-seeding-starts.R
# ---------------------------------------------------------------------------
# The seeding functions' partition arguments: `kmeans_nstart = NULL` and
# `nstart = NULL` (the kind's own number of starts), which partition
# arguments a call is taken to name beside scored centres and beside
# method = "random" or "provided", and the same on copies of the two
# functions whose own kind is the density kind, the edit a later step makes
# to their formals.
#
# Every call names `partition` and the number of starts where its result
# depends on them.  A test of what a call that names no kind is told is made
# on copies of the functions whose own kind is written into their formals,
# the k-means kind or the density kind, so that no test here reads a default.
# ---------------------------------------------------------------------------

.ss_pts <- function()
  partition_points_sf(partition_clustered_xy(n_bg = 100L, n_cl = 300L, seed = 29L))
.ss_bnd <- function(pts) sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(pts)))

# The four profiles of scored centres, and a fifth of the k-means kind from
# random starts: two levels, two starts, no response.
.ss_store <- new.env(parent = emptyenv())
.ss_profiles <- function() {
  if (is.null(.ss_store$p)) {
    pts <- .ss_pts()
    prof <- function(...)
      resolution_profile(pts, levels = c(6L, 14L), nstart = 2, estimate = NULL,
                         obs_per_cell = NULL, sac = NULL, ...)
    .ss_store$p <- list(
      kmeans  = prof(partition = "kmeans"),
      density = prof(partition = "density"),
      power   = prof(partition = "density", density_power = 0.5),
      k_nn    = prof(partition = "density", density_k = 12),
      random  = prof(partition = "kmeans", init = "random"))
  }
  .ss_store$p
}
# A selection at 14 cells that carries the level's centres and nothing a
# warning about the count could be raised for.
.ss_sel <- function(p, L = 14L) structure(
  list(best = L, seeds = spatialkit:::.profile_seeds(p, L), criterion = "cp",
       edge = NA_character_, structured = TRUE), class = "resolution_selection")

# A copy of a seeding function whose own kind is the one named: the k-means
# kind, or the density kind.
.ss_own <- function(f, kind) {
  formals(f)$partition <- if (identical(kind, "density")) quote(c("density", "kmeans"))
                          else quote(c("kmeans", "density"))
  f
}
.ss_text <- function(r) vapply(r$warnings, conditionMessage, character(1))


# ---- NULL starts ---------------------------------------------------------------

test_that("kmeans_nstart = NULL is 10 starts for k-means and 3 for the density kind", {
  pts <- .ss_pts()
  seeds <- function(...) partition_quietly(
    get_voronoi_seeds(method = "kmeans", n = 25, sample_points = pts, set_seed = 1, ...))
  # The layer and the seed are ones on which 3 starts and 10 starts keep
  # different centres, for either kind and either starting rule: a NULL read
  # as the other kind's number would show.
  for (case in list(list(partition = "kmeans"),
                    list(partition = "kmeans", init = "kmeans++"),
                    list(partition = "density"))) {
    own   <- if (identical(case$partition, "density")) 3 else 10
    other <- if (identical(case$partition, "density")) 10 else 3
    with_null  <- do.call(seeds, c(case, list(kmeans_nstart = NULL)))
    with_own   <- do.call(seeds, c(case, list(kmeans_nstart = own)))
    with_other <- do.call(seeds, c(case, list(kmeans_nstart = other)))
    info <- paste(unlist(case), collapse = ", ")
    expect_identical(with_null, with_own, info = info)
    expect_identical(attr(with_null, "kmeans")$nstart, as.integer(own), info = info)
    expect_false(isTRUE(all.equal(partition_xy(with_null), partition_xy(with_other))),
                 info = info)
  }
})

test_that("nstart = NULL in voronoi_seeds_kmeans() is 10 starts for k-means and 3 for the density kind", {
  pts <- .ss_pts()
  seeds <- function(...) partition_quietly(voronoi_seeds_kmeans(pts, 25, set_seed = 1, ...))
  for (case in list(list(partition = "kmeans"),
                    list(partition = "kmeans", init = "kmeans++"),
                    list(partition = "density"))) {
    own   <- if (identical(case$partition, "density")) 3 else 10
    other <- if (identical(case$partition, "density")) 10 else 3
    with_null  <- do.call(seeds, c(case, list(nstart = NULL)))
    with_own   <- do.call(seeds, c(case, list(nstart = own)))
    with_other <- do.call(seeds, c(case, list(nstart = other)))
    info <- paste(unlist(case), collapse = ", ")
    # The function keeps no record of the number, so the identity is the test.
    expect_identical(with_null, with_own, info = info)
    expect_null(attr(with_null, "kmeans"), info = info)
    expect_false("nstart" %in% names(attr(with_null, "partition")), info = info)
    expect_false(isTRUE(all.equal(partition_xy(with_null), partition_xy(with_other))),
                 info = info)
  }
})

test_that("a number of starts that is refused is refused in words that name NULL", {
  pts <- .ss_pts()
  # voronoi_seeds_kmeans() checks the number on every path.
  for (kind in c("kmeans", "density")) {
    expect_error(voronoi_seeds_kmeans(pts, 5, partition = kind, nstart = 0),
                 paste0("voronoi_seeds_kmeans(): `nstart` must be NULL or a single positive ",
                        "number of k-means restarts and at least 1; got 0."), fixed = TRUE)
    expect_error(voronoi_seeds_kmeans(pts, 5, partition = kind, nstart = "three"),
                 paste0("voronoi_seeds_kmeans(): `nstart` must be NULL or a single positive ",
                        "number of k-means restarts; got character."), fixed = TRUE)
    expect_error(voronoi_seeds_kmeans(pts, 5, partition = kind, nstart = c(2, 3)),
                 "`nstart` must be NULL or a single positive number of k-means restarts; got numeric of length 2.",
                 fixed = TRUE)
    expect_error(voronoi_seeds_kmeans(pts, 5, partition = kind, nstart = NA_real_),
                 "`nstart` must be NULL or a single positive number of k-means restarts; got NA.",
                 fixed = TRUE)
  }
  # get_voronoi_seeds() checks it on the paths that are not the run it has
  # always made: the density kind, and k-means++ starts.
  for (extra in list(list(partition = "density"), list(partition = "kmeans", init = "kmeans++"))) {
    f <- function(...) do.call(get_voronoi_seeds,
                               c(list(method = "kmeans", n = 5, sample_points = pts, set_seed = 1),
                                 extra, list(...)))
    expect_error(f(kmeans_nstart = 0),
                 paste0("get_voronoi_seeds(): `kmeans_nstart` must be NULL or a single positive ",
                        "whole number of starts and at least 1; got 0."), fixed = TRUE)
    expect_error(f(kmeans_nstart = 2.5),
                 paste0("get_voronoi_seeds(): `kmeans_nstart` must be NULL or a single positive ",
                        "whole number of starts; got 2.5."), fixed = TRUE)
    expect_error(f(kmeans_nstart = NA),
                 "`kmeans_nstart` must be NULL or a single positive whole number of starts; got logical.",
                 fixed = TRUE)
    # The iteration limit takes no NULL, and its text has none.
    expect_error(f(kmeans_iter = 0),
                 paste0("get_voronoi_seeds(): `kmeans_iter` must be a single positive whole ",
                        "number of iterations and at least 1; got 0."), fixed = TRUE)
  }
})

test_that("method = \"random\" and \"provided\" do not read the number of starts, NULL or not", {
  pts <- .ss_pts(); bnd <- .ss_bnd(pts)
  base <- get_voronoi_seeds(bnd, method = "random", n = 5, set_seed = 2)
  expect_identical(get_voronoi_seeds(bnd, method = "random", n = 5, set_seed = 2,
                                     kmeans_nstart = NULL), base)
  expect_identical(get_voronoi_seeds(bnd, method = "random", n = 5, set_seed = 2,
                                     kmeans_nstart = 0), base)
  given <- get_voronoi_seeds(method = "provided", seeds = pts[1:4, ])
  expect_identical(get_voronoi_seeds(method = "provided", seeds = pts[1:4, ],
                                     kmeans_nstart = NULL), given)
})

test_that("a call on a path that never reads the number of starts does not evaluate it", {
  pts <- .ss_pts(); bnd <- .ss_bnd(pts)
  sel <- .ss_sel(.ss_profiles()$kmeans)
  # An argument that stops when it is evaluated, and a wrapper whose own
  # argument for the number of starts the caller left out.
  never <- quote(stop("the number of starts was evaluated"))
  gvs_via <- function(f) function(starts, ...) f(..., kmeans_nstart = starts)
  vsk_via <- function(f) function(starts, ...) f(..., nstart = starts)

  for (kind in c("kmeans", "density")) {
    gvs <- .ss_own(get_voronoi_seeds, kind)
    vsk <- .ss_own(voronoi_seeds_kmeans, kind)
    # method = "random" and "provided".
    base <- gvs(bnd, method = "random", n = 5, set_seed = 2)
    expect_identical(eval(bquote(gvs(bnd, method = "random", n = 5, set_seed = 2,
                                     kmeans_nstart = .(never)))), base, info = kind)
    expect_identical(gvs_via(gvs)(, bnd, method = "random", n = 5, set_seed = 2), base,
                     info = kind)
    given <- gvs(method = "provided", seeds = pts[1:4, ])
    expect_identical(eval(bquote(gvs(method = "provided", seeds = pts[1:4, ],
                                     kmeans_nstart = .(never)))), given, info = kind)
    expect_identical(gvs_via(gvs)(, method = "provided", seeds = pts[1:4, ]), given, info = kind)
    # Scored centres from voronoi_seeds_kmeans(), which returns them before
    # its check of the number, whichever kind is in force.
    scored <- partition_quietly(vsk(pts, sel, partition = "kmeans"))
    expect_identical(nrow(scored), 14L)
    expect_identical(partition_quietly(eval(bquote(vsk(pts, sel, partition = "kmeans",
                                                       nstart = .(never))))), scored, info = kind)
    expect_identical(partition_quietly(vsk_via(vsk)(, pts, sel, partition = "kmeans")), scored,
                     info = kind)
  }

  # get_voronoi_seeds() beside scored centres: its own run of k-means from
  # random starts checks nothing, and the centres are returned unread.
  own_kmeans <- .ss_own(get_voronoi_seeds, "kmeans")
  scored <- own_kmeans(method = "kmeans", n = sel, sample_points = pts)
  expect_identical(eval(bquote(own_kmeans(method = "kmeans", n = sel, sample_points = pts,
                                          kmeans_nstart = .(never)))), scored)
  expect_identical(gvs_via(own_kmeans)(, method = "kmeans", n = sel, sample_points = pts), scored)
  # The density kind and k-means++ starts check the number before anything
  # else, scored centres or not, as they always have.
  for (args in list(list(partition = "density"), list(init = "kmeans++")))
    expect_error(partition_quietly(do.call(own_kmeans, c(
      list(method = "kmeans", n = sel, sample_points = pts, kmeans_nstart = never), args))),
      "the number of starts was evaluated", fixed = TRUE)
  expect_error(partition_quietly(eval(bquote(
    .ss_own(get_voronoi_seeds, "density")(method = "kmeans", n = sel, sample_points = pts,
                                          kmeans_nstart = .(never))))),
    "the number of starts was evaluated", fixed = TRUE)

  # A fresh run of k-means from random starts reads the number where the fit
  # first needs it and no sooner: after the checks of the other arguments,
  # and after the stream is seeded.
  expect_error(eval(bquote(own_kmeans(method = "kmeans", n = 5, set_seed = 1,
                                      kmeans_nstart = .(never)))),
               "k-means seeds require 'boundary'", fixed = TRUE)
  expect_error(eval(bquote(own_kmeans(method = "kmeans", sample_points = pts, set_seed = 1,
                                      kmeans_nstart = .(never)))),
               "argument 'n' is required", fixed = TRUE)
  expect_error(eval(bquote(own_kmeans(method = "kmeans", n = 5, sample_points = pts,
                                      set_seed = 1, kmeans_nstart = .(never)))),
               "the number of starts was evaluated", fixed = TRUE)
  drawn <- NULL
  with_draw <- own_kmeans(method = "kmeans", n = 5, sample_points = pts, set_seed = 1,
                          kmeans_nstart = { drawn <- stats::runif(1); 4 })
  expect_identical(drawn, withr::with_seed(1, stats::runif(1)))
  expect_identical(attr(with_draw, "kmeans")$nstart, 4L)
  # voronoi_seeds_kmeans() reads it at its check, after `k` is checked.
  own_vsk <- .ss_own(voronoi_seeds_kmeans, "kmeans")
  expect_error(eval(bquote(own_vsk(pts, 0, nstart = .(never)))),
               "`k` must be a single positive number of seeds", fixed = TRUE)
  expect_error(eval(bquote(own_vsk(pts, 5, nstart = .(never)))),
               "the number of starts was evaluated", fixed = TRUE)
})


# ---- The partition arguments a call names --------------------------------------

test_that("with the k-means kind as its own, a random or provided seeding is warned of a kind or a starting rule it names, and of no other", {
  pts <- .ss_pts(); bnd <- .ss_bnd(pts)
  gvs <- .ss_own(get_voronoi_seeds, "kmeans")
  random   <- function(...) gvs(bnd, method = "random", n = 5, set_seed = 2, ...)
  provided <- function(...) gvs(method = "provided", seeds = pts[1:4, ], ...)
  named <- "get_voronoi_seeds(): `partition` and `init` are read by method = \"kmeans\" only; method = \"%s\" ignores them."
  for (case in list(list(random, "random"), list(provided, "provided"))) {
    f <- case[[1L]]; method <- case[[2L]]
    base <- partition_warnings(f())
    expect_length(base$warnings, 0L)
    # Naming the kind the function has anyway asks for nothing it does not get.
    expect_length(partition_warnings(f(partition = "kmeans"))$warnings, 0L)
    for (args in list(list(partition = "density"), list(init = "kmeans++"),
                      list(init = "random"), list(partition = "kmeans", init = "random"))) {
      r <- partition_warnings(do.call(f, args))
      expect_identical(.ss_text(r), sprintf(named, method))
      expect_identical(r$value, base$value)
    }
    # A density setting with no kind named: the present warning that the
    # k-means kind does not read it.
    r <- partition_warnings(f(density_power = 0.5))
    expect_identical(.ss_text(r), paste0(
      "get_voronoi_seeds(): `density_power` is read by partition = \"density\" only, and ",
      "partition = \"kmeans\" ignores it."))
    r <- partition_warnings(f(density_power = 0.5, density_k = 12))
    expect_identical(.ss_text(r), paste0(
      "get_voronoi_seeds(): `density_power` and `density_k` are read by partition = ",
      "\"density\" only, and partition = \"kmeans\" ignores them."))
    # With the density kind named too, one warning: the kind is not read.
    r <- partition_warnings(f(partition = "density", density_power = 0.5))
    expect_identical(.ss_text(r), sprintf(named, method))
  }
  # The functions themselves, whichever kind is their own: nothing named,
  # nothing said; the density kind or a starting rule named, the one warning.
  expect_length(partition_warnings(
    get_voronoi_seeds(bnd, method = "random", n = 5, set_seed = 2))$warnings, 0L)
  expect_length(partition_warnings(
    get_voronoi_seeds(method = "provided", seeds = pts[1:4, ]))$warnings, 0L)
  expect_identical(.ss_text(partition_warnings(
    get_voronoi_seeds(bnd, method = "random", n = 5, set_seed = 2, partition = "density"))),
    sprintf(named, "random"))
  expect_identical(.ss_text(partition_warnings(
    get_voronoi_seeds(method = "provided", seeds = pts[1:4, ], init = "kmeans++"))),
    sprintf(named, "provided"))
})

test_that("a selection passed as the count by a call that names no partition argument raises nothing", {
  pts <- .ss_pts()
  profiles <- .ss_profiles()
  for (nm in c("kmeans", "density", "power", "k_nn")) {
    sel <- .ss_sel(profiles[[nm]])
    cen <- attr(profiles[[nm]], "centres")[["14"]]
    for (r in list(partition_warnings(voronoi_seeds_kmeans(pts, sel)),
                   partition_warnings(get_voronoi_seeds(method = "kmeans", n = sel,
                                                        sample_points = pts)))) {
      expect_length(r$warnings, 0L)
      expect_equal(partition_xy(r$value), cen, tolerance = 1e-12)
      expect_identical(attr(r$value, "partition"),
                       attr(profiles[[nm]], "partition")[c("kind", "init", "density_power",
                                                           "density_k")])
    }
  }
  # The kind named beside centres of that kind built at another power: both
  # density settings are named with it, at their defaults.
  sel <- .ss_sel(profiles$power)
  for (r in list(partition_warnings(voronoi_seeds_kmeans(pts, sel, partition = "density")),
                 partition_warnings(get_voronoi_seeds(method = "kmeans", n = sel,
                                                      sample_points = pts,
                                                      partition = "density")))) {
    expect_length(partition_of_class(r, "spatialkit_partition_unused"), 1L)
    expect_length(r$warnings, 1L)
  }
})

# The rule as it was before a call's supplied density settings were read:
# both settings named whenever the density kind is in force.
.ss_rule_before <- function(seeds, part, init, partition_given, arg, caller) {
  rec <- attr(seeds, "partition", exact = TRUE)
  kind <- if (is.null(rec)) "kmeans" else rec$kind
  rule <- if (is.null(rec)) "kmeans++" else rec$init
  named <- if (identical(part$kind, "density")) c("density_power", "density_k")
           else part$unused
  differs <- (isTRUE(partition_given) && !identical(part$kind, kind)) ||
    (!is.null(init) && !identical(init, rule)) ||
    (identical(kind, "density") &&
       (("density_power" %in% named &&
           !isTRUE(all.equal(part$asked$density_power, rec$density_power))) ||
          ("density_k" %in% named &&
             !identical(as.integer(part$asked$density_k), as.integer(rec$density_k)))))
  if (differs) "spatialkit_partition_unused"
  else if (!identical(kind, "density") && length(part$unused)) "density setting unused"
  else "nothing"
}

test_that("with the k-means kind as the functions' own, every call is told what it was told before", {
  # Five kinds of profile by the 81 ways of leaving each of the four
  # arguments unset or setting it to one of two values: 405 calls.
  profiles <- .ss_profiles()
  values <- list(partition = list(NULL, "kmeans", "density"),
                 init = list(NULL, "kmeans++", "random"),
                 density_power = list(NULL, 0.5, 0.7),
                 density_k = list(NULL, 12, 4))
  grid <- expand.grid(partition = 1:3, init = 1:3, density_power = 1:3, density_k = 1:3)
  calls <- 0L
  for (nm in names(profiles)) {
    seeds <- spatialkit:::.profile_seeds(profiles[[nm]], 14L)
    for (i in seq_len(nrow(grid))) {
      a <- lapply(names(values), function(v) values[[v]][[grid[[v]][i]]])
      names(a) <- names(values)
      given <- !vapply(a, is.null, logical(1))
      kind  <- if (is.null(a$partition)) "kmeans" else a$partition      # the functions' own
      part  <- .seed_partition_args(kind, a$init, a$density_power %||% 1, a$density_k %||% 8,
                                    "f")
      r <- partition_warnings(.scored_partition_unused(
        seeds, part, part$init_given, given[["partition"]], arg = "n", caller = "f",
        density_given = names(given)[given & names(given) %in% c("density_power", "density_k")]))
      got <- if (!length(r$warnings)) "nothing"
             else if (inherits(r$warnings[[1L]], "spatialkit_partition_unused"))
               "spatialkit_partition_unused"
             else "density setting unused"
      expect_lte(length(r$warnings), 1L)
      expect_identical(got, .ss_rule_before(seeds, part, part$init_given, given[["partition"]],
                                            "n", "f"),
                       info = paste(nm, paste(names(a)[given], unlist(a[given]), collapse = " ")))
      calls <- calls + 1L
    }
  }
  expect_identical(calls, 405L)
})


# ---- The same with the functions' own kind flipped -----------------------------

test_that("with the density kind as the functions' own, a call that names nothing is told nothing", {
  pts <- .ss_pts()
  profiles <- .ss_profiles()
  gvs <- .ss_own(get_voronoi_seeds, "density")
  vsk <- .ss_own(voronoi_seeds_kmeans, "density")
  expect_identical(eval(formals(gvs)$partition), c("density", "kmeans"))
  expect_identical(eval(formals(.ss_own(gvs, "kmeans"))$partition), c("kmeans", "density"))
  via <- list(function(sel, ...) vsk(pts, sel, ...),
              function(sel, ...) gvs(method = "kmeans", n = sel, sample_points = pts, ...))
  unused <- function(r) length(partition_of_class(r, "spatialkit_partition_unused"))

  for (nm in c("kmeans", "density", "power", "k_nn")) {
    sel <- .ss_sel(profiles[[nm]])
    for (f in via) {
      r <- partition_warnings(f(sel))
      expect_length(r$warnings, 0L)
      # The scored centres, of the kind the profile was built with.
      expect_equal(partition_xy(r$value), attr(profiles[[nm]], "centres")[["14"]],
                   tolerance = 1e-12)
      expect_identical(attr(r$value, "partition"),
                       attr(profiles[[nm]], "partition")[c("kind", "init", "density_power",
                                                           "density_k")])
    }
  }
  # With the rule as it was, the two profiles built at other settings than
  # the defaults were warned about on a call that names nothing.
  part <- .seed_partition_args("density", NULL, 1, 8, "f")
  for (nm in c("power", "k_nn"))
    expect_identical(.ss_rule_before(spatialkit:::.profile_seeds(profiles[[nm]], 14L), part,
                                     NULL, FALSE, "n", "f"), "spatialkit_partition_unused")

  # Beside the profile built at power 0.5.
  sel <- .ss_sel(profiles$power)
  for (f in via) {
    expect_length(partition_warnings(f(sel, density_power = 0.5))$warnings, 0L)
    for (other in list(list(density_power = 0.7), list(density_k = 4))) {
      r <- partition_warnings(do.call(f, c(list(sel), other)))
      expect_identical(unused(r), 1L)
      expect_length(r$warnings, 1L)
      expect_match(conditionMessage(r$warnings[[1L]]),
                   "built with partition = \"density\" from k-means++ starts (power 0.5, 8 neighbours)",
                   fixed = TRUE)
    }
    # The neighbour count the profile has, supplied: nothing; and the kind
    # named names both settings at the values in force, as it always has.
    expect_length(partition_warnings(f(sel, density_k = 8))$warnings, 0L)
    expect_identical(unused(partition_warnings(f(sel, partition = "density"))), 1L)
    expect_length(partition_warnings(f(sel, partition = "density",
                                       density_power = 0.5))$warnings, 0L)
  }
  # Beside the profile built with 12 neighbours.
  sel <- .ss_sel(profiles$k_nn)
  for (f in via) {
    expect_length(partition_warnings(f(sel, density_k = 12))$warnings, 0L)
    expect_identical(unused(partition_warnings(f(sel, density_k = 8))), 1L)
    expect_length(partition_warnings(f(sel, density_power = 1))$warnings, 0L)
  }
  # Beside the k-means profile: a density setting supplied is read by
  # nothing, and the kind named as built is nothing to say.
  sel <- .ss_sel(profiles$kmeans)
  for (f in via) {
    r <- partition_warnings(f(sel, density_power = 0.5))
    expect_identical(unused(r), 1L)
    expect_length(r$warnings, 1L)
    expect_match(conditionMessage(r$warnings[[1L]]),
                 "built with partition = \"kmeans\" from k-means++ starts, and those are returned",
                 fixed = TRUE)
    # Supplied at its default counts as supplied.
    expect_identical(unused(partition_warnings(f(sel, density_k = 8))), 1L)
    expect_length(partition_warnings(f(sel, partition = "kmeans"))$warnings, 0L)
    expect_identical(unused(partition_warnings(f(sel, partition = "density"))), 1L)
  }
})

test_that("with the density kind as its own, a random or provided seeding is warned of what it names", {
  pts <- .ss_pts(); bnd <- .ss_bnd(pts)
  gvs <- .ss_own(get_voronoi_seeds, "density")
  random   <- function(...) gvs(bnd, method = "random", n = 5, set_seed = 2, ...)
  provided <- function(...) gvs(method = "provided", seeds = pts[1:4, ], ...)
  for (case in list(list(random, "random"), list(provided, "provided"))) {
    f <- case[[1L]]; method <- case[[2L]]
    base <- partition_warnings(f())
    expect_length(base$warnings, 0L)
    r <- partition_warnings(f(density_power = 0.5))
    expect_identical(.ss_text(r), sprintf(paste0(
      "get_voronoi_seeds(): `density_power` is read by method = \"kmeans\" only; method = ",
      "\"%s\" ignores it."), method))
    expect_identical(r$value, base$value)
    r <- partition_warnings(f(density_k = 12))
    expect_identical(.ss_text(r), sprintf(paste0(
      "get_voronoi_seeds(): `density_k` is read by method = \"kmeans\" only; method = ",
      "\"%s\" ignores it."), method))
    r <- partition_warnings(f(density_power = 0.5, density_k = 12))
    expect_identical(.ss_text(r), sprintf(paste0(
      "get_voronoi_seeds(): `density_power` and `density_k` are read by method = \"kmeans\" ",
      "only; method = \"%s\" ignores them."), method))
    # A kind or a starting rule named is the present warning, and one only.
    for (args in list(list(partition = "density"), list(init = "random"),
                      list(partition = "density", density_power = 0.5))) {
      r <- partition_warnings(do.call(f, args))
      expect_identical(.ss_text(r), sprintf(paste0(
        "get_voronoi_seeds(): `partition` and `init` are read by method = \"kmeans\" only; ",
        "method = \"%s\" ignores them."), method))
    }
    # The k-means kind named: nothing of the kind, and a density setting
    # beside it is the present warning of that kind.
    expect_length(partition_warnings(f(partition = "kmeans"))$warnings, 0L)
    r <- partition_warnings(f(partition = "kmeans", density_power = 0.5))
    expect_identical(.ss_text(r), paste0(
      "get_voronoi_seeds(): `density_power` is read by partition = \"density\" only, and ",
      "partition = \"kmeans\" ignores it."))
  }
  # The line is for a random or provided seeding alone.  A fresh run of the
  # density kind as the function's own reads the settings, and is told
  # nothing: it is the run of the same call with the kind named.
  vsk <- .ss_own(voronoi_seeds_kmeans, "density")
  for (setting in list(list(density_power = 0.5), list(density_k = 12),
                       list(density_power = 0.5, density_k = 12))) {
    for (f in list(function(...) gvs(method = "kmeans", n = 25, sample_points = pts,
                                     set_seed = 1, kmeans_nstart = 3, ...),
                   function(...) vsk(pts, 25, set_seed = 1, nstart = 3, ...))) {
      info <- paste(names(setting), collapse = ", ")
      r <- partition_warnings(do.call(f, setting))
      expect_identical(.ss_text(r), character(0), info = info)
      named <- partition_warnings(do.call(f, c(list(partition = "density"), setting)))
      expect_identical(.ss_text(named), character(0), info = info)
      expect_identical(r$value, named$value, info = info)
      expect_identical(attr(r$value, "partition")[names(setting)],
                       lapply(setting, function(v) if (v == 12) 12L else v), info = info)
      # With the k-means kind named beside it, the present warning that the
      # kind does not read the setting, and no other.
      r <- partition_warnings(do.call(f, c(list(partition = "kmeans"), setting)))
      expect_length(r$warnings, 1L)
      expect_match(.ss_text(r), "read by partition = \"density\" only, and partition = \"kmeans\" ignores",
                   fixed = TRUE, info = info)
    }
  }
  # With the k-means kind as the function's own the added line cannot be
  # reached: the density kind is in force only where the call named it.
  own_kmeans <- .ss_own(get_voronoi_seeds, "kmeans")
  for (args in list(list(density_power = 0.5), list(partition = "density", density_power = 0.5),
                    list(partition = "kmeans", density_k = 12))) {
    r <- partition_warnings(do.call(own_kmeans, c(list(bnd, method = "random", n = 5,
                                                       set_seed = 2), args)))
    expect_false(any(grepl("is read by method = \"kmeans\" only; method", .ss_text(r),
                           fixed = TRUE)))
  }
})

test_that("with the kind flipped and the starts NULL, a count as a number is 3 density starts", {
  pts <- .ss_pts()
  gvs <- .ss_own(get_voronoi_seeds, "density")
  formals(gvs)["kmeans_nstart"] <- list(NULL)
  vsk <- .ss_own(voronoi_seeds_kmeans, "density")
  formals(vsk)["nstart"] <- list(NULL)
  expect_null(formals(gvs)$kmeans_nstart)
  a <- partition_quietly(gvs(method = "kmeans", n = 25, sample_points = pts, set_seed = 1))
  expect_identical(a, partition_quietly(
    get_voronoi_seeds(method = "kmeans", n = 25, sample_points = pts, set_seed = 1,
                      partition = "density", kmeans_nstart = 3)))
  expect_identical(attr(a, "kmeans")$nstart, 3L)
  expect_identical(partition_quietly(vsk(pts, 25, set_seed = 1)),
                   partition_quietly(voronoi_seeds_kmeans(pts, 25, set_seed = 1,
                                                          partition = "density", nstart = 3)))
  # And the k-means kind named there is the run of 10 random starts.
  expect_identical(gvs(method = "kmeans", n = 25, sample_points = pts, set_seed = 1,
                       partition = "kmeans"),
                   get_voronoi_seeds(method = "kmeans", n = 25, sample_points = pts, set_seed = 1,
                                     partition = "kmeans", kmeans_nstart = 10))
  expect_identical(vsk(pts, 25, partition = "kmeans"),
                   voronoi_seeds_kmeans(pts, 25, partition = "kmeans", nstart = 10))
  # With the k-means kind as their own and the starts NULL: the same run.
  gvk <- .ss_own(get_voronoi_seeds, "kmeans")
  formals(gvk)["kmeans_nstart"] <- list(NULL)
  expect_identical(gvk(method = "kmeans", n = 25, sample_points = pts, set_seed = 1),
                   get_voronoi_seeds(method = "kmeans", n = 25, sample_points = pts, set_seed = 1,
                                     partition = "kmeans", kmeans_nstart = 10))
})
