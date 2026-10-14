# tests/testthat/test-core-defaults.R
# ---------------------------------------------------------------------------
# The defaults that one later step switches together, at their present
# values: eleven formal defaults in five functions, and the value of one
# internal function.  Until that step a call that names none of the new
# arguments returns what it returned before they existed, and these are the
# values that make it so.  The step that switches them edits the expected
# values here in the same change.
# ---------------------------------------------------------------------------

test_that("the eleven defaults that switch together hold their present values", {
  profile <- formals(resolution_profile)
  seeds   <- formals(get_voronoi_seeds)
  kmeans  <- formals(voronoi_seeds_kmeans)
  select  <- formals(select_resolution)
  values  <- formals(summarize_by_cell)

  # 1 to 3: no estimate is scored, the ladder is the automatic one, and the
  # variogram is estimated from the response.
  expect_true("estimate" %in% names(profile))
  expect_null(profile$estimate)
  expect_true("obs_per_cell" %in% names(profile))
  expect_null(profile$obs_per_cell)
  expect_true("sac" %in% names(profile))
  expect_null(profile$sac)

  # 4: a profile is read by cp.
  expect_identical(eval(select$criterion),
                   c("cp", "reliability", "elbow", "moran_z", "loo", "size", "auto"))
  expect_identical(eval(select$criterion)[1L], "cp")

  # 5: the summary reports plain means.
  expect_identical(eval(values$estimate), c("mean", "pooled"))
  expect_identical(eval(values$estimate)[1L], "mean")

  # 6 to 8: the regions are k-means regions in all three functions.
  for (f in list(profile, seeds, kmeans)) {
    expect_identical(eval(f$partition), c("kmeans", "density"))
    expect_identical(eval(f$partition)[1L], "kmeans")
  }

  # 9 to 11: the number of starts is a number, whatever the kind.
  expect_identical(profile$nstart, 25L)
  expect_identical(seeds$kmeans_nstart, 10)
  expect_identical(kmeans$nstart, 10)
})

test_that("Moran's z is left out without its packages on a ladder of sizes only", {
  expect_identical(.moran_skip_any_ladder(), FALSE)
  expect_identical(formals(.moran_skip_any_ladder), NULL)
})

test_that("the values the defaults switch to are accepted where this version reads them", {
  # Each value one of the eleven switches to is a value its argument takes
  # today, so the switch edits defaults and nothing else.
  expect_true("auto" %in% eval(formals(select_resolution)$criterion))
  expect_true("pooled" %in% eval(formals(summarize_by_cell)$estimate))
  for (f in list(resolution_profile, get_voronoi_seeds, voronoi_seeds_kmeans))
    expect_true("density" %in% eval(formals(f)$partition))
  expect_identical(.NSTART_DENSITY, 3L)
  expect_identical(.OBS_PER_CELL_DEFAULT, 4.5)
  # The arguments this stage's rules read sit after the present ones.
  expect_identical(names(formals(select_resolution)),
                   c("profile", "criterion", "tol", "range_floor", "diagnostic", "obs_per_cell"))
  expect_null(formals(select_resolution)$obs_per_cell)
  expect_identical(formals(select_resolution)$tol, 0.02)
  expect_identical(formals(select_resolution)$range_floor, TRUE)
  expect_identical(formals(select_resolution)$diagnostic, FALSE)
})
