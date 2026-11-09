# tests/testthat/helper-level-streams.R
# ---------------------------------------------------------------------------
# The stream resolution_profile() runs a level's partition on.  After the
# subsample it draws one integer from the stream of `seed`, and each level L
# runs under set.seed() of that integer plus L, the stream put back after
# it.  On a layer of `n` located rows with no subsample to draw (n at most
# `sample_n`, or a ladder of sizes),
#   set.seed(profile_level_seed(seed, n, L))
# before .partition_centres() gives the starts the profile's level L draws.
# ---------------------------------------------------------------------------

profile_level_seed <- function(seed, n, level) {
  cleanup <- spatialkit:::.with_seed(seed)
  on.exit(cleanup(), add = TRUE)
  sample.int(.Machine$integer.max - as.integer(n), 1L) + as.integer(level)
}
