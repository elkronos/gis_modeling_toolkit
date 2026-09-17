# =============================================================================
# 00  Run the whole tour
# =============================================================================
#   source(system.file("scripts", "00-run-all.R", package = "spatialkit"))
#
# Runs scripts 01 to 10 in order, reporting how long each took and what it
# skipped. Scripts 09 (GWR) and 10 (Bayesian) need optional packages and skip
# themselves cleanly when those are missing.
#
# Scripts 01 to 09 run in about a minute between them. Script 10 takes several
# more on its own, because Stan compiles the model separately for each of its
# four fits. Skip it by sourcing 01 to 09 individually if you are in a hurry.
#
# To write every figure to files instead of drawing them:
#   Sys.setenv(SPATIALKIT_TOUR_OUTPUT = "~/spatialkit-tour")
# To run one script on its own, source it directly -- each is self-contained.
# =============================================================================
.tour_dir <- if (nzchar(system.file("scripts", package = "spatialkit")))
  system.file("scripts", package = "spatialkit") else "."

# Non-interactive by default: with ten scripts the per-figure pause is a chore.
# Unset SPATIALKIT_TOUR_PAUSE before sourcing if you want it back.
if (!nzchar(Sys.getenv("SPATIALKIT_TOUR_PAUSE")))
  Sys.setenv(SPATIALKIT_TOUR_PAUSE = "no")

scripts <- sort(grep("^[0-9]{2}-", list.files(.tour_dir, pattern = "[.]R$"),
                     value = TRUE))
scripts <- setdiff(scripts, "00-run-all.R")

results <- data.frame(script = scripts, seconds = NA_real_,
                      status = NA_character_, stringsAsFactors = FALSE)

for (i in seq_along(scripts)) {
  cat("\n\n", strrep("#", 78), "\n# ", scripts[i], "\n", strrep("#", 78), "\n",
      sep = "")
  t0 <- Sys.time()
  ok <- tryCatch({
    # local() so one script's objects cannot leak into the next.
    local(source(file.path(.tour_dir, scripts[i]), local = TRUE, echo = FALSE))
    "ok"
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  results$seconds[i] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  results$status[i]  <- ok
}

cat("\n\n", strrep("=", 78), "\n", sep = "")
results$seconds <- round(results$seconds, 1)
print(results, row.names = FALSE)
failed <- sum(results$status != "ok")
cat(sprintf("\n%d script(s) ran, %d failed, %.0f seconds in total.\n",
            nrow(results), failed, sum(results$seconds)))
out <- Sys.getenv("SPATIALKIT_TOUR_OUTPUT", unset = "")
if (nzchar(out)) cat("Figures are in", out, "\n")
