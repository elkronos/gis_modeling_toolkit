# dev/baseline-core-frozen.R
# ---------------------------------------------------------------------------
# The frozen baseline of the cell pipeline.
#
# Runs the calls of tests/testthat/helper-core-frozen.R on the baseline code
# and saves what each returns, prints, raises and logs, and where it leaves
# the random-number stream, as tests/testthat/fixtures/core-frozen.rds.
# tests/testthat/test-core-frozen.R runs the same calls on the code under
# test and compares, so a change to what a default call returns fails that
# test.
#
# THE BASELINE is the package code with the two partition kinds
# (partition = "kmeans" | "density", R/partition.R) and nothing built on top
# of them.  It is named by a git TREE, not by a commit: when the baseline was
# taken that code was staged and not yet committed, and BASELINE_TREE is the
# tree of the index then (`git write-tree`), on top of commit c7a2a10.  The
# commit that records this tree (`git rev-parse <commit>^{tree}` prints
# BASELINE_TREE) is the baseline commit.
#
# The script refuses to run on any other code: it hashes R/, DESCRIPTION and
# NAMESPACE of the package directory as git does and compares them with the
# entries of BASELINE_TREE (`git ls-tree <BASELINE_TREE> R DESCRIPTION
# NAMESPACE` prints the three hashes below).
#
# Usage, from the root of the working tree:
#   Rscript --vanilla dev/baseline-core-frozen.R [<baseline dir> [<fixture file>]]
# <baseline dir> is a directory holding the baseline code, "." by default.
# Once the working tree has moved on, export the baseline first, by its tree
# or by the commit that records it:
#   mkdir <baseline dir>
#   git archive c64fcac7be7ae3c83027198e1176d25368727633 | tar -x -C <baseline dir>
# The calls and the recorder are always read from the working tree
# (tests/testthat/helper-core-frozen.R): they are what the test replays.
# When a fixture is already at <fixture file>, the script says how many of
# its records the new run gives again, so that a call added to the list, or
# a fixture taken again on another machine, can be seen to move no record.
#
# The test holds a run to identical() only on the platform the fixture was
# recorded on.  When R, sf, GEOS, PROJ or one of the packages the fixture
# names changes on the machine the work is done on, take the fixture again
# from an export of the baseline; do not take it from code that has moved.
#
# The fixture stores the text of every call, and the test compares it with
# the helper's list.  When the calls are rewritten to name the defaults they
# were recorded under (the step that changes those defaults), the records
# stay and only the stored text is replaced:
#   Rscript --vanilla dev/baseline-core-frozen.R --calls-only [<fixture file>]
# That mode runs nothing and needs no baseline code.  It refuses a list with
# other names than the fixture's records, prints every call whose text
# changed, for the commit to list, and keeps the text the records were taken
# with as `calls_recorded`.  After it the records can no longer be taken
# again from this list: the baseline code does not know the arguments the
# rewritten calls name.
#
# `dev/` is in .Rbuildignore, so nothing here ships.
# ---------------------------------------------------------------------------

BASELINE_TREE <- "c64fcac7be7ae3c83027198e1176d25368727633"
BASELINE_CODE <- c(R           = "2df204570e3e0378e1919c6bd1f9a1c5be24d8ca",   # tree
                   DESCRIPTION = "645da4e0c2ef52c044b40e7885104ab1a17ce45b",   # blob
                   NAMESPACE   = "43ea3214d0eebda8f3a1b2792075b4d61d909ef1")   # blob

args     <- commandArgs(trailingOnly = TRUE)
calls_only <- length(args) >= 1L && identical(args[[1L]], "--calls-only")
if (calls_only) args <- c(".", args[-1L])
base_dir <- if (length(args) >= 1L) args[[1L]] else "."
fixture  <- if (length(args) >= 2L) args[[2L]] else
  file.path("tests", "testthat", "fixtures", "core-frozen.rds")
helper   <- file.path("tests", "testthat", "helper-core-frozen.R")
if (!file.exists(helper))
  stop("Run from the root of the working tree: ", helper, " not found.")
if (!file.exists(file.path(base_dir, "DESCRIPTION")))
  stop("No package in '", base_dir, "'.")

# ---- The stored text of the calls alone --------------------------------------------

if (calls_only) {
  source(helper, local = FALSE)
  if (!file.exists(fixture)) stop("No fixture at '", fixture, "'.")
  baseline <- readRDS(fixture)
  now <- core_frozen_call_text()
  if (!identical(names(now), names(baseline$records)))
    stop("The helper's calls are not the fixture's records by name: only in the helper, ",
         paste(setdiff(names(now), names(baseline$records)), collapse = ", "),
         "; only in the fixture, ",
         paste(setdiff(names(baseline$records), names(now)), collapse = ", "),
         "; or their order differs. A call is added or removed on the baseline code only.")
  changed <- names(now)[!mapply(identical, now, baseline$calls[names(now)])]
  for (nm in changed)
    cat(sprintf("%s\n  was: %s\n  now: %s\n", nm, baseline$calls[[nm]], now[[nm]]))
  if (is.null(baseline$calls_recorded)) baseline$calls_recorded <- baseline$calls
  baseline$calls <- now
  saveRDS(baseline, fixture, compress = "xz")
  cat(sprintf("\n%d of %d calls have another text; %s holds the new text and the records it had.\n",
              length(changed), length(now), fixture))
  quit(save = "no", status = 0L)
}

# ---- The code is the baseline's ------------------------------------------------

# The hash git gives a file, and the one it gives a directory of plain files.
git_blob <- function(path) {
  bytes <- readBin(path, "raw", n = file.size(path))
  digest::digest(c(charToRaw(sprintf("blob %d", length(bytes))), as.raw(0L), bytes),
                 algo = "sha1", serialize = FALSE)
}
git_tree <- function(dir) {
  files <- sort(list.files(dir, all.files = TRUE, no.. = TRUE), method = "radix")
  if (any(dir.exists(file.path(dir, files))))
    stop("git_tree(): '", dir, "' holds a directory; only plain files are hashed.")
  entry <- function(f) {
    sha <- git_blob(file.path(dir, f))
    c(charToRaw(paste("100644", f)), as.raw(0L),
      as.raw(strtoi(substring(sha, seq(1L, 39L, 2L), seq(2L, 40L, 2L)), 16L)))
  }
  body <- unlist(lapply(files, entry))
  digest::digest(c(charToRaw(sprintf("tree %d", length(body))), as.raw(0L), body),
                 algo = "sha1", serialize = FALSE)
}
found <- c(R           = git_tree(file.path(base_dir, "R")),
           DESCRIPTION = git_blob(file.path(base_dir, "DESCRIPTION")),
           NAMESPACE   = git_blob(file.path(base_dir, "NAMESPACE")))
if (!identical(found, BASELINE_CODE)) {
  off <- names(found)[found != BASELINE_CODE]
  stop("The code in '", base_dir, "' is not the baseline's: ",
       paste(sprintf("%s is %s, the baseline's is %s", off, found[off], BASELINE_CODE[off]),
             collapse = "; "),
       ". Export the baseline (see the header) and pass its directory.")
}

# ---- Load it, and the calls -----------------------------------------------------

suppressPackageStartupMessages(library(sf))
pkgload::load_all(base_dir, quiet = TRUE)
loaded_from <- normalizePath(getNamespaceInfo(asNamespace("spatialkit"), "path"))
if (!identical(loaded_from, normalizePath(base_dir)))
  stop("spatialkit is loaded from ", loaded_from, ", not from ", normalizePath(base_dir), ".")
message("spatialkit loaded from : ", loaded_from)
source(helper, local = FALSE)

# ---- The data --------------------------------------------------------------------

# The layers, as data frames of coordinates and columns.  The first three are
# in EPSG:32632.
#   pts   400 points spread over a 1 km square.  `z` is a Gaussian field with
#         an exponential covariance (range parameter 90 m) and a nugget;
#         `z_na` is `z` with 50 values missing; `z_bad` is `z` with missing
#         and infinite values; `flag` is logical; `none` has no value;
#         `x1` and `x2` are predictors with 20 values missing each; `noise`
#         has no spatial structure.
#   clus  500 points, 470 of them in four tight clusters, with a response `w`.
#   big   1,800 points over a 3 km square, more than resolution_profile()'s
#         `sample_n`, with a smooth response `z` plus noise.
# Three more were added after the first list of calls.  They are drawn after
# the three above, so those are the numbers they were.
#   ll       `pts` with its coordinates in longitude and latitude (EPSG:4326).
#   twice    the first 150 locations of `pts`, each held by two rows, with a
#            response `z` that differs between the two.
#   lattice  a 12 x 12 lattice at 80 m, with a smooth response `z` plus
#            noise, a constant column `const` and an integer column `count`.
core_frozen_data <- function() {
  set.seed(20261012, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  n <- 400L
  x <- 5e5 + runif(n, 0, 1000)
  y <- 5e6 + runif(n, 0, 1000)
  D <- as.matrix(dist(cbind(x, y)))
  z <- 10 + as.numeric(t(chol(exp(-D / 90) + diag(0.3, n))) %*% rnorm(n))
  z_na <- z
  z_na[sample.int(n, 50L)] <- NA
  z_bad <- z
  z_bad[sample.int(n, 12L)] <- c(NA, NaN, Inf, -Inf)
  x1 <- rnorm(n)
  x2 <- 0.5 * z + rnorm(n)
  x1[sample.int(n, 20L)] <- NA
  x2[sample.int(n, 20L)] <- NA
  pts <- data.frame(x = x, y = y, z = z, z_na = z_na, z_bad = z_bad,
                    flag = z > stats::median(z), none = NA_real_, x1 = x1, x2 = x2,
                    noise = rnorm(n))

  n_bg <- 30L; n_cl <- 470L; k <- 4L
  bg  <- cbind(runif(n_bg, 0, 1000), runif(n_bg, 0, 1000))
  ctr <- cbind(runif(k, 150, 850), runif(k, 150, 850))
  who <- sample.int(k, n_cl, replace = TRUE)
  cl  <- ctr[who, , drop = FALSE] + matrix(rnorm(2L * n_cl, sd = 20), n_cl, 2L)
  xy  <- rbind(bg, cl)
  xy  <- xy[sample.int(nrow(xy)), , drop = FALSE]
  clus <- data.frame(x = 5e5 + xy[, 1L], y = 5e6 + xy[, 2L])
  clus$w <- 5 + xy[, 1L] / 500 + sin(xy[, 2L] / 150) + rnorm(nrow(xy), sd = 0.4)

  m  <- 1800L
  bx <- runif(m, 0, 3000)
  by <- runif(m, 0, 3000)
  big <- data.frame(x = 5e5 + bx, y = 5e6 + by,
                    z = sin(bx / 400) + cos(by / 300) + rnorm(m, sd = 0.5))

  lonlat <- sf::st_coordinates(sf::st_transform(
    sf::st_as_sf(pts[, c("x", "y")], coords = c("x", "y"), crs = 32632), 4326))
  ll <- pts
  ll$x <- unname(lonlat[, 1L])
  ll$y <- unname(lonlat[, 2L])

  both  <- c(1:150, 1:150)
  twice <- data.frame(x = pts$x[both], y = pts$y[both],
                      z = pts$z[both] + rnorm(length(both), sd = 0.2))

  g <- expand.grid(gx = seq(40, 920, by = 80), gy = seq(40, 920, by = 80))
  lattice <- data.frame(x = 5e5 + g$gx, y = 5e6 + g$gy,
                        z = sin(g$gx / 200) + cos(g$gy / 150) + rnorm(nrow(g), sd = 0.3),
                        const = 7)
  lattice$count <- as.integer(round(10 * lattice$z))
  list(pts = pts, clus = clus, big = big, ll = ll, twice = twice, lattice = lattice)
}

# ---- Record ------------------------------------------------------------------------

data  <- core_frozen_data()
calls <- core_frozen_calls()
t0 <- proc.time()[[3L]]
records <- core_frozen_run(data, calls)
message(sprintf("%d calls recorded in %.1f s", length(records), proc.time()[[3L]] - t0))
# A second run in the same session has to give the same records: a call
# whose result moved with the session's state could not be replayed.
again <- core_frozen_run(data, calls)
moved <- names(records)[!mapply(identical, records, again)]
if (length(moved))
  stop("These calls do not return the same record twice: ", paste(moved, collapse = ", "))
# The same calls where FNN and Matrix are not installed.  The test replays
# the list that way too.  A call whose record is another without them has
# that record saved beside the first, and the replay is held to it.
bare <- core_frozen_run(data, calls, without = c("FNN", "Matrix"))
bare_again <- core_frozen_run(data, calls, without = c("FNN", "Matrix"))
moved <- names(bare)[!mapply(identical, bare, bare_again)]
if (length(moved))
  stop("These calls do not return the same record twice without FNN and Matrix: ",
       paste(moved, collapse = ", "))
other <- names(records)[!mapply(identical, records, bare)]
message(sprintf("without FNN and Matrix %d of %d records are another%s", length(other),
                length(records),
                if (length(other)) paste0(": ", paste(other, collapse = ", ")) else ""))

# What a fixture already there holds of these calls.
if (file.exists(fixture)) {
  had <- readRDS(fixture)$records
  old_parts <- c("value", "output", "printed", "log")
  same <- function(nm) {
    a <- had[[nm]]; b <- records[[nm]]
    said <- function(r) lapply(r$conditions, function(c) c[c("when", "kind", "class", "message")])
    identical(a[old_parts], b[old_parts]) && identical(said(a), said(b))
  }
  both <- intersect(names(had), names(records))
  off  <- both[!vapply(both, same, logical(1))]
  message(sprintf("%s held %d records: %d are also in this list, %d of them with another value, condition, output, print or log%s",
                  fixture, length(had), length(both), length(off),
                  if (length(off)) paste0(" (", paste(off, collapse = ", "), ")") else ""))
}

baseline <- list(
  info = c(list(tree = BASELINE_TREE,
                tree_is = paste("a git tree (`git write-tree` of the index), not a commit:",
                                "the code was staged on top of commit c7a2a10"),
                code = BASELINE_CODE,
                recorded = format(Sys.Date())),
           core_frozen_platform()),
  data = data,
  calls = core_frozen_call_text(calls),
  records = records,
  # What the calls return where FNN and Matrix are not installed: the records
  # above, but for these.
  without = list(packages = c("FNN", "Matrix"), records = bare[other]))
dir.create(dirname(fixture), recursive = TRUE, showWarnings = FALSE)
saveRDS(baseline, fixture, compress = "xz")

n_cond <- function(kind) sum(vapply(records, function(r)
  sum(vapply(r$conditions, function(c) identical(c$kind, kind), logical(1))), numeric(1)))
cat(sprintf("\nWrote %s (%.0f kB)\n", fixture, file.size(fixture) / 1024))
cat(sprintf("  %d calls; %d with a value, %d that stop\n", length(records),
            sum(!vapply(records, function(r) is.null(r$value), logical(1))),
            n_cond("error")))
cat(sprintf("  %d warnings, %d messages, %d logged lines, %d prints\n",
            n_cond("warning"), n_cond("message"),
            sum(lengths(lapply(records, `[[`, "log"))),
            sum(!vapply(records, function(r) is.null(r$printed), logical(1)))))
cat(sprintf("  %d calls leave the random-number stream moved\n",
            sum(!vapply(records, function(r) r$rng$untouched, logical(1)))))
cat(sprintf("  %d records are another without FNN and Matrix\n", length(other)))
p <- baseline$info
cat(sprintf("  %s, %s; sf %s, GEOS %s, PROJ %s\n", p$platform, p$r, p$sf, p$geos, p$proj))
