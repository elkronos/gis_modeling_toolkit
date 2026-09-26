# ===========================================================================
# Regressions from the adversarial review of the tessellation slice:
#
#   * coerce_to_points(mode = "auto") segfaulted R on an EMPTY
#     MULTILINESTRING (st_cast() gives one empty part, and st_line_sample()
#     on it crashes in sf 1.0.x), and refused an EMPTY LINESTRING outright.
# ===========================================================================


# ---------------------------------------------------------------------------
# Empty lines in coerce_to_points()
# ---------------------------------------------------------------------------

# Run `code` in a fresh R process that loads this same copy of spatialkit and
# return its exit status, its output and whatever it saved as `res`.  The
# regression being guarded is a SEGFAULT: run in this process, a relapse
# would kill the whole test run instead of failing one test.
.rt_run_child <- function(code) {
  ns_path   <- getNamespaceInfo(asNamespace("spatialkit"), "path")
  installed <- file.exists(file.path(ns_path, "Meta", "package.rds"))
  load_line <- if (installed) {
    sprintf("suppressPackageStartupMessages(library(spatialkit, lib.loc = %s))",
            deparse(dirname(ns_path)))
  } else {
    sprintf("suppressMessages(pkgload::load_all(%s, quiet = TRUE))",
            deparse(ns_path))
  }
  script <- tempfile(fileext = ".R")
  out    <- tempfile(fileext = ".rds")
  on.exit(unlink(c(script, out)), add = TRUE)
  writeLines(c(sprintf(".libPaths(%s)", paste(deparse(.libPaths()), collapse = "")),
               load_line,
               code,
               sprintf("saveRDS(res, %s)", deparse(out))),
             script)
  output <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
                                     c("--vanilla", shQuote(script)),
                                     stdout = TRUE, stderr = TRUE))
  status <- attr(output, "status")
  list(status = if (is.null(status)) 0L else as.integer(status),
       output = output,
       res    = if (file.exists(out)) readRDS(out) else NULL)
}


test_that("coerce_to_points turns empty lines into empty POINTs instead of crashing", {
  child <- .rt_run_child(c(
    "seg <- function(x0, y0, x1, y1) rbind(c(x0, y0), c(x1, y1))",
    "mls <- function(...) sf::st_multilinestring(list(...))",
    "try_ <- function(expr) tryCatch(expr, error = function(e) conditionMessage(e))",
    "res <- list()",
    # An empty MULTILINESTRING between two real ones, projected.
    "x <- sf::st_sf(id = 1:3, geometry = sf::st_sfc(",
    "  mls(seg(0, 0, 10, 0)), sf::st_multilinestring(), mls(seg(0, 0, 0, 20)),",
    "  crs = 32617))",
    "res$mls_proj <- try_(coerce_to_points(x, 'auto'))",
    # The same in lon/lat, projected temporarily (the default).
    "x_ll <- sf::st_sf(id = 1:3, geometry = sf::st_sfc(",
    "  mls(seg(-80, 35, -79, 35)), sf::st_multilinestring(),",
    "  mls(seg(-80, 35, -80, 36)), crs = 4326))",
    "res$mls_ll <- try_(suppressWarnings(coerce_to_points(x_ll, 'auto')))",
    # A layer of one empty row.
    "res$mls_one <- try_(coerce_to_points(x[2, ], 'auto'))",
    # An empty PART beside a real one: the empty part used to be sampled
    # whenever it came first among parts of equal (zero) length.
    "x_part <- sf::st_sf(id = 1:2, geometry = sf::st_sfc(",
    "  sf::st_multilinestring(list(matrix(numeric(0), 0, 2), seg(0, 0, 10, 0))),",
    "  sf::st_multilinestring(list(matrix(numeric(0), 0, 2), seg(3, 4, 3, 4))),",
    "  crs = 32617))",
    "res$mls_part <- try_(coerce_to_points(x_part, 'auto'))",
    # An empty LINESTRING, in both modes that sample lines.
    "x_ls <- sf::st_sf(id = 1:3, geometry = sf::st_sfc(",
    "  sf::st_linestring(seg(0, 0, 10, 0)), sf::st_linestring(),",
    "  sf::st_linestring(seg(0, 0, 0, 20)), crs = 32617))",
    "res$ls_auto <- try_(coerce_to_points(x_ls, 'auto'))",
    "res$ls_mid  <- try_(coerce_to_points(x_ls, 'line_midpoint'))",
    # The documented cleaning path: prep_model_data() drops empty rows, but
    # it pointizes first, so this is where a GeoPackage null geometry crashed.
    "d <- x; d$y <- c(1, 2, 3); d$x1 <- c(0.1, 0.5, 0.9)",
    "res$prep <- try_(prep_model_data(d, 'y', 'x1'))"
  ))

  expect_equal(child$status, 0L,
               info = paste(utils::tail(child$output, 25), collapse = "\n"))
  res <- child$res
  expect_type(res, "list")
  if (!is.list(res)) return(invisible())

  xy <- function(p) unname(sf::st_coordinates(p)[, 1:2, drop = FALSE])
  for (nm in c("mls_proj", "mls_ll", "mls_one", "mls_part", "ls_auto",
               "ls_mid", "prep")) {
    expect_s3_class(res[[nm]], "sf")
  }
  if (!all(vapply(res, inherits, logical(1), "sf"))) return(invisible())

  # Row for row with the input: the empty feature is an empty POINT in its
  # own row (st_coordinates() gives it an NA row), and its neighbours get
  # their true midpoints.
  for (nm in c("mls_proj", "ls_auto", "ls_mid")) {
    got <- res[[nm]]
    expect_equal(got$id, 1:3, info = nm)
    expect_true(all(sf::st_geometry_type(got) == "POINT"), info = nm)
    expect_equal(sf::st_is_empty(got), c(FALSE, TRUE, FALSE), info = nm)
    expect_equal(xy(got)[c(1, 3), ], rbind(c(5, 0), c(0, 10)), info = nm)
  }

  ll <- res$mls_ll
  expect_equal(sf::st_is_empty(ll), c(FALSE, TRUE, FALSE))
  expect_equal(sf::st_crs(ll), sf::st_crs(4326))
  expect_equal(xy(ll)[c(1, 3), ], rbind(c(-79.5, 35), c(-80, 35.5)),
               tolerance = 1e-4)

  expect_equal(nrow(res$mls_one), 1L)
  expect_true(sf::st_is_empty(res$mls_one))

  # The real part is the one sampled; a zero-length part gives its point.
  expect_false(any(sf::st_is_empty(res$mls_part)))
  expect_equal(xy(res$mls_part), rbind(c(5, 0), c(3, 4)))

  # prep_model_data() drops the empty row and says why.
  expect_equal(nrow(res$prep), 2L)
  expect_equal(res$prep$y, c(1, 3))
  expect_equal(attr(res$prep, "dropped")$n_geometry, 1L)
  expect_equal(attr(res$prep, "dropped")$which, 2L)
})
