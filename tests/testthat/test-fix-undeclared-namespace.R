# R CMD check warned "'loadNamespace' or 'requireNamespace' call not
# declared from: 'lwgeom'": .sf_box_from_coordinates() asked
# requireNamespace("lwgeom") with the name written out, and lwgeom is not in
# DESCRIPTION. It asks the same question with the name passed as a value.
# These tests read the package's functions, so they run on the installed
# package as R CMD check runs them.

# Every package a function of the namespace names in a literal
# requireNamespace(), loadNamespace(), library() or require() call, or
# before :: or :::.
.un_named_packages <- function() {
  ns <- asNamespace("spatialkit")
  found <- character(0)
  walk <- function(e) {
    if (is.call(e)) {
      head <- e[[1L]]
      if (is.name(head)) {
        h <- as.character(head)
        if (h %in% c("requireNamespace", "loadNamespace", "library", "require") &&
            length(e) >= 2L && is.character(e[[2L]]))
          found <<- c(found, e[[2L]])
        if (h %in% c("::", ":::") && length(e) >= 2L)
          found <<- c(found, as.character(e[[2L]]))
      }
      for (a in as.list(e)) if (!missing(a)) walk(a)
    } else if (is.function(e)) {
      walk(body(e))
      for (f in formals(e)) if (!missing(f)) walk(f)
    }
  }
  for (nm in ls(ns, all.names = TRUE)) {
    f <- get(nm, envir = ns)
    if (is.function(f)) walk(f)
  }
  unique(found)
}

.un_declared_packages <- function() {
  # The DESCRIPTION of the package under test, installed or loaded from
  # its source.
  d <- read.dcf(system.file("DESCRIPTION", package = "spatialkit"))
  fields <- d[1L, intersect(c("Depends", "Imports", "Suggests", "LinkingTo"), colnames(d))]
  pk <- trimws(sub("\\(.*$", "", unlist(strsplit(paste(fields, collapse = ","), ","))))
  base <- rownames(utils::installed.packages(priority = "base"))
  unique(c(pk[nzchar(pk)], base, "spatialkit"))
}

test_that("every package the code names by literal is declared in DESCRIPTION", {
  named <- .un_named_packages()
  # The walk sees the calls it is meant to see: sf is named before ::.
  expect_true("sf" %in% named)
  expect_identical(setdiff(named, .un_declared_packages()), character(0))
})

test_that("the lwgeom test asks what sf asks, without naming lwgeom by literal", {
  expect_false("lwgeom" %in% .un_named_packages())
  # The name is passed as a value here too: R CMD check reads a literal
  # in the tests as an undeclared use as well.
  pkg <- "lwgeom"
  expect_identical(spatialkit:::.sf_box_from_coordinates(),
                   !requireNamespace(pkg, quietly = TRUE))
})
