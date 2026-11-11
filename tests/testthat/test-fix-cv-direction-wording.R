# tests/testthat/test-fix-cv-direction-wording.R
# ---------------------------------------------------------------------------
# Which way a cross-validated score errs, and what folds = NULL estimates.
#
# A cross-validated error is the error at the distance the folds put between
# the held-out points and their training data.  It reads low for predictions
# made farther from the data than that and high for predictions made closer
# (?fold_separation), so no page may say which way a score errs without
# naming where the model will predict.  The help pages, the vignettes, the
# README and the tour scripts are read as text, sentence by sentence:
#
#   * the statements that gave a direction without a target are gone from
#     every page, each page that held one says instead what the score
#     estimates, and make_folds()'s five warnings on blocks smaller than the
#     range no longer say "may leak";
#   * every sentence that holds a direction word (flatter, optimistic,
#     pessimistic, leak, honest, overstate, understate, conservative, the
#     safe side, reads low or high) names a prediction target in the same
#     sentence, or the word stands in one of the phrases listed below, which
#     use it of something else: an in-sample fit, a selection, the response
#     among its own predictors, intervals and standard errors, a test, the
#     halves of a split;
#   * no sentence pairs a direction with the wrong distance (held-out points
#     that sit farther, read low), and the sentences that give a direction
#     are pinned word for word where they stand;
#   * each cv_*() help says what folds = NULL estimates, and
#     ?cv_block_size_sweep what the rise of its curve measures, and no page
#     presents the sweep as the way to choose a block size;
#   * the texts the code prints or logs are held to the same rules: the
#     print of a forest, the cautions of fit_rf_model(include_coords = TRUE)
#     and of select_features_forward(method = "random_kfold"), and the
#     warning of the fallback to random folds, each line of a print's
#     header as a sentence of its own.
#
# From the source tree the pages are read from man/, vignettes/, README.Rmd
# and inst/scripts/.  Under R CMD check the help comes from the installed Rd
# database, the vignettes from their installed sources and the scripts from
# their installed folder; the README is not installed and is not read there,
# and the vignettes are read only when their sources were installed.
# ---------------------------------------------------------------------------

# ---- The readers --------------------------------------------------------------------

.cvw_source <- function() {
  root <- testthat::test_path("..", "..")
  if (file.exists(file.path(root, "DESCRIPTION")) && dir.exists(file.path(root, "man")))
    normalizePath(root)
  else NULL
}

.cvw_flat <- function(x)
  paste(unlist(lapply(x, function(el) if (is.character(el)) el else .cvw_flat(el))),
        collapse = "")

# Text on one line: no backticks, runs of white space as one space.
.cvw_squash <- function(s) trimws(gsub("[[:space:]]+", " ", gsub("`", "", s, fixed = TRUE)))

# A text cut at its blank lines, each piece on one line.
.cvw_paragraphs <- function(text) {
  p <- strsplit(text, "\n[ \t]*\n", perl = TRUE)[[1L]]
  p <- vapply(p, .cvw_squash, character(1), USE.NAMES = FALSE)
  p[nzchar(p)]
}

# A paragraph cut into sentences: after a full stop, question or exclamation
# mark that is followed by a space and anything but a lower-case letter.
.cvw_sentences <- function(paragraph)
  unlist(strsplit(paragraph, "(?<=[.!?])\\s+(?=[^a-z])", perl = TRUE), use.names = FALSE)

.cvw_cache <- new.env(parent = emptyenv())

# Every help page as parsed Rd, named by its file.
.cvw_rd_all <- function() {
  if (is.null(.cvw_cache$rd)) {
    src <- .cvw_source()
    .cvw_cache$rd <- if (is.null(src)) tools::Rd_db("spatialkit") else {
      files <- list.files(file.path(src, "man"), "\\.Rd$", full.names = TRUE)
      stats::setNames(lapply(files, tools::parse_Rd), basename(files))
    }
  }
  .cvw_cache$rd
}
.cvw_rd <- function(page) .cvw_rd_all()[[paste0(page, ".Rd")]]

.cvw_tag <- function(node) attr(node, "Rd_tag") %||% ""

# The \arguments entries of a page: their text, named by label.
.cvw_rd_arguments <- function(rd) {
  out <- character(0)
  for (node in rd[vapply(rd, .cvw_tag, character(1)) == "\\arguments"])
    for (it in node)
      if (identical(.cvw_tag(it), "\\item") && length(it) == 2L)
        out[[.cvw_squash(.cvw_flat(it[[1L]]))]] <- .cvw_squash(.cvw_flat(it[[2L]]))
  out
}

# The text of the top-level parts of a page with one tag ("\\description").
.cvw_rd_part <- function(rd, tag)
  .cvw_squash(paste(vapply(rd[vapply(rd, .cvw_tag, character(1)) == tag], .cvw_flat,
                           character(1)), collapse = " "))

# The lines of a vignette's source, or of the README in the source tree.
.cvw_vignette_dir <- function() {
  src <- .cvw_source()
  if (!is.null(src)) file.path(src, "vignettes") else system.file("doc", package = "spatialkit")
}
# Whether the vignettes' sources are there to read: under a check that does
# not build the vignettes they are not installed.
.cvw_have_vignettes <- function()
  length(list.files(.cvw_vignette_dir(), "\\.Rmd$")) > 0L
.cvw_vignette <- function(name) {
  file <- file.path(.cvw_vignette_dir(), paste0(name, ".Rmd"))
  if (!file.exists(file)) return(NULL)
  paste(readLines(file, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

# The tour scripts.  A comment loses its marker and the pieces of one printed
# sentence (cat("...\n", "...")) are joined, so that a sentence reads as one.
.cvw_scripts_dir <- function() {
  src <- .cvw_source()
  if (!is.null(src)) file.path(src, "inst", "scripts") else system.file("scripts", package = "spatialkit")
}
.cvw_script <- function(file, dir = .cvw_scripts_dir()) {
  x <- readLines(file.path(dir, file), warn = FALSE, encoding = "UTF-8")
  x <- paste(sub("^[[:space:]]*#+[[:space:]]?", "", x), collapse = "\n")
  gsub("(\\\\n)?\"[[:space:]]*,[[:space:]]*\"", " ", x, perl = TRUE)
}

# Every page as raw text, named "?page", "vignette name", "README" or
# "script file".
.cvw_pages <- function() {
  if (is.null(.cvw_cache$pages)) {
    rd  <- .cvw_rd_all()
    out <- stats::setNames(lapply(rd, .cvw_flat), paste0("?", sub("\\.Rd$", "", names(rd))))
    for (v in sub("\\.Rmd$", "", list.files(.cvw_vignette_dir(), "\\.Rmd$")))
      out[[paste("vignette", v)]] <- .cvw_vignette(v)
    src <- .cvw_source()
    if (!is.null(src) && file.exists(file.path(src, "README.Rmd")))
      out[["README"]] <- paste(readLines(file.path(src, "README.Rmd"), warn = FALSE,
                                         encoding = "UTF-8"), collapse = "\n")
    for (f in list.files(.cvw_scripts_dir(), "\\.R$"))
      out[[paste("script", f)]] <- .cvw_script(f)
    .cvw_cache$pages <- out
  }
  .cvw_cache$pages
}

# The lines of one chunk of a vignette, by its label.
.cvw_chunk <- function(text, label) {
  lines <- strsplit(text, "\n", fixed = TRUE)[[1L]]
  start <- grep(paste0("^```\\{r ", label, "[,}]"), lines)
  if (length(start) != 1L) return(NULL)
  end <- start + which(grepl("^```[[:space:]]*$", lines[-seq_len(start)]))[1L]
  lines[(start + 1L):(end - 1L)]
}


# ---- The rules ----------------------------------------------------------------------

# Words by which a text says which way a score errs.
.CVW_STEMS <- paste0(
  "(?i)\\b(flatter\\w*|optimisti\\w*|pessimisti\\w*|leak\\w*|honest\\w*|",
  "overstat\\w*|understat\\w*|(anti)?conservativ\\w*|safe side|reads? (low|high))\\b")

# Words by which a text names where the model will predict.
.CVW_TARGET <- paste0(
  "(?i)(prediction_points|prediction points|prediction target|",
  "\\b(the|a|this|that) target\\b|will predict|\\bpredicts? (at|over|into|onto)\\b|",
  "\\bpredictions? made\\b|\\bmaps? (of|over|whose)\\b|\\b(the|a|this|that) map\\b)")

# The sentences that hold a direction word and name no target, once the
# phrases that use such a word of something else are taken out.
.cvw_untargeted <- function(paragraphs, other = character(0)) {
  rest <- paragraphs
  for (ph in other) rest <- gsub(ph, " ", rest, fixed = TRUE)
  s <- unlist(lapply(rest, .cvw_sentences), use.names = FALSE)
  s[grepl(.CVW_STEMS, s, perl = TRUE) & !grepl(.CVW_TARGET, s, perl = TRUE)]
}

# A direction paired with the wrong distance: held-out points that sit
# farther than the map's locations make the score read high, not low, and a
# score reads low, not high, for predictions made farther from the data.
.CVW_WRONG_WAY <- c(
  "(?i)\\bread(s)? high for (predictions|locations|a map|maps|one) (made |that predicts? )?farther",
  "(?i)\\bread(s)? low for (predictions|locations|a map|maps|one) (made |that predicts? )?closer",
  "(?i)\\b(sits?|stays?|lies?) farther\\b[^.;]{0,160}?\\bread(s)? low\\b",
  "(?i)\\b(sits?|stays?|lies?) closer\\b[^.;]{0,160}?\\bread(s)? high\\b",
  "(?i)\\bunderstates? the error of a map where the samples cover",
  # A reading named in quotes and given its direction: "farther" reads high,
  # "closer" low.
  "(?i)\\\\?\"farther\\\\?\"[^.;\"]{0,60}?\\bread(s)? low\\b",
  "(?i)\\\\?\"closer\\\\?\"[^.;\"]{0,60}?\\bread(s)? high\\b")
.cvw_wrong_way <- function(text) {
  text <- .cvw_squash(text)
  unlist(lapply(.CVW_WRONG_WAY, function(re)
    regmatches(text, gregexpr(re, text, perl = TRUE))[[1L]]), use.names = FALSE)
}

# The phrases that use a direction word of something other than which way a
# cross-validated score errs for a map, by page.  Each stands on its page.
.CVW_OTHER_SENSES <- list(
  # An in-sample fit, and a score chosen as the best of many.
  "?fitted.gwr_fit" = "summary() on a gwr_fit reports an optimistic fit",
  "?gwr_model_selection" = c(
    "an in-sample criterion is optimistic in a way that a spatially blocked estimate is not",
    # AICc counts the rows at one location as independent.
    "so the gaps overstate what a predictor adds"),
  "?plot.feature_selection" = "optimistically biased by the selection",
  "?select_features_forward" = c(
    "optimistically biased because it was chosen as the best of many",
    "optimistically biased by construction",
    "An honest performance estimate needs data the selection never saw",
    "automatically nested and leak-free",
    # A standard error over the folds.
    "a standard error over rows would understate the noise"),
  # A selection made inside each fold, and folds that failed.
  "?cv_spatial" = c(
    "already nested and leak-free",
    "the fold that fails is often the hardest to predict, so it may flatter the model"),
  # The response among its own predictors.
  "?prep_model_data" = "Using it as its own predictor is leakage no backend catches",
  "README" = c("it is leakage in the forest",
               "the criteria that disagree, the leak, and the default size"),
  # Intervals, standard errors, variances, tests, a comparison of estimands,
  # a print, seeds, distances in degrees, the area of applicability.
  "?cv_bayes" = "whether those intervals are honest at locations the model has not seen",
  "?fit_bayesian_spatial_model" = c(
    "the honest diagnostic is to report both side by side",
    "On a fit with nothing left over the test is conservative"),
  "?summary.resolution_profile" = "a row subset of the result reads honestly",
  "?voronoi_seeds_random" = "It is also the honest choice for a null or sensitivity comparison",
  "?area_of_applicability" = c(
    "giving a small threshold and a conservative AOA",
    "which overstates dissimilarity along directions the model barely uses",
    "silently understate the distances"),
  "?kriging_adequacy" = c(
    "under spatial correlation it understates the error",
    "in that unit overstate the disagreement",
    "how badly the kriging variance is understated (or overstated)",
    "With the nugget understated tenfold"),
  "?plot.spatial_fit" = "residuals of a fitted trend understate the variogram",
  "?residual_morans_i" = c(
    "The random forest is anticonservative under both",
    "the variance is understated whichever is used",
    "GWR is conservative under both"),
  "?summarize_by_cell" = c(
    "the naive SE is anticonservative (too small)",
    "understates sigma^2",
    "it is anticonservative for the block average too",
    "their standard errors are understated"),
  "vignette diagnostics" = c("says the naive standard error was optimistic",
                             "above 1 means it is understated"),
  "vignette spatialkit_nc_demo" =
    "it is anticonservative when points inside a cell are spatially correlated",
  # The selection and estimation halves of select_on = "split".
  "vignette resolution" = c(
    "on an autocorrelated field that leaks the structure you are choosing against",
    "Blocking reduces that leak without removing it"),
  # The tour scripts: R objects, a scatter, kriging variances, a selection,
  # intervals, standard errors.
  "script 00-run-all.R" = "so one script's objects cannot leak into the next",
  "script 02-resolution.R" = c(
    "it is the more honest measure of how much the data pins this down",
    "The leak: choosing the resolution on data you will test on"),
  "script 05-fit-diagnose.R" = c(
    "A flatter-than-1:1 cloud",
    "are the kriging uncertainties honest?",
    "the model understates its own error",
    "the model overstates its own error",
    "The interesting case is honest uncertainties AND a large gain"),
  "script 06-cv-compare.R" = "The honest write-up either (a) reports the selection",
  "script 08-feature-select.R" = c(
    "Choosing predictors without letting the choice flatter the score",
    "The honest score for a selected model"),
  "script 10-bayes.R" = "only worth the wait if the distribution is honest",
  "script example_nc_demo.R" = "which is anticonservative when points inside a cell are")

# What the pages said, each a statement of which way a score errs with no
# target, a name for it, or the sweep presented as the way to choose a block
# size.  None may stand on any page.
.CVW_BARRED <- c(
  "flatter autocorrelated", "without flattering", "comes back flattering",
  "the one to report", "optimistic bound", "still optimistic", "stays optimistic",
  "may leak", "leaky reference", "fully leaky", "below it leak",
  "become optimistic", "overstate performance", "overstate its accuracy",
  "overstates its accuracy", "overstated the model", "number overstated",
  "were hiding", "the general case", "two ways to leak", "ways a score can leak",
  "honest map-accuracy", "honest cv folds", "scored honestly",
  "for the folds to be honest", "therefore leaking", "whether the scheme leaks",
  "the direction that leaks", "which leaks", "avoid leakage",
  "cross-validation pessimistic", "validation pessimistic",
  "leakage mechanism", "predictors without leaking", "spatial structure left to leak",
  "stops being about leakage", "leaky inner loop", "dishonestly chosen",
  "question a map is used for", "question most maps are used for",
  "random folds overstated", "the safe side", "less conservative than it claims",
  "size blocks conservatively", "flatters a spatial model", "the optimism random folds",
  "the leak is gone", "honest blocking", "single honest number",
  "what the sweep is for", "defending a block size", "let the data choose the block size",
  "sizing it by measurement", "sizing a block by measurement", "blocks by measurement",
  "choosing a block size by cross-validated error", "should be setting your block size",
  "makes random cross-validation optimistic",
  # The run-time texts: the forest's print, its include_coords caution, the
  # random_kfold caution of select_features_forward() and the random-folds
  # fallback.
  "optimistic under spatial", "leak spatial autocorrelation", "honest-looking",
  "optimistic performance")

# The sentences each cv_*() help gives on what folds = NULL estimates.
.CVW_FOLDS_NULL <- paste(
  "With folds = NULL the score estimates the error at the distance the block",
  "grid built here puts between the held-out points and their training data,",
  "which overstates the error of a map where the samples cover the area it",
  "predicts over. The excess is how much the model's error grows from the",
  "map's distances to the data to the block grid's, and is small for a model",
  "whose error changes little with the distance from the data.")

# The sentences that give a direction, or say what a score estimates, word
# for word where they stand (squashed: no backticks, single spaces).
.CVW_PINNED <- list(
  "?spatialkit-package" = c(
    paste("Fold. make_folds() builds spatial cross-validation folds: blocked,",
          "buffered, leave-location-out or nearest-neighbour distance-matched.",
          "A cross-validated error is the error at the distance the folds put",
          "between the held-out points and their training data: it reads low",
          "for predictions made farther from the data and high for predictions",
          "made closer, so the folds to use are the ones matched to where the",
          "model will predict. fold_separation() with prediction_points",
          "compares the two distances."),
    "make_folds() (split the data into cross-validation folds)",
    paste("Sizing blocks from it is one way to set how far the held-out points",
          "sit from their training data, and suits a map that has to reach that",
          "far from the data.")),
  "?make_folds" = c(
    paste("A cross-validated error is the error at the distances the folds put",
          "between the held-out points and their training data. It reads low",
          "for predictions made farther from the data than the held-out points",
          "sit, and high for predictions made closer"),
    paste("\"random_kfold\" ignores geography: a held-out point keeps most of its",
          "neighbours in the training data, so the score is the error of",
          "predicting among the samples. That is the error of a map over an area",
          "the samples cover evenly, and it reads low for predictions made",
          "farther from the data, such as a map of a whole region from clustered",
          "samples."),
    paste("\"block_kfold\" holds out whole blocks, so the held-out points sit",
          "farther from their training data the wider the blocks are, and the",
          "score is the error of predicting into unsampled gaps of that size;",
          "for a map of an area the samples cover, it reads high."),
    "erring towards held-out distances longer than the target."),
  "?cv_gwr" = paste(
    "A random split leaves a held-out point's neighbours in the training",
    "set, so its score is the error of predicting among the samples;",
    "blocks hold out whole areas"),
  "?cv_rf" = "Built via block_kfold when NULL, as in cv_spatial().",
  "?compare_models_cv" = paste(
    "When it is NULL, one set of block_kfold folds is built and every model",
    "is scored on it."),
  "?fit_rf_model" = c(
    paste("Random folds, and the out-of-bag error, a random hold-out too,",
          "measure it among the samples: the error of a map over an area the",
          "samples cover evenly, which reads low for predictions made farther",
          "from the data."),
    paste("so under spatial autocorrelation it is the error of predicting among",
          "the samples, as random k-fold is, and reads low for predictions made",
          "farther from the data than the samples sit from each other")),
  "?print.rf_fit" = paste(
    "its out-of-bag error, a random hold-out, is its error among the samples",
    "and reads low for predictions made farther from the data;"),
  "?estimate_sac_range" = c(
    paste("On fields smoother than exponential that makes the range long, and",
          "blocks sized from it come out larger: they still hold the test points",
          "at least the field's range from their training data, which blocks",
          "sized from a short range do not."),
    paste("Blocks sized from that range are then narrower than the field's",
          "range, and hold the test points less far from their training data",
          "than range-sized blocks are meant to."),
    "or size its blocks with a margin above the range returned, to make up the shortfall.",
    paste("Treat the range from a count response as no more than an order of",
          "magnitude, since it can come out long or short")),
  "?fold_separation" = paste(
    "cv_block_size_sweep() for the cross-validated error at each block size.",
    "The sweep does not say which size matches a map; this function does,",
    "read on the folds of each size with the map's locations as",
    "prediction_points."),
  "vignette spatial-cross-validation" = c(
    "## The error at each block size",
    "{r separation}",
    paste("| random_kfold | a random 1/k of rows | the map predicts among the",
          "samples (an area they cover evenly), or the data are not",
          "autocorrelated |"),
    "fold_separation(fr, pts, prediction_points = cells)",
    paste("\"matched\" says the held-out points sit about as far from their",
          "training data as the map's cells do, so the score is the error of",
          "that map; \"farther\" says they sit farther, and the score reads high",
          "for it; \"closer\", the reading of random folds on clustered samples",
          "with a map of the whole region, says it reads low."),
    paste("The curve does not say which size to report; the size to read is the",
          "one whose folds match the map")),
  "vignette spatialkit_nc_demo" = c(
    "Which of the two to report depends on where the map will predict.",
    "{r cv-separation}",
    "fold_separation(folds_random, points_sf, prediction_points = cells)",
    paste("\"matched\" is the scheme whose score is the error of that map,",
          "\"farther\" one that reads high for it, and \"closer\" one that reads",
          "low."),
    "cv_block_size_sweep() gives the cross-validated",
    paste("What 5c measured was how much worse the coordinate forest predicts",
          "whole blocks of the state than it predicts among the points;")),
  "vignette reporting" = c(
    "{r block-vs-range, echo=FALSE, results='asis', eval=has_gstat}",
    "that is the error of the map depends on where the map predicts.",
    "fold_separation(folds, pts, sac = sac, prediction_points = surface)",
    "m(folds$params$block_size), m(sac), reading, says[[reading]]",
    paste("Read together, those lines say the run failed two of its own checks,",
          "the score and the residuals, which is what a report is for.")),
  "README" = "| check which way a CV score errs for my map | fold_separation(prediction_points = ) |",
  "script 03-folds.R" = c(
    "fold_separation(f, pts, prediction_points = map_pts)",
    paste("Neither number is the error of every map. Report the one whose folds",
          "match where the map will predict:")),
  "script 04-block-size.R" = c(
    "fold_separation(res[[nm]]$folds, pts, prediction_points = map_pts)",
    paste("The curve does not say which size to report: that is the size whose",
          "test points sit as far from the data as the locations the map will",
          "predict at")),
  "script example_nc_demo.R" = c(
    "fold_separation(f, points_sf, prediction_points = cells)",
    "Which of the two to report depends on where the map will predict"))


# The sentences of make_folds()'s five warnings after their first.
.CVW_READING <- paste(
  "Held-out points then have training data within the range, so the score is",
  "the error of predicting that close to the data: it reads low for",
  "predictions made farther from the data, and is the error of a map whose",
  "locations sit as close (fold_separation(prediction_points = ) compares the",
  "held-out distances with the map's).")

.cvw_warnings <- function(expr) {
  w <- character(0)
  withCallingHandlers(suppressMessages(expr), warning = function(cnd) {
    w <<- c(w, conditionMessage(cnd))
    invokeRestart("muffleWarning")
  })
  w
}

# An exponential field of effective range 210 on a 1000-unit square: the
# range estimated from it (about 200) lies between the default block edge at
# k = 5 (250) and at k = 10 (about 167).
.cvw_field <- function(n = 300, seed = 8) {
  set.seed(seed)
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000)
  d <- as.matrix(stats::dist(cbind(x, y)))
  z <- as.numeric(t(chol(exp(-d / 70) + diag(1e-8, n))) %*% rnorm(n)) + rnorm(n, sd = 0.3)
  sf::st_as_sf(data.frame(x = x, y = y, z = z, w = rnorm(n)), coords = c("x", "y"),
               crs = 32632)
}


# ---- The rules themselves -------------------------------------------------------------

test_that("the rules find a direction without a target, and pass one with a target", {
  bare <- "Random folds flatter autocorrelated data; these do not."
  expect_identical(.cvw_untargeted(bare), bare)
  expect_identical(.cvw_untargeted("Blocks below the range leak."), "Blocks below the range leak.")
  expect_identical(.cvw_untargeted("The score is still optimistic."), "The score is still optimistic.")
  expect_identical(.cvw_untargeted("Use cv_rf() for an honest figure."), "Use cv_rf() for an honest figure.")
  expect_identical(.cvw_untargeted("A pessimistic score."), "A pessimistic score.")
  for (bad in c("Random folds overstate the forest's accuracy.",
                "Blocks narrower than the range understate the error.",
                "Size blocks conservatively.",
                "A long range is the safe side.",
                "The blocked score reads high."))
    expect_identical(.cvw_untargeted(bad), bad, info = bad)
  # A target in the sentence: the locations the model will predict at.
  for (ok in c("It reads optimistic for predictions made farther from the data.",
               "Random folds flatter a map of the whole region.",
               "The folds stay optimistic for the prediction points.",
               "Blocks leak where the model will predict at the samples.",
               "It is pessimistic against the target.",
               "The score reads high for this map."))
    expect_identical(.cvw_untargeted(ok), character(0), info = ok)
  # A target in another sentence of the paragraph does not count.
  expect_identical(.cvw_untargeted("A map of the region is the target. Random folds overstate the accuracy."),
                   "Random folds overstate the accuracy.")
  # Letters inside another word are not the word, and no direction word is
  # no finding.
  expect_identical(.cvw_untargeted(c("An unleaked, dishonestly framed flat curve.",
                                     "The error grows with the distance from the data.")),
                   character(0))
  # A phrase that uses the word of something else is taken out first.
  expect_identical(.cvw_untargeted("It is leak-free; nested.", other = "leak-free"),
                   character(0))
  # The sentence found is returned with the phrase taken out.
  expect_identical(.cvw_untargeted("It is leak-free and leaky.", other = "leak-free"),
                   "It is   and leaky.")
  expect_identical(.cvw_paragraphs("One `a`\nline.\n\n  \nTwo\n\nThree"),
                   c("One a line.", "Two", "Three"))
  expect_identical(.cvw_sentences("One (see f()). Two, e.g. three. \"Four\" ends."),
                   c("One (see f()).", "Two, e.g. three.", "\"Four\" ends."))
  # A direction paired with the wrong distance.
  expect_length(.cvw_wrong_way("It reads low for predictions made farther from the data, and high for predictions made closer."), 0L)
  expect_length(.cvw_wrong_way("It reads high for predictions made farther from the data."), 1L)
  expect_length(.cvw_wrong_way("It reads low for a map that predicts closer to the data."), 1L)
  expect_length(.cvw_wrong_way("The held-out points sit farther than the cells do, so the score tends to read high."), 0L)
  expect_length(.cvw_wrong_way("The held-out points sit farther than the cells do, so the score reads low."), 1L)
  expect_length(.cvw_wrong_way("They stay closer to their training data, and the error reads high for them."), 1L)
  expect_length(.cvw_wrong_way("which understates the error of a map where the samples cover the area"), 1L)
  expect_length(.cvw_wrong_way("\"farther\" one that reads high for it, and \"closer\" one that reads low."), 0L)
  expect_length(.cvw_wrong_way("\"farther\" one that reads low for it, and \"closer\" one that reads high."), 2L)
  expect_length(.cvw_wrong_way("cat(\"\\\"farther\\\" tends to read low for it\")"), 1L)
  # A script's comment loses its marker, and the pieces of a printed
  # sentence are joined.
  dir <- tempfile("cvw"); dir.create(dir)
  writeLines(c("# A comment", "  #   indented", "cat(\"  one piece\\n\",", "    \"  and the next.\\n\", sep = \"\")"),
             file.path(dir, "s.R"))
  expect_identical(.cvw_squash(.cvw_script("s.R", dir)),
                   "A comment indented cat(\" one piece and the next.\\n\", sep = \"\")")
})


# ---- Every page ------------------------------------------------------------------------

test_that("no page says which way a cross-validated score errs without naming where the model will predict", {
  pages <- .cvw_pages()
  expect_gt(length(grep("^\\?", names(pages))), 75L)
  expect_true(all(paste("script", c("03-folds.R", "04-block-size.R", "08-feature-select.R",
                                    "example_nc_demo.R")) %in% names(pages)))
  if (.cvw_have_vignettes())
    expect_true(all(paste("vignette", c("spatial-cross-validation", "reporting",
                                        "spatialkit_nc_demo", "getting-started",
                                        "diagnostics", "resolution")) %in% names(pages)))
  for (what in names(pages)) {
    expect_false(is.null(pages[[what]]), info = what)
    paras <- .cvw_paragraphs(pages[[what]])
    expect_identical(.cvw_untargeted(paras, .CVW_OTHER_SENSES[[what]]), character(0),
                     info = what)
    expect_identical(.cvw_wrong_way(pages[[what]]), character(0), info = what)
    low <- tolower(.cvw_squash(pages[[what]]))
    for (phrase in .CVW_BARRED)
      expect_false(grepl(phrase, low, fixed = TRUE), info = paste(what, phrase))
  }
})

test_that("each phrase that uses a direction word of something else stands on its page", {
  pages <- .cvw_pages()
  for (what in names(.CVW_OTHER_SENSES)) {
    if (identical(what, "README") && is.null(pages[["README"]])) next
    if (startsWith(what, "vignette ") && !.cvw_have_vignettes()) next
    expect_false(is.null(pages[[what]]), info = what)
    text <- .cvw_squash(pages[[what]])
    for (phrase in .CVW_OTHER_SENSES[[what]])
      expect_true(grepl(phrase, text, fixed = TRUE), info = paste(what, phrase))
  }
})


# ---- The places the statements stood ----------------------------------------------------

test_that("each page that gave a direction without a target says instead what the score estimates, in the direction pinned", {
  pages <- .cvw_pages()
  for (what in names(.CVW_PINNED)) {
    if (identical(what, "README") && is.null(pages[["README"]])) next
    if (startsWith(what, "vignette ") && !.cvw_have_vignettes()) next
    expect_false(is.null(pages[[what]]), info = what)
    text <- .cvw_squash(pages[[what]])
    for (s in .CVW_PINNED[[what]])
      expect_identical(lengths(regmatches(text, gregexpr(s, text, fixed = TRUE))), 1L,
                       info = paste(what, s))
  }
  if (!.cvw_have_vignettes()) skip("the vignettes' sources are not installed")
  # The reporting vignette words its sentence by fold_separation()'s reading:
  # each reading has its own sentence, and each sentence its direction.
  code <- .cvw_chunk(.cvw_vignette("reporting"), "block-vs-range")
  expect_false(is.null(code))
  exprs <- parse(text = code)
  at <- which(vapply(exprs, function(e) is.call(e) && identical(e[[1L]], as.name("<-")) &&
                       identical(e[[2L]], as.name("says")), logical(1)))
  expect_length(at, 1L)
  says <- eval(exprs[[at]][[3L]], baseenv())
  expect_identical(names(says), c("matched", "closer", "farther", "mixed"))
  expect_match(says[["matched"]], "sit about as far from their training data as the cells do, so the score is the error of this map", fixed = TRUE)
  expect_match(says[["closer"]], "sit closer to their training data than the cells do, so the score tends to read low for this map; wider blocks", fixed = TRUE)
  expect_match(says[["farther"]], "sit farther from their training data than the cells do, so the score tends to read high for this map; narrower blocks", fixed = TRUE)
  expect_match(says[["mixed"]], "no single direction", fixed = TRUE)
  for (r in c("closer", "farther"))
    expect_match(says[[r]], "make_folds(method = \"knndm\", prediction_points = surface)", fixed = TRUE)
})

test_that("each cv_*() help says what folds = NULL estimates", {
  for (page in c("cv_spatial", "cv_gwr", "cv_bayes", "cv_rf", "compare_models_cv")) {
    folds <- .cvw_rd_arguments(.cvw_rd(page))[["folds"]]
    expect_false(is.null(folds), info = page)
    expect_match(folds, .CVW_FOLDS_NULL, fixed = TRUE, info = page)
    expect_match(folds, "fold_separation() with prediction_points", fixed = TRUE, info = page)
    # Once, in the entry for folds and nowhere else on the page.
    page_text <- .cvw_squash(.cvw_flat(.cvw_rd(page)))
    expect_identical(lengths(regmatches(page_text, gregexpr(.CVW_FOLDS_NULL, page_text,
                                                            fixed = TRUE))), 1L, info = page)
  }
})

test_that("?cv_block_size_sweep says what the rise of its curve measures", {
  d <- .cvw_rd_part(.cvw_rd("cv_block_size_sweep"), "\\description")
  expect_match(d, paste("The rise of the curve with block size measures the error at that",
                        "distance from the data: how much worse the model predicts a location",
                        "the farther it lies from the observations. It is the error of a map",
                        "only at locations that far from the data, so the curve does not by",
                        "itself say which size to report."), fixed = TRUE)
  expect_false(grepl("plateau", d, fixed = TRUE))
  expect_match(.cvw_rd_arguments(.cvw_rd("cv_block_size_sweep"))[["include_random"]],
               "Also cross-validate with random folds, as the reference: the error of predicting among the samples",
               fixed = TRUE)
  p <- .cvw_rd_part(.cvw_rd("plot.block_size_sweep"), "\\description")
  expect_match(p, paste("its height at a size is the error at that distance from the data,",
                        "and the error of a map only where the locations it predicts at sit",
                        "that far from the observations."), fixed = TRUE)
})


# ---- What the code says at run time ---------------------------------------------------

test_that("make_folds()'s five warnings on blocks below the range say what the score estimates, not that it may leak", {
  skip_if_not_installed("gstat")
  pts <- .cvw_field()
  r <- as.numeric(estimate_sac_range(pts, "z"))
  expect_gt(r, 170); expect_lt(r, 250)   # the geometry below depends on this
  blocks <- sf::st_sf(geometry = sf::st_make_grid(pts, n = c(8, 8)))
  # One call for each of the five warnings, run one at a time.
  calls <- list(
    "auto_range and a smaller block_size" = function()
      make_folds(pts, k = 5, method = "block_kfold", response_var = "z",
                 auto_range = TRUE, block_size = 100, seed = 1),
    "a smaller block_size" = function()
      make_folds(pts, k = 5, method = "block_kfold", response_var = "z",
                 block_size = 100, seed = 1),
    "the default grid" = function()
      make_folds(pts, k = 10, method = "block_kfold", response_var = "z", seed = 1),
    "block_nx and block_ny" = function()
      make_folds(pts, k = 5, method = "block_kfold", response_var = "z",
                 block_nx = 8, block_ny = 8, seed = 1),
    "supplied blocks" = function()
      make_folds(pts, k = 5, method = "block_kfold", response_var = "z",
                 blocks = blocks, seed = 1))
  heads <- c("make_folds\\(\\): block_size \\(100\\) < estimated autocorrelation range",
             "make_folds\\(\\): block_size \\(100\\) < estimated autocorrelation range",
             "make_folds\\(\\): block dimension \\([0-9.]+\\) < autocorrelation range",
             "make_folds\\(\\): block_nx/block_ny yield blocks smaller than autocorrelation range",
             "make_folds\\(\\): the supplied blocks \\(median scale [0-9.]+\\) are smaller than the autocorrelation range")
  for (i in seq_along(calls)) {
    w <- .cvw_warnings(calls[[i]]())
    w <- grep("autocorrelation range \\(", w, value = TRUE)
    expect_length(w, 1L)
    expect_match(w, paste0("^", heads[i], " \\([0-9.]+\\)\\. "), info = names(calls)[i])
    expect_false(grepl("leak", w, fixed = TRUE), info = names(calls)[i])
    expect_match(w, .CVW_READING, fixed = TRUE, info = names(calls)[i])
    if (i < 5L) {
      # The advice on the block size follows, as before.
      expect_match(w, paste0(.CVW_READING, " For the error farther from the data, widen the ",
                             "blocks. Pass block_size = "), fixed = TRUE, info = names(calls)[i])
    } else {
      expect_match(w, paste0(.CVW_READING, " For the error farther from the data, supply ",
                             "blocks at least "), fixed = TRUE, info = names(calls)[i])
      expect_match(w, "supply blocks at least [0-9]+ units across\\.$", info = names(calls)[i])
    }
  }
})

test_that("make_folds() logs the check of the blocks against the range as skipped, not a leakage check", {
  skip_if_not_installed("gstat")
  set.seed(4)
  pts <- sf::st_as_sf(data.frame(x = runif(25, 0, 1000), y = runif(25, 0, 1000), z = rnorm(25)),
                      coords = c("x", "y"), crs = 32632)
  lines <- capture_spatialkit_log(
    make_folds(pts, k = 4, method = "block_kfold", response_var = "z", block_size = 300,
               seed = 1),
    level = logger::INFO)
  expect_true(any(grepl(paste("make_folds(block_kfold): the check of the blocks against the",
                              "autocorrelation range is skipped: fewer than 30 points"),
                        lines, fixed = TRUE)))
  expect_false(any(grepl("leak", lines, fixed = TRUE)))
})

test_that("the plots of the folds and of the sweep no longer say that blocks below the range leak", {
  skip_if_not_installed("ggplot2")
  set.seed(21)
  n <- 160
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000)
  d <- as.matrix(stats::dist(cbind(x, y)))
  pts <- sf::st_as_sf(data.frame(x = x, y = y, a = rnorm(n)), coords = c("x", "y"),
                      crs = 32632)
  pts$z <- as.numeric(t(chol(exp(-d / 60) + diag(1e-6, n))) %*% rnorm(n)) +
    0.5 * pts$a + rnorm(n, 0, 0.2)
  msgs <- character(0)
  sw <- withCallingHandlers(
    cv_block_size_sweep(pts, "z", "a", fit_fn = function(tr) lm_spatial_fit(tr, "z", "a"),
                        k = 3, n_sizes = 2, sac = structure(120, class = "sac_range")),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  expect_true(any(grepl("cv_block_size_sweep(): random folds (the reference) ...", msgs,
                        fixed = TRUE)))
  expect_false(any(grepl("leak", msgs, fixed = TRUE)))
  p <- plot(sw)
  expect_identical(p$labels$subtitle, "Estimated autocorrelation range 120 (dotted)")
  expect_match(p$labels$caption, "^Dashed: random folds, RMSE = [0-9.e-]+\nBand: ")
  expect_false(grepl("leak", paste(p$labels$subtitle, p$labels$caption), fixed = TRUE))

  fb <- make_folds(pts, k = 4, method = "block_kfold", block_size = 300, seed = 1)
  expect_identical(strsplit(plot_folds(fb, pts)$labels$subtitle, "\n", fixed = TRUE)[[1L]][2L],
                   "Wider blocks hold points out farther from the data")
  fb$params$sac_range <- 250
  expect_identical(strsplit(plot_folds(fb, pts)$labels$subtitle, "\n", fixed = TRUE)[[1L]][2L],
                   "Estimated autocorrelation range 250")
})

# The texts the code prints or logs at run time: the forest's print, the
# caution of fit_rf_model(include_coords = TRUE), the caution of
# select_features_forward(method = "random_kfold") and the warning of the
# fallback to random folds.  Each is judged by the rules the pages are.

# A logged line starts with its level and time.
.cvw_drop_stamp <- function(x) sub("^[A-Z]+ \\[[^]]*\\] ", "", x)

# The texts that need ranger: the print of a forest without and with the
# coordinates and of one with no row out of bag, and the include_coords
# caution.
.cvw_forest_texts <- function() {
  pts <- .cvw_field(n = 120)
  plain <- fit_rf_model(pts, "z", "w", num_trees = 50L)
  rm(list = ls(spatialkit:::.warned_once), envir = spatialkit:::.warned_once)
  coords_log <- capture_spatialkit_log(
    with_xy <- fit_rf_model(pts, "z", "w", num_trees = 50L, include_coords = TRUE))
  # Every tree grown on every row: the fit warns, and logs it.
  no_oob_log <- capture_spatialkit_log(
    no_oob <- suppressWarnings(fit_rf_model(pts, "z", "w", num_trees = 50L,
                                            replace = FALSE, sample_fraction = 1)))
  list(
    "print of a forest" = utils::capture.output(print(plain)),
    "print of a forest with coordinates" = utils::capture.output(print(with_xy)),
    "print of a forest with no row out of bag" = utils::capture.output(print(no_oob)),
    "include_coords caution" = .cvw_drop_stamp(grep("include_coords = TRUE", coords_log,
                                                    value = TRUE)))
}

# The texts that need no forest: the random_kfold caution, with a linear
# fit, and the fallback to random folds.
.cvw_other_texts <- function() {
  pts <- .cvw_field(n = 120)
  select_log <- capture_spatialkit_log(
    select_features_forward(pts, "z", "w", function(tr, vars) lm_spatial_fit(tr, "z", vars),
                            k = 4, method = "random_kfold", seed = 1, quiet = TRUE))
  # The fallback is raised as an R warning and logged; the log is caught so
  # that it does not echo.
  fallback_log <- capture_spatialkit_log(
    fallback <- .cvw_warnings(spatialkit:::.remap_folds(NULL, 1:20)))
  list(
    "random_kfold caution" = .cvw_drop_stamp(grep("random_kfold inner folds", select_log,
                                                  value = TRUE)),
    "fallback to random folds" = fallback)
}

# What a text is judged in.  Above its first blank line a print is a header
# of fields, one to a line and with no full stop, so each line is a sentence
# of its own: run together, the header and the first sentence below it read
# as one sentence, and a target named there would pass a direction word on a
# header line.  Below the blank line the lines are paragraphs, as on a page.
# A logged line or a warning is one paragraph.
.cvw_runtime_units <- function(lines) {
  blank <- which(!nzchar(trimws(lines)))
  if (!length(blank)) return(.cvw_squash(paste(lines, collapse = " ")))
  head <- lines[seq_len(blank[1L] - 1L)]
  c(vapply(head, .cvw_squash, character(1), USE.NAMES = FALSE),
    .cvw_paragraphs(paste(lines[-seq_len(blank[1L])], collapse = "\n")))
}

# The rules of the pages, applied to each text.
.cvw_judge_runtime <- function(texts) {
  for (what in names(texts)) {
    # Each text is there to judge: a print, or one logged line or warning.
    if (startsWith(what, "print")) expect_gt(length(texts[[what]]), 5L, label = what)
    else expect_length(texts[[what]], 1L)
    units <- .cvw_runtime_units(texts[[what]])
    expect_identical(.cvw_untargeted(units), character(0), info = what)
    txt <- paste(units, collapse = " ")
    expect_identical(.cvw_wrong_way(txt), character(0), info = what)
    low <- tolower(txt)
    for (phrase in .CVW_BARRED)
      expect_false(grepl(phrase, low, fixed = TRUE), info = paste(what, phrase))
  }
}

test_that("a print is judged line by line above its closing paragraph", {
  # A direction word on a header line is found although the closing
  # paragraph names a target; read as one paragraph, the print passed it.
  probe <- c("<Random Forest (ranger)> spatial model fit",
             "  Coords as predictors: YES - see ?fit_rf_model (its OOB error is optimistic)",
             "  OOB RMSE: 1.2900   OOB R^2: 0.3500 (both flatter the forest)",
             "",
             "  OOB is a random hold-out: it reads low for predictions made",
             "  farther from the data. Score it with cv_rf().")
  expect_identical(.cvw_untargeted(.cvw_squash(paste(probe, collapse = " "))), character(0))
  units <- .cvw_runtime_units(probe)
  expect_identical(units,
                   c("<Random Forest (ranger)> spatial model fit",
                     "Coords as predictors: YES - see ?fit_rf_model (its OOB error is optimistic)",
                     "OOB RMSE: 1.2900 OOB R^2: 0.3500 (both flatter the forest)",
                     paste("OOB is a random hold-out: it reads low for predictions made farther",
                           "from the data. Score it with cv_rf().")))
  expect_identical(.cvw_untargeted(units), units[2:3])
  # Below the blank line a direction word needs its target in its own
  # sentence, as on a page.
  expect_identical(.cvw_untargeted(.cvw_runtime_units(c(probe[1L], "", "  Predictions made",
                                                        "  here. It reads low."))),
                   "It reads low.")
  # A logged line or a warning is one paragraph.
  expect_identical(.cvw_runtime_units("  One   line. It reads low.  "), "One line. It reads low.")
})

test_that("the forest's prints and its include_coords caution name where the model will predict", {
  skip_if_not_installed("ranger")
  texts <- .cvw_forest_texts()
  .cvw_judge_runtime(texts)
  # What each says instead, word for word.
  for (what in c("print of a forest", "print of a forest with coordinates"))
    expect_match(.cvw_squash(paste(texts[[what]], collapse = " ")),
                 paste("OOB is a random hold-out: under spatial autocorrelation it is the",
                       "error of predicting among the samples, and reads low for predictions",
                       "made farther from the data than the samples sit from each other.",
                       "Score a map that predicts farther with cv_rf() on folds at its",
                       "distances (?fold_separation)."), fixed = TRUE, info = what)
  # A forest with no row out of bag has no out-of-bag error to describe: its
  # print says how to score it instead.
  nob <- texts[["print of a forest with no row out of bag"]]
  expect_match(.cvw_squash(paste(nob, collapse = " ")),
               "OOB RMSE: undefined (no row is out of bag)", fixed = TRUE)
  expect_no_match(.cvw_squash(paste(nob, collapse = " ")), "OOB is a random hold-out",
                  fixed = TRUE)
  expect_identical(utils::tail(.cvw_runtime_units(nob), 1L),
                   paste("With no out-of-bag error to read, score the forest with cv_rf() on",
                         "folds at the distances the map will predict at (?fold_separation)."))
  expect_match(texts[["include_coords caution"]],
               paste("so its error grows quickly with the distance from the data. The",
                     "out-of-bag error and random cross-validation measure it among the",
                     "samples: the error of a map over an area the samples cover evenly,",
                     "which reads low for predictions made farther from the data."), fixed = TRUE)
  expect_match(texts[["include_coords caution"]],
               "Score the forest with cv_rf() on blocked or distance-matched folds at the distances the map will predict at",
               fixed = TRUE)
})

test_that("the random_kfold caution and the random-folds fallback name where the model will predict", {
  # Neither needs a forest, so both are judged without ranger.
  texts <- .cvw_other_texts()
  .cvw_judge_runtime(texts)
  expect_match(texts[["random_kfold caution"]],
               paste("select_features_forward(): random_kfold inner folds score each candidate",
                     "set by its error among the samples, so the sweep chooses the variables",
                     "that predict best there, where one that only tracks nearby values of the",
                     "response can win. A blocked outer loop then scores the chosen set at the",
                     "blocks' distance, for which it was not chosen. Use inner folds that hold",
                     "points out as the outer ones do: 'block_kfold', the default, inside",
                     "blocked outer folds."), fixed = TRUE)
  expect_identical(texts[["fallback to random folds"]],
                   paste("cross-validation: falling back to random k-fold CV, whose score is the",
                         "error of predicting among the samples. For a map that predicts",
                         "farther from the data, build folds at its distances with make_folds()",
                         "(fold_separation(prediction_points = ) compares them)."))
})


# ---- What the new sentences rest on ----------------------------------------------------

test_that("on the fixture of the cross-validation vignette, random folds match a map of the square and blocks read high for it", {
  skip_on_cran()
  # 400 sites on a 1000-unit square, an exponential field of effective range
  # 180 and noise of sd 0.3, the vignette's five-nearest-neighbour smoother.
  # The truth is the error of the map made from all 400 sites at the 900
  # centres of a 30 x 30 grid, drawn with the sites from one field.
  knn <- function(train_sf, k = 5L)
    new_spatial_fit(subclass = "cvw_knn",
                    engine = list(xy = sf::st_coordinates(train_sf), y = train_sf$z, k = k),
                    formula = z ~ 1, response_var = "z", predictor_vars = character(0),
                    data_sf = train_sf)
  pred <- function(object, newdata = NULL, ...) {
    e <- object$engine
    xy <- if (is.null(newdata)) e$xy else sf::st_coordinates(newdata)
    vapply(seq_len(nrow(xy)), function(i) {
      d <- (e$xy[, 1] - xy[i, 1])^2 + (e$xy[, 2] - xy[i, 2])^2
      mean(e$y[order(d)[seq_len(e$k)]])
    }, numeric(1))
  }
  registerS3method("predict", "cvw_knn", pred)
  registerS3method("fitted", "cvw_knn", function(object, ...) pred(object))
  registerS3method("residuals", "cvw_knn", function(object, ...) object$engine$y - pred(object))
  off <- (seq_len(30) - 0.5) * 1000 / 30
  grid_xy <- expand.grid(x = 5e5 + off, y = 5e6 + off)
  out <- lapply(c(7, 101, 102, 103), function(seed) {
    set.seed(seed)
    s  <- data.frame(x = 5e5 + runif(400, 0, 1000), y = 5e6 + runif(400, 0, 1000))
    xy <- rbind(s, grid_xy)
    z  <- as.numeric(t(chol(exp(-as.matrix(stats::dist(xy)) / 60) + diag(1e-8, nrow(xy)))) %*%
                       rnorm(nrow(xy))) + rnorm(nrow(xy), sd = 0.3)
    pts   <- sf::st_as_sf(cbind(s, z = z[1:400]), coords = c("x", "y"), crs = 32632)
    cells <- sf::st_as_sf(cbind(grid_xy, z = z[-(1:400)]), coords = c("x", "y"), crs = 32632)
    truth <- sqrt(mean((cells$z - pred(knn(pts), cells))^2))
    fr <- make_folds(pts, k = 5, method = "random_kfold", seed = 1)
    fb <- make_folds(pts, k = 5, method = "block_kfold", block_size = 200, seed = 1)
    score <- function(...) suppressWarnings(suppressMessages(
      cv_spatial(pts, "z", character(0), fit_fn = function(tr, ...) knn(tr), ...)))$overall$RMSE
    reading <- function(f) attr(fold_separation(f, pts, prediction_points = cells),
                                "target")$reading
    c(random = score(folds = fr) / truth - 1, blocked = score(folds = fb) / truth - 1,
      default = score(k = 5, seed = 1) / truth - 1,
      random_matched = reading(fr) == "matched", blocked_farther = reading(fb) == "farther")
  })
  out <- do.call(rbind, out)
  # Random folds are the error of this map; the blocks the vignette uses, and
  # the folds cv_spatial() builds with folds = NULL, overstate it.
  expect_lt(abs(mean(out[, "random"])), 0.10)
  expect_gt(mean(out[, "blocked"]), 0.15)
  expect_gt(mean(out[, "default"]), 0.15)
  expect_true(all(out[, "random_matched"] == 1))
  expect_true(all(out[, "blocked_farther"] == 1))
})
