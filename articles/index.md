# Articles

### Get started

- [Get
  started](https://elkronos.github.io/gis_modeling_toolkit/articles/getting-started.md):

  Take a layer of points to cells with standard errors, a score from
  spatial cross-validation and a map that marks where its predictions
  apply.

### Step by step

One step of the pipeline per article, in the order you take them.

- [Choosing a
  resolution](https://elkronos.github.io/gis_modeling_toolkit/articles/resolution.md):

  Decide how many cells the data support, and pass that count to the
  tessellation.

- [Building cells and summarising
  them](https://elkronos.github.io/gis_modeling_toolkit/articles/cells.md):

  Draw the cells, put each observation in one, and report each cell’s
  value with the right standard error.

- [Spatial
  cross-validation](https://elkronos.github.io/gis_modeling_toolkit/articles/spatial-cross-validation.md):

  Choose folds whose score matches the map you will make, size the
  blocks, and read every part of the result.

- [Fitting and checking
  models](https://elkronos.github.io/gis_modeling_toolkit/articles/models.md):

  Pick a backend or write your own, compare backends on the same folds,
  check the residuals, and map where the predictions apply.

- [Exporting and reporting
  results](https://elkronos.github.io/gis_modeling_toolkit/articles/reporting.md):

  Write the cells to a file others can join to, and report the numbers
  that say how far a cell value or a map can be trusted.

### Worked example

- [Worked example: North
  Carolina](https://elkronos.github.io/gis_modeling_toolkit/articles/spatialkit_nc_demo.md):

  The pipeline on a simulated field over North Carolina, comparing four
  tessellations, two fold schemes and two backends at cell counts set by
  hand.

### Help

Pages to consult when a step fails or a term is unclear.

- [Troubleshooting and
  performance](https://elkronos.github.io/gis_modeling_toolkit/articles/troubleshooting.md):

  The errors and warnings you are most likely to meet, what each means,
  and how to handle large layers, parallel folds, caches and the log.

- [Terms used on these
  pages](https://elkronos.github.io/gis_modeling_toolkit/articles/glossary.md):

  Short definitions of the terms the articles and help pages use.
