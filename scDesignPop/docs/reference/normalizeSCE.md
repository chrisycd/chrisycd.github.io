# Normalize data in a SingleCellExperiment object

Normalizes the data in a SingleCellExperiment object with specified
method and returns a SingleCellExperiment object.

## Usage

``` r
normalizeSCE(
  sce,
  method = c("log1p", "seurat", "scran", "none"),
  overwrite = FALSE,
  ...
)
```

## Arguments

- sce:

  A SingleCellExperiment object.

- method:

  A string value specifying the normalization method. Must be one of
  'log1p', seurat', 'scran', or 'none'. The default is 'log1p'. See
  Normalization details.

- overwrite:

  A logical value on whether to overwrite the "logcounts" slot in input
  sce object. Must be set to `overwrite = TRUE` to overwrite. The
  default is FALSE.

- ...:

  Additional arguments passed to internal functions.

## Value

A SingleCellExperiment object

## Details

### Normalization options

If 'seurat' is used, the default options from `NormalizeData` from the
`Seurat` package is used. Seurat v4 is currently supported; If 'scran'
is used, the default options from `computePooledFactors`, followed by
`logNormCounts` from the `scuttle` package is used. If 'log1p' is used,
the `log1p` function from base R is used. In general, the computational
time from fast to slow is 'log1p', 'seurat', and 'scran'.

## Examples

``` r
data("example_sce")
example_sce <- normalizeSCE(example_sce, method = "log1p", overwrite = TRUE)
```
