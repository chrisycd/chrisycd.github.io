# Extract data from SingleCellExperiment object

Extract data from SingleCellExperiment object

## Usage

``` r
extractFromSCE(
  sce_obj,
  features = NULL,
  slot_name = "counts",
  indiv_colname = "indiv",
  ...
)
```

## Arguments

- sce_obj:

  a SingleCellExperiment object

- features:

  string vector to filter for features

- slot_name:

  string scalar to specify the type of input assay data

- indiv_colname:

  string scalar to specify the column that has the indiv ids

- ...:

  other options to specify extracted outputs

## Value

a list of extracted data

## Examples

``` r
NULL
#> NULL
```
