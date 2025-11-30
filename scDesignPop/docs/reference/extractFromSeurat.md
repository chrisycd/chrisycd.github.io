# Extract data from Seurat object

Extract data from Seurat object

## Usage

``` r
extractFromSeurat(
  seurat_obj,
  features = NULL,
  assay_name = "RNA",
  slot_name = "counts",
  indiv_colname = "indiv",
  ...
)
```

## Arguments

- seurat_obj:

  a Seurat object

- features:

  string vector to filter for features

- assay_name:

  string scalar to specify the name of input assay data

- slot_name:

  string scalar to specify the type of input assay data

- indiv_colname:

  string scalar to specify the column that has the indiv ids

- ...:

  other options to specify extracted outputs. Currently available:

  - `col_class` : logical value

  - `feat_names` : logical value

  - `cell_names` : logical value

  - `sc_indiv` : logical value

  - `expr_mat` : logical value

  - `sc_cov` : logical value

  - `feat_cov` : logical value

## Value

a list of extracted data with following optional outputs:

- `col_class` : a vector of the column classes for cell covariate
  variables

- `feat_names` : a vector of feature names

- `cell_names` : a vector of cell names

- `sc_indiv` : a vector of distinct sample ids extracted from cell
  covariate using `indiv_colname`

- `expr_mat` : an expression matrix extracted using `assay_name` and
  `slot_name` options

- `sc_cov` : a dataframe for cell covariates

- `feat_cov` : a dataframe for feature covariates

## Examples

``` r
NULL
#> NULL
```
