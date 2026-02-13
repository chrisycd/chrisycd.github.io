# Construct a SingleCellExperiment object with specified covariates.

Construct a SingleCellExperiment object with specified covariates.

## Usage

``` r
constructSCE(
  data_obj,
  cellcov_df = NULL,
  featcov_df = NULL,
  cellcov_names,
  sampid_vec = NULL,
  overlap_features = NULL,
  assay_name = "RNA",
  slot_name = "counts",
  sce_name = slot_name,
  indiv_colname = "indiv",
  celltype_colname = "cell_type",
  factor_colnames = NULL,
  cellcov_renames = cellcov_names
)
```

## Arguments

- data_obj:

  a Seurat, SingleCellExperiment, or expression matrix (gene-by-cell)
  used to construct the new SingleCellExperiment object.

- cellcov_df:

  a dataframe containing covariates for all cells. Row names must have
  cell names that match `data_obj`. Only used when `data_obj` is a
  matrix object.

- featcov_df:

  a dataframe containing covariates for all features (ie genes). Row
  names must have feature names that match `data_obj`. Only used when
  `data_obj` is a matrix object.

- cellcov_names:

  a string vector of all columns to extract from cell covariates.

- sampid_vec:

  an optional string vector of sample ids (individuals). Default is
  NULL. If provided, it must match sample ids found in `data_obj`.

- overlap_features:

  an optional string vector for features to filter `data_obj`. Default
  is NULL. If provided, all features must be present in the `data_obj`.

- assay_name:

  a string scalar of name of the Seurat object slot of interest. Only
  used if `data_obj` is a Seurat object.

- slot_name:

  a string scalar for the type of assay data (either 'logcounts' or
  'counts'). Default is 'counts'. Only used if `data_obj` is a Seurat or
  SingleCellExperiment object.

- sce_name:

  a string scalar for the type of assay data for output
  SingleCellExperiment.

- indiv_colname:

  a string scalar for the sample id colname. Must be present in
  `cellcov_df` or `data_obj`. Default is 'indiv'.

- celltype_colname:

  a string scalar for the cell type colname. Must be present in
  `cellcov_df` or `data_obj`. Default is 'cell_type'.

- factor_colnames:

  an optional string vector or scalar for the columns to coerce into
  factor variables. Default is NULL. If provided, it must be present in
  `cellcov_df` or `data_obj`.

- cellcov_renames:

  an optional string vector to rename the cell covariates. Must be same
  length as `cellcov_names`. Default is same as `cellcov_names`.

## Value

a SingleCellExperiment object

## Examples

``` r
NULL
#> NULL
```
