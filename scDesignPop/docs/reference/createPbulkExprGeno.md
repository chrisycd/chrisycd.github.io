# Creates pseudobulk expression grouped by eQTL genotypes.

Computes pseudobulk based on normalization and aggregation options, and
combines with eQTL genotypes to output dataframe of pseudobulk
expression vs genotypes.

## Usage

``` r
createPbulkExprGeno(
  sce_list,
  eqtlgeno,
  feature_sel,
  celltype_sel = "allct",
  eqtl_indx = 1L,
  eqtl_snp = NULL,
  normalize_type = c("log1p", "seurat", "scran", "none"),
  aggregate_type = c("mean", "sum"),
  slot_name = c("counts", "logcounts"),
  rescale = FALSE,
  overwrite = FALSE,
  feature_colname = "gene_id",
  snp_colname = "snp_id",
  indiv_colname = "indiv",
  celltype_colname = "cell_type",
  loc_colname = "POS",
  if_plot = TRUE,
  ...
)
```

## Arguments

- sce_list:

  A list of named SingleCellExperiment (SCE) objects. Each SCE object
  must have the `indiv_colname` and `celltype_colname` columns in its
  colData. If no names are provided in the list, the objects will be
  named "sce_1", "sce_2", etc, by default.

- eqtlgeno:

  A dataframe of eQTL with genotypes for samples. Must have the columns
  corresponding to `feature_colname`, `snp_colname`, `celltype_colname`,
  and `loc_colname`. The `loc_colname` variable must be the last column
  that precedes the genotypes of the individuals (samples).

- feature_sel:

  A string scalar for the feature name (eg. a gene id) to select. Must
  exist in rowname of every SCE object in `sce_list`.

- celltype_sel:

  A string scalar or vector of the cell states (eg. cell types) to
  select. Must be one of the cell types in the `celltype_colname` of
  every SCE object in `sce_list`.

- eqtl_indx:

  An integer scalar or vector to select the eQTL from `eqtlgeno` for a
  given celltype and feature. Only used if `eqtl_snp = NULL`. The
  default is 1. The values correspond to `celltype_sel` and are recycled
  if the length of elements is shorter.

- eqtl_snp:

  A string scalar or vector to select the eQTL from `eqtlgeno`.
  Overrides `eqtl_indx` option if both `eqtl_indx` and `eqtl_snp` are
  specified. Can select any SNP in the `eqtlgeno` (even from another
  feature). Must exist in the `snp_colname` column of `eqtlgeno`. The
  values correspond to `celltype_sel` and are recycled if the length of
  elements is shorter.

- normalize_type:

  A string scalar for normalization method applied to `sce_list`.

- aggregate_type:

  A string scalar for pseudobulk aggregation applied to `sce_list` after
  subsetting by `celltype_sel`.

- slot_name:

  A string scalar for the SCE slot used. The default is "counts".

- rescale:

  A logical scalar for whether to rescale the response values based on
  maximum expression value.

- overwrite:

  A logical value on whether to overwrite the "logcounts" slot in input
  sce object. Must be set to `overwrite = TRUE` to overwrite. The
  default is FALSE.

- feature_colname:

  A string scalar for column name of feature variable (i.e. gene). The
  default is "gene_id".

- snp_colname:

  A string scalar for column name of the SNP variable. The default is
  "snp_id".

- indiv_colname:

  A string scalar for column name of individuals (samples). The default
  is "indiv".

- celltype_colname:

  A string scalar for column name of cell type. The default is
  "cell_type".

- loc_colname:

  A string scalar for column name of SNP position variable. The default
  is "POS".

- if_plot:

  A logical scalar specifying whether to return a pseudobulk plot. The
  default is TRUE.

- ...:

  Additional arguments passed to internal functions.

## Value

Outputs a dataframe with following columns:

- `indiv`:

  Name of individual's sample ID.

- `response`:

  Pseudobulk expression of given feature.

- `genotype`:

  Genotype values of individual for selected eQTL SNP.

- `feature`:

  Name of selected feature.

- `celltype`:

  Name of selected cell type.

- `snp`:

  Name of selected eQTL SNP.

- `sce_name`:

  Name of the SingleCellExperiment object.

- `normalization`:

  Type of normalization method used.

- `aggregate_type`:

  Type of pseudobulk aggregation used.

- `slot_used`:

  Slot used from input SingleCellExperiment object.

## Examples

``` r
NULL
#> NULL
```
