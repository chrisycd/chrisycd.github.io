# Construct a list of input data

Function extracts an expression matrix, cell covariates, and filters the
SNPs in the eQTL genotype dataframe.

## Usage

``` r
constructDataPop(
  sce,
  eqtlgeno_df,
  new_eqtlgeno_df = NULL,
  new_covariate = NULL,
  overlap_features = NULL,
  sampid_vec = NULL,
  copula_variable = celltype_colname,
  n_quantiles = 10,
  slot_name = "counts",
  snp_mode = c("single", "multi"),
  time_colname = NULL,
  celltype_colname = "cell_type",
  feature_colname = "gene_id",
  snp_colname = "snp_id",
  loc_colname = "POS",
  chrom_colname = "CHR",
  indiv_colname = "indiv",
  prune_thres = 0.9,
  ...
)
```

## Arguments

- sce:

  a SingleCellExperiment object.

- eqtlgeno_df:

  a dataframe with eQTL annotations and SNP genotypes for each gene.

- new_eqtlgeno_df:

  a dataframe with eQTL annotations and SNP genotypes for each gene in
  new individuals. The default is NULL.

- new_covariate:

  a cell covariate dataframe for which to simulate. The default is NULL.

- overlap_features:

  an optional string vector to filter for features (ie. genes). The
  default is NULL.

- sampid_vec:

  an optional string vector to filter for sample ids. The default is
  NULL.

- copula_variable:

  a string scalar specifying the column name in cell covariate of `sce`
  object for which to fit the copula by. If NULL, then copula will be
  fitted for all cells. if it's the same as celltype_colname, then
  copula will be fitted per cell type, while if it's the same as
  time_colname, copula will be fitted per time quantile.

- n_quantiles:

  an integer value used for how many time quantiles of continuous cell
  states will be assigned for copula fitting. The default is 10.

- slot_name:

  a string scalar specifying the slot to use in input `sce`. The default
  is "counts".

- snp_mode:

  a string scalar specifying the type of SNP mode used in the marginal
  models. Options are either "single" for single-SNP mode, or "multi"
  for multi-SNP mode. If multi-SNP mode is used, each gene will have one
  or more SNPs in the model upon SNP pruning. The default is "single".

- time_colname:

  a string scalar specifying the time/pseudotime variable in the cell
  covariate of `sce` object.

- celltype_colname:

  a string scalar specifying the cell type variable in `eqtlgeno_df` and
  cell covariate of `sce` object. The default is "cell_type".

- feature_colname:

  a string scalar specifying the feature variable (ie. genes) in
  `eqtlgeno_df`. The default is "gene_id".

- snp_colname:

  a string scalar for the SNP variable in `eqtlgeno_df`. The default is
  "snp_id".

- loc_colname:

  a string scalar for the last column of eQTL annotation in
  `eqtlgeno_df`. The default is "POS".

- chrom_colname:

  a string scalar of the chromosome variable in `eqtlgeno_df`. The
  default is "CHR".

- indiv_colname:

  a string scalar of the sample ID variable in cell covariate of `sce`.
  The default is "indiv".

- prune_thres:

  a numerical value between 0 and 1 used to threshold the pairwise
  correlations of eQTLs' genotypes for each feature. The default value
  is 0.9.

- ...:

  additional arguments passed to internal functions.

## Value

outputs a list with following elements:

- `count_mat`:

  a cell-by-gene matrix of response values.

- `covariate`:

  a cell-by-covariate data frame used for fit marginal.

- `new_covariate`:

  an optional cell-by-covariate data frame used in simulation.

- `important_features`:

  a string vector of gene ids.

- `noeqtl_features`:

  a string vector of gene ids that do not have eQTLs.

- `eqtl_geno_list`:

  a list of eQTL genotype dataframes for each gene (for fit marginal).

- `new_eqtl_geno_list`:

  a optional list of eQTL genotype dataframes for each gene (for new
  individual simulation).

- `filtered_gene`:

  string vector of features QC filtered.

## Examples

``` r
NULL
#> NULL
```
