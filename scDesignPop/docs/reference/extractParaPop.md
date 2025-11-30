# Extract parameter matrix for new covariate df

This is the main function.

## Usage

``` r
extractParaPop(
  sce,
  assay_use = "counts",
  marginal_list,
  n_cores,
  family_use,
  new_covariate,
  new_eqtl_geno_list,
  indiv_colname = "indiv",
  snp_colname = "snp_id",
  loc_colname = "POS",
  parallelization = "mcmapply",
  BPPARAM = NULL,
  data
)
```

## Arguments

- sce:

  add later

- assay_use:

  add later

- marginal_list:

  add later

- n_cores:

  add later

- family_use:

  a string scalar or vector of marginal distribution used.

- new_covariate:

  a cell-by-feature covariate dataframe (from construct_data.R) plus
  corr_group.

- new_eqtl_geno_list:

  a list of eQTL genotype dataframes for each gene (to be predicted).

- indiv_colname:

  add later

- snp_colname:

  add later

- loc_colname:

  add later

- parallelization:

  add later

- BPPARAM:

  add later

- data:

  a cell-by-feature covariate dataframe (from construct_data.R) plus
  corr_group. Used only in gamlss fits.

## Value

a list of mean, sigma, and zero parameter cell by feature matrices:

## Examples

``` r
NULL
#> NULL
```
