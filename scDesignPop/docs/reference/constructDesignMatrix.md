# Construct a Design Matrix Dataframe

A helper function that constructs a design matrix using a given
feature's expression vector for every cell, eQTL genotype dataframe
(optional), and cell covariate dataframe.

## Usage

``` r
constructDesignMatrix(
  response_vec,
  cellcov_df,
  eqtlgeno_df,
  loc_colname = "POS",
  snp_colname = "snp_id",
  indiv_colname = "indiv",
  filter_snps = TRUE,
  snpvar_thres = 0,
  cleanup = TRUE
)
```

## Arguments

- response_vec:

  a vector of values for the response variable.

- cellcov_df:

  a cell-by-covariates dataframe containing the covariates (explanatory
  variables) for all cells.

- eqtlgeno_df:

  a SNP-by-sample genotype dataframe containing a feature's eQTL
  annotations and SNP genotypes (explanatory variables) for all samples
  (ie. individuals).

- loc_colname:

  a string scalar for column name of SNP position variable.

- snp_colname:

  a string scalar for column name of SNP id variable.

- indiv_colname:

  a string scalar for column name of individuals (samples).

- filter_snps:

  a logical scalar for whether to filter out SNP covariates with either
  low-variance or with only 1 distinct genotype (ie. all 1's) prior to
  fitting the model.

- snpvar_thres:

  a numeric scalar (between 0 to 1) used to filter out SNPs whose
  variance of genotypes across samples are below this threshold. Used
  together when `filter_snps = TRUE`.

- cleanup:

  a logical scalar for whether to clean up variables after constructing
  `dmat_df`.

## Value

a list containing the following:

- `dmat_df`:

  a dataframe of the design matrix containing all covariates (both cell
  covariates and eQTL genotype covariates) for a given feature.

- `snp_cov`:

  a string scalar or vector of SNP ids in the design matrix.

## Details

When `response_vec` is provided, the function constructs a full design
matrix dataframe (ie. with response variable), whereas only a covariate
dataframe is constructed when `response_vec = NULL`.

## Examples

``` r
NULL
#> NULL
```
