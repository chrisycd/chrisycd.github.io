# Fit a Marginal Model

Fits a specified parametric model for a feature using a response
variable, and eQTL genotype and cell covariates as explanatory
variables.

## Usage

``` r
fitModel(
  feature_name,
  response_vec,
  cellcov_df,
  eqtlgeno_df,
  mu_formula,
  model_family = "nb",
  interact_colnames = NULL,
  loc_colname = "POS",
  snp_colname = "snp_id",
  cellstate_colname = "cell_type",
  indiv_colname = "indiv",
  filter_snps = TRUE,
  snpvar_thres = 0,
  force_formula = FALSE
)
```

## Arguments

- feature_name:

  a string scalar of a feature's name (ie. gene id).

- response_vec:

  a vector of values for the response variable.

- cellcov_df:

  a cell-by-covariates dataframe containing the covariates (explanatory
  variables) for all cells.

- eqtlgeno_df:

  a SNP-by-sample genotype dataframe containing a feature's eQTL
  annotations and SNP genotypes (explanatory variables) for all samples
  (ie. individuals).

- mu_formula:

  a string scalar to specify the mean formula, including random effects
  (if any) but without SNP genotypes or SNP genotype interaction
  effects.

- model_family:

  a string scalar to specify model fitting used.

- interact_colnames:

  a string scalar or vector for the variable names that have first-order
  interaction with SNP genotypes.

- loc_colname:

  a string scalar for column name of SNP position variable.

- snp_colname:

  a string scalar for column name of SNP id variable.

- cellstate_colname:

  a string scalar for column name of cell state (ie. cell type).

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

- force_formula:

  a logical scalar for whether to bypass model parsimony check. If
  `force_formula = TRUE`, interaction terms whose covariates are not
  main effects in the model would be permitted. Results in error if
  `force_formula = FALSE` and `length(geno_interact_names) > 0`.

## Value

a list containing the fitted model object, elapsed time, SNP ids of
covariates, and removed cells (NOT currently implemented.)

## Examples

``` r
NULL
#> NULL
```
