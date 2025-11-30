# Constructs a Model Formula

A helper function that specifies the formula for the marginal model.

## Usage

``` r
constructFormula(
  model_formula,
  interact_colnames = NULL,
  snp_cov = NULL,
  force_formula = FALSE
)
```

## Arguments

- model_formula:

  a string scalar to specify the model formula, including random effects
  (if any) but without SNP genotypes or SNP genotype interaction
  effects.

- interact_colnames:

  a string scalar or vector for the variable names that have first-order
  interaction with SNP genotypes.

- snp_cov:

  a string scalar or vector that specifies the SNP covariates used in
  the model.

- force_formula:

  a logical scalar for whether to bypass model parsimony check. If
  `force_formula = TRUE`, interaction terms whose covariates are not
  main effects in the model would be permitted. Results in error if
  `force_formula = FALSE` and `length(geno_interact_names) > 0`.

## Value

a formula object

## Examples

``` r
NULL
#> NULL
```
