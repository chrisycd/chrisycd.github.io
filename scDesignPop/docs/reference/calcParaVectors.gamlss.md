# A calcParaVectors Method for gamlss Objects

A calcParaVectors Method for gamlss Objects

## Usage

``` r
# S3 method for class 'gamlss'
calcParaVectors(fit, family_use, new_covariate, data, ...)
```

## Arguments

- fit:

  a fitted object in the marginal_list.

- family_use:

  a string scalar or vector of marginal distribution used.

- new_covariate:

  a cell-by-covariate data frame obtained in the list output from
  [`constructDataPop()`](https://chrisycd.github.io/scDesignPop/reference/constructDataPop.md).
  It must have a corr_group variable.

- data:

  a cell-by-covariate data frame used to fit marginal models in
  [`fitMarginalPop()`](https://chrisycd.github.io/scDesignPop/reference/fitMarginalPop.md).
  It must have a corr_group variable.

- ...:

  additional arguments passed to calcParaVectors S3 method functions.
