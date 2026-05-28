# A calcParaVectors Method for glmmTMB Objects

A calcParaVectors Method for glmmTMB Objects

## Usage

``` r
# S3 method for class 'glmmTMB'
calcParaVectors(fit, family_use, new_covariate, ...)
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

- ...:

  additional arguments passed to calcParaVectors S3 method functions.
