# Generic function to compute model parameter vectors

A S3 generic function for computing model parameters for with or without
new covariate of a feature

## Usage

``` r
calcParaVectors(fit, family_use, new_covariate, total_cells, data, ...)
```

## Arguments

- fit:

  a fitted object in the marginal_list.

- family_use:

  a string scalar or vector of marginal distribution used.

- new_covariate:

  a cell-by-covariate data frame obtained in the list output from
  [`constructDataPop`](https://github.com/chrisycd/scDesignPop/reference/constructDataPop.md).
  It must have a corr_group variable.

- total_cells:

  a positive integer for the number of total cells to simulate.

- data:

  a cell-by-covariate data frame obtained in the list output from
  [`constructDataPop`](https://github.com/chrisycd/scDesignPop/reference/constructDataPop.md).
  It must have a corr_group variable. Used only in gamlss fits.

- ...:

  additional arguments passed to calcParaVectors S3 method functions.
