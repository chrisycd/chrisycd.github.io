# Generic function to compute model parameter vectors

A S3 generic function for computing the mean, theta, zero parameter
vectors with covariates for a feature, given the parametric family of
the marginal model.

## Usage

``` r
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

## Details

It also uses either train or new data, and has option to compute
parameter vectors for new individuals options. Note that there is
randomness introduced when computing the parameter vectors for new
individuals, as there are random effects.
