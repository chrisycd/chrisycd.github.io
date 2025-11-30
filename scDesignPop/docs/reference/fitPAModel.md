# Fit a marginal model for power analysis

Fit a marginal model for power analysis

## Usage

``` r
fitPAModel(
  df,
  model_formula,
  idx,
  method = c("nb", "poisson", "gaussian", "pseudoBulkLinear"),
  snpid,
  indiv_colname
)
```

## Arguments

- df:

  a data frame object contains the design matrix.

- model_formula:

  a stats::formula object contains the model formula for power analysis.

- idx:

  a numeric value recording the serial number of the simulation

- method:

  a character object specifying the method that will be analyzed for
  power. (Options: nb,poisson,gaussian,pseudoBulkLinear).

- snpid:

  a character object contains snpid.

- indiv_colname:

  a string scalar of the sample ID variable in cell covariate of
  `marginal_list[[geneid]]$frame`. The default is "indiv".

## Value

a fitted stats::model object.

## Examples

``` r
NULL
#> NULL
```
