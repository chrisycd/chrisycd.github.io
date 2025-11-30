# Simulate new design matrix for power analysis

Simulate new design matrix for power analysis

## Usage

``` r
simulatePADesignMatrix(
  fit,
  df_sel,
  nindiv_total,
  model = c("nb", "poisson", "gaussian"),
  snpid,
  nindiv,
  ncell,
  celltype_colname,
  indiv_colname
)
```

## Arguments

- fit:

  a fitted stats::model object.

- df_sel:

  a data frame contains the new design matrix.

- nindiv_total:

  a vector contains the number of individuals for each genotype.

- model:

  a character showing the model types of the full marginal model.

- snpid:

  a character object contains snpid.

- nindiv:

  a numeric value showing the number of individuals that user wants to
  simulate.

- ncell:

  a numeric value showing the number of cells per each individual that
  user wants to simulate.

- celltype_colname:

  a string scalar specifying the cell state variable in `df_sel`. The
  default is "cell_type".

- indiv_colname:

  a string scalar of the sample ID variable in cell covariate of
  `marginal_list[[geneid]]$frame`. The default is "indiv".

## Value

a new data frame contains the design matrix with simulated response.

## Examples

``` r
NULL
#> NULL
```
