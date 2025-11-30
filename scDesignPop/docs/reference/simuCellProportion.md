# Simulate cell proportions using genotype principal components and population-level covariates

Function fits a multinomial regression model using the input genotype
principal components (PCs) and other population-level covariates as
training data and outputs simulated cell proportions.

## Usage

``` r
simuCellProportion(
  sce,
  genoPC,
  new_genoPC,
  new_othercov,
  PCnum = 5L,
  cov_colnames = NULL,
  indiv_colname = "indiv",
  cellstate_colname = "cell_type",
  cn_model_family = "lognormal",
  cn_meanlog = NULL,
  cn_sdlog = NULL,
  cp_model_family = "MN",
  cp_intercept = TRUE,
  ...
)
```

## Arguments

- sce:

  a SingleCellExperiment object.

- genoPC:

  a data frame of individual by genotype principal components for `sce`
  input. The first column must be the variable for individual with same
  name as `indiv_colname`.

- new_genoPC:

  a data frame of individual by genotype principal components for
  simulated individuals. The first column must be the variable for
  individual with same name as `indiv_colname`, followed by "PC1",
  "PC2", etc.

- new_othercov:

  a data frame of the test data containing same additional covariates as
  in colData of `sce`.

- PCnum:

  an integer scalar specifying the number of principal components used
  in multinomial regression.

- cov_colnames:

  an optional string vector or scalar for the variable names to include
  in the cell proportion model. Variables must exist in both
  `new_othercov` and colData of `sce`.

- indiv_colname:

  a string scalar to specify the variable in `sce` containing
  individuals.

- cellstate_colname:

  a string scalar to specify the variable in `sce` containing cell
  states (ie. cell types).

- cn_model_family:

  a string scalar to specify the model family used for total cell
  modeling. Currently only 'lognormal' from
  [fitdistr](https://rdrr.io/pkg/MASS/man/fitdistr.html) is supported.

- cn_meanlog:

  a numeric scalar for the mean parameter (on log scale) of the total
  cell number model. When `cn_meanlog = NULL`, the parameter is
  estimated from input data.

- cn_sdlog:

  a numeric scalar for the standard deviation parameter (on log scale)
  of the total cell number model. When `n_sdlog = NULL`, the parameter
  is estimated from input data.

- cp_model_family:

  a string scalar to specify the model family used for cell proportion
  modeling. Currently only 'MN' from
  [dist](https://rdrr.io/pkg/MGLM/man/dist.html) is supported.

- cp_intercept:

  a logical scalar for whether to include an intercept in the cell
  proportion model.

- ...:

  additional optional arguments.

## Value

outputs a list with following elements:

- `simu_cov`:

  a cell-by-covariate data frame of simulated cell types and
  corresponding individual.

- `cp_simu_df`:

  a cell type-by-covariate data frame summarizing the simulate cell
  proportions, total cell numbers, and cells per cell types.

- `cp_modelfit`:

  a fitted model object for the cell proportion model.

- `cn_modelfit`:

  a fitted model object for the cell number model.

## Examples

``` r
NULL
#> NULL
```
