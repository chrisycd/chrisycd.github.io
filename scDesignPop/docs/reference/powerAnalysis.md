# Perform a cell-type-specific power analysis on eQTL effects

Perform a cell-type-specific power analysis on eQTL effects

## Usage

``` r
powerAnalysis(
  marginal_list,
  marginal_model = NULL,
  refit_formula = NULL,
  geneid = NULL,
  snpid = NULL,
  celltype_colname = "cell_type",
  celltype_vector = NULL,
  celltype_specific_ES_vector = NULL,
  indiv_colname = "indiv",
  method = c("nb", "poisson", "gaussian", "pseudoBulkLinear"),
  nindivs = NULL,
  ncells = NULL,
  nPool = NULL,
  nIndivPerPool = NULL,
  nCellPerPool = NULL,
  alpha = 0.05,
  nsims = 100,
  ncores = 1
)
```

## Arguments

- marginal_list:

  the output of function fitMarginalPop().

- marginal_model:

  a character showing the model types of the full marginal model.

- refit_formula:

  the formula used to refit the marginal full model if user wants to.
  Default is null.

- geneid:

  a character object contains geneid.

- snpid:

  a character object contains snpid.

- celltype_colname:

  a string scalar specifying the cell state variable in
  `marginal_list[[geneid]]$frame`. The default is "cell_type".

- celltype_vector:

  a vector object specifies the cell type that will be tested

- celltype_specific_ES_vector:

  a vector object specifies the genotype effect size (ES) for each cell
  type

- indiv_colname:

  a string scalar of the sample ID variable in cell covariate of
  `marginal_list[[geneid]]$frame`. The default is "indiv".

- method:

  a character object specifying the method that will be analyzed for
  power. (Options: nb,poisson,gaussian,pseudoBulkLinear).

- nindivs:

  a vector of numeric values showing the numbers of individuals that
  user wants to simulate.

- ncells:

  a vector of numeric values showing the numbers of cells per each
  individual that user wants to simulate.

- nPool:

  a vector of numeric values showing how many pools of sequencing has
  been performed.

- nIndivPerPool:

  a numerical value showing how many individuals are sequenced in one
  pool.

- nCellPerPool:

  a vector of numeric values showing how many cells are sequenced in one
  pool.

- alpha:

  the p value threshold for rejecting the H0 hypothesis.

- nsims:

  number of simulations for calculating the power. This parameter will
  affect the resolution of the power value.

- ncores:

  number of CPU cores user wants to use.

## Value

a list of named features, each containing a list with the following
items:

- `intercept`:

  the intercept for the genotype effect in the specified type/level.

- `slope`:

  the slope for the genotype effect in the specified type/level.

- `power`:

  a data frame contains power values in different parameter settings.

- `data`:

  a data frame contains both the H1 and H0 genotype effect estimates in
  different parameter settings and simulation times.

## Examples

``` r
NULL
#> NULL
```
