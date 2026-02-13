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
  ncores = 2L,
  parallelization = c("pbmcapply", "future.apply", "parallel", "biocparallel"),
  BPPARAM = NULL,
  future.seed = FALSE,
  data_maxsize = 1
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

  a positive integer value (greater or equal to 1) to specify the number
  of CPU cores used in parallelization. The default is 2.

- parallelization:

  a string scalar specifying the parallelization backend used when
  simulating data. Must be one of "parallel", "future.apply",
  "biocparallel", or "pbmcapply". The default value is "pbmcapply". See
  details.

- BPPARAM:

  a BiocParallelParam class object (from `BiocParallel` R package) that
  must be specified when using `parallelization = "biocparallel"`.
  Either
  [`BiocParallel::SnowParam()`](https://rdrr.io/pkg/BiocParallel/man/SnowParam-class.html)
  or
  [`BiocParallel::MulticoreParam()`](https://rdrr.io/pkg/BiocParallel/man/MulticoreParam-class.html)
  can be used to initialize, depending on the operating system. BPPARAM
  is not used in other parallelization options. The default is NULL.

- future.seed:

  a logical or an integer (of length one or seven), or a list of
  length(X) with pre-generated random seeds that can be specified when
  using `parallelization = "future.apply"`. See
  [`future.apply::future_eapply`](https://future.apply.futureverse.org/reference/future_lapply.html)
  documentation for more details on its usage. future.seed is not used
  in other parallelization options. The default is FALSE.

- data_maxsize:

  a positive numeric value used to set max marginal_list size in GiB
  increments. Used only when `parallelization = "future.apply"`. The
  default is 1.

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
