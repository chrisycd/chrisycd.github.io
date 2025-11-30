# The wrapper function for power analysis

The wrapper function for power analysis

## Usage

``` r
runPowerAnalysis(
  marginal_list,
  marginal_model = "nb",
  refit_formula = NULL,
  geneid = NULL,
  snpid = NULL,
  celltype_colname = "cell_type",
  celltype_vector = NULL,
  celltype_specific_ES_list = NULL,
  indiv_colname = "indiv",
  methods = NULL,
  nindivs = NULL,
  ncells = NULL,
  nPool = NULL,
  nIndivPerPool = NULL,
  nCellPerPool = NULL,
  alpha = 0.05,
  power_nsim = 100,
  snp_number = 10,
  gene_number = 800,
  CI_nsim = 1000,
  CI_conf = 0.05,
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

- celltype_specific_ES_list:

  a list object specifies different vectors of the genotype effect size
  (ES) for each cell type

- indiv_colname:

  a string scalar of the sample ID variable in cell covariate of
  `marginal_list[[geneid]]$frame`. The default is "indiv".

- methods:

  a vector of character objects specifying the methods that will be
  analyzed for power. (Options: nb,poisson,gaussian,pseudoBulkLinear).

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

- power_nsim:

  a number of simulations for calculating the power. This parameter will
  affect the resolution of the power value.

- snp_number:

  the number of SNPs for multiple testing correction.

- gene_number:

  the number of genes for multiple testing correction.

- CI_nsim:

  number of simulations for calculating the Bootstrap CI.

- CI_conf:

  Bootstrap CI interval.

- ncores:

  number of CPU cores user wants to use.

## Value

a data frame contains power analysis result in different parameter
settings.

## Examples

``` r
NULL
#> NULL
```
