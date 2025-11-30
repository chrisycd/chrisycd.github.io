# Simulate new data

`simuNewPop` generates new simulated data based on fitted marginal and
copula models. This function is adapted from simu_new function in
scDesign3 v0.99.7

## Usage

``` r
simuNewPop(
  sce,
  assay_use = "counts",
  mean_mat,
  sigma_mat,
  zero_mat,
  quantile_mat = NULL,
  copula_list,
  n_cores,
  fastmvn = FALSE,
  family_use,
  nonnegative = TRUE,
  nonzerovar = FALSE,
  input_data,
  new_covariate,
  important_feature = "all",
  parallelization = "mcmapply",
  BPPARAM = NULL,
  filtered_gene,
  mean_limit = 1e+15,
  debug = FALSE
)
```

## Arguments

- sce:

  A `SingleCellExperiment` object.

- assay_use:

  A string which indicates the assay you will use in the sce. Default is
  'counts'.

- mean_mat:

  A cell by feature matrix of the mean parameter.

- sigma_mat:

  A cell by feature matrix of the sigma parameter.

- zero_mat:

  A cell by feature matrix of the zero-inflation parameter.

- quantile_mat:

  A cell by feature matrix of the multivariate quantile.

- copula_list:

  A list of copulas for generating the multivariate quantile matrix. If
  provided, the `quantile_mat` must be NULL.

- n_cores:

  An integer. The number of cores to use.

- fastmvn:

  An logical variable. If TRUE, the sampling of multivariate Gaussian is
  done by `mvnfast`, otherwise by `mvtnorm`. Default is FALSE.

- family_use:

  A string of the marginal distribution. Must be one of 'poisson', 'nb',
  or 'gaussian'.

- nonnegative:

  A logical variable. If TRUE, values \< 0 in the synthetic data will be
  converted to 0. Default is TRUE (since the expression matrix is
  nonnegative).

- nonzerovar:

  A logical variable. If TRUE, for any gene with zero variance, a cell
  will be replaced with 1. This is designed for avoiding potential
  errors, for example, PCA.

- input_data:

  A input count matrix.

- new_covariate:

  A data.frame which contains covariates of targeted simulated data from
  [`constructDataPop`](https://github.com/chrisycd/scDesignPop/reference/constructDataPop.md).

- important_feature:

  important_feature A string or vector which indicates whether a gene
  will be used in correlation estimation or not. If this is a string,
  then this string must be either "all" (using all genes) or "auto",
  which indicates that the genes will be automatically selected based on
  the proportion of zero expression across cells for each gene. Gene
  with zero proportion greater than 0.8 will be excluded form gene-gene
  correlation estimation. If this is a vector, then this should be a
  logical vector with length equal to the number of genes in `sce`.
  `TRUE` in the logical vector means the corresponding gene will be
  included in gene-gene correlation estimation and `FALSE` in the
  logical vector means the corresponding gene will be excluded from the
  gene-gene correlation estimation. The default value for is "all".

- parallelization:

  A string indicating the specific parallelization function to use. Must
  be one of 'mcmapply', 'bpmapply', or 'pbmcmapply', which corresponds
  to the parallelization function in the package
  `parallel`,`BiocParallel`, and `pbmcapply` respectively. The default
  value is 'mcmapply'.

- BPPARAM:

  A `MulticoreParam` object or NULL. When the parameter parallelization
  = 'mcmapply' or 'pbmcmapply', this parameter must be NULL. When the
  parameter parallelization = 'bpmapply', this parameter must be one of
  the `MulticoreParam` object offered by the package 'BiocParallel. The
  default value is NULL.

- filtered_gene:

  A vector or NULL which contains genes that are excluded in the
  marginal and copula fitting

- mean_limit:

  A numeric scalar to filter genes which has cells that exceed the limit
  in the `mean_mat`. The default value is 1e15. This is to avoid
  features that have extremely high and unreasonable means.

- debug:

  A logical scalar for whether to return a list of variables in addition
  to simulated count matrix. The default is FALSE.

## Value

A feature by cell matrix of the new simulated count (expression) matrix
or sparse matrix.

## Details

The function takes the new covariate (if use) from
[`constructDataPop`](https://github.com/chrisycd/scDesignPop/reference/constructDataPop.md),
parameter matrices from
[`extractParaPop`](https://github.com/chrisycd/scDesignPop/reference/extractParaPop.md)
and multivariate Unifs from
[`fitCopulaPop`](https://github.com/chrisycd/scDesignPop/reference/fitCopulaPop.md).

## Examples

``` r
NULL
#> NULL
```
