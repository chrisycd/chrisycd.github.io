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
  n_cores = 2L,
  fastmvn = FALSE,
  family_use,
  nonnegative = TRUE,
  nonzerovar = FALSE,
  input_data,
  new_covariate,
  important_feature = "all",
  parallelization = c("pbmcapply", "future.apply", "parallel", "biocparallel"),
  BPPARAM = NULL,
  future.seed = FALSE,
  data_maxsize = 1,
  filtered_gene,
  mean_limit = 1e+15,
  debug = FALSE,
  ...
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

  a positive integer value (greater or equal to 1) to specify the number
  of CPU cores used in parallelization. The default is 2.

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

  a string scalar specifying the parallelization backend used when
  simulating data. Must be one of "parallel", "future.apply",
  "biocparallel", or "pbmcapply". The default value is "parallel". See
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

- ...:

  additional arguments passed to internal functions.

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

### Parallelization options

If "parallel" is used then `mcmapply` is called from the `parallel`
package; if "biocparallel" is used, then `bpmapply` is called from the
`BiocParallel` package; if "future.apply" is used, then `future_mapply`
is called from the `future.apply` package; if "pbmcapply" is used, then
`pbmcmapply` is called from the `pbmcapply` package.

## Examples

``` r
NULL
#> NULL
```
