# Extract parameter matrix for a new covariate data frame

This is the main function for extracting parameter matrices.

## Usage

``` r
extractParaPop(
  sce,
  assay_use = "counts",
  marginal_list,
  n_cores = 2L,
  family_use,
  new_covariate,
  new_eqtl_geno_list,
  indiv_colname = "indiv",
  snp_colname = "snp_id",
  loc_colname = "POS",
  parallelization = c("pbmcapply", "future.apply", "parallel", "biocparallel"),
  BPPARAM = NULL,
  future.seed = FALSE,
  data_maxsize = 1,
  data,
  ...
)
```

## Arguments

- sce:

  a SingleCellExperiment object.

- assay_use:

  a string scalar specifying the slot to use in input `sce`. The default
  is "counts".

- marginal_list:

  a list of named features, each with the fitted object and other
  variables as output from
  [`fitMarginalPop`](https://github.com/chrisycd/scDesignPop/reference/fitMarginalPop.md).

- n_cores:

  a positive integer value (greater or equal to 1) to specify the number
  of CPU cores used in parallelization. The default is 2.

- family_use:

  a string scalar or vector of marginal distribution used.

- new_covariate:

  a cell-by-covariate data frame obtained in the list output from
  [`constructDataPop`](https://github.com/chrisycd/scDesignPop/reference/constructDataPop.md).
  It must have a corr_group variable.

- new_eqtl_geno_list:

  a list of eQTL genotype data frames for each gene to be simulated. If
  using same list as in
  [`fitMarginalPop`](https://github.com/chrisycd/scDesignPop/reference/fitMarginalPop.md),
  then the in those samples

- indiv_colname:

  a string scalar of the sample ID variable in cell covariate of `sce`.
  The default is "indiv".

- snp_colname:

  a string scalar for the SNP variable in `eqtlgeno_df` used in
  [`constructDataPop`](https://github.com/chrisycd/scDesignPop/reference/constructDataPop.md).
  The default is "snp_id".

- loc_colname:

  a string scalar for the last column of eQTL annotation in
  `eqtlgeno_df`. The default is "POS".

- parallelization:

  a string scalar specifying the parallelization backend used when
  extracting parameters. Must be one of "parallel", "future.apply",
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

- data:

  a cell-by-covariate data frame obtained in the list output from
  [`constructDataPop`](https://github.com/chrisycd/scDesignPop/reference/constructDataPop.md).
  It must have a corr_group variable. Used only in gamlss fits.

- ...:

  additional arguments passed to internal functions.

## Value

a list of mean, sigma, and zero parameter cell by feature matrices:

- `mean_mat`:

  a cell by feature matrix containing the conditional mean values.

- `sigma_mat`:

  a cell by feature matrix containing the gene specific dispersion
  values.

- `zero_mat`:

  a cell by feature matrix containing the gene specific zero probability
  values (for zip and zinb models).

## Details

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
