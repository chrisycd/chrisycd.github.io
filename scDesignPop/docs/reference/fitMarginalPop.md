# Fit marginal models for every feature

Fits a specified parametric model using various input parameters.

## Usage

``` r
fitMarginalPop(
  data_list,
  mean_formula,
  model_family = "nb",
  interact_colnames = NULL,
  parallelization = c("pbmcapply", "future.apply", "parallel", "biocparallel"),
  n_cores = 2L,
  loc_colname = "POS",
  snp_colname = "snp_id",
  celltype_colname = "cell_type",
  indiv_colname = "indiv",
  filter_snps = TRUE,
  snpvar_thres = 0,
  force_formula = FALSE,
  keep_cellnames = FALSE,
  BPPARAM = NULL,
  future.seed = FALSE,
  data_maxsize = 1,
  ...
)
```

## Arguments

- data_list:

  a list of input data.

- mean_formula:

  a string scalar to specify the mean formula, including random effects
  (if any) but without SNP genotypes or SNP genotype interaction
  effects.

- model_family:

  a string scalar to specify model fitting used.

- interact_colnames:

  a string scalar or vector for the variable names that have first-order
  interaction with SNP genotypes.

- parallelization:

  a string scalar specifying the parallelization backend used during
  marginal fitting. Must be one of either "parallel", "future.apply",
  "biocparallel", or "pbmcapply". The default value is "pbmcapply". See
  details.

- n_cores:

  positive integer value (greater or equal to 1) to specify the number
  of CPU cores used in parallelization. The default is 2.

- loc_colname:

  a string scalar for column name of SNP position variable.

- snp_colname:

  a string scalar for column name of SNP id variable.

- celltype_colname:

  a string scalar for column name of cell type.

- indiv_colname:

  a string scalar for column name of individuals (samples).

- filter_snps:

  a logical scalar for whether to filter out SNP covariates with either
  low-variance or with only 1 distinct genotype (ie. all 1's) prior to
  fitting the model.

- snpvar_thres:

  a numeric scalar (between 0 to 1) used to filter out SNPs whose
  variance of genotypes across samples are below this threshold. Used
  together when `filter_snps = TRUE`.

- force_formula:

  a logical scalar for whether to bypass model parsimony check. If
  `force_formula = TRUE`, interaction terms whose covariates are not
  main effects in the model would be permitted. Results in error if
  `force_formula = FALSE` and `length(geno_interact_names) > 0`.

- keep_cellnames:

  a logical scalar for whether to keep cell barcode names. If
  `keep_cellnames = TRUE`, the memory will be larger. The default is
  FALSE.

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

  a positive numeric value used to set max data_list size in GiB
  increments. Used only when `parallelization = "future.apply"`. The
  default is 1.

- ...:

  additional arguments passed to internal functions.

## Value

a list of named features, each containing a list with the following
items:

- `fit`:

  a `glmmTMB` fit object.

- `time`:

  a numeric scalar for the elapsed time to fit the given feature.

- `snp_cov`:

  a string scalar or vector of SNP ids used in the fit for given
  feature.

- `model_attr`:

  a list of attributes extracted from each model. NOT currently
  implemented.

- `removed_cell`:

  a string scalar or vector of the cell names removed due to
  low-variance (NOT currently implemented).

## Details

### Parallelization options

If "parallel" is used then `mclapply` is called from the `parallel`
package; if "biocparallel" is used, then `bplapply` is called from the
`BiocParallel` package; if "future.apply" is used, then `future_lapply`
is called from the `future.apply` package; if "pbmcapply" is used, then
`pbmclapply` is called from the `pbmcapply` package.

## Examples

``` r
NULL
#> NULL
```
