# Modify parameters of a glmmTMB model object

Modify parameters of a glmmTMB model object

## Usage

``` r
modifyModelPara(
  model_obj,
  eqtlgeno,
  celltype,
  neg_ctrl = FALSE,
  mean_log2fc = 0,
  eqtl_log2fc = mean_log2fc,
  eqtl_reverse = FALSE,
  mean_baseline = NULL,
  eqtl_baseline = NULL,
  mean_baseline_only = FALSE,
  eqtl_baseline_only = FALSE,
  disp_scaling = "linear",
  cellstate_colname = "cell_type",
  snp_colname = "snp_id",
  verbose = TRUE,
  debug = FALSE,
  log_tol = 1e-04,
  ...
)
```

## Arguments

- model_obj:

  A marginal model object for a gene.

- eqtlgeno:

  A dataframe with eQTL annotations and samples' genotype for a gene.

- celltype:

  A string to specify the cell type in which to make the modification.

- neg_ctrl:

  A logical value for whether to set a negative control eQTL (ie. a
  non-eGene). This option sets the conditional means to be identical
  across genotypes (0, 1, 2). If `neg_ctrl = TRUE`, the `mean_log2fc`
  option will still be applied if set, but eqtl_log2fc will be overidden
  and have no impact. Default is `FALSE`.

- mean_log2fc:

  A numeric value for the log2 fold-change parameter to increase or
  decrease the conditional mean at genotype 1 \\\mu\_{1}\\ in a cell
  type. Default is `mean_log2fc = 0` (no parameters are modified and
  uses estimated parameters from the fitted marginal model).

- eqtl_log2fc:

  A numeric value for the log2 fold-change parameter to increase or
  decrease the slope of eQTL effect in a celltype. The eQTL slope is
  defined as the difference between the conditional mean at genotype 1
  and genotype 0 (\\\mu\_{1}\\ - \\\mu\_{0}\\). Default is
  `eqtl_log2fc = mean_log2fc` (eQTL slope is scaled the same as the
  conditional mean log2 fold-change).

- eqtl_reverse:

  A logical value to determine whether the eQTL slope trends in the
  reverse direction (TRUE) or same (FALSE). Default is `FALSE`.

- mean_baseline:

  A numeric value to specify the minimum conditional mean at genotype 1
  \\\mu\_{1}\\. If `mean_baseline_only = FALSE`, then the conditional
  mean will be the maximum of the fitted (estimated from marginal model)
  and the `mean_baseline` value. Otherwise, the conditional mean will be
  set to the `mean_baseline` value. Default value is `NULL`.

- eqtl_baseline:

  A numeric value to specify the minimum eQTL slope between genotype 1
  and 0 (\\\mu\_{1}\\ - \\\mu\_{0}\\). If `eqtl_baseline_only = FALSE`,
  then the eQTL slope will be the maximum of the slope of fitted
  (estimated from marginal model) and the `eqtl_baseline` value.
  Otherwise, the eQTL slope will be set to the `eqtl_baseline` value.
  Default value is `NULL`.

- mean_baseline_only:

  A logical value to force the conditional mean (in linear prediction)
  at genotype 1 \\\mu\_{1}\\. Default is `FALSE`.

- eqtl_baseline_only:

  A logical value to force the eQTL slope between genotype 1 and 0
  (\\\mu\_{1}\\ - \\\mu\_{0}\\). Default is `FALSE`.

- disp_scaling:

  A string value to specify the dispersion-mean scaling for certain
  parametric models. Current options are either `"linear"`,
  `"quadratic"`, or `"none"`. (NOTE: currently only applicable to the
  negative binomial model.)

- cellstate_colname:

  A string for cell state variable name.

- snp_colname:

  A string for SNP id variable name.

- verbose:

  A logical value for whether to output messages related to modified
  parameters. Default is `TRUE`.

- debug:

  A logical value for whether to output intermediate objects used for
  debugging purposes. Default is `FALSE`.

- log_tol:

  A numeric value used as tolerance in log computation. Default value is
  \\1e-4\\.

- ...:

  Additional options.

## Value

A list of dataframe of coefficients, model objects, and optional outputs
if debugging is enabled.

## Examples

``` r
NULL
#> NULL
```
