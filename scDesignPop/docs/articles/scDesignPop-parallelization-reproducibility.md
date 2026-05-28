# Parallel computing and reproducibility in scDesignPop

## Introduction

scDesignPop currently supports several different backend options for
parallel computing in its computationally intensive functions. This
tutorial is split into two parts: first, we explain the differences
between the parallelization options; second, we explain where parallel
computing is offered and where random number generation (RNG) occurs in
scDesignPop’s simulation workflow.

## Parallelization options

Currently, scDesignPop leverages the `parallel`, `pbmcapply`,
`BiocParallel`, and `future.apply` R packages for parallel computing.
These backends differ in their operating system support, reproducibility
guarantees, and usability.

### parallel

The `parallel` package offers parallelization using its apply-like
functions such as
[`mclapply()`](https://rdrr.io/r/parallel/mclapply.html) and `mcmapply`.
It uses forking for parallel processing, which is availabe on *Linux*
and *macOS* (M-series processors) systems, but not on *Windows* unless
`mc.cores = 1`.

By default, using [`set.seed()`](https://rdrr.io/r/base/Random.html)
**does not guarantee reproducibility** under parallel execution. This
can be seen in the example below, where two sets of $10$ random numbers
generated from a standard Gaussian distribution produce different
results despite using the same seed.

``` r
library(parallel)

set.seed(123)

unlist(mclapply(1:10, function(i) rnorm(1), mc.cores = 5))
#>  [1]  1.4528822 -0.1468824 -0.9429165 -1.0865478  0.1332380 -0.6302264
#>  [7]  1.8086562  0.8421334 -0.8377302  1.1807826
unlist(mclapply(1:10, function(i) rnorm(1), mc.cores = 5))
#>  [1] -0.7385392  0.5248758  0.9804686 -0.7200886  0.5263955 -1.2640428
#>  [7] -0.5579965  0.5095406 -1.7711345 -0.3073335
```

To ensure reproducibility, the RNG method must be set to
`"L'Ecuyer-CMRG"`, which supports parallel-safe random number streams.
Since this modifies the global RNG settings, it is good practice
afterwards to reset it to the default, which is `"Mersenne-Twister"` in
R.

``` r
set.seed(123, kind = "L'Ecuyer-CMRG")

unlist(mclapply(1:10, function(i) rnorm(1), mc.cores = 5))
#>  [1] -0.4094454 -0.4890608 -1.0388664  0.7613014 -1.1488680  0.8909694
#>  [7]  0.4330424  1.5745125  2.2994158  1.0644774
unlist(mclapply(1:10, function(i) rnorm(1), mc.cores = 5))
#>  [1] -0.4094454 -0.4890608 -1.0388664  0.7613014 -1.1488680  0.8909694
#>  [7]  0.4330424  1.5745125  2.2994158  1.0644774

RNGkind("Mersenne-Twister")
```

### pbmcapply

The `pbmcapply` package extends `parallel` by adding a **progress bar**
to the
[`pbmclapply()`](https://rdrr.io/pkg/pbmcapply/man/pbmclapply.html) and
[`pbmcmapply()`](https://rdrr.io/pkg/pbmcapply/man/pbmcmapply.html)
functions. This provides an estimated completion time, which can be
useful for long-running simulations. Like the `parallel` package, it
relies on forking, thus only supports parallel execution on *Linux* and
*macOS* (M-series processors) when `mc.cores > 1`.

Another key difference is that `pbmcapply` currently **does not generate
reproducible code** when RNG is present within the apply code. This is
demonstrated in the example below when using the `"L'Ecuyer-CMRG"` RNG
method. This behavior has been noted as an issue in various discussions,
noted
[here](https://stackoverflow.com/questions/67655726/parallel-processing-in-r-setting-seed-with-mclapply-vs-pbmclapply)
and [here](https://github.com/kvnkuang/pbmcapply/issues/47).

``` r
library(pbmcapply)

set.seed(123, kind = "L'Ecuyer-CMRG")

unlist(pbmclapply(1:10, function(i) rnorm(1), mc.cores = 5))
#>  [1] -0.4094454 -0.4890608 -1.0388664  0.7613014 -1.1488680  0.8909694
#>  [7]  0.4330424  1.5745125  2.2994158  1.0644774
unlist(pbmclapply(1:10, function(i) rnorm(1), mc.cores = 2))
#>  [1] -0.40944544 -0.48906078  0.89096942  0.43304237 -0.86537045 -0.03195349
#>  [7]  1.46427111  0.14670372  1.26748453 -1.75239095

RNGkind("Mersenne-Twister")
```

### BiocParallel

The `BiocParallel` package provides parallel computing across different
operating system platforms. In practice, the two relevant backends used
in scDesignPop are `MulticoreParam` and `SnowParam`. `MulticoreParam`
uses forking and is therefore available only on *Linux* and *macOS*,
whereas `SnowParam` uses separate R worker processes and is available on
*Windows*, *Linux*, and *macOS*. This makes `BiocParallel` a convenient
cross-operating system option.

`BiocParallel` has built-in support for reproducible random number
generation via its `RNGseed` argument instead of using
[`set.seed()`](https://rdrr.io/r/base/Random.html). Given the same seed,
identical results are returned across repeated runs, invariant to the
backend or number of workers used.

The toy example below illustrates this using `SnowParam` and
`MulticoreParam`, as well as using different number of workers.

``` r
library(BiocParallel)

unlist(bplapply(1:10, function(i) rnorm(1),
                BPPARAM = SnowParam(workers = 2, RNGseed = 123)))
#>  [1]  0.4254817 -2.0340375 -0.1219214  2.3644493  0.6723667  0.2474520
#>  [7]  0.2592861 -0.7682602  0.3982019 -0.8227634
unlist(bplapply(1:10, function(i) rnorm(1),
                BPPARAM = SnowParam(workers = 4, RNGseed = 123)))
#>  [1]  0.4254817 -2.0340375 -0.1219214  2.3644493  0.6723667  0.2474520
#>  [7]  0.2592861 -0.7682602  0.3982019 -0.8227634
unlist(bplapply(1:10, function(i) rnorm(1),
                BPPARAM = MulticoreParam(workers = 3, RNGseed = 123)))
#>  [1]  0.4254817 -2.0340375 -0.1219214  2.3644493  0.6723667  0.2474520
#>  [7]  0.2592861 -0.7682602  0.3982019 -0.8227634
unlist(bplapply(1:10, function(i) rnorm(1),
                BPPARAM = MulticoreParam(workers = 6, RNGseed = 123)))
#>  [1]  0.4254817 -2.0340375 -0.1219214  2.3644493  0.6723667  0.2474520
#>  [7]  0.2592861 -0.7682602  0.3982019 -0.8227634
```

Lastly, a progress bar can be enabled using the `progressbar = TRUE`
option.

### future.apply

The `future.apply` package offers parallel computing using the `future`
framework. Its two parallel computing options, `multisession` and
`multicore`, can be specified using
[`plan()`](https://future.futureverse.org/reference/plan.html).
`multisession` works across *Windows*, *Linux*, and *macOS* by launching
separate R sessions, whereas `multicore` uses forking and is therefore
limited to *Linux*, and *macOS* systems.

`future.apply` does support reproducible RNG, but a warning is produced
when RNG is detected when using its functions with default options.

``` r
library(future)
library(future.apply)

unlist(future_lapply(1:10, function(i) rnorm(1)))
#> Warning: UNRELIABLE VALUE: One of the 'future.apply' iterations
#> ('future_lapply-1') unexpectedly generated random numbers without declaring so.
#> There is a risk that those random numbers are not statistically sound and the
#> overall results might be invalid. To fix this, specify 'future.seed=TRUE'. This
#> ensures that proper, parallel-safe random numbers are produced via the
#> L'Ecuyer-CMRG method. To disable this check, use 'future.seed = NULL', or set
#> option 'future.rng.onMisuse' to "ignore".
#>  [1] -0.1207144 -0.7934806  0.2565556  0.4045028  1.9082073  1.6040098
#>  [7] -1.7511067 -0.5945040  0.4369179  0.5360741
```

To generate reproducible RNG, its `future.seed` argument must be
specified. This option ensures the same results are returned across
repeated runs given the same initial seed, regardless of the future
backend used or the number of workers. Since
[`plan()`](https://future.futureverse.org/reference/plan.html) changes
the global environment in R, it’s usually a good idea to reset this back
to the default.

The following example generates reproducible random numbers in parallel:

``` r
plan(multisession, workers = 5)

unlist(future_lapply(1:10, function(i) rnorm(1), future.seed = 123))
#>  [1] -1.01425025 -0.03074715  0.08281673 -0.96674384  0.30474400  0.37906458
#>  [7]  0.01244892  0.83755002  0.05058956  0.10957834
unlist(future_lapply(1:10, function(i) rnorm(1), future.seed = 123))
#>  [1] -1.01425025 -0.03074715  0.08281673 -0.96674384  0.30474400  0.37906458
#>  [7]  0.01244892  0.83755002  0.05058956  0.10957834

plan("default")
```

One caveat is that reproducibility in `future.apply` can introduce some
overhead, because seeds are generated in advance for all iterations.
This may increase computation time or memory use.

## Parallel computing in scDesignpop

Parallel computing is currently supported in scDesignPop’s
[`fitMarginalPop()`](https://chrisycd.github.io/scDesignPop/reference/fitMarginalPop.md),
[`fitCopulaPop()`](https://chrisycd.github.io/scDesignPop/reference/fitCopulaPop.md),
[`extractParaPop()`](https://chrisycd.github.io/scDesignPop/reference/extractParaPop.md),
[`simuNewPop()`](https://chrisycd.github.io/scDesignPop/reference/simuNewPop.md),
and
[`powerAnalysis()`](https://chrisycd.github.io/scDesignPop/reference/powerAnalysis.md)
functions. These are the most computationally intensive parts in
scDesignPop’s workflow. The four different parallel backend options are
supported as follows:

- For `parallel` and `pbmcparallel`, parallel computing is supported in
  scDesignPop when using `parallelization = "parallel"` or
  `parallelization = "pbmcapply"`.

- For `BiocParallel`, this backend option is enabled in scDesignPop when
  using `parallelization = "biocparallel"`. Note that this requires
  initializing either a class object using `BPPARAM = MulticoreParam()`
  or `BPPARAM = SnowParam()`.

- For `future.apply`, this backend option is enabled in scDesignPop when
  using `parallelization = "future.apply"`. Currently, only the
  `multisession` option is supported as it offers the best
  cross-operating system compatibility.

- The number of CPU cores used is controlled using the `n_cores`
  argument for all backend parallelization options.

As an example, we show how `parallel` and `BiocParallel` can be used and
where they impact data reproducibility in a typical scDesignPop
workflow.

### Loading data and constructing input

``` r
library(scDesignPop)
library(SingleCellExperiment)
library(SummarizedExperiment)

data(example_sce)
data(example_eqtlgeno)
```

### Parallelization in marginal fitting

We use a subset of $200$ genes from the `example_sce` data to do
marginal fitting. Here, we specify a negative binomial mixed model using
`model_family = "nb"`, with individuals as the random intercept and cell
type as the fixed effect `(1|indiv) + cell_type`. To specify
cell-type-specific eQTLs, we use the single-SNP mode via
`snp_mode = "single"` with genotype-cell-type interactions via
`interact_colnames = "cell_type"`.

scDesignPop’s marginal fitting process is largely deterministic (i.e.,
no RNG is present), so using `pbmcapply` still generates reproducible
marginal fittings. Below, we confirm each fitted gene for the $2$
repeated marginal fittings are the same.

``` r
gene_subset <- rownames(example_sce)[1:200]

data_list <- constructDataPop(
    sce = example_sce,
    eqtlgeno_df = example_eqtlgeno,
    new_covariate = as.data.frame(colData(example_sce)),
    copula_variable = "cell_type",
    slot_name = "counts",
    snp_mode = "single",
    overlap_features = gene_subset
    )
#> slot_name argument is deprecated and will be removed in a future update.Please use assay_use instead.
#> Constructing eqtlgeno list...

marginal_list <- fitMarginalPop(
    data_list = data_list,
    mean_formula = "(1|indiv) + cell_type",
    model_family = "nb",
    interact_colnames = "cell_type",
    parallelization = "pbmcapply",
    n_cores = 10L
)

marginal_list2 <- fitMarginalPop(
    data_list = data_list,
    mean_formula = "(1|indiv) + cell_type",
    model_family = "nb",
    interact_colnames = "cell_type",
    parallelization = "pbmcapply",
    n_cores = 10L
)

# reproducibility check
all(unlist(lapply(1:length(marginal_list), function(x) { 
    all.equal(marginal_list[[x]]$fit, marginal_list2[[x]]$fit)
    })))
#> Warning in all(unlist(lapply(1:length(marginal_list), function(x) {: coercing
#> argument of type 'character' to logical
#> [1] NA
```

Here, we use
[`all.equal()`](https://rdrr.io/pkg/Matrix/man/all.equal-methods.html)
to compare each fitted gene between the two marginal list objects
because there may be negligible differences
($< 2.220446 \times 10^{- 16}$) during the numerical optimization
process.

### Parallelization in copula fitting

Since we used a count model (NB mixed), a distributional transformation
is applied, which introduces RNG during the copula fitting process.
Thus, in order to have the same fitted copula, we can use the `parallel`
package. As discussed in the first part of this vignette, the RNG method
must be set to `"L'Ecuyer-CMRG"`. Using this, the example below shows
that the two copula fittings are identical.

``` r
set.seed(123, kind = "L'Ecuyer-CMRG")

copula_fit <- fitCopulaPop(
    sce = example_sce,
    assay_use = "counts",
    input_data = data_list$covariate,
    marginal_list = marginal_list,
    family_use = "nb",
    n_cores = 10L,
    parallelization = "parallel"
)
#> Convert Residuals to Multivariate Gaussian
#> Converting End
#> Copula group 1 starts
#> Copula group 3 starts
#> Copula group 4 starts
#> Copula group 2 starts

copula_fit2 <- fitCopulaPop(
    sce = example_sce,
    assay_use = "counts",
    input_data = data_list$covariate,
    marginal_list = marginal_list2,
    family_use = "nb",
    n_cores = 10L,
    parallelization = "parallel"
)
#> Convert Residuals to Multivariate Gaussian
#> Converting End
#> Copula group 1 starts
#> Copula group 3 starts
#> Copula group 4 starts
#> Copula group 2 starts

RNGkind("Mersenne-Twister")  # reset

# reproducibility check
identical(copula_fit, copula_fit2)
#> [1] FALSE
```

### Parallelization in other scDesignPop functions

Besides
[`fitCopulaPop()`](https://chrisycd.github.io/scDesignPop/reference/fitCopulaPop.md),
RNG occurs also in
[`extractParaPop()`](https://chrisycd.github.io/scDesignPop/reference/extractParaPop.md)
when simulating for new individuals not present in the training data, as
well as in
[`simuNewPop()`](https://chrisycd.github.io/scDesignPop/reference/simuNewPop.md)
and
[`powerAnalysis()`](https://chrisycd.github.io/scDesignPop/reference/powerAnalysis.md).
We will not cover usage of these in this vignette. In general, when
reproducibility is desired, appropriate seed options should be used
depending on the parallelization option.

In the future, we plan to offer user recommendations in terms of memory
usage and computational speed after we extensively test all combinations
of various parallel backend options.

## Session information

``` r
sessionInfo()
#> R version 4.2.3 (2023-03-15)
#> Platform: x86_64-pc-linux-gnu (64-bit)
#> Running under: Ubuntu 22.04.5 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> 
#> attached base packages:
#> [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
#> [8] methods   base     
#> 
#> other attached packages:
#>  [1] SingleCellExperiment_1.20.1 SummarizedExperiment_1.28.0
#>  [3] Biobase_2.58.0              GenomicRanges_1.50.2       
#>  [5] GenomeInfoDb_1.34.9         IRanges_2.32.0             
#>  [7] S4Vectors_0.36.2            BiocGenerics_0.44.0        
#>  [9] MatrixGenerics_1.10.0       matrixStats_1.1.0          
#> [11] scDesignPop_0.0.0.9012      future.apply_1.11.3        
#> [13] future_1.49.0               BiocParallel_1.32.6        
#> [15] pbmcapply_1.5.1             BiocStyle_2.26.0           
#> 
#> loaded via a namespace (and not attached):
#>  [1] sass_0.4.10            jsonlite_2.0.0         splines_4.2.3         
#>  [4] bslib_0.9.0            assertthat_0.2.1       BiocManager_1.30.25   
#>  [7] GenomeInfoDbData_1.2.9 yaml_2.3.10            globals_0.18.0        
#> [10] numDeriv_2016.8-1.1    pillar_1.10.2          lattice_0.22-6        
#> [13] glue_1.8.0             digest_0.6.37          XVector_0.38.0        
#> [16] RColorBrewer_1.1-3     glmmTMB_1.1.9          minqa_1.2.8           
#> [19] htmltools_0.5.8.1      Matrix_1.6-5           pkgconfig_2.0.3       
#> [22] listenv_0.9.1          zlibbioc_1.44.0        bookdown_0.43         
#> [25] mvtnorm_1.3-3          scales_1.4.0           snow_0.4-4            
#> [28] lme4_1.1-35.3          tibble_3.2.1           mgcv_1.9-1            
#> [31] generics_0.1.4         farver_2.1.2           ggplot2_3.5.2         
#> [34] withr_3.0.2            cachem_1.1.0           pbapply_1.7-2         
#> [37] TMB_1.9.11             cli_3.6.5              magrittr_2.0.3        
#> [40] evaluate_1.0.3         fs_1.6.6               parallelly_1.44.0     
#> [43] nlme_3.1-164           MASS_7.3-58.2          textshaping_0.4.0     
#> [46] tools_4.2.3            lifecycle_1.0.4        DelayedArray_0.24.0   
#> [49] irlba_2.3.5.1          compiler_4.2.3         pkgdown_2.2.0         
#> [52] jquerylib_0.1.4        systemfonts_1.2.3      rlang_1.1.6           
#> [55] grid_4.2.3             RCurl_1.98-1.17        nloptr_2.2.1          
#> [58] rstudioapi_0.17.1      htmlwidgets_1.6.4      bitops_1.0-9          
#> [61] rmarkdown_2.27         boot_1.3-30            gtable_0.3.6          
#> [64] codetools_0.2-20       R6_2.6.1               knitr_1.50            
#> [67] dplyr_1.1.4            fastmap_1.2.0          uwot_0.2.3            
#> [70] ragg_1.5.0             desc_1.4.3             Rcpp_1.0.14           
#> [73] vctrs_0.6.5            tidyselect_1.2.1       xfun_0.52
```
