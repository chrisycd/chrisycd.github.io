# Modeling batch effect in population-scale scRNA-seq data

``` r
library(scDesignPop)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(ggplot2)
theme_set(theme_bw())
```

## Introduction

scDesignPop can also model population-scale scRNA-seq data with batch
information. The detailed tutorial is shown below.

## Fitting batch effects

### Step 1: construct a data list

To run scDesignPop, a list of data is required as input. This is done
using the `constructDataPop` function. A `SingleCellExperiment` object
and an `eqtlgeno` dataframe are the two main inputs needed. The
`eqtlgeno` dataframe consists of eQTL annotations (it must have cell
state, gene, SNP, chromosome, and position columns at a minimum), and
genotypes across individuals (columns) for every SNP (rows). The
structure of an example `eqtlgeno` dataframe is given below.

``` r
data("example_sce")
data("example_eqtlgeno")
```

Here, the `example_sce` also has a batch column.

``` r
colData(example_sce)
#> DataFrame with 7811 rows and 5 columns
#>            indiv cell_type       sex       age    batch
#>      <character>  <factor> <integer> <integer> <factor>
#> 1          SAMP1     cd4nc         2        65        3
#> 2          SAMP1     cd4nc         2        65        3
#> 3          SAMP1     cd4nc         2        65        3
#> 4          SAMP1     cd4nc         2        65        3
#> 5          SAMP1     cd4nc         2        65        3
#> ...          ...       ...       ...       ...      ...
#> 7807      SAMP14     monoc         1        87        1
#> 7808      SAMP14     bmem          1        87        1
#> 7809      SAMP14     bmem          1        87        1
#> 7810      SAMP14     bmem          1        87        1
#> 7811      SAMP14     bmem          1        87        1
```

``` r
data_list <- constructDataPop(
    sce = example_sce,
    eqtlgeno_df = example_eqtlgeno,
    new_covariate = as.data.frame(colData(example_sce)),
    overlap_features = NULL,
    sampid_vec = NULL,
    copula_variable = "cell_type",
    slot_name = "counts",
    snp_model = "single",
    celltype_colname = "cell_type",
    feature_colname = "gene_id",
    snp_colname = "snp_id",
    loc_colname = "POS",
    chrom_colname = "CHR",
    indiv_colname = "indiv",
    prune_thres = 0.9
    )
```

### Step 2: fit marginal model

Next, a marginal model is specified to fit each gene using the
`fitMarginalPop` function.  
Here we use a Negative Binominal as the parametric model using `"nb"`.
Since we want to include the batch effect modeling, we specify
`mean_formula = "(1|indiv) + cell_type + batch`.

``` r
marginal_list <- fitMarginalPop(
    data_list = data_list,
    mean_formula = "(1|indiv) + cell_type + batch",
    model_family = "nb",
    interact_colnames = "cell_type",
    parallelization = "pbmcapply",
    n_threads = 20L,
    loc_colname = "POS",
    snp_colname = "snp_id",
    celltype_colname = "cell_type",
    indiv_colname = "indiv",
    filter_snps = TRUE,
    snpvar_thres = 0,
    force_formula = FALSE,
    data_maxsize = 1
    )
```

### Step 3: fit a Gaussian copula

The third step is to fit a Gaussian copula using the `fitCopulaPop`
function.

``` r
set.seed(123, kind = "L'Ecuyer-CMRG")

copula_fit <- fitCopulaPop(
    sce = example_sce,
    assay_use = "counts",
    input_data = data_list[["new_covariate"]],
    marginal_list = marginal_list,
    family_use = "nb",
    copula = "gaussian",
    n_cores = 2L,
    parallelization = "mcmapply"
    )

RNGkind("Mersenne-Twister")  # reset
```

### Step 4: extract parameters

The fourth step is to compute the mean, sigma, and zero probability
parameters using the `extractParaPop` function.

``` r
para_new <- extractParaPop(
    sce = example_sce,
    assay_use = "counts",
    marginal_list = marginal_list,
    n_cores = 2L,
    family_use = "nb",
    indiv_colname = "indiv",
    new_covariate = data_list[["new_covariate"]],
    new_eqtl_geno_list = data_list[["eqtl_geno_list"]],
    data = data_list[["covariate"]],
    parallelization = "pbmcmapply"
    )
```

### Step 5: simulate counts

The fifth step is to simulate counts using the `simuNewPop` function.

``` r
set.seed(123)

newcount_mat <- simuNewPop(
    sce = example_sce,
    mean_mat = para_new[["mean_mat"]],
    sigma_mat = para_new[["sigma_mat"]],
    zero_mat = para_new[["zero_mat"]],
    quantile_mat = NULL,
    copula_list = copula_fit[["copula_list"]],
    n_cores = 2L,
    family_use = "nb",
    nonnegative = TRUE,
    input_data = data_list[["covariate"]],
    new_covariate = data_list[["new_covariate"]],
    important_feature = copula_fit[["important_feature"]],
    filtered_gene = data_list[["filtered_gene"]],
    parallelization = "pbmcmapply"
    )
```

## Removing the batch effects

### Step 1: set zero batch effect sizes in marginal models

Check the rank of the batch effect sizes in the marginal distributions.
Here, we take the gene `ENSG00000023902` as example since every gene are
fitted using the same formula.

``` r
marginal_list$ENSG00000023902$fit
#> Formula:          response ~ (1 | indiv) + cell_type + batch + `1:150133323` +  
#>     `1:150133323`:cell_type
#> Data: res_list[["dmat_df"]]
#>       AIC       BIC    logLik -2*log(L)  df.resid 
#>  8048.368  8138.891 -4011.184  8022.368      7798 
#> Random-effects (co)variances:
#> 
#> Conditional model:
#>  Groups Name        Std.Dev.
#>  indiv  (Intercept) 0.1807  
#> 
#> Number of obs: 7811 / Conditional model: indiv, 40
#> 
#> Dispersion parameter for nbinom2 family (): 1.63 
#> 
#> Fixed Effects:
#> 
#> Conditional model:
#>                   (Intercept)                cell_typemononc  
#>                      -0.87229                        0.33417  
#>                 cell_typebmem                 cell_typecd4nc  
#>                      -0.69781                       -1.65119  
#>                        batch2                         batch3  
#>                      -0.34148                       -0.12831  
#>                        batch4                  `1:150133323`  
#>                      -0.11636                        0.18264  
#> cell_typemononc:`1:150133323`    cell_typebmem:`1:150133323`  
#>                      -0.16217                       -0.04084  
#>  cell_typecd4nc:`1:150133323`  
#>                      -0.11806
```

After the ranks of the batch effects are obtained, set zero to the gene
marginal models and check the summary report again.

``` r
marginal_list_null <- lapply(marginal_list, function(x) {
  x$fit$fit$par[5:7] <- 0
  x$fit$fit$parfull[5:7] <- 0
  x$fit$sdr$par.fixed[5:7] <- 0
  x
})

marginal_list_null$ENSG00000023902$fit
#> Formula:          response ~ (1 | indiv) + cell_type + batch + `1:150133323` +  
#>     `1:150133323`:cell_type
#> Data: res_list[["dmat_df"]]
#>       AIC       BIC    logLik -2*log(L)  df.resid 
#>  8048.368  8138.891 -4011.184  8022.368      7798 
#> Random-effects (co)variances:
#> 
#> Conditional model:
#>  Groups Name        Std.Dev.
#>  indiv  (Intercept) 0.1807  
#> 
#> Number of obs: 7811 / Conditional model: indiv, 40
#> 
#> Dispersion parameter for nbinom2 family (): 1.63 
#> 
#> Fixed Effects:
#> 
#> Conditional model:
#>                   (Intercept)                cell_typemononc  
#>                      -0.87229                        0.33417  
#>                 cell_typebmem                 cell_typecd4nc  
#>                      -0.69781                       -1.65119  
#>                        batch2                         batch3  
#>                       0.00000                        0.00000  
#>                        batch4                  `1:150133323`  
#>                       0.00000                        0.18264  
#> cell_typemononc:`1:150133323`    cell_typebmem:`1:150133323`  
#>                      -0.16217                       -0.04084  
#>  cell_typecd4nc:`1:150133323`  
#>                      -0.11806
```

### Step 2: re-extract the parameters

``` r
para_new_null <- extractParaPop(
    sce = example_sce,
    assay_use = "counts",
    marginal_list = marginal_list_null,
    n_cores = 2L,
    family_use = "nb",
    indiv_colname = "indiv",
    new_covariate = data_list[["new_covariate"]],
    new_eqtl_geno_list = data_list[["eqtl_geno_list"]],
    data = data_list[["covariate"]],
    parallelization = "pbmcmapply"
    )
```

### Step 3: simulate counts with no batch effects

``` r
set.seed(123)

newcount_mat_null <- simuNewPop(
    sce = example_sce,
    mean_mat = para_new_null[["mean_mat"]],
    sigma_mat = para_new_null[["sigma_mat"]],
    zero_mat = para_new_null[["zero_mat"]],
    quantile_mat = NULL,
    copula_list = copula_fit[["copula_list"]],
    n_cores = 2L,
    family_use = "nb",
    nonnegative = TRUE,
    input_data = data_list[["covariate"]],
    new_covariate = data_list[["new_covariate"]],
    important_feature = copula_fit[["important_feature"]],
    filtered_gene = data_list[["filtered_gene"]],
    parallelization = "pbmcmapply"
    )
```

## Set new batch effects

### Step 1: set new batch effect sizes in marginal models

Check the rank of the batch effect sizes in the marginal distributions.
Here, we take the gene `ENSG00000023902` as example since every gene are
fitted using the same formula.

``` r
marginal_list$ENSG00000023902$fit
#> Formula:          response ~ (1 | indiv) + cell_type + batch + `1:150133323` +  
#>     `1:150133323`:cell_type
#> Data: res_list[["dmat_df"]]
#>       AIC       BIC    logLik -2*log(L)  df.resid 
#>  8048.368  8138.891 -4011.184  8022.368      7798 
#> Random-effects (co)variances:
#> 
#> Conditional model:
#>  Groups Name        Std.Dev.
#>  indiv  (Intercept) 0.1807  
#> 
#> Number of obs: 7811 / Conditional model: indiv, 40
#> 
#> Dispersion parameter for nbinom2 family (): 1.63 
#> 
#> Fixed Effects:
#> 
#> Conditional model:
#>                   (Intercept)                cell_typemononc  
#>                      -0.87229                        0.33417  
#>                 cell_typebmem                 cell_typecd4nc  
#>                      -0.69781                       -1.65119  
#>                        batch2                         batch3  
#>                      -0.34148                       -0.12831  
#>                        batch4                  `1:150133323`  
#>                      -0.11636                        0.18264  
#> cell_typemononc:`1:150133323`    cell_typebmem:`1:150133323`  
#>                      -0.16217                       -0.04084  
#>  cell_typecd4nc:`1:150133323`  
#>                      -0.11806
```

After the ranks of the batch effects are obtained, set a new value to
the gene marginal models and check the summary report again.

``` r
marginal_list_diff <- lapply(marginal_list, function(x) {
  val <- rnorm(3, mean = 5, sd = 2)
  x$fit$fit$par[5:7] <- val
  x$fit$fit$parfull[5:7] <- val
  x$fit$sdr$par.fixed[5:7] <- val
  x
})

marginal_list_diff$ENSG00000023902$fit
#> Formula:          response ~ (1 | indiv) + cell_type + batch + `1:150133323` +  
#>     `1:150133323`:cell_type
#> Data: res_list[["dmat_df"]]
#>       AIC       BIC    logLik -2*log(L)  df.resid 
#>  8048.368  8138.891 -4011.184  8022.368      7798 
#> Random-effects (co)variances:
#> 
#> Conditional model:
#>  Groups Name        Std.Dev.
#>  indiv  (Intercept) 0.1807  
#> 
#> Number of obs: 7811 / Conditional model: indiv, 40
#> 
#> Dispersion parameter for nbinom2 family (): 1.63 
#> 
#> Fixed Effects:
#> 
#> Conditional model:
#>                   (Intercept)                cell_typemononc  
#>                      -0.87229                        0.33417  
#>                 cell_typebmem                 cell_typecd4nc  
#>                      -0.69781                       -1.65119  
#>                        batch2                         batch3  
#>                       1.00619                        5.68744  
#>                        batch4                  `1:150133323`  
#>                       3.15227                        0.18264  
#> cell_typemononc:`1:150133323`    cell_typebmem:`1:150133323`  
#>                      -0.16217                       -0.04084  
#>  cell_typecd4nc:`1:150133323`  
#>                      -0.11806
```

### Step 2: re-extract the parameters

``` r
para_new_diff <- extractParaPop(
    sce = example_sce,
    assay_use = "counts",
    marginal_list = marginal_list_diff,
    n_cores = 2L,
    family_use = "nb",
    indiv_colname = "indiv",
    new_covariate = data_list[["new_covariate"]],
    new_eqtl_geno_list = data_list[["eqtl_geno_list"]],
    data = data_list[["covariate"]],
    parallelization = "pbmcmapply"
    )
```

### Step 3: simulate counts with no batch effects

``` r
set.seed(123)

newcount_mat_diff <- simuNewPop(
    sce = example_sce,
    mean_mat = para_new_diff[["mean_mat"]],
    sigma_mat = para_new_diff[["sigma_mat"]],
    zero_mat = para_new_diff[["zero_mat"]],
    quantile_mat = NULL,
    copula_list = copula_fit[["copula_list"]],
    n_cores = 2L,
    family_use = "nb",
    nonnegative = TRUE,
    input_data = data_list[["covariate"]],
    new_covariate = data_list[["new_covariate"]],
    important_feature = copula_fit[["important_feature"]],
    filtered_gene = data_list[["filtered_gene"]],
    parallelization = "pbmcmapply"
    )
```

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
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] ggplot2_3.5.2               SingleCellExperiment_1.20.1
#>  [3] SummarizedExperiment_1.28.0 Biobase_2.58.0             
#>  [5] GenomicRanges_1.50.2        GenomeInfoDb_1.34.9        
#>  [7] IRanges_2.32.0              S4Vectors_0.36.2           
#>  [9] BiocGenerics_0.44.0         MatrixGenerics_1.10.0      
#> [11] matrixStats_1.1.0           scDesignPop_0.0.0.9009     
#> [13] BiocStyle_2.26.0           
#> 
#> loaded via a namespace (and not attached):
#>  [1] sass_0.4.10            jsonlite_2.0.0         splines_4.2.3         
#>  [4] bslib_0.9.0            Rdpack_2.6.4           assertthat_0.2.1      
#>  [7] BiocManager_1.30.25    GenomeInfoDbData_1.2.9 yaml_2.3.10           
#> [10] numDeriv_2016.8-1.1    pillar_1.10.2          lattice_0.22-6        
#> [13] glue_1.8.0             reformulas_0.4.1       digest_0.6.37         
#> [16] RColorBrewer_1.1-3     XVector_0.38.0         rbibutils_2.3         
#> [19] glmmTMB_1.1.13         minqa_1.2.8            sandwich_3.1-1        
#> [22] htmltools_0.5.8.1      Matrix_1.6-5           pkgconfig_2.0.3       
#> [25] bookdown_0.43          zlibbioc_1.44.0        mvtnorm_1.3-3         
#> [28] scales_1.4.0           lme4_1.1-35.3          tibble_3.2.1          
#> [31] mgcv_1.9-1             generics_0.1.4         farver_2.1.2          
#> [34] cachem_1.1.0           withr_3.0.2            pbapply_1.7-2         
#> [37] TMB_1.9.11             cli_3.6.5              magrittr_2.0.3        
#> [40] evaluate_1.0.3         fs_1.6.6               nlme_3.1-164          
#> [43] MASS_7.3-58.2          textshaping_0.4.0      tools_4.2.3           
#> [46] lifecycle_1.0.4        DelayedArray_0.24.0    compiler_4.2.3        
#> [49] pkgdown_2.2.0          jquerylib_0.1.4        pbmcapply_1.5.1       
#> [52] systemfonts_1.2.3      rlang_1.1.6            grid_4.2.3            
#> [55] RCurl_1.98-1.17        nloptr_2.2.1           rstudioapi_0.17.1     
#> [58] htmlwidgets_1.6.4      bitops_1.0-9           rmarkdown_2.27        
#> [61] boot_1.3-30            gtable_0.3.6           R6_2.6.1              
#> [64] zoo_1.8-14             knitr_1.50             dplyr_1.1.4           
#> [67] fastmap_1.2.0          ragg_1.5.0             desc_1.4.3            
#> [70] parallel_4.2.3         Rcpp_1.0.14            vctrs_0.6.5           
#> [73] tidyselect_1.2.1       xfun_0.52
```
