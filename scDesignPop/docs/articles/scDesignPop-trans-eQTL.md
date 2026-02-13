# Modeling trans-eQTL effect using public data

``` r
library(scDesignPop)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(ggplot2)
theme_set(theme_bw())
```

## Introduction

scDesignPop supports user-defined eQTL pairs, making it straightforward
to extend the framework to trans-eQTL studies with known ground truth.

## Library and data preparation

### prepare the eQTL genotype dataframe as input

Here, we take the trans-eQTL list from eQTLGen phase 1
(<https://www.eqtlgen.org/trans-eqtls.html>) as an example to show the
scDesignPop’s flexibility in eQTL ground truth settings. The download
link of the trans-eQTL list is
<https://download.gcc.rug.nl/downloads/eqtlgen/trans-eqtl/2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded.txt.gz>.

``` bash
wget https://download.gcc.rug.nl/downloads/eqtlgen/trans-eqtl/2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded.txt.gz
#> --2026-02-13 01:24:58--  https://download.gcc.rug.nl/downloads/eqtlgen/trans-eqtl/2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded.txt.gz
#> Resolving download.gcc.rug.nl (download.gcc.rug.nl)... 195.169.22.66
#> Connecting to download.gcc.rug.nl (download.gcc.rug.nl)|195.169.22.66|:443... connected.
#> HTTP request sent, awaiting response... 200 OK
#> Length: 1839918 (1.8M) [application/octet-stream]
#> Saving to: ‘2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded.txt.gz.12’
#> 
#>      0K .......... .......... .......... .......... ..........  2%  177K 10s
#>     50K .......... .......... .......... .......... ..........  5%  329K 7s
#>    100K .......... .......... .......... .......... ..........  8% 61.5M 5s
#>    150K .......... .......... .......... .......... .......... 11% 2.37M 4s
#>    200K .......... .......... .......... .......... .......... 13%  382K 4s
#>    250K .......... .......... .......... .......... .......... 16% 68.6M 3s
#>    300K .......... .......... .......... .......... .......... 19% 85.3M 2s
#>    350K .......... .......... .......... .......... .......... 22%  107M 2s
#>    400K .......... .......... .......... .......... .......... 25% 2.48M 2s
#>    450K .......... .......... .......... .......... .......... 27%  382K 2s
#>    500K .......... .......... .......... .......... .......... 30% 84.8M 2s
#>    550K .......... .......... .......... .......... .......... 33%  141M 1s
#>    600K .......... .......... .......... .......... .......... 36% 84.5M 1s
#>    650K .......... .......... .......... .......... .......... 38%  103M 1s
#>    700K .......... .......... .......... .......... .......... 41%  112M 1s
#>    750K .......... .......... .......... .......... .......... 44%  129M 1s
#>    800K .......... .......... .......... .......... .......... 47% 2.69M 1s
#>    850K .......... .......... .......... .......... .......... 50%  122M 1s
#>    900K .......... .......... .......... .......... .......... 52%  383K 1s
#>    950K .......... .......... .......... .......... .......... 55%  123M 1s
#>   1000K .......... .......... .......... .......... .......... 58% 89.0M 1s
#>   1050K .......... .......... .......... .......... .......... 61%  136M 1s
#>   1100K .......... .......... .......... .......... .......... 64%  131M 1s
#>   1150K .......... .......... .......... .......... .......... 66% 94.2M 0s
#>   1200K .......... .......... .......... .......... .......... 69%  119M 0s
#>   1250K .......... .......... .......... .......... .......... 72%  136M 0s
#>   1300K .......... .......... .......... .......... .......... 75%  130M 0s
#>   1350K .......... .......... .......... .......... .......... 77% 80.2M 0s
#>   1400K .......... .......... .......... .......... .......... 80%  143K 0s
#>   1450K .......... .......... .......... .......... .......... 83%  110K 0s
#>   1500K .......... .......... .......... .......... .......... 86%  329K 0s
#>   1550K .......... .......... .......... .......... .......... 89%  329K 0s
#>   1600K .......... .......... .......... .......... .......... 91%  329K 0s
#>   1650K .......... .......... .......... .......... .......... 94%  329K 0s
#>   1700K .......... .......... .......... .......... .......... 97% 93.6M 0s
#>   1750K .......... .......... .......... .......... ......    100%  170M=2.3s
#> 
#> 2026-02-13 01:25:01 (778 KB/s) - ‘2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded.txt.gz.12’ saved [1839918/1839918]
```

Then, we reformat this trans-eQTL list into the required eQTL dataframe
as follows,

``` r
library(scDesignPop)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(ggplot2)
library(data.table)

transeqtl <- fread("2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded.txt.gz")
transeqtl <- transeqtl[,c("Gene","SNP","SNPChr","SNPPos")]
transeqtl <- cbind("bin", transeqtl)
colnames(transeqtl) <- c("cell_type","gene_id","snp_id","CHR","POS")
transeqtl <- transeqtl[!duplicated(transeqtl$gene_id), ]
transeqtl$snp_id <- paste(transeqtl$CHR,transeqtl$POS,sep = ":")
```

Since the genotype data are not provided and are often under restricted
access, here we extract the SNPs that we have genotype data in the
`example_eqtlgeno`

``` r
data("example_eqtlgeno")
transeqtl <- transeqtl[which(transeqtl$snp_id%in%example_eqtlgeno$snp_id),]

geno_cols <- grep("^SAMP", names(example_eqtlgeno), value = TRUE)
geno_by_snp <- example_eqtlgeno[!duplicated(example_eqtlgeno$snp_id), 
                                c("snp_id", geno_cols)]
idx <- match(transeqtl$snp_id, geno_by_snp$snp_id)
trans_eqtlgeno <- cbind(
  transeqtl,
  geno_by_snp[idx, geno_cols, drop = FALSE]
)
trans_eqtlgeno
#>      cell_type         gene_id      snp_id   CHR       POS SAMP1 SAMP2 SAMP3
#>         <char>          <char>      <char> <int>     <int> <num> <num> <num>
#>   1:       bin ENSG00000172262 12:69732105    12  69732105     2     0     1
#>   2:       bin ENSG00000258730 12:69732105    12  69732105     2     0     1
#>   3:       bin ENSG00000117862   11:327143    11    327143     0     1     0
#>   4:       bin ENSG00000171101 10:94450233    10  94450233     1     2     1
#>   5:       bin ENSG00000197057 10:94450233    10  94450233     1     2     1
#>  ---                                                                        
#> 129:       bin ENSG00000138433 3:127474030     3 127474030     0     2     2
#> 130:       bin ENSG00000182405  5:64355060     5  64355060     2     2     2
#> 131:       bin ENSG00000105085 12:56401085    12  56401085     1     0     2
#> 132:       bin ENSG00000214955 3:194403578     3 194403578     0     1     1
#> 133:       bin ENSG00000120093 10:94450233    10  94450233     1     2     1
#>      SAMP4 SAMP5 SAMP6 SAMP7 SAMP8 SAMP9 SAMP10 SAMP11 SAMP12 SAMP13 SAMP14
#>      <num> <num> <num> <num> <num> <num>  <num>  <num>  <num>  <num>  <num>
#>   1:     2     1     1     1     1     0      2      2      1      2      1
#>   2:     2     1     1     1     1     0      2      2      1      2      1
#>   3:     1     0     1     0     2     1      0      1      1      1      0
#>   4:     0     2     2     1     0     0      1      1      0      2      2
#>   5:     0     2     2     1     0     0      1      1      0      2      2
#>  ---                                                                       
#> 129:     2     2     0     1     2     0      2      2      2      2      2
#> 130:     0     2     2     1     1     2      1      0      0      0      1
#> 131:     2     1     0     2     1     1      1      1      0      1      1
#> 132:     2     2     0     1     0     1      0      1      1      0      1
#> 133:     0     2     2     1     0     0      1      1      0      2      2
#>      SAMP15 SAMP16 SAMP17 SAMP18 SAMP19 SAMP20 SAMP21 SAMP22 SAMP23 SAMP24
#>       <num>  <num>  <num>  <num>  <num>  <num>  <num>  <num>  <num>  <num>
#>   1:      1      1      0      1      0      1      2      1      1      1
#>   2:      1      1      0      1      0      1      2      1      1      1
#>   3:      0      0      1      0      0      1      0      0      0      0
#>   4:      1      0      1      0      0      2      1      2      1      0
#>   5:      1      0      1      0      0      2      1      2      1      0
#>  ---                                                                      
#> 129:      2      1      1      1      0      2      1      0      2      1
#> 130:      0      2      2      2      1      1      1      1      2      1
#> 131:      1      1      2      0      1      1      2      2      0      0
#> 132:      2      1      1      2      1      0      2      1      2      1
#> 133:      1      0      1      0      0      2      1      2      1      0
#>      SAMP25 SAMP26 SAMP27 SAMP28 SAMP29 SAMP30 SAMP31 SAMP32 SAMP33 SAMP34
#>       <num>  <num>  <num>  <num>  <num>  <num>  <num>  <num>  <num>  <num>
#>   1:      0      0      1      1      1      2      1      1      0      0
#>   2:      0      0      1      1      1      2      1      1      0      0
#>   3:      0      1      1      1      0      0      0      2      1      0
#>   4:      0      1      1      0      0      1      2      0      1      1
#>   5:      0      1      1      0      0      1      2      0      1      1
#>  ---                                                                      
#> 129:      2      1      2      0      1      1      2      1      2      1
#> 130:      2      0      2      1      0      0      2      2      0      1
#> 131:      1      2      2      0      2      0      2      2      1      2
#> 132:      1      2      2      1      1      0      0      0      0      0
#> 133:      0      1      1      0      0      1      2      0      1      1
#>      SAMP35 SAMP36 SAMP37 SAMP38 SAMP39 SAMP40
#>       <num>  <num>  <num>  <num>  <num>  <num>
#>   1:      1      1      1      1      2      1
#>   2:      1      1      1      1      2      1
#>   3:      0      0      0      0      0      0
#>   4:      1      1      0      0      2      2
#>   5:      1      1      0      0      2      2
#>  ---                                          
#> 129:      2      2      2      2      1      2
#> 130:      1      2      2      2      2      1
#> 131:      2      1      2      1      2      1
#> 132:      1      1      2      1      0      1
#> 133:      1      1      0      0      2      2
```

## Modeling and simulation

### Step 1: construct a data list

To run scDesignPop, a list of data is required as input. This is done
using the `constructDataPop` function. A `SingleCellExperiment` object
and an `eqtlgeno` dataframe are the two main inputs needed. The
`eqtlgeno` dataframe consists of eQTL annotations (it must have cell
state, gene, SNP, chromosome, and position columns at a minimum), and
genotypes across individuals (columns) for every SNP (rows). Here, we
use `example_sce` to filter for overlapped genes in \`trans_eqtlgeno’
object created

``` r
data("example_sce")
overlap_genes <- intersect(rownames(example_sce),trans_eqtlgeno$gene_id)
trans_eqtlgeno <- trans_eqtlgeno[which(trans_eqtlgeno$gene_id%in%overlap_genes),]
example_sce <- example_sce[which(rownames(example_sce)%in%overlap_genes),]
```

``` r
data_list <- constructDataPop(
    sce = example_sce,
    eqtlgeno_df = trans_eqtlgeno,
    new_covariate = as.data.frame(colData(example_sce)),
    overlap_features = NULL,
    sampid_vec = NULL,
    copula_variable = "cell_type",
    slot_name = "counts",
    snp_mode = "single",
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

``` r
marginal_list <- fitMarginalPop(
    data_list = data_list,
    mean_formula = "(1|indiv) + cell_type",
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

### Step 6: create SingleCellExperiment object using simulated data

After simulating the data, we can create a `SingleCellExperiment` object
as follows.

``` r
simu_sce <- SingleCellExperiment(list(counts = newcount_mat), 
                                 colData = data_list[["new_covariate"]])
names(assays(simu_sce)) <- "counts"

rowData(simu_sce) <- rowData(example_sce)
simu_sce
#> class: SingleCellExperiment 
#> dim: 22 7811 
#> metadata(0):
#> assays(1): counts
#> rownames(22): ENSG00000198574 ENSG00000189292 ... ENSG00000100079
#>   ENSG00000225783
#> rowData names(0):
#> colnames: NULL
#> colData names(6): indiv cell_type ... batch corr_group
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
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
#>  [1] data.table_1.17.4           ggplot2_3.5.2              
#>  [3] SingleCellExperiment_1.20.1 SummarizedExperiment_1.28.0
#>  [5] Biobase_2.58.0              GenomicRanges_1.50.2       
#>  [7] GenomeInfoDb_1.34.9         IRanges_2.32.0             
#>  [9] S4Vectors_0.36.2            BiocGenerics_0.44.0        
#> [11] MatrixGenerics_1.10.0       matrixStats_1.1.0          
#> [13] scDesignPop_0.0.0.9010      BiocStyle_2.26.0           
#> 
#> loaded via a namespace (and not attached):
#>  [1] viridis_0.6.5          sass_0.4.10            jsonlite_2.0.0        
#>  [4] viridisLite_0.4.2      splines_4.2.3          R.utils_2.13.0        
#>  [7] RhpcBLASctl_0.23-42    bslib_0.9.0            assertthat_0.2.1      
#> [10] BiocManager_1.30.25    GenomeInfoDbData_1.2.9 yaml_2.3.10           
#> [13] numDeriv_2016.8-1.1    pillar_1.10.2          lattice_0.22-6        
#> [16] glue_1.8.0             digest_0.6.37          RColorBrewer_1.1-3    
#> [19] XVector_0.38.0         glmmTMB_1.1.9          minqa_1.2.8           
#> [22] R.oo_1.27.1            htmltools_0.5.8.1      Matrix_1.6-5          
#> [25] pkgconfig_2.0.3        bookdown_0.43          zlibbioc_1.44.0       
#> [28] mvtnorm_1.3-3          scales_1.4.0           lme4_1.1-35.3         
#> [31] tibble_3.2.1           mgcv_1.9-1             generics_0.1.4        
#> [34] farver_2.1.2           withr_3.0.2            cachem_1.1.0          
#> [37] pbapply_1.7-2          TMB_1.9.11             cli_3.6.5             
#> [40] magrittr_2.0.3         evaluate_1.0.3         R.methodsS3_1.8.2     
#> [43] fs_1.6.6               nlme_3.1-164           MASS_7.3-58.2         
#> [46] textshaping_0.4.0      tools_4.2.3            lifecycle_1.0.4       
#> [49] DelayedArray_0.24.0    irlba_2.3.5.1          compiler_4.2.3        
#> [52] pkgdown_2.2.0          jquerylib_0.1.4        pbmcapply_1.5.1       
#> [55] systemfonts_1.2.3      rlang_1.1.6            grid_4.2.3            
#> [58] RCurl_1.98-1.17        nloptr_2.2.1           rstudioapi_0.17.1     
#> [61] htmlwidgets_1.6.4      bitops_1.0-9           rmarkdown_2.27        
#> [64] boot_1.3-30            gtable_0.3.6           R6_2.6.1              
#> [67] gridExtra_2.3          knitr_1.50             dplyr_1.1.4           
#> [70] fastmap_1.2.0          uwot_0.2.3             ragg_1.5.0            
#> [73] desc_1.4.3             parallel_4.2.3         Rcpp_1.0.14           
#> [76] vctrs_0.6.5            tidyselect_1.2.1       xfun_0.52
```
