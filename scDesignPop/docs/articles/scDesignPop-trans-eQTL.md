# Modeling trans-eQTL effect using public data

``` r
library(scDesignPop)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(ggplot2)
theme_set(theme_bw())
```

## Introduction

scDesignPop allows users to specify eQTL pairs, thus providing the
capability to also provide ground truth for trans eQTL studies.

## Step 0: prepare the eQTL genotype dataframe as input

Here, we take the trans-eQTL list from eQTLGen phase 1
(<https://www.eqtlgen.org/trans-eqtls.html>) as an example to show the
scDesignPop’s flexibility in eQTL ground truth settings. The download
link of the trans-eQTL list is
<https://download.gcc.rug.nl/downloads/eqtlgen/trans-eqtl/2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded.txt.gz>.

``` bash
wget https://download.gcc.rug.nl/downloads/eqtlgen/trans-eqtl/2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded.txt.gz
#> --2026-01-23 17:58:42--  https://download.gcc.rug.nl/downloads/eqtlgen/trans-eqtl/2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded.txt.gz
#> Resolving download.gcc.rug.nl (download.gcc.rug.nl)... 195.169.22.66
#> Connecting to download.gcc.rug.nl (download.gcc.rug.nl)|195.169.22.66|:443... connected.
#> HTTP request sent, awaiting response... 200 OK
#> Length: 1839918 (1.8M) [application/octet-stream]
#> Saving to: ‘2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded.txt.gz.7’
#> 
#>      0K .......... .......... .......... .......... ..........  2%  165K 11s
#>     50K .......... .......... .......... .......... ..........  5%  330K 8s
#>    100K .......... .......... .......... .......... ..........  8% 74.6M 5s
#>    150K .......... .......... .......... .......... .......... 11% 1.20M 4s
#>    200K .......... .......... .......... .......... .......... 13%  449K 4s
#>    250K .......... .......... .......... .......... .......... 16% 91.4M 3s
#>    300K .......... .......... .......... .......... .......... 19% 81.3M 3s
#>    350K .......... .......... .......... .......... .......... 22% 77.6M 2s
#>    400K .......... .......... .......... .......... .......... 25% 1.24M 2s
#>    450K .......... .......... .......... .......... .......... 27%  451K 2s
#>    500K .......... .......... .......... .......... .......... 30% 70.0M 2s
#>    550K .......... .......... .......... .......... .......... 33%  106M 2s
#>    600K .......... .......... .......... .......... .......... 36%  109M 1s
#>    650K .......... .......... .......... .......... .......... 38% 94.3M 1s
#>    700K .......... .......... .......... .......... .......... 41%  126M 1s
#>    750K .......... .......... .......... .......... .......... 44%  130M 1s
#>    800K .......... .......... .......... .......... .......... 47% 1.29M 1s
#>    850K .......... .......... .......... .......... .......... 50% 99.2M 1s
#>    900K .......... .......... .......... .......... .......... 52%  454K 1s
#>    950K .......... .......... .......... .......... .......... 55% 72.0M 1s
#>   1000K .......... .......... .......... .......... .......... 58%  102M 1s
#>   1050K .......... .......... .......... .......... .......... 61% 76.8M 1s
#>   1100K .......... .......... .......... .......... .......... 64%  151M 1s
#>   1150K .......... .......... .......... .......... .......... 66% 99.1M 0s
#>   1200K .......... .......... .......... .......... .......... 69% 57.0M 0s
#>   1250K .......... .......... .......... .......... .......... 72%  171M 0s
#>   1300K .......... .......... .......... .......... .......... 75%  196M 0s
#>   1350K .......... .......... .......... .......... .......... 77%  176M 0s
#>   1400K .......... .......... .......... .......... .......... 80% 55.5M 0s
#>   1450K .......... .......... .......... .......... .......... 83%  327K 0s
#>   1500K .......... .......... .......... .......... .......... 86% 82.5M 0s
#>   1550K .......... .......... .......... .......... .......... 89%  255M 0s
#>   1600K .......... .......... .......... .......... .......... 91% 1.40M 0s
#>   1650K .......... .......... .......... .......... .......... 94%  131M 0s
#>   1700K .......... .......... .......... .......... .......... 97%  117M 0s
#>   1750K .......... .......... .......... .......... ......    100%  277M=1.1s
#> 
#> 2026-01-23 17:58:44 (1.59 MB/s) - ‘2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded.txt.gz.7’ saved [1839918/1839918]
```

Then, we reformat this trans-eQTL list into the required eQTL dataframe
as follows,

``` r
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

## Step 1: construct a data list

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

## Step 2: fit marginal model

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

## Step 3: fit a Gaussian copula

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

## Step 4: extract parameters

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

## Step 5: simulate counts

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

## Step 6: create SingleCellExperiment object using simulated data

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
#> [13] scDesignPop_0.0.0.9009      BiocStyle_2.26.0           
#> 
#> loaded via a namespace (and not attached):
#>  [1] sass_0.4.10            jsonlite_2.0.0         splines_4.2.3         
#>  [4] R.utils_2.13.0         bslib_0.9.0            Rdpack_2.6.4          
#>  [7] assertthat_0.2.1       BiocManager_1.30.25    GenomeInfoDbData_1.2.9
#> [10] yaml_2.3.10            numDeriv_2016.8-1.1    pillar_1.10.2         
#> [13] lattice_0.22-6         glue_1.8.0             reformulas_0.4.1      
#> [16] digest_0.6.37          rbibutils_2.3          RColorBrewer_1.1-3    
#> [19] XVector_0.38.0         glmmTMB_1.1.13         minqa_1.2.8           
#> [22] sandwich_3.1-1         htmltools_0.5.8.1      Matrix_1.6-5          
#> [25] R.oo_1.27.1            pkgconfig_2.0.3        bookdown_0.43         
#> [28] zlibbioc_1.44.0        mvtnorm_1.3-3          scales_1.4.0          
#> [31] lme4_1.1-35.3          tibble_3.2.1           mgcv_1.9-1            
#> [34] generics_0.1.4         farver_2.1.2           cachem_1.1.0          
#> [37] withr_3.0.2            pbapply_1.7-2          TMB_1.9.11            
#> [40] cli_3.6.5              magrittr_2.0.3         evaluate_1.0.3        
#> [43] R.methodsS3_1.8.2      fs_1.6.6               nlme_3.1-164          
#> [46] MASS_7.3-58.2          textshaping_0.4.0      tools_4.2.3           
#> [49] lifecycle_1.0.4        DelayedArray_0.24.0    compiler_4.2.3        
#> [52] pkgdown_2.2.0          jquerylib_0.1.4        pbmcapply_1.5.1       
#> [55] systemfonts_1.2.3      rlang_1.1.6            grid_4.2.3            
#> [58] RCurl_1.98-1.17        nloptr_2.2.1           rstudioapi_0.17.1     
#> [61] htmlwidgets_1.6.4      bitops_1.0-9           rmarkdown_2.27        
#> [64] boot_1.3-30            gtable_0.3.6           R6_2.6.1              
#> [67] zoo_1.8-14             knitr_1.50             dplyr_1.1.4           
#> [70] fastmap_1.2.0          ragg_1.5.0             desc_1.4.3            
#> [73] parallel_4.2.3         Rcpp_1.0.14            vctrs_0.6.5           
#> [76] tidyselect_1.2.1       xfun_0.52
```
