# Example list of scDesignPop's marginal models

A named list of fitted marginal models for selected gene–SNP pairs used
in the examples of scDesignPop. Each element corresponds to one gene and
contains the fitted model object, the SNP used as covariate, and
additional model attributes.

## Usage

``` r
data("marginal_list_sel")
```

## Format

A named list. Each element corresponds to one gene (e.g.,
`"ENSG00000163221"`, `"ENSG00000135218"`) and is a list with the
following components:

- `fit`:

  A fitted generalized linear mixed model, typically a negative binomial
  mixed model (e.g., from glmmTMB) with random intercepts for
  individuals and fixed effects for cell type, SNP genotype, and
  SNP–cell-type interactions.

- `time`:

  Computation time (in seconds) required to fit the model.

- `snp_cov`:

  Character string giving the SNP covariate used in the model (e.g.,
  `"1:153337943"`).

- `model_attr`:

  A list of auxiliary model attributes, including the `terms` object,
  variable names, and data classes used for constructing the model
  matrix.

- `removed_cell`:

  Information about removed cells (if any) during model fitting; `NA` if
  no cells were removed.

## Details

This dataset is useful for demonstrating how to directly perform power
analysis using the fitted marginal models of the scDesignPop.

## Examples

``` r
data("marginal_list_sel")
names(marginal_list_sel)
#> [1] "ENSG00000163221" "ENSG00000135218"
str(marginal_list_sel[[1]]$fit)
#> List of 7
#>  $ obj      :List of 10
#>   ..$ par     : Named num [1:10] 0 0 0 0 0 0 0 0 0 0
#>   .. ..- attr(*, "names")= chr [1:10] "beta" "beta" "beta" "beta" ...
#>   ..$ fn      :function (x = last.par[lfixed()], ...)  
#>   ..$ gr      :function (x = last.par[lfixed()], ...)  
#>   ..$ he      :function (x = last.par[lfixed()], atomic = usingAtomics())  
#>   ..$ hessian : logi FALSE
#>   ..$ method  : chr "BFGS"
#>   ..$ retape  :function (set.defaults = TRUE)  
#>   ..$ env     :<environment: 0x5d9b212f8b70> 
#>   ..$ report  :function (par = last.par)  
#>   ..$ simulate:function (par = last.par, complete = FALSE)  
#>  $ fit      :List of 7
#>   ..$ par        : Named num [1:10] 1.8039 -4.6539 -5.9649 -6.4468 -0.0227 ...
#>   .. ..- attr(*, "names")= chr [1:10] "beta" "beta" "beta" "beta" ...
#>   ..$ objective  : num 5534
#>   ..$ convergence: int 0
#>   ..$ iterations : int 41
#>   ..$ evaluations: Named int [1:2] 52 42
#>   .. ..- attr(*, "names")= chr [1:2] "function" "gradient"
#>   ..$ message    : chr "relative convergence (4)"
#>   ..$ parfull    : Named num [1:50] 1.8039 -4.6539 -5.9649 -6.4468 -0.0227 ...
#>   .. ..- attr(*, "names")= chr [1:50] "beta" "beta" "beta" "beta" ...
#>   ..- attr(*, "optTime")= 'proc_time' Named num [1:5] 0.584 0.014 0.599 0 0
#>   .. ..- attr(*, "names")= chr [1:5] "user.self" "sys.self" "elapsed" "user.child" ...
#>  $ sdr      :List of 10
#>   ..$ value          : num(0) 
#>   ..$ sd             : num(0) 
#>   ..$ cov            : logi[0 , 0 ] 
#>   ..$ par.fixed      : Named num [1:10] 1.8039 -4.6539 -5.9649 -6.4468 -0.0227 ...
#>   .. ..- attr(*, "names")= chr [1:10] "beta" "beta" "beta" "beta" ...
#>   ..$ cov.fixed      : num [1:10, 1:10] 0.00603 -0.0014 -0.00241 -0.00191 -0.00452 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. .. ..$ : chr [1:10] "beta" "beta" "beta" "beta" ...
#>   .. .. ..$ : chr [1:10] "beta" "beta" "beta" "beta" ...
#>   ..$ pdHess         : logi TRUE
#>   ..$ gradient.fixed : num [1:10] -2.43e-03 2.80e-04 9.69e-05 1.48e-03 -1.47e-03 ...
#>   ..$ par.random     : Named num [1:40] -0.0836 -0.1549 -0.3842 0.1453 0.005 ...
#>   .. ..- attr(*, "names")= chr [1:40] "b" "b" "b" "b" ...
#>   ..$ diag.cov.random: num [1:40] 0.0613 0.0166 0.0123 0.0432 0.0538 ...
#>   ..$ env            :<environment: 0x5d9b20fce240> 
#>   ..- attr(*, "class")= chr "sdreport"
#>  $ call     : language glmmTMB::glmmTMB(formula = response ~ (1 | indiv) + cell_type + `1:153337943` +      `1:153337943`:cell_type, dat| __truncated__ ...
#>  $ frame    :'data.frame':   7811 obs. of  4 variables:
#>   ..$ response   : num [1:7811] 0 0 0 0 0 0 0 0 0 0 ...
#>   ..$ indiv      : Factor w/ 40 levels "SAMP1","SAMP10",..: 1 1 1 1 1 1 1 1 1 1 ...
#>   ..$ cell_type  : Factor w/ 4 levels "monoc","mononc",..: 4 4 4 4 4 4 4 4 4 4 ...
#>   ..$ 1:153337943: num [1:7811] 1 1 1 1 1 1 1 1 1 1 ...
#>   ..- attr(*, "terms")=Classes 'terms', 'formula'  language response ~ (1 + indiv) + cell_type + `1:153337943` + `1:153337943`:cell_type +      0 + 1
#>   .. .. ..- attr(*, "variables")= language list(response, indiv, cell_type, `1:153337943`)
#>   .. .. ..- attr(*, "factors")= int [1:4, 1:4] 0 1 0 0 0 0 1 0 0 0 ...
#>   .. .. .. ..- attr(*, "dimnames")=List of 2
#>   .. .. .. .. ..$ : chr [1:4] "response" "indiv" "cell_type" "`1:153337943`"
#>   .. .. .. .. ..$ : chr [1:4] "indiv" "cell_type" "`1:153337943`" "cell_type:`1:153337943`"
#>   .. .. ..- attr(*, "term.labels")= chr [1:4] "indiv" "cell_type" "`1:153337943`" "cell_type:`1:153337943`"
#>   .. .. ..- attr(*, "order")= int [1:4] 1 1 1 2
#>   .. .. ..- attr(*, "intercept")= int 1
#>   .. .. ..- attr(*, "response")= int 1
#>   .. .. ..- attr(*, ".Environment")=<environment: 0x5d9b20fce7b8> 
#>   .. .. ..- attr(*, "predvars")= language list(response, indiv, cell_type, `1:153337943`)
#>   .. .. ..- attr(*, "dataClasses")= Named chr [1:4] "numeric" "factor" "factor" "numeric"
#>   .. .. .. ..- attr(*, "names")= chr [1:4] "response" "indiv" "cell_type" "1:153337943"
#>  $ modelInfo:List of 15
#>   ..$ nobs          : int 7811
#>   ..$ respCol       : Named int 1
#>   .. ..- attr(*, "names")= chr "response"
#>   ..$ grpVar        : chr "indiv"
#>   ..$ family        :List of 11
#>   .. ..$ family    : chr "nbinom2"
#>   .. ..$ variance  :function (mu, theta = NULL)  
#>   .. ..$ initialize:  expression({  if (any(y < 0))  stop("negative values not allowed for the negative binomial family")  n <- rep(1, | __truncated__
#>   .. ..$ dev.resids:function (y, mu, wt, theta = NULL)  
#>   .. ..$ link      : chr "log"
#>   .. ..$ linkfun   :function (mu)  
#>   .. ..$ linkinv   :function (eta)  
#>   .. ..$ mu.eta    :function (eta)  
#>   .. ..$ valideta  :function (eta)  
#>   .. ..$ name      : chr "log"
#>   .. ..$ aic       :function (...)  
#>   .. ..- attr(*, "class")= chr "family"
#>   ..$ contrasts     : NULL
#>   ..$ reTrms        :List of 3
#>   .. ..$ cond:List of 3
#>   .. .. ..$ cnms :List of 1
#>   .. .. .. ..$ indiv: chr "(Intercept)"
#>   .. .. ..$ flist:List of 1
#>   .. .. .. ..$ indiv: Factor w/ 40 levels "SAMP1","SAMP10",..: 1 1 1 1 1 1 1 1 1 1 ...
#>   .. .. .. ..- attr(*, "assign")= int 1
#>   .. .. ..$ terms:List of 1
#>   .. .. .. ..$ fixed:Classes 'terms', 'formula'  language response ~ cell_type + `1:153337943` + `1:153337943`:cell_type
#>   .. .. .. .. .. ..- attr(*, "variables")= language list(response, cell_type, `1:153337943`)
#>   .. .. .. .. .. ..- attr(*, "factors")= int [1:3, 1:3] 0 1 0 0 0 1 0 1 1
#>   .. .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
#>   .. .. .. .. .. .. .. ..$ : chr [1:3] "response" "cell_type" "`1:153337943`"
#>   .. .. .. .. .. .. .. ..$ : chr [1:3] "cell_type" "`1:153337943`" "cell_type:`1:153337943`"
#>   .. .. .. .. .. ..- attr(*, "term.labels")= chr [1:3] "cell_type" "`1:153337943`" "cell_type:`1:153337943`"
#>   .. .. .. .. .. ..- attr(*, "order")= int [1:3] 1 1 2
#>   .. .. .. .. .. ..- attr(*, "intercept")= int 1
#>   .. .. .. .. .. ..- attr(*, "response")= int 1
#>   .. .. .. .. .. ..- attr(*, ".Environment")=<environment: 0x5d9b20fce7b8> 
#>   .. .. .. .. .. ..- attr(*, "predvars")= language list(response, cell_type, `1:153337943`)
#>   .. .. .. .. .. ..- attr(*, "dataClasses")= Named chr [1:3] "numeric" "factor" "numeric"
#>   .. .. .. .. .. .. ..- attr(*, "names")= chr [1:3] "response" "cell_type" "1:153337943"
#>   .. ..$ zi  :List of 1
#>   .. .. ..$ terms: NULL
#>   .. ..$ disp:List of 1
#>   .. .. ..$ terms:List of 1
#>   .. .. .. ..$ fixed:Classes 'terms', 'formula'  language ~1
#>   .. .. .. .. .. ..- attr(*, "variables")= language list()
#>   .. .. .. .. .. ..- attr(*, "factors")= int(0) 
#>   .. .. .. .. .. ..- attr(*, "term.labels")= chr(0) 
#>   .. .. .. .. .. ..- attr(*, "order")= int(0) 
#>   .. .. .. .. .. ..- attr(*, "intercept")= int 1
#>   .. .. .. .. .. ..- attr(*, "response")= int 0
#>   .. .. .. .. .. ..- attr(*, ".Environment")=<environment: 0x5d9b20fce7b8> 
#>   .. .. .. .. .. ..- attr(*, "predvars")= language list()
#>   .. .. .. .. .. ..- attr(*, "dataClasses")= Named chr(0) 
#>   .. .. .. .. .. .. ..- attr(*, "names")= chr(0) 
#>   ..$ terms         :List of 3
#>   .. ..$ cond:List of 1
#>   .. .. ..$ fixed:Classes 'terms', 'formula'  language response ~ cell_type + `1:153337943` + `1:153337943`:cell_type
#>   .. .. .. .. ..- attr(*, "variables")= language list(response, cell_type, `1:153337943`)
#>   .. .. .. .. ..- attr(*, "factors")= int [1:3, 1:3] 0 1 0 0 0 1 0 1 1
#>   .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
#>   .. .. .. .. .. .. ..$ : chr [1:3] "response" "cell_type" "`1:153337943`"
#>   .. .. .. .. .. .. ..$ : chr [1:3] "cell_type" "`1:153337943`" "cell_type:`1:153337943`"
#>   .. .. .. .. ..- attr(*, "term.labels")= chr [1:3] "cell_type" "`1:153337943`" "cell_type:`1:153337943`"
#>   .. .. .. .. ..- attr(*, "order")= int [1:3] 1 1 2
#>   .. .. .. .. ..- attr(*, "intercept")= int 1
#>   .. .. .. .. ..- attr(*, "response")= int 1
#>   .. .. .. .. ..- attr(*, ".Environment")=<environment: 0x5d9b20fce7b8> 
#>   .. .. .. .. ..- attr(*, "predvars")= language list(response, cell_type, `1:153337943`)
#>   .. .. .. .. ..- attr(*, "dataClasses")= Named chr [1:3] "numeric" "factor" "numeric"
#>   .. .. .. .. .. ..- attr(*, "names")= chr [1:3] "response" "cell_type" "1:153337943"
#>   .. ..$ zi  : NULL
#>   .. ..$ disp:List of 1
#>   .. .. ..$ fixed:Classes 'terms', 'formula'  language ~1
#>   .. .. .. .. ..- attr(*, "variables")= language list()
#>   .. .. .. .. ..- attr(*, "factors")= int(0) 
#>   .. .. .. .. ..- attr(*, "term.labels")= chr(0) 
#>   .. .. .. .. ..- attr(*, "order")= int(0) 
#>   .. .. .. .. ..- attr(*, "intercept")= int 1
#>   .. .. .. .. ..- attr(*, "response")= int 0
#>   .. .. .. .. ..- attr(*, ".Environment")=<environment: 0x5d9b20fce7b8> 
#>   .. .. .. .. ..- attr(*, "predvars")= language list()
#>   .. .. .. .. ..- attr(*, "dataClasses")= Named chr(0) 
#>   .. .. .. .. .. ..- attr(*, "names")= chr(0) 
#>   ..$ reStruc       :List of 3
#>   .. ..$ condReStruc:List of 1
#>   .. .. ..$ 1 | indiv:List of 6
#>   .. .. .. ..$ blockReps    : num 40
#>   .. .. .. ..$ blockSize    : num 1
#>   .. .. .. ..$ blockNumTheta: num 1
#>   .. .. .. ..$ blockCode    : Named num 1
#>   .. .. .. .. ..- attr(*, "names")= chr "us"
#>   .. .. .. ..$ simCode      : num 2
#>   .. .. .. ..$ fullCor      : int 1
#>   .. ..$ ziReStruc  : list()
#>   .. ..$ dispReStruc: list()
#>   ..$ allForm       :List of 4
#>   .. ..$ combForm   :Class 'formula'  language response ~ (1 + indiv) + cell_type + `1:153337943` + `1:153337943`:cell_type +      0 + 1
#>   .. .. .. ..- attr(*, ".Environment")=<environment: 0x5d9b20fce7b8> 
#>   .. ..$ formula    :Class 'formula'  language response ~ (1 | indiv) + cell_type + `1:153337943` + `1:153337943`:cell_type
#>   .. .. .. ..- attr(*, ".Environment")=<environment: 0x5d9b20fce7b8> 
#>   .. ..$ ziformula  :Class 'formula'  language ~0
#>   .. .. .. ..- attr(*, ".Environment")=<environment: 0x5d9b20fce7b8> 
#>   .. ..$ dispformula:Class 'formula'  language ~1
#>   .. .. .. ..- attr(*, ".Environment")=<environment: 0x5d9b20fce7b8> 
#>   ..$ REML          : logi FALSE
#>   ..$ map           : NULL
#>   ..$ sparseX       : Named logi [1:3] FALSE FALSE FALSE
#>   .. ..- attr(*, "names")= chr [1:3] "cond" "zi" "disp"
#>   ..$ parallel      :List of 2
#>   .. ..$ n      : int 1
#>   .. ..$ autopar: logi FALSE
#>   ..$ priors        : NULL
#>   ..$ packageVersion:Classes 'package_version', 'numeric_version'  hidden list of 1
#>   .. ..$ : int [1:3] 1 1 13
#>  $ fitted   : NULL
#>  - attr(*, "class")= chr "glmmTMB"
```
