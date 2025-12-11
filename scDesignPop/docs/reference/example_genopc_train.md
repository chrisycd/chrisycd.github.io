# Genotype principal components for training individuals

Principal component (PC) scores computed from genotype data for a set of
training individuals. These PCs are typically used as covariates to
account for population structure in the cell type proportion modeling.

## Usage

``` r
data("example_genopc_train")
```

## Format

A tibble with 40 rows and 31 columns:

- `indiv`:

  Individual identifier (e.g., `"SAMP1"`, `"SAMP2"`).

- `PC1`, `PC2`, ..., `PC30`:

  Principal component scores summarizing genotype variation across
  individuals.

## Examples

``` r
data("example_genopc_train")
head(example_genopc_train)
#> # A tibble: 6 × 31
#>   indiv      PC1      PC2       PC3      PC4       PC5     PC6      PC7      PC8
#>   <chr>    <dbl>    <dbl>     <dbl>    <dbl>     <dbl>   <dbl>    <dbl>    <dbl>
#> 1 SAMP1  0.00988  0.0240   0.0552   -0.00774 -0.0239    0.0380  0.00775 -0.0179 
#> 2 SAMP2 -0.0217  -0.00495  0.00717  -0.0191   0.000858  0.0210 -0.00205 -0.0159 
#> 3 SAMP3  0.00503  0.0103  -0.00238  -0.0363   0.00505  -0.0245  0.0645  -0.0149 
#> 4 SAMP4 -0.0410   0.0297  -0.0341    0.0263  -0.0249    0.0171 -0.00949 -0.00863
#> 5 SAMP5 -0.00978  0.0294   0.0504   -0.0194  -0.0179    0.0368  0.0135   0.0115 
#> 6 SAMP6 -0.0212  -0.0106  -0.000811 -0.0106   0.0300    0.0547  0.00962 -0.0380 
#> # ℹ 22 more variables: PC9 <dbl>, PC10 <dbl>, PC11 <dbl>, PC12 <dbl>,
#> #   PC13 <dbl>, PC14 <dbl>, PC15 <dbl>, PC16 <dbl>, PC17 <dbl>, PC18 <dbl>,
#> #   PC19 <dbl>, PC20 <dbl>, PC21 <dbl>, PC22 <dbl>, PC23 <dbl>, PC24 <dbl>,
#> #   PC25 <dbl>, PC26 <dbl>, PC27 <dbl>, PC28 <dbl>, PC29 <dbl>, PC30 <dbl>
```
