# Genotype principal components for new individuals

Principal component (PC) scores computed from genotype data for a set of
new individuals. These individuals are projected into the same PC space
as the training set to be used in the cell type proportion modeling.

## Usage

``` r
data("example_genopc_new")
```

## Format

A tibble with 982 rows and 31 columns:

- `indiv`:

  Individual identifier (e.g., `"NEW_SAMP1"`).

- `PC1`, `PC2`, ..., `PC30`:

  Principal component scores for new individuals projected into the PC
  space derived from training individuals.

## Examples

``` r
data("example_genopc_new")
head(example_genopc_new)
#> # A tibble: 6 × 31
#>   indiv        PC1      PC2      PC3      PC4      PC5      PC6      PC7     PC8
#>   <chr>      <dbl>    <dbl>    <dbl>    <dbl>    <dbl>    <dbl>    <dbl>   <dbl>
#> 1 NEW_SA… -0.0167  -4.34e-2 -0.0174  -2.63e-2 -0.0391   1.36e-2  0.00141 -0.0226
#> 2 NEW_SA… -0.0182   2.96e-2 -0.0114  -5.91e-2 -0.0170   6.63e-3 -0.00878 -0.0453
#> 3 NEW_SA… -0.0247  -1.73e-2 -0.00227 -8.02e-4  0.0589  -1.50e-2 -0.0375   0.0260
#> 4 NEW_SA…  0.0199  -4.28e-2 -0.0490   5.99e-2 -0.0323   1.22e-2  0.0627   0.0562
#> 5 NEW_SA…  0.0743  -3.98e-2 -0.0151   3.53e-2 -0.00486  4.28e-2 -0.0277  -0.0236
#> 6 NEW_SA… -0.00289 -8.21e-4 -0.0423   6.06e-4  0.0245  -9.70e-4 -0.0200  -0.0487
#> # ℹ 22 more variables: PC9 <dbl>, PC10 <dbl>, PC11 <dbl>, PC12 <dbl>,
#> #   PC13 <dbl>, PC14 <dbl>, PC15 <dbl>, PC16 <dbl>, PC17 <dbl>, PC18 <dbl>,
#> #   PC19 <dbl>, PC20 <dbl>, PC21 <dbl>, PC22 <dbl>, PC23 <dbl>, PC24 <dbl>,
#> #   PC25 <dbl>, PC26 <dbl>, PC27 <dbl>, PC28 <dbl>, PC29 <dbl>, PC30 <dbl>
```
