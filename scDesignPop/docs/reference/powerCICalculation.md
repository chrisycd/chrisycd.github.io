# Calculate a bootstrap confidence interval for each power

Calculate a bootstrap confidence interval for each power

## Usage

``` r
powerCICalculation(
  res,
  celltype_vector,
  nindivs,
  ncells,
  snp_number = 10,
  gene_number = 800,
  alpha = 0.05,
  nsim = 1000,
  conf = 0.05
)
```

## Arguments

- res:

  the output of function powerAnalysis()

- celltype_vector:

  a vector object specifies the cell type that will be tested

- nindivs:

  a vector of numeric values showing the numbers of individuals that
  user wants to simulate.

- ncells:

  a vector of numeric values showing the numbers of cells per each
  individual that user wants to simulate.

- snp_number:

  the number of SNPs for multiple testing correction.

- gene_number:

  the number of genes for multiple testing correction.

- alpha:

  the p value threshold for rejecting the H0 hypothesis.

- nsim:

  number of simulations for calculating the Bootstrap CI.

- conf:

  Bootstrap CI interval.

## Value

a data frame contains average power with standard deviations in each
parameter settings.

## Examples

``` r
NULL
#> NULL
```
