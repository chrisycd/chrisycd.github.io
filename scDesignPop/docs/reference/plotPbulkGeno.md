# Plot pseudobulk expression vs. genotype

Creates a faceted violin plot of pseudobulk expression values across
genotype groups, with optional boxplots and mean summary points.

## Usage

``` r
plotPbulkGeno(
  res_df,
  celltype_colname,
  snp_colname,
  show_boxplot = TRUE,
  show_meanstat = TRUE,
  show_snp = FALSE,
  ...
)
```

## Arguments

- res_df:

  A dataframe

- celltype_colname:

  A string scalar for column name of cell type. The default is
  "cell_type".

- snp_colname:

  A string scalar for column name of the SNP variable. The default is
  "snp_id".

- show_boxplot:

  A logical scalar for whether to overlay boxplots.

- show_meanstat:

  A logical scalar for whether to display the mean marked by a red dot.

- show_snp:

  A logical scalar for whether to display selected SNP for each y-axis
  facet.

- ...:

  Additional arguments passed to internal functions.

## Value

A ggplot object.

## Examples

``` r
NULL
#> NULL
```
