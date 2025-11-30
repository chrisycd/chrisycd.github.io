# Construct eQTL annotation and genotype dataframe

This function creates a combined dataframe using input eQTL annotations,
and sample genotype dataframe.

## Usage

``` r
constructEqtlGeno(
  eqtl_annot_df,
  geno_df,
  sampid_vec,
  name = NULL,
  feature_colname = "gene_name",
  cellstate_colname = "cell_type",
  snp_colname = "snp_id",
  loc_colname = "POS"
)
```

## Arguments

- eqtl_annot_df:

  eQTL annotation dataframe

- geno_df:

  eQTL by sample genotype dataframe

- sampid_vec:

  sample id vector

- name:

  string or integer scalar to identify input

- feature_colname:

  feature variable name

- cellstate_colname:

  cell state variable name

- snp_colname:

  SNP id variable name

- loc_colname:

  SNP position variable name

## Value

a dataframe with eQTL annotations and genotype of SNPs

## Examples

``` r
NULL
#> NULL
```
