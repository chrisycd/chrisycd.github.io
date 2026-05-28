# Example B cell eQTL genotype data

A tibble containing example genotype data for selected cis-eQTL SNPs for
the selected B cells from OneK1K. This dataset can be used together with
`example_sce_Bcell` to show dynamic or cell-type–specific eQTL modeling.
No eQTL effect size results are included. Genotypes are all permuted.

## Usage

``` r
data("example_eqtlgeno_Bcell")
```

## Format

A tibble with 2,332 rows and 105 columns:

- `cell_type`:

  Cell type.

- `gene_id`:

  Ensembl gene ID.

- `snp_id`:

  SNP identifier in `CHR:POS` format.

- `CHR`:

  Chromosome number.

- `POS`:

  Genomic position (base-pair coordinate).

- `SAMP1`, `SAMP2`, ...:

  Genotype dosage (0/1/2) for each individual. The remaining columns
  `SAMPk` store genotypes for all samples used in the example.

## Examples

``` r
data("example_eqtlgeno_Bcell")
example_eqtlgeno_Bcell
#> # A tibble: 2,332 × 105
#>    cell_type gene_id     snp_id   CHR    POS SAMP1 SAMP2 SAMP3 SAMP4 SAMP5 SAMP6
#>    <chr>     <chr>       <chr>  <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#>  1 monoc     ENSG000001… 1:101…     1 1.02e7     0     1     0     0     0     0
#>  2 cd4et     ENSG000001… 1:110…     1 1.11e7     2     2     2     2     2     2
#>  3 cd8et     ENSG000000… 1:111…     1 1.12e8     1     2     2     1     2     2
#>  4 nk        ENSG000000… 1:111…     1 1.12e8     1     2     2     1     2     2
#>  5 nkr       ENSG000000… 1:111…     1 1.12e8     1     2     2     1     2     2
#>  6 bin       ENSG000000… 1:111…     1 1.12e8     2     1     0     1     0     1
#>  7 bmem      ENSG000000… 1:111…     1 1.12e8     2     1     0     1     0     1
#>  8 cd8nc     ENSG000000… 1:111…     1 1.12e8     1     1     2     2     1     2
#>  9 cd4et     ENSG000000… 1:111…     1 1.12e8     2     1     1     1     1     2
#> 10 cd4nc     ENSG000000… 1:111…     1 1.12e8     2     2     0     1     1     1
#> # ℹ 2,322 more rows
#> # ℹ 94 more variables: SAMP7 <dbl>, SAMP8 <dbl>, SAMP9 <dbl>, SAMP10 <dbl>,
#> #   SAMP11 <dbl>, SAMP12 <dbl>, SAMP13 <dbl>, SAMP14 <dbl>, SAMP15 <dbl>,
#> #   SAMP16 <dbl>, SAMP17 <dbl>, SAMP18 <dbl>, SAMP19 <dbl>, SAMP20 <dbl>,
#> #   SAMP21 <dbl>, SAMP22 <dbl>, SAMP23 <dbl>, SAMP24 <dbl>, SAMP25 <dbl>,
#> #   SAMP26 <dbl>, SAMP27 <dbl>, SAMP28 <dbl>, SAMP29 <dbl>, SAMP30 <dbl>,
#> #   SAMP31 <dbl>, SAMP32 <dbl>, SAMP33 <dbl>, SAMP34 <dbl>, SAMP35 <dbl>, …
dplyr::count(example_eqtlgeno_Bcell, gene_id)
#> # A tibble: 753 × 2
#>    gene_id             n
#>    <chr>           <int>
#>  1 ENSG00000002549     3
#>  2 ENSG00000002933     1
#>  3 ENSG00000003147     1
#>  4 ENSG00000004468     1
#>  5 ENSG00000004809     1
#>  6 ENSG00000005020     6
#>  7 ENSG00000006075     3
#>  8 ENSG00000007264     2
#>  9 ENSG00000007312     6
#> 10 ENSG00000008517     8
#> # ℹ 743 more rows
```
