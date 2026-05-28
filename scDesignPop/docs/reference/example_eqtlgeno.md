# Example genotype data for cell-type-specific cis-eQTLs

A tibble containing example genotype data for lead cis-eQTL SNPs across
multiple cell types and genes in OneK1K cohort. Each row corresponds to
a cell-type–gene–SNP combination, and each sample column stores genotype
dosages (0/1/2) for one individual. No eQTL effect size results are
included. Genotypes are all permuted.

## Usage

``` r
data("example_eqtlgeno")
```

## Format

A tibble with 2,792 rows and 45 columns:

- `cell_type`:

  Cell type (e.g., `"cd4nc"`, `"cd8nc"`, `"nk"`, `"bin"`, `"bmem"`).

- `gene_id`:

  Ensembl gene ID (e.g., `"ENSG00000023902"`).

- `snp_id`:

  SNP identifier in `CHR:POS` format (e.g., `"1:150123456"`).

- `CHR`:

  Chromosome number.

- `POS`:

  Genomic position (base-pair coordinate).

- `SAMP1`, `SAMP2`, ..., `SAMP40`:

  Genotype dosage for each individual (typically encoded as 0, 1, or 2).

## Examples

``` r
data("example_eqtlgeno")
example_eqtlgeno
#> # A tibble: 2,792 × 45
#>    cell_type gene_id     snp_id   CHR    POS SAMP1 SAMP2 SAMP3 SAMP4 SAMP5 SAMP6
#>    <chr>     <chr>       <chr>  <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#>  1 monoc     ENSG000001… 1:101…     1 1.02e7     0     0     0     0     0     0
#>  2 monoc     ENSG000002… 1:108…     1 1.09e8     1     0     0     0     1     0
#>  3 cd4et     ENSG000001… 1:110…     1 1.11e7     2     2     2     2     2     2
#>  4 cd8et     ENSG000000… 1:111…     1 1.12e8     2     2     2     2     2     2
#>  5 nk        ENSG000000… 1:111…     1 1.12e8     2     2     2     2     2     2
#>  6 nkr       ENSG000000… 1:111…     1 1.12e8     2     2     2     2     2     2
#>  7 bin       ENSG000000… 1:111…     1 1.12e8     0     0     1     2     0     0
#>  8 bmem      ENSG000000… 1:111…     1 1.12e8     0     0     1     2     0     0
#>  9 cd8nc     ENSG000000… 1:111…     1 1.12e8     2     1     1     2     1     2
#> 10 cd4et     ENSG000000… 1:111…     1 1.12e8     1     1     2     1     1     2
#> # ℹ 2,782 more rows
#> # ℹ 34 more variables: SAMP7 <dbl>, SAMP8 <dbl>, SAMP9 <dbl>, SAMP10 <dbl>,
#> #   SAMP11 <dbl>, SAMP12 <dbl>, SAMP13 <dbl>, SAMP14 <dbl>, SAMP15 <dbl>,
#> #   SAMP16 <dbl>, SAMP17 <dbl>, SAMP18 <dbl>, SAMP19 <dbl>, SAMP20 <dbl>,
#> #   SAMP21 <dbl>, SAMP22 <dbl>, SAMP23 <dbl>, SAMP24 <dbl>, SAMP25 <dbl>,
#> #   SAMP26 <dbl>, SAMP27 <dbl>, SAMP28 <dbl>, SAMP29 <dbl>, SAMP30 <dbl>,
#> #   SAMP31 <dbl>, SAMP32 <dbl>, SAMP33 <dbl>, SAMP34 <dbl>, SAMP35 <dbl>, …
dplyr::count(example_eqtlgeno, cell_type)
#> # A tibble: 14 × 2
#>    cell_type     n
#>    <chr>     <int>
#>  1 bin         203
#>  2 bmem        179
#>  3 cd4et       196
#>  4 cd4nc       510
#>  5 cd4sox4      33
#>  6 cd8et       358
#>  7 cd8nc       271
#>  8 cd8s100b    121
#>  9 dc           77
#> 10 monoc       199
#> 11 mononc      177
#> 12 nk          361
#> 13 nkr          57
#> 14 plasma       50
```
