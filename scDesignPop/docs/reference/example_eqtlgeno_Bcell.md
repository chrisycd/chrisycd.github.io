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

A tibble with 2,406 rows and 106 columns:

- `cell_type`:

  Cell type (e.g., `"cd4nc"`, `"cd8nc"`, `"nk"`, `"cd4et"`, etc.).

- `gene_name`:

  Gene symbol (e.g., `"PLEKHO1"`, `"TNFRSF1B"`, `"UTS2"`).

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
#> # A tibble: 2,406 × 106
#>    cell_type gene_name gene_id snp_id   CHR    POS SAMP1 SAMP2 SAMP3 SAMP4 SAMP5
#>    <chr>     <chr>     <chr>   <chr>  <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#>  1 cd4nc     PLEKHO1   ENSG00… 1:150…     1 1.50e8     0     1     0     0     0
#>  2 cd8nc     PLEKHO1   ENSG00… 1:150…     1 1.50e8     1     1     2     2     2
#>  3 cd4nc     TNFRSF1B  ENSG00… 1:121…     1 1.22e7     1     1     2     1     2
#>  4 cd8nc     TNFRSF1B  ENSG00… 1:122…     1 1.23e7     2     0     1     1     0
#>  5 nk        TNFRSF1B  ENSG00… 1:122…     1 1.23e7     2     2     1     1     1
#>  6 cd4et     UTS2      ENSG00… 1:817…     1 8.18e6     0     1     1     1     1
#>  7 cd4nc     UTS2      ENSG00… 1:798…     1 7.98e6     2     2     1     1     0
#>  8 cd8et     UTS2      ENSG00… 1:787…     1 7.87e6     1     2     2     2     2
#>  9 cd8nc     UTS2      ENSG00… 1:798…     1 7.98e6     2     2     1     2     1
#> 10 nk        UTS2      ENSG00… 1:797…     1 7.97e6     0     0     0     0     0
#> # ℹ 2,396 more rows
#> # ℹ 95 more variables: SAMP6 <dbl>, SAMP7 <dbl>, SAMP8 <dbl>, SAMP9 <dbl>,
#> #   SAMP10 <dbl>, SAMP11 <dbl>, SAMP12 <dbl>, SAMP13 <dbl>, SAMP14 <dbl>,
#> #   SAMP15 <dbl>, SAMP16 <dbl>, SAMP17 <dbl>, SAMP18 <dbl>, SAMP19 <dbl>,
#> #   SAMP20 <dbl>, SAMP21 <dbl>, SAMP22 <dbl>, SAMP23 <dbl>, SAMP24 <dbl>,
#> #   SAMP25 <dbl>, SAMP26 <dbl>, SAMP27 <dbl>, SAMP28 <dbl>, SAMP29 <dbl>,
#> #   SAMP30 <dbl>, SAMP31 <dbl>, SAMP32 <dbl>, SAMP33 <dbl>, SAMP34 <dbl>, …
dplyr::count(example_eqtlgeno_Bcell, gene_name)
#> # A tibble: 817 × 2
#>    gene_name      n
#>    <chr>      <int>
#>  1 ABCB9          1
#>  2 ABCC3          2
#>  3 AC002331.1     2
#>  4 AC006129.4     2
#>  5 AC013264.2     7
#>  6 AC016629.8     1
#>  7 AC018816.3     5
#>  8 AC069277.2     1
#>  9 AC079767.4     2
#> 10 AC092580.4     1
#> # ℹ 807 more rows
```
