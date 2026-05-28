# Example genotype data for cell-type-specific trans-eQTLs

A tibble containing example genotype data for lead cis-eQTL SNPs across
multiple cell types and genes in OneK1K cohort. Each row corresponds to
a cell-type–gene–SNP combination, and each sample column stores genotype
dosages (0/1/2) for one individual. No eQTL effect size results are
included. Genotypes are all permuted.

## Usage

``` r
data("example_eqtlgeno_trans")
```

## Format

A tibble with 27 rows and 45 columns:

- `cell_type`:

  Cell type (`"bulk"`).

- `gene_id`:

  Ensembl gene ID.

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
data("example_eqtlgeno_trans")
example_eqtlgeno_trans
#> # A tibble: 27 × 45
#>    cell_type gene_id     snp_id   CHR    POS SAMP1 SAMP2 SAMP3 SAMP4 SAMP5 SAMP6
#>    <chr>     <chr>       <chr>  <int>  <int> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#>  1 bulk      ENSG000001… 10:94…    10 9.45e7     2     2     1     2     1     2
#>  2 bulk      ENSG000001… 2:182…     2 1.82e8     2     1     1     2     1     2
#>  3 bulk      ENSG000001… 10:94…    10 9.45e7     2     2     1     2     1     2
#>  4 bulk      ENSG000001… 10:94…    10 9.45e7     2     2     1     2     1     2
#>  5 bulk      ENSG000000… 10:94…    10 9.45e7     2     2     1     2     1     2
#>  6 bulk      ENSG000001… 10:94…    10 9.45e7     2     2     1     2     1     2
#>  7 bulk      ENSG000001… 10:94…    10 9.45e7     2     2     1     2     1     2
#>  8 bulk      ENSG000000… 2:182…     2 1.82e8     2     1     1     2     1     2
#>  9 bulk      ENSG000000… 10:94…    10 9.45e7     2     2     1     2     1     2
#> 10 bulk      ENSG000001… 10:94…    10 9.45e7     2     2     1     2     1     2
#> # ℹ 17 more rows
#> # ℹ 34 more variables: SAMP7 <dbl>, SAMP8 <dbl>, SAMP9 <dbl>, SAMP10 <dbl>,
#> #   SAMP11 <dbl>, SAMP12 <dbl>, SAMP13 <dbl>, SAMP14 <dbl>, SAMP15 <dbl>,
#> #   SAMP16 <dbl>, SAMP17 <dbl>, SAMP18 <dbl>, SAMP19 <dbl>, SAMP20 <dbl>,
#> #   SAMP21 <dbl>, SAMP22 <dbl>, SAMP23 <dbl>, SAMP24 <dbl>, SAMP25 <dbl>,
#> #   SAMP26 <dbl>, SAMP27 <dbl>, SAMP28 <dbl>, SAMP29 <dbl>, SAMP30 <dbl>,
#> #   SAMP31 <dbl>, SAMP32 <dbl>, SAMP33 <dbl>, SAMP34 <dbl>, SAMP35 <dbl>, …
```
