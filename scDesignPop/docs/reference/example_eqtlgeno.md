# Example genotype data for cell-type-specific eQTLs

A tibble containing example genotype data for lead cis-eQTL SNPs across
multiple cell types and genes in OneK1K cohort. Each row corresponds to
a cell-type–gene–SNP combination, and each sample column stores genotype
dosages (0/1/2) for one individual. No eQTL effect size results are
included. Genotypes are all permuted.

A subset of data from the OneK1K dataset whose individual ids are
anonymized, and with simulated genotypes.

## Usage

``` r
data("example_eqtlgeno")

data("example_eqtlgeno")
```

## Format

A tibble with 2,826 rows and 45 columns:

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

### `example_eqtlgeno`

A data frame with 270 cell-type specific eQTLs (rows), 5 eQTL relevant
annotations (cols) consisting of `cell_type`, `gene_id`, `snp_id`,
`CHR`, and `POS`, and 20 samples' genotypes (cols).

- cell_type:

  Abbreviated cell type

- gene_id:

  Gene ID

- snp_id:

  SNP id of the cell-type specific eQTL

- CHR:

  Chromosome

- POS:

  Chromosome location of SNP

## Examples

``` r
data("example_eqtlgeno")
example_eqtlgeno
#> # A tibble: 2,826 × 45
#>    cell_type gene_id     snp_id   CHR    POS SAMP1 SAMP2 SAMP3 SAMP4 SAMP5 SAMP6
#>    <chr>     <chr>       <chr>  <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#>  1 cd4nc     ENSG000000… 1:150…     1 1.50e8     0     0     0     0     0     1
#>  2 cd8nc     ENSG000000… 1:150…     1 1.50e8     2     2     1     1     1     2
#>  3 cd4nc     ENSG000000… 1:156…     1 1.57e8     0     0     0     0     1     0
#>  4 cd4nc     ENSG000000… 1:121…     1 1.22e7     2     2     1     1     1     0
#>  5 cd8nc     ENSG000000… 1:122…     1 1.23e7     2     1     1     2     1     2
#>  6 nk        ENSG000000… 1:122…     1 1.23e7     1     0     1     0     0     1
#>  7 bin       ENSG000000… 1:111…     1 1.12e8     0     1     1     0     1     0
#>  8 bmem      ENSG000000… 1:111…     1 1.12e8     1     0     0     0     0     0
#>  9 cd4et     ENSG000000… 1:111…     1 1.12e8     1     2     0     2     1     1
#> 10 cd4nc     ENSG000000… 1:111…     1 1.12e8     2     2     2     2     1     2
#> # ℹ 2,816 more rows
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
#>  1 bin         216
#>  2 bmem        186
#>  3 cd4et       186
#>  4 cd4nc       545
#>  5 cd4sox4      28
#>  6 cd8et       342
#>  7 cd8nc       289
#>  8 cd8s100b    126
#>  9 dc           71
#> 10 monoc       207
#> 11 mononc      189
#> 12 nk          342
#> 13 nkr          57
#> 14 plasma       42
```
