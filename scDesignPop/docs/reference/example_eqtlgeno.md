# A dataframe containing both eQTL annotations and genotypes of individuals

A subset of data from the OneK1K dataset whose individual ids are
anonymized, and with simulated genotypes.

## Usage

``` r
data("example_eqtlgeno")
```

## Format

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
