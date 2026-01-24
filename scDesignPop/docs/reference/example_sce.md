# Example single-cell RNA-seq data (multi–cell type)

A small `SingleCellExperiment` object used as example input for
scDesignPop functions. It contains a subset of genes and cells from the
OneK1K cohort with associated individual and cell-type metadata.

A subset of the scRNA-seq data from the OneK1K cohort which has been
anonymized.

## Usage

``` r
data("example_sce")

data("example_sce")
```

## Format

A `SingleCellExperiment` object with 1,000 genes and 7,811 cells.

- assays:

  - `counts`: raw UMI count matrix (genes x cells).

  - `logcounts`: log-normalized expression values.

- rowData:

  (empty in this toy example).

- colData:

  A `DataFrame` with 4 columns:

  `indiv`

  :   Individual identifier (e.g., `"SAMP1"`, `"SAMP2"`).

  `cell_type`

  :   Cell-type annotation (e.g., T cells, B cells, NK cells).

  `sex`

  :   Sex of the donor.

  `age`

  :   Age of the donor.

  `batch`

  :   Batch of the sequencing (Simulated).

- reducedDimNames:

  None in this example.

### `example_sce`

A data frame with 100 genes (rows) and 17,918 cells (cols) for 20
individuals :

- cell_type:

  Cell type

- indiv:

  Individual id

## Examples

``` r
data("example_sce")
example_sce
#> class: SingleCellExperiment 
#> dim: 1000 7811 
#> metadata(0):
#> assays(2): counts logcounts
#> rownames(1000): ENSG00000023902 ENSG00000027869 ... ENSG00000254709
#>   ENSG00000272216
#> rowData names(0):
#> colnames: NULL
#> colData names(5): indiv cell_type sex age batch
#> reducedDimNames(0):
#> mainExpName: RNA
#> altExpNames(0):
table(SingleCellExperiment::colData(example_sce)$cell_type)
#> 
#>  monoc mononc   bmem  cd4nc 
#>   1750    689   1766   3606 
```
