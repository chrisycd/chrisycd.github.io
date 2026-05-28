# Example single-cell RNA-seq data (multi–cell type)

A small `SingleCellExperiment` object used as example input for
scDesignPop functions. It contains a subset of genes and cells from the
OneK1K cohort with associated individual and cell-type metadata.

## Usage

``` r
data("example_sce")
```

## Format

A `SingleCellExperiment` object with 982 genes and 7,998 cells.

- assays:

  - `counts`: raw UMI count matrix (genes x cells).

- rowData:

  (empty in this toy example).

- colData:

  A `DataFrame` with 4 columns:

  `indiv`

  :   Individual identifier (e.g., `"SAMP1"`, `"SAMP2"`).

  `pool`

  :   Batch of the donor.

  `cell_type`

  :   Cell-type annotation (e.g., T cells, B cells, NK cells).

  `sex`

  :   Sex of the donor.

  `age`

  :   Age of the donor.

- reducedDimNames:

  None in this example.

## Examples

``` r
data("example_sce")
example_sce
#> class: SingleCellExperiment 
#> dim: 982 7998 
#> metadata(0):
#> assays(1): counts
#> rownames(982): ENSG00000023902 ENSG00000028137 ... ENSG00000254709
#>   ENSG00000273272
#> rowData names(0):
#> colnames(7998): Cell1 Cell2 ... Cell7997 Cell7998
#> colData names(5): indiv pool cell_type sex age
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
table(SingleCellExperiment::colData(example_sce)$cell_type)
#> 
#>   bmem  cd4nc  monoc mononc 
#>   1640   3731   1952    675 
```
