# Example B cell single-cell RNA-seq data with pseudotime

A `SingleCellExperiment` object containing a subset of B cells from the
OneK1K cohort with associated individual metadata and a pseudotime
trajectory. This dataset is useful for illustrating dynamic eQTL
modeling along a continuous trajectory.

## Usage

``` r
data("example_sce_Bcell")
```

## Format

A `SingleCellExperiment` object with 753 genes and 3,785 cells.

- assays:

  - `counts`: raw UMI count matrix.

- rowData:

  A `DataFrame` with 2 columns:

  `GeneSymbol`

  :   Gene symbol.

  `features`

  :   Feature identifier (e.g., Ensembl ID).

- colData:

  A `DataFrame` with 5 columns:

  `cell_type`

  :   Cell-type annotation within the B-cell compartment (e.g.,
      immature/naive B cells, memory B cells).

  `sex`

  :   Batch of the donor.

  `indiv`

  :   Individual identifier.

  `sex`

  :   Sex of the donor.

  `age`

  :   Age of the donor.

  `slingPseudotime_1`

  :   Inferred pseudotime along a trajectory (e.g., from immature to
      memory B cells).

- reducedDimNames:

  `"PHATE"`: PHATE embedding used for trajectory inference and
  visualization.

## Examples

``` r
data("example_sce_Bcell")
example_sce_Bcell
#> class: SingleCellExperiment 
#> dim: 753 3785 
#> metadata(0):
#> assays(1): counts
#> rownames(753): ENSG00000023902 ENSG00000028137 ... ENSG00000244509
#>   ENSG00000254709
#> rowData names(2): GeneSymbol features
#> colnames(3785): Cell1 Cell2 ... Cell3784 Cell3785
#> colData names(5): indiv cell_type sex age slingPseudotime_1
#> reducedDimNames(1): PHATE
#> mainExpName: NULL
#> altExpNames(0):
summary(SingleCellExperiment::colData(example_sce_Bcell)$slingPseudotime_1)
#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#> 0.003699 0.317517 0.465263 0.474149 0.624296 1.000000 
```
