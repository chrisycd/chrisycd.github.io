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

A `SingleCellExperiment` object with 817 genes and 3,726 cells.

- assays:

  - `counts`: raw UMI count matrix.

  - `norm`: normalized expression values.

  - `logcounts`: log-transformed normalized expression values.

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
#> dim: 817 3726 
#> metadata(0):
#> assays(3): counts norm logcounts
#> rownames(817): ENSG00000023902 ENSG00000028137 ... ENSG00000254709
#>   ENSG00000273272
#> rowData names(2): GeneSymbol features
#> colnames(3726): AAAGATGGTTATGCGT-1 AACTCAGGTCCGCTGA-1 ...
#>   TTCGAAGGTGATGCCC-9 TTGCCGTGTTCGTCTC-9
#> colData names(5): cell_type indiv sex age slingPseudotime_1
#> reducedDimNames(1): PHATE
#> mainExpName: NULL
#> altExpNames(0):
summary(SingleCellExperiment::colData(example_sce_Bcell)$slingPseudotime_1)
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#> 0.0006717 0.2882431 0.4254734 0.4389588 0.5794676 0.9997564 
```
