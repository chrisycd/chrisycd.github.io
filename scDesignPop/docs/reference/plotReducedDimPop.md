# Dimensionality reduction and visualization for population-scale data

This function takes a reference sce and a list of new sces, performs the
dimensionality reduction on the reference data, projects the synthetic
datasets on the same low dimensional space (PCA and UMAP) of the
reference sce, and then visualize the results.

## Usage

``` r
plotReducedDimPop(
  ref_sce,
  sce_list,
  name_vec,
  assay_use = "logcounts",
  pc_umap = TRUE,
  n_pc = 50,
  center = TRUE,
  scale. = TRUE,
  if_plot = TRUE,
  shape_by = NULL,
  color_by,
  point_size = 1
)
```

## Arguments

- ref_sce:

  a reference sce object which synthetic sce objects will be projected
  to.

- sce_list:

  a list of synthetic sce objects.

- name_vec:

  a string vector specifiying the names of each dataset. The length
  should be `length(sce_list) + 1`, where the first name is for
  `ref_sce`.

- assay_use:

  a string scalar which indicates the assay you will use in the sce.
  Default is 'logcounts'.

- pc_umap:

  a boolean value specifying whether using PCs as the input of UMAP.
  Default is TRUE.

- n_pc:

  an integer specifying the number of PCs.

- center:

  a boolean value specifying whether centering the data before PCA.
  Default is TRUE.

- scale.:

  a boolean value specifying whether scaling the data before PCA.
  Default is TRUE.

- if_plot:

  a boolean value specifying whether returning the plot. If FALSE,
  return the reduced dimensions of each dataset.

- shape_by:

  a string scalar which indicates the column in `colData` used for
  shape.

- color_by:

  a string scalar which indicates the column in `colData` used for
  color.

- point_size:

  an numeric scalar specifying of the point size in the final plot.
  Default is 1.

## Value

A data frame of the reduced dimensions (both PCA and UMAP) or a list
contains both the data frame and two ggplot2 object of PCA plot and UMAP
plot.

## Examples

``` r
NULL
#> NULL
```
