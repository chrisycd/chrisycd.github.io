# Visualize the power analysis result

Visualize the power analysis result

## Usage

``` r
visualizePowerResult(
  power_result,
  celltype_vector,
  x_axis = "nindiv",
  x_facet = "celltype",
  y_facet = "ncell",
  col_group = "method",
  geneid,
  snpid
)
```

## Arguments

- power_result:

  a data frame contains power analysis result in different parameter
  settings.

- celltype_vector:

  a vector of cell types that will be selected for visualization.

- x_axis:

  a character that specifies the x axis. Default is the number of
  individuals.

- x_facet:

  a character that specifies the x axis facet. Default is the cell types
  given in the celltype_vector.

- y_facet:

  a character that specifies the y axis facet. Default is the number of
  cells per individual.

- col_group:

  a character that specifies the color groups. Default is the eQTL
  model.

- geneid:

  a character object contains geneid to be part of the plot title.

- snpid:

  a character object contains snpid to be part of the plot title.

## Value

a ggplot object

## Examples

``` r
NULL
#> NULL
```
