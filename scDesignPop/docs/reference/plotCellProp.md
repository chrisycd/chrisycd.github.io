# Visualize the cell type proportions across individuals

Function visiualizes the simulated output of the cell type proportion
modelling function

## Usage

``` r
plotCellProp(
  col_data,
  title = NULL,
  color_vec = NULL,
  celltype_colname = "cell_type",
  indiv_colname = "indiv",
  width = 1,
  linewidth = 0.01
)
```

## Arguments

- col_data:

  a data frame with two columns: individual ids with column name
  `indiv_colname`, and cell types with column name `celltype_colname`.

- title:

  a string scalar specifiying the title of the output barplot.

- color_vec:

  a named vector of color hex codes named by the ordered cell type
  names. If NULL, the cell types contained in `col_data` will be ordered
  by alphabetical order with R default color schemes.

- celltype_colname:

  a string scalar to specify the cell type variable in `col_data`

- indiv_colname:

  a string scalar to specify the individual id variable in `col_data`.

- width:

  an numeric scalar specifying the width of each bar in the output
  barplot.

- linewidth:

  an numeric scalar specifying the line width of each bar in the output
  barplot.

## Value

outputs a ggplot2 object

## Examples

``` r
NULL
#> NULL
```
