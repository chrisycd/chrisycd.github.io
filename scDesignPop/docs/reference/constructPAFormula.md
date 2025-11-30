# Construct a model formula for power analysis

Construct a model formula for power analysis

## Usage

``` r
constructPAFormula(
  fm,
  method = c("nb", "poisson", "gaussian", "pseudoBulkLinear"),
  snpid = NULL,
  celltype_colname = NULL
)
```

## Arguments

- fm:

  a stats::formula object from the full marginal model.

- method:

  a character object specifying the method that will be analyzed for
  power. (Options: nb,poisson,gaussian,pseudoBulkLinear).

- snpid:

  a character object contains snpid.

- celltype_colname:

  a string scalar specifying the cell state variable The default is
  "cell_type".

## Value

a stats::formula object for power analysis in each specified type.

## Examples

``` r
NULL
#> NULL
```
