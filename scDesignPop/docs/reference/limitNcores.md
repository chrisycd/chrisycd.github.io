# Limit number of cores for parallel computation

This is an internal helper function to determine a safe number of CPU
cores for parallel execution. Under `R CMD check` or CRAN-like
environments, the number of cores is automatically restricted to 1 to
avoid check failures. Otherwise, the number of cores is capped by the
available cores on the system.

## Usage

``` r
limitNcores(ncores = 1L, max_cores = Inf)
```

## Arguments

- ncores:

  integer scalar specifying the requested number of cores.

- max_cores:

  optional integer scalar specifying an upper bound on cores.

## Value

an integer scalar representing the number of cores to use.
