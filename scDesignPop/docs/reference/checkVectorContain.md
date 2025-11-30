# Check membership of first vector compared to other vectors

This is an internal helper function to check if elements of the first
vector are subsets of one or more other vectors. Ordering of elements is
ignored.

## Usage

``` r
checkVectorContain(..., ignore_dups = FALSE)
```

## Arguments

- ...:

  arguments of two or more vectors.

- ignore_dups:

  logical scalar to disregard duplicate elements during comparison.
  Default is FALSE.

## Value

a logical scalar

## Examples

``` r
vec1 <- c("cherry", "apple", "cherry")
vec2 <- c("apple", "cherry", "banana", "cherry")
vec3 <- c("banana", "kiwi", "apple", "cherry")
vec4 <- c("banana", "apple", "apple", "cherry")

checkVectorContain(vec1, vec2, vec3, ignore_dups = TRUE)   # returns TRUE
#> [1] TRUE
checkVectorContain(vec1, vec2, vec4, ignore_dups = FALSE)  # returns FALSE
#> [1] TRUE
```
