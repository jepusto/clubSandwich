# Detect cluster structure of an rma.mv object

`findCluster.rma.mv` returns a vector of ID variables for the highest
level of clustering in a fitted `rma.mv` model.

## Usage

``` r
findCluster.rma.mv(obj)
```

## Arguments

- obj:

  A fitted `rma.mv` object.

## Value

A a vector of ID variables for the highest level of clustering in `obj`.

## Examples

``` r
if (requireNamespace("metafor", quietly = TRUE)) {

library(metafor)
data(dat.assink2016, package = "metadat")

mfor_fit <- rma.mv(yi ~ year + deltype, 
                 V = vi, random = ~ 1 | study / esid,
                 data = dat.assink2016)
                 
findCluster.rma.mv(mfor_fit)

}
#> Loading required package: Matrix
#> Loading required package: metadat
#> Loading required package: numDeriv
#> 
#> Loading the 'metafor' package (version 5.0-1). For an
#> introduction to the package please type: help(metafor)
#>   [1] 1  1  1  1  1  1  2  2  2  3  3  3  3  3  3  3  3  3  3  4  5  6  6  6  6 
#>  [26] 6  7  7  7  7  7  7  8  9  9  10 10 10 10 10 11 11 11 11 11 11 11 11 11 11
#>  [51] 11 11 11 11 11 11 11 11 11 11 11 11 12 12 12 12 12 13 13 14 14 14 14 14 14
#>  [76] 15 15 15 16 16 16 16 16 16 16 16 16 16 16 16 16 16 16 16 17 17 17 17 17 17
#> Levels: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17
```
