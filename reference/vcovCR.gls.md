# Cluster-robust variance-covariance matrix for a gls object.

`vcovCR` returns a sandwich estimate of the variance-covariance matrix
of a set of regression coefficient estimates from a
[`gls`](https://rdrr.io/pkg/nlme/man/gls.html) object.

## Usage

``` r
# S3 method for class 'gls'
vcovCR(obj, cluster, type, target, inverse_var, form = "sandwich", ...)
```

## Arguments

- obj:

  Fitted model for which to calculate the variance-covariance matrix

- cluster:

  Optional expression or vector indicating which observations belong to
  the same cluster. If not specified, will be set to `getGroups(obj)`.

- type:

  Character string specifying which small-sample adjustment should be
  used, with available options `"CR0"`, `"CR1"`, `"CR1p"`, `"CR1S"`,
  `"CR2"`, or `"CR3"`. See "Details" section of
  [`vcovCR`](http://jepusto.github.io/clubSandwich/reference/vcovCR.md)
  for further information.

- target:

  Optional matrix or vector describing the working variance-covariance
  model used to calculate the `CR2` and `CR4` adjustment matrices. If
  not specified, the target is taken to be the estimated
  variance-covariance structure of the `gls` object.

- inverse_var:

  Optional logical indicating whether the weights used in fitting the
  model are inverse-variance. If not specified, `vcovCR` will attempt to
  infer a value.

- form:

  Controls the form of the returned matrix. The default `"sandwich"`
  will return the sandwich variance-covariance matrix. Alternately,
  setting `form = "meat"` will return only the meat of the sandwich and
  setting `form = B`, where `B` is a matrix of appropriate dimension,
  will return the sandwich variance-covariance matrix calculated using
  `B` as the bread. `form = "estfun"` will return the (appropriately
  scaled) estimating function, the transposed crossproduct of which is
  equal to the sandwich variance-covariance matrix.

- ...:

  Additional arguments available for some classes of objects.

## Value

An object of class `c("vcovCR","clubSandwich")`, which consists of a
matrix of the estimated variance of and covariances between the
regression coefficient estimates.

## See also

[`vcovCR`](http://jepusto.github.io/clubSandwich/reference/vcovCR.md)

## Examples

``` r

if (requireNamespace("nlme", quietly = TRUE)) {

  library(nlme)
  data(Ovary, package = "nlme")
  Ovary$time_int <- 1:nrow(Ovary)
  lm_AR1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), data = Ovary, 
                correlation = corAR1(form = ~ time_int | Mare))
  vcovCR(lm_AR1, type = "CR2")

}
#>                    (Intercept) sin(2 * pi * Time) cos(2 * pi * Time)
#> (Intercept)          1.0319747        -0.25050263        -0.28914347
#> sin(2 * pi * Time)  -0.2505026         0.31620081         0.02677139
#> cos(2 * pi * Time)  -0.2891435         0.02677139         0.16709681
    
```
