# Cluster-robust variance-covariance matrix for an `estimatr::iv_robust` object.

`vcovCR` returns a sandwich estimate of the variance-covariance matrix
of a set of regression coefficient estimates from an
[`iv_robust`](https://declaredesign.org/r/estimatr/reference/iv_robust.html)
object.

## Usage

``` r
# S3 method for class 'iv_robust'
vcovCR(
  obj,
  cluster,
  type,
  target = NULL,
  inverse_var = FALSE,
  form = "sandwich",
  ...
)
```

## Arguments

- obj:

  Fitted model for which to calculate the variance-covariance matrix

- cluster:

  Vector indicating which observations belong to the same cluster.
  Required for `iv_robust` objects unless the model was fitted with a
  `clusters` argument.

- type:

  Character string specifying which small-sample adjustment should be
  used, with available options `"CR0"`, `"CR1"`, `"CR1p"`, `"CR1S"`,
  `"CR2"`, or `"CR3"`. See "Details" section of
  [`vcovCR`](http://jepusto.github.io/clubSandwich/reference/vcovCR.md)
  for further information.

- target:

  Optional matrix or vector describing the working variance-covariance
  model used to calculate the `CR2` and `CR4` adjustment matrices. If a
  vector, the target matrix is assumed to be diagonal. If not specified,
  the target is taken to be an identity matrix.

- inverse_var:

  Not used for `iv_robust` objects.

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
if (requireNamespace("estimatr", quietly = TRUE)) withAutoprint({
  library(estimatr)

  set.seed(20261201)
  N <- 200
  cluster <- factor(rep(1:20, each = 10))
  z1 <- rnorm(N)
  z2 <- rnorm(N)
  x_exog <- rnorm(N)
  x_endog <- 0.5 * z1 + 0.3 * z2 + rnorm(N)
  y <- 1 + 2 * x_endog + 0.5 * x_exog + rnorm(N)
  dat <- data.frame(y, x_endog, x_exog, z1, z2, cluster)

  # Basic 2SLS with clusters and CR2
  iv_fit <- iv_robust(
    y ~ x_endog + x_exog | x_exog + z1 + z2,
    data = dat, 
    clusters = cluster, se_type = "CR2"
  )
  vcovCR(iv_fit)
  conf_int(iv_fit, vcov = "CR2")

  # 2SLS with fixed effects
  iv_fe <- iv_robust(
    y ~ x_endog + x_exog | x_exog + z1 + z2,
    fixed_effects = ~ cluster, 
    data = dat, 
    clusters = cluster,
    se_type = "CR2"
  )
  vcovCR(iv_fe)
  conf_int(iv_fe, vcov = "CR2")
})
#> > library(estimatr)
#> > set.seed(20261201)
#> > N <- 200
#> > cluster <- factor(rep(1:20, each = 10))
#> > z1 <- rnorm(N)
#> > z2 <- rnorm(N)
#> > x_exog <- rnorm(N)
#> > x_endog <- 0.5 * z1 + 0.3 * z2 + rnorm(N)
#> > y <- 1 + 2 * x_endog + 0.5 * x_exog + rnorm(N)
#> > dat <- data.frame(y, x_endog, x_exog, z1, z2, cluster)
#> > iv_fit <- iv_robust(y ~ x_endog + x_exog | x_exog + z1 + z2, data = dat, 
#> +     clusters = cluster, se_type = "CR2")
#> > vcovCR(iv_fit)
#>               (Intercept)      x_endog       x_exog
#> (Intercept)  0.0028574796 -0.001938895 0.0004459813
#> x_endog     -0.0019388949  0.012634349 0.0020154301
#> x_exog       0.0004459813  0.002015430 0.0059362835
#> > conf_int(iv_fit, vcov = "CR2")
#>        Coef. Estimate     SE d.f. Lower 95% CI Upper 95% CI
#>  (Intercept)    1.006 0.0535 18.8        0.894        1.118
#>      x_endog    2.126 0.1124 16.1        1.888        2.364
#>       x_exog    0.458 0.0770 14.3        0.293        0.623
#> > iv_fe <- iv_robust(y ~ x_endog + x_exog | x_exog + z1 + z2, fixed_effects = ~cluster, 
#> +     data = dat, clusters = cluster, se_type = "CR2")
#> > vcovCR(iv_fe)
#>             x_endog      x_exog
#> x_endog 0.013025078 0.002932599
#> x_exog  0.002932599 0.007977294
#> > conf_int(iv_fe, vcov = "CR2")
#>    Coef. Estimate     SE d.f. Lower 95% CI Upper 95% CI
#>  x_endog    2.115 0.1141 15.9         1.87        2.357
#>   x_exog    0.491 0.0893 14.1         0.30        0.683
```
