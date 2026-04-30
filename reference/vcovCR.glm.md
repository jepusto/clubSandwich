# Cluster-robust variance-covariance matrix for a glm object.

`vcovCR` returns a sandwich estimate of the variance-covariance matrix
of a set of regression coefficient estimates from an
[`glm`](https://rdrr.io/r/stats/glm.html) object.

## Usage

``` r
# S3 method for class 'glm'
vcovCR(
  obj,
  cluster,
  type,
  target = NULL,
  inverse_var = NULL,
  form = "sandwich",
  ...
)
```

## Arguments

- obj:

  Fitted model for which to calculate the variance-covariance matrix

- cluster:

  Expression or vector indicating which observations belong to the same
  cluster. Required for `glm` objects.

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
  the target is taken to be the estimated variance function.

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

if (requireNamespace("geepack", quietly = TRUE)) {

  data(dietox, package = "geepack")
  dietox$Cu <- as.factor(dietox$Cu)
  weight_fit <- glm(Weight ~ Cu * poly(Time, 3), data=dietox, family = "quasipoisson")
  V_CR <- vcovCR(weight_fit, cluster = dietox$Pig, type = "CR2")
  coef_test(weight_fit, vcov = V_CR, test = "Satterthwaite")
  
}
#> Alternative hypothesis: two-sided 
#>                   Coef. Estimate     SE Null value  t-stat d.f. (Satt)
#>             (Intercept)   4.0124 0.0190          0 211.193        22.0
#>                 CuCu035  -0.0134 0.0286          0  -0.469        45.7
#>                 CuCu175   0.0330 0.0333          0   0.993        44.8
#>          poly(Time, 3)1  12.7115 0.2414          0  52.655        22.0
#>          poly(Time, 3)2  -1.6810 0.1456          0 -11.545        22.0
#>          poly(Time, 3)3   0.0292 0.0566          0   0.517        21.9
#>  CuCu035:poly(Time, 3)1  -0.0823 0.3120          0  -0.264        45.6
#>  CuCu175:poly(Time, 3)1  -0.3242 0.3433          0  -0.944        44.8
#>  CuCu035:poly(Time, 3)2   0.0927 0.2113          0   0.439        45.6
#>  CuCu175:poly(Time, 3)2  -0.1777 0.1656          0  -1.073        44.8
#>  CuCu035:poly(Time, 3)3  -0.1010 0.1013          0  -0.997        45.5
#>  CuCu175:poly(Time, 3)3   0.1146 0.0998          0   1.149        44.7
#>  p-val (Satt) Sig.
#>        <0.001  ***
#>         0.641     
#>         0.326     
#>        <0.001  ***
#>        <0.001  ***
#>         0.611     
#>         0.793     
#>         0.350     
#>         0.663     
#>         0.289     
#>         0.324     
#>         0.257     
```
