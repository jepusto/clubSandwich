# Cluster-robust variance-covariance matrix for an lme object.

`vcovCR` returns a sandwich estimate of the variance-covariance matrix
of a set of regression coefficient estimates from a
[`lme`](https://rdrr.io/pkg/nlme/man/lme.html) object.

## Usage

``` r
# S3 method for class 'lme'
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
  variance-covariance structure of the `lme` object.

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
  rat_weight <- lme(weight ~ Time * Diet, data=BodyWeight, ~ Time | Rat) 
  vcovCR(rat_weight, type = "CR2")

}
#>              (Intercept)          Time         Diet2        Diet3    Time:Diet2
#> (Intercept)  20.50523572 -0.0586895628  -20.50523572 -20.50523572  0.0586895628
#> Time         -0.05868956  0.0008918201    0.05868956   0.05868956 -0.0008918201
#> Diet2       -20.50523572  0.0586895628 1192.99481298  20.50523572  0.3728534881
#> Diet3       -20.50523572  0.0586895628   20.50523572 238.25477649 -0.0586895628
#> Time:Diet2    0.05868956 -0.0008918201    0.37285349  -0.05868956  0.0458806483
#> Time:Diet3    0.05868956 -0.0008918201   -0.05868956  -1.87063320  0.0008918201
#>                Time:Diet3
#> (Intercept)  0.0586895628
#> Time        -0.0008918201
#> Diet2       -0.0586895628
#> Diet3       -1.8706331980
#> Time:Diet2   0.0008918201
#> Time:Diet3   0.0237304192

pkgs_available <- 
  requireNamespace("nlme", quietly = TRUE) & 
  requireNamespace("mlmRev", quietly = TRUE)

if (pkgs_available) {

  data(egsingle, package = "mlmRev")
  subset_ids <- levels(egsingle$schoolid)[1:10]
  egsingle_subset <- subset(egsingle, schoolid %in% subset_ids)
  
  math_model <- lme(math ~ year * size + female + black + hispanic, 
                    random = list(~ year | schoolid, ~ 1 | childid), 
                    data = egsingle_subset)
                    
  vcovCR(math_model, type = "CR2")
  
}
#>               (Intercept)          year          size    femaleMale
#> (Intercept)  4.398970e-02  1.018504e-02 -5.251581e-05 -7.784067e-03
#> year         1.018504e-02  2.000100e-02 -6.243677e-06 -3.116308e-03
#> size        -5.251581e-05 -6.243677e-06  9.759313e-08  9.072820e-06
#> femaleMale  -7.784067e-03 -3.116308e-03  9.072820e-06  8.999071e-03
#> black1      -6.030535e-03 -5.258604e-03 -3.438467e-06 -5.458323e-03
#> hispanic1    1.347498e-02  4.667589e-03 -1.564763e-05 -4.631514e-03
#> year:size   -1.206983e-05 -2.477655e-05  9.626510e-09  1.923760e-06
#>                    black1     hispanic1     year:size
#> (Intercept) -6.030535e-03  1.347498e-02 -1.206983e-05
#> year        -5.258604e-03  4.667589e-03 -2.477655e-05
#> size        -3.438467e-06 -1.564763e-05  9.626510e-09
#> femaleMale  -5.458323e-03 -4.631514e-03  1.923760e-06
#> black1       2.304693e-02  5.559583e-03  8.910563e-06
#> hispanic1    5.559583e-03  1.063522e-02 -4.438542e-06
#> year:size    8.910563e-06 -4.438542e-06  3.239673e-08
```
