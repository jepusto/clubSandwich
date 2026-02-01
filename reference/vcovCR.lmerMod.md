# Cluster-robust variance-covariance matrix for an lmerMod object.

`vcovCR` returns a sandwich estimate of the variance-covariance matrix
of a set of regression coefficient estimates from
[`merMod`](https://rdrr.io/pkg/lme4/man/merMod-class.html) object.

## Usage

``` r
# S3 method for class 'lmerMod'
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
  variance-covariance structure of the `lmerMod` object.

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
if (requireNamespace("lme4", quietly = TRUE)) {

library(lme4)
sleep_fit <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
vcovCR(sleep_fit, type = "CR2")

}
#> 
#> Attaching package: ‘lme4’
#> The following object is masked from ‘package:nlme’:
#> 
#>     lmList
#>             (Intercept)      Days
#> (Intercept)   46.574572 -1.451096
#> Days          -1.451096  2.389463

pkgs_available <- 
  requireNamespace("lme4", quietly = TRUE) & 
  requireNamespace("mlmRev", quietly = TRUE)

if (pkgs_available) {

data(egsingle, package = "mlmRev")
subset_ids <- levels(egsingle$schoolid)[1:10]
math_model <- lmer(math ~ year * size + female + black + hispanic 
                   + (1 | schoolid) + (1 | childid), 
                   data = egsingle, subset = schoolid %in% subset_ids)
vcovCR(math_model, type = "CR2")
}
#> Warning: Some predictor variables are on very different scales: consider rescaling. 
#> You may also use (g)lmerControl(autoscale = TRUE) to improve numerical stability.
#>               (Intercept)          year          size    femaleMale
#> (Intercept)  3.783508e-02  4.761442e-03 -4.339287e-05 -7.507242e-03
#> year         4.761442e-03  1.055863e-02 -2.769454e-06 -3.659318e-03
#> size        -4.339287e-05 -2.769454e-06  8.544310e-08  6.516505e-06
#> femaleMale  -7.507242e-03 -3.659318e-03  6.516505e-06  9.417935e-03
#> black1      -7.037052e-03  3.942216e-04 -1.389629e-06 -4.086171e-03
#> hispanic1    1.603053e-02  6.423362e-03 -1.644442e-05 -4.676828e-03
#> year:size   -6.313817e-06 -1.201793e-05  8.181346e-09  3.644258e-06
#>                    black1     hispanic1     year:size
#> (Intercept) -7.037052e-03  1.603053e-02 -6.313817e-06
#> year         3.942216e-04  6.423362e-03 -1.201793e-05
#> size        -1.389629e-06 -1.644442e-05  8.181346e-09
#> femaleMale  -4.086171e-03 -4.676828e-03  3.644258e-06
#> black1       2.086581e-02  3.779801e-03  9.423158e-07
#> hispanic1    3.779801e-03  1.375137e-02 -6.919894e-06
#> year:size    9.423158e-07 -6.919894e-06  1.553385e-08
```
