# Cluster-robust variance-covariance matrix for a geeglm object.

`vcovCR` returns a sandwich estimate of the variance-covariance matrix
of a set of regression coefficient estimates from an
[`geeglm`](https://rdrr.io/pkg/geepack/man/geeglm.html) object.

## Usage

``` r
# S3 method for class 'geeglm'
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
  cluster. Required for `geeglm` objects.

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

  library(geepack)
  data(dietox, package = "geepack")
  dietox$Cu <- as.factor(dietox$Cu)
  mf <- formula(Weight ~ Cu * (Time + I(Time^2) + I(Time^3)))
  gee1 <- geeglm(mf, data=dietox, id=Pig, family=poisson("identity"), corstr="ar1")
  V_CR <- vcovCR(gee1, cluster = dietox$Pig, type = "CR2")
  coef_test(gee1, vcov = V_CR, test = "Satterthwaite")
  
}
#> Warning: non-integer x = 26.500000
#> Warning: non-integer x = 27.599990
#> Warning: non-integer x = 36.500000
#> Warning: non-integer x = 40.299990
#> Warning: non-integer x = 49.099980
#> Warning: non-integer x = 55.399990
#> Warning: non-integer x = 59.599980
#> Warning: non-integer x = 76.599980
#> Warning: non-integer x = 86.500000
#> Warning: non-integer x = 91.599980
#> Warning: non-integer x = 98.599980
#> Warning: non-integer x = 28.299990
#> Warning: non-integer x = 30.099990
#> Warning: non-integer x = 38.299990
#> Warning: non-integer x = 44.500000
#> Warning: non-integer x = 51.599980
#> Warning: non-integer x = 57.599980
#> Warning: non-integer x = 82.299990
#> Warning: non-integer x = 99.699950
#> Warning: non-integer x = 106.699950
#> Warning: non-integer x = 27.599990
#> Warning: non-integer x = 30.599990
#> Warning: non-integer x = 38.699980
#> Warning: non-integer x = 47.199980
#> Warning: non-integer x = 54.099980
#> Warning: non-integer x = 61.500000
#> Warning: non-integer x = 68.500000
#> Warning: non-integer x = 75.199950
#> Warning: non-integer x = 81.699950
#> Warning: non-integer x = 90.199950
#> Warning: non-integer x = 98.399960
#> Warning: non-integer x = 105.399960
#> Warning: non-integer x = 31.500000
#> Warning: non-integer x = 34.799990
#> Warning: non-integer x = 40.699980
#> Warning: non-integer x = 47.699980
#> Warning: non-integer x = 55.899990
#> Warning: non-integer x = 62.199980
#> Warning: non-integer x = 70.699950
#> Warning: non-integer x = 86.199950
#> Warning: non-integer x = 94.699950
#> Warning: non-integer x = 109.199950
#> Warning: non-integer x = 27.099990
#> Warning: non-integer x = 42.500000
#> Warning: non-integer x = 50.099980
#> Warning: non-integer x = 56.500000
#> Warning: non-integer x = 72.500000
#> Warning: non-integer x = 80.500000
#> Warning: non-integer x = 108.399960
#> Warning: non-integer x = 31.799990
#> Warning: non-integer x = 44.799990
#> Warning: non-integer x = 50.899990
#> Warning: non-integer x = 57.399990
#> Warning: non-integer x = 62.500000
#> Warning: non-integer x = 71.699950
#> Warning: non-integer x = 78.199950
#> Warning: non-integer x = 85.799990
#> Warning: non-integer x = 91.799990
#> Warning: non-integer x = 100.199950
#> Warning: non-integer x = 27.700000
#> Warning: non-integer x = 33.599980
#> Warning: non-integer x = 46.699980
#> Warning: non-integer x = 54.799990
#> Warning: non-integer x = 69.799990
#> Warning: non-integer x = 83.599980
#> Warning: non-integer x = 88.199950
#> Warning: non-integer x = 96.799990
#> Warning: non-integer x = 102.299990
#> Warning: non-integer x = 23.599990
#> Warning: non-integer x = 36.299990
#> Warning: non-integer x = 42.399990
#> Warning: non-integer x = 50.500000
#> Warning: non-integer x = 57.599980
#> Warning: non-integer x = 66.599980
#> Warning: non-integer x = 74.199950
#> Warning: non-integer x = 89.299990
#> Warning: non-integer x = 26.899990
#> Warning: non-integer x = 32.299990
#> Warning: non-integer x = 38.599980
#> Warning: non-integer x = 57.799990
#> Warning: non-integer x = 58.599980
#> Warning: non-integer x = 65.899960
#> Warning: non-integer x = 71.199950
#> Warning: non-integer x = 79.399960
#> Warning: non-integer x = 91.699950
#> Warning: non-integer x = 22.599990
#> Warning: non-integer x = 28.500000
#> Warning: non-integer x = 32.699980
#> Warning: non-integer x = 39.899990
#> Warning: non-integer x = 47.599980
#> Warning: non-integer x = 53.799990
#> Warning: non-integer x = 67.199950
#> Warning: non-integer x = 74.099980
#> Warning: non-integer x = 89.399960
#> Warning: non-integer x = 93.500000
#> Warning: non-integer x = 19.399990
#> Warning: non-integer x = 21.399990
#> Warning: non-integer x = 25.299990
#> Warning: non-integer x = 31.899990
#> Warning: non-integer x = 37.599980
#> Warning: non-integer x = 43.699980
#> Warning: non-integer x = 48.899990
#> Warning: non-integer x = 53.699980
#> Warning: non-integer x = 59.299990
#> Warning: non-integer x = 64.899960
#> Warning: non-integer x = 65.799990
#> Warning: non-integer x = 27.299990
#> Warning: non-integer x = 32.399990
#> Warning: non-integer x = 38.299990
#> Warning: non-integer x = 45.799990
#> Warning: non-integer x = 53.799990
#> Warning: non-integer x = 68.699950
#> Warning: non-integer x = 75.799990
#> Warning: non-integer x = 80.799990
#> Warning: non-integer x = 86.799990
#> Warning: non-integer x = 91.399960
#> Warning: non-integer x = 99.799990
#> Warning: non-integer x = 26.200000
#> Warning: non-integer x = 31.099990
#> Warning: non-integer x = 34.599980
#> Warning: non-integer x = 41.199980
#> Warning: non-integer x = 49.199980
#> Warning: non-integer x = 55.599980
#> Warning: non-integer x = 62.599980
#> Warning: non-integer x = 72.099980
#> Warning: non-integer x = 77.500000
#> Warning: non-integer x = 83.799990
#> Warning: non-integer x = 93.399960
#> Warning: non-integer x = 98.399960
#> Warning: non-integer x = 24.899990
#> Warning: non-integer x = 29.700000
#> Warning: non-integer x = 33.799990
#> Warning: non-integer x = 38.799990
#> Warning: non-integer x = 48.299990
#> Warning: non-integer x = 52.899990
#> Warning: non-integer x = 70.399960
#> Warning: non-integer x = 76.799990
#> Warning: non-integer x = 83.299990
#> Warning: non-integer x = 102.699950
#> Warning: non-integer x = 24.599990
#> Warning: non-integer x = 34.599980
#> Warning: non-integer x = 41.500000
#> Warning: non-integer x = 49.699980
#> Warning: non-integer x = 56.500000
#> Warning: non-integer x = 66.899960
#> Warning: non-integer x = 75.799990
#> Warning: non-integer x = 84.599980
#> Warning: non-integer x = 98.799990
#> Warning: non-integer x = 107.199950
#> Warning: non-integer x = 32.699980
#> Warning: non-integer x = 38.299990
#> Warning: non-integer x = 45.500000
#> Warning: non-integer x = 48.599980
#> Warning: non-integer x = 62.599980
#> Warning: non-integer x = 70.199950
#> Warning: non-integer x = 73.799990
#> Warning: non-integer x = 80.299990
#> Warning: non-integer x = 85.500000
#> Warning: non-integer x = 93.099980
#> Warning: non-integer x = 35.199980
#> Warning: non-integer x = 41.500000
#> Warning: non-integer x = 47.699980
#> Warning: non-integer x = 53.799990
#> Warning: non-integer x = 61.599980
#> Warning: non-integer x = 79.099980
#> Warning: non-integer x = 85.699950
#> Warning: non-integer x = 92.500000
#> Warning: non-integer x = 101.099980
#> Warning: non-integer x = 34.500000
#> Warning: non-integer x = 40.399990
#> Warning: non-integer x = 53.500000
#> Warning: non-integer x = 60.199980
#> Warning: non-integer x = 68.599980
#> Warning: non-integer x = 73.399960
#> Warning: non-integer x = 78.500000
#> Warning: non-integer x = 90.599980
#> Warning: non-integer x = 97.199950
#> Warning: non-integer x = 31.200000
#> Warning: non-integer x = 38.799990
#> Warning: non-integer x = 44.699980
#> Warning: non-integer x = 53.599980
#> Warning: non-integer x = 59.199980
#> Warning: non-integer x = 67.099980
#> Warning: non-integer x = 76.399960
#> Warning: non-integer x = 83.799990
#> Warning: non-integer x = 71.599980
#> Warning: non-integer x = 97.399960
#> Warning: non-integer x = 107.099980
#> Warning: non-integer x = 24.299990
#> Warning: non-integer x = 28.399990
#> Warning: non-integer x = 40.299990
#> Warning: non-integer x = 46.199980
#> Warning: non-integer x = 51.899990
#> Warning: non-integer x = 58.099980
#> Warning: non-integer x = 62.599980
#> Warning: non-integer x = 68.199950
#> Warning: non-integer x = 78.399960
#> Warning: non-integer x = 87.500000
#> Warning: non-integer x = 90.399960
#> Warning: non-integer x = 24.599990
#> Warning: non-integer x = 25.799990
#> Warning: non-integer x = 28.599990
#> Warning: non-integer x = 30.399990
#> Warning: non-integer x = 36.099980
#> Warning: non-integer x = 36.099980
#> Warning: non-integer x = 49.299990
#> Warning: non-integer x = 57.899990
#> Warning: non-integer x = 75.699950
#> Warning: non-integer x = 81.399960
#> Warning: non-integer x = 88.699950
#> Warning: non-integer x = 21.200000
#> Warning: non-integer x = 25.799990
#> Warning: non-integer x = 29.399990
#> Warning: non-integer x = 36.799990
#> Warning: non-integer x = 43.199980
#> Warning: non-integer x = 51.099980
#> Warning: non-integer x = 59.299990
#> Warning: non-integer x = 72.799990
#> Warning: non-integer x = 78.599980
#> Warning: non-integer x = 93.799990
#> Warning: non-integer x = 95.199950
#> Warning: non-integer x = 29.399990
#> Warning: non-integer x = 34.599980
#> Warning: non-integer x = 40.599980
#> Warning: non-integer x = 46.599980
#> Warning: non-integer x = 62.500000
#> Warning: non-integer x = 78.299990
#> Warning: non-integer x = 83.199950
#> Warning: non-integer x = 34.699980
#> Warning: non-integer x = 42.799990
#> Warning: non-integer x = 48.399990
#> Warning: non-integer x = 53.299990
#> Warning: non-integer x = 71.099980
#> Warning: non-integer x = 77.399960
#> Warning: non-integer x = 94.500000
#> Warning: non-integer x = 104.599980
#> Warning: non-integer x = 29.599990
#> Warning: non-integer x = 34.199980
#> Warning: non-integer x = 63.500000
#> Warning: non-integer x = 72.699950
#> Warning: non-integer x = 78.699950
#> Warning: non-integer x = 84.899960
#> Warning: non-integer x = 96.399960
#> Warning: non-integer x = 100.199950
#> Warning: non-integer x = 22.399990
#> Warning: non-integer x = 25.200000
#> Warning: non-integer x = 36.500000
#> Warning: non-integer x = 43.299990
#> Warning: non-integer x = 58.099980
#> Warning: non-integer x = 70.699950
#> Warning: non-integer x = 78.099980
#> Warning: non-integer x = 83.500000
#> Warning: non-integer x = 26.599990
#> Warning: non-integer x = 28.200000
#> Warning: non-integer x = 38.399990
#> Warning: non-integer x = 45.799990
#> Warning: non-integer x = 54.899990
#> Warning: non-integer x = 60.099980
#> Warning: non-integer x = 56.899990
#> Warning: non-integer x = 64.899960
#> Warning: non-integer x = 81.899960
#> Warning: non-integer x = 90.899960
#> Warning: non-integer x = 28.200000
#> Warning: non-integer x = 35.299990
#> Warning: non-integer x = 42.699980
#> Warning: non-integer x = 48.199980
#> Warning: non-integer x = 53.199980
#> Warning: non-integer x = 64.799990
#> Warning: non-integer x = 71.500000
#> Warning: non-integer x = 76.099980
#> Warning: non-integer x = 86.099980
#> Warning: non-integer x = 92.500000
#> Warning: non-integer x = 27.099990
#> Warning: non-integer x = 34.699980
#> Warning: non-integer x = 40.699980
#> Warning: non-integer x = 47.599980
#> Warning: non-integer x = 56.500000
#> Warning: non-integer x = 72.299990
#> Warning: non-integer x = 95.099980
#> Warning: non-integer x = 105.699950
#> Warning: non-integer x = 37.799990
#> Warning: non-integer x = 43.299990
#> Warning: non-integer x = 51.500000
#> Warning: non-integer x = 63.699980
#> Warning: non-integer x = 70.500000
#> Warning: non-integer x = 80.399960
#> Warning: non-integer x = 106.099980
#> Warning: non-integer x = 29.500000
#> Warning: non-integer x = 31.399990
#> Warning: non-integer x = 39.399990
#> Warning: non-integer x = 44.799990
#> Warning: non-integer x = 50.299990
#> Warning: non-integer x = 65.599980
#> Warning: non-integer x = 72.399960
#> Warning: non-integer x = 80.199950
#> Warning: non-integer x = 90.399960
#> Warning: non-integer x = 99.099980
#> Warning: non-integer x = 22.799990
#> Warning: non-integer x = 33.899990
#> Warning: non-integer x = 37.699980
#> Warning: non-integer x = 41.399990
#> Warning: non-integer x = 48.199980
#> Warning: non-integer x = 53.699980
#> Warning: non-integer x = 62.599980
#> Warning: non-integer x = 70.399960
#> Warning: non-integer x = 78.500000
#> Warning: non-integer x = 85.899960
#> Warning: non-integer x = 22.299990
#> Warning: non-integer x = 24.700000
#> Warning: non-integer x = 30.399990
#> Warning: non-integer x = 36.899990
#> Warning: non-integer x = 44.099980
#> Warning: non-integer x = 52.799990
#> Warning: non-integer x = 60.299990
#> Warning: non-integer x = 67.799990
#> Warning: non-integer x = 75.399960
#> Warning: non-integer x = 87.599980
#> Warning: non-integer x = 101.399960
#> Warning: non-integer x = 23.799990
#> Warning: non-integer x = 37.099980
#> Warning: non-integer x = 43.199980
#> Warning: non-integer x = 50.399990
#> Warning: non-integer x = 60.199980
#> Warning: non-integer x = 82.399960
#> Warning: non-integer x = 91.199950
#> Warning: non-integer x = 101.799990
#> Warning: non-integer x = 109.599980
#> Warning: non-integer x = 25.799990
#> Warning: non-integer x = 29.200000
#> Warning: non-integer x = 35.799990
#> Warning: non-integer x = 41.099980
#> Warning: non-integer x = 55.199980
#> Warning: non-integer x = 62.199980
#> Warning: non-integer x = 69.799990
#> Warning: non-integer x = 76.399960
#> Warning: non-integer x = 83.199950
#> Warning: non-integer x = 93.199950
#> Warning: non-integer x = 32.899990
#> Warning: non-integer x = 46.199980
#> Warning: non-integer x = 53.500000
#> Warning: non-integer x = 58.599980
#> Warning: non-integer x = 68.199950
#> Warning: non-integer x = 80.799990
#> Warning: non-integer x = 84.500000
#> Warning: non-integer x = 24.799990
#> Warning: non-integer x = 35.799990
#> Warning: non-integer x = 42.299990
#> Warning: non-integer x = 50.899990
#> Warning: non-integer x = 57.799990
#> Warning: non-integer x = 64.799990
#> Warning: non-integer x = 71.799990
#> Warning: non-integer x = 83.799990
#> Warning: non-integer x = 92.099980
#> Warning: non-integer x = 26.200000
#> Warning: non-integer x = 30.399990
#> Warning: non-integer x = 37.399990
#> Warning: non-integer x = 49.500000
#> Warning: non-integer x = 55.099980
#> Warning: non-integer x = 64.199950
#> Warning: non-integer x = 75.599980
#> Warning: non-integer x = 83.299990
#> Warning: non-integer x = 87.399960
#> Warning: non-integer x = 94.699950
#> Warning: non-integer x = 98.399960
#> Warning: non-integer x = 32.500000
#> Warning: non-integer x = 38.799990
#> Warning: non-integer x = 47.799990
#> Warning: non-integer x = 54.500000
#> Warning: non-integer x = 70.099980
#> Warning: non-integer x = 76.599980
#> Warning: non-integer x = 86.199950
#> Warning: non-integer x = 95.799990
#> Warning: non-integer x = 101.399960
#> Warning: non-integer x = 112.299990
#> Warning: non-integer x = 32.399990
#> Warning: non-integer x = 38.199980
#> Warning: non-integer x = 44.699980
#> Warning: non-integer x = 53.199980
#> Warning: non-integer x = 62.399990
#> Warning: non-integer x = 70.500000
#> Warning: non-integer x = 72.799990
#> Warning: non-integer x = 86.599980
#> Warning: non-integer x = 93.299990
#> Warning: non-integer x = 104.199950
#> Warning: non-integer x = 27.399990
#> Warning: non-integer x = 32.199980
#> Warning: non-integer x = 37.899990
#> Warning: non-integer x = 45.899990
#> Warning: non-integer x = 51.199980
#> Warning: non-integer x = 57.099980
#> Warning: non-integer x = 76.799990
#> Warning: non-integer x = 83.099980
#> Warning: non-integer x = 91.599980
#> Warning: non-integer x = 100.500000
#> Warning: non-integer x = 105.299990
#> Warning: non-integer x = 27.099990
#> Warning: non-integer x = 30.299990
#> Warning: non-integer x = 38.199980
#> Warning: non-integer x = 43.199980
#> Warning: non-integer x = 52.500000
#> Warning: non-integer x = 67.299990
#> Warning: non-integer x = 74.199950
#> Warning: non-integer x = 81.199950
#> Warning: non-integer x = 87.500000
#> Warning: non-integer x = 94.399960
#> Warning: non-integer x = 101.399960
#> Warning: non-integer x = 26.799990
#> Warning: non-integer x = 31.099990
#> Warning: non-integer x = 37.799990
#> Warning: non-integer x = 44.599980
#> Warning: non-integer x = 50.500000
#> Warning: non-integer x = 58.299990
#> Warning: non-integer x = 68.299990
#> Warning: non-integer x = 75.099980
#> Warning: non-integer x = 84.399960
#> Warning: non-integer x = 87.299990
#> Warning: non-integer x = 96.399960
#> Warning: non-integer x = 104.399960
#> Warning: non-integer x = 24.500000
#> Warning: non-integer x = 31.599990
#> Warning: non-integer x = 37.399990
#> Warning: non-integer x = 46.399990
#> Warning: non-integer x = 53.599980
#> Warning: non-integer x = 61.199980
#> Warning: non-integer x = 70.099980
#> Warning: non-integer x = 77.500000
#> Warning: non-integer x = 86.199950
#> Warning: non-integer x = 89.199950
#> Warning: non-integer x = 96.899960
#> Warning: non-integer x = 104.500000
#> Warning: non-integer x = 23.099990
#> Warning: non-integer x = 30.099990
#> Warning: non-integer x = 40.699980
#> Warning: non-integer x = 46.599980
#> Warning: non-integer x = 52.099980
#> Warning: non-integer x = 60.500000
#> Warning: non-integer x = 73.500000
#> Warning: non-integer x = 77.299990
#> Warning: non-integer x = 86.899960
#> Warning: non-integer x = 96.399960
#> Warning: non-integer x = 21.500000
#> Warning: non-integer x = 24.599990
#> Warning: non-integer x = 28.500000
#> Warning: non-integer x = 39.099980
#> Warning: non-integer x = 45.299990
#> Warning: non-integer x = 59.599980
#> Warning: non-integer x = 69.500000
#> Warning: non-integer x = 75.199950
#> Warning: non-integer x = 84.899960
#> Warning: non-integer x = 93.399960
#> Warning: non-integer x = 99.299990
#> Warning: non-integer x = 24.099990
#> Warning: non-integer x = 28.099990
#> Warning: non-integer x = 32.500000
#> Warning: non-integer x = 40.599980
#> Warning: non-integer x = 53.799990
#> Warning: non-integer x = 61.199980
#> Warning: non-integer x = 71.899960
#> Warning: non-integer x = 79.500000
#> Warning: non-integer x = 87.399960
#> Warning: non-integer x = 93.299990
#> Warning: non-integer x = 101.599980
#> Warning: non-integer x = 34.199980
#> Warning: non-integer x = 42.099980
#> Warning: non-integer x = 50.099980
#> Warning: non-integer x = 57.599980
#> Warning: non-integer x = 62.299990
#> Warning: non-integer x = 76.099980
#> Warning: non-integer x = 81.500000
#> Warning: non-integer x = 84.299990
#> Warning: non-integer x = 97.899960
#> Warning: non-integer x = 108.599980
#> Warning: non-integer x = 25.099990
#> Warning: non-integer x = 23.200000
#> Warning: non-integer x = 28.899990
#> Warning: non-integer x = 36.799990
#> Warning: non-integer x = 44.199980
#> Warning: non-integer x = 52.799990
#> Warning: non-integer x = 62.199980
#> Warning: non-integer x = 69.399960
#> Warning: non-integer x = 78.099980
#> Warning: non-integer x = 88.799990
#> Warning: non-integer x = 102.399960
#> Warning: non-integer x = 106.799990
#> Warning: non-integer x = 32.199980
#> Warning: non-integer x = 36.399990
#> Warning: non-integer x = 41.799990
#> Warning: non-integer x = 48.199980
#> Warning: non-integer x = 55.899990
#> Warning: non-integer x = 60.099980
#> Warning: non-integer x = 68.500000
#> Warning: non-integer x = 77.799990
#> Warning: non-integer x = 86.399960
#> Warning: non-integer x = 96.699950
#> Warning: non-integer x = 111.099980
#> Warning: non-integer x = 115.399960
#> Warning: non-integer x = 24.700000
#> Warning: non-integer x = 28.899990
#> Warning: non-integer x = 34.599980
#> Warning: non-integer x = 42.199980
#> Warning: non-integer x = 48.199980
#> Warning: non-integer x = 54.399990
#> Warning: non-integer x = 62.500000
#> Warning: non-integer x = 69.799990
#> Warning: non-integer x = 78.599980
#> Warning: non-integer x = 96.399960
#> Warning: non-integer x = 103.199950
#> Warning: non-integer x = 24.200000
#> Warning: non-integer x = 29.200000
#> Warning: non-integer x = 36.199980
#> Warning: non-integer x = 43.099980
#> Warning: non-integer x = 49.199980
#> Warning: non-integer x = 70.899960
#> Warning: non-integer x = 87.199950
#> Warning: non-integer x = 94.699950
#> Warning: non-integer x = 100.599980
#> Warning: non-integer x = 24.500000
#> Warning: non-integer x = 28.399990
#> Warning: non-integer x = 41.899990
#> Warning: non-integer x = 55.599980
#> Warning: non-integer x = 64.099980
#> Warning: non-integer x = 71.199950
#> Warning: non-integer x = 80.500000
#> Warning: non-integer x = 90.199950
#> Warning: non-integer x = 98.899960
#> Warning: non-integer x = 103.799990
#> Warning: non-integer x = 26.599990
#> Warning: non-integer x = 33.099980
#> Warning: non-integer x = 38.399990
#> Warning: non-integer x = 46.599980
#> Warning: non-integer x = 52.299990
#> Warning: non-integer x = 61.500000
#> Warning: non-integer x = 69.299990
#> Warning: non-integer x = 74.199950
#> Warning: non-integer x = 80.199950
#> Warning: non-integer x = 87.199950
#> Warning: non-integer x = 94.799990
#> Warning: non-integer x = 26.799990
#> Warning: non-integer x = 33.599980
#> Warning: non-integer x = 42.099980
#> Warning: non-integer x = 47.899990
#> Warning: non-integer x = 56.199980
#> Warning: non-integer x = 63.399990
#> Warning: non-integer x = 71.799990
#> Warning: non-integer x = 87.599980
#> Warning: non-integer x = 94.899960
#> Warning: non-integer x = 30.599990
#> Warning: non-integer x = 37.599980
#> Warning: non-integer x = 42.799990
#> Warning: non-integer x = 47.399990
#> Warning: non-integer x = 56.599980
#> Warning: non-integer x = 62.399990
#> Warning: non-integer x = 71.799990
#> Warning: non-integer x = 86.799990
#> Warning: non-integer x = 94.199950
#> Warning: non-integer x = 101.899960
#> Warning: non-integer x = 112.500000
#> Warning: non-integer x = 21.899990
#> Warning: non-integer x = 26.799990
#> Warning: non-integer x = 34.699980
#> Warning: non-integer x = 41.299990
#> Warning: non-integer x = 48.699980
#> Warning: non-integer x = 57.899990
#> Warning: non-integer x = 65.399960
#> Warning: non-integer x = 71.299990
#> Warning: non-integer x = 80.599980
#> Warning: non-integer x = 88.399960
#> Warning: non-integer x = 96.299990
#> Warning: non-integer x = 103.500000
#> Warning: non-integer x = 26.799990
#> Warning: non-integer x = 39.099980
#> Warning: non-integer x = 46.399990
#> Warning: non-integer x = 61.599980
#> Warning: non-integer x = 68.799990
#> Warning: non-integer x = 72.199950
#> Warning: non-integer x = 77.199950
#> Warning: non-integer x = 83.799990
#> Warning: non-integer x = 90.799990
#> Warning: non-integer x = 24.200000
#> Warning: non-integer x = 25.399990
#> Warning: non-integer x = 30.700000
#> Warning: non-integer x = 35.599980
#> Warning: non-integer x = 42.799990
#> Warning: non-integer x = 49.500000
#> Warning: non-integer x = 57.199980
#> Warning: non-integer x = 65.599980
#> Warning: non-integer x = 72.399960
#> Warning: non-integer x = 81.199950
#> Warning: non-integer x = 85.199950
#> Warning: non-integer x = 90.199950
#> Warning: non-integer x = 25.200000
#> Warning: non-integer x = 29.799990
#> Warning: non-integer x = 35.699980
#> Warning: non-integer x = 55.899990
#> Warning: non-integer x = 63.500000
#> Warning: non-integer x = 79.500000
#> Warning: non-integer x = 83.599980
#> Warning: non-integer x = 89.399960
#> Warning: non-integer x = 98.500000
#> Warning: non-integer x = 23.500000
#> Warning: non-integer x = 25.899990
#> Warning: non-integer x = 32.199980
#> Warning: non-integer x = 32.399990
#> Warning: non-integer x = 29.799990
#> Warning: non-integer x = 32.799990
#> Warning: non-integer x = 56.899990
#> Warning: non-integer x = 64.799990
#> Warning: non-integer x = 72.099980
#> Warning: non-integer x = 81.399960
#> Warning: non-integer x = 26.599990
#> Warning: non-integer x = 38.199980
#> Warning: non-integer x = 46.399990
#> Warning: non-integer x = 53.199980
#> Warning: non-integer x = 59.599980
#> Warning: non-integer x = 66.199950
#> Warning: non-integer x = 72.599980
#> Warning: non-integer x = 80.099980
#> Warning: non-integer x = 97.699950
#> Warning: non-integer x = 35.399990
#> Warning: non-integer x = 42.299990
#> Warning: non-integer x = 49.500000
#> Warning: non-integer x = 58.699980
#> Warning: non-integer x = 66.599980
#> Warning: non-integer x = 71.399960
#> Warning: non-integer x = 78.799990
#> Warning: non-integer x = 85.500000
#> Warning: non-integer x = 91.299990
#> Warning: non-integer x = 97.799990
#> Warning: non-integer x = 29.299990
#> Warning: non-integer x = 45.399990
#> Warning: non-integer x = 51.799990
#> Warning: non-integer x = 58.799990
#> Warning: non-integer x = 62.599980
#> Warning: non-integer x = 71.599980
#> Warning: non-integer x = 76.699950
#> Warning: non-integer x = 81.699950
#> Warning: non-integer x = 86.399960
#> Warning: non-integer x = 92.099980
#> Warning: non-integer x = 25.299990
#> Warning: non-integer x = 31.599990
#> Warning: non-integer x = 35.199980
#> Warning: non-integer x = 37.299990
#> Warning: non-integer x = 40.799990
#> Warning: non-integer x = 47.699980
#> Warning: non-integer x = 64.599980
#> Warning: non-integer x = 71.199950
#> Warning: non-integer x = 80.500000
#> Warning: non-integer x = 85.899960
#> Warning: non-integer x = 22.099990
#> Warning: non-integer x = 26.299990
#> Warning: non-integer x = 33.299990
#> Warning: non-integer x = 39.599980
#> Warning: non-integer x = 41.399990
#> Warning: non-integer x = 46.299990
#> Warning: non-integer x = 52.399990
#> Warning: non-integer x = 59.399990
#> Warning: non-integer x = 66.699950
#> Warning: non-integer x = 73.899960
#> Warning: non-integer x = 76.799990
#> Warning: non-integer x = 23.700000
#> Warning: non-integer x = 28.799990
#> Warning: non-integer x = 35.099980
#> Warning: non-integer x = 49.099980
#> Warning: non-integer x = 58.799990
#> Warning: non-integer x = 65.399960
#> Warning: non-integer x = 68.500000
#> Warning: non-integer x = 90.699950
#> Warning: non-integer x = 28.899990
#> Warning: non-integer x = 35.799990
#> Warning: non-integer x = 47.500000
#> Warning: non-integer x = 56.399990
#> Warning: non-integer x = 63.699980
#> Warning: non-integer x = 70.299990
#> Warning: non-integer x = 89.199950
#> Warning: non-integer x = 96.299990
#> Warning: non-integer x = 99.799990
#> Warning: non-integer x = 30.099990
#> Warning: non-integer x = 37.099980
#> Warning: non-integer x = 53.599980
#> Warning: non-integer x = 61.099980
#> Warning: non-integer x = 78.199950
#> Warning: non-integer x = 87.099980
#> Warning: non-integer x = 94.399960
#> Warning: non-integer x = 102.500000
#> Warning: non-integer x = 107.500000
#> Warning: non-integer x = 35.399990
#> Warning: non-integer x = 48.599980
#> Warning: non-integer x = 56.399990
#> Warning: non-integer x = 75.500000
#> Warning: non-integer x = 77.799990
#> Warning: non-integer x = 85.799990
#> Warning: non-integer x = 101.399960
#> Warning: non-integer x = 27.299990
#> Warning: non-integer x = 31.200000
#> Warning: non-integer x = 43.899990
#> Warning: non-integer x = 49.799990
#> Warning: non-integer x = 57.199980
#> Warning: non-integer x = 65.199950
#> Warning: non-integer x = 73.199950
#> Warning: non-integer x = 78.799990
#> Warning: non-integer x = 87.399960
#> Warning: non-integer x = 100.500000
#> Warning: non-integer x = 25.700000
#> Warning: non-integer x = 28.700000
#> Warning: non-integer x = 33.399990
#> Warning: non-integer x = 46.699980
#> Warning: non-integer x = 56.599980
#> Warning: non-integer x = 65.199950
#> Warning: non-integer x = 73.199950
#> Warning: non-integer x = 81.699950
#> Warning: non-integer x = 90.299990
#> Warning: non-integer x = 103.500000
#> Alternative hypothesis: two-sided 
#>              Coef.  Estimate      SE Null value  t-stat d.f. (Satt)
#>        (Intercept) 21.857868 0.70806          0 30.8702        22.0
#>            CuCu035  0.526734 0.96051          0  0.5484        45.7
#>            CuCu175  0.042784 1.04101          0  0.0411        44.9
#>               Time  2.885157 0.36697          0  7.8622        22.0
#>          I(Time^2)  0.614042 0.07137          0  8.6035        21.9
#>          I(Time^3) -0.026295 0.00389          0 -6.7577        21.9
#>       CuCu035:Time -0.405424 0.64915          0 -0.6245        45.6
#>       CuCu175:Time  0.856941 0.62800          0  1.3645        45.0
#>  CuCu035:I(Time^2)  0.018286 0.11746          0  0.1557        45.5
#>  CuCu175:I(Time^2) -0.096071 0.11301          0 -0.8501        44.9
#>  CuCu035:I(Time^3)  0.000547 0.00610          0  0.0897        45.3
#>  CuCu175:I(Time^3)  0.002700 0.00577          0  0.4682        44.7
#>  p-val (Satt) Sig.
#>        <0.001  ***
#>         0.586     
#>         0.967     
#>        <0.001  ***
#>        <0.001  ***
#>        <0.001  ***
#>         0.535     
#>         0.179     
#>         0.877     
#>         0.400     
#>         0.929     
#>         0.642     
```
