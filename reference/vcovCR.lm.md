# Cluster-robust variance-covariance matrix for an lm object.

`vcovCR` returns a sandwich estimate of the variance-covariance matrix
of a set of regression coefficient estimates from an
[`lm`](https://rdrr.io/r/stats/lm.html) object.

## Usage

``` r
# S3 method for class 'lm'
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
  cluster. Required for `lm` objects.

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

data("ChickWeight", package = "datasets")
lm_fit <- lm(weight ~ Time + Diet:Time, data = ChickWeight)
vcovCR(lm_fit, cluster = ChickWeight$Chick, type = "CR2")
#>             (Intercept)       Time Time:Diet2 Time:Diet3 Time:Diet4
#> (Intercept)   3.9317840 -0.7851508 -0.2560762 -0.0750498  0.2695861
#> Time         -0.7851508  0.4360726 -0.3314209 -0.3436398 -0.3665016
#> Time:Diet2   -0.2560762 -0.3314209  1.3680862  0.3475327  0.3482061
#> Time:Diet3   -0.0750498 -0.3436398  0.3475327  1.1195632  0.3483146
#> Time:Diet4    0.2695861 -0.3665016  0.3482061  0.3483146  0.5416587

if (requireNamespace("plm", quietly = TRUE)) withAutoprint({

  data("Produc", package = "plm")
  lm_individual <- lm(log(gsp) ~ 0 + state + log(pcap) + log(pc) + log(emp) + unemp, data = Produc)
  individual_index <- !grepl("state", names(coef(lm_individual)))
  vcovCR(lm_individual, cluster = Produc$state, type = "CR2")[individual_index,individual_index]

  # compare to plm()
  plm_FE <- plm::plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, 
                     data = Produc, index = c("state","year"), 
                     effect = "individual", model = "within")
  vcovCR(plm_FE, type="CR2")
  
})
#> > data("Produc", package = "plm")
#> > lm_individual <- lm(log(gsp) ~ 0 + state + log(pcap) + log(pc) + log(emp) + 
#> +     unemp, data = Produc)
#> > individual_index <- !grepl("state", names(coef(lm_individual)))
#> > vcovCR(lm_individual, cluster = Produc$state, type = "CR2")[individual_index, 
#> +     individual_index]
#>               log(pcap)       log(pc)      log(emp)         unemp
#> log(pcap)  3.900840e-03 -7.382313e-04 -0.0024837779 -7.825846e-05
#> log(pc)   -7.382313e-04  4.177767e-03 -0.0040101610 -9.960967e-05
#> log(emp)  -2.483778e-03 -4.010161e-03  0.0073152051  1.736024e-04
#> unemp     -7.825846e-05 -9.960967e-05  0.0001736024  6.745842e-06
#> > plm_FE <- plm::plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, 
#> +     data = Produc, index = c("state", "year"), effect = "individual", model = "within")
#> > vcovCR(plm_FE, type = "CR2")
#>               log(pcap)       log(pc)      log(emp)         unemp
#> log(pcap)  3.900840e-03 -7.382313e-04 -0.0024837779 -7.825846e-05
#> log(pc)   -7.382313e-04  4.177767e-03 -0.0040101610 -9.960967e-05
#> log(emp)  -2.483778e-03 -4.010161e-03  0.0073152051  1.736024e-04
#> unemp     -7.825846e-05 -9.960967e-05  0.0001736024  6.745842e-06
```
