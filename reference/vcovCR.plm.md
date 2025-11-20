# Cluster-robust variance-covariance matrix for a plm object.

`vcovCR` returns a sandwich estimate of the variance-covariance matrix
of a set of regression coefficient estimates from a
[`plm`](https://rdrr.io/pkg/plm/man/plm.html) object.

## Usage

``` r
# S3 method for class 'plm'
vcovCR(
  obj,
  cluster,
  type,
  target,
  inverse_var,
  form = "sandwich",
  ignore_FE = FALSE,
  ...
)
```

## Arguments

- obj:

  Fitted model for which to calculate the variance-covariance matrix

- cluster:

  Optional character string, expression, or vector indicating which
  observations belong to the same cluster. For fixed-effect models that
  include individual effects or time effects (but not both), the cluster
  will be taken equal to the included fixed effects if not otherwise
  specified. Clustering on individuals can also be obtained by
  specifying the name of the individual index (e.g.,
  `cluster = "state"`) or `cluster = "individual"`; clustering on time
  periods can be obtained by specifying the name of the time index
  (e.g., `cluster = "year"`) or `cluster = "time"`; if a group index is
  specified, clustering on groups (in which individuals are nested) can
  be obtained by specifying the name of the group index or
  `cluster = "group"`. For random-effects models, the cluster will be
  taken equal to the included random effect identifier if not otherwise
  specified.

- type:

  Character string specifying which small-sample adjustment should be
  used, with available options `"CR0"`, `"CR1"`, `"CR1p"`, `"CR1S"`,
  `"CR2"`, or `"CR3"`. See "Details" section of
  [`vcovCR`](http://jepusto.github.io/clubSandwich/reference/vcovCR.md)
  for further information.

- target:

  Optional matrix or vector describing the working variance-covariance
  model used to calculate the `CR2` and `CR4` adjustment matrices. By
  default, the target is taken to be an identity matrix for fixed effect
  models or the estimated compound-symmetric covariance matrix for
  random effects models.

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

- ignore_FE:

  Optional logical controlling whether fixed effects are ignored when
  calculating small-sample adjustments in models where fixed effects are
  estimated through absorption.

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
if (requireNamespace("plm", quietly = TRUE)) withAutoprint({

  library(plm)
  # fixed effects
  data("Produc", package = "plm")
  plm_FE <- plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp,
                data = Produc, index = c("state","year","region"),
                effect = "individual", model = "within")
  vcovCR(plm_FE, type="CR2")
  vcovCR(plm_FE, type = "CR2", cluster = Produc$region) # clustering on region
  
  # random effects
  plm_RE <- update(plm_FE, model = "random")
  vcovCR(plm_RE, type = "CR2")
  vcovCR(plm_RE, type = "CR2", cluster = Produc$region) # clustering on region
  
  # nested random effects
  plm_nested <- update(plm_FE, effect = "nested", model = "random")
  vcovCR(plm_nested, type = "CR2") # clustering on region
})
#> > library(plm)
#> > data("Produc", package = "plm")
#> > plm_FE <- plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, data = Produc, 
#> +     index = c("state", "year", "region"), effect = "individual", model = "within")
#> > vcovCR(plm_FE, type = "CR2")
#>               log(pcap)       log(pc)      log(emp)         unemp
#> log(pcap)  3.900840e-03 -7.382313e-04 -0.0024837779 -7.825846e-05
#> log(pc)   -7.382313e-04  4.177767e-03 -0.0040101610 -9.960967e-05
#> log(emp)  -2.483778e-03 -4.010161e-03  0.0073152051  1.736024e-04
#> unemp     -7.825846e-05 -9.960967e-05  0.0001736024  6.745842e-06
#> > vcovCR(plm_FE, type = "CR2", cluster = Produc$region)
#>               log(pcap)       log(pc)      log(emp)         unemp
#> log(pcap)  6.973505e-03 -0.0022279395 -0.0044833594 -6.786254e-05
#> log(pc)   -2.227940e-03  0.0054359010 -0.0043767642 -1.297104e-04
#> log(emp)  -4.483359e-03 -0.0043767642  0.0108260634  2.011267e-04
#> unemp     -6.786254e-05 -0.0001297104  0.0002011267  1.170089e-05
#> > plm_RE <- update(plm_FE, model = "random")
#> > vcovCR(plm_RE, type = "CR2")
#>               (Intercept)     log(pcap)       log(pc)      log(emp)
#> (Intercept)  0.0608153880 -8.772256e-03 -7.015754e-03  0.0136354448
#> log(pcap)   -0.0087722560  3.221991e-03  3.731861e-05 -0.0031502393
#> log(pc)     -0.0070157539  3.731861e-05  2.060851e-03 -0.0020993882
#> log(emp)     0.0136354448 -3.150239e-03 -2.099388e-03  0.0054047466
#> unemp        0.0004200323 -7.947309e-05 -6.416279e-05  0.0001390737
#>                     unemp
#> (Intercept)  4.200323e-04
#> log(pcap)   -7.947309e-05
#> log(pc)     -6.416279e-05
#> log(emp)     1.390737e-04
#> unemp        5.848168e-06
#> > vcovCR(plm_RE, type = "CR2", cluster = Produc$region)
#>               (Intercept)     log(pcap)       log(pc)      log(emp)
#> (Intercept)  0.0565942863 -0.0138696773 -0.0056768267  0.0190768009
#> log(pcap)   -0.0138696773  0.0065345184 -0.0003786724 -0.0063601303
#> log(pc)     -0.0056768267 -0.0003786724  0.0024373929 -0.0022193705
#> log(emp)     0.0190768009 -0.0063601303 -0.0022193705  0.0091711531
#> unemp        0.0004716871 -0.0000926249 -0.0000839537  0.0001715663
#>                     unemp
#> (Intercept)  4.716871e-04
#> log(pcap)   -9.262490e-05
#> log(pc)     -8.395370e-05
#> log(emp)     1.715663e-04
#> unemp        1.057574e-05
#> > plm_nested <- update(plm_FE, effect = "nested", model = "random")
#> > vcovCR(plm_nested, type = "CR2")
#>               (Intercept)     log(pcap)       log(pc)      log(emp)
#> (Intercept)  0.0537315947 -1.210951e-02 -0.0073309652  0.0193853190
#> log(pcap)   -0.0121095110  6.078432e-03 -0.0003141431 -0.0060636185
#> log(pc)     -0.0073309652 -3.141431e-04  0.0027593587 -0.0025387198
#> log(emp)     0.0193853190 -6.063618e-03 -0.0025387198  0.0091748282
#> unemp        0.0005354663 -9.267671e-05 -0.0001000830  0.0001852566
#>                     unemp
#> (Intercept)  5.354663e-04
#> log(pcap)   -9.267671e-05
#> log(pc)     -1.000830e-04
#> log(emp)     1.852566e-04
#> unemp        1.135799e-05

pkgs_available <- requireNamespace("plm", quietly = TRUE) & requireNamespace("AER", quietly = TRUE)

if (pkgs_available) withAutoprint({
  # first differencing
  data(Fatalities, package = "AER")
  Fatalities <- within(Fatalities, {
    frate <- 10000 * fatal / pop
    drinkagec <- cut(drinkage, breaks = 18:22, include.lowest = TRUE, right = FALSE)
    drinkagec <- relevel(drinkagec, ref = 4)
  })

  plm_FD <- plm(frate ~ beertax + drinkagec + miles + unemp + log(income),
                data = Fatalities, index = c("state", "year"),
                model = "fd")
  vcovHC(plm_FD, method="arellano", type = "sss", cluster = "group")
  vcovCR(plm_FD, type = "CR1S")
  vcovCR(plm_FD, type = "CR2")
  
})
#> > data(Fatalities, package = "AER")
#> > Fatalities <- within(Fatalities, {
#> +     frate <- 10000 * fatal/pop
#> +     drinkagec <- cut(drinkage, breaks = 18:22, include.lowest = TRUE, right = FALSE)
#> +     drinkagec <- relevel(drinkagec, ref = 4)
#> + })
#> > plm_FD <- plm(frate ~ beertax + drinkagec + miles + unemp + log(income), 
#> +     data = Fatalities, index = c("state", "year"), model = "fd")
#> > vcovHC(plm_FD, method = "arellano", type = "sss", cluster = "group")
#>                    (Intercept)       beertax drinkagec[18,19) drinkagec[19,20)
#> (Intercept)       4.491252e-04  1.312107e-03     3.177175e-05    -3.224139e-05
#> beertax           1.312107e-03  6.112682e-02    -1.374164e-03    -2.155192e-03
#> drinkagec[18,19)  3.177175e-05 -1.374164e-03     5.978249e-03     2.531421e-03
#> drinkagec[19,20) -3.224139e-05 -2.155192e-03     2.531421e-03     2.130959e-03
#> drinkagec[20,21) -2.531288e-04 -9.279081e-04     1.738413e-03     1.185214e-03
#> miles             4.688676e-09 -1.642830e-07    -8.440044e-08    -5.492005e-08
#> unemp            -7.197958e-05 -3.079273e-04    -4.794408e-05    -3.200981e-05
#> log(income)      -1.853129e-02 -2.489411e-02     3.443431e-03     4.787025e-03
#>                  drinkagec[20,21)         miles         unemp   log(income)
#> (Intercept)         -2.531288e-04  4.688676e-09 -7.197958e-05 -1.853129e-02
#> beertax             -9.279081e-04 -1.642830e-07 -3.079273e-04 -2.489411e-02
#> drinkagec[18,19)     1.738413e-03 -8.440044e-08 -4.794408e-05  3.443431e-03
#> drinkagec[19,20)     1.185214e-03 -5.492005e-08 -3.200981e-05  4.787025e-03
#> drinkagec[20,21)     2.529571e-03 -9.131630e-08  1.900456e-04  1.890722e-02
#> miles               -9.131630e-08  1.570809e-11  1.928974e-09 -6.148657e-07
#> unemp                1.900456e-04  1.928974e-09  1.970444e-04  8.234312e-03
#> log(income)          1.890722e-02 -6.148657e-07  8.234312e-03  1.034150e+00
#> attr(,"cluster")
#> [1] "group"
#> > vcovCR(plm_FD, type = "CR1S")
#>                    (Intercept)       beertax drinkagec[18,19) drinkagec[19,20)
#> (Intercept)       4.491252e-04  1.312107e-03     3.177175e-05    -3.224139e-05
#> beertax           1.312107e-03  6.112682e-02    -1.374164e-03    -2.155192e-03
#> drinkagec[18,19)  3.177175e-05 -1.374164e-03     5.978249e-03     2.531421e-03
#> drinkagec[19,20) -3.224139e-05 -2.155192e-03     2.531421e-03     2.130959e-03
#> drinkagec[20,21) -2.531288e-04 -9.279081e-04     1.738413e-03     1.185214e-03
#> miles             4.688676e-09 -1.642830e-07    -8.440044e-08    -5.492005e-08
#> unemp            -7.197958e-05 -3.079273e-04    -4.794408e-05    -3.200981e-05
#> log(income)      -1.853129e-02 -2.489411e-02     3.443431e-03     4.787025e-03
#>                  drinkagec[20,21)         miles         unemp   log(income)
#> (Intercept)         -2.531288e-04  4.688676e-09 -7.197958e-05 -1.853129e-02
#> beertax             -9.279081e-04 -1.642830e-07 -3.079273e-04 -2.489411e-02
#> drinkagec[18,19)     1.738413e-03 -8.440044e-08 -4.794408e-05  3.443431e-03
#> drinkagec[19,20)     1.185214e-03 -5.492005e-08 -3.200981e-05  4.787025e-03
#> drinkagec[20,21)     2.529571e-03 -9.131630e-08  1.900456e-04  1.890722e-02
#> miles               -9.131630e-08  1.570809e-11  1.928974e-09 -6.148657e-07
#> unemp                1.900456e-04  1.928974e-09  1.970444e-04  8.234312e-03
#> log(income)          1.890722e-02 -6.148657e-07  8.234312e-03  1.034150e+00
#> > vcovCR(plm_FD, type = "CR2")
#>                    (Intercept)       beertax drinkagec[18,19) drinkagec[19,20)
#> (Intercept)       4.926777e-04  0.0014903060     5.194011e-05    -1.917433e-05
#> beertax           1.490306e-03  0.0639227452    -1.545309e-03    -2.356755e-03
#> drinkagec[18,19)  5.194011e-05 -0.0015453092     6.783986e-03     2.627033e-03
#> drinkagec[19,20) -1.917433e-05 -0.0023567549     2.627033e-03     2.182684e-03
#> drinkagec[20,21) -2.585076e-04 -0.0011673142     1.811511e-03     1.193471e-03
#> miles            -1.036421e-07 -0.0000011056    -1.901583e-07    -1.298187e-07
#> unemp            -7.270549e-05 -0.0002748930    -3.606990e-05    -3.739358e-05
#> log(income)      -1.903417e-02 -0.0209995682     3.953875e-03     4.817831e-03
#>                  drinkagec[20,21)         miles         unemp   log(income)
#> (Intercept)         -2.585076e-04 -1.036421e-07 -7.270549e-05 -1.903417e-02
#> beertax             -1.167314e-03 -1.105600e-06 -2.748930e-04 -2.099957e-02
#> drinkagec[18,19)     1.811511e-03 -1.901583e-07 -3.606990e-05  3.953875e-03
#> drinkagec[19,20)     1.193471e-03 -1.298187e-07 -3.739358e-05  4.817831e-03
#> drinkagec[20,21)     2.626823e-03 -1.514746e-07  1.955850e-04  1.969001e-02
#> miles               -1.514746e-07  4.431806e-10 -5.460872e-10 -1.329604e-06
#> unemp                1.955850e-04 -5.460872e-10  2.005515e-04  8.377011e-03
#> log(income)          1.969001e-02 -1.329604e-06  8.377011e-03  1.065829e+00
```
