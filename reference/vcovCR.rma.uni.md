# Cluster-robust variance-covariance matrix for a rma.uni object.

`vcovCR` returns a sandwich estimate of the variance-covariance matrix
of a set of regression coefficient estimates from a
[`rma.uni`](https://wviechtb.github.io/metafor/reference/rma.uni.html)
object.

## Usage

``` r
# S3 method for class 'rma.uni'
vcovCR(obj, cluster, type, target, inverse_var, form = "sandwich", ...)
```

## Arguments

- obj:

  Fitted model for which to calculate the variance-covariance matrix

- cluster:

  Expression or vector indicating which observations belong to the same
  cluster. Required for `rma.uni` objects.

- type:

  Character string specifying which small-sample adjustment should be
  used, with available options `"CR0"`, `"CR1"`, `"CR1p"`, `"CR1S"`,
  `"CR2"`, or `"CR3"`. See "Details" section of
  [`vcovCR`](http://jepusto.github.io/clubSandwich/reference/vcovCR.md)
  for further information.

- target:

  Optional matrix or vector describing the working variance-covariance
  model used to calculate the `CR2` and `CR4` adjustment matrices. If
  not specified, the target is taken to be diagonal with entries equal
  to the estimated marginal variance of the effect sizes.

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

pkgs_available <- 
  requireNamespace("metafor", quietly = TRUE) & 
  requireNamespace("metadat", quietly = TRUE)
  
if (pkgs_available) withAutoprint({

library(metafor)
data(dat.assink2016, package = "metadat")

mfor_fit <- rma.uni(yi ~ year + deltype, vi = vi, 
                 data = dat.assink2016)
mfor_fit

mfor_CR2 <- vcovCR(mfor_fit, type = "CR2", cluster = dat.assink2016$study)
mfor_CR2
coef_test(mfor_fit, vcov = mfor_CR2, test = c("Satterthwaite", "saddlepoint"))
Wald_test(mfor_fit, constraints = constrain_zero(2:4), vcov = mfor_CR2)

})
#> > library(metafor)
#> > data(dat.assink2016, package = "metadat")
#> > mfor_fit <- rma.uni(yi ~ year + deltype, vi = vi, data = dat.assink2016)
#> > mfor_fit
#> 
#> Mixed-Effects Model (k = 100; tau^2 estimator: REML)
#> 
#> tau^2 (estimated amount of residual heterogeneity):     0.2016 (SE = 0.0373)
#> tau (square root of estimated tau^2 value):             0.4490
#> I^2 (residual heterogeneity / unaccounted variability): 89.68%
#> H^2 (unaccounted variability / sampling variability):   9.69
#> R^2 (amount of heterogeneity accounted for):            47.74%
#> 
#> Test for Residual Heterogeneity:
#> QE(df = 96) = 610.2644, p-val < .0001
#> 
#> Test of Moderators (coefficients 2:4):
#> QM(df = 3) = 66.1189, p-val < .0001
#> 
#> Model Results:
#> 
#>                 estimate      se     zval    pval    ci.lb    ci.ub      
#> intrcpt           0.0197  0.1724   0.1143  0.9090  -0.3181   0.3575      
#> year             -0.0682  0.0098  -6.9476  <.0001  -0.0874  -0.0490  *** 
#> deltypegeneral    0.5499  0.1827   3.0097  0.0026   0.1918   0.9079   ** 
#> deltypeovert      0.5606  0.2226   2.5186  0.0118   0.1244   0.9969    * 
#> 
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> > mfor_CR2 <- vcovCR(mfor_fit, type = "CR2", cluster = dat.assink2016$study)
#> > mfor_CR2
#>                     intrcpt          year deltypegeneral  deltypeovert
#> intrcpt         0.004914493 -0.0017836453   0.0019509336 -0.0022872378
#> year           -0.001783645  0.0007185062  -0.0008482905  0.0004549436
#> deltypegeneral  0.001950934 -0.0008482905   0.0083677591  0.0016321051
#> deltypeovert   -0.002287238  0.0004549436   0.0016321051  0.0060609926
#> > coef_test(mfor_fit, vcov = mfor_CR2, test = c("Satterthwaite", "saddlepoint"))
#> Alternative hypothesis: two-sided 
#>           Coef. Estimate     SE Null value t-stat d.f. (Satt) p-val (Satt) Sig.
#>         intrcpt   0.0197 0.0701          0  0.281        6.44       0.7875     
#>            year  -0.0682 0.0268          0 -2.544        6.00       0.0438    *
#>  deltypegeneral   0.5499 0.0915          0  6.011        8.97       <0.001  ***
#>    deltypeovert   0.5606 0.0779          0  7.201        1.95       0.0200    *
#>    s.p. p-val (Saddle) Sig.
#>  -5.009         0.7774     
#>   0.375         0.0380    *
#>   0.446         <0.001  ***
#>   0.375         0.0088   **
#> > Wald_test(mfor_fit, constraints = constrain_zero(2:4), vcov = mfor_CR2)
#>  test Fstat df_num df_denom p_val sig
#>   HTZ  14.9      3     2.63 0.036   *
```
