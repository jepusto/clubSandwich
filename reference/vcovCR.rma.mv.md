# Cluster-robust variance-covariance matrix for a rma.mv object.

`vcovCR` returns a sandwich estimate of the variance-covariance matrix
of a set of regression coefficient estimates from a
[`rma.mv`](https://wviechtb.github.io/metafor/reference/rma.mv.html)
object.

## Usage

``` r
# S3 method for class 'rma.mv'
vcovCR(obj, cluster, type, target, inverse_var, form = "sandwich", ...)
```

## Arguments

- obj:

  Fitted model for which to calculate the variance-covariance matrix

- cluster:

  Optional expression or vector indicating which observations belong to
  the same cluster. If not specified, will be set to the factor in the
  random-effects structure with the fewest distinct levels. Caveat
  emptor: the function does not check that the random effects are
  nested.

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
  variance-covariance structure of the `rma.mv` object.

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

mfor_fit <- rma.mv(yi ~ year + deltype, 
                 V = vi, random = ~ 1 | study / esid,
                 data = dat.assink2016)
mfor_fit

mfor_CR2 <- vcovCR(mfor_fit, type = "CR2")
mfor_CR2

coef_test(mfor_fit, vcov = mfor_CR2, test = c("Satterthwaite", "saddlepoint"))
Wald_test(mfor_fit, constraints = constrain_zero(3:4), vcov = mfor_CR2)

})
#> > library(metafor)
#> > data(dat.assink2016, package = "metadat")
#> > mfor_fit <- rma.mv(yi ~ year + deltype, V = vi, random = ~1 | study/esid, 
#> +     data = dat.assink2016)
#> > mfor_fit
#> 
#> Multivariate Meta-Analysis Model (k = 100; method: REML)
#> 
#> Variance Components:
#> 
#>             estim    sqrt  nlvls  fixed      factor 
#> sigma^2.1  0.1493  0.3863     17     no       study 
#> sigma^2.2  0.0853  0.2920    100     no  study/esid 
#> 
#> Test for Residual Heterogeneity:
#> QE(df = 96) = 610.2644, p-val < .0001
#> 
#> Test of Moderators (coefficients 2:4):
#> QM(df = 3) = 19.2399, p-val = 0.0002
#> 
#> Model Results:
#> 
#>                 estimate      se     zval    pval    ci.lb    ci.ub      
#> intrcpt          -0.2438  0.2101  -1.1605  0.2458  -0.6556   0.1680      
#> year             -0.0380  0.0183  -2.0773  0.0378  -0.0738  -0.0021    * 
#> deltypegeneral    0.7094  0.1914   3.7069  0.0002   0.3343   1.0845  *** 
#> deltypeovert      0.5054  0.2099   2.4078  0.0160   0.0940   0.9168    * 
#> 
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> > mfor_CR2 <- vcovCR(mfor_fit, type = "CR2")
#> > mfor_CR2
#>                     intrcpt          year deltypegeneral deltypeovert
#> intrcpt         0.015081059 -0.0020131555  -0.0036579995 -0.007130740
#> year           -0.002013156  0.0007603265  -0.0002607379  0.000218784
#> deltypegeneral -0.003658000 -0.0002607379   0.0042972555  0.004769307
#> deltypeovert   -0.007130740  0.0002187840   0.0047693066  0.010075586
#> > coef_test(mfor_fit, vcov = mfor_CR2, test = c("Satterthwaite", "saddlepoint"))
#> Alternative hypothesis: two-sided 
#>           Coef. Estimate     SE Null value t-stat d.f. (Satt) p-val (Satt) Sig.
#>         intrcpt   -0.244 0.1228          0  -1.99        5.35      0.10015     
#>            year   -0.038 0.0276          0  -1.38        7.22      0.20977     
#>  deltypegeneral    0.709 0.0656          0  10.82        2.25      0.00555   **
#>    deltypeovert    0.505 0.1004          0   5.03        2.03      0.03614    *
#>   s.p. p-val (Saddle) Sig.
#>  0.327         0.0915    .
#>  0.209         0.2103     
#>  0.437         <0.001  ***
#>  0.344         0.0278    *
#> > Wald_test(mfor_fit, constraints = constrain_zero(3:4), vcov = mfor_CR2)
#>  test Fstat df_num df_denom p_val sig
#>   HTZ  30.6      2    0.842 0.164    
```
