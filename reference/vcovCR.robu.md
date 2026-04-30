# Cluster-robust variance-covariance matrix for a robu object.

`vcovCR` returns a sandwich estimate of the variance-covariance matrix
of a set of regression coefficient estimates from a
[`robu`](https://rdrr.io/pkg/robumeta/man/robu.html) object.

## Usage

``` r
# S3 method for class 'robu'
vcovCR(obj, cluster, type, target, inverse_var, form = "sandwich", ...)
```

## Arguments

- obj:

  Fitted model for which to calculate the variance-covariance matrix

- cluster:

  Optional expression or vector indicating which observations belong to
  the same cluster. If not specified, will be set to the `studynum` used
  in fitting the [`robu`](https://rdrr.io/pkg/robumeta/man/robu.html)
  object.

- type:

  Character string specifying which small-sample adjustment should be
  used, with available options `"CR0"`, `"CR1"`, `"CR1p"`, `"CR1S"`,
  `"CR2"`, or `"CR3"`. See "Details" section of
  [`vcovCR`](http://jepusto.github.io/clubSandwich/reference/vcovCR.md)
  for further information.

- target:

  Optional matrix or vector describing the working variance-covariance
  model used to calculate the `CR2` and `CR4` adjustment matrices. If
  not specified, the target is taken to be the inverse of the estimated
  weights used in fitting the
  [`robu`](https://rdrr.io/pkg/robumeta/man/robu.html) object.

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

if (requireNamespace("robumeta", quietly = TRUE)) withAutoprint({
library(robumeta)
data(hierdat)

robu_fit <- robu(effectsize ~ binge + followup + sreport + age, 
                 data = hierdat, studynum = studyid, 
                 var.eff.size = var, modelweights = "HIER")
robu_fit

robu_CR2 <- vcovCR(robu_fit, type = "CR2")
robu_CR2
coef_test(robu_fit, vcov = robu_CR2, test = c("Satterthwaite", "saddlepoint"))

Wald_test(robu_fit, constraints = constrain_zero(c(2,4)), vcov = robu_CR2)
Wald_test(robu_fit, constraints = constrain_zero(2:5), vcov = robu_CR2)

})
#> > library(robumeta)
#> > data(hierdat)
#> > robu_fit <- robu(effectsize ~ binge + followup + sreport + age, data = hierdat, 
#> +     studynum = studyid, var.eff.size = var, modelweights = "HIER")
#> > robu_fit
#> RVE: Hierarchical Effects Model with Small-Sample Corrections 
#> 
#> Model: effectsize ~ binge + followup + sreport + age 
#> 
#> Number of clusters = 15 
#> Number of outcomes = 68 (min = 1 , mean = 4.53 , median = 2 , max = 29 )
#> Omega.sq = 0.1086551 
#> Tau.sq = 0.02362071 
#> 
#>                Estimate   StdErr t-value  dfs P(|t|>) 95% CI.L 95% CI.U Sig
#> 1 X.Intercept.  0.39695 0.658006   0.603 3.06  0.5882 -1.67534  2.46924    
#> 2        binge  0.45158 0.101641   4.443 3.59  0.0144  0.15628  0.74689  **
#> 3     followup  0.00133 0.000723   1.842 2.03  0.2048 -0.00173  0.00439    
#> 4      sreport  0.53876 0.143398   3.757 4.36  0.0170  0.15324  0.92428  **
#> 5          age -0.04371 0.037874  -1.154 2.69  0.3404 -0.17235  0.08492    
#> ---
#> Signif. codes: < .01 *** < .05 ** < .10 *
#> ---
#> Note: If df < 4, do not trust the results> robu_CR2 <- vcovCR(robu_fit, type = "CR2")
#> > robu_CR2
#>               X.Intercept.         binge      followup       sreport
#> X.Intercept.  0.4329720570 -8.055322e-03 -2.957783e-04 -7.180868e-04
#> binge        -0.0080553221  1.033098e-02 -1.580973e-05  3.890863e-03
#> followup     -0.0002957783 -1.580973e-05  5.221116e-07 -1.951302e-05
#> sreport      -0.0007180868  3.890863e-03 -1.951302e-05  2.056284e-02
#> age          -0.0245235730  3.757319e-04  1.683702e-05 -8.892270e-04
#>                        age
#> X.Intercept. -2.452357e-02
#> binge         3.757319e-04
#> followup      1.683702e-05
#> sreport      -8.892270e-04
#> age           1.434475e-03
#> > coef_test(robu_fit, vcov = robu_CR2, test = c("Satterthwaite", "saddlepoint"))
#> Alternative hypothesis: two-sided 
#>         Coef. Estimate       SE Null value t-stat d.f. (Satt) p-val (Satt) Sig.
#>  X.Intercept.  0.39695 0.658006          0  0.603        3.06       0.5882     
#>         binge  0.45158 0.101641          0  4.443        3.59       0.0144    *
#>      followup  0.00133 0.000723          0  1.842        2.03       0.2048     
#>       sreport  0.53876 0.143398          0  3.757        4.36       0.0170    *
#>           age -0.04371 0.037874          0 -1.154        2.69       0.3404     
#>     s.p. p-val (Saddle) Sig.
#>  -0.6440        0.59343     
#>   0.4170        0.00424   **
#>   0.2537        0.18900     
#>   0.4112        0.00767   **
#>   0.0919        0.34521     
#> > Wald_test(robu_fit, constraints = constrain_zero(c(2, 4)), vcov = robu_CR2)
#>  test Fstat df_num df_denom  p_val sig
#>   HTZ  10.4      2     3.48 0.0339   *
#> > Wald_test(robu_fit, constraints = constrain_zero(2:5), vcov = robu_CR2)
#>  test Fstat df_num df_denom p_val sig
#>   HTZ 0.176      4   0.0388 0.962    
```
