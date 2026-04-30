# Cluster-robust variance-covariance matrix for an ivreg object.

`vcovCR` returns a sandwich estimate of the variance-covariance matrix
of a set of regression coefficient estimates from an ivreg object fitted
from the [AER](https://CRAN.R-project.org/package=AER) package or the
[ivreg](https://CRAN.R-project.org/package=ivreg) package.

## Usage

``` r
# S3 method for class 'ivreg'
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

  Expression or vector indicating which observations belong to the same
  cluster. Required for `ivreg` objects.

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

  Not used for `ivreg` objects.

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

## Details

For any "ivreg" objects fitted via the
[`ivreg`](https://zeileis.github.io/ivreg/reference/ivreg.html) function
from the [ivreg](https://CRAN.R-project.org/package=ivreg) package, only
traditional 2SLS regression method (method = "OLS") is supported.
clubSandwich currently cannot support robust-regression methods such as
M-estimation (method = "M") or MM-estimation (method = "MM").

## See also

[`vcovCR`](http://jepusto.github.io/clubSandwich/reference/vcovCR.md)

## Examples

``` r

if (requireNamespace("AER", quietly = TRUE)) withAutoprint({

  library(AER)
  data("CigarettesSW")
  Cigs <- within(CigarettesSW, {
    rprice <- price/cpi
    rincome <- income/population/cpi
    tdiff <- (taxs - tax)/cpi
  })

  iv_fit_AER <- AER::ivreg(log(packs) ~ log(rprice) + log(rincome) | 
                  log(rincome) + tdiff + I(tax/cpi), data = Cigs)
  vcovCR(iv_fit_AER, cluster = Cigs$state, type = "CR2")
  coef_test(iv_fit_AER, vcov = "CR2", cluster = Cigs$state)

})
#> > library(AER)
#> Loading required package: car
#> Loading required package: carData
#> 
#> Attaching package: ‘car’
#> The following object is masked from ‘package:metafor’:
#> 
#>     vif
#> Loading required package: lmtest
#> Loading required package: zoo
#> 
#> Attaching package: ‘zoo’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.Date, as.Date.numeric
#> Loading required package: sandwich
#> Loading required package: survival
#> > data("CigarettesSW")
#> > Cigs <- within(CigarettesSW, {
#> +     rprice <- price/cpi
#> +     rincome <- income/population/cpi
#> +     tdiff <- (taxs - tax)/cpi
#> + })
#> > iv_fit_AER <- AER::ivreg(log(packs) ~ log(rprice) + log(rincome) | log(rincome) + 
#> +     tdiff + I(tax/cpi), data = Cigs)
#> > vcovCR(iv_fit_AER, cluster = Cigs$state, type = "CR2")
#>              (Intercept) log(rprice) log(rincome)
#> (Intercept)   0.31782384 -0.08285646   0.02774789
#> log(rprice)  -0.08285646  0.03453204  -0.03009923
#> log(rincome)  0.02774789 -0.03009923   0.04292233
#> > coef_test(iv_fit_AER, vcov = "CR2", cluster = Cigs$state)
#> Alternative hypothesis: two-sided 
#>         Coef. Estimate    SE Null value t-stat d.f. (Satt) p-val (Satt) Sig.
#>   (Intercept)    9.736 0.564          0  17.27        22.0       <0.001  ***
#>   log(rprice)   -1.229 0.186          0  -6.61        21.1       <0.001  ***
#>  log(rincome)    0.257 0.207          0   1.24        23.5        0.227     

pkgs_available <- 
  requireNamespace("AER", quietly = TRUE) & 
  requireNamespace("ivreg", quietly = TRUE)
#> Registered S3 methods overwritten by 'ivreg':
#>   method              from
#>   anova.ivreg         AER 
#>   hatvalues.ivreg     AER 
#>   model.matrix.ivreg  AER 
#>   predict.ivreg       AER 
#>   print.ivreg         AER 
#>   print.summary.ivreg AER 
#>   summary.ivreg       AER 
#>   terms.ivreg         AER 
#>   update.ivreg        AER 
#>   vcov.ivreg          AER 

if (pkgs_available) withAutoprint ({

data("CigarettesSW")
  Cigs <- within(CigarettesSW, {
    rprice <- price/cpi
    rincome <- income/population/cpi
    tdiff <- (taxs - tax)/cpi
  })
iv_fit_ivreg <- ivreg::ivreg(log(packs) ~ log(rprice) + log(rincome) | 
                  log(rincome) + tdiff + I(tax/cpi), data = Cigs)
  vcovCR(iv_fit_ivreg, cluster = Cigs$state, type = "CR2")
  coef_test(iv_fit_ivreg, vcov = "CR2", cluster = Cigs$state)
})
#> > data("CigarettesSW")
#> > Cigs <- within(CigarettesSW, {
#> +     rprice <- price/cpi
#> +     rincome <- income/population/cpi
#> +     tdiff <- (taxs - tax)/cpi
#> + })
#> > iv_fit_ivreg <- ivreg::ivreg(log(packs) ~ log(rprice) + log(rincome) | 
#> +     log(rincome) + tdiff + I(tax/cpi), data = Cigs)
#> > vcovCR(iv_fit_ivreg, cluster = Cigs$state, type = "CR2")
#>              (Intercept) log(rprice) log(rincome)
#> (Intercept)   0.31782384 -0.08285646   0.02774789
#> log(rprice)  -0.08285646  0.03453204  -0.03009923
#> log(rincome)  0.02774789 -0.03009923   0.04292233
#> > coef_test(iv_fit_ivreg, vcov = "CR2", cluster = Cigs$state)
#> Alternative hypothesis: two-sided 
#>         Coef. Estimate    SE Null value t-stat d.f. (Satt) p-val (Satt) Sig.
#>   (Intercept)    9.736 0.564          0  17.27        22.0       <0.001  ***
#>   log(rprice)   -1.229 0.186          0  -6.61        21.1       <0.001  ***
#>  log(rincome)    0.257 0.207          0   1.24        23.5        0.227     
```
