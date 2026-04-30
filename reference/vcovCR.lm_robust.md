# Cluster-robust variance-covariance matrix for an `estimatr::lm_robust` object.

`vcovCR` returns a sandwich estimate of the variance-covariance matrix
of a set of regression coefficient estimates from an
[`lm_robust`](https://declaredesign.org/r/estimatr/reference/lm_robust.html)
object.

## Usage

``` r
# S3 method for class 'lm_robust'
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
  cluster. If not specified, will be detected from the `clusters`
  argument of `obj`.

- type:

  Character string specifying which small-sample adjustment should be
  used, with available options `"CR0"`, `"CR1"`, `"CR1p"`, `"CR1S"`,
  `"CR2"`, or `"CR3"`. If not specified, will be detected from the
  `se_type` argument of `obj`. See "Details" section of
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
ChickWeight$Chick <- factor(ChickWeight$Chick, ordered = FALSE)

if (requireNamespace("estimatr", quietly = TRUE)) withAutoprint({
  library(estimatr)

  lm_fit <- lm_robust(
    weight ~ Time + Diet:Time,
    data = ChickWeight
   )
  vcovCR(lm_fit, cluster = ChickWeight$Chick, type = "CR2")

  lm_fit_clust <- lm_robust(
    weight ~ Time + Diet:Time, data = ChickWeight,
    clusters = Chick
   )
  conf_int(lm_fit_clust, vcov = "CR2")

  # similar model via lm_lin()
  lin_fit_clust <- lm_lin(
    weight ~ Diet, 
    covariates = ~ Time,
    data = ChickWeight,
    clusters = Chick
  )
  conf_int(lin_fit_clust, vcov = "CR2")
  
  lm_fit_fe <- lm_robust(
    weight ~ Time:Diet, data = ChickWeight,
    clusters = Chick,
    fixed_effects = ~ Chick
   )
  vcovCR(lm_fit_fe)
  
  # two-way fixed effects model
  data("MortalityRates")
  MortalityRates <- subset(MortalityRates, cause == "Motor Vehicle")
  MortalityRates$state <- factor(MortalityRates$state)
  MortalityRates$year <- factor(MortalityRates$year)
  MLDA_fit <- lm_robust(
    mrate ~ legal + beertaxa + beerpercap + winepercap + spiritpercap,
    fixed_effects = ~ year + state,
    data = MortalityRates,
    cluster = state
  )
  conf_int(MLDA_fit, vcov = "CR2")

  if (requireNamespace("plm", quietly = TRUE)) withAutoprint({

    data("Produc", package = "plm")
    lm_individual <- lm_robust(
      log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp,
      data = Produc,
      fixed_effects = ~ state,
      cluster = state
     )
    vcovCR(lm_individual, type = "CR2")

  })

})
#> > library(estimatr)
#> > lm_fit <- lm_robust(weight ~ Time + Diet:Time, data = ChickWeight)
#> > vcovCR(lm_fit, cluster = ChickWeight$Chick, type = "CR2")
#>             (Intercept)       Time Time:Diet2 Time:Diet3 Time:Diet4
#> (Intercept)   3.9317840 -0.7851508 -0.2560762 -0.0750498  0.2695861
#> Time         -0.7851508  0.4360726 -0.3314209 -0.3436398 -0.3665016
#> Time:Diet2   -0.2560762 -0.3314209  1.3680862  0.3475327  0.3482061
#> Time:Diet3   -0.0750498 -0.3436398  0.3475327  1.1195632  0.3483146
#> Time:Diet4    0.2695861 -0.3665016  0.3482061  0.3483146  0.5416587
#> > lm_fit_clust <- lm_robust(weight ~ Time + Diet:Time, data = ChickWeight, 
#> +     clusters = Chick)
#> > conf_int(lm_fit_clust, vcov = "CR2")
#>        Coef. Estimate    SE d.f. Lower 95% CI Upper 95% CI
#>  (Intercept)    27.86 1.983 48.7       23.874        31.84
#>         Time     7.05 0.660 31.4        5.703         8.40
#>   Time:Diet2     1.61 1.170 19.0       -0.837         4.06
#>   Time:Diet3     3.74 1.058 19.0        1.524         5.95
#>   Time:Diet4     2.86 0.736 18.4        1.317         4.41
#> > lin_fit_clust <- lm_lin(weight ~ Diet, covariates = ~Time, data = ChickWeight, 
#> +     clusters = Chick)
#> > conf_int(lin_fit_clust, vcov = "CR2")
#>         Coef. Estimate     SE d.f. Lower 95% CI Upper 95% CI
#>   (Intercept)   104.26  5.727 18.0       92.227       116.30
#>         Diet2    16.64 11.257 18.7       -6.940        40.23
#>         Diet3    36.42 10.177 18.7       15.095        57.74
#>         Diet4    30.65  6.957 18.5       16.064        45.23
#>        Time_c     6.84  0.759 18.0        5.247         8.44
#>  Diet2:Time_c     1.77  1.488 18.8       -1.349         4.88
#>  Diet3:Time_c     4.58  1.351 18.8        1.751         7.41
#>  Diet4:Time_c     2.87  1.008 18.3        0.757         4.99
#> > lm_fit_fe <- lm_robust(weight ~ Time:Diet, data = ChickWeight, clusters = Chick, 
#> +     fixed_effects = ~Chick)
#> > vcovCR(lm_fit_fe)
#>            Time:Diet1 Time:Diet2 Time:Diet3 Time:Diet4
#> Time:Diet1  0.5644892   0.000000   0.000000  0.0000000
#> Time:Diet2  0.0000000   1.638116   0.000000  0.0000000
#> Time:Diet3  0.0000000   0.000000   1.249162  0.0000000
#> Time:Diet4  0.0000000   0.000000   0.000000  0.4523152
#> > data("MortalityRates")
#> > MortalityRates <- subset(MortalityRates, cause == "Motor Vehicle")
#> > MortalityRates$state <- factor(MortalityRates$state)
#> > MortalityRates$year <- factor(MortalityRates$year)
#> > MLDA_fit <- lm_robust(mrate ~ legal + beertaxa + beerpercap + winepercap + 
#> +     spiritpercap, fixed_effects = ~year + state, data = MortalityRates, cluster = state)
#> > conf_int(MLDA_fit, vcov = "CR2")
#>         Coef. Estimate    SE  d.f. Lower 95% CI Upper 95% CI
#>         legal    0.224  2.59 40.57        -5.00         5.45
#>      beertaxa   -1.433  5.15  8.35       -13.22        10.35
#>    beerpercap   32.387 12.83 28.57         6.12        58.65
#>    winepercap   11.000 11.75 19.12       -13.58        35.59
#>  spiritpercap   -1.353 11.88  3.95       -34.51        31.80
#> > if (requireNamespace("plm", quietly = TRUE)) withAutoprint({
#> +     data("Produc", package = "plm")
#> +     lm_individual <- lm_robust(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, 
#> +         data = Produc, fixed_effects = ~state, cluster = state)
#> +     vcovCR(lm_individual, type = "CR2")
#> + })
#> > data("Produc", package = "plm")
#> > lm_individual <- lm_robust(log(gsp) ~ log(pcap) + log(pc) + log(emp) + 
#> +     unemp, data = Produc, fixed_effects = ~state, cluster = state)
#> > vcovCR(lm_individual, type = "CR2")
#>               log(pcap)       log(pc)      log(emp)         unemp
#> log(pcap)  3.900840e-03 -7.382313e-04 -0.0024837779 -7.825846e-05
#> log(pc)   -7.382313e-04  4.177767e-03 -0.0040101610 -9.960967e-05
#> log(emp)  -2.483778e-03 -4.010161e-03  0.0073152051  1.736024e-04
#> unemp     -7.825846e-05 -9.960967e-05  0.0001736024  6.745842e-06
```
