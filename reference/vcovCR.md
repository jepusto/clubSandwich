# Cluster-robust variance-covariance matrix

This is a generic function, with specific methods defined for
[`lm`](https://rdrr.io/r/stats/lm.html),
[`plm`](https://rdrr.io/pkg/plm/man/plm.html),
[`glm`](https://rdrr.io/r/stats/glm.html),
[`gls`](https://rdrr.io/pkg/nlme/man/gls.html),
[`lme`](https://rdrr.io/pkg/nlme/man/lme.html),
[`robu`](https://rdrr.io/pkg/robumeta/man/robu.html),
[`rma.uni`](https://wviechtb.github.io/metafor/reference/rma.uni.html),
[`rma.mv`](https://wviechtb.github.io/metafor/reference/rma.mv.html),
and
[`lm_robust`](https://declaredesign.org/r/estimatr/reference/lm_robust.html)
objects.

`vcovCR` returns a sandwich estimate of the variance-covariance matrix
of a set of regression coefficient estimates.

## Usage

``` r
vcovCR(obj, cluster, type, target, inverse_var, form, ...)

# Default S3 method
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
  cluster. For some classes, the cluster will be detected automatically
  if not specified.

- type:

  Character string specifying which small-sample adjustment should be
  used, with available options `"CR0"`, `"CR1"`, `"CR1p"`, `"CR1S"`,
  `"CR2"`, or `"CR3"`. See "Details" section of `vcovCR` for further
  information.

- target:

  Optional matrix or vector describing the working variance-covariance
  model used to calculate the `CR2` and `CR4` adjustment matrices. If a
  vector, the target matrix is assumed to be diagonal. If not specified,
  `vcovCR` will attempt to infer a value.

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
regression coefficient estimates. The matrix has several attributes:

- type:

  indicates which small-sample adjustment was used

- cluster:

  contains the factor vector that defines independent clusters

- bread:

  contains the bread matrix

- v_scale:

  constant used in scaling the sandwich estimator

- est_mats:

  contains a list of estimating matrices used to calculate the sandwich
  estimator

- adjustments:

  contains a list of adjustment matrices used to calculate the sandwich
  estimator

- target:

  contains the working variance-covariance model used to calculate the
  adjustment matrices. This is needed for calculating small-sample
  corrections for Wald tests.

## Details

`vcovCR` returns a sandwich estimate of the variance-covariance matrix
of a set of regression coefficient estimates.

Several different small sample corrections are available, which run
parallel with the "HC" corrections for heteroskedasticity-consistent
variance estimators, as implemented in
[`vcovHC`](https://sandwich.R-Forge.R-project.org/reference/vcovHC.html).
The "CR2" adjustment is recommended (Pustejovsky & Tipton, 2017; Imbens
& Kolesar, 2016). See Pustejovsky and Tipton (2017) and Cameron and
Miller (2015) for further technical details. Available options include:

- "CR0":

  is the original form of the sandwich estimator (Liang & Zeger, 1986),
  which does not make any small-sample correction.

- "CR1":

  multiplies CR0 by `m / (m - 1)`, where `m` is the number of clusters.

- "CR1p":

  multiplies CR0 by `m / (m - p)`, where `m` is the number of clusters
  and `p` is the number of covariates.

- "CR1S":

  multiplies CR0 by `(m (N-1)) / [(m - 1)(N - p)]`, where `m` is the
  number of clusters, `N` is the total number of observations, and `p`
  is the number of covariates. Some Stata commands use this correction
  by default.

- "CR2":

  is the "bias-reduced linearization" adjustment proposed by Bell and
  McCaffrey (2002) and further developed in Pustejovsky and Tipton
  (2017). The adjustment is chosen so that the variance-covariance
  estimator is exactly unbiased under a user-specified working model.

- "CR3":

  approximates the leave-one-cluster-out jackknife variance estimator
  (Bell & McCaffrey, 2002).

## References

Bell, R. M., & McCaffrey, D. F. (2002). Bias reduction in standard
errors for linear regression with multi-stage samples. Survey
Methodology, 28(2), 169-181.

Cameron, A. C., & Miller, D. L. (2015). A Practitioner's Guide to
Cluster-Robust Inference. *Journal of Human Resources, 50*(2), 317-372.
[doi:10.3368/jhr.50.2.317](https://doi.org/10.3368/jhr.50.2.317)

Imbens, G. W., & Kolesar, M. (2016). Robust standard errors in small
samples: Some practical advice. *Review of Economics and Statistics,
98*(4), 701-712.
[doi:10.1162/rest_a_00552](https://doi.org/10.1162/rest_a_00552)

Liang, K.-Y., & Zeger, S. L. (1986). Longitudinal data analysis using
generalized linear models. *Biometrika, 73*(1), 13-22.
[doi:10.1093/biomet/73.1.13](https://doi.org/10.1093/biomet/73.1.13)

Pustejovsky, J. E. & Tipton, E. (2018). Small sample methods for
cluster-robust variance estimation and hypothesis testing in fixed
effects models. *Journal of Business and Economic Statistics, 36*(4),
672-683.
[doi:10.1080/07350015.2016.1247004](https://doi.org/10.1080/07350015.2016.1247004)

## See also

[`vcovCR.lm`](http://jepusto.github.io/clubSandwich/reference/vcovCR.lm.md),
[`vcovCR.plm`](http://jepusto.github.io/clubSandwich/reference/vcovCR.plm.md),
[`vcovCR.glm`](http://jepusto.github.io/clubSandwich/reference/vcovCR.glm.md),
[`vcovCR.gls`](http://jepusto.github.io/clubSandwich/reference/vcovCR.gls.md),
[`vcovCR.lme`](http://jepusto.github.io/clubSandwich/reference/vcovCR.lme.md),
[`vcovCR.lmerMod`](http://jepusto.github.io/clubSandwich/reference/vcovCR.lmerMod.md),
[`vcovCR.robu`](http://jepusto.github.io/clubSandwich/reference/vcovCR.robu.md),
[`vcovCR.rma.uni`](http://jepusto.github.io/clubSandwich/reference/vcovCR.rma.uni.md),
[`vcovCR.rma.mv`](http://jepusto.github.io/clubSandwich/reference/vcovCR.rma.mv.md),
[`vcovCR.lm_robust`](http://jepusto.github.io/clubSandwich/reference/vcovCR.lm_robust.md)

## Examples

``` r

# simulate design with cluster-dependence
m <- 8
cluster <- factor(rep(LETTERS[1:m], 3 + rpois(m, 5)))
n <- length(cluster)
X <- matrix(rnorm(3 * n), n, 3)
nu <- rnorm(m)[cluster]
e <- rnorm(n)
y <- X %*% c(.4, .3, -.3) + nu + e
dat <- data.frame(y, X, cluster, row = 1:n)

# fit linear model
lm_fit <- lm(y ~ X1 + X2 + X3, data = dat)
vcov(lm_fit)
#>              (Intercept)            X1            X2           X3
#> (Intercept)  0.024119296 -0.0033629563  0.0030753483  0.001678683
#> X1          -0.003362956  0.0175124039 -0.0001519176 -0.002512630
#> X2           0.003075348 -0.0001519176  0.0166577695  0.003156097
#> X3           0.001678683 -0.0025126298  0.0031560972  0.024070803

# cluster-robust variance estimator with CR2 small-sample correction
vcovCR(lm_fit, cluster = dat$cluster, type = "CR2")
#>              (Intercept)           X1           X2          X3
#> (Intercept)  0.118258531 -0.058199675 -0.019174051 0.002266691
#> X1          -0.058199675  0.043149660  0.013660753 0.002410578
#> X2          -0.019174051  0.013660753  0.011604149 0.001719445
#> X3           0.002266691  0.002410578  0.001719445 0.001876608

# compare small-sample adjustments
CR_types <- paste0("CR",c("0","1","1S","2","3"))
sapply(CR_types, function(type) 
       sqrt(diag(vcovCR(lm_fit, cluster = dat$cluster, type = type))))
#>                    CR0       CR1       CR1S        CR2        CR3
#> (Intercept) 0.31083475 0.3322963 0.34011613 0.34388738 0.38123090
#> X1          0.18390490 0.1966026 0.20122919 0.20772496 0.23542350
#> X2          0.09753920 0.1042738 0.10672762 0.10772256 0.11989761
#> X3          0.04068856 0.0434979 0.04452152 0.04331983 0.04774566
```
