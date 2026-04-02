# Cluster-robust variance-covariance matrix for an mmrm object.

`vcovCR` returns a sandwich estimate of the variance-covariance matrix
of a set of regression coefficient estimates from an
[`mmrm`](https://openpharma.github.io/mmrm/latest-tag/reference/mmrm.html)
object.

## Usage

``` r
# S3 method for class 'mmrm'
vcovCR(obj, cluster, type, target, inverse_var, form = "sandwich", ...)
```

## Arguments

- obj:

  Fitted model for which to calculate the variance-covariance matrix

- cluster:

  Optional expression or vector indicating which observations belong to
  the same cluster. If not specified, will be set to the subject
  variable from the `mmrm` model formula.

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
  variance-covariance structure of the `mmrm` object, with weights
  applied if present.

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
if (requireNamespace("mmrm", quietly = TRUE)) withAutoprint({

  library(mmrm)
  data(fev_data, package = "mmrm")

  # Fit an mmrm model with unstructured covariance
  mmrm_fit <- mmrm(
    FEV1 ~ RACE + SEX + ARMCD * AVISIT + us(AVISIT | USUBJID),
    data = fev_data
  )

  # CR2 cluster-robust variance estimator (cluster auto-detected)
  vcovCR(mmrm_fit, type = "CR2")

  # Coefficient tests with Satterthwaite degrees of freedom
  coef_test(mmrm_fit, vcov = "CR2", test = "Satterthwaite")

  # Fit a weighted mmrm model
  fev_data$wt <- 1 + rpois(nrow(fev_data), lambda = 3)
  mmrm_wt <- mmrm(
    FEV1 ~ RACE + SEX + ARMCD + us(AVISIT | USUBJID),
    data = fev_data,
    weights = fev_data$wt
  )

  # CR2 works with weighted models
  vcovCR(mmrm_wt, type = "CR2")

})
#> > library(mmrm)
#> > data(fev_data, package = "mmrm")
#> > mmrm_fit <- mmrm(FEV1 ~ RACE + SEX + ARMCD * AVISIT + us(AVISIT | USUBJID), 
#> +     data = fev_data)
#> > vcovCR(mmrm_fit, type = "CR2")
#>                               (Intercept) RACEBlack or African American
#> (Intercept)                     0.7486699                  -0.275377814
#> RACEBlack or African American  -0.2753778                   0.399818190
#> RACEWhite                      -0.2193198                   0.222373089
#> SEXFemale                      -0.1744407                   0.029046368
#> ARMCDTRT                       -0.4509784                   0.008250065
#> AVISITVIS2                     -0.3227883                  -0.038293391
#> AVISITVIS3                     -0.4547670                   0.021404213
#> AVISITVIS4                     -0.4201979                   0.033844036
#> ARMCDTRT:AVISITVIS2             0.2848088                   0.111528886
#> ARMCDTRT:AVISITVIS3             0.4401008                   0.024412918
#> ARMCDTRT:AVISITVIS4             0.3785846                   0.057210773
#>                                  RACEWhite    SEXFemale     ARMCDTRT
#> (Intercept)                   -0.219319791 -0.174440682 -0.450978374
#> RACEBlack or African American  0.222373089  0.029046368  0.008250065
#> RACEWhite                      0.489112634  0.080650472 -0.141644549
#> SEXFemale                      0.080650472  0.289953729 -0.027945701
#> ARMCDTRT                      -0.141644549 -0.027945701  1.184979675
#> AVISITVIS2                    -0.056550656  0.019401793  0.343856320
#> AVISITVIS3                    -0.005112526 -0.009524999  0.454656341
#> AVISITVIS4                     0.017644969  0.037878067  0.384949958
#> ARMCDTRT:AVISITVIS2            0.126475338 -0.035829368 -0.850549586
#> ARMCDTRT:AVISITVIS3            0.027960034 -0.016350149 -1.046822315
#> ARMCDTRT:AVISITVIS4            0.013297790 -0.046038291 -0.868455239
#>                                AVISITVIS2   AVISITVIS3  AVISITVIS4
#> (Intercept)                   -0.32278835 -0.454767019 -0.42019794
#> RACEBlack or African American -0.03829339  0.021404213  0.03384404
#> RACEWhite                     -0.05655066 -0.005112526  0.01764497
#> SEXFemale                      0.01940179 -0.009524999  0.03787807
#> ARMCDTRT                       0.34385632  0.454656341  0.38494996
#> AVISITVIS2                     0.52294746  0.321243882  0.21005177
#> AVISITVIS3                     0.32124388  0.628571724  0.35259632
#> AVISITVIS4                     0.21005177  0.352596325  1.52761435
#> ARMCDTRT:AVISITVIS2           -0.52558755 -0.319810333 -0.21130872
#> ARMCDTRT:AVISITVIS3           -0.32449443 -0.626073105 -0.35340385
#> ARMCDTRT:AVISITVIS4           -0.20868869 -0.350202589 -1.52838250
#>                               ARMCDTRT:AVISITVIS2 ARMCDTRT:AVISITVIS3
#> (Intercept)                            0.28480875          0.44010080
#> RACEBlack or African American          0.11152889          0.02441292
#> RACEWhite                              0.12647534          0.02796003
#> SEXFemale                             -0.03582937         -0.01635015
#> ARMCDTRT                              -0.85054959         -1.04682232
#> AVISITVIS2                            -0.52558755         -0.32449443
#> AVISITVIS3                            -0.31981033         -0.62607311
#> AVISITVIS4                            -0.21130872         -0.35340385
#> ARMCDTRT:AVISITVIS2                    1.28044209          0.81192593
#> ARMCDTRT:AVISITVIS3                    0.81192593          1.40680254
#> ARMCDTRT:AVISITVIS4                    0.72261379          0.79805570
#>                               ARMCDTRT:AVISITVIS4
#> (Intercept)                            0.37858456
#> RACEBlack or African American          0.05721077
#> RACEWhite                              0.01329779
#> SEXFemale                             -0.04603829
#> ARMCDTRT                              -0.86845524
#> AVISITVIS2                            -0.20868869
#> AVISITVIS3                            -0.35020259
#> AVISITVIS4                            -1.52838250
#> ARMCDTRT:AVISITVIS2                    0.72261379
#> ARMCDTRT:AVISITVIS3                    0.79805570
#> ARMCDTRT:AVISITVIS4                    3.41530964
#> > coef_test(mmrm_fit, vcov = "CR2", test = "Satterthwaite")
#> Alternative hypothesis: two-sided 
#>                          Coef. Estimate    SE Null value  t-stat d.f. (Satt)
#>                    (Intercept)  30.7775 0.865          0 35.5703        90.2
#>  RACEBlack or African American   1.5305 0.632          0  2.4205       117.9
#>                      RACEWhite   5.6436 0.699          0  8.0695        94.4
#>                      SEXFemale   0.3261 0.538          0  0.6055       153.8
#>                       ARMCDTRT   3.7742 1.089          0  3.4672       142.4
#>                     AVISITVIS2   4.8396 0.723          0  6.6924        68.9
#>                     AVISITVIS3  10.3421 0.793          0 13.0446        77.3
#>                     AVISITVIS4  15.0539 1.236          0 12.1799        73.9
#>            ARMCDTRT:AVISITVIS2  -0.0419 1.132          0 -0.0371       136.1
#>            ARMCDTRT:AVISITVIS3  -0.6937 1.186          0 -0.5849       149.3
#>            ARMCDTRT:AVISITVIS4   0.6242 1.848          0  0.3378       144.8
#>  p-val (Satt) Sig.
#>        <0.001  ***
#>         0.017    *
#>        <0.001  ***
#>         0.546     
#>        <0.001  ***
#>        <0.001  ***
#>        <0.001  ***
#>        <0.001  ***
#>         0.970     
#>         0.560     
#>         0.736     
#> > fev_data$wt <- 1 + rpois(nrow(fev_data), lambda = 3)
#> > mmrm_wt <- mmrm(FEV1 ~ RACE + SEX + ARMCD + us(AVISIT | USUBJID), data = fev_data, 
#> +     weights = fev_data$wt)
#> > vcovCR(mmrm_wt, type = "CR2")
#>                               (Intercept) RACEBlack or African American
#> (Intercept)                     0.4280773                   -0.29394827
#> RACEBlack or African American  -0.2939483                    0.50963863
#> RACEWhite                      -0.2703860                    0.24747987
#> SEXFemale                      -0.2040653                    0.02679204
#> ARMCDTRT                       -0.1428009                    0.06503542
#>                                 RACEWhite   SEXFemale    ARMCDTRT
#> (Intercept)                   -0.27038600 -0.20406529 -0.14280092
#> RACEBlack or African American  0.24747987  0.02679204  0.06503542
#> RACEWhite                      0.57879567  0.08783448 -0.07998619
#> SEXFemale                      0.08783448  0.35577426 -0.02560111
#> ARMCDTRT                      -0.07998619 -0.02560111  0.39340382
```
