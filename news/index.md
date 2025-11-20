# Changelog

## clubSandwich 0.6.1.9999

- Fixed bug in internal function get_cor_grouping(), which identifies
  the grouping variable in nlme::lme() models with correlation
  structures.
- Fixed bug in bread.rma.mv() occurring in models with test = ‘knha’ (an
  undocumented, experimental option).

## clubSandwich 0.6.1

CRAN release: 2025-07-30

- Added support for
  [`estimatr::lm_robust()`](https://declaredesign.org/r/estimatr/reference/lm_robust.html)
  and
  [`estimatr::lm_lin()`](https://declaredesign.org/r/estimatr/reference/lm_lin.html)
  models.

## clubSandwich 0.6.0

CRAN release: 2025-04-01

- Added option to Wald_test() and linear_contrast() to correct
  hypothesis tests for multiple comparisons.
- Added option to coef_test() to specify one- or two-sided alternative
  hypotheses.
- Added options to Wald_test() and coef_test() to specify non-zero
  values for null hypotheses.
- Fixed a bug in methods for
  [`geepack::geeglm()`](https://rdrr.io/pkg/geepack/man/geeglm.html)
  models that occurred for models with nonlinear link functions.
- Corrected typo in deprecation messages for
  [`impute_covariance_matrix()`](http://jepusto.github.io/clubSandwich/reference/impute_covariance_matrix.md)
  and
  [`pattern_covariance_matrix()`](http://jepusto.github.io/clubSandwich/reference/pattern_covariance_matrix.md).

## clubSandwich 0.5.11

CRAN release: 2024-06-20

- Corrected a unit test related to the plm package, for compatibility
  with upcoming release of plm.
- Deprecated
  [`impute_covariance_matrix()`](http://jepusto.github.io/clubSandwich/reference/impute_covariance_matrix.md)
  and
  [`pattern_covariance_matrix()`](http://jepusto.github.io/clubSandwich/reference/pattern_covariance_matrix.md),
  because they have been superseded by
  [`metafor::vcalc()`](https://wviechtb.github.io/metafor/reference/vcalc.html).

## clubSandwich 0.5.10

CRAN release: 2023-07-20

- Fixed another bug in
  [`linear_contrast()`](http://jepusto.github.io/clubSandwich/reference/linear_contrast.md)
  to handle specified contrasts that are scalars when
  variance-covariance matrix is computed with a working model that is
  not inverse-variance.
- Fixed formatting of package version numbers in unit tests to conform
  to changes in
  [`packageVersion()`](https://rdrr.io/r/utils/packageDescription.html)
  in R-devel.

## clubSandwich 0.5.9

CRAN release: 2023-07-12

- Added support for
  [`geepack::geeglm()`](https://rdrr.io/pkg/geepack/man/geeglm.html)
  models.
- Added support for `rma.ls` models (location-scale models estimated
  using `metafor::rma.uni(scale = )`).
- Improved error handling of
  [`Wald_test()`](http://jepusto.github.io/clubSandwich/reference/Wald_test.md)
  when vcov of contrasts is not positive definite.
- Fixed a bug in
  [`linear_contrast()`](http://jepusto.github.io/clubSandwich/reference/linear_contrast.md)
  to handle specified contrasts that are scalars.
- Improved internal `get_data` function for `gls` and `lme` objects to
  allow use of expressions in addition to object names.

## clubSandwich 0.5.8

CRAN release: 2022-08-15

- Added support for
  [`ivreg::ivreg`](https://zeileis.github.io/ivreg/reference/ivreg.html)
  objects when estimated by ordinary least squares (support for objects
  estimated by 2SM and 2SMM is not yet implemented).
- Updated unit tests for
  [`plm::plm()`](https://rdrr.io/pkg/plm/man/plm.html) when
  `method = "FD"` to account for bug fixes in version 2.6-2 of plm.

## clubSandwich 0.5.7

CRAN release: 2022-06-15

- Fixed bug in methods for multi-variate multi-level models estimated
  with lme().
- Updated vignettes, examples, and unit tests so that the package can be
  compiled without any packages from SUGGESTS.

## clubSandwich 0.5.6

CRAN release: 2022-04-23

- Corrected bug in methods for `plm` objects estimated by random
  effects, which occurred when a user-specified clustering variable was
  at a higher level than the random effects.
- Added support for `plm` objects with nested random effects
  (`effects = "nested"`).
- Added additional syntactic options for specifying clustering variable
  with `plm` objects. See
  [`?plm`](https://rdrr.io/pkg/plm/man/plm.html).
- Corrected bug in how
  [`Wald_test()`](http://jepusto.github.io/clubSandwich/reference/Wald_test.md)
  labeled results when `test = "Naive-Fp"`.

## clubSandwich 0.5.5

CRAN release: 2022-01-18

- New function
  [`linear_contrast()`](http://jepusto.github.io/clubSandwich/reference/linear_contrast.md)
  calculates robust confidence intervals and p-values for linear
  contrasts of regression coefficients from a fitted model. Works with
  [`constrain_pairwise()`](http://jepusto.github.io/clubSandwich/reference/constraint_matrices.md)
  and other `constrain_*()` helper functions.

- Corrected precision of unit test leading to error on M1mac.

## clubSandwich 0.5.4

CRAN release: 2022-01-09

### New features

- [`Wald_test()`](http://jepusto.github.io/clubSandwich/reference/Wald_test.md)
  gains an option for `test = "Naive-Fp"`, which uses denominator
  degrees of freedom equal to the number of clusters minus the number of
  coefficients in the fitted model.
- [`coef_test()`](http://jepusto.github.io/clubSandwich/reference/coef_test.md)
  and
  [`conf_int()`](http://jepusto.github.io/clubSandwich/reference/conf_int.md)
  gain an option for `test = "naive-tp"`, which uses denominator degrees
  of freedom equal to the number of clusters minus the number of
  coefficients in the fitted model.

### Minor improvements and bug fixes

- Corrected a bug in the Satterthwaite degrees of freedom calculations
  for models that include only an intercept.

- Output of
  [`coef_test()`](http://jepusto.github.io/clubSandwich/reference/coef_test.md)
  and
  [`conf_int()`](http://jepusto.github.io/clubSandwich/reference/conf_int.md)
  now include a variable containing the coefficient names, so that the
  results are “tidy.”

- [`conf_int()`](http://jepusto.github.io/clubSandwich/reference/conf_int.md)
  now includes an option to report a p-value for each coefficient.

- [`coef_test()`](http://jepusto.github.io/clubSandwich/reference/coef_test.md)
  now reports degrees of freedom for `test = 'z'` and
  `test = 'naive-t'`.

- [`vcovCR()`](http://jepusto.github.io/clubSandwich/reference/vcovCR.md)
  now provides a more informative error message when the clustering
  variable is a constant.

- [`vcovCR()`](http://jepusto.github.io/clubSandwich/reference/vcovCR.md)
  now handles models estimated using analytic weights, where some
  weights are equal to zero. Results are consistent with omitting
  observations with weights of zero.

- Added more informative error messages for
  [`Wald_test()`](http://jepusto.github.io/clubSandwich/reference/Wald_test.md)
  and
  [`conf_int()`](http://jepusto.github.io/clubSandwich/reference/conf_int.md),
  triggered if the test argument does not match any of the available
  tests.

- Corrected a bug in
  [`findCluster.rma.mv()`](http://jepusto.github.io/clubSandwich/reference/findCluster.rma.mv.md),
  which threw an error if a random effects factor in the rma.mv model
  had unobserved levels.

- Corrected a bug in
  [`Wald_test()`](http://jepusto.github.io/clubSandwich/reference/Wald_test.md),
  which threw an error for tests of intercept-only models.

- Fixed a minor bug in print method for
  [`Wald_test()`](http://jepusto.github.io/clubSandwich/reference/Wald_test.md)
  results, which threw an error when the p-value was `NA`.

## clubSandwich 0.5.3

CRAN release: 2021-01-23

- Removed dependency on mathjaxr

## clubSandwich 0.5.2

CRAN release: 2020-11-14

- Added mathjaxr to Imports

## clubSandwich 0.5.1

CRAN release: 2020-10-12

### New features

- New functionality for
  [`impute_covariance_matrix()`](http://jepusto.github.io/clubSandwich/reference/impute_covariance_matrix.md):
  - Compute covariance matrices with AR1 correlation structure or with a
    combination of constant correlation and AR1 correlation structure.
  - Compute covariance matrices that are blocked by subgroup.
  - Average the variance estimates by cluster before computing
    covariance matrices.
- New function
  [`pattern_covariance_matrix()`](http://jepusto.github.io/clubSandwich/reference/pattern_covariance_matrix.md),
  which creates a covariance matrix based on a specified pattern of
  correlations between different categories of effects.

### Minor improvements and bug fixes

- Corrected bug in methods for `rma.mv` objects, which previously led to
  incorrect identification of clustering variables in models with
  multiple levels of random effects, where at least one set of random
  effects has inner \| outer structure.

## clubSandwich 0.5.0

CRAN release: 2020-09-01

### New features: a major update to `Wald_test()`

- [`Wald_test()`](http://jepusto.github.io/clubSandwich/reference/Wald_test.md)
  now uses new helper functions
  [`constrain_zero()`](http://jepusto.github.io/clubSandwich/reference/constraint_matrices.md),
  [`constrain_equal()`](http://jepusto.github.io/clubSandwich/reference/constraint_matrices.md),
  and
  [`constrain_pairwise()`](http://jepusto.github.io/clubSandwich/reference/constraint_matrices.md)
  to specify constraint matrices.

- [`Wald_test()`](http://jepusto.github.io/clubSandwich/reference/Wald_test.md)
  gains an argument `tidy`. When `TRUE`, results for a list of tests
  will be tidied into a single data.frame.

- Output of
  [`Wald_test()`](http://jepusto.github.io/clubSandwich/reference/Wald_test.md)
  now includes both numerator and denominator degrees of freedom.

### Minor improvements and bug fixes

- Corrected bug in methods for `plm` objects, which occurred when
  “within” models included cluster-level interactions. Previously main
  effects of cluster-level variables were not getting dropped from
  `model_matrix.plm()`.

- Corrected bugs in methods for `robu` objects

  - Corrected a bug that previously led to errors for models with only
    one column in the model matrix (i.e., intercept-only models).
  - Corrected a bug in an internal function that previously led to
    errors in
    [`constrain_equal()`](http://jepusto.github.io/clubSandwich/reference/constraint_matrices.md)
    and
    [`constrain_zero()`](http://jepusto.github.io/clubSandwich/reference/constraint_matrices.md)
    when called on robu objects.

## clubSandwich 0.4.2

CRAN release: 2020-04-17

- Updated and streamlined unit tests for R 4.0.0.

## clubSandwich 0.4.1

CRAN release: 2020-01-28

- Updated unit tests to satisfy obscure CRAN checks.

## clubSandwich 0.4.0

CRAN release: 2019-12-18

- Added methods for `lmerMod` objects fitted by
  [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html).

- Updated internals to use
  [`inherits()`](https://rdrr.io/r/base/class.html) instead of checking
  [`class()`](https://rdrr.io/r/base/class.html) directly.

## clubSandwich 0.3.5

CRAN release: 2019-05-14

- Added t statistics to output of
  [`coef_test()`](http://jepusto.github.io/clubSandwich/reference/coef_test.md).

- Fixed a bug in `get_index_order()`, an internal function used with plm
  objects. Previously, the function assumed that both individual and
  time indices were specified in the `plm` call. The new function works
  even when zero or one indices are specified.

## clubSandwich 0.3.3

CRAN release: 2019-01-24

- [`impute_covariance_matrix()`](http://jepusto.github.io/clubSandwich/reference/impute_covariance_matrix.md)
  now drops unobserved factor levels.

- updated method for handling residuals from `rma.uni` and `rma.mv`
  objects, for consistency with metafor 2.1-0.

## clubSandwich 0.3.2

CRAN release: 2018-05-21

- Added
  [`conf_int()`](http://jepusto.github.io/clubSandwich/reference/conf_int.md)
  to provide easy cluster-robust confidence intervals.

- Added examples to documentation for
  [`conf_int()`](http://jepusto.github.io/clubSandwich/reference/conf_int.md)
  and
  [`coef_test()`](http://jepusto.github.io/clubSandwich/reference/coef_test.md).

## clubSandwich 0.3.1

CRAN release: 2018-04-04

- Added `coefs` option to
  [`coef_test()`](http://jepusto.github.io/clubSandwich/reference/coef_test.md)
  to allow testing of subsets of coefficients.

- Updated tests to use `carData` instead of car package.

## clubSandwich 0.3.0

CRAN release: 2017-11-13

- Added methods for `ivreg` objects.

- Added methods for `mlm` objects.

- Updated `residuals_CS.plm` to account for changes in plm 1.6-6.

## clubSandwich 0.2.3

CRAN release: 2017-08-10

### New features

- Added methods for `glm` objects.

- Provide facility to cluster at higher level than highest random
  effects for `lme` and `gls` objects.

- Added
  [`impute_covariance_matrix()`](http://jepusto.github.io/clubSandwich/reference/impute_covariance_matrix.md)
  utility function for multivariate meta-analysis.

### Minor improvements and bug fixes

- Updated methods for plm objects to account for changes in plm 1.6-6.

- Added documentation of `type` options in
  [`vcovCR()`](http://jepusto.github.io/clubSandwich/reference/vcovCR.md).

- Added examples for all
  [`vcovCR()`](http://jepusto.github.io/clubSandwich/reference/vcovCR.md)
  methods.

## clubSandwich 0.2.2

CRAN release: 2016-12-01

### New features

- Added
  [`bread()`](https://sandwich.R-Forge.R-project.org/reference/bread.html)
  methods for all supported model classes.

- [`vcovCR()`](http://jepusto.github.io/clubSandwich/reference/vcovCR.md)
  is now calculated using
  [`bread()`](https://sandwich.R-Forge.R-project.org/reference/bread.html),
  and carries attributes for `bread`, `est_mat`, and adjustment
  matrices.

- [`vcovCR()`](http://jepusto.github.io/clubSandwich/reference/vcovCR.md)
  gains a `form` argument to obtain just the meat of the sandwich, or to
  use a user-specified bread matrix.

### Minor improvements and bug fixes

- Refactored internal functions for degrees of freedom calculation to
  improve speed and memory usage.

- Bug fixes:

  - updated `nobs.plm()` method to handle first-differenced models

## clubSandwich 0.2.1

CRAN release: 2016-07-23

- First version released on CRAN.
