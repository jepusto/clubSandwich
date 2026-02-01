# clubSandwich 0.6.2

* Fixed bug in internal function `get_cor_grouping()`, which identifies the grouping variable in `nlme::lme()` models with correlation structures.
* Fixed bug in `bread.rma.mv()` occurring in models with `test = 'knha'` (an undocumented, experimental option).
* Fixed bug in `linear_contrast()` occurring for certain linear contrasts in models with weights that are not inverse variance. 

# clubSandwich 0.6.1

* Added support for `estimatr::lm_robust()` and `estimatr::lm_lin()` models.

# clubSandwich 0.6.0

* Added option to Wald_test() and linear_contrast() to correct hypothesis tests for multiple comparisons.
* Added option to coef_test() to specify one- or two-sided alternative hypotheses.
* Added options to Wald_test() and coef_test() to specify non-zero values for null hypotheses.
* Fixed a bug in methods for `geepack::geeglm()` models that occurred for models with nonlinear link functions.
* Corrected typo in deprecation messages for `impute_covariance_matrix()` and `pattern_covariance_matrix()`.

# clubSandwich 0.5.11

* Corrected a unit test related to the plm package, for compatibility with upcoming release of plm.
* Deprecated `impute_covariance_matrix()` and `pattern_covariance_matrix()`, because they have been superseded by `metafor::vcalc()`.

# clubSandwich 0.5.10

* Fixed another bug in `linear_contrast()` to handle specified contrasts that are scalars when variance-covariance matrix is computed with a working model that is not inverse-variance.
* Fixed formatting of package version numbers in unit tests to conform to changes in `packageVersion()` in R-devel.

# clubSandwich 0.5.9

* Added support for `geepack::geeglm()` models. 
* Added support for `rma.ls` models (location-scale models estimated using `metafor::rma.uni(scale = )`).
* Improved error handling of `Wald_test()` when vcov of contrasts is not positive definite.
* Fixed a bug in `linear_contrast()` to handle specified contrasts that are scalars.
* Improved internal `get_data` function for `gls` and `lme` objects to allow use of expressions in addition to object names.

# clubSandwich 0.5.8

* Added support for `ivreg::ivreg` objects when estimated by ordinary least squares (support for objects estimated by 2SM and 2SMM is not yet implemented).
* Updated unit tests for `plm::plm()` when `method = "FD"` to account for bug fixes in version 2.6-2 of plm. 

# clubSandwich 0.5.7

* Fixed bug in methods for multi-variate multi-level models estimated with lme().
* Updated vignettes, examples, and unit tests so that the package can be compiled without any packages from SUGGESTS.

# clubSandwich 0.5.6

* Corrected bug in methods for `plm` objects estimated by random effects, which occurred when a user-specified clustering variable was at a higher level than the random effects.
* Added support for `plm` objects with nested random effects (`effects = "nested"`).
* Added additional syntactic options for specifying clustering variable with `plm` objects. See `?plm`.
* Corrected bug in how `Wald_test()` labeled results when `test = "Naive-Fp"`.

# clubSandwich 0.5.5

* New function `linear_contrast()` calculates robust confidence intervals and p-values for linear contrasts of regression coefficients from a fitted model. Works with `constrain_pairwise()` and other `constrain_*()` helper functions.

* Corrected precision of unit test leading to error on M1mac.

# clubSandwich 0.5.4

## New features

* `Wald_test()` gains an option for `test = "Naive-Fp"`, which uses denominator degrees of freedom equal to the number of clusters minus the number of coefficients in the fitted model. 
* `coef_test()` and `conf_int()` gain an option for `test = "naive-tp"`, which uses denominator degrees of freedom equal to the number of clusters minus the number of coefficients in the fitted model.

## Minor improvements and bug fixes

* Corrected a bug in the Satterthwaite degrees of freedom calculations for models that include only an intercept.

* Output of `coef_test()` and `conf_int()` now include a variable containing the coefficient names, so that the results are "tidy."

* `conf_int()` now includes an option to report a p-value for each coefficient.

* `coef_test()` now reports degrees of freedom for `test = 'z'` and `test = 'naive-t'`.

* `vcovCR()` now provides a more informative error message when the clustering variable is a constant.

* `vcovCR()` now handles models estimated using analytic weights, where some weights are equal to zero. Results are consistent with omitting observations with weights of zero.

* Added more informative error messages for `Wald_test()` and `conf_int()`, triggered if the test argument does not match any of the available tests.

* Corrected a bug in `findCluster.rma.mv()`, which threw an error if a random effects factor in the rma.mv model had unobserved levels.

* Corrected a bug in `Wald_test()`, which threw an error for tests of intercept-only models.

* Fixed a minor bug in print method for `Wald_test()` results, which threw an error when the p-value was `NA`.

# clubSandwich 0.5.3

* Removed dependency on mathjaxr

# clubSandwich 0.5.2

* Added mathjaxr to Imports

# clubSandwich 0.5.1

## New features

* New functionality for `impute_covariance_matrix()`: 
    * Compute covariance matrices with AR1 correlation structure or with a combination of constant correlation and AR1 correlation structure. 
    * Compute covariance matrices that are blocked by subgroup.
    * Average the variance estimates by cluster before computing covariance matrices.

* New function `pattern_covariance_matrix()`, which creates a covariance matrix based on a specified pattern of correlations between different categories of effects. 

## Minor improvements and bug fixes

* Corrected bug in methods for `rma.mv` objects, which previously led to incorrect identification of clustering variables in models with multiple levels of random effects, where at least one set of random effects has inner | outer structure.


# clubSandwich 0.5.0

## New features: a major update to `Wald_test()`

* `Wald_test()` now uses new helper functions `constrain_zero()`, `constrain_equal()`, and `constrain_pairwise()` to specify constraint matrices.

* `Wald_test()` gains an argument `tidy`. When `TRUE`, results for a list of tests will be tidied into a single data.frame. 

* Output of `Wald_test()` now includes both numerator and denominator degrees of freedom.

## Minor improvements and bug fixes

* Corrected bug in methods for `plm` objects, which occurred when "within" models included cluster-level interactions. Previously main effects of cluster-level variables were not getting dropped from `model_matrix.plm()`.

* Corrected bugs in methods for `robu` objects
    * Corrected a bug that previously led to errors for models with only one column in the model matrix (i.e., intercept-only models). 
    * Corrected a bug in an internal function that previously led to errors in `constrain_equal()` and `constrain_zero()` when called on robu objects.


# clubSandwich 0.4.2

* Updated and streamlined unit tests for R 4.0.0.

# clubSandwich 0.4.1

* Updated unit tests to satisfy obscure CRAN checks.

# clubSandwich 0.4.0

* Added methods for `lmerMod` objects fitted by `lme4::lmer()`.

* Updated internals to use `inherits()` instead of checking `class()` directly.

# clubSandwich 0.3.5

* Added t statistics to output of `coef_test()`. 

* Fixed a bug in `get_index_order()`, an internal function used with plm objects. Previously, the function assumed that both individual and time indices were specified in the `plm` call. The new function works even when zero or one indices are specified.

# clubSandwich 0.3.3

* `impute_covariance_matrix()` now drops unobserved factor levels.

* updated method for handling residuals from `rma.uni` and `rma.mv` objects, for consistency with metafor 2.1-0.

# clubSandwich 0.3.2

* Added `conf_int()` to provide easy cluster-robust confidence intervals.

* Added examples to documentation for `conf_int()` and `coef_test()`.

# clubSandwich 0.3.1

* Added `coefs` option to `coef_test()` to allow testing of subsets of coefficients.

* Updated tests to use `carData` instead of car package.

# clubSandwich 0.3.0

* Added methods for `ivreg` objects.

* Added methods for `mlm` objects.

* Updated `residuals_CS.plm` to account for changes in plm 1.6-6.

# clubSandwich 0.2.3

## New features

* Added methods for `glm` objects.

* Provide facility to cluster at higher level than highest random effects for `lme` and `gls` objects.

* Added `impute_covariance_matrix()` utility function for multivariate meta-analysis.

## Minor improvements and bug fixes

* Updated methods for plm objects to account for changes in plm 1.6-6.

* Added documentation of `type` options in `vcovCR()`.

* Added examples for all `vcovCR()` methods.

# clubSandwich 0.2.2

## New features

* Added `bread()` methods for all supported model classes.

* `vcovCR()` is now calculated using `bread()`, and carries attributes for `bread`, `est_mat`, and adjustment matrices.

* `vcovCR()` gains a `form` argument to obtain just the meat of the sandwich, or to use a user-specified bread matrix. 


## Minor improvements and bug fixes

* Refactored internal functions for degrees of freedom calculation to improve speed and memory usage.

* Bug fixes:
  - updated `nobs.plm()` method to handle first-differenced models
  

# clubSandwich 0.2.1

* First version released on CRAN.
