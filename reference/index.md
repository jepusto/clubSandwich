# Package index

## Cluster-robust variance covariance estimation

S3 methods for calculating cluster-robust variance-covariance matrices
for various fitted model objects.

- [`vcovCR()`](http://jepusto.github.io/clubSandwich/reference/vcovCR.md)
  : Cluster-robust variance-covariance matrix

- [`vcovCR(`*`<geeglm>`*`)`](http://jepusto.github.io/clubSandwich/reference/vcovCR.geeglm.md)
  : Cluster-robust variance-covariance matrix for a geeglm object.

- [`vcovCR(`*`<glm>`*`)`](http://jepusto.github.io/clubSandwich/reference/vcovCR.glm.md)
  : Cluster-robust variance-covariance matrix for a glm object.

- [`vcovCR(`*`<gls>`*`)`](http://jepusto.github.io/clubSandwich/reference/vcovCR.gls.md)
  : Cluster-robust variance-covariance matrix for a gls object.

- [`vcovCR(`*`<ivreg>`*`)`](http://jepusto.github.io/clubSandwich/reference/vcovCR.ivreg.md)
  : Cluster-robust variance-covariance matrix for an ivreg object.

- [`vcovCR(`*`<lm>`*`)`](http://jepusto.github.io/clubSandwich/reference/vcovCR.lm.md)
  : Cluster-robust variance-covariance matrix for an lm object.

- [`vcovCR(`*`<lm_robust>`*`)`](http://jepusto.github.io/clubSandwich/reference/vcovCR.lm_robust.md)
  :

  Cluster-robust variance-covariance matrix for an
  [`estimatr::lm_robust`](https://declaredesign.org/r/estimatr/reference/lm_robust.html)
  object.

- [`vcovCR(`*`<lme>`*`)`](http://jepusto.github.io/clubSandwich/reference/vcovCR.lme.md)
  : Cluster-robust variance-covariance matrix for an lme object.

- [`vcovCR(`*`<lmerMod>`*`)`](http://jepusto.github.io/clubSandwich/reference/vcovCR.lmerMod.md)
  : Cluster-robust variance-covariance matrix for an lmerMod object.

- [`vcovCR(`*`<mlm>`*`)`](http://jepusto.github.io/clubSandwich/reference/vcovCR.mlm.md)
  : Cluster-robust variance-covariance matrix for an mlm object.

- [`vcovCR(`*`<mmrm>`*`)`](http://jepusto.github.io/clubSandwich/reference/vcovCR.mmrm.md)
  : Cluster-robust variance-covariance matrix for an mmrm object.

- [`vcovCR(`*`<plm>`*`)`](http://jepusto.github.io/clubSandwich/reference/vcovCR.plm.md)
  : Cluster-robust variance-covariance matrix for a plm object.

- [`vcovCR(`*`<rma.mv>`*`)`](http://jepusto.github.io/clubSandwich/reference/vcovCR.rma.mv.md)
  : Cluster-robust variance-covariance matrix for a rma.mv object.

- [`vcovCR(`*`<rma.uni>`*`)`](http://jepusto.github.io/clubSandwich/reference/vcovCR.rma.uni.md)
  : Cluster-robust variance-covariance matrix for a rma.uni object.

- [`vcovCR(`*`<robu>`*`)`](http://jepusto.github.io/clubSandwich/reference/vcovCR.robu.md)
  : Cluster-robust variance-covariance matrix for a robu object.

## Inference functions

Functions for inference based on cluster-robust variance-covariance
matrices.

- [`coef_test()`](http://jepusto.github.io/clubSandwich/reference/coef_test.md)
  : Test all or selected regression coefficients in a fitted model
- [`conf_int()`](http://jepusto.github.io/clubSandwich/reference/conf_int.md)
  : Calculate confidence intervals for all or selected regression
  coefficients in a fitted model
- [`linear_contrast()`](http://jepusto.github.io/clubSandwich/reference/linear_contrast.md)
  : Calculate confidence intervals and p-values for linear contrasts of
  regression coefficients in a fitted model
- [`Wald_test()`](http://jepusto.github.io/clubSandwich/reference/Wald_test.md)
  : Test parameter constraints in a fitted linear regression model
- [`constrain_zero()`](http://jepusto.github.io/clubSandwich/reference/constraint_matrices.md)
  [`constrain_equal()`](http://jepusto.github.io/clubSandwich/reference/constraint_matrices.md)
  [`constrain_pairwise()`](http://jepusto.github.io/clubSandwich/reference/constraint_matrices.md)
  : Create constraint matrices

## Helper functions for multivariate meta-analysis

Helper functions for use in estimating multivariate meta-analysis models
with
[`metafor::rma.mv()`](https://wviechtb.github.io/metafor/reference/rma.mv.html).

- [`impute_covariance_matrix()`](http://jepusto.github.io/clubSandwich/reference/impute_covariance_matrix.md)
  : Impute a block-diagonal covariance matrix
- [`pattern_covariance_matrix()`](http://jepusto.github.io/clubSandwich/reference/pattern_covariance_matrix.md)
  : Impute a patterned block-diagonal covariance matrix
- [`findCluster.rma.mv()`](http://jepusto.github.io/clubSandwich/reference/findCluster.rma.mv.md)
  : Detect cluster structure of an rma.mv object

## Data

Example datasets

- [`AchievementAwardsRCT`](http://jepusto.github.io/clubSandwich/reference/AchievementAwardsRCT.md)
  : Achievement Awards Demonstration program
- [`dropoutPrevention`](http://jepusto.github.io/clubSandwich/reference/dropoutPrevention.md)
  : Dropout prevention/intervention program effects
- [`MortalityRates`](http://jepusto.github.io/clubSandwich/reference/MortalityRates.md)
  : State-level annual mortality rates by cause among 18-20 year-olds
- [`SATcoaching`](http://jepusto.github.io/clubSandwich/reference/SATcoaching.md)
  : Randomized experiments on SAT coaching
