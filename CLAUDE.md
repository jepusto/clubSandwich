# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

clubSandwich is an R package providing cluster-robust variance estimators (sandwich estimators) with small-sample corrections for a wide range of regression models. The key methodological contribution is the bias-reduced linearization (BRL) estimator (CR2) from Bell & McCaffrey (2002) / Pustejovsky & Tipton (2017). Version 0.6.2.9999 (dev), authored by James E. Pustejovsky, licensed GPL-3.

## Build and test commands

```bash
# Full package check (runs R CMD check with all standard checks)
R CMD build . && R CMD check clubSandwich_*.tar.gz

# Run all tests
Rscript -e 'devtools::test()'

# Run a single test file
Rscript -e 'testthat::test_file("tests/testthat/test_lm.R")'

# Build documentation (roxygen2)
Rscript -e 'devtools::document()'

# Install from local source
Rscript -e 'devtools::install()'

# Code coverage
Rscript -e 'covr::package_coverage()'
```

## Architecture

### S3 dispatch chain

The core entry point is `vcovCR()`, a generic that dispatches to model-specific methods:

```
vcovCR(obj, cluster, type, target, inverse_var, form)
  └─ vcovCR.lm() / vcovCR.glm() / vcovCR.lme() / ... (13 model classes)
      └─ vcov_CR()  [core implementation in clubSandwich.R]
```

`vcov_CR()` orchestrates:
1. Process cluster variable → `droplevels(as.factor(cluster))`
2. Get design matrix → `model_matrix(obj)`
3. Get weights → `weightMatrix(obj, cluster)`
4. Get working variance target → `targetVariance(obj, cluster)`
5. Compute small-sample adjustment → `CR0()`/`CR1()`/`CR2()`/`CR3()`/`CR4()` (in `CR-adjustments.R`)
6. Sandwich construction → `bread %*% meat %*% bread / v_scale`

### Key source files

- **`clubSandwich.R`** — `vcov_CR()` core implementation + `vcovCR()` generic
- **`S3-methods.R`** — Base S3 generics: `targetVariance()`, `weightMatrix()`, `model_matrix()`, `residuals_CS()`, `coef_CS()`, `v_scale()`
- **`CR-adjustments.R`** — Small-sample correction functions (CR0–CR4), `IH_jj_list()` for projection matrices
- **`get_arrays.R`** — `get_GH()`, `get_P_array()`, `get_S_array()` for Satterthwaite/saddlepoint computations
- **`matrix-functions.R`** — Block matrix utilities: `matrix_list()`, `matrix_power()`, `chol_psd()`
- **Model files** (`lm.R`, `glm.R`, `gls.R`, `lme.R`, `lmer.R`, `plm.R`, `robu.R`, `rma-uni.R`, `rma-mv.R`, etc.) — Each implements the S3 methods above for its model class

### Adding support for a new model class

Create a new file `R/newmodel.R` implementing these S3 methods:

1. `vcovCR.NewModel()` — Entry point, extracts cluster if auto-detectable, sets `inverse_var` default, calls `vcov_CR()`
2. `model_matrix.NewModel()` — Return the design matrix X
3. `residuals_CS.NewModel()` — Return working residuals
4. `targetVariance.NewModel()` — Return list of working variance-covariance blocks per cluster
5. `weightMatrix.NewModel()` — Return list of weight matrices per cluster
6. `v_scale.NewModel()` — Return variance scaling constant (typically N or N-p)
7. `bread.NewModel()` — Only if not already provided by the `sandwich` package

### Inference functions

- **`coef_test()`** (`coef_test.R`) — t-tests with Satterthwaite or saddlepoint df corrections
- **`conf_int()`** (`conf_int.R`) — Confidence intervals with small-sample corrections
- **`Wald_test()`** (`Wald_test.R`) — Multi-contrast F-tests using Hotelling's T² approximation
- **`linear_contrast()`** — Linear contrast interface wrapping `coef_test()`
- **`constrain_zero()`**, **`constrain_equal()`**, **`constrain_pairwise()`** — Constraint matrix builders for `Wald_test()`

## Testing conventions

Tests use `testthat` in `tests/testthat/`. Every model-type test file follows this pattern:

1. **Setup**: `set.seed(20190513)`, `skip_if_not_installed()` for optional deps, simulate/load data, fit models
2. **Standard checks** (in `R/utilities.R`):
   - `check_bread(obj, cluster, y)` — Reconstructs (X'WX)^{-1} and verifies coefficients match
   - `check_CR(obj, vcov, cluster)` — Verifies E(V^CR_j) equals the target (core unbiasedness property)
   - `check_sort_order(obj, dat, cluster)` — Scrambles rows, refits, confirms identical results
   - `compare_ttests(a, b)` / `compare_Waldtests(a, b)` — Compare result sets
3. **Recurring test categories**: bread works, vcovCR options for CR2/CR4, target-unbiasedness, equivalence to `sandwich::vcovHC` for size-1 clusters, dropped observations, aliased predictors, weight scaling invariance, single-cluster error

## Dependencies

**Imports**: `stats`, `sandwich` (for `bread()`), `lifecycle`

**Suggests** (test/vignette only): `testthat`, `nlme`, `lme4`, `plm`, `metafor`, `metadat`, `robumeta`, `geepack`, `AER`, `ivreg`, `estimatr`, `carData`, `mlmRev`, `Matrix`, `zoo`, `fixest`, `Formula`, `knitr`, `rmarkdown`, `covr`
