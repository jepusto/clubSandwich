<!-- badges: start -->
[![R-CMD-check](https://github.com/jepusto/clubSandwich/workflows/R-CMD-check/badge.svg)](https://github.com/jepusto/clubSandwich/actions)
[![Codecov Status](https://codecov.io/gh/jepusto/clubSandwich/branch/main/graph/badge.svg)](https://codecov.io/gh/jepusto/clubSandwich?branch=main)
[![CRAN status](https://www.r-pkg.org/badges/version/clubSandwich)](https://CRAN.R-project.org/package=clubSandwich)
[![](https://cranlogs.r-pkg.org/badges/grand-total/clubSandwich)](https://CRAN.R-project.org/package=clubSandwich)
[![](https://cranlogs.r-pkg.org/badges/last-month/clubSandwich)](https://CRAN.R-project.org/package=clubSandwich)
<!-- badges: end -->

# clubSandwich

`clubSandwich` provides several cluster-robust variance estimators 
(i.e., sandwich estimators) for ordinary and weighted least squares linear regression models, two-stage least squares regression models, and generalized linear models. 
Several adjustments are incorporated to improve small-sample performance. 
The package includes functions for estimating the variance-covariance matrix and 
for testing single- and multiple-contrast hypotheses based on Wald test statistics. 
Tests of single regression coefficients use Satterthwaite or saddlepoint corrections.
Tests of multiple-contrast hypotheses use an approximation to Hotelling's T-squared distribution. 
Methods are provided for a variety of fitted models, including:

- `lm()`
- `mlm()`
- `glm()` 
- `geeglm()` (from package `geepack`)
- `ivreg()` (from package `ivreg`, when estimated using `method = "OLS"`)
- `ivreg()` (from package `AER`)
- `plm()` (from package `plm`), 
- `gls()` and `lme()` (from package `nlme`)
- `lmer()` (from package `lme4`)
- `robu()` (from package `robumeta`)
- `rma.uni()` and `rma.mv()` (from package `metafor`) 

# Installing clubSandwich

The package is available on the Comprehensive R Archive Network. To install it, type 
```{r}
install.packages("clubSandwich")
```

To install the latest development version directly from Github, type:
```{r}
install.packages("remotes")
remotes::install_github("jepusto/clubSandwich")
```

Once installed, have a look at the available vignettes by typing:
```{r}
browseVignettes(package="clubSandwich")
```
