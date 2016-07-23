[![Travis-CI Build Status](https://travis-ci.org/jepusto/clubSandwich.svg?branch=master)](https://travis-ci.org/jepusto/clubSandwich)
[![Coverage Status](https://img.shields.io/codecov/c/github/jepusto/clubSandwich/master.svg)](https://codecov.io/github/jepusto/clubSandwich?branch=master)
[![](http://www.r-pkg.org/badges/version/clubSandwich)](http://cran.rstudio.com/web/packages/clubSandwich/index.html)
[![](http://cranlogs.r-pkg.org/badges/grand-total/clubSandwich)](http://cran.rstudio.com/web/packages/clubSandwich/index.html)

# clubSandwich

`clubSandwich` provides several cluster-robust variance estimators 
(i.e., sandwich estimators) for ordinary and weighted least squares linear regression models. 
Several adjustments are incorporated to improve small-sample performance. 
The package includes functions for estimating the variance-covariance matrix and 
for testing single- and multiple-contrast hypotheses based on Wald test statistics. 
Tests of single regression coefficients use Satterthwaite or saddlepoint corrections.
Tests of multiple-contrast hypotheses use an approximation to Hotelling's T-squared distribution. 
Methods are provided for a variety of fitted models, including  `lm`, `plm` (from package `plm`), `gls` and `lme` (from `nlme`), `robu` (from `robumeta`), and `rma.uni` and `rma.mv` (from `metafor`). 

# Installing clubSandwich

The package is available on the Comprehensive R Archive Network. To install it, type 
```{r}
install.packages("clubSandwich")
```

To install the latest development version directly from Github, type:
```{r}
install.packages("devtools")
devtools::install_github("jepusto/clubSandwich")
```
