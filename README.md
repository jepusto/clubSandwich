# clubSandwich

`clubSandwich` provides several cluster-robust variance estimators 
(i.e., sandwich estimators) for ordinary and weighted least squares linear regression models. 
Several adjustments are incorporated to improve small-sample performance. 
The package includes functions for estimating the variance-covariance matrix and 
for testing single- and multiple-contrast hypotheses based on Wald test statistics. 
Tests of single regression coefficients use Satterthwaite or saddlepoint corrections.
Tests of multiple-contrast hypotheses use adjustments based on a Hotelling's T-squared
approximation. Methods are provided for a variety of fitted models, including 
`lm`, `plm` (from package `plm`), `gls` and `lme` (from `nlme`), `robu` (from `robumeta`), and 
`rma.uni` and `rma.mv` (from `metafor`). 

# Installing clubSandwich

Currently, the package is only available here on Github. To install it, type the following commands in the R console:
```{r}
install.packages("devtools")
library(devtools)
install_github("jepusto/clubSandwich")
```