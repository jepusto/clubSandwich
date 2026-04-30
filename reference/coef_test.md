# Test all or selected regression coefficients in a fitted model

`coef_test` reports one- or two-sided t-tests for each coefficient
estimate in a fitted linear regression model, using a sandwich estimator
for the standard errors and (optionally) a small sample correction for
the p-value. Available small-sample corrections include Satterthwaite
approximation or a saddlepoint approximation. Coefficients can be tested
against non-zero null values by specifying `null_constants`.

## Usage

``` r
coef_test(
  obj,
  vcov,
  test = "Satterthwaite",
  alternative = c("two-sided", "greater", "less"),
  coefs = "All",
  null_constants = 0,
  p_values = TRUE,
  ...
)
```

## Arguments

- obj:

  Fitted model for which to calculate t-tests.

- vcov:

  Variance covariance matrix estimated using `vcovCR` or a character
  string specifying which small-sample adjustment should be used to
  calculate the variance-covariance.

- test:

  Character vector specifying which small-sample corrections to
  calculate. `"z"` returns a z test (i.e., using a standard normal
  reference distribution). `"naive-t"` returns a t test with `m - 1`
  degrees of freedom, where `m` is the number of unique clusters.
  `"naive-tp"` returns a t test with `m - p` degrees of freedom, where
  `p` is the number of regression coefficients in `obj`.
  `"Satterthwaite"` returns a Satterthwaite correction. `"saddlepoint"`
  returns a saddlepoint correction. Default is `"Satterthwaite"`.

- alternative:

  Character string specifying the alternative hypothesis, with options
  "two-sided" (the default), "greater" or "less".

- coefs:

  Character, integer, or logical vector specifying which coefficients
  should be tested. The default value `"All"` will test all estimated
  coefficients.

- null_constants:

  vector of null values for each coefficient to test. Must have length
  equal to the number of coefficients specified in `coefs`. Default is
  `0`, in which case the null values are taken to be zero.

- p_values:

  Logical indicating whether to report p-values. The default value is
  `TRUE`.

- ...:

  Further arguments passed to
  [`vcovCR`](http://jepusto.github.io/clubSandwich/reference/vcovCR.md),
  which are only needed if `vcov` is a character string.

## Value

A data frame containing estimated regression coefficients, standard
errors, specified values of null hypotheses, and test results. For the
Satterthwaite approximation, degrees of freedom and a p-value are
reported. For the saddlepoint approximation, the saddlepoint and a
p-value are reported.

## See also

[`vcovCR`](http://jepusto.github.io/clubSandwich/reference/vcovCR.md)

## Examples

``` r

data("ChickWeight", package = "datasets")
lm_fit <- lm(weight ~ Diet  * Time, data = ChickWeight)
diet_index <- grepl("Diet.:Time", names(coef(lm_fit)))
coef_test(lm_fit, vcov = "CR2", cluster = ChickWeight$Chick, coefs = diet_index)
#> Alternative hypothesis: two-sided 
#>       Coef. Estimate   SE Null value t-stat d.f. (Satt) p-val (Satt) Sig.
#>  Diet2:Time     1.77 1.49          0   1.19        18.8       0.2497     
#>  Diet3:Time     4.58 1.35          0   3.39        18.8       0.0031   **
#>  Diet4:Time     2.87 1.01          0   2.85        18.3       0.0105    *

V_CR2 <- vcovCR(lm_fit, cluster = ChickWeight$Chick, type = "CR2")
coef_test(lm_fit, vcov = V_CR2, coefs = diet_index)
#> Alternative hypothesis: two-sided 
#>       Coef. Estimate   SE Null value t-stat d.f. (Satt) p-val (Satt) Sig.
#>  Diet2:Time     1.77 1.49          0   1.19        18.8       0.2497     
#>  Diet3:Time     4.58 1.35          0   3.39        18.8       0.0031   **
#>  Diet4:Time     2.87 1.01          0   2.85        18.3       0.0105    *

# non-inferiority test whether time-by-diet interaction effects are 2 or greater
coef_test(lm_fit, vcov = V_CR2, coefs = diet_index, null_constants = 2, alternative = "greater")
#> Alternative hypothesis: greater 
#>       Coef. Estimate   SE Null value t-stat d.f. (Satt) p-val (Satt) Sig.
#>  Diet2:Time     1.77 1.49          2 -0.156        18.8       0.5613     
#>  Diet3:Time     4.58 1.35          2  1.911        18.8       0.0357    *
#>  Diet4:Time     2.87 1.01          2  0.866        18.3       0.1990     
```
