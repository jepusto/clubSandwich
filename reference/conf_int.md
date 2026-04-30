# Calculate confidence intervals for all or selected regression coefficients in a fitted model

`conf_int` reports confidence intervals for each coefficient estimate in
a fitted linear regression model, using a sandwich estimator for the
standard errors and a small sample correction for the critical values.
The small-sample correction is based on a Satterthwaite approximation.

## Usage

``` r
conf_int(
  obj,
  vcov,
  level = 0.95,
  test = "Satterthwaite",
  coefs = "All",
  ...,
  p_values = FALSE
)
```

## Arguments

- obj:

  Fitted model for which to calculate confidence intervals.

- vcov:

  Variance covariance matrix estimated using `vcovCR` or a character
  string specifying which small-sample adjustment should be used to
  calculate the variance-covariance.

- level:

  Desired coverage level for confidence intervals.

- test:

  Character vector specifying which small-sample corrections to
  calculate. `"z"` returns a z test (i.e., using a standard normal
  reference distribution). `"naive-t"` returns a t test with `m - 1`
  degrees of freedom, where `m` is the number of unique clusters.
  `"naive-tp"` returns a t test with `m - p` degrees of freedom, where
  `p` is the number of regression coefficients in `obj`.
  `"Satterthwaite"` returns a Satterthwaite correction. Unlike in
  [`coef_test()`](http://jepusto.github.io/clubSandwich/reference/coef_test.md),
  `"saddlepoint"` is not currently supported in `conf_int()` because
  saddlepoint confidence intervals do not have a closed-form solution.

- coefs:

  Character, integer, or logical vector specifying which coefficients
  should be tested. The default value `"All"` will test all estimated
  coefficients.

- ...:

  Further arguments passed to
  [`vcovCR`](http://jepusto.github.io/clubSandwich/reference/vcovCR.md),
  which are only needed if `vcov` is a character string.

- p_values:

  Logical indicating whether to report p-values. The default value is
  `TRUE`.

## Value

A data frame containing estimated regression coefficients, standard
errors, confidence intervals, and (optionally) p-values.

## See also

[`vcovCR`](http://jepusto.github.io/clubSandwich/reference/vcovCR.md)

## Examples

``` r

data("ChickWeight", package = "datasets")
lm_fit <- lm(weight ~ Diet  * Time, data = ChickWeight)
diet_index <- grepl("Diet.:Time", names(coef(lm_fit)))
conf_int(lm_fit, vcov = "CR2", cluster = ChickWeight$Chick, coefs = diet_index)
#>       Coef. Estimate   SE d.f. Lower 95% CI Upper 95% CI
#>  Diet2:Time     1.77 1.49 18.8       -1.349         4.88
#>  Diet3:Time     4.58 1.35 18.8        1.751         7.41
#>  Diet4:Time     2.87 1.01 18.3        0.757         4.99

V_CR2 <- vcovCR(lm_fit, cluster = ChickWeight$Chick, type = "CR2")
conf_int(lm_fit, vcov = V_CR2, level = .99, coefs = diet_index)
#>       Coef. Estimate   SE d.f. Lower 99% CI Upper 99% CI
#>  Diet2:Time     1.77 1.49 18.8      -2.4946         6.03
#>  Diet3:Time     4.58 1.35 18.8       0.7115         8.45
#>  Diet4:Time     2.87 1.01 18.3      -0.0237         5.77
```
