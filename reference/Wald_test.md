# Test parameter constraints in a fitted linear regression model

`Wald_test` reports Wald-type tests of linear contrasts from a fitted
linear regression model, using a sandwich estimator for the
variance-covariance matrix and a small sample correction for the
p-value. Several different small-sample corrections are available.

## Usage

``` r
Wald_test(
  obj,
  constraints,
  vcov,
  null_constant = 0,
  test = "HTZ",
  tidy = FALSE,
  adjustment_method = "none",
  ...
)
```

## Arguments

- obj:

  Fitted model for which to calculate Wald tests.

- constraints:

  constraint or list of multiple constraints to test. See details and
  examples.

- vcov:

  Variance covariance matrix estimated using `vcovCR` or a character
  string specifying which small-sample adjustment should be used to
  calculate the variance-covariance.

- null_constant:

  vector of null values or list of such vectors for each set of
  constraints to test. For a single `constraint`, the null values must
  have length equal to the number of rows in the constraint. For lists
  of null values, each entry must have length equal to the number of
  rows in the corresponding entry of `constraints`. Default is `0`, in
  which case the null values are taken to be zero (for every entry, if
  `constraints` is a list).

- test:

  Character vector specifying which small-sample correction(s) to
  calculate. The following corrections are available: `"chi-sq"`,
  `"Naive-F"`, `"Naive-Fp"`, `"HTA"`, `"HTB"`, `"HTZ"`, `"EDF"`,
  `"EDT"`. Default is `"HTZ"`.

- tidy:

  Logical value controlling whether to tidy the test results. If
  `constraints` is a list with multiple constraints, the result will be
  coerced into a data frame when `tidy = TRUE`.

- adjustment_method:

  A character string indicating a multiple comparisons correction to
  apply to p-values in instances where multiple tests are run. Possible
  options are from
  [`p.adjust.methods`](https://rdrr.io/r/stats/p.adjust.html), which is
  passed to [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html) to
  correct p-values for multiple comparisons. Defaults to `"none"`.

- ...:

  Further arguments passed to
  [`vcovCR`](http://jepusto.github.io/clubSandwich/reference/vcovCR.md),
  which are only needed if `vcov` is a character string.

## Value

A list of test results.

## Details

Constraints can be specified directly as q X p matrices or indirectly
through
[`constrain_equal`](http://jepusto.github.io/clubSandwich/reference/constraint_matrices.md),
[`constrain_zero`](http://jepusto.github.io/clubSandwich/reference/constraint_matrices.md),
or
[`constrain_pairwise`](http://jepusto.github.io/clubSandwich/reference/constraint_matrices.md).
By default, each constraint will be tested against the null hypothesis
that it equal to a zero vector. Non-zero values for null-hypotheses can
be specified using the `null_constant` argument.

## See also

[`vcovCR`](http://jepusto.github.io/clubSandwich/reference/vcovCR.md),
[`constrain_equal`](http://jepusto.github.io/clubSandwich/reference/constraint_matrices.md),
[`constrain_zero`](http://jepusto.github.io/clubSandwich/reference/constraint_matrices.md),
[`constrain_pairwise`](http://jepusto.github.io/clubSandwich/reference/constraint_matrices.md)

## Examples

``` r


if (requireNamespace("carData", quietly = TRUE)) withAutoprint({

data(Duncan, package = "carData")
Duncan$cluster <- sample(LETTERS[1:8], size = nrow(Duncan), replace = TRUE)

Duncan_fit <- lm(prestige ~ 0 + type + income + type:income + type:education, data=Duncan)
# Note that type:income terms are interactions because main effect of income is included
# but type:education terms are separate slopes for each unique level of type

# Test equality of intercepts
Wald_test(Duncan_fit,
          constraints = constrain_equal(1:3),
          vcov = "CR2", cluster = Duncan$cluster)

# Test equality of type-by-education slopes
Wald_test(Duncan_fit,
          constraints = constrain_equal(":education", reg_ex = TRUE),
          vcov = "CR2", cluster = Duncan$cluster)

# Pairwise comparisons of type-by-education slopes
Wald_test(Duncan_fit,
          constraints = constrain_pairwise(":education", reg_ex = TRUE),
          vcov = "CR2", cluster = Duncan$cluster)

# Test type-by-income interactions
Wald_test(Duncan_fit,
          constraints = constrain_zero(":income", reg_ex = TRUE),
          vcov = "CR2", cluster = Duncan$cluster)

# Pairwise comparisons of type-by-income interactions
Wald_test(Duncan_fit,
          constraints = constrain_pairwise(":income", reg_ex = TRUE, with_zero = TRUE),
          vcov = "CR2", cluster = Duncan$cluster)

# Pairwise comparisons of type-by-education slopes, with two tests and 
# multiple comparisons p-value adjustment
Wald_test(Duncan_fit,
          constraints = constrain_pairwise(":education", reg_ex = TRUE),
          vcov = "CR2",
          cluster = Duncan$cluster,
          test = c("HTZ","chi-sq"),
          adjustment_method = "holm")

})
#> > data(Duncan, package = "carData")
#> > Duncan$cluster <- sample(LETTERS[1:8], size = nrow(Duncan), replace = TRUE)
#> > Duncan_fit <- lm(prestige ~ 0 + type + income + type:income + type:education, 
#> +     data = Duncan)
#> > Wald_test(Duncan_fit, constraints = constrain_equal(1:3), vcov = "CR2", 
#> +     cluster = Duncan$cluster)
#>  test Fstat df_num df_denom p_val sig
#>   HTZ  1.47      2      1.7 0.426    
#> > Wald_test(Duncan_fit, constraints = constrain_equal(":education", reg_ex = TRUE), 
#> +     vcov = "CR2", cluster = Duncan$cluster)
#>  test  Fstat df_num df_denom p_val sig
#>   HTZ 0.0648      2     2.72 0.939    
#> > Wald_test(Duncan_fit, constraints = constrain_pairwise(":education", reg_ex = TRUE), 
#> +     vcov = "CR2", cluster = Duncan$cluster)
#> $`typeprof:education - typebc:education`
#>  test   Fstat df_num df_denom p_val sig
#>   HTZ 0.00289      1      4.3  0.96    
#> 
#> $`typewc:education - typebc:education`
#>  test Fstat df_num df_denom p_val sig
#>   HTZ 0.136      1      4.2 0.731    
#> 
#> $`typewc:education - typeprof:education`
#>  test  Fstat df_num df_denom p_val sig
#>   HTZ 0.0918      1     2.63 0.784    
#> 
#> > Wald_test(Duncan_fit, constraints = constrain_zero(":income", reg_ex = TRUE), 
#> +     vcov = "CR2", cluster = Duncan$cluster)
#>  test Fstat df_num df_denom p_val sig
#>   HTZ  1.07      2     1.93 0.486    
#> > Wald_test(Duncan_fit, constraints = constrain_pairwise(":income", reg_ex = TRUE, 
#> +     with_zero = TRUE), vcov = "CR2", cluster = Duncan$cluster)
#> $`typeprof:income`
#>  test Fstat df_num df_denom p_val sig
#>   HTZ   1.5      1     2.77 0.314    
#> 
#> $`typewc:income`
#>  test Fstat df_num df_denom p_val sig
#>   HTZ  2.15      1     2.28 0.265    
#> 
#> $`typewc:income - typeprof:income`
#>  test    Fstat df_num df_denom p_val sig
#>   HTZ 0.000588      1     2.41 0.983    
#> 
#> > Wald_test(Duncan_fit, constraints = constrain_pairwise(":education", reg_ex = TRUE), 
#> +     vcov = "CR2", cluster = Duncan$cluster, test = c("HTZ", "chi-sq"), adjustment_method = "holm")
#> $`typeprof:education - typebc:education`
#>    test   Fstat df_num df_denom p_val sig
#>  chi-sq 0.00289      1      Inf     1    
#>     HTZ 0.00289      1      4.3     1    
#> 
#> $`typewc:education - typebc:education`
#>    test Fstat df_num df_denom p_val sig
#>  chi-sq 0.136      1      Inf     1    
#>     HTZ 0.136      1      4.2     1    
#> 
#> $`typewc:education - typeprof:education`
#>    test  Fstat df_num df_denom p_val sig
#>  chi-sq 0.0918      1      Inf     1    
#>     HTZ 0.0918      1     2.63     1    
#> 
```
