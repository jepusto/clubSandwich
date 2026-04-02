# Calculate confidence intervals and p-values for linear contrasts of regression coefficients in a fitted model

`linear_contrast` reports confidence intervals and (optionally) p-values
for linear contrasts of regression coefficients from a fitted model,
using a sandwich estimator for the standard errors and (optionally) a
small sample correction for the critical values. The default
small-sample correction is based on a Satterthwaite approximation.

## Usage

``` r
linear_contrast(
  obj,
  vcov,
  contrasts,
  level = 0.95,
  test = "Satterthwaite",
  ...,
  p_values = FALSE,
  adjustment_method = "none"
)
```

## Arguments

- obj:

  Fitted model for which to calculate confidence intervals.

- vcov:

  Variance covariance matrix estimated using `vcovCR` or a character
  string specifying which small-sample adjustment should be used to
  calculate the variance-covariance.

- contrasts:

  A contrast matrix, or a list of multiple contrast matrices to test.
  See details and examples.

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
  `"saddlepoint"` is not currently supported in
  [`conf_int()`](http://jepusto.github.io/clubSandwich/reference/conf_int.md)
  because saddlepoint confidence intervals do not have a closed-form
  solution.

- ...:

  Further arguments passed to
  [`vcovCR`](http://jepusto.github.io/clubSandwich/reference/vcovCR.md),
  which are only needed if `vcov` is a character string.

- p_values:

  Logical indicating whether to report p-values. The default value is
  `FALSE`.

- adjustment_method:

  Correction method, a
  [`character`](https://rdrr.io/r/base/character.html) string from
  `p.adjust.methods`, which is passed to
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html) to correct
  p-values in the case of multiple comparisons. Defaults to `"none"`.
  Any value given will be ignored if `p_values = FALSE`.

## Value

A data frame containing estimated contrasts, standard errors, confidence
intervals, and (optionally) p-values.

## Details

Constraints can be specified directly as q X p matrices or indirectly
through
[`constrain_pairwise`](http://jepusto.github.io/clubSandwich/reference/constraint_matrices.md),
[`constrain_equal`](http://jepusto.github.io/clubSandwich/reference/constraint_matrices.md),
or
[`constrain_zero`](http://jepusto.github.io/clubSandwich/reference/constraint_matrices.md).

## See also

[`vcovCR`](http://jepusto.github.io/clubSandwich/reference/vcovCR.md)

## Examples

``` r
data("ChickWeight", package = "datasets")
lm_fit <- lm(weight ~ 0 + Diet + Time:Diet, data = ChickWeight)

# Pairwise comparisons of diet-by-time slopes
linear_contrast(lm_fit, vcov = "CR2", cluster = ChickWeight$Chick, 
                contrasts = constrain_pairwise("Diet.:Time", reg_ex = TRUE))
#>                    Coef. Estimate   SE d.f. Lower 95% CI Upper 95% CI
#>  Diet2:Time - Diet1:Time     1.77 1.49 18.8       -1.349         4.88
#>  Diet3:Time - Diet1:Time     4.58 1.35 18.8        1.751         7.41
#>  Diet4:Time - Diet1:Time     2.87 1.01 18.3        0.757         4.99
#>  Diet3:Time - Diet2:Time     2.81 1.70 18.0       -0.756         6.38
#>  Diet4:Time - Diet2:Time     1.11 1.44 17.9       -1.925         4.14
#>  Diet4:Time - Diet3:Time    -1.71 1.30 17.9       -4.441         1.02


if (requireNamespace("carData", quietly = TRUE)) withAutoprint({

  data(Duncan, package = "carData")
  Duncan$cluster <- sample(LETTERS[1:8], size = nrow(Duncan), replace = TRUE)

  Duncan_fit <- lm(prestige ~ 0 + type + income + type:income + type:education, data=Duncan)
  # Note that type:income terms are interactions because main effect of income is included
  # but type:education terms are separate slopes for each unique level of type

  # Pairwise comparisons of type-by-education slopes
  linear_contrast(Duncan_fit, vcov = "CR2", cluster = Duncan$cluster,
                  contrasts = constrain_pairwise(":education", reg_ex = TRUE),
                  test = "Satterthwaite")
 
  # Repeat the above, but with p-values displayed and multiple comparisons adjustment applied.
  linear_contrast(Duncan_fit, vcov = "CR2", cluster = Duncan$cluster,
                  contrasts = constrain_pairwise(":education", reg_ex = TRUE),
                  test = "Satterthwaite", p_values = TRUE, adjustment_method = "hochberg")
 
  # Pairwise comparisons of type-by-income interactions
  linear_contrast(Duncan_fit, vcov = "CR2", cluster = Duncan$cluster,
                  contrasts = constrain_pairwise(":income", reg_ex = TRUE, with_zero = TRUE),
                  test = "Satterthwaite")
                  
})
#> > data(Duncan, package = "carData")
#> > Duncan$cluster <- sample(LETTERS[1:8], size = nrow(Duncan), replace = TRUE)
#> > Duncan_fit <- lm(prestige ~ 0 + type + income + type:income + type:education, 
#> +     data = Duncan)
#> > linear_contrast(Duncan_fit, vcov = "CR2", cluster = Duncan$cluster, contrasts = constrain_pairwise(":education", 
#> +     reg_ex = TRUE), test = "Satterthwaite")
#>                                  Coef. Estimate    SE d.f. Lower 95% CI
#>  typeprof:education - typebc:education   0.0186 0.390 4.42       -1.025
#>    typewc:education - typebc:education   0.1068 0.319 3.44       -0.838
#>  typewc:education - typeprof:education   0.0882 0.365 2.15       -1.382
#>  Upper 95% CI
#>          1.06
#>          1.05
#>          1.56
#> > linear_contrast(Duncan_fit, vcov = "CR2", cluster = Duncan$cluster, contrasts = constrain_pairwise(":education", 
#> +     reg_ex = TRUE), test = "Satterthwaite", p_values = TRUE, adjustment_method = "hochberg")
#>                                  Coef. Estimate    SE d.f. Lower 95% CI
#>  typeprof:education - typebc:education   0.0186 0.390 4.42       -1.025
#>    typewc:education - typebc:education   0.1068 0.319 3.44       -0.838
#>  typewc:education - typeprof:education   0.0882 0.365 2.15       -1.382
#>  Upper 95% CI p-value Sig.
#>          1.06   0.964     
#>          1.05   0.964     
#>          1.56   0.964     
#> > linear_contrast(Duncan_fit, vcov = "CR2", cluster = Duncan$cluster, contrasts = constrain_pairwise(":income", 
#> +     reg_ex = TRUE, with_zero = TRUE), test = "Satterthwaite")
#>                            Coef. Estimate    SE d.f. Lower 95% CI Upper 95% CI
#>                  typeprof:income -0.36914 0.367 3.56        -1.44        0.703
#>                    typewc:income -0.36031 0.296 2.23        -1.52        0.797
#>  typewc:income - typeprof:income  0.00883 0.363 2.64        -1.24        1.257
```
