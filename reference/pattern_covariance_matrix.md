# Impute a patterned block-diagonal covariance matrix

\`r lifecycle::badge("superseded")\`

This function is superseded by the
[`vcalc`](https://wviechtb.github.io/metafor/reference/vcalc.html)
provided by the `metafor` package. Compared to
`pattern_covariance_matrix`,
[`vcalc`](https://wviechtb.github.io/metafor/reference/vcalc.html)
provides many further features, includes a `data` argument, and uses
syntax that is consistent with other functions in `metafor`.

`pattern_covariance_matrix` calculates a block-diagonal covariance
matrix, given the marginal variances, the block structure, and an
assumed correlation structure defined by a patterned correlation matrix.

## Usage

``` r
pattern_covariance_matrix(
  vi,
  cluster,
  pattern_level,
  r_pattern,
  r,
  smooth_vi = FALSE,
  subgroup = NULL,
  return_list = identical(as.factor(cluster), sort(as.factor(cluster))),
  check_PD = TRUE
)
```

## Arguments

- vi:

  Vector of variances

- cluster:

  Vector indicating which effects belong to the same cluster. Effects
  with the same value of \`cluster\` will be treated as correlated.

- pattern_level:

  Vector of categories for each effect size, used to determine which
  entry of the pattern matrix will be used to impute a correlation.

- r_pattern:

  Patterned correlation matrix with row and column names corresponding
  to the levels of `pattern`.

- r:

  Vector or numeric value of assumed constant correlation(s) between
  effect size estimates from each study.

- smooth_vi:

  Logical indicating whether to smooth the marginal variances by taking
  the average `vi` within each cluster. Defaults to `FALSE`.

- subgroup:

  Vector of category labels describing sub-groups of effects. If
  non-null, effects that share the same category label and the same
  cluster will be treated as correlated, but effects with different
  category labels will be treated as uncorrelated, even if they come
  from the same cluster.

- return_list:

  Optional logical indicating whether to return a list of matrices (with
  one entry per block) or the full variance-covariance matrix.

- check_PD:

  Optional logical indicating whether to check whether each covariance
  matrix is positive definite. If `TRUE` (the default), the function
  will display a warning if any covariance matrix is not positive
  definite.

## Value

If `cluster` is appropriately sorted, then a list of matrices, with one
entry per cluster, will be returned by default. If `cluster` is out of
order, then the full variance-covariance matrix will be returned by
default. The output structure can be controlled with the optional
`return_list` argument.

## Details

A block-diagonal variance-covariance matrix (possibly represented as a
list of matrices) with a specified correlation structure, defined by a
patterned correlation matrix. Let \\v\_{ij}\\ denote the specified
variance for effect \\i\\ in cluster \\j\\ and \\C\_{hij}\\ be the
covariance between effects \\h\\ and \\i\\ in cluster \\j\\. Let
\\p\_{ij}\\ be the level of the pattern variable for effect \\i\\ in
cluster \\j\\, taking a value in \\1,...,C\\. A patterned correlation
matrix is defined as a set of correlations between pairs of effects
taking each possible combination of patterns. Formally, let \\r\_{cd}\\
be the correlation between effects in categories \\c\\ and \\d\\,
respectively, where \\r\_{cd} = r\_{dc}\\. Then the covariance between
effects \\h\\ and \\i\\ in cluster \\j\\ is taken to be \$\$C\_{hij} =
\sqrt{v\_{hj} v\_{ij}} \times r\_{p\_{hj} p\_{ij}}.\$\$

Correlations between effect sizes within the same category are defined
by the diagonal values of the pattern matrix, which may take values less
than one.

Combinations of pattern levels that do not occur in the patterned
correlation matrix will be set equal to `r`.

If `smooth_vi = TRUE`, then all of the variances within cluster \\j\\
will be set equal to the average variance of cluster \\j\\, i.e.,
\$\$v'\_{ij} = \frac{1}{n_j} \sum\_{i=1}^{n_j} v\_{ij}\$\$ for
\\i=1,...,n_j\\ and \\j=1,...,k\\.

## Examples

``` r

pkgs_available <- 
  requireNamespace("metafor", quietly = TRUE) & 
  requireNamespace("robumeta", quietly = TRUE)
  
if (pkgs_available) {
library(metafor)

data(oswald2013, package = "robumeta")
dat <- escalc(data = oswald2013, measure = "ZCOR", ri = R, ni = N)
subset_ids <- unique(dat$Study)[1:20]
dat <- subset(dat, Study %in% subset_ids)

# make a patterned correlation matrix 

p_levels <- levels(dat$Crit.Cat)
r_pattern <- 0.7^as.matrix(dist(1:length(p_levels)))
diag(r_pattern) <- seq(0.75, 0.95, length.out = 6)
rownames(r_pattern) <- colnames(r_pattern) <- p_levels

# impute the covariance matrix using patterned correlations
V_list <- pattern_covariance_matrix(vi = dat$vi, 
                                    cluster = dat$Study, 
                                    pattern_level = dat$Crit.Cat,
                                    r_pattern = r_pattern,
                                    smooth_vi = TRUE)
                                    
# fit a model using imputed covariance matrix

MVFE <- rma.mv(yi ~ 0 + Crit.Cat, V = V_list, 
               random = ~ Crit.Cat | Study,
               data = dat)
               
conf_int(MVFE, vcov = "CR2")

}
#>                           Coef. Estimate       SE d.f. Lower 95% CI
#>          Crit.CatBrain Activity   0.1425 4.49e-01 2.84     -1.33431
#>  Crit.CatInterpersonal Behavior   0.0167 3.65e-03 1.63     -0.00307
#>           Crit.CatMicrobehavior   0.1308 3.43e-02 3.64      0.03178
#>       Crit.CatPerson Perception   0.1544 1.94e-02 8.08      0.10963
#>       Crit.CatPolicy Preference   0.0864 5.10e-02 1.98     -0.13561
#>           Crit.CatResponse Time   0.2988 1.42e-13 1.00      0.29882
#>  Upper 95% CI
#>        1.6193
#>        0.0364
#>        0.2298
#>        0.1991
#>        0.3084
#>        0.2988
```
