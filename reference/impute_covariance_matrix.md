# Impute a block-diagonal covariance matrix

\`r lifecycle::badge("superseded")\`

This function is superseded by the
[`vcalc`](https://wviechtb.github.io/metafor/reference/vcalc.html)
provided by the `metafor` package. Compared to
`impute_covariance_matrix`,
[`vcalc`](https://wviechtb.github.io/metafor/reference/vcalc.html)
provides many further features, includes a `data` argument, and uses
syntax that is consistent with other functions in `metafor`.

`impute_covariance_matrix` calculates a block-diagonal covariance
matrix, given the marginal variances, the block structure, and an
assumed correlation structure. Can be used to create compound-symmetric
structures, AR(1) auto-correlated structures, or combinations thereof.

## Usage

``` r
impute_covariance_matrix(
  vi,
  cluster,
  r,
  ti,
  ar1,
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

- r:

  Vector or numeric value of assumed constant correlation(s) between
  effect size estimates from each study.

- ti:

  Vector of time-points describing temporal spacing of effects, for use
  with auto-regressive correlation structures.

- ar1:

  Vector or numeric value of assumed AR(1) auto-correlation(s) between
  effect size estimates from each study. If specified, then `ti`
  argument must be specified.

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
list of matrices) with a specified structure. The structure depends on
whether the `r` argument, `ar1` argument, or both arguments are
specified. Let \\v\_{ij}\\ denote the specified variance for effect
\\i\\ in cluster \\j\\ and \\C\_{hij}\\ be the covariance between
effects \\h\\ and \\i\\ in cluster \\j\\.

- If only `r` is specified, each block of the variance-covariance matrix
  will have a constant (compound symmetric) correlation, so that
  \$\$C\_{hij} = r_j \sqrt{v\_{hj} v\_{ij},}\$\$ where \\r_j\\ is the
  specified correlation for cluster \\j\\. If only a single value is
  given in `r`, then it will be used for every cluster.

- If only `ar1` is specified, each block of the variance-covariance
  matrix will have an AR(1) auto-correlation structure, so that
  \$\$C\_{hij} = \phi_j^{\|t\_{hj}- t\_{ij}\|} \sqrt{v\_{hj}
  v\_{ij},}\$\$ where \\\phi_j\\ is the specified auto-correlation for
  cluster \\j\\ and \\t\_{hj}\\ and \\t\_{ij}\\ are specified
  time-points corresponding to effects \\h\\ and \\i\\ in cluster \\j\\.
  If only a single value is given in `ar1`, then it will be used for
  every cluster.

- If both `r` and `ar1` are specified, each block of the
  variance-covariance matrix will have combination of compound symmetric
  and an AR(1) auto-correlation structures, so that \$\$C\_{hij} =
  \left\[r_j + (1 - r_j)\phi_j^{\|t\_{hj} - t\_{ij}\|}\right\]
  \sqrt{v\_{hj} v\_{ij},}\$\$ where \\r_j\\ is the specified constant
  correlation for cluster \\j\\, \\\phi_j\\ is the specified
  auto-correlation for cluster \\j\\ and \\t\_{hj}\\ and \\t\_{ij}\\ are
  specified time-points corresponding to effects \\h\\ and \\i\\ in
  cluster \\j\\. If only single values are given in `r` or `ar1`, they
  will be used for every cluster.

If `smooth_vi = TRUE`, then all of the variances within cluster \\j\\
will be set equal to the average variance of cluster \\j\\, i.e.,
\$\$v'\_{ij} = \frac{1}{n_j} \sum\_{i=1}^{n_j} v\_{ij}\$\$ for
\\i=1,...,n_j\\ and \\j=1,...,k\\.

## Examples

``` r

if (requireNamespace("metafor", quietly = TRUE)) {

library(metafor)

# Constant correlation
data(SATcoaching)
V_list <- impute_covariance_matrix(vi = SATcoaching$V, cluster = SATcoaching$study, r = 0.66)
MVFE <- rma.mv(d ~ 0 + test, V = V_list, data = SATcoaching)
conf_int(MVFE, vcov = "CR2", cluster = SATcoaching$study)

}
#> Warning: `impute_covariance_matrix()` was deprecated in clubSandwich 0.5.11.
#> ℹ Please use `metafor::vcalc()` instead.
#>       Coef. Estimate     SE d.f. Lower 95% CI Upper 95% CI
#>    testMath    0.132 0.0376 13.2       0.0506        0.213
#>  testVerbal    0.122 0.0250 16.9       0.0686        0.174
```
