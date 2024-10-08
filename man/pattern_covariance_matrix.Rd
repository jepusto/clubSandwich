% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rma-mv.R
\name{pattern_covariance_matrix}
\alias{pattern_covariance_matrix}
\title{Impute a patterned block-diagonal covariance matrix}
\usage{
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
}
\arguments{
\item{vi}{Vector of variances}

\item{cluster}{Vector indicating which effects belong to the same cluster.
Effects with the same value of `cluster` will be treated as correlated.}

\item{pattern_level}{Vector of categories for each effect size, used to
determine which entry of the pattern matrix will be used to impute a
correlation.}

\item{r_pattern}{Patterned correlation matrix with row and column names
corresponding to the levels of \code{pattern}.}

\item{r}{Vector or numeric value of assumed constant correlation(s) between
effect size estimates from each study.}

\item{smooth_vi}{Logical indicating whether to smooth the marginal variances
by taking the average \code{vi} within each cluster. Defaults to
\code{FALSE}.}

\item{subgroup}{Vector of category labels describing sub-groups of effects.
If non-null, effects that share the same category label and the same
cluster will be treated as correlated, but effects with different category
labels will be treated as uncorrelated, even if they come from the same
cluster.}

\item{return_list}{Optional logical indicating whether to return a list of
matrices (with one entry per block) or the full variance-covariance matrix.}

\item{check_PD}{Optional logical indicating whether to check whether each
covariance matrix is positive definite. If \code{TRUE} (the default), the
function will display a warning if any covariance matrix is not positive
definite.}
}
\value{
If \code{cluster} is appropriately sorted, then a list of matrices,
  with one entry per cluster, will be returned by default. If \code{cluster}
  is out of order, then the full variance-covariance matrix will be returned
  by default. The output structure can be controlled with the optional
  \code{return_list} argument.
}
\description{
`r lifecycle::badge("superseded")`

  This function is superseded by the \code{\link[metafor]{vcalc}} provided by
  the \code{metafor} package. Compared to \code{pattern_covariance_matrix},
  \code{\link[metafor]{vcalc}} provides many further features, includes a
  \code{data} argument, and uses syntax that is consistent with other
  functions in \code{metafor}.

\code{pattern_covariance_matrix} calculates a
  block-diagonal covariance matrix, given the marginal variances, the block
  structure, and an assumed correlation structure defined by a patterned
  correlation matrix.
}
\details{
A block-diagonal variance-covariance matrix (possibly represented as
  a list of matrices) with a specified correlation structure, defined by a
  patterned correlation matrix. Let \eqn{v_{ij}}{v-ij} denote the specified
  variance for effect \eqn{i}{i} in cluster \eqn{j}{j} and
  \eqn{C_{hij}}{C-hij} be the covariance between effects \eqn{h}{h} and
  \eqn{i}{i} in cluster \eqn{j}{j}. Let \eqn{p_{ij}}{p-ij} be the level
  of the pattern variable for effect \eqn{i}{i} in cluster \eqn{j}{j},
  taking a value in \eqn{1,...,C}{1,...,C}. A patterned correlation matrix
  is defined as a set of correlations between pairs of effects taking each
  possible combination of patterns. Formally, let \eqn{r_{cd}}{r-cd} be the
  correlation between effects in categories \eqn{c}{c} and \eqn{d}{d},
  respectively, where \eqn{r_{cd} = r_{dc}}{r-cd = r-dc}. Then the
  covariance between effects \eqn{h}{h} and \eqn{i}{i} in cluster
  \eqn{j}{j} is taken to be \deqn{C_{hij} = \sqrt{v_{hj} v_{ij}} \times
  r_{p_{hj} p_{ij}}.}{C-hij = sqrt(v-hj v-ij) * r[p-hj, p-ij].} 
  
  Correlations between effect sizes within the same category are defined by the diagonal
  values of the pattern matrix, which may take values less than one. 
  
  Combinations of pattern levels that do not occur in the patterned correlation matrix will be set equal to \code{r}.
  
  If \code{smooth_vi = TRUE}, then all of the variances within cluster
  \eqn{j}{j} will be set equal to the average variance of cluster
  \eqn{j}{j}, i.e., \deqn{v'_{ij} = \frac{1}{n_j} \sum_{i=1}^{n_j}
  v_{ij}}{v-ij' = (v-1j + ... + v-nj,j) / n-j} for
  \eqn{i=1,...,n_j}{i=1,...,n-j} and \eqn{j=1,...,k}{j=1,...,k}.
}
\examples{

pkgs_available <- 
  requireNamespace("metafor", quietly = TRUE) & 
  requireNamespace("robumeta", quietly = TRUE)
  
if (pkgs_available) {
library(metafor)

data(oswald2013, package = "robumeta")
dat <- escalc(data = oswald2013, measure = "ZCOR", ri = R, ni = N)
subset_ids <- unique(dat$Study)[1:20]
dat <- subset(dat, Study \%in\% subset_ids)

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

}
