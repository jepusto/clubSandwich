% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plm.R
\name{vcovCR.plm}
\alias{vcovCR.plm}
\title{Cluster-robust variance-covariance matrix for a plm object.}
\usage{
\method{vcovCR}{plm}(
  obj,
  cluster,
  type,
  target,
  inverse_var,
  form = "sandwich",
  ignore_FE = FALSE,
  ...
)
}
\arguments{
\item{obj}{Fitted model for which to calculate the variance-covariance matrix}

\item{cluster}{Optional character string, expression, or vector indicating
which observations belong to the same cluster. For fixed-effect models that
include individual effects or time effects (but not both), the cluster will
be taken equal to the included fixed effects if not otherwise specified.
Clustering on individuals can also be obtained by specifying the name of
the individual index (e.g., \code{cluster = "state"}) or \code{cluster =
"individual"}; clustering on time periods can be obtained by specifying the
name of the time index (e.g., \code{cluster = "year"}) or \code{cluster =
"time"}; if a group index is specified, clustering on groups (in which
individuals are nested) can be obtained by specifying the name of the group
index or \code{cluster = "group"}. For random-effects models, the cluster
will be taken equal to the included random effect identifier if not
otherwise specified.}

\item{type}{Character string specifying which small-sample adjustment should 
be used, with available options \code{"CR0"}, \code{"CR1"}, \code{"CR1p"},
\code{"CR1S"}, \code{"CR2"}, or \code{"CR3"}. See "Details" section of 
\code{\link{vcovCR}} for further information.}

\item{target}{Optional matrix or vector describing the working
variance-covariance model used to calculate the \code{CR2} and \code{CR4}
adjustment matrices. By default, the target is taken to be an identity
matrix for fixed effect models or the estimated compound-symmetric
covariance matrix for random effects models.}

\item{inverse_var}{Optional logical indicating whether the weights used in 
fitting the model are inverse-variance. If not specified, \code{vcovCR} 
will attempt to infer a value.}

\item{form}{Controls the form of the returned matrix. The default 
\code{"sandwich"} will return the sandwich variance-covariance matrix. 
Alternately, setting \code{form = "meat"} will return only the meat of the 
sandwich and setting \code{form = B}, where \code{B} is a matrix of 
appropriate dimension, will return the sandwich variance-covariance matrix 
calculated using \code{B} as the bread. \code{form = "estfun"} will return the 
(appropriately scaled) estimating function, the transposed crossproduct of 
which is equal to the sandwich variance-covariance matrix.}

\item{ignore_FE}{Optional logical controlling whether fixed effects are
ignored when calculating small-sample adjustments in models where fixed
effects are estimated through absorption.}

\item{...}{Additional arguments available for some classes of objects.}
}
\value{
An object of class \code{c("vcovCR","clubSandwich")}, which consists
  of a matrix of the estimated variance of and covariances between the
  regression coefficient estimates.
}
\description{
\code{vcovCR} returns a sandwich estimate of the variance-covariance matrix
of a set of regression coefficient estimates from a \code{\link[plm]{plm}}
object.
}
\examples{

if (requireNamespace("plm", quietly = TRUE)) withAutoprint({

  library(plm)
  # fixed effects
  data("Produc", package = "plm")
  plm_FE <- plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp,
                data = Produc, index = c("state","year","region"),
                effect = "individual", model = "within")
  vcovCR(plm_FE, type="CR2")
  vcovCR(plm_FE, type = "CR2", cluster = Produc$region) # clustering on region
  
  # random effects
  plm_RE <- update(plm_FE, model = "random")
  vcovCR(plm_RE, type = "CR2")
  vcovCR(plm_RE, type = "CR2", cluster = Produc$region) # clustering on region
  
  # nested random effects
  plm_nested <- update(plm_FE, effect = "nested", model = "random")
  vcovCR(plm_nested, type = "CR2") # clustering on region
})

pkgs_available <- requireNamespace("plm", quietly = TRUE) & requireNamespace("AER", quietly = TRUE)

if (pkgs_available) withAutoprint({
  # first differencing
  data(Fatalities, package = "AER")
  Fatalities <- within(Fatalities, {
    frate <- 10000 * fatal / pop
    drinkagec <- cut(drinkage, breaks = 18:22, include.lowest = TRUE, right = FALSE)
    drinkagec <- relevel(drinkagec, ref = 4)
  })

  plm_FD <- plm(frate ~ beertax + drinkagec + miles + unemp + log(income),
                data = Fatalities, index = c("state", "year"),
                model = "fd")
  vcovHC(plm_FD, method="arellano", type = "sss", cluster = "group")
  vcovCR(plm_FD, type = "CR1S")
  vcovCR(plm_FD, type = "CR2")
  
})

}
\seealso{
\code{\link{vcovCR}}
}
