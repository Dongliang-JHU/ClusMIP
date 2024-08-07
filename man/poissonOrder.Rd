% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gsf_modified.r
\name{poissonOrder}
\alias{poissonOrder}
\title{Estimate the Number of Components in a Poisson Finite Mixture Model.}
\usage{
poissonOrder(y, lambdas, K = NULL, ...)
}
\arguments{
\item{y}{Vector consisting of the data.}

\item{lambdas}{Vector of tuning parameter values.}

\item{K}{Upper bound on the true number of components.
If \code{K} is \code{NULL}, at least one of \code{theta} and
\code{pii} must be non-\code{NULL}, and K is inferred from their
number of columns.}

\item{...}{Additional control parameters. See the \strong{Details} section.}
}
\value{
An object with S3 classes \code{gsf} and \code{poissonGsf},
consisting of a list with the estimates produced for every tuning
parameter in \code{lambdas}.
}
\description{
Estimate the order of a finite mixture of Poisson distributions.
}
\details{

}
\examples{
 data(notices)

 # Run the GSF with the Adaptive Lasso penalty.
 set.seed(1)
 out <- poissonOrder(notices, lambdas=c(0, .001, .01, .5, 1, 2),
                     K=12, penalty="ADAPTIVE-LASSO", maxMem=1000)
 # Requires ggplot2.
 plot(out, eta=FALSE)

 # Run the GSF with the SCAD penalty.
 set.seed(1)
 out <- poissonOrder(notices, lambdas=c(.00005, .001, .005, .01, .1, .5, 1, 2, 5),
                     K=12, uBound=0.01, penalty="SCAD", maxMem=1000)
 plot(out, eta=FALSE)

 # Select a tuning parameter using the BIC.
 bicTuning(notices, out)
}
\references{
Manole, T., Khalili, A. 2019. "Estimating the Number of Components in Finite Mixture Models
via the Group-Sort-Fuse Procedure".
}
