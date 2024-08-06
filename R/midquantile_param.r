# Author: Dongliang Zhang (^)
#
# Supervisors: Professor Masoud Asgharian (*)
#              Professor Martin Lindquist (^)
#              Professor Mei-Cheng Wang (^)
#
# Current Affiliation: (*) Department of Mathematics and Statistics, McGill University, Montreal, Quebec, Canada
#                      (^) Department of Biostatistics, Johns Hopkins University, Baltimore, Maryland, USA
#
# Title of Project: High-dimensional Influential Diagnosis on Variable Selection: Revival of the Exchangeability
#
# R Script Purpose: Mid-quantile Computation
#
# Created on: October 24, 2023
#
# Modified on: October 24, 2023
#              October 25, 2023
#              October 26, 2023
#              October 28, 2023
#
# Dependent R Scripts:
#
# Functions: midcdf.param, midquantile.param

################################################################################
# Function: midcdf.param
#
# Args:
#
# x: tau vector
# p: number of predictors
# param.string: name of parametric approximating distribution
# param.val: MLE estimates for parameters
#
# Returns:
# a midcdf object

#' Title
#'
#' @param x tau vector
#' @param p number of predictors
#' @param param.string name of parametric approximating distribution
#' @param param.val MLE estimates for parameters
#'
#' @return a midcdf object
#' @export
#'
midcdf.param <- function(x, p, param.string, param.val){

  xo <- unique(x)

  if(param.string=="cmb"){
    pmf <- sapply(xo, function(x) d_cmb(x, m=p, p=param.val[1], nu=param.val[2], take_log=FALSE, normalize=TRUE))
    cdf.cmb <- function(x){
      summand <- c()
      for(i in 0:x){summand <- c(summand, d_cmb(i, m=p, p=param.val[1], nu=param.val[2], take_log=FALSE, normalize=TRUE))}
      return(sum(summand))}
    cdf <- sapply(xo, cdf.cmb)
  }else if(param.string=="cmp"){
    pmf <- dcmp(xo, lambda=param.val[1], nu=param.val[2], log=FALSE, control=NULL)
    cdf <- pcmp(xo, lambda=param.val[1], nu=param.val[2], control=NULL)
  }else if(param.string=="bb"){
    pmf <- dbb(xo, N=p, u=param.val[1], v=param.val[2], log=FALSE)
    cdf <- pbb(xo, N=p, u=param.val[1], v=param.val[2])
  }else if(param.string=="gp"){
    pmf <- dgenpois0(xo, theta=param.val[1], lambda=param.val[2])
    cdf <- pgenpois0(xo, theta=param.val[1], lambda=param.val[2])
  }else if(param.string=="mb"){
    pmf.mb <- as.matrix(mapply(function(psi,success) psi*dbinom(xo, size=p, prob=success), psi=param.val[1,], success=param.val[2,]))
    pmf <- rowSums(pmf.mb)
    cdf.mb <- function(x){
      pmf.mb.tmp <- as.matrix(mapply(function(psi,success) psi*dbinom(0:x, size=p, prob=success), psi=param.val[1,], success=param.val[2,]))
      return(sum(rowSums(pmf.mb.tmp)))
    }
    cdf <- sapply(xo,cdf.mb)
  }else if(param.string=="mp"){
    pmf <- dmixpois(xo, alpha=param.val[1,], lambda=param.val[2,], log=FALSE)
    cdf <- pmixpois(xo, alpha=param.val[1,], lambda=param.val[2,], lower.tail=TRUE, log.p=FALSE)
  }

  val <- list()
  val$call <- match.call()
  val$x <- xo
  val$y <- cdf - 0.5*pmf
  val$fn <- approxfun(val$x, val$y, method = "linear", rule = 1)
  val$data <- x
  class(val) <- "midcdf"
  return(val)

}

################################################################################
# Function: midquantile.param
#
# Args:
#
# x            : tau vector
# p            : number of predictors
# param.string : name of parametric approximating distribution
# param.val    : MLE estimates for parameters
# probs        : percentile
#
# Returns:
#
# a midquantile object
#
# Dependent R functions: midcdf.param

#' Title
#'
#' @param x tau vector
#' @param p number of predictors
#' @param param.string name of parametric approximating distribution
#' @param param.val MLE estimates for parameters
#' @param probs percentile
#'
#' @return a midquantile object
#' @export
#'
midquantile.param <- function(x, p, param.string, param.val, probs){

  Fn <- midcdf.param(x, p, param.string, param.val)
  Qn <- approxfun(Fn$y, Fn$x, method = "linear", rule = 2)
  val <- list()
  val$call <- match.call()
  val$x <- probs
  val$y <- Qn(probs)
  val$fn <- Qn
  val$data <- x
  class(val) <- "midquantile"
  return(val)
}

