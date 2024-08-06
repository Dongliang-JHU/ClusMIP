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
# R Script Purpose: Non-parametric (Bootstrap) Approach
#
# Created on: October 25, 2023
#
# Modified on: October 25, 2023
#              October 26, 2023
#
# R Packages: Qtools
#
# Total Functions: boot.approx
#
################################################################################
# Function: boot.approx
#
# Args:
#
# tau.vec: tau vector
# n.boot: number of bootstrap samples
# alpha: FDR
#
# Returns:
#
# 3 thresholds based on the three bootstrap sampling schemes
#
# References:
# 1. Muliere and Secchi (1992), Exhangeability, Predictive Sufficiency and Bayesian Bootstrap, Journal of Italian Statistical Society.
#
# 2. Jentsch and Leucht (2016), Bootstrapping Sample Quantiles of Discrete Data, Annals of the Institute of Statistical Mathematics.

#' Threshold computation from three bootstrap schemes
#'
#' @param tau.vec vector of tau
#' @param n.boot number bootstrap copies
#' @param alpha FDR
#'
#' @return a vector of length 3 indicating 3 thresholds
#' @export
#'

boot.approx <- function(tau.vec, n.boot, alpha){

  n.tau <- length(tau.vec)

  # m-out-of-n low-intensity bootstrap
  m <- ceiling(n.tau^(3/4))

  mean.boot <- c()
  quantile.boot <- c()
  midquantile.boot <- c()

  ##############################################################################

  for(i in 1:n.boot){

    tau.boot <- sample(tau.vec, size=n.tau, replace=TRUE)

    # Boot I: Muliere and Secchi (1992)
    mean.boot <- c(mean.boot, mean(tau.boot))

    # Boot II: Jentsch and Leucht (2016)
    tau.m.outof.n <- sample(tau.vec, size=m, replace=TRUE)
    quantile.boot <- c(quantile.boot, as.numeric(quantile(tau.m.outof.n, 1-alpha)))

    # Boot III: Jentsch and Leucht (2016)
    skip_to_next_outer <- FALSE

    tryCatch({
      midquant <- Qtools::midquantile(tau.m.outof.n, probs=1-alpha)$y},
    error = function(e){skip_to_next_outer <<- TRUE})

    if(skip_to_next_outer==TRUE){midquant <- as.numeric(quantile(tau.m.outof.n,1-alpha))}

    midquantile.boot <- c(midquantile.boot, midquant)
  }

  ##############################################################################

  bootI.threshold <- ceiling(as.numeric(quantile(mean.boot,1-alpha)))

  bootII.threshold <- ceiling(mean(quantile.boot))

  bootIII.threshold <- ceiling(mean(midquantile.boot))

  return(c(bootI.threshold, bootII.threshold, bootIII.threshold))
}
