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
# R Script Purpose: Parametric Approximation Approach
#
# Created on: October 23, 2023
#
# Modified on: October 23, 2023
#              October 24, 2023
#              October 25, 2023
#              October 26, 2023
#              October 27, 2023
#              October 28, 2023
#              October 29, 2023
#
# Dependent R Scripts: midquantile_param.r, supplementary.r
#
# R Packages:
#
# Total Functions: param.approx
#
################################################################################
# Function: param.approx
#
# Args:
#
# tau.vec     : tau vector
# p           : number of predictors
# mb.maxComp  : upper bound of number of components for binomial mixture
# mp.maxComp  : upper bound of number of components for Poisson mixture
# n.boot      : number of bootstrap samples for Poisson mixutre order selection (likelihood-based approach)
# alpha       : FDR
#
# Returns:
#
# 6 thresholds based on the 6 parametric approximating distributions
#
# Dependent R functions: supplementary.r, midquantile_param.r

#' Title
#'
#' @param tau.vec tau vector
#' @param p number of predictors
#' @param mb.maxComp upper bound of number of components for binomial mixture
#' @param mp.maxComp upper bound of number of components for Poisson mixture
#' @param n.boot number of bootstrap samples for Poisson mixutre order selection (likelihood-based approach)
#' @param alpha FDR
#'
#' @return 6 thresholds based on the 6 parametric approximating distributions
#' @export
#'

param.approx <- function(tau.vec,p,mb.maxComp,mp.maxComp,n.boot,alpha){



  ##############################################################################
  # Conway-Maxwell-Binomial (CMB)
  # Checked

  skip_to_next_outer <- FALSE

  tryCatch({
    tau.tmp <- as.matrix(cbind(tau.vec, p-tau.vec))
    cmp.dat <- data.frame(y1=tau.tmp[,1], y2=tau.tmp[,2], x=1, w=1)
    k <- 2
    cmb.fit <- cmm_reg(formula_x=tau.tmp~x-1, formula_w=~1, data=cmp.dat)
    cmb.p <- 1/(as.numeric(coef(cmb.fit)$beta)+1) # log((1-p)/p)=beta
    cmb.nu <- as.numeric(coef(cmb.fit)$gamma) # identity link
    cmb.threshold <- midquantile.param(x=tau.vec, p, param.string="cmb", param.val=c(cmb.p,cmb.nu), probs=1-alpha)$y},
  error = function(e){skip_to_next_outer <<- TRUE})

  if(skip_to_next_outer==TRUE){cmb.threshold <- as.numeric(quantile(tau.vec,1-alpha))}

  ##############################################################################
  # Conway-Maxwell-Poisson (CMP)
  # Checked

  skip_to_next_outer <- FALSE

  tryCatch({
    cmp.fit <- COMPoissonReg::glm.cmp(formula.lambda=tau.vec~1, formula.nu=~1)
    cmp.lambda <- exp(as.numeric(cmp.fit$beta)) # log link
    cmp.nu <- exp(as.numeric(cmp.fit$gamma)) # log link
    cmp.threshold <- midquantile.param(x=tau.vec, p, param.string="cmp", param.val=c(cmp.lambda,cmp.nu), probs=1-alpha)$y},
  error = function(e){skip_to_next_outer <<- TRUE})

  if(skip_to_next_outer==TRUE){cmp.threshold <- as.numeric(quantile(tau.vec,1-alpha))}

  ##############################################################################
  # Beta-Binomial (BB)
  # Checked

  skip_to_next_outer <- FALSE

  tryCatch({
    bb.fit <- bb.mle(x=tau.vec, n=p, alpha1=1, alpha2=1) #(1,1) initial values
    bb.alpha1 <- as.numeric(bb.fit[1,2])
    bb.alpha2 <- as.numeric(bb.fit[1,3])
    bb.threshold <- midquantile.param(x=tau.vec, p, param.string="bb", param.val=c(bb.alpha1,bb.alpha2), probs=1-alpha)$y},
  error = function(e){skip_to_next_outer <<- TRUE})

  if(skip_to_next_outer==TRUE){bb.threshold <- as.numeric(quantile(tau.vec,1-alpha))}

  ##############################################################################
  # Generalized Poisson (GP)
  # Checked

  skip_to_next_outer <- FALSE

  tryCatch({
    gp.dat <- data.frame(x=1, y=tau.vec)
    suppressWarnings({gp.fit <- vglm(y~1, genpoisson0, data=gp.dat, trace=FALSE)}) # summary: log link for theta, logit link for lambda
    gp.coef <- as.numeric(coef(gp.fit, matrix=TRUE))
    gp.theta <- exp(gp.coef[1]) # log link
    gp.lambda <- exp(gp.coef[2]) / ( 1 + exp(gp.coef[2]) ) # logit link
    gp.threshold <- midquantile.param(x=tau.vec, p, param.string="gp", param.val=c(gp.theta,gp.lambda), probs=1-alpha)$y},
  error = function(e){skip_to_next_outer <<- TRUE})

  if(skip_to_next_outer==TRUE){gp.threshold <- as.numeric(quantile(tau.vec,1-alpha))}

  ##############################################################################
  # Mixture of Binomial (MB)
  # Checked

  skip_to_next_outer <- FALSE

  tryCatch({
    mb.dat <- matrix(0,nrow=length(tau.vec),ncol=2)
    mb.dat[,1] <- tau.vec
    mb.dat[,2] <- p-tau.vec
    mb.mix <- BinomialMixtures(maxNumComponents=mb.maxComp, phi="default")
    mb.sbic <- sBIC(mb.dat, mb.mix)
    mb.n.comp <- (1:mb.maxComp)[which.max(mb.sbic$sBIC)]},
  error = function(e){skip_to_next_outer <<- TRUE})

  if(skip_to_next_outer==TRUE){mb.n.comp <- mixbinom.OrderSelection.BIC(tau.vec,p,mb.maxComp)}

  skip_to_next_outer <- FALSE

  tryCatch({
    mb.fit <- mixBinom(k=tau.vec, n=rep(p,length(tau.vec)), n_components=mb.n.comp)
    mb.success <- mb.fit$p # success rates of components
    mb.psi <- mb.fit$psi # fractions of components
    param.mat <- rbind(mb.psi, mb.success)
    mb.threshold <- midquantile.param(x=tau.vec, p, param.string="mb", param.val=param.mat, probs=1-alpha)$y},
  error = function(e){skip_to_next_outer <<- TRUE})

  if(skip_to_next_outer==TRUE){mb.threshold <- as.numeric(quantile(tau.vec,1-alpha))}

  ##############################################################################
  # Mixture of Poisson (MP)
  # Checked

  skip_to_next_outer <- FALSE

  tryCatch({
    poisson.order <- poissonOrder(tau.vec,K=mp.maxComp,penalty="ADAPTIVE-LASSO",lambdas=c(0, 0.001, 0.01, 0.1, 0.5, 1, 1.5, 2), verbose=F)
    gsf.fit <- bicTuning.modified(tau.vec, poisson.order)

    if(typeof(gsf.fit)=="double"){
      mp.n.comp <- pois.OrderSelection(tau.vec,mp.maxComp,n.boot,alpha)
    }else{
      mp.n.comp <- gsf.fit$result$order
    }

    mp.dat <- data.frame(tau.vec)
    mp.fit <- flexmix(tau.vec~1, data=mp.dat, k=mp.n.comp, model=FLXMRglm(family="poisson"))
    mp.lambda <- as.numeric(exp(parameters(mp.fit)))
    mp.psi <- prior(mp.fit)
    param.mat <- rbind(mp.psi, mp.lambda)
    mp.threshold <- midquantile.param(x=tau.vec, p, param.string="mp", param.val=param.mat, probs=1-alpha)$y},

  error = function(e){skip_to_next_outer <<- TRUE})

  if(skip_to_next_outer==TRUE){mp.threshold <- as.numeric(quantile(tau.vec,1-alpha))}

  ##############################################################################

  return(c(cmb.threshold, cmp.threshold, bb.threshold, gp.threshold, mb.threshold, mp.threshold))

}
