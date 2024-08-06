# Author: Dongliang Zhang (^)
#
# Supervisors: Professor Masoud Asgharian (*)
#              Professor Martin Lindquist (^)
#
# Current Affiliation: (*) Department of Mathematics and Statistics, McGill University, Montreal, Quebec, Canada
#                      (^) Department of Biostatistics, Johns Hopkins University, Baltimore, Maryland, USA
#
# Title of Project: Detection of Multiple Influential Observations on Model Selection
#
# R Script Purpose: Main function for the ClusMIP package
#
# Created on: August 05, 2024
#
# Modified on: August 05, 2024
#
# Dependent R Scripts: cluster.r, tau.r, df_measure.r, param_approx.r, boot_approx.r, midquantile_param.r, supplementary.r
#
# R Packages :
#
################################################################################
# Args:
#
#   X                :  design matrix
#   y                :  response vector
#   cluster.type     :  five type of clustering from "kmeans", "kmeans++", "tsne", "spectral" and "rkmeans"
#   dat.type         :  "gaussian" or "binomial"
#   nfolds           :  number of folds for cross-validation
#   n.boot           :  number of bootstrap samples
#   alpha            :  FDR
#   mb.maxComp       :  upper bound of number of components for binomial mixture
#   mp.maxComp       :  upper bound of number of components for Poisson mixture
#   selector.switch  :  a binary vector of length 4 indicating which selector is used
#
# Returns:
#   infl.pos         :  position of potentially influential observations
#   clean.pos        :  positions of potentially clean observations
#   mip              :  binary decision vector of length n for MIP
#   df.lasso         :  binary decision vector of length n for the DF(LASSO)
#   clusmip          :  a list of length of "infl.pos"
#                       each element of the list contains a matrix of 0 and 1
#                       the rows correspond to which selectors are used
#                       the columns correspond to the 10 approximation methods

#' Main function for assessing influential observations by different detection procedures.
#'
#'
#' @param X design matrix
#' @param y response vector
#' @param cluster.type five type of clustering from "kmeans", "kmeans++", "tsne", "spectral" and "rkmeans"
#' @param dat.type "gaussian" or "binomial"
#' @param nfolds number of folds for cross-validation
#' @param n.boot number of bootstrap samples
#' @param alpha FDR
#' @param mb.maxComp upper bound of number of components for binomial mixture
#' @param mp.maxComp upper bound of number of components for Poisson mixture
#' @param selector.switch a binary vector of length 4 indicating which selector is used
#'
#' @return A list indicating positions of the potentially influential and clean observations, decision vector from MIP and DFLASSO, and decision from the ClusMIP method for each model selector and approximation method.
#'
#' @export
#'

clusmip_main <- function(X, y, cluster.type, dat.type, nfolds, n.boot, alpha, mb.maxComp, mp.maxComp, selector.switch){

  n <- nrow(X)
  p <- ncol(X)

  # clustering
  Xy.dat <- as.data.frame(cbind(X, y))
  names(Xy.dat)<-c(paste("X", 1:p, sep=""),"y")
  cluster.fit <- cluster.est(df=Xy.dat, nclust=2, cluster.type)
  infl.pos <- cluster.fit$infl.pos
  clean.pos <- cluster.fit$clean.pos

  ##############################################################################
  # MIP
  if(dat.type=="gaussian"){
    detection.mip <- rep(0,n)
    infl.pos.MIP <- MIP(X=as.matrix(X), Y=t(as.matrix(y)), n, p, q=1, n_subset=100, subset_vol=n/2, ep=0.1, alpha)$inf_setfinal
    infl.pos.MIP <- sort(infl.pos.MIP, decreasing=FALSE)
    detection.mip[infl.pos.MIP] <- 1
  }else{
    detection.mip <- NA
  }

  ##############################################################################
  # DF(LASSO)
  detection.dflasso <- df.measure(X=X, Y=y, dat.type)

  ##############################################################################
  # ClusMIP (revised)
  # 4 model selectors (LASSO, SLASSO, SCAD, MCP)
  # 10 approximation methods (CLT, 6 parametric and 3 non-parametric)

  dec.clusmip.list <- vector(mode="list", length=length(infl.pos))

  names(dec.clusmip.list) <- as.character(infl.pos)

  for(i in 1:length(infl.pos)){

    infl.pos.single <- infl.pos[i]

    dec.mat <- matrix(0,nrow=4,ncol=10)
    rownames(dec.mat) <- c("LASSO", "SLASSO", "SCAD", "MCP")
    colnames(dec.mat) <- c("CLT", "CMB", "CMP", "BB", "GP", "MB", "MP", "BootI", "BootII", "BootIII")

    # LASSO
    if(selector.switch[1]==1){

      # tau computation
      tau.lasso <- tau.compute.real(X, y, nfolds, dat.type, infl.pos=infl.pos.single, clean.pos, selector.switch=c(1,0,0,0))

       # CLT
      tau.clusmip.lasso <- tau.lasso$tau.infl.mat
      dec.mat[1,1] <- as.integer(clt.tau(tau.clusmip.lasso,alpha))

      # Parametric approximation
      tau.lasso.clean <- tau.lasso$tau.clean
      tau.lasso.infl <- tau.lasso$tau.infl
      lasso.param.threshold <- floor(param.approx(tau.vec=tau.lasso.clean,p,mb.maxComp,mp.maxComp,n.boot,alpha))
      dec.mat[1,2:7] <- as.integer(tau.lasso.infl > lasso.param.threshold)

      # Non-parametric (Bootstrap) Approximation
      lasso.boot.threshold <- boot.approx(tau.lasso.clean, n.boot, alpha)
      dec.mat[1,8:10] <- as.integer(tau.lasso.infl > lasso.boot.threshold)
    }

    # SLASSO
    if(selector.switch[2]==1){
      if(dat.type=="gaussian"){

        # tau computation
        tau.slasso <- tau.compute.real(X, y, nfolds, dat.type, infl.pos=infl.pos.single, clean.pos, selector.switch=c(0,1,0,0))

        # CLT
        tau.clusmip.slasso <- tau.slasso$tau.infl.mat
        dec.mat[2,1] <- as.integer(clt.tau(tau.clusmip.slasso,alpha))

        # Parametric approximation
        tau.slasso.clean <- tau.slasso$tau.clean
        tau.slasso.infl <- tau.slasso$tau.infl
        slasso.param.threshold <- floor(param.approx(tau.vec=tau.slasso.clean,p,mb.maxComp,mp.maxComp,n.boot,alpha))
        dec.mat[2,2:7] <- as.integer(tau.slasso.infl > slasso.param.threshold)

        # Non-parametric (Bootstrap) Approximation
        slasso.boot.threshold <- boot.approx(tau.slasso.clean, n.boot, alpha)
        dec.mat[2,8:10] <- as.integer(tau.slasso.infl > slasso.boot.threshold)
      }else{
        dec.mat[2,] <- rep(NA,10)
      }
    }

    # SCAD
    if(selector.switch[3]==1){
      # tau computation
      tau.scad <- tau.compute.real(X, y, nfolds, dat.type, infl.pos=infl.pos.single, clean.pos, selector.switch=c(0,0,1,0))

      # CLT
      tau.clusmip.scad <- tau.scad$tau.infl.mat
      dec.mat[3,1] <- as.integer(clt.tau(tau.clusmip.scad,alpha))

      # Parametric approximation
      tau.scad.clean <- tau.scad$tau.clean
      tau.scad.infl <- tau.scad$tau.infl
      scad.param.threshold <- floor(param.approx(tau.vec=tau.scad.clean,p,mb.maxComp,mp.maxComp,n.boot,alpha))
      dec.mat[3,2:7] <- as.integer(tau.scad.infl > scad.param.threshold)

      # Non-parametric (Bootstrap) Approximation
      scad.boot.threshold <- boot.approx(tau.scad.clean, n.boot, alpha)
      dec.mat[3,8:10] <- as.integer(tau.scad.infl > scad.boot.threshold)
    }

    # MCP
    if(selector.switch[4]==1){
      # tau computation
      tau.mcp <- tau.compute.real(X, y, nfolds, dat.type, infl.pos=infl.pos.single, clean.pos, selector.switch=c(0,0,0,1))

      # CLT
      tau.clusmip.mcp <- tau.mcp$tau.infl.mat
      dec.mat[4,1] <- as.integer(clt.tau(tau.clusmip.mcp,alpha))

      # Parametric approximation
      tau.mcp.clean <- tau.mcp$tau.clean
      tau.mcp.infl <- tau.mcp$tau.infl
      mcp.param.threshold <- floor(param.approx(tau.vec=tau.mcp.clean,p,mb.maxComp,mp.maxComp,n.boot,alpha))
      dec.mat[4,2:7] <- as.integer(tau.mcp.infl > mcp.param.threshold)

      # Non-parametric (Bootstrap) Approximation
      mcp.boot.threshold <- boot.approx(tau.mcp.clean, n.boot, alpha)
      dec.mat[4,8:10] <- as.integer(tau.mcp.infl > mcp.boot.threshold)

    }

    dec.clusmip.list[[i]] <- dec.mat[which(selector.switch==1),]

  }

  ##############################################################################
  # return

  returnList <- list("infl.pos" = infl.pos,
                     "clean.pos" = clean.pos,
                     "mip" = detection.mip,
                     "dflasso" = detection.dflasso,
                     "clusmip" = dec.clusmip.list)

  return(returnList)

}
