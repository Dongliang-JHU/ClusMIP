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
# R Script Purpose: to compute the tau's for clean and contaminated influential observations for real data   
#
# Created on: November 03, 2023  
#
# Modified on: November 03, 2023
#              December 02, 2023 
#
# R Packages: glmnet, scalreg, ncvreg 
#
# Functions: tau.compute, tau.one.out
# 
################################################################################ 
# Function: tau.compute
#
# Args:
#
#   X                :  potentially contaminated design matrx 
#   y                :  potentially contaminated response vector  
#   nfolds           :  number of folds for cross-validation  
#   dat.type         :  "gaussian" or "binomial" 
#   infl.pos.single  :  the index of a single potentially influential point 
#   clean.pos        :  indices of estimated clean dataset from clustering  
#   selector.switch  :  a binary ( 0 or 1) vector indicating LASSO, SLASSO, SCAD and MCP 
#                       selector.switch=1: LASSO
#                       selector.switch=2: SLASSO
#                       selector.switch=3: SCAD
#                       selector.switch=4: MCP
#
# Returns: 
#
#   tau.lasso.clean   
#   tau.slasso.clean 
#   tau.scad.clean 
#   tau.mcp.clean
#   tau.lasso.infl 
#   tau.slasso.infl 
#   tau.scad.infl 
#   tau.mcp.infl
#   tau.lasso.infl.mat 
#   tau.slasso.infl.mat 
#   tau.scad.infl.mat 
#   tau.mcp.infl.mat 
#
# Dependent R functions:   
#

tau.compute.real <- function(X, y, nfolds, dat.type, infl.pos.single, clean.pos, selector.switch){
  
  n <- nrow(X)
  p <- ncol(X)
  
  X.clean <- as.matrix(X[clean.pos,])
  y.clean <- y[clean.pos]
  
  X.infl <- rbind(X[infl.pos.single,], X[clean.pos,])
  y.infl <- c(y[infl.pos.single], y[clean.pos])

  ##############################################################################
  # LASSO
  if(selector.switch[1]==1){
    
    # clean 
    lasso.clean <- cv.glmnet(as.matrix(X.clean), y.clean, alpha=1, family=dat.type, nfolds=nfolds)   
    lasso.clean.coef <- as.numeric(coef(lasso.clean, s="lambda.min")[-1]) 
    lasso.clean.zero <- which(lasso.clean.coef==0) 
    lasso.clean.ind <- rep(0,p) 
    lasso.clean.ind[lasso.clean.zero] <- 1  
    
    tau.clean <- tau.one.out(X=X.clean,y=y.clean,selector.switch=c(1,0,0,0),dat.type,nfolds,indicator.complete=lasso.clean.ind)  
    
    lasso.infl <- cv.glmnet(as.matrix(X.infl), y.infl, alpha=1, family=dat.type, nfolds=nfolds)   
    lasso.infl.coef <- as.numeric(coef(lasso.infl, s="lambda.min")[-1]) 
    lasso.infl.zero <- which(lasso.infl.coef==0) 
    lasso.infl.ind <- rep(0,p) 
    lasso.infl.ind[lasso.infl.zero] <- 1  
    
    tau.lasso.inf <- tau.one.out(X=X.infl,y=y.infl,selector.switch=c(1,0,0,0),dat.type,nfolds,indicator.complete=lasso.infl.ind)
    
    tau.infl.vec <- tau.lasso.inf[1]
    tau.infl.mat <- tau.lasso.inf

  }
  
  ##############################################################################
  # SLASSO
  if(selector.switch[2]==1){
    
    # clean 
    slasso.clean <- scalreg(as.matrix(X.clean), y.clean, lam0=NULL, LSE=FALSE)
    slasso.clean.coef <- as.numeric(slasso.clean$coefficients)  
    slasso.clean.zero <- which(slasso.clean.coef==0) 
    slasso.clean.ind <- rep(0,p) 
    slasso.clean.ind[slasso.clean.zero] <- 1    
    
    tau.clean <- tau.one.out(X=X.clean,y=y.clean,selector.switch=c(0,1,0,0),dat.type,nfolds,indicator.complete=slasso.clean.ind)
    
    slasso.infl <- scalreg(as.matrix(X.infl), y.infl, lam0 = NULL, LSE = FALSE) 
    slasso.infl.coef <- as.numeric(slasso.infl$coefficients)  
    slasso.infl.zero <- which(slasso.infl.coef==0) 
    slasso.infl.ind <- rep(0,p) 
    slasso.infl.ind[slasso.infl.zero] <- 1  
    
    tau.slasso.infl <- tau.one.out(X=X.infl,y=y.infl,selector.switch=c(0,1,0,0),dat.type,nfolds,indicator.complete=slasso.infl.ind)
    
    tau.infl.vec <- tau.slasso.infl[1]
    tau.infl.mat <- tau.slasso.infl
    
  }
  
  ##############################################################################
  # SCAD
  
  if(selector.switch[3]==1){
    
    # clean 
    scad.clean.obj <- cv.ncvreg(as.matrix(X.clean), y.clean, family=dat.type, penalty="SCAD") 
    scad.clean <- scad.clean.obj$fit 
    scad.clean.coef <- as.numeric(scad.clean$beta[,scad.clean.obj$min][-1])
    scad.clean.zero <- which(scad.clean.coef==0) 
    scad.clean.ind <- rep(0,p) 
    scad.clean.ind[scad.clean.zero] <- 1      
    
    tau.clean <- tau.one.out(X=X.clean,y=y.clean,selector.switch=c(0,0,1,0),dat.type,nfolds,indicator.complete=scad.clean.ind)
    
    scad.infl.obj <- cv.ncvreg(as.matrix(X.infl), y.infl, family=dat.type, penalty="SCAD")
    scad.infl <- scad.infl.obj$fit 
    scad.infl.coef <- as.numeric(scad.infl$beta[,scad.infl.obj$min][-1])
    scad.infl.zero <- which(scad.infl.coef==0) 
    scad.infl.ind <- rep(0,p) 
    scad.infl.ind[scad.infl.zero] <- 1  
    
    tau.scad.infl <- tau.one.out(X=X.infl,y=y.infl,selector.switch=c(0,0,1,0),dat.type,nfolds,indicator.complete=scad.infl.ind)
    
    tau.infl.vec <- tau.scad.infl[1]
    tau.infl.mat <- tau.scad.infl
  }
  
  ##############################################################################
  # MCP
  
  if(selector.switch[4]==1){
    
    # clean 
    mcp.clean.obj <- cv.ncvreg(as.matrix(X.clean), y.clean, family=dat.type, penalty="MCP") 
    mcp.clean <- mcp.clean.obj$fit 
    mcp.clean.coef <- as.numeric(mcp.clean$beta[,mcp.clean.obj$min][-1])
    mcp.clean.zero <- which(mcp.clean.coef==0) 
    mcp.clean.ind <- rep(0,p) 
    mcp.clean.ind[mcp.clean.zero] <- 1   
    
    tau.clean <- tau.one.out(X=X.clean,y=y.clean,selector.switch=c(0,0,0,1),dat.type,nfolds,indicator.complete=mcp.clean.ind)
    
    mcp.infl.obj <- cv.ncvreg(as.matrix(X.infl), y.infl, family=dat.type, penalty="MCP")
    mcp.infl <- mcp.infl.obj$fit 
    mcp.infl.coef <- as.numeric(mcp.infl$beta[,mcp.infl.obj$min][-1])
    mcp.infl.zero <- which(mcp.infl.coef==0) 
    mcp.infl.ind <- rep(0,p) 
    mcp.infl.ind[mcp.infl.zero] <- 1  
    
    tau.mcp.infl <- tau.one.out(X=X.infl,y=y.infl,selector.switch=c(0,0,0,1),dat.type,nfolds,indicator.complete=mcp.infl.ind)
    
    tau.infl.vec <- tau.mcp.infl[1]
    tau.infl.mat <- tau.mcp.infl
  }
  
  ##############################################################################
  
  returnList <- list("tau.clean" = tau.clean, 
                     "tau.infl" = tau.infl.vec, 
                     "tau.infl.mat" = tau.infl.mat)
  
  return(returnList)  
  
}

################################################################################ 
# Function: tau.one.out
#
# Args:
# X                  : design matrix
# y                  : response vector 
# selector.switch    : a binary (0 or 1) vector:LASSO, SLASSO, SCAD and MCP
# dat.type           : gaussian or binomial       :  
# nfolds             : number of folds for cross-validation  
# indicator.complete : a binary (0 or 1) vector indicating zeros of beta on the complete dataset  
#
# Returns: 
# tau : tau for assessing influence 

tau.one.out <- function(X,y,selector.switch,dat.type,nfolds,indicator.complete){
  
  p <- ncol(X)
  
  tau <- sapply(1:nrow(X), function(j){
    
    X.sub <- X[-j,]
    y.sub <- y[-j] 
    
    # selector.switch=1: LASSO
    # selector.switch=2: SLASSO
    # selector.switch=3: SCAD
    # selector.switch=4: MCP
    
    if(selector.switch[1]==1){
      fit <- cv.glmnet(as.matrix(X.sub), y.sub, alpha=1, family=dat.type, nfolds=nfolds)  
      coef <- as.numeric(coef(fit, s="lambda.min")[-1])  
    }else if(selector.switch[2]==1){
      fit <- scalreg(as.matrix(X.sub), y.sub, lam0 = NULL, LSE = FALSE) 
      coef <- as.numeric(fit$coefficients)
    }else if(selector.switch[3]==1){
      scad <- cv.ncvreg(as.matrix(X.sub), y.sub, family=dat.type, penalty="SCAD") 
      fit <- scad$fit
      coef <- as.numeric(fit$beta[,scad$min][-1])
    }else if(selector.switch[4]==1){
      mcp <- cv.ncvreg(as.matrix(X.sub), y.sub, family=dat.type, penalty="MCP")
      fit <- mcp$fit 
      coef <- as.numeric(fit$beta[,mcp$min][-1])
    }else{
      return("Choose a model selector from LASSO, SLASSO, SCAD and MCP")
      break
    }
    
    zero.pos <- which(coef==0) 
    indicator.sub <- rep(0,p) 
    indicator.sub[zero.pos] <- 1  
    sum((indicator.complete-indicator.sub)^2) 	  
  })
  
  return(tau)
  
}
