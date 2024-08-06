# Author: Dongliang Zhang (^)
#
# Supervisors: Professor Masoud Asgharian (*)
#              Professor Martin Lindquist (^)
#
# Current Affiliation: (1) Department of Mathematics and Statistics, McGill University, Montreal, Quebec, Canada (*)
#                      (2) Department of Biostatistics, Johns Hopkins University, Baltimore, Maryland, USA (^)
#
# Title of Project: Detection of Influential Cases on Variable Selection
#
# R Script Purpose: Difference in Model Selected by LASSO (DF measure)
#
# References: B. Rajaratnam, S. Roberts, D. Sparks and H. Yu, Inflence Diagnostics for High-Dimensional Lasso Regression,
#             Journal of Computational and Graphical Statistics, 2019.
#
# Created on      : April 05, 2022
#
# Modified on     : April 05, 2022
#                   November 01, 2023
#
# Dependent R Scripts:
#
# R Packages         : glmnet;
#
#############################################################################################################################
# Function: df.measure
#
# Args:
#   X                :  the design matrix
#   Y                :  response vector
#   dat.type.        :  type of data
#
# Returns:

#
# R Packages         :  glmnet;
#

#' DF(LASSO)
#'
#' @param X design matrix
#' @param Y response vector
#' @param dat.type types of data "gaussian" or "binomial"
#'
#' @return a decision vector indicating the positions of clean and influential observations
#' @export
#'

df.measure <- function(X,Y,dat.type){

  n <- nrow(X);
  p <- ncol(X);

  dec.vec <- rep(0,n)

  tau.vec <- vector(mode="numeric", length=n)

  #LASSO
  glm.obj.full <- glmnet::cv.glmnet(X, Y, alpha=1, family=dat.type, nfolds=10, intercept=FALSE);
  coeff.full <- as.numeric(coef(glm.obj.full, s="lambda.min")[-1]);
  boo.coeff <- which(coeff.full==0);
  boo.org <- rep(0,p);
  boo.org[boo.coeff] <- 1;

  #"for" loop for leave-one-out algorithm
  for(j in 1:n){
    Y.red <- Y[-j]
    X.red <- X[-j,]
    glm.obj.red <- glmnet::cv.glmnet(X.red, Y.red, alpha=1, family=dat.type, nfolds=10, intercept=FALSE)
    coeff.red <- as.numeric(coef(glm.obj.red, s="lambda.min")[-1])
    boo.coeff.red <- which(coeff.red==0)
    boo.red <- rep(0,p)
    boo.red[boo.coeff.red] <- 1
    tau.vec[j] <- sum((boo.org-boo.red)^2)
  }

  vec.tmp <- abs((tau.vec-mean(tau.vec)) / sd(tau.vec))

  pos.infl <- which(vec.tmp >= 2)

  if(length(pos.infl) > 0){
    dec.vec[pos.infl] <- 1
  }

  return(dec.vec)

}
