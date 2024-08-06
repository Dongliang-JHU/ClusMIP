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
# R Script Purpose: Supplementary Functions
#
# Created on: October 26, 2023 
#
# Modified on:  October 26, 2023 
#               October 27, 2023 
#               October 28, 2023 
#
# R Packages :  
#
# Functions: decision.fun, 
#            clt.tau, 
#            pois.OrderSelection, 
#            pois.OrderSelection.getSample,
#            mixbinom.OrderSelection.BIC
#                      
################################################################################ 
# Function: decision.fun
#
# Args:
#
# Returns: 
#

decision.fun <- function(list.pos, mat){
  
  l <- length(list.pos)
  
  for(i in 1:l){
    mat[i,list.pos[[i]]] <- 1
  }
  
  return(mat)
}

################################################################################ 
# Function: clt.tau 
#
# Args:
# 
# vec   : vector of tau 
# alpha : FDR
#
# Returns: 
#
# TRUE or FALSE indicating whether the first element is influential or not based on the 
# Central Limit Theorem 
#

clt.tau <- function(vec,alpha){
  tmp <- abs( (vec[1]-mean(vec))/sd(vec) )
  normal.quantile <- qnorm(1-alpha/2,mean=0,sd=1)
  return(tmp >= normal.quantile)
}

################################################################################ 
# Function: pois.OrderSelection
#
# Args:
#
# tau.vec     : tau vector    
# mp.maxComp  : upper bound on the number of Poisson mixture components 
# n.boot      : number of bootstrap samples 
# alpha       : FDR 
#
# Returns: 
#
# estimate on the order of Poisson mixtures based on the reference below 
#
# References: 
# 1. Karlis and Xekalaki (1999), On Testing for the Number of Components in a Mixed Poisson Model  
#                                Annals of the Institute of Statistical Mathematics 
# 

pois.OrderSelection <- function(tau.vec,mp.maxComp,n.boot,alpha){

  n.order <- 1
  
  ##############################################################################
  # iterative sampling scheme 
  
  mp.fit.pre <- flexmix(tau.vec~1, data=data.frame(tau.vec), k=1, model=FLXMRglm(family="poisson"))
  
  mp.fit.next <- flexmix(tau.vec~1, data=data.frame(tau.vec), k=2, model=FLXMRglm(family="poisson"))
  
  ratio <- (-2) * ( as.numeric(logLik(mp.fit.pre)) - as.numeric(logLik(mp.fit.next)) )
  
  param.pre <- matrix(c(as.numeric(exp(parameters(mp.fit.pre))), prior(mp.fit.pre)), ncol=length(prior(mp.fit.pre)), byrow=T)
  
  cut.off <- pois.OrderSelection.getSample(n=length(tau.vec),
                                           param.null=param.pre,
                                           n.boot,
                                           loglike.next=as.numeric(logLik(mp.fit.next)),
                                           alpha)
  
  
  while(ratio>=cut.off & n.order<=mp.maxComp){
    
    mp.fit.pre <- flexmix(tau.vec~1, data=data.frame(tau.vec), k=n.order+1, model=FLXMRglm(family="poisson"))
    
    mp.fit.next <- flexmix(tau.vec~1, data=data.frame(tau.vec), k=n.order+2, model=FLXMRglm(family="poisson"))
    
    ratio <- (-2) * ( as.numeric(logLik(mp.fit.pre)) - as.numeric(logLik(mp.fit.next)) )
    
    param.pre <- matrix(c(as.numeric(exp(parameters(mp.fit.pre))), prior(mp.fit.pre)), ncol=length(prior(mp.fit.pre)), byrow=T)
    
    cut.off <- pois.OrderSelection.getSample(n=length(tau.vec),
                                             param.null=param.pre,
                                             n.boot,
                                             loglike.next=as.numeric(logLik(mp.fit.next)),
                                             alpha)
    
    n.order <- n.order + 1 
  }
  
  ##############################################################################
  
  return(n.order)
}

################################################################################ 
# Function: pois.OrderSelection.getSample
#
# Args:
#
# n            : sample size  
# param.null   : MLE estimates under the null distribution 
# n.boot       : number of bootstrap samples 
# loglike.next : log-likelihood for the next iteration 
# alpha        : FDR
#
# Returns: 
#
# quantile of the LRT ratio statistics 
# 
pois.OrderSelection.getSample <- function(n,param.null,n.boot,loglike.next,alpha){
  
  ratio.boot <- c()
  
  for(i in 1:n.boot){
    mixpois.sample.boot <- rmixpois(n, lambda=param.null[1,], alpha=param.null[2,]) # under the null hypothesis 
    loglike.boot.null <- logLikePoisMix(y=as.matrix(mixpois.sample.boot), mean=param.null[1,], pi=param.null[2,])$ll
    ratio.boot <- c(ratio.boot, (-2)*( loglike.boot.null - loglike.next ))
  }
  
  return(as.numeric(quantile(ratio.boot,1-alpha)))
  
}

################################################################################ 
# Function: mixbinom.OrderSelection.BIC
#
# Args:
#
# tau.vec.    : tau vector
# mp.maxComp  : upper bound of the number of Poisson mixture components  
# n.boot      : number of bootstrap samples 
#
# Returns: 
#
# estimate on the order of binomial mixtures based on BIC 

mixbinom.OrderSelection.BIC <- function(tau.vec,p,mb.maxComp){
  
  bic.vec <- c()
  
  comp.n.vec <- 1:mb.maxComp
  
  for(i in 1:mb.maxComp){
    mb.fit <- mixBinom(k=tau.vec, n=rep(p,length(tau.vec)), n_components=i)
    k <- (i-1) + i 
    bic.mb <- (-2)*mb.fit$logLik + k*log(length(tau.vec))
    bic.vec <- c(bic.vec,bic.mb)
  }
  
  return(which.min(bic.vec))

}
