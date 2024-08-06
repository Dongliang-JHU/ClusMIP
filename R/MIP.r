library(Matrix)
library(MASS)

###############################
##main function
##MIP: function to detect multiple influential points
##Usage: MIP(X,Y,n,p,q,n_subset,subset_vol,ep=0.1,alpha)
##input:
# X: the data of predictors with dimension n by p
# Y: the data of response with dimension  n by q
# n: the sample size
# p: the dimension of predictor
# q: the dimension of response
# n_subset: the number of subsets  chosen at random to compute the Min and Max statistics
# subset_vol: the  samples size in each subset
# ep: the proportion of maximum number of  rejected null hypothesis in the Min-step. The defaulted value is set at  0.1.
# alpha: significant level used in FDR procedure
##output
#inf_setfinal: the indices of the influential points detected by MIP algorithm
####################################################




####explanation  for the  functions used
##########################################################
##fun_pv: function to comupute the max-statistics and the min-statistics
##Usage: fun_pv(X,Y,n,p,q,n_subset,subset_vol,clean_setv)
##input:
#clean_setv: the estimated clean set obtained by Min/Max step
##output:
#Tmax: the values of the max-statistics
#Tmin: the values of the min-statistics
##########################################################



##########################################################
##fun_masking: function to detect the influential points  using the max-statistics
##Usage:fun_masking(X,Y,n,p,q,n_subset,subset_vol,clean_setv,alpha)
##input:
#clean_setv: an input value of estimated clean set  obtained by Min/Max step
##output:
#clean_set: the indices of the estimated clean set obtained by Max-step
##########################################################



##########################################################
##fun_swamping: function to detect the influential points by using the min-statistics
##Usage: fun_swamping(X,Y,n,p,q,n_subset,subset_vol,clean_setv,ep=0.1,alpha)
##input:
#clean_setv: an input value of estimated clean set  obtained by Min/Max step
##output:
#clean_setv: the indices of the  estimated clean set  obtained  by the Min-step (obtained by  the min-statistics)
##########################################################



##########################################################
##fun_checking:function for checking whether there are noninfluentila points being identified as influential ones
##Usage: fun_checking(X,Y,n,p,q,inf_t,clean_t,alpha)
##input:
#inf_t: the estimated indices of influential poins found by Min-Max algorithm
#clean_t: the estimated  indices of clean poins found by Min-Max algorithm
##output:
#inf_setfinal: the estimated indices of influential points  obtained by MIP algorithm, after applying the checking algorithm to   the potential influential point inf_t.
##########################################################


fun_pv=function(X,Y,n,p,q,n_subset,subset_vol,clean_setv)
{
  rob_sd=apply(X,2,mad)
  X=X-rep(1,n)%o%apply(X,2,median);
  X=X%*%diag(1/rob_sd);
  if (q==1)
  {Y=(Y-median(Y))/mad(Y)}
  if (q>1)
  {
    y_sig=apply(Y,1,mad)
    Y_sigma_D=diag(y_sig,q,q)
    Y_sigma_R=sin(pi/2*(cor(t(Y),method='kendall')))
    Y_sigma=Y_sigma_D%*%Y_sigma_R%*%Y_sigma_D
    sY_sigma=eigen(Y_sigma)
    YY_sigma=sY_sigma$vectors%*%diag(1/(sqrt(sY_sigma$values)))%*%t(sY_sigma$vectors)
    Y=YY_sigma%*%(Y-median(Y))
  }

  TT=rep(0,n_subset)
  Tmax=rep(0,length(clean_setv))
  Tmin=rep(0,length(clean_setv))
  for (i in 1:length(clean_setv))
  {
    S_i=setdiff(clean_setv,clean_setv[i])

    for (m in 1: n_subset)
    {
      I=sample(S_i,size =subset_vol,replace=FALSE,prob=NULL)
      X1=X[c(clean_setv[i],I),]
      Y1=Y[,c(clean_setv[i],I)]
      X2=X[I,]
      Y2=Y[,I]
      rhat1= Y1%*%X1/(subset_vol+1)
      rhat2=Y2%*%X2/subset_vol
      TT[m]=((subset_vol+1)^2)*(sum((rhat1-rhat2)^2)/p)
    }
    Tmax[i]=max(TT)
    Tmin[i]=min(TT)
  }

  list(Tmax=Tmax,Tmin=Tmin)
}


fun_masking=function(X,Y,n,p,q,n_subset,subset_vol,clean_setv,alpha)
{
  Tmax=(fun_pv(X,Y,n,p,q,n_subset,subset_vol,clean_setv))$Tmax
  pv=1-pchisq(Tmax,q)
  Spv=sort.int(pv,index.return=TRUE)  # sorted p value
  Si=Spv$ix
  dp=Spv$x-alpha*c(1:length(Tmax))/length(Tmax)

  In=which(dp<=0)    # BH procedure to control the error rate
  if (length(In)==0)
  {clean_set=clean_setv}
  else
  {
    rin=max(In)
    inf_set=clean_setv[Si[1:rin]]
    clean_set=setdiff(clean_setv,inf_set)
  }

  list(clean_set=clean_set)
}

fun_swamping=function(X,Y,n,p,q,n_subset,subset_vol,clean_setv,ep=0.1,alpha)
{
  Tmin=(fun_pv(X,Y,n,p,q,n_subset,subset_vol,clean_setv))$Tmin
  pvv=1-pchisq(Tmin,q)
  Spvv=sort.int(pvv,index.return=TRUE)
  Sii=Spvv$ix
  dpv=Spvv$x-alpha*c(1:length(Tmin))/length(Tmin)

  In=which(dpv<=0)
  if (length(In)==0)
  {clean_set=clean_setv}
  else
  {
    rin=max(In)
    inf_setv=clean_setv[Sii[1:min(floor(ep*n),rin)]]
    clean_set=setdiff(clean_setv,inf_setv)
  }

  clean_setv=clean_set
  list(clean_setv=clean_setv)
}


fun_checking=function(X,Y,n,p,q,inf_t,clean_t,alpha)
{
  rob_sd=apply(X,2,mad)
  X=X-rep(1,n)%o%apply(X,2,median);
  X=X%*%diag(1/rob_sd);
  if (q==1)
  {Y=(Y-median(Y))/mad(Y)}
  if (q>1)
  {
    y_sig=apply(Y,1,mad)
    Y_sigma_D=diag(y_sig,q,q)
    Y_sigma_R=sin(pi/2*(cor(t(Y),method='kendall')))
    Y_sigma=Y_sigma_D%*%Y_sigma_R%*%Y_sigma_D
    sY_sigma=eigen(Y_sigma)
    YY_sigma=sY_sigma$vectors%*%diag(1/(sqrt(sY_sigma$values)))%*%t(sY_sigma$vectors)
    Y=YY_sigma%*%(Y-median(Y))
  }
  T=rep(0,length(inf_t))
  for (i in 1:length(inf_t))
  {
    X1=X[c(inf_t[i],clean_t),]
    Y1=Y[,c(inf_t[i],clean_t)]
    X2=X[clean_t,]
    Y2=Y[,clean_t]
    rhat1=Y1%*%X1/(length(clean_t)+1)
    rhat2=Y2%*%X2/length(clean_t)
    T[i]=((length(clean_t)+1)^2)*(sum((rhat1-rhat2)^2)/p)
  }
	 pv_inf=1-pchisq(T,q)
  Spv_inf=sort.int(pv_inf,index.return=TRUE)
  Si=Spv_inf$ix
  dp=Spv_inf$x-alpha*c(1:length(inf_t))/length(inf_t)

  In=which(dp<=0)
  if (length(In)==0)
  {clean_setfinal=c(1:n)
  inf_setfinal=setdiff(c(1:n),clean_setfinal)
  }

  else{
    rin=max(In)
    inf_setfinal=inf_t[Si[1:rin]]
    clean_setfinal=setdiff(c(1:n),inf_setfinal)
  }

  list(inf_setfinal=inf_setfinal)
}

MIP=function(X,Y,n,p,q,n_subset,subset_vol,ep=0.1,alpha)
{
  clean_setv=c(1:n)
  clean_setv=fun_swamping(X,Y,n,p,q,n_subset,subset_vol,clean_setv,ep,alpha)$clean_setv
  clean_set=fun_masking(X,Y,n,p,q,n_subset,subset_vol,clean_setv,alpha)$clean_set
  while (length(clean_set)<n/2)
  {
  clean_setv=fun_swamping(X,Y,n,p,q,n_subset,subset_vol,clean_setv,ep,alpha)$clean_setv
  clean_set=fun_masking(X,Y,n,p,q,n_subset,subset_vol,clean_setv,alpha)$clean_set
  }
  clean_t=clean_set
  inf_t=setdiff(c(1:n),clean_t)
  inf_setfinal=fun_checking(X,Y,n,p,q,inf_t,clean_t,alpha)$inf_setfinal
  list(inf_setfinal=inf_setfinal)

}








