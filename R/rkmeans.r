
# Wei Sun and Junhui Wang (2012)
# Regularized K-means Clustering of High-dimensional Data and Its Asymptotic Consistency
# Electronic Journal of Statistics, Volume 6, 148-167

eps=1e-3

################################################################################
###Functions######
mode=function(x) return(-as.numeric(names((sort(-table(-x)))[1])))

max.loc=function(x,rev=T) {if (rev) return(rev(which(x==max(x)))[1]) else return((which(x==max(x)))[1])}

min.loc=function(x) {return(rev((which(x==min(x))))[1])}

std=function(x) { return((x-mean(x))/sd(x)) }

ksearch=function(vec,ctr,k) {
	# find the cluster label of an input vec, given cluster center ctr and number of clusters k.

    ctr.mat=matrix(ctr,nrow=k,byrow=F)
    ctr.new=as.vector(t(ctr.mat))
    sq.diff=(rep(vec,k)-ctr.new)^2
    p=length(ctr)/k
    ind=rep(1:k, rep(p,k))
    kdist=tapply(sq.diff, ind, sum)
    return(order(kdist)[1])
}


clus.dist=function(clus1,clus2) {
	#Stability Distance
	if (length(clus1)!=length(clus2)) return("cluster sizes don't match!")
	n=length(clus1)
	s=0
	for ( i in 2:n ) {
		for ( j in 1:(i-1)) {
			s=s+as.numeric(clus1[i]==clus1[j] & clus2[i]!=clus2[j])
			s=s+as.numeric(clus1[i]!=clus1[j] & clus2[i]==clus2[j])
		}
	}
	return(2*s/(n*(n-1)))
}

#function of f(x)-normal
f.density=function(x.vec,c.vec,sigma.vec)
  {
    p=length(x.vec)
    V=diag(sigma.vec)
    return(as.numeric(1/(((2*pi)^(p/2))*(det(V))^2)*exp((-1/2*t(x.vec-c.vec)%*%solve(V)%*%(x.vec-c.vec)))))
  }

##log-likelihood function for complete data
log.ld=function(x,k,l.mat,pi.vec,sigma.vec,c.mat,c.mat.km,lambda)
  {
   n=nrow(x)
   p=ncol(x)
   l1=matrix(0,n,k)
   l2=rep(0,p)
   for(j in 1:p)  {
      for(k.ind in 1:k)  {
        for(i in 1:n)  l1[i,k.ind]=l.mat[i,k.ind]*(log(pi.vec[k.ind])+log(f.density(x[i,],c.mat[k.ind,],sigma.vec)))
       }
      l2[j]=sqrt(t(c.mat[,j])%*%c.mat[,j])/(sqrt(t(c.mat.km[,j])%*%c.mat.km[,j]))
     }
   return(sum(l1)-lambda*sum(l2))
  }



################################################################################
##Standard Kmeans###
my.kmeans=function(data,k) {
       p=ncol(data)
       kkm=kmeans(data,k,nstart=20)
       return( list( center=kkm$centers, activeset=1:p, membership=kkm$cluster))
}

################################################################################
###Regularized Kmeans with adaptive group lasso#####

adp.skmeans=function(data,k,lambda,itn=5)
    {
       n=nrow(data)
       p=ncol(data)
       x=data

       kkm=kmeans(x,k,nstart=100)
       c.mat=kkm$centers
       c.mat.old=matrix(0,k,p)
       lambda.old=lambda
       ii=0
       while(sum(abs(c.mat-c.mat.old))/(sum(abs(c.mat))+eps)>=eps&ii<itn)
            {  ii=ii+1
               c.mat.old=c.mat
               l.vec=apply(x,1,ksearch,k=k,ctr=as.vector(c.mat))
               if(length(unique(l.vec))<k)  break
               ll.mat=matrix(0,n,k)
               for(jj in 1:n) ll.mat[jj, l.vec[jj]]=1
               for(j in 1:p) {
                     tmp.dem=sum(abs(c.mat[,j]))
                     if (tmp.dem<=eps) tmp.dem=eps
                     m.temp=grplasso::grplasso(x=ll.mat,y=x[,j],index=rep(1,k),model=LinReg(),lambda=n*lambda.old/tmp.dem,center=F,control=grpl.control(max.iter=10,trace=0))
                     c.mat[,j]=m.temp$coef
                 }
               if (sum(abs(c.mat))<eps) break
            }
    activeset=which(apply(abs(c.mat),2,sum)>=eps)
    return( list( center=c.mat, activeset=activeset, membership=l.vec))
}


################################################################################
##Model based clustering: EM algorithm###
## In Sun et al. (EJS, 2012), there is a typo on page 9, the updating for $\widehat{C}_{(j)}^{(t+1)}$.
## The derivative of Q w.r.t. C_{kj} should be \sum_{i=1}^n L_{ik}^{(t)} [ X_{ij} - C_{kj} ] / \sigma_{j}^2 - \lambda_j C_{kj} = 0
## Therefore, C_{kj} = [ \sum_{i=1}^n L_{ik}^{(t)} X_{ij} ] / [ \sum_{i=1}^n L_{ik}^{(t)} / \sigma_{j}^2 + \lambda_j ]
## The following "mbc" needs to be updated
################################################################################

mbc=function(data,k,lambda,itn=5)
    {
       n=nrow(data)
       p=ncol(data)
       x=data

       kkm=kmeans(x,k,nstart=100)
       c.mat=kkm$centers
       c.mat.km=c.mat
       l.vec=apply(x,1,ksearch,k=k,ctr=as.vector(c.mat))
       ll.mat=matrix(0,n,k)
       for(jj in 1:n) ll.mat[jj, l.vec[jj]]=1
       sigma.vec=diag(var(x))
       pi.vec=rep(1/k,k)

       c.mat.old=matrix(0,k,p)
       ll.mat.old=matrix(0,n,k)
       sigma.vec.old=rep(0,p)
       pi.vec.old=rep(0,k)
       temp=matrix(0,n,k)
       temp.f=matrix(0,n,k)
       temp2=matrix(0,p,k)

       ii=0
       while(sum(abs(c.mat-c.mat.old))/sum(abs(c.mat))>=eps&ii<itn)
            {  ii=ii+1
               l.vec=apply(x,1,ksearch,k=k,ctr=as.vector(c.mat))
               if(length(unique(l.vec))<k)  break
               c.mat.old=c.mat
               ll.mat.old=ll.mat
               sigma.vec.old=sigma.vec
               pi.vec.old=pi.vec



               for(k.ind in 1:k) {
                    for(i in 1:n) {
                     temp.f[i,k.ind]=f.density(x[i,],c.mat.old[k.ind,],sigma.vec.old)
                     }
                 }

              for(i in 1:n){ if(sum(temp.f[i,])==0) temp.f[i,]=ll.mat.old[i,]}


              for(j in 1:p){
                   for(k.ind in 1:k) {
                      for(i in 1:n) {
                       if(pi.vec.old%*%temp.f[i,]==0) ll.mat[i,]=ll.mat.old[i,] else
                       ll.mat[i,k.ind]=pi.vec.old[k.ind]*temp.f[i,k.ind]/pi.vec.old%*%temp.f[i,]
                       temp[i,k.ind]=ll.mat[i,k.ind]*(x[i,j]-c.mat.old[k.ind,j])^2/n
                       }
                   pi.vec[k.ind]=mean(ll.mat[,k.ind])
                   temp2[j,k.ind]=(ll.mat[,k.ind]%*%x[,j]/sum(ll.mat[,k.ind]))^2
                   }
                sigma.vec[j]=sum(temp)
               }


              for(j in 1:p){
                  lambda.j=n*lambda/sqrt(c.mat.km[,j]%*%c.mat.km[,j])
                   for(k.ind in 1:k) {
                     if(lambda.j>=sqrt(sum(temp2[j,]))*sum(ll.mat[,k.ind])/sigma.vec[j]) c.mat[k.ind,j]=0 else
                    c.mat[k.ind,j]=(1-lambda.j*sigma.vec[j]/(sqrt(sum(temp2[j,]))*sum(ll.mat[,k.ind])))*(ll.mat[,k.ind]%*%x[,j]/sum(ll.mat[,k.ind]))

                   }
               }

              if (sum(abs(c.mat))<eps) break


            }


    activeset=which(apply(abs(c.mat),2,sum)>=eps)
    return( list( center=c.mat, activeset=activeset, membership=l.vec))
}


################################################################################
##Choice of K and \lambda for adaptive regularized kmeans####

clus.bt=function(data,clus=adp.skmeans,k.max=10,bt.num=10,lambda.num=20){
	x=data.matrix(data)
	n=nrow(x)
	k.s=rep(0,bt.num)
	tmp.s=array(0,c(bt.num,k.max-1,lambda.num))
	res.k=matrix(0,bt.num,k.max-1)
	lambda.vec=10^seq(-2,2,length=lambda.num)

	for(centers in 2:k.max) {
        for(l.ind in 1:lambda.num){
           tr1.clus=clus(x,k=centers,lambda=lambda.vec[l.ind])

           for(bt.i in 1:bt.num){
				x.ind=sample(1:n,n,replace=T)
				x.bt=x[x.ind,]

				tr2.clus=clus(x.bt,k=centers,lambda=lambda.vec[l.ind])

				tr1.val=apply(x.bt,1,ksearch,ctr=tr1.clus$center,k=centers)
				tr2.val=apply(x.bt,1,ksearch,ctr=tr2.clus$center,k=centers)

				if(length(unique(tr1.val))+length(unique(tr2.val))<2*centers)  {tmp.s[bt.i,centers-1,l.ind]==1;quit}
				if(min(table(tr1.val))<2||min(table(tr2.val))<2)  {tmp.s[bt.i,centers-1,l.ind]==1;quit}
				tmp.s[bt.i,centers-1,l.ind]=ifelse(sum(abs(tr1.clus$center))*sum(abs(tr2.clus$center))==0, 1, clus.dist(tr1.val,tr2.val))
      	    }
       }
       for (bt.i in 1:bt.num) res.k[bt.i,centers-1]=min(tmp.s[bt.i,centers-1,])
	}
	k.best=mode(apply(res.k,1,min.loc))
	l.best=mode(apply(tmp.s[,k.best,],1,min.loc))

	return(list(k=k.best+1,lambda=lambda.vec[l.best]))
}


