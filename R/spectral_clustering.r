# spectral clustering 

################################################################################
# affinity matrix A 

affinity.mat <- function(dat){
  
  n <- nrow(dat)
  
  A <- matrix(0, nrow=n, ncol=n)
  
  for(i in 1:n) {
    for(j in 1:n) {
      if(i!=j){
        A[i,j] <- exp( (-1) * norm(as.matrix(dat[i,]-dat[j,]), type="F"))
      }else{
        A[i,j] <- 0
      }
    }
  }
  return(A)
}

################################################################################
# formation of the matrix L=D^(-1/2) %*% A %*% D^(-1/2)

spectral_clustering <- function(dat,k){
  
  dat <- as.matrix(dat)
  
  A <- affinity.mat(dat)
  
  D <- diag(apply(A, 1, sum))
  
  Q <- diag(diag(D)^(-1/2))

  L <- Q %*% A %*% Q 

  spec.dat <- as.matrix(scale(eigen(L)$vector[,1:k], center=F, scale=TRUE))
  
  return(spec.dat)
  
}
  

