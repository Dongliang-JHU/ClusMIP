# Author: Dongliang Zhang (^)
#
# Supervisors: Professor Masoud Asgharian (*)
#              Professor Martin Lindquist (^)
#              Professor Mei-Cheng Wang (^)
#
# Current Affiliation: (1) Department of Mathematics and Statistics, McGill University, Montreal, Quebec, Canada (*)
#                      (2) Department of Biostatistics, Johns Hopkins University, Baltimore, Maryland, USA (^)
#
# Title of Project: Detection of Influential Observations on Variable Selection
#
# R Script Purpose: High-dimensional Clustering Assessment and Comparison
#
# Created on  : April 10, 2022
#
# Modified on : April 10, 2022
#               April 11, 2022
#               April 19, 2022
#               November 18, 2023
#               November 19, 2023
#               November 20, 2023
#               August 05, 2024
#
# Dependent R Scripts: spectral_clustering.r
#                      rkmeans.r (Sun et al. 2012)
#
# R Packages : stats (K-means)
#              LICORS (K-means++)
#              Rtsne (t-SNE)
#              grplasso (regularized K-means)
#
##########################################################################################################################################
# Function: cluster.est
#
# Reference:
#
# Args:
#   df           :  data frame (X and y)
#   nclust       :  number of cluster centers
#   cluster.type :  five type of clustering from "kmeans", "kmeans++", "tsne", "spectral" and "rkmeans"
#

#' Clustering
#'
#' @param df data frame
#' @param nclust number of cluster centers
#' @param cluster.type types of clustering schemes
#'
#' @return A list of positions of potentially clean and influential observations
#' @export
#'

cluster.est <- function(df, nclust, cluster.type){

  n <- nrow(df)

  p <- ncol(df)

  df <- as.matrix(df)

  ##############################################################################
  # K-means
  if(cluster.type=="kmeans"){
    fit.kmeans <- stats::kmeans(df, centers=nclust)
    membership <- fit.kmeans$cluster
  }

  ##############################################################################
  # K-means++
  if(cluster.type=="kmeans++"){
    #fit.kmeanspp <- flexclust::kcca(x=df, k=nclust, family=kccaFamily("kmeans"), control=list(initcent="kmeanspp"))
    #membership <- clusters(fit.kmeanspp)

    fit.kmeanspp <- LICORS::kmeanspp(df, k=nclust)
    membership <- fit.kmeanspp$cluster

  }

  ##############################################################################
  # t-SNE
  if(cluster.type=="tsne"){
    fit.tsne <- Rtsne::Rtsne(df,dims=nclust,perplexity=10)
    fit.tsne.kmeans <- stats::kmeans(fit.tsne$Y, centers=nclust)
    membership <- fit.tsne.kmeans$cluster
  }

  ##############################################################################
  # spectral clustering
  if(cluster.type=="spectral"){
    fit.sc <- stats::kmeans(spectral_clustering(df,k=2),nclust)
    membership <- fit.sc$cluster
  }

  ##############################################################################
  # regularized K-means
  if(cluster.type=="rkmeans"){
    tmp <- clus.bt(data=df,clus=adp.skmeans,k.max=nclust,bt.num=10,lambda.num=10)
    fit.rkmeans <- adp.skmeans(data=df,k=nclust,lambda=tmp$lambda,itn=5)
    membership <- fit.rkmeans$membership
  }

  ##############################################################################
  # sparse convex clustering
  #if(cluster.type=="scc"){
  #  n <- nrow(df)
  #  g1 <- 6
  #  g2 <- 0
  #  Gamma2.weight <- c(rep(0.5, true.p), rep(1,p-true.p))
  #  k_w <- 5    # Number of nearest neighbors
  #  phi <- 0.5  # scale of the kernel
  #  w <- dist_weight( t(df) / sqrt(p),phi, dist.type = "euclidean", p=2)
  #  w <- knn_weights(w,k=k_w,n)
  #  nu <- AMA_step_size(w,n) /2
  #  fit.scvxclust <- scvxclust(X=df, w=w, Gamma1=g1, Gamma2=g2, Gamma2_weight=Gamma2.weight, method="ama", nu=nu, max_iter=10000, tol_abs=1e-5)
  #  membership <- find_clusters(create_adjacency(fit.scvxclust$V[[1]],w,n))$size
  #}

  ##############################################################################

  cluster.level <- unique(membership)

  l.vec <- vector(length=length(cluster.level))

  for(i in 1:length(cluster.level)){
    l.vec[i] <- length(which(membership==cluster.level[i]))
  }

  clean.pos <- sort(which(membership==cluster.level[which.max(l.vec)]), decreasing=FALSE)

  infl.pos <- sort(setdiff(1:n,clean.pos), decreasing=FALSE)

  ##############################################################################

  returnList <- list("clean.pos"=clean.pos,
                     "infl.pos"=infl.pos)

  return(returnList)
}
