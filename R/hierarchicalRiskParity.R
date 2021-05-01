getIVP <- function(Sigma) {
  if(is.null(nrow(Sigma))) return(1) 
  ivp <- 1./diag(Sigma)
  ivp <- ivp/sum(ivp)
  return(ivp)
}

getClusterVar <- function(Sigma, cluster_idx) {
  Sigma_ <- Sigma[cluster_idx, cluster_idx]
  w_ <- getIVP(Sigma_)
  cluster_var <- t(w_) %*% Sigma_ %*% w_
  return(cluster_var[1])
}

bisectClusters <- function(clusters){
  clusters_new <- list()
  for(cl in clusters){
    if(length(cl)<2) next
    n <- length(cl) %/% 2
    len <- length(cl)
    clusters_new <- c(clusters_new, list(cl[1:n]), list(cl[(n+1):len]))
  }
  return(clusters_new)
}

getRecBipart <- function(Sigma, sorted_idx, w_min, w_max, lam) {
  N <- length(sorted_idx)
  w <- rep(1., N)
  clusters <- bisectClusters(list(sorted_idx))

  while(length(clusters) > 0) {
    for(i in 1:length(clusters)){
      if(i%%2==0 || length(clusters)==1) next

      cl0 <- unlist(clusters[i])
      cl1 <- unlist(clusters[i+1])

      cl_var0 <- lam * getClusterVar(Sigma, cl0)
      cl_var1 <- getClusterVar(Sigma, cl1)
      alpha <- 1 - cl_var0 / (cl_var0 + cl_var1)

      alpha <- min(sum(w_max[cl0]) / w[cl0[1]],
                   max(sum(w_min[cl0]) / w[cl0[1]],
                       alpha))

      alpha <- 1- min(sum(w_max[cl1]) / w[cl1[1]],
                   max(sum(w_min[cl1]) / w[cl1[1]],
                       1-alpha))

      w[cl0] <- w[cl0] * alpha
      w[cl1] <- w[cl1] * (1-alpha)
    }
    clusters <- bisectClusters(clusters)

  }
  names(w) <- rownames(Sigma)
  return(w)
}

#' @title Design of hierarchical risk parity portfolios
#'
#' @description This function designs a hierarchical risk parity portfolio
#' based on Lopez de Prado's paper "Building Diversified Portfolios that
#' Outperform Out-of-Sample" (2015).
#' 
#' @details This portfolio allocation method is a heuristic method which makes
#' use of hierarchical clustering to seriate the correlation matrix of asset
#' returns, then recursively allocates the portfolio weights via naive risk
#' parity and "recursive bisection".
#'
#' @param asset_prices An XTS object of the asset prices.
#' @param asset_returns An XTS object of the asset returns. 
#' @param Sigma Covariance matrix of returns. If none is provided, the
#'        covariance matrix will be computed from the returns.
#' @param method String indicating the desired hierarchical clustering method.
#'        Must be one of c("single", "complete", "average" ,"ward.D", "ward.D2",
#'        "divisive"). If method="divisive", divisive clustering (or the DIANA
#'        algorithm)is used, otherwise agglomerative clustering is used with
#'        method referring to the desired linkage function.
#' @param w_min Scalar or vector with values between [0,1] to control the
#'        minimum value of weights.
#' @param w_max Scalar or vector with values between [0,1] to control the
#'        maximum value of weights.
#' @param lam Non-negative tuning parameter to control the concentration into 
#'        different clusters.
#'
#' @export
hierarchicalRiskParity <- function(asset_prices=NULL, asset_returns=NULL,
                                   Sigma=NULL, method='single', w_min=NULL,
                                   w_max=NULL, lam=1) {
  if(is.null(Sigma) && is.null(asset_prices) && is.null(asset_returns))
    stop("Invalid input. Please provide either a covariance matrix, asset returns or asset prices.")
  if(is.null(asset_returns))
    asset_returns <- PerformanceAnalytics::Return.calculate(asset_prices,
                                                            method='arithmetic')[-1]
  if(is.null(Sigma))
    Sigma <- cov(asset_returns)
  rho <- cov2cor(Sigma)
  distance <- as.dist(sqrt((1-rho)/2))

  if(method=='divisive'){
    hcluster <- cluster::diana(distance)
  }else{
    hcluster <- hclust(distance, method=method)
  }

  sorted_idx <- hcluster$order

  if(is.null(w_min))
    w_min = rep(0, ncol(Sigma))
  if(is.null(w_max))
    w_max = rep(1, ncol(Sigma))

  if(!is.vector(w_min) | length(w_min)==1)
    w_min <- rep(w_min, ncol(Sigma))

  if(!is.vector(w_max) | length(w_max)==1)
    w_max <- rep(w_max, ncol(Sigma))

  if(sum(w_min>1))
    stop("Invalid minimum weights.")

  w <- getRecBipart(Sigma, sorted_idx, w_min, w_max, lam)

  return(list('w'=w))
}
