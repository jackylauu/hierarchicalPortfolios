getIVP <- function(Sigma) {
  if(is.null(nrow(Sigma))) return(1) 
  ivp <- 1./diag(Sigma)
  ivp <- ivp/sum(ivp)
  return(ivp)
}

getInverseCVaR <- function(asset_returns) {
  if(ncol(asset_returns)==1) return(1)
  w <- PerformanceAnalytics::CVaR(asset_returns, method="historical")
  w <- 1/w
  w <- w / sum(w)
  return(w)
}

getInverseCDaR <- function(asset_returns) {
  if(ncol(asset_returns)==1) return(1)
  w <- PerformanceAnalytics::CDD(asset_returns)
  w <- 1/w
  w <- w / sum(w)
  return(w)
}

getClusterVariance <- function(Sigma, cluster_idx) {
  Sigma_ <- Sigma[cluster_idx, cluster_idx]
  w_ <- getIVP(Sigma_)
  cluster_var <- t(w_) %*% Sigma_ %*% w_
  return(cluster_var[1])
}

getClusterCVaR <- function(asset_returns, cluster_idx) {
  conf <- 0.05
  cluster_returns <- asset_returns[,cluster_idx]

  w <- getInverseCVaR(cluster_returns)

  # portfolio returns is then an xts
  portfolio_returns <- w %*% t(cluster_returns) 
  portfolio_returns <- t(portfolio_returns)

  cluster_cvar <- PerformanceAnalytics::CVaR(portfolio_returns, p=conf, method="historical")

  # weights doesn't seem to work for PerformanceAnalytics
  # cluster_cvar <- PerformanceAnalytics::CVaR(cluster_returns, weights=w, p=conf, method="historical")
  return(cluster_cvar[1])
}

getClusterCDaR <- function(asset_returns, cluster_idx) {
  conf <- 0.05
  cluster_returns <- asset_returns[,cluster_idx]

  w <- PerformanceAnalytics::CDD(cluster_returns, p=conf)
  w <- 1/w
  w <- w / sum(w)

  # portfolio returns is then an xts
  portfolio_returns <- w %*% t(cluster_returns) 
  portfolio_returns <- t(portfolio_returns)

  cluster_cdd <- PerformanceAnalytics::CDD(portfolio_returns, p=conf)

  return(cluster_cdd[[1]])
}

getClusterRisk <- function(Sigma, asset_returns, inds, risk_measure){
  switch(risk_measure,
         'variance'={
           return(getClusterVariance(Sigma, inds))
         },
         'standard-deviation'={
           return(sqrt(getClusterVariance(Sigma, inds)))
         },
         'CVaR'={
           return(getClusterCVaR(asset_returns, inds))
         },
         'CDaR'={
           return(getClusterCDaR(asset_returns, inds))
         }
  )
}

getRecursiveBisection <- function(dend, Sigma, asset_returns,
                                  risk_measure, w_min, w_max, lam){
  w <- rep(1,nrow(Sigma))

  computeChildClusterWeights <- function(dend, Sigma, asset_returns,
                                         risk_measure, w_min, w_max,
                                         lam) {
    left <- dend[[1]]
    right <- dend[[2]]
    cl0 <- unlist(left)
    cl1 <- unlist(right)

    left_contrib <-  lam * getClusterRisk(Sigma, asset_returns, cl0, risk_measure)        
    right_contrib <- getClusterRisk(Sigma, asset_returns, cl1, risk_measure)

    alpha <- right_contrib /(right_contrib+left_contrib)
    alpha <- min(sum(w_max[cl0]) / w[cl0[1]],
                 max(sum(w_min[cl0]) / w[cl0[1]], 
                     alpha))

    alpha <- 1- min(sum(w_max[cl1]) / w[cl1[1]],
                    max(sum(w_min[cl1]) / w[cl1[1]], 
                        1-alpha))
    w[cl0] <<- w[cl0] * alpha
    w[cl1] <<- w[cl1] * (1-alpha)

    if(nobs(left)>1)
      computeChildClusterWeights(left, Sigma, asset_returns, risk_measure,
                                 w_min, w_max, lam)
    if(nobs(right)>1)
      computeChildClusterWeights(right, Sigma, asset_returns, risk_measure,
                                 w_min, w_max, lam) }
  computeChildClusterWeights(dend, Sigma, asset_returns, risk_measure, w_min, w_max, lam)

  names(w) <- rownames(Sigma)
  return(w)
}

#' @title Nodal Hierarchical Risk Parity
#'
#' @description A modified version of Lopez de Prado's (2015) Hierarchical Risk
#' Parity (HRP) method which incorporates the hierarchical structure of the minimum
#' spanning tree obtained from hierarchical clustering.
#' 
#' @details This portfolio allocation method is an extension of the
#' Hierarchical Risk Parity (HRP) method which learns the minimum spanning tree
#' from historical asset returns, then recursively allocates the portfolio
#' weights via naive risk parity and "recursive bisection" according to the
#' tree structure.
#'
#' @param asset_prices An XTS object of the asset prices.
#' @param asset_returns An XTS object of the asset returns. 
#' @param Sigma Covariance matrix of returns. If none is provided, the
#'        covariance matrix will be computed from the returns.
#' @param risk_measure String indicating the desired risk measure for assigning
#'        portfolio weights. Must be one of c('variance', 'standard-deviation',
#'        'CVaR', 'CDaR') 
#' @param method String indicating the desired hierarchical clustering method.
#'        Must be one of c("complete", "single", "average" ,"ward.D", "ward.D2",
#'        "divisive"). If method="divisive", divisive clustering (or the DIANA
#'        algorithm)is used, otherwise agglomerative clustering is used with
#'        method referring to the desired linkage function. The default is
#'        "complete".
#' @param w_min Scalar or vector with values between [0,1] to control the
#'        minimum value of weights.
#' @param w_max Scalar or vector with values between [0,1] to control the
#'        maximum value of weights.
#' @param lam Non-negative tuning parameter to control the concentration into 
#'        different clusters.
#'
#' @export
NHRP <- function(asset_prices=NULL, asset_returns=NULL, Sigma=NULL,
                     risk_measure=c('variance', 'standard-deviation',
                                    'CVaR', 'CDaR'),
                     method='complete', w_min=NULL, w_max=NULL, lam=1) {
  risk_measure <- match.arg(risk_measure)
  if(is.null(asset_prices) && is.null(asset_returns))
    stop("Asset prices and returns not given.")
  if(is.null(asset_returns))
    asset_returns <- PerformanceAnalytics::Return.calculate(asset_prices, method='arithmetic')[-1]

  if(is.null(Sigma))
    Sigma <- cov(asset_returns)

  rho <- cov2cor(Sigma)
  distance <- as.dist(sqrt(2 * (1-rho)))

  hcluster <- hclust(distance, method=method)
 
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

  w <- getRecursiveBisection(as.dendrogram(hcluster), Sigma, asset_returns,
                             risk_measure, w_min, w_max, lam)

  return(list('w'=w))
}
