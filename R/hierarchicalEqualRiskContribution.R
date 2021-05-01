computeClusterInertia <- function(labels, asset_returns) {
  inertia <- 0
  n <- nrow(asset_returns)
  for(label in unique(labels)){
    # removing `as.matrix` increases the speed significantly
    # inertia <- inertia + mean(as.matrix(dist(asset_returns[,labels==label])))
    inertia <- inertia + sum(dist(asset_returns[,labels==label])) * 2 / n^2
  }
  return(log(inertia))
} 

computeExpectedInertia <- function(N, D, num_clusters, method,
                                   num_reference_datasets=5) {
  reference_inertia <- 0
  for(i in 1:num_reference_datasets){
    reference_asset_returns = replicate(N, runif(D))
    reference_corr <- cor(reference_asset_returns)
    reference_dist <- as.dist(sqrt(2 * (1 - reference_corr)))

    if(method=='divisive'){
      reference_clust <- cluster::diana(reference_dist)
    }else{
      reference_clust<- hclust(reference_dist, method=method)
    }

    reference_cluster_assignments <- dendextend::cutree(reference_clust,
                                                        k=num_clusters,
                                                        order_clusters_as_data=F)

    inertia <- computeClusterInertia(reference_cluster_assignments,
                                     reference_asset_returns)
    reference_inertia <- reference_inertia + inertia
  }
  return(reference_inertia/num_reference_datasets)
}

computeNumClusters <- function(rho, method, asset_returns, distance) {
  gap_values <- list()
  num_clusters <- 1
  max_clusters <- -Inf
  N <- ncol(asset_returns)
  D <- nrow(asset_returns)

  if(method=='divisive'){
    original_clusters <- cluster::diana(distance)
  }else{
    original_clusters <- hclust(distance, method=method)
  }

  for(i in 1:nrow(rho)){
    original_cluster_assignments <- dendextend::cutree(original_clusters,
                                                       k=num_clusters,
                                                       order_clusters_as_data=F)

    if(max(original_cluster_assignments)==max_clusters |
       max(original_cluster_assignments)>10){
      break
    }
    max_clusters = max(original_cluster_assignments)
    inertia <- computeClusterInertia(original_cluster_assignments, asset_returns)
    expected_inertia <- computeExpectedInertia(N, D, num_clusters, method)

    gap_values <- c(gap_values, expected_inertia - inertia)
    num_clusters <- num_clusters + 1
  }
  return(which.max(gap_values))
}

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

calculateClusterRiskContribution <- function(risk_measure, Sigma, asset_returns, 
                                             num_clusters, cluster_assignments){
  cluster_contributions <- rep(1, num_clusters)

  for(i in 1:num_clusters){
    cluster_asset_ind <- which(cluster_assignments==i)

    switch(risk_measure,
           'variance'={
             cluster_contributions[i] <- getClusterVariance(Sigma, cluster_asset_ind)
           },
           'standard-deviation'={
             cluster_contributions[i] <- sqrt(getClusterVariance(Sigma, cluster_asset_ind))
           },
           'CVaR'={
             cluster_contributions[i] <- getClusterCVaR(asset_returns, cluster_asset_ind)
           },
           'CDaR'={
             cluster_contributions[i] <- getClusterCDaR(asset_returns, cluster_asset_ind)
           }
    )
  }
  return(cluster_contributions)
}

getRecursiveBisection <- function(dend, cluster_contributions, risk_measure){
  num_clusters <- nobs(dend)
  weights <- rep(1,num_clusters)

  computeChildClusterWeights <- function(dend, cluster_contributions, risk_measure) {
    left <- dend[[1]]
    right <- dend[[2]]
    left_order <- order.dendrogram(left)
    right_order <- order.dendrogram(right)

    if(risk_measure=='equal-weighting'){
      alpha <- 0.5
    }else {
      left_contrib <- sum(cluster_contributions[left_order])
      right_contrib <- sum(cluster_contributions[right_order])

      alpha <- right_contrib /(right_contrib+left_contrib)
    }

    weights[left_order] <<- weights[left_order] * alpha
    weights[right_order] <<- weights[right_order] * (1-alpha)

    if(nobs(left)>1)
      computeChildClusterWeights(left, cluster_contributions, risk_measure)
    if(nobs(right)>1)
      computeChildClusterWeights(right, cluster_contributions, risk_measure)
  }
  computeChildClusterWeights(dend, cluster_contributions, risk_measure)
  return(weights)
}


calculateFinalPortfolioWeights <- function(risk_measure, clusters_weights, Sigma,
                                           asset_returns, num_clusters, cluster_assignments) {
  w <- rep(1, nrow(Sigma))
  for(i in 1:num_clusters){
    cluster_inds <- cluster_assignments==i
    Sigma_ <- Sigma[cluster_inds, cluster_inds]
    cluster_asset_returns <- asset_returns[,cluster_inds]

    switch(risk_measure,
           'variance'={
             parity_weights <- getIVP(Sigma_)
           },
           'standard-deviation'={
             parity_weights <- getIVP(Sigma_)
           },
           'equal-weighting'={
             n <- sum(cluster_assignments==i)
             parity_weights <- rep(1/n, n)
           },
           'CVaR'={
             parity_weights <- getInverseCVaR(cluster_asset_returns)
           },
           'CDaR'={
             parity_weights <- getInverseCDaR(cluster_asset_returns)
           }
    )

    w[cluster_inds] <- parity_weights * clusters_weights[i]
  }
  return(w)
}

getPortfolioWeights <- function(hcluster, asset_returns, Sigma, num_clusters,
                                risk_measure, cluster_assignments) {
  N <- ncol(Sigma)

  clusters_weights <- rep(1, num_clusters)
  cluster_contributions <- calculateClusterRiskContribution(risk_measure,Sigma,
                                                            asset_returns,num_clusters,
                                                            cluster_assignments)

  dend_depth <- hcluster$height[N - num_clusters]
  clusters_dend <- cut(as.dendrogram(hcluster), h=dend_depth)$upper

  clusters_weights <- getRecursiveBisection(clusters_dend, cluster_contributions, risk_measure)

  w <- calculateFinalPortfolioWeights(risk_measure, clusters_weights, Sigma, 
                                      asset_returns, num_clusters, cluster_assignments)

  names(w) <- names(asset_returns)
  return(w)
}


#' @title Design of hierarchical equal risk contribution portfolios
#'
#' @description This function designs hierarchical equal risk contribution
#'              portfolios based on the method proposed by Raffinot (2018).
#'
#' @details This portfolio allocation method makes use of hierarchical clustering
#'          to assign portfolio weights.
#'
#' @param asset_prices An XTS object of the asset prices.
#' @param asset_returns An XTS object of the asset returns. 
#' @param Sigma Covariance matrix of returns. If none is provided, the
#'        covariance matrix will be computed from the returns.
#' @param risk_measure String indicating the desired risk measure for assigning
#'        portfolio weights. Must be one of c('variance', 'standard-deviation',
#'        'equal-weighting', 'CVaR', 'CDaR') 
#' @param method String indicating the desired hierarchical clustering method.
#'        Must be one of c("single", "complete", "average" ,"ward.D", "ward.D2",
#'        "divisive"). If method="divisive", divisive clustering (or the DIANA
#'        algorithm)is used, otherwise agglomerative clustering is used with
#'        method referring to the desired linkage function.
#' @param num_clusters Integer value representing the optimal number of clusters.
#'        If no value is given, the optimal number of clusters will be computed 
#'        automatically.
#'
#' @export
hierarchicalEqualRiskContribution <- function(asset_prices=NULL,
                                              asset_returns=NULL, Sigma=NULL,
                                              risk_measure=c('variance', 
                                                             'standard-deviation', 
                                                             'equal-weighting',
                                                             'CVaR', 'CDaR'), 
                                              method='ward.D2', num_clusters=NULL) {

  risk_measure <- match.arg(risk_measure)
  if(is.null(asset_prices) && is.null(asset_returns))
    stop("Asset prices and returns not given.")
  if(is.null(asset_returns))
    asset_returns <- diff(log(asset_prices))[-1]
  if(is.null(Sigma))
    Sigma <- cov(asset_returns)
  rho <- cov2cor(Sigma)
  distance <- as.dist(sqrt(2 * (1-rho)))

  if(method=='divisive'){
    hcluster <- cluster::diana(distance)
  }else{
    hcluster <- hclust(distance, method=method)
  }

  if(is.null(num_clusters))
    num_clusters <- computeNumClusters(rho, method, asset_returns, distance)

  cluster_assignments <- dendextend::cutree(hcluster, k=num_clusters, order_clusters_as_data=F)
  asset_names <- names(asset_prices)
  cluster_assignments <- cluster_assignments[asset_names]

  w <- getPortfolioWeights(hcluster, asset_returns, Sigma, num_clusters,
                           risk_measure, cluster_assignments)

  return(list('w'=w))
}

