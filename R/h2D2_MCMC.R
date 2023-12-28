#' @title MCMC sampling under h2-D2 prior.
#' 
#' @description
#' Using MCMC sampler to fit posterior distribution under h2-D2 prior.
#' 
#' @usage
#' h2D2 = h2D2_MCMC(h2D2, mcmc_n = 10000, burn_in = 5000, 
#'                  thin = 1, stepsize = 2, seed = 428)
#' 
#' @param h2D2 An h2D2 object.
#' @param mcmc_n Number of posterior samples to generate.
#' @param burn_in Number of early samples to discard.
#' @param thin Thinning parameter of the MCMC chain. The default is 1, 
#' indicating no thinning.
#' @param stepsize Step size for proposal.
#' @param seed Random seed for MCMC sampling.
#' 
#' @return An h2D2 object with MCMC samples. See \code{\link{h2D2-class}}.
#' 
#' @import Matrix
#' @export

h2D2_MCMC = function(h2D2, mcmc_n = 100, burn_in = 0, thin = 1,
                     stepsize = 2, seed = 428)
{
  if(!is.numeric(mcmc_n) | mcmc_n < 0)
  {
    stop("mcmc_n should be a non-negative integer.")
  }
  mcmc_n = round(mcmc_n)
  
  if(!is.numeric(burn_in) | burn_in < 0)
  {
    stop("burn_in should be a non-negative integer.")
  }
  burn_in = round(burn_in)
  
  if(!is.numeric(thin) | thin <= 0)
  {
    stop("thin should be a positive integer.")
  }
  thin = ceiling(thin)
  
  if(!is.numeric(stepsize) | stepsize <= 0)
  {
    stop("stepsize should be a positive number.")
  }
  
  if(!is.numeric(seed) | seed < 0)
  {
    stop("seed should be a non-negative integer.")
  }
  seed = round(seed)
  
  if(h2D2@trait == "quantitative")
  {
    scale2 = 1
    S_2 = h2D2@N
  } else if(h2D2@trait == "binary") {
    scale2 = h2D2@N1 * h2D2@N0 / h2D2@N^2
    S_2 = (h2D2@N1 * h2D2@N0) / h2D2@N
  }
  
  W = h2D2@R * S_2
  W = as(W, "dgCMatrix")
  
  S_2betaHat = h2D2@betaHat * S_2
  
  h2D2_sampling(h2D2, list(W, rep(S_2, h2D2@M), S_2betaHat), 
                mcmc_n, thin, stepsize, scale2, seed)
  
  #burn in
  h2D2@mcmc_samples[["n_burnin"]] = h2D2@mcmc_samples[["n_burnin"]] + burn_in
  if(burn_in > 0)
  {
    for(l in c("beta", "log_sigma2", "psi"))
    {
      h2D2@mcmc_samples[[l]] = h2D2@mcmc_samples[[l]][-c(1:burn_in),]
    }
    
    h2D2@mcmc_samples$h2 = h2D2@mcmc_samples$h2[-c(1:burn_in)]
    h2D2@mcmc_samples$h2_beta = h2D2@mcmc_samples$h2_beta[-c(1:burn_in)]
  }
  for(l in c("beta", "log_sigma2", "psi"))
  {
    colnames(h2D2@mcmc_samples[[l]]) = h2D2@SNP_ID
  }
  h2D2@mcmc_samples$n_samples = nrow(h2D2@mcmc_samples$beta)
  
  h2D2@mcmc_samples$h2_beta = h2D2@mcmc_samples$h2_beta * scale2
  
  #summary statistics
  sigma2 = exp(h2D2@mcmc_samples$log_sigma2) / scale2
  h2D2@mcmc_mean$sigma2 = colMeans(sigma2)
  h2D2@mcmc_sd$sigma2 = apply(sigma2, 2, sd)
  h2D2@mcmc_quantile$sigma2 = apply(sigma2, 2, quantile, 
                                    probs = c(0.025,0.25,0.5,0.75,0.975))
  names(h2D2@mcmc_mean$sigma2) = names(h2D2@mcmc_sd$sigma2) = 
    colnames(h2D2@mcmc_quantile$sigma2) = h2D2@SNP_ID
  
  h2D2@mcmc_mean$beta = colMeans(h2D2@mcmc_samples$beta)
  h2D2@mcmc_sd$beta = apply(h2D2@mcmc_samples$beta, 2, sd)
  h2D2@mcmc_quantile$beta = apply(h2D2@mcmc_samples$beta, 2, quantile, 
                                  probs = c(0.025,0.25,0.5,0.75,0.975))
  names(h2D2@mcmc_mean$beta) = names(h2D2@mcmc_sd$beta) = 
    colnames(h2D2@mcmc_quantile$beta) = h2D2@SNP_ID
  
  for(l in c("h2", "h2_beta"))
  {
    h2D2@mcmc_mean[[l]] = mean(h2D2@mcmc_samples[[l]])
    h2D2@mcmc_sd[[l]] = sd(h2D2@mcmc_samples[[l]])
    h2D2@mcmc_quantile[[l]] = quantile(h2D2@mcmc_samples[[l]], 
                                       probs = c(0.025,0.25,0.5,0.75,0.975))
  }
  
  h2D2@CL = apply(h2D2@mcmc_samples$beta, 2, function(x){
    abs(sum(x > 0) - sum(x < 0)) / length(x)})
  h2D2@CS = h2D2_CS(h2D2, h2D2@coverage, h2D2@purity)
  
  # Convergence
  n = h2D2@mcmc_samples$n_samples
  n_2 = floor(n/2)
  h2D2@mcmc_samples$PSRF_beta = foreach(j = 1:h2D2@M, .combine = "c") %do%
  {
    s1 = h2D2@mcmc_samples$beta[1:n_2,j]
    s2 = h2D2@mcmc_samples$beta[(n-n_2+1):n,j]
    B = n_2 * (mean(s1) - mean(s2))^2 / 2
    W = (var(s1) + var(s2))/2
    sqrt(1 - 1/n_2 + B/(W+5e-324) / n_2)
  }
  if(any(h2D2@mcmc_samples$PSRF_beta > 1.1))
  {
    warning(sprintf("The MCMC chain may not converge well for beta %s.",
                    paste(which(h2D2@mcmc_samples$PSRF_beta > 1.1), 
                          collapse = ',')))
  }
  
  return(h2D2)
}