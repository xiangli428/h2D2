#' @title MCMC sampling under h2-D2 prior.
#' 
#' @description
#' Using MCMC sampler to fit posterior distribution under h2-D2 prior.
#' 
#' @usage
#' h2D2 = h2D2_MCMC(h2D2, mcmc_n = 5500, burn_in = 500, thin = 1, 
#'                  stepsize = 2, seed = 428, get_CS = TRUE, 
#'                  n_chain = 3, fold = 1.1,
#'                  pre_mcmc_n = 400, pre_burn_in = 200, pre_p = 0.05,
#'                  pre_maxiter = 10, pre_miniter = 2)
#' 
#' @param h2D2 An h2D2 object.
#' @param mcmc_n Number of posterior samples to generate.
#' @param burn_in Number of early samples to discard.
#' @param thin Thinning parameter of the MCMC chain. The default is 1, 
#' indicating no thinning.
#' @param stepsize Step size for proposal.
#' @param seed Random seed for MCMC sampling.
#' @param get_CS Logical. If get_CS = TRUE, detect CS after MCMC. Otherwise, 
#' CS will not be detected.
#' @param n_chain Number of parallel chains with different temperature.
#' @param fold Fold change of temperature between two neighbored parallel chains.
#' For instance, if \code{n_chain = 3} and \code{fold = 1.1}, the temperature
#' for the three parallel chains will be 1, 1.1, 1.21.
#' @param pre_mcmc_n Number of posterior samples to generate in each iteration
#' of pretraining. If b is given, the pretrain will not run.
#' @param pre_burn_in Number of early samples to discard in each iteration of
#' pretraining.
#' @param pre_p p-value threshold for the pretraining.
#' @param pre_maxiter Maximum number of pretraining iterations allowed.
#' @param pre_miniter Minimum number of pretraining iterations required.
#' 
#' @return An h2D2 object with MCMC samples. See \code{\link{h2D2-class}}.
#' 
#' @import Matrix
#' @export

h2D2_MCMC = function(h2D2, mcmc_n = 5500, burn_in = 500, thin = 1,
                     stepsize = 2, seed = 428, get_CS = T, 
                     n_chain = 3, fold = 1.1,
                     pre_mcmc_n = 400, pre_burn_in = 200, pre_p = 0.05,
                     pre_maxiter = 10, pre_miniter = 2)
{
  if(!is.numeric(mcmc_n) | mcmc_n < 0)
  {
    stop("'mcmc_n' should be a non-negative integer.")
  }
  mcmc_n %<>% round()

  if(!is.numeric(burn_in) | burn_in < 0)
  {
    stop("'burn_in' should be a non-negative integer.")
  }
  burn_in %<>% round()

  if(!is.numeric(thin) | thin <= 0)
  {
    stop("'thin' should be a positive integer.")
  }
  thin %<>% round()
  thin %<>% max(1)

  if(!is.numeric(stepsize) | stepsize <= 0)
  {
    stop("'stepsize' should be a positive number.")
  }

  if(!is.numeric(seed) | seed < 0)
  {
    stop("'seed' should be a non-negative integer.")
  }
  seed %<>% round()
  
  if(!is.numeric(pre_mcmc_n) | mcmc_n < 0)
  {
    stop("'pre_mcmc_n' should be a non-negative integer.")
  }
  pre_mcmc_n %<>% round()
  
  if(!is.numeric(pre_burn_in) | pre_burn_in < 0)
  {
    stop("'pre_burn_in' should be a non-negative integer.")
  }
  pre_burn_in %<>% round()
  
  if(pre_mcmc_n <= pre_burn_in)
  {
    stop("'pre_mcmc_n' should be larger than 'pre_burn_in'.")
  }
  
  N = h2D2@N
  M = h2D2@M
  z = h2D2@z
  lambda = h2D2@lambda
  if(lambda == 0)
  {
    W = h2D2@R * N
    W %<>% as("generalMatrix")
    mu = z * sqrt(N)
    dW = rep(N, M)
  } else {
    values = h2D2@values
    tvectors = h2D2@tvectors
    mu = t(tvectors) %*% ((tvectors %*% z) * values / (values + lambda)) *
      sqrt(N)
    W = t(tvectors) %*% (tvectors * values^2 / (values + lambda)) * N
    W = (W+t(W)) / 2
    dW = diag(W)
    diag(W) = 0
    W %<>% as("dgCMatrix")
  }
  
  LD_pairs = as.matrix(h2D2@R)^2
  LD_pairs[LD_pairs < h2D2@purity^2] = 0
  LD_pairs = t(LD_pairs / (rowSums(LD_pairs) + 5e-324))
  LD_pairs %<>% as("dgCMatrix")
  
  sample = list()
  
  if(h2D2@mcmc_samples$n_samples == 0)
  {
    sample[[1]] = matrix(0, M, n_chain)
    sample[[2]] = matrix(log(1e-4/M), M, n_chain)
    sample[[3]] = matrix(1, M, n_chain)
  } else {
    sample = h2D2@mcmc_samples$last_sample
    if(n_chain != ncol(sample[[1]]))
    {
      warning("'n_chain' should be consistent.")
      n_chain = ncol(sample[[1]])
    }
  }
  
  temp = fold^c(1:n_chain - 1)
  
  if(is.null(h2D2@b))
  {
    sample = h2D2_pretrain(h2D2, sample, dW, W, mu, LD_pairs, n_chain, temp, 
                           pre_mcmc_n, pre_burn_in, pre_p, 
                           pre_maxiter, pre_miniter, stepsize, seed)
  }
  
  if(mcmc_n > 0)
  {
    samples = h2D2_sampling(h2D2, sample, dW, W, mu, LD_pairs, n_chain, temp, 
                            mcmc_n, thin, stepsize, seed)
    
    # Add new MCMC samples
    
    rownames(samples[[1]]) = rownames(samples[[2]]) = rownames(samples[[3]]) =
      h2D2@SNP_ID
    
    colnames(samples[[1]]) = colnames(samples[[2]]) = colnames(samples[[3]]) =
      names(samples[[4]]) = names(samples[[5]]) = paste(
        "mcmc", (h2D2@mcmc_samples$n_samples+h2D2@mcmc_samples$n_burnin+1):
          (h2D2@mcmc_samples$n_samples+h2D2@mcmc_samples$n_burnin+mcmc_n),
        sep = '_')
    
    if(h2D2@mcmc_samples$n_samples == 0)
    {
      h2D2@mcmc_samples$beta = samples[[1]]
      h2D2@mcmc_samples$tau = samples[[2]]
      h2D2@mcmc_samples$psi = samples[[3]]
      h2D2@mcmc_samples$h2_beta = samples[[4]]
      h2D2@mcmc_samples$h2 = samples[[5]]
    } else {
      h2D2@mcmc_samples$beta %<>% cbind(samples[[1]])
      h2D2@mcmc_samples$tau %<>% cbind(samples[[2]])
      h2D2@mcmc_samples$psi %<>% cbind(samples[[3]])
      h2D2@mcmc_samples$h2_beta %<>% c(samples[[4]])
      h2D2@mcmc_samples$h2 %<>% c(samples[[5]])
    }
    
    h2D2@mcmc_samples$last_sample = samples[[6]]
  }
  
  #burn in
  
  h2D2@mcmc_samples$n_burnin %<>% add(burn_in)
  
  if(burn_in > 0)
  {
    for(l in c("beta", "tau", "psi"))
    {
      h2D2@mcmc_samples[[l]] = h2D2@mcmc_samples[[l]][,-c(1:burn_in)]
    }
    h2D2@mcmc_samples$h2_beta = h2D2@mcmc_samples$h2_beta[-c(1:burn_in)]
    h2D2@mcmc_samples$h2 = h2D2@mcmc_samples$h2[-c(1:burn_in)]
  }
  
  h2D2@mcmc_samples$n_samples = ncol(h2D2@mcmc_samples$beta)
  
  #summary statistics
  h2D2@mcmc_mean$beta = rowMeans(h2D2@mcmc_samples$beta)
  h2D2@mcmc_sd$beta = apply(h2D2@mcmc_samples$beta, 1, sd)
  h2D2@mcmc_quantile$beta = apply(h2D2@mcmc_samples$beta, 1, quantile, 
                                  probs = c(0.025,0.25,0.5,0.75,0.975))
  
  sigma2 = exp(h2D2@mcmc_samples$tau)
  h2D2@mcmc_mean$sigma2 = rowMeans(sigma2)
  h2D2@mcmc_sd$sigma2 = apply(sigma2, 1, sd)
  h2D2@mcmc_quantile$sigma2 = apply(sigma2, 1, quantile, 
                                    probs = c(0.025,0.25,0.5,0.75,0.975))
  
  for(l in c("h2", "h2_beta"))
  {
    h2D2@mcmc_mean[[l]] = mean(h2D2@mcmc_samples[[l]])
    h2D2@mcmc_sd[[l]] = sd(h2D2@mcmc_samples[[l]])
    h2D2@mcmc_quantile[[l]] = quantile(h2D2@mcmc_samples[[l]], 
                                       probs = c(0.025,0.25,0.5,0.75,0.975))
  }
  
  h2D2@CL = apply(h2D2@mcmc_samples$beta, 1, function(x){
    abs(sum(x > 0) - sum(x < 0)) / length(x)})
  if(get_CS)
  {
    h2D2@CS = h2D2_CS(h2D2, h2D2@coverage, h2D2@purity)
  }
  
  # Convergence
  n = h2D2@mcmc_samples$n_samples
  n_2 = floor(n/2)
  h2D2@mcmc_samples$PSRF_beta = foreach(j = 1:h2D2@M, .combine = "c") %do%
  {
    s1 = h2D2@mcmc_samples$beta[j,1:n_2]
    s2 = h2D2@mcmc_samples$beta[j,(n-n_2+1):n]
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