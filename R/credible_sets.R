#' @title Identify credible sets for h2D2 object.
#' 
#' @description 
#' A greedy search algorithm to identify all credible sets achieving target
#' "coverage" with minimum absolute correlation not less than "purity".
#' 
#' @param h2D2 An h2D2 object with MCMC samples.
#' @param coverage A number between 0 and 1 specifying the "coverage"
#' of the credible sets.
#' @param purity A number between 0 and 1 specifying the minimum 
#' absolute correlation allowed in a credible set.
#' 
#' @return A list with the following elements:
#' \describe{
#'   \item{sets}{A list in which each element is a vector containing the
#'   indices of the variables in the CS.}
#'   \item{CL}{The credible level of each CS.}
#'   \item{purity}{The purity of each CS.}
#' }
#' 
#' @export

credible_sets = function(h2D2, coverage = 0.95, purity = 0.5)
{
  M = h2D2@M
  n = h2D2@mcmc_samples$n_samples
  LD_pairs = as(h2D2@LD_pairs, "dgTMatrix")
  
  CL = data.frame("index" = 1:M,
                  "CL" = h2D2@CL)
  CL = arrange(CL, -CL, index)
  
  CS = list("sets" = list(), "CL" = c())
  n1 = sum(CL$CL > coverage)
  if(n1 > 0)
  {
    for(r in 1:n1)
    {
      CS$sets[[r]] = CL$index[r]
      CS$CL[r] = CL$CL[r]
    }
  }
  r = length(CS$sets) + 1
  
  candidate = CL$index[(n1+1):M]
  discard = c()
  CL = arrange(CL, index)
  
  while(length(candidate) > 0)
  {
    j = candidate[1]
    tmp_CL = CL$CL[j]
    tests = intersect(LD_pairs@i[LD_pairs@j == j-1]+1, 
                      c(candidate[-1], discard))
    set = c(j)
    findCS = F
    
    while((length(tests) > 0) & (tmp_CL + sum(CL$CL[tests]) > coverage))
    {
      tests_CL = foreach(k = tests, .combine = "c") %do%
      {
        idx = c(set, k)
        eig = eigen(h2D2@R[idx,idx])
        v = h2D2@mcmc_samples$beta[,idx] %*% 
          eig$vectors[,which.max(abs(eig$values + 1))]
        abs(sum(v > 0) - sum(v < 0)) / n
      }
      k = which.max(tests_CL)
      if(max(tests_CL) > tmp_CL)
      {
        set = c(set, tests[k])
        tmp_CL = max(tests_CL)
        if(tmp_CL > coverage)
        {
          CS$sets[[r]] = sort(set)
          CS$CL[r] = tmp_CL
          r = r + 1
          candidate = setdiff(candidate, set)
          discard = setdiff(discard, set)
          findCS = T
          break
        } else {
          tests = intersect(tests, LD_pairs@i[LD_pairs@j == tests[k]-1]+1)
        }
      } else {
        break
      }
    }
    
    if(!findCS)
    {
      candidate = candidate[-1]
      discard = c(discard,j)
    }
  }
  
  CS$purity = foreach(set = CS$sets, .combine = "rbind") %do%
  {
    if(length(set) > 1)
    {
      R = abs(h2D2@R[set, set]@x)
      data.frame("min.abs.corr" = min(R),
                 "mean.abs.corr" = mean(R),
                 "median.abs.corr" = median(R))
    } else {
      data.frame("min.abs.corr" = 1,
                 "mean.abs.corr" = 1,
                 "median.abs.corr" = 1)
    }
  }
  
  return(CS)
}