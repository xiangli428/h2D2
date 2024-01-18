#' @title Estimate local heritability using high-definition likelihood method.
#' 
#' @description 
#' Obtain maximum likelihood estimate of local heritability using z-scores and
#' LD matrix.
#' 
#' @usage 
#' HDL(z, N, R_eig = R_eig)
#' HDL(z, N, R = R)
#'
#' @param z M-vector of z-scores.
#' @param N GWAS sample size. For binary traits, 'N' is the sum of the
#' number of cases and the number of controls.
#' @param R_eig Spectral decomposition of LD matrix. A list return by 
#' \code{\link[eigen]{base}} with components 'values' and 'vectors'. The
#' diagonal elements of LD matrix are zeros.
#' @param R An M-by-M LD matrix with all diagonal elements zeros.
#' The diagonal elements will be coerced to zeros.
#' @param tol A numeric tolerance. Eigenvalues smaller than this value will be
#' removed.
#' 
#' @return A numeric value of local heritability estimate.
#' 
#' @references https://www.nature.com/articles/s41588-020-0653-y
#' 
#' @export

HDL = function(z, N, R = NULL, R_eig = NULL, tol = 1e-8)
{
  if(is.null(R_eig))
  {
    if(is.null(R))
    {
      stop("Either 'R' or 'R_eig' should be input.")
    } else {
      R_eig = eigen(R, symmetric = T)
    }
  }
  
  d = R_eig$values + 1
  u2 = ((z %*% R_eig$vectors)[1,])^2
  
  idx = which(d > tol)
  d = d[idx]
  d2 = d^2
  u2 = u2[idx]
  
  P = length(idx)
  
  fn = function(h2)
  {
    denom = d + N * h2 * d2 / P
    sum(log(denom)) + sum(u2 / denom)
  }
  
  opt = optimize(fn, c(0,1))
  opt$minimum
}
