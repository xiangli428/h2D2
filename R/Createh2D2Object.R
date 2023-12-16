#' @title Create an h2D2 object.
#' 
#' @description 
#' Create an h2D2 object from GWAS summary statistics and an LD matrix estimate.
#' GWAS summary data, h2-D2 hyper-parameters, and useful variables are stored.
#' 
#' @usage 
#' For quantitative traits,
#' h2D2 = Createh2D2Object(z,
#'                         R,
#'                         N,
#'                         SNP_ID = NULL,
#'                         trait = "quantitative",
#'                         in_sample_LD = F,
#'                         a = 0.005,
#'                         coverage = 0.95,
#'                         purity = 0.5)
#' 
#' For binary traits,
#' h2D2 = Createh2D2Object(z,
#'                         R,
#'                         N1,
#'                         N0,
#'                         SNP_ID = NULL,
#'                         trait = "binary",
#'                         in_sample_LD = F,
#'                         a = 0.005,
#'                         coverage = 0.95,
#'                         purity = 0.5)
#'
#' @param z M-vector of z-scores.
#' @param N The sample size. Required for quantitative traits.
#' @param N1 Number of cases. Required for binary traits.
#' @param N0 Number of controls. Required for binary traits.
#' @param R An M-by-M LD matrix that can be coerced to 
#' \code{\link[Matrix]{dsCMatrix-class}}.
#' The diagonal elements should be set as 0s.
#' @param SNP_ID Identifiers of SNPs. The default is c("SNP_1", ...).
#' @param trait One of "quantitative" or "binary".
#' @param in_sample_LD Whether the LD matrix is in-sample.
#' @param a Shape parameters for sigma^2. Either a positive real number or an
#' M-vector of positive real numbers.
#' @param b Shape parameter for 1-h^2. A positive real number.
#' @param coverage A number between 0 and 1 specifying the "coverage"
#' of the credible sets.
#' @param purity A number between 0 and 1 specifying the minimum 
#' absolute correlation allowed in a credible set.
#' 
#' @return An h2D2 object. See \code{\link{h2D2-class}}.
#' 
#' @export

Createh2D2Object <- function(z,
                             R,
                             N = NULL,
                             N1 = NULL,
                             N0 = NULL,
                             SNP_ID = NULL,
                             trait = "quantitative",
                             in_sample_LD = F,
                             a = 0.005,
                             b = NULL,
                             coverage = 0.95,
                             purity = 0.5)
{
  # Check z.
  M = length(z)
  if(M <= 1)
  {
    stop("Should contain at least 2 SNPs.")
  }
  if(any(is.na(z)))
  {
    stop("z cannot contain missing values.")
  }
  
  # Check R.
  R = as.matrix(R)
  if(nrow(R) != M | ncol(R) != M)
  {
    stop("Dimensions of LD matrix are not equal to the length of effect sizes.")
  }
  if(any(abs(R) > 1))
  {
    stop("There cannot be a value greater than 1 in the LD matrix.")
  }
  if(!isSymmetric(R))
  {
    stop("LD matrix must be symmetric.")
  }
  if(any(diag(R) != 0))
  {
    diag(R) = 0
    warning("The diagonal elements of LD matrix are coerced to zeros.")
  }
  R_eig = eigen(R, symmetric = T)
  if(R_eig$values[M] < -1)
  {
    warning("LD matrix is not positive semidefinite.")
  }
  
  # Check trait.
  if(trait == "quantitative")
  {
    if(is.null(N))
    {
      stop("N is required for quantitative traits.")
    }
    betaHat = z / sqrt(N + z^2)
  } else if(trait == "binary")
  {
    if(is.null(N1) | is.null(N0))
    {
      stop("N1 and N0 are required for binary traits.")
    }
    N = N1 + N0
    betaHat = z * sqrt(N/N1/N0)
  } else {
    stop("'trait' must be either 'quantitative' or 'binary'.")
  }
  
  # Check SNP ID.
  if(is.null(SNP_ID))
  {
    SNP_ID = paste("SNP", 1:M, sep = "_")
  }
  if(length(SNP_ID) != M)
  {
    warning("The length of SNP_ID is not equal to the length of effect sizes. Rename SNPs.")
    SNP_ID = paste("SNP", 1:M, sep = "_")
  }
  
  # Hyper-parameters
  if(length(a) == 1)
  {
    a = rep(a, M)
  }
  if(length(a) != M)
  {
    stop("The length of 'a' is not equal to the length of effect sizes.")
  }
  if(any(a <= 0))
  {
    stop("Shape parameters a should be positive.")
  }
  
  if(is.null(b))
  {
    if(in_sample_LD)
    {
      # compute HESS estimater
      order_SNP = order(abs(z), decreasing = T)
      keep_SNPs = order_SNP[1]
      for(j in order_SNP[2:M])
      {
        if(all(abs(R[j,keep_SNPs]) < 0.5))
        {
          keep_SNPs = c(keep_SNPs, j)
        }
      }
      keep_SNPs = sort(keep_SNPs)
      p = length(keep_SNPs)
      bVb = crossprod(solve(R[keep_SNPs,keep_SNPs] + diag(p),
                            betaHat[keep_SNPs], system = "LDLt"),
                      betaHat[keep_SNPs])[1]
      if(trait == "quantitative")
      {
        h2_hess = (N * bVb - p) / (N - p)
      } else {
        h2_hess = N1 * N0 / N^2 * bVb - p / N
      }
      b = sum(a) * (1 - h2_hess) / h2_hess
      
      if(b <= 1)
      {
        warning("Setting b as default value 1e5.")
      }
    } else {
      b = 1e5
    }
  } else {
    if(b < 1)
    {
      warning("Setting b<=1 may lead to divergent result.")
    }
  }
  
  R = as(R, "dsCMatrix")
  
  # Check coverage, purity, and rho
  if(coverage < 0 | coverage > 1)
  {
    stop("Coverage must be in a range between 0 and 1.")
  }
  if(purity < 0 | purity > 1)
  {
    stop("Purity must be in a range between 0 and 1.")
  }
  
  LD_pairs = as.matrix(R^2)
  LD_pairs[LD_pairs < purity^2] = 0
  LD_pairs = t(LD_pairs / (rowSums(LD_pairs) + 5e-324))
  LD_pairs = as(LD_pairs, "dgCMatrix")
  
  return(new("h2D2",
             M = M,
             N = N,
             N1 = N1,
             N0 = N0,
             SNP_ID = SNP_ID,
             z = z,
             betaHat = betaHat,
             R = R,
             trait = trait,
             LD_pairs = LD_pairs,
             a = a,
             b = b,
             mcmc_samples = list(),
             mcmc_mean = list(),
             mcmc_sd = list(),
             mcmc_quantile = list(),
             CL = rep(0,M),
             coverage = coverage,
             purity = purity,
             CS = list()
  ))
}