#' @title Create an h2D2 object.
#' 
#' @description 
#' Create an h2D2 object from GWAS summary statistics and an LD matrix estimate.
#' GWAS summary data, h2-D2 hyper-parameters, and useful variables are stored.
#'
#' @param betaHat M-vector of standardized effect sizes. Per-allele effect sizes
#' should be multiplied by sqrt(2\*MAF\*(1-MAF)) before input.
#' @param sigmaHat M-vector of standard errors of standardized effect sizes.
#' Standard errors of per-allele effect sizes should be multiplied by 
#' sqrt(2\*MAF\*(1-MAF)) before input.
#' @param R An M-by-M LD matrix that can be coerced to "dsCMatrix".
#' The diagonal elements should be set as 0s.
#' @param SNP_ID Identifiers of SNPs. The default is c("SNP_1", ...).
#' @param trait One of "quantitative" or "binary".
#' @param N The sample size. Required for quantitative traits.
#' @param N1 Number of cases. Required for binary traits.
#' @param N0 Number of controls. Required for binary traits.
#' @param a Shape parameters for sigma^2. Either a positive real number or an
#' M-vector of positive real numbers.
#' @param b Shape parameters for 1-h^2. A positive real number.
#' @param coverage A number between 0 and 1 specifying the "coverage"
#' of the credible sets.
#' @param purity A number between 0 and 1 specifying the minimum 
#' absolute correlation allowed in a credible set.
#' 
#' @export

Createh2D2Object <- function(betaHat,
                             sigmaHat,
                             R,
                             SNP_ID = NULL,
                             trait = "quantitative",
                             N = NULL,
                             N1 = NULL,
                             N0 = NULL,
                             a = 0.005,
                             b = 2e5,
                             coverage = 0.95,
                             purity = 0.5)
{
  # Check betaHat and sigmaHat.
  M = length(betaHat)
  if(M <= 1)
  {
    stop("Should contain at least 2 SNPs.")
  }
  if(length(sigmaHat) != M)
  {
    stop("The length of standard errors is not equal to the length of effect sizes.")
  }
  if(any(is.na(betaHat)) | any(is.na(sigmaHat)))
  {
    stop("betaHat and sigmaHat cannot contain missing values.")
  }
  if(any(sigmaHat <= 0))
  {
    stop("sigmaHat should be positive.")
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
  if(R_eig$values[M] <= -1)
  {
    warning("LD matrix is nearly singular.")
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
  
  if(length(a) == 1)
  {
    a = rep(a, M)
  }
  if(length(a) != M)
  {
    stop("The length of 'a' is not equal to the length of effect sizes.")
  }
  
  if(any(a <= 0) | b <= 0)
  {
    stop("Shape parameter should be positive.")
  }
  
  if(b <= 1)
  {
    warning("Setting b<=1 may lead to divergent result.")
  }
  
  # Check trait.
  if(trait == "quantitative")
  {
    S = sqrt(sigmaHat^2 + betaHat^2 / N)
  } else if(trait == "binary")
  {
    S = sigmaHat
  } else {
    stop("'trait' must be either 'quantitative' or 'binary'.")
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
             betaHat = betaHat,
             S = S,
             z = betaHat / sigmaHat,
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