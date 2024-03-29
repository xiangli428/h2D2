#' @title Create an h2D2 object.
#' 
#' @description 
#' Create an h2D2 object from GWAS summary statistics and an LD matrix estimate.
#' GWAS z-scores, h2-D2 hyper-parameters, and useful variables are stored.
#' 
#' @usage 
#' h2D2 = Createh2D2Object(z,
#'                         R,
#'                         N,
#'                         SNP_ID = NULL,
#'                         trait = c("quantitative","binary"),
#'                         in_sample_LD = F,
#'                         a = 0.005,
#'                         b = 1e4,
#'                         coverage = 0.95,
#'                         purity = 0.5)
#'
#' @param z M-vector of z-scores.
#' @param R An M-by-M LD matrix that can be coerced to 
#' \code{\link[Matrix]{dsCMatrix-class}}.
#' The diagonal elements will be coerced to zeros.
#' @param N GWAS sample size. For binary traits, 'N' is the sum of the
#' number of cases and the number of controls.
#' @param SNP_ID Identifiers of SNPs. The default is c("SNP_1", ...).
#' @param trait Either "quantitative" or "binary".
#' @param in_sample_LD Whether the LD matrix is in-sample.
#' @param a Shape parameters for sigma^2. Either a positive real number or an
#' M-vector of positive real numbers.
#' @param b Shape parameter for 1-h^2. Must be a positive real number.
#' If the in-sample LD matrix is provided, we recommend setting b = NULL and
#' b will be estimated by a pre-training process before MCMC.
#' @param coverage A number between 0 and 1 specifying the "coverage"
#' of the credible sets.
#' @param purity A number between 0 and 1 specifying the minimum 
#' absolute correlation allowed in a credible set.
#' 
#' @return An h2D2 object. See \code{\link{h2D2-class}}.
#' For quantitative traits, the input z-scores will be modified to
#' z / sqrt(1 + z^2/N - 1/N).
#' 
#' @export

Createh2D2Object <- function(z,
                             R,
                             N,
                             SNP_ID = NULL,
                             trait = "quantitative",
                             in_sample_LD = F,
                             a = 0.005,
                             b = 1e4,
                             coverage = 0.95,
                             purity = 0.5,
                             tol = 1e-8)
{
  # Check z.
  
  M = length(z)
  if(M <= 1)
  {
    stop("Should contain at least 2 SNPs.")
  }
  if(any(is.na(z)))
  {
    stop("'z' cannot contain missing values.")
  }
  
  # Check SNP ID.
  
  if(is.null(SNP_ID))
  {
    SNP_ID = paste("SNP", 1:M, sep = "_")
  }
  if(length(SNP_ID) != M)
  {
    warning("The length of 'SNP_ID' is not equal to the length of 'z'.")
    warning("Rename SNPs.")
    SNP_ID = paste("SNP", 1:M, sep = "_")
  }
  
  # Check R.
  
  R = as.matrix(R)
  
  if(nrow(R) != M | ncol(R) != M)
  {
    stop("Dimensions of LD matrix are not equal to the length of effect sizes.")
  }
  
  rownames(R) = colnames(R) = SNP_ID
  
  if(any(abs(R) > 1))
  {
    warning("There cannot be a value greater than 1 in the LD matrix.")
    warning("Values larger than 1 or smaller than -1 are coerced to be 1 or -1.")
    R[R > 1] = 1
    R[R < -1] = -1
  }
  if(!isSymmetric(R, check.attributes = F))
  {
    warning("LD matrix must be symmetric. Coerce to be symmetric.")
    R = (R + t(R)) / 2
  }
  if(any(diag(R) != 0))
  {
    diag(R) = 0
    warning("The diagonal elements of LD matrix are coerced to zeros.")
  }
  
  # Check trait.
  
  if(trait == "quantitative")
  {
    z = z / sqrt(1 + z^2 / N - 1 / N)
  } else if(trait != "binary") {
    stop("'trait' must be either 'quantitative' or 'binary'.")
  }
  
  # Projection
  
  if(!in_sample_LD)
  {
    R_eig = eigen(R, symmetric = T)
    if(sum(R_eig$values + 1 < tol) > 0)
    {
      idx = which(R_eig$values + 1 >= tol)
      z = (R_eig$vectors[,idx] %*% (t(R_eig$vectors[,idx]) %*% z))[,1]
    }
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
    stop("Shape parameters 'a' should be positive.")
  }
  
  if(!is.null(b))
  {
    if(b <= 1)
    {
      warning("Setting b<=1 may lead to a divergent result.")
    }
  }
  
  # Check coverage, purity, and rho
  
  if(coverage < 0 | coverage > 1)
  {
    stop("Coverage must be in a range between 0 and 1.")
  }
  if(purity < 0 | purity > 1)
  {
    stop("Purity must be in a range between 0 and 1.")
  }
  
  LD_pairs = R^2
  LD_pairs[LD_pairs < purity^2] = 0
  LD_pairs = t(LD_pairs / (rowSums(LD_pairs) + 5e-324))
  LD_pairs = as(LD_pairs, "dgCMatrix")
  
  R = as(R, "dsCMatrix")
  
  return(new("h2D2",
             M = M,
             N = N,
             SNP_ID = SNP_ID,
             z = z,
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
             CS = list()))
}