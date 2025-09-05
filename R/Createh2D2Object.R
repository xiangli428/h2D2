#' @title Create an h2D2 object.
#' 
#' @description 
#' Create an h2D2 object from GWAS summary statistics and an LD matrix estimate.
#' GWAS z-scores, LD matrix, h2-D2 hyper-parameters, MCMC samples, and outputs 
#' are stored.
#' 
#' @usage 
#' h2D2 = Createh2D2Object(z,
#'                         R,
#'                         N,
#'                         SNP_ID = NULL,
#'                         in_sample_LD = FALSE,
#'                         lambda = NULL,
#'                         R_eig = NULL,
#'                         a = 0.005,
#'                         b = NULL,
#'                         coverage = 0.95,
#'                         purity = 0.5,
#'                         tol = 1e-8)
#'
#' @param z An M-vector of z-scores.
#' @param R An M-by-M LD matrix that can be coerced to 
#' \code{\link[Matrix]{dsCMatrix-class}}.
#' The diagonal elements of 'R' will be coerced to zeros.
#' @param N GWAS sample size. For binary traits, 'N' is the sum of the
#' number of cases and the number of controls.
#' @param SNP_ID Identifiers of SNPs. The default is c("SNP_1", ...).
#' @param in_sample_LD Logical. Setting in_sample_LD = TRUE if in-sample LD
#' matrix is used.
#' @param lambda Hyper-parameter quantifying the discrepancy between z-scores
#' and R. If in_sample_LD == TRUE, lambda will be set as 0. Otherwise, if
#' lambda is not provided, lambda will be estimated.
#' @param R_eig An output of "eigen" function. A list with two components
#' "values" and "vector". If in_sample_LD == TRUE, don't need to provide it.
#' If provided, it will be used to create h2D2 object directly. Otherwise, 
#' it will be computed.
#' @param a Shape parameters for sigma^2. Either a positive real number or an
#' M-vector of positive real numbers.
#' @param b Shape parameter for 1-h^2. A positive real number. If not specified,
#' it will be estimated by an empirical Bayesian approach.
#' @param coverage A number between 0 and 1 specifying the required coverage
#' of the credible sets.
#' @param purity A number between 0 and 1 specifying the minimum 
#' absolute correlation allowed in a credible set.
#' @param tol Eigenvalues less than tol will be coerced to zeros.
#' 
#' @return An h2D2 object. See \code{\link{h2D2-class}}.
#' The input z-scores will be modified to z / sqrt(1 + z^2/N - 1/N).
#' 
#' @export

Createh2D2Object <- function(z,
                             R,
                             N,
                             SNP_ID = NULL,
                             in_sample_LD = F,
                             lambda = NULL,
                             R_eig = NULL,
                             a = 0.005,
                             b = NULL,
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
  
  # Check sample size.
  if(length(N) > 0)
  {
    N = N[1]
  }
  if(!is.numeric(N) | N < 0)
  {
    stop("'N' must be a positive number.")
  }
  
  z = z / sqrt(1 + z^2 / N - 1 / N)
  
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
  
  R %<>% as.matrix()

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
  
  # Compute eigen decomposition of R.
  
  if(!in_sample_LD)
  {
    if(is.null(R_eig))
    {
      R_eig = eigen(R, symmetric = T)
      R_eig$values = R_eig$values + 1
    }
    values = R_eig$values
    values[values < tol] = 0
    tvectors = t(R_eig$vectors)
    
    if(is.null(lambda))
    {
      j = which.max(abs(z))
      dz = z - R[,j] * z[j]
      dz[j] = 0
      Udz2 = (tvectors %*% dz)^2
      
      gr = function(lambda)
      {
        denom = lambda + values
        sum(1 / denom - Udz2 / denom^2)
      }
      
      if(gr(1e-100) > 0)
      {
        opt = list(root = 0)
      } else {
        opt = uniroot(gr, interval = c(1e-100, 1e100), tol = 1e-10)
      }
      
      if(opt$root <= 1e-100)
      {
        lambda = max(mean(dz^2) - 1, 0)
      } else {
        lambda = max(mean(dz^2) - 1, opt$root)
      }
    }
  } else {
    values = numeric()
    tvectors = matrix(0,0,0)
    lambda = 0
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
  if(any(a < 0))
  {
    stop("Shape parameters 'a' should be positive.")
  }
  
  if(!is.null(b))
  {
    if(b <= 1)
    {
      warning("Setting b<=1 may lead to divergent result.")
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
  
  R = as(as(as(R, "dMatrix"), "symmetricMatrix"), "CsparseMatrix")
  
  return(new("h2D2",
             M = M,
             N = N,
             z = z,
             R = R,
             SNP_ID = SNP_ID,
             values = values,
             tvectors = tvectors,
             lambda = lambda,
             a = a,
             b = b,
             mcmc_samples = list("n_samples" = 0, "n_burnin" = 0),
             mcmc_mean = list(),
             mcmc_sd = list(),
             mcmc_quantile = list(),
             CL = rep(0,M),
             coverage = coverage,
             purity = purity,
             CS = list()))
}