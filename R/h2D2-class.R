#' @title The h2D2 class
#' 
#' @description 
#' The h2D2 object is used to store GWAS summary data, LD matrix, variables
#' useful for MCMC sampling, h2-D2 hyper-parameters, posterior samples, and
#' credible sets. Objects can be created by Createh2D2Object.
#' 
#' @slot M Number of SNPs.
#' @slot N GWAS sample size. For binary traits, 'N' is the sum of the number 
#' of cases and the number of controls.
#' @slot z An M-vector of z-scores.
#' @slot R An M-by-M correlation matrix in \code{\link[Matrix]{dsCMatrix-class}} 
#' format. The diagonal elements are set as 0.
#' @slot SNP_ID IDs of SNPs. A character vector of length M.
#' @slot values Eigenvalues of R (with diagonal elements as 1). An M-vector
#' of nonnegative values.
#' @slot tvectors An M-by-M matrix. Each column of its transpose is an
#' eigenvector of R.
#' @slot lambda Hyper-parameter quantifying the discrepancy between z-scores
#' and R.
#' @slot a An M-vector of shape parameters for sigma^2.
#' @slot b Shape parameter for 1-h^2.
#' @slot mcmc_samples A list storing MCMC samples.
#' @slot mcmc_mean A list storing posterior sample means.
#' @slot mcmc_sd A list storing posterior sample standard errors.
#' @slot mcmc_quantile A list storing posterior sample 2.5\%, 25\%, 50\%, 75\%, 
#' and 97.5\% quantiles.
#' @slot CL An M-vector of credible levels of variables.
#' @slot coverage A number between 0 and 1 specifying the "coverage" of the 
#' credible sets.
#' @slot purity A number between 0 and 1 specifying the minimum absolute 
#' correlation allowed in a credible set.
#' @slot CS A list of identified credible sets. See \code{\link{h2D2_CS}}.
#' 
#' @export

setClass("h2D2", 
         slots = c(
           M = "numeric", #Number of SNPs
           N = "numeric", #Number of individuals
           z = "numeric", #z-scores
           R = "dsCMatrix", #Genotype correlation, symmetric sparse, triu
           SNP_ID = "character",
           values = "numeric",
           tvectors = "matrix",
           lambda = "numeric",
           a = "numeric", #Shape parameter for sigma2
           b = "ANY", #Shape parameter for 1-h2
           mcmc_samples = "list",
           mcmc_mean = "list",
           mcmc_sd = "list",
           mcmc_quantile = "list",
           CL = "numeric",
           coverage = "numeric",
           purity = "numeric",
           CS = "list"
         ))