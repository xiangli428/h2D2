#' @title h2D2 package
#' @name h2D2
#' @description h2-D2: heritability-induced Dirichlet decomposition
#' 
#' @import Rcpp
#' @import RcppEigen
#' @import Matrix
#' @import foreach
#' @import dplyr
#' 
#' @details 
#' This package implements a fine-mapping method based on the 
#' "heritability-induced Dirichlet Decomposition" (h2-D2) prior.
#' It is a novel fine-mapping method that utilizes a continuous global-local 
#' shrinkage prior.
#' This method can be applied to both quantitative and binary traits.
#' An MCMC algorithm is employed to obtain samples from the posterior 
#' distribution.
#' This method can also provide "credible sets" of candidate causal variants, 
#' which are generally provided by fine-mapping methods based on discrete-mixture priors.
#' 
#' @keywords internal
#' 
#' @references https://www.medrxiv.org/content/10.1101/2023.08.04.23293456v1

"_PACKAGE"

