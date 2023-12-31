% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/h2D2-package.R
\docType{package}
\name{h2D2}
\alias{h2D2}
\alias{h2D2-package}
\title{h2D2 package}
\description{
h2-D2: heritability-induced Dirichlet decomposition
}
\details{
This package implements a fine-mapping method based on the 
"heritability-induced Dirichlet Decomposition" (h2-D2) prior.
It is a novel fine-mapping method that utilizes a continuous global-local 
shrinkage prior.
This method can be applied to both quantitative and binary traits.
An MCMC algorithm is employed to obtain samples from the posterior 
distribution.
This method can also provide "credible sets" of candidate causal variants, 
which are generally provided by fine-mapping methods based on discrete-mixture priors.
}
\references{
https://doi.org/10.1016/j.ajhg.2023.12.007
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/xiangli428/h2D2}
}

}
\author{
\strong{Maintainer}: Xiang Li \email{freddyl@connect.hku.hk}

Authors:
\itemize{
  \item Yan Dora Zhang \email{doraz@hku.hk}
}

}
\keyword{internal}
