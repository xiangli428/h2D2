# h2D2 v2.0

The `h2D2` package implements a fine-mapping method based on the 
"heritability-induced Dirichlet Decomposition" (h2-D2) prior.
It is a novel fine-mapping method that utilizes a continuous global-local 
shrinkage prior.
This method can be applied to both quantitative and binary traits.
An MCMC algorithm is employed to obtain samples from the posterior 
distribution.
This method can also provide "credible sets" of candidate
causal variants, which are generally provided by fine-mapping methods
based on discrete-mixture priors.
In version 2.0, we implement an approach for robust fine mapping with
out-of-sample LD matrix.

## Quick Start

### 1. Install `h2D2` from GitHub
```
devtools::install_github("https://github.com/xiangli428/h2D2")
library(h2D2)
```

### 2. load test dataset
```
data(h2D2_test_data)
```
When out-of-sample LD matrix is used, eigendecomposition of `R` will be 
performed during the creation of `h2D2` object. By default, set `R_eig = NULL`.
It can also be precomputed by `eigen`.
```
R_eig = eigen(R, symmetric = TRUE)

# If the diagonal elements of R are already set as 0, run
R_eig$values = R_eig$values + 1 
```

### 3. Create an h2D2 object with summary data

```
h2D2 = Createh2D2Object(z,
                        R,
                        N,
                        SNP_ID,
                        in_sample_LD = FALSE,
                        R_eig = R_eig,
                        a = 0.005)
```
where `a` specify the hyper-parameters of the prior.
In v2.0, by default, when both `in_sample = TRUE` and `in_sample = FALSE`, we 
recommend setting `b = NULL` and `b` will be estimated by a pre-training 
process before MCMC. 

### 4. MCMC sampling

```
h2D2 = h2D2_MCMC(h2D2, mcmc_n = 5500, burn_in = 500, thin = 2, stepsize = 2, 
                 seed = 428, get_CS = T, n_chain = 3, fold = 1.1,
                 pre_mcmc_n = 400, pre_burn_in = 200, pre_p = 0.05,
                 pre_maxiter = 10, pre_miniter = 2)
```
In v2.0, we introduce parallel tempering to improve mixing of MCMC. Set
`n_chain = 1` can turn off parallel tempering.

### 5. Results

Credible levels of SNPs:
```
h2D2@CL
```
Credible sets:
```
h2D2@CS
```

### 6. Citation

Li, X., Sham, P. C., & Zhang, Y. D. (2024). A Bayesian fine-mapping model using 
a continuous global-local shrinkage prior with applications in prostate cancer 
analysis. The American Journal of Human Genetics.
https://doi.org/10.1016/j.ajhg.2023.12.007
