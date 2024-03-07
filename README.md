# h2D2 v1.1

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

## Quick Start

### 1. Install `h2D2` from GitHub
```
devtools::install_github("https://github.com/xiangli428/h2D2")
library(h2D2)
```
### 2. Create an h2D2 object with summary data

```
h2D2 = Createh2D2Object(z,
                        R,
                        N,
                        SNP_ID = NULL,
                        trait = "quantitative",
                        in_sample_LD = F,
                        a = 0.005,
                        b = 1e4,
                        coverage = 0.95,
                        purity = 0.5)
```
where `a` and `b` specify the hyperparameters of the prior.
In v1.1, if `in_sample_LD = T`, we recommend setting `b = NULL` and `b` will be estimated by a pre-training process before MCMC.

### 3. MCMC sampling

```
h2D2 = h2D2_MCMC(h2D2, mcmc_n = 10000, burn_in = 5000, 
                 thin = 1, stepsize = 2, seed = 428)
```

### 4. Results

Credible levels of SNPs:
```
h2D2@CL
```
95% Credible sets:
```
h2D2@CS
```

### 5. Citation

Li, X., Sham, P. C., & Zhang, Y. D. (2024). A Bayesian fine-mapping model using a continuous global-local shrinkage prior with applications in prostate cancer analysis. The American Journal of Human Genetics.
https://doi.org/10.1016/j.ajhg.2023.12.007
