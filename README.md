# h2D2

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

For quantitative traits,
```
h2D2 = Createh2D2Object(betaHat,
                        sigmaHat,
                        R,
                        SNP_ID = NULL,
                        trait = "quantitative",
                        N,
                        a = 0.005,
                        b = 2e5,
                        coverage = 0.95,
                        purity = 0.5)
```
For binary traits,
```
h2D2 = Createh2D2Object(betaHat,
                        sigmaHat,
                        R,
                        SNP_ID = NULL,
                        trait = "binary",
                        N1,
                        N0,
                        a = 0.005,
                        b = 2e5,
                        coverage = 0.95,
                        purity = 0.5)
```
where `a` and `b` specify the hyperparameters of the prior.

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
