# freq_vs_baye

A learing note for side-by-side comparison of frequentist and bayesian models

```
library(brms)
library(easystats)
```

### Ordinary linear regression 

`lm.R`

- Compare `lm` and `brm`
- Hand calculated R^2 and F statistics from `lm` object


### Linear mixed model

`lmm.R`

- Compared `lme4::lmer` and `brm`
- Perform variance component analysis. Use `VarCorr` for `lmer` object and compute per-sample group variance over total variance for `brm`

### Generalized linear model

`glm.R`

- Compare `glmmTMB` and `brm`
- Simulate 6 types of response distributions
    - gaussian
    - logistic
    - poisson
    - negative binomial
    - zero-inflated poisson
    - zero-inflated negative binomial
- Estimate R2 and Bayesian R2
- Estimate probablity of zero inflation


- For GLM, the variance is usually calculated at the latent scale. It's tricky to calculate the R2 in glm because of the latent vs response scales



### Generalized linear mixed model

`glmm.R`

- Compare `glmmTMB` and `brm`
- Simulate a poisson  binomial
- Random effect only includes one predictor (`group`)
- Hand calculate Pearson residual
- Bayesian R2
- Estimate random effect of each group by the mean across posterior samples
- Estimate the posterior distribution of proportion of variance explained by random effects

### Phylogenetic linear mixed model 

`plmm.R`

- Only `brm`
- Random effect includes phylogenetic effect (`phylo_cov`)
- Phylogenetic effect is sampled from Multivariate Gaussian Distribution based on the phylo covariance matrix
- Compare one measurement per species and multiple measurements per species
- Sampling from a multivariate Gaussian based on a cov is equivalent to sampling from a gaussian multiplied by the Cholesky decomposition of the cov

- Well it seems that the estimate of phylo variance is only trustworthy when there are 50-100 species

### GLMM with covariance matrix 

`glmm_cov.R`

- Testing simple covariance matrix


### Phylogenetic generalized linear mixed model

`pglmm.R`

- Only `brm`
- Random effect includes phylogenetic effect (`phylo_cov`) and species effect that is not captured by shared phylogeny (`species`)
- Phylogenetic effect is sampled from Multivariate Normal Distribution based on the phylo covariance matrix
- Generate counts by negative binomial
- Multiple measurements per species
- Compare two models: species mean or raw species value

- If the number of species is too small, the intrinsic variation in NB can mask phylogenetic and species effect.

### Phylogenetic generalized linear mixed model with different distribution

`pglmm2.R`

- Generate counts by zero-inflated negative binomial

