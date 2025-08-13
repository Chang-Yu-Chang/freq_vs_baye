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
- Compute bayesian R2, using var_fit / (var_fit + var_res), where var_fit is the variance among the expectations of the new data; var_res is the empirical residual variance

### Linear mixed model

`lmm.R`

- Compared `lme4::lmer` and `brm`
- Perform variance component analysis. Use `VarCorr` for `lmer` object and compute per-sample group variance over total variance for `brm`
- Compute bayesian R2, using var_fit / (var_fit + var_res), where var_fit is the variance among the expectations of the new data; var_res is the empirical residual variance

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
- Hand calculate bayesian R2, basically the empirical R2 = var_mu / (var_mu + var_res), whereas var_mu is the var of predicted response, var_res is the empirical residual variance
    - Bayesian R2 is estaimted at the response (data, or count) scale
- Estimate probability of zero inflation
- For GLM, it's tricky to calculate the R2 in glm because of the latent vs response scales


### Generalized linear mixed model

`glmm.R`

- Compare `glmmTMB` and `brm`
- Simulate 6 types of response distributions as in glm.R
- Random effect only includes one predictor (`group`)
- Hand calculate Pearson residuals
- Calculate R2_GLMM
    - Use `insight::get_variance` to obtain the variance of fixed, random effects, and residuals
    - Marginal R2 only includes fixed effects. var(fix) / (var(fix) + var(ran) + var(res))
    - Conditional R2 includes both fixed and random effects. (var(fix) + var(ran)) / (var(fix) + var(ran) + var(res))
- Calculate Bayesian R2 from empirical draws R2 = var_mu / (var_mu + var_res)

### GLMM with covariance matrix 

`glmm_cov.R`

- GLMM with a group factor and a covariance matrix as random effects
- Test 3 responses: gaussian, negative binomial, and ZINB
- Compute three R2:
    - GLMM R2 (nagakawa) on the link (linear) scale. It allows variance partitioning
    - Bayesian R2 on the response scale. It allows "how much var in the response is explained by the model"
    - semi-partial "part R2" on bayesian R2. 

### GLMM with nested random effects

`glmm_nes.R`

- GLMM with two random effects, and they are nested
- test brms and glmmTMB
- Test 3 responses: gaussian, negative binomial, and ZINB
- Variance components of random effects
    - Laten scale
    - Respsone scale
- Compute three R2:
    - GLMM R2 (nagakawa) on the link (linear) scale. It allows variance partitioning. This does not work for zero inflation
    - Bayesian R2 on the response scale. It allows "how much var in the response is explained by the model"
    - semi-partial "part R2" on bayesian R2.  -> no it's for fixed effects


### Phylogenetic linear mixed model 

`plmm.R`

- Only `brm`
- Random effect includes phylogenetic effect (`phylo_cov`)
- Phylogenetic effect is sampled from Multivariate Gaussian Distribution based on the phylo covariance matrix
- Compare one measurement per species and multiple measurements per species
- Sampling from a multivariate Gaussian based on a cov is equivalent to sampling from a gaussian multiplied by the Cholesky decomposition of the cov

- Well it seems that the estimate of phylo variance is only trustworthy when there are 50-100 species


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

