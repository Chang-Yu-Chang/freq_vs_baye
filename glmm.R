#' Compare glmm between glmmTMB and brms

library(tidyverse)
library(easystats)
library(glmmTMB)
library(brms)

# Simulate data with a poisson distribution ----
set.seed(123)
n_groups <- 10
n_per_group <- 100
group <- factor(rep(1:n_groups, each = n_per_group))
X <- rep(rnorm(n_per_group, 0, 1), n_groups)

# Random effect, group effect
b_count <- rnorm(n_groups, 0, 1)
re_count <- b_count[group]

# Linear predictor
# for poisson: log(lambda) = eta = beta * X. lambda = exp(eta)
# for nbinom: log(mu) = eta = beta * X. mu = exp(eta)
eta <- 0.01 + 0.08 * X + re_count
Y <- rpois(n_groups*n_per_group, lambda = exp(eta))
tb <- tibble(Y, X, eta, group)
hist(Y, breaks = 100)

# Frequentist ----
mod <- glmmTMB(Y ~ X + (1|group), family = poisson, data = tb)
summary(mod)

# Residual at the response scale
var(resid(mod, type = "response"))
var(Y - predict(mod, type = "response"))

## Pearson residual = (y - y_hat) / sqrt(var(y_hat))
var(resid(mod, type = "pearson"))
lambda_hat <- predict(mod, type = "response")
var((Y - lambda_hat) / sqrt(lambda_hat)) # poisson has same mean and variance

# Random effects
summary(mod)
VarCorr(mod)[[1]][[1]][1] # Model based theoretical estimate of the entire random effect distribution
as_tibble(ranef(mod))[["condval"]] %>% var # Empirical BLUPs (Best Linear Unbiased Predictors) for each group


# Variance on the response scale
## Total variance minus those explained by residuals
var_res <- var(Y - predict(mod, type = "response")) # residual variance. it's build in poisson
var_tot <- var(Y) # total variance
1 - var_res/var_tot

## Variance explined by random effect
var_ran <- VarCorr(mod)[[1]][[1]][1] # Model based theoretical estimate of the entire random effect distribution
var_ran / var_tot



# Bayesian ----
mod2 <- brm(Y ~ X + (1|group), data = tb, family = poisson, chains = 2, cores = 2, iter = 10000, thin = 10, seed = 123)

# Fixed effect estimates
fixef(mod2)

# Posterior predictive checks
pp_check(mod2)

# Compute Bayesian R2
bayes_R2(mod2)

# Variance components
set.seed(1)
pp <- posterior_predict(mod2)
pp_mean <- colMeans(pp) # mean posterior predictor
var_res2 <- var(Y - pp_mean)

## Proportion of total variance explained not by residual
1 - var_res2 / var(Y) # should be similar to bayes R2

## variance explained by random effects
ran_samples <- as_draws_df(mod2) %>%
    as_tibble() %>%
    select(starts_with("r_group"))

apply(ran_samples, 2, mean) # group effect averaged over posterior samples
ranef(mod2)$group[,,1]

var_ran <- apply(ran_samples, 1, var) # posterior distribution of variance of random effects
(var_ran / var_tot) %>% hist
mean(var_ran / var_tot)
