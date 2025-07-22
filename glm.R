# GLM

library(tidyverse)
library(easystats)
library(glmmTMB)
library(brms)

# Simulate data with a poisson distribution ----
set.seed(123)
n <- 200
X <- rnorm(n, 0, 1)

# Linear predictor
# for poisson: log(lambda) = eta = beta * X. lambda = exp(eta)
# for nbinom: log(mu) = eta = beta * X. mu = exp(eta)
eta <- 1 + 0.8 * X
Y <- rpois(n, lambda = exp(eta))
# dispersion parameter k = 2
#Y <- rnbinom(n_groups*n_per_group, mu = exp(eta), size = 2)
tb <- tibble(Y, X, eta)
hist(Y, breaks = 100)


## Frequentist ----
mod <- glmmTMB(Y ~ X, data = tb, family = poisson)
summary(mod)

# Prediction. three ways
exp(predict(mod))
fitted(mod) # predicted at the response scale
predict(mod, type = "response")

# Residuals
resid(mod)
Y - exp(predict(mod))
check_model(mod)

# pseudo R2
## Based on variance of residuals that reduce the explanation of total variance
var_res <- var(Y - predict(mod, type = "response"))
var_tot <- var(Y)
1 - var_res/var_tot


## Bayesian ----
mod2 <- brm(Y ~ X, data = tb, family = poisson, chains = 2, cores = 2, iter = 10000, thin = 10, seed = 123)

# Fixed effect estimates
fixef(mod2)

# Posterior predictive checks
pp_check(mod2)

# Compute Bayesian R2
bayes_R2(mod2)

# Variance components
set.seed(1)
fit <- predict(mod2, summary = F) # predicted value
apply(fit, 1, var) %>% hist # estimated variance in response
var(Y) # observed variance in response

# Residual variance per sample
var_res_samples <- (matrix(rep(Y, 1000), ncol = length(Y), byrow = T) - fit) %>% apply(1, var)
range(var_res_samples)

# proportion of total variance explained
(1 - var_res_samples / var(Y)) %>% hist




# Simulate data with a negative binomial  distribution ----
set.seed(123)
n <- 200
X <- rnorm(n, 0, 3)

# Linear predictor
# for poisson: log(lambda) = eta = beta * X. lambda = exp(eta)
# for nbinom: log(mu) = eta = beta * X. mu = exp(eta)
# dispersion parameter k = 2
eta <- 1 + 0.8 * X
Y <- rnbinom(n, mu = exp(eta), size = 10) # size is the dispersion parameter
tb <- tibble(Y, X, eta)
hist(Y, breaks = 100)


## Frequentist ----
mod <- glmmTMB(Y ~ X, data = tb, family = nbinom2)
summary(mod)

# Prediction. three ways
exp(predict(mod))
fitted(mod) # predicted at the response scale
predict(mod, type = "response")

# Residuals
resid(mod)
Y - exp(predict(mod))

# pseudo R2
## Based on variance of residuals that reduce the explanation of total variance
var_res <- var(Y - predict(mod, type = "response"))
var_tot <- var(Y)
1 - var_res/var_tot


## Bayesian ----
mod2 <- brm(Y ~ X, data = tb, family = "negbinomial2", chains = 2, cores = 2, iter = 10000, thin = 10, seed = 123)

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

# proportion of total variance explained
(1 - var_res2 / var(Y))

# Compute R2 for each posterior sample
## By loop
R2 <- list()
for (i in 1:1000) R2[i] <- 1 - var(Y - pp[i,]) / var(Y)
unlist(R2)

## matrix
x <- pp - matrix(Y, nrow = nrow(pp), ncol = ncol(pp), byrow = T)
(1 - (apply(x, 1, var) / var(Y)))

