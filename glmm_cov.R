#' Simulate covariance matrix as the random effect in glmm

library(tidyverse)
library(easystats)
library(brms)

# Simulation ----
# Example data with grouping factor
set.seed(123)
n_groups <- 10
n_per_group <- 10
group_names <- paste0("g", 1:n_groups)
group <- factor(rep(group_names, each = n_per_group))
X <- rep(1:n_per_group, times = n_groups)

# Random intercepts for cov effects
#cov <- matrix(c(1, .8, .8, 1), nrow = n_groups)
A <- matrix(rnorm(n_groups*n_groups), nrow=n_groups)
cov <- t(A) %*% A # Create a covariance matrix as A * A'
colnames(cov) <- group_names
rownames(cov) <- group_names
cov_scaled <- cov/sum(diag(cov)) # rescaled the cov to have a total variance = 1
#cov_effects <- MASS::mvrnorm(n = 1, mu = rep(0, nrow(cov)), Sigma = cov_scaled)
cov_effects <- MASS::mvrnorm(n = 1, mu = rep(0, nrow(cov)), Sigma = cov)
names(cov_effects) <- group_names

#
Y <- 1 + 1.6 * X + cov_effects[group] + rnorm(n_groups * n_per_group, 0, 1)
tb <- tibble(cov_group = group, X = X, Y = Y)

#
mod <- brm(
    Y ~ X + (1|gr(cov_group, cov = cov)),
    data = tb,
    family = gaussian(),
    data2 = list(cov = cov_scaled),
    chains = 2, cores = 2, iter = 10000, thin = 10, seed = 123,
    control = list(adapt_delta = 0.95)
)

pp_check(mod)

hyp <- hypothesis(mod, class = NULL, "sd_cov_group__Intercept^2 / (sd_cov_group__Intercept^2 + sigma^2) = 0")
plot(hyp)
var(cov_effects)



all_samples <- MASS::mvrnorm(n=1000, mu=rep(0, 2), Sigma=sigma)
cov_estimate <- cov(all_samples)
sum(diag(cov_estimate))
