#' Compare lmm between lmer and brms

library(tidyverse)
library(easystats)
library(brms)
library(lme4)

# Example data with grouping factor
set.seed(123)
n_groups <- 10
n_per_group <- 100
group <- factor(rep(1:n_groups, each = n_per_group))
X <- rep(1:n_per_group, times = n_groups)

# Random intercepts for groups
b0 <- rnorm(n_groups, 10, 3) # group variance
Y <- 2 + 1.6 * X + b0[group] + rnorm(n_groups * n_per_group, 0, 1)

tb <- tibble(group = group, X = X, Y = Y)

tb %>%
    ggplot(aes(x = X, y = Y, color = group)) +
    geom_point() +
    geom_smooth(method = "lm") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme() +
    guides() +
    labs()

# --- Frequentist Linear Mixed Model ---
# Model with random intercept for groups
lmm <- lmer(Y ~ X + (1 | group), data = tb)
su_lmm <- summary(lmm)

# Fixed effects estimates
fixef(lmm)

# Random effects
ranef(lmm)$group[,1]
var(ranef(lmm)$group[,1])
var(b0)

# Residuals, or sigma
resid(lmm) %>% var

# Variance components
Var_corr <- VarCorr(lmm)
var_ra <- as_tibble(Var_corr)$vcov[1] # random effect
var_re <- as_tibble(Var_corr)$vcov[2] # residual

var_ra / (var_ra+var_re) # Proportion of variance due to random intercept
var_re / (var_ra+var_re) # Proportion of variance due to residual variance


# --- Bayesian Hierarchical Model ---
# Same model specification in brm
brm_model <- brm(Y ~ X + (1 | group), data = tb, chains = 2, cores = 2, iter = 10000, thin = 10, seed = 123)

# Fixed effect estimates
fixef(brm_model)

# Random effects estimates
ranef(brm_model)

# Posterior predictive checks
pp_check(brm_model)

# Compute Bayesian R2
bayes_R2(brm_model)

# Variance components
## Theoretical variance components

## Variance random effect per each posterior sample
var_ra_samples <- as_draws_df(brm_model) %>%
    as_tibble() %>%
    select(starts_with("r_group")) %>%
    apply(1, var)

## Variance residuals per each posterior sample
var_re_samples <- as_draws_df(brm_model) %>%
    as_tibble() %>%
    pull(sigma) %>%
    `^`(2)

mean(var_ra_samples + var_re_samples) # total variance
mean(var_ra_samples) # Random intercept variance by group
mean(var_re_samples) # Residual variance
mean(var_ra_samples/(var_ra_samples + var_re_samples)) # proportion of total variance explained by random intecepts

hypothesis(brm_model, class = NULL, hypothesis = "sd_group__Intercept^2 / (sd_group__Intercept^2 + sigma^2)=0") %>%
    plot
