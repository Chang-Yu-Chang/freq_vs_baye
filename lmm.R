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
b0 <- rnorm(n_groups, 10, 3) # group variance  the true RE variance is 9
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
var(ranef(lmm)$group[,1]) # conditional modes (BLUPs) and are shrunk toward 0, so it's always smaller than the true/group‐level variance
var(b0)

# Residuals, or sigma
resid(lmm) %>% var
sigma(lmm)

# Variance components
vc <- VarCorr(lmm)
var_u  <- attr(vc$group, "stddev")^2   # random intercept variance
var_e  <- sigma(lmm)^2                 # residual variance
c(var_u = var_u, var_e = var_e, prop_u = var_u/(var_u+var_e))

var_ra <- as_tibble(vc)$vcov[1] # random effect
var_re <- as_tibble(vc)$vcov[2] # residual

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

# Variance components
## Theoretical variance components
## Variance random effect per each posterior sample. This is using the parametric sigma^2 instead of empirical residuals
draws <- as_draws_df(brm_model)
var_ra_draws <- draws$sd_group__Intercept^2   # random-intercept variance
var_re_draws <- draws$sigma^2                 # residual variance

mean(var_ra_draws + var_re_draws) # total variance
mean(var_ra_draws) # Random intercept variance by group
mean(var_re_draws) # Residual variance
mean(var_ra_draws/(var_ra_draws + var_re_draws)) # proportion of total variance explained by random intecepts

hypothesis(brm_model, class = NULL, hypothesis = "sd_group__Intercept^2 / (sd_group__Intercept^2 + sigma^2)=0")

## Compute Bayesian R2. Empirical R2
mu <- posterior_epred(brm_model)
y  <- model.response(model.frame(brm_model)) # prediction
var_mu   <- apply(mu, 1, var)
var_res  <- apply(mu, 1, function(m) var(y - m))     # empirical residual var
r2_draws <- var_mu / (var_mu + var_res)

c(mean = mean(r2_draws), sd = sd(r2_draws),
  q2.5 = quantile(r2_draws, .025), q50 = quantile(r2_draws, .5),
  q97.5 = quantile(r2_draws, .975))
bayes_R2(brm_model)  # should match (mean vs median may differ slightly)

