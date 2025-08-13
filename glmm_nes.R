# =========================================================
# Simulation + Fits + R2 (Nakagawa, Bayes/Empirical) + part R2
# Nested REs: site nested in region  =>  (1|region/site)
# Families: gaussian, NB2, ZINB
# =========================================================

library(tidyverse)
library(easystats)
library(glmmTMB)
library(brms)
library(partR2) # for part R2, basically the R2 for fixed effect predictors
set.seed(1)

# Simulation ----
R = 8 # number of regions
S_per_R = 3 # number of sites per regions
n_per_site = 10 # number of replicates per site
theta = 3               # NB/ZI-NB size
zi_logit = -1.2        # zero-inflation logit (≈ 0.23 zeros)

# Factors
region <- factor(rep(paste0("r", seq_len(R)), each = S_per_R * n_per_site))
site_id <- rep(rep(seq_len(S_per_R), each = n_per_site), times = R)
site <- interaction(region, site_id, sep = ":", drop = TRUE)
N <- length(region)

# Random intercepts
b_region <- rnorm(R, 0, 0.6) # region random intercept
names(b_region) <- levels(region)
b_site <- rnorm(nlevels(site), 0, 0.9) # site random intercept
names(b_site) <- levels(site)

# Linear predictor
X <- rnorm(N, 0, 1)
X2 <- rnorm
eta <- 0.5 + 0.8*X + b_region[region] + b_site[site]

# Link
mu <- exp(eta)
pi0 <- plogis(-1.2)                      # zero-inflation logit (≈ 0.23 zeros)
inflation <- rbinom(N, 1, 1-pi0)

tb <- tibble(X, region, site,
       Y1 = rnorm(N, mean = eta, sd = 1), # Gaussian
       Y2 = rnbinom(N, size = theta, mu = mu), # NB
       Y3 = inflation * Y2 # ZINB
)


# Gaussian ----
## Frequentist ----
mod_g_f <- lme4::lmer(Y1 ~ X + (1|region) + (1|region:site), data = tb)
partR2(mod_g_f, partvars = c("X")) # part R2 only suppory merMod

mod_g_f <- glmmTMB(Y1 ~ X + (1|region) + (1|region:site), data = tb)

# Variance components and glmm R2
var_fix <- get_variance_fixed(mod_g_f)
var_ran <- get_variance_random(mod_g_f)
var_res <- get_variance_residual(mod_g_f)
(var_fix + var_ran) / (var_fix + var_ran + var_res) # conditional R2
(var_fix) / (var_fix + var_ran + var_res) # marginal R2
r2_nakagawa(mod_g_f)

# Empirical R2 (analogous to bayesian R2)
y_pred <- predict(mod_g_f, type = "response") # prediction
y <- get_response(mod_g_f) # data
var(y_pred) / (var(y_pred) + var(y - y_pred))


## Bayesian ----
mod_g_b <- brm(Y1 ~ X + (1|region/site), data = tb, family = gaussian, chains = 2, iter = 2000, cores = 2, seed = 1, control = list(adapt_delta = 0.95, max_treedepth = 12))

# Relative contribution of random effects
dr <- as_tibble(as_draws_df(mod_g_b))
## 1. Latent scale
var_r <- dr$sd_region__Intercept^2
var_s <- dr$`sd_region:site__Intercept`^2
vv <- var_r / (var_r + var_s)
c(mean(vv), quantile(vv, probs = c(0.025, 0.975)))
hypothesis(dr, hypothesis = "sd_region__Intercept^2 / (sd_region__Intercept^2 + sd_region:site__Intercept^2) =0")

## 2. response scale. link-adjusted: exp(b) variance ratio. For gaussan it's identity
var_r_adj <- var_r
var_s_adj <- var_s
vv <- var_r_adj / (var_r_adj + var_s_adj)
c(mean(vv), quantile(vv, probs = c(0.025, 0.975)))

## R2
# 1. Variance components and glmm R2
var_fix <- get_variance_fixed(mod_g_b)
var_ran <- get_variance_random(mod_g_b)
var_res <- get_variance_residual(mod_g_b)
(var_fix + var_ran) / (var_fix + var_ran + var_res) # conditional R2
(var_fix) / (var_fix + var_ran + var_res) # marginal R2
r2_nakagawa(mod_g_b)

# 2. Bayesian r2
y <- as.numeric(get_response(mod_g_b)) # data
y_draws <- posterior_linpred(mod_g_b) # pp predicted y
res_draws <- -sweep(y_draws, 2, y)
## variance
var_ypred <- apply(y_draws, 1, var) # var in  predicted
var_res <- apply(res_draws, 1, var) # var in residuals. Not sigma (parametric residuals)
R2_bayes_hand <- var_ypred / (var_ypred + var_res)
c(mean(R2_bayes_hand), sd(R2_bayes_hand), quantile(target_var, probs = c(0.025, 0.975)))
bayes_R2(mod_g_b)

# # 3. Shapley r2. Additive and robust
# # Example for two nested REs: A = (1|region), B = (1|region:site)
# # Full model: Y ~ X + (1|region) + (1|region:site)
# R2_full <- bayes_R2(mod_g_b, summary = F)[,1]
# R2_Aonly <- bayes_R2(mod_g_b, re_formula = ~(1|region), summary = F)[,1]
# R2_Bonly <- bayes_R2(mod_g_b, re_formula = ~(1|region:site), summary = F)[,1]
# R2_none  <- bayes_R2(mod_g_b, re_formula = NA, summary = F)[,1]
#
# # Shapley for two REs (additive, no refit needed for RE toggling)
# phi_A <- 0.5 * ((R2_Aonly - R2_none) + (R2_full - R2_Bonly))
# phi_B <- 0.5 * ((R2_Bonly - R2_none) + (R2_full - R2_Aonly))
# mean(phi_A / (phi_A + phi_B))
# 0.6^2 / (0.6^2 + 0.9^2)
#
# phi_A / (phi_A + phi_B)




# Negative binomial ----
## Frequentist ----
mod_nb_f <- glmmTMB(Y2 ~ X + (1|region/site), data = tb, family = nbinom2)

# Variance components and glmm R2
var_fix <- get_variance_fixed(mod_nb_f)
var_ran <- get_variance_random(mod_nb_f)
var_res <- get_variance_residual(mod_nb_f)
(var_fix + var_ran) / (var_fix + var_ran + var_res) # conditional R2
(var_fix) / (var_fix + var_ran + var_res) # marginal R2
r2_nakagawa(mod_nb_f)

# Empirical R2 (analogous to bayesian R2) at response scale
y_pred <- predict(mod_nb_f, type = "response") # prediction
y <- get_response(mod_nb_f) # data
var(y_pred) / (var(y_pred) + var(y - y_pred))

## Bayesian ----
mod_nb_b <- brm(Y2 ~ X + (1|region/site), data = tb, family = negbinomial, chains = 2, iter = 2000, cores = 2, seed = 1, control = list(adapt_delta = 0.95, max_treedepth = 12))

# Relative contribution of random effects
dr <- as_tibble(as_draws_df(mod_nb_b))
## 1. Latent scale
var_r <- dr$sd_region__Intercept^2
var_s <- dr$`sd_region:site__Intercept`^2
vv <- var_r / (var_r + var_s)
c(mean(vv), quantile(vv, probs = c(0.025, 0.975)))
hypothesis(dr, hypothesis = "sd_region__Intercept^2 / (sd_region__Intercept^2 + sd_region:site__Intercept^2) =0")

## 2. response scale. link-adjusted: exp(b) variance ratio
var_r_adj <- exp(var_r) * (exp(var_r) - 1) # for log link function, VA <- exp(sdA^2) * (exp(sdA^2) - 1)
var_s_adj <- exp(var_s) * (exp(var_s) - 1)
vv <- var_r_adj / (var_r_adj + var_s_adj)
c(mean(vv), quantile(vv, probs = c(0.025, 0.975)))

## R2
# 1. Variance components and glmm R2
var_fix <- get_variance_fixed(mod_nb_b)
var_ran <- get_variance_random(mod_nb_b)
var_res <- get_variance_residual(mod_nb_b)
(var_fix + var_ran) / (var_fix + var_ran + var_res) # conditional R2
(var_fix) / (var_fix + var_ran + var_res) # marginal R2
r2_nakagawa(mod_nb_b)

# 2. Bayesian r2 at response scale
y <- as.numeric(get_response(mod_nb_b)) # data
y_draws <- posterior_epred(mod_nb_b) # pp predicted y at response scale
res_draws <- -sweep(y_draws, 2, y)
## variance
var_ypred <- apply(y_draws, 1, var) # var in  predicted
var_res <- apply(res_draws, 1, var) # var in residuals. Not sigma (parametric residuals)
R2_bayes_hand <- var_ypred / (var_ypred + var_res)
c(mean(R2_bayes_hand), sd(R2_bayes_hand), quantile(target_var, probs = c(0.025, 0.975)))
bayes_R2(mod_nb_b)


# ZINB ----
## Frequenties ----
mod_zinb_f <- glmmTMB(Y3 ~ X + (1|region/site), data = tb, family = nbinom2, ziformula = ~1)

# Variance components and glmm R2
var_fix <- get_variance_fixed(mod_zinb_f)
var_ran <- get_variance_random(mod_zinb_f)
var_res <- get_variance_residual(mod_zinb_f)
(var_fix + var_ran) / (var_fix + var_ran + var_res) # conditional R2
(var_fix) / (var_fix + var_ran + var_res) # marginal R2
r2_nakagawa(mod_zinb_f)

# Empirical R2 (analogous to bayesian R2) at response scale
y_pred <- predict(mod_zinb_f, type = "response") # prediction
y <- get_response(mod_zinb_f) # data
var(y_pred) / (var(y_pred) + var(y - y_pred))

r2_zeroinflated(mod_zinb_f)

## Bayesian ----
mod_zinb_b <- brm(bf(Y3 ~ X + (1|region/site), zi ~ 1), data = tb, family = zero_inflated_negbinomial, chains = 2, iter = 2000, cores = 2, seed = 1, control = list(adapt_delta = 0.95, max_treedepth = 12))

# Relative contribution of random effects
dr <- as_tibble(as_draws_df(mod_zinb_b))
## 1. Latent scale
var_r <- dr$sd_region__Intercept^2
var_s <- dr$`sd_region:site__Intercept`^2
vv <- var_r / (var_r + var_s)
c(mean(vv), quantile(vv, probs = c(0.025, 0.975)))
hypothesis(dr, hypothesis = "sd_region__Intercept^2 / (sd_region__Intercept^2 + sd_region:site__Intercept^2) =0")

## 2. response scale. link-adjusted: exp(b) variance ratio
var_r_adj <- exp(var_r) * (exp(var_r) - 1) # for log link function, VA <- exp(sdA^2) * (exp(sdA^2) - 1)
var_s_adj <- exp(var_s) * (exp(var_s) - 1)
vv <- var_r_adj / (var_r_adj + var_s_adj)
c(mean(vv), quantile(vv, probs = c(0.025, 0.975)))

## R2
# 1. Variance components and glmm R2 -> does not work with zero inflation
var_fix <- get_variance_fixed(mod_zinb_b) # does not work
var_ran <- get_variance_random(mod_zinb_b)
var_res <- get_variance_residual(mod_zinb_b)# does not work
(var_fix + var_ran) / (var_fix + var_ran + var_res) # conditional R2
(var_fix) / (var_fix + var_ran + var_res) # marginal R2
r2_nakagawa(mod_zinb_b)
r2(mod_zinb_b)

# 2. Bayesian r2 at response scale
y <- as.numeric(get_response(mod_zinb_b)) # data
y_draws <- posterior_epred(mod_zinb_b) # pp predicted y at response scale
res_draws <- -sweep(y_draws, 2, y)
## variance
var_ypred <- apply(y_draws, 1, var) # var in  predicted
var_res <- apply(res_draws, 1, var) # var in residuals. Not sigma (parametric residuals)
R2_bayes_hand <- var_ypred / (var_ypred + var_res)
c(mean(R2_bayes_hand), sd(R2_bayes_hand), quantile(target_var, probs = c(0.025, 0.975)))
bayes_R2(mod_zinb_b)
