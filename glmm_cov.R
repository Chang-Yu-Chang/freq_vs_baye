#' Simulate glmm with two random effects: groups and covariance matrix
#' do it in gaussian and negative binomial, and ZINB

library(tidyverse)
library(easystats)
library(MASS)
library(brms)

# Simulation ----
n_groups    <- 10
n_per_group <- 10
n_total     <- n_groups * n_per_group

group_names <- paste0("g", seq_len(n_groups))
sp_group    <- factor(rep(group_names, each = n_per_group), levels = group_names)

X <- rnorm(n_total, 0, 1)  # predictor

# --- correlation matrix C for the structured RE ---
A      <- matrix(rnorm(n_groups^2), n_groups, n_groups)
Sigma0 <- crossprod(A)                           # PD covariance
D      <- diag(1 / sqrt(diag(Sigma0)))
C      <- D %*% Sigma0 %*% D                     # CORRELATION matrix
dimnames(C) <- list(group_names, group_names)
C_ord <- C[levels(sp_group), levels(sp_group)]   # ensure matching order

# --- simulate two RE components on same levels ---
sd_cov <- 0.9                                    # structured RE SD
sd_iid <- 0.7                                    # i.i.d. RE SD

b_iid <- rnorm(n_groups, 0, sd_iid)              # i.i.d. random intercepts
names(b_iid) <- group_names

SigmaRE <- (sd_cov^2) * C                        # structured covariance
b_cov   <- drop(mvrnorm(1, mu = rep(0, n_groups), Sigma = SigmaRE))
names(b_cov) <- group_names

# --- linear predictor (your spec) ---
eta <- 1 + 0.8 * X + b_iid[sp_group] + b_cov[sp_group]


# ===== Gaussian LMM =====
sigma_eps <- 1.0
Y_g <- rnorm(n_total, mean = eta, sd = sigma_eps)
tb_g <- tibble(sp_group, X, Y = Y_g, cov_group = sp_group)

mod_g <- brm(
    Y ~ X + (1 | sp_group) + (1 | gr(cov_group, cov = C_ord)),
    data  = tb_g,
    data2 = list(C_ord = C_ord),      # name must match 'cov = C_ord' in formula
    family = gaussian(),
    chains = 2, cores = 2, iter = 4000, seed = 123,
    control = list(adapt_delta = 0.95, max_treedepth = 12)
)

# Model estimates
mod_g
pp_check(mod_g)
model_parameters(mod_g)
fixef(mod_g)
ranef(mod_g) # BLUPs
icc(mod_g) # ICC

# Estimate marginal means and contrasts
estimate_means(mod_g, by = "sp_group")
estimate_contrasts(mod_g, contrast = "sp_group", comparison = "reference")


# Variance components
VarCorr(mod_g) # this only provides random effects and residuals

# Use the posterior draws
dr <- as_tibble(as_draws_df(mod_g))
## fixed effect variance on link scale
eta_fix_draws <- posterior_linpred(mod_g, re_formula = NA, transform = FALSE)
var_fix <- apply(eta_fix_draws, 1, var)
## random effect variance
var_iid <- dr[["sd_sp_group__Intercept"]]^2 # group
var_cov <- dr[["sd_cov_group__Intercept"]]^2 # cov
## Residual (distribution-specific) variance
var_res <- dr[["sigma"]]^2
## total variance
var_total <- var_fix + var_iid + var_cov + var_res

## Relative var of random effects
target_var <- var_iid / (var_iid + var_cov)
c(mean(target_var), sd(target_var), quantile(target_var, probs = c(0.025, 0.975)))
hypothesis(mod_g, class = NULL, "sd_sp_group__Intercept^2 / (sd_sp_group__Intercept^2 + sd_cov_group__Intercept^2) =0")

## nagakawa R2
R2_marg_draws <- var_fix / var_total
R2_cond_draws <- (var_fix + var_iid + var_cov) / var_total
mean(R2_cond_draws)
mean(R2_marg_draws)
r2_nakagawa(mod_g)

## Bayes_r2. Expected values per draw (includes fixed + random effects)
mu_draws <- posterior_epred(mod_g)         # at response scale
var_mu   <- apply(mu_draws, 1, var)        # Var(mu) per draw
target_var <- var_mu / (var_mu + var_res)
c(mean(target_var), sd(target_var), quantile(target_var, probs = c(0.025, 0.975)))
bayes_R2(mod_g)

# ===== Negative Binomial GLMM (log link) =====
# Here, eta is on the log scale: mu = exp(eta)
theta <- 3                           # NB 'size' (Var = mu + mu^2/theta)
mu <- exp(eta)                       # eta is on log scale here
Y_nb <- rnbinom(n_total, mu = mu, size = theta)
tb_nb <- tibble(sp_group, X, Y = Y_nb, cov_group = sp_group)

mod_nb <- brm(
    Y ~ X + (1 | sp_group) + (1 | gr(cov_group, cov = C_ord)),
    data  = tb_nb,
    data2 = list(C_ord = C_ord),
    family = negbinomial(),            # NB2 in brms
    chains = 2, cores = 2, iter = 4000, seed = 123,
    control = list(adapt_delta = 0.95, max_treedepth = 12)
)

mod_nb
model_parameters(mod_nb)
pp_check(mod_nb)
fixef(mod_nb)
ranef(mod_nb) # BLUPs
icc(mod_nb) # ICC

# Estimate marginal means and contrasts
estimate_means(mod_nb, by = "sp_group")
estimate_contrasts(mod_nb, contrast = "sp_group", comparison = "reference")

# Variance components
VarCorr(mod_nb) # this only provides random effects. Note that it does not have residuals

# Use the posterior draws
dr <- as_tibble(as_draws_df(mod_nb))

## fixed effect variance on link scale
eta_fix_draws <- posterior_linpred(mod_nb, re_formula = NA, transform = FALSE)
var_fix <- apply(eta_fix_draws, 1, var)
## random effect variance
var_iid <- dr[["sd_sp_group__Intercept"]]^2 # group
var_cov <- dr[["sd_cov_group__Intercept"]]^2 # cov
## DISTRIBUTION-SPECIFIC ("RESIDUAL") VARIANCE on link scale (Nakagawa, lognormal) ---
# NB2 in brms uses 'shape' (theta in Var = mu + mu^2/theta)
theta_draw <- dr[["shape"]]
b0_draw    <- dr[["b_Intercept"]] # null model
mu0_draw <- exp(b0_draw)
# Null-model mean on response scale per draw, adjusted by 0.5 * sum(tau^2)
#re_int_var_sum <- var_iid + var_cov
# link-scale residual variance per draw (lognormal approximation)
var_res <- log1p(1 / mu0_draw + 1 / theta_draw)
# --- TOTAL VARIANCE (link scale) ---
var_total <- var_fix + var_iid + var_cov + var_res

## Relative var of random effects
target_var <- var_iid / (var_iid + var_cov)
c(mean(target_var), sd(target_var), quantile(target_var, probs = c(0.025, 0.975)))
hypothesis(mod_nb, class = NULL, "sd_sp_group__Intercept^2 / (sd_sp_group__Intercept^2 + sd_cov_group__Intercept^2) =0")

## Bayes_r2. Response scale
# Observations
y <- as.numeric(get_response(mod_nb))
# Expected values per draw
y_draws <- posterior_epred(mod_nb)
# Empirical residuals per draw (response scale)
e_draws <- -sweep(y_draws, 2, y)
# per-draw variances across observations
var_ypred <- apply(y_draws, 1, var)       # matches brms' .bayes_R2()
var_e     <- apply(e_draws, 1, var)

R2_bayes_hand <- var_ypred / (var_ypred + var_e)
c(mean(R2_bayes_hand), sd(R2_bayes_hand), quantile(target_var, probs = c(0.025, 0.975)))
bayes_R2(mod_nb)



# ===== Zero-inflated Negative Binomial GLMM (log link) =====
p_zero <- 0.2 # probability of zeros
inflate <- rbinom(n_total, size = 1, prob = 1-p_zero)
theta <- 3                           # NB 'size' (Var = mu + mu^2/theta)
mu <- exp(eta)                       # eta is on log scale here
Y_nb <- rnbinom(n_total, mu = exp(eta), size = theta)
Y_zinb <- Y_nb*inflate
tb_zinb <- tibble(sp_group, X, Y = Y_zinb, cov_group = sp_group)

mod_zinb <- brm(
    bf(Y ~ X + (1 | sp_group) + (1 | gr(cov_group, cov = C_ord)),
       zi ~ 1),
    data  = tb_zinb,
    data2 = list(C_ord = C_ord),
    family = zero_inflated_negbinomial,            # NB2 in brms
    chains = 2, cores = 2, iter = 4000, seed = 123,
    control = list(adapt_delta = 0.95, max_treedepth = 12)
)


mod_zinb
model_parameters(mod_zinb)
pp_check(mod_zinb)
fixef(mod_zinb)
ranef(mod_zinb) # BLUPs
icc(mod_zinb) # ICC # can't do it in ZI

# Estimate marginal means and contrasts
estimate_means(mod_zinb, by = "sp_group")
estimate_contrasts(mod_zinb, contrast = "sp_group", comparison = "reference")

# Variance components
VarCorr(mod_zinb) # this only provides random effects. Note that it does not have residuals

# Use the posterior draws
dr <- as_tibble(as_draws_df(mod_zinb))
names(dr)

## fixed effect variance on link scale
eta_fix_draws <- posterior_linpred(mod_zinb, re_formula = NA, transform = FALSE)
var_fix <- apply(eta_fix_draws, 1, var)
## random effect variance
var_iid <- dr[["sd_sp_group__Intercept"]]^2 # group
var_cov <- dr[["sd_cov_group__Intercept"]]^2 # cov
## DISTRIBUTION-SPECIFIC ("RESIDUAL") VARIANCE on link scale (Nakagawa, lognormal) ---
# NB2 in brms uses 'shape' (theta in Var = mu + mu^2/theta)
theta_draw <- dr[["shape"]]
b0_draw    <- dr[["b_Intercept"]] # null model
mu0_draw <- exp(b0_draw)
# Null-model mean on response scale per draw, adjusted by 0.5 * sum(tau^2)
#re_int_var_sum <- var_iid + var_cov
# link-scale residual variance per draw (lognormal approximation)
var_res <- log1p(1 / mu0_draw + 1 / theta_draw)
# --- TOTAL VARIANCE (link scale) ---
var_total <- var_fix + var_iid + var_cov + var_res

## Relative var of random effects
target_var <- var_iid / (var_iid + var_cov)
c(mean(target_var), sd(target_var), quantile(target_var, probs = c(0.025, 0.975)))
hypothesis(mod_zinb, class = NULL, "sd_sp_group__Intercept^2 / (sd_sp_group__Intercept^2 + sd_cov_group__Intercept^2) =0")

## Bayes_r2. Response scale
# Observations
y <- as.numeric(get_response(mod_zinb))
# Expected values per draw
y_draws <- posterior_epred(mod_zinb)
# Empirical residuals per draw (response scale)
e_draws <- -sweep(y_draws, 2, y)
# per-draw variances across observations
var_ypred <- apply(y_draws, 1, var)       # matches brms' .bayes_R2()
var_e     <- apply(e_draws, 1, var)

R2_bayes_hand <- var_ypred / (var_ypred + var_e)
c(mean(R2_bayes_hand), sd(R2_bayes_hand), quantile(target_var, probs = c(0.025, 0.975)))
bayes_R2(mod_zinb)
