# Simulate phylogenetic with repeated measurements

library(tidyverse)
library(easystats)     # r2()
library(brms)
library(ape)           # phylogeny simulation
library(MASS)          # mvrnorm
library(posterior)
library(matrixStats)


set.seed(1)

# Simulation ----
S <- 20            # number of species (tips)
G_per_S <- 4       # number of genotypes per species
n_rep <- 8         # replicates per genotype (=> each genotype has multiple measurements)
theta <- 3         # NB/ZI-NB size
zi_logit <- -1.2   # zero-inflation logit (≈ 0.23 zeros)

# Simulate a phylogeny and build species correlation matrix (A_phylo)
tree <- ape::rtree(S)
tree <- chronos(tree, lambda = 1) # the modern penalized-likelihood chronogrammer in ape. It doesn’t require calibrations and is a clean drop-in for “make this ultrametric” workflows. lambda controls rate smoothing (higher → closer to a strict clock).
A_phylo <- cov2cor(ape::vcv(tree))                  # correlation matrix (1s on diag)
sp_levels <- rownames(A_phylo)

# Factors: species and nested genotype (species:genotype_id), with replicates
species <- factor(rep(sp_levels, each = G_per_S * n_rep), levels = sp_levels)
geno_id <- rep(rep(seq_len(G_per_S), each = n_rep), times = S)
genotype <- interaction(species, geno_id, sep = ":", drop = TRUE)
N <- length(species)

# Random intercepts (on link scale)
sd_phylo <- 0.5
sd_geno  <- 0.6
b_phylo <- drop(MASS::mvrnorm(1, mu = rep(0, S), Sigma = (sd_phylo^2) * A_phylo))
names(b_phylo) <- sp_levels
b_geno <- rnorm(nlevels(genotype), 0, sd_geno)
names(b_geno) <- levels(genotype)

# Linear predictor
beta0 <- 0.3
beta1 <- 0.7
X <- rnorm(N, 0, 1)
eta <- beta0 + beta1 * X + b_phylo[species] + b_geno[genotype]

# Responses
mu <- exp(eta)                                # for GLMMs with log link
pi0 <- plogis(zi_logit)                       # ZI probability

Y1 <- rnorm(N, mean = eta, sd = 1.0)          # Gaussian (identity link)
Y2 <- rnbinom(N, size = theta, mu = mu)       # NB2 (log link)
is_zero <- rbinom(N, 1, pi0)
Y3 <- ifelse(is_zero == 1, 0L, rnbinom(N, size = theta, mu = mu))  # ZINB

tb <- tibble(X, species, genotype, Y1, Y2, Y3, phylo = species)

# Ensure matrix row/col order matches species factor
A_ord <- A_phylo[levels(tb$species), levels(tb$species)]

# Gaussian ----
mod_g <- brm(
    Y1 ~ X + (1|gr(phylo, cov = A_ord)) + (1|species) + (1|species:genotype),
    data = tb, data2 = list(A_ord = A_ord),
    family = gaussian(),
    chains = 2, iter = 2000, cores = 2, seed = 1,
    control = list(adapt_delta = 0.95, max_treedepth = 12)
)
pp_check(mod_g)
summary(mod_g)

# r2
r2_nakagawa(mod_g) # Nakagawa R2 (link = response for Gaussian)
r2_bayes(mod_g) # bayes r2

# RE variance contributions
dr <- as_draws_df(mod_g)
sd_phy  <- dr$sd_phylo__Intercept^2
sd_sp <- dr$sd_species__Intercept^2
sd_gen <- dr$`sd_species:genotype__Intercept`^2
## 1) Latent scale ratio (identity link => tau^2)
lat_ratio <- (sd_sp^2) / (sd_phy^2 + sd_sp^2 + sd_gen^2)
hist(lat_ratio)
c(mean = mean(lat_ratio), q025 = quantile(lat_ratio, .025), q975 = quantile(lat_ratio, .975))
hyp <- hypothesis(mod_g, class= NULL, hypothesis = "sd_species__Intercept^2 / (sd_phylo__Intercept^2 + sd_species__Intercept^2 + sd_species:genotype__Intercept^2)=0")
plot(hyp)

## 2) Response-scale, link-adjusted ratio (identity => same as latent)
resp_ratio <- lat_ratio
c(mean = mean(resp_ratio), q025 = quantile(resp_ratio, .025), q975 = quantile(resp_ratio, .975))


# Negative binomial (NB2) ----
mod_nb <- brm(
    Y2 ~ X + (1|gr(phylo, cov = A_ord)) + (1|species) + (1|species:genotype),
    data = tb, data2 = list(A_ord = A_ord),
    family = negbinomial(),
    chains = 2, iter = 2000, cores = 2, seed = 1,
    control = list(adapt_delta = 0.95, max_treedepth = 12)
)
pp_check(mod_nb)
summary(mod_nb)

# R2
## Link scale
r2_nakagawa(mod_nb)
# Bayesian R2 (response scale) + hand check
y_nb <- as.numeric(get_response(mod_nb))
yhat_nb <- posterior_epred(mod_nb)        # response scale
res_nb  <- -sweep(yhat_nb, 2, y_nb)
var_mu_nb  <- rowVars(yhat_nb)
var_res_nb <- rowVars(res_nb)
R2_bayes_hand_nb <- var_mu_nb / (var_mu_nb + var_res_nb)
c(mean = mean(R2_bayes_hand_nb), q025 = quantile(R2_bayes_hand_nb, .025), q975 = quantile(R2_bayes_hand_nb, .975))
bayes_R2(mod_nb) # response scale

# RE variance contributions
dr <- as_draws_df(mod_nb)
sd_phy  <- dr$sd_phylo__Intercept^2
sd_sp <- dr$sd_species__Intercept^2
sd_gen <- dr$`sd_species:genotype__Intercept`^2

## 1) Latent scale ratio
lat_ratio <- (sd_sp^2) / (sd_phy^2 + sd_sp^2 + sd_gen^2)
hist(lat_ratio)
c(mean = mean(lat_ratio), q025 = quantile(lat_ratio, .025), q975 = quantile(lat_ratio, .975))
hyp <- hypothesis(mod_g, class= NULL, hypothesis = "sd_species__Intercept^2 / (sd_phylo__Intercept^2 + sd_species__Intercept^2 + sd_species:genotype__Intercept^2)=0")
plot(hyp)

## 2) Response-scale, link-adjusted ratio (log link => Var(exp(b)) = e^{tau^2}(e^{tau^2}-1))
var_exp_normal <- function(sd_draws) exp(sd_draws^2) * (exp(sd_draws^2) - 1)
VA <- var_exp_normal(sd_phy)
VB <- var_exp_normal(sd_sp)
VC <- var_exp_normal(sd_gen)
resp_ratio_nb <- VA / (VA + VB + VC)
c(mean = mean(resp_ratio_nb), q025 = quantile(resp_ratio_nb, .025), q975 = quantile(resp_ratio_nb, .975))

# Zero-inflated NB (ZINB) ----
mod_zi <- brm(
    bf(Y3 ~ X + (1|gr(phylo, cov = A_ord)) + (1|species) + (1|species:genotype), zi ~ 1),
    data = tb, data2 = list(A_ord = A_ord),
    family = zero_inflated_negbinomial(),
    chains = 2, iter = 2000, cores = 2, seed = 1,
    control = list(adapt_delta = 0.95, max_treedepth = 12)
)
pp_check(mod_zi)
summary(mod_zi)

# Nakagawa R2 (link scale). r2() supports brms; for ZI it returns overall R2 on link scale
r2_nakagawa(mod_zi) # not available
r2_zeroinflated(mod_zi)

# Bayesian R2 (response scale) + hand check
y_zi <- as.numeric(get_response(mod_zi))
yhat_zi <- posterior_epred(mod_zi)        # mixture expectation E[Y|X] on response scale
res_zi  <- -sweep(yhat_zi, 2, y_zi)
var_mu_zi  <- rowVars(yhat_zi)
var_res_zi <- rowVars(res_zi)
R2_bayes_hand_zi <- var_mu_zi / (var_mu_zi + var_res_zi)
c(mean = mean(R2_bayes_hand_zi), q025 = quantile(R2_bayes_hand_zi, .025), q975 = quantile(R2_bayes_hand_zi, .975))
bayes_R2(mod_zi)

# RE variance contributions
dr <- as_draws_df(mod_zi)
sd_phy  <- dr$sd_phylo__Intercept^2
sd_sp <- dr$sd_species__Intercept^2
sd_gen <- dr$`sd_species:genotype__Intercept`^2

## 1) Latent scale ratio
lat_ratio <- (sd_sp^2) / (sd_phy^2 + sd_sp^2 + sd_gen^2)
hist(lat_ratio)
c(mean = mean(lat_ratio), q025 = quantile(lat_ratio, .025), q975 = quantile(lat_ratio, .975))
hyp <- hypothesis(mod_g, class= NULL, hypothesis = "sd_species__Intercept^2 / (sd_phylo__Intercept^2 + sd_species__Intercept^2 + sd_species:genotype__Intercept^2)=0")
plot(hyp)

## 2) Response-scale, link-adjusted ratio (log link => Var(exp(b)) = e^{tau^2}(e^{tau^2}-1))
var_exp_normal <- function(sd_draws) exp(sd_draws^2) * (exp(sd_draws^2) - 1)
VA <- var_exp_normal(sd_phy)
VB <- var_exp_normal(sd_sp)
VC <- var_exp_normal(sd_gen)
resp_ratio_nb <- VA / (VA + VB + VC)
hist(resp_ratio_nb)
c(mean = mean(resp_ratio_nb), q025 = quantile(resp_ratio_nb, .025), q975 = quantile(resp_ratio_nb, .975))



if (F) {

# Create phylogeny and covariance matrix
## The diagonals of the cov matrix are the total branch lengths from the root to tips for each species
## The off-diagonals are the total branch lengths from the root to the most recent common ancestor (MRCA) of two species
phylo <- rtree(n = n_species)
phylo <- chronopl(phylo, lambda = 1) # turn the genetic distance based tree to a time calibrated tree
species <- phylo$tip.label
plot(phylo, main = "Mock Phylogeny")
phylo_cov <- vcv(phylo)
#cov_scaled <- phylo_cov/sum(diag(phylo_cov)) # rescaled the cov to have a total variance = 1

# Simulate phylogenetic effects
# mvrnorm: Simulate from a Multivariate Normal Distribution
phylo_effects <- MASS::mvrnorm(n = 1, mu = rep(0, nrow(phylo_cov)), Sigma = phylo_cov)
names(phylo_effects) <- species
#apply(phylo_effects, 2, var) %>% sum # this should be close to 1 if n = 1000 in m

# Simulate species effect
sp_group <- rep(species, each = measurements_per_species)
species_effects <- rnorm(n_species, 0, 1)
names(species_effects) <- species

# Simulate response in negative binomial
X <- rnorm(n_species*measurements_per_species, 0, 1) # predictor
eta <- 1 + 0.8 * X + phylo_effects[sp_group] + 0.1*species_effects[sp_group]
# Set a high dispersion parameter to reduce randomness generated by NB sampling
Y <- rnbinom(n_species * measurements_per_species, mu = exp(eta), size = 100) # dispersion parameter k. var(NB) = mu + mu^2/k

tb <- tibble(species = sp_group, phylo = sp_group, X = X, Y = Y)
hist(Y, breaks = 50)

# species mean
tb <- tb %>%
    group_by(phylo) %>%
    mutate(X_bet = mean(X), X_wit = X - X_bet)


# First model using the species mean ----
mod <- brm(
    Y ~ X_bet + X_wit + (1|species) + (1|gr(phylo, cov = phylo_cov)),
    data = tb,
    family = "negbinomial2",
    data2 = list(phylo_cov = phylo_cov),
    chains = 4, cores = 4, iter = 20000, thin = 10, seed = 123,
    control = list(adapt_delta = 0.99)
)

pp_check(mod)
describe_posterior(mod)

# phylogenetic signal
# In NB, there is a sigma estimate to describe variability beyond poisson (which equals mean)
hyp <- hypothesis(mod, class = NULL, "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sd_species__Intercept^2 + sigma^2) = 0")
plot(hyp)

hyp <- hypothesis(mod, class = NULL, "sd_species__Intercept^2 / (sd_phylo__Intercept^2 + sd_species__Intercept^2 + sigma^2) = 0")
plot(hyp)

# Variance components
set.seed(1)
pp <- posterior_predict(mod)
pp_mean <- colMeans(pp) # mean posterior predictor
var_res <- var(Y - pp_mean)

## Proportion of total variance explained not by residual
1 - var_res / var(Y) # should be similar to bayes R2
bayes_R2(mod)

## variance explained by random effects (phylogeny)
pp <- as_draws_df(mod) %>% as_tibble()
ran_samples_phylo <- select(pp, starts_with("r_phylo"))
ran_samples_sp <- select(pp, starts_with("r_species"))
var_ran_phylo <- apply(ran_samples_phylo, 1, var) # posterior distribution of variance of random effects
var_ran_sp <- apply(ran_samples_sp, 1, var) # posterior distribution of variance of random effects

var_tot <- var_ran_phylo + var_ran_sp + pp$sigma^2 # total variance
ps <- var_ran_phylo / var_tot
hist(ps); mean(ps); quantile(ps, c(0.05, 0.5, 0.95)) # this should be similar to hypothesis results

# Second model using the raw values ----
mod2 <- brm(
    Y ~ X + (1|species) + (1|gr(phylo, cov = phylo_cov)),
    data = tb,
    family = "negbinomial2",
    data2 = list(phylo_cov = phylo_cov),
    chains = 4, cores = 4, iter = 20000, thin = 10, seed = 123,
    control = list(adapt_delta = 0.99)
)

pp_check(mod2)
describe_posterior(mod2)

# In NB, there is a sigma estimate to describe variability beyond poisson (which equals mean)
hyp <- hypothesis(mod2, class = NULL, "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sd_species__Intercept^2 + sigma^2) = 0")
plot(hyp)
}


