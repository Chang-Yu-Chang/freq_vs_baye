# Simulate phylogenetic mixed model with repeated measurements

library(tidyverse)
library(easystats)
library(ape)
library(brms)
#library(MASS) # for multivariate gaussian

# One measurement per species ----
# Parameters
set.seed(123)
n_species <- 20
#measurements_per_species <- 15

# Create phylogeny
phylo <- rtree(n = n_species)
species <- phylo$tip.label
plot(phylo, main = "Mock Phylogeny")
phylo_cov <- vcv(phylo)

# Simulate phylogenetic effects
# mvrnorm: Simulate from a Multivariate Normal Distribution
phylo_effects <- MASS::mvrnorm(n = 1, mu = rep(0, nrow(phylo_cov)), Sigma = phylo_cov)
names(phylo_effects) <- species
sd(phylo_effects) # the phylo sd

# Generate data: one measurements per species
X <- rnorm(n_species, 0, 1)
Y <- 2 + 5 * X + 10*phylo_effects[species] + rnorm(n_species, 0, 1)
tb <- tibble(phylo = species, X = X, Y = Y)

# Fit the model with phylocentric random effect
mod <- brm(
    Y ~ X + (1|gr(phylo, cov = phylo_cov)),
    data = tb,
    data2 = list(phylo_cov = phylo_cov),
    chains = 2, cores = 2, iter = 10000, thin = 10, seed = 123
)

pp_check(mod)
describe_posterior(mod)

# phylogenetic signal
hyp <- hypothesis(mod, class = NULL, "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sigma^2) = 0")
plot(hyp)

# Variance components
set.seed(1)
pp <- posterior_predict(mod)
pp_mean <- colMeans(pp) # mean posterior predictor
var_res2 <- var(Y - pp_mean)

## Proportion of total variance explained not by residual
1 - var_res2 / var(Y) # should be similar to bayes R2
bayes_R2(mod)

## variance explained by random effects (phylogeny)
ran_samples <- as_draws_df(mod) %>%
    as_tibble() %>%
    select(starts_with("r_phylo"))

var_ran <- apply(ran_samples, 1, var) # posterior distribution of variance of random effects
var_tot <- var(Y) # total variance
ps <- (var_ran / var_tot)
hist(ps)
mean(ps)
quantile(ps, c(0.05, 0.5, 0.95))


# Repeated measurements per species ----
# Parameters
set.seed(123)
n_species <- 10
measurements_per_species <- 10

# Create phylogeny
phylo <- rtree(n = n_species)
species <- phylo$tip.label
plot(phylo, main = "Mock Phylogeny")
phylo_cov <- vcv(phylo)

# Simulate phylogenetic effects
# mvrnorm: Simulate from a Multivariate Normal Distribution
phylo_effects <- MASS::mvrnorm(n = 1, mu = rep(0, nrow(phylo_cov)), Sigma = phylo_cov)
names(phylo_effects) <- species

# Simulate species effect
# Simulate response
sp_group <- rep(species, each = measurements_per_species)
X <- rnorm(n_species*measurements_per_species, 0, 1)
Y <- 2 + 5 * X + 10*phylo_effects[sp_group] + rnorm(n_species*measurements_per_species, 0, 1)
tb <- tibble(phylo = sp_group, X = X, Y = Y)

# species mean
tb <- tb %>%
    group_by(phylo) %>%
    mutate(X_bet = mean(X), X_wit = X - X_bet)

# Fit the model with phylocentric random effect
mod2 <- brm(
    Y ~ X_bet + X_wit + (1|gr(phylo, cov = phylo_cov)),
    data = tb,
    data2 = list(phylo_cov = phylo_cov),
    chains = 2, cores = 2, iter = 10000, thin = 10, seed = 123
)

pp_check(mod2)
describe_posterior(mod2)

# phylogenetic signal
hyp <- hypothesis(mod2, class = NULL, "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sigma^2) = 0")
plot(hyp)

# Variance components
set.seed(1)
pp <- posterior_predict(mod2)
pp_mean <- colMeans(pp) # mean posterior predictor
var_res2 <- var(Y - pp_mean)

## Proportion of total variance explained not by residual
1 - var_res2 / var(Y) # should be similar to bayes R2
bayes_R2(mod2)

## variance explained by random effects (phylogeny)
pp <- as_draws_df(mod2) %>%
    as_tibble()
ran_samples <- select(pp, starts_with("r_phylo"))
var_ran <- apply(ran_samples, 1, var) # posterior distribution of variance of random effects
var_tot <- var_ran + pp$sigma^2 # total variance
ps <- var_ran / var_tot
hist(ps)
mean(ps)
quantile(ps, c(0.05, 0.5, 0.95))



plot(phylo)

tb %>%
    mutate(phylo = factor(phylo, species)) %>%
    ggplot() +
    geom_col(data = distinct(tb, phylo, X_bet), aes(x = phylo, y = X_bet)) +
    geom_jitter(aes(x = phylo, y = X), width = .1) +
    coord_flip(clip = "off") +
    theme_bw() +
    theme() +
    guides() +
    labs()
