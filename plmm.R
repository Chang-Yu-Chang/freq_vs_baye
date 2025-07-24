# Simulate phylogenetic mixed model with repeated measurements

library(tidyverse)
library(easystats)
library(ape)
library(brms)
#library(MASS) # for sampling from multivariate gaussian

# One measurement per species ----
set.seed(123)
n_species <- 100

# Create phylogeny
phylo <- rtree(n = n_species)
phylo <- chronopl(phylo, lambda = 1) # so the tree is ultrametric
species <- phylo$tip.label
plot(phylo, main = "Mock Phylogeny")

# Simulate phylogenetic effects
phylo_cov <- vcv(phylo) # the diagonals should all be 1, meaning the variance is 1
sd_phylo <- 5 # the sd of brownian motion
phylo_effects <- MASS::mvrnorm(n = 1, mu = rep(0, n_species), Sigma = phylo_cov * sd_phylo^2)

# Generate data: one measurements per species
Y <- 5 + phylo_effects + rnorm(n_species, 0, 1)
tb <- tibble(phylo = species, Y = Y)

# Fit the model with phylocentric random effect
mod <- brm(
    Y ~ 1 + (1|gr(phylo, cov = phylo_cov)),
    data = tb,
    data2 = list(phylo_cov = phylo_cov),
    chains = 2, cores = 2, iter = 10000, thin = 10, seed = 123
)
mod
pp_check(mod)

# phylogenetic signal
hyp <- hypothesis(mod, class = NULL, "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sigma^2) = 0")
plot(hyp)

if (F) {
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

}


# Repeated measurements per species ----
set.seed(123)
n_species <- 50
measurements_per_species <- 10
n_total <- n_species*measurements_per_species

# Create phylogeny
phylo <- rtree(n_species)
phylo <- chronopl(phylo, lambda = 1)
species <- rep(phylo$tip.label, each = measurements_per_species)
plot(phylo, main = "Mock Phylogeny")

# Simulate phylogenetic effects
phylo_cov <- vcv(phylo) # the diagonals should all be 1, meaning the variance is 1
sd_phylo <- 5 # the sd of brownian motion
phylo_effects <- MASS::mvrnorm(n = 1, mu = rep(0, n_species), Sigma = phylo_cov * sd_phylo^2)


# Generate data: repeated measurements per species
X <- rnorm(n_total, 0, 1)
Y <- 5 + 2*X + phylo_effects[species]
tb <- tibble(phylo = species, X = X, Y = Y) %>%
    # species mean
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
mod2

# phylogenetic signal
hyp <- hypothesis(mod2, class = NULL, "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sigma^2) = 0")
plot(hyp)


if (F) {

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
}
