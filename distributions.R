#' Just to learn some statistic distributions

library(tidyverse)

# R functions ----
# In R, every distribution comes as a family of four functions with the same suffix:
# dname: density / mass function → “what’s the height at x?”
# pname: cumulative distribution function (CDF) → “what’s P(X ≤ x)?”
# qname: quantile function (inverse CDF) → “what x gives me P(X ≤ x) = p?”
# rname: random generation → “give me n samples”
# Density at x = 1.96 for N(0,1)
dnorm(1.96)
# P(|Z| <= 1.96)
pnorm(1.96) # 0.975
# Two-sided z-test tail area: P(|Z| ≥ 1.96)
pnorm(1.96, lower.tail = FALSE) # 0.025
2*pnorm(1.96, lower.tail = FALSE) # this gives right tail prob without 1-p
# 97.5th percentile (useful for 95% CI endpoints)
qnorm(0.975)
# Simulate
set.seed(1); rnorm(5, mean=10, sd=2)

# Counts
lambda = 2
# P(X = 0) when X ~ Pois(λ=2)
dpois(0, lambda)
# P(X ≤ 3) when X ~ Pois(λ=2)
ppois(3, lambda=2)
# Upper tail: P(X ≥ 4) is P(X > 3) = 1 - P(X ≤ 3)
ppois(3, lambda, lower.tail = FALSE)
1-ppois(3, lambda)
# Interval probability: P(4 ≤ X ≤ 8)
ppois(8, lambda) - ppois(3, lambda)
# Quantiles
qpois(0.95, lambda)   # 95th percentile
# random draws
rpois(10, lambda)



# Plot distributions ----
# Gaussian
rnorm(1000, mean = 0, sd = 1) %>% hist(100)

# Multivariate gaussion
## For example, 2-D
set.seed(1)
mu <- c(0,0) # means of the variables
cov <- matrix(c(1, 0.8, 0.8, 1), nrow = 2) # covariance matrix. a postive-definite symmetric matrix
samples <- MASS::mvrnorm(1000, mu = mu, Sigma = cov*4) # variance is 1*4
plot(samples)
contour(MASS::kde2d(samples[,1], samples[,2]), add = T, col = "red") # density
var(samples[,1]) # should be close to 4
var(samples[,1])

## Cholesky decomposition
m <- t(chol(cov)) %*% matrix(rnorm(2*1000, 0, 2), nrow = 2) %>% t() # with sd = 2
plot(m)
contour(MASS::kde2d(m[,1], m[,2]), add = T, col = "red") # density
var(m[,1]) # should be close to 2^2
var(m[,2])

# Poisson: number of events that happens over a fixed period of time
rpois(1000, lambda = 10) %>% hist(100)

# Binomial: number of "successes" in a fixed number of tries
rbinom(1000, size = 10, prob = 0.5) %>% hist(100)

# Bernoulli or the response in logistic regression
p = 0.5
rbinom(1000, size = 1, prob = p) %>% hist(100)
# var = p * (1-p); mean = p
var(rbinom(1000, size = 1, prob = 0.5))
p * (1-p)

# Negative binomial: number of "failures" until you get a certain number of successes (k)
## Overdispersion: var > u
## It has extra variability controlled by the dispersion parameter.
## The dispersion parameter k (or size in the rnbinom) determines the number of successes you want
## probability p is the success probability at each toss
## mu is the expected number of failures before reaching k

## Small k, large overdispersion: var >> mean
rnbinom(1000, mu = 10, size = 2) %>% hist(100)

## Large k, var approximates mean, so becomes poisson
rnbinom(1000, mu = 10, size = 100) %>% hist(100)
rpois(1000, lambda = 10) %>% hist(100)

# mu = k (1-p)/p
# var = mu + mu^2/k
k = 100
p = 0.5
mu = k * (1-p) / p
rnbinom(10000, prob = p, size = k) %>% hist(100)
rnbinom(10000, mu = mu, size = k) %>% hist(100)


# zero-truncated, zero inflated, zero-altered (hurdle)
n = 10000
pi0 <- 0.2     # 20% structural zeros
lambda = 10

## Zero truncated
rztpois <- function(n, lambda) {
    u <- runif(n, min = .Machine$double.eps, max = 1)
    qpois(u * (1 - exp(-lambda)) + exp(-lambda), lambda)
}
rztpois(1000, lambda) %>% hist

## Zero-Inflated Poisson (ZIP): with prob pi_zero -> forced zero; else Poisson(lambda)
z <- rbinom(n, 1, pi0)               # 1 = structural zero
y <- rpois(n, lambda)                     # sampling zeros also possible
y[z == 1] <- 0L
hist(y)

## Hurdle Poisson: with prob pi_zero -> zero; else positive via zero-trunc Poisson
y <- integer(n) # all 0s
is_pos <- rbinom(n, 1, 1 - pi0)      # 1 = crosses the hurdle
k <- sum(is_pos == 1)
if (k > 0) y[is_pos == 1] <- rztpois(k, lambda)
hist(y)
table(y == 0) # roughly 20% is 0s
