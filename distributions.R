#

library(tidyverse)

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
