#

library(tidyverse)

# Gaussian
rnorm(1000, mean = 0, sd = 1) %>% hist(100)

# Poisson: number of events that happens over a fixed period of time
rpois(1000, lambda = 10) %>% hist(100)

# Binomial: number of "successes" in a fixed number of tries
rbinom(1000, size = 10, prob = 0.5) %>% hist(100)


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

