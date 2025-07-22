#' Compare lm results from lm and brms

library(tidyverse)
library(brms)

# Example data
set.seed(123)
X <- 1:10
Y <- 2 + 1.9 * X + rnorm(10, 0, 1)  # True slope=1.9, intercept=2, noise=1
df <- data.frame(X=X, Y=Y)

# Frequentist ----
mod <- lm(Y~X, data = df)
su <- summary(mod)
su

Y-predict(mod) # residuals
residuals(mod)

sst <- sum((Y-mean(Y))^2)              # SST: sum of squares total
ssr <- sum((predict(mod) - mean(Y))^2) # SSR: sum of squares regression
sse <- sum((Y-predict(mod))^2)         # SSE: sum of squares error
n = length(X) # number of samples
k = 2         # number of parameters
df1 = 1      # n of predictors
df2 = n-k    # residual degree of freedom:  n of samples - n of parameters
msr = ssr / df1     # MSR: mean square for regression
mse = sse / (n-2)   # MSE: mean square for error df for residual is n of samples - n of predicators - 1

# R squared = SSR/SST or 1-SSR-SST
ssr/sst
1 - sse/sst
su$r.squared

# adjusted R squared penalized by n of parameters
1 - (sse/(n-k)) / (sst/(n-1))
su$adj.r.squared

# F statistics
msr/mse
su$fstatistic[1]

# F test
pf(msr/mse, df1, df2, lower.tail = F)
su$coefficients["X", "Pr(>|t|)"]

# t test for the coefficient
beta_hat <- su$coefficients["X", "Estimate"]
se_beta <- su$coefficients["X", "Std. Error"]
t_value = beta_hat/se_beta # t = beta_hat/se_beta
t_value
su$coefficients["X", "t value"]
2 * (1-pt(abs(t_value), df2))
su$coefficients["X", "Pr(>|t|)"]

# RSE: Residual standard error. Or sigma
sqrt(sse/(n-k))
su$sigma


# Bayesian ----
mod2 <- brm(Y ~ X, data=df, chains=2, cores=2,iter=2000,seed=123)

# Fixed effect estimates
fixef(mod2)
plot(mod2) # ok it shows the average across chains and only uses iterations after burnins

# Estiamte distribution
as_draws_df(mod2) %>%
    as_tibble %>%
    #filter(.chain == 2) %>%
    ggplot() +
    geom_histogram(aes(x = sigma), binwidth = 0.1, color = 1, fill = "white") +
    theme_classic()


# R2
fitted_samples <- fitted(mod2, summary = F)
dim(fitted_samples)

n_samples <- ndraws(mod2)
R2_samples <- numeric(n_samples)

for (s in 1:n_samples) {
    Y_pred_s <- fitted_samples[s,]
    sst <- sum((Y - mean(Y))^2)
    sse <- sum((Y - Y_pred_s)^2)
    R2_samples[s] <- 1-sse/sst
}

median(R2_samples)
sd(R2_samples)
quantile(R2_samples, probs = c(0.025, 0.5, 0.975))
bayes_R2(mod2)
