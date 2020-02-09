# source("functions.R")

r <- 30
pf1 <- 0.1
beta <- 2
eta <- 1

pf2 <- 0.2
t_w <- qweibull(pf2, beta, eta)

# simulate dataset
dat <- generate_censored_data(r, pf1, beta, eta)

# find MLEs
(mles <- find_mle2_with_backup(dat))

# find the 95% and 90% prediction interval
lik_ratio_pred(0.9, dat, t_w)
lik_ratio_pred(0.8, dat, t_w)
