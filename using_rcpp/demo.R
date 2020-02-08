source("functions.R")

r <- 10
pf1 <- 0.1
beta <- 2
eta <- 1

pf2 <- 0.2
t_w <- qweibull(pf2, beta, eta)

dat <- generate_censored_data(r, pf1, beta, eta)

y_n <- 20
eval_y(y_n, t_w, mles, dat)