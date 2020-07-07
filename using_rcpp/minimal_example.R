source("functions.R")
library(parallel)

beta <- 2
pf1 <- 0.1
delta <- 0.1
r <- 10
eta <- 1

t_w <- qweibull(pf1+delta, beta, eta)

dat <- generate_censored_data(r, pf1, beta, eta)
prediction_main_paper(dat, t_w, beta, eta, B = 5000)