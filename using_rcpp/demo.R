source("functions.R")

Sys.time() -> t0
Er <- 50
pf1 <- 0.01
beta <- 2
eta <- 1
pf2 <- 0.2

dat <- generate_censored_data(Er, pf1, beta, eta)
t_c <- qweibull(pf1, beta, eta)
t_w <- qweibull(pf2, beta, eta)

prediction_six_methods(dat, t_w, beta, eta)

Sys.time() - t0