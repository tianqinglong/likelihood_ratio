source("functions.R")

Er <- 15
pf1 <- 0.1
beta <- 4
eta <- 1
pf2 <- 0.2

dat <- generate_censored_data(Er, pf1, beta, eta)
t_c <- qweibull(pf1, beta, eta)
t_w <- qweibull(pf2, beta, eta)
( mles <- find_mle2_with_backup(dat) )

# list_mles_r <- generate_bootstrap_draws(dat, 5000)
# get_p_star(list_mles_r, t_w, t_c)

prediction_five_methods(dat, t_w, beta, eta)
