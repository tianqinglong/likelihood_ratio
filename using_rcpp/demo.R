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

dat <- generate_censored_data(5, 0.1, 2, 1)
r <- dat[[1]]
n <- dat[[4]]
t_c <- qweibull(0.1, 2, 1)
t_w <- qweibull(0.2, 2, 1)
mles <- find_mle2_with_backup(dat)
list_mles_r <- generate_bootstrap_draws(dat)
p_ast <- get_p_star(list_mles_r, t_w, t_c)
p_astast <- get_p_starstar(list_mles_r, mles, t_w, t_c)
boot_solve_discrete(0.05, p_ast, n-r)