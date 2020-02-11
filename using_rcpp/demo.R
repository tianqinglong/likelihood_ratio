# source("functions.R")

r <- 30
pf1 <- 0.1
beta <- 2
eta <- 1
n <- r/pf1

pf2 <- 0.2
t_c <- qweibull(pf1, beta, eta)
t_w <- qweibull(pf2, beta, eta)

# simulate dataset
dat <- generate_censored_data(r, pf1, beta, eta)
mles <- find_mle2_with_backup(dat)

# find MLEs
B <- 300
num_y <- 20
list_bootstrap_samples <- lapply(1:B, function(x) {bootstrap_sample(mles, t_c, n)})

# for a single bootstrap sample
bs <- list_bootstrap_samples[[1]]
r_star <- bs$Number_of_Failures
y_array <- generate_y_n(r_star, n, t_c, t_w, mles[1], mles[2], num_y)
eval_y_array <- double(length = num_y)
for (i in 1:num_y)
{
  mles_star <- find_mle2_with_backup(bs)
  eval_y_array[i] <- eval_y(y_array[i], t_w, mles_star, bs)
}

#
ratio_emp <- generate_ratio_array(mles, t_c, t_w, n, num_per_sample = 30, B = 3000)
qt <- quantile(ratio_emp, probs = c(0.8,0.9))
qch <- qt[1];
mid <- find_mid_boot(qch, dat, t_w)
solve_discrete_root_boot(qch, 0, mid, dat, mles, t_w)

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

t0 <- Sys.time()
Er <- 5
pf1 <- 0.01
beta <- 4
eta <- 1
pf2 <- 0.05
dat <- generate_censored_data(Er, pf1, beta, eta)
t_w <- qweibull(pf2, beta, eta)
prediction_five_methods(dat, t_w, 2, 1)
Sys.time() - t0