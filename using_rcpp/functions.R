Rcpp::sourceCpp("simulate_data.cpp")
Rcpp::sourceCpp("solve_mle_root.cpp")
Rcpp::sourceCpp("solve_mle_optim_backup.cpp")
Rcpp::sourceCpp("solve_mle_3_params.cpp")
Rcpp::sourceCpp("solve_mle_3_params_optim_backup.cpp")

find_mle2_with_backup <- function(dat)
{
  mles <- mle_solve_root(dat)
  
  # if the root finding routine breaks down use the optimization
  if(mles[1] < 0)
  {
    mles <- mle_solve_backup_2(dat)
  }
  
  return (mles)
}

gn_loglik <- function(y_n, mles3, dat, t_w)
{
	r <- dat[[1]]
	t_c <- dat[[2]]
	ft <- dat[[3]]
	n <- dat[[4]]

	loglik <- sum( dweibull(ft, mles3[1], mles3[2], log = TRUE) )+
		y_n*log(pweibull(t_w, mles3[1], mles3[2]) - pweibull(t_c, mles3[1], mles3[2]))+
		(n-y_n-r)*pweibull(t_w, mles3[1], mles3[2], lower.tail = FALSE, log.p = TRUE)

	return (loglik)
}

gd_loglik <- function(y_n, dat, mles)
{
	r <- dat[[1]]
	t_c <- dat[[2]]
	ft <- dat[[3]]
	n <- dat[[4]]

	loglik <- sum( dweibull(ft, mles[1], mles[2], log = TRUE) )+
		(n-r)*pweibull(t_c, mles[1], mles[2], lower.tail = FALSE, log.p = TRUE)+
		y_n*log(y_n/(n-r))+(n-r-y_n)*log(1-y_n/(n-r))

	return (loglik);
}

eval_y <- function(y_n, t_w, mles, dat)
{
	mles_3param <- mle_solve_3param(dat, y_n, t_w)
	log_ratio <- -2*( gn_loglik(y_n, mles_3param, dat, t_w)-gd_loglik(y_n, dat, mles) )

	return(log_ratio)
}