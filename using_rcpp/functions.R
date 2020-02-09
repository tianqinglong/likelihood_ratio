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

# likelihood-ratio test
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

	if (y_n != 0 && y_n !=(n-r))
	{
		loglik <- sum( dweibull(ft, mles[1], mles[2], log = TRUE) )+
			(n-r)*pweibull(t_c, mles[1], mles[2], lower.tail = FALSE, log.p = TRUE)+
			y_n*log(y_n/(n-r))+(n-r-y_n)*log(1-y_n/(n-r))
	}
	else
	{
		loglik <- sum( dweibull(ft, mles[1], mles[2], log = TRUE) )+
			(n-r)*pweibull(t_c, mles[1], mles[2], lower.tail = FALSE, log.p = TRUE)
	}

	return (loglik);
}

eval_y <- function(y_n, t_w, mles, dat)
{
	mles_3param <- mle_solve_3_params_backup(dat, y_n, t_w)
	log_ratio <- -2*( gn_loglik(y_n, mles_3param, dat, t_w)-gd_loglik(y_n, dat, mles) )

	return(log_ratio)
}

find_mid <- function(p, dat, t_w)
{
  r <- dat[[1]]
  t_c <- dat[[2]]
  n <- dat[[4]]
  mles <- find_mle2_with_backup(dat)
  delta <- pweibull(t_w, mles[1], mles[2])-pweibull(t_c, mles[1], mles[2])
  eF <- n*delta
  qch <- qchisq(p, df = 1)
  
  if (eval_y(eF, t_w, mles, dat) < qch) {return(eF)}
  
  sample_points <- seq(from = 0, to = n-r, length.out = min(5, n-r))
  for (i in 1:length(sample_points[c(-1, -5)]))
  {
    x <- sample_points[i]
    if (eval_y(x, t_w, mles, dat) < qch) {return(x)}
  }

  sample_points <- seq(from = 0, to = n-r, length.out = min(20, n-r))
  for (i in 1:length(sample_points[c(-1, -20)]))
  {
  	x <- sample_points[i]
  	if (eval_y(x, t_w, mles, dat) < qch) {return(x)}
  }

  sample_points <- seq(from = 0, to = n-r, length.out = min(50, n-r))
  for (i in 1:length(sample_points[c(-1, -50)]))
  {
	  x <- sample_points[i]
	  if (eval_y(x, t_w, mles, dat) < qch) {return(x)}
  }

  sample_points <- seq(from = 0, to = n-r, length.out = min(100, n-r))
  for (i in 1:length(sample_points[c(-1, -100)]))
  {
	  x <- sample_points[i]
	  if (eval_y(x, t_w, mles, dat) < qch) {return(x)}
  }

  sample_points <- seq(from = 0, to = n-r, length.out = min(200, n-r))
  for (i in 1:length(sample_points[c(-1, -200)]))
  {
	  x <- sample_points[i]
	  if (eval_y(x, t_w, mles, dat) < qch) {return(x)}
  }

  sample_points <- seq(from = 0, to = n-r, length.out = min(500, n-r))
  for (i in 1:length(sample_points[c(-1, -500)]))
  {
	  x <- sample_points[i]
	  if (eval_y(x, t_w, mles, dat) < qch) {return(x)}
  }

  sample_points <- seq(from = 0, to = n-r, length.out = min(1000, n-r))
  for (i in 1:length(sample_points[c(-1, -1000)]))
  {
	  x <- sample_points[i]
	  if (eval_y(x, t_w, mles, dat) < qch) {return(x)}
  }

  stop("cannot find mid point!")
}

solve_discrete_root <- function(p, lp, sp, dat, mles, t_w)
{
  mid <- round((lp+sp)/2)
  
  while(abs(lp-sp) > 1)
  {
    y_val <- eval_y(mid, t_w, mles, dat)
    if(y_val > qchisq(p, df = 1))
    {
      lp <- mid
    }
    else
    {
      sp <- mid
    }
    
    mid <- round((lp+sp)/2)
  }
  
  pred <- ifelse(lp < sp, mid, mid+1)
  
  return(pred)
}

lik_ratio_pred <- function(p, dat, t_w)
{
  r <- dat[[1]]
  t_c <- dat[[2]]
  n <- dat[[4]]
  
  mles <- find_mle2_with_backup(dat)
  midpoint <- find_mid(p, dat, t_w)
  
  lwb <- solve_discrete_root(p, 0, midpoint, dat, mles, t_w)
  upb <- solve_discrete_root(p, n-r, midpoint, dat, mles, t_w)
  
  return(c(lwb, upb))
}

# bootstrap and gpq
bootstrap_sample <- function(mles, t_c, n)
{
  complete_data <- rweibull(n, mles[1], mles[2])
  index <- complete_data < t_c
  r <- sum(index)
  if (r < 2)
  {
    return( bootstrap_sample(mles, t_c, n) )
  }
  ft <- complete_data[index]
  return(list(Number_of_Failures = r, Censor_Time = t_c, Failure_Times = ft, Total_Number = n))
}

generate_bootstrap_draws <- function(dat,B = 1000)
{
  mles <- find_mle2_with_backup(dat)

  list_mle_r <- lapply(1:B, function(x) {
    bsample <- bootstrap_sample(mles, dat[[2]], dat[[4]])
    return(list(MLEs = find_mle2_with_backup(bsample), R = bsample[[1]]))
  })
  
  return(list_mle_r)
}

get_p_star <- function(list_mle_r, t_w, t_c)
{
  lapply(list_mle_r, function(x) {
    compute_p(t_c, t_w, x$MLEs[1], x$MLEs[2])
  }) -> p_array
  
  return(unlist(p_array))
}

get_p_starstar <- function(list_mle_r, mles, t_w, t_c)
{
  lapply(list_mle_r, function(x) {
    compute_p_starstar(mles[1], mles[2], x$MLEs[1], x$MLEs[2], t_c, t_w)
  }) -> pp_array
  
  return(unlist(pp_array))
}

boot_solve_discrete <- function(p, p_array, n)
{
  lp <- 0
  up <- n
  
  if(pred_dist(lp, p_array, n) >= p)
  {
    return(lp)
  }
  if(pred_dist(up, p_array, n) <= p)
  {
    return(up)
  }
  
  mid <- (lp+up)/2
  while(abs(lp-up)>1)
  {
    pmid <- pred_dist(mid, p_array, n)
    if (pmid < p)
    {
      lp <- mid
    }
    else if(pmid > p)
    {
      up <- mid
    }
    else if(abs(pmid-p) <= 1e-5)
    {
      return(mid)
    }
    
    mid <- round((lp+up)/2)
  }
  
  pred <- ifelse(p>0.5, lp, max(lp-1, 0))
  
  return(lp)
}

bootstrap_and_GPQ_and_calibration_pred <- function(p, dat, t_w)
{
  t_c <- dat$Censor_Time
  mles <- find_mle2_with_backup(dat)
  list_mle_r <- generate_bootstrap_draws(dat)
  
  p_ast <- get_p_star(list_mle_r, t_w, t_c)
  p_astast <- get_p_starstar(list_mle_r, mles, t_w, t_c)
  
  bp <- boot_solve_discrete(p, p_ast, dat[[4]]-dat[[1]])
  gp <- boot_solve_discrete(p, p_astast, dat[[4]]-dat[[1]])
  
  return(c(bp, gp))
}
