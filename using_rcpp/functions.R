Rcpp::sourceCpp("simulate_data.cpp")
Rcpp::sourceCpp("solve_mle_root.cpp")
Rcpp::sourceCpp("solve_mle_optim_backup.cpp")
Rcpp::sourceCpp("solve_mle_3_params.cpp")
Rcpp::sourceCpp("solve_mle_3_params_optim_backup.cpp")

# compute_p <- function(t_c, t_w, beta, eta){
#   cond_p <- (pweibull(t_w, beta, eta)-pweibull(t_c, beta, eta))/(1-pweibull(t_c, beta, eta))
#   return(cond_p)
# }

find_mle2_with_backup <- function(dat)
{
  mles <- mle_solve_root(dat)
  
  # if the root finding routine breaks down use the optimization
  if(mles[1] < 1e-6 || mles[2] < 1e-6 || mles[1] > 20 || mles[2] > 20)
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

gd_loglik <- function(y_n, dat)
{
	r <- dat[[1]]
	t_c <- dat[[2]]
	ft <- dat[[3]]
	n <- dat[[4]]

	mles <- find_mle2_with_backup(dat)

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
	log_ratio <- -2*( gn_loglik(y_n, mles_3param, dat, t_w)-gd_loglik(y_n, dat) )

	return(log_ratio)
}

find_mid <- function(qch, dat, t_w)
{
  r <- dat[[1]]
  t_c <- dat[[2]]
  n <- dat[[4]]
  mles <- find_mle2_with_backup(dat)
  delta <- pweibull(t_w, mles[1], mles[2])-pweibull(t_c, mles[1], mles[2])
  eF <- round(n*delta); eF <- max(1, eF); eF <- min(n-r-1, eF)
  
  if (eval_y(eF, t_w, mles, dat) < qch) {return(eF)}
  
  # first backup
  sample_points <- round(seq(from = 0, to = n-r, length.out = min(30, n-r)))
  for (i in 1:length(sample_points))
  {
    x <- sample_points[i]
    if (eval_y(x, t_w, mles, dat) < qch) {return(x)}
  }
  
  # last backup
  for (x in 0:(n-r))
  {
  	if (eval_y(x, t_w, mles, dat) < qch) {return(x)}
  }

  stop("cannot find mid point!")
}

solve_discrete_root <- function(qch, lp, sp, dat, mles, t_w)
{
  mid <- ceiling((lp+sp)/2)
  
  while(mid!=lp & mid!=sp)
  {
    y_val <- eval_y(mid, t_w, mles, dat)
    if(y_val > qch)
    {
      lp <- mid
    }
    else if(y_val < qch)
    {
      sp <- mid
    }
    
    mid <- ceiling((lp+sp)/2)
  }

  return(lp)
}

lik_ratio_pred <- function(p, dat, t_w)
{
  r <- dat[[1]]
  t_c <- dat[[2]]
  n <- dat[[4]]
  qch <- qchisq(p, df = 1)
  
  mles <- find_mle2_with_backup(dat)
  midpoint <- tryCatch(find_mid(qch, dat, t_w), error = function(e) {-1})

  if (midpoint < 0)
  {
  	return(c(NA, NA))
  }
  
  if (eval_y(0, t_w, mles, dat) > qchisq(p, df = 1))
  {
  	lwb <- solve_discrete_root(qch, 0, midpoint, dat, mles, t_w)
  }
  else
  {
  	lwb <- 0
  }
  
  if (eval_y(n-r, t_w, mles, dat) > qchisq(p, df = 1))
  {
  	upb <- solve_discrete_root(qch, n-r, midpoint, dat, mles, t_w)
  }
  else
  {
  	upb <- n-r
  }
  
  return(c(lwb, upb))
}

generate_bootstrap_draws <- function(dat,B = 5000)
{
  mles <- find_mle2_with_backup(dat)

  iter <- 1
  output <- NULL
  while(iter <= B)
  {
  	bsample <- bootstrap_sample(mles, dat[[2]], dat[[4]])
  	MLE_raw <- find_mle2_with_backup(bsample)

  	if (MLE_raw[1] < 0)
  	{
  		next
  	}

  	output[[iter]] <- list(MLEs = MLE_raw, R = bsample[[1]], BData = bsample)
  	iter <- iter + 1
  }
  
  return(output)
}

get_p_star <- function(list_mle_r, t_w, t_c)
{
  p_array <- numeric(length(list_mle_r))
  for(i in 1:length(list_mle_r))
  {
    mle_r <- list_mle_r[[i]]
    MLE <- mle_r[[1]]
    
    p_array[i] <- compute_p(t_c, t_w, MLE[1], MLE[2])
  }
  
  p_array[p_array > 1] <- 1
  return(p_array)
}

get_p_starstar <- function(list_mle_r, MLEs, t_w, t_c)
{
  p_array <- numeric(length(list_mle_r))
  for(i in 1:length(list_mle_r))
  {
    mle_r <- list_mle_r[[i]]
    BT_MLEs <- mle_r[[1]]
    
    MLE_mu <- log (MLEs[[2]])
    MLE_sigma <- 1/MLEs[[1]]
    
    BT_mu <- log (BT_MLEs[[2]])
    BT_sigma <- 1/BT_MLEs[[1]]
    
    GPQ_Beta <- 1 / ( MLE_sigma*MLE_sigma / BT_sigma )
    GPQ_Eta <- exp( MLE_mu + (MLE_mu - BT_mu)/BT_sigma*MLE_sigma )
    
    p_array[i] <- compute_p(t_c, t_w, GPQ_Beta, GPQ_Eta)
  }
  
  p_array[p_array > 1] <- 1
  return(p_array)
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
  
  mid <- ceiling((lp+up)/2)
  while(mid!=lp & mid!=up)
  {
    pmid <- pred_dist(mid, p_array, n)
    if (pmid < p)
    {
      lp <- mid
    }
    else
    {
      up <- mid
    }
    mid <- ceiling((lp+up)/2)
  }
  pred <- ifelse(p>0.5, up, lp)
  
  return(pred)
}

# likelihood ratio + bootstrap
generate_ratio_array <- function(mles, t_c, t_w, n, list_mles_r, num_per_sample)
{
	sapply(list_mles_r, function(x){
		x <- x$BData
		r_star <- x$Number_of_Failures
		y_array <- generate_y_n(r_star, n, t_c, t_w, mles[1], mles[2], num_per_sample)
		mles_star <- find_mle2_with_backup(x)

		eval_y_array <- double(length = num_per_sample)
		for (i in 1:num_per_sample)
		{
			eval_y_array[i] <- eval_y(y_array[i], t_w, mles_star, x)
		}
		return(eval_y_array)
	}
	) -> ratio_emp

  return(as.vector(ratio_emp))
}

lik_ratio_pred_boot <- function(dat, t_w, list_mles_r, num_of_samples = 50)
# default 90% 95% prediction bonuds
{
  r <- dat[[1]]
  t_c <- dat[[2]]
  n <- dat[[4]]

  mles <- find_mle2_with_backup(dat)
  ratio_emp <- generate_ratio_array(mles, t_c, t_w, n, list_mles_r, num_of_samples)

  qt <- quantile(ratio_emp, probs = c(0.8, 0.9), na.rm = T)
  qch <- qt[1]
  mid <- tryCatch(find_mid(qch, dat, t_w), error = function(e) {-1})
  if (mid < 0)
  {
	L90 <- NA
	U90 <- NA
  }
  else
  {
	L90 <- solve_discrete_root(qch, 0, mid, dat, mles, t_w)
	U90 <- solve_discrete_root(qch, n-r, mid, dat, mles, t_w)
  } 
  qch <- qt[2]
  mid <- tryCatch(find_mid(qch, dat, t_w), error = function(e) {-1})
  if (mid < 0)
  {
	L95 <- NA
	U95 <- NA
  }
  else
  {
	L95 <- solve_discrete_root(qch, 0, mid, dat, mles, t_w)
	U95 <- solve_discrete_root(qch, n-r, mid, dat, mles, t_w)
  }

  return(list(bounds = c(L95, L90, U90, U95), emean = mean(ratio_emp)))
}

# Bartlett Corrections

lik_ratio_pred_bartlett_correct <- function(p, dat, t_w, emean)
{
  r <- dat[[1]]
  n <- dat[[4]]
  qch <- emean * qchisq(p, df = 1)
  mid <- tryCatch(find_mid(qch, dat, t_w), error = function(e) {-1})
  if (mid < 0)
  {
	LB <- NA
	UB <- NA
  }
  else
  {
	LB <- solve_discrete_root(qch, 0, mid, dat, mles, t_w)
	UB <- solve_discrete_root(qch, n-r, mid, dat, mles, t_w)
  }

  return(c(LB, UB))
}

# calibration method
pred_root_empirical <- function(list_mles_r, mles, t_c, t_w, n, num_per_sample = 30)
{
	phat <- compute_p(t_c, t_w, mles[1], mles[2])
	sapply(list_mles_r, function(x) {
		p_star <- compute_p(t_c, t_w, x$MLEs[1], x$MLEs[2])
		p_star <- ifelse(p_star > 1, 1, p_star)
		p_star <- ifelse(p_star < 0, 0, p_star)

		ystar <- rbinom(num_per_sample, n-x$R, phat)
		u_array <- pbinom(ystar, n-x$R, p_star)

		return(u_array)
	}
	) -> u_emp

	u_emp <- as.vector(u_emp)
	four_quantile <- quantile(u_emp, probs = c(0.05, 0.1, 0.9, 0.95), na.rm = T)

	return(four_quantile)
}

compute_cp <- function(pb, t_c, t_w, beta, eta, n, r)
# four bounds
{
  p <- compute_p(t_c, t_w, beta, eta)
  cp1 <- pbinom(pb[1], n-r, p, lower.tail = F)+dbinom(pb[1], n-r, p)
  cp2 <- pbinom(pb[2], n-r, p, lower.tail = F)+dbinom(pb[2], n-r, p)
  cp3 <- pbinom(pb[3], n-r, p)
  cp4 <- pbinom(pb[4], n-r, p)
  
  return(c(cp1, cp2, cp3, cp4))
}

pb2cp <- function(pb_mat, t_c, t_w, beta, eta, n, r)
{
  cp_mat <- matrix(nrow = nrow(pb_mat), ncol = 4)
  rownames(cp_mat) <- rownames(pb_mat)
  colnames(cp_mat) <- colnames(pb_mat)
  
  for(i in 1:nrow(pb_mat))
  {
    cp_mat[i,] <- compute_cp(pb_mat[i,], t_c, t_w, beta, eta, n, r)
  }
  
  return(cp_mat)
}

# predictions using four methods
prediction_four_methods <- function(dat, t_w, beta, eta)
{
  n <- dat$Total_Number
  r <- dat$Number_of_Failures
  t_c <- dat$Censor_Time
  mles <- find_mle2_with_backup(dat)
  list_mles_r <- generate_bootstrap_draws(dat)
  p_ast <- get_p_star(list_mles_r, t_w, t_c)
  p_astast <- get_p_starstar(list_mles_r, mles, t_w, t_c)
  
  pb_mat <- matrix(nrow = 4, ncol = 4)
  
  rownames(pb_mat) <- c("Ratio", "Bootstrap", "GPQ", "Calibration")
  colnames(pb_mat) <- c("Lower95", "Lower90", "Upper90", "Upper95")

  # likelihood ratio based prediction
  LU95 <- lik_ratio_pred(0.9, dat, t_w)
  LU90 <- lik_ratio_pred(0.8, dat, t_w)
  pb_mat[1, c(1, 4)] <- LU95
  pb_mat[1, c(2, 3)] <- LU90
  
  # bootstrap
  L95 <- boot_solve_discrete(0.05, p_ast, n-r)
  L90 <- boot_solve_discrete(0.1, p_ast, n-r)
  U90 <- boot_solve_discrete(0.9, p_ast, n-r)
  U95 <- boot_solve_discrete(0.95, p_ast, n-r)
  pb_mat[2,] <- c(L95, L90, U90, U95)

  # gpq
  L95 <- boot_solve_discrete(0.05, p_astast, n-r)
  L90 <- boot_solve_discrete(0.1, p_astast, n-r)
  U90 <- boot_solve_discrete(0.9, p_astast, n-r)
  U95 <- boot_solve_discrete(0.95, p_astast, n-r)
  pb_mat[3,] <- c(L95, L90, U90, U95)

  # calibration
  alpha_cali <- pred_root_empirical(list_mles_r, mles, t_c, t_w, n)
  cap <- qbinom(alpha_cali, n-r, compute_p(t_c, t_w, mles[1], mles[2]))
  cap[1] <- max(0, cap[1]-1)
  cap[2] <- max(0, cap[2]-1)
  pb_mat[4,] <- cap
  
  cp_mat <- pb2cp(pb_mat, t_c, t_w, beta, eta, n, r)
  
  return(list(Prediction_Bounds = pb_mat, Coverage_Probability = cp_mat))
}

# to fix the supplementary of main paper
prediction_main_paper <- function(dat, t_w, beta, eta, B = 5000)
{
  n <- dat$Total_Number
  r <- dat$Number_of_Failures
  t_c <- dat$Censor_Time
  mles <- find_mle2_with_backup(dat)
  list_mles_r <- generate_bootstrap_draws(dat, B)
  p_ast <- get_p_star(list_mles_r, t_w, t_c)
  p_astast <- get_p_starstar(list_mles_r, mles, t_w, t_c)
  
  pb_mat <- matrix(nrow = 4, ncol = 4)
  rownames(pb_mat) <- c("Bootstrap", "GPQ", "Calibration", "Plug-in")
  colnames(pb_mat) <- c("Lower95", "Lower90", "Upper90", "Upper95")
  
  # bootstrap
  L95 <- boot_solve_discrete(0.05, p_ast, n-r)
  L90 <- boot_solve_discrete(0.1, p_ast, n-r)
  U90 <- boot_solve_discrete(0.9, p_ast, n-r)
  U95 <- boot_solve_discrete(0.95, p_ast, n-r)
  pb_mat[1,] <- c(L95, L90, U90, U95)
  
  # gpq
  L95 <- boot_solve_discrete(0.05, p_astast, n-r)
  L90 <- boot_solve_discrete(0.1, p_astast, n-r)
  U90 <- boot_solve_discrete(0.9, p_astast, n-r)
  U95 <- boot_solve_discrete(0.95, p_astast, n-r)
  pb_mat[2,] <- c(L95, L90, U90, U95)
  
  # calibration
  alpha_cali <- pred_root_empirical(list_mles_r, mles, t_c, t_w, n)
  cap <- qbinom(alpha_cali, n-r, compute_p(t_c, t_w, mles[1], mles[2]))
  cap[1] <- max(0, cap[1]-1)
  cap[2] <- max(0, cap[2]-1)
  pb_mat[3,] <- cap
  
  # plug-in
  p_hat <- compute_p(t_c, t_w, mles[1], mles[2])
  L95 <- qbinom(0.05, n-r, p_hat)
  L90 <- qbinom(0.1, n-r, p_hat)
  U90 <- qbinom(0.9, n-r, p_hat)
  U95 <- qbinom(0.95, n-r, p_hat)
  pb_mat[4,] <- c(L95, L90, U90, U95)
  
  cp_mat <- pb2cp(pb_mat, t_c, t_w, beta, eta, n, r)
  
  return(list(Prediction_Bounds = pb_mat, Coverage_Probability = cp_mat))
}

# add the calibrated likelihood method
prediction_six_methods <- function(dat, t_w, beta, eta, B = 5000)
{
  n <- dat$Total_Number
  r <- dat$Number_of_Failures
  t_c <- dat$Censor_Time
  mles <- find_mle2_with_backup(dat)
  list_mles_r <- generate_bootstrap_draws(dat, B)
  p_ast <- get_p_star(list_mles_r, t_w, t_c)
  p_astast <- get_p_starstar(list_mles_r, mles, t_w, t_c)
  
  pb_mat <- matrix(nrow = 6, ncol = 4)
  
  rownames(pb_mat) <- c("LRT", "Bootstrap", "GPQ", "Calibration", "C-LRT", "BC-LRT")
  colnames(pb_mat) <- c("Lower95", "Lower90", "Upper90", "Upper95")

  # likelihood ratio based prediction
  pb_mat[1, c(1, 4)] <- lik_ratio_pred(0.9, dat, t_w)
  pb_mat[1, c(2, 3)] <- lik_ratio_pred(0.8, dat, t_w)
  
  # bootstrap
  L95 <- boot_solve_discrete(0.05, p_ast, n-r)
  L90 <- boot_solve_discrete(0.1, p_ast, n-r)
  U90 <- boot_solve_discrete(0.9, p_ast, n-r)
  U95 <- boot_solve_discrete(0.95, p_ast, n-r)
  pb_mat[2,] <- c(L95, L90, U90, U95)

  # gpq
  L95 <- boot_solve_discrete(0.05, p_astast, n-r)
  L90 <- boot_solve_discrete(0.1, p_astast, n-r)
  U90 <- boot_solve_discrete(0.9, p_astast, n-r)
  U95 <- boot_solve_discrete(0.95, p_astast, n-r)
  pb_mat[3,] <- c(L95, L90, U90, U95)

  # calibration
  alpha_cali <- pred_root_empirical(list_mles_r, mles, t_c, t_w, n)
  cap <- qbinom(alpha_cali, n-r, compute_p(t_c, t_w, mles[1], mles[2]))
  cap[1] <- max(0, cap[1]-1)
  cap[2] <- max(0, cap[2]-1)
  pb_mat[4,] <- cap
  
  # calibrated-likelihood ratio
  clr <- lik_ratio_pred_boot(dat, t_w, list_mles_r)
  pb_mat[5,] <- clr$bounds

  # bartlett_correction
  emean <- clr$emean
  pb_mat[6, c(1, 4)] <- lik_ratio_pred_bartlett_correct(0.9, dat, t_w, emean)
  pb_mat[6, c(2, 3)] <- lik_ratio_pred_bartlett_correct(0.8, dat, t_w, emean)

  cp_mat <- pb2cp(pb_mat, t_c, t_w, beta, eta, n, r)
  
  return(list(Prediction_Bounds = pb_mat, Coverage_Probability = cp_mat, Data = dat))
}
