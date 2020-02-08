#include <RcppGSL.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

extern "C" {
	#include<math.h>
	#include "df_fn.h"
	#include "df_fn.c"
}

// [[Rcpp::depends(RcppGSL)]]

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
Rcpp::NumericVector mle_solve_root(Rcpp::List dat)
{
  int r = dat[0], n = dat[3];
  double t_c = dat[1];
  Rcpp::NumericVector ft = dat[2];
  double *fp = ft.begin();
  
  int status;
  int iter = 0, max_iter = 500;
  Rcpp::NumericVector MLEs (2);

  gsl_set_error_handler_off();

  double sol = 0;
  double x_lo = 0.01, x_hi = 25;

  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  gsl_function F;
  
  struct weibull_data_params dat_param = {r, n, t_c, fp};
  F.function = &deriv;
  F.params = &dat_param;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc(T);
  gsl_root_fsolver_set(s, &F, x_lo, x_hi);

  do
  {
  	iter++;
  	status = gsl_root_fsolver_iterate(s);

  	if(status) {break;}

  	sol = gsl_root_fsolver_root(s);
  	x_lo = gsl_root_fsolver_x_lower(s);
  	x_hi = gsl_root_fsolver_x_upper(s);

  	status = gsl_root_test_interval(x_lo, x_hi, 0, 0.001);
  }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fsolver_free(s);

  if(status == GSL_SUCCESS)
  {
  	double lambda = 0;

  	for(int i=0; i<r; i++)
  	{
  		lambda += pow(ft[i], sol);
  	}

  	lambda += (n-r)*pow(t_c, sol);
  	lambda = lambda/r;
  	lambda = pow(lambda, 1/sol);

  	MLEs[0] = sol;
  	MLEs[1] = lambda;
  }
  else
  {
  	MLEs[0] = -1;
  	MLEs[1] = -1;
  }

  return MLEs;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
