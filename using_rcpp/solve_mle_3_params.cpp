#include <RcppGSL.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

extern "C" {
	#include <math.h>
	#include <gsl/gsl_vector.h>
	#include <gsl/gsl_multiroots.h>
	#include "df_fn_2.h"
	#include "df_fn_2.c"
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
Rcpp::NumericVector mle_solve_3param(Rcpp::List dat, int y_n, double t_w)
{
	int r = dat[0], n = dat[3];
	double t_c = dat[1];
	Rcpp::NumericVector ft = dat[2];
	double ft_log[r];

	int i;
	for (i = 0; i < r; i++)
	{
		ft_log[i] = log(ft[i]);
	}

	struct rparams p = {y_n, r, n, log(t_c), log(t_w), ft_log};

	const int np = 2;
	int status, iter = 0;
	
	const gsl_multiroot_fsolver_type *T;
	gsl_multiroot_fsolver *s;

	gsl_multiroot_function f = {&loglik_f, np, &p};
	double gev_init[2] = {0, log(0.25)};
	gsl_vector *x = gsl_vector_alloc(np);

	gsl_vector_set(x, 0, gev_init[0]);
	gsl_vector_set(x, 1, gev_init[1]);

	T = gsl_multiroot_fsolver_hybrids;
	s = gsl_multiroot_fsolver_alloc (T, np);
	gsl_multiroot_fsolver_set (s, &f, x);

	do
	{
		iter++;
		status = gsl_multiroot_fsolver_iterate (s);

		// print_state (iter, s);

		if (status)
			break;

		status = gsl_multiroot_test_residual (s->f, 1e-7);
	}
	while (status == GSL_CONTINUE && gsl_strerror(status));

	Rcpp::NumericVector MLEs (2);
	MLEs[0] = 1/exp(gsl_vector_get (s->x, 1));
	MLEs[1] = exp( gsl_vector_get (s->x, 0) );

	printf ("status = %s\n", gsl_strerror (status));
	gsl_multiroot_fsolver_free (s);
	gsl_vector_free (x);

	return MLEs;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
