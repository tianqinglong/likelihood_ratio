#include <RcppGSL.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>

extern "C" {
	#include <math.h>
	#include <gsl/gsl_vector.h>
	#include <gsl/gsl_multimin.h>
	#include "df_fn_4.h"
	#include "df_fn_4.c"
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

Rcpp::NumericVector mle_solve_3_params_backup(Rcpp::List dat, int y_n, double t_w)
{
	int r = dat[0], n = dat[3];
	double t_c = dat[1];
	Rcpp::NumericVector ft = dat[2];
	double *fp = ft.begin();

	struct rparams par = {y_n, r, n, t_c, t_w, fp};

	const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer *s = NULL;
	gsl_vector *ss, *x;
	gsl_multimin_function minex_func;

	int iter = 0;
	int status;
	double size;

	x = gsl_vector_alloc (2);

	gsl_vector_set (x, 0, log(2));
	gsl_vector_set (x, 1, 0);

	/* Set initial step sizes to 1 */
	ss = gsl_vector_alloc (2);
	gsl_vector_set_all (ss, 1);

	/* Initialize method and iterate */
	minex_func.n = 2;
	minex_func.f = my_f;
	minex_func.params = &par;

	s = gsl_multimin_fminimizer_alloc (T, 2);
	gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

	do
    {
    	iter++;
    	status = gsl_multimin_fminimizer_iterate(s);

    	if (status)
    		break;

    	size = gsl_multimin_fminimizer_size (s);
    	status = gsl_multimin_test_size (size, 1e-2);
    }
    while (status == GSL_CONTINUE && iter < 100);

    Rcpp::NumericVector MLEs (2);
    if (status == GSL_SUCCESS)
    {
    	MLEs[0] = exp(gsl_vector_get (s->x, 0));
    	MLEs[1] = exp(gsl_vector_get (s->x, 1));
    }
    else
    {
    	MLEs[0] = -1;
    	MLEs[1] = -1;
    }

    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);

    return MLEs;
}