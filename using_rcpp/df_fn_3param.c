#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

struct rparams
{
	int y_n;
	int r;
	int n;
	double t_c;
	double t_w;
	double *ft;
};

int
loglik_f (const gsl_vector * x, void * params, gsl_vector * f)
{
	int r = ((struct rparams *) params)->r;
	int n = ((struct rparams *) params)->n;
	int y_n = ((struct rparams *) params)->y_n;

	double t_c = ((struct rparams *) params)->t_c;
	double t_w = ((struct rparams *) params)->t_w;
	double *ft = ((struct rparams *) params)->ft;

	const double u = gsl_vector_get(x, 0);
	const double b = gsl_vector_get(x, 1);

	int i;
	double dfdu = 0, dfdb = 0;
	for(i=0; i<r; i++)
	{
		dfdu += exp((ft[i]-u)/b)/b-1/b;
		dfdb += ( -b*exp(exp((ft[i] - u)/b) - (ft[i] - u)/b)*(exp((ft[i] - u)/b - exp((ft[i] - u)/b))/pow(b,2) + (exp((ft[i] - u)/b - exp((ft[i] - u)/b))*((ft[i] - u)/pow(b,2) - (exp((ft[i] - u)/b)*(ft[i] - u))/pow(b,2)))/b) );
	}
	dfdu += y_n * ( ((exp((t_c - u)/b)*exp(-exp((t_c - u)/b)))/b - (exp((t_w - u)/b)*exp(-exp((t_w - u)/b)))/b)/(exp(-exp((t_c - u)/b)) - exp(-exp((t_w - u)/b))) );
	dfdb += y_n * ( ((exp((t_c - u)/b)*exp(-exp((t_c - u)/b))*(t_c - u))/pow(b,2) - (exp((t_w - u)/b)*exp(-exp((t_w - u)/b))*(t_w - u))/pow(b,2))/(exp(-exp((t_c - u)/b)) - exp(-exp((t_w - u)/b))) );

	dfdu += (n-r-y_n)*( exp((t_w - u)/b)/b );
	dfdb += (n-r-y_n)*((exp((t_w - u)/b)*(t_w - u))/pow(b, 2));

	gsl_vector_set (f, 0, dfdu);
	gsl_vector_set (f, 1, dfdb);

	return GSL_SUCCESS;
}

int
print_state (int iter, gsl_multiroot_fsolver * s)
{
  printf ("iter = %3d x = % .3f % .3f "
          "f(x) = % .3e % .3e\n",
          iter,
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->f, 0),
          gsl_vector_get (s->f, 1));
}

int
main (void)
{
	const gsl_multiroot_fsolver_type *T;
	gsl_multiroot_fsolver *s;

	int status;
	int i, iter = 0;

	const int n = 2;

	double ftime[5] = {log(0.03400426), log(0.04005004), log(0.12370353), log(0.15395537), log(0.16553593)};
	struct rparams p = {9, 5, 100, log(0.18), log(0.3245928), ftime};
	gsl_multiroot_function f = {&loglik_f, n, &p};

	double gev_init[2] = {0, 0.5};
	gsl_vector *x = gsl_vector_alloc(n);

	gsl_vector_set(x, 0, gev_init[0]);
	gsl_vector_set(x, 1, gev_init[1]);

	T = gsl_multiroot_fsolver_hybrids;
	s = gsl_multiroot_fsolver_alloc (T, 2);
	gsl_multiroot_fsolver_set (s, &f, x);

	print_state(iter, s);

	do
	{
		iter++;
		status = gsl_multiroot_fsolver_iterate (s);

		print_state (iter, s);

		if (status)
			break;

		status = gsl_multiroot_test_residual (s->f, 1e-7);
	}
	while (status == GSL_CONTINUE && gsl_strerror(status));

	printf ("status = %s\n", gsl_strerror (status));
	printf("shape=%0.3f; scale= %0.3f\n", 1/gsl_vector_get (s->x, 1), exp(gsl_vector_get (s->x, 0)));
	gsl_multiroot_fsolver_free (s);
	gsl_vector_free (x);

	return 0;
}
