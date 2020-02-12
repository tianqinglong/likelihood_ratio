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

  return 0;
}

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
	const double log_b = gsl_vector_get(x, 1);

	int i;
	double b = exp(log_b);
	double dfdu = 0, dfdb = 0;
	for(i=0; i<r; i++)
	{
		dfdu += exp((ft[i]-u)/b)/b-1/b;
		dfdb += ( -b*exp(exp((ft[i] - u)/b) - (ft[i] - u)/b)*(exp((ft[i] - u)/b - exp((ft[i] - u)/b))/pow(b,2) + (exp((ft[i] - u)/b - exp((ft[i] - u)/b))*((ft[i] - u)/pow(b,2) - (exp((ft[i] - u)/b)*(ft[i] - u))/pow(b,2)))/b) );
	}
	dfdu += y_n * ( ((exp((t_c - u)/b)*exp(-exp((t_c - u)/b)))/b - (exp((t_w - u)/b)*exp(-exp((t_w - u)/b)))/b)/(exp(-exp((t_c - u)/b)) - exp(-exp((t_w - u)/b))) );
	dfdb += y_n * ( ((exp((t_c - u)/b)*exp(-exp((t_c - u)/b))*(t_c - u))/pow(b,2) - (exp((t_w - u)/b)*exp(-exp((t_w - u)/b))*(t_w - u))/pow(b,2))/(exp(-exp((t_c - u)/b)) - exp(-exp((t_w - u)/b))) );

	dfdu += (n-r-y_n)*( exp((t_w - u)/b)/b );
	dfdb += (n-r-y_n)*( (exp((t_w - u)/b)*(t_w - u))/pow(b, 2) );

	gsl_vector_set (f, 0, dfdu);
	gsl_vector_set (f, 1, dfdb);

	return GSL_SUCCESS;
}
