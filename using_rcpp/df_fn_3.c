double
my_f (const gsl_vector *v, void *params)
{
	double a, b, beta, eta, loglik=0;
	struct weibull_data_params *dat = (struct weibull_data_params *) params;

	int r = dat->r, n = dat->n;
	double t_c = dat->t_c;
	double* t_f = dat->t_f;

	a = gsl_vector_get(v, 0);
	b = gsl_vector_get(v, 1);

	beta = exp(a);
	eta = exp(b);

	for (int i=0; i < r; i++)
	{
		loglik += log(beta/eta*pow(t_f[i]/eta, beta-1))-pow(t_f[i]/eta, beta);
	}
	loglik += (n-r)*(-pow(t_c/eta, beta));

	double negloglik = -loglik;

	return negloglik;
}
