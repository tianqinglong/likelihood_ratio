double
my_f (const gsl_vector *v, void *params)
{
	double a, b, beta, eta, loglik=0;
	struct rparams *dat = (struct rparams *) params;

	int r = dat->r, n = dat->n, y_n = dat->y_n;
	double t_c = dat->t_c, t_w = dat->t_w;
	double* ft = dat->ft;

	a = gsl_vector_get(v, 0);
	b = gsl_vector_get(v, 1);

	beta = exp(a);
	eta = exp(b);

	for (int i=0; i < r; i++)
	{
		loglik += log(beta/eta*pow(ft[i]/eta, beta-1))-pow(ft[i]/eta, beta);
	}
	loglik += (n-r-y_n)*(-pow(t_c/eta, beta));
	loglik += y_n * log(exp(-pow(t_c/eta, beta)) - exp(-pow(t_w/eta, beta)));

	double negloglik = -loglik;

	return negloglik;
}
