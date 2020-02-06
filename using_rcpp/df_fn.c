double deriv(double beta, void *params)
{
	struct weibull_data_params *dat = (struct weibull_data_params *) params;

	int r = dat->r, n = dat->n;
	double t_c = dat->t_c;
	double* t_f = dat->t_f;

	double df1;
	double lambda = 0, df;

	int i;
	for (i=0; i<r; i++)
	{
		lambda += pow(t_f[i], beta);
	}
	lambda += (n-r)*pow(t_c, beta);

	lambda = lambda/r;
	lambda = pow(lambda, 1/beta);

	df1 = r/beta - r*log(lambda);
	for(i=0; i<r; i++)
	{
		df1 += log(t_f[i]);
		df1 -= pow(t_f[i]/lambda, beta) * (log(t_f[i])-log(lambda));
	}

	df1 -= (n-r)*pow(t_c/lambda, beta)*(log(t_c)-log(lambda));

	return df1;
}
