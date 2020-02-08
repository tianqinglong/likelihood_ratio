struct weibull_data_params
{
	int r;
	int n;
	double t_c;
	double* t_f;
};

double
my_f (const gsl_vector *v, void *params);