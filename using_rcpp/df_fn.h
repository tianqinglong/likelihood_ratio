struct weibull_data_params
{
	int r;
	int n;
	double t_c;
	double* t_f;
};

double deriv(double beta, void *params);
