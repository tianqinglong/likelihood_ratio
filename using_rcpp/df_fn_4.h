struct rparams
{
	int y_n;
	int r;
	int n;
	double t_c;
	double t_w;
	double *ft;
};

double
my_f (const gsl_vector *v, void *params);