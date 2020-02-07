struct rparams
{
	int y_n;
	int r;
	int n;
	double t_c;
	double t_w;
	double *ft;
};

int print_state (int iter, gsl_multiroot_fsolver * s);
int loglik_f (const gsl_vector * x, void * params, gsl_vector * f);