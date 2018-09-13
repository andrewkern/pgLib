int rmvnorm(const gsl_rng *r, const int n, const gsl_vector *mean, const gsl_matrix *var, gsl_vector *result);
double dmvnorm(const int n, const gsl_vector *x, const gsl_vector *mean, const gsl_matrix *var);
int rmvt(const gsl_rng *r, const int n, const gsl_vector *location, const gsl_matrix *scale, const int dof, gsl_vector *result);
double dmvt(const int n, const gsl_vector *x, const gsl_vector *location, const gsl_matrix *scale, const int dof);
int rwishart(const gsl_rng *r, const int n, const int dof, const gsl_matrix *scale, gsl_matrix *work, gsl_matrix *result);
double mv_gamma_func(double a, double d);

/// mine
int rmvnorm_prealloc(const gsl_rng *r, const int n, const gsl_vector *mean, const gsl_matrix *var, gsl_matrix *work, gsl_vector *result);
double dmvnorm_prealloc(const int n, const gsl_vector *x, const gsl_vector *mean, const gsl_matrix *var, gsl_matrix *work, 
	gsl_matrix *winv, gsl_vector *ym, gsl_vector *xm, gsl_permutation *p);

double dInverseWishart(const gsl_matrix *x, const gsl_matrix *sigma, const double dof, gsl_matrix *work, 
	gsl_matrix *winv,  gsl_permutation *p);
	
void rwishartGelman(const gsl_rng *r, const int n, const int dof, const gsl_matrix *scale, 
	gsl_matrix *work,gsl_matrix *work2, gsl_vector *mean,gsl_vector *xm, gsl_matrix *output);
	
void rwishartOdellFeiveson(const gsl_rng *r, const int n, const int dof, const gsl_matrix *scale, 
	gsl_matrix *work,gsl_matrix *work2, gsl_vector *mean,gsl_vector *xm, gsl_matrix *output);