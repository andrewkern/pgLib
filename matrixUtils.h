gsl_matrix fill_rand_matrix(gsl_matrix *X, gsl_rng *r_num);
gsl_matrix fill_matrix(gsl_matrix *X);
double matrix_determ(gsl_matrix *X);
double matrix_trace(gsl_matrix *X);
double mv_gamma(double a, double d);
double iwishpdf(gsl_matrix *X, gsl_matrix *Scale, gsl_matrix *inv, double dof);
gsl_matrix inv_matrix(gsl_matrix *X, gsl_matrix *inv);
void print_matrix(gsl_matrix *X);
gsl_vector array_to_gsl_vec(gsl_vector *dest, double src[3]);
gsl_matrix mdarray_to_gsl_mat(gsl_matrix *dest, double src[3][3]);
void print_vector(gsl_vector *src);
gsl_matrix get_chol(gsl_matrix *src, gsl_matrix *dest);
gsl_matrix mat_noise(gsl_matrix *src, gsl_rng *r_num);
gsl_vector vec_noise(gsl_rng *r_num, gsl_vector *src, gsl_vector *dest);
void fill_chol_vec(gsl_matrix *src, gsl_vector *dest);
gsl_matrix chol_mult(gsl_matrix *src, gsl_matrix *dest);
gsl_matrix gsl_matrix_lower_tri(gsl_matrix *src);

