#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_permutation.h>
#include "matrixUtils.h"
#include "../mcmc/mvnormMCMC/mvrandist.h"

/*A collection of matrix functions*/
/*Alexander G. Shanku*/
/*Wed 12/05/2012 12:44:41 EST*/

//// Makes a matrix filled with random U(0,1) elements ////
gsl_matrix fill_rand_matrix(gsl_matrix *X, gsl_rng *r_num){
	int i,j;
	/*int step = 0;*/
	int n = X->size1;
	int m = X->size2;
	for (i = 0; i < m; ++i){
		for (j = 0; j < n; ++j){
			gsl_matrix_set(X, i, j, gsl_rng_uniform(r_num));						
		}
	}
	return(*X);	
}

//// Makes a matrix filled with sequential elements, starting with 0 ////
gsl_matrix fill_matrix(gsl_matrix *X){
	int i,j;
	int n = X->size1;
	int m = X->size2;
	for (i = 0; i < m; ++i){
		for (j = 0; j < n; ++j){
			gsl_matrix_set(X, i, j, i+(j*m));
		}
	}
	return(*X);	
}

//// Do LU decomp, get determinant and keep input matrix untouched ////
double matrix_determ(gsl_matrix *X){
	int s;
	int n = X->size1;
	int m = X->size2;
	gsl_matrix *a_copy = gsl_matrix_alloc(m, n);
	gsl_matrix_memcpy(a_copy, X );
	gsl_permutation *P = gsl_permutation_alloc(n);
	gsl_linalg_LU_decomp(a_copy, P, &s);
	double my_det = gsl_linalg_LU_det (a_copy, s);
	return(my_det);
}

//// Get trace of given matrix ////
double matrix_trace(gsl_matrix *X){
	int i, m;
	m = X->size1;
	double trace = 0.0;
	for(i=0;i<m;i++){
		trace += gsl_matrix_get(X,i,i);
	}
	return(trace);
}

//// Multivariate Gamma function ////
double mv_gamma(double a, double d){
	double val = 1.0;
	int i;
		for(i = 1; i <= d; i++){
		 val *= gsl_sf_gamma(a - (0.5 * (i - 1)));
		}
	val *=  pow(M_PI, (d * (d - 1) / 4.0));
	return(val);
}

//// Get Inverse Wishart PDF ////
double iwishpdf(gsl_matrix *X, gsl_matrix *Scale, gsl_matrix *inv, double dof){
	// From http://en.wikipedia.org/wiki/Inverse-Wishart_distribution
	double X_det, scale_det, denom, pdf, trace, numerator;	
	int m = X->size1;
	int n = X->size2;
	
	// Allocate matrix for inv(X)
	gsl_matrix *for_mult = gsl_matrix_alloc(m, n);
	
	// Get determinant of X and Scale matrix
	X_det = matrix_determ(X);
	scale_det = matrix_determ(Scale);
	
	// Invert X
	inv_matrix(X,inv);
	
	// Multiple Scale * inv(X)
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Scale, inv, 1.0, for_mult);
	
	// Get trace of above.
	trace = matrix_trace(for_mult);

	numerator =  pow(scale_det, dof / 2.0) * pow(X_det, (-dof-m-1)/ 2.0) * exp(-0.5 * trace);
	denom = pow(2,dof * m / 2) * mv_gamma(dof/2, m);

	pdf = (numerator/denom);

	return(pdf);
}

//// Invert a matrix and keep input matrix untouched ////
gsl_matrix inv_matrix(gsl_matrix *X, gsl_matrix *inv){
	int s;
	int n = X->size1;
	int m = X->size2;
	gsl_matrix *a_copy = gsl_matrix_alloc(m, n);
	gsl_matrix_memcpy( a_copy, X );
	gsl_permutation *P = gsl_permutation_alloc(n);
	gsl_linalg_LU_decomp(a_copy, P, &s);
	gsl_linalg_LU_invert(a_copy, P, inv);
	return(*inv);
}

//// Print a matrix to stdout ////
void print_matrix(gsl_matrix *X){
	int i, j;
	int n = X->size1;
	int m = X->size2;
	printf("\n");
	for(i=0; i<n; i++){
		for(j=0; j<m; j++){
			printf("%05.2f ", gsl_matrix_get(X, i, j));
		}
		printf("\n");
	}
}


///////// Added Sat 01/19/2013 15:24:41 EST /////////		

gsl_vector array_to_gsl_vec(gsl_vector *dest, double src[3]){
	int i;
	int n = dest->size;
		for (i = 0; i < n; ++i){
			
				gsl_vector_set(dest, i, src[i]);
			}
	return(*dest);	
}	

///// Transform user supplied matrix into a gsl matrix ////
gsl_matrix mdarray_to_gsl_mat(gsl_matrix *dest, double src[3][3]){
	int i, j;
	int n = dest->size1;
	int m = dest->size2;
		for (i = 0; i < m; ++i){
			for (j = 0; j < n; ++j){
				gsl_matrix_set(dest, i, j, src[i][j]);
			}
		}
	return(*dest);	
}	

//// Print a vector to stdout ////
void print_vector(gsl_vector *src){
	int i;
	int n = src->size;
	printf("\n");
	for (i = 0; i < n; ++i){
		printf("%05.4f ", gsl_vector_get(src,i));
	}
	printf("\n");
}

//// Do cholesky decomp to matrix ////
gsl_matrix get_chol(gsl_matrix *src, gsl_matrix *dest){
	gsl_matrix_memcpy(dest, src);
	gsl_linalg_cholesky_decomp(dest);
	return(*dest);
}

//// Add gaussian noise to chol decomp matrix ////
gsl_matrix mat_noise(gsl_matrix *src, gsl_rng *r_num){
	int i,j;
	int n = src->size1;
	int m = src->size2;
	double tmp;

	for(i = 0; i < n; i++){
		for(j = 0; j < m; j++){
			if(i >= j){
				tmp = gsl_matrix_get(src,i,j) + gsl_ran_gaussian(r_num,0.05);
				gsl_matrix_set(src, i, j, tmp);
			}
			
		}
	}	
	return(*src);
}

//// Transform cholesky decomp matrix into vector ////
void fill_chol_vec(gsl_matrix *src, gsl_vector *dest){
	int i, j, count = 0;
	
	for(i = 0; i < src->size1; i++){
		for(j=0;j<src->size2;j++){
			if(i<=j){
				 gsl_vector_set(dest,count++, gsl_matrix_get(src,i,j));
			}
			
		}
	}	
}

//// Add gaussian noise to each element in vector ////
gsl_vector vec_noise(gsl_rng *r_num, gsl_vector *src, gsl_vector *dest){
	int i; 
	int n = src->size;
	double tmp;
	
	for(i=0;i<n;i++){
			tmp = gsl_vector_get(src,i) + gsl_ran_gaussian(r_num, 0.05);
			gsl_vector_set(dest,i,tmp);
		}
		return(*dest);
}

//// Multiply cholesky decomp matrix by its transpose ////
gsl_matrix chol_mult(gsl_matrix *src, gsl_matrix *dest){
	gsl_blas_dgemm(CblasNoTrans, CblasTrans,1.0, src, src,1.0, dest);
	return(*dest);
}

//// Zeros out everything but lower tri-- in place ////
//// From adk's adkGSL.c ////
gsl_matrix gsl_matrix_lower_tri(gsl_matrix *src){
	int i, j;
	for(i = 0; i < src->size1; i++){
		for(j = 0; j < src->size2; j++){
			if(i<j)	gsl_matrix_set(src,i,j,0.0);
		}
	}
	return(*src);
}


// Random Wishart matrix.  Taken from Andy Kern's mvnorm_gsl.c
/*void rwishartGelman(gsl_matrix *Scale, double dof){
*//* Wishart distribution random number generator */
/*
*	n	 gives the dimension of the random matrix
*	dof	 degrees of freedom
*	scale	 scale matrix of dimension n x n
*	result	 output variable with a single random matrix Wishart distribution generation
*/
/*	int i;
	gsl_vector gsl_vector_alloc(ndim)

	gsl_vector_set_all(mean,0.0);
	gsl_matrix_set_all(output,0.0);
	for(i = 0; i < dof; i++){
		gsl_vector_set_all(xm,0);
		rmvnorm_prealloc(r, n, mean, scale, work, xm);
		gsl_vector_outer_product(xm,xm,work2);
		gsl_matrix_add(output,work2);
	}
}
*/


/*void rwishartGelman(const gsl_rng *r, const int n, const int dof, const gsl_matrix *scale, 
	gsl_matrix *work,gsl_matrix *work2, gsl_vector *mean,gsl_vector *xm, gsl_matrix *output)

rwishartGelman(r, ndim, nsnps+ndim, p2->winv, 
			p2->work,p2->work2, p2->means, p2->xm, p2->covMat);

*/





























