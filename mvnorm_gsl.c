/***************************************************************************************
 *  Multivariate Normal density function and random number generator
 *  Multivariate Student t density function and random number generator
 *  Wishart random number generator
 *  Using GSL -> www.gnu.org/software/gsl
 *
 *  Copyright (C) 2006  Ralph dos Santos Silva
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *  AUTHOR 
 *     Ralph dos Santos Silva,  [EMAIL PROTECTED]
 *     March, 2006       
***************************************************************************************/
#include <math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include "../hmm/adkGSL.h"
#include "mvnorm_gsl.h"
/*****************************************************************************************************************/
/*****************************************************************************************************************/
int rmvnorm(const gsl_rng *r, const int n, const gsl_vector *mean, const gsl_matrix *var, gsl_vector *result){
/* multivariate normal distribution random number generator */
/*
*	n	dimension of the random vetor
*	mean	vector of means of size n
*	var	variance matrix of dimension n x n
*	result	output variable with a sigle random vector normal distribution generation
*/
int k;
gsl_matrix *work = gsl_matrix_alloc(n,n);

gsl_matrix_memcpy(work,var);
gsl_linalg_cholesky_decomp(work);

for(k=0; k<n; k++)
	gsl_vector_set( result, k, gsl_ran_ugaussian(r) );

gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, work, result);
gsl_vector_add(result,mean);

gsl_matrix_free(work);

return 0;
}

/*****************************************************************************************************************/
/*****************************************************************************************************************/
int rmvnorm_prealloc(const gsl_rng *r, const int n, const gsl_vector *mean, const gsl_matrix *var, gsl_matrix *work, gsl_vector *result){
/* multivariate normal distribution random number generator */
/*
*	n	dimension of the random vetor
*	mean	vector of means of size n
*	var	variance matrix of dimension n x n
*	result	output variable with a sigle random vector normal distribution generation
*/
int k;

gsl_matrix_memcpy(work,var);
gsl_linalg_cholesky_decomp(work);

for(k=0; k<n; k++)
	gsl_vector_set( result, k, gsl_ran_ugaussian(r) );

gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, work, result);
gsl_vector_add(result,mean);


return 0;
}

/*****************************************************************************************************************/
/*****************************************************************************************************************/
double dmvnorm(const int n, const gsl_vector *x, const gsl_vector *mean, const gsl_matrix *var){
/* multivariate normal density function    */
/*
*	n	dimension of the random vetor
*	mean	vector of means of size n
*	var	variance matrix of dimension n x n
*/
int s;
double ax,ay;
gsl_vector *ym, *xm;
gsl_matrix *work = gsl_matrix_alloc(n,n), 
           *winv = gsl_matrix_alloc(n,n);
gsl_permutation *p = gsl_permutation_alloc(n);

gsl_matrix_memcpy( work, var );
gsl_linalg_LU_decomp( work, p, &s );
gsl_linalg_LU_invert( work, p, winv );
ax = gsl_linalg_LU_det( work, s );
gsl_matrix_free( work );
gsl_permutation_free( p );

xm = gsl_vector_alloc(n);
gsl_vector_memcpy( xm, x);
gsl_vector_sub( xm, mean );
ym = gsl_vector_alloc(n);
gsl_blas_dsymv(CblasUpper,1.0,winv,xm,0.0,ym);
gsl_matrix_free( winv );
gsl_blas_ddot( xm, ym, &ay);
gsl_vector_free(xm);
gsl_vector_free(ym);
ay = exp(-0.5*ay)/sqrt( pow((2*M_PI),n)*ax );

return ay;
}

/*****************************************************************************************************************/
/*****************************************************************************************************************/
int rmvt(const gsl_rng *r, const int n, const gsl_vector *location, const gsl_matrix *scale, const int dof, gsl_vector *result){
/* multivariate Student t distribution random number generator */
/*
*	n	 dimension of the random vetor
*	location vector of locations of size n
*	scale	 scale matrix of dimension n x n
*	dof	 degrees of freedom
*	result	 output variable with a single random vector normal distribution generation
*/
int k;
gsl_matrix *work = gsl_matrix_alloc(n,n);
double ax = 0.5*dof; 

ax = gsl_ran_gamma(r,ax,(1/ax));     /* gamma distribution */

gsl_matrix_memcpy(work,scale);
gsl_matrix_scale(work,(1/ax));       /* scaling the matrix */
gsl_linalg_cholesky_decomp(work);

for(k=0; k<n; k++)
	gsl_vector_set( result, k, gsl_ran_ugaussian(r) );

gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, work, result);
gsl_vector_add(result, location);

gsl_matrix_free(work);

return 0;
}

/*****************************************************************************************************************/
/*****************************************************************************************************************/
double dmvt(const int n, const gsl_vector *x, const gsl_vector *location, const gsl_matrix *scale, const int dof){
/* multivariate Student t density function */
/*
*	n	 dimension of the random vetor
*	location vector of locations of size n
*	scale	 scale matrix of dimension n x n
*	dof	 degrees of freedom
*/
int s;
double ax,ay,az=0.5*(dof + n);
gsl_vector *ym, *xm;
gsl_matrix *work = gsl_matrix_alloc(n,n), 
           *winv = gsl_matrix_alloc(n,n);
gsl_permutation *p = gsl_permutation_alloc(n);

gsl_matrix_memcpy( work, scale );
gsl_linalg_LU_decomp( work, p, &s );
gsl_linalg_LU_invert( work, p, winv );
ax = gsl_linalg_LU_det( work, s );
gsl_matrix_free( work );
gsl_permutation_free( p );

xm = gsl_vector_alloc(n);
gsl_vector_memcpy( xm, x);
gsl_vector_sub( xm, location );
ym = gsl_vector_alloc(n);
gsl_blas_dsymv(CblasUpper,1.0,winv,xm,0.0,ym);
gsl_matrix_free( winv );
gsl_blas_ddot( xm, ym, &ay);
gsl_vector_free(xm);
gsl_vector_free(ym);

ay = pow((1+ay/dof),-az)*gsl_sf_gamma(az)/(gsl_sf_gamma(0.5*dof)*sqrt( pow((dof*M_PI),n)*ax ));

return ay;
}
/*****************************************************************************************************************/
/*****************************************************************************************************************/
int rwishart(const gsl_rng *r, const int n, const int dof, const gsl_matrix *scale, gsl_matrix *work, gsl_matrix *output){
/* Wishart distribution random number generator */
/*
*	n	 gives the dimension of the random matrix
*	dof	 degrees of freedom
*	scale	 scale matrix of dimension n x n
*	result	 output variable with a single random matrix Wishart distribution generation
*/
int k,l;

for(k=0; k<scale->size1; k++){
	gsl_matrix_set( work, k, k, sqrt( gsl_ran_chisq( r, (dof-k) ) ) );
	for(l=0; l<k; l++){
		gsl_matrix_set( work, k, l, gsl_ran_ugaussian(r) );
	}
}
gsl_matrix_memcpy(output,scale);
gsl_linalg_cholesky_decomp(output);
gsl_blas_dtrmm(CblasLeft,CblasLower,CblasNoTrans,CblasNonUnit,1.0,output,work);
gsl_blas_dsyrk(CblasUpper,CblasNoTrans,1.0,work,0.0,output);
return 0;
}

// ADK additions from here on

/*****************************************************************************************************************/
/*****************************************************************************************************************/
double dmvnorm_prealloc(const int n, const gsl_vector *x, const gsl_vector *mean, const gsl_matrix *var, gsl_matrix *work, 
	 gsl_matrix *winv,  gsl_vector *ym,  gsl_vector *xm, gsl_permutation *p){
/* multivariate normal density function    */
/*
*	n	dimension of the random vetor
*	mean	vector of means of size n
*	var	variance matrix of dimension n x n
*
*   added this prealloc version to cut down on memory calls
*/

int s;
double ax,ay;

gsl_matrix_memcpy( work, var );
gsl_linalg_LU_decomp( work, p, &s );
gsl_linalg_LU_invert( work, p, winv );
ax = gsl_linalg_LU_det( work, s );
gsl_vector_memcpy( xm, x);
gsl_vector_sub( xm, mean );
gsl_blas_dsymv(CblasUpper,1.0,winv,xm,0.0,ym);
gsl_blas_ddot( xm, ym, &ay);
ay = exp(-0.5*ay)/sqrt( pow((2*M_PI),n)*ax );

return ay;
}

/*****************************************************************************************************************/
/*****************************************************************************************************************/
double dInverseWishart(const gsl_matrix *x, const gsl_matrix *sigma, const double dof, gsl_matrix *work, 
	 gsl_matrix *winv,  gsl_permutation *p){
	
	double denom, num;
	int s, d, i;
	double detSig,detX,trace;
	
	//get det of both matrices, take inverse of X	
	gsl_matrix_memcpy( work, sigma );
	gsl_linalg_LU_decomp( work, p, &s );
	detSig = gsl_linalg_LU_det( work, s );
	gsl_matrix_memcpy( work, x );
	gsl_linalg_LU_decomp( work, p, &s );
	detX = gsl_linalg_LU_det( work, s );
	gsl_linalg_LU_invert( work, p, winv );
	//get product of gamma's for denominator; calc denom
	d = x->size1;	
	denom = 1.0;
	denom *= pow(2,(d*dof/2.0)) * mv_gamma_func(dof/2.0, d);
	
	//multiply x by inverse of sigma and take trace
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,1.0, sigma, winv,0.0, work);
	trace = gsl_matrix_trace(work);
	num =  pow(detSig, dof / 2.0) * pow(detX , (-dof-d-1)/ 2.0) * exp(-0.5 * trace);
	return(num/denom);
}

/*****************************************************************************************************************/
//rwishartGelman-- uses alg described in Gelman et al appendix
void rwishartGelman(const gsl_rng *r, const int n, const int dof, const gsl_matrix *scale, 
	gsl_matrix *work,gsl_matrix *work2, gsl_vector *mean,gsl_vector *xm, gsl_matrix *output){
/* Wishart distribution random number generator */
/*
*	n	 gives the dimension of the random matrix
*	dof	 degrees of freedom
*	scale	 scale matrix of dimension n x n
*	result	 output variable with a single random matrix Wishart distribution generation
*/
	int i;

	gsl_vector_set_all(mean,0.0);
	gsl_matrix_set_all(output,0.0);
	for(i = 0; i < dof; i++){
		gsl_vector_set_all(xm,0);
		rmvnorm_prealloc(r, n, mean, scale, work, xm);
		gsl_vector_outer_product(xm,xm,work2);
		gsl_matrix_add(output,work2);
	}
}

//rwishartOdellFeiveson-- uses alg of Odell and Feiveson
void rwishartOdellFeiveson(const gsl_rng *r, const int n, const int dof, const gsl_matrix *scale, 
	gsl_matrix *work,gsl_matrix *work2, gsl_vector *mean,gsl_vector *xm, gsl_matrix *output){
/* Wishart distribution random number generator */
/*
*	n	 gives the dimension of the random matrix
*	dof	 degrees of freedom
*	scale	 scale matrix of dimension n x n
*	result	 output variable with a single random matrix Wishart distribution generation
*/
	int i,j, k;
	double sum, tmp;

	gsl_vector_set_all(mean,0.0);
	gsl_matrix_set_all(output,0.0);

	for(j=0;j<n;j++){
		gsl_vector_set(mean,j,gsl_ran_chisq(r,dof-j));
		for(i=0;i<j;i++){
			gsl_matrix_set(work,i,j,gsl_ran_ugaussian(r));
		}
	}
	
	//set up B matrix in work2
	gsl_matrix_set(work2,0,0,gsl_vector_get(mean,0));
	for(j=0;j<n;j++){
		for(i=0;i<j;i++){
			if(i==0){
				tmp = sqrt(gsl_vector_get(mean,j)) * gsl_matrix_get(work,0,j);
				gsl_matrix_set(work2,0,j,tmp);
			}
			else{
				tmp = sqrt(gsl_vector_get(mean,i)) * gsl_matrix_get(work,i,j);
				for(k=0;k<i;k++){
					tmp += gsl_matrix_get(work,k,i) * gsl_matrix_get(work,k,j);
				}
				gsl_matrix_set(work2,i,j,tmp);
			}
		}
		tmp = gsl_vector_get(mean,j);
		for(k=0;k<j;k++){
			tmp += gsl_matrix_get(work,k,j) * gsl_matrix_get(work,k,j);
		}
		gsl_matrix_set(work2,j,j,tmp);
	}
	gsl_matrix_memcpy(output,scale);
	gsl_linalg_cholesky_decomp(output);
	gsl_blas_dtrmm(CblasLeft,CblasLower,CblasNoTrans,CblasNonUnit,1.0,output,work2);
	gsl_blas_dtrmm(CblasRight,CblasUpper,CblasNoTrans,CblasNonUnit,1.0,work2,output);
	
}

//mv_gamma_func -- returns the multivariate gamma function
double mv_gamma_func(double a, double d){
	double val = 1.0;
	int i;
	
	for(i=1;i<=d;i++){
		 val *= gsl_sf_gamma(a- (0.5 * (i-1)));
		}
	val *=  pow(M_PI,(d*(d-1)/4.0));
	return(val);
}

