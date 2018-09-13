//popGen stuff for Next Gen Sequencing Data

#include "pgSummaryStats.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "vector.h"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>

//Sampling Probability stuff


//prob_i returns the probability under SNM of a snp freq i
double prob_i(int i, int n){
	double tmp;
	tmp = 1.0 / (double) i;
	return(tmp/ a1f(n));
}


double prob_km(int k, int m, int n){
	double sum = 0.0;
	int i;
	double tmpFreq;
	
	for(i=1;i<n;i++){
		tmpFreq = (double) i / (double) n;
		sum += gsl_ran_binomial_pdf(k,tmpFreq,m) * prob_i(i,n);
	}
	return(sum);
}



//achaz system stuff
double achazSFS(gsl_matrix *sfsCount, gsl_matrix *weightMat){
	int i, j;
	double theta = 0.0;
	double weightSum, tmpTheta;
	
	//rows = sampleSizes
	for(i = 2; i < sfsCount->size1; i++){
		weightSum = 0.0;
		tmpTheta = 0.0;
		for(j = 1; j < i; j++){
			//size i freq j
			tmpTheta += gsl_matrix_get(sfsCount, i, j) * j * gsl_matrix_get(weightMat,i,j);
			weightSum += gsl_matrix_get(weightMat,i,j);
		}
		if(weightSum > 0){ theta += tmpTheta / weightSum; }
	}
	return(theta);
}

double achazSFSCorrected(gsl_matrix *sfsCount, gsl_matrix *weightMat,int nPool){
	int i, j;
	double theta = 0.0;
	double weightSum, tmpTheta;
	
	//rows = sampleSizes
	for(i = 2; i < sfsCount->size1; i++){
		weightSum = 0.0;
		tmpTheta = 0.0;
		for(j = 1; j < i; j++){
			//size i freq j
			tmpTheta += gsl_matrix_get(sfsCount, i, j) * gsl_matrix_get(weightMat,i,j) * (1.0/prob_km(j,i,nPool));
			weightSum += gsl_matrix_get(weightMat,i,j);
		}
		if(weightSum > 0){ 
		
			theta += tmpTheta / weightSum / a1f(nPool); 
		}
	}
	return(theta);
}


// fillPiWeights -- returns the weight matrix for theta_pi
void fillPiWeights(gsl_matrix *weightMat, int minFreq){
	int i, j;
	gsl_matrix_set_zero(weightMat);
	for(i=minFreq+1; i < weightMat->size1; i++){
		for(j=minFreq;j<i;j++){
			gsl_matrix_set(weightMat,i,j,i-j);
		}
	}
}

// fillHWeights -- returns the weight matrix for theta_H
void fillHWeights(gsl_matrix *weightMat, int minFreq){
	int i, j;
	gsl_matrix_set_zero(weightMat);
	for(i=minFreq+1; i < weightMat->size1; i++){
		for(j=minFreq;j<i;j++){
			gsl_matrix_set(weightMat,i,j,j);
		}
	}
}

// fillWWeights -- returns the weight matrix for theta_w
void fillWWeights(gsl_matrix *weightMat, int minFreq){
	int i, j;
	gsl_matrix_set_zero(weightMat);
	for(i=minFreq+1; i < weightMat->size1; i++){
		for(j=minFreq;j<i;j++){
			gsl_matrix_set(weightMat,i,j,1/(double)j);
		}
	}
}
