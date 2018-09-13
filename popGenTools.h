/* popGenTools.h
/
/
/ Andrew Kern
/
/
*/
#ifndef PGTOOLS_INC
#define PGTOOLS_INC


//#define MAXSNPS 1000000
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_min.h>
#include "mpi.h"


	
struct my_F_params{
	int i, n, ascSize, ascFreq, derivedFlag;
	gsl_integration_workspace * w;
};

struct trans_params{
	int xi, i, n, ascSize, ascFreq, derivedFlag;
	double t, q, p, tau,v;	
	gsl_integration_workspace * w;
};

struct prf_params{
	int xi, i, n, ascSize, ascFreq, derivedFlag;	
};

struct my_snpProb_params{
  int i, n, ascSize;
  double denom;
	gsl_integration_workspace * w;
};

struct my_f_params{
  int i,n, ascSize;
  double beta;
};

struct prfLik_params{
  int snpNumber, maxSampleSize, ascSize, runMode;
  gsl_vector *weights, *sampleSizeVector, *betas;
  struct snp *data;
  gsl_matrix *sfsCounts;		
};

struct my_lik_params{
  int snpNumber, maxSampleSize, ascSize, runMode;
  gsl_vector *weights, *sampleSizeVector, *betas;
  gsl_integration_workspace * w;
  struct snp *data;
  gsl_matrix *sfsBools;
  gsl_min_fminimizer *min;

};
struct my_likCI_params{
  int snpNumber, maxSampleSize;
  gsl_vector *weights, *sampleSizeVector;
  struct snp *data;
  gsl_matrix *sfsBools;
  double lMax, logUnits, beta_hat;
};


struct my_jointLik_params{
  int snpNumber, maxSampleSize, ascSize, nstates;
  gsl_vector *weights, *sampleSizeVector, *betas, *startingPoint;
  struct snp *data;
  gsl_matrix *sfsBools, *posts;

};

double my_f(double q, void * params);
double my_F(double beta, void * params);
double my_F2(double beta, void * p);
double my_FAdapt(double beta, void * p);
double my_lik(double beta, void * params);
double bustProb(double beta, double theta, void *p);
double sfsLikBetaVector(gsl_vector *beta, void * p);
double prfLik(double beta, void * p);
double weightedLik(double beta, void * params);
double weightedLikLook(double beta, void * params);
double weightedLikLookCI(double beta, void * params);
double likWrap(double beta);
double snpProb(double beta, void * params);
double snpProbDenomLookup(double beta, void * params);
double snpProbDenom(double beta, int sampleSize);
gsl_vector *makeSnpProbDenomVector(double beta, int maxSampleSize, gsl_vector *sampleSizeVec);
double ml_est(double * lik_beta_hat);
double ml_prf(double * lik_beta_hat, void * p);
double weighted_ml_est(double * lik_beta_hat, void * p);
double weighted_ml_est_lookup(double * lik_beta_hat, void * p);
double weighted_ml_est_lookup_errReport(double * lik_beta_hat, void * p);
double weighted_ml_CILower_lookup(double * lik_beta_hat, void * p);
double weighted_ml_CIUpper_lookup(double * lik_beta_hat, void * p);
double wlikWrap(double beta, void * p);
double  wlikWrapLook(double beta, void * p);
double  wlikCIWrapLook(double beta, void * p);
double  prfLikWrap(double beta, void * p);
gsl_vector *sampleSizeVector(struct snp *data, int snpNumber, int maxSampleSize);
void sampleSizeVectorFill(struct snp *data, int snpNumber, int maxSampleSize, gsl_vector *bools);
gsl_matrix *summarizeSFS(int maxSampleSize, struct snp *data, int snpNumber);
gsl_matrix *summarizeSFSBool(int max, struct snp *data, int snpNumber);
void summarizeSFSBoolFill(gsl_matrix *bool, struct snp *data, int snpNumber);
void summarizeSFSFill(gsl_matrix *counts, struct snp *data, int snpNumber);
int summarizeSFSFillFromTo(gsl_matrix *counts, struct snp *data, int snpNumber, int from, int to, int startPos);
gsl_matrix *snpProbMatrix(double beta,  int maxSampleSize, gsl_vector *sampleSizeVector, gsl_matrix *sfsBools);
gsl_matrix *snpProbMatrixNotLog(double beta,  int maxSampleSize, gsl_vector *sampleSizeVector, gsl_matrix *sfsBools);
gsl_matrix *snpProbMatrixNotLogFull(double beta,  int maxSampleSize,gsl_vector *sampleSizeVector);
gsl_vector *snpProbVectorNotLog(double beta,  int sampleSize);
int simulateSiteFreq(int sampleSize, double alpha, void *r);
double simpleAscertain(double sampleFreq, void *p);
double probAscertainmentGivenModel(double beta, void *p);
gsl_matrix *snpAscMatrix(int maxSampleSize, gsl_vector *sampleSizeVector, gsl_matrix *sfsBools, int ascSize);
gsl_matrix *estimateAscSFS(int maxSampleSize, gsl_vector *sampleSizeVector,  gsl_matrix *sfsBools, int ascSize, gsl_matrix *sfsSummary);
double weightedLikLookAsc(double beta, void * p);
double  wlikWrapLookAsc(double beta, void * p);
double weighted_ml_est_lookup_asc(double * lik_beta_hat, void * p);
double probAscertainmentGivenModelLookup(double beta, int sampSize, gsl_matrix *snpProbs, gsl_matrix *ascProbs);
double probAscertainmentGivenModelHemiLookup(double beta, int sampSize, gsl_vector *snpProbs, gsl_matrix *ascProbs);
double outgroupAscertain(double sampleFreq, void * p);
gsl_matrix *snpOutgroupAscMatrix(int maxSampleSize, gsl_vector *sampleSizeVector, int ascSize, int derived, int runMode);
double weightedLikLookOutgroupAsc(double beta, void * p);
double  wlikWrapLookOutgroupAsc(double beta, void * p);
double weighted_ml_est_lookup_outgroup_asc(double * lik_beta_hat, void * p);
double probOutgroupAscertainmentGivenModel(double beta, void *p);
double sfsLikBetaVectorOutgroupAsc(gsl_vector *betas, void * p);
double sfsLikBetaVectorOutgroupAscPre(gsl_vector *betas, gsl_matrix *ascProbsAncestral,gsl_matrix *ascProbsDerived, gsl_vector *snpProbs, void * p);
//pre-alloced versions
gsl_matrix *snpOutgroupAscMatrixPre(int maxSampleSize, gsl_vector *sampleSizeVector, int ascSize, int derived, \
	int runMode, gsl_matrix *probs);
gsl_vector *snpProbVectorNotLogPre(double beta,  int sampleSize, gsl_vector *probs);
//MPI stuff
double sfsLikBetaVectorMPI(gsl_vector *betas, void * p, int myrank, int mprocs, MPI_Comm comm);
double sfsLikBetaVectorOutgroupAscPreMPI(gsl_vector *betas, gsl_matrix *ascProbsAncestral,gsl_matrix *ascProbsDerived, \
	gsl_vector *snpProbs, void * p,int myrank, int mprocs, MPI_Comm comm);

double weightedLikLookFolded(double beta, void * p);
double weighted_ml_est_lookup_folded(double * lik_beta_hat, void * p);
double weighted_ml_CILower_lookupFolded(double * lik_beta_hat, void * p);
double weighted_ml_CIUpper_lookupFolded(double * lik_beta_hat, void * p);
double  wlikWrapLookFolded(double beta, void * p);


int segSites(struct snp *data, int sites);
double theta_pi(struct snp *data, int sites);
double a1f(int n);
double theta_w(int segSites, int n);
double a2f(int);
double b1f(int);
double b2f(int);
double c1f(double, double);
double c2f(int, double, double, double);
double e1f(double, double);
double e2f(double, double, double);
double tajdSNP(int, int, double) ;



#endif


