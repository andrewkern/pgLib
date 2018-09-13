/* population genetics tools- mainly for working with site frequency spectrum
/
/
/ Andrew Kern
*/


#include "stdio.h"
#include <math.h>
#include "snpFile.h"
#include "popGenTools.h"
#include "numerical.h"
#include "assert.h"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h> 
#include <gsl/gsl_sf.h>

/* function for sfs generation -- stationary solution for snp frequency density*/
double my_f(double q, void * p){
  struct my_f_params * params = (struct my_f_params *) p;
  double beta = params->beta;
  int i = params->i;
  int n = params->n;
	
  return ((1.0 - exp(-2.0 * beta * (1.0 - q))) / (1.0 - exp(-2.0 * beta))) * (2.0 / (q * (1.0 - q))) * gsl_ran_binomial_pdf(i, q, n);
}



/* integral of f over allele freq */
double my_F(double beta, void * p){
  struct my_F_params * params = (struct my_F_params *) p;
  struct my_f_params fParams;
  gsl_function f;
  //gsl_integration_workspace * w = gsl_integration_workspace_alloc (100000);
  double result,error;
  size_t size= 1000;
  
  if (beta == 0){
    //gsl_integration_workspace_free(w);
    return 1.0 / params->i;
  }
  else{
    fParams.beta = beta;
    fParams.i = params->i;
    fParams.n = params->n;
    f.function = &(my_f);
    f.params = &fParams;
    //adaptive integration-- slow and steady
//	gsl_integration_qags(&f,0.0,1.0,1e-5,1e-5, 100000,w, &result,&error);
	//non-adaptive version is rougher and faster
	gsl_integration_qng(&f, 0.0, 1.0, 1e-2,1e-2,&result,&error,&size);
//    gsl_integration_workspace_free(w);
    return result;
  }
}


/* integral of f over allele freq using numerical recipes integration */
double my_F2(double beta, void * p){
  struct my_F_params * params = (struct my_F_params *) p;
  struct my_f_params fParams;
  double result;
  
  if (beta == 0){
    return 1.0 / params->i;
  }
  else{
    fParams.beta = beta;
    fParams.i = params->i;
    fParams.n = params->n;
    result = d_qromb2(&my_f, 0.0, 1.0, &fParams);
    return result;
  }
}

/* integral of f over allele freq-- using adaptive integration */
double my_FAdapt(double beta, void * p){
  struct my_F_params * params = (struct my_F_params *) p;
  struct my_f_params fParams;
  gsl_function f;
  double result,error;
  size_t size= 100;
  
  if (beta == 0){
    return 1.0 / params->i;
  }
  else{
    fParams.beta = beta;
    fParams.i = params->i;
    fParams.n = params->n;
    f.function = &(my_f);
    f.params = &fParams;
    //adaptive integration-- slow and steady
	gsl_integration_qags(&f,0.0,1.0,1e-3,1e-3,size,params->w, &result,&error);
    return result;
  }
}


/* this returns the prob of SNP freq i given n and beta */
double snpProb(double beta, void * p){			
  struct my_F_params * params = (struct my_F_params *) p;
  struct my_F_params tempParams;
  long double tot;
  int i;
  gsl_function F;
  
  /* add up total prob space */
  tot = 0.0;
  tempParams.n = params->n;

  for(i = 1; i < params->n ; i++){
    tempParams.i = i;
    tempParams.w = params->w;
    F.function =  &my_F;
    F.params = &tempParams;
    tot += GSL_FN_EVAL(&F, beta);
  }
  
  /* reset parameters */	
  F.function =  &my_FAdapt;
  F.params = params;
  
  return   GSL_FN_EVAL(&F, beta) / tot ;
}		

/* trying to speedup the snpProb routine with a lookup of the denominators (in snpProb_params) */
double snpProbDenomLookup(double beta, void * p){			
	struct my_snpProb_params * params = (struct my_snpProb_params *) p;
	struct my_F_params tempParams;
	gsl_function F;

/* set parameters */	
	tempParams.n = params->n;
	tempParams.i = params->i;
	tempParams.w = params->w;
	F.function =  &my_FAdapt;
	F.params = &tempParams;

	return   GSL_FN_EVAL(&F, beta) / params->denom ;
}	

/* sampleSize vector - returns a vector of 1s or 0s depending on sample sizes in the data */
gsl_vector *sampleSizeVector(struct snp data[], int snpNumber, int maxSampleSize){
  gsl_vector *bools;
  int i;

  //alloc and initialize bools vector
  bools = gsl_vector_alloc(maxSampleSize + 1);
  gsl_vector_set_zero(bools);

  //go through data, set each sampleSize value to 1
  for(i = 0; i < snpNumber; i++){
    gsl_vector_set(bools,data[i].n, 1);
  }
  return(bools);
}


/* sampleSize vector - returns a vector of 1s or 0s depending on sample sizes in the data */
void sampleSizeVectorFill(struct snp data[], int snpNumber, int maxSampleSize, gsl_vector *bools){
  int i;

  //initialize bools vector
  gsl_vector_set_zero(bools);

  //go through data, set each sampleSize value to 1
  for(i = 0; i < snpNumber; i++){
    gsl_vector_set(bools,data[i].n, 1);
  }
}
/* returns a vector of snpProbDenominators from i = 2 to n, 1 indexed */
gsl_vector *makeSnpProbDenomVector(double beta, int maxSampleSize, gsl_vector *sampleSizeVector){
	gsl_vector *probs;
	int i, j, test;
	double tot;
	struct my_F_params tempParams;
	gsl_function F;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	probs = gsl_vector_alloc(maxSampleSize + 1);

	tempParams.w = w;
/* go from n = 2 to maxSampleSize, calc total prob space for each, put in vector probs */
	for(j = 2; j < maxSampleSize + 1; j++){
	//check if sampleSize is included
		test = gsl_vector_get(sampleSizeVector,j);
		if (test){
			tot = 0.0;
			tempParams.n = j;
			for(i = 1; i < j ; i++){
				tempParams.i = i;
				F.function =  &my_FAdapt;
				F.params = &tempParams;
				tot += GSL_FN_EVAL(&F, beta);
			}
			gsl_vector_set(probs, j, tot);
		}
	}
	gsl_integration_workspace_free(w);
	return probs;
}
  
/* snpProbDenom - returns denom for snpProb */
double snpProbDenom(double beta, int sampleSize){
  int i;
  double tot;
  struct my_F_params tempParams;
  gsl_function F;

  tot = 0.0;
  tempParams.n = sampleSize;
  for(i = 1; i < sampleSize ; i++){
	tempParams.i = i;
	F.function =  &my_F;
	F.params = &tempParams;
	tot += GSL_FN_EVAL(&F, beta);
      }
  return(tot);
}


/* snpProbMatrix- returns a matrix of log snpProbs. matrix is maxSampleSize+1 by maxSampleSize. rows represent
   variable sampleSizes (2 to max), columns represent freqs (1 to max -1) */
gsl_matrix *snpProbMatrix(double beta,  int maxSampleSize, gsl_vector *sampleSizeVector, gsl_matrix *sfsBools){
  gsl_vector *denoms;
  gsl_matrix *probs;
  struct my_snpProb_params FParams;
  gsl_function f;
  int i, j;
  
  //alloc probs
  probs = gsl_matrix_alloc(maxSampleSize + 1, maxSampleSize);
  gsl_matrix_set_zero(probs);
  //get denoms
  denoms = makeSnpProbDenomVector(beta,maxSampleSize, sampleSizeVector);
  FParams.w = gsl_integration_workspace_alloc(1000);
  //go through sampleSizes
  for(i = 2; i <= maxSampleSize; i++){
    //have sample size i?
    if (gsl_vector_get(sampleSizeVector, i)){
      //go through freqs
      for(j = 1; j < maxSampleSize; j++){
	//have freq j?
	if(gsl_matrix_get(sfsBools, i, j)){
	  //calc prob
	  FParams.i = j;
	  FParams.n = i;
	  FParams.denom = gsl_vector_get(denoms,i);
	  f.function = &snpProbDenomLookup;
	  f.params = &FParams;
	  gsl_matrix_set(probs, i, j, log(GSL_FN_EVAL(&f, beta)));
	}
      }
    }
  }
  gsl_vector_free(denoms);
  gsl_integration_workspace_free(FParams.w);
  return(probs);
}
/* snpProbMatrix- returns a matrix of  snpProbs. matrix is maxSampleSize+1 by maxSampleSize. rows represent
variable sampleSizes (2 to max), columns represent freqs (1 to max -1) */
gsl_matrix *snpProbMatrixNotLog(double beta,  int maxSampleSize, gsl_vector *sampleSizeVector, gsl_matrix *sfsBools){
	gsl_vector *denoms;
	gsl_matrix *probs;
	struct my_snpProb_params FParams;
	gsl_function f;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (100000);
	int i, j;

	FParams.w = w;
//alloc probs
	probs = gsl_matrix_alloc(maxSampleSize + 1, maxSampleSize);
	gsl_matrix_set_zero(probs); 
//get denoms
	denoms = makeSnpProbDenomVector(beta,maxSampleSize, sampleSizeVector);

//go through sampleSizes
	for(i = 2; i <= maxSampleSize; i++){
	//have sample size i?
		if (gsl_vector_get(sampleSizeVector, i)){
	//go through freqs
			for(j = 1; j < maxSampleSize; j++){
	//have freq j?
				if(gsl_matrix_get(sfsBools, i, j)){
	//calc prob
					FParams.i = j;
					FParams.n = i;
					FParams.denom = gsl_vector_get(denoms,i);
					f.function = &snpProbDenomLookup;
					f.params = &FParams;
					gsl_matrix_set(probs, i, j, GSL_FN_EVAL(&f, beta));
				}
			}
		}
	}
	gsl_vector_free(denoms);
	return(probs);
}
/* snpProbMatrixNotLogFull- returns a matrix of  snpProbs. matrix is maxSampleSize+1 by maxSampleSize. rows represent
   variable sampleSizes (2 to max), columns represent freqs (1 to max -1). All freqs/samplesizes included */
gsl_matrix *snpProbMatrixNotLogFull(double beta,  int maxSampleSize, gsl_vector *sampleSizeVector){
  gsl_vector *denoms;
  gsl_matrix *probs;
  struct my_snpProb_params FParams;
  gsl_function f;
  int i, j;
  
  //alloc probs
  probs = gsl_matrix_alloc(maxSampleSize + 1, maxSampleSize);
  gsl_matrix_set_zero(probs); 

  //get denoms

  denoms = makeSnpProbDenomVector(beta,maxSampleSize, sampleSizeVector);

  //go through sampleSizes
  for(i = 2; i <= maxSampleSize; i++){
    //check for sampleSize
    if (gsl_vector_get(sampleSizeVector,i)){
      //go through freqs
      for(j = 1; j < maxSampleSize; j++){
	//calc prob
	FParams.i = j;
	FParams.n = i;
	FParams.denom = gsl_vector_get(denoms,i);
	f.function = &snpProbDenomLookup;
	f.params = &FParams;
	gsl_matrix_set(probs, i, j, GSL_FN_EVAL(&f, beta));
      }
    }
  }
  gsl_vector_free(denoms);
  return(probs);
}


/* snpProbVectorNotLog- returns a vector of log snpProbs at specified sampleSize columns represent freqs (1 to max -1) */
gsl_vector *snpProbVectorNotLog(double beta,  int sampleSize){
  double denom;
  gsl_vector *probs;
  struct my_snpProb_params FParams;
  gsl_function f;
  int i;
  
  //alloc probs
  probs = gsl_vector_alloc(sampleSize);
  gsl_vector_set_zero(probs);
  //get denoms
  denom = snpProbDenom(beta, sampleSize);

  //go through freqs
  for(i = 1; i < sampleSize; i++){
    FParams.i = i;
    FParams.n = sampleSize;
    FParams.denom = denom;
    f.function = &snpProbDenomLookup;
    f.params = &FParams;
    gsl_vector_set(probs, i, GSL_FN_EVAL(&f, beta));
  }
  return(probs);
}

gsl_vector *snpProbVectorNotLogPre(double beta,  int sampleSize, gsl_vector *probs){
  double denom;
  struct my_snpProb_params FParams;
  gsl_function f;
  int i;
  
   gsl_vector_set_zero(probs);
  //get denoms
  denom = snpProbDenom(beta, sampleSize);

  //go through freqs
  for(i = 1; i < sampleSize; i++){
    FParams.i = i;
    FParams.n = sampleSize;
    FParams.denom = denom;
    f.function = &snpProbDenomLookup;
    f.params = &FParams;
    gsl_vector_set(probs, i, GSL_FN_EVAL(&f, beta));
  }
  return(probs);
}

/* likelihood function for SFS given the data; p here contains the data and the weights */
double my_lik(double beta, void * p){
	struct my_lik_params * params = (struct my_lik_params *) p;
	struct my_F_params FParams;
	gsl_function f;
	double lik;
	int i;
	
	lik = 0;
	for(i = 0; i < params->snpNumber; i++){
		FParams.i = params->data[i].i;
		FParams.n = params->data[i].n;
		f.function = &snpProb;
		f.params = &FParams;
		lik += log(GSL_FN_EVAL(&f, beta));
		}
	return -lik;
}

/* likelihood function for SFS where each SNP has independent selection coeff */
double sfsLikBetaVector(gsl_vector *betas, void * p){
	struct my_lik_params * params = (struct my_lik_params *) p;
	struct my_F_params FParams;
	gsl_function f;
	double lik;
	int i;
	
	lik = 0;
	for(i = 0; i < params->snpNumber; i++){
		FParams.i = params->data[i].i;
		FParams.n = params->data[i].n;
		FParams.w = params->w;
		f.function = &snpProb;
		f.params = &FParams;
		//	printf("here lik:%f beta:%f i: %d n: %d\n",lik,gsl_vector_get(betas, i),FParams.i,FParams.n);
		
		lik += log(GSL_FN_EVAL(&f, gsl_vector_get(betas, i)));
		}
	return -lik;
}

/*likelihood from Bustamante et al. 2001 */
double bustaLik(double beta, void * p){
	struct my_lik_params * params = (struct my_lik_params *) p;
	struct my_F_params FParams;
	gsl_function f;
	double lik;
	int i;

	lik = 0;
	for(i = 0; i < params->snpNumber; i++){
		FParams.i = params->data[i].i;
		FParams.n = params->data[i].n;
		f.function = &snpProb;
		f.params = &FParams;
		lik += log(GSL_FN_EVAL(&f, beta));
		}
	return -lik;
}

/* weighted likelihood function for SFS given the data; here contains the data and the weights */
double weightedLik(double beta, void * p){
  struct my_lik_params * params = (struct my_lik_params *) p;
  struct my_F_params FParams;
  gsl_function f;
  double lik;
  int i;
	
  lik = 0;
  for(i = 0; i < params->snpNumber; i++){
    FParams.i = params->data[i].i;
    FParams.n = params->data[i].n;
    f.function = &snpProb;
    f.params = &FParams;
    lik += gsl_vector_get(params->weights, i) * log(GSL_FN_EVAL(&f, beta)); 
  }
  return -lik;
}

/* weighted likelihood function for SFS given the data; here contains the data and the weights and uses the snpProbLookup routine */
double weightedLikLook(double beta, void * p){
  struct my_lik_params * params = (struct my_lik_params *) p;
  double lik;
  int i;
  gsl_matrix *probs;
  
  //make prob matrix (note this is in logs)
  probs = snpProbMatrix(beta, params->maxSampleSize, params->sampleSizeVector, params->sfsBools);
  //now go and tally likelihood
  lik = 0;
  for(i = 0; i < params->snpNumber; i++){
    lik += gsl_matrix_get(probs,params->data[i].n,params->data[i].i)  * gsl_vector_get(params->weights,i);
  }
  gsl_matrix_free(probs);
  return -lik;
}
 
/* Folded SFS version of above; weighted likelihood function for SFS given the data; here contains the data and the weights and uses the snpProbLookup routine */
double weightedLikLookFolded(double beta, void * p){
  struct my_lik_params * params = (struct my_lik_params *) p;
  double lik;
  int i;
  gsl_matrix *probs;
  
  //make prob matrix (note this is in logs)
  probs = snpProbMatrix(beta, params->maxSampleSize, params->sampleSizeVector, params->sfsBools);
  //now go and tally likelihood
  lik = 0;
  for(i = 0; i < params->snpNumber; i++){
    lik += gsl_matrix_get(probs,params->data[i].n,params->data[i].i)  * gsl_vector_get(params->weights,i);
	lik += gsl_matrix_get(probs,params->data[i].n, params->data[i].n - params->data[i].i)  * gsl_vector_get(params->weights,i);
  }
  gsl_matrix_free(probs);
  return -lik;
}

/* PRF log Likelihood  */
double prfLik(double beta, void * p){
	 struct prfLik_params * params = (struct prfLik_params *) p;
	double lik;
	double term1, term2,theta;
  	int i;	
	struct my_F_params tempParams;
	gsl_function F;
	double fTerm = 0.;
	
	theta = gsl_vector_get(params->betas,1);
	term1 = term2 = 0;
	for(i = 1; i < params->maxSampleSize;i++){
		tempParams.i = i;
		tempParams.n = params->maxSampleSize;
		F.function = &my_F;
		F.params = &tempParams;
		fTerm = GSL_FN_EVAL(&F,beta);
		term1 += gsl_matrix_get(params->sfsCounts,params->maxSampleSize,i) * log( theta* fTerm) ;
		term2 += fTerm;
	}
	return((term1 - (theta * term2)) * -1);
}	
  



double weightedLikLookCI(double beta, void * p){
  struct my_likCI_params *params = (struct my_likCI_params *) p;
  struct my_lik_params likParams;
  gsl_function l;
  double result;

  likParams.data = params->data;
  likParams.snpNumber = params->snpNumber;
  likParams.sampleSizeVector = params->sampleSizeVector;
  likParams.maxSampleSize = params->maxSampleSize;
  likParams.weights = params->weights;
  likParams.sfsBools = params->sfsBools;
  l.function = &weightedLikLook;
  l.params = &likParams;
  result = GSL_FN_EVAL(&l, beta) - (params->lMax) - (params->logUnits) ;
  //  printf("%f\t%f\t%f\n",params->lMax,GSL_FN_EVAL(&l, beta),result);
  return result;
}


double weightedLikLookCIFolded(double beta, void * p){
  struct my_likCI_params *params = (struct my_likCI_params *) p;
  struct my_lik_params likParams;
  gsl_function l;
  double result;

  likParams.data = params->data;
  likParams.snpNumber = params->snpNumber;
  likParams.sampleSizeVector = params->sampleSizeVector;
  likParams.maxSampleSize = params->maxSampleSize;
  likParams.weights = params->weights;
  likParams.sfsBools = params->sfsBools;
  l.function = &weightedLikLookFolded;
  l.params = &likParams;
  result = GSL_FN_EVAL(&l, beta) - (params->lMax) - (params->logUnits) ;
  //  printf("%f\t%f\t%f\n",params->lMax,GSL_FN_EVAL(&l, beta),result);
  return result;
}


/* wrapper for my_lik for use in mnbrak */
double likWrap(double beta){
	gsl_function f;
	
	f.function = &my_lik;
	f.params = NULL;
	return GSL_FN_EVAL(&f, beta);
}

/* wrapper for my_lik for use in mnbrak */
double  wlikWrap(double beta, void * p){
	gsl_function f;
	struct my_lik_params * params = (struct my_lik_params *) p;
	
	f.function = &weightedLik;
	f.params = params;
	return GSL_FN_EVAL(&f, beta);
}

/* wrapper for my_lik for use in mnbrak */
double  wlikWrapLook(double beta, void * p){
	gsl_function f;
	struct my_lik_params * params = (struct my_lik_params *) p;
	
	f.function = &weightedLikLook;
	f.params = params;
	return GSL_FN_EVAL(&f, beta);
}

double  wlikWrapLookFolded(double beta, void * p){
	gsl_function f;
	struct my_lik_params * params = (struct my_lik_params *) p;
	
	f.function = &weightedLikLookFolded;
	f.params = params;
	return GSL_FN_EVAL(&f, beta);
}

double  wlikCIWrapLook(double beta, void * p){
	gsl_function f;
	struct my_likCI_params * params = (struct my_likCI_params *) p;
	
	f.function = &weightedLikLookCI;
	f.params = params;
	return GSL_FN_EVAL(&f, beta);
}
 

/* wrapper for prfLik for use in mnbrak */
double  prfLikWrap(double beta, void * p){
	gsl_function f;
	struct prfLik_params * params = (struct prfLik_params *) p;
	
	f.function = &prfLik;
	f.params = params;
	return GSL_FN_EVAL(&f, beta);
} 

/*error checking function for likelihood function */
void printProfile(double min, double max, void * p){
  gsl_function f;
  struct my_lik_params *params =  (struct my_lik_params *) p;
  int i;

  f.function = &weightedLikLook;
  f.params = params;
  for(i = min; i < max; i++){
    printf("%d\t%f\n",i, GSL_FN_EVAL(&f, i));
  }
 
}
    
/* get ML estimate of beta using likelihood function; BRENT method for gsl */
double ml_est(double * lik_beta_hat){
	int status;
	int iter = 0;
	int max_iter = 100;
	const gsl_min_fminimizer_type *T;
	gsl_min_fminimizer *s;
	double m, a, b, ax, bx, cx, fa, fb, fc;
	gsl_function L;
	
	
	
	/* get minimum bracket */
	ax = -10.0;
	bx = 0.0;
	cx = 10.0;
	
	mnbrak(&ax, &bx, &cx, &fa, &fb, &fc, &my_lik);
	
	/* check if there is a minimum, if not return inf */
	if (fa == fb || fc == fb || ax > cx){
		*lik_beta_hat = fa;
		return 666.0;
		}
	else{
		
		/* do mimization */
		
		L.function = &my_lik;
		//should have param line here but no params....
		
		T = gsl_min_fminimizer_brent;
		s = gsl_min_fminimizer_alloc(T);
		gsl_min_fminimizer_set(s, &L, bx, ax, cx);
		
		do{
			iter++;
			status = gsl_min_fminimizer_iterate(s);
			 m = gsl_min_fminimizer_x_minimum (s);
			  a = gsl_min_fminimizer_x_lower (s);
			  b = gsl_min_fminimizer_x_upper (s);
		
			  status = gsl_min_test_interval (a, b, 0.001, 0.0);
		
		}
	  while (status == GSL_CONTINUE && iter < max_iter);
	  gsl_min_fminimizer_free(s);
	  
	  return m;
	}
}

/* get weighted ML estimate of beta using likelihood function; BRENT method for gsl */
double weighted_ml_est(double * lik_beta_hat, void * p){
	int status;
	int iter = 0;
	int max_iter = 100;
	const gsl_min_fminimizer_type *T;
	gsl_min_fminimizer *s;
	double m, a, b, ax, bx, cx, fa, fb, fc, min;
	gsl_function L;
	struct my_lik_params * params = (struct my_lik_params *) p;	
	
	
	// get minimum bracket 
	ax = -10.0;
	bx = 1.0;
	cx = 10;
	L.function = &weightedLik;
	L.params = params;
	mnbrak2(&ax, &bx, &cx, &fa, &fb, &fc, &wlikWrap, bx, params);

	//check if there is a minimum, if not return 666
	if (fa == fb || fc == fb || ax > cx){
	  min = fa;
	  if(fb < min){
	    min = fb;
	  }
	  if(fc < min){
	    min = fc;
	  }
	  *lik_beta_hat = min;
 	
	  return 666.0;
	}
	else{
	  // do mimization 
		
	  //initialize lik function
	  L.function = &weightedLik;
	  L.params = params;
	  
	  //min routine
	  T = gsl_min_fminimizer_brent;
	  s = gsl_min_fminimizer_alloc(T);
	  gsl_min_fminimizer_set(s, &L, bx, ax, cx);
	  do{
	    iter++;
	    status = gsl_min_fminimizer_iterate(s);
	    m = gsl_min_fminimizer_x_minimum (s);
	    a = gsl_min_fminimizer_x_lower (s);
	    b = gsl_min_fminimizer_x_upper (s);
	    status = gsl_min_test_interval (a, b, 0.001, 0.0);
	  }
	  while (status == GSL_CONTINUE && iter < max_iter);
	  *lik_beta_hat = gsl_min_fminimizer_f_minimum(s);
	  gsl_min_fminimizer_free(s);
	  return m;
	}
}

double weighted_ml_est_lookup(double * lik_beta_hat, void * p){
	int status;
	int iter = 0;
	int max_iter = 1000;
	const gsl_min_fminimizer_type *T;
	gsl_min_fminimizer *s;
	double m, a, b, ax, bx, cx, fa, fb, fc, dummy;
	gsl_function L;
	struct my_lik_params * params = (struct my_lik_params *) p;	
	
	// get minimum bracket 
	ax = -20.0;
	bx = -1.0;
	cx = 20000;
	L.function = &weightedLikLook;
	L.params = params;
	mnbrak2(&ax, &bx, &cx, &fa, &fb, &fc, &wlikWrapLook, bx, params);
	//swap bounds if needed
	if (cx < ax){
	  dummy = cx;
	  cx = ax;
	  ax = dummy;
	}

	if (fb == fc){
		return HUGE_VAL;
	}
	// do mimization 
		
	//initialize lik function
	L.function = &weightedLikLook;
	L.params = params;
	  
	//min routine
	T = gsl_min_fminimizer_brent;
	s = gsl_min_fminimizer_alloc(T);
	gsl_min_fminimizer_set(s, &L, bx, ax, cx);
	do{
	  iter++;
	  status = gsl_min_fminimizer_iterate(s);
	  m = gsl_min_fminimizer_x_minimum (s);
	  a = gsl_min_fminimizer_x_lower (s);
	  b = gsl_min_fminimizer_x_upper (s);
	  status = gsl_min_test_interval (a, b, 0.01, 0.01);
	}
	while (status == GSL_CONTINUE && iter < max_iter);
	*lik_beta_hat = gsl_min_fminimizer_f_minimum(s);
	gsl_min_fminimizer_free(s);
	return m;
}

/* same as above but outputs 666 on minimization error */
double weighted_ml_est_lookup_errReport(double * lik_beta_hat, void * p){
	int status;
	int iter = 0;
	int max_iter = 100;
	gsl_min_fminimizer *s;
	double m, a, b, ax, bx, cx, fa, fb, fc, dummy;
	gsl_function L;
	struct my_lik_params * params = (struct my_lik_params *) p;	
	
	
	// get minimum bracket 
	ax = -20.0;
	bx = -1.0;
	cx = 20;
	L.function = &weightedLikLook;
	L.params = params;
	mnbrak2(&ax, &bx, &cx, &fa, &fb, &fc, &wlikWrapLook, bx, params);
	//swap bounds if needed
	if (cx < ax){
	  dummy = cx;
	  cx = ax;
	  ax = dummy;
	}
	//check for error
	if (fb >= fc || fb >= fa){
	  return(666.0);
	}
	
	// do mimization 
		
	//initialize lik function
	L.function = &weightedLikLook;
	L.params = params;
	  
	//min routine
	s = params->min;
	gsl_set_error_handler_off ();
	status = gsl_min_fminimizer_set(s, &L, bx, ax, cx);
	if (status == GSL_EINVAL){
	  return(666);
	}
	do{
	  iter++;
	  status = gsl_min_fminimizer_iterate(s);
	  m = gsl_min_fminimizer_x_minimum (s);
	  a = gsl_min_fminimizer_x_lower (s);
	  b = gsl_min_fminimizer_x_upper (s);
	  status = gsl_min_test_interval (a, b, 0.01, 0.01);
	}
	while (status == GSL_CONTINUE && iter < max_iter);
	*lik_beta_hat = gsl_min_fminimizer_f_minimum(s);
	//gsl_min_fminimizer_free(s);
	return m;
}

/* get weighted ML lower CI beta using weightedLikLookCI method; BRENT method for gsl */
double weighted_ml_CILower_lookup(double * lik_beta_hat, void * p){
	int status;
	int iter = 0;
	int max_iter = 100;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	double m, a, b, ax, cx;
	gsl_function L;
	struct my_likCI_params * params = (struct my_likCI_params *) p;	
	
	
	// get minimum bracket 
	ax = -50.0;
	cx = params->beta_hat;
	 
	// do root finding 		
	//initialize lik function
	L.function = &weightedLikLookCI;
	L.params = params;
	  
	//min routine
	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc(T);
	gsl_root_fsolver_set(s, &L, ax, cx);
	do{
	  iter++;
	  status = gsl_root_fsolver_iterate(s);
	  m = gsl_root_fsolver_root(s);
	  a = gsl_root_fsolver_x_lower(s);
	  b = gsl_root_fsolver_x_upper(s);
	  status = gsl_min_test_interval (a, b, 0.001, 0.0);
	}
	while (status == GSL_CONTINUE && iter < max_iter);
	gsl_root_fsolver_free(s);
	return m;
}

/* get weighted ML upper CI beta using weightedLikLookCI method; BRENT method for gsl */
double weighted_ml_CIUpper_lookup(double * lik_beta_hat, void * p){
	int status;
	int iter = 0;
	int max_iter = 100;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	double m, a, b, ax, cx;
	gsl_function L;
	struct my_likCI_params * params = (struct my_likCI_params *) p;	
	
	
	// get minimum bracket 
	cx = 100.0;
	ax = params->beta_hat;
	 
	// do root finding 		
	//initialize lik function
	L.function = &weightedLikLookCI;
	L.params = params;
	  
	//min routine
	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc(T);
	gsl_root_fsolver_set(s, &L, ax, cx);
	do{
	  iter++;
	  status = gsl_root_fsolver_iterate(s);
	  m = gsl_root_fsolver_root(s);
	  a = gsl_root_fsolver_x_lower(s);
	  b = gsl_root_fsolver_x_upper(s);
	  status = gsl_min_test_interval (a, b, 0.001, 0.0);
	}
	while (status == GSL_CONTINUE && iter < max_iter);
	gsl_root_fsolver_free(s);
	return m;
}

//Folded Version
double weighted_ml_est_lookup_folded(double * lik_beta_hat, void * p){
	int status;
	int iter = 0;
	int max_iter = 1000;
	const gsl_min_fminimizer_type *T;
	gsl_min_fminimizer *s;
	double m, a, b, ax, bx, cx, fa, fb, fc, dummy;
	gsl_function L;
	struct my_lik_params * params = (struct my_lik_params *) p;	
	
	// get minimum bracket 
	ax = -20.0;
	bx = -1.0;
	cx = 20000;
	L.function = &weightedLikLook;
	L.params = params;
	mnbrak2(&ax, &bx, &cx, &fa, &fb, &fc, &wlikWrapLookFolded, bx, params);
	//swap bounds if needed
	if (cx < ax){
	  dummy = cx;
	  cx = ax;
	  ax = dummy;
	}

	if (fb == fc){
		return HUGE_VAL;
	}
	// do mimization 
		
	//initialize lik function
	L.function = &weightedLikLook;
	L.params = params;
	  
	//min routine
	T = gsl_min_fminimizer_brent;
	s = gsl_min_fminimizer_alloc(T);
	gsl_min_fminimizer_set(s, &L, bx, ax, cx);
	do{
	  iter++;
	  status = gsl_min_fminimizer_iterate(s);
	  m = gsl_min_fminimizer_x_minimum (s);
	  a = gsl_min_fminimizer_x_lower (s);
	  b = gsl_min_fminimizer_x_upper (s);
	  status = gsl_min_test_interval (a, b, 0.01, 0.01);
	}
	while (status == GSL_CONTINUE && iter < max_iter);
	*lik_beta_hat = gsl_min_fminimizer_f_minimum(s);
	gsl_min_fminimizer_free(s);
	return m;
}

/* get weighted ML lower CI beta using weightedLikLookCI method; BRENT method for gsl */
double weighted_ml_CILower_lookupFolded(double * lik_beta_hat, void * p){
	int status;
	int iter = 0;
	int max_iter = 100;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	double m, a, b, ax, cx;
	gsl_function L;
	struct my_likCI_params * params = (struct my_likCI_params *) p;	
	
	
	// get minimum bracket 
	ax = -50.0;
	cx = params->beta_hat;
	 
	// do root finding 		
	//initialize lik function
	L.function = &weightedLikLookCIFolded;
	L.params = params;
	  
	//min routine
	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc(T);
	gsl_root_fsolver_set(s, &L, ax, cx);
	do{
	  iter++;
	  status = gsl_root_fsolver_iterate(s);
	  m = gsl_root_fsolver_root(s);
	  a = gsl_root_fsolver_x_lower(s);
	  b = gsl_root_fsolver_x_upper(s);
	  status = gsl_min_test_interval (a, b, 0.001, 0.0);
	}
	while (status == GSL_CONTINUE && iter < max_iter);
	gsl_root_fsolver_free(s);
	return m;
}

/* get weighted ML upper CI beta using weightedLikLookCI method; BRENT method for gsl */
double weighted_ml_CIUpper_lookupFolded(double * lik_beta_hat, void * p){
	int status;
	int iter = 0;
	int max_iter = 100;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	double m, a, b, ax, cx;
	gsl_function L;
	struct my_likCI_params * params = (struct my_likCI_params *) p;	
	
	
	// get minimum bracket 
	cx = 100.0;
	ax = params->beta_hat;
	 
	// do root finding 		
	//initialize lik function
	L.function = &weightedLikLookCIFolded;
	L.params = params;
	  
	//min routine
	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc(T);
	gsl_root_fsolver_set(s, &L, ax, cx);
	do{
	  iter++;
	  status = gsl_root_fsolver_iterate(s);
	  m = gsl_root_fsolver_root(s);
	  a = gsl_root_fsolver_x_lower(s);
	  b = gsl_root_fsolver_x_upper(s);
	  status = gsl_min_test_interval (a, b, 0.001, 0.0);
	}
	while (status == GSL_CONTINUE && iter < max_iter);
	gsl_root_fsolver_free(s);
	return m;
}

//ML estimation using PRF likelihood model
double ml_prf(double * lik_beta_hat, void * p){
	int status;
	int iter = 0;
	int max_iter = 100;
	const gsl_min_fminimizer_type *T;
	gsl_min_fminimizer *s;
	double m, a, b, ax, bx, cx, fa, fb, fc, dummy;
	gsl_function L;
	struct prfLik_params * params = (struct prfLik_params *) p;	
	
	
	// get minimum bracket 
	ax = -10.0;
	bx = 0.0;
	cx = 10;

	mnbrak2(&ax, &bx, &cx, &fa, &fb, &fc, &prfLikWrap, bx, params);
	//swap bounds if needed
	if (cx < ax){
	  dummy = cx;
	  cx = ax;
	  ax = dummy;
	}

	if (fb == fc){
		return HUGE_VAL;
	}
	// do mimization 
		
	//initialize lik function
	L.function = &prfLik;
	L.params = params;
	  
	//min routine
	T = gsl_min_fminimizer_brent;
	s = gsl_min_fminimizer_alloc(T);
	gsl_min_fminimizer_set(s, &L, bx, ax, cx);
	do{
	  iter++;
	  status = gsl_min_fminimizer_iterate(s);
	  m = gsl_min_fminimizer_x_minimum (s);
	  a = gsl_min_fminimizer_x_lower (s);
	  b = gsl_min_fminimizer_x_upper (s);
	  status = gsl_min_test_interval (a, b, 0.001, 0.0);
	}
	while (status == GSL_CONTINUE && iter < max_iter);
	*lik_beta_hat = gsl_min_fminimizer_f_minimum(s);
	gsl_min_fminimizer_free(s);
	return m;
}

/* summarizeSFS- returns a matrix of counts which summarize the SFS conditional
   upon sampleSize. rows = sampleSize, column = freqs */

gsl_matrix *summarizeSFS(int maxSampleSize, struct snp data[], int snpNumber){
  gsl_matrix *counts;
  int i;

  //alloc matrix of size maxSampleSize by maxSampleSize
  counts = gsl_matrix_alloc(maxSampleSize + 1, maxSampleSize + 1);
  gsl_matrix_set_zero(counts);
  assert(counts != NULL);
    
  //gsl_matrix_set_zero(counts);
  
  //go through data, fill up matrix
  for(i = 0; i < snpNumber; i++){
    gsl_matrix_set(counts, data[i].n, data[i].i, gsl_matrix_get(counts, data[i].n, data[i].i) + 1);
  }
  return(counts);

}

/* summarizeSFSFill- fills a matrix of counts which summarize the SFS conditional
   upon sampleSize. rows = sampleSize, column = freqs */

void summarizeSFSFill(gsl_matrix *counts, struct snp data[], int snpNumber){
  int i;

  //fill matrix of size maxSampleSize by maxSampleSize
  gsl_matrix_set_zero(counts);
  assert(counts != NULL);
  gsl_matrix_set_zero(counts);
  
  //go through data, fill up matrix
  for(i = 0; i < snpNumber; i++){
    gsl_matrix_set(counts, data[i].n, data[i].i, gsl_matrix_get(counts, data[i].n, data[i].i) + 1);
  }
}

/* summarizeSFSFill- fills a matrix of counts which summarize the SFS conditional
   upon sampleSize. rows = sampleSize, column = freqs */

int summarizeSFSFillFromTo(gsl_matrix *counts, struct snp data[], int snpNumber, int from, int to, int startPos){
  int i;

  //fill matrix of size maxSampleSize by maxSampleSize
  gsl_matrix_set_zero(counts);
  assert(counts != NULL);
  gsl_matrix_set_zero(counts);
  
  //go through data, fill up matrix
	i = startPos;
	while(data[i].pos < to){
		if(data[i].pos >= from)
    		gsl_matrix_set(counts, data[i].n, data[i].i, gsl_matrix_get(counts, data[i].n, data[i].i) + 1);
		i++;
  	}
	return(i);
}



/* summarizeSFS- returns a matrix of booleans which summarize the SFS conditional
   upon sampleSize. rows = sampleSize, column = freqs */

gsl_matrix *summarizeSFSBool(int maxSampleSize, struct snp data[], int snpNumber){
  gsl_matrix *counts;
  int i;

  //alloc matrix of size maxSampleSize by maxSampleSize
  counts = gsl_matrix_alloc(maxSampleSize + 1, maxSampleSize + 1);
  gsl_matrix_set_zero(counts);
  assert(counts != NULL);
    
  gsl_matrix_set_zero(counts);
  
  //go through data, fill up matrix
  for(i = 0; i < snpNumber; i++){
    gsl_matrix_set(counts, data[i].n, data[i].i, 1);
  }
  return(counts);

}

/* summarizeSFSBoolFill- fills a matrix of booleans which summarize the SFS conditional
   upon sampleSize. rows = sampleSize, column = freqs */

void summarizeSFSBoolFill(gsl_matrix *counts, struct snp data[], int snpNumber){
  int i;

  //fill matrix of size maxSampleSize by maxSampleSize
  //counts = gsl_matrix_alloc(maxSampleSize + 1, maxSampleSize + 1);
  gsl_matrix_set_zero(counts);
  assert(counts != NULL);
    
  gsl_matrix_set_zero(counts);
  
  //go through data, fill up matrix
  for(i = 0; i < snpNumber; i++){
    gsl_matrix_set(counts, data[i].n, data[i].i, 1);
  }
}


/* simulateSiteFreq- returns an int representing the frequency of
   of a site (unfolded), i, in sampleSize n, conditional on sel. 
   coeff. beta , also needs a pointer to a gsl random number generator */

int simulateSiteFreq(int sampleSize, double alpha, void *r){
  int i;
  double sum, probs[sampleSize - 1], rand;
  gsl_function p;
  struct my_F_params simParams;
  gsl_rng * rn = (gsl_rng *) r;
  
  simParams.n = sampleSize;
  /* put probs in vector */
  for(i = 0; i < sampleSize - 1; i++){
    simParams.i = i + 1;
    p.function = &snpProb;
    p.params = &simParams;
    probs[i] = GSL_FN_EVAL(&p, alpha);  
  }

  /* choose frequency */
  rand = gsl_rng_uniform(rn);
  sum = 0;
  for(i = 0; i < sampleSize - 1; i++){
    sum += probs[i];
    if (rand <= sum){
      return(i+1);
    }
  }
  return(666);
}

/* ascertainment stuff -- following Nielsen et al. 2004  */

/* simpleAscertain--  Nielsen et al. 2004 eqn 2 *, returns the prob of ascertainment 
given a sample freq */
double simpleAscertain(double sampleFreq, void * p){
  struct my_F_params * params = (struct my_F_params *) p;
  double sampleSize, ascSize,tmpSampleFreq;
  double num1 = 0.0;
  double num2 = 0.0;
 
  sampleSize = params->n;
  ascSize = params->ascSize;
  tmpSampleFreq = (int) sampleFreq;

  //test for weirdness
  
  if (sampleFreq < ascSize){
    num1 = 0;
  }
  else{
    num1 = gsl_sf_choose(sampleFreq, ascSize);
  }

  if ((sampleSize -sampleFreq) < ascSize){
    num2 = 0;
  }
  else{
    num2 =  gsl_sf_choose(sampleSize - sampleFreq, ascSize);
  }

  return(1.0 - ((num1 + num2) / gsl_sf_choose(sampleSize, ascSize)));
}

/* outgroupAscertain-- equivalent to Nielsen et al. 2004 eqn 2,
 returns the prob of ascertainment given a sample freq, but
is for use when ascertainment was based on single outgroup comparisons
(e.g. divergence between reference genomes). specifically this
calculates the prob of sampling only ancestral alleles.
IMPORTANT- sampleFreq here is derived allele freq */
double outgroupAscertain(double sampleFreq, void * p){
	struct my_F_params * params = (struct my_F_params *) p;
	int sampleSize, ascSize, derived;
	double num = 0.0;

	sampleSize = params->n;
	ascSize = params->ascSize;
	derived = params->derivedFlag;
	if(derived){	    
	    num =  xchoosey( sampleFreq, ascSize);
	  }
	  else{
	    num =  xchoosey((sampleSize - sampleFreq), ascSize);
	  }
	  return(num /xchoosey(sampleSize, ascSize));
}

/*same as above for case where ascSamp is not subset of sample */
double outgroupAscertain2(double sampleFreq, void * p){
	struct my_F_params * params = (struct my_F_params *) p;
	int sampleSize, ascSize,ascFreq, derived, N, j;
	double num = 0.0, denom = 0.0;

	sampleSize = params->n;
	ascSize = params->ascSize;
	ascFreq = params->ascFreq;
	derived = params->derivedFlag;
	
	N = sampleSize + ascSize;
	j = sampleFreq + ascFreq;
  	//was ascertainment for derived?
  	if(derived == 1){
		num =  xchoosey(j,sampleFreq) * xchoosey((N-j), sampleSize - sampleFreq);
  	}
  	else{
		num =  xchoosey((N-j), (sampleSize - sampleFreq)) * xchoosey(j,  sampleFreq); 
  	}
//	printf("num: %f\n", num );
//	printf("N: %d j: %d sampleSize: %d ascSize: %d derivedFlag: %d\n",N,j, sampleSize,ascSize,derived);
	denom = xchoosey(N, sampleSize);
	if (denom == 0 || num == 0){
		return(0);
	}
//	printf("oaVal: %f\n",num)/denom;
  	return(num / denom);
}

/* probAscertainmentGivenModel-- Nielsen et al. 2005 eqn 8
   returns the prob of ascertainment given the model */
double probAscertainmentGivenModel(double beta, void *p){
  struct my_F_params * params = (struct my_F_params *) p;
  struct my_F_params tempParams;
 
  long double tot, tmpres, tmpres2;
  int i=1;
  gsl_function F, G;
  
  /* add up total prob space */
  tot = 0.0;
  tempParams.n = params->n;
  tempParams.ascSize = params->ascSize;
  
  for(i = 1; i < params->n ; i++){
    tempParams.i = i;
  
    F.function =  &(snpProb);
    F.params = &tempParams;
    tmpres = 0.0;
    tmpres = GSL_FN_EVAL(&F, beta);
    G.function = &simpleAscertain;
    G.params = &tempParams;
    tmpres2 = GSL_FN_EVAL(&G, i) ;
    tot += tmpres * tmpres2;
    
  }
  
  
  return(tot) ;
}

/* probAscertainmentGivenModel-- Nielsen et al. 2005 eqn 8
   returns the prob of ascertainment given the model */
double probOutgroupAscertainmentGivenModel(double beta, void *p){
  struct my_F_params * params = (struct my_F_params *) p;
  struct my_F_params tempParams;
 
  long double tot, tmpres, tmpres2;
  int i=1;
  gsl_function F, G;
  
  /* add up total prob space */
  tot = 0.0;
  tempParams.n = params->n;
  tempParams.ascSize = params->ascSize;
  
  for(i = 1; i < params->n ; i++){
    tempParams.i = i;
  
    F.function =  &(snpProb);
    F.params = &tempParams;
    tmpres = 0.0;
    tmpres = GSL_FN_EVAL(&F, beta);
    G.function = &outgroupAscertain;
    G.params = &tempParams;
    tmpres2 = GSL_FN_EVAL(&G, i) ;
    tot += tmpres * tmpres2;
    
  }
  
  
  return(tot) ;
}

/* probAscertainmentGivenModelLookup-- Nielsen et al. 2005 eqn 8
   returns the prob of ascertainment given the model, lookup version
   requires snpProb matrix and ascProb matrix */
double probAscertainmentGivenModelLookup(double beta, int sampSize, gsl_matrix *snpProbs, gsl_matrix *ascProbs){
 
  long double tot;
  int i=1;
   
  /* add up total prob space */
  tot = 0.0;
  
  for(i = 1; i < sampSize ; i++){
    tot += gsl_matrix_get(snpProbs, sampSize, i) * gsl_matrix_get(ascProbs, sampSize, i);
  }
  return(tot) ;
}

/* probAscertainmentGivenModelHemiLookup-- same
   as above but expects vector not matrix of snpProbs */
double probAscertainmentGivenModelHemiLookup(double beta, int sampSize, gsl_vector *snpProbs, gsl_matrix *ascProbs){
 
  long double tot;
  int i=1;
   
  /* add up total prob space */
  tot = 0.0;
  
  for(i = 1; i < sampSize ; i++){
    tot += gsl_vector_get(snpProbs, i) * gsl_matrix_get(ascProbs, sampSize, i);
  }
  return(tot) ;
}


/* ascSize vector - returns a vector of 1s or 0s depending on ascertainment sample sizes in the data */
gsl_vector *ascSizeVector(struct snp data[], int snpNumber, int maxSampleSize){
  gsl_vector *bools;
  int i;

  //alloc and initialize bools vector
  bools = gsl_vector_alloc(maxSampleSize + 1);
  gsl_vector_set_zero(bools);

  //go through data, set each sampleSize value to 1
  for(i = 0; i < snpNumber; i++){
    gsl_vector_set(bools,data[i].ascSize, 1);
  }
  return(bools);
}

/* snpAscMatrix- returns a matrix of simpleAscertainment probs. matrix is maxSampleSize+1 by maxSampleSize. rows represent
   variable sampleSizes (2 to max), columns represent freqs (1 to max -1). Currently this only supports a single ascertainment size  */
gsl_matrix *snpAscMatrix(int maxSampleSize, gsl_vector *sampleSizeVector, gsl_matrix *sfsBools, \
		 int ascSize){
  gsl_matrix *probs;
  struct my_F_params FParams;
  gsl_function f;
  int i, j;
  
  //alloc probs
  probs = gsl_matrix_alloc(maxSampleSize + 1, maxSampleSize);
  gsl_matrix_set_zero(probs);
  //go through sampleSizes
  for(i = 2; i <= maxSampleSize; i++){
    //have sample size i?
    if (gsl_vector_get(sampleSizeVector, i)){
      //go through freqs
      for(j = 1; j < maxSampleSize; j++){
	//have freq j?
	if(gsl_matrix_get(sfsBools, i, j)){
	  //calc prob
	  FParams.i = j;
	  FParams.n = i;
	  FParams.ascSize = ascSize;
	  f.function = &simpleAscertain;
	  f.params = &FParams;
	  gsl_matrix_set(probs, i, j, GSL_FN_EVAL(&f, j));
	}
      }
    }
  }
  return(probs);
}

/* snpOutgroupAscMatrix- returns a matrix of outgroupAscertainment probs. 
matrix is maxSampleSize+1 by maxSampleSize. rows represent
   variable sampleSizes (2 to max), columns represent freqs (1 to max -1). 
Currently this only supports a single ascertainment size. Also has inelegance
of calculating derived or anscertral matrix  */
gsl_matrix *snpOutgroupAscMatrix(int maxSampleSize, gsl_vector *sampleSizeVector, int ascSize, int derived, \
int runMode){
	gsl_matrix *probs;
	struct my_F_params FParams;
	gsl_function f;
	int i, j;

//alloc probs
	probs = gsl_matrix_alloc(maxSampleSize + 1, maxSampleSize);
	gsl_matrix_set_zero(probs);
//go through sampleSizes
	for(i = 2; i <= maxSampleSize; i++){
	//have sample size i?
		if (gsl_vector_get(sampleSizeVector, i)){
	//go through freqs
			for(j = 1; j < i; j++){
	//calc prob
				FParams.i = j;
				FParams.n = i;
				FParams.ascSize = ascSize;
				FParams.ascFreq = derived;
				FParams.derivedFlag = derived;
				if(runMode == 0){
					f.function = &outgroupAscertain;
				}
				else{
					f.function = &outgroupAscertain2;
				}
				f.params = &FParams;
				gsl_matrix_set(probs, i, j, GSL_FN_EVAL(&f, j));
			}
		}
	}
	return(probs);
}

gsl_matrix *snpOutgroupAscMatrixPre(int maxSampleSize, gsl_vector *sampleSizeVector, int ascSize, int derived, \
int runMode, gsl_matrix *probs){
	struct my_F_params FParams;
	gsl_function f;
	int i, j;

	gsl_matrix_set_zero(probs);
//go through sampleSizes
	for(i = 2; i <= maxSampleSize; i++){
	//have sample size i?
		if (gsl_vector_get(sampleSizeVector, i)){
	//go through freqs
			for(j = 1; j < i; j++){
	//calc prob
				FParams.i = j;
				FParams.n = i;
				FParams.ascSize = ascSize;
				FParams.ascFreq = derived;
				FParams.derivedFlag = derived;
				if(runMode == 0){
					f.function = &outgroupAscertain;
				}
				else{
					f.function = &outgroupAscertain2;
				}
				f.params = &FParams;
				gsl_matrix_set(probs, i, j, GSL_FN_EVAL(&f, j));
			}
		}
	}
	return(probs);
}


/*estimateAscSFS-- ML estimates of SFS, given ascertainment. returns 
 a matrix of probs. matrix is maxSampleSize+1 by maxSampleSize. rows represent
variable sampleSizes (2 to max), columns represent freqs (1 to max -1). Currently this only supports a single ascertainment size */

gsl_matrix *estimateAscSFS(int maxSampleSize, gsl_vector *sampleSizeVector, gsl_matrix *sfsBools, \
			   int ascSize, gsl_matrix *sfsSummary){
  gsl_matrix *probs, *ascProbs;
    int i, j;
  double ssTmp, p_hat; 
  double xij;

  //alloc probs
  probs = gsl_matrix_alloc(maxSampleSize + 1, maxSampleSize);
  gsl_matrix_set_zero(probs);
  //calculate ascProbs
  ascProbs = snpAscMatrix(maxSampleSize, sampleSizeVector, sfsBools, \
			  ascSize);
  

  //go through sampleSizes 
  for(i = 2; i <= maxSampleSize; i++){
    //have sample size i?
    if (gsl_vector_get(sampleSizeVector, i)){
      //go through freqs to tally up denom at sampleSize
      ssTmp = 0.0;
      for(j = 1; j < maxSampleSize; j++){
	xij = gsl_matrix_get(sfsSummary, i, j);
	if (xij > 0){
	  ssTmp += xij / gsl_matrix_get(ascProbs, i, j);
	}
      }
      //go back through freqs to assign p_hats at sampleSize
      for(j = 1; j < maxSampleSize; j++){
	xij = gsl_matrix_get(sfsSummary, i, j);
	p_hat = (double) xij / gsl_matrix_get(ascProbs, i, j);
	gsl_matrix_set(probs, i, j, p_hat / ssTmp);
      }
    }
  }
  return(probs);
}

/* weighted likelihood function for SFS given the data and ascertainment; here contains the data and the weights and uses the snpProbLookup routine */
double weightedLikLookAsc(double beta, void * p){
  struct my_lik_params * params = (struct my_lik_params *) p;
  double lik, num;
  int i;
  gsl_matrix *probs, *ascProbs;
  gsl_vector *ascDenoms;

  
  //make prob matrix (note this is in logs)
  probs = snpProbMatrix(beta, params->maxSampleSize, params->sampleSizeVector, params->sfsBools);
  //make ascProb matrix (NOT in logs!)
  ascProbs = snpAscMatrix(params->maxSampleSize, params->sampleSizeVector, params->sfsBools, params->ascSize);

  //make vector of denoms
  ascDenoms = gsl_vector_alloc(params->maxSampleSize + 1);
  gsl_vector_set_zero(ascDenoms);
  for(i = params->ascSize; i <= params->maxSampleSize; i++){
     gsl_vector_set(ascDenoms,i,probAscertainmentGivenModelLookup(beta, i, probs, ascProbs));
  }

  //now go and tally likelihood
  lik = 0;
  num = 0;
  for(i = 0; i < params->snpNumber; i++){
    lik += gsl_matrix_get(probs,params->data[i].n,params->data[i].i)	\
      + log(gsl_matrix_get(ascProbs,params->data[i].n,params->data[i].i)) \
      - log(gsl_vector_get(ascDenoms, params->data[i].n));
     }
  gsl_matrix_free(probs);
  gsl_matrix_free(ascProbs);
  gsl_vector_free(ascDenoms);
  
  return -lik;
}
 

/* weighted likelihood function for SFS given the data and outgroup ascertainment; here contains the data and the weights and uses the snpProbLookup routine */
double weightedLikLookOutgroupAsc(double beta, void * p){
	struct my_lik_params * params = (struct my_lik_params *) p;
	double lik, num;
	int i;
	gsl_matrix *probs, *ascProbsDerived, *ascProbsAncestral;
	gsl_vector *ascDenomsDerived, *ascDenomsAncestral;

	//make prob matrix 
	probs = snpProbMatrixNotLogFull(beta, params->maxSampleSize, params->sampleSizeVector);
	//make ascProb matrices, one for derived one for ancestral
	ascProbsAncestral = snpOutgroupAscMatrix(params->maxSampleSize, params->sampleSizeVector, \
		params->ascSize, 0, params->runMode);

	ascProbsDerived = snpOutgroupAscMatrix(params->maxSampleSize, params->sampleSizeVector, \
		params->ascSize, 1, params->runMode);
	
	//make vectors of denoms- derived, ancestral
	ascDenomsDerived = gsl_vector_alloc(params->maxSampleSize + 1);
	ascDenomsAncestral = gsl_vector_alloc(params->maxSampleSize + 1);
	gsl_vector_set_zero(ascDenomsDerived);
	gsl_vector_set_zero(ascDenomsAncestral);
	for(i = params->ascSize; i <= params->maxSampleSize; i++){
		gsl_vector_set(ascDenomsDerived,i,probAscertainmentGivenModelLookup(beta, i, probs, ascProbsDerived));
		gsl_vector_set(ascDenomsAncestral,i,probAscertainmentGivenModelLookup(beta, i, probs, ascProbsAncestral));
	}

//now go and tally likelihood
	lik = 0;
	num = 0;
	for(i = 0; i < params->snpNumber; i++){
		if(params->data[i].derivedFlag == 1){
			num = (gsl_matrix_get(probs,params->data[i].n,params->data[i].i) * gsl_matrix_get(ascProbsDerived,params->data[i].n,params->data[i].i)) \
				/ gsl_vector_get(ascDenomsDerived, params->data[i].n); 
			}
		else{
			num = (gsl_matrix_get(probs,params->data[i].n,params->data[i].i) * gsl_matrix_get(ascProbsAncestral,params->data[i].n,params->data[i].i)) \
				/ gsl_vector_get(ascDenomsAncestral, params->data[i].n);
		}	
		lik += log(num);
	}

	gsl_matrix_free(probs);
	gsl_matrix_free(ascProbsDerived);
	gsl_matrix_free(ascProbsAncestral);
	gsl_vector_free(ascDenomsDerived);
	gsl_vector_free(ascDenomsAncestral);
	return -lik;
}


/* wrapper for my_lik for use in mnbrak */
double  wlikWrapLookAsc(double beta, void * p){
	gsl_function f;
	struct my_lik_params * params = (struct my_lik_params *) p;
	
	f.function = &weightedLikLookAsc;
	f.params = params;
	return GSL_FN_EVAL(&f, beta);
}

/* wrapper for my_lik for use in mnbrak */
double  wlikWrapLookOutgroupAsc(double beta, void * p){
	gsl_function f;
	struct my_lik_params * params = (struct my_lik_params *) p;
	
	f.function = &weightedLikLookOutgroupAsc;
	f.params = params;
	return GSL_FN_EVAL(&f, beta);
}

/* get weighted ML estimate of beta using likelihood function (lookup) conditional on ascertainment
; BRENT method for gsl */
double weighted_ml_est_lookup_asc(double * lik_beta_hat, void * p){
	int status;
	int iter = 0;
	int max_iter = 100;
	const gsl_min_fminimizer_type *T;
	gsl_min_fminimizer *s;
	double m, a, b, ax, bx, cx, fa, fb, fc, dummy;
	gsl_function L;
	struct my_lik_params * params = (struct my_lik_params *) p;	
	
	
	// get minimum bracket 
	ax = -20.0;
	bx = -1.0;
	cx = 20;
	L.function = &weightedLikLookAsc;
	L.params = params;
	mnbrak2(&ax, &bx, &cx, &fa, &fb, &fc, &wlikWrapLookAsc, bx, params);
	//swap bounds if needed
	if (cx < ax){
	  dummy = cx;
	  cx = ax;
	  ax = dummy;
	}
	// do mimization 
		
	//initialize lik function
	L.function = &weightedLikLookAsc;
	L.params = params;
	  
	//min routine
	T = gsl_min_fminimizer_brent;
	s = gsl_min_fminimizer_alloc(T);
	gsl_min_fminimizer_set(s, &L, bx, ax, cx);
	do{
	  iter++;
	  status = gsl_min_fminimizer_iterate(s);
	  m = gsl_min_fminimizer_x_minimum (s);
	  a = gsl_min_fminimizer_x_lower (s);
	  b = gsl_min_fminimizer_x_upper (s);
	  status = gsl_min_test_interval (a, b, 0.001, 0.0);
	}
	while (status == GSL_CONTINUE && iter < max_iter);
	*lik_beta_hat = gsl_min_fminimizer_f_minimum(s);
	gsl_min_fminimizer_free(s);
	return m;
}

/* get weighted ML estimate of beta using likelihood function (lookup) conditional on outgroup ascertainment
; BRENT method for gsl */
double weighted_ml_est_lookup_outgroup_asc(double * lik_beta_hat, void * p){
	int status;
	int iter = 0;
	int max_iter = 100;
	const gsl_min_fminimizer_type *T;
	gsl_min_fminimizer *s;
	double m, a, b, ax, bx, cx, fa, fb, fc, dummy;
	gsl_function L;
	struct my_lik_params * params = (struct my_lik_params *) p;	
	
	
	// get minimum bracket 
	ax = -20.0;
	bx = -1.0;
	cx = 20;
	L.function = &weightedLikLookOutgroupAsc;
	L.params = params;
	mnbrak2(&ax, &bx, &cx, &fa, &fb, &fc, &wlikWrapLookOutgroupAsc, bx, params);
	//swap bounds if needed
	if (cx < ax){
	  dummy = cx;
	  cx = ax;
	  ax = dummy;
	}
	// do mimization 
		
	//initialize lik function
	L.function = &weightedLikLookOutgroupAsc;
	L.params = params;
	  
	//min routine
	T = gsl_min_fminimizer_brent;
	s = gsl_min_fminimizer_alloc(T);
	gsl_min_fminimizer_set(s, &L, bx, ax, cx);
	do{
	  iter++;
	  status = gsl_min_fminimizer_iterate(s);
	  m = gsl_min_fminimizer_x_minimum (s);
	  a = gsl_min_fminimizer_x_lower (s);
	  b = gsl_min_fminimizer_x_upper (s);
	  status = gsl_min_test_interval (a, b, 0.001, 0.0);
	}
	while (status == GSL_CONTINUE && iter < max_iter);
	*lik_beta_hat = gsl_min_fminimizer_f_minimum(s);
	gsl_min_fminimizer_free(s);
	return m;
}

/* likelihood function for SFS where each SNP has
 independent selection coeff and outgroup ascertainment */
double sfsLikBetaVectorOutgroupAsc(gsl_vector *betas, void * p){
	struct my_lik_params * params = (struct my_lik_params *) p;
	double pSNP, pAsc, pAscModel,lik, tmp;
	int i;
	gsl_matrix *ascProbsDerived, *ascProbsAncestral;
	gsl_vector *snpProbs;

	lik = 0;
//make ascProb matrix 
	ascProbsAncestral = snpOutgroupAscMatrix(params->maxSampleSize, params->sampleSizeVector, \
		params->ascSize, 0, params->runMode);

	ascProbsDerived = snpOutgroupAscMatrix(params->maxSampleSize, params->sampleSizeVector, \
		params->ascSize, 1, params->runMode);
//go through snps
	for(i = 0; i < params->snpNumber; i++){
	//create vector of snpProbs, find pSNP_i
		snpProbs = snpProbVectorNotLog(gsl_vector_get(betas, i), params->data[i].n);
		pSNP = gsl_vector_get(snpProbs, params->data[i].i);

	//lookup pAsc_i
		if(params->data[i].derivedFlag){
			pAsc = gsl_matrix_get(ascProbsDerived, params->data[i].n, params->data[i].i);
			//get pAscModel
			pAscModel = probAscertainmentGivenModelHemiLookup(gsl_vector_get(betas, i),params->data[i].n, \
					snpProbs, ascProbsDerived);
		}
		else{
			pAsc = gsl_matrix_get(ascProbsAncestral, params->data[i].n, params->data[i].i);
			//get pAscModel
			pAscModel = probAscertainmentGivenModelHemiLookup(gsl_vector_get(betas, i),params->data[i].n, \
					snpProbs, ascProbsAncestral);
		}
		tmp = pSNP * pAsc / pAscModel;
		lik += log(tmp);
		gsl_vector_free(snpProbs);
	}
	gsl_matrix_free(ascProbsAncestral);
	gsl_matrix_free(ascProbsDerived);
	return -lik;
}

double sfsLikBetaVectorOutgroupAscPre(gsl_vector *betas, gsl_matrix *ascProbsAncestral,gsl_matrix *ascProbsDerived, gsl_vector *snpProbs, void * p){
	struct my_lik_params * params = (struct my_lik_params *) p;
	double pSNP, pAsc, pAscModel,lik, tmp;
	int i;

	lik = 0.0;
//make ascProb matrix 
	ascProbsAncestral = snpOutgroupAscMatrixPre(params->maxSampleSize, params->sampleSizeVector, \
		params->ascSize, 0, params->runMode,ascProbsAncestral);

	ascProbsDerived = snpOutgroupAscMatrixPre(params->maxSampleSize, params->sampleSizeVector, \
		params->ascSize, 1, params->runMode,ascProbsDerived);
//go through snps
	for(i = 0; i < params->snpNumber; i++){
	//create vector of snpProbs, find pSNP_i
		snpProbs = snpProbVectorNotLogPre(gsl_vector_get(betas, i), params->data[i].n,snpProbs);
		pSNP = gsl_vector_get(snpProbs, params->data[i].i);

	//lookup pAsc_i
		if(params->data[i].derivedFlag){
			pAsc = gsl_matrix_get(ascProbsDerived, params->data[i].n, params->data[i].i);
			//get pAscModel
			pAscModel = probAscertainmentGivenModelHemiLookup(gsl_vector_get(betas, i),params->data[i].n, \
					snpProbs, ascProbsDerived);
		}
		else{
			pAsc = gsl_matrix_get(ascProbsAncestral, params->data[i].n, params->data[i].i);
			//get pAscModel
			pAscModel = probAscertainmentGivenModelHemiLookup(gsl_vector_get(betas, i),params->data[i].n, \
					snpProbs, ascProbsAncestral);
		}
		tmp = pSNP * pAsc / pAscModel;
		lik += log(tmp);
//		printf("%f\n",lik);
	}
	return(-lik);
}
//classical summary stats

//tajimas theta_pi
double theta_pi(struct snp *data, int sites){
	int i;
	double n,p, sum;
	
	sum = 0.0;
	for(i=0;i<sites;i++){
		n = (double) data[i].n;
		p = (double) data[i].i / n;
		sum += (n / (n-1.0)) * (2.0 * p * (1.0 - p));
	}
	
	return(sum);
}

//this simply counts size with 0 < i < sampleSize
int segSites(struct snp *data, int sites){
	int i, count = 0;
	//go thru sites
	for(i=0;i<sites;i++){
		//check if monomorphic ##note I'm doing zero and n
		if  (data[i].i != 0 && data[i].i != data[i].n) 
			count++;
	}
	return(count);
}	
/*
double theta_w(int segSites, int n){
	return((double)segSites / a1f(n));
}
	
	
//harmonic sum for Waterson's theta
double a1f(int n){
	int i;
	double sum = 0;
		
	for(i=1;i<n;i++){
		sum += 1.0 / (double) i;	
	}
	return(sum);
} 
/*
//tajd
double tajd(int nsam, int segsites, double sumk){
  double  a1, a2, b1, b2, c1, c2, e1, e2; 
  
  if( segsites == 0 ) return( 0.0) ;
  
  a1 = a1f(nsam);
  a2 = a2f(nsam);
  b1 = b1f(nsam);
  b2 = b2f(nsam);
  c1 = c1f(a1, b1);
  c2 = c2f(nsam, a1, a2, b2);
  e1 = e1f(a1, c1);
  e2 = e2f(a1, a2, c2);

  return( (sumk - (segsites/a1))/sqrt((e1*segsites) + ((e2*segsites)*(segsites
								      -1))) ) ;
}

double a2f(int nsam) {
  double a2;
  int i;
  a2 = 0.0;
  for (i=1; i<=nsam-1; i++) a2 += 1.0/(i*i);
  return (a2);
}


double b1f(int nsam){
  double b1;
  b1 = (nsam + 1.0)/(3.0*(nsam-1.0));
  return (b1);
}


double b2f(int nsam){
  double b2;
  b2 = (2*(nsam*nsam + nsam + 3.0))/(9*nsam*(nsam - 1));
  return (b2);
}


double e1f(double a1f, double c1){
  double e1;
  e1 = c1/a1f;
  return (e1);
}

double e2f(double a1f, double a2, double c2){ 
  double e2;
  e2 = c2/((a1f*a1f)+a2);
  return (e2);
}

double c1f(double a1f, double b1){
  double c1;
  c1 = b1 - (1/a1f);
  return (c1);
}

double c2f(int nsam, double a1f, double a2, double b2){
  double c2;
  c2 = b2 - ((nsam+2)/(a1f*nsam)) + (a2/(a1f * a1f));
  return (c2);
}


*/

/////////
/////
////  MPI stuff


//MPI version of above for parallel processing
double sfsLikBetaVectorMPI(gsl_vector *betas, void * p, int myrank, int mprocs, MPI_Comm comm){
	struct my_lik_params * params = (struct my_lik_params *) p;
	struct my_F_params FParams;
	gsl_function f;
	double lik, final;
	int snpNumber = params->snpNumber;
	int i, start, end, window;
	
	//initialize; where am i?
	lik = 0.;
	window = ceil((float)snpNumber / mprocs);
	start = myrank * window;
	end = start + window;
	FParams.w = params->w;
	f.function = &snpProb;
	for(i = start; i < end && i < snpNumber; i++){
		FParams.i = params->data[i].i;
		FParams.n = params->data[i].n;

		f.params = &FParams;
	//	printf("here lik:%f beta:%f i: %d n: %d\n",lik,gsl_vector_get(betas, i),FParams.i,FParams.n);
		
		lik += log(GSL_FN_EVAL(&f, gsl_vector_get(betas, i)));
		}
	MPI_Reduce(&lik,&final,1, MPI_DOUBLE,MPI_SUM,0,comm);
	return -final;
}

double sfsLikBetaVectorOutgroupAscPreMPI(gsl_vector *betas, gsl_matrix *ascProbsAncestral,gsl_matrix *ascProbsDerived, \
	gsl_vector *snpProbs, void * p,int myrank, int mprocs, MPI_Comm comm){
		
	struct my_lik_params * params = (struct my_lik_params *) p;
	double pSNP, pAsc, pAscModel,lik, tmp, final;
	int i,start, end, window,snpNumber;
	
	//initialization; window for proc to process
	lik = 0;
	snpNumber = params->snpNumber;
	window = ceil((float)snpNumber / mprocs);
	start = myrank * window;
	end = start + window;
	
	//currently every process calculates the matrices below-- could be faster
	//make ascProb matrix 
	ascProbsAncestral = snpOutgroupAscMatrixPre(params->maxSampleSize, params->sampleSizeVector, \
		params->ascSize, 0, params->runMode,ascProbsAncestral);

	ascProbsDerived = snpOutgroupAscMatrixPre(params->maxSampleSize, params->sampleSizeVector, \
		params->ascSize, 1, params->runMode,ascProbsDerived);
		
	//go through snps
	for(i = start; i < end && i < snpNumber; i++){
	//create vector of snpProbs, find pSNP_i
		snpProbs = snpProbVectorNotLogPre(gsl_vector_get(betas, i), params->data[i].n,snpProbs);
		pSNP = gsl_vector_get(snpProbs, params->data[i].i);

	//lookup pAsc_i
		if(params->data[i].derivedFlag){
			pAsc = gsl_matrix_get(ascProbsDerived, params->data[i].n, params->data[i].i);
			//get pAscModel
			pAscModel = probAscertainmentGivenModelHemiLookup(gsl_vector_get(betas, i),params->data[i].n, \
					snpProbs, ascProbsDerived);
		}
		else{
			pAsc = gsl_matrix_get(ascProbsAncestral, params->data[i].n, params->data[i].i);
			//get pAscModel
			pAscModel = probAscertainmentGivenModelHemiLookup(gsl_vector_get(betas, i),params->data[i].n, \
					snpProbs, ascProbsAncestral);
		}
		tmp = pSNP * pAsc / pAscModel;
		lik += log(tmp);
	}
	MPI_Reduce(&lik,&final,1, MPI_DOUBLE,MPI_SUM,0,comm);
	return -final;
}


