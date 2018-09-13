/* popGenTools.h
/
/
/ Andrew Kern
/
/
*/
#ifndef NGSPGTOOLS_INC
#define NGSPGTOOLS_INC


#include "stringWrap.h"
#include "sequenceMatrix.h"
#include "vector.h"
#include <stdlib.h>
#include <stdio.h>


double prob_i(int i, int n);
double prob_km(int k, int m, int n);

double achazSFS(gsl_matrix *sfsCount, gsl_matrix *weightMat);
double achazSFSCorrected(gsl_matrix *sfsCount, gsl_matrix *weightMat,int nPool);
void fillPiWeights(gsl_matrix *weightMat, int minFreq);
void fillHWeights(gsl_matrix *weightMat, int minFreq);
void fillWWeights(gsl_matrix *weightMat, int minFreq);


#endif
