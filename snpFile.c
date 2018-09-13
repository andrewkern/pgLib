/* snpFile utilities
/
/
/ Andrew Kern
*/

#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "assert.h"
#include "snpFile.h"


struct snpFile *snpFile_new(){
	struct snpFile *newSnp;
	
	newSnp = malloc(sizeof(snpFile));
	if(newSnp == NULL){
		fprintf(stderr,"didn't make snpFile malloc\n");
		exit(1);
	}
	newSnp->next = NULL;

	return(newSnp);
}

void snpFile_free(struct snpFile *aSnpFile){
	free(aSnpFile->data);
	aSnpFile->next = NULL;
	free(aSnpFile);
}

void addSnpFile(struct snpFile *parent,struct snpFile *new){
	parent->next = new;
}

//counts number of separate snp samples in snpFile list
int snpFileCount(struct snpFile *aSnpFile){
	int i = 0;
	struct snpFile *curPtr;
	
	curPtr = aSnpFile;
	while(curPtr != NULL){
		i++;
		curPtr = curPtr->next;
	}
	return(i);
}

//snpFileImport-- reads a file and stores info into pre-alloc'd data for 5 column
//returns snpNumber
int snpFileImport5(char *fileName, struct snp *data){
	FILE *infile;
	int pos, i, j, n, flag, ascSize;

	/* open file, errors? */
	infile = fopen(fileName, "r");
	if (infile == NULL){
		fprintf(stderr,"Error opening infile! ARRRRR!!!!\n");
		exit(1);
	}
	/* go through infile and get snp info; expect pos, i, n */
	j = 0;
	while (fscanf(infile, "%d %d %d %d %d", &pos, &i, &n, &flag, &ascSize) != EOF){
		data[j].pos = pos;
		data[j].i = i;
		data[j].n = n;
		data[j].derivedFlag = flag;
		data[j].ascSize = ascSize;
		j += 1;
	}
	fclose(infile);
	return(j);
}

//snpFileImport-- reads a file and stores info into pre-alloc'd data for 3 column
//returns snpNumber
int snpFileImport3(char *fileName, struct snp *data){
	FILE *infile;
	int pos, i, j, n;

	/* open file, errors? */
	infile = fopen(fileName, "r");
	if (infile == NULL){
		fprintf(stderr,"Error opening infile! ARRRRR!!!!\n");
		exit(1);
	}
	/* go through infile and get snp info; expect pos, i, n */
	j = 0;
	while (fscanf(infile, "%d %d %d", &pos, &i, &n) != EOF){
		data[j].pos = pos;
		data[j].i = i;
		data[j].n = n;
		j += 1;
	}
	fclose(infile);
	return(j);
}

//snpFileImport-- reads a file and stores info into pre-alloc'd data for 1 column
//returns snpNumber
int snpFileImport1(char *fileName, struct snp *data){
	FILE *infile;
	int pos,j;

	/* open file, errors? */
	infile = fopen(fileName, "r");
	if (infile == NULL){
		fprintf(stderr,"Error opening infile! ARRRRR!!!!\n");
		exit(1);
	}
	/* go through infile and get snp info; expect pos, i, n */
	j = 0;
	while (fscanf(infile, "%d", &pos) != EOF){
		data[j].pos = pos;
		j += 1;
	}
	fclose(infile);
	return(j);
}


struct snpFile *snpFileMultiple(char *fileName){
	FILE *infile;
	int pos, i, j, n;
	struct snpFile *rootFile, *lastFile, *newFile=NULL;
	char line[1001];

	/* open file, errors? */
	infile = fopen(fileName, "r");
	if (infile == NULL){
		fprintf(stderr,"Error opening infile! ARRRRR!!!!\n");
		exit(1);
	}
	//setup snpFile
	rootFile = NULL;
	lastFile = NULL;
	j= 0;
	/* go through infile and get snp info; expect pos, i, n */
	while (fgets( line, 1000, infile)){
		if (line[0] == '/'){
			if(rootFile == NULL){
				rootFile = snpFile_new();
				newFile = rootFile;
			}
			else{
			if(lastFile != NULL){
				addSnpFile(lastFile, newFile);
			}
			lastFile = newFile;
			newFile = snpFile_new();
			}
			j = 0;
		}
		else{
			sscanf(line, "%d %d %d", &pos, &i, &n);
			newFile->data[j].pos = pos;
			newFile->data[j].i = i;
			newFile->data[j].n = n;
			j += 1;
			newFile->snpNumber = j;
			assert(j<MAXSNPS);
		}
	}
	addSnpFile(lastFile, newFile);
	fclose(infile);
	return(rootFile);	
	
}


//returns averageSampleSize of snpFile
double averageSampleSizeSNPs(struct snp *data, int siteNumber){
	int i, sum = 0;
	
	for(i=0; i < siteNumber; i++){
		sum += data[i].n;
	}
	return((double) sum / siteNumber);
}
//returns max sampleSize of snpFile
int maxSampleSizeSNPs(struct snp *data, int siteNumber){
	int i, max = 0;
	
	for(i=0; i < siteNumber; i++){
		if(data[i].n > max){
			max = data[i].n;
		}
	}
	return(max);
}



