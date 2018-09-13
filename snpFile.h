/* snpFile.h -- for dealing with the snpFile, file format 
/
/
/
/ Also contains the definition of a snp */


#ifndef SNP_INC
#define SNP_INC


#define MAXSNPS 2000000
#define MAXALLELES 50

 /*SNP struct definition- meant to hold position, derived freq, sample size */
struct snp{
  	int pos, i, n, ascSize, ascFreq, derivedFlag;
	char *array;
	double freq;
}snp;

struct snpFile{
	int snpNumber;
	struct snp data[MAXSNPS];
	struct snpFile *next;
}snpFile;

struct snpFile *snpFile_new();
void snpFile_free(struct snpFile *aSnpFile);
void addSnpFile(struct snpFile *parent, struct snpFile *new);
int snpFileCount(struct snpFile *aSnpFile);

struct snpFile *snpFileMultiple(char *fileName);
int snpFileImport3(char *fileName, struct snp *data);
int snpFileImport5(char *fileName, struct snp *data);
int snpFileImport1(char *fileName, struct snp *data);

double averageSampleSizeSNPs(struct snp *data, int siteNumber);
int maxSampleSizeSNPs(struct snp *data, int siteNumber);

#endif

