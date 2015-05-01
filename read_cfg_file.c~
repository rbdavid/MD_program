/* file read_cfg_file.c */
// c code for the read_cfg_file subroutine

#include "read_cfg_file.h"
#include "stringlib.h"

void read_cfg_file(char *trajFileName, char *logFileName, int *nAtoms, int *nIter, int *deltaWrite, int *reevalvel, int *ig, double *temp, double *Ar_eps, double *Ar_sigma, double *cutoff) {	

	char buffer[1024];
	char tempBuffer[1024];
	char check[15];
	char *firstWord;

	while (fgets(buffer,1024,stdin) != NULL) {
		strncpy(tempBuffer,buffer,1024);
		firstWord=string_firstword(tempBuffer);
		if (strncmp(firstWord,"trajFile",8)==0) {
			strcpy(trajFileName,string_secondword(buffer));
		} else if (strncmp(firstWord,"logFile",7)==0) {
			strcpy(logFileName,string_secondword(buffer));
		} else if (strncmp(firstWord,"temperature",11)==0) {
			*temp = atof(string_secondword(buffer));
		} else if (strncmp(firstWord,"nAtoms",6)==0) {
			*nAtoms = atoi(string_secondword(buffer));
		} else if (strncmp(firstWord,"nIter",5)==0) {
			*nIter = atoi(string_secondword(buffer));
		} else if (strncmp(firstWord,"deltaWrite",10)==0) {
			*deltaWrite = atoi(string_secondword(buffer));
		} else if (strncmp(firstWord,"reevalvel",9)==0) {
			*reevalvel = atoi(string_secondword(buffer));
		} else if (strncmp(firstWord,"ig",2)==0) {
			*ig = atoi(string_secondword(buffer));
		} else if (strncmp(firstWord,"Ar_eps",6)==0) {
			*Ar_eps = atof(string_secondword(buffer));
		} else if (strncmp(firstWord,"Ar_sigma",8)==0) {
			*Ar_sigma = atof(string_secondword(buffer));
		} else if (strncmp(firstWord,"cutoff",6)==0) {
			*cutoff = atof(string_secondword(buffer));
		}
	}
}

