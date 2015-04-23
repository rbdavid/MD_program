
#include <stdio.h>                                                    
#include <math.h>                                                     
#include <time.h>                                                     
#include "stringlib.h"                                                
                                                                      
// Declare Subroutines                                                
void read_cfg_file(char *, char *, double *, int *, int *, int *, int *, int *, double *, double *, double *);

//Main Program                                                        

int main() {
	
	// variable declaration for reading config file
	
	char trajFileName[1024];  	// output trajectory file name
	char logFileName[1024];   	// log file name
	
	double temp;    		// temperature
	int nAtoms; 			// Number of atoms
//	double box;      		// cubic box size

	int nIter;      		// number of MC iterations
	int deltaWrite; 		// how often to write coordinates and log info in MC
	double cutoff2;			// nonbonding interaction cutoff distance squared (Angstrom^2)
	int reevalvel;			// number of steps between reevaluation of velocities...
	int ig;				// random number for velocity initialization

	double Ar_eps;			// epsilon parameter for LJ energy calculation
//	double Ar_sigma;		// sigma parameter for LJ energy calculation
	double Ar_sigma6;		// sigma^6 for quick LJ energy calculation...


	// read config data from standard in
	read_cfg_file(trajFileName, logFileName, &temp, &nAtoms, &nIter, &deltaWrite, &reevalvel, &ig, &Ar_eps, &Ar_sigma6, &cutoff2); 


}

// Subroutines

void read_cfg_file(char *trajFileName, char *logFileName, double *temp, int *nAtoms, int *nIter, int *deltaWrite, int *reevalvel, int *ig, double *Ar_eps, double *Ar_sigma6, double *cutoff2) {

	char buffer[1024];
	char tempBuffer[1024];
	char check[15];
	char *firstWord;
	double Ar_sigma;
	double cutoff;

	while (fgets(buffer,1024,stdin) != NULL) {
		strncpy(tempBuffer,buffer,1024);
		firstWord=string_firstword(tempBuffer);
//              printf ("First word = %s\n",firstWord);
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
			Ar_sigma = atof(string_secondword(buffer));
			*Ar_sigma6 = Ar_sigma*Ar_sigma*Ar_sigma*Ar_sigma*Ar_sigma*Ar_sigma;
		} else if (strncmp(firstWord,"cutoff",6)==0) {
			cutoff = atof(string_secondword(buffer));
			*cutoff2 = cutoff*cutoff;
		}
	}

// 	Print log file info
	printf("Molecular Dynamics Simulation of Argon; developed by RBD\n");
	printf("Name of trajectory file: %s\n",trajFileName);
	printf("Name of log file: %s\n",logFileName);
	printf("Number of atoms: %d\n",*nAtoms);
	printf("Temperature: %f\n",*temp);
	printf("Number of steps: %d\n",*nIter);
	printf("Write to trajectory and log files every %d\n",*deltaWrite);
//	printf("box dimension: %f\n", *box);
	printf("Nonbonding distance cutoff: %f\n", cutoff);
	printf("Re-evaluate velocities every %d\n", *reevalvel);
	printf("Argon Parameters: \n");
	printf("Epsilon (kcal/mol): %f\n", *Ar_eps);
	printf("Sigma (Angstrom): %f\n", Ar_sigma);
	printf("Sigma^6 (Angstrom^6): %f\n", *Ar_sigma6);

}


