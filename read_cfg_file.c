/* file read_cfg_file.c */
// c code for the read_cfg_file subroutine

#include "read_cfg_file.h"
#include "stringlib.h"

void read_cfg_file(char *trajFileName, char *logFileName, double *temp, int *nAtoms, int *nIter, int *deltaWrite, int *reevalvel, int *ig, double *Ar_eps, double *Ar_sigma6, double *cutoff2, FILE *xyzOut, FILE *logOut) {

	char buffer[1024];
	char tempBuffer[1024];
	char check[15];
	char *firstWord;
	double Ar_sigma;
	double cutoff;

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
			Ar_sigma = atof(string_secondword(buffer));
			*Ar_sigma6 = Ar_sigma*Ar_sigma*Ar_sigma*Ar_sigma*Ar_sigma*Ar_sigma;
		} else if (strncmp(firstWord,"cutoff",6)==0) {
			cutoff = atof(string_secondword(buffer));
			*cutoff2 = cutoff*cutoff;
		}
	}

//  ---------------------------------
//  open new files for traj, log, forces, velocities, etc...
//  ---------------------------------   

	xyzOut = fopen(trajFileName, "w");
	logOut = fopen(logFileName, "w");
//	velocitiesOut = fopen(velocityFileName, "w");
//	forcesOut = fopen(forcesFileName, "w");

//  ---------------------------------
//      Print log file info
//  ---------------------------------

	fprintf(logOut, "Molecular Dynamics Simulation of Argon; developed by RBD\n");
	fprintf(logOut, "Name of trajectory file: %s\n",trajFileName);
	fprintf(logOut, "Name of log file: %s\n",logFileName);
	fprintf(logOut, "Number of atoms: %d\n",*nAtoms);
	fprintf(logOut, "Temperature: %f\n",*temp);
	fprintf(logOut, "Number of steps: %d\n",*nIter);
	fprintf(logOut, "Write to trajectory and log files every %d\n",*deltaWrite);
	fprintf(logOut, "Nonbonding distance cutoff: %f\n", cutoff);
	fprintf(logOut, "Re-evaluate velocities every %d\n", *reevalvel);
	fprintf(logOut, "Argon Parameters: \n");
	fprintf(logOut, "Epsilon (kcal/mol): %f\n", *Ar_eps);
	fprintf(logOut, "Sigma (Angstrom): %f\n", Ar_sigma);
	fprintf(logOut, "Sigma^6 (Angstrom^6): %f\n", *Ar_sigma6);

}

