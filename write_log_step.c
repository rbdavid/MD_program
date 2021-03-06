/* file write_log_step.c */
// c code for write_log_step subroutine

#include "stringlib.h"
#include "write_log_step.h"

void write_log_step(int iter, FILE *logOut, char *logFileName, char *trajFileName, char *velFileName, char *forFileName, int nAtoms, double temp, int nIter, double delta_t, int deltaWrite, double cutoff, double Ar_m, double Ar_mmass, double Ar_eps, double Ar_sigma, double kBT, double box, double Tot_en, double Tot_potential_en, double Tot_kinetic_en, double int_temp) {

	if(iter==0) {
		fprintf(logOut, "Molecular Dynamics Simulation of Argon; developed by RBD\n");
		fprintf(logOut, "Name of log file: %s\n", logFileName);
		fprintf(logOut, "Name of trajectory file: %s\n", trajFileName);
		fprintf(logOut, "Name of velocity file: %s\n", velFileName);
		fprintf(logOut, "Name of force file: %s\n", forFileName);
		fprintf(logOut, "Number of atoms: %d\n", nAtoms);
		fprintf(logOut, "Temperature: %f\n", temp);
		fprintf(logOut, "Number of steps: %d\n", nIter);
		fprintf(logOut, "Time step: %f pecoseconds\n", delta_t);
		fprintf(logOut, "Write to trajectory and log files every %d steps\n", deltaWrite);
		fprintf(logOut, "Nonbonding distance cutoff: %f Angstroms\n", cutoff);
//		fprintf(logOut, "Re-evaluate velocities every %d steps\n\n", reevalvel);
		fprintf(logOut, "\nArgon Parameters: \n");
		fprintf(logOut, "Argon Mass: %E kg \n", Ar_m);
		fprintf(logOut, "Argon Molar Mass: %f kg mol^-1 \n", Ar_mmass);
		fprintf(logOut, "Epsilon (kcal/mol): %f\n", Ar_eps);
		fprintf(logOut, "Sigma (Angstrom): %f\n\n", Ar_sigma);
		fprintf(logOut, "kB*T = %E  \n", kBT);  // Print kB*T value to log file; units!!!
		fprintf(logOut, "Box dimension: %f Angstroms \n \n", box);

		fprintf(logOut, "Step %d\n", iter);
//		fprintf(logOut, "Total Energy (kcal mol^-1) = %E \n", Tot_en);
		fprintf(logOut, "Total Potential Energy (kcal mol^-1) = %E \n", Tot_potential_en);
//		fprintf(logOut, "Total Kinetic Energy (kcal mol^-1) = %E \n", Tot_kinetic_en);
//		fprintf(logOut, "Instantaneous Temperature (K) = %E \n", int_temp);
	} else {
		fprintf(logOut, "Step %d\n", iter);
		fprintf(logOut, "Total Energy (kcal mol^-1) = %E \n", Tot_en);
		fprintf(logOut, "Total Potential Energy (kcal mol^-1) = %E \n", Tot_potential_en);
		fprintf(logOut, "Total Kinetic Energy (kcal mol^-1) = %E \n", Tot_kinetic_en);
		fprintf(logOut, "Instantaneous Temperature (K) = %E \n", int_temp);
	}
	
//	if(iter==nIter) {
	//Timing information
//	}
}
