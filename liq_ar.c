
#include <stdio.h>                                                    
#include <math.h>                                                     
#include <time.h>                                                     
#include "stringlib.h"                                                
#include "read_cfg_file.h"
#include "init_positions.h"
#include "force_energy_calc.h"
#include "write_xyz_step.h"

//  ---------------------------------
// Declare Subroutines                                                
//  ---------------------------------

//  ---------------------------------
// Declare Global Variables
//  ---------------------------------

static double kB = 1.9872041E-3; // actually R in units of kcal/mol/K

//  ---------------------------------
//Main Program                                                        
//  ---------------------------------

int main() {

//  ---------------------------------	
// variable declaration for reading config file
//  ---------------------------------
	
	char trajFileName[1024];  	// output trajectory file name
	char logFileName[1024];   	// log file name
	
	double temp;    		// temperature
	int nAtoms; 			// Number of atoms

	int nIter;      		// number of MC iterations
	int deltaWrite; 		// how often to write coordinates and log info in MC
	double cutoff;			// cutoff distance (Angstrom)
	double cutoff2;			// nonbonding interaction cutoff distance squared (Angstrom^2)
	int reevalvel;			// number of steps between reevaluation of velocities...
	int ig;				// random number for velocity initialization

	double Ar_eps;			// epsilon parameter for LJ energy calculation
	double Ar_sigma;		// sigma value...
	double Ar_sigma6;		// sigma^6 for quick LJ energy calculation...

	double kBT;

	int i, iter;				// generic indeces

//  ---------------------------------
// variable declaration for initializing particle positions, velocities
//  ---------------------------------

	double box; 			// cubic box dimension

	double **coord;			// variable declaration for the coordinates of particles array

//	double **atomVelocities;	// variable declaration for the velocities of particles array

	double **atomForces;		// variable declaration for force acting on atom array

	double Tot_potential_en;	// total potential energy of the system

//  ---------------------------------
// variable declaration for writing to log, xyz, force, velocity, etc files...
//  ---------------------------------

	FILE *xyzOut;
	FILE *logOut;

//  ---------------------------------
// read config data from standard in
//  ---------------------------------

	read_cfg_file(trajFileName, logFileName, &nAtoms, &nIter, &deltaWrite, &reevalvel, &ig, &temp, &Ar_eps, &Ar_sigma, &cutoff); 

	Ar_sigma6 = Ar_sigma*Ar_sigma*Ar_sigma*Ar_sigma*Ar_sigma*Ar_sigma;

	cutoff2 = cutoff*cutoff;

//  ---------------------------------
// open new files for traj, log, forces, velocities, etc...
//  ---------------------------------

	xyzOut = fopen(trajFileName, "w");
	logOut = fopen(logFileName, "w");
//      velocitiesOut = fopen(velocityFileName, "w");
//      forcesOut = fopen(forcesFileName, "w");

	fprintf(logOut, "Molecular Dynamics Simulation of Argon; developed by RBD\n");
	fprintf(logOut, "Name of trajectory file: %s\n", trajFileName);
	fprintf(logOut, "Name of log file: %s\n", logFileName);
	fprintf(logOut, "Number of atoms: %d\n", nAtoms);
	fprintf(logOut, "Temperature: %f\n", temp);
	fprintf(logOut, "Number of steps: %d\n", nIter);
	fprintf(logOut, "Write to trajectory and log files every %d steps\n", deltaWrite);
	fprintf(logOut, "Nonbonding distance cutoff: %f Angstroms\n", cutoff);
	fprintf(logOut, "Re-evaluate velocities every %d steps\n", reevalvel);
	fprintf(logOut, "Argon Parameters: \n");
	fprintf(logOut, "Epsilon (kcal/mol): %f\n", Ar_eps);
	fprintf(logOut, "Sigma (Angstrom): %f\n", Ar_sigma);
	
//  ---------------------------------
// array memory assignment ???
//  ---------------------------------

        coord = (double**) malloc(nAtoms*sizeof(double*));              // allocate coordinate array memory
	atomForces = (double**) calloc(nAtoms,sizeof(double*));         // allocate forces array memory; was calloc(nAtoms,sizeof(double*));
//      atomVelocities = (double**) malloc(nAtoms,sizeof(double*));     // allocate velocity array memory

	for (i=0; i<nAtoms; i++) {
		coord[i] = (double*) malloc(3*sizeof(double));
		atomForces[i] = (double*) malloc(3*sizeof(double));
//		atomVelocities[i] = (double*) malloc(3*sizeof(double));
	}

//  ---------------------------------
//  open new files for traj, log, forces, velocities, etc...
//  ---------------------------------	

	kBT = kB*temp;			// solve for kB*T value
	fprintf(logOut, "kB*T = %f kcal/mol \n", kBT);	// Print kB*T value to log file;

//  ---------------------------------	
// initialize particle positions
//  ---------------------------------   
	
	init_positions(coord,nAtoms,&box);

	fprintf(logOut, "Box dimension: %f Angstroms \n", box);

	// need to initialize the velocities for each particle... 

	Tot_potential_en = force_energy_calc(nAtoms, box, cutoff2, Ar_eps, Ar_sigma6, coord, atomForces);

	fprintf(logOut, "Initial Total potential energy = %f\n", Tot_potential_en);


	fflush(xyzOut);
	fflush(logOut);

	fclose(xyzOut);
	fclose(logOut);

}


