
#include <stdio.h>                                                    
#include <math.h>                                                     
#include <time.h>                                                     
#include "stringlib.h"                                                
#include "read_cfg_file.h"
#include "init_positions.h"
#include "init_velocities.h"
#include "init_old_coord.h"
#include "force_energy_calc.h"
#include "write_xyz_step.h"
#include "write_vel_step.h"
#include "write_force_step.h"
#include "write_log_step.h"
#include "verlet.h"
#include "velocity_calc.h"
#include "positions_calc.h"

//  ---------------------------------
// Declare Subroutines                                                
//  ---------------------------------

//  ---------------------------------
// Declare Global Variables
//  ---------------------------------

static double R = 0.0019872;		// units: kcal mol^-1 K^-1

//  ---------------------------------
// Main Program                                                        
//  ---------------------------------

int main() {

//  ---------------------------------	
// variable declaration for reading config file
//  ---------------------------------
	
	char trajFileName[1024];  	// output trajectory file name
	char logFileName[1024];   	// log file name
	char velFileName[1024];
	char forFileName[1024];

	double temp;    		// temperature
	int nAtoms; 			// Number of atoms

	int nIter;      		// number of MD iterations
	int deltaWrite; 		// how often to write coordinates and log info 
	double cutoff;			// cutoff distance (Angstrom)
	double cutoff2;			// nonbonding interaction cutoff distance squared (Angstrom^2)
	int ig;				// random number for velocity initialization
	double delta_t;			// delta t value; MD time step; ps 

	double Ar_m;			// mass of Ar in kg
	double Ar_mmass;		// molar mass of Ar in kg mol^-1
	double Ar_eps;			// epsilon parameter for LJ energy calculation
	double Ar_sigma;		// sigma value...
	double Ar_sigma6;		// sigma^6 for quick LJ energy calculation...

	double RT;
	double delta_t2;

	int i, iter;			// generic indeces

//  ---------------------------------
// variable declaration for initializing particle positions, velocities
//  ---------------------------------

	double box; 			// cubic box dimension

	double **old_coord;

	double **coord;

	double **atomVelocities;	// variable declaration for the velocities of particles array

	double **atomForces;		// variable declaration for force acting on atom array

	double **old_atomForces;

	double Tot_potential_en;	// total potential energy of the system

	double int_temp;

	double Tot_kinetic_en;

	double Tot_en;

//  ---------------------------------
// Declaration of Timing Variables
//  ---------------------------------

	time_t startTime;   // initial clock time
	time_t stopTime;    // final clock time
	time_t routineStartTime; // start time for a routine
	time_t routineStopTime;  // stop time for a routine
	double timeSpent; // amount of time in seconds 
	double MDLoopTime; // 
	
//  ---------------------------------
// variable declaration for writing to log, xyz, force, velocity, etc files...
//  ---------------------------------

	FILE *xyzOut;
	FILE *logOut;
	FILE *velOut;
	FILE *forOut;

//  ---------------------------------
// initializing job timing	
//  ---------------------------------

	startTime = clock();

//  ---------------------------------
// read config data from standard in
//  ---------------------------------

	read_cfg_file(trajFileName, logFileName, velFileName, forFileName, &nAtoms, &nIter, &deltaWrite, &ig, &temp, &Ar_eps, &Ar_sigma, &cutoff, &delta_t, &Ar_m, &Ar_mmass); 

	Ar_sigma6 = Ar_sigma*Ar_sigma*Ar_sigma*Ar_sigma*Ar_sigma*Ar_sigma;		// Units: Angstrom^6

	cutoff2 = cutoff*cutoff;	// Units: Angstrom^2

	RT = R*temp;

	delta_t2 = delta_t*delta_t;	// delta_t squared; used often in the position calculation calc; ps^2

//  ---------------------------------
// array memory assignment ???
//  ---------------------------------

	old_coord = (double**) calloc(nAtoms,sizeof(double*));
	coord = (double**) calloc(nAtoms,sizeof(double*));
	old_atomForces = (double**) calloc(nAtoms,sizeof(double*));
	atomForces = (double**) calloc(nAtoms,sizeof(double*));         // allocate forces array memory; was calloc(nAtoms,sizeof(double*));
	atomVelocities = (double**) calloc(nAtoms,sizeof(double*));     // allocate velocity array memory

	for (i=0; i<nAtoms; i++) {
		old_coord[i] = (double*) calloc(3,sizeof(double));
		coord[i] = (double*) calloc(3,sizeof(double));
		old_atomForces[i] = (double*) calloc(3,sizeof(double));
		atomForces[i] = (double*) calloc(3,sizeof(double));
		atomVelocities[i] = (double*) calloc(3,sizeof(double));
	}

//  ---------------------------------
// open new files for traj, log, forces, velocities, etc...
//  ---------------------------------

	xyzOut = fopen(trajFileName, "w");
	logOut = fopen(logFileName, "w");
	velOut = fopen(velFileName, "w");
	forOut = fopen(forFileName, "w");

//  ---------------------------------	
// initialize particle positions
//  ---------------------------------   

	MDLoopTime=0;

	iter =0;
	Tot_kinetic_en =0.0;
	Tot_en = 0.0;
	
	init_positions(coord,nAtoms,&box);		// Units of coord: Angstrom

	init_velocities(nAtoms, ig, Ar_mmass, RT, R, &int_temp, &Tot_kinetic_en, atomVelocities);

//	init_old_coord(nAtoms, delta_t, atomVelocities, coord, old_coord);

	force_energy_calc(nAtoms, iter, deltaWrite, box, cutoff2, Ar_eps, Ar_sigma6, &Tot_potential_en, coord, atomForces, old_atomForces);

	Tot_en = Tot_kinetic_en + Tot_potential_en;

	write_log_step(iter, logOut, logFileName, trajFileName, velFileName, forFileName, nAtoms, temp, nIter, delta_t, deltaWrite, cutoff, Ar_mmass, Ar_eps, Ar_sigma, RT, box, Tot_en, Tot_potential_en, Tot_kinetic_en, int_temp);

	write_xyz_step(nAtoms, iter, box, coord, xyzOut);

	write_vel_step(nAtoms, iter, atomVelocities, velOut);

	write_force_step(nAtoms, iter, atomForces, forOut);

	fflush(xyzOut);
	fflush(logOut);
	fflush(velOut);
	fflush(forOut);

//  ---------------------------------   
//	MD LOOP
//  ---------------------------------   

	routineStartTime =clock();

	for(iter=1;iter<=nIter;iter++) {
		Tot_en = 0.0;
		
		positions_calc(nAtoms, Ar_mmass, delta_t, delta_t2, box, coord, atomVelocities, atomForces);

		routineStartTime =clock();
		
		force_energy_calc(nAtoms, iter, deltaWrite, box, cutoff2, Ar_eps, Ar_sigma6, &Tot_potential_en, coord, atomForces, old_atomForces);
		routineStopTime = clock();
		
		MDLoopTime += (double)(routineStopTime-routineStartTime)/CLOCKS_PER_SEC;

		velocity_calc(nAtoms, iter, deltaWrite, Ar_mmass, delta_t, R, &int_temp, &Tot_kinetic_en, atomVelocities, atomForces, old_atomForces);
//		verlet(nAtoms, delta_t2, Ar_mmass, box, coord, old_coord, atomForces);

		if(iter%deltaWrite==0) {

			Tot_en = Tot_kinetic_en + Tot_potential_en;
			
			write_log_step(iter, logOut, logFileName, trajFileName, velFileName, forFileName, nAtoms, temp, nIter, delta_t, deltaWrite, cutoff, Ar_mmass, Ar_eps, Ar_sigma, RT, box, Tot_en, Tot_potential_en, Tot_kinetic_en, int_temp);

			write_xyz_step(nAtoms, iter, box, coord, xyzOut);
			
			write_vel_step(nAtoms, iter, atomVelocities, velOut);
			
			write_force_step(nAtoms, iter, atomForces, forOut);

			fflush(xyzOut);
			fflush(logOut);
			fflush(velOut);
			fflush(forOut);

		}	

	}

	fprintf(logOut, "Total time to compute forces (seconds): %f\n", MDLoopTime);

	stopTime = clock();
	timeSpent = (double)(stopTime-startTime)/CLOCKS_PER_SEC;
	fprintf(logOut, "Total job time (seconds): %f\n", timeSpent);
	fprintf(logOut, "Nanoseconds per day: %f \n", nIter*delta_t*86.4/timeSpent);
	
	fclose(xyzOut);
	fclose(logOut);
	fclose(velOut);
	fclose(forOut);

}

