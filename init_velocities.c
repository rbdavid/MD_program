/* file init_velocities.c */
// c code for init_velocities subroutine

#include <stdio.h>
#include <stdlib.h>
#include "init_velocities.h"

void init_velocities(int nAtoms, int ig, double Ar_m, double kBT, double **atomVelocities) {
	
	double v[3] = {0.0, 0.0, 0.0};
	double sumv[3] = {0.0, 0.0, 0.0};
	double sumv2;
	int i, j;
	double scale_factor;
	double dof;
	double convert;

	double sqrt(double x);

	srand((unsigned) ig);

	dof = nAtoms*3;

	convert = 10E-2; 		// conversion factor to turn m s^-1 to Ang ps^-1

	for(i=0; i<nAtoms; i++) {
		for(j=0; j<3; j++) {
			atomVelocities[i][j] = (rand()/((double) RAND_MAX) -0.5);
			sumv[j] += atomVelocities[i][j];
			sumv2 += atomVelocities[i][j]*atomVelocities[i][j];
		}
	}

	sumv2 = sumv2/dof;		// mean squared velocity for randomly assigned velocities (around -0.5 and 0.5)

	scale_factor = sqrt(3*kBT/(Ar_m*sumv2)); 

	for(i=0; i<nAtoms; i++) {
		for(j=0; j<3; j++) {
			sumv[j] = sumv[j]/nAtoms;
			atomVelocities[i][j] = (atomVelocities[i][j] - sumv[j])*scale_factor*convert;
		}
	}
}

