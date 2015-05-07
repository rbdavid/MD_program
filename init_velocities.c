/* file init_velocities.c */
// c code for init_velocities subroutine

#include <stdio.h>
#include <stdlib.h>
#include "stringlib.h"
#include "init_velocities.h"

void init_velocities(int nAtoms, int ig, double Ar_m, double kBT, double **atomVelocities) {
	
	double v[3] = {0.0, 0.0, 0.0};
	double sumv[3] = {0.0, 0.0, 0.0};
	double sumv2;
	int i, j;
	double scale_factor;
	double dof;
	double convert;
	double convert2;

	double sqrt(double x);

	srand((unsigned) ig);

	dof = nAtoms*3;

	sumv2 = 0.0;

	convert = 0.0001; 		// conversion factor to turn m^2 s^-2 to Ang^2 ps^-2
//	convert2 = 0.01;		// conversion factor to turn m^1 s^-1 to Ang^1 ps^-1
//	convert2 = 1E-2;

	for(i=0; i<nAtoms; i++) {
		for(j=0; j<3; j++) {
			atomVelocities[i][j] = (rand()/((double) RAND_MAX) -0.5);
			sumv[j] += atomVelocities[i][j];
			sumv2 += atomVelocities[i][j]*atomVelocities[i][j];
		}
	}
	
	printf("first particle, unscaled velocity: %f\n", atomVelocities[0][0]);

	sumv2 = sumv2/dof;		// mean squared velocity for randomly assigned velocities (around -0.5 and 0.5)
//	sumv2 = sumv2/nAtoms;		// mean squared velocity for randomly assigned velocities (around -0.5 and 0.5)
	printf("average velocity squared is %f\n", sumv2);
	
	scale_factor = sqrt(3*kBT*convert/(Ar_m*sumv2)); 
//	scale_factor = sqrt(3*kBT/(Ar_m*sumv2)); 

	printf("scale factor equals %f", scale_factor);

	for(j=0; j<3; j++) {
		sumv[j] = sumv[j]/nAtoms;
	}

	printf("average velocity is %f\n", sumv);

	for(i=0; i<nAtoms; i++) {
		for(j=0; j<3; j++) {
			atomVelocities[i][j] = (atomVelocities[i][j] - sumv[j])*scale_factor;
//			atomVelocities[i][j] = scale_factor*convert2*(atomVelocities[i][j] - sumv[j]);
		}
	}
	
	printf("first particle, scaled velocity: %f\n", atomVelocities[0][0]);
}

