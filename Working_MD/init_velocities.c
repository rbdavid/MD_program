/* file init_velocities.c */
// c code for init_velocities subroutine

#include <stdio.h>
#include <stdlib.h>
#include "stringlib.h"
#include "init_velocities.h"

void init_velocities(int nAtoms, int ig, double Ar_mmass, double RT, double R, double *int_temp, double *Tot_kinetic_en, double **atomVelocities) {
	
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

	convert2 = 4184*convert;

	*Tot_kinetic_en = 0.0;
	*int_temp = 0.0;

	for(j=0; j<3; j++) {
		sumv[j] = 0.0;
		for(i=0; i<nAtoms; i++) {
			atomVelocities[i][j] = (rand()/((double) RAND_MAX) -0.5);
			sumv[j] += atomVelocities[i][j];
			sumv2 += atomVelocities[i][j]*atomVelocities[i][j];
		}
		sumv[j] = sumv[j]/nAtoms;
	}

	printf("average velocity (Ang ps^-1): %12.6f%12.6f%12.6f\n", sumv[0], sumv[1], sumv[2]);

	printf("Total velocity squared (Ang^2 ps^-2): %12.6f\n", sumv2);

	sumv2 = sumv2/dof;		// mean squared velocity for randomly assigned velocities (around -0.5 and 0.5)
	
	printf("average velocity squared (Ang^2 ps^-2): %12.6f\n", sumv2);

	scale_factor = sqrt(RT*convert2/(Ar_mmass*sumv2)); 

	printf("scale factor = %12.6f\n", scale_factor);

	sumv2 = 0.0;

	for(i=0; i<nAtoms; i++) {
		for(j=0; j<3; j++) {
			sumv[j] = 0.0;
			atomVelocities[i][j] = (atomVelocities[i][j] - sumv[j])*scale_factor;
			sumv[j] += atomVelocities[i][j];
			sumv[j] = sumv[j]/nAtoms;
			sumv2 += atomVelocities[i][j]*atomVelocities[i][j];
		}
	}

	printf("After velocity scaling:\n");
	printf("average velocity (Ang ps^-1): %12.6f%12.6f%12.6f\n", sumv[0], sumv[1], sumv[2]);
	printf("Total velocity squared (Ang^2 ps^-2): %12.6f\n", sumv2);
	
	*Tot_kinetic_en = 0.5*Ar_mmass*sumv2/convert2;
	sumv2 = sumv2/dof;
	*int_temp = (Ar_mmass*sumv2)/(R*convert2);

	printf("Total Kinetic Energy (kcal mol^-1): %12.6f\n", *Tot_kinetic_en);
	printf("average velocity squared (Ang^2 ps^-2): %12.6f\n", sumv2);
	printf("Instantaneous Temperature (K): %12.6f\n", *int_temp);

}

