/* file velocity_calc.c */
// c code for velocity_calc subroutines

#include <math.h>
#include "velocity_calc.h"

void velocity_calc(int nAtoms, int iter, int deltaWrite, double Ar_m, double Ar_mmass, double delta_t, double kB, double *int_temp, double *Tot_kinetic_en, double **atomVelocities, double **atomForces, double **old_atomForces) {
	
	int i, j;
	double component;
	double sumv2;
	double dof;
	double convert;
	double convert2;
	
	convert = 4184*0.0001;		// conversion factor for kcal m^2 s^-2 to J Ang^2 ps^2

	convert2 = 0.0001;

	dof = nAtoms*3.0;

	sumv2 = 0.0;

	for(i=0; i<nAtoms; i++) {
		for(j=0; j<3; j++) {
			component = atomVelocities[i][j] + 0.5*convert*delta_t*(atomForces[i][j] + old_atomForces[i][j])/Ar_mmass;
			atomVelocities[i][j] = component;
			
			if(iter%deltaWrite==0) {
				sumv2 += atomVelocities[i][j]*atomVelocities[i][j];
			}
		}
	}
	
	*Tot_kinetic_en = 0.5*Ar_mmass*sumv2/convert;
	sumv2 = sumv2/dof;
	*int_temp = (Ar_m*sumv2)/(kB*convert2);

}	
