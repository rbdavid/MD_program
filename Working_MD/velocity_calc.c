/* file velocity_calc.c */
// c code for velocity_calc subroutines

#include <math.h>
#include "velocity_calc.h"

void velocity_calc(int nAtoms, int iter, int deltaWrite, double Ar_mmass, double delta_t, double R, double *int_temp, double *Tot_kinetic_en, double **atomVelocities, double **atomForces, double **old_atomForces) {
	
	int i, j;
	double component;
	double sumv2;
	double dof;
	double convert;
	
	convert = 0.0001*4184;

	dof = nAtoms*3.0;

	sumv2 = 0.0;

	*Tot_kinetic_en = 0.0;
	*int_temp = 0.0;

	for(j=0;j<3;j++) {
		for(i=0;i<nAtoms;i++) {
			atomVelocities[i][j] += 0.5*convert*delta_t*(atomForces[i][j] + old_atomForces[i][j])/Ar_mmass;

			if(iter%deltaWrite==0) {
				sumv2 += atomVelocities[i][j]*atomVelocities[i][j];
			}
		}
	}

	*Tot_kinetic_en = 0.5*Ar_mmass*sumv2/convert;
	sumv2 = sumv2/dof;
	*int_temp = (Ar_mmass*sumv2)/(R*convert);

}	
