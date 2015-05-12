/* file positions_calc.c */
// c code for the positions_calc subroutine

#include <math.h>
#include "positions_calc.h"

void positions_calc(int nAtoms, double Ar_mmass, double delta_t, double delta_t2, double box, double **coord, double **atomVelocities, double **atomForces) {
	
	int i, j;
	double component;
	double convert;

	convert = 4184*0.0001;

	for(i=0; i<nAtoms; i++) {
		for(j=0; j<3; j++) {
			
			component = 0.0;
			// units of component: Angstrom
			component = coord[i][j] + delta_t*atomVelocities[i][j] + 0.5*delta_t2*convert*atomForces[i][j]/Ar_mmass;
			// wrapping the particles within the box...
			if(component<0) {
				component += box;
			} else if(component > box) {
				component -= box;
			}
			coord[i][j] = component;
		}
	}
}

