/* file positions_calc.c */
// c code for the positions_calc subroutine

#include <math.h>
#include "positions_calc.h"

void positions_calc(int nAtoms, double Ar_m_na2d, double delta_t, double delta_t2, double box, double **coord, double **atomVelocities, double **atomForces) {
	
	int i, j;
	double component;
	double convert;

	convert = 10.E-4;

	for(i=0; i<nAtoms; i++) {
		for(j=0; j<3; j++) {

			// units of component: Angstrom
			component = coord[i][j] + delta_t*atomVelocities[i][j] + delta_t2*Ar_m_na2d*convert*atomForces[i][j];
			// wrapping the particles within the box...
			if(component<0) {
				component += box;
			} else if(component > box) {
				component -= box;
			}
			coord[i][j] = component;
			component = 0.0;
		}
	}
}

