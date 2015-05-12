/* file verlet.c */
// c code for the verlet subroutine

#include <math.h>
#include "verlet.h"

void verlet(int nAtoms, double delta_t2, double Ar_mmass, double box, double **coord, double **old_coord, double **atomForces) {

	int atom1;
	int j;

	double component;
	double convert;

	convert = 4184*0.0001;

	for(atom1=0;atom1<nAtoms;atom1++) {
		for(j=0;j<3;j++) {
			component = 0.0;
			component = 2*coord[atom1][j] - old_coord[atom1][j] + 0.5*convert*delta_t2*atomForces[atom1][j]/Ar_mmass;
			
			if(component<0) {
				component += box;
			} else if(component > box) {
				component -= box;
			}
			
			old_coord[atom1][j] = coord[atom1][j];
			coord[atom1][j] = component;
		}
	}
}

