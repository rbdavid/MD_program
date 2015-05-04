/* file init_old_coord.c */
// c code for init_old_coord subroutine

#include "init_old_coord.h"

void init_old_coord(double nAtoms, double delta_t, double **coord, double **atomVelocities, double **old_coord) {

	int i, j;

	for(i=0; i<nAtoms; i++) {
		for(j=0; j<3; j++) {
		old_coord[i][j] = coord[i][j] - atomVelocities[i][j]*delta_t;	
		}
	}
}	

