/* file write_vel_step.c */
// c code for write_vel_step subroutine

#include "stringlib.h"
#include "write_vel_step.h"

void write_vel_step(int nAtoms, int iter, double **atomVelocities, FILE *velOut) {
	int atom;

	fprintf(velOut, "Step %d\n", iter);
	for(atom=0; atom<nAtoms; atom++) {
		fprintf(velOut, "Ar%16.6f%16.6f%16.6f\n", atomVelocities[atom][0], atomVelocities[atom][1], atomVelocities[atom][2]);
	}
}

