/* file write_force_step.c */
// c code for write_force_step subroutine

#include "stringlib.h"
#include "write_force_step.h"

void write_force_step(int nAtoms, int iter, double **atomForces, FILE *forOut) {

	int atom;

	fprintf(forOut, "Step %d\n", iter);
	for(atom=0; atom<nAtoms; atom++) {
		fprintf(forOut, "Ar%16.6E%16.6E%16.6E\n", atomForces[atom][0], atomForces[atom][1], atomForces[atom][2]);
	}
}
