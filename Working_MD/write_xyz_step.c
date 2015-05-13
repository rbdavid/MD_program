/* file write_xyz_step.c */
// c code for the write_xyz_step subroutine

#include "stringlib.h"
#include "write_xyz_step.h"

void write_xyz_step(int nAtoms, int iter, double box, double **coord, FILE *xyzOut) {

	int atom;
	
	fprintf(xyzOut, "%d \n", nAtoms);
	fprintf(xyzOut, "Step %d , box: %8.3f  %8.3f  %8.3f\n", iter, box, box, box);
	for (atom=0;atom<nAtoms;atom++) {
		fprintf(xyzOut, "Ar%12.6f%12.6f%12.6f\n",coord[atom][0],coord[atom][1],coord[atom][2]);

	}
}

