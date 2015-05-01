/* file init_positions.c */
// c code for init_positions subroutine

#include <stdio.h>
#include "init_positions.h"

void init_positions(double **coord, int nAtoms, double *box) {

	double cbrt(double x);          // cube root function
	int iBoxD;                      // integer box dimension
	double fBoxD;                   // float box dimension
	
	int x, y, z;
	double xPos, yPos, zPos;
	int atomCount;

	// determine how many bins to divide the box into
	iBoxD = (int) cbrt((double) nAtoms);
	if (iBoxD*iBoxD*iBoxD < nAtoms) {
		iBoxD++;
	}

	// determine the size of the bins
	fBoxD = 3.55;           // Density??
	*box = iBoxD*fBoxD;     // calculate the box dimension

	// add a particle in each bin
	atomCount = 0;
	for(x=0;x<iBoxD;x++) {
		xPos = (x+0.5)*fBoxD;
		for(y=0;y<iBoxD;y++) {
			yPos = (y+0.5)*fBoxD;
			for(z=0; z<iBoxD; z++) {
				if (atomCount < nAtoms) {
					zPos = (z+0.5)*fBoxD;
					coord[atomCount][0]=xPos;
					coord[atomCount][1]=yPos;
					coord[atomCount][2]=zPos;
					atomCount++;
				} else {
					break;
				}
			}

			if (atomCount>=nAtoms) {
				break;
			}
		}
	}
}

