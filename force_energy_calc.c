/* file force_energy_calc.c */
// c code for force_energy_calc subroutine

#include <math.h>                                                     
#include "force_energy_calc.h"

double force_energy_calc(int nAtoms, double box, double cutoff2, double Ar_eps, double Ar_sigma6, double **coord, double **atomForces) {

	int atom1;
	int atom2;
	double dist2;
	double dist2d;
	double dist6d;
	double sigma_dist_6;
	int i;
	double component[3] = {0.0, 0.0, 0.0};
//	int atomCount;

//	component = (double**) malloc(3*sizeof(double*));
	
	double pot_energy;
	double ff;

	pot_energy = 0;

	for(atom1=0; atom1<nAtoms-1; atom1++){
		for(atom2=atom1+1;atom2<nAtoms;atom2++) {
			dist2 = 0;
			dist2d = 0;
			dist6d = 0;
			for(i=0;i<3;i++) {
				component[i] = 0.0;
				component[i] = coord[atom1][i] - coord[atom2][i];
				if(component[i] < -box/2.0) {
					component[i] += box;
				} else if(component[i] > box/2.0) {
					component[i] -= box;
				}
				dist2 += component[i]*component[i]; 		// r^2 = x^2 + y^2 + z^2; 
			}
			if(dist2 < cutoff2) {
				dist2d = 1.0/dist2;
				dist6d = dist2d*dist2d*dist2d;
				sigma_dist_6 = Ar_sigma6*dist6d;
				ff = 48*Ar_eps*dist2d*sigma_dist_6*(sigma_dist_6-0.5);
				for(i=0; i<3; i++) {
					atomForces[atom1][i]+= ff*component[i];
					atomForces[atom2][i]+= ff*component[i];
				}
				pot_energy += sigma_dist_6*(sigma_dist_6 - 1.0);
			}
		}
	}
	
	pot_energy = 4.0*Ar_eps*pot_energy;

//	IDEA for increasing performance; move the 48*Ar_eps calculation outside of cutoff loop; I can't determine if this will actually be quicker though...
//	
//	for(atomCount=0; atomCount<nAtoms; atomCount++) {
//		for(i=0;i<3;i++) {
//			
//		}
//	}

	return pot_energy;
}



