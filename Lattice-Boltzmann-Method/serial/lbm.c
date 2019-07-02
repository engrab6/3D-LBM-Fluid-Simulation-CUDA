#include "lbm.h"
#include <petscmat.h>
#include "lbm_subroutines.h"


int lattice_boltzmann_method(size_t n, int n_iter, double tau, double (*f)[27]) // (*f)[27] = f[0][27]
{
	double** 	f_temp = malloc(n*n*n*sizeof(double*));
	double* 	rho = malloc(n*n*n*sizeof(double));
	double** 	vel = malloc(n*n*n*sizeof(double*));
	for (int i = 0; i < n*n*n; i++) {
		f_temp[i] = malloc(27*sizeof(double)); 
		vel[i] = malloc(3*sizeof(double));
	}
	//memcpy(f_temp, f, n*n*n*27*sizeof(double)); // set f_temp = f;
	
	for (int i = 0; i < n_iter; i++) { // O(n_iter)
		combined(f, f_temp, rho, vel, n, tau);

		// swap_deepCopy(f, f_temp, n); // set f_temp = f; O(n^3)
		// momentUpdate(f_temp, rho, vel, n); // gives macroscopic properties (rho & vel); need these for collision calc
		// collision(f_temp, rho, vel, n, tau); //f goes in, comes out as f* (f_temp = f*) 		
		// stream(f_temp, f, n); // propagate f* at each location to new location in f
		//writeToFile(&filename, f_temp, vel, rho, n); 
	}
	free(f_temp); 
	free(rho); 
	free(vel);
	return 0;
}