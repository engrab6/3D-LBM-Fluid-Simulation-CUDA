void stream(double** f_temp, double (*f)[27], size_t n) {
	for (size_t z = 0, pos_index = 0; z < n; z++) {
		for (size_t y = 0; y < n; y++) {
			for (size_t x = 0; x < n; x++, pos_index++) {
				// For each direction at this location, f_i = f_i_star(previous_location, previous_time)
				for (int dir = 0; dir < 27; dir++) {
					// Calculate old positions, and apply periodic boundary condition:
					size_t x_periodic = (n + x - lbm_velocities[dir][0]) % n; 
					size_t y_periodic = (n + y - lbm_velocities[dir][1]) % n;
					size_t z_periodic = (n + z - lbm_velocities[dir][2]) % n;
					// Convert Cartesian coordinates to n*n*n row index: TODO: CHECK THIS!!!
					size_t old_pos = x_periodic + y_periodic * n + z_periodic * n * n; 

					// f_temp at each position is 
					f[pos_index][dir] = f_temp[old_pos][dir]; 
				}
			}
		}
	}
}

void momentUpdate(double** f, double* rho, double** vel, size_t n) {
	for (size_t z = 0, pos_index = 0; z < n; z++) {
		for (size_t y = 0; y < n; y++) {
			for (size_t x = 0; x < n; x++, pos_index++) {
				double u = 0, v = 0, w = 0;
				rho[pos_index] = f[pos_index][0] + f[pos_index][1] + f[pos_index][2] + f[pos_index][3] + f[pos_index][4] +
							 	 f[pos_index][5] + f[pos_index][6] + f[pos_index][7] + f[pos_index][8] + f[pos_index][9] +
							 	 f[pos_index][10] + f[pos_index][11] + f[pos_index][12] + f[pos_index][13] + f[pos_index][14] +
							 	 f[pos_index][15] + f[pos_index][16] + f[pos_index][17] + f[pos_index][18] + f[pos_index][19] +
							 	 f[pos_index][20] + f[pos_index][21] + f[pos_index][22] + f[pos_index][23] + f[pos_index][24] +
							 	 f[pos_index][25] + f[pos_index][26];
				for (int dir = 0; dir < 27; dir++) {				
					u += f[pos_index][dir] * lbm_velocities[dir][0]; 
					v += f[pos_index][dir] * lbm_velocities[dir][1]; 
					w += f[pos_index][dir] * lbm_velocities[dir][2];
				}
				vel[pos_index][0] = u/rho[pos_index]; 
				vel[pos_index][1] = v/rho[pos_index]; 
				vel[pos_index][2] = w/rho[pos_index];
			}
		}
	}
}

void collision(double** f, double* rho, double** vel, size_t n, double tau) {
	const double tau_reciprocal = 1/tau; 
	const double one_minus_tau_reciprocal = 1 - tau_reciprocal;
	for (size_t z = 0, counter = 0; z < n; z++) {
		for (size_t y = 0; y < n; y++) {
			for (size_t x = 0; x < n; x++, counter++) {
				for (int dir = 0; dir < 27; dir++) {
					// calculate u_alpha (the dot product)
					double u = vel[counter][0]; 
					double v = vel[counter][1];
					double w = vel[counter][2];
					double u_alpha = lbm_velocities[dir][0] * u +
						lbm_velocities[dir][1] * v + lbm_velocities[dir][2] * w;
					// calculate equilibrium distribution
					double f_equil_i = lbm_weights[dir] * rho[counter] * 
						(1. + 3. * u_alpha + 4.5 * u_alpha * u_alpha - 1.5*(u*u + v*v + w*w));
					f[counter][dir] = tau_reciprocal * f_equil_i + one_minus_tau_reciprocal * f[counter][dir];
				}
				// set rho and velocity to zero since we're not using them anymore
				rho[counter] = vel[counter][0] = vel[counter][1] = vel[counter][2] = 0;
			}
		}
	}

}
size_t getOldPos(size_t x, size_t y, size_t z, size_t n, size_t dir) {
	// Calculate old positions, and apply periodic boundary condition:
	size_t x_periodic = (n + x - lbm_velocities[dir-1][0]) % n; 
	size_t y_periodic = (n + y - lbm_velocities[dir-1][1]) % n;
	size_t z_periodic = (n + z - lbm_velocities[dir-1][2]) % n;
	// Convert Cartesian coordinates to n*n*n row index: TODO: CHECK THIS!!!
	return x_periodic + y_periodic * n + z_periodic * n * n; 	
}

void combined(double (*f)[27], double **f_temp, double* rho, double** vel, size_t n, double tau)
{
	const double tau_reciprocal = 1/tau; 
	const double one_minus_tau_reciprocal = 1 - tau_reciprocal;  
	const double tau_reciprocal = 1/tau; 
	const double one_minus_tau_reciprocal = 1 - tau_reciprocal;
	for (size_t z = 0, counter = 0, pos_index=0; z < n; z++) {
		for (size_t y = 0; y < n; y++) {
			for (size_t x = 0; x < n; x++, counter++, pos_index++) {
  for (int dir = 0; dir < 27; dir++) {
    size_t x_periodic = (n + x - lbm_velocities[dir][0]) % n; 
    size_t y_periodic = (n + y - lbm_velocities[dir][1]) % n;
    size_t z_periodic = (n + z - lbm_velocities[dir][2]) % n;
    // Convert Cartesian coordinates to n*n*n row index: TODO: CHECK THIS!!!
    size_t old_pos = x_periodic + y_periodic * n + z_periodic * n * n;         
    f_temp[dir] = f[old_pos][dir]; 
    rho += f_temp[dir];
  } 
  double rhoinv = 1.0/rho;
  double u = 0, v = 0, w = 0;
  for (int dir = 0; dir < 27; dir++) {
    // f_temp at each position is 
    f[pos_index][dir] = f_temp[dir];    
    u += f[pos_index][dir] * lbm_velocities[dir][0]; 
    v += f[pos_index][dir] * lbm_velocities[dir][1]; 
    w += f[pos_index][dir] * lbm_velocities[dir][2];         
  }
  u *= rhoinv; v *= rhoinv; w *= rhoinv;
  for (int dir = 0; dir < 27; dir++) {
    // calculate u_alpha (the dot product)
    double u_alpha = lbm_velocities[dir][0] * u +
      lbm_velocities[dir][1] * v + lbm_velocities[dir][2] * w;
    // calculate equilibrium distribution
    double f_equil_i = lbm_weights[dir] * rho * (1. + 3. * u_alpha + 4.5 * u_alpha * u_alpha - 1.5*(u*u + v*v + w*w));
   	f[pos_index][dir] = tau_reciprocal * f_equil_i + one_minus_tau_reciprocal * f[pos_index][dir];
  }
			}
		}
	}

}


void swap_deepCopy(double (*f)[27], double** f_temp, size_t n) {
	for (int i = 0; i < n*n*n; i++) {
		for (int j = 0; j < 27; j++) {
			f_temp[i][j] = f[i][j];
		}
	}
}