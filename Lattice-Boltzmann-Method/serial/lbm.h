#if !defined(LBM_H)
#define      LBM_H

#include <stddef.h>
#include <petscmat.h>

static const int lbm_velocities[27][3] = {{-1, -1, -1}, { 0, -1, -1}, {+1, -1, -1},
                                          {-1,  0, -1}, { 0,  0, -1}, {+1,  0, -1},
                                          {-1,  1, -1}, { 0,  1, -1}, {+1,  1, -1},

                                          {-1, -1,  0}, { 0, -1,  0}, {+1, -1,  0},
                                          {-1,  0,  0}, { 0,  0,  0}, {+1,  0,  0},
                                          {-1,  1,  0}, { 0,  1,  0}, {+1,  1,  0},

                                          {-1, -1, +1}, { 0, -1, +1}, {+1, -1, +1},
                                          {-1,  0, +1}, { 0,  0, +1}, {+1,  0, +1},
                                          {-1,  1, +1}, { 0,  1, +1}, {+1,  1, +1}};

#if 1
#define lbm_wc (8./27.)
#define lbm_wf (2./27.)
#define lbm_we (1./54.)
#define lbm_wv (1./216.)
#else
#define lbm_wc (2./9.)
#define lbm_wf (1./9.)
#define lbm_we (0.)
#define lbm_wv (1./72.)
#endif

static double lbm_weights[27] = {lbm_wv, lbm_we, lbm_wv,
                                 lbm_we, lbm_wf, lbm_we,
                                 lbm_wv, lbm_we, lbm_wv,

                                 lbm_we, lbm_wf, lbm_we,
                                 lbm_wf, lbm_wc, lbm_wf,
                                 lbm_we, lbm_wf, lbm_we,

                                 lbm_wv, lbm_we, lbm_wv,
                                 lbm_we, lbm_wf, lbm_we,
                                 lbm_wv, lbm_we, lbm_wv};

/** Lattice-Boltzmann BGK Method for Newtonian fluid without forcing.

  There are n x n x n grid points.  The grid points (i, j, k) has index
  $p(i,j,k) = i + n * (j + n * k)$, so that the x-coordinate varies first, then the
  y-coordinate, then the z-coordinate.  The grid is periodic.

  Each grid point has 27 densities `f[p][alpha]` associated with 27 discrete
  velocities.  The velocities are given in `lbm_velocities`.  We consider the
  time step to be 1, so that the velocities also correspond to the
  displacements between grid points:

      (i, j, k) + lbm_velocities[alpha] * dt =
      (i + lbm_velocities[alpha][0], j + lbm_velocities[alpha][1], k + lbm_velocities[alpha][2]).

      (Let's abbreviate this calculation p + v[alpha])

  At each point, those densities define macroscopic density and velocities by

  rho[p] = \sum_{alpha} f[p][alpha]
  u[p][0] = (1/rho[p]) * \sum_{alpha} f[p][alpha] * lbm_velocities[alpha][0]
  u[p][1] = (1/rho[p]) * \sum_{alpha} f[p][alpha] * lbm_velocities[alpha][1]
  u[p][2] = (1/rho[p]) * \sum_{alpha} f[p][alpha] * lbm_velocities[alpha][2]

  The macroscopic velocity distributes to the quantized directions by

  u_d[p][alpha] = sum_{d=0}^3 u[p][d] * lbm_velocities[alpha][d]

  These values determine the equilibrium densities f_0[p][alpha]:

  f_0[p][alpha] = rho[p] * lbm_weights[alpha] * (1 + 3 u_d[p][alpha] + \frac{9}{2} u_d[p][alpha]^2 - \frac{3}{2} (u[p][0]^2 + u[p][1]^2 + u[p][2]^2)).

  Finally, the new densities at the point downstream of

  f_new[p + v[alpha]][alpha] = (1/tau) * f_0[p][alpha] + (1 - 1/tau) * f[p][alpha]

  \param[in] n        The number of grid points in each direction (the domain is periodic in each direction).
  \param[in] n_iter   The number of time steps to simulate.
  \param[in] tau      The relaxation time.
  \param[in/out] f    [n x n x n x 27 vector] The densities in each of the 27
                      discrete velocities.  The grid points a ordered lexicographically as described above.
 */
int lattice_boltzmann_method(size_t n, int n_iter, double tau, double (*f)[27]);

#endif