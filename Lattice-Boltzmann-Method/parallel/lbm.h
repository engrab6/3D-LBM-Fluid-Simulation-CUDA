#if !defined(LBM_H)
#define      LBM_H

#include <petscmat.h>
#include <petscsys.h>
#include <stddef.h>
#include <mpi.h>
#include <driver_types.h>


#if defined(__cplusplus)
extern "C" {
#endif

typedef struct _lbm *LBM;

/** Each device gets a 3D block of the indices.  N is the number of points per
 * dimension, so each of the *Start and *End pairs should be a subinterval of
 * [0, N], and indicates that this block goes from coordinates
 * [*Start, ..., *End - 1] in that direction.

 The indexing is in lexicographic ordering, so the global index of (i, j, k)
 is

 i + N * (j + (N * k),

 whereas the *local* index of the same (i, j, k) is

 i - ib.xStart + (ib.xEnd - ib.xStart) * (j - ib.yStart + (ib.yEnd - yStart) * (k - ib.zStart)

 You are responsible for understanding how the to convert between local and
 global indexes for your choice of partition.
 */
struct IndexBlock
{
  size_t xStart, xEnd;
  size_t yStart, yEnd;
  size_t zStart, zEnd;
};

__constant__ int gpu_lbm_velocities[27][3] = {{-1, -1, -1}, { 0, -1, -1}, {+1, -1, -1},
                                          {-1,  0, -1}, { 0,  0, -1}, {+1,  0, -1},
                                          {-1,  1, -1}, { 0,  1, -1}, {+1,  1, -1},

                                          {-1, -1,  0}, { 0, -1,  0}, {+1, -1,  0},
                                          {-1,  0,  0}, { 0,  0,  0}, {+1,  0,  0},
                                          {-1,  1,  0}, { 0,  1,  0}, {+1,  1,  0},

                                          {-1, -1, +1}, { 0, -1, +1}, {+1, -1, +1},
                                          {-1,  0, +1}, { 0,  0, +1}, {+1,  0, +1},
                                          {-1,  1, +1}, { 0,  1, +1}, {+1,  1, +1}};

static int const lbm_velocities[27][3] = {{-1, -1, -1}, { 0, -1, -1}, {+1, -1, -1},
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

__constant__ double gpu_lbm_weights[27] = {lbm_wv, lbm_we, lbm_wv,
                                 lbm_we, lbm_wf, lbm_we,
                                 lbm_wv, lbm_we, lbm_wv,

                                 lbm_we, lbm_wf, lbm_we,
                                 lbm_wf, lbm_wc, lbm_wf,
                                 lbm_we, lbm_wf, lbm_we,

                                 lbm_wv, lbm_we, lbm_wv,
                                 lbm_we, lbm_wf, lbm_we,
                                 lbm_wv, lbm_we, lbm_wv};

static double const lbm_weights[27] = {lbm_wv, lbm_we, lbm_wv,
                                 lbm_we, lbm_wf, lbm_we,
                                 lbm_wv, lbm_we, lbm_wv,

                                 lbm_we, lbm_wf, lbm_we,
                                 lbm_wf, lbm_wc, lbm_wf,
                                 lbm_we, lbm_wf, lbm_we,

                                 lbm_wv, lbm_we, lbm_wv,
                                 lbm_we, lbm_wf, lbm_we,
                                 lbm_wv, lbm_we, lbm_wv};

/** Create a Lattice-Boltzmann method iteration context object
 *
 * This object encapsulates all the decisions about how best to implement LBM
 * for various number of processors and problem sizes
 *
 * comm: input, MPI communicator describing the parallel layout.
 * out:  lbm_p, the object that will define the graph on which we will
 * search.
 */
int LBMCreate(MPI_Comm comm, LBM *lbm_p);
/** Destroy the object */
int LBMDestroy(LBM *lbm_p);

/** Get arrays for the fields to be filled.

  \param[in]     N             The number of points per dimension
  \param[in]     numDevices    The number of cuda devices available
  \param[out]    iBlock        A description of the index block for each cuda device
  \param[out]    fLocal        For each cuda device, this should be made to
                               point to an array, allocated by you *on the device*, where you would like
                               the local solution vector to be stored (can be NULL if the device is empty)
*/
int LBMGetFieldArrayDevice(LBM lbm, size_t N, int numDevices, struct IndexBlock iBlock[], PetscReal (*fLocal_p[])[27]);

/** Restore arrays for the fields to be filled.

  \param[in]     N             The number of points per dimension
  \param[in]     numDevices    The number of cuda devices available
  \param[in]     iBlock        A description of the index block for each cuda device
  \param[out]    fLocal        For each cuda device, you should free the device data that you allocated
*/
int LBMRestoreFieldArrayDevice(LBM lbm, size_t N, int numDevices, struct IndexBlock iBlock[], PetscReal (*fLocal_p[])[27]);

/** Lattice-Boltzmann BGK Method for Newtonian fluid without forcing.

  The grid layout is described above in the definition of indexBlock.  The grid is periodic.

  See the serial problem for the description about what this calculation means and the correct order to perform the steps in.

  \param[in] N          The number of grid points in each direction (the domain is periodic in each direction).
  \param[in] num_iter   The number of time steps to simulate.
  \param[in] tau        The relaxation time.
  \param[in] numDevices The number of cuda devices.
  \param[in] iBlock     A description of the layout on each device
  \param[in/out] f      The densities in each of the 27 discrete velocities on each device.
 */
int LBMIterateDevice(LBM lbm, size_t N, int num_iter, double tau, int numDevices, const struct IndexBlock iBlock[], PetscReal (*fLocal[])[27], int block_size);

#if defined(__cplusplus)
}
#endif
#endif
__global__ void set_f_temp(PetscReal (*f)[27], double (*f_temp)[27], size_t n);
__global__ void run_lbm(PetscReal (*f)[27], double (*f_temp)[27], size_t n, double tau, struct IndexBlock iBlock);
void printDistribution(PetscReal (*f)[27], size_t n);