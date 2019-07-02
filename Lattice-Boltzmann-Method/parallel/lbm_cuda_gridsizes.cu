
#include <petscsys.h>
#include <driver_types.h>
#include "lbm.h"
#include <math.h>
#include <time.h>
#include "seconds.h" // seconds.cpp function sourced from Krüger Timm, et al. The Lattice Boltzmann Method: Principles and Practice. Springer International Publishing, 2017.


#define CUDA_CHK(cerr) do {cudaError_t _cerr = (cerr); if ((_cerr) != cudaSuccess) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_LIB,"Cuda error");} while(0)
// #define PRINT_OUTPUT 1
  // __constant__ int gpu_lbm_velocities[27][3]; 
  // __constant__ double gpu_lbm_weights[27];

int LBMGetFieldArrayDevice(LBM lbm, size_t N, int numDevices,
                           struct IndexBlock iBlock[],
                           PetscReal (*fLocal_p[])[27])
{
  cudaError_t       cerr;

  PetscFunctionBeginUser;

  iBlock[0].xStart = 0;
  iBlock[0].xEnd   = N;
  iBlock[0].yStart = 0;
  iBlock[0].yEnd   = N;
  iBlock[0].zStart = 0;
  iBlock[0].zEnd   = N;
  cerr = cudaMalloc(&fLocal_p[0], N * N * N * sizeof(PetscReal[27])); CUDA_CHK(cerr);

  for (int i = 1; i < numDevices; i++) {
    iBlock[i].xStart = N;
    iBlock[i].xEnd   = N;
    iBlock[i].yStart = N;
    iBlock[i].yEnd   = N;
    iBlock[i].zStart = N;
    iBlock[i].zEnd   = N;
    fLocal_p[i] = NULL;
  }

  PetscFunctionReturn(0);
}

int LBMRestoreFieldArrayDevice(LBM lbm, size_t N, int numDevices,
                               struct IndexBlock iBlock[],
                               PetscReal (*fLocal_p[])[27])
{
  cudaError_t    cerr;

  PetscFunctionBeginUser;
  cerr = cudaFree(fLocal_p[0]); CUDA_CHK(cerr);
  PetscFunctionReturn(0);
}

void printDistribution(PetscReal (*f)[27], size_t n){ 
  for (size_t i = 0; i < n; i++) {
    for (int dir = 0; dir < 27; dir++) {
      printf("%3.2f ", f[i][dir]); 
    }
    printf("\n");
  }
}

int LBMIterateDevice(LBM lbm, size_t N, int num_iter, double tau,
                     int numDevices,
                     const struct IndexBlock iBlock[],
                     PetscReal (*fLocal[])[27], int threads_per_block) {
  
  PetscFunctionBeginUser;
  cudaError_t cerr, ierr;  
  struct cudaDeviceProp prop;  
  ierr = cudaGetDeviceProperties(&prop, 0); CHKERRQ(ierr);

  printf("    Initializing LBM calculation...\n");   
  printf("    Total points = %lu * %lu * %lu = %lu\n", N, N, N, N*N*N);

  // dim3 grid_of_blocks(N, 1, 1); 
  // dim3 block_of_threads(N, N, 1);

  // int threads_per_block = PetscMin(prop.maxThreadsPerBlock, N);
  threads_per_block = PetscMin(prop.maxThreadsPerBlock, threads_per_block);
  int numBlocks = ceil(N*N*N / threads_per_block); 
  dim3 grid_of_blocks(numBlocks, 1, 1); 
  dim3 block_of_threads(threads_per_block, 1, 1); 
  printf("    Allocated grid of blocks with dimensions: (%d, %d, %d); %d blocks total\n", grid_of_blocks.x, grid_of_blocks.y, grid_of_blocks.z, grid_of_blocks.x * grid_of_blocks.y * grid_of_blocks.z);
  printf("    Allocated %d thread blocks with dimensions: (%d, %d, %d); %d threads per block\n", grid_of_blocks.x * grid_of_blocks.y * grid_of_blocks.z, block_of_threads.x, block_of_threads.y, block_of_threads.z, block_of_threads.x * block_of_threads.y * block_of_threads.z);
  printf("    %d blocks * %d threads per block = %d threads total\n", grid_of_blocks.x * grid_of_blocks.y * grid_of_blocks.z, block_of_threads.x * block_of_threads.y * block_of_threads.z, (block_of_threads.x * block_of_threads.y * block_of_threads.z)*(grid_of_blocks.x * grid_of_blocks.y * grid_of_blocks.z)); 
  
  double (*f_temp)[27];
  cerr = cudaMalloc(&f_temp, N * N * N * sizeof(double[27])); CUDA_CHK(cerr);
  double bytesPerMiB = 1024.0*1024.0;
  double bytesPerGiB = 1024.0*1024.0*1024.0;
  

  size_t gpu_free_mem, gpu_total_mem;
  (cudaMemGetInfo(&gpu_free_mem,&gpu_total_mem));
  
  printf("Remote CUDA information: \n");
  printf("Device name: %s\n",prop.name);
  printf("# of multiprocessors: %d\n",prop.multiProcessorCount);
  printf("compute capability: %d.%d\n",prop.major,prop.minor);
  printf("GPU global memory available: %.1f MiB\n\n",prop.totalGlobalMem/bytesPerMiB);


  // MAIN SIMULATION LOOP
  double begin = seconds();
  for (int i = 0; i < num_iter; i++) {
    set_f_temp <<< grid_of_blocks, block_of_threads >>> (fLocal[0], f_temp, N); // set f_temp = f //TODO: move this into run_lbm
    run_lbm <<<grid_of_blocks, block_of_threads>>> (fLocal[0], f_temp, N, tau, iBlock[0]); // *fLocal[0] is an array of length 27, with n^3 rows
    cudaDeviceSynchronize();
  }
    double end = seconds();
    double runtime = end - begin;
    // double gpu_runtime = 0.001*0.0f;

    size_t doubles_read = 27; // per node every time step
    size_t doubles_written = 27;
    size_t doubles_saved = 0; // per node every NSAVE time steps
    
    // note NX*NY overflows when NX=NY=65536
    size_t nodes_updated = num_iter*size_t(N*N*N);
    size_t nodes_saved   = 0;
    double speed = nodes_updated/(1e6*runtime);
    
    double bandwidth = (nodes_updated*(doubles_read + doubles_written)+nodes_saved*(doubles_saved))*sizeof(double)/(runtime*bytesPerGiB);
    
    printf("Performance Statistics: \n");
    printf("GPU Memory Allocated: %.1f (MiB)\n",432*N*N*N/bytesPerMiB);
    printf("Runtime: %.3f (s)\n",runtime);
    printf("Speed: %.2f (Mlups)\n",speed);
    printf("Bandwidth: %.1f (GiB/s)\n",bandwidth);  
#if PRINT_OUTPUT==1   
    printf("Printing number_of_points,runtime (s) ,speed (Mlups),bandwidth (GiB/s)\n");
    char filename[100];
    sprintf(filename, "output_N=%d_%d_iterations.txt", N, num_iter); 
    printf("Filename = %s\n", filename);
    FILE        *fp = NULL;
    fp = fopen(filename,"a");
    fprintf(fp,"%d,%.3f,%.4f,%.4f\n", N*N*N, runtime, speed, bandwidth);
    fclose(fp);
#endif


  PetscFunctionReturn(0);
}

__global__ void set_f_temp(PetscReal (*f)[27], double (*f_temp)[27], size_t n) {
  // int x = threadIdx.x; 
  // int y = blockIdx.x; 
  // int z = blockIdx.y; 
  // int threadNum = threadIdx.x;
  // int blockNum = blockIdx.x; 
  int pos_index = threadIdx.x + blockIdx.x * blockDim.x;
  // int y = x + 1; 
  // int z = x + 2;
  // int pos_index = x + n * (y + n * z); // current (x, y, z)   
  for (int i = 0; i < 27; i++) {
    f_temp[pos_index][i] = f[pos_index][i]; 
  }
}

__global__ void run_lbm(PetscReal (*f)[27], double (*f_temp)[27], size_t n, double tau, struct IndexBlock iBlock) 
{
  const double tau_reciprocal = 1/tau; 
  const double one_minus_tau_reciprocal = 1 - tau_reciprocal;  
  // int x = threadIdx.x; 
  // int y = blockIdx.x; 
  // int z = blockIdx.y; 
  int pos_index = threadIdx.x + blockIdx.x * blockDim.x;

  // int pos_index = x + n * (y + n * z); // current (x, y, z) corresponds to a pos_index in the global f array

  double rho = 0, u = 0, v = 0, w = 0;

  for (int dir = 0; dir < 27; dir++) {
    size_t x_periodic = (n + x - gpu_lbm_velocities[dir][0]) % n; 
    size_t y_periodic = (n + y - gpu_lbm_velocities[dir][1]) % n;
    size_t z_periodic = (n + z - gpu_lbm_velocities[dir][2]) % n;
    // Convert Cartesian coordinates to n*n*n row index:
    size_t old_pos = x_periodic + y_periodic * n + z_periodic * n * n;         

    f[pos_index][dir] = f_temp[old_pos][dir]; 
    rho += f[pos_index][dir];
    u += f[pos_index][dir] * gpu_lbm_velocities[dir][0]; 
    v += f[pos_index][dir] * gpu_lbm_velocities[dir][1]; 
    w += f[pos_index][dir] * gpu_lbm_velocities[dir][2];     
  }

  double rhoinv = 1.0/rho; 

  u *= rhoinv; v *= rhoinv; w *= rhoinv;
  for (int dir = 0; dir < 27; dir++) {
    // calculate u_alpha (the dot product)
    double u_alpha = gpu_lbm_velocities[dir][0] * u +
      gpu_lbm_velocities[dir][1] * v + gpu_lbm_velocities[dir][2] * w;
    // calculate equilibrium distribution
    double f_equil_i = gpu_lbm_weights[dir] * rho * 
      (1. + 3. * u_alpha + 4.5 * u_alpha * u_alpha - 1.5*(u*u + v*v + w*w));
    f[pos_index][dir] = tau_reciprocal * f_equil_i + one_minus_tau_reciprocal * f[pos_index][dir];
  }
}