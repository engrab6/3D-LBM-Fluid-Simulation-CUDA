const char help[] = "Test driver for the correctness of Lattice-Boltzmann Method implementation";

#include <petscmat.h>
#include "lbm.h"

static PetscErrorCode TaylorGreen(PetscInt n, PetscReal t, PetscReal scale, PetscReal nu, PetscInt i, PetscInt j, PetscReal *u, PetscReal *v)
{
  PetscFunctionBegin;
  *u =   PetscExpReal(-2. * t * nu) * scale * PetscSinReal(2 * PETSC_PI * (i + 0.5) / n) * PetscCosReal(2 * PETSC_PI * (j + 0.5) / n);
  *v = - PetscExpReal(-2. * t * nu) * scale * PetscCosReal(2 * PETSC_PI * (i + 0.5) / n) * PetscSinReal(2 * PETSC_PI * (j + 0.5) / n);
  PetscFunctionReturn(0);
}

int main(int argc, char **argv)
{
  PetscInt       test, numTests = 10;
  PetscInt       scale = 12;
  PetscRandom    rand;
  MPI_Comm       comm;
  PetscViewer    viewer;
  PetscBool      write_pre, write_post;
  char           pre_vtk[BUFSIZ] = {0};
  char           post_vtk[BUFSIZ] = {0};
  PetscErrorCode ierr;

  ierr = PetscInitialize(&argc, &argv, NULL, help); if (ierr) return ierr;
  comm = PETSC_COMM_WORLD;
  viewer = PETSC_VIEWER_STDOUT_WORLD;
  ierr = PetscOptionsBegin(comm, NULL, "Lattice Boltzmann Method Test Options", "test_lbm.c");CHKERRQ(ierr);
  ierr = PetscOptionsInt("-num_tests", "Number of tests to run", "test_lbm.c", numTests, &numTests, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-scale", "Scale (log2) of the array in the test", "test_lbm.c", scale, &scale, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsString("-pre_vtk", "File name for vtk output of Taylor-Green initial condition", "test_lbm.c", pre_vtk, pre_vtk, BUFSIZ, &write_pre);CHKERRQ(ierr);
  ierr = PetscOptionsString("-post_vtk", "File name for vtk output of Taylor-Green final condition", "test_lbm.c", post_vtk, post_vtk, BUFSIZ, &write_post);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);

  ierr = PetscRandomCreate(comm, &rand);CHKERRQ(ierr);
  ierr = PetscRandomSetFromOptions(rand);CHKERRQ(ierr);
  ierr = PetscRandomSeed(rand);CHKERRQ(ierr);

  ierr = PetscViewerASCIIPrintf(viewer, "Running %D tests of gauss_seidel_iteration()\n", numTests);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPushTab(PETSC_VIEWER_STDOUT_WORLD);
  for (test = 0; test < numTests; test++) {
    PetscReal   diff, norm, nreal, dt, dx, nu, tau, tol;
    PetscInt    dir, n, nits, i, j, k, p, a, d;
    PetscReal   (*f)[27];
    FILE        *fp = NULL;

    ierr = PetscViewerASCIIPrintf(viewer, "Test %D:\n", test);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPushTab(PETSC_VIEWER_STDOUT_WORLD);
    ierr = PetscRandomSetInterval(rand, PetscPowReal(2., scale / 3. - 1), PetscPowReal(2., scale / 3.));CHKERRQ(ierr);
    ierr = PetscRandomGetValueReal(rand, &nreal);CHKERRQ(ierr);

    n = (PetscInt) nreal;
    nits = n * n;
    nu = 1.;
    dx = (2. * PETSC_PI) / n;
    dt = 0.1 * dx * dx;
    tau = (1./2.) + 3 * nu * dt / (dx * dx);

    tol = 200. / (n * n);

    ierr = PetscViewerASCIIPrintf(viewer, "Test dimensions: [%D x %D x %D], %D iterations, tolerance %g\n", n, n, n, nits, (double) tol);CHKERRQ(ierr);

    ierr = PetscMalloc1(n*n*n, &f);CHKERRQ(ierr);

    /* Create a Taylor-Green instability oriented in varying directions */
    dir = test % 3;
    if (write_pre) {
      fp = fopen(pre_vtk,"w");
      if (!fp) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_LIB,"Could not open %s for writing\n",pre_vtk);
      fprintf(fp,"# vtk DataFile Version 2.0\n2D Taylor-Green vortex\nASCII\nDATASET STRUCTURED_POINTS\n");
      fprintf(fp,"DIMENSIONS %d %d %d\nASPECT_RATIO 1 1 1\nORIGIN 0.5 0.5 0.5\nPOINT_DATA %d\n",(int)n,(int)n,(int)n,(int)(n*n*n));
      fprintf(fp,"VECTORS u double\n");
    }
    for (k = 0, p = 0; k < n; k++) {
      for (j = 0; j < n; j++) {
        for (i = 0; i < n; i++, p++) {
          PetscReal U[3] = {0.}, u2;
          int       ijk[3];
          PetscReal u, v;

          ijk[0] = i;
          ijk[1] = j;
          ijk[2] = k;

          ierr = TaylorGreen(n, 0., 0.01, nu, ijk[dir], ijk[(dir + 1) % 3], &u, &v);CHKERRQ(ierr);
          U[ dir       ] = u;
          U[(dir+1) % 3] = v;

          if (fp) fprintf(fp,"%e %e %e\n", U[0], U[1], U[2]);

          /* Convert to lattice boltzmann dimensions */
          U[0] *= dt / dx;
          U[1] *= dt / dx;
          U[2] *= dt / dx;

          u2 = U[0]*U[0] + U[1]*U[1] + U[2]*U[2];

          for (a = 0; a < 27; a++) {
            PetscReal u_d = 0.;

            for (d = 0; d < 3; d++) {
              u_d += U[d] * lbm_velocities[a][d];
            }
            f[p][a] = lbm_weights[a] * (1. + 3. * u_d + (9./2.) * (u_d * u_d) - (3./2.) * u2);
            if (f[p][a] < 0.) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_LIB, "Negative weight");
          }
        }
      }
    }
    if (fp) fclose(fp);
    fp = NULL;

    ierr = lattice_boltzmann_method((size_t) n, nits, tau, f);CHKERRQ(ierr);

    diff = 0.;
    norm = 0.;
    dir = test % 3;
    if (write_post) {
      fp = fopen(post_vtk,"w");
      if (!fp) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_LIB,"Could not open %s for writing\n",post_vtk);
      fprintf(fp,"# vtk DataFile Version 2.0\n2D Taylor-Green vortex\nASCII\nDATASET STRUCTURED_POINTS\n");
      fprintf(fp,"DIMENSIONS %d %d %d\nASPECT_RATIO 1 1 1\nORIGIN 0.5 0.5 0.5\nPOINT_DATA %d\n",(int)n,(int)n,(int)n,(int)(n*n*n));
      fprintf(fp,"VECTORS u double\n");
    }
    for (k = 0, p = 0; k < n; k++) {
      for (j = 0; j < n; j++) {
        for (i = 0; i < n; i++, p++) {
          PetscReal rho  = 0.;
          PetscReal U[3] = {0.};
          PetscReal Ucheck[3] = {0.};
          int       ijk[3];
          PetscReal u, v;

          for (a = 0; a < 27; a++) {
            rho += f[p][a];
            for (d = 0; d < 3; d++) {
              U[d] += f[p][a] * lbm_velocities[a][d];
            }
          }
          if (rho <= 0.) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_LIB,"non-positive density");
          for (d = 0; d < 3; d++) {
            U[d] /= rho;
          }

          /* Convert to physical dimensions */
          U[0] *= dx / dt;
          U[1] *= dx / dt;
          U[2] *= dx / dt;

          if (fp) fprintf(fp,"%e %e %e\n", U[0], U[1], U[2]);

          ijk[0] = i;
          ijk[1] = j;
          ijk[2] = k;

          ierr = TaylorGreen(n, nits * dt, 0.01, nu, ijk[dir], ijk[(dir + 1) % 3], &u, &v);CHKERRQ(ierr);
          Ucheck[ dir       ] = u;
          Ucheck[(dir+1) % 3] = v;

          for (d = 0; d < 3; d++) {
            diff += PetscSqr(Ucheck[d] - U[d]);
            norm += PetscSqr(Ucheck[d]);
          }
        }
      }
    }
    if (fp) fclose(fp);
    ierr = PetscViewerASCIIPushTab(PETSC_VIEWER_STDOUT_WORLD);

    diff = PetscSqrtReal(diff/norm);
    ierr = PetscViewerASCIIPrintf(viewer, "Error %g\n", (double) diff);CHKERRQ(ierr);
    if (diff > tol) SETERRQ3(comm, PETSC_ERR_LIB, "Test %D failed residual test at threshold %g with value %g\n", test, (double) tol, (double) diff);

    ierr = PetscViewerASCIIPrintf(viewer, "Passed.\n");CHKERRQ(ierr);

    ierr = PetscFree(f);CHKERRQ(ierr);

    ierr = PetscViewerASCIIPopTab(PETSC_VIEWER_STDOUT_WORLD);
    ierr = PetscViewerASCIIPopTab(PETSC_VIEWER_STDOUT_WORLD);
  }
  ierr = PetscViewerASCIIPopTab(PETSC_VIEWER_STDOUT_WORLD);
  ierr = PetscRandomDestroy(&rand);CHKERRQ(ierr);

  ierr = PetscFinalize();
  return ierr;
}
