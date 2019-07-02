#include "lbm.h"
#include <petscsys.h>
#include <driver_types.h>

struct _lbm
{
  MPI_Comm comm;
};

int LBMCreate(MPI_Comm comm, LBM *lbm_p)
{
  LBM lbm = NULL;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = PetscCalloc1(1,&lbm); CHKERRQ(ierr);

  lbm->comm = comm;

  *lbm_p = lbm;
  PetscFunctionReturn(0);
}

__global__

int LBMDestroy(LBM *lbm_p)
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;

  ierr = PetscFree(*lbm_p); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
