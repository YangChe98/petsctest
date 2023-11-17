#include </usr/lib/petscdir/petsc64-3.15/x86_64-linux-gnu-real/include/petscksp.h>
int min(int a, int b) {
  if (a < b) {
    return a;
  } else {
    return b;
  }
}
int main(int argc, char **args) {
  Vec x, c;       // x=exp(A*deltat)*c
  Mat A, Adeltat; // A , Adeltat=A*deltat
  PetscScalar deltat = 0.0025;
  PetscInt Nbasis = 16;
  PetscInt Istart, Iend, Ii, i, j, J, n = 16, m = 16;
  PetscScalar v;
  PetscInitialize(&argc, &args, (char *)0, (char *)0);
  PetscOptionsGetInt(NULL, NULL, "-Nbasis", &Nbasis, NULL);
  PetscPrintf(PETSC_COMM_WORLD, "18");

  MatCreate(PETSC_COMM_WORLD, &A);
  MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, Nbasis, Nbasis);
  MatMPIAIJSetPreallocation(A, 5, NULL, 5, NULL);

  MatSeqAIJSetPreallocation(A, 5, NULL);
  PetscPrintf(PETSC_COMM_WORLD, "24");
  // MatSeqSBAIJSetPreallocation(A, 1, 5, NULL);
  // MatMPISBAIJSetPreallocation(A, 1, 5, NULL, 5, NULL);
  // MatMPISELLSetPreallocation(A, 5, NULL, 5, NULL);
  // MatSeqSELLSetPreallocation(A, 5, NULL);
  MatGetOwnershipRange(A, &Istart, &Iend);

  PetscPrintf(PETSC_COMM_WORLD, "30");
  for (Ii = Istart; Ii < Iend; Ii++) {
    v = Ii * Ii;
    MatSetValues(A, 1, &Ii, 1, &Ii, &v, INSERT_VALUES);
  }
  PetscPrintf(PETSC_COMM_WORLD, "34");
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

  PetscPrintf(PETSC_COMM_WORLD, "38");

  VecCreate(PETSC_COMM_WORLD, &c);

  VecSetSizes(c, PETSC_DECIDE, Nbasis);

  VecSet(c, 1.0);
  /*
  MatAXPY(Adeltat, deltat, A, SAME_NONZERO_PATTERN);
  PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,
  PETSC_VIEWER_ASCII_MATLAB); MatView(Adeltat, PETSC_VIEWER_STDOUT_WORLD);
  PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);
  MatMult(Adeltat, c, x);
  */
  VecDestroy(&c);
  VecDestroy(&x);
  MatDestroy(&A);
  MatDestroy(&Adeltat);
  PetscFinalize();
  return 0;
}