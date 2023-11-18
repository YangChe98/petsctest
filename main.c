#include </usr/lib/petscdir/petsc64-3.15/x86_64-linux-gnu-real/include/petscksp.h>
#include <malloc.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

int min(int a, int b) {
  if (a < b) {
    return a;
  } else {
    return b;
  }
}
void matrix_expontenial(Mat A, Vec c, double tdeltat, int m, int mmax, int smax,
                        double tol, Vec *Ac) {
  PetscScalar deltat = tdeltat;
  Mat Adeltat;
  int j;
  Vec c_tmp, c1, c1_tmp;
  VecDuplicate(c, &c_tmp);
  VecDuplicate(c, Ac);
  VecDuplicate(c, &c1_tmp);
  VecDuplicate(c, &c1);
  VecCopy(c, c_tmp);

  MatDuplicate(A, MAT_COPY_VALUES, &Adeltat);
  MatScale(Adeltat, deltat);

  MatMult(Adeltat, c_tmp, c1);

  // VecCopy(c1, c1_tmp);
  for (j = 2; j < (m + 2); j++) {
    MatMult(Adeltat, c1, c1_tmp);
    VecCopy(c1_tmp, c1);
  }

  int s;
  double nfac[mmax];
  nfac[0] = 1;
  for (j = 1; j < mmax; j++) {

    nfac[j] = nfac[j - 1] / (j + 1);
  }

  PetscScalar c1norm1;
  VecNorm(c1, NORM_1, &c1norm1);

  double c1norm = (double)c1norm1;
  // double tt = 0.9;
  //  PetscPrintf(PETSC_COMM_WORLD, "c1norm1=%.16f\n", c1norm);
  //  PetscPrintf(PETSC_COMM_WORLD, "c1norm1=%d\n", (int)tt);
  //  PetscPrintf(PETSC_COMM_WORLD, "c1norm1=%.200f\n", c1norm1 * nfac[m] /
  //  tol); PetscPrintf(PETSC_COMM_WORLD, "c1norm1=%.200f\n",
  //             PetscPowReal(c1norm1 * nfac[m] / tol, 1. / (m + 1)));
  s = (int)ceil(pow(c1norm * nfac[m] / tol, 1. / (m + 1)));
  float powtmp = pow(c1norm * nfac[m] / tol, 1. / (m + 1));
  double p = m * s;
  double s1, p1;
  int f = 0;
  while (f == 0 && m < mmax) {
    m++;
    MatMult(Adeltat, c1, c1_tmp);
    VecCopy(c1_tmp, c1);
    VecNorm(c1, NORM_1, &c1norm1);
    c1norm = (double)c1norm1;
    s1 = (int)ceil(pow(c1norm * nfac[m] / tol, 1. / (m + 1)));

    p1 = m * s1;
    if (p1 <= p) {
      p = p1;
      s = s1;
    } else {
      m = m - 1;
      f = 1;
    }
  }

  s = min(s, smax);
  p = m * s;
  PetscPrintf(PETSC_COMM_WORLD, "c1norm1=%d\n", s);
  Vec w;
  VecDuplicate(c, &w);
  VecCopy(c, w);

  VecCopy(c, c1);
  for (j = 0; j < m; j++) {
    // PetscPrintf(PETSC_COMM_WORLD, "c1norm1=%d\n", j);
    MatMult(Adeltat, c1, c1_tmp);
    VecCopy(c1_tmp, c1);
    VecAXPY(w, nfac[j] / pow(s, (j + 1)), c1);
  }

  Mat Adeltats;
  MatDuplicate(Adeltat, MAT_COPY_VALUES, &Adeltats);
  MatScale(Adeltats, 1. / s);

  int i;
  Vec v2, v2tmp;
  VecDuplicate(w, &v2);
  VecCopy(w, v2);
  VecDuplicate(w, &v2tmp);
  VecCopy(w, v2tmp);

  for (i = 0; i < (s - 1); i++) {
    VecCopy(w, v2);

    for (j = 0; j < m; j++) {
      MatMult(Adeltats, v2, v2tmp);
      VecCopy(v2tmp, v2);
      VecAXPY(w, nfac[j], v2);
    }
  }
  VecCopy(w, *Ac);

  MatDestroy(&Adeltat);
  VecDestroy(&c_tmp);
  VecDestroy(&c1);
  VecDestroy(&c1_tmp);
}
int main(int argc, char **args) {
  Vec Ac, c; // Ac=exp(A*deltat)*c
  Mat A;     // A , Adeltat=A*deltat
  PetscScalar deltat = 0.0025;
  int fd;
  PetscInt Nbasis = 16;
  PetscInt Istart, Iend, Ii, j, J, n = 16;
  int m = 50;
  int mmax = 100;
  int smax = 100;
  double tol = 1e-16;

  PetscScalar v;
  PetscViewer view_out, view_in;
  PetscViewer viewer, cviewer;
  PetscInitialize(&argc, &args, (char *)0, (char *)0);
  PetscOptionsGetInt(NULL, NULL, "-Nbasis", &Nbasis, NULL);

  MatCreate(PETSC_COMM_WORLD, &A);

  // MatSetType(A, MATSEQAIJ);
  MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, Nbasis, Nbasis);
  MatSetFromOptions(A);
  MatMPIAIJSetPreallocation(A, 5, NULL, 5, NULL);
  MatSeqAIJSetPreallocation(A, 5, NULL);
  //  MatSeqSBAIJSetPreallocation(A, 1, 5, NULL);
  //  MatMPISBAIJSetPreallocation(A, 1, 5, NULL, 5, NULL);
  //  MatMPISELLSetPreallocation(A, 5, NULL, 5, NULL);
  //  MatSeqSELLSetPreallocation(A, 5, NULL);

  // MatSeqSBAIJSetPreallocation(A, 1, 5, NULL);
  // MatMPISBAIJSetPreallocation(A, 1, 5, NULL, 5, NULL);
  // MatMPISELLSetPreallocation(A, 5, NULL, 5, NULL);
  // MatSeqSELLSetPreallocation(A, 5, NULL);
  MatGetOwnershipRange(A, &Istart, &Iend);

  for (Ii = Istart; Ii < Iend; Ii++) {
    v = (Ii + 1) * (Ii + 1);
    MatSetValues(A, 1, &Ii, 1, &Ii, &v, ADD_VALUES);
  }

  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

  VecCreate(PETSC_COMM_WORLD, &c);

  VecSetSizes(c, PETSC_DECIDE, Nbasis);
  VecSetFromOptions(c);
  VecSet(c, 1.0);
  matrix_expontenial(A, c, deltat, m, mmax, smax, tol, &Ac);
  PetscViewerASCIIOpen(PETSC_COMM_WORLD, "outputc.m", &cviewer);
  PetscViewerPushFormat(cviewer, PETSC_VIEWER_ASCII_MATLAB);
  VecView(Ac, cviewer);
  PetscViewerDestroy(&cviewer);
  /*
  float powtmp = pow(c1norm * nfac[m] / tol, 1. / (m + 1));
  double p = m * s;
  double s1, p1;

  int f = 0;
 */

  PetscViewerASCIIOpen(PETSC_COMM_WORLD, "outputA.m", &viewer);
  PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
  // MatView(Adeltat, viewer);
  MatView(A, viewer);
  PetscViewerDestroy(&viewer);

  PetscPrintf(PETSC_COMM_WORLD, "38\n");

  VecDestroy(&c);
  VecDestroy(&Ac);
  MatDestroy(&A);

  PetscFinalize();
  return 0;
}