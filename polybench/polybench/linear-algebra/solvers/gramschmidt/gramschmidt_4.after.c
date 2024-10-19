/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* gramschmidt.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "gramschmidt.h"

#include <omp.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y) ((x) > (y)? (x) : (y))
#define min(x,y) ((x) < (y)? (x) : (y))
#ifndef POLYBENCH_DUMP_ARRAYS
#define DUMP 0
#else
#define DUMP 1
#endif
#ifndef POLYBENCH_CHECKSUM_ARRAYS
#define CHECKSUM 0
#else
#define CHECKSUM 1
#endif


/* Array initialization. */
static
void init_array(int m, int n,
		DATA_TYPE POLYBENCH_2D(A,M,N,m,n),
		DATA_TYPE POLYBENCH_2D(R,N,N,n,n),
		DATA_TYPE POLYBENCH_2D(Q,M,N,m,n))
{
  int i, j;

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {
      A[i][j] = (((DATA_TYPE) ((i*j) % m) / m )*100) + 10;
      Q[i][j] = 0.0;
    }
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      R[i][j] = 0.0;
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int m, int n,
		 DATA_TYPE POLYBENCH_2D(A,M,N,m,n),
		 DATA_TYPE POLYBENCH_2D(R,N,N,n,n),
		 DATA_TYPE POLYBENCH_2D(Q,M,N,m,n))
{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("R");
  DATA_TYPE tmp_R = 0;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      if (DUMP) {
  	if ((i*n+j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
  	fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, R[i][j]);
      }
      if (CHECKSUM) tmp_R += R[i][j];
    }
  if (CHECKSUM) {
    fprintf(POLYBENCH_DUMP_TARGET,"\nchecksum: ");
    fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, tmp_R);
  }
  POLYBENCH_DUMP_END("R");

  POLYBENCH_DUMP_BEGIN("Q");
  DATA_TYPE tmp_Q = 0;
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {
      if (DUMP) {
  	if ((i*n+j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
  	fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, Q[i][j]);
      }
      if (CHECKSUM) tmp_Q += Q[i][j];
    }
  if (CHECKSUM) {
    fprintf(POLYBENCH_DUMP_TARGET,"\nchecksum: ");
    fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, tmp_Q);
  }
  POLYBENCH_DUMP_END("Q");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
/* QR Decomposition with Modified Gram Schmidt:
 http://www.inf.ethz.ch/personal/gander/ */
static
void kernel_gramschmidt(int m, int n,
			DATA_TYPE POLYBENCH_2D(A,M,N,m,n),
			DATA_TYPE POLYBENCH_2D(R,N,N,n,n),
			DATA_TYPE POLYBENCH_2D(Q,M,N,m,n))
{
  int i, j, k;

  DATA_TYPE nrm;

#pragma scop
for (k = 0; k < _PB_N; k++) {
    double nrm = SCALAR_VAL(0.0);
    for (i = 0; i < _PB_M; i += 4) {
        double A_ik1 = A[i][k];
        double A_ik2 = A[i + 1][k];
        double A_ik3 = A[i + 2][k];
        double A_ik4 = A[i + 3][k];
        nrm += A_ik1 * A_ik1 + A_ik2 * A_ik2 + A_ik3 * A_ik3 + A_ik4 * A_ik4;
        Q[i][k] = A_ik1;
        Q[i + 1][k] = A_ik2;
        Q[i + 2][k] = A_ik3;
        Q[i + 3][k] = A_ik4;
    }
    if (_PB_M % 4 != 0) {
        for (i = _PB_M - (_PB_M % 4); i < _PB_M; i++) {
            double A_ik = A[i][k];
            nrm += A_ik * A_ik;
            Q[i][k] = A_ik;
        }
    }
    double R_kk = SQRT_FUN(nrm);
    R[k][k] = R_kk;
    for (i = 0; i < _PB_M; i += 4) {
        Q[i][k] /= R_kk;
        Q[i + 1][k] /= R_kk;
        Q[i + 2][k] /= R_kk;
        Q[i + 3][k] /= R_kk;
    }
    if (_PB_M % 4 != 0) {
        for (i = _PB_M - (_PB_M % 4); i < _PB_M; i++) {
            Q[i][k] /= R_kk;
        }
    }
    for (j = k + 1; j < _PB_N; j++) {
        double R_kj = SCALAR_VAL(0.0);
        for (i = 0; i < _PB_M; i += 4) {
            R_kj += Q[i][k] * A[i][j] + Q[i + 1][k] * A[i + 1][j] + Q[i + 2][k] * A[i + 2][j] + Q[i + 3][k] * A[i + 3][j];
        }
        if (_PB_M % 4 != 0) {
            for (i = _PB_M - (_PB_M % 4); i < _PB_M; i++) {
                R_kj += Q[i][k] * A[i][j];
            }
        }
        R[k][j] = R_kj;
        for (i = 0; i < _PB_M; i += 4) {
            A[i][j] -= Q[i][k] * R_kj;
            A[i + 1][j] -= Q[i + 1][k] * R_kj;
            A[i + 2][j] -= Q[i + 2][k] * R_kj;
            A[i + 3][j] -= Q[i + 3][k] * R_kj;
        }
        if (_PB_M % 4 != 0) {
            for (i = _PB_M - (_PB_M % 4); i < _PB_M; i++) {
                A[i][j] -= Q[i][k] * R_kj;
            }
        }
    }
}
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int m = M;
  int n = N;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(A,DATA_TYPE,M,N,m,n);
  POLYBENCH_2D_ARRAY_DECL(R,DATA_TYPE,N,N,n,n);
  POLYBENCH_2D_ARRAY_DECL(Q,DATA_TYPE,M,N,m,n);

  /* Initialize array(s). */
  init_array (m, n,
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(R),
	      POLYBENCH_ARRAY(Q));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_gramschmidt (m, n,
		      POLYBENCH_ARRAY(A),
		      POLYBENCH_ARRAY(R),
		      POLYBENCH_ARRAY(Q));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(m, n, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(R), POLYBENCH_ARRAY(Q)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(R);
  POLYBENCH_FREE_ARRAY(Q);

  return 0;
}
