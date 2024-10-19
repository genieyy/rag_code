#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <polybench.h>
#include "gramschmidt.h"
#include <omp.h>
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


/* Include polybench common header. */

/* Include benchmark-specific header. */

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
  int t1, t2, t3, t4, t5, t6, t7, t8, t9;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
lbp=0;
ubp=floord(_PB_N-2,32);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6,t7,t8,t9)
for (t2=lbp;t2<=ubp;t2++) {
  for (t4=t2;t4<=floord(_PB_N-1,32);t4++) {
    for (t5=32*t2;t5<=min(min(_PB_N-2,32*t2+31),32*t4+30);t5++) {
      lbv=max(32*t4,t5+1);
      ubv=min(_PB_N-1,32*t4+31);
#pragma ivdep
#pragma vector always
      for (t7=lbv;t7<=ubv;t7++) {
        R[t5][t7] = SCALAR_VAL(0.0);;
      }
    }
  }
}
for (t2=0;t2<=_PB_N-1;t2++) {
  nrm = SCALAR_VAL(0.0);;
  for (t4=0;t4<=_PB_M-1;t4++) {
    nrm += A[t4][t2] * A[t4][t2];;
  }
  R[t2][t2] = SQRT_FUN(nrm);;
  lbp=0;
  ubp=floord(_PB_M-1,32);
#pragma omp parallel for private(lbv,ubv,t5,t6,t7,t8,t9)
  for (t4=lbp;t4<=ubp;t4++) {
    lbv=32*t4;
    ubv=min(_PB_M-1,32*t4+31);
#pragma ivdep
#pragma vector always
    for (t5=lbv;t5<=ubv;t5++) {
      Q[t5][t2] = A[t5][t2] / R[t2][t2];;
    }
  }
  if (t2 <= _PB_N-2) {
    lbp=ceild(t2-30,32);
    ubp=floord(_PB_N-1,32);
#pragma omp parallel for private(lbv,ubv,t5,t6,t7,t8,t9)
    for (t4=lbp;t4<=ubp;t4++) {
      for (t6=0;t6<=floord(_PB_M-1,32);t6++) {
        for (t7=max(32*t4,t2+1);t7<=min(_PB_N-1,32*t4+31);t7++) {
          for (t9=32*t6;t9<=min(_PB_M-1,32*t6+31);t9++) {
            R[t2][t7] += Q[t9][t2] * A[t9][t7];;
          }
        }
      }
      for (t6=0;t6<=floord(_PB_M-1,32);t6++) {
        for (t7=max(32*t4,t2+1);t7<=min(_PB_N-1,32*t4+31);t7++) {
          lbv=32*t6;
          ubv=min(_PB_M-1,32*t6+31);
#pragma ivdep
#pragma vector always
          for (t9=lbv;t9<=ubv;t9++) {
            A[t9][t7] = A[t9][t7] - Q[t9][t2] * R[t2][t7];;
          }
        }
      }
    }
  }
}
/* End of CLooG code */

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
