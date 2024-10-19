#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <polybench.h>
#include "jacobi-2d.h"
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
/* jacobi-2d.c: this file is part of PolyBench/C */


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
void init_array (int n,
		 DATA_TYPE POLYBENCH_2D(A,N,N,n,n),
		 DATA_TYPE POLYBENCH_2D(B,N,N,n,n))
{
  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      {
	A[i][j] = ((DATA_TYPE) i*(j+2) + 2) / n;
	B[i][j] = ((DATA_TYPE) i*(j+3) + 3) / n;
      }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_2D(A,N,N,n,n))

{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("A");
  DATA_TYPE tmp_A = 0;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      if (DUMP) {
        if ((i * n + j) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
        fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, A[i][j]);
      }
      if (CHECKSUM) tmp_A += A[i][j];
    }
  if (CHECKSUM) {
    fprintf(POLYBENCH_DUMP_TARGET,"\nchecksum: ");
    fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, tmp_A);
  }
  POLYBENCH_DUMP_END("A");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_jacobi_2d(int tsteps,
			    int n,
			    DATA_TYPE POLYBENCH_2D(A,N,N,n,n),
			    DATA_TYPE POLYBENCH_2D(B,N,N,n,n))
{
  int t, i, j;
  int t1, t2, t3, t4, t5, t6;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
for (t1=0;t1<=floord(3*_PB_TSTEPS+_PB_N-4,32);t1++) {
  lbp=max(ceild(2*t1,3),ceild(32*t1-_PB_TSTEPS+1,32));
  ubp=min(min(floord(2*_PB_TSTEPS+_PB_N-3,32),floord(64*t1+_PB_N+61,96)),t1);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6)
  for (t2=lbp;t2<=ubp;t2++) {
    for (t3=max(ceild(32*t2-_PB_N-28,32),2*t1-2*t2);t3<=min(min(floord(2*_PB_TSTEPS+_PB_N-3,32),floord(32*t2+_PB_N+28,32)),floord(64*t1-64*t2+_PB_N+61,32));t3++) {
      if ((_PB_N == 1301) && (t1 <= floord(2*t2+t3-41,2)) && (t2 <= t3-1)) {
        for (t5=max(32*t2,32*t3-1298);t5<=32*t2+31;t5++) {
          A[(-32*t3+t5+1299)][1299] = SCALAR_VAL(0.2) * (B[(-32*t3+t5+1299)][1299] + B[(-32*t3+t5+1299)][1299 -1] + B[(-32*t3+t5+1299)][1+1299] + B[1+(-32*t3+t5+1299)][1299] + B[(-32*t3+t5+1299)-1][1299]);;
        }
      }
      if ((_PB_N == 1301) && (t1 <= floord(3*t2-41,2)) && (t2 >= t3)) {
        for (t6=max(32*t3,32*t2-1298);t6<=min(32*t2,32*t3+31);t6++) {
          A[1299][(-32*t2+t6+1299)] = SCALAR_VAL(0.2) * (B[1299][(-32*t2+t6+1299)] + B[1299][(-32*t2+t6+1299)-1] + B[1299][1+(-32*t2+t6+1299)] + B[1+1299][(-32*t2+t6+1299)] + B[1299 -1][(-32*t2+t6+1299)]);;
        }
      }
      for (t4=max(max(ceild(32*t2-_PB_N+2,2),ceild(32*t3-_PB_N+2,2)),32*t1-32*t2);t4<=min(min(min(_PB_TSTEPS-1,16*t2+14),16*t3+14),32*t1-32*t2+31);t4++) {
        if (t2 <= floord(t4,16)) {
          for (t6=max(32*t3,2*t4+1);t6<=min(32*t3+31,2*t4+_PB_N-2);t6++) {
            B[1][(-2*t4+t6)] = SCALAR_VAL(0.2) * (A[1][(-2*t4+t6)] + A[1][(-2*t4+t6)-1] + A[1][1+(-2*t4+t6)] + A[1+1][(-2*t4+t6)] + A[1 -1][(-2*t4+t6)]);;
          }
        }
        for (t5=max(32*t2,2*t4+2);t5<=min(32*t2+31,2*t4+_PB_N-2);t5++) {
          if (t3 <= floord(t4,16)) {
            B[(-2*t4+t5)][1] = SCALAR_VAL(0.2) * (A[(-2*t4+t5)][1] + A[(-2*t4+t5)][1 -1] + A[(-2*t4+t5)][1+1] + A[1+(-2*t4+t5)][1] + A[(-2*t4+t5)-1][1]);;
          }
          for (t6=max(32*t3,2*t4+2);t6<=min(32*t3+31,2*t4+_PB_N-2);t6++) {
            B[(-2*t4+t5)][(-2*t4+t6)] = SCALAR_VAL(0.2) * (A[(-2*t4+t5)][(-2*t4+t6)] + A[(-2*t4+t5)][(-2*t4+t6)-1] + A[(-2*t4+t5)][1+(-2*t4+t6)] + A[1+(-2*t4+t5)][(-2*t4+t6)] + A[(-2*t4+t5)-1][(-2*t4+t6)]);;
            A[(-2*t4+t5-1)][(-2*t4+t6-1)] = SCALAR_VAL(0.2) * (B[(-2*t4+t5-1)][(-2*t4+t6-1)] + B[(-2*t4+t5-1)][(-2*t4+t6-1)-1] + B[(-2*t4+t5-1)][1+(-2*t4+t6-1)] + B[1+(-2*t4+t5-1)][(-2*t4+t6-1)] + B[(-2*t4+t5-1)-1][(-2*t4+t6-1)]);;
          }
          if (t3 >= ceild(2*t4+_PB_N-32,32)) {
            A[(-2*t4+t5-1)][(_PB_N-2)] = SCALAR_VAL(0.2) * (B[(-2*t4+t5-1)][(_PB_N-2)] + B[(-2*t4+t5-1)][(_PB_N-2)-1] + B[(-2*t4+t5-1)][1+(_PB_N-2)] + B[1+(-2*t4+t5-1)][(_PB_N-2)] + B[(-2*t4+t5-1)-1][(_PB_N-2)]);;
          }
        }
        if (t2 >= ceild(2*t4+_PB_N-32,32)) {
          for (t6=max(32*t3,2*t4+2);t6<=min(32*t3+31,2*t4+_PB_N-1);t6++) {
            A[(_PB_N-2)][(-2*t4+t6-1)] = SCALAR_VAL(0.2) * (B[(_PB_N-2)][(-2*t4+t6-1)] + B[(_PB_N-2)][(-2*t4+t6-1)-1] + B[(_PB_N-2)][1+(-2*t4+t6-1)] + B[1+(_PB_N-2)][(-2*t4+t6-1)] + B[(_PB_N-2)-1][(-2*t4+t6-1)]);;
          }
        }
      }
      if ((t1 >= ceild(3*t2-1,2)) && (t2 <= min(floord(_PB_TSTEPS-16,16),t3-1))) {
        for (t6=32*t3;t6<=min(32*t3+31,32*t2+_PB_N+28);t6++) {
          B[1][(-32*t2+t6-30)] = SCALAR_VAL(0.2) * (A[1][(-32*t2+t6-30)] + A[1][(-32*t2+t6-30)-1] + A[1][1+(-32*t2+t6-30)] + A[1+1][(-32*t2+t6-30)] + A[1 -1][(-32*t2+t6-30)]);;
        }
      }
      if ((t1 >= ceild(2*t2+t3-1,2)) && (t2 >= t3) && (t3 <= floord(_PB_TSTEPS-16,16))) {
        for (t5=max(32*t2,32*t3+31);t5<=min(32*t2+31,32*t3+_PB_N+28);t5++) {
          B[(-32*t3+t5-30)][1] = SCALAR_VAL(0.2) * (A[(-32*t3+t5-30)][1] + A[(-32*t3+t5-30)][1 -1] + A[(-32*t3+t5-30)][1+1] + A[1+(-32*t3+t5-30)][1] + A[(-32*t3+t5-30)-1][1]);;
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
  int n = N;
  int tsteps = TSTEPS;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(B, DATA_TYPE, N, N, n, n);


  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(B));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_jacobi_2d(tsteps, n, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(B));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(A)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(B);

  return 0;
}
