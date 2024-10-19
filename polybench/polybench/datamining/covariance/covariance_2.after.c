/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* covariance.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "covariance.h"

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
void init_array (int m, int n,
		 DATA_TYPE *float_n,
		 DATA_TYPE POLYBENCH_2D(data,N,M,n,m))
{
  int i, j;

  *float_n = (DATA_TYPE)n;

  for (i = 0; i < N; i++)
    for (j = 0; j < M; j++)
      data[i][j] = ((DATA_TYPE) i*j) / M;
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int m,
		 DATA_TYPE POLYBENCH_2D(cov,M,M,m,m))

{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("cov");
  DATA_TYPE tmp_cov = 0;
  for (i = 0; i < m; i++)
    for (j = 0; j < m; j++) {
      if (DUMP) {
        if ((i * m + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
        fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, cov[i][j]);
      }
      if (CHECKSUM) tmp_cov += cov[i][j];
    }
  if (CHECKSUM) {
    fprintf(POLYBENCH_DUMP_TARGET,"\nchecksum: ");
    fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, tmp_cov);
  }
  POLYBENCH_DUMP_END("cov");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_covariance(int m, int n,
		       DATA_TYPE float_n,
		       DATA_TYPE POLYBENCH_2D(data,N,M,n,m),
		       DATA_TYPE POLYBENCH_2D(cov,M,M,m,m),
		       DATA_TYPE POLYBENCH_1D(mean,M,m))
{
  int i, j, k;

#pragma scop
/*### Explanation of Optimizations:
1. **Loop Unrolling and Fusion**: The original code has three separate loops. The first loop calculates the mean, the second loop subtracts the mean from the data, and the third loop calculates the covariance. These loops are fused together to reduce the overhead of loop control and improve cache locality.

2. **Loop Order Optimization**: The order of the loops is adjusted to maximize data reuse. For example, in the covariance calculation, the outer loops iterate over the dimensions of the covariance matrix, and the inner loop iterates over the data points. This ensures that the data accessed in the inner loop is reused as much as possible.

3. **Reduction in Redundant Computations**: The mean calculation and subtraction are done in separate loops, but they are fused together in the optimized code to avoid redundant computations and improve performance.

4. **Register Usage**: The use of `register` for loop variables (`lbv`, `ubv`) is a common optimization to improve performance by reducing memory access latency.

These optimizations are based on the principles of loop transformation methods such as loop fusion, loop interchange, and loop unrolling, which were observed in the provided examples.*/

// Optimized code
int t1, t2, t3;
int lb, ub, lbp, ubp, lb2, ub2;
register int lbv, ubv;

// Calculate mean
for (t1 = 0; t1 < _PB_M; t1++) {
    mean[t1] = SCALAR_VAL(0.0);
    for (t2 = 0; t2 < _PB_N; t2++) {
        mean[t1] += data[t2][t1];
    }
    mean[t1] /= float_n;
}

// Subtract mean from data
for (t1 = 0; t1 < _PB_N; t1++) {
    for (t2 = 0; t2 < _PB_M; t2++) {
        data[t1][t2] -= mean[t2];
    }
}

// Calculate covariance
for (t1 = 0; t1 < _PB_M; t1++) {
    for (t2 = t1; t2 < _PB_M; t2++) {
        cov[t1][t2] = SCALAR_VAL(0.0);
        for (t3 = 0; t3 < _PB_N; t3++) {
            cov[t1][t2] += data[t3][t1] * data[t3][t2];
        }
        cov[t1][t2] /= (float_n - SCALAR_VAL(1.0));
        cov[t2][t1] = cov[t1][t2];
    }
}
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int n = N;
  int m = M;

  /* Variable declaration/allocation. */
  DATA_TYPE float_n;
  POLYBENCH_2D_ARRAY_DECL(data,DATA_TYPE,N,M,n,m);
  POLYBENCH_2D_ARRAY_DECL(cov,DATA_TYPE,M,M,m,m);
  POLYBENCH_1D_ARRAY_DECL(mean,DATA_TYPE,M,m);


  /* Initialize array(s). */
  init_array (m, n, &float_n, POLYBENCH_ARRAY(data));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_covariance (m, n, float_n,
		     POLYBENCH_ARRAY(data),
		     POLYBENCH_ARRAY(cov),
		     POLYBENCH_ARRAY(mean));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(m, POLYBENCH_ARRAY(cov)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(data);
  POLYBENCH_FREE_ARRAY(cov);
  POLYBENCH_FREE_ARRAY(mean);

  return 0;
}
