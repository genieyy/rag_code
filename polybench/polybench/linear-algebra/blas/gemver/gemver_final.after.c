/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* gemver.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "gemver.h"

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
void init_array (int n,
		 DATA_TYPE *alpha,
		 DATA_TYPE *beta,
		 DATA_TYPE POLYBENCH_2D(A,N,N,n,n),
		 DATA_TYPE POLYBENCH_1D(u1,N,n),
		 DATA_TYPE POLYBENCH_1D(v1,N,n),
		 DATA_TYPE POLYBENCH_1D(u2,N,n),
		 DATA_TYPE POLYBENCH_1D(v2,N,n),
		 DATA_TYPE POLYBENCH_1D(w,N,n),
		 DATA_TYPE POLYBENCH_1D(x,N,n),
		 DATA_TYPE POLYBENCH_1D(y,N,n),
		 DATA_TYPE POLYBENCH_1D(z,N,n))
{
  int i, j;

  *alpha = 1.5;
  *beta = 1.2;

  DATA_TYPE fn = (DATA_TYPE)n;

  for (i = 0; i < n; i++)
    {
      u1[i] = i;
      u2[i] = ((i+1)/fn)/2.0;
      v1[i] = ((i+1)/fn)/4.0;
      v2[i] = ((i+1)/fn)/6.0;
      y[i] = ((i+1)/fn)/8.0;
      z[i] = ((i+1)/fn)/9.0;
      x[i] = 0.0;
      w[i] = 0.0;
      for (j = 0; j < n; j++)
        A[i][j] = (DATA_TYPE) (i*j % n) / n;
    }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_1D(w,N,n))
{
  int i;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("w");
  DATA_TYPE tmp_w = 0;
  for (i = 0; i < n; i++) {
    if (DUMP) {
      if (i % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
      fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, w[i]);
    }
    if (CHECKSUM) tmp_w += w[i];
  }
  if (CHECKSUM) {
    fprintf(POLYBENCH_DUMP_TARGET,"\nchecksum: ");
    fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, tmp_w);
  }
  POLYBENCH_DUMP_END("w");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_gemver(int n,
		   DATA_TYPE alpha,
		   DATA_TYPE beta,
		   DATA_TYPE POLYBENCH_2D(A,N,N,n,n),
		   DATA_TYPE POLYBENCH_1D(u1,N,n),
		   DATA_TYPE POLYBENCH_1D(v1,N,n),
		   DATA_TYPE POLYBENCH_1D(u2,N,n),
		   DATA_TYPE POLYBENCH_1D(v2,N,n),
		   DATA_TYPE POLYBENCH_1D(w,N,n),
		   DATA_TYPE POLYBENCH_1D(x,N,n),
		   DATA_TYPE POLYBENCH_1D(y,N,n),
		   DATA_TYPE POLYBENCH_1D(z,N,n))
{
  int i, j;

#pragma scop
/*### Explanation:
1. **Loop Unrolling and Tiling**: The loops are already tiled with a tile size of 32, which helps in better cache utilization and parallelization.
2. **Reduction in Memory Accesses**: By storing frequently accessed variables (`u1[t2]`, `u2[t2]`, `x[t2]`, `w[t2]`) in temporary variables (`temp_u1`, `temp_u2`, `temp_x`, `temp_w`), we reduce the number of memory accesses, which can significantly improve performance.
3. **OpenMP Parallelization**: The loops are parallelized using OpenMP, which allows for multi-threaded execution. The `private(t2)` clause ensures that each thread has its own copy of `t2`, preventing race conditions.
4. **Avoiding Redundant Computations**: By using temporary variables, we avoid redundant computations and memory accesses, which can be costly, especially for large values of `_PB_N`.

This optimized version should provide better performance by reducing memory latency and improving parallel execution efficiency.*/

#include <omp.h>
#include <math.h>

// Helper macro to compute floord
#define floord(n, d) ((n) / (d) - ((n) % (d) < 0))

// Optimized Version:
int t1, t2;
register int lbv, ubv;

// First loop transformation
for (int t1 = 0; t1 <= floord(_PB_N - 1, 32); t1++) {
    lbv = max(0, 32 * t1);
    ubv = min(_PB_N - 1, 32 * t1 + 31);
#pragma omp parallel for private(t2)
    for (int t2 = lbv; t2 <= ubv; t2++) {
        double temp_u1 = u1[t2];
        double temp_u2 = u2[t2];
        for (int j = 0; j < _PB_N; j++) {
            A[t2][j] = A[t2][j] + temp_u1 * v1[j] + temp_u2 * v2[j];
        }
    }
}

// Second loop transformation
for (int t1 = 0; t1 <= floord(_PB_N - 1, 32); t1++) {
    lbv = max(0, 32 * t1);
    ubv = min(_PB_N - 1, 32 * t1 + 31);
#pragma omp parallel for private(t2)
    for (int t2 = lbv; t2 <= ubv; t2++) {
        double temp_x = x[t2];
        for (int j = 0; j < _PB_N; j++) {
            temp_x += beta * A[j][t2] * y[j];
        }
        x[t2] = temp_x;
    }
}

// Third loop transformation
for (int t1 = 0; t1 <= floord(_PB_N - 1, 32); t1++) {
    lbv = max(0, 32 * t1);
    ubv = min(_PB_N - 1, 32 * t1 + 31);
#pragma omp parallel for
    for (int t2 = lbv; t2 <= ubv; t2++) {
        x[t2] += z[t2];
    }
}

// Fourth loop transformation
for (int t1 = 0; t1 <= floord(_PB_N - 1, 32); t1++) {
    lbv = max(0, 32 * t1);
    ubv = min(_PB_N - 1, 32 * t1 + 31);
#pragma omp parallel for private(t2)
    for (int t2 = lbv; t2 <= ubv; t2++) {
        double temp_w = w[t2];
        for (int j = 0; j < _PB_N; j++) {
            temp_w += alpha * A[t2][j] * x[j];
        }
        w[t2] = temp_w;
    }
}
#pragma endscop
}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int n = N;

  /* Variable declaration/allocation. */
  DATA_TYPE alpha;
  DATA_TYPE beta;
  POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, N, N, n, n);
  POLYBENCH_1D_ARRAY_DECL(u1, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(v1, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(u2, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(v2, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(w, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(x, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(y, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(z, DATA_TYPE, N, n);


  /* Initialize array(s). */
  init_array (n, &alpha, &beta,
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(u1),
	      POLYBENCH_ARRAY(v1),
	      POLYBENCH_ARRAY(u2),
	      POLYBENCH_ARRAY(v2),
	      POLYBENCH_ARRAY(w),
	      POLYBENCH_ARRAY(x),
	      POLYBENCH_ARRAY(y),
	      POLYBENCH_ARRAY(z));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_gemver (n, alpha, beta,
		 POLYBENCH_ARRAY(A),
		 POLYBENCH_ARRAY(u1),
		 POLYBENCH_ARRAY(v1),
		 POLYBENCH_ARRAY(u2),
		 POLYBENCH_ARRAY(v2),
		 POLYBENCH_ARRAY(w),
		 POLYBENCH_ARRAY(x),
		 POLYBENCH_ARRAY(y),
		 POLYBENCH_ARRAY(z));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(w)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(u1);
  POLYBENCH_FREE_ARRAY(v1);
  POLYBENCH_FREE_ARRAY(u2);
  POLYBENCH_FREE_ARRAY(v2);
  POLYBENCH_FREE_ARRAY(w);
  POLYBENCH_FREE_ARRAY(x);
  POLYBENCH_FREE_ARRAY(y);
  POLYBENCH_FREE_ARRAY(z);

  return 0;
}
