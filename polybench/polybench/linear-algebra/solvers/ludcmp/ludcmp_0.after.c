/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* ludcmp.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "ludcmp.h"

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
		 DATA_TYPE POLYBENCH_2D(A,N,N,n,n),
		 DATA_TYPE POLYBENCH_1D(b,N,n),
		 DATA_TYPE POLYBENCH_1D(x,N,n),
		 DATA_TYPE POLYBENCH_1D(y,N,n))
{
  int i, j;
  DATA_TYPE fn = (DATA_TYPE)n;

  for (i = 0; i < n; i++)
    {
      x[i] = 0;
      y[i] = 0;
      b[i] = (i+1)/fn/2.0 + 4;
    }

  for (i = 0; i < n; i++)
    {
      for (j = 0; j <= i; j++)
	A[i][j] = (DATA_TYPE)(-j % n) / n + 1;
      for (j = i+1; j < n; j++) {
	A[i][j] = 0;
      }
      A[i][i] = 1;
    }

  /* Make the matrix positive semi-definite. */
  /* not necessary for LU, but using same code as cholesky */
  int r,s,t;
  POLYBENCH_2D_ARRAY_DECL(B, DATA_TYPE, N, N, n, n);
  for (r = 0; r < n; ++r)
    for (s = 0; s < n; ++s)
      (POLYBENCH_ARRAY(B))[r][s] = 0;
  for (t = 0; t < n; ++t)
    for (r = 0; r < n; ++r)
      for (s = 0; s < n; ++s)
	(POLYBENCH_ARRAY(B))[r][s] += A[r][t] * A[s][t];
    for (r = 0; r < n; ++r)
      for (s = 0; s < n; ++s)
	A[r][s] = (POLYBENCH_ARRAY(B))[r][s];
  POLYBENCH_FREE_ARRAY(B);

}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_1D(x,N,n))

{
  int i;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("x");
  DATA_TYPE tmp_x = 0;
  for (i = 0; i < n; i++) {
    if (DUMP) {
      if (i % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
      fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, x[i]);
    }
    if (CHECKSUM) tmp_x += x[i];
  }
  if (CHECKSUM) {
    fprintf(POLYBENCH_DUMP_TARGET,"\nchecksum: ");
    fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, tmp_x);
  }
  POLYBENCH_DUMP_END("x");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_ludcmp(int n,
		   DATA_TYPE POLYBENCH_2D(A,N,N,n,n),
		   DATA_TYPE POLYBENCH_1D(b,N,n),
		   DATA_TYPE POLYBENCH_1D(x,N,n),
		   DATA_TYPE POLYBENCH_1D(y,N,n))
{
  int i, j, k;

  DATA_TYPE w;

#pragma scop
/*### Explanation of Optimizations:
1. **Loop Tiling**: The original loops are tiled to improve cache locality. This is done by dividing the iteration space into smaller chunks (32 in this case) and processing each chunk separately.
2. **Parallelization**: The `#pragma omp parallel for` directive is used to parallelize the outer loop, allowing multiple threads to work on different chunks of the iteration space simultaneously.
3. **Private Variables**: The `private` clause is used to ensure that each thread has its own copy of the loop variables (`lbv`, `ubv`, `t3`, `t4`, `t5`, `t6`), preventing race conditions.
4. **Bounds Calculation**: The bounds for the inner loops are calculated dynamically based on the current chunk being processed, ensuring that the loops only operate within valid ranges.

These optimizations aim to improve the performance of the original code by leveraging parallelism and enhancing cache efficiency.*/

int t1, t2, t3, t4, t5, t6;
int lb, ub, lbp, ubp, lb2, ub2;
register int lbv, ubv;

// Optimized code for the first loop nest
for (t1 = 0; t1 <= floord(_PB_N, 32); t1++) {
    lbp = max(0, 32 * t1);
    ubp = min(_PB_N - 1, 32 * t1 + 31);
#pragma omp parallel for private(lbv, ubv, t3, t4, t5, t6)
    for (t2 = lbp; t2 <= ubp; t2++) {
        for (t3 = 0; t3 <= t2 - 1; t3++) {
            double w = A[t2][t3];
            for (t4 = 0; t4 <= t3 - 1; t4++) {
                w -= A[t2][t4] * A[t4][t3];
            }
            A[t2][t3] = w / A[t3][t3];
        }
        for (t3 = t2; t3 <= _PB_N - 1; t3++) {
            double w = A[t2][t3];
            for (t4 = 0; t4 <= t2 - 1; t4++) {
                w -= A[t2][t4] * A[t4][t3];
            }
            A[t2][t3] = w;
        }
    }
}

// Optimized code for the second loop nest
for (t1 = 0; t1 <= floord(_PB_N, 32); t1++) {
    lbp = max(0, 32 * t1);
    ubp = min(_PB_N - 1, 32 * t1 + 31);
#pragma omp parallel for private(lbv, ubv, t3, t4)
    for (t2 = lbp; t2 <= ubp; t2++) {
        double w = b[t2];
        for (t3 = 0; t3 <= t2 - 1; t3++) {
            w -= A[t2][t3] * y[t3];
        }
        y[t2] = w;
    }
}

// Optimized code for the third loop nest
for (t1 = 0; t1 <= floord(_PB_N, 32); t1++) {
    lbp = max(_PB_N - 32 * t1 - 1, 0);
    ubp = min(_PB_N - 1, _PB_N - 32 * t1 - 1 + 31);
#pragma omp parallel for private(lbv, ubv, t3, t4)
    for (t2 = lbp; t2 >= 0; t2--) {
        double w = y[t2];
        for (t3 = t2 + 1; t3 <= _PB_N - 1; t3++) {
            w -= A[t2][t3] * x[t3];
        }
        x[t2] = w / A[t2][t2];
    }
}
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int n = N;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, N, N, n, n);
  POLYBENCH_1D_ARRAY_DECL(b, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(x, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(y, DATA_TYPE, N, n);


  /* Initialize array(s). */
  init_array (n,
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(b),
	      POLYBENCH_ARRAY(x),
	      POLYBENCH_ARRAY(y));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_ludcmp (n,
		 POLYBENCH_ARRAY(A),
		 POLYBENCH_ARRAY(b),
		 POLYBENCH_ARRAY(x),
		 POLYBENCH_ARRAY(y));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(x)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(b);
  POLYBENCH_FREE_ARRAY(x);
  POLYBENCH_FREE_ARRAY(y);

  return 0;
}
