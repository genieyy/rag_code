/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* gesummv.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "gesummv.h"

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
void init_array(int n,
		DATA_TYPE *alpha,
		DATA_TYPE *beta,
		DATA_TYPE POLYBENCH_2D(A,N,N,n,n),
		DATA_TYPE POLYBENCH_2D(B,N,N,n,n),
		DATA_TYPE POLYBENCH_1D(x,N,n))
{
  int i, j;

  *alpha = 1.5;
  *beta = 1.2;
  for (i = 0; i < n; i++)
    {
      x[i] = (DATA_TYPE)( i % n) / n;
      for (j = 0; j < n; j++) {
	A[i][j] = (DATA_TYPE) ((i*j+1) % n) / n;
	B[i][j] = (DATA_TYPE) ((i*j+2) % n) / n;
      }
    }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_1D(y,N,n))

{
  int i;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("y");
  DATA_TYPE tmp_y = 0;
  for (i = 0; i < n; i++) {
    if (DUMP) {
      if (i % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
      fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, y[i]);
    }
    if (CHECKSUM) tmp_y += y[i];
  }
  if (CHECKSUM) {
    fprintf(POLYBENCH_DUMP_TARGET,"\nchecksum: ");
    fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, tmp_y);
  }
  POLYBENCH_DUMP_END("y");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_gesummv(int n,
		    DATA_TYPE alpha,
		    DATA_TYPE beta,
		    DATA_TYPE POLYBENCH_2D(A,N,N,n,n),
		    DATA_TYPE POLYBENCH_2D(B,N,N,n,n),
		    DATA_TYPE POLYBENCH_1D(tmp,N,n),
		    DATA_TYPE POLYBENCH_1D(x,N,n),
		    DATA_TYPE POLYBENCH_1D(y,N,n))
{
  int i, j;

#pragma scop
/*### Explanation:
1. **Loop Unrolling**: The inner loop is partially unrolled by a factor of 4 to reduce loop control overhead and improve instruction-level parallelism. This is particularly beneficial for architectures that can execute multiple instructions in parallel.
2. **Reduction in Array Accesses**: Temporary variables `tmp_val` and `y_val` are used to accumulate the results, reducing the number of array accesses and improving cache performance.
3. **Loop Fusion**: The two inner loop operations (`tmp[i]` and `y[i]`) are fused into a single loop, reducing loop control overhead.
4. **Scalar Replacement**: Temporary variables replace direct accesses to the arrays within the inner loop, improving register allocation and reducing memory bandwidth usage.
5. **OpenMP Parallelization**: The outer loop is parallelized using OpenMP to leverage multi-core processors, which can significantly improve performance for large values of `_PB_N`.*/

/**
 * Further Optimized Version:
 * 
 * ### Explanation of Transformations:
 * 1. **Loop Unrolling**: The inner loop is partially unrolled by a factor of 4 to reduce loop control overhead and improve instruction-level parallelism.
 * 2. **Reduction in Array Accesses**: By using temporary variables `tmp_val` and `y_val` to accumulate the results of the inner loop, we reduce the number of array accesses to `tmp` and `y`. This can improve performance by reducing cache misses.
 * 3. **Loop Fusion**: The two inner loop operations (`tmp[i]` and `y[i]`) are fused into a single loop, which reduces the overhead of loop control.
 * 4. **Scalar Replacement**: The temporary variables `tmp_val` and `y_val` replace direct accesses to the arrays `tmp` and `y` within the inner loop, which can help with register allocation and reduce memory bandwidth usage.
 * 5. **OpenMP Parallelization**: The outer loop is parallelized using OpenMP to leverage multi-core processors, which can significantly improve performance for large values of `_PB_N`.
 * 
 * These transformations are designed to maximize performance by reducing loop overhead, improving cache locality, and leveraging parallel processing capabilities.
 */

int t1, t2;
#pragma omp parallel for private(t1, t2)
for (t1 = 0; t1 < _PB_N; t1++) {
    double tmp_val = SCALAR_VAL(0.0);
    double y_val = SCALAR_VAL(0.0);
    for (t2 = 0; t2 < _PB_N - 3; t2 += 4) {
        tmp_val += A[t1][t2] * x[t2] + A[t1][t2 + 1] * x[t2 + 1] + A[t1][t2 + 2] * x[t2 + 2] + A[t1][t2 + 3] * x[t2 + 3];
        y_val += B[t1][t2] * x[t2] + B[t1][t2 + 1] * x[t2 + 1] + B[t1][t2 + 2] * x[t2 + 2] + B[t1][t2 + 3] * x[t2 + 3];
    }
    for (; t2 < _PB_N; t2++) {
        tmp_val += A[t1][t2] * x[t2];
        y_val += B[t1][t2] * x[t2];
    }
    y[t1] = alpha * tmp_val + beta * y_val;
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
  POLYBENCH_2D_ARRAY_DECL(B, DATA_TYPE, N, N, n, n);
  POLYBENCH_1D_ARRAY_DECL(tmp, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(x, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(y, DATA_TYPE, N, n);


  /* Initialize array(s). */
  init_array (n, &alpha, &beta,
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(B),
	      POLYBENCH_ARRAY(x));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_gesummv (n, alpha, beta,
		  POLYBENCH_ARRAY(A),
		  POLYBENCH_ARRAY(B),
		  POLYBENCH_ARRAY(tmp),
		  POLYBENCH_ARRAY(x),
		  POLYBENCH_ARRAY(y));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(y)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(B);
  POLYBENCH_FREE_ARRAY(tmp);
  POLYBENCH_FREE_ARRAY(x);
  POLYBENCH_FREE_ARRAY(y);

  return 0;
}
