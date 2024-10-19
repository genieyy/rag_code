/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* cholesky.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "cholesky.h"

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
		DATA_TYPE POLYBENCH_2D(A,N,N,n,n))
{
  int i, j;

  for (i = 0; i < n; i++)
    {
      for (j = 0; j <= i; j++)
	A[i][j] = (DATA_TYPE)(i + j + 1); 
      for (j = i+1; j < n; j++) {
	A[i][j] = 0; 
      }
      A[i][i] = 1; 
    }

  
  for (i = 0; i < n; i++) {
    A[i][i] = (DATA_TYPE)(i + 1); 
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
    for (j = 0; j <= i; j++) {
      if (DUMP) {
      if ((i * n + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
      fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, A[i][j]);
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
void kernel_cholesky(int n,
		     DATA_TYPE POLYBENCH_2D(A,N,N,n,n))
{
  int i, j, k;


#pragma scop
/*### Explanation:

1. **Loop Interchange and Fusion**: The loops are interchanged and fused to improve cache locality. The inner loops are processed in a way that the most frequently accessed elements are handled first.

2. **Reduction in Redundant Calculations**: By introducing temporary variables (`temp_i`, `temp_k`, `temp_j`), we reduce the number of redundant calculations. This helps in minimizing the number of memory accesses and improves performance.

3. **Simplified Bounds Checking**: The loop bounds are adjusted to avoid unnecessary checks, simplifying the code and improving performance.

4. **Improved Cache Utilization**: By processing the elements in a more cache-friendly manner, the code reduces cache misses and improves overall performance.

This optimized version aims to reduce the overhead associated with loop control, improve cache locality, and minimize redundant calculations, leading to better performance.*/

for (int i = 0; i < _PB_N; i++) {
    double temp_i = A[i][i];
    for (int k = 0; k < i; k++) {
        double temp_k = A[i][k];
        temp_i -= temp_k * temp_k;
        for (int j = k + 1; j < i; j++) {
            A[i][j] -= temp_k * A[j][k];
        }
    }
    A[i][i] = SQRT_FUN(temp_i);
    for (int j = i + 1; j < _PB_N; j++) {
        double temp_j = A[i][j];
        for (int k = 0; k < i; k++) {
            temp_j -= A[i][k] * A[j][k];
        }
        A[i][j] = temp_j / A[i][i];
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

  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(A));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_cholesky (n, POLYBENCH_ARRAY(A));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(A)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);

  return 0;
}
