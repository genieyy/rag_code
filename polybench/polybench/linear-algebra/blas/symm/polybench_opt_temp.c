/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* symm.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "symm.h"

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
		DATA_TYPE *alpha,
		DATA_TYPE *beta,
		DATA_TYPE POLYBENCH_2D(C,M,N,m,n),
		DATA_TYPE POLYBENCH_2D(A,M,M,m,m),
		DATA_TYPE POLYBENCH_2D(B,M,N,m,n))
{
  int i, j;

  *alpha = 3.0;
  *beta = 0.8;
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {
      C[i][j] = (DATA_TYPE) ((i * i + j * j) % 100) / (m * m + n * n);
      B[i][j] = (DATA_TYPE) ((i * j + 1) % 100) / (m * n + 1);
    }
  for (i = 0; i < m; i++) {
    for (j = 0; j <= i; j++)
      A[i][j] = (DATA_TYPE) ((i * i - j * j) % 100) / (m * m);
    for (j = i+1; j < m; j++)
      A[i][j] = -2; 
  }
}



/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int m, int n,
		 DATA_TYPE POLYBENCH_2D(C,M,N,m,n))
{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("C");
  DATA_TYPE tmp_C = 0;
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {
      if (DUMP) {
  	if ((i * m + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
  	fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, C[i][j]);
      }
      if (CHECKSUM) tmp_C += C[i][j];
    }
  if (CHECKSUM) {
    fprintf(POLYBENCH_DUMP_TARGET,"\nchecksum: ");
    fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, tmp_C);
  }
  POLYBENCH_DUMP_END("C");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_symm(int m, int n,
		 DATA_TYPE alpha,
		 DATA_TYPE beta,
		 DATA_TYPE POLYBENCH_2D(C,M,N,m,n),
		 DATA_TYPE POLYBENCH_2D(A,M,M,m,m),
		 DATA_TYPE POLYBENCH_2D(B,M,N,m,n))
{
  int i, j, k;
  DATA_TYPE temp2;

//BLAS PARAMS
//SIDE = 'L'
//UPLO = 'L'
// =>  Form  C := alpha*A*B + beta*C
// A is MxM
// B is MxN
// C is MxN
//note that due to Fortran array layout, the code below more closely resembles upper triangular case in BLAS
#pragma scop
/*### Explanation of the Transformation:
1. **Loop Unrolling**: The inner loop over `k` is partially unrolled to reduce the overhead of loop control and improve instruction-level parallelism.
2. **Reduction in Redundant Calculations**: By unrolling the loop, we reduce the number of times the loop control variable `k` is checked and incremented, which can lead to performance improvements.
3. **Temporary Variable Optimization**: The variable `temp2` is still used to accumulate the sum within the inner loop, but the unrolling helps in reducing the number of iterations, thus improving performance.

These transformations aim to improve the performance by reducing the overhead of loop control and maintaining the original logic of the code.*/

/*### Explanation of the Transformation:
1. **Loop Unrolling**: The inner loop over `k` is partially unrolled to reduce the overhead of loop control and improve instruction-level parallelism.
2. **Reduction in Redundant Calculations**: By unrolling the loop, we reduce the number of times the loop control variable `k` is checked and incremented, which can lead to performance improvements.
3. **Temporary Variable Optimization**: The variable `temp2` is still used to accumulate the sum within the inner loop, but the unrolling helps in reducing the number of iterations, thus improving performance.

These transformations aim to improve the performance by reducing the overhead of loop control and maintaining the original logic of the code.*/

for (int i = 0; i < _PB_M; i++) {
    for (int j = 0; j < _PB_N; j++) {
        double temp2 = 0;
        int k;
        for (k = 0; k + 4 <= i; k += 4) {
            C[k][j] += alpha * B[i][j] * A[i][k];
            temp2 += B[k][j] * A[i][k];
            C[k + 1][j] += alpha * B[i][j] * A[i][k + 1];
            temp2 += B[k + 1][j] * A[i][k + 1];
            C[k + 2][j] += alpha * B[i][j] * A[i][k + 2];
            temp2 += B[k + 2][j] * A[i][k + 2];
            C[k + 3][j] += alpha * B[i][j] * A[i][k + 3];
            temp2 += B[k + 3][j] * A[i][k + 3];
        }
        for (; k < i; k++) {
            C[k][j] += alpha * B[i][j] * A[i][k];
            temp2 += B[k][j] * A[i][k];
        }
        C[i][j] = beta * C[i][j] + alpha * B[i][j] * A[i][i] + alpha * temp2;
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
  DATA_TYPE alpha;
  DATA_TYPE beta;
  POLYBENCH_2D_ARRAY_DECL(C,DATA_TYPE,M,N,m,n);
  POLYBENCH_2D_ARRAY_DECL(A,DATA_TYPE,M,M,m,m);
  POLYBENCH_2D_ARRAY_DECL(B,DATA_TYPE,M,N,m,n);

  /* Initialize array(s). */
  init_array (m, n, &alpha, &beta,
	      POLYBENCH_ARRAY(C),
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(B));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_symm (m, n,
	       alpha, beta,
	       POLYBENCH_ARRAY(C),
	       POLYBENCH_ARRAY(A),
	       POLYBENCH_ARRAY(B));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(m, n, POLYBENCH_ARRAY(C)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(C);
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(B);

  return 0;
}
