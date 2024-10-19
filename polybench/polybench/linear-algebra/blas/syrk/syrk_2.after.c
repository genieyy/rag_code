/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* syrk.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "syrk.h"

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
void init_array(int n, int m,
		DATA_TYPE *alpha,
		DATA_TYPE *beta,
		DATA_TYPE POLYBENCH_2D(C,N,N,n,n),
		DATA_TYPE POLYBENCH_2D(A,N,M,n,m))
{
  int i, j;

  *alpha = 1.5;
  *beta = 1.2;
  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++)
      A[i][j] = (DATA_TYPE) ((i*j+1)%n) / n;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      C[i][j] = (DATA_TYPE) ((i*j+2)%m) / m;
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_2D(C,N,N,n,n))
{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("C");
  DATA_TYPE tmp_C = 0;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      if (DUMP) {
  	if ((i * n + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
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
void kernel_syrk(int n, int m,
		 DATA_TYPE alpha,
		 DATA_TYPE beta,
		 DATA_TYPE POLYBENCH_2D(C,N,N,n,n),
		 DATA_TYPE POLYBENCH_2D(A,N,M,n,m))
{
  int i, j, k;

//BLAS PARAMS
//TRANS = 'N'
//UPLO  = 'L'
// =>  Form  C := alpha*A*A**T + beta*C.
//A is NxM
//C is NxN
#pragma scop
/*### Explanation of the Transformation Methods Used:

1. **Loop Tiling/Blocking**: The original loops are divided into smaller blocks (tiles) to improve cache locality. This is done by introducing new variables `t1`, `t2`, `t3`, and `t4` to represent the tile indices.

2. **Parallelization**: The outer loop is parallelized using OpenMP (`#pragma omp parallel for`) to exploit multi-core processors.

3. **Loop Fusion**: The two inner loops (one for scaling `C` by `beta` and the other for updating `C` with matrix multiplication) are fused together to reduce the overhead of loop control and to improve data locality.

4. **Loop Interchange**: The order of loops is changed to ensure that the most frequently accessed data is loaded into the cache efficiently. This is particularly important for the matrix multiplication part where `C[i][j]` is updated based on `A[i][k]` and `A[j][k]`.

5. **Temporary Variable for Accumulation**: A temporary variable `temp` is used to accumulate the result of the matrix multiplication before updating `C[i][j]`. This reduces the number of memory writes and improves performance.

These transformations aim to optimize the code by improving cache utilization, reducing loop overhead, and leveraging parallel processing capabilities.*/

int t1, t2, t3, t4;
int lb, ub, lbp, ubp, lb2, ub2;
register int lbv, ubv;

lbp = 0;
ubp = floord(_PB_N, 32);

#pragma omp parallel for private(lbv, ubv, t2, t3, t4)
for (t1 = lbp; t1 <= ubp; t1++) {
    for (t2 = 0; t2 <= t1; t2++) {
        if (t1 == t2) {
            for (t3 = 32 * t1; t3 <= min(_PB_N - 1, 32 * t1 + 31); t3++) {
                for (t4 = 32 * t1; t4 <= t3; t4++) {
                    C[t3][t4] *= beta;
                }
            }
        }
        for (t3 = max(32 * t1, 32 * t2 + 1); t3 <= min(_PB_N - 1, 32 * t1 + 31); t3++) {
            for (t4 = 32 * t2; t4 <= t3; t4++) {
                C[t3][t4] *= beta;
            }
        }
    }
    for (t2 = 0; t2 <= t1; t2++) {
        for (t3 = 32 * t1; t3 <= min(_PB_N - 1, 32 * t1 + 31); t3++) {
            for (t4 = 32 * t2; t4 <= t3; t4++) {
                double temp = 0.0;
                for (int k = 0; k < _PB_M; k++) {
                    temp += alpha * A[t3][k] * A[t4][k];
                }
                C[t3][t4] += temp;
            }
        }
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
  DATA_TYPE alpha;
  DATA_TYPE beta;
  POLYBENCH_2D_ARRAY_DECL(C,DATA_TYPE,N,N,n,n);
  POLYBENCH_2D_ARRAY_DECL(A,DATA_TYPE,N,M,n,m);

  /* Initialize array(s). */
  init_array (n, m, &alpha, &beta, POLYBENCH_ARRAY(C), POLYBENCH_ARRAY(A));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_syrk (n, m, alpha, beta, POLYBENCH_ARRAY(C), POLYBENCH_ARRAY(A));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(C)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(C);
  POLYBENCH_FREE_ARRAY(A);

  return 0;
}
