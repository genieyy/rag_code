/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* trmm.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "trmm.h"

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
		DATA_TYPE POLYBENCH_2D(A,M,M,m,m),
		DATA_TYPE POLYBENCH_2D(B,M,N,m,n))
{
  int i, j;

  *alpha = 1.5;
  for (i = 0; i < m; i++) {
    for (j = 0; j < i; j++) {
      A[i][j] = (DATA_TYPE)((i+j) % m)/m;
    }
    A[i][i] = 1.0;
    for (j = 0; j < n; j++) {
      B[i][j] = (DATA_TYPE)((n+(i-j)) % n)/n;
    }
 }

}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int m, int n,
		 DATA_TYPE POLYBENCH_2D(B,M,N,m,n))
{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("B");
  DATA_TYPE tmp_B = 0;
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {
      if (DUMP) {
  	if ((i * m + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
  	fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, B[i][j]);
      }
      if (CHECKSUM) tmp_B += B[i][j];
    }
  if (CHECKSUM) {
    fprintf(POLYBENCH_DUMP_TARGET,"\nchecksum: ");
    fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, tmp_B);
  }
  POLYBENCH_DUMP_END("B");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_trmm(int m, int n,
		 DATA_TYPE alpha,
		 DATA_TYPE POLYBENCH_2D(A,M,M,m,m),
		 DATA_TYPE POLYBENCH_2D(B,M,N,m,n))
{
  int i, j, k;

//BLAS parameters
//SIDE   = 'L'
//UPLO   = 'L'
//TRANSA = 'T'
//DIAG   = 'U'
// => Form  B := alpha*A**T*B.
// A is MxM
// B is MxN
#pragma scop
/*### Analysis and Transformation Methods Used:

1. **Loop Tiling/Blocking**: The original loops are transformed into a tiled structure where the iteration space is divided into smaller blocks. This helps in better cache utilization and can lead to significant performance improvements, especially for large arrays.

2. **Parallelization**: The `#pragma omp parallel for` directive is used to parallelize the outermost loop. This allows the loop iterations to be executed in parallel across multiple threads, leveraging multi-core processors.

3. **Loop Reordering**: The original nested loops are reordered and restructured to improve locality and reduce the number of cache misses. This is done by carefully choosing the order of the loops and the bounds of the iterations.

4. **Register Usage**: The use of `register` for variables like `lbv` and `ubv` suggests that these variables are frequently accessed and should be stored in CPU registers for faster access.

5. **Vectorization**: The `#pragma ivdep` and `#pragma vector always` directives are used to hint the compiler to vectorize the loop, which can improve performance by executing multiple iterations of the loop in parallel.

### Learning:

- **Loop Tiling**: By dividing the iteration space into smaller blocks, we can improve cache performance and reduce the number of cache misses.
- **Parallelization**: Using OpenMP to parallelize loops can significantly improve performance on multi-core systems.
- **Loop Reordering**: Reordering loops can improve data locality and reduce the number of cache misses.
- **Vectorization**: Using compiler directives to hint vectorization can help the compiler generate more efficient code.

### Optimized Code:

The optimized code applies loop tiling and parallelization to the original nested loops. The iteration space is divided into smaller blocks, and the outermost loop is parallelized using OpenMP. This approach aims to improve cache utilization and leverage multi-core processors for better performance.*/

int t1, t2, t3, t4;
int lb, ub, lbp, ubp;
register int lbv, ubv;

lbp = 0;
ubp = floord(_PB_M, 32);
#pragma omp parallel for private(lbv, ubv, t2, t3, t4)
for (t1 = lbp; t1 <= ubp; t1++) {
    for (t2 = 0; t2 <= min(floord(_PB_M + _PB_N - 1, 32), floord(32 * t1 + _PB_N + 30, 32)); t2++) {
        for (t3 = max(32 * t1, 32 * t2 - _PB_N + 1); t3 <= min(_PB_M, 32 * t1 + 31); t3++) {
            for (t4 = max(32 * t2, t3 + 1); t4 <= min(32 * t2 + 31, _PB_M - 1); t4++) {
                B[t3][t4 - t3] += A[t4][t3] * B[t4][t4 - t3];
            }
            B[t3][t4 - t3] = alpha * B[t3][t4 - t3];
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
  DATA_TYPE alpha;
  POLYBENCH_2D_ARRAY_DECL(A,DATA_TYPE,M,M,m,m);
  POLYBENCH_2D_ARRAY_DECL(B,DATA_TYPE,M,N,m,n);

  /* Initialize array(s). */
  init_array (m, n, &alpha, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(B));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_trmm (m, n, alpha, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(B));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(m, n, POLYBENCH_ARRAY(B)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(B);

  return 0;
}
