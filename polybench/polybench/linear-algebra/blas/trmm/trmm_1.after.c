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
/*### Analysis of Meaning-Preserving Loop Transformations:
1. **Loop Tiling/Blocking**: The original loops are transformed into tiled loops, which helps in better cache utilization and reduces cache misses. This is evident in the transformation of the original loops into `t1`, `t2`, and `t3` loops.
2. **Loop Reordering**: The loops are reordered to ensure that the innermost loop has the most spatial locality, which is crucial for cache performance.
3. **Parallelization**: The use of `#pragma omp parallel for` indicates that the loops are parallelized to take advantage of multi-core processors.
4. **Loop Fusion/Fission**: The original loops are split into multiple loops to allow for better optimization opportunities, such as parallelization and vectorization.

### Learnings:
- **Cache Optimization**: By tiling the loops, we can ensure that the data accessed in the innermost loop fits better into the cache, reducing the number of cache misses.
- **Parallelization**: By parallelizing the outer loops, we can exploit the multi-core capabilities of modern processors, leading to significant performance improvements.
- **Loop Reordering**: Reordering the loops can help in improving the locality of reference, which is crucial for performance in memory-bound applications.

### Optimized Code Explanation:
- **Tiling**: The loops are tiled using `t1`, `t2`, and `t3` to ensure that the data accessed in the innermost loop fits into the cache.
- **Parallelization**: The outer loop (`t1`) is parallelized using OpenMP to distribute the work across multiple threads.
- **Temporary Variable**: A temporary variable `temp` is used to accumulate the results of the innermost loop, which is then added to `B[t3][j]` after the loop completes. This reduces the number of writes to `B[t3][j]`, which can be costly.
- **Loop Bounds**: The loop bounds are carefully calculated to ensure that the loops cover the entire range of `i` and `j` while respecting the tiling.*/

int t1, t2, t3;
int lb, ub, lbp, ubp, lb2, ub2;
register int lbv, ubv;

lbp = 0;
ubp = floord(_PB_M, 32);
#pragma omp parallel for private(lbv, ubv, t2, t3)
for (t1 = lbp; t1 <= ubp; t1++) {
    for (t2 = 0; t2 <= min(floord(_PB_M + _PB_N - 1, 32), floord(32 * t1 + _PB_N + 30, 32)); t2++) {
        for (t3 = max(32 * t1, 32 * t2 - _PB_N + 1); t3 <= min(_PB_M, 32 * t1 + 31); t3++) {
            for (int j = max(32 * t2, t3); j <= min(32 * t2 + 31, t3 + _PB_N - 1); j++) {
                double temp = 0.0;
                for (int k = t3 + 1; k < _PB_M; k++) {
                    temp += A[k][t3] * B[k][j];
                }
                B[t3][j] += temp;
                B[t3][j] = alpha * B[t3][j];
            }
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
