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

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "jacobi-2d.h"

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

#pragma scop
/*### Analysis of Meaning-Preserving Loop Transformation Methods

1. **Loop Tiling/Blocking**:
   - The original loops are transformed into tiled loops where the iteration space is divided into smaller blocks. This allows for better cache utilization and can reduce the number of cache misses.
   - Example: The transformation from `for (i = 1; i < _PB_N - 1; i++)` to `for (t3 = max(max(1, 32 * t1 - 32 * t2), 32 * t2 + 1); t3 <= min(_PB_N - 2, 32 * t1 - 32 * t2 + 31); t3++)` introduces tiling with a block size of 32.

2. **Loop Fusion/Fission**:
   - The original loops are fused into a single loop to reduce the overhead of loop control and to allow for better parallelism.
   - Example: The two loops updating `B` and `A` are fused into a single loop in the optimized code.

3. **Parallelization**:
   - The use of `#pragma omp parallel for` indicates that the loop iterations are executed in parallel, which can significantly improve performance on multi-core processors.
   - Example: The outer loop is parallelized to exploit multiple CPU cores.

4. **Loop Interchange**:
   - The order of nested loops is changed to improve locality of reference and to better align with the memory hierarchy.
   - Example: The original nested loops are reordered to ensure that the innermost loop accesses memory in a contiguous manner.

### Learning and Application

- **Tiling**: By dividing the iteration space into smaller blocks (tiles), we can improve cache performance. This is applied by introducing new loop variables (`t1`, `t2`, `t3`, `t4`) that control the tile sizes.
- **Parallelization**: Using OpenMP directives to parallelize the outer loop can leverage multi-core processors effectively.
- **Loop Fusion**: Combining the two loops that update `B` and `A` into a single loop reduces loop overhead and improves parallelism.

The optimized code above applies these techniques to improve the performance of the original nested loops.*/

int t1, t2, t3, t4;
int lb, ub, lbp, ubp, lb2, ub2;
register int lbv, ubv;

for (t1 = 0; t1 <= floord(_PB_TSTEPS + _PB_N - 3, 32); t1++) {
    lbp = max(0, ceild(32 * t1 - _PB_N + 2, 32));
    ubp = min(floord(t1, 2), floord(_PB_N - 1, 32));
#pragma omp parallel for private(lbv, ubv, t3, t4)
    for (t2 = lbp; t2 <= ubp; t2++) {
        for (t3 = max(max(1, 32 * t1 - 32 * t2), 32 * t2 + 1); t3 <= min(_PB_N - 2, 32 * t1 - 32 * t2 + 31); t3++) {
            for (t4 = max(1, 32 * t2); t4 <= min(32 * t2 + 31, _PB_N - 2); t4++) {
                B[t3][t4] = SCALAR_VAL(0.2) * (A[t3][t4] + A[t3][t4 - 1] + A[t3][t4 + 1] + A[t3 + 1][t4] + A[t3 - 1][t4]);
                A[t3][t4] = SCALAR_VAL(0.2) * (B[t3][t4] + B[t3][t4 - 1] + B[t3][t4 + 1] + B[t3 + 1][t4] + B[t3 - 1][t4]);
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
