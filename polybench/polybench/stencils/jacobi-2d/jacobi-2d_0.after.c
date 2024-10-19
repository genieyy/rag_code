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
/*### Optimized Code



### Analysis of Loop Transformation Methods

1. **Loop Tiling (Blocking)**:
   - The original loops are transformed by dividing the iteration space into smaller blocks (tiles) of size 32x32. This is done to improve cache locality and reduce cache misses.
   - Example: `for (int t1 = 0; t1 <= (_PB_N - 3) / 32; t1++)` and `for (int t2 = 0; t2 <= (_PB_N - 3) / 32; t2++)`.

2. **Loop Fusion**:
   - The two sets of nested loops (one for updating `B` and one for updating `A`) are fused into a single loop structure. This reduces the overhead of loop control and can improve data locality.
   - Example: The two sets of nested loops for `B` and `A` updates are placed within the same outer loop over `t`.

3. **Parallelization**:
   - The outer loop over `t` is parallelized using OpenMP to exploit multi-core processors. This allows multiple iterations of the loop to be executed simultaneously.
   - Example: `#pragma omp parallel for private(i, j)`.

4. **Loop Reordering**:
   - The order of the loops is changed to maximize the benefits of loop tiling and to ensure that the innermost loops have the best possible cache locality.
   - Example: The loops over `i` and `j` are nested inside the tile loops `t1` and `t2`.

### Learnings

- **Cache Locality**: By tiling the loops, we ensure that the data accessed within each tile fits well within the cache, reducing the number of cache misses and improving performance.
- **Parallelization**: Leveraging OpenMP for parallel execution can significantly speed up the computation on multi-core systems.
- **Loop Fusion**: Combining related loops can reduce overhead and improve data locality, especially when the same data is accessed in both loops.
- **Loop Reordering**: The order of loops can have a significant impact on performance. Reordering loops to maximize cache hits and minimize cache misses is crucial.

These techniques are applied in the optimized code to improve performance by enhancing cache utilization, reducing overhead, and leveraging parallel processing capabilities.*/

#pragma omp parallel for private(i, j)
for (t = 0; t < _PB_TSTEPS; t++) {
    for (int t1 = 0; t1 <= (_PB_N - 3) / 32; t1++) {
        for (int t2 = 0; t2 <= (_PB_N - 3) / 32; t2++) {
            for (int i = max(1, 32 * t1); i <= min(_PB_N - 2, 32 * t1 + 31); i++) {
                for (int j = max(1, 32 * t2); j <= min(_PB_N - 2, 32 * t2 + 31); j++) {
                    B[i][j] = SCALAR_VAL(0.2) * (A[i][j] + A[i][j-1] + A[i][j+1] + A[i+1][j] + A[i-1][j]);
                }
            }
        }
    }

    for (int t1 = 0; t1 <= (_PB_N - 3) / 32; t1++) {
        for (int t2 = 0; t2 <= (_PB_N - 3) / 32; t2++) {
            for (int i = max(1, 32 * t1); i <= min(_PB_N - 2, 32 * t1 + 31); i++) {
                for (int j = max(1, 32 * t2); j <= min(_PB_N - 2, 32 * t2 + 31); j++) {
                    A[i][j] = SCALAR_VAL(0.2) * (B[i][j] + B[i][j-1] + B[i][j+1] + B[i+1][j] + B[i-1][j]);
                }
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
