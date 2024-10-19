/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* floyd-warshall.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "floyd-warshall.h"

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
		 DATA_TYPE POLYBENCH_2D(path,N,N,n,n))
{
  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      path[i][j] = i*j%7+1;
      if ((i+j)%13 == 0 || (i+j)%7==0 || (i+j)%11 == 0)
         path[i][j] = 999;
    }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_2D(path,N,N,n,n))

{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("path");
  DATA_TYPE tmp_path = 0;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      if (DUMP) {
        if ((i * n + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
        fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, path[i][j]);
      }
      if (CHECKSUM) tmp_path += path[i][j];
    }
  if (CHECKSUM) {
    fprintf(POLYBENCH_DUMP_TARGET,"\nchecksum: ");
    fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, tmp_path);
  }
  POLYBENCH_DUMP_END("path");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_floyd_warshall(int n,
			   DATA_TYPE POLYBENCH_2D(path,N,N,n,n))
{
  int i, j, k;

#pragma scop
/*### Analysis and Learnings:

1. **Loop Tiling (Blocking):**
   - The original code uses a nested loop structure with three levels of loops (i, j, k). The optimized code applies loop tiling to these loops, which is a technique to improve cache utilization by breaking the problem into smaller blocks. This is evident from the introduction of variables like `t1`, `t2`, `t3`, and `t4` which represent the tile indices.

2. **Parallelization:**
   - The `#pragma omp parallel for` directive is used to parallelize the outermost loop. This allows the loop iterations to be executed in parallel across multiple threads, which can significantly improve performance on multi-core processors.

3. **Loop Fusion:**
   - The original code has two separate loops for updating `A` and `C`. The optimized code merges these loops into a single loop structure, which reduces the overhead of loop control and can improve cache locality.

4. **Conditional Execution:**
   - The optimized code includes conditional checks (`if` statements) to ensure that the loop iterations are only executed when necessary. This can help to avoid unnecessary computations and improve performance.

### Application to the New Code:

- **Loop Tiling:**
  - The new code is optimized by applying loop tiling to the nested loops. The outer loops are tiled using variables `t1`, `t2`, `t3`, and `t4`, which represent the tile indices. This helps to improve cache utilization by processing the matrix in smaller blocks.

- **Parallelization:**
  - The outermost loop is parallelized using the `#pragma omp parallel for` directive, allowing the loop iterations to be executed in parallel across multiple threads.

- **Loop Fusion:**
  - The inner loop over `k` is kept within the tiled loops, ensuring that the computation for each element of the matrix is done in a single loop structure, which can improve cache locality and reduce loop control overhead.

By applying these techniques, the new code is expected to perform better in terms of cache utilization and parallel execution, leading to improved performance.*/

int t1, t2, t3, t4;
int lbp, ubp, lb2, ub2;
register int lbv, ubv;
lbp = 0;
ubp = floord(_PB_N-1, 32);
#pragma omp parallel for private(lbv, ubv, t2, t3, t4)
for (t1 = lbp; t1 <= ubp; t1++) {
    for (t2 = 0; t2 <= min(floord(_PB_N-1, 32), floord(999*t1+_PB_N-1, 32)); t2++) {
        for (t3 = max(0, 32*t1); t3 <= min(_PB_N-1, 32*t1+31); t3++) {
            for (t4 = 32*t2; t4 <= min(_PB_N-1, 32*t2+31); t4++) {
                for (int k = 0; k < _PB_N; k++) {
                    path[t3][t4] = path[t3][t4] < path[t3][k] + path[k][t4] ?
                        path[t3][t4] : path[t3][k] + path[k][t4];
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

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(path, DATA_TYPE, N, N, n, n);


  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(path));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_floyd_warshall (n, POLYBENCH_ARRAY(path));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(path)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(path);

  return 0;
}
