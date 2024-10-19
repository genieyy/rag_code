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
   - The original code uses a simple triple-nested loop. The optimized code introduces loop tiling (also known as blocking) by dividing the loops into chunks of size 32. This technique improves cache locality by ensuring that the data accessed within each chunk fits better into the cache, reducing cache misses.

2. **Parallelization:**
   - The `#pragma omp parallel for` directive is used to parallelize the outermost loop. This allows multiple threads to work on different chunks of the data simultaneously, leveraging multi-core processors to improve performance.

3. **Loop Fusion:**
   - Although not explicitly shown in the provided code, the original example demonstrates loop fusion by combining two loops into a single loop structure. This can be beneficial for reducing loop overhead and improving data locality.

4. **Conditional Execution:**
   - The original example uses conditional checks within the loop to handle specific cases (e.g., `if (t2 == 187)`). This can be useful for optimizing specific scenarios but should be used judiciously to avoid complicating the code.

### Application to the New Code:

- **Loop Tiling:** The new code applies loop tiling by dividing the `k`, `i`, and `j` loops into chunks of size 32. This ensures that the data accessed within each chunk is more likely to be in the cache, improving performance.
  
- **Parallelization:** The outermost loop is parallelized using OpenMP, allowing multiple threads to work on different chunks of the data simultaneously.

- **Nested Loops:** The inner loops iterate over the smaller chunks, ensuring that the data accessed within each chunk is processed efficiently.

This approach leverages the lessons learned from the original example to improve the performance of the new code by optimizing cache usage and parallelizing the computation.*/

#pragma omp parallel for private(i, j, k)
for (int k = 0; k < _PB_N; k += 32) {
    for (int i = 0; i < _PB_N; i += 32) {
        for (int j = 0; j < _PB_N; j += 32) {
            for (int kk = k; kk < min(k + 32, _PB_N); kk++) {
                for (int ii = i; ii < min(i + 32, _PB_N); ii++) {
                    for (int jj = j; jj < min(j + 32, _PB_N); jj++) {
                        path[ii][jj] = path[ii][jj] < path[ii][kk] + path[kk][jj] ?
                                       path[ii][jj] : path[ii][kk] + path[kk][jj];
                    }
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
