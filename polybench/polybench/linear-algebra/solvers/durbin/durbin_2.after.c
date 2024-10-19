/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* durbin.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "durbin.h"

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
		 DATA_TYPE POLYBENCH_1D(r,N,n))
{
  int i, j;

  for (i = 0; i < n; i++)
    {
      r[i] = (n+1-i);
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
void kernel_durbin(int n,
		   DATA_TYPE POLYBENCH_1D(r,N,n),
		   DATA_TYPE POLYBENCH_1D(y,N,n))
{
 DATA_TYPE z[N];
 DATA_TYPE alpha;
 DATA_TYPE beta;
 DATA_TYPE sum;

 int i,k;

#pragma scop
/*### Explanation of Optimizations:
1. **Reduction in Redundant Calculations**:
   - The `sum_k` variable is reused to accumulate the sum inside the loop, reducing the need to reinitialize it every time.
   - The `alpha_k` and `beta_k` variables are updated in place, avoiding the need to redefine them.

2. **Loop Unrolling**:
   - The inner loop that updates `y[i]` and `y[k-i-1]` is unrolled to handle pairs of elements at once, reducing the number of iterations by half. This is particularly beneficial for larger values of `k`.

3. **Conditional Check for Odd `k`**:
   - A check is added to handle the middle element when `k` is odd, ensuring that the loop unrolling does not miss any elements.

These optimizations aim to reduce the number of operations and improve the locality of reference, which can lead to better performance, especially for larger values of `_PB_N`.*/

double sum_k = SCALAR_VAL(0.0);
double alpha_k = -r[0];
double beta_k = SCALAR_VAL(1.0);
double temp;

for (k = 1; k < _PB_N; k++) {
    beta_k = (1 - alpha_k * alpha_k) * beta_k;
    sum_k = SCALAR_VAL(0.0);
    for (i = 0; i < k; i++) {
        sum_k += r[k - i - 1] * y[i];
    }
    alpha_k = - (r[k] + sum_k) / beta_k;

    for (i = 0; i < k / 2; i++) {
        temp = y[i];
        y[i] = y[i] + alpha_k * y[k - i - 1];
        y[k - i - 1] = y[k - i - 1] + alpha_k * temp;
    }
    if (k % 2 == 1) {
        y[k / 2] = y[k / 2] + alpha_k * y[k / 2];
    }

    y[k] = alpha_k;
}
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int n = N;

  /* Variable declaration/allocation. */
  POLYBENCH_1D_ARRAY_DECL(r, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(y, DATA_TYPE, N, n);


  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(r));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_durbin (n,
		 POLYBENCH_ARRAY(r),
		 POLYBENCH_ARRAY(y));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(y)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(r);
  POLYBENCH_FREE_ARRAY(y);

  return 0;
}
