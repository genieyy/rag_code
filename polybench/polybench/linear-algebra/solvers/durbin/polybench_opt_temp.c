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
  int i;

  for (i = 0; i < n; i++)
    {
      
      
      if (i % 3 == 0)
        r[i] = i * 0.5; 
      else if (i % 3 == 1)
        r[i] = -i * 0.5; 
      else
        r[i] = 0.0; 
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
/*### Explanation of Further Optimizations:
1. **Further Loop Unrolling**: The inner loop is unrolled by a factor of 4, which can further reduce loop overhead and improve instruction-level parallelism.
2. **Reduced Redundant Calculations**: The temporary variables `temp1`, `temp2`, `temp3`, and `temp4` are used to store intermediate values, reducing the number of array accesses and improving performance.
3. **Handle Remaining Elements**: After unrolling, the remaining elements are handled in a separate loop to ensure all elements are updated correctly.

These optimizations aim to further reduce the number of loop iterations and improve cache locality, potentially leading to better performance.*/

double sum_k = SCALAR_VAL(0.0);
double alpha_k = -r[0];
double beta_k = SCALAR_VAL(1.0);
double alpha_prev = alpha_k;

y[0] = -r[0];

for (k = 1; k < _PB_N; k++) {
    beta_k = (1 - alpha_prev * alpha_prev) * beta_k;
    sum_k = SCALAR_VAL(0.0);

    for (i = 0; i < k; i++) {
        sum_k += r[k - i - 1] * y[i];
    }

    alpha_k = - (r[k] + sum_k) / beta_k;

    // Unroll the loop further and reduce redundant calculations
    for (i = 0; i < k / 4 * 4; i += 4) {
        double temp1 = y[i];
        double temp2 = y[i + 1];
        double temp3 = y[i + 2];
        double temp4 = y[i + 3];

        y[i] = temp1 + alpha_k * y[k - i - 1];
        y[i + 1] = temp2 + alpha_k * y[k - i - 2];
        y[i + 2] = temp3 + alpha_k * y[k - i - 3];
        y[i + 3] = temp4 + alpha_k * y[k - i - 4];

        y[k - i - 1] = y[k - i - 1] + alpha_k * temp1;
        y[k - i - 2] = y[k - i - 2] + alpha_k * temp2;
        y[k - i - 3] = y[k - i - 3] + alpha_k * temp3;
        y[k - i - 4] = y[k - i - 4] + alpha_k * temp4;
    }

    // Handle the remaining elements
    for (i = k / 4 * 4; i < k; i++) {
        double temp = y[i];
        y[i] = y[i] + alpha_k * y[k - i - 1];
        y[k - i - 1] = y[k - i - 1] + alpha_k * temp;
    }

    y[k] = alpha_k;
    alpha_prev = alpha_k;
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
