/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* bicg.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "bicg.h"

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
void init_array (int m, int n,
                 DATA_TYPE POLYBENCH_2D(A,N,M,n,m),
                 DATA_TYPE POLYBENCH_1D(r,N,n),
                 DATA_TYPE POLYBENCH_1D(p,M,m))
{
  int i, j;

  for (i = 0; i < m; i++)
    p[i] = (DATA_TYPE)(i % (m + 5)) / (m + 5);
  for (i = 0; i < n; i++) {
    r[i] = (DATA_TYPE)(i % (n + 5)) / (n + 5);
    for (j = 0; j < m; j++)
      A[i][j] = (DATA_TYPE) (i*(j+6) % (n + 5))/(n + 5);
  }
}



/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int m, int n,
		 DATA_TYPE POLYBENCH_1D(s,M,m),
		 DATA_TYPE POLYBENCH_1D(q,N,n))

{
  int i;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("s");
  DATA_TYPE tmp_s = 0;
  for (i = 0; i < m; i++) {
    if (DUMP) {
      if (i % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
      fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, s[i]);
    }
    if (CHECKSUM) tmp_s += s[i];
  }
  if (CHECKSUM) {
    fprintf(POLYBENCH_DUMP_TARGET,"\nchecksum: ");
    fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, tmp_s);
  }
  POLYBENCH_DUMP_END("s");
  POLYBENCH_DUMP_BEGIN("q");
  DATA_TYPE tmp_q = 0;
  for (i = 0; i < n; i++) {
    if (DUMP) {
      if (i % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
      fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, q[i]);
    }
    if (CHECKSUM) tmp_q += q[i];
  }
  if (CHECKSUM) {
    fprintf(POLYBENCH_DUMP_TARGET,"\nchecksum: ");
    fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, tmp_q);
  }
  POLYBENCH_DUMP_END("q");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_bicg(int m, int n,
		 DATA_TYPE POLYBENCH_2D(A,N,M,n,m),
		 DATA_TYPE POLYBENCH_1D(s,M,m),
		 DATA_TYPE POLYBENCH_1D(q,N,n),
		 DATA_TYPE POLYBENCH_1D(p,M,m),
		 DATA_TYPE POLYBENCH_1D(r,N,n))
{
  int i, j;

#pragma scop
/*### Explanation:
1. **Reduction in Array Accesses**: The variable `temp_q` is introduced to accumulate the value of `q[i]` inside the inner loop. This reduces the number of array accesses to `q[i]` from `_PB_M` to 1.
2. **Loop Fusion**: The operations inside the inner loop are fused together, reducing the number of iterations and improving cache locality.
3. **Avoiding Redundant Initialization**: The initialization of `q[i]` is moved inside the loop and only assigned once after the inner loop completes, reducing redundant operations.*/

for (i = 0; i < _PB_M; i++)
    s[i] = 0;
for (i = 0; i < _PB_N; i++)
{
    double temp_q = SCALAR_VAL(0.0);
    for (j = 0; j < _PB_M; j++)
    {
        s[j] += r[i] * A[i][j];
        temp_q += A[i][j] * p[j];
    }
    q[i] = temp_q;
}
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int n = N;
  int m = M;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, N, M, n, m);
  POLYBENCH_1D_ARRAY_DECL(s, DATA_TYPE, M, m);
  POLYBENCH_1D_ARRAY_DECL(q, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(p, DATA_TYPE, M, m);
  POLYBENCH_1D_ARRAY_DECL(r, DATA_TYPE, N, n);

  /* Initialize array(s). */
  init_array (m, n,
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(r),
	      POLYBENCH_ARRAY(p));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_bicg (m, n,
	       POLYBENCH_ARRAY(A),
	       POLYBENCH_ARRAY(s),
	       POLYBENCH_ARRAY(q),
	       POLYBENCH_ARRAY(p),
	       POLYBENCH_ARRAY(r));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(m, n, POLYBENCH_ARRAY(s), POLYBENCH_ARRAY(q)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(s);
  POLYBENCH_FREE_ARRAY(q);
  POLYBENCH_FREE_ARRAY(p);
  POLYBENCH_FREE_ARRAY(r);

  return 0;
}
