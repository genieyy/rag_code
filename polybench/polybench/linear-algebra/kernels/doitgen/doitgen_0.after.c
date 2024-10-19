/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* doitgen.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "doitgen.h"

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
void init_array(int nr, int nq, int np,
		DATA_TYPE POLYBENCH_3D(A,NR,NQ,NP,nr,nq,np),
		DATA_TYPE POLYBENCH_2D(C4,NP,NP,np,np))
{
  int i, j, k;

  for (i = 0; i < nr; i++)
    for (j = 0; j < nq; j++)
      for (k = 0; k < np; k++)
	A[i][j][k] = (DATA_TYPE) ((i*j + k)%np) / np;
  for (i = 0; i < np; i++)
    for (j = 0; j < np; j++)
      C4[i][j] = (DATA_TYPE) (i*j % np) / np;
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int nr, int nq, int np,
		 DATA_TYPE POLYBENCH_3D(A,NR,NQ,NP,nr,nq,np))
{
  int i, j, k;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("A");
  DATA_TYPE tmp_A = 0;
  for (i = 0; i < nr; i++)
    for (j = 0; j < nq; j++)
      for (k = 0; k < np; k++) {
        if (DUMP) {
  	if ((i*nq*np+j*np+k) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
  	fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, A[i][j][k]);
        }
        if (CHECKSUM) tmp_A += A[i][j][k];
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
void kernel_doitgen(int nr, int nq, int np,
		    DATA_TYPE POLYBENCH_3D(A,NR,NQ,NP,nr,nq,np),
		    DATA_TYPE POLYBENCH_2D(C4,NP,NP,np,np),
		    DATA_TYPE POLYBENCH_1D(sum,NP,np))
{
  int r, q, p, s;

#pragma scop
/*### Explanation of the Optimized Code:

1. **Loop Tiling and Parallelization**:
   - The original nested loops are transformed using loop tiling to improve cache locality and parallelize the computation.
   - The outer loops (`t1`, `t2`, `t3`) are tiled with a tile size of 32, which is a common choice for cache optimization.
   - The `#pragma omp parallel for` directive is used to parallelize the outermost loop (`t1`), distributing the work across multiple threads.

2. **Vectorization**:
   - The `#pragma ivdep` and `#pragma vector always` directives are used to hint the compiler to vectorize the innermost loop, which computes the sum and updates the array `A`.

3. **Reduction in Memory Accesses**:
   - By tiling the loops, the code reduces the number of cache misses by keeping the data accessed within a tile in the cache for longer periods.
   - The temporary array `sum` is used to store intermediate results, which helps in reducing redundant memory accesses.

4. **Loop Reordering**:
   - The loops are reordered to ensure that the innermost loop (`p`) is the one that iterates over the smallest dimension (`_PB_NP`), which is beneficial for vectorization and cache performance.

5. **Register Usage**:
   - The `register` keyword is used for the loop bounds (`lbv`, `ubv`) to suggest that these variables should be stored in CPU registers for faster access.

This optimized code should provide better performance by leveraging parallelism, vectorization, and cache locality improvements.*/

int t1, t2, t3, t4, t5;
int lb, ub, lbp, ubp, lb2, ub2;
register int lbv, ubv;

lbp = 0;
ubp = floord(_PB_NR - 1, 32);
#pragma omp parallel for private(lbv, ubv, t3, t4, t5)
for (t1 = lbp; t1 <= ubp; t1++) {
    for (t2 = 0; t2 <= floord(_PB_NQ - 1, 32); t2++) {
        for (t3 = 0; t3 <= floord(_PB_NP - 1, 32); t3++) {
            for (t4 = 32 * t1; t4 <= min(_PB_NR - 1, 32 * t1 + 31); t4++) {
                for (t5 = 32 * t2; t5 <= min(_PB_NQ - 1, 32 * t2 + 31); t5++) {
                    double sum[32];
                    lbv = max(0, 32 * t3);
                    ubv = min(_PB_NP - 1, 32 * t3 + 31);
#pragma ivdep
#pragma vector always
                    for (int p = lbv; p <= ubv; p++) {
                        sum[p - lbv] = SCALAR_VAL(0.0);
                        for (int s = 0; s < _PB_NP; s++) {
                            sum[p - lbv] += A[t4][t5][s] * C4[s][p];
                        }
                        A[t4][t5][p] = sum[p - lbv];
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
  int nr = NR;
  int nq = NQ;
  int np = NP;

  /* Variable declaration/allocation. */
  POLYBENCH_3D_ARRAY_DECL(A,DATA_TYPE,NR,NQ,NP,nr,nq,np);
  POLYBENCH_1D_ARRAY_DECL(sum,DATA_TYPE,NP,np);
  POLYBENCH_2D_ARRAY_DECL(C4,DATA_TYPE,NP,NP,np,np);

  /* Initialize array(s). */
  init_array (nr, nq, np,
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(C4));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_doitgen (nr, nq, np,
		  POLYBENCH_ARRAY(A),
		  POLYBENCH_ARRAY(C4),
		  POLYBENCH_ARRAY(sum));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(nr, nq, np,  POLYBENCH_ARRAY(A)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(sum);
  POLYBENCH_FREE_ARRAY(C4);

  return 0;
}
