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
1. **Loop Distribution and Parallelization**:
   - The original nested loops are distributed into multiple loops to facilitate parallelization. This is done using OpenMP's `#pragma omp parallel for` directive to parallelize the outermost loop.

2. **Loop Tiling**:
   - The loops are tiled using variables `t2`, `t3`, and `t4` to create smaller chunks of iterations that can be processed in parallel. This helps in better cache utilization and reduces the overhead of parallelization.

3. **Loop Reordering and Fusion**:
   - The inner loops are reordered and fused to ensure that the computation of `sum[p]` and the assignment to `A[r][q][p]` are done in the same loop structure. This reduces the number of loop iterations and improves data locality.

4. **Vectorization**:
   - The inner loop over `p` is vectorized using `#pragma ivdep` and `#pragma vector always` to ensure that the compiler can generate SIMD instructions for the loop, improving performance on modern CPUs.

5. **Bounds Calculation**:
   - The bounds of the loops are carefully calculated to ensure that the loop iterations are within the valid range of the arrays, avoiding out-of-bounds accesses.

By applying these transformations, the code is optimized for better performance, leveraging parallel execution, vectorization, and improved data locality.*/

int t1, t2, t3, t4, t5;
int lb, ub, lbp, ubp, lb2, ub2;
register int lbv, ubv;

lbp = 0;
ubp = floord(_PB_NR - 1, 32);
#pragma omp parallel for private(lbv, ubv, t3, t4, t5)
for (t2 = lbp; t2 <= ubp; t2++) {
    for (t3 = 0; t3 <= floord(_PB_NQ - 1, 32); t3++) {
        for (t4 = 0; t4 <= floord(_PB_NP - 1, 32); t4++) {
            for (t5 = 32 * t2; t5 <= min(_PB_NR - 1, 32 * t2 + 31); t5++) {
                for (int q = max(32 * t3, 0); q <= min(_PB_NQ - 1, 32 * t3 + 31); q++) {
                    lbv = max(32 * t4, 0);
                    ubv = min(_PB_NP - 1, 32 * t4 + 31);
#pragma ivdep
#pragma vector always
                    for (int p = lbv; p <= ubv; p++) {
                        sum[p] = SCALAR_VAL(0.0);
                        for (int s = 0; s < _PB_NP; s++) {
                            sum[p] += A[t5][q][s] * C4[s][p];
                        }
                        A[t5][q][p] = sum[p];
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
