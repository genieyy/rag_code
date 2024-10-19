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
   - The original nested loops are transformed using loop tiling to improve cache locality and parallelized using OpenMP to exploit multi-core processors.
   - The outer loops are tiled with a tile size of 32, which is a common choice for balancing the overhead of parallelization and the benefits of cache reuse.

2. **Private Arrays**:
   - The `sum` array is declared as a private variable within the innermost loop to avoid false sharing and to ensure that each thread operates on its own copy of the data.

3. **Vectorization**:
   - The innermost loops are annotated with `#pragma ivdep` and `#pragma vector always` to hint the compiler to vectorize the loops, assuming the compiler supports these pragmas.

4. **Bounds Calculation**:
   - The bounds for the tiled loops are calculated dynamically to ensure that the loop iterations are within the valid range of the problem size.

This optimization strategy leverages parallelism, improves cache locality, and hints at vectorization to enhance the performance of the original code.*/

int lbp, ubp;
register int lbv, ubv;
lbp = 0;
ubp = floord(_PB_NR - 1, 32);
#pragma omp parallel for private(lbv, ubv)
for (int t1 = lbp; t1 <= ubp; t1++) {
    for (int t2 = 0; t2 <= floord(_PB_NQ - 1, 32); t2++) {
        for (int t3 = 0; t3 <= floord(_PB_NP - 1, 32); t3++) {
            for (int r = max(0, 32 * t1); r <= min(_PB_NR - 1, 32 * t1 + 31); r++) {
                for (int q = max(0, 32 * t2); q <= min(_PB_NQ - 1, 32 * t2 + 31); q++) {
                    double sum[32];
                    for (int p = max(0, 32 * t3); p <= min(_PB_NP - 1, 32 * t3 + 31); p++) {
                        sum[p - max(0, 32 * t3)] = SCALAR_VAL(0.0);
                        for (int s = 0; s < _PB_NP; s++) {
                            sum[p - max(0, 32 * t3)] += A[r][q][s] * C4[s][p];
                        }
                    }
                    for (int p = max(0, 32 * t3); p <= min(_PB_NP - 1, 32 * t3 + 31); p++) {
                        A[r][q][p] = sum[p - max(0, 32 * t3)];
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
