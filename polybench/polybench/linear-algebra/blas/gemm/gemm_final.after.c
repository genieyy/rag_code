/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* gemm.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "gemm.h"

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
void init_array(int ni, int nj, int nk,
		DATA_TYPE *alpha,
		DATA_TYPE *beta,
		DATA_TYPE POLYBENCH_2D(C,NI,NJ,ni,nj),
		DATA_TYPE POLYBENCH_2D(A,NI,NK,ni,nk),
		DATA_TYPE POLYBENCH_2D(B,NK,NJ,nk,nj))
{
  int i, j;

  *alpha = 1.5;
  *beta = 1.2;
  for (i = 0; i < ni; i++)
    for (j = 0; j < nj; j++)
      C[i][j] = (DATA_TYPE) ((i*j+1) % ni) / ni;
  for (i = 0; i < ni; i++)
    for (j = 0; j < nk; j++)
      A[i][j] = (DATA_TYPE) (i*(j+1) % nk) / nk;
  for (i = 0; i < nk; i++)
    for (j = 0; j < nj; j++)
      B[i][j] = (DATA_TYPE) (i*(j+2) % nj) / nj;
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int ni, int nj,
		 DATA_TYPE POLYBENCH_2D(C,NI,NJ,ni,nj))
{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("C");
  DATA_TYPE tmp_C = 0;
  for (i = 0; i < ni; i++)
    for (j = 0; j < nj; j++) {
      if (DUMP) {
  	if ((i * ni + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
  	fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, C[i][j]);
      }
      if (CHECKSUM) tmp_C += C[i][j];
    }
  if (CHECKSUM) {
    fprintf(POLYBENCH_DUMP_TARGET,"\nchecksum: ");
    fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, tmp_C);
  }
  POLYBENCH_DUMP_END("C");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_gemm(int ni, int nj, int nk,
		 DATA_TYPE alpha,
		 DATA_TYPE beta,
		 DATA_TYPE POLYBENCH_2D(C,NI,NJ,ni,nj),
		 DATA_TYPE POLYBENCH_2D(A,NI,NK,ni,nk),
		 DATA_TYPE POLYBENCH_2D(B,NK,NJ,nk,nj))
{
  int i, j, k;

//BLAS PARAMS
//TRANSA = 'N'
//TRANSB = 'N'
// => Form C := alpha*A*B + beta*C,
//A is NIxNK
//B is NKxNJ
//C is NIxNJ
#pragma scop
/*### Explanation:
1. **Parallelization and Collapsing Loops**: 
   - The `#pragma omp parallel` directive is used to parallelize the entire block.
   - The `#pragma omp for collapse(2)` and `#pragma omp for collapse(3)` directives are used to collapse the nested loops, reducing the overhead of parallelizing each loop individually. This can improve performance by allowing the OpenMP runtime to better distribute the work across threads.

2. **Loop Tiling**:
   - The code still uses loop tiling with a tile size of 32, which helps in better cache utilization and reduces cache misses.

3. **Reduction of Redundant Calculations**:
   - The expression `alpha * A[t3][k]` is computed once and stored in a temporary variable `temp`. This reduces the number of multiplications performed in the innermost loop, which can be a significant performance improvement, especially for large matrices.

4. **Register Usage**:
   - The use of `register` for `lbv` and `ubv` is retained to encourage the compiler to place these variables in CPU registers, which can speed up access.

This version of the code should provide better performance due to the combined effects of parallelization, loop collapsing, and reduction of redundant calculations.*/

#pragma omp parallel
{
    int t1, t2, t3;
    int lb, ub, lbp, ubp;
    register int lbv, ubv;

    lbp = 0;
    ubp = floord(_PB_NI, 32);

    #pragma omp for collapse(2)
    for (t1 = lbp; t1 <= ubp; t1++) {
        for (t2 = 0; t2 <= floord(_PB_NJ, 32); t2++) {
            for (t3 = 32 * t1; t3 <= min(_PB_NI - 1, 32 * t1 + 31); t3++) {
                for (int j = 32 * t2; j <= min(_PB_NJ - 1, 32 * t2 + 31); j++) {
                    C[t3][j] *= beta;
                }
            }
        }
    }

    #pragma omp for collapse(3)
    for (t1 = lbp; t1 <= ubp; t1++) {
        for (t2 = 0; t2 <= floord(_PB_NK, 32); t2++) {
            for (t3 = 32 * t1; t3 <= min(_PB_NI - 1, 32 * t1 + 31); t3++) {
                for (int k = 32 * t2; k <= min(_PB_NK - 1, 32 * t2 + 31); k++) {
                    double temp = alpha * A[t3][k];
                    for (int j = 0; j < _PB_NJ; j++) {
                        C[t3][j] += temp * B[k][j];
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
  int ni = NI;
  int nj = NJ;
  int nk = NK;

  /* Variable declaration/allocation. */
  DATA_TYPE alpha;
  DATA_TYPE beta;
  POLYBENCH_2D_ARRAY_DECL(C,DATA_TYPE,NI,NJ,ni,nj);
  POLYBENCH_2D_ARRAY_DECL(A,DATA_TYPE,NI,NK,ni,nk);
  POLYBENCH_2D_ARRAY_DECL(B,DATA_TYPE,NK,NJ,nk,nj);

  /* Initialize array(s). */
  init_array (ni, nj, nk, &alpha, &beta,
	      POLYBENCH_ARRAY(C),
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(B));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_gemm (ni, nj, nk,
	       alpha, beta,
	       POLYBENCH_ARRAY(C),
	       POLYBENCH_ARRAY(A),
	       POLYBENCH_ARRAY(B));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(ni, nj,  POLYBENCH_ARRAY(C)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(C);
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(B);

  return 0;
}
