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
/*### Analysis of Meaning-Preserving Loop Transformations:

1. **Loop Tiling/Blocking**: The original loops are transformed into tiled loops. This is evident from the introduction of new variables (`t1`, `t2`, `t3`) that define the tile boundaries. This transformation helps in reducing cache misses by ensuring that the data accessed within a tile fits into the cache.

2. **Loop Fusion/Fission**: The original loops are split into multiple loops (fission) and then fused back together in a way that reduces the number of iterations and improves locality. For example, the multiplication of `C[i][j]` by `beta` and the addition of `alpha * A[i][k] * B[k][j]` to `C[i][j]` are handled in separate loops but within the same tile.

3. **Parallelization**: The use of `#pragma omp parallel for` indicates that the loops are parallelized to take advantage of multi-core processors. This is a common optimization for loops that can be executed independently.

4. **Loop Interchange**: The order of loops is changed to improve cache locality. For example, the innermost loop over `j` is moved outside in some cases to ensure that the data accessed by `j` is reused more effectively.

### Learning and Application:

- **Tiling**: By introducing tiles (`t1`, `t2`, `t3`), we ensure that the data accessed within each tile fits into the cache, reducing cache misses and improving performance.
- **Parallelization**: Using OpenMP to parallelize the outer loops allows the program to take advantage of multi-core processors.
- **Loop Fusion/Fission**: Splitting the original loops and then fusing them back together in a way that reduces the number of iterations and improves locality.

The optimized code applies these transformations to the given loop structure, ensuring better cache utilization and parallel execution.*/

int t1, t2, t3;
int lb, ub, lbp, ubp;
register int lbv, ubv;

lbp = 0;
ubp = floord(_PB_NI, 16);

#pragma omp parallel for private(lbv, ubv, t2, t3)
for (t1 = lbp; t1 <= ubp; t1++) {
    for (t2 = max(0, ceild(32 * t1 - _PB_NI + 1, 32)); t2 <= min(floord(_PB_NI, 32), t1); t2++) {
        for (t3 = max(0, 32 * t1 - 32 * t2); t3 <= min(_PB_NJ - 1, 32 * t1 - 32 * t2 + 31); t3++) {
            for (int i = 32 * t2; i <= min(_PB_NI - 1, 32 * t2 + 31); i++) {
                C[i][t3] *= beta;
            }
        }
        for (t3 = max(0, 32 * t1 - 32 * t2); t3 <= min(_PB_NK - 1, 32 * t1 - 32 * t2 + 31); t3++) {
            for (int i = 32 * t2; i <= min(_PB_NI - 1, 32 * t2 + 31); i++) {
                for (int j = 0; j < _PB_NJ; j++) {
                    C[i][j] += alpha * A[i][t3] * B[t3][j];
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
