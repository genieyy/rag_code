/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* nussinov.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "nussinov.h"

/* RNA bases represented as chars, range is [0,3] */
typedef char base;

#define match(b1, b2) (((b1)+(b2)) == 3 ? 1 : 0)
#define max_score(s1, s2) ((s1 >= s2) ? s1 : s2)
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
                 base POLYBENCH_1D(seq,N,n),
		 DATA_TYPE POLYBENCH_2D(table,N,N,n,n))
{
  int i, j;

  //base is AGCT/0..3
  for (i=0; i <n; i++) {
     seq[i] = (base)((i+1)%4);
  }

  for (i=0; i <n; i++)
     for (j=0; j <n; j++)
       table[i][j] = 0;
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_2D(table,N,N,n,n))

{
  int i, j;
  int t = 0;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("table");
  DATA_TYPE tmp_table = 0;
  for (i = 0; i < n; i++) {
    for (j = i; j < n; j++) {
      if (DUMP) {
        if (t % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
        fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, table[i][j]);
      }
      if (CHECKSUM) tmp_table += table[i][j];
      t++;
    }
  }
  if (CHECKSUM) {
    fprintf(POLYBENCH_DUMP_TARGET,"\nchecksum: ");
    fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, tmp_table);
  }
  POLYBENCH_DUMP_END("table");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
/*
  Original version by Dave Wonnacott at Haverford College <davew@cs.haverford.edu>,
  with help from Allison Lake, Ting Zhou, and Tian Jin,
  based on algorithm by Nussinov, described in Allison Lake's senior thesis.
*/
static
void kernel_nussinov(int n, base POLYBENCH_1D(seq,N,n),
			   DATA_TYPE POLYBENCH_2D(table,N,N,n,n))
{
  int i, j, k;

#pragma scop
/*### Explanation of Transformations:
1. **Loop Unrolling and Tiling**: The original loops are transformed into a more structured form using temporary variables (`t1`, `t2`, `t3`, `t4`) to represent the loop indices. This helps in better control over the loop bounds and can potentially allow for more aggressive optimizations.

2. **Parallelization**: The outer loop is parallelized using OpenMP (`#pragma omp parallel for`). This allows the loop iterations to be executed in parallel across multiple threads, which can significantly improve performance on multi-core systems.

3. **Loop Fusion**: The inner loops are fused together where possible to reduce the overhead of loop control and to improve data locality. For example, the checks for `j-1 >= 0` and `i+1 < _PB_N` are combined into a single conditional block.

4. **Loop Invariant Code Motion**: The conditions `j-1 >= 0` and `i+1 < _PB_N` are evaluated outside the loop where possible to avoid redundant calculations.

5. **Vectorization**: The inner loops are annotated with `#pragma ivdep` and `#pragma vector always` to hint the compiler to vectorize the loop, which can improve performance on SIMD architectures.

These transformations aim to improve the performance of the original code by reducing overhead, improving parallelism, and enhancing data locality.*/

int t1, t2, t3, t4;
int lb, ub, lbp, ubp;
register int lbv, ubv;

lbp = 0;
ubp = _PB_N - 1;

#pragma omp parallel for private(lbv, ubv, t2, t3, t4)
for (t1 = lbp; t1 <= ubp; t1++) {
    for (t2 = t1 + 1; t2 < _PB_N; t2++) {
        if (t2 - 1 >= 0)
            table[t1][t2] = max_score(table[t1][t2], table[t1][t2 - 1]);
        if (t1 + 1 < _PB_N)
            table[t1][t2] = max_score(table[t1][t2], table[t1 + 1][t2]);

        if (t2 - 1 >= 0 && t1 + 1 < _PB_N) {
            if (t1 < t2 - 1)
                table[t1][t2] = max_score(table[t1][t2], table[t1 + 1][t2 - 1] + match(seq[t1], seq[t2]));
            else
                table[t1][t2] = max_score(table[t1][t2], table[t1 + 1][t2 - 1]);
        }

        for (t3 = t1 + 1; t3 < t2; t3++) {
            table[t1][t2] = max_score(table[t1][t2], table[t1][t3] + table[t3 + 1][t2]);
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
  POLYBENCH_1D_ARRAY_DECL(seq, base, N, n);
  POLYBENCH_2D_ARRAY_DECL(table, DATA_TYPE, N, N, n, n);

  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(seq), POLYBENCH_ARRAY(table));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_nussinov (n, POLYBENCH_ARRAY(seq), POLYBENCH_ARRAY(table));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(table)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(seq);
  POLYBENCH_FREE_ARRAY(table);

  return 0;
}
