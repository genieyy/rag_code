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

  
  for (i = 0; i < n; i++) {
    seq[i] = (base)((i % 2 == 0) ? 0 : 1); 
  }

  
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (i == j) {
        table[i][j] = 0; 
      } else if (i > j) {
        table[i][j] = -1; 
      } else {
        
        if ((i + j) % 2 == 0) {
          table[i][j] = (DATA_TYPE)(i + j); 
        } else {
          table[i][j] = (DATA_TYPE)(i * j); 
        }
      }
    }
  }
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
/*### Analysis and Transformation Methods Used:

1. **Precompute Common Conditions**:
   - The conditions `t2-1 >= 0` and `t1+1 < _PB_N` are precomputed to avoid redundant checks inside the loop. This reduces the number of conditional checks and improves performance.

2. **Loop Unrolling**:
   - The inner loop is unrolled by a factor of 4 to reduce the loop overhead. This can improve performance by reducing the number of loop iterations and allowing the CPU to execute more instructions per cycle.

3. **Reduced Redundant Checks**:
   - By precomputing the conditions, the code avoids redundant checks, which can improve performance by reducing the number of conditional branches.

### Performance Improvement:
- The optimized code reduces the number of conditional checks by precomputing common conditions, which can improve performance by reducing the number of branch mispredictions.
- The inner loop is unrolled to reduce the loop overhead, which can improve performance by allowing the CPU to execute more instructions per cycle.
- Overall, these optimizations can lead to a significant performance improvement, especially for large values of `_PB_N`.*/

/**/

int t1, t2, t3;
int lb, ub, lbp, ubp, lb2, ub2;
register int lbv, ubv;

for (int t1 = _PB_N-1; t1 >= 0; t1--) {
    for (int t2 = t1+1; t2 < _PB_N; t2++) {
        // Precompute common conditions to reduce redundant checks
        int t2_minus_1_ge_0 = (t2-1 >= 0);
        int t1_plus_1_lt_N = (t1+1 < _PB_N);

        if (t2_minus_1_ge_0)
            table[t1][t2] = max_score(table[t1][t2], table[t1][t2-1]);
        if (t1_plus_1_lt_N)
            table[t1][t2] = max_score(table[t1][t2], table[t1+1][t2]);

        if (t2_minus_1_ge_0 && t1_plus_1_lt_N) {
            if (t1 < t2-1)
                table[t1][t2] = max_score(table[t1][t2], table[t1+1][t2-1] + match(seq[t1], seq[t2]));
            else
                table[t1][t2] = max_score(table[t1][t2], table[t1+1][t2-1]);
        }

        // Unroll the inner loop slightly to reduce loop overhead
        int t3;
        for (t3 = t1+1; t3 + 4 < t2; t3 += 4) {
            table[t1][t2] = max_score(table[t1][t2], table[t1][t3] + table[t3+1][t2]);
            table[t1][t2] = max_score(table[t1][t2], table[t1][t3+1] + table[t3+2][t2]);
            table[t1][t2] = max_score(table[t1][t2], table[t1][t3+2] + table[t3+3][t2]);
            table[t1][t2] = max_score(table[t1][t2], table[t1][t3+3] + table[t3+4][t2]);
        }
        // Handle the remaining iterations
        for (; t3 < t2; t3++) {
            table[t1][t2] = max_score(table[t1][t2], table[t1][t3] + table[t3+1][t2]);
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
