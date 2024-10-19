#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <polybench.h>
#include "nussinov.h"
#include <omp.h>
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


/* Include polybench common header. */

/* Include benchmark-specific header. */

/* RNA bases represented as chars, range is [0,3] */
typedef char base;

#define match(b1, b2) (((b1)+(b2)) == 3 ? 1 : 0)
#define max_score(s1, s2) ((s1 >= s2) ? s1 : s2)
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
  int t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
for (t2=-1;t2<=floord(_PB_N-1,32);t2++) {
  lbp=max(0,t2);
  ubp=floord(_PB_N-1,32);
#pragma omp parallel for private(lbv,ubv,t5,t6,t7,t8,t9,t10)
  for (t4=lbp;t4<=ubp;t4++) {
    if ((t2 == -1) && (t4 == 78)) {
      table[(_PB_N-2)][(_PB_N-1)] = max_score(table[(_PB_N-2)][(_PB_N-1)], table[(_PB_N-2)][(_PB_N-1)-1]);;
      table[(_PB_N-2)][(_PB_N-1)] = max_score(table[(_PB_N-2)][(_PB_N-1)], table[(_PB_N-2)+1][(_PB_N-1)]);;
      table[(_PB_N-2)][(_PB_N-1)] = max_score(table[(_PB_N-2)][(_PB_N-1)], table[(_PB_N-2)+1][(_PB_N-1)-1]);;
    }
    if ((t2 == -1) && (t4 <= floord(_PB_N-32,32))) {
      table[(32*t4+30)][(32*t4+31)] = max_score(table[(32*t4+30)][(32*t4+31)], table[(32*t4+30)][(32*t4+31)-1]);;
      table[(32*t4+30)][(32*t4+31)] = max_score(table[(32*t4+30)][(32*t4+31)], table[(32*t4+30)+1][(32*t4+31)]);;
      table[(32*t4+30)][(32*t4+31)] = max_score(table[(32*t4+30)][(32*t4+31)], table[(32*t4+30)+1][(32*t4+31)-1]);;
    }
    for (t5=max(max(-_PB_N+3,32*t2-32*t4),-32*t4-29);t5<=min(min(0,-32*t4+1),32*t2-32*t4+31);t5++) {
      table[-t5][(-t5+1)] = max_score(table[-t5][(-t5+1)], table[-t5][(-t5+1)-1]);;
      table[-t5][(-t5+1)] = max_score(table[-t5][(-t5+1)], table[-t5+1][(-t5+1)]);;
      table[-t5][(-t5+1)] = max_score(table[-t5][(-t5+1)], table[-t5+1][(-t5+1)-1]);;
      for (t7=-t5+2;t7<=min(_PB_N-1,32*t4+31);t7++) {
        table[-t5][t7] = max_score(table[-t5][t7], table[-t5][t7-1]);;
        table[-t5][t7] = max_score(table[-t5][t7], table[-t5+1][t7]);;
        table[-t5][t7] = max_score(table[-t5][t7], table[-t5+1][t7-1]+match(seq[-t5], seq[t7]));;
        for (t9=-t5+1;t9<=t7-1;t9++) {
          table[-t5][t7] = max_score(table[-t5][t7], table[-t5][t9] + table[t9+1][t7]);;
        }
      }
    }
    for (t5=max(32*t2-32*t4,-32*t4+2);t5<=min(0,32*t2-32*t4+31);t5++) {
      for (t7=32*t4;t7<=min(_PB_N-1,32*t4+31);t7++) {
        table[-t5][t7] = max_score(table[-t5][t7], table[-t5][t7-1]);;
        table[-t5][t7] = max_score(table[-t5][t7], table[-t5+1][t7]);;
        table[-t5][t7] = max_score(table[-t5][t7], table[-t5+1][t7-1]+match(seq[-t5], seq[t7]));;
        for (t9=-t5+1;t9<=t7-1;t9++) {
          table[-t5][t7] = max_score(table[-t5][t7], table[-t5][t9] + table[t9+1][t7]);;
        }
      }
    }
  }
}
/* End of CLooG code */

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
