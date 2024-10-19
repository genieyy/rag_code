/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* syr2k.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#include "syr2k.h"

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
void init_array(int n, int m,
		DATA_TYPE *alpha,
		DATA_TYPE *beta,
		DATA_TYPE POLYBENCH_2D(C,N,N,n,n),
		DATA_TYPE POLYBENCH_2D(A,N,M,n,m),
		DATA_TYPE POLYBENCH_2D(B,N,M,n,m))
{
  int i, j;

  *alpha = 1.5;
  *beta = 1.2;
  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++) {
      A[i][j] = (DATA_TYPE) ((i*j+1)%n) / n;
      B[i][j] = (DATA_TYPE) ((i*j+2)%m) / m;
    }
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      C[i][j] = (DATA_TYPE) ((i*j+3)%n) / m;
    }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_2D(C,N,N,n,n))
{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("C");
  DATA_TYPE tmp_C = 0;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      if (DUMP) {
  	if ((i * n + j) % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
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
void kernel_syr2k(int n, int m,
		  DATA_TYPE alpha,
		  DATA_TYPE beta,
		  DATA_TYPE POLYBENCH_2D(C,N,N,n,n),
		  DATA_TYPE POLYBENCH_2D(A,N,M,n,m),
		  DATA_TYPE POLYBENCH_2D(B,N,M,n,m))
{
  int i, j, k;

//BLAS PARAMS
//UPLO  = 'L'
//TRANS = 'N'
//A is NxM
//B is NxM
//C is NxN
#pragma scop
/**/

int t1, t2, t3, t4;
int lb, ub, lbp, ubp, lb2, ub2;
register int lbv, ubv;

for (int t1 = 0; t1 <= floord(_PB_N - 1, 32); t1++) {
    lbp = max(0, ceild(32 * t1 - _PB_N + 1, 32));
    ubp = min(min(floord(_PB_N - 1, 16), floord(32 * t1 + _PB_N + 30, 64)), t1);
#pragma omp parallel for private(lbv, ubv, t3, t4)
    for (int t2 = lbp; t2 <= ubp; t2++) {
        if (t2 <= floord(_PB_N - 1, 32)) {
            for (int t3 = max(0, 32 * t1 - 32 * t2); t3 <= min(32 * t1 - 32 * t2 + 31, t2); t3++) {
                for (int t4 = 32 * t2; t4 <= min(t3, 32 * t2 + 31); t4++) {
                    C[t3][t4] *= beta;
                }
            }
        }
        if (t2 >= ceild(_PB_N - 1, 32)) {
            for (int t3 = max(0, 32 * t1 - 32 * t2); t3 <= min(32 * t1 - 32 * t2 + 31, t2); t3++) {
                for (int t4 = 32 * t2; t4 <= min(t3, 32 * t2 + 31); t4++) {
                    C[t3][t4] *= beta;
                }
            }
        }
    }
}

for (int t1 = 0; t1 <= floord(_PB_N - 1, 32); t1++) {
    lbp = max(0, ceild(32 * t1 - _PB_N + 1, 32));
    ubp = min(min(floord(_PB_N - 1, 16), floord(32 * t1 + _PB_N + 30, 64)), t1);
#pragma omp parallel for private(lbv, ubv, t3, t4)
    for (int t2 = lbp; t2 <= ubp; t2++) {
        if (t2 <= floord(_PB_N - 1, 32)) {
            for (int t3 = max(0, 32 * t1 - 32 * t2); t3 <= min(32 * t1 - 32 * t2 + 31, t2); t3++) {
                for (int t4 = 32 * t2; t4 <= min(t3, 32 * t2 + 31); t4++) {
                    for (int k = 0; k < _PB_M; k++) {
                        C[t3][t4] += A[t4][k] * alpha * B[t3][k] + B[t4][k] * alpha * A[t3][k];
                    }
                }
            }
        }
        if (t2 >= ceild(_PB_N - 1, 32)) {
            for (int t3 = max(0, 32 * t1 - 32 * t2); t3 <= min(32 * t1 - 32 * t2 + 31, t2); t3++) {
                for (int t4 = 32 * t2; t4 <= min(t3, 32 * t2 + 31); t4++) {
                    for (int k = 0; k < _PB_M; k++) {
                        C[t3][t4] += A[t4][k] * alpha * B[t3][k] + B[t4][k] * alpha * A[t3][k];
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
  int n = N;
  int m = M;

  /* Variable declaration/allocation. */
  DATA_TYPE alpha;
  DATA_TYPE beta;
  POLYBENCH_2D_ARRAY_DECL(C,DATA_TYPE,N,N,n,n);
  POLYBENCH_2D_ARRAY_DECL(A,DATA_TYPE,N,M,n,m);
  POLYBENCH_2D_ARRAY_DECL(B,DATA_TYPE,N,M,n,m);

  /* Initialize array(s). */
  init_array (n, m, &alpha, &beta,
	      POLYBENCH_ARRAY(C),
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(B));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_syr2k (n, m,
		alpha, beta,
		POLYBENCH_ARRAY(C),
		POLYBENCH_ARRAY(A),
		POLYBENCH_ARRAY(B));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(C)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(C);
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(B);

  return 0;
}
