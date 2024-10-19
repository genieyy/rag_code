#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <polybench.h>
#include "adi.h"
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
/* adi.c: this file is part of PolyBench/C */


/* Include polybench common header. */

/* Include benchmark-specific header. */

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
		 DATA_TYPE POLYBENCH_2D(u,N,N,n,n))
{
  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      {
	u[i][j] =  (DATA_TYPE)(i + n-j) / n;
      }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_2D(u,N,N,n,n))

{
  int i, j;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("u");
  DATA_TYPE tmp_u = 0;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      if (DUMP) {
        if ((i * n + j) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
        fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, u[i][j]);
      }
      if (CHECKSUM) tmp_u += u[i][j];
    }
  if (CHECKSUM) {
    fprintf(POLYBENCH_DUMP_TARGET,"\nchecksum: ");
    fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, tmp_u);
  }
  POLYBENCH_DUMP_END("u");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
/* Based on a Fortran code fragment from Figure 5 of
 * "Automatic Data and Computation Decomposition on Distributed Memory Parallel Computers"
 * by Peizong Lee and Zvi Meir Kedem, TOPLAS, 2002
 */
static
void kernel_adi(int tsteps, int n,
		DATA_TYPE POLYBENCH_2D(u,N,N,n,n),
		DATA_TYPE POLYBENCH_2D(v,N,N,n,n),
		DATA_TYPE POLYBENCH_2D(p,N,N,n,n),
		DATA_TYPE POLYBENCH_2D(q,N,N,n,n))
{
  int t, i, j;
  DATA_TYPE DX, DY, DT;
  DATA_TYPE B1, B2;
  DATA_TYPE mul1, mul2;
  DATA_TYPE a, b, c, d, e, f;


  DX = SCALAR_VAL(1.0)/(DATA_TYPE)_PB_N;
  DY = SCALAR_VAL(1.0)/(DATA_TYPE)_PB_N;
  DT = SCALAR_VAL(1.0)/(DATA_TYPE)_PB_TSTEPS;
  B1 = SCALAR_VAL(2.0);
  B2 = SCALAR_VAL(1.0);
  mul1 = B1 * DT / (DX * DX);
  mul2 = B2 * DT / (DY * DY);

  a = -mul1 /  SCALAR_VAL(2.0);
  b = SCALAR_VAL(1.0)+mul1;
  c = a;
  d = -mul2 / SCALAR_VAL(2.0);
  e = SCALAR_VAL(1.0)+mul2;
  f = d;
  
  int t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
for (t2=1;t2<=_PB_TSTEPS;t2++) {
  lbp=0;
  ubp=floord(_PB_N-2,32);
#pragma omp parallel for private(lbv,ubv,t5,t6,t7,t8,t9,t10)
  for (t4=lbp;t4<=ubp;t4++) {
    lbv=max(1,32*t4);
    ubv=min(_PB_N-2,32*t4+31);
#pragma ivdep
#pragma vector always
    for (t7=lbv;t7<=ubv;t7++) {
      v[0][t7] = SCALAR_VAL(1.0);;
    }
    lbv=max(1,32*t4);
    ubv=min(_PB_N-2,32*t4+31);
#pragma ivdep
#pragma vector always
    for (t7=lbv;t7<=ubv;t7++) {
      p[t7][0] = SCALAR_VAL(0.0);;
    }
    lbv=max(1,32*t4);
    ubv=min(_PB_N-2,32*t4+31);
#pragma ivdep
#pragma vector always
    for (t7=lbv;t7<=ubv;t7++) {
      q[t7][0] = v[0][t7];;
    }
    for (t6=0;t6<=floord(_PB_N-2,32);t6++) {
      for (t7=max(1,32*t4);t7<=min(_PB_N-2,32*t4+31);t7++) {
        for (t9=max(1,32*t6);t9<=min(_PB_N-2,32*t6+31);t9++) {
          p[t7][t9] = -c / (a*p[t7][t9-1]+b);;
          q[t7][t9] = (-d*u[t9][t7-1]+(SCALAR_VAL(1.0)+SCALAR_VAL(2.0)*d)*u[t9][t7] - f*u[t9][t7+1]-a*q[t7][t9-1])/(a*p[t7][t9-1]+b);;
        }
      }
    }
    lbv=max(1,32*t4);
    ubv=min(_PB_N-2,32*t4+31);
#pragma ivdep
#pragma vector always
    for (t7=lbv;t7<=ubv;t7++) {
      v[_PB_N-1][t7] = SCALAR_VAL(1.0);;
    }
    for (t6=ceild(-_PB_N-29,32);t6<=-1;t6++) {
      for (t7=max(1,32*t4);t7<=min(_PB_N-2,32*t4+31);t7++) {
        for (t9=max(32*t6,-_PB_N+2);t9<=32*t6+31;t9++) {
          v[-t9][t7] = p[t7][-t9] * v[-t9+1][t7] + q[t7][-t9];;
        }
      }
    }
  }
  lbp=0;
  ubp=floord(_PB_N-2,32);
#pragma omp parallel for private(lbv,ubv,t5,t6,t7,t8,t9,t10)
  for (t4=lbp;t4<=ubp;t4++) {
    lbv=max(1,32*t4);
    ubv=min(_PB_N-2,32*t4+31);
#pragma ivdep
#pragma vector always
    for (t7=lbv;t7<=ubv;t7++) {
      u[t7][0] = SCALAR_VAL(1.0);;
    }
    lbv=max(1,32*t4);
    ubv=min(_PB_N-2,32*t4+31);
#pragma ivdep
#pragma vector always
    for (t7=lbv;t7<=ubv;t7++) {
      p[t7][0] = SCALAR_VAL(0.0);;
    }
    lbv=max(1,32*t4);
    ubv=min(_PB_N-2,32*t4+31);
#pragma ivdep
#pragma vector always
    for (t7=lbv;t7<=ubv;t7++) {
      q[t7][0] = u[t7][0];;
    }
    for (t6=0;t6<=floord(_PB_N-2,32);t6++) {
      for (t7=max(1,32*t4);t7<=min(_PB_N-2,32*t4+31);t7++) {
        for (t9=max(1,32*t6);t9<=min(_PB_N-2,32*t6+31);t9++) {
          p[t7][t9] = -f / (d*p[t7][t9-1]+e);;
          q[t7][t9] = (-a*v[t7-1][t9]+(SCALAR_VAL(1.0)+SCALAR_VAL(2.0)*a)*v[t7][t9] - c*v[t7+1][t9]-d*q[t7][t9-1])/(d*p[t7][t9-1]+e);;
        }
      }
    }
    lbv=max(1,32*t4);
    ubv=min(_PB_N-2,32*t4+31);
#pragma ivdep
#pragma vector always
    for (t7=lbv;t7<=ubv;t7++) {
      u[t7][_PB_N-1] = SCALAR_VAL(1.0);;
    }
    for (t6=ceild(-_PB_N-29,32);t6<=-1;t6++) {
      for (t7=max(1,32*t4);t7<=min(_PB_N-2,32*t4+31);t7++) {
        for (t9=max(32*t6,-_PB_N+2);t9<=32*t6+31;t9++) {
          u[t7][-t9] = p[t7][-t9] * u[t7][-t9+1] + q[t7][-t9];;
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
  int tsteps = TSTEPS;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(u, DATA_TYPE, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(v, DATA_TYPE, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(p, DATA_TYPE, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(q, DATA_TYPE, N, N, n, n);


  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(u));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_adi (tsteps, n, POLYBENCH_ARRAY(u), POLYBENCH_ARRAY(v), POLYBENCH_ARRAY(p), POLYBENCH_ARRAY(q));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(u)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(u);
  POLYBENCH_FREE_ARRAY(v);
  POLYBENCH_FREE_ARRAY(p);
  POLYBENCH_FREE_ARRAY(q);

  return 0;
}
