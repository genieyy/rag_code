#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <polybench.h>
#include "heat-3d.h"
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
/* heat-3d.c: this file is part of PolyBench/C */


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
		 DATA_TYPE POLYBENCH_3D(A,N,N,N,n,n,n),
		 DATA_TYPE POLYBENCH_3D(B,N,N,N,n,n,n))
{
  int i, j, k;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      for (k = 0; k < n; k++)
        A[i][j][k] = B[i][j][k] = (DATA_TYPE) (i + j + (n-k))* 10 / (n);
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_3D(A,N,N,N,n,n,n))

{
  int i, j, k;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("A");
  DATA_TYPE tmp_A = 0;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      for (k = 0; k < n; k++) {
        if (DUMP) {
           if ((i * n * n + j * n + k) % 20 == 0) fprintf(POLYBENCH_DUMP_TARGET, "\n");
           fprintf(POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, A[i][j][k]);
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
static
void kernel_heat_3d(int tsteps,
		      int n,
		      DATA_TYPE POLYBENCH_3D(A,N,N,N,n,n,n),
		      DATA_TYPE POLYBENCH_3D(B,N,N,N,n,n,n))
{
  int t, i, j, k;
  int t1, t2, t3, t4, t5, t6, t7, t8;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
for (t1=0;t1<=floord(3*_PB_TSTEPS+_PB_N-1,32);t1++) {
  lbp=max(ceild(2*t1,3),ceild(32*t1-_PB_TSTEPS,32));
  ubp=min(min(floord(2*_PB_TSTEPS+_PB_N-1,32),floord(64*t1+_PB_N+61,96)),t1);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6,t7,t8)
  for (t2=lbp;t2<=ubp;t2++) {
    for (t3=max(ceild(32*t2-_PB_N-28,32),2*t1-2*t2);t3<=min(min(floord(2*_PB_TSTEPS+_PB_N-1,32),floord(32*t2+_PB_N+28,32)),floord(64*t1-64*t2+_PB_N+61,32));t3++) {
      for (t4=max(max(ceild(32*t2-_PB_N-28,32),ceild(32*t3-_PB_N-28,32)),2*t1-2*t2);t4<=min(min(min(floord(2*_PB_TSTEPS+_PB_N-1,32),floord(32*t2+_PB_N+28,32)),floord(32*t3+_PB_N+28,32)),floord(64*t1-64*t2+_PB_N+61,32));t4++) {
        if ((_PB_N == 121) && (t1 <= floord(2*t2+t4-4,2)) && (t2 <= t4-1) && (t3 <= t4-1)) {
          for (t6=max(32*t2,32*t4-118);t6<=32*t2+31;t6++) {
            for (t7=max(32*t3,32*t4-118);t7<=32*t3+31;t7++) {
              A[(-32*t4+t6+119)][(-32*t4+t7+119)][119] = SCALAR_VAL(0.125) * (B[(-32*t4+t6+119)+1][(-32*t4+t7+119)][119] - SCALAR_VAL(2.0) * B[(-32*t4+t6+119)][(-32*t4+t7+119)][119] + B[(-32*t4+t6+119)-1][(-32*t4+t7+119)][119]) + SCALAR_VAL(0.125) * (B[(-32*t4+t6+119)][(-32*t4+t7+119)+1][119] - SCALAR_VAL(2.0) * B[(-32*t4+t6+119)][(-32*t4+t7+119)][119] + B[(-32*t4+t6+119)][(-32*t4+t7+119)-1][119]) + SCALAR_VAL(0.125) * (B[(-32*t4+t6+119)][(-32*t4+t7+119)][119 +1] - SCALAR_VAL(2.0) * B[(-32*t4+t6+119)][(-32*t4+t7+119)][119] + B[(-32*t4+t6+119)][(-32*t4+t7+119)][119 -1]) + B[(-32*t4+t6+119)][(-32*t4+t7+119)][119];;
            }
          }
        }
        if ((_PB_N == 121) && (t1 <= floord(2*t2+t3-4,2)) && (t2 <= t3-1) && (t3 >= t4)) {
          for (t6=max(32*t2,32*t3-118);t6<=32*t2+31;t6++) {
            lbv=max(32*t4,32*t3-118);
            ubv=min(32*t3,32*t4+31);
#pragma ivdep
#pragma vector always
            for (t8=lbv;t8<=ubv;t8++) {
              A[(-32*t3+t6+119)][119][(-32*t3+t8+119)] = SCALAR_VAL(0.125) * (B[(-32*t3+t6+119)+1][119][(-32*t3+t8+119)] - SCALAR_VAL(2.0) * B[(-32*t3+t6+119)][119][(-32*t3+t8+119)] + B[(-32*t3+t6+119)-1][119][(-32*t3+t8+119)]) + SCALAR_VAL(0.125) * (B[(-32*t3+t6+119)][119 +1][(-32*t3+t8+119)] - SCALAR_VAL(2.0) * B[(-32*t3+t6+119)][119][(-32*t3+t8+119)] + B[(-32*t3+t6+119)][119 -1][(-32*t3+t8+119)]) + SCALAR_VAL(0.125) * (B[(-32*t3+t6+119)][119][(-32*t3+t8+119)+1] - SCALAR_VAL(2.0) * B[(-32*t3+t6+119)][119][(-32*t3+t8+119)] + B[(-32*t3+t6+119)][119][(-32*t3+t8+119)-1]) + B[(-32*t3+t6+119)][119][(-32*t3+t8+119)];;
            }
          }
        }
        if ((_PB_N == 121) && (t1 <= floord(3*t2-4,2)) && (t2 >= max(t3,t4))) {
          for (t7=max(32*t3,32*t2-118);t7<=min(32*t2,32*t3+31);t7++) {
            lbv=max(32*t4,32*t2-118);
            ubv=min(32*t2,32*t4+31);
#pragma ivdep
#pragma vector always
            for (t8=lbv;t8<=ubv;t8++) {
              A[119][(-32*t2+t7+119)][(-32*t2+t8+119)] = SCALAR_VAL(0.125) * (B[119 +1][(-32*t2+t7+119)][(-32*t2+t8+119)] - SCALAR_VAL(2.0) * B[119][(-32*t2+t7+119)][(-32*t2+t8+119)] + B[119 -1][(-32*t2+t7+119)][(-32*t2+t8+119)]) + SCALAR_VAL(0.125) * (B[119][(-32*t2+t7+119)+1][(-32*t2+t8+119)] - SCALAR_VAL(2.0) * B[119][(-32*t2+t7+119)][(-32*t2+t8+119)] + B[119][(-32*t2+t7+119)-1][(-32*t2+t8+119)]) + SCALAR_VAL(0.125) * (B[119][(-32*t2+t7+119)][(-32*t2+t8+119)+1] - SCALAR_VAL(2.0) * B[119][(-32*t2+t7+119)][(-32*t2+t8+119)] + B[119][(-32*t2+t7+119)][(-32*t2+t8+119)-1]) + B[119][(-32*t2+t7+119)][(-32*t2+t8+119)];;
            }
          }
        }
        for (t5=max(max(max(max(1,ceild(32*t2-_PB_N+2,2)),ceild(32*t3-_PB_N+2,2)),ceild(32*t4-_PB_N+2,2)),32*t1-32*t2);t5<=min(min(min(min(_PB_TSTEPS,16*t2+14),16*t3+14),16*t4+14),32*t1-32*t2+31);t5++) {
          if (t2 <= floord(t5,16)) {
            for (t7=max(32*t3,2*t5+1);t7<=min(32*t3+31,2*t5+_PB_N-2);t7++) {
              lbv=max(32*t4,2*t5+1);
              ubv=min(32*t4+31,2*t5+_PB_N-2);
#pragma ivdep
#pragma vector always
              for (t8=lbv;t8<=ubv;t8++) {
                B[1][(-2*t5+t7)][(-2*t5+t8)] = SCALAR_VAL(0.125) * (A[1 +1][(-2*t5+t7)][(-2*t5+t8)] - SCALAR_VAL(2.0) * A[1][(-2*t5+t7)][(-2*t5+t8)] + A[1 -1][(-2*t5+t7)][(-2*t5+t8)]) + SCALAR_VAL(0.125) * (A[1][(-2*t5+t7)+1][(-2*t5+t8)] - SCALAR_VAL(2.0) * A[1][(-2*t5+t7)][(-2*t5+t8)] + A[1][(-2*t5+t7)-1][(-2*t5+t8)]) + SCALAR_VAL(0.125) * (A[1][(-2*t5+t7)][(-2*t5+t8)+1] - SCALAR_VAL(2.0) * A[1][(-2*t5+t7)][(-2*t5+t8)] + A[1][(-2*t5+t7)][(-2*t5+t8)-1]) + A[1][(-2*t5+t7)][(-2*t5+t8)];;
              }
            }
          }
          for (t6=max(32*t2,2*t5+2);t6<=min(32*t2+31,2*t5+_PB_N-2);t6++) {
            if (t3 <= floord(t5,16)) {
              lbv=max(32*t4,2*t5+1);
              ubv=min(32*t4+31,2*t5+_PB_N-2);
#pragma ivdep
#pragma vector always
              for (t8=lbv;t8<=ubv;t8++) {
                B[(-2*t5+t6)][1][(-2*t5+t8)] = SCALAR_VAL(0.125) * (A[(-2*t5+t6)+1][1][(-2*t5+t8)] - SCALAR_VAL(2.0) * A[(-2*t5+t6)][1][(-2*t5+t8)] + A[(-2*t5+t6)-1][1][(-2*t5+t8)]) + SCALAR_VAL(0.125) * (A[(-2*t5+t6)][1 +1][(-2*t5+t8)] - SCALAR_VAL(2.0) * A[(-2*t5+t6)][1][(-2*t5+t8)] + A[(-2*t5+t6)][1 -1][(-2*t5+t8)]) + SCALAR_VAL(0.125) * (A[(-2*t5+t6)][1][(-2*t5+t8)+1] - SCALAR_VAL(2.0) * A[(-2*t5+t6)][1][(-2*t5+t8)] + A[(-2*t5+t6)][1][(-2*t5+t8)-1]) + A[(-2*t5+t6)][1][(-2*t5+t8)];;
              }
            }
            for (t7=max(32*t3,2*t5+2);t7<=min(32*t3+31,2*t5+_PB_N-2);t7++) {
              if (t4 <= floord(t5,16)) {
                B[(-2*t5+t6)][(-2*t5+t7)][1] = SCALAR_VAL(0.125) * (A[(-2*t5+t6)+1][(-2*t5+t7)][1] - SCALAR_VAL(2.0) * A[(-2*t5+t6)][(-2*t5+t7)][1] + A[(-2*t5+t6)-1][(-2*t5+t7)][1]) + SCALAR_VAL(0.125) * (A[(-2*t5+t6)][(-2*t5+t7)+1][1] - SCALAR_VAL(2.0) * A[(-2*t5+t6)][(-2*t5+t7)][1] + A[(-2*t5+t6)][(-2*t5+t7)-1][1]) + SCALAR_VAL(0.125) * (A[(-2*t5+t6)][(-2*t5+t7)][1 +1] - SCALAR_VAL(2.0) * A[(-2*t5+t6)][(-2*t5+t7)][1] + A[(-2*t5+t6)][(-2*t5+t7)][1 -1]) + A[(-2*t5+t6)][(-2*t5+t7)][1];;
              }
              lbv=max(32*t4,2*t5+2);
              ubv=min(32*t4+31,2*t5+_PB_N-2);
#pragma ivdep
#pragma vector always
              for (t8=lbv;t8<=ubv;t8++) {
                B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = SCALAR_VAL(0.125) * (A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] - SCALAR_VAL(2.0) * A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)]) + SCALAR_VAL(0.125) * (A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] - SCALAR_VAL(2.0) * A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)]) + SCALAR_VAL(0.125) * (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] - SCALAR_VAL(2.0) * A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1]) + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)];;
                A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = SCALAR_VAL(0.125) * (B[(-2*t5+t6-1)+1][(-2*t5+t7-1)][(-2*t5+t8-1)] - SCALAR_VAL(2.0) * B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] + B[(-2*t5+t6-1)-1][(-2*t5+t7-1)][(-2*t5+t8-1)]) + SCALAR_VAL(0.125) * (B[(-2*t5+t6-1)][(-2*t5+t7-1)+1][(-2*t5+t8-1)] - SCALAR_VAL(2.0) * B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] + B[(-2*t5+t6-1)][(-2*t5+t7-1)-1][(-2*t5+t8-1)]) + SCALAR_VAL(0.125) * (B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)+1] - SCALAR_VAL(2.0) * B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] + B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)-1]) + B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
              }
              if (t4 >= ceild(2*t5+_PB_N-32,32)) {
                A[(-2*t5+t6-1)][(-2*t5+t7-1)][(_PB_N-2)] = SCALAR_VAL(0.125) * (B[(-2*t5+t6-1)+1][(-2*t5+t7-1)][(_PB_N-2)] - SCALAR_VAL(2.0) * B[(-2*t5+t6-1)][(-2*t5+t7-1)][(_PB_N-2)] + B[(-2*t5+t6-1)-1][(-2*t5+t7-1)][(_PB_N-2)]) + SCALAR_VAL(0.125) * (B[(-2*t5+t6-1)][(-2*t5+t7-1)+1][(_PB_N-2)] - SCALAR_VAL(2.0) * B[(-2*t5+t6-1)][(-2*t5+t7-1)][(_PB_N-2)] + B[(-2*t5+t6-1)][(-2*t5+t7-1)-1][(_PB_N-2)]) + SCALAR_VAL(0.125) * (B[(-2*t5+t6-1)][(-2*t5+t7-1)][(_PB_N-2)+1] - SCALAR_VAL(2.0) * B[(-2*t5+t6-1)][(-2*t5+t7-1)][(_PB_N-2)] + B[(-2*t5+t6-1)][(-2*t5+t7-1)][(_PB_N-2)-1]) + B[(-2*t5+t6-1)][(-2*t5+t7-1)][(_PB_N-2)];;
              }
            }
            if (t3 >= ceild(2*t5+_PB_N-32,32)) {
              lbv=max(32*t4,2*t5+2);
              ubv=min(32*t4+31,2*t5+_PB_N-1);
#pragma ivdep
#pragma vector always
              for (t8=lbv;t8<=ubv;t8++) {
                A[(-2*t5+t6-1)][(_PB_N-2)][(-2*t5+t8-1)] = SCALAR_VAL(0.125) * (B[(-2*t5+t6-1)+1][(_PB_N-2)][(-2*t5+t8-1)] - SCALAR_VAL(2.0) * B[(-2*t5+t6-1)][(_PB_N-2)][(-2*t5+t8-1)] + B[(-2*t5+t6-1)-1][(_PB_N-2)][(-2*t5+t8-1)]) + SCALAR_VAL(0.125) * (B[(-2*t5+t6-1)][(_PB_N-2)+1][(-2*t5+t8-1)] - SCALAR_VAL(2.0) * B[(-2*t5+t6-1)][(_PB_N-2)][(-2*t5+t8-1)] + B[(-2*t5+t6-1)][(_PB_N-2)-1][(-2*t5+t8-1)]) + SCALAR_VAL(0.125) * (B[(-2*t5+t6-1)][(_PB_N-2)][(-2*t5+t8-1)+1] - SCALAR_VAL(2.0) * B[(-2*t5+t6-1)][(_PB_N-2)][(-2*t5+t8-1)] + B[(-2*t5+t6-1)][(_PB_N-2)][(-2*t5+t8-1)-1]) + B[(-2*t5+t6-1)][(_PB_N-2)][(-2*t5+t8-1)];;
              }
            }
          }
          if (t2 >= ceild(2*t5+_PB_N-32,32)) {
            for (t7=max(32*t3,2*t5+2);t7<=min(32*t3+31,2*t5+_PB_N-1);t7++) {
              lbv=max(32*t4,2*t5+2);
              ubv=min(32*t4+31,2*t5+_PB_N-1);
#pragma ivdep
#pragma vector always
              for (t8=lbv;t8<=ubv;t8++) {
                A[(_PB_N-2)][(-2*t5+t7-1)][(-2*t5+t8-1)] = SCALAR_VAL(0.125) * (B[(_PB_N-2)+1][(-2*t5+t7-1)][(-2*t5+t8-1)] - SCALAR_VAL(2.0) * B[(_PB_N-2)][(-2*t5+t7-1)][(-2*t5+t8-1)] + B[(_PB_N-2)-1][(-2*t5+t7-1)][(-2*t5+t8-1)]) + SCALAR_VAL(0.125) * (B[(_PB_N-2)][(-2*t5+t7-1)+1][(-2*t5+t8-1)] - SCALAR_VAL(2.0) * B[(_PB_N-2)][(-2*t5+t7-1)][(-2*t5+t8-1)] + B[(_PB_N-2)][(-2*t5+t7-1)-1][(-2*t5+t8-1)]) + SCALAR_VAL(0.125) * (B[(_PB_N-2)][(-2*t5+t7-1)][(-2*t5+t8-1)+1] - SCALAR_VAL(2.0) * B[(_PB_N-2)][(-2*t5+t7-1)][(-2*t5+t8-1)] + B[(_PB_N-2)][(-2*t5+t7-1)][(-2*t5+t8-1)-1]) + B[(_PB_N-2)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
              }
            }
          }
        }
        if ((t1 >= ceild(3*t2-1,2)) && (t2 <= min(min(floord(_PB_TSTEPS-15,16),t3-1),t4-1))) {
          for (t7=32*t3;t7<=min(32*t3+31,32*t2+_PB_N+28);t7++) {
            lbv=32*t4;
            ubv=min(32*t4+31,32*t2+_PB_N+28);
#pragma ivdep
#pragma vector always
            for (t8=lbv;t8<=ubv;t8++) {
              B[1][(-32*t2+t7-30)][(-32*t2+t8-30)] = SCALAR_VAL(0.125) * (A[1 +1][(-32*t2+t7-30)][(-32*t2+t8-30)] - SCALAR_VAL(2.0) * A[1][(-32*t2+t7-30)][(-32*t2+t8-30)] + A[1 -1][(-32*t2+t7-30)][(-32*t2+t8-30)]) + SCALAR_VAL(0.125) * (A[1][(-32*t2+t7-30)+1][(-32*t2+t8-30)] - SCALAR_VAL(2.0) * A[1][(-32*t2+t7-30)][(-32*t2+t8-30)] + A[1][(-32*t2+t7-30)-1][(-32*t2+t8-30)]) + SCALAR_VAL(0.125) * (A[1][(-32*t2+t7-30)][(-32*t2+t8-30)+1] - SCALAR_VAL(2.0) * A[1][(-32*t2+t7-30)][(-32*t2+t8-30)] + A[1][(-32*t2+t7-30)][(-32*t2+t8-30)-1]) + A[1][(-32*t2+t7-30)][(-32*t2+t8-30)];;
            }
          }
        }
        if ((t1 >= ceild(2*t2+t3-1,2)) && (t2 >= t3) && (t3 <= min(floord(_PB_TSTEPS-15,16),t4-1))) {
          for (t6=max(32*t2,32*t3+31);t6<=min(32*t2+31,32*t3+_PB_N+28);t6++) {
            lbv=32*t4;
            ubv=min(32*t4+31,32*t3+_PB_N+28);
#pragma ivdep
#pragma vector always
            for (t8=lbv;t8<=ubv;t8++) {
              B[(-32*t3+t6-30)][1][(-32*t3+t8-30)] = SCALAR_VAL(0.125) * (A[(-32*t3+t6-30)+1][1][(-32*t3+t8-30)] - SCALAR_VAL(2.0) * A[(-32*t3+t6-30)][1][(-32*t3+t8-30)] + A[(-32*t3+t6-30)-1][1][(-32*t3+t8-30)]) + SCALAR_VAL(0.125) * (A[(-32*t3+t6-30)][1 +1][(-32*t3+t8-30)] - SCALAR_VAL(2.0) * A[(-32*t3+t6-30)][1][(-32*t3+t8-30)] + A[(-32*t3+t6-30)][1 -1][(-32*t3+t8-30)]) + SCALAR_VAL(0.125) * (A[(-32*t3+t6-30)][1][(-32*t3+t8-30)+1] - SCALAR_VAL(2.0) * A[(-32*t3+t6-30)][1][(-32*t3+t8-30)] + A[(-32*t3+t6-30)][1][(-32*t3+t8-30)-1]) + A[(-32*t3+t6-30)][1][(-32*t3+t8-30)];;
            }
          }
        }
        if ((t1 >= ceild(2*t2+t4-1,2)) && (t2 >= t4) && (t3 >= t4) && (t4 <= floord(_PB_TSTEPS-15,16))) {
          for (t6=max(32*t2,32*t4+31);t6<=min(32*t2+31,32*t4+_PB_N+28);t6++) {
            for (t7=max(32*t3,32*t4+31);t7<=min(32*t3+31,32*t4+_PB_N+28);t7++) {
              B[(-32*t4+t6-30)][(-32*t4+t7-30)][1] = SCALAR_VAL(0.125) * (A[(-32*t4+t6-30)+1][(-32*t4+t7-30)][1] - SCALAR_VAL(2.0) * A[(-32*t4+t6-30)][(-32*t4+t7-30)][1] + A[(-32*t4+t6-30)-1][(-32*t4+t7-30)][1]) + SCALAR_VAL(0.125) * (A[(-32*t4+t6-30)][(-32*t4+t7-30)+1][1] - SCALAR_VAL(2.0) * A[(-32*t4+t6-30)][(-32*t4+t7-30)][1] + A[(-32*t4+t6-30)][(-32*t4+t7-30)-1][1]) + SCALAR_VAL(0.125) * (A[(-32*t4+t6-30)][(-32*t4+t7-30)][1 +1] - SCALAR_VAL(2.0) * A[(-32*t4+t6-30)][(-32*t4+t7-30)][1] + A[(-32*t4+t6-30)][(-32*t4+t7-30)][1 -1]) + A[(-32*t4+t6-30)][(-32*t4+t7-30)][1];;
            }
          }
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
  POLYBENCH_3D_ARRAY_DECL(A, DATA_TYPE, N, N, N, n, n, n);
  POLYBENCH_3D_ARRAY_DECL(B, DATA_TYPE, N, N, N, n, n, n);


  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(B));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_heat_3d (tsteps, n, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(B));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(A)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);

  return 0;
}
