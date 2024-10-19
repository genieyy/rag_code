#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lore.h"
#include <omp.h>
/*
ALPBench/Face_Rec/src/csuCommonMatrix.c
line 243 - line 269
*/

#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y) ((x) > (y)? (x) : (y))
#define min(x,y) ((x) < (y)? (x) : (y))

/* start param define */ 
#define N 60
#define M 70
#define L 80
/* end parameters define */

/* start kernel func */
void ALPBench_transposeMultiplyMatrixR(double **P, double **A, double **B, int B_row_dim, int A_row_dim, int A_col_dim) {
    int i, j, k;                                   /* Variables used as indices */

    double time_start = omp_get_wtime();
  int t1, t2, t3, t4, t5, t6, t7, t8;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
lbp=0;
ubp=floord(A_row_dim-1,32);
#pragma omp parallel for private(lbv,ubv,t2,t3,t4,t5,t6,t7,t8)
for (t1=lbp;t1<=ubp;t1++) {
  for (t2=0;t2<=floord(B_row_dim-1,32);t2++) {
    for (t3=0;t3<=floord(ITERATIONS-1,32);t3++) {
      for (t4=32*t1;t4<=min(A_row_dim-1,32*t1+31);t4++) {
        for (t5=32*t2;t5<=min(B_row_dim-1,32*t2+31);t5++) {
          for (t6=32*t3;t6<=min(ITERATIONS-1,32*t3+31);t6++) {
            P[t5][t4] = 0;;
            for (t7=0;t7<=A_col_dim-1;t7++) {
              P[t5][t4] += A[t7][t4] * B[t7][t5];;
            }
          }
        }
      }
    }
  }
}
/* End of CLooG code */
    double time_end = omp_get_wtime();
    printf("%f\n", time_end - time_start);
}
/* end kernel func */

int main(int argc, char *argv[])
{
    ARRAY_PREPARATION_2D(array_0, L+1, M+1);
    ARRAY_PREPARATION_2D(array_1, N+1, M+1);
    ARRAY_PREPARATION_2D(array_2, N+1, L+1);

    ALPBench_transposeMultiplyMatrixR(array_0, array_1, array_2, L, M, N);
    
    print_array_2d(array_0, L+1, M+1);

    free_array_2d(array_0, L+1, M+1);
    free_array_2d(array_1, N+1, M+1);
    free_array_2d(array_2, N+1, L+1);
    return 0;
}
