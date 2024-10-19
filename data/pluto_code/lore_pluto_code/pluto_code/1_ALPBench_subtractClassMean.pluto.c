#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lore.h"
#include <omp.h>
/*
ALPBench/Face_Rec/src/csuSubspaceFisher.c
line 89 - line 106
*/

#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y) ((x) > (y)? (x) : (y))
#define min(x,y) ((x) < (y)? (x) : (y))

/* start param define */ 
#define N 200
#define M 300
/* end parameters define */

/* start kernel func */
void ALPBench_subtractClassMean(double **mean, double **classmatrix, int classmatrix_row_dim, int classEnd) {
    int i, j;

    double time_start = omp_get_wtime();
  int t1, t2, t3, t4, t5, t6;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
lbp=0;
ubp=floord(classmatrix_row_dim-1,32);
#pragma omp parallel for private(lbv,ubv,t2,t3,t4,t5,t6)
for (t1=lbp;t1<=ubp;t1++) {
  for (t2=0;t2<=floord(ITERATIONS-1,32);t2++) {
    for (t3=32*t1;t3<=min(classmatrix_row_dim-1,32*t1+31);t3++) {
      for (t4=32*t2;t4<=min(ITERATIONS-1,32*t2+31);t4++) {
        mean[0][t3] = 0.0;;
        for (t6=0;t6<=classEnd;t6++) {
          mean[0][t3] += classmatrix[t6][t3];;
        }
        mean[0][t3] = mean[0][t3] / (classEnd + 1);;
        lbv=0;
        ubv=classEnd;
#pragma ivdep
#pragma vector always
        for (t6=lbv;t6<=ubv;t6++) {
          classmatrix[t6][t3] -= mean[0][t3];;
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
    ARRAY_PREPARATION_2D(array_0, 1, M+1);
    ARRAY_PREPARATION_2D(array_1, N+1, M+1);

    ALPBench_subtractClassMean(array_0, array_1, M, N);
    
    print_array_2d(array_1, N+1, M+1);

    free_array_2d(array_0, 1, M+1);
    free_array_2d(array_1, N+1, M+1);
    
    return 0;
}
