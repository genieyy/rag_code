/*
ALPBench/Face_Rec/src/csuCommonMatrix.c
line 449 - line 465
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lore.h"

#include <omp.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y) ((x) > (y)? (x) : (y))
#define min(x,y) ((x) < (y)? (x) : (y))

/* start param define */ 
#define N 200
#define M 300
/* end parameters define */

/* start kernel func */
void ALPBench_subtractMatrix(double **diff, double **A, double **B, int row_dim, int col_dim) {
    int i, j;

    double time_start = omp_get_wtime();
#pragma scop
/**/

int lbp, ubp;
register int lbv, ubv;
lbp = 0;
ubp = floord(row_dim - 1, 32);
#pragma omp parallel for private(lbv, ubv)
for (int iter = 0; iter < ITERATIONS; iter++) {
    for (int t1 = lbp; t1 <= ubp; t1++) {
        for (int t2 = 0; t2 <= floord(col_dim - 1, 32); t2++) {
            for (int t3 = 32 * t1; t3 <= min(row_dim - 1, 32 * t1 + 31); t3++) {
                lbv = 32 * t2;
                ubv = min(col_dim - 1, 32 * t2 + 31);
#pragma ivdep
#pragma vector always
                for (int t4 = lbv; t4 <= ubv; t4++) {
                    diff[t4][t3] = A[t4][t3] - B[t4][t3];
                }
            }
        }
    }
}
#pragma endscop
    double time_end = omp_get_wtime();
    printf("%f\n", time_end - time_start);
}
/* end kernel func */

int main(int argc, char *argv[])
{
    ARRAY_PREPARATION_2D(array_0, N+1, M+1);
    ARRAY_PREPARATION_2D(array_1, N+1, M+1);
    ARRAY_PREPARATION_2D(array_2, N+1, M+1);

    ALPBench_subtractMatrix(array_0, array_1, array_2, M, N);
    
    print_array_2d(array_0, N+1, M+1);

    free_array_2d(array_0, N+1, M+1);
    free_array_2d(array_1, N+1, M+1);
    free_array_2d(array_2, N+1, M+1);
    
    return 0;
}