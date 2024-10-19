/*
ALPBench/Face_Rec/src/csuCommonMatrix.c
line 128 - line 140
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
#define N 400
#define M 600
/* end parameters define */

/* start kernel func */
void ALPBench_makeZeroMatrix(double **A, int row_dim, int col_dim) {
    int i, j;

    double time_start = omp_get_wtime();
#pragma scop
/**/

int lbp, ubp;
register int lbv, ubv;
lbp = 0;
ubp = floord(ITERATIONS - 1, 32);
#pragma omp parallel for private(lbv, ubv)
for (int t1 = lbp; t1 <= ubp; t1++) {
    for (int t2 = 0; t2 <= floord(row_dim - 1, 32); t2++) {
        for (int t3 = 0; t3 <= floord(col_dim - 1, 32); t3++) {
            lbv = 32 * t1;
            ubv = min(ITERATIONS - 1, 32 * t1 + 31);
            for (int iter = lbv; iter <= ubv; iter++) {
                for (int i = 32 * t2; i <= min(row_dim - 1, 32 * t2 + 31); i++) {
                    #pragma ivdep
                    #pragma vector always
                    for (int j = 32 * t3; j <= min(col_dim - 1, 32 * t3 + 31); j++) {
                        A[j][i] = 0;
                    }
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
    ARRAY_PREPARATION_2D(array_0, M+1, N+1);

    ALPBench_makeZeroMatrix(array_0, N, M);
    
    print_array_2d(array_0, M+1, N+1);

    free_array_2d(array_0, M+1, N+1);

    return 0;
}