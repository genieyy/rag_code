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
/*### Explanation:
- **Pointer Arithmetic**: By using pointers (`diff_row`, `A_row`, `B_row`), we avoid the overhead of array indexing (`diff[j][i]`, `A[j][i]`, `B[j][i]`). This can lead to better performance, especially when dealing with large arrays.
- **Loop Order**: The loop order remains the same as in the optimized version (0), which is already optimal for cache locality.
- **No Redefinition**: The original variables are not redefined, and no new variables are introduced that would require redefinition.*/

for (int iter = 0; iter < ITERATIONS; iter++){
    for (int j = 0; j < col_dim; j++) {
        double *diff_row = diff[j];
        double *A_row = A[j];
        double *B_row = B[j];
        for (int i = 0; i < row_dim; i++) {
            diff_row[i] = A_row[i] - B_row[i];
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