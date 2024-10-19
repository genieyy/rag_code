/*
ALPBench/Face_Rec/src/csuCommonMatrix.c
line 486 - line 499
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
void ALPBench_addMatrixEquals(double **A, double **B, int row_dim, int col_dim)
{
    int i, j;

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation:
- **Pointer Arithmetic**: By using pointers (`A_row` and `B_row`), we avoid the overhead of array indexing (`A[j][i]` and `B[j][i]`). This can lead to better performance, especially when dealing with large arrays.
- **Loop Order**: The loop order remains the same as in the best-ranked version (0), which is optimal for cache locality.
- **Variable Type**: The new variables (`A_row` and `B_row`) are of type `double*`, which matches the type of the array elements.

This version should provide a further performance improvement by reducing the overhead of array indexing and leveraging pointer arithmetic.*/

for (int iter = 0; iter < ITERATIONS; iter++) {
    for (int j = 0; j < col_dim; j++) {
        double *A_row = A[j];
        double *B_row = B[j];
        for (int i = 0; i < row_dim; i++) {
            A_row[i] += B_row[i];
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
    ARRAY_PREPARATION_2D(array_0, N, M);
    ARRAY_PREPARATION_2D(array_1, N, M);

    ALPBench_addMatrixEquals(array_0, array_1, M, N);

    print_array_2d(array_0, N, M);

    free_array_2d(array_0, N, M);
    free_array_2d(array_1, N, M);
    
    return 0;
}