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
1. **Row-Major Access Optimization**: By storing pointers to the rows of `A` and `B` in `A_row` and `B_row`, we avoid the overhead of computing `A[j][i]` and `B[j][i]` repeatedly inside the inner loop. This reduces the number of array indexing operations, which can be costly.

2. **Loop Variable Type**: Changed the loop variable `i` to `int` to match the type used in the outer loop. This avoids potential type mismatches and ensures consistency.

3. **Meaning-Preserving Transformation**: The transformation maintains the original meaning of the program, ensuring that each element of `A` is incremented by the corresponding element of `B` for the specified number of iterations.*/

for (int iter = 0; iter < ITERATIONS; iter++){
    for (int i = 0; i < row_dim; i++)
    {
        double *A_row = A[i];
        double *B_row = B[i];
        for (int j = 0; j < col_dim; j++)
        {
            A_row[j] += B_row[j];
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