/*
ALPBench/Face_Rec/src/csuSubspaceFisher.c
line 89 - line 106
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
void ALPBench_subtractClassMean(double **mean, double **classmatrix, int classmatrix_row_dim, int classEnd) {
    int i, j;

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation of Optimizations:
1. **Combined Mean Calculation and Subtraction**: The calculation of the mean and the subtraction of the mean from `classmatrix` are combined into a single loop for each row. This reduces the number of iterations over the rows of `classmatrix` from two to one, improving performance.
2. **Temporary Variable for Mean**: Introduced a temporary variable `mean_val` to store the mean value for each row. This reduces the number of accesses to the `mean` array within the inner loop.
3. **Reduction Variable**: Introduced a local variable `sum` to accumulate the sum of elements in `classmatrix` for each row. This avoids repeatedly accessing and modifying `mean[0][i]` within the inner loop.
4. **Loop Invariant Code Motion**: Moved the calculation of `mean[0][i]` outside the inner loop where it is invariant, reducing redundant calculations.*/

for (int iter = 0; iter < ITERATIONS; iter++) {
    for (int i = 0; i < classmatrix_row_dim; i++) {
        double sum = 0.0;
        for (int j = 0; j <= classEnd; j++) {
            sum += classmatrix[j][i];
        }
        double mean_val = sum / (classEnd + 1);
        mean[0][i] = mean_val;
        for (int j = 0; j <= classEnd; j++) {
            classmatrix[j][i] -= mean_val;
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
    ARRAY_PREPARATION_2D(array_0, 1, M+1);
    ARRAY_PREPARATION_2D(array_1, N+1, M+1);

    ALPBench_subtractClassMean(array_0, array_1, M, N);
    
    print_array_2d(array_1, N+1, M+1);

    free_array_2d(array_0, 1, M+1);
    free_array_2d(array_1, N+1, M+1);
    
    return 0;
}