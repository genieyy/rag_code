/*
ALPBench/Face_Rec/src/csuCommonSubspace.c
line 761 - line 776
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
void ALPBench_get_mean_image(double **mean, double **images, int row_dim, int col_dim) {
    int i, j;

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation of Optimizations:
1. **Reduction Variable (`temp_sum`)**: Introduced a temporary variable `temp_sum` to accumulate the sum of `images[j][i]` within the inner loop. This reduces the number of memory accesses to `mean[0][i]` from `col_dim + 1` to just 2 (one for initialization and one for assignment after the loop).

2. **Loop Order**: The loop order remains the same, but the use of `temp_sum` ensures that the computation is more efficient by reducing redundant memory operations.

3. **Type Consistency**: The new variable `temp_sum` is of type `double` to match the type of the elements being summed and the final mean calculation.

These optimizations are inspired by the examples provided, where temporary variables were used to reduce redundant operations and improve performance.*/

double temp_sum;
for (int iter = 0; iter < ITERATIONS; iter++) {
    for (int i = 0; i < row_dim; i++) {
        temp_sum = 0.0;
        for (int j = 0; j < col_dim; j++) {
            temp_sum += images[j][i];
        }
        mean[0][i] = temp_sum / col_dim;
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

    ALPBench_get_mean_image(array_0, array_1, M, N);
    
    print_array_2d(array_0, 1, M+1);

    free_array_2d(array_0, 1, M+1);
    free_array_2d(array_1, N+1, M+1);
    
    return 0;
}