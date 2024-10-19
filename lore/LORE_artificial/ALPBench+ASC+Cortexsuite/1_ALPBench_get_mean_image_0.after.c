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
/*### Explanation of the Optimization:
1. **Reduction of Array Accesses**: By introducing a temporary variable `temp_sum`, we reduce the number of times we access the `mean` array. This can improve performance by minimizing cache misses and reducing memory bandwidth usage.
2. **Loop Ordering**: The loop ordering remains the same, which is optimal for this type of computation where we accumulate values across columns for each row.
3. **Data Locality**: The inner loop accesses `images[j][i]`, which is a column-major access pattern. This is efficient if the `images` array is stored in column-major order, which is typical for many numerical libraries.

This optimization focuses on reducing redundant memory accesses and improving data locality, which are common techniques used in the provided examples.*/

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