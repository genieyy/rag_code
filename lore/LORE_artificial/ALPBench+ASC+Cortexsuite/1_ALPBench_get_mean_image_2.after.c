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
1. **Reduction Variable (`temp_sum`)**: 
   - Introduced a temporary variable `temp_sum` to accumulate the sum of `images[j][i]` within the inner loop. This avoids repeatedly accessing and modifying `mean[0][i]` within the inner loop, which can be more efficient.

2. **Loop Order**:
   - The loop order remains the same, iterating over `i` first and then `j`. This is generally efficient for memory access patterns, assuming `images` is stored in row-major order.

3. **Division Outside Inner Loop**:
   - The division by `col_dim` is moved outside the inner loop, reducing the number of division operations from `col_dim` times to just once per row.

These optimizations are based on common techniques such as loop unrolling, reduction variables, and minimizing redundant operations, which are often used to improve performance in nested loops.*/

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