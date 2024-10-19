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

1. **Chunking**: The loop is divided into chunks of 32 elements (`chunk_size * t1` and `chunk_size * t2`). This helps in better cache utilization.
2. **Parallelization**: The outer loop is parallelized using `#pragma omp parallel for` to distribute the work across multiple threads.
3. **Vectorization**: The inner loop is vectorized using `#pragma ivdep` and `#pragma vector always` to ensure that the loop is optimized for SIMD architectures.
4. **Loop Distribution**: The loops are distributed across different levels to optimize cache usage and reduce the number of iterations.

This optimized code should provide better performance by leveraging parallelization, vectorization, and loop distribution techniques.*/

#include <omp.h>
#include <math.h>

#define min(a, b) ((a) < (b) ? (a) : (b))
#define floord(n, d) ((n) / (d))

int t1, t2, t3, t4;
int lb, ub, lbp, ubp, lb2, ub2;
register int lbv, ubv;
double chunk_size = 32.0;

lbp = 0;
ubp = floord(row_dim - 1, chunk_size);
#pragma omp parallel for private(lbv, ubv, t2, t3, t4)
for (int t1 = lbp; t1 <= ubp; t1++) {
    for (int t2 = 0; t2 <= floord(col_dim - 1, chunk_size); t2++) {
        for (int t3 = chunk_size * t1; t3 <= min(row_dim - 1, chunk_size * t1 + chunk_size - 1); t3++) {
            lbv = chunk_size * t2;
            ubv = min(col_dim - 1, chunk_size * t2 + chunk_size - 1);
#pragma ivdep
#pragma vector always
            for (int t4 = lbv; t4 <= ubv; t4++) {
                diff[t4][t3] = A[t4][t3] - B[t4][t3];
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