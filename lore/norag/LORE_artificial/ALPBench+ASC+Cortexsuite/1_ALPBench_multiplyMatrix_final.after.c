/*
ALPBench/Face_Rec/src/csuCommonMatrix.c
line 190 - line 211
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
#define N 60
#define M 70
#define L 80
/* end parameters define */

/* start kernel func */
void ALPBench_multiplyMatrix(double **P, double **A, double **B, int B_col_dim, int A_row_dim, int A_col_dim) {
    int i, j, k;                                   /* Variables used as indices */

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation:
- **Loop Unrolling**: The inner loop is unrolled by a factor of 4 to reduce the number of loop iterations and improve instruction-level parallelism. This can help the CPU to better utilize its execution units.
- **Reduced Memory Access**: By unrolling the loop, we reduce the number of times the loop condition and increment operations are performed, which can lead to fewer memory accesses and better cache utilization.
- **Vectorization Potential**: This version of the code is more amenable to vectorization by modern compilers, which can further optimize the performance by using SIMD (Single Instruction, Multiple Data) instructions.*/

for (int iter = 0; iter < ITERATIONS; iter++){
    /** Build the product matrix P */
    for ( j = 0; j < B_col_dim; j++) {
        for ( i = 0; i < A_row_dim; i++) {
            double temp = 0.0;
            int k = 0;
            for (; k <= A_col_dim - 4; k += 4) {
                temp += A[k][i] * B[j][k] +
                        A[k+1][i] * B[j][k+1] +
                        A[k+2][i] * B[j][k+2] +
                        A[k+3][i] * B[j][k+3];
            }
            for (; k < A_col_dim; k++) {
                temp += A[k][i] * B[j][k];
            }
            P[j][i] = temp;
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
    ARRAY_PREPARATION_2D(array_0, L+1, M+1);
    ARRAY_PREPARATION_2D(array_1, N+1, M+1);
    ARRAY_PREPARATION_2D(array_2, L+1, N+1);

    ALPBench_multiplyMatrix(array_0, array_1, array_2, L, M, N);
    
    print_array_2d(array_0, L+1, M+1);

    free_array_2d(array_0, L+1, M+1);
    free_array_2d(array_1, N+1, M+1);
    free_array_2d(array_2, L+1, N+1);
    
    return 0;
}