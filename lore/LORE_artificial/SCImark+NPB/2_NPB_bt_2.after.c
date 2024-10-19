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
#define N 20
#define M 30
#define Q 40
#define P 6
/* end parameters define */

/* start kernel func */
void NPB_bt(double ****u, double ****rhs, int n, int m, int q, int p) {
    int i, j, k, l;
    double time_start = omp_get_wtime();
#pragma scop
/*### Analysis and Transformation Methods Used:

1. **Loop Unrolling**: The original code has multiple nested loops, which can be unrolled to reduce loop overhead and improve instruction-level parallelism. However, in the provided optimized code, loop unrolling is not explicitly applied, but the structure is optimized for vectorization.

2. **Vectorization**: The use of `#pragma ivdep` and `#pragma vector always` in the original code indicates that the loops are being vectorized to take advantage of SIMD (Single Instruction, Multiple Data) instructions. This is done to process multiple elements of the array in parallel.

3. **Parallelization**: The use of `#pragma omp parallel for` in the optimized code indicates that the outer loop is being parallelized using OpenMP. This allows the iterations of the loop to be executed in parallel across multiple threads, which can significantly improve performance on multi-core processors.

4. **Pointer Arithmetic**: In the optimized code, pointer arithmetic is used to access the elements of the arrays `u` and `rhs`. This can sometimes be more efficient than using multi-dimensional array indexing, especially when dealing with large arrays.

5. **Loop Fusion**: The original code has multiple loops that are fused together in the optimized code. This reduces the overhead of loop control and can improve cache locality by processing related data in a single loop.

### Optimized Code Explanation:

- **Parallelization**: The outer loop over `iter` is parallelized using OpenMP, allowing multiple iterations to be executed concurrently.
- **Pointer Arithmetic**: The inner loops use pointer arithmetic to access the elements of the arrays `u` and `rhs`, which can be more efficient than using multi-dimensional array indexing.
- **Vectorization**: Although not explicitly shown, the inner loop over `l` is structured in a way that could be vectorized by the compiler, potentially improving performance.

This optimized code should perform better than the original, especially on multi-core processors, due to the parallelization and the use of pointer arithmetic.*/

#pragma omp parallel for private(i, j, k, l)
for (int iter = 0; iter < ITERATIONS; iter++) {
    for (i = 1; i < n - 1; i++) {
        for (j = 1; j < m - 1; j++) {
            for (k = 1; k < q - 1; k++) {
                double *u_ptr = &u[i][j][k][0];
                double *rhs_ptr = &rhs[i][j][k][0];
                for (l = 0; l < p - 1; l++) {
                    u_ptr[l] += rhs_ptr[l];
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

int main() {
    ARRAY_PREPARATION_4D(array_0, N + 1, M + 1, Q + 1, P + 1);
    ARRAY_PREPARATION_4D(array_1, N + 1, M + 1, Q + 1, P + 1);

    NPB_bt(array_0, array_1, N, M, Q, P);
    
    print_array_4d(array_0, N, M, Q, P);

    free_array_4d(array_0, N + 1, M + 1, Q + 1, P + 1);
    free_array_4d(array_1, N + 1, M + 1, Q + 1, P + 1);

    return 0;
}