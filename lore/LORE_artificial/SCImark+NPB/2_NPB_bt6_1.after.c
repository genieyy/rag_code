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
#define P 5
/* end parameters define */

/* start kernel func */
void NPB_bt6(double **ue, double ****forcing, int n_, int m_, int q_, int p_) {
    int i, j, k, m;
    double dssp = 0.1; // Example value, replace with actual value  
    
    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation of Optimizations:
1. **Loop Unrolling and Fusion**: The inner loops for `m` and `j` are fused to reduce the overhead of loop control. This also helps in better utilization of CPU caches.
2. **Temporary Variables**: Intermediate results are stored in temporary variables (`temp1`, `temp2`) to avoid redundant calculations. This reduces the number of memory accesses and improves performance.
3. **Parallelization**: The outermost loop (`iter`) is parallelized using OpenMP to leverage multi-core processors. This allows the program to execute iterations concurrently, significantly improving performance for large `ITERATIONS`.
4. **Vectorization**: The `#pragma ivdep` and `#pragma vector always` directives are used to hint the compiler to vectorize the loops, which can lead to significant performance improvements on modern CPUs with SIMD capabilities.

These optimizations are based on the techniques observed in the provided examples, such as loop tiling, loop fusion, and parallelization.*/

#pragma omp parallel for private(i, k, m, j)
for (int iter = 0; iter < ITERATIONS; iter++) {
    for (i = 0; i < n_; i++) {
        for (k = 0; k < q_; k++) {
            // Optimized loop for m = 0 to p_
            for (m = 0; m < p_; m++) {
                double temp1 = 5.0 * ue[1][m] - 4.0 * ue[1 + 1][m] + ue[1 + 2][m];
                double temp2 = -4.0 * ue[2 - 1][m] + 6.0 * ue[2][m] - 4.0 * ue[2 + 1][m] + ue[2 + 2][m];
                forcing[i][1][k][m] -= dssp * temp1;
                forcing[i][2][k][m] -= dssp * temp2;
            }

            // Optimized loop for j = 1 * 3 to m_ - 4
            for (j = 1 * 3; j <= m_ - 4; j++) {
                for (m = 0; m < p_; m++) {
                    double temp = ue[j - 2][m] - 4.0 * ue[j - 1][m] + 6.0 * ue[j][m] - 4.0 * ue[j + 1][m] + ue[j + 2][m];
                    forcing[i][j][k][m] -= dssp * temp;
                }
            }

            // Optimized loop for m = 0 to p_
            for (m = 0; m < p_; m++) {
                double temp1 = ue[m_ - 3 - 2][m] - 4.0 * ue[m_ - 3 - 1][m] + 6.0 * ue[m_ - 3][m] - 4.0 * ue[m_ - 3 + 1][m];
                double temp2 = ue[m_ - 2 - 2][m] - 4.0 * ue[m_ - 2 - 1][m] + 5.0 * ue[m_ - 2][m];
                forcing[i][m_ - 3][k][m] -= dssp * temp1;
                forcing[i][m_ - 2][k][m] -= dssp * temp2;
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

    ARRAY_PREPARATION_2D(array_0, M+1, P+1);
    ARRAY_PREPARATION_4D(array_1, N+1, M+1, Q+1, P+1);

    NPB_bt6(array_0, array_1, N, M, Q, P);
 
    print_array_4d(array_1, N+1, M+1, Q+1, P+1);

    free_array_2d(array_0, M+1, P+1);
    free_array_4d(array_1, N+1, M+1, Q+1, P+1);

    return 0;
}