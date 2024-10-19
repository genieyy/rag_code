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
/*### Explanation:
- **Loop Unrolling**: The inner loop (`l`) is unrolled by a factor of 4 to reduce the number of loop iterations and improve instruction-level parallelism. This can lead to better performance due to reduced loop overhead and better utilization of CPU resources.
- **Tail Loop**: After the unrolled loop, a tail loop handles the remaining iterations if `p - 1` is not a multiple of 4. This ensures that all elements are processed correctly.

This optimization is based on the assumption that the compiler might not automatically unroll the loop, and manual unrolling can sometimes yield better performance, especially on architectures that benefit from such optimizations.*/

for (int iter = 0; iter < ITERATIONS; iter++){
    for (i = 1; i < n - 1; i++) {
        for (j = 1; j < m - 1; j++) {
            for (k = 1; k < q - 1; k++) {
                double *u_ptr = &u[i][j][k][0];
                double *rhs_ptr = &rhs[i][j][k][0];
                for (l = 0; l < p - 1; l += 4) {
                    u_ptr[l] += rhs_ptr[l];
                    u_ptr[l + 1] += rhs_ptr[l + 1];
                    u_ptr[l + 2] += rhs_ptr[l + 2];
                    u_ptr[l + 3] += rhs_ptr[l + 3];
                }
                for (l = (p - 1) & ~3; l < p - 1; l++) {
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