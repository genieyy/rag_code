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
1. **Pointer Arithmetic**: By using pointer arithmetic, we avoid the repeated indexing calculations, which can be costly in deeply nested loops. This reduces the overhead of accessing elements in the arrays.
2. **Loop Unrolling**: Although not explicitly unrolled here, the use of pointers can sometimes lead to better optimization by the compiler, potentially enabling loop unrolling or other optimizations.
3. **Local Variable for Pointers**: The `u_ptr` and `rhs_ptr` variables are used to store the base addresses of the current slice of `u` and `rhs`, reducing the number of times the base address needs to be recalculated.

This transformation should improve the performance of the original code by reducing the overhead of array indexing and potentially enabling better optimizations by the compiler.*/

for (int iter = 0; iter < ITERATIONS; iter++){
    for (int i = 1; i < n - 1; i++) {
        for (int j = 1; j < m - 1; j++) {
            for (int k = 1; k < q - 1; k++) {
                double *u_ptr = &u[i][j][k][0];
                double *rhs_ptr = &rhs[i][j][k][0];
                for (int l = 0; l < p - 1; l++) {
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