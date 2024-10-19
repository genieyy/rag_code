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
#define N 40
#define M 50
#define Q 60
#define P 5
/* end parameters define */

/* start kernel func */
void NPB_bt8(double ****forcing, int n_, int m_, int q_, int p_) {
    int i, j, k, m;
    double dssp = 0.5; // Example value, replace with actual value

    double time_start = omp_get_wtime();
 #pragma scop
/*### Explanation:
1. **Pointer Arithmetic**: The pointer `forcing_ptr` is used to access the elements of the `forcing` array, which avoids the overhead of multi-dimensional array indexing.
2. **Constant Multiplication**: The multiplication by `-1.0` is replaced with a multiplication by a constant `neg_one`, which is precomputed outside the innermost loop. This reduces the number of floating-point operations inside the loop.
3. **Loop Unrolling**: Although not explicitly unrolled, the use of pointer arithmetic and constant multiplication helps the compiler to potentially optimize the loop further through automatic loop unrolling and other optimizations.

This version should be more efficient due to the reduced number of operations inside the innermost loop and the use of pointer arithmetic.*/

for (int iter = 0; iter < ITERATIONS; iter++){
    for (i = 1; i < n_ - 1; i++) {
        for (j = 1; j < m_ - 1; j++) {
            for (k = 1; k < q_ - 1; k++) {
                double *forcing_ptr = &forcing[i][j][k][0];
                double neg_one = -1.0;
                for (m = 0; m < p_; m++) {
                    forcing_ptr[m] *= neg_one;
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
    ARRAY_PREPARATION_4D(array_0, N+1, M+1, Q+1, P+1);

    NPB_bt8(array_0, N, M, Q, P);

    print_array_4d(array_0, N+1, M+1, Q+1, P+1);

    free_array_4d(array_0, N+1, M+1, Q+1, P+1);
    
    return 0;
}