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
1. **Loop Unrolling**: The innermost loop is not unrolled in this example, but it could be considered if `p_` is a small constant. However, since `p_` is not specified as a constant, unrolling might not be beneficial without further information.

2. **Pointer Arithmetic**: By using a pointer (`forcing_ptr`) to access the elements of the `forcing` array, we reduce the overhead of array indexing. This can improve cache performance and reduce the number of memory accesses.

3. **Loop Fusion**: The loops are already fused in the sense that they are nested, which is optimal for this type of access pattern.

4. **Loop Invariant Code Motion**: The multiplication by `-1.0` is loop-invariant and is applied inside the innermost loop, which is correct but could be optimized further if the loop structure allowed for more aggressive code motion.

5. **Cache Optimization**: By accessing memory sequentially using a pointer, we improve cache locality, which can lead to better performance, especially for large arrays.

This optimized code should perform better by reducing the overhead of array indexing and improving cache utilization.*/

for (int iter = 0; iter < ITERATIONS; iter++){
    for (i = 1; i < n_ - 1; i++) {
        for (j = 1; j < m_ - 1; j++) {
            for (k = 1; k < q_ - 1; k++) {
                double *forcing_ptr = &forcing[i][j][k][0];
                for (m = 0; m < p_; m++) {
                    forcing_ptr[m] = -1.0 * forcing_ptr[m];
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