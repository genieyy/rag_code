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
- **Partial Unrolling**: The innermost loop is partially unrolled by a factor of 4. This reduces the loop overhead and allows the CPU to perform more operations in parallel, potentially improving performance.
- **Pointer Arithmetic**: The use of `forcing_ptr` continues to reduce the overhead of array indexing and improve cache locality.
- **Loop Fission**: The loop is split into two parts: one for the unrolled iterations and one for the remaining iterations. This ensures that the unrolling does not cause any out-of-bounds accesses and handles the case where `p_` is not a multiple of 4.*/

/*### Explanation of Further Optimizations:
1. **Loop Unrolling**: While full unrolling of the innermost loop is not practical due to the unknown size of `p_`, we can partially unroll the loop by a small constant factor (e.g., 4) to reduce the loop overhead and improve instruction-level parallelism.

2. **Pointer Arithmetic**: We continue to use pointer arithmetic to reduce the overhead of array indexing and improve cache locality.

3. **Reduction in Indirection**: By using a pointer (`forcing_ptr`), we reduce the number of indirections, which can lead to better performance, especially if the array is large and the memory hierarchy is deep.

4. **Partial Unrolling**: The innermost loop is partially unrolled by a factor of 4, which can help the compiler generate more efficient code by reducing the number of loop control instructions and enabling better use of CPU registers.

These optimizations are based on the principles of reducing redundant computations, improving memory access patterns, and leveraging instruction-level parallelism, which are common techniques in loop transformation for performance improvement.*/

for (int iter = 0; iter < ITERATIONS; iter++){
    for (i = 1; i < n_ - 1; i++) {
        for (j = 1; j < m_ - 1; j++) {
            for (k = 1; k < q_ - 1; k++) {
                double *forcing_ptr = &forcing[i][j][k][0];
                int m;
                for (m = 0; m + 3 < p_; m += 4) {
                    forcing_ptr[m] = -1.0 * forcing_ptr[m];
                    forcing_ptr[m + 1] = -1.0 * forcing_ptr[m + 1];
                    forcing_ptr[m + 2] = -1.0 * forcing_ptr[m + 2];
                    forcing_ptr[m + 3] = -1.0 * forcing_ptr[m + 3];
                }
                for (; m < p_; m++) {
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