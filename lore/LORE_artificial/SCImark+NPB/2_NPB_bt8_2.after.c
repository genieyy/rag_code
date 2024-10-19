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
/*### Explanation of Optimizations:
1. **Loop Unrolling**: The innermost loop over `m` is not unrolled here, but it could be considered if `p_` is a small constant. However, since `p_` is not specified as a constant, unrolling might not be beneficial without further information.

2. **Pointer Arithmetic**: By using a pointer `forcing_ijk` to point to the start of the `forcing[i][j][k]` array, we reduce the number of array indexing operations. This can improve cache performance and reduce the overhead of multiple array indexing operations.

3. **Loop Order**: The loop order is maintained as it is, which is optimal for spatial locality in memory access patterns. Changing the order of the loops could potentially degrade performance due to less efficient cache usage.

4. **Constant Folding**: The multiplication by `-1.0` is kept as is, as it is a simple and efficient operation. No need to optimize further unless there are specific constraints or requirements.

These optimizations focus on reducing the overhead of array indexing and improving cache locality, which are common techniques to enhance the performance of nested loops in numerical computations.*/

for (int iter = 0; iter < ITERATIONS; iter++){
    for (int i = 1; i < n_ - 1; i++) {
        for (int j = 1; j < m_ - 1; j++) {
            for (int k = 1; k < q_ - 1; k++) {
                double *forcing_ijk = forcing[i][j][k];
                for (int m = 0; m < p_; m++) {
                    forcing_ijk[m] = -1.0 * forcing_ijk[m];
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