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
#define L 40
#define P 5
/* end parameters define */

/* start kernel func */
void NPB_bt14(double ****rhs, double ****forcing, int nz, int mz, int q, int p)
{
    int i, j, k, m;

    double time_start = omp_get_wtime();
#pragma scop
/*### Analysis of Loop Transformations

The provided optimized code demonstrates several loop transformation techniques:

1. **Loop Tiling/Blocking**: The original nested loops are transformed into a tiled structure, which helps in improving cache locality. This is evident from the introduction of new variables like `t1`, `t3`, `t4`, etc., which are used to create smaller blocks of iterations.

2. **Parallelization**: The use of `#pragma omp parallel for` indicates that the loop iterations are being parallelized across multiple threads, which can significantly improve performance on multi-core processors.

3. **Vectorization**: The use of `#pragma ivdep` and `#pragma vector always` suggests that the loop is being vectorized, which allows the compiler to generate SIMD (Single Instruction, Multiple Data) instructions, improving performance on modern CPUs with vector units.

4. **Loop Fusion/Fission**: The original loops are split and fused in a way that reduces the number of loop nests and combines operations that can be performed together, reducing the overhead of loop control.

### Applying These Techniques to the New Code

To optimize the given code, we can apply similar techniques:

1. **Loop Tiling**: We can tile the loops to improve cache locality.
2. **Parallelization**: We can parallelize the outermost loop to leverage multi-core processors.
3. **Vectorization**: We can vectorize the innermost loop to take advantage of SIMD instructions.

Here is the optimized code:



### Explanation of the Optimized Code

1. **Loop Tiling**: The loops are tiled using `tile_size_i`, `tile_size_j`, `tile_size_k`, and `tile_size_m` to create smaller blocks of iterations. This improves cache locality by ensuring that the data accessed within each tile fits into the cache.

2. **Parallelization**: The outermost loop (`iter`) is parallelized using `#pragma omp parallel for collapse(4)`. The `collapse(4)` directive ensures that the four nested loops (`ti`, `tj`, `tk`, `tm`) are parallelized together, which can improve load balancing across threads.

3. **Vectorization**: The innermost loop (`m`) is vectorized using `#pragma ivdep` and `#pragma vector always`. This allows the compiler to generate SIMD instructions for the assignment operation, improving performance on CPUs with vector units.

By applying these techniques, the optimized code should achieve better performance in terms of cache utilization, parallel execution, and vectorization.*/

int tile_size_i = 32;
int tile_size_j = 32;
int tile_size_k = 32;
int tile_size_m = 32;

for (int iter = 0; iter < ITERATIONS; iter++) {
    #pragma omp parallel for collapse(4)
    for (int ti = 0; ti < nz; ti += tile_size_i) {
        for (int tj = 0; tj < mz; tj += tile_size_j) {
            for (int tk = 0; tk < q; tk += tile_size_k) {
                for (int tm = 0; tm < p; tm += tile_size_m) {
                    for (int i = ti; i < min(ti + tile_size_i, nz); i++) {
                        for (int j = tj; j < min(tj + tile_size_j, mz); j++) {
                            for (int k = tk; k < min(tk + tile_size_k, q); k++) {
                                #pragma ivdep
                                #pragma vector always
                                for (int m = tm; m < min(tm + tile_size_m, p); m++) {
                                    rhs[i][j][k][m] = forcing[i][j][k][m];
                                }
                            }
                        }
                    }
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

int main(int argc, char *argv[])
{
    ARRAY_PREPARATION_4D(array_0, N + 1, M + 1, L + 1, P + 1);
    ARRAY_PREPARATION_4D(array_1, N + 1, M + 1, L + 1, P + 1);

    NPB_bt14(array_0, array_1, N, M, L, P);

    print_array_4d(array_0, N + 1, M + 1, L + 1, P + 1);

    free_array_4d(array_0, N + 1, M + 1, L + 1, P + 1);
    free_array_4d(array_1, N + 1, M + 1, L + 1, P + 1);

    return 0;
}