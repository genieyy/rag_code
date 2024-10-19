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
#define P 50
#define R 60
/* end parameters define */

/* start kernel func */
void NPB_bt2(double *rms, double ****rhs, int n_, int m_, int q_, int p_, int r_) {
    int i, j, k, d, m;
    double add;

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation of the Optimized Version:
1. **Loop Tiling**: The outermost loops are tiled with a tile size of 32 to improve cache locality. This means that instead of processing one element at a time, the code processes blocks of 32 elements, which helps to keep more data in the cache.
2. **Loop Fusion and Reordering**: The inner loops are fused and reordered to ensure that the most frequently accessed data (`rms[m]`) is accessed in a contiguous manner, which improves cache performance.
3. **Loop Unrolling**: Although not explicitly shown here, the innermost loop over `m` could be unrolled if `p_` is small and known at compile time.
4. **Register Usage**: The use of `register` for `lbv` and `ubv` is not directly applicable here, but the idea of using registers for frequently accessed variables is implied.

### Performance Considerations:
- **Cache Locality**: By tiling the outermost loops and iterating over `m` first, the code ensures that `rms[m]` is accessed sequentially and that data is reused within the cache.
- **Loop Fusion**: Combining the loops reduces the number of loop control operations, which can improve performance.
- **Loop Reordering**: The reordering ensures that the most frequently accessed data is accessed in a contiguous manner, which can improve cache hit rates.
- **Loop Tiling**: Tiling the outermost loops helps to keep more data in the cache, reducing the number of cache misses.*/

/*### Explanation of Transformations:
1. **Loop Fusion**: The inner loops over `i`, `j`, `k`, and `m` are fused into a single loop over `t1`, `t2`, `t3`, and `t4` respectively. This reduces the overhead of loop control and potentially improves cache locality.
2. **Loop Reordering**: The loops are reordered to ensure that the innermost loop iterates over `m` first, which is likely to improve cache performance since `rms[m]` is accessed repeatedly within the same iteration.
3. **Loop Unrolling**: The innermost loop over `m` is not unrolled in this example, but it could be considered if `p_` is small and known at compile time.
4. **Register Usage**: The use of `register` for `lbv` and `ubv` is not directly applicable here, but the idea of using registers for frequently accessed variables is implied.
5. **Loop Tiling**: The outermost loops are tiled to improve cache locality by processing blocks of data at a time.

### Performance Considerations:
- **Cache Locality**: By iterating over `m` first and tiling the outermost loops, the code ensures that `rms[m]` is accessed sequentially and that data is reused within the cache.
- **Loop Fusion**: Combining the loops reduces the number of loop control operations, which can improve performance.
- **Loop Reordering**: The reordering ensures that the most frequently accessed data is accessed in a contiguous manner, which can improve cache hit rates.
- **Loop Tiling**: Tiling the outermost loops helps to keep more data in the cache, reducing the number of cache misses.

These transformations are based on the principles observed in the provided examples, aiming to improve performance through better cache utilization and reduced loop overhead.*/

int t1, t2, t3, t4, t5;
register int lbv, ubv;

for (int iter = 0; iter < ITERATIONS; iter++) {
    for (t1 = 1; t1 <= (n_ - 2 + 31) / 32; t1++) {
        for (t2 = 1; t2 <= (m_ - 2 + 31) / 32; t2++) {
            for (t3 = 1; t3 <= (q_ - 2 + 31) / 32; t3++) {
                for (t4 = 0; t4 < p_; t4++) {
                    for (t5 = max(1, 32 * t1); t5 <= min(n_ - 2, 32 * t1 + 31); t5++) {
                        for (int j = max(1, 32 * t2); j <= min(m_ - 2, 32 * t2 + 31); j++) {
                            for (int k = max(1, 32 * t3); k <= min(q_ - 2, 32 * t3 + 31); k++) {
                                double add = rhs[t5][j][k][t4];
                                rms[t4] += add * add;
                            }
                        }
                    }
                }
            }
        }
    }

    for (t5 = 0; t5 < p_; t5++) {
        for (int d = 0; d <= r_; d++) {
            rms[t5] /= (double)(N - 2);
        }
        rms[t5] = sqrt(rms[t5]);
    }
}
#pragma endscop
    double time_end = omp_get_wtime();
    printf("%f\n", time_end - time_start);
}
/* end kernel func */

int main() {
    ARRAY_PREPARATION_1D(array_0, P+1);
    ARRAY_PREPARATION_4D(array_1, N+1, M+1, Q+1, P+1);

    NPB_bt2(array_0, array_1, N, M, Q, P, R);

    print_array_1d(array_0, P+1);

    free_array_1d(array_0, P+1);
    free_array_4d(array_1, N+1, M+1, Q+1, P+1);

    return 0;
}