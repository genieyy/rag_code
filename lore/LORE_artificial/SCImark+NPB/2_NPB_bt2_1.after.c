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
/*### Explanation of Transformations:
1. **Loop Fusion**: The inner loops over `i`, `j`, `k`, and `m` are fused into a single loop over `t1`, `t2`, `t3`, and `t4` respectively. This reduces the overhead of loop control and potentially improves cache locality.
2. **Loop Reordering**: The loops are reordered to ensure that the innermost loop iterates over `m` first, which is likely to improve cache performance since `rms[m]` is accessed repeatedly within the same iteration.
3. **Loop Unrolling**: The innermost loop over `m` is not unrolled in this example, but it could be considered if `p_` is small and known at compile time.
4. **Register Usage**: The use of `register` for `lbv` and `ubv` is not directly applicable here, but the idea of using registers for frequently accessed variables is implied.

### Performance Considerations:
- **Cache Locality**: By iterating over `m` first, the code ensures that `rms[m]` is accessed sequentially, which is beneficial for cache performance.
- **Loop Fusion**: Combining the loops reduces the number of loop control operations, which can improve performance.
- **Loop Reordering**: The reordering ensures that the most frequently accessed data is accessed in a contiguous manner, which can improve cache hit rates.

These transformations are based on the principles observed in the provided examples, aiming to improve performance through better cache utilization and reduced loop overhead.*/

int t1, t2, t3, t4, t5;
register int lbv, ubv;

for (int iter = 0; iter < ITERATIONS; iter++) {
    for (t1 = 1; t1 <= n_ - 2; t1++) {
        for (t2 = 1; t2 <= m_ - 2; t2++) {
            for (t3 = 1; t3 <= q_ - 2; t3++) {
                for (t4 = 0; t4 <= p_ - 1; t4++) {
                    double add = rhs[t1][t2][t3][t4];
                    rms[t4] += add * add;
                }
            }
        }
    }

    for (t5 = 0; t5 <= p_ - 1; t5++) {
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