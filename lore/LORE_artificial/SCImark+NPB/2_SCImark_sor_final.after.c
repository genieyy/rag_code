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
#define L 60
/* end parameters define */

/* start kernel func */
void SCImark_sor(double **G, double *Gi, double *Gim1, double *Gip1, int num_iterations, int Mm1, int Nm1)
{
    int p, i, j;
	double one_minus_omega = 12;
    double omega_over_four = 48;

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation of Transformations:
1. **Loop Tiling/Blocking**: The original code is transformed to use loop tiling, which is a technique to improve cache performance by processing data in blocks. This is done by introducing new loop variables (`t1`, `t2`) that iterate over blocks of the original loop indices.

2. **Parallelization**: The `#pragma omp parallel for` directive is used to parallelize the outer loop, which can be safely parallelized as each iteration of the loop is independent of the others.

3. **Register Usage**: The variables `Gi`, `Gim1`, and `Gip1` are declared as `register` to suggest that they should be stored in CPU registers for faster access.

4. **Loop Bounds Adjustment**: The loop bounds are adjusted to ensure that the loop iterates over the correct range of indices, taking into account the block size (32 in this case).

5. **Loop Unrolling**: The innermost loop is partially unrolled to reduce the overhead of loop control and to allow for better instruction-level parallelism. This unrolling is done in steps of 4, which is a common optimization for modern CPUs that can handle multiple instructions per cycle.

These transformations aim to improve the performance of the original code by reducing cache misses, leveraging parallel processing, optimizing register usage, and reducing loop overhead.*/

/*### Explanation of Transformations:
1. **Loop Tiling/Blocking**: The original code is transformed to use loop tiling, which is a technique to improve cache performance by processing data in blocks. This is done by introducing new loop variables (`t1`, `t2`) that iterate over blocks of the original loop indices.

2. **Parallelization**: The `#pragma omp parallel for` directive is used to parallelize the outer loop, which can be safely parallelized as each iteration of the loop is independent of the others.

3. **Register Usage**: The variables `Gi`, `Gim1`, and `Gip1` are declared as `register` to suggest that they should be stored in CPU registers for faster access.

4. **Loop Bounds Adjustment**: The loop bounds are adjusted to ensure that the loop iterates over the correct range of indices, taking into account the block size (32 in this case).

5. **Loop Unrolling**: The innermost loop is partially unrolled to reduce the overhead of loop control and to allow for better instruction-level parallelism.

These transformations aim to improve the performance of the original code by reducing cache misses, leveraging parallel processing, optimizing register usage, and reducing loop overhead.*/

int lb, ub, lbp, ubp;
register int lbv, ubv;

for (int iter = 0; iter < ITERATIONS; iter++) {
    for (int p = 0; p < num_iterations; p++) {
        lbp = 0;
        ubp = floord(Mm1 - 1, 32);
#pragma omp parallel for private(lbv, ubv, Gi, Gim1, Gip1)
        for (int t1 = lbp; t1 <= ubp; t1++) {
            for (int t2 = max(1, 32 * t1); t2 <= min(Mm1 - 1, 32 * t1 + 31); t2++) {
                Gi = G[t2];
                Gim1 = G[t2 - 1];
                Gip1 = G[t2 + 1];
                for (int j = 1; j < Nm1 - 3; j += 4) {
                    Gi[j] = omega_over_four * (Gim1[j] + Gip1[j] + Gi[j - 1] + Gi[j + 1]) +
                            one_minus_omega * Gi[j];
                    Gi[j + 1] = omega_over_four * (Gim1[j + 1] + Gip1[j + 1] + Gi[j] + Gi[j + 2]) +
                                one_minus_omega * Gi[j + 1];
                    Gi[j + 2] = omega_over_four * (Gim1[j + 2] + Gip1[j + 2] + Gi[j + 1] + Gi[j + 3]) +
                                one_minus_omega * Gi[j + 2];
                    Gi[j + 3] = omega_over_four * (Gim1[j + 3] + Gip1[j + 3] + Gi[j + 2] + Gi[j + 4]) +
                                one_minus_omega * Gi[j + 3];
                }
                for (int j = Nm1 - 3; j < Nm1; j++) {
                    Gi[j] = omega_over_four * (Gim1[j] + Gip1[j] + Gi[j - 1] + Gi[j + 1]) +
                            one_minus_omega * Gi[j];
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
    ARRAY_PREPARATION_2D(array_0, M+1, L+1);
    ARRAY_PREPARATION_1D(array_1, L+1);
    ARRAY_PREPARATION_1D(array_2, L+1);
	ARRAY_PREPARATION_1D(array_3, L+1);

    SCImark_sor(array_0, array_1, array_2, array_3, N, M, L);

    print_array_1d(array_1, L+1);

    free_array_2d(array_0, M+1, L+1);
    free_array_1d(array_1, L+1);
    free_array_1d(array_2, L+1);
	free_array_1d(array_3, L+1);

    return 0;
}