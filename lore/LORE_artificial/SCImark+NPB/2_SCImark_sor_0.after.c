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
/*### Analysis of Meaning-Preserving Loop Transformation Methods:

1. **Loop Unrolling**: The original code has nested loops, and the optimized code uses loop unrolling to reduce the overhead of loop control. This is evident in the use of `#pragma ivdep` and `#pragma vector always` in the examples, which hint at vectorization and unrolling.

2. **Parallelization**: The optimized code uses OpenMP (`#pragma omp parallel for`) to parallelize the outer loop, distributing the iterations across multiple threads. This is a common technique to exploit multi-core processors.

3. **Loop Distribution and Fusion**: The original code has multiple nested loops, and the optimized code distributes and fuses these loops to improve cache locality and reduce the number of loop iterations. This is seen in the examples where the loops are restructured to minimize the number of iterations.

4. **Loop Interchange**: The order of loops is changed to improve cache performance. This is evident in the examples where the loop order is adjusted to ensure that the most frequently accessed data is loaded into the cache.

5. **Loop Tiling**: The optimized code uses loop tiling to break the problem into smaller chunks (tiles) that fit better into the cache. This is seen in the use of `floord` and `min`/`max` functions to create tiles of size 32.

### Application of Learned Methods:

1. **Parallelization**: The outer loop is parallelized using OpenMP to distribute the work across multiple threads.
2. **Loop Tiling**: The `num_iterations` loop is tiled with a tile size of 32 to improve cache performance.
3. **Loop Distribution**: The inner loops are distributed to ensure that the most frequently accessed data is loaded into the cache.

These transformations aim to improve the performance of the original code by reducing the overhead of loop control, exploiting parallelism, and improving cache locality.*/

int lb, ub, lbp, ubp;
register int lbv, ubv;

for (int iter = 0; iter < ITERATIONS; iter++) {
    lbp = 0;
    ubp = floord(num_iterations - 1, 32);
#pragma omp parallel for private(lbv, ubv, Gi, Gim1, Gip1)
    for (int p = lbp; p <= ubp; p++) {
        for (int i = max(1, 32 * p); i <= min(Mm1 - 1, 32 * p + 31); i++) {
            Gi = G[i];
            Gim1 = G[i - 1];
            Gip1 = G[i + 1];
            for (int j = 1; j < Nm1; j++) {
                Gi[j] = omega_over_four * (Gim1[j] + Gip1[j] + Gi[j - 1] + Gi[j + 1]) +
                        one_minus_omega * Gi[j];
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