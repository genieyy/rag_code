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
/*### Explanation of the Optimization:
1. **Loop Tiling and Parallelization**:
   - The outer loop over `p` is tiled and parallelized using OpenMP to distribute the work across multiple threads. This is done by calculating the bounds `lbp` and `ubp` to ensure that each thread processes a chunk of the loop iterations.
   - The `#pragma omp parallel for` directive is used to parallelize the loop, with private variables `lbv`, `ubv`, `Gi`, `Gim1`, and `Gip1` to avoid race conditions.

2. **Loop Unrolling**:
   - The inner loops are not unrolled in this example, but the tiling approach helps in reducing cache misses and improving locality of reference.

3. **Reduction of Redundant Computations**:
   - The variables `Gi`, `Gim1`, and `Gip1` are precomputed outside the innermost loop to avoid redundant memory accesses.

This optimization leverages the techniques observed in the provided examples, such as loop tiling and parallelization, to improve the performance of the given code.*/

int lb, ub, lbp, ubp;
register int lbv, ubv;

for (int iter = 0; iter < ITERATIONS; iter++) {
    lbp = 0;
    ubp = floord(Mm1 - 2, 32);
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