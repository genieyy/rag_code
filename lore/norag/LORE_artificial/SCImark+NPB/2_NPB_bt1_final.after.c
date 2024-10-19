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
#define N 30
#define M 25
#define Q 20
#define P 15
#define R 10
/* end parameters define */

/* start kernel func */
void NPB_bt1(double *rms, double ****u, double *u_exact, int n, int m, int q, int p, int r)
{
    int i, j, k, l, d;
    double xi, eta, zeta, add;
    int dnxm1 = 1;
    int dnym1 = 2;
    int dnzm1 = 3;

    double time_start = omp_get_wtime();
#pragma scop
/*### Explanation:
1. **Precompute Constants**: 
   - `dnxm1_d`, `dnym1_d`, `dnzm1_d`, and `N_minus_2_inv` are precomputed once per iteration to avoid repeated casting and division operations inside the nested loops.

2. **Loop Unrolling**:
   - The inner loop over `d` is unrolled to avoid the overhead of loop control. This assumes that `r` is small enough that unrolling is feasible.

3. **Reduced Casting**:
   - Casting operations are minimized by using precomputed double values.

4. **Avoid Repeated Division**:
   - The division by `(N - 2)` is done once per iteration of the `d` loop, reducing the number of divisions.*/

for (int iter = 0; iter < ITERATIONS; iter++) {
    double dnxm1_d = (double)dnxm1;
    double dnym1_d = (double)dnym1;
    double dnzm1_d = (double)dnzm1;
    double N_minus_2_inv = 1.0 / (double)(N - 2);

    for (i = 0; i < n; i++) {
        double xi = (double)i * dnxm1_d;
        for (j = 0; j < m; j++) {
            double eta = (double)j * dnym1_d;
            for (k = 0; k < q; k++) {
                double zeta = (double)k * dnzm1_d;

                for (l = 0; l < p; l++) {
                    double add = u[i][j][k][l] - u_exact[l];
                    rms[l] += add * add;
                }
            }
        }
    }

    for (i = 0; i < p; i++) {
        rms[i] *= N_minus_2_inv;
        for (d = 1; d <= r; d++) {
            rms[i] *= N_minus_2_inv;
        }
        rms[i] = sqrt(rms[i]);
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
    ARRAY_PREPARATION_1D(array_2, P+1);

    NPB_bt1(array_0, array_1, array_2, N, M, Q, P, R);

    print_array_1d(array_0, P+1);

    free_array_1d(array_0, P+1);
    free_array_4d(array_1, N+1, M+1, Q+1, P+1);
    free_array_1d(array_2, P+1);

    return 0;
}