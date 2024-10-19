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
#define N 10
#define M 12
#define L 14
#define P 16
/* end parameters define */

/* start kernel func */
void NPB_bt18(double *****lhs, double ****rhs, int nz, int mz, int q, int BLOCK_SIZE)
{
    int i, j, k, m, n;

    double time_start = omp_get_wtime();
#pragma scop
/**/

int t1, t2, t3, t4, t5, t6, t7, t8;
int lb, ub, lbp, ubp, lb2, ub2;
register int lbv, ubv;

lbp = 0;
ubp = floord(ITERATIONS - 1, 32);
#pragma omp parallel for private(lbv, ubv, t3, t4, t5, t6, t7, t8)
for (int t1 = lbp; t1 <= ubp; t1++) {
    for (int t2 = 0; t2 <= floord(nz - 3, 32); t2++) {
        for (int t3 = 0; t3 <= floord(mz - 3, 32); t3++) {
            for (int t4 = 0; t4 <= floord(q - 2, 32); t4++) {
                for (int t5 = max(1, 32 * t2); t5 <= min(nz - 2, 32 * t2 + 31); t5++) {
                    for (int t6 = max(1, 32 * t3); t6 <= min(mz - 2, 32 * t3 + 31); t6++) {
                        for (int t7 = max(0, 32 * t4); t7 <= min(q - 3, 32 * t4 + 31); t7++) {
                            for (int t8 = 0; t8 < BLOCK_SIZE; t8++) {
                                lbv = 0;
                                ubv = BLOCK_SIZE - 1;
#pragma ivdep
#pragma vector always
                                for (int n = lbv; n <= ubv; n++) {
                                    rhs[t5][t6][t7][t8] = rhs[t5][t6][t7][t8] - lhs[t5][t6][t7][t8][n] * rhs[t5][t6][t7 + 1][n];
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
    ARRAY_PREPARATION_5D(array_0, N + 1, M + 1, L + 1, P + 1, P + 1);
    ARRAY_PREPARATION_4D(array_1, N + 1, M + 1, L + 1, P + 1);

    NPB_bt18(array_0, array_1, N, M, L, P);

    print_array_4d(array_1, N + 1, M + 1, L + 1, P + 1);

    free_array_5d(array_0, N + 1, M + 1, L + 1, P + 1, P + 1);
    free_array_4d(array_1, N + 1, M + 1, L + 1, P + 1);

    return 0;
}