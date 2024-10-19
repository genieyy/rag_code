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
#define N 400
#define M 600
/* end parameters define */

/* start kernel func */
void SCImark_lu1(double **A, double *Aii, double *Aj, int m, int n)
{
    int j = 0;
    int ii, jj;
    double AiiJ;

    double time_start = omp_get_wtime();
#pragma scop
/**/

int t1, t2, t3, t4;
int lb, ub, lbp, ubp;
register int lbv, ubv;

for (int t1 = 0; t1 <= floord(ITERATIONS + m - j - 2, 32); t1++) {
    lbp = max(0, ceild(32 * t1 - m + j + 1, 32));
    ubp = min(floord(ITERATIONS - 1, 32), t1);
#pragma omp parallel for private(lbv, ubv, t3, t4)
    for (int t2 = lbp; t2 <= ubp; t2++) {
        for (int t3 = max(32 * t1 - 32 * t2, j + 1); t3 <= min(m - 1, 32 * t1 - 32 * t2 + 31); t3++) {
            double *Aii = A[t3];
            double *Aj = A[j];
            double AiiJ = Aii[j];
            lbv = j + 1;
            ubv = n - 1;
#pragma ivdep
#pragma vector always
            for (int t4 = lbv; t4 <= ubv; t4++) {
                Aii[t4] -= AiiJ * Aj[t4];
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
    ARRAY_PREPARATION_2D(array_0, N+1, M+1);
    ARRAY_PREPARATION_1D(array_1, M+1);
    ARRAY_PREPARATION_1D(array_2, M+1);

    SCImark_lu1(array_0, array_1, array_2, N, M);

    print_array_1d(array_1, M+1);

    free_array_2d(array_0, N+1, M+1);
    free_array_1d(array_1, M+1);
    free_array_1d(array_2, M+1);

    return 0;
}
